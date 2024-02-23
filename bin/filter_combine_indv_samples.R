#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(SeuratObject)
library(intrinsicDimension)
library(plyr)
library(scuttle)
library(reticulate)
use_condaenv("/nemo/stp/babs/working/software/anaconda/envs/scvelo-0.2.2")


GEN.dir<-"/nemo/stp/babs/working/schneid/projects/"
#"/home/schneid/Documents/nemo/"
#"/home/deborah/Documents/Crick/Projects/"

dir_all_sub=paste0(GEN.dir,
                   "/haydaya/leticia.monin/lm510/")

dir_all=paste0(dir_all_sub,"/Data_interim")

fileRDS<-paste0(dir_all,"/lm510_Umapped_seurat_object.RDS")

seurat_list<-readRDS(fileRDS)

thrsh.mt<-20
thrsh.dbl<-50

cc.genes.mm10<-read_rds(paste0(dir_all_sub,'NF_scRNAseq/assets/cc_genes_mm10.rds'))
s.genes.mm10<-cc.genes.mm10[['s']]
g2m.genes.mm10<-cc.genes.mm10[['g2m']]

################################
################################## Filtering - mt content
################################

###### Data driven thresholding 
# 
# sum.stats<-map(seurat_list, function(seur_obj){
#   df.dat<-FetchData(object=seur_obj, vars=c("orig.ident","nFeature_RNA", "nCount_RNA",
#                                             "percent_ribosomal","percent.mt")) %>% as.data.frame()
# })
# mad.stats<-map(seurat_list, function(seur_obj){
#   df.dat<-FetchData(object=seur_obj, vars=c("orig.ident","nFeature_RNA", "nCount_RNA",
#                                             "percent_ribosomal","percent.mt")) %>% as.data.frame() %>%
#     summarise_if(is.numeric, ~mad(.x))
# })
# median.stats<-map(seurat_list, function(seur_obj){
#   df.dat<-FetchData(object=seur_obj, vars=c("orig.ident","nFeature_RNA", "nCount_RNA",
#                                             "percent_ribosomal","percent.mt")) %>% as.data.frame() %>%
#     summarise_if(is.numeric, ~median(.x))
# })
# 
# Df.summary <- do.call("rbind", sum.stats)
# Df.mad<- do.call("rbind", mad.stats)
# Df.median<- do.call("rbind", median.stats)
# 
# Df.summary %>% dplyr::select(-orig.ident) %>% 
#   summarise_if(is.numeric, ~quantile(.x, probs = c(0,0.025,0.05,0.1,0.5,0.9,0.95,1))) -> table1
# 
# th.nFeature.down<-min(200,round_any(table1$nFeature_RNA[2], 100, f = floor))
# th.nCount<-round_any(table1$nCount_RNA[2], 100, f = floor)
# thrsh.mt<-min(thrsh.mt,table1$percent.mt[6])
# 
# cat('nFeature thrsh.: ',th.nFeature.down,' \n')
# cat('nCount thrsh.: ',th.nCount,' \n')
# cat('Mito. % thrsh.: ',thrsh.mt,' \n')

Umap_seuratFilt_list<-map(seurat_list, function(seur_obj) {
  df.dat<-FetchData(object=seur_obj, vars=c('RNA_snn_res.0.5',"nFeature_RNA", "nCount_RNA",
                                            "percent_ribosomal","percent.mt","scDblFinder.class"))
  df.dat %>% summarise_if(is.numeric, ~mad(.x)) ->Df.mad
  df.dat %>% summarise_if(is.numeric, ~median(.x)) ->Df.median
  df.dat %>% summarise_if(is.numeric, ~quantile(.x, probs = c(0,0.025))) ->Df.quant
  
  Thr.high<-ceiling(Df.mad$nFeature_RNA*3+Df.median$nFeature_RNA)
  Thr.low<-floor(-Df.mad$nFeature_RNA*3+Df.median$nFeature_RNA) %>% max(.,Df.quant$nFeature_RNA[2])
  
  df.display<-df.dat %>% group_by(RNA_snn_res.0.5) %>% 
    dplyr::summarise(doublet.perc=sum(scDblFinder.class == 'doublet')/n()*100, 
                     perc_Mito = mean(percent.mt)) %>%
    as.data.frame()
  colnames(df.display)[1]<-'Cluster'
  mt_excl<-levels(droplevels(df.display[df.display$perc_Mito>thrsh.mt,]$Cluster))
  dlt_excl<-levels(droplevels(df.display[df.display$doublet.perc>thrsh.dbl,]$Cluster))
  Umap_filtered<-subset(seur_obj, subset = seurat_clusters %in% mt_excl, invert= TRUE)
  Umap_filtered<-subset(Umap_filtered, subset = seurat_clusters %in% dlt_excl, invert= TRUE)
  Umap_filtered2<-subset(Umap_filtered, subset = percent.mt < thrsh.mt )
  Umap_filtered2<-subset(Umap_filtered2, subset = scDblFinder.class == "singlet" )
  Umap_filtered2<-subset(Umap_filtered2, subset = (nFeature_RNA > Thr.low)&(nFeature_RNA < Thr.high) )
  Umap_filtered2@meta.data %<>% dplyr::select("orig.ident","nCount_RNA","nFeature_RNA","group",
                                              "sorting_batch","mouse_id","percent.mt","percent_ribosomal","percent.eGFP",
                                              "S.Score","G2M.Score","Phase","RNA_snn_res.0.5","seurat_clusters","scDblFinder.class")%>%
    mutate(nFeature_high=Thr.high,nFeature_low=Thr.low )
  Umap_filtered2
})

do_the_things<-function(seur_obj){
  seur_obj<-NormalizeData(seur_obj)
  seur_obj<-FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 3000)
  seur_obj <- CellCycleScoring(seur_obj, s.features = s.genes.mm10,
                               g2m.features = g2m.genes.mm10, set.ident = FALSE)
  all.genes <- rownames(seur_obj)
  seur_obj <- ScaleData(seur_obj, features = all.genes)
  seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj ))
  int_dim <- intrinsicDimension::maxLikGlobalDimEst(seur_obj@reductions$pca@cell.embeddings, k = 15)
  D<-ceiling(int_dim$dim.est)
  seur_obj <- FindNeighbors(seur_obj, dims = 1:D)
  seur_obj <- FindClusters(seur_obj, resolution = seq(0.3,1.1,0.2))#, method = "igraph", algorithm = 4, verbose=TRUE)
  seur_obj <- RunUMAP(seur_obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated.QCed")
  return(seur_obj)
}


Umap_seuratFilt_list<-lapply(Umap_seuratFilt_list, function(X) do_the_things(X))


################################## Save objects 2
saveRDS(
  Umap_seuratFilt_list,
  file = #"SC21137_Umapped_Filteredseurat_object.RDS"
    paste0(dir_all,"/lm510_UmappedFiltered_seurat_object.RDS")
)
