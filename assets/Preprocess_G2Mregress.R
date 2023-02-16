
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(intrinsicDimension)
library(plyr)


args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]

seurat_list <-read_rds(fileRDS)
# fileRDS=SC21137_raw_scrubletted.RDS
thrsh.mt=20

###########

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
convertHumanGeneList <- function(x){
  require("biomaRt")
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
  mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}
s.genes.mm10<-convertHumanGeneList(s.genes)
g2m.genes.mm10<-convertHumanGeneList(g2m.genes)


################################### Normalisation 
normalised_seurat_list<-  map(seurat_list, function(seur_obj) {
  tmp.1 <- ifelse(seur_obj[["percent.tdT"]]>0,"Mouse_host","Cancer_inject")
  seur_obj@meta.data["Cell.orig"] <- tmp.1
  seur_obj<-NormalizeData(seur_obj)
  seur_obj<-FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 3000)
  seur_obj
})

################################### Scale + PCA
# normalised_seurat_list<-  map(normalised_seurat_list, function(seur_obj) {
#   all.genes <- rownames(seur_obj)
#   seur_obj <- ScaleData(seur_obj, features = all.genes)
#   seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj ))
#   seur_obj
# })

################################### Cell-cycle scoring and regression

normalised_seurat_list<-map(normalised_seurat_list,function(seur_obj){
  
  seur_obj <- CellCycleScoring(seur_obj, s.features = s.genes.mm10,
                                          g2m.features = g2m.genes.mm10, set.ident = FALSE)
  all.genes <- rownames(seur_obj)
  seur_obj  <- ScaleData(seur_obj , vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
  seur_obj <- RunPCA(seur_obj, features = VariableFeatures(seur_obj))
 # seur_obj <- RunPCA(seur_obj, features = c(s.genes.mm10, g2m.genes.mm10))
  seur_obj
})



# intrin_dim_list<-map(normalised_seurat_list, function(seur_obj, seur_name) {
#   int_dim <-
#     intrinsicDimension::maxLikGlobalDimEst(seur_obj@reductions$pca@cell.embeddings, k = 15)
#   ceiling(int_dim$dim.est)
# })

################################### Clustering + UMAP 
Umap_seurat_list<-map(normalised_seurat_list, function(pbmc) {
  int_dim <-
    intrinsicDimension::maxLikGlobalDimEst(pbmc@reductions$pca@cell.embeddings, k = 10)
  D<-ceiling(int_dim$dim.est)
  pbmc <- FindNeighbors(pbmc, dims = 1:D)
  pbmc <- FindClusters(pbmc, resolution = 0.5)# seq(0.3,1.1,0.2))
  pbmc <- RunUMAP(pbmc, dims = 1:30)
  pbmc
})



################################## Save objects 1
saveRDS(
  Umap_seurat_list,
  file = "SC21137_Umapped_seurat_object.RDS"
  #paste0(dir_all,"/SC21137_Umapped_Filteredseurat_object.RDS")
)

################################
################################## Filtering - mt content
################################

###### Data driven thresholding 

sum.stats<-map(seurat_list, function(seur_obj){
  df.dat<-FetchData(object=seur_obj, vars=c("orig.ident","nFeature_RNA", "nCount_RNA","percent.tdT",
                                            "percent_ribosomal","percent.mt")) %>% as.data.frame()
})
Df.summary <- do.call("rbind", sum.stats)
Df.summary %>% select(-orig.ident) %>% 
  summarise_if(is.numeric, ~quantile(.x, probs = c(0,0.025,0.05,0.1,0.5,0.9,0.95,1))) -> table1

th.nFeature<-round_any(table1$nFeature_RNA[2], 100, f = floor)
th.nCount<-round_any(table1$nCount_RNA[2], 100, f = floor)
thrsh.mt<-min(thrsh.mt,table1$percent.mt[6])

cat('nFeature thrsh.: ',th.nFeature,' \n')
cat('nCount thrsh.: ',th.nCount,' \n')
cat('Mito. % thrsh.: ',thrsh.mt,' \n')

seurat_list<-map(seurat_list, function(seur_obj) {
  df.dat<-FetchData(object=seur_obj, vars=c('RNA_snn_res.0.5',"nFeature_RNA", "nCount_RNA","percent.tdT",
                                            "percent_ribosomal","percent.mt","scrublet_pred")) #"scrublet_score"))
  ##
  df.display<-df.dat %>% group_by(RNA_snn_res.0.5) %>% 
    dplyr::summarise(perc_ribosomal = mean(percent_ribosomal), 
                     perc_Mito = mean(percent.mt)) %>%
    # perc_tdT = mean(percent.tdT)) %>%
    as.data.frame()
  colnames(df.display)[1]<-'Cluster'
  mt_excl<-levels(droplevels(df.display[df.display$perc_Mito>thrsh.mt,]$Cluster))
  Umap_filtered<-subset(seur_obj, subset = seurat_clusters %in% mt_excl, invert= TRUE)
  Umap_filtered2<-subset(Umap_filtered, subset = percent.mt < thrsh.mt )
  Umap_filtered2<-subset(Umap_filtered2, subset = scrublet_pred == FALSE )
  Umap_filtered2<-subset(Umap_filtered2, subset = nFeature_RNA > th.nFeature )
  Umap_filtered2<-subset(Umap_filtered2, subset = nCount_RNA > th.nCount )
  Umap_filtered2
})

################################## Save objects 2
saveRDS(
  Umap_seuratFilt_list,
  file = "SC21137_Umapped_Filteredseurat_object.RDS"
  #paste0(dir_all,"/SC21137_Umapped_Filteredseurat_object.RDS")
)
