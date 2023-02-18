#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(SeuratObject)
library(intrinsicDimension)
library(plyr)



Cells_file<-list.files(pattern="SC21137_Umapped_(.*).RDS")

cell.names<-lapply(Cells_file, function(cell) str_match(cell,"SC21137_Umapped_(.*).RDS")[,2])

seurat_list<-lapply(Cells_file, function(cell) read_rds(cell))
names(seurat_list)<-cell.names

thrsh.mt<-20

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

Umap_seuratFilt_list<-map(seurat_list, function(seur_obj) {
  df.dat<-FetchData(object=seur_obj, vars=c('RNA_snn_res.0.5',"nFeature_RNA", "nCount_RNA","percent.tdT",
                                            "percent_ribosomal","percent.mt","scrublet_pred")) #"scrublet_score"))
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
