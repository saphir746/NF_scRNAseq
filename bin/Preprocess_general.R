#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(SeuratObject)
library(intrinsicDimension)
library(plyr)
library(SingleCellExperiment)
library(scDblFinder)


# args = commandArgs(trailingOnly=TRUE)
# 
# fileRDS<-args[1]

dir_all_sub<-"/nemo/stp/babs/working/schneid/projects/haydaya/leticia.monin/lm510/"
dir_all<-paste0(dir_all_sub,'Data_interim/')

fileRDS<-paste0(dir_all,"/lm510_raw_seurat_object.RDS")

seurat_list <-read_rds(fileRDS)
# fileRDS=SC21137_raw_scrubletted.RDS
thrsh.mt=20


cc.genes.mm10<-read_rds(paste0(dir_all_sub,'NF_scRNAseq/assets/cc_genes_mm10.rds'))

doublet_ident<-function(seur_obj){
  Obj.sce<-as.SingleCellExperiment(seur_obj,assay='RNA')
  sce <- scDblFinder(Obj.sce,dbr=0.05,dbr.sd=0, clusters=seur_obj@meta.data$seurat_clusters)
  sce@colData %>% as.data.frame() ->metadata.new
  seur_obj@meta.data<-metadata.new
  return(seur_obj)
}

################################### Normalisation + Scale + PCA
normalised_seurat_list<-  map(seurat_list, function(seur_obj) {
  seur_obj<-NormalizeData(seur_obj)
  seur_obj<-FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(seur_obj)
  seur_obj <- ScaleData(seur_obj, features = all.genes)
  seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj ))
  seur_obj <- CellCycleScoring(seur_obj, s.features = cc.genes.mm10$s,
                               g2m.features = cc.genes.mm10$g2m, set.ident = FALSE)
  #all.genes <- rownames(seur_obj)
  #seur_obj  <- ScaleData(seur_obj , vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
  #seur_obj <- RunPCA(seur_obj, features = VariableFeatures(seur_obj))
  # seur_obj <- RunPCA(seur_obj, features = c(s.genes.mm10, g2m.genes.mm10))
  seur_obj
})

################################### Clustering + UMAP 
Umap_seurat_list<-map(normalised_seurat_list, function(pbmc) {
  int_dim <-
    intrinsicDimension::maxLikGlobalDimEst(pbmc@reductions$pca@cell.embeddings, k = 10)
  D<-ceiling(int_dim$dim.est)
  pbmc <- FindNeighbors(pbmc, dims = 1:D)
  pbmc <- FindClusters(pbmc, resolution = 0.5)# seq(0.3,1.1,0.2))
  pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  pbmc<-doublet_ident(pbmc)
  pbmc
})



################################## Save objects 1
saveRDS(
  Umap_seurat_list,
  file = #"SC21137_Umapped_seurat_object.RDS"
  paste0(dir_all,"/lm510_Umapped_seurat_object.RDS")
)
