#!/usr/bin/env Rscript

library(tidyverse)
library(ggpubr)
#library(SingleCellExperiment)
#library(scDblFinder)
library(Seurat)
library(SeuratObject)
library(intrinsicDimension)

thrsh.mt<-20

# GEN.dir<-"/camp/stp/babs/working/schneid/"
# #"/home/schneid/Documents/CAMP/"
# #"/home/deborah/Documents/Crick/Projects/"
# 
# dir_all_sub=paste0(GEN.dir,
#                    "projects/sahaie/giovanni.giangreco/",
#                    "Characterisation_of_CAF_in_HPV_cancer_scrnaseq")
# dir_all=paste0(dir_all_sub,"/Data_interim_files")

args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]

######################################

Umap_seurat_list<-read_rds(fileRDS)
# paste0(dir_all,"/SC21137_Umapped_Filteredseurat_object.RDS"))
# exclude "MOC2_3" due to poor QC everything
#Umap_seurat_list.bis<-Umap_seurat_list[names(Umap_seurat_list)!="MOC2_3"]

######################################

do_the_things<-function(seur_obj){
  seur_obj<-NormalizeData(seur_obj)
  seur_obj<-FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(seur_obj)
  seur_obj <- ScaleData(seur_obj, features = all.genes)
  seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj ))
  int_dim <- intrinsicDimension::maxLikGlobalDimEst(seur_obj@reductions$pca@cell.embeddings, k = 15)
  D<-ceiling(int_dim$dim.est)
  seur_obj <- FindNeighbors(seur_obj, dims = 1:D)
  seur_obj <- FindClusters(seur_obj, resolution = 0.5)#seq(0.3,1.1,0.2))
  seur_obj <- RunUMAP(seur_obj, dims = 1:D)
  return(seur_obj)
}

convertHumanGeneList <- function(x){
  require("biomaRt")
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
  mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes.mm10<-convertHumanGeneList(s.genes)
g2m.genes.mm10<-convertHumanGeneList(g2m.genes)

do_all_of_it<-function(Umap_seurat_list){
  normalised_seurat_list<-  map(Umap_seurat_list, function(seur_obj) { do_the_things(seur_obj)})
  
  #### Integration round 1 
  
  features <- SelectIntegrationFeatures(object.list = normalised_seurat_list)
  ## remove cc genes from potential anchors 
  features <- features[!(features %in% c(s.genes.mm10,g2m.genes.mm10))]
  #
  ### Integration_large_datasets =https://satijalab.org/seurat/articles/integration_large_datasets.html 
  normalised_seurat_list <- lapply(X = normalised_seurat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  Anchors <- FindIntegrationAnchors(object.list = normalised_seurat_list, reference=c(1,2,7),  reduction = "rpca",
                                    dims = 1:50)
  Everything.combined <- IntegrateData(Anchors, normalization.method = "LogNormalize",  dims = 1:50, verbose = TRUE)
  DefaultAssay(Everything.combined) <- "integrated"
  
  cell_type<-Everything.combined[[]]$orig.ident %>% gsub('_\\d','',.)
  Everything.combined@meta.data$Cell.type<-cell_type
  
  # Run the standard workflow for visualization and clustering
  all.genes <- rownames(Everything.combined)
  Everything.combined <- CellCycleScoring(Everything.combined, s.features = s.genes.mm10,
                                          g2m.features = g2m.genes.mm10, set.ident = FALSE)
  Everything.combined   <- ScaleData(Everything.combined,
                                     vars.to.regress = c("S.Score", "G2M.Score"), 
                                     features = all.genes)
  #
  #Everything.combined <- ScaleData(Everything.combined, verbose = FALSE)
  stuff<-VariableFeatures(object = Everything.combined)
  Everything.combined <- RunPCA(Everything.combined, features = stuff, verbose = FALSE)
  #
  int_dim <- intrinsicDimension::maxLikGlobalDimEst(Everything.combined@reductions$pca@cell.embeddings, k = 30)
  D<-ceiling(int_dim$dim.est)
  Everything.combined <- FindNeighbors(Everything.combined, dims = 1:D)
  Everything.combined <- RunUMAP(Everything.combined, reduction = "pca", dims = 1:D)
  Everything.combined <- FindNeighbors(Everything.combined, reduction = "pca", dims = 1:D)
  Everything.combined <- FindClusters(Everything.combined, resolution = seq(0.1,1.1,0.2))#0.5)
  
  ##
  
  Everything.combined[["percent.ptprc"]] <- PercentageFeatureSet(Everything.combined, pattern = "Ptprc", assay = 'RNA')
  Everything.combined@meta.data$Ptprc_detect <-ifelse(Everything.combined@meta.data$percent.ptprc>0,1,0)
  Everything.combined@meta.data$Ptprc_tdT_pos<-ifelse((Everything.combined@meta.data$Ptprc_detect==1)&
                                                        (Everything.combined@meta.data$Cell.orig=="Mouse_host"),'yes','no')
  Everything.combined@meta.data$Ptprc_tdT_pos<-as.factor(Everything.combined@meta.data$Ptprc_tdT_pos)
  
  ###### Cell -cycle sorting
  # ### done in previous scripts ###############
  
  return(Everything.combined)
}

# Everything.combined<-do_all_of_it(Umap_seurat_list.bis)
# saveRDS(
#   Everything.combined,
#   file = paste0(dir_all,"/SC21137_Integrated_Filtered.RDS")
# )

Everything.combined.all<-do_all_of_it(Umap_seurat_list)
saveRDS(
  Everything.combined.all,
  file = "SC21137_Integrated_Filtered.RDS"
  #paste0(dir_all,"/SC21137_Integrated_Filtered_withMOC2_3.RDS")
)