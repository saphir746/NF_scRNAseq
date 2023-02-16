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
suffix<-str_match(fileRDS,'SC21137_Integrated_(.*).RDS')[,2]

######################################

Everything.combined <-read_rds(fileRDS)

######################################

  # Run the standard workflow for visualization and clustering

Everything.combined <- ScaleData(Everything.combined, verbose = FALSE)
stuff<-VariableFeatures(object = Everything.combined)
Everything.combined <- RunPCA(Everything.combined, features = stuff, verbose = FALSE)
int_dim <- intrinsicDimension::maxLikGlobalDimEst(Everything.combined@reductions$pca@cell.embeddings, k = 30)
D<-ceiling(int_dim$dim.est)
Everything.combined <- FindNeighbors(Everything.combined, dims = 1:D)
Everything.combined <- RunUMAP(Everything.combined, reduction = "pca", dims = 1:D)
Everything.combined <- FindNeighbors(Everything.combined, reduction = "pca", dims = 1:D)
Everything.combined <- FindClusters(Everything.combined, resolution = seq(0.1,1.1,0.2))#0.5)
  

final.name<-paste0('SC21137_Integrated_Filtered_',suffix,'.RDS')
saveRDS(
  Everything.combined,
  file = final.name
  #paste0(dir_all,"/",final.name)
)