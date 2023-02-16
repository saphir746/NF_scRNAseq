library(sceasy)
library(reticulate)
library(tidyverse)
library(ggpubr)
library(Seurat)

sc <- import("scanpy")

#files_dir="/camp/stp/babs/working/schneid/projects/sahaie/giovanni.giangreco/Characterisation_of_CAF_in_HPV_cancer_scrnaseq/Data_interim_files/"

Files=list.files(, pattern = '_scrubletted.h5ad')


seurat_list <-	map(Files, function(f) {
  atlas.data <- sc$read_h5ad(f)
  counts <- t(as.matrix(atlas.data$X))
  colnames(counts) <-  atlas.data$obs_names$to_list()
  rownames(counts) <-  atlas.data$var_names$to_list()
  counts <- Matrix::Matrix(as.matrix(counts), sparse = T)
  # 
  seur_obj <- CreateSeuratObject(counts)
  seur_obj <- AddMetaData(seur_obj,  atlas.data$obs)
  id<-gsub('_scrubletted.h5ad','',f)
  seur_obj@meta.data$orig.ident<-id
  seur_obj
})

names(seurat_list)<-gsub('_scrubletted.h5ad','',Files)

saveRDS(
  seurat_list,
  file = "SC21137_raw_scrubletted.RDS"
)
