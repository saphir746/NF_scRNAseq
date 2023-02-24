#!/usr/bin/env Rscript

library(sceasy)
library(reticulate)
#library(tidyverse)
#library(ggpubr)
library(Seurat)


# GEN.dir<-"/camp/stp/babs/working/schneid/"
#   #"/home/schneid/Documents/CAMP/"
# #"/home/deborah/Documents/Crick/Projects/"
# 
# dir_all_sub=paste0(GEN.dir,
#                    "projects/sahaie/giovanni.giangreco/",
#                    "Characterisation_of_CAF_in_HPV_cancer_scrnaseq")
# #
# dir_all=paste0(dir_all_sub,"/Data_interim_files")
# # # 
# seurat_list <-
#   read_rds(paste0(dir_all,"/SC21137_raw_seurat_object.RDS"))


args = commandArgs(trailingOnly=TRUE)

seurat_list <-readRDS(args[1])

lapply(names(seurat_list), function(x){
seurat_obj<-seurat_list[[x]]
sceasy::convertFormat(seurat_obj, from="seurat", to="anndata",
                       outFile=paste0(x,'_raw.h5ad'))
			#paste0(dir_all,'/',x,'_raw.h5ad'))
                       })
