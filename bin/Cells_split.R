#!/usr/bin/env Rscript

library(tidyverse)
library(ggpubr)
library(Seurat)
library(SeuratObject)


# GEN.dir<-"/camp/stp/babs/working/schneid/"
# #"/home/schneid/Documents/CAMP/"
# #"/home/deborah/Documents/Crick/Projects/"
# 
# dir_all_sub=paste0(GEN.dir,
#                    "projects/sahaie/giovanni.giangreco/",
#                    "Characterisation_of_CAF_in_HPV_cancer_scrnaseq")
# dir_all=paste0(dir_all_sub,"/Data_interim_files")
# Everything.combined <-
#   read_rds(paste0(dir_all,"/SC21137_Integrated_Filtered_augmented.RDS"))

args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]

Everything.combined<-read_rds(fileRDS)

#####

Cancer_cells<-subset(Everything.combined, subset = Cell.ident == 'Cancer_cells')
Cancer_cells %>% saveRDS(.,
  file = "SC21137_Integrated_cancer_cells.RDS"
  #paste0(dir_all,"/SC21137_Integrated_Cancer_cells.RDS")
)
#
Immune_cells<-subset(Everything.combined, subset = Cell.ident == 'Immune_cells')
Immune_cells %>% saveRDS(.,
                         file = "SC21137_Integrated_Immune_cells.RDS"
                         #paste0(dir_all,"/SC21137_Integrated_Immune_cells.RDS")
)
#
Stromal_cells<-subset(Everything.combined, subset = Cell.ident == 'Stromal_cells')
Stromal_cells %>% saveRDS(.,
                         file = "SC21137_Integrated_Stromal_cells.RDS"
                         #paste0(dir_all,"/SC21137_Integrated_Stromal_cells.RDS")
)

#
CAFs_cells<-subset(Everything.combined, subset = Cell.ident == 'CAFs')
CAFs_cells %>% saveRDS(.,
                          file = "SC21137_Integrated_CAFs_cells.RDS"
                          #paste0(dir_all,"/SC21137_Integrated_Stromal_cells.RDS")
)