#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)

# 
# GEN.dir<-"/camp/stp/babs/working/schneid/"
# #"/home/schneid/Documents/CAMP/"
# #"/home/deborah/Documents/Crick/Projects/"
# 
# dir_all_sub=paste0(GEN.dir,
#                    "projects/sahaie/giovanni.giangreco/",
#                    "Characterisation_of_CAF_in_HPV_cancer_scrnaseq")
# 
# dir_all=paste0(dir_all_sub,"/Data_interim_files")
# 
# Everything.combined <-
#   read_rds(paste0(dir_all,"/SC21137_Integrated_Filtered_withMOC2_3.RDS"))


args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]
resol<-args[2]

fileName<-str_match(fileRDS,'(.*).RDS')[,2]

R<-paste0('integrated_snn_res.',resol)
r<-gsub('0.','',resol)

Everything.combined <-read_rds(fileRDS)

DefaultAssay(Everything.combined)<-'RNA'


Everything.combined <- SetIdent(Everything.combined,
                               value = Everything.combined@meta.data[,R])

pbmc.markers <- FindAllMarkers(Everything.combined,
                              only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25,
                              #test.use="DESeq2",
                              random.seed = 1, return.thresh = 0.01)

pbmc.markers %>% write.csv(.,file = paste0(fileName,"_lead_markers_0",r,".csv"))
                             #paste0(dir_all,"/lead_markers_03.csv"))


# Everything.combined <- SetIdent(Everything.combined, 
#                                 value = Everything.combined@meta.data$integrated_snn_res.0.9)
# 
# pbmc.markers <- FindAllMarkers(Everything.combined, 
#                                only.pos = TRUE, min.pct = 0.25, 
#                                logfc.threshold = 0.25,
#                               # test.use="DESeq2",
#                                random.seed = 1, return.thresh = 0.01)
# 
# pbmc.markers %>% write.csv(.,file = "lead_markers_09.csv")
#                              #paste0(dir_all,"/lead_markers_09.csv"))

