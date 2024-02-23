#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)

GEN.dir<-"/nemo/stp/babs/working/schneid/projects/"
#"/home/schneid/Documents/nemo/"
#"/home/deborah/Documents/Crick/Projects/"

dir_all_sub=paste0(GEN.dir,
                   "/haydaya/leticia.monin/lm510/")

dir_all=paste0(dir_all_sub,"/Data_interim")
dir_all_save=paste0(dir_all,"/LMG_ind")

Seur_list<-
  read_rds(paste0(dir_all,"/lm510_UmappedFiltered_seurat_object.RDS"))

find_LMGs<-function(seur_obj,sample_name){
  
  DefaultAssay(seur_obj)<-'RNA'
  R <- "RNA_snn_res.1.1"
  seur_obj <- SetIdent(seur_obj, value = seur_obj@meta.data[,R])
  
  pbmc.markers <- FindAllMarkers(seur_obj,
                                 only.pos = TRUE, min.pct = 0.25,
                                 logfc.threshold = 0.25,
                                 #test.use="DESeq2",
                                 random.seed = 1, return.thresh = 0.01)
  
  pbmc.markers %>% write.csv(.,file = paste0(dir_all_save,"/LMG_",sample_name,"_res.1.1.csv"))
}

map2(Seur_list,names(Seur_list), function(X,n){ find_LMGs(X,n) })