#!/usr/bin/env Rscript

library(tidyverse)
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
# Everything.combined.moc23 <-
#   read_rds(paste0(dir_all,"/SC21137_Integrated_Filtered_augmented.RDS"))

args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]

Everything.combined <-read_rds(fileRDS)

Everything.combined[[]] %>%
  as.data.frame() -> meta.df# %>%
# rownames_to_column('cell_id') -> meta.df

##

meta.df$tdT_detect<-ifelse(meta.df$percent.tdT>0,1,0)

meta.df$Ptprc_detect <-ifelse(meta.df$percent.ptprc>0,1,0)

meta.df$Ptprc_tdT_pos<-ifelse((meta.df$Ptprc_detect==1)&(meta.df$tdT_detect==1),1,0)

meta.df$Ptprc_tdT_pos<-as.factor(meta.df$Ptprc_tdT_pos)


##################################

GEN.dir<-"/camp/stp/babs/working/schneid/"
#"/home/schneid/Documents/CAMP/"
#"/home/deborah/Documents/Crick/Projects/"

dir_all_sub=paste0(GEN.dir,
                   "projects/sahaie/giovanni.giangreco/",
                   "Characterisation_of_CAF_in_HPV_cancer_scrnaseq")

dir_all=paste0(dir_all_sub,"/Data_interim_files")

Everything.combined.moc23 <-
  read_rds(paste0(dir_all,"/SC21137_Integrated_Filtered_augmented.RDS"))


Everything.combined.moc23@meta.data %>% 
  select(cell_id,integrated_snn_res.0.3,Cell.ident,Cell.ident.new) -> cell_IDed

meta.df %>% 
  mutate(cell_id=rownames(.)) %>% 
  select(cell_id,integrated_snn_res.0.3) %>% 
  rename(Newsnn_res.0.3=integrated_snn_res.0.3) -> cell_nonIDed

cell_IDs<-merge(cell_nonIDed,cell_IDed,by='cell_id')
merge(meta.df %>% mutate(cell_id=rownames(.)),
      cell_IDs %>% select(cell_id,Cell.ident,Cell.ident.new),by='cell_id', all.x = TRUE) -> meta.new

##############################

meta.new %>% mutate(Cell.ident=as.factor(Cell.ident))->meta.new
rownames(meta.new)<-meta.new$cell_id
Everything.combined<-AddMetaData(Everything.combined, metadata=meta.new)

#

saveRDS(
  Everything.combined.moc23,
  file = "SC21137_Integrated_Filtered_augmented.RDS"
  #  paste0(dir_all,"/SC21137_Integrated_Filtered_augmented.RDS")
)

