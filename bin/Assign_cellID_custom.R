#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(SeuratObject)

GEN.dir<-"/camp/stp/babs/working/schneid/"
#"/home/schneid/Documents/CAMP/"
#"/home/deborah/Documents/Crick/Projects/"

# dir_all_sub=paste0(GEN.dir,
#                    "projects/sahaie/giovanni.giangreco/",
#                    "Characterisation_of_CAF_in_HPV_cancer_scrnaseq")
# dir_all=paste0(dir_all_sub,"/Data_interim_files_2")
# Everything.combined.moc23 <-
#   read_rds(paste0(dir_all,"/SC21137_Integrated_Filtered.RDS"))
# annot.df<-read.csv(paste0(dir_all_sub,"/custom_annot_3.txt"),sep='\t')

args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]
annot<-args[2]

Everything.combined.moc23 <-read_rds(fileRDS)
Everything.combined.moc23[["percent.tdT"]] <- PercentageFeatureSet(Everything.combined.moc23, pattern = "tdTomato", assay = 'RNA')
Everything.combined.moc23[["percent.ptprc"]] <- PercentageFeatureSet(Everything.combined.moc23,
                                                                     pattern = "Ptprc", assay = 'RNA')

Everything.combined.moc23[[]] %>%
  as.data.frame() -> meta.df# %>%
# rownames_to_column('cell_id') -> meta.df

annot.df<-read.csv(annot)
annot.df$Cluster <- as.factor(annot.df$Cluster)

##

meta.df$tdT_detect<-ifelse(meta.df$percent.tdT>0,1,0)

meta.df$Ptprc_detect <-ifelse(meta.df$percent.ptprc>0,1,0)

meta.df$Ptprc_tdT_pos<-ifelse((meta.df$Ptprc_detect==1)&(meta.df$tdT_detect==1),1,0)

meta.df$Ptprc_tdT_pos<-as.factor(meta.df$Ptprc_tdT_pos)

### assign custom made cell-type identities

meta.df$Cluster<-meta.df$integrated_snn_res.0.9
meta.df<-meta.df %>% left_join(.,annot.df,by='Cluster') %>% rename(.,'Cell.ident.new'='Annot')

rownames(meta.df)<-meta.df$cell_id

Everything.combined.moc23<-AddMetaData(Everything.combined.moc23, metadata=meta.df)

saveRDS(
  Everything.combined.moc23,
  file = "SC21137_Integrated_Filtered_augmented.RDS"
  #  paste0(dir_all,"/SC21137_Integrated_Filtered_augmented.RDS")
)
