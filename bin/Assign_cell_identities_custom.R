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

Everything.combined.moc23 <-read_rds(fileRDS)
Everything.combined.moc23[["percent.tdT"]] <- PercentageFeatureSet(Everything.combined.moc23, pattern = "tdTomato", assay = 'RNA')
Everything.combined.moc23[["percent.ptprc"]] <- PercentageFeatureSet(Everything.combined.moc23,
                                                                     pattern = "Ptprc", assay = 'RNA')

Everything.combined.moc23[[]] %>%
  as.data.frame() -> meta.df# %>%
 # rownames_to_column('cell_id') -> meta.df

##

meta.df$tdT_detect<-ifelse(meta.df$percent.tdT>0,1,0)

meta.df$Ptprc_detect <-ifelse(meta.df$percent.ptprc>0,1,0)

meta.df$Ptprc_tdT_pos<-ifelse((meta.df$Ptprc_detect==1)&(meta.df$tdT_detect==1),1,0)

meta.df$Ptprc_tdT_pos<-as.factor(meta.df$Ptprc_tdT_pos)

### assign custom made cell-type identities

meta.df<-meta.df %>%
  mutate(Cell.ident = recode(integrated_snn_res.0.3,"0"="Cancer_cells",
                             "1"="Cancer_cells",
                             "2"="Immune_cells",
                             "3"="CAFs",
                             "4"="Cancer_cells",
                             "5"="Stromal_cells",
                             "6"="Immune_cells",
                             "7"="Immune_cells",
                             "8"="CAFs",
                             "9"="Stromal_cells",
                             "10"="Immune_cells",
                             "11"="Stromal_cells",
                             "12"="Immune_cells",
                             "13"="Stromal_cells",
                             "14"="Immune_cells",
                             "15"="Stromal_cells",
                             "16"="Stromal_cells"))
meta.df<-meta.df %>%
  mutate(Cell.ident = ifelse((Cell.ident=="Cancer_cells")&(Ptprc_detect==0)&(tdT_detect==0),"Cancer_cells",
                             ifelse((Cell.ident=="Immune_cells")|(percent.ptprc>0),"Immune_cells",
                                    ifelse(Cell.ident=="CAFs","CAFs","Stromal_cells")
                                    )
                             ))

############### 12/2022 - custom new cell identities ##########

# 0,1,4 - Cancer cells
# 2,7,12 - Monocytes / Macrophages / Dendritic cells
# 3 - myoCAFs
# 5,11,13 - Skin epithelial cells
# 6 - T cells / NK cells
# 8 - iCAFs
# 9 - Endothelial cells
# 10 - RBC / Myeloid cells
# 14 - Basophil / Mast cells
# 15 - Perycites
# 16 - Myoepithelial cells
### based on identities of integrated_snn_res.0.3


meta.df<-meta.df %>%
  mutate(Cell.ident.new = recode(integrated_snn_res.0.3,"0"="Cancer cells",
                                 "1"="Cancer cells",
                                 "2"="Mon Mac Den cells",
                                 "3"="myoCAFs",
                                 "4"="Cancer cells",
                                 "5"="Skin epithelial cells",
                                 "6"="T cells / NK cells",
                                 "7"="Mon Mac Den cells",
                                 "8"="iCAFs",
                                 "9"="Endothelial cells",
                                 "10"="RBC / Myeloid cells",
                                 "11"="Skin epithelial cells",
                                 "12"="Mon Mac Den cells",
                                 "13"="Skin epithelial cells",
                                 "14"="Basophil / Mast cells",
                                 "15"="Perycites",
                                 "16"="Myoepithelial cells"))
rownames(meta.df)<-meta.df$cell_id

Everything.combined.moc23<-AddMetaData(Everything.combined.moc23, metadata=meta.df)

saveRDS(
  Everything.combined.moc23,
  file = #"SC21137_Integrated_Filtered_augmented.RDS"
  paste0(dir_all,"/SC21137_Integrated_Filtered_augmented.RDS")
)