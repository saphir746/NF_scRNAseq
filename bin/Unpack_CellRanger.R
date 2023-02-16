#!/usr/bin/env Rscript


library(tidyverse)
#library(ggpubr)
#library(SingleCellExperiment)
#library(scDblFinder)
library(Seurat)
library(SeuratObject)
#library(sceasy)


args = commandArgs(trailingOnly=TRUE)

dir_all_sub<-args[1]

###################################################

cellranger_quantifications <- paste0(dir_all_sub,'/output_cellRanger')

sample_meta <-
  read_csv(paste0(dir_all_sub,'/SC21137_samplesheet_ASF.csv'))

sample_meta<-sample_meta[,c("Sample Name","Sample Replicate Group")]
colnames(sample_meta)<-gsub(' ','_',colnames(sample_meta))
rep<-unlist(lapply(sample_meta$Sample_Replicate_Group, function(x) str_split(x,' ')[[1]][1]))
genotype<-unlist(lapply(sample_meta$Sample_Replicate_Group, function(x) str_split(x,' ')[[1]][3]))
rep<-gsub('Replicate','',rep)
sample_meta$Replicate<-rep
sample_meta$Genotype<-genotype
sample_meta<-sample_meta[,c("Sample_Name","Replicate","Genotype")]


###################################################

sample_path_10x <-
  dir(path=cellranger_quantifications,
      recursive = TRUE,
      pattern = "filtered_feature_bc_matrix$",
      include.dirs = TRUE
  )

sample_path_10x <-
  set_names(sample_path_10x,
            str_extract(sample_path_10x, "^[a-zA-Z0-9_]+"))

###################################################
############################### Load cellRanger stuff and filter out obvious crap

seurat_list <-
  map2(sample_path_10x, names(sample_path_10x), function(file_path, proj_name) {
    data_10x <- Read10X(data.dir = paste0(cellranger_quantifications,'/',file_path))
    # data_10x
    CreateSeuratObject(
      counts = data_10x,
      project = proj_name,
      min.cells = 2, #50,
      # include features detected in at least this many cells
      min.features = 100#400
    ) # include cells where at least this many features are detected
  })

################

seurat_list <-
  map(seurat_list, function(seur_obj) {
    m_dat <- seur_obj[[]] %>%
      as_tibble(rownames = "cell_id") %>%
      left_join(sample_meta, by = c("orig.ident" = "Sample_Name")) %>%
      dplyr::select(-orig.ident,-nCount_RNA,-nFeature_RNA) %>%
      column_to_rownames(var = "cell_id")
    seur_obj <- AddMetaData(seur_obj, metadata = m_dat)
    seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^mt-")
    seur_obj[["percent_ribosomal"]] <- PercentageFeatureSet(seur_obj, pattern = "^Rp[sl]")
    seur_obj[["percent.tdT"]] <- PercentageFeatureSet(seur_obj, pattern = "tdTomato")
   # seur_obj[["Cell.orig"]] <- ifelse(seur_obj[["percent.tdT"]]>0,'Mouse_host','Cancer_inject')
    seur_obj
  })


################

saveRDS(
  seurat_list,
  file = "SC21137_raw_seurat_object.RDS"
    #paste0(dir_all,"/SC21137_raw_seurat_object.RDS")
)
