#!/usr/bin/env Rscript


library(tidyverse)
#library(ggpubr)
#library(SingleCellExperiment)
#library(scDblFinder)
library(Seurat)
library(SeuratObject)
#library(sceasy)
options(Seurat.object.assay.version = 'v5')

###


args = commandArgs(trailingOnly=TRUE)

dir_all_sub<-args[1]
#"/nemo/stp/babs/working/schneid/projects/haydaya/leticia.monin/lm510/"
dir_all<-paste0(dir_all_sub,'Data_interim/')

###################################################

cellranger_quantifications <- paste0(dir_all_sub,'/output_cellRanger')

sample_meta <-
  read_csv(paste0(dir_all_sub,'/sampleSheet_design.csv'))

sample_meta<-sample_meta[,c("Sample_name","group","sorting_batch","mouse_id")]
colnames(sample_meta)<-gsub(' ','_',colnames(sample_meta))


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
      left_join(sample_meta, by = c("orig.ident" = "Sample_name")) %>%
      dplyr::select(-orig.ident,-nCount_RNA,-nFeature_RNA) %>%
      column_to_rownames(var = "cell_id")
    seur_obj <- AddMetaData(seur_obj, metadata = m_dat)
    seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^mt-")
    seur_obj[["percent_ribosomal"]] <- PercentageFeatureSet(seur_obj, pattern = "^Rp[sl]")
    seur_obj[["percent.eGFP"]] <- PercentageFeatureSet(seur_obj, pattern = "eGFP")
    seur_obj
  })


################

saveRDS(
  seurat_list,
  file = #"SC21137_raw_seurat_object.RDS"
    paste0(dir_all,"/lm510_raw_seurat_object.RDS")
)
