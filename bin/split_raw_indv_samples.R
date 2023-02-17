#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(SeuratObject)
library(intrinsicDimension)
library(plyr)


args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]

seurat_list <-read_rds(fileRDS)
#fileRDS=SC21137_raw_scrubletted.RDS

lapply(seurat_list, function(seur_obj) {
  sample.name<-seur_obj@meta.data$orig.ident %>% unique()
  new.name<-paste0("SC21137_raw_scrubletted_",sample.name,".RDS")
  saveRDS(
    seur_obj,
    file = new.name
  )
 }
)