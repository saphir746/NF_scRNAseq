#!/usr/bin/env Rscript

library(tidyverse)
library(Seurat)
library(SeuratObject)
library(intrinsicDimension)
library(plyr)
library(future)

args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]

seur_obj <-read_rds(fileRDS)

sample.name<-str_match(fileRDS,'SC21137_raw_scrubletted_(.*).RDS')[,2]

###########

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
convertHumanGeneList <- function(x){
  require("biomaRt")
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
  mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", version = 105)
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}
s.genes.mm10<-convertHumanGeneList(s.genes)
g2m.genes.mm10<-convertHumanGeneList(g2m.genes)

plan()

################################### Normalisation 

tmp.1 <- ifelse(seur_obj[["percent.tdT"]]>0,"Mouse_host","Cancer_inject")
seur_obj@meta.data["Cell.orig"] <- tmp.1
seur_obj<-NormalizeData(seur_obj)
seur_obj<-FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 3000)

################################### Cell-cycle scoring and regression
options(future.globals.maxSize = 8000 * 1024^2)

plan("multicore", workers = 4) 
seur_obj <- CellCycleScoring(seur_obj, s.features = s.genes.mm10,
                               g2m.features = g2m.genes.mm10, set.ident = FALSE)
all.genes <- rownames(seur_obj)
seur_obj  <- ScaleData(seur_obj , vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
seur_obj <- RunPCA(seur_obj, features = VariableFeatures(seur_obj))

################################### Clustering + UMAP 
int_dim <-
    intrinsicDimension::maxLikGlobalDimEst(seur_obj@reductions$pca@cell.embeddings, k = 10)
D<-ceiling(int_dim$dim.est)
seur_obj <- FindNeighbors(seur_obj, dims = 1:D)
seur_obj <- FindClusters(seur_obj, resolution = 0.5)# seq(0.3,1.1,0.2))
seur_obj <- RunUMAP(seur_obj, dims = 1:D)

################################## Save object

final.name<-paste0("SC21137_Umapped_",sample.name,".RDS")

saveRDS(
  seur_obj,
  file = final.name
  #paste0(dir_all,final.name)
)
