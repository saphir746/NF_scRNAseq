

library(Seurat)
options(Seurat.object.assay.version = "v5")
#library(BPCells)
library(ggplot2)
library(SeuratObject)
#library(biomaRt)
library(tidyverse)
library(Seurat)
library(SeuratObject)
library(reticulate)
library(data.table)
library(intrinsicDimension)
library(SeuratWrappers)
library(magrittr)

options(future.globals.maxSize = 3e9)

scviPyEnv="/nemo/stp/babs/working/schneid/conda/envs/scviPy"
use_condaenv(scviPyEnv)

GEN.dir<-"/nemo/stp/babs/working/schneid/projects/"
#"/home/schneid/Documents/nemo/"
#"/home/deborah/Documents/Crick/Projects/"

dir_all_sub=paste0(GEN.dir,
                   "/haydaya/leticia.monin/lm510/")

dir_all=paste0(dir_all_sub,"/Data_interim")
dir_all_save=paste0(dir_all,"/LMG_ind")

Seur_list<-
  read_rds(paste0(dir_all,"/lm510_UmappedFiltered_seurat_object.RDS"))


##### rename cells


Seur_list <- map(Seur_list, function(Seur){
  new.names <- Seur[[]] %>% rownames_to_column(.,var='Cell_ID') %>%
    mutate(Cell_IDnew = paste0(Cell_ID,'_',gsub('_','',orig.ident))) %>%
    select(Cell_ID,Cell_IDnew)
  Seur<-RenameCells(Seur,old.names=new.names$Cell_ID,new.names=new.names$Cell_IDnew)
  Seur
})


##### concat ind sample data into "layers"

data.list <- lapply(Seur_list, function(Seur) Seur[["RNA"]]$counts)

metadata.true <- lapply(Seur_list, function(Seur) Seur@meta.data %>% 
                     rownames_to_column(.,var='Cell_ID') %>%
                     dplyr::select(-contains('RNA_snn_res'),-seurat_clusters,-contains('nFeature_'))) %>% rbindlist(.) %>%
  #distinct(Cell_ID, .keep_all = TRUE) %>% 
  column_to_rownames(.,var='Cell_ID')

names(data.list) <- names(Seur_list)
merged.object <- CreateSeuratObject(counts = data.list, meta.data = metadata.true)

# 
# merged.object[["RNA"]] <- split(merged.object[["RNA"]], f = merged.object$orig.ident)

merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object, verbose = FALSE)
merged.object <- ScaleData(merged.object)
merged.object <- RunPCA(merged.object)

############################

do_the_things<-function(seur_obj,reduc_name){
  seur_obj<-NormalizeData(seur_obj)
  seur_obj<-FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 3000)
  # seur_obj <- CellCycleScoring(seur_obj, s.features = s.genes.mm10,
  #                              g2m.features = g2m.genes.mm10, set.ident = FALSE)
  all.genes <- rownames(seur_obj)
  seur_obj <- ScaleData(seur_obj, features = all.genes)
  seur_obj  <- RunPCA(seur_obj , features = VariableFeatures(object = seur_obj ))
  int_dim <- intrinsicDimension::maxLikGlobalDimEst(seur_obj@reductions$pca@cell.embeddings, k = 15)
  D<-ceiling(int_dim$dim.est)
  clustPrefix<- reduc_name %>% gsub("integrated.","",.) %>% paste0(.,"_snn_")
  seur_obj <- FindNeighbors(seur_obj, dims = 1:D)
  seur_obj <- FindClusters(seur_obj, resolution = seq(0.3,1.1,0.2),
                           method = "igraph", algorithm = 4,  verbose=TRUE)
  seur_obj <- RunUMAP(seur_obj, dims = 1:30, reduction = "pca", reduction.name = reduc_name)
  return(seur_obj)
}

############################


Everything.combined <- IntegrateLayers(
  object = merged.object, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

Everything.combined[["RNA"]]  <- JoinLayers(Everything.combined[['RNA']] )
Everything.combined %<>% do_the_things(.,"integrated.rpca")

saveRDS(
  Everything.combined ,
  file = #"SC21137_Umapped_Filteredseurat_object.RDS"
    paste0(dir_all,"/lm510_IntegratedRPCA_seurat_object.RDS")
)
########
Everything.combined <- IntegrateLayers(
  object = merged.object, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE)
  # object = merged.object, method = scVIIntegration,
  # new.reduction = "integrated.scvi", orig.reduction = "pca",
  #  conda_env = scviPyEnv, verbose = FALSE)

Everything.combined[["RNA"]]  <- JoinLayers(Everything.combined[['RNA']] )
Everything.combined %<>% do_the_things(.,"integrated.mnn")

saveRDS(
  Everything.combined ,
  file = #"SC21137_Umapped_Filteredseurat_object.RDS"
    paste0(dir_all,"/lm510_IntegratedFastMNN_seurat_object.RDS")
)

##########

Everything.combined <- IntegrateLayers(
object = merged.object, method = scVIIntegration,
new.reduction = "integrated.scvi", orig.reduction = "pca",
 conda_env = scviPyEnv, verbose = FALSE)

Everything.combined[["RNA"]]  <- JoinLayers(Everything.combined[['RNA']] )
Everything.combined %<>% do_the_things(.,"integrated.scvi")



saveRDS(
  Everything.combined ,
  file = #"SC21137_Umapped_Filteredseurat_object.RDS"
    paste0(dir_all,"/lm510_Integratedscvi_seurat_object.RDS")
)
