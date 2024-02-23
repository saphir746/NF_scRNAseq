#!/usr/bin/env Rscript

library(scviewer)
library(Seurat)
library(gtools)
library(formattable)
library(dplyr)

GEN.dir<-"/nemo/stp/babs/working/schneid/projects/"
# # #"/home/schneid/Documents/CAMP/"
# # #"/home/deborah/Documents/Crick/Projects/"
# #
# 
dir_all_sub=paste0(GEN.dir,
                   "/haydaya/leticia.monin/lm510/")

dir_all=paste0(dir_all_sub,"/Data_interim")
# Seur_list <-
#   read_rds(paste0(dir_all,"/lm510_UmappedFiltered_seurat_object.RDS"))

#args = commandArgs(trailingOnly=TRUE)

#fileRDS<-args[1]
Seur_list_QCed <-
  read_rds(paste0(dir_all,"/lm510_UmappedFiltered_seurat_object.RDS"))

Seur_obj_FastMNN <-
  read_rds(paste0(dir_all,"/lm510_IntegratedFastMNN_seurat_object.RDS"))

Seur_obj_FastMNN_annot <-
  read_rds(paste0(dir_all,"/lm510_IntegratedFastMnn_cellAnnotated.RDS"))

keep<-c("nCount_RNA","nFeature_RNA","percent.mt","percent_ribosomal")

collate_feature_matrix<-function(seurat){
  cbind({
    seurat[['RNA']]$data %>% Matrix::as.matrix() %>% 
      Matrix::t()
  }, {
     seurat@meta.data[,keep] %>% Matrix::as.matrix()
  }) %>% Matrix::as.matrix()
}

Do_everyhting<-function(Everything.combined,sample_name){
### for visual display only - downsample by a lot 
  if(nrow(Everything.combined[[]])>50000){
    N.cells<-30000
    Everything.combined = subset(Everything.combined, cells = sample(Cells(Everything.combined), N.cells))
  }
  file_name<-sample_name
  ###
  Everything.combined <- SetIdent(Everything.combined,
                                  value = Everything.combined@meta.data$RNA_snn_res.0.5)
  DefaultAssay(Everything.combined)<-'RNA'

##############

#h5_file <- paste0(file_name,".scv",sep='')
  h5_file <- paste0(dir_all_sub,"/scviewer_files/",file_name,".scv",sep='')
  if(!file.exists(h5_file)){
    create_h5_scv(h5_file=h5_file)
  }
  npc=30
  reductions <- scviewer:::guess_reductions(Everything.combined)
  Reductions <- list(
    pca     = reductions$pca,
    pca_3d  = reductions$pca_3d,
    umap    = as.data.frame(Embeddings(RunUMAP(Everything.combined, reduction = "pca", dims = 1:npc,
                                               verbose = FALSE, seed.use = 42, n.components=2),"umap"))[,1:2],
    umap_3d = as.data.frame(Embeddings(RunUMAP(Everything.combined, reduction = "pca", dims = 1:npc,
                                               verbose = FALSE, seed.use = 42, n.components=3),"umap"))[,1:3],
    tsne    = as.data.frame(Embeddings(RunTSNE(Everything.combined, reduction = "pca", dims = 1:npc,
                                               verbose = FALSE, seed.use = 42, dim.embed=2),"tsne"))[,1:2],
    tsne_3d = as.data.frame(Embeddings(RunTSNE(Everything.combined, reduction = "pca", dims = 1:npc,
                                               verbose = FALSE, seed.use = 42, dim.embed=3),"tsne"))[,1:3]
  )
  Reductions <- lapply(Reductions,function(x){
    reduct <- x
    reduct <- cbind(rownames(reduct),reduct)
    if (ncol(reduct)==3) { colnames(reduct) <- c("cell_id","x","y")}
    if (ncol(reduct)==4) { colnames(reduct) <- c("cell_id","x","y","z")}
    rownames(reduct) <- NULL
    reduct
  })
#lapply(reductions, head, n=3)
  write_reductions(h5_file=h5_file, reductions=Reductions)

  ### features
  features_matrix <- collate_feature_matrix(Everything.combined)
  write_features(h5_file=h5_file, features_matrix=features_matrix)


### meta data

  clust.identites<- Everything.combined[[]] %>%
    as.data.frame() %>% dplyr::select(contains('RNA_snn_res.')) %>% colnames()
  
  Everything.combined[[]] %>%
    as.data.frame() %>%
    rownames_to_column('cell_id') %>%
    mutate(datasets_filter=str_replace_all(orig.ident, '_', ' ')) %>%
    dplyr::select(datasets_filter, cell_id , orig.ident, group , sorting_batch , Cell.type, scDblFinder.class, contains('RNA_snn_res.'), Phase) %>% #
    group_by(orig.ident) %>%
    dplyr::mutate(N=n()) %>%
    ungroup() %>%
    mutate(datasets_filter={sprintf(fmt='%s (n=%s)', datasets_filter, comma(N)) %>% factor() %>% fct_relevel({levels(.) %>% mixedsort()})}) %>%
    mutate_at(clust.identites, function(x) x %>% fct_relevel({levels(.) %>% mixedsort()})) %>%
    mutate_at(c('orig.ident', 'group', 'sorting_batch', 'Phase', 'Cell.type','scDblFinder.class'), as.factor) %>%
    dplyr::select(-N) -> metadata
  write_metadata(h5_file=h5_file, metadata=metadata) 

### Cell clusters
# cluster_identity_sets <- scviewer:::guess_cluster_identity_sets(seurat)
# ID_keep<-grep('integrated_',names(cluster_identity_sets))
  cluster_identity_sets <- list()
  
  meta<-Everything.combined[[]]
  meta<- meta %>%
    mutate_at(c('orig.ident','Phase','scDblFinder.class','Cell.type'), as.factor)
  
  for (m in grep('RNA_|orig.ident|Phase|group|sorting_batch|Cell.type|scDblFinder.class',colnames(meta))) {  #
    cluster_identity_sets[[colnames(meta)[m]]] <- list(
      var      = colnames(Everything.combined[[]])[m],
      name     = gsub('_',' ',colnames(meta)[m]),
      selected = sort(as.character(unique(meta[,m])))
    )
  }
  write_cluster_identity_sets(h5_file=h5_file, cluster_identity_sets=cluster_identity_sets)
}

map2(Seur_list_QCed,names(Seur_list_QCed), function(Obj,n){ Do_everyhting(Obj,n) })
#
Do_everyhting(Seur_obj_FastMNN,"IntegratedFastMNN")
