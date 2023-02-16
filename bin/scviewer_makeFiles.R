#!/usr/bin/env Rscript

library(scviewer)
library(Seurat)
library(gtools)
library(formattable)

# GEN.dir<-"/camp/stp/babs/working/schneid/"
# # #"/home/schneid/Documents/CAMP/"
# # #"/home/deborah/Documents/Crick/Projects/"
# #
# 
# dir_all_sub=paste0(GEN.dir,
#                    "projects/sahaie/giovanni.giangreco/",
#                    "Characterisation_of_CAF_in_HPV_cancer_scrnaseq")
# 
# dir_all=paste0(dir_all_sub,"/Data_interim_files")
# Everything.combined <-
#   read_rds(paste0(dir_all,"/SC21137_Integrated_Filtered_withMOC2_3.RDS"))

args = commandArgs(trailingOnly=TRUE)

fileRDS<-args[1]
Everything.combined <-read_rds(fileRDS)


### for visual display only - downsample by a lot 
if(nrow(Everything.combined[[]])>40000){
 N.cells<-30000
  Everything.combined = subset(Everything.combined, cells = sample(Cells(Everything.combined), N.cells))
  file_name<-'SC21137_Integrated_Filtered'
}else{
  file_name<-str_match(fileRDS,'(.*).RDS')[,2]
}

###
Everything.combined <- SetIdent(Everything.combined,
                                      value = Everything.combined@meta.data$integrated_snn_res.0.3)
DefaultAssay(Everything.combined)<-'RNA'

##############

h5_file <- paste0(file_name,".scv",sep='')
#h5_file <- paste0(dir_all_sub,"/scviewer_files/",file_name,".scv",sep='')
if(!file.exists(h5_file)){
  create_h5_scv(h5_file=h5_file)
}
npc=30
reductions <- list(
  pca     = as.data.frame(Embeddings(Everything.combined,"pca"))[,1:2],
  pca_3d  = as.data.frame(Embeddings(Everything.combined,"pca"))[,1:3],
  umap    = as.data.frame(Embeddings(Everything.combined,"umap")),
  umap_3d = as.data.frame(Embeddings(RunUMAP(Everything.combined, reduction = "pca", dims = 1:npc, 
                                              verbose = FALSE, seed.use = 42, n.components=3),"umap"))[,1:3],
  tsne    = as.data.frame(Embeddings(RunTSNE(Everything.combined, reduction = "pca", dims = 1:npc, 
                                             verbose = FALSE, seed.use = 42, dim.embed=2),"tsne"))[,1:2],
  tsne_3d = as.data.frame(Embeddings(RunTSNE(Everything.combined, reduction = "pca", dims = 1:npc, 
                                              verbose = FALSE, seed.use = 42, dim.embed=3),"tsne"))[,1:3]
)
reductions <- lapply(reductions,function(x){
  reduct <- x
  reduct <- cbind(rownames(reduct),reduct)
  if (ncol(reduct)==3) { colnames(reduct) <- c("cell_id","x","y")}
  if (ncol(reduct)==4) { colnames(reduct) <- c("cell_id","x","y","z")}
  rownames(reduct) <- NULL
  reduct
})
#lapply(reductions, head, n=3)
write_reductions(h5_file=h5_file, reductions=reductions)

### features
features_matrix <- scviewer:::guess_features_matrix(Everything.combined)
keep<-c("nCount_RNA","nFeature_RNA","percent.mt","percent_ribosomal")
write_features(h5_file=h5_file, features_matrix=features_matrix[,c(rownames(Everything.combined),keep)])
#write_features(h5_file=h5_file, features_matrix=features_matrix[,keep])


### meta data
library(dplyr)

clust.identites<- Everything.combined[[]] %>%
       as.data.frame() %>% dplyr::select(contains('integrated_snn_res.')) %>% colnames()

Everything.combined[[]] %>%
  as.data.frame() %>%
  #rownames_to_column('cell_id') %>%
  mutate(datasets_filter=str_replace(orig.ident, '_', ' ')) %>%
  dplyr::select(datasets_filter, cell_id , contains('integrated_snn_res.'), Cell.ident, Cell.ident.new, Cell.type, orig.ident, Phase, Ptprc_tdT_pos) %>% #
  group_by(datasets_filter) %>%
  mutate(N=n()) %>%
  ungroup() %>%
  mutate(datasets_filter={sprintf(fmt='%s (n=%s)', datasets_filter, comma(N)) %>% factor() %>% fct_relevel({levels(.) %>% mixedsort()})}) %>%
  mutate_at(clust.identites, function(x) x %>% fct_relevel({levels(.) %>% mixedsort()})) %>%
  mutate_at(c('Cell.type', 'Cell.ident','Cell.ident.new','orig.ident', 'Phase', 'Ptprc_tdT_pos'), as.factor) %>%
  dplyr::select(-N) -> metadata
write_metadata(h5_file=h5_file, metadata=metadata) 

### Cell clusters
# cluster_identity_sets <- scviewer:::guess_cluster_identity_sets(seurat)
# ID_keep<-grep('integrated_',names(cluster_identity_sets))
cluster_identity_sets <- list()

meta<-Everything.combined[[]]
meta<- meta %>%
  mutate_at(c('Cell.type', 'orig.ident', 'Phase', 'Cell.ident', 'Cell.ident.new','Ptprc_tdT_pos'), as.factor)

for (m in grep('integrated_|Cell.type|orig.ident|Phase|Cell.ident|Cell.ident.new|Ptprc_tdT_pos',colnames(meta))) {  #
  cluster_identity_sets[[colnames(meta)[m]]] <- list(
    var      = colnames(Everything.combined[[]])[m],
    name     = gsub('_',' ',colnames(meta)[m]),
    selected = sort(as.character(unique(meta[,m])))
  )
}
write_cluster_identity_sets(h5_file=h5_file, cluster_identity_sets=cluster_identity_sets)
