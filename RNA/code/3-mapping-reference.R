reference <- subset(seu, subset = cell_type_main == "B/Plasma")
reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference, npcs = 50)
reference <- RunUMAP(reference, dims = 1:40)
umap_new_model <- list()
umap_new_model$n_epochs <- 500
umap_new_model$alpha <-1
umap_new_model$method <- "umap"
umap_new_model$negative_sample_rate <- 5
umap_new_model$gamma <- 1
umap_new_model$approx_pow <- 0
umap_new_model$n_neighbors <- 30
umap_new_model$metric$cosine <- list()
umap_new_model$embedding <- reference[["umap"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap_new_model$a <- ab_param["a"]
umap_new_model$b <- ab_param["b"]
reference[["umap"]]@misc$model <- umap_new_model
#map_obj_ref <- function(seurat_obj, reference, resolution = 10){
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 50)
seurat_obj <- seurat_obj %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:40, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = 10)

anchors <- FindTransferAnchors(reference = reference, query = seurat_obj, dims = 1:30, 
                               reference.reduction = "pca")
seurat_obj <- MapQuery(
  anchorset = anchors,
  query = seurat_obj,
  reference = reference,
  refdata = list(
    celltypemain = "cell_type_main",
    celltypefine = "cell_type_fine"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

seurat_obj$corrected.major_cluster <- correct_cluster(seurat_obj, "predicted.celltypemain")
seurat_obj$corrected.minor_cluster <- correct_cluster(seurat_obj, "predicted.celltypefine")

correct_table <- table(seurat_obj$corrected.major_cluster,seurat_obj$corrected.minor_cluster)
seurat_obj$corrected.major_cluster <- sapply(seurat_obj$corrected.minor_cluster, function(x){
  rownames(correct_table)[which.max(correct_table[,x])]
})

#return(seurat_obj)
#}
saveRDS()
correct_cluster <- function(seurat_obj, correct_group){
  # Extract the metadata
  metadata <- seurat_obj@meta.data
  # Get the most frequent cell type for each cluster
  table_cluster <- table(seurat_obj$seurat_clusters, as.vector(seurat_obj[[correct_group]])[[1]])
  cluster_celltype <- apply(table_cluster,1,function(x) colnames(table_cluster)[which.max(x)])
  cluster_vector <- sapply(seurat_obj$seurat_clusters, function(x) cluster_celltype[x])
  names(cluster_vector) <- sub("(.*)\\..*","\\1",names(cluster_vector))
  return(cluster_vector)
}
DimPlot(seurat_obj,group.by = "corrected.minor_cluster")
reference$cell_type_main
reference$cluster <- case_when(reference$major_cluster %in% c("T/NK cells") ~ "Immune")
seurat_obj$predicted.celltypefine
seurat_obj$celltype_bped_main
storeTNK$
  DimPlot(storeTNK,group.by = "seurat_clusters",label = T)

corrected_minor_cluster_B <- seurat_obj$corrected.minor_cluster
cell_annotations_B_Plasma <- data.frame(barcode = names(corrected_minor_cluster_B),cell_type = unname(corrected_minor_cluster_B))
ro10w_filtered <- AddMetaData(
  object = ro10w_filtered,
  metadata = setNames(combined_annotations$cell_type, combined_annotations$barcode),
  col.name = "corrected_minor_cluster"
)

p1 <- VlnPlot(Myeloid_sub, features = markergenes, stack = TRUE, sort = FALSE, 
              group.by = "corrected.minor_cluster", flip = TRUE) + 
  theme(legend.position = "none")
markergenes <- c(              
  "CD19","CD79A","MS4A1", # B                
  "MZB1","IGKC","DERL3"#Plasma
)