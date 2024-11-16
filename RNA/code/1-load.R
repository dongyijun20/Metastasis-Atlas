library(Seurat)
library(tidyverse)
library(harmony)
library(DoubletFinder)
minFeature <- 500
maxFeature <- 10000
minCount <- 1000
maxCount <- 60000
maxMT <- 10
datasource <- c("2291206","2484697","44956007","45831556","48253637","4993157","BM1778LF","BM1778RO","BM2049332","BM6284RF","BMC377923","ThBM2441")
sub_dirs <- paste0("~/mestastasis_atlas/RNA/",datasource,"/count/sample_feature_bc_matrix")
seurat_list <- list()
for (dataset in datasource) {
  seurat_obj <- Read10X(data.dir = paste0("~/mestastasis_atlas/RNA/",dataset,"/count/sample_feature_bc_matrix"))
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = dataset)  
  seurat_list[[dataset]] <- seurat_obj
}
for(i in datasource) {
  pat <- i
  seurat_list[[i]] <- PercentageFeatureSet(seurat_list[[i]], pattern = '^MT-', col.name = 'percent.mt')
  seurat_list[[i]] <- PercentageFeatureSet(seurat_list[[i]], pattern = '^RPS', col.name = 'percent.rps')
  seurat_list[[i]] <- PercentageFeatureSet(seurat_list[[i]], pattern = '^RPL', col.name = 'percent.rpl')
  seurat_list[[i]]$percent.rp <- seurat_list[[i]]$percent.rps + seurat_list[[i]]$percent.rpl
  seurat_list[[i]]$barcode <- rownames(seurat_list[[i]]@meta.data)
  seurat_list[[i]]$barcode_pat <- paste0(rownames(seurat_list[[i]]@meta.data), '_', pat)
  doublet_rate_tmp <- ncol(seurat_list[[i]])*8*1e-6 
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]])
  seurat_list[[i]] <- ScaleData(seurat_list[[i]])
  seurat_list[[i]] <- RunPCA(seurat_list[[i]])
  seurat_list[[i]] <- RunUMAP(seurat_list[[i]], dims = 1:40)
  seurat_list[[i]] <- FindNeighbors(seurat_list[[i]], dims = 1:40)
  seurat_list[[i]] <- FindClusters(seurat_list[[i]], resolution = 0.5)
  sweep.list <- paramSweep(seurat_list[[i]], PCs = 1:40, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>%
    arrange(desc(BCmetric))
  pK <- pK[1, 2]
  pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
  nExp <- round(doublet_rate_tmp * dim(seurat_list[[i]]@assays$RNA@layers$counts)[2])
  seurat_list[[i]] <- doubletFinder(seurat_list[[i]], PCs = 1:40, pK = pK, nExp = nExp)
  seurat_list[[i]]$doublet <- seurat_list[[i]]@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)]
  seurat_list[[i]]$DF_score <- seurat_list[[i]]@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)]
  seurat_list[[i]] <- subset(seurat_list[[i]], subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature &
                  nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & 
                  doublet == 'Singlet')
  print(paste("data", i, "is finished alreadyyyyyyy"))
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]])
  seurat_list[[i]] <- ScaleData(seurat_list[[i]])
  seurat_list[[i]] <- RunPCA(seurat_list[[i]], npcs = 50)
  seurat_list[[i]] <- RunUMAP(seurat_list[[i]], dims = 1:40)
  seurat_list[[i]] <- FindNeighbors(seurat_list[[i]], dims = 1:40)
  seurat_list[[i]] <- FindClusters(seurat_list[[i]])
  saveRDS(seurat_list[[i]],file = paste0("~/mestastasis_atlas/RNA/",i,"/",i,".rds"))
  p1 <- VlnPlot(seurat_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  pdf(file = paste0("~/mestastasis_atlas/RNA/",i,"/QC_plots.pdf"), width = 10, height = 8)
  print(p1)
  dev.off()
}
combined <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
)
saveRDS(combined,"~/mestastasis_atlas/RNA/combinedV5.rds")