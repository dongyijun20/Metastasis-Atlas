library(dplyr)
library(Seurat)
library(purrr)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(DoubletFinder)
pat <- "LUAD"
seu.data <- Read10X(data.dir = "45831556/count/sample_feature_bc_matrix/")
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^MT-', col.name = 'percent.mt')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPS', col.name = 'percent.rps')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPL', col.name = 'percent.rpl')
seu_raw$percent.rp <- seu_raw$percent.rps + seu_raw$percent.rpl
seu_raw$barcode <- rownames(seu_raw@meta.data)
seu_raw$barcode_pat <- paste0(rownames(seu_raw@meta.data), '_', pat)

# ### Identify doublets using scrublet
doublet_rate_tmp <- ncol(seu_raw)*8*1e-6 
#workflow
seu_raw <- NormalizeData(seu_raw)
seu_raw <- FindVariableFeatures(seu_raw)
seu_raw <- ScaleData(seu_raw)
seu_raw <- RunPCA(seu_raw)
seu_raw <- RunUMAP(seu_raw, dims = 1:40)
seu_raw <- FindNeighbors(seu_raw, dims = 1:40)
seu_raw <- FindClusters(seu_raw, resolution = 0.5)

### Identify doublets using DF
sweep.list <- paramSweep(seu_raw, PCs = 1:40, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn %>%
  arrange(desc(BCmetric))
pK <- pK[1, 2]
pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
nExp <- round(doublet_rate_tmp * dim(seu_raw@assays$RNA@layers$counts)[2])
seu_raw <- doubletFinder(seu_raw, PCs = 1:40, pK = pK, nExp = nExp)
seu_raw$doublet <- seu_raw@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)]
seu_raw$DF_score <- seu_raw@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)]

### Filtering steps 
pdf(paste0("plot_new/",pat,"_QC.pdf"))
VlnPlot(seu_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
# table(seu_raw$doublet,seu_raw$predicted_doublets)
DimPlot(seu_raw, cells.highlight = colnames(seu_raw)[seu_raw$doublet=="Doublet"], sizes.highlight = 0.1)
# DimPlot(seu_raw, cells.highlight = colnames(seu_raw)[seu_raw$predicted_doublets==T], sizes.highlight = 0.1)
dev.off()
minFeature <- 500
maxFeature <- 10000
minCount <- 1000
maxCount <- 60000
maxMT <- 10

seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature &
                nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & 
                doublet == 'Singlet')
ncol(seu_raw)
ncol(seu)

### Seurat workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
seu <- RunUMAP(seu, dims = 1:40)
seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu)

### Preliminary cell type identification using SingleR
seu_sce <- as.SingleCellExperiment(seu)

bped <- BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels

seu@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)] <- NULL
seu@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)] <- NULL
saveRDS(seu, file = paste0(pat, ".rds"))
