library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)

counts <- Read10X(data.dir = "SISA137.lung/filtered_mtx", gene.column = 1,cell.column = 1)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = 'SISA137.lung/fragments_filtered.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc.137 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#import annotations dataset

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(test.137) <- annotations

# compute nucleosome signal score per cell
test.137 <- NucleosomeSignal(object = test.137)

# compute TSS enrichment score per cell
test.137 <- TSSEnrichment(object = test.137, fast = FALSE)

# # add blacklist ratio and fraction of reads in peaks
# pbmc.137.137$pct_reads_in_peaks <- pbmc.137.137$peak_region_fragments / pbmc.137.137$passed_filters * 100
# pbmc.137.137$blacklist_ratio <- pbmc.137.137$blacklist_region_fragments / pbmc.137.137$peak_region_fragments

# DensityScatter(pbmc.137.137, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

pbmc.137$high.tss <- ifelse(pbmc.137$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(pbmc.137, group.by = 'high.tss') + NoLegend()

pbmc.137$nucleosome_group <- ifelse(pbmc.137$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc.137, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc.137,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

pbmc.137 <- subset(
  x = pbmc.137,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 40000 &
    # pct_reads_in_peaks > 15 &
    # blacklist_ratio < 0.05 &
    nucleosome_signal < 5 &
    TSS.enrichment > 3
)

#归一化、SVD分解、降维、聚类
pbmc.137 <- RunTFIDF(pbmc.137)
pbmc.137 <- FindTopFeatures(pbmc.137, min.cutoff = 'q0')
pbmc.137 <- RunSVD(pbmc.137)

pbmc.137 <- RunUMAP(object = pbmc.137, reduction = 'lsi', dims = 2:30)
pbmc.137 <- FindNeighbors(object = pbmc.137, reduction = 'lsi', dims = 2:30)
pbmc.137 <- FindClusters(object = pbmc.137, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc.137, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(pbmc.137)

# add the gene activity matrix(important) to the Seurat object as a new assay and normalize it
pbmc.137[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc.137 <- NormalizeData(
  object = pbmc.137,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc.137$nCount_RNA)
)

DefaultAssay(pbmc.137) <- 'RNA'

markerGenes  <- c(
  "CLDN5","PECAM1", #endothelial
  "RGS5", "ACTA2", #vascular
  "ISLR","CTHRC1", # mesenchymal stromal cell
  "JCHAIN", "MZB1", #B-Cell
  "AIF1", "LYZ", #metastasis-associated macrophages
  "CD1C", "CLEC10A", #DC cell
  "GFAP", "S100B", #reactive astrocyte
  "CD3D", "IL7R" #T-Cells
)

FeaturePlot(
  object = pbmc.137,
  features =markerGenes,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

#import SISA81
counts <- Read10X(data.dir = "SISA81.lung/filtered_mtx", gene.column = 1,cell.column = 1)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = 'SISA81.lung/fragments_filtered.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc.81 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

# add the gene information to the object
Annotation(pbmc.81) <- annotations

# compute nucleosome signal score per cell
pbmc.81 <- NucleosomeSignal(object = pbmc.81)

# compute TSS enrichment score per cell
pbmc.81 <- TSSEnrichment(object = pbmc.81, fast = FALSE)

# DensityScatter(pbmc.81, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

pbmc.81$high.tss <- ifelse(pbmc.81$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(pbmc.81, group.by = 'high.tss') + NoLegend()

pbmc.81$nucleosome_group <- ifelse(pbmc.81$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc.81, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc.81,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

pbmc.81 <- subset(
  x = pbmc.81,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 40000 &
    # pct_reads_in_peaks > 15 &
    # blacklist_ratio < 0.05 &
    nucleosome_signal < 5 &
    TSS.enrichment > 3
)

#归一化、SVD分解、降维、聚类
pbmc.81 <- RunTFIDF(pbmc.81)
pbmc.81 <- FindTopFeatures(pbmc.81, min.cutoff = 'q0')
pbmc.81 <- RunSVD(pbmc.81)

pbmc.81 <- RunUMAP(object = pbmc.81, reduction = 'lsi', dims = 2:30)
pbmc.81 <- FindNeighbors(object = pbmc.81, reduction = 'lsi', dims = 2:30)
pbmc.81 <- FindClusters(object = pbmc.81, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc.81, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(pbmc.81)

# add the gene activity matrix(important) to the Seurat object as a new assay and normalize it
pbmc.81[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc.81 <- NormalizeData(
  object = pbmc.81,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc.81$nCount_RNA)
)

DefaultAssay(pbmc.81) <- 'RNA'

markerGenes  <- c(
  "CLDN5","PECAM1", #endothelial
  "RGS5", "ACTA2", #vascular
  "ISLR","CTHRC1", # mesenchymal stromal cell
  "JCHAIN", "MZB1", #B-Cell
  "AIF1", "LYZ", #metastasis-associated macrophages
  "CD1C", "CLEC10A", #DC cell
  "GFAP", "S100B", #reactive astrocyte
  "CD3D", "IL7R" #T-Cells
)

FeaturePlot(
  object = pbmc.81,
  features =markerGenes,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
