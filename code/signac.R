library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)

counts <- Read10X(data.dir = "SISA80/filtered_mtx/", gene.column = 1,cell.column = 1)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = 'SISA80/fragments_filtered.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(pbmc) <- annotations

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# # add blacklist ratio and fraction of reads in peaks
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
# pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

pbmc <- subset(
  x = pbmc,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 40000 &
    # pct_reads_in_peaks > 15 &
    # blacklist_ratio < 0.05 &
    nucleosome_signal < 5 &
    TSS.enrichment > 3
)
pbmc

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(pbmc)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

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
  object = pbmc,
  features =markerGenes,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("../neurosurgery/metastasis/landscape/GSM5645891_seurat.rds")
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$cell_type,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
saveRDS(pbmc, file = "result/signac_BC2.rds")

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'cell_type',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot2

# replace each label with its most likely prediction
for(i in levels(pbmc)) {
  cells_to_reid <- WhichCells(pbmc, idents = i)
  newid <- names(which.max(table(pbmc$predicted.id[cells_to_reid])))
  Idents(pbmc, cells = cells_to_reid) <- newid
}
table(pbmc$predicted.id,pbmc$seurat_clusters)

pdf("plots/celltype_BC2.pdf")
DimPlot(
  object = pbmc,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
DimPlot(
  object = pbmc, group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
dev.off()

# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
head(da_peaks)

plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14+ Monocytes")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14+ Monocytes")
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
head(fc)

open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

head(closest_genes_cd4naive)

# set plotting order
levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 effector","Double negative T cell","NK dim", "NK bright", "pre-B cell",'B cell progenitor',"pDC","CD14+ Monocytes",'CD16+ Monocytes')

CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)

# motif
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

mouse_brain <- BC1
features.keep <- as.character(seqnames(granges(mouse_brain))) %in% standardChromosomes(granges(mouse_brain))
mouse_brain <- mouse_brain[features.keep, ] # if you have multiple assays you'll need to adjust this to keep features from the different assays
# add motif information
mouse_brain <- AddMotifs(
  object = mouse_brain,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

da_peaks <- FindAllMarkers(
  object = mouse_brain,
  ident.1 = 'Pvalb',
  ident.2 = 'Sst',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

# merge
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

# convert to genomic ranges
gr.500 <- makeGRangesFromDataFrame(peaks.500)
gr.1k <- makeGRangesFromDataFrame(peaks.1k)
gr.5k <- makeGRangesFromDataFrame(peaks.5k)
gr.10k <- makeGRangesFromDataFrame(peaks.10k)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.500, gr.1k, gr.5k, gr.10k))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

