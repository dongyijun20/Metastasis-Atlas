library(ArchR)
library(parallel)
addArchRThreads(threads = 4)
addArchRGenome("hg38")
h5disableFileLocking()
sampleID_fileID<-paste0(c("SISA137.lung/","SISA81.lung/"),"fragments_filtered.tsv.gz")
names(sampleID_fileID) <- c("SISA137","SISA81")
ArrowFiles_sample <- createArrowFiles(
  inputFiles = sampleID_fileID,
  outputNames = names(sampleID_fileID),
  sampleNames = names(sampleID_fileID),
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  minTSS = 4,
  minFrags = 1000,
  maxFrags = 100000,
  addTileMat = TRUE, 
  addGeneScoreMat = TRUE
)
proj_name <- ArchRProject(
  ArrowFiles = ArrowFiles_sample,
  outputDirectory = "NewArchROutput1_LUAD",
  copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usag.sdo
)

#保证与signac的细胞一致，共3461个细胞(SISA81)
ArchRsub.lungnew <- subsetArchRProject(
  ArchRProj = proj_name,
  cells = cells_to_keep_true, #来源于signac的barcodes
  outputDirectory = "ArchROutput_lung_filtered_with_signac_new",
  force = TRUE
)
genescore<-getMatrixFromProject(ArchRsub.lungnew,useMatrix = "GeneScoreMatrix")
save(genescore,file="result/LUADgenescore.rdata")
h5disableFileLocking()

ArchRsub.lungnew<-addIterativeLSI(
  ArchRProj = ArchRsub.lungnew,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 8000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
)
ArchRsub.lungnew <- addHarmony(
  ArchRProj = ArchRsub.lungnew,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force =TRUE
)
h5disableFileLocking()
ArchRsub.lungnew <- addClusters(
  input = ArchRsub.lungnew,
  reducedDims = "IterativeLSI",
  method = "seurat",
  name = "Clusters",
  resolution =2,
  force =TRUE
)
h5disableFileLocking()
ArchRsub.lungnew <- addUMAP(
  ArchRProj = ArchRsub.lungnew,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force =TRUE
)
#saveArchRProject(ArchRsub.lungnew,outputDirectory ="result/all_proj_QC", dropCells = T)

p1 <- ArchR::plotEmbedding(ArchRProj = ArchRsub.lungnew, colorBy = "cellColData", name = "Sample", embedding = "UMAP", randomize = T)
p2 <- ArchR::plotEmbedding(ArchRProj = ArchRsub.lungnew, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", randomize = T)
ggAlignPlots(p1, p2, type = "h")
#plotPDF(p1,p2, name = "LUDA-Plot-UMAP-Sample-Clusters.pdf", ArchRProj = ArchRsub.lungnew, addDOC = FALSE, width = 5, height = 5)
ArchRsub.lungnew <- addUMAP(
  ArchRProj = ArchRsub.lungnew, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
p3 <- plotEmbedding(ArchRProj = ArchRsub.lungnew, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = ArchRsub.lungnew, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")

p1 <- plotGroups(
  ArchRProj = ArchRsub.lungnew, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p2 <- plotGroups(
  ArchRProj = ArchRsub.lungnew, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = ArchRsub.lungnew, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p4 <- plotGroups(
  ArchRProj = ArchRsub.lungnew, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
#plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = ArchRsub.lungnew, addDOC = FALSE, width = 4, height = 4)

#ArchRsub.lungnew<-loadArchRProject("result/all_proj_QC")

h5disableFileLocking()
markersGS <- getMarkerFeatures(
  ArchRProj = ArchRsub.lungnew, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

ArchRsub.lungnew <- addImputeWeights(ArchRsub.lungnew)


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

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = ArchRsub.lungnew, addDOC = FALSE)
p <- plotEmbedding(
  ArchRProj = ArchRsub.lungnew, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = ArchRsub.lungnew, 
        addDOC = FALSE, width = 5, height = 5
)
#Rearrange for grid plotting
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = ArchRsub.lungnew, 
        addDOC = FALSE, width = 5, height = 5)

markerGenes <- sapply(markerList, function(x) x$name[1:2])
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
ArchRsub.lungnew <- addGeneIntegrationMatrix(
  ArchRProj = ArchRsub.lungnew, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "celltype_bped_main",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
