library(ArchR)
library(parallel)
library(pheatmap)
addArchRThreads(threads = 4)

#Make a ridge plot for each sample for the TSS enrichment scores.
pp1 <-plotGroups(  
  ArchRProj = ArchRsub.lung,  
  groupBy = "Sample",  
  colorBy = "cellColData",  
  name = "TSSEnrichment",  
  plotAs = "ridges"
)

#Make a violin plot for each sample for the TSS enrichment scores.
pp2 <- plotGroups(
  ArchRProj = ArchRsub.lung, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

#Make a ridge plot for each sample for the log10(unique nuclear fragments).
pp3 <- plotGroups(
  ArchRProj = ArchRsub.lung, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

#Make a violin plot for each sample for the log10(unique nuclear fragments).
pp4 <- plotGroups(
  ArchRProj = ArchRsub.lung, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotPDF(pp1,pp2,pp3,pp4, name = "QC-Sample-Statistics.pdf", ArchRProj = ArchRsub.lung, addDOC = FALSE, width = 4, height = 4)

saveArchRProject(ArchRProj = ArchRsub.lung, outputDirectory = "ArchRsub.lung", load = FALSE)

###画图的位置
#LSI降维，迭代了2次/3次
ArchRsub.lung<- addIterativeLSI(  
  ArchRProj = ArchRsub.lung,  
  useMatrix = "TileMatrix",  
  name = "NewIterativeLSI.2",  
  iterations = 4,  
  clusterParams = list(    
  resolution = c(0.2),    
  sampleCells = 8000,   
  n.start = 10  
  ),  
  varFeatures = 25000,  
  dimsToUse = 1:30,  
)
###画图的位置

#harmony去batch
ArchRsub.lung <- addHarmony(
  ArchRProj = ArchRsub.lung,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)
ArchRsub.lung <- addHarmony(
  ArchRProj = ArchRsub.lung,
  reducedDims = "NewIterativeLSI",
  name = "NewHarmony",
  groupBy = "Sample"
)
#聚类
ArchRsub.lung <- addClusters(
  input = ArchRsub.lung,
  reducedDims = "IterativeLSI",
  method = "Seurat", #采用的是Seurat的方法
  name = "Clusters",
  resolution = 0.8
)
ArchRsub.lung <- addClusters(
  input = ArchRsub.lung,
  reducedDims = "NewIterativeLSI",
  method = "Seurat", #采用的是Seurat的方法
  name = "NewClusters",
  resolution = 0.8
)
ArchRsub.lung <- addClusters(
  input = ArchRsub.lung,
  reducedDims = "NewIterativeLSI.2",
  method = "Seurat", #采用的是Seurat的方法
  name = "NewClusters.2",
  resolution = 0.8
)
ArchRsub.lung <- addClusters(
  input = ArchRsub.lung,
  reducedDims = "NewHarmony",
  method = "Seurat", #采用的是Seurat的方法
  name = "NewClusters_harmony",
  resolution = 0.8
)
ArchRsub.lung <- addClusters(
  input = ArchRsub.lung,
  reducedDims = "IterativeLSI",
  method = "scran", #采用的是scran的方法
  name = "ScranClusters",
  k = 15
)
cM <- confusionMatrix(paste0(ArchRsub.lung$Clusters), paste0(ArchRsub.lung$Sample))
cM <- confusionMatrix(paste0(ArchRsub.lung$NewClusters_harmony), paste0(ArchRsub.lung$Sample_true))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
       mat = as.matrix(cM), 
       color = paletteContinuous("whiteBlue"), 
       border_color = "black"
)
#UMAP
ArchRsub.lung <- addUMAP(
  ArchRProj = ArchRsub.lung, 
  reducedDims = "NewIterativeLSI.2", 
  name = "NewUMAP.2", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
ArchRsub.lung <- addUMAP(
  ArchRProj = ArchRsub.lung, 
  reducedDims = "NewIterativeLSI", 
  name = "NewUMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
p1 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Sample_true", embedding = "NewUMAP")
p2.4 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "NewClusters.2", embedding = "NewUMAP")
#t-SNE
ArchRsub.lung <- addTSNE(
  ArchRProj = ArchRsub.lung, 
  reducedDims = "IterativeLSI", 
  name = "TSNE", 
  perplexity = 30
)
p1 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
#基于batch后的重新UMAP/t-SNE
ArchRsub.lung <- addUMAP(
  ArchRProj = ArchRsub.lung, 
  reducedDims = "NewHarmony", 
  name = "NewUMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
ArchRsub.lung <- addTSNE(
  ArchRProj = ArchRsub.lung, 
  reducedDims = "Harmony", 
  name = "TSNEHarmony", 
  perplexity = 30
)
p3.1 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Sample", embedding = "NewUMAPHarmony")
p4.1 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "NewClusters", embedding = "NewUMAPHarmony")
p3 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")

#Identifying Marker Genes
markersGS <- getMarkerFeatures(
  ArchRProj = ArchRsub.lung, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "NewClusters.2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6

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
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = ArchRsub.lung, addDOC = FALSE)
#embed on UMAP
p <- plotEmbedding(
  ArchRProj = ArchRsub.lung, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "NewUMAP.2",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p <- plotBrowserTrack(
  ArchRProj = ArchRsub.lung, 
  groupBy = "NewClusters.2", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$)