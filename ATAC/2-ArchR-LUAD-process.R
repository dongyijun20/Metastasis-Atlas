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

#LSI降维，迭代了2次
ArchRsub.lung<- addIterativeLSI(  
  ArchRProj = ArchRsub.lung,  
  useMatrix = "TileMatrix",  
  name = "IterativeLSI.3times",  
  iterations = 3,  
  clusterParams = list(    
  resolution = c(0.2),    
  sampleCells = 8000,   
  n.start = 10  
  ),  
  varFeatures = 25000,  
  dimsToUse = 1:30,  
)
#harmony去batch
ArchRsub.lung <- addHarmony(
  ArchRProj = ArchRsub.lung,
  reducedDims = "IterativeLSI",
  name = "Harmony",
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
  reducedDims = "IterativeLSI",
  method = "scran", #采用的是scran的方法
  name = "ScranClusters",
  k = 15
)
cM <- confusionMatrix(paste0(ArchRsub.lung$Clusters), paste0(ArchRsub.lung$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
       mat = as.matrix(cM), 
       color = paletteContinuous("whiteBlue"), 
       border_color = "black"
)
#UMAP
ArchRsub.lung <- addUMAP(
  ArchRProj = ArchRsub.lung, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
p1 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
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
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
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
p3 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
p4 <- plotEmbedding(ArchRProj = ArchRsub.lung, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
