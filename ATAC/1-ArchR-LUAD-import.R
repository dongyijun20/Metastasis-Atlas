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

#ArchR总细胞
ArchR.lung <- ArchRProject(  
  ArrowFiles = ArrowFiles_sample,  
  outputDirectory = "ArchROutputraw.lung"
)


#保证与signac的细胞一致，共3461个细胞(SISA81)
ArchRsub.lung<- subsetArchRProject(
  ArchRProj = ArchR.lung,
  cells = cells_to_keep_true, #来源于signac的barcodes
  outputDirectory = "ArchROutput_lung_filtered_with_signac"
)
#双细胞目前未成功
#doubScores <- addDoubletScores(  
  #input = ArrowFiles_sample,  
  #k = 10, #Refers to how many cells near a "pseudo-doublet" to count.  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.  LSIMethod = 1
#)
