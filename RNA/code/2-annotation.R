markergenes <- c(
  "CD3D","CD3E","CD3G",#T  
  "KLRD1","GNLY","NKG7",#NK                 
  "CD19","CD79A","MS4A1", # B                 
  "EPCAM","KRT8","KRT19",#Epithelial                 
  "VWF", "PECAM1", "ENG",#Endothelial                 
  "COL1A1","COL1A2","DCN",#Stromal
  "CD68","CD163","FCGR3A","MS4A7",#Myeloid       
  "MZB1","IGKC","DERL3",#Plasma
  "TPSAB1","TPSB2","CPA3","MS4A2" #Mast
)
markergenes <- c(  "CD68","CD163","FCGR3A","MS4A7"#Myeloid
                   )
markergenes <- c("RGS5","CSPG4","ABCC9")  #pericytes
p2 = DotPlot(seu,features = markergenes, dot.scale = 8) + RotatedAxis()
p2
p1 <- VlnPlot(seu,markergenes,stack=TRUE,sort=TRUE)+  theme(legend.position = "none")
p1
seu$celltype_main <- recode(seu@meta.data$seurat_clusters,
                            "0" = "T/NK",
                            "1" = "Epithelial",
                            "2" = "Epithelial",
                            "3" = "Epithelial",
                            "4" = "Melanocytes",
                            "5" = "Myeloid",
                            "6" = "NK",
                            "7" = "Epithelial",
                            "8" = "Epithelial",
                            "9" = "T/NK",
                            "10" = "T/NK",
                            "11" = "T/NK",
                            "12" = "T/NK",
                            "13" = "T/NK",
                            "14" = "Melanocytes",
                            "15" = "Epithelial",
                            "16" = "Astrocyte",
                            "17" = "Myeloid",
                            "18" = "B/Plasma",
                            "19" = "T/NK",
                            "20" = "T/NK",
                            "21" = "B/Plasma",
                            "22" = "B/Plasma",
                            "23" = "Epithelial",
                            "24" = "Epithelial",
                            "25" = "Myeloid", #都不明显
                            "26" = "Epithelial",
                            "27" = "Epithelial",
                            "28" = "Stromal", #不明显
                            "29" = "B/Plasma",
                            "30" = "Stromal",
                            "31" = "Epithelial",
                            "32" = "B/Plasma",
                            "33" = "B/Plasma",
                            "34" = "Endothelial",
                            "35" = "Mast",
                            "36" = "Myeloid",
                            "37" = "B/Plasma",
                            "38" = "Epithelial",
                            "39" = "B/Plasma",
                            "40" = "B/Plasma",
                            "41" = "B/Plasma",
                            "42" = "Endothelial"#不明显
)
DimPlot(seu, reduction = "umap",group.by = "celltype_main", label = T,raster=FALSE)
markergenes <- c("CD8A","FOXP3","IL7R","TOX","HAVCR2") #supplementary
markergenes <- c("CD4","IL7R","SELL","CCR7","TCF7","LEF1","CD8A","TOX","LAG3","HAVCR2","PDCD1",
                 "ENTPD1","TIGIT","GZMA","CCL4","CCL5","NCAM1","GNLY","KLRF1","CD28","CTLA4",
                 "BATF","TOX2","FOXP3","IL2RA")
markergenes <- c("KLRD1","GNLY","NKG7","CD56","CD16")#NK  
markergenes <- c("CLEC9A","XCR1","CLNK","SLC38A1","RTN1","SIRPA","IDO1","LAMP3","CD200","MS4A2","KIT","MRC1",
                 "MERTK","FCGR1A","CD38","CD163L1","SELENOP","F13A1","DAB2","SIGLEC1","FTL","FTH1","NAV3",
                 "P2RY12","VCAN","FCN1","LYZ")
markergenes <- c("ABCA1","ACSL6","AGT","ANLN","APOE") #Astrocyte
p2 = DotPlot(seu, features = markergenes, dot.scale = 8) + RotatedAxis()
p2
p1 <- VlnPlot(seu,markergenes,stack=TRUE,sort=TRUE)+  theme(legend.position = "none")
p1

FeaturePlot(storeTNK,features = markergenes) 
storeTNK$corrected.minor_cluster