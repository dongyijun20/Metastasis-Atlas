library(GeneNMF)
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
seu <- FindVariableFeatures(seu, nfeatures = 1000)
seu <- runNMF(seu, k = ndim, assay="SCT")
seu <- RunUMAP(seu, reduction = "NMF", dims=1:ndim, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")
DimPlot(seu, reduction = "NMF_UMAP", group.by = "Metastasis_Type", label=T) + theme(aspect.ratio = 1,
                                                                                axis.text = element_blank(),
                                                                                axis.title = element_blank(),
                                                                                axis.ticks = element_blank()) + ggtitle("NMF UMAP") + NoLegend()

seu.list <- SplitObject(seu, split.by = "Metastasis_Type")

geneNMF.programs <- multiNMF(seu.list, assay="RNA", slot="data", k=4:9, nfeatures = 1000)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=8,
                                        weight.explained = 0.7,
                                        max.genes=100)
ph <- plotMetaPrograms(geneNMF.metaprograms)
geneNMF.metaprograms$metaprograms.metrics
lapply(geneNMF.metaprograms$metaprograms.genes, head)
library(msigdbr)
library(fgsea)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
       runGSEA(program, universe=rownames(seu), category = "H")
   })
head(top_p$MP3)
head(top_p$MP1)
gsea_result <- top_p$MP10  # 选择其中一个结果（例如 MP3）
gsea_result$logP <- -log10(gsea_result$pval)
gsea_result <- gsea_result[order(gsea_result$logP, decreasing = TRUE), ]
n_top <- 10
gsea_top <- gsea_result[1:n_top, ]
barplot <- ggplot(gsea_top, aes(x = reorder(pathway, logP), y = logP, fill = logP)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Top Enriched Pathways for MP10",
    x = "Pathway",
    y = "-log10(p-value)"
  ) +
  scale_fill_gradient(low = "steelblue", high = "firebrick") +
  theme_minimal()
print(barplot)
dev.off()
mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(seu, features = mp.genes, assay="RNA", ncores=4, name = "")
VlnPlot(seu, features=names(mp.genes), group.by = "Metastasis_Type",
        pt.size = 0, ncol=5)

matrix <- seu@meta.data[,names(mp.genes)]

#dimred <- scale(matrix)
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
seu@reductions[["MPsignatures"]] <- new("DimReduc",
                                        cell.embeddings = dimred,
                                        assay.used = "RNA",
                                        key = "MP_",
                                        global = FALSE)
seu <- RunUMAP(seu, reduction="MPsignatures", dims=1:length(seu@reductions[["MPsignatures"]]),
                               metric = "euclidean", reduction.name = "umap_MP")
FeaturePlot(seu, features = names(mp.genes), reduction = "umap_MP", ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())

annotation_col$Metastasis_Type <- "other"
annotation_col$Metastasis_Type[
  grepl("^Melanoma", rownames(annotation_col))
] <- "Melanoma"
annotation_col$Metastasis_Type[
  grepl("^CervicalSquamous", rownames(annotation_col))
] <- "CervicalSquamous"
annotation_col$Metastasis_Type[
  grepl("^RCC", rownames(annotation_col))
] <- "RCC"
annotation_col$Metastasis_Type[
  grepl("^SCLC", rownames(annotation_col))
] <- "SCLC"
annotation_col$Metastasis_Type[
  grepl("^Colon.Adenocarcinoma", rownames(annotation_col))
] <- "Colon.Adenocarcinoma"
annotation_col$Metastasis_Type[
  grepl("^Serous.Ovarian.Cancer", rownames(annotation_col))
] <- "Serous.Ovarian.Cancer"
annotation_col$Metastasis_Type[
  grepl("^BRCA", rownames(annotation_col))
] <- "BRCA"
annotation_col$Metastasis_Type[
  grepl("^Eso.Squamous", rownames(annotation_col))
] <- "Eso.Squamous"
annotation_col$Metastasis_Type[
  grepl("^LUAD", rownames(annotation_col))
] <- "LUAD"
annotation_col$Metastasis_Type[
  grepl("^KRAS mut LC", rownames(annotation_col))
] <- "KRAS mut LC"
annotation_colors <- list(
  Metastasis_Type = c(
    Melanoma = "skyblue", 
    CervicalSquamous = "tomato",
    RCC = "#ece399",
    SCLC = "#408444",
    Colon.Adenocarcinoma = "#6c408e",
    Serous.Ovarian.Cancer = "#2d3462",
    BRCA = "#bc9a7f",
    Eso.Squamous = "#efd2c9",
    LUAD = "#4b6aa8",
    `KRAS mut LC` = "#b7deea"
  )
)
color_database <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                    '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
                    '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
                    '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
                    '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
                    '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
                    '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
                    '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
                    '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
                    '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
                    '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
                    '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')
ph <- pheatmap(J,
               scale = scale,
               color = palette,
               main = main,
               cluster_rows = tree,
               cluster_cols = tree,
               cutree_rows = nMP,
               cutree_cols = nMP,
               gaps_row = gaps,
               gaps_col = gaps,
               annotation_col = annotation_col,
               annotation_row = annotation_col,
               annotation_colors = annotation_colors,
               annotation_names_col = FALSE,
               annotation_names_row = FALSE,
               show_rownames = FALSE,
               show_colnames = FALSE
)


###KEGG###
module_1_genes <- geneNMF.metaprograms$metaprograms.genes$MP1
entrez_genes <- bitr(module_1_genes, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
kegg_result <- enrichKEGG(gene = entrez_genes$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.05)
kegg_result <- setReadable(kegg_result,
                           OrgDb= org.Hs.eg.db,
                           keyType= "ENTREZID")
kegg_result_MP6<- kegg_result
view(kegg_result_MP4@result)
pdf("~/fengqian/metastasis_atlasnew/RNA/plot/NMF/KEGGsub.pdf",width=10,height=8)
barplot(
  kegg_result_MP10, 
  x= "Count", #or "GeneRatio"
  color= "pvalue", #or "p.adjust", "qvalue"
  showCategory= 15, #显示pathway的数量
  font.size = 12, #字号
  title = "KEGG enrichment barplot MP10", #标题
  label_format = 30 #pathway标签长度超过30个字符串换行
)


###GO###
module_1_genes <- geneNMF.metaprograms$metaprograms.genes$MP10
entrez_genes <- bitr(module_1_genes, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
ego <- enrichGO(
  gene          = entrez_genes$ENTREZID,          # 基因列表
  OrgDb         = org.Hs.eg.db,       # 物种数据库
  keyType       = "ENTREZID",         # 基因ID类型（可以是 "SYMBOL", "ENSEMBL", 等）
  ont           = "BP",               # GO类别: BP（生物过程），MF（分子功能），CC（细胞组分）
  pAdjustMethod = "BH",               # p值校正方法
  pvalueCutoff  = 0.05,               # p值阈值
  qvalueCutoff  = 0.2                 # q值阈值
)
barplot(ego, showCategory = 10, title = "GO Enrichment Analysis MP10")

pdf("~/fengqian/metastasis_atlasnew/RNA/plot/NMF/KEGGsub.pdf",width=10,height=8)
barplot(ego, showCategory = 10, title = "GO Enrichment Analysis MP1")
dev.off()


markergenes <- c("MT1X","MT1E","MT2A","MT1F","S100A8","EGR1","DNAJB1","PLOD2")
p1 <- VlnPlot(seu,markergenes,group.by = "Metastasis_Type",stack=TRUE,sort=TRUE)+  theme(legend.position = "none")
pdf("~/test2.pdf",width = 15,height = 15)
FeaturePlot(seu,features = markergenes) 
DimPlot(seu, group.by = "Metastasis_Type",label = T,raster = FALSE)
p1
p2
dev.off()


###multiNMF###
geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05)