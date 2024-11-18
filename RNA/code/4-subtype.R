library(forcats)
library(ggplot2)
library(gridExtra)
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
subtype_color <- c("T/NK" = "red","Epithelial" = "orange", "Melanocytes" = "yellow","Myeloid" = "green","Neuron" = "blue","B/Plasma" = "pink","Stromal" = "grey","Endothelial" = "black")
subtype_color <- c("T/NK" = "#a9c2cb","Epithelial" = "#c376a7", "Melanocytes" = "#cea5c7","Myeloid" = "#3ca0cf","Neuron" = "#405993","B/Plasma" = "pink","Stromal" = "grey","Endothelial" = "black")
subtype_color <- c("CD4+ T cells" = "#a9c2cb",
                   "CD8+ T cells TCF7+" = "#c4daec",
                   "CD8+ T cells TOX+" = "#b7deea",
                   "Tfh-like cells" = "#696a6c",
                   "T-cell doublets" = "#3674a2",
                   "Cycling cells" = "#76a2be",
                   "Tregs" = "#a5a9b0",
                   "NK cells" = "#c8c7e1",
                   "Epithelial" = "#c376a7", 
                   "Melanocytes" = "#cea5c7",
                   "MDM" = "#3ca0cf",
                   "cDC1" = "#76a2be",
                   "DC3" = "#83ab8e",
                   "Mast" = "#53738c",
                   "MDM FTL+" = "#63a3b8",
                   "Microglia" = "#76a2be",
                   "Monocytes" = "#4b6aa8",
                   "Neuron" = "#405993",
                   "Activated B cells" = "pink",
                   "Plasma cells" = "#d69a55",
                   "CAFs" = "grey",
                   "Endothelial" = "black"
                   #"Undetermined" = "#696a6c"
)
CellInfo <- ro10w_filtered@meta.data
P1=CellInfo %>% ggplot(aes(x=Metastasis_Type, fill=fct_rev(corrected_minor_cluster))) +scale_fill_manual(values = subtype_color)+
  geom_bar(color="black",position = "fill",width = 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), legend.text = element_text(color="black",size=13,face = "bold"),legend.title = element_text(color="black",size=13,face = "bold"),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15, face="bold"), axis.title=element_text(size=15,face="bold"), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+labs(y ="Composition (percentage of cells)", x= NULL)+ 
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)+
  theme(legend.title = element_blank())
P1g=ggplotGrob(P1)
P2=ggplot(CellInfo, aes(corrected_minor_cluster , fill=corrected_minor_cluster))+geom_bar(stat="count",colour = "black",width = 0.7)+  scale_fill_manual(values = subtype_color)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5, face="bold"), axis.title=element_text(size=15,face="bold"),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2",limits=c(256,131072),oob =  scales::squish)+coord_flip()+
  theme(legend.position = "none")
P2g=ggplotGrob(P2)
p3 <- DimPlot(ro10w_filtered, group.by = "corrected_minor_cluster",label = T,raster = FALSE) +
  scale_color_manual(values = subtype_color)

pdf("~/fengqian/metastasis_atlasnew/RNA/plot/ro10-4-subtype/proportion_subtypesample.pdf",width=15.5,height=13.5)
grid.arrange(grobs=list(P1g,P2g), widths = c(0.8,0.75),heights=c(1.09,0.8),layout_matrix = rbind(c(1, NA),c(1,2))) 
print(p3)
dev.off()
pdf("~/fengqian/metastasis_atlasnew/RNA/plot/ro10-4-subtype/proportion_subtypesample2.pdf",width=15.5,height=13.5)
grid.arrange(grobs=list(P1g,P2g), widths = c(0.8,0.75),heights=c(0.8,0.8),layout_matrix = rbind(c(1, NA),c(1,2))) 
print(p3)
dev.off()
#添加癌种信息
ro10w_filtered$Metastasis_Type <- "Other"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "BMC377923"] <- "CervicalSquamous"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "ThBM2441"] <- "RCC"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "4993157"] <- "SCLC"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "2484697"] <- "Colon.Adenocarcinoma"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "2291206"] <- "Serous.Ovarian.Cancer"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "44956007"] <- "BRCA"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "48253637"] <- "Eso.Squamous"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "45831556"] <- "LUAD"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "BM1778LF"] <- "KRAS mut LC"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "BM1778RO"] <- "KRAS mut LC"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "MBM8352T"] <- "Melanoma"
ro10w_filtered$Metastasis_Type[ro10w_filtered$orig.ident == "BM6284RF"] <- "BRCA"



P2 = ggplot(CellInfo, aes(corrected_minor_cluster, fill=corrected_minor_cluster)) +
  geom_bar(stat="count", colour="black", width=0.7) +
  geom_text(stat="count", aes(label=corrected_minor_cluster), hjust=-0.2, size=3.5, fontface="bold") +
  scale_fill_manual(values=subtype_color) +
  theme(
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(),
    panel.background=element_blank(),
    axis.line=element_line(colour="black"), 
    axis.text.y=element_blank(),
    axis.text.x=element_text(color="black", angle=45, hjust=1, size=10),
    axis.text=element_text(size=6.5, face="bold"), 
    axis.title=element_text(size=15, face="bold"),
    plot.margin=unit(c(-0.1, 0.5, 2.5, -0.5), "cm")
  ) +
  labs(y="number of cells", x=NULL) +
  scale_y_continuous(
    expand=c(0, 0), trans="log2", limits=c(256, 131072), oob=scales::squish
  ) +
  coord_flip() +
  theme(legend.position="none")

DimPlot(ro10w_filtered, group.by = "corrected_minor_cluster",label = T,raster = FALSE) +
  scale_color_manual(values = subtype_color)
CellInfo_filtered <- CellInfo %>%
  filter(!corrected_minor_cluster %in% c("Undetermined", "DC3", "cDC1"))

P1 = CellInfo %>% 
  ggplot(aes(x=Metastasis_Type, fill=fct_rev(corrected_minor_cluster))) +
  scale_fill_manual(values=subtype_color) +
  geom_bar(color="black", position="fill", width=0.7) +
  theme(
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(),
    panel.background=element_blank(), 
    legend.text=element_text(color="black", size=13, face="bold"),
    legend.title=element_blank(),
    axis.line=element_line(colour="black"), 
    axis.text.y=element_text(color="black", size=12),
    axis.text.x=element_text(color="black", angle=45, hjust=1, size=15),
    axis.text=element_text(size=15, face="bold"), 
    axis.title=element_text(size=15, face="bold"), 
    plot.margin=unit(c(0.2, 0.5, 0.2, 0.5), "cm")
  ) +
  labs(y="Composition (percentage of cells)", x=NULL) +
  scale_y_continuous(
    expand=c(0, 0), limits=c(0, 1), 
    breaks=c(0, 0.20, 0.40, 0.60, 0.80, 1), 
    labels=scales::percent
  ) +
  guides(fill=guide_legend(ncol=1)) + # 设置图例为单列
  theme(
    legend.position="right" # 图例在右侧
  )