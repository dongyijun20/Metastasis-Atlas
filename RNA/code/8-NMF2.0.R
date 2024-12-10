# 初始化列表
Genes_nmf_w_basis <- list()

# 获取要合并的矩阵名称
k_values <- paste0("Serous.Ovarian.Cancer.k", 3:10)

# 初始化合并矩阵为第一个矩阵
merged_matrix <- geneNMF.programs[[k_values[1]]][["w"]]

# 循环合并后续矩阵
for (k in k_values[-1]) {
  # 当前矩阵
  current_matrix <- geneNMF.programs[[k]][["w"]]
  
  # 按行名对齐并合并
  merged_matrix <- cbind(merged_matrix, current_matrix[rownames(merged_matrix), , drop = FALSE])
}

# 将最终合并的矩阵赋值给 Genes_nmf_w_basis[[1]]
Genes_nmf_w_basis[[1]] <- merged_matrix
# 假设每个对象的来源信息存储在一个向量中
sources <- c("Serous.Ovarian.Cancer_program_rank3_11_nrun10.RDS","Colon.Adenocarcinoma_program_rank3_11_nrun10.RDS","BRCA_program_rank3_11_nrun10.RDS", 
             "LUAD_program_rank3_11_nrun10.RDS", "Eso.Squamous_program_rank3_11_nrun10.RDS", "SCLC_program_rank3_11_nrun10.RDS", 
             "KRAS mut LC_program_rank3_11_nrun10.RDS", "CervicalSquamous_program_rank3_11_nrun10.RDS", "Melanoma_program_rank3_11_nrun10.RDS", "RCC_program_rank3_11_nrun10.RDS")

# 动态为每个对象添加注释
for (i in seq_along(Genes_nmf_w_basis)) {
  attr(Genes_nmf_w_basis[[i]], "source") <- sources[i]
}
colnames(nmf_programs[["Serous.Ovarian.Cancer"]]) <- paste0("Serous.Ovarian.Cancer_", colnames(nmf_programs[["Serous.Ovarian.Cancer"]]))
# 更新正则表达式，适配所有名称格式，包括有空格的
updated_names <- gsub(
  "^(.*?_program)(\\d+\\.\\d+)$",        # 捕获 '名称_program' 和 '数字.数字' 的部分
  "\\1_rank3_11_nrun10.RDS.\\2",          # 插入 '_rank3_11_nrun10.RDS.'
  colnames(nmf_programs[[10]])
)
colnames(nmf_programs[[10]]) <- updated_names
# 检查结果
updated_names[1:5]


#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#filter robust program
#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) {
  
  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  
  final_filter <- NULL 
  for(i in names(nmf_sel)) {
    #a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
    a <- inter_intersect[grep(substr(i, 1, 5), substr(colnames(inter_intersect), 1, 5), invert = TRUE), 
                         grep(substr(i, 1, 5), substr(colnames(inter_intersect), 1, 5))]
    b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
    if(inter_filter==T) b <- b[b>=inter_min] # selects programs with a maximum intersection of at least 10
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)                                                      
}

nmf_programs          <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs          <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 
##names(nmf_programs) <- sources
nmf_filter_ccle       <- robust_nmf_programs(nmf_programs, intra_min = 35, intra_max = 20, inter_filter=T, inter_min = 5)  
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)

# calculate similarity between programs
nmf_intersect         <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

#pheatmap(inter_intersect, 
#         color = colorRampPalette(c("blue", "white", "red"))(50), # 颜色渐变
#         clustering_distance_rows = "euclidean", # 行的聚类距离度量
#         clustering_distance_cols = "euclidean", # 列的聚类距离度量
#         clustering_method = "complete",
#         colnames = FALSE,
#         rownames = FALSE,# 聚类方法
#         display_numbers = FALSE) # 是否显示数值


#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#generate MPs with robust program cluster
#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#nmf_programs <- store_nmf_programs
### Parameters for clustering
Min_intersect_initial <- 20    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 21    # the minimal intersection cutoff for adding a new NMF to the forming cluster 
Min_group_size        <- 3     # the minimal group size to consider for defining the first NMF program in a cluster 

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list              <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first NMF. Genes_MP always has only 50 genes and evolves during the formation of the cluster
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]  # remove selected NMF
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                 <- Genes_MP  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF 
  
  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.    
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {  
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste( strsplit(i , "[.]")[[1]][1 : which(strsplit(i , "[.]")[[1]] == "RDS")]   , collapse = "."  )
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      
      Genes_MP_temp             <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history     <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP        <- Genes_MP_temp[1:50]
    
    nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1
  
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # Remove current chosen cluster
  
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
  
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered

#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#plot heatmap of NMF
#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new]) 

annotation_colors <- c("red", "blue", "green", "purple", "orange", 
                       "yellow", "cyan", "pink", "brown", "grey")

annotation_heatmap <- ggplot(data = annotation_data, 
                             aes(x = Var1, y = 1, fill = Annotation)) + 
  geom_tile() +
  scale_fill_manual(values = annotation_colors, name = "Annotation") + 
  theme_void() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "right") +
    coord_flip()
ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))


MP_list <-  do.call(cbind, MP_list)

#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#adjust the parameter
#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

nmf_programs  <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs  <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 
intra_min <- c(10,20,30,40)
intra_max <- c(10,15,20,30)
inter_min <- c(1,3,9)
for (i in intra_min){
  for (j in intra_max){
    for (k in inter_min){
      print(paste("intra_min:", i, "；intra_max:", j, "；inter_min:", k))
      nmf_filter_ccle  <- robust_nmf_programs(nmf_programs, intra_min = i, intra_max = j, inter_filter=T, inter_min = k)
      print(length(nmf_filter_ccle))
      }
  }
}
