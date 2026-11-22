

################################################################################################################################################################
#---------------------------------------------------------part1: prepare data---------------------------------------------------------------#
################################################################################################################################################################

 
library(Seurat)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dbplyr)
library(VennDiagram)
library(CellChat)
library(nichenetr)
library(tibble)
library(tidyr )
library(patchwork)
library(edgeR)
library(RColorBrewer )
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]

setwd('~/2021_NRBC_chlyu/')
NRBC_altas_seu=readRDS('~/2021_NRBC_chlyu/primitive_NRBC_project/all_NRBC_altas_seu.rds')


# modified the function netAnalysis_signalingRole_heatmap 
computeCentralityLocal <- function(net) {
  centr <- vector("list")
  G <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  centr$outdeg_unweighted <- rowSums(net > 0)
  centr$indeg_unweighted <- colSums(net > 0)
  centr$outdeg <- igraph::strength(G, mode="out")
  centr$indeg <- igraph::strength(G, mode="in")
  centr$hub <- igraph::hub_score(G)$vector
  centr$authority <- igraph::authority_score(G)$vector # A node has high authority when it is linked by many other nodes that are linking many other nodes.
  centr$eigen <- igraph::eigen_centrality(G)$vector # A measure of influence in the network that takes into account second-order connections
  centr$page_rank <- igraph::page_rank(G)$vector
  igraph::E(G)$weight <- 1/igraph::E(G)$weight
  centr$betweenness <- igraph::betweenness(G)
  #centr$flowbet <- try(sna::flowbet(net)) # a measure of its role as a gatekeeper for the flow of communication between any two cells; the total maximum flow (aggregated across all pairs of third parties) mediated by v.
  #centr$info <- try(sna::infocent(net)) # actors with higher information centrality are predicted to have greater control over the flow of information within a network; highly information-central individuals tend to have a large number of short paths to many others within the social structure.
  centr$flowbet <- tryCatch({
    sna::flowbet(net)
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  centr$info <- tryCatch({
    sna::infocent(net, diag = T, rescale = T, cmode = "lower")
    # sna::infocent(net, diag = T, rescale = T, cmode = "weak")
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net), ncol = 1))
  })
  return(centr)
}
netAnalysis_signalingRole_heatmap2 <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                               color.use = NULL, color.heatmap = "BuGn",title = NULL, width = 10, height = 8,
                                               font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat=mat[,colSums(mat)>0]
  mat=mat[rowSums(mat)>0,]
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA
  
  
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = c('#FFFFFF',grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 2, name = color.heatmap)))(100))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  
  #row_annotation_df=data.frame(row.names = slot(object, slot.name)$pathways,annotation=object@DB$interaction[match(slot(object, slot.name)$pathways,object@DB$interaction$pathway_name),'annotation'])
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)
  
  
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  mat[is.na(mat)]=0
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.cols,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = legend.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}


netAnalysis_signalingRole_heatmap3 <- function(object, signaling = NULL, pattern = c("outgoing", "incoming","all"), slot.name = "netP",
                                               color.use = NULL, color.heatmap = "BuGn",title = NULL, width = 10, height = 8,
                                               font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE){
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]]$outdeg
    incoming[,i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  } else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  } else if (pattern == "all") {
    mat <- t(outgoing+ incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  } else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ",title)
  }
  
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat=mat[,colSums(mat)>0]
  # mat=mat[rowSums(mat)>0,]
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  mat[mat == 0] <- NA
  
  
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = c('#FFFFFF',grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 2, name = color.heatmap)))(100))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  
  #row_annotation_df=data.frame(row.names = slot(object, slot.name)$pathways,annotation=object@DB$interaction[match(slot(object, slot.name)$pathways,object@DB$interaction$pathway_name),'annotation'])
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), show_annotation_name = FALSE)
  
  
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  } else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), round(max(mat, na.rm = T), digits = 1))
  }
  mat[is.na(mat)]=0
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative strength",
                bottom_annotation = col_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.cols,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = legend.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}


#-------------------------------- read data----------------------------------#

YS_subcelltype_ccd=readRDS('NRBC_YS_altas/YS_subcelltype_cellchat.rds')

object.list=list(YS=YS_subcelltype_ccd)
levels(YS_subcelltype_ccd@idents)


stage_type_list=list()
stage_type_list[['YS']][['sources.use']]=c(1:2,7:length(levels(YS_subcelltype_ccd@idents))) # no Ery
stage_type_list[['YS']][['targets.use']]=3:6    # Ery celltype


NRBC_LR_df=data.frame()
Other2Ery_df=data.frame()
Ery2Other_df=data.frame()
Ery2Ery_df=data.frame()


for ( stage in names(object.list)){
  tmp_Other2Ery_df=netVisual_bubble(object.list[[stage]],sources.use = stage_type_list[[stage]][['sources.use']],targets.use =stage_type_list[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Other2Ery_df$stage=stage;rownames(tmp_Other2Ery_df)=NULL
  tmp_Other2Ery_df=tmp_Other2Ery_df[!is.na(tmp_Other2Ery_df$prob),]
  
  Other2Ery_df=rbind(Other2Ery_df,tmp_Other2Ery_df)
  
  tmp_Ery2Other_df=netVisual_bubble(object.list[[stage]],sources.use = stage_type_list[[stage]][['targets.use']],targets.use =stage_type_list[[stage]][['sources.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Other_df$stage=stage;rownames(tmp_Ery2Other_df)=NULL
  Ery2Other_df=rbind(Ery2Other_df,tmp_Ery2Other_df)
  Ery2Other_df=Ery2Other_df[!is.na(Ery2Other_df$prob),]
  
  
  tmp_Ery2Ery_df=netVisual_bubble(object.list[[stage]],sources.use = stage_type_list[[stage]][['targets.use']],targets.use =stage_type_list[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Ery_df$stage=stage;rownames(tmp_Ery2Ery_df)=NULL
  Ery2Ery_df=rbind(Ery2Ery_df,tmp_Ery2Ery_df)
  Ery2Ery_df=Ery2Ery_df[!is.na(Ery2Ery_df$prob),]
  
}

Other2Ery_df$target_type='Other2Ery'
Ery2Other_df$target_type='Ery2Other'
Ery2Ery_df$target_type='Ery2Ery'

YS_NRBC_altas_LR_df=rbind(rbind(Other2Ery_df,Ery2Other_df),Ery2Ery_df)
YS_NRBC_altas_LR_df=YS_NRBC_altas_LR_df[!is.na(YS_NRBC_altas_LR_df$prob),]
YS_NRBC_altas_LR_df=YS_NRBC_altas_LR_df[-grep(pattern = 'ERYTHROID',YS_NRBC_altas_LR_df$source.target)]
YS_NRBC_altas_LR_df=YS_NRBC_altas_LR_df[!is.na(YS_NRBC_altas_LR_df$pathway_name),]
YS_NRBC_altas_LR_df=YS_NRBC_altas_LR_df[!is.na(YS_NRBC_altas_LR_df$interaction_name),]

#---------------------在更细的细胞亚类水平subcelltype-----------------------------------#

#-------------------------------- read data----------------------------------#
YS_allsubcelltype_ccd=readRDS('NRBC_YS_altas/YS_all_subcelltype_cellchat.rds')
levels(YS_allsubcelltype_ccd@idents)

#------------------------------------重新统计 source:target, Other celltype: NRBC celltye--------------------------------------------#
# 准备，构建所有stage list 
object.list2 <- list(YS = YS_allsubcelltype_ccd)
stage_type_list2=list()

# 未考虑YS NRBC 在FL、FBM中接收到的细胞互作
stage_type_list2[['YS']][['sources.use']]=c(5:length(levels(YS_allsubcelltype_ccd@idents))) # no Ery
stage_type_list2[['YS']][['targets.use']]=1:4  # Ery celltype

rm(YS_allsubcelltype_ccd);gc()


NRBC_LR_df2=data.frame()
Other2Ery_df2=data.frame()
Ery2Other_df2=data.frame()
Ery2Ery_df2=data.frame()


for ( stage in names(object.list2)){
  tmp_Other2Ery_df=netVisual_bubble(object.list2[[stage]],sources.use = stage_type_list2[[stage]][['sources.use']],targets.use =stage_type_list2[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Other2Ery_df$stage=stage;rownames(tmp_Other2Ery_df)=NULL
  tmp_Other2Ery_df=tmp_Other2Ery_df[!is.na(tmp_Other2Ery_df$prob),]
  
  Other2Ery_df2=rbind(Other2Ery_df2,tmp_Other2Ery_df)
  
  tmp_Ery2Other_df=netVisual_bubble(object.list2[[stage]],sources.use = stage_type_list2[[stage]][['targets.use']],targets.use =stage_type_list2[[stage]][['sources.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Other_df$stage=stage;rownames(tmp_Ery2Other_df)=NULL
  Ery2Other_df2=rbind(Ery2Other_df2,tmp_Ery2Other_df)
  Ery2Other_df2=Ery2Other_df2[!is.na(Ery2Other_df2$prob),]
  
  
  tmp_Ery2Ery_df=netVisual_bubble(object.list2[[stage]],sources.use = stage_type_list2[[stage]][['targets.use']],targets.use =stage_type_list2[[stage]][['targets.use']],title.name = stage,font.size.title = 20,sort.by.target  = T)$data
  tmp_Ery2Ery_df$stage=stage;rownames(tmp_Ery2Ery_df)=NULL
  Ery2Ery_df2=rbind(Ery2Ery_df2,tmp_Ery2Ery_df)
  Ery2Ery_df2=Ery2Ery_df2[!is.na(Ery2Ery_df2$prob),]
  
}

Other2Ery_df2$target_type='Other2Ery'
Ery2Other_df2$target_type='Ery2Other'
Ery2Ery_df2$target_type='Ery2Ery'

YS_NRBC_altas_LR_df2=rbind(rbind(Other2Ery_df2,Ery2Other_df2),Ery2Ery_df2)
YS_NRBC_altas_LR_df2=YS_NRBC_altas_LR_df2[!is.na(YS_NRBC_altas_LR_df2$prob),]
YS_NRBC_altas_LR_df2=YS_NRBC_altas_LR_df2[-grep(pattern = 'ERYTHROID',YS_NRBC_altas_LR_df2$source.target)]
YS_NRBC_altas_LR_df2=YS_NRBC_altas_LR_df2[!is.na(YS_NRBC_altas_LR_df2$pathway_name),]
YS_NRBC_altas_LR_df2=YS_NRBC_altas_LR_df2[!is.na(YS_NRBC_altas_LR_df2$interaction_name),]


rownames(YS_NRBC_altas_LR_df)=NULL
write.table(Other2Ery_df,file = 'res_data/Other2Ery_df_new.csv',quote = F,sep = "\t")
write.table(YS_NRBC_altas_LR_df,file = 'res_data/YS_NRBC_altas_LR_df_new.csv',quote = F,sep = "\t")

write.table(Other2Ery_df2,file = 'res_data/Other2Ery_df2.csv',quote = F,sep = "\t")
write.table(YS_NRBC_altas_LR_df2,file = 'res_data/YS_NRBC_altas_LR_df2.csv',quote = F,sep = "\t")


all_incoming_LRs=unique(YS_NRBC_altas_LR_df$interaction_name[YS_NRBC_altas_LR_df$target_type=='Other2Ery'])
all_outgoing_LRs=unique(YS_NRBC_altas_LR_df$interaction_name[YS_NRBC_altas_LR_df$target_type=='Ery2Other'])
all_Ery2Ery_LRs=unique(YS_NRBC_altas_LR_df$interaction_name[YS_NRBC_altas_LR_df$target_type=='Ery2Ery'])

all_incoming_LRs2=unique(YS_NRBC_altas_LR_df2$interaction_name[YS_NRBC_altas_LR_df2$target_type=='Other2Ery'])
all_outgoing_LRs2=unique(YS_NRBC_altas_LR_df2$interaction_name[YS_NRBC_altas_LR_df2$target_type=='Ery2Other'])
all_Ery2Ery_LRs2=unique(YS_NRBC_altas_LR_df2$interaction_name[YS_NRBC_altas_LR_df2$target_type=='Ery2Ery'])

length(all_incoming_LRs) # 28
table(all_incoming_LRs2 %in% all_incoming_LRs) # F/T:1/28 
table(all_incoming_LRs %in% all_incoming_LRs2) # F/T:0/28 
length(all_outgoing_LRs)# 53
table(all_outgoing_LRs2 %in% all_outgoing_LRs) # F/T:5/49
table(all_outgoing_LRs %in% all_outgoing_LRs2) # F/T:4/49

length(all_Ery2Ery_LRs) # 5 
table(all_Ery2Ery_LRs2 %in% all_Ery2Ery_LRs) # F/T:0/5 
table(all_Ery2Ery_LRs %in% all_Ery2Ery_LRs2) # F/T:0/5 

rm(object.list2);gc()

# 粗颗粒细胞分类在更细颗粒细胞分类中的与NRBC细胞互作绝大部分都存在，少量未有发现的，基于表达水平,可以不考虑
all_outgoing_LRs[!all_outgoing_LRs %in%  all_outgoing_LRs2]#  SELPLG_SELP SELPLG_SELL HLA-A_CD8A  HLA-C_CD8A
# 来自NRBC 到免疫细胞，由于免疫细胞很少，
temp_df=YS_NRBC_altas_LR_df[YS_NRBC_altas_LR_df$interaction_name %in% c('SELPLG_SELP', 'SELPLG_SELL', 'HLA-A_CD8A',  'HLA-C_CD8A'),]# MESOTHELIUM target BFUE/CFUE, HEPATOCYTE target YS_Bas/Poly
VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features = c('SELPLG','HLA-A','HLA-C'),stack=T) #  增加到 YS_NRBC_altas_LR_df2

YS_NRBC_altas_LR_df2=rbind(YS_NRBC_altas_LR_df2,temp_df)

rm(stage_type_list,object.list);gc()

# 更细的类可以得到哪类细胞是主要ligand的提供者抛弃粗分类，直接采用细分类进行
YS_NRBC_altas_LR_df=YS_NRBC_altas_LR_df2
YS_NRBC_altas_LR_df$main_type='no'
YS_NRBC_altas_LR_df$main_type[YS_NRBC_altas_LR_df$target_type=='Other2Ery' & YS_NRBC_altas_LR_df$interaction_name %in% all_incoming_LRs]='yes'
YS_NRBC_altas_LR_df$main_type[YS_NRBC_altas_LR_df$target_type=='Ery2Other' & YS_NRBC_altas_LR_df$interaction_name %in% all_outgoing_LRs]='yes'
YS_NRBC_altas_LR_df$main_type[YS_NRBC_altas_LR_df$target_type=='Ery2Ery' & YS_NRBC_altas_LR_df$interaction_name %in% all_Ery2Ery_LRs]='yes'

YS_NRBC_altas_LR_df$celltype=as.character(YS_NRBC_altas_LR_df$target)
YS_NRBC_altas_LR_df$celltype[YS_NRBC_altas_LR_df$target_type=='Ery2Other']=as.character(YS_NRBC_altas_LR_df$source[YS_NRBC_altas_LR_df$target_type=='Ery2Other'])

NRBC_subcelltype=c("EryP-CFC","Bas","Poly","Orth" )
YS_NRBC_altas_LR_df$celltype=factor(YS_NRBC_altas_LR_df$celltype,levels = NRBC_subcelltype)

NRBC_celltype_mexp_df=data.frame(AverageExpression(subset(NRBC_altas_seu,tissue_stage=='YS'),  layer = 'data',group.by = 'source_celltype')$RNA)
# 去掉不同niche中 NRBC LR 表达均值低于0.1的基因 对应的LR
temp_LRs=YS_NRBC_altas_LR_df[YS_NRBC_altas_LR_df$stage==stage,]
LR_genes=unique(c(CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name  %in% unique(temp_LRs[temp_LRs$target_type=='Other2Ery','interaction_name']),'receptor.symbol'],
                  CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name  %in% unique(temp_LRs[temp_LRs$target_type=='Ery2Other','interaction_name']),'ligand.symbol']))

LR_genes=unique(unlist(strsplit(LR_genes,split=', ')))
LR_gene_exp_test=rowSums(NRBC_celltype_mexp_df[LR_genes,grep(stage,colnames(NRBC_celltype_mexp_df))] >0.1)
print(names(LR_gene_exp_test)[LR_gene_exp_test <1])
# 无

Other2Ery_df=YS_NRBC_altas_LR_df[YS_NRBC_altas_LR_df$target_type=='Other2Ery',]

write.table(YS_NRBC_altas_LR_df,file = 'primitive_NRBC_project/res_data/filt_YS_NRBC_altas_LR_df.csv',quote = F,sep="\t")
write.table(Other2Ery_df,file = 'primitive_NRBC_project/res_data/filt_YS_Other2Ery_df_new.csv',quote = F,sep="\t")

rm(list=ls())
gc()

Idents(NRBC_altas_seu)='source_celltype'
 
##############################################################################################################################################################################
#------------------------------------part2: analysis   interactioin  across niches & primitive vs definitiive----------------------#
##############################################################################################################################################################################

#################-----------------------2.1 DEGs-LR -------------------------------------------------------#

setwd('primitive_NRBC_project')
NRBC_altas_seu=readRDS('all_NRBC_altas_seu.rds')

definitive_NRBC_df=read.csv('../NRBC_altas_CC/res_data/def_filt_NRBC_altas_LR_df.csv',sep="\t")
NRBC_altas_LR_df=rbind(YS_NRBC_altas_LR_df,definitive_NRBC_df)

NRBC_subcelltype=c("BFUE/CFUE","EryP-CFC","ProE","Bas","Poly","Orth" )
NRBC_altas_LR_df$stage=factor(NRBC_altas_LR_df$stage,levels = c('YS','FL','FBM','ABM'))
NRBC_altas_LR_df$celltype=factor(NRBC_altas_LR_df$celltype,levels = NRBC_subcelltype)

 
Idents(NRBC_altas_seu)='source_celltype'
NRBC_altas_seu$pd_celltype='definitive'
NRBC_altas_seu$pd_celltype[grep('YS',NRBC_altas_seu$tissue_stage)]='primitive'


all_NRBC_receptor_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']),'receptor.symbol'] # 取receptor gene

all_NRBC_ligand_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in%unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']),'ligand.symbol']# 取ligand gene
Ery2Ery_ligands=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Ery']),'ligand.symbol'] # 取receptor gene
Ery2Ery_ligands[!Ery2Ery_ligands %in% all_NRBC_ligand_genes] # 全部在all_NRBC_ligand_genes


all_NRBC_receptor_genes=sort(unique(unlist(strsplit(all_NRBC_receptor_genes,split=','))))
all_NRBC_ligand_genes=sort(unique(unlist(strsplit(all_NRBC_ligand_genes,split=','))))
length(all_NRBC_receptor_genes);length(all_NRBC_ligand_genes)# 49， 76

all_NRBC_LR_genes=unique(c(all_NRBC_receptor_genes,all_NRBC_ligand_genes));length(all_NRBC_LR_genes) # 116


# ----------------whole level------------------#
pd_whole_level_markers_degs_df=readRDS('res_data/pd_whole_level_markers.rds')

receptor_pd_whole_level_markers_degs_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$gene %in% all_NRBC_receptor_genes ,]
receptor_pd_whole_level_markers_degs_df=receptor_pd_whole_level_markers_degs_df[receptor_pd_whole_level_markers_degs_df$cluster=='primitive',]
receptor_pd_whole_level_markers_degs_df$qval=-log10(receptor_pd_whole_level_markers_degs_df$p_val_adj)
receptor_pd_whole_level_markers_degs_df[receptor_pd_whole_level_markers_degs_df$qval >150,'qval']=150
receptor_pd_whole_level_markers_degs_df$sig="Not significant"
receptor_pd_whole_level_markers_degs_df$sig[receptor_pd_whole_level_markers_degs_df$avg_log2FC >0]="Up"
receptor_pd_whole_level_markers_degs_df$sig[receptor_pd_whole_level_markers_degs_df$avg_log2FC < 0]="Dn"

dir.create('res_data/main_figure4')

# 火山图
p1=ggplot(receptor_pd_whole_level_markers_degs_df, aes(x = avg_log2FC, y = qval, color = sig)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#E64B35", "Dn" = "#4DBBD5")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(aes(label = gene), size = 3, max.overlaps = 20, color = "black", fontface = "bold") +
  labs(title = "Primitive vs Definitive Receptors", x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "top")
p1

ggsave(p1,filename='res_pic/main_figure4/NRBC_receptor_DE_whole_level_dotplot.pdf',width = 6,height =6 )

ligand_pd_whole_level_markers_degs_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$gene %in% all_NRBC_ligand_genes ,]
ligand_pd_whole_level_markers_degs_df=ligand_pd_whole_level_markers_degs_df[ligand_pd_whole_level_markers_degs_df$cluster=='primitive',]
ligand_pd_whole_level_markers_degs_df$qval=-log10(ligand_pd_whole_level_markers_degs_df$p_val_adj)
ligand_pd_whole_level_markers_degs_df[ligand_pd_whole_level_markers_degs_df$qval >150,'qval']=150
ligand_pd_whole_level_markers_degs_df$sig="Not significant"
ligand_pd_whole_level_markers_degs_df$sig[ligand_pd_whole_level_markers_degs_df$avg_log2FC >0]="Up"
ligand_pd_whole_level_markers_degs_df$sig[ligand_pd_whole_level_markers_degs_df$avg_log2FC < 0]="Dn"

# 火山图
p2=ggplot(ligand_pd_whole_level_markers_degs_df, aes(x = avg_log2FC, y = qval, color = sig)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("Up" = "#E64B35", "Dn" = "#4DBBD5")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(aes(label = gene), size =3, max.overlaps = 20, color = "black", fontface = "bold") +
  labs(title = "Primitive vs Definitive Ligands", x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "top")
p2

ggsave(p2,filename='res_pic/main_figure4/NRBC_ligand_DE_whole_level_dotplot.pdf',width = 6,height =6 )


# -----------early stage--------------#

# subcelltype level
sub_pd_all_Ery_tissue_markers=read.csv('res_data/primitive_definitive_all_Ery_RNA_markers.csv')
sub_pd_all_Ery_tissue_markers=sub_pd_all_Ery_tissue_markers[,-1]

sub_pd_pic_list=list()
for( celltype in unique( sub_pd_all_Ery_tissue_markers$celltype)){
  
  receptor_pd_sub_level_markers_degs_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$gene %in% all_NRBC_receptor_genes & sub_pd_all_Ery_tissue_markers$celltype==celltype ,]
  receptor_pd_sub_level_markers_degs_df=receptor_pd_sub_level_markers_degs_df[receptor_pd_sub_level_markers_degs_df$cluster=='primitive',]
  receptor_pd_sub_level_markers_degs_df$qval=-log10(receptor_pd_sub_level_markers_degs_df$p_val_adj)
  receptor_pd_sub_level_markers_degs_df[receptor_pd_sub_level_markers_degs_df$qval >150,'qval']=150
  receptor_pd_sub_level_markers_degs_df$sig="Not significant"
  receptor_pd_sub_level_markers_degs_df$sig[receptor_pd_sub_level_markers_degs_df$avg_log2FC >0]="Up"
  receptor_pd_sub_level_markers_degs_df$sig[receptor_pd_sub_level_markers_degs_df$avg_log2FC < 0]="Dn"
  
  # 火山图
  p1=ggplot(receptor_pd_sub_level_markers_degs_df, aes(x = avg_log2FC, y = qval, color = sig)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#E64B35", "Dn" = "#4DBBD5")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(aes(label = gene), size = 3, max.overlaps = 20, color = "black", fontface = "bold") +
    labs(title = "Primitive vs Definitive Receptors", x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "top")
  sub_pd_pic_list[[celltype]][['receptor']]=p1
  
  ligand_pd_sub_level_markers_degs_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$gene %in% all_NRBC_ligand_genes & sub_pd_all_Ery_tissue_markers$celltype==celltype ,]
  ligand_pd_sub_level_markers_degs_df=ligand_pd_sub_level_markers_degs_df[ligand_pd_sub_level_markers_degs_df$cluster=='primitive',]
  ligand_pd_sub_level_markers_degs_df$qval=-log10(ligand_pd_sub_level_markers_degs_df$p_val_adj)
  ligand_pd_sub_level_markers_degs_df[ligand_pd_sub_level_markers_degs_df$qval >150,'qval']=150
  ligand_pd_sub_level_markers_degs_df$sig="Not significant"
  ligand_pd_sub_level_markers_degs_df$sig[ligand_pd_sub_level_markers_degs_df$avg_log2FC >0]="Up"
  ligand_pd_sub_level_markers_degs_df$sig[ligand_pd_sub_level_markers_degs_df$avg_log2FC < 0]="Dn"
  
  # 火山图
  p2=ggplot(ligand_pd_sub_level_markers_degs_df, aes(x = avg_log2FC, y = qval, color = sig)) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#E64B35", "Dn" = "#4DBBD5")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(aes(label = gene), size =3, max.overlaps = 20, color = "black", fontface = "bold") +
    labs(title = "Primitive vs Definitive Ligands", x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "top")
  sub_pd_pic_list[[celltype]][['ligand']]=p2
  
  
}

p=sub_pd_pic_list[['early_Ery']][[1]]+sub_pd_pic_list[['early_Ery']][[2]];p
ggsave(p2,filename='res_pic/main_figure4/NRBC_receptor_ligand_DE_sub_early_dotplot.pdf',width = 12,height =5 )
p=sub_pd_pic_list[['mid_Ery']][[1]]+sub_pd_pic_list[['mid_Ery']][[2]];p
ggsave(p2,filename='res_pic/main_figure4/NRBC_receptor_ligand_DE_sub_mid_dotplot.pdf',width = 12,height =5 )
p=sub_pd_pic_list[['late_Ery']][[1]]+sub_pd_pic_list[['late_Ery']][[2]];p
ggsave(p2,filename='res_pic/main_figure4/NRBC_receptor_ligand_DE_sub_late_dotplot.pdf',width = 12,height =5 )




#----------------------recepted LR--------------------------#
LR_max_prob_func=function(LR_df){
  temp_max_prob_LRs_df=data.frame() # interaction,stage,celltype,target_type
  for(LR in unique(LR_df$interaction_name) ){
    temp_df=LR_df[LR_df$interaction_name ==LR ,]
    temp_df$stage_celltype=paste(as.character(temp_df$stage),as.character(temp_df$celltype),sep='_')
    if(dim(temp_df)[1]>1){
      for (celltype  in unique(temp_df$stage_celltype)) {
        temp_df1=temp_df[temp_df$stage_celltype==celltype,]
        if(length(unique(temp_df1$target_type)) >1){
          for( target_type in unique(temp_df1$target_type) ){
            temp_df2=temp_df1[temp_df1$target_type==target_type,]
            temp_df2=temp_df2[which.max(temp_df2$prob),]
            temp_max_prob_LRs_df=rbind(temp_max_prob_LRs_df,temp_df2)
          }
        }else{
          temp_df2=temp_df1[which.max(temp_df1$prob),]
          temp_max_prob_LRs_df=rbind(temp_max_prob_LRs_df,temp_df2)
          
        }
      }
    }else{
      temp_df2=temp_df[which.max(temp_df$prob),]
      temp_max_prob_LRs_df=rbind(temp_max_prob_LRs_df,temp_df2)
      
    }
    
  }
  colnames(temp_max_prob_LRs_df)=gsub(pattern ='prob' ,replacement = 'max.prob',colnames(temp_max_prob_LRs_df))
  return(temp_max_prob_LRs_df)
  
}

recepted_NRBC_altas_LR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$target_type!='Ery2Other',]
recepted_LR_list <- split(recepted_NRBC_altas_LR_df$interaction_name,   f = recepted_NRBC_altas_LR_df$stage)
recepted_LR_list <- lapply(recepted_LR_list, unique)
recepted_LR_list <- lapply(recepted_LR_list, as.character)
common_recepted_interactions <- Reduce(intersect, recepted_LR_list)
max_common_recepted_LR_df=LR_max_prob_func(LR_df =recepted_NRBC_altas_LR_df[recepted_NRBC_altas_LR_df$interaction_name %in% common_recepted_interactions,] )
p1=ggplot(max_common_recepted_LR_df,aes(x=celltype,y=interaction_name,size=max.prob,shape=target_type,color=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'NRBCs recepted commone interactions')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])
p1
ggsave(p1,filename='res_pic/main_figure4/NRBC_common_recepted.pdf',width = 10,height =8 )


#
outgoing_NRBC_altas_LR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$target_type=='Ery2Other',]
outgoing_LR_list <- split(outgoing_NRBC_altas_LR_df$interaction_name,   f = outgoing_NRBC_altas_LR_df$stage)
outgoing_LR_list <- lapply(outgoing_LR_list, unique)
outgoing_LR_list <- lapply(outgoing_LR_list, as.character)
common_outgoing_interactions <- Reduce(intersect, outgoing_LR_list)
max_common_outgoing_LR_df=LR_max_prob_func(LR_df =outgoing_NRBC_altas_LR_df[outgoing_NRBC_altas_LR_df$interaction_name %in% common_outgoing_interactions,] )

p2=ggplot(max_common_outgoing_LR_df,aes(x=celltype,y=interaction_name,size=max.prob,shape=target_type,color=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'NRBCs outgoing common interactions')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])
p2
ggsave(p1,filename='res_pic/main_figure4/NRBC_common_outgoing.pdf',width = 8,height =8 )


def_common_recepted_interactions=Reduce(intersect, recepted_LR_list[c('FL','FBM','ABM')])
specific_recepted_interactions <- c(recepted_LR_list[['YS']][! recepted_LR_list[['YS']] %in% c(recepted_LR_list[['FL']],recepted_LR_list['FBM'],recepted_LR_list[['ABM']])],
                                    def_common_recepted_interactions[!def_common_recepted_interactions %in% recepted_LR_list[['YS']]]
                                    )

max_specific_recepted_LR_df=LR_max_prob_func(LR_df =recepted_NRBC_altas_LR_df[recepted_NRBC_altas_LR_df$interaction_name %in% specific_recepted_interactions,] )

p21=ggplot(max_specific_recepted_LR_df,aes(x=celltype,y=interaction_name,size=max.prob,shape=target_type,color=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'NRBCs recepted specific interactions')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])
p21
ggsave(p21,filename='res_pic/main_figure4/NRBC_specific_recepted.pdf',width = 10,height =10 )

#

def_specific_outgoing_interactions=Reduce(intersect, outgoing_LR_list[c('FL','FBM','ABM')])
specific_outgoing_interactions <- c(outgoing_LR_list[['YS']][! outgoing_LR_list[['YS']] %in% c(outgoing_LR_list[['FL']],outgoing_LR_list[['FBM']],outgoing_LR_list[['ABM']])],
                                    def_common_outgoing_interactions[!def_common_outgoing_interactions %in% outgoing_LR_list[['YS']]]
)


max_specific_outgoing_LR_df=LR_max_prob_func(LR_df =outgoing_NRBC_altas_LR_df[outgoing_NRBC_altas_LR_df$interaction_name %in% specific_outgoing_interactions,] )

p22=ggplot(max_specific_outgoing_LR_df,aes(x=celltype,y=interaction_name,size=max.prob,shape=target_type,color=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'NRBCs outgoing specific interactions')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])
p22
ggsave(p22,filename='res_pic/main_figure4/NRBC_specific_recepted.pdf',width = 10,height =10 )


# key receptor or ligands

focused_receptor_genes=c('RORA','TMEM219','NR1H2','CD55','EPOR','CD44','CD47','CD36','ITGA4','CD99','PLXNC1','NOTCH2')

p=VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features =focused_receptor_genes ,
        stack = T,split.by = 'tissue_stage',cols = cols)  +NoLegend()
p

ggsave(p,filename='res_pic/main_figure4/NRBC_focused_receptor_vlnplot.pdf',width = 10,height =8 )


temp_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$gene %in%focused_receptor_genes & sub_pd_all_Ery_tissue_markers$cluster=='primitive', ]
temp_df$gene=factor(temp_df$gene,levels =focused_receptor_genes )
p1=ggplot(temp_df[temp_df$celltype=='early_Ery',],aes(x=gene,y=avg_log2FC ))+geom_bar(stat = 'identity',fill=alpha('firebrick3',0.8))+theme_classic()
p2=ggplot(temp_df[temp_df$celltype=='mid_Ery',],aes(x=gene,y=avg_log2FC ))+geom_bar(stat = 'identity',fill=alpha('firebrick3',0.8))+theme_classic()
p3=ggplot(temp_df[temp_df$celltype=='late_Ery',],aes(x=gene,y=avg_log2FC ))+geom_bar(stat = 'identity',fill=alpha('firebrick3',0.8))+theme_classic()

p22=p1/p2/p3;p22




############################################################################################################################################################
#-----------------------------part3: using the nichnet to analysis the key LR to target genes-------------------------#
############################################################################################################################################################
library(nichenetr)
library(network)
library(ggnetwork)


lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS("../NRBC_altas_CC/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("../NRBC_altas_CC/weighted_networks_nsga2r_final.rds")

sig_network <- readRDS('../NRBC_altas_CC/signaling_network_human_21122021.rds')
gr_network <- readRDS('../NRBC_altas_CC/gr_network_human_21122021.rds')
ligand_tf_matrix <- readRDS('../NRBC_altas_CC/ligand_tf_matrix_nsga2r_final.rds')

signaling_pic_func=function(ligands_oi,targets_oi,weighted_networks=weighted_networks,ligand_tf_matrix=ligand_tf_matrix,top_n_regulators=4){
  
  active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,targets_all = targets_oi,
                                                        weighted_networks = weighted_networks,ligand_tf_matrix = ligand_tf_matrix,
                                                        top_n_regulators = top_n_regulators,minmax_scaling = TRUE) 
  
  graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                     ligands_all = ligands_oi, targets_all = targets_oi,
                                                     sig_color = "indianred", gr_color = "steelblue")
  
  
  #构建 network 对象
  # 使用 edges_df 中的 from 和 to 列构建有向图
  edges_df <- graph_min_max$edges_df
  nodes_df <- graph_min_max$nodes_df
  net_obj <- network(as.matrix(edges_df[, c("from", "to")]), directed = TRUE)
  
  # 将节点属性合并到 network 对象中
  # 这一步至关重要，它让 ggnetwork 能读取节点的颜色、类型等信息
  # 确保节点名称匹配
  net_vertices <- network.vertex.names(net_obj)
  
  # 将属性赋值给 network 对象的顶点
  # 使用 match 确保顺序一致
  net_obj %v% "label" <- nodes_df$label[match(net_vertices, nodes_df$id)]
  net_obj %v% "type" <- nodes_df$type[match(net_vertices, nodes_df$id)] # 节点类型（配体/靶点/中间分子）
  net_obj %v% "fillcolor" <- nodes_df$fillcolor[match(net_vertices, nodes_df$id)] # 节点颜色
  
  # 为了映射边的粗细，我们需要将边的权重信息也加入到 network 对象中
  # 注意：network 对象构建时边的顺序可能与 edges_df 不完全一致，需要重新匹配
  # 获取 network 对象中的边列表
  net_edges <- as.data.frame(net_obj, edge.names = TRUE)
  # 匹配 from 和 to 来合并权重
  net_edges$weight <- edges_df$penwidth[match(paste(net_edges$.tail, net_edges$.head, sep = "-"), 
                                              paste(edges_df$from, edges_df$to, sep = "-"))]
  # 将权重赋值给 network 对象的边属性
  net_obj %e% "weight" <- net_edges$weight
  
  # 开始绘图
  p <- ggplot(net_obj, aes(x = x, y = y, xend = xend, yend = yend))+
    # 绘制边 使用 weight 映射线条粗细
    geom_edges(aes(size = weight),alpha = 0.4,  color = "grey70", arrow = arrow(length = unit(1, "mm"), type = "closed")) +
    # 绘制节点 使用 fillcolor 映射填充颜色
    geom_nodes(aes(color = fillcolor), size = 5) +
    # 添加节点标签
    geom_nodetext(aes(label = label),   fontface = "bold", color = "black", size = 3, vjust = 1.5) + # vjust 调整标签位置
    # 主题调整, 移除背景网格
    theme_blank() + labs(title = "NicheNet Signaling Pathway") +
    # 手动设置颜色（保持 NicheNet 原有的配色）
    scale_color_identity() + scale_size_continuous(range = c(0.5, 2)) # 调整线条粗细范围
  
  print(unique(net_obj$label))
  return(p)
}

draw_gonetcwork_pic_func=function(res=fa_fetal_LR_targetgene_enrichGO_res,showCategory=20,xlimits = c(-0.2, 2.5)){
  
  res <- as.data.frame(res)
  res <- head(res,showCategory)
  
  # 构建连接
  links <- do.call(rbind, lapply(1:nrow(res), function(i) {
    data.frame(gene = strsplit(res$geneID[i], "/")[[1]], 
               pathway = res$Description[i])}))
  
  # 计算位置
  links$y_gene <- match(links$gene, names(sort(table(links$gene), T))) / (length(unique(links$gene)) + 1)
  links$y_pathway <- match(links$pathway, res$Description[order(res$p.adjust)]) / (nrow(res) + 1)
  
  # 绘图
  p=ggplot() +  geom_segment(aes(x = 0, xend = 0.8, y = y_gene, yend = y_pathway), data = links,   color = "gray70", alpha = 0.4, size = 0.4) +
    geom_point(aes(x = 0, y = y_gene), data = distinct(links, gene, y_gene), size = 2, color = "#4878D0") +
    geom_text(aes(x = -0.05, y = y_gene, label = gene), data = distinct(links, gene, y_gene), hjust = 1, size = 3, color = "black") +
    geom_point(aes(x = 0.8, y = y_pathway), data = distinct(links, pathway, y_pathway), size = 2, color = "#EE854A") +
    geom_text(aes(x = 0.85, y = y_pathway, label = pathway), data = distinct(links, pathway, y_pathway),   hjust = 0, size = 2.5, color = "black") +
    scale_x_continuous(limits =xlimits, breaks = c(0, 1), labels = c("Genes", "Pathways")) +
    scale_y_continuous(breaks = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 10),
          plot.margin = margin(5, 15, 25, 5))
  print(p)
  return(p)
}



LR_df=CellChatDB.human$interaction
LR_df$ligand.symbol[LR_df$ligand.symbol==""]=LR_df$ligand[LR_df$ligand.symbol==""]

# --------compare the LR on primitive and definitve NRBC -------------#
primtive_NBRC_ligands=unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='YS' & NRBC_altas_LR_df$target_type!='Ery2Other','ligand'])
primtive_NBRC_ligands= unique(unlist(strsplit(as.character(t(data.frame(strsplit(primtive_NBRC_ligands,split = '_')))[,1]),split = '-')))
primtive_NBRC_ligands=primtive_NBRC_ligands[primtive_NBRC_ligands %in% rownames(NRBC_altas_seu)]  



nichenet_predict_func=function(potential_ligands,target_genes,ligand_target_matrix,background_expressed_genes,legend_title='primitive vs definitve score'){
  predicted_pd_ligand_activites=predict_ligand_activities(geneset =target_genes,ligand_target_matrix = ligand_target_matrix,
                                                          potential_ligands = potential_ligands, background_expressed_genes =background_expressed_genes )
  
  predicted_ligand_activites_mtx=predicted_pd_ligand_activites[order(predicted_pd_ligand_activites$aupr_corrected),]  %>% 
    column_to_rownames('test_ligand') %>% dplyr::select(aupr_corrected) %>% as.matrix(ncol = 1)                                                        
  p=make_heatmap_ggplot(matrix =predicted_ligand_activites_mtx,y_name ='Ligand activaty' ,x_name = "Prioritized ligands",legend_title = "AUPR", color = "darkorange")+
    theme(axis.text.x.top = element_blank())  
  
  # target the DEGs
  ligand_target_gene_link_df <-predicted_pd_ligand_activites$test_ligand%>% lapply(get_weighted_ligand_target_links,geneset =target_genes,ligand_target_matrix = ligand_target_matrix, n = 100) %>%bind_rows() %>% drop_na()
  ligand_target_gene_link_vis=prepare_ligand_target_visualization(ligand_target_df =ligand_target_gene_link_df,ligand_target_matrix = ligand_target_matrix,cutoff = 0.25)
  p2=make_heatmap_ggplot(matrix =t(ligand_target_gene_link_vis[,rownames(predicted_ligand_activites_mtx)[rownames(predicted_ligand_activites_mtx) %in% colnames(ligand_target_gene_link_vis)]]),y_name = 'ligand',x_name = 'target',legend_title =legend_title )
  
  return(list(p1=p,p2=p2,link_vis=ligand_target_gene_link_vis))
  
}
NRBC_altas_seu$pd_celltype2=paste(as.character(NRBC_altas_seu$pd_celltype),as.character(NRBC_altas_seu$final_celltype),sep = "_")
all_mexp_df=as.matrix(AverageExpression(NRBC_altas_seu,group.by = 'pd_celltype2',features =rownames(NRBC_altas_seu) )$RNA)

#primtive_NBRC_ligands=primtive_NBRC_ligands[-grep('CXCL|CCL',primtive_NBRC_ligands)]   

# cho the target geneset-#
# primitive vs definitive
top_pd_degs=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$avg_log2FC >1 & pd_whole_level_markers_degs_df$pct.1 >0.05, ] %>% group_by(cluster) %>% top_n(wt =avg_log2FC,n = 800 ) #& pd_whole_level_markers$pct.2 <0.3
table(top_pd_degs[,c('cluster')])

sub_pd_all_Ery_tissue_maconserved_LR_genesrkers=sub_pd_all_Ery_tissue_markers[,-1]
top_pd_sub_degs=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >1 & sub_pd_all_Ery_tissue_markers$pct.1>0.1  ,]  %>% group_by(cluster,celltype) %>% top_n(wt =avg_log2FC,n = 500 )
table(top_pd_sub_degs$cluster)
table(top_pd_sub_degs[,c('cluster','celltype')])

primitive_positive_target_genes=unique(c(top_pd_degs$gene[top_pd_degs$cluster=='primitive'],top_pd_sub_degs$gene[top_pd_sub_degs$cluster=='primitive']))
definitive_positive_target_genes=unique(c(top_pd_degs$gene[top_pd_degs$cluster=='definitive'],top_pd_sub_degs$gene[top_pd_sub_degs$cluster=='definitive']))
definitive_positive_target_genes=definitive_positive_target_genes[-grep('^HIS',definitive_positive_target_genes)]


known_key_regulator_erythropoiesis_genes=c('GATA1','GATA2','KLF1','TAL1','LMO2','FOG1','BCL11A','EPOR','HIF1A','HIF2A','STAT5A','ALAS2','SLC11A2','FOXO3')
known_key_regulator_erythropoiesis_genes[known_key_regulator_erythropoiesis_genes %in% pd_whole_level_markers_degs_df$gene]
known_key_regulator_erythropoiesis_genes=known_key_regulator_erythropoiesis_genes
#known_key_regulator_erythropoiesis_genes[known_key_regulator_erythropoiesis_genes %in%primitive_positive_target_genes ] #  "KLF1"  "ALAS2"
#known_key_regulator_erythropoiesis_genes[known_key_regulator_erythropoiesis_genes %in%definitive_positive_target_genes ] #  "GATA2"  "BCL11A" "STAT5A"

potential_ligands=primtive_NBRC_ligands;potential_ligands=potential_ligands[potential_ligands %in% colnames(ligand_target_matrix)]

# 排除低表达基因
target_genes=unique(c(primitive_positive_target_genes)) ;target_genes =unique(c(target_genes,known_key_regulator_erythropoiesis_genes));length(target_genes)
target_genes=target_genes[target_genes %in% rownames(all_mexp_df)]
target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('primi',colnames(all_mexp_df))])>0.1] ;length(target_genes)

background_expressed_genes=rownames(NRBC_altas_seu)[!rownames(NRBC_altas_seu)%in%  target_genes]

primitive_ligand_target_res_list=nichenet_predict_func(potential_ligands = potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
p=primitive_ligand_target_res_list[[2]];p
ggsave(p,filename='res_pic/main_figure4/pd_primitive_ligand_target_gene_heatmap_score.pdf',width = 10,height = 6)
saveRDS(primitive_ligand_target_res_list,file = 'res_data/pd_primitive_ligand_target_res_list.rds')

primitive_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(primitive_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(primitive_LR_targetgene_enrichGO_res,file = 'res_data/primitive_LR_targetgene_enrichGO_res.rds')

p=draw_gonetcwork_pic_func(res =primitive_LR_targetgene_enrichGO_res,showCategory =20,xlimits = c(-0.2,2.2) );p
ggsave(p,filename='res_pic/main_figure4/pd_primitive_target_gene_enrichGO_res.pdf',width =6,height = 6)

p=VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features = 'CDKN1A',raster =T,pt.size = 0.5)+NoLegend()
ggsave(p,filename='res_pic/main_figure4/CDKN1A_vlnplot.pdf',width =10,height = 6)

cho_LR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='YS' & NRBC_altas_LR_df$receptor %in% focused_receptor_genes,]
cho_LR_df=cho_LR_df[cho_LR_df$target_type!='Ery2Other',]
unique(cho_LR_df$ligand)
# "ADGRE5"             "Desmosterol-DHCR7"  "EPO"                "IGFBP3"             "THBS1"              "Cholesterol-DHCR24" "Cholesterol-LIPA"   "Cholesterol-DHCR7" 
cho_LR_d2f=NRBC_altas_LR_df[NRBC_altas_LR_df$stage=='YS' & NRBC_altas_LR_df$interaction_name %in% specific_recepted_interactions,]
cho_LR_df2=cho_LR_d2f[cho_LR_d2f$target_type!='Ery2Other',]
unique(cho_LR_df2$ligand)
#"CXCL1"              "Cholesterol-DHCR24"

cho_potential_ligands=c('IGFBP3','EPO');
cho_potential_ligands=cho_potential_ligands[cho_potential_ligands %in% colnames(ligand_target_matrix)]

# -------------------up DEG----------------------#
primitive_ligand_target_res_list2=nichenet_predict_func(potential_ligands = cho_potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
p=primitive_ligand_target_res_list2[[2]];p
ggsave(p,filename='res_pic/main_figure4/EPO_IGFBP3_pd_primitive_ligand_target_gene_heatmap_score.pdf',width = 8,height = 3)

VlnPlot(NRBC_altas_seu,features = rownames(primitive_ligand_target_res_list2[[3]]),group.by = 'source_celltype',stack = T,split.by = 'tissue_stage')+NoLegend()
saveRDS(primitive_ligand_target_res_list2,file = 'res_data/IGFBP3_pd_primitive_ligand_target_res_up_list.rds')


#----------------------- Dn DEGs------------------------------#
cho_potential_ligands=c('IGFBP3','EPO');
cho_potential_ligands=cho_potential_ligands[cho_potential_ligands %in% colnames(ligand_target_matrix)]

target_genes=unique(c(definitive_positive_target_genes)) ;length(target_genes)
target_genes=target_genes[target_genes %in% rownames(all_mexp_df)]
target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('def',colnames(all_mexp_df))])>0.1] ;length(target_genes)
background_expressed_genes=rownames(NRBC_altas_seu)[!rownames(NRBC_altas_seu)%in%  target_genes]

primitive_ligand_target_res_list3=nichenet_predict_func(potential_ligands = cho_potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
p=primitive_ligand_target_res_list3[[2]];p
ggsave(p,filename='res_pic/main_figure4/EPO_IGFBP3_pd_primitive_ligand_target_dn_gene_heatmap_score.pdf',width = 6,height = 3)
saveRDS(primitive_ligand_target_res_list3,file = 'res_data/IGFBP3_pd_primitive_ligand_target_res_dn_list.rds')


cho_primitive_LR_targetgene_enrichGO_res=enrichGO(gene =unique(c(rownames(primitive_ligand_target_res_list2[[3]]), rownames(primitive_ligand_target_res_list3[[3]]))),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(cho_primitive_LR_targetgene_enrichGO_res,file = 'res_data/IGFBP3_EPO_primitive_LR_targetgene_enrichGO_res.rds')

p=draw_gonetcwork_pic_func(res =cho_primitive_LR_targetgene_enrichGO_res,showCategory =20,xlimits = c(-0.2,2.2) );p
ggsave(p,filename='res_pic/main_figure4/IGFBP3_EPO_pd_primitive_target_gene_enrichGO_res.pdf',width =8,height = 6)

{
  
  # YS celltype DEGs 
  filt_YS_altas_Ery_seu_celltype_Markers=readRDS('res_data/filt_YS_altas_Ery_seu_celltype_Markers.rds')
  top500_filt_YS_altas_Ery_seu_celltype_Markers =filt_YS_altas_Ery_seu_celltype_Markers  %>% filter(p_val_adj < 0.0001, avg_log2FC > 0.5)  %>% group_by(cluster) %>% 
    arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 500) %>% ungroup()
  length(unique(top500_filt_YS_altas_Ery_seu_celltype_Markers$gene))
  target_genes=unique(top500_filt_YS_altas_Ery_seu_celltype_Markers$gene)
  
  target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('prim',colnames(all_mexp_df))])>0.1] ;length(target_genes)
  background_expressed_genes=rownames(NRBC_altas_seu)[!rownames(NRBC_altas_seu)%in%  target_genes]
  
  primitive_ligand_target_res_list4=nichenet_predict_func(potential_ligands = cho_potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
  p=primitive_ligand_target_res_list4[[2]];p
  
  # Definitive celltype DEGs
  df_celltype_Markers=FindAllMarkers(subset(NRBC_altas_seu,pd_celltype=='definitive'),group.by = 'final_celltype')
  saveRDS(df_celltype_Markers,filename='res_data/df_celltype_Markers.rds')
  
  top500_filt_YS_altas_Ery_seu_celltype_Markers =filt_YS_altas_Ery_seu_celltype_Markers  %>% filter(p_val_adj < 0.0001, avg_log2FC > 0.5)  %>% group_by(cluster) %>% 
    arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 500) %>% ungroup()
  length(unique(top500_filt_YS_altas_Ery_seu_celltype_Markers$gene))
  target_genes=unique(top500_filt_YS_altas_Ery_seu_celltype_Markers$gene)
  
  target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('prim',colnames(all_mexp_df))])>0.1] ;length(target_genes)
  background_expressed_genes=rownames(NRBC_altas_seu)[!rownames(NRBC_altas_seu)%in%  target_genes]
  
  primitive_ligand_target_res_list5=nichenet_predict_func(potential_ligands = cho_potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
  p=primitive_ligand_target_res_list5[[2]];p
  
  
}

ligands_oi <- c('IGFBP3') # this can be a list of multiple ligands if required
targets_oi <- names(primitive_ligand_target_res_list2[[3]][,'IGFBP3'][primitive_ligand_target_res_list2[[3]][,'IGFBP3'] >0.01])#"IFITM2
IGFBP3_signaling=signaling_pic_func(ligands_oi =ligands_oi,targets_oi = targets_oi,weighted_networks =weighted_networks ,ligand_tf_matrix =ligand_tf_matrix ,top_n_regulators =4  )
all_related_genes=unique(IGFBP3_signaling$data$label)
all_related_TF_genes=all_related_genes[!all_related_genes %in% c(targets_oi,ligands_oi)]
DotPlot(NRBC_altas_seu,features = all_related_TF_genes,group.by = 'source_celltype',scale = T)+RotatedAxis()


temp_df=GetAssayData(object =subset(NRBC_altas_seu, TMEM219 >0 & EPOR>0),layer='data' )
temp_df=data.frame(t(temp_df[c('TMEM219','EPOR'),]))
temp_df$tissue_stage=NRBC_altas_seu@meta.data[rownames(temp_df),'tissue_stage']
temp_df$celltype=NRBC_altas_seu@meta.data[rownames(temp_df),'final_celltype']

p=ggplot(temp_df,aes(x=TMEM219,y=EPOR))+geom_point(aes(color=celltype))+facet_grid(~tissue_stage)+scale_color_manual(values = cols)+
  geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
  stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
p

ggsave(p,filename='main_figure5/YS_TMEM219_EPOR_correlation.pdf',height = 5,width = 12)



############################################################################################################################################################
#---------------------------------------definitive ligand to target genes------------------------------#
############################################################################################################################################################

definitive_NBRC_ligands=unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage!='YS' & NRBC_altas_LR_df$target_type!='Ery2Other','ligand'])
definitive_NBRC_ligands=unique(LR_df[LR_df$ligand %in%  definitive_NBRC_ligands,'ligand.symbol'])
definitive_NBRC_ligands=unique(unlist(strsplit(definitive_NBRC_ligands,', ')))


potential_ligands=definitive_NBRC_ligands;potential_ligands=potential_ligands[potential_ligands %in% colnames(ligand_target_matrix)]
target_genes=unique(c(definitive_positive_target_genes))
target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('defi',colnames(all_mexp_df))])>0.1]
background_expressed_genes=rownames(NRBC_altas_seu)[!rownames(NRBC_altas_seu)%in%  target_genes]

definitive_ligand_target_res_list=nichenet_predict_func(potential_ligands =potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes,legend_title = 'definitve ligand to target genes' )
p=definitive_ligand_target_res_list[[2]];p
ggsave(p,filename='res_pic/main_figure4/raw_pd_definitive_ligand_target_gene_heatmap_score.pdf',width = 6,height = 8)


definitive_NBRC_ligands=unique(LR_df[LR_df$interaction_name %in%  def_common_recepted_interactions,'ligand.symbol'])
potential_ligands=c('EPO',definitive_NBRC_ligands);potential_ligands=potential_ligands[potential_ligands %in% colnames(ligand_target_matrix)]
target_genes=unique(c(definitive_positive_target_genes))
target_genes=target_genes[rowMax(all_mexp_df[target_genes,grep('defi',colnames(all_mexp_df))])>0.1]
background_expressed_genes=rownames(NRBC_altas_seu)[!rownames(NRBC_altas_seu)%in%  target_genes]

definitive_ligand_target_res_list=nichenet_predict_func(potential_ligands =potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes,legend_title = 'definitve ligand to target genes' )
p=definitive_ligand_target_res_list[[2]];p
ggsave(p,filename='res_pic/main_figure4/pd_definitive_ligand_target_gene_heatmap_score.pdf',width = 8,height = 8)


#
an_df=data.frame(definitive_ligand_target_res_list[[1]]$data[definitive_ligand_target_res_list[[1]]$data$y %in% c('EPO','MDK','APP'),])
an_df=data.frame(row.names = an_df$y,aupr_score=an_df$score)

p=pheatmap(t(definitive_ligand_target_res_list[[3]][,c('EPO','APP','MDK')]),annotation_row = an_df,border_color = 'white',cluster_rows = F,cluster_cols =F ,color =colorRampPalette(colors =  c("white" ,"#EF6548"))(100))
ggsave(as.ggplot(p),file='res_pic/main_figure4/definitive_key_ligand_target_genes_heatmap.pdf',width = 12,height = 3,dpi = 300)

saveRDS(definitive_ligand_target_res_list,file = 'res_data/pd_definitive_ligand_target_res_list.rds')

DotPlot(NRBC_altas_seu,features = rownames(definitive_ligand_target_res_list[[3]]),scale = F)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')

definitive_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(definitive_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(definitive_LR_targetgene_enrichGO_res,file = 'res_data/definitive_LR_targetgene_enrichGO_res.rds')
cnetplot(definitive_LR_targetgene_enrichGO_res,20)+ggtitle('the  ligands of definitive NRBC to target genes')
dotplot(definitive_LR_targetgene_enrichGO_res)
p=draw_gonetcwork_pic_func(res =definitive_LR_targetgene_enrichGO_res,showCategory =20,xlimits = c(-0.2,2) );p
ggsave(p,filename='res_pic/main_figure4/pd_definitive_target_gene_enrichGO_res.pdf',width =6,height = 6)


temp_df=pd_whole_level_markers_degs_df[pd_whole_level_markers_degs_df$gene %in% rownames(definitive_ligand_target_res_list[[3]]) & pd_whole_level_markers_degs_df$avg_log2FC>1,]
#temp_df$gene=factor(temp_df$gene,levels = rownames(fetal_ligand_target_res_list[[3]]))
temp_df=temp_df[temp_df$pct.2<0.1 & temp_df$pct.1>0.08,]
temp_df=temp_df[order(temp_df$avg_log2FC,decreasing = T),]
temp_df$gene=factor(temp_df$gene,levels = temp_df$gene)
p1=ggplot(temp_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()+ggtitle('the hDEGs of the target genes by definitive NRBC ligand')

p2=DotPlot(NRBC_altas_seu,features =as.character(temp_df$gene),scale = F)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
p=p1+p2+plot_layout(heights = c(0.4,1.2),ncol = 1);p
ggsave(p,filename='res_pic/main_figure4/pd_definitive_hDEG_targetgene_expression.pdf',width = 8,height = 12)



ligands_oi <- c('EPO') # this can be a list of multiple ligands  if required
targets_oi <- c( names(filt_primitive_ligand_target_res_list[[3]][,'EPO'][filt_primitive_ligand_target_res_list[[3]][,'EPO'] >0.05]),
                 names(definitive_ligand_target_res_list[[3]][,'EPO'][definitive_ligand_target_res_list[[3]][,'EPO'] >0.05]))
EPO_signaling=signaling_pic_func(ligands_oi =ligands_oi,targets_oi = targets_oi,weighted_networks =weighted_networks ,ligand_tf_matrix =ligand_tf_matrix ,top_n_regulators =4  )
# EPO_signaling_network.pdf,5 X6
key_network_genes[['EPO']]= unique(EPO_signaling$data$label)[!unique(EPO_signaling$data$label) %in% c(targets_oi,ligands_oi)]
targets_gene_list[['EPO']]=targets_oi

saveRDS(key_network_genes,file = 'res_data/key_network_genes.rds')
saveRDS(targets_gene_list,file = 'res_data/targets_gene_list.rds')


p=DotPlot(NRBC_altas_seu,group.by = 'pd_celltype1',cols = c('gray','firebrick3'),scale = F,
        features =unique(c(key_network_genes[['EPO']],key_network_genes[['IGFBP3']],
                           key_network_genes[['MDK']],key_network_genes[['APP']])
                    ) )+RotatedAxis()
p
ggsave(p,filename='res_pic/main_figure4/key_signaling_gene_expression_dotplot.pdf',width =12 ,height = 5)

p=DotPlot(NRBC_altas_seu,group.by = 'pd_celltype1',cols = c('gray','firebrick3'),scale = F,
          features = c('KRT19','HMOX1', 'NDRG1','SOX4','MYB','LMNA','CMBL','ISG20','ITGA2B'))+RotatedAxis()
p
ggsave(p,filename='res_pic/main_figure4/pd_specific_target_gene_expression_dotplot.pdf',width =8 ,height = 5)

#
#基因	主要功能	与EPO信号的关系
#BCL2L1	抗凋亡、细胞存活	已确认的STAT5直接靶基因
#PIM1	细胞存活、增殖	EPO快速响应基因
#CISH	细胞因子信号负反馈	STAT5直接靶基因
#SOCS3	抑制JAK/STAT信号	STAT5直接靶基因
#CLINT1	EPOR内吞及受体周转	EPO/STAT5响应靶基因
#SUV420H2	染色质调控、H4K20甲基化	STAT5直接靶基因
#RBM25	RNA剪接调控	EPO/STAT5响应靶基因
#GYPC	红细胞膜结构	延迟性EPO响应基因
#PODXL	红系细胞膜及细胞结构	延迟性EPO响应基因

EPO_response <- c("BCL2L1","PIM1","CISH","SOCS3","CLINT1","SUV420H2","RBM25")
EPO_response=EPO_response[EPO_response %in% rownames(NRBC_altas_seu)]
VlnPlot(NRBC_altas_seu,features = EPO_response,stack = T,group.by ='source_celltype',split.by = 'tissue_stage',cols = cols)

#EPOR下游核心信号转导基因
#通路	核心基因	主要功能
#JAK2–STAT5	JAK2, STAT5A, STAT5B	EPO信号转导及转录调控
#PI3K–AKT	PIK3CA, PIK3R1, AKT1, AKT2	细胞存活、增殖
#MAPK	GRB2, SOS1, RAF1, MAPK1, MAPK3	增殖及分化
#信号调节	LYN, FYN, SRC, CRKL	调节JAK2/STAT5活化
EPOR_signaling <- c("JAK2","STAT5A","STAT5B","PIK3CA", "PIK3R1","AKT1","AKT2","MAPK1", "MAPK3", "GRB2", "RAF1")
EPOR_signaling=EPOR_signaling[EPOR_signaling %in% rownames(NRBC_altas_seu)]

#AL1和GATA1不仅参与EPO响应，还能直接参与EPOR基因的转录调控；KLF1则主要参与EPO下游红系转录程序的协同调控。
#TAL1	直接结合EPOR调控区域，促进EPOR转录	较强，有启动子结合及敲低实验
#GATA1	直接结合EPOR启动子中的GATA motif，促进EPOR表达	较强，有启动子功能实验
#KLF1	参与EPOR相关的红系转录调控网络，并与STAT5协同调控下游基因	有功能关联，但不是目前证据最明确的EPOR直接调控因子
EPOR_regulator_genes=c('GATA1','TAL1','KLF1')

VlnPlot(NRBC_altas_seu,features = EPOR_signaling,stack = T,group.by ='source_celltype',split.by = 'tissue_stage',cols = cols)
temp_df=AverageExpression(NRBC_altas_seu,features =c(EPOR_regulator_genes,EPO_response,EPOR_signaling),group.by = 'source_celltype' )
temp_df=temp_df$RNA;temp_df=temp_df[c(EPOR_regulator_genes,EPO_response,EPOR_signaling),]
p=pheatmap(t(temp_df),cluster_rows = F,cluster_cols = F,scale = 'column',color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))
ggsave(as.ggplot(p),filename='res_pic/main_figure4/EPO_EPOR_response_genes_heatmap.pdf',width = 8,height = 8)


####################################################################################################################################################################
#----------------------------------------part4: expression of key ligand in altas---------------------------------------#
####################################################################################################################################################################

p=VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features =c( 'CDKN1A'),raster =T,pt.size = 0.5)+NoLegend()
ggsave(p,filename='res_pic/main_figure4/CDKN1A_vlnplot.pdf',width =10,height = 6)

VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features =c( 'CCND2','MYC'),raster =T,pt.size = 0,stack = T,cols = cols)+NoLegend()


cho_feature=c('IGFBP3',"EPO")

YS_altas_seu= readRDS('../NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' )
YS_altas_seu=NormalizeData(YS_altas_seu)
levels(YS_altas_seu$subcelltype)
p1=VlnPlot(subset(YS_altas_seu,subcelltype  %in% c('MESOTHELIUM','SMOOTH_MUSCLE','ENDODERM','ENDOTHELIUM','FIBROBLAST')),
          group.by = 'subcelltype',features = c('IGFBP3'), cols = cols,pt.size = 0)+NoLegend()
p2= VlnPlot(subset(YS_altas_seu,subcelltype %in% c('ENDODERM')),features = c('EPO'),group.by = 'subcelltype', cols = cols,pt.size = 0) +NoLegend()
  

p=(p1|p2)+plot_layout(widths = c(5,1),ncol = 2)
ggsave(p,filename='res_pic/main_figure4/EPO_IGFBP3_YS_celltype_expression_vlnplot.pdf',width = 4,height = 6)



p1=VlnPlot(subset(YS_altas_seu,subcelltype %in% c('MESOTHELIUM','SMOOTH_MUSCLE','ENDODERM','ENDOTHELIUM','FIBROBLAST')),features = c('IGFBP3'), group.by ='stage',cols = cols,pt.size = 0)+
  VlnPlot(subset(YS_altas_seu,subcelltype %in% c('ENDODERM')),features = c('EPO'), group.by ='stage',cols = cols,pt.size = 0)
p1
ggsave(p1,filename='res_pic/main_figure4/EPO_IGFBP3_YS_expression_vlnplot.pdf',width = 6,height = 6)


YS_cho_gene_aggregated_exp=AggregateExpression(YS_altas_seu,group.by = 'id')$RNA[cho_feature,]
YS_cho_gene_aggregated_exp=log2(cpm(YS_cho_gene_aggregated_exp+1))
colnames(YS_cho_gene_aggregated_exp)=paste('YS',colnames(YS_cho_gene_aggregated_exp),sep='_')
saveRDS(YS_cho_gene_aggregated_exp,file='res_data/YS_IGFBP3_EPO_aggregated_exp.rds')
rm(YS_altas_seu);gc()

FL_altas_seu=readRDS('../NRBC_FL_altas/tmp_FL_altas_seu.rds')
p1=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("ENDOTHELIUM","FIBROBLASTS","SMOOTH MUSCLE"  , "SKELETAL MUSCLE" ,"MESOTHELIUM","NEURON")),
           group.by = 'subcelltype',features =c( 'IGFBP3'),pt.size = 0,cols = cols)+NoLegend()
p2=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("HEPATOCYTE")),
           group.by = 'subcelltype',features =c( 'EPO'),pt.size = 0,cols = cols)+NoLegend()

p=(p1|p2) + plot_layout(widths = c(6, 1))
ggsave(p,filename='res_pic/main_figure4/FL_IGFBP3_EPO_expression.pdf',width = 8,height = 5)

FL_altas_seu$age=factor(FL_altas_seu$age,levels = c("CS14_4PCW" ,"CS15_5PCW", "CS17_6PCW", "CS18","CS22","CS23",  "8PCW","8.1PCW","9.1PCW","9.7PCW",  "11PCW","11.4PCW" ,"12PCW","13.9PCW","14.4PCW","15PCW","16.3PCW","16PCW","17PCW"  ))
p1=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("ENDOTHELIUM","FIBROBLASTS","SMOOTH MUSCLE"  , "SKELETAL MUSCLE" ,"MESOTHELIUM","NEURON")),features = c('IGFBP3'), group.by ='age',ncol = 1,cols = cols,pt.size = 0)+NoLegend()
p1
p2=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("HEPATOCYTE")),features = c('EPO'), group.by ='age',ncol = 1,cols = cols,pt.size = 0)+NoLegend()
p2
p=p1/p2;p

ggsave(p,filename='res_pic/main_figure4/FL_IGFBP3_key_ligand_expression_vlnplot.pdf',height = 12,width = 8)
rm(p);gc()


FL_altas_seu$id=FL_altas_seu$donor
FL_altas_seu$id[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]]=paste(FL_altas_seu$donor[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]],
                FL_altas_seu$age[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]],sep='_')
FL_cho_gene_aggregated_exp=AggregateExpression(FL_altas_seu,group.by = 'id')$RNA
FL_cho_gene_aggregated_exp=log2(cpm(FL_cho_gene_aggregated_exp)+1)
colnames(FL_cho_gene_aggregated_exp)=gsub(pattern = 'wk','WPC',colnames(FL_cho_gene_aggregated_exp))

order_sampleid=c( "FL-4WPC" ,"FL-5WPC" ,"FL-6WPC", "F61-CS18","F35-CS22","F32-CS22","F34-CS23", "FL-8WPC"  ,"F16-8.1PCW","F17-9.1PCW","F22-9.7PCW","F33-9.7PCW","FL-11WPC" ,
                  "F23-11.4PCW","F30-14.4PCW","F38-12PCW" , "F45-13.9PCW" ,"F30-14.4PCW", "FL-15WPC"  ,"F41-16PCW","F21-16.3PCW","F29-17PCW")

colnames(FL_cho_gene_aggregated_exp)[!colnames(FL_cho_gene_aggregated_exp) %in% order_sampleid]
FL_cho_gene_aggregated_exp=FL_cho_gene_aggregated_exp[,order_sampleid]
colnames(FL_cho_gene_aggregated_exp)=paste('FL',colnames(FL_cho_gene_aggregated_exp),sep='_')
FL_cho_gene_aggregated_exp=FL_cho_gene_aggregated_exp[cho_feature,]
saveRDS(FL_cho_gene_aggregated_exp,file='res_data/FL_IGFBP3_EPO_aggregated_exp.rds')
rm(FL_altas_seu);gc()


BM_altas_seu=readRDS('../NRBC_BM_altas/BM_altas_seu_v2.rds')#
VlnPlot(BM_altas_seu,group.by = 'age',features = cho_feature,stack = T)

p2=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% c(  "FIBROBLASTS","ENDOTHELIUM", "OSTEOBLAST")),cols=cols,features = 'IGFBP3',group.by = 'new_celltype')+NoLegend()+ggtitle('FBM ALTAS')
p2
ggsave(p2,filename='res_pic/main_figure4/FBM_celltype_IGFBP3_expression_vlnplot.pdf',height =5,width = 6)

p2=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% c(  "FIBROBLASTS")),features = c('IGFBP3'), group.by ='age',ncol = 1,cols = cols,pt.size = 0)+NoLegend()
p2

p3=VlnPlot(subset(subset(BM_altas_seu,stage=='ABM'),new_celltype %in% c("OSTEOBLAST","VSMC" , "ENDOTHELIUM","Adipo-MSC","APOD+ MSC","Fibro-MSC","Osteo-MSC","THY1+ MSC","RNAlo MSC")),
           cols=cols,features = 'IGFBP3',group.by = 'new_celltype')+NoLegend()+ggtitle('ABM ALTAS')
p3
ggsave(p3,filename='res_pic/main_figure4/ABM_celltype_IGFBP3_expression_vlnplot.pdf',height =5,width = 6)


BM_altas_seu$id=paste(BM_altas_seu$age,BM_altas_seu$sample,sep='_')
BM_cho_gene_aggregated_exp=AggregateExpression(BM_altas_seu,group.by ='id' )$RNA
BM_cho_gene_aggregated_exp=log2(cpm(BM_cho_gene_aggregated_exp)+1)

FBM_cho_gene_aggregated_exp_df=BM_cho_gene_aggregated_exp[,grep('PCW',colnames(BM_cho_gene_aggregated_exp))]
ABM_cho_gene_aggregated_exp_df=BM_cho_gene_aggregated_exp[,grep('y',colnames(BM_cho_gene_aggregated_exp))]

colnames(FBM_cho_gene_aggregated_exp_df)=paste('FBM',colnames(FBM_cho_gene_aggregated_exp_df),sep='_')
colnames(FBM_cho_gene_aggregated_exp_df)=gsub(pattern = 'g',colnames(FBM_cho_gene_aggregated_exp_df),replacement = '')
colnames(ABM_cho_gene_aggregated_exp_df)=paste('ABM',colnames(ABM_cho_gene_aggregated_exp_df),sep='_')
colnames(ABM_cho_gene_aggregated_exp_df)=gsub(pattern = 'g',colnames(ABM_cho_gene_aggregated_exp_df),replacement = '')

FBM_cho_gene_aggregated_exp_df=FBM_cho_gene_aggregated_exp_df[cho_feature,]
saveRDS(FBM_cho_gene_aggregated_exp_df,file='res_data/FBM_IGFBP3_EPO_aggregated_exp.rds')
ABM_cho_gene_aggregated_exp_df=ABM_cho_gene_aggregated_exp_df[cho_feature,]
saveRDS(ABM_cho_gene_aggregated_exp_df,file='res_data/ABM_IGFBP3_EPO_aggregated_exp.rds')
rm(BM_altas_seu);gc()


EPO_list=list('YS'=YS_cho_gene_aggregated_exp['IGFBP3',],
                 'FL'=FL_cho_gene_aggregated_exp['EPO',])
                 
EPO_df_igf <- data.frame(
  group = rep(names(EPO_list), times = sapply(EPO_list, length)),
  expr  = unlist(EPO_list, use.names = FALSE)
)

df_n <- EPO_df_igf %>%
  group_by(group) %>%
  summarise(n = n(), y = max(expr) * 1.08, .groups = "drop") %>%
  mutate(label = paste0("n = ", n))

EPO_df_igf$group=factor(EPO_df_igf$group,levels = c('YS','FL'))
p <- ggplot(EPO_df_igf, aes(x = group, y = expr, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.3) +
  geom_text(data = df_n, aes(x = group, y = 0, label = label),  inherit.aes = FALSE, size = 6) +
  theme_classic() +
  labs(x = NULL, y = "log2(CPM+1)", title = "EPO expression")+NoLegend()

# 两两 Wilcoxon 比较（自动加显著性星号）
library(ggpubr)
p1=p + stat_compare_means(
  comparisons = list(c("YS", "FL")),
  method = "t.test",
  p.adjust.method = "BH"
)

p1
ggsave(p3,filename='res_pic/main_figure4/niche_EPO_expression_vlnplot.pdf',height =6,width =2)


IGFBP3_list=list('YS'=YS_cho_gene_aggregated_exp['IGFBP3',],
               'FL'=FL_cho_gene_aggregated_exp['IGFBP3',],
              'FBM'=FBM_cho_gene_aggregated_exp_df['IGFBP3',],
              'ABM'=ABM_cho_gene_aggregated_exp_df['IGFBP3',])
saveRDS(IGFBP3_list,file = 'res_data/niche_exp_IGFBP3_list.rds')

df_igf <- data.frame(
  group = rep(names(IGFBP3_list), times = sapply(IGFBP3_list, length)),
  expr  = unlist(IGFBP3_list, use.names = FALSE)
)

df_igf$group=factor(df_igf$group,levels = c('YS','FL','FBM','ABM'))

df_n <- df_igf %>%
  group_by(group) %>%
  summarise(n = n(), y = max(expr) * 1.08, .groups = "drop") %>%
  mutate(label = paste0("n = ", n))


p <- ggplot(df_igf, aes(x = group, y = expr, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.3) +
  geom_text(data = df_n, aes(x = group, y = 0, label = label),  inherit.aes = FALSE, size = 6) +
  theme_classic() +
  labs(x = NULL, y = "log2(CPM +1) ", title = "IGFBP3 expression")+NoLegend()

# 两两 Wilcoxon 比较（自动加显著性星号）
library(ggpubr)
p2=p + stat_compare_means(
  comparisons = list(c("YS", "FL"), c("YS", "FBM"), c("YS", "ABM"),
                     c("FL", "FBM"), c("FL", "ABM"), c("FBM", "ABM")),
  method = "t.test",
  p.adjust.method = "BH"
)

p2

p=(p1|p2)+plot_layout(widths = c(1,3))
ggsave(p,filename='res_pic/main_figure4/niche_IGFBP3_EPO_expression_compare_ttest.pdf',height =6,width =8)



rm(BM_altas_seu);gc()



##################################################################################################################################################################
# ----------------------------------------查看重点关注通路中的信号流-----------------------------------------------------#c()
##################################################################################################################################################################
YS_subcelltype_ccd=readRDS('../NRBC_YS_altas/YS_subcelltype_cellchat.rds')

pathways.show='IGFBP';LR.show='IGFBP3_TMEM219'
p=netVisual_individual(YS_subcelltype_ccd, sources.use =c(1,2,8:length(YS_subcelltype_ccd@idents)) ,targets.use =3:6 ,signaling = pathways.show, pairLR.use = LR.show, layout = "circle",remove.isolate = T)[[1]]
p
#ggsave(p,filename='res_pic/main_figure4/YS_IGFBP3_interaction.pdf',width = 5,height = 5)

p=netVisual_individual(YS_subcelltype_ccd, sources.use =c(1,2,8:length(YS_subcelltype_ccd@idents)) ,targets.use =3:6 ,signaling = 'EPO', pairLR.use = 'EPO_EPOR', layout = "circle",remove.isolate = T)[[1]]
p
#ggsave(p,filename='res_pic/main_figure4/YS_EPO_interaction.pdf',width = 5,height = 5)

YS_altas_seu= readRDS('../NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' )




