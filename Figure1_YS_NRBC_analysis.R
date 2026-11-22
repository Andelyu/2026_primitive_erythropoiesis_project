
############################################################################################################################### 
#------------------------------------YS Ery altas data analysis, version2 yuchunlai_whu@163.com------------------------------------
############################################################################################################################### 
source('/home/gibh/2021_NRBC_chlyu/zx_lab_NRBC/scripts/scRNAseq_pipline/scRNAseq_analysis_model.R')

library(SeuratDisk)
library(anndata)
library(Seurat)
library(RColorBrewer)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
setwd('/home/gibh/2021_NRBC_chlyu/')
set.seed(1234)


if(F){
  #----------------------- refer gene symbol ------------------------#
  all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  
  all_shared_ensembl_id_info$org_symbol=as.character(mapIds(x = org.Hs.eg.db,keys = all_shared_ensembl_id_info$X,keytype = 'ENSEMBL',column = 'SYMBOL'))
  all_shared_ensembl_id_info$org_symbol[is.na(all_shared_ensembl_id_info$org_symbol)]=all_shared_ensembl_id_info$ref_symbol[is.na(all_shared_ensembl_id_info$org_symbol)]
  table(all_shared_ensembl_id_info$ref_symbol==all_shared_ensembl_id_info$org_symbol) # F:2482,T:33248 /23824
  
  all_shared_ensembl_id_info$ref_symbol1=all_shared_ensembl_id_info$fl_symbol
  all_shared_ensembl_id_info$ref_symbol1[is.na(all_shared_ensembl_id_info$ref_symbol1)]=all_shared_ensembl_id_info$org_symbol[is.na(all_shared_ensembl_id_info$ref_symbol1)]
  table(is.na(all_shared_ensembl_id_info$ref_symbol1))
  table(all_shared_ensembl_id_info$ref_symbol1==all_shared_ensembl_id_info$org_symbol) # F:11906 ,T: 23824
  rownames(all_shared_ensembl_id_info)=all_shared_ensembl_id_info$X
  
  write.csv(all_shared_ensembl_id_info,file = 'ref_data/all_shared_ensembl_id_info.csv')
  
  # ---------------------prepare refer bulk data data---------------#
  load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T)  
  table(rownames(YS_altas_Ery_seu@assays$RNA) %in% rownames(nrbc_ref_se2@assays@data$logcounts)) #FALSE/TRUE :12958 16594 
  
  logcounts=nrbc_ref_se2@assays@data$logcounts
  symbol_ToEn=as.character(mapIds(x = org.Hs.eg.db,keys = rownames(nrbc_ref_se2@assays@data$logcounts),keytype = 'ALIAS',column = 'ENSEMBL' ))
  table(symbol_ToEn %in% all_shared_ensembl_id_info$X) # F/T: 3487 18407
  symbol_ToEn=all_shared_ensembl_id_info[match(x =as.character(symbol_ToEn) ,table = all_shared_ensembl_id_info$X),'ref_symbol1']
  table(is.na(symbol_ToEn))# TRUE:3487, FALSE:18407
  symbol_ToEn[is.na(symbol_ToEn)]=rownames(logcounts)[is.na(symbol_ToEn)]
  du_num=which(symbol_ToEn %in% symbol_ToEn[duplicated(symbol_ToEn)] )
  du_logcounts=logcounts[du_num,]
  du_logcounts=aggregate(du_logcounts,list(symbol_ToEn[du_num]),sum);du_logcounts=column_to_rownames(du_logcounts,'Group.1')
  logcounts=rbind(logcounts[-du_num,],du_logcounts)
  
  pheatmap(cor(logcounts))
  nrbc_ref_se2=SummarizedExperiment(assays =list(logcounts=logcounts),colData = data.frame(nrbc_ref_se2@colData) )
  
  
  logcounts=nrbc_ref_se2@assays@data$logcounts[,-1:-4]
  pheatmap(cor(logcounts))
  nrbc_ref_se=SummarizedExperiment(assays =list(logcounts=logcounts),colData = data.frame(nrbc_ref_se@colData) )
  
  rm(logcounts);gc()
  
  logcounts=nrbc_ref_se2@assays@data$logcounts[,-grep('HSPC|CD34',colnames(nrbc_ref_se2@assays@data$logcounts))]
  colnames(logcounts)=c("BFUE","CFUE",colnames(logcounts)[3:length(colnames(logcounts))])
  colData=data.frame(celltype=as.character(t(data.frame(strsplit(colnames(logcounts),split = '_')))[,1]),
                     sourceis=c(rep('GSE128268',3),rep('GSE53983',15),rep('GSE61566',6)))
  rownames(colData)= colnames(logcounts)
  
  colnames(logcounts)=c(colnames(logcounts)[1:2],'ProE',colnames(logcounts)[4:length(colnames(logcounts))])
  annotation_df=data.frame(row.names =  colnames(logcounts),sourceid=c(rep('GSE128268',3),rep('GSE53983',15),rep('GSE61566',6)))
  p=pheatmap(annotation_row = annotation_df,cor(logcounts),color = colorRampPalette(colors = c('navy','white','firebrick3'))(100),main = 'singleR:refence bulkRNAseq',cluster_cols = F)
  
  # 最后还是剔除HSPC
  
  nrbc_ref_se2=SummarizedExperiment(assays =list(logcounts=logcounts),colData =colData )# 去掉HSCP细胞
  
  logcounts=nrbc_ref_se@assays@data$logcounts[,-grep('HSPC|CD34',colnames(nrbc_ref_se@assays@data$logcounts))]
  nrbc_ref_se=SummarizedExperiment(assays =list(logcounts=logcounts),colData =nrbc_ref_se@colData[-grep('HSPC|CD34',rownames(nrbc_ref_se@colData)),] )# 去掉HSCP细胞
  colnames(nrbc_ref_se@colData)='celltype'
  
  save(nrbc_ref_se,nrbc_ref_se2,file = '/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata')
  
  library(ggplotify)
  ggsave(as.ggplot(p),filename='ref_data/NRBC_ref_bullk_RNAseq_cor_heatmap.pdf',width = 10,height = 8)
  
}else{
  load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T)  
  all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  
}


   
############################################################################################################################### 
#-------------------------------------prepare  data------------------------------------
############################################################################################################################### 

# ----------------------------------- YS  data information ------------------------------#
# 多个数据来源，导致gene symbol并没有统一，需要统一注释不同来源gene symbol
# 这个样本数据也是集合了多篇文章的数据，包括： 
# Popescu_et_al_Nature_2019:Decoding human fetal liver haematopoiesis 的YS：F32,F35(CS23(8 PWC)，74 个Ery,太少，进行细胞比例的时候删掉了), F37，数据在验证集中已经使用
# Wang_et_al_Cell_Stem Cell_2021 ：Decoding Human Megakaryocyte Development， Wang_CS10, Wang_CS11， 验证集中已经使用，
# HM_YS：MAPPING HUMAN HAEMATOPOIETIC STEM CELLS FROM HAEMOGENIC ENDOTHELIUM TO BIRTH ，2022 nature, :mikola_1，只有一个样本:CS14
# WE_YS：F138:CS17，一个样本,
# Mapping the developing human immune system across organs*,Science,"YS_F37_4PCW"  "LI_F61_7PCW"  "YS_F61_7PCW"  "LI_F19_10PCW"；YS_F37_4PCW在上述数据中已经在验证集被使用
# De_novo：F61,F79,F80，三个样本，CS18，CS18，CS15,其中F61在 接下来的数据中已经被使用
# 总结：F79,F80，F138，mikola_1四个样本没有使用过
# YS-FL match sample id: F32,F35,F61 三个样本
# F开头的样本中，同一个样本，有不同区域来源的NRBC， 需要再次优化
library(rhdf5)
library(Matrix)

if(F){
  #------------------------这个YS数据，来自多个文章-----------------#
  Convert(source ='NRBC_YS_altas/raw_ref_data/ys_portal_object_20221208.h5ad', dest="h5seurat" )
  YS_altas_seu=LoadH5Seurat('NRBC_YS_altas/raw_ref_data/ys_portal_object_20221208.h5seurat', meta.data = FALSE)
  meta_df=read.csv('NRBC_YS_altas/raw_ref_data/ys_portal_object_20221208_meta.csv',header = T)
  # 查定，修改F138-->F158
  rownames(meta_df)=meta_df$X;meta_df=meta_df[,-1]
  meta_df$id=paste(meta_df$stage,meta_df$fetal.ids,sep = '_')
  
  YS_altas_seu@meta.data=meta_df
  
  table(YS_altas_seu$LVL1);table(YS_altas_seu$LVL2);table(YS_altas_seu$LVL3)
  t(table(YS_altas_seu@meta.data[,c('LVL1','LVL2')]))
  
  
  # M. Haniffa, Decoding human fetal liver haematopoiesis. Nature574, 365–371 (2019).Dorin-Mirel Popescu 1，E-MTAB-7407
  # M. Haniffa, S. A. Teichmann, Mapping the developing human immune system across organs. Science376, eabo0510 (2022).Chenqu Suo 1，E-MTAB-11343
  # H. K. A. Mikkola, Mapping human haematopoietic stem cells from haemogenic endothelium to birth. Nature604, 534–540 (2022).Vincenzo Calvanese 1， GSE162950. Data from published：GSE135202(无feature 信息).
  # Decoding Human Megakaryocyte Development. Cell Stem Cell28, 535–549.e8 (2021). GSE144024
  # The role of the yolk sac in human fetal development and identification of a hepatocyte-like cell in the human yolk sac: E-MTAB-10552, Michael Mather 1Muzz Haniffa 
  # S. Webb, M. Haniffa, E. Stephenson, “Human fetal yolk sac scRNA-seq data (sample ID: F158 for Haniffa Lab; 16099 for HDBR),” BioStudies, E-MTAB-11673 (2022);
  # We :(CellRanger version 6.0.1 using human reference genome GRCh38-2020-A ,https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11673
  
  # Smart2-seq M. Haniffa, M. Mather, R. Botting, “The role of the yolk sac in human fetal development and identification of a hepatocyte-like cell in the human yolk sac (SS2),” BioStudies, E-MTAB-10888 (2023)
  
  c('GNB2L1','ATP5E','RACK1','ATP5F1E') %in% rownames(YS_altas_seu)
  object_list=SplitObject(YS_altas_seu,split.by = 'orig.dataset')
  object_list[[1]]=CreateSeuratObject(counts = GetAssayData(object_list[[1]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data =object_list[[1]]@meta.data );dim(object_list[[1]])
  object_list[[2]]=CreateSeuratObject(counts = GetAssayData(object_list[[2]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[2]]@meta.data);dim(object_list[[2]])
  object_list[[3]]=CreateSeuratObject(counts = GetAssayData(object_list[[3]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[3]]@meta.data);dim(object_list[[3]])
  object_list[[4]]=CreateSeuratObject(counts = GetAssayData(object_list[[4]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[4]]@meta.data);dim(object_list[[4]])
  object_list[[5]]=CreateSeuratObject(counts = GetAssayData(object_list[[5]],assay = 'RNA'),min.cells = 10,min.features = 200,meta.data = object_list[[5]]@meta.data);dim(object_list[[5]])
  
  all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  rownames(all_shared_ensembl_id_info)=all_shared_ensembl_id_info$X
  all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  table(duplicated(all_shared_ensembl_id_info$ref_symbol1)) #497 symbol replicated 
  
  
  GSE162950_feature_info=read.table('NRBC_YS_altas/raw_ref_data/GSE162950_features.tsv',sep="\t")
  YS_feature_info=read.table('NRBC_YS_altas/YS_10x_genes.tsv',sep="\t",header = T)
  
  We_feature_info=read.table('NRBC_YS_altas/raw_ref_data/gene_info.txt',sep="\t",header = F) # cellranger_genome_2020A
  We_feature_info=We_feature_info[,-1];We_feature_info$V2=gsub(pattern = ' ',replacement = '',We_feature_info$V2)
  table(We_feature_info$V2 %in% all_shared_ensembl_id_info$X)# F/T:3695/32905 
  table(We_feature_info[!We_feature_info$V2 %in% all_shared_ensembl_id_info$ref_symbol1,'V3'])
  #lncRNA  protein_coding
  # 16561            19394  
  
  names(object_list)
  #  Wang_et_al_Cell_Stem_Cell_2021:zx_lab,Popescu_et_al_Nature_2019 & De_novo: FL , Popescu_跟FL数据是同一来源, HM_YS、 WE_YS注释信息各不同, 需要找到原文数据
  # 以 Popescu_et_al_Nature_2019 数据注释为参考，统一全文，需要重新注解
  
  # 统计不同来源数据注释基因差异
  # 发现CS10_Wang_CS10、CS11_Wang_CS11、CS14_mikola_1、CS17_F158: Wang_et_al_Cell_Stem_Cell_2021: zx_lab,
  # Popescu_et_al_Nature_2019+De_novo+WE_YS:gene feature 注释信息比较一致，F*样本的gene feature 注释信息比较一致，
  ## 5ad文件无feature信息, Wang_et_al_Cell_Stem: GSE144024
  
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$abm_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$zxlab_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$fl_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$ref_symbol)})
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$ref_symbol1)})
  sapply(object_list, function(x){table(rownames(x) %in% GSE162950_feature_info$V2)}) # 同Wang_et_al_Cell_Stem
  sapply(object_list, function(x){table(rownames(x) %in% as.character(YS_feature_info$X))}) # WE_YS:2174
  sapply(object_list, function(x){table(rownames(x) %in% We_feature_info$V4)}) # WE_YS:F/T:9/23814
  
  
  new_symbol_list=list()
  refId="Wang_et_al_Cell_Stem Cell_2021"
  new_symbol_list[[refId]]=all_shared_ensembl_id_info[match( x = rownames(object_list[[refId]]),table = all_shared_ensembl_id_info$zxlab_symbol),'ref_symbol1']
  out_re_symbols=rownames(object_list[[refId]])[! rownames(object_list[[refId]]) %in% all_shared_ensembl_id_info$zxlab_symbol ];out_re_symbols
  out_re_symbols=as.character(t(data.frame(strsplit(out_re_symbols,split = '.',fixed = T)))[,1]);out_re_symbols
  new_symbol_list[[refId]][! rownames(object_list[[refId]]) %in% all_shared_ensembl_id_info$zxlab_symbol ]=out_re_symbols
  
  
  refId="WE_YS"
  new_symbol_list[[refId]]=all_shared_ensembl_id_info[We_feature_info$V2[match( x = rownames(object_list[[refId]]),table = We_feature_info$V4)],'ref_symbol1']
  table(is.na(new_symbol_list[[refId]])) # F/T:22217/1606, lncRNA 占绝大多数，可能跟其测序方法有关
  new_symbol_list[[refId]][is.na(new_symbol_list[[refId]])]=rownames(object_list[[refId]])[is.na(new_symbol_list[[refId]])] 
  out_re_symbols=rownames(object_list[[refId]])[! rownames(object_list[[refId]]) %in% We_feature_info$V4 ];out_re_symbols
  out_re_symbols=as.character(t(data.frame(strsplit(out_re_symbols,split = '-1',fixed = T)))[,1]);out_re_symbols
  out_re_symbols[!out_re_symbols %in% all_shared_ensembl_id_info$ref_symbol1] # 3 不在
  trans_out_re_symbols=all_shared_ensembl_id_info[We_feature_info$V2[match( x = out_re_symbols,table = We_feature_info$V4)],'ref_symbol1']
  trans_out_re_symbols[is.na(trans_out_re_symbols)]=out_re_symbols[is.na(trans_out_re_symbols)]
  new_symbol_list[[refId]][! rownames(object_list[[refId]]) %in% We_feature_info$V4 ]=trans_out_re_symbols
  
  HM_to_refsymbol=mapIds(x = org.Hs.eg.db,keys =rownames( object_list[['HM_YS']]),keytype = 'ALIAS',column = 'ENSEMBL')
  HM_to_refsymbol=all_shared_ensembl_id_info[match(x = as.character(HM_to_refsymbol),table =all_shared_ensembl_id_info$X),'ref_symbol1']
  HM_to_refsymbol[is.na(HM_to_refsymbol)]=rownames( object_list[['HM_YS']])[is.na(HM_to_refsymbol)]
  new_symbol_list[['HM_YS']]=HM_to_refsymbol;rm(HM_to_refsymbol)
  table(new_symbol_list[['HM_YS']] %in% all_shared_ensembl_id_info$ref_symbol1) # F/T:368/16301
  
  new_object_list=list();library(tidyverse)
  for(refId in names(new_symbol_list)){
    print(refId)
    du_symbol=new_symbol_list[[refId]][duplicated(new_symbol_list[[refId]])];
    if( length(du_symbol) == 0){print(refId);next}
    tmp_assay=GetAssayData(object_list[[refId]],assay = 'RNA')
    rownames(tmp_assay)=new_symbol_list[[refId]]
    du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% du_symbol,]
    tmp_assay=tmp_assay[!rownames(tmp_assay) %in% du_symbol,]
    
    du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum) 
    du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
    du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
    
    tmp_assay=rbind(tmp_assay,du_tmp_assay)
    object_list[[refId]]=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =object_list[[refId]]@meta.data )
  }
  
  rm(du_tmp_assay,tmp_assay);gc()
  
  sapply(object_list, function(x){table(rownames(x) %in% all_shared_ensembl_id_info$ref_symbol1)})
  #Wang_et_al_Cell_Stem Cell_2021 Popescu_et_al_Nature_2019 De_novo HM_YS WE_YS
  #FALSE                              3                        13      21   371  1556
  #TRUE                           18397                     19716   23062 16032 22240
  
  YS_altas_seu_umap=YS_altas_seu[['umap']]
  YS_altas_seu=merge(object_list[[1]],c(object_list[[2]],object_list[[3]],object_list[[4]],object_list[[5]]))
  YS_altas_seu[['RNA']]=JoinLayers(YS_altas_seu[['RNA']])
  table(rownames(YS_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol1)# F/T:1939/24374 
  
  YS_altas_seu[['umap']]= YS_altas_seu_umap
  saveRDS(YS_altas_seu,file = 'NRBC_YS_altas/raw_ref_data/YS_altas_seu.rds')
  rm(object_list,YS_altas_seu_umap);gc()
  
  
  DimPlot(YS_altas_seu,group.by = c('orig.dataset','stage'),cols = col,raster=FALSE)
  DimPlot(YS_altas_seu,group.by = c('LVL1','LVL2'),cols = col,raster=FALSE)
  
  
  # 双细胞分析
  YS_altas_seu[["RNA"]]=split(YS_altas_seu[['RNA']], f= YS_altas_seu$id)
  # run standard analysis workflow
  YS_altas_seu <- NormalizeData(YS_altas_seu)
  YS_altas_seu <- FindVariableFeatures(YS_altas_seu)
  YS_altas_seu <- ScaleData(YS_altas_seu)
  YS_altas_seu=RunPCA(YS_altas_seu)
  YS_altas_seu[['RNA']]=JoinLayers(YS_altas_seu[['RNA']])
  
  
  library(DoubletFinder)
  library(scDblFinder)
  
  sce <- as.SingleCellExperiment(YS_altas_seu)
  sce <- scDblFinder(sce,clusters = 'LVL2',samples = 'id',dims = 30,dbr = 0.06,threshold = TRUE)
  
  YS_altas_seu$scDblFinder.class <- sce$scDblFinder.class
  YS_altas_seu$scDblFinder.score <- sce$scDblFinder.score
  saveRDS(YS_altas_seu@meta.data,file = 'tmp_YS_altas_seu_metadat.rds')
  
  sce <- scDblFinder(sce,clusters = 'LVL2',samples = 'lanes',dims = 30,dbr = 0.06,threshold = TRUE)
  YS_altas_seu$scDblFinder.class2 <- sce$scDblFinder.class
  YS_altas_seu$scDblFinder.score2 <- sce$scDblFinder.score
  saveRDS(YS_altas_seu,file = '../NRBC_YS_altas/raw_ref_data/YS_altas_seu.rds')
  YS_altas_seu_metadata=YS_altas_seu@meta.data
  rm(YS_altas_seu)
  rm(sce);
  gc()
  
}else{
  YS_altas_seu=readRDS('NRBC_YS_altas/raw_ref_data/YS_altas_seu.rds')
}

if(F)else{
  load(file = '../ref_data/all_NI_Ery_meta.Rdata',verbose = T)
  CS10_YS=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/CS10_YS.loom');CS10_YS$orig.ident='CS10_YS';head(CS10_YS) # 与zx_lab feuature inf 高度重合，
  CS11_YS=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/CS11_YS.loom');CS11_YS$orig.ident='CS11_YS';head(CS11_YS)
  CS15_YS=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/CS15_YS.loom');CS15_YS$orig.ident='CS15_YS';head(CS15_YS)
  CS11_YS=RenameCells(CS11_YS,new.names = gsub(pattern = 'YS_R2',replacement = 'YS',rownames(CS11_YS@meta.data)))
  CS15_YS=RenameCells(CS15_YS,new.names = gsub(pattern = 'YS_R1',replacement = 'YS',rownames(CS15_YS@meta.data)))
  
  fetal_YS=merge(CS10_YS,c(CS11_YS,CS15_YS));
  fetal_YS$type='YS'
  
  dim(fetal_YS)
  rm(list = c('CS10_YS','CS11_YS','CS15_YS'))
  
  fetal_YS[['percent.mt']] <- PercentageFeatureSet(fetal_YS,pattern = '^MT-' )
  fetal_YS$GS_score=fetal_YS$G2M.Score-fetal_YS$S.Score
  
  all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  rownames(all_shared_ensembl_id_info)=all_shared_ensembl_id_info$X
  all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
  table(duplicated(all_shared_ensembl_id_info$ref_symbol1)) #497 symbol replicated 
  
  table(rownames(fetal_YS) %in% all_shared_ensembl_id_info$zxlab_symbol)# F/T:24/ 33514
  table(rownames(fetal_YS) %in% all_shared_ensembl_id_info$ref_symbol)# F/T:24/ 33514
  
  
  # 先整理不在本身注释基因组信息中，因为添加了后缀
  symbol_addpostfix=rownames(fetal_YS)[!rownames(fetal_YS) %in% all_shared_ensembl_id_info$ref_symbol]
  symbol_addpostfix=data.frame(strsplit2(symbol_addpostfix,split = '.',fix=T))[,1] 
  symbol_addpostfix_inf=all_shared_ensembl_id_info[all_shared_ensembl_id_info$ref_symbol %in% symbol_addpostfix,]
  symbol_addpostfix_inf=symbol_addpostfix_inf[order(symbol_addpostfix_inf$ref_symbol),]
  symbol_addpostfix_inf$number=1:dim(symbol_addpostfix_inf)[1]
  du_symbol_addpostfix=symbol_addpostfix_inf[match(symbol_addpostfix,symbol_addpostfix_inf$ref_symbol)+1,'ref_symbol1']
  names(du_symbol_addpostfix)=paste(symbol_addpostfix,'1',sep = '.')
  
  # 再处理在范围内的symbol
  uniq_symbols=rownames(fetal_YS)[rownames(fetal_YS) %in% all_shared_ensembl_id_info$ref_symbol]
  new_symbol_list=all_shared_ensembl_id_info[match(uniq_symbols,all_shared_ensembl_id_info$ref_symbol),'ref_symbol1']
  names(new_symbol_list)=uniq_symbols
  new_symbol_list=c(du_symbol_addpostfix,new_symbol_list)
  
  # 得到总list
  du_symbol=as.character(new_symbol_list)[duplicated(as.character(new_symbol_list))];length(du_symbol)
  tmp_assay=GetAssayData(fetal_YS,assay = 'RNA',layer = 'counts')
  rownames(tmp_assay)=as.character(new_symbol_list[rownames(tmp_assay)])
  # 将symbol进行映射处理
  du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% du_symbol,]
  tmp_assay=tmp_assay[!rownames(tmp_assay) %in% du_symbol,]
  
  # 再处理映射后重复symbol 
  du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum) 
  du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
  du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
  
  #合并得到最后symbol处理矩阵
  tmp_assay=rbind(tmp_assay,du_tmp_assay)
  fetal_YS=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =fetal_YS@meta.data )
  rm(tmp_assay,du_tmp_assay)
  
  fetal_YS[['RNA']]=split(fetal_YS[['RNA']],f = fetal_YS$orig.ident) 
  fetal_YS=NormalizeData(fetal_YS)
  fetal_YS=FindVariableFeatures(fetal_YS)
  fetal_YS=ScaleData(fetal_YS,vars.to.regress = c('percent.mt','nCount_RNA',''))
  fetal_YS=RunPCA(fetal_YS)
  ElbowPlot(fetal_YS,ndims = 50)
  fetal_YS=IntegrateLayers(fetal_YS,method = CCAIntegration,orig.reduction = 'pca',new.reduction='cca')
  fetal_YS=RunUMAP(fetal_YS,dims = 1:30,reduction = 'cca',reduction.name = 'cca_umap')
  DimPlot(fetal_YS,cols = col,group.by = 'orig.ident',reduction = 'cca_umap')+FeaturePlot(fetal_YS,features = 'GYPA',reduction = 'cca_umap')
  
  fetal_YS=IntegrateLayers(fetal_YS,method = HarmonyIntegration,orig.reduction = 'pca',new.reduction='harmony')
  fetal_YS=RunUMAP(fetal_YS,dims = 1:30,reduction = 'harmony')
  fetal_YS[['RNA']]=JoinLayers(fetal_YS[['RNA']])
  fetal_YS=FindNeighbors(fetal_YS,reduction = 'harmony',dims = 1:30);fetal_YS=FindClusters(fetal_YS,resolution = 0.2)
  DimPlot(fetal_YS,reduction = 'umap',group.by = 'seurat_clusters',cols = cols,label = T)+
    FeaturePlot(fetal_YS,reduction = 'umap',features = c('GYPA'))
  
  fetal_YS$celltype='others'
  fetal_YS@meta.data[fetal_YS$seurat_clusters %in% c('0','3') ,'celltype']='Ery'
  DimPlot(fetal_YS,cols = col,group.by = c('orig.ident','celltype','seurat_clusters'),reduction = 'umap')
    FeaturePlot(fetal_YS,features =c ('GYPA','SOX6','CD69','CD14','SPP1','CD34'),reduction = 'umap')
  
  library(DoubletFinder)
  library(scDblFinder)
    
  sce <- as.SingleCellExperiment(fetal_YS)
  sce <- scDblFinder(sce,clusters = 'seurat_clusters',samples = 'orig.ident',dims = 30,dbr = 0.06,threshold = TRUE)
    
  fetal_YS$scDblFinder.class <- sce$scDblFinder.class
  fetal_YS$scDblFinder.score <- sce$scDblFinder.score
  saveRDS(fetal_YS,file = 'ref_data/NI_fetal_YS_test.rds')
  
  
  fetal_YS_Ery_seu=subset(fetal_YS,celltype=='Ery' & fetal_YS$scDblFinder.class=='singlet' )
  fetal_YS_Ery_seu=RunUMAP(fetal_YS_Ery_seu,reduction = 'harmony',dims = 1:10,reduction.name = 'umap1')
  DimPlot(fetal_YS_Ery_seu,cols = col,group.by = c('orig.ident'),reduction = 'umap1')
  
  FeaturePlot(fetal_YS_Ery_seu,features = c('GYPA','HBE1','CD63','SPP1','HBB'),reduction = 'umap1')
  fetal_YS_Ery_seu$id=paste('NI_',fetal_YS_Ery_seu$orig.ident) 
   DimPlot(fetal_YS_Ery_seu,cols = col,group.by = c('orig.ident'),reduction = 'umap')
  
  fetal_YS_Ery_seu=singleR_analysis_func(refdata =nrbc_ref_se2,test =as.matrix(fetal_YS_Ery_seu@assays$RNA$data)[VariableFeatures(fetal_YS_Ery_seu),],outdata = fetal_YS_Ery_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype')
  fetal_YS_Ery_seu=fetal_YS_Ery_seu[[1]]
  fetal_YS_Ery_seu=singleR_analysis_func(refdata =nrbc_ref_se,test =as.matrix(fetal_YS_Ery_seu@assays$RNA$data)[VariableFeatures(fetal_YS_Ery_seu),],outdata = fetal_YS_Ery_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype2')
  fetal_YS_Ery_seu=fetal_YS_Ery_seu[[1]]
  table(is.na(fetal_YS_Ery_seu$Pre_celltype2))
  table(is.na(fetal_YS_Ery_seu$Pre_celltype))
  
  fetal_YS_Ery_seu$Pre_celltype=factor(fetal_YS_Ery_seu$Pre_celltype,levels = c( "BFUE","CFUE","ProE" ,"eBas" ,"lBas" ,"Poly","Orth"))
  fetal_YS_Ery_seu$Pre_celltype2=factor(fetal_YS_Ery_seu$Pre_celltype2,levels = c( "BFUE","CFUE","ProE" ,"eBas" ,"lBas" ,"Poly","Orth"))
  
  DimPlot(fetal_YS_Ery_seu,cols = col,group.by = c('orig.ident','Pre_celltype','Pre_celltype2'),reduction = 'umap1')
   
  fetal_YS_Ery_seu@meta.data=fetal_YS_Ery_seu@meta.data[,-grep('X_',colnames(fetal_YS_Ery_seu@meta.data))]
  
   
  fetal_YS_Ery_seu@meta.data$id=paste0('NI_',fetal_YS_Ery_seu@meta.data$orig.ident)
  colnames( fetal_YS_Ery_seu@meta.data)=gsub(pattern = 'old.ident',replacement = 'stage',x = colnames( fetal_YS_Ery_seu@meta.data))
  fetal_YS_Ery_seu$orig.dataset='Xu_et_al_Nature_immune_2022'
  saveRDS(fetal_YS_Ery_seu,file = 'primitive_NRBC_project/NI_ref_YS_ery_seu.rds')
 
}else{
  fetal_YS_Ery_seu=readRDS('primitive_NRBC_project/NI_ref_YS_ery_seu.rds')
}

# CS11,46NRBC，太少，CS15有全局性巨噬细胞RNA轻微污染,舍弃，GSE137010 舍弃
if(F){
  CS11_ys_seu=Read10X(data.dir = 'primitive_NRBC_project/raw_data/GSE137010/CS11_YS')
  CS15_ys_seu=Read10X(data.dir = 'primitive_NRBC_project/raw_data/GSE137010/CS15_YS')
  CS11_ys_seu=CreateSeuratObject(counts = CS11_ys_seu,assay = 'RNA',project = 'YS_CS11',min.cells = 10,min.features = 200)
  CS15_ys_seu=CreateSeuratObject(counts = CS15_ys_seu,assay = 'RNA',project = 'YS_CS15',min.cells = 10,min.features = 200)
  YS_seu3=merge(CS11_ys_seu,CS15_ys_seu)
  YS_seu3$orig.dataset='GSE137010'
  YS_seu3$id=paste('Bian',YS_seu3$orig.ident,sep = '_')
  YS_seu3$stage='CS11'
  YS_seu3@meta.data[YS_seu3$orig.ident=='YS_CS15','stage']='CS15'
  
  
  YS_seu3=NormalizeData(YS_seu3) %>% FindVariableFeatures()  %>% ScaleData()
  YS_seu3=RunPCA(YS_seu3)
  ElbowPlot(YS_seu3,ndims = 50)
  
  YS_seu3=RunHarmony(YS_seu3,group.by.vars='id')
  YS_seu3=RunUMAP(YS_seu3,dims = 1:15,reduction = 'harmony')
  DimPlot(YS_seu3,cols = col,group.by = c('orig.ident'),reduction = 'umap')+
    FeaturePlot(YS_seu3,features = c('GYPA','HBE1','CD69','CD14','TNF','SPP1','HBB','EPO','COL3A1'),reduction = 'umap')
  YS_seu3=FindNeighbors(YS_seu3,reduction = 'harmony',dims = 1:30)
  YS_seu3=FindClusters(YS_seu3,resolution = 0.3)
  
  DimPlot(YS_seu3,cols = col,group.by = c('seurat_clusters','id'),reduction = 'umap',label = T)
  YS_seu3$type='other'
  YS_seu3@meta.data[YS_seu3$seurat_clusters %in% c('6','7','10'),'type']='Ery' # 主要来自CS15, NRBC 同样存在Mac污染：SPP1，CS14
  YS_Ery_seu3=subset(YS_seu3,type=='Ery')
  table(YS_Ery_seu3$id)
 #Bian_YS_CS11 Bian_YS_CS15 
  #46         1931 
  DimPlot(YS_Ery_seu3,cols = col,group.by = c('seurat_clusters','id'),reduction = 'umap',label = T)
  YS_Ery_seu3=RunUMAP(YS_Ery_seu3,reduction = 'harmony',dims = 1:10,reduction.name = 'umap1')
  DimPlot(YS_Ery_seu3,cols = col,group.by = c('seurat_clusters','id'),reduction = 'umap1',label = T)
  FeaturePlot(YS_Ery_seu3,features = c('GYPA','HBE1','CD69','CD14','SPP1','HBB','COL3A1'),reduction = 'umap1')
  
  library(decontX)
  YS_seu3[['RNA']]=JoinLayers(YS_seu3[['RNA']])
  raw_counts <- GetAssayData(object = YS_seu3, assay = "RNA", layer = "counts")
  # 准备基因和细胞标签 (传入聚类标签能显著提高计算精度)
  # 如果有明确的批次信息，请将 batch 也加入 z 参数 (选项有: 细胞聚类标签, 或仅用基因名)
  z_labels <- YS_seu3$seurat_clusters  # 这里的标签使用你之前的聚类结果
  
  # 运行 DecontX
  decontX_res <- decontX(counts = raw_counts, z = z_labels,batch ='id')
  YS_seu3$contamination <- decontX_res$contamination
  
  saveRDS(YS_seu3,file = 'raw_data/GSE137010_YS_seu3.rds')
  
  saveRDS(YS_Ery_seu3,file = 'raw_data/GSE137010_YS_Ery_seu3.rds')
  
  FeaturePlot(subset(YS_seu3,stage=='CS15'),features = c('SPP1','CD14'),reduction ='umap') #有全局性巨噬细胞RNA轻微污染,舍弃
}


setwd('primitive_NRBC_project')

#----------- QC ------------#
YS_altas_seu=NormalizeData(YS_altas_seu)
FeaturePlot(subset(YS_altas_seu,id=='CS14_F37'),features = c('SPP1','CD14'),reduction ='umap') #有全局性巨噬细胞RNA 验证污染，剔除该样本
YS_altas_seu=subset(YS_altas_seu,id!='CS14_F37')

# Ery
YS_altas_Ery_seu=subset(YS_altas_seu,LVL1 %in% 'ERYTHROID'  );rm(YS_altas_seu);gc() 
YS_altas_Ery_seu=subset(YS_altas_Ery_seu,scDblFinder.class2=='singlet' & scDblFinder.class2=='singlet')

YS_altas_Ery_seu=NormalizeData(YS_altas_Ery_seu) %>%FindVariableFeatures() %>% ScaleData() %>%RunPCA()
YS_altas_Ery_seu=CellCycleScoring(YS_altas_Ery_seu,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes,set.ident = F)

#------------------------singleR NRBC 亚类注释----------------------------#
# 自动化注释效果在bas，提示primitive NRBC 与defintive NRBC 转录组层面差异
YS_altas_Ery_seu=FindVariableFeatures(YS_altas_Ery_seu) %>% ScaleData() %>% RunPCA()
YS_altas_Ery_seu=singleR_analysis_func(refdata =nrbc_ref_se2,test =as.matrix(YS_altas_Ery_seu@assays$RNA$data)[VariableFeatures(YS_altas_Ery_seu),],outdata = YS_altas_Ery_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype')
YS_altas_Ery_seu=YS_altas_Ery_seu[[1]]
YS_altas_Ery_seu=singleR_analysis_func(refdata =nrbc_ref_se,test =as.matrix(YS_altas_Ery_seu@assays$RNA$data)[VariableFeatures(YS_altas_Ery_seu),],outdata = YS_altas_Ery_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype2')
YS_altas_Ery_seu=YS_altas_Ery_seu[[1]]

table(is.na(YS_altas_Ery_seu$Pre_celltype)) # F/T:60851 
table(is.na(YS_altas_Ery_seu$Pre_celltype2)) # F/T:60851
DimPlot(YS_altas_Ery_seu,group.by = c('Pre_celltype','Pre_celltype2'),cols = cols)


YS_altas_Ery_seu$Pre_celltype=factor(YS_altas_Ery_seu$Pre_celltype,levels = c( "BFUE","CFUE","ProE" ,"eBas" ,"lBas" ,"Poly","Orth"))
YS_altas_Ery_seu$Pre_celltype2=factor(YS_altas_Ery_seu$Pre_celltype2,levels = c( "BFUE","CFUE","ProE" ,"eBas" ,"lBas" ,"Poly","Orth"))

colnames(fetal_YS_Ery_seu@meta.data)
colnames(YS_altas_Ery_seu@meta.data)
fetal_YS_Ery_seu$lanes=fetal_YS_Ery_seu$id
cho_cols=intersect(colnames(YS_altas_Ery_seu@meta.data),colnames(fetal_YS_Ery_seu@meta.data))

YS_altas_Ery_seu=merge(YS_altas_Ery_seu,fetal_YS_Ery_seu)
YS_altas_Ery_seu@meta.data=YS_altas_Ery_seu@meta.data[,cho_cols]
saveRDS(YS_altas_Ery_seu,file = 'new_YS_altas_Ery_seu.rds')


################################  #####################################################################################################################
#----------------------------------defined the YS_NRBC subcelltype --------------------------------------------------------#
##################################################################################################################################################### 
#-----------------------------------------CS10,CS15, CS18 as reference umap --------------------------------#
#-----------------------------------------test1 --------------------------------#
table(YS_altas_Ery_seu$id)
YS_altas_Ery_seu[['RNA']]=JoinLayers(YS_altas_Ery_seu[['RNA']])
YS_altas_Ery_seu =RunPCA(YS_altas_Ery_seu)
YS_altas_refer_NRBC_seu=subset(YS_altas_Ery_seu,id %in% c('CS18_F61','CS15_F80','CS10_Wang_CS10')) #
ElbowPlot(YS_altas_refer_NRBC_seu)
YS_altas_refer_NRBC_seu=FindVariableFeatures(YS_altas_refer_NRBC_seu) %>%ScaleData() %>%RunPCA() %>%RunHarmony(group.by.vars=c('id'))

# test pca number
YS_altas_refer_NRBC_seu=RunUMAP(YS_altas_refer_NRBC_seu,reduction = 'harmony',dims = 1:15,reduction.name ='umap_pc15' )
YS_altas_refer_NRBC_seu=RunUMAP(YS_altas_refer_NRBC_seu,reduction = 'harmony',dims = 1:8,reduction.name ='umap_pc8')

DimPlot(YS_altas_refer_NRBC_seu,reduction = 'umap_pc15',group.by = 'Pre_celltype',cols=cols,label = T)+
  DimPlot(YS_altas_refer_NRBC_seu,reduction = 'umap_pc8',group.by = 'Pre_celltype',cols=cols,label = T)

YS_altas_refer_NRBC_seu=RunUMAP(YS_altas_refer_NRBC_seu,reduction = 'harmony',dims = 1:15,reduction.name ='umap1',return.model = T,spread =0.5 ,seed.use = 42 )
DimPlot(YS_altas_refer_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('id','Pre_celltype','celltype'))

p=DimPlot(YS_altas_refer_NRBC_seu,group.by = c('id','Pre_celltype2'),cols = col,reduction = 'umap1',label = F);p
ggsave(p,filename='res_pic/main_figure1/raw_YS_altas_refer_NRBC_sample_umap.pdf',width = 12,height = 5)


YS_altas_refer_NRBC_seu=FindNeighbors(YS_altas_refer_NRBC_seu,reduction = 'harmony',dims = 1:30)
YS_altas_refer_NRBC_seu=FindClusters(YS_altas_refer_NRBC_seu,resolution = c(0.3,0.2,0.1))
DimPlot(YS_altas_refer_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('RNA_snn_res.0.3','RNA_snn_res.0.2','RNA_snn_res.0.1'))


celltype_marker=c("TFRC","GYPA",'HBE1','HBB',"MALAT1",'MKI67')
DimPlot(YS_altas_refer_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('RNA_snn_res.0.3'))+
  VlnPlot(YS_altas_refer_NRBC_seu,features = celltype_marker,group.by = 'RNA_snn_res.0.3',stack = T)+NoLegend()

Idents(YS_altas_refer_NRBC_seu)='RNA_snn_res.0.3'
YS_altas_refer_NRBC_seu=FindSubCluster(YS_altas_refer_NRBC_seu,cluster = '1',graph.name = 'RNA_snn',resolution = 0.2,subcluster.name = 'sub.cluster')
DimPlot(YS_altas_refer_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('sub.cluster'))

DimPlot(YS_altas_refer_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('sub.cluster'))+
  VlnPlot(YS_altas_refer_NRBC_seu,features = celltype_marker,group.by = 'sub.cluster',stack = T)+NoLegend()

# define the celltype
YS_altas_refer_NRBC_seu$celltype='Bas'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('7'),'celltype']='EryP-CFC'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('4','2'),'celltype']='eBas'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('3'),'celltype']='mBas'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('0'),'celltype']='lBas'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('1_0'),'celltype']='Poly'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('1_1'),'celltype']='Orth'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('8'),'celltype']='def_like_Ery'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('5'),'celltype']='pEry1'
YS_altas_refer_NRBC_seu@meta.data[YS_altas_refer_NRBC_seu@meta.data$sub.cluster %in% c('6'),'celltype']='pEry2'

YS_altas_refer_NRBC_seu$celltype=factor(YS_altas_refer_NRBC_seu$celltype,levels = c('EryP-CFC','eBas','mBas','lBas','Poly','Orth','def_like_Ery','pEry1','pEry2'))

p=DimPlot(YS_altas_refer_NRBC_seu,group.by =c('sub.cluster','celltype'),cols = col,reduction = 'umap1',label = F);p
ggsave(p,filename='res_pic/main_figure1/raw_YS_altas_refer_celltype_NRBC_sample_umap.pdf',width = 12,height = 5)

celltype_marker=c("TFRC","GYPA",'HBE1','HBB','MKI67')
p=VlnPlot(YS_altas_refer_NRBC_seu,features = celltype_marker,group.by = 'celltype',stack = T,cols = cols)+NoLegend();p
ggsave(p,filename='res_pic/main_figure1/raw_YS_altas_refer_celltype_maker_vlnplot.pdf',width = 6,height = 5)


saveRDS(YS_altas_refer_NRBC_seu,file = 'YS_altas_refer_NRBC_seu.rds')

Other_YS_altas_Ery_seu=subset(YS_altas_Ery_seu,!id %in% unique(YS_altas_refer_NRBC_seu$id))
anchors=FindTransferAnchors(reference =YS_altas_refer_NRBC_seu,query = Other_YS_altas_Ery_seu, reference.reduction='pca')
Other_YS_altas_Ery_seu=MapQuery(anchorset =anchors,query = Other_YS_altas_Ery_seu,reference =YS_altas_refer_NRBC_seu,reference.reduction = 'pca',
                                reduction.model = 'umap1',refdata =  list(celltype = "sub.cluster")  )
DimPlot(Other_YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,group.by = c('predicted.celltype'))

YS_altas_Ery_seu[['ref.pca']]=merge(CreateDimReducObject(embeddings = YS_altas_refer_NRBC_seu[['pca']]@cell.embeddings[,1:30],
                                                         loadings = YS_altas_refer_NRBC_seu[['pca']]@feature.loadings[,1:30]),Other_YS_altas_Ery_seu[['ref.pca']])
YS_altas_Ery_seu[['ref.umap']]=merge(YS_altas_refer_NRBC_seu[['umap1']],Other_YS_altas_Ery_seu[['ref.umap']])
YS_altas_Ery_seu@meta.data[c(rownames(YS_altas_refer_NRBC_seu@meta.data),rownames(Other_YS_altas_Ery_seu@meta.data)),
                           'predicted.cluster']= c(as.character(YS_altas_refer_NRBC_seu$sub.cluster),as.character(Other_YS_altas_Ery_seu$predicted.celltype))

Other_YS_altas_Ery_seu=MapQuery(anchorset =anchors,query = Other_YS_altas_Ery_seu,reference =YS_altas_refer_NRBC_seu,
                                reference.reduction = 'pca',reduction.model = 'umap1',refdata =  list(celltype = "celltype")  )
YS_altas_Ery_seu@meta.data[c(rownames(YS_altas_refer_NRBC_seu@meta.data),rownames(Other_YS_altas_Ery_seu@meta.data)),
                           'predicted.celltype']= c(as.character(YS_altas_refer_NRBC_seu$celltype),as.character(Other_YS_altas_Ery_seu$predicted.celltype))
YS_altas_Ery_seu$predicted.celltype=factor(YS_altas_Ery_seu$predicted.celltype,levels = levels(YS_altas_refer_NRBC_seu$celltype))
DimPlot(YS_altas_Ery_seu,reduction = 'ref.umap',group.by = c('predicted.cluster','predicted.celltype'),cols = cols,label = F);p

YS_altas_Ery_seu$predicted.celltype=factor(YS_altas_Ery_seu$predicted.celltype,levels = levels(YS_altas_refer_NRBC_seu$celltype))
p=DimPlot(YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,group.by = c('id','predicted.celltype'));p
ggsave(p,filename='res_pic/main_figure1/raw_YS_altas_predicted_celltype_NRBC_sample_umap.pdf',width = 12,height = 5)

p=VlnPlot(YS_altas_Ery_seu,features = celltype_marker,group.by = 'predicted.celltype',stack = T)+NoLegend();p
ggsave(p,filename='res_pic/main_figure1/raw_marker_vlnplot.pdf',width = 6,height = 5)

YS_altas_Ery_seu_celltype_Markers=FindAllMarkers(YS_altas_refer_NRBC_seu,group.by = 'predicted.celltype')
saveRDS(YS_altas_Ery_seu_celltype_Markers,file = 'res_data/raw_YS_altas_Ery_seu_celltype_Markers.rds')
top_marker_celltype=YS_altas_Ery_seu_celltype_Markers[YS_altas_Ery_seu_celltype_Markers$pct.1>0.1 & YS_altas_Ery_seu_celltype_Markers$pct.2<0.1 & YS_altas_Ery_seu_celltype_Markers$avg_log2FC>3,] %>% group_by(cluster) %>%top_n(n = 10,wt = pct.1-pct.2)
top_marker_celltype=top_marker_celltype[-grep('^CTC-|^CTA-|^AC104|^RP',top_marker_celltype$gene),]
p=DotPlot(YS_altas_Ery_seu,group.by = 'predicted.celltype',scale = F,features = unique(top_marker_celltype$gene),cols = c('white','firebrick3'))+coord_flip()+RotatedAxis()
p
ggsave(p,filename='res_pic/main_figure1/raw_top_marker_dotplot.pdf',width = 6,height = 10)

rm(Other_YS_altas_Ery_seu,anchors);gc()


# pEry2 :KRT8 和 KRT18 是编码细胞角蛋白（Keratin）的基因，主要参与细胞骨架形成、维持细胞形态和机械稳定性，并调节细胞应激、凋亡及分化
# 'EPCAM','KRT19','KRT7','CDH1', 上皮marker， 等内胚层基因：AFP、TTR、FOXA2 
VlnPlot(YS_altas_Ery_seu,features = c('KRT8','KRT18','EPCAM','KRT19','KRT7','CDH1','VIM','AFP','TTR','FOXA2','GYPA','HBE1', 'HBZ'),stack = T,group.by = 'predicted.celltype')
table(YS_altas_Ery_seu@meta.data[,c('predicted.celltype','id')])

YS_altas_seu=NormalizeData(YS_altas_seu)
FeaturePlot(YS_altas_seu,features = c('KRT8','KRT18'),reduction = 'umap')

DimPlot(YS_altas_seu,group.by = 'LVL1',cols = cols,label = T)
DimPlot(YS_altas_seu,group.by = 'LVL2',cols = cols,label = F)

order_ids=c("CS10_Wang_CS10",  "NI_CS10_YS","CS11_Wang_CS11","NI_CS11_YS" , "CS14_mikola_1", "CS15_F80" ,"NI_CS15_YS",
            "CS17_F138" , "CS18_F61","CS18_F79", "CS22_F32","CS23_F35")
table(unique(YS_altas_Ery_seu$id) %in% order_ids)
YS_altas_Ery_seu$id=factor(YS_altas_Ery_seu$id,levels = order_ids)

table(YS_altas_Ery_seu@meta.data[,c('id','predicted.celltype')])
celltype_df=data.frame(table(YS_altas_Ery_seu@meta.data[,c('id','predicted.celltype')]))
celltype_fre_Pro_df3=data.frame(table(YS_altas_Ery_seu@meta.data[,c('id')]))
colnames(celltype_fre_Pro_df3)=c('id','Count')

p1=ggplot(celltype_df,aes(x = id,y=Freq,fill=predicted.celltype  ))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+scale_fill_manual(values =  col[1:11])+
  theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold')) +ggtitle('primitive NRBCs')+
  geom_label(data = celltype_fre_Pro_df3,aes(x=id,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
p1
ggsave(p1,filename='res_pic/main_figure1/raw_celltype_distribution_barplot.pdf',width = 8,height = 6)

# pEry1:endoderm siganature & stroma signature
VlnPlot(subset(YS_altas_Ery_seu,predicted.celltype=='pEry1'),group.by = 'id',features = c('GYPA','HBB','GATA1','KLF1','ALAS2', top_marker_celltype[top_marker_celltype$cluster %in% c('pEry1','pEry2'),'gene']$gene),stack = T)

saveRDS(YS_altas_Ery_seu,file = 'new_YS_altas_Ery_seu.rds')

def_YS_YS_altas_Ery_seu=subset(YS_altas_Ery_seu,predicted.celltype=='def_like_Ery')
DimPlot(def_YS_YS_altas_Ery_seu,reduction = 'ref.umap',group.by = 'Pre_celltype2',cols = col)
table(def_YS_YS_altas_Ery_seu$id) # 主要来自CS22_F32,其次为CS18_F79 样本 
#CS10_Wang_CS10     NI_CS10_YS CS11_Wang_CS11     NI_CS11_YS  CS14_mikola_1       CS15_F80     NI_CS15_YS      CS17_F138       CS18_F61       CS18_F79       CS22_F32       CS23_F35 
#1              5              0              5              0              5             26              8            149              3            460              0 

def_YS_YS_altas_Ery_seu=FindVariableFeatures(def_YS_YS_altas_Ery_seu) %>% ScaleData()
def_YS_YS_altas_Ery_seu= RunPCA(def_YS_YS_altas_Ery_seu) %>% RunUMAP(dims = 1:10)
DimPlot(def_YS_YS_altas_Ery_seu,reduction = 'umap',group.by = 'Pre_celltype2',cols = col)

def_YS_YS_altas_Ery_seu=singleR_analysis_func(refdata =nrbc_ref_se2,test =as.matrix(def_YS_YS_altas_Ery_seu@assays$RNA$data)[VariableFeatures(def_YS_YS_altas_Ery_seu),],
                                              outdata = def_YS_YS_altas_Ery_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype3')
def_YS_YS_altas_Ery_seu=def_YS_YS_altas_Ery_seu[[1]]
def_YS_YS_altas_Ery_seu$Pre_celltype3=factor(def_YS_YS_altas_Ery_seu$Pre_celltype3,levels =  c('CFUE','ProE', 'eBas', 'lBas' ,'Poly' ,'Orth' ))
table(def_YS_YS_altas_Ery_seu$Pre_celltype3) # Orth :436/662=0.66, 这也挺神奇的，是因为YS niche 原因导致的吗？
# 有批次效应问题，要处理， 
FeaturePlot(def_YS_YS_altas_Ery_seu,features = c('TFRC','CD63','HBE1','HBG1'),reduction = 'umap')
# 高表达 HBF、HBE1
DimPlot(def_YS_YS_altas_Ery_seu,reduction = 'umap',group.by = 'Pre_celltype3',cols = col)
saveRDS(def_YS_YS_altas_Ery_seu,file = 'YS_def_like_Ery_seu.rds')


##############################################################################################################################################
#--------------------------------filt YS NRBC altas------------------#
##############################################################################################################################################
filt_YS_altas_Ery_seu=subset(YS_altas_Ery_seu,predicted.celltype %in% c("EryP-CFC","eBas","mBas","lBas","Poly","Orth" ))
filt_YS_altas_Ery_seu$predicted.celltype=factor(filt_YS_altas_Ery_seu$predicted.celltype,levels = c("EryP-CFC","eBas","mBas","lBas","Poly","Orth" ))
DimPlot(filt_YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,label = T,group.by = c('predicted.celltype'))
# Visualize QC metrics as a violin plot
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
filt_YS_altas_Ery_seu[["percent.mt"]] <- PercentageFeatureSet(filt_YS_altas_Ery_seu, pattern = "^MT-")
VlnPlot(YS_altas_Ery_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = 'predicted.celltype')
filt_YS_altas_Ery_seu <- subset(filt_YS_altas_Ery_seu, subset = nFeature_RNA > 200 & nCount_RNA < 60000 & percent.mt < 10)
rm(new_YS_altas_Ery_seu.rds);gc()

filt_YS_altas_refer_NRBC_seu=subset(filt_YS_altas_Ery_seu,id %in% c('CS18_F61','CS15_F80','CS10_Wang_CS10')) #
ElbowPlot(filt_YS_altas_refer_NRBC_seu)
filt_YS_altas_refer_NRBC_seu=FindVariableFeatures(filt_YS_altas_refer_NRBC_seu) %>%ScaleData() %>%RunPCA() %>%RunHarmony(group.by.vars=c('id'))
filt_YS_altas_refer_NRBC_seu=RunUMAP(filt_YS_altas_refer_NRBC_seu,reduction = 'harmony',dims = 1:15,reduction.name ='umap1',return.model = T,spread =0.5 ,seed.use = 42 )
p=DimPlot(filt_YS_altas_refer_NRBC_seu,reduction = 'umap1',cols=cols,label = T,group.by = c('id','predicted.celltype'));p
ggsave(p,filename='res_pic/main_figure1/filt_YS_altas_refer_celltype_NRBC_sample_umap.pdf',width = 12,height = 5)

Other_YS_altas_Ery_seu=subset(YS_altas_Ery_seu,!id %in% unique(filt_YS_altas_refer_NRBC_seu$id))
anchors=FindTransferAnchors(reference =filt_YS_altas_refer_NRBC_seu,query = Other_YS_altas_Ery_seu, reference.reduction='pca')
Other_YS_altas_Ery_seu=MapQuery(anchorset =anchors,query = Other_YS_altas_Ery_seu,reference =filt_YS_altas_refer_NRBC_seu,reference.reduction = 'pca',reduction.model = 'umap1' )
DimPlot(Other_YS_altas_Ery_seu,reduction = 'ref.umap',cols=cols,group.by = c('id','predicted.celltype'))

ref_umap = filt_YS_altas_refer_NRBC_seu[["umap1"]]
query_umap= Other_YS_altas_Ery_seu[["ref.umap"]]
Key(ref_umap) = "refUMAP_"
Key(query_umap) = "refUMAP_"
combined_umap = merge(ref_umap, query_umap)
filt_YS_altas_Ery_seu[["ref.umap2"]]= combined_umap

filt_YS_altas_Ery_seu$celltype=filt_YS_altas_Ery_seu$predicted.celltype

p1=DimPlot(filt_YS_altas_Ery_seu,reduction = 'ref.umap2',group.by = c('id','celltype'),cols = col)

celltype_marker=c("TFRC","GYPA",'HBE1','CD63','MKI67')
p2=VlnPlot(filt_YS_altas_Ery_seu,features = celltype_marker,group.by = 'celltype',stack = T,cols = cols)+NoLegend()
p2
p=p1+p2;p
ggsave(p,filename='res_pic/main_figure1/filt_umap_cellype_classical_marker.pdf',width = 18,height = 5)

saveRDS(filt_YS_altas_refer_NRBC_seu,file = 'filt_YS_altas_refer_NRBC_seu.rds')
rm(Other_YS_altas_Ery_seu,anchors,filt_YS_altas_refer_NRBC_seu);gc()


filt_YS_altas_Ery_seu$celltype=as.character(filt_YS_altas_Ery_seu$predicted.celltype)
filt_YS_altas_Ery_seu$celltype=factor(filt_YS_altas_Ery_seu$celltype,levels =levels(filt_YS_altas_Ery_seu$predicted.celltype) )
order_ids=c("CS10_Wang_CS10",  "NI_CS10_YS","CS11_Wang_CS11","NI_CS11_YS" , "CS14_mikola_1", "CS15_F80" ,"NI_CS15_YS",
            "CS17_F138" , "CS18_F61","CS18_F79", "CS22_F32","CS23_F35")
table(unique(filt_YS_altas_Ery_seu$id) %in% order_ids)
filt_YS_altas_Ery_seu$id=factor(filt_YS_altas_Ery_seu$id,levels = order_ids)

table(filt_YS_altas_Ery_seu@meta.data[,c('id','celltype')])
celltype_df=data.frame(table(filt_YS_altas_Ery_seu@meta.data[,c('id','celltype')]))
celltype_fre_Pro_df3=data.frame(table(filt_YS_altas_Ery_seu@meta.data[,c('id')]))
colnames(celltype_fre_Pro_df3)=c('id','Count')

p1=ggplot(celltype_df,aes(x = id,y=Freq,fill=celltype  ))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+scale_fill_manual(values =  col[1:11])+
  theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold')) +ggtitle('primitive NRBCs')+
  geom_label(data = celltype_fre_Pro_df3,aes(x=id,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
p1
ggsave(p1,filename='res_pic/main_figure1/filt_celltype_distribution_barplot.pdf',width = 8,height = 6)


filt_YS_altas_Ery_seu_celltype_Markers=FindAllMarkers(filt_YS_altas_Ery_seu,group.by = 'celltype')
saveRDS(filt_YS_altas_Ery_seu_celltype_Markers,file = 'res_data/filt_YS_altas_Ery_seu_celltype_Markers.rds')


top_marker_celltype=filt_YS_altas_Ery_seu_celltype_Markers[filt_YS_altas_Ery_seu_celltype_Markers$pct.1>0.1 & filt_YS_altas_Ery_seu_celltype_Markers$pct.2<0.3 & filt_YS_altas_Ery_seu_celltype_Markers$avg_log2FC>2,] %>% group_by(cluster) %>%top_n(n = 10,wt = pct.1-pct.2)
top_marker_celltype=top_marker_celltype[-grep('^CTC-|^CTA-|^AC104|^AC01|^RP11|CMB9',top_marker_celltype$gene),]
p=DotPlot(filt_YS_altas_Ery_seu,group.by = 'celltype',scale = F,features = unique(top_marker_celltype$gene),cols = c('white','firebrick3'))+coord_flip()+RotatedAxis()
p
# m/l Bas 缺少marker gene

# main celltype: e/m/lBas:Bas
filt_YS_altas_Ery_seu$main_celltype=as.character(filt_YS_altas_Ery_seu$celltype)
filt_YS_altas_Ery_seu@meta.data[filt_YS_altas_Ery_seu$main_celltype %in% c('eBas','mBas','lBas'),'main_celltype']='Bas'
filt_YS_altas_Ery_seu$main_celltype=factor(filt_YS_altas_Ery_seu$main_celltype,levels = c("EryP-CFC",'Bas','Poly','Orth'))
filt_YS_altas_Ery_seu_mcelltype_Markers=FindAllMarkers(filt_YS_altas_Ery_seu,group.by = 'main_celltype')
saveRDS(filt_YS_altas_Ery_seu_mcelltype_Markers,file = 'res_data/filt_YS_altas_Ery_seu_mcelltype_Markers.rds')

top_mmarker_celltype=filt_YS_altas_Ery_seu_mcelltype_Markers[filt_YS_altas_Ery_seu_mcelltype_Markers$pct.1>0.1 & filt_YS_altas_Ery_seu_mcelltype_Markers$pct.2<0.3 & filt_YS_altas_Ery_seu_mcelltype_Markers$avg_log2FC>2,] %>% group_by(cluster) %>%top_n(n = 10,wt = pct.1-pct.2)
top_mmarker_celltype=top_mmarker_celltype[-grep('^CTC-|^CTA-|^AC104|^AC01|^RP11',top_mmarker_celltype$gene),]
p=DotPlot(filt_YS_altas_Ery_seu,group.by = 'main_celltype',scale = F,features = unique(top_mmarker_celltype$gene),cols = c('white','firebrick3'))+coord_flip()+RotatedAxis()
p

# 两个对比：early Bas vs mid Bas, mid Bas vs late Bas;
Idents(filt_YS_altas_Ery_seu)='celltype'
eBas_mBas_celltype_Markers=FindMarkers(filt_YS_altas_Ery_seu,ident.1 = 'mBas',ident.2 = 'eBas')
mBas_lBas_celltype_Markers=FindMarkers(filt_YS_altas_Ery_seu,ident.1 = 'lBas',ident.2 = 'mBas')

top1_Bas_marker=eBas_mBas_celltype_Markers[eBas_mBas_celltype_Markers$avg_log2FC>=1 & eBas_mBas_celltype_Markers$pct.1>0.2  , ]
top2_Bas_marker=mBas_lBas_celltype_Markers[mBas_lBas_celltype_Markers$avg_log2FC>=1 & mBas_lBas_celltype_Markers$pct.1>0.2  , ]
dim(top2_Bas_marker)
cho_bas_genes= intersect(rownames(top1_Bas_marker),rownames(top2_Bas_marker))
DotPlot(filt_YS_altas_Ery_seu,group.by = 'celltype',scale = F,features =cho_bas_genes,cols = c('white','firebrick3'))+coord_flip()+RotatedAxis()

top1_Bas_marker=eBas_mBas_celltype_Markers[eBas_mBas_celltype_Markers$avg_log2FC>=0.9 & eBas_mBas_celltype_Markers$pct.1>0.2  , ]
top2_Bas_marker=mBas_lBas_celltype_Markers[mBas_lBas_celltype_Markers$avg_log2FC>=0.9 & mBas_lBas_celltype_Markers$pct.1>0.2  , ]
dim(top2_Bas_marker)
cho_bas_genes= intersect(rownames(top1_Bas_marker),rownames(top2_Bas_marker))
DotPlot(filt_YS_altas_Ery_seu,group.by = 'celltype',scale = F,features =cho_bas_genes,cols = c('white','firebrick3'))+coord_flip()+RotatedAxis()
# chose the 5 marker increasing from early to late Bas;
cho_bas_genes=c("MT1M","ABTB1","YPEL4","CYB5R1","OPTN","TUBA1A"  )

p=DotPlot(filt_YS_altas_Ery_seu,group.by = 'celltype',scale = F,features = unique(c(unique(top_marker_celltype$gene)[1:22],cho_bas_genes,unique(top_marker_celltype$gene)[23:30])),cols = c('white','firebrick3'))+coord_flip()+RotatedAxis()
p
ggsave(p,filename='res_pic/main_figure1/filt_celltype_top_marker.pdf',width = 6,height = 10)


YS_altas_Ery_seu_celltype_Markers_list=list()
YS_altas_Ery_seu_celltype_gseGO_list=list()
YS_altas_Ery_seu_celltype_gseGO_df=data.frame()

for( cluster in levels(filt_YS_altas_Ery_seu_celltype_Markers$cluster)){
  YS_altas_Ery_seu_celltype_Markers_list[[cluster]]=filt_YS_altas_Ery_seu_celltype_Markers[filt_YS_altas_Ery_seu_celltype_Markers$cluster==cluster,'avg_log2FC']
  names( YS_altas_Ery_seu_celltype_Markers_list[[cluster]])=filt_YS_altas_Ery_seu_celltype_Markers[filt_YS_altas_Ery_seu_celltype_Markers$cluster==cluster,'gene']
  YS_altas_Ery_seu_celltype_Markers_list[[cluster]]=sort(YS_altas_Ery_seu_celltype_Markers_list[[cluster]],decreasing = T)
  YS_altas_Ery_seu_celltype_gseGO_list[[cluster]]= gseGO(geneList = YS_altas_Ery_seu_celltype_Markers_list[[cluster]],ont = 'BP',OrgDb = org.Hs.eg.db,keyType = 'SYMBOL')
  tmp=YS_altas_Ery_seu_celltype_gseGO_list[[cluster]]@result
  tmp$cluster=cluster
  YS_altas_Ery_seu_celltype_gseGO_df=rbind(YS_altas_Ery_seu_celltype_gseGO_df,tmp)
  
}

saveRDS(YS_altas_Ery_seu_celltype_gseGO_df,file = 'res_data/filt_YS_altas_Ery_seu_celltype_gseGO_df.rds')

top_YS_altas_Ery_seu_celltype_gseGO_df=YS_altas_Ery_seu_celltype_gseGO_df[YS_altas_Ery_seu_celltype_gseGO_df$NES>0 & YS_altas_Ery_seu_celltype_gseGO_df$p.adjust <0.01,] %>%group_by(cluster) %>%  slice_head(n = 5)
#top_YS_altas_Ery_seu_celltype_gseGO_df2=YS_altas_Ery_seu_celltype_gseGO_df[YS_altas_Ery_seu_celltype_gseGO_df$NES<0 & YS_altas_Ery_seu_celltype_gseGO_df$p.adjust <0.01,] %>%group_by(cluster) %>%  slice_head(n = 5)
#top_YS_altas_Ery_seu_celltype_gseGO_df=rbind(top_YS_altas_Ery_seu_celltype_gseGO_df1,top_YS_altas_Ery_seu_celltype_gseGO_df2)
#top_YS_altas_Ery_seu_celltype_gseGO_df=YS_altas_Ery_seu_celltype_gseGO_df[YS_altas_Ery_seu_celltype_gseGO_df$Description %in% unique(top_YS_altas_Ery_seu_celltype_gseGO_df$Description),]
top_YS_altas_Ery_seu_celltype_gseGO_df$cluster=factor(top_YS_altas_Ery_seu_celltype_gseGO_df$cluster,levels = levels(filt_YS_altas_Ery_seu$celltype))
top_YS_altas_Ery_seu_celltype_gseGO_df=top_YS_altas_Ery_seu_celltype_gseGO_df[order(top_YS_altas_Ery_seu_celltype_gseGO_df$cluster),]
top_YS_altas_Ery_seu_celltype_gseGO_df$Description=factor(top_YS_altas_Ery_seu_celltype_gseGO_df$Description,levels =unique( top_YS_altas_Ery_seu_celltype_gseGO_df$Description))
p=ggplot(top_YS_altas_Ery_seu_celltype_gseGO_df,aes(x=cluster,y=Description,color=NES,size=-log10(p.adjust)))+geom_point()+theme_classic()+
  scale_color_gradient(low ='gray',high = 'firebrick3')+RotatedAxis()
p
ggsave(p,filename='res_pic/main_figure1/filt_enrichedGO_dotplot.pdf',width = 8,height = 6)

# check the gene enriched in immune related pathways
Idents(filt_YS_altas_Ery_seu)='celltype'
VlnPlot(filt_YS_altas_Ery_seu,features = unique(unlist(strsplit(top_YS_altas_Ery_seu_celltype_gseGO_df$core_enrichment[1:2],split = '/',fixed = T))),stack = T)+NoLegend()
VlnPlot(filt_YS_altas_Ery_seu,features = c('GAL','BID','MIF','CD81','IFGBP2','CD320'),stack = T)

saveRDS(filt_YS_altas_Ery_seu,file = 'filt_YS_altas_Ery_seu.rds')
rm(list=ls());gc()


#####################################################################################################################################################################
###########-------------------------------------------HBB--------------------------------#
#####################################################################################################################################################################
p1=FeaturePlot(subset(filt_YS_altas_Ery_seu,HBB >0),features = c('HBB'),reduction = 'ref.umap2',cols = c('gray','firebrick3'));p1
ggsave(p1,filename='res_pic/main_figure1/filt_HBB_expression_umap.pdf',width = 6,height = 5)

p2=DimPlot(subset(filt_YS_altas_Ery_seu,HBB >0),group.by = 'celltype',cols = cols,reduction = 'ref.umap2');p2
ggsave(p2,filename='res_pic/main_figure1/filt_HBB_pos_celltype_umap.pdf',width = 6,height =5)



hbb_expr <- GetAssayData(filt_YS_altas_Ery_seu, assay = "RNA", layer = "counts")["HBB", ]
filt_YS_altas_Ery_seu@meta.data$HBB_status <- ifelse(hbb_expr >=1, "positive", "negtive")
tmp_df=table(filt_YS_altas_Ery_seu@meta.data[,c('id','HBB_status')])
tmp_df=data.frame(tmp_df)

celltype_fre_Pro_df=data.frame(table(filt_YS_altas_Ery_seu$id));colnames(celltype_fre_Pro_df)=c('id','Count')

p=ggplot(data =tmp_df,aes(x=id,y = Freq,fill=HBB_status  ))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+RotatedAxis()+
  geom_label(data = celltype_fre_Pro_df,aes(x=id,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
p
ggsave(p,filename='res_pic/main_figure1/filt_HBB_ratio_distribution_barplot.pdf',width = 8,height = 6)


#############################################################################################################################
#---------add  the subcelltype of nRBC to YS ALTAS----------#
#############################################################################################################################
YS_altas_seu=readRDS('../NRBC_YS_altas/raw_ref_data/YS_altas_seu.rds')
YS_altas_seu$subcelltype=YS_altas_seu$LVL2
sort(table(YS_altas_seu$subcelltype))# MONO_Mac_DC 细胞很少
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('EOSINOPHIL_BASOPHIL','MAST_CELL')]='EO/BASO/MAST'
#DimPlot(subset(YS_altas_seu,subcelltype %in% c('MONOCYTE_MACROPHAGE','MONO MAC PRE DC2','PDC PRECURSOR','MONO MAC DC2','PROMONOCYTE','MONOCYTE_YS_1', 'MONOCYTE','MOP','DC')),group.by = 'LVL2',cols = cols)
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('MONO MAC PRE DC2','PDC PRECURSOR','MONO MAC DC2','DC')]='MOMO_MAC_DC'
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('PROMONOCYTE')]='MOP'
YS_altas_seu$subcelltype[YS_altas_seu$subcelltype %in% c('MONOCYTE_YS_1')]='MONOCYTE'
YS_altas_seu$subcelltype[YS_altas_seu$LVL3 %in% c('DEF_HSPC_1','DEF_HSPC_2')]='DEF_HSPC'
YS_altas_seu$subcelltype[YS_altas_seu$LVL3 %in% c('HSPC_1','HSPC_2')]='HSPC'
YS_altas_seu$subcelltype[YS_altas_seu$LVL3 %in% c('PRIM_HSPC_1','PRIM_HSPC_2')]='PRIM_HSPC'
YS_altas_seu=subset(YS_altas_seu,id!='CS14_F37')

cho_cells=intersect(rownames(filt_YS_altas_Ery_seu@meta.data),rownames(YS_altas_seu@meta.data))
YS_altas_seu@meta.data[cho_cells,'subcelltype']=as.character(filt_YS_altas_Ery_seu@meta.data[cho_cells,'main_celltype'])
sort(table(YS_altas_seu$subcelltype)) # ERYTHROID: 6514  , definitive NRBC，
YS_altas_seu$subcelltype=gsub(pattern = 'MK',replacement = 'MEGAKARYOCYTE',YS_altas_seu$subcelltype) 

celltype_levels=c('EO/BASO/MAST','MEGAKARYOCYTE','EryP-CFC','Bas','Poly','Orth', 'ERYTHROID','MACROPHAGE','PRIM_HSPC','DEF_HSPC','HSPC','CMP','LMPP','MEMP','MOP','MONOCYTE','MOMO_MAC_DC','ELP','ILC','NK','B_CELL','NEUTROPHIL_PRECURSOR','MONOCYTE_MACROPHAGE',                 
                  'MESOTHELIUM','SMOOTH_MUSCLE', 'ENDODERM','ENDOTHELIUM','FIBROBLAST' )
table(unique(YS_altas_seu$subcelltype )%in% celltype_levels)
YS_altas_seu$subcelltype=factor(YS_altas_seu$subcelltype,levels = celltype_levels)
DimPlot(subset(YS_altas_seu,subcelltype!='ERYTHROID'),group.by = 'subcelltype',cols = cols,raster = F) # 169,494-599=168,895 cells
saveRDS(YS_altas_seu,file ='../NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' )


######################################################################################################################################################################
# ----------------------validated in other YS /embryo scRNAseq data--------------#
######################################################################################################################################################################


# CD12-CS15 时期血液系统单细胞数据,此时，只有 YS 产生红细胞，FL还未产生 dNRBC

# -----------------------process-------------------------------------------------------------------#
if(F){
  library(Matrix)
  early_organogenesis_count=readMM(file = gzfile('../ref_data/ref_scRNAseq_data/GSE157329_early_organogenesis_in_human_embryos/GSE157329_raw_counts.mtx.gz'),header = T)
  dim(early_organogenesis_count)
  
  gene_info=read.table(gzfile('../ref_data/ref_scRNAseq_data/GSE157329_early_organogenesis_in_human_embryos/GSE157329_gene_annotate.txt.gz'),header = T)
  head(gene_info)
  rownames(early_organogenesis_count)=gene_info$gene_short_name
  dim(gene_info)
  
  
  early_organogenesis_meta=read.table(file = gzfile('../ref_data/ref_scRNAseq_data/GSE157329_early_organogenesis_in_human_embryos/GSE157329_cell_annotate.txt.gz'),sep="\t",header = T)
  head(early_organogenesis_meta)
  rownames(early_organogenesis_meta)=early_organogenesis_meta$cell_id;early_organogenesis_meta=early_organogenesis_meta[,-1]
  colnames(early_organogenesis_count)=rownames(early_organogenesis_meta)
  early_organogenesis_count=early_organogenesis_count[rowSums(early_organogenesis_count)>5,]
  dup_num=which(rownames(early_organogenesis_count) %in% rownames(early_organogenesis_count)[duplicated(rownames(early_organogenesis_count))])
  dup_early_organogenesis_count=early_organogenesis_count[dup_num,]
  
  dup_early_organogenesis_count=data.frame(LINC01238=colSums(dup_early_organogenesis_count[1:2,]),
                                           CYB561D2=colSums(dup_early_organogenesis_count[3:4,]),
                                           MATR3=colSums(dup_early_organogenesis_count[5:6,]),
                                           TMSB15B=colSums(dup_early_organogenesis_count[7:8,]),
                                           LINC01505=colSums(dup_early_organogenesis_count[9:10,]),
                                           GOLGA8M=colSums(dup_early_organogenesis_count[9:10,])
  )
  dup_early_organogenesis_count=t(dup_early_organogenesis_count)
  
  
  early_organogenesis_count=early_organogenesis_count[-dup_num,]
  early_organogenesis_count=rbind(early_organogenesis_count,dup_early_organogenesis_count)
  early_organogenesis_count=as(early_organogenesis_count,'dgTMatrix')
  
  early_organogenesis_seu=CreateSeuratObject(counts =early_organogenesis_count, meta.data =early_organogenesis_meta,project = 'early_organogenesis' )
  rm(early_organogenesis_count,early_organogenesis_meta)
  
  
  
  early_organogenesis_seu=NormalizeData(early_organogenesis_seu)
  early_organogenesis_seu[['RNA']]=split(early_organogenesis_seu[['RNA']],f = early_organogenesis_seu$sample)
  early_organogenesis_seu=subset(early_organogenesis_seu,developmental.system=='blood')
  
  early_organogenesis_seu=FindVariableFeatures(early_organogenesis_seu)
  early_organogenesis_seu=ScaleData(early_organogenesis_seu)
  early_organogenesis_seu=RunPCA(early_organogenesis_seu)
  ElbowPlot(early_organogenesis_seu)
  early_organogenesis_seu=RunUMAP(early_organogenesis_seu,dims = 1:15)
  DimPlot(early_organogenesis_seu,group.by =c('embryo','dissection_part','stage','annotation'),cols = col)
  
  FeaturePlot(early_organogenesis_seu,features = c('HBE1','HBZ','HBG1','HBG2'),cols = c('white','firebrick3'))
  
  early_organogenesis_seu=IntegrateLayers(early_organogenesis_seu,method = HarmonyIntegration,orig.reduction = "pca",new.reduction='harmony')
  early_organogenesis_seu=RunUMAP(early_organogenesis_seu,reduction = 'harmony',reduction.name = 'har_umap',dims = 1:15)
  #viscera: 内脏, 内脏下干:viscera-lowerTrunk
  DimPlot(early_organogenesis_seu,group.by =c('embryo','dissection_part','stage','annotation'),cols = col,reduction = 'har_umap')
  early_organogenesis_seu[['RNA']]=JoinLayers(early_organogenesis_seu[['RNA']])
  
  FeaturePlot(early_organogenesis_seu,features = c('HBE1','HBZ','HBG1','HBG2'),cols = c('gray','firebrick3'),reduction = 'har_umap')
  FeaturePlot(early_organogenesis_seu,features = c('HBE1','HBZ','HBG1','HBG2'),cols = c('gray','firebrick3'),reduction = 'har_umap',split.by = 'stage')
  
  early_organogenesis_seu=RunUMAP(early_organogenesis_seu,reduction = 'harmony',reduction.name = 'har_umap1',dims = 1:6)
  DimPlot(early_organogenesis_seu,group.by =c('embryo','dissection_part','stage','annotation'),cols = col,reduction = 'har_umap1')
  FeaturePlot(early_organogenesis_seu,features = c('KIT','TFRC','GYPA','CD63'),cols = c('gray','firebrick3'),reduction = 'har_umap1')
  
  
  VlnPlot(early_organogenesis_seu,features = 'HBD',split.by = 'annotation',cols = col)+
    FeaturePlot(early_organogenesis_seu,features = 'HBD',cols = c('gray','firebrick3'),reduction = 'har_umap')
  
  DotPlot(subset(early_organogenesis_seu,annotation=='erythroid'),group.by = 'stage',features = heamoglobin_genes)+RotatedAxis()
  saveRDS(early_organogenesis_seu,file = '../ref_data/early_organogenesis_seu.rds')
  
  meta_data=read.csv(gzfile('ref_data/ref_scRNAseq_data/GSE157329_early_organogenesis_in_human_embryos/GSE156793_S1_metadata_cells.txt.gz'))
  table(meta_data$Organ ) # 包含了多个组织的单细胞数据,非空间组数据
  #  Adrenal Cerebellum   Cerebrum        Eye      Heart  Intestine     Kidney      Liver       Lung     Muscle   Pancreas   Placenta     Spleen    Stomach     Thymus 
  ggplot(meta_data[meta_data$Organ=='Placenta',],aes(x=Main_cluster_umap_1 ,y=Main_cluster_umap_2))+geom_point()+theme_bw()
  VlnPlot(early_organogenesis_seu,features = c('CXCL12'))+NoLegend()
  
   
}


early_organogenesis_seu=readRDS('../ref_data/early_organogenesis_seu.rds')
Ery_early_organogenesis_seu=subset(early_organogenesis_seu,annotation=='erythroid')
Ery_early_organogenesis_seu=RunUMAP(Ery_early_organogenesis_seu,reduction = 'harmony',dims = 1:8,reduction.name = 'umap1')
Ery_early_organogenesis_seu=FindVariableFeatures(Ery_early_organogenesis_seu) %>%NormalizeData()
Ery_early_organogenesis_seu=singleR_analysis_func(refdata =nrbc_ref_se,test =as.matrix(Ery_early_organogenesis_seu@assays$RNA$data)[VariableFeatures(Ery_early_organogenesis_seu),],outdata = Ery_early_organogenesis_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype')
Ery_early_organogenesis_seu=Ery_early_organogenesis_seu[[1]]
Ery_early_organogenesis_seu$Pre_celltype=factor(Ery_early_organogenesis_seu$Pre_celltype,levels =c('BFUE','CFUE','ProE','eBas','mBas','lBas','Poly','Orth'))
FeaturePlot(Ery_early_organogenesis_seu,features = c('HBE1'),reduction = 'umap1')
Ery_early_organogenesis_seu <- CellCycleScoring(Ery_early_organogenesis_seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

Ery_early_organogenesis_seu=FindNeighbors(Ery_early_organogenesis_seu,reduction = 'harmony',dims = 1:20)
Ery_early_organogenesis_seu=FindClusters(Ery_early_organogenesis_seu,resolution = 0.2)
p=DimPlot(Ery_early_organogenesis_seu,reduction = 'umap1',group.by = c('seurat_clusters','Pre_celltype'),cols = col)
p

filt_YS_altas_refer_NRBC_seu=readRDS('filt_YS_altas_refer_NRBC_seu.rds')
anchors=FindTransferAnchors(reference =filt_YS_altas_refer_NRBC_seu,query = Ery_early_organogenesis_seu, reference.reduction='pca')
Ery_early_organogenesis_seu=MapQuery(anchorset =anchors,query = Ery_early_organogenesis_seu,reference =filt_YS_altas_refer_NRBC_seu,reference.reduction = 'pca',reduction.model = 'umap1',refdata = 'predicted.celltype' )
Ery_early_organogenesis_seu$predicted.id=factor(Ery_early_organogenesis_seu$predicted.id,levels = levels(filt_YS_altas_refer_NRBC_seu$predicted.celltype))

DimPlot(Ery_early_organogenesis_seu,reduction = 'umap1',cols=cols,group.by = c('predicted.id'))+
  DimPlot(Ery_early_organogenesis_seu,reduction = 'ref.umap',cols=cols,group.by = c('predicted.id','seurat_clusters'))

p=DimPlot(Ery_early_organogenesis_seu,reduction = 'ref.umap',group.by = c('stage','predicted.id'),cols = col);p
ggsave(p,filename = 'res_pic/main_figure1/YS_embryo_blood_systerm_umap_stage_celltype.pdf',width = 5,height = 12,dpi = 300)

Ery_early_organogenesis_seu$celltype=Ery_early_organogenesis_seu$predicted.id # reference 更准确

celltype_df=data.frame(table(Ery_early_organogenesis_seu@meta.data[,c('stage','celltype')]))
celltype_fre_Pro_df3=data.frame(table(Ery_early_organogenesis_seu@meta.data[,c('stage')]))
colnames(celltype_fre_Pro_df3)=c('stage','Count')

p=ggplot(celltype_df,aes(x = stage,y=Freq,fill=celltype))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+scale_fill_manual(values =  col[1:11])+
  theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold')) +
  geom_label(data = celltype_fre_Pro_df3,aes(x=stage,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  ggtitle('embryo blood system')
p
ggsave(p,filename = 'res_pic/main_figure1/YS_embryo_blood_systerm_nRBC_analysis_barplot.pdf',width = 4,height = 8,dpi = 300)

FeaturePlot(Ery_early_organogenesis_seu,reduction = 'ref.umap',features = c('HBE1','HBB'))

saveRDS(Ery_early_organogenesis_seu,file = '../NRBC_YS_altas/raw_ref_data/Ery_early_organogenesis_seu.rds')

p1=DimPlot(Ery_early_organogenesis_seu,reduction = 'ref.umap',group.by = c('predicted.id'),cols = col)+ggtitle('Blood')
p1

############################################## YS-derived NRBC in FL altas #############################
FL_altas_Ery_seu=readRDS('../NRBC_FL_altas/tmp_FL_altas_Ery_seu.rds')
table(FL_altas_Ery_seu$subcelltype)
primitive_NRBC_inFL_seu=subset(FL_altas_Ery_seu,subcelltype %in% c('YS_Bas/Poly','YS_Orth'))
primitive_NRBC_inFL_seu=subset(primitive_NRBC_inFL_seu,age %in% c('CS14_4PCW' ,'CS15_5PCW' ,'CS17_6PCW','CS18' ) )
DimPlot(primitive_NRBC_inFL_seu,reduction = 'umap2',group.by = 'subcelltype')


filt_YS_altas_refer_NRBC_seu=readRDS('filt_YS_altas_refer_NRBC_seu.rds')
anchors=FindTransferAnchors(reference =filt_YS_altas_refer_NRBC_seu,query = primitive_NRBC_inFL_seu, reference.reduction='pca')
primitive_NRBC_inFL_seu=MapQuery(anchorset =anchors,query = primitive_NRBC_inFL_seu,reference =filt_YS_altas_refer_NRBC_seu,reference.reduction = 'pca',reduction.model = 'umap1',refdata = 'predicted.celltype' )
primitive_NRBC_inFL_seu$predicted.id=factor(primitive_NRBC_inFL_seu$predicted.id,levels = levels(filt_YS_altas_refer_NRBC_seu$predicted.celltype))
DimPlot(primitive_NRBC_inFL_seu,reduction = 'ref.umap',group.by = 'predicted.id',cols = cols)+ggtitle('FL ')

#----------------YS NRBC in EBM ----------------------#
fetal_bm_seu=readRDS('../NRBC_BM_altas/tmp_embryo_bm_seu.rds')
YS_NRBC_inEBM_seu=subset(fetal_bm_seu,celltype  %in% c('YS_Bas'))

anchors=FindTransferAnchors(reference =filt_YS_altas_refer_NRBC_seu,query = YS_NRBC_inEBM_seu, reference.reduction='pca')
YS_NRBC_inEBM_seu=MapQuery(anchorset =anchors,query = fetal_bm_seu,reference =filt_YS_altas_refer_NRBC_seu,reference.reduction = 'pca',reduction.model = 'umap1',refdata = 'predicted.celltype' )
YS_NRBC_inEBM_seu$predicted.id=factor(YS_NRBC_inEBM_seu$predicted.id,levels = levels(filt_YS_altas_refer_NRBC_seu$predicted.celltype))

DimPlot(YS_NRBC_inEBM_seu,reduction = 'ref.umap',group.by = 'predicted.id',cols = cols)+ggtitle('EBM')

Ery_early_organogenesis_seu$tissue='Blood:CS12-CS16'
primitive_NRBC_inFL_seu$tissue='FL:CS14-CS18'
YS_NRBC_inEBM_seu$tissue='EBM:CS13-CS22'

primitive_NRBC_inFL_seu$stage='CS14-CS17'
primitive_NRBC_inFL_seu$stage[primitive_NRBC_inFL_seu$age=='CS18']='CS18'
YS_NRBC_inEBM_seu$stage='CS13-CS15'
YS_NRBC_inEBM_seu$stage[YS_NRBC_inEBM_seu$age %in% c('CS20' , 'CS22' ,'CS22b') ]='CS20-CS22'

other_tissue_YS_NRBC_seu=merge(Ery_early_organogenesis_seu,c(primitive_NRBC_inFL_seu,YS_NRBC_inEBM_seu))
other_tissue_YS_NRBC_seu[['ref.umap']]=merge(Ery_early_organogenesis_seu[['ref.umap']],c(primitive_NRBC_inFL_seu[['ref.umap']],YS_NRBC_inEBM_seu[['ref.umap']]))
other_tissue_YS_NRBC_seu$predicted.id=factor(other_tissue_YS_NRBC_seu$predicted.id,levels = levels(filt_YS_altas_refer_NRBC_seu$predicted.celltype))
rm(Ery_early_organogenesis_seu,primitive_NRBC_inFL_seu,YS_NRBC_inEBM_seu);gc()

p=DimPlot(other_tissue_YS_NRBC_seu,reduction = 'ref.umap',group.by = 'tissue',cols = cols)+
  DimPlot(other_tissue_YS_NRBC_seu,reduction = 'ref.umap',group.by = 'predicted.id',cols = cols)


other_tissue_YS_NRBC_seu$stage=paste(strsplit2(other_tissue_YS_NRBC_seu$tissue,split = ':')[,1],other_tissue_YS_NRBC_seu$stage)


celltype_df=data.frame(table(other_tissue_YS_NRBC_seu@meta.data[,c('stage','predicted.id')]))
celltype_fre_Pro_df3=data.frame(table(other_tissue_YS_NRBC_seu@meta.data[,c('stage')]))
colnames(celltype_fre_Pro_df3)=c('stage','Count')

p2=ggplot(celltype_df,aes(x = stage,y=Freq,fill=predicted.id))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+scale_fill_manual(values =  col[1:11])+
  theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold')) +
  geom_label(data = celltype_fre_Pro_df3,aes(x=stage,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  ggtitle('YS-derived NRBC in other tissues')
p2
p=p+p2
ggsave(p,filename = 'res_pic/main_figure1/other_tissue_YSnRBC_analysis_barplot.pdf',width = 18,height = 6,dpi = 300)

FeaturePlot(other_tissue_YS_NRBC_seu,reduction = 'ref.umap',features = c('HBE1','HBB'),cols =  c('gray','firebrick3'))

saveRDS(other_tissue_YS_NRBC_seu,file = 'other_tissue_YS_NRBC_seu.rds')

#--------------------------------------  HBE1  & HBB-----------------------------------------#
FL_altas_Ery_seu=subset(FL_altas_Ery_seu, subcelltype %in% c('BFUE/CFUE','ProE','Bas','Poly','Orth'))
p=FeaturePlot(subset(FL_altas_Ery_seu, HBE1 >0),cols = c('gray','firebrick3'),features = c('HBE1'),reduction = 'umap2')+
   DimPlot(subset(FL_altas_Ery_seu, HBE1 >0),group.by = 'subcelltype',cols = cols,reduction = 'umap2')+
  DimPlot(subset(FL_altas_Ery_seu, HBE1 >0),group.by = 'age',cols = cols,reduction = 'umap2')

p

ggsave(p,filename = 'res_pic/main_figure1/HBE1_FL_dNRBC.pdf',width = 6,height = 15,dpi = 300)


hbe1_expr <- GetAssayData(FL_altas_Ery_seu, assay = "RNA", layer = "counts")["HBE1", ]
FL_altas_Ery_seu@meta.data$HBE1_status <- ifelse(hbe1_expr >=1, "positive", "negtive")
tmp_df=table(FL_altas_Ery_seu@meta.data[,c('id','HBE1_status')])
tmp_df=data.frame(tmp_df)
rm(hbe1_expr);gc()                  

celltype_fre_Pro_df=data.frame(table(FL_altas_Ery_seu$id));colnames(celltype_fre_Pro_df)=c('id','Count')

p=ggplot(data =tmp_df,aes(x=id,y = Freq,fill=HBE1_status  ))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+RotatedAxis()+
  geom_label(data = celltype_fre_Pro_df,aes(x=id,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
p
ggsave(p,filename='res_pic/main_figure1/FL_dNRBC_HBE1_ratio_distribution_barplot.pdf',width = 15,height = 6)


#--------------------------- FBM  NRBC--------------------------#
BM_NRBC_altas_seu=readRDS('../NRBC_BM_altas/res_data/BM_NRBC_altas_seu.rds')
BM_NRBC_altas_seu=subset(BM_NRBC_altas_seu,stage=='FBM')
table(BM_NRBC_altas_seu$subcelltype)

BM_NRBC_altas_seu=subset(BM_NRBC_altas_seu, subcelltype %in% levels(BM_NRBC_altas_seu$subcelltype)[1:5] )
p=FeaturePlot(subset(BM_NRBC_altas_seu, HBE1 >0),cols = c('gray','firebrick3'),features = c('HBE1'),reduction = 'ref.umap')+
  DimPlot(subset(BM_NRBC_altas_seu, HBE1 >0),group.by = c('subcelltype'),cols = cols,reduction = 'ref.umap')+
  DimPlot(subset(BM_NRBC_altas_seu, HBE1 >0),group.by = c('age'),cols = cols,reduction = 'ref.umap')
p
ggsave(p,filename = 'res_pic/main_figure1/HBE1_FBM_dNRBC.pdf',width = 6,height = 15,dpi = 300)


hbe1_expr <- GetAssayData(BM_NRBC_altas_seu, assay = "RNA", layer = "counts")["HBE1", ]
BM_NRBC_altas_seu@meta.data$HBE1_status <- ifelse(hbe1_expr >=1, "positive", "negtive")
tmp_df=table(BM_NRBC_altas_seu@meta.data[,c('id','HBE1_status')])
tmp_df=data.frame(tmp_df)
tmp_df=tmp_df[tmp_df$Freq >0,]
rm(hbe1_expr);gc()                  

celltype_fre_Pro_df=data.frame(table(BM_NRBC_altas_seu$id));colnames(celltype_fre_Pro_df)=c('id','Count')
celltype_fre_Pro_df=celltype_fre_Pro_df[celltype_fre_Pro_df$id %in% tmp_df$id,]
p=ggplot(data =tmp_df,aes(x=id,y = Freq,fill=HBE1_status  ))+geom_bar(position = "fill", stat = "identity",alpha=0.9)+theme_classic()+RotatedAxis()+
  geom_label(data = celltype_fre_Pro_df,aes(x=id,y=1.16,label=Count,fill=NULL),vjust ='top',show.legend = F, label.size = NA)+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
p
ggsave(p,filename='res_pic/main_figure1/FBM_dNRBC_HBE1_ratio_distribution_barplot.pdf',width = 8,height = 6)






