#---------------------------------------------------tissue RNA marker ----------------------------------------------------#

################################################################################################################################################################
#---------------------------------------------------prepare pkg and function ----------------------------------------------------#
################################################################################################################################################################
library(ggrepel)
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(SingleCellExperiment)
library(RColorBrewer )
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tibble)
library(pheatmap)
library(clusterProfiler)
library(ggplotify)
library(ggplot2)
library(reshape2)
library(corrplot)
library(APL)
library(patchwork)
library(yulab.utils)

#install.packages("patchwork")
library(RColorBrewer )
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
getwd()
setwd('/home/gibh/2021_NRBC_chlyu')

find_mDEGs_func=function(seu,group,sfile){
  early_Ery_tissue_markers=FindAllMarkers(subset(seu,Ery_stage=='early_Ery'))
  mid_Ery_tissue_markers=  FindAllMarkers(subset(seu,Ery_stage=='mid_Ery'))
  late_Ery_tissue_markers= FindAllMarkers(subset(seu,Ery_stage=='late_Ery'))
  
  early_Ery_tissue_markers$celltype='early_Ery'
  mid_Ery_tissue_markers$celltype='mid_Ery'
  late_Ery_tissue_markers$celltype='late_Ery'
  all_Ery_tissue_markers=rbind(early_Ery_tissue_markers,rbind(mid_Ery_tissue_markers,late_Ery_tissue_markers))
  rownames(all_Ery_tissue_markers)=NULL
  all_Ery_tissue_markers$group=group
  
  write.csv(all_Ery_tissue_markers,file = sfile,quote = F)
  rm(early_Ery_tissue_markers,mid_Ery_tissue_markers,late_Ery_tissue_markers);gc()
  
  pos_all_Ery_tissue_markers=all_Ery_tissue_markers[all_Ery_tissue_markers$avg_log2FC >0  & all_Ery_tissue_markers$p_val_adj <0.01,]
  if(dim(pos_all_Ery_tissue_markers)[2] <1 ){ stop('no DEGs are found!!!')}
  
  count_tmp_df=cbind( data.frame(table(pos_all_Ery_tissue_markers[,c('cluster','celltype')])[1,]), data.frame(table(pos_all_Ery_tissue_markers[,c('cluster','celltype')])[2,]))
  colnames(count_tmp_df)= strsplit(group,'_')[[1]]
  type1=strsplit(group,split = '_')[[1]][1]
  type2=strsplit(group,split = '_')[[1]][2]
  
  count_tmp_df['all',type1]=length(unique(pos_all_Ery_tissue_markers$gene[pos_all_Ery_tissue_markers$cluster==type1]))
  count_tmp_df['all',type2]=length(unique(pos_all_Ery_tissue_markers$gene[pos_all_Ery_tissue_markers$cluster==type2]))
  count_tmp_df['add_all',]=length(unique(pos_all_Ery_tissue_markers$gene))
  
  
  return(list(all_Ery_tissue_markers,count_tmp_df))
}

subcelltype_gseGO_func=function(RNA_markers,keyType='SYMBOL'){
  subcelltype_gseGO_glist=list()
  subcelltype_gseGO_list=list()
  tmp_df=data.frame()
  for (type in unique(RNA_markers$celltype)) {
    print(type)
    subcelltype_gseGO_glist[[type]]=RNA_markers[RNA_markers$celltype==type,'avg_log2FC']
    names(subcelltype_gseGO_glist[[type]])=RNA_markers[RNA_markers$celltype==type,'gene']
    subcelltype_gseGO_glist[[type]]=sort(subcelltype_gseGO_glist[[type]],decreasing = T)
    subcelltype_gseGO_list[[type]] =gseGO(geneList = subcelltype_gseGO_glist[[type]],ont ='all' ,OrgDb = org.Hs.eg.db,keyType = keyType)
    
    if(dim(subcelltype_gseGO_list[[type]]@result)[1] >0){
      subcelltype_gseGO_list[[type]]@result$celltype=type
      tmp_df=rbind(tmp_df,subcelltype_gseGO_list[[type]]@result)}
  }
  return(list(subcelltype_gseGO_list,tmp_df))
}
subcelltype_enrichGO_func=function(RNA_markers){
  subcelltype_enrichGO_glist=list()
  subcelltype_enrichGO_list=list()
  RNA_markers=RNA_markers[RNA_markers$avg_log2FC>0,]
  for (type in unique(RNA_markers$celltype)) {
    cluster=as.character(unique(RNA_markers$cluster))
    subcelltype_enrichGO_glist[[type]][[ cluster[1] ]]=unique(RNA_markers[RNA_markers$cluster==cluster[1] & RNA_markers$celltype==type ,'gene'])
    subcelltype_enrichGO_glist[[type]][[ cluster[2] ]]=unique(RNA_markers[RNA_markers$cluster==cluster[2] & RNA_markers$celltype==type ,'gene'])
    subcelltype_enrichGO_list[[type]][[ cluster[1] ]] =enrichGO(gene =subcelltype_enrichGO_glist[[type]][[ cluster[1] ]] ,OrgDb =org.Hs.eg.db ,keyType = 'SYMBOL',ont = 'all')
    subcelltype_enrichGO_list[[type]][[ cluster[2] ]] =enrichGO(gene =subcelltype_enrichGO_glist[[type]][[ cluster[2] ]] ,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'all')
  }
  
  return(subcelltype_enrichGO_list)
}



#################################################################################################################################################################
#--------------------------------------------------prepare data ----------------------------------------------------#
#################################################################################################################################################################

#--------------------------prepare NRBC data-----------------------#
#--------------------------prepare NRBC data-----------------------#
if(F){
  YS_altas_Ery_seu=readRDS('filt_YS_altas_Ery_seu.rds')
  defintive_Ery_seu=readRDS('definitive_nRBC_marker_project/defintive_filt_NRBC_altas_seu.rds')
  
  
  p=DimPlot(defintive_Ery_seu,group.by = c('tissue_stage','final_celltype'),cols = cols);p
  ggsave(p,filename='res_pic/main_figure2/defintive_NRBC_tissue_celltype.pdf',width = 12,height = 5)
 
  p=DimPlot(subset(NRBC_altas_seu,tissue_stage!='YS'),group.by = c('tissue_stage','final_celltype'),cols = cols);p
  ggsave(p,filename='res_pic/main_figure2/defintive_NRBC_tissue_celltype.pdf',width = 12,height = 5)
  
  
  NRBC_altas_seu=merge(YS_altas_Ery_seu,defintive_Ery_seu);
  NRBC_altas_seu[['RNA']]=JoinLayers(NRBC_altas_seu[['RNA']])
  NRBC_altas_seu[['umap']]=merge(YS_altas_Ery_seu[['ref.umap2']],defintive_Ery_seu[['umap']])
  
  NRBC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('resource')]=as.character(YS_altas_Ery_seu$orig.dataset)
  NRBC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('age')]=as.character(YS_altas_Ery_seu$stage)
  
  NRBC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('tissue_stage')]='YS'
  NRBC_altas_seu$tissue_stage=factor(NRBC_altas_seu$tissue_stage,levels = c('YS','FL','FBM','ABM'))
  
  NRBC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('final_celltype')]=as.character(YS_altas_Ery_seu$main_celltype)
  NRBC_altas_seu$final_celltype=factor(NRBC_altas_seu$final_celltype,levels = c("BFUE/CFUE","ProE",'EryP-CFC','',"Bas","Poly" ,"Orth"))
  
  NRBC_altas_seu$Ery_stage='late_Ery'
  NRBC_altas_seu$Ery_stage[NRBC_altas_seu$final_celltype %in% c('BFUE/CFUE','ProE','EryP-CFC')]='early_Ery'
  NRBC_altas_seu$Ery_stage[NRBC_altas_seu$final_celltype %in% c('Bas')]='mid_Ery'
  NRBC_altas_seu$Ery_stage=factor(NRBC_altas_seu$Ery_stage,levels = c('early_Ery','mid_Ery','late_Ery'))
  
  
  DimPlot(NRBC_altas_seu,group.by = 'final_celltype',split.by = 'tissue_stage',cols = cols,raster=FALSE)
  
  table(NRBC_altas_seu@meta.data[,c('Ery_stage','final_celltype')])
  
  rm(YS_altas_Ery_seu,defintive_Ery_seu);gc()
  
  
  NRBC_altas_seu_meta=NRBC_altas_seu@meta.data
  cho_meta_an_cols=c('orig.ident', 'nCount_RNA', 'nFeature_RNA','age',  'sex', 'id','resource', 'S.Score' ,  'G2M.Score', 'Phase',
                     'tissue_stage', 'final_celltype', 'Ery_stage')
  NRBC_altas_seu_meta=NRBC_altas_seu_meta[,cho_meta_an_cols]
  
  
  NRBC_altas_seu@meta.data=NRBC_altas_seu_meta
  NRBC_altas_seu$source_celltype=paste(NRBC_altas_seu$tissue_stage,NRBC_altas_seu$final_celltype,sep = '_')
  NRBC_altas_seu$source_celltype=factor(NRBC_altas_seu$source_celltype,levels = c( "YS_EryP-CFC","YS_Bas","YS_Poly","YS_Orth", 
                                                                                   "FL_BFUE/CFUE" ,"FL_ProE", "FL_Bas","FL_Poly","FL_Orth",
                                                                                   "FBM_BFUE/CFUE" ,"FBM_ProE","FBM_Bas","FBM_Poly","FBM_Orth", 
                                                                                   "ABM_BFUE/CFUE","ABM_ProE","ABM_Bas","ABM_Poly","ABM_Orth"))
  NRBC_altas_seu$sourceid=paste(NRBC_altas_seu$tissue_stage,NRBC_altas_seu$id,sep = '_')
  Idents(NRBC_altas_seu)='tissue_stage'
  saveRDS(NRBC_altas_seu,file = 'primitive_NRBC_project/all_NRBC_altas_seu.rds')
  
}else{
  NRBC_altas_seu=readRDS('primitive_NRBC_project/all_NRBC_altas_seu.rds')
}



setwd('primitive_NRBC_project')



#################################################################################################################################################################
#--------------------------------------------------figure 1: hematopoietic niche shape the differenation of nRBC---------------------------------------------------#
#################################################################################################################################################################

setwd('primitive_NRBC_project')


# nRBC from different niche, are quite different and 异质性的
dir.create(recursive =T ,'res_pic/main_figure1')
old_var_genes=VariableFeatures(NRBC_altas_seu)
NRBC_altas_seu



p=DimPlot(NRBC_altas_seu,group.by = c('tissue_stage','final_celltype'),cols = cols,raster = F);p
ggsave(p,filename='res_pic/main_figure2/all_NRBC_tissue_celltype.pdf',width = 12,height = 5)

p=pheatmap(cor(data.frame(AverageExpression(NRBC_altas_seu,group.by = 'tissue_stage',layer = 'data',features =VariableFeatures(NRBC_altas_seu) )$RNA)),
           color =  colorRampPalette(colors = c('navy','white','firebrick3'))(100),main = 'based on variable genes')
ggsave(as.ggplot(p),filename = 'res_pic/main_figure2/tissue_basedonvariable_cor_heatmap.pdf',width = 5,height = 5,dpi = 300)

p=pheatmap(cor(data.frame(AverageExpression(NRBC_altas_seu,group.by = 'source_celltype',layer = 'data',features =VariableFeatures(NRBC_altas_seu) )$RNA)),
           colorRampPalette(colors = c('navy','white','firebrick3'))(100),main = 'based on variable genes')
ggsave(as.ggplot(p),filename = 'res_pic/main_figure2/tissue_subcelltype_basedonvariable_cor_heatmap.pdf',width = 8,height =8,dpi = 300)


p=pheatmap(cor(data.frame(AverageExpression(NRBC_altas_seu,group.by = 'source_celltype',layer = 'data',features =rownames(NRBC_altas_seu) )$RNA)),
           colorRampPalette(colors = c('navy','white','firebrick3'))(100),main = 'based on all genes')
ggsave(as.ggplot(p),filename = 'res_pic/main_figure2/tissue_subcelltype_basedonall_cor_heatmap.pdf',width = 8,height =8,dpi = 300)


# NRBC reference marker
p=VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features =c('GATA1','KLF1','TFRC','GYPA','CCNB1','NCL'),stack = T,cols = col)+NoLegend();p
ggsave(p,filename='res_pic/main_figure2/classical_nRBC_marker_expression_in_niches_subceltlype_vlnplot.pdf',width = 8,height = 8)

# 计算不同组织微环境稳定造血时期
NRBC_altas_seu_meta=NRBC_altas_seu@meta.data

YS_CS15_CS18_sample_df=NRBC_altas_seu_meta[NRBC_altas_seu_meta$tissue_stage=='YS' & NRBC_altas_seu_meta$age %in% c('CS15', "CS17" ,"CS18" ),]
YS_celltype_mratio_df=data.frame( colMeans(prop.table(table(YS_CS15_CS18_sample_df[,c('id','final_celltype')]),margin = 1)))
colnames(YS_celltype_mratio_df)='YS'

FL_CS18post_sample_df=NRBC_altas_seu_meta[NRBC_altas_seu_meta$tissue_stage=='FL' & !NRBC_altas_seu_meta$age %in% c('CS14_4PCW', "CS15_5PCW" ,"CS17_6PCW" ),]
FL_celltype_mratio_df=data.frame( colMeans(prop.table(table(FL_CS18post_sample_df[,c('id','final_celltype')]),margin = 1)))
colnames(FL_celltype_mratio_df)='FL'

FBM_sample_df=NRBC_altas_seu_meta[NRBC_altas_seu_meta$tissue_stage=='FBM',]
FBM_celltype_mratio_df=data.frame( colMeans(prop.table(table(FBM_sample_df[,c('id','final_celltype')]),margin = 1)))
colnames(FBM_celltype_mratio_df)='FBM'

ABM_sample_df=NRBC_altas_seu_meta[NRBC_altas_seu_meta$tissue_stage=='ABM',]
table(ABM_sample_df$source)# CD34- CD235+: GSE150774:CD34+,GSE133181,GSE135194,GSE169426,  Mononuclear cells:GSE165645, Mononuclear:CD34+=4:1:GSE181989
ABM_celltype_mratio_CD34pos_df=data.frame( round(colMeans(table(ABM_sample_df[!ABM_sample_df$source %in% c('GSE150774','GSE181989','GSE165645'),c('id','final_celltype')])),digits = 0) )
colnames(ABM_celltype_mratio_CD34pos_df)='mcount'
ABM_celltype_mratio_GYPAposCD34neg_df=data.frame( round(colMeans(table(ABM_sample_df[ABM_sample_df$source=='GSE150774',c('id','final_celltype')])),digits = 0))
colnames(ABM_celltype_mratio_GYPAposCD34neg_df)='mcount'
ABM_celltype_mratio_GYPAposCD34neg_df['BFUE/CFUE','mcount']=0

ABM_NRBC_ratio_df=prop.table(ABM_celltype_mratio_CD34pos_df+ ABM_celltype_mratio_GYPAposCD34neg_df[rownames(ABM_celltype_mratio_CD34pos_df),])
colnames(ABM_NRBC_ratio_df)='ABM'

NRBC_ratio_df= cbind(YS_celltype_mratio_df,cbind(cbind(FL_celltype_mratio_df,FBM_celltype_mratio_df),ABM_NRBC_ratio_df)[rownames(YS_celltype_mratio_df),])
NRBC_ratio_df$celltytpe=rownames(NRBC_ratio_df)
NRBC_ratio_df= melt(NRBC_ratio_df)
NRBC_ratio_df=NRBC_ratio_df[!is.na(NRBC_ratio_df$value),]
NRBC_ratio_df=NRBC_ratio_df[!is.na(NRBC_ratio_df$celltytpe),]

NRBC_ratio_df$celltytpe=factor(NRBC_ratio_df$celltytpe,levels = c("BFUE/CFUE" ,"ProE","EryP-CFC","Bas","Poly","Orth"  ))
p=ggplot(NRBC_ratio_df ,aes(x=variable,y=value,fill=celltytpe))+geom_bar(position = "fill",stat = 'identity')+
  scale_fill_manual(values =  cols)+theme_classic()+theme(axis.text.x = element_text(hjust =1 ,angle = 45,face = 'bold'))
p
ggsave(p,file='res_pic/main_figure2/tissue_celltype_ratio_barplot.pdf',width =3 ,height = 6)


#################################################################################################################################################################
##----------------------------------------------------Figure2 primitive vs definitive nRBC-------------------------# 
#################################################################################################################################################################

NRBC_altas_seu$NRBC_type='definitive'
NRBC_altas_seu$NRBC_type[NRBC_altas_seu$tissue_stage=='YS']='primitive'

Idents(NRBC_altas_seu)='source_celltype'

#------------------------------ 排除nRBC亚类细胞比例影响 ，统一细胞亚类比例-------------------------------------#
if(T){
  # fetal NRBC: YS、FL、FBM，三个阶段各抽取1:1:1的细胞，构成early、mid、late NRBC， 与整体不抽样分析，结果差异很小
  
  YS_early_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='YS',downsample =3000)
  YS_mid_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='YS',downsample =3000)
  YS_late_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='YS',downsample =3000)
  YS_NRBC_altas_seu=merge(YS_early_Ery_NRBC_altas_seu,c(YS_mid_Ery_NRBC_altas_seu,YS_late_Ery_NRBC_altas_seu))
  rm(YS_early_Ery_NRBC_altas_seu,YS_mid_Ery_NRBC_altas_seu,YS_late_Ery_NRBC_altas_seu)
  
  FL_early_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FL',downsample =1000)
  FL_mid_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FL',downsample =1000)
  FL_late_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FL',downsample =1000)
  FL_NRBC_altas_seu=merge(FL_early_Ery_NRBC_altas_seu,c(FL_mid_Ery_NRBC_altas_seu,FL_late_Ery_NRBC_altas_seu))
  rm(FL_early_Ery_NRBC_altas_seu,FL_mid_Ery_NRBC_altas_seu,FL_late_Ery_NRBC_altas_seu)
  
  FBM_early_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FBM',downsample =1000)
  FBM_mid_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FBM',downsample =1000)
  FBM_late_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FBM',downsample =1000)
  FBM_NRBC_altas_seu=merge(FBM_early_Ery_NRBC_altas_seu,c(FBM_mid_Ery_NRBC_altas_seu,FBM_late_Ery_NRBC_altas_seu))
  rm(FBM_early_Ery_NRBC_altas_seu,FBM_mid_Ery_NRBC_altas_seu,FBM_late_Ery_NRBC_altas_seu)
  
  ABM_early_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='ABM',downsample =1000)
  ABM_mid_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='ABM',downsample =1000)
  ABM_late_Ery_NRBC_altas_seu=subset(NRBC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='ABM',downsample =1000)
  ABM_NRBC_altas_seu=merge(ABM_early_Ery_NRBC_altas_seu,c(ABM_mid_Ery_NRBC_altas_seu,ABM_late_Ery_NRBC_altas_seu))
  rm(ABM_early_Ery_NRBC_altas_seu,ABM_mid_Ery_NRBC_altas_seu,ABM_late_Ery_NRBC_altas_seu)
  
  pd_subset_NRBC_altas_seu=merge(ABM_NRBC_altas_seu,c(FBM_NRBC_altas_seu,FL_NRBC_altas_seu,YS_NRBC_altas_seu))
  rm(ABM_NRBC_altas_seu,FBM_NRBC_altas_seu,FL_NRBC_altas_seu,YS_NRBC_altas_seu)
  gc()
  
  pd_subset_NRBC_altas_seu <- JoinLayers(pd_subset_NRBC_altas_seu)
  pd_subset_NRBC_altas_seu$type_stage='definitive'
  pd_subset_NRBC_altas_seu$type_stage[pd_subset_NRBC_altas_seu$tissue_stage=='YS']='primitive'
  Idents(pd_subset_NRBC_altas_seu)='type_stage'
  saveRDS(pd_subset_NRBC_altas_seu,file = 'res_data/temp_pd_subset_NRBC_altas_seu.rds')
  
}else{
  pd_subset_NRBC_altas_seu=readRDS('res_data/temp_pd_subset_NRBC_altas_seu.rds')
  pd_subset_NRBC_altas_seu$source_celltype=factor(pd_subset_NRBC_altas_seu$source_celltype,levels = levels(NRBC_altas_seu$source_celltype))
}

#-------------------------WHOLE LEVEL--------------------------------------#
Idents(pd_subset_NRBC_altas_seu)='type_stage'
pd_whole_level_markers=FindAllMarkers(pd_subset_NRBC_altas_seu)
pd_whole_level_markers=pd_whole_level_markers[pd_whole_level_markers$p_val_adj <0.05,]
saveRDS(pd_whole_level_markers,file = 'res_data/pd_whole_level_markers.rds')

temp_df=pd_whole_level_markers[pd_whole_level_markers$avg_log2FC >1 & pd_whole_level_markers$pct.2<0.2  ,]
top_pd_whole_level_markers=temp_df[-grep('AC09',temp_df$gene),] %>% group_by(cluster) %>% do(head(.,10))
p=DotPlot(NRBC_altas_seu,features = c(top_pd_whole_level_markers$gene[11:20],top_pd_whole_level_markers$gene[1:10]),cols = c('gray','firebrick3'),scale = F,group.by = 'source_celltype')+RotatedAxis() ;p
ggsave(p,filename = 'res_pic/main_figure2/pd_top10_marker_dotplot.pdf',width = 8,height = 10)


whole_pd_enrichgo_res=compareCluster(geneClusters =list('primitive'=unique(pd_whole_level_markers[pd_whole_level_markers$cluster=='primitive' & pd_whole_level_markers$avg_log2FC >0,'gene']),
                                                        'definitive'=unique(pd_whole_level_markers[pd_whole_level_markers$cluster=='definitive' & pd_whole_level_markers$avg_log2FC>0,'gene'])),
                                     keyType = 'SYMBOL',ont = "BP",fun = 'enrichGO',  OrgDb='org.Hs.eg.db' )

saveRDS(whole_pd_enrichgo_res,file = 'res_data/whole_tissue_pd_NRBC_marker_enrichgo_res1.rds')




pd_enrichgo_res=pairwise_termsim(whole_pd_enrichgo_res)

# ggplot美化, 不同niche来源nRBC具有不同的特征状态
top10_enrichgo_res_df=pd_enrichgo_res@compareClusterResult %>% group_by(Cluster) %>% do(head(.,15))
top10_enrichgo_res_df$ratio=as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[1,])/as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[2,])
top10_enrichgo_res_df$Description=factor(top10_enrichgo_res_df$Description,levels =unique(top10_enrichgo_res_df$Description) )
p=ggplot(top10_enrichgo_res_df,aes(x=Cluster,y=Description,color=-log10(p.adjust),size=Count))+geom_point()+theme_bw()+
  scale_color_gradient(low = '#4387B5',high = 'firebrick3')+theme(text = element_text(face = 'bold'))
p 
ggsave(p,filename = 'res_pic/main_figure2/pd_whole_enrichGO_dotplot.pdf',width = 8,height = 8)

emapplot(pd_enrichgo_res, node_label = "group")
p=emapplot(pd_enrichgo_res, node_label = "category",showCategory =30);p
ggsave(p,filename = 'res_pic/main_figure2/pd_whole_enrichGO_emapplot.pdf',width = 8,height = 8)

genelist=pd_whole_level_markers[pd_whole_level_markers$cluster=='primitive','avg_log2FC'];names(genelist)=pd_whole_level_markers[pd_whole_level_markers$cluster=='primitive','gene']
genelist=sort(genelist,decreasing = T)
pd_whole_level_marker_gseGO_res=data.frame(gseGO(geneList =genelist,ont = 'BP',OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',pvalueCutoff = 0.05 ))
pd_whole_level_marker_gseGO_res$res='up';pd_whole_level_marker_gseGO_res$res[pd_whole_level_marker_gseGO_res$NES<0]='down'
table(pd_whole_level_marker_gseGO_res$res)#down/up: 188   61 
saveRDS(pd_whole_level_marker_gseGO_res,file = 'res_data/pd_whole_level_marker_gseGO_res.rds')


top_pd_whole_level_marker_gseGO_res=pd_whole_level_marker_gseGO_res %>% group_by( res)  %>% do(head(.,15))
top_pd_whole_level_marker_gseGO_res=top_pd_whole_level_marker_gseGO_res[order(top_pd_whole_level_marker_gseGO_res$res,top_pd_whole_level_marker_gseGO_res$NES,decreasing = T),]
top_pd_whole_level_marker_gseGO_res$Description=factor(top_pd_whole_level_marker_gseGO_res$Description,levels = unique(top_pd_whole_level_marker_gseGO_res$Description))
p1=ggplot(top_pd_whole_level_marker_gseGO_res,aes(x =res ,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+
  scale_color_gradient2(low = 'navy',mid = 'gray',high = 'firebrick3')+
  theme(text = element_text(face = 'bold'))+ggtitle("primitive vs definitive  in whole level")
p1
ggsave(p1,filename = 'res_pic/main_figure2/pd_whole_deg_gsego_dotplot.pdf',width = 6,height = 8)



#-------------------------substage LEVEL--------------------------------------#
group='primitive_definitive'
sfile='res_data/primitive_definitive_all_Ery_RNA_markers.csv'
res=find_mDEGs_func(seu = pd_subset_NRBC_altas_seu,group = group,sfile = sfile)
sub_pd_all_Ery_tissue_markers=res[[1]]
sub_pd_count_df=res[[2]]
rm(res);gc()

sub_pd_all_Ery_tissue_markers=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$p_val_adj <0.05,]
length(unique(sub_pd_all_Ery_tissue_markers$gene) )
table(sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >0,c('celltype','cluster')])
length(unique(sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >0 & sub_pd_all_Ery_tissue_markers$cluster=='primitive' ,'gene']))
length(unique(sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$avg_log2FC >0 & sub_pd_all_Ery_tissue_markers$cluster=='definitive' ,'gene']))

table(unique(sub_pd_all_Ery_tissue_markers$gene) %in% unique(pd_whole_level_markers$gene)) 
#FALSE  TRUE 
#1740  7222 
table(unique(pd_whole_level_markers$gene) %in% unique(sub_pd_all_Ery_tissue_markers$gene)) 
#FALSE  TRUE 
#2806  7222 

# 查看在亚群中鉴定到的新marker中top marker 表达情况，效果比较差
new_sub_pd_all_Ery_tissue_markers=sub_pd_all_Ery_tissue_markers[!sub_pd_all_Ery_tissue_markers$gene %in% unique(pd_whole_level_markers$gene) & sub_pd_all_Ery_tissue_markers$avg_log2FC >0 ,]
new_sub_pd_all_Ery_tissue_markers=new_sub_pd_all_Ery_tissue_markers[new_sub_pd_all_Ery_tissue_markers$avg_log2FC >1 & new_sub_pd_all_Ery_tissue_markers$pct.1 >0.1 & new_sub_pd_all_Ery_tissue_markers$pct.2 <0.3,]
length(unique(new_sub_pd_all_Ery_tissue_markers$gene))# 600 
top5_new_sub_pd_all_Ery_tissue_markers=new_sub_pd_all_Ery_tissue_markers %>% group_by(celltype) %>%top_n(wt = avg_log2FC,n = 5) 
VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features =top5_new_sub_pd_all_Ery_tissue_markers$gene,stack = T)+NoLegend()

#查看在整体水平中鉴定到的marker，而在亚群中未鉴定到，中top marker 表达情况，效果比较好,ABM nBRC pct低，而 表达大水平接近
out_pd_whole_level_markers=pd_whole_level_markers[!pd_whole_level_markers$gene %in% unique(sub_pd_all_Ery_tissue_markers$gene),]
table(out_pd_whole_level_markers$avg_log2FC >1 &  out_pd_whole_level_markers$pct.2 <0.2 & out_pd_whole_level_markers$pct.1 >0.1)
#FALSE  TRUE 
#5582         30 ,几乎不存在特异性，除了KRCC1， 候选筛选还是应该考虑sub_pd_all_Ery_tissue_markers
out_pd_whole_level_markers=out_pd_whole_level_markers[out_pd_whole_level_markers$avg_log2FC >1 &  out_pd_whole_level_markers$pct.2 <0.3 & out_pd_whole_level_markers$pct.1 >0.1,]
VlnPlot(pd_subset_NRBC_altas_seu,group.by = 'source_celltype',features =out_pd_whole_level_markers$gene,stack = T,cols = cols)
DotPlot(pd_subset_NRBC_altas_seu,group.by = 'source_celltype',features =out_pd_whole_level_markers$gene,scale = T)+RotatedAxis()

sub_pd_all_Ery_tissue_markers_gseGO_res=subcelltype_gseGO_func(RNA_markers = sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$cluster=='primitive',])
sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$res='up';sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$res[sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$NES <0]='down'
saveRDS(sub_pd_all_Ery_tissue_markers_gseGO_res,file = 'res_data/sub_pd_all_Ery_tissue_markers_gseGO_res.rds')

top_sub_pd_all_Ery_tissue_markers_gseGO_res=sub_pd_all_Ery_tissue_markers_gseGO_res[[2]][sub_pd_all_Ery_tissue_markers_gseGO_res[[2]]$ONTOLOGY=='BP',] %>% group_by(celltype,res) %>% do(head(.,10))
top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype=factor(top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))
top_sub_pd_all_Ery_tissue_markers_gseGO_res=top_sub_pd_all_Ery_tissue_markers_gseGO_res[order(top_sub_pd_all_Ery_tissue_markers_gseGO_res$celltype,top_sub_pd_all_Ery_tissue_markers_gseGO_res$res,top_sub_pd_all_Ery_tissue_markers_gseGO_res$NES),]
top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description=factor(top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description,levels = unique(top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description))
p=ggplot(top_sub_pd_all_Ery_tissue_markers_gseGO_res,aes(x =celltype ,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+
  scale_color_gradient2(low = 'navy',mid = 'gray',high = 'firebrick3')+theme(text = element_text(face = 'bold'))+ggtitle(' primitive vs definitive in substage')+RotatedAxis()
#scale_y_discrete(label=function(x)str_wrap(string = x,width = 50))

p
ggsave(p,filename = 'res_pic/main_figure2/pd_erysubstage_gsetop10_dotplot.pdf',width = 8,height = 12)


# 利用emapplot 查看top enriched pathway 之间的关系
test=pairwise_termsim(sub_pd_all_Ery_tissue_markers_gseGO_res[[1]]$early_Ery)
emapplot(test,color='NES',showCategory=30)+ggtitle('pd_early_Ery_top30')
test=pairwise_termsim(sub_pd_all_Ery_tissue_markers_gseGO_res[[1]]$mid_Ery)
emapplot(test,color='NES',showCategory=30)+ggtitle('pd_mid_Ery_top30')
test=pairwise_termsim(sub_pd_all_Ery_tissue_markers_gseGO_res[[1]]$late_Ery)
emapplot(test,color='NES')+ggtitle('pd_late_Ery_top30')


# 比较enrichGO 得到的结果
pd_all_Ery_tissue_markers_enrichGO_list=subcelltype_enrichGO_func(RNA_markers = sub_pd_all_Ery_tissue_markers[ sub_pd_all_Ery_tissue_markers$avg_log2FC>0,] )#no term enriched under specific pvalueCutoff.

pd_all_Ery_tissue_markers_enrichGO_df=data.frame()
for (celltype in c('early_Ery','mid_Ery','late_Ery')) {
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['primitive']]@result$celltype=celltype
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['definitive']]@result$celltype=celltype
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['primitive']]@result$type='primitive'
  pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['definitive']]@result$type='definitive'
  temp_df=rbind(pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['primitive']]@result,pd_all_Ery_tissue_markers_enrichGO_list[[celltype]][['definitive']]@result)
  pd_all_Ery_tissue_markers_enrichGO_df=rbind(pd_all_Ery_tissue_markers_enrichGO_df,temp_df)
}

pd_all_Ery_tissue_markers_enrichGO_df$ratio=as.numeric(t(data.frame(strsplit(pd_all_Ery_tissue_markers_enrichGO_df$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit( pd_all_Ery_tissue_markers_enrichGO_df$GeneRatio,split = '/')))[,2])
pd_all_Ery_tissue_markers_enrichGO_df$celltype=factor(pd_all_Ery_tissue_markers_enrichGO_df$celltype,levels =c('early_Ery','mid_Ery','late_Ery') )
saveRDS(pd_all_Ery_tissue_markers_enrichGO_df,file = 'res_data/pd_all_Ery_tissue_markers_enrichGO_df.rds')

top_pd_all_Ery_tissue_markers_enrichGO_df=pd_all_Ery_tissue_markers_enrichGO_df %>%group_by(type,celltype) %>% do(head(.,10))
top_pd_all_Ery_tissue_markers_enrichGO_df=top_pd_all_Ery_tissue_markers_enrichGO_df[order(top_pd_all_Ery_tissue_markers_enrichGO_df$type,top_pd_all_Ery_tissue_markers_enrichGO_df$celltype),]
top_pd_all_Ery_tissue_markers_enrichGO_df$Description =factor(top_pd_all_Ery_tissue_markers_enrichGO_df$Description,levels = unique(top_pd_all_Ery_tissue_markers_enrichGO_df$Description))
p=ggplot(top_pd_all_Ery_tissue_markers_enrichGO_df,aes(x=type,y=Description,color=-log10(p.adjust),size=ratio))+geom_point()+theme_bw()+
  scale_color_gradient2(low = 'navy',mid='gray',high = 'firebrick3')+facet_grid(~celltype)+
  theme( axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'enrichGO')#+scale_y_discrete(label=function(x)str_wrap(x,width = 50))
p
ggsave(p,filename='res_pic/main_figure2/pd_erystage_enrichtop10_dotplot.pdf',width = 8,height = 12)


# 逻辑假设：强代谢、强金属离子压力以及ROS---->加速功能成熟和分化导致

############################################################################################################################################
#--------------------------------------character of Primitive NRBC-------------------------------#
############################################################################################################################################

# marker 表现明显的压力：金属离子，ROS：antioxidant enzymes(FAM213A,SOD3),the stress marker GDF15
copper_metabolism_core <- c( 'SLC31A1', "ATOX1", "CCS", 'COX17', "SCO1")  # 胞内运输
VlnPlot(NRBC_altas_seu,features =copper_metabolism_core,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()

metal_detox_genes <- c(
  "ATOX1", "CCS",'SCO1', # 铜离子
  "SLC39A8",        # 离子转运
  "CALM1", "CALR", # 钙信号/结合
  'SLC40A1', # 铁例子转出
  'MTF1', # 金属离子应激转录因子
  "SOD1",'GPX1', "GPX4", "CAT"   # 抗氧化
)
p=VlnPlot(NRBC_altas_seu,features =metal_detox_genes,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend();p
ggsave(p,filename='res_pic/main_figure2/metal_antioxidant_gene_lvnplot.pdf',width = 7,height = 8)

# MTF1	经典金属响应转录因子，很少细胞表达
DotPlot(NRBC_altas_seu,features=metal_detox_genes ,group.by = 'source_celltype',scale = F,cols =c( "gray", "firebrick3") )+RotatedAxis()


# ----------------------------------- cell cycle analysis-----------------------------------#
YS_id_phase_ratio=prop.table(table(NRBC_altas_seu@meta.data[NRBC_altas_seu$tissue_stage=='YS' & NRBC_altas_seu$final_celltype %in% c(''),c('sourceid','Phase')]),margin = 1)

p=as.ggplot(pheatmap(t(prop.table(table(NRBC_altas_seu@meta.data[,c('source_celltype','Phase')]),margin = 1)),display_numbers = T,cluster_rows = F,cluster_cols = F,main = 'ratio of Phase',
                     color=colorRampPalette(colors = brewer.pal(11, "PiYG")[c(8,6,3)])(100)))
p
ggsave(p,width =8 ,height =4,filename='res_pic/main_figure2/niche_NRBC_altas_CellCycle_phase_ratio_heatmap.pdf' )

p=DimPlot(NRBC_altas_seu,reduction = 'umap',group.by = 'Phase',raster = F);p
ggsave(p,width =6 ,height =6,filename='res_pic/main_figure2/NRBC_CellCycle_phase_umap.pdf' )

filt_YS_altas_Ery_seu=readRDS('filt_YS_altas_Ery_seu.rds')
temp_df=t(prop.table(table(filt_YS_altas_Ery_seu@meta.data[filt_YS_altas_Ery_seu$main_celltype=='Bas',c('celltype','Phase')]),margin = 1))
temp_df
p=as.ggplot(pheatmap(temp_df[,2:4],display_numbers = T,cluster_rows = F,cluster_cols = F,main = 'ratio of Phase',
                     color=colorRampPalette(colors = brewer.pal(11, "PiYG")[c(8,6,3)])(100)))
p
ggsave(p,width =4 ,height =3,filename='res_pic/main_figure2/NRBC_Bas_CellCycle_phase_umap.pdf' )


# DNA 复制
#PCNA	DNA 复制滑动夹	增殖最经典、最广泛使用的标志物
#MCM2	DNA 复制许可因子	代表细胞处于增殖周期或具备增殖能力
#CDK1	细胞周期核心激酶	代表 G2/M 期推进，真正进入分裂
#细胞周期退出	CDKN1A、CDKN1B、CDKN2C ，抑制增殖BTG1 / BTG2 / BTG3
proliferation_regulate_markers <- c("PCNA", "MCM2", "CDK1",'MKI67','CDKN1A','CDKN1B','BTG1' , 'BTG2' , 'BTG3') 

p=VlnPlot(NRBC_altas_seu,features =proliferation_regulate_markers,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend();p
ggsave(p,width =6,height =6,filename='res_pic/main_figure2/NRBC_CellCycle_gene_vlnplot.pdf' )


# -----------------------RNA and Protein metabolism-------------------------------------#
#
RNA_spicing_core <- c( "ZMAT2", "PQBP1", "RALY", "RBM42", "CIRBP")   # RNA 剪接
VlnPlot(NRBC_altas_seu,features =RNA_spicing_core,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()

#ARF1	高尔基体 COPI 出芽	高尔基体运输的起始
#ARF5	高尔基体运输	高尔基体运输执行
#RAB6A	高尔基体 → ER 逆向运输	高尔基体回路完整
#RAB7A	晚期内体 → 溶酶体	降解方向
#RAB11B	循环内体 → 质膜	回收方向
#EHD1	循环内体 → 质膜	回收执行
#这 3 个基因代表内体系统的两大去向：降解和循环。
#DNM2	膜分裂	所有囊泡形成都需要
#NAPA	SNARE 融合	所有囊泡融合都需要
#这 2 个基因是通用执行元件，放在最后。
golgi_centered_trafficking <- c("ARF1", "ARF5", "RAB6A", "RAB7A", "RAB11B", "EHD1",  "DNM2", "NAPA" )
VlnPlot(NRBC_altas_seu,features =golgi_centered_trafficking,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()

# Ubiquitin–proteasome system：PSMB5、PSMD1、PSMD11、PSMC4
# Proteostasis regulation	NFE2L1、HSPA8、DNAJB1	蛋白酶体表达调节及分子伴侣相关质量控制
Ubiquitin_proteasome_genes
VlnPlot(NRBC_altas_seu,features =Ubiquitin_proteasome_genes,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()

p=VlnPlot(NRBC_altas_seu,features =c(RNA_spicing_core,golgi_centered_trafficking,Ubiquitin_proteasome_genes),group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()
p
ggsave(p,width =6,height =6,filename='res_pic/main_figure2/NRBC_CellCycle_gene_vlnplot.pdf' )

cho_kn_genes=c(metal_detox_genes,proliferation_regulate_markers,RNA_spicing_core,golgi_centered_trafficking,Ubiquitin_proteasome_genes)
an_df2=data.frame(row.names = cho_kn_genes,type=c(rep('metal_response',length(metal_detox_genes)),rep('proliferation_regulator',length(proliferation_regulate_markers)),rep('RNA_spicing',length(RNA_spicing_core)),rep('golgi_centered_trafficking',length(golgi_centered_trafficking)),
                                                   rep('proteasome_genes',length(Ubiquitin_proteasome_genes)) ) )

temp_df=pd_whole_level_markers[pd_whole_level_markers$gene %in% cho_kn_genes & pd_whole_level_markers$cluster=='primitive',]
rownames(temp_df)=temp_df$gene
an_df2$global_avg_log2FC=temp_df[rownames(an_df2),'avg_log2FC']

temp_df=sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$gene %in% cho_kn_genes & sub_pd_all_Ery_tissue_markers$cluster=='primitive',]
early_temp_df=temp_df[temp_df$celltype=='early_Ery',]
rownames(early_temp_df)=early_temp_df$gene
an_df2$early_avg_log2FC=early_temp_df[rownames(an_df2),'avg_log2FC']
mid_temp_df=temp_df[temp_df$celltype=='mid_Ery',]
rownames(mid_temp_df)=mid_temp_df$gene
an_df2$mid_avg_log2FC=mid_temp_df[rownames(an_df2),'avg_log2FC']
late_temp_df=temp_df[temp_df$celltype=='late_Ery',]
rownames(late_temp_df)=late_temp_df$gene
an_df2$late_avg_log2FC=late_temp_df[rownames(an_df2),'avg_log2FC']
an_df2[is.na(an_df2)]=0

ann_colors <- list(
  type=setNames(col[c(1,3:6)], c('metal_response','proliferation_regulator',"RNA_spicing", "golgi_centered_trafficking", "proteasome_genes")),
  global_avg_log2FC = colorRampPalette(c( "white", "firebrick3"))(100),
  early_avg_log2FC = colorRampPalette(c( "white", "firebrick3"))(100),
  mid_avg_log2FC = colorRampPalette(c( "white", "firebrick3"))(100),
  late_avg_log2FC = colorRampPalette(c( "white", "firebrick3"))(100)
)
p=pheatmap(t(AverageExpression(NRBC_altas_seu,layer = 'data',features =cho_kn_genes ,group.by = 'source_celltype')$RNA[cho_kn_genes,]),
           cluster_rows = F,cluster_cols = F,
           annotation_col= an_df2,annotation_colors = ann_colors,colorRampPalette(colors = c('white','firebrick3'))(100))
p
ggsave(as.ggplot(p),filename='res_pic/main_figure2/metabolism_genes_heatmap.pdf',width = 8,height = 6)

p=DotPlot(NRBC_altas_seu,features=cho_kn_genes ,group.by = 'source_celltype',scale = F,cols =c( "gray", "firebrick3") )+RotatedAxis()
ggsave(p,filename='res_pic/main_figure2/metabolism_genes_heatmap_dotplot.pdf',width = 16,height = 6)


VlnPlot(NRBC_altas_seu,features=cho_kn_genes ,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols,stack = T)+RotatedAxis()


an_df1=data.frame(row.names = cho_kn_genes,type=c(rep('metal_response',length(metal_detox_genes)),rep('proliferation_regulator',length(proliferation_regulate_markers)),rep('RNA_spicing',length(RNA_spicing_core)),rep('golgi_centered_trafficking',length(golgi_centered_trafficking)),
                                                  rep('proteasome_genes',length(Ubiquitin_proteasome_genes)) ) )
an_df1$type=factor(an_df1$type,levels = unique(an_df1$type))
p=pheatmap(t(an_df2[,-1]),cluster_rows = F,annotation_col =an_df1 ,cluster_cols = F,color =  colorRampPalette(c( 'navy',"white", "firebrick3"))(100),display_numbers = T)
ggsave(as.ggplot(p),filename='res_pic/main_figure2/metabolism_genes_heatmap_FC.pdf',width = 16,height = 4)


############################################################################################################################################
# ----------------------------------- Figure 3 fast maturation：提前功能实现，主要用词-----------------------------------#
# 提前退出细胞周期，提前表达分化基因，提前表达红细胞功能基因，提前启动线粒体自噬
############################################################################################################################################


# Erythroid functional maturation	ALAS2、AHSP、SLC4A1、GYPA、ANK1、EPB42、SLC25A37、FECH
#血红素合成与铁代谢 ALAS2 · SLC25A37 · FECH
#珠蛋白稳定与蛋白质稳态:AHSP
#红细胞膜结构与气体运输:SLC4A1 · GYPA · ANK1 · EPB42
maturation_differentiation_genes=c( "HBA1", "HBA2",'AHSP','ALAS2','FECH','HMBS','SLC25A37','SLC4A1','ANK1','EPB42') 
p=VlnPlot(NRBC_altas_seu,features =maturation_differentiation_genes,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()
p

NRBC_altas_seu1=NRBC_altas_seu
NRBC_altas_seu1$source_celltype=factor(NRBC_altas_seu$source_celltype,levels = c('YS_EryP-CFC','FL_BFUE/CFUE','FBM_BFUE/CFUE','ABM_BFUE/CFUE',
                                                                                 'FL_ProE','FBM_ProE','ABM_ProE'))
p=VlnPlot(subset(NRBC_altas_seu1,final_celltype %in% c( "BFUE/CFUE", "EryP-CFC",'ProE')),features =maturation_differentiation_genes,
          group.by = 'source_celltype',stack = T,cols = cols,split.by = 'tissue_stage' ) +NoLegend()
p

Idents(NRBC_altas_seu)='final_celltype'
temp_df=FindMarkers(NRBC_altas_seu,ident.1 = 'EryP-CFC' ,ident.2 = 'BFUE/CFUE')
cho_temp_df=temp_df[ maturation_differentiation_genes,]
cho_temp_df$gene=rownames(cho_temp_df)
cho_temp_df$gene=factor(cho_temp_df$gene,levels =maturation_differentiation_genes )

p1=ggplot(cho_temp_df,aes(x=gene,y=avg_log2FC))+geom_bar(fill=alpha('firebrick3',0.8),stat = 'identity')+theme_classic()

p=p1/p+plot_layout(heights = c(1, 2),)
ggsave(p,width =8 ,height =6,filename='res_pic/main_figure3/NRBC_function_maturation_gene_vlnplot.pdf' )

temp_df2=FindMarkers(NRBC_altas_seu,ident.1 = 'EryP-CFC' ,ident.2 = 'ProE')
cho_temp_df2=temp_df2[ maturation_differentiation_genes,]
cho_temp_df2$gene=rownames(cho_temp_df2)
cho_temp_df2$gene=factor(cho_temp_df2$gene,levels =maturation_differentiation_genes )
p2=ggplot(cho_temp_df2,aes(x=gene,y=avg_log2FC))+geom_bar(fill=alpha('firebrick3',0.8),stat = 'identity')+theme_classic()
p2

p=p2/p1/p+plot_layout(heights = c(1,1,3))
p
ggsave(p,width =8 ,height =10,filename='res_pic/main_figure3/NRBC_function_maturation_gene_vlnplot.pdf' )

# ISR
# stress 调控：抗氧化响应，金属离子响应，整合应激反应（Integrated Stress Response，ISR）？
#ATF4	整合应激反应（ISR）核心转录因子
#NFE2L2（NRF2）	抗氧化响应
#ISR 由 4 种 eIF2α 激酶启动，它们感知不同压力：
#激酶	全名	感知的应激
#EIF2AK1	HRI	血红素缺乏、氧化应激、红系特异压力
#EIF2AK2	PKR	双链 RNA、病毒感染
#EIF2AK3	PERK	内质网应激、蛋白折叠压力
#EIF2AK4	GCN2	氨基酸饥饿、空载 tRNA


# ASNS和DDIT3作为ISR下游响应基因
#ASNS	ATF4介导的氨基酸代谢适应
#PPP1R15A	ISR负反馈
#TRIB3	ATF4相关应激反馈
#SLC7A11	氨基酸转运及抗氧化适应
#DDIT3	ISR及持续应激相关响应
#CHAC1	应激相关谷胱甘肽代谢
#ATF5	应激转录调控
#DDIT4	应激相关生长调节
#EIF2AK1	红系HRI通路上游激酶
#TXNRD1	氧化还原稳态
#HSPD1	线粒体蛋白稳态
#HSPA9	线粒体蛋白折叠与红系功能
#应激相关	DDIT4、SESN2、MTOR、EIF4EBP1
VlnPlot(NRBC_altas_seu,features =c('DDIT4','SESN2','MTOR','EIF4EBP1'),group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()
VlnPlot(NRBC_altas_seu,features =c('ATF4','ASNS', 'DDIT3','EIF2AK1','EIF2AK2','EIF2AK3','EIF2AK4'),group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()

p=VlnPlot(NRBC_altas_seu,features =c('ATF4','EIF2AK1','ASNS', 'DDIT3','TRIB3','PPP1R15A','CHAC1','DDIT4','EIF4EBP1'),group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()
p

cho_temp_df=pd_whole_level_markers[ pd_whole_level_markers$gene %in% c('ATF4','EIF2AK1','ASNS', 'DDIT3','TRIB3','PPP1R15A','CHAC1','DDIT4','EIF4EBP1'),]
cho_temp_df=cho_temp_df[cho_temp_df$cluster=='primitive',]
rownames(cho_temp_df)=cho_temp_df$gene
cho_temp_df=cho_temp_df[c('ATF4','EIF2AK1','ASNS', 'DDIT3','TRIB3','PPP1R15A','CHAC1','DDIT4','EIF4EBP1'),]
cho_temp_df[3,'gene']='ASNS'
rownames(cho_temp_df)=cho_temp_df$gene
cho_temp_df[is.na(cho_temp_df)]=0
cho_temp_df$gene=factor(cho_temp_df$gene,levels = cho_temp_df$gene)
p1=ggplot(cho_temp_df,aes(x=gene,y=avg_log2FC))+geom_bar(fill=alpha('firebrick3',0.8),stat = 'identity')+theme_classic()

p=p1/p+plot_layout(heights = c(1, 3))
p
ggsave(p,width =8 ,height =8,filename='res_pic/main_figure3/NRBC_ISR_gene_vlnplot.pdf' )

# early stage
sub_pd_all_Ery_tissue_markers[sub_pd_all_Ery_tissue_markers$gene=='ASNS' & sub_pd_all_Ery_tissue_markers$celltype=='early_Ery' ,][2,]

# NRBC最后阶段发生线粒体自噬清除
#BNIP3、BNIP3L： 线粒体自噬受体，负责识别并招募自噬相关机器。
#SQSTM1： 选择性自噬适配蛋白，连接泛素化底物与自噬体。
#MAP1LC3B： 自噬体膜蛋白，参与自噬底物的包裹。
mitophagy <- c("BNIP3L", "BNIP3", "SQSTM1", "MAP1LC3B") # 线粒体自噬相关转录程序增强
p=VlnPlot(NRBC_altas_seu,features =mitophagy,group.by = 'source_celltype',stack = T,cols = cols ) +NoLegend()
p


cho_temp_df=pd_whole_level_markers[ pd_whole_level_markers$gene %in%  c("BNIP3L", "BNIP3", "SQSTM1", "MAP1LC3B"),]
cho_temp_df=cho_temp_df[cho_temp_df$cluster=='primitive',]
rownames(cho_temp_df)=cho_temp_df$gene
cho_temp_df=cho_temp_df[ c("BNIP3L", "BNIP3", "SQSTM1", "MAP1LC3B"),]
cho_temp_df$gene=factor(cho_temp_df$gene,levels = cho_temp_df$gene)
p1=ggplot(cho_temp_df,aes(x=gene,y=avg_log2FC))+geom_bar(fill=alpha('firebrick3',0.8),stat = 'identity')+theme_classic()

p=p1/p+plot_layout(heights = c(1, 2))
p

ggsave(p,width =6 ,height =8,filename='res_pic/main_figure2/NRBC_mitophagy_gene_vlnplot.pdf' )



cho_kn_genes2=c(hemoglobin_synthesis_core,Fe_iron_metabolism_core,known_prodifferention)
an_df2=data.frame(row.names = cho_kn_genes2,type=c(rep('hemoglobin_synthesis',length(hemoglobin_synthesis_core)),rep('Fe_iron_metabolism',length(Fe_iron_metabolism_core)),
                                                   rep('prodifferention',length(known_prodifferention)) ) )
p=pheatmap(AverageExpression(NRBC_altas_seu,layer = 'data',features =cho_kn_genes2 ,group.by = 'source_celltype')$RNA[cho_kn_genes2,],cluster_rows = F,cluster_cols = F,
           annotation_row = an_df2,colorRampPalette(colors = c('navy','white','firebrick3'))(100))
ggsave(as.ggplot(p),filename='res_pic/main_figure2/regulte_stress_genes_heatmap.pdf',width = 8,height = 8)


# 
# 血红素中铁转运相关基因
NRBC_altas_seu=readRDS('~/2021_NRBC_chlyu/primitive_NRBC_project/all_NRBC_altas_seu.rds')

#功能	代表基因
#铁摄取	TFRC、SLC11A2
#铁储存与调节	FTH1、FTL
#铁输出	SLC40A1
#线粒体铁转运	SLC25A37、SLC25A28
#血红素合成	ALAS2、FECH、UROS、UROD
Fe_key_genes=c('TFRC','SLC11A2','FTH1','FTL','SLC40A1')
VlnPlot(NRBC_altas_seu,features = Fe_key_genes,group.by = 'source_celltype',stack = T,cols = cols)+NoLegend()

###########################################################################################################################
####------------- pathway activity score analysis-----------------#
###########################################################################################################################
library(UCell)


library(GO.db)
library(AnnotationDbi)


#cellular oxidant detoxification	GO:0098869	清除细胞内超氧阴离子和过氧化氢
#Response to oxidative stress	GO:0006979	对氧化应激的响应
#Cell redox homeostasis	GO:0045454	维持细胞氧化还原稳态

cho_gene_list=mapIds(x = org.Hs.eg.db,keys =c( 'GO:0098869','GO:0006979','GO:0045454'),column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
names(cho_gene_list)=c('cellular_oxidant_detoxification','Response_to_oxidative_stress','Cell_redox_homeostasis')
cho_gene_list=sapply(cho_gene_list, unique)
gene_sets=c(cho_gene_list)
set.seed(123)
NRBC_altas_seu <- AddModuleScore_UCell(NRBC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快


# response to metal ion :GO:0010038
# stress response to metal ion: GO:0097501 
cho_gene_list=mapIds(x = org.Hs.eg.db,keys =c( 'GO:0010038','GO:0097501'),column = 'SYMBOL',keytype = 'GOALL',multiVals = 'list') # gene limited 
names(cho_gene_list)=c('response_to_metal_ion','stress_response_to_metal_ion')
gene_sets=c(cho_gene_list)
NRBC_altas_seu <- AddModuleScore_UCell(NRBC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快


#Protein homeostasis：Proteasome-mediated ubiquitin-dependent protein catabolic process	GO:0043161	蛋白酶体介导的泛素依赖性蛋白质降解
# RNA homeostasis ： RNA splicing	GO:0008380	RNA 剪接
# Intracellular trafficking / vesicular transport： Golgi vesicle transport	GO:0048193	高尔基体囊泡运输

cho_gene_list=mapIds(x = org.Hs.eg.db,keys =c( 'GO:0043161','GO:0008380','GO:0048193'),column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
names(cho_gene_list)=c('Proteasome-mediated_protein_catabolic_process','RNA_splicing','Golgi_vesicle_transport')
gene_sets=c(cho_gene_list)
NRBC_altas_seu <- AddModuleScore_UCell(NRBC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快


VlnPlot(NRBC_altas_seu,features = c('cellular_oxidant_detoxification_UCell', 'Response_to_oxidative_stress_UCell',
                                    'Proteasome-mediated_protein_catabolic_process_UCell'),
        group.by = 'final_celltype',split.by = 'tissue_stage',cols = cols,stack = T)


#GO:0140467 — integrated stress response signaling：用于概括整体ISR相关转录特征。
#GO:0140468 — HRI-mediated signaling：用于探索与红系血红素稳态相关的特定分支。
# GO:0039585 — PKR/eIFalpha signaling 更适合作为补充分析，尤其在EIF2AK2及其相关基因确实表现出差异时。无 gene symbol 信息
cho_gene_list=mapIds(x = org.Hs.eg.db,keys =c( 'GO:0140467','GO:0140468'),column = 'SYMBOL',keytype = 'GOALL',multiVals = 'list')
names(cho_gene_list)=c('integrated_stress_response_signaling','HRI_mediated_signaling')
gene_sets=c(cho_gene_list)
set.seed(123)
NRBC_altas_seu <- AddModuleScore_UCell(NRBC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快


#negative regulation of erythrocyte differentiation	GO:0045649
# heme biosynthetic process	GO:0006783
# hemoglobin biosynthetic process	GO:0042541
# oxygen transport	GO:0015671
#GO ID	GO Biological Process	适用方向
#GO:0030218	erythrocyte differentiation	红系分化
#GO:0043249	erythrocyte maturation	红细胞成熟

#膜骨架组织与重塑
#GO:0030036
#actin cytoskeleton organization
#GO:0007009
#plasma membrane organization
#GO:0030866
#cortical actin cytoskeleton organization


cho_gene_list=mapIds(x = org.Hs.eg.db,keys =c( 'GO:0030218','GO:0006783','GO:0042541','GO:0015671','GO:0030866'),column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
names(cho_gene_list)=c('erythrocyte_differentiation','heme_biosynthetic_process','hemoglobin_biosynthetic_process','oxygen_transport','cortical_actin_cytoskeleton_organization')
gene_sets=c(cho_gene_list)
set.seed(123)
NRBC_altas_seu <- AddModuleScore_UCell(NRBC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快

# Mitophagy：GO:0000423
# Regulation of mitophagy ：GO:1901524
# Positive regulation of mitophagy：GO:1901526
mitophagy_go <- c( "GO:0000423",  "GO:1901524", "GO:1901526" )
mitophagy_gene_list=mapIds(x = org.Hs.eg.db,keys =mitophagy_go,column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
names(mitophagy_gene_list)=c('Mitophagy','Regulation of mitophagy','Positive regulation of mitophagy')
mitophagy_gene_list=mitophagy_gene_list[1]

NRBC_altas_seu <- AddModuleScore_UCell(NRBC_altas_seu, features = mitophagy_gene_list,ncores = 2) # 计算速度很快
colnames(NRBC_altas_seu@meta.data) =gsub(pattern = 'GO:0000423_UCell',replacement ='Mitophagy_UCell' ,colnames(NRBC_altas_seu@meta.data))


VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features = c('Mitophagy_UCell'),pt.size = 0)+NoLegend()
NRBC_altas_seu$final_celltype=factor(NRBC_altas_seu$final_celltype,levels =c( "BFUE/CFUE","EryP-CFC","ProE","Bas","Poly","Orth") )

cho_pathways=c('heme_biosynthetic_process_UCell' ,'hemoglobin_biosynthetic_process_UCell' ,'oxygen_transport_UCell','Mitophagy_UCell','HRI_mediated_signaling_UCell') # cortical_actin_cytoskeleton_organization 不明显
p=VlnPlot(NRBC_altas_seu,group.by = 'final_celltype',split.by = 'tissue_stage',stack = T,cols = col,flip = T,
        features =cho_pathways )
p
ggsave(p,filename='res_pic/main_figure2/pd_maturation_score_stage_vlnplot.pdf',width = 8,height = 8)

saveRDS(NRBC_altas_seu@meta.data,file = 'NRBC_altas_seu_metadata_UCell_score.rds')




test_p_list=list()
# 执行 t 检验
for( cho_pathway in cho_pathways ){
  print(cho_pathway)
  for( celltype in levels(NRBC_altas_seu$final_celltype)[4:6]){
    group1_expr <- NRBC_altas_seu@meta.data[NRBC_altas_seu@meta.data$NRBC_type == "definitive"& NRBC_altas_seu@meta.data$final_celltype == celltype, cho_pathway]
    group2_expr <- NRBC_altas_seu@meta.data[NRBC_altas_seu@meta.data$NRBC_type == "primitive" & NRBC_altas_seu@meta.data$final_celltype == celltype, cho_pathway]
    t_result <- t.test(group1_expr, group2_expr)
    test_p_list[[cho_pathway]][[celltype]] <- t_result$p.value
    
  }
}

for( cho_pathway in cho_pathways ){
  group1_expr <- NRBC_altas_seu@meta.data[NRBC_altas_seu@meta.data$NRBC_type == "definitive"& NRBC_altas_seu@meta.data$final_celltype == 'BFUE/CFUE', cho_pathway]
  group2_expr <- NRBC_altas_seu@meta.data[NRBC_altas_seu@meta.data$NRBC_type == "primitive" & NRBC_altas_seu@meta.data$final_celltype == 'EryP-CFC', cho_pathway]
  t_result <- t.test(group1_expr, group2_expr)
  test_p_list[[cho_pathway]][['progenitor']] <- t_result$p.value
}

#$oxygen_transport_UCell$Orth
# [1] 0.6818525
#$Mitophagy_UCell$progenitor
#[1] 0.8305787
# 其他 都是极端显著，接近0 

cho_gene_list=mapIds(x = org.Hs.eg.db,keys =c( 'GO:0030036','GO:0030866'),column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
names(cho_gene_list)=c('actin_cytoskeleton_organization','plasma_membrane_organization')
gene_sets=c(cho_gene_list)
set.seed(123)
NRBC_altas_seu <- AddModuleScore_UCell(NRBC_altas_seu, features = gene_sets,ncores = 6) # 计算速度很快
cho_pathways2=c('actin_cytoskeleton_organization_UCell' ,'plasma_membrane_organization_UCell' ) # cortical_actin_cytoskeleton_organization 不明显


p=VlnPlot(NRBC_altas_seu,group.by = 'final_celltype',split.by = 'tissue_stage',stack = T,cols = col,
          features =cho_pathways )
p


saveRDS(NRBC_altas_seu@meta.data,file = 'NRBC_altas_seu_metadata_UCell_score.rds')





# HSPC-derive NRBC 最显著特的特征就是免疫调控和末期的细胞复制
pathway=c('chromosome segregation','DNA-templated DNA replication','leukocyte proliferation','positive regulation of cell-cell adhesion','leukocyte cell-cell adhesion','B cell activation','positive regulation of T cell activation')
gene_list=unique(unlist(strsplit(pd_enrichgo_res@compareClusterResult[pd_enrichgo_res@compareClusterResult$Description %in% pathway & pd_enrichgo_res@compareClusterResult$Cluster=='definitive','geneID'],split = '/',fixed = T)))
gene_list_df=pd_whole_level_markers[pd_whole_level_markers$gene %in% gene_list,]

DoHeatmap(NRBC_altas_seu,features = gene_list_df[gene_list_df$avg_log2FC >1 & gene_list_df$pct.2 <0.1,'gene'])
pheatmap(celltype_mexp_df[gene_list_df[gene_list_df$avg_log2FC >1 & gene_list_df$pct.2 <0.1,'gene'],],cluster_cols = F,scale = 'row',color =colorRampPalette(colors = c('navy','white','#C31E1F'))(100))


gene_sets <- list( macroautophagy = unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='macroautophagy'],split = '/')),
                   'regulation_of_autophagy'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='regulation of autophagy'],split = '/')),
                   'vacuole organization'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='vacuole organization'],split = '/')),
                   'proteasome-mediated_ubiquitin-dependent_protein_catabolic_process'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='proteasome-mediated ubiquitin-dependent protein catabolic process'],split = '/')),
                   'DNA replication'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='DNA replication'],split = '/')),
                   'positive_regulation_of_T_cell_activation'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='positive regulation of T cell activation'],split = '/')),
                   'regulation_of_leukocyte_proliferation'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='regulation of leukocyte proliferation'],split = '/')),
                   'positive_regulation_of_leukocyte_cell-cell adhesion'=unlist(strsplit(top10_enrichgo_res_df$geneID[top10_enrichgo_res_df$Description=='positive regulation of leukocyte cell-cell adhesion'],split = '/')),
                   'stress response_to_metal_ion'= unlist(strsplit(top_pd_whole_level_marker_gseGO_res$core_enrichment[top_pd_whole_level_marker_gseGO_res$Description=='stress response to metal ion'],split = '/')),
                   'detoxification_of_copper_ion'= unlist(strsplit(top_pd_whole_level_marker_gseGO_res$core_enrichment[top_pd_whole_level_marker_gseGO_res$Description=='detoxification of copper ion'],split = '/')),
                   'intracellular_zinc_ion_homeostasis'= unlist(strsplit(top_pd_whole_level_marker_gseGO_res$core_enrichment[top_pd_whole_level_marker_gseGO_res$Description=='intracellular zinc ion homeostasis'],split = '/'))
)



########################################################################################################################################################################
# ---------------HSPC-derived NRBC unique signature :ABM BFUE/CFUE-----------------#
########################################################################################################################################################################

if(F){
  #TAPBP	TAP结合蛋白（亦称 Tapasin）	抗原肽装载的“质量控制员”与“桥梁”	专一服务于 MHC-I类分子。
  #它在内质网中物理性桥接抗原肽转运体 和新合成的MHC-I重链-β2M复合物，确保将合适长度和亲和力的抗原肽高效、精准地装载到MHC-I分子的肽结合槽中，是MHC-I成熟的关键步骤
  #B2M	β2-微球蛋白	MHC-I类分子的结构基石与稳定伴侣	是 所有MHC-I类分子（如HLA-A, B, C）的不变轻链。
  #它与MHC-I重链非共价结合，对其正确折叠、稳定结构和转运至细胞表面绝对必需。没有B2M，MHC-I重链无法稳定存在。它与MHC-II类分子无关。
  
  cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep('MHC|antigen',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
  MHC_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
  sort(MHC_genes)
  VlnPlot(NRBC_altas_seu ,features= MHC_genes,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()
  MHC_genes=c("B2M","HLA-A","HLA-B","HLA-C","HLA-E","CD74","HLA-DMA","HLA-DPA1","HLA-DPB1","HLA-DRA" ,"HLA-DRB1", "HLA-DRB5") # order, and cho
  VlnPlot(NRBC_altas_seu ,features= MHC_genes,stack = T,cols = cols,split.by = 'tissue_stage',group.by = 'source_celltype')+NoLegend()
  
  peptide_load=c('TAPBP','TAP1','TAP2','ERAP1')
  immune_Proteasome_genes=c(  "PSMB8", "PSMB9","PSMB10"  )
  VlnPlot(NRBC_altas_seu ,features =c(peptide_load,immune_Proteasome_genes),stack = T,group.by = 'source_celltype',cols = cols,split.by = 'tissue_stage')
  
  cho_enrich_pathways=top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description[grep('monocyte|leukocyte|T cell',top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description)]
  immune_remodle_genes=unique(unlist(strsplit(top_sub_pd_all_Ery_tissue_markers_gseGO_res[top_sub_pd_all_Ery_tissue_markers_gseGO_res$Description %in% cho_enrich_pathways,]$core_enrichment,split = '/')))
  sort(immune_remodle_genes)
  immune_remodle_genes=immune_remodle_genes[!immune_remodle_genes %in% MHC_genes ]
  VlnPlot(NRBC_altas_seu ,features =immune_remodle_genes,stack = T,group.by = 'source_celltype',cols = cols)+NoLegend()
  #转录调控,信号转导，细胞因子以及受体，细胞粘附
  cho_immune_remodle_genes=c( "SOX4",'FOXP1',"ZBTB16", "SMARCA2", "PTPRC","INPP5D", "RHOH","TESPA1","NCKAP1L", 'IL1B',"LGALS3","LGALS9","CSF2RB","IL2RG","ITGA4","CD44" ,'SPN') # PNP:嘌呤代谢酶，ATAD5：	DNA损伤应答/周期调控
  VlnPlot(NRBC_altas_seu ,features =cho_immune_remodle_genes,stack = T,group.by = 'source_celltype',cols = cols)+NoLegend()
  
  
  #SOX4	SRY-box转录因子4	免疫细胞发育的“命运调控器”。在B细胞和T细胞的早期发育、谱系决定与存活中不可或缺。它调控多种靶基因，影响淋巴细胞的增殖、分化和功能。
  #TESPA1:T细胞受体信号转导的“专用适配器”。在胸腺T细胞阳性选择阶段，特异性介导TCR信号向下游（如NFAT通路）的高效传递，对功能性T细胞库的形成至关重要。
  #CD74	恒定链（Ii链）	MHC-II类分子抗原呈递的“专职伴侣与导航员”。作为伴侣蛋白，防止MHC-II在内质网中过早结合抗原；作为靶向信号，引导MHC-II至内吞系统装载外源抗原肽。
  #DOCK8	胞质因子DOCK8	淋巴细胞迁移、存活与功能的“支架”。一种鸟嘌呤核苷酸交换因子（GEF），通过激活CDC42等小G蛋白，调控淋巴细胞（尤其CD8⁺ T和B细胞）的骨架重组、免疫突触形成、存活和效应功能。
  #INPP5D	SHIP1磷酸酶	髓系与B细胞的负向调控“刹车”。通过降解PIP3，负向调控PI3K-Akt-mTOR这一关键的促存活、增殖和代谢通路，防止免疫细胞过度活化，在免疫耐受和稳态中起关键作用。
  #LEF1	淋巴样增强结合因子1	T/B淋巴细胞发育的“命运决定者”。是Wnt/β-catenin信号通路下游的关键转录因子，对早期T细胞发育、CD4⁺ T细胞谱系决定以及B细胞生成不可或缺。
  #RHOH	Rho GTP酶H	T细胞发育与信号转导的“特化调控器”。一种缺乏GTP酶活性的特殊Rho蛋白，主要在T细胞中表达。在胸腺细胞发育的β选择阶段和成熟T细胞受体（TCR）信号传导中起关键调控作用。
  # SOX4,TESPA1, CD74
  # 显著性检验
  immune_pathway_related_genes=c('SOX4','TESPA1','CD74','DOCK8','INPP5D','LEF1','RHOH') #
  VlnPlot(NRBC_altas_seu ,features=c(MHC_genes,immune_pathway_related_genes) ,group.by = 'source_celltype',split.by = 'tissue_stage',cols = cols,pt.size = 0,raster=FALSE,alpha = 0.6,stack = T)+NoLegend()
  
}


