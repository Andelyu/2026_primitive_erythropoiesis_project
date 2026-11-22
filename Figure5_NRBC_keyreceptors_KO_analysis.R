
library(ggplot2)
library(patchwork)
library(ggplotify)
library(ggnetwork)
library(dplyr)
library(stringr)
library(limma)
library(ggplot2)
library(ggbreak) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggrepel)

library(RColorBrewer)
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
setwd('/home/gibh/2021_NRBC_chlyu/')
set.seed(1234)



cells=sample(rownames(filt_YS_altas_Ery_seu@meta.data),size = 5000,replace = F)
sub_filt_YS_altas_Ery_seu=subset(filt_YS_altas_Ery_seu,cells=cells)
saveRDS(sub_filt_YS_altas_Ery_seu,file = 'primitive_NRBC_project/sub_5k_filt_YS_altas_Ery_seu.rds')


draw_gonetcwork_pic_func=function(res=enrichGO_res,showCategory=20,xlimits = c(-0.2, 2.5)){
  
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


setwd('~/2021_NRBC_chlyu/keygenes_KOres')
dir.create('main_figure5')
dir.create('res_data')

library(ggrepel)



library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 提取 LogNormalized 数据 (千万不要用 raw counts)
df <- FetchData(NRBC_altas_seu, vars = c("TMEM219", "EPOR", "tissue_stage"), layer = "data") 

# 2. 过滤掉双 0 的细胞（关键！因为双 0 会严重拉低相关性）
df_filtered <- df %>% filter(TMEM219 > 0 & EPOR > 0)

# 3. 按细胞亚群计算 Spearman 相关性
cor_res <- df_filtered %>%
  group_by(tissue_stage) %>%
  summarise(
    Spearman_r = cor(TMEM219, EPOR, method = "spearman"),
    p_value = cor.test(TMEM219, EPOR, method = "spearman")$p.value,
    n_cells = n()
  )
print(cor_res)

# 4. 画图 (带有真正的拟合线和相关性系数)
p=ggplot(df_filtered, aes(x = TMEM219, y = EPOR, color = tissue_stage)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", color = "black", se = FALSE) + # 加一条黑色的总体拟合线
  facet_wrap(~tissue_stage) +  # 按细胞类型分面
  theme_bw() +
  labs(title = "TMEM219 vs EPOR (Filtered & LogNormalized)")

ggsave(p,filename='main_figure5/correlation_dotplot.pdf',width = 8,height = 8)



focused_gene_sets <- list(
  
  # 1. Cell cycle regulation
  Cell_Cycle_Progression = c(
    "CDK1", "CDK2", "CDK4", "CDK6",
    "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2",
    "CCNA2", "CCNB1", "CCNB2",
    "MYC", "MKI67", "PCNA", "TOP2A",
    "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7"
  ),
  
  Cell_Cycle_Inhibition = c(
    "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B",
    "BTG1", "BTG2", "GADD45A", "GADD45B",
    "RBL1", "RBL2"
  ),
  
  # 2. Integrated stress response and metabolic stress
  ISR_HRI_ATF4 = c(
    "EIF2AK1", "EIF2S1", "ATF4", "DDIT3", "ATF3",
    "GDF15", "PPP1R15A", "TRIB3", "ASNS", "CHAC1"
  ),
  
  Metabolic_Stress = c(
    "DDIT4", "SESN2", "SLC7A11",
    "SLC7A5", "EIF4EBP1", "HMOX1"
  ),
  
  # 3. Cell survival and apoptosis
  Anti_Apoptosis = c(
    "BCL2L1", "BCL2", "MCL1",
    "BCL2L2", "XIAP", "BIRC5"
  ),
  
  Pro_Apoptosis = c(
    "BAX", "BAK1", "BCL2L11", "BBC3", "PMAIP1",
    "CASP3", "CASP8", "CASP9", "APAF1"
  ),
  
  # 4. Erythroid functional maturation
  Globin_Production = c(
    "HBA1", "HBA2", "HBB", "HBE1",
    "HBZ", "HBG1", "HBG2", "AHSP"
  ),
  
  Heme_Biosynthesis = c(
    "ALAS2", "ALAD", "HMBS", "UROS",
    "UROD", "CPOX", "PPOX", "FECH"
  ),
  
  Iron_Homeostasis = c(
    "SLC25A37", "SLC25A28", "TFRC",
    "SLC11A2", "FTH1", "FTL", "SLC40A1"
  ),
  
  Erythrocyte_Membrane = c(
    "SLC4A1", "ANK1", "EPB42", "EPB41",
    "SPTA1", "SPTB", "GYPA", "GYPC", "RHAG"
  ),
  
  Antioxidant_Defense = c(
    "SOD1", "SOD2", "SOD3", "GPX1", "GPX4",
    "PRDX1", "PRDX2", "CAT", "G6PD",
    "TXN", "TXNRD1"
  ),
  
  Mitophagy = c(
    "BNIP3L", "BNIP3", "ULK1", "SQSTM1",
    "MAP1LC3B", "PINK1", "PRKN"
  ),
  
  Erythroid_Transcription_Factors = c(
    "GATA1", "KLF1", "NFE2",
    "SOX6", "TAL1", "GFI1B"
  )
  
)

####################################################### ##################################################################
# --------------------------figure5  IGFBP3-TMEM219 、EPOR ----------------------#
####################################################### ##################################################################

def_YS_YS_altas_Ery_seu=readRDS('../primitive_NRBC_project/YS_def_like_Ery_seu.rds') # 662  细胞
VlnPlot(def_YS_YS_altas_Ery_seu,features = c('TMEM219','EPOR'),group.by = 'Pre_celltype2')
DimPlot(def_YS_YS_altas_Ery_seu,reduction = 'umap',group.by = 'Phase',cols = cols)


####################################################### ##################################################################
#-----------------------------------虚拟扰动分析-----------------------------------#
####################################################### ##################################################################
plot_KOres_func=function(df,type='primitive ',KO_gene='TMEM219'){
  top_df=df[1:10,]
  p=ggplot(df, aes(y = -log10(p.value), x = log2(FC), color = log2(FC))) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_gradient(low = "gray", high = "firebrick3") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_df, aes(label = gene), size = 4, max.overlaps = 20, color = "black", fontface = "bold") +
    labs(title = paste(type,"NRBC KO ",KO_gene), x = "log2 Fold Change", y = "-log10(adjusted p-value)", color = "") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "top")
  print(p)
  return(p)
}


##  --------------------------------TMEM219 -------------------------#
#---------------------raw： 保留血红蛋白：HBB，HBA1，HBA2，HBE1，HBG1,HBG2-----------------------#

TEME219_knk_gene_res_list=list()
TEME219_knk_gene_res_list[['YS']]=readRDS('scTenifoldKnk/KO_data/YS_TMEM219_KO_res_raw.rds')
TEME219_knk_gene_res_list[['def_fetal']]=readRDS('scTenifoldKnk/KO_data/fetal_TMEM219_KO_res_raw.rds')
TEME219_knk_gene_res_list[['def_abm']]=readRDS('scTenifoldKnk/KO_data/ABM_TMEM219_res.rds')

# 太少了，说明血红蛋白干扰比较大
table(TEME219_knk_gene_res_list[['YS']]$diffRegulation$p.adj <0.05)# 2
table(TEME219_knk_gene_res_list[['def_fetal']]$diffRegulation$p.adj <0.05) # 6
table(TEME219_knk_gene_res_list[['def_abm']]$diffRegulation$p.adj <0.05) # 6 

TMEM219_diffRegulation_list=list('YS'=TEME219_knk_gene_res_list[['YS']]$diffRegulation,
                         'fetal'=TEME219_knk_gene_res_list[['def_fetal']]$diffRegulation,
                         'abm'=TEME219_knk_gene_res_list[['def_abm']]$diffRegulation
)


df <- TMEM219_diffRegulation_list[['YS']][-1,]
p=plot_KOres_func(df = df)
ggsave(p,width = 6,height = 6,filename='main_figure5/TMEM219_YS_KO_diffRegulation.pdf')

df <- TMEM219_diffRegulation_list[['fetal']][-1,]
p=plot_KOres_func(df = df,type = 'definitive fetal')
ggsave(p,width = 6,height = 6,filename='main_figure5/TMEM219_fetal_KO_diffRegulation.pdf')

df <- TMEM219_diffRegulation_list[['abm']][-1,]
p=plot_KOres_func(df = df,type = 'definitive adult')
ggsave(p,width = 6,height = 6,filename='main_figure5/TMEM219_adult_KO_diffRegulation.pdf')

#---------------------filt： 剔除血红蛋白：HBB，HBA1，HBA2，HBE1，HBG1,HBG2-----------------------#
TEME219_knk_gene_res_list=list()
TEME219_knk_gene_res_list[['YS']]=readRDS('scTenifoldKnk/KO_data/YS_TMEM219_KO_res_new.rds')
TEME219_knk_gene_res_list[['def_fetal']]=readRDS('scTenifoldKnk/KO_data/fetal_TMEM219_KO_res_filt.rds')
TEME219_knk_gene_res_list[['def_abm']]=readRDS('scTenifoldKnk/KO_data/adult_TMEM219_KO_res_filt.rds')


table(TEME219_knk_gene_res_list[['YS']]$diffRegulation$p.adj <0.05)# 47
table(TEME219_knk_gene_res_list[['def_fetal']]$diffRegulation$p.adj <0.05) # 6
table(TEME219_knk_gene_res_list[['def_abm']]$diffRegulation$p.adj <0.05) # 6

TMEM219_diffRegulation_list=list('YS'=TEME219_knk_gene_res_list[['YS']]$diffRegulation,
                         'fetal'=TEME219_knk_gene_res_list[['def_fetal']]$diffRegulation,
                         'abm'=TEME219_knk_gene_res_list[['def_abm']]$diffRegulation
)


df <- TMEM219_diffRegulation_list[['YS']][-1,]
p=plot_KOres_func(df = df)
ggsave(p,width = 6,height = 6,filename='main_figure5/TMEM219_YS_KO_diffRegulation.pdf')

df <- TMEM219_diffRegulation_list[['fetal']][-1,]
p=plot_KOres_func(df = df)
ggsave(p,width = 6,height = 6,filename='main_figure5/TMEM219_fetal_KO_diffRegulation.pdf')
 
df <- TMEM219_diffRegulation_list[['abm']][-1,]
p=plot_KOres_func(df = df)
ggsave(p,width = 6,height = 6,filename='main_figure5/TMEM219_adult_KO_diffRegulation.pdf')



Ys_TMEM9_enrichgo_res=enrichGO(OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP',gene = TEME219_knk_gene_res_list[['YS']]$diffRegulation[TEME219_knk_gene_res_list[['YS']]$diffRegulation$p.adj <0.05,'gene'][-1])  
Ys_TMEM9_enrichgo_res@result$type='TEME219_KO' 
Ys_TMEM9_enrichgo_res@result$ratio=round(as.numeric(data.frame(strsplit2(Ys_TMEM9_enrichgo_res@result$GeneRatio,split = '/'))[,1])/as.numeric(data.frame(strsplit2(Ys_TMEM9_enrichgo_res@result$GeneRatio,split = '/'))[,2]),3)

p=draw_gonetcwork_pic_func(res = Ys_TMEM9_enrichgo_res)
p
ggsave(p,width = 8,height = 6,filename='main_figure5/EPOR_YS_KO_enrichres.pdf')


##################################################################################################################################################    
#--------------------------------k562 细胞系扰动结果，真实实验结果数据------------------------------#
##################################################################################################################################################    
file_list=list.files('GPSA/primitive')
file_list=file_list[grep('csv',file_list)]

# 只取第 1 和第 2 类（最核心的红系祖细胞和红白血病
core_erythroid_cells <- c("HUDEP-2","CD34 positive cells", "HSPCs","K-562","HEL","MEL-JUSO"  )

# 主要来自 K-562 
#CRISPRi      ko   shRNA   siRNA 

#-----------------------KO------------------------------------------#
keygene_list=list()
keygene_ko_res_list=list()
for( id in file_list){
  new_id=gsub(pattern = '_gpsa_query_table.csv',replacement = '',id)
  keygene_list[[new_id]]=read.csv(file = paste0('GPSA/primitive/',id))
  keygene_list[[new_id]]= keygene_list[[new_id]][ keygene_list[[new_id]]$CellLineName %in%  'K-562',]
  keygene_list[[new_id]] <- subset( keygene_list[[new_id]], method %in% c("ko"))
  print( table(keygene_list[[new_id]]$CellLineName ))
  keygene_list[[new_id]]$type='up'
  keygene_list[[new_id]]$type[keygene_list[[new_id]]$logFC<0]='dn'
  temp_df=keygene_list[[new_id]][,c('pbgene','logFC')]
  temp_df=temp_df[!duplicated(temp_df$pbgene),];
  print(new_id);print(dim(temp_df))
  genelist= temp_df$logFC
  names(genelist)= temp_df$pbgene
  genelist=sort(genelist,decreasing = T)
  #keygene_ko_res_list[[new_id]][['gsea']]=gseGO(geneList = genelist,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP') # 基因太少，不适合
  keygene_ko_res_list[[new_id]][['go_up']]=enrichGO(gene =names(genelist)[genelist >0],OrgDb =  org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
  keygene_ko_res_list[[new_id]][['go_dn']]=enrichGO(gene =names(genelist)[genelist <0],OrgDb =  org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
  keygene_ko_res_list[[new_id]][['go_up']]@result$ratio= as.numeric(t(data.frame(strsplit(keygene_ko_res_list[[new_id]][['go_up']]@result$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit(keygene_ko_res_list[[new_id]][['go_up']]@result$GeneRatio,split = '/')))[,2])
  keygene_ko_res_list[[new_id]][['go_dn']]@result$ratio= as.numeric(t(data.frame(strsplit(keygene_ko_res_list[[new_id]][['go_dn']]@result$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit(keygene_ko_res_list[[new_id]][['go_dn']]@result$GeneRatio,split = '/')))[,2])
  keygene_ko_res_list[[new_id]][['go_up']]@result$restype='up'
  keygene_ko_res_list[[new_id]][['go_dn']]@result$restype='dn'
}

saveRDS(keygene_list,file = 'res_data/K562_keygene_ko_list.rds')
saveRDS(keygene_ko_res_list,file = 'res_data/K562_keygene_ko_res_list.rds')



unique(keygene_list[['EPOR']]$CellLineName)
unique(keygene_list[['TMEM219']]$CellLineName)

table(keygene_list[['EPOR']]$method) 
#29      57      85       4 
table(keygene_list[['TMEM219']]$method) 
#CRISPRi      ko   shRNA   siRNA 
#25      63      60       2 

length(intersect( keygene_ko_res_list[['EPOR']][['go_dn']]@result$Description[keygene_ko_res_list[['EPOR']][['go_dn']]@result$p.adjust <0.05],
                  keygene_ko_res_list[['TMEM219']][['go_dn']]@result$Description[keygene_ko_res_list[['TMEM219']][['go_dn']]@result$p.adjust <0.05]))
# 0
length(intersect( keygene_ko_res_list[['EPOR']][['go_up']]@result$Description[keygene_ko_res_list[['EPOR']][['go_up']]@result$p.adjust <0.05],
                  keygene_ko_res_list[['TMEM219']][['go_up']]@result$Description[keygene_ko_res_list[['TMEM219']][['go_up']]@result$p.adjust <0.05]))
# 5


ggplot(keygene_list[['TMEM219']][keygene_list[['TMEM219']]$CellLineName=='K-562' & keygene_list[['TMEM219']]$pbgene!='TMEM219', ],aes(x=logFC,y=-log10(adj.P.Val)))+geom_point()+theme_bw()

top10=keygene_list[['EPOR']][keygene_list[['EPOR']]$CellLineName=='K-562' & keygene_list[['EPOR']]$pbgene!='EPOR', ]
top10=top10[order(top10$adj.P.Val),][1:10,]
ggplot(keygene_list[['EPOR']][keygene_list[['EPOR']]$CellLineName=='K-562' & keygene_list[['EPOR']]$pbgene!='EPOR', ],aes(x=logFC,y=-log10(adj.P.Val)))+geom_point()+theme_bw()+
  geom_text_repel(data = top10, aes(label = pbgene), size = 4, max.overlaps = 20, color = "black", fontface = "bold") +
  scale_color_gradient(low = 'gray',high = 'firebrick3')+ggtitle(label = paste(new_id,'KO genelist'))

draw_top_marker_topgores_pic=function(new_id,keygene_list=keygene_list,keygene_ko_res_list=keygene_ko_res_list,topn=5,topg=10){
  
  top_markers=keygene_list[[new_id]] %>% group_by(type) %>% top_n(n = topn,wt = -log10(adj.P.Val+10e-10))
  top_markers
  p1=ggplot(keygene_list[[new_id]],aes(x=logFC,y=-log10(adj.P.Val),color=-log10(adj.P.Val+10e-10)))+geom_point()+theme_classic()+
    geom_vline(xintercept = c(1,-1),linetype = "dashed", color = "red", linewidth = 0.4)+
    geom_text(data = top_markers,aes(x=logFC,y=-log10(adj.P.Val),label=pbgene),color='black',vjust = -0.5)+
    scale_color_gradient(low = 'gray',high = 'firebrick3')+ggtitle(label = paste(new_id,'KO genelist'))
  
  up_df=keygene_ko_res_list[[new_id]][['go_up']]@result[keygene_ko_res_list[[new_id]][['go_up']]@result$p.adjust <0.05,]
  dn_df=keygene_ko_res_list[[new_id]][['go_dn']]@result[keygene_ko_res_list[[new_id]][['go_dn']]@result$p.adjust <0.05,]
  
  temp_df=rbind(keygene_ko_res_list[[new_id]][['go_up']]@result[1:topg,],keygene_ko_res_list[[new_id]][['go_dn']]@result[1:topg,])
  temp_df$Description=factor(temp_df$Description,unique(temp_df$Description))
  p2=ggplot(temp_df,aes(x=restype,y=Description ,color= -log10(p.adjust),size= ratio))+geom_point()+
    theme_classic()+scale_color_gradient(low = 'gray',high = 'firebrick3')+ scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+theme(text = element_text(face = 'bold'))
  
  
  p=p1+p2+plot_layout(ncol = 2,widths = c(1,0.6))
  return(p)
}



p_EPOR=draw_top_marker_topgores_pic(new_id = 'EPOR',keygene_list =keygene_list,keygene_ko_res_list=keygene_ko_res_list ,topg=15)
p_TMEM219=draw_top_marker_topgores_pic(new_id = 'TMEM219',keygene_list =keygene_list ,keygene_ko_res_list=keygene_ko_res_list,topg=15)

ggsave(p_EPOR,width = 12,height = 8,filename='main_figure5/EPOR_K562_KO_enrichres.pdf')
ggsave(p_TMEM219,width = 12,height = 8,filename='main_figure5/TMEM219_K562_KO_enrichres.pdf')

keygene_list[['TMEM219']][ keygene_list[['TMEM219']]$pbgene %in% as.character(unlist(focused_gene_sets)),]
keygene_list[['EPOR']][ keygene_list[['EPOR']]$pbgene %in% as.character(unlist(focused_gene_sets)),]




#-----------------------KD------------------------------------------#
kn_keygene_list=list()
keygene_kn_res_list=list()
for( id in file_list){
  new_id=gsub(pattern = '_gpsa_query_table.csv',replacement = '',id)
  kn_keygene_list[[new_id]]=read.csv(file = paste0('GPSA/primitive/',id))
  kn_keygene_list[[new_id]]= kn_keygene_list[[new_id]][ kn_keygene_list[[new_id]]$CellLineName %in%  'K-562',]
  kn_keygene_list[[new_id]] <- subset( kn_keygene_list[[new_id]], method %in% c("shRNA"))
  print( table(kn_keygene_list[[new_id]]$CellLineName ))
  kn_keygene_list[[new_id]]$type='up'
  kn_keygene_list[[new_id]]$type[kn_keygene_list[[new_id]]$logFC<0]='dn'
  temp_df=kn_keygene_list[[new_id]][,c('pbgene','logFC')]
  temp_df=temp_df[!duplicated(temp_df$pbgene),];
  print(new_id);print(dim(temp_df))
  genelist= temp_df$logFC
  names(genelist)= temp_df$pbgene
  genelist=sort(genelist,decreasing = T)
  keygene_kn_res_list[[new_id]][['gsea']]=gseGO(geneList = genelist,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP')
  keygene_kn_res_list[[new_id]][['go_up']]=enrichGO(gene =names(genelist)[genelist >0],OrgDb =  org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
  keygene_kn_res_list[[new_id]][['go_dn']]=enrichGO(gene =names(genelist)[genelist <0],OrgDb =  org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
  keygene_kn_res_list[[new_id]][['go_up']]@result$ratio= as.numeric(t(data.frame(strsplit(keygene_kn_res_list[[new_id]][['go_up']]@result$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit(keygene_kn_res_list[[new_id]][['go_up']]@result$GeneRatio,split = '/')))[,2])
  keygene_kn_res_list[[new_id]][['go_dn']]@result$ratio= as.numeric(t(data.frame(strsplit(keygene_kn_res_list[[new_id]][['go_dn']]@result$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit(keygene_kn_res_list[[new_id]][['go_dn']]@result$GeneRatio,split = '/')))[,2])
  keygene_kn_res_list[[new_id]][['go_up']]@result$restype='up'
  keygene_kn_res_list[[new_id]][['go_dn']]@result$restype='dn'
}

saveRDS(kn_keygene_list,file = 'res_data/K562_kn_keygene_list.rds')
saveRDS(keygene_kn_res_list,file = 'res_data/K562_keygene_kn_res_list.rds')




length(intersect( keygene_kn_res_list[['EPOR']][['go_dn']]@result$Description[keygene_kn_res_list[['EPOR']][['go_dn']]@result$p.adjust <0.05],
                  keygene_kn_res_list[['TMEM219']][['go_dn']]@result$Description[keygene_kn_res_list[['TMEM219']][['go_dn']]@result$p.adjust <0.05]))
# 14
length(intersect( keygene_kn_res_list[['EPOR']][['go_up']]@result$Description[keygene_kn_res_list[['EPOR']][['go_up']]@result$p.adjust <0.05],
                  keygene_kn_res_list[['TMEM219']][['go_up']]@result$Description[keygene_kn_res_list[['TMEM219']][['go_up']]@result$p.adjust <0.05]))
# 79

cho_gene='TMEM219'
top10=kn_keygene_list[[cho_gene]][kn_keygene_list[[cho_gene]]$CellLineName=='K-562' & kn_keygene_list[[cho_gene]]$pbgene!=cho_gene, ]
top10=top10[order(top10$adj.P.Val),][1:10,]
ggplot(kn_keygene_list[[cho_gene]][kn_keygene_list[[cho_gene]]$CellLineName=='K-562' & kn_keygene_list[[cho_gene]]$pbgene!=cho_gene, ],aes(x=logFC,y=-log10(adj.P.Val),color=-log10(adj.P.Val)) )+geom_point()+theme_bw()+
  geom_text_repel(data = top10, aes(label = pbgene), size = 4, max.overlaps = 20, color = "black", fontface = "bold") +
  scale_color_gradient(low = 'gray',high = 'firebrick3')+ggtitle(label = paste(new_id,'KO genelist'))


cho_gene='EPOR'
top10=kn_keygene_list[[cho_gene]][kn_keygene_list[[cho_gene]]$CellLineName=='K-562' & kn_keygene_list[[cho_gene]]$pbgene!=cho_gene, ]
top10=top10[order(top10$adj.P.Val),][1:10,]
ggplot(kn_keygene_list[[cho_gene]][kn_keygene_list[[cho_gene]]$CellLineName=='K-562' & kn_keygene_list[[cho_gene]]$pbgene!=cho_gene, ],aes(x=logFC,y=-log10(adj.P.Val),color=-log10(adj.P.Val) ) )+geom_point()+theme_bw()+
  geom_text_repel(data = top10, aes(label = pbgene), size = 4, max.overlaps = 20, color = "black", fontface = "bold") +
  scale_color_gradient(low = 'gray',high = 'firebrick3')+ggtitle(label = paste(new_id,'KO genelist'))


p_EPOR=draw_top_marker_topgores_pic(new_id = 'EPOR',keygene_list =kn_keygene_list,keygene_ko_res_list =keygene_kn_res_list ,topg=15)
p_TMEM219=draw_top_marker_topgores_pic(new_id = 'TMEM219',keygene_list =kn_keygene_list ,keygene_ko_res_list =keygene_kn_res_list,topg=15)

ggsave(p_EPOR,width = 12,height = 8,filename='main_figure5/EPOR_K562_KD_enrichres.pdf')
ggsave(p_TMEM219,width = 12,height = 8,filename='main_figure5/TMEM219_K562_KD_enrichres.pdf')

kn_keygene_list[['TMEM219']][ kn_keygene_list[['TMEM219']]$pbgene %in% as.character(unlist(focused_gene_sets)),]
kn_keygene_list[['EPOR']][ kn_keygene_list[['EPOR']]$pbgene %in% as.character(unlist(focused_gene_sets)),]



#-----------------------CRISPRi------------------------------------------#
CRISPRi_keygene_list=list()
keygene_CRISPRi_res_list=list()
for( id in file_list){
  new_id=gsub(pattern = '_gpsa_query_table.csv',replacement = '',id)
  CRISPRi_keygene_list[[new_id]]=read.csv(file = paste0('GPSA/primitive/',id))
  CRISPRi_keygene_list[[new_id]]= CRISPRi_keygene_list[[new_id]][ CRISPRi_keygene_list[[new_id]]$CellLineName %in%  'K-562',]
  CRISPRi_keygene_list[[new_id]] <- subset( CRISPRi_keygene_list[[new_id]], method %in% c("CRISPRi"))
  print( table(CRISPRi_keygene_list[[new_id]]$CellLineName ))
  CRISPRi_keygene_list[[new_id]]$type='up'
  CRISPRi_keygene_list[[new_id]]$type[CRISPRi_keygene_list[[new_id]]$logFC<0]='dn'
  temp_df=CRISPRi_keygene_list[[new_id]][,c('pbgene','logFC')]
  temp_df=temp_df[!duplicated(temp_df$pbgene),];
  print(new_id);print(dim(temp_df))
  genelist= temp_df$logFC
  names(genelist)= temp_df$pbgene
  genelist=sort(genelist,decreasing = T)
  #keygene_CRISPRi_res_list[[new_id]][['gsea']]=gseGO(geneList = genelist,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP')
  keygene_CRISPRi_res_list[[new_id]][['go_up']]=enrichGO(gene =names(genelist)[genelist >0],OrgDb =  org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
  keygene_CRISPRi_res_list[[new_id]][['go_dn']]=enrichGO(gene =names(genelist)[genelist <0],OrgDb =  org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
  keygene_CRISPRi_res_list[[new_id]][['go_up']]@result$ratio= as.numeric(t(data.frame(strsplit(keygene_CRISPRi_res_list[[new_id]][['go_up']]@result$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit(keygene_CRISPRi_res_list[[new_id]][['go_up']]@result$GeneRatio,split = '/')))[,2])
  keygene_CRISPRi_res_list[[new_id]][['go_dn']]@result$ratio= as.numeric(t(data.frame(strsplit(keygene_CRISPRi_res_list[[new_id]][['go_dn']]@result$GeneRatio,split = '/')))[,1])/as.numeric(t(data.frame(strsplit(keygene_CRISPRi_res_list[[new_id]][['go_dn']]@result$GeneRatio,split = '/')))[,2])
  keygene_CRISPRi_res_list[[new_id]][['go_up']]@result$restype='up'
  keygene_CRISPRi_res_list[[new_id]][['go_dn']]@result$restype='dn'
}

saveRDS(CRISPRi_keygene_list,file = 'res_data/K562_CRISPRi_keygene_list.rds')
saveRDS(keygene_CRISPRi_res_list,file = 'res_data/K562_keygene_CRISPRi_res_list.rds')



length(intersect( keygene_CRISPRi_res_list[['EPOR']][['go_dn']]@result$Description[keygene_CRISPRi_res_list[['EPOR']][['go_dn']]@result$p.adjust <0.05],
                  keygene_CRISPRi_res_list[['TMEM219']][['go_dn']]@result$Description[keygene_CRISPRi_res_list[['TMEM219']][['go_dn']]@result$p.adjust <0.05]))
# 4
length(intersect( keygene_CRISPRi_res_list[['EPOR']][['go_up']]@result$Description[keygene_CRISPRi_res_list[['EPOR']][['go_up']]@result$p.adjust <0.05],
                  keygene_CRISPRi_res_list[['TMEM219']][['go_up']]@result$Description[keygene_CRISPRi_res_list[['TMEM219']][['go_up']]@result$p.adjust <0.05]))
# 14

cho_gene='TMEM219'
top10=CRISPRi_keygene_list[[cho_gene]][CRISPRi_keygene_list[[cho_gene]]$CellLineName=='K-562' & CRISPRi_keygene_list[[cho_gene]]$pbgene!=cho_gene, ]
top10=top10[order(top10$adj.P.Val),][1:10,]
ggplot(CRISPRi_keygene_list[[cho_gene]][CRISPRi_keygene_list[[cho_gene]]$CellLineName=='K-562' & CRISPRi_keygene_list[[cho_gene]]$pbgene!=cho_gene, ],aes(x=logFC,y=-log10(adj.P.Val),color=-log10(adj.P.Val)) )+geom_point()+theme_bw()+
  geom_text_repel(data = top10, aes(label = pbgene), size = 4, max.overlaps = 20, color = "black", fontface = "bold") +
  scale_color_gradient(low = 'gray',high = 'firebrick3')+ggtitle(label = paste(new_id,'KO genelist'))


cho_gene='EPOR'
top10=CRISPRi_keygene_list[[cho_gene]][CRISPRi_keygene_list[[cho_gene]]$CellLineName=='K-562' & CRISPRi_keygene_list[[cho_gene]]$pbgene!=cho_gene, ]
top10=top10[order(top10$adj.P.Val),][1:10,]
ggplot(CRISPRi_keygene_list[[cho_gene]][CRISPRi_keygene_list[[cho_gene]]$CellLineName=='K-562' & CRISPRi_keygene_list[[cho_gene]]$pbgene!=cho_gene, ],aes(x=logFC,y=-log10(adj.P.Val),color=-log10(adj.P.Val) ) )+geom_point()+theme_bw()+
  geom_text_repel(data = top10, aes(label = pbgene), size = 4, max.overlaps = 20, color = "black", fontface = "bold") +
  scale_color_gradient(low = 'gray',high = 'firebrick3')+ggtitle(label = paste(new_id,'KO genelist'))


p_EPOR=draw_top_marker_topgores_pic(new_id = 'EPOR',keygene_list =CRISPRi_keygene_list,keygene_ko_res_list =keygene_CRISPRi_res_list ,topg=15)
p_TMEM219=draw_top_marker_topgores_pic(new_id = 'TMEM219',keygene_list =CRISPRi_keygene_list ,keygene_ko_res_list =keygene_CRISPRi_res_list,topg=15)

ggsave(p_EPOR,width = 12,height = 8,filename='main_figure5/EPOR_K562_CRISPRi_enrichres.pdf')
ggsave(p_TMEM219,width = 12,height = 10,filename='main_figure5/TMEM219_K562_CRISPRi_enrichres.pdf')

CRISPRi_keygene_list[['TMEM219']][ CRISPRi_keygene_list[['TMEM219']]$pbgene %in% as.character(unlist(focused_gene_sets)),]
CRISPRi_keygene_list[['EPOR']][ CRISPRi_keygene_list[['EPOR']]$pbgene %in% as.character(unlist(focused_gene_sets)),]




#---------------------ISR 与TMEM219 的相关性分析--------------------------------#
NRBC_altas_seu=readRDS('~/2021_NRBC_chlyu/primitive_NRBC_project/all_NRBC_altas_seu.rds')
VlnPlot(NRBC_altas_seu,group.by = 'source_celltype',features = c('GPX1','CAT'),stack = T)+NoLegend()

# 计算EPOR 与TMEM219 的相关性， 共表达网路聚类分析
p1=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & ATF4>0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "ATF4", cols = cols)
p2=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & EIF2AK1>0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "EIF2AK1", cols = cols)
p3=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & GPX1>0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "GPX1", cols = cols)
p=p1+p2+p3+plot_layout(ncol = 1,heights = c(1,1,1))

cor.test()
ggsave(p,width = 5,height = 12,filename='main_figure5/YS_TMEM219_ATF4_EIF2AK1_GPX1_correlation.pdf')

p1=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & EIF2S1>0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "EIF2S1", cols = cols)
p2=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & DDIT3>0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "DDIT3", cols = cols)
p3=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & CAT>0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "CAT", cols = cols)
p=p1+p2+p3+plot_layout(ncol = 1,heights = c(1,1,1))
p
ggsave(p,width = 5,height = 12,filename='main_figure5/YS_TMEM219_EIF2S1_DDIT3_CAT_correlation.pdf')



p1=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & SLC40A1 >0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "SLC40A1", cols = cols) # 细胞铁输出的主要转运蛋白
p2=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & TFRC >0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "TFRC", cols = cols) # 铁摄取
p=p2+p1+plot_layout(ncol = 2);p
ggsave(p,width = 10,height = 4,filename='main_figure5/YS_TMEM219_SLC40A1_TFRC_correlation.pdf')

FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & STEAP3 >0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "STEAP3", cols = cols) # STEAP3


p=FeatureScatter(subset(subset(NRBC_altas_seu, TMEM219 >0 & EPOR>0),tissue_stage=='YS'), feature1 = "TMEM219", feature2 = "EPOR", cols = cols)
ggsave(p,width = 5,height = 12,filename='main_figure5/YS_TMEM219_EPOR_correlation.pdf')

temp_df=GetAssayData(object =subset(NRBC_altas_seu, TMEM219 >0 & EPOR>0),layer='data' )
temp_df=data.frame(t(temp_df[c('TMEM219','EPOR'),]))
temp_df$tissue_stage=NRBC_altas_seu@meta.data[rownames(temp_df),'tissue_stage']
temp_df$celltype=NRBC_altas_seu@meta.data[rownames(temp_df),'final_celltype']

p=ggplot(temp_df,aes(x=TMEM219,y=EPOR))+geom_point(aes(color=celltype))+facet_grid(~tissue_stage)+scale_color_manual(values = cols)+
  geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
  stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
p

ggsave(p,filename='main_figure5/YS_TMEM219_EPOR_correlation.pdf',height = 5,width = 12)





#功能位置	代表基因
#ISR上游激酶	EIF2AK1、EIF2AK2、EIF2AK3、EIF2AK4
#核心信号	EIF2S1
#核心转录因子	ATF4
#ISR下游靶点	ATF3、DDIT3、ASNS、CHAC1、TRIB3
#负反馈调控	PPP1R15A、PPP1R15B

ISR_genes <- c("ATF4",'EIF2S1', "DDIT3", "ATF3", "ASNS", "CHAC1", "TRIB3", "GDF15",  "PPP1R15A",'PPP1R15B', "DDIT4", "SESN2",'EIF2AK1','EIF2AK2','EIF2AK3','EIF2AK4')
Redox_genes <- c( "HMOX1",'SOD1', "SOD2", "GPX1", "CAT", "NQO1", "GCLC", "GCLM")
Iron_genes <- c( "FTH1", "FTL", "SLC40A1","TFRC") # Fe


aggre_df=AggregateExpression(subset(NRBC_altas_seu,tissue_stage =='YS'),group.by = 'id')$RNA
aggre_df=log2(cpm(aggre_df)+1)
cho_aggre_df=aggre_df[c('TMEM219',ISR_genes,Redox_genes,Iron_genes),]
p=pheatmap(cor((t(cho_aggre_df))),annotation_row = an_df ,color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))
ggsave(as.ggplot(p),width = 10,height = 8,filename='main_figure5/YS_TMEM219_cor_ISR_Redox_metal_ion.pdf')

p=pheatmap(cor((t(cho_aggre_df)),method = 'spearman'),annotation_row = an_df , color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))
ggsave(as.ggplot(p),width = 10,height = 8,filename='main_figure5/YS_TMEM219_spearman_ISR_Redox_metal_ion.pdf')

# 
UCell_score_metad_df=readRDS('../primitive_NRBC_project/NBRC_altas_seu_metadata_UCell_score.rds')
data=GetAssayData(NRBC_altas_seu,layer = 'data')
UCell_score_metad_df$TMEM219=data['TMEM219',rownames(UCell_score_metad_df)]
rm(data);gc()

p=ggplot(UCell_score_metad_df[UCell_score_metad_df$tissue_stage=='YS',],aes(x=HRI_mediated_signaling_UCell,y=TMEM219))+#geom_point()+
  geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
  stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
p

ggsave(p,filename='main_figure6/AML_TMEM219_ISR_cor.pdf',width = 5,height = 5)



######################################################################################################################################## 
################################### ISR 与TMEM219 的在癌症中的拓展分析##########################################
######################################################################################################################################## 

# Figure 6
dir.create('main_figure6')
# M6 样本类型太少，不适用，2 个， 确实表现更高EPOR，IGFBP3 表达倾向
{
  library(edgeR)
  AML_count_df=read.csv('../primitive_NRBC_project/raw_data/AML_data/LAML.symbol_count.csv',sep="\t",header = T)
  AML_ph_df=read.csv('../primitive_NRBC_project/raw_data/AML_data/TCGA-LAML.GDC_phenotype.tsv.gz',sep="\t",header = T)
  AML_cl_df=read.csv('../primitive_NRBC_project/raw_data/AML_data/TCGA-LAML.survival.tsv.gz',sep="\t",header = T)
  
  rownames(AML_count_df)=AML_count_df$Symbol;AML_count_df=AML_count_df[,-1]
  AML_count_df= AML_count_df[rowMeans(AML_count_df)>1,]
  dim(AML_count_df)
  
  length(unique(rownames(AML_count_df)));dim(AML_count_df)
  rownames(AML_count_df)=AML_count_df$Symbol;AML_count_df=AML_count_df[,-1]
  rownames(AML_ph_df)=AML_ph_df$submitter_id.samples
  
  table(rownames(AML_ph_df) %in% colnames(AML_count_df))
  head( colnames(AML_count_df))
  head(rownames(AML_ph_df))
  rownames(AML_ph_df)=gsub(pattern = '-',rownames(AML_ph_df),replacement = '.')
  rownames(AML_cl_df)=gsub(pattern = '-',rownames(AML_cl_df),replacement = '.')
  
  
  table(rownames(AML_ph_df) %in% colnames(AML_count_df))
  AML_ph_df=AML_ph_df[ colnames(AML_count_df),]
  AML_cl_df=AML_cl_df[ colnames(AML_count_df),]
  
  
  
  AML_count_nr_df=2**AML_count_df-1
  summary(as.numeric(as.matrix(AML_count_df)))
  head(AML_count_df[rowMeans(AML_count_df) >25,1:4])
  
  AML_count_nr_df=AML_count_nr_df[! rownames(AML_count_nr_df) %in% rownames(AML_count_df[rowMeans(AML_count_df) >25,1:4]),]
  # 去掉5S_sRNA
  AML_count_nr_df=log2(cpm(AML_count_nr_df)+1) # 数据有些问题
  head(AML_count_nr_df[,1:4])
  
  cor.test(t(AML_count_nr_df)[,c('TMEM219')],t(AML_count_nr_df)[,c('EPOR')]) #
  

  library(GSVA)
  library(edgeR)
  
  

  gsvapar <- gsvaParam(as.matrix(AML_count_nr_df), list(ISR_genes=ISR_genes))
  ISR_score <- gsva(gsvapar)
  
  gsvapar <- gsvaParam(as.matrix(AML_count_nr_df), list(Redox_genes=Redox_genes))
  Redox_score <- gsva(gsvapar) # Redox 是氧化还原（Reduction–Oxidation）
  
  gsvapar <- gsvaParam(as.matrix(AML_count_nr_df), list(Iron_genes=Iron_genes))
  Iron_score <- gsva(gsvapar)
  
  LAML_df=data.frame(t(AML_count_nr_df[c('TMEM219',unique(c(ISR_genes,Redox_genes,Iron_genes))),]))
  LAML_df$ISR_score=as.numeric(ISR_score)
  LAML_df$Redox_score=as.numeric(Redox_score)
  LAML_df$Iron_score=as.numeric(Iron_score)
 
  cor.test(AML_count_nr_df['TMEM219',],AML_count_nr_df['EPOR',])
  p=ggplot(t(AML_count_nr_df),aes(x=EPOR,y=TMEM219))+geom_point()+
    geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
    stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
  p
  
  cor.test( LAML_df$TMEM219, LAML_df$ISR_score )
  cor.test( LAML_df$TMEM219, LAML_df$Redox_score )
  cor.test( LAML_df$TMEM219, LAML_df$Iron_score )
  
  an_df=data.frame(row.names = c('TMEM219',ISR_genes,Redox_genes,Iron_genes),type=c('focused',rep('ISR',length(ISR_genes)),
                                                                                  rep('Redox',length(Redox_genes)),
                                                                                  rep('Fe_iron',length(Iron_genes))))
  pheatmap(cor(LAML_df[,c('TMEM219','ISR_score','Redox_score','Iron_score')]),display_numbers = T)
  
  p=pheatmap(cor(LAML_df[,c('TMEM219',ISR_genes,Redox_genes,Iron_genes)]),display_numbers = F,annotation_row = an_df,
           color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))
  ggsave(as.ggplot(p),width = 10,height = 8,filename='main_figure5/AML_TMEM219_cor_ISR_Redox_metal_ion.pdf')

  
  p=ggplot(LAML_df,aes(x=ATF4,y=TMEM219))+geom_point()+
    geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
    stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
  p
  ggsave(p,filename='main_figure6/AML_TMEM219_ATF_cor.pdf',width = 5,height = 5)
  
  p=ggplot(LAML_df,aes(x=EIF2S1,y=TMEM219))+geom_point()+
    geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
    stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
  p
  ggsave(p,filename='main_figure6/AML_TMEM219_EIF2S1_cor.pdf',width = 5,height = 5)
  
  
  pval=cor.test(LAML_df$GPX1,LAML_df$TMEM219,method = 'spearman')
  pval=cor.test(LAML_df$GPX1,LAML_df$TMEM219)
  p=ggplot(LAML_df,aes(x=GPX1,y=TMEM219))+geom_point()+
    geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
    stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4)
  p
  ggsave(p,filename='main_figure6/AML_TMEM219_GPX1_cor.pdf',width = 5,height = 5)
  
  

  p=ggplot(LAML_df,aes(x=ISR_score,y=TMEM219))+geom_point()+
    geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
    stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
  
  p
  
  ggsave(p,filename='main_figure6/AML_TMEM219_ISR_cor.pdf',width = 5,height = 5)
  
  p=ggplot(LAML_df,aes(x=Redox_score,y=TMEM219))+geom_point()+
    geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
    stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
  p
  
  ggsave(p,filename='main_figure6/AML_TMEM219_Redox_score_cor.pdf',width = 5,height = 5)
  
  
  p=ggplot(LAML_df,aes(x=Iron_score,y=Redox_score))+geom_point()+
    geom_smooth(method = "lm", se = T,linewidth = 1,col='firebrick3')+theme_bw()+
    stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4) 
  p
  
  ggsave(p,filename='main_figure6/AML_ISR_Redox_score_Redox_score_cor.pdf',width = 5,height = 5)
  
   rownames(AML_cl_df)=gsub(pattern = '-',replacement = '.',x = AML_cl_df$sample)
  share_cli_sample_ids<- intersect(colnames(AML_count_nr_df),rownames(AML_cl_df));length(share_cli_sample_ids)
  cli_count_inf_df <- cbind(t(AML_count_nr_df)[share_cli_sample_ids,],AML_cl_df[share_cli_sample_ids,c('OS','OS.time')])
  cli_count_inf_df$cancer='LAML'
  cli_count_inf_df$OS.time=cli_count_inf_df$OS.time/30
  TMEM219_res_pic=single_gene_os_cox_func(geneid = 'TMEM219',data_df =cli_count_inf_df,pic = T,cancer='LAML' )
  TMEM219_res=single_gene_os_cox_func(geneid = 'TMEM219',data_df =cli_count_inf_df,pic = F,cancer='LAML' )
  
  ATF4_res_pic=single_gene_os_cox_func(geneid = 'ATF4',data_df =cli_count_inf_df,pic = T,cancer='LAML' )
  ATF4_res=single_gene_os_cox_func(geneid = 'ATF4',data_df =cli_count_inf_df,pic = F,cancer='LAML' )
  
  
  single_gene_os_cox_func<-function(geneid,data_df,pic=F,os_type='OS',cancer){
    # COX analysis 
    #print(data_df[1:5,])
    os_type_time=paste0(os_type,'.time')
    data_df<- data_df[!is.na(data_df[,os_type_time]),]
    #temp_formula<- as.formula(paste0('Surv(OS.time,OS)~',geneid))
    temp_formula<- as.formula(paste0('Surv(',os_type_time,',',os_type,')~',geneid))
    surv_scox_res <- summary(coxph(temp_formula,data =data_df))
    surv_scox_ph_p <- cox.zph(coxph(temp_formula,data =data_df ))$table[1,3]
    
    # KM survival analysis
    #data_df$q3 <- 'low';data_df$q4 <- 'low'
    data_df$q3 <- 'low';data_df$q4 <- 'low'
    data_df$q3[data_df[,geneid] >= quantile(data_df[,geneid])[3]]='high'
    data_df$q4[data_df[,geneid] >= quantile(data_df[,geneid])[4]]='high'
    
    #https://cloud.tencent.com/developer/article/1725994
    # surv_fit是survfit的扩展，更有利于解析Surv(OS.time,OS)公式
    
    #q3_km_p <-surv_pvalue( survfit(Surv(OS.time,OS)~q3,data = data_df),data = data_df)$pval
    #q4_km_p <-surv_pvalue( survfit(Surv(OS.time,OS)~q4,data = data_df),data = data_df)$pval
    
    temp_formula1=as.formula(paste0('Surv(',os_type_time,',',os_type,')~q3'))
    q3_km_p <-surv_pvalue(surv_fit(temp_formula1,data = data_df),data = data_df)$pval
    temp_formula2=as.formula(paste0('Surv(',os_type_time,',',os_type,')~q4'))
    q4_km_p <-surv_pvalue( surv_fit(temp_formula2,data = data_df),data = data_df)$pval
    km_pvalue <- min(q3_km_p,q4_km_p,na.rm = T) 
    
    # ROC survival analysis
    
    roc_res1= survivalROC(Stime=data_df[,os_type_time], status=data_df[,os_type], marker = data_df[,geneid], predict.time = 60,lambda = 0.1 )
    roc_res2= survivalROC(Stime=data_df[,os_type_time], status=data_df[,os_type], marker = data_df[,geneid], predict.time = 60,method = 'KM' )
    roc_res3= survivalROC(Stime=data_df[,os_type_time], status=data_df[,os_type], marker = data_df[,c('q3','q4')[which.min(c(q3_km_p,q4_km_p))]], predict.time = max(data_df[,os_type_time]),method = 'KM' )
    
    singe_gene_cox_res <- data.frame('gene'  =geneid,
                                     'coef'  =surv_scox_res$coefficients[,1],
                                     'HR'    =exp(surv_scox_res$coefficients[,1]),
                                     'C_index'=as.numeric(surv_scox_res$concordance[1]),
                                     'roc_AUC_NNK'=round(roc_res1$AUC,3),'roc_AUC_KM'=round(roc_res2$AUC,3),'roc_AUC_KMgroup'=round(roc_res3$AUC,3),
                                     'zph.p' =surv_scox_ph_p,
                                     'logtest.p'   =surv_scox_res$logtest[3],
                                     'waldtst.p'   =surv_scox_res$waldtest[3],
                                     'Scoretest.p'   =surv_scox_res$sctest[3],
                                     'km.p'=km_pvalue)
    
    if(pic){
      # 参考画图网址 ：
      p=ggplot(data = data.frame('FP'=roc_res2$FP,'TP'=roc_res2$TP)) +  geom_line(aes(x = FP, y = TP),color = "firebrick3",size=1) + 
        theme_bw()+labs(title = cancer)+ scale_x_continuous(name  = "FP")+ scale_y_continuous(name = "TP")+ geom_abline(slope=1, intercept=0,linetype='dotted')+
        annotate("text",x = .75, y = .25, label = paste("AUC:",round(roc_res2$AUC,2)),color = "blue",size=6)
      
      if(q3_km_p > km_pvalue){
        pic_out<- ggsurvplot( surv_fit(temp_formula2, data_df),data = data_df,pval =T,conf.int = F,risk.table = T,surv.median.line = 'hv',title =(paste(cancer,'HR:',round(singe_gene_cox_res$HR,2))))
      }else{
        pic_out<- ggsurvplot( surv_fit(temp_formula1, data_df),data = data_df,pval =T,conf.int = F,risk.table = T,surv.median.line = 'hv', title =(paste(cancer,'HR:',round(singe_gene_cox_res$HR,2))))
      } 
      return(list(pic_out,p))
    }else{
      return(singe_gene_cox_res)
    }
    
  }
  
  
  
}

