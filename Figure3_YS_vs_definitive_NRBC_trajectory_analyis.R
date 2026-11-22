######################################################################################################################################################################
########## ----------------------trajectory ananlysis---------------##################################################################################################
######################################################################################################################################################################
 

set.seed(123)
# 设置全局变量内存大小，默认是500M，我们数据超过，需要重新设置
options(future.globals.maxSize=50000 * 1024 ^ 2 )

library(corrplot)
library(clusterProfiler)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(limma)# 使用strsplit2函数
#  变量 col，颜色
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(SeuratDisk)
library(monocle)
library(RColorBrewer)
library(Seurat)
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
setwd('~/2021_NRBC_chlyu/primitive_NRBC_project')

NRBC_altas_seu=readRDS('~/2021_NRBC_chlyu/primitive_NRBC_project/all_NRBC_altas_seu.rds')
Idents(NBRC_altas_seu)='final_celltype'
celltype_markers=FindAllMarkers(NBRC_altas_seu)
cho_genes=celltype_markers %>% filter(p_val_adj<0.000001) %>% group_by(cluster) %>% top_n(n = 300,wt = avg_log2FC); 
NBRC_altas_seu$final_celltype=factor(NBRC_altas_seu$final_celltype,levels = c("BFUE/CFUE" ,"EryP-CFC" ,"ProE", "Bas","Poly","Orth"))

write.csv(celltype_markers,file = 'res_data/all_nrbc_celltype_markers_df.csv')

# 各类别中要YS,FBM,ABM,UCB,FL各抽取200各细胞；
sample_cellids=c(sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='EryP-CFC' & NBRC_altas_seu$tissue_stage=='YS',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Bas' & NBRC_altas_seu$tissue_stage=='YS',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Poly' & NBRC_altas_seu$tissue_stage=='YS',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Orth' & NBRC_altas_seu$tissue_stage=='YS',]),size = 200),# 8 个
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='BFUE/CFUE' & NBRC_altas_seu$tissue_stage=='FL',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='ProE' & NBRC_altas_seu$tissue_stage=='FL',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Bas' &  NBRC_altas_seu$tissue_stage=='FL',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Poly' & NBRC_altas_seu$tissue_stage=='FL',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Orth' & NBRC_altas_seu$tissue_stage=='FL',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='BFUE/CFUE' & NBRC_altas_seu$tissue_stage=='FBM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='ProE' & NBRC_altas_seu$tissue_stage=='FBM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Bas' &  NBRC_altas_seu$tissue_stage=='FBM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Poly' & NBRC_altas_seu$tissue_stage=='FBM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Orth' & NBRC_altas_seu$tissue_stage=='FBM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='BFUE/CFUE' & NBRC_altas_seu$tissue_stage=='ABM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='ProE' & NBRC_altas_seu$tissue_stage=='ABM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Bas' &  NBRC_altas_seu$tissue_stage=='ABM',]),size = 200),
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Poly' & NBRC_altas_seu$tissue_stage=='ABM',]),size = 200), 
                 sample(rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$final_celltype=='Orth' & NBRC_altas_seu$tissue_stage=='ABM',]),size = 200))         

nrbc_sampled_seu=subset(NBRC_altas_seu,cells=sample_cellids)
pd <- new("AnnotatedDataFrame", data = nrbc_sampled_seu@meta.data)
fd=data.frame(gene_short_name=rownames(nrbc_sampled_seu),row.names = rownames(nrbc_sampled_seu))
fd <- new("AnnotatedDataFrame", data = fd)

nrbc_sampled_cds <- newCellDataSet(GetAssayData(nrbc_sampled_seu), phenoData = pd, featureData = fd)
nrbc_sampled_cds=estimateSizeFactors(nrbc_sampled_cds)
nrbc_sampled_cds=estimateDispersions(nrbc_sampled_cds)
nrbc_sampled_cds=setOrderingFilter(nrbc_sampled_cds,ordering_genes = unique(cho_genes$gene)) # 这一步很重要，基因越多，在细胞比较多的情况下，速度也慢
plot_ordering_genes(nrbc_sampled_cds)

nrbc_sampled_cds=reduceDimension(nrbc_sampled_cds,reduction_method = 'DDRTree')
nrbc_sampled_cds=orderCells(nrbc_sampled_cds)

p1=plot_cell_trajectory(nrbc_sampled_cds,color_by='final_celltype',cell_size = 0.5)+
  scale_colour_manual(values = alpha(brewer.pal(name = 'Set3',n = 9),alpha = 0.6) )+theme(text = element_text(face = 'bold'))

p2=plot_cell_trajectory(nrbc_sampled_cds,color_by='tissue_stage',)+
  scale_colour_manual(values = alpha(brewer.pal(name = 'Set3',n = 9),alpha = 0.6) )+theme(text = element_text(face = 'bold'))
p=p1+p2
p
ggsave(p,width = 8,height = 5,filename='res_pic/main_figure3/monocle_all_tissue_NRBC_trajectory_celltype_tissue_info.pdf')

p=plot_cell_trajectory( nrbc_sampled_cds, color_by = "Pseudotime",cell_size = 0.5)+
   plot_cell_trajectory( nrbc_sampled_cds, color_by = "State",cell_size = 0.5)
p
ggsave(p,width = 8,height = 5,filename='res_pic/main_figure3/monocle_all_tissue_NRBC_trajectory_Pseudotime_State_info.pdf')

maturation_differentiation_genes=c( "HBA1", "HBA2",'ALAS2','FECH','HMBS','SLC25A37','AHSP','SLC4A1','ANK1','EPB42') 
cds_subset <- nrbc_sampled_cds[maturation_differentiation_genes[2],]
plot_genes_in_pseudotime(cds_subset, color_by = "final_celltype")+facet_grid(~tissue_stage)+scale_color_manual(values = cols[1:7])

saveRDS(nrbc_sampled_cds,file = 'all_sub_nrbc_sampled_cds.rds')


nrbc_sampled_cds=readRDS('all_sub_nrbc_sampled_cds.rds')
UCell_score_metad_df=readRDS('NRBC_altas_seu_metadata_UCell_score.rds')
UCell_score_metad_df=UCell_score_metad_df[,-grep('cortical',colnames(UCell_score_metad_df))]
head(UCell_score_metad_df)

cho_pathways=c('erythrocyte_differentiation_UCell','oxygen_transport_UCell','hemoglobin_biosynthetic_process_UCell','heme_biosynthetic_process_UCell',
               'Cell_redox_homeostasis_UCell','response_to_metal_ion_UCell','stress_response_to_metal_ion_UCell','RNA_splicing_UCell','Golgi_vesicle_transport_UCell','Proteasome-mediated_protein_catabolic_process_UCell')
p=pheatmap(cor(UCell_score_metad_df[,cho_pathways]),display_numbers = T,
           color =colorRampPalette(colors = c('navy','white','firebrick3'))(100) )

ggsave(as.ggplot(p),filename='res_pic/main_figure3/all_NRBC_pathway_correlation_heatmap.pdf',width = 12,height = 12)


p=pheatmap(cor(UCell_score_metad_df[UCell_score_metad_df$tissue_stage=='YS',cho_pathways]),display_numbers = T,
           main = 'YS NRBC',color =colorRampPalette(colors = c('navy','white','firebrick3'))(100) )

ggsave(as.ggplot(p),filename='res_pic/main_figure3/YS_NRBC_pathway_correlation_heatmap.pdf',width = 12,height = 12)

p=pheatmap(cor(UCell_score_metad_df[UCell_score_metad_df$tissue_stage=='FL',cho_pathways]),display_numbers = T,
           main = 'FL NRBC',color =colorRampPalette(colors = c('navy','white','firebrick3'))(100) )

ggsave(as.ggplot(p),filename='res_pic/main_figure3/FL_NRBC_pathway_correlation_heatmap.pdf',width = 12,height = 12)

p=pheatmap(cor(UCell_score_metad_df[UCell_score_metad_df$tissue_stage=='FBM',cho_pathways]),display_numbers = T,
           main = 'FBM NRBC',color =colorRampPalette(colors = c('navy','white','firebrick3'))(100) )

ggsave(as.ggplot(p),filename='res_pic/main_figure3/FBM_NRBC_pathway_correlation_heatmap.pdf',width = 12,height = 12)

p=pheatmap(cor(UCell_score_metad_df[UCell_score_metad_df$tissue_stage=='ABM',cho_pathways]),display_numbers = T,
           main = 'ABM NRBC',color =colorRampPalette(colors = c('navy','white','firebrick3'))(100) )

ggsave(as.ggplot(p),filename='res_pic/main_figure3/ABM_NRBC_pathway_correlation_heatmap.pdf',width = 12,height = 12)

ave_df=AverageExpression(NRBC_altas_seu,group.by = 'source_celltype')$RNA
gsvapar <- gsvaParam(as.matrix(AML_count_nr_df), list(ISR_genes=ISR_genes))
ISR_score <- gsva(gsvapar)


cho_cols=c('erythrocyte_differentiation_UCell','oxygen_transport_UCell','heme_biosynthetic_process_UCell','hemoglobin_biosynthetic_process_UCell',
           'stress_response_to_metal_ion_UCell','cellular_oxidant_detoxification_UCell',
           'integrated_stress_response_signaling_UCell','HRI_mediated_signaling_UCell','Mitophagy_UCell')

pheatmap(cor(UCell_score_metad_df[UCell_score_metad_df$tissue_stage=='YS',cho_cols]),display_numbers = T)


cho_cols=c('oxygen_transport_UCell','heme_biosynthetic_process_UCell','hemoglobin_biosynthetic_process_UCell',
           'response_to_metal_ion_UCell','stress_response_to_metal_ion_UCell','response_to_metal_ion_UCell','cellular_oxidant_detoxification_UCell',
           'integrated_stress_response_signaling_UCell','HRI_mediated_signaling_UCell','Mitophagy_UCell')

cho_UCell_score_metad_df=UCell_score_metad_df[,cho_cols]

nrbc_sampled_cds@phenoData@data[,cho_cols]=UCell_score_metad_df[rownames(nrbc_sampled_cds@phenoData@data),cho_cols]


plot_cell_trajectory(nrbc_sampled_cds,color_by = 'oxygen_transport_UCell',use_color_gradient = T,cell_size = 0.5)+facet_grid(~tissue_stage)
plot_genes_in_pseudotime(nrbc_sampled_cds, color_by = "oxygen_transport_UCell")

cho_UCell_score_metad_df=UCell_score_metad_df[rownames(nrbc_sampled_cds@phenoData@data),]
cho_UCell_score_metad_df$Pseudotime=nrbc_sampled_cds@phenoData@data[,'Pseudotime']
cho_UCell_score_metad_df$DDRT_X=nrbc_sampled_cds@reducedDimS[1,]
cho_UCell_score_metad_df$DDRT_Y=nrbc_sampled_cds@reducedDimS[2,]


p1=ggplot(cho_UCell_score_metad_df,aes(x=Pseudotime,y=oxygen_transport_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)
p2=ggplot(cho_UCell_score_metad_df,aes(x=Pseudotime,y=heme_biosynthetic_process_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)
p3=ggplot(cho_UCell_score_metad_df,aes(x=Pseudotime,y=hemoglobin_biosynthetic_process_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)

p=p1+p2+p3
p
ggsave(p,width =18 ,height =5 ,filename='res_pic/main_figure3/monocle_all_tissue_ucell_score_pseudotime.pdf')

p4=ggplot(cho_UCell_score_metad_df,aes(x=Pseudotime,y=stress_response_to_metal_ion_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)
p5=ggplot(cho_UCell_score_metad_df,aes(x=Pseudotime,y=cellular_oxidant_detoxification_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)
p6=ggplot(cho_UCell_score_metad_df,aes(x=Pseudotime,y= HRI_mediated_signaling_UCell ,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)
p7=ggplot(cho_UCell_score_metad_df,aes(x=Pseudotime,y= Mitophagy_UCell ,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)

p=p4+p5+p6+p7
p
ggsave(p,width =12 ,height =10 ,filename='res_pic/main_figure3/monocle_all_tissue_tress_ucell_score_pseudotime_stress.pdf')



stress_df=cho_UCell_score_metad_df[,c('tissue_stage','Pseudotime','HRI_mediated_signaling_UCell','integrated_stress_response_signaling_UCell')]
stress_df=stress_df[order(stress_df$tissue_stage,stress_df$Pseudotime),]

cho_UCell_score_metad_df=cho_UCell_score_metad_df[order(cho_UCell_score_metad_df$tissue_stage,cho_UCell_score_metad_df$Pseudotime),]
p=ggplot(cho_UCell_score_metad_df,aes(x=heme_biosynthetic_process_UCell,y=cellular_oxidant_detoxification_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)+stat_cor()
p
ggsave(p,width =5 ,height =5 ,filename='res_pic/main_figure3/stress_score_pseudotime_correlation.pdf')


p=ggplot(cho_UCell_score_metad_df,aes(x=HRI_mediated_signaling_UCell,y=cellular_oxidant_detoxification_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)+stat_cor()
p

p=ggplot(cho_UCell_score_metad_df,aes(x=stress_response_to_metal_ion_UCell,y=cellular_oxidant_detoxification_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)+stat_cor()
p





ggsave(p,width =5 ,height =5 ,filename='res_pic/main_figure3/stress_score_pseudotime_correlation.pdf')

p1=ggplot(cho_UCell_score_metad_df,aes(x=HRI_mediated_signaling_UCell,y=oxygen_transport_UCell ,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)+
  stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4)
p1

p2=ggplot(cho_UCell_score_metad_df,aes(x=HRI_mediated_signaling_UCell,y=stress_response_to_metal_ion_UCell,fill=tissue_stage))+theme_classic()+
  scale_color_manual(values = cols)+ geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,color=tissue_stage ),linewidth = 1)+
  stat_cor( method = "pearson",label.x.npc = "left",label.y.npc = "top",   size = 4)
p=p1/p2;p
ggsave(p,width =5 ,height =10 ,filename='res_pic/main_figure3/stress_score_pseudotime_correlation.pdf')



YS_pval=cor.test(cho_UCell_score_metad_df$HRI_mediated_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='YS'],
         cho_UCell_score_metad_df$integrated_stress_response_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='YS'])
FL_pval=cor.test(cho_UCell_score_metad_df$HRI_mediated_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='FL'],
         cho_UCell_score_metad_df$integrated_stress_response_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='FL'])
FBM_pval=cor.test(cho_UCell_score_metad_df$HRI_mediated_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='FBM'],
         cho_UCell_score_metad_df$integrated_stress_response_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='FL'])
ABM_pval=cor.test(cho_UCell_score_metad_df$HRI_mediated_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='ABM'],
         cho_UCell_score_metad_df$integrated_stress_response_signaling_UCell[cho_UCell_score_metad_df$tissue_stage=='FL'])




p=ggplot(cho_UCell_score_metad_df[cho_UCell_score_metad_df$tissue_stage=='YS',],aes(x=Pseudotime,fill=tissue_stage))+
  theme_classic()+ 
  geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,y=stress_response_to_metal_ion_UCell, ),linewidth = 1)+
  geom_smooth(method = "loess", se = FALSE, aes(group = tissue_stage,y=HRI_mediated_signaling_UCell ),linewidth = 1)

p

ggsave(p,width =5 ,height =5 ,filename='res_pic/main_figure3/stress_score_along_pseudotime.pdf')


##############################################################################################################################################################
#---------------------------------------------BEAM analysis compare primitive Bas-Orth with definitive Bas-Orth---------------------------#
##############################################################################################################################################################
fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch"
reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)"
branch_states = NULL
branch_point=1
relative_expr = TRUE

cds=buildBranchCellDataSet2(cds =nrbc_sampled_cds,progenitor_method ='duplicate',branch_point = 1   )
branchTest_res <- differentialGeneTest(cds, 
                                       fullModelFormulaStr = fullModelFormulaStr, 
                                       reducedModelFormulaStr = reducedModelFormulaStr, 
                                       cores = 4, 
                                       relative_expr = T, 
                                       verbose=T)

branchTest_res <- branchTest(nrbc_sampled_cds, branch_point=1,    cores = 4)

cmbn_df <- branchTest_res[, 1:4] 


#make a newCellDataSet object with the smoothed data? 

fd <- fData(cds)[row.names(cmbn_df),]

#combined dataframe: 
cmbn_df <- cbind(cmbn_df, fd)
write.csv(cmbn_df,'res_data/BEAM_res.csv')

cmbn_df=cmbn_df[order(cmbn_df$qval),]
BEAM_resBEAM_res <- cmbn_df[order(cmbn_df$qval),]
BEAM_res <- BEAM_res[BEAM_res$qval <0.001,]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
dim(BEAM_res) # 3233

library(clusterProfiler)
library(enrichplot)

pd_whole_level_markers=readRDS('res_data/pd_whole_level_markers.rds')
BEAM_pUp_enrichGO_res=enrichGO(gene = cmbn_df$gene_short_name[ cmbn_df$gene_short_name %in% pd_whole_level_markers$gene[pd_whole_level_markers$cluster=='primitive' & pd_whole_level_markers$avg_log2FC >0] ],
                           OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP',pvalueCutoff = 0.05)
BEAM_pDn_enrichGO_res=enrichGO(gene = cmbn_df$gene_short_name[ cmbn_df$gene_short_name %in% pd_whole_level_markers$gene[pd_whole_level_markers$cluster=='primitive' & pd_whole_level_markers$avg_log2FC <0] ],
                               OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP',pvalueCutoff = 0.05)

dotplot(BEAM_pUp_enrichGO_res)
dotplot(BEAM_pDn_enrichGO_res)
# 同DEGs富集分析结果基本一致

#1. 红系分化与成熟（Erythroid differentiation and maturation）
#GO ID	GO Biological Process	适用方向
#GO:0030218	erythrocyte differentiation	红系分化
#GO:0043249	erythrocyte maturation	红细胞成熟
#GO:0048821	erythrocyte development	红系发育全过程
#GO:0043362	nucleate erythrocyte maturation	有核红细胞成熟
#
#2. 血红蛋白及血红素合成（Hemoglobin and heme biosynthesis）
#GO ID	GO Biological Process	适用方向
#GO:0042541	hemoglobin biosynthetic process	血红蛋白合成
#GO:0006783	heme biosynthetic process	血红素合成
#GO:0006785	heme B biosynthetic process	血红蛋白所用血红素类型的合成
#GO:0070453	regulation of heme biosynthetic process	血红素合成调控

#3. 红系去核（Erythroblast enucleation）
#GO ID	GO Biological Process	适用方向
#GO:0043131	erythrocyte enucleation	红系细胞去核
#GO:0061930	regulation of erythrocyte enucleation	去核过程调控
#GO:0061931	positive regulation of erythrocyte enucleation	促进去核
#GO:0090601	enucleation	广义细胞去核

#4. 红细胞膜骨架组织与重塑（Membrane and cytoskeletal remodeling）
#GO ID	GO术语	GO类型
#GO:0030036	actin cytoskeleton organization	BP
#GO:0007009	plasma membrane organization	BP
#GO:0030866	cortical actin cytoskeleton organization	BP
#GO:0007015	actin filament organization	BP
#GO:0032956	regulation of actin cytoskeleton organization	BP
#GO:0007010	cytoskeleton organization	BP
#GO:0061024	membrane organization	BP
#GO:0014731	spectrin-associated cytoskeleton	CC ***
#
#红系分化与成熟
#GO:0030218
#erythrocyte differentiation
#GO:0043249
#erythrocyte maturation
#GO:0043362 ，空
#nucleate erythrocyte maturation
#功能成熟与血红蛋白合成
#GO:0015671
#oxygen transport
#GO:0042541
#hemoglobin biosynthetic process
#GO:0006783
#heme biosynthetic process
#红系去核
#GO:0043131
#erythrocyte enucleation
#GO:0061930
#regulation of erythrocyte enucleation
#膜骨架组织与重塑
#GO:0030036
#actin cytoskeleton organization
#GO:0007009
#plasma membrane organization
#GO:0030866
#cortical actin cytoskeleton organization
#

cho_GOids=list('erythrocyte_maturation'=c('GO:0030218','GO:0043249'),
               'oxygen_function'=c( 'GO:0015671','GO:0042541','GO:0006783'),
               'erythrocyte_enucleation'=c('GO:0043131','GO:0061930'),
               'membrane_organization'=c('GO:0030036','GO:0007009','GO:0030866'))

cho_gene_list=mapIds(x = org.Hs.eg.db,keys =cho_GOids[['oxygen_function']],column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
plot_genes_branched_heatmap(nrbc_sampled_cds[row.names(BEAM_res)[row.names(BEAM_res) %in% unlist(cho_gene_list)],],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)
#pd_monocle_BEAM_oxygen_transport_geneset_heatmap.pdf

cho_gene_list=mapIds(x = org.Hs.eg.db,keys =cho_GOids[['erythrocyte_maturation']],column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
plot_genes_branched_heatmap(nrbc_sampled_cds[row.names(BEAM_res)[row.names(BEAM_res) %in% unlist(cho_gene_list)],],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)

# pd_monocle_BEAM_mature_geneset_heatmap.pdf, 5x5


cho_gene_list=mapIds(x = org.Hs.eg.db,keys =cho_GOids[['erythrocyte_enucleation']],column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
plot_genes_branched_heatmap(nrbc_sampled_cds[row.names(BEAM_res)[row.names(BEAM_res) %in% unlist(cho_gene_list)],],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)
#pd_monocle_BEAM_erythrocyte_enucleation_geneset_heatmap.pdf


cho_gene_list=mapIds(x = org.Hs.eg.db,keys =cho_GOids[['membrane_organization']],column = 'SYMBOL',keytype = 'GO',multiVals = 'list')
plot_genes_branched_heatmap(nrbc_sampled_cds[row.names(BEAM_res)[row.names(BEAM_res) %in% unlist(cho_gene_list)],],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)
#pd_monocle_BEAM_membrane_organization_heatmap.pdf

plot_cell_trajectory( nrbc_sampled_cds, markers  = "HIPK2",cell_size = 0.5,use_color_gradient = T) +  scale_color_gradient(low = "lightgrey", high = "red")+facet_grid(~tissue_stage)

VlnPlot(NBRC_altas_seu,features = 'STAT5B',group.by = 'source_celltype',pt.size = 0,split.by = 'tissue_stage')+NoLegend()





