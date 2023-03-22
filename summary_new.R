library(Seurat)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(multtest)
library(gridExtra)
library(scHCL)
library(scMCA)
library(clusterProfiler)
library(org.Mm.eg.db)
merge_rds <- readRDS('/Users/kangyimei/Desktop/scRNA/cluster_after.rds')
library(ggsci)
cors <- c(pal_npg()(6))
cluster_mac <- readRDS('/Users/kangyimei/Desktop/scRNA/rds/cluster_mac.rds')

# subsets markers
cluster_mac$cell_type_orig <- paste0(cluster_mac$seurat_clusters, '_', cluster_mac$orig.ident)
mac_marker <- FindAllMarkers(cluster_mac)
cluster_0 <- mac_marker[mac_marker$cluster=='0',]
cluster_0 <- cluster_0[cluster_0$avg_log2FC>0, ]
c0_ID <- bitr(rownames(cluster_0), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Mm.eg.db)
c0_go <- enrichGO(c0_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'All',pvalueCutoff = 0.5, qvalueCutoff = 0.5)
c0_marker <- rownames(cluster_0)[1:10]
VlnPlot(cluster_mac,features = c0_marker ,pt.size = 0, group.by = 'seurat_clusters')

cluster_1 <- mac_marker[mac_marker$cluster=='1',]
cluster_1 <- cluster_1[cluster_1$avg_log2FC>0, ]
c1_ID <- bitr(cluster_1$gene, fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Mm.eg.db)
c1_go <- enrichGO(c1_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP',pvalueCutoff = 0.5, qvalueCutoff = 0.5)

cluster_2 <- mac_marker[mac_marker$cluster=='2',]
cluster_2 <- cluster_2[cluster_2$avg_log2FC>0, ]
c2_ID <- bitr(cluster_2$gene, fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Mm.eg.db)
c2_go <- enrichGO(c2_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP',pvalueCutoff = 0.5, qvalueCutoff = 0.5)

cluster_5 <- mac_marker[mac_marker$cluster=='5',]
cluster_5 <- cluster_5[cluster_5$avg_log2FC>0, ]
c5_ID <- bitr(cluster_5$gene, fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Mm.eg.db)
c5_go <- enrichGO(c5_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP',pvalueCutoff = 0.5, qvalueCutoff = 0.5)

cluster_6 <- mac_marker[mac_marker$cluster=='6',]
cluster_6 <- cluster_6[cluster_6$avg_log2FC>0, ]
c6_ID <- bitr(cluster_6$gene, fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Mm.eg.db)
c6_go <- enrichGO(c6_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP',pvalueCutoff = 0.5, qvalueCutoff = 0.5)

cluster_11 <- mac_marker[mac_marker$cluster=='11',]
cluster_11 <- cluster_11[cluster_11$avg_log2FC>0, ]
c11_ID <- bitr(cluster_11$gene, fromType = "SYMBOL",
              toType = c( "ENTREZID"),
              OrgDb = org.Mm.eg.db)
c11_go <- enrichGO(c11_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP',pvalueCutoff = 0.5, qvalueCutoff = 0.5)

# cluster ratio

table(cluster_mac$orig.ident)
table(Idents(cluster_mac), cluster_mac$orig.ident)
cellratio <- prop.table(table(Idents(cluster_mac),cluster_mac$orig.ident), margin = 2)
#plot
cellratio <- as.data.frame(cellratio)
ggplot(cellratio) + geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = 'identity', width = 0.7, size = 0.5) +
  theme_classic() + labs(x = 'Samples', y = 'Ratio') + 
  coord_flip() + theme(panel.border = element_rect(fill = NA, color = 'black', size = 0.5, linetype = 'solid'))


# split into clusters

sample_list <- SplitObject(cluster_mac, split.by = 'seurat_clusters')
cluster_0 <- sample_list$'0'
cluster_1 <- sample_list$'1'
cluster_2 <- sample_list$'2'
cluster_5 <- sample_list$'5'
cluster_6 <- sample_list$'6'
cluster_11 <- sample_list$'11'


#------------------------------------------------------------
#------------------------------------------------------------

#rhythmic genes total
Idents(cluster_mac) <- cluster_mac$orig.ident
All_rhythmic_genes <- FindMarkers(cluster_mac, ident.1 = 'WT4', ident.2 = 'WT16', logfc.threshold = 0.15)
All_rhythmic_genes <- All_rhythmic_genes[All_rhythmic_genes$p_val<0.05, ]
All_RG_up4 <- All_rhythmic_genes[All_rhythmic_genes$avg_log2FC>0, ]
All_RG_up4 <- All_RG_up4[All_RG_up4$p_val<0.05, ]
All_RG_up4_ID <- bitr(rownames(All_RG_up4), fromType = "SYMBOL",
                   toType = c( "ENTREZID"),
                   OrgDb = org.Mm.eg.db)
All_RG_up4_GO <- enrichGO(All_RG_up4_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP', pvalueCutoff = 0.5, qvalueCutoff = 0.5)
All_RG_up16 <- All_rhythmic_genes[All_rhythmic_genes$avg_log2FC<0, ]
All_RG_up16 <- All_RG_up16[All_RG_up16$p_val<0.05, ]
All_RG_up16_ID <- bitr(rownames(All_RG_up16), fromType = "SYMBOL",
                      toType = c( "ENTREZID"),
                      OrgDb = org.Mm.eg.db)
All_RG_up16_GO <- enrichGO(All_RG_up16_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP', pvalueCutoff = 0.5, qvalueCutoff = 0.5)

#RG in IRKO

AllKO_rhythmic_genes <- FindMarkers(cluster_mac, ident.1 = 'KO4', ident.2 = 'KO16', logfc.threshold = 0.15)
AllKO_rhythmic_genes <- AllKO_rhythmic_genes[AllKO_rhythmic_genes$p_val<0.05, ]
All_RGKO_up4 <- AllKO_rhythmic_genes[AllKO_rhythmic_genes$avg_log2FC>0, ]
All_RGKO_up4 <- All_RGKO_up4[All_RGKO_up4$p_val<0.05, ]
All_RGKO_up4_ID <- bitr(rownames(All_RGKO_up4), fromType = "SYMBOL",
                      toType = c( "ENTREZID"),
                      OrgDb = org.Mm.eg.db)
All_RGKO_up4_GO <- enrichGO(All_RGKO_up4_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP', pvalueCutoff = 0.5, qvalueCutoff = 0.5)
All_RGKO_up16 <- AllKO_rhythmic_genes[AllKO_rhythmic_genes$avg_log2FC<0, ]
All_RGKO_up16 <- All_RGKO_up16[All_RGKO_up16$p_val<0.05, ]
All_RGKO_up16_ID <- bitr(rownames(All_RGKO_up16), fromType = "SYMBOL",
                       toType = c( "ENTREZID"),
                       OrgDb = org.Mm.eg.db)
All_RGKO_up16_GO <- enrichGO(All_RGKO_up16_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP', pvalueCutoff = 0.5, qvalueCutoff = 0.5)

only <- c()
for(gene in rownames(All_rhythmic_genes)){
  if(gene %in% rownames(AllKO_rhythmic_genes)){
    only <- only
  }
  else{
    only <- c(only, gene)
  }
}

# RG in cluster0
WT_DEG0 <- FindMarkers(cluster_mac, ident.1 = '0_WT4', ident.2 = '0_WT16', logfc.threshold = 0.15)
WT_DEG0 <- WT_DEG0[WT_DEG0$p_val<0.05,]
KO_DEG0 <- FindMarkers(cluster_mac, ident.1 = '0_KO4', ident.2 = '0_KO16', logfc.threshold = 0.15)
KO_DEG0 <- KO_DEG0[KO_DEG0$p_val<0.05,]
only_DEG0 <- c()
for (gene in rownames(WT_DEG0)){
  if (gene %in% rownames(KO_DEG0)){
    only_DEG0 <- only_DEG0
  }
  else{
    only_DEG0 <- c(only_DEG0, gene)
  }
}

#RG in c1
WT_DEG1 <- FindMarkers(cluster_mac, ident.1 = '1_WT4', ident.2 = '1_WT16', logfc.threshold = 0.15)
WT_DEG1 <- WT_DEG1[WT_DEG1$p_val<0.05,]
KO_DEG1 <- FindMarkers(cluster_mac, ident.1 = '1_KO4', ident.2 = '1_KO16', logfc.threshold = 0.15)
KO_DEG1 <- KO_DEG1[KO_DEG1$p_val<0.05,]
only_DEG1 <- c()
for (gene in rownames(WT_DEG1)){
  if (gene %in% rownames(KO_DEG1)){
    only_DEG1 <- only_DEG1
  }
  else{
    only_DEG1 <- c(only_DEG1, gene)
  }
}

# rhythmic genes in cluster2
WT_DEG2 <- FindMarkers(cluster_mac, ident.1 = '2_WT4', ident.2 = '2_WT16', logfc.threshold = 0.15)
WT_DEG2 <- WT_DEG2[WT_DEG2$p_val<0.05,]
#----------------
up16_c2 <- WT_DEG2[WT_DEG2$avg_log2FC<0,]
up16_c2_ID <- bitr(rownames(up16_c2), fromType = "SYMBOL",
                   toType = c( "ENTREZID"),
                   OrgDb = org.Mm.eg.db)
up16_c2_go <- enrichGO(up16_c2_ID$ENTREZID, OrgDb = org.Mm.eg.db, ont = 'BP',pvalueCutoff = 0.5, qvalueCutoff = 0.5)
ego2 <- pairwise_termsim(up16_c2_go)
p1 <- emapplot(ego2)
WTKO16_DEG2 <- FindMarkers(cluster_mac, ident.1 = '2_WT16', ident.2 = '2_KO16', logfc.threshold = 0.15)
#-------------------
KO_DEG2 <- FindMarkers(cluster_mac, ident.1 = '2_KO4', ident.2 = '2_KO16', logfc.threshold = 0.15)
KO_DEG2 <- KO_DEG2[KO_DEG2$p_val<0.05, ]
only_DEG2 <- c()
for (gene in rownames(WT_DEG2)){
  if (gene %in% rownames(KO_DEG2)){
    only_DEG2 <- only_DEG2
  }
  else{
    only_DEG2 <- c(only_DEG2, gene)
  }
}

#RG in c5
WT_DEG5 <- FindMarkers(cluster_mac, ident.1 = '5_WT4', ident.2 = '5_WT16', logfc.threshold = 0.15)
WT_DEG5 <- WT_DEG5[WT_DEG5$p_val<0.05,]
WT_DEG5 <- WT_DEG5[WT_DEG5$p_val<0.05,]
KO_DEG5 <- FindMarkers(cluster_mac, ident.1 = '5_KO4', ident.2 = '5_KO16', logfc.threshold = 0.15)
KO_DEG5 <- KO_DEG5[KO_DEG5$p_val<0.05, ]
only_DEG5 <- c()
for (gene in rownames(WT_DEG5)){
  if (gene %in% rownames(KO_DEG5)){
    only_DEG5 <- only_DEG5
  }
  else{
    only_DEG5 <- c(only_DEG5, gene)
  }
}

#RG in c6
WT_DEG6 <- FindMarkers(cluster_mac, ident.1 = '6_WT4', ident.2 = '6_WT16', logfc.threshold = 0.15)
WT_DEG6 <- WT_DEG6[WT_DEG6$p_val<0.05,]
KO_DEG6 <- FindMarkers(cluster_mac, ident.1 = '6_KO4', ident.2 = '6_KO16', logfc.threshold = 0.15)
KO_DEG6 <- KO_DEG6[KO_DEG6$p_val<0.05, ]
only_DEG6 <- c()
for (gene in rownames(WT_DEG6)){
  if (gene %in% rownames(KO_DEG6)){
    only_DEG6 <- only_DEG6
  }
  else{
    only_DEG6 <- c(only_DEG6, gene)
  }
}

#RG in c11
WT_DEG11 <- FindMarkers(cluster_mac, ident.1 = '11_WT4', ident.2 = '11_WT16', logfc.threshold = 0.15)
WT_DEG11 <- WT_DEG11[WT_DEG11$p_val<0.05,]
KO_DEG11 <- FindMarkers(cluster_mac, ident.1 = '11_KO4', ident.2 = '11_KO16', logfc.threshold = 0.15)
KO_DEG11 <- KO_DEG11[KO_DEG11$p_val<0.05, ]
only_DEG11 <- c()
for (gene in rownames(WT_DEG11)){
  if (gene %in% rownames(KO_DEG11)){
    only_DEG11 <- only_DEG11
  }
  else{
    only_DEG11 <- c(only_DEG11, gene)
  }
}

# RG ratio
RG_ratio <- read.csv(file = '/Users/kangyimei/Desktop/scRNA/RG_ratio.csv', header = T)
RG_ratio <- as.data.frame(RG_ratio)
ggplot(RG_ratio) + geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = 'identity', width = 0.7, size = 0.5) +
  theme_classic() + labs(x = 'Samples', y = 'Ratio') + 
  coord_flip() + theme(panel.border = element_rect(fill = NA, color = 'black', size = 0.5, linetype = 'solid'))

#B target ratio
B_ratio <- read.csv(file = '/Users/kangyimei/Desktop/scRNA/B_target.csv', header = T)
B_ratio <- as.data.frame(B_ratio)
ggplot(B_ratio) + geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = 'identity', width = 0.7, size = 0.5) +
  theme_classic() + labs(x = 'Samples', y = 'Ratio') + 
  coord_flip() + theme(panel.border = element_rect(fill = NA, color = 'black', size = 0.5, linetype = 'solid'))

#----------------------------------------------------------------------------------
# bmal putative targets
B_targets <- read.csv(file = '/Users/kangyimei/Desktop/scRNA/sci-reports/lps0.csv', header = T)
B_targets <- B_targets[1:10000,]
genes <- c()
list <- B_targets$symbol
for (gene in B_targets$symbol){
  if(gene %in% genes){
    genes <- genes
  }
  else{
    genes <- c(genes, gene)
  }
}

#---------------------------------------------------------------------------------
#elife bulk seq
BKO_data <- read.csv('/Users/kangyimei/Desktop/scRNA/elife_bulkseq/elife_bulkseq_BKO.csv', header = T, row.names = 1)
tpm_only_IFN <- data.frame(WT_IFN_1 = BKO_data$WT_Ifng_only_S1...TPM, WT_IFN_2 = BKO_data$WT_Ifng_only_S2...TPM, WT_IFN_3 = BKO_data$WT_Ifng_only_S3...TPM,
                           BKO_IFN_1 = BKO_data$M.BKO_Ifng_only_S1...TPM, BKO_IFN_2 = BKO_data$M.BKO_Ifng_only_S2...TPM, BKO_IFN_3 = BKO_data$M.BKO_Ifng_only_S3...TPM)
rownames(tpm_only_IFN) <- rownames(BKO_data)

group_list=c(rep('WT_IFN',3),rep('BKO_IFN',3))
## 强制限定顺序
group_list <- factor(group_list,levels = c("WT_IFN","BKO_IFN"),ordered = F)
#表达矩阵数据校正
exprSet <- tpm_only_IFN
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
#差异分析：
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
#save(deg,file = 'deg.Rdata')

## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
if(T){
  logFC_t=0.15
  deg$condition=ifelse(deg$logFC>0,'up_BKO','up_WT')
  deg$g=ifelse(deg$P.Value>0.05,'stable',
               ifelse(deg$logFC > logFC_t,'UP',
                      ifelse(deg$logFC < -logFC_t,'DOWN','stable') )
  )
  table(deg$g)
  head(deg)
  deg$symbol=rownames(deg)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Mm.eg.db)
  head(df)
  DEG=deg
  head(DEG)
  
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  head(DEG)
  
  save(DEG,file = 'anno_DEG.Rdata')
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
}

# 最简单的超几何分布检验
###这里就拿KEGG数据库举例吧，拿自己判定好的上调基因集进行超几何分布检验，如下
if(T){
  gene_down
  gene_up
  enrichKK <- enrichKEGG(gene         =  gene_up,
                         organism     = 'mmu',
                         #universe     = gene_all,
                         pvalueCutoff = 0.05,
                         qvalueCutoff =0.05)
  head(enrichKK)[,1:6] 
  browseKEGG(enrichKK, 'hsa04512')
  dotplot(enrichKK)
  ggsave("enrichKK.png")
  enrichKK=DOSE::setReadable(enrichKK, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
  enrichKK 
}
##最基础的条形图和点图
#条带图
barplot(enrichKK,showCategory=20)
#气泡图
dotplot(enrichKK)
# GO 分析
GO_up <- enrichGO(gene_up, OrgDb = org.Mm.eg.db, ont = 'BP',pvalueCutoff = 0.5, qvalueCutoff = 0.5)
DEG_BK <- DEG[DEG$g != 'stable', ]
