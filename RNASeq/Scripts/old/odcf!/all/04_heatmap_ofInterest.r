#plot arid1a figures
---
title: 'Ali''s RNAseq: All samples together:04_DiffExpr'
output:
  html_document:
    theme: lumen
    highlight: pygments
    toc_float: true
    toc: true
    fig_align: left
---
Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(knitr)
library(DT)
library(clusterProfiler)
library(Glimma)
library(limma)
library(edgeR)
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Read in Data
dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

anno <- colData(dds_genes)

#get log trans
vst <-vst(dds_genes)

#extract samples of interest
#wt vs arid1a
#filter samples
anno_sub <- anno[which(anno$genotype %in% c("ARID1A", "WT") & anno$tamoxifen_treatment =="untreated" & anno$hdac_treatment=="untreated"  
& anno$irradiation_treatment=="untreated"), ]
#filter DEG
DEG_results_list_sub <- DEG_results_list$WT_Tam.untreated_Irr.untreated_Hdac.untreated_vs_ARID1A_Tam.untreated_Irr.untreated_Hdac.untreated
DEG_results_list_sub_sub <- DEG_results_list_sub[which(DEG_results_list_sub$padj<0.05 & abs(DEG_results_list_sub$log2FoldChange)>1), ] 
dim(DEG_results_list_sub_sub)
#filter expression
vst_sub <- assay(vst)[rownames(DEG_results_list_sub_sub),rownames(anno_sub)]

#plot heatmap
annovst <- as.data.frame(anno_sub)[, c("genotype", "tamoxifen_treatment")] 
pheatmap(vst_sub, scale="row", show_colnames=F,
                       annotation_col=annovst, show_rownames=F, file=file.path(PostDE.dir, "WT_Tam.untreated_Irr.untreated_Hdac.untreated_vs_ARID1A_Tam.untreated_Irr.untreated_Hdac.untreated", "WT_vs_Arid1A_heat.pdf" )) 


#wt treated vs arid1 treated
#filter samples
anno_sub <- anno[which(anno$genotype %in% c("ARID1A", "WT") & anno$tamoxifen_treatment =="treated" & anno$hdac_treatment=="untreated" & 
 anno$irradiation_treatment=="untreated"), ]
#filter DEG
DEG_results_list_sub <- DEG_results_list$WT_Tam.treated_Irr.untreated_Hdac.untreated_vs_ARID1A_Tam.treated_Irr.untreated_Hdac.untreated
DEG_results_list_sub_sub <- DEG_results_list_sub[which(DEG_results_list_sub$padj<0.05 & abs(DEG_results_list_sub$log2FoldChange)>1), ] 
dim(DEG_results_list_sub_sub)
#filter expression
vst_sub <- assay(vst)[rownames(DEG_results_list_sub_sub),rownames(anno_sub)]

#plot heatmap
annovst <- as.data.frame(anno_sub)[, c("genotype", "tamoxifen_treatment")] 
pheatmap(vst_sub, scale="row", show_colnames=F,
                       annotation_col=annovst, show_rownames=F, file=file.path(PostDE.dir, "WT_Tam.untreated_Irr.untreated_Hdac.untreated_vs_ARID1A_Tam.untreated_Irr.untreated_Hdac.untreated", "WT_vs_Arid1A_TamTreated_heat.pdf" )) 
