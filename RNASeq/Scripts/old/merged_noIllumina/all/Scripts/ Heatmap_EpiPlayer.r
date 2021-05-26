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
library(dplyr)

#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190531_odcf_NoIll/"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
rownames(dds_genes)<- sapply(strsplit(rownames(dds_genes), ".", fixed=TRUE), "[", 1)
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
#DEG_results_list_sign <- lapply(DEG_results_list, function(x)x[which(x$padj < 0.05 & abs(x$log2FoldChange)>1),])
DEG_results_list_sign <- lapply(DEG_results_list, function(x)x[which(x$padj < 0.05),])

vst_genes <- vst(dds_genes)


#load epigenetic players
epi_player <- read.csv(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data","data", "HumanEpiPlayer.csv"), header=TRUE, stringsAsFactors = FALSE)
epi_player$ensembl<- mapIds(org.Hs.eg.db, keys=epi_player$gene.name, keytype ="SYMBOL", column = "ENSEMBL", multiVals = "first" )

#subset of dna damage genes
vst_genes <- assay(vst_genes)
vst_damage <-vst_genes[rownames(vst_genes) %in% epi_player$ensembl,]
vst_damage<-as.data.frame(vst_damage)
vst_damage$ensembl <- rownames(vst_damage)
vst_damage <- left_join(vst_damage, epi_player)
rownames(vst_damage)<- vst_damage$ensembl
dim(vst_damage)

#subset samples to plot
anno <- colData(dds_genes)
anno_sub <- anno[which(anno$genotype %in% c("WT","ARID1A" )& 
    anno$tamoxifen_treatment=="untreated" &
    anno$hdac_treatment=="untreated" &
    anno$irradiation_treatment=="untreated"),]
#vst_damage <- vst_damage[,rownames(anno_sub)]

#subset sign genes
vst_sub_sign <- vst_damage[which(rownames(vst_damage) %in% c(DEG_results_list_sign$WT_Tam.untreated_Irr.untreated_Hdac.untreated_vs_ARID1A_Tam.untreated_Irr.untreated_Hdac.untreated$ensembl)),]
#,
# DEG_results_list_sign$WT_Tam.untreated_Irr.treated_Hdac.untreated_vs_ID3_Tam.untreated_Irr.treated_Hdac.untreated$ensembl ,
# DEG_results_list_sign$ID3_Tam.untreated_Irr.treated_Hdac.untreated_vs_ID3_rescue_Tam.untreated_Irr.untreated_Hdac.untreated$ensembl
dim(vst_sub_sign)
vst_sub_sign <- as.data.frame(vst_sub_sign)

#heatmap annotation
#columns
anno_vst <- anno_sub[,c("genotype", "irradiation_treatment")]
#rows
rownames(vst_sub_sign) <- vst_sub_sign$gene.name
#anno_row <- vst_sub_sign[,c("Pathway"), drop=F]
#anno_names<- ifelse(is.na(anno_names),vst_sub_sign$ensembl, anno_names )

#plot heatmap
dir.create(file.path(PreDE.dir , "EpiPlayer"))
pheatmap(vst_sub_sign[, rownames(anno_sub)], annotation_col=as.data.frame(anno_vst) ,scale="row", clustering_distance_cols="correlation",#annotation_row=anno_row,
show_colnames=F, show_rownames=T,file=file.path(PreDE.dir , "EpiPlayer","EpiPlayer_correlation_WtUT_vs_Arid1aUT_DEG.pdf"),fontsize_row=4)
pheatmap(vst_sub_sign[, rownames(anno_sub)], annotation_col=as.data.frame(anno_vst) ,scale="row", clustering_distance_cols="euclidean",#annotation_row=anno_row,
show_colnames=F, show_rownames=T,file=file.path(PreDE.dir , "EpiPlayer","EpiPlayer_euclidean_WtUT_vs_Arid1aUT_DEG.pdf"),fontsize_row=4)




#subset sign genes
vst_sub_sign <- vst_damage[which(rownames(vst_damage) %in% c(DEG_results_list_sign$WT_Tam.untreated_Irr.treated_Hdac.untreated_vs_ID3_Tam.untreated_Irr.treated_Hdac.untreated$ensembl )),]
dim(vst_sub_sign)
vst_sub_sign <- as.data.frame(vst_sub_sign)

#heatmap annotation
#columns
anno_vst <- anno_sub[,c("genotype", "irradiation_treatment")]
#rows
rownames(vst_sub_sign) <- vst_sub_sign$symbol
#anno_row <- vst_sub_sign[,c("Pathway"), drop=F]
#anno_names<- ifelse(is.na(anno_names),vst_sub_sign$ensembl, anno_names )

#plot heatmap
dir.create(file.path(PreDE.dir , "DNA_damage"))
pheatmap(vst_sub_sign[, rownames(anno_sub)], annotation_col=as.data.frame(anno_vst) ,scale="row", clustering_distance_cols="correlation",#annotation_row=anno_row,
show_colnames=F, show_rownames=T,file=file.path(PreDE.dir , "DNA_damage","DNAdamageGenes_correlation_WtIrr_ID3Irr_ID3rescueIrr_compWTIrr_vsID3Irr_onlypadj.pdf"),fontsize_row=8)
pheatmap(vst_sub_sign[, rownames(anno_sub)], annotation_col=as.data.frame(anno_vst) ,scale="row", clustering_distance_cols="euclidean",#annotation_row=anno_row,
show_colnames=F, show_rownames=T,file=file.path(PreDE.dir , "DNA_damage","DNAdamageGenes_euclidean_WtIrr_ID3Irr_ID3rescueIrr_compWTIrr_vsID3Irr_onlypadj.pdf"),fontsize_row=8)



