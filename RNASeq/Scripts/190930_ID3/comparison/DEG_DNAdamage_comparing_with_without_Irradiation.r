library(DESeq2)
library(openxlsx)
library(dplyr)
library(ggpubr)
library("rafalib")
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(randomcoloR)
library(RColorBrewer)

#Epigenetic Players 
#RNA analysis
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
DEG_results_list_Irr<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
rld_b_Irr<- readRDS(file.path(results.dir,"rld_b_replicate.rds"))

#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
DEG_results_list_Untr<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
rld_b_Untr<- readRDS(file.path(results.dir,"rld_b_replicate.rds"))


#load dna damage lists
#load term2 gene list of damage
t2g <- read.csv2(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190531_odcf_NoIll/","data", "updated_DNA_damage_player.csv"))
t2g$ENSEMBL<-as.character(t2g$ensembl) 
t2g$Pathway<-as.character(t2g$Pathway) 
rownames(t2g)<- t2g$ensembl
t2g$symbol<- mapIds(org.Hs.eg.db, keys= t2g$ENSEMBL , keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )

#subset list
lapply(DEG_results_list_Irr,function(x)dim(x))
DEG_results_list_Irr_sub <- lapply(DEG_results_list_Irr, function(x){
    x <- x[which(rownames(x) %in% t2g$ENSEMBL),]
    x
})
lapply(DEG_results_list_Irr_sub,function(x)dim(x))

for(i in names(DEG_results_list_Irr_sub)){
    write.table(DEG_results_list_Irr_sub[[i]], file.path(PostDE.dir, i , paste0("all_DNAdamage_sub.tsv")), sep="\t", row.names=TRUE, quote=FALSE)
}

DEG_results_list_Irr_sub_sig <- lapply(DEG_results_list_Irr_sub, function(x){
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange)>0.5),]
    x
})
lapply(DEG_results_list_Irr_sub_sig,function(x)dim(x))

DEG_results_list_Irr_sub_sig_onlyFDR <- lapply(DEG_results_list_Irr_sub, function(x){
    x <- x[which(x$padj < 0.05) ,]
    x
})
lapply(DEG_results_list_Irr_sub_sig_onlyFDR,function(x)dim(x))

#merge list which will be plotted
#for irr data
pheno_irr <- colData(rld_b_Irr)
pheno_irr_sub <- pheno_irr[which(pheno_irr$genotype %in% c("ID3", "WT")),]
rld_b_Irr_sub <- rld_b_Irr[, rownames(pheno_irr_sub) ]
#for untr data
pheno_Untr <- colData(rld_b_Untr)
pheno_Untr_sub <- pheno_Untr[which(pheno_Untr$genotype %in% c("ID3", "WT", "ID3_rescue")),]
rld_b_Untr_sub <- rld_b_Untr[, rownames(pheno_Untr_sub) ]
#combine
rld_b_Irr_sub_as  <- as.data.frame(assay(rld_b_Irr_sub))
rld_b_Irr_sub_as$ensembl <- rownames(rld_b_Irr_sub_as)
rld_b_Untr_sub_as <- as.data.frame(assay(rld_b_Untr_sub))
rld_b_Untr_sub_as$ensembl <- rownames(rld_b_Untr_sub_as)
merged_rld_b <- left_join(rld_b_Irr_sub_as, rld_b_Untr_sub_as)
rownames(merged_rld_b)<- merged_rld_b$ensembl
merged_rld_b$ensembl <- NULL

pheno <- rbind(pheno_irr, pheno_Untr)
pheno_sub <- pheno[rownames(pheno) %in% colnames(merged_rld_b),]
#subset genes of interest
merged_rld_b_sub <- merged_rld_b[rownames(merged_rld_b) %in% rownames(DEG_results_list_Irr_sub_sig$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED),]

#plot data
anno_vst <- pheno_sub[,c("genotype","Replicate","irradiation_treatment")]
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(pheno_sub$genotype)), name = 'Dark2')
names(color_genotype)<- as.character(unique(pheno_sub$genotype))
color_replicate<-randomColor(length(unique(pheno_sub$Replicate)))
names(color_replicate)<- as.character(unique(pheno_sub$Replicate))
color_irradiation<-randomColor(length(unique(pheno_sub$irradiation_treatment)))
names(color_irradiation)<- as.character(unique(pheno_sub$irradiation_treatment))
anno_colors <- list(genotype=color_genotype,  Replicate=color_replicate, irradiation_treatment=color_irradiation)

heat<- pheatmap(merged_rld_b_sub,
                        scale="row", show_colnames=F,
                       annotation_col=as.data.frame(anno_vst),fontsize_row=3,
                       annotation_colors=anno_colors, clustering_distance_rows="correlation",
                       show_rownames=F) 
pdf(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/","Heat_sig_IrrComp_DNAdamage.pdf"))
print(heat)
dev.off() 

#seperated heatmap
#irradiation
anno_vst <- pheno_irr_sub[,c("genotype","Replicate","irradiation_treatment")]
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(pheno_irr_sub$genotype)), name = 'Dark2')
names(color_genotype)<- as.character(unique(pheno_irr_sub$genotype))
color_replicate<-randomColor(length(unique(pheno_irr_sub$Replicate)))
names(color_replicate)<- as.character(unique(pheno_irr_sub$Replicate))
color_irradiation<-randomColor(length(unique(pheno_irr_sub$irradiation_treatment)))
names(color_irradiation)<- as.character(unique(pheno_irr_sub$irradiation_treatment))
anno_colors <- list(genotype=color_genotype,  Replicate=color_replicate, irradiation_treatment=color_irradiation)

heat<- pheatmap(merged_rld_b_sub[,rownames(anno_vst)],
                        scale="row", show_colnames=F,
                       annotation_col=as.data.frame(anno_vst),fontsize_row=3,
                       annotation_colors=anno_colors, clustering_distance_rows="correlation",
                       show_rownames=T) 

pdf(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/","Heat_sig_IrrComp_DNAdamage_irradiation_only.pdf"))
print(heat)
dev.off() 

#non irradiated
anno_vst <- pheno_Untr_sub[,c("genotype","Replicate","irradiation_treatment")]
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(pheno_Untr_sub$genotype)), name = 'Dark2')
names(color_genotype)<- as.character(unique(pheno_Untr_sub$genotype))
color_replicate<-randomColor(length(unique(pheno_Untr_sub$Replicate)))
names(color_replicate)<- as.character(unique(pheno_Untr_sub$Replicate))
color_irradiation<-randomColor(length(unique(pheno_Untr_sub$irradiation_treatment)))
names(color_irradiation)<- as.character(unique(pheno_Untr_sub$irradiation_treatment))
anno_colors <- list(genotype=color_genotype,  Replicate=color_replicate, irradiation_treatment=color_irradiation)


temp <- merged_rld_b_sub[heat$tree_row$order,rownames(anno_vst)]
temp2 <- DEG_results_list_Irr_sub_sig$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED
temp2 <- temp2[rownames(temp),]
rownames(temp)<- temp2$symbol
heat2<- pheatmap(temp,
                        scale="row", show_colnames=F,cluster_rows=FALSE,
                       annotation_col=as.data.frame(anno_vst),#fontsize_row=3,
                       annotation_colors=anno_colors, clustering_distance_rows="correlation",
                       filename="/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/Heat_sig_IrrComp_DNAdamage_Untreated_only.pdf",
                       show_rownames=T) 
#second version
anno_vst <- pheno_Untr_sub[,c("genotype"),drop=F]
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(pheno_Untr_sub$genotype)), name = 'Dark2')
names(color_genotype)<- as.character(unique(pheno_Untr_sub$genotype))
anno_colors <- list(genotype=color_genotype)
temp <- merged_rld_b_sub[heat$tree_row$order,rownames(anno_vst)]
temp2 <- DEG_results_list_Irr_sub_sig$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED
temp2 <- temp2[rownames(temp),]
rownames(temp)<- temp2$symbol
heat2<- pheatmap(temp,
                        scale="row", show_colnames=F,cluster_rows=FALSE,
                       annotation_col=as.data.frame(anno_vst),#fontsize_row=3,
                       annotation_colors=anno_colors, clustering_distance_rows="correlation",
                       filename="/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/Heat_sig_IrrComp_DNAdamage_Untreated_only2.pdf",
                       show_rownames=T) 

#plot log2fold changes in an heatmap
DEG_sub_treated <- DEG_results_list_Irr_sub_sig$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED
DEG_sub_Untrtreated <- DEG_results_list_Untr$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED
l2fc <- DEG_sub_treated[,c("ensembl","log2FoldChange" )]
colnames(l2fc) <- c("ensembl","log2FoldChange_Irr" )
l2fc <- left_join(l2fc, DEG_sub_Untrtreated[,c("ensembl","log2FoldChange")])
rownames(l2fc)<- l2fc$ensembl
l2fc$ensembl <- NULL
colnames(l2fc)<- c("log2FoldChange_Irr", "log2FoldChange_Untr")
l2fc <- l2fc[order(l2fc$log2FoldChange_Irr, decreasing=TRUE),]
heat_lf<- pheatmap(l2fc,color=RColorBrewer::brewer.pal(n=10, "RdBu"),
                scale="none", show_colnames=T,cluster_rows=F,cluster_cols=F,
                clustering_distance_rows="correlation",
                show_rownames=T) 

pdf(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/","Heat_sig_IrrComp_DNAdamage_l2fcomp.pdf"), width=4)
print(heat_lf)
dev.off() 
