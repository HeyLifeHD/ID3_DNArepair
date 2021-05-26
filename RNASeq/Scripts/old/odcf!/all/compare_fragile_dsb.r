
#libraries
library(DESeq2)
library(ggpubr)
library(pheatmap)
library(org.Hs.eg.db)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(reshape2)
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
data.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/external_data/"
#Read in Data
dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
anno <- colData(dds_genes)
#expression
vst_genes<- readRDS(file = file.path(results.dir,"vst_genes_batchRem_rep.rds"))
mat_genes<- assay(vst_genes)
rownames(mat_genes)<- sapply(strsplit(rownames(mat_genes), ".", fixed=TRUE), "[", 1)
#sites of interest
fragile_dsb_hg18 <-import.bed(file.path(data.dir,"GSE93038_Dellino_BLISS_Processed_Data_Tier1_DSBs.bed.txt"))

#lift over
ch <- import.chain(file.path(data.dir, "hg18ToHg19.over.chain"))
seqlevelsStyle(fragile_dsb_hg18) = "UCSC"
fragile_dsb_hg19<-unlist(liftOver(fragile_dsb_hg18, ch))

#annotate those regions
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
fragile_dsb_hg19_anno <- annotatePeak(peak= fragile_dsb_hg19, tssRegion=c(-1500, 1500),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
fragile_dsb_hg19_anno_df <- as.data.frame(fragile_dsb_hg19_anno)

#start by looking at promoters
promoter_dsb <- fragile_dsb_hg19_anno_df[grep("Promoter", fragile_dsb_hg19_anno_df$annotation),]$ENSEMBL

#subset interesting samples: ID3
anno_sub <- anno[which(anno$genotype %in% c("WT", "ID3") & anno$replicate != 0 &anno$irradiation_treatment=="untreated" & anno$hdac_treatment=="untreated" & anno$tamoxifen_treatment =="untreated"),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% promoter_dsb,]
mat_genes_sub_fil <- mat_genes_sub[apply(mat_genes_sub, MARGIN = 1, FUN = function(x) sd(x) != 0),]

length(promoter_dsb)
dim(mat_genes_sub_fil)

#plot heatmap of those
#anno_vis <- anno_sub[, c("replicate","genotype","irradiation_treatment", "tamoxifen_treatment", "hdac_treatment")]
anno_vis <- anno_sub[, c("replicate","genotype")]
#anno_vis$irradiation_treatment <- droplevels(as.factor(anno_vis$irradiation_treatment))
#anno_vis$hdac_treatment <- droplevels(as.factor(anno_vis$hdac_treatment))
anno_vis <- as.data.frame(anno_vis)
pheatmap(mat_genes_sub_fil, annotation_col=anno_vis, scale="row",show_rownames=FALSE, show_colnames=FALSE,
 file=file.path(base_results.dir, "fragile_dsb", "promoter_dsb_all_id3_vs_wt_untreated.pdf"))


#subset interesting samples: Arid1a
anno_sub <- anno[which(anno$genotype %in% c("WT", "ARID1A") & anno$replicate != 0 & anno$irradiation_treatment=="untreated" & anno$hdac_treatment=="untreated" & anno$tamoxifen_treatment =="untreated"),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% promoter_dsb,]
mat_genes_sub_fil <- mat_genes_sub[apply(mat_genes_sub, MARGIN = 1, FUN = function(x) sd(x) != 0),]

length(promoter_dsb)
dim(mat_genes_sub_fil)

#plot heatmap of those
#anno_vis <- anno_sub[, c("replicate","genotype","irradiation_treatment", "tamoxifen_treatment", "hdac_treatment")]
anno_vis <- anno_sub[, c("replicate","genotype")]
#anno_vis$irradiation_treatment <- droplevels(as.factor(anno_vis$irradiation_treatment))
#anno_vis$hdac_treatment <- droplevels(as.factor(anno_vis$hdac_treatment))
anno_vis <- as.data.frame(anno_vis)
pheatmap(mat_genes_sub_fil, annotation_col=anno_vis, scale="row",show_rownames=FALSE, show_colnames=FALSE,
 file=file.path(base_results.dir, "fragile_dsb", "promoter_dsb_all_arid1a_vs_wt_untreated.pdf"))


#plot all
anno_sub <- anno[which(anno$replicate != 0 ),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% promoter_dsb,]
mat_genes_sub_fil <- mat_genes_sub[apply(mat_genes_sub, MARGIN = 1, FUN = function(x) sd(x) != 0),]

length(promoter_dsb)
dim(mat_genes_sub_fil)

#plot heatmap of those
anno_vis <- anno_sub[, c("replicate","genotype","irradiation_treatment", "tamoxifen_treatment", "hdac_treatment")]
#anno_vis <- anno_sub[, c("replicate","genotype")]
#anno_vis$irradiation_treatment <- droplevels(as.factor(anno_vis$irradiation_treatment))
#anno_vis$hdac_treatment <- droplevels(as.factor(anno_vis$hdac_treatment))
anno_vis <- as.data.frame(anno_vis)
pheatmap(mat_genes_sub_fil, annotation_col=anno_vis, scale="row",show_rownames=FALSE, show_colnames=FALSE,
 file=file.path(base_results.dir, "fragile_dsb", "promoter_dsb_all_allSamples.pdf"))











#start by looking at promoters and gene bodies
fragile_dsb_hg19_anno_df_sub <- fragile_dsb_hg19_anno_df[grep("Inter", fragile_dsb_hg19_anno_df$annotation, invert=TRUE),]
fragile_dsb_hg19_anno_df_sub <- fragile_dsb_hg19_anno_df_sub[grep("Down", fragile_dsb_hg19_anno_df_sub$annotation, invert=TRUE),]
DSBInterest <- fragile_dsb_hg19_anno_df_sub[grep("UTR", fragile_dsb_hg19_anno_df_sub$annotation, invert=TRUE),]$ENSEMBL

#subset interesting samples: ID3
anno_sub <- anno[which(anno$genotype %in% c("WT", "ID3") & anno$replicate != 0 &anno$irradiation_treatment=="untreated" & anno$hdac_treatment=="untreated" & anno$tamoxifen_treatment =="untreated"),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% DSBInterest,]
mat_genes_sub_fil <- mat_genes_sub[apply(mat_genes_sub, MARGIN = 1, FUN = function(x) sd(x) != 0),]

length(DSBInterest)
dim(mat_genes_sub_fil)

#plot heatmap of those
#anno_vis <- anno_sub[, c("replicate","genotype","irradiation_treatment", "tamoxifen_treatment", "hdac_treatment")]
anno_vis <- anno_sub[, c("replicate","genotype")]
#anno_vis$irradiation_treatment <- droplevels(as.factor(anno_vis$irradiation_treatment))
#anno_vis$hdac_treatment <- droplevels(as.factor(anno_vis$hdac_treatment))
anno_vis <- as.data.frame(anno_vis)
pheatmap(mat_genes_sub_fil, annotation_col=anno_vis, scale="row",show_rownames=FALSE, show_colnames=FALSE,
 file=file.path(base_results.dir, "fragile_dsb", "PromoterGenebody_dsb_all_id3_vs_wt_untreated.pdf"))

#get average expression for genotypes
mat_genes_sub_genoSum <- data.frame(ID3= rowMeans(mat_genes_sub[,anno_sub$genotype=="ID3"]),WT= rowMeans(mat_genes_sub[,anno_sub$genotype=="WT"]))

mat_genes_sub_genoSum_1 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) < quantile(rowMeans(mat_genes_sub_genoSum))[2]), ]
mat_genes_sub_genoSum_1$cluster <- "cluster1"
mat_genes_sub_genoSum_2 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) > quantile(rowMeans(mat_genes_sub_genoSum))[2] & rowMeans(mat_genes_sub_genoSum) < quantile(rowMeans(mat_genes_sub_genoSum))[3]), ]
mat_genes_sub_genoSum_2$cluster <- "cluster2"
mat_genes_sub_genoSum_3 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) > quantile(rowMeans(mat_genes_sub_genoSum))[3] & rowMeans(mat_genes_sub_genoSum) < quantile(rowMeans(mat_genes_sub_genoSum))[4]), ]
mat_genes_sub_genoSum_3$cluster <- "cluster3"
mat_genes_sub_genoSum_4 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) > quantile(rowMeans(mat_genes_sub_genoSum))[4]), ]
mat_genes_sub_genoSum_4$cluster <- "cluster4"
mat_genes_sub_genoSum <- rbind(mat_genes_sub_genoSum_1, mat_genes_sub_genoSum_2, mat_genes_sub_genoSum_3, mat_genes_sub_genoSum_4)

mat_genes_sub_genoSum<- melt(mat_genes_sub_genoSum)
pdf(file.path(base_results.dir, "fragile_dsb", "PromoterGenebody_dsb_all_arid1a_vs_wt_untreated_averageExpr.pdf"))
ggboxplot(mat_genes_sub_genoSum, x="variable", "value", col="variable", facet.by="cluster")+ stat_compare_means(method = "wilcox.test")
dev.off()

counts <- counts(dds_genes, normalized=TRUE)
rownames(counts)<- sapply(strsplit(rownames(counts), ".", fixed=TRUE), "[", 1)
anno_sub <- anno[which(anno$genotype %in% c("WT", "ID3") & anno$replicate != 0 &anno$irradiation_treatment=="untreated" & anno$hdac_treatment=="untreated" & anno$tamoxifen_treatment =="untreated"),]
counts_sub <- counts[, rownames(anno_sub)]
counts_sub <- counts_sub[rownames(counts_sub) %in% DSBInterest,]
counts_sub_genoSum <- data.frame(ID3= rowMeans(counts_sub[,anno_sub$genotype=="ID3"]),WT= rowMeans(counts_sub[,anno_sub$genotype=="WT"]))

counts_sub_genoSum_1 <- counts_sub_genoSum[which(rowMeans(counts_sub_genoSum) < quantile(rowMeans(counts_sub_genoSum))[2]), ]
counts_sub_genoSum_1$cluster <- "cluster1"
counts_sub_genoSum_2 <- counts_sub_genoSum[which(rowMeans(counts_sub_genoSum) > quantile(rowMeans(counts_sub_genoSum))[2] & rowMeans(counts_sub_genoSum) < quantile(rowMeans(counts_sub_genoSum))[3]), ]
counts_sub_genoSum_2$cluster <- "cluster2"
counts_sub_genoSum_3 <- counts_sub_genoSum[which(rowMeans(counts_sub_genoSum) > quantile(rowMeans(counts_sub_genoSum))[3] & rowMeans(counts_sub_genoSum) < quantile(rowMeans(counts_sub_genoSum))[4]), ]
counts_sub_genoSum_3$cluster <- "cluster3"
counts_sub_genoSum_4 <- counts_sub_genoSum[which(rowMeans(counts_sub_genoSum) > quantile(rowMeans(counts_sub_genoSum))[4]), ]
counts_sub_genoSum_4$cluster <- "cluster4"

counts_sub_genoSum <- rbind(counts_sub_genoSum_1, counts_sub_genoSum_2, counts_sub_genoSum_3, counts_sub_genoSum_4)

counts_sub_genoSum<- melt(counts_sub_genoSum)
pdf(file.path(base_results.dir, "fragile_dsb", "PromoterGenebody_dsb_all_arid1a_vs_wt_untreated_averageExprCounts.pdf"))
ggboxplot(counts_sub_genoSum, x="variable", "value", facet.by="cluster")
dev.off()

#subset interesting samples: Arid1a
anno_sub <- anno[which(anno$genotype %in% c("WT", "ARID1A") & anno$replicate != 0 & anno$irradiation_treatment=="untreated" & anno$hdac_treatment=="untreated" & anno$tamoxifen_treatment =="untreated"),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% promoter_dsb,]
mat_genes_sub_fil <- mat_genes_sub[apply(mat_genes_sub, MARGIN = 1, FUN = function(x) sd(x) != 0),]

length(DSBInterest)
dim(mat_genes_sub_fil)

#plot heatmap of those
#anno_vis <- anno_sub[, c("replicate","genotype","irradiation_treatment", "tamoxifen_treatment", "hdac_treatment")]
anno_vis <- anno_sub[, c("replicate","genotype")]
#anno_vis$irradiation_treatment <- droplevels(as.factor(anno_vis$irradiation_treatment))
#anno_vis$hdac_treatment <- droplevels(as.factor(anno_vis$hdac_treatment))
anno_vis <- as.data.frame(anno_vis)
pheatmap(mat_genes_sub_fil, annotation_col=anno_vis, scale="row",show_rownames=FALSE, show_colnames=FALSE,
 file=file.path(base_results.dir, "fragile_dsb", "PromoterGenebody_dsb_all_arid1a_vs_wt_untreated.pdf"))


#plot all
anno_sub <- anno[which(anno$replicate != 0 ),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% promoter_dsb,]
mat_genes_sub_fil <- mat_genes_sub[apply(mat_genes_sub, MARGIN = 1, FUN = function(x) sd(x) != 0),]

length(DSBInterest)
dim(mat_genes_sub_fil)

#plot heatmap of those
anno_vis <- anno_sub[, c("replicate","genotype","irradiation_treatment", "tamoxifen_treatment", "hdac_treatment")]
#anno_vis <- anno_sub[, c("replicate","genotype")]
#anno_vis$irradiation_treatment <- droplevels(as.factor(anno_vis$irradiation_treatment))
#anno_vis$hdac_treatment <- droplevels(as.factor(anno_vis$hdac_treatment))
anno_vis <- as.data.frame(anno_vis)
pheatmap(mat_genes_sub_fil, annotation_col=anno_vis, scale="row",show_rownames=FALSE, show_colnames=FALSE,
 file=file.path(base_results.dir, "fragile_dsb", "PromoterGenebody_dsb_all_allSamples.pdf"))


#subset interesting samples: no hdac , irradiation and tamoxifen treatment
anno_sub <- anno[which(anno$replicate != 0 &anno$irradiation_treatment=="untreated" & anno$hdac_treatment=="untreated" & anno$tamoxifen_treatment =="untreated"),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% DSBInterest,]
mat_genes_sub_fil <- mat_genes_sub[apply(mat_genes_sub, MARGIN = 1, FUN = function(x) sd(x) != 0),]

length(DSBInterest)
dim(mat_genes_sub_fil)

#plot heatmap of those
#anno_vis <- anno_sub[, c("replicate","genotype","irradiation_treatment", "tamoxifen_treatment", "hdac_treatment")]
anno_vis <- anno_sub[, c("replicate","genotype")]
#anno_vis$irradiation_treatment <- droplevels(as.factor(anno_vis$irradiation_treatment))
#anno_vis$hdac_treatment <- droplevels(as.factor(anno_vis$hdac_treatment))
anno_vis <- as.data.frame(anno_vis)
pheatmap(mat_genes_sub_fil, annotation_col=anno_vis, scale="row",show_rownames=FALSE, show_colnames=FALSE,
 file=file.path(base_results.dir, "fragile_dsb", "PromoterGenebody_dsb_all_Alluntreated.pdf"))




#get average expression for genotypes
anno_sub <- anno[which(anno$replicate != 0 &anno$irradiation_treatment=="untreated" & anno$hdac_treatment=="untreated" & anno$tamoxifen_treatment =="untreated"),]
mat_genes_sub <- mat_genes[, rownames(anno_sub)]
mat_genes_sub <- mat_genes_sub[rownames(mat_genes_sub) %in% DSBInterest,]
mat_genes_sub_genoSum <- data.frame(ID3= rowMeans(mat_genes_sub[,anno_sub$genotype=="ID3"]),WT= rowMeans(mat_genes_sub[,anno_sub$genotype=="WT"]), ARID1A= rowMeans(mat_genes_sub[,anno_sub$genotype=="ARID1A"]),
 ID3rescue= rowMeans(mat_genes_sub[,anno_sub$genotype=="ID3_rescue"]))

mat_genes_sub_genoSum_1 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) < quantile(rowMeans(mat_genes_sub_genoSum))[2]), ]
mat_genes_sub_genoSum_1$cluster <- "cluster1"
mat_genes_sub_genoSum_2 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) > quantile(rowMeans(mat_genes_sub_genoSum))[2] & rowMeans(mat_genes_sub_genoSum) < quantile(rowMeans(mat_genes_sub_genoSum))[3]), ]
mat_genes_sub_genoSum_2$cluster <- "cluster2"
mat_genes_sub_genoSum_3 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) > quantile(rowMeans(mat_genes_sub_genoSum))[3] & rowMeans(mat_genes_sub_genoSum) < quantile(rowMeans(mat_genes_sub_genoSum))[4]), ]
mat_genes_sub_genoSum_3$cluster <- "cluster3"
mat_genes_sub_genoSum_4 <- mat_genes_sub_genoSum[which(rowMeans(mat_genes_sub_genoSum) > quantile(rowMeans(mat_genes_sub_genoSum))[4]), ]
mat_genes_sub_genoSum_4$cluster <- "cluster4"
mat_genes_sub_genoSum <- rbind(mat_genes_sub_genoSum_1, mat_genes_sub_genoSum_2, mat_genes_sub_genoSum_3, mat_genes_sub_genoSum_4)

mat_genes_sub_genoSum<- melt(mat_genes_sub_genoSum)
pdf(file.path(base_results.dir, "fragile_dsb", "PromoterGenebody_dsb_all_allGeno_untreated_averageExpr.pdf"))
ggboxplot(mat_genes_sub_genoSum, x="variable", "value", col="variable", facet.by="cluster")+ stat_compare_means(method = "anova")
dev.off()
compare_means(value~variable,data=mat_genes_sub_genoSum, method="t.test")