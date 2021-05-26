#Directories
#NGS-machine
.libPaths("/home/c010-ngs/R/x86_64-pc-linux-gnu-library/3.4_alt")
base.dir <- "/home/c010-ngs/c010-datasetshey/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")

#VM
base.dir <- "/C010-datasets/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")

#PC-09
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")
temp.dir <- "/home/heyj/rTemp"


#PC-09
base.dir <- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")
temp.dir <- "/home/epicwl/rTemp"

#library
library(data.table)
library(matrixStats)
#load counts
counts <- read.table("c010-datasets/Internal/2018-Ali/RNASeq/181114_analysis_rep1/repeat_counts.txt")

sample_anno <- data.frame(Sample_ID=as.character(c("WT_n4OHT","ARID1A-KO_n4OHT","ID3-KO_n4OHT","WTsiMDC1_n4OHT","WT_4OHT","ARID1A-KO_4OHT","ID3-KO_4OHT","WTsiMDC1_4OHT")))
sample_anno$Sample_ID<- paste0(sample_anno$Sample_ID, "_rep1")
sample_anno$Sample_ID <- as.character(sample_anno$Sample_ID)
sample_anno$Replicate <- 1
sample_anno$Treatment <- sapply(strsplit(sample_anno$Sample_ID ,"_", fixed=TRUE),`[`, 2)
sample_anno$Genotype <- sapply(strsplit(sample_anno$Sample_ID ,"_", fixed=TRUE),`[`, 1)

colnames(counts)<- c("repeat_id", "seqnames", "start", "end","strand","length",sample_anno$Sample_ID)

#Load repeat regions and subset LTR12 repeats
repeats=readRDS(c("/bigdisk/Nanopore/raw_data/repeats/repeats.RDS"))

#normalize
counts_norm <- counts
counts_norm[,  (ncol(counts_norm)-8):ncol(counts_norm)]<- sapply(counts_norm[,  (ncol(counts_norm)-8):ncol(counts_norm)],function(x) as.numeric(x))
#calculate cpm
counts_norm[,  (ncol(counts_norm)-8):ncol(counts_norm)]<- counts_norm[,  (ncol(counts_norm)-8):ncol(counts_norm)]/colSums(counts_norm[,  (ncol(counts_norm)-8):ncol(counts_norm)])
counts_norm[,  (ncol(counts_norm)-8):ncol(counts_norm)] <- counts_norm[,  (ncol(counts_norm)-8):ncol(counts_norm)]*10^6
#combine 
counts_norm<- cbind(repeats,counts_norm[-1,  (ncol(counts_norm)-8):ncol(counts_norm)])

#as data.table
#setDT(counts_norm)
#summarize repeat
#counts_melt_repFamily <- data.table::melt(counts_norm[,c(samples, 'repFamily'), with = FALSE])[,sum(value), by=,.(repFamily)]
#sum_counts_repFamily <- (data.table::dcast(counts_melt_repFamily, variable ~ repFamily))

counts_sum_repFamily <- aggregate(cbind(counts_norm[,(27-8):27]), by=list(repFamily=counts_norm$repFamily), FUN=sum)
counts_sum_repClass <- aggregate(cbind(counts_norm[,(27-8):27]), by=list(repClass=counts_norm$repClass), FUN=sum)
counts_sum_repName <- aggregate(cbind(counts_norm[,(27-8):27]), by=list(repName=counts_norm$repName), FUN=sum)

#subset ltr
counts_sum_repFamily_erv <- counts_sum_repFamily[grep("ERV", counts_sum_repFamily$repFamily),]
counts_sum_repClass_ltr <- counts_sum_repClass[grep("LTR", counts_sum_repClass$repClass),]
counts_sum_repName_ltr12 <- counts_sum_repName[grep("LTR12", counts_sum_repName$repName),]

#Plot LTR12
pdf("c010-datasets/Internal/2018-Ali/RNASeq/181114_analysis_rep1/repName_Heatmap_log2_scaled.pdf")
pheatmap(log2(counts_sum_repName_ltr12[,sample_anno$Sample_ID]), 
                scale="row",
                #annotation_col=anno, 
                labels_row=counts_sum_repName_ltr12$repName,
                show_colnames = TRUE,show_rownames = TRUE,
                #annotation_colors = anno_colors,
                cluster_rows=F,fontsize = 5,
                #color=col_palette,
                border_color = "grey"
                #cellwidth=5,
)
dev.off()

#get most variant rep
topVar<- head(order(rowVars(as.matrix(counts_sum_repName[, sample_anno$Sample_ID])), decreasing=TRUE),100)

pdf("c010-datasets/Internal/2018-Ali/RNASeq/181114_analysis_rep1/all_Heatmap.pdf")
pheatmap(counts_sum_repName[topVar,sample_anno$Sample_ID], 
                scale="none",
                #annotation_col=anno, 
                labels_row=counts_sum_repName[topVar,"repName"],
                show_colnames = TRUE,show_rownames = TRUE,
                #annotation_colors = anno_colors,
                cluster_rows=F,fontsize = 5,
                #color=col_palette,
                border_color = "grey"
                #cellwidth=5,
)
dev.off()

#all no subset
pdf("c010-datasets/Internal/2018-Ali/RNASeq/181114_analysis_rep1/all_nosub_Heatmap_scaled.pdf")
pheatmap(counts_sum_repName[,sample_anno$Sample_ID], 
                scale="row",
                #annotation_col=anno, 
                labels_row=counts_sum_repName[,"repName"],
                show_colnames = TRUE,show_rownames = TRUE,
                #annotation_colors = anno_colors,
                cluster_rows=T,fontsize = 5,
                #color=col_palette,
                border_color = "grey"
                #cellwidth=5,
)
dev.off()