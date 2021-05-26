#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
#folder

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/repeats"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
featureCounts.dir <-"/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/results/repeats/"

#Data Read-in
#DESeq2 Analysis
#dds <- readRDS(file = file.path("/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/results/tables","dds_all.rds"))
anno <- readRDS(file = file.path("/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/data","anno.rds"))
rownames(anno)<- anno$sample_name
anno$group <- as.factor(paste0(anno$genotype,"_Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.",anno$hdac_treatment))
anno$replicate <- as.factor(anno$replicate)

#set index for analysis
#subset data
idx <- c("control-untreated_C010_DivA_WT_HDACi_R1", "control-untreated_C010_DivA_WT_HDACi_R2",
        "control-treated_C010_DivA_WT_R1", "control-treated_C010_DivA_WT_R2", "control-treated_C010_DivA_WT",
        "control-untreated_C010_DivA_WT_R1", "control-untreated_C010_DivA_WT_R2", "control-untreated_C010_DivA_WT",
        "knockout-treated_C010_DivA_ARID1A_R1", "knockout-treated_C010_DivA_ARID1A_R2", "knockout-treated_C010_DivA_ARID1A", 
        "knockout-untreated_C010_DivA_ARID1A_R1", "knockout-untreated_C010_DivA_ARID1A_R2", "knockout-untreated_C010_DivA_ARID1A",
        "control-irradiated_C010_DivA_WT_IR_R1", "control-irradiated_C010_DivA_WT_IR_R2")

#load repeat counts
#raw data for counts and repeats
repeats=readRDS(c("/bigdisk/Nanopore/raw_data/repeats/repeats.RDS"))
count_repeats<- readRDS(file.path(featureCounts.dir, "merged_repeat_counts.rds"))
#sum up per rep name
count_repeats<- sapply(count_repeats,function(x) as.integer(x))
#combine 
count_repeats<- cbind(as.data.frame(repeats),count_repeats)
#sum up
counts_sum_repFamily <- aggregate(cbind(count_repeats[,19:ncol(count_repeats)]), by=list(repFamily=count_repeats$repFamily), FUN=sum)
counts_sum_repClass <- aggregate(cbind(count_repeats[,19:ncol(count_repeats)]), by=list(repClass=count_repeats$repClass), FUN=sum)
counts_sum_repName <- aggregate(cbind(count_repeats[,19:ncol(count_repeats)]), by=list(repName=count_repeats$repName), FUN=sum)
repeat_list <- list(counts_sum_repFamily,counts_sum_repClass, counts_sum_repName )
names(repeat_list)<-  c("repFamily", "repClass", "repName")
saveRDS(repeat_list,file.path(results.dir, "repeat_list.rds") )
#create df
count_repeats <- counts_sum_repName[,rownames(anno[anno$sample_name %in% idx,])]
rownames(count_repeats)<- counts_sum_repName$repName

#create dds
#dds_repeats <- DESeqDataSetFromMatrix(countData = count_repeats, 
#                              colData = anno[anno$sample_name %in% idx,], 
#                              design = ~ replicate + sequencer + genotype + hdac_treatment + tamoxifen_treatment+irradiation_treatment)

dds_repeats <- DESeqDataSetFromMatrix(countData = count_repeats, 
                              colData = anno[anno$sample_name %in% idx,], 
                              design = ~ replicate + group)




#Estimating size factors
dds_repeats <- estimateSizeFactors(dds_repeats)

#Filter genes which are only expressed in 1 sample
length(dds_repeats)
idx_dds_repeats <- rowSums( counts(dds_repeats, normalized=TRUE) >= 1 ) >= 1
table(idx_dds_repeats)
#dds <- dds[idx,]
#dim(dds)

#Running the differential expression 
dds_repeats <- DESeq(dds_repeats)
saveRDS(dds_repeats, file.path(results.dir, "dds_repeats.rds"))















#counts genes
counts<- readRDS(file.path("/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/data", "counts.rds"))

#subset counts
counts_genes <- counts[,rownames(anno[anno$sample_name %in% idx,])]

#create dds sets
dds_genes <- DESeqDataSetFromMatrix(countData = counts_genes, 
                              colData = anno[anno$sample_name %in% idx,], 
                              design = ~ replicate + genotype + hdac_treatment + tamoxifen_treatment)

#Estimating size factors
dds_genes <- estimateSizeFactors(dds_genes)

#Filter genes which are only expressed in 1 sample
length(dds_genes)
idx_dds_genes <- rowSums( counts(dds_genes, normalized=TRUE) >= 1 ) >= 1
table(idx_dds_genes)
dds_genes <- dds_genes[idx_dds_genes,]
dim(dds_genes)

#Running the differential expression 
dds_genes <- DESeq(dds_genes)
saveRDS(dds_genes, file.path(results.dir, "dds_genes.rds"))
