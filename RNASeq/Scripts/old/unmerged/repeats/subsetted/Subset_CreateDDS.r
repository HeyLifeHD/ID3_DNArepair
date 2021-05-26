#subset hddac sample set
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(rtracklayer)

#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)

#load data
counts_genes<- readRDS(file.path(data.dir, "gene_counts.rds"))          

#load repeat counts
#raw data for counts and repeats
repeats=readRDS(c("/bigdisk/Nanopore/raw_data/repeats/repeats.RDS"))
count_repeats<- readRDS("c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged/results/repeats/merged_counts_repeats.rds")
#sum up per rep name
count_repeats<- sapply(count_repeats,function(x) as.integer(x))
#combine 
count_repeats<- cbind(as.data.frame(repeats),count_repeats)

counts_sum_repFamily <- aggregate(cbind(count_repeats[,19:ncol(count_repeats)]), by=list(repFamily=count_repeats$repFamily), FUN=sum)
counts_sum_repClass <- aggregate(cbind(count_repeats[,19:ncol(count_repeats)]), by=list(repClass=count_repeats$repClass), FUN=sum)
counts_sum_repName <- aggregate(cbind(count_repeats[,19:ncol(count_repeats)]), by=list(repName=count_repeats$repName), FUN=sum)

#subset data
idx <- c("ARID1A_knockout-treated", "ARID1A_R1_knockout-treated", "ARID1A_R2_knockout-treated", 
            "ARID1A_knockout-untreated", "ARID1A_R1_knockout-untreated", "ARID1A_R2_knockout-untreated",
            "WT_R1_control-treated", "WT_R2_control-treated", "WT_control-treated", 
            "WT_R1_control-untreated", "WT_R2_control-untreated", 
            "WT_HDACi_R1_control-untreated", "WT_HDACi_R2_control-untreated"
            )
#repeats
count_repeats <- counts_sum_repName[,rownames(anno[anno$group %in% idx,])]
rownames(count_repeats)<- counts_sum_repName$repName
#counts
#subset counts
counts_genes <- counts_genes[,rownames(anno[anno$group %in% idx,])]

#create dds sets
dds_genes <- DESeqDataSetFromMatrix(countData = counts_genes, 
                              colData = anno[anno$group %in% idx,], 
                              design = ~ replicate + genotype + hdac_treatment + dox_treatment)
dds_repeats <- DESeqDataSetFromMatrix(countData = count_repeats, 
                              colData = anno[anno$group %in% idx,], 
                              design = ~ replicate + genotype + hdac_treatment + dox_treatment)

#Estimating size factors
dds_genes <- estimateSizeFactors(dds_genes)
dds_repeats <- estimateSizeFactors(dds_repeats)

#Filter genes which are only expressed in 1 sample
length(dds_genes)
idx_dds_genes <- rowSums( counts(dds_genes, normalized=TRUE) >= 1 ) >= 1
length(idx_dds_genes)
length(dds_repeats)
idx_dds_repeats <- rowSums( counts(dds_repeats, normalized=TRUE) >= 1 ) >= 1
length(idx_dds_repeats)
#dds <- dds[idx,]
#dim(dds)

#Running the differential expression 
dds_genes <- DESeq(dds_genes)
saveRDS(dds_genes, file.path(base.dir, "results",  "repeats", "subsetted", "dds_genes.rds"))
dds_repeats <- DESeq(dds_repeats)
saveRDS(dds_repeats, file.path(base.dir, "results", "repeats", "subsetted", "dds_repeats.rds"))

