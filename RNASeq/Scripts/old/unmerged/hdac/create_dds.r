#subset hddac sample set
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)

#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
#dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
#dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
#dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
#dir.create(PostDE.dir)

#load data
counts<- readRDS(file.path(data.dir, "gene_counts.rds"))

#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)

#get necessary samples
#hdac_pos <- unique(c(anno[anno$hdac_treatment=="treated",]$sample_name, 
#anno[anno$genotype=="ID3" & anno$irr_treatment=="untreated",]$sample_name, 
#anno[anno$genotype=="WT" & anno$irr_treatment=="untreated",]$sample_name))

hdac_pos <- unique(c(anno[anno$hdac_treatment=="treated",]$sample_name))

#subset data
counts_hdac <- counts[,hdac_pos]
anno_hdac <- anno[hdac_pos,]
saveRDS(counts_hdac, file.path(results.dir, "counts_hdac.rds"))
saveRDS(anno_hdac, file.path(results.dir, "anno_hdac.rds"))

#create dds
#dds_hdac <- DESeqDataSetFromMatrix(countData = counts_hdac, 
#                              colData = anno_hdac, 
#                              design = ~ replicate + genotype + hdac_treatment + dox_treatment)
#
dds_hdac <- DESeqDataSetFromMatrix(countData = counts_hdac, 
                              colData = anno_hdac, 
                              design = ~ replicate + genotype  + dox_treatment)

#Estimating size factors
#for Genotype comparison
dds_hdac <- estimateSizeFactors(dds_hdac)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds_hdac, normalized=TRUE) >= 1 ) >= 1
dim(dds_hdac[idx,])
dds_hdac <- dds_hdac[idx,]
dim(dds_hdac)

#Running the differential expression 
dds_hdac <- DESeq(dds_hdac)
saveRDS(dds_hdac,file = file.path(results.dir,"dds_hdac.rds"))
