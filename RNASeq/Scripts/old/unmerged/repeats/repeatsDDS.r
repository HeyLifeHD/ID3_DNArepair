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
repeats.dir <- file.path(base_results.dir, "repeats")
#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno_dds <- colData(dds)

#raw data for counts and repeats
repeats=readRDS(c("/bigdisk/Nanopore/raw_data/repeats/repeats.RDS"))
count_list<- readRDS("c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged/results/repeats/merged_counts_repeats.rds")

#sum up per rep name
counts<- count_list
counts<- sapply(counts,function(x) as.integer(x))

#combine 
counts<- cbind(as.data.frame(repeats),counts)

counts_sum_repFamily <- aggregate(cbind(counts[,19:ncol(counts)]), by=list(repFamily=counts$repFamily), FUN=sum)
saveRDS(counts_sum_repFamily, file.path(repeats.dir,"counts_sum_repFamily.rds" ))
counts_sum_repClass <- aggregate(cbind(counts[,19:ncol(counts)]), by=list(repClass=counts$repClass), FUN=sum)
saveRDS(counts_sum_repClass, file.path(repeats.dir,"counts_sum_repClass.rds" ))
counts_sum_repName <- aggregate(cbind(counts[,19:ncol(counts)]), by=list(repName=counts$repName), FUN=sum)
saveRDS(counts_sum_repName, file.path(repeats.dir,"counts_sum_repName.rds" ))



#create dds from counts_sum_repFamily
dds <- DESeqDataSetFromMatrix(countData = counts_sum_repName[,anno_dds$sample_name], 
                              colData = anno_dds, 
                              design = ~ replicate + genotype + hdac_treatment + dox_treatment +irr_treatment)

#Estimating size factors
#for Genotype comparison
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
length(idx)
#dds <- dds[idx,]
#dim(dds)

#Plot dispersion
#pdf(file.path(PreDE.dir, "Dispersion(>1).pdf"))
#plotDispEsts(dds, main="Dispersion plot")
#dev.off()

#Running the differential expression 
dds <- DESeq(dds)
saveRDS(dds,file = file.path(results.dir,"dds_repeats.rds"))
