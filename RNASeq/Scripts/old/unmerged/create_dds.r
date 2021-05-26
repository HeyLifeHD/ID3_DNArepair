#libraries
library(DESeq2)
#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)


#load data
counts<- readRDS(file.path(data.dir, "gene_counts.rds"))
anno <- readRDS(file.path(data.dir, "anno.rds"))

anno$irr_treatment <- "untreated"
anno[anno$dox_treatment =="irradiated",]$irr_treatment <- "treated"
anno$dox_treatment <- as.character(anno$dox_treatment)
anno[anno$dox_treatment =="irradiated",]$dox_treatment <- "untreated"


#create dds
#dds <- DESeqDataSetFromMatrix(countData = counts, 
#                              colData = anno, 
#                              design = ~ replicate + genotype + hdac_treatment + dox_treatment + genotype:hdac_treatment + genotype:dox_treatment)
#
#model.matrix(~ replicate + genotype + hdac_treatment + dox_treatment + genotype:hdac_treatment + genotype:dox_treatment, 
#anno[,c("replicate", "genotype","hdac_treatment","dox_treatment"  )])

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = anno, 
                              design = ~ replicate + genotype + hdac_treatment + dox_treatment +irr_treatment)

#Estimating size factors
#for Genotype comparison
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
length(idx)
dds <- dds[idx,]
dim(dds)

#Plot dispersion
pdf(file.path(PreDE.dir, "Dispersion(>1).pdf"))
plotDispEsts(dds, main="Dispersion plot")
dev.off()

#Running the differential expression 
dds <- DESeq(dds)
saveRDS(dds,file = file.path(results.dir,"dds.rds"))
