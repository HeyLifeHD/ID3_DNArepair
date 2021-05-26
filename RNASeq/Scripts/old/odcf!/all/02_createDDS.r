#libraries
library(DESeq2)
#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
counts<- readRDS(file.path(data.dir, "counts.rds"))
#counts_df <- as.data.frame(counts)
#rownames(counts_df)<- rownames(counts)
#counts <- counts_df
#saveRDS(counts,file.path(data.dir, "counts.rds"))
anno <- readRDS(file.path(data.dir, "anno.rds"))
anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment)

#create dds
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = anno, 
                              design = ~ replicate + genotype + tamoxifen_treatment + irradiation_treatment +hdac_treatment )

#Estimating size factors
#for Genotype comparison
dds <- estimateSizeFactors(dds)

#Filter genes which are only expressed in 1 sample
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
table(idx)
dds <- dds[idx,]
dim(dds)

#Running the DESEQ
dds <- DESeq(dds)
#Plot dispersion
pdf(file.path(PreDE.dir, "Dispersion(1).pdf"))
DESeq2::plotDispEsts(dds, main="Dispersion plot")
dev.off()
#save object
saveRDS(dds,file = file.path(results.dir,"dds_all.rds"))
