#libraries
library(DESeq2)
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/"
data.dir <- file.path(base.dir, "data")
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/test_wo4"
dir.create(base.dir)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)

#load data
counts<- readRDS(file.path(data.dir , "counts.rds"))
sample_anno <- readRDS(file.path(data.dir, "sample_anno.rds"))
#subset
sample_anno_sub <- sample_anno[sample_anno$genotype %in% c("WT", "ARID1A", "ID3"), ]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$tamoxifen_treatment =="UNTREATED",]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$hdac_treatment. =="UNTREATED",]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$Plattform =="Illumina NovaSeq 6000 S1",]
sample_anno_sub <- sample_anno_sub[!sample_anno_sub$Replicate %in% c(8),]

counts_sub <- counts[,rownames(sample_anno_sub)]
sample_anno_sub$Replicate <-as.character(sample_anno_sub$Replicate )

#create dds
dds <- DESeqDataSetFromMatrix(countData = counts_sub, 
                              colData = sample_anno_sub, 
                              design = ~ Replicate+group )

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
saveRDS(dds,file = file.path(results.dir,"dds.rds"))
