#libraries
library(DESeq2)
#folder
data_or.dir <- file.path("/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data", "data")

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190531_odcf_NoIll/"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

dir.create(base.dir)
dir.create(data.dir)
dir.create(base_results.dir)
dir.create(results.dir)
dir.create(PreDE.dir)
dir.create(PostDE.dir)

#load data
counts<- as.data.frame(readRDS(file.path(data_or.dir, "counts.rds")))
anno <- readRDS(file.path(data_or.dir, "anno.rds"))
counts <- counts[,rownames(anno[anno$replicate!=0,])]
saveRDS(counts,file.path(data.dir, "counts.rds"))
anno <- readRDS(file.path(data_or.dir, "anno.rds"))
anno<- anno[anno$replicate!=0,]
saveRDS(anno, file.path(data.dir, "anno.rds"))

anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment)

#create group annotation
anno$group <- as.factor(paste0(anno$genotype,"_Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.",anno$hdac_treatment))
anno$replicate <- as.factor(anno$replicate)
#create dds
dds_genes <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = anno, 
                              design = ~ replicate + group)
#Estimating size factors
dds_genes <- estimateSizeFactors(dds_genes)

#Filter genes which are only expressed in 1 sample
length(dds_genes)
idx_dds_genes <- rowSums( counts(dds_genes, normalized=TRUE) >= 1 ) >= 1
table(idx_dds_genes)
dds_genes <- dds_genes
dim(dds_genes)

#Running the differential expression 
dds_genes <- DESeq(dds_genes)
saveRDS(dds_genes, file.path(results.dir, "dds_group_genes.rds"))
