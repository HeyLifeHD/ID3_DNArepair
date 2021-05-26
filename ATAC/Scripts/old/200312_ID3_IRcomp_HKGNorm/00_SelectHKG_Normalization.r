#create common deseq object
#libraries
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/"
data.dir <- file.path(base.dir, "data")
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
dir.create(base.dir)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)

#load data for irradiated
counts<- readRDS(file.path(data.dir , "counts.rds"))
sample_anno <- readRDS(file.path(data.dir, "sample_anno.rds"))
#select subset
sample_anno_sub <- sample_anno[sample_anno$id3_irradiation_project==TRUE,]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$tamoxifen_treatment =="UNTREATED", ]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$irradiation_treatment =="TREATED", ]
#sample_anno_sub <- sample_anno[sample_anno$genotype %in% c("WT", "ARID1A", "ID3"), ]
#sample_anno_sub <- sample_anno_sub[sample_anno_sub$tamoxifen_treatment =="UNTREATED",]
#sample_anno_sub <- sample_anno_sub[sample_anno_sub$hdac_treatment. =="UNTREATED",]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$Plattform =="Illumina NovaSeq 6000 S1",]
dim(sample_anno_sub)
sample_anno_sub$Replicate <- as.factor(sample_anno_sub$Replicate)
sample_anno_sub$Replicate <- droplevels(sample_anno_sub$Replicate)
sample_anno_sub$group <- droplevels(sample_anno_sub$group)

#subset
counts_sub_irr <- counts[,rownames(sample_anno_sub)]
sample_anno_sub$Replicate <-as.character(sample_anno_sub$Replicate )
sample_anno_sub_irr <- sample_anno_sub

#load data for unirradiated 
#libraries
library(DESeq2)
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/"
data.dir <- file.path(base.dir, "data")
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
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
#select subset
sample_anno_sub <- sample_anno[sample_anno$id3_irradiation_project==TRUE,]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$tamoxifen_treatment =="UNTREATED", ]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$irradiation_treatment =="UNTREATED", ]
#sample_anno_sub <- sample_anno[sample_anno$genotype %in% c("WT", "ARID1A", "ID3"), ]
#sample_anno_sub <- sample_anno_sub[sample_anno_sub$tamoxifen_treatment =="UNTREATED",]
#sample_anno_sub <- sample_anno_sub[sample_anno_sub$hdac_treatment. =="UNTREATED",]
sample_anno_sub <- sample_anno_sub[sample_anno_sub$Plattform =="Illumina NovaSeq 6000 S1",]
dim(sample_anno_sub)
sample_anno_sub$Replicate <- as.factor(sample_anno_sub$Replicate)
sample_anno_sub$Replicate <- droplevels(sample_anno_sub$Replicate)
sample_anno_sub$group <- droplevels(sample_anno_sub$group)

#subset
counts_sub_uir <- counts[,rownames(sample_anno_sub)]
sample_anno_sub$Replicate <-as.character(sample_anno_sub$Replicate )
sample_anno_sub_uir <- sample_anno_sub


#combine data
sample_anno <- rbind(sample_anno_sub_irr,sample_anno_sub_uir)
counts <- cbind(counts_sub_irr, counts_sub_uir)

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = sample_anno, 
                              design = ~ Replicate+group )

#Estimating size factors
dds <- estimateSizeFactors(dds)

#extract normalized counts
counts_norm <-counts(dds, normalized=TRUE)

#calculate rowmeans and their top 3 quantiles
counts_norm_mean <- rowMeans(counts_norm)
counts_norm_sub <- counts_norm[counts_norm_mean>quantile(counts_norm_mean)[2],]
counts_norm_sub <- as.data.frame(counts_norm_sub)

#get  annotation
data.dir <- "icgc/dkfzlsdf/project/OE0219/id3/sequencing/rna_sequencing/view-by-pid/"
#libraries
library(data.table)
#load counts
files <-list.files(path=data.dir,pattern=".featureCounts.tsv$", recursive = TRUE)
anno<- fread(file.path(data.dir,files[1]))
anno$ensembl <- sapply(strsplit(anno$gene_id, ".",fixed=TRUE),"[",1)
anno_sub <- anno[, c("ensembl","name" )]
counts_norm_sub$ensembl<- rownames(counts_norm_sub)
library(dplyr)
counts_norm_sub <- left_join(counts_norm_sub, anno_sub)

#select housekeeping genes
hkg <- readRDS("c010-datasets/Internal/ID3/analysis/data/all_hk_transcripts.rds")
counts_norm_sub_hkg <- counts_norm_sub[counts_norm_sub$name %in% hkg$gene_,]
dim(counts_norm_sub_hkg)

#order counts by sd
counts_norm_sub_hkg_ord <- counts_norm_sub_hkg[order(rowVars(as.matrix(counts_norm_sub_hkg[,rownames(sample_anno)])), decreasing=FALSE),]
head(counts_norm_sub_hkg_ord, 50)$name
rowVars(as.matrix(counts_norm_sub_hkg_ord[,rownames(sample_anno)]))

#get promoter regions and subset
require(EnsDb.Hsapiens.v75)
library(rtracklayer)
txdb <- EnsDb.Hsapiens.v75
gene_list <- genes(txdb)
gene_list <- promoters(gene_list,  upstream=1000,    downstream=1000 )
counts_norm_sub_hkg_ord_regions <- gene_list[gene_list$gene_id %in% head(counts_norm_sub_hkg_ord, 50)$ensembl,]
length(counts_norm_sub_hkg_ord_regions)
saveRDS(counts_norm_sub_hkg_ord_regions,file.path("c010-datasets/Internal/ID3/analysis/data/", "top50_stable_hkg.rds"))
#remove duplciated rows
mcols(counts_norm_sub_hkg_ord_regions) <- mcols(counts_norm_sub_hkg_ord_regions)[,-1]
counts_norm_sub_hkg_ord_regions[distinct(mcols(counts_norm_sub_hkg_ord_regions))),]
temp<-counts_norm_sub_hkg_ord_regions$symbol[duplicated(counts_norm_sub_hkg_ord_regions$symbol)]
C010-datasets/Internal/ID3/analysis/data/all_hk_transcripts.rds
rowVars(counts_norm_sub_ord)
#top 100
counts_norm_sub_hkg_ord_regions <- gene_list[gene_list$gene_id %in% head(counts_norm_sub_hkg_ord, 100)$ensembl,]
length(counts_norm_sub_hkg_ord_regions)
saveRDS(counts_norm_sub_hkg_ord_regions,file.path("c010-datasets/Internal/ID3/analysis/data/", "top100_stable_hkg.rds"))