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
base.dir <- "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")
temp.dir <- "/home/heyj/rTemp"

#Find out AsiSI cut sites int the hg19 genome
#in silico prediction of sites
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(data.table)

dnaseq <- DNAString("GCGATCGC")

AsiSI<- vmatchPattern(dnaseq, BSgenome.Hsapiens.UCSC.hg19)
chr <- c(paste0("chr", 1:21), "chrY", "chrX")
AsiSI[seqnames(AsiSI)%in% chr,]
sub <- AsiSI[strand(AsiSI)=="+",]

AsiSIsub <- sub[seqnames(sub)%in% chr,]
export.bed(AsiSIsub, file.path(data.dir, "AsiSi_insilico.bed"))


#100 best cut sites
AsiSI_100<- import.bed(file.path(data.dir,"ASIsites_100BEST_hg19.bed"))


#sample annotation sheet 
sample_anno<- readRDS(file.path(data.dir, "sample_anno.rds"))


#Profile Plot 
library(peakSeason)
##100 most sign sites
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


#Profile plot of all TSS
genes <- genes(txdb)
prom <- promoters(genes, upstream=1, downstream=1)
prom <- prom[seqnames(prom) %in%  c(paste0("chr", 1:21),"chrY", "chrX")]
export.bed(prom, file.path(data.dir, "Promoters.bed"))

col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= file.path(data.dir, "Promoters.bed"), genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir,"new", "all_TSS_Genotype_profileplot.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE)
dev.off()

#stratify top 100 cutting sites
library(ChIPseeker)
#annotate DE regions
anno_AsiSi_100 <- annotatePeak(peak= AsiSI_100, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

#Promoter AsiSi
anno_AsiSi_100_df <- as.data.frame(anno_AsiSi_100)
anno_AsiSi_100_Promoter1kb <- anno_AsiSi_100_df[anno_AsiSi_100_df$annotation == "Promoter (<=1kb)", ]
anno_AsiSi_100_Promoter1kb <- unlist(makeGRangesListFromDataFrame(anno_AsiSi_100_Promoter1kb))
export.bed(anno_AsiSi_100_Promoter1kb, file.path(data.dir,"ASIsites_100BEST_Promoter.bed"))

for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, file.path(data.dir,"ASIsites_100BEST_Promoter.bed"), 
genome=hg19, up=500, down=500, binSize=5)
#draw profile plot 
pdf(file.path(analysis.dir,"new", paste0(unique(sample_anno$Genotype)[Geno], "Genotype_AsiSiPromoter_500_profile_500.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()
}

#Distal intergenic
anno_AsiSi_100_DG <- anno_AsiSi_100_df[anno_AsiSi_100_df$annotation == "Distal Intergenic", ]
anno_AsiSi_100_DG <- unlist(makeGRangesListFromDataFrame(anno_AsiSi_100_DG))
export.bed(anno_AsiSi_100_DG, file.path(data.dir,"ASIsites_100BEST_DistalIntergenic.bed"))

for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, file.path(data.dir,"ASIsites_100BEST_DistalIntergenic.bed"), 
genome=hg19, up=500, down=500, binSize=5)
#draw profile plot 
pdf(file.path(analysis.dir,"new", paste0(unique(sample_anno$Genotype)[Geno], "Genotype_AsiSiDistalIntergenic_500_profile_500.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()
}