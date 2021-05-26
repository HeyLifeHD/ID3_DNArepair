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
library(peakSeason)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#100 best cut sites
AsiSI_100<- import.bed(file.path(data.dir,"ASIsites_100BEST_hg19.bed"))

#sample annotation sheet 
sample_anno<- readRDS(file.path(data.dir, "sample_anno.rds"))


#Profile Plot 
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
#for genotype comparison
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= file.path(data.dir, "Promoters.bed"), genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir,"TSS", paste0(unique(sample_anno$Genotype)[Geno], "GenotypeComp_TSS_2500.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()
}

#stratify top 100 cutting sites based on promoter location
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
genome=hg19, up=500, down=500, binSize=2)
#draw profile plot 
pdf(file.path(analysis.dir,"Promoter_Damage", paste0(unique(sample_anno$Genotype)[Geno], "GenotypeComp_PromoterDamage_500.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()
}

#plot promoter of damage sites next to promoter
anno_AsiSi_100 <- annotatePeak(peak= AsiSI_100, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

#Promoter AsiSi
anno_AsiSi_100_df <- as.data.frame(anno_AsiSi_100)
anno_AsiSi_100_Promoter1kb <- anno_AsiSi_100_df[anno_AsiSi_100_df$annotation == "Promoter (<=1kb)", ]

#get promoters
genes <- genes(txdb)
prom <- promoters(genes, upstream=1, downstream=1)
prom <- prom[seqnames(prom) %in%  c(paste0("chr", 1:21),"chrY", "chrX")]

#subset prmoters with a damage
prom_damage <- prom[prom$gene_id %in% anno_AsiSi_100_Promoter1kb$geneId, ]
export.bed(prom_damage, file.path(data.dir,"Promoter_CloseToDamageAsiSi100best.bed"))
#plot genotype comparison
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, file.path(data.dir,"Promoter_CloseToDamageAsiSi100best.bed"), 
genome=hg19, up=500, down=500, binSize=2)
#draw profile plot 
pdf(file.path(analysis.dir,"Promoter_CloseToDamage", paste0(unique(sample_anno$Genotype)[Geno], "GenotypeComp_PromoterCloseToDamage_500.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()
}