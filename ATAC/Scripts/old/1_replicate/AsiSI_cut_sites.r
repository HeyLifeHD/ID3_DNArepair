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
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/"
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
pdf(file.path(analysis.dir, "all_TSS_Genotype_profileplot.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSEf)
dev.off()


library(ChIPseeker)
#annotate DE regions
anno_AsiSi_100 <- annotatePeak(peak= AsiSI_100, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

pdf(file.path(analysis.dir,   "AsiSI_100_DistanceTSSPC.pdf"))
plotDistToTSS(anno_AsiSi_100 , 
              title="Distribution of DMR relative to TSS")
dev.off()  

pdf(file.path(analysis.dir, "AsiSI_100_AnnoDist_Pie.pdf"))
plotAnnoPie(anno_AsiSi_100)
dev.off()

pdf(file.path(analysis.dir,"AsiSI_100_AnnoDist_Upset.pdf"))
upsetplot(anno_AsiSi_100, vennpie=TRUE)
dev.off()

#all
#create column data
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed", genome=hg19, up=10000, down=10000)
#draw profile plot 
pdf(file.path(analysis.dir,"new", "all_Treatment_AsiSi_100_10000_profile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()

pdf(file.path(analysis.dir,"new", "all_Genotype_AsiSi_100_10000_profile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE)
dev.off()

#for each subgroup
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed", 
genome=hg19, up=500, down=500, binSize=5)
#draw profile plot 
pdf(file.path(analysis.dir,"new", paste0(unique(sample_anno$Genotype)[Geno], "Genotype_AsiSi_100_profile_500.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()
}
#heatmap
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed", 
genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "_AsiSi_100_heatmap.pdf")))
plot_heatmap(matrix_data_ASItop100,zmaxs=15)
dev.off()
}

#heatmap 
#AsiSI 100Ja
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed", genome=hg19, up=2500, down=2500)
pdf(file.path(analysis.dir, "all_Genotype_AsiSi_100_profile.pdf"))
plot_heatmap(matrix_data_ASItop100, sortBy = "mean", rcb_pal = "Blues",
        title_size = 0.8, top_profile = FALSE,
       zmins = NULL, zmaxs = NULL, scale = FALSE, file_name = NULL,
       hm_width = NULL, hm_height = 12, ...)
dev.off()
## all in silico predicted sites
AsiSIsub<- import.bed(file.path(data.dir, "AsiSi_insilico.bed"))

#annotate DE regions
anno_AsiSIsub <- annotatePeak(peak= AsiSIsub, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

pdf(file.path(analysis.dir,   "AsiSIsub_DistanceTSSPC.pdf"))
plotDistToTSS(anno_AsiSIsub , 
              title="Distribution of DMR relative to TSS")
dev.off()  

pdf(file.path(analysis.dir, "AsiSIsub_AnnoDist_Pie.pdf"))
plotAnnoPie(anno_AsiSIsub)
dev.off()

pdf(file.path(analysis.dir,"AsiSIsub_AnnoDist_Upset.pdf"))
upsetplot(anno_AsiSIsub, vennpie=TRUE)
dev.off()

#all
#create column data
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC//data/AsiSi_insilico.bed", genome=hg19, up=100, down=100)
#draw profile plot 
pdf(file.path(analysis.dir,"new", "all_AsiSIsub_profile_100.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()

pdf(file.path(analysis.dir,"new", "all_Genotype_AsiSIsub_profile_100.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE)
dev.off()


#for each subgroup
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_AsiSIsub <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC//data/AsiSi_insilico.bed", 
genome=hg19, up=100, down=100, binSize=5)
#draw profile plot 
pdf(file.path(analysis.dir, "new",paste0(unique(sample_anno$Genotype)[Geno], "_AsiSIsub_profile_100.pdf")))
plot_profile(mat_list=matrix_data_AsiSIsub, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE)
dev.off()
}


#extrac area under the curve for every cut site
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
dir.create(file.path(analysis.dir, "summaries"))
sum_ASItop100<- extract_summary(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed",  
up=500, down=500, , startFrom="center",rmAfter=FALSE ,op_dir=file.path(analysis.dir, "summaries"))
#since function is not working get the summary lists
temp = list.files(path=file.path(analysis.dir, "summaries"), pattern="*.summary")
temp<- paste0(file.path(analysis.dir, "summaries/"), temp)
myfiles = lapply(temp, fread)
#extract sums
sums <-lapply(myfiles, "[", , "sum")
sums<- do.call("cbind", sums)
sums <- as.data.frame(sums)
colnames(sums)<- list.files(path=file.path(analysis.dir, "summaries"), pattern="*.summary")

pheatmap(sums, scale="row",filename=file.path(analysis.dir,"AsiSi_100_Heatmap_2500.pdf") )


col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
sum_ASItop100<- extract_summary(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed",  
up=100, down=100, , startFrom="center",rmAfter=FALSE ,remove_dups=FALSE)
#extract sums
sums <- as.data.frame(sum_ASItop100$summaries[,-c(1:8)])
#log2 transform
sums <- log2(sums+0.1)
sums <- sums[complete.cases(sums),]
anno_sub <- as.data.frame(sample_anno[,c(2:4)])
rownames(anno_sub)<- sample_anno$SampleID

pheatmap(sums, scale="none",annotation_col=anno_sub,filename=file.path(analysis.dir,"AsiSi_100_Heatmap_100.pdf"))

scale_rows <- function(x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}
sums_scaled <- scale_rows(sums)
sums_scaled <- sums_scaled[complete.cases(sums_scaled),]

pheatmap(sums_scaled, annotation_col=anno_sub,scale="none",filename=file.path(analysis.dir,"AsiSi_100_Heatmap_100_scaled.pdf"))



#plot pca
sample_anno$SampleID2 <- sample_anno$SampleID
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
sum_ASItop100<- extract_summary(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed",  
up=100, down=100, , startFrom="center",rmAfter=FALSE ,remove_dups=FALSE)

pdf(file.path(analysis.dir, "ASIsite100_pca_100bp.pdf"))
bw_pca(sum_ASItop100, condition="Genotype")
dev.off()