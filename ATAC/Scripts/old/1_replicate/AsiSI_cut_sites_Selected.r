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
base.dir <- "c010-datasets/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")
temp.dir <- "/home/heyj/rTemp"

#Find out AsiSI cut sites int the hg19 genome
#in silico prediction of sites
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(data.table)

sample_anno<- readRDS(file.path(data.dir, "sample_anno.rds"))
sample_anno$bw <- str_replace(sample_anno$bw, "heyj", "epicwl")

#HS prone sites
HS_sites <- read.csv(file.path(data.dir, "HS_AsiSi_sites.csv"), sep=";")
HS_sites <- data.frame(lapply(HS_sites, function(x) {
                  gsub("Chr", "chr", x)
              }))
HS_sites_gr <- makeGRangesFromDataFrame(HS_sites, keep.extra.columns=TRUE)
chain<- import.chain(file.path(data.dir, "hg18ToHg19.over.chain"))
HS_sites_gr_lifted <- unlist(liftOver(HS_sites_gr, chain))
HS_sites_gr_lifted <- HS_sites_gr_lifted
export.bed(HS_sites_gr_lifted, file.path(data.dir, "HS_sites.bed"))


#all
#create column data
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", genome=hg19, up=100, down=100)
#draw profile plot 
pdf(file.path(analysis.dir, "all_Genotype_HS_profile_100.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10))
dev.off()

#all only positive
#create column data
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos",], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", genome=hg19, up=100, down=100)
#draw profile plot 
pdf(file.path(analysis.dir, "all_POS_Genotype_HS_profile_100.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10))
dev.off()

#for each subgroup
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", 
genome=hg19, up=500, down=500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "500_HS_profile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE, y_lim=c(0,10) )
dev.off()
}

#for each subgroup only positive
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos" & sample_anno$Genotype %in% c(as.character(unique(sample_anno$Genotype)[Geno]),"wt"),], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/epicwl/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", 
genome=hg19, up=500, down=500, binSize=5)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "andwt_POS_500_HS_profile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10) )
dev.off()
}


#NHEJ prone sites
NHEJ_sites <- read.csv(file.path(data.dir, "NHEJ_AsiSi_sites.csv"), sep=";")
NHEJ_sites <- data.frame(lapply(NHEJ_sites, function(x) {
                  gsub("Chr", "chr", x)
              }))
NHEJ_sites_gr <- makeGRangesFromDataFrame(NHEJ_sites, keep.extra.columns=TRUE)
chain<- import.chain(file.path(data.dir, "hg18ToHg19.over.chain"))
NHEJ_sites_gr_lifted <- unlist(liftOver(NHEJ_sites_gr, chain))
NHEJ_sites_gr_lifted <- NHEJ_sites_gr_lifted
export.bed(NHEJ_sites_gr_lifted, file.path(data.dir, "NHEJ_sites.bed"))

#all
#create column data
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, "all_Genotype_NHEJ_profile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10))
dev.off()

#all only positive
#create column data
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos",], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, "all_POS_Genotype_NHEJ_profile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE,  y_lim=c(0,10))
dev.off()

#for each subgroup
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", 
genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "2500_NHEJ_profile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE , y_lim=c(0,10))
dev.off()
}


#for each subgroup only positive
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos" & sample_anno$Genotype %in% c(as.character(unique(sample_anno$Genotype)[Geno]),"wt"),], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", 
genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "andwt_POS_2500_NHEJ_profile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Genotype", collapse_reps=FALSE , y_lim=c(0,10))
dev.off()
}












### only top 5
HS_sites_gr <- makeGRangesFromDataFrame(HS_sites, keep.extra.columns=TRUE)
chain<- import.chain(file.path(data.dir, "hg18ToHg19.over.chain"))
HS_sites_gr_lifted <- unlist(liftOver(HS_sites_gr, chain))
HS_sites_gr_lifted <- HS_sites_gr_lifted[1:5,]
export.bed(HS_sites_gr_lifted, file.path(data.dir, "HS_sites.bed"))

#all
#create column data
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, "all_Genotype_HS_top5_profile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10))
dev.off()

#all only positive
#create column data
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos",], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, "all_POS_Genotype_HS_top5_pprofile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10))
dev.off()

#for each subgroup
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", 
genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "2500_HS_top5_pprofile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE, y_lim=c(0,10) )
dev.off()
}

#for each subgroup only positive
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos" & sample_anno$Genotype %in% c(as.character(unique(sample_anno$Genotype)[Geno]),"wt"),], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/HS_sites.bed", 
genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "andwt_POS_2500_HS_top5_pprofile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10) )
dev.off()
}


#NHEJ prone sites
NHEJ_sites <- read.csv(file.path(data.dir, "NHEJ_AsiSi_sites.csv"), sep=";")
NHEJ_sites <- data.frame(lapply(NHEJ_sites, function(x) {
                  gsub("Chr", "chr", x)
              }))
NHEJ_sites_gr <- makeGRangesFromDataFrame(NHEJ_sites, keep.extra.columns=TRUE)
chain<- import.chain(file.path(data.dir, "hg18ToHg19.over.chain"))
NHEJ_sites_gr_lifted <- unlist(liftOver(NHEJ_sites_gr, chain))
NHEJ_sites_gr_lifted <- NHEJ_sites_gr_lifted[1:5,]
export.bed(NHEJ_sites_gr_lifted, file.path(data.dir, "NHEJ_sites.bed"))

#all
#create column data
col_data <- read_coldata(coldata=sample_anno, files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, "all_Genotype_NHEJ_top5_pprofile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE, y_lim=c(0,10))
dev.off()

#all only positive
#create column data
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos",], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, "all_POS_Genotype_NHEJ_top5_pprofile.pdf"))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=FALSE, condition="Genotype", collapse_reps=FALSE,  y_lim=c(0,10))
dev.off()

#for each subgroup
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Genotype==unique(sample_anno$Genotype)[Geno]], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", 
genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "2500_NHEJ_top5_pprofile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Treatment", collapse_reps=FALSE , y_lim=c(0,10))
dev.off()
}


#for each subgroup only positive
#create column data
for (Geno in 1:length(unique(sample_anno$Genotype))){
col_data <- read_coldata(coldata=sample_anno[sample_anno$Treatment =="pos" & sample_anno$Genotype %in% c(as.character(unique(sample_anno$Genotype)[Geno]),"wt"),], files_idx=9, sample_idx=1)
#create matrix list
matrix_data_ASItop100 <- extract_matrix(coldata=col_data, bed= "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/data/NHEJ_sites.bed", 
genome=hg19, up=2500, down=2500)
#draw profile plot 
pdf(file.path(analysis.dir, paste0(unique(sample_anno$Genotype)[Geno], "andwt_POS_2500_NHEJ_top5_pprofile.pdf")))
plot_profile(mat_list=matrix_data_ASItop100, summarizeBy="mean", ci=TRUE, condition="Genotype", collapse_reps=FALSE , y_lim=c(0,10))
dev.off()
}
