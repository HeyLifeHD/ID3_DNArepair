
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

library(Gviz)
library("EnsDb.Hsapiens.v86")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(org.Hs.eg.db)



#load Data
sample_anno<- readRDS(file.path(data.dir, "sample_anno.rds"))
AsiSI_100<- import.bed(file.path(data.dir,"ASIsites_100BEST_hg19.bed"))
HS_sites <- import.bed(file.path(data.dir, "HS_sites.bed")
HS_sites<- HS_sites[1:5,]
NHEJ_sites <- import.bed(file.path(data.dir, "NHEJ_sites.bed"))
NHEJ_sites<- NHEJ_sites[1:5,]
#get gene annotation
#ENS<-EnsDb.Hsapiens.v86
#seqlevelsStyle(ENS) <- "UCSC"
#GENE <- genes(ENS)
TXDB <- TxDb.Hsapiens.UCSC.hg19.knownGene

#choose the sites that should be investigated
SoI <- NHEJ_sites
#for loop for all regions
for (i in 1:length(SoI)) {
#get limits
lim <- c(start(SoI[i,]), end(SoI[i,]))
Chr<- unique(as.character(seqnames(SoI[i,])))
ext <- 2500

#gene annotation
grtrack <- GeneRegionTrack(TXDB, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, geneSymbol=TRUE, name="GeneModel")
symbols <- unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(grtrack) <- symbols[gene(grtrack)]

#get annotation tracks of cut sites
sites <- AnnotationTrack(SoI, name="AsiSI",    chromosome=Chr,genome = "hg19", id=paste0(seqnames(SoI[i,]), start(SoI[i,]), end(SoI[i,])))

#get ideogramm tracks
itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr)

#get data track from bam
bamFile <- sample_anno[7, "bw"]
dataTrack1neg <- DataTrack(range = bamFile, genome = "hg19", type = c("l","a"), name = "pos_wt",  chromosome =Chr, size=15)

bamFile <- sample_anno[1, "bw"]
dataTrack1pos <- DataTrack(range = bamFile, genome = "hg19",type =c("l","a"), name = "neg_wt",  chromosome =Chr, size=15)

bamFile <- sample_anno[8, "bw"]
dataTrack2neg <- DataTrack(range = bamFile, genome = "hg19",type =c("l","a") ,name = "pos_ARID1A",  chromosome =Chr, size=15)

bamFile <- sample_anno[2, "bw"]
dataTrack2pos <- DataTrack(range = bamFile, genome = "hg19",type = c("l","a"), name = "neg_ARID1A",  chromosome =Chr, size=15)

bamFile <- sample_anno[9, "bw"]
dataTrack3neg <- DataTrack(range = bamFile, genome = "hg19",type = c("l","a"), name = "pos_ARID1B",  chromosome =Chr, size=15)

bamFile <- sample_anno[3, "bw"]
dataTrack3pos <- DataTrack(range = bamFile, genome = "hg19",type = c("l","a"), name = "neg_ARID1B",  chromosome =Chr, size=15)

bamFile <- sample_anno[10, "bw"]
dataTrack4neg <- DataTrack(range = bamFile, genome = "hg19",type = c("l","a"), name = "pos_IDR3",  chromosome =Chr, size=15)

bamFile <- sample_anno[4, "bw"]
dataTrack4pos <- DataTrack(range = bamFile, genome = "hg19",type =c("l","a"), name = "neg_IDR3",  chromosome =Chr, size=15)

bamFile <- sample_anno[11, "bw"]
dataTrack5neg <- DataTrack(range = bamFile, genome = "hg19",type = c("l","a"), name = "pos_siMDC1",  chromosome =Chr, size=15)

bamFile <- sample_anno[5, "bw"]
dataTrack5pos <- DataTrack(range = bamFile, genome = "hg19",type =c("l","a"), name = "neg_siMDC1",  chromosome =Chr, size=15)

bamFile <- sample_anno[12, "bw"]
dataTrack6neg <- DataTrack(range = bamFile, genome = "hg19",type = c("l","a"), name = "pos_IDR3_siMDC1",  chromosome =Chr, size=15)

bamFile <- sample_anno[6, "bw"]
dataTrack6pos <- DataTrack(range = bamFile, genome = "hg19",type = c("l","a"), name = "neg_IDR3_siMDC1",  chromosome =Chr, size=15)

pdf(file.path(base.dir, "locus_plots_HS", paste0(seqnames(SoI[i,]),":", start(SoI[i,]), "-",end(SoI[i,]), ".pdf")), height=20)
plotTracks(list(itrack, grtrack, sites, dataTrack1neg, dataTrack1pos, dataTrack2neg, dataTrack2pos, dataTrack3neg, dataTrack3pos,  dataTrack4neg, dataTrack4pos, dataTrack5neg, dataTrack5pos, dataTrack6neg, dataTrack6pos),
from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
dev.off()
}

