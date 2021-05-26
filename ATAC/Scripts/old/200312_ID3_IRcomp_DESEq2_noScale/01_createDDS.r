
#Directories
base_old.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp_HKGNorm"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis_old.dir <- file.path(base_old.dir, "analysis", "DiffBind")
dir.create(analysis.dir, recursive=TRUE)

#Libraries
library(DiffBind)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ChIPseeker)
library(rtracklayer)
library(DESeq2)

#load data
counts_raw <- readRDS(file.path(analysis_old.dir,  "data","allPeaks_raw_gr.rds"))
sample_anno <- readRDS(file.path(analysis_old.dir, "data","sample_anno.rds" ))
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp_DESeq_noScale"
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")
dir.create(analysis.dir, recursive=TRUE)
#annotate peaks
counts_raw_anno<- annotatePeak(peak= counts_raw, tssRegion=c(-3000, 3000),
                        TxDb=txdb, annoDb="org.Hs.eg.db")
                        #plot anno
#plot anno
dir.create(file.path(analysis.dir, "Anno"))
pdf(file.path(analysis.dir, "Anno",  paste0("Anno_", "allpeaks", ".pdf")))
plotAnnoPie(counts_raw_anno)
dev.off()  
pdf(file.path(analysis.dir, "Anno",   paste0("DistTSS_", "allpeaks", ".pdf")))
plotDistToTSS(counts_raw_anno)
dev.off()
pdf(file.path(analysis.dir, "Anno",   paste0("UpsetAnno_", "allpeaks", ".pdf")))
upsetplot(counts_raw_anno)
dev.off()
#make df
counts_raw_anno_df <- as.data.frame(counts_raw_anno)
counts_raw_gr <- makeGRangesFromDataFrame(counts_raw_anno_df, keep.extra.columns=TRUE)
saveRDS(counts_raw_gr,file.path(analysis.dir,  "data","allPeaks_raw_gr.rds"))
#create dds
#create dds
dds <- DESeqDataSetFromMatrix(countData = mcols(counts_raw_gr)[, rownames(sample_anno)], 
                              colData = sample_anno, 
                              design = ~ Condition, rowRanges=counts_raw_gr)

#Running the DESEQ
dds <- DESeq(dds)

#Plot dispersion
dir.create(file.path(analysis.dir, "DESEQ"))
pdf(file.path(analysis.dir, "DESEQ","Dispersion(1).pdf"))
DESeq2::plotDispEsts(dds, main="Dispersion plot")
dev.off()
saveRDS(dds,file = file.path(analysis.dir, "DESEQ","dds.rds"))
