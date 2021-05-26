
#Directories
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp_HKGNorm"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")
dir.create(analysis.dir, recursive=TRUE)

#Libraries
library(DiffBind)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ChIPseeker)
library(rtracklayer)
#Create Sample sheet
#get files 
bam_files <- list.files(data.dir, pattern="_tn5_center_73bp.bam$", full.names=TRUE, recursive=TRUE)
bam_files <- bam_files[grep("cache", bam_files, invert=TRUE)]
peak_files <- list.files(data.dir, pattern="_tn5_center_73bp.macs2_peaks.broadPeak$", full.names=TRUE, recursive=TRUE)
peak_files <- peak_files[grep("cache", peak_files, invert=TRUE)]
#GET NAMES
file_names <- sapply(strsplit(peak_files, "/", fixed=TRUE), "[", 11)
file_names <-gsub("ID3_KO", "ID3-ko", file_names)
sample_anno <- data.frame("SampleID"=file_names,
    "Factor"= sapply(strsplit(file_names, "_", fixed=TRUE), "[", 1),
    "Treatment"= sapply(strsplit(file_names, "_", fixed=TRUE), "[", 2),
    "Replicate"= paste0("rep",sapply(strsplit(file_names, "_", fixed=TRUE), "[", 3)),
    "PeakCaller"= rep("macs", length(file_names)),
    "Peaks"=peak_files,
    "bamReads"=bam_files)
sample_anno$Condition <- paste0(sample_anno$Factor, "_", sample_anno$Treatment)
sample_anno$SampleID <- gsub("ID3-ko", "ID3.ko",sample_anno$SampleID)
rownames(sample_anno)<- sample_anno$SampleID

dir.create(file.path(analysis.dir, "data"))
saveRDS(sample_anno,file.path(analysis.dir, "data","sample_anno.rds" ))
sample_anno <- readRDS(file.path(analysis.dir, "data","sample_anno.rds" ))
#create dataset
dataset <-dba(sampleSheet=sample_anno)
saveRDS(dataset, file.path(analysis.dir, "data","Dataset_peak.rds"))

#correlation Heatmap
dir.create(file.path(analysis.dir, "Expl_Diffbind"))
pdf(file.path(analysis.dir, "Expl_Diffbind", "correlation_heatmap_peaks.pdf"))
plot(dataset)
dev.off()

#count reads for each peak
dataset<- dba.count(dataset,  bParallel=T,  score=DBA_SCORE_READS)
saveRDS(dataset, file.path(analysis.dir, "data", "Dataset_count.rds"))
#dataset<- readRDS( file.path(analysis.dir, "data","Dataset_count.rds"))

#correlation Heatmap
pdf(file.path(analysis.dir, "Expl_Diffbind", "correlation_heatmap_affinityscore.pdf"))
plot(dataset)
dev.off()

#PCA analysis
pdf(file.path(analysis.dir, "Expl_Diffbind", "pca_affinity_score.pdf"))
dba.plotPCA(dataset, label=DBA_CONDITION)
dev.off()

#extract counts
counts_raw <- dba.peakset(dataset, bRetrieve=TRUE)
counts_raw <- dba.peakset(dataset, bRetrieve=TRUE)

saveRDS(counts_raw, file.path(analysis.dir,  "data","allPeaks_raw_all.rds"))
counts_raw <- counts_raw[seqnames(counts_raw) %in% paste0("chr",c(1:19, "X","Y")),]
write.table(as.data.frame(counts_raw), file.path(analysis.dir, "data", "allPeaks_raw.txt"),quote=F)
export.bed(counts_raw, file.path(analysis.dir,  "data","allPeaks_raw.bed"))
saveRDS(counts_raw, file.path(analysis.dir,  "data","allPeaks_raw.rds"))