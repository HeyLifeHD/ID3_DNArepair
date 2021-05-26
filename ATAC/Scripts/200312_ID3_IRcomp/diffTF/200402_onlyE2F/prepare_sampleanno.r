
#Directories
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp"
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

#adapt for diffTF
sample_anno <- sample_anno[,c("SampleID", "Peaks", "bamReads", "Condition")]
sample_anno$conditionSummary <- sample_anno$Condition
sample_anno$Peaks <- paste0("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data/",sapply(strsplit(peak_files,"/",fixed=TRUE ), "[", 12))
sample_anno$bamReads <- paste0("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data/",sapply(strsplit(bam_files,"/",fixed=TRUE ), "[", 12))
sample_anno_sub <- sample_anno[sample_anno$Condition %in% c("ID3-ko_Irr", "wt_Irr"), ]
write.table(sample_anno_sub, file.path("c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data","sample_anno.txt" ), row.names=FALSE, sep="\t", quote=FALSE)

sample_anno_sub2 <- sample_anno[sample_anno$Condition %in% c("ID3-ko_unIrr", "wt_unIrr"), ]
write.table(sample_anno_sub2, file.path("c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data","sample_anno_UnIrr.txt" ), row.names=FALSE, sep="\t", quote=FALSE)
