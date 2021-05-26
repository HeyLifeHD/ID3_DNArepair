
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
#Create Sample sheetich
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
write.table(sample_anno, file.path(analysis.dir,"sample_anno.txt" ), row.names=FALSE, sep="\t")
rownames(sample_anno)<- sample_anno$SampleID
saveRDS(sample_anno,file.path(analysis.dir,"sample_anno.rds" ))
sample_anno <- readRDS(file.path(analysis.dir,"sample_anno.rds" ))
#create dataset
dataset <-dba(sampleSheet=sample_anno)
saveRDS(dataset, file.path(analysis.dir,"Dataset_peak.rds"))

#correlation Heatmap
pdf(file.path(analysis.dir, "correlation_heatmap_peaks.pdf"))
plot(dataset)
dev.off()
#count reads for each peak
dataset<- dba.count(dataset,  bParallel=T)
saveRDS(dataset, file.path(analysis.dir,"Dataset_count.rds"))
#dataset<- readRDS( file.path(analysis.dir,"Dataset_count.rds"))

#correlation Heatmap
pdf(file.path(analysis.dir,  "correlation_heatmap_affinityscore.pdf"))
plot(dataset)
dev.off()

#PCA analysis
pdf(file.path(analysis.dir,  "pca_affinity_score.pdf"))
dba.plotPCA(dataset, label=DBA_CONDITION)
dev.off()

#set up the contrast DAR analysis
dataset<- dba.contrast(dataset, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)

#run DAR analysis
dataset <- dba.analyze(dataset,method=DBA_ALL_METHODS, bParallel=F)
saveRDS(dataset, file.path(analysis.dir, "diffbind.rds"))
#dataset<- readRDS(file.path(analysis.dir, "diffbind.rds"))

#extact DAR regions for each of the comparisons
contrasts <- data.frame(contrasts=1:length(unlist(lapply(dataset$contrasts, function(x)x$name1))), group1=c(unlist(lapply(dataset$contrasts, function(x)x$name1))),
group2= c(unlist(lapply(dataset$contrasts, function(x)x$name2))))
contrasts$comparison <- paste0(contrasts$group1, "_vs_", contrasts$group2)
saveRDS(contrasts, file.path(analysis.dir, "contrasts"))

#write for loop to extract all possible comparisons in different plots
dir.create(file.path(analysis.dir, "Anno"))
dir.create(file.path(analysis.dir, "MAs"))
dir.create(file.path(analysis.dir, "Volcanos"))
dir.create(file.path(analysis.dir, "Boxplots"))
dir.create(file.path(analysis.dir, "Heatmaps"))

DE_list <- list(NULL)
for (i in 1:nrow(contrasts)){
    #extract DARs
    #specific for this mixed analysis
    if(i <=2){
        temp <- dba.report(dataset, contrast=i ,method =DBA_DESEQ2_BLOCK, th=0.05)
    } else{
        temp <- dba.report(dataset, contrast=i ,method =DBA_DESEQ2_BLOCK, th=0.05)
    }
    if(is.null(temp)) {}
    if(length(temp)==1) {
    temp_anno<- annotatePeak(peak= temp, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
    #make df
    temp_anno_df <- as.data.frame(temp_anno)
    temp_anno_gr <- makeGRangesFromDataFrame(temp_anno_df, keep.extra.columns=TRUE)
    #make a list
    DE_list[[i]]<- temp_anno_gr
    #plot anno
    pdf(file.path(analysis.dir, "Anno",  paste0("Anno_", contrasts[i, 4], ".pdf")))
    print(plotAnnoPie(temp_anno))
    dev.off()  
    #plot MA
    pdf(file.path(analysis.dir, "MAs",  paste0("MA_", contrasts[i, 4], ".pdf")))
    print(dba.plotMA(dataset,contrast=i, method= DBA_DESEQ2_BLOCK))
    dev.off()   
    #volcano plot
    pdf(file.path(analysis.dir, "Volcanos",  paste0("Vol_", contrasts[i, 4], ".pdf")))
    print(dba.plotVolcano(dataset,contrast=i, method= DBA_DESEQ2_BLOCK))
    dev.off()
    }
    if(length(temp)>1){
    #annotate DE regions
    temp_anno<- annotatePeak(peak= temp, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
    #make df
    temp_anno_df <- as.data.frame(temp_anno)
    temp_anno_gr <- makeGRangesFromDataFrame(temp_anno_df, keep.extra.columns=TRUE)
    #make a list
    DE_list[[i]]<- temp_anno_gr
    #plot anno
    pdf(file.path(analysis.dir, "Anno",  paste0("Anno_", contrasts[i, 4], ".pdf")))
    print(plotAnnoPie(temp_anno))
    dev.off()  
    #plot distance to TSS
    pdf(file.path(analysis.dir, "Anno",  paste0("DisTSS_", contrasts[i, 4], ".pdf")))
    print(print(plotDistToTSS(temp_anno)))
    dev.off()
    #plot MA
    pdf(file.path(analysis.dir, "MAs",  paste0("MA_", contrasts[i, 4], ".pdf")))
    print(dba.plotMA(dataset,contrast=i, method= DBA_DESEQ2_BLOCK))
    dev.off()   
    #volcano plot
    pdf(file.path(analysis.dir, "Volcanos",  paste0("Vol_", contrasts[i, 4], ".pdf")))
    print(dba.plotVolcano(dataset,contrast=i, method= DBA_DESEQ2_BLOCK))
    dev.off()
    #Heatmap
    pdf(file.path(analysis.dir, "Heatmaps",  paste0("Heat_", contrasts[i, 4], ".pdf")))
    print(dba.plotHeatmap(dataset, contrast=i, correlations=FALSE, scale="row", method=DBA_DESEQ2_BLOCK))
    dev.off()
    #Boxplot
    pdf(file.path(analysis.dir, "Boxplots",  paste0("Box_", contrasts[i, 4], ".pdf")))
    print(dba.plotBox(dataset,contrast=i, method=DBA_DESEQ2_BLOCK))
    dev.off()
    }
}
lapply(DE_list, function(x)length(x))
#subset results only for canonical chromosomes
DE_list_sub <- lapply(DE_list, function(x){
    x <- x[seqnames(x) %in% paste0("chr",c(1:19, "X","Y")),]
})
lapply(DE_list_sub, function(x)length(x))

names(DE_list) <- contrasts$comparison
saveRDS(DE_list, file.path(analysis.dir, "de_list.rds"))
DE_list <- readRDS(file.path(analysis.dir, "de_list.rds"))

# retrieve complete list
#get complete peak set TPM
#remove non canonical
counts <- dba.peakset(dataset, bRetrieve=TRUE)
counts <- counts[seqnames(counts) %in% paste0("chr",c(1:19, "X","Y")),]
saveRDS(counts, file.path(analysis.dir, "countsperpeak.rds"))
#counts<- readRDS(file.path(analysis.dir, "countsperpeak.rds"))
write.table(as.data.frame(counts), file.path(analysis.dir, "allPeaks.txt"),quote=F)
export.bed(counts, file.path(analysis.dir, "allPeaks.bed"))

#get raw couts
dataset_raw<- readRDS( file.path(analysis.dir,"Dataset_count.rds"))
counts_raw <- dba.peakset(dataset_raw, bRetrieve=TRUE,score=DBA_SCORE_READS)
counts_raw <- counts_raw[seqnames(counts_raw) %in% paste0("chr",c(1:19, "X","Y")),]
write.table(as.data.frame(counts_raw), file.path(analysis.dir, "allPeaks_raw.txt"),quote=F)
export.bed(counts_raw, file.path(analysis.dir, "allPeaks_raw.bed"))
saveRDS(counts_raw, file.path(analysis.dir, "allPeaks_raw.rds"))

#get all peaks with values
DE_list_complete <- list(NULL)
for (i in 1:nrow(contrasts)){
    #specific for this mixed analysis
    if(i <=2){
        temp <- dba.report(dataset, contrast=i ,method =DBA_DESEQ2_BLOCK, th=1)
    } else{
        temp <- dba.report(dataset, contrast=i ,method =DBA_DESEQ2_BLOCK, th=1)
    }    
    DE_list_complete[[i]]<- temp
}
names(DE_list_complete)<- contrasts$comparison

#get regions on canonical chromosomes
DE_list_complete <- lapply(DE_list_complete, function(x){
    x <- x[seqnames(x) %in% paste0("chr",c(1:19, "X","Y")),]
})

#order and join with counts
counts <- sortSeqlevels(counts)
counts <- sort(counts)
DE_list_complete <- lapply(DE_list_complete, function(x){
    x<- sortSeqlevels(x)
    x <- sort(x)
    mcols(x)<- cbind(mcols(x), mcols(counts))
    x <- x[order(x$FDR, decreasing=F),]
    x
})'
#annotate them
DE_list_complete <- lapply(DE_list_complete, function(x){
        temp_anno<- annotatePeak(peak= x, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
        temp_anno_df <- as.data.frame(temp_anno)
        x <- makeGRangesFromDataFrame(temp_anno_df, keep.extra.columns=TRUE)
        x
}) 
names(DE_list_complete)<- contrasts$comparison
saveRDS(DE_list_complete, file.path(analysis.dir, "de_list_complete.rds"))
#DE_list_complete<- readRDS(file.path(analysis.dir, "de_list_complete.rds"))
for (i in names(DE_list_complete)){
    write.table(as.data.frame(DE_list_complete[[i]]), file.path(analysis.dir, paste0("All_peaks_list_",i,".txt")),quote=F, sep="\t", row.names=FALSE)
    write.table(as.data.frame(DE_list_complete[[i]][which(DE_list_complete[[i]]$FDR < 0.05),]), file.path(analysis.dir, paste0("All_Sig_peaks_list_",i,".txt")),quote=F, sep="\t", row.names=FALSE)
    write.table(as.data.frame(DE_list_complete[[i]][which(DE_list_complete[[i]]$FDR < 0.05 & abs(DE_list_complete[[i]]$Fold)>1 ),]), file.path(analysis.dir, paste0("All_SigFC_peaks_list_",i,".txt")),quote=F, sep="\t", row.names=FALSE)
    write.table(as.data.frame(DE_list_complete[[i]][DE_list_complete[[i]]$Fold >0,]), file.path(analysis.dir, paste0("Open_peaks_list_",i,".txt")),quote=F, sep="\t", row.names=FALSE)
    write.table(as.data.frame(DE_list_complete[[i]][DE_list_complete[[i]]$Fold <0,]), file.path(analysis.dir, paste0("Closed_peaks_list_",i,".txt")),quote=F, sep="\t", row.names=FALSE)
}

# complete heatmap of all DARs
DE_list_df <-lapply(DE_list, function(x) as.data.frame(x))
DE_list_df <-DE_list_df[c(which(!unlist(lapply(DE_list_df, function(x)nrow(x)))==0))]
for (i in 1:length(DE_list_df)){    DE_list_df[[i]]<- DE_list_df[[i]][abs(DE_list_df[[i]]$Fold)>1,]
}
DE_list_df <-lapply(DE_list_df, "[" ,1:5)
DE_list_df <- do.call("rbind", DE_list_df)
DE_list_gr <- makeGRangesFromDataFrame(DE_list_df)
length(DE_list_gr)
DE_list_gr_unique <- unique(DE_list_gr)
length(DE_list_gr_unique)

#find overlap between DAR and peak set
overlaps <- findOverlaps(DE_list_gr_unique, counts)
#subset counts based on it being dmrs
counts_DARs <- counts[subjectHits(overlaps), ]
length(counts_DARs)
plot <- as.data.frame(mcols(counts_DARs))
#topVar<- head(order(rowVars(as.matrix(plot)), decreasing=TRUE),20000)
#plot <- plot[topVar,]
#plot heatmap
annorld<- sample_anno[,c("Factor","Treatment")]
library(pheatmap)
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(annorld$Factor)), name = 'Dark2')
names(color_genotype)<- as.character(unique(annorld$Factor))
color_treatment<-RColorBrewer::brewer.pal(n = length(unique(annorld$Treatment)), name = 'Accent')
names(color_treatment)<- as.character(unique(annorld$Treatment))
anno_colors <- list(Factor=color_genotype,  Treatment=color_treatment)

pheatmap(plot, scale="none",annotation_col=annorld, clustering_distance_cols="correlation", annotation_colors=anno_colors, 
    filename=file.path(analysis.dir,"Heatmaps","allDAR_Heatmap.pdf"),colnames=FALSE)
pheatmap(plot, scale="row",annotation_col=annorld,clustering_distance_cols="correlation", annotation_colors=anno_colors,  
    filename=file.path(analysis.dir,"Heatmaps","allDAR_Heatmap_zscaled.pdf"),colnames=FALSE)

'