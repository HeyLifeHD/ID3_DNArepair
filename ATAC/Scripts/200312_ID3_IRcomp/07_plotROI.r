#directories
#atac
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")
#RNA 
#Irradiation
base_Irr.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data_Irr.dir <- file.path(base_Irr.dir, "data")
base_results_Irr.dir <- file.path(base_Irr.dir, "results")
results_Irr.dir<- file.path(base_results_Irr.dir , "tables")
PreDE_Irr.dir <- file.path(base_results_Irr.dir,"PreDE")
PostDE_Irr.dir <- file.path(base_results_Irr.dir,"PostDE")

#Unirradiated
base_UIrr.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data_UIrr.dir <- file.path(base_UIrr.dir, "data")
base_results_UIrr.dir <- file.path(base_UIrr.dir, "results")
results_UIrr.dir<- file.path(base_results_UIrr.dir , "tables")
PreDE_UIrr.dir <- file.path(base_results_UIrr.dir,"PreDE")
PostDE_UIrr.dir <- file.path(base_results_UIrr.dir,"PostDE")

#output
output.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/Locusplots/"
dir.create(output.dir,recursive=TRUE)

#libraries#Libraries
library(parallel)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(rtracklayer)
library(rtracklayer)
library(Gviz)

#load data
#ATAC
#DE_list <- readRDS(file.path(analysis.dir, "de_list.rds"))
#counts<- readRDS(file.path(analysis.dir, "countsperpeak.rds"))
#counts_raw<- readRDS(file.path(analysis.dir, "allPeaks_raw.rds"))
#sample_anno<- readRDS(file.path(analysis.dir, "sample_anno.rds"))
DAR_list_complete<- readRDS(file.path(analysis.dir, "de_list_complete.rds"))

#RNA
#Irr
DEG_results_list_Irr<- readRDS(file.path(PostDE_Irr.dir,  "DEG_results_group_list.rds"))

#UnIrr
DEG_results_list_UIrr<- readRDS(file.path(PostDE_UIrr.dir,  "DEG_results_group_list.rds"))

#join list of relevant rna seq results
DEG_results_list <- list(DEG_results_list_Irr$"ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED",
    DEG_results_list_UIrr$"ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED")
names(DEG_results_list) <- c("ID3-ko_Irr_vs_wt_Irr", "ID3-ko_unIrr_vs_wt_unIrr")

#subset ATAC comparisons
DAR_list_complete_sub <- DAR_list_complete[names(DEG_results_list)]

#subset significant ATAC DARs
DAR_list_complete_sign <- lapply(DAR_list_complete_sub, function(x){
    x <- x[x$FDR<0.05,]
    x
})
lapply(DAR_list_complete_sign, function(x)length(x))
#subset fold change
DAR_list_complete_sign_fc <- lapply(DAR_list_complete_sign, function(x){
    x <- x[abs(x$Fold)>1,]
    x
})
lapply(DAR_list_complete_sign_fc, function(x)length(x))

dars <- GRangesList(DAR_list_complete_sign_fc$"ID3-ko_Irr_vs_wt_Irr" ,DAR_list_complete_sign_fc$"ID3-ko_unIrr_vs_wt_unIrr")
dars <- unlist(dars)
dars <- reduce(dars)

#Create Sample sheet
#get files 
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")

bigwig_files <- list.files(data.dir, pattern="_tn5_center_73bp.bigwig$", full.names=TRUE, recursive=TRUE)
bigwig_files <- bigwig_files[grep("cache", bigwig_files, invert=TRUE)]

#GET NAMES
file_names <- sapply(strsplit(bigwig_files, "/", fixed=TRUE), "[", 11)
file_names <-gsub("ID3_KO", "ID3-ko", file_names)
sample_anno <- data.frame("SampleID"=file_names,
    "Factor"= sapply(strsplit(file_names, "_", fixed=TRUE), "[", 1),
    "Treatment"= sapply(strsplit(file_names, "_", fixed=TRUE), "[", 2),
    "Replicate"= paste0("rep",sapply(strsplit(file_names, "_", fixed=TRUE), "[", 3)),
    "PeakCaller"= rep("macs", length(file_names)),
    "bigwig"=bigwig_files)
sample_anno$Condition <- as.factor(paste0(sample_anno$Factor, "_", sample_anno$Treatment))
sample_anno$SampleID <- gsub("ID3-ko", "ID3.ko",sample_anno$SampleID)
rownames(sample_anno)<- sample_anno$SampleID
sample_anno$bigwig <- as.character(sample_anno$bigwig)

#get atac track
atac_bw <- lapply(sample_anno$bigwig,import.bw)
names(atac_bw)<- sample_anno$SampleID
saveRDS(atac_bw,file.path(analysis.dir, "all_atac_bw_gr.rds"))
atac_bw <- readRDS(file.path(analysis.dir, "all_atac_bw_gr.rds"))
#remove non canonical chromosomes
atac_bw_sub <- lapply(atac_bw, function(x){
    x <- x[seqnames(x) %in% paste0("chr",c(1:21, "X","Y")),]
    seqlevels(x)<- paste0("chr",c(1:21, "X","Y"))
    x
})
atac_bw_cov <-  mclapply(atac_bw_sub, function(x)GenomicRanges::coverage(x, weight="score"),mc.cores=1)
saveRDS(atac_bw_cov, file.path(analysis.dir, "atac_bw_cov.rds"))
atac_bw_cov <- readRDS(file.path(analysis.dir, "atac_bw_cov.rds"))
#summarize bedgraphs
chrlengths <- seqlengths(atac_bw_sub[[1]])
chrlengths <- GRanges(
  seqnames = names(chrlengths),
 ranges = IRanges(start = rep(1, length(chrlengths)),
                   end = chrlengths
  )
)
tiles <- unlist(tile(range(chrlengths), width=20))
seqlevels(tiles)<- names(atac_bw_cov[[1]])
atac_bw_cov_tiled <-   mclapply(atac_bw_cov, function(x)GenomicRanges::binnedAverage(tiles,x,"score"),mc.cores=1)
saveRDS(atac_bw_cov_tiled, file.path(analysis.dir, "atac_bw_cov_tiled.rds"))
atac_bw_cov_tiled <- readRDS(file.path(analysis.dir, "atac_bw_cov_tiled.rds"))
scores <- lapply(atac_bw_cov_tiled, function(x)mcols(x))
scores <- do.call("cbind", scores)
colnames(scores) <- as.character(sample_anno$SampleID)
scores <- as.matrix(scores)*20
saveRDS(scores, file.path(analysis.dir, "scores.rds"))
scores <- readRDS(file.path(analysis.dir, "scores.rds"))
atac_gr <- tiles
mcols(atac_gr)<- scores
saveRDS(atac_gr,file.path(analysis.dir, "allATAC_gr.rds"))

atac_gr <-readRDS(file.path(analysis.dir, "allATAC_gr.rds"))



#Gene gene annotation
library(EnsDb.Hsapiens.v75)
ENS<-EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"
GENE <- genes(ENS)

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
TXDB <- TxDb.Hsapiens.UCSC.hg19.knownGene

#select genes
dir.create(file.path(output.dir,"selectedGenes"))
#goi <- unlist(lapply(c("Irf"), function(x) GENE[grep(x, GENE$symbol),]$symbol))
#
goi <- c("XAB2","EXO1","RBBP8","FANCM","BRCA2","RAD51","POLQ","RCF4")

#plot regions
for (i in goi){
    #Define Location
    if(i %in% GENE$symbol){

        lim <- c(start(GENE[GENE$symbol==i,])[1], end(GENE[GENE$symbol==i,])[1])
        Chr<- unique(as.character(seqnames(GENE[GENE$symbol==i,])[1]))
        ext <- 2000 
        range<- GRanges(
        seqnames = Chr,
        ranges = IRanges(start = lim[1]-ext,
                    end = lim[2]+ext))

        #gene annotation
        grtrack <- GeneRegionTrack(TXDB, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, geneSymbol=TRUE, 
            name="Gene Model EnsDB v.79",collapseTrack=TRUE,
            start =lim[1]-ext, end = lim[2]+ext,)
        symbols <- unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
        symbol(grtrack) <- symbols[gene(grtrack)] 
        
        #get dmr annotation
        dars_track <- AnnotationTrack(dars[seqnames(dars)==Chr,], name="red.\nDARs", chromosome=Chr,genome = "hg19", fill="darkgray",
            start =lim[1]-ext, end = lim[2]+ext,)

        #get ideogramm tracks
        itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr)

        #genome axis track
        getrack <- GenomeAxisTrack()
        #Data Track
        col <- c("darkorange4", "darkorange", "#7570B3",  "cyan3")

        trackall <- DataTrack(range = atac_gr, genome = "hg19",start =lim[1]-ext, end = lim[2]+ext,
        type = c("smooth"), chromosome = Chr, name = "Chromatin accessibility",groups=as.character(sample_anno$Condition),col=col, span=.05)

        #Plot track
        pdf(file.path(output.dir,"selectedGenes",paste0("locusPlot_Gene_",i,"_v1.pdf")), width = 7, height = 6)
        plotTracks(list(itrack ,grtrack,getrack,dars_track, trackall),from =lim[1]-ext, to = lim[2]+ext, 
        chromosome=Chr)
        dev.off()
        print(i)
    }
}


library(ChIPseeker)
dars_anno<- annotatePeak(peak= dars, tssRegion=c(-3000, 3000),
                         TxDb=TXDB, annoDb="org.Hs.eg.db")
dars_anno_df <- as.data.frame(dars_anno)
dars <- makeGRangesFromDataFrame(dars_anno_df, keep.extra.columns=TRUE)
    
dir.create(file.path(output.dir,"selectedDARs"))
"EXO1","RBBP8",
goi <- c("FANCM","BRCA2","RAD51","POLQ","RCF4","XAB2")
#plot dars regions
for (i in goi){
    #Define Location
    if(i %in% dars$SYMBOL){
        lim <- c(start(dars[which(dars$SYMBOL==i),])[1], end(dars[which(dars$SYMBOL==i),])[1])
        Chr<- unique(as.character(seqnames(dars[which(dars$SYMBOL==i),])[1]))
        ext <- 5000 
        range<- GRanges(
        seqnames = Chr,
        ranges = IRanges(start = lim[1]-ext,
                    end = lim[2]+ext))

        #gene annotation
        grtrack <- GeneRegionTrack(TXDB, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, geneSymbol=TRUE, 
            name="Gene Model EnsDB v.79",collapseTrack=TRUE)
        symbols <- unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
        symbol(grtrack) <- symbols[gene(grtrack)] 
    
        #get dmr annotation
        dars_track <- AnnotationTrack(dars[seqnames(dars)==Chr,], name="red.\nDARs", chromosome=Chr,genome = "hg19", fill="darkgray",
            start =lim[1]-ext, end = lim[2]+ext,)

        #get ideogramm tracks
        itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr)

        #genome axis track
        getrack <- GenomeAxisTrack()
        #Data Track
        col <- c("darkorange4", "darkorange", "#7570B3",  "cyan3")

        trackall <- DataTrack(range = atac_gr, genome = "hg19",start =lim[1]-ext, end = lim[2]+ext,
        type = c("smooth"), chromosome = Chr, name = "Chromatin accessibility",groups=as.character(sample_anno$Condition),col=col, span=.05)

        #Plot track
        pdf(file.path(output.dir,"selectedDARs",paste0("locusPlot_Gene_",i,"_v1.pdf")), width = 7, height = 6)
        plotTracks(list(itrack ,grtrack,getrack,dars_track, trackall),from =lim[1]-ext, to = lim[2]+ext, 
        chromosome=Chr)
        dev.off()
        print(i)
    }
}

