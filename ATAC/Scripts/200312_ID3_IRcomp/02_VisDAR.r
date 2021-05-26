#directories
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")

#libraries
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ChIPseeker)
library(ChIPpeakAnno)
library(reshape2)
library(ggpubr)

#load data
DE_list <- readRDS(file.path(analysis.dir, "de_list.rds"))
counts<- readRDS(file.path(analysis.dir, "countsperpeak.rds"))
counts_raw<- readRDS(file.path(analysis.dir, "allPeaks_raw.rds"))
sample_anno<- readRDS(file.path(analysis.dir, "sample_anno.rds"))

#create common object
DE_list_complete<- readRDS(file.path(analysis.dir, "de_list_complete.rds"))
#subset significant
DE_list_complete_sign <- lapply(DE_list_complete, function(x){
    x <- x[x$FDR<0.05,]
    x
})
lapply(DE_list_complete_sign, function(x)length(x))
#subset fold change
DE_list_complete_sign_fc <- lapply(DE_list_complete_sign, function(x){
    x <- x[abs(x$Fold)>1,]
    x
})
lapply(DE_list_complete_sign_fc, function(x)length(x))


#Do Venn Diagramm of sign
#all important
temp <- DE_list_complete_sign[c("ID3-ko_Irr_vs_ID3-ko_unIrr","wt_Irr_vs_wt_unIrr")]
pdf(file.path(analysis.dir, "Venn_diagramm_TreatmentComp_onlySig.pdf"))
makeVennDiagram(temp, names(temp), TXDB=txdb)
dev.off()
#LPS vs medium
temp <- DE_list_complete_sign[c("ID3-ko_unIrr_vs_wt_unIrr", "ID3-ko_Irr_vs_wt_Irr")]
pdf(file.path(analysis.dir, "Venn_diagramm_GenotypeComp_onlySig.pdf"))
makeVennDiagram(temp, names(temp), TXDB=txdb)
dev.off()
#tg vs medium
temp <- DE_list_complete_sign[c("ID3-ko_Irr_vs_ID3-ko_unIrr","wt_Irr_vs_wt_unIrr","ID3-ko_unIrr_vs_wt_unIrr", "ID3-ko_Irr_vs_wt_Irr")]
pdf(file.path(analysis.dir, "Venn_diagramm_tgvswt_allCompn_onlySig.pdf"))
makeVennDiagram(temp, names(temp), TXDB=txdb)
dev.off()

#Do Venn Diagramm of sign and fc
#all important
temp <- DE_list_complete_sign_fc[c("ID3-ko_Irr_vs_ID3-ko_unIrr","wt_Irr_vs_wt_unIrr")]
pdf(file.path(analysis.dir, "Venn_diagramm_TreatmentComp_sigFC.pdf"))
makeVennDiagram(temp, names(temp), TXDB=txdb)
dev.off()
#LPS vs medium
temp <- DE_list_complete_sign_fc[c("ID3-ko_unIrr_vs_wt_unIrr", "ID3-ko_Irr_vs_wt_Irr")]
pdf(file.path(analysis.dir, "Venn_diagramm_GenotypeComp_sigFC.pdf"))
makeVennDiagram(temp, names(temp), TXDB=txdb)
dev.off()
#tg vs medium
temp <- DE_list_complete_sign[c("ID3-ko_Irr_vs_ID3-ko_unIrr","wt_Irr_vs_wt_unIrr","ID3-ko_unIrr_vs_wt_unIrr", "ID3-ko_Irr_vs_wt_Irr")]
pdf(file.path(analysis.dir, "Venn_diagramm_tgvswt_allCompn_sigFC.pdf"))
makeVennDiagram(temp, names(temp), TXDB=txdb)
dev.off()




#piechart of peaks
#for sign
DE_list_complete_sign<- lapply(DE_list_complete_sign, function(x){
    x$direction <- ifelse(x$Fold>0, "open", "closed")
    x
})
table <- lapply(DE_list_complete_sign, function(x){
    x <- as.data.frame(table(x$direction))
    x
})
labs<- lapply(table, function(x){
    x <- paste0(x$Var1,"\n(", round((x$Freq/sum(x$Freq)*100)),"%)\n", x$Freq)
    x
})
pie1 <- list()
for(i in names(table)){
pie1[[i]] <- ggpie(table[[i]], "Freq", fill="Var1", palette=c("#7AA6DCFF","#CD534CFF"), title=paste0("DAR analysis:\n", i),
    label=labs[[i]],lab.font = "white",  lab.pos = "in", #main="Differential accessible regions", 
    submain="p.adj < 0.05\n ") + rremove("legend")
}
#for sign +fc
DE_list_complete_sign_fc<- lapply(DE_list_complete_sign_fc, function(x){
    x$direction <- ifelse(x$Fold>0, "open", "closed")
    x
})
table <- lapply(DE_list_complete_sign_fc, function(x){
    x <- as.data.frame(table(x$direction))
    x
})
labs<- lapply(table, function(x){
    x <- paste0(x$Var1,"\n(", round((x$Freq/sum(x$Freq)*100)),"%)\n", x$Freq)
    x
})
pie2 <- list()
for(i in names(table)){
pie2[[i]] <- ggpie(table[[i]], "Freq", fill="Var1", palette=c("#7AA6DCFF","#CD534CFF"), label=labs[[i]],lab.font = "white",  lab.pos = "in", #main="Differential accessible regions", 
submain="\np.adj < 0.05\nabs(log2Fold) > 1") + rremove("legend")
}
#plot both pies
dir.create(file.path(analysis.dir, "Pies"))
for (i in names(table)){
    fig <- ggarrange(pie1[[i]], pie2[[i]], ncol = 2, nrow = 1)

    pdf(file.path(analysis.dir,"Pies",paste0("directionofDiffPeaks_", i,".pdf" )), height= 7, width= 7)
    #annotate_figure(fig, top= text_grob("Differential accessible regions", face="bold", size=14))
    print(fig)
    dev.off()
}

#distance to TSS
hist1 <- lapply(DE_list_complete_sign, function(x){
    x <- mcols(x[,"distanceToTSS"])
    x$distanceToTSS <- x$distanceToTSS/1000
    x <- as.data.frame(x)
    x <- gghistogram(x, x="distanceToTSS", fill="lightgray", bins=40, rug = T, xlab="Distance from TSS (kb)", ylab="# of DAR",
    main= paste0("\np.adj < 0.05"))  + geom_vline(xintercept=0, linetype = 1, color="red", size=0.5 )
    x
})
hist2 <- lapply(DE_list_complete_sign_fc, function(x){
    x <- mcols(x[,"distanceToTSS"])
    x$distanceToTSS <- x$distanceToTSS/1000
    x <- as.data.frame(x)
    x <- gghistogram(x, x="distanceToTSS", fill="lightgray", bins=40, rug = T, xlab="Distance from TSS (kb)", ylab="# of DAR",
    main= paste0("p.adj < 0.05\nabs(log2Fold) > 1"))  + geom_vline(xintercept=0, linetype = 1, color="red", size=0.5 )
    x
})
dir.create(file.path(analysis.dir, "Distance"))
for(i in names(DE_list_complete_sign_fc)){
    fig <- ggarrange(hist1[[i]], hist2[[i]], ncol = 2, nrow = 1)
    fig <- annotate_figure(fig,
               top = text_grob(paste0("Distance to TSS: ",i),  face = "bold", size = 9))
    pdf(file.path(analysis.dir,"Distance",paste0("distance_to_TSS_", i,".pdf" )), height= 3.5, width= 7)
    print(fig)
    dev.off()
}

#volcano plot
plot <- lapply(DE_list_complete, function(x){
    x <- as.data.frame(x)
    x$col <- "black"
    #x[which(x$FDR < 0.05 & x$Fold>1),"col" ]<- "blue"
    #x[which(x$FDR < 0.05 & x$Fold<(-1)),"col" ]<- "red"
    x[which(x$FDR < 0.05 & x$Fold>1),"col" ]<- "red"
    x[which(x$FDR < 0.05 & x$Fold<(-1)),"col" ]<- "blue"
    x$logpadj <- -(log10(x$FDR))

    x$peakNr <- rownames(x)
    x <- x[order(x$FDR, decreasing=FALSE), ]
    names <- c(head(x[x$Fold>0,]$peakNr, 10))
    names <- c(names, head(x[x$Fold<0,]$peakNr, 10))
    x$label <- NA
    x$label <- ifelse(x$peakNr %in% names, x$SYMBOL, NA)
    x
})

dir.create(file.path(analysis.dir, "Volcanos2"))
dir.create(file.path(analysis.dir, "MAs2"))
for(i in names(plot)){
    pdf(file.path(analysis.dir,"Volcanos2",paste0("Volcano_nolabel_noFC_", i,".pdf" )), height= 3.5, width= 3.5)
    print(ggscatter(plot[[i]], x ="Fold", y = "logpadj",xlab="Chromatin accessibility (log2 fold)",#xscale="log10" ,
          ylab="- log10(P.adjusted)", title=i, 
          color = "col" ,shape = 16, size = 1, palette=c("grey","#7AA6DCFF","#CD534CFF"), # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance",# Customize reg. line
          #label="label", repel = T, font.label = c( "black")#, legend = "right",
            )  +geom_hline(yintercept=1.30103, linetype = 2) + rremove("legend")+ geom_vline(xintercept=c(-1, 1), linetype = 2))#+ geom_vline(xintercept=c(-1,1), linetype = 2)
    dev.off()
    png(file.path(analysis.dir,"Volcanos2",paste0("Volcano_nolabel_noFC_", i,".png" )))
    print(ggscatter(plot[[i]], x ="Fold", y = "logpadj",xlab="Chromatin accessibility (log2 fold)",#xscale="log10" ,
          ylab="- log10(P.adjusted)", title=i, 
          color = "col" ,shape = 16, size = 1, palette=c("grey","#7AA6DCFF","#CD534CFF"), # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance",# Customize reg. line
          #label="label", repel = T, font.label = c( "black")#, legend = "right",
            )  +geom_hline(yintercept=1.30103, linetype = 2) + rremove("legend")+ geom_vline(xintercept=c(-1, 1), linetype = 2))#+ geom_vline(xintercept=c(-1,1), linetype = 2)
    dev.off()

    pdf(file.path(analysis.dir,"MAs2",paste0("MAs2_noFC_", i,".pdf" )), height= 3.5, width= 3.5)
    print(ggscatter(plot[[i]], x ="Conc", y =  "Fold",xlab="Mean of log2 read counts",xscale="log10" ,
          ylab="Chromatin accessibility (log2 fold)", #xlim=c(0, max(plot[[i]]$Conc)),
          color = "col",palette=c("grey","#7AA6DCFF","#CD534CFF"), ,shape = 16, size = 1, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), #legend.title="Significance",# Customize reg. line
        ) +geom_hline(yintercept=c(-1,1), linetype = 2) + rremove("legend"))
    dev.off()
}

#overlap with gene regulatory regions
#load regions
ctcf <- read.table("c010-datasets/Internal/COPD/enrichment_databases//mm10_genomic_regions/ensembl_regulatory/regions/multi-cell_fe-ctcf-binding-site_st-undef.bed", sep="\t")
enhancer <- read.table("c010-datasets/Internal/COPD/enrichment_databases//mm10_genomic_regions/ensembl_regulatory/regions/multi-cell_fe-enhancer_st-undef.bed", sep="\t")
openchrom <- read.table("c010-datasets/Internal/COPD/enrichment_databases//mm10_genomic_regions/ensembl_regulatory/regions/multi-cell_fe-open-chromatin-region_st-undef.bed", sep="\t")
promoterflank <- read.table("c010-datasets/Internal/COPD/enrichment_databases//mm10_genomic_regions/ensembl_regulatory/regions/multi-cell_fe-promoter-flanking-region_st-undef.bed", sep="\t")
promoter <- read.table("c010-datasets/Internal/COPD/enrichment_databases//mm10_genomic_regions/ensembl_regulatory/regions/multi-cell_fe-promoter_st-undef.bed", sep="\t")
tfbind <-read.table("c010-datasets/Internal/COPD/enrichment_databases//mm10_genomic_regions/ensembl_regulatory/regions/multi-cell_fe-tf-binding-site_st-undef.bed", sep="\t")
RegulatoryRegions<- list(ctcf, enhancer, openchrom,promoterflank, promoter,tfbind)
for (i in 1:length(RegulatoryRegions)) {
  RegulatoryRegions[[i]]<- RegulatoryRegions[[i]][,1:4]
  colnames(RegulatoryRegions[[i]])<-c("seqnames", "start", "end", "class")
  RegulatoryRegions[[i]]$strand <- "*"
  RegulatoryRegions[[i]]<- makeGRangesFromDataFrame(RegulatoryRegions[[i]], keep.extra.columns = TRUE)

}
names(RegulatoryRegions)<- c("ctcf binding site", "enhancer", "open chromatin","promoter flanking regions", "promoter","tf binding site")

overlaps_pos <- list()
Total_pos <- list()
overlaps_neg <- list()
Total_neg <- list()
sum_regions<- list()
sum_regions_melt <- list()
dir.create(file.path(analysis.dir,"GeneReg_Anno"))
for(i in names(DE_list_complete_sign)){
    overlaps_pos[[i]] <- lapply(RegulatoryRegions, function(x) findOverlaps(makeGRangesFromDataFrame(DE_list_complete_sign[[i]][DE_list_complete_sign[[i]]$Fold>0,]), x))
    Total_pos[[i]] <- (unlist(lapply(overlaps_pos[[i]], function(x) length(unique(queryHits(x)))))/length(makeGRangesFromDataFrame(DE_list_complete_sign[[i]][DE_list_complete_sign[[i]]$Fold>0,])))*100

    overlaps_neg[[i]] <- lapply(RegulatoryRegions, function(x) findOverlaps(makeGRangesFromDataFrame(DE_list_complete_sign[[i]][DE_list_complete_sign[[i]]$Fold<0,]), x))
    Total_neg[[i]] <- (unlist(lapply(overlaps_pos[[i]], function(x) length(unique(queryHits(x)))))/length(makeGRangesFromDataFrame(DE_list_complete_sign[[i]][DE_list_complete_sign[[i]]$Fold<0,])))*100

    sum_regions[[i]] <- data.frame("open DARs"=Total_pos[[i]], "closed DARs"= Total_neg[[i]], Regions=names(RegulatoryRegions))
    sum_regions_melt[[i]] <- melt(sum_regions[[i]])
    sum_regions_melt[[i]] <- melt(sum_regions[[i]])

    pdf(file.path(analysis.dir,"GeneReg_Anno",paste0("Anno_", i,".pdf" )), height= 5, width= 5)
    print(ggbarplot(sum_regions_melt[[i]], x="Regions", y="value", fill="variable", palette=rev(c("#7AA6DCFF","#CD534CFF")), position = position_dodge(0.9), ylim=c(0,100),ylab= "% of DARs")+ rotate_x_text(45)+coord_flip()+rremove("legend.title")+rremove("y.title"))
    dev.off()
}

#for sign +fc
overlaps_pos <- list()
Total_pos <- list()
overlaps_neg <- list()
Total_neg <- list()
sum_regions<- list()
sum_regions_melt <- list()
dir.create(file.path(analysis.dir,"GeneReg_Anno"))
for(i in names(DE_list_complete_sign_fc)){
    overlaps_pos[[i]] <- lapply(RegulatoryRegions, function(x) findOverlaps(makeGRangesFromDataFrame(DE_list_complete_sign_fc[[i]][DE_list_complete_sign_fc[[i]]$Fold>0,]), x))
    Total_pos[[i]] <- (unlist(lapply(overlaps_pos[[i]], function(x) length(unique(queryHits(x)))))/length(makeGRangesFromDataFrame(DE_list_complete_sign_fc[[i]][DE_list_complete_sign_fc[[i]]$Fold>0,])))*100

    overlaps_neg[[i]] <- lapply(RegulatoryRegions, function(x) findOverlaps(makeGRangesFromDataFrame(DE_list_complete_sign_fc[[i]][DE_list_complete_sign_fc[[i]]$Fold<0,]), x))
    Total_neg[[i]] <- (unlist(lapply(overlaps_pos[[i]], function(x) length(unique(queryHits(x)))))/length(makeGRangesFromDataFrame(DE_list_complete_sign_fc[[i]][DE_list_complete_sign_fc[[i]]$Fold<0,])))*100

    sum_regions[[i]] <- data.frame("open DARs"=Total_pos[[i]], "closed DARs"= Total_neg[[i]], Regions=names(RegulatoryRegions))
    sum_regions_melt[[i]] <- melt(sum_regions[[i]])
    sum_regions_melt[[i]] <- melt(sum_regions[[i]])

    pdf(file.path(analysis.dir,"GeneReg_Anno",paste0("Anno_Sign_FC_", i,".pdf" )), height= 5, width= 5)
    print(ggbarplot(sum_regions_melt[[i]], x="Regions", y="value", fill="variable", palette=rev(c("#7AA6DCFF","#CD534CFF")), position = position_dodge(0.9), ylim=c(0,100),ylab= "% of DARs")+ rotate_x_text(45)+coord_flip()+rremove("legend.title")+rremove("y.title"))
    dev.off()
}

 
#plot heatmaps of DARs
#get contrast
contrasts <- readRDS(file.path(analysis.dir, "contrasts"))
rownames(contrasts)<- contrasts$comparison
#get annotation
annorld<- sample_anno[,c("Factor","Treatment")]
library(pheatmap)
color_genotype<-c("#D95F02","#1B9E77")
names(color_genotype)<- as.character(unique(annorld$Factor))
color_treatment<-c("#7FC97F", "#BEAED4")
names(color_treatment)<- as.character(unique(as.character(annorld$Treatment)))
anno_colors <- list(Factor=color_genotype,  Treatment=color_treatment)

heat <- list()
heat_sub <-list()
for (i in names(DE_list_complete_sign_fc)[5:6]){
    plot <- mcols(DE_list_complete_sign_fc[[i]])[, rownames(sample_anno)]
    heat[[i]]<- pheatmap(plot, 
    scale="row", show_colnames=F, annotation_col=annorld,show_rownames=F,
        annotation_colors=anno_colors, clustering_distance_rows="correlation") 
    pdf(file.path(analysis.dir,"heatmaps",paste0("Heat_correlation_",i,"_DAR_SigFC_allGroups.pdf")))
    print(heat[[i]])
    dev.off()
    plot_sub <- plot[, c(rownames(sample_anno[sample_anno$Condition==as.character(contrasts[i,2]),]), rownames(sample_anno[sample_anno$Condition==as.character(contrasts[i,3]),]))]
    heat_sub[[i]]<- pheatmap(plot_sub, 
    scale="row", show_colnames=F, annotation_col=annorld,show_rownames=F,
        annotation_colors=anno_colors, clustering_distance_rows="correlation") 
    pdf(file.path(analysis.dir,"heatmaps",paste0("Heat_correlation_",i,"_DAR_SigFC_subgroups.pdf")))
    print(heat_sub[[i]])
    dev.off()

    plot <- mcols(DE_list_complete_sign[[i]])[, rownames(sample_anno)]
    heat<- pheatmap(plot, 
    scale="row", show_colnames=F, annotation_col=annorld,show_rownames=F,
        annotation_colors=anno_colors, clustering_distance_rows="correlation") 
    pdf(file.path(analysis.dir,"heatmaps",paste0("Heat_correlation_",i,"_DAR_Sig_allGroups.pdf")))
    print(heat)
    dev.off()
    plot_sub <- plot[, c(rownames(sample_anno[sample_anno$Condition==as.character(contrasts[i,2]),]), rownames(sample_anno[sample_anno$Condition==as.character(contrasts[i,3]),]))]
    heat_sub<- pheatmap(plot_sub, 
    scale="row", show_colnames=F, annotation_col=annorld,show_rownames=F,
        annotation_colors=anno_colors, clustering_distance_rows="correlation") 
    pdf(file.path(analysis.dir,"heatmaps",paste0("Heat_correlation_",i,"_DAR_Sig_subgroups.pdf")))
    print(heat_sub)
    dev.off()
}

#do all of them combined
#sign and fc
DE_list_df <-lapply(DE_list_complete_sign_fc, function(x) as.data.frame(x))
DE_list_df <-lapply(DE_list_df,function(x)x[,c("seqnames","start", "end", "width", "strand", rownames(sample_anno))])
DE_list_df <- do.call("rbind", DE_list_df)
DE_list_df <- unique(DE_list_df)
DE_list_gr <- makeGRangesFromDataFrame(DE_list_df, keep.extra.columns=TRUE)
length(DE_list_gr)
plot <- mcols(DE_list_gr)

heat<- pheatmap(plot, 
    scale="row", show_colnames=F, annotation_col=annorld,show_rownames=F,
    annotation_colors=anno_colors, clustering_distance_rows="correlation") 
pdf(file.path(analysis.dir,"heatmaps",paste0("Heat_correlation_","All","_DAR_SigFC_allGroups.pdf")))
print(heat)
dev.off()

#only sig
DE_list_df <-lapply(DE_list_complete_sign, function(x) as.data.frame(x))
DE_list_df <-lapply(DE_list_df,function(x)x[,c("seqnames","start", "end", "width", "strand", rownames(sample_anno))])
DE_list_df <- do.call("rbind", DE_list_df)
DE_list_df <- unique(DE_list_df)
DE_list_gr <- makeGRangesFromDataFrame(DE_list_df, keep.extra.columns=TRUE)
length(DE_list_gr)
plot <- mcols(DE_list_gr)

heat<- pheatmap(plot, 
    scale="row", show_colnames=F, annotation_col=annorld,show_rownames=F,
    annotation_colors=anno_colors, clustering_distance_rows="correlation") 
pdf(file.path(analysis.dir,"heatmaps",paste0("Heat_correlation_","All","_DAR_Sig_allGroups.pdf")))
print(heat)
dev.off()