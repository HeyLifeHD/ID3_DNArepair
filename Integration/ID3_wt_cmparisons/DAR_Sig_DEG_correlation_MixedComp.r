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
output.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/Integration"
dir.create(output.dir,recursive=TRUE)

#libraries
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ChIPseeker)
library(ChIPpeakAnno)
library(reshape2)
library(ggpubr)
library(dplyr)
library(LSD)

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

#edit deg table
DEG_results_list <- lapply(DEG_results_list, function(x){
    x$ENSEMBL <- x$ensembl
    x
})

#merge data
merged <- list()
temp <- inner_join(as.data.frame(DAR_list_complete_sign[["ID3-ko_unIrr_vs_wt_unIrr"]]),as.data.frame(DEG_results_list[["ID3-ko_Irr_vs_wt_Irr"]]), by="ENSEMBL")
temp$Col <- NA
temp$Col <- ifelse(temp$padj < 0.05 & temp$log2FoldChange>0.5, "sign. upregulated",ifelse(temp$padj < 0.1 & temp$log2FoldChange<(-0.5), "sign. downregulated", "no sign. change"))
temp$Col <- ifelse(is.na(temp$Col), "no sign. change", temp$Col)
temp$Col<- as.factor(temp$Col)
temp$Col2 <- NA
temp$Col2 <- ifelse(temp$padj < 0.05 & temp$log2FoldChange>0.5, "sign. upregulated",ifelse(temp$padj < 0.1 & temp$log2FoldChange<(-0.5), "sign. downregulated", "no sign. change"))
temp$Col2 <- ifelse(is.na(temp$Col2), "no sign. change", temp$Col2)
temp$Col2<- as.factor(temp$Col2)
temp$Plot<-NA
temp[temp$Col2!="no sign. change",]$Plot <- temp[temp$Col2!="no sign. change",]$symbol
#ind <-  which(temp$SYMBOL %in% c("Pdpn","Ogn", "Sfrp4",  "Col3a1", "Fbn1", "Cxcl12", "Gsn", "Acta2", "Saa3", "Saa3", "Il6", "Spp1", "S100a4", "Cd74"))
#temp[ind, "Plot"] <-temp[ind, "SYMBOL"]
merged[["DAR_UIrrComp_DEG_IrrComp"]] <- temp


lapply(merged, function(x)dim(x))

#write tables
for(i in names(merged)){
    dir.create(file.path(output.dir,i))
   write.table(merged[[i]], file.path(output.dir, i , paste0("DAR_RNA_merged_DAR",".txt")),quote=F, sep="\t", row.names=FALSE)
}

#all DARs
merged_sub <- merged
#create folders
#plot
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged_sub)){
#name <- gsub("_", " ", i)
name <- i
g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_allDARs_allGenes_nolabel.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

pdf(file.path(output.dir, i , "DAR_RNA_Heat_allDARs_allGenes_nolabel.pdf"), width=6,height=6)
heatscatter(merged_sub[[i]]$log2FoldChange, merged_sub[[i]]$Fold,
    cexplot = .2,
    xlab=paste0("\nGene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    cor=TRUE,
    add.contour=FALSE) +
    abline(lm(Fold ~0 +log2FoldChange, merged_sub[[i]]), 
        col="black", lwd=2, lty=1) +
    abline(h=0, lwd=1, lty=2, col="black")  +
    abline(v=-0.1,, lwd=1, lty=2, col="black")
dev.off()

print(i)
}

#subset based on promoter  (all) location
merged_sub<- lapply(merged, function(x){
    x <- x[grep("Promoter", x$annotation), ]
    #x <- x[which(x$padj < 0.05),]
    #x <- x[which(x$annotation =="Promoter (<=1kb)"),]
    #x <-  x[which(x$SYMBOL %in% c("Pdpn","Ogn", "Sfrp4",  "Col3a1", "Fbn1", "Cxcl12", "Gsn", "Acta2", "Saa3", "Saa3", "Il6", "Spp1", "S100a4", "Cd74")),]
    x
})
lapply(merged_sub, function(x)dim(x))
#create folders
#plot
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged_sub)){
#name <- gsub("_", " ", i)
name <- i
g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_PromoterDARs_allGenes_nolabel.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

pdf(file.path(output.dir, i , "DAR_RNA_Heat_PromoterDARs_allGenes_nolabel.pdf"), width=6,height=6)
heatscatter(merged_sub[[i]]$log2FoldChange, merged_sub[[i]]$Fold,
    cexplot = .2,
    xlab=paste0("\nGene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    cor=TRUE,
    add.contour=FALSE) +
    abline(lm(Fold ~0 +log2FoldChange, merged_sub[[i]]), 
        col="black", lwd=2, lty=1) +
    abline(h=0, lwd=1, lty=2, col="black")  +
    abline(v=-0.1,, lwd=1, lty=2, col="black")
dev.off()

print(i)
}

#same for only sign DEGs
#all DARs
merged_sub<- lapply(merged, function(x){
    #x <- x[grep("Promoter", x$annotation), ]
    x <- x[which(x$padj < 0.05),]
    #x <- x[which(x$annotation =="Promoter (<=1kb)"),]
    #x <-  x[which(x$SYMBOL %in% c("Pdpn","Ogn", "Sfrp4",  "Col3a1", "Fbn1", "Cxcl12", "Gsn", "Acta2", "Saa3", "Saa3", "Il6", "Spp1", "S100a4", "Cd74")),]
    x
})#create folders
lapply(merged_sub, function(x)dim(x))

#plot
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged_sub)){
#name <- gsub("_", " ", i)
name <- i
g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_allDARs_allSigGenes_nolabel.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

pdf(file.path(output.dir, i , "DAR_RNA_Heat_allDARs_allSigGenes_nolabel.pdf"), width=6,height=6)
heatscatter(merged_sub[[i]]$log2FoldChange, merged_sub[[i]]$Fold,
    cexplot = .2,
    xlab=paste0("\nGene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    cor=TRUE,
    add.contour=FALSE) +
    abline(lm(Fold ~0 +log2FoldChange, merged_sub[[i]]), 
        col="black", lwd=2, lty=1) +
    abline(h=0, lwd=1, lty=2, col="black")  +
    abline(v=-0.1,, lwd=1, lty=2, col="black")
dev.off()

print(i)
}

#subset based on promoter  (all) location
merged_sub<- lapply(merged, function(x){
    x <- x[grep("Promoter", x$annotation), ]
    x <- x[which(x$padj < 0.05),]
    #x <- x[which(x$annotation =="Promoter (<=1kb)"),]
    #x <-  x[which(x$SYMBOL %in% c("Pdpn","Ogn", "Sfrp4",  "Col3a1", "Fbn1", "Cxcl12", "Gsn", "Acta2", "Saa3", "Saa3", "Il6", "Spp1", "S100a4", "Cd74")),]
    x
})
lapply(merged_sub, function(x)dim(x))
#create folders
#plot
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged_sub)){
#name <- gsub("_", " ", i)
name <- i
g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_PromoterDARs_allSigGenes_nolabel.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    label="Plot",  font.label = c(5, "plain"), repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_PromoterDARs_allSigGenes_label.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

pdf(file.path(output.dir, i , "DAR_RNA_Heat_PromoterDARs_allSigGenes_nolabel.pdf"), width=6,height=6)
heatscatter(merged_sub[[i]]$log2FoldChange, merged_sub[[i]]$Fold,
    cexplot = .2,
    xlab=paste0("\nGene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    cor=TRUE,
    add.contour=FALSE) +
    abline(lm(Fold ~0 +log2FoldChange, merged_sub[[i]]), 
        col="black", lwd=2, lty=1) +
    abline(h=0, lwd=1, lty=2, col="black")  +
    abline(v=-0.1,, lwd=1, lty=2, col="black")
dev.off()

print(i)
}



#same for only sign and fcDEGs
#all DARs
merged_sub<- lapply(merged, function(x){
    #x <- x[grep("Promoter", x$annotation), ]
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange) >0.5),]
    #x <- x[which(x$annotation =="Promoter (<=1kb)"),]
    #x <-  x[which(x$SYMBOL %in% c("Pdpn","Ogn", "Sfrp4",  "Col3a1", "Fbn1", "Cxcl12", "Gsn", "Acta2", "Saa3", "Saa3", "Il6", "Spp1", "S100a4", "Cd74")),]
    x
})#create folders
lapply(merged_sub, function(x)dim(x))

#plot
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged_sub)){
#name <- gsub("_", " ", i)
name <- i
g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_allDARs_allSigFCGenes_nolabel.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

pdf(file.path(output.dir, i , "DAR_RNA_Heat_allDARs_allSigFCGenes_nolabel.pdf"), width=6,height=6)
heatscatter(merged_sub[[i]]$log2FoldChange, merged_sub[[i]]$Fold,
    cexplot = .2,
    xlab=paste0("\nGene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    cor=TRUE,
    add.contour=FALSE) +
    abline(lm(Fold ~0 +log2FoldChange, merged_sub[[i]]), 
        col="black", lwd=2, lty=1) +
    abline(h=0, lwd=1, lty=2, col="black")  +
    abline(v=-0.1,, lwd=1, lty=2, col="black")
dev.off()

print(i)
}

#subset based on promoter  (all) location
merged_sub<- lapply(merged, function(x){
    x <- x[grep("Promoter", x$annotation), ]
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange) >0.5),]
    #x <- x[which(x$annotation =="Promoter (<=1kb)"),]
    #x <-  x[which(x$SYMBOL %in% c("Pdpn","Ogn", "Sfrp4",  "Col3a1", "Fbn1", "Cxcl12", "Gsn", "Acta2", "Saa3", "Saa3", "Il6", "Spp1", "S100a4", "Cd74")),]
    x
})
lapply(merged_sub, function(x)dim(x))
#create folders
#plot
col <- c("gray","#7AA6DCFF","#CD534CFF")
g <- list()
for(i in names(merged_sub)){
#name <- gsub("_", " ", i)
name <- i
g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    #label="Plot", repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_PromoterDARs_allSigFCGenes_nolabel.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

g[[i]] <- ggscatter(merged_sub[[i]], x = "log2FoldChange", y = "Fold",
    xlab=paste0("Gene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    color = "Col" ,
    shape = 16, size = 1,
    palette=col, # Points color, shape and size
    add = "reg.line",  # Add regressin line
    conf.int = T, # Add confidence interval
    add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
    label="Plot",  font.label = c(5, "plain"),repel = TRUE, #legend.title="P.adjusted",legend = "right",
    cor.coef = F, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)+  stat_cor( label.y.npc = "top", label.x.npc = "left")   + 
    geom_vline(xintercept=0, linetype = 2) +
    geom_hline(yintercept=0, linetype = 2) + rremove("legend")

dir.create(file.path(output.dir,i))
pdf(file.path(output.dir, i , "DAR_RNA_PromoterDARs_allSigFCGenes_label.pdf"), width=4,height=4)
    print(g[[i]])
dev.off()

pdf(file.path(output.dir, i , "DAR_RNA_Heat_PromoterDARs_allSigFCGenes_nolabel.pdf"), width=6,height=6)
heatscatter(merged_sub[[i]]$log2FoldChange, merged_sub[[i]]$Fold,
    cexplot = .2,
    xlab=paste0("\nGene expression change\n[log2foldchange]","\n",name),
    ylab=paste0("Chromatin accessibility change\n[log2foldchange]","\n",name),
    cor=TRUE,
    add.contour=FALSE) +
    abline(lm(Fold ~0 +log2FoldChange, merged_sub[[i]]), 
        col="black", lwd=2, lty=1) +
    abline(h=0, lwd=1, lty=2, col="black")  +
    abline(v=-0.1,, lwd=1, lty=2, col="black")
dev.off()

print(i)
}








