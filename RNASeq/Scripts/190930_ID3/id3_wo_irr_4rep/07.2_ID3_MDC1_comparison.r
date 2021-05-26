#Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(knitr)
library(clusterProfiler)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)

#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)
rld <- readRDS(file =file.path(results.dir, "rld_b_replicate.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#Take a look at design and annotation of samples
design(dds)


alpha <- 0.05 #set FDR cutoff
lfc <- 0##set logfold2 cutoff

#Filter genes which are only expressed in 1 sample
length(dds)
idx_dds<- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
table(idx_dds)
#dds <- dds[idx,]
#dim(dds)
#heatmaps
#define padj cutoff for plotting
cutoff <- 0.05
genes2plot <- lapply(DEG_results_list, function(x){
  x <- x[which(x$padj < cutoff & abs(x$log2FoldChange)>0.5),]
  x <- rownames(x)
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]

#heatmap
mat_genes<- assay(rld)
#subset groups of interest
mat_genes_sub <- mat_genes[, rownames(colData(rld)[colData(rld)$genotype %in% c("ID3", "WT", "MDC1"),])]
annovst <- as.data.frame(colData(rld)[colData(rld)$genotype %in% c("ID3", "WT", "MDC1"),])[, c("genotype"), drop=F] 
color_genotype<-RColorBrewer::brewer.pal(n =4, name = 'Dark2')[c(2,4,1)]

names(color_genotype)<- as.character(unique(colData(rld)[colData(rld)$genotype %in% c("ID3", "WT", "MDC1"),]$genotype))
#color_treatment<-randomColor(length(unique(colData(rld)$irradiation_treatment)))
#names(color_treatment)<- as.character(unique(colData(rld)$irradiation_treatment))
color_replicate<-randomColor(length(unique(colData(rld)$Replicate)))
names(color_replicate)<- as.character(unique(colData(rld)$Replicate))
#anno_colors <- list(genotype=color_genotype, irradiation_treatment=color_treatment, Replicate=color_replicate)
anno_colors <- list(genotype=color_genotype,  Replicate=color_replicate)

heat<-list()

  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[["ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED"]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
annorow <- rownames(plot)
pheatmap(plot[, rownames(annovst)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors,
                       filename=file.path(PostDE.dir,
                       "ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED","Heat_DEG_woID3rescue.pdf"))
