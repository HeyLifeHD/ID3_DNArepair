#directories
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")

#Packages
library(DESeq2)
library(RnBeads)
library(org.Mm.eg.db)
library(ChIPseeker)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(topGO)
library(gage)
library(biomaRt)
library(GenomicRanges)
library(rafalib)
library(rtracklayer)
library(pheatmap)
library(ggpubr)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ggpubr)
library(Rsubread)
library("EnsDb.Hsapiens.v75")
ENS<-EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"
GENE = genes(ENS)
Promoter <- promoters(GENE, upstream = 3000, downstream = 3000)


#ATAC
db.anno_all<- readRDS(file.path(analysis.dir, "DiffBind_DB_anno_all.rds"))


#load data
sample_anno<- readRDS(file.path(analysis.dir, "sample_anno.rds"))
counts_raw<- readRDS(file.path(analysis.dir, "allPeaks_raw.rds"))

#with interactive term
dds <- DESeqDataSetFromMatrix(countData = mcols(counts_raw), 
                              colData = sample_anno, 
                              design = ~  Replicate + Factor+ Treatment, rowRanges=counts_raw)
#for Genotype comparison
dds <- estimateSizeFactors(dds)
dim(dds)
#Running the differential expression 
dds <- DESeq(dds)
saveRDS(dds,file = file.path(analysis.dir,"dds_group.rds"))
dds <- readRDS( file.path(analysis.dir,"dds_group.rds"))
#extract transformed values
rld <- rlog(dds, blind=FALSE)
saveRDS(rld,file = file.path(analysis.dir,"rld.rds"))
rld<- readRDS(file.path(analysis.dir,"rld.rds"))


#plot PC1 and 2
pcaDatarld<- plotPCA(rld, intgroup =  colnames(colData(rld)), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))
dir.create(file.path(analysis.dir, "Explor"))
pdf(file.path(analysis.dir, "Explor", "PCA12_allGenes.pdf"), height = 4, width = 4)
ggscatter(pcaDatarld, x="PC1", y="PC2",size=3,
          color = "Factor", shape = "Treatment", legend = "right",
          ellipse = F ,mean.point = FALSE,palette = "jco",#label="SampleID",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC1: ", percentvarrld[1], "% variance")), ylab=(paste0("PC2: ", percentvarrld[2], "% variance")) )
dev.off()

#function to plot PC2 and 3
plotPCA.alt <- function (object, intgroup = "condition", ntop = Inf, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + 
    geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}

pcaDatarld.alt<- plotPCA.alt(rld, intgroup = colnames(sample_anno), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))


pdf(file.path(analysis.dir, "Explor","PCA23_allGenes.pdf"), height = 4, width = 4)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",size=3,
          color = "Factor", shape = "Treatment", legend = "right",
          ellipse = F ,mean.point = FALSE,palette="jco",
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()


#Extract genes behind PC1
ir.pca <- prcomp(t(assay(rld)),
                 center = T,
                 scale. = T) 

summary(ir.pca)
x <- ir.pca$x
x<- as.data.frame(cbind(x , colData(rld)))

color_genotype<-RColorBrewer::brewer.pal(n = length(unique(annorld$Factor)), name = 'Dark2')

pdf(file.path(analysis.dir, "Explor", "PCA12_allGenes_2.pdf"), height=4, width=5)
ggscatter(x, x="PC1", y="PC2",
          color = "Factor", shape = "Treatment",
          ellipse = F , mean.point = FALSE, palette = color_genotype,
          star.plot = F, xlab=(paste0("PC1: ", round(summary(ir.pca)$importance[2,1]*100,2), "% variance")), ylab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance")))+theme(legend.position="right")
dev.off()
pdf(file.path(analysis.dir, "Explor", "PCA23_allGenes_2.pdf"), height=4, width=5)
ggscatter(x, x="PC2", y="PC3",
          color = "Factor", shape = "Treatment",
          ellipse = F , mean.point = FALSE, palette =color_genotype,
          star.plot = F, xlab=(paste0("PC2: ", round(summary(ir.pca)$importance[2,2]*100,2), "% variance")), ylab=(paste0("PC3: ", round(summary(ir.pca)$importance[2,3]*100,2), "% variance")))+theme(legend.position="right")
dev.off()

#prop <- as.data.frame(ir.pca$x)
#top200both <- rownames(prop[head(order(abs(prop$PC1), decreasing=TRUE), 200),])
#as.factor(top200both)
#top200pos <- rownames(prop[head(order((prop$PC3), decreasing=TRUE), 200),])
#as.factor(top200pos)

#Sample Clustering
library(dendextend)

dvst <- dist(t(assay(rld)))
hvst <- hclust(dvst)
dend<- dist(t(assay(rld))) %>% hclust %>% as.dendrogram 
#dend<-dend %>% set("labels_col", c( "#EFC00099", "#EFC00099", "#EFC00099",  "#86868699","#86868699", "#86868699", "#86868699","#0073C299","#0073C299","#0073C299")) %>% 
#set("labels_cex", .5) 
 dend %>% plot
pdf(file.path(analysis.dir, "Explor", "Clustering_rld.pdf"), height = 7, width = 4)
#myplclust(hvst, labels=colData(rld)$SampleID, lab.col= as.numeric(as.factor(colData(rld)$Celltype),   tip.color  = c("#0073C299", "#EFC00099", "#86868699")), main="", cex=0.5)
 dend %>% plot
dev.off()

#Gene clustering: Heatmaps
#500 most variable Genes
topVarGenesvst<- head(order(rowVars(assay(rld)), decreasing=TRUE),20000)

annorld<- sample_anno[,c("Factor","Treatment")]
library(pheatmap)
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(annorld$Factor)), name = 'Dark2')
names(color_genotype)<- as.character(unique(annorld$Factor))
color_treatment<-c("#7FC97F", "#BEAED4")
names(color_treatment)<- as.character(unique(annorld$Treatment))
anno_colors <- list(Factor=color_genotype,  Treatment=color_treatment)

pheatmap(assay(rld)[topVarGenesvst,], annotation_col=annorld ,scale="row", annotation_colors=anno_colors, 
show_colnames=F, show_rownames=F,file=file.path(analysis.dir, "Explor","Heatmap_20000mv_Rld.pdf"))
pheatmap(assay(rld)[topVarGenesvst,], annotation_col=annorld ,scale="none", annotation_colors=anno_colors, 
show_colnames=F, show_rownames=F,file=file.path(analysis.dir, "Explor","Heatmap_20000mv_Rld_noScale.pdf"))
