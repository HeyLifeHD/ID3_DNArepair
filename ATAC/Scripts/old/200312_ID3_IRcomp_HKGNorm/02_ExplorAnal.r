
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
library(DESeq2)
library(ggpubr)
library(rafalib)
#load data
counts_raw <- readRDS(file.path(analysis.dir,  "data","allPeaks_raw.rds"))
sample_anno <- readRDS(file.path(analysis.dir, "data","sample_anno.rds" ))
dds <- readRDS(file = file.path(analysis.dir, "DESEQ","dds.rds"))



#Extracting transformed values
rld <- rlog(dds)
anno <- colData(rld)
saveRDS(rld,file = file.path(analysis.dir, "DESEQ","rld.rds"))
rld <- readRDS(file.path(analysis.dir, "DESEQ","rld.rds"))

#function to plot PC1 and 2
pcaDatarld<- plotPCA(rld, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))
col <- c("#D95F02","#1B9E77")
dir.create(file.path(analysis.dir,"Explor"))
pdf(file.path(analysis.dir,"Explor", "PCA12.pdf"), height = 5, width = 6)
ggscatter(pcaDatarld, x="PC1", y="PC2",
            size=5,label="Replicate",
            color = "Factor", shape = "Treatment",
            repel=TRUE,
            legend = "right",palette=col,
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarrld[1], "% variance")), ylab=(paste0("PC2: ", percentvarrld[2], "% variance")) )
dev.off()

#function
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
#function to plot PC2 and 3
pcaDatarld.alt<- plotPCA.alt(rld, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
#pcaDatarld.alt$Replicate <- as.factor(pcaDatarld.alt$Replicate )
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))

pdf(file.path(analysis.dir,"Explor", "PCA23.pdf"), height = 5, width = 6)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",            
            size=5,label="Replicate",
            color = "Factor", shape = "Treatment",
            repel=TRUE,
            legend = "right",palette=col,
            ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#Sample Clustering
drld <- dist(t(assay(rld)))
hrld<- hclust(drld)

pdf(file.path(analysis.dir,"Explor", "Clustering_rld.pdf"), height = 7, width = 7)
myplclust(hrld, labels= rownames(anno), lab.col=as.numeric((anno$Condition)), main="", cex=0.2)
dev.off()


