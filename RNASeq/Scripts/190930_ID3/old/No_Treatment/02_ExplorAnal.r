#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(Rtsne)
#folder

base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/No_Irr"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)

#functions
plotPCA.alt <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
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

#Extracting transformed values
vst_genes <- vst(dds)
anno <- colData(vst_genes)
anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment)
colData(vst_genes)<- anno 
saveRDS(vst_genes,file = file.path(results.dir,"vst_genes.rds"))

#function to plot PC1 and 2
pcaDatarld<- plotPCA(vst_genes, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_genes.pdf"), height = 10, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",size=3,
          color = "group.1", shape = "Replicate", legend = "right",
          ellipse = F ,mean.point = FALSE,repel=TRUE,
          star.plot = F,  xlab=(paste0("PC1: ", percentvarrld[1], "% variance")), ylab=(paste0("PC2: ", percentvarrld[2], "% variance")) )
dev.off()

#function to plot PC2 and 3
pcaDatarld.alt<- plotPCA.alt(vst_genes, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
pcaDatarld.alt$Replicate <- as.factor(pcaDatarld.alt$Replicate )
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))

pdf(file.path(PreDE.dir, "PCA23_genes.pdf"), height = 10, width = 10)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",size=3,
          color = "group.1", shape = "Replicate", legend = "right",
          ellipse = F ,mean.point = FALSE,repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#Sample Clustering
dvst <- dist(t(assay(vst_genes)))
hvst <- hclust(dvst)

pdf(file.path(PreDE.dir, "Clustering_vst_genes.pdf"), height = 7, width = 7)
myplclust(hvst, labels= rownames(anno), lab.col=as.numeric((anno$group)), main="", cex=0.2)
dev.off()


#Remove batch from transformed Data 
library(limma)
vst_genes <- vst(dds)
BR<-removeBatchEffect(assay(vst_genes), batch=colData(vst_genes)$Replicate)
vst_genes_b <- SummarizedExperiment(assays =BR , colData=colData(vst_genes))
vst_genes_b <- DESeqTransform(vst_genes_b)
saveRDS(vst_genes_b,file =file.path(results.dir, "vst_genes_b_replicate.rds"))


#PCA
pcaDatarld<- plotPCA(vst_genes_b, intgroup = c("group", "Replicate" ), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_afB.pdf"), height = 10, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",
          color = "group.1", shape = "Replicate", legend = "right",
          ellipse = F ,mean.point = FALSE,repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#PC2 and 3
pcaDatarld.alt<- plotPCA.alt(vst_genes_b, intgroup = c("group", "Replicate" ), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))

pdf(file.path(PreDE.dir, "PCA23_afB.pdf"), height = 10, width = 10)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",
          color = "group.1", shape = "Replicate", legend = "right",
          ellipse = F ,mean.point = FALSE,repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()





