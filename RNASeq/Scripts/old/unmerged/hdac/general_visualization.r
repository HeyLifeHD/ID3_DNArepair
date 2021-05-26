#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)

#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
counts_hdac<- readRDS(file.path(data.dir, "counts_hdac.rds"))

#Data Read-in
#DESeq2 Analysis
dds_hdac <- readRDS(file = file.path(results.dir,"dds_hdac.rds"))
anno_hdac <- colData(dds_hdac)

#Extracting transformed values
#rld_hdac <- rlog(dds_hdac, blind=FALSE)
#saveRDS(rld_hdac,file = file.path(results.dir,"rld_hdac.rds"))

vst_hdac <- vst(dds_hdac)
saveRDS(vst_hdac,file = file.path(results.dir,"vst_hdac.rds"))

#function to plot PC2 and 3
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
pcaDatarld.alt<- plotPCA.alt(vst_hdac, intgroup = c("sample_name", "replicate", "genotype", "group", "hdac_treatment", "dox_treatment"), returnData =TRUE, ntop=10000)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))


pdf(file.path(PreDE.dir, "hdac_PCA23.pdf"), height = 7, width = 7)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",
          color = "genotype", shape = "dox_treatment",
          ellipse = F ,mean.point = FALSE,label="group.1",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()


#Sample Clustering
dvst <- dist(t(assay(vst)))
hvst <- hclust(dvst)

pdf(file.path(PreDE.dir, "Clustering_vst.pdf"), height = 12, width = 7)
myplclust(hvst, labels=anno$group, lab.col= as.numeric(colData(vst)$genotype), main="", cex=0.5)
dev.off()

#Gene clustering: Heatmaps
#500 most variable Genes
topVarGenesvst<- head(order(rowVars(assay(vst)), decreasing=TRUE),2000)

matvst <- assay(vst)[topVarGenesvst,]
annovst <- as.data.frame(colData(vst)[,c( "replicate", "genotype", "dox_treatment", "hdac_treatment", "irr_treatment")])
annovst$replicate <- as.factor(annovst$replicate)
pdf(file.path(PreDE.dir,"Heatmap2000vst.pdf"),height= 25)
pheatmap(matvst, annotation_col=annovst ,scale="row")
dev.off()


