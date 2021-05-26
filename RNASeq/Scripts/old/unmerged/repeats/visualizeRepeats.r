#subset hddac sample set
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(rtracklayer)

#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
repeats.dir <- file.path(base_results.dir, "repeats")
#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds_repeats.rds"))
anno <- colData(dds)

#Extracting transformed values
rld <- rlog(dds, blind=FALSE)
saveRDS(rld,file = file.path(results.dir,"rld_repeats.rds"))
#rld<-readRDS(file.path(results.dir,"rld.rds"))


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
pcaDatarld.alt<- plotPCA.alt(rld, intgroup = c("sample_name", "replicate", "genotype", "dox_treatment", "hdac_treatment", "irr_treatment", "group"), returnData =TRUE, ntop=10000)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))


pdf(file.path(repeats.dir, "PCA23_repeats.pdf"), height = 7, width = 7)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",
          color = "genotype", shape = "hdac_treatment",
          ellipse = F ,mean.point = FALSE,label="group.1",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()


#Sample Clustering
drld <- dist(t(assay(rld)))
hrld <- hclust(drld)

pdf(file.path(repeats.dir, "Clustering_rld_repeats.pdf"), height = 12, width = 7)
myplclust(hrld, labels=anno$group, lab.col= as.numeric(colData(rld)$genotype), main="", cex=0.5)
dev.off()

#repeat clustering: Heatmaps
#100 most variable repeats
topVarGenesrld<- head(order(rowVars(assay(rld)), decreasing=TRUE),100)

matrld <- assay(rld)[topVarGenesrld,]
annorld <- as.data.frame(colData(rld)[,c( "replicate", "genotype", "dox_treatment", "hdac_treatment", "irr_treatment")])
annorld$replicate <- as.factor(annorld$replicate)
pdf(file.path(repeats.dir,"Heatmap_100mv_repeatst.pdf"),height= 25)
pheatmap(matrld, annotation_col=annorld ,scale="row",
labels_row=counts_sum_repName[topVarGenesrld,"repName"],)
dev.off()

pdf(file.path(repeats.dir,"Heatmap_noScale_100mv_repeatst.pdf"),height= 25)
pheatmap(matrld, annotation_col=annorld ,scale="none",
labels_row=counts_sum_repName[topVarGenesrld,"repName"],)
dev.off()

#repeat clustering:Heatmaps
#LTR
ltr_rld <- assay(rld)[grep("LTR12", counts_sum_repName$repName),]
annorld <- as.data.frame(colData(rld)[,c( "replicate", "genotype", "dox_treatment", "hdac_treatment", "irr_treatment")])
annorld$replicate <- as.factor(annorld$replicate)

pdf(file.path(repeats.dir,"Heatmap_LTR_repeatst.pdf"))
pheatmap(ltr_rld, annotation_col=annorld ,scale="row",
labels_row=counts_sum_repName[grep("LTR12", counts_sum_repName$repName),"repName"],)
dev.off()