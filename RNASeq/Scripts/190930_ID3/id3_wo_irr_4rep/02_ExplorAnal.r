#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)

#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Data Read-in
#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)

#Extracting transformed values
vst <- vst(dds)
anno <- colData(vst)
saveRDS(vst,file = file.path(results.dir,"vst.rds"))

rld <- rlog(dds)
anno <- colData(rld)
saveRDS(rld,file = file.path(results.dir,"rld.rds"))

#function to plot PC1 and 2
pcaDatarld<- plotPCA(rld, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
pcaDatarld$Replicate <- as.factor(pcaDatarld$Replicate )
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12.pdf"), height = 4, width = 7)
ggscatter(pcaDatarld, x="PC1", y="PC2",
            size=5,
            color = "genotype", 
            #shape = "irradiation_treatment", 
            label="Replicate",
            repel=TRUE,
            legend = "right",
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC1: ", percentvarrld[1], "% variance")), ylab=(paste0("PC2: ", percentvarrld[2], "% variance")) )
dev.off()

#function
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
#function to plot PC2 and 3
pcaDatarld.alt<- plotPCA.alt(rld, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
pcaDatarld.alt$Replicate <- as.factor(pcaDatarld.alt$Replicate )
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))

pdf(file.path(PreDE.dir, "PCA23.pdf"), height = 4, width = 7)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",            
            size=5,
            color = "genotype", 
            shape = "irradiation_treatment", 
            label="Replicate",
            repel=TRUE,
            legend = "right",
          ellipse = F ,mean.point = FALSE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#Sample Clustering
drld <- dist(t(assay(rld)))
hrld<- hclust(drld)

pdf(file.path(PreDE.dir, "Clustering_rld.pdf"), height = 7, width = 7)
myplclust(hrld, labels= rownames(anno), lab.col=as.numeric((anno$group)), main="", cex=0.2)
dev.off()







#Remove batch from transformed Data 
library(limma)
#BR<-removeBatchEffect(assay(rld), batch=colData(rlöd)$Replicate, batch2= colData(rld)$Plattform)
BR<-removeBatchEffect(assay(rld), batch=colData(rld)$Replicate)
rld_b <- SummarizedExperiment(assays =BR , colData=colData(rld))
rld_b <- DESeqTransform(rld_b)
saveRDS(rld_b,file =file.path(results.dir, "rld_b_replicate.rds"))


#PCA
pcaDatarld<- plotPCA(rld_b, intgroup =colnames(colData(rld_b)), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_afB.pdf"), height = 7, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",
            size=5,
            color = "genotype", 
            #shape = "irradiation_treatment", 
            label="Replicate",
            repel=TRUE,
            legend = "right",
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#PC2 and 3
pcaDatarld.alt<- plotPCA.alt(rld_b, intgroup = colnames(colData(rld_b)), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))

pdf(file.path(PreDE.dir, "PCA23_afB.pdf"), height = 4, width = 5)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",
            size=5,
            color = "genotype", 
            shape = "irradiation_treatment", 
            label="Replicate",
            repel=TRUE,
            legend = "right",
            ellipse = F ,mean.point = FALSE,
            star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()


#Sample Clustering
topVarGenesRld_b <- head(order(rowVars(assay(rld_b)), decreasing=TRUE),Inf)

drld_b <- dist(t(assay(rld_b[topVarGenesRld_b,])))
hcrld_b <- hclust(drld_b)

pdf(file.path(PreDE.dir, "clustering_RLD_afB.pdf"), height = 7, width = 7)
myplclust(hcrld_b, labels=paste0(colData(rld_b)$group),  main="",lab.col= as.numeric(colData(rld_b)$group),cex=0.2)
dev.off()

dissimilarity <- 1 - cor(assay(rld_b), use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)

pdf(file.path(PreDE.dir, "Correlation_RLD_afB.pdf"), height = 7, width = 5)
myplclust(hclust(distance),hang=0.2,cex=0.3,labels=paste0(colData(rld_b)$group),  main="Distance",lab.col= as.numeric(colData(rld_b)$group))
dev.off()

#Scaled clustering
scale_rows <- function(x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}
scale_rld_b <-scale_rows(assay(rld_b))

d <- dist(t(scale_rld_b), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 

pdf(file.path(PreDE.dir, "Scaled_cllustering_RLD_afB.pdf"), height = 4, width = 5)
myplclust(fit,hang=0.01,cex=0.5,labels=paste0(colData(rld_b)$SampleID),  main="Distance",lab.col= as.numeric(colData(rld_b)$group))# display dendogram
dev.off()




