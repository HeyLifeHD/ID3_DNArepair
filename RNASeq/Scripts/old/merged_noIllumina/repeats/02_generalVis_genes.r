#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(Rtsne)
#folder

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190531_odcf_NoIll/repeats"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Data Read-in
#Data Read-in
#DESeq2 Analysis
dds_genes <- readRDS(file = file.path(results.dir,"dds_genes.rds"))
anno <- colData(dds_genes)

#Extracting transformed values
#rld <- rlog(dds, blind=FALSE)
#saveRDS(rld,file = file.path(results.dir,"rld.rds"))
#rld<-readRDS(file.path(results.dir,"rld.rds"))
vst_genes <- vst(dds_genes)
anno <- colData(vst_genes)
anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment)
colData(vst_genes)<- anno
saveRDS(vst_genes,file = file.path(results.dir,"vst_genes.rds"))

#function to plot PC1 and 2
pcaDatarld<- plotPCA(vst_genes, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_genes.pdf"), height = 4, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",size=3,
          color = "group_uncomplete2", shape = "genotype", legend = "right",
          ellipse = F ,mean.point = FALSE,label="replicate",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC1: ", percentvarrld[1], "% variance")), ylab=(paste0("PC2: ", percentvarrld[2], "% variance")) )
dev.off()


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
pcaDatarld.alt<- plotPCA.alt(vst_genes, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))


pdf(file.path(PreDE.dir, "PCA23_genes.pdf"), height = 7, width = 10)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",size=3,
          color = "group_uncomplete2", shape = "genotype", legend = "right",
          ellipse = F ,mean.point = FALSE,label="replicate",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#tsne
#tsne_model_1 = Rtsne(as.matrix(t(matvst[,1:500])), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
#d_tsne_1 = as.data.frame(tsne_model_1$Y)  


#extract eigenvectors of top 100 genes
ir.pca <- prcomp(assay(vst_genes),
                 center = TRUE,
                 scale. = TRUE) 
ev<- as.data.frame(ir.pca$x)
mv_ev_PC2 <- head(ev[order(abs(ev$PC2), decreasing=TRUE),],100)
as.factor(rownames(mv_ev_PC2))

#Sample Clustering
dvst <- dist(t(assay(vst_genes)))
hvst <- hclust(dvst)

pdf(file.path(PreDE.dir, "Clustering_vst_genes.pdf"), height = 12, width = 7)
myplclust(hvst, labels=anno$group_complete, lab.col= as.numeric(as.factor(anno$group_uncomplete2)), main="", cex=0.5)
dev.off()

#Gene clustering: Heatmaps
#500 most variable Genes
topVarGenesvst<- head(order(rowVars(assay(vst_genes)), decreasing=TRUE),10000)

matvst <- assay(vst_genes)[topVarGenesvst,]
annovst <- as.data.frame(colData(vst_genes)[,colnames(anno)])
annovst$replicate <- as.factor(annovst$replicate)
annovst <- annovst[, c("genotype", "tamoxifen_treatment","hdac_treatment", "replicate")]
pdf(file.path(PreDE.dir,"Heatmap1000vst_genes.pdf"),height= 25)
pheatmap(matvst, annotation_col=annovst ,scale="row", show_colnames=F, show_rownames=F)
dev.off()


