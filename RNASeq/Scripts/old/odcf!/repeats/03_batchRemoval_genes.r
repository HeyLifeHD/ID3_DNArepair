#Remove batch effects
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(limma)
#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/repeats"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

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

#load vst files for each of the data sets
vst_repeats_or <- readRDS(file = file.path(results.dir,"vst_repeats.rds"))
vst_genes_or<- readRDS(file = file.path(results.dir,"vst_genes.rds"))

anno <- colData(vst_repeats)

#remove batch effect of replicates
mat_genes <- assay(vst_genes_or)
mat_genes <- limma::removeBatchEffect(mat_genes, vst_genes$replicate)
assay(vst_genes) <- mat_genes
#plot batch removed data
#pc 1 and 2
pcaDatarld<- plotPCA(vst_genes, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_batchRem_rep_genes.pdf"), height = 4, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",size=3,
          color = "group_uncomplete2", shape = "genotype", legend = "right",
          ellipse = F ,mean.point = FALSE,label="replicate",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC1: ", percentvarrld[1], "% variance")), ylab=(paste0("PC2: ", percentvarrld[2], "% variance")) )
dev.off()

#pc2 and 3
pcaDatarld.alt<- plotPCA.alt(vst_genes, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))


pdf(file.path(PreDE.dir, "PCA23_batchRem_rep_genes.pdf"), height = 4, width = 10)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",size=3,
          color = "group_uncomplete2", shape = "genotype", legend = "right",
          ellipse = F ,mean.point = FALSE,label="replicate",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()



#Sample Clustering
dvst <- dist(t(assay(vst_genes)))
hvst <- hclust(dvst)

pdf(file.path(PreDE.dir, "Clustering_batchRem_rep_genes.pdf"), height = 12, width = 7)
myplclust(hvst, labels=anno$group_complete, lab.col= as.numeric(as.factor(anno$group_uncomplete2)), main="", cex=0.5)
dev.off()

#remove batch effect of sequencer
mat_genes <- assay(vst_genes_or)
mat_genes <- limma::removeBatchEffect(mat_genes, vst_genes$sequencer)
assay(vst_genes) <- mat_genes

#plot PCA of batch removed data
#pc 1 and 2
pcaDatarld<- plotPCA(vst_genes, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_batchRem_seq_genes.pdf"), height = 4, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",size=3,
          color = "group_uncomplete2", shape = "genotype", legend = "right",
          ellipse = F ,mean.point = FALSE,label="replicate",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC1: ", percentvarrld[1], "% variance")), ylab=(paste0("PC2: ", percentvarrld[2], "% variance")) )
dev.off()

#pc 2 and 3
pcaDatarld.alt<- plotPCA.alt(vst_genes, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))


pdf(file.path(PreDE.dir, "PCA23_batchRem_seq_genes.pdf"), height = 4, width = 10)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",size=3,
          color = "group_uncomplete2", shape = "genotype", legend = "right",
          ellipse = F ,mean.point = FALSE,label="replicate",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#Sample Clustering
dvst <- dist(t(assay(vst_genes)))
hvst <- hclust(dvst)

pdf(file.path(PreDE.dir, "Clustering_batchRem_seq_genes.pdf"), height = 12, width = 7)
myplclust(hvst, labels=anno$group_complete, lab.col= as.numeric(as.factor(anno$group_uncomplete2)), main="", cex=0.5)
dev.off()

#decide on which batch removal should be chosen --> replicate
mat_genes <- assay(vst_genes_or)
mat_genes <- limma::removeBatchEffect(mat_genes, vst_genes$replicate)
assay(vst_genes) <- mat_genes
saveRDS(vst_genes, file = file.path(results.dir,"vst_genes_batchRem_rep.rds"))



#residuals function by Reka; try it tomorrow
residuals <- function(meth, pheno, adjustment_columns=c()){
  
  residual_values <- apply(meth, 1, function(x) {
    pheno$x <- x
    form <- as.formula(paste0("x ~ ", paste(adjustment_columns, collapse="+")))
    mod<-lm(form, data = pheno)
    pval<-pf(summary(mod)$fstatistic[1], summary(mod)$fstatistic[2], summary(mod)$fstatistic[3], lower.tail = FALSE)
    if(pval<0.01){ 
      stats::residuals(mod) 
    }else{
      scale(x, scale = FALSE)
    }
  })
  return(residual_values)
} 
pheno <- colData(vst_genes_or)
pheno <- pheno[, c("replicate","sequencer","genotype","hdac_treatment","tamoxifen_treatment", "irradiation_treatment"), drop=FALSE]
pheno$replicate <- as.factor(pheno$replicate)
residuals<- residuals(as.data.frame(assay(vst_genes_or)), pheno, "replicate")