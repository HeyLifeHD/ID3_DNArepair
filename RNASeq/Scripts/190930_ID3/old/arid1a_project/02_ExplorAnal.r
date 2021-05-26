#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(Rtsne)
#folder

base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/arid1a_project"
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
vst_genes <- vst(dds)
anno <- colData(vst_genes)
anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment)
colData(vst_genes)<- anno 
saveRDS(vst_genes,file = file.path(results.dir,"vst_genes.rds"))

#function to plot PC1 and 2
pcaDatarld<- plotPCA(vst_genes, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
pcaDatarld$Replicate <- as.factor(pcaDatarld$Replicate )
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_genes.pdf"), height = 4, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",size=3,
          color = "group.1", shape = "Replicate", legend = "right",
          ellipse = F ,mean.point = FALSE,repel=TRUE,
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
pcaDatarld.alt$Replicate <- as.factor(pcaDatarld.alt$Replicate )
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))

pdf(file.path(PreDE.dir, "PCA23_genes.pdf"), height = 7, width = 10)
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
BR<-removeBatchEffect(assay(vst_genes), batch=colData(vst_genes)$Replicate, batch2= colData(vst_genes)$Plattform)
vst_genes_b <- SummarizedExperiment(assays =BR , colData=colData(vst_genes))
vst_genes_b <- DESeqTransform(vst_genes_b)
saveRDS(vst_genes_b,file =file.path(results.dir, "vst_genes_b_replicate.rds"))


#PCA
pcaDatarld<- plotPCA(vst_genes_b, intgroup = c("group", "Replicate" ), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA_afB.pdf"), height = 7, width = 10)
ggscatter(pcaDatarld, x="PC1", y="PC2",
          color = "group.1", shape = "Replicate", legend = "right",
          ellipse = F ,mean.point = FALSE,repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#PC2 and 3
pcaDatarld.alt<- plotPCA.alt(vst_genes_b, intgroup = c("group", "Replicate" ), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))

pdf(file.path(PreDE.dir, "PCA2_3_afB.pdf"), height = 4, width = 5)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",
          color = "group.1", shape = "Replicate", legend = "right",
          ellipse = F ,mean.point = FALSE,repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()



#Sample Clustering
topVarGenesRld_b <- head(order(rowVars(assay(vst_genes_b)), decreasing=TRUE),Inf)

drld_b <- dist(t(assay(vst_genes_b[topVarGenesRld_b,])))
hcrld_b <- hclust(drld_b)

pdf(file.path(PreDE.dir, "clustering_RLD_afB2.pdf"), height = 3, width = 3)
myplclust(hcrld_b, labels=paste0(colData(vst_genes_b)$group),  main="",lab.col= as.numeric(colData(vst_genes_b)$group))
dev.off()



dissimilarity <- 1 - cor(assay(rld_b), use="pairwise.complete.obs")
distance <- as.dist(dissimilarity)

pdf(file.path(PreDE.dir, "Correlation_RLD_afB2_replicate.pdf"), height = 4, width = 5)
myplclust(hclust(distance),hang=0.2,cex=0.5,labels=paste0(colData(rld_b)$group),  main="Distance",lab.col= as.numeric(colData(rld_b)$group))
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

pdf(file.path(PreDE.dir, "Scaled_cllustering_RLD_afB_replicate.pdf"), height = 4, width = 5)
myplclust(fit,hang=0.01,cex=0.5,labels=paste0(colData(rld_b)$SampleID),  main="Distance",lab.col= as.numeric(colData(rld_b)$group))# display dendogram
dev.off()
pdf(file.path(PreDE.dir, "Scaled_cllustering_RLD_afB2_replicate.pdf"), height = 4, width = 5)
myplclust(fit,labels=paste0(colData(rld_b)$SampleID),  main="Distance",lab.col= as.numeric(colData(rld_b)$group))# display dendogram
dev.off()
#Gene clustering: Heatmaps
#500 most variable Genes
topVarGenesRld_b <- head(order(rowVars(assay(rld_b)), decreasing=TRUE),500)
#####BETTER with MAD
topVarGenesRld_b <- head(order(rowMads(assay(rld_b)), decreasing=TRUE),10000)

matrld_b <- assay(rld_b)[topVarGenesRld_b,]

#annorld_b <- as.data.frame(colData(rld_b)[,c("Genotype", "batch", "EosinophilFactor",  "Eosinophil")])
annorld_b <- as.data.frame(colData(rld_b)[,c("Genotype", "replicate", "Treatment")])
pdf(file.path(PreDE.dir, "Heatmap500mV_RLD_afB_replicate.pdf"),height= 25)
pheatmap(matrld_b,show_colnames= FALSE, annotation_col=annorld_b , scale="row",fontsize_row=5)
dev.off()













subresult <- subset(result ,abs(log2FoldChange)>0.5 &padj <0.1)
subresult <- as.data.frame(subresult)
#topVarGenesRld_b <- head(order(rowVars(assay(rld_b)[rownames(rld_b)%in%rownames(subresult)]), decreasing=TRUE),1000)


annorld_b <- as.data.frame(colData(rld_b)[,c("Genotype", "batch", "Eosinophil")])
matrld_b <- assay(rld_b)[rownames(assay(rld_b))%in%rownames(subresult) ,]
#annorld_b <- as.data.frame(colData(rld_b)[,c("Genotype", "batch")])
png(file.path(PreDE.dir,  "Heatmap_DE.png"))
pheatmap(matrld_b, annotation_col=annorld_b , labels_row = ensembl2symbol[row.names(matrld_b)],
           fontsize_row=5,  scale="row",  show_colnames= FALSE, show_rownames = TRUE,  color=brewer.pal(11,"RdBu")
) 
dev.off()

##DE_Genes
###Import DE genes

Sig <- subset(result, padj <0.1)
Sig <- result[sapply(result$padj<0.1, isTRUE),]
ind <- order(Sig$padj, decreasing=FALSE)
Sig <- Sig [ind,]
matrld_b <- assay(rld_b)
matrld_b<-matrld_b[rownames(matrld_b) %in% rownames(Sig),]
annorld_b <- as.data.frame(colData(rld_b)[,c("Genotype", "batch")])

pdf(paste(dir=PreDE.dir, "HeatmapDE_nocluster_RLD_afB_", ct,".pdf", sep=""),height= 15)
pheatmap(matrld_b, scale = "row",annotation_col=annorld_b , labels_row = ensembl2symbol[row.names(matrld_b)],fontsize_row=5, cluster_rows = FALSE,show_colnames= FALSE)
dev.off()

png(file.path(PreDE.dir,  "2HeatmapDE_RLD_afB.png"))
pheatmap(matrld_b,scale="row", annotation_col=annorld_b , labels_row = ensembl2symbol[row.names(matrld_b)],fontsize_row=5, show_colnames= FALSE, show_rownames = FALSE)
dev.off()


#Scattertplot
CT=colData(rld_b)
WT <-  as.vector(CT[CT$Genotype == "wt",]$Sample_name)
TG <- as.vector(CT[CT$Genotype == "tg",]$Sample_name)

bin<-hexbin(rowMeans(assay(rld_b)[ ,WT]), rowMeans(assay(rld_b)[ ,TG]), xbins=75)

pdf(file.path(PreDE.dir,  "SmoothScatter.pdf"))
my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
par(mgp=c(3, 1, 0), cex.axis=1.5)
plot(bin, main="" , colramp=my_colors ,
     colorcut = 20,
     legend=1.2 , xlab="Rlog transformed Counts, WT", ylab="Rlog transformed Counts, TG", clip="off")

dev.off()


