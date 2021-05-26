---
title: 'Ali''s RNAseq: All samples together:04_DiffExpr'
output:
  html_document:
    theme: lumen
    highlight: pygments
    toc_float: true
    toc: true
    fig_align: left
---
Perform differential epression analysis of different contrasts
#libraries

library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(knitr)
library(DT)
library(clusterProfiler)
library(Glimma)
library(limma)
library(edgeR)

#folder

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190531_odcf_NoIll/"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data

dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
anno <- colData(dds_genes)

#Take a look at design and annotation of samples

design(dds_genes)

#extract with results
#Set specifications

alpha <- 0.05 #set FDR cutoff
lfc <- 0 #set logfold2 cutoff


#set up all possible contrasts
contrasts <- as.data.frame(combn(as.character(unique(colData(dds_genes)$group)), 2))

results <- list()

for (i in 1:length(contrasts)) {
  results[i]<- results(dds_genes, contrast=c("group",as.character(contrasts[1,i]),  as.character(contrasts[2,i])), alpha = alpha, lfcThreshold = lfc) #extract results of wanted comparison
  print(paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i])))
  print(summary(results[i]))
}

nam<-vector()
for (i in 1:length(contrasts)) {
nam[i]<- paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i]))
}
names(results)<- nam


#annotate samples
#Make a dataframe out of it
DEG_results_list <- lapply(results, function(x){
  x<- as.data.frame(x)
  x
})
#create ensembl annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$ensembl <- sapply(strsplit(rownames(x) ,".", fixed=TRUE),`[`, 1)
  x
})
#add symbol annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$symbol<- mapIds(org.Hs.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )
  x
})
#add entrezgene annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$entrezgene<- mapIds(org.Hs.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "ENTREZID", multiVals = "first" )
  x
})
#add gene name annotation
DEG_results_list <- lapply(DEG_results_list, function(x){
  x$genename<- mapIds(org.Hs.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "GENENAME", multiVals = "first" )
  x
})
#order by padjusted value
DEG_results_list <- lapply(DEG_results_list, function(x){
  x <- x[order(x$padj),]
  x
})
names(DEG_results_list)<-names(results)

#save lists
saveRDS(DEG_results_list, file.path(PostDE.dir, "DEG_results_group_list.rds"))

for (i in names(DEG_results_list)){
  dir.create(file.path(PostDE.dir, i))
}
for (i in names(DEG_results_list)){
  write.table(DEG_results_list[[i]], file=file.path(PostDE.dir,i,  paste0(i, "_DEG.txt")),row.names =T, col.names = T  )
}

#Plot results in Volcano and MA plot {.tabset}

DEG_results_list_plot <- DEG_results_list
DEG_results_list_plot <- lapply(DEG_results_list_plot, function(x){
  x <- as.data.frame(x)
  x$col <- "black"
  x$col <- ifelse(test = x$padj < 0.05 & x$log2FoldChange>1, yes = "red", 
                  no = ifelse(test = x$padj < 0.05  & x$log2FoldChange<(-1), yes = "blue", "black"))
  x$pvalue <-(-(log10(x$pvalue)))
  x$padj <- (-(log10(x$padj)))
  x
})

#Volcano Plot
Volcanos<- list()
for (i in names(DEG_results_list_plot)){
  Volcanos[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="log2FoldChange", y = "padj",
          xlab="Gene expression (log2 fold change)",#xscale="log10" ,
          ylab="- log10(P.adjusted)", #xlim=c(-0.5,0.5),
          color = "col" ,shape = 16, size = 1.5, palette=c("grey","red","blue"), # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance", main=i# Customize reg. line
          #label="SYMBOL", repel = T, font.label = c( "black"), 
) + rremove("legend") + geom_vline(xintercept=c(-1,1), linetype = 2) +geom_hline(yintercept=1.3, linetype = 2)
  pdf(file.path(PostDE.dir,i,"volcano.pdf"))
  print(Volcanos[[i]])
  dev.off()
}
#print(Volcanos)

#MA plot
MAs <- list()
for (i in names(DEG_results_list_plot)){
  MAs[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="baseMean", y =  "log2FoldChange",xlab="Mean of normalized counts",xscale="log10" ,
          ylab="Gene expression (log2 fold change)",
          color = "col",palette=c("grey","red","blue") ,shape = 16, size = 1.5, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray") #legend.title="Significance",# Customize reg. line
          #label="SYMBOL", repel = TRUE#, legend = "right",
) +geom_hline(yintercept=c(-1, 1), linetype = 2) + rremove("legend")
  pdf(file.path(PostDE.dir,i,"MA.pdf"))
  print(MAs[[i]])
  dev.off()
}
#print(MAs)

#heatmaps

#define padj cutoff for plotting
cutoff <- 0.1
genes2plot <- lapply(DEG_results_list, function(x){
  x <- x[which(x$padj < cutoff),]
  x <- x$ensembl
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]

#heatmap
vst_genes<- readRDS(file = file.path(results.dir,"vst_genes_batchRem_rep.rds"))
mat_genes<- assay(vst_genes)
rownames(mat_genes)<- sapply(strsplit(rownames(mat_genes) ,".", fixed=TRUE),`[`, 1)
#mat_genes$symbol<- mapIds(org.Hs.eg.db, keys=rownames(mat_genes), keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )

annovst <- colData(dds_genes)
annovst <- annovst[, c("genotype", "tamoxifen_treatment","hdac_treatment", "replicate" )]
annovst$replicate <- as.factor(annovst$replicate)
heat<-list()
for (i in names(genes2plot)){
  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[[i]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
  annorow <- rownames(plot)
  heat[[i]]<- pheatmap(plot[, rownames(annovst)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),main=i) 
  pdf(file.path(PostDE.dir,i,"Heat.pdf"))
  print(heat[[i]])
  dev.off()
}

#GSEA
library(clusterProfiler)
GSEA <- list()
GSKEGG<- list()
for (i in names(DEG_results_list)){
  geneList <- DEG_results_list[[i]][,"log2FoldChange"]
  names(geneList) <- DEG_results_list[[i]]$ensembl
  geneList <- sort(geneList, decreasing=TRUE)
  GSEA[[i]] <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType = "ENSEMBL",
              ont          = "ALL",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 1,
              verbose      = FALSE)
  GSEA[[i]] <- setReadable(GSEA[[i]], OrgDb = org.Hs.eg.db)

  pdf(file.path(PostDE.dir,i,"GSEA_Dotplot.pdf"))
  dotplot(GSEA[[i]],showCategory = 20)
  dev.off()
  print(i)
#GSKEGG[[i]] <- gseKEGG(geneList     = geneList,
#               organism     = "hs",
#               nPerm        = 1000,
#               minGSSize    = 120,
#               pvalueCutoff = 0.05,
#               verbose      = FALSE)
#GSKEGG[[i]] <- setReadable(GSKEGG[[i]], OrgDb = org.Mm.eg.db)

#  pdf(file.path(PostDE.dir,i,"GSEKEGG_Dotplot.pdf"))
#  dotplot(GSKEGG[[i]],showCategory = 20)
#  dev.off()
}
plot <- list()
for (i in names(GSEA)){
  plot[[i]] <- dotplot(GSEA[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"GSEA.pdf"))
  print(plot[[i]])
  dev.off()
}

saveRDS(GSEA, file.path(PostDE.dir, "GSEA.rds"))




#Glimma plots

#base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
#data.dir <- file.path(base.dir, "data")
#ase_results.dir <- file.path(base.dir, "results")
#esults.dir<- file.path(base_results.dir , "tables")
#PreDE.dir <- file.path(base_results.dir,"PreDE")
#PostDE.dir <- file.path(base_results.dir,"PostDE")

#dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
#vst_genes<- readRDS(file = file.path(results.dir,"vst_genes_batchRem_rep.rds"))

#DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

# convert result table
res.df <- lapply(DEG_results_list, function(x) {
x <- as.data.frame(x)
x$log10MeanNormCount <- log10(x$baseMean)
x
})
# remove very low counts / it would be sufficent to remove "zeros"
idx<-rowSums(DESeq2::counts(dds_genes))>0
res.df.fil<- list()
for (i in names(DEG_results_list)){
res.df.fil[[i]] <- res.df[[i]][idx,]
res.df.fil[[i]]$padj[is.na(res.df.fil[[i]]$padj)]<- 1
res.df.fil[[i]]$pvalue[is.na(res.df.fil[[i]]$pvalue)]<- 1
}

# to color siginif genes
res.sig <- lapply(res.df.fil, function(x) {
x <- subset(x, padj<0.1 & abs(log2FoldChange)>0.5) 
x
})
wx <- list()
stat <- list()
for (i in names(res.sig)){
wx[[i]] <- which(rownames(res.df.fil[[i]]) %in% rownames(res.sig[[i]]))
stat[[i]]<-rep(0,nrow(res.df.fil[[i]]))
stat[[i]][wx[[i]]]<-rep(1,length(wx[[i]]))
}

#plot
for (i in names(DEG_results_list)){
dir.create(path=file.path(PostDE.dir,i,"Glimma"))
glMDPlot(res.df.fil[[i]], 
         xval = "log10MeanNormCount",
         yval="log2FoldChange",
         counts=DESeq2::counts(dds_genes)[idx,],
         anno=data.frame(GeneID=rownames(dds_genes)[idx]),
         groups=colData(dds_genes)$group, samples=colnames(dds_genes),
         display.columns = c("ensembl","symbol", "padj","pvalue"),
         status = stat[[i]], path=file.path(PostDE.dir,i,"Glimma"), las=2)
}

for (i in names(DEG_results_list)) {
  dir.create(file.path(PostDE.dir,i,"Glimma_new"))
  glMDPlot(DEG_results_list[[i]],
           status=sign[[i]], 
           counts=counts(dds_plot[[i]], normalized=TRUE), 
           groups=as.character(colData(dds_plot[[i]])$group), 
           transform=FALSE, 
           samples=as.character(colnames(counts(dds_plot[[i]]))), 
           anno=rownames(counts(dds_plot[[i]])),
           path=file.path(PostDE.dir,i,"Glimma_new"))
}


#which groups to plot 
anno <- colData(dds_genes)
group1 <- sapply(strsplit(names(DEG_results_list) ,"_vs_", fixed=TRUE),`[`, 1)
group2 <- sapply(strsplit(names(DEG_results_list) ,"_vs_", fixed=TRUE),`[`, 2)

counts_plot <- list()
sample_list <-list()
group_list <- list()
anno_list <- list()

for (i in 1:length(DEG_results_list)){
counts_plot[[i]] <- DESeq2::counts(dds_genes, normalized=TRUE)[idx,rownames(anno[which(anno$group %in% c(group1[i], group2[i])), ])]
sample_list[[i]] <- rownames(anno[which(anno$group %in% c(group1[i], group2[i])), ])
group_list[[i]] <- anno[which(anno$group %in% c(group1[i], group2[i])), ]$group
}

names(counts_plot)<- names(DEG_results_list)
names(sample_list)<- names(DEG_results_list)
names(group_list)<- names(DEG_results_list)


for (i in names(DEG_results_list)){
dir.create(path=file.path(PostDE.dir,i,"Glimma_sub"))
glMDPlot(res.df.fil[[i]], 
         xval = "log10MeanNormCount",
         yval="log2FoldChange",
         counts=counts_plot[[i]],
         anno=data.frame(GeneID=rownames(dds_genes)[idx]),
         groups=group_list[[i]], samples=sample_list[[i]],
         display.columns = c("ensembl","symbol", "padj","pvalue"),
         status = stat[[i]], path=file.path(PostDE.dir,i,"Glimma_sub"), las=2)
}

dds_plot <- list()
dt.deseq2 <- list()
anno <- colData(dds_genes)
group1 <- sapply(strsplit(names(DEG_results_list) ,"_vs_", fixed=TRUE),`[`, 1)
group2 <- sapply(strsplit(names(DEG_results_list) ,"_vs_", fixed=TRUE),`[`, 2)

for (i in names(DEG_results_list)){
dds_plot[[i]] <- dds_plot[,rownames(anno[which(anno$group %in% c(group1[i], group2[i])), ])]
dt.deseq2[[i]] <- as.numeric(DEG_results_list[[i]]$padj<0.05)
}

for (i in names(DEG_results_list)){
dir.create(path=file.path(PostDE.dir,i,"Glimma_new"))
glMDPlot(DEG_results_list[[i]], status=dt.deseq2[[i]] ,path=file.path(PostDE.dir,i,"Glimma_new"), counts=DESeq2::counts(dds_genes, normalized=TRUE),
groups=group_list[[i]], transform=FALSE,samples=sample_list[[i]], anno=rownames(DESeq2::counts(dds_genes, normalized=TRUE)))
}

library(DESeq2)
rnaseq.deseq2 <- DESeqDataSetFromMatrix(
rnaseq$counts, colData=rnaseq$samples, design=~group)
mcols(rnaseq.deseq2) <- DataFrame(mcols(rnaseq.deseq2), rnaseq$genes)
rnaseq.deseq2 <- DESeq(rnaseq.deseq2)
fit.deseq2 <- results(rnaseq.deseq2, contrast=c("group", "Smchd1-null", "WT"))
 <- as.numeric(fit.deseq2$padj<0.05)
glMDPlot(fit.deseq2, status=dt.deseq2, counts=rnaseq, groups=groups, transform=FALSE,
samples=colnames(rnaseq), anno=rnaseq$genes)


