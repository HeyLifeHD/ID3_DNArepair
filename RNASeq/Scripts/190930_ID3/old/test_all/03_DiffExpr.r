#Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(knitr)
library(clusterProfiler)
library(Glimma)
library(limma)
library(edgeR)

#folder
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/test_all"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)

#Take a look at design and annotation of samples
design(dds)

#extract with results
#Set specifications

alpha <- 0.05 #set FDR cutoff
lfc <- 1##set logfold2 cutoff

#Filter genes which are only expressed in 1 sample
length(dds)
idx_dds<- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 1
table(idx_dds)
#dds <- dds[idx,]
#dim(dds)

#Running the differential expression 
#set up all possible contrasts
contrasts <- as.data.frame(combn(as.character(unique(colData(dds)$group)), 2))

results <- list()

for (i in 1:length(contrasts)) {
  results[[i]]<- results(dds, contrast=c("group",as.character(contrasts[1,i]),  as.character(contrasts[2,i])), alpha = alpha, lfcThreshold = lfc) #extract results of wanted comparison
  print(paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i])))
  print(summary(results[[i]]))
}

nam<-vector()
for (i in 1:length(contrasts)) {
nam[i]<- paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i]))
}
names(results)<- nam

#resultshrink <- lfcShrink(dds,contrast=c("group",as.character(contrasts[1,5]),  as.character(contrasts[2,5])), res=results$wt_LPS_vs_tg_LPS )


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
  x$genename<- mapIds(org.Mm.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "GENENAME", multiVals = "first" )
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
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
for (i in names(DEG_results_list)){
  dir.create(file.path(PostDE.dir, i))
}
for (i in names(DEG_results_list)){
  #DEG_results_list[[i]][which(is.na(DEG_results_list[[i]]$padj)),]$padj<-1
  DEG_results_list[[i]]$padj <- round(DEG_results_list[[i]]$padj, 4)
    DEG_results_list[[i]]$log2FoldChange <- round(DEG_results_list[[i]]$log2FoldChange, 4)

    DEG_results_list[[i]] <- DEG_results_list[[i]][!is.na(DEG_results_list[[i]]$padj),]
  write.table(DEG_results_list[[i]], file=file.path(PostDE.dir,i,  paste0(i, "_DEG.txt")),row.names =F, col.names = T , sep="\t", quote=FALSE)
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
rld_b<- readRDS(file = file.path(results.dir,"rld_b.rds"))
mat_genes<- assay(rld_b)
#mat_genes$symbol<- mapIds(org.Hs.eg.db, keys=rownames(mat_genes), keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )

annovst <- colData(dds)
annovst <- annovst[, c("group", "Replicate" )]
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
#GSKEGG<- list()
for (i in names(DEG_results_list)){
  geneList <- DEG_results_list[[i]][,"log2FoldChange"]
  names(geneList) <- DEG_results_list[[i]]$ensembl
  geneList <- sort(geneList, decreasing=TRUE)
  GSEA[[i]] <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              keyType = "ENSEMBL",
              ont          = "ALL",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 1,
              verbose      = FALSE)
  GSEA[[i]] <- setReadable(GSEA[[i]], OrgDb = org.Mm.eg.db)
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



