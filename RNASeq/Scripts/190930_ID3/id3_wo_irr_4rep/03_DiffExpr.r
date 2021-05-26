#Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(p)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(knitr)
library(clusterProfiler)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)

#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)
rld <- readRDS(file =file.path(results.dir, "rld_b_replicate.rds"))

#Take a look at design and annotation of samples
design(dds)

#extract with results
#Set specifications

alpha <- 0.05 #set FDR cutoff
lfc <- 0##set logfold2 cutoff

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
names(contrasts)<- nam
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
  x$genename<- mapIds(org.Hs.eg.db, keys=x$ensembl, keytype ="ENSEMBL", column = "GENENAME", multiVals = "first" )
  x
})
#order by padjusted value
DEG_results_list <- lapply(DEG_results_list, function(x){
  x <- x[order(x$padj),]
  x
})
names(DEG_results_list) <- names(results)

#save lists
saveRDS(DEG_results_list, file.path(PostDE.dir, "DEG_results_group_list.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))


for (i in names(DEG_results_list)){
    dir.create(file.path(PostDE.dir, i))
  DEG_results_list[[i]][which(is.na(DEG_results_list[[i]]$padj)),]$padj<-1
  DEG_results_list[[i]]$padj <- round(DEG_results_list[[i]]$padj, 4)
    DEG_results_list[[i]]$log2FoldChange <- round(DEG_results_list[[i]]$log2FoldChange, 4)

   # DEG_results_list[[i]] <- DEG_results_list[[i]][!is.na(DEG_results_list[[i]]$padj),]
  write.table(DEG_results_list[[i]], file=file.path(PostDE.dir,i,  paste0( "DEG.txt")),row.names =F, col.names = T , sep="\t", quote=FALSE)
}




#Plot results in Volcano and MA plot {.tabset}
DEG_results_list_plot <- DEG_results_list
DEG_results_list_plot <- lapply(DEG_results_list_plot, function(x){
  x <- as.data.frame(x)
  x$col <- "black"
  x$col <- ifelse(test = x$padj < 0.05 & x$log2FoldChange>0.5, yes = "red", 
                  no = ifelse(test = x$padj < 0.05  & x$log2FoldChange<(-0.5), yes = "blue", "black"))
  x$pvalue <-(-(log10(x$pvalue)))
  x$padj <- (-(log10(x$padj)))
  x
}) 

#Volcano Plot
library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- c("grey", col[3], col[1])

Volcanos<- list() 
for (i in names(DEG_results_list_plot)){
  Volcanos[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="log2FoldChange", y = "padj",
          xlab="Gene expression (log2 fold change)",#xscale="log10" ,
          ylab="- log10(P.adjusted)", #xlim=c(-0.5,0.5),
          color = "col" ,shape = 16, size = 1.5, palette=col, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance",# main=i# Customize reg. line
          #label="SYMBOL", repel = T, font.label = c( "black"), 
) + rremove("legend") + geom_vline(xintercept=c(-0.5,0.5), linetype = 2) +geom_hline(yintercept=1.3, linetype = 2)
  pdf(file.path(PostDE.dir,i,"volcano.pdf"))
  print(Volcanos[[i]])
  dev.off()
  pdf(file.path(PostDE.dir,i,"volcano2.pdf"), height=3.5, width=3.5)
  print(Volcanos[[i]])
  dev.off()
}

#MA plot
MAs <- list()
for (i in names(DEG_results_list_plot)){
  MAs[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="baseMean", y =  "log2FoldChange",xlab="Mean of normalized counts",xscale="log10" ,
          ylab="Gene expression (log2 fold change)",
          color = "col",palette=col,shape = 16, size = 1.5, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray") #legend.title="Significance",# Customize reg. line
          #label="SYMBOL", repel = TRUE#, legend = "right",
) +geom_hline(yintercept=c(-0.5, 0.5), linetype = 2) + rremove("legend")
  pdf(file.path(PostDE.dir,i,"MA.pdf"))
  print(MAs[[i]])
  dev.off()
  pdf(file.path(PostDE.dir,i,"MA2.pdf"), height=3.5, width=3.5)
  print(MAs[[i]])
  dev.off()
}

#s
#define padj cutoff for plotting
cutoff <- 0.05
genes2plot <- lapply(DEG_results_list, function(x){
  x <- x[which(x$padj < cutoff & abs(x$log2FoldChange)>0.5),]
  x <- rownames(x)
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]

#
mat_genes<- assay(rld)
annovst <- as.data.frame(colData(rld))[, c("genotype"), drop=FALSE] 
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(colData(rld)$genotype)), name = 'Dark2')
names(color_genotype)<- as.character(unique(colData(rld)$genotype))
#color_treatment<-randomColor(length(unique(colData(rld)$irradiation_treatment)))
#names(color_treatment)<- as.character(unique(colData(rld)$irradiation_treatment))
color_replicate<-randomColor(length(unique(colData(rld)$Replicate)))
names(color_replicate)<- as.character(unique(colData(rld)$Replicate))
#anno_colors <- list(genotype=color_genotype, irradiation_treatment=color_treatment, Replicate=color_replicate)
anno_colors <- list(genotype=color_genotype)

heat<-list()

for (i in names(genes2plot)){
  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[[i]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
  annorow <- rownames(plot)
  heat[[i]]<- pheatmap(plot[, rownames(annovst)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors
                       ) 
  pdf(file.path(PostDE.dir,i,"Heat_DEG.pdf"))
  print(heat[[i]])
  dev.off()
    heat[[i]]<- pheatmap(plot[, rownames(annovst)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors,
                       clustering_distance_rows="correlation") 
  pdf(file.path(PostDE.dir,i,"Heat_correlation_DEG.pdf"))
  print(heat[[i]])
  dev.off()
}

#subset only comparison groups
anno <- colData(rld)
for (i in names(genes2plot)){
  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[[i]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
  annorow <- rownames(plot)
  heat[[i]]<- p(plot[, c(rownames(anno[anno$group==as.character(contrasts[1,i]),]), rownames(anno[anno$group==as.character(contrasts[2,i]),]))], 
  scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors
                       ) 
  pdf(file.path(PostDE.dir,i,"Heat_DEG_subgroups.pdf"))
  print(heat[[i]])
  dev.off()
    heat[[i]]<- p(plot[, c(rownames(anno[anno$group==as.character(contrasts[1,i]),]), rownames(anno[anno$group==as.character(contrasts[2,i]),]))], 
    scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),show_rownames=F,annotation_colors=anno_colors,
                       clustering_distance_rows="correlation") 
  pdf(file.path(PostDE.dir,i,"Heat_correlation_DEG_subgroups.pdf"))
  print(heat[[i]])
  dev.off()
}
#subset for manuscript
pheno <- colData(rld)
rld_sub <- rld[,rownames(pheno[pheno$genotype %in% c("ID3", "WT", "MDC1"),])]
mat_genes<- assay(rld_sub)
annovst <- as.data.frame(colData(rld_sub))[, c("genotype"), drop=FALSE] 
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(colData(rld_sub)$genotype)), name = 'Dark2')
names(color_genotype)<- as.character(unique(colData(rld_sub)$genotype))
anno_colors <- list(genotype=c(ID3="#D95F02",WT="#7570B3",MDC1="#E7298A" ))

heat<-list()

for (i in names(genes2plot)){
  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[[i]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
  annorow <- rownames(plot)
  heat[[i]]<- pheatmap(plot[, rownames(annovst)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors
                       ) 
  pdf(file.path(PostDE.dir,i,"Heat_DEG_new.pdf"))
  print(heat[[i]])
  dev.off()
    heat[[i]]<- pheatmap(plot[, rownames(annovst)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors,
                       clustering_distance_rows="correlation") 
  pdf(file.path(PostDE.dir,i,"Heat_correlation_DEG_new.pdf"))
  print(heat[[i]])
  dev.off()
}

#common 
DEG <- unique(unlist(genes2plot))
mat_genes<- assay(rld)
plot <- mat_genes[which(rownames(mat_genes) %in% DEG ),]
set.seed(123)
p <- p(plot, scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),
                       annotation_colors=anno_colors, show_rownames=F, clustering_distance_rows="correlation" ,
                       #cutree_rows =5, 
                       file=file.path(PostDE.dir, "Heat_allDEG_correlation.pdf")) 

#clust_genes <- sort(cutree(p$tree_row, k=5))
#for (i in unique(clust_genes)){
#temp <- names(clust_genes[clust_genes==i])
#write.table(temp, file.path(PostDE.dir, paste0("cluster_",i, "_genes.txt")), quote=FALSE)
#print(length(temp))
#}

#clust_genes<- as.data.frame(clust_genes)
#enrich for each cluster
#GO<- list()
#for (i in 1:2){
#temp <- names(clust_genes[clust_genes==i])
#  GO[[i]] <- enrichGO(temp, OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
#                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "All")
#
# GO[[i]] <- setReadable(GO[[i]], OrgDb = org.Hs.eg.db)
#}
#plot <- list()
#plot2 <- list()
#for (i in 1:2){
#  plot[[i]] <- dotplot(GO[[i]],showCategory = 20)
#  pdf(file.path(PostDE.dir,paste0("Cluster",i,"_GO_All_dotplot.pdf")), width=12)
#  print(plot[[i]])
#  dev.off()
#  plot2[[i]] <- barplot(GO[[i]],showCategory = 10)
#  pdf(file.path(PostDE.dir,paste0("Cluster",i,"_GO_All_barplot.pdf")), height=4, width=10)
#  print(plot2[[i]])
#  dev.off()
#}

#pies Distribution
DEG_results_list<- lapply(DEG_results_list, function(x){
  x$direction <- NA
  x$direction <- ifelse(x$log2FoldChange>0, "up", "down")
  x
})
for(i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]][which(abs(DEG_results_list[[i]]$log2FoldChange)>0.5 & DEG_results_list[[i]]$padj <0.05),]
  temp <- DEG_results_list[[i]][which(DEG_results_list[[i]]$padj <0.05),]
  temp <- temp[abs(temp$log2FoldChange)>0.5 ,]

  table <- as.data.frame(table(temp$direction))
  labs <- (paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq))

  pdf(file.path(PostDE.dir,i,paste0("Direction_DEGs_sign_fc.pdf")), height=3.5, width=3.5)
  print(ggpie(table, "Freq", fill="Var1", palette=col[c(2,3)], label=labs, lab.pos = "in", main="DEG Analysis\nadj. p-value < 0.05 &\nabs. log2 fold-change > 0.5",  lab.font = c(5, "bold", "white")) + rremove("legend"))
  dev.off()

    temp <- DEG_results_list[[i]][which(DEG_results_list[[i]]$padj <0.05),]
  table <- as.data.frame(table(temp$direction))
  labs <- (paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq))

  pdf(file.path(PostDE.dir,i,paste0("Direction_DEGs_sign.pdf")), height=3.5, width=3.5)
  print(ggpie(table, "Freq", fill="Var1", palette=col[c(2,3)], label=labs, lab.pos = "in", main="DEG Analysis\nadj. p-value < 0.05",  lab.font = c(5, "bold", "white")) + rremove("legend"))
  dev.off()
}

#GSEA
library(clusterProfiler)
#All
GSEA <- list()
for (i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]]
  temp <- temp[!is.na(temp$padj),]
  geneList <-temp[,"log2FoldChange"]
  names(geneList) <- rownames(temp)
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

}
saveRDS(GSEA, file.path(PostDE.dir, "GSEA.rds"))
plot <- list()
for (i in names(GSEA)){
  plot[[i]] <- clusterProfiler::dotplot(GSEA[[i]],showCategory = 20,orderBy = "p.adjust")
  pdf(file.path(PostDE.dir,i,"GSEA_dotplot.pdf"), width=12)
  print(plot[[i]])
  dev.off()
}
#BP
GSEA <- list()
for (i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]]
  temp <- temp[!is.na(temp$padj),]
  geneList <-temp[,"log2FoldChange"]
  names(geneList) <- rownames(temp)
  geneList <- sort(geneList, decreasing=TRUE)
  GSEA[[i]] <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType = "ENSEMBL",
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 1,
              verbose      = FALSE)
  GSEA[[i]] <- setReadable(GSEA[[i]], OrgDb = org.Hs.eg.db)

}
saveRDS(GSEA, file.path(PostDE.dir, "GSEA_BP.rds"))
plot <- list()
for (i in names(GSEA)){
  plot[[i]] <- clusterProfiler::dotplot(GSEA[[i]],showCategory = 20,orderBy = "p.adjust")
  pdf(file.path(PostDE.dir,i,"GSEA_dotplot_BP.pdf"), width=12)
  print(plot[[i]])
  dev.off()
}
#MF
GSEA <- list()
for (i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]]
  temp <- temp[!is.na(temp$padj),]
  geneList <-temp[,"log2FoldChange"]
  names(geneList) <- rownames(temp)
  geneList <- sort(geneList, decreasing=TRUE)
  GSEA[[i]] <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType = "ENSEMBL",
              ont          = "MF",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 1,
              verbose      = FALSE)
  GSEA[[i]] <- setReadable(GSEA[[i]], OrgDb = org.Hs.eg.db)

}
saveRDS(GSEA, file.path(PostDE.dir, "GSEA_MF.rds"))
plot <- list()
for (i in names(GSEA)){
  plot[[i]] <- clusterProfiler::dotplot(GSEA[[i]],showCategory = 20,orderBy = "p.adjust")
  pdf(file.path(PostDE.dir,i,"GSEA_dotplot_MF.pdf"), width=12)
  print(plot[[i]])
  dev.off()
}



#GO analysis
#All
GO_up<- list()
GO_down<- list()

for (i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]]
  temp <- temp[which(temp$padj < 0.05 & temp$log2FoldChange>1),]
  GO_up[[i]] <- enrichGO(rownames(temp), OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "All")
  GO_up[[i]] <- setReadable(GO_up[[i]], OrgDb = org.Hs.eg.db)

  temp <- DEG_results_list[[i]]
  temp <- temp[which(temp$padj < 0.05 & temp$log2FoldChange<(-1)),]
  GO_down[[i]] <- enrichGO(rownames(temp), OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "All")
  GO_down[[i]] <- setReadable(GO_down[[i]], OrgDb = org.Hs.eg.db)
}
saveRDS(GO_up, file.path(PostDE.dir, "GO_sigUp_All.rds"))
saveRDS(GO_down, file.path(PostDE.dir, "GO_sigDown_All.rds"))

plot <- list()
for (i in names(GO_up)){
  plot[[i]] <- dotplot(GO_up[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"GO_sigUp_All_dotplot.pdf"), width=16)
  print(plot[[i]])
  dev.off()

  plot[[i]] <- dotplot(GO_down[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"GO_sigDown_All_dotplot.pdf"), width=16)
  print(plot[[i]])
  dev.off()
}

#MF
GO_up<- list()
GO_down<- list()
for (i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]]
  temp <- temp[which(temp$padj < 0.05 & temp$log2FoldChange>1),]
  GO_up[[i]] <- enrichGO(rownames(temp), OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "MF")
  GO_up[[i]] <- setReadable(GO_up[[i]], OrgDb = org.Hs.eg.db)

  temp <- DEG_results_list[[i]]
  temp <- temp[which(temp$padj < 0.05 & temp$log2FoldChange<(-1)),]
  GO_down[[i]] <- enrichGO(rownames(temp), OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "MF")
  GO_down[[i]] <- setReadable(GO_down[[i]], OrgDb = org.Hs.eg.db)
}
saveRDS(GO_up, file.path(PostDE.dir, "GO_sigUp_MF.rds"))
saveRDS(GO_down, file.path(PostDE.dir, "GO_sigDown_MF.rds"))

plot <- list()
for (i in names(GO_down)){
  plot[[i]] <- dotplot(GO_up[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"GO_sigUp_MF_dotplot.pdf"), width=16)
  print(plot[[i]])
  dev.off()

  plot[[i]] <- dotplot(GO_down[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"GO_sigDown_MF_dotplot.pdf"), width=16)
  print(plot[[i]])
  dev.off()
}

#BP
GO_up<- list()
GO_down<- list()
for (i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]]
  temp <- temp[which(temp$padj < 0.05 & temp$log2FoldChange>1),]
  GO_up[[i]] <- enrichGO(rownames(temp), OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "BP")
  GO_up[[i]] <- setReadable(GO_up[[i]], OrgDb = org.Hs.eg.db)

  temp <- DEG_results_list[[i]]
  temp <- temp[which(temp$padj < 0.05 & temp$log2FoldChange<(-1)),]
  GO_down[[i]] <- enrichGO(rownames(temp), OrgDb = org.Hs.eg.db,keyType="ENSEMBL",  
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "BP")
  GO_down[[i]] <- setReadable(GO_down[[i]], OrgDb = org.Hs.eg.db)
}
saveRDS(GO_up, file.path(PostDE.dir, "GO_sigUp_BP.rds"))
saveRDS(GO_down, file.path(PostDE.dir, "GO_sigDown_BP.rds"))

plot <- list()
for (i in names(GO_down)){
  plot[[i]] <- dotplot(GO_up[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"GO_sigUp_BP_dotplot.pdf"), width=16)
  print(plot[[i]])
  dev.off()

  plot[[i]] <- dotplot(GO_down[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"GO_sigDown_BP_dotplot.pdf"), width=16)
  print(plot[[i]])
  dev.off()
}

