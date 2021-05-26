
#Directories
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp_HKGNorm"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")
dir.create(analysis.dir, recursive=TRUE)

#Libraries
library(DiffBind)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(ChIPseeker)
library(rtracklayer)
library(DESeq2)
library(ggpubr)
library(rafalib)
#load data
counts_raw_gr <- readRDS(file.path(analysis.dir,  "data","allPeaks_raw_gr.rds"))
sample_anno <- readRDS(file.path(analysis.dir, "data","sample_anno.rds" ))
dds <- readRDS(file = file.path(analysis.dir, "DESEQ","dds.rds"))

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
#saveRDS(dds,file.path(results.dir,"dds.rds") )

#Running the differential expression 
#set up all possible contrasts
contrasts <- as.data.frame(combn(as.character(unique(colData(dds)$Condition)), 2))

results <- list()
for (i in 1:length(contrasts)) {
  results[[i]]<- results(dds, contrast=c("Condition",as.character(contrasts[1,i]),  
  as.character(contrasts[2,i])), alpha = alpha, lfcThreshold = lfc, format="GRanges") #extract results of wanted comparison
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

#merge with annotated data
for(i in names(results)){
    mcols(results[[i]])<- cbind(mcols(results[[i]]), mcols(counts_raw_gr))
}

#save lists
saveRDS(results, file.path(analysis.dir, "DESEQ", "results.rds"))

DEG_results_list_exp<- list()
for (i in names(results)){
    DEG_results_list_exp[[i]] <- as.data.frame(results[[i]])
    DEG_results_list_exp[[i]]$padj <- round(DEG_results_list_exp[[i]]$padj, 4)
    DEG_results_list_exp[[i]]$log2FoldChange <- round(DEG_results_list_exp[[i]]$log2FoldChange, 4)
    dir.create(file.path(analysis.dir, "DESEQ", i))
    write.table(DEG_results_list_exp[[i]], file=file.path(analysis.dir, "DESEQ", i, "DEG.txt"),row.names =F, col.names = T , sep="\t", quote=FALSE)
}

#Plot results in Volcano and MA plot {.tabset}
DEG_results_list_plot <- results
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
) + rremove("legend") + geom_vline(xintercept=c(-1,1), linetype = 2) +geom_hline(yintercept=1.3, linetype = 2)
  dir.create(file.path(analysis.dir,i))
  pdf(file.path(analysis.dir,i,"volcano.pdf"))
  print(Volcanos[[i]])
  dev.off()
  pdf(file.path(analysis.dir,i,"volcano2.pdf"), height=3.5, width=3.5)
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
) +geom_hline(yintercept=c(-1, 1), linetype = 2) + rremove("legend")
  pdf(file.path(analysis.dir,i,"MA.pdf"))
  print(MAs[[i]])
  dev.off()
  pdf(file.path(analysis.dir,i,"MA2.pdf"), height=3.5, width=3.5)
  print(MAs[[i]])
  dev.off()
}

#heatmaps
#define padj cutoff for plotting
cutoff <- 0.05
genes2plot <- lapply(DEG_results_list, function(x){
  x <- x[which(x$padj < cutoff),]
  x <- rownames(x)
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]
lapply(genes2plot, function(x)length(x))

#heatmap
mat_genes<- assay(vst_b)
annovst <- as.data.frame(colData(vst_b))[, c("group_new", "CellType", "Origin", "model", "seq")] 
color_group<-randomColor(length(unique(colData(vst_b)$group_new)))
names(color_group)<- as.character(unique(colData(vst_b)$group_new))
color_CellType<-randomColor(length(unique(colData(vst_b)$CellType)))
names(color_CellType)<- as.character(unique(colData(vst_b)$CellType))
color_Origin<-randomColor(length(unique(colData(vst_b)$Origin)))
names(color_Origin)<- as.character(unique(colData(vst_b)$Origin))
color_model<-randomColor(length(unique(colData(vst_b)$model)))
names(color_model)<- as.character(unique(colData(vst_b)$model))
color_seq<-randomColor(length(unique(colData(vst_b)$seq)))
names(color_seq)<- as.character(unique(colData(vst_b)$seq))
anno_colors <- list(group_new=color_group,  CellType=color_CellType, Origin=color_Origin, model=color_model, seq=color_seq)

heat<-list()
mat_genes <- assay(vst_b)
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

#subset only groups of interest for heatmap
anno <- colData(vst_b_b)
for (i in names(genes2plot)){
  plot <- mat_genes[which(rownames(mat_genes) %in% genes2plot[[i]] ),]
  #annorow <- ifelse(is.na(plot$symbol),rownames(plot), plot$symbol)
  annorow <- rownames(plot)
  heat[[i]]<- pheatmap(plot[, c(rownames(anno[anno$group_new==as.character(contrasts[1,i]),]), rownames(anno[anno$group_new==as.character(contrasts[2,i]),]))], 
  scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst), show_rownames=F,annotation_colors=anno_colors
                       ) 
  pdf(file.path(PostDE.dir,i,"Heat_DEG_subgroups.pdf"))
  print(heat[[i]])
  dev.off()
    heat[[i]]<- pheatmap(plot[, c(rownames(anno[anno$group_new==as.character(contrasts[1,i]),]), rownames(anno[anno$group_new==as.character(contrasts[2,i]),]))], 
    scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),show_rownames=F,annotation_colors=anno_colors,
                       clustering_distance_rows="correlation") 
  pdf(file.path(PostDE.dir,i,"Heat_correlation_DEG_subgroups.pdf"))
  print(heat[[i]])
  dev.off()
}

#common heatmap
DEG <- unique(unlist(genes2plot))
mat_genes<- assay(vst_b)
plot <- mat_genes[which(rownames(mat_genes) %in% DEG ),]
set.seed(123)
p <- pheatmap(plot, scale="row", show_colnames=F,#labels_row=annorow,
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
#  GO[[i]] <- enrichGO(temp, OrgDb = org.Mm.eg.db,keyType="ENSEMBL",  
#                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "All")
#
# GO[[i]] <- setReadable(GO[[i]], OrgDb = org.Mm.eg.db)
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
  temp <- DEG_results_list[[i]][which(abs(DEG_results_list[[i]]$log2FoldChange)>1 & DEG_results_list[[i]]$padj <0.05),]
  table <- as.data.frame(table(temp$direction))
  labs <- rev(paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq))
i
  pdf(file.path(PostDE.dir,i,paste0("Direction_DEGs_sign_fc.pdf")), height=3.5, width=3.5)
  print(ggpie(table, "Freq", fill="Var1", palette=col[c(2,3)], label=labs, lab.pos = "in", main="DEG Analysis\nadj. p-value < 0.05 &\nabs. log2 fold-change > 1",  lab.font = c(5, "bold", "white")) + rremove("legend"))
  dev.off()

    temp <- DEG_results_list[[i]][which(DEG_results_list[[i]]$padj <0.05),]
  table <- as.data.frame(table(temp$direction))
  labs <- rev(paste0(table$Var1,"\n(", round((table$Freq/sum(table$Freq)*100)),"%)\n", table$Freq))

  pdf(file.path(PostDE.dir,i,paste0("Direction_DEGs_sign.pdf")), height=3.5, width=3.5)
  print(ggpie(table, "Freq", fill="Var1", palette=col[c(2,3)], label=labs, lab.pos = "in", main="DEG Analysis\nadj. p-value < 0.05",  lab.font = c(5, "bold", "white")) + rremove("legend"))
  dev.off()
}




