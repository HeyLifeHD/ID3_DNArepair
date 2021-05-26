#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(Rtsne)
#folder

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/repeats"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Data Read-in
#Data Read-in
#DESeq2 Analysis
dds_repeats <- readRDS(file = file.path(results.dir,"dds_repeats.rds"))
anno <- colData(dds_repeats)

#Extracting transformed values
rld <- rlog(dds_repeats, blind=FALSE)
#saveRDS(rld,file = file.path(results.dir,"rld.rds"))
#rld<-readRDS(file.path(results.dir,"rld.rds"))
vst_repeats <- vst(dds_repeats)
anno <- colData(vst_repeats)
anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment)
colData(vst_repeats)<- anno
saveRDS(vst_repeats,file = file.path(results.dir,"vst_repeats.rds"))

#function to plot PC1 and 2
pcaDatarld<- plotPCA(vst_repeats, intgroup =  colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld <- round(100* attr(pcaDatarld, "percentVar"))

pdf(file.path(PreDE.dir, "PCA12_repeats.pdf"), height = 4, width = 10)
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
pcaDatarld.alt<- plotPCA.alt(vst_repeats, intgroup = colnames(anno), returnData =TRUE, ntop=Inf)
percentvarrld.alt <- round(100* attr(pcaDatarld.alt, "percentVar"))


pdf(file.path(PreDE.dir, "PCA23_repeats.pdf"), height = 7, width = 10)
ggscatter(pcaDatarld.alt, x="PC2", y="PC3",size=3,
          color = "group_uncomplete2", shape = "genotype", legend = "right",
          ellipse = F ,mean.point = FALSE,label="replicate",repel=TRUE,
          star.plot = F,  xlab=(paste0("PC2: ", percentvarrld.alt[1], "% variance")), ylab=(paste0("PC3: ", percentvarrld.alt[2], "% variance")) )
dev.off()

#tsne
#tsne_model_1 = Rtsne(as.matrix(t(matvst[,1:500])), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
#d_tsne_1 = as.data.frame(tsne_model_1$Y)  


#extract eigenvectors of top 100 repeats
#ir.pca <- prcomp(assay(vst_repeats),
#                 center = TRUE,
                 scale. = TRUE) 
#ev<- as.data.frame(ir.pca$x)
#mv_ev_PC2 <- head(ev[order(abs(ev$PC2), decreasing=TRUE),],100)
#as.factor(rownames(mv_ev_PC2))

#Sample Clustering
dvst <- dist(t(assay(vst_repeats)))
hvst <- hclust(dvst)

pdf(file.path(PreDE.dir, "Clustering_vst_repeats.pdf"), height = 12, width = 7)
myplclust(hvst, labels=anno$group_complete, lab.col= as.numeric(as.factor(anno$group_uncomplete2)), main="", cex=0.5)
dev.off()

#Gene clustering: Heatmaps
#500 most variable repeats
topVarrepeatsvst<- head(order(rowVars(assay(vst_repeats)), decreasing=TRUE),10000)

matvst <- assay(vst_repeats)[topVarrepeatsvst,]
annovst <- as.data.frame(colData(vst_repeats)[,colnames(anno)])
annovst$replicate <- as.factor(annovst$replicate)
annovst <- annovst[, c("genotype", "tamoxifen_treatment","hdac_treatment", "replicate" )]
pdf(file.path(PreDE.dir,"Heatmap1000vst_noScale_repeats.pdf"),height= 7)
pheatmap(matvst, annotation_col=annovst ,scale="none", show_colnames=F, show_rownames=F)
dev.off()



#subset ltr
#repeat_list<- readRDS(file.path(results.dir, "repeat_list.rds"))
#counts_sum_repFamily <- repeat_list[["repFamily"]]
#counts_sum_repClass <- repeat_list[["repClass"]]
#counts_sum_repName <- repeat_list[["repName"]]
#counts_sum_repFamily_erv <- counts_sum_repFamily[grep("ERV", counts_sum_repFamily$repFamily),]
#counts_sum_repClass_ltr <- counts_sum_repClass[grep("LTR", counts_sum_repClass$repClass),]
#counts_sum_repName_ltr12 <- counts_sum_repName[grep("LTR12", counts_sum_repName$repName),]
vst_as <- assay(vst_repeats)
vst_LTR <- vst_as[grep("LTR12",rownames(vst_as)),]
#Plot LTR12
annovst <- as.data.frame(colData(vst_repeats)[anno$sample_name,])
annovst$replicate <- as.factor(annovst$replicate)
annovst <- annovst[, c("genotype", "tamoxifen_treatment","hdac_treatment", "replicate", "irradiation_treatment" )]

#pdf(file.path(PreDE.dir, "repName_Heatmap_log2_scaled.pdf"))
pheatmap(vst_LTR, 
                scale="row",
                annotation_col=annovst, 
                labels_row=rownames(vst_LTR),
                show_colnames = TRUE,show_rownames = TRUE,
                #annotation_colors = anno_colors,
                cluster_rows=F,fontsize = 5,
                #color=col_palette,
                filename= file.path(PreDE.dir, "repName_LTR_Heatmap_Scaled.pdf"),
                border_color = "grey"
                #cellwidth=5,
)
#dev.off()

