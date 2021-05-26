#Comparison of ID3 knockdown and MDC1 siRNA
#Epigenetic Players 
#RNA analysis
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

library(DESeq2)
library(openxlsx)
library(dplyr)
library(ggpubr)
library("rafalib")
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(randomcoloR)
library(RColorBrewer)
library(dplyr)

#load datax
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
rld_b<- readRDS(file.path(results.dir,"rld_b_replicate.rds"))

#look at overlap of differentially expressed genes
DEG_results_list<- lapply(DEG_results_list, function(x){
    x <- as.data.frame(x)
    x$direction <- ifelse(x$log2FoldChange>0, "up", "down") 
    x$ENSEMBL <- NA
    x$ENSEMBL <- rownames(x)
    x
})
DEG_results_list_sign<- lapply(DEG_results_list, function(x){
    x <- x[which(x$padj < 0.05), ]
    x
})
DEG_results_list_sign_fc<- lapply(DEG_results_list, function(x){
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange)>0.5), ]
    x
})

#look at overlap between DEG
dir.create(file.path(base.dir,"MDC1_ID3_overlap"))
temp <- list(DEG_results_list_sign$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
        DEG_results_list_sign$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",],
        DEG_results_list_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
        DEG_results_list_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",]
)
names(temp)<- c("MDC1_vs_WT_sign_up", "MDC1_vs_WT_sign_down","ID3_vs_WT_sign_up", "ID3_vs_WT_sign_down")
temp <- lapply(temp, function(x){rownames(x)})
library(UpSetR)
pdf(file.path(base.dir,"MDC1_ID3_overlap", "Upset_overlap_sign_DEG.pdf"))
upset(fromList(temp), order.by = "freq")
dev.off()

#for fc and sign
dir.create(file.path(base.dir,"MDC1_ID3_overlap"))
temp <- list(DEG_results_list_sign_fc$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign_fc$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
        DEG_results_list_sign_fc$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign_fc$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",],
        DEG_results_list_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
        DEG_results_list_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_results_list_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",]
)
names(temp)<- c("MDC1_vs_WT_sign_up", "MDC1_vs_WT_sign_down","ID3_vs_WT_sign_up", "ID3_vs_WT_sign_down")
temp <- lapply(temp, function(x){rownames(x)})
library(UpSetR)
pdf(file.path(base.dir,"MDC1_ID3_overlap", "Upset_overlap_sign_fc_DEG.pdf"), height=4, width=4)
upset(fromList(temp), order.by = "freq")
dev.off()


#look at correlation
merge <- left_join(DEG_results_list$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED, DEG_results_list$MDC1_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED, by="ENSEMBL" )
merge$log2foldMean <- rowMeans(merge[,c("log2FoldChange.x", "log2FoldChange.y")])
merge<- merge[order(abs(merge$log2foldMean),decreasing=TRUE ),]
merge$label <- NA
merge[1:50,]$label<- merge[1:50,]$symbol.x

pdf(file.path(base.dir, "MDC1_ID3_overlap","MDC1vsWT_ID3_vsWT_Correlation.pdf"), height=3.5, width=3.5)
ggscatter(merge, x = "log2FoldChange.x", y = "log2FoldChange.y",xlab="Gene expression (log2 fold)\nID3-KD vs WT", 
ylab="Gene expression (log2 fold)\nsiMDC1 vs WT",shape = 16, size = .5,
   #color = "col" , palette=c("grey","#CD534CFF","#7AA6DCFF"), # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   #label="label", repel = TRUE, legend.title="P.adjusted",legend = "right",
  # cor.coef = TRUE, #font.label = c( "black"),# Add correlation coefficient. see ?stat_cor
   #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
   )+  stat_cor( label.y.npc = "top", label.x.npc = "left") + geom_vline(xintercept=0, linetype = 2) +geom_hline(yintercept=0, linetype = 2) + rremove("legend")
dev.off()