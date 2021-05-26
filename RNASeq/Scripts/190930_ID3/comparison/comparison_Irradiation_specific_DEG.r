#comparison of DEG with and without irradiation
#folders
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/comparison"
dir.create(base.dir)
data_wIRR.dir <- file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep", "results", "PostDE")
data_woIRR.dir <- file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep", "results", "PostDE")

#load data
DEG_results_list_wIrr<- readRDS(file.path(data_wIRR.dir, "DEG_results_group_list.rds"))
DEG_results_list_woIrr<- readRDS(file.path(data_woIRR.dir, "DEG_results_group_list.rds"))

#create new list to compare
DEG_oI <- list(DEG_results_list_wIrr$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED,
DEG_results_list_woIrr$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED)
names(DEG_oI) <- c("ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED",
"ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED")
#subset list
DEG_oI<- lapply(DEG_oI, function(x){
    x <- as.data.frame(x)
    x$direction <- ifelse(x$log2FoldChange>0, "up", "down") 
    x$ENSEMBL <- NA
    x$ENSEMBL <- rownames(x)
    x
})
DEG_oI_sign<- lapply(DEG_oI, function(x){
    x <- x[which(x$padj < 0.05), ]
    x
})
DEG_oI_sign_fc<- lapply(DEG_oI, function(x){
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange)>1), ]
    x
})

#look at overlap between DEG
temp <- list(DEG_oI_sign$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_oI_sign$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
            DEG_oI_sign$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_oI_sign$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",],        
            DEG_oI_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[ DEG_oI_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
            DEG_oI_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[ DEG_oI_sign$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",]
)
names(temp)<- c("ID3_vs_WT_Irr_sign_up", "ID3_vs_WT_Irr_sign_down","ID3_vs_WT_UnIrr_sign_up", "ID3_vs_WT_UnIrr_sign_down")
temp <- lapply(temp, function(x){rownames(x)})
library(UpSetR)
pdf(file.path(base.dir, "Upset_overlap_sign_DEG.pdf"))
upset(fromList(temp), order.by = "freq")
dev.off()

#for fc and sign
temp <- list(DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
            DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",],        
            DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[ DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
            DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[ DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",]
)
names(temp)<- c("ID3_vs_WT_Irr_sign_up", "ID3_vs_WT_Irr_sign_down","ID3_vs_WT_UnIrr_sign_up", "ID3_vs_WT_UnIrr_sign_down")
temp <- lapply(temp, function(x){rownames(x)})
library(UpSetR)
pdf(file.path(base.dir, "Upset_overlap_sign_fc_DEG.pdf"))
upset(fromList(temp), order.by = "freq")
dev.off()


#look at correlation
merge <- left_join(DEG_oI$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED, DEG_oI$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED, by="ENSEMBL" )
merge$log2foldMean <- rowMeans(merge[,c("log2FoldChange.x", "log2FoldChange.y")])
merge<- merge[order(abs(merge$log2foldMean),decreasing=TRUE ),]
merge$label <- NA
merge[1:50,]$label<- merge[1:50,]$symbol.x

pdf(file.path(base.dir,"ID3_Irr_noIrr_Correlation_label.pdf"), height=3.5, width=3.5)
ggscatter(merge, x = "log2FoldChange.x", y = "log2FoldChange.y",xlab="Gene expression (log2 fold)\nID3-KD vs WT", 
ylab="Gene expression (log2 fold)\nsiMDC1 vs WT",shape = 16, size = .5,
   #color = "col" , palette=c("grey","#CD534CFF","#7AA6DCFF"), # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   label="label", repel = TRUE, legend.title="P.adjusted",legend = "right",
  # cor.coef = TRUE, #font.label = c( "black"),# Add correlation coefficient. see ?stat_cor
   #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
   )+  stat_cor( label.y.npc = "top", label.x.npc = "left") + geom_vline(xintercept=0, linetype = 2) +geom_hline(yintercept=0, linetype = 2) + rremove("legend")
dev.off()



#look at enrichment of specific genes
temp <- list(DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
            DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED[DEG_oI_sign_fc$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",],        
            DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[ DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="up",],
            DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED[ DEG_oI_sign_fc$ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED$direction=="down",]
)
names(temp)<- c("ID3_vs_WT_Irr_sign_up", "ID3_vs_WT_Irr_sign_down","ID3_vs_WT_UnIrr_sign_up", "ID3_vs_WT_UnIrr_sign_down")
temp <- lapply(temp, function(x){rownames(x)})
as.factor(temp$ID3_vs_WT_Irr_sign_up[!temp$ID3_vs_WT_Irr_sign_up %in% temp$ID3_vs_WT_UnIrr_sign_up])
as.factor(temp$ID3_vs_WT_Irr_sign_down[!temp$ID3_vs_WT_Irr_sign_down %in% temp$ID3_vs_WT_UnIrr_sign_down])
as.factor(temp$ID3_vs_WT_Irr_sign_down[temp$ID3_vs_WT_Irr_sign_down %in% temp$ID3_vs_WT_UnIrr_sign_down])
as.factor(temp$ID3_vs_WT_Irr_sign_up[temp$ID3_vs_WT_Irr_sign_up %in% temp$ID3_vs_WT_UnIrr_sign_up])
