#look at expression of tf target genes
#RNA 
#Irradiation
base_Irr.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data_Irr.dir <- file.path(base_Irr.dir, "data")
base_results_Irr.dir <- file.path(base_Irr.dir, "results")
results_Irr.dir<- file.path(base_results_Irr.dir , "tables")
PreDE_Irr.dir <- file.path(base_results_Irr.dir,"PreDE")
PostDE_Irr.dir <- file.path(base_results_Irr.dir,"PostDE")


#libraries
library(data.table)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(ggpubr)
library(dplyr)
library(ggrepel)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#load data
#RNA
#Irr
DEG_results_list<- readRDS(file.path(PostDE_Irr.dir,  "DEG_results_group_list.rds"))

#load tf binding site table
tf_binding <- fread("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_anl/FINAL_OUTPUT/extension100/wt_IrrvsID3-ko_Irr.all.allMotifs.tsv.gz")

#load translation table
translation_table <- fread("tools/diffTF/src/TF_Gene_TranslationTables/HOCOMOCO_v10/translationTable_hg19.csv")

#load diff tf final results
diffTF_res_Irr <- fread("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_anl/FINAL_OUTPUT/extension100/wt_IrrvsID3-ko_Irr.all.summary.tsv.gz")

#subset significant results
diffTF_res_Irr_sub <- diffTF_res_Irr[which(diffTF_res_Irr$pvalueAdj < 0.05),]
dim(diffTF_res_Irr_sub)
tf_binding_sub <- tf_binding[tf_binding$TF %in% diffTF_res_Irr_sub$TF]
length(unique(tf_binding_sub$TF))
#split by transcription factor
tf_binding_sub_split <- split(tf_binding_sub, tf_binding_sub$TF)
length(tf_binding_sub_split)
#make granges out of it 
tf_binding_sub_split_gr <- lapply(tf_binding_sub_split, function(x){
   y <- GRanges(
  seqnames = paste0(x$chr),
  ranges = IRanges(start = x$MSS,
                   end = x$MES
  )
)
mcols(y)<- x 
y
})
#annotate target genes
tf_binding_sub_split_anno <- lapply(tf_binding_sub_split_gr, function(x){
x<-annotatePeak(peak= x, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
x <- as.data.frame(x)
x
})

#prepare expression data
ID3KO_Irr_vs_WT_Irr <- DEG_results_list$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED



#final plotting 
#--->Promoter sub
#get average log2fold change
tf_target_gene_promoter_list <- list()
#subset list
tf_binding_sub_split_anno_promoter <- lapply(tf_binding_sub_split_anno, function(x){
    x <- x[grep("Promoter", x$annotation),]
    x
})
for(i in names(tf_binding_sub_split_anno_promoter)){
tf_target_gene_promoter_list[[i]] <-mean(ID3KO_Irr_vs_WT_Irr[ID3KO_Irr_vs_WT_Irr$ensembl %in% tf_binding_sub_split_anno_promoter[[i]]$ENSEMBL,]$log2FoldChange, na.rm=  TRUE) 
}
tf_target_gene_promoter <- as.data.frame(unlist(tf_target_gene_promoter_list))
tf_target_gene_promoter$TF <- rownames(tf_target_gene_promoter)
colnames(tf_target_gene_promoter)<- c("MeanTarget_log2FoldChange", "TF")

#merge tf fold and target gene expression
merge <- inner_join(tf_target_gene_promoter, diffTF_res_Irr, by="TF")
#merge_sub <- merge[merge$classification_q0.01_final %in% c("activator", "repressor", "undetermined"),]
merge$TF_class<- "undetermined"
# Currently hard-coded parameters
par.l = list()
# 1. Misc
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"
# 2. Statistical thresholds and values
par.l$significanceThresholds  = c(0.001, 0.01, 0.05,0.1,0.2) # p-value thresholds
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.1, 0.5, 0.8)
# 3. RNA-Seq specific
par.l$corMethod = "pearson" # Expression-peak count correlation method. As we quantile normalize now, should be pearson
par.l$regressionMethod = "glm" # for correlating RNA-Seq classification with TF activity
par.l$filter_minCountsPerCondition = 5 # For filtering RNA-seq genes, see Documentation
# 4. Volcano plot settings
par.l$maxTFsToLabel = 150 # Maximum Tfs to label in the Volcano plot
par.l$volcanoPlot_minDimensions  = 12
par.l$minPointSize = 0.3
par.l$plot_grayColor = "grey50"
par.l$colorCategories = c("activator" = "#4daf4a", "undetermined" = "black", "repressor" = "#e41a1c", "not-expressed" = "Snow3") # diverging, modified
par.l$colorConditions = c("#ef8a62", "#67a9cf") # Colors of the theme

#plot
g = ggplot() + geom_point(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference", fill="TF_class"),
    shape=21, stroke = 0.5) +
        scale_fill_manual("TF_class", values = par.l$colorCategories) + 
        geom_rect(aes(xmin = -Inf,xmax = 0,ymin = 0, ymax = Inf),
                    alpha = .3, fill = "gray", size = 0) +
        geom_rect(aes(xmin = 0,xmax = Inf,ymin = 0, ymax = -Inf),
                    alpha = .3, fill = "gray", size = 0) +
        geom_rect(aes(xmin = 0,xmax = Inf,ymin = 0, ymax = Inf),
                    alpha = .3, fill = par.l$colorConditions[1], size = 0) +
        geom_rect(aes(xmin = -Inf,xmax = 0,ymin = 0, ymax =-Inf),
                    alpha = .3, fill = par.l$colorConditions[2], size = 0) +
        geom_vline(xintercept=0, linetype = 2) +
        geom_hline(yintercept=0, linetype = 2) +
        ylab(paste0("Transcription factor activity","\nID3-KO Irr. vs WT Irr")) + 
        xlab(paste0("Mean target gene expression change\n[log2foldchange]","\nID3-KO Irr. vs WT Irr")) +
        # scale_alpha_manual(paste0("TF", " < ", significanceThresholdCur), values = c(alphaValueNonSign, 1), labels = c("no", "yes")) + 
        geom_label_repel(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference", label = "TF", fill="TF_class"),
                                    fontface = 'bold', color = 'white',size=2,
                                    segment.size = 0.3,# box.padding = unit(0.2, "lines"), max.iter = 5000,
                                    label.padding = unit(0.2, "lines"), # how thick is connectin line
                                    #nudge_y = 0.05, nudge_x = 0.05,  # how far from center points
                                    segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE) +
        theme_bw() + 
            theme(axis.text.x = element_text(size=rel(1)), axis.text.y = element_text(size=rel(1)), 
                axis.title.x = element_text(size=rel(1)), axis.title.y = element_text(size=rel(1)),
                legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(1))) +
        stat_cor(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference")) +
        geom_smooth(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference"),method = "lm",color='darkgray',se = FALSE)

pdf(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_anl/FINAL_OUTPUT/",  "MeanTargetGeneExpression_TFactivity_correlation_onlyPromoter.pdf"), width=7,height=7)
print(g)
dev.off()


#subset plot only e2f
#subset significant results
diffTF_res_Irr_sub <- diffTF_res_Irr[grep("E2F", diffTF_res_Irr$TF),]
diffTF_res_Irr_sub <- diffTF_res_Irr_sub[which(diffTF_res_Irr_sub$pvalueAdj < 0.05),]


dim(diffTF_res_Irr_sub)
tf_binding_sub <- tf_binding[tf_binding$TF %in% diffTF_res_Irr_sub$TF]
length(unique(tf_binding_sub$TF))
#split by transcription factor
tf_binding_sub_split <- split(tf_binding_sub, tf_binding_sub$TF)
length(tf_binding_sub_split)
#make granges out of it 
tf_binding_sub_split_gr <- lapply(tf_binding_sub_split, function(x){
   y <- GRanges(
  seqnames = paste0(x$chr),
  ranges = IRanges(start = x$MSS,
                   end = x$MES
  )
)
mcols(y)<- x 
y
})
#annotate target genes
tf_binding_sub_split_anno <- lapply(tf_binding_sub_split_gr, function(x){
x<-annotatePeak(peak= x, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
x <- as.data.frame(x)
x
})

#prepare expression data
ID3KO_Irr_vs_WT_Irr <- DEG_results_list$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED



#final plotting 
#--->Promoter sub
#get average log2fold change
tf_target_gene_promoter_list <- list()
#subset list
tf_binding_sub_split_anno_promoter <- lapply(tf_binding_sub_split_anno, function(x){
    x <- x[grep("Promoter", x$annotation),]
    x
})
for(i in names(tf_binding_sub_split_anno_promoter)){
tf_target_gene_promoter_list[[i]] <-mean(ID3KO_Irr_vs_WT_Irr[ID3KO_Irr_vs_WT_Irr$ensembl %in% tf_binding_sub_split_anno_promoter[[i]]$ENSEMBL,]$log2FoldChange, na.rm=  TRUE) 
}
tf_target_gene_promoter <- as.data.frame(unlist(tf_target_gene_promoter_list))
tf_target_gene_promoter$TF <- rownames(tf_target_gene_promoter)
colnames(tf_target_gene_promoter)<- c("MeanTarget_log2FoldChange", "TF")

#merge tf fold and target gene expression
merge <- inner_join(tf_target_gene_promoter, diffTF_res_Irr, by="TF")
#merge_sub <- merge[merge$classification_q0.01_final %in% c("activator", "repressor", "undetermined"),]
merge$TF_class<- "undetermined"
# Currently hard-coded parameters
par.l = list()
# 1. Misc
par.l$verbose = TRUE
par.l$log_minlevel = "INFO"
# 2. Statistical thresholds and values
par.l$significanceThresholds  = c(0.001, 0.01, 0.05,0.1,0.2) # p-value thresholds
par.l$classes_CohensD = c("small", "medium", "large", "very large")
par.l$thresholds_CohensD = c(0.1, 0.5, 0.8)
# 3. RNA-Seq specific
par.l$corMethod = "pearson" # Expression-peak count correlation method. As we quantile normalize now, should be pearson
par.l$regressionMethod = "glm" # for correlating RNA-Seq classification with TF activity
par.l$filter_minCountsPerCondition = 5 # For filtering RNA-seq genes, see Documentation
# 4. Volcano plot settings
par.l$maxTFsToLabel = 150 # Maximum Tfs to label in the Volcano plot
par.l$volcanoPlot_minDimensions  = 12
par.l$minPointSize = 0.3
par.l$plot_grayColor = "grey50"
par.l$colorCategories = c("activator" = "#4daf4a", "undetermined" = "black", "repressor" = "#e41a1c", "not-expressed" = "Snow3") # diverging, modified
par.l$colorConditions = c("#ef8a62", "#67a9cf") # Colors of the theme

#plot
g = ggplot() + geom_point(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference", fill="TF_class"),
    shape=21, stroke = 0.5) +
        scale_fill_manual("TF_class", values = par.l$colorCategories) + 
        geom_rect(aes(xmin = -Inf,xmax = 0,ymin = 0, ymax = Inf),
                    alpha = .3, fill = "gray", size = 0) +
        geom_rect(aes(xmin = 0,xmax = Inf,ymin = 0, ymax = -Inf),
                    alpha = .3, fill = "gray", size = 0) +
        geom_rect(aes(xmin = 0,xmax = Inf,ymin = 0, ymax = Inf),
                    alpha = .3, fill = par.l$colorConditions[1], size = 0) +
        geom_rect(aes(xmin = -Inf,xmax = 0,ymin = 0, ymax =-Inf),
                    alpha = .3, fill = par.l$colorConditions[2], size = 0) +
        geom_vline(xintercept=0, linetype = 2) +
        geom_hline(yintercept=0, linetype = 2) +
        ylab(paste0("Transcription factor activity","\nID3-KO Irr. vs WT Irr")) + 
        xlab(paste0("Mean target gene expression change\n[log2foldchange]","\nID3-KO Irr. vs WT Irr")) +
        # scale_alpha_manual(paste0("TF", " < ", significanceThresholdCur), values = c(alphaValueNonSign, 1), labels = c("no", "yes")) + 
        geom_label_repel(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference", label = "TF", fill="TF_class"),
                                    fontface = 'bold', color = 'white',size=2,
                                    segment.size = 0.3,# box.padding = unit(0.2, "lines"), max.iter = 5000,
                                    label.padding = unit(0.2, "lines"), # how thick is connectin line
                                    #nudge_y = 0.05, nudge_x = 0.05,  # how far from center points
                                    segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE) +
        theme_bw() + 
            theme(axis.text.x = element_text(size=rel(1)), axis.text.y = element_text(size=rel(1)), 
                axis.title.x = element_text(size=rel(1)), axis.title.y = element_text(size=rel(1)),
                legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(1))) +
        stat_cor(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference")) +
        geom_smooth(data = merge, aes_string("MeanTarget_log2FoldChange", "weighted_meanDifference"),method = "lm",color='darkgray',se = FALSE)

pdf(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_anl/FINAL_OUTPUT/",  "MeanTargetGeneExpression_TFactivity_correlation_onlyPromoter_onlyE2F.pdf"), width=5,height=5)
print(g)
dev.off()









#further investigate those genes
tf_target_gene_promoter_list_comp <- list()
tf_target_gene_promoter_list_comp_sub<- list()
for(i in names(tf_binding_sub_split_anno_promoter)){
tf_target_gene_promoter_list_comp[[i]] <-ID3KO_Irr_vs_WT_Irr[ID3KO_Irr_vs_WT_Irr$ensembl %in% tf_binding_sub_split_anno_promoter[[i]]$ENSEMBL,]
tf_target_gene_promoter_list_comp_sub[[i]] <- tf_target_gene_promoter_list_comp[[i]][which(tf_target_gene_promoter_list_comp[[i]]$padj < 0.05 & abs(tf_target_gene_promoter_list_comp[[i]]$log2FoldChange) > 0.5),]
}

write.table(tf_target_gene_promoter_list_comp_sub$E2F4, file.path("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_anl/FINAL_OUTPUT/",  "DEG_SigFC_E2F4_targets.txt"),
    quote=FALSE,row.names=FALSE, sep="\t")



#further investigate those genes
tf_target_gene_promoter_list_comp <- list()
tf_target_gene_promoter_list_comp_sub<- list()
for(i in names(tf_binding_sub_split_anno_promoter)){
tf_target_gene_promoter_list_comp[[i]] <-ID3KO_Irr_vs_WT_Irr[ID3KO_Irr_vs_WT_Irr$ensembl %in% tf_binding_sub_split_anno_promoter[[i]]$ENSEMBL,]
tf_target_gene_promoter_list_comp_sub[[i]] <- tf_target_gene_promoter_list_comp[[i]][which(tf_target_gene_promoter_list_comp[[i]]$padj < 0.05),]
}

write.table(tf_target_gene_promoter_list_comp_sub$E2F4, file.path("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_onlyE2F_anal/FINAL_OUTPUT/",  "DEG_Sig_E2F4_targets.txt"),
    quote=FALSE,row.names=FALSE, sep="\t")

