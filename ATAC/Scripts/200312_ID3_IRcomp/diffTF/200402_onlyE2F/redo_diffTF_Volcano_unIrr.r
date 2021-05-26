
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
tf_binding <- fread("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_unIrr_anal/FINAL_OUTPUT/extension100/wt_unIrrvsID3-ko_unIrr.all.allMotifs.tsv.gz")

#load translation table
translation_table <- fread("tools/diffTF/src/TF_Gene_TranslationTables/HOCOMOCO_v10/translationTable_hg19.csv")

#load diff tf final results
diffTF_res_Irr <- fread("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_unIrr_anal/FINAL_OUTPUT/extension100/wt_unIrrvsID3-ko_unIrr.all.summary.tsv.gz")

#select genes for highlighting
diffTF_res_Irr_plot <- as.data.frame(diffTF_res_Irr)
diffTF_res_Irr_plot$Plot <- NA
diffTF_res_Irr_plot[grep("E2F", diffTF_res_Irr_plot$TF),]$Plot  <- diffTF_res_Irr_plot[grep("E2F", diffTF_res_Irr$TF),]$TF 
diffTF_res_Irr_plot$logPadj <- -log10(diffTF_res_Irr_plot$pvalueAdj)
#merge_sub <- merge[merge$classification_q0.01_final %in% c("activator", "repressor", "undetermined"),]
diffTF_res_Irr_plot$TF_class<- "undetermined"
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
g = ggplot() + geom_point(data = diffTF_res_Irr_plot, aes_string("weighted_meanDifference", "logPadj", fill="TF_class"),
    shape=21, stroke = 0.3) +
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
        ylab(paste0("-log10(adjusted p-value)")) + 
        xlab(paste0("Transcription factor activity","\nID3-KO UnIrr. vs WT UnIrr")) +
        # scale_alpha_manual(paste0("TF", " < ", significanceThresholdCur), values = c(alphaValueNonSign, 1), labels = c("no", "yes")) + 
        geom_label_repel(data = diffTF_res_Irr_plot, aes_string("weighted_meanDifference", "logPadj", label = "Plot", fill="TF_class"),
                                    fontface = 'bold', color = 'white',size=2,
                                    segment.size = 0.3,# box.padding = unit(0.2, "lines"), max.iter = 5000,
                                    label.padding = unit(0.2, "lines"), # how thick is connectin line
                                    #nudge_y = 0.05, nudge_x = 0.05,  # how far from center points
                                    segment.alpha = .8, segment.color = par.l$plot_grayColor, show.legend = FALSE) +
        theme_bw() + 
            theme(axis.text.x = element_text(size=rel(1)), axis.text.y = element_text(size=rel(1)), 
                axis.title.x = element_text(size=rel(1)), axis.title.y = element_text(size=rel(1)),
                legend.title=element_text(size=rel(1)), legend.text=element_text(size=rel(1))) +
                rremove("legend")
        #stat_cor(data = diffTF_res_Irr_plot, aes_string("weighted_meanDifference", "logPadj")) +
        #geom_smooth(data = diffTF_res_Irr_plot, aes_string("weighted_meanDifference", "logPadj"),method = "lm",color='darkgray',se = FALSE)

pdf(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_unIrr_anal/FINAL_OUTPUT/",  "Redone_volcanoe_E2F_highlight.pdf"), width=5,height=5)
print(g)
dev.off()
