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

#Epigenetic Players 
#RNA analysis
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
DEG_results_list_Irr_ID3<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
#rld_b_Irr<- readRDS(file.path(results.dir,"rld_b_replicate.rds"))

#folder
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/old/arid1a_project"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#Read in Data
DEG_results_list_Irr_arid1a<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#select data
deg_id3 <- DEG_results_list_Irr_ID3$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED
deg_arid1a <- DEG_results_list_Irr_arid1a$ARID1A_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED

#join data
merge <- left_join(deg_id3, deg_arid1a, by="ensembl")

#plot scatter
dir.create(file.path(base.dir, "ID3_Arid1a_comparison"))
pdf(file.path(base.dir, "ID3_Arid1a_comparison","ID3vsWT_Arid1AvsWT_Irr_Correlation.pdf"), height=3.5, width=3.5)
ggscatter(merge, x = "log2FoldChange.x", y = "log2FoldChange.y",xlab="Gene expression (log2 fold)\nID3-KD vs WT", 
ylab="Gene expression (log2 fold)\nARID1A-KD vs WT",shape = 16, size = .5,
   #color = "col" , palette=c("grey","#CD534CFF","#7AA6DCFF"), # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   #label="label", repel = TRUE, legend.title="P.adjusted",legend = "right",
  # cor.coef = TRUE, #font.label = c( "black"),# Add correlation coefficient. see ?stat_cor
   #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
   )+  stat_cor( label.y.npc = "top", label.x.npc = "left") + geom_vline(xintercept=0, linetype = 2) +geom_hline(yintercept=0, linetype = 2) + rremove("legend")
dev.off()

#look at overlap of sig.
deg_id3_sig <- deg_id3[which(deg_id3$padj < 0.05 & deg_id3$log2FoldChange >0.5),]
deg_arid1a_sig <- deg_arid1a[which(deg_arid1a$padj < 0.05& deg_arid1a$log2FoldChange >0.5),]

deg_id3_sig[deg_id3_sig$ensembl %in% deg_arid1a_sig$ensembl,]