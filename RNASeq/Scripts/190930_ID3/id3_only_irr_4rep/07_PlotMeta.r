#Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
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
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load files
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
rld_b<- readRDS(file = file.path(results.dir,"rld_b_replicate.rds"))

DEG_results_list<- readRDS(file.path(PostDE.dir,  "DEG_results_group_list.rds"))

#subset list
DEG_results_list_sig <- lapply(DEG_results_list, function(x){
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange)>0.5),]
    x
})
lapply(DEG_results_list_sig,function(x)dim(x))

DEG_results_list_sig_onlyFDR <- lapply(DEG_results_list, function(x){
    x <- x[which(x$padj < 0.05) ,]
    x
})
lapply(DEG_results_list_sig_onlyFDR,function(x)dim(x))

#look at sig, fc up and down in id3 vs wt
temp <-DEG_results_list_sig[["ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED"]]
up <-temp[which(temp$log2FoldChange>0),]
dim(up)
as.factor(up$ensembl)
down <-temp[which(temp$log2FoldChange<0),]
dim(down)
as.factor(down$ensembl)





#Plot Metascape
library(readxl)
library(ggpubr)
library(RColorBrewer)
#function
plotMeta <- function(excel_path, n=10, output_path, color="#CD534CFF", height=3.5, width=3.5){
dataset <- read_excel(excel_path, sheet=2)
datasetsub <- dataset[grep("Summary", dataset$GroupID), ]
datasetsub <- head(datasetsub, n)
datasetsub$LogP<- abs(datasetsub$LogP)
datasetsub$Description<- rev(datasetsub$Description)
datasetsub$LogP<- rev(datasetsub$LogP)
datasetsub$Term<- rev(datasetsub$Term)
datasetsub$description <- as.factor(1:nrow(datasetsub))
 p <- ggplot(datasetsub, aes(x=description)) +
geom_bar(aes(x=description, y=LogP, fill = "-log10(p.value)"), stat = "identity", width=0.25)
 p<-p+ annotate("text",x=as.integer(datasetsub$description)+0.35, y=0.01, 
                label=datasetsub$Description, size=4, hjust = 0)
 p <- p + scale_y_continuous( expand = c(0, 0), name = "-log10(p-value)")
 p <- p + scale_fill_manual(values = c( color))
 p<-p+coord_flip() 
 p<-p+ ggpubr:::theme_pubr()
 p<-p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=15))
 p<-p+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +rremove(c("legend"))+rremove(c("ylab"))
 p$labels$fill <- ""
 pdf(output_path, height=height, width=width)
 print(p)
 dev.off()
}
#directory
input.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep/results/PostDE/Metascape/ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED"

plotMeta(file.path(input.dir, "up", "metascape_result.xlsx"), file.path(input.dir, "up", "metascape_vis3.pdf"), n=20, width=7, height=7)
plotMeta(file.path(input.dir, "down", "metascape_result.xlsx"), file.path(input.dir, "down", "metascape_vis3.pdf"),color= "#67A9CF", n=20, width=7, height=7)


plotMeta(file.path(input.dir, "upregulated_EpiPlayers", "metascape_result.xlsx"), file.path(input.dir, "upregulated_EpiPlayers", "metascape_vis_top20_1.pdf"), n=20, width=5, height=5)
plotMeta(file.path(input.dir, "downregulated_EpiPlayers", "metascape_result.xlsx"), file.path(input.dir, "downregulated_EpiPlayers", "metascape_vis_top20_1.pdf"),color= "#67A9CF", n=20, width=5, height=5)


 
 
 
 