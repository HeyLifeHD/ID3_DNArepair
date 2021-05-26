#look at overlaps

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
library(UpSetR)


#folder
#folder
#libraries
library(openxlsx)
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#read data
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
rld_b<- readRDS(file.path(results.dir,"rld_b_replicate.rds"))

#select comparison 
name <- "IID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED"
result.df <- DEG_results_list[["ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED"]]



gene_name <- "MCPH1"
goi_in <- result.df[grep(gene_name,result.df$symbol),]$ensembl
names(goi_in)<- result.df[grep(gene_name,result.df$symbol),]$symbol
goi_in

result.df$padj <- round(result.df$padj, 3)
result.df$padj <- ifelse(result.df$padj < 0.001, "<0.001", result.df$padj)

#Gene of interest
dir.create(file.path(base_results.dir,"Genecounts", name), recursive=TRUE)
for(gene in names(goi_in)){ 
  pdf(file.path(base_results.dir,"Genecounts",name,paste0( names(goi_in[gene]),".pdf")), height = 5, width = 5)
  Count <- plotCounts(dds, gene = goi_in[gene],intgroup = c("genotype", "irradiation_treatment"), returnData = TRUE, normalized=TRUE)
  Count <- Count[Count$genotype !="MDC1",]
  print(ggboxplot(Count, x="genotype", y="count",merge=FALSE,combine=TRUE,fill="genotype",#yscale = "log2",# main = sapply(strsplit(names(goi_in[gene]), ".", fixed=TRUE), "[",2),
                  xlab="genotype",shape="genotype", ylab=paste0(names(goi_in[gene]),"\nnormalized gene counts"),  palette= rev(c("#D95F02",  "#1B9E77" ,  "#7570B3")),order=c("WT", "ID3_rescue", "ID3"),
                  #add = "jitter",
                  subtitle=paste0("ID3-ko UnIrr. vs wt UnIrr\n","log2fc = ", round(result.df[goi_in[gene],]$log2FoldChange,2), 
                              "; adj. p-value = ",(result.df[goi_in[gene],]$padj)#, #"(pvalue= ~",result.df[factors[gene],]$pvalue,")"  ,sep=" "),
                  #title=
                  #font.label = list(size = 11,color = "black")
                 ))
         +  rotate_x_text(angle = 45)+rremove("legend")+rremove("xlab"))
  dev.off()
}