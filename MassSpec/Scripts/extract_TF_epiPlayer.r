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
library(readxl)
library(data.table)

#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#data
#epiplayer
epi_human <- fread("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/external_data/HumanEpiPlayer.csv")
#TF
TF_list<-read.table(file.path("c010-datasets/External/2018-07-20-Coral/" , "TF_list.txt"), sep="\t", header=TRUE)
TF_list <- TF_list[, c("Symbol", "Gene.name")]
TF_list$Symbol<-as.character(TF_list$Symbol)
TF_list$Symbol <- toupper(TF_list$Symbol)

#read in excel sheets
Irr_vs_UT_15min <- as.data.frame(read_excel(file.path(base.dir, "IP_MS", "volcano", "200130_15min_vs_UT_Table.xlsx" ), sheet=1))
Irr_vs_UT_1h <- as.data.frame(read_excel(file.path(base.dir, "IP_MS", "volcano", "200130_1h_vs_UT_Table.xlsx" ), sheet=1))

ms <- list(Irr_vs_UT_15min, Irr_vs_UT_1h )
names(ms)<- c("Irr_vs_UT_15min", "Irr_vs_UT_1h" )
#adjust names
ms <- lapply(ms, function(x){
    x$gene_name <- sapply(strsplit(x$"Gene names" , ";", fixed=TRUE), "[", 1)
    x
})

#subset based on present in epi
ms_sub_epi <- lapply(ms, function(x){
    x <- x[which(x$gene_name %in% epi_human$gene.name), ]
    x
})

for(i in names(ms_sub_epi)){
    write.table(ms_sub_epi[[i]], file.path(base.dir, "IP_MS",paste0("EpiPlayer_overlap_",i, ".txt")), sep="\t")
}

#subset based on present in TF
ms_sub_tf <- lapply(ms, function(x){
    x <- x[which(x$gene_name %in% TF_list$Symbol), ]
    x
})
for(i in names(ms_sub_epi)){
    write.table(ms_sub_epi[[i]], file.path(base.dir, "IP_MS",paste0("TF_overlap_",i, ".txt")), sep="\t")
}
