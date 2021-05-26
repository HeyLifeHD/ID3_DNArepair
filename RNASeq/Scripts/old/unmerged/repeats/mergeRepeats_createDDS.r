#subset hddac sample set
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(rtracklayer)

#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
counts<- readRDS(file.path(data.dir, "gene_counts.rds"))

#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno_dds <- colData(dds)


#load repeat counts
repeat_counts <- read.table("c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged/results/repeats/repeat_counts.txt", header=TRUE, stringsAsFactors=FALSE)
repeat_anno <- repeat_counts[, c(1:6 )]
repeat_counts<- repeat_counts[, (-c(1:6 ))]

#annotation
anno <- read.table(file.path(data.dir, "annotation_unmerged.csv"), sep=";", stringsAsFactors=FALSE)
anno$sample_name <- paste0(anno$V1, "_", anno$V2)
colnames(anno)<- c("PID", "Treatment", "ID", "sample_name")
#modify anno
anno$celltype <- sapply(strsplit(anno$sample_name ,"_", fixed=TRUE),`[`, 2)
anno$genotype <- sapply(strsplit(anno$sample_name ,"_", fixed=TRUE),`[`, 3)
anno$replicate <- sapply(strsplit(anno$sample_name ,"_", fixed=TRUE),`[`, 4)


#data
colnames(repeat_counts)<- sapply(strsplit(colnames(repeat_counts) ,"STAR.", fixed=TRUE),`[`, 2)
colnames(repeat_counts)<- sapply(strsplit(colnames(repeat_counts) ,"Aligned.", fixed=TRUE),`[`, 1)
colnames(repeat_counts) <- gsub(".", "-", colnames(repeat_counts), fixed=TRUE)
#sum up unmerged samples
count_list <- list()
for (i in unique(anno$sample_name)){
    temp <- repeat_counts[,anno[anno$sample_name == i,]$ID ]
    if(class(temp)=="integer"){
        count_list[[i]]<- data.frame(temp)
       
    } else {
         count_list[[i]]<- data.frame(rowSums(temp))
    }
}
#count_list <- count_list[-1]
count_list <- do.call( "cbind", count_list)
colnames(count_list)<- unique(anno$sample_name)
saveRDS(count_list, "c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged/results/repeats/merged_counts_repeats.rds")
