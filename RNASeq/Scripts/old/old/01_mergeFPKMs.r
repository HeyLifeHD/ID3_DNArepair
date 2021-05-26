
#Directories
#NGS-machine
.libPaths("/home/c010-ngs/R/x86_64-pc-linux-gnu-library/3.4_alt")
base.dir <- "/home/c010-ngs/c010-datasetshey/Internal/2018-Ali/RNASeq/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")

#VM
base.dir <- "/C010-datasets/Internal/2018-Ali/RNASeq/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")

#PC-09
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")
temp.dir <- "/home/heyj/rTemp"

#libraries
library(dplyr)
#load FPKM files
files <- list.files(file.path(base.dir, "181114_analysis_rep1", "stringtieFPKM"), pattern=".txt")
file_list <- as.vector(file.path(base.dir, "181114_analysis_rep1", "stringtieFPKM", file_list))
FPKM_list <- list(NULL)
for (i in 1:length(file_list)) {
    FPKM_list[[i]]<- read.table(file_list[i], sep="\t", header=TRUE)
}
#subset columns
FPKM_list_sub <- lapply(FPKM_list, function(x) x[, -c(2:8)])
for(i in 1:length(FPKM_list_sub)){
    FPKM_list_sub[[i]]$Gene.ID<-  sapply(strsplit(as.character(FPKM_list_sub[[i]]$Gene.ID) ,".", fixed=TRUE),`[`, 1)
}
#merge files
FPKM_df_sub<-Reduce(function(x, y) full_join(x, y, by="Gene.ID"), FPKM_list_sub)
#rename columns
colnames(FPKM_df_sub)[c(2:9)]<-  sapply(strsplit(as.character(files) ,".", fixed=TRUE),`[`, 1)
#write table
write.table(FPKM_df_sub, file.path(base.dir, "181114_analysis_rep1", "stringtieFPKM", "FPKM_merged_rep1.txt"), quote=TRUE, sep="\t", col.names=TRUE )