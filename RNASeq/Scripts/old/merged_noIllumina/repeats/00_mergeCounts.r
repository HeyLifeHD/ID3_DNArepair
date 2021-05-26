#libraries
library(data.table)

#folder

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/repeats"
data.dir <- file.path(base.dir, "data")
featureCounts.dir <-"/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/results/repeats/"

base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#load data
files <- dir(path = featureCounts.dir, full.names = T,
             pattern = ".mdup.bam_repeat_counts.txt")

featureCounts <- list()
for (cov_file in files){
  featureCounts[[cov_file]] <- fread(input = sprintf(cov_file))
}
#check for same length of files
lapply(featureCounts, function(x)length(x))
featureCounts_sub<- featureCounts[-grep("summary",names(featureCounts))]
lapply(featureCounts_sub, function(x)length(x))

#subset counts and create dataframe
counts <- lapply(featureCounts_sub, function(x){
    y <- as.data.frame(x[,c(7)])
    y 
})
counts <- do.call("cbind", counts)
#rename counts
sample_names <- sapply(strsplit(colnames(counts) ,"_merged", fixed=TRUE),`[`, 1)
colnames(counts)<- sample_names
#reorder counts
anno <- readRDS(file = file.path("/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/data","anno.rds"))
anno$sample_name <- as.character(anno$sample_name)
counts<- counts[,anno$sample_name]

colnames(counts) == anno$sample_name

saveRDS(counts , file.path(featureCounts.dir, "merged_repeat_counts.rds"))
saveRDS(as.data.frame(featureCounts[[1]][,1:6]) , file.path(featureCounts.dir, "merged_repeat_anno.rds"))