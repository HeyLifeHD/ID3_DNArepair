#libraries
library(DESeq2)
#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")

base_results.dir <- file.path(base.dir, "results")
#dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
#dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
#dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
#dir.create(PostDE.dir)

intern.dir <- file.path(base.dir, "intern")
dir.create(intern.dir)

#save data for interns
#load data
counts<- readRDS(file.path(data.dir, "gene_counts.rds"))

#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno <- colData(dds)


#seperate  projects
id3_project <- c("ID3", "WT")
arid1a_project <- c("ARID1A", "WT")
#id3
counts_id3 <- counts[, rownames(anno[anno$genotype %in% id3_project,])]
dds_id3 <- dds[, rownames(anno[anno$genotype %in% id3_project,])]
anno_id3 <- anno[ rownames(anno[anno$genotype %in% id3_project,]),]

saveRDS(counts_id3,file.path(intern.dir, "counts_id3.rds"))
saveRDS(dds_id3,file.path(intern.dir, "dds_id3.rds"))
saveRDS(anno_id3,file.path(intern.dir, "anno_id3.rds"))

#arid1a
counts_arid1a <- counts[, rownames(anno[anno$genotype %in% arid1a_project,])]
dds_arid1a <- dds[, rownames(anno[anno$genotype %in% arid1a_project,])]
anno_arid1a <- anno[ rownames(anno[anno$genotype %in% arid1a_project,]),]

saveRDS(counts_arid1a,file.path(intern.dir, "counts_arid1a.rds"))
saveRDS(dds_arid1a,file.path(intern.dir, "dds_arid1a.rds"))
saveRDS(anno_arid1a,file.path(intern.dir, "anno_arid1a.rds"))