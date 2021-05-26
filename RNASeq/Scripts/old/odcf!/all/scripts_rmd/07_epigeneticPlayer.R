#overlap with epigenetic players
#directories
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
library(org.Hs.eg.db)
library(biomaRt)
#packages
library(DESeq2)
library(clusterProfiler)

#load files
dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
vst_genes<- readRDS(file = file.path(results.dir,"vst_genes_batchRem_rep.rds"))

DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#load epigenetic players
epi_player <- read.csv(file.path(data.dir, "HumanEpiPlayer.csv"), header=TRUE, stringsAsFactors = FALSE)
epi_player$ensembl<- mapIds(org.Hs.eg.db, keys=epi_player$gene.name, keytype ="SYMBOL", column = "ENSEMBL", multiVals = "first" )

#find overlap between epi player and genes
epi_player_ol <- lapply(DEG_results_list, function(x){
  ol <- x[which(x$ensembl %in% epi_player$ensembl),]
  ol
})

#find genes which have an sign change in gene expression
epi_player_ol_sign <- lapply(epi_player_ol, function(x){
  x <- as.data.frame(x)
  x <- x[which(x$padj<0.1), ]
  x
})

#add reader writer information
epi_player_ol_sign<- lapply(epi_player_ol_sign, function(x){
  x <- left_join(x, epi_player, by=c("ensembl"="ensembl"))
  x
})

#export them as a table
for (i in names(epi_player_ol_sign)){
  if (nrow(epi_player_ol_sign[[i]])>0){
    write.table(epi_player_ol_sign[[i]], file.path(PostDE.dir,i,"epigenetic_player.txt"),row.names=FALSE)
  }
}











#analysis
cutoff <- 0.05
genes2plot <- lapply(DEG_results_list, function(x){
  x <- x[which(x$padj < cutoff),]
  x <- x$ensembl
  x
})
# find non-complete elements
ids.to.remove <- sapply(genes2plot, function(i) length(i) <= 0)
# remove found elements
genes2plot <- genes2plot[!ids.to.remove]
#do overrepresentation analysis with subsetted genes
res <- list()
for (i in names(genes2plot)){
  temp <- DEG_results_list[[i]]
  temp <- temp[!is.na(temp$padj),]
  res[[i]] <- goseq_wrapper(temp$ensembl, genes2plot[[i]])
  plotgo(res[[i]], bn= paste0(PostDE.dir, "/", i,"/OR"))
}
