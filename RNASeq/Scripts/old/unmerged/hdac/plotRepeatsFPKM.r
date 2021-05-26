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
repeats.dir <- file.path(base_results.dir, "repeats")
#Data Read-in
#DESeq2 Analysis
dds <- readRDS(file = file.path(results.dir,"dds.rds"))
anno_dds <- colData(dds)

repeats=readRDS(c("/bigdisk/Nanopore/raw_data/repeats/repeats.RDS"))

count_list<- readRDS("c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged/results/repeats/merged_counts_repeats.rds")

#normalize
counts_norm <- count_list
counts_norm<- sapply(counts_norm,function(x) as.numeric(x))
#calculate cpm
counts_norm<- counts_norm/colSums(counts_norm)
counts_norm <- counts_norm*10^6
#combine 
counts_norm<- cbind(as.data.frame(repeats),counts_norm)

#as data.table
#setDT(counts_norm)
#summarize repeat
#counts_melt_repFamily <- data.table::melt(counts_norm[,c(samples, 'repFamily'), with = FALSE])[,sum(value), by=,.(repFamily)]
#sum_counts_repFamily <- (data.table::dcast(counts_melt_repFamily, variable ~ repFamily))

counts_sum_repFamily <- aggregate(cbind(counts_norm[,19:ncol(counts_norm)]), by=list(repFamily=counts_norm$repFamily), FUN=sum)
saveRDS(counts_sum_repFamily, file.path(repeats.dir,"cpm_sum_repFamily.rds" ))
counts_sum_repClass <- aggregate(cbind(counts_norm[,19:ncol(counts_norm)]), by=list(repClass=counts_norm$repClass), FUN=sum)
saveRDS(counts_sum_repClass, file.path(repeats.dir,"cpm_sum_repClass.rds" ))
counts_sum_repName <- aggregate(cbind(counts_norm[,19:ncol(counts_norm)]), by=list(repName=counts_norm$repName), FUN=sum)
saveRDS(counts_sum_repName, file.path(repeats.dir,"cpm_sum_repName.rds" ))

#subset ltr
counts_sum_repFamily_erv <- counts_sum_repFamily[grep("ERV", counts_sum_repFamily$repFamily),]
counts_sum_repClass_ltr <- counts_sum_repClass[grep("LTR", counts_sum_repClass$repClass),]
counts_sum_repName_ltr12 <- counts_sum_repName[grep("LTR12", counts_sum_repName$repName),]

#Plot LTR12
anno_heat <- as.data.frame(colData(dds)[,c( "replicate", "genotype", "dox_treatment", "hdac_treatment", "irr_treatment")])

pdf(file.path(repeats.dir, "LTR12_heatmap_noScale_log2FPKM.pdf"))
pheatmap(log2(counts_sum_repName_ltr12[,anno$sample_name]),
                 scale="none",
                 annotation_col=anno_heat,
                 labels_row=counts_sum_repName_ltr12$repName,
                 show_colnames = TRUE,show_rownames = TRUE,
                 #annotation_colors = anno_colors,
                 cluster_rows=F,fontsize = 5,
                 #color=col_palette,
                 border_color = "grey"
                 #cellwidth=5,
 )
dev.off()

pdf(file.path(repeats.dir, "LTR12_heatmap_Scale_log2FPKM.pdf"))
pheatmap(log2(counts_sum_repName_ltr12[,anno$sample_name]),
                 scale="row",
                 annotation_col=anno_heat,
                 labels_row=counts_sum_repName_ltr12$repName,
                 show_colnames = TRUE,show_rownames = TRUE,
                 #annotation_colors = anno_colors,
                 cluster_rows=F,fontsize = 5,
                 #color=col_palette,
                 border_color = "grey"
                 #cellwidth=5,
 )
dev.off()
#get most variant rep
topVar<- head(order(rowVars(as.matrix(counts_sum_repName[, anno$sample_name])), decreasing=TRUE),100)

pdf(file.path(repeats.dir, "repeat_top100_heatmap_noScale_log2FPKM.pdf"))
pheatmap(counts_sum_repName[topVar,anno$sample_name], 
                scale="none",
                 annotation_col=anno_heat,
                labels_row=counts_sum_repName[topVar,"repName"],
                show_colnames = TRUE,show_rownames = TRUE,
                #annotation_colors = anno_colors,
                cluster_rows=F,fontsize = 5,
                #color=col_palette,
                border_color = "grey"
                #cellwidth=5,
)
dev.off()

pdf(file.path(repeats.dir, "repeat_top100_heatmap_Scale_log2FPKM.pdf"))
pheatmap(counts_sum_repName[topVar,anno$sample_name], 
                scale="row",
                 annotation_col=anno_heat,
                labels_row=counts_sum_repName[topVar,"repName"],
                show_colnames = TRUE,show_rownames = TRUE,
                #annotation_colors = anno_colors,
                cluster_rows=F,fontsize = 5,
                #color=col_palette,
                border_color = "grey"
                #cellwidth=5,
)
dev.off()

#all no subset
pdf(file.path(repeats.dir, "repeat_all_heatmap_noScale_log2FPKM.pdf"))
pheatmap(counts_sum_repName[,anno$sample_name], 
                scale="none",
                annotation_col=anno_heat, 
                labels_row=counts_sum_repName[,"repName"],
                show_colnames = TRUE,show_rownames = TRUE,
                #annotation_colors = anno_colors,
                cluster_rows=T,fontsize = 5,
                #color=col_palette,
                border_color = "grey"
                #cellwidth=5,
)
dev.off()