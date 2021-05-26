#Epigenetic Players 
#RNA analysis
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_4"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

library(DESeq2)
library(openxlsx)
library(dplyr)
library(ggpubr)
library("rafalib")
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(randomcoloR)

#load data
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
rld_b<- readRDS(file.path(results.dir,"rld_b_replicate.rds"))
epi_mouse <- readRDS("c010-datasets/External/2018-07-20-Coral/3_rep/data/MouseEpiPlayer.rds")

#make epi player "human"
epi_human <- epi_mouse
epi_human$Symbol <- toupper(rownames(epi_mouse))
rownames(epi_human)<- epi_human$Symbol 
epi_human$ENSEMBL<- mapIds(org.Hs.eg.db, keys= epi_human$Symbol , keytype ="SYMBOL", column = "ENSEMBL", multiVals = "first" )
epi_human$ENTREZ<- mapIds(org.Hs.eg.db, keys= epi_human$Symbol , keytype ="SYMBOL", column = "ENTREZID", multiVals = "first" )

#subset list
lapply(DEG_results_list,function(x)dim(x))
DEG_results_list_sub <- lapply(DEG_results_list, function(x){
    x <- x[which(rownames(x) %in% epi_human$ENSEMBL),]
    x
})
lapply(DEG_results_list_sub,function(x)dim(x))

for(i in names(DEG_results_list_sub)){
    write.table(DEG_results_list_sub[[i]], file.path(PostDE.dir, i , paste0("all_EpiPlayer_sub.tsv")), sep="\t", row.names=TRUE, quote=FALSE)
}

DEG_results_list_sub_sig <- lapply(DEG_results_list_sub, function(x){
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange)>1),]
    x
})
lapply(DEG_results_list_sub_sig,function(x)dim(x))

DEG_results_list_sub_sig_onlyFDR <- lapply(DEG_results_list_sub, function(x){
    x <- x[which(x$padj < 0.05) ,]
    x
})
lapply(DEG_results_list_sub_sig_onlyFDR,function(x)dim(x))

#plot heatmaps of those genes
rld_expr <- assay(rld_b)
rld_sub_sig <- lapply(DEG_results_list_sub_sig_onlyFDR, function(x){
    x <- rld_expr[which(rownames(rld_expr) %in% rownames(x)),]
})
annovst <- as.data.frame(colData(rld_b))[, c("genotype", "irradiation_treatment", "Replicate")] 
color_genotype<-randomColor(length(unique(colData(rld_b)$genotype)))
names(color_genotype)<- as.character(unique(colData(rld_b)$genotype))
color_treatment<-randomColor(length(unique(colData(rld_b)$irradiation_treatment)))
names(color_treatment)<- as.character(unique(colData(rld_b)$irradiation_treatment))
color_replicate<-randomColor(length(unique(colData(rld_b)$Replicate)))
names(color_replicate)<- as.character(unique(colData(rld_b)$Replicate))
anno_colors <- list(genotype=color_genotype, irradiation_treatment=color_treatment, Replicate=color_replicate)


# find non-complete elements
ids.to.remove <- sapply(rld_sub_sig, function(i) nrow(i) <= 0)
# remove found elements
rld_sub_sig <- rld_sub_sig[!ids.to.remove]

#get subset factor
contrasts <- as.data.frame(combn(as.character(unique(colData(dds)$group)), 2))
nam<-vector()
for (i in 1:length(contrasts)) {
nam[i]<- paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i]))
}
names(results)<- nam
names(contrasts)<- nam
#get symbol naming
lapply(rld_sub_sig, function(x){
    x$SYMBOL <-  <- mapIds(org.Hs.eg.db, keys= epi_human$Symbol , keytype ="SYMBOL", column = "ENTREZID", multiVals = "first" )
})
rld_sub_sig$Symbol


heat<-list()
for (i in names(rld_sub_sig)){
    heat[[i]]<- pheatmap(rld_sub_sig[[i]], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),
                       annotation_colors=anno_colors, 
                       clustering_distance_rows="correlation",
                       show_rownames=T) 
    pdf(file.path(PostDE.dir,i,"Heat_sig_EpiPlayer.pdf"))
    print(heat[[i]])
    dev.off() 

    heat[[i]]<- pheatmap(rld_sub_sig[[i]][, c(rownames(anno[anno$group==as.character(contrasts[1,i]),]), rownames(anno[anno$group==as.character(contrasts[2,i]),]))],
                        scale="row", show_colnames=F,
                       annotation_col=as.data.frame(annovst),
                       annotation_colors=anno_colors, clustering_distance_rows="correlation",
                       show_rownames=F) 
    pdf(file.path(PostDE.dir,i,"Heat_sig_EpiPlayer_subset.pdf"))
    print(heat[[i]])
    dev.off() 
}

#As barplot
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]
for (i in names(DEG_results_list_sub_sig)){
    DEG_results_list_sub_sig[[i]]$Symbol <- rownames(DEG_results_list_sub_sig[[i]])
    DEG_results_list_sub_sig[[i]]<-     DEG_results_list_sub_sig[[i]][order(DEG_results_list_sub_sig[[i]]$log2FoldChange, decreasing=TRUE),]
    DEG_results_list_sub_sig[[i]]$Change <- ifelse(DEG_results_list_sub_sig[[i]]$log2FoldChange >0, "blue", "red")
    pdf(file.path(PostDE.dir,i,"Bar_sig_EpiPlayer.pdf"))
    print(ggbarplot(DEG_results_list_sub_sig[[i]], x="Symbol", y="log2FoldChange",fill="Change",palette=col, main=i)
    + rremove("xlab")+ rremove("legend")+rotate_x_text(angle = 45))
    dev.off()
}









#look for specific genes
search <- "Dnmt"
lapply(DEG_results_list,function(x){
    x<- x[grep(search, rownames(x)), ]
    x
} )
assay(rld)[which(rownames(assay(rld))==search), ]

#plot all eigenetic players
rld_expr <- assay(rld)
rld_sub <-rld_expr[which(rownames(rld_expr) %in% rownames(epi_mouse)),]

annovst <- as.data.frame(colData(rld))[, c("Celltype", "Replicate")] 
color<-c("#0073C299", "#EFC00099", "#86868699")
names(color)<- c("pCAF", "sCAF", "NMF")
anno_colors <- list(Celltype=color)


temp<- pheatmap(rld_sub[,rownames(annovst[annovst$Celltype %in% c("pCAF","NMF"),])], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst)[annovst$Celltype %in% c("pCAF", "NMF"),],main=i, annotation_colors=anno_colors, show_rownames=T) 
pdf(file.path(PreDE.dir,"Heat_all_EpiMouse_pCAFvsNMF.pdf"))
print(temp)
dev.off() 
