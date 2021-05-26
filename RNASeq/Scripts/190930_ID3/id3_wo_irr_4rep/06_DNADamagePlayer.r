#Epigenetic Players 
#RNA analysis
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
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

#load dna damage lists
#load term2 gene list of damage
t2g <- read.csv2(file.path("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190531_odcf_NoIll/","data", "updated_DNA_damage_player.csv"))
t2g$ENSEMBL<-as.character(t2g$ensembl) 
t2g$Pathway<-as.character(t2g$Pathway) 
rownames(t2g)<- t2g$ensembl
t2g$symbol<- mapIds(org.Hs.eg.db, keys= t2g$ENSEMBL , keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )

#subset list
lapply(DEG_results_list,function(x)dim(x))
DEG_results_list_sub <- lapply(DEG_results_list, function(x){
    x <- x[which(rownames(x) %in% t2g$ENSEMBL),]
    x
})
lapply(DEG_results_list_sub,function(x)dim(x))

for(i in names(DEG_results_list_sub)){
    write.table(DEG_results_list_sub[[i]], file.path(PostDE.dir, i , paste0("all_DNAdamage_sub.tsv")), sep="\t", row.names=TRUE, quote=FALSE)
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
annovst <- as.data.frame(colData(rld_b))[, c("genotype", "Replicate")] 
color_genotype<-RColorBrewer::brewer.pal(n = length(unique(colData(rld_b)$genotype)), name = 'Dark2')
names(color_genotype)<- as.character(unique(colData(rld_b)$genotype))
color_replicate<-randomColor(length(unique(colData(rld_b)$Replicate)))
names(color_replicate)<- as.character(unique(colData(rld_b)$Replicate))
anno_colors <- list(genotype=color_genotype,  Replicate=color_replicate)


# find non-complete elements
ids.to.remove <- sapply(rld_sub_sig, function(i) nrow(i) <= 1)
# remove found elements
rld_sub_sig <- rld_sub_sig[!ids.to.remove]
#remove numbber one !!!!!!!!!!!!!
rld_sub_sig[[1]]<- NULL

#get subset factor
contrasts <- as.data.frame(combn(as.character(unique(colData(dds)$group)), 2))
nam<-vector()
for (i in 1:length(contrasts)) {
nam[i]<- paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i]))
}
names(results)<- nam
names(contrasts)<- nam
#get symbol naming
rld_sub_sig<-lapply(rld_sub_sig, function(x){
    x <- as.data.frame(x)
    x$ENSEMBL <- rownames(x)
    x$SYMBOL <- mapIds(org.Hs.eg.db, keys= x$ENSEMBL , keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )
    rownames(x)<- x$SYMBOL
    x
})
#get pathway annotaiton
rld_sub_sig<-lapply(rld_sub_sig, function(x){
    x<-left_join(x, t2g)
    rownames(x)<- x$SYMBOL
    x
})

#plot
anno <- colData(rld_b)
heat<-list()
for (i in names(rld_sub_sig)){
    anno_row <- rld_sub_sig[[i]][,c("Pathway"), drop=F]

    heat[[i]]<- pheatmap(rld_sub_sig[[i]][,rownames(anno)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),fontsize_row=3,annotation_row=anno_row,
                       annotation_colors=anno_colors, 
                       clustering_distance_rows="correlation",
                       show_rownames=T) 
    pdf(file.path(PostDE.dir,i,"Heat_sig_DNAdamage.pdf"))
    print(heat[[i]])
    dev.off() 

    heat[[i]]<- pheatmap(rld_sub_sig[[i]][, c(rownames(anno[anno$group==as.character(contrasts[1,i]),]), rownames(anno[anno$group==as.character(contrasts[2,i]),]))],
                        scale="row", show_colnames=F,
                       annotation_col=as.data.frame(annovst),fontsize_row=3,annotation_row=anno_row,
                       annotation_colors=anno_colors, clustering_distance_rows="correlation",
                       show_rownames=T) 
    pdf(file.path(PostDE.dir,i,"Heat_sig_DNAdamage_subset.pdf"))
    print(heat[[i]])
    dev.off() 
}

#As barplot
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]
for (i in names(DEG_results_list_sub_sig_onlyFDR)){
    DEG_results_list_sub_sig_onlyFDR[[i]]$Symbol <- DEG_results_list_sub_sig_onlyFDR[[i]]$symbol
    DEG_results_list_sub_sig_onlyFDR[[i]]<-     DEG_results_list_sub_sig_onlyFDR[[i]][order(DEG_results_list_sub_sig_onlyFDR[[i]]$log2FoldChange, decreasing=TRUE),]
    DEG_results_list_sub_sig_onlyFDR[[i]]$Change <- ifelse(DEG_results_list_sub_sig_onlyFDR[[i]]$log2FoldChange >0, "blue", "red")
    pdf(file.path(PostDE.dir,i,"Bar_sig_DNAdamage.pdf"), width=16)
    print(ggbarplot(DEG_results_list_sub_sig_onlyFDR[[i]], x="Symbol", y="log2FoldChange",fill="Change",palette=col)
    + rremove("xlab")+ rremove("legend")+rotate_x_text(angle = 45))
    dev.off()
}








#Same analysis for smaller DNA damage list
#load dna damage lists
#load term2 gene list of damage
t2g <- read.table(file.path( "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data","data", "DRGTerm2Gene.csv"), sep=";", stringsAsFactors = F)
colnames(t2g)<- c("Pathway", "symbol")
t2g$ensembl<- mapIds(org.Hs.eg.db, keys= t2g$symbol , keytype ="SYMBOL", column = "ENSEMBL", multiVals = "first" )
t2g$ENSEMBL<-as.character(t2g$ensembl) 
t2g$Pathway<-as.character(t2g$Pathway) 
#rownames(t2g)<- t2g$ensembl

#subset list
lapply(DEG_results_list,function(x)dim(x))
DEG_results_list_sub <- lapply(DEG_results_list, function(x){
    x <- x[which(rownames(x) %in% t2g$ENSEMBL),]
    x
})
lapply(DEG_results_list_sub,function(x)dim(x))

for(i in names(DEG_results_list_sub)){
    write.table(DEG_results_list_sub[[i]], file.path(PostDE.dir, i , paste0("all_DNAdamage_sub.tsv")), sep="\t", row.names=TRUE, quote=FALSE)
}

DEG_results_list_sub_sig <- lapply(DEG_results_list_sub, function(x){
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange)>.5),]
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
annovst <- as.data.frame(colData(rld_b))[, c("genotype", "Replicate")] 
color_genotype<-randomColor(length(unique(colData(rld_b)$genotype)))
names(color_genotype)<- as.character(unique(colData(rld_b)$genotype))
color_replicate<-randomColor(length(unique(colData(rld_b)$Replicate)))
names(color_replicate)<- as.character(unique(colData(rld_b)$Replicate))
anno_colors <- list(genotype=color_genotype,  Replicate=color_replicate)


# find non-complete elements
ids.to.remove <- sapply(rld_sub_sig, function(i) nrow(i) <= 1)
# remove found elements
rld_sub_sig <- rld_sub_sig[!ids.to.remove]
#remove numbber one !!!!!!!!!!!!!
#rld_sub_sig[[1]]<- NULL

#get subset factor
contrasts <- as.data.frame(combn(as.character(unique(colData(rld_b)$group)), 2))
nam<-vector()
for (i in 1:length(contrasts)) {
nam[i]<- paste0( as.character(contrasts[1,i]), "_vs_",  as.character(contrasts[2,i]))
}
names(contrasts)<- nam
#get symbol naming
rld_sub_sig<-lapply(rld_sub_sig, function(x){
    x <- as.data.frame(x)
    x$ENSEMBL <- rownames(x)
    x$SYMBOL <- mapIds(org.Hs.eg.db, keys= x$ENSEMBL , keytype ="ENSEMBL", column = "SYMBOL", multiVals = "first" )
    rownames(x)<- x$SYMBOL
    x
})
#get pathway annotaiton
rld_sub_sig<-lapply(rld_sub_sig, function(x){
    x<-left_join(x, t2g)
    x <- x[!duplicated(x$SYMBOL),]
    rownames(x)<- x$SYMBOL
    x
})

#plot
anno <- colData(rld_b)
heat<-list()
for (i in names(rld_sub_sig)){
    anno_row <- rld_sub_sig[[i]][,c("Pathway"), drop=F]

    heat[[i]]<- pheatmap(rld_sub_sig[[i]][,rownames(anno)], scale="row", show_colnames=F,#labels_row=annorow,
                       annotation_col=as.data.frame(annovst),fontsize_row=5,#annotation_row=anno_row,
                       annotation_colors=anno_colors, 
                       clustering_distance_rows="correlation",
                       show_rownames=T) 
    pdf(file.path(PostDE.dir,i,"Heat_sig_DNAdamageShortList.pdf"))
    print(heat[[i]])
    dev.off() 

    heat[[i]]<- pheatmap(rld_sub_sig[[i]][, c(rownames(anno[anno$group==as.character(contrasts[1,i]),]), rownames(anno[anno$group==as.character(contrasts[2,i]),]))],
                        scale="row", show_colnames=F,
                       annotation_col=as.data.frame(annovst),fontsize_row=5,#annotation_row=anno_row,
                       annotation_colors=anno_colors, clustering_distance_rows="correlation",
                       show_rownames=T) 
    pdf(file.path(PostDE.dir,i,"Heat_sig_DNAdamageShortList_subset.pdf"))
    print(heat[[i]])
    dev.off() 
}

#As barplot
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]
for (i in names(DEG_results_list_sub_sig_onlyFDR)){
    DEG_results_list_sub_sig_onlyFDR[[i]]$Symbol <- DEG_results_list_sub_sig_onlyFDR[[i]]$symbol
    DEG_results_list_sub_sig_onlyFDR[[i]]<-     DEG_results_list_sub_sig_onlyFDR[[i]][order(DEG_results_list_sub_sig_onlyFDR[[i]]$log2FoldChange, decreasing=TRUE),]
    DEG_results_list_sub_sig_onlyFDR[[i]]$Change <- ifelse(DEG_results_list_sub_sig_onlyFDR[[i]]$log2FoldChange >0, "blue", "red")
    pdf(file.path(PostDE.dir,i,"Bar_sig_DNAdamageShortList.pdf"), width=16)
    print(ggbarplot(DEG_results_list_sub_sig_onlyFDR[[i]], x="Symbol", y="log2FoldChange",fill="Change",palette=col)
    + rremove("xlab")+ rremove("legend")+rotate_x_text(angle = 45))
    dev.off()
}

#look at overlap
library(UpSetR)
t2g_split <- split(t2g, t2g$Pathway)
t2g_split <- lapply(t2g_split, function(x){
    x <- x$symbol
    x
})
pdf(file.path(PostDE.dir, "Overlap_DNAdamageShortList.pdf"))
upset(fromList(t2g_split), order.by = "freq")
dev.off()

