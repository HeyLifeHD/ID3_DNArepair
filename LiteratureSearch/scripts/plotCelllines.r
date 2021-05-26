#plot cell lines of interest
#directories
base.dir<- "c010-datasets/Internal/2018-Ali/LiteratureSearch"
analysis.dir <- file.path(base.dir, "analysis")
data.dir <- file.path(base.dir, "data")

#libraries 
library(ggpubr)

#load data
cl_id3<-read.table(file.path(data.dir, "mRNA expression (RNAseq)_ ID3.txt"), sep="\t",header=T)
#subset only ones with data
cl_id3 <- cl_id3[!is.na(cl_id3$ID3_expr),]
cl_id3$Cancer<- as.character(cl_id3$Cancer)
dim(cl_id3)
cl_id3$Cell_line_name <- sapply(strsplit(cl_id3$Cancer ,"_", fixed=TRUE),`[`, 1)
temp <- strsplit(cl_id3$Cancer ,"_", fixed=TRUE)
cl_id3$Cancer_type <- unlist(lapply(temp,FUN=function(x){paste(x[2:length(x)],collapse="_")}))

#subset cancer types of interest
cancer_oi <- c("BREAST", "PROSTATE", "PANCREAS")
cl_id3_sub<- cl_id3[cl_id3$Cancer_type %in% cancer_oi,]
dim(cl_id3_sub)

#select cell-lines for highlighting
cl_id3_sub$Cell_line_name[grep("BT", cl_id3_sub$Cell_line_name)]
cellline_oi <- c("MCF7", "MDAMB231", "MDAMB468","BT483","DU145", "LNCAPCLONEFGC", "MIAPACA2", "PSN1")
cl_id3_sub$label<-NA
cl_id3_sub$label<-ifelse(cl_id3_sub$Cell_line_name %in% cellline_oi, cl_id3_sub$Cell_line_name , NA)

col <- RColorBrewer::brewer.pal(n = length(unique(cl_id3_sub$Cancer_type)), name = 'Dark2')

pdf(file.path(analysis.dir, "boxplot_cancer_celllines_ofInterest2.pdf"), height=5, width=5)
ggboxplot(cl_id3_sub, x="Cancer_type", y="ID3_expr", label="label", repel=TRUE,fill="Cancer_type", col="grey", add="jitter",ylab="normalized ID3 expression", palette=col)+rremove("xlab")+rremove("legend")
dev.off()

pdf(file.path(analysis.dir, "dotplot_cancer_celllines_ofInterest2.pdf"), height=5, width=5)
ggdotplot(cl_id3_sub, x="Cancer_type", y="ID3_expr", label="label", repel=TRUE,fill="Cancer_type", col="grey",ylab="normalized ID3 expression", palette=col)+rremove("xlab")+rremove("legend")
dev.off()

cl_id3_sub[cl_id3_sub$Cell_line_name %in% c("MCF7", "MDAMB231"),]
