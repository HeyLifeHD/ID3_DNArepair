#GSEA with DNA damage list
#directories
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#packages
library(DESeq2)
library(clusterProfiler)

#load files
dds_genes <- readRDS(file = file.path(results.dir,"dds_group_genes.rds"))
vst_genes<- readRDS(file = file.path(results.dir,"vst_genes_batchRem_rep.rds"))

DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#load term2 gene list
t2g <- read.table(file.path(data.dir, "DRGTerm2Gene.csv"), sep=";", stringsAsFactors = F)
t2g$ensembl <-  mapIds(org.Hs.eg.db, keys=t2g$V2, keytype ="SYMBOL", column = "ENSEMBL", multiVals = "first" )
t2g <- t2g[!is.na(t2g$ensembl),]
t2g<- data.frame(TERM= t2g$V1, GENE=t2g$ensembl)

#Run GSEA
GSEA_DSS<- list()
for (i in names(DEG_results_list)){
  temp <- DEG_results_list[[i]]
  temp <- temp[!is.na(temp$padj),]
  geneList <-temp[,"log2FoldChange"]
  names(geneList) <- temp$ensembl
  geneList <- sort(geneList, decreasing=TRUE)
  
  GSEA_DSS[[i]] <- GSEA(geneList,TERM2GENE = t2g , pvalueCutoff=1)
}
#plot dotplot
plot2 <- list()
for (i in names(GSEA_DSS)){
  plot[[i]] <- dotplot(GSEA_DSS[[i]],showCategory = 20)
  pdf(file.path(PostDE.dir,i,"DNA_damage_GSEA.pdf"))
  print(plot[[i]])
  dev.off()
  plot2[[i]] <- dotplot(GSEA_DSS[[i]],showCategory = 20, color="pvalue")
  pdf(file.path(PostDE.dir,i,"DNA_damage_GSEA_pval.pdf"))
  print(plot2[[i]])
  dev.off()
}
#plot barcode plot
for (i in names(GSEA_DSS)){
temp <- GSEA_DSS[[i]]
  for (term in 1:length(temp$Description)) {
    pdf(file.path(PostDE.dir,i,paste0("SainiGSEA_Barcode_",temp$Description[term], ".pdf")),width = 7,height=3.5)
    print(gseaplot(temp, geneSetID =temp$ID[term]))
    dev.off()
  }
}
#plot barplot
GSEAPlot <- function(GSEARes){
  data <- as.data.frame(GSEARes)
  data <- data[order(abs(data$NES), decreasing = TRUE),]
  data <- data[, c("NES", "p.adjust","ID")]
  #data<- head(data, 10)
  
  data$p.adjust <- abs(log10(data$p.adjust ))
  data$factor <- as.factor(c(1:nrow(data)))
  data$factor<- rev(data$factor)
  p<-ggplot(data, aes(x=factor)) +
    geom_bar(aes(x=as.integer(as.factor(data$factor))-0.15, y=NES, fill = "Normalized Enrichment Score"), stat = "identity", width=0.25)+ 
    geom_bar(aes(x=as.integer(as.factor(data$factor))+0.15, y=p.adjust, fill = "P.adjusted"), stat = "identity", width=0.25)
  
  p<-p+ annotate("text",x=as.integer(as.factor(data$factor))+0.5, y=0.01, 
                 label=data$ID, size=4, hjust = 0)
  
  p <- p + scale_y_continuous(sec.axis = sec_axis(~.*2, name = "-log10(P.adjusted)"), expand = c(0, 0))
  p <- p + scale_fill_manual(values = c("darkolivegreen", "grey"))
  
  p<-p+coord_flip() 
  p<-p+ ggpubr:::theme_pubr()
  p<-p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=15))
  p<-p+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
  
  
  p <- p + labs(y = "Absolute enrichment score",
                x = NULL,
                colour = "Parameter")
  p$labels$fill <- ""
  print(p)
}
for (i in names(GSEA_DSS)){
  pdf(file.path(PostDE.dir,i,"DNA_damage_GSEA_barl.pdf"))
  GSEAPlot(GSEA_DSS[[i]])
  dev.off()
}

for (i in names(GSEA_DSS)){
  write.table(as.data.frame(GSEA_DSS[[i]]), file.path(PostDE.dir,i, "DNA_damage_GSEA.txt"),row.names = TRUE,col.names = TRUE)
}

GSEA_DSS[["WT_Tam.untreated_Irr.treated_Hdac.untreated_vs_ID3_Tam.untreated_Irr.treated_Hdac.untreated"]]
