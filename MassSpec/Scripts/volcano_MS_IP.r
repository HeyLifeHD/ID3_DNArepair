#Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(knitr)
library(clusterProfiler)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)
library(readxl)


#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data.dir <- file.path(base.dir, "data")
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")

#read in excel sheets
Irr_vs_UT_15min <- as.data.frame(read_excel(file.path(base.dir, "IP_MS", "volcano", "200130_15min_vs_UT_Table.xlsx" ), sheet=1))
Irr_vs_UT_1h <- as.data.frame(read_excel(file.path(base.dir, "IP_MS", "volcano", "200130_1h_vs_UT_Table.xlsx" ), sheet=1))


DEG_results_list_plot<- list(Irr_vs_UT_15min, Irr_vs_UT_1h)
names(DEG_results_list_plot)<- c("Irr_vs_UT_15min", "Irr_vs_UT_1h")

#Plot results in Volcano and MA plot {.tabset}
DEG_results_list_plot <- lapply(DEG_results_list_plot, function(x){
  x <- as.data.frame(x)
  x$col <- "black"
  x$col <- ifelse(test = x$"-LOG(P-value)" > 1.3 & x$Difference>0.5, yes = "red", 
                  no = ifelse(test = x$"-LOG(P-value)" >1.3  & x$Difference<(-0.5), yes = "blue", "black"))
    x$plot <- NA
    #x$plot <- ifelse(x$"Gene names"  %in% c("RAD50", "c", "NBN", "PRMT1;HRMT1L2", "CDK2;CDK3") ,x$"Gene names", NA)
    x$plot <- ifelse(x$"Gene names"  %in% c("PRMT1;HRMT1L2", "CDK2;CDK3") ,x$"Gene names", NA)
    x$plot <- gsub("CDK2;CDK3", "CDK2",x$plot  )
    x$plot <- gsub("PRMT1;HRMT1L2", "PRMT1",x$plot  )
  x$pvalue <-x$"-LOG(P-value)"
  #x$padj <- (-(log10(x$padj)))
  x
}) 


#Volcano Plot
library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- c("grey", col[3], col[1])

Volcanos<- list() 
for (i in names(DEG_results_list_plot)){
  Volcanos[[i]] <- ggscatter(DEG_results_list_plot[[i]], x ="Difference", y = "pvalue",
          xlab="Normalized fold change",#xscale="log10" ,
          ylab="- log10(p-value)", #xlim=c(-0.5,0.5),
          color = "col" ,shape = 16, size = 1.5, palette=col, # Points color, shape and size
          add.params = list(color = "grey", fill = "lightgray"), legend.title="Significance",# main=i# Customize reg. line
          label="plot", repel = T, font.label = c("black","bold"), 
) + rremove("legend") + geom_vline(xintercept=c(-0.5,0.5), linetype = 2) +geom_hline(yintercept=1.3, linetype = 2)
  pdf(file.path(base.dir, "IP_MS", "volcano",paste0(i,"volcano.pdf")))
  print(Volcanos[[i]])
  dev.off()
  pdf(file.path(base.dir, "IP_MS", "volcano",paste0(i,"volcano4.pdf")), height=3.5, width=3.5)
  print(Volcanos[[i]])
  dev.off()
}
