#libraries
library(openxlsx)
#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/"
data.dir <- file.path(base.dir, "data")
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
dir.create(base.dir)
base_results.dir <- file.path(base.dir, "results")
dir.create(base_results.dir)
results.dir<- file.path(base_results.dir , "tables")
dir.create(results.dir)
PreDE.dir <- file.path(base_results.dir,"PreDE")
dir.create(PreDE.dir)
PostDE.dir <- file.path(base_results.dir,"PostDE")
dir.create(PostDE.dir)

#load data
MS_data <- read.xlsx(file.path(data.dir, "MS_ID3_IP.xlsx"))
MS_data_short <- MS_data[, c("Gene.names", "Log2_15min/UT", "Log2_1h/UT")]
MS_data_short$"Log2_15min/UT" <- as.numeric(MS_data_short$"Log2_15min/UT")
MS_data_short$"Log2_1h/UT" <- as.numeric(MS_data_short$"Log2_1h/UT")
#remove double gene names
new_df <- list()
for (i in 1:nrow(MS_data_short)){
    temp <- MS_data_short[i,]
    genes <- unlist(strsplit(temp$Gene.names, ";", fixed=TRUE))
    new_df[[i]] <- temp[rep(1, each = length(genes)), ]
    new_df[[i]]$Gene.names <- genes
}
new_df <- do.call("rbind",new_df)
saveRDS(new_df,file.path(data.dir, "MS_ID3_IP.rds") )

new_df_sub <- new_df[which(abs(new_df$"Log2_15min/UT")>1 & abs(new_df$"Log2_1h/UT")>1),]
dim(new_df_sub)
as.factor(new_df_sub[new_df_sub$"Log2_15min/UT">0,]$Gene.names)

#load expression results
DEG_results_list<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))
DEG_results_list_sig <- lapply(DEG_results_list, function(x){
    x <- x[which(x$padj < 0.05 & abs(x$log2FoldChange)>1),]
    x
})
lapply(DEG_results_list_sig,function(x)dim(x))
DEG_results_list_onlyFDR <- lapply(DEG_results_list, function(x){
    x <- x[which(x$padj < 0.05) ,]
    x
})
lapply(DEG_results_list_onlyFDR,function(x)dim(x))

results_FDR <- DEG_results_list_onlyFDR$ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED

#look for overlap of ip stuffed and deg
temp <- results_FDR[results_FDR$symbol  %in% new_df_sub$Gene.names,]
dim(temp)

dir.create(base.dir, "IP_MS")
write.table(temp, file.path(base.dir, "IP_MS", "ID3_WT_Irr_FDR_overlap_IP.txt"), row.names=FALSE, sep="\t")



#Plot Metascape
library(readxl)
library(ggpubr)
library(RColorBrewer)
#function
plotMeta <- function(excel_path, n=10, output_path, color="#CD534CFF", height=3.5, width=3.5){
dataset <- read_excel(excel_path, sheet=2)
datasetsub <- dataset[grep("Summary", dataset$GroupID), ]
datasetsub <- head(datasetsub, n)
datasetsub$LogP<- abs(datasetsub$LogP)
datasetsub$Description<- rev(datasetsub$Description)
datasetsub$LogP<- rev(datasetsub$LogP)
datasetsub$Term<- rev(datasetsub$Term)
datasetsub$description <- as.factor(1:nrow(datasetsub))
 p <- ggplot(datasetsub, aes(x=description)) +
geom_bar(aes(x=description, y=LogP, fill = "-log10(p.value)"), stat = "identity", width=0.25)
 p<-p+ annotate("text",x=as.integer(datasetsub$description)+0.35, y=0.01, 
                label=datasetsub$Description, size=4, hjust = 0)
 p <- p + scale_y_continuous( expand = c(0, 0), name = "-log10(p-value)")
 p <- p + scale_fill_manual(values = c( color))
 p<-p+coord_flip() 
 p<-p+ ggpubr:::theme_pubr()
 p<-p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=15))
 p<-p+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +rremove(c("legend"))+rremove(c("ylab"))
 p$labels$fill <- ""
 pdf(output_path, height=height, width=width)
 print(p)
 dev.off()
}

plotMeta( file.path(base.dir, "IP_MS","IP_Up", "metascape_result.xlsx"), file.path(base.dir, "IP_MS", "IP_Up","metascape_vis_20.pdf"), n=20, width=7, height=7)
plotMeta( file.path(base.dir, "IP_MS","IP_Up", "metascape_result.xlsx"), file.path(base.dir, "IP_MS", "IP_Up","metascape_vis_10.pdf"), n=10, width=5, height=5)
plotMeta( file.path(base.dir, "IP_MS","IP_Up", "metascape_result.xlsx"), file.path(base.dir, "IP_MS","IP_Up", "metascape_vis_5.pdf"), n=5, width=5, height=5)

plotMeta(file.path(input.dir, "down", "metascape_result.xlsx"), file.path(input.dir, "down", "metascape_vis3.pdf"),color= "#67A9CF", n=20, width=7, height=7)


#plot one subterm
dataset <- read_excel(file.path(base.dir, "IP_MS","IP_Up", "metascape_result.xlsx"), sheet=2)
dataset <- as.data.frame(dataset)
datasetsub <- dataset[grep("6_Member", dataset$GroupID), ]
datasetsub <- head(datasetsub, 10)
datasetsub$LogP<- abs(datasetsub$LogP)
datasetsub$Description<- rev(datasetsub$Description)
datasetsub$LogP<- rev(datasetsub$LogP)
datasetsub$Term<- rev(datasetsub$Term)
datasetsub$description <- as.factor(1:nrow(datasetsub))
 p <- ggplot(datasetsub, aes(x=description)) +
geom_bar(aes(x=description, y=LogP, fill = "-log10(p.value)"), stat = "identity", width=0.25)
 p<-p+ annotate("text",x=as.integer(datasetsub$description)+0.35, y=0.01, 
                label=datasetsub$Description, size=4, hjust = 0)
 p <- p + scale_y_continuous( expand = c(0, 0), name = "-log10(p-value)")
 p <- p + scale_fill_manual(values = c( "#CD534CFF"))
 p<-p+coord_flip() 
 p<-p+ ggpubr:::theme_pubr()
 p<-p+ theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(size=15))
 p<-p+theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +rremove(c("legend"))+rremove(c("ylab"))
 p$labels$fill <- ""
 pdf(file.path(base.dir, "IP_MS","IP_Up", "metascape_Vis_DNAreplication.pdf"), height=5, width=5)
 print(p)
 dev.off()
