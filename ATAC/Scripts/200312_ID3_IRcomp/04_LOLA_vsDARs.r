##### Joschka Hey 
##### 
#Libraries
library(DSS)
library(bsseq)
library(doParallel)
library(ChIPseeker)
library(foreach)
library(SummarizedExperiment)
library(doMC)
library(rtracklayer)
library(HDF5Array)
library(pheatmap)
library(randomcoloR)
library(MeDeCom)
library(LOLA)
library(ggpubr)

#directories
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")

#load data
#create common object
DE_list_complete<- readRDS(file.path(analysis.dir, "de_list_complete.rds"))
#subset significant
DE_list_complete_sign <- lapply(DE_list_complete, function(x){
    x <- x[x$FDR<0.05,]
    x
})
lapply(DE_list_complete_sign, function(x)length(x))
#subset fold change
DE_list_complete_sign_fc <- lapply(DE_list_complete_sign, function(x){
    x <- x[abs(x$Fold)>1,]
    x
})
lapply(DE_list_complete_sign_fc, function(x)length(x))

#load LOLA
datasets.dir <- "c010-datasets/Internal/COPD/enrichment_databases/"
region_list <- readRDS(file.path(datasets.dir, "region_list_LOLA.rds"))
regionDB_genomicRegions <- region_list[[1]]
regionDB_homer <- region_list[[2]]
regionDB_chipSeq <- region_list[[3]]
regionDB_msigdb <- region_list[[4]]

#Run Enrichment
results_Core <- list()
results_genomicRegions <- list()
results_homer <- list()
results_chipSeq <- list()
results_msigdb <- list() 

for (i in names(DE_list_complete_sign_fc)){
    open <- DE_list_complete_sign_fc[[i]][which( DE_list_complete_sign_fc[[i]]$Fold>1),]
    closed <- DE_list_complete_sign_fc[[i]][which( DE_list_complete_sign_fc[[i]]$Fold<(-1)),]
    userSets<- list(open, closed)
    names(userSets)<- c("open","closed")
    dir.create(file.path(analysis.dir, "LOLA_vsDARs",i))
    #set  Universe
    userUnisverse <-DE_list_complete_sign_fc[[i]]
    #Run analysis
    #results_Core[[i]]= runLOLA(userSets, userUniverse, regionDB_Core, cores=4)
    results_genomicRegions[[i]]= runLOLA(userSets, userUnisverse, regionDB_genomicRegions, cores=3)
    results_homer[[i]]= runLOLA(userSets, userUnisverse, regionDB_homer, cores=3)
    results_chipSeq[[i]]= runLOLA(userSets, userUnisverse, regionDB_chipSeq, cores=3)
    results_msigdb[[i]]= runLOLA(userSets, userUnisverse, regionDB_msigdb, cores=3)
    print(i)
}

dir.create(file.path(analysis.dir,"LOLA_resultTables_vsDARs"))
saveRDS(results_genomicRegions, file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_genomicRegions.rds"))
saveRDS(results_homer, file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_homer.rds"))
saveRDS(results_chipSeq, file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_chipSeq.rds"))
saveRDS(results_msigdb, file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_msigdb.rds"))

results_genomicRegions<- readRDS(file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_genomicRegions.rds"))
results_homer<- readRDS(file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_homer.rds"))
results_chipSeq <- readRDS(file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_chipSeq.rds"))
results_msigdb <- readRDS(file.path(analysis.dir, "LOLA_resultTables_vsDARs", "results_msigdb.rds"))

#Plot genomic regions in bubble plot
#function 
label_func <- function(x){
    breaks <- x
    breaks[breaks>=200] <- ">=200"
    breaks
}
#subset multi-cell
results <- lapply(results_genomicRegions, function(x) x[grep("multi-cell", x$filename),])
#plotting
g<- list(NULL)

library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_vsDARs",i,"genomicRegions"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,"_"),`[`, 2)
combined_data$filename<-sapply(strsplit(combined_data$filename,"fe-"),`[`, 2)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_vsDARs",i,"genomicRegions", paste0("EnrichLOLA_","genomicRegions.pdf")))
print(g[i])
dev.off()
}



#Plot homer results in bubble plot
#subset multi-cell
results <- results_homer
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="hypo",]$filename, 20), head(x[x$userSet=="hyper",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_vsDARs",i,"homer"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_vsDARs",i,"homer", paste0("EnrichLOLA", "_homer.pdf")))
print(g[i])
dev.off()
}



#Plot regionDB_chipSeq results in bubble plot
#subset multi-cell
results <- results_chipSeq
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="open",]$filename, 20), head(x[x$userSet=="closed",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_vsDARs",i,"ChipSeq"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_vsDARs",i,"ChipSeq", paste0("EnrichLOLA_ChipSeq.pdf")), width=20)
print(g[i])
dev.off()
}




#For MsigDB
results <- results_msigdb
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="open",]$filename, 20), head(x[x$userSet=="closed",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_vsDARs",i,"MSigDB"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_vsDARs",i,"MSigDB", paste0("EnrichLOLA_MSigDB.pdf")))
print(g[i])
dev.off()
}


#For MsigDB only hallmarks
results <- results_msigdb
results <- lapply(results, function(x){
    x <- x[x$collection =="hallmarks",]
    temp <- unique(c(head(x[x$userSet=="open",]$filename, 20), head(x[x$userSet=="closed",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})



#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_vsDARs",i,"MSigDBHallmarks"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_vsDARs",i,"MSigDBHallmarks", paste0("EnrichLOLAHallmarks_MSigDB.pdf")))
print(g[i])
dev.off()
}



#same for only sign
for (i in names(DE_list_complete_sign)){
    open <- DE_list_complete_sign[[i]][which( DE_list_complete_sign[[i]]$Fold>1),]
    closed <- DE_list_complete_sign[[i]][which( DE_list_complete_sign[[i]]$Fold<(-1)),]
    userSets<- list(open, closed)
    names(userSets)<- c("open","closed")
    dir.create(file.path(analysis.dir, "LOLA_onlySig_vsDARs",i))
    #set  Universe
    userUnisverse <-DE_list_complete_sign[[i]]
    #Run analysis
    #results_Core[[i]]= runLOLA_onlySig(userSets, userUniverse, regionDB_Core, cores=4)
    results_genomicRegions[[i]]= runLOLA(userSets, userUnisverse, regionDB_genomicRegions, cores=3)
    results_homer[[i]]= runLOLA(userSets, userUnisverse, regionDB_homer, cores=3)
    results_chipSeq[[i]]= runLOLA(userSets, userUnisverse, regionDB_chipSeq, cores=3)
    results_msigdb[[i]]= runLOLA(userSets, userUnisverse, regionDB_msigdb, cores=3)
    print(i)
}

dir.create(file.path(analysis.dir,"LOLA_vsDARs_onlySig_resultTables"))
saveRDS(results_genomicRegions, file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_genomicRegions.rds"))
saveRDS(results_homer, file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_homer.rds"))
saveRDS(results_chipSeq, file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_chipSeq.rds"))
saveRDS(results_msigdb, file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_msigdb.rds"))

results_genomicRegions<- readRDS(file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_genomicRegions.rds"))
results_homer<- readRDS(file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_homer.rds"))
results_chipSeq <- readRDS(file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_chipSeq.rds"))
results_msigdb <- readRDS(file.path(analysis.dir, "LOLA_vsDARs_onlySig_resultTables", "results_msigdb.rds"))

#Plot genomic regions in bubble plot
#function 
label_func <- function(x){
    breaks <- x
    breaks[breaks>=200] <- ">=200"
    breaks
}
#subset multi-cell
results <- lapply(results_genomicRegions, function(x) x[grep("multi-cell", x$filename),])
#plotting
g<- list(NULL)

library(RColorBrewer)
col <- brewer.pal(3,"RdBu")
col <- col[c(1,3)]

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"genomicRegions"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LOLA_onlySigCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,"_"),`[`, 2)
combined_data$filename<-sapply(strsplit(combined_data$filename,"fe-"),`[`, 2)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"genomicRegions", paste0("EnrichLOLA_onlySig_","genomicRegions.pdf")))
print(g[i])
dev.off()
}



#Plot homer results in bubble plot
#subset multi-cell
results <- results_homer
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="open",]$filename, 20), head(x[x$userSet=="closed",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"homer"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LOLA_onlySigCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"homer", paste0("EnrichLOLA_onlySig", "_homer.pdf")))
print(g[i])
dev.off()
}



#Plot regionDB_chipSeq results in bubble plot
#subset multi-cell
results <- results_chipSeq
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="open",]$filename, 20), head(x[x$userSet=="closed",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"ChipSeq"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LOLA_onlySigCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"ChipSeq", paste0("EnrichLOLA_onlySig_ChipSeq.pdf")), width=20)
print(g[i])
dev.off()
}




#For MsigDB
results <- results_msigdb
results <- lapply(results, function(x){
    temp <- unique(c(head(x[x$userSet=="open",]$filename, 20), head(x[x$userSet=="closed",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})

#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"MSigDB"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LOLA_onlySigCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"MSigDB", paste0("EnrichLOLA_onlySig_MSigDB.pdf")))
print(g[i])
dev.off()
}


#For MsigDB only hallmarks
results <- results_msigdb
results <- lapply(results, function(x){
    x <- x[x$collection =="hallmarks",]
    temp <- unique(c(head(x[x$userSet=="open",]$filename, 20), head(x[x$userSet=="closed",]$filename, 20)))
    x <- x[x$filename %in% temp,]
    x
})



#plotting
g<- list(NULL)

for(i in names(results)){
    dir.create(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"MSigDBHallmarks"), recursive=TRUE)

  #data preparation
combined_data <- results[[i]][,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename", "qValue")]#[userSet=="closed",]
#combined_data$significant<- ifelse(locResults.LOLA_onlySigCore$qValue> 0.05, "No", "Yes" )
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
#combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)

g[[i]] <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
    geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
    scale_fill_gradient2( midpoint = 1, low="darkblue", high="darkred", name = "Odds Ratio")+
    scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="P-value\n(-log10)", labels = label_func) +
    scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
    theme(text =element_text(size=14, color="black", family = "sans"),
          axis.ticks = element_blank(), axis.line = element_blank(), 
          axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
          axis.text.y=element_text(size=12, family="sans", color="black"))+
    scale_x_discrete(name=NULL)+
    theme(legend.text=element_text(size=12, family="sans"), 
          legend.title=element_text(size=12, family= "sans"),
          legend.background = element_rect(fill="white", color="white"),
          panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
          legend.key = element_rect(fill="white"))+rremove("ylab")

pdf(file.path(analysis.dir,"LOLA_onlySig_vsDARs",i,"MSigDBHallmarks", paste0("EnrichLOLA_onlySigHallmarks_MSigDB.pdf")))
print(g[i])
dev.off()
}
