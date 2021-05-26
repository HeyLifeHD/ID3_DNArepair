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
regionDB = loadRegionDB(file.path(datasets.dir,"scratch/ns5bc/resources/regions/LOLACore/hg19/"))


#Run Enrichment
results_regionDB <- list() 
for (i in names(DE_list_complete_sign_fc)){
    open <- DE_list_complete_sign_fc[[i]][which( DE_list_complete_sign_fc[[i]]$Fold>1),]
    closed <- DE_list_complete_sign_fc[[i]][which( DE_list_complete_sign_fc[[i]]$Fold<(-1)),]
    userSets<- list(open, closed)
    names(userSets)<- c("open","closed")
    dir.create(file.path(analysis.dir, "LOLA",i))
    #set  Universe
    userUnisverse <-DE_list_complete[[i]]
    #Run analysis
    #results_Core[[i]]= runLOLA(userSets, userUniverse, regionDB_Core, cores=4)
    results_regionDB[[i]]= runLOLA(userSets, userUnisverse, regionDB, cores=3)

    print(i)
}

dir.create(file.path(analysis.dir,"LOLA_resultTables"))
saveRDS(results_regionDB, file.path(analysis.dir, "LOLA_resultTables", "results_regionDB.rds"))
results_regionDB<- readRDS(file.path(analysis.dir, "LOLA_resultTables", "results_regionDB.rds"))

#function
label_func <- function(x){
    breaks <- x
    breaks[breaks>=200] <- ">=200"
    breaks
}
#run plotting
for(i in names(results_regionDB)){
    result <- results_regionDB[[i]]
    for(j in unique(unique(result$collection))){
        result_sub<- result[result$collection == j,]
        temp <- unique(c(head(result_sub[result_sub$userSet=="open",]$filename, 20), head(result_sub[result_sub$userSet=="dmrcloseds_hyper",]$filename, 20)))
        result_sub2 <- result_sub[result_sub$filename %in% temp, ]


        #data preparation
        combined_data <- result_sub2[,c("userSet","dbSet", "pValueLog", "oddsRatio" ,"filename", "qValue")]#[userSet=="closed",]
        #combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
        combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
        #change infinite values
        #combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
        combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
        combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)
        #plot
        g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
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
        dir.create(file.path(analysis.dir,"LOLA",i), recursive=TRUE)
        pdf(file.path(analysis.dir,"LOLA",i, paste0("EnrichLOLA_", j, ".pdf")), width=15, height=20)
        print(g)
        dev.off()
    }
}






#Run Enrichmen for only Sig
results_regionDB <- list() 
for (i in names(DE_list_complete_sign)){
    open <- DE_list_complete_sign[[i]][which( DE_list_complete_sign[[i]]$Fold>1),]
    closed <- DE_list_complete_sign[[i]][which( DE_list_complete_sign[[i]]$Fold<(-1)),]
    userSets<- list(open, closed)
    names(userSets)<- c("open","closed")
    dir.create(file.path(analysis.dir, "LOLA_onlySig",i))
    #set  Universe
    userUnisverse <-DE_list_complete[[i]]
    #Run analysis
    #results_Core[[i]]= runLOLA(userSets, userUniverse, regionDB_Core, cores=4)
    results_regionDB[[i]]= runLOLA(userSets, userUnisverse, regionDB, cores=3)

    print(i)
}

dir.create(file.path(analysis.dir,"LOLA_onlySig_resultTables"))
saveRDS(results_regionDB, file.path(analysis.dir, "LOLA_onlySig_resultTables", "results_regionDB.rds"))
results_regionDB<- readRDS(file.path(analysis.dir, "LOLA_onlySig_resultTables", "results_regionDB.rds"))

#function
label_func <- function(x){
    breaks <- x
    breaks[breaks>=200] <- ">=200"
    breaks
}
#run plotting
for(i in names(results_regionDB)){
    result <- results_regionDB[[i]]
    for(j in unique(unique(result$collection))){
        result_sub<- result[result$collection == j,]
        temp <- unique(c(head(result_sub[result_sub$userSet=="open",]$filename, 20), head(result_sub[result_sub$userSet=="dmrcloseds_hyper",]$filename, 20)))
        result_sub2 <- result_sub[result_sub$filename %in% temp, ]


        #data preparation
        combined_data <- result_sub2[,c("userSet","dbSet", "pValueLog", "oddsRatio" ,"filename", "qValue")]#[userSet=="closed",]
        #combined_data$significant<- ifelse(locResults.LolaCore$qValue> 0.05, "No", "Yes" )
        combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
        #change infinite values
        #combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 10*(max(combined_data$pValueLog[!is.infinite(combined_data$pValueLog)],na.rm=TRUE))
        combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 200
        combined_data$filename<-sapply(strsplit(combined_data$filename,".bed", fixed=TRUE),`[`, 1)
        #plot
        g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
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
        dir.create(file.path(analysis.dir,"LOLA_onlySig",i), recursive=TRUE)
        pdf(file.path(analysis.dir,"LOLA_onlySig",i, paste0("EnrichLOLA_", j, ".pdf")), width=15, height=20)
        print(g)
        dev.off()
    }
}





