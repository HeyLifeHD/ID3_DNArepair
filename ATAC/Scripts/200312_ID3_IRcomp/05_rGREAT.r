##### Joschka Hey 
##### 
#Libraries
library(dplyr)
library(ggpubr)
library(LSD)
library(rGREAT)
library(forcats)
library(ggpubr)


#functions

#function
#plotting function
bubblPlotEnr<- function(data, n, set){
enr <- data[order(data$Hyper_Adjp_BH, decreasing=FALSE),]
datasetsub<-enr[1:n,]
datasetsub$LogP<- abs(log10(datasetsub$Hyper_Adjp_BH))
datasetsub$order <- 1:nrow(datasetsub)
datasetsub$significant<- ifelse(datasetsub$LogP< -log10(0.05), "No", "Yes" )
datasetsub$set <- set
datasetsub<- mutate(datasetsub, name = fct_reorder(name, LogP))

label_func <- function(x){
    breaks <- x
    breaks[breaks>=300] <- ">=300"
    breaks
}
ggplot(data = datasetsub, aes(y=name, x=set))+coord_fixed()+
    geom_point(aes(size=LogP, fill=Hyper_Fold_Enrichment), pch=21)+
    scale_fill_gradient2( midpoint = 0, low="white", high="darkred", name = "Hyper fold enrichment")+
    #scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
    scale_size(name="Q-value\n(-log10)", labels = label_func) +
    #scale_y_discrete(limits=rev(levels(as.factor(datasetsub$name))))+
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
}

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

#create lists to compare
test <- list()
for(i in names(DE_list_complete_sign_fc)){
    temp <- DE_list_complete_sign_fc[[i]]
    open <- temp[temp$Fold >1, ]
    closed <- temp[temp$Fold < (-1), ]

    test[[i]]<- list(open, closed)
    names(test[[i]])<- c("open", "closed")
}

universe <- DE_list_complete[[1]]


#rgreat enrichment
#rules for rGreat
Dist<- 10000
rule<-"basalPlusExt"


for (i in names(test)){
job_open <- submitGreatJob(test[[i]]$open, bg=universe,
    species= "hg19", includeCuratedRegDoms = TRUE,
    rule= c("basalPlusExt"),adv_upstream= 5.0,adv_downstream= 1.0,
    adv_span= 1000.0,adv_twoDistance= 1000.0,
    adv_oneDistance= 1000.0, request_interval = 300, max_tries = 10, version = "default")
job_closed <- submitGreatJob(test[[i]]$closed, bg=universe,
    species= "hg19", includeCuratedRegDoms = TRUE,
    rule= c("basalPlusExt"),adv_upstream= 5.0,adv_downstream= 1.0,
    adv_span= 1000.0,adv_twoDistance= 1000.0,
    adv_oneDistance= 1000.0, request_interval = 300, max_tries = 10, version = "default")


#hypo DARs enrichments
enrichments <- list()
dir.create(file.path(analysis.dir,"rGreat_SignFC" ,i, "openDARs"), recursive=TRUE)
for(j in availableCategories(job_open)){
    enrichments[[j]] <- getEnrichmentTables(job_open, category=j)
    print(j)
    for (k in names(enrichments[[j]])){
        table <- enrichments[[j]][[k]]
        #p <- plotGreat(tb_pathways_MSigDB)
        #pdf(file.path(analysis.dir, "rGreat", "pCAFvsNMF",paste0("allDMRs_",i,j, ".pdf"), height=6, width=6)
        #print(p)
        #dev.off()
    
        pdf(file.path(analysis.dir,"rGreat_SignFC" ,i, "openDARs",paste0("openDARs_bubble","_",j,"_",k, ".pdf")), height=7, width=14)
        print(bubblPlotEnr(table, 30, "open DARs"))
        dev.off()
        write.table(table, file.path(analysis.dir,"rGreat_SignFC" ,i, "openDARs",paste0("openDARs_","_",j,"_",k, ".txt")))
    }
}

#cloed DARs enrichments
enrichments <- list()
dir.create(file.path(analysis.dir,"rGreat_SignFC" ,i, "closedDARs"), recursive=TRUE)
for(j in availableCategories(job_closed)){
    enrichments[[j]] <- getEnrichmentTables(job_closed, category=j)
    print(j)
    for (k in names(enrichments[[j]])){
        table <- enrichments[[j]][[k]]
        #p <- plotGreat(tb_pathways_MSigDB)
        #pdf(file.path(analysis.dir, "rGreat", "pCAFvsNMF",paste0("allDMRs_",i,j, ".pdf"), height=6, width=6)
        #print(p)
        #dev.off()
    
        pdf(file.path(analysis.dir,"rGreat_SignFC" ,i, "closedDARs",paste0("closedDARs_bubble","_",j,"_",k, ".pdf")), height=7, width=14)
        print(bubblPlotEnr(table, 30, "closed DARs"))
        dev.off()
        write.table(table, file.path(analysis.dir,"rGreat_SignFC" ,i, "closedDARs",paste0("closedDARs","_",j,"_",k, ".txt")))
    }
}
}
