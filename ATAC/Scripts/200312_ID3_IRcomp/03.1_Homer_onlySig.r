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

#create bg
DE_list_complete_df <- as.data.frame(DE_list_complete[[1]])
DE_list_complete_df$peak_id <- paste(DE_list_complete_df$seqnames,DE_list_complete_df$start,
    DE_list_complete_df$end, DE_list_complete_df$strand, sep="_")
dir.create(file.path(analysis.dir,"homer_onlySig"))
write.table(data.frame(chr=DE_list_complete_df$seqnames, 
        start=DE_list_complete_df$start, end=DE_list_complete_df$end, 
        id=DE_list_complete_df$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir,"homer_onlySig", "all_peaks_df.txt"),
        ,sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

#export lists for open 
#annotate peaks
DE_list_complete_sign <- lapply(DE_list_complete_sign, function(x){
    x <- as.data.frame(x)
    x$peak_id <- paste(x$seqnames,x$start,
    x$end, x$strand, sep="_")
    x
})

dars_final_df<-list()
for(i in names(DE_list_complete_sign)){
    dars_final_df[[i]] <- DE_list_complete_sign[[i]]
    dars_final_df[[i]]<- dars_final_df[[i]][which(dars_final_df[[i]]$Fold >1),]
    dir.create(file.path(analysis.dir, "homer_onlySig",i, "open"), recursive=TRUE)
    #open
    write.table(data.frame(chr=dars_final_df[[i]]$seqnames, 
        start=dars_final_df[[i]]$start, end=dars_final_df[[i]]$end, 
        id=dars_final_df[[i]]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer_onlySig",i, "open", paste0("dars.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
    #closed
    dars_final_df[[i]] <- DE_list_complete_sign[[i]]
    dars_final_df[[i]]<- dars_final_df[[i]][which(dars_final_df[[i]]$Fold <(-1)),]
    dir.create(file.path(analysis.dir, "homer_onlySig",i, "closed"), recursive=TRUE)
    write.table(data.frame(chr=dars_final_df[[i]]$seqnames, 
        start=dars_final_df[[i]]$start, end=dars_final_df[[i]]$end, 
        id=dars_final_df[[i]]$peak_id, notUsed=NA, strand="+"),
        file.path(analysis.dir, "homer_onlySig",i, "closed", paste0("dars.bed")),sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)

}

#run in command line
cd /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/DiffBind/homer_onlySig
conda activate homer2

for file in `ls */*/dars.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path} -size given -preparsedDir ${path}/ -p 6
    echo ${path}
done

for file in `ls */*/dars.bed`
do
    echo ${file}
    path=`dirname ${file}`
    findMotifsGenome.pl ${file} hg19 ${path}_bg -size given -preparsedDir ${path}/ -p 6 \
    -bg all_peaks_df.txt
    echo ${path}
done
