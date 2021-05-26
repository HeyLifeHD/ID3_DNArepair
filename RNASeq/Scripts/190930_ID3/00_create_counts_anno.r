#folder
base.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/"
dir.create(base.dir)
data.dir <- "icgc/dkfzlsdf/project/OE0219/id3/sequencing/rna_sequencing/view-by-pid/"

#libraries
library(data.table)

#load counts
files <-list.files(path=data.dir,pattern=".featureCounts.tsv$", recursive = TRUE)
files <- paste0(data.dir, files)
count_list<-list()
for (cov_file in files){
  count_list[[cov_file]] <- fread(input =  cov_file)
}

#get real ensembl id
for(i in names(count_list)){
    count_list[[i]]$gene_id<- sapply(strsplit(count_list[[i]]$gene_id, ".", fixed=TRUE), "[",1)
    }
#subset counts
counts <- lapply(count_list , function(x){
    x <- x[, "num_reads"]
    x
})
counts <- do.call("cbind", counts)
rownames(counts)<- count_list[[1]]$gene_id
#subset fpkm
FPKM <- lapply(count_list , function(x){
    x <- x[, "FPKM"]
    x
})
FPKM <- do.call("cbind", FPKM)
rownames(FPKM)<- count_list[[1]]$gene_id
#subset tpm
TPM <- lapply(count_list , function(x){
    x <- x[, "TPM"]
    x
})
TPM <- do.call("cbind", TPM)
rownames(TPM)<- count_list[[1]]$gene_id


#load data and annotation#
sample <-read.table("/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/data/190930_ali_sample_sheet.tsv", sep="\t",header=TRUE, stringsAsFactors=FALSE)[,1:18]
sample$genotype <- ifelse(sample$genotype=="RESCUE", "ID3_rescue",  sapply(strsplit(sample$Patient.ID,"_", fixed=TRUE),`[`, 3))
#annotation
sample$match <- paste0(tolower(sample$Sample.Type),"_", sample$Patient.ID, ".fpkm_tpm.featureCounts.tsv.num_reads")
colnames(counts) <- sapply(strsplit(colnames(counts),"/", fixed=TRUE),`[`, 14)
sample$match %in% colnames(counts)
counts <- as.data.frame(counts)
counts <- counts[,sample$match]
colnames(counts)==sample$match

#subsetting and renaming
sample_anno <- sample[,c("genotype", "irradiation_treatment", "tamoxifen_treatment", "hdac_treatment.","Replicate" ,"Plattform", "id3_irradiation_project","id3_hdac_project", "arid1a_project")]
rownames(sample_anno)<- paste0(sample_anno$genotype, "_Irr.", sample_anno$irradiation_treatment, "_Tam.",sample_anno$tamoxifen_treatment, "_Hdac.",sample_anno$hdac_treatment, "_Rep.",sample_anno$Replicate)
colnames(counts)<- rownames(sample_anno)
sample_anno$group <- paste0(sample_anno$genotype, "_Irr.", sample_anno$irradiation_treatment, "_Tam.",sample_anno$tamoxifen_treatment, "_Hdac.",sample_anno$hdac_treatment)
rownames(counts) <- count_list[[1]]$gene_id
#subset projects
sample_anno_irr <- sample_anno[sample_anno$id3_irradiation_project==TRUE,]
counts_irr <- counts[,rownames(sample_anno_irr)]

sample_anno_arid <- sample_anno[sample_anno$arid1a_project==TRUE,]
counts_arid <- counts[,rownames(sample_anno_arid)]

sample_anno_hdac <- sample_anno[sample_anno$hdac_treatment==TRUE,]
counts_hdac <- counts[,rownames(sample_anno_hdac)] 

saveRDS(sample_anno ,file.path(base.dir,"data", "sample_anno.rds"))
write.table(sample_anno,file.path(base.dir,"data", "sample_anno.txt"), sep="\t",row.names=TRUE, col.names=TRUE, quote=FALSE)
saveRDS(counts ,file.path(base.dir,"data", "counts.rds"))

#save subprojects
counts_list <- list(counts_irr, counts_arid, counts_hdac)
names(counts_list)<- c("id3_irradiation_project","arid1a_project", "id3_hdac_project")
sample_anno_list <- list(sample_anno_irr, sample_anno_arid, sample_anno_hdac)
names(sample_anno_list)<- c("id3_irradiation_project","arid1a_project", "id3_hdac_project")
saveRDS(counts_list ,file.path(base.dir,"data", "counts_list.rds"))
saveRDS(sample_anno_list ,file.path(base.dir,"data", "sample_anno_list.rds"))


#new sample_anno
sample_anno <-read.delim(file.path(base.dir,"data","sample_anno_updated.tsv"))
rownames(sample_anno)<- sample_anno$SampleID
sample_anno$date_sequenced <- NA
sample_anno[sample_anno$first.sequencing ==TRUE,]$date_sequenced<-1
sample_anno[sample_anno$second.sequencing ==TRUE,]$date_sequenced<-2
sample_anno[sample_anno$third.sequencing ==TRUE,]$date_sequenced<-3
saveRDS(sample_anno ,file.path(base.dir,"data", "sample_anno.rds"))
write.table(sample_anno,file.path(base.dir,"data", "sample_anno.txt"), sep="\t",row.names=TRUE, col.names=TRUE, quote=FALSE)

#subset projects
sample_anno_irr <- sample_anno[sample_anno$id3_irradiation_project==TRUE,]

sample_anno_arid <- sample_anno[sample_anno$arid1a_project==TRUE,]

sample_anno_hdac <- sample_anno[sample_anno$id3_hdac_project==TRUE,]

sample_anno_NoTreat <- sample_anno[which(sample_anno$irradiation_treatment=="UNTREATED" & sample_anno$tamoxifen_treatment=="UNTREATED" &sample_anno$hdac_treatment=="UNTREATED"),]

sample_anno_Irradiated <- sample_anno[which(sample_anno$irradiation_treatment=="TREATED" & sample_anno$tamoxifen_treatment=="UNTREATED" &sample_anno$hdac_treatment=="UNTREATED"),]

sample_anno_list <- list(sample_anno_irr, sample_anno_arid, sample_anno_hdac,sample_anno_NoTreat,sample_anno_Irradiated)
names(sample_anno_list)<- c("id3_irradiation_project","arid1a_project", "id3_hdac_project", "NoTreatment", "OnlyIrradiated")
saveRDS(sample_anno_list ,file.path(base.dir,"data", "sample_anno_list.rds"))

#subset projects
counts<- readRDS(file.path(base.dir,"data", "counts.rds"))
counts_irr <- counts[,rownames(sample_anno_irr)]
counts_arid <- counts[,rownames(sample_anno_arid)]
counts_hdac <- counts[,rownames(sample_anno_hdac)]
counts_noTreat <- counts[,rownames(sample_anno_NoTreat)]
counts_Irradiated <- counts[,rownames(sample_anno_Irradiated)]

counts_list <- list(counts_irr, counts_arid, counts_hdac,counts_noTreat,counts_Irradiated)
names(counts_list)<- c("id3_irradiation_project","arid1a_project", "id3_hdac_project", "NoTreatment", "OnlyIrradiated")
saveRDS(counts_list ,file.path(base.dir,"data", "counts_list.rds"))