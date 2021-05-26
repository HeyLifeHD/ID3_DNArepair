#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data"
data.dir <- file.path("icgc/dkfzlsdf/analysis/C010/id3/odcf_featureCounts/")

library(data.table)

#load counts
files <- dir(path = data.dir, full.names = T,
             pattern = ".tsv")

count_list <- list()
for (cov_file in files){
  count_list[[cov_file]] <- fread(input =  cov_file)
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




#load data and annotation
#annotation
sample_names <- sapply(strsplit(files ,"//", fixed=TRUE),`[`, 2)
sample_names <- sapply(strsplit(sample_names ,".fpkm", fixed=TRUE),`[`, 1)

knockout_status <- sapply(strsplit(sample_names ,"-", fixed=TRUE),`[`, 1)

treatment_status <- sapply(strsplit(sample_names ,"-", fixed=TRUE),`[`, 2)
treatment_status<- sapply(strsplit(treatment_status ,"_", fixed=TRUE),`[`, 1)

genotype <- sapply(strsplit(sample_names ,"_", fixed=TRUE),`[`,4)

anno <- data.frame(sample_name = as.character(sample_names), knockout_status= as.character(knockout_status), tamoxifen_treatment=as.character(treatment_status), genotype=as.character(genotype ))
anno[grep("FL", anno$sample_name),"genotype"] <- "ID3_rescue"

anno[grep("irradiated", anno$tamoxifen_treatment),]$tamoxifen_treatment <- "untreated"

anno$hdac_treatment <- "untreated"
anno[grep("HDACi", anno$sample_name),]$hdac_treatment <- "treated"

anno$irradiation_treatment <- "untreated"
anno[grep("irradiated", anno$sample_name),]$irradiation_treatment <- "treated"

anno$replicate <- 0
anno[grep("R1", anno$sample_name),]$replicate <- 1
anno[grep("R2", anno$sample_name),]$replicate <- 2


anno$group_complete<- paste0(anno$genotype, "_Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment, "_Rep.", anno$replicate)
anno$group_uncomplete<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment, "_Rep.", anno$replicate)
anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment)

rownames(anno)<- anno$sample_name

saveRDS(anno ,file.path(base.dir,"data", "anno.rds"))

write.table(anno[,1:7], file.path(base.dir, "data", "anno.txt"), quote=FALSE, sep=";")
colnames(counts)<- as.character(anno$sample_name)
colnames(FPKM)<-  as.character(anno$sample_name)
colnames(TPM)<-  as.character(anno$sample_name)

saveRDS(counts, file.path(base.dir,"data", "counts.rds"))
saveRDS(FPKM, file.path(base.dir,"data", "FPKM.rds"))
saveRDS(TPM, file.path(base.dir, "data","TPM.rds"))
write.table(FPKM, file.path(base.dir,"data", "FPKM.txt"), quote=F)