#folder

base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190415_merged/all_data"
data.dir <- file.path("icgc/dkfzlsdf/analysis/C010/id3/190412_RNAseq_processing/")

library(data.table)

#load counts
counts <- fread(file.path(data.dir, "featureCounts/merged_gene_counts.txt"))

#subset counts
counts<- as.data.frame(counts)
rownames(counts)<- counts$gene_id

#read annotation table 
sample_anno <- read.table(file.path(base.dir, "sampleAnnotation.csv"), stringsAsFactors=FALSE, header=TRUE, sep=";")
#load data and annotation
#annotation

knockout_status <- sapply(strsplit(sample_anno$Patient.ID ,"_", fixed=TRUE),`[`, 3)
knockout_status[grep("FL", sample_anno$Patient.ID)]<- "ID3_rescue"

treatment_status <- sapply(strsplit(sample_anno$Sample.Type ,"-", fixed=TRUE),`[`, 2)
treatment_status<- sapply(strsplit(treatment_status ,"_", fixed=TRUE),`[`, 1)

anno <- data.frame(sample_name = sample_anno$ASID, genotype= as.character(knockout_status))
anno$tamoxifen_treatment <- "untreated"
anno[grep("+T", sample_anno$identifiere.name, fixed=TRUE),]$tamoxifen_treatment <- "treated"
anno[grep("+4OHT", sample_anno$identifiere.name, fixed=TRUE),]$tamoxifen_treatment <- "treated"

anno$hdac_treatment <- "untreated"
anno[grep("+HDi", sample_anno$identifiere.name, fixed=TRUE),]$hdac_treatment <- "treated"
anno$irradiation_treatment <- "untreated"
anno[grep("IRRADIATED", sample_anno$Sample.Type, fixed=TRUE),]$irradiation_treatment <- "treated"
anno$replicate <- 0
anno[grep("R1", sample_anno$Patient.ID),]$replicate <- 1
anno[grep("R2", sample_anno$Patient.ID),]$replicate <- 2

anno$group_complete<- paste0(anno$genotype, "_Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment, "_Rep.", anno$replicate)
anno$group_uncomplete<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment, "_Rep.", anno$replicate)
anno$group_uncomplete2<- paste0("Tam.",anno$tamoxifen_treatment ,"_Irr.", anno$irradiation_treatment, "_Hdac.", anno$hdac_treatment)

rownames(anno)<- anno$sample_name

dir.create(file.path(base.dir,"data"))
saveRDS(anno ,file.path(base.dir,"data", "anno.rds"))

write.table(anno[,1:7], file.path(base.dir, "data", "anno.txt"), quote=FALSE, sep=";")
colnames(counts)<- as.character(anno$sample_name)

colnames(counts) <- sapply(strsplit(colnames(counts) ,"-.", fixed=TRUE),`[`, 1)

counts_new <- counts[, rownames(anno)]
rownames(counts_new) <- counts$Geneid

saveRDS(counts_new, file.path(base.dir,"data", "counts.rds"))

