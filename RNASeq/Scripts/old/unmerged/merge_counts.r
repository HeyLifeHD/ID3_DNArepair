#folder
base.dir<- "/home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged"
data.dir <- file.path(base.dir, "data")

#load data and annotation
#annotation
anno <- read.table(file.path(data.dir, "annotation_unmerged.csv"), sep=";", stringsAsFactors=FALSE)
anno$sample_name <- paste0(anno$V1, "_", anno$V2)
colnames(anno)<- c("PID", "Treatment", "ID", "sample_name")
#modify anno
anno$celltype <- sapply(strsplit(anno$sample_name ,"_", fixed=TRUE),`[`, 2)
anno$genotype <- sapply(strsplit(anno$sample_name ,"_", fixed=TRUE),`[`, 3)
anno$replicate <- sapply(strsplit(anno$sample_name ,"_", fixed=TRUE),`[`, 4)


#data
counts <- read.table(file.path(data.dir, "merged_gene_counts.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(counts)<- counts$Geneid
counts <-  counts[,-c(1:2)]
colnames(counts)<- sapply(strsplit(colnames(counts) ,"Ali", fixed=TRUE),`[`, 1)
colnames(counts) <- gsub(".", "-", colnames(counts), fixed=TRUE)
#sum up unmerged samples
count_list <- list(NULL)
for (i in unique(anno$sample_name)){
    temp <- counts[,anno[anno$sample_name == i,]$ID ]
    if(class(temp)=="integer"){
        count_list[[i]]<- data.frame(temp)
       
    } else {
         count_list[[i]]<- data.frame(rowSums(temp))
    }
}
count_list <- count_list[-1]
count_list <- do.call( "cbind", count_list)
colnames(count_list)<- unique(anno$sample_name)

#create updated anno
anno_sub <- anno[!duplicated(anno$sample_name),]
anno_sub$replicate <- c(0,0,1,1,2,2,0,0,1,1,2,2,0,0,1,1,2,2,1,2,1,2,1,1,2,2,0,0,0,0,1,2,1,2,1,1,2,2)
anno_sub$dox_treatment <- sapply(strsplit(anno_sub$Treatment ,"-", fixed=TRUE),`[`, 2)
anno_sub$hdac_treatment <- "untreated"
anno_sub[grep("HDACi", anno_sub$sample_name),]$hdac_treatment <-"treated"
saveRDS(anno_sub, file.path(data.dir, "anno.rds"))

colnames(count_list)<- anno$sample_name
saveRDS(count_list, file.path(data.dir, "gene_counts.rds"))