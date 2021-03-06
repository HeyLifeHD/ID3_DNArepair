#directories
#atac
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp"
data.dir <- file.path("/icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/")
analysis.dir <- file.path(base.dir, "analysis", "DiffBind")
#RNA 
#Irradiation
base_Irr.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_only_irr_4rep"
data_Irr.dir <- file.path(base_Irr.dir, "data")
base_results_Irr.dir <- file.path(base_Irr.dir, "results")
results_Irr.dir<- file.path(base_results_Irr.dir , "tables")
PreDE_Irr.dir <- file.path(base_results_Irr.dir,"PreDE")
PostDE_Irr.dir <- file.path(base_results_Irr.dir,"PostDE")

#Unirradiated
base_UIrr.dir<- "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190930_odcf/id3_wo_irr_4rep"
data_UIrr.dir <- file.path(base_UIrr.dir, "data")
base_results_UIrr.dir <- file.path(base_UIrr.dir, "results")
results_UIrr.dir<- file.path(base_results_UIrr.dir , "tables")
PreDE_UIrr.dir <- file.path(base_results_UIrr.dir,"PreDE")
PostDE_UIrr.dir <- file.path(base_results_UIrr.dir,"PostDE")


#libraries
library(rtracklayer)
library(ChIPseeker)
library(ChIPpeakAnno)
library(reshape2)
library(ggpubr)
library(dplyr)
library(LSD)
library(peakSeason)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(EnsDb.Hsapiens.v75)
ENS<-EnsDb.Hsapiens.v75
seqlevelsStyle(ENS) <- "UCSC"


#load data
#ATAC
#DE_list <- readRDS(file.path(analysis.dir, "de_list.rds"))
#counts<- readRDS(file.path(analysis.dir, "countsperpeak.rds"))
#counts_raw<- readRDS(file.path(analysis.dir, "allPeaks_raw.rds"))
#sample_anno<- readRDS(file.path(analysis.dir, "sample_anno.rds"))
DAR_list_complete<- readRDS(file.path(analysis.dir, "de_list_complete.rds"))

#RNA
#Irr
DEG_results_list_Irr<- readRDS(file.path(PostDE_Irr.dir,  "DEG_results_group_list.rds"))

#UnIrr
DEG_results_list_UIrr<- readRDS(file.path(PostDE_UIrr.dir,  "DEG_results_group_list.rds"))

#Create Sample sheet
#get files 
bigwig_files <- list.files(data.dir, pattern="_tn5_center_73bp.bigwig$", full.names=TRUE, recursive=TRUE)
bigwig_files <- bigwig_files[grep("cache", bigwig_files, invert=TRUE)]

#GET NAMES
file_names <- sapply(strsplit(bigwig_files, "/", fixed=TRUE), "[", 11)
file_names <-gsub("ID3_KO", "ID3-ko", file_names)
sample_anno <- data.frame("SampleID"=file_names,
    "Factor"= sapply(strsplit(file_names, "_", fixed=TRUE), "[", 1),
    "Treatment"= sapply(strsplit(file_names, "_", fixed=TRUE), "[", 2),
    "Replicate"= paste0("rep",sapply(strsplit(file_names, "_", fixed=TRUE), "[", 3)),
    "PeakCaller"= rep("macs", length(file_names)),
    "bigwig"=bigwig_files)
sample_anno$Condition <- as.factor(paste0(sample_anno$Factor, "_", sample_anno$Treatment))
sample_anno$SampleID <- gsub("ID3-ko", "ID3.ko",sample_anno$SampleID)
rownames(sample_anno)<- sample_anno$SampleID


#create column data
col_data <- read_coldata(coldata=sample_anno, files_idx=6, sample_idx=1)
col <- c("darkorange4", "darkorange", "#7570B3",  "cyan3")

#Profile plot of all TSS
genes <- genes(ENS)
prom <- promoters(genes, upstream=1, downstream=1)
prom <- prom[seqnames(prom) %in%  c(paste0("chr", 1:21),"chrY", "chrX")]
prom_protein_coding <- prom[prom$gene_biotype=="protein_coding",]
dir.create(file.path(analysis.dir, "data","ProfilePlot"),recursive=TRUE)
export.bed(prom_protein_coding, file.path(analysis.dir, "data","ProfilePlot","Promoters.bed"))
#create matrix list
regions_matrix <- extract_matrix(coldata=col_data, 
    bed= file.path(analysis.dir, "data","ProfilePlot","Promoters.bed"), 
    genome=hg19, up=500, down=500,binSize=10, 
    nthreads=6)
#draw profile plot 
dir.create(file.path(analysis.dir, "ProfilePlots"))
pdf(file.path(analysis.dir, "ProfilePlots","all_TSS_allGroups_profileplot_500bp.pdf"), height=3,5, width=3.5)
plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=TRUE, condition="Condition", collapse_reps=TRUE, color=col)
dev.off()
#run statistical test
mean_list <- list()
for( i in rownames(sample_anno)){
mean_list[[i]] <- regions_matrix[[i]]
mean_list[[i]] <- mean(rowMeans(mean_list[[i]]))
}
shapiro.test( unlist(mean_list))
res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
test <- compare_means(value ~ group, data = res, method="t.test")
write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots","all_TSS_allGroups_profileplot_500bp_test.txt"), 
    row.names=FALSE, quote=FALSE, sep="\t")
mean_list <- list()
#test for center
for( i in rownames(sample_anno)){
mean_list[[i]] <- regions_matrix[[i]]
mean_list[[i]] <- mean(mean_list[[i]]$V50)
}
shapiro.test( unlist(mean_list))
res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
test <- compare_means(value ~ group, data = res, method="t.test")
write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots","all_TSS_allGroups_profileplot_500bp_test_center.txt"), 
    row.names=FALSE, quote=FALSE, sep="\t")
#create matrix list
regions_matrix <- extract_matrix(coldata=col_data, 
    bed= file.path(analysis.dir, "data","ProfilePlot","Promoters.bed"), 
    genome=hg19, up=1500, down=1500,binSize=10, 
    nthreads=6)
#draw profile plot 
dir.create(file.path(analysis.dir, "ProfilePlots"))
pdf(file.path(analysis.dir, "ProfilePlots","all_TSS_allGroups_profileplot_1500bp.pdf"), height=3,5, width=3.5)
plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=TRUE, condition="Condition", collapse_reps=TRUE, color=col)
dev.off()



#Profile Plot of all DNA damage player
t2g <- read.table(file.path( "/home/heyj/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data","data", "DRGTerm2Gene.csv"), sep=";", stringsAsFactors = F)
t2g$ensembl <-  mapIds(org.Hs.eg.db, keys=t2g$V2, keytype ="SYMBOL", column = "ENSEMBL", multiVals = "first" )
t2g <- t2g[!is.na(t2g$ensembl),]
t2g<- data.frame(TERM= t2g$V1, GENE=t2g$ensembl)
DNA_damage_promoter_all <- prom[prom$gene_id %in% t2g$GENE, ]
export.bed(DNA_damage_promoter_all, file.path(analysis.dir, "data","ProfilePlot","DNA_damage_promoter_all.bed"))
#500bp
#create matrix list
regions_matrix <- extract_matrix(coldata=col_data, 
    bed= file.path(analysis.dir, "data","ProfilePlot","DNA_damage_promoter_all.bed"), 
    genome=hg19, up=500, down=500,binSize=10, 
    nthreads=6)
#draw profile plot 
dir.create(file.path(analysis.dir, "ProfilePlots"))
pdf(file.path(analysis.dir, "ProfilePlots","DNA_damage_promoter_all_allGroups_profileplot_500bp.pdf"), height=3,5, width=3.5)
plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col)
dev.off()
#run statistical test
mean_list <- list()
for( i in rownames(sample_anno)){
mean_list[[i]] <- regions_matrix[[i]]
mean_list[[i]] <- mean(rowMeans(mean_list[[i]]))
}
shapiro.test( unlist(mean_list))
res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
test <- compare_means(value ~ group, data = res, method="t.test")
write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots","DNA_damage_promoter_all_allGroups_profileplot_500bp_test.txt"), 
    row.names=FALSE, quote=FALSE, sep="\t")
#test for center
for( i in rownames(sample_anno)){
mean_list[[i]] <- regions_matrix[[i]]
mean_list[[i]] <- mean(mean_list[[i]]$V50)
}
shapiro.test( unlist(mean_list))
res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
test <- compare_means(value ~ group, data = res, method="t.test")
write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots","DNA_damage_promoter_all_allGroups_profileplot_500bp_test_center.txt"), 
    row.names=FALSE, quote=FALSE, sep="\t")

#15000bp
#create matrix list
regions_matrix <- extract_matrix(coldata=col_data, 
    bed= file.path(analysis.dir, "data","ProfilePlot","DNA_damage_promoter_all.bed"), 
    genome=hg19, up=1500, down=1500,binSize=10, 
    nthreads=6)
#draw profile plot 
dir.create(file.path(analysis.dir, "ProfilePlots"))
pdf(file.path(analysis.dir, "ProfilePlots","DNA_damage_promoter_all_allGroups_profileplot_1500bp.pdf"), height=3,5, width=3.5)
plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col)
dev.off()


#loop over different dna damage player subtypes
t2g$TERM <- as.character(t2g$TERM )
for(i in unique(t2g$TERM)){
    t2g_sub <- t2g[t2g$TERM == i, ]
    region_oi <- prom[prom$gene_id %in% t2g_sub$GENE, ]
    name <- i 
    name <- gsub(" ", "", name)
    name <- gsub("(", "", name, fixed=TRUE)
    name <- gsub(")", "", name, fixed=TRUE)

    export.bed(region_oi, file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_",name,".bed")))

    regions_matrix <- extract_matrix(coldata=col_data, 
        bed= file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_",name,".bed")), 
        genome=hg19, up=500, down=500,binSize=10, 
        nthreads=6)
    pdf(file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_",name,"allGroups_profileplot_500bp.pdf")), height=3,5, width=3.5)
    print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
    dev.off()
    #test
    mean_list <- list()
    for( j in rownames(sample_anno)){
    mean_list[[j]] <- regions_matrix[[j]]
    mean_list[[j]] <- mean(rowMeans(mean_list[[j]]))
    }
    shapiro.test( unlist(mean_list))
    res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
    test <- compare_means(value ~ group, data = res, method="t.test")
    write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_",name,"allGroups_profileplot_500bp_test.txt")), 
    row.names=FALSE, quote=FALSE, sep="\t")
#test for center
    for( j in rownames(sample_anno)){
    mean_list[[j]] <- regions_matrix[[j]]
    mean_list[[j]] <- mean(mean_list[[j]]$V50)
    }
    shapiro.test( unlist(mean_list))
    res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
    test <- compare_means(value ~ group, data = res, method="t.test")
    write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_",name,"allGroups_profileplot_500bp_test_center.txt")), 
        row.names=FALSE, quote=FALSE, sep="\t")


    regions_matrix <- extract_matrix(coldata=col_data, 
        bed= file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_",name,".bed")), 
        genome=hg19, up=100, down=100,binSize=10, 
        nthreads=6)
    pdf(file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_",name,"allGroups_profileplot_100bp.pdf")), height=3,5, width=3.5)
    print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
    dev.off()

    regions_matrix <- extract_matrix(coldata=col_data, 
        bed= file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_",name,".bed")), 
        genome=hg19, up=1500, down=1500,binSize=10, 
        nthreads=6)
    pdf(file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_",name,"allGroups_profileplot_1500bp.pdf")), height=3,5, width=3.5)
    print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
    dev.off()
    print(name)
}

#E2F transcription factors
prom_e2f <-  prom[grep("E2F",prom$symbol),]
prom_e2f <-  prom_e2f[grep("UB",prom_e2f$symbol, invert=TRUE),]
prom_e2f <- prom_e2f[prom_e2f$gene_biotype=="protein_coding",]
export.bed(prom_e2f, file.path(analysis.dir, "data","ProfilePlot","prom_e2f_all.bed"))
#500bp
#create matrix list
regions_matrix <- extract_matrix(coldata=col_data, 
    bed= file.path(analysis.dir, "data","ProfilePlot","prom_e2f_all.bed"), 
    genome=hg19, up=500, down=500,binSize=10, 
    nthreads=6)
#draw profile plot 
dir.create(file.path(analysis.dir, "ProfilePlots"))
pdf(file.path(analysis.dir, "ProfilePlots","E2FPromoter_all_allGroups_profileplot_500bp.pdf"), height=3,5, width=3.5)
plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col)
dev.off()
#run statistical test
mean_list <- list()
for( i in rownames(sample_anno)){
mean_list[[i]] <- regions_matrix[[i]]
mean_list[[i]] <- mean(rowMeans(mean_list[[i]]))
}
shapiro.test( unlist(mean_list))
res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
test <- compare_means(value ~ group, data = res, method="t.test")
write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots","E2FPromoter_all_allGroups_profileplot_500bp_test.txt"), 
    row.names=FALSE, quote=FALSE, sep="\t")


#15000bp
#create matrix list
regions_matrix <- extract_matrix(coldata=col_data, 
    bed= file.path(analysis.dir, "data","ProfilePlot","prom_e2f_all.bed"), 
    genome=hg19, up=1500, down=1500, binSize=10, 
    nthreads=6)
#draw profile plot 
dir.create(file.path(analysis.dir, "ProfilePlots"))
pdf(file.path(analysis.dir, "ProfilePlots","E2FPromoter_all_allGroups_profileplot_1500bp.pdf"), height=3,5, width=3.5)
plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col)
dev.off()


#differentially expressed genes
DEG_results_list <- list(DEG_results_list_Irr$"ID3_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.TREATED_Tam.UNTREATED_Hdac.UNTREATED",
    DEG_results_list_UIrr$"ID3_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED_vs_WT_Irr.UNTREATED_Tam.UNTREATED_Hdac.UNTREATED")
names(DEG_results_list) <- c("ID3-ko_Irr_vs_wt_Irr", "ID3-ko_unIrr_vs_wt_unIrr")

#subset significant ATAC DARs
DEG_results_list_sign <- lapply(DEG_results_list, function(x){
    x <- x[which(x$padj<0.05),]
    x
})
lapply(DEG_results_list_sign, function(x)dim(x))
#subset fold change
DEG_results_list_sign_fc <- lapply(DEG_results_list_sign, function(x){
    x <- x[abs(x$log2FoldChange)>0.5,]
    x
})
lapply(DEG_results_list_sign_fc, function(x)dim(x))

for(i in names(DEG_results_list_sign_fc)){
    up <- DEG_results_list_sign_fc[[i]][DEG_results_list_sign_fc[[i]]$log2FoldChange>0, ]
    down <- DEG_results_list_sign_fc[[i]][DEG_results_list_sign_fc[[i]]$log2FoldChange<0, ]
    region_up <- prom[prom$gene_id %in% up$ensembl, ]
    region_down <- prom[prom$gene_id %in% down$ensembl, ]

    export.bed(region_up, file.path(analysis.dir, "data","ProfilePlot",paste0("Upregulated_DEGSigFC_",i,".bed")))
    export.bed(region_down, file.path(analysis.dir, "data","ProfilePlot",paste0("Downregulated_DEGSigFC_",i,".bed")))

    regions_matrix <- extract_matrix(coldata=col_data, 
        bed= file.path(analysis.dir, "data","ProfilePlot",paste0("Upregulated_DEGSigFC_",i,".bed")), 
        genome=hg19, up=500, down=500,binSize=10, 
        nthreads=6)
    pdf(file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Upregulated_DEGSigFC_",i,"_500bp.pdf")), height=3,5, width=3.5)
    print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
    dev.off()
    #test
    mean_list <- list()
    for( j in rownames(sample_anno)){
    mean_list[[j]] <- regions_matrix[[j]]
    mean_list[[j]] <- mean(rowMeans(mean_list[[j]]))
    }
    shapiro.test( unlist(mean_list))
    res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
    test <- compare_means(value ~ group, data = res, method="t.test")
    write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Upregulated_DEGSigFC_",i,"allGroups_profileplot_500bp_test.txt")), 
    row.names=FALSE, quote=FALSE, sep="\t")
    #test for center
    for( j in rownames(sample_anno)){
    mean_list[[j]] <- regions_matrix[[j]]
    mean_list[[j]] <- mean(mean_list[[j]]$V50)
    }
    shapiro.test( unlist(mean_list))
    res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
    test <- compare_means(value ~ group, data = res, method="t.test")
    write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Upregulated_DEGSigFC_",i,"allGroups_profileplot_500bp_test_center.txt")), 
        row.names=FALSE, quote=FALSE, sep="\t")


    regions_matrix <- extract_matrix(coldata=col_data, 
        bed= file.path(analysis.dir, "data","ProfilePlot",paste0("Upregulated_DEGSigFC_",i,".bed")), 
        genome=hg19, up=1500, down=1500,binSize=10, 
        nthreads=6)
    pdf(file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Upregulated_DEGSigFC_",i,"_1500bp.pdf")), height=3,5, width=3.5)
    print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
    dev.off()

    regions_matrix <- extract_matrix(coldata=col_data, 
        bed= file.path(analysis.dir, "data","ProfilePlot",paste0("Downregulated_DEGSigFC_",i,".bed")), 
        genome=hg19, up=500, down=500,binSize=10, 
        nthreads=6)
    pdf(file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Downregulated_DEGSigFC_",i,"_500bp.pdf")), height=3,5, width=3.5)
    print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
    dev.off()
    #test
    mean_list <- list()
    for( j in rownames(sample_anno)){
    mean_list[[j]] <- regions_matrix[[j]]
    mean_list[[j]] <- mean(rowMeans(mean_list[[j]]))
    }
    shapiro.test( unlist(mean_list))
    res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
    test <- compare_means(value ~ group, data = res, method="t.test")
    write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Downregulated_DEGSigFC_",i,"allGroups_profileplot_500bp_test.txt")), 
    row.names=FALSE, quote=FALSE, sep="\t")
    #test for center
    for( j in rownames(sample_anno)){
    mean_list[[j]] <- regions_matrix[[j]]
    mean_list[[j]] <- mean(mean_list[[j]]$V50)
    }
    shapiro.test( unlist(mean_list))
    res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
    test <- compare_means(value ~ group, data = res, method="t.test")
    write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Downregulated_DEGSigFC_",i,"allGroups_profileplot_500bp_test_center.txt")), 
        row.names=FALSE, quote=FALSE, sep="\t")
    
    regions_matrix <- extract_matrix(coldata=col_data, 
        bed= file.path(analysis.dir, "data","ProfilePlot",paste0("Downregulated_DEGSigFC_",i,".bed")), 
        genome=hg19, up=1500, down=1500,binSize=10, 
        nthreads=6)
    pdf(file.path(analysis.dir, "ProfilePlots",paste0("allGroups_Downregulated_DEGSigFC_",i,"_1500bp.pdf")), height=3,5, width=3.5)
    print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
    dev.off()

}


#differentially expressed DNA damage repair genes
#loop over different dna damage player subtypes
t2g$TERM <- as.character(t2g$TERM )
for(i in unique(t2g$TERM)){
    t2g_sub <- t2g[t2g$TERM == i, ]
    name <- i 
    name <- gsub(" ", "", name)
    name <- gsub("(", "", name, fixed=TRUE)
    name <- gsub(")", "", name, fixed=TRUE)

    t2g_sub_DEG <- t2g_sub[t2g_sub$GENE %in% DEG_results_list_sign[["ID3-ko_Irr_vs_wt_Irr"]]$ensembl,]
    if(nrow(t2g_sub_DEG)>1){
        region_oi <- prom[prom$gene_id %in% t2g_sub_DEG$GENE, ]
        export.bed(region_oi, file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_DEGSigFC_IrrComp_",name,".bed")))

        regions_matrix <- extract_matrix(coldata=col_data, 
            bed= file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_DEGSigFC_IrrComp_",name,".bed")), 
            genome=hg19, up=500, down=500,binSize=10, 
            nthreads=6)
        pdf(file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_DEGSig_IrrComp_",name,"allGroups_profileplot_500bp.pdf")), height=3,5, width=3.5)
        print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
        dev.off()
        #test
        mean_list <- list()
        for( j in rownames(sample_anno)){
        mean_list[[j]] <- regions_matrix[[j]]
        mean_list[[j]] <- mean(rowMeans(mean_list[[j]]))
        }
        shapiro.test( unlist(mean_list))
        res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
        test <- compare_means(value ~ group, data = res, method="t.test")
        write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_DEGSig_IrrComp_",name,"allGroups_profileplot_500bp_test.txt")), 
        row.names=FALSE, quote=FALSE, sep="\t")
        #test for center
        for( j in rownames(sample_anno)){
        mean_list[[j]] <- regions_matrix[[j]]
        mean_list[[j]] <- mean(mean_list[[j]]$V50)
        }
        shapiro.test( unlist(mean_list))
        res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
        test <- compare_means(value ~ group, data = res, method="t.test")
        write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_DEGSig_IrrComp_",name,"allGroups_profileplot_500bp_test_center.txt")), 
            row.names=FALSE, quote=FALSE, sep="\t")
    }


    t2g_sub_DEG <- t2g_sub[t2g_sub$GENE %in% DEG_results_list_sign[["ID3-ID3-ko_unIrr_vs_wt_unIrr"]]$ensembl,]
    if(nrow(t2g_sub_DEG)>1){
        region_oi <- prom[prom$gene_id %in% t2g_sub_DEG$GENE, ]
        export.bed(region_oi, file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_DEGSigFC_UIrrComp_",name,".bed")))

        regions_matrix <- extract_matrix(coldata=col_data, 
            bed= file.path(analysis.dir, "data","ProfilePlot",paste0("DNA_damage_promoter_DEGSigFC_UIrrComp_",name,".bed")), 
            genome=hg19, up=500, down=500,binSize=10, 
            nthreads=6)
        pdf(file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_DEGSig_UIrrComp_",name,"allGroups_profileplot_500bp.pdf")), height=3,5, width=3.5)
        print(plot_profile(mat_list=regions_matrix, summarizeBy="mean", ci=FALSE, condition="Condition", collapse_reps=TRUE, color=col))
        dev.off()
        #test
        mean_list <- list()
        for( j in rownames(sample_anno)){
        mean_list[[j]] <- regions_matrix[[j]]
        mean_list[[j]] <- mean(rowMeans(mean_list[[j]]))
        }
        shapiro.test( unlist(mean_list))
        res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
        test <- compare_means(value ~ group, data = res, method="t.test")
        write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_DEGSig_UIrrComp_",name,"allGroups_profileplot_500bp_test.txt")), 
        row.names=FALSE, quote=FALSE, sep="\t")
        #test for center
        for( j in rownames(sample_anno)){
        mean_list[[j]] <- regions_matrix[[j]]
        mean_list[[j]] <- mean(mean_list[[j]]$V50)
        }
        shapiro.test( unlist(mean_list))
        res <-data.frame(value=unlist(mean_list), group=sample_anno$Condition)
        test <- compare_means(value ~ group, data = res, method="t.test")
        write.table(as.data.frame(test), file.path(analysis.dir, "ProfilePlots",paste0("DNA_damage_promoter_DEGSig_UIrrComp_",name,"allGroups_profileplot_500bp_test_center.txt")), 
            row.names=FALSE, quote=FALSE, sep="\t")
    }
}
