
#Directories
#NGS-machine
.libPaths("/home/c010-ngs/R/x86_64-pc-linux-gnu-library/3.4_alt")
base.dir <- "/home/c010-ngs/c010-datasetshey/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")

#VM
base.dir <- "/C010-datasets/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")

#PC-09
base.dir <- "/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/"
data.dir <- file.path(base.dir, "data")
analysis.dir <- file.path(base.dir, "181002_1rep")
temp.dir <- "/home/heyj/rTemp"

#Libraries
library(DiffBind)
#Create Sample sheet
sample_anno <- data.frame("SampleID"=NA,
"Treatment"=c(rep("neg",6,),rep("pos",6,)), "Genotype"=rep(c("wt","ARID1A_KD","ARID1B_KD", "ID3_KD", "siMDC1", "ID3_KD_siMDC1"),2),
"Replicate"=rep("1",12), "lane"=rep("1",12),"Peaks"=NA, "PeakCaller"=rep("macs", 4) , "bamReads"=NA)
sample_anno$PeakCaller<- as.character(sample_anno$PeakCaller)
sample_anno$SampleID<- paste0(sample_anno$Treatment, "_4OHT_", sample_anno$Genotype, "_rep_", sample_anno$Replicate)
rownames(sample_anno)<- sample_anno$SampleID

#Peaks
sample_anno[1,6]<- file.path(base.dir,"processing", "set1", "neg_4OHT_wt", "peak", "macs2","rep1", "neg_4OHT_wt_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[2,6]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ARID1A_KO", "peak", "macs2","rep1", "neg_4OHT_ARID1A_KO_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[3,6]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ARID1B_KO", "peak", "macs2","rep1", "neg_4OHT_ARID1B_KO_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[4,6]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ID3_KO", "peak", "macs2","rep1", "neg_4OHT_ID3_KO_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[5,6]<- file.path(base.dir,"processing", "set1", "neg_4OHT_siMDC1", "peak", "macs2","rep1", "neg_4OHT_siMDC1_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[6,6]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ID3_KO_siMD", "peak", "macs2","rep1", "neg_4OHT_ID3_KO_siMD_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")

sample_anno[7,6]<- file.path(base.dir,"processing", "set1", "pos_4OHT_wt", "peak", "macs2","rep1", "pos_4OHT_wt_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[8,6]<- file.path(base.dir,"processing", "set1", "pos_4OHT_ARID1A_KO", "peak", "macs2","rep1", "pos_4OHT_ARID1A_KO_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[9,6]<- file.path(base.dir,"processing", "set1", "pos_4OHT_ARIDB_KO", "peak", "macs2","rep1", "pos_4OHT_ARIDB_KO_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[10,6]<- file.path(base.dir,"processing", "set1", "pos_40HT_IDR_KO", "peak", "macs2","rep1", "pos_40HT_IDR_KO_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[11,6]<- file.path(base.dir,"processing", "set1", "pos_4OHT_siMDC1", "peak", "macs2","rep1", "pos_4OHT_siMDC1_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")
sample_anno[12,6]<- file.path(base.dir,"processing", "set1", "pos_4OHT_IDR3_KO_siM", "peak", "macs2","rep1", "pos_4OHT_IDR3_KO_siM_R1.trim.PE2SE.nodup.tn5.pf_peaks.xls")

#BAM
sample_anno[1,8]<- file.path(base.dir,"processing", "set1", "neg_4OHT_wt", "align","rep1", "neg_4OHT_wt_R1.trim.PE2SE.nodup.bam")
sample_anno[2,8]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ARID1A_KO", "align","rep1", "neg_4OHT_ARID1A_KO_R1.trim.PE2SE.nodup.bam")
sample_anno[3,8]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ARID1B_KO", "align","rep1", "neg_4OHT_ARID1B_KO_R1.trim.PE2SE.nodup.bam")
sample_anno[4,8]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ID3_KO", "align","rep1", "neg_4OHT_ID3_KO_R1.trim.PE2SE.nodup.bam")
sample_anno[5,8]<- file.path(base.dir,"processing", "set1", "neg_4OHT_siMDC1", "align","rep1", "neg_4OHT_siMDC1_R1.trim.PE2SE.nodup.bam")
sample_anno[6,8]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ID3_KO_siMD","align","rep1", "neg_4OHT_ID3_KO_siMD_R1.trim.PE2SE.nodup.bam")

sample_anno[7,8]<- file.path(base.dir,"processing", "set1", "pos_4OHT_wt", "align","rep1", "pos_4OHT_wt_R1.trim.PE2SE.nodup.bam")
sample_anno[8,8]<- file.path(base.dir,"processing", "set1", "pos_4OHT_ARID1A_KO", "align","rep1", "pos_4OHT_ARID1A_KO_R1.trim.PE2SE.nodup.bam")
sample_anno[9,8]<- file.path(base.dir,"processing", "set1", "pos_4OHT_ARIDB_KO", "align","rep1", "pos_4OHT_ARIDB_KO_R1.trim.PE2SE.nodup.bam")
sample_anno[10,8]<- file.path(base.dir,"processing", "set1", "pos_40HT_IDR_KO", "align","rep1", "pos_40HT_IDR_KO_R1.trim.PE2SE.nodup.bam")
sample_anno[11,8]<- file.path(base.dir,"processing", "set1", "pos_4OHT_siMDC1", "align","rep1", "pos_4OHT_siMDC1_R1.trim.PE2SE.nodup.bam")
sample_anno[12,8]<- file.path(base.dir,"processing", "set1", "pos_4OHT_IDR3_KO_siM", "align","rep1", "pos_4OHT_IDR3_KO_siM_R1.trim.PE2SE.nodup.bam")

#bw 
sample_anno$bw <- NA

sample_anno[1,9]<- file.path(base.dir,"processing", "set1", "neg_4OHT_wt", "signal", "macs2","rep1", "neg_4OHT_wt_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[2,9]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ARID1A_KO", "signal", "macs2","rep1", "neg_4OHT_ARID1A_KO_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[3,9]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ARID1B_KO", "signal", "macs2","rep1", "neg_4OHT_ARID1B_KO_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[4,9]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ID3_KO", "signal", "macs2","rep1", "neg_4OHT_ID3_KO_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[5,9]<- file.path(base.dir,"processing", "set1", "neg_4OHT_siMDC1", "signal", "macs2","rep1", "neg_4OHT_siMDC1_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[6,9]<- file.path(base.dir,"processing", "set1", "neg_4OHT_ID3_KO_siMD", "signal", "macs2","rep1", "neg_4OHT_ID3_KO_siMD_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")

sample_anno[7,9]<- file.path(base.dir,"processing", "set1", "pos_4OHT_wt", "signal", "macs2","rep1", "pos_4OHT_wt_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[8,9]<- file.path(base.dir,"processing", "set1", "pos_4OHT_ARID1A_KO", "signal", "macs2","rep1", "pos_4OHT_ARID1A_KO_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[9,9]<- file.path(base.dir,"processing", "set1", "pos_4OHT_ARIDB_KO", "signal", "macs2","rep1", "pos_4OHT_ARIDB_KO_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[10,9]<- file.path(base.dir,"processing", "set1", "pos_40HT_IDR_KO", "signal", "macs2","rep1", "pos_40HT_IDR_KO_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[11,9]<- file.path(base.dir,"processing", "set1", "pos_4OHT_siMDC1", "signal", "macs2","rep1", "pos_4OHT_siMDC1_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
sample_anno[12,9]<- file.path(base.dir,"processing", "set1", "pos_4OHT_IDR3_KO_siM", "signal", "macs2","rep1", "pos_4OHT_IDR3_KO_siM_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")


#add total reads
#sample_anno$totalReads <- c("76889662","59725664","?????", "61569492", "80869376","?????")
#save
saveRDS(sample_anno, file.path(data.dir, "sample_anno.rds"))

#create dataset
dataset <-dba(sampleSheet=sample_anno, peakCaller="macs")

#correlation Heatmap
pdf(file.path(analysis.dir, "DiffBind","correlation_heatmap.pdf"))
plot(dataset)
dev.off()

#count reads for each peak
dataset<- dba.count(dataset)

#correlation Heatmap
pdf(file.path(analysis.dir, "DiffBind", "correlation_heatmap_affinityscore.pdf"))
plot(dataset)
dev.off()

#setup contrast
dataset <- dba.contrast(dataset, categories=DBA_CONDITION,minMembers=2,block=DBA_REPLICATE)

#DE Analysis
dataset<- dba.analyze(dataset)
saveRDS(dataset, file.path(analysis.dir, "DiffBind","Dataset.rds"))

pdf(file.path(analysis.dir, "DiffBind", "correlation_heatmap_DE.pdf"))
plot(dataset, contrast=1)
dev.off()

#get DE regions
datasets.db <- dba.report(dataset, contrast=1,method =DBA_DESEQ2_BLOCK, th=0.1)

pdf(file.path(analysis.dir, "DiffBind", "pca_affinity_score.pdf"))
dba.plotPCA(dataset,DBA_CONDITION, label=DBA_CONDITION)
dev.off()

pdf(file.path(analysis.dir, "DiffBind", "MA_plot.pdf"))
dba.plotMA(dataset, method= DBA_DESEQ2_BLOCK)
dev.off()

#annotation
require(TxDb.Mmusculus.UCSC.hg19.knownGene)
txdb <- TxDb.Mmusculus.UCSC.hg19.knownGene
library(ChIPseeker)
#annotate DE regions
db.anno <- annotatePeak(peak= datasets.db, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
db.annodf <- as.data.frame(db.anno)
saveRDS(db.annodf, file.path(analysis.dir, "DiffBind", "DiffBind_DB_anno.rds"))
write.table(db.annodf, file=file.path(analysis.dir, "DiffBind","DiffBind_DB_anno.txt"), sep="\t", quote=FALSE, row.names = FALSE)
db.annodf<- readRDS(file.path(analysis.dir, "DiffBind", "DiffBind_DB_anno.rds"))
db.annogr <- makeGRangesFromDataFrame(db.annodf)

db.annodf<- readRDS(file.path(analysis.dir, "DiffBind", "DiffBind_DB_anno.rds"))
pdf(file.path(analysis.dir,  "DiffBind", "DB_DistanceTSSPC.pdf"))
plotDistToTSS(db.anno , 
              title="Distribution of DMR relative to TSS")
dev.off()  

pdf(file.path(analysis.dir, "DiffBind","DB_AnnoDist_Pie.pdf"))
plotAnnoPie(db.anno )
dev.off()

pdf(file.path(analysis.dir, "DiffBind","DB_AnnoDist_Upset.pdf"))
upsetplot(db.anno , vennpie=TRUE)
dev.off()

#functional enrichment
library(clusterProfiler)
genes =db.annodf$geneId
names(genes) = sub("_", "\n", names(genes))

GO_all<- enrichGO(gene= genes, OrgDb = org.Mm.eg.db,  
                 pvalueCutoff = 1, qvalueCutoff = 1, ont = "All")
GO_all <- setReadable(GO_all, OrgDb = org.Mm.eg.db)

GO_all.df <- as.data.frame(GO_all)

write.table(GO_all.df, file=file.path(analysis.dir, "DiffBind","GO_all.txt"), sep="\t", quote=FALSE, row.names = FALSE)

pdf(file.path(analysis.dir, "DiffBind","GO_dot.pdf"))
dotplot(GO_all, showCategory=30, colorBy = "pvalue",  font.size=8)
dev.off()

pdf(file.path(analysis.dir, "DiffBind","GO_bar.pdf"))
barplot(GO_all, drop=TRUE, showCategory = 30)
dev.off()

KEGG<-enrichKEGG(gene= genes, organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1)

pdf(file.path(analysis.dir, "DiffBind","KEGG_dot.pdf"))
dotplot(KEGG , showCategory=30, colorBy = "pvalue",  font.size=8)
dev.off()

pdf(file.path(analysis.dir, "DiffBind","KEGG_bar.pdf"))
barplot(KEGG, drop=TRUE, showCategory = 30)
dev.off()

KEGG.df <- as.data.frame(KEGG)

write.table(KEGG.df, file=file.path(analysis.dir, "DiffBind","KEGG.txt"), sep="\t", quote=FALSE, row.names = FALSE)










#get overlap
dataset_new <-dba(sampleSheet=sample_anno, peakCaller="macs")

dataset_new_1 <- dba.peakset(dataset_new, consensus=DBA_CONDITION, minOverlap=0.33)


#plot unique 
pdf(file.path(analysis.dir, "dba_unique_peaks_venn.pdf"))
dba.plotVenn(dataset_new_1,dataset_new_1$masks$Consensus)
dev.off()

#get consensus to find unique
dataset_new_1.OL <- dba.overlap(dataset_new_1, dataset_new_1$masks$Consensus)

wt_peak_set <-dataset_new_1.OL$onlyA
tg_peak_set <-dataset_new_1.OL$onlyB
consens_peak_set <- dataset_new_1.OL$inAll

dataset_new_1.OL.rep <- dba.report(dataset_new_1 ,bCalled=TRUE,th=1)j