#!/bin/sh
#merge fasta files
cd /media/epicwl/Toshiba/fasta/
#neg_4OHT_wt
cat AS-262605-LR-38364_R1.fastq.gz AS-262605-LR-38365_R1.fastq.gz > neg_4OHT_wt_R1.fastq.gz	
cat AS-262605-LR-38364_R2.fastq.gz AS-262605-LR-38365_R2.fastq.gz > neg_4OHT_wt_R2.fastq.gz

#neg_4OHT_ARID1A_KO
cat AS-262607-LR-38364_R1.fastq.gz AS-262607-LR-38365_R1.fastq.gz > neg_4OHT_ARID1A_KO_R1.fastq.gz
cat AS-262607-LR-38364_R2.fastq.gz AS-262607-LR-38365_R2.fastq.gz > neg_4OHT_ARID1A_KO_R2.fastq.gz

#neg_4OHT_ARIDB_KO
cat AS-262609-LR-38364_R1.fastq.gz AS-262609-LR-38365_R1.fastq.gz > neg_4OHT_ARIDB_KO_R1.fastq.gz
cat AS-262609-LR-38364_R2.fastq.gz AS-262609-LR-38365_R2.fastq.gz > neg_4OHT_ARIDB_KO_R2.fastq.gz

#neg_4OHT_ID3_KO
cat AS-262611-LR-38364_R1.fastq.gz AS-262611-LR-38365_R1.fastq.gz > neg_4OHT_ID3_KO_R1.fastq.gz
cat AS-262611-LR-38364_R2.fastq.gz AS-262611-LR-38365_R2.fastq.gz > neg_4OHT_ID3_KO_R2.fastq.gz

#neg_4OHT_siMDC1
cat AS-262613-LR-38364_R1.fastq.gz AS-262613-LR-38365_R1.fastq.gz > neg_4OHT_siMDC1_R1.fastq.gz
cat AS-262613-LR-38364_R2.fastq.gz AS-262613-LR-38365_R2.fastq.gz > neg_4OHT_siMDC1_R2.fastq.gz	

#neg_4OHT_ID3_KO_siMD
cat AS-262615-LR-38364_R1.fastq.gz AS-262615-LR-38365_R1.fastq.gz > neg_4OHT_ID3_KO_siMD_R1.fastq.gz
cat AS-262615-LR-38364_R2.fastq.gz AS-262615-LR-38365_R2.fastq.gz > neg_4OHT_ID3_KO_siMD_R2.fastq.gz

#pos_4OHT_wt
cat AS-262617-LR-38364_R1.fastq.gz AS-262617-LR-38365_R1.fastq.gz > pos_4OHT_wt_R1.fastq.gz	
cat AS-262617-LR-38364_R2.fastq.gz AS-262617-LR-38365_R2.fastq.gz > pos_4OHT_wt_R2.fastq.gz

#pos_4OHT_ARID1A_KO
cat AS-262619-LR-38364_R1.fastq.gz AS-262619-LR-38365_R1.fastq.gz > pos_4OHT_ARID1A_KO_R1.fastq.gz
cat AS-262619-LR-38364_R2.fastq.gz AS-262619-LR-38365_R2.fastq.gz > pos_4OHT_ARID1A_KO_R2.fastq.gz

#pos_4OHT_ARIDB_KO
cat AS-262621-LR-38364_R1.fastq.gz AS-262621-LR-38365_R1.fastq.gz > pos_4OHT_ARIDB_KO_R1.fastq.gz
cat AS-262621-LR-38364_R2.fastq.gz AS-262621-LR-38365_R2.fastq.gz > pos_4OHT_ARIDB_KO_R2.fastq.gz	

#pos_40HT_IDR_KO
cat AS-262623-LR-38364_R1.fastq.gz AS-262623-LR-38365_R1.fastq.gz > pos_40HT_IDR_KO_R1.fastq.gz
cat AS-262623-LR-38364_R2.fastq.gz AS-262623-LR-38365_R2.fastq.gz > pos_40HT_IDR_KO_R2.fastq.gz

#pos_4OHT_siMDC1
cat AS-262625-LR-38364_R1.fastq.gz AS-262625-LR-38365_R1.fastq.gz > pos_4OHT_siMDC1_R1.fastq.gz
cat AS-262625-LR-38364_R2.fastq.gz AS-262625-LR-38365_R2.fastq.gz > pos_4OHT_siMDC1_R2.fastq.gz

#pos_4OHT_IDR3_KO_siM
cat AS-262627-LR-38364_R1.fastq.gz AS-262627-LR-38365_R1.fastq.gz > pos_4OHT_IDR3_KO_siM_R1.fastq.gz	
cat AS-262627-LR-38364_R2.fastq.gz AS-262627-LR-38365_R2.fastq.gz > pos_4OHT_IDR3_KO_siM_R2.fastq.gz
#!/bin/sh




#neg_4OHT_wt
mkdir ATAC_Processing/Ali/neg_4OHT_wt

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/neg_4OHT_wt_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/neg_4OHT_wt_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/neg_4OHT_wt \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#neg_4OHT_ARID1A_KO
mkdir ATAC_Processing/Ali/neg_4OHT_ARID1A_KO

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/neg_4OHT_ARID1A_KO_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/neg_4OHT_ARID1A_KO_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/neg_4OHT_ARID1A_KO \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#neg_4OHT_ARID1B_KO
mkdir ATAC_Processing/Ali/neg_4OHT_ARID1V_KO

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/neg_4OHT_ARID1B_KO_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/neg_4OHT_ARID1B_KO_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/neg_4OHT_ARID1B_KO \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#neg_4OHT_ID3_KO
mkdir ATAC_Processing/Ali/neg_4OHT_ID3_KO

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/neg_4OHT_ID3_KO_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/neg_4OHT_ID3_KO_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/neg_4OHT_ID3_KO \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr


#neg_4OHT_siMDC1
mkdir ATAC_Processing/Ali/neg_4OHT_siMDC1

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/neg_4OHT_siMDC1_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/neg_4OHT_siMDC1_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/neg_4OHT_siMDC1 \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#neg_4OHT_ID3_KO_siMD
mkdir ATAC_Processing/Ali/neg_4OHT_ID3_KO_siMD

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/neg_4OHT_ID3_KO_siMD_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/neg_4OHT_ID3_KO_siMD_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/neg_4OHT_ID3_KO_siMD \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#pos_4OHT_wt
mkdir ATAC_Processing/Ali/pos_4OHT_wt

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/pos_4OHT_wt_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/pos_4OHT_wt_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/pos_4OHT_wt \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#pos_4OHT_ARID1A_KO
mkdir ATAC_Processing/Ali/pos_4OHT_ARID1A_KO

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/pos_4OHT_ARID1A_KO_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/pos_4OHT_ARID1A_KO_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/pos_4OHT_ARID1A_KO \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#pos_4OHT_ARIDB_KO
mkdir ATAC_Processing/Ali/pos_4OHT_ARIDB_KO

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/pos_4OHT_ARIDB_KO_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/pos_4OHT_ARIDB_KO_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/pos_4OHT_ARIDB_KO \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#pos_40HT_IDR_KO
mkdir ATAC_Processing/Ali/pos_40HT_IDR_KO

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/pos_40HT_IDR_KO_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/pos_40HT_IDR_KO_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/pos_40HT_IDR_KO \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr

#pos_4OHT_siMDC1
mkdir ATAC_Processing/Ali/pos_4OHT_siMDC1

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/pos_4OHT_siMDC1_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/pos_4OHT_siMDC1_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/pos_4OHT_siMDC1 \
-auto_detect_adapter  \
-no_ataqc \
#-enable_idr

#pos_4OHT_IDR3_KO_siM
mkdir ATAC_Processing/Ali/pos_4OHT_IDR3_KO_siM

bds atac_dnase_pipelines/atac.bds \
-species hg19 \
-fastq1_1  /media/epicwl/Toshiba/fasta/pos_4OHT_IDR3_KO_siM_R1.fastq.gz \
-fastq1_2  /media/epicwl/Toshiba/fasta/pos_4OHT_IDR3_KO_siM_R2.fastq.gz \
-bwt2_idx atac_dnase_pipelines/genomes/hg19/bowtie2_index/male.hg19.fa \
-gensz hs \
-chrsz  atac_dnase_pipelines/genomes/hg19/hg19.chrom.sizes \
-out_dir ATAC_Processing/Ali/pos_4OHT_IDR3_KO_siM \
-auto_detect_adapter \
-no_ataqc \
#-enable_idr
