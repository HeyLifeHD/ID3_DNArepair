cd c010-datasets/Internal/2018-Ali/ATAC/

fluff profile -i  chr1:68602071-68612071 -d \
processing/set1/pos_4OHT_wt/align/rep1/pos_4OHT_wt_R1.trim.PE2SE.nodup.bam \
processing/set1/pos_4OHT_wt/align/rep1/pos_4OHT_wt_R1.trim.PE2SE.nodup.bam \
-o locus_plots/lo

cd ~/
fluff heatmap -f c010-datasets/Internal/2018-Ali/ATAC/data/ASIsites_100BEST_hg19.bed \
-d c010-datasets/Internal/2018-Ali/ATAC/processing/set1/neg_4OHT_wt/align/rep1/neg_4OHT_wt_R1.trim.PE2SE.nodup.bam \
c010-datasets/Internal/2018-Ali/ATAC/processing/set1/pos_4OHT_wt/align/rep1/pos_4OHT_wt_R1.trim.PE2SE.nodup.bam \
-o c010-datasets/Internal/2018-Ali/ATAC/locus_plots/heatmap_ex.pdf
