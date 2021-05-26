##Count overlap or reads with repeats
files=`ls icgc/analysis/C010/id3/190228_RNAseq_processing/STAR/*.bam`
featureCounts -p -f -F 'SAF' -a /bigdisk/Nanopore/raw_data/repeats/repeats.gtf \
-o c010-datasets/Internal/2018-Ali/RNASeq/190305_unmerged/results/repeats/repeat_counts.txt $files 

