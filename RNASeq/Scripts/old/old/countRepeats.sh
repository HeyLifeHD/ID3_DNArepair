##Count overlap or reads with repeats
files=`ls c010-datasets/Internal/2018-Ali/RNASeq/181114_analysis_rep1/markDuplicates/*.bam`
featureCounts -p -f -F 'SAF' -a /bigdisk/Nanopore/raw_data/repeats/repeats.gtf \
-o c010-datasets/Internal/2018-Ali/RNASeq/181114_analysis_rep1/repeat_counts.txt $files 

