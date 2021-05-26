##Count overlap or reads with repeats

#files=`ls icgc/dkfzlsdf/analysis/C010/id3/odcf_bams/*.bam`
#for files 
#featureCounts -p -f -F 'SAF' -a /bigdisk/Nanopore/raw_data/repeats/repeats.gtf \
#-o /home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/results/repeats/repeat_counts.txt $files 


cd icgc/dkfzlsdf/analysis/C010/id3_new/odcf_bams/
for file in ls *.bam;do
featureCounts -p   --tmpDir  ./ -f -F 'SAF' -a /bigdisk/Nanopore/raw_data/repeats/repeats.gtf \
-o /home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/results/repeats/${file}_repeat_counts.txt $file
echo  /home/epicwl/c010-datasets/Internal/2018-Ali/RNASeq/190308_odcf/all_data/results/repeats/${file}_repeat_counts.txt
done