#load anaconda module
module load anaconda
#load pipelin environment
source activate nf-core-rnaseq-1.2
#go to working directory
cd /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master
#start pipeline
/icgc/dkfzlsdf/analysis/C010/nextflow/nextflow run /icgc/dkfzlsdf/analysis/C010/nextflow/rnaseq-master/main.nf  \
 --saveTrimmed --saveAlignedIntermediates --genome Hg19 -profile cluster_odcf  --skip_genebody_coverage

#copy transcript files into their own folder named by sample
#cd /icgc/dkfzlsdf/analysis/C010/id3/190228_RNAseq_processing/stringtieFPKM/transcripts
#mkdir ../../transcriptsls
#for file in ls *.gtf;do
#filename="${file%.*}"
#mkdir  ../../transcripts/$filename
#cp $file ../../transcripts/$filename/
#done

#run python script to create counts
#cd ../../transcripts
#/icgc/dkfzlsdf/analysis/C010/nextflow/prepDE.py -i ./ -g ./gene_count_matrix.csv -t ./transcript_count_matrix.csv -l 100 


