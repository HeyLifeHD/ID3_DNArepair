
for file in `ls /icgc/dkfzlsdf/project/C010/id3/sequencing/rna_sequencing/core/*/ | tr "LR" "\t" | cut -f1 | uniq`; do
#echo /icgc/dkfzlsdf/project/C010/id3/sequencing/rna_sequencing/core/*/$file*R1.fastq.gz 
cat /icgc/dkfzlsdf/project/C010/id3/sequencing/rna_sequencing/core/*/$file*R1.fastq.gz > /icgc/dkfzlsdf/analysis/C010/id3/merged_fastq/$file.R1.fastq.gz
cat /icgc/dkfzlsdf/project/C010/id3/sequencing/rna_sequencing/core/*/$file*R2.fastq.gz > /icgc/dkfzlsdf/analysis/C010/id3/merged_fastq/$file.R2.fastq.gz
echo $file
done