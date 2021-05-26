#copy feature counts
cp  /icgc/dkfzlsdf/project/C010/id3/sequencing/rna_sequencing/view-by-pid/*/*/paired/merged-alignment/featureCounts/*.fpkm_tpm.featureCounts.tsv \
/icgc/dkfzlsdf/analysis/C010/id3/odcf_featureCounts/

#copy bam and bai
cp /icgc/dkfzlsdf/project/C010/id3/sequencing/rna_sequencing/view-by-pid/*/*/paired/merged-alignment/*_merged.mdup.bam \
/icgc/dkfzlsdf/analysis/C010/id3/odcf_bams/
cp /icgc/dkfzlsdf/project/C010/id3/sequencing/rna_sequencing/view-by-pid/*/*/paired/merged-alignment/*_merged.mdup.bam.bai   \
/icgc/dkfzlsdf/analysis/C010/id3/odcf_bams/

rm -f /icgc/dkfzlsdf/analysis/C010/id3/odcf_bams/*chimeric*