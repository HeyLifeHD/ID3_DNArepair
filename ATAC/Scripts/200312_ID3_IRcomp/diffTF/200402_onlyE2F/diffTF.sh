#copy atac data do c010 datasets
cd /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF
mkdir data
cd data
for file in `ls /icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/*/*_tn5_center_73bp.bam | grep "cache" -v`; do
    cp $file .
    echo $file
done 
for file in `ls /icgc/dkfzlsdf/analysis/C010/cwlab_processing_out/canepi-srv1/20200310_107_ID3_ATAC/runs_out/*/*_tn5_center_73bp.macs2_peaks.broadPeak| grep "cache" -v`; do
    cp $file .
    echo $file
done 

mkdir /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_onlyE2F

conda activate base
snakemake --snakefile /home/heyj/tools/diffTF/src/Snakefile  --cores 5 --configfile ./200402_input_onlyE2F/config.json \
--use-singularity --singularity-args "--bind ./,/home/heyj/tools/diffTF/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_input_onlyE2F,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_onlyE2F" \
--singularity-prefix /home/heyj/tools/singularity_containers
snakemake --snakefile /home/heyj/tools/diffTF/src/Snakefile  --cores 5 --configfile ./200402_input_onlyE2F/config.json \
--use-singularity --singularity-args "--bind ./,/home/heyj/tools/diffTF/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_input_onlyE2F,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_onlyE2F_anal" \
--singularity-prefix /home/heyj/tools/singularity_containers


snakemake --snakefile /home/heyj/tools/diffTF/src/Snakefile  --cores 5 --configfile ./200402_input_all/config.json \
--use-singularity --singularity-args "--bind ./,/home/heyj/tools/diffTF/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_input_onlyE2F,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all" \
--singularity-prefix /home/heyj/tools/singularity_containers


mkdir /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_anl
snakemake --snakefile /home/heyj/tools/diffTF/src/Snakefile  --cores 5 --configfile ./200402_input_all/config.json \
--use-singularity --singularity-args "--bind /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/,/home/heyj/tools/diffTF/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_input_onlyE2F,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_anl/" \
--singularity-prefix /home/heyj/tools/singularity_containers

#all
snakemake --snakefile /home/heyj/tools/diffTF/src/Snakefile  --cores 5 --configfile ./200402_input_all/config.json \
--use-singularity --singularity-args "--bind /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/,/home/heyj/tools/diffTF/" \
--singularity-prefix /home/heyj/tools/singularity_containers

#same for unirradiated
snakemake --snakefile /home/heyj/tools/diffTF/src/Snakefile  --cores 5 --configfile ./200402_input_onlyE2F_unIrr/config.json \
--use-singularity --singularity-args "--bind ./,/home/heyj/tools/diffTF/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_input_onlyE2F,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_onlyE2F_unIrr" \
--singularity-prefix /home/heyj/tools/singularity_containers

#all
mkdir /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_unIrr_anal
snakemake --snakefile /home/heyj/tools/diffTF/src/Snakefile  --cores 5 --configfile ./200402_input_all_unIrr/config.json \
--use-singularity --singularity-args "--bind /home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/,/home/heyj/tools/diffTF/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_input_all_unIrr/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/data/,/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_all_unIrr_anal/" \
--singularity-prefix /home/heyj/tools/singularity_containers

/home/heyj/c010-datasets/Internal/2018-Ali/ATAC/200312_ID3_IRcomp/analysis/diffTF/200402_output_onlyE2F_anal/TEMP/extension100/wt_IrrvsID3-ko_Irr.all.allMotifs_log2FC_perm20.tsv.gz