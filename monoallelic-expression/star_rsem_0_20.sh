#!/usr/bin/env bash
#SBATCH -c 8
#SBATCH -t 1:0:0
#SBATCH --mem-per-cpu=16G
#SBATCH -p short
#SBATCH -o /home/am717/slurm_runs/StarRsem.%A.%a.out
#SBATCH -e /home/am717/slurm_runs/StarRsem.%A.%a.err

module load gcc/6.2.0 star/2.7.3a rsem/1.3.0

j=${SLURM_ARRAY_TASK_ID}
# --array=0-20

d=/n/scratch2/am717/sandbox_and_testing/opnsc/
dd=$d/sim_$j
cd $dd

mkdir -p $dd/ref/star
mkdir -p $dd/ref/rsem
mkdir $dd/rsem

cat $d/$j"_mat_transcriptome.fa" $d/$j"_pat_transcriptome.fa" > $dd/ref/star/$j"_transcriptome.fa"
cp $dd/ref/star/$j"_transcriptome.fa" $dd/ref/rsem/$j"_transcriptome.fa"
# reads are here: $dd/reads/$j"_reads.fa"


# 1. STAR:

STAR --runThreadN 8 --runMode genomeGenerate \
     --genomeDir $dd/ref/star/ \
     --genomeFastaFiles $dd/ref/star/$j"_transcriptome.fa"

STAR --readFilesIn $dd/reads/$j"_reads.fa" \
     --outFileNamePrefix $dd/reads/sim_$j"." \
     --runThreadN 4 --outSAMtype BAM Unsorted \
     --genomeDir $dd/ref/star/ \
     --alignEndsType EndToEnd --alignIntronMax 1 --alignIntronMin 2 --scoreDelOpen -10000 --scoreInsOpen -10000 


# 2. RSEM:

rsem-prepare-reference $dd/ref/rsem/$j"_transcriptome.fa" $dd/ref/rsem/$j"_transcriptome"

rsem-calculate-expression --no-qualities --bam $dd/reads/sim_$j".Aligned.out.bam" $dd/ref/rsem/$j"_transcriptome" $dd/rsem/sim_$j".rsem"  




