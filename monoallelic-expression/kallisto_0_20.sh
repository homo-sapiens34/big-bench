#!/usr/bin/env bash
#SBATCH -c 8
#SBATCH -t 10:0:0
#SBATCH --mem-per-cpu=8G
#SBATCH -p short
#SBATCH -o /home/am717/slurm_runs/Kallisto.%A.%a.out
#SBATCH -e /home/am717/slurm_runs/Kallisto.%A.%a.err

module load gcc/6.2.0 kallisto/0.45.1

j=${SLURM_ARRAY_TASK_ID}
# --array=0-20

d=/n/scratch2/am717/sandbox_and_testing/opnsc/
dd=$d/sim_$j
mkdir $dd

cd $dd

mkdir -p $dd/ref/kallisto
mkdir $dd/reads
mkdir $dd/kallisto

cat $d/$j"_mat_transcriptome.fa" $d/$j"_pat_transcriptome.fa" > $dd/ref/kallisto/$j"_transcriptome.fa"
mv $d/$j"_reads.fa" $dd/reads/$j"_reads.fa"


kallisto index -i $dd/ref/kallisto/$j"_transcriptome.31.idx" $dd/ref/kallisto/$j"_transcriptome.fa"


kallisto quant --single -i $dd/ref/kallisto/$j"_transcriptome.31.idx" -o $dd/kallisto $dd/reads/$j"_reads.fa" -l 100 -s 1

cp $dd/kallisto/abundance.tsv $dd/kallisto/sim_$j".kallisto.abundance.tsv"



