#!/bin/bash
#SBATCH -J canFam_trim
#SBATCH -o ./logs/bbduk/%x_%a_%A.log
#SBATCH -e ./logs/bbduk/%x_%a_%A.err
#SBATCH -p cpu-preempt
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH --mem=32000
#SBATCH --array=1-3

module load BBMap/38.90
module load fastqc/0.11.9
module load sratoolkit/3.0.7

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./sample_files/felCat_s)
SPP=felCat
TMPDIR=/scratch/workspace/bbentley_smith_edu-EVE

gzip $TMPDIR/raw_data/${SPP}/${SAMPLE}.fastq

fastqc $TMPDIR/raw_data/${SPP}/${SAMPLE}.fastq.gz -o ./FastQC/raw/

bbduk.sh \
ref=./adapters.fa \
in=$TMPDIR/raw_data/${SPP}/${SAMPLE}.fastq.gz \
out=./trimmed_reads/${SPP}/${SAMPLE}.trimmed.fastq \
trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo

gzip ./trimmed_reads/${SPP}/${SAMPLE}.trimmed.fastq

fastqc ./trimmed_reads/${SPP}/${SAMPLE}.trimmed.fastq.gz -o ./FastQC/trimmed/
