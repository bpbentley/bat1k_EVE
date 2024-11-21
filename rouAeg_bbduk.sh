#!/bin/bash
#SBATCH -J rouAeg_bbduk
#SBATCH -o ./logs/bbduk/%x_%a_%A.log
#SBATCH -e ./logs/bbduk/%x_%a_%A.err
#SBATCH -p cpu-preempt
#SBATCH -t 02:00:00
#SBATCH -c 8
#SBATCH --mem=32000
#SBATCH --array=1-4

module load BBMap/38.90
module load fastqc/0.11.9
module load sratoolkit/3.0.7

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./sample_files/rouAeg_s)
SPP=rouAeg
TMPDIR=/scratch/workspace/bbentley_smith_edu-EVE

#fasterq-dump ${SAMPLE} -o ${TMPDIR}/raw_data/${SAMPLE} -t $TMPDIR

#gzip $TMPDIR/raw_data/${SPP}/${SAMPLE}_1.fastq
#gzip $TMPDIR/raw_data/${SPP}/${SAMPLE}_2.fastq

#fastqc $TMPDIR/raw_data/${SPP}/${SAMPLE}_1.fastq.gz -o ./FastQC/raw/
#fastqc $TMPDIR/raw_data/${SPP}/${SAMPLE}_2.fastq.gz -o ./FastQC/raw/

bbduk.sh \
ref=./adapters.fa \
in1=$TMPDIR/raw_data/${SPP}/${SAMPLE}_1.fastq.gz \
in2=$TMPDIR/raw_data/${SPP}/${SAMPLE}_2.fastq.gz \
out1=./trimmed_reads/${SPP}/${SAMPLE}_R1.trimmed.fastq \
out2=./trimmed_reads/${SPP}/${SAMPLE}_R2.trimmed.fastq \
trimpolya=10 ktrim=r k=23 mink=11 hdist=1 tpe tbo

gzip ./trimmed_reads/${SPP}/${SAMPLE}_R1.trimmed.fastq
gzip ./trimmed_reads/${SPP}/${SAMPLE}_R2.trimmed.fastq

fastqc ./trimmed_reads/${SPP}/${SAMPLE}_R1.trimmed.fastq.gz -o ./FastQC/trimmed/
fastqc ./trimmed_reads/${SPP}/${SAMPLE}_R2.trimmed.fastq.gz -o ./FastQC/trimmed/
