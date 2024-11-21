#!/bin/bash
#SBATCH -J HTSeq_homSap
#SBATCH -o ./logs/HTSeq/homSap/%x_%a_%A.log
#SBATCH -e ./logs/HTSeq/homSap/%x_%a_%A.err
#SBATCH -p cpu-long
#SBATCH -t 24:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=64000  # Requested Memory
#SBATCH --array=2

NUM=5
SPP=$(sed -n ${NUM}p spp | cut -f1)
REFDIR=./refs/${SPP}
REF=$(sed -n ${NUM}p spp | cut -f3)
OUTDIR=./htseq/${SPP}
SCRATCH=/scratch/workspace/bbentley_smith_edu-EVE/temp
STR=$(sed -n ${NUM}p spp | cut -f4)

#pip install HTSeq

#for file in alignments/${SPP}/*.bam; do
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p sample_files/homSap_s)
file=alignments/${SPP}/${SAMPLE}_Aligned.sortedByCoord.out.bam
SAMPLE=$(echo "$file" | sed "s|alignments/${SPP}/||; s|_Aligned.sortedByCoord.out.bam||")
/home/bbentley_smith_edu/.local/bin/htseq-count -f bam -r pos -s ${STR} --idattr=gene_id \
    ./alignments/${SPP}/${SAMPLE}_Aligned.sortedByCoord.out.bam \
    ${REFDIR}/${REF}_genomic.gtf > htseq/$SPP/${SAMPLE}_GeneCount.txt
#done
