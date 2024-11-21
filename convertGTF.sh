#!/bin/bash
#SBATCH -J gtf2bed
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 05:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=64000  # Requested Memory
#SBATCH --array=11

SPP=$(sed -n ${SLURM_ARRAY_TASK_ID}p spp | cut -f1)
REFDIR=./refs/${SPP}
REF=$(sed -n ${SLURM_ARRAY_TASK_ID}p spp | cut -f3)

module load miniconda/22.11.1-1

conda activate bedops
gtf2bed < ${REFDIR}/${REF}_genomic.gtf > ${REFDIR}/${REF}_genomic.bed
conda deactivate bedops

infer_experiment.py -i alignments/${SPP}/*.bam -r ${REFDIR}/${REF}_genomic.bed 
