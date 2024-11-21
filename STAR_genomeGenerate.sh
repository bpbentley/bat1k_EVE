#!/bin/bash
#SBATCH -J STARgenome
#SBATCH -o ./logs/%x_%a.log
#SBATCH -e ./logs/%x_%a.err
#SBATCH -p cpu
#SBATCH -t 05:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=180GB  # Requested Memory
#SBATCH --array=11

module load star/2.7.10b

SPP=$(sed -n ${SLURM_ARRAY_TASK_ID}p spp | cut -f1)
echo $SPP
LEN=$(sed -n ${SLURM_ARRAY_TASK_ID}p spp | cut -f2)
REFDIR=./refs/${SPP}

#gunzip $REFDIR/*.fna.gz
#gunzip $REFDIR/*.gtf.gz

REF=${REFDIR}/*.fna
GTF=${REFDIR}/*.gtf

STAR --runMode genomeGenerate --genomeDir $REFDIR --genomeFastaFiles $REF --outTmpDir /scratch/workspace/bbentley_smith_edu-EVE/temp/$SPP --sjdbGTFfile $GTF --sjdbOverhang $LEN --runThreadN 16 --limitGenomeGenerateRAM 164436989536
