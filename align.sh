#!/bin/bash
#SBATCH -J STAR_Align
#SBATCH -o ./logs/align/%x_%a_%A.log
#SBATCH -e ./logs/align/%x_%a_%A.err
#SBATCH -p cpu-long
#SBATCH -t 24:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=64G  # Requested Memory
#SBATCH --array=5

module load star/2.7.10b

SPP=$(sed -n ${SLURM_ARRAY_TASK_ID}p spp | cut -f1)
REFDIR=./refs/${SPP}
REF=$(sed -n ${SLURM_ARRAY_TASK_ID}p spp | cut -f3)
GTF=${REF}_genomic.gtf

for file in trimmed/${SPP}/*_R1.trimmed.fastq.gz; do
	prefix=$(echo "$file" | cut -d'_' -f1 | cut -d'/' -f3)
	echo $prefix

	f1=trimmed/${SPP}/${prefix}_R1.trimmed.fastq.gz
	f2=trimmed/${SPP}/${prefix}_R2.trimmed.fastq.gz

	rm -r /scratch/workspace/bbentley_smith_edu-EVE/temp/${prefix}

	STAR --runThreadN 16 \
        --readFilesCommand zcat \
        --genomeDir ${REFDIR} \
	--sjdbGTFfile ${REFDIR}/${REF}_genomic.gtf \
        --outFileNamePrefix alignments/${SPP}/${prefix}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outTmpDir /scratch/workspace/bbentley_smith_edu-EVE/temp/${prefix} \
	--readFilesIn "$f1" "$f2"
done

for file in trimmed/${SPP}/*.SE.trimmed.fastq.gz; do
        prefix=$(echo "$file" | cut -d'.' -f1 | cut -d'/' -f3)
        echo $prefix

        f1=trimmed/${SPP}/${prefix}.SE.trimmed.fastq.gz
       
	rm -r /scratch/workspace/bbentley_smith_edu-EVE/temp/${prefix}

	STAR --runThreadN 16 \
        --readFilesCommand zcat \
        --genomeDir ${REFDIR} \
        --sjdbGTFfile $REFDIR/${REF}_genomic.gtf \
        --outFileNamePrefix alignments/${SPP}/${prefix}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outTmpDir /scratch/workspace/bbentley_smith_edu-EVE/temp/${prefix} \
        --readFilesIn "$f1"
done
