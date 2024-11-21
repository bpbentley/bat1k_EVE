#!/bin/bash
#SBATCH -J Orthofinder_primary_trancripts
#SBATCH -o ./logs/%x_%a_%A.log
#SBATCH -e ./logs/%x_%a_%A.err
#SBATCH -p cpu-long
#SBATCH -t 24:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --mem=64000  # Requested Memory

ORTHODIR=/home/bbentley_smith_edu/OrthoFinder
REFDIR=./refs/primary_transcripts
OUTDIR=./Ortho_out
SCRATCH=./temp

$ORTHODIR/orthofinder.py -p $SCRATCH -f $REFDIR
