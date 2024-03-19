#!/bin/bash
#SBATCH --job-name=prepyscenic
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=70G
#SBATCH --time=12:00:00
#SBATCH --array=0-7
#SBATCH --output=pysceniclogs/job.%A.%a.out

samples=("Aguirre_EGSFR0074" "Aguirre_EGSFR0128" "Aguirre_EGSFR0148" "Aguirre_EGSFR1982" "Aguirre_EGSFR1938" "CCG1153_4411" "Aguirre_EGSFR2218" "Aguirre_EGSFR1732")

python run-pre-pyscenic-script.py --samplename ${samples[${SLURM_ARRAY_TASK_ID}]} --ncpus 8
