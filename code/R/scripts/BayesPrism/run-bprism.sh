#!/bin/bash
 
#SBATCH --job-name=bprism
#SBATCH --output=output_eac.txt
#SBATCH --cpus-per-task=25
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=10G
 
source ~/.bashrc
enable_modules
module load python/3.10.2 scipy-stack/2022a rstudio-server/4.1
VENV=/path/to/rstudio_env
source $VENV/bin/activate
export R_LIBS_USER=$VENV
Rscript /path/to/run_BPrism.R
