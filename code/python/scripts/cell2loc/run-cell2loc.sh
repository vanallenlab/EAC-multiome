#!/bin/bash
#SBATCH --job-name=cell2Loc
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=70G
#SBATCH --time=12:00:00
#SBATCH --array=0-0
#SBATCH --output=cell2loc_logs/job.%A.%a.out

#patients=("EGSFR0074_A" "EGSFR0148" "EGSFR1938_A" "EGSFR1938_B" "EGSFR1938_C" "EGSFR1982")
patients=("EGSFR1982")

python script-cell2loc.py --patient ${patients[${SLURM_ARRAY_TASK_ID}]}
