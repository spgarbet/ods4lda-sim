#!/bin/bash

#SBATCH --mail-user=shawn.garbett@vumc.org
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:30:00
#SBATCH --mem=4G
#SBATCH --output=status/job%a.out
#SBATCH --array=1-2000
#SBATCH --constrain=haswell

module purge
module load Intel/2016.3.210-GCC-5.4.0-2.26
module load IntelMPI/5.1.3.181
module load R/3.4.0-X11-20160819

Rscript --no-save sim-accre.R ${SLURM_ARRAY_TASK_ID}
