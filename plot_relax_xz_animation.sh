#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-15:00                  # wall time (D-HH:MM)
#SBATCH -q public -p highmem
#SBATCH --mem=200G
#SBATCH -A tpatel28            # Account to pull cpu hours from (commented out)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#ulimit -l unlimited

python proj_xz_plot.py rpot


