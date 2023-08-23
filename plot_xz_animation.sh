#!/bin/bash

#SBATCH -c 1
#SBATCH -t 1-00:00                  # wall time (D-HH:MM)
#SBATCH -q public -p general
#SBATCH --mem=100G
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#ulimit -l unlimited

python xz_indiv_animations.py E $1 0
python xz_indiv_animations.py E $1 1
python xz_indiv_animations.py E $1 2
python xz_indiv_animations.py E $1 3
python xz_indiv_animations.py E $1 4
python xz_indiv_animations.py E $1 5
python xz_indiv_animations.py E $1 6
python xz_indiv_animations.py E $1 7
python xz_indiv_animations.py E $1 8
python xz_indiv_animations.py E $1 9
python xz_indiv_animations.py E $1 10


