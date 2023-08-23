#!/bin/bash

#SBATCH -c 125
#SBATCH -t 4-0:00                  # wall time (D-HH:MM)
#SBATCH -p general -q public
#SBATCH --mem=0
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#ulimit -l unlimited

#module load openmpi-4.1.3-gcc-11.2.0
time srun -n 125 ./EW_dumbbell_relax.out 00 $1 $2
#module unload openmpi-4.1.3-gcc-11.2.0


#Above openmpi module and sbatch options are those that 
#work on ASU SOL but have been tested for
#intel mpi-fortran distributions like parallel-studio and one-api
#For intel mpi-fortran distributions use distribution specific mpirun instead of srun
