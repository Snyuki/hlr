#!/bin/sh

#SBATCH --job-name=aintnothing
#SBATCH --output=SlurmOut/mpi-24-2-%j
#SBATCH -p west
#SBATCH --nodes=2

mpirun -np 24 ./partdiff 12 2 512 2 2 1024
# srun ./partdiff 12 2 512 2 2 1024