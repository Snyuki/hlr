#!/bin/sh

#SBATCH --job-name=aintnothing
#SBATCH --output=SlurmOut/mpi-%j
#SBATCH -p west

mpirun -np 1 ./partdiff 1 2 256 2 1 1e-12
# mpirun -np 1 ./partdiff 12 1 2  2 2 1024