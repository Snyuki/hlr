#!/bin/sh

#SBATCH --job-name=mpi
#SBATCH --output=SlurmOut/mpi-%j
#SBATCH --ntasks-per-node=24
#SBATCH -p west
#SBATCH --nodes=7


mpirun -np 168 ./timempi