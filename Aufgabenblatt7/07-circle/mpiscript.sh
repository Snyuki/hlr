#!/bin/sh

#SBATCH --job-name=mpi-circle
#SBATCH --output=SlurmOut/mpi-%j
#SBATCH -p west

# spack load scorep

# srun make
mpirun -np 5 ./circle 13
# mpirun -np 5 ./a.out