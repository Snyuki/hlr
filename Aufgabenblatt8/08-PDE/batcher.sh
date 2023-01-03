#!/bin/bash
#SBATCH -p west
#SBATCH --job-name=butwhy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=parallel-%j

mpirun -np 8 ./partdiff 12 2 512 2 2 256
