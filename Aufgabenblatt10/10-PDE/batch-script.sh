#!/bin/sh

#SBATCH --job-name=10/11
#SBATCH -p west


# Do in python script call --output=SlurmOut/mpi-%j


# check if the output directory exists
if [ ! -d "SlurmOut" ]; then
    mkdir SlurmOut
fi

# check if there are exactly 7 arguments
if [ $# -ne 7 ]; then
    echo "Invalid number of arguments"
    echo "Usage: $0 <number of processes> <number of threads> <method> <interlines> <funciton> <end condition> <term it/prec>"
    exit 1
fi

# print the call with the given arguments
echo "running: \"mpirun -np $1 ./partdiff $2 $3 $4 $5 $6 $7\" ..."
echo ""

# call partdiff with the given arguments
mpirun -np $1 ./partdiff $2 $3 $4 $5 $6 $7

echo ""