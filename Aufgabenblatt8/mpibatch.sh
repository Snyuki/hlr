#!/bin/sh

#SBATCH --job-name=aintnothing
#SBATCH --output=SlurmOut/mpi-%j
#SBATCH -p west

echo "Starting Gauss-Seidel\n"

# Gauss Seidel all processes
mpirun -np 1 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 2 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 3 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 4 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 5 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 6 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 7 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 8 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 9 ./partdiff 12 1 2 2 2 1024
echo "\n"
mpirun -np 10 ./partdiff 12 1 2 2 2 1024


echo "\n\nStarting Jacobi\n"

# Jacobi all processes
echo "\n\n1 Processes\n"
mpirun -np 1 ./partdiff 12 2 2 2 2 1024
echo "\n\n2 Processes\n"
mpirun -np 2 ./partdiff 12 2 2 2 2 1024
echo "\n\n3 Processes\n"
mpirun -np 3 ./partdiff 12 2 2 2 2 1024
echo "\n\n4 Processes\n"
mpirun -np 4 ./partdiff 12 2 2 2 2 1024
echo "\n\n5 Processes\n"
mpirun -np 5 ./partdiff 12 2 2 2 2 1024
echo "\n\n6 Processes\n"
mpirun -np 6 ./partdiff 12 2 2 2 2 1024
echo "\n\n7 Processes\n"
mpirun -np 7 ./partdiff 12 2 2 2 2 1024
echo "\n\n8 Processes\n"
mpirun -np 8 ./partdiff 12 2 2 2 2 1024
echo "\n\n9 Processes\n"
mpirun -np 9 ./partdiff 12 2 2 2 2 1024
echo "\n\n10 Processes\n"
mpirun -np 10 ./partdiff 12 2 2 2 2 1024

echo "\n\nStarting Jacobi with FUNC_F0\n"

# Jacobi all processes
echo "\n\n1 Processes\n"
mpirun -np 1 ./partdiff 12 2 2 1 2 1024
echo "\n\n2 Processes\n"
mpirun -np 2 ./partdiff 12 2 2 1 2 1024
echo "\n\n3 Processes\n"
mpirun -np 3 ./partdiff 12 2 2 1 2 1024
echo "\n\n4 Processes\n"
mpirun -np 4 ./partdiff 12 2 2 1 2 1024
echo "\n\n5 Processes\n"
mpirun -np 5 ./partdiff 12 2 2 1 2 1024
echo "\n\n6 Processes\n"
mpirun -np 6 ./partdiff 12 2 2 1 2 1024
echo "\n\n7 Processes\n"
mpirun -np 7 ./partdiff 12 2 2 1 2 1024
echo "\n\n8 Processes\n"
mpirun -np 8 ./partdiff 12 2 2 1 2 1024
echo "\n\n9 Processes\n"
mpirun -np 9 ./partdiff 12 2 2 1 2 1024
echo "\n\n10 Processes\n"
mpirun -np 10 ./partdiff 12 2 2 1 2 1024

echo "\n\nStarting Jacobi complex \n"

# Jacobi complex x processes
mpirun -np 24 ./partdiff 12 2 512 2 2 1024