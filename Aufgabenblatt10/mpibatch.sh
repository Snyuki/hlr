#!/bin/sh

#SBATCH --job-name=deadlock.exe
#SBATCH --output=SlurmOut/mpi-%j
#SBATCH -p west

# echo "Starting Gauss-Seidel with Term_Iteration\n"

# Gauss Seidel all processes
# mpirun -np 1 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 2 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 3 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 4 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 5 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 6 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 7 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 8 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 9 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 10 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 11 ./partdiff 12 1 256 2 2 1024
# echo "\n"
# mpirun -np 12 ./partdiff 12 1 256 2 2 1024
# echo "\n"


echo "Starting Gauss-Seidel with Term_Precision 4 Processes \n"
mpirun -np 4 ./partdiff 1 1 256 2 1 3e-4
echo "\n"
mpirun -np 4 ./partdiff 1 1 256 2 1 7e-4
echo "\n"
mpirun -np 4 ./partdiff 1 1 256 2 1 3e-5
echo "\n"
mpirun -np 4 ./partdiff 1 1 256 2 1 7e-5
echo "\n"
mpirun -np 4 ./partdiff 1 1 256 2 1 3e-6
echo "\n"
mpirun -np 4 ./partdiff 1 1 256 2 1 7e-6


echo "Starting Gauss-Seidel with Term_Precision all Processes \n"
mpirun -np 1 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 2 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 3 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 4 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 5 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 6 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 7 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 8 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 9 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 10 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 11 ./partdiff 1 1 256 2 1 7e-6
echo "\n"
mpirun -np 12 ./partdiff 1 1 256 2 1 7e-6

# # Jacobi complex function
# echo "\n\nStarting Jacobi\n"

# mpirun -np 1 ./partdiff 12 2 2 2 2 1024

# echo "\n\nStarting Jacobi with FUNC_F0\n"

# # Jacobi Function F0
# echo "\n\n \n"
# mpirun -np 1 ./partdiff 12 2 2 1 2 1024

# echo "\n\nStarting Jacobi complex \n"

# # Jacobi complex x processes
# mpirun -np 24 ./partdiff 12 2 512 2 2 1024