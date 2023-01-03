#!/bin/sh

#SBATCH --job-name=aintnothing-leit
#SBATCH --output=SlurmOut/mpi-hybrid-Leistungsanalyse-%j
#SBATCH -p west
#SBATCH --nodes=3

# echo "Starting Gauss-Seidel\n"

# Gauss Seidel all processes
# mpirun -np 1 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 2 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 3 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 4 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 5 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 6 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 7 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 8 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 9 ./partdiff-par-hybrid 12 1 512 2 2 1024
# echo "\n"
# mpirun -np 10 ./partdiff-par-hybrid 12 1 512 2 2 1024


# echo "\n\nStarting Jacobi\n"

# Jacobi all processes
# echo "\n\n1 Processes\n"
# mpirun -np 1 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n2 Processes\n"
# mpirun -np 2 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n3 Processes\n"
# mpirun -np 3 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n4 Processes\n"
# mpirun -np 4 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n5 Processes\n"
# mpirun -np 5 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n6 Processes\n"
# mpirun -np 6 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n7 Processes\n"
# mpirun -np 7 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n8 Processes\n"
# mpirun -np 8 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n9 Processes\n"
# mpirun -np 9 ./partdiff-par-hybrid 12 2 2 2 2 1024
# echo "\n\n10 Processes\n"
# mpirun -np 10 ./partdiff-par-hybrid 12 2 2 2 2 1024

# echo "\n\nStarting Jacobi with FUNC_F0\n"

# # Jacobi all processes
# echo "\n\n1 Processes\n"
# mpirun -np 1 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n2 Processes\n"
# mpirun -np 2 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n3 Processes\n"
# mpirun -np 3 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n4 Processes\n"
# mpirun -np 4 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n5 Processes\n"
# mpirun -np 5 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n6 Processes\n"
# mpirun -np 6 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n7 Processes\n"
# mpirun -np 7 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n8 Processes\n"
# mpirun -np 8 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n9 Processes\n"
# mpirun -np 9 ./partdiff-par-hybrid 12 2 2 1 2 1024
# echo "\n\n10 Processes\n"
# mpirun -np 10 ./partdiff-par-hybrid 12 2 2 1 2 1024

# echo "\n\nStarting Jacobi complex \n"

# # Jacobi complex x processes
# mpirun -np 24 ./partdiff-par-hybrid 12 2 512 2 2 1024


echo "\n\nStarting Leistungsanalyse Jacobi\n"
mpirun -np 12 ./partdiff-par-hybrid 1 2 512 2 2 1024
mpirun -np 24 ./partdiff-par-hybrid 1 2 512 2 2 1024
mpirun -np 1 ./partdiff-par-hybrid 12 2 512 2 2 1024
mpirun -np 1 ./partdiff-par-hybrid 24 2 512 2 2 1024
mpirun -np 2 ./partdiff-par-hybrid 6 2 512 2 2 1024
mpirun -np 2 ./partdiff-par-hybrid 12 2 512 2 2 1024
mpirun -np 12 ./partdiff-par-hybrid 2 2 512 2 2 1024