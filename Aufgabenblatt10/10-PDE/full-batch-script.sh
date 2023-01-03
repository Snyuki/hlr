#!/bin/sh

#SBATCH --job-name=running-test
#SBATCH -p west


# check if the output directory exists
if [ ! -d "SlurmOut" ]; then
    mkdir SlurmOut
fi

# run partdiff with all all possible argument combinations:
# 1. number of processes (1-24)
# 2. number of threads (1-24)
# 3. method (1) -> only Gauss-Seidel
# 4. interlines (1-512) in powers of 2
# 5. function (1)
# 6. end condition (1-2)
# 7. term iterations (1-1024) in powers of 2 or term precision (1e-4 - 1e-20)

CounterIt=0
CounterPrec=0
TotalTime=0

process_tries=($(seq 1 1 12))
Interlines=(256 512)
term_iterations=(256 512 1024)
# term_precisions=(1e-12 1e-16 1e-20)
term_precisions=(7e-4)
thread_tries=(1)

# if any argument equals '--sequential', or '-s', then run the sequential reference
if [[ " ${@} " =~ " --sequential " ]] || [[ " ${@} " =~ " -s " ]]; then
    echo "Running sequential version"
    for interline in ${Interlines[@]}
    do
        for term_condition in ${term_iterations[@]}
        do
            echo "srun ./backup/partdiff 1 1 $interline 2 2 $term_condition"
            CounterIt=$((CounterIt+1))
            srun ./backup/partdiff 1 1 $interline 2 2 $term_condition
        done
    done

    for interline in ${Interlines[@]}
    do
        for term_condition in ${term_precisions[@]}
        do
            echo "srun ./backup/partdiff 1 1 $interline 2 1 $term_condition"
            CounterIt=$((CounterPrec+1))
            srun ./backup/partdiff 1 1 $interline 2 1 $term_condition
        done
    done

    echo "Iteration runs: $CounterIt"
    echo "Precision runs: $CounterPrec"
    echo "Total runs: $((CounterIt+CounterPrec))"

    exit 0
fi

# if any argument equals '--use-threads', or '-t', then try all thread numbers from 1 to 24
if [[ " ${@} " =~ " --use-threads " ]] || [[ " ${@} " =~ " -t " ]]; then
    thread_tries=($(seq 1 1 24)) 
fi

# if any argument equals '--use-all-processes' or '-ap', then try all process numbers from 1 to 24
if [[ " ${@} " =~ " --use-all-processes " ]] || [[ " ${@} " =~ " -ap " ]]; then
    process_tries=($(seq 1 1 24))
fi

# if any argument equals '--use-reduced-processes' or '-rp', then try only process numbers 1 and from 2 to 12 in steps of 2
if [[ " ${@} " =~ " --use-reduced-processes " ]] || [[ " ${@} " =~ " -rp " ]]; then
    process_tries=(1 2 4 6 8 10 12) 
fi


# Do for termination iterations
for processes in ${process_tries[@]}
do
    for threads in ${thread_tries[@]}
    do
        for interline in ${Interlines[@]}
        do
            for term_condition in ${term_iterations[@]}
            do
                echo "mpirun -np $processes ./partdiff $threads 1 $interline 2 2 $term_condition"
                CounterIt=$((CounterIt+1))
                # TotalTime=$(( TotalTime+(400/($processes*$threads)/4) ))
                mpirun -np $processes ./partdiff $threads 1 $interline 2 2 $term_condition
            done
        done
    done
done

# Do for termination precision
for processes in ${process_tries[@]}
do
    for threads in ${thread_tries[@]}
    do
        for interline in ${Interlines[@]}
        do
            for term_condition in ${term_precisions[@]}
            do
                echo "mpirun -np $processes ./partdiff $threads 1 $interline 2 1 $term_condition"
                CounterPrec=$((CounterPrec+1))
                # TotalTime=$(( TotalTime+(400/($processes*$threads)/4) ))
                mpirun -np $processes ./partdiff $threads 1 $interline 2 1 $term_condition
            done
        done
    done
done

echo "Iteration runs: $CounterIt"
echo "Precision runs: $CounterPrec"
echo "Total runs: $((CounterIt+CounterPrec))"
# echo "Total time: $TotalTime"
