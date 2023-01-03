#!/bin/bash

## do not run on login node!
# -> this script doesn't use srun, but will instead simply run on the current node.

## this example of a test-script is incomplete. You will need to add useful test-cases yourself.
# change it at your discression and make it prettier!

## the general idea:
# ssh over to an allocated west-node and run the script there.
# the difference of the tested program's and the reference program's output is stored in diff.txt
# result.txt holds an elaborate log of each run


## some syntax hints:
# -  ' ' and '\n' are (sometimes) very important in bashscript syntax and easily lead to errors :/
# -  don't forget to reference variables with the '$' prefix
# -  a string like 'this' will be interpreted as-is. a string like "this" will evaluate the variables inside.


# names of compiled files that are ready to be executed by mpirun
experimental_program='partdiff'
reference_program='partdiff-ref'

# these names could also be passed in by the 1. and 2. parameter when calling ./test.sh:
#experimental_program=$1
#reference_program=$2

# remove old results
rm result.txt
rm diff.txt

# this function could be nicer..
# temporary files like actual.txt and expected.txt might not be necessary for the 'diff' command.
call_func(){
    mpirun -n $7 ./$experimental_program $1 $2 $3 $4 $5 $6 > actual.txt
    mpirun -n 1 ./$reference_program $1 $2 $3 $4 $5 $6 > expected.txt

    # write out everything that happens; both outputs and their difference:
    echo '###############   test run    ########################################' >> result.txt
    echo "args: $1 $2 $3 $4 $5 $6  procnum: $7" >> result.txt
    cat actual.txt >> result.txt
    echo '############### reference run ########################################' >> result.txt
    cat expected.txt >> result.txt
    echo '###############      diff     ########################################' >> result.txt
    diff actual.txt expected.txt >> result.txt

    # write out only the difference to a seperate file, to be able to quickly look it over:
    echo "################################ args: $1 $2 $3 $4 $5 $6  procnum: $7 :"  >> diff.txt
    diff actual.txt expected.txt >> diff.txt
}

# example on how to call the function with different parameters:

term_iter=2
term_prec=1

method=2

thread_list=(1 2)

# for-each loop over thread_list
for threads in ${thread_list[@]}
do
    call_func 1 $method 0 1 $term_prec 1e-5 $threads
done

# remove temporary files
rm actual.txt
rm expected.txt

echo 'done'
