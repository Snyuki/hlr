import os
import subprocess
from time import sleep
from datetime import timedelta


# TODO Matrix is empty on all runs for some reason
# TODO parallel deadlocks on 'mpirun -np 1 ./partdiff 1 1 256 1 1 1e-12'
# TODO sequential deadlocks on 'srun ./backup/partdiff 1 1 256 1 1 1e-12'

# TODO Set output file name in args
# TODO Set sequential reference file path in args
# TODO Set mode in args (sequential, parallel, both)
# TODO Test the comparison of the results with dummy numbers
# TODO reduce redundancy of parameter listing in full-batch-script.sh and test-script.py

# plan

# run partdiff with all possible argument combinations with threads and processes (parallel).

# run partdiff with all possible argument combinations without threads and processes (sequential).
# But only if no sequential run with this specifications was done before.

# compare the results of the sequential run with the results of the parallel runs.


# number of threads and processes used (has to be updated when changed in full-batch-script.sh)
num_threads = 1
num_processes = 12

# number of lines of one output
output_length = 20

time_steps = 10
sequential_reference_path = "SlurmOut/blatt10-reference-155940"
parallel_test_path = ""


def get_time_elapsed(seconds):
    return timedelta(seconds=seconds)


def run_sequential_reference(username="mueller2"):
    print("Starting Sequential-Reference...")

    seconds_elapsed = 0

    # Run Sequential-Reference
    ret_value = subprocess.check_output("sbatch --output=SlurmOut/blatt10-reference-%j ./full-batch-script.sh -s", shell=True)
    
    # Get Job-ID
    job_id = ret_value.decode("utf-8").split(" ")[-1].strip()
    print("Job-ID: {}".format(job_id))

    # Wait for Batch-Script to finish
    is_still_running = job_id in subprocess.check_output("squeue -u {}".format(username), shell=True)
    while(is_still_running):
        sleep(time_steps)
        seconds_elapsed += time_steps
        print("Sequential-Reference is still running... (time elapsed: {})".format(get_time_elapsed(seconds_elapsed)))
        is_still_running = job_id in subprocess.check_output("squeue -u mueller2", shell=True)

    print("Sequential-Reference finished! (Total Time elapsed: {})".format(get_time_elapsed(seconds_elapsed)))


def run_batch_script(username="mueller2"):

    print("Starting Batch-Script...")

    # Run Batch-Script
    ret_value = subprocess.check_output("sbatch --output=SlurmOut/blatt10-parallel-test-%j ./full-batch-script.sh", shell=True)
    
    # Get Job-ID
    job_id = ret_value.decode("utf-8").split(" ")[-1].strip()
    print("Job-ID: {}".format(job_id))

    # Set path to output of file
    parallel_test_path = "SlurmOut/blatt10-parallel-test-{}".format(job_id)

    seconds_elapsed = 0

    # Wait for Batch-Script to finish
    is_still_running = job_id in subprocess.check_output("squeue -u {}".format(username), shell=True)
    while(is_still_running):
        sleep(time_steps)
        seconds_elapsed += time_steps
        print("Test is still running... (time elapsed: {})".format(get_time_elapsed(seconds_elapsed)))
        is_still_running = job_id in subprocess.check_output("squeue -u {}".format(username), shell=True)

    print("Batch-Script finished! (Total Time elapsed: {})".format(get_time_elapsed(seconds_elapsed)))

    return parallel_test_path


def compare_results(parallel_test_path):
    print("Comparing results...")

    # if the parallel path was not set, abort
    if parallel_test_path == "":
        print("Error: parallel_test_path is not set!")
        return

    # Compare results of sequential run with results of parallel runs
    sequential_reference = []
    parallel_test = []

    # read files
    with open(sequential_reference_path, "r") as sequential_reference_file:
        sequential_reference = sequential_reference_file.readlines()

    with open(parallel_test_path, "r") as parallel_test_file:
        parallel_test = parallel_test_file.readlines()
    
    # compare line by line
    sequential_length = len(sequential_reference)
    parallel_length = len(parallel_test)

    # calculate how often one output is repeated in the parallel output compared to the sequential output
    factor_parallel_repeat = sequential_length
    print("seq/par: ({}/{}) = {}".format(parallel_length, sequential_length, parallel_length / sequential_length))

    current_index_sequential = 0
    current_index_parallel = 0

    print("current_index_sequential_parallel: ({}/{}), parallel_seq_length: ({}/{})".format(current_index_sequential, current_index_parallel, sequential_length, parallel_length))

    while (current_index_parallel < parallel_length):
        # setback sequential_length factor_parallel_repeat - 1 times
        # print("factor_parallel_repeat: {}".format(factor_parallel_repeat))
        for repeat_count in range(factor_parallel_repeat):
            print("repeat_count: {}, factor_parallel_repeat: {}".format(repeat_count, factor_parallel_repeat))
            diff_string = ""
            for local_index in range(output_length):
                print(local_index, current_index_sequential, current_index_parallel)
                line_sequential = sequential_reference[current_index_sequential + local_index]
                line_parallel = parallel_test[current_index_parallel + local_index]

                if line_sequential.startswith('srun') or line_parallel.startswith('mpirun'):
                    continue

                if line_sequential != line_parallel:
                    diff_string += "Line seq:{}/par:{} differs: Sequential: {} Parallel: {}".format(current_index_sequential, current_index_parallel, line_sequential, line_parallel)

                # increment indicies
                current_index_sequential += 1
                current_index_parallel += 1

            # setback sequential_length by output_length
            if repeat_count < factor_parallel_repeat:
                current_index_sequential -= output_length

            # write diff to output file
            with open("diff.txt", "a") as diff_file:
                # get call with parameters
                current_sequential_call_index = current_index_sequential - output_length if (repeat_count == current_index_sequential) else current_index_sequential
                call_sequential = sequential_reference[current_sequential_call_index]
                call_parallel = parallel_test[current_index_parallel - output_length]

                diff_file.write("Sequential call: {} vs. Parallel call: {}\n".format(call_sequential, call_parallel))
                diff_file.write(diff_string)
                diff_file.write("\n\n")






    print("Comparing completed!")



if __name__ == "__main__":
    # run_sequential_reference()
    parallel_test_path = run_batch_script()
    compare_results(parallel_test_path)
