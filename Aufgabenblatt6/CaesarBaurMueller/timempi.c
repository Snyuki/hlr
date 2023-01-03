#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

// Include header files
#include <mpi.h>



int main(int argc, char *argv[]) {

    struct timeval tv;
    time_t time;
    int micro_sec;
    char time_string[30];
    char output[80];
    char hostname[30];
    int rank, size;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Arrays for Task 3
    // int onedarray[24] = {3, 4, 2, 3, 0, -3, 9, 11, 23, 12, 23, 2, 13, 4, 56, 3, 5, 9, 3, 5, 3, 1, 4, 9};
    // int twodarray[6][4] = {{3, 4, 2, 3}, {0, -3, 9, 11}, {23, 12, 23, 2}, {13, 4, 56, 3}, {5, 9, 3, 5}, {3, 1, 4, 9}};
    // int threedarray[2][3][4] = {{{3, 4, 2, 3}, {0, -3, 9, 11}, {23, 12, 23, 2}}, {{13, 4, 56, 3}, {5, 9, 3, 5}, {3, 1, 4, 9}}};

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // if only one process is existing, print the same message a slave would print
    if (size == 1)
    {
        gethostname(hostname, 30);
        gettimeofday(&tv, NULL);

        time = tv.tv_sec;
        micro_sec = tv.tv_usec;

        strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
        snprintf(output, 80, "[%i] %s : %s.%d\n", rank, hostname, time_string, (int)micro_sec);

        printf("%s", output);
        printf("[%i] beendet jetzt!\n", rank);
        
        MPI_Finalize();

        return 0; 
    }


    // If master process
    if (rank == (size - 1))
    {
        // receive output from all other processes
        for (int i = 0; i < size - 1; i++)
        {
            MPI_Recv(output, 80, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("%s", output);
        }

        // collect microsecond time from all other processes with collective operation
        int micro_sec_array[size - 1];
        MPI_Gather(&micro_sec, 1, MPI_INT, micro_sec_array, 1, MPI_INT, size - 1, MPI_COMM_WORLD);

        // calculate the minimum microsecond time and the maximum difference between two microsecond times of all threads
        int min_micro_sec = micro_sec_array[0];
        int max_micro_sec = micro_sec_array[0];
        int max_diff = 0;
        for (int i = 0; i < size - 1; i++)
        {
            if (micro_sec_array[i] < min_micro_sec)
            {
                min_micro_sec = micro_sec_array[i];
            }

            if (micro_sec_array[i] > max_micro_sec)
            {
                max_micro_sec = micro_sec_array[i];
            }
        }
        max_diff = max_micro_sec - min_micro_sec;

        // print the minimum microsecond time and the maximum difference between two microsecond times
        printf("[%i] Kleinster MS-Differenz: %d\n", rank, min_micro_sec);
        printf("[%i] Größte Differenz: %d\n", rank, max_diff);      

        fflush(stdout); // Flush output to make sure all output is printed before the sleep
        // Sleep short time to make sure there aren't any abnormalities in the output
        // sleep(1);
        // Allow the other processes to pass the barrier
        MPI_Barrier(MPI_COMM_WORLD);
        printf("[%i] beendet jetzt!\n", rank);
    }
    // Else slave process
    else
    {
        gethostname(hostname, 30);
        gettimeofday(&tv, NULL);

        time = tv.tv_sec;
        micro_sec = tv.tv_usec;

        strftime(time_string, 30, "%Y-%m-%d %T", localtime(&time));
        snprintf(output, 80, "[%i] %s : %s.%d\n", rank, hostname, time_string, (int)micro_sec);
        // printf("%d\n", (int)micro_sec);

        // send output to master process
        MPI_Send(output, 80, MPI_CHAR, size - 1, 0, MPI_COMM_WORLD);

        // Send microsecond time to master process with collective operation
        MPI_Gather(&micro_sec, 1, MPI_INT, NULL, 0, MPI_INT, size - 1, MPI_COMM_WORLD);

        // Barrier until all processes are done
        MPI_Barrier(MPI_COMM_WORLD);
        printf("[%i] beendet jetzt!\n", rank);
    }

    MPI_Finalize(); 

    return 0;
}