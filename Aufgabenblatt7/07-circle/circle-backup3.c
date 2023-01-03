#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// tmp for sleep
#include <unistd.h>

#include <mpi.h>

void
print_array(int* array, char* name, int size, int rank)
{
	printf("\n");
	printf("\n[%i] %s = ", rank, name);
	for (int i = 0; i < size; i++)
	{
		printf("%d, ", array[i]);
	}
	printf("\n");
}

int*
init (int bufsize, int true_bufsize, int perfect_distribution, int rank)
{
	int max_elements_in_buf = (perfect_distribution) ? bufsize : bufsize + 1;
	int* buf = (int*)malloc(sizeof(int) * (max_elements_in_buf));

	srand(time(NULL) + rank);

	for (int i = 0; i < true_bufsize; i++)
	{
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	}
	if (!perfect_distribution && true_bufsize == bufsize)
	{
		buf[bufsize] = -1;
	}

	return buf;
}

int*
circle (int* buf, int rank, int size, int bufsize)
{
	int package_size = bufsize * sizeof(int);

	int destination = (rank + 1) % size;
	int source = (rank - 1 + size) % size;

	int first_element = -1;
    int finished = 0;
    if (rank == 0)
    {
        MPI_Send(&buf[0], 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
    }
    if (rank == size - 1)
    {
        MPI_Recv(&first_element, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


	// Start by sending the buffer off to the next process if the rank is 0
	// Else start by receiving the buffer from the previous process to kick of the cycle
	// and prevent deadlock.

	for (int i = 0; i < size; i++)
    {
        if (rank == 0)
        {
            MPI_Ssend(buf, package_size, MPI_INT, destination, 0, MPI_COMM_WORLD);
            MPI_Recv(buf, package_size, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            int* tmp_buf = (int*)calloc(package_size, sizeof(int));
            for (int i = 0; i < package_size; i++)
            {
                tmp_buf[i] = buf[i];
            }
            MPI_Recv(buf, package_size, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Ssend(tmp_buf, package_size, MPI_INT, destination, 0, MPI_COMM_WORLD);
            free(tmp_buf);
        }

        if (rank == size - 1)
        {
            finished = (buf[0] != first_element) ? 0 : 1;
        }
        MPI_Bcast(&finished, 1, MPI_INT, size - 1, MPI_COMM_WORLD);
        if (finished == 1)
        {
            break;
        }
    }
	
	return buf;
}

int
main (int argc, char** argv)
{
	// Initialize MPI
	MPI_Init(&argc, &argv);

	int N;
	int rank;
	int* buf;
	int size;
	
	// Get the number of processes and rank
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc < 2)
	{
		printf("Arguments error!\nPlease specify a buffer size.\n");
		return EXIT_FAILURE;
	}
	
	// Array length
	N = atoi(argv[1]);

	// calculate the buffer size
	int bufsize = (int) N / size;
	int true_bufsize = (rank < (N % size)) ? bufsize + 1 : bufsize;
	int perfect_distribution = (N % size == 0) ? 1 : 0;

	int to_display =  (perfect_distribution) ? bufsize : bufsize + 1;
	
	buf = init(bufsize, true_bufsize, perfect_distribution, rank);

	printf("\nBEFORE\n");

	for (int i = 0; i < to_display; i++)
	{
		printf("rank %d: %d\n", rank, buf[i]);
	}

	buf = circle(buf, rank, size, bufsize);

	// flush stdout
	fflush(stdout);

	// sleep 1 sec to synchronize all Processes better
	sleep(1);

	// Wait for all Processes to finish the circles
	MPI_Barrier(MPI_COMM_WORLD);

	printf("\nAFTER\n");

	for (int j = 0; j < to_display; j++)
	{
		printf("rank %d: %d\n", rank, buf[j]);
	}

	// free the buffer
	free(buf);

	// Finalize MPI
	MPI_Finalize();
	return EXIT_SUCCESS;
}
