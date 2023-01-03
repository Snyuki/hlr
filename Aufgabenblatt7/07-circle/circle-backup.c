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
init (int bufsize, int rank)
{
	int* buf = (int*)malloc(sizeof(int) * bufsize);

	srand(time(NULL) + rank);

	for (int i = 0; i < bufsize; i++)
	{
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	}

	return buf;
}

int*
circle (int* buf, int rank, int size, int bufsize)
{
	int package_size = bufsize * sizeof(int);

	int* recvbuf = (int*)calloc(bufsize, sizeof(int));

	int destination = (rank + 1) % size;
	int source = (rank - 1 + size) % size;

	// Start by sending the buffer off to the next process if the rank is 0
	// Else start by receiving the buffer from the previous process to kick of the cycle
	// and prevent deadlock.
	if (rank == 0)
	{
		MPI_Ssend(buf, package_size, MPI_INT, destination, 0, MPI_COMM_WORLD);
		MPI_Recv(recvbuf, package_size, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else
	{
		MPI_Recv(recvbuf, package_size, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Ssend(buf, package_size, MPI_INT, destination, 0, MPI_COMM_WORLD);
		
	}

	for (int i = 0; i < bufsize; i++)
	{
		buf[i] = recvbuf[i];
	}
	free(recvbuf);

	// print_array(recvbuf, "recvbuf", bufsize, rank);
	// print_array(buf, "buf", bufsize, rank);
	
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
	// if (rank < (N % size)) bufsize++;
	bufsize++;
	
	buf = init(bufsize, rank);

	printf("\nBEFORE\n");

	for (int i = 0; i < bufsize; i++)
	{
		printf("rank %d: %d\n", rank, buf[i]);
	}

	// Send first element of first process to last process
	if (rank == 0)
	{
		MPI_Send(&buf[0], 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
	}
	
	// Only important for the last process
	int first_element = -1;

	// Receive first element of first process from first process
	if (rank == (size - 1))
	{
		MPI_Recv(&first_element, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	// int finished = 1;

	// while (finished)
	// {
	// 	buf = circle(buf, rank, size, bufsize);

	// 	if (rank == size - 1)
	// 	{
	// 		int finished = (first_element != buf[0]) ? 1 : 0;

	// 		// Send to all if finished
	// 		MPI_Bcast(&finished, 1, MPI_INT, size - 1, MPI_COMM_WORLD);
	// 	}
	// 	else
	// 	{
	// 		MPI_Bcast(&finished, 1, MPI_INT, size - 1, MPI_COMM_WORLD);
	// 	}

	// 	MPI_Barrier(MPI_COMM_WORLD);
		
	// }
	

	// Some kind of while loop with end condition
	// -------
	buf = circle(buf, rank, size, bufsize);
	buf = circle(buf, rank, size, bufsize);
	// -------

	// flush stdout
	fflush(stdout);

	// sleep 1 sec to synchronize all Processes better
	sleep(1);

	// Wait for all Processes to finish the circles
	MPI_Barrier(MPI_COMM_WORLD);

	printf("\nAFTER\n");

	for (int j = 0; j < bufsize; j++)
	{
		printf("rank %d: %d\n", rank, buf[j]);
	}

	// free the buffer
	free(buf);

	// Finalize MPI
	MPI_Finalize();
	return EXIT_SUCCESS;
}
