#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// tmp for sleep
#include <unistd.h>

#include <mpi.h>


/**
 * Return value for circle function
*/
typedef struct circle_return
{
	int* buf;						// Buffer for the section of this process
	int iterations;					// Number of iterations done until completed
	int last_element_placeholder;	// If 1 the last element is a placeholder, else its not
	int end_condition_element; 		// The element that represents the end condition (only needed by the first and the last process)
} circle_return_t;

/**
 * Print an array
 * 
 * @param array The array to print
 * @param name The name of the array (displayed in the output)
 * @param size The size of the array
 * @param rank The rank of the process
*/
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

/**
 * Initializes the array for this process and allocates the correct amount of memory.
 * If there is not a perfect distribution of the elements, the arrays of the last x processes
 * that dont hold a extra element are equalized in size (appended a '-1').
 * 
 * @param bufsize The minimum size the array of each process needs to have (N / size)
 * @param true_bufsize The actual size the array of this process needs to have (N / size or N / size + 1)
 * @param perfect_distribution If 1 the distribution between elements and processes is perfect, else its not
 * @param rank The rank of the process
*/
int*
init (int bufsize, int true_bufsize, int perfect_distribution, int rank)
{
	int max_elements_in_buf = (perfect_distribution) ? bufsize : bufsize + 1;			// necessary size to hold each array of each process
	int* buf = (int*)malloc(sizeof(int) * (max_elements_in_buf));						// allocate memory for the array of this process

	srand(time(NULL) + rank);

	for (int i = 0; i < true_bufsize; i++)
	{
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	}

	// if the distribution is not perfect and this process holds a extra element, add a placeholder
	if (!perfect_distribution && true_bufsize == bufsize)
	{
		buf[bufsize] = -1;
	}

	return buf;
}

circle_return_t*
circle (int* buf, int rank, int size, int bufsize)
{
	int package_size = bufsize * sizeof(int);					// size of the package to send

	int destination = (rank + 1) % size;						// rank of the process to send the package to
	int source = (rank - 1 + size) % size;						// rank of the process to receive the package from

	int first_element = -1;										// first element of the array of the first process (used to check if the circle is completed)
    int finished = 0;											// if 1 the algorithm is finished, else its not
    
	// Tell the last process what the first element of the first process is (end_condition) (communication between last and first process)
	if (rank == 0)
    {
		first_element = buf[0];
        MPI_Send(&buf[0], 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
    }
    if (rank == size - 1)
    {
        MPI_Recv(&first_element, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


	// Start by sending the buffer off to the next process if the rank is 0
	// Else start by receiving the buffer from the previous process to kick of the cycle
	// and prevent deadlock.
	// Do this until either the end_condition is met or the circle is completed.
	int i;
	for (i = 0; i < size; i++)
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

		// Check if the circle is completed (only done by the last process)
        if (rank == size - 1)
        {
            finished = (buf[0] != first_element) ? 0 : 1;
        }

		// The last process notifies the other processes if the circle is completed and all break if so
        MPI_Bcast(&finished, 1, MPI_INT, size - 1, MPI_COMM_WORLD);
        if (finished == 1)
        {
			i++;
            break;
        }
    }
	
	circle_return_t* ret = malloc(sizeof(circle_return_t));
	ret->buf = buf;
	ret->iterations = i;
	ret->last_element_placeholder = (buf[bufsize] == -1) ? 1 : 0;
	ret->end_condition_element = first_element;

	return ret;
}

int
main (int argc, char** argv)
{
	// Initialize MPI
	MPI_Init(&argc, &argv);

	int N;
	int rank;
	int* buf;
	int iterations;
	int end_condition_element;
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

	// If too few processes are used or too many processes for elements exist, exit
	if (size < 2 || size > N)
	{
		fprintf(stderr, "Too few or too many processes! (%i processes for %i elements)\n", size, N);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// calculate the buffer size
	int bufsize = (int) N / size;
	int true_bufsize = (rank < (N % size)) ? bufsize + 1 : bufsize;
	int perfect_distribution = (N % size == 0) ? 1 : 0;

	buf = init(bufsize, true_bufsize, perfect_distribution, rank);

	/* ----------------------- Print the initial array ----------------------- */

	int ok = 1;
    if(rank != 0) 
	{
        MPI_Recv(&ok, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    printf("\nBEFORE\n");
    for (int i = 0; i < true_bufsize; i++)
    {
        printf("rank %d: %d\n", rank, buf[i]);
    }
    if(rank != size-1) 
	{
        MPI_Ssend(&ok, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }

	/* ----------------------- Print the initial array ----------------------- */

	circle_return_t* ret = circle(buf, rank, size, bufsize);
	buf = ret->buf;																						// get the new buffer
	iterations = ret->iterations;																		// number of iterations the circle took
	true_bufsize = (ret->last_element_placeholder || perfect_distribution) ? bufsize : bufsize + 1;		// Recalculate the true_bufsize for this process after circle (dont want to print '-1')
	end_condition_element = ret->end_condition_element;													// get the end_condition_element
	free(ret);


	// flush stdout
	fflush(stdout);

	// sleep 0.5 sec to synchronize all Processes better
	sleep(0.5);

	// Wait for all Processes to finish the circles
	MPI_Barrier(MPI_COMM_WORLD);


	/* ----------------------- Print the final array ----------------------- */

	if(rank != 0)
	{
        MPI_Recv(&ok, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
	printf("\nAFTER\n");
    for (int i = 0; i < true_bufsize; i++)
    {
        printf("rank %d: %d\n", rank, buf[i]);
    }
    if(rank != size-1)
	{
        MPI_Ssend(&ok, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }

	/* ----------------------- Print the final array ----------------------- */

	// free the buffer
	free(buf);

	// Sleep to synchronize the output better (Barrier does not work properly because too fast)
	sleep(1);
	// Synchronize all Processes to end simultaneously
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Print the number of iterations the circle took
	if (rank == 0)
	{
		printf("\nIterations done: %d (%i skipped), Abbruchwert (condition_element): %i\n", iterations, size - 1 - iterations, end_condition_element);
	}

	// Finalize MPI
	MPI_Finalize();
	return EXIT_SUCCESS;
}
