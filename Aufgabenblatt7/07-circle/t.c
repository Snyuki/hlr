#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// tmp for sleep
#include <unistd.h>

#include <mpi.h>


void
print_array(int* array, char* name, int size, int rank)
{
	printf("[%i] %s = ", rank, name);
	for (int i = 0; i < size; i++)
	{
		printf("%d, ", array[i]);
	}
	printf("\n");
}

int
main (int argc, char** argv)
{
    int rank, size;
	// Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int** array = (int**)malloc(sizeof(int) * 5);
    array[1][5];

    if (rank == 0)
    {
        for (int i = 0; i < 5; i++)
        {
            array[0][i] = i;
        }
    }
    else if (rank == 2)
    {
        array[0][0] = 0;
        for (int i = 1; i < 5; i++)
        {
            array[0][i] = i + 100;
        }
    }
    else
    {
        for (int i = 0; i < 5; i++)
        {
            array[0][i] = i + (5 * rank);
        }
    }

    int first_element = -1;
    int finished = 1;
    if (rank == 0)
    {
        MPI_Send(&array[0][0], 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
    }
    if (rank == size - 1)
    {
        MPI_Recv(&first_element, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    int source = (rank - 1 + size) % size;
    int destination = (rank + 1) % size;

    for (int i = 0; i < size; i++)
    {
        if (rank == 0)
        {
            MPI_Ssend(array[0], 5, MPI_INT, destination, 0, MPI_COMM_WORLD);
            MPI_Recv(array[0], 5, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            int* tmp_array = (int*)calloc(5, sizeof(int));
            // tmp_array[5];
            for (int i = 0; i < 5; i++)
            {
                tmp_array[i] = array[0][i];
            }
            MPI_Recv(array[0], 5, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Ssend(tmp_array, 5, MPI_INT, destination, 0, MPI_COMM_WORLD);
            free(tmp_array);
        }

        if (rank == size - 1)
        {
            finished = (array[0][0] != first_element) ? 1 : 0;
        }
        MPI_Bcast(&finished, 1, MPI_INT, size - 1, MPI_COMM_WORLD);
        if (finished == 0)
        {
            break;
        }
    }

    print_array(array[0], "array[0]", 5, rank);

    free(array);

    MPI_Finalize();
    return 0;
}
