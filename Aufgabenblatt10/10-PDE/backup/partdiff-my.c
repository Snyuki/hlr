/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <mpi.h>
#include <unistd.h>

#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     					*/
	uint64_t  num_matrices;   /* number of matrices                             					*/
	double    h;              /* length of a space between two lines            					*/
	double    ***Matrix;      /* index matrix used for addressing M             					*/
	double    *M;             /* two matrices with real values                 						*/
	uint64_t  num_rows;       /* number of rows prior to the exchange of rows for this process      */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

typedef struct start_end_results
{
	int start;				/* start index of the rows in memory for this process */
	int end;				/* end index of the rows in memory for this process   */
} start_end_results_t;


/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options, int rank, int size)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	uint64_t const N = arguments->N;
	// calculate number of rows in memory for this process
	uint64_t rows_in_memory = (N + 1) / size;
	int rest = (N + 1) % size;
	if (rank < rest)
	{
		rows_in_memory++;
	}
	arguments->num_rows = rows_in_memory;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments, int rank, int size)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;
	uint64_t const num_rows = arguments->num_rows;

	// Need to allocate more memory because we need to know the border row of the neighbor processes
	int is_border_process = (rank == 0 || rank == size - 1) ? 1 : 0;
	int rows_from_neighbor_processes = (is_border_process) ? 1 : 2;

	arguments->M = allocateMemory(arguments->num_matrices * (num_rows + rows_from_neighbor_processes) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	// printf("[%d] Allocated %ld with %ld rows and %d neighbors\n", rank, (num_rows + rows_from_neighbor_processes), num_rows, rows_from_neighbor_processes);

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((num_rows + rows_from_neighbor_processes) * sizeof(double*));

		for (j = 0; j <= num_rows; j++)
		{
			int actual_row = j;
			if (rank != 0)
			{
				actual_row += rows_from_neighbor_processes;
			}

			arguments->Matrix[i][j] = arguments->M + (i * (num_rows + rows_from_neighbor_processes) * (N + 1)) + (j * (N + 1));
		}
	}
}

static
void
allocateMatrices2 (struct calculation_arguments* arguments, int start, int end, int rank, int size)
{
    uint64_t i, j;

    uint64_t range = end-start;

    uint64_t const N = arguments->N;


    if(rank==0) {
        arguments->M = allocateMemory(arguments->num_matrices * (range+1+1) * (N + 1) * sizeof(double));
        arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));
        printf("Im first rank %d, and I contain %lu lines\n", rank, range+1+1);

        for (i = 0; i < arguments->num_matrices; i++)
        {
            arguments->Matrix[i] = allocateMemory((range+1+1) * sizeof(double*));

            for (j = 0; j <= range+1; j++)
            {
                arguments->Matrix[i][j] = arguments->M + (i * (range+1+1) * (N + 1)) + (j * (N + 1));
            }
        }
    }
    if(rank==size-1) {
        arguments->M = allocateMemory(arguments->num_matrices * (range+1+1) * (N + 1) * sizeof(double));
        arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));
        printf("Im last rank %d, and I contain %lu lines\n", rank, range+1+1);

        for (i = 0; i < arguments->num_matrices; i++)
        {
            arguments->Matrix[i] = allocateMemory((range+1+1) * sizeof(double*));

            for (j = 0; j <= range+1; j++)
            {
                arguments->Matrix[i][j] = arguments->M + (i * (range+1+1) * (N + 1)) + (j * (N + 1));
            }
        }
    }
    else if(rank!=0 && rank!=size-1) {
        arguments->M = allocateMemory(arguments->num_matrices * (range+1+2) * (N + 1) * sizeof(double));
        arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));
        printf("Im middle rank %d, and I contain %lu lines\n", rank, range+1+2);

        for (i = 0; i < arguments->num_matrices; i++)
        {
            arguments->Matrix[i] = allocateMemory((range+1+2) * sizeof(double*));

            for (j = 0; j <= range+2; j++)
            {
                arguments->Matrix[i][j] = arguments->M + (i * (range+1+2) * (N + 1)) + (j * (N + 1));
            }
        }
    }

} 

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options, int rank, int size)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;
	uint64_t num_rows = arguments->num_rows;

	// printf("Rank %d: rows_in_memory = %ld\n", rank, num_rows);

	// If this is not the first process, we need to initialize from the second [1] row
	// If this is not the last process, we need to initialize to the second last row
	if (rank != size - 1) num_rows--;


	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= num_rows; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	// printf("Rank %d: Matrix initialized\n", rank);

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		// Only initialize own lines/borders
		for(i = ((rank != 0) ? 1 : 0); i <= num_rows; i++)
		{
			for (g = 0; g <= N; g++)
			{
				for (j = 0; j < arguments->num_matrices; j++)
				{
					if (rank == 0)
					{
						Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
					}
					if (rank == size - 1)
					{
						Matrix[j][num_rows][i] = 1 - (h * i); // Untere Kante
					}
					Matrix[j][i][0] = 1 + (1 - (h * i)); // Linke Kante
					Matrix[j][num_rows - i][N] = h * i; // Rechte Kante
				}
			}
		}
	}
	// printf("Rank %d: Borders initialized\n", rank);
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, int start_row, int end_row, int rank, int size)
{
	uint64_t i;			/* local variable for loops */
	int j;           	/* local variable for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;
	uint64_t const num_rows = arguments->num_rows;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	int range = end_row - start_row;
	int bufferlines = (rank == 0 || rank == size - 1) ? 1 : 2;

	// printf("[%d] Starting calculation with %ld iterations\n", rank, options->term_iteration);

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		/* exchange borders with neighbor processes */
		if (rank == 0)
		{
			// Send last row to next process
			MPI_Ssend(Matrix_In[num_rows - 1], N + 1, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
			// Receive first row from next process
			MPI_Recv(Matrix_In[num_rows], N + 1, MPI_DOUBLE, rank + 1, rank * 777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else if (rank != size - 1)
		{
			// Receive first row from previous process
			MPI_Recv(Matrix_In[0], N + 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// Send first row to previous process
			MPI_Ssend(Matrix_In[1], N + 1, MPI_DOUBLE, rank - 1, (rank - 1) * 777, MPI_COMM_WORLD);
			// Send last row to next process
			MPI_Ssend(Matrix_In[num_rows - 1], N + 1, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD);
			// Receive first row from next process
			MPI_Recv(Matrix_In[num_rows], N + 1, MPI_DOUBLE, rank + 1, rank * 777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			// Receive first row from previous process
			MPI_Recv(Matrix_In[0], N + 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// Send first row to previous process
			MPI_Ssend(Matrix_In[1], N + 1, MPI_DOUBLE, rank - 1, (rank - 1) * 777, MPI_COMM_WORLD);
		}

		maxResiduum = 0;

		/* over all rows */
		for (i = 1; i < range + bufferlines; i++)
		{
			// printf("[%d] I = %ld of %ld rows\n", rank, i, num_rows);
			double fpisin_i = 0.0;
			int realrow = start_row + i;;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)realrow);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/**
 * Display Matrix for multiple MPI Threads
*/
static
void
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int start, int end, int rank, int size)
{
    int x, y;

    double** Matrix = arguments->Matrix[results->m];

    int const interlines = options->interlines;

    // printf("Matrix: \n");

	int range = end - start;

	int add = 0;
	if (rank != 0)
	{
		add = 1;
	}

    int ok = 1;
    if(rank > 0) 
    {
		// printf("[%d] Waiting for rank %d\n", rank, rank-1);
        MPI_Recv(&ok, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    for (y = 0; y < 9; y++)
    {
        for (x = 0; x < 9; x++)
        {
            int global_index = y * (interlines + 1);
            if(global_index - start >= 0 && global_index - start <= end-start)
			{
                printf("%7.8f ", Matrix[global_index - start + add][x * (interlines + 1)]);
            }
        }
        int i_index = y * (interlines + 1);
        if(i_index - start >= 0 && i_index - start <= end-start) 
		{
            printf ("\n");
        }
    }
    if(rank != size - 1) 
    {
		// printf("[%d] Sending to rank %d\n", rank, rank+1);
        MPI_Ssend(&ok, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }

	// Print matrix off all processes
	MPI_Barrier(MPI_COMM_WORLD);
    sleep(1);
    if(rank == 0) 
	{
        printf ("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int ok2 = 1;
    if(rank != 0) 
    {
        MPI_Recv(&ok2, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(rank != 0 && rank != size-1) 
	{
        range += 1;
    }

    for (y = 0; y <= range + 1; y++)
    {
        for (x = 0; x <= (int)arguments->N; x++)
        {
            printf("%7.8f ", Matrix[y][x]);
        }
        printf ("\n");
    }
    printf ("\n");
    if(rank != size-1) 
    {
        sleep(1);
        MPI_Ssend(&ok2, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }

    fflush (stdout);
}

start_end_results_t*
calcStartEnd (int N, int rank, int size)
{
    int lines = (N+1)/size;
    int extra = (N+1)%size;

    int startline = 0;
    for (int i = 0; i < rank; i++) 
	{
        startline += lines;
        if (i < extra)
		{
            startline += 1;
        }
    }
    int endline = startline + lines - 1;
    if (rank < extra) 
	{
            endline += 1;
    }

    printf("Lines = %d, Startline = %d, Endline = %d\n", endline-startline+1, startline, endline);

	start_end_results_t* results = malloc(sizeof(start_end_results_t));
	results->start = startline;
	results->end = endline;

    return results;
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	// Init MPI Stuff
	MPI_Init(&argc, &argv);

	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	askParams(&options, argc, argv, rank);
	
	initVariables(&arguments, &results, &options, rank, size);

	// calc start and end
	start_end_results_t* start_end = calcStartEnd(arguments.N, rank, size);
	int start_row = start_end->start;
	int end_row = start_end->end;
	free(start_end);
	
	// allocateMatrices(&arguments, rank, size);
	allocateMatrices2(&arguments, start_row, end_row, rank, size);

	initMatrices(&arguments, &options, rank, size);


	// displayMatrix(&arguments, &results, &options, start_row, end_row, rank, size);
	// sleep(1);
	MPI_Barrier(MPI_COMM_WORLD);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options, start_row, end_row, rank, size);
	gettimeofday(&comp_time, NULL);


	// displayStatistics(&arguments, &results, &options);
	// DisplayMatrix(&arguments, &results, &options, rank, size, 1, arguments.num_rows - 2);
	displayMatrix(&arguments, &results, &options, start_row, end_row, rank, size);

	sleep(1);
	MPI_Barrier(MPI_COMM_WORLD);

	// printf("[%d] Displayed Statistics successfully...\n", rank);
	freeMatrices(&arguments);
	// printf("[%d] Freed Matrices successfully...\n", rank);

	// Finalize MPI Stuff
	MPI_Finalize();

	return 0;
}
