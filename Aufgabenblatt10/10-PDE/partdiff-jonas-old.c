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
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
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
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* sequential variant of initVariables				                        */
/* ************************************************************************ */
static
void
initVariablesSequential (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

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
allocateMatrices (struct calculation_arguments* arguments, int start, int end, int rank, int size)
{
	uint64_t i, j;

	uint64_t range = end - start + 1;
	uint64_t bufferlines = 0;

	if(rank == 0 || rank == size - 1) 
	{
		bufferlines = 1;
	}
	else if(rank != 0 && rank != size - 1)
	{
		bufferlines = 2;
	}
	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (range + bufferlines) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));
	// printf("Rank %d allocates %lu lines, with %lu bufferlines\n", rank, range + bufferlines, bufferlines);

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((range + bufferlines) * sizeof(double*));

		for (j = 0; j < range + bufferlines; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (range + bufferlines) * (N + 1)) + (j * (N + 1));
		}
	}
	
}

/**
 * sequential variant of allocateMatrices
*/
static
void
allocateMatricesSequential (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}


/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options, int start, int end, int rank, int size)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	uint64_t range = end - start + 1;
	uint64_t bufferlines = 0;

	if(rank == 0 || rank == size - 1) 
	{
		bufferlines = 1;
	}
	else if(rank != 0 && rank != size - 1)
	{
		bufferlines = 2;
	}

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < range + bufferlines; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0)
    {
		// Be aware that the 'i' variable of the initial program is named 'r' here
        for(uint64_t r = 0; r < range; r++)
        {
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < arguments->num_matrices; j++)
                {
					int true_row = r + start;	// (== end - r) The actual row in the combined matrix 

					// h * (N - true_row) invertiert die Ausgabe Elemente jedes Prozesses
					// Index [range - 1 - (range - r - 1)] invertiert die lesweise (nachher: von oben nach unten)
					/**
					 * [4] = 2 			-> 			[0] = 0
					 * [3] = 1 			-> 			[1] = 1
					 * [2] = 0 			-> 			[2] = 2
					 * [1] = 4 			-> 			[3] = 3
					 * [0] = 3 			-> 			[4] = 4
					*/

                    if(rank == 0) {
                        Matrix[j][r][0] = 1 + (1 - (h * (true_row))); // Linke Kante
                        Matrix[j][range - 1 - (range - r - 1)][N] = (h * (N - true_row)); // Rechte Kante
						// if(i == 0 && j == 0)
						// {
						// 	printf("rank %i: (%f * (%li - %i)) = %f\n", rank, h, N, true_row, (h * (N - true_row)));
						// }
						
                    }
                    if(rank == size - 1) {
                        Matrix[j][r+1][0] = 1 + (1 - (h * (true_row))); // Linke Kante
                        Matrix[j][range - (range - r) + 1][N] = (h * (N - true_row)); // Rechte Kante
						// if(i == 0 && j == 0)
						// {
						// 	printf("rank %i: (%f * (%li - %i)) = %f\n", rank, h, N, true_row, (h * (N - true_row)));
						// }
                    }
                    else if(rank != 0 && rank != size - 1) {
                        Matrix[j][r+1][0] = 1 + (1 - (h * (true_row))); // Linke Kante
                        Matrix[j][range - (range - r) + 1][N] = (h * (N - true_row)); // Rechte Kante
						// if(i == 0 && j == 0)
						// {
						// 	printf("rank %i: (%f * (%li - %i)) = %f\n", rank, h, N, true_row, (h * (N - true_row)));
						// }
                    }

                    if(rank == 0) {
                        Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
                    }
                    if(rank == size-1) {
                        Matrix[j][range + bufferlines - 1][i] = 1 - (h * i); // Untere Kante
                    }
                }
            }
        }
    }
	sleep(1);
}


/* ************************************************************************ */
/* initMatricesSequential: Sequential variant of initMatrices       		*/
/* ************************************************************************ */
static
void
initMatricesSequential (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for(i = 0; i < N; i++)
		{
			for (j = 0; j < arguments->num_matrices; j++)
			{
				Matrix[j][i][0] = 1 + (1 - (h * i)); // Linke Kante
				Matrix[j][N][i] = 1 - (h * i); // Untere Kante
				Matrix[j][N - i][N] = h * i; // Rechte Kante
				Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
			}
		}
	}
}

start_end_results_t*
calcStartEnd (int N, int rank, int size)
{
    int lines = (N + 1) / size;
    int extra = (N + 1) % size;

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

    // printf("Lines = %d, Startline = %d, Endline = %d\n", endline-startline+1, startline, endline);

	start_end_results_t* results = malloc(sizeof(start_end_results_t));
	results->start = startline;
	results->end = endline;

    return results;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, int start, int end, int rank, int size)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	// Einfach als arg->range - 1 speichern
	int range = end - start;
	int bufferlines = 0;

	if(rank == 0 || rank == size - 1) 
	{
		bufferlines = 1;
	}
	else if(rank != 0 && rank != size - 1)
	{
		bufferlines = 2;
	}

	// printf("Rank %d calculates %d lines, from local index %d to local index %d\n", rank, range + bufferlines - 1, 1, range + bufferlines - 1);

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

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

		maxResiduum = 0;

		if(rank != 0)
		{
			 // Get last line from previous process
			MPI_Recv(Matrix_In[0], N+1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if(rank != size-1) 
		{
			// Send last line to next process
			MPI_Ssend(Matrix_In[range + bufferlines - 1], N+1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			// TODO die beiden ifs zusammenfassen? Ja
		}
		if(rank != size-1)
		{
			 // Get first line from next process
			MPI_Recv(Matrix_In[range + bufferlines], N+1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if(rank != 0) 
		{
			// Send first line to previous process
			MPI_Ssend(Matrix_In[1], N+1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);

		/* over all rows */
		for (i = 1; i < range + bufferlines; i++)
		{
			double fpisin_i = 0.0;

			int subtract = 0;
			if(rank != 0) 
			{
				subtract = -1;
			}
			
			int true_line = start + i + subtract;

			if(term_iteration == 1) {
			}
			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)true_line);
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

		// Sync residuum
		double global_maxResiduum;

        MPI_Reduce(&maxResiduum, &global_maxResiduum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        MPI_Bcast(&global_maxResiduum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        results->stat_iteration++;
        results->stat_precision = global_maxResiduum;

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

/**
 * sequential variant of calculate
*/
static
void
calculateSequential (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

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

		maxResiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
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
	}

	results->m = m2;
}

/**
 * Parallel variant of calculate for GS
*/
static
void
calculateGSparallel (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, int start, int end, int rank, int size)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	// Einfach als arg->range - 1 speichern
	int range = end - start;
	int bufferlines = 0;

	if(rank == 0 || rank == size - 1) 
	{
		bufferlines = 1;
	}
	else if(rank != 0 && rank != size - 1)
	{
		bufferlines = 2;
	}

	// printf("Rank %d calculates %d lines, from local index %d to local index %d\n", rank, range + bufferlines - 1, 1, range + bufferlines - 1);

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

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

		maxResiduum = 0;

		if(rank != 0)
		{
			 // Get last line from previous process
			MPI_Recv(Matrix_In[0], N+1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if(rank != size-1) 
		{
			// Send last line to next process
			MPI_Ssend(Matrix_In[range + bufferlines - 1], N+1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			// TODO die beiden ifs zusammenfassen? Ja
		}
		if(rank != size-1)
		{
			 // Get first line from next process
			MPI_Recv(Matrix_In[range + bufferlines], N+1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		if(rank != 0) 
		{
			// Send first line to previous process
			MPI_Ssend(Matrix_In[1], N+1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);

		/* over all rows */
		for (i = 1; i < range + bufferlines; i++)
		{
			double fpisin_i = 0.0;

			int subtract = 0;
			if(rank != 0) 
			{
				subtract = -1;
			}
			
			int true_line = start + i + subtract;

			if(term_iteration == 1) {
			}
			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)true_line);
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

		// Sync residuum
		double global_maxResiduum;

        MPI_Reduce(&maxResiduum, &global_maxResiduum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        MPI_Bcast(&global_maxResiduum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        results->stat_iteration++;
        results->stat_precision = global_maxResiduum;

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

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int start, int end, int rank, int size)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	//printf("Matrix: ");
	int range = end - start;
	int add = 0;

	if(rank != 0)
	{
		add = 1;
	}

	int ok = 1;

	// Print in order
    if(rank != 0) 
	{
        MPI_Recv(&ok, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

	// Do the Printing
    for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			// global index == true_line ?
			int global_index = y * (interlines + 1);
			if(global_index - start >= 0 && global_index - start <= range)
			{
				printf("%7.8f ", Matrix[global_index - start + add][x * (interlines + 1)]);
			}
		}

		int i_index = y * (interlines + 1);
		if(i_index - start >= 0 && i_index - start <= range) 
		{
			printf ("\n");
		}
	}
	// Allow next Process to start
    if(rank != size-1) 
	{
        MPI_Ssend(&ok, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);
	if(rank == 0) {
		printf ("\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// print detailed matrix
	// int ok2 = 1;
	// if(rank != 0) 
	// {
    //     MPI_Recv(&ok2, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // }
	// if(rank!=0 && rank!=size-1) {
	// 	range+=1;
	// }
    // for (y = 0; y <= range+1; y++)
	// {
	// 	for (x = 0; x <= (int)arguments->N; x++)
	// 	{
	// 		printf("%7.8f ", Matrix[y][x]);
	// 	}
	// 	printf ("\n");
	// }
	// printf ("\n");
    // if(rank != size-1) 
	// {
	// 	sleep(1);
    //     MPI_Ssend(&ok2, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    // }

	fflush (stdout);
}

/**
 * Sequential variant of displayMatrix
*/
static
void
displayMatrixSequential (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.8f ", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}



/****************************************************************************/
/* Main for one process -> Sequential.										*/
/****************************************************************************/
int
mainOneProcess(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	initVariablesSequential(arguments, results, options);

	allocateMatricesSequential(arguments);
	initMatricesSequential(arguments, options);

	gettimeofday(&start_time, NULL);
	calculateSequential(arguments, results, options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(arguments, results, options);
	displayMatrixSequential(arguments, results, options);

	freeMatrices(arguments);

	return 0;
}

/****************************************************************************/
/* Main for GS Parallel														*/
/****************************************************************************/
int
mainGSParallel(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size)
{
	initVariables(arguments, results, options);
	int N = arguments->N;

	start_end_results_t* start_end_results = calcStartEnd(N, rank, size);
	int start = start_end_results->start;
	int end = start_end_results->end;
	

	allocateMatrices(arguments, start, end, rank, size);
	initMatrices(arguments, options, start, end, rank, size);

	gettimeofday(&start_time, NULL);
	calculateGSparallel(arguments, results, options, start, end, rank, size);
	gettimeofday(&comp_time, NULL);

    if(rank == 0) 
	{
        displayStatistics(arguments, results, options);
    }
	
	MPI_Barrier(MPI_COMM_WORLD);
	displayMatrix(arguments, results, options, start, end, rank, size);

	freeMatrices(arguments);

	return 0;
}

int
mainJacobiParallel(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size)
{
	initVariables(arguments, results, options);
	int N = arguments->N;

	start_end_results_t* start_end_results = calcStartEnd(N, rank, size);
	int start = start_end_results->start;
	int end = start_end_results->end;
	

	// TODO only  calculate range etc once
	// TODO Do OpenMP stuff
	// TODO calcStartEnd aufhübschen
	// TODO prints entfernen

	allocateMatrices(arguments, start, end, rank, size);
	initMatrices(arguments, options, start, end, rank, size);

	gettimeofday(&start_time, NULL);
	calculate(arguments, results, options, start, end, rank, size);
	gettimeofday(&comp_time, NULL);

    if(rank == 0) 
	{
        displayStatistics(arguments, results, options);
    }
	
	MPI_Barrier(MPI_COMM_WORLD);
	displayMatrix(arguments, results, options, start, end, rank, size);

	freeMatrices(arguments);

	return 0;
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	askParams(&options, argc, argv);

	// Check if we want to use sequential version (only one process or Gauss-Seidel)
	int use_sequential = (size == 1) ? 1 : 0;

	// If sequential should be used, call mainOneProcess on one process
	if(use_sequential)
	{
		// Only one process has to do anything, everyone else can exit
		if (rank == 0)
		{
			printf("Using sequential version\n");
			mainOneProcess(&arguments, &results, &options);
		}
	}
	else if(options.method == METH_GAUSS_SEIDEL)
	{
		// If we use Gauss-Seidel, call mainGSParallel on all processes
		printf("Using parallel version for GS\n");
		mainGSParallel(&arguments, &results, &options, rank, size);
	}
	else
	{
		// If we use parallel version, call mainJacobiParallel on all processes
		printf("Using parallel version for Jacobi\n");
		mainJacobiParallel(&arguments, &results, &options, rank, size);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}