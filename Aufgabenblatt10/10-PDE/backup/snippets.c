// Synchronize the starting row for each process
	int start_row = malloc(sizeof(MPI_INT) * 1);
	int end_row = malloc(sizeof(MPI_INT) * 1);
	if (rank == 0)
	{
		start_row = 1 + (rank * (arguments.num_rows));
		end_row = ((rank + 1) * (arguments.num_rows));
		MPI_Send(end_row, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	else if (rank == size - 1)
	{
		MPI_Recv(&start_row, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		start_row++;
		end_row = arguments.num_rows;
	}
	else
	{
		MPI_Recv(&start_row, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		start_row++;
		end_row = start_row + arguments.num_rows;
		MPI_Send(&end_row, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
	}


/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
	printf("[%d] from: %d, to: %d\n", rank, from, to);

  	int const elements = 8 * options->interlines + 9;

	int x, y;
	double** Matrix = arguments->Matrix[results->m];
	MPI_Status status;

  /* first line belongs to rank 0 */
	if (rank == 0) from--;

  /* last line belongs to rank size - 1 */
	if (rank + 1 == size) to++;

	if (rank == 0) printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		int line = y * (options->interlines + 1);

		if (rank == 0)
		{
			/* check whether this line belongs to rank 0 */
			if (line < from || line > to)
			{
				/* use the tag to receive the lines in the correct order the line is stored in Matrix[0], because we do not need it anymore */
				MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
			}
		}
		else
		{
			if (line >= from && line <= to)
			{
				/* if the line belongs to this process, send it to rank 0 (line - from + 1) is used to calculate the correct local address */
				MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
			}
		}
		// printf("[%d] Send Recv done...!\n", rank);

		if (rank == 0)
		{
			for (x = 0; x < 9; x++)
			{
				int col = x * (options->interlines + 1);

				if (line >= from && line <= to)
				{
					/* this line belongs to rank 0 */
					printf("%7.4f", Matrix[line][col]);
				}
				else
				{
					/* this line belongs to another rank and was received above */
					printf("%7.8f", Matrix[0][col]);
				}
			}

			printf("\n");
		}
		printf("[%d] DisplayMatrix barrier reached...!\n", rank);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	printf("[%d] DisplayMatrix done...!\n", rank);
	fflush(stdout);
}
