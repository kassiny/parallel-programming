//программа считает скалярное произведение векторов 
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "mpi.h"

int main(int argc, char** argv)
{
	int num, rank, i, blocksize;
	int size, oldSize;
	double *v1, *v2, *v11, *v22;
	double prod, scalarProd;
	MPI_Status status;

	srand(time(NULL));
	v1 = v2 = NULL;
	size = atoi(argv[1]);
	prod = 0;

	/*if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
	{
		printf("pmi can't init");
	}*/
	MPI_Init(0, 0);

	MPI_Comm_size(MPI_COMM_WORLD, &num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (size % num == 0)
	{
		
		v1 = (double*)malloc(sizeof(double) * size);
		v2 = (double*)malloc(sizeof(double) * size);
		oldSize = size;
	}
	else
	{
		v1 = (double*)malloc(sizeof(double) * (size + (num - size % num)));
		v2 = (double*)malloc(sizeof(double) * (size + (num - size % num)));
		oldSize = size;
		size = (size + (num - size % num));
	}

	if (rank == 0) 
	{
		
		if (oldSize % num == 0)
		{
			// генерируем v1
			printf("v1: ");
			for (i = 0; i < size; i++) {
				v1[i] = (rand() % 100)/10.0;
				printf("%f ", v1[i]);
			}
			printf("\n");
			// генерируем v2
			printf("v2: ");
			for (int i = 0; i < size; i++) {
				v2[i] = (rand() % 100) / 10.0;
				printf("%f ", v2[i]);
			}
			printf("\n");
			size = size + size%num;
		}
		else
		{
			// генерируем v1
			printf("v1: ");
			for (i = 0; i < oldSize; i++) {
				v1[i] = (rand() % 100) / 10.0;
				printf("%f ", v1[i]);
			}
			for (; i < oldSize + num - oldSize%num; i++)
			{
				v1[i] = 0;
				printf("%f ", v1[i]);
			}
			printf("\n");
			// генерируем v2
			printf("v2: ");
			for (int i = 0; i < oldSize; i++) {
				v2[i] = (rand() % 100) / 10.0;
				printf("%f ", v2[i]);
			}
			for (; i < oldSize + num - oldSize%num; i++)
			{
				v2[i] = 0;
				printf("%f ", v2[i]);
			}
			printf("\n");
		}
	}
	MPI_Bcast(v1, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(v2, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	blocksize = size / num;
	//MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);4
	for (i = blocksize*rank; i < (rank + 1)*blocksize; i++)
	{
		prod += v1[i] * v2[i];
	}
	if (rank == 0)
	{
		for (i = 1; i < num; i++)
		{
			MPI_Recv(&scalarProd, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			prod += scalarProd;
		}
	}
	else
	{
		MPI_Send(&prod, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	if (rank == 0)
	{
		printf("scalar product is %f", prod);
	}

	MPI_Finalize();
	free(v1);
	free(v2);
	return 0;
}