#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

void dividirVec(double *vec, double norma, int tam)
{
	for(int i = 0; i < tam; i++)
	{
		vec[i] = vec[i] / norma; 
	}
}

void calcularVec(double *vec2, double *vec1, double **mat, int tam, int filaRep)
{
	for(int i = 0; i < filaRep; i++)
	{
		vec2[i] = 0.0;
		for(int j = 0; j < tam; j++)
		{
			vec2[i] += mat[i][j] * vec1[j];
			// printf("%f\n", mat[i][j]);
		}
		
	}
}

double calcularNorma(double *vec, int tam, int inicio) 
{
    double norma = 0.0;
    for (int i = 0; i < tam; i++) 
	{
        norma += vec[inicio + i] * vec[inicio + i];
    }
    return norma;
}

int main(int argc, char *argv[])
{
	double **mat, *vec1, *vec2, *vec3, norma;
    int myrank, nproces, tam, numIteraciones, filaRep;
	int *vecProc = (int*)malloc(nproces * sizeof(int));
	int *vecInic = (int*)malloc(nproces * sizeof(int));	
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproces);
	MPI_Status status;
	
	//Reparto trabajo
	if(myrank == 0)
	{
		tam = atoi(argv[2]);
		numIteraciones = atoi(argv[3]);
		filaRep = tam / nproces;
		vecInic[0] = 0;
		vecProc[0] = filaRep;
		for(int i = 0; i < nproces; i++)
		{
			vecProc[i] = (i == nproces - 1) ? filaRep + (tam % nproces) : filaRep;
			vecInic[i] = vecInic[i - 1] + vecProc[i - 1];
		}
	}
	MPI_Bcast(vecProc, nproces, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vecInic, nproces, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(vecProc, 1, MPI_INT, &filaRep, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&tam, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numIteraciones, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//Reservo memoria
	if(myrank == 0)
	{
		mat = (double**)malloc(tam * sizeof(double*));
		for(int i = 0; i < tam; i++)
		{
			mat[i] = (double*)malloc(tam * sizeof(double));
		}		
	}else{
		mat = (double**)malloc(filaRep * sizeof(double*));
		for(int i = 0; i < filaRep; i++)
		{
			mat[i] = (double*)malloc(tam * sizeof(double));
		}	
	}
	vec1 = (double*)malloc(tam * sizeof(double));
	vec2 = (double*)calloc(filaRep, sizeof(double));
	
	//Recibo datos y Reparto datos
	if(myrank == 0)
	{
		char *nom = argv[1];
		FILE *fich = fopen(nom, "rb");
		if(fich == NULL)
		{
			for(int i = 0; i < tam; i++)
			{
				for(int j = 0; j < tam;j++)
				{
					if(i == j)
					{
						mat[i][j] = 1;
					}else{
						mat[i][j] = ((double)rand() / RAND_MAX) * 0.02 - 0.01;
					}
				}
			}
		}else{
			for(int i = 0; i < tam; i++)
			{
				fread(mat[i], sizeof(double), tam, fich);
				vec1[i] = 1;
			}
		}
		for(int i = 1; i < nproces; i++)
		{
			for(int j = 0; j < vecProc[i]; j++)
			{
				MPI_Send(mat[vecInic[i] + j], tam, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
	}else{
		for(int i = 0; i < filaRep; i++)
		{
			MPI_Recv(mat[i], tam, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			// printf("%d : %f\n", i, mat[i][1]);
		}
	}
	MPI_Bcast(vec1, tam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//Procesamiento
	norma = 1;
	for(int i = 0; i < numIteraciones; i++)
	{		
		calcularVec(vec2, vec1, mat, tam, vecProc[myrank]);
		dividirVec(vec2, norma, vecProc[myrank]);
		double auxNorm = calcularNorma(vec1, filaRep, vecInic[myrank]);	
		MPI_Allreduce(&auxNorm, &norma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		norma = sqrt(norma);
		MPI_Allgatherv(vec2, vecProc[myrank], MPI_DOUBLE, vec1, vecProc, vecInic, MPI_DOUBLE, MPI_COMM_WORLD);
		if(myrank == 0)
		{
			printf("%f\n", norma);
		}		
	}
	
	if(myrank == 0)
	{
		for(int i = =; i < tam; i++)
		{
			free(mat[i]);
		}
		free(mat);
	}else{
		for(int i = 0; i < filaRep; i++)
		{
			free(mat[i]);
		}
		free(mat);
	}
	free(vec1);
	free(vec2);
	
	MPI_Finalize();
}