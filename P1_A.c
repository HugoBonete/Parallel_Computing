#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	int myrank, nproces, col, filaRep, fil, resto;
	double tiempo_inicio, tiempo_final;
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproces);
	MPI_Status status;
	uint8_t **mat;
	int histograma[256] = {0};
	
	//Reparto trabajo
	if(myrank == 0)
	{
		fil = atoi(argv[2]);
		col = atoi(argv[3]);
		resto = fil % nproces;
		filaRep = fil / nproces;
		
		for(int i = 1; i < nproces; i++)
		{
			int auxFila = (i == nproces - 1) ? filaRep + resto : filaRep;
			
			MPI_Send(&auxFila, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&col, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
		}
	}else{
		MPI_Recv(&filaRep, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&col, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
	}
	
	//Reservo memoria y obtengo datos
	if(myrank == 0)
	{
		char *nom = argv[1];
		FILE *fich = fopen(nom, "rb");
		mat = (uint8_t**)malloc(fil * sizeof(uint8_t*));
		for(int i = 0; i < fil; i++)
		{
			mat[i] = (uint8_t*)malloc(col * sizeof(uint8_t));
			fread(mat[i], sizeof(uint8_t), col, fich);
		}	
	}else{
		mat = (uint8_t**)malloc(filaRep * sizeof(uint8_t*));
		for(int i = 0; i < filaRep; i++)
		{
			mat[i] = (uint8_t*)malloc(col * sizeof(uint8_t));
		}	
	}
	
	//Reparto datos
	if(myrank == 0)
	{
		tiempo_inicio = MPI_Wtime();
		for(int i = 1; i < nproces; i++)
		{
			int auxFila = (i == nproces - 1) ? filaRep + resto : filaRep;
			for(int j = 0; j < auxFila; j++)
			{
				MPI_Send(mat[(i * filaRep) + j], col, MPI_UINT8_T, i, 2, MPI_COMM_WORLD);
			}
		}
	}else{
		for(int i = 0; i < filaRep; i++)
		{
			MPI_Recv(mat[i], col, MPI_UINT8_T, 0, 2, MPI_COMM_WORLD, &status);
		}
	}
	
	//Procesamiento
	if(myrank == 0)
	{
		for(int i = 0; i < filaRep; i++)
		{
			for(int j = 0; j < col; j++)
			{
				histograma[mat[i][j]]++;
			}
		}
	}else{
		for(int i = 0; i < filaRep; i++)
		{
			for(int j = 0; j < col; j++)
			{
				histograma[mat[i][j]]++;
			}
		}
	}
	
	//Recibo resultados
	if(myrank != 0)
	{
		MPI_Send(&histograma, 256, MPI_INT, 0, 3, MPI_COMM_WORLD);
	}else{
		for(int i = 1; i < nproces; i++)
		{
			int histAux[256] = {0};
			MPI_Recv(&histAux, 256, MPI_INT, i, 3, MPI_COMM_WORLD, &status);
			for(int j = 0; j < 256; j++)
			{
				histograma[j] += histAux[j];
			}
		}
		FILE *fich2 = fopen("hola.txt", "w");
		for(int i = 0; i < 256; i++)
		{
			fprintf(fich2, "El color %d se ha repetido %d veces\n", i, histograma[i]);
		}
		tiempo_final = MPI_Wtime();
		printf("%f", tiempo_final-tiempo_inicio);		
	}
	
	if(myrank == 0)
	{
		for (int i = 0; i < fil; i++) {
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
	
	MPI_Finalize();
}
