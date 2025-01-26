#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>

//Funcion donde usamos el buble metodo para sacar la mediana de la ventana de 9 elementosa
uint8_t calcMediana(uint8_t ventana[9]) {
    for (int i = 0; i < 8; i++) {
        for (int j = i + 1; j < 9; j++) {
            if (ventana[i] > ventana[j]) {
                uint8_t temp = ventana[i];
                ventana[i] = ventana[j];
                ventana[j] = temp;
            }
        }
    }
    return ventana[4];
}

//Funcion para calcular las filas y columnas extendidads de la matriz Extendida 
void calcMatExt(uint8_t** mat, int fila, int col)
{
    mat[0][0] = mat[2][2];
    mat[fila - 1][0] = mat[fila - 3][2];
    mat[0][col - 1] = mat[2][col - 3];
    mat[fila - 1][col - 1] = mat[fila - 3][col - 3];

    for (int j = 1; j < col - 1; j++) 
    {
        mat[0][j] = mat[2][j];
		mat[fila - 1][j] = mat[fila - 3][j];
    }
	
    for (int i = 1; i < fila; i++) 
    {
        mat[i][0] = mat[i][2];
		mat[i][col - 1] = mat[i][col - 3];
    }
}

int main(int argc, char *argv[])
{
	int myrank, nproces, col, filaRep, fil, resto, filExt, colExt, *vecInic;
	double tiempo_inicio, tiempo_final;
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproces);
	MPI_Status status;	
	uint8_t **matExt, **matFiltrada;
	int *vecRep = (int*)malloc(nproces * sizeof(int));
	vecInic = (int*)malloc(nproces * sizeof(int));
	
	//Reparto trabajo
	if(myrank == 0)
	{
		char *nom = argv[1];
		fil = atoi(argv[2]);
		col = atoi(argv[3]);
		filExt = fil + 2, colExt = col + 2;
		filaRep = fil / nproces;
		vecRep[0] = filaRep + 2;
		vecInic[0] = 0;		
		for(int i = 1; i < nproces; i++)
		{
			int aux = (i == nproces - 1) ? filaRep + (fil % nproces) + 2 : filaRep + 2;
			vecRep[i] = aux;
			vecInic[i] = (filaRep * i);
			// printf("%d\n", vecRep[i]);
			printf("%d\n", vecInic[i]);
			if(i > 0)
			{
				MPI_Send(&vecRep[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&colExt, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			}
		}
	}else{
		MPI_Recv(&filExt, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&colExt, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
	}
	
	//Reservo Memoria
	if(myrank == 0)
	{
		matExt = (uint8_t**)malloc((filExt) * sizeof(uint8_t*));
		for(int i = 0; i < filExt; i++)
		{
			matExt[i] = (uint8_t*)malloc((colExt) * sizeof(uint8_t));
		}
		matFiltrada = (uint8_t**)malloc(fil * sizeof(uint8_t*));
		for(int i = 0; i < fil; i++)
		{
			matFiltrada[i] = (uint8_t*)malloc(col * sizeof(uint8_t));
		}
	}else{
		matExt = (uint8_t**)malloc((filExt) * sizeof(uint8_t*));
		for(int i = 0; i < filExt; i++)
		{
			matExt[i] = (uint8_t*)malloc((colExt) * sizeof(uint8_t));
		}
		matFiltrada = (uint8_t**)malloc((filExt - 2) * sizeof(uint8_t*));
		for(int i = 0; i < filExt - 2; i++)
		{
			matFiltrada[i] = (uint8_t*)malloc((colExt - 2) * sizeof(uint8_t));
		}
	}
	
	//Obtengo datos y los envio
	if(myrank == 0)
	{
		char *nom = argv[1];
		FILE *fich = fopen(nom, "rb");
		for(int i = 1; i < filExt - 1; i++)
		{
			fread(&matExt[i][1], sizeof(uint8_t), col, fich);
		}
		calcMatExt(matExt, filExt, colExt);
		
		for(int i = 1; i < nproces; i++)
		{
			for(int j = 0; j < vecRep[i]; j++)
			{
				MPI_Send(matExt[vecInic[i] + j], colExt, MPI_UINT8_T, i, 2, MPI_COMM_WORLD);
			}
		}
	}else{
		for(int i = 0; i < filExt; i++)
		{
			MPI_Recv(matExt[i], colExt, MPI_UINT8_T, 0, 2, MPI_COMM_WORLD, &status);
		}
	}
	
	//Procesamiento
	for(int i = 1; i < filExt - 1; i++) {
		for (int j = 1; j < colExt - 1; j++) {
		uint8_t ventana[9] = {
			matExt[i - 1][j - 1], matExt[i - 1][j], matExt[i - 1][j + 1],
			matExt[i][j - 1],     matExt[i][j],     matExt[i][j + 1],
			matExt[i + 1][j - 1], matExt[i + 1][j], matExt[i + 1][j + 1]
		};
		matFiltrada[i - 1][j - 1] = calcMediana(ventana);
		}
	}
	
	//Envio Resultados
	if(myrank != 0)
	{
		for(int i = 0; i < filExt - 2; i++)
		{
			MPI_Send(matFiltrada[i], colExt - 2, MPI_UINT8_T, 0, 3, MPI_COMM_WORLD);
		}
	}else{
		for(int i = 1; i < nproces; i++)
		{
			for(int j = 0; j < vecRep[i] - 2; j++)
			{
				MPI_Recv(matFiltrada[vecInic[i] + j], col, MPI_UINT8_T, i, 3, MPI_COMM_WORLD, &status);
			}
		}
		FILE *fich2 = fopen("resultado.raw", "w");
		
		for(int i = 0; i < fil; i++)
		{
			fwrite(matFiltrada[i], col, sizeof(uint8_t), fich2); 
		}
		fclose(fich2);		
	}
	
	//Libero memorio
	for(int i = 0; i < filExt; i++)
	{
		free(matExt[i]);
	}
	for(int i = 0; i < filExt - 2; i ++)
	{
		free(matFiltrada[i]);
	}
	free(matExt);
	free(matFiltrada);
	
	MPI_Finalize();
}