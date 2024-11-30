#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>
#include <limits.h>

void crearMatriz(uint8_t** mat, FILE *fich, int fil, int col) {
    for (int i = 0; i < fil; i++) {
        fread(mat[i], sizeof(uint8_t), col, fich);
    }
}



int main(int argc, char *argv[]) 
{
    int myrank, nproces, fil, col, *histograma, auxRepHisto, max = 0, min;
	int *histoRepar;
	double tiempo_inicio, tiempo_final;
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproces);
	MPI_Status status;
	

    if (myrank == 0) {
		char *nom_fich = argv[1];
        FILE *fich = fopen(nom_fich, "rb");
		fil = atoi(argv[2]);
		col = atoi(argv[3]);
	
		tiempo_inicio = MPI_Wtime();
		int filaRepartida = fil / nproces;
		int restoFila = fil % nproces;

        // Reservo la memoria para la matriz
        uint8_t **mat = (uint8_t**) malloc(fil * sizeof(uint8_t*));
        for (int i = 0; i < fil; i++) 
		{
            mat[i] = (uint8_t*) malloc(col * sizeof(uint8_t));
        }

        // Creo la matriz leyendo el archivo
        crearMatriz(mat, fich, fil, col);
		
        // Enviar partes de la matriz a otros procesos
        for (int i = 1; i < nproces; i++) {
            int auxRepartida = (i == nproces - 1) ? (filaRepartida + restoFila) : filaRepartida;
			for (int j = 0; j < auxRepartida; j++) 
			{
				MPI_Send(mat[filaRepartida * i + j], col, MPI_BYTE, i, 99, MPI_COMM_WORLD);
			}
			MPI_Send(&col, 1, MPI_INT, i, 98, MPI_COMM_WORLD);
			MPI_Send(&auxRepartida, 1, MPI_INT, i, 97, MPI_COMM_WORLD);			
        }
		
        // Inicializo el histograma en el proceso 0
        histograma = (int*)calloc(256, sizeof(int));
		for(int i = 0; i < filaRepartida; i++)
		{
			for(int j = 0; j < col; j++)
			{
				histograma[mat[i][j]]++;
			}
		}
		
        // Recibo los histogramas de los otros procesos
        for (int i = 1; i < nproces; i++) {
            int vecAux[256];
            MPI_Recv(vecAux, 256, MPI_INT, i, 0, MPI_COMM_WORLD, &status);

            for (int j = 0; j < 256; j++) {
                histograma[j] += vecAux[j];
            }
        }
        // Liberar memoria de la matriz en el proceso 0
        for (int i = 0; i < fil; i++) {
            free(mat[i]);
        }
        free(mat);
    }
	else 
	{
		int auxFilaRec;
		MPI_Recv(&auxFilaRec, 1, MPI_INT, 0, 97, MPI_COMM_WORLD, &status);
		MPI_Recv(&col, 1, MPI_INT, 0, 98, MPI_COMM_WORLD, &status);
        uint8_t **matAux = (uint8_t**) malloc(auxFilaRec * sizeof(uint8_t*));
        for (int i = 0; i < auxFilaRec; i++) 
		{
            matAux[i] = (uint8_t*) malloc(col * sizeof(uint8_t));
        }
		
		int *vecHisto = (int*)calloc(256, sizeof(int));
		
        // Recibo el bloque de la matriz
		for(int i = 0; i < auxFilaRec; i++)
		{
			MPI_Recv((matAux[i]), col, MPI_BYTE, 0, 99, MPI_COMM_WORLD, &status);
			
			//Cuento los colores del histograma
			for (int j = 0; j < col; j++) {
                vecHisto[matAux[i][j]]++;
            }
			
		}
        // Envío el histograma al proceso 0
        MPI_Send(vecHisto, 256, MPI_INT, 0, 0, MPI_COMM_WORLD);
		
		free(vecHisto);
		for(int i = 0; i < auxFilaRec; i++)
		{
			free(matAux[i]);
		}
    }
	//ENTREGA B
	if(myrank == 0)
	{
		int repartHisto = 256 / nproces;
		int restoHisto = 256 % nproces;
		for(int i = 1; i < nproces; i++)
		{
			auxRepHisto = (i == nproces - 1) ? (repartHisto + restoHisto) : repartHisto;
			int inicioAux = auxRepHisto * i;
			MPI_Send(&auxRepHisto, 1, MPI_INT, i, 51, MPI_COMM_WORLD);
			MPI_Send(&histograma[inicioAux], auxRepHisto, MPI_INT, i, 50, MPI_COMM_WORLD);
		}
		histoRepar = (int*)calloc(auxRepHisto, sizeof(int));
		for(int i = 0; i < auxRepHisto; i++)
		{
			histoRepar[i] = histograma[i];
		}
	}else{
		MPI_Recv(&auxRepHisto, 1, MPI_INT, 0, 51, MPI_COMM_WORLD, &status);
		histoRepar = (int*)calloc(auxRepHisto, sizeof(int));
		MPI_Recv(histoRepar, auxRepHisto, MPI_INT, 0, 50, MPI_COMM_WORLD, &status);
	}
	
	if(myrank != 0)
	{
		min = INT_MAX;
		for(int i = 0; i < auxRepHisto; i++)
		{
			if(i == 0)
			{
				max = histoRepar[i];
				if(histoRepar[i] != 0)
				{
					min = histoRepar[i];
				}else{
					min = min;
				}
			}
			else{
				if(max < histoRepar[i])
				{
					max = histoRepar[i];
				}
				if(min > histoRepar[i] && histoRepar[i] != 0)
				{
					min = histoRepar[i];
				}
			}
			
		}
		
		MPI_Send(&max, 1, MPI_INT, 0, 52, MPI_COMM_WORLD);
		MPI_Send(&min, 1, MPI_INT, 0, 53, MPI_COMM_WORLD);
	}else{
		FILE *fich2 = fopen("histograma.txt", "w");
		for(int i = 0; i < auxRepHisto; i++)
		{
			min = INT_MAX;
			if(i == 0)
			{
				max = histoRepar[i];
				if(histoRepar[i] != 0)
				{
					min = histoRepar[i];
				}else{
					min = min;
				}
			}
			else{
				if(max < histoRepar[i])
				{
					max = histoRepar[i];
				}
				if(min > histoRepar[i])
				{
					min = histoRepar[i];
				}
			}
		}
		for(int i = 1; i < nproces; i++)
		{
			int auxMax, auxMin;
			MPI_Recv(&auxMax, 1, MPI_INT, i, 52, MPI_COMM_WORLD, &status);
			MPI_Recv(&auxMin, 1, MPI_INT, i, 53, MPI_COMM_WORLD, &status);
			max = (auxMax > max) ? auxMax : max;
			min = (auxMin < min) ? auxMin : min;
		}
		fprintf(fich2, "El maximo es : %d.\n El minimo es: %d\n", max, min);
		tiempo_final = MPI_Wtime();
		printf("%f", tiempo_final-tiempo_inicio);
		for(int i = 0; i < 256; i++)
		{
			fprintf(fich2, "El color %d se ha repetido %d veces\n", i, histograma[i]);
		}
		
		fclose(fich2);
		free(histoRepar);
	}
    MPI_Finalize();
    return EXIT_SUCCESS;
}
