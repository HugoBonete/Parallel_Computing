z#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>

void crearMatriz(uint8_t** mat, FILE *fich, int fil, int col) {
    for (int i = 0; i < fil; i++) {
        fread(mat[i], sizeof(uint8_t), col, fich);
    }
}

int main(int argc, char *argv[]) 
{
    int myrank, nproces, fil, col;
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproces);
	MPI_Status status;
	

    if (myrank == 0) {
		
		double tiempo_inicio, tiempo_final;
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
        int histograma[256] = {0};
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
		
		FILE *fich2 = fopen("histograma.txt", "w");
		
		for(int i = 0; i < 256; i++)
		{
			fprintf(fich2, "El color %d se ha repetido %d veces\n", i, histograma[i]);
		}
		
		tiempo_final = MPI_Wtime();
		printf("%f", tiempo_final-tiempo_inicio);

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
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
