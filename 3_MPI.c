#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <mpi.h>

//Funcion que compara los numeros de la ventana
int compare(const void *a, const void *b)
{
    return (*(uint8_t *)a - *(uint8_t *)b);
}

//Funcion donde usamos el qsort para sacar la mediana de la ventana de 9 elementosa
uint8_t calcMediana(uint8_t ventana[9])
{
    qsort(ventana, 9, sizeof(uint8_t), compare);
    return ventana[4];
}

//Funcion que Realiza la matriz Extendida a partir de la matriz normal
void calcMatExt(uint8_t** mat, uint8_t** matExt, int fila, int col)
{
    matExt[0][0] = mat[1][1];
    matExt[fila + 1][0] = mat[fila - 1][1];
    matExt[0][col + 1] = mat[1][col - 1];
    matExt[fila + 1][col + 1] = mat[fila - 1][col - 1];

    for (int j = 0; j < col; j++) 
    {
        matExt[0][j + 1] = mat[0][j];
    }

    for (int j = 0; j < col; j++) 
    {
        matExt[fila + 1][j + 1] = mat[fila - 1][j];
    }

    for (int i = 0; i < fila; i++) 
    {
        matExt[i + 1][0] = mat[i][0];
    }

    for (int i = 0; i < fila; i++) 
    {
        matExt[i + 1][col + 1] = mat[i][col - 1];
    }

    for(int i = 0; i < fila; i++)
    {
        for(int j = 0; j < col; j++)
        {
            matExt[i + 1][j + 1] = mat[i][j];
        }
    }

}

int main(int argc, char *argv[]) 
{
	double tiempo_inicio, tiempo_final;
    int myrank, nproces, fil, col, filaRep, filaInic;
	uint8_t **mat, **matExt, **matRep, **matAux;
	
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproces);
	MPI_Status status;

	if(myrank == 0)
	{
		char *nom_fich = argv[1];
        FILE *fich = fopen(nom_fich, "rb");
		fil = atoi(argv[2]);
		col = atoi(argv[3]);
		tiempo_inicio = MPI_Wtime();
		//Reservo memoria
		mat = (uint8_t**)calloc(fil, sizeof(uint8_t*));
		for(int i = 0; i < fil; i++)
		{
			mat[i] = (uint8_t*)calloc(col, sizeof(uint8_t));
		}
		
		//Creo la matriz con los datos del archivo
		for(int i = 0; i < fil; i++)
		{
			fread(mat[i], sizeof(uint8_t),  col, fich);
		}
		
		//Reservo la memoria de la matriz extendida
		matExt = (uint8_t**)calloc((fil + 2), sizeof(uint8_t*));
		for(int i = 0; i < fil + 2; i++)
		{
			matExt[i] = (uint8_t*)calloc((col + 2), sizeof(uint8_t));
		}
		
		//Calculo la matriz Extendida con la funcion creada
		calcMatExt(mat, matExt, fil, col);
		
		//Calculo las filas a repartir
		filaRep = ((fil) / nproces);
		int restoRep = (fil) % nproces;
		for(int i = 1; i < nproces; i++)
		{
			int filAux;
			if (i == nproces - 1) {
				filAux = (filaRep + restoRep) + 2;
			} else {
				filAux = filaRep + 2;
			}
			//La fila inicial le sumo 1 para saltarme el borde de la matriz Extendida
			filaInic = (filaRep * i) + 1;
			int colAux = col + 2;
			
			//Envio las filas y columnas que voy a enviar a cada proceso
			MPI_Send(&filAux, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&colAux, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
			
			//Envio la MatrizExtendida fila por fila a cada proceso
			for(int j = 0; j < filAux; j++)
			{	
				MPI_Send(matExt[(filaInic - 1) + j], colAux, MPI_UINT8_T, i, 0, MPI_COMM_WORLD);
				// printf("Enviando %d fila al proceso %d\n", j, i);
			}
		}
		fclose(fich);
		
	}else{
		
		//Recivo las filas y columnas
		MPI_Recv(&fil, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&col, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
		// printf("Proceso: %d, fila %d, col %d\n",myrank,  fil, col);
		printf("filas %d\n", fil);
		//Reservo memoria para la matriz que voy a recibir
		matExt = (uint8_t**)calloc((fil), sizeof(uint8_t*));
		for(int i = 0; i < fil; i++)
		{
			matExt[i] = (uint8_t*)calloc((col), sizeof(uint8_t));
		}
		
		//Recibo la matriz
		for(int i = 0; i < fil; i++)
		{
			MPI_Recv((matExt[i]), (col), MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, &status);
			// printf("recibiendo: fila %d, Procesp %d\n", i, myrank);
		}
	}
	
	
	if(myrank != 0) {
		
		//Reservo memoria para la matriz procesada
		matRep = (uint8_t**)malloc((fil - 2) * sizeof(uint8_t*));
		for (int i = 0; i < fil - 2; i++) {
			matRep[i] = (uint8_t*)malloc((col - 2) * sizeof(uint8_t));
		}

		//Proceso la matriz creando ventanas y ejecutando la funcion para calcular la mediana
		for (int i = 1; i < fil - 1; i++) 
		{
			for (int j = 1; j < col - 1; j++) 
			{
				uint8_t ventana[9] = 
				{
					matExt[i-1][j-1], matExt[i-1][j], matExt[i-1][j+1],
					matExt[i][j-1],   matExt[i][j],   matExt[i][j+1],
					matExt[i+1][j-1], matExt[i+1][j], matExt[i+1][j+1]
				};
				matRep[i-1][j-1] = calcMediana(ventana);
				
			}
		}

		//Envio la matriz procesada al proceso 0
		for(int i = 0; i < fil - 2; i++) 
		{
			MPI_Send(matRep[i], col - 2, MPI_UINT8_T, 0, 3, MPI_COMM_WORLD);
			//printf("Proceso %d: Enviando fila %d\n", myrank, i);
		}
		
		//LIBERO MEMORIA
		for(int i = 0; i < fil - 2; i++)
		{
			free(matRep[i]);
		}
		for(int i = 0; i < fil; i++)
		{
			free(matExt[i]);
		}
	} 
	else {
		//Reservo memoria para la matriz que almacenara todo el archivo procesado
		uint8_t **matAux = (uint8_t**)malloc((fil) * sizeof(uint8_t*));
		for (int i = 0; i < fil; i++) 
		{
			matAux[i] = (uint8_t*)malloc((col) * sizeof(uint8_t));
		}
		
		//Proceso la parte a procesar del proceso 0 y lo guardo en la matriz procesada
		for(int i = 1; i < filaRep + 1; i++)
		{
			for(int j = 1; j < col + 1; j++)
			{
				uint8_t ventana[9] = {
					matExt[i-1][j-1], matExt[i-1][j], matExt[i-1][j+1],
					matExt[i][j-1],   matExt[i][j],   matExt[i][j+1],
					matExt[i+1][j-1], matExt[i+1][j], matExt[i+1][j+1]
				};
				matAux[i-1][j-1] = calcMediana(ventana);
			}
		}
		
		//Recivo las matrices procesadas de cada processo y las junto en la matAux
		for (int i = 1; i <	 nproces; i++) 
		{
			
			filaRep = ((fil) / nproces);
			int restoRep = (fil) % nproces;
			int filAux;
			if (i == nproces - 1) 
			{
				filAux = (filaRep + restoRep);
			} else {
				filAux = filaRep;
			}
			filaInic = (filaRep * i);

			for (int j = 0; j < filAux; j++) 
			{
				MPI_Recv(matAux[filaInic + j], col, MPI_UINT8_T, i, 3, MPI_COMM_WORLD, &status);
				//printf("Proceso 0: Recibida fila %d de proceso %d: %d\n", j, i, filaRep);
			}
		}
		
		//Abro un archivo de escritura
		FILE *fich2 = fopen("resultado.raw", "w");
		
		//Copio la matriz procesada en el archivo creado.
		for(int i = 0; i < fil; i++)
		{
			fwrite(matAux[i], col, sizeof(uint8_t), fich2); 
		}
		fclose(fich2);
		
		tiempo_final = MPI_Wtime();
		printf("%f", tiempo_final-tiempo_inicio);
		
		//LIBERO MEMORIA
		for(int i = 0; i < fil; i++)
		{
			free(matAux[i]);
			free(mat[i]);
		}
		for(int i = 0; i < fil + 2; i++)
		{
			free(matExt[i]);
		}
	}

    MPI_Finalize();
    return EXIT_SUCCESS;
}
