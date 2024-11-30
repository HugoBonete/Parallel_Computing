#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

//Creo una funcion que cree una matriz
void crearMatriz(double** mat, FILE *fich, int tam) {
    for (int i = 0; i < tam; i++) {
        fread(mat[i], sizeof(double), tam, fich);
    }
}

//Funcion para calcular el vector resultado despues de la division
void dividirVec(double *vec, double norma, int tam)
{
	for(int i = 0; i < tam; i++)
	{
		vec[i] = vec[i] / norma; 
	}
}

//Funcion para calcular los vectores resultados de la multiplicacion en cada proceso
void calcularVec(double *vec2, double *vec, double **mat, int tam, int filRep)
{
	for(int i = 0; i < filRep; i++)
	{
		vec2[i] = 0.0;
		for(int j = 0; j < tam; j++)
		{
			vec2[i] += mat[i][j] * vec[j];
		}
	}
}

//Funcion para calcular la norma 
double calcularNorma(double *vec, int tam) {
    double norma = 0.0;
    for (int i = 0; i < tam; i++) 
	{
        norma += vec[i] * vec[i];
    }
    return sqrt(norma);
}

int main(int argc, char *argv[]) 
{
	double **mat, *vec1, *vec2, *vec3, norma;
	double tiempo_inicio, tiempo_final;
    int myrank, nproces, tam, numIteraciones, filRep, filaInic, *desplazamiento, *aRecibir;
	
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproces);
	MPI_Status status;
	if(myrank == 0)
	{
		
		char *nom_fich = argv[1];
		tam = atoi(argv[2]);
		numIteraciones = atoi(argv[3]);
		FILE *fich;
		tiempo_inicio = MPI_Wtime();
		//Reservo memoria de los arrays de desplazamiento y datos recibidos para el Gatherv
		desplazamiento = (int*)malloc(nproces * sizeof(int));
		aRecibir = (int*)malloc(nproces * sizeof(int));
		
		// printf("%d",tam);
		
		//Reservo memoria
		mat = (double**)malloc(tam * sizeof(double*));
		for(int i = 0; i < tam; i++)
		{
			mat[i] = (double*)malloc(tam * sizeof(double));
		}
		fich = fopen(nom_fich, "rb");
		//Si el fichero no existe me crea una matriz con 1 en la diagonal y lo demas aleatorio y que sea menor 0.1 y mayor que -0.1
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
			//Si si que existe el fichero me crea una matriz con los datos del fichero
			crearMatriz(mat, fich, tam);
		}
		
		//Calculo los arrays para el Gatherv y los datos a enviar a cada proceso
		filRep = tam / nproces;
		int resto = tam % nproces;
		aRecibir[0] = filRep;
		desplazamiento[0] = 0;
		
		//Envio los datos a cada proceso
		for(int i = 1; i < nproces; i++)
		{
			int filAux = filRep;
			filaInic = filRep * i;
			desplazamiento[i] = filaInic;
			if(i == nproces - 1)
			{
				filAux += resto;
			}else{
				filAux = filRep;
			}
			aRecibir[i] = filAux;
			for(int j = 0; j < filAux; j++)
			{
				MPI_Send(mat[filaInic + j], tam, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
			MPI_Send(&filAux, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
		}
	}
	//Realizo un Bcast sobre el numero de iteraciones y el tamaño
	MPI_Bcast(&tam, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numIteraciones, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(myrank != 0)
	{
		MPI_Recv(&filRep, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		
		mat = (double**)malloc(filRep * sizeof(double*));
		for(int i = 0; i < filRep; i++)
		{
			mat[i] = (double*)malloc(tam * sizeof(double));
		}
		
		for(int i = 0; i < filRep; i++)
		{
			MPI_Recv(mat[i], tam, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}
	}
	
	//Reservo memoria en todos los procesos para los vectores que voy a usar para los calculos
	vec1 = (double*)malloc(tam * sizeof(double));
	vec2 = (double*)calloc(filRep, sizeof(double));
	
	if(myrank == 0)
	{
		//Reservo memoria para todos 
		vec3 = (double*)malloc(tam * sizeof(double));
		
		//Inicializo el vector 1 como el vector unitario
		for(int i = 0; i < tam; i++)
		{
			vec1[i] = 1;
		}
	}
	
		
	for(int i = 0; i < numIteraciones; i++)
	{
		
		if(myrank == 0)
		{
			// for(int j = 0; j < tam; j++)
			// {
				// printf("%f, %f\n", vec1[j], vec3[j]);
			// }
				
			//Calculo la norma en el proceso 0
			norma = calcularNorma(vec1, tam);
			//Si estamos en una iteracion diferente a la 0 igualo el vec1 al vector resultado del anterior iteracion
			if(i != 0)
			{
				for (int j = 0; j < tam; j++) 
				{
					vec1[j] = vec3[j];
				}
			}
		}
		// printf("%d: %f\n",myrank, vec1[3]);
		
		//Envio el vec1 a todos los procesos con Bcast
		MPI_Bcast(vec1, tam, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//Calculo el vec2 en todos los procesos
		calcularVec(vec2, vec1, mat, tam, filRep);
		//Recivo todos los vec2 en el vector 0 y los recibo en el vec3 con Gatherv
		MPI_Gatherv(vec2, filRep, MPI_DOUBLE, vec3, aRecibir, desplazamiento, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if(myrank == 0)
		{
			//Si estamos e una iteracion diferente de 0 realizo la division
			if(i != 0)
			{
				dividirVec(vec3, norma, tam);
			}
			printf("%d: %f\n",i, norma);
			// for(int j = 0; j < tam; j++){
				// printf("%f, %f\n", vec1[j], vec3[j]);
			// }
		}
	}
	
	//Paro el contador y libero memoria en los procesos correspondientes
	if (myrank == 0) {
		tiempo_final = MPI_Wtime();
		printf("%f", tiempo_final-tiempo_inicio);
		
		free(vec3);
		free(aRecibir);
		free(desplazamiento);
	}
	for (int i = 0; i < filRep; i++) 
	{
		free(mat[i]);
	}
	free(mat);
	free(vec1);
	free(vec2);
	
	MPI_Finalize();
    return EXIT_SUCCESS;
} 