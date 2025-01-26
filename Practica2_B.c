#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <math.h>


//Funcion para crear una matriz aleatoria si no encuentra el fichero especificado
FILE *crearMatAleatoria(double **mat, int tam)
{
	FILE *fich = fopen("random.bin", "w");
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
		fwrite(mat[i], sizeof(double), tam, fich);
	}
	return fich;
}

//Funcion que hace la multuplicacion de la matriz por el vector teniendo en cuenta el inicio y el final de cada proceso
void calcularVec(double *vec2, double *vec1, double **mat, int tam, int filRep, int inicio)
{
	for(int i = inicio; i < filRep + inicio; i++)
	{
		vec2[i] = 0.0;
		for(int j = 0; j < tam; j++)
		{
			vec2[i] += mat[i][j] * vec1[j];
		}
		
	}
}

//Funcion para leer un archivo y almacenarlo en la matriz especificada
void crearMatriz(double** mat, FILE *fich, int tam) {
    for (int i = 0; i < tam; i++) {
        fread(mat[i], sizeof(double), tam, fich);
    }
}

int main(int argc, char *argv[])
{
	//Inicializo las variables
	int iam, np = atoi(argv[1]), tam = atoi(argv[3]), m = atoi(argv[4]), vecRep[np], vecInic[np];
	char *nom = argv[2];
	double **mat, *vec1, *vec2, norma = 0, inicio, final, total;
	
	FILE *fich = fopen(nom, "rb");
	
	//Reservo la memoria de la matriz y los vectores utilizados
	mat = (double**)malloc(tam * sizeof(double*));
	for(int i = 0; i < tam; i++)
	{
		mat[i] = (double*)malloc(tam * sizeof(double));
	}
	
	vec1 = (double*)malloc(tam * sizeof(double));
	for(int i = 0; i < tam; i++)
	{
		vec1[i] = 1;
	}
	vec2 = (double*)malloc(tam * sizeof(double));
	
	//Realizo la funcion de leer el archivo y si es null le asigno a la matriz valores aleatorios
	if(fich == NULL)
	{
		FILE *fich = crearMatAleatoria(mat, tam);
	}else{
		//Si si que existe el fichero me crea una matriz con los datos del fichero
		crearMatriz(mat, fich, tam);
	}
	
	//Creo los vectores de reparticion de trabajo
	int tamRep = tam/np;
	for(int i = 0; i < np; i++)
	{
		vecRep[i] = (i == 0) ? (tamRep) + (tam % np) : tamRep;
		vecInic[i] = (i == 0) ? 0 : vecRep[i - 1] + vecInic[i - 1];
	}
	
	//Empiezo las iteraciones
	if(m > 0)
	{
		inicio = omp_get_wtime();
		//Primera iteracion sin automatizadores utilizando los vectores de reparticion de trabajo
		#pragma omp parallel num_threads(np) shared(mat, vec1, vec2, tam, vecRep, vecInic) private(iam)
		{
			iam = omp_get_thread_num();
			calcularVec(vec2, vec1, mat, tam, vecRep[iam], vecInic[iam]);
		}
		//Las demas iteraciones se realizan con automatizadores en esta region paralela
		#pragma omp parallel num_threads(np) shared(mat, vec1, vec2, tam, norma)
		{
			for(int l = 1; l < m; l++)
			{
				//Asigno la norma a 0 en un solo proceso
				#pragma omp single
				{
					norma = 0;
				}
				
				#pragma omp for schedule(static, tamRep) reduction(+ : norma)
				for(int i = 0; i < tam; i++)
				{
					norma += (vec1[i] * vec1[i]);
				}
				#pragma omp single
				{
					norma = sqrt(norma);
					// printf("%.6e\n", norma);
					double *temp = vec1;
					vec1 = vec2;
					vec2 = temp;
				}
				
				#pragma omp for schedule(static, tamRep)
				for(int i = 0; i < tam; i++)
				{
					vec2[i] = 0.0;
					for(int j = 0; j < tam; j++)
					{
						vec2[i] += mat[i][j] * vec1[j];
					}
					vec2[i] = vec2[i] / norma;
				}
			}
		}
		final = omp_get_wtime();
		total = final - inicio;
		printf("%f \n", total);
	}
	for (int i = 0; i < tam; i++) 
	{
		free(mat[i]);
	}
	free(mat);
	free(vec1);
	free(vec2);
	
}