#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

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
	//Inicializo todas las variables
   int iam = 0, np = atoi(argv[1]);
   char *nom = argv[2];
   
   double inicio, final, total;
   
   int fil = atoi(argv[3]), col = atoi(argv[4]), filExt = fil + 2, colExt = col + 2, min = INT_MAX,
	max = INT_MIN, *histograma, tamRep = (fil) / np, hisRep = 256/np;
	
	//Reservo memoria
   uint8_t **mat, **matFil;
   histograma = (int*)calloc(256, sizeof(int));
   
   int *vecDesp = (int*)malloc(np * sizeof(int));   
   int *vecProc = (int*)malloc(np * sizeof(int));
   
   //Creo el vector de desplazamientos y elementos procesados.
   int filRep = fil / np;
   int resto = fil % np;
   for (int i = 0; i < np; i++) 
   {
	   vecDesp[i] = (i == 0) ? 1 : vecDesp[i - 1] + vecProc[i - 1];
	   vecProc[i] = (i == np - 1) ? filRep + resto : filRep;
   }
   
   mat = (uint8_t**)malloc(filExt * sizeof(uint8_t*));
   for(int i = 0; i < filExt; i++)
   {
	   mat[i] = (uint8_t*)malloc(colExt * sizeof(uint8_t));
   }
   
   matFil = (uint8_t**)malloc(fil * sizeof(uint8_t*));
   for(int i = 0; i < fil; i++)
   {
	   matFil[i] = (uint8_t*)malloc(col * sizeof(uint8_t));
   }
   
   FILE *fich = fopen(nom, "rb");
   for(int i = 1; i < filExt - 1; i++)
   {
	   fread(&mat[i][1], sizeof(uint8_t), col, fich);
   }
   
   //Calculo la matriz extendidad
   calcMatExt(mat, filExt, colExt);
   
   inicio = omp_get_wtime();
   //printf("%d\n", tamRep);
   
#pragma omp parallel num_threads(np) shared(mat, filExt, colExt)
{
	//Genero un histograma local en cada proceso
	int *histoLocal = (int*)calloc(256, sizeof(int));
	
	//Con un for paralelo calculo el histograma local de cada proceso
	#pragma omp for schedule(static, tamRep)
	for(int i = 1; i < filExt - 1; i++)
	{	
		for(int j = 1; j < colExt - 1; j++)
		{
			histoLocal[mat[i][j]]++;
		}
	}
	
	//Actualizo los valores en el histograma local con un critical
	#pragma omp critical
	{
		for(int i = 0; i < 256; i++)
		{
			histograma[i] += histoLocal[i];
		}
	}
	//Añado una barrera para esperar a que terminen todos los procesos
	#pragma omp barrier
	
	//Con otro for paralelo calculo el maximo y minimo del histograma global para repartirlo entre todos los procesos.
	#pragma omp for schedule(static, hisRep) reduction(max : max) reduction(min : min)
	for(int i = 0; i < 256; i++)
	{
		if(max < histograma[i])
		{
			max = histograma[i];
		}
		if(min > histograma[i] && histograma[i] != 0)
		{
			min = histograma[i];
		}
	}
	free(histoLocal);
}//parallel
final = omp_get_wtime();
total = final - inicio;
printf("%f \n", total);

inicio = omp_get_wtime();
#pragma omp parallel num_threads(np) shared(mat, matFil, filExt, colExt, vecDesp, vecProc) private(iam)
{
	
	//Calculo la mediana usando el vector de elementos repartidos para cada elemento y donde inician
	iam = omp_get_thread_num();
	for (int i = vecDesp[iam]; i < vecDesp[iam] + vecProc[iam]; i++) {
		for (int j = 1; j < colExt - 1; j++) {
		uint8_t ventana[9] = {
			mat[i - 1][j - 1], mat[i - 1][j], mat[i - 1][j + 1],
			mat[i][j - 1],     mat[i][j],     mat[i][j + 1],
			mat[i + 1][j - 1], mat[i + 1][j], mat[i + 1][j + 1]
		};
		matFil[i - 1][j - 1] = calcMediana(ventana);
		}
	}
}

final = omp_get_wtime();
total = final - inicio;
printf("%f \n", total);

//Escribo el resultado en dos ficheros uno el histograma y otro la foto filtrada
FILE *fich2 = fopen("resultado.raw", "wb");

for(int i = 0; i < fil; i++)
{
	fwrite(matFil[i], sizeof(uint8_t), col	, fich2);
}

FILE *fich3 = fopen("histograma.txt", "w");
fprintf( fich3, "max : %d, min : %d\n", max, min);
for(int i = 0; i < 256; i++)
{
	fprintf(fich3, "%d: %d\n", i, histograma[i]);
}

//Libero memoria
for(int i = 0; i < filExt; i++)
{
	free(mat[i]);
}

for(int i = 0; i < fil; i++)
{
	free(matFil[i]);
}

free(vecDesp);
free(vecProc);

// for(int i = 0; i < np; i++)
// {
	// printf("%d, %d\n", vecDesp[i], vecProc[i]);
// }

}