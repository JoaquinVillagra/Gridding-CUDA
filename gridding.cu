#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <time.h>
#include <math.h>
#include <cuda_runtime.h>
#define PI 3.14159265358979323846
#define FactorArcosegRad 0.00000484814
#define BLOQUESIZE 32

clock_t timestart, timeend;

/**
@brief Función que transforma un valor en arco segundo a radianes
@param deltax: Valor numérico a transformar
@returns Valor correspondiente a la entrada en radianes */
float arcoseg_radian(float deltax){
	return FactorArcosegRad*deltax;
}

/**
@brief Función que lee el archivo de entrada
@param archivo: puntero al archivo a leer
@param tamano: Numero de visibilidades del archivo a leer
@returns  */
double* readFile(FILE* archivo, int tamano){
	double* elementos =(double*) malloc(sizeof(double)*4*tamano);
	fread(elementos, tamano*4, sizeof(double), archivo);
	return elementos;
}


/**
@brief Función ejecuta el proceso de gridding
@param U: Valores de la coordenada U en el plano de Fourier
@param V: Valores de la coordenada V en el plano de Fourier
@param R: Valores reales de la visibilidad en el plano de Fourier
@param I: Valores imaginarios la visibilidad en el plano de Fourier
@param num_datos: Cantidad de visibilidades ingresadas o dimensión de los vectores anteriores
@param tamano: Lado de la matriz a construir, si tamano es 512 se construye una matriz de 512X512
@param V: Valores de la coordenada V en el plano de Fourier
@param deltaU: Valor delta necesario para determinar la vecindad de cada pixel de la grilla regular
@param r: vector de valores reales de la salida del proceso de gridding
@param k: vector de valores imaginarios de la salida del proceso de gridding
@returns  */
__global__ void gridding_process(float *U, float *V, float *R, float *I, int num_datos, int tamano, float deltaU, float *r, float *k)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i<num_datos)
	{
		float x, y, modx, mody;
		x = U[i]/deltaU+tamano/2;
		y = V[i]/deltaU+tamano/2;
		modx = U[i] - x*deltaU;
		mody = V[i] - y*deltaU;

		if(modx>deltaU/2){	
			x+=1;
		}

		if (mody>deltaU/2)
		{
			y+=1;
		}

		if ((int)x<tamano && (int)y<tamano)
		{
			atomicAdd(&r[(int)y*tamano+(int)x], R[i]);
			atomicAdd(&k[(int)y*tamano+(int)x], I[i]);
		}
	}
}

__host__ unsigned long upper_power_of_two(unsigned long v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;

}

int main(int argc, char * const argv[])
{
	int tamano;//tamaño de imagen
	int numdatos;//número de pasos
	float deltaX_arcoseg, deltaX_radian;
	float deltaU; 
	char* archivo_entrada=NULL;
	char* archivo_salida=NULL;
	int i, c;
	opterr = 0;
	while ((c = getopt (argc, argv, "i:z:d:N:o:")) != -1)
		switch (c)
			{
			case 'i':
				archivo_entrada = optarg;
				break;
			case 'z':
				numdatos = atoi(optarg);
				break;
			case 'd':
				deltaX_arcoseg = atof(optarg);
				break;
			case 'N':
				tamano = atoi(optarg);
				break;
			case 'o':
				archivo_salida = optarg;
				break;
			case '?':
				if (optopt == 'i' ||optopt == 'z' ||optopt == 'd'||optopt == 'N' ||optopt == 'o')
					fprintf (stderr, "Opcion -%c requiere un argumento.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Opcion desconocida `-%c'.\n", optopt);
				else
					fprintf (stderr,
									"Carater opcion desconocido `\\x%x'.\n",
									optopt);
				return 1;
			default:
				abort ();
			}
	/**
		Comprobación de Inputs
			- Valores mayores que cero
			- Cadenas no nulas
	**/
	if(tamano<=0){
		printf("El parametro -N debe existir y ser mayor que 0\n");
		exit(1);
	}
	if(numdatos==0){
		printf("El parametro -z debe existir y ser mayor que 0\n");
		exit(1);
	}
	if(deltaX_arcoseg==0){
		printf("El parametro -d debe existir y ser mayor que 0\n");
		exit(1);
	}
	if(archivo_entrada==NULL){
		printf("Debe especificarse un archivo de entrada\n");
	}
	if(archivo_salida==NULL){
		printf("Debe especificarse un archivo de salida\n");
	}
	//Transformacion de unidades necesaria para calcular delta U
	deltaX_radian = arcoseg_radian(deltaX_arcoseg);

	//Determina delta U/V a utilizar
	deltaU = 1/(tamano*deltaX_radian);

	//Medición de tiempo de computo
	timestart = clock(); 

	//Lectura de entrada
	FILE *entrada = fopen(archivo_entrada,"r");
	double* data = readFile(entrada,numdatos);
	fclose(entrada);

	//Creando arrays para coordenada X, Y, R e I
	float *X = (float*)malloc(sizeof(float)*numdatos); 
	float *Y = (float*)malloc(sizeof(float)*numdatos); 
	float *R = (float*)malloc(sizeof(float)*numdatos); 
	float *I = (float*)malloc(sizeof(float)*numdatos);	
	//Quizas necesite dos vectores adicionales para el gridding [matrices desenroyadas]
	float *r = (float*)malloc(sizeof(float)*tamano*tamano);
	float *k = (float*)malloc(sizeof(float)*tamano*tamano);
	//Se asigan los valores correspondientes de la lectura
	for (i = 0; i < numdatos; i++)
	{
		X[i] = (float)data[i];
		Y[i] = (float)data[i+numdatos];
		R[i] = (float)data[i+2*numdatos];
		I[i] = (float)data[i+3*numdatos];

	}
	for (i = 0; i < tamano*tamano; ++i)
	{
		r[i] = 0;
		k[i] = 0;
	}
	//se declaran las variables CUDA
	float *C_X;
	float *C_Y;
	float *C_R;
	float *C_I;
	float *C_r;
	float *C_k;
	//Se reserva memoria CUDA
	cudaMalloc( (void**)&C_X, numdatos*sizeof(float)); 
	cudaMalloc( (void**)&C_Y, numdatos*sizeof(float)); 
	cudaMalloc( (void**)&C_R, numdatos*sizeof(float)); 
	cudaMalloc( (void**)&C_I, numdatos*sizeof(float)); 
	cudaMalloc( (void**)&C_r, tamano*tamano*sizeof(float)); 
	cudaMalloc( (void**)&C_k, tamano*tamano*sizeof(float)); 
	//se copia la matriz iniciada en las matrices de trabajo en memoria global GPU
	cudaMemcpy( C_X, X, numdatos*sizeof(float), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_Y, Y, numdatos*sizeof(float), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_R, R, numdatos*sizeof(float), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_I, I, numdatos*sizeof(float), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_r, r, tamano*tamano*sizeof(float), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_k, k, tamano*tamano*sizeof(float), cudaMemcpyHostToDevice); 
	//determino dimension para el kernel
	long data_size_2 = upper_power_of_two(numdatos);
	//Se declaran las dimenciones
	dim3 dimBlock(BLOQUESIZE, 1);
	dim3 dimGrid(data_size_2/BLOQUESIZE, 1);
	//se ejecuta el kernel en la GPU
	//printf("%d - %d - %d\n", numdatos, kernel_size, kernel_size/BLOQUESIZE);

	gridding_process<<<dimGrid, dimBlock>>>(C_X, C_Y, C_R, C_I, numdatos, tamano, deltaU, C_r, C_k);
	//se espera a que terminen
	cudaDeviceSynchronize();
	//se obtiene la memoria de regreso
	cudaMemcpy( r, C_r, tamano*tamano*sizeof(float), cudaMemcpyDeviceToHost); 
	cudaMemcpy( k, C_k, tamano*tamano*sizeof(float), cudaMemcpyDeviceToHost); 
	//se libera la memoria global CUDA para que pueda ser usada por otro proceso
	cudaFree( C_X );
	cudaFree( C_Y );
	cudaFree( C_R );
	cudaFree( C_I );
	cudaFree( C_r );
	cudaFree( C_k );
	//Se imprime salida
	FILE *f = fopen(strcat(archivo_salida, "real.raw"),"wb");
	FILE *g = fopen(strcat(archivo_salida, "img.raw"),"wb");
	fwrite(r, tamano*tamano, sizeof(float),f);
	fwrite(k, tamano*tamano, sizeof(float),g);
	fclose(f);
	fclose(g);
	//Se mide el tiempo utilizado
	timeend = clock();
	printf("Total = %f\n", (double) (timeend-timestart)/(double)CLOCKS_PER_SEC);
	return EXIT_SUCCESS;
}