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

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* a, double b) { return b; }
#endif

/**
@brief Función que transforma un valor en arco segundo a radianes
@param deltax: Valor numérico a transformar
@returns Valor correspondiente a la entrada en radianes */
double arcoseg_radian(double deltax){
	return FactorArcosegRad*deltax;
}

/**
@brief Función que lee el archivo de entrada
@param archivo: puntero al archivo a leer
@param archivo: puntero al archivo a leer
@returns  */
double* readFile(FILE* archivo, int tamano){
	double* elementos =(double*) malloc(sizeof(double)*4*tamano);
	fread(elementos, tamano*4, sizeof(double), archivo);
	return elementos;
}

__global__ 
void gridding_process(float *X, float *Y, float *R, float *I, int num_datos, int tamano, float deltaU, float *r, float *k)
{
	long i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i<num_datos)
	{
		double x, y, modx, mody;
		x = X[i]/deltaU+tamano/2;
		y = Y[i]/deltaU+tamano/2;
		modx = X[i] - x*deltaU;
		mody = Y[i] - y*deltaU;
		if(modx>deltaU/2){	
			x+=1;
		}
		if (mody>deltaU/2)
		{
			y+=1;
		}
		atomicAdd(&r[(int)y*tamano+(int)x], R[i]);
		atomicAdd(&k[(int)y*tamano+(int)x], I[i]);
	}
}

int main(int argc, char * const argv[])
{
	int tamano;//tamaño de imagen
	int numdatos;//número de pasos
	double deltaX_arcoseg, deltaX_radian;
	double deltaU; 
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
		printf("El parametro -N debe estár y ser mayor que 0\n");
		exit(1);
	}
	if(numdatos==0){
		printf("El parametro -z debe estár y ser mayor que 0\n");
		exit(1);
	}
	if(deltaX_arcoseg==0){
		printf("El parametro -d debe estár y ser mayor que 0\n");
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
	//determino dimension para el kernel
	long kernel_size, raiz;
	raiz = (int)sqrt(numdatos);;
	if(raiz % 2 != 0)
		kernel_size = (raiz + 2)*(raiz + 2);
	kernel_size = (raiz+1)*(raiz+1); 
	
	//Se declaran las dimenciones
	dim3 dimBlock(BLOQUESIZE, 1);
	dim3 dimGrid(kernel_size, 1);
	//se ejecuta el kernel en la GPU
	gridding_process<<<dimGrid, dimBlock>>>(C_X, C_Y, C_R, C_Y, numdatos, tamano, (float)deltaU, C_r, C_k);
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
	FILE *f = fopen("salida_real","wb");
	FILE *g = fopen("salida_imaginaria","wb");

	fwrite(r, tamano*tamano, sizeof(float),f);
	fwrite(k, tamano*tamano, sizeof(float),g);

	timeend = clock(); // registramos el tiempo hasta el final
	printf("Total = %f\n", (double) (timeend-timestart)/(double)CLOCKS_PER_SEC);
	return EXIT_SUCCESS;
}