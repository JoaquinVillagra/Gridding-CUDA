#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <time.h>


clock_t timestart, timeend;

/**
@brief Función que transforma un valor en arco segundo a radianes
@param deltax: Valor numérico a transformar
@returns Valor correspondiente a la entrada en radianes */
float arcoseg_radian(float deltax){
	return 0.000004848140*deltax;
}


int main(int argc, char * const argv[])
{
	int tamano;//tamaño de imagen
	int numdatos;//número de pasos
	float deltaX_arcoseg, deltaX_radian;
	float deltaU; 
	int threads_num_X;//numero de hebras
	int threads_num_Y;//numero de hebras
	char* archivo_entrada=NULL;
	char* archivo_salida=NULL;
	int c;

	opterr = 0;
	//GETOPT
	while ((c = getopt (argc, argv, "i:z:d:N:o")) != -1)
		switch (c)
			{
			case 'i':
				archivo_entrada = optarg;
				break;
			case 'z':
				numdatos = atoi(optarg);
				break;
			case 'd':
				deltaX_arcoseg = atoi(optarg);
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
	if(tamano==0){
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
	deltaX_radian = arcoseg_radian(delta_arcoseg);

	timestart = clock(); 

	//se calculan las dimenciones
	const long n = Grilla;
	const long N = n*n; 

	float* b= (float*) malloc(sizeof(float)*N);
	long i,j;
	//se inicia la matriz, es decir H(t=0)
	for(i =0;i<n;i++){
		for(j=0;j<n;j++){
			if(i>0.4*n && i<0.6*n && j>0.4*n && j<0.6*n)
				b[j+n*i]=20;
			else
				b[j+n*i]=0;
		}
	}
	
	







	//se declaran las variables CUDA
	float *H;
	float *H_1;
	float *H_2;
	float *H_t;
	//Se reserva memoria CUDA
	cudaMalloc( (void**)&H, N*sizeof(float) ); 
	cudaMalloc( (void**)&H_1, N*sizeof(float) ); 
	cudaMalloc( (void**)&H_2, N*sizeof(float) ); 
	cudaMalloc( (void**)&H_t, N*sizeof(float) ); 
	//se copia la matriz iniciada en las matrices de trabajo en memoria global GPU
	cudaMemcpy( H, b, N*sizeof(float), cudaMemcpyHostToDevice ); 
	cudaMemcpy( H_1, b, N*sizeof(float), cudaMemcpyHostToDevice ); 
	cudaMemcpy( H_2, b, N*sizeof(float), cudaMemcpyHostToDevice ); 
	cudaMemcpy( H_t, b, N*sizeof(float), cudaMemcpyHostToDevice ); 

	//Se declaran las dimenciones
	dim3 dimBlock( threads_num_X, threads_num_Y );
	dim3 dimGrid( 1, 1 );
		
	for(i=0;i<Steps;i++){
	//se ejecuta el kernel en la GPU
	schrodinger<<<dimGrid, dimBlock>>>(H,H_1,H_2,H_t,i,Steps,n);
	//se espera a que terminen
	cudaDeviceSynchronize();
	}
	//se obtiene la memoria de regreso
	cudaMemcpy( b, H_t, N*sizeof(float), cudaMemcpyDeviceToHost ); 
	//se libera la memoria global CUDA para que pueda ser usada por otro proceso
	cudaFree( H );
	cudaFree( H_1 );
	cudaFree( H_2 );
	cudaFree( H_t );
	//Se escribe el archivo
	FILE *f = fopen(salida_arg,"wb");
	fwrite(b,N, sizeof(float),f);

	timeend = clock(); // registramos el tiempo hasta el final
	printf("Total = %f\n", (double) (timeend-timestart)/(double)CLOCKS_PER_SEC);
	return EXIT_SUCCESS;
}