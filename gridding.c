#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <time.h>
#define PI 3.14159265358979323846
#define FactorArcosegRad 0.00000484814

clock_t timestart, timeend;

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
	double* elementos = malloc(sizeof(double)*4*tamano);
	fread(elementos, tamano*4, sizeof(double), archivo);
	return elementos;
}

int main(int argc, char * const argv[])
{
	printf("0!!!!\n");
	int tamano;//tamaño de imagen
	int numdatos;//número de pasos
	double deltaX_arcoseg, deltaX_radian;
	double deltaU; 
	char* archivo_entrada=NULL;
	char* archivo_salida=NULL;
	int i, j, c;

	printf("1!!!!\n");
	opterr = 0;
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
	printf("ACAAA1\n");
	//Transformacion de unidades necesaria para calcular delta U
	deltaX_radian = arcoseg_radian(deltaX_arcoseg);
	printf("Delta X en radian %lf\n",  deltaX_radian);

	printf("ACAAA2\n");
	//Determina delta U/V a utilizar
	deltaU = 1/(tamano*deltaX_radian);

	//Medición de tiempo de computo
	timestart = clock(); 

	//Lectura de entrada
	FILE *entrada = fopen(archivo_entrada,"r");
	//double* data = readFile(entrada, 4*numdatos);

	//double *regularReal = malloc(sizeof(double)*tamano*tamano);
	//double *regularImaginario = malloc(sizeof(double)*tamano*tamano);

	double x, y, modx, mody;
	double matriz_real[2048][2048];
	double matriz_i[2048][2048];
	for (int i = 0; i < 2048; ++i)
	{
		for (int j = 0; j < 2048; ++j)
		{
			matriz_real[i][j]=0;
			matriz_i[i][j]=0;
		}
	}
	printf("AQUI despues de set(0)\n");
	double* data = readFile(entrada,numdatos);
	for (int i = 0; i < numdatos; i++)
	{
		printf("[%lf,%lf,%lf,%lf]\n",data[i],data[numdatos+i],data[2*numdatos+i],data[3*numdatos+i] );
		x = data[i]/deltaU+tamano/2;
		y = data[numdatos+i]/deltaU+tamano/2;
		modx = data[i] - x*deltaU;
		mody = data[numdatos+i] - y*deltaU;
		if(modx>deltaU/2){	
			x+=1;
		}
		if (mody>deltaU/2)
		{
			y+=1;
		}
		printf("x es: %d e y es: %d\n", (int)x, (int)y);
		//matriz_real[(int)y][(int)x] += data[2*numdatos+i];
		//matriz_i[(int)y][(int)x] += data[3*numdatos+i];
	}
	fclose(entrada);
	//printf("Delta U %lf\n", deltaU );
	FILE *f = fopen("salida_real","wb");
	FILE *g = fopen("salida_imaginaria","wb");
	for (i = 0; i < tamano; i++)
	{
		fwrite(matriz_real[i],2048, sizeof(double),f);
		fwrite(matriz_i[i],2048, sizeof(double),g);
	}
	printf("OK!	\n");
	return(0);
	
/*


	//se declaran las variables CUDA
	double *H;
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
	<<<dimGrid, dimBlock>>>(H,H_1,H_2,H_t,i,Steps,n);
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
*/
}