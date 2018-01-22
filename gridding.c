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

__global__ 
void gridding_process(){

}

int main(int argc, char * const argv[])
{
	int tamano;//tamaño de imagen
	int numdatos;//número de pasos
	double deltaX_arcoseg, deltaX_radian;
	double deltaU; 
	char* archivo_entrada=NULL;
	char* archivo_salida=NULL;
	int i, j, c;
	
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

	double x, y, modx, mody;
	//Crea matrices con memoria dinamica
	double **matriz_real = (double**)malloc(sizeof(double*)*tamano);
	double **matriz_i = (double**)malloc(sizeof(double*)*tamano);
	for (i = 0; i < tamano; ++i)
	{
		matriz_real[i] = (double*)malloc(sizeof(double)*tamano);
		matriz_i[i] = (double*)malloc(sizeof(double)*tamano);
	}
	//Inicializa matrices en cero
	for (i = 0; i < tamano; ++i)
	{
		for (int j = 0; j < tamano; ++j)
		{
			matriz_real[i][j]=0;
			matriz_i[i][j]=0;
		}
	}
	//Creando arrays para coordenada X, Y, R e I
	double *X = (double*)malloc(sizeof(double)*numdatos); 
	double *Y = (double*)malloc(sizeof(double)*numdatos); 
	double *R = (double*)malloc(sizeof(double)*numdatos); 
	double *I = (double*)malloc(sizeof(double)*numdatos);	
	//Quizas necesite dos vectores adicionales para el gridding [matrices desenroyadas]
	double *r = (double*)malloc(sizeof(double)*tamano*tamano);
	double *k = (double*)malloc(sizeof(double)*tamano*tamano);
	//Se asigan los valores correspondientes de la lectura
	for (i = 0; i < numdatos; i++)
	{
		X[i] = data[i];
		Y[i] = data[i+numdatos];
		R[i] = data[i+2*numdatos];
		I[i] = data[i+3*numdatos];

	}
	//se declaran las variables CUDA
	double *C_X;
	double *C_Y;
	double *C_R;
	double *C_I;
	//Se reserva memoria CUDA
	cudaMalloc( (void*)&C_X, numdatos*sizeof(double)); 
	cudaMalloc( (void*)&C_Y, numdatos*sizeof(double)); 
	cudaMalloc( (void*)&C_R, numdatos*sizeof(double)); 
	cudaMalloc( (void*)&C_I, numdatos*sizeof(double)); 
	//se copia la matriz iniciada en las matrices de trabajo en memoria global GPU
	cudaMemcpy( C_X, X, numdatos*sizeof(double), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_Y, Y, numdatos*sizeof(double), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_R, R, numdatos*sizeof(double), cudaMemcpyHostToDevice); 
	cudaMemcpy( C_I, I, numdatos*sizeof(double), cudaMemcpyHostToDevice); 

	//Se declaran las dimenciones
	//dim3 dimBlock( ? , ? );
	//dim3 dimGrid( 1, 1 );

	for (i = 0; i < numdatos; i++)
	{
		//printf("[%lf,%lf,%lf,%lf]\n",data[i],data[numdatos+i],data[2*numdatos+i],data[3*numdatos+i] );
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
		//printf("x es: %d e y es: %d\n", (int)x, (int)y);
		matriz_real[(int)y][(int)x] += R[i];
		matriz_i[(int)y][(int)x] += I[i];
	}
	fclose(entrada);
	//printf("Delta U %lf\n", deltaU );
	FILE *f = fopen("salida_real","wb");
	FILE *g = fopen("salida_imaginaria","wb");
	for (i = 0; i < tamano; i++)
	{
		fwrite(matriz_real[i],tamano, sizeof(double),f);
		fwrite(matriz_i[i],tamano, sizeof(double),g);
	}
	timeend = clock(); // registramos el tiempo hasta el final
	printf("Total = %f\n", (double) (timeend-timestart)/(double)CLOCKS_PER_SEC);
	return EXIT_SUCCESS;
}