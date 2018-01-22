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
	double **matriz_real = (double**)malloc(sizeof(double*)*tamano);
	double **matriz_i = (double**)malloc(sizeof(double*)*tamano);
	for (i = 0; i < tamano; ++i)
	{
		matriz_real[i] = (double*)malloc(sizeof(double)*tamano);
		matriz_i[i] = (double*)malloc(sizeof(double)*tamano);
	}
	for (int i = 0; i < tamano; ++i)
	{
		for (int j = 0; j < tamano; ++j)
		{
			matriz_real[i][j]=0;
			matriz_i[i][j]=0;
		}
	}
	for (int i = 0; i < numdatos; i++)
	{
		//printf("[%lf,%lf,%lf,%lf]\n",data[i],data[numdatos+i],data[2*numdatos+i],data[3*numdatos+i] );
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
		//printf("x es: %d e y es: %d\n", (int)x, (int)y);
		matriz_real[(int)y][(int)x] += data[2*numdatos+i];
		matriz_i[(int)y][(int)x] += data[3*numdatos+i];
	}
	fclose(entrada);
	//printf("Delta U %lf\n", deltaU );
	FILE *f = fopen("salida_real","wb");
	FILE *g = fopen("salida_imaginaria","wb");
	for (i = 0; i < tamano; i++)
	{
		fwrite(matriz_real[i],tamano, sizeof(double),f);
		fwrite(matriz_i[i],tamno, sizeof(double),g);
	}
	timeend = clock(); // registramos el tiempo hasta el final
	printf("Total = %f\n", (double) (timeend-timestart)/(double)CLOCKS_PER_SEC);
	return EXIT_SUCCESS;
}