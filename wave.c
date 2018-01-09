#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/times.h>
#include <time.h>
 

#define K1 0.0025
#define K2 0.00125

clock_t timestart, timeend;
 
__global__ 
void schrodinger(float *H,float *H_1,float *H_2,float *H_t,int index,int t,long N) 
{
	// memoria para darle orden a las memorias globales, sirve para almacenar los 3 buffer requeridos
	__shared__ float* Hs[3];
	long i,j;
	
	//el hilo 0,0 es quien da el orden de los buffer para que estos roten
	if(threadIdx.x==0 && threadIdx.y==0){
		if(index%3==0){
			Hs[0]=H;
			Hs[1]=H_1;
			Hs[2]=H_2;
		}else if(index%3==1){
			Hs[0]=H_2;
			Hs[1]=H;
			Hs[2]=H_1;
		}else{
			Hs[0]=H_1;
			Hs[1]=H_2;
			Hs[2]=H;
		}
	}
	//mientras se realiza eso las demas hebras esperan
	__syncthreads();
	//cada hebra calcula una parte de la matriz, de modo que sea proporcionalmente muy similar la cantidad de trabajo
	for (j = (N/gridDim.x)*blockIdx.x+threadIdx.y; j < N-1; j+=blockDim.y)//Y
	{
		if(j == 0) continue;
		for (i = threadIdx.x; i < N-1; i+=blockDim.x)//X
		{
			if( i == 0 )continue;
			if (index == 1)
			{	
				// en la primera iteración se usa la funcion para H(t=1)
				Hs[0][j*N+i] = Hs[1][j*N+i] + K2 * (Hs[1][j*N+i+1] + Hs[1][j*N+i-1]+Hs[1][(j+1)*N+i]+Hs[1][(j-1)*N+i]-4*Hs[1][j*N+i]);
				if(index==t)
					H_t[j*N+i]=Hs[0][j*N+i];
			}
			else
			{	
				//finalmente para el resto de los t se calcula usando la segunda función
				Hs[0][j*N+i] = 2*Hs[1][j*N+i] - Hs[2][j*N+i] + K1 * (Hs[1][j*N+i+1] + Hs[1][j*N+i-1]+Hs[1][(j+1)*N+i]+Hs[1][(j-1)*N+i]-4*Hs[1][j*N+i]);
				if(index==t)
					H_t[j*N+i]=Hs[0][j*N+i];
			}
		}
	}
	__syncthreads();
	
}

int main(int argc, char * const argv[])
{
	int Grilla;//tamaño grilla
	int Steps;//número de pasos
	int threads_num_X;//numero de hebras
	int threads_num_Y;//numero de hebras
	char* salida_arg=NULL;
	int Steps;//iteración de salida
	int c;

	opterr = 0;
	//GETOPT
	while ((c = getopt (argc, argv, "N:T:X:Y:f:t:")) != -1)
		switch (c)
			{
			case 'N':
				Grilla = atoi(optarg);
				break;
			case 'T':
				Steps = atoi(optarg);
				break;
			case 'X':
				threads_num_X = atoi(optarg);
				break;
			case 'Y':
				threads_num_Y = atoi(optarg);
				break;
			case 'f':
				salida_arg = optarg;
				break;
			case 't':
				Steps = atoi(optarg);
				break;
			case '?':
				if (optopt == 'N' ||optopt == 'T' ||optopt == 'X'||optopt == 'Y' ||optopt == 'f' ||optopt == 't' )
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
	//comprobar inputs
	if(Grilla==0){
		printf("El parametro -N debe estár y ser mayor que 0\n");
		exit(1);
	}
	if(Steps==0){
		printf("El parametro -T debe estár y ser mayor que 0\n");
		exit(1);
	}
	if(threads_num_X==0){
		printf("El parametro -X debe estár y ser mayor que 0\n");
		exit(1);
	}
	if(threads_num_Y==0){
		printf("El parametro -Y debe estár y ser mayor que 0\n");
		exit(1);
	}
	if(Steps==0){
		printf("El parametro -t debe estár y ser mayor que 0\n");
		exit(1);
	}
	if(salida_arg==NULL){
		printf("Debe especificarse un archivo de salida\n");
	}

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