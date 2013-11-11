/********************************************************************************/
/*Cannon's Algoritm								*/
/*Input parameters 	M :: dim of the square matrices A,B,C			*/
/*			N :: dim of the square cartesian process's topology	*/
/*Output		Computation Time in secs for scalability studies	*/
/*Compile		gcc -o openmpcannon openmpcannon.c -Wall -fopenmp	*/
/*Usage			./openmpcannon  <M> <N>					*/
/********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "omp.h"

int main(int argc,char *argv[])
{
	double ***A,***B,***C;							//A,B::source matrices ,C::target matrix
	double t0,t1;								//timer
	int i,j,k;								//counters
	int M,N,m;								

	srand(time(NULL));							//seed random generator with the value of time in seconds since the Epoch
	if ( (argc != 3) || (M=atoi(argv[1])) < 4 || (N=atoi(argv[2])) < 1 || ((int)(m=((int)M/N)) != (float)((float) M/N)) ) 
	{
		printf("%c[%d;%dmPlease provide legal parameters!%c[%dm\n",27,1,31,27,0);
		return -1;
	}

	if ( (A=malloc(N*N*sizeof(double**))) == NULL ) 			//allocation
	{
		printf("Can not allocate mem\n");
		return -1;
	}
	if ( (B=malloc(N*N*sizeof(double**))) == NULL )
	{
		printf("Can not allocate mem\n");
		free(A);
		return -1;
	}
	if ( (C=malloc(N*N*sizeof(double**))) == NULL )
	{
		printf("Can not allocate mem\n");
		free(A);free(B);
		return -1;
	}
	for(i=0;i<N*N;i++)
	{									//allocation
		if ( (A[i]=malloc(m*sizeof(double*))) == NULL ) {printf("Can not allocate mem\n");return -1;}					
		if ( (B[i]=malloc(m*sizeof(double*))) == NULL ) {printf("Can not allocate mem\n");return -1;}	
		if ( (C[i]=malloc(m*sizeof(double*))) == NULL ) {printf("Can not allocate mem\n");return -1;}						
		for(j=0;j<m;j++)						
		{								//allocation
			if ( (A[i][j]=malloc(m*sizeof(double))) == NULL ) {printf("Can not allocate mem\n");return -1;}
			if ( (B[i][j]=malloc(m*sizeof(double))) == NULL ) {printf("Can not allocate mem\n");return -1;}
			if ( (C[i][j]=malloc(m*sizeof(double))) == NULL ) {printf("Can not allocate mem\n");return -1;}
			for(k=0;k<m;k++)					//initialize A,B,C
			{
				A[i][j][k]=(double)rand()/(double)RAND_MAX;	
				B[i][j][k]=(double)rand()/(double)RAND_MAX;	
				C[i][j][k]=0;
			}
		}	
	}

	t0=omp_get_wtime();							//get time
	#pragma omp parallel shared(A,B,C,M,N) num_threads(N*N)			
	{
		int vida,vidb,id,i,j,k,p;
        	id=omp_get_thread_num();					//determinate thread id
		i=(int)id/N;							//determinate initial shift for a
		vida=(id-i<i*N?id-i+N:id-i);					//get a view
		i=(int)id%N;							//determinate initial shift for b
		vidb=( id -i*N < i? id+(N-i)*N : id-i*N);			//get b view

		#pragma omp barrier
		//calculation
		for(p=0; p<N; p++)
		{
			for(i=0; i<m;i++)					//multiplication
				for(j=0;j<m;j++)	
					for(k=0;k<m;k++)
				(C[id])[i][j]+=((A[vida])[i][k])*((B[vidb])[k][j]);

		#pragma omp barrier
			vida=(vida-1  <(vida/N)*N? vida+N-1     :vida-1);	//shifting left
			vidb=(vidb-1*N<(vidb%N)  ? vidb+(N-1)*N :vidb-1*N);	//shifting up
		}
	}
	t1=omp_get_wtime();

	for(i=0;i<N*N;i++)							//deallocation
	{
		for(j=0;j<m;j++) {free(A[i][j]);free(B[i][j]);free(C[i][j]);}
		free(A[i]);free(B[i]);free(C[i]);
	}
	free(A);free(B);free(C);
	printf("A , B , C  	<-- %d x %d matrices\n",M,M);		//output
	printf("Toplogy  	<-- %d x %d threads\n",N,N);
	printf("Computation time  : %.5f secs\n",(t1-t0));
	return 0;
}
