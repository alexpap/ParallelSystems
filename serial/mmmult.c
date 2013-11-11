/********************************************************************************/
/*The conventional serial algorithm for multiplication of two n x n matrices.	*/
/*Input	parameter	N :: dim of the square matrices A,B,C			*/
/*Output		Computation Time in secs for scalability studies	*/
/*Compile		gcc -o mmmult mmmult.c -Wall				*/
/*Usage			./mmmult <N>						*/
/********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc,char *argv[])
{
	clock_t c0,c1;								//timers
	int i,j,k,n=0;								//i,j,k :: counters ,n::input parameter
	double *A=NULL,*B=NULL,*C=NULL;						//A,B::source matrices ,C::target matrix

	(void)srand(time(NULL));						//seed random generator with the value of time in seconds since the Epoch
	if ( argc != 2 || (n=atoi(argv[1])) < 1 ) 				//Check/initialize cmd line arguments
	{
		printf("%c[%d;%dmPlease provide legal parameters!%c[%dm\n",27,1,31,27,0);
		return 0;
	}
	if ( (A=malloc(n*n*sizeof(double))) == NULL ) 				//Allocate memory for the matrices
	{
		printf("Can not allocate mem\n");
		return -1;
	}
	if ( (B=malloc(n*n*sizeof(double))) == NULL) 
	{
		printf("Can not allocate mem\n");
		free(A);
		return -1;
	}
	if ( (C=malloc(n*n*sizeof(double))) == NULL) 
	{
		printf("Can not allocate mem\n");
		free(A);free(B);
		return -1;
	}
	for( i=0; i<n ;i++)							//generate values for matrices A,B
		for( j=0; j<n; j++)
		{
			A[i*n+j]=(double)rand()/(double)RAND_MAX;	
			B[i*n+j]=(double)rand()/(double)RAND_MAX;
			C[i*n+j]=0.0;						//initialize matrix C
		}
	if ( (c0=clock()) == -1 ) {printf("clock failure\n");return -1;}	//get time

	for( i=0; i<n ;i++)							//Calculate C=A*B
		for( j=0; j<n; j++)
			for( k=0; k<n; k++)
				C[i*n+j]+=A[i*n+k]*B[k*n+j];

	if ( (c1=clock()) == -1 ) {printf("clock failure\n");return -1;}	//get time
	printf("A , B , C 	<-- %d x %d \n",n,n);				//output	
	printf("Computation time  : %.5f secs\n",(double)(c1-c0)/CLOCKS_PER_SEC);				
	free(A);free(B);free(C);						//deallocate mem
	return 0;
}
