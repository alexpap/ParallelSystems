/********************************************************************************/
/*DNS Algorithm									*/
/*Input	parameter	N dimesion of square matrices A,B,C			*/
/*Output		Computation time					*/
/*Compile		nvcc -o cudadns cudadns.c				*/
/*Usage			./cudadns <N>						*/
/********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

__global__ void dns(float *A,float *B,float *C,int N) 
{
	extern __shared__ float sum[];

	sum[threadIdx.x]=A[blockIdx.x*N+threadIdx.x]*B[threadIdx.x*N+blockIdx.y]; // multiplication 

	__syncthreads();

	for(unsigned int s=N/2; s>0; s>>=1) 					// reduction
	{
        	if (threadIdx.x < s) 
        	{
        	    sum[threadIdx.x] += sum[threadIdx.x + s];
        	}
        	__syncthreads();
    	}
	if ( !threadIdx.x ) C[blockIdx.x*N+blockIdx.y]=sum[threadIdx.x];	// result
}

int main(int argc, char* argv[])
{

	float *A, *B, *C;							// A,B,C matrices@host
	float *devA,*devB,*devC;						// A,B,C matrices @device 
	int i,j,N=0;								//i,j: counters , N:dims of A,B,C
	float elapsedTime;							//timers
	cudaEvent_t start, stop;
	
	if ( argc != 2 || (N=atoi(argv[1])) < 1 )				//check cmd line args
	{
		printf("Please provide legal args\n");
		return -1;
	}
	(void)srand(time(NULL));				//seed random generator with the value of time in seconds since the Epoch

	A = (float*)malloc( N*N*sizeof(float) );				//Allocate host mem
	B = (float*)malloc( N*N*sizeof(float) );
	C = (float*)malloc( N*N*sizeof(float) );

	for( i=0; i<N ;i++)							//generate values for matrices A,B
		for( j=0; j<N; j++)
		{
			A[i*N+j]=(float)rand()/(float)RAND_MAX;	
			B[i*N+j]=(float)rand()/(float)RAND_MAX;
			C[i*N+j]=0.0;						//initialize matrix C
		}
	
	cudaMalloc((void**)&devA,N*N*sizeof(float));					//Allocate  dev mem 
	cudaMalloc((void**)&devB,N*N*sizeof(float));	
	cudaMalloc((void**)&devC,N*N*sizeof(float));	
	cudaMemcpy((void*)devA,A,N*N*sizeof(float),cudaMemcpyHostToDevice);		// Copy A,B,C to device 
	cudaMemcpy((void*)devB,B,N*N*sizeof(float),cudaMemcpyHostToDevice);		
	cudaMemcpy((void*)devC,C,N*N*sizeof(float),cudaMemcpyHostToDevice);  	

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);						// start event	
	
	dns<<<dim3(N,N),N,N*sizeof(float)>>>(devA,devB,devC,N);		//run dns kernel with N*N blocks of N threads and a table sizeof N

	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
	{ 
		printf("Error: %s\n", cudaGetErrorString(err));
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		cudaFree(devA);cudaFree(devB);cudaFree(devC);				//deallocate device mem
		free(A);free(B);free(C);						//deallocate host mem
		return -1;
	}

	cudaEventRecord(stop, 0);						// Stop event
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop); 			// get elapsed time
	cudaMemcpy(C, devC, N*N*sizeof(float), cudaMemcpyDeviceToHost);		// get matrix C back from device
	cudaFree(devA);cudaFree(devB);cudaFree(devC);				//deallocate device mem
	free(A);free(B);free(C);						//deallocate host mem
	
	printf("A, B, C 	<-- %d x %d matrices\n",N,N);					//output				
	printf("Computation time  : %f\n", elapsedTime);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	return 0;
}
