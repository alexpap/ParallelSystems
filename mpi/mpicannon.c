/********************************************************************************/
/*Cannon's Algoritm								*/
/*Input parameters 	M :: dim of the square matrices A,B,C			*/
/*			N :: dim of the square cartesian process's topology	*/
/*Output		Computation Time in secs for scalability studies	*/
/*Mpd conf		mpdboot -c -v -r rsh -f mpd.hosts 			*/
/*				--ncpus=<#cpus> -n <#machines> 			*/
/*Compile		mpicc -o mpicannon mpicannon.c -Wall -lm		*/
/*Usage			mpiexec -machinefile machines -np <N*N> ./mpicannon <M>	*/
/********************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<time.h>

int main(int argc, char *argv[])
{
	int np,rank;								//np : number of processes ,rank : rank of current process
	int M,N;								//M  : dims of A,B,C square matrices ,N : dims of square Cartesian Topology
	int m;									//m  : dims of square submatrices a,b,c
	float *A=NULL,*B=NULL,*C=NULL;						//A,B,C: serialized square matrices
	float *a=NULL,*b=NULL,*c=NULL;						//a,b,c: serialized square submatrices
	double elapsed=0,elocal;						//timers
	int p,i,j,k;								//p,i,j,k : counters
	MPI_Comm cartComm;							//cartComm : Topology Communicator
	MPI_Status status;							//status :Recv's Status
	int dims[2],periods[2],coords[2];					//Characteristics of square cartesian Topology
	int shiftsource,shiftdest;		
	int up,down,left,right;							

	MPI_Init(&argc,&argv);							//init MPI envi
	MPI_Comm_size(MPI_COMM_WORLD,&np);					//get envi info
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
										//determinate the dims of matrices A,B,C & the dims of local submatrices
	if (argc != 2 || (M=atoi(argv[1])) < 4 || (N=sqrt(np)) < 1 || (int)(m=((int)M/N)) != (float)((float) M/N) ) 				
	{
		if (!rank ) printf("%c[%d;%dmPlease provide legal parameters!%c[%dm\n",27,1,31,27,0);
		MPI_Abort(MPI_COMM_WORLD,MPI_ERR_ARG);
		MPI_Finalize();
		return -1;
	}

	a=malloc(m*m*sizeof(float));						//allocate block for every process
	b=malloc(m*m*sizeof(float));		
	c=malloc(m*m*sizeof(float));

	if ( !rank )								//root process
	{
		(void)srand(time(NULL));					//seed number generator
		A=malloc(M*M*sizeof(float));					//allocate mem for the 'total' matrices
		B=malloc(M*M*sizeof(float));		
		C=malloc(M*M*sizeof(float));
		for( i=0; i<M*M; i++)						//initialize matirices
		{
			A[i]=(float)rand()/(float)RAND_MAX;	
			B[i]=(float)rand()/(float)RAND_MAX;	
			C[i]=0;
		}
	}

	MPI_Scatter(A,m*m,MPI_FLOAT,a,m*m,MPI_FLOAT,0,MPI_COMM_WORLD);	//spread/partition the values of A,B
	MPI_Scatter(B,m*m,MPI_FLOAT,b,m*m,MPI_FLOAT,0,MPI_COMM_WORLD);

	dims[0]=dims[1]=N;							//initialize square cart topology characteristics
	periods[0] = periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&cartComm);		//create topology
	MPI_Comm_rank(cartComm,&rank);						//get topology info 
	MPI_Cart_coords(cartComm,rank,2,coords);
	MPI_Cart_shift(cartComm,1,-1,&right,&left);				
	MPI_Cart_shift(cartComm,0,-1,&down,&up);
	MPI_Cart_shift(cartComm,1,-coords[0],&shiftsource,&shiftdest);		//initial stage
	MPI_Sendrecv_replace(a,m*m,MPI_FLOAT,shiftdest,1,shiftsource,1,cartComm,&status);
	MPI_Cart_shift(cartComm,0,-coords[1],&shiftsource,&shiftdest);
	MPI_Sendrecv_replace(b,m*m,MPI_FLOAT,shiftdest,1,shiftsource,1,cartComm,&status);

	MPI_Barrier(MPI_COMM_WORLD);						//start timer
	elocal=MPI_Wtime();
	for (p=0; p<dims[0]; p++)						//main calulation
	{
		for (i=0; i<m; i++)
			for (j=0; j<m; j++)
				for (k=0; k<m; k++)
					c[i*m+j] += a[i*m+k]*b[k*m+j];

		MPI_Sendrecv_replace(a,m*m, MPI_FLOAT,left,1,right,1,cartComm,&status);
		MPI_Sendrecv_replace(b,m*m,MPI_FLOAT,up,1,down,1,cartComm,&status);
	}
	MPI_Barrier(MPI_COMM_WORLD);						//stop timer
	elocal=MPI_Wtime()-elocal;
	//MPI_Cart_shift(cartComm,1,+coords[0],&shiftsource,&shiftdest);		//'undo' initial stage for A,B
	//MPI_Sendrecv_replace(a,m*m,MPI_FLOAT,shiftdest,1,shiftsource,1,cartComm,&status);
	//MPI_Cart_shift(cartComm,0,+coords[1],&shiftsource,&shiftdest);
	//MPI_Sendrecv_replace(b,m*m,MPI_FLOAT,shiftdest,1,shiftsource,1,cartComm,&status);
	MPI_Comm_free(&cartComm);						//destroy topology
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);					//gather results from every process
	//MPI_Gather(c,m*m,MPI_FLOAT,C,m*m,MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Reduce(&elocal,&elapsed,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	if ( !rank )
	{
		printf("A , B , C  	<-- %d x %d matrices\n",M,M);		//output
		printf("Toplogy  	<-- %d x %d threads\n",N,N);
		printf("Computation time  : %.5f sec\n",elapsed);		//report execution time		
		free(A);free(B);free(C);					//deallocate mem of A,B,B
	}
	free(a);free(b);free(c);						//deallocate mem of local a,b,c
	MPI_Finalize();
	return 0;
}
