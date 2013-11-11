Matrix-Matrix Multiplication with openMP - Cannon Algorithm 

	+Working Directory includes :

		-readme.txt	 project description
		-openmpcannon.c	 source code , algorithm implementation
		-openMpResults.pdf result's report of sampling@linuxLab@DIT

	+Usage	: M dim of square matrices A,B,C , N dim of square process's topology 

		gcc -o openmpcannon openmpcannon.c -Wall -fopenmp
		./openmpcannon <M> <N>

	+Assumptions

		-M mod N = 0 
		-Output time = Computation Time ( Execution Time without initialization , partitioning,deallocation )
		-The matrices A,B,C are logical M x M matrices allocated as N*N array of M/N x M/N submatrices ( => M x M = N x N x M/N x M/N ) 
