Matrix-Matrix Multiplication with MPI - Cannon's Algorithm 

	+Working Directory includes :

		-readme.txt	project description
		-machines	file mapping procs to machines
		-mpd.hosts	file specifies the file of machines to start the rest of the mpds on
		-mpicannon.c	source code , algorithm implementation
		-MpiResults.pdf result's report of sampling@linux'sLab@DIT

	+Usage	: N - dim of square process's topology, M - dim of square matrices A,B,C

		mpdboot -f mpd.hosts -r rsh -c -v -n <#machines>
		mpiexec -np <N*N> ./mpicannon M

	+Assumptions

		-M mod N = 0
		-Output time = Calculation Time ( Execution Time without initialization , partitioning , deallocation )
		-The matrices A,B,C = { {block(0,0)},{block(0,1)},...,{block(M/N,M/N)} } , where {block(i,j)} :: M/N x M/N 2D-Matrix 
		-Τhe shaded code in source code ,can be used for concentration of the matrices A,B,C which in our approach is unnecessary.
