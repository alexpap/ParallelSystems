Matrix-Matrix Multiplication with CUDA - DNS Algorithm 

	+Working Directory includes :

		-readme.txt	 project description
		-cudadns.cu	 source code , algorithm implementation
		-cudaResults.pdf result's report of sampling@GTX-480

	+Usage	: N dim of square matrices A,B,C

		nvcc -o cudadns cudadns.cu 
		./cudadns <N>

	+Assumptions

		-N <= 1024 , has the form 2^i where i= 2,3,4... ( N is bounded from the card's capabilities(#threads/block) , 1024 is for nvidia GTX-480)
		-Output time in milisecs = Computation Time ( Execution Time without initialization , deallocation )
		-The matrices A,B,C in CPU main memory are 2D square tables
		-GPU kernel creates a 3D logical process's topology (as DNS required) with N x N grid of blocks ,where every block has N threads as z-dimension
		-In GPU main memory we use the thread's indexing((x,y),z) to point to the  apropiate(as DNS required) element of 2D matrices A,B,C without further allocation
		-Each Block of processe calculate an element of C with reduction to spread the workload among the threads (4reduction T(n)=log(n)) 

