#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cuda_runtime_api.h"
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>

#define result_tag 11

const int INF = 10000000;
const int V = 10010;
void input(char *inFileName);
void output(char *outFileName);

void block_FW(int B,int numdevs);
int ceil(int a, int b);
void cal(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height);
int init_device();

int n, m;	// Number of vertices, edges
static int Dist[V][V];
int *gpudist0;
int *gpudist1;
int *hostdist;
int *hostdist2;


//MPI
int size, rank;
MPI_Status status;
MPI_Request status_req[2];

//time
/*double comtime=0;
double cmitime=0;
double memcpy_t=0;*/



__global__ void calculat(int B, int n, int Round, int block_start_x, int block_start_y, int block_width, int block_height, int *gpudist){
	//for (int b_i =  block_start_x; b_i < block_end_x; ++b_i) {
	//	for (int b_j = block_start_y; b_j < block_end_y; ++b_j) {
			// To calculate B*B elements in the block (b_i, b_j)
			// For each block, it need to compute B times
	//		for (int k = Round * B; k < (Round +1) * B && k < n; ++k) {
				// To calculate original index of elements in the block (b_i, b_j)
				// For instance, original index of (0,0) in block (1,2) is (2,5) for V=6,B=2
	int b_i = blockIdx.x+block_start_x;
	int b_j = blockIdx.y+block_start_y;
	int distij,distik,distkj;// register value
	int block_internal_start_x = b_i * B;
	int block_internal_start_y = b_j * B; 
	int j = block_internal_start_y + threadIdx.x; // set column
	if(j > n-1) j=n-1;
	for (int k = Round * B; k < (Round +1) * B && k < n; ++k) {
		distkj = gpudist[k*n+j];	
		for (int blk_idx = 0;blk_idx<B;blk_idx++){
				//int block_internal_start_x = b_i * B;
	//			int block_internal_end_x   = (b_i +1) * B;
				//int block_internal_start_y = b_j * B; 
	//			int block_internal_end_y   = (b_j +1) * B;

	//			if (block_internal_end_x > n)	block_internal_end_x = n;
	//			if (block_internal_end_y > n)	block_internal_end_y = n;

	//			for (int i = block_internal_start_x; i < block_internal_end_x; ++i) {
	//				for (int j = block_internal_start_y; j < block_internal_end_y; ++j) {
				int i = block_internal_start_x + blk_idx;
				//int j = block_internal_start_y + threadIdx.x; // set column
				if(i > n-1) i=n-1;
				
				
				distij = gpudist[i*n+j];
				distik = gpudist[i*n+k];
					if (distik + distkj < distij)
						gpudist[i*n+j] = distik + distkj;
				__syncthreads();
		}
	}
	
			//}
		//}
	//}
}
static __global__ void calculat32(int B, int n, int Round, int block_start_x, int block_start_y, int block_width, int block_height, int *dist){

		int b_i = blockIdx.x+block_start_x;
        int b_j = blockIdx.y+block_start_y;
        int distij,distik,distkj;// register value
        int block_internal_start_x = b_i * B;
        int block_internal_start_y = b_j * B;
        int j = block_internal_start_y + threadIdx.y;
        int i = block_internal_start_x + threadIdx.x;
        if(i > n-1) i=n-1;
        if(j > n-1) j=n-1;
            distij = dist[i*n+j];
                for (int k = Round * B; k < (Round +1) * B && k < n; ++k) {
                        distik = dist[i*n+k];
                        distkj = dist[k*n+j];
                        if (distik + distkj < distij){
                                distij = distik + distkj;
                                dist[i*n+j]=distij;
                        }
                        __syncthreads();
                }

}


double timer(void)
{	struct timeval tv;
	struct timezone tz;
  
  	double t;

  	gettimeofday(&tv, &tz);

  	t = (double)tv.tv_sec*1000;
  	t += ((double)tv.tv_usec)/1000.0;

  	return t/1000;
}

int main(int argc, char* argv[])
{
	//MPI init
	
	MPI_Init(&argc,&argv);
	MPI_Comm_size (MPI_COMM_WORLD, &size); 
	MPI_Comm_rank (MPI_COMM_WORLD, &rank); 
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*
	//time
	double t_st,t_end;
	struct timeval tv;
	struct timezone tz;
	clock_t t0,t1;
	t_st = timer();
	t0 = clock();
	double memcpy_start_t = 0, memcpy_end_t = 0;
	double IO_start_t = 0, IO_end_t = 0;
	double IO_t=0;
	//time
	*/
	
	cudaError_t err;
	int numdevs; // get devices number
	//printf("go~%d\n",argc);
	//IO_start_t = timer();
	input(argv[1]);
	/*IO_end_t = timer();
	IO_t = IO_t + IO_end_t -IO_start_t;
	*/
	int B = atoi(argv[3]);
	//printf("%d",B);
	
	
	//init gpu && openMP threads
	numdevs = init_device();
	//omp_set_num_threads(numdevs);
	
	//gpudist = (int *)malloc(sizeof(int)*2);
	//allocate GPU memory
	//err = cudaMalloc((void**)&gpudist0,sizeof(int)*n*n);
	//printf("cudaMallocgpudist0 %s\n",err);
	//if(size>1 && rank ==1){
	cudaSetDevice(rank);
	//memcpy_start_t = timer();
	cudaMalloc((void**)&gpudist0,sizeof(int)*n*n);
	//printf("cudaMallocgpudist1 %s\n",err);
	//}
	//copy DIST to Device(GPU)
	//cudaSetDevice(0);
	cudaMemcpy(gpudist0,hostdist,sizeof(int)*n*n,cudaMemcpyHostToDevice);
	/*memcpy_end_t = timer();
	memcpy_t = memcpy_t +memcpy_end_t - memcpy_start_t;*/
	
	
	block_FW(B,numdevs);
	
	
	//memcpy_start_t = timer();
	cudaMemcpy(hostdist,gpudist0,sizeof(int)*n*n,cudaMemcpyDeviceToHost);
	/*memcpy_end_t = timer();
	memcpy_t = memcpy_t +memcpy_end_t - memcpy_start_t;*/
	
	cudaFree(gpudist0);
	//cudaFree(gpudist1);
	//IO_start_t = timer();
	if(rank==0)
		output(argv[2]);
	//IO_end_t = timer();
	//IO_t = IO_t + IO_end_t -IO_start_t;
	
	//free(gpudist);
	free(hostdist);
	free(hostdist2);
	
	/*t_end = timer();
	t1 = clock();
	//printf("#b:%d Total Cuda[%d] time = %lf sec (%lf)\n" , B,rank, (t_end - t_st),(t1-t0)/(double)(CLOCKS_PER_SEC));
	//printf("Cuda[%d] computation time = %f   IO = %f  communication = %f  memcpy = %f\n",rank, comtime, IO_t, cmitime, memcpy_t);
	*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

void input(char *inFileName)
{
	FILE *infile = fopen(inFileName, "r");
	fscanf(infile, "%d %d", &n, &m);
	//allocate nxn matrix
	hostdist = (int *)malloc(sizeof(int)*n*n);
	hostdist2 = (int *)malloc(sizeof(int)*n*n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j){	
					hostdist2[i*n+j] = 0;
					hostdist[i*n+j]=0;
			}
			else{		
					hostdist2[i*n+j] = INF;
					hostdist[i*n+j]=INF;
			}
		}
	}

	while (--m >= 0) {
		int a, b, v;
		fscanf(infile, "%d %d %d", &a, &b, &v);
		--a, --b;
		hostdist2[a*n+b] = v;
		hostdist[a*n+b] = v;
	}
}

void output(char *outFileName)
{
	FILE *outfile = fopen(outFileName, "w");
	if(rank == 0){
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (hostdist[i*n+j] >= INF)	fprintf(outfile, "INF ");
				else					fprintf(outfile, "%d ", hostdist[i*n+j]);
			}
			fprintf(outfile, "\n");
		}		
	}
}

int ceil(int a, int b)
{
	return (a + b -1)/b;
}

void block_FW(int B, int numdevs)
{
	//omp_set_num_threads(B);
	double str;



	int round = ceil(n, B);//?–n/B?„ä???	
	//printf("round %d\n",round);
	for (int r = 0; r < round; ++r) {
		/* Phase 1*/
		//str=timer();
		cal(B,	r,	r,	r,	1,	1);
		//printf("phase1\n");
		/* Phase 2*/
		cal(B, r,     r,     0,             r,             1);
		cal(B, r,     r,  r +1,  round - r -1,             1);
		cal(B, r,     0,     r,             1,             r);
		cal(B, r,  r +1,     r,             1,  round - r -1);
		//printf("phase2\n");
		cudaDeviceSynchronize();
		//comtime +=(timer()-str);
		//if(numdevs >1 ){
		//	cudaMemcpy(gpudist1,gpudist0,sizeof(int)*n*n,cudaMemcpyDeviceToDevice);
		//	printf("D2D\n");
		//}
		//cudaMemcpy(hostdist,gpudist0,sizeof(int)*n*n,cudaMemcpyDeviceToHost);
		//output("iniout");
		/* Phase 3*/
		//#pragma omp parallel
		//{
		//str=timer();
			if(size==1){
			cal(B, r,     0,     0,            r,             r);
			cal(B, r,     0,  r +1,  round -r -1,             r);
			cal(B, r,  r +1,     0,            r,  round - r -1);
			cal(B, r,  r +1,  r +1,  round -r -1,  round - r -1);
			}
			else{
				cudaSetDevice(rank);
				cal(B, r,  (r +1)*rank,     0,            r,  (round - r -1)*rank-r*(rank-1));
				cal(B, r,  (r +1)*rank,  r +1,  round -r -1,  (round - r -1)*rank-r*(rank-1));
				/*
				(round - r -1)*rank-r*(rank-1)
				if rank ==0 then (round - r -1)*0-r*(0-1) => r
				if rank ==1 then (round - r -1)*1-r*(1-1) => round - r -1
				*/
				cudaDeviceSynchronize();
			}
			//comtime +=(timer()-str);
			
		//}
	//str = timer();
		if(size >1 ){
			cudaSetDevice(rank);
			//cudaMemcpy(hostdist2,gpudist1,sizeof(int)*n*n,cudaMemcpyDeviceToHost);
			//cudaSetDevice(0);
			cudaMemcpy(hostdist,gpudist0,sizeof(int)*n*n,cudaMemcpyDeviceToHost);
			int j=0;
			int i=0;
			memset(hostdist2,0,sizeof(int)*n*n);
			//MPI_Barrier(MPI_COMM_WORLD);
			//printf("will send");
			/*if(rank==1){
				MPI_Irecv(hostdist2,n*n,MPI_INT,1,result_tag,MPI_COMM_WORLD,&status_req[0]);
				MPI_Isend(hostdist,n*n,MPI_INT,0,result_tag,MPI_COMM_WORLD,&status_req[1]);
				printf("will recive");
			}
			if(rank==0){
				MPI_Irecv(hostdist2,n*n,MPI_INT,0,result_tag,MPI_COMM_WORLD,&status_req[0]);
				MPI_Isend(hostdist,n*n,MPI_INT,1,result_tag,MPI_COMM_WORLD,&status_req[1]);
				printf("will recive");
			}
			MPI_Waitall(1,status_req,&status);*/
			if(rank==1){
				MPI_Send(hostdist,n*n,MPI_INT,size-1-rank,result_tag,MPI_COMM_WORLD);
				MPI_Recv(hostdist2,n*n,MPI_INT,MPI_ANY_SOURCE,result_tag,MPI_COMM_WORLD,&status);
			}
			else{
				MPI_Recv(hostdist2,n*n,MPI_INT,MPI_ANY_SOURCE,result_tag,MPI_COMM_WORLD,&status);
				MPI_Send(hostdist,n*n,MPI_INT,size-1-rank,result_tag,MPI_COMM_WORLD);
			}
			
			//cmitime = cmitime +timer() - str;
			
			//			str = timer();
			#pragma omp parallel private(j) 
			{
				#pragma omp for schedule(dynamic)
				for(i=0;i<n;i++)
					for(j=0;j<n;j++){
						if(hostdist[i*n+j]>hostdist2[i*n+j]){
							hostdist[i*n+j]=hostdist2[i*n+j];
							//printf("1dist[%d][%d]= %d  ",i,j,hostdist[i*n+j]);
							//printf("2dist[%d][%d]= %d  ",i,j,hostdist2[i*n+j]);
						}
						//if(hostdist[i*n+j]==hostdist2[i*n+j])
							//printf("fuck\n");
					}
			}
			
			//comtime+=(timer()-str);
			//str = timer();
			cudaMemcpy(gpudist0,hostdist,sizeof(int)*n*n,cudaMemcpyHostToDevice);
			//memcpy_t = memcpy_t +timer() - str;
		}
	}
}

void cal(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height)
{
	//int block_end_x = block_start_x + block_height;
	//int block_end_y = block_start_y + block_width;
	dim3 guid_size = dim3(block_height, block_width);
	//dim3 block_size = dim3(B,B);
	dim3 block_size;
	if(B<=32)
		block_size = dim3(B,B);
	else
		block_size = dim3(B,1);
		
		if(B<=32)
				calculat32<<<guid_size,block_size>>>(B,n,Round,block_start_x,block_start_y,block_width,block_height,gpudist0);
		else	calculat<<<guid_size,block_size>>>(B,n,Round,block_start_x,block_start_y,block_width,block_height,gpudist0);

	


}

int init_device(){
	cudaError_t err;
	int numdevs;
	//printf("go");
	cudaGetDeviceCount(&numdevs);
	if(numdevs > 1){
	err = cudaSetDevice(rank);
	//printf("suda set =%s\n",err);
	}
	else
		cudaSetDevice(0);
	//printf("numdevs=%d\n",numdevs);
	
	return numdevs;
}

