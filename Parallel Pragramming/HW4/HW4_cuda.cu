#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "cuda_runtime_api.h"
#include <time.h>
#include <sys/time.h>

const int INF = 10000000;
const int V = 10010;
void input(char *inFileName);
void output(char *outFileName);

void block_FW(int B);
int ceil(int a, int b);
void cal(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height);
void init_device();

int n, m;	// Number of vertices, edges
static int Dist[V][V];
int *gpudist;
int *hostdist;

//time
//double comtime=0;
//double cmitime=0;

//cudaError_t cudaEventElapsedTime(float* ms,cudaEvent_t start,cudaEvent_t end)

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
	for (int k = Round * B; k < (Round +1) * B && k < n; ++k) {
		for (int blk_idx = 0;blk_idx<B;blk_idx++){
				int block_internal_start_x = b_i * B;
	//			int block_internal_end_x   = (b_i +1) * B;
				int block_internal_start_y = b_j * B; 
	//			int block_internal_end_y   = (b_j +1) * B;

	//			if (block_internal_end_x > n)	block_internal_end_x = n;
	//			if (block_internal_end_y > n)	block_internal_end_y = n;

	//			for (int i = block_internal_start_x; i < block_internal_end_x; ++i) {
	//				for (int j = block_internal_start_y; j < block_internal_end_y; ++j) {
				int i = block_internal_start_x + blk_idx;
				int j = block_internal_start_y + threadIdx.x; // set column
				if(i > n-1) i=n-1;
				if(j > n-1) j=n-1;
				
				distij = gpudist[i*n+j];
				distik = gpudist[i*n+k];
				distkj = gpudist[k*n+j];	
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
{	cudaError_t err;
	
	//time
	/*double t_st,t_end;
	struct timeval tv;
	struct timezone tz;
	clock_t t0,t1;
	t_st = timer();
	t0 = clock();
	double memcpy_start_t = 0, memcpy_end_t = 0;
	double memcpy_t=0;
	double IO_start_t = 0, IO_end_t = 0;
	double IO_t=0;
	
	//IO
	IO_start_t = timer();*/
	input(argv[1]);
	//IO_end_t = timer();
	//IO_t = IO_t + IO_end_t -IO_start_t;
	
	int B = atoi(argv[3]);
	//printf("%d",B);
	
	
	//init gpu
	init_device();
	//printf("init_device %d \n",B);
	
	
	//allocate GPU memory
	//memcpy_start_t = timer();
	cudaMalloc((void**)&gpudist,sizeof(int)*n*n);
	//printf("cudaMalloc %s \n",err);
	//copy DIST to Device(GPU)
	cudaMemcpy(gpudist,hostdist,sizeof(int)*n*n,cudaMemcpyHostToDevice);
	//memcpy_end_t = timer();
	//memcpy_t = memcpy_t +memcpy_end_t - memcpy_start_t;
	//printf("cudaMemcpy %d \n",B);
	
	
	block_FW(B);
	
	
	//memcpy_start_t = timer();
	cudaMemcpy(hostdist,gpudist,sizeof(int)*n*n,cudaMemcpyDeviceToHost);
	//memcpy_end_t = timer();
	//memcpy_t = memcpy_t +memcpy_end_t - memcpy_start_t;
	
	
	
	cudaFree(gpudist);
	
	//IO
	//IO_start_t = timer();
	output(argv[2]);
	//IO_end_t = timer();
	//IO_t = IO_t + IO_end_t -IO_start_t;
	
	
	free(hostdist);
	
	//gettimeofday(&tv, &tz);
	//t_end = timer();
	//t1 = clock();
	//printf("#b:%d Total Cuda time = %lf sec (%lf)\n" , B, (t_end - t_st),(t1-t0)/(double)(CLOCKS_PER_SEC));
	//printf("computation time = %f  IO = %f  communication = %f  memcpy = %f\n", comtime/1000, IO_t, cmitime, memcpy_t);
	
	return 0;
}

void input(char *inFileName)
{
	FILE *infile = fopen(inFileName, "r");
	fscanf(infile, "%d %d", &n, &m);
	//allocate nxn matrix
	hostdist = (int *)malloc(sizeof(int)*n*n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j){	
					hostdist[i*n+j]=0;}
			else{		
					hostdist[i*n+j]=INF;}
		}
	}

	while (--m >= 0) {
		int a, b, v;
		fscanf(infile, "%d %d %d", &a, &b, &v);
		--a, --b;
		
		hostdist[a*n+b] = v;
	}
}

void output(char *outFileName)
{
	FILE *outfile = fopen(outFileName, "w");
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (hostdist[i*n+j] >= INF)	fprintf(outfile, "INF ");
			else					fprintf(outfile, "%d ", hostdist[i*n+j]);
		}
		fprintf(outfile, "\n");
	}		
}

int ceil(int a, int b)
{
	return (a + b -1)/b;
}

void block_FW(int B)
{	
	//cudaEvent_t start, stop;
	//cudaEventCreate(&start);
	//cudaEventCreate(&stop);
	//float time;
	
	
	int round = ceil(n, B);//?–n/B?„ä???	
	for (int r = 0; r < round; ++r) {
		/* Phase 1*/
		//cudaEventRecord(start);
		cal(B,	r,	r,	r,	1,	1);
		//cudaEventRecord(stop);
		//cudaEventSynchronize(stop);
		//cudaEventElapsedTime(&time, start, stop);
		//comtime +=time;
		
		/* Phase 2*/
		//cudaEventRecord(start);
		cal(B, r,     r,     0,             r,             1);
		cal(B, r,     r,  r +1,  round - r -1,             1);
		cal(B, r,     0,     r,             1,             r);
		cal(B, r,  r +1,     r,             1,  round - r -1);
		//cudaDeviceSynchronize();
		//cudaEventRecord(stop);
		//cudaEventSynchronize(stop);
		//cudaEventElapsedTime(&time, start, stop);
		//comtime +=time;
		cudaDeviceSynchronize();
		/* Phase 3*/
		//cudaEventRecord(start);
		cal(B, r,     0,     0,            r,             r);
		cal(B, r,     0,  r +1,  round -r -1,             r);
		cal(B, r,  r +1,     0,            r,  round - r -1);
		cal(B, r,  r +1,  r +1,  round -r -1,  round - r -1);
		cudaDeviceSynchronize();
		//cudaDeviceSynchronize();
		//cudaDeviceSynchronize();
		//cudaEventRecord(stop);
		//cudaEventSynchronize(stop);
		//cudaEventElapsedTime(&time, start, stop);
		//comtime +=time;
	}
}

void cal(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height)
{
	//int block_end_x = block_start_x + block_height;
	//int block_end_y = block_start_y + block_width;
	dim3 guid_size = dim3(block_height, block_width);
	dim3 block_size;
	if(B<=32)
		block_size = dim3(B,B);
	else
		block_size = dim3(B,1);
		//printf("a=%d\n",k);
	if(B<=32)	
		calculat32<<<guid_size,block_size>>>(B,n,Round,block_start_x,block_start_y,block_width,block_height,gpudist);
	else
		calculat<<<guid_size,block_size>>>(B,n,Round,block_start_x,block_start_y,block_width,block_height,gpudist);
}

void init_device(){
	int numdevs;
	//printf("go");
	cudaGetDeviceCount(&numdevs);
	if(numdevs > 0){
	cudaSetDevice(0);
	}
	//printf("numdevs=%d\n",numdevs);
}


