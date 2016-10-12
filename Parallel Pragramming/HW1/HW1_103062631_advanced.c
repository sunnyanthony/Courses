#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define ROOT 0


int compare(const void * a, const void * b){
	int *a1=(int *)a;
	int *a2=(int *)b;
	if(*a1>*a2) return 1;
	if(*a1==*a2) return 0;
	if(*a1<*a2) return -1;
}

int Merge(int *local_new_arr,int *local_arr,int leftCount,int *local_comp,int rightCount) {
	int i,j,k;
	i = 0; j = 0; k =0;
	int sorted=0;
//	printf("comp=%d",local_comp[0]);
//	printf("**%d*%d*%d*",i,j,k);
	while(i<leftCount && j< rightCount) {
		if((local_arr[i]  < local_comp[j]) || (local_arr[i]==local_comp[j])) local_new_arr[k++] = local_arr[i++];
		else {
			local_new_arr[k++] = local_comp[j++];
			sorted=1;
		}	
	}
	//printf("--%d-%d--%d----",i,local_new_arr[0],local_new_arr[1]);
	while(i < leftCount) local_new_arr[k++] = local_arr[i++];
	while(j < rightCount) local_new_arr[k++] = local_comp[j++];
	//printf("XX %d XX %d XX",local_new_arr[0],local_new_arr[1]);
	if(sorted){
		for(i=0;i<leftCount;i++)
			local_arr[i]=local_new_arr[i];
	}
//	for(j=0;j<leftCount*2;j++)
//		printf("h(%d)ihi=%d",j,local_new_arr[i]);
	return sorted;
}

int main (int argc, char *argv[]) { //argv is file to read
	int rank, size;	//rank is taskID of communicator
	printf("start");	//size is the number of processes in the group
	MPI_Status status;
	MPI_Init(&argc, &argv);//initializes MPI resource
	MPI_Comm_size(MPI_COMM_WORLD, &size);//get size by address
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//get taskID
	double cput=MPI_Wtime();
	if (argc < 3) {//argc must 3 , so less 3 that is error case
		if (rank == ROOT) {	//defined taskID==ROOT to be master
			fprintf(stderr, "Insufficient args\n");
			fprintf(stderr, "Usage: %s N input_file", argv[0]);
		}
		MPI_Barrier(MPI_COMM_WORLD);//all processes be Blocked until synchronized
		MPI_Finalize();
		return 0;
	}
	const int N = atoi(argv[1]);//N is number of file size
	const char *inName = argv[2];//isName is file name

	int *root_arr; // for root process (which rank == 0) only
	
	int remain = N % size; //
	int odd = rank & 1;
	int even = !odd;
	int num_per_node = N / size;
	MPI_File f;
	MPI_Offset offset;
	if(remain>0)num_per_node++;
	int *local_arr=malloc((num_per_node)*sizeof(int));
   	MPI_File_open(MPI_COMM_WORLD,inName,MPI_MODE_RDONLY,MPI_INFO_NULL,&f);

	offset=rank*num_per_node*sizeof(int);
		MPI_File_read_at(f,offset,local_arr,num_per_node,MPI_INT,&status);//root read file

		int d;
		int x=num_per_node-size+remain;//if x>=0 i=x and this only size-1
		int y=(-x)/num_per_node;//if x<0 size-1-(y+1) ~size-1
		int z=(-x)%num_per_node;//if x<0 size-2-y i=num_per_node-z
		if(remain != 0 && x>=0 && rank==size-1){
			
			for(d=x;d<num_per_node;d++)
				local_arr[d]=2147483647;
		}else if(remain !=0 && x<0){
			if(rank>=size-1-y)
				for(d=0;d<num_per_node;d++)
					local_arr[d]=2147483647;
			if(z>0 && rank==size-2-y)
				for(d=num_per_node-z;d<num_per_node;d++)
					local_arr[d]=2147483647;
		}
	
	MPI_Barrier(MPI_COMM_WORLD);
	int r=1;
	int *local_comp=malloc((num_per_node)*sizeof(int));
	int *local_new_arr=malloc((num_per_node*2)*sizeof(int));
    	int recgs;
	int rcom=0;
	int sorted;
	int i;
	int count=0;
	int lnum=0;
	int rnum=0;
    	int gsorted=0;
    	qsort(local_arr,num_per_node,sizeof(int),compare);

	while(!gsorted){

		gsorted=1;
		sorted=0;//1 , continu . 0,out
//		printf("OK");
			if((even && r)||(odd && !r)){
//			printf("++_%d=_%d_%d_%d__",rank,local_arr[0],local_arr[1],local_arr[2]);

				if(rank != ROOT){
//				printf("rank=%d ",rank);
				MPI_Send(local_arr,num_per_node,MPI_INT,rank-1,0,MPI_COMM_WORLD);
//				 printf("++_%d=_%d_%d_%d__",rank,local_arr[0],local_arr[1],local_arr[2]);
				MPI_Recv(local_arr,num_per_node,MPI_INT,rank-1,0,MPI_COMM_WORLD,&status);
				}					
			}
			else if((odd && r) || (even && !r) ){
				if(rank != size-1){
//				printf("rankk=%d ",rank);
					MPI_Recv(local_comp,num_per_node,MPI_INT,rank+1,0,MPI_COMM_WORLD,&status);
//					printf("%d=_%d_%d_%d__",rank,local_comp[0],local_comp[1],local_comp[2]);
					sorted=Merge(local_new_arr,local_arr,num_per_node,local_comp,num_per_node);
			                MPI_Send(&local_new_arr[num_per_node],num_per_node,MPI_INT,rank+1,0,MPI_COMM_WORLD);
				}
				if(sorted)gsorted=0;
			}

//		printf(" rank(%d=%d,%d) ",rank,local_arr[0],local_new_arr[1]);
//	printf("yo");
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce (&gsorted,&recgs,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

		if(recgs==1) count=!(count^1)+1;
        	else count=0;

		if(count==2)gsorted=1;
		else gsorted=0;	
		r=!r;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_CREATE + MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);//write file
		if(remain !=0 && rank==size-1 && x>=0)
			num_per_node=x;
		else if(remain !=0 && x<0 && rank>=size-1-y)
			num_per_node=0;
		else if (remain !=0 && x<0 && rank==size-2-y)
			num_per_node=num_per_node-z;
		MPI_File_write_at(fh,offset,local_arr,num_per_node,MPI_INT,&status);
		MPI_File_close(&fh);
		MPI_File_close(&f);
	free(local_arr);
	MPI_Barrier(MPI_COMM_WORLD);
	cput=MPI_Wtime()-cput;
	if(rank==ROOT)printf("%lf",cput);
	MPI_Finalize();

	return 0;
}
