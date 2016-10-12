// Sample MPI Program for Parallel Programming 2015 Lab1
// You can reuse this code as a template for homework 1.
// Or you can write your own from scratch!
// It's up to you.

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define ROOT 0

void swap(int *local_arr1,int *local_arr2){
	int temp;
	temp=*local_arr1;
	*local_arr1=*local_arr2;
	*local_arr2=temp;

}
/*
void odevsort(int *local_arr,int num_per_node){
	int sorted = 0;
	int i;
	int j;
	while (!sorted) {
		sorted = 1;
		
                for(i = 1;i < num_per_node-1;i+=2){   //odd
                	if(local_arr[i] > local_arr[i+1]){
                        	swap(&local_arr[i],&local_arr[i+1]);
                                sorted = 0;
                        }
                }
                for(j = 0;j < num_per_node-1;j+=2){   //even
                	if(local_arr[j] > local_arr[j+1]){
                        	swap(&local_arr[j],&local_arr[j+1]);
                                sorted = 0;
                        }
                }
	}
}
*/
int odevsort(int *local_arr,int num_per_node,int r,int rank){
	int sorted = 0;
	int odd = rank & 1;
	int i= !(odd^r);
		for(i;i<num_per_node-1;i+=2)
			if(local_arr[i]>local_arr[i+1]){
				swap(&local_arr[i],&local_arr[i+1]);
				sorted=1;
			}
	return sorted;
}

int main (int argc, char *argv[]) { //argv is file to read
	int rank, size;	//rank is taskID of communicator
	printf("start");	//size is the number of processes in the group
	MPI_Status status;
	MPI_Init(&argc, &argv);//initializes MPI resource
	MPI_Comm_size(MPI_COMM_WORLD, &size);//get size by address
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//get taskID
//	MPI_Status status;
//	printf("\tMPI %d ",rank);
	if (argc < 3) {//argc must 3 , so less 3 that is error case
		if (rank == ROOT) {	//defined taskID==ROOT to be master
			fprintf(stderr, "Insufficient args\n");
			fprintf(stderr, "Usage: %s N input_file", argv[0]);
		}
		MPI_Barrier(MPI_COMM_WORLD);//all processes be Blocked until synchronized
		MPI_Finalize();
		return 0;
	}
//	printf("\twhat\n");
	const int N = atoi(argv[1]);//N is number of file size
	const char *inName = argv[2];//isName is file name

	int *root_arr; // for root process (which rank == 0) only
	
	// Part 1: Read file
	/* Note: You should deal with cases where (N < size) in Homework 1 */
	/*if (rank == ROOT) {
		root_arr = malloc(N * sizeof(int));//root allocate memory for file size
		FILE *f = fopen(inName, "rb");
		fread(root_arr, sizeof(int), N, f);//root read file
		fclose(f);
	}
	*/
	// Part 2: Scatter
	/* Note: You should deal with cases where (N % size != 0) in Homework 1 */
	int remain = N % size; //
	int odd = rank & 1;
	int even = !odd;
	int num_per_node = N / size;
	MPI_File f;
	MPI_Offset offset;
	if(remain>0)num_per_node++;
	//if (remain>0 && ((num_per_node+1)*size-N) < (num_per_node+1))num_per_node++;
	int *local_arr=malloc((num_per_node)*sizeof(int));
/*	if(rank == ROOT){
	
	root_arr = malloc((num_per_node)*size * sizeof(int));//root allocate memory for file size
	printf("\n%d\t%d\t\n",num_per_node,size);	
	     //File *f = fopen(inName, "rb");
	 }
*/	double t1=MPI_Wtime();
   	MPI_File_open(MPI_COMM_WORLD,inName,MPI_MODE_RDONLY,MPI_INFO_NULL,&f);
//	printf("\noh ya\n");
	offset=rank*num_per_node*sizeof(int);
//	MPI_File_set_view(f,offset,
	//if(num_per_node*size-N < num_per_node)
		MPI_File_read_at(f,offset,local_arr,num_per_node,MPI_INT,&status);//root read file
	double t_IO=MPI_Wtime()-t1;
	//printf("file OK");
		//fclose(f);
	/*if(((num_per_node+1)*size-N) > (num_per_node+1)){
		MPI_File_read_at(f,size*num_per_node*sizeof(int)+rank*sizeof(int),&local_arr[num_per_node],1,MPI_INT,&status);
		if(size*num_per_node+rank>=N)local_arr[num_per_node]=2147483647;
		num_per_node++;
	}*/
				//free(local_arr);
				//MPI_File_close(&f);
				//MPI_Finalize();
				//return 0;
		//	}else{
		//		size=N;
		//		remain=0;
		//	}
		//}
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
		
	//int *local_arr= malloc (num_per_node * sizeof(int));//allocate memory for recieve
//	printf("\t local");
	// Ref:
	// MPI_Scatter (send_buf, send_count, send_type, recv_buf, recv_count, recv_type, root, comm)
	//MPI_Scatter    (root_arr/*file data*/, num_per_node/*send number to node*/, MPI_INT/*integer type*/, local_arr/*receive*/, num_per_node/*receive number*/, MPI_INT/*receive type*/, ROOT/*sender*/, MPI_COMM_WORLD/*communicator*/);
	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Scatter(root_arr, num_per_node, MPI_INT, local_arr, num_per_node, MPI_INT, ROOT, MPI_COMM_WORLD);
//	if(rank==ROOT)printf("\n   %d scatter OK\n ",rank);
/*
*	if (rank == ROOT) {
*		free (root_arr);
*	}
*/
//	odevsort(local_arr,num_per_node);
	// Part 3: sorting ~ ~ ~ ~
	int r=0;
    	int gsorted; //all process sorted?
    	int comp=0;
    	int recgs;
	int rcom=0;
	int sorted;
	int num=num_per_node&1;
	int i;
	int count=0;
    	gsorted=0;
/*	for(i=0;i<num_per_node;i++)
		printf(" a[%d]=%d ",rank,local_arr[i]);    	
	printf("\n");
	
*/	double t_commu=0;
	while(!gsorted){

//	printf("  %d hi %d \n",rank,r);		
		
		
		gsorted=1;
	
		if(!num){
			if(!r){
			//printf(" r=oe %d ",rank);
				sorted=odevsort(local_arr,num_per_node,r,rank);
				rcom=local_arr[0];
				t1=MPI_Wtime();
				if(rank != ROOT)
					MPI_Send(&local_arr[0],1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
				if(rank != size-1){
					MPI_Recv(&comp,1,MPI_INT,rank+1,0,MPI_COMM_WORLD,&status);
					if(comp<local_arr[num_per_node-1])swap(&comp,&local_arr[num_per_node-1]);
					MPI_Send(&comp,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
				}
				if(rank != ROOT){
				MPI_Recv(&local_arr[0],1,MPI_INT,rank-1,0,MPI_COMM_WORLD,&status);            
				if(rcom != local_arr[0] || sorted) gsorted=0;
				}else if(sorted) gsorted=0;
				
				t_commu=MPI_Wtime()-t1+t_commu;
			}
			else{
				sorted=odevsort(local_arr,num_per_node,r,rank);
				if(sorted)gsorted=0;
			}
		}else{ 
	//printf("***%d***",rank);
			if((even && !r)||(odd && r)){
				sorted=odevsort(local_arr,num_per_node,r,rank);
				t1=MPI_Wtime();
				if(rank != ROOT){
				MPI_Send(&local_arr[0],1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
				comp=local_arr[0];
				MPI_Recv(&local_arr[0],1,MPI_INT,rank-1,0,MPI_COMM_WORLD,&status);
				
				if(comp!=local_arr[0] || sorted)  gsorted=0;
				}
				else if(sorted) gsorted=0;
				t_commu=MPI_Wtime()-t1+t_commu;
				
			}
			else if((odd && !r) || (even && r) ){
				sorted=odevsort(local_arr,num_per_node,r,rank);
                                t1=MPI_Wtime();
				if(rank != size-1){
                                MPI_Recv(&comp,1,MPI_INT,rank+1,0,MPI_COMM_WORLD,&status);
				//rcom=local_arr[num_per_node-1];
				if(comp<local_arr[num_per_node-1])swap(&comp,&local_arr[num_per_node-1]);
                                MPI_Send(&comp,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
				}
				if(sorted)gsorted=0;
				t_commu=MPI_Wtime()-t1+t_commu;
			}

		}
		
		//printf("%d ",rank);
		//MPI_Barrier(MPI_COMM_WORLD);
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		t1=MPI_Wtime();
		MPI_Allreduce (&gsorted,&recgs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		t_commu=MPI_Wtime()-t1+t_commu;
		//printf("\n oooo");
		if(recgs == size) count=!(count^1)+1;
        	else count=0;
		
		//printf("%d",recgs);
		if(count==2)gsorted=1;
		else gsorted=0;	
		r=!r;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	//free (local_arr);

	// Part 4: Accumulate the result to root process
	//int accumulated;
	// Ref:
	// MPI_Reduce (send_buf, recv_buf, count, data_type, op, root, comm)
	//MPI_Reduce    (&sum, &accumulated, 1, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);
	//MPI_Gather(&local_arr,num_per_node,MPI_INT,&root_arr,num_per_node,MPI_INT,ROOT,MPI_COMM_WORLD);
	//if (rank == ROOT) {
		t1=MPI_Wtime();
		MPI_File fh;
		MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_CREATE + MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);//write file
		if(remain !=0 && rank==size-1 && x>=0)
			num_per_node=x;
		else if(remain !=0 && x<0 && rank>=size-1-y)
			num_per_node=0;
		else if (remain !=0 && x<0 && rank==size-2-y)
			num_per_node=num_per_node-z;
		MPI_File_write_at(fh,offset,local_arr,num_per_node,MPI_INT,&status);
		//else
		//	MPI_File_write_at(fh,offset,local_arr,remain,MPI_INT,&status);
		MPI_File_close(&fh);
		MPI_File_close(&f);
		t_IO=MPI_Wtime()-t1+t_IO;
	//}
	free(local_arr);
	if(rank==ROOT){
		printf("IO=%f\n",t_IO);
		printf("Commu=%f",t_commu);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
