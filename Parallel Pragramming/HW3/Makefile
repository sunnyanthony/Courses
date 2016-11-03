all:	mpis mpid omps ompd hybs hybd
mpis:
	mpicc -o MS_MPI_static -lX11 MS_MPI_static.c
mpid:
	mpicc -o MS_MPI_dynamic -lX11 MS_MPI_dynamic.c
omps:
	gcc -o MS_OpenMP_static -fopenmp -lX11 MS_OpenMP_static.c
ompd:
	gcc -o MS_OpenMP_dynamic -fopenmp -lX11 MS_OpenMP_dynamic.c
hybs:
	mpicc -o MS_Hybrid_static -fopenmp -lX11 MS_Hybrid_static.c
hybd:
	mpicc -o MS_Hybrid_dynamic -fopenmp -lX11 MS_Hybrid_dynamic.c
