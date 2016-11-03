 /* 
   Hybrid static Mandelbort sort
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <mpi.h>


struct timeval tv1, tv2;
double t = 0;
typedef struct complextype
{
	double real, imag;
} Compl;

int main(int argc, char *argv[])
{
	// ----------input---------------
	int size, rank; 
	int number_thread = atoi(argv[1]);
	double leftR = atof(argv[2]);
	double rightR = atof(argv[3]);
	double lowerR = atof(argv[4]);
	double upperR = atof(argv[5]);
	int width = atoi(argv[6]);
	int height = atoi(argv[7]);
	char *en = argv[8];
	int num_r;
	int dis;
	int remainder=0;
	if( *en == 'e')dis = 1;
	else dis = 0;
	
// ------------MPI initial----------------
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size (MPI_COMM_WORLD, &size); 
	MPI_Comm_rank (MPI_COMM_WORLD, &rank); 
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank == 0)gettimeofday(&tv1, NULL);
	double t_start[24], total[24]={0};
	double t=0;
	double t_win = omp_get_wtime();
// -----------node of number row
	if(height % size == 0)
		num_r = height / size;
	else
		num_r = height / size +1;
	
	if(rank == size-1 && height % size!=0)
		remainder = height % size;
	
// ------display && root-----
	if(rank==0){
		Display *display;
		Window window;      //initialization for a window
		int screen;         //which screen 
		GC gc;
		XGCValues values;
		long valuemask = 0;
		if(dis ==1){
			/* open connection with the server */ 
			display = XOpenDisplay(NULL);
			if(display == NULL) {
				fprintf(stderr, "cannot open display\n");
				return 0;
			}
		
			screen = DefaultScreen(display);
		
			/* set window size */
			//int width = 800;
			//int height = 800;
		
			/* set window position */
			int x = 0;
			int y = 0;
		
			/* border width in pixels */
			int border_width = 0;
		
			/* create window */
			window = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, width, height, border_width,
							BlackPixel(display, screen), WhitePixel(display, screen));
			
			/* create graph */
			//GC gc;
			//XGCValues values;
			//long valuemask = 0;
			
			gc = XCreateGC(display, window, valuemask, &values);
			//XSetBackground (display, gc, WhitePixel (display, screen));
			XSetForeground (display, gc, BlackPixel (display, screen));
			XSetBackground(display, gc, 0X0000FF00);
			XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
			
			/* map(show) the window */
			XMapWindow(display, window);
			XSync(display, 0);
		}
		
		
		Compl z, c;
		int repeats;
		double temp, lengthsq;
		int i, j, k;
		int Pi[num_r][width];
		omp_set_num_threads(number_thread);
		#pragma omp parallel default(shared) private(j, z, c, repeats, lengthsq, temp) 
		{
			#pragma omp for schedule(auto)
			for(i=0; i<width; i++) {
				for(j=rank*num_r; j<rank*num_r+num_r && j<height; j++) {
					z.real = 0.0;
					z.imag = 0.0;
					c.real = leftR + (double)((double)i *((rightR-leftR)/(double)width)); /* Theorem : If c belongs to M(Mandelbrot set), then |c| <= 2 */
					c.imag = lowerR + (double)((double)j * ((upperR-lowerR)/(double)height)); // imag is y axis
					repeats = 0;
					lengthsq = 0.0;
					
					while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
						temp = z.real*z.real - z.imag*z.imag + c.real;
						z.imag = 2*z.real*z.imag + c.imag;
						z.real = temp;
						lengthsq = z.real*z.real + z.imag*z.imag; 
						repeats++;
					}
					//Pi[j-(rank*num_r)][i] = repeats;
					
					if(dis ==1 ){
						#pragma omp critical
						{
							XSetForeground (display, gc,  1024 * 1024 * (repeats% 256));
							XDrawPoint (display, window, gc, i, j); 
						}
					}
				}
			}				
		}
				if(dis == 1)
			XFlush(display);
		
		printf("this");
		// master recive
		#pragma omp critical
		{
			for(k=0;k<size-1;k++){
				MPI_Recv(Pi,num_r*width,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				
				for(i=(status.MPI_SOURCE*num_r);i<(status.MPI_SOURCE*num_r+num_r) && i<height ;i++)
					for(j=0;j<width;j++){
						if(dis == 1){
							XSetForeground (display, gc,  1024 * 1024 * (Pi[i-status.MPI_SOURCE*num_r][j]% 256));
							XDrawPoint (display, window, gc, j, i);
							//printf("\n n %d",status.MPI_SOURCE);
						}
					}
			}
		}
		t_win = omp_get_wtime()-t_win;
		if(dis == 1)
			XFlush(display);
		
	}else{

		/* draw points */
		Compl z, c;
		int repeats;
		double temp, lengthsq;
		int i, j;
		int Pi[num_r][width];
		t_win = omp_get_wtime()-t_win;
		omp_set_num_threads(number_thread);
		#pragma omp parallel default(shared) private(i,j, z, c, repeats, lengthsq, temp) 
		{
			t_start[omp_get_thread_num()] = 0;
			#pragma omp for schedule(static)
			for(j=0; j<num_r ; j++) {
				t_start[omp_get_thread_num()] = omp_get_wtime();
				if(j+rank*num_r<height){
					for(i=0; i<width; i++) {
						z.real = 0.0;
						z.imag = 0.0;
						c.real = leftR + (double)((double)i *((rightR-leftR)/(double)width)); /* Theorem : If c belongs to M(Mandelbrot set), then |c| <= 2 */
						c.imag = lowerR + (double)((double)(j+rank*num_r )* ((upperR-lowerR)/(double)height)); // imag is y axis
						repeats = 0;
						lengthsq = 0.0;
						
						while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
							temp = z.real*z.real - z.imag*z.imag + c.real;
							z.imag = 2*z.real*z.imag + c.imag;
							z.real = temp;
							lengthsq = z.real*z.real + z.imag*z.imag; 
							repeats++;
						}
						Pi[j][i] = repeats;
					}
				}else
					Pi[j][i]=0;
				total[omp_get_thread_num()]+=(omp_get_wtime()-t_start[omp_get_thread_num()]);
			}
			//printf("q=%d \n",rank);				
		}


		#pragma omp critical
		{	
			MPI_Send(Pi,num_r*width , MPI_INT,0,rank,MPI_COMM_WORLD);
		}
	}
		int i;
		for(i=0;i<number_thread;i++){
              printf("%d %f\n",i+rank*number_thread,total[i]+t_win);
        }
	if(rank == 0){
	gettimeofday(&tv2, NULL);
	t += (double)(tv2.tv_sec - tv1.tv_sec)+(double)(tv2.tv_usec - tv1.tv_usec)/1000000.0;
	//printf("hys=%lf\n", t);
	}

	sleep(5);	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	

	return 0;
}
