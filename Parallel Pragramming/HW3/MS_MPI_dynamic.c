 /* 
   MPI Dynamic Mandelbort sort
 */
#include <stdlib.h>
#include <mpi.h>
#include <X11/Xlib.h>
#include <stdio.h>

#define data_tag 100
#define result_tag 101
#define terminate_tag 110

typedef struct complextype
{
	double real, imag;
} Compl;

int main(int argc, char *argv[])
{
// ----------input---------------
	int size, rank; 
	double leftR = atof(argv[2]);
	double rightR = atof(argv[3]);
	double lowerR = atof(argv[4]);
	double upperR = atof(argv[5]);
	int width = atoi(argv[6]);
	int height = atoi(argv[7]);
	char *en = argv[8];
	int num_r;
	int dis;
	if( *en == 'e')dis = 1;
	else dis = 0;

// ------------MPI initial----------------
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size (MPI_COMM_WORLD, &size); 
	MPI_Comm_rank (MPI_COMM_WORLD, &rank); 
	MPI_Barrier(MPI_COMM_WORLD);

	int t=0;
// ------- size =1
	if(size == 1){
		Display *display;
		Window window;      //initialization for a window
		int screen;         //which screen 
		/* open connection with the server */ 
		GC gc;
		XGCValues values;
		long valuemask = 0;
		if(dis == 1){
			display = XOpenDisplay(NULL);
			if(display == NULL) {
				fprintf(stderr, "cannot open display\n");
				exit(1);
			}
    
			screen = DefaultScreen(display);
    
			/* set window position */
			int x = 0;
			int y = 0;
    
			/* border width in pixels */
			int border_width = 0;
    
			/* create window */
			window = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, width, height, border_width,
							BlackPixel(display, screen), WhitePixel(display, screen));
			
			/* create graph */
			
			
			gc = XCreateGC(display, window, valuemask, &values);
			//XSetBackground (display, gc, WhitePixel (display, screen));
			XSetForeground (display, gc, BlackPixel (display, screen));
			XSetBackground(display, gc, 0X0000FF00);
			XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
			
			/* map(show) the window */
			XMapWindow(display, window);
			XSync(display, 0);
		}
	
		/* draw points */
		Compl z, c;
		int repeats;
		double temp, lengthsq;
		int i, j;
		for(i=0; i<width; i++) {
			for(j=0; j<height; j++) {
				z.real = 0.0;
				z.imag = 0.0;
				c.real = leftR + (double)((double)i *((rightR-leftR)/(double)width)); 
				c.imag = lowerR + (double)((double)j * ((upperR-lowerR)/(double)height));
				repeats = 0;
				lengthsq = 0.0;
		
				while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
					temp = z.real*z.real - z.imag*z.imag + c.real;
					z.imag = 2*z.real*z.imag + c.imag;
					z.real = temp;
					lengthsq = z.real*z.real + z.imag*z.imag; 
					repeats++;
				}
				if(dis == 1){
					XSetForeground (display, gc,  1024 * 1024 * (repeats % 256));		
					XDrawPoint (display, window, gc, i, j);
				}
			}
		}
		if(dis == 1)XFlush(display);
	}else{
// --------------display-------------
		GC gc;
		Display *display;
		Window window;
		if(rank==0){
			if(dis){
			//Display *display;
			//Window window;      //initialization for a window
			int screen;         //which screen 
		
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
			XGCValues values;
			long valuemask = 0;
			
			gc = XCreateGC(display, window, valuemask, &values);
			//XSetBackground (display, gc, WhitePixel (display, screen));
			XSetForeground (display, gc, BlackPixel (display, screen));
			XSetBackground(display, gc, 0X0000FF00);
			XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
			
			/* map(show) the window */
			XMapWindow(display, window);
			XSync(display, 0);
			} // display-----------------------------------------
		
			/* draw points */
			//Compl z, c;
			//int repeats;
			//double temp, lengthsq;
			int i;
			int count = 0;
			int *slave = malloc((width+1)*sizeof(int));
			num_r = 0;
			//send chunk = 1
			for(i=0;i<size-1;i++){
				MPI_Send(&num_r,1,MPI_INT,i+1,data_tag,MPI_COMM_WORLD);
				count++;
				num_r++;
			}
		
			while(count >0){
				MPI_Recv(slave,width+1,MPI_INT,MPI_ANY_SOURCE,result_tag,MPI_COMM_WORLD,&status);
				count--;
				if(dis == 1){ // ---------draw row
							for(i=0; i<width; i++){
								XSetForeground (display, gc,  1024 * 1024 * (slave[i+1] % 256));
								XDrawPoint (display, window, gc, i, *slave);               			
					}
				}
				if(num_r<height){
					MPI_Send(&num_r,1,MPI_INT,status.MPI_SOURCE,data_tag,MPI_COMM_WORLD);//CHUNK=1
					count++;
					num_r++;
				}else{
					MPI_Send(&num_r,1,MPI_INT,status.MPI_SOURCE,terminate_tag,MPI_COMM_WORLD);//CHUNK=1
				}
				
			}
			if(dis == 1)XFlush(display);
		
		}else{
			Compl z,c;
			int repeats;
			double temp, lengthsq;
			int i, j;
			int *Pi = malloc((width+1)*sizeof(int));
			MPI_Recv(&num_r,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status); //
			
			while(status.MPI_TAG == data_tag){
				*Pi = num_r; //Pi[0] save this cpu do number row
				c.imag = lowerR + (double)((double)num_r * ((upperR-lowerR)/(double)height)); // imag is y axis
				for(i=0; i<width; i++) {
				//	for(j=0; j<height; j++) {
					t++;
						z.real = 0.0;
						z.imag = 0.0;
						c.real = leftR + (double)((double)i *((rightR-leftR)/(double)width)); /* Theorem : If c belongs to M(Mandelbrot set), then |c| <= 2 */
						//c.imag = ((double)j - 400.0)/200.0; /* So needs to scale the window */
						repeats = 0;
						lengthsq = 0.0;
			
						while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
							temp = z.real*z.real - z.imag*z.imag + c.real;
							z.imag = 2*z.real*z.imag + c.imag;
							z.real = temp;
							lengthsq = z.real*z.real + z.imag*z.imag; 
							repeats++;
						}
					Pi[i+1] = repeats;
			//			XSetForeground (display, gc,  1024 * 1024 * (repeats % 256));		
			//			XDrawPoint (display, window, gc, i, j);
				}
				MPI_Send(Pi,width+1,MPI_INT,0,result_tag,MPI_COMM_WORLD);
				MPI_Recv(&num_r,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			}
		}
	}
//	XFlush(display);

		//t=MPI_Wtime()-t;
			printf("mpis[%d]=%d\n",rank, t);

	sleep(5);
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}
