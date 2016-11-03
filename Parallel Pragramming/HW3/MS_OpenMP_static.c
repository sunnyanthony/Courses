 /* 
 *    ompd Mandelbort sort
 *     
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
typedef struct complexType
{
	double real, imag;
} Compl;

int main(int argc, char *argv[])
{
// ----------input---------------
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
	if( *en == 'e')dis = 1;
	else dis = 0;
	double t_start[24], total[24]={0};
	double t = omp_get_wtime();
	double t_win = omp_get_wtime();
// --------display----------------------	
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

	t_win = omp_get_wtime()-t_win;
	omp_set_num_threads(number_thread);
	#pragma omp parallel default(shared) private(j, z, c, repeats, lengthsq, temp) 
	{
		t_start[omp_get_thread_num()] = 0;
		#pragma omp for schedule(static)
		for(i=0; i<width; i++) {
			t_start[omp_get_thread_num()] = omp_get_wtime();
			c.real = leftR + (double)((double)i *((rightR-leftR)/(double)width)); 
			for(j=0; j<height; j++) {
				z.real = 0.0;
				z.imag = 0.0;
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
				
					if(dis==1){
						#pragma omp critical 
						{
						XSetForeground (display, gc,  1024 * 1024 * (repeats % 256));
						XDrawPoint (display, window, gc, i, j);
						}
					}
			}
			total[omp_get_thread_num()]+=(omp_get_wtime()-t_start[omp_get_thread_num()]);
		}
	
	}

	for(i=0;i<number_thread;i++){
        printf("%d %f\n",i,total[i]+t_win);
    }
    t = omp_get_wtime()-t;	
    //printf("omps=%f\n",t);
    
	if(dis==1)	
		XFlush(display);
	sleep(5);	
	return 0;
}
	
