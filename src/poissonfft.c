#include "poissonfft.h"
#include <math.h>


void fft2d_init(FFTContext *c, double *y, double *b, int w, int h) {
	c->y = y;
	c->b = b;
	c->w = w;
	c->h = h;
	c->plan = fftw_plan_r2r_2d(h,w,b,b,FFTW_RODFT10,FFTW_RODFT00,0);
	c->planinv = fftw_plan_r2r_2d(h,w,y,y,FFTW_RODFT01,FFTW_RODFT00,0);	
}

void poissonfft2d(FFTContext *c) {
	fftw_execute(c->plan);
	for(int k = 0; k < c->h; k++) {
		for(int l = 0; l < c->w; l++) {
			c->y[k*c->w+l] = c->b[k*c->w+l]/(2*(cos(M_PI*(l+1.)/c->w)+cos(M_PI*(k+1.)/c->h)-2))/4/(c->w)/(c->h);
		}
	}
	fftw_execute(c->planinv);
}

void diffusionfft2d(FFTContext *c, double D, double dt, double dx) {
	fftw_execute(c->plan);
	double alpha = D*dt;
	for(int k = 0; k < c->h; k++) {
		for(int l = 0; l < c->w; l++) {
			c->y[k*c->w+l] = -0.9995*c->b[k*c->w+l]/(2*(alpha*cos(M_PI*(l+1.)/c->w)+alpha*cos(M_PI*(k+1.)/c->h)-2*(alpha+0.25*dx*dx)))/4/c->w/c->h*dx*dx;
		}
	}
	fftw_execute(c->planinv);
}

void fft2d_cleanup(FFTContext *c) {
	fftw_destroy_plan(c->plan);
	fftw_destroy_plan(c->planinv);
}

void poissonfft2d_test() {
	FFTContext c;
	int n = 10;
	double arr[n*n];
	fft2d_init(&c, arr,arr,n,n);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			arr[i*n+j] = (i*n+j*1.)/(4*(n+1)*(n+1));
			printf("%g ",arr[i*n+j]);
		}
	}
	printf("\n");
	fftw_execute(c.plan);
	fftw_execute(c.planinv);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			printf("%g ",arr[i*n+j]);
		}
	}
	printf("\n");
	fft2d_cleanup(&c);
}
