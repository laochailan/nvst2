#include "cg.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// taken from https://de.wikipedia.org/wiki/CG-Verfahren
int conjugate_gradient(double *x0, void (*op)(double *res, double *x, int n, void *arg), double *b, int n, void *arg, double tol, int maxiter, double *tmp1, double *tmp2) {
	double *r = b;
	double *z = tmp1;
	double *d = tmp2;
	memset(z,0,sizeof(double)*n);

	op(z,x0,n,arg);
	
	double resnorm = 0;
	double bnorm = 0;
	for(int i = 0; i < n; i++) {
		bnorm += b[i]*b[i];
		r[i] -= z[i];
		d[i] = r[i];

		resnorm += r[i]*r[i];
	}
	
	int step = 0;
	for(; step < maxiter && sqrt(resnorm/bnorm) > tol; step++) {
		op(z,d,n,arg);

		double dz = 0;
		for(int i = 0; i < n; i++) {
			dz += d[i]*z[i];
		}

		double alpha = resnorm/dz;
		double newresnorm = 0;
		dz = 0;
		#pragma omp parallel for
		for(int i = 0; i < n; i++) {
			x0[i] += alpha*d[i];
			r[i] -= alpha*z[i];
			newresnorm += r[i]*r[i];
			dz = r[i]*z[i];
		}
		double beta = -alpha*dz/resnorm;
		#pragma omp parallel for
		for(int i = 0; i < n; i++) {
			d[i] = r[i]+beta*d[i];
		}
		//if(step > maxiter/2 && newresnorm > resnorm)
		//	return step;
		resnorm = newresnorm;
		//printf("%d: %g\n",step,sqrt(resnorm/bnorm));
	}
	printf("cg - %d: %g\n",step,sqrt(resnorm/bnorm));

	//assert(step < maxiter);
	return step;
}

void laplace2D(double *res, double *f, int n, void *arg) { 
	struct laplace2Dargs *a = (struct laplace2Dargs*) arg;

	int w = a->w;
	int h = n/w;
	assert(w*h == n);

#pragma omp parallel for
	for(int y = 1; y < h-1; y++) {
		for(int x = 1; x < w-1; x++) {
			res[y*w+x] = f[y*w+x+1]+f[y*w+x-1]+f[(y+1)*w+x]+f[(y-1)*w+x]-4*f[y*w+x];
		}
	}
}

void laplace2D5(double *res, double *f, int n, void *arg) { 
	struct laplace2Dargs *a = (struct laplace2Dargs*) arg;

	int w = a->w;
	int h = n/w;
	assert(w*h == n);

	for(int y = 2; y < h-2; y++) {
		for(int x = 2; x < w-2; x++) {
			double l2 = -(f[y*w+x+2]+f[y*w+x-2]+f[(y+2)*w+x]+f[(y-2)*w+x]);
			double l1 = 16*(f[y*w+x+1]+f[y*w+x-1]+f[(y+1)*w+x]+f[(y-1)*w+x]);
			double l0 = -2*30*f[y*w+x];
			res[y*w+x] = (l2+l1+l0)/12;
		}
	}
}


static void mat_op(double *res, double *x, int n, void *arg) {
	double *A = (double *)arg;
	for(int i = 0; i < n; i++) {
		res[i] = 0;
		for(int j = 0; j < n; j++) {
			res[i] += A[i*n+j]*x[j];
		}
	}
}

void test_conjugate_gradient2(void) {
	int n = 2;
	
	double b[] = {1,2};
	double A[] = {4,1,1,3};
	double x0[] = {2,1};
	double tmp1[2],tmp2[2];
	
	int step = conjugate_gradient(x0,mat_op,b,n,A,1e-8,100,tmp1,tmp2);
	printf("%d steps, (%g, %g)\n",step,x0[0],x0[1]);
}

void test_conjugate_gradient(void) {
	int w = 101;
	int h = 101;

	double *b = calloc(w*h,sizeof(double));
	double *x0 = calloc(w*h,sizeof(double));
	double *tmp1 = calloc(w*h,sizeof(double));
	double *tmp2 = calloc(w*h,sizeof(double));

	for(int y = 0; y < h; y++) {
		for(int x = 0; x < w; x++) {
			x0[y*w+x] = sin(M_PI*x/(w-1))*sin(5*M_PI*y/(h-1))*0.0;
			b[y*w+x] = -2*M_PI*M_PI/(w-1)/(h-1)*sin(M_PI*x/(w-1))*sin(M_PI*y/(h-1));
		}
	}


	struct laplace2Dargs a = { w };
	
	int step = conjugate_gradient(x0,laplace2D,b,w*h,&a,1e-15,6*w*h,tmp1,tmp2);

	double err = 0;
	for(int y = 0; y < h; y++) {
		for(int x = 0; x < w; x++) {
			double cx = x - w/2;
			double cy = y - h/2;

			if(cx*cx+cy*cy > 4)
				err += fabs(x0[y*w+x]-sin(M_PI*x/(w-1))*sin(M_PI*y/(h-1)));
		}
	}
	err/=w*h;

	printf("CG test:\nConvergence in %d steps.\nError: %g\n",step,err);
	/*for(int x = 0; x < w; x++) {
		double cx = x - w/2;
		double cy = 0;
		printf("%d: %g %g\n",x,x0[h/2*w+x], sin(M_PI*x/(w-1))*sin(M_PI*h/2/(h-1)));
	}*/

/*
	for(int y = 0; y < h; y++) {
		for(int x = 0; x < w; x++) {
			double cx = x - w/2;
			double cy = y - h/2;
			x0[y*w+x] = dx*dx*(cx*cx+cy*cy);
		}
	}
	laplace2D5(b,x0,w*h,&a);
	for(int x = 0; x < w; x++) {
		double cx = x - w/2;
		double cy = 0;
		printf("%d: %g\t%g\n",x, b[h/2*w/2+x],sin(20*dx*x));
	}*/
}

