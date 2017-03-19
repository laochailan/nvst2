#ifndef POISSONFFT_H
#define POISSONFFT_H

#include <fftw3.h>

// Solve poisson equation using FFT (like in Numerical Recipes)

struct FFTContext {
	int w;
       	int h;
	double *y;
	double *b;
	fftw_plan plan;
	fftw_plan planinv;
};
typedef struct FFTContext FFTContext;


void fft2d_init(FFTContext *c, double *y, double *b, int w, int h);

// poissonfft2d solves the equation Δy = b with dirichlet boundary conditions.
void poissonfft2d(FFTContext *c);

// diffusionfft2d does a diffusion timestep where y = u(t+dt) and b = u(t).
// D is the diffusion constant, dt the timestep and dx the space step.
void diffusionfft2d(FFTContext *c, double D, double dt, double dx);

void fft2d_cleanup(FFTContext *c);

void poissonfft2d_test();

#endif
