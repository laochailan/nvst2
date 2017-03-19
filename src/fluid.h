#ifndef PHYS_H
#define PHYS_H

#include "poissonfft.h"

struct Particle {
	double x;
	double y;
};
typedef struct Particle Particle;

struct Fluid {
	int w; // number of horizontal cells
	int h; // number of vertical cells

	double t; // time

	double dx; // grid spacing

	double nu; // dynamic viscosity
	
	double *v; // 2*w*h velocity field
	double *dvdt; // 2*w*h change in velocity field
	double *p; // w*h pressure field
	
	double *tmp; // n*w*h auxiliary arrays

	int npart;
	Particle *parts;
	FFTContext pressurefft;
	FFTContext diffusionfft;
};
typedef struct Fluid Fluid;

void fluid_new(Fluid *f, int w, int h, double dx, double nu, int npart);
void fluid_step(Fluid *f, double dt);
void fluid_remove_divergence(Fluid *f);
void fluid_cleanup(Fluid *f);

#endif
