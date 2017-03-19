#include "fluid.h"
#include "cg.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void fluid_new(Fluid *f, int w, int h, double dx, double nu, int npart) {
	f->w = w;
	f->h = h;
	f->dx = dx;

	f->nu = nu;

	f->t = 0;

	f->v = calloc(2*w*h,sizeof(double));
	f->dvdt = calloc(2*w*h,sizeof(double));
	f->p = calloc(w*h,sizeof(double));
	f->tmp = calloc(5*w*h,sizeof(double));

	fft2d_init(&f->pressurefft,f->p,f->tmp+4*f->w*f->h,f->w-2,f->h-2);
	fft2d_init(&f->diffusionfft,f->tmp,f->tmp+f->w*f->h,f->w,f->h);

	// for an initial velocity distribution
	/*for(int y = 0; y < f->h; y++) {
		for(int x = 0; x < f->w; x++) {
			double cx = x - w/2;
			double cy = y - h/2;

			int idx = y*f->w+x;
			double r = sqrt(cx*cx + cy*cy+0.0001);
			//f->v[2*idx] = cy*exp(-r);
			//f->v[2*idx+1] = -cx*exp(-r);
		}
	}*/

	f->npart = npart;
	f->parts = calloc(npart,sizeof(Particle));
	for(int i = 0; i < f->npart; i++) {
		f->parts[i].x = (f->w-1)*dx*random()/RAND_MAX;
		f->parts[i].y = (f->h-1)*dx*random()/RAND_MAX;
	}
			
}

double force(int d, double x, double y, double t, int w, int h) {
	t/=2;

	// the commented out part is for some rotating jet which also produces a nice pattern
	/*double cx = x - w/2+w/6*sin(t/5);
	double cy = y - h/2+w/6*cos(t/5);

	double r = sqrt((cx*cx+cy*cy))/w;
	//return 10*(d == 1 ? cx/w : -cy/h)*sin(12*3*r)*exp(-0.05*(cx*cx+cy*cy));
	if(d == 0)
		return 400*cos(3*t)*exp(-0.1*(cx*cx+cy*cy));
	else 
		return 400*sin(3*t)*exp(-0.1*(cx*cx+cy*cy));
	*/

	double cx = x-w/2;
	double cy = y-h/4;
	if(d == 1)
		return 400*exp(-0.1*(cx*cx+cy*cy));
	return 0;

}


void open_bc(Fluid *f) {
	for(int y = 1; y < f->h-1; y++) {
		f->v[2*y*f->w] = f->v[2*(y*f->w+1)];
		f->v[2*y*f->w+1] = f->v[2*(y*f->w+1)+1];
		f->v[2*(y*f->w+f->w-1)] = f->v[2*(y*f->w+f->w-2)];
		f->v[2*(y*f->w+f->w-1)+1] = f->v[2*(y*f->w+f->w-2)+1];
	}
	for(int x = 1; x < f->w-1; x++) {
		f->v[2*x] = f->v[2*(x+f->w)];
		f->v[2*x+1] = f->v[2*(x+f->w)+1];
		f->v[2*(x+f->w*(f->h-1))] = f->v[2*(x+f->w*(f->h-2))];
		f->v[2*(x+f->w*(f->h-1))+1] = f->v[2*(x+f->w*(f->h-2))+1];
	}

}

void closed_bc(Fluid *f) {
#pragma omp parallel for
	for(int y = 0; y < f->h; y++) {
		for(int x = 0; x < f->w; x++) {
			int idx = y * f->w+x;
			const int border = 5;
			if(x > border && y > border && x < f->w-border && y < f->h-border)
				continue;
			f->v[2*idx] = 0;
			f->v[2*idx+1] = 0;
		}
	}

}
void periodic_bc(Fluid *f) {
	for(int y = 0; y < f->h; y++) {
		for(int x = 0; x < f->w; x++) {
			int idx = y * f->w+x;
			const int border = 5;
			if(x > border && y > border && x < f->w-border && y < f->h-border)
				continue;
			int effw = f->w-2*border, effh = f->h-2*border;
			int nx=x, ny=y;
			if(x <= border)
				nx+=effw;
			if(x >= f->w-border)
				nx-=effw;
			if(y <= border)
				ny+=effh;
			if(y >= f->h-border)
				ny-=effh;
			int idx2 = ny * f->w+nx;
			f->v[2*idx] = f->v[2*idx2];
			f->v[2*idx+1] = f->v[2*idx2+1];
		}
	}
}


void fluid_step(Fluid *f, double dt) {
	periodic_bc(f);

	poissonfft2d(&f->pressurefft);
	// calculate dv/dt
	for(int y = 1; y < f->h-1; y++) {
		for(int x = 1; x < f->w-1; x++) {
			int idx = y*f->w+x;
			for(int d = 0; d < 2; d++) {
				f->dvdt[2*idx+d] =
					-f->v[2*idx]*(f->v[2*(idx+1)+d]-f->v[2*(idx-1)+d])/2/f->dx
					-f->v[2*idx+1]*(f->v[2*(idx+f->w)+d]-f->v[2*(idx-f->w)+d])/2/f->dx
					+force(d,x,y,f->t,f->w,f->h);
			}
		}
	}
	periodic_bc(f);
	// explicit iteration for the convection term
	#pragma omp parallel for
	for(int y = 1; y < f->h-1; y++) {
		for(int x = 1; x < f->w-1; x++) {
			int idx = y*f->w+x;
			f->v[2*idx] += dt*0.25*(f->dvdt[2*(idx+1)]+f->dvdt[2*(idx-1)]+f->dvdt[2*(idx+f->w)]+f->dvdt[2*(idx-f->w)]);
			f->v[2*idx+1] += dt*0.25*(f->dvdt[2*(idx+1)+1]+f->dvdt[2*(idx-1)+1]+f->dvdt[2*(idx+f->w)+1]+f->dvdt[2*(idx-f->w)+1]);
			
			double vx = f->v[2*idx];
			double vy = f->v[2*idx+1];
			double v = sqrt(vx*vx+vy*vy);
			double vmax = f->dx/dt/10;
			if(v > vmax) {
				f->v[2*idx] *= vmax/v;
				f->v[2*idx+1] *= vmax/v;
			}
		}
	}
	periodic_bc(f);

	// implicit diffusion step (could probably be more efficient using operator splitting)
	int offset = f->w*f->h;
	
	#pragma omp parallel for
	for(int y = 0; y < f->h; y++) {
		for(int x = 0; x < f->w; x++) {
			int idx = y*f->w+x;
			f->tmp[idx+offset] = f->v[2*idx];
		}
	}
	diffusionfft2d(&f->diffusionfft,f->nu,dt,f->dx);
	#pragma omp parallel for
	for(int y = 0; y < f->h; y++) {
		for(int x = 0; x < f->w; x++) {
			int idx = y*f->w+x;
			f->v[2*idx] = f->tmp[idx];
			f->tmp[idx+offset] = f->v[2*idx+1];
		}
	}
	diffusionfft2d(&f->diffusionfft,f->nu,dt,f->dx);
	#pragma omp parallel for
	for(int y = 0; y < f->h; y++) {
		for(int x = 0; x < f->w; x++) {
			int idx = y*f->w+x;
			f->v[2*idx+1] = f->tmp[idx];
		}
	}
	periodic_bc(f);


	fluid_remove_divergence(f);
	
	// advance particles
#pragma omp parallel for
	for(int i = 0; i < f->npart; i++) {
		int x = f->parts[i].x/f->dx+0.5;
		int y = f->parts[i].y/f->dx+0.5;

		f->parts[i].x += dt*f->v[2*(y*f->w+x)];
		f->parts[i].y += dt*f->v[2*(y*f->w+x)+1];

		if(f->parts[i].x < 0)
			f->parts[i].x = (f->w-2)*f->dx;
		if(f->parts[i].y < 0)
			f->parts[i].y = (f->h-2)*f->dx;
		if(f->parts[i].x > (f->w-2)*f->dx)
			f->parts[i].x = 2*f->dx;
		if(f->parts[i].y > (f->h-2)*f->dx)
			f->parts[i].y = 2*f->dx;
	}
	
	
	f->t += dt;
	periodic_bc(f);
}

void fluid_remove_divergence(Fluid *f) {
	int offset = f->w*f->h;
	#pragma omp parallel for
	for(int y = 1; y < f->h-1; y++) {
		for(int x = 1; x < f->w-1; x++) {
			int idx = y*f->w+x;
			int idx2 = (y-1)*(f->w-2)+x-1;
			f->tmp[4*offset+idx2] = 0.5*f->dx*(f->v[2*(idx+1)]-f->v[2*(idx-1)]
				+f->v[2*(idx+f->w)+1]-f->v[2*(idx-f->w)+1]);
		}
	}
	poissonfft2d(&f->pressurefft);
	// add pressure contribution
	#pragma omp parallel for
	for(int y = 2; y < f->h-2; y++) {
		for(int x = 2; x < f->w-2; x++) {
			int idx = y*f->w+x;
			int idx2 = (y-1)*(f->w-2)+x-1;
			f->v[2*idx] -= (f->p[idx2+1]-f->p[idx2-1])/2/f->dx;
			f->v[2*idx+1] -= (f->p[idx2+(f->w-2)]-f->p[idx2-(f->w-2)])/2/f->dx;
		}
	}
}

void fluid_cleanup(Fluid *f) {
	free(f->v);
	free(f->dvdt);
	free(f->p);
	free(f->tmp);
	free(f->parts);
	fft2d_cleanup(&f->pressurefft);
}
