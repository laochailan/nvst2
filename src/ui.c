#include "ui.h"
#include <SDL.h>

int ui_init(UI *ui, int w, int h, int scale) {
	ui->w = w;
	ui->h = h;
	ui->scale = scale;
	ui->mode = DISPLAY_VELDIR;

	if(SDL_Init(SDL_INIT_VIDEO) != 0)
		return 1;

	ui->win = SDL_CreateWindow("nvst2",SDL_WINDOWPOS_UNDEFINED,SDL_WINDOWPOS_UNDEFINED,w*scale,h*scale,0);
	if(ui->win == 0) {
		SDL_Quit();
		return 2;
	}

	ui->ren = SDL_CreateRenderer(ui->win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
	if(ui->ren == 0) {
		SDL_DestroyWindow(ui->win);
		SDL_Quit();
		return 3;
	}

	ui->tex = SDL_CreateTexture(ui->ren,SDL_PIXELFORMAT_RGB24,SDL_TEXTUREACCESS_STREAMING,w,h);
	if(ui->tex == 0) {
		SDL_DestroyWindow(ui->win);
		SDL_DestroyRenderer(ui->ren);
		SDL_Quit();
		return 3;
	}


	return 0;
}

	
int ui_cleanup(UI *ui) {
	SDL_DestroyTexture(ui->tex);
	SDL_DestroyRenderer(ui->ren);
	SDL_DestroyWindow(ui->win);
	SDL_Quit();
	return 0;
}

static double map(double x) {
	return 0.5+0.5*tanh(x);
}


static void add_clr(int clr, int x, int y, int pitch, uint8_t *pixels) {
	for(int d = 0; d < 3; d++) {
		if(pixels[y*pitch+3*x+d]>255-clr)
			pixels[y*pitch+3*x+d] = 255;
		else 
			pixels[y*pitch+3*x+d] += clr;
	}
}

static void draw_particle(double x, double y, double clr, int pitch, uint8_t *pixels) {
	if(!isnormal(x) || !isnormal(y))
		return;

	int fx = x+0.5;
	int fy = y+0.5;
	double rx = x+0.5-fx;
	double ry = y+0.5-fy;
	add_clr(clr*rx*ry,fx+1,fy+1,pitch,pixels);
	add_clr(clr*(1-rx)*ry,fx,fy+1,pitch,pixels);
	add_clr(clr*(1-rx)*(1-ry),fx,fy,pitch,pixels);
	add_clr(clr*rx*(1-ry),fx+1,fy,pitch,pixels);
}

void ui_draw_pixels(Fluid *f, uint8_t *pixels, int pitch, DisplayMode mode) {
	for(int y = 0; y < f->h; y++) {
		for(int d = 0; d < 2; d++) {
			pixels[y*pitch+3*d*(f->w-1)+0] = 0;
			pixels[y*pitch+3*d*(f->w-1)+1] = 0;
			pixels[y*pitch+3*d*(f->w-1)+2] = 0;
		}
	}
	for(int d = 0; d < 2; d++) {
		for(int x = 0; x < f->w; x++) {
			pixels[d*(f->h-1)*pitch+0] = 0;
			pixels[d*(f->h-1)*pitch+1] = 0;
			pixels[d*(f->h-1)*pitch+2] = 0;
		}
	}

	if(mode == DISPLAY_VELDIR) {
	#pragma omp parallel for
		for(int y = 0; y < f->h; y++) {
			for(int x = 0; x < f->w; x++) {
				double vx = f->v[2*(y*f->w+x)];
				double vy = f->v[2*(y*f->w+x)+1];
				double vr = sqrt(vx*vx+vy*vy);
				pixels[y*pitch+3*x+0] = 255*map(vx/vr);
				pixels[y*pitch+3*x+1] = 255*map(vy/vr);
				pixels[y*pitch+3*x+2] = 128;
			}
		}
	} else if(mode == DISPLAY_VEL) {
	#pragma omp parallel for
		for(int y = 0; y < f->h; y++) {
			for(int x = 0; x < f->w; x++) {
				double vx = f->v[2*(y*f->w+x)];
				double vy = f->v[2*(y*f->w+x)+1];
				pixels[y*pitch+3*x+0] = 255*map(2*vx);
				pixels[y*pitch+3*x+1] = 255*map(2*vy);
				pixels[y*pitch+3*x+2] = 128;
			}
		}
	} else if(mode == DISPLAY_PART) {
	#pragma omp parallel for
		for(int y = 0; y < f->h; y++) {
			for(int x = 0; x < f->w; x++) {
				pixels[y*pitch+3*x+0] = 0;
				pixels[y*pitch+3*x+1] = 0;
				pixels[y*pitch+3*x+2] = 0;
			}
		}
	} else if(mode == DISPLAY_DIVVEL) {
	#pragma omp parallel for
		for(int y = 1; y < f->h-1; y++) {
			for(int x = 1; x < f->w-1; x++) {
				double div = (f->v[2*(y*f->w+x+1)+0]-f->v[2*(y*f->w+x-1)+0]
						+f->v[2*((y+1)*f->w+x)+1]-f->v[2*((y-1)*f->w+x)+1])/2/f->dx;

				div = 255*map(20*fabs(div));
				pixels[y*pitch+3*x+0] = div;
				pixels[y*pitch+3*x+1] = div;
				pixels[y*pitch+3*x+2] = div;

			}

		}
	}

	for(int i = 0; i < f->npart; i++) {
		if(!isnormal(f->parts[i].x) || !isnormal(f->parts[i].y))
			continue;
		int x = f->parts[i].x/f->dx+0.5;
		int y = f->parts[i].y/f->dx+0.5;
		double vx = f->v[2*(y*f->w+x)];
		double vy = f->v[2*(y*f->w+x)+1];
		double vr = sqrt(vx*vx+vy*vy);
		int div = 50*vr;
		
		draw_particle(f->parts[i].x/f->dx,f->parts[i].y/f->dx,div,pitch,pixels);
	}

}

void ui_draw(UI *ui, Fluid *f) {
	uint8_t *pixels;
	int pitch;
	SDL_LockTexture(ui->tex,0,(void **)&pixels,&pitch);
	ui_draw_pixels(f,pixels,pitch,ui->mode);
	pixels[0] = 128;
	pixels[1] = 128;
	pixels[2] = 128;
	pixels[ui->mode] = 255;
	SDL_UnlockTexture(ui->tex);

	SDL_RenderCopy(ui->ren, ui->tex, NULL, NULL);
	SDL_RenderPresent(ui->ren);
}

void ui_poll_event(UI *ui, int *quit) {
	SDL_Event e;
	while (SDL_PollEvent(&e)){
		if (e.type == SDL_QUIT){
			*quit = 1;
		}
		if(e.type == SDL_KEYDOWN) {
			switch(e.key.keysym.scancode) {
				case SDL_SCANCODE_SPACE:
					ui->mode++;
					if(ui->mode >= _DISPLAY_MODES)
						ui->mode = 0;
					break;
				case SDL_SCANCODE_BACKSPACE:
					ui->mode--;
					if(ui->mode < 0)
						ui->mode = _DISPLAY_MODES-1;
					break;
				default:
					break;
			}
		}

	}
}
