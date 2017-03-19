#ifndef UI_H
#define UI_H

#include <SDL.h>
#include "fluid.h"
#include <stdint.h>

// different ways to visualize the velocity field.
enum DisplayMode {
	DISPLAY_NONE= -1,
	DISPLAY_VELDIR = 0,	// velocity direction
	DISPLAY_VEL,		// velocity direction + absolute value
	DISPLAY_DIVVEL,		// velocity divergence
	DISPLAY_PART,		// only particles
	_DISPLAY_MODES = 4
};
typedef enum DisplayMode DisplayMode;

struct UI {
	int w;
	int h;
	int scale;

	DisplayMode mode; // display mode

	SDL_Window *win;
	SDL_Renderer *ren;
	SDL_Texture *tex;
};
typedef struct UI UI;

int ui_init(UI *ui, int w, int h, int scale);
int ui_cleanup(UI *ui);
void ui_draw_pixels(Fluid *f, uint8_t *pixels, int pitch, DisplayMode mode);
void ui_draw(UI *ui, Fluid *f);

void ui_poll_event(UI *ui, int *quit);

#endif
