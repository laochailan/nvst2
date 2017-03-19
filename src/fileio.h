#ifndef FILEIO_H
#define FILEIO_H

#include "fluid.h"
#include "ui.h"
#include <stdint.h>

typedef struct {
	FILE *ffmpeg;

	int w, h;
	DisplayMode mode;
	uint8_t *pixels;
} FileIO;

int fileio_init(FileIO *f, const char *filename, int w, int h, DisplayMode mode);
void fileio_write_frame(FileIO *fio, Fluid *f);
void fileio_cleanup(FileIO *f);

#endif
