#ifndef VIDEO_H
#define VIDEO_H

#include <stdio.h>
#include "fluid.h"
#include "ui.h"

// This module implements recording to a video file on the disk via ffmpeg.

struct Video {
	FILE *ffmpeg;

	DisplayMode mode;

	int w;
	int h;
	uint8_t *pixels;
};
typedef struct Video Video;

int video_init(Video *v, char *filename, int w, int h, DisplayMode mode);
int video_write_frame(Video *v, Fluid *v);
void video_close(Video *v);

#endif
