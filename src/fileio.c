#include "fileio.h"
#include "ui.h"
#include <stdio.h>

int fileio_init(FileIO *f, const char *filename, int w, int h, DisplayMode mode) {
	char *cmd;
	asprintf(&cmd,"ffmpeg -y -s %dx%d -f rawvideo -pix_fmt rgb24 -i pipe: -codec:v libx264 -crf 0 -pix_fmt yuv420p %s",w,h,filename);
	f->w = w;
	f->h = h;
	f->mode = mode;
	f->ffmpeg = popen(cmd,"w");
	f->pixels = calloc(w*h,3);
	free(cmd);

	return 0;
}

void fileio_write_frame(FileIO *fio, Fluid *f) {
	ui_draw_pixels(f, fio->pixels, 3*fio->w, fio->mode);
	fwrite(fio->pixels,fio->w*fio->h*3,1,fio->ffmpeg);
}

void fileio_cleanup(FileIO *f) {
	free(f->pixels);
	pclose(f->ffmpeg);
}
