#include "video.h"


int video_init(Video *v, char *filename, int w, int h, DisplayMode mode) {
	v->w = w;
	v->h = h;
	v->mode = mode;
	v->pixels = malloc(3*w*h);

	const char *ffmpeg_cmd_fmt = "ffmpeg -s %dx%d -f rawvideo -pix_fmt rgb24 -i pipe: -codec:v libvpx -crf 10 -y -pix_fmt yuv420p %s";


	v->ffmpeg = popen(ffmpegcmd,"w");
}

int video_write_frame(Video *v, Fluid *v) {
}

void video_close(Video *v) {

}

