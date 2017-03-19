#include "fluid.h"
#include "fileio.h"
#include "ui.h"
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>

int main(int argc, char **argv) {
	Fluid f;
	UI ui;
	FileIO fio;

	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
	
	fluid_new(&f,214,214,0.05,0.005,10000);

	fileio_init(&fio, "out.mkv",  f.w, f.h, DISPLAY_PART);
	int rc = ui_init(&ui, f.w, f.h, 1);
	if(rc != 0) {
		printf("Error setting up UI.\n");
		return 1;
	}
	
	while(1) {
		fluid_step(&f,0.005);
		fluid_step(&f,0.005);
		fluid_step(&f,0.005);
		//printf("step %g\n", f.t);
		
		ui_draw(&ui,&f);
		fileio_write_frame(&fio, &f);

		int quit = 0;
		ui_poll_event(&ui, &quit);
		if(quit)
			break;
	}

	fileio_cleanup(&fio);
	ui_cleanup(&ui);
	fluid_cleanup(&f);
	return 0;
}
