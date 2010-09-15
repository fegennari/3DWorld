// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 5/10/02

#include "3DWorld.h"
#include "tree_3dw.h"


#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN 1
#include <windows.h>
#include "pthread.h"

#else

#include <pthread.h>

//#define Sleep sleep
//#define Sleep(x) sleep(unsigned(0.001*x))
#define Sleep(x) usleep(1000*x)

#endif


int calcs_running(0), kill_vis_calc(0), disable_threads(0);


struct tdata {

	char light;
	pthread_t td;
};


extern bool double_buffer_shadow;
extern int visibility_thread, world_mode;



void wait_for_vis_calc_finish() {

	while (calcs_running > 0) {
		kill_vis_calc = 1;
		Sleep(10); // should replace this later
	}
	kill_vis_calc = 0;
}


void *run_visibility(void *data) {

	++calcs_running;
	calc_visibility2(((tdata *)data)->light);
	double_buffer_shadow = 0;
	--calcs_running;
	//pthread_join((tdata *)data)->td, NULL);
	return NULL;
}


void calc_visibility(char light_sources) {

	if (world_mode == WMODE_UNIVERSE) return;
	if (!(light_sources & TREE_ONLY)) clear_all_lightmaps(1);

	if (disable_threads) {
		calc_visibility2(light_sources);
		return;
	}
	pthread_attr_t attr;
	static tdata data;

	if (!visibility_thread) {
		calc_visibility2(light_sources);
		double_buffer_shadow = 0;
		return;
	}
	wait_for_vis_calc_finish();
	data.light = light_sources;
	pthread_attr_init(&attr);
	pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	data.td = pthread_self();
	int err = pthread_create(&data.td, &attr, run_visibility, &data);

	if (err) {
		const char* estr;
		fprintf (stderr, "Error creating shadowing thread.\n");
		estr = strerror(err);
		if (estr) fprintf (stderr, "%s\n", estr);
		else fprintf (stderr, "%d\n", err);
		exit(1);
	}
}


void post_delete_trees(vector<tree> &t_trees) {

	wait_for_vis_calc_finish();
	delete_trees(t_trees);
}


// *************************************************************



