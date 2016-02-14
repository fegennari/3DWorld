// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 2/13/03

#ifndef _MAIN_H_
#define _MAIN_H_


#include "3DWorld.h"
#include "u_event.h"
#include "tree_3dw.h"


// Global Variables
extern bool underwater;
extern int xoff, yoff, xoff2, yoff2, camera_change, animate, init_x, game_mode, map_mode;
extern int display_framerate, is_cloudy, recreated;
extern int displayed, show_framerate, pause_frame, show_fog, camera_view, camera_mode, camera_surf_collide;
extern int camera_coll_id, display_mode, frame_counter, tree_mode;
extern float temperature, water_plane_z, sm_green_int;
extern float zcenter, leaf_color_coherence, tree_color_coherence, sun_radius, moon_radius, earth_radius, map_zoom;
extern vector3d wind;
extern colorRGBA leaf_base_color;
extern lightning l_strike;
extern obj_group obj_groups[];

// camera variables
extern double c_radius, c_theta, c_phi, up_theta, camera_y;
extern float sun_rot, moon_rot;
extern point cpos2;


#endif

