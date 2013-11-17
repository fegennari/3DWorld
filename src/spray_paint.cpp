// 3D World - OpenGL CS184 Computer Graphics Project - Spray Paint/Graffiti Code
// by Frank Gennari
// 9/14/13

#include "collision_detect.h"
#include "gameplay.h"
#include "tree_3dw.h"
#include "openal_wrap.h"

using namespace std;


colorRGBA const NEAR_BLACK(0.01, 0.01, 0.01, 1.0); // not quite black, so that it uses the regular lighting shader, and draw order is correct
unsigned const NUM_PAINT_COLORS = 10;
string const paint_color_names[NUM_PAINT_COLORS] = {"WHITE", "RED", "GREEN", "BLUE", "YELLOW", "PINK", "ORANGE", "PURPLE", "BROWN", "BLACK"};
colorRGBA const paint_colors  [NUM_PAINT_COLORS] = {WHITE, RED, GREEN, BLUE, YELLOW, PINK, ORANGE, PURPLE, BROWN, NEAR_BLACK};

bool spraypaint_mode(0);
unsigned paint_color_ix(0);

extern float CAMERA_RADIUS;
extern coll_obj_group coll_objects;


colorRGBA get_cur_paint_color() {

	assert(paint_color_ix < NUM_PAINT_COLORS);
	return paint_colors[paint_color_ix];
}


void show_cur_spraypaint_mode() {

	string const str(paint_color_names[paint_color_ix] + " Spray Paint");
	print_text_onscreen(str, get_cur_paint_color(), 1.0, TICKS_PER_SECOND, 1); // 1 second
}


void toggle_spraypaint_mode() {

	if (world_mode != WMODE_GROUND) return;
	spraypaint_mode ^= 1;
	show_cur_spraypaint_mode();
}


void change_spraypaint_color(int val) {

	if (world_mode != WMODE_GROUND) return;
	paint_color_ix = (paint_color_ix + NUM_PAINT_COLORS + val) % NUM_PAINT_COLORS;
	show_cur_spraypaint_mode();
}


void draw_spraypaint_crosshair() {

	if (world_mode != WMODE_GROUND) return;
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	enable_blend();
	colorRGBA(get_cur_paint_color(), 0.5).do_glColor();
	draw_circle_normal(0.0009, 0.0010, 64, 1, -0.05);
	disable_blend();
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
}


float get_spray_radius(point const &pos, float &alpha) {

	float const dist(distance_to_camera(pos)), radius(min(0.1, max(0.001, 0.05*dist)));
	if (radius > 0.05) {alpha = 1.0 - 10.0*(radius - 0.05);} // 0.5 - 1.0
	return radius;
}


void spray_paint(bool mode) {

	// spray paint affects flat cobjs (cubes and polygons), mesh, grass, and tree leaves - also adds/removes volume to voxels
	// play sound?
	point const pos(get_camera_pos());
	colorRGBA color(get_cur_paint_color());
	int xpos(0), ypos(0), cindex(-1);
	point coll_pos;
	vector3d coll_norm(plus_z);
	float range(FAR_CLIP);
	bool const mesh_int(get_range_to_mesh(pos, cview_dir, coll_pos, xpos, ypos) == 1);

	if (mesh_int) { // mesh (not ice) intersection
		range    = p2p_dist(pos, coll_pos);
		coll_pos = pos + cview_dir*(range - 0.01); // simple and inexact, but seems OK
		coll_pos.z += SMALL_NUMBER;
	}
	if (check_coll_line_exact(pos, (pos + cview_dir*range), coll_pos, coll_norm, cindex, 0.0, -1, 0, 0, 1, 0)) { // hit cobjs (skip_dynamic=1), ignore voxels
		assert(cindex >= 0 && unsigned(cindex) < coll_objects.size());
		coll_obj &cobj(coll_objects[cindex]);
		float const radius(get_spray_radius(coll_pos, color.alpha));

		if (cobj.status == COLL_STATIC && (!cobj.no_draw() || (cobj.cp.cobj_type != COBJ_TYPE_STD))) { // similar to cobj.can_be_scorched()
			if (decal_contained_in_cobj(cobj, coll_pos, coll_norm, radius, get_max_dim(coll_norm))) {
				int const lifetime(((mode & 1) ? 3600 : 60)*TICKS_PER_SECOND); // 1 min / 1 hour
				gen_decal(coll_pos, radius, coll_norm, BLUR_CENT_TEX, cindex, color, 0, 0, lifetime);
			}
		}
		else if (cobj.is_billboard) { // Note: currently, only tree leaves use billboard cobj
			spraypaint_tree_leaves(coll_pos, 1.5*radius, cindex, color);
		}
	}
	else if (mesh_int) { // mesh intersection
		float const radius(get_spray_radius(coll_pos, color.alpha));
		add_color_to_landscape_texture(color, coll_pos.x, coll_pos.y, 1.5*radius);
		modify_grass_at(coll_pos, 1.5*radius, 0, 0, 0, 1, 1, 0, color);
	}
	gen_sound(SOUND_SPRAY, (pos + CAMERA_RADIUS*cview_dir), 0.2, 1.0);
}

