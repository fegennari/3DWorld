// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/10/02

#include "3DWorld.h"
#include "mesh.h"


int   const DEF_MESH_X_SIZE  = 128;
int   const DEF_MESH_Y_SIZE  = 128;
int   const DEF_MESH_Z_SIZE  = 1;
float const DEF_X_SCENE_SIZE = 4.0;
float const DEF_Y_SCENE_SIZE = 4.0;
float const DEF_Z_SCENE_SIZE = 4.0;
int const INTERPOLATE_METHOD = 1;


// Global Variables
int matrix_alloced(0);
int MESH_X_SIZE(DEF_MESH_X_SIZE), MESH_Y_SIZE(DEF_MESH_Y_SIZE), MESH_Z_SIZE(DEF_MESH_Z_SIZE);
float X_SCENE_SIZE(DEF_X_SCENE_SIZE), Y_SCENE_SIZE(DEF_Y_SCENE_SIZE), Z_SCENE_SIZE(DEF_Z_SCENE_SIZE);
int MESH_SIZE[3], MAX_XY_SIZE, XY_MULT_SIZE, XY_SUM_SIZE, MAX_RUN_DIST, I_TIMESCALE;
float SCENE_SIZE[3], MESH_HEIGHT, OCEAN_DEPTH, XY_SCENE_SIZE, TWO_XSS, TWO_YSS;
float DX_VAL, DY_VAL, HALF_DXY, DX_VAL_INV, DY_VAL_INV, DZ_VAL, dxdy, CLOUD_CEILING, LARGE_ZVAL;


// global arrays dependent on mesh size
valley_w  **watershed_matrix = NULL; // inside: 0 = outside mesh, 1 = inside mesh, 2 = under water level
char      **wminside = NULL;
vector3d  **wat_surf_normals = NULL;
vector3d  **wat_vert_normals = NULL;
float     **mesh_height = NULL;
float     **z_min_matrix = NULL;
float     **accumulation_matrix = NULL;
float     **h_collision_matrix = NULL;
coll_cell **v_collision_matrix = NULL; // *** should use 3X 1D arrays ***
float     **water_matrix = NULL;
short     **spillway_matrix = NULL;
surf_adv  **w_motion_matrix = NULL;
vector3d  **surface_normals = NULL;
vector3d  **vertex_normals = NULL;
float     **charge_dist = NULL;
float     **surface_damage = NULL;
ripple_state **ripples = NULL;
char      **mesh_draw = NULL;
char      **water_enabled = NULL;
short     ***volume_matrix = NULL;
unsigned char ***shadow_mask = NULL;

extern bool last_int, mesh_invalidated;
extern int world_mode, MAX_RUN_DIST, island, xoff, yoff, I_TIMESCALE2, DISABLE_WATER;
extern float zmax, zmin, water_plane_z, def_water_level, temperature, max_obj_radius;
extern point ocean;


void update_motion_zmin_matrices(int xpos, int ypos);


void set_scene_constants() {

	MESH_SIZE[0]  = MESH_X_SIZE;
	MESH_SIZE[1]  = MESH_Y_SIZE;
	MESH_SIZE[2]  = MESH_Z_SIZE;
	SCENE_SIZE[0] = X_SCENE_SIZE;
	SCENE_SIZE[1] = Y_SCENE_SIZE;
	SCENE_SIZE[2] = Z_SCENE_SIZE;
	MAX_XY_SIZE   = max(MESH_X_SIZE, MESH_Y_SIZE);
	XY_MULT_SIZE  = MESH_X_SIZE*MESH_Y_SIZE;
	XY_SUM_SIZE   = MESH_X_SIZE + MESH_Y_SIZE;
	I_TIMESCALE   = min(MAX_I_TIMESCALE, max(1, int(XY_SUM_SIZE/128)));
	I_TIMESCALE2  = I_TIMESCALE;
	MESH_HEIGHT   = 0.10*Z_SCENE_SIZE;
	OCEAN_DEPTH   = 0.05*Z_SCENE_SIZE;
	XY_SCENE_SIZE = 0.5*(X_SCENE_SIZE + Y_SCENE_SIZE);
	TWO_XSS       = 2.0*X_SCENE_SIZE;
	TWO_YSS       = 2.0*Y_SCENE_SIZE;
	DX_VAL        = TWO_XSS/(float)MESH_X_SIZE;
	DY_VAL        = TWO_YSS/(float)MESH_Y_SIZE;
	HALF_DXY      = 0.5*(DX_VAL + DY_VAL);
	DX_VAL_INV    = 1.0/DX_VAL;
	DY_VAL_INV    = 1.0/DY_VAL;
	DZ_VAL        = float(2.0*Z_SCENE_SIZE)/(float)max(MESH_Z_SIZE, 1);
	dxdy          = DX_VAL*DY_VAL;
	MAX_RUN_DIST  = min(MESH_X_SIZE, MESH_Y_SIZE)/2;
	CLOUD_CEILING = CLOUD_CEILING0*Z_SCENE_SIZE;
	LARGE_ZVAL    = 100.0*CLOUD_CEILING;
}


// This file is mostly preprocessing.
void alloc_matrices() { // called at the beginning of main()

	RESET_TIME;
	assert(!matrix_alloced && MESH_Y_SIZE >= 4 && MESH_Y_SIZE <= 4096 && MESH_X_SIZE >= 4 && MESH_X_SIZE <= 4096);
	assert(X_SCENE_SIZE > 0.0 && Y_SCENE_SIZE > 0.0 && Z_SCENE_SIZE > 0.0);
	cout << "mesh = " << MESH_X_SIZE << "x" << MESH_Y_SIZE << ", scene = " << X_SCENE_SIZE << "x" << Y_SCENE_SIZE << endl;

	// reset parameters in case size has changed
	set_scene_constants();
	set_coll_rmax(max_obj_radius); // recompute
	matrix_gen_2d(watershed_matrix);
	matrix_gen_2d(wminside);
	matrix_gen_2d(wat_vert_normals);
	matrix_gen_2d(mesh_height);
	matrix_gen_2d(z_min_matrix);
	matrix_gen_2d(accumulation_matrix);
	matrix_gen_2d(h_collision_matrix);
	matrix_gen_2d(v_collision_matrix);
	matrix_gen_2d(water_matrix);
	matrix_gen_2d(spillway_matrix);
	matrix_gen_2d(w_motion_matrix);
	matrix_gen_2d(surface_normals);
	matrix_gen_2d(vertex_normals);
	matrix_gen_2d(charge_dist);
	matrix_gen_2d(surface_damage);
	matrix_gen_2d(ripples);
	matrix_gen_3d(volume_matrix, MESH_Z_SIZE);
	matrix_gen_3d(shadow_mask,   NUM_LIGHT_SRC);
	matrix_gen_2d(wat_surf_normals, MESH_X_SIZE, 2); // only two rows
	matrix_alloced = 1;
	PRINT_TIME("matrix alloc");
}


void delete_matrices() { // called at the end of main()

	RESET_TIME;
	if (!matrix_alloced) return;
	matrix_delete_2d(watershed_matrix);
	matrix_delete_2d(wminside);
	matrix_delete_2d(wat_surf_normals);
	matrix_delete_2d(wat_vert_normals);
	matrix_delete_2d(mesh_height);
	matrix_delete_2d(z_min_matrix);
	matrix_delete_2d(accumulation_matrix);
	matrix_delete_2d(h_collision_matrix);
	matrix_delete_2d(v_collision_matrix);
	matrix_delete_2d(water_matrix);
	matrix_delete_2d(spillway_matrix);
	matrix_delete_2d(w_motion_matrix);
	matrix_delete_2d(surface_normals);
	matrix_delete_2d(vertex_normals);
	matrix_delete_2d(charge_dist);
	matrix_delete_2d(surface_damage);
	matrix_delete_2d(ripples);
	matrix_delete_3d(volume_matrix, MESH_Z_SIZE);
	matrix_delete_3d(shadow_mask,   NUM_LIGHT_SRC);
	matrix_alloced = 0;
	PRINT_TIME("matrix delete");
}


void compute_matrices() {

	last_int = 0;

	// initialize objects
	reset_other_objects_status();
	matrix_clear_2d(accumulation_matrix);
	matrix_clear_2d(surface_damage);
	matrix_clear_2d(ripples);
	matrix_clear_2d(spillway_matrix);
	
	for (int i = 0; i < NUM_LIGHT_SRC; ++i) {
		matrix_clear_2d(shadow_mask[i]);
	}
	remove_all_coll_obj();

	for (int y = 0; y < MESH_Y_SIZE; ++y) {
		for (int x = 0; x < MESH_X_SIZE; ++x) {
			update_matrix_element(x, y);
		}
	}
	gen_mesh_bsp_tree();
}


void update_matrix_element(int xpos, int ypos) {

	water_matrix    [ypos][xpos] = ((world_mode == WMODE_INF_TERRAIN) ? water_plane_z : def_water_level);
	wat_vert_normals[ypos][xpos] = plus_z;
	charge_dist     [ypos][xpos] = 0.4 + 0.6*rand_float();
	calc_matrix_normal_at(mesh_height, vertex_normals, surface_normals, mesh_draw, MESH_X_SIZE, MESH_Y_SIZE, xpos, ypos);
}


struct mesh_update_t {
	int x, y;
	float old_mh;

	mesh_update_t(int x_, int y_, float mh) : x(x_), y(y_), old_mh(mh) {}
};


// mode is currently: 0=crater, 1=erosion
void update_mesh_height(int xpos, int ypos, int rad, float scale, float offset, int mode, bool is_large_change) {

	//RESET_TIME;
	assert(rad >= 0);
	int const x1(max(0, xpos-rad)), y1(max(0, ypos-rad));
	int const x2(min(MESH_X_SIZE-1, xpos+rad)), y2(min(MESH_Y_SIZE-1, ypos+rad));
	float const zbot(island ? (ocean.z + 0.01) : (zbottom - MESH_LOWEST_DZ));
	vector<mesh_update_t> to_update; // {x, y}
	set<pair<int, int> > grass_update;

	// first pass to update mesh
	for (int i = y1; i <= y2; ++i) {
		for (int j = x1; j <= x2; ++j) {
			float const dh(sqrt(float((i - ypos)*(i - ypos) + (j - xpos)*(j - xpos))));
			if (dh > rad)     continue;
			float const mh(mesh_height[i][j]);
			if (mh < zbottom) continue;
			float delta_h;

			if (mode == 0) { // crater
				delta_h = offset + 0.9*rad - dh; // push the mesh up along the lip of the crater
			}
			else { // erosion
				delta_h = 1.0/(offset + dh);
			}
			float const mh2(max(zbot, (mh - scale*delta_h)));
			mesh_height[i][j] = mh2;//min(mh, mh2);
			if (h_collision_matrix[i][j] == mh) {h_collision_matrix[i][j] = mesh_height[i][j];} // hcm was determined by mh
			to_update.push_back(mesh_update_t(j, i, mh));
			for (int d = 0; d < 4; ++d) {grass_update.insert(make_pair(max(0, j-(d&1)), max(0, i-((d>>1)&1))));}
		}
	}

	// second pass to update adjacency data
	for (vector<mesh_update_t>::const_iterator i = to_update.begin(); i != to_update.end(); ++i) {
		update_matrix_element(i->x, i->y); // requires mesh_height
		update_motion_zmin_matrices(i->x, i->y); // requires mesh_height
	}

	// third pass to update water, which depends on w_motion_matrix
	for (vector<mesh_update_t>::const_iterator i = to_update.begin(); i != to_update.end(); ++i) {
		update_water_zval(i->x, i->y, i->old_mh);
	}
	bool cobjs_updated(update_scenery_zvals(x1, y1, x2, y2));

	if (is_large_change) {
		cobjs_updated |= update_decid_tree_zvals(x1, y1, x2, y2);
		cobjs_updated |= update_small_tree_zvals(x1, y1, x2, y2);
		if (cobjs_updated) {build_cobj_tree(0, 0);} // slow, but probably necessary
		//calc_visibility(SUN_SHADOW | MOON_SHADOW); // too slow
	}

	// fourth pass to update grass, after cobjs/shadows have been updated
	for (set<pair<int, int> >::const_iterator i = grass_update.begin(); i != grass_update.end(); ++i) {
		grass_mesh_height_change(i->first, i->second, (is_large_change && 0));
	}
	// update waypoints?
	mesh_invalidated = 1;
	//PRINT_TIME("Mesh Height Update");
}


// not always correct if !enabled[i][j]
vector3d get_matrix_surf_norm(float **matrix, char **enabled, int xsize, int ysize, int x, int y) {

	assert(matrix);
	float nx(0.0), ny(0.0);
	float const mhij(matrix[y][x]);
	bool const test_md(enabled != NULL);
	assert(y >= 0 && x >= 0 && y < ysize && x < xsize);

	if (!test_md || enabled[y][x]) {
		if (y < ysize-1) {
			if (!test_md || enabled[y+1][x]) {ny =  DX_VAL*(mhij - matrix[y+1][x]);}
		}
		else {
			if (!test_md || enabled[y-1][x]) {ny = -DX_VAL*(mhij - matrix[y-1][x]);}
		}
		if (x < xsize-1) {
			if (!test_md || enabled[y][x+1]) {nx =  DY_VAL*(mhij - matrix[y][x+1]);}
		}
		else {
			if (!test_md || enabled[y][x-1]) {nx = -DY_VAL*(mhij - matrix[y][x-1]);}
		}
	}
	return vector3d(nx, ny, dxdy).get_norm();
}


void calc_matrix_normal_at(float **matrix, vector3d **vn, vector3d **sn, char **enabled, int xsize, int ysize, int xpos, int ypos) {

	vector3d const norm(get_matrix_surf_norm(matrix, enabled, xsize, ysize, xpos, ypos));
	sn[ypos][xpos] = norm;
	vn[ypos][xpos] = (norm + sn[max(ypos-1, 0)][xpos] + sn[max(ypos-1, 0)][max(xpos-1, 0)] + sn[ypos][max(xpos-1, 0)])*0.25;
}


void calc_matrix_normals(float **matrix, vector3d **vn, vector3d **sn, char **enabled, int xsize, int ysize) {

	assert(matrix && vn && sn);

	for (int y = 0; y < ysize; ++y) {
		for (int x = 0; x < xsize; ++x) {
			calc_matrix_normal_at(matrix, vn, sn, enabled, xsize, ysize, x, y);
		}
	}
}


void get_matrix_point(int xpos, int ypos, point &pt) {

	if (mesh_height != NULL && !point_outside_mesh(xpos, ypos)) {
		pt.assign(get_xval(xpos), get_yval(ypos), mesh_height[ypos][xpos]);
	}
}


int is_in_ice(int xpos, int ypos) {

	if (DISABLE_WATER || temperature > W_FREEZE_POINT || point_outside_mesh(xpos, ypos)) return 0;
	return (wminside[ypos][xpos] && (!island || mesh_height[ypos][xpos] > (ocean.z + SMALL_NUMBER)));
}


float get_exact_zval(float xval, float yval) {

	return eval_one_surface_point(((xval + X_SCENE_SIZE)*DX_VAL_INV + 0.5), ((yval + Y_SCENE_SIZE)*DY_VAL_INV + 0.5));
}


// This really solves for the z value of an x/y point on the mesh plane, but it is somewhat like interpolation.
float interpolate_mesh_zval(float xval, float yval, float rad, int use_real_equation, int ignore_ice) {

	if (world_mode == WMODE_INF_TERRAIN) {return get_exact_zval(xval, yval);}
	int const xpos(get_xpos(xval)), ypos(get_ypos(yval));
	if (point_outside_mesh(xpos, ypos)) {return (use_real_equation ? get_exact_zval(xval, yval) : (zbottom - MESH_LOWEST_DZ));}
	assert(mesh_height != NULL);
	float const xp((xval + X_SCENE_SIZE)*DX_VAL_INV), yp((yval + Y_SCENE_SIZE)*DY_VAL_INV);
	int const x0((int)xp), y0((int)yp);
	bool const xy0_bad(x0 < 0 || y0 < 0 || x0 >= MESH_X_SIZE-1 || y0 >= MESH_Y_SIZE-1);
	float zval;

	if (INTERPOLATE_METHOD == 0 || xy0_bad) {
		vector3d const &norm(surface_normals[ypos][xpos]); // 4 points doesn't define a plane - need two faces/planes
		assert(fabs(norm.z) > TOLERANCE); // can't have a vertical mesh quad
		zval = (-norm.x*(xval - get_xval(xpos)) - norm.y*(yval - get_yval(ypos)) + norm.z*mesh_height[ypos][xpos])/norm.z;
	}
	else {
		float const xpi(fabs(xp - (float)x0)), ypi(fabs(yp - (float)y0));
		zval = (1.0 - xpi)*((1.0 - ypi)*mesh_height[y0][x0] + ypi*mesh_height[y0+1][x0]) + xpi*((1.0 - ypi)*mesh_height[y0][x0+1] + ypi*mesh_height[y0+1][x0+1]);
	}
	if (rad > 0.0 && !xy0_bad) {
		float hcm(min(h_collision_matrix[y0][x0], h_collision_matrix[y0+1][x0+1]));
		hcm = min(hcm, min(h_collision_matrix[y0][x0+1], h_collision_matrix[y0+1][x0]));
		if (zval + 0.5*rad + SMALL_NUMBER < hcm) return h_collision_matrix[ypos][xpos];
	}
	if (!ignore_ice && is_in_ice(xpos, ypos)) zval = max(zval, water_matrix[ypos][xpos]); // on ice
	return zval;
}


float int_mesh_zval_pt_off(point const &pos, int use_real_equation, int ignore_ice) {
	
	return interpolate_mesh_zval((pos.x-DX_VAL*xoff), (pos.y-DY_VAL*yoff), 0.0, use_real_equation, ignore_ice);
}


void update_motion_zmin_matrices(int xpos, int ypos) {

	float const old_z(mesh_height[ypos][xpos]);
	float new_z(old_z);
	int new_x(xpos), new_y(ypos), xlevel(xpos), ylevel(ypos);

	// this part used for calc_rest_pos() and draw_spillover() (water)
	for (int dy = -1; dy <= 1; ++dy) {
		for (int dx = -1; dx <= 1; ++dx) {
			if (dx == 0 && dy == 0)         continue;
			int const nx(xpos+dx), ny(ypos+dy);
			if (point_outside_mesh(nx, ny)) continue;
			float const this_z(mesh_height[ny][nx]);
			
			if (this_z < new_z) {
				new_z = this_z;
				new_x = nx;
				new_y = ny;
			}
			else if (dx == 1 && dy == 1 && new_z == old_z && this_z == new_z) { // NE only
				xlevel = nx;
				ylevel = ny;
			}
		}
	}
	if (new_z < old_z) {
		w_motion_matrix[ypos][xpos].assign(new_x, new_y); // was new_x
	}
	else { // force flat areas to still have flow: default flow is NE
		w_motion_matrix[ypos][xpos].assign(xlevel, ylevel);
	}
	float z_min(zmin);

	// this part used to calculate water height/surface intersection in draw_water()
	if (point_interior_to_mesh(xpos, ypos)) {
		z_min = zmax;

		for (int y = -1; y <= 1; ++y) {
			for (int x = ((y == -1) ? 0 : -1); x <= 1; ++x) {
				z_min = min(z_min, mesh_height[ypos+y][xpos+x]); // (-1,-1) intentionally skipped
			}
		}
	}
	z_min_matrix[ypos][xpos] = z_min;
}


void calc_motion_direction() {

	for (int ypos = 0; ypos < MESH_Y_SIZE; ++ypos) {
		for (int xpos = 0; xpos < MESH_X_SIZE; ++xpos) {
			update_motion_zmin_matrices(xpos, ypos);
		}
	}
}


typedef float (*comp_func)(float, float);
inline float comp_min(float A, float B) {return min(A, B);}
inline float comp_max(float A, float B) {return max(A, B);}


float proc_mesh_point(point const &pt, float radius, comp_func cf) {

	int const xpos(get_xpos(pt.x)), ypos(get_ypos(pt.y));
	//if (point_outside_mesh(xpos, ypos)) return zmin;
	float zval(interpolate_mesh_zval(pt.x, pt.y, 0.0, 1, 0));
	if (point_outside_mesh(xpos, ypos)) return zval;
	int const dx((int)ceil(radius*DX_VAL_INV)), dy((int)ceil(radius*DY_VAL_INV));
	int const xv[2] = {max(0, (xpos - dx)), min(MESH_X_SIZE-1, (xpos + dx))};
	int const yv[2] = {max(0, (ypos - dy)), min(MESH_Y_SIZE-1, (ypos + dy))};

	for (unsigned y = 0; y < 2; ++y) {
		for (unsigned x = 0; x < 2; ++x) {
			if (xv[x] != xpos && yv[y] != ypos && !point_outside_mesh(xv[x], yv[y])) {
				zval = cf(zval, mesh_height[yv[y]][xv[x]]);
			}
		}
	}
	return zval;
}


float lowest_mesh_point(point const &pt, float radius) {
	return proc_mesh_point(pt, radius, comp_min);
}

float highest_mesh_point(point const &pt, float radius) {
	return proc_mesh_point(pt, radius, comp_max);
}







