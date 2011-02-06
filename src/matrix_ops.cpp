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
int const USE_REAL_EQ_WM0    = 0;
int const INTERPOLATE_METHOD = 1;


// Global Variables
int matrix_alloced(0);
int MESH_X_SIZE(DEF_MESH_X_SIZE), MESH_Y_SIZE(DEF_MESH_Y_SIZE), MESH_Z_SIZE(DEF_MESH_Z_SIZE);
float X_SCENE_SIZE(DEF_X_SCENE_SIZE), Y_SCENE_SIZE(DEF_Y_SCENE_SIZE), Z_SCENE_SIZE(DEF_Z_SCENE_SIZE);
int MESH_SIZE[3], MAX_XY_SIZE, XY_MULT_SIZE, XY_SUM_SIZE, MAX_RUN_DIST, I_TIMESCALE;
float SCENE_SIZE[3], MESH_HEIGHT, OCEAN_DEPTH, XY_SCENE_SIZE, TWO_XSS, TWO_YSS;
float DX_VAL, DY_VAL, HALF_DXY, DX_VAL_INV, DY_VAL_INV, DZ_VAL, dxdy, CLOUD_CEILING, LARGE_ZVAL;
float FOG_DIST1, FOG_DIST_W, FOG_DIST_UW0, FOG_DIST_UW3;


// global arrays dependent on mesh size
unsigned char **num_obj_on_mesh;
valley_w  **watershed_matrix; // inside: 0 = outside mesh, 1 = inside mesh, 2 = under water level
char      **wminside;
vector3d  **wat_surf_normals;
vector3d  **wat_vert_normals;
float     **mesh_height = NULL;
float     **z_min_matrix;
float     **accumulation_matrix;
float     **h_collision_matrix;
coll_cell **v_collision_matrix; // *** should use 3X 1D arrays ***
float     **water_matrix;
short     **spillway_matrix;
surf_adv  **w_motion_matrix;
vector3d  **surface_normals;
vector3d  **vertex_normals;
float     **charge_dist;
float     **surface_damage;
ripple_state **ripples;
char      **mesh_draw;
char      **water_enabled;
short     ***volume_matrix;
unsigned char ***shadow_mask;

extern bool last_int, mesh_invalidated;
extern int world_mode, MAX_RUN_DIST, island, xoff, yoff, I_TIMESCALE2, DISABLE_WATER;
extern float zmax, zmin, czmin, czmax, zbottom, water_plane_z, def_water_level, temperature, max_obj_radius;
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
	MESH_HEIGHT   = MESH_HEIGHT0*Z_SCENE_SIZE;
	OCEAN_DEPTH   = OCEAN_DEPTH0*Z_SCENE_SIZE;
	XY_SCENE_SIZE = 0.5*(X_SCENE_SIZE + Y_SCENE_SIZE);
	TWO_XSS       = 2.0*X_SCENE_SIZE;
	TWO_YSS       = 2.0*Y_SCENE_SIZE;
	DX_VAL        = TWO_XSS/(float)MESH_X_SIZE;
	DY_VAL        = TWO_YSS/(float)MESH_Y_SIZE;
	HALF_DXY      = 0.5*(DX_VAL + DY_VAL);
	DX_VAL_INV    = 1.0/DX_VAL;
	DY_VAL_INV    = 1.0/DY_VAL;
	DZ_VAL        = float(2.0*Z_SCENE_SIZE)/(float)MESH_Z_SIZE;
	dxdy          = DX_VAL*DY_VAL;
	MAX_RUN_DIST  = min(MESH_X_SIZE, MESH_Y_SIZE)/2;
	CLOUD_CEILING = CLOUD_CEILING0*Z_SCENE_SIZE;
	LARGE_ZVAL    = 100.0*CLOUD_CEILING;
	FOG_DIST1     = 2.5*Z_SCENE_SIZE;
	FOG_DIST_W    = 2.5*Z_SCENE_SIZE;
	FOG_DIST_UW0  = 1.5*Z_SCENE_SIZE;
	FOG_DIST_UW3  = 1.5*Z_SCENE_SIZE;
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
	init_draw_stats();

	matrix_gen_2d(num_obj_on_mesh);
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
	assert(matrix_alloced);
	matrix_delete_2d(num_obj_on_mesh);
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
	czmin    =  FAR_CLIP;
	czmax    = -FAR_CLIP;

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

	for (int i = 0; i < MESH_Y_SIZE; ++i) {
		for (int j = 0; j < MESH_X_SIZE; ++j) {
			update_matrix_element(i, j);
		}
	}
	gen_mesh_bsp_tree();
}


void update_matrix_element(int i, int j) {

	float const mh(mesh_height[i][j]);
	water_matrix[i][j]     = ((world_mode == WMODE_INF_TERRAIN) ? water_plane_z : def_water_level);
	wat_vert_normals[i][j] = plus_z;
	charge_dist[i][j]      = 0.4 + 0.6*rand_float();
	calc_matrix_normal_at(mesh_height, vertex_normals, surface_normals, mesh_draw, MESH_X_SIZE, MESH_Y_SIZE, i, j);
}


void update_mesh_height(int xpos, int ypos, int rad, float scale, float offset, int mode) {

	assert(rad >= 0);
	int const x1(max(0, xpos-rad)), y1(max(0, ypos-rad));
	int const x2(min(MESH_X_SIZE-1, xpos+rad)), y2(min(MESH_Y_SIZE-1, ypos+rad));
	float const zbot(island ? (ocean.z + 0.01) : (zbottom - 0.04));

	for (int i = y1; i <= y2; ++i) {
		for (int j = x1; j <= x2; ++j) {
			float const dh(sqrt(float((i - ypos)*(i - ypos) + (j - xpos)*(j - xpos))));
			if (dh > rad)     continue;
			float const mh(mesh_height[i][j]);
			if (mh < zbottom) continue;
			float const mh2(max(zbot, (mh - scale*((mode == 0) ? (offset + rad - dh) : 1.0f/(offset + dh)))));
			mesh_height[i][j] = min(mh, mh2);
			if (h_collision_matrix[i][j] == mh) h_collision_matrix[i][j] = mesh_height[i][j]; // hcm was determined by mh
		}
	}
	for (int i = y1; i <= y2; ++i) {
		for (int j = x1; j <= x2; ++j) {
			if ((sqrt(float((i - ypos)*(i - ypos) + (j - xpos)*(j - xpos)))) > rad) continue;
			update_matrix_element(i, j); // requires mesh_height
			update_motion_zmin_matrices(j, i); // requires mesh_height
		}
	}
	update_scenery_zvals(x1, y1, x2, y2);
	update_water_zvals(x1, y1, x2, y2);
	mesh_invalidated = 1;
}


// not always correct if !enabled[i][j]
vector3d get_matrix_surf_norm(float **matrix, char **enabled, int xsize, int ysize, int i, int j) {

	assert(matrix);
	float nx(0.0), ny(0.0);
	float const mhij(matrix[i][j]);
	bool const test_md(enabled != NULL);
	assert(i >= 0 && j >= 0 && i < ysize && j < xsize);

	if (!test_md || enabled[i][j]) {
		if (i < ysize-1) {
			if (!test_md || enabled[i+1][j]) ny =  DX_VAL*(mhij - matrix[i+1][j]);
		}
		else {
			if (!test_md || enabled[i-1][j]) ny = -DX_VAL*(mhij - matrix[i-1][j]);
		}
		if (j < xsize-1) {
			if (!test_md || enabled[i][j+1]) nx =  DY_VAL*(mhij - matrix[i][j+1]);
		}
		else {
			if (!test_md || enabled[i][j-1]) nx = -DY_VAL*(mhij - matrix[i][j-1]);
		}
	}
	return vector3d(nx, ny, dxdy).get_norm();
}


void calc_matrix_normal_at(float **matrix, vector3d **vn, vector3d **sn, char **enabled, int xsize, int ysize, int i, int j) {

	vector3d norm(get_matrix_surf_norm(matrix, enabled, xsize, ysize, i, j));
	sn[i][j] = norm;
	vn[i][j] = (norm + sn[max(i-1, 0)][j] + sn[max(i-1, 0)][max(j-1, 0)] + sn[i][max(j-1, 0)])*0.25;
}


void calc_matrix_normals(float **matrix, vector3d **vn, vector3d **sn, char **enabled, int xsize, int ysize) {

	assert(matrix && vn && sn);

	for (int y = 0; y < ysize; ++y) {
		for (int x = 0; x < xsize; ++x) {
			calc_matrix_normal_at(matrix, vn, sn, enabled, xsize, ysize, y, x);
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

	int const xpos(get_xpos(xval)), ypos(get_ypos(yval));
	
	if ((USE_REAL_EQ_WM0 || world_mode == WMODE_INF_TERRAIN || point_outside_mesh(xpos, ypos)) && use_real_equation) {
		return get_exact_zval(xval, yval);
	}
	if (point_outside_mesh(xpos, ypos)) return (zbottom - SMALL_NUMBER);
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







