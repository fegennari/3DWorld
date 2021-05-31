// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/11/03
#pragma once

#include "3DWorld.h"
#include "collision_detect.h"

float const MESH_MIN_Z = -1.0E6; // special mesh height that is guaranteed to be smaller than any mesh zval

extern float sthresh[2][2];


struct ripple_state {
	float rval, acc;
};


class compute_shader_t;
class compute_shader_comp_t;

class mesh_xy_grid_cache_t {

	vector<float> xyterms, sine_mag_terms, cached_vals;
	unsigned cur_nx, cur_ny, yterms_start, tid;
	float mx0, my0, mdx, mdy, sine_offset;
	int gen_mode, gen_shape;
	bool do_glaciate;

	// compute_shader_t or compute_shader_comp_t, but only compute_shader_t works for tiled terrain (size not a multiple of block_size=16)
	typedef compute_shader_t grid_gen_shader_t;
	//typedef compute_shader_comp_t grid_gen_shader_t;
	grid_gen_shader_t *cshader;

	void run_gpu_simplex();
	void cache_gpu_simplex_vals();

public:
	mesh_xy_grid_cache_t() : cur_nx(0), cur_ny(0), yterms_start(0), tid(0), mx0(0.0), my0(0.0), mdx(0.0), mdy(0.0), sine_offset(0.0),
		gen_mode(MGEN_SINE), gen_shape(0), do_glaciate(0), cshader(nullptr) {}
	~mesh_xy_grid_cache_t() {clear_context();}
	bool build_arrays(float x0, float y0, float dx, float dy, unsigned nx, unsigned ny, bool cache_values=0, bool force_sine_mode=0, bool no_wait=0);
	void enable_glaciate();
	float eval_index(unsigned x, unsigned y, int min_start_sin=0, bool use_cache=1) const;
	void clear_context();
	void free_cshader();
};


struct valley { // size = 70

	struct spill_func { // size = 16

		short index, i, j, si, sj, spill;
		float z_over;

		spill_func() : index(-1), i(0), j(0), si(0), sj(0), spill(0), z_over(0.0) {}
		spill_func(short ix, short i_, short j_, short si_, short sj_, short s, float z)
			: index(ix), i(i_), j(j_), si(si_), sj(sj_), spill(s), z_over(z) {}
	};

	short x, y, spill_index;
	bool has_spilled;
	float w_volume, spill_vol, lwv, zval, min_zval, dz, area, fvol, depth, blood_mix, mud_mix, spill_integral;
	spill_func sf;

	valley(short x_=0, short y_=0) : x(x_), y(y_), spill_index(-1), has_spilled(0),
		w_volume(0.0), spill_vol(0.0), lwv(0.0), zval(0.0), min_zval(0.0), dz(0.0),
		area(0.0), fvol(0.0), depth(0.0), blood_mix(0.0), mud_mix(0.0), spill_integral(0.0) {}
	void copy_state_from(valley const &v);
	void create(int wsi);
	float get_volume() const;
};


struct valley_w { // size = 8
	short wsi, x, y, inside8;
	valley_w() : wsi(0), x(0), y(0), inside8(0) {}
};


struct surf_adv { // size = 4
	short x, y;
	surf_adv() : x(0), y(0) {}
	void assign(int x_, int y_) {x = x_; y = y_;}
};


float const hmap_large_zval = 1000.0;

struct hmap_params_t {
	//int mode, shape;
	float plat_bot, plat_h, plat_s, plat_max, crat_h, crat_s, crack_lo, crack_hi, crack_d, sine_mag, sine_freq, sine_bias, volcano_width, volcano_height;
	hmap_params_t() : plat_bot(hmap_large_zval), plat_h(0), plat_s(0), plat_max(0), crat_h(hmap_large_zval), crat_s(0), crack_lo(0), crack_hi(0), crack_d(0),
		sine_mag(0), sine_freq(0), sine_bias(0), volcano_width(0), volcano_height(0) {}
	bool need_postproc() const {return (plat_bot < hmap_large_zval || crat_h < hmap_large_zval || crack_lo < crack_hi);}
};


// should be const, but depend on mesh size
extern int MESH_X_SIZE, MESH_Y_SIZE, MESH_Z_SIZE, MAX_XY_SIZE, XY_MULT_SIZE, XY_SUM_SIZE, I_TIMESCALE;
extern int MESH_SIZE[];
extern float DX_VAL, DY_VAL, DZ_VAL, HALF_DXY, DX_VAL_INV, DY_VAL_INV, dxdy, CLOUD_CEILING, LARGE_ZVAL;
extern float SCENE_SIZE[];
extern float czmin, DZ_VAL2, DZ_VAL_INV2;


// extern global arrays dependent on mesh size
extern valley_w  **watershed_matrix;
extern char      **wminside;
extern vector3d  **wat_surf_normals;
extern vector3d  **wat_vert_normals;
extern float     **mesh_height;
extern float     **z_min_matrix;
extern float     **accumulation_matrix;
extern float     **h_collision_matrix;
extern coll_cell **v_collision_matrix;
extern float     **water_matrix;
extern short     **spillway_matrix;
extern surf_adv  **w_motion_matrix;
extern vector3d  **surface_normals;
extern vector3d  **vertex_normals;
extern float     **charge_dist;
extern float     **surface_damage;
extern ripple_state **ripples;
extern unsigned char **mesh_draw;
extern unsigned char **water_enabled;
extern unsigned char **flower_weight;


inline float get_xval(int xpos)  {return -X_SCENE_SIZE + DX_VAL*xpos;}
inline float get_yval(int ypos)  {return -Y_SCENE_SIZE + DY_VAL*ypos;}
inline float get_zval(int zpos)  {return czmin + DZ_VAL2*zpos;}
inline float get_zval2(int zpos) {return -Z_SCENE_SIZE + DZ_VAL*zpos;}
inline float get_zval_min()      {return get_zval(0);}
inline float get_zval_max()      {return get_zval(max(MESH_SIZE[2], 1));}

inline int get_xpos(float xval) {return int((xval + X_SCENE_SIZE)*DX_VAL_INV + 0.5);}
inline int get_ypos(float yval) {return int((yval + Y_SCENE_SIZE)*DY_VAL_INV + 0.5);}
inline int get_xpos_floor(float xval) {return int((xval + X_SCENE_SIZE)*DX_VAL_INV);}
inline int get_ypos_floor(float yval) {return int((yval + Y_SCENE_SIZE)*DY_VAL_INV);}
inline int get_zpos(float z)    {return int((z - czmin)*DZ_VAL_INV2);}
inline int get_zpos2(float z)   {return int((z + Z_SCENE_SIZE)/DZ_VAL + 0.5);}

inline int get_xpos_clamp(float xval) {return max(0, min(MESH_X_SIZE-1, get_xpos(xval)));}
inline int get_ypos_clamp(float yval) {return max(0, min(MESH_Y_SIZE-1, get_ypos(yval)));}
inline int get_xpos_round_down(float xval) {return int((xval + X_SCENE_SIZE)*DX_VAL_INV);}
inline int get_ypos_round_down(float yval) {return int((yval + Y_SCENE_SIZE)*DY_VAL_INV);}


inline float get_dim_val(int val, unsigned dim) {
	switch (dim) {
		case 0: return get_xval(val);
		case 1: return get_yval(val);
		case 2: return get_zval(val);
		default: assert(0);
	}
	return 0;
}

inline int get_dim_pos(float val, unsigned dim) {
	switch (dim) {
		case 0: return get_xpos(val);
		case 1: return get_ypos(val);
		case 2: return get_zpos(val);
		default: assert(0);
	}
	return 0;
}

inline point get_xyz_pos(int x, int y, int z) {
	return point(get_xval(x), get_yval(y), get_zval(z));
}
inline point get_mesh_xyz_pos(int x, int y) { // Note: not bounds checked
	return point(get_xval(x), get_yval(y), mesh_height[y][x]);
}

