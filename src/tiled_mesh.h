// 3D World - tiled terrain classes
// by Frank Gennari
// 4/4/13

#ifndef _TILED_TERRAIN_H_
#define _TILED_TERRAIN_H_

#include "3DWorld.h"
#include "function_registry.h"
#include "inlines.h"
#include "mesh.h"
#include "small_tree.h"
#include "scenery.h"
#include "grass.h"
#include "tree_3dw.h"
#include "shadow_map.h"


bool const ENABLE_TREE_LOD    = 1; // faster but has popping artifacts
bool const ENABLE_TERRAIN_ENV = 1;
bool const GRASS_CLOUD_SHADOWS= 1; // slow, but looks nice
bool const USE_TREE_BILLBOARDS= 1; // decidious trees: faster but lower quality
int  const TILE_RADIUS        = 6; // in mesh sizes
unsigned const NUM_LODS       = 5; // > 0
float const TREE_LOD_THRESH   = 6.0;
float const GEOMORPH_THRESH   = 6.0;
float const SCENERY_THRESH_REF= 2.0;
float const SCENERY_THRESH    = 5.0;
float const GRASS_THRESH      = 1.6;
float const SMAP_NEW_THRESH   = 1.2;
float const SMAP_DEL_THRESH   = 1.3;
float const BCUBE_ZTOLER      = 1.0E-6;


extern int xoff, yoff, frame_counter;
extern float grass_length, water_plane_z;



class lightning_strike_t {

	int time;
	line3d path;

public:
	lightning_strike_t() : time(0) {}
	bool enabled() const {return (time > 0);}
	void clear() {path.points.clear(); time = 0;}
	point get_pos() const;
	void gen();
	void update();
	void draw() const;
	void end_draw() const;
};


struct ix_sz_pair {
	unsigned short ix, sz;
	ix_sz_pair() : ix(0), sz(0) {}
};

class crack_ibuf_t {
	vector<ix_sz_pair> offsets;

public:
	void gen_offsets(vector<unsigned> &indices, unsigned size);
	unsigned get_index(unsigned dim, unsigned dir, unsigned cur_lod, unsigned adj_lod) const;
	ix_sz_pair const &lookup(unsigned ix) const {assert(ix < offsets.size()); return offsets[ix];}
};


class tile_t;

struct tile_smap_data_t : public smap_data_t {

	int dxoff, dyoff;
	tile_t *tile;

	tile_smap_data_t(unsigned tu_id_, unsigned smap_sz_, tile_t *tile_, smap_data_state_t const &init_state=smap_data_state_t())
		: smap_data_t(tu_id_, smap_sz_, init_state), dxoff(0), dyoff(0), tile(tile_) {}
	virtual void render_scene_shadow_pass(point const &lpos);
	virtual bool needs_update(point const &lpos);
};


class tile_shadow_map_manager {

	vector<smap_data_state_t> free_list[NUM_LIGHT_SRC];
public:
	tile_smap_data_t new_smap_data(unsigned tu_id, tile_t *tile, unsigned light);
	void release_smap_data(tile_smap_data_t &smd, unsigned light);
	void clear_context();
};


inline float get_tile_width        () {return (X_SCENE_SIZE + Y_SCENE_SIZE);}
inline float get_scaled_tile_radius() {return TILE_RADIUS*get_tile_width();}
inline float get_tree_scale_denom  () {return max(1.0f, TREE_LOD_THRESH*calc_tree_size());}


struct tile_xy_pair {

	int x, y;
	tile_xy_pair(int x_=0, int y_=0) : x(x_), y(y_) {}
	bool operator<(tile_xy_pair const &t) const {return ((y == t.y) ? (x < t.x) : (y < t.y));}
	void operator+=(tile_xy_pair const &tp) {x += tp.x; y += tp.y;}
	void operator-=(tile_xy_pair const &tp) {x -= tp.x; y -= tp.y;}
	tile_xy_pair operator+(tile_xy_pair const &tp) const {return tile_xy_pair(x+tp.x, y+tp.y);}
	tile_xy_pair operator-(tile_xy_pair const &tp) const {return tile_xy_pair(x-tp.x, y-tp.y);}
};

tile_t *get_tile_from_xy(tile_xy_pair const &tp);


struct tile_cloud_t : public volume_part_cloud {

	point pos;
	vector3d size; // {x, y, z}

	float get_rmax() const {return size.get_max_val();}
	void draw(vpc_shader_t &s, vector3d const &xlate) const;
};

class tile_cloud_manager_t : public vector<tile_cloud_t> {

	bool generated;
	cube_t bcube;
	vector<pair<float, unsigned>> sorted;
public:
	tile_cloud_manager_t() : generated(0) {}
	void gen(int x1, int y1, int x2, int y2);
	bool any_visible(vector3d const &xlate) const;
	void draw(vpc_shader_t &s, vector3d const &xlate);
};


class tile_t {

public:
	struct offset_t {
		int dxoff, dyoff;

		offset_t(int dxoff_=0, int dyoff_=0) : dxoff(dxoff_), dyoff(dyoff_) {}
		void set_from_xyoff2() {dxoff = -xoff2; dyoff = -yoff2;}
		vector3d get_xlate() const {return vector3d(((xoff - xoff2) - dxoff)*DX_VAL, ((yoff - yoff2) - dyoff)*DY_VAL, 0.0);}
		vector3d subtract_from(offset_t const &o) const {return vector3d((o.dxoff - dxoff)*DX_VAL, (o.dyoff - dyoff)*DY_VAL, 0.0);}
	};
	struct tree_map_val {
		unsigned char ao, sh;
		tree_map_val() : ao(255), sh(255) {}
	};

private:
	int x1, y1, x2, y2, wx1, wy1, wx2, wy2, last_occluded_frame;
	unsigned weight_tid, height_tid, normal_tid, shadow_tid;
	unsigned size, stride, zvsize, base_tsize, gen_tsize;
	float radius, mzmin, mzmax, ptzmax, dtzmax, trmax, xstart, ystart, min_normal_z, deltax, deltay;
	bool shadows_invalid, recalc_tree_grass_weights, mesh_height_invalid, in_queue, last_occluded, has_any_grass, is_distant, no_trees;
	offset_t mesh_off, ptree_off, dtree_off, scenery_off;
	float sub_zmin[4][4], sub_zmax[4][4];
	vector<float> zvals;
	vector<tree_map_val> tree_map;
	vector<unsigned char> mesh_weight_data, weight_data, ao_lighting;
	vector<unsigned char> smask[NUM_LIGHT_SRC];
	vector<float> sh_out[NUM_LIGHT_SRC][2];
	vect_smap_t<tile_smap_data_t> smap_data;
	small_tree_group pine_trees;
	scenery_group scenery;
	tree_cont_t decid_trees;
	flower_tile_manager_t flowers;
	tile_cloud_manager_t clouds;

	struct grass_block_t {
		unsigned ix; // 0 is unused
		float zmin, zmax;
		grass_block_t() : ix(0), zmin(0.0), zmax(0.0) {}
	};

	vector<grass_block_t> grass_blocks;

	struct terrain_params_t {
		float hoff, hscale, veg, grass, dirt;
		terrain_params_t() : hoff(0.0), hscale(1.0), veg(1.0), grass(1.0), dirt(0.0) {}
	};

	terrain_params_t params[2][2]; // {ylo,yhi} x {xlo,xhi}

	void update_terrain_params();
	unsigned get_lod_level(bool reflection_pass) const;

public:
	tile_t();
	tile_t(unsigned size_, int x, int y);
	// can't free in the destructor because the gl context may be destroyed before this point
	//~tile_t() {clear_vbo_tid();}
	float calc_radius() const {return 0.5*sqrt(deltax*deltax + deltay*deltay)*size;} // approximate (lower bound)
	float get_zmin() const {return mzmin;}
	float get_zmax() const {return mzmax;}
	float get_tile_zmax() const {return max((mzmax + (grass_blocks.empty() ? 0.0f : grass_length)), max(ptzmax, dtzmax));}
	float get_zval(int x, int y) const {assert(!zvals.empty()); assert(x >= 0 && y >= 0 && x < (int)zvsize && y < (int)zvsize); return zvals[y*zvsize + x];}
	bool has_water() const {return (mzmin < water_plane_z);}
	bool all_water() const {return (mzmax < water_plane_z);} // get_tile_zmax()? - grass and trees should not be underwater
	bool can_have_trees() const {return (!no_trees && !is_distant && !all_water());}
	bool pine_trees_generated() const {return pine_trees.generated;}
	bool has_pine_trees() const {return (pine_trees_generated() && !pine_trees.empty());}
	void invalidate_mesh_height() {mesh_height_invalid = 1;}
	float get_avg_veg() const {return 0.25*(params[0][0].veg + params[0][1].veg + params[1][0].veg + params[1][1].veg);}
	void set_last_occluded(bool val) {last_occluded = val; last_occluded_frame = frame_counter;}
	bool was_last_occluded  () const {return (last_occluded_frame == frame_counter &&  last_occluded);}
	bool was_last_unoccluded() const {return (last_occluded_frame == frame_counter && !last_occluded);}

	point get_center() const {
		return point(get_xval(((x1+x2)>>1) + (xoff - xoff2)), get_yval(((y1+y2)>>1) + (yoff - yoff2)), 0.5*(mzmin + mzmax));
	}
	cube_t get_bcube() const {
		// Note: here we include the water plane's contribution to the upper z bound, since we draw the water as part of the tile contents
		float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2)), z2(max(get_tile_zmax(), water_plane_z));
		return cube_t(xv1-trmax, xv1+(x2-x1)*deltax+trmax, yv1-trmax, yv1+(y2-y1)*deltay+trmax, mzmin-BCUBE_ZTOLER, z2+BCUBE_ZTOLER);
	}
	cube_t get_mesh_bcube() const {
		float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2));
		return cube_t(xv1, xv1+(x2-x1)*deltax, yv1, yv1+(y2-y1)*deltay, mzmin-BCUBE_ZTOLER, mzmax+BCUBE_ZTOLER); // Note: bias by BCUBE_ZTOLER so dz != 0 when mzmin == mzmax
	}
	cube_t get_mesh_sub_bcube(unsigned x, unsigned y) const {
		assert(x < 4 && y < 4);
		float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2));
		float const xv2(xv1+(x2-x1)*deltax), yv2(yv1+(y2-y1)*deltay), dx((xv2 - xv1)/4), dy((yv2 - yv1)/4);
		return cube_t((xv1 + x*dx), (xv1 + (x+1)*dx), (yv1 + y*dy), (yv1 + (y+1)*dy), sub_zmin[y][x], sub_zmax[y][x]);
	}
	cube_t get_water_bcube() const {
		float const xv1(get_xval(wx1 + xoff - xoff2)), yv1(get_yval(wy1 + yoff - yoff2));
		return cube_t(xv1, xv1+(wx2-wx1)*deltax, yv1, yv1+(wy2-wy1)*deltay, water_plane_z, water_plane_z); // zero area in z
	}
	void fill_adj_mask(bool mask[3][3], int x, int y) const;
	float get_min_dist_to_pt(point const &pt, bool xy_only=0, bool mesh_only=1) const;
	float get_max_xy_dist_to_pt(point const &pt) const;
	bool contains_point(point const &pos) const {return get_bcube().contains_pt_xy(pos);}
	bool contains_camera() const {return contains_point(get_camera_pos());}
	unsigned get_gpu_mem() const;

	unsigned get_tree_mem() const { // only accounts for top-level class memory
		return (pine_trees.capacity()*sizeof(small_tree) + decid_trees.capacity()*sizeof(tree));
	}
	void clear();
	void clear_flowers() {flowers.clear();}
	void clear_shadows();
	void clear_shadow_map(tile_shadow_map_manager *smap_manager);
	void clear_vbo_tid(tile_shadow_map_manager *smap_manager);
	void invalidate_shadows() {shadows_invalid = 1;}
	void create_zvals(mesh_xy_grid_cache_t &height_gen);

	vector3d get_norm_not_normalized(unsigned ix) const {
		return vector3d(DY_VAL*(zvals[ix] - zvals[ix + 1]), DX_VAL*(zvals[ix] - zvals[ix + zvsize]), dxdy);
	}
	vector3d get_norm(unsigned ix) const {
		return get_norm_not_normalized(ix).get_norm();
	}

	// *** shadows ***
	void calc_mesh_ao_lighting();
	void calc_shadows_for_light(unsigned l);
	static void proc_tile_queue(tile_t *init_tile, unsigned l);
	void calc_shadows(bool calc_sun, bool calc_moon, bool no_push=0);

	tile_xy_pair get_tile_xy_pair(int dx=0, int dy=0) const {
		return tile_xy_pair((x1/(int)size)+dx, (y1/(int)size)+dy);
	}
	tile_t *get_adj_tile(int dx, int dy) const {
		return get_tile_from_xy(get_tile_xy_pair(dx, dy));
	}
	tile_t *get_adj_tile_smap(int dx, int dy) const {
		tile_t *adj_tile(get_adj_tile(dx, dy));
		return ((adj_tile && !adj_tile->tree_map.empty()) ? adj_tile : NULL);
	}
	void push_tree_ao_shadow(int dx, int dy, point const &pos, float tradius, float center_height) const;
	void add_tree_ao_shadow(point const &pos, float tradius, float center_height, bool no_adj_test);
	void apply_ao_shadows_for_trees(tile_t const *const tile, bool no_adj_test);
	void apply_tree_ao_shadows();
	void check_shadow_map_and_normal_texture();
	void upload_normal_texture(bool tid_is_valid);
	void upload_shadow_map_texture(bool tid_is_valid);
	void setup_shadow_maps(tile_shadow_map_manager &smap_manager);
	bool using_shadow_maps() const {return !smap_data.empty();}

	// *** mesh creation ***
	void ensure_height_tid();
	unsigned get_grass_block_dim() const {return (1+(size-1)/GRASS_BLOCK_SZ);} // ceil
	void create_texture(mesh_xy_grid_cache_t &height_gen);

	float get_rel_dist_to_camera(bool xy_dist=1) const {
		return max(0.0f, (xy_dist ? p2p_dist_xy(get_camera_pos(), get_center()) : p2p_dist(get_camera_pos(), get_center())) - radius)/get_scaled_tile_radius();
	}
	float get_bsphere_radius_inc_water() const;
	bool update_range(tile_shadow_map_manager &smap_manager);
	bool is_visible() const {return camera_pdu.sphere_and_cube_visible_test(get_center(), get_bsphere_radius_inc_water(), get_bcube());}
	float get_dist_to_camera_in_tiles(bool xy_dist=1) const {return get_rel_dist_to_camera(xy_dist)*TILE_RADIUS;}
	float get_scenery_thresh    (bool reflection_pass) const {return (reflection_pass ? SCENERY_THRESH_REF : SCENERY_THRESH);}
	float get_scenery_dist_scale(bool reflection_pass) const {return get_dist_to_camera_in_tiles(0)/get_scenery_thresh(reflection_pass);}
	float get_tree_dist_scale () const {return get_dist_to_camera_in_tiles()/get_tree_scale_denom();}
	float get_tree_far_weight () const {return (ENABLE_TREE_LOD ? CLIP_TO_01(GEOMORPH_THRESH*(get_tree_dist_scale() - 1.0f)) : 0.0);}

	// *** trees ***
	template <typename T> void postproc_trees(T const &trees, float &tzmax) { // pine/decidious trees
		tzmax  = mzmin;
		trees.update_zmax(tzmax);
		trmax  = max(trmax, trees.get_rmax());
		radius = max(radius, (calc_radius() + trmax)); // is this really needed?
	}
	void init_pine_tree_draw();
	void update_pine_tree_state(bool upload_if_needed);
	unsigned num_pine_trees() const {return pine_trees.size();}
	void draw_tree_leaves_lod(vector3d const &xlate, bool low_detail, int xlate_loc);
	void draw_pine_trees(shader_t &s, vector<vert_wrap_t> &trunk_pts, bool draw_branches, bool draw_near_leaves, bool draw_far_leaves, bool reflection_pass, int xlate_loc=-1);
	unsigned num_decid_trees() const {return decid_trees.size();}
	void gen_decid_trees_if_needed();
	void draw_decid_trees(shader_t &s, tree_lod_render_t &lod_renderer, bool draw_branches, bool draw_leaves, bool reflection_pass, bool shadow_pass, bool enable_smap);
	void update_decid_trees();

	// *** scenery/grass ***
	void update_scenery();
	void draw_scenery(shader_t &s, bool draw_opaque, bool draw_leaves, bool reflection_pass);
	void pre_draw_grass_flowers(shader_t &s, bool use_cloud_shadows) const;
	void draw_grass(shader_t &s, vector<vector<vector2d> > *insts, bool use_cloud_shadows, int lt_loc);
	void draw_flowers(shader_t &s, bool use_cloud_shadows);

	// *** clouds ***
	bool any_clouds_visible() const;
	void draw_tile_clouds(vpc_shader_t &s, bool reflection_pass);
	void gen_tile_clouds();

	// *** rendering ***
	void pre_draw(mesh_xy_grid_cache_t &height_gen);
	void shader_shadow_map_setup(shader_t &s, xform_matrix const *const mvm=nullptr) const;
	void bind_and_setup_shadow_map(shader_t &s) const;
	void bind_textures() const;
	void draw(shader_t &s, indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1], crack_ibuf_t const &crack_ibuf, bool reflection_pass, int shader_locs[2]) const;
	void draw_water_cap(shader_t &s, bool textures_already_set) const;
	void draw_water(shader_t &s, float z) const;
	bool is_water_visible() const;
	bool check_sphere_collision(point &pos, float radius) const;
	int get_tid_under_point(point const &pos) const;
	bool line_intersect_mesh(point const &v1, point const &v2, float &t, int &xpos, int &ypos) const;
}; // tile_t


class tile_draw_t : public indexed_vbo_manager_t {

	typedef map<tile_xy_pair, std::unique_ptr<tile_t> > tile_map;
	typedef set<tile_xy_pair> tile_set_t;
	typedef vector<pair<float, tile_t *> > draw_vect_t;

	tile_map tiles;
	unsigned ivbo_ixs[NUM_LODS+1];
	float terrain_zmin;
	draw_vect_t to_draw;
	vector<tile_t *> occluded_tiles;
	vector<vert_wrap_t> tree_trunk_pts;
	mesh_xy_grid_cache_t height_gen;
	lightning_strike_t lightning_strike;
	tree_lod_render_t lod_renderer;
	crack_ibuf_t crack_ibuf;
	tile_shadow_map_manager smap_manager;

	struct occluder_pts_t {
		point cube_pts[4];
		void calc_cube_top_points(cube_t const &bcube);
	};

public:
	tile_draw_t();
	~tile_draw_t() {/*clear();*/ height_gen.free_cshader();}
	void clear();
	void free_compute_shader() {height_gen.clear_context();}
	float update(float &min_camera_dist);
	static void setup_terrain_textures(shader_t &s, unsigned start_tu_id);
	static void shared_shader_lighting_setup(shader_t &s, unsigned lighting_shader);
	static void lighting_with_cloud_shadows_setup(shader_t &s, unsigned lighting_shader, bool cloud_shadows);
	void setup_mesh_draw_shaders(shader_t &s, bool reflection_pass, bool enable_shadow_map) const;
	bool can_have_reflection_recur(tile_t const *const tile, point const corners[3], tile_set_t &tile_set, unsigned dim_ix);
	bool can_have_reflection(tile_t const *const tile, tile_set_t &tile_set);
	void pre_draw();
	void draw(bool reflection_pass);
	void draw_tiles(bool reflection_pass, bool enable_shadow_map) const;
	void draw_shadow_pass(point const &lpos, tile_t *tile);
	void draw_water(shader_t &s, float zval) const;
	void end_lightning() const;
	static void set_noise_tex(shader_t &s, unsigned tu_id);
	static void set_tree_dither_noise_tex(shader_t &s, unsigned tu_id);
	static void set_pine_tree_shader(shader_t &s, string const &vs);
	void draw_pine_tree_bl(shader_t &s, bool branches, bool near_leaves, bool far_leaves, bool reflection_pass, int xlate_loc=-1);
	void draw_pine_trees(bool reflection_pass, bool shadow_pass=0);
	void draw_decid_tree_bl(shader_t &s, tree_lod_render_t &lod_renderer, bool branches, bool leaves, bool reflection_pass, bool shadow_pass, bool enable_smap);
	static void billboard_tree_shader_setup(shader_t &s);
	static void tree_branch_shader_setup(shader_t &s, bool enable_shadow_maps, bool enable_opacity);
	void draw_decid_trees(bool reflection_pass, bool shadow_pass=0);
	void draw_scenery(bool reflection_pass);
	static void setup_grass_flower_shader(shader_t &s, bool enable_wind, bool use_smap, float dist_const_mult);
	void draw_grass(bool reflection_pass);
	void draw_tile_clouds(bool reflection_pass);
	void update_lightning(bool reflection_pass);
	void clear_vbos_tids();
	void clear_flowers();
	tile_t *get_tile_from_xy(tile_xy_pair const &tp) const;
	tile_t *get_tile_containing_point(point const &pos) const;
	bool check_sphere_collision(point &pos, float radius) const;
	bool check_player_collision() const;
	int get_tid_under_point(point const &pos) const;
	bool line_intersect_mesh(point const &v1, point const &v2, float &t, tile_t *&intersected_tile, int &xpos, int &ypos) const;
	float get_actual_zmin() const;
}; // tile_draw_t


#endif // _TILED_TERRAIN_H_
