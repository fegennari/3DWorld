// 3D World - tiled terrain classes
// by Frank Gennari
// 4/4/13
#pragma once

#include "3DWorld.h"
#include "function_registry.h"
#include "inlines.h"
#include "mesh.h"
#include "small_tree.h"
#include "scenery.h"
#include "grass.h"
#include "tree_3dw.h"
#include "shadow_map.h"
#include "animals.h"
#include <unordered_map>
#include <unordered_set>


bool const ENABLE_TREE_LOD    = 1; // faster but has popping artifacts
bool const ENABLE_TERRAIN_ENV = 1;
bool const GRASS_CLOUD_SHADOWS= 1; // slow, but looks nice
bool const USE_TREE_BILLBOARDS= 1; // decidious trees: faster but lower quality
int  const TILE_RADIUS        = 6; // in mesh sizes
unsigned const NUM_LODS       = 5; // > 0
unsigned const NUM_SMAP_LODS  = 4;
float const TREE_LOD_THRESH   = 6.0;
float const GEOMORPH_THRESH   = 6.0;
float const PALM_DIST_SCALE   = 0.75;
float const SCENERY_THRESH_REF= 2.0;
float const SCENERY_THRESH    = 5.0;
float const BCUBE_ZTOLER      = 1.0E-6;


extern int frame_counter;
extern float grass_length, water_plane_z, tree_scale, far_clip_ratio;


class lightning_strike_t {

	int time=0;
	line3d path;
public:
	bool enabled() const {return (time > 0);}
	void clear() {path.points.clear(); time = 0;}
	point get_pos() const;
	void gen();
	void update();
	void draw() const;
	void end_draw() const;
};


struct ix_sz_pair {
	unsigned short ix=0, sz=0;
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

	int dxoff=0, dyoff=0;
	unsigned lod_level;
	tile_t *tile;

	tile_smap_data_t(unsigned tu_id_, unsigned smap_sz_, unsigned lod_level_, tile_t *tile_, smap_data_state_t const &init_state=smap_data_state_t())
		: smap_data_t(tu_id_, smap_sz_, init_state), lod_level(lod_level_), tile(tile_) {}
	virtual void render_scene_shadow_pass(point const &lpos);
	virtual bool needs_update(point const &lpos);
};


class tile_shadow_map_manager {

	vector<smap_data_state_t> free_list[NUM_LIGHT_SRC][NUM_SMAP_LODS];
public:
	tile_smap_data_t new_smap_data(unsigned tu_id, tile_t *tile, unsigned light, unsigned lod_level);
	void release_smap_data(tile_smap_data_t &smd, unsigned light);
	void clear_context();
	unsigned get_free_list_mem_usage() const;
};


inline float get_tile_width        () {return (X_SCENE_SIZE + Y_SCENE_SIZE);}
inline float get_scaled_tile_radius() {return TILE_RADIUS*get_tile_width();} // *far_clip_ratio
inline float get_tree_scale_denom  () {return max(1.0f, TREE_LOD_THRESH*calc_tree_size());}


struct tile_xy_pair {

	int x, y;
	tile_xy_pair(int x_=0, int y_=0) : x(x_), y(y_) {}
	bool operator==(tile_xy_pair const &tp) const {return (x == tp.x && y == tp.y);}
	bool operator< (tile_xy_pair const &tp) const {return ((y == tp.y) ? (x < tp.x) : (y < tp.y));}
	void operator+=(tile_xy_pair const &tp) {x += tp.x; y += tp.y;}
	void operator-=(tile_xy_pair const &tp) {x -= tp.x; y -= tp.y;}
	tile_xy_pair operator+(tile_xy_pair const &tp) const {return tile_xy_pair(x+tp.x, y+tp.y);}
	tile_xy_pair operator-(tile_xy_pair const &tp) const {return tile_xy_pair(x-tp.x, y-tp.y);}
};
struct hash_tile_xy_pair {
	uint32_t operator()(tile_xy_pair const &tp) const {return ((tp.x * 0x1F1F1F1F) ^ tp.y);}
};

tile_t *get_tile_from_xy(tile_xy_pair const &tp);


struct tile_cloud_t : public volume_part_cloud {

	float pos_hash=0.0;
	point pos;
	vector3d size; // {x, y, z}

	float get_rmax() const {return size.get_max_val();}
	cube_t get_bcube() const {cube_t bcube(pos, pos); bcube.expand_by(size); return bcube;}
	void draw(vpc_shader_t &s, vector3d const &xlate, float alpha_mult=1.0) const;
};

struct cloud_inst_t {
	tile_cloud_t const *cloud;
	float dist_sq, alpha;
	cloud_inst_t(tile_cloud_t const &cloud_, float dist_sq_, float alpha_=1.0) : cloud(&cloud_), dist_sq(dist_sq_), alpha(alpha_) {}
	bool operator<(cloud_inst_t const &ci) const {return (dist_sq < ci.dist_sq);}
};

typedef vector<cloud_inst_t> cloud_draw_list_t;

class tile_cloud_manager_t : public vector<tile_cloud_t> {

	bool generated=0;
	unsigned num_clouds=0;
	float cur_move_dist=0.0;
	cube_t bcube, range;
	rand_gen_t rgen;

	void update_bcube(tile_cloud_t const &c);
	void populate_clouds();
	void choose_num_clouds();
public:
	void gen_new_cloud();
	void gen(int x1, int y1, int x2, int y2);
	void move_by_wind(tile_t const &tile);
	void try_add_cloud(tile_cloud_t const &cloud);
	void get_draw_list(cloud_draw_list_t &clouds_to_draw, float mesh_zmin, float mesh_zmax) const;
};


class tile_t {
public:
	struct tree_map_val {
		unsigned char ao, sh;
		tree_map_val() : ao(255), sh(255) {}
	};
	mutable bool vis_ref_call=0; // cached visited flag for tile_draw_t::can_have_reflection_recur()
private:
	int x1=0, y1=0, x2=0, y2=0, wx1=0, wy1=0, wx2=0, wy2=0, last_occluded_frame=0;
	unsigned weight_tid=0, height_tid=0, normal_tid=0, shadow_tid=0;
	unsigned size=0, stride=0, zvsize=0, base_tsize=0, gen_tsize=0, smap_lod_level=0;
	float radius=0, mzmin=0, mzmax=0, mesh_dz=0, ptzmax=0, dtzmax=0, trmax=0, xstart=0, ystart=0, min_normal_z=0, deltax=0, deltay=0;
	bool sun_shadows_invalid=1, moon_shadows_invalid=1, recalc_tree_grass_weights=1, mesh_height_invalid=0, in_queue=0, last_occluded=0, has_any_grass=0;
	bool is_distant=0, no_trees=0, just_cleared=0, has_tunnel=0;
	colorRGB avg_mesh_tex_color;
	tile_offset_t mesh_off, ptree_off, dtree_off, scenery_off;
	float sub_zmin[4][4] = {0}, sub_zmax[4][4] = {0};
	vector<float> zvals, ao_zvals;
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
	vect_fish_t fish;
	vect_bird_t birds;
	vect_butterfly_t bflies;

	struct grass_block_t {
		unsigned ix; // 0 is unused
		float zmin, zmax;
		grass_block_t() : ix(0), zmin(0.0), zmax(0.0) {}
	};

	vector<grass_block_t> grass_blocks;

	struct terrain_params_t { // settings for different biomes
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
	float get_tile_zmax() const {return max((mzmax + (has_grass() ? grass_length : 0.0f)), max(ptzmax, dtzmax));}
	float get_zval(int x, int y) const {assert(!zvals.empty()); assert(x >= 0 && y >= 0 && x < (int)zvsize && y < (int)zvsize); return zvals[y*zvsize + x];}
	bool has_water() const {return (mzmin < water_plane_z);}
	bool all_water() const {return (mzmax < water_plane_z);} // get_tile_zmax()? - grass and trees should not be underwater
	bool can_have_trees() const {return (!no_trees && !is_distant && !all_water());}
	bool can_have_pine_palm_trees() const;
	bool can_have_decid_trees    () const;
	bool pine_trees_generated() const {return pine_trees.generated;}
	bool has_pine_trees() const {return (pine_trees_generated() && !pine_trees.empty());}
	bool has_valid_shadow_map() const {return !smap_data.empty();}
	bool has_grass() const {return !grass_blocks.empty();}
	bool get_checkerboard_bit() const {return (((x1/128) + (y1/128)) & 1);}
	void invalidate_mesh_height() {mesh_height_invalid = 1;}
	float get_avg_veg() const {return 0.25f*(params[0][0].veg + params[0][1].veg + params[1][0].veg + params[1][1].veg);}
	void set_last_occluded(bool val) {last_occluded = val; last_occluded_frame = frame_counter;}
	bool was_last_occluded  () const {return (last_occluded_frame == frame_counter &&  last_occluded);}
	bool was_last_unoccluded() const {return (last_occluded_frame == frame_counter && !last_occluded);}
	vect_bird_t      &get_birds () {return birds; } // for flocking
	vect_butterfly_t &get_bflies() {return bflies;} // for mating

	// all of these are in the current camera's local coordinate space (based on xoff/yoff/xoff2/yoff2)
	point get_center() const {
		return point(get_xval(((x1+x2)>>1) + (xoff - xoff2)), get_yval(((y1+y2)>>1) + (yoff - yoff2)), 0.5f*(mzmin + mzmax));
	}
	cube_t get_bcube() const {
		float const xv1(get_xval(x1 + xoff - xoff2)), yv1(get_yval(y1 + yoff - yoff2));
		float const z2(max(get_tile_zmax()+BCUBE_ZTOLER, water_plane_z)); // include the water plane's contribution, since we draw the water as part of the tile contents
		return cube_t(xv1-trmax, xv1+(x2-x1)*deltax+trmax, yv1-trmax, yv1+(y2-y1)*deltay+trmax, mzmin-BCUBE_ZTOLER, z2);
	}
	cube_t get_shadow_bcube() const;
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
	cube_t get_sub_bcube(unsigned x, unsigned y) const {
		cube_t bcube(get_mesh_sub_bcube(x, y));
		bcube.d[2][1] += (get_tile_zmax() - mzmax); // approximate accounting for grass and tree height
		return bcube;
	}
	cube_t get_water_bcube() const {
		float const xv1(get_xval(wx1 + xoff - xoff2)), yv1(get_yval(wy1 + yoff - yoff2));
		return cube_t(xv1, xv1+(wx2-wx1)*deltax, yv1, yv1+(wy2-wy1)*deltay, water_plane_z, water_plane_z); // zero area in z
	}
	// this is in global space
	cube_t get_mesh_bcube_global() const {
		float const xv1(get_xval(x1)), yv1(get_yval(y1));
		return cube_t(xv1, xv1+(x2-x1)*deltax, yv1, yv1+(y2-y1)*deltay, mzmin, mzmax);
	}
	void fill_adj_mask(bool mask[3][3], int x, int y) const;
	float get_min_dist_to_pt(point const &pt, bool xy_only=0, bool mesh_only=1) const;
	float get_max_xy_dist_to_pt(point const &pt) const;
	bool contains_point(point const &pos) const {return get_bcube().contains_pt_xy(pos);} // XY only
	bool contains_camera() const {return contains_point(get_camera_pos());}
	unsigned get_gpu_mem() const;
	unsigned get_smap_mem() const;
	unsigned count_shadow_maps() const;

	unsigned get_tree_mem() const { // only accounts for top-level class memory + palm verts
		return (get_cont_mem_usage(pine_trees) + get_cont_mem_usage(decid_trees) + pine_trees.palm_vbo_mem);
	}
	void clear();
	void clear_flowers() {flowers.clear();}
	void clear_shadows(bool clear_sun=1, bool clear_moon=1, bool no_clear_adj=0);
	void clear_shadow_map(tile_shadow_map_manager *smap_manager);
	void clear_vbo_tid(tile_shadow_map_manager *smap_manager);
	void clear_pine_tree_vbos() {pine_trees.clear_vbos();}
	bool create_zvals(mesh_xy_grid_cache_t &height_gen, bool no_wait);
	void get_z_minmax_for_area(point const &pos, float radius, float &zmin, float &zmax) const;
	float get_zval_at(float x, float y, bool in_global_space) const;

	vector3d get_norm_not_normalized(unsigned ix) const {
		return vector3d(DY_VAL*(zvals[ix] - zvals[ix + 1]), DX_VAL*(zvals[ix] - zvals[ix + zvsize]), dxdy);
	}
	vector3d get_norm(unsigned ix) const {return get_norm_not_normalized(ix).get_norm();}
	vector3d get_mesh_xlate() const {return mesh_off.get_xlate() + vector3d(xstart, ystart, 0.0);}

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
	void add_cloud(tile_cloud_t const &cloud) {clouds.try_add_cloud(cloud);}
	void push_tree_ao_shadow(int dx, int dy, point const &pos, float tradius) const;
	void add_tree_ao_shadow(point const &pos, float tradius, bool no_adj_test);
	template<typename T> void apply_ao_shadows_for_tree_group(T const &trees, tile_offset_t const &toff, bool no_adj_test, float rscale);
	void apply_ao_shadows_for_trees(tile_t const *const tile, bool no_adj_test);
	void apply_tree_ao_shadows();
	void check_shadow_map_and_normal_texture(bool no_push=0);
	void upload_normal_texture(bool tid_is_valid);
	void upload_shadow_map_texture(bool tid_is_valid);
	void draw_smap_debug_vis(shader_t &s) const;
	void setup_shadow_maps(tile_shadow_map_manager &smap_manager, bool cleanup_only);
	bool shadow_maps_allocated() const;
	bool using_shadow_maps() const {return !smap_data.empty();}

	// *** mesh creation ***
	void ensure_height_tid();
	unsigned get_grass_block_dim() const {return (1+(size-1)/GRASS_BLOCK_SZ);} // ceil
	void create_texture(mesh_xy_grid_cache_t &height_gen);
	void add_grass_block_at(unsigned x, unsigned y, float mhmin, float mhmax, unsigned grass_block_dim);
	void create_or_update_weight_tex();
	void calc_avg_mesh_color();

	float get_rel_dist_to_camera(bool xy_dist=1) const {
		return max(0.0f, (xy_dist ? p2p_dist_xy(get_camera_pos(), get_center()) : p2p_dist(get_camera_pos(), get_center())) - radius)/get_scaled_tile_radius();
	}
	bool rel_dist_to_camera_xy_lt(float rel_dist) const {
		return dist_xy_less_than(get_camera_pos(), get_center(), (rel_dist*get_scaled_tile_radius() + radius));
	}
	float get_bsphere_radius_inc_water() const;
	bool use_as_occluder() const;
	bool mesh_sphere_intersect(point const &pos, float rradius) const;
	bool update_range(tile_shadow_map_manager &smap_manager);
	bool is_visible() const {return camera_pdu.sphere_and_cube_visible_test(get_center(), get_bsphere_radius_inc_water(), get_bcube());}
	bool is_smap_bounds_visible() const {return camera_pdu.cube_visible(get_shadow_bcube());}
	float get_dist_to_camera_in_tiles(bool xy_dist=1) const {return get_rel_dist_to_camera(xy_dist)*TILE_RADIUS;}
	float get_scenery_thresh    (bool reflection_pass) const {return (reflection_pass ? SCENERY_THRESH_REF : SCENERY_THRESH);}
	float get_scenery_dist_scale(bool reflection_pass) const {return tree_scale*get_dist_to_camera_in_tiles(0)/get_scenery_thresh(reflection_pass);}
	float get_tree_dist_scale(bool has_palm=0) const {return (has_palm ? PALM_DIST_SCALE : 1.0f)*get_dist_to_camera_in_tiles()/get_tree_scale_denom();}
	float get_tree_far_weight(bool force_high_detail, bool has_palm=0) const {
		return ((ENABLE_TREE_LOD && !force_high_detail) ? CLIP_TO_01(GEOMORPH_THRESH*(get_tree_dist_scale(has_palm) - 1.0f)) : 0.0);
	}
	float get_draw_priority() const;

	// *** trees ***
	template <typename T> void postproc_trees(T const &trees, float &tzmax) { // pine/decidious trees
		tzmax  = mzmin;
		trees.update_zmax(tzmax);
		trmax  = max(trmax, trees.get_rmax());
		radius = max(radius, (calc_radius() + trmax)); // is this really needed?
	}
	void init_pine_tree_draw();
	void update_pine_tree_state(bool upload_if_needed, bool force_high_detail=0);
	unsigned num_pine_trees() const {return pine_trees.size();}
	void draw_tree_leaves_lod(shader_t &s, vector3d const &xlate, bool low_detail, int xlate_loc);
	void draw_pine_trees(shader_t &s, vector<tile_t *> &to_draw_trunk_pts, bool draw_trunks, bool draw_near_leaves, bool draw_far_leaves,
		bool shadow_pass, bool reflection_pass, bool enable_smap, int xlate_loc);
	void draw_trunk_pts(shader_t &s);
	unsigned num_decid_trees() const {return decid_trees.size();}
	void gen_decid_trees_if_needed();
	void set_mesh_ambient_color(shader_t &s) const;
	void draw_decid_trees(shader_t &s, tree_lod_render_t &lod_renderer, bool draw_branches, bool draw_leaves, bool reflection_pass, bool shadow_pass, bool enable_smap);
	void update_decid_trees();
	void register_tree_change(tile_shadow_map_manager &smap_manager);
	template <typename T> bool add_new_trees(T &trees, tile_offset_t const &toff, cube_t &update_bcube, float &tzmax, point const &tpos, float rradius, bool is_square);
	int  add_or_remove_trees_at(point const &pos, float rradius, bool add_trees, int brush_shape, tile_shadow_map_manager &smap_manager, cube_t &update_bcube);
	bool add_or_remove_grass_at(point const &pos, float rradius, bool add_grass, int brush_shape, float brush_weight);

	// *** scenery/grass ***
	void update_scenery();
	void draw_scenery(shader_t &s, shader_t &vrs, bool draw_opaque, bool draw_leaves, bool reflection_pass, bool shadow_pass=0, bool enable_shadow_maps=0);
	void pre_draw_grass_flowers(shader_t &s, bool use_cloud_shadows) const;
	unsigned draw_grass(shader_t &s, vector<vector<vector2d> > *insts, bool use_cloud_shadows, bool enable_tess, int lt_loc);
	unsigned draw_flowers(shader_t &s, bool use_cloud_shadows);
	bool choose_butterfly_dest(point &dest, sphere_t &plant_bsphere, rand_gen_t &rgen) const;

	// *** clouds ***
	void get_cloud_draw_list(cloud_draw_list_t &clouds_to_draw) const {clouds.get_draw_list(clouds_to_draw, mzmin, mzmax);}
	unsigned update_tile_clouds();

	// *** animals ***
	void add_animal(fish_t      const &f) {fish.push_back  (f);}
	void add_animal(bird_t      const &b) {birds.push_back (b);}
	void add_animal(butterfly_t const &b) {bflies.push_back(b);}
	template<typename A> void propagate_animals_to_neighbor_tiles(animal_group_t<A> &animals);
	void update_animals();
	void clear_animals() {fish.clear(); birds.clear(); bflies.clear();}
	void draw_birds (shader_t &s) const {birds.draw_animals (s, this);}
	void draw_fish  (shader_t &s) const {fish.draw_animals  (s, this);}
	void draw_bflies(shader_t &s) const {bflies.draw_animals(s, this);}

	// *** rendering ***
	void pre_draw(mesh_xy_grid_cache_t &height_gen);
	void shader_shadow_map_setup(shader_t &s, xform_matrix const *const mvm=nullptr) const;
	void bind_and_setup_shadow_map(shader_t &s) const;
	bool try_bind_shadow_map(shader_t &s, bool check_only) const;
	void bind_textures() const;
	void draw_mesh_vbo(indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1], unsigned lod_level) const;
	void draw(shader_t &s, indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1], crack_ibuf_t const &crack_ibuf, int reflection_pass, int shader_locs[2]) const;
	void draw_shadow_pass(shader_t &s, indexed_vbo_manager_t const &vbo_mgr, unsigned const ivbo_ixs[NUM_LODS+1]);
	void draw_water_cap(shader_t &s, bool textures_already_set) const;
	void draw_water(shader_t &s, float z) const;
	bool is_water_visible() const;
	bool check_sphere_collision(point &pos, float sradius, bool inc_dtrees=1, bool inc_ptrees=1, bool inc_scenery=1) const;
	bool check_cube_int_trees(cube_t const &c) const;
	int get_tid_under_point(point const &pos) const;
	bool line_intersect_mesh(point const &v1, point const &v2, float &t, int &xpos, int &ypos, float inc_trees) const;
}; // tile_t


class tile_draw_t : public indexed_vbo_manager_t {

	typedef unordered_map<tile_xy_pair, unique_ptr<tile_t>, hash_tile_xy_pair> tile_map;
	typedef vector<pair<float, tile_t *> > draw_vect_t;

	tile_map tiles;
	bool buildings_valid=0;
	unsigned ivbo_ixs[NUM_LODS+1] = {0};
	unsigned tiles_gen_prev_frame=0;
	float terrain_zmin=0.0;
	draw_vect_t to_draw, to_gen_zvals;
	vector<tile_t *> occluded_tiles, to_draw_trunk_pts;
	cloud_draw_list_t to_draw_clouds;
	vector<mesh_xy_grid_cache_t> height_gens;
	lightning_strike_t lightning_strike;
	tree_lod_render_t lod_renderer;
	crack_ibuf_t crack_ibuf;
	tile_shadow_map_manager smap_manager;
	vector<pair<float, tile_xy_pair>> shadow_recomp_queue;
	vect_cube_t occluder_cubes;
	vector<vector<vector2d> > grass_insts[NUM_GRASS_LODS];

	struct occluder_pts_t {
		point cube_pts[4];
		void calc_cube_top_points(cube_t const &bcube);
	};
	struct occluder_cubes_t {
		tile_t const *const tile;
		cube_t bcube, sub_cubes[16];
		occluder_cubes_t(tile_t const *const tile_);
	};
	vector<occluder_cubes_t> occluders; // reused across draw calls
	vector<unsigned> occluder_ixs; // reused across draw calls
	void insert_tile(tile_t *tile);

public:
	tile_draw_t();
	~tile_draw_t() {/*clear();*/}
	void clear(bool no_regen_buildings);
	void free_compute_shader();
	float update(float &min_camera_dist);
private:
	static void setup_terrain_textures(shader_t &s, unsigned start_tu_id);
	static void shared_shader_lighting_setup(shader_t &s, unsigned lighting_shader);
	static void lighting_with_cloud_shadows_setup(shader_t &s, unsigned lighting_shader, bool cloud_shadows);
	void setup_mesh_draw_shaders(shader_t &s, bool reflection_pass, bool enable_shadow_map) const;
	bool can_have_reflection_recur(tile_t const *const tile, point const corners[3], unsigned dim_ix);
	bool can_have_reflection(tile_t const *const tile);
public:
	uint64_t show_debug_stats(bool calc_mem_only) const;
	void pre_draw();
	void draw(int reflection_pass);
	void draw_shadow_pass(point const &lpos, tile_t *tile, bool decid_trees_only=0);
	void draw_smap_debug_vis() const;
	void draw_decid_tree_shadows() {draw_shadow_pass(camera_pdu.pos, nullptr, 1);}
	void draw_water(shader_t &s, float zval) const;
	static void billboard_tree_shader_setup(shader_t &s);
	static void tree_branch_shader_setup(shader_t &s, bool enable_shadow_maps, bool enable_opacity, bool shadow_only, bool enable_dlights=0);
private:
	void draw_tiles(int reflection_pass, bool enable_shadow_map) const;
	void draw_tiles_shadow_pass(point const &lpos, tile_t const *const tile);
	bool find_and_bind_any_valid_shadow_map(shader_t &s) const;
	static void set_noise_tex(shader_t &s, unsigned tu_id);
	static void set_tree_dither_noise_tex(shader_t &s, unsigned tu_id);
	static void set_pine_tree_shader(shader_t &s, string const &vs, bool use_texgen=1);
	static void set_pine_tree_shader_post(shader_t &s);
	void draw_pine_tree_bl(shader_t &s, bool branches, bool near_leaves, bool far_leaves, bool shadow_pass, bool reflection_pass, bool enable_smap, int xlate_loc);
	void draw_pine_trees(bool reflection_pass, bool shadow_pass=0);
	void draw_decid_tree_bl(shader_t &s, tree_lod_render_t &lod_renderer, bool branches, bool leaves, bool reflection_pass, bool shadow_pass, bool enable_smap);
	void draw_decid_trees(bool reflection_pass, bool shadow_pass=0);
	void draw_scenery(bool reflection_pass, bool shadow_pass=0);
	static void setup_grass_flower_shader(shader_t &s, bool enable_wind, bool use_smap, float dist_const_mult);
	void draw_grass(bool reflection_pass);
	void draw_animals(bool reflection_pass);
public:
	void draw_tile_clouds(bool reflection_pass);
	void update_lightning(bool reflection_pass);
	void end_lightning() const;
	void clear_vbos_tids();
	void clear_flowers();
	bool remove_tile(tile_xy_pair const &tp) {return (tiles.erase(tp) > 0);} // okay if tile doesn't exist; unused
	tile_t *get_tile_from_xy(tile_xy_pair const &tp) const;
	tile_t *get_tile_containing_point(point const &pos) const;
	void invalidate_tile_smap_in_region(cube_t region);
	bool try_bind_tile_smap_at_point(point const &pos, shader_t &s, bool check_only) const;
	bool check_sphere_collision(point &pos, float radius) const;
	bool check_cube_int_trees(cube_t const &c) const;
	bool check_player_collision() const;
	int get_tid_under_point(point const &pos) const;
	bool line_intersect_mesh(point const &v1, point const &v2, float &t, tile_t *&intersected_tile, int &xpos, int &ypos, float inc_trees) const;
	float get_actual_zmin() const;
	void add_or_remove_trees_at(point const &pos, float radius, bool add_trees, int brush_shape);
	void add_or_remove_grass_at(point const &pos, float radius, bool add_grass, int brush_shape, float brush_weight);
}; // tile_draw_t

