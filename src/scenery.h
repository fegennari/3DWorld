// 3D World - Scenery Classes Header (plants, rocks, logs stumps)
// by Frank Gennari
// 6/26/12
#pragma once

#include "3DWorld.h"
#include "shape_line3d.h"
#include "upsurface.h"
#include "draw_utils.h"
#include "voxels.h"


enum {PLANT_MJ=0, PLANT1, PLANT2, PLANT3, PLANT4, COFFEE, SEAWEED, NUM_PLANT_TYPES};
enum {LEAFY_PLANT_UW=0, LEAFY_PLANT_DIRT, LEAFY_PLANT_GRASS, LEAFY_PLANT_ROCK, NUM_LEAFY_PLANT_TYPES};
unsigned const NUM_LAND_PLANT_TYPES  = 6;
unsigned const NUM_WATER_PLANT_TYPES = NUM_PLANT_TYPES - NUM_LAND_PLANT_TYPES;

struct plant_type {
	int tid;
	float leaf_length, leaf_width_base, leaf_width_end;
	colorRGBA stemc, leafc, berryc;

	plant_type(int tid_, colorRGBA const &sc, colorRGBA const &lc, colorRGBA const &bc, float leaf_length_=1.0, float leaf_width_base_=1.0, float leaf_width_end_=1.0)
		: tid(tid_), leaf_length(leaf_length_), leaf_width_base(leaf_width_base_), leaf_width_end(leaf_width_end_), stemc(sc), leafc(lc), berryc(bc) {}
};


struct surface_cache {
	typedef pair<long, long> seed_pair;
	typedef map<seed_pair, p_upsurface> surface_map;
	surface_map scache;

	p_upsurface get_surface(bool fixed_sz_rock_cache);
	void clear() {scache.clear();}
	void clear_unref();
};

class surface_rock : public scenery_obj { // size = 1456+
	int vbo_mgr_ix=-1;
	float scale=0.0;
	vector3d dir;
	p_upsurface surface;
public:
	void create(int x, int y, int use_xy, bool fixed_sz_rock_cache);
	void gen_points(vbo_vnt_block_manager_t &vbo_manager);
	unsigned get_num_verts() const;
	void add_cobjs();
	void draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val, vbo_vnt_block_manager_t &vbo_manager) const;
	void destroy();
};

class s_rock : public scenery_obj { // size = 48
	float size=0.0, angle=0.0;
	vector3d scale, dir;
public:
	void create(int x, int y, int use_xy);
	void add_cobjs();
	void draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val) const;
	void add_bounds_to_bcube(cube_t &bcube) const {scenery_obj::add_bounds_to_bcube(bcube, 1.3*radius);}
};

class voxel_rock_manager_t {
	vector<unique_ptr<voxel_model_rock>> models;
	set<unsigned> to_gen;
	noise_texture_manager_t ntg;
public:
	unsigned gen_model_ix(int rseed);
	void build_models(unsigned num_lod_levels);
	voxel_model_rock       &get_model(unsigned ix)       {assert(ix < models.size()); assert(models[ix] != nullptr); return *models[ix];}
	voxel_model_rock const &get_model(unsigned ix) const {assert(ix < models.size()); assert(models[ix] != nullptr); return *models[ix];}
	void free_context();
};

class voxel_rock : public scenery_obj {
	unsigned model_ix=0;
	int rseed=0;

	unsigned get_tid() const;
public:
	void create(int x, int y, int use_xy);
	void build_model();
	void add_cobjs();
	void draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val, shader_t &s, bool use_model_texgen);
};

class burnable_scenery_obj : public scenery_obj {
protected:
	float fire_amt=0.0, burn_amt=0.0;
public:
	virtual ~burnable_scenery_obj() {}
	virtual float get_bsphere_radius() const = 0;
	virtual point get_center() const = 0;
	void next_frame();
	void draw_fire(fire_drawer_t &fire_drawer, float rscale, unsigned ix) const;
};

class wood_scenery_obj : public burnable_scenery_obj {
protected:
	int closest_tree_type=-1; // cached from drawing functions
	void calc_type();
	bool is_from_large_trees() const;
	int get_tid() const;
	colorRGBA get_bark_color(vector3d const &xlate=zero_vector) const;
public:
	void cache_closest_tree_type(tree_cont_t const &trees);
};

class s_log : public wood_scenery_obj { // size = 57 (60)
	float length=0.0, radius2=0.0;
	point pt2;
	vector3d dir;
public:
	virtual point get_center() const {return (pos + pt2)*0.5;}
	bool check_sphere_coll(point &center, float sphere_radius) const {return 0;} // no collisions
	void shift_by(vector3d const &vd);
	int create(int x, int y, int use_xy, float minz);
	void add_cobjs();
	void draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val) const;
	bool update_zvals(int x1, int y1, int x2, int y2);
	virtual float get_bsphere_radius() const {return max(length, max(radius, radius2));}
	void add_bounds_to_bcube(cube_t &bcube) const {scenery_obj::add_bounds_to_bcube(bcube, get_bsphere_radius());}
};

class s_stump : public wood_scenery_obj { // size = 29 (32)
	float radius2=0.0, height=0.0;
public:
	virtual point get_center() const {return (pos + point(0.0, 0.0, 0.5*height));}
	int create(int x, int y, int use_xy, float minz);
	void add_cobjs();
	bool check_sphere_coll(point &center, float sphere_radius) const;
	void draw(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate, float scale_val) const;
	virtual float get_bsphere_radius() const {return max(height, max(radius, radius2));}
	void add_bounds_to_bcube(cube_t &bcube) const {scenery_obj::add_bounds_to_bcube(bcube, get_bsphere_radius());}
};

struct texture_binder_t {
	int cur_tid=-1;
	void do_bind(int tid);
};

struct plant_base : public burnable_scenery_obj { // size = 32
	struct shader_state_t {
		int color_scale_loc, normal_scale_loc, wind_scale_loc, wind_add_loc;
		float wind_scale;
		texture_binder_t texture_binder;
		shader_state_t() : color_scale_loc(-1), normal_scale_loc(-1), wind_scale_loc(-1), wind_add_loc(-1), wind_scale(1.0) {}
		void set_color_scale(shader_t &s, colorRGBA const &color);
		void set_normal_scale(shader_t &s, float normal_scale);
		void set_wind_scale(shader_t &s, float wscale);
		void set_wind_add(shader_t &s, float w_add);
	};
	int vbo_mgr_ix=-1;

	bool operator<(plant_base const &p) const {return (type < p.type);}
	int create(int x, int y, int use_xy, float minz);
	void next_frame();
	colorRGBA get_plant_color(vector3d const &xlate) const;
	virtual point get_center() const {return pos;}
};

class s_plant : public plant_base { // size = 56
	int coll_id2=-1;
	float height=1.0;
	vector<vert_wrap_t> berries;
public:
	virtual float get_bsphere_radius() const {return 0.5f*(height + radius);}
	point get_top_pt() const {return pos + point(0.0, 0.0, height);}
	bool is_water_plant() const;
	colorRGBA const &get_leaf_color() const;
	colorRGBA const &get_stem_color() const;
	int get_leaf_tid() const;
	int create(int x, int y, int use_xy, float minz);
	void create2(point const &pos_, float height_, float radius_, int type_, int calc_z);
	void create_no_verts(point const &pos_, float height_, float radius_, int type_, int calc_z=0, bool land_plants_only=0, bool water_plants_only=0);
	void place_in_pond(cube_t const &pond);
	void add_cobjs();
	bool check_sphere_coll(point &center, float sphere_radius) const;
	void create_leaf_points(vector<vert_norm_comp> &points, float plant_scale, float nlevels_scale=1.0, unsigned nrings=3) const;
	void gen_points(vbo_vnc_block_manager_t &vbo_manager, vector<vert_norm_comp> &pts);
	void update_points_vbo(vbo_vnc_block_manager_t &vbo_manager);
	bool update_zvals(int x1, int y1, int x2, int y2, vbo_vnc_block_manager_t &vbo_manager);
	bool is_shadowed() const;
	void draw_stem(float sscale, bool shadow_only, bool reflection_pass, vector3d const &xlate) const;
	void draw_leaves(shader_t &s, vbo_vnc_block_manager_t &vbo_manager, bool shadow_only, bool reflection_pass, vector3d const &xlate, shader_state_t &state) const;
	void draw_berries(shader_t &s, vector3d const &xlate) const;
	void remove_cobjs();
	void write_to_cobj_file(std::ostream &out) const;
	void add_bounds_to_bcube(cube_t &bcube) const {scenery_obj::add_bounds_to_bcube(bcube, (height + radius));}
};

class leafy_plant : public plant_base {
	unsigned plant_ix=0;
	float delta_z=0.0, motion_amt=0.0, cur_motion_energy=0.0, prev_motion_energy=0.0;
	struct plant_leaf {xform_matrix m;};
	vector<plant_leaf> leaves;

	void gen_leaves();
public:
	virtual float get_bsphere_radius() const {return radius;}
	int create(int x, int y, int use_xy, float minz, unsigned plant_ix_);
	void create2(point const &pos_, float radius_, int type_, int calc_z, unsigned plant_ix_);
	unsigned num_leaves() const {return leaves.size();}
	point get_top_pt() const {return pos + point(0.0, 0.0, radius);}
	void gen_points(vbo_vnt_block_manager_t &vbo_manager, vector<vert_norm_tc> const &sphere_verts);
	void add_cobjs();
	void obj_collision(float energy) {cur_motion_energy += energy;}
	void next_frame();
	bool update_zvals(int x1, int y1, int x2, int y2, vbo_vnt_block_manager_t &vbo_manager);
	int get_tid() const;
	void draw_leaves(shader_t &s, bool shadow_only, bool reflection_pass, vector3d const &xlate, shader_state_t &state, vbo_vnt_block_manager_t &vbo_manager) const;
};


class scenery_group {
	vector<rock_shape3d> rock_shapes;
	vector<surface_rock> surface_rocks;
	vector<voxel_rock>   voxel_rocks;
	vector<s_rock>       rocks;
	vector<s_log>        logs;
	vector<s_stump>      stumps;
	vector<s_plant>      plants;
	vector<leafy_plant>  leafy_plants;
	vbo_vnc_block_manager_t plant_vbo_manager;
	vbo_vnt_block_manager_t rock_vbo_manager;
	vbo_vnt_block_manager_t leafy_vbo_manager;
	cube_t all_bcube;
	vector<vert_norm_comp> temp_pts;
public:
	bool generated=0;

	void clear_vbos();
	void clear();
	void free_cobjs();
	void add_cobjs();
	bool check_sphere_coll(point &center, float radius) const;
	void shift(vector3d const &vd);
	void calc_bcube();
	bool update_zvals(int x1, int y1, int x2, int y2);
	void do_rock_damage(point const &pos, float radius, float damage);
	void add_plant(point const &pos, float height, float radius, int type, int calc_z);
	void add_leafy_plant(point const &pos, float radius, int type, int calc_z);
	void gen(int x1, int y1, int x2, int y2, float vegetation_, bool fixed_sz_rock_cache, tree_cont_t const &trees);
	void post_gen_setup(tree_cont_t const &trees);
	void draw_plant_leaves(shader_t &s, bool shadow_only, vector3d const &xlate, bool reflection_pass=0);
	void draw_opaque_objects(shader_t &s, shader_t &vrs, bool shadow_only, vector3d const &xlate, bool draw_pld, float scale_val=0.0, bool reflection_pass=0);
	bool setup_voxel_rocks_shader(shader_t &vrs, bool shadow_only) const;
	void draw(bool shadow_only, vector3d const &xlate=zero_vector);
	void draw_fires(shader_t &s) const;
	void leafy_plant_coll(unsigned plant_ix, float energy);
	bool choose_butterfly_dest(point &dest, sphere_t &plant_bsphere, rand_gen_t &rgen) const;
	void write_plants_to_cobj_file(std::ostream &out) const;
	size_t get_gpu_mem() const {return (plant_vbo_manager.get_gpu_mem() + rock_vbo_manager.get_gpu_mem() + leafy_vbo_manager.get_gpu_mem());} // only accounts for part of memory
};

