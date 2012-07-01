// 3D World - Scenery Classes Header (plants, rocks, logs stumps)
// by Frank Gennari
// 6/26/12

#ifndef _SCENERY_H_
#define _SCENERY_H_

#include "3DWorld.h"
#include "shape_line3d.h"
#include "upsurface.h"
#include "draw_utils.h"


int get_bark_tex_for_tree_type(int type);


enum {PLANT_MJ = 0, PLANT1, PLANT2, PLANT3, PLANT4, COFFEE, NUM_PLANT_TYPES};

struct plant_type {

	int tid;
	colorRGBA stemc, leafc;

	plant_type(int tid_, colorRGBA const &sc, colorRGBA const &lc) : tid(tid_), stemc(sc), leafc(lc) {}
};


struct surface_cache {

	typedef pair<long, long> seed_pair;
	typedef map<seed_pair, upsurface *> surface_map;
	surface_map scache;

	upsurface *get_surface();
	void clear();
	void clear_unref();
};


class surface_rock : public scenery_obj { // size = 1456+

	int vbo_mgr_ix;
	float scale;
	vector3d dir;
	upsurface *surface;

public:
	surface_rock() : vbo_mgr_ix(-1), surface(NULL) {}
	void create(int x, int y, int use_xy, vbo_quad_block_manager_t &vbo_manager);
	void add_cobjs();
	void draw(float sscale, bool shadow_only, vector3d const &xlate, vbo_quad_block_manager_t &vbo_manager) const;
	void destroy();
};


class s_rock : public scenery_obj { // size = 48

	float size, angle;
	vector3d scale, dir;

public:
	void create(int x, int y, int use_xy);
	void add_cobjs();
	void draw(float sscale, bool shadow_only, vector3d const &xlate) const;
};


class s_log : public scenery_obj { // size = 57 (60)

	float length, radius2;
	point pt2;
	vector3d dir;

	int get_tid() const {return get_bark_tex_for_tree_type(type);}

public:
	void shift_by(vector3d const &vd);
	int create(int x, int y, int use_xy, float minz);
	void add_cobjs();
	void draw(float sscale, bool shadow_only, vector3d const &xlate) const;
	bool update_zvals(int x1, int y1, int x2, int y2);
};


class s_stump : public scenery_obj { // size = 29 (32)

	float radius2, height;

	int get_tid() const {return get_bark_tex_for_tree_type(type);}

public:
	int create(int x, int y, int use_xy, float minz);
	void add_cobjs();
	void draw(float sscale, bool shadow_only, vector3d const &xlate) const;
};


class s_plant : public scenery_obj { // size = 40

	int coll_id2, vbo_mgr_ix;
	float height;

public:
	s_plant() : coll_id2(-1), vbo_mgr_ix(-1), height(1.0) {}
	bool operator<(s_plant const &p) const {return (type < p.type);}
	int create(int x, int y, int use_xy, float minz, vbo_quad_block_manager_t &vbo_manager);
	void create2(point const &pos_, float height_, float radius_, int type_, int calc_z, vbo_quad_block_manager_t &vbo_manager);
	void add_cobjs();
	void gen_points(vbo_quad_block_manager_t &vbo_manager);
	bool is_shadowed() const;
	void draw_stem(float sscale, bool shadow_only, vector3d const &xlate) const;
	void draw_leaves(shader_t &s, vbo_quad_block_manager_t &vbo_manager, bool shadow_only, vector3d const &xlate) const;
	void remove_cobjs();
	void destroy();
};


class scenery_group {

	vector<rock_shape3d> rock_shapes;
	vector<surface_rock> surface_rocks;
	vector<s_rock>       rocks;
	vector<s_log>        logs;
	vector<s_stump>      stumps;
	vector<s_plant>      plants;
	vbo_quad_block_manager_t vbo_manager;

public:
	bool generated;

	scenery_group() : generated(0) {}
	void clear_vbos_and_dlists();
	void clear();
	void free();
	void add_cobjs();
	void shift(vector3d const &vd);
	void update_zvals(int x1, int y1, int x2, int y2);
	void do_rock_damage(point const &pos, float radius, float damage);
	void add_plant(point const &pos, float height, float radius, int type, int calc_z);
	void gen(int x1, int y1, int x2, int y2);
	void draw_plant_leaves(shader_t &s, bool shadow_only, vector3d const &xlate);
	void draw_opaque_objects(bool shadow_only, vector3d const &xlate, bool draw_pld);
	void draw(bool draw_opaque, bool draw_transparent, bool shadow_only, vector3d const &xlate=zero_vector);
};


#endif // _SCENERY_H_

