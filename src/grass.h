// 3D World - grass_t and grass_manager classes
// by Frank Gennari
// 7/9/12

#ifndef _GRASS_H_
#define _GRASS_H_

#include "3DWorld.h"


unsigned const NUM_GRASS_LODS       = 6;
unsigned const NUM_RND_GRASS_BLOCKS = 16;
unsigned const GRASS_BLOCK_SZ       = 4;
float const TT_GRASS_COLOR_SCALE    = 0.5;


class detail_scenery_t {
protected:
	unsigned vbo;

public:
	detail_scenery_t() : vbo(0) {}
	void free_vbo();
	static void setup_shaders_pre(shader_t &s);
	static void setup_shaders_post(shader_t &s);
};


class grass_manager_t : public detail_scenery_t {

protected:
	struct grass_t { // size = 44
		point p;
		vector3d dir, n;
		unsigned char c[3];
		unsigned char on_mesh : 1;
		float w;

		grass_t() {} // optimization
		grass_t(point const &p_, vector3d const &dir_, vector3d const &n_, unsigned char const *const c_, float w_, bool on_mesh_)
			: p(p_), dir(dir_), n(n_), on_mesh(on_mesh_), w(w_) {c[0] = c_[0]; c[1] = c_[1]; c[2] = c_[2];}
		void merge(grass_t const &g);
	};

	vector<grass_t> grass;
	bool data_valid;
	rand_gen_pregen_t rgen;
	typedef vert_norm_comp_color grass_data_t;

	vector3d interpolate_mesh_normal(point const &pos) const;

public:
	grass_manager_t() : data_valid(0) {}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~grass_manager_t() {clear();}
	size_t size() const {return grass.size ();} // 2 points per grass blade
	bool empty()  const {return grass.empty();}
	void clear();
	void add_grass_blade(point const &pos, float cscale, bool on_mesh);
	void create_new_vbo();
	void add_to_vbo_data(grass_t const &g, vector<grass_data_t> &data, unsigned &ix, vector3d &norm) const;
	void begin_draw() const;
	void end_draw() const;
};


class grass_tile_manager_t : public grass_manager_t {

	vector<unsigned> vbo_offsets[NUM_GRASS_LODS];
	unsigned start_render_ix, end_render_ix;

	void gen_block(unsigned bix);
	void gen_lod_block(unsigned bix, unsigned lod);

public:
	grass_tile_manager_t() : start_render_ix(0), end_render_ix(0) {}
	void clear();
	unsigned get_gpu_mem() const {return (vbo ? 3*size()*sizeof(grass_data_t) : 0);}
	void upload_data();
	void gen_grass();
	void update();
	void render_block(unsigned block_ix, unsigned lod, float density=1.0, unsigned num_instances=1);
};


class mesh_xy_grid_cache_t;

class flower_manager_t : public detail_scenery_t {

protected:
	struct flower_t { // size = 44
		point pos;
		vector3d normal;
		float radius;
		colorRGBA color;

		flower_t() {}
		flower_t(point const &p, vector3d const &n, float r, colorRGBA const &c) : pos(p), normal(n), radius(r), color(c) {}
	};

	vector<flower_t> flowers;
	rand_gen_t rgen;
	bool generated;

public:
	flower_manager_t() : generated(0) {}
	unsigned get_vertex_count() const {return 4*flowers.size();} // one quad (4 verts) per flower
	size_t size() const {return flowers.size ();}
	bool empty () const {return flowers.empty();}
	bool skip_generate() const;
	void clear() {free_vbo(); flowers.clear(); generated = 0;}
	void check_vbo();
	static void setup_flower_shader_post(shader_t &shader);
	void draw_triangles(shader_t &shader) const;
	void add_flowers(mesh_xy_grid_cache_t const density_gen[2], float grass_den, float hthresh, float dx, float dy, int xpos, int ypos, bool gen_zval);
	void gen_density_cache(mesh_xy_grid_cache_t density_gen[2], int x1, int y1);
};


class flower_tile_manager_t : public flower_manager_t {
public:
	void gen_flowers(vector<unsigned char> const &weight_data, unsigned wd_stride, int x1, int y1);
};


#endif // _GRASS_H_
