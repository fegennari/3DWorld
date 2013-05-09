// 3D World - grass_t and grass_manager classes
// by Frank Gennari
// 7/9/12

#ifndef _GRASS_H_
#define _GRASS_H_

#include "3DWorld.h"


unsigned const NUM_GRASS_LODS    = 6;
unsigned const GRASS_BLOCK_SZ    = 4;
float const TT_GRASS_COLOR_SCALE = 0.5;


class grass_manager_t {

protected:
	struct grass_t { // size = 48
		point p;
		vector3d dir, n;
		unsigned char c[3];
		unsigned char shadowed;
		float w;

		grass_t() {} // optimization
		grass_t(point const &p_, vector3d const &dir_, vector3d const &n_, unsigned char const *const c_, float w_)
			: p(p_), dir(dir_), n(n_), shadowed(0), w(w_) {c[0] = c_[0]; c[1] = c_[1]; c[2] = c_[2];}
		void merge(grass_t const &g);
	};

	vector<grass_t> grass;
	unsigned vbo;
	bool data_valid;
	rand_gen_t rgen;
	typedef vert_norm_comp_tc_color grass_data_t; // could use vert_norm_comp_tc_comp_color, it takes less vmem but is slightly slower

public:
	grass_manager_t() : vbo(0), data_valid(0) {}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~grass_manager_t() {clear();}
	size_t size() const {return grass.size ();} // 2 points per grass blade
	bool empty()  const {return grass.empty();}
	void free_vbo();
	void clear();
	void add_grass_blade(point const &pos, float cscale);
	void create_new_vbo();
	void add_to_vbo_data(grass_t const &g, vector<grass_data_t> &data, unsigned &ix, vector3d &norm) const;
	void begin_draw(float spec_weight) const;
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
	void render_block(unsigned block_ix, unsigned lod);
};


#endif // _GRASS_H_
