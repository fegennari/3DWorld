// 3D World - Heightmap Texture Managment
// by Frank Gennari
// 10/19/13

#ifndef _HEIGHTMAP_H_
#define _HEIGHTMAP_H_

#include "3DWorld.h"


enum {BSHAPE_CONST_SQ=0, BSHAPE_CNST_CIR, BSHAPE_LINEAR, BSHAPE_QUADRATIC, BSHAPE_COSINE, BSHAPE_SINE, BSHAPE_FLAT_SQ, BSHAPE_FLAT_CIR, NUM_BSHAPES};

float scale_mh_texture_val(float val);


class heightmap_t : public texture_t {

	unsigned get_pixel_ix(unsigned x, unsigned y) const;

public:
	heightmap_t() {}
	heightmap_t(char t, char f, int w, int h, std::string const &n, bool inv) :
	texture_t(t, f, w, h, 0, 1, 0, n, inv) {}
	unsigned get_pixel_value (unsigned x, unsigned y) const;
	float get_heightmap_value(unsigned x, unsigned y) const;
	void modify_heightmap_value(unsigned x, unsigned y, int val, bool val_is_delta);
};


class tex_mod_map_manager_t {

public:
	typedef unsigned short tex_ix_t;
	typedef int hmap_val_t; // really needs to be from -2^16 to 2^16, so we need 17 bits
	unsigned max_tex_ix() const {return (1 << (sizeof(tex_ix_t) << 3));} // use a constant?

	struct tex_xy_t {
		tex_ix_t x, y;

		tex_xy_t() {x = y = 0;}
		tex_xy_t(tex_ix_t x_, tex_ix_t y_) : x(x_), y(y_) {}
		bool operator==(tex_xy_t const &t) const {return (x == t.x && y == t.y);}
		bool operator< (tex_xy_t const &t) const {return ((x == t.x) ? (y < t.y) : (x < t.x));}
	};

	struct mod_map_val_t {
		hmap_val_t val;
		mod_map_val_t(hmap_val_t v=0) : val(v) {} // init to zero
	};

	struct mod_elem_t : public tex_xy_t {
		hmap_val_t delta;
		mod_elem_t() {}
		mod_elem_t(tex_ix_t x_, tex_ix_t y_, hmap_val_t d) : tex_xy_t(x_, y_), delta(d) {}
		mod_elem_t(pair<tex_xy_t, mod_map_val_t> const &v) : tex_xy_t(v.first), delta(v.second.val) {}
	};

	struct tex_mod_map_t : public map<tex_xy_t, mod_map_val_t> { // for uniquing/combining modifications to the same xy point
		// Note: this isn't entirely correct due to the clamping in the height texture update
		void add(mod_elem_t const &elem) {operator[](elem).val += elem.delta;}
	};

	struct hmap_brush_t {
		int x, y;
		unsigned radius;
		hmap_val_t delta;
		short shape;

		hmap_brush_t() : x(0), y(0), radius(0), delta(0), shape(0) {}
		hmap_brush_t(int x_, int y_, hmap_val_t d, unsigned r, short s) : x(x_), y(y_), delta(d), radius(r), shape(s) {assert(shape < NUM_BSHAPES);}
		bool is_flatten_brush() const {return (shape == BSHAPE_FLAT_SQ || shape == BSHAPE_FLAT_CIR);}
		void apply(tex_mod_map_manager_t *tmmm, int step_sz=1, unsigned num_steps=1) const;
	};

	typedef vector<mod_elem_t> tex_mod_vect_t;
	typedef vector<hmap_brush_t> brush_vect_t;

protected:
	tex_mod_map_t mod_map;
	brush_vect_t brush_vect;

public:
	void add_mod(mod_elem_t const &elem) {mod_map.add(elem);}
	void add_mod(tex_mod_vect_t const &mod);
	void add_mod(tex_mod_map_t const &mod);
	void apply_brush(hmap_brush_t const &brush, int step_sz=1, unsigned num_steps=1) {brush.apply(this, step_sz, num_steps);}

	void apply_and_cache_brush(hmap_brush_t const &brush, int step_sz=1, unsigned num_steps=1) {
		brush_vect.push_back(brush);
		apply_brush(brush, step_sz, num_steps);
	}
	bool pop_last_brush(hmap_brush_t &last_brush);
	bool undo_last_brush(); // unused
	bool read_mod(std::string const &fn);
	bool write_mod(std::string const &fn) const;

	virtual bool modify_height_value(int x, int y, hmap_val_t val, bool is_delta, float fract_x=0.0, float fract_y=0.0) = 0;
	virtual ~tex_mod_map_manager_t() {}
};


class terrain_hmap_manager_t : public tex_mod_map_manager_t {

	heightmap_t hmap;

public:
	void load(char const *const fn, bool invert_y=0);
	bool maybe_load(char const *const fn, bool invert_y=0);
	bool clamp_xy(int &x, int &y, float fract_x=0.0, float fract_y=0.0) const;
	bool clamp_no_scale(int &x, int &y) const;
	hmap_val_t get_clamped_pixel_value(int x, int y) const;
	float get_raw_height(int x, int y) const {return scale_mh_texture_val(hmap.get_heightmap_value(x, y));}
	float get_clamped_height(int x, int y) const;
	float interpolate_height(float x, float y) const;
	vector3d get_norm(int x, int y) const;

	virtual bool modify_height_value(int x, int y, hmap_val_t val, bool is_delta, float fract_x=0.0, float fract_y=0.0) { // unused
		assert(fract_x == 0.0 && fract_y == 0.0);
		modify_height(mod_elem_t(x, y, val), is_delta);
		return 1;
	}
	void modify_height(mod_elem_t const &elem, bool is_delta);
	void modify_and_cache_height(mod_elem_t const &elem, bool is_delta) {modify_height(elem, is_delta); add_mod(elem);} // unused
	hmap_val_t scale_delta(float delta) const;
	bool read_mod(std::string const &fn);
	bool enabled() const {return hmap.is_allocated();}
	~terrain_hmap_manager_t() {hmap.free_data();}
};


struct hmap_brush_param_t {

	unsigned delay; // placement delay in ticks
	int shape;
	int radius_exp;
	unsigned delta_exp;

	hmap_brush_param_t() : delay(1), shape(BSHAPE_COSINE), radius_exp(5), delta_exp(4) {}
	unsigned get_radius() const {return ((radius_exp < 0) ? 0 : (1 << radius_exp));}
	float get_delta_mag() const {return 0.001*pow(2.0f, (float)delta_exp);}
};


#endif // _HEIGHTMAP_H_
