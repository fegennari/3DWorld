// 3D World - Heightmap Texture Managment
// by Frank Gennari
// 10/19/13

#ifndef _HEIGHTMAP_H_
#define _HEIGHTMAP_H_

#include "3DWorld.h"


class heightmap_t : public texture_t {
public:
	heightmap_t() {}
	heightmap_t(char t, char f, int w, int h, std::string const &n, bool inv) :
	texture_t(t, f, w, h, 0, 1, 0, n, inv) {}
	float get_heightmap_value(unsigned x, unsigned y) const;
	void modify_heightmap_value(unsigned x, unsigned y, int val, bool val_is_delta);
};


class tex_mod_map_manager_t {

public:
	typedef unsigned short tex_ix_t;
	typedef unsigned short hmap_val_t;
	unsigned max_tex_ix() const {return (1 << (sizeof(tex_ix_t) << 3));} // use a constant?

	struct tex_xy_t {
		tex_ix_t x, y;

		tex_xy_t() {x = y = 0;}
		tex_xy_t(tex_ix_t x_, tex_ix_t y_) {x = x; y = y;}
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
		mod_elem_t(pair<tex_xy_t, mod_map_val_t> const &v) : tex_xy_t(v.first), delta(v.second.val) {}
	};

	struct tex_mod_map_t : public map<tex_xy_t, mod_map_val_t> { // for uniquing/combining modifications to the same xy point
		// Note: this isn't entirely correct due to the clamping in the height texture update
		void add(mod_elem_t const &elem) {operator[](elem).val += elem.delta;}
	};

	typedef vector<mod_elem_t> tex_mod_vect_t;

protected:
	tex_mod_map_t mod_map;

public:
	void add_mod(mod_elem_t const &elem) {mod_map.add(elem);}
	void add_mod(tex_mod_vect_t const &mod);
	void add_mod(tex_mod_map_t const &mod);
	bool read_mod(std::string const &fn);
	bool write_mod(std::string const &fn) const;
};


class terrain_hmap_manager_t : public tex_mod_map_manager_t {

	heightmap_t hmap;

	bool clamp_xy(int &x, int &y) const;
	bool clamp_elem(mod_elem_t &elem);

public:
	void load(char const *const fn, bool invert_y=0);
	void maybe_load(char const *const fn, bool invert_y=0) {
		if (fn != NULL && !enabled()) {load(fn, invert_y);}
	}
	float get_clamped_height(int x, int y) const;
	float interpolate_height(float x, float y) const;
	vector3d get_norm(int x, int y) const;
	bool modify_height(mod_elem_t &elem);
	bool modify_and_cache_height(mod_elem_t &elem);
	bool read_mod(std::string const &fn);
	bool enabled() const {return hmap.is_allocated();}
	~terrain_hmap_manager_t() {hmap.free_data();}
};


#endif // _HEIGHTMAP_H_
