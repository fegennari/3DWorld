// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16
#pragma once

#include "3DWorld.h"
#include "function_registry.h"
#include "inlines.h"

class tile_t;

struct tile_offset_t {
	int dxoff, dyoff;

	tile_offset_t(int dxoff_=0, int dyoff_=0) : dxoff(dxoff_), dyoff(dyoff_) {}
	void set_from_xyoff2() {dxoff = -xoff2; dyoff = -yoff2;}
	int get_delta_xoff() const {return ((xoff - xoff2) - dxoff);}
	int get_delta_yoff() const {return ((yoff - yoff2) - dyoff);}
	vector3d get_xlate() const {return vector3d(get_delta_xoff()*DX_VAL, get_delta_yoff()*DY_VAL, 0.0);}
	vector3d subtract_from(tile_offset_t const &o) const {return vector3d((o.dxoff - dxoff)*DX_VAL, (o.dyoff - dyoff)*DY_VAL, 0.0);}
};


class animal_t : public sphere_t {

public:
	vector3d velocity;

protected:
	bool enabled;
	vector3d dir;
	colorRGBA color;

	int get_ndiv(point const &pos_) const;
	void gen_dir_vel(rand_gen_t &rgen, float speed);
public:
	animal_t() : enabled(0), color(BLACK) {}
	void apply_force(vector3d const &force) {velocity += force;}
	void apply_force_xy(vector3d const &force) {velocity.x += force.x; velocity.y += force.y;}
	bool is_enabled() const {return enabled;}
	bool is_visible(point const &pos_, float vis_dist_scale=1.0) const;
	point get_draw_pos() const;
};

class fish_t : public animal_t {

	float get_mesh_zval_at_pos(tile_t const *const tile) const;
	float get_half_height() const {return 0.4*radius;} // approximate
public:
	static bool type_enabled();
	bool gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile);
	bool update(rand_gen_t &rgen, tile_t const *const tile);
	void draw(shader_t &s) const;
};

class bird_t : public animal_t {

	bool flocking;
	float time;
public:
	bird_t() : flocking(0), time(0) {}
	static bool type_enabled() {return 1;} // no model, always enabled
	bool gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile);
	bool update(rand_gen_t &rgen, tile_t const *const tile);
	void apply_force_xy_const_vel(vector3d const &force);
	void draw(shader_t &s) const;
};

class butterfly_t : public animal_t {

	float time;
public:
	butterfly_t() : time(0) {}
	static bool type_enabled();
	bool gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile);
	bool update(rand_gen_t &rgen, tile_t const *const tile);
	void draw(shader_t &s) const;
};


class animal_group_base_t {
protected:
	rand_gen_t rgen;
	bool generated;
public:
	animal_group_base_t() : generated(0) {}
	bool was_generated() const {return generated;}
};

template<typename A> class animal_group_t : public vector<A>, public animal_group_base_t {

	cube_t bcube;
public:
	void gen(unsigned num, cube_t const &range, tile_t const *const tile);
	void update(tile_t const *const tile);
	void remove(unsigned ix);
	void remove_disabled();
	void draw_animals(shader_t &s) const;
	void clear() {vector<A>::clear(); generated = 0;}
};

struct vect_fish_t : public animal_group_t<fish_t> {
	static void begin_draw(shader_t &s);
	static void end_draw  (shader_t &s);
};

struct vect_bird_t : public animal_group_t<bird_t> {
	void flock(tile_t const *const tile);
	static void begin_draw(shader_t &s);
	static void end_draw  (shader_t &s);
};

struct vect_butterfly_t : public animal_group_t<butterfly_t> {
	static void begin_draw(shader_t &s);
	static void end_draw  (shader_t &s);
};

bool birds_active();

