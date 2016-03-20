// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#ifndef _ANIMALS_H_
#define _ANIMALS_H_

#include "3DWorld.h"
#include "function_registry.h"
#include "inlines.h"


struct tile_offset_t {
	int dxoff, dyoff;

	tile_offset_t(int dxoff_=0, int dyoff_=0) : dxoff(dxoff_), dyoff(dyoff_) {}
	void set_from_xyoff2() {dxoff = -xoff2; dyoff = -yoff2;}
	vector3d get_xlate() const {return vector3d(((xoff - xoff2) - dxoff)*DX_VAL, ((yoff - yoff2) - dyoff)*DY_VAL, 0.0);}
	vector3d subtract_from(tile_offset_t const &o) const {return vector3d((o.dxoff - dxoff)*DX_VAL, (o.dyoff - dyoff)*DY_VAL, 0.0);}
};


class animal_t : public sphere_t {

protected:
	bool enabled;
	vector3d dir;
	vector3d velocity;
	colorRGBA color;
	//tile_offset_t tile_off;

	int get_ndiv(point const &pos_) const;

public:
	animal_t() : enabled(0) {}
	bool is_enabled() const {return enabled;}
	bool is_visible(point const &pos_) const;
	point get_draw_pos() const;
};

class fish_t : public animal_t {

	vector3d scale; // x = width, y = height, z = length

public:
	bool gen(rand_gen_t &rgen, cube_t const &range);
	bool update(rand_gen_t &rgen);
	void draw(shader_t &s) const;
};

class bird_t : public animal_t {

public:
	bool gen(rand_gen_t &rgen, cube_t const &range);
	bool update(rand_gen_t &rgen);
	void draw(shader_t &s) const;
};


class animal_group_base_t {
protected:
	rand_gen_t rgen;
	bool generated;

public:
	animal_group_base_t() : generated(0) {}
	bool was_generated() const {return generated;}
	static void begin_draw(shader_t &s);
	static void end_draw(shader_t &s);
};

template<typename A> class animal_group_t : public vector<A>, public animal_group_base_t {
public:
	void gen(unsigned num, cube_t const &range);
	void update();
	void remove(unsigned ix);
	void remove_disabled();
	void draw_animals(shader_t &s) const;
	void clear() {vector<A>::clear(); generated = 0;}
};

struct vect_fish_t : public animal_group_t<fish_t> {
	void draw() const;
};

struct vect_bird_t : public animal_group_t<bird_t> {
	void draw() const;
};


#endif // _ANIMALS_H_
