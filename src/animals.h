// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#ifndef _ANIMALS_H_
#define _ANIMALS_H_

#include "3DWorld.h"


class animal_t : public sphere_t {

protected:
	bool enabled;
	vector3d dir;
	vector3d velocity;
	colorRGBA color;

	int get_ndiv() const;

public:
	animal_t() : enabled(0) {}
	bool is_visible() const;
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


class animal_group_t {

	rand_gen_t rgen;

public:
	template<typename A> void update_animals(vector<A> &animals) {
		for (auto i = animals.begin(); i != animals.end(); ++i) {i->update(rgen);}
	}
	template<typename A> void draw_animals(vector<A> const &animals, shader_t &s) const {
		for (auto i = animals.begin(); i != animals.end(); ++i) {i->draw(s);}
	}
};

struct vect_fish_t : public vector<fish_t>, public animal_group_t {
	void update() {update_animals(*this);}
	void draw() const;
};

struct vect_bird_t : public vector<bird_t>, public animal_group_t {
	void update() {update_animals(*this);}
	void draw() const;
};


#endif // _ANIMALS_H_
