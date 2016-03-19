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
	bool is_enabled() const {return enabled;}
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
	void gen(unsigned num, cube_t const &range) { // Note: okay if nonempty
		reserve(size() + num); // optional

		for (unsigned n = 0; n < num; ++n) {
			A animal;
			if (animal.gen(rgen, range)) {push_back(animal);} // only add if generation was successful
		}
	}
	void update() {
		for (iterator i = begin(); i != end(); ++i) {i->update(rgen);}
	}
	void remove(unsigned ix) {
		assert(ix < size());
		std::swap(operator[](ix), back());
		pop_back();
	}
	void remove_disabled() {
		iterator i(begin()), o(i);
		for (; i != end(); ++i) {if (i->is_enabled()) {*(o++) = *i;}}
		erase(o, end());
	}
	void draw_animals(shader_t &s) const {
		for (const_iterator i = begin(); i != end(); ++i) {i->draw(s);}
	}
	void clear() {vector<A>::clear(); generated = 0;}
};

struct vect_fish_t : public animal_group_t<fish_t> {
	void draw() const;
};

struct vect_bird_t : public animal_group_t<bird_t> {
	void draw() const;
};


#endif // _ANIMALS_H_
