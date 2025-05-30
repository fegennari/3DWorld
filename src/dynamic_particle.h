// 3D World Dynamic Particle class declaration
// by Frank Gennari
// 7/17/06
#pragma once

#include "cube_map_shadow_manager.h"


struct dpart_params_t {
	float rmin, rmax, vmin, vmax, imin, imax; // {min, max}x{radius, velocity, intensity}
	vector3d sdist[2]; // low, high spawn distances

	dpart_params_t() : rmin(0.03), rmax(0.07), vmin(0.6), vmax(3.0), imin(0.12), imax(0.26) {
		sdist[0] = sdist[1] = vector3d(1,1,1);
	}
};


class dynamic_particle : public sphere_t { // size = 68

	bool moves, lighted, collides, chdir, gravity, shadows_setup;
	int tid, cid;
	float intensity, bwidth;
	vector3d velocity;
	colorRGBA color;

public:
	dynamic_particle();
	~dynamic_particle() {remove_cobj();}
	void gen_pos();
	void draw() const; // lights, color, texture, shadowed
	void apply_physics(float stepsize, int index); // begin_motion, move, random dir change, collision (mesh and cobjs), forces applied to?
	void add_light(cube_map_shadow_manager &smgr, int index); // dynamic lights
	void add_cobj_shadows() const; // cobjs, dynamic objects
	void add_cobj();
	void remove_cobj();
};


class dynamic_particle_system : public cube_map_shadow_manager { // size = 16

	vector<dynamic_particle>  particles;
	vector<vector<unsigned> > bins;
	int valid_frame;

public:
	dynamic_particle_system() : valid_frame(-1) {}
	void add_particle(dynamic_particle const &p) {particles.push_back(p);}
	size_t size() const {return particles.size();}
	void clear() {remove_lights(); particles.clear();}
	void create_particles(unsigned num, bool only_if_empty);
	void draw() const;
	void apply_physics(float stepsize=1.0);
	void add_lights();
	void remove_lights();
	void add_cobj_shadows() const;
	void build_lookup_matrix();
};


extern dynamic_particle_system d_part_sys;

