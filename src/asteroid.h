// 3D World - Universe Asteroid Header
// by Frank Gennari
// 12/15/12

#ifndef _ASTEROID_H_
#define _ASTEROID_H_

#include "universe.h"

unsigned const AF_GRID_SZ = 8;


class uasteroid : public uobject, public rotated_obj {
	
	unsigned inst_id;
	vector3d scale, velocity;

public:
	int last_coll_id;

	uasteroid() : inst_id(0), last_coll_id(-1) {}
	void gen(upos_point_type const &pos_offset, float max_dist, float max_radius);
	void apply_physics(point const &af_pos, float af_radius);
	void draw(point_d const &pos_, point const &camera, shader_t &s) const;
	void destroy();
	void set_velocity(vector3d const &v) {velocity = v;}
	vector3d const &get_scale()    const {return scale;}
	vector3d const &get_velocity() const {return velocity;}
	float get_rel_mass()           const {return scale.x*scale.y*scale.z*radius*radius*radius;} // mass is proportional to volume which is proportional to radius^3

	virtual std::string get_name() const {return "Asteroid";}
	virtual bool rename(std::string const &name_) {return 0;} // not renaemable
	virtual int  get_owner() const {return NO_OWNER;}
	virtual int  get_fragment_tid(point const &hit_pos) const;
};


class uasteroid_field : public uobject_base, public vector<uasteroid> {

	vector<unsigned short> grid[AF_GRID_SZ][AF_GRID_SZ][AF_GRID_SZ];
	int rseed;

public:
	uasteroid_field() : rseed(0) {}
	void init(point const &pos, float radius);
	static void begin_render(shader_t &shader);
	static void end_render(shader_t &shader);
	void gen_asteroids();
	void apply_physics(point_d const &pos_, point const &camera);
	void draw(point_d const &pos_, point const &camera, shader_t &s);
	void destroy_asteroid(unsigned ix);
	void free_uobj() {clear();}
};


#endif // _ASTEROID_H_


