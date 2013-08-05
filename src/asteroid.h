// 3D World - Universe Asteroid Header
// by Frank Gennari
// 12/15/12

#ifndef _ASTEROID_H_
#define _ASTEROID_H_

#include "universe.h"

unsigned const AF_GRID_SZ = 12;


class uasteroid : public uobject, public rotated_obj {

	unsigned inst_id;
	vector3d scale, velocity;

public:
	int last_coll_id;

	uasteroid() : inst_id(0), last_coll_id(-1) {}
	void gen_base(float max_radius);
	void gen_spherical(upos_point_type const &pos_offset, float max_dist, float max_radius);
	float gen_belt(upos_point_type const &pos_offset, vector3d const &orbital_plane_normal,
		float belt_radius, float belt_width, float belt_thickness, float max_radius);
	void apply_field_physics(point const &af_pos, float af_radius);
	void apply_belt_physics(upos_point_type const &af_pos, vector3d const &op_normal, float af_radius);
	void draw(point_d const &pos_, point const &camera, shader_t &s, pt_line_drawer &pld) const;
	void destroy();
	void set_velocity(vector3d const &v) {velocity = v;}
	vector3d const &get_scale()    const {return scale;}
	vector3d const &get_velocity() const {return velocity;}
	float get_rel_mass()           const {return scale.x*scale.y*scale.z*radius*radius*radius;} // mass is proportional to volume which is proportional to radius^3
	bool operator<(uasteroid const &a) const {return (inst_id < a.inst_id);} // for sorting by inst_id
	bool line_intersection(point const &p1, vector3d const &v12, float line_length, float line_radius, float &ldist) const;

	virtual std::string get_name() const {return "Asteroid";}
	virtual bool rename(std::string const &name_) {return 0;} // not renaemable
	virtual int  get_owner() const {return NO_OWNER;}
	virtual int  get_fragment_tid(point const &hit_pos) const;
};


class uasteroid_cont : public uobject_base, public vector<uasteroid> {

	int rseed;

public:
	uasteroid_cont() : rseed(0) {}
	void init(point const &pos, float radius);
	void gen_asteroids();
	void draw(point_d const &pos_, point const &camera, shader_t &s);
	void destroy_asteroid(unsigned ix);
	void free_uobj() {clear();}
	virtual void gen_asteroid_placements() = 0;

	static void begin_render(shader_t &shader);
	static void end_render(shader_t &shader);
};


class uasteroid_field : public uasteroid_cont {

	vector<unsigned short> grid[AF_GRID_SZ][AF_GRID_SZ][AF_GRID_SZ];

public:
	void apply_physics(point_d const &pos_, point const &camera);
	virtual void gen_asteroid_placements();
};


class uasteroid_belt : public uasteroid_cont {

	vector3d orbital_plane_normal;
	float inner_radius, outer_radius, max_asteroid_radius;

public:
	uasteroid_belt(vector3d const &opn) : orbital_plane_normal(opn), inner_radius(0.0), outer_radius(0.0), max_asteroid_radius(0.0) {}
	void apply_physics(point_d const &pos_, point const &camera);
	virtual void gen_asteroid_placements();
	bool line_might_intersect(point const &p1, point const &p2) const;
	bool sphere_might_intersect(point const &sc, float sr) const;
	float get_max_asteroid_radius() const {return max_asteroid_radius;}
};


#endif // _ASTEROID_H_


