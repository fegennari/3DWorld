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
	float orbital_dist; // for asteroid_belt asteroids

public:
	int last_coll_id;

	uasteroid() : inst_id(0), orbital_dist(0.0), last_coll_id(-1) {}
	void gen_base(float max_radius);
	void gen_spherical(upos_point_type const &pos_offset, float max_dist, float max_radius);
	void gen_belt(upos_point_type const &pos_offset, vector3d const &orbital_plane_normal, vector3d const vxy[2],
		float belt_radius, float belt_width, float belt_thickness, float max_radius, float &ri_max, float &plane_dmax);
	void apply_field_physics(point const &af_pos, float af_radius);
	void apply_belt_physics(upos_point_type const &af_pos, upos_point_type const &op_normal);
	void draw(point_d const &pos_, point const &camera, shader_t &s, pt_line_drawer &pld) const;
	void destroy();
	void set_velocity(vector3d const &v) {velocity = v;}
	unsigned get_rseed()           const {return inst_id;}
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
protected:
	vector<sphere_t> shadow_casters;
	sphere_t sun_pos_radius;

	virtual void gen_asteroid_placements() = 0;
	void remove_asteroid(unsigned ix);

public:
	uasteroid_cont() : rseed(0) {}
	void init(point const &pos, float radius);
	void gen_asteroids();
	void draw(point_d const &pos_, point const &camera, shader_t &s, bool sun_light_already_set);
	void detatch_asteroid(unsigned ix);
	void destroy_asteroid(unsigned ix);
	void free_uobj() {clear();}
	void begin_render(shader_t &shader, bool custom_lighting) {begin_render(shader, shadow_casters.size(), custom_lighting);}

	static void begin_render(shader_t &shader, unsigned num_shadow_casters, bool custom_lighting);
	static void end_render(shader_t &shader);
};


class uasteroid_field : public uasteroid_cont {

	vector<unsigned short> grid[AF_GRID_SZ][AF_GRID_SZ][AF_GRID_SZ];

public:
	void apply_physics(point_d const &pos_, point const &camera);
	virtual void gen_asteroid_placements();
};


class uasteroid_belt : public uasteroid_cont {

protected:
	vector3d orbital_plane_normal, scale;
	float max_asteroid_radius, inner_radius, outer_radius;

	void xform_to_local_torus_coord_space(point &pt) const;
	void gen_belt_placements(unsigned max_num, float belt_width, float belt_thickness, float max_ast_radius);

public:
	uasteroid_belt(vector3d const &opn, vector3d const &scale_) :
	  orbital_plane_normal(opn), inner_radius(0.0), outer_radius(0.0), max_asteroid_radius(0.0), scale(scale_) {}
	virtual void apply_physics(point_d const &pos_, point const &camera) = 0;
	bool line_might_intersect(point const &p1, point const &p2, float line_radius) const;
	bool sphere_might_intersect(point const &sc, float sr) const;
	float get_dist_to_boundary(point const &pt) const;
	float get_max_asteroid_radius() const {return max_asteroid_radius;}
	void draw_bounding_torus(point const &pos_, colorRGBA const &color=GREEN) const;
};


class uasteroid_belt_system : public uasteroid_belt {

	ussystem *system;

	virtual void gen_asteroid_placements();
	bool might_cast_shadow(uobject const &uobj) const;
	void calc_shadowers();

public:
	uasteroid_belt_system(vector3d const &opn, ussystem *system_) :
	  uasteroid_belt(opn, vector3d(1,1,1)), system(system_) {}
	virtual void apply_physics(upos_point_type const &pos_, point const &camera);
	void draw_detail(point_d const &pos_, point const &camera) const;
};


class uasteroid_belt_planet : public uasteroid_belt {

	float bwidth;
	uplanet *planet;

	virtual void gen_asteroid_placements();
	void calc_shadowers();

public:
	uasteroid_belt_planet(vector3d const &opn, uplanet *planet_) : uasteroid_belt(opn, planet_->rscale), bwidth(0.0), planet(planet_) {}
	void init_rings(point const &pos);
	virtual void apply_physics(upos_point_type const &pos_, point const &camera);
};


#endif // _ASTEROID_H_


