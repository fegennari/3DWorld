// 3D World - Universe Asteroid Header
// by Frank Gennari
// 12/15/12
#pragma once

#include "universe.h"

unsigned const AF_GRID_SZ = 12;


struct asteroid_belt_cloud : public volume_part_cloud {

	point pos;
	float radius;
	unsigned vbo_pos;

	asteroid_belt_cloud() : radius(0.0f), vbo_pos(0) {}
	void gen(rand_gen_t &rgen, float def_radius);
	static void pre_draw(vpc_shader_t &s, colorRGBA const &color, float noise_scale);
	static void post_draw(vpc_shader_t &s);
	void draw(vpc_shader_t &s, point_d const &pos_, float def_cloud_radius, float shadow_atten) const;
};


class uasteroid : public uobject, public rotated_obj {

	unsigned inst_id;
	vector3d scale, velocity;
	float orbital_dist; // for asteroid_belt asteroids

public:
	int last_coll_id;
	bool is_ice;

	uasteroid() : inst_id(0), orbital_dist(0.0), last_coll_id(-1), is_ice(0) {}
	void gen_base(float max_radius);
	void gen_spherical(upos_point_type const &pos_offset, float max_dist, float max_radius);
	void gen_belt(upos_point_type const &pos_offset, vector3d const &orbital_plane_normal, vector3d const vxy[2],
		float belt_radius, float belt_width, float belt_thickness, float max_radius, float &ri_max, float &plane_dmax);
	void apply_field_physics(point const &af_pos, float af_radius);
	void apply_belt_physics(upos_point_type const &af_pos, upos_point_type const &op_normal, vector3d const &orbit_scale, vector<sphere_t> const &colliders);
	void draw(point_d const &pos_, point const &camera, shader_t &s, pt_line_drawer &pld) const;
	void destroy();
	void set_velocity(vector3d const &v) {velocity = v;}
	unsigned get_rseed()           const {return inst_id;}
	vector3d const &get_scale()    const {return scale;}
	vector3d const &get_velocity() const {return velocity;}
	float get_rel_mass()           const {return scale.x*scale.y*scale.z*radius*radius*radius;} // mass is proportional to volume which is proportional to radius^3
	bool operator<(uasteroid const &a) const {return (inst_id < a.inst_id);} // for sorting by inst_id
	bool line_intersection(point const &p1, vector3d const &v12, float line_length, float line_radius, float &ldist) const;
	virtual bool sphere_intersection(point const &c, float r) const;

	virtual std::string get_name() const {return (is_ice ? "Ice Fragment" : "Asteroid");}
	virtual bool rename(std::string const &name_) {return 0;} // not renameable
	virtual int  get_owner() const {return NO_OWNER;}
	virtual int  get_fragment_tid(point const &hit_pos) const;
};


class shadowed_uobject {
protected:
	vector<sphere_t> shadow_casters;
	sphere_t sun_pos_radius;

public:
	void calc_shadowers_for_planet(uplanet const &planet);
	void upload_shadow_casters(shader_t &s) const;
};


class uasteroid_cont : public uobject_base, public shadowed_uobject, public vector<uasteroid> {

	int rseed;
protected:
	pt_line_drawer pld; // for drawing

	virtual void gen_asteroid_placements() = 0;
	virtual void remove_asteroid(unsigned ix);

public:
	uasteroid_cont() : rseed(0) {}
	virtual ~uasteroid_cont() {}
	void init(point const &pos, float radius);
	virtual bool get_is_ice() const {return 0;}
	virtual void gen_asteroids();
	void draw(point_d const &pos_, point const &camera, shader_t &s, bool sun_light_already_set);
	void detach_asteroid(unsigned ix);
	void destroy_asteroid(unsigned ix);
	void free_uobj() {clear();}
	void begin_render(shader_t &shader, bool custom_lighting) {begin_render(shader, shadow_casters.size(), custom_lighting);}
	float calc_shadow_atten(point const &cpos) const;

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
	struct cloud_inst {
		unsigned asteroid_id, cloud_id;
		cloud_inst(unsigned aid=0, unsigned cid=0) : asteroid_id(aid), cloud_id(cid) {}
	};
	struct cloud_dist_cmp {
		bool operator()(pair<float, cloud_inst> const &a, pair<float, cloud_inst> const &b) const {return (a.first < b.first);}
	};
	vector3d orbital_plane_normal, orbit_scale;
	float max_asteroid_radius, inner_radius, outer_radius, temperature;
	vector<cloud_inst> cloud_insts;
	mutable vector<pair<float, cloud_inst>> clouds_to_draw;

	void xform_to_local_torus_coord_space(point &pt) const;
	void xform_from_local_torus_coord_space(point &pt) const;
	void gen_belt_placements(unsigned max_num, float belt_width, float belt_thickness, float max_ast_radius);

public:
	uasteroid_belt(vector3d const &opn, vector3d const &scale_) :
	orbital_plane_normal(opn), orbit_scale(scale_), max_asteroid_radius(0.0), inner_radius(0.0), outer_radius(0.0), temperature(0.0) {}
	virtual bool is_planet_ab() const {return 0;}
	virtual void gen_asteroids();
	virtual void apply_physics(point_d const &pos_, point const &camera) = 0;
	virtual void remove_asteroid(unsigned ix);
	bool line_might_intersect(point const &p1, point const &p2, float line_radius, point *p_int=nullptr) const;
	bool sphere_might_intersect(point const &sc, float sr) const;
	float get_line_sphere_int_radius_scale() const;
	float get_dist_to_boundary(point const &pt) const;
	float get_max_asteroid_radius() const {return max_asteroid_radius;}
	void draw_detail(point_d const &pos_, point const &camera, bool no_asteroid_dust, bool draw_dust, float density) const;
};


class uasteroid_belt_system : public uasteroid_belt {

	ussystem *system;
	vector<sphere_t> colliders;

	virtual void gen_asteroid_placements();
	void add_potential_collider(point const &cpos, float cradius);
	bool might_cast_shadow(uobject const &uobj) const;
	void calc_colliders();
	void calc_shadowers();

public:
	uasteroid_belt_system(vector3d const &opn, ussystem *system_) : uasteroid_belt(opn, system_->orbit_scale), system(system_) {}
	virtual bool get_is_ice() const {return (temperature < 6.0);} // 50% of FREEZE_TEMP
	virtual void apply_physics(upos_point_type const &pos_, point const &camera);
};


class uasteroid_belt_planet : public uasteroid_belt {

	float bwidth;
	uplanet *planet;

	virtual void gen_asteroid_placements();
	void calc_shadowers();

public:
	uasteroid_belt_planet(vector3d const &opn, uplanet *planet_) : uasteroid_belt(opn, planet_->rscale), bwidth(0.0), planet(planet_) {}
	virtual bool is_planet_ab() const {return 1;}
	virtual bool get_is_ice  () const {return planet->has_ice_debris();}
	void init_rings(point const &pos);
	virtual void apply_physics(upos_point_type const &pos_, point const &camera);
};

