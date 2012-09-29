// 3D World - asteroid classes universe mode
// by Frank Gennari
// 9/26/12

#include "ship.h"
#include "voxels.h"
#include "shape_line3d.h" // for rock_shape3d
#include "upsurface.h"
#include "shaders.h"
#include "gl_ext_arb.h"


unsigned const ASTEROID_NDIV   = 32; // for sphere model, better if a power of 2
unsigned const ASTEROID_VOX_SZ = 64; // for voxel model
float    const AST_COLL_RAD    = 0.25;


extern vector<us_weapon> us_weapons;


class uobj_asteroid_sphere : public uobj_asteroid {
public:
	uobj_asteroid_sphere(point const &pos_, float radius_, unsigned lt) : uobj_asteroid(pos_, radius_, lt) {}
	virtual void draw_obj(uobj_draw_data &ddata) const {ddata.draw_asteroid();}
};


class uobj_asteroid_rock3d : public uobj_asteroid {

	rock_shape3d model3d;

public:
	uobj_asteroid_rock3d(point const &pos_, float radius_, unsigned lt, int type) : uobj_asteroid(pos_, radius_, lt) {
		static int obj_id(0); // for random seed
		model3d.gen_rock(24, 1.0, ++obj_id, type); // pos starts at and stays at all_zeros
		model3d.set_texture(MOON_TEX, 0.2);
	}
	~uobj_asteroid_rock3d() {model3d.destroy();}

	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(); return;}
		WHITE.do_glColor();
		model3d.draw();
		end_texture();
		enable_blend();
	}
};


class uobj_asteroid_destroyable : public uobj_asteroid {

public:
	uobj_asteroid_destroyable(point const &pos_, float radius_, unsigned lt) : uobj_asteroid(pos_, radius_, lt) {}
	virtual void apply_damage(float damage, point const &hit_pos) = 0;

	virtual float damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {
		// similar to add_damage_to_smiley_surface()
		bool gen_fragments(type == DAMAGE_EXP && val >= 20.0 && hit_pos != all_zeros && hit_pos != pos);
	
		if (gen_fragments && wc >= 0) {
			assert(unsigned(wc) < us_weapons.size());
			if (us_weapons[wc].const_dam) gen_fragments = 0;
		}
		if (gen_fragments) {
			float const damage(min(0.5, 0.02*val));
			apply_damage(damage, hit_pos);
			gen_moving_fragments(hit_pos, min(25U, max(1U, unsigned(20*damage))), 0.25);
		}
		return uobj_asteroid::damage(val, type, hit_pos, source, wc);
	}

	virtual bool has_detailed_coll(free_obj const *const other_obj) const {
		assert(other_obj);
		return (!other_obj->is_ship() && other_obj->get_radius() <= AST_COLL_RAD*radius);
	}

	virtual bool ship_int_obj(u_ship const *const ship,  intersect_params &ip=intersect_params()) const {
		assert(ship);
		if (ship->has_detailed_coll(this)) return uobj_asteroid::ship_int_obj(ship, ip);
		return sphere_int_obj(ship->get_pos(), ship->get_c_radius(), ip); // has simple cobjs, use mesh model
	}

	virtual bool obj_int_obj (free_obj const *const obj, intersect_params &ip=intersect_params()) const {
		assert(obj);
		if (obj->has_detailed_coll(this)) return uobj_asteroid::obj_int_obj(obj, ip);
		return sphere_int_obj(obj->get_pos(), obj->get_c_radius(), ip); // has simple cobjs, use mesh model
	}
};


class uobj_asteroid_hmap : public uobj_asteroid_destroyable {

	float scale_val;
	upsurface surface;
	static vector<float> pmap_vector;

public:
	uobj_asteroid_hmap(point const &pos_, float radius_, unsigned lt) : uobj_asteroid_destroyable(pos_, radius_, lt) {
		static int obj_id(0); // for random seed
		set_rand2_state(++obj_id, 1);
		surface.gen(0.15, 2.0, 10, 1.0);
		surface.setup(ASTEROID_NDIV, 0.0, 0);
		surface.setup_draw_sphere(all_zeros, 1.0, 0.0, ASTEROID_NDIV, NULL);
		surface.calc_rmax();
		scale_val = 1.0/surface.rmax;
	}

	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(); return;}
		uniform_scale(scale_val);
		WHITE.do_glColor();
		//select_texture(MOON_TEX);
		select_texture(ROCK_SPHERE_TEX);
		surface.sd.draw_ndiv_pow2(ddata.ndiv); // use dlist?
		end_texture();
	}

	virtual void apply_damage(float damage, point const &hit_pos) {
		int tx, ty;
		get_tex_coords_at(hit_pos, tx, ty);
		int const tsize(int(0.05*damage*ASTEROID_NDIV + 1.0)), radsq(4*tsize*tsize);
		int x1(tx - tsize), y1(ty - 2*tsize), x2(tx + tsize), y2(ty + 2*tsize);
		point **points = surface.sd.get_points();
		assert(points);

		for (int yy = y1; yy < y2; ++yy) { // allow texture wrap
			int const y((yy + ASTEROID_NDIV) % ASTEROID_NDIV), yterm((yy-ty)*(yy-ty));

			for (int xx = x1; xx < x2; ++xx) { // allow texture wrap
				int const x((xx + ASTEROID_NDIV) % ASTEROID_NDIV);
				if (((xx-tx)*(xx-tx) << 2) + yterm > radsq) continue;
				float const dist(sqrt(float((xx-tx)*(xx-tx) + yterm)));
				point &pt(points[x][y]);
				float const pt_mag(pt.mag());
				if (pt_mag <= 0.2) continue;
				pt -= pt*(damage/((dist + 2.0)*pt_mag));
				float const pt_mag2(pt.mag());
				if (pt_mag2 < 0.19) pt *= 0.19/pt_mag2;
			}
		}
	}

	virtual float const *get_sphere_shadow_pmap(point const &sun_pos, point const &obj_pos, int ndiv) const {
		assert(ndiv >= 3);
		float const dist_to_sun(p2p_dist(pos, sun_pos)), scale_val((dist_to_sun + p2p_dist(pos, obj_pos))/dist_to_sun);
		pmap_vector.resize(ndiv);
		point const ce[2] = {pos, sun_pos};
		vector3d v12; // unused
		vector_point_norm const &vpn(gen_cylinder_data(ce, c_radius, 0.0, ndiv, v12));

		for (int i = 0; i < ndiv; ++i) { // assumes the cylinder is more or less constant radius
			pmap_vector[i] = scale_val*(get_radius_at(vpn.p[i<<1]) - c_radius);
		}
		return &pmap_vector.front();
	}

	virtual bool sphere_int_obj(point const &c, float r, intersect_params &ip=intersect_params()) const {
		if (r > AST_COLL_RAD*radius) return uobj_asteroid_destroyable::sphere_int_obj(c, r, ip); // use default sphere collision
		float const asteroid_radius(get_radius_at(c));
		if (!dist_less_than(pos, c, (r + asteroid_radius))) return 0;
	
		if (ip.calc_int) {
			get_sphere_mov_sphere_int_pt(pos, c, (c - ip.p_last), 1.01*(r + asteroid_radius), ip.p_int);
			ip.norm = (ip.p_int - pos).get_norm();
			//ip.norm  = (c - pos).get_norm();
			//ip.p_int = pos + ip.norm*(1.01*(r + asteroid_radius));
		}
		return 1;
	}

	private:
	float get_radius_at(point const &pt) const {
		int tx, ty;
		get_tex_coords_at(pt, tx, ty);
		point **points = surface.sd.get_points();
		assert(points);
		return radius*scale_val*points[tx][ty].mag();
	}

	void get_tex_coords_at(point const &query_pos, int &tx, int &ty) const {
		vector3d query_dir(query_pos - pos);
		rotate_point(query_dir);
		query_dir.normalize();
		get_tex_coord(query_dir, vector3d(0.0, 1.0, 0.0), ASTEROID_NDIV, ASTEROID_NDIV, tx, ty, 1);
	}
};

vector<float> uobj_asteroid_hmap::pmap_vector; // static


class uobj_asteroid_voxel : public uobj_asteroid_destroyable {

	mutable voxel_model model; // FIXME: const problems

public:
	uobj_asteroid_voxel(point const &pos_, float radius_, unsigned lt) : uobj_asteroid_destroyable(pos_, radius_, lt) {
		static int obj_id(0); // for random seed
		RESET_TIME;
		gen_voxel_asteroid(model, all_zeros, 1.0, ASTEROID_VOX_SZ, ++obj_id); // will be translated to pos and scaled by radius during rendering
		model.build(0, 0, 0); // no cobjs
		PRINT_TIME("Create Asteroid");
	}
	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(); return;}
		if (ddata.shader.is_setup()) {ddata.shader.disable();}
		WHITE.do_glColor();
		shader_t s;
		// FIXME: write
		camera_pdu.valid = 0; // disable view frustum culling because it's not correct (due to transform matrices)
		model.core_render(s, 0);
		camera_pdu.valid = 1;
		//model.render(0);
		s.end_shader();
		if (ddata.shader.is_setup()) {ddata.shader.enable();}
	}
	virtual void apply_damage(float damage, point const &hit_pos) {
		// FIXME: write
	}
	virtual bool sphere_int_obj(point const &c, float r, intersect_params &ip=intersect_params()) const {
		if (r > AST_COLL_RAD*radius) return uobj_asteroid_destroyable::sphere_int_obj(c, r, ip); // use default sphere collision
		// FIXME: write
		return 1;
	}
	virtual bool line_int_obj(point const &p1, point const &p2, point *p_int=NULL, float *dscale=NULL) const {
		// FIXME: write
		return 1;
	}
	virtual void clear_context() {model.free_context();}
	// FIXME: shadow function?
};


uobj_asteroid *uobj_asteroid::create(point const &pos, float radius, unsigned model, unsigned lt) {

	switch (model) {
	case AS_MODEL_SPHERE: return new uobj_asteroid_sphere(pos, radius, lt);
	case AS_MODEL_ROCK1:  return new uobj_asteroid_rock3d(pos, radius, lt, 0);
	case AS_MODEL_ROCK2:  return new uobj_asteroid_rock3d(pos, radius, lt, 1);
	case AS_MODEL_HMAP:   return new uobj_asteroid_hmap  (pos, radius, lt);
	case AS_MODEL_VOXEL:  return new uobj_asteroid_voxel (pos, radius, lt);
	default: assert(0);
	}
	return NULL; // never gets here
}


void uobj_asteroid::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass,
							int align, unsigned eflags, free_obj const *parent_)
{
	gen_fragments();
	uobject::explode(damage, bradius, etype, edir, exp_time, wclass, align, eflags, parent_);
}


