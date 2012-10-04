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
float    const AST_COLL_RAD    = 0.25; // limit collisions of large objects for accuracy (heightmap) or performance (voxel)


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
			float const damage_val(min(0.5, 0.02*val));
			apply_damage(damage_val, hit_pos);
			gen_moving_fragments(hit_pos, min(25U, max(1U, unsigned(20*damage_val))), get_fragment_tid(hit_pos), 0.25);
		}
		return uobj_asteroid::damage(val, type, hit_pos, source, wc);
	}
	virtual int get_fragment_tid(point const &hit_pos) const = 0;

	virtual bool has_detailed_coll(free_obj const *const other_obj) const {
		assert(other_obj);
		return (!other_obj->is_ship() && other_obj->get_radius() <= AST_COLL_RAD*radius);
	}
	bool check_sphere_int(free_obj const *const obj, intersect_params &ip) const {
		assert(obj);
		return sphere_int_obj(obj->get_pos(), obj->get_c_radius(), ip);
	}
	virtual bool ship_int_obj(u_ship const *const ship,  intersect_params &ip=intersect_params()) const {
		if (!check_sphere_int(ship, ip))    return 0; // use mesh model
		if (!ship->has_detailed_coll(this)) return 1; // simple intersection
		return uobj_asteroid::ship_int_obj(ship, ip); // has detailed cobjs, do detailed intersection
	}
	virtual bool obj_int_obj (free_obj const *const obj, intersect_params &ip=intersect_params()) const {
		if (!check_sphere_int(obj, ip))    return 0; // use mesh model
		if (!obj->has_detailed_coll(this)) return 1; // simple intersection
		return uobj_asteroid::obj_int_obj(obj, ip);  // has detailed cobjs, do detailed intersection
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
	virtual int get_fragment_tid(point const &hit_pos) const {return ROCK_SPHERE_TEX;}

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


// FIXME:
// sphere collision normal
// AO lighting
class uobj_asteroid_voxel : public uobj_asteroid_destroyable {

	mutable voxel_model_space model; // FIXME: const problems
	bool have_sun_pos;

public:
	uobj_asteroid_voxel(point const &pos_, float radius_, unsigned lt) : uobj_asteroid_destroyable(pos_, radius_, lt), have_sun_pos(0) {
		static int obj_id(0); // for random seed
		RESET_TIME;
		gen_voxel_asteroid(model, all_zeros, 1.0, ASTEROID_VOX_SZ, ++obj_id); // will be translated to pos and scaled by radius during rendering
		model.build(0, 0);
		float const gen_radius(model.get_bsphere().radius);
		assert(gen_radius > 0.0);
		radius /= gen_radius;
		PRINT_TIME("Create Asteroid");
	}

	virtual void apply_physics() {
		point sun_pos;

		if (!have_sun_pos && get_universe_sun_pos(pos, sun_pos)) { // too slow if there is no sun?
			dir = (pos - sun_pos).get_norm(); // orient toward the sun
			have_sun_pos = 1;
		}
		uobj_asteroid_destroyable::apply_physics();
		model.proc_pending_updates();
	}

	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(); return;}
		if (ddata.shader.is_setup()) {ddata.shader.disable();}
		
		shader_t s;
		s.set_int_prefix("num_lights", min(8U, exp_lights.size()+2U), 1); // FS
		s.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
		if (!glIsEnabled(GL_FOG)) s.set_prefix("#define NO_FOG", 1); // FS
		s.set_vert_shader("asteroid");
		s.set_frag_shader("linear_fog.part+ads_lighting.part*+triplanar_texture.part+procedural_texture.part+voxel_texture.part+asteroid");
		s.begin_shader();
		s.setup_fog_scale();
		s.add_uniform_int("tex0", 0);
		s.add_uniform_int("tex1", 8);
		s.add_uniform_int("noise_tex", 5);
		s.add_uniform_float("tex_scale", 1.0);
		s.add_uniform_float("noise_scale", 0.1);
		s.add_uniform_float("tex_mix_saturate", 5.0);

		model.setup_tex_gen_for_rendering(s);
		WHITE.do_glColor();
		set_specular(0.0, 1.0); // ???
		glEnable(GL_CULL_FACE);
		camera_pdu.valid = 0; // disable view frustum culling because it's not correct (due to transform matrices)
		model.core_render(s, 0);
		camera_pdu.valid = 1;
		glDisable(GL_CULL_FACE);
		s.end_shader();
		if (ddata.shader.is_setup()) {ddata.shader.enable();}
		select_texture(WHITE_TEX, 0);
	}

	virtual void apply_damage(float damage, point const &hit_pos) {
		float const damage_radius(min(0.2, 0.1*damage));
		point center(hit_pos);
		xform_point(center);
		model.update_voxel_sphere_region(center, damage_radius, -1.0);
	}

	virtual int get_fragment_tid(point const &hit_pos) const {
		point p(hit_pos);
		xform_point(p);
		return model.get_texture_at(p);
	}

	virtual bool sphere_int_obj(point const &c, float r, intersect_params &ip=intersect_params()) const {
		if (r > AST_COLL_RAD*radius) return uobj_asteroid_destroyable::sphere_int_obj(c, r, ip); // use default sphere collision
		r /= radius; // scale to 1.0
		point p(c);
		xform_point(p);
		if (!model.sphere_intersect(p, r, (ip.calc_int ? &ip.p_int : NULL))) return 0;

		if (ip.calc_int) {
			ip.norm = (p - pos).get_norm(); // we can't actually calculate the normal, so we use the direction from asteroid center to object center
			if (ip.p_int == p) {ip.p_int = p - ip.norm*r;} // intersecting at the center, determine actual pos based on normal
			xform_point_inv(ip.p_int);
			rotate_point_inv(ip.norm); // ip.norm will be normalized
		}
		return 1;
	}
	virtual bool line_int_obj(point const &p1, point const &p2, point *p_int=NULL, float *dscale=NULL) const { // Note: dscale is ignored
		point p[2] = {p1, p2};
		xform_point_x2(p[0], p[1]);
		if (!model.line_intersect(p[0], p[1], p_int)) return 0;
		if (p_int) {xform_point_inv(*p_int);}
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


