// 3D World - asteroid classes universe mode
// by Frank Gennari
// 9/26/12

#include "ship.h"
#include "voxels.h"
#include "shape_line3d.h" // for rock_shape3d
#include "upsurface.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "explosion.h"
#include "asteroid.h"


unsigned const ASTEROID_NDIV   = 32; // for sphere model, better if a power of 2
unsigned const ASTEROID_VOX_SZ = 64; // for voxel model
float    const AST_COLL_RAD    = 0.25; // limit collisions of large objects for accuracy (heightmap)
float    const AST_PROC_HEIGHT = 0.1; // height values of procedural shader asteroids


extern vector<us_weapon> us_weapons;

shader_t cached_voxel_shaders[8]; // one for each value of num_lights
shader_t cached_proc_shaders [8];


void clear_cached_shaders() {

	for (unsigned i = 0; i < 8; ++i) {
		cached_voxel_shaders[i].end_shader();
		cached_proc_shaders [i].end_shader();
	}
}


class uobj_asteroid_sphere : public uobj_asteroid {
	int tex_id;

public:
	uobj_asteroid_sphere(point const &pos_, float radius_, int tex_id_, unsigned lt) : uobj_asteroid(pos_, radius_, lt), tex_id(tex_id_) {}
	virtual void draw_obj(uobj_draw_data &ddata) const {ddata.draw_asteroid(tex_id);}
	virtual int get_fragment_tid(point const &hit_pos) const {return tex_id;}
};


class uobj_asteroid_rock3d : public uobj_asteroid {

	rock_shape3d model3d;

public:
	uobj_asteroid_rock3d(point const &pos_, float radius_, unsigned rseed_ix, unsigned lt, int type) : uobj_asteroid(pos_, radius_, lt) {
		model3d.gen_rock(24, 1.0, rseed_ix, type); // pos starts at and stays at all_zeros
		model3d.set_texture(MOON_TEX, 0.2);
	}
	~uobj_asteroid_rock3d() {model3d.destroy();}

	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(MOON_TEX); return;}
		WHITE.do_glColor();
		model3d.draw();
		end_texture();
		enable_blend();
	}
};


class uobj_asteroid_shader : public uobj_asteroid {
	int rseed_ix;

public:
	uobj_asteroid_shader(point const &pos_, float radius_, int rseed_ix_, unsigned lt) : uobj_asteroid(pos_, radius_, lt) {
		c_radius = (1.0 + AST_PROC_HEIGHT)*radius;
	}
	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.shader.is_setup()) {ddata.shader.disable();}
		unsigned const num_lights(min(8U, exp_lights.size()+2U));
		shader_t &s(cached_proc_shaders[num_lights]);
		
		if (s.is_setup()) { // already setup
			s.enable();
		}
		else {
			bind_3d_texture(get_noise_tex_3d(64, 1)); // grayscale noise
			s.set_int_prefix("num_lights", num_lights, 1); // FS
			s.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
			s.set_prefix("#define NO_SPECULAR",      1); // FS (optional/optimization)
			s.set_prefix("#define NUM_OCTAVES 8",    0); // VS
			s.set_vert_shader("perlin_clouds_3d.part*+procedural_rock");
			s.set_frag_shader("linear_fog.part+ads_lighting.part*+procedural_rock");
			s.begin_shader();
			s.setup_fog_scale();
			s.add_uniform_int("cloud_noise_tex", 0);
			s.add_uniform_float("time", float(rseed_ix));
			s.add_uniform_float("noise_scale",  0.1);
			s.add_uniform_float("height_scale", AST_PROC_HEIGHT);
		}
		colorRGBA(0.5, 0.45, 0.4, 1.0).do_glColor();
		select_texture(WHITE_TEX);
		draw_sphere_dlist(all_zeros, 1.0, 3*ddata.ndiv/2, 1);
		s.disable();

		if (ddata.final_pass) {
			if (ddata.shader.is_setup()) {ddata.shader.enable();}
			end_texture();
		}
	}
};


class uobj_asteroid_destroyable : public uobj_asteroid {

public:
	uobj_asteroid_destroyable(point const &pos_, float radius_, unsigned lt) : uobj_asteroid(pos_, radius_, lt) {}
	virtual bool apply_damage(float damage, point &hit_pos) = 0;

	virtual float damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {
		// similar to add_damage_to_smiley_surface()
		bool gen_fragments(type == DAMAGE_EXP && val >= 20.0 && hit_pos != all_zeros && hit_pos != pos);
	
		if (gen_fragments && wc >= 0) {
			assert(unsigned(wc) < us_weapons.size());
			if (us_weapons[wc].const_dam) gen_fragments = 0;
		}
		if (gen_fragments) {
			float const damage_val(min(0.5, 0.02*val));
			point mod_hit_pos(hit_pos);
			
			if (apply_damage(((wc == WCLASS_EXPLODE) ? 10.0 : 1.0)*damage_val, mod_hit_pos)) { // ship explosions are more damaging
				gen_moving_fragments(mod_hit_pos, min(25U, max(1U, unsigned(20*damage_val))), get_fragment_tid(mod_hit_pos), 0.25);
			}
		}
		return uobj_asteroid::damage(val, type, hit_pos, source, wc);
	}

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
	uobj_asteroid_hmap(point const &pos_, float radius_, unsigned rseed_ix, unsigned lt) : uobj_asteroid_destroyable(pos_, radius_, lt) {
		surface.rgen.set_state(rseed_ix, 1);
		surface.gen(0.15, 2.0, 10, 1.0);
		surface.setup(ASTEROID_NDIV, 0.0, 0);
		surface.setup_draw_sphere(all_zeros, 1.0, 0.0, ASTEROID_NDIV, NULL);
		surface.calc_rmax();
		scale_val = 1.0/surface.rmax;
	}

	virtual void draw_obj(uobj_draw_data &ddata) const {
		int const tex_id(ROCK_SPHERE_TEX); // MOON_TEX?
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(tex_id); return;}
		if (scale_val != 1.0) {uniform_scale(scale_val);}
		WHITE.do_glColor();
		select_texture(tex_id);
		surface.sd.draw_ndiv_pow2(ddata.ndiv); // use dlist?
		end_texture();
	}

	virtual bool apply_damage(float damage, point &hit_pos) {
		int tx, ty;
		get_tex_coords_at(hit_pos, tx, ty);
		int const tsize(int(0.05*damage*ASTEROID_NDIV + 1.0)), radsq(4*tsize*tsize);
		int x1(tx - tsize), y1(ty - 2*tsize), x2(tx + tsize), y2(ty + 2*tsize);
		point **points = surface.sd.get_points();
		assert(points);
		bool damaged(0);

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
				damaged = 1;
			}
		}
		return damaged;
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


// FIXME: sphere collision normal / pos
class uobj_asteroid_voxel : public uobj_asteroid_destroyable {

	mutable voxel_model_space model; // FIXME: const problems with draw()
	bool have_sun_pos;

public:
	uobj_asteroid_voxel(point const &pos_, float radius_, unsigned rseed_ix, unsigned lt) : uobj_asteroid_destroyable(pos_, radius_, lt), have_sun_pos(0) {
		//RESET_TIME;
		float const gen_radius(gen_voxel_rock(model, all_zeros, 1.0, ASTEROID_VOX_SZ, 2, rseed_ix)); // ndiv=2x2, will be translated to pos and scaled by radius during rendering
		assert(gen_radius > 0.0);
		radius /= gen_radius;
		//PRINT_TIME("Create Asteroid");
	}

	virtual void first_frame_hook() {
		point sun_pos;

		if (get_universe_sun_pos(pos, sun_pos)) { // too slow if there is no sun?
			dir = (sun_pos - pos).get_norm(); // orient toward the sun
			force_calc_rotation_vectors();
		}
	}

	virtual void apply_physics() {
		uobj_asteroid_destroyable::apply_physics();
		model.proc_pending_updates();

		if (!model.has_triangles()) { // completely destroyed (center anchor point is gone)
			explode(0.0, radius, ETYPE_NONE, zero_vector, 0, WCLASS_EXPLODE, ALIGN_NEUTRAL, 0, NULL);
		}
	}

	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(model.get_params().tids[0]); return;}
		if (ddata.shader.is_setup()) {ddata.shader.disable();}
		unsigned const num_lights(min(8U, exp_lights.size()+2U));
		shader_t &s(cached_voxel_shaders[num_lights]);
		
		if (s.is_setup()) { // already setup
			s.enable();
		}
		else {
			s.set_int_prefix("num_lights", num_lights, 1); // FS
			s.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
			s.set_prefix("#define NO_SPECULAR", 1); // FS (optional/optimization)
			s.set_vert_shader("asteroid");
			s.set_frag_shader("linear_fog.part+ads_lighting.part*+triplanar_texture.part+procedural_texture.part+voxel_texture.part+asteroid");
			s.begin_shader();
			s.setup_fog_scale();
			s.add_uniform_int("tex0", 0);
			s.add_uniform_int("tex1", 8);
			s.add_uniform_int("noise_tex", 5);
			s.add_uniform_int("ao_tex", 9);
			s.add_uniform_int("shadow_tex", 10);
			s.add_uniform_float("tex_scale", 1.0);
			s.add_uniform_float("noise_scale", 0.1);
			s.add_uniform_float("tex_mix_saturate", 5.0);
			s.add_uniform_vector3d("tex_eval_offset", zero_vector);
		}
		if (ddata.first_pass) {model.setup_tex_gen_for_rendering(s);}
		WHITE.do_glColor();
		glEnable(GL_CULL_FACE);
		model.core_render(s, 0, 1); // disable view frustum culling because it's incorrect (due to transform matrices)
		glDisable(GL_CULL_FACE);
		s.disable();

		if (ddata.final_pass) {
			if (ddata.shader.is_setup()) {ddata.shader.enable();}
			end_texture();
		}
	}

	virtual bool apply_damage(float damage, point &hit_pos) {
		float const damage_radius(min(0.5, 0.1*damage));
		xform_point(hit_pos);
		point const center(hit_pos);
		bool const damaged(model.update_voxel_sphere_region(center, damage_radius, -1.0, &hit_pos));
		xform_point_inv(hit_pos);
		hit_pos += 0.5*(radius/ASTEROID_VOX_SZ)*(hit_pos - pos).get_norm(); // move slightly away from asteroid center
		return damaged;
	}

	virtual int get_fragment_tid(point const &hit_pos) const {
		point p(hit_pos);
		xform_point(p);
		return model.get_texture_at(p);
	}

	virtual bool ship_int_obj(u_ship const *const ship, intersect_params &ip=intersect_params()) const {
		if (!uobj_asteroid_destroyable::ship_int_obj(ship, ip)) return 0;
		if (!ship->has_detailed_coll(this)) return 1; // simple intersection
		assert(ship);
		cobj_vector_t const &cobjs(ship->get_cobjs());
		assert(!cobjs.empty());
		point center(ship->get_pos()), p_last(ip.p_last);
		xform_point(center); // global to current local
		ship->xform_point(p_last); // global to ship local
		float const obj_radius(ship->get_radius()), sphere_radius(obj_radius/radius), sr(sphere_radius/ASTEROID_VOX_SZ);
		cube_t bcube;
		bcube.set_from_sphere(center, sphere_radius);
		int llc[3], urc[3];
		model.get_xyz(bcube.get_llc(), llc);
		model.get_xyz(bcube.get_urc(), urc);

		for (int y = max(0, llc[1]); y <= min(int(model.ny)-1, urc[1]); ++y) {
			for (int x = max(0, llc[0]); x <= min(int(model.nx)-1, urc[0]); ++x) {
				for (int z = max(0, llc[2]); z <= min(int(model.nz)-1, urc[2]); ++z) {
					point p(model.get_pt_at(x, y, z));
					if (!dist_less_than(p, center, sphere_radius) || model.is_outside(model.get_ix(x, y, z))) continue;
					xform_point_inv(p); // local to global
					ship->xform_point(p); // global to ship local

					for (cobj_vector_t::const_iterator c = cobjs.begin(); c != cobjs.end(); ++c) {
						assert(*c);
						
						if ((*c)->sphere_intersect(p, sr, p_last, ip.p_int, ip.norm, ip.calc_int)) {
							if (ip.calc_int) {
								ship->xform_point_inv(ip.p_int);
								ship->rotate_point_inv(ip.norm);
							}
							return 1;
						}
					} // for c
				} // for z
			} // for x
		} // for y
		return 0;
	}

	virtual bool sphere_int_obj(point const &c, float r, intersect_params &ip=intersect_params()) const { // Note: no size check
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

	virtual void draw_shadow_volumes(point const &targ_pos, float cur_radius, point const &sun_pos, int ndiv, bool test) const {
		ushadow_triangle_mesh(model.get_shadow_edge_tris(), targ_pos, cur_radius, sun_pos, radius, pos, this).draw_geom(targ_pos, test);
	}
	virtual bool casts_detailed_shadow() const {return !model.get_shadow_edge_tris().empty();}
	virtual void clear_context() {model.free_context();}
};


uobj_asteroid *uobj_asteroid::create(point const &pos, float radius, unsigned model, int tex_id, unsigned rseed_ix, unsigned lt) {

	switch (model) {
	case AS_MODEL_SPHERE: return new uobj_asteroid_sphere(pos, radius, tex_id, lt);
	case AS_MODEL_ROCK1:  return new uobj_asteroid_rock3d(pos, radius, rseed_ix, lt, 0);
	case AS_MODEL_ROCK2:  return new uobj_asteroid_rock3d(pos, radius, rseed_ix, lt, 1);
	case AS_MODEL_HMAP:   return new uobj_asteroid_hmap  (pos, radius, rseed_ix, lt);
	case AS_MODEL_VOXEL:  return new uobj_asteroid_voxel (pos, radius, rseed_ix, lt);
	case AS_MODEL_SHADER: return new uobj_asteroid_shader(pos, radius, rseed_ix, lt);
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


// *** asteroid field ***


unsigned const AST_FIELD_MODEL   = AS_MODEL_HMAP;
unsigned const NUM_AST_MODELS    = 100;
unsigned const AST_FLD_MAX_NUM   = 1000;
float    const AST_RADIUS_SCALE  = 0.04;
float    const AST_AMBIENT_SCALE = 20.0;
float    const AST_AMBIENT_VAL   = 0.15;
float    const NDIV_SCALE_AST    = 500.0;


class asteroid_model_gen_t {

	vector<uobj_asteroid *> asteroids; // FIXME: a case for boost::shared_ptr<>?

public:
	//~asteroid_model_gen_t() {clear();} // can't free after GL context has been destroyed
	bool empty() const {return asteroids.empty();}
	
	void gen(unsigned num, unsigned model) {
		RESET_TIME;
		assert(asteroids.empty());
		asteroids.resize(num);

		for (unsigned i = 0; i < num; ++i) { // create at the origin with radius=1
			asteroids[i] = uobj_asteroid::create(all_zeros, 1.0, model, ROCK_SPHERE_TEX, i);
		}
		PRINT_TIME("Asteroid Model Gen");
	}
	void draw(unsigned ix, point const &pos, vector3d const &scale, point const &camera, vector3d const &rot_axis, float rot_ang, shader_t &s) const {
		float const radius(scale.mag()); // approximate (upper bound)
		float const dist(p2p_dist(pos, camera)), dscale(NDIV_SCALE_AST*(radius/(dist + 0.1*radius)));
		if (dscale < 0.5) return; // too far/small - clip it
		glPushMatrix();
		global_translate(pos);
		scale_by(scale);
		if (rot_ang != 0.0) {rotate_about(rot_ang, rot_axis);}
		assert(ix < asteroids.size());
		assert(asteroids[ix]);
		int ndiv(max(3, min((int)ASTEROID_NDIV, int(sqrt(10.0*dscale)))));
		uobj_draw_data ddata(asteroids[ix], s, ndiv, 0, 0, 0, 0, pos, zero_vector, plus_z, plus_y, dist, radius, 1.0, 0, 1, 1, 1, 1);
		asteroids[ix]->draw_obj(ddata);
		glPopMatrix();
	}
	void clear() {
		for (vector<uobj_asteroid *>::iterator i = asteroids.begin(); i != asteroids.end(); ++i) {
			delete *i;
		}
		asteroids.clear();
	}
};

asteroid_model_gen_t asteroid_model_gen;


void uasteroid_field::gen_asteroids() {

	if (asteroid_model_gen.empty()) {asteroid_model_gen.gen(NUM_AST_MODELS, AST_FIELD_MODEL);}
	clear();
	resize((rand2() % AST_FLD_MAX_NUM) + 1);

	for (vector<uasteroid>::iterator i = begin(); i != end(); ++i) {
		i->gen(radius, AST_RADIUS_SCALE*radius);
	}
}


void setup_ship_draw_shader(shader_t &s);


void uasteroid_field::begin_render(shader_t &shader) {

	if (!shader.is_setup()) {setup_ship_draw_shader(shader);} // ambient only? not sure if this is what we want in the end
	shader.enable();
	BLACK.do_glColor();
	glEnable(GL_LIGHTING);
	colorRGBA const acolor(AST_AMBIENT_VAL, AST_AMBIENT_VAL, AST_AMBIENT_VAL, 1.0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, &acolor.R);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, &BLACK.R);
}


void uasteroid_field::end_render(shader_t &shader) {

	glDisable(GL_LIGHTING);
	shader.disable();
	glDisable(GL_TEXTURE_2D);
}


void move_in_front_of_far_clip(point_d &pos, point const &camera, float &size, float dist, float dscale);


void uasteroid_field::draw(point_d const &pos_, point const &camera, shader_t &s) {

	point_d const afpos(pos + pos_);
	if (!univ_sphere_vis(afpos, radius)) return;
	if (calc_sphere_size(afpos, camera, AST_RADIUS_SCALE*radius) < 1.0) return; // asteroids are too small/far away
	if (empty()) {gen_asteroids();}
	/*set_color(WHITE);
	WHITE.do_glColor();
	draw_sphere_at(make_pt_global(afpos), radius, N_SPHERE_DIV);*/

	for (vector<uasteroid>::const_iterator i = begin(); i != end(); ++i) {
		i->draw(afpos, camera, s);
	}
}


void uasteroid::gen(float max_dist, float max_radius) {

	assert(max_radius > 0.0 && max_radius < max_dist && max_radius);
	rgen_values();
	radius  = max_radius*rand_uniform(0.1, 1.0);
	pos     = signed_rand_vector2_spherical(max_dist - radius);
	inst_id = rand2() % NUM_AST_MODELS;
	UNROLL_3X(scale[i_] = rand_uniform2(0.5, 1.0);)
}


void uasteroid::draw(point_d const &pos_, point const &camera, shader_t &s) const {

	point_d const apos(pos_ + pos);
	if (!univ_sphere_vis(apos, radius)) return;
	if (calc_sphere_size(apos, camera, radius) < 1.0) return; // too small/far away
	asteroid_model_gen.draw(inst_id, apos, radius*scale, camera, rot_axis, rot_ang, s);
}


