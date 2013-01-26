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
#include "ship_util.h" // for gen_particle


unsigned const ASTEROID_NDIV   = 32; // for sphere model, better if a power of 2
unsigned const ASTEROID_VOX_SZ = 64; // for voxel model
float    const AST_COLL_RAD    = 0.25; // limit collisions of large objects for accuracy (heightmap)
float    const AST_PROC_HEIGHT = 0.1; // height values of procedural shader asteroids


extern int animate2, iticks;
extern float fticks;
extern colorRGBA sun_color;
extern s_object clobj0;
extern vector<us_weapon> us_weapons;
extern vector<usw_ray> t_wrays;

shader_t cached_voxel_shaders[8]; // one for each value of num_lights
shader_t cached_proc_shaders [8];


void set_uniform_atten_lighting(int light);


void clear_cached_shaders() {

	for (unsigned i = 0; i < 8; ++i) {
		cached_voxel_shaders[i].end_shader();
		cached_proc_shaders [i].end_shader();
	}
}


class uobj_asteroid_sphere : public uobj_asteroid {

public:
	uobj_asteroid_sphere(point const &pos_, float radius_, int tid, unsigned lt) : uobj_asteroid(pos_, radius_, tid, lt) {}
	virtual void draw_obj(uobj_draw_data &ddata) const {ddata.draw_asteroid(tex_id);}
};


class uobj_asteroid_rock3d : public uobj_asteroid {

	rock_shape3d model3d;

public:
	uobj_asteroid_rock3d(point const &pos_, float radius_, unsigned rseed_ix, int tid, unsigned lt, int type)
		: uobj_asteroid(pos_, radius_, tid, lt)
	{
		model3d.gen_rock(24, 1.0, rseed_ix, type); // pos starts at and stays at all_zeros
		model3d.set_texture(tex_id, 0.2);
	}
	~uobj_asteroid_rock3d() {model3d.destroy();}

	virtual void draw_obj(uobj_draw_data &ddata) const {
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(tex_id); return;}
		ddata.color_a.do_glColor();
		model3d.draw();
		end_texture();
		enable_blend();
	}
};


class uobj_asteroid_shader : public uobj_asteroid {
	int rseed_ix;

public:
	uobj_asteroid_shader(point const &pos_, float radius_, int rseed_ix_, int tid, unsigned lt)
		: uobj_asteroid(pos_, radius_, tid, lt)
	{
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
			s.set_frag_shader("ads_lighting.part*+procedural_rock");
			s.begin_shader();
			s.add_uniform_int("cloud_noise_tex", 0);
			s.add_uniform_float("time", float(rseed_ix));
			s.add_uniform_float("noise_scale",  0.1);
			s.add_uniform_float("height_scale", AST_PROC_HEIGHT);
		}
		colorRGBA(0.5, 0.45, 0.4, 1.0).do_glColor(); // Note: ignores color_a
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
	uobj_asteroid_destroyable(point const &pos_, float radius_, int tid, unsigned lt) : uobj_asteroid(pos_, radius_, tid, lt) {}
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
	uobj_asteroid_hmap(point const &pos_, float radius_, unsigned rseed_ix, int tid, unsigned lt)
		: uobj_asteroid_destroyable(pos_, radius_, tid, lt)
	{
		surface.rgen.set_state(rseed_ix, 1);
		surface.gen(0.15, 2.0, 10, 1.0);
		surface.setup(ASTEROID_NDIV, 0.0, 0);
		surface.setup_draw_sphere(all_zeros, 1.0, 0.0, ASTEROID_NDIV, NULL);
		surface.calc_rmax();
		scale_val = 1.0/surface.rmax;
	}

	// Note: this class overrides draw_with_texture() because it's used instanced
	virtual void draw_with_texture(uobj_draw_data &ddata, int force_tex_id) const { // to allow overriding the texture id
		if (ddata.ndiv <= 4) {ddata.draw_asteroid(force_tex_id); return;}
		if (scale_val != 1.0) {uniform_scale(scale_val);}
		ddata.color_a.do_glColor();
		select_texture(force_tex_id);
		surface.sd.draw_ndiv_pow2(ddata.ndiv); // use dlist?
		end_texture();
	}
	virtual void draw_obj(uobj_draw_data &ddata) const {
		draw_with_texture(ddata, tex_id);
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
	uobj_asteroid_voxel(point const &pos_, float radius_, unsigned rseed_ix, int tid, unsigned lt)
		: uobj_asteroid_destroyable(pos_, radius_, tid, lt), have_sun_pos(0)
	{
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
			s.set_frag_shader("ads_lighting.part*+triplanar_texture.part+procedural_texture.part+voxel_texture.part+voxel_asteroid");
			s.begin_shader();
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
		ddata.color_a.do_glColor();
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
		bool const contains(r > radius && dist_less_than(c, pos, r-radius)); // sphere contains asteroid
		if (contains && !ip.calc_int) return 1;
		r /= radius; // scale to 1.0
		point p(c);
		xform_point(p);

		if (contains) {
			ip.p_int = all_zeros; // report asteroid center as collision point
		}
		else if (!model.sphere_intersect(p, r, (ip.calc_int ? &ip.p_int : NULL))) {
			return 0;
		}
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
	case AS_MODEL_ROCK1:  return new uobj_asteroid_rock3d(pos, radius, rseed_ix, tex_id, lt, 0);
	case AS_MODEL_ROCK2:  return new uobj_asteroid_rock3d(pos, radius, rseed_ix, tex_id, lt, 1);
	case AS_MODEL_HMAP:   return new uobj_asteroid_hmap  (pos, radius, rseed_ix, tex_id, lt);
	case AS_MODEL_VOXEL:  return new uobj_asteroid_voxel (pos, radius, rseed_ix, tex_id, lt);
	case AS_MODEL_SHADER: return new uobj_asteroid_shader(pos, radius, rseed_ix, tex_id, lt);
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
float    const AST_VEL_SCALE     = 0.0002;
float    const NDIV_SCALE_AST    = 800.0;


float get_eq_vol_scale(vector3d const &scale) {return pow(scale.x*scale.y*scale.z, 1.0f/3.0f);}


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
			asteroids[i] = uobj_asteroid::create(all_zeros, 1.0, model, MOON_TEX/*ROCK_TEX*/, i);
		}
		PRINT_TIME("Asteroid Model Gen");
	}
	uobj_asteroid *get_asteroid(unsigned ix) const {
		assert(ix < asteroids.size());
		assert(asteroids[ix]);
		return asteroids[ix];
	}
	void draw(unsigned ix, point const &pos, vector3d const &scale, point const &camera, vector3d const &rot_axis, float rot_ang, shader_t &s) const {
		float const radius(max(scale.x, max(scale.y, scale.z)));
		float const dist(p2p_dist(pos, camera)), dscale(NDIV_SCALE_AST*(radius/(dist + 0.1*radius)));
		if (dscale < 0.5) return; // too far/small - clip it
		glPushMatrix();
		global_translate(pos);
		if (rot_ang != 0.0) {rotate_about(rot_ang, rot_axis);}
		scale_by(scale);
		int ndiv(max(3, min((int)ASTEROID_NDIV, int(sqrt(4.0*dscale)))));
		uobj_asteroid const *const asteroid(get_asteroid(ix));
		uobj_draw_data ddata(asteroid, s, ndiv, 0, 0, 0, 0, pos, zero_vector, plus_z, plus_y, dist, radius, 1.0, 0, 1, 1, 1, 1);
		asteroid->draw_obj(ddata);
		glPopMatrix();
	}
	void destroy_inst(unsigned ix, point const &pos, vector3d const &scale) {
		get_asteroid(ix)->gen_fragments(pos, get_eq_vol_scale(scale));
	}
	int get_fragment_tid(unsigned ix, point const &hit_pos) const {
		return get_asteroid(ix)->get_fragment_tid(hit_pos);
	}
	void clear() {
		for (vector<uobj_asteroid *>::iterator i = asteroids.begin(); i != asteroids.end(); ++i) {
			delete *i;
		}
		asteroids.clear();
	}
};

asteroid_model_gen_t asteroid_model_gen;


void ensure_asteroid_models() {

	if (asteroid_model_gen.empty()) {asteroid_model_gen.gen(NUM_AST_MODELS, AST_FIELD_MODEL);}
}


void uasteroid_field::init(point const &pos_, float radius_) {

	pos    = pos_;
	radius = radius_;
	rseed  = rand2();
}


void uasteroid_field::gen_asteroids() {

	global_rand_gen.set_state(rseed, 123);
	ensure_asteroid_models();
	clear();
	resize((rand2() % AST_FLD_MAX_NUM) + 1);

	for (iterator i = begin(); i != end(); ++i) {
		i->gen(pos, radius, AST_RADIUS_SCALE*radius);
	}
}


void uasteroid_field::apply_physics(point_d const &pos_, point const &camera) { // only needs to be called when visible

	if (empty() || calc_sphere_size((pos + pos_), camera, AST_RADIUS_SCALE*radius) < 1.0) return; // asteroids are too small/far away

	for (iterator i = begin(); i != end(); ++i) {
		i->apply_physics(pos, radius);
	}

	// check for collisions between asteroids
	float const mult(0.5*AF_GRID_SZ/radius);

	for (unsigned z = 0; z < AF_GRID_SZ; ++z) {
		for (unsigned y = 0; y < AF_GRID_SZ; ++y) {
			for (unsigned x = 0; x < AF_GRID_SZ; ++x) {
				grid[z][y][x].resize(0);
			}
		}
	}
	for (iterator i = begin(); i != end(); ++i) {
		unsigned const ix(i - begin());
		unsigned bnds[3][2];

		for (unsigned d = 0; d < 3; ++d) {
			bnds[d][0] = max(0, min((int)AF_GRID_SZ-1, int((i->pos[d] - i->radius - (pos[d] - radius))*mult)));
			bnds[d][1] = max(0, min((int)AF_GRID_SZ-1, int((i->pos[d] + i->radius - (pos[d] - radius))*mult)));
		}
		for (unsigned z = bnds[2][0]; z <= bnds[2][1]; ++z) {
			for (unsigned y = bnds[1][0]; y <= bnds[1][1]; ++y) {
				for (unsigned x = bnds[0][0]; x <= bnds[0][1]; ++x) {
					vector<unsigned short> &gv(grid[z][y][x]);

					for (vector<unsigned short>::const_iterator g = gv.begin(); g != gv.end(); ++g) {
						uasteroid &j(at(*g));
						if (j.last_coll_id == ix) continue; // already collided with this object
						float const dmin(i->radius + j.radius);
						if (!dist_less_than(i->pos, j.pos, dmin)) continue;
						vector3d norm_dir(i->pos - j.pos);
						UNROLL_3X(norm_dir[i_] /= (i->get_scale()[i_]*j.get_scale()[i_]);)
						if (norm_dir.mag_sq() < dmin*dmin) continue;
						// see free_obj::coll_physics(): v1' = v1*(m1 - m2)/(m1 + m2) + v2*2*m2/(m1 + m2)
						float const mi(i->get_rel_mass()), mj(j.get_rel_mass());
						vector3d const &vi(i->get_velocity()), &vj(j.get_velocity());
						vector3d const vin(vi*(mi - mj)/(mi + mj) + vj*2*mj/(mi + mj));
						vector3d const vjn(vj*(mj - mi)/(mj + mi) + vi*2*mi/(mj + mi));
						i->set_velocity(vin); i->last_coll_id = *g;
						j.set_velocity (vjn); j.last_coll_id  = ix; // FIXME: move so they don't collide?
					}
					gv.push_back(ix);
				}
			}
		}
	}
}


void uasteroid_field::begin_render(shader_t &shader) {

	if (!shader.is_setup()) {
		shader.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
		shader.set_vert_shader("asteroid");
		shader.set_frag_shader("ads_lighting.part*+triplanar_texture.part+asteroid");
		shader.begin_shader();
		shader.add_uniform_int("tex0", 0);
		shader.add_uniform_float("tex_scale", 0.5);
	}
	shader.enable();
	BLACK.do_glColor();
	glEnable(GL_LIGHTING);
	colorRGBA const acolor(AST_AMBIENT_VAL, AST_AMBIENT_VAL, AST_AMBIENT_VAL, 1.0);
	int const light(GL_LIGHT0);
	glLightfv(light, GL_AMBIENT, &acolor.R);
	glLightfv(light, GL_DIFFUSE, &BLACK.R);
	set_uniform_atten_lighting(light);
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

	for (const_iterator i = begin(); i != end(); ++i) {
		i->draw(pos_, camera, s);
	}
}


void uasteroid_field::destroy_asteroid(unsigned ix) {

	assert(ix < size());
	operator[](ix).destroy();
	erase(begin()+ix); // probably okay if empty after this call
}


void uasteroid::gen(upos_point_type const &pos_offset, float max_dist, float max_radius) {

	assert(max_radius > 0.0 && max_radius < max_dist && max_radius);
	rgen_values();
	radius   = max_radius*rand_uniform(0.1, 1.0);
	pos      = pos_offset + signed_rand_vector2_spherical(max_dist - radius);
	rot_ang0 = 0.5*fabs(rand_gaussian2(0.0, 1.0)); // rotation rate
	UNROLL_3X(velocity[i_] = AST_VEL_SCALE*rand_gaussian2(0.0, 1.0);)
	inst_id  = rand2() % NUM_AST_MODELS;
	UNROLL_3X(scale[i_] = rand_uniform2(0.5, 1.0);)
}


void uasteroid::apply_physics(point const &af_pos, float af_radius) {

	last_coll_id = -1;
	float const vmag(velocity.mag()), vmax(10.0*AST_VEL_SCALE);
	if (vmag > vmax) {velocity *= vmax/vmag;} // clamp max velocity (from collisions)
	rot_ang += fticks*rot_ang0;
	pos     += velocity;
	if (!dist_less_than((pos - af_pos), all_zeros, (af_radius - radius))) {velocity *= -1.0;} // outside asteroid field bounds, reverse
}


void uasteroid::draw(point_d const &pos_, point const &camera, shader_t &s) const {

	point_d const apos(pos_ + pos);
	if (!univ_sphere_vis(apos, radius)) return;
	if (calc_sphere_size(apos, camera, radius) < 1.0) return; // too small/far away
	//((last_coll_id >= 0) ? RED : WHITE).do_glColor(); // testing
	asteroid_model_gen.draw(inst_id, apos, radius*scale, camera, rot_axis, rot_ang, s);
}


void uasteroid::destroy() {

	//int const tid(get_fragment_tid(all_zeros));
	int const tid(ROCK_SPHERE_TEX);
	def_explode(u_exp_size[UTYPE_ASTEROID], ETYPE_ANIM_FIRE, signed_rand_vector());
	gen_moving_fragments(pos, (40 + (rand()%20)), tid, 2.0*get_eq_vol_scale(scale), 0.5, velocity, WHITE*2.5); // scale up the color to increase ambient
	//asteroid_model_gen.destroy_inst(inst_id, pos, radius*scale); // fragments don't have correct velocity and are very dark (no added ambient)
}


int uasteroid::get_fragment_tid(point const &hit_pos) const {
	return asteroid_model_gen.get_fragment_tid(inst_id, pos+hit_pos); // Note: asteroid_field pos is not added to hit_pos here
}


// *** uobject_rand_spawn_t / ucomet ***


unsigned const comet_tids[2] = {ROCK_SPHERE_TEX, ICE_TEX};


uobject_rand_spawn_t *uobject_rand_spawn_t::create(unsigned type, float radius_, float dmax, float vmag) {
	
	switch (type) {
	case SPO_COMET: {return new ucomet(radius_, dmax, vmag);}
	default: assert(0);
	}
	return NULL;
}


uobject_rand_spawn_t::uobject_rand_spawn_t(float radius_, float dmax, float vmag) : first_pos(1), pos_valid(0), max_cdist(dmax) {

	assert(radius_ > 0.0 && dmax > 0.0 && vmag >= 0.0); // sanity checks
	radius   = c_radius = radius_;
	velocity = ((vmag == 0.0) ? zero_vector : signed_rand_vector(vmag));
	pos      = all_zeros; // temporary
}


void uobject_rand_spawn_t::mark_pos_invalid() {

	pos_valid = 0; // respawn
	flags    |= OBJ_FLAGS_NCOL; // disable collisions until respawned
}


void uobject_rand_spawn_t::gen_pos() {

	if (clobj0.galaxy < 0) return; // player not within a galaxy, so don't respawn
	if (player_ship().get_velocity().mag() > 0.2) return; // player ship moving too quickly, don't respawn

	for (unsigned iter_count = 0; iter_count < 10; ++iter_count) { // limit the number of iterations
		vector3d const dir(first_pos ? signed_rand_vector(max_cdist) : signed_rand_vector_norm(max_cdist));
		pos = get_player_pos() + dir;
		if (!player_pdu.valid || !univ_sphere_vis(pos, radius)) break; // don't spawn in player's view
	}
	if (dot_product(velocity, dir) > 0.0) {velocity *= -1.0;} // make it approach the camera
	first_pos = 0;
	pos_valid = 1;
	flags    &= ~OBJ_FLAGS_NCOL; // clear the flag
}


void uobject_rand_spawn_t::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass, int align, unsigned eflags, free_obj const *parent_) {

	free_obj::explode(0.2*damage, 0.2*bradius, etype, edir, exp_time, wclass, align, eflags, parent_);

	if (status == 1) { // actually exploded and was destroyed
		status    = 0;
		mark_pos_invalid(); // respawn
	}
}


void uobject_rand_spawn_t::advance_time(float timestep) {

	if (!pos_valid) {gen_pos();}
	free_obj::advance_time(timestep);

	if (!dist_less_than(pos, get_player_pos(), (univ_sphere_vis(pos, radius) ? 2.0 : 1.01)*max_cdist)) { // further if visible
		mark_pos_invalid(); // if too far from the player, respawn at a different location
	}
}


ucomet::ucomet(float radius_, float dmax, float vmag) : uobject_rand_spawn_t(radius_, dmax, vmag), sun_pos(all_zeros) {

	dir = signed_rand_vector_norm();
	gen_inst_ids();
}


void ucomet::gen_pos() {

	gen_inst_ids();
	uobject_rand_spawn_t::gen_pos();
}


void ucomet::gen_inst_ids() {

	for (unsigned d = 0; d < 2; ++d) { // need two instances, ice and rock
		inst_ids[d] = rand2() % NUM_AST_MODELS;
	}
	if (inst_ids[0] == inst_ids[1] && NUM_AST_MODELS > 2) { // same inst_ids
		inst_ids[1] = (inst_ids[0] + 1) % NUM_AST_MODELS; // make them differ by 1
	}
}


void ucomet::set_temp(float temp, point const &tcenter, free_obj const *source) {

	sun_pos = tcenter;
	free_obj::set_temp(temp, tcenter, source);
}


float ucomet::damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {

	if (rand_float() > 1.0E-5*val/radius) {return 0.0;} // higher damage = higher chance of destroying the comet

	for (unsigned i = 0; i < 2; ++i) { // mixed rock and ice
		gen_moving_fragments(pos, (20 + (rand()%12)), comet_tids[i], 1.0, 1.0, velocity, WHITE);
	}
	return free_obj::damage(val, type, hit_pos, source, wc); // will be destroyed
}


void ucomet::draw_obj(uobj_draw_data &ddata) const {

	if (!pos_valid) return;
	ensure_asteroid_models();

	for (unsigned i = 0; i < 2; ++i) { // mixed rock and ice
		glPushMatrix();
		ddata.color_a = (i ? colorRGBA(2.0, 1.2, 1.0, 1.0) : WHITE); // less blue for ice
		if (i == 1) {set_specular(0.8, 50.0);} // not sure if this actually works
		asteroid_model_gen.get_asteroid((inst_ids[i] + i) % NUM_AST_MODELS)->draw_with_texture(ddata, comet_tids[i]);
		if (i == 1) {set_specular(0.0, 1.0);}
		glPopMatrix();
	}
	if (temperature > 1.0) {
		float const glow_weight(CLIP_TO_01(get_true_temp()/40.0f)), z_offset(0.0); // 1.0 if camera is facing the lit side?
		colorRGBA color(sun_color), color2(color);
		color.alpha  = glow_weight;
		color2.alpha = 0.0;
		ddata.enable_ship_flares(color);
		ddata.draw_engine(color, all_zeros, 4.0, 1.0, all_zeros, z_offset);
		ddata.disable_ship_flares();

		if (animate2) { // create tails
			color.alpha *= 0.5;

			if (temperature > 4.0 && sun_pos != all_zeros) { // ion tail points away from the sun
				rand_gen_t rgen;
				rgen.set_state(inst_ids[0], inst_ids[1]);

				for (unsigned i = 0; i < 10; ++i) {
					vector3d const dir(radius*rgen.signed_rand_vector());
					point const pos2(pos + 30.0*radius*rgen.rand_uniform(0.75, 1.0)*(pos - sun_pos).get_norm() + 2.0*dir);
					float const width(rgen.rand_uniform(0.5, 1.0));
					t_wrays.push_back(usw_ray(1.0*width*radius, 2.5*width*radius, (pos + 0.3*dir), pos2, color, color2));
				}
			}
			if (temperature > 2.0 && ddata.ndiv > 6) { // dust tail follows velocity/path
				// FIXME: iterate and use iticks?
				vector3d const delta(signed_rand_vector()), pvel(0.2*velocity.mag()*delta);
				gen_particle(PTYPE_GLOW, color, color2, unsigned(2.0*(3.0 - delta.mag())*TICKS_PER_SECOND),
					(pos + 0.75*delta*radius), pvel, 0.3*radius, 0.0, ALIGN_NEUTRAL, 0);
			}
		}
	}
}


