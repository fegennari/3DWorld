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
#include "transform_obj.h"


bool     const ENABLE_AF_INSTS  = 1; // more efficient on large asteroid fields, but less efficient when close/sparse (due to overhead), and normals are incorrect
bool     const ENABLE_CRATERS   = 1;
bool     const ENABLE_SHADOWS   = 1;
unsigned const ASTEROID_NDIV    = 32; // for sphere model, better if a power of 2
unsigned const ASTEROID_VOX_SZ  = 64; // for voxel model
unsigned const AST_VOX_NUM_BLK  = 2; // ndiv=2x2
unsigned const NUM_VOX_AST_LODS = 3;
float    const AST_COLL_RAD     = 0.25; // limit collisions of large objects for accuracy (heightmap)
float    const AST_PROC_HEIGHT  = 0.1; // height values of procedural shader asteroids

unsigned const DEFAULT_AST_TEX  = MOON_TEX; // ROCK_TEX or MOON_TEX
unsigned const comet_tids[2]    = {ROCK_SPHERE_TEX, ICE_TEX};


extern bool no_asteroid_dust;
extern int animate2, display_mode, frame_counter, window_width, window_height;
extern float fticks;
extern colorRGBA sun_color;
extern s_object clobj0;
extern vector<us_weapon> us_weapons;
extern vector<usw_ray> t_wrays;

shader_t cached_voxel_shaders[9]; // one for each value of num_lights (0-8)
shader_t cached_proc_shaders [9];


unsigned calc_lod_pow2(unsigned max_ndiv, unsigned ndiv);


int get_spherical_texture(int tid) {
	//return ((tid == MOON_TEX) ? ROCK_SPHERE_TEX : tid);
	return tid; // ROCK_SPHERE_TEX no longer matches MOON_TEX
}


void clear_cached_shaders() {

	for (unsigned i = 0; i < 9; ++i) {
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
		ddata.set_color(ddata.color_a);
		model3d.draw_using_vbo();
		end_texture();
		enable_blend();
	}
};


class uobj_asteroid_shader : public uobj_asteroid { // unused
	int rseed_ix;

public:
	uobj_asteroid_shader(point const &pos_, float radius_, int rseed_ix_, int tid, unsigned lt)
		: uobj_asteroid(pos_, radius_, tid, lt)
	{
		c_radius = (1.0 + AST_PROC_HEIGHT)*radius;
	}
	virtual void draw_obj(uobj_draw_data &ddata) const {
		unsigned const num_lights(min(8U, num_exp_lights+2U));
		shader_t &s(cached_proc_shaders[num_lights]);
		
		if (s.is_setup()) { // already setup
			s.enable();
		}
		else {
			bind_3d_texture(get_noise_tex_3d(64, 1)); // grayscale noise
			s.set_int_prefix("num_lights", num_lights, 1); // FS
			s.set_prefix("#define NO_SPECULAR",   1); // FS (optional/optimization)
			s.set_prefix("#define NUM_OCTAVES 8", 0); // VS
			s.set_vert_shader("perlin_clouds_3d.part*+procedural_rock");
			s.set_frag_shader("ads_lighting.part*+procedural_rock");
			s.begin_shader();
			s.add_uniform_int("cloud_noise_tex", 0);
			s.add_uniform_float("time", float(rseed_ix));
			s.add_uniform_float("noise_scale",  0.1);
			s.add_uniform_float("height_scale", AST_PROC_HEIGHT);
		}
		s.set_cur_color(colorRGBA(0.5, 0.45, 0.4, 1.0)); // Note: ignores color_a
		end_texture();
		draw_sphere_vbo(all_zeros, 1.0, 3*ddata.ndiv/2, 1); // ndiv may be too large to use a vbo
		s.disable();
		if (ddata.final_pass && ddata.shader->is_setup()) {ddata.shader->enable();}
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
				int const fragment_tid(get_fragment_tid(mod_hit_pos));
				unsigned const num_fragments(min(25U, max(1U, unsigned(20*damage_val))));
				gen_moving_fragments(mod_hit_pos, num_fragments, fragment_tid, 0.25);
				if (gen_more_small_fragments()) {gen_moving_fragments(mod_hit_pos, num_fragments, fragment_tid, 0.15);}
			}
		}
		return uobj_asteroid::damage(val, type, hit_pos, source, wc);
	}
	virtual bool gen_more_small_fragments() const {return 0;}

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

	class ast_instance_render_t : public instance_render_t { // so we can add max_ndiv without needing to make it mutable as well
	public:
		unsigned max_ndiv;
		ast_instance_render_t() : max_ndiv(0) {}
	};

	float scale_val;
	vector3d xyz_scale;
	mutable upsurface surface; // FIXME: mutable so that the contained sd_sphere_vbo_d can modify its vbo indexes
	ast_instance_render_t inst_render; // could be per-LOD level (4)

public:
	uobj_asteroid_hmap(point const &pos_, float radius_, unsigned rseed_ix, int tid, unsigned lt)
		: uobj_asteroid_destroyable(pos_, radius_, tid, lt), xyz_scale(1.0, 1.0, 1.0)
	{
		surface.rgen.set_state(rseed_ix, 1);
		surface.gen(0.15, 2.0, 10, 1.0);
		surface.setup(ASTEROID_NDIV, 0.0, 0);
		surface.setup_draw_sphere(all_zeros, 1.0, 0.0, ASTEROID_NDIV, NULL);
		surface.calc_rmax();
		scale_val = 1.0/surface.rmax;
	}
	virtual void set_scale(vector3d const &scale) {xyz_scale = scale;}
	virtual bool has_custom_shadow_profile() const {return 1;}

	// Note: this class overrides draw_with_texture() because it's used instanced
	virtual void draw_with_texture(uobj_draw_data &ddata, int force_tex_id, bool no_reset_texture=0) const { // to allow overriding the texture id
		unsigned const ndiv(3*ddata.ndiv/2); // increase ndiv because we want higher resolution to capture details
		scale_by(scale_val*xyz_scale);
		ddata.set_color(ddata.color_a);
		select_texture((force_tex_id >= 0) ? force_tex_id : tex_id);
		surface.sd.draw_ndiv_pow2_vbo(ndiv); // use a vbo
		if (!no_reset_texture) {end_texture();}
	}
	virtual void draw_obj(uobj_draw_data &ddata) const {
		draw_with_texture(ddata, -1);
	}
	virtual bool draw_instanced(unsigned ndiv) { // non-const because it caches instance transforms
		ndiv = 3*ndiv/2; // increase ndiv because we want higher resolution to capture details
		scale_by(scale_val*xyz_scale);
		//unsigned const lod(calc_lod_pow2(ASTEROID_NDIV, ndiv));
		inst_render.add_cur_inst();
		inst_render.max_ndiv = max(inst_render.max_ndiv, ndiv);
		return 1;
	}
	virtual void final_draw(int xfm_shader_loc) { // non-const because it clears instance transforms
		select_texture(tex_id);
		inst_render.set_loc(xfm_shader_loc);
		surface.sd.draw_instances(inst_render.max_ndiv, inst_render);
		inst_render.max_ndiv = 0;
		// end_texture() will be called by the caller
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

	virtual void clear_context() {surface.free_context();}

	private:
	virtual float get_radius_at(point const &pt, bool exact=0) const {
		int tx, ty;
		get_tex_coords_at(pt, tx, ty);
		point **points = surface.sd.get_points();
		assert(points);
		return radius*scale_val*(xyz_scale*points[tx][ty]).mag();
	}

	void get_tex_coords_at(point const &query_pos, int &tx, int &ty) const {
		vector3d query_dir(query_pos - pos);
		rotate_point(query_dir);
		query_dir.normalize();
		get_tex_coord(query_dir, vector3d(0.0, 1.0, 0.0), ASTEROID_NDIV, ASTEROID_NDIV, tx, ty, 1);
	}
};


// FIXME: sphere collision normal / pos
class uobj_asteroid_voxel : public uobj_asteroid_destroyable {

	mutable voxel_model_space model; // FIXME: const problems with draw()
	bool have_sun_pos;
	static noise_texture_manager_t global_asteroid_ntg;

public:
	uobj_asteroid_voxel(point const &pos_, float radius_, unsigned rseed_ix, int tid, unsigned lt)
		: uobj_asteroid_destroyable(pos_, radius_, tid, lt), model(&global_asteroid_ntg, NUM_VOX_AST_LODS), have_sun_pos(0)
	{
		//RESET_TIME;
		float const gen_radius(gen_voxel_rock(model, all_zeros, 1.0, ASTEROID_VOX_SZ, AST_VOX_NUM_BLK, rseed_ix)); // will be translated to pos and scaled by radius during rendering
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
		unsigned const num_lights(min(8U, num_exp_lights+2U));
		unsigned const lod_level(min(16U/ddata.ndiv, NUM_VOX_AST_LODS-1));
		shader_t &s(cached_voxel_shaders[num_lights]);
		
		if (s.is_setup()) { // already setup
			s.enable();
		}
		else {
			s.set_int_prefix("num_lights", num_lights, 1); // FS
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
		s.add_uniform_color("color", ddata.color_a);
		glEnable(GL_CULL_FACE);
		model.core_render(s, lod_level, 0, 1); // disable view frustum culling because it's incorrect (due to transform matrices)
		glDisable(GL_CULL_FACE);
		s.disable();

		if (ddata.final_pass) {
			if (ddata.shader->is_setup()) {ddata.shader->enable();}
			end_texture();
		}
	}

	virtual bool apply_damage(float damage, point &hit_pos) {
		float const damage_radius(min(0.5, 0.1*damage));
		xform_point(hit_pos);
		point const center(hit_pos);
		bool const damaged(model.update_voxel_sphere_region(center, damage_radius, -1.0, 1, 1, &hit_pos));
		xform_point_inv(hit_pos);
		hit_pos += 0.5*(radius/ASTEROID_VOX_SZ)*(hit_pos - pos).get_norm(); // move slightly away from asteroid center
		return damaged;
	}

	virtual int get_fragment_tid(point const &hit_pos) const {
		point p(hit_pos);
		xform_point(p);
		return model.get_texture_at(p);
	}
	virtual bool gen_more_small_fragments() const {return 1;}
	virtual bool check_fragment_self_coll() const {return is_ok();}

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

noise_texture_manager_t uobj_asteroid_voxel::global_asteroid_ntg;


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
	clear_context();
}


// *** asteroid instance management ***


unsigned const AST_FIELD_MODEL   = AS_MODEL_HMAP;
unsigned const NUM_AST_MODELS    = 40;
unsigned const AST_FLD_MAX_NUM   = 1200;
unsigned const AST_BELT_MAX_NS   = 10000;
unsigned const AST_BELT_MAX_NP   = 4000;
float    const AST_RADIUS_SCALE  = 0.04;
float    const AST_AMBIENT_S     = 2.5;
float    const AST_AMBIENT_NO_S  = 10.0;
float    const AST_AMBIENT_VAL   = 0.15;
float    const AST_VEL_SCALE     = 0.0002;
float    const NDIV_SCALE_AST    = 800.0;


float get_eq_vol_scale(vector3d const &scale) {return pow(scale.x*scale.y*scale.z, 1.0f/3.0f);}


class asteroid_model_gen_t {

	vector<uobj_asteroid *> asteroids; // FIXME: a case for boost::shared_ptr<>?
	colorRGBA tex_color;

public:
	//~asteroid_model_gen_t() {clear();} // can't free after GL context has been destroyed
	bool empty() const {return asteroids.empty();}
	
	void gen(unsigned num, unsigned model) {
		RESET_TIME;
		assert(asteroids.empty());
		asteroids.resize(num);
		tex_color = texture_color(DEFAULT_AST_TEX);

		for (unsigned i = 0; i < num; ++i) { // create at the origin with radius=1
			asteroids[i] = uobj_asteroid::create(all_zeros, 1.0, model, DEFAULT_AST_TEX, i);
		}
		PRINT_TIME("Asteroid Model Gen");
	}
	uobj_asteroid const *get_asteroid(unsigned ix) const {
		assert(ix < asteroids.size());
		assert(asteroids[ix]);
		return asteroids[ix];
	}
	void draw(unsigned ix, point_d const &pos, vector3d const &scale, point const &camera, vector3d const &rot_axis, float rot_ang, shader_t &s, pt_line_drawer &pld) {
		assert(ix < asteroids.size());
		float const radius(max(scale.x, max(scale.y, scale.z)));
		float const dist(p2p_dist(pos, camera)), dscale(NDIV_SCALE_AST*(radius/(dist + 0.1*radius)));
		if (dscale < 0.5) return; // too far/small - clip it

		if (dscale < 1.0) {
			point const global_pos(make_pt_global(pos));
			pld.add_pt(global_pos, (get_player_pos() - global_pos).get_norm(), tex_color);
			return;
		}
		int ndiv(max(3, min((int)ASTEROID_NDIV, int(sqrt(4.0*dscale)))));
		fgPushMatrix();
		global_translate(pos);
		if (rot_ang != 0.0 && dscale > 2.0) {rotate_about(rot_ang, rot_axis);}
		scale_by(scale);

		// try to draw instanced if enabled, but do a normal draw without resetting the texture if that fails
		if (!asteroids[ix]->draw_instanced(ndiv)) { // maybe this case shouldn't exist and we can only create instanceable asteroids?
			uobj_draw_data ddata(asteroids[ix], &s, ndiv, 0, 0, 0, 0, pos, zero_vector, plus_z, plus_y, dist, radius, 1.0, 0, 1, 1, 1, 1);
			asteroids[ix]->draw_with_texture(ddata, -1, 1);
		}
		fgPopMatrix();
	}
	void destroy_inst(unsigned ix, point const &pos, vector3d const &scale) {
		get_asteroid(ix)->gen_fragments(pos, get_eq_vol_scale(scale));
	}
	int get_fragment_tid(unsigned ix, point const &hit_pos) const {
		return get_asteroid(ix)->get_fragment_tid(hit_pos);
	}
	void final_draw(int loc) {
		for (vector<uobj_asteroid *>::iterator i = asteroids.begin(); i != asteroids.end(); ++i) {
			(*i)->final_draw(loc); // could call only on asteroids that have been rendered at least once, but probably okay to always call
		}
	}
	void clear_contexts() {
		for (vector<uobj_asteroid *>::iterator i = asteroids.begin(); i != asteroids.end(); ++i) {
			(*i)->clear_context();
		}
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


// *** asteroid belt particles ***


float rand_gaussian_limited(rand_gen_t &rgen, float mean, float std_dev, float max_val) {

	while (1) {
		float const v(rgen.rand_gaussian(mean, std_dev));
		if (fabs(v - mean) < max_val) return v;
	}
	assert(0); // never gets here
}


class ast_belt_part_manager_t {

	vector<point> pts;
	unsigned vbo;
	sphere_t bsphere;

public:
	ast_belt_part_manager_t() : vbo(0) {}
	void clear() {pts.clear(); clear_context();}
	unsigned size () const {return pts.size();}
	bool     empty() const {return pts.empty();}

	// generates a toriodal section in the z=0 plane centered at (0,0,0)
	void gen_torus_section(unsigned npts, float ro, float ri, float max_angle) {
		pts.resize(npts);
		rand_gen_t rgen;

		if (max_angle >= PI) { // at least half the torus, use a central bounding sphere
			bsphere.pos = all_zeros;
		}
		else { // torus subsection, use a local bounding sphere
			bsphere.pos = ro*vector3d(sinf(0.5*max_angle), cosf(0.5*max_angle), 0.0);
		}
		for (unsigned i = 0; i < npts; ++i) {
			float const theta(rgen.rand_uniform(0.0, max_angle));
			float const dval(rand_gaussian_limited(rgen, 0.0, 1.0, 2.0));
			float const dplane(rand_gaussian_limited(rgen, 0.0, ri, 2.5*ri));
			vector3d const dir(sinf(theta), cosf(theta), 0.0);
			pts[i] = (ro + ri*dval)*dir;
			pts[i].z += sqrt(1.0 - 0.25*dval*dval)*dplane; // elliptical distribution
			bsphere.radius = max(bsphere.radius, p2p_dist(bsphere.pos, pts[i])); // no point radius
		}
	}
	void clear_context() {
		delete_vbo(vbo);
		vbo = 0;
	}
	void draw_vbo() const {
		assert(vbo);
		bind_vbo(vbo);
		draw_verts<vert_wrap_t>(NULL, pts.size(), GL_POINTS);
		bind_vbo(0);
	}
	bool draw(point const &center, vector3d const &rot_axis, float rot_degrees, vector3d const &size) {
		if (empty()) return 0;
		float const bradius(bsphere.radius*max(max(size.x, size.y), size.z));
		point spos(bsphere.pos*size);
		if (rot_degrees != 0.0) {rotate_vector3d(plus_z, -TO_RADIANS*rot_degrees, spos);} // rotation is backwards from GL
		rotate_vector3d_by_vr(plus_z, rot_axis, spos);
		spos += center;
		if (!univ_sphere_vis(spos, bradius)) return 0; // VFC
		if ((distance_to_camera(spos) - bradius) > 0.8) return 0; // distance/size culling
		
		create_vbo_and_upload(vbo, pts); // non-const due to this call
		fgPushMatrix();
		global_translate(center);
		rotate_into_plus_z(rot_axis);
		if (rot_degrees != 0.0) {rotate_about(rot_degrees, plus_z);}
		scale_by(size);
		draw_vbo();
		fgPopMatrix();
		return 1;
	}
};

ast_belt_part_manager_t ast_belt_part[2]; // full, partial segment


void clear_asteroid_contexts() {

	asteroid_model_gen.clear_contexts();
	for (unsigned d = 0; d < 2; ++d) {ast_belt_part[d].clear_context();}
}


unsigned const AB_NUM_PARTS_F  = 100000;
unsigned const AB_NUM_PARTS_S  = 15000;
float const AB_WIDTH_TO_RADIUS = 0.035;
float const AB_THICK_TO_WIDTH  = 0.22;
float const AST_PARTICLE_SIZE  = 0.1;
float const AST_BELT_ROT_RATE  = 0.2;
unsigned const AB_NUM_PART_SEG = 50;


// *** asteroid fields and belts ***


void set_shader_prefix_for_shadow_casters(shader_t &shader, unsigned num_shadow_casters) {

	if (num_shadow_casters > 0) {
		assert(num_shadow_casters <= 8);
		for (unsigned d = 0; d < 2; ++d) {shader.set_prefix("#define ENABLE_SHADOWS", d);} // VS/FS
	}
	shader.set_int_prefix("num_shadow_casters", num_shadow_casters, 1); // FS
}


bool set_af_color_from_system(point_d const &afpos, float radius, shader_t *shader) {

	point sun_pos; // unused
	uobject const *sobj(NULL); // unused
	return (set_uobj_color(afpos, radius, 0, 1, sun_pos, sobj, AST_AMBIENT_S, AST_AMBIENT_NO_S, shader) >= 0);
}


void uasteroid_belt_system::draw_detail(point_d const &pos_, point const &camera) const {

	point_d const afpos(pos_ + pos);
	bool const has_sun(set_af_color_from_system(afpos, radius, NULL));
	enable_blend();
	shader_t shader;

	if (AB_NUM_PARTS_F > 0 && world_mode == WMODE_UNIVERSE) { // global asteroid dust (points), only in universe mode
		if (ast_belt_part[0].empty()) {ast_belt_part[0].gen_torus_section(AB_NUM_PARTS_F, 1.0, AB_WIDTH_TO_RADIUS, TWO_PI);}
		shader.set_vert_shader("asteroid_dust");
		shader.set_frag_shader("ads_lighting.part*+asteroid_dust"); // +sphere_shadow.part*
		shader.begin_shader();
		shader.add_uniform_float("alpha_scale", 2.0);
		shader.add_uniform_color("color", texture_color(DEFAULT_AST_TEX));
		ast_belt_part[0].draw(afpos, orbital_plane_normal, 0.0, outer_radius*scale); // full/sparse
		shader.end_shader();
	}
	if (AB_NUM_PARTS_S > 0 && !no_asteroid_dust) { // local small asteroid bits (spheres)
		if (ast_belt_part[1].empty()) {ast_belt_part[1].gen_torus_section(AB_NUM_PARTS_S, 1.0, AB_WIDTH_TO_RADIUS, TWO_PI/AB_NUM_PART_SEG);}
		for (unsigned i = 0; i < 2; ++i) {shader.set_prefix("#define DRAW_AS_SPHERES", i);} // VS/FS
		if (ENABLE_SHADOWS && has_sun) {set_shader_prefix_for_shadow_casters(shader, shadow_casters.size());}
		shader.set_vert_shader("asteroid_dust");
		shader.set_frag_shader("ads_lighting.part*+sphere_shadow.part*+sphere_shadow_casters.part+asteroid_dust"); // +sphere_shadow.part*
		shader.begin_shader();
		shader.add_uniform_float("alpha_scale", 5.0);
		shader.add_uniform_float("sphere_size", AST_PARTICLE_SIZE*window_height*max_asteroid_radius);
		shader.add_uniform_int("tex0", 0);
		shader.add_uniform_color("color", WHITE);
		select_texture(DEFAULT_AST_TEX);
		if (ENABLE_SHADOWS && has_sun) {upload_shader_casters(shader);}
		set_point_sprite_mode(1);

		for (unsigned i = 0; i < AB_NUM_PART_SEG; ++i) {
			ast_belt_part[1].draw(afpos, orbital_plane_normal, (360.0*i/AB_NUM_PART_SEG), outer_radius*scale);
		}
		set_point_sprite_mode(0);
		shader.end_shader();
	}
	disable_blend();
}


void uasteroid_cont::init(point const &pos_, float radius_) {

	pos    = pos_;
	radius = radius_;
	rseed  = rand2();
}


void uasteroid_cont::gen_asteroids() {

	global_rand_gen.set_state(rseed, 123);
	ensure_asteroid_models();
	clear();
	gen_asteroid_placements();
	sort(begin(), end()); // sort by inst_id to help reduce rendering context switch time (probably irrelevant when instancing is enabled)
}


void uasteroid_field::gen_asteroid_placements() {

	resize((rand2() % AST_FLD_MAX_NUM) + 1);

	for (iterator i = begin(); i != end(); ++i) {
		i->gen_spherical(pos, radius, AST_RADIUS_SCALE*radius);
	}
}


void uasteroid_belt::gen_belt_placements(unsigned max_num, float belt_width, float belt_thickness, float max_ast_radius) {

	float rmax(0.0), plane_dmax(0.0);
	inner_radius = 0.0;
	outer_radius = radius;
	max_asteroid_radius = 0.0;
	resize((rand2() % max_num/2) + max_num/2); // 50% to 100% of max
	vector3d vxy[2] = {plus_x, plus_y};
	
	for (unsigned d = 0; d < 2; ++d) {
		rotate_vector3d_by_vr(plus_z, orbital_plane_normal, vxy[d]);
		vxy[d] *= scale[d];
	}
	for (iterator i = begin(); i != end(); ++i) {
		i->gen_belt(pos, orbital_plane_normal, vxy, outer_radius, belt_width, belt_thickness, max_ast_radius, inner_radius, plane_dmax);
		rmax = max(rmax, (p2p_dist(pos, i->pos) + i->radius));
		max_asteroid_radius = max(max_asteroid_radius, i->radius);
	}
	radius  = rmax; // update with the bounding radius of the generated asteroids
	scale.z = plane_dmax / inner_radius;
}


void uasteroid_belt_system::gen_asteroid_placements() { // radius is the asteroid belt distance from the sun

	//RESET_TIME;
	float const belt_width(AB_WIDTH_TO_RADIUS*rand_uniform2(0.9, 1.1)*radius);
	float const belt_thickness(AB_THICK_TO_WIDTH*rand_uniform2(0.9, 1.1)*belt_width);
	gen_belt_placements(AST_BELT_MAX_NS, belt_width, belt_thickness, 0.002*radius); // circular orbit, animated
	//PRINT_TIME("Asteroid Belt"); // 4ms
}


void uasteroid_belt_planet::gen_asteroid_placements() { // radius is the asteroid belt distance from the planet

	assert(planet);
	float const belt_thickness(rand_uniform2(0.08, 0.10)*bwidth);
	gen_belt_placements(AST_BELT_MAX_NP, bwidth, belt_thickness, 0.005*radius); // elliptical orbit, static
}


void uasteroid_belt_planet::init_rings(point const &pos) {
	
	assert(planet);
	//float const rscale(planet->rscale.xy_mag()/SQRT2), ri(planet->ring_ri*rscale), ro(planet->ring_ro*rscale);
	float const ri(planet->ring_ri), ro(planet->ring_ro);
	bwidth = 0.25*(ro - ri); // divide by 4 to account for the clamping of the gaussian distance function to 2*radius, and for radius vs. diameter
	init(pos, 0.5*(ro + ri)); // center of the rings
}


void uasteroid_belt::draw_bounding_torus(point const &pos_, colorRGBA const &color) const { // unused

	shader_t s;
	s.begin_color_only_shader(color);
	enable_blend();
	fgPushMatrix();
	global_translate(pos_ + pos);
	rotate_into_plus_z(orbital_plane_normal);
	scale_by(scale);
	draw_torus(inner_radius, outer_radius, 32, 32);
	fgPopMatrix();
	disable_blend();
	s.end_shader();
}


void uasteroid_belt::xform_to_local_torus_coord_space(point &pt) const {

	pt -= pos;
	rotate_vector3d_by_vr(orbital_plane_normal, plus_z, pt);
	UNROLL_3X(pt[i_] /= scale[i_];) // account for squished/elliptical torus in orbital plane
}


bool uasteroid_belt::line_might_intersect(point const &p1, point const &p2, float line_radius) const {

	if (empty()) return 0;
	if (!line_sphere_intersect(p1, p2, pos, (radius + line_radius))) return 0; // optional optimization, may not be useful
	if (sphere_might_intersect(p1, line_radius)) return 1; // first point inside torus
	point pt[2] = {p1, p2};

	for (unsigned d = 0; d < 2; ++d) {
		xform_to_local_torus_coord_space(pt[d]);
	}
	float t(0.0); // unused
	// FIXME: line_radius is incorrect in the z-dimension, should scale by <scale> to compensate for torus elliptical scaling
	float const ri(min((inner_radius + line_radius), outer_radius)), scale(1.0/outer_radius);
	// line_intersect_torus() seems to have some bug where small ro causes it to always return false, so we rescale so that ro == 1.0
	return line_torus_intersect(scale*pt[0], scale*pt[1], all_zeros, scale*ri, 1.0, t);
}


bool uasteroid_belt::sphere_might_intersect(point const &sc, float sr) const {

	if (empty()) return 0;
	if (!dist_less_than(sc, pos, (radius + sr))) return 0;
	point pt(sc);
	xform_to_local_torus_coord_space(pt);
	point p_int; // unused
	vector3d norm; // unused
	return sphere_torus_intersect(pt, sr, all_zeros, inner_radius, outer_radius, p_int, norm, 0);
}


float uasteroid_belt::get_dist_to_boundary(point const &pt) const {

	//if (sphere_might_intersect(pt, 0.0)) return 0.0;
	float const dp(dot_product((pt - pos), orbital_plane_normal));
	point const vproj(pt - dp*orbital_plane_normal);
	float const zd(max(0.0f, (fabs(dp) - inner_radius*scale.z)));
	float const xyd(max(0.0f, (fabs(outer_radius - p2p_dist(pos, vproj)) - inner_radius))); // Note: ignores scale.x and scale.y
	return sqrt(xyd*xyd + zd*zd);
}


void uasteroid_field::apply_physics(point_d const &pos_, point const &camera) { // only needs to be called when visible

	if (empty()) return;
	float const sphere_size(calc_sphere_size((pos + pos_), camera, AST_RADIUS_SCALE*radius));
	if (sphere_size < 2.0) return; // asteroids are too small/far away

	for (iterator i = begin(); i != end(); ++i) {
		i->apply_field_physics(pos, radius);
	}
	if (sphere_size < 8.0) return; // asteroids are too small/far away

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


void uasteroid_belt_system::apply_physics(upos_point_type const &pos_, point const &camera) { // only needs to be called when visible

	if (empty()) return;
	//RESET_TIME;
	calc_colliders();
	upos_point_type opn(orbital_plane_normal);

	for (iterator i = begin(); i != end(); ++i) {
		i->apply_belt_physics(pos, opn, colliders);
	}
	calc_shadowers();
	//PRINT_TIME("Physics");
	// no collision detection between asteroids as it's rare and too slow
}


void uasteroid_belt_planet::apply_physics(upos_point_type const &pos_, point const &camera) { // only needs to be called when visible

	if (empty()) return;

	if (planet) { // move all asteroids along the planet's orbit
		upos_point_type const delta_pos(planet->pos - pos);
		pos = planet->pos;

		for (iterator i = begin(); i != end(); ++i) {
			i->rot_ang += fticks*i->rot_ang0; // rotation
			i->pos     += delta_pos;
		}
	}
	calc_shadowers();
}


bool uasteroid_belt_system::is_potential_collider(uobject const &uobj) const {

	if (!dist_less_than(uobj.pos, get_camera_pos(), 250*max_asteroid_radius)) return 0; // too far away to notice
	if (!univ_sphere_vis(uobj.pos, 2.0*(uobj.radius + max_asteroid_radius)))  return 0; // expand radius somewhat
	return sphere_might_intersect(uobj.pos, uobj.radius);
}


void uasteroid_belt_system::calc_colliders() {

	colliders.clear();
	if (!animate2 || !system) return;

	for (vector<uplanet>::const_iterator p = system->planets.begin(); p != system->planets.end(); ++p) {
		if (is_potential_collider(*p)) {colliders.push_back(sphere_t(p->pos, p->radius));}

		for (vector<umoon>::const_iterator m = p->moons.begin(); m != p->moons.end(); ++m) {
			if (is_potential_collider(*m)) {colliders.push_back(sphere_t(m->pos, m->radius));}
		}
	}
}


bool uasteroid_belt_system::might_cast_shadow(uobject const &uobj) const {

	assert(system);
	if (!dist_less_than(system->sun.pos, uobj.pos, (uobj.radius + inner_radius + outer_radius))) return 0;
	if (sphere_might_intersect(uobj.pos, uobj.radius)) return 1; // uobj within asteroid belt
	float const projected_r(uobj.radius*(inner_radius + outer_radius)/p2p_dist(system->sun.pos, uobj.pos));
	return line_might_intersect(uobj.pos, system->sun.pos, projected_r);
}


void uasteroid_belt_system::calc_shadowers() {

	shadow_casters.clear();
	if (!ENABLE_SHADOWS || !system || !system->sun.is_ok()) return;

	for (vector<uplanet>::const_iterator p = system->planets.begin(); p != system->planets.end(); ++p) {
		if (might_cast_shadow(*p)) {shadow_casters.push_back(sphere_t(p->pos, p->radius));}

		for (vector<umoon>::const_iterator m = p->moons.begin(); m != p->moons.end(); ++m) {
			if (might_cast_shadow(*m)) {shadow_casters.push_back(sphere_t(m->pos, m->radius));}
		}
	}
	if (shadow_casters.size() > 8) {shadow_casters.resize(8);} // max number
	sun_pos_radius.pos    = system->sun.pos;
	sun_pos_radius.radius = system->sun.radius;
}


void uasteroid_belt_planet::calc_shadowers() {

	shadow_casters.clear();
	if (!ENABLE_SHADOWS || !planet || !planet->is_ok() || !planet->system || !planet->system->sun.is_ok()) return;
	shadow_casters.push_back(sphere_t(planet->pos, planet->radius));
	sun_pos_radius.pos    = planet->system->sun.pos;
	sun_pos_radius.radius = planet->system->sun.radius;
}


void uasteroid_cont::begin_render(shader_t &shader, unsigned num_shadow_casters, bool custom_lighting) {

	if (!shader.is_setup()) {
		set_shader_prefix_for_shadow_casters(shader, num_shadow_casters);
		if (ENABLE_CRATERS ) {shader.set_prefix("#define HAS_CRATERS",      1);} // FS
		if (ENABLE_AF_INSTS) {shader.set_prefix("#define USE_CUSTOM_XFORM", 0);} // VS
		shader.set_vert_shader("asteroid");
		string frag_shader_str("ads_lighting.part*+triplanar_texture.part+sphere_shadow.part*+sphere_shadow_casters.part");
		if (ENABLE_CRATERS) {frag_shader_str += "+rand_gen.part*+craters.part";}
		shader.set_frag_shader(frag_shader_str + "+asteroid");
		shader.begin_shader();
		shader.add_uniform_int("tex0", 0);
		shader.add_uniform_float("tex_scale", 0.5);
		shader.add_uniform_color("color", WHITE);
	}
	shader.enable();
	set_fill_mode();
	glEnable(GL_CULL_FACE);

	if (custom_lighting) {
		colorRGBA const acolor(AST_AMBIENT_VAL, AST_AMBIENT_VAL, AST_AMBIENT_VAL, 1.0);
		set_light_colors(0, acolor, BLACK, &shader);
		setup_gl_light_atten(0, 1.0, 0.0, 0.0, &shader);
	}
}


void uasteroid_cont::end_render(shader_t &shader) {

	shader.disable();
	end_texture();
	glDisable(GL_CULL_FACE);
}


void uasteroid_cont::upload_shader_casters(shader_t &s) const {

	if (!shadow_casters.empty()) {
		int const sc_loc(s.get_uniform_loc("shadow_casters"));
		assert(sc_loc >= 0); // must be available

		for (vector<sphere_t>::const_iterator sc = shadow_casters.begin(); sc != shadow_casters.end(); ++sc) {
			point const gpos(make_pt_global(sc->pos));
			float const vals[4] = {gpos.x, gpos.y, gpos.z, sc->radius};
			unsigned const ix(sc - shadow_casters.begin());
			glUniform4fv((sc_loc + ix), 1, vals);
		}
	}
	s.add_uniform_int("num_shadow_casters", shadow_casters.size());
	s.add_uniform_vector3d("sun_pos", make_pt_global(sun_pos_radius.pos));
	s.add_uniform_float("sun_radius", sun_pos_radius.radius);
	upload_mvm_to_shader(s, "fg_ViewMatrix");
}


void uasteroid_cont::draw(point_d const &pos_, point const &camera, shader_t &s, bool sun_light_already_set) {

	point_d const afpos(pos + pos_);
	if (!univ_sphere_vis(afpos, radius)) return;
	if (sphere_size_less_than(afpos, camera, AST_RADIUS_SCALE*radius, 1.0)) return; // asteroids are too small/far away
	if (empty()) {gen_asteroids();}

	// Note: can be made more efficient for asteroid_belt, since we know what the current star is, but probably not worth the complexity
	point sun_pos; // unused
	uobject const *sobj(NULL); // unused
	bool const has_sun(sun_light_already_set || set_af_color_from_system(afpos, radius, &s));

	// Note: this block and associated variables could be moved to uasteroid_belt, but we may want to use them for asteriod fields near within systems/near stars later
	if (ENABLE_SHADOWS) { // setup shadow casters
		if (!has_sun) {shadow_casters.clear();} // optional, may never get here/not make a difference
		upload_shader_casters(s);
	}
	s.add_uniform_float("crater_scale", (has_sun ? 1.0 : 0.0));
	int const loc(s.get_attrib_loc("inst_xform_matrix", 1)); // shader should include: attribute mat4 inst_xform_matrix;
	pt_line_drawer pld;

	for (const_iterator i = begin(); i != end(); ++i) {
		i->draw(pos_, camera, s, pld); // move in front of far clipping plane?
	}
	asteroid_model_gen.final_draw(loc); // flush and drawing buffers/state (will do the actual rendering here in instanced mode)

	if (!pld.empty()) {
		end_texture();
		s.set_attrib_float_array(loc, fgGetMVM().get_ptr(), 16);
		pld.draw();
	}
}


void uasteroid_cont::remove_asteroid(unsigned ix) {

	//std::swap(at(ix), back()); pop_back();
	erase(begin()+ix); // probably okay if empty after this call
}


void uasteroid_cont::detatch_asteroid(unsigned ix) {

	assert(ix < size());
	cout << "Detatch asteroid " << ix << " of " << size() << endl; // TESTING
	// create a new asteroid from the instance and copy all the parameters
	uasteroid const &inst(operator[](ix));
	uobj_asteroid *asteroid(uobj_asteroid::create(inst.pos, inst.radius, AST_FIELD_MODEL, inst.get_fragment_tid(inst.pos), inst.get_rseed(), 0)); // lt=0
	asteroid->set_vel(inst.get_velocity());
	asteroid->set_scale(inst.get_scale());
	// FIXME: rotation
	add_uobj(asteroid);
	remove_asteroid(ix);
}


void uasteroid_cont::destroy_asteroid(unsigned ix) {

	assert(ix < size());
	operator[](ix).destroy();
	remove_asteroid(ix);
}


void uasteroid::gen_base(float max_radius) {

	assert(max_radius > 0.0);
	rgen_values(); // sets rot_axis and rot_ang
	UNROLL_3X(scale[i_] = rand_uniform2(0.5, 1.0);)
	inst_id  = rand2() % NUM_AST_MODELS;
	radius   = max_radius*rand_uniform2(0.2, 1.0);
	rot_ang0 = 0.5*fabs(rand_gaussian2(0.0, 1.0)); // rotation rate
}


void uasteroid::gen_spherical(upos_point_type const &pos_offset, float max_dist, float max_radius) {

	assert(max_radius < max_dist);
	gen_base(max_radius);
	pos = pos_offset + signed_rand_vector2_spherical(max_dist - radius);
	UNROLL_3X(velocity[i_] = AST_VEL_SCALE*rand_gaussian2(0.0, 1.0);)
}


// returns the distance from the asteroid to the asteroid belt center line (torus inner radius)
void uasteroid::gen_belt(upos_point_type const &pos_offset, vector3d const &orbital_plane_normal, vector3d const vxy[2],
	float belt_radius, float belt_width, float belt_thickness, float max_radius, float &ri_max, float &plane_dmax)
{
	gen_base(max_radius);
	float const theta(TWO_PI*rand_float2());
	float const dval(rand_gaussian_limited(global_rand_gen, 0.0, 1.0, 2.0));
	float const delta_dist_to_sun(belt_width*dval);
	float const dplane(rand_gaussian_limited(global_rand_gen, 0.0, belt_thickness, 2.5*belt_thickness));
	float const dist_from_plane(sqrt(1.0 - 0.25*dval*dval)*dplane); // elliptical distribution
	vector3d const dir(vxy[0]*sinf(theta) + vxy[1]*cosf(theta)); // Note: only normalized if xscale=yscale=1
	orbital_dist = delta_dist_to_sun + belt_radius;
	pos          = pos_offset; // start at center of system/sun
	pos         += orbital_dist*dir; // move out to belt radius
	pos         += dist_from_plane*orbital_plane_normal;
	velocity     = zero_vector; // for now
	float const aoR(orbital_dist/radius), rev_rate(AST_BELT_ROT_RATE/(aoR*sqrt(aoR))); // see urev_body::gen_rotrev()
	rev_ang0     = rev_rate/(orbital_dist*orbital_dist);
	ri_max       = max(ri_max, float(radius + sqrt(delta_dist_to_sun*delta_dist_to_sun + dist_from_plane*dist_from_plane)));
	plane_dmax   = max(plane_dmax, (radius + dplane)); // approximate
}


void uasteroid::apply_field_physics(point const &af_pos, float af_radius) {

	last_coll_id = -1;
	float const vmag_sq(velocity.mag_sq()), vmax(10.0*AST_VEL_SCALE);
	if (vmag_sq > vmax*vmax) {velocity *= 0.99*vmax/sqrt(vmag_sq);} // clamp max velocity (from collisions)
	rot_ang += fticks*rot_ang0;
	pos     += velocity;
	
	if (!dist_less_than((pos - af_pos), all_zeros, (af_radius - radius))) {  // outside asteroid field bounds
		calc_reflection_angle(velocity, velocity, (af_pos - pos).get_norm()); // reflect
	}
}


void uasteroid::apply_belt_physics(upos_point_type const &af_pos, upos_point_type const &op_normal, vector<sphere_t> const &colliders) {

	upos_point_type const dir(pos - af_pos);
	rot_ang += 0.25*fticks*rot_ang0; // slow rotation
	// adjust velocity so asteroids revolve around the sun
	// Note: slightly off for asteroids not in the plane, should be cross_product(dir, op_normal).get_norm() but that's slower
	velocity = rev_ang0*cross_product(dir, op_normal);
	pos      = af_pos + orbital_dist*(pos + fticks*velocity - af_pos).get_norm(); // renormalize for constant distance

	for (vector<sphere_t>::const_iterator i = colliders.begin(); i != colliders.end(); ++i) {
		if (dist_less_than(pos, i->pos, (radius + i->radius))) {
			vector3d const coll_normal((pos - i->pos).get_norm());
			pos = i->pos + (radius + i->radius)*coll_normal;
		}
	}
}


void uasteroid::draw(point_d const &pos_, point const &camera, shader_t &s, pt_line_drawer &pld) const {

	point_d const apos(pos_ + pos);
	if (!univ_sphere_vis(apos, radius)) return;
	if (sphere_size_less_than(apos, camera, radius, 1.0)) return; // too small/far away
	asteroid_model_gen.draw(inst_id, apos, radius*scale, camera, rot_axis, rot_ang, s, pld);
}


void uasteroid::destroy() {

	int const tid(get_spherical_texture(get_fragment_tid(all_zeros)));
	float const color_scale(has_sun_lighting(pos) ? 1.0 : 2.5); // scale up the color to increase ambient
	def_explode(u_exp_size[UTYPE_ASTEROID], ETYPE_ANIM_FIRE, signed_rand_vector());
	gen_moving_fragments(pos, (40 + (rand()%20)), tid, 2.0*get_eq_vol_scale(scale), 0.5, velocity, WHITE*color_scale);
	//asteroid_model_gen.destroy_inst(inst_id, pos, radius*scale); // fragments don't have correct velocity and are very dark (no added ambient)
}


int uasteroid::get_fragment_tid(point const &hit_pos) const {
	return asteroid_model_gen.get_fragment_tid(inst_id, pos+hit_pos); // Note: asteroid_field pos is not added to hit_pos here
}


bool uasteroid::line_intersection(point const &p1, vector3d const &v12, float line_length, float line_radius, float &ldist) const {

	if (!dist_less_than(p1, pos, (radius + line_length))) return 0;
	float t, rdist;
	if (!line_intersect_sphere(p1, v12, pos, (radius+line_radius), rdist, ldist, t)) return 0; // line intersects asteroid bounding sphere
	// transform line into asteroid's translated and scaled coord space
	point p1b;
	vector3d v12b;
	UNROLL_3X(p1b[i_] = (p1[i_] - pos[i_])/scale[i_]; v12b[i_] = v12[i_]/scale[i_];)
	v12b.normalize();
	float rdist2, ldist2, t2; // unused
	if (!line_intersect_sphere(p1b, v12b, all_zeros, (radius+line_radius), rdist2, ldist2, t2)) return 0; // somewhat more accurate test
	return (t > 0.0 && ldist <= line_length);
}


// *** rand_spawn_mixin / uobject_rand_spawn_t / ucomet ***


void rand_spawn_mixin::gen_rand_pos() {

	for (unsigned iter_count = 0; iter_count < 10; ++iter_count) { // limit the number of iterations
		vector3d const dir(first_pos ? signed_rand_vector(max_cdist) : signed_rand_vector_norm(max_cdist));
		obj_pos = get_player_pos() + dir;
		if (!player_pdu.valid || !univ_sphere_vis(obj_pos, obj_radius)) break; // don't spawn in player's view
	}
	first_pos = 0;
	pos_valid = 1;
}


bool rand_spawn_mixin::okay_to_respawn() const {

	if (clobj0.system < 0) return 0; // player not near a system, so don't respawn
	if (player_ship().get_velocity().mag() > 0.2) return 0; // player ship moving too quickly, don't respawn
	return 1;
}


bool rand_spawn_mixin::needs_respawned() const {

	return !dist_less_than(obj_pos, get_player_pos(), (univ_sphere_vis(obj_pos, obj_radius) ? 2.0 : 1.1)*max_cdist); // further if visible
}


uobject_rand_spawn_t *uobject_rand_spawn_t::create(unsigned type, float radius_, float dmax, float vmag) {
	
	switch (type) {
	case SPO_COMET: {return new ucomet(radius_, dmax, vmag);}
	default: assert(0);
	}
	return NULL;
}


uobject_rand_spawn_t::uobject_rand_spawn_t(float radius_, float dmax, float vmag) : rand_spawn_mixin(pos, radius_, dmax) {

	assert(radius_ > 0.0 && dmax > 0.0 && vmag >= 0.0); // sanity checks
	radius   = c_radius = radius_;
	velocity = ((vmag == 0.0) ? zero_vector : signed_rand_vector_norm(vmag));
	pos      = all_zeros; // temporary
}


void uobject_rand_spawn_t::mark_pos_invalid() {

	pos_valid = 0; // respawn
	flags    |= OBJ_FLAGS_NCOL; // disable collisions until respawned
}


void uobject_rand_spawn_t::gen_pos() {

	if (!okay_to_respawn()) return;
	gen_rand_pos();
	if (dot_product(velocity, dir) > 0.0) {velocity *= -1.0;} // make it approach the camera
	flags &= ~OBJ_FLAGS_NCOL; // clear the flag
}


void uobject_rand_spawn_t::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass, int align, unsigned eflags, free_obj const *parent_) {

	free_obj::explode(0.2*damage, 0.2*bradius, etype, edir, exp_time, wclass, align, eflags, parent_);

	if (status == 1) { // actually exploded and was destroyed
		status = 0;
		mark_pos_invalid(); // respawn
	}
}


void uobject_rand_spawn_t::advance_time(float timestep) {

	if (!pos_valid) {gen_pos();}
	free_obj::advance_time(timestep);
	if (needs_respawned()) {mark_pos_invalid();} // if too far from the player, respawn at a different location
}


ucomet::ucomet(float radius_, float dmax, float vmag) : uobject_rand_spawn_t(radius_, dmax, vmag), sun_pos(all_zeros) {

	draw_rscale = 10.0; // increase draw radius so that tails are visible even when the comet body is not
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
		fgPushMatrix();
		ddata.color_a = (i ? colorRGBA(2.0, 1.2, 1.0, 1.0) : WHITE); // less blue for ice
		if (i == 1) {ddata.set_uobj_specular(0.8, 50.0);} // not sure if this actually works
		asteroid_model_gen.get_asteroid((inst_ids[i] + i) % NUM_AST_MODELS)->draw_with_texture(ddata, comet_tids[i], 1);
		if (i == 1) {ddata.end_specular();}
		fgPopMatrix();
	}
	end_texture();

	if (temperature > 1.0) {
		float const glow_weight(CLIP_TO_01(get_true_temp()/40.0f)); // 1.0 if camera is facing the lit side?
		colorRGBA color(sun_color, glow_weight), color2(color, 0.0);
		ddata.draw_engine(color, all_zeros, 4.0, 1.0, all_zeros); // coma
		ddata.draw_ship_flares(color);

		if (animate2) { // create tails
			color.alpha *= 0.5;

			if (temperature > 4.0 && sun_pos != all_zeros) { // ion tail points away from the sun
				rand_gen_t rgen;
				rgen.set_state(inst_ids[0], inst_ids[1]);

				for (unsigned i = 0; i < 10; ++i) { // Note: could use a procedural 1D texture similar to planet rings instead of creating multiple rays
					vector3d const dir(radius*rgen.signed_rand_vector());
					point const pos2(pos + 30.0*radius*rgen.rand_uniform(0.75, 1.0)*(pos - sun_pos).get_norm() + 2.0*dir);
					float const width(rgen.rand_uniform(0.5, 1.0));
					t_wrays.push_back(usw_ray(1.0*width*radius, 3.0*width*radius, (pos + 0.3*dir), pos2, color, color2));
				}
			}
			if (temperature > 2.0 && ddata.ndiv > 6) { // dust tail follows velocity/path
				float const vmag(velocity.mag()), fnum(8000.0*fticks*vmag);
				unsigned num(fnum);
				if (rand_float() < (fnum - num)) {++num;}

				for (unsigned i = 0; i < num; ++i) {
					vector3d const delta(signed_rand_vector()), pvel(0.2*vmag*delta);
					gen_particle(PTYPE_GLOW, color, color2, unsigned(2.0*(2.5 - delta.mag())*TICKS_PER_SECOND),
						(pos + 0.75*delta*radius), pvel, 0.25*radius, 0.0, ALIGN_NEUTRAL, 0);
				}
			}
		}
	}
}


