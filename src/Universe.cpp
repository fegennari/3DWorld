// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 9/4/02
#include "function_registry.h"
#include "universe.h"
#include "explosion.h"
#include "textures_3dw.h"
#include "draw_utils.h"
#include "shaders.h"
#include "gl_ext_arb.h"


// temperatures
float const CGAS_TEMP        = 5.00;
float const MIN_LAND_TEMP    = 5.50;
float const MIN_COLONY_TEMP  = 6.00;
float const MIN_PLANT_TEMP   = 7.00;
float const MIN_LIVE_TEMP    = 9.00;
float const FREEZE_TEMP      = 12.0;
float const MAX_LIVE_TEMP    = 20.0;
float const MAX_PLANT_TEMP   = 25.0;
float const MAX_COLONY_TEMP  = 28.0;
float const MAX_LAND_TEMP    = 29.0;
float const BOIL_TEMP        = 30.0;
float const NO_AIR_TEMP      = 32.0;
float const NDIV_SIZE_SCALE  = 12.0;

bool const USE_HEIGHTMAP     = 1; // for planets/moons/(stars)
bool const CACHE_SPHERE_DATA = 1; // more memory but much faster rendering
bool const SHOW_SPHERE_TIME  = 0; // debugging
bool const USE_SPHERE_DLISTS = 1; // sometimes faster on newer video cards, slower on older ones

unsigned const MAX_TRIES     = 100;
unsigned const SPHERE_MAX_ND = 256;
unsigned const RING_TEX_SZ   = 256;
int   const RAND_CONST       = 1;
float const ROTREV_TIMESCALE = 1.0;
float const ROT_RATE_CONST   = 0.5*ROTREV_TIMESCALE;
float const REV_RATE_CONST   = 1.0*ROTREV_TIMESCALE;
float const STAR_BRIGHTNESS  = 1.4;
float const MIN_TEX_OBJ_SZ   = 3.0;
float const MAX_WATER        = 0.75;
float const GLOBAL_AMBIENT   = 0.25;


bool have_sun(1);
int uxyz[3] = {0, 0, 0};
unsigned char water_c[3] = {0}, ice_c[3] = {0};
unsigned char const *const wic[2] = {water_c, ice_c};
float univ_sun_rad(AVG_STAR_SIZE), univ_temp(0.0);
point last_camera(all_zeros), univ_sun_pos(all_zeros);
colorRGBA sun_color(SUN_LT_C);
s_object current;
universe_t universe; // the top level universe
pt_line_drawer universe_pld;


extern bool univ_planet_lod, disable_shaders;
extern int window_width, window_height, animate2, display_mode, onscreen_display, iticks, frame_counter;
extern float tan_term, sin_term, fticks, tfticks;
extern colorRGBA bkg_color;
extern exp_type_params et_params[];
extern GLUquadricObj* quadric;
extern rand_gen_t global_rand_gen;


void draw_cell(int const cxyz[3], camera_mv_speed const &cmvs, s_object const &clobj, unsigned pass, bool no_move);
void draw_cell_contents(ucell &cell, camera_mv_speed const &cmvs, s_object const &clobj, unsigned pass, bool no_move);

int gen_rand_seed1(point const &center);
int gen_rand_seed2(point const &center);
bool is_shadowed(point const &pos, float radius, bool expand, ussystem const &sol, uobject const *&sobj);
unsigned get_texture_size(float psize);
int moon_shadowed_by_planet(umoon const &moon);
bool get_gravity(s_object &result, point pos, vector3d &gravity, int offset);
void set_sun_loc_color(point const &pos, colorRGBA const &color, float radius, bool shadowed, bool no_ambient, float a_scale, float d_scale);
void set_light_galaxy_ambient_only();
void set_ambient_color(colorRGBA const &color);
void set_lighting_params();
void get_point_of_collision(s_object const &result, point const &pos, point &cpos);



point get_scaled_upt() {
	return point(CELL_SIZE*uxyz[0], CELL_SIZE*uxyz[1], CELL_SIZE*uxyz[2]);
}

inline void offset_pos(point &pos) { // unrolled
	pos += get_scaled_upt();
}

inline float temp_to_deg_c(float temp) {
	return (100.0*(temp - FREEZE_TEMP)/(BOIL_TEMP - FREEZE_TEMP)); // 10 => 0, 20 => 100
}

inline float temp_to_deg_f(float temp) {
	return (1.8*temp_to_deg_c(temp) + 32.0);
}

inline float calc_sphere_size(point const &pos, point const &camera, float radius, float d_adj=0.0) {
	return ((float)window_width)*radius/max(TOLERANCE, (p2p_dist(pos, camera) + d_adj)); // approx. in pixels
}


s_object get_shifted_sobj(s_object const &sobj) {

	s_object sobj2(sobj);
	UNROLL_3X(sobj2.cellxyz[i_] += uxyz[i_];)
	return sobj2;
}


void setup_universe_fog(s_object const &closest, bool damaged) {

	if (damaged) {
		glEnable(GL_FOG);
		glFogfv(GL_FOG_COLOR, (float *)&RED);
		glFogf(GL_FOG_END, 0.02); // interesting effect - shows red color towards star
	}
	if (closest.type == UTYPE_PLANET) {
		assert(closest.object != NULL);

		if (p2p_dist(get_player_pos(), closest.object->pos) < PLANET_ATM_RSCALE*closest.object->radius) {
			float const atmos(closest.get_planet().atmos);
			
			if (atmos > 0.0) {
				colorRGBA color(0.8, 0.8, 0.8, 0.33*atmos);
				blend_color(color, color, closest.get_star().color, 0.7, 0);
				add_camera_filter(color, 4, -1, CAM_FILT_FOG); // ***texture?, GL_FOG? ***
			}
		}
	}
}


uobject *line_intersect_universe(point const &start, vector3d const &dir, float length, float line_radius, float &dist) {

	point coll;
	s_object target;

	if (universe.get_trajectory_collisions(target, coll, dir, start, length, line_radius)) { // destroy and query
		if (target.is_solid() && target.object != NULL) {
			dist = p2p_dist(start, coll);
			return (uobject *)target.object;
		}
	}
	dist = 0.0;
	return NULL;
}


void destroy_sobj(s_object const &target) {

	//point coll_pos; get_point_of_collision(target, ship_pos, coll_pos);
	int const etype(ETYPE_ANIM_FIRE);

	switch (target.type) {
	case UTYPE_STAR:
		target.get_star().def_explode(u_exp_size[target.type], etype, signed_rand_vector());
		target.get_galaxy().calc_color(); // recompute color (what about when galaxies overlap?)
		break;
	case UTYPE_PLANET:
		target.get_planet().def_explode(u_exp_size[target.type], etype, signed_rand_vector());
		break;
	case UTYPE_MOON:
		target.get_moon().def_explode(u_exp_size[target.type], etype, signed_rand_vector());
		break;
	}
	target.register_destroyed_sobj();
}


// *** DRAW CODE ***


float get_light_scale(unsigned light) {return (glIsEnabled(light) ? 1.0 : 0.0);}


class ushader_group {
	shader_t planet_shader, star_shader, ring_shader, cloud_shader, atmospheric_shader;

	void set_light_scale(shader_t const &shader) const {
		vector2d const light_scale(get_light_scale(GL_LIGHT0), get_light_scale(GL_LIGHT1));
		shader.add_uniform_vector2d("light_scale", light_scale);
	}
	void shared_shader_setup(shader_t &shader, const char *texture_name="tex0") const {
		shader.begin_shader();
		shader.add_uniform_int(texture_name, 0);
		shader.setup_fog_scale(); // fog not needed?
	}

public:
	bool enable_planet_shader() {
		if (disable_shaders) return 0;

		if (!planet_shader.is_setup()) {
			planet_shader.set_vert_shader("per_pixel_lighting");
			planet_shader.set_frag_shader("linear_fog.part+ads_lighting.part*+planet_draw");
			shared_shader_setup(planet_shader);
		}
		planet_shader.enable();
		set_light_scale(planet_shader);
		set_specular(1.0, 80.0);
		return 1;
	}
	void disable_planet_shader() {
		if (planet_shader.is_setup()) {
			planet_shader.disable();
			set_specular(0.0, 1.0);
		}
	}

	bool enable_star_shader(colorRGBA const &colorA, colorRGBA const &colorB) { // no lighting
		if (disable_shaders) return 0;

		if (!star_shader.is_setup()) {
			set_active_texture(1);
			bind_3d_texture(get_noise_tex_3d());
			set_active_texture(0);
			star_shader.set_prefix("#define NUM_OCTAVES 8", 1); // FS
			star_shader.set_vert_shader("star_draw");
			star_shader.set_frag_shader("perlin_clouds_3d.part*+star_draw");
			star_shader.begin_shader();
			star_shader.add_uniform_int("tex0", 0);
			star_shader.add_uniform_int("cloud_noise_tex", 1);
			star_shader.add_uniform_float("time", 5.0E-5*frame_counter);
			star_shader.add_uniform_float("noise_scale", 3.5);
		}
		star_shader.enable();
		star_shader.add_uniform_color("colorA", colorA);
		star_shader.add_uniform_color("colorB", colorB);
		return 1;
	}
	void disable_star_shader() {
		if (star_shader.is_setup()) {star_shader.disable();}
	}

	bool enable_ring_shader(point const &planet_pos, float planet_radius, float ring_ri, float ring_ro, point const &sun_pos, float sun_radius) {
		if (disable_shaders) return 0;
		
		if (!ring_shader.is_setup()) {
			ring_shader.set_vert_shader("planet_rings");
			ring_shader.set_frag_shader("linear_fog.part+ads_lighting.part*+planet_rings");
			shared_shader_setup(ring_shader, "ring_tex");
			ring_shader.add_uniform_int("noise_tex", 1);
		}
		set_active_texture(1);
		select_texture(NOISE_GEN_MIPMAP_TEX, 0);
		set_active_texture(0);
		ring_shader.enable();
		ring_shader.add_uniform_vector3d("planet_pos", planet_pos);
		ring_shader.add_uniform_float("planet_radius", planet_radius);
		ring_shader.add_uniform_float("ring_ri",       ring_ri);
		ring_shader.add_uniform_float("ring_ro",       ring_ro);
		ring_shader.add_uniform_vector3d("sun_pos",    sun_pos);
		ring_shader.add_uniform_float("sun_radius",    sun_radius);
		upload_mvm_to_shader(ring_shader, "world_space_mvm");
		set_specular(0.5, 50.0);
		return 1;
	}
	void disable_ring_shader() {
		if (ring_shader.is_setup()) {
			set_specular(0.0, 1.0);
			ring_shader.disable();
		}
	}

	bool enable_cloud_shader(int cloud_tex_repeat) {
		if (disable_shaders) return 0;
		
		if (!cloud_shader.is_setup()) {
			bind_3d_texture(get_noise_tex_3d());
			cloud_shader.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
			cloud_shader.set_prefix("#define NUM_OCTAVES 8",    1); // FS
			cloud_shader.set_vert_shader("planet_clouds");
			cloud_shader.set_frag_shader("linear_fog.part+ads_lighting.part*+perlin_clouds_3d.part*+planet_clouds");
			shared_shader_setup(cloud_shader, "cloud_noise_tex");
			cloud_shader.add_uniform_float("time", 4.0E-6*frame_counter);
		}
		cloud_shader.enable();
		set_light_scale(cloud_shader);
		cloud_shader.add_uniform_float("noise_scale", 2.5*cloud_tex_repeat);
		return 1;
	}
	void disable_cloud_shader() {
		if (cloud_shader.is_setup()) {cloud_shader.disable();}
	}

	bool enable_atmospheric_shader(float atmosphere) {
		if (disable_shaders) return 0;
		
		if (!atmospheric_shader.is_setup()) {
			atmospheric_shader.set_prefix("#define USE_LIGHT_COLORS", 1); // FS
			atmospheric_shader.set_vert_shader("planet_clouds");
			atmospheric_shader.set_frag_shader("linear_fog.part+ads_lighting.part*+atmosphere");
			//shared_shader_setup(atmospheric_shader, "tex0");
			atmospheric_shader.begin_shader();
			atmospheric_shader.setup_fog_scale();
		}
		atmospheric_shader.enable();
		atmospheric_shader.add_uniform_float("atmosphere", atmosphere);
		set_light_scale(atmospheric_shader);
		return 1;
	}
	void disable_atmospheric_shader() {
		if (atmospheric_shader.is_setup()) {atmospheric_shader.disable();}
	}
};


void universe_t::draw_all_cells(s_object const &clobj, bool skip_closest, bool no_move, bool no_distant) {

	//RESET_TIME;
	unpack_color(water_c, P_WATER_C); // recalculate every time
	unpack_color(ice_c,   P_ICE_C  );
	int cxyz[3];
	point const &camera(get_player_pos());
	vector3d const motion_vector(camera, last_camera);
	float const speed(motion_vector.mag());
	camera_mv_speed const cmvs(camera, motion_vector, speed);
	ushader_group usg;
	last_camera = camera;
	glDisable(GL_LIGHTING);

	if (!no_distant || clobj.type < UTYPE_SYSTEM) { // drawing pass 1
		for (cxyz[2] = 0; cxyz[2] < int(U_BLOCKS); ++cxyz[2]) { // z
			for (cxyz[1] = 0; cxyz[1] < int(U_BLOCKS); ++cxyz[1]) { // y
				for (cxyz[0] = 0; cxyz[0] < int(U_BLOCKS); ++cxyz[0]) { // x
					draw_cell(cxyz, cmvs, usg, clobj, 0, no_move);
				}
			}
		}
	}
	if (clobj.type >= UTYPE_SYSTEM) { // in a system
		unsigned const npasses((skip_closest && (clobj.type == UTYPE_PLANET || clobj.type == UTYPE_MOON)) ? 2 : 3);

		for (unsigned pass = 1; pass < npasses; ++pass) { // drawing passes 2-3
			draw_cell(clobj.cellxyz, cmvs, usg, clobj, pass, no_move);
		}
	}
	universe_pld.draw_and_clear();
	glEnable(GL_LIGHTING);
	//PRINT_TIME("Draw Cells");
	//exit(0);
}


void set_current_system_light(s_object const &clobj, point const &pspos, float a_scale, float d_scale) {

	if (clobj.type >= UTYPE_SYSTEM && clobj.get_star().is_ok()) { // in a system
		ustar const &sun(clobj.get_star());
		point const pos(clobj.get_ucell().rel_center + sun.pos);
		set_sun_loc_color(pos, sun.color, sun.radius, 0, 0, a_scale, d_scale); // ending light from current system
		univ_sun_pos = pos;
		univ_sun_rad = sun.radius;
		univ_temp    = sun.get_energy()/max(TOLERANCE, p2p_dist_sq(sun.pos, pspos));
		have_sun     = 1;

		if (clobj.type == UTYPE_PLANET) { // planet (what about moon?)
			// calculate temperature based on distance of pspos to the poles
			uplanet const &planet(clobj.get_planet());
			vector3d const dir(pspos - planet.pos);
			float const pole_align(fabs(dot_product(planet.rot_axis, dir)/dir.mag()));
			// +/- 2 depending on distance to pole (should this be a percentage instead?)
			univ_temp += 2.0 - 4.0*pole_align;
			univ_temp  = max(0.0f, univ_temp);
		}
	}
	else {
		//set_light_galaxy_ambient_only();
		set_lighting_params(); // ending light default
		univ_sun_pos = point(0.0, 0.0, 1.0);
		univ_sun_rad = AVG_STAR_SIZE;
		univ_temp    = 0.0;
		sun_color    = SUN_LT_C;
		have_sun     = 0;
	}
	univ_temp = temp_to_deg_c(univ_temp);
}


bool is_shadowed(point const &pos, float radius, bool expand, ussystem const &sol, uobject const *&sobj) {

	float const exp(expand ? radius : 0.0);
	vector3d const dir(sol.sun.pos, pos);
	point const pos2(pos + dir*(1.1*radius)); // extend beyond the object radius so that we don't get a collision with it
	point const spos(sol.sun.pos - dir*(1.1*sol.sun.radius));

	for (unsigned k = 0; k < sol.planets.size(); ++k) {
		uplanet const &planet(sol.planets[k]);
		
		if (line_sphere_int_cont(pos2, spos, planet.pos, (planet.radius + exp))) {
			sobj = &planet;
			return 1; // shadowed by planet
		}
		if (planet.moons.empty()) continue;
		if (!line_sphere_int_cont(pos2, spos, planet.pos, (planet.mosize + exp))) continue; // not shadowed by planet or its moons

		for (unsigned l = 0; l < planet.moons.size(); ++l) {
			if (line_sphere_int_cont(pos2, spos, planet.moons[l].pos, (planet.moons[l].radius + exp))) {
				sobj = &planet.moons[l];
				return 1; // shadowed by moon
			}
		}
	}
	return 0;
}


void atten_color(colorRGBA &color, point const &p1, point const &p2, float radius, float expand) {

	assert(radius > 0.0 && expand > 0.0 && expand != 1.0);
	float const dist(p2p_dist(p1, p2));
	assert(!is_nan(dist) && dist > 0.0);
	//assert(dist >= radius && dist <= expand*radius); // might be too strong (failed after running for 10 hours)
	float const atten(CLIP_TO_01((expand - dist/radius)/(expand - 1.0f)));
	color *= atten;
}


bool get_universe_sun_pos(point const &pos, point &spos) {

	s_object result;

	if (universe.get_close_system(pos, result, 1.0)) {
		ustar const &star(result.get_star());
		
		if (star.is_ok()) {
			spos = star.pos;
			return 1;
		}
	}
	return 0;
}


void clear_ambient_color() {

	float const uambient[4] = {0.0, 0.0, 0.0, 0.0};
	glLightfv(get_universe_ambient_light(), GL_AMBIENT, uambient); // light is not enabled or disabled
}


// return values: -1: no sun, 0: not shadowed, 1: < half shadowed, 2: > half shadowed, 3: fully shadowed
// caches sobj and only determines if there is a shadow if NULL
int set_uobj_color(point const &pos, float radius, bool known_shadowed, int shadow_thresh, point &sun_pos,
				   uobject const *&sobj, bool no_ambient) // based on current star and current galaxy
{
	assert(radius < CELL_SIZE);
	float const expand(2.0);
	s_object result;
	bool const found(universe.get_close_system(pos, result, 1.0) != 0);
	bool const blend(!found && universe.get_close_system(pos, result, expand));
	ussystem *sol(NULL);
	point pos2(pos);
	offset_pos(pos2);

	if (result.galaxy >= 0) { // close to a galaxy
		pos2 -= result.get_ucell().pos;
		colorRGBA color;

		if (result.system >= 0) { // close to a solar system
			sol   = &result.get_system();
			color = sol->get_galaxy_color(); // non-const
		}
		else {
			ugalaxy const &galaxy(result.get_galaxy()); // what about overlapping galaxies?
			color = galaxy.color;
			//atten_color(color, pos2, galaxy.pos, (galaxy.radius + MAX_SYSTEM_EXTENT), expand);
		}
		if (no_ambient) {
			clear_ambient_color();
		}
		else {
			set_ambient_color(color*3.2); // hack to increase light on ships (isn't ambient reset below in some cases?)
		}
	}
	else { // not near any galaxies
		clear_ambient_color();
	}
	if (!found && !blend) { // no sun
		set_light_galaxy_ambient_only();
		return -1;
	}
	assert(sol);
	int shadow_val(0);
	ustar const &sun(sol->sun);

	if (!sun.is_ok()) { // no sun
		set_light_galaxy_ambient_only();
		return -1;
	}
	point const spos(result.get_ucell().rel_center + sun.pos);
	sun_pos = sun.pos; // spos?
	colorRGBA color(sun.color);
	if (blend) atten_color(color, pos2, sun.pos, (sol->radius + MAX_PLANET_EXTENT), expand);
	
	if (sun.is_ok() && (sobj != NULL || is_shadowed(pos2, radius, 1, *sol, sobj))) { // check for planet/moon shadowing
		++shadow_val;
		assert(sobj != NULL);
		float const sr(sobj->get_radius());

		if (line_sphere_int_cont(pos2, sun.pos, sobj->get_pos(), sr)) {
			++shadow_val;
			if (sr > radius && line_sphere_int_cont(pos2, sun.pos, sobj->get_pos(), (sr - radius))) ++shadow_val;
		}
	}
	set_sun_loc_color(spos, color, sun.radius, (known_shadowed || shadow_val > shadow_thresh), no_ambient, BASE_AMBIENT, BASE_DIFFUSE); // check for sun
	return shadow_val;
}


uobject const *get_shadowing_object(uobject const &obj, ustar const &sun) { // unused

	float dist(0.0); // unused
	vector3d const sun_dir(sun.pos - obj.pos);
	float const sun_dist(sun_dir.mag()), length(sun_dist - 1.5*(obj.radius + sun.radius));
	point const start(obj.pos + sun_dir*(1.5*obj.radius));
	uobject *res(line_intersect_universe(start, sun_dir, length, obj.radius, dist));
	assert(res != &obj && res != &sun);
	return res;
}


void universe_t::draw_cell(int const cxyz[3], camera_mv_speed const &cmvs, ushader_group &usg, s_object const &clobj, unsigned pass, bool no_move) {

	ucell &cell(get_cell(cxyz));
	if (!univ_sphere_vis(cell.rel_center, CELL_SPHERE_RAD)) return; // could use player_pdu.cube_visible()
	UNROLL_3X(current.cellxyz[i_] = cxyz[i_] + uxyz[i_];)
	draw_cell_contents(cell, cmvs, usg, clobj, pass, no_move);
}


// pass:
//  0. Draw all except for the player's system
//  If player's system is none then stop
//  1. Draw the player's system except for the player's sobj
//  2. Draw the player's sobj
void universe_t::draw_cell_contents(ucell &cell, camera_mv_speed const &cmvs, ushader_group &usg, s_object const &clobj, unsigned pass, bool no_move) {

	if (cell.galaxies == NULL) return; // galaxies not yet allocated
	point const &camera(cmvs.camera);
	assert(!is_nan(camera));
	bool const p_system(clobj.type >= UTYPE_SYSTEM);
	point_d const pos(cell.rel_center);
	float const wwsq((float)window_width*(float)window_width);

	for (unsigned i = 0; i < cell.galaxies->size(); ++i) { // remember, galaxies can overlap
		bool const sel_g((int)i == clobj.galaxy), sclip(!sel_g || (display_mode & 0x01) != 0);
		if (p_system && pass > 0 && !sel_g) continue; // drawn in previous pass
		ugalaxy &galaxy((*cell.galaxies)[i]);
		point_d const pos2(pos + galaxy.pos);
		if (calc_sphere_size(pos2, camera, STAR_MAX_SIZE, -galaxy.radius) < 0.18) continue; // too far away
		if (!univ_sphere_vis(pos2, galaxy.radius)) continue; // conservative, since galaxies are not spherical
		current.galaxy = i;
		galaxy.process(cell);
		float max_size(0.0);
		bool visible(0);

		for (unsigned c = 0; c < galaxy.clusters.size() && !visible; ++c) { // is there a better way to use clusters?
			point const pos_(pos + galaxy.clusters[c].center);
			visible  = univ_sphere_vis(pos_, galaxy.clusters[c].bounds);
			max_size = max(max_size, calc_sphere_size(pos_, camera, STAR_MAX_SIZE));
		}
		if (!visible)        continue; // galaxy is not visible
		if (max_size < 0.18) continue; // too small (normally will be freed earlier)
		set_ambient_color(galaxy.color);

		for (unsigned j = 0; j < galaxy.sols.size(); ++j) {
			bool const sel_s(sel_g && (int)j == clobj.system);
			if (p_system && ((pass == 0 && sel_s) || (pass != 0 && !sel_s))) continue; // draw in another pass
			ussystem &sol(galaxy.sols[j]);
			float const sradius(sol.sun.radius);

			if (max_size*sradius < 0.1*STAR_MAX_SIZE) {
				if (!sel_g) sol.free();
				continue;
			}
			point_d const pos3(pos + sol.pos);
			float const sizes(calc_sphere_size(pos3, camera, sradius));

			if (sizes < 0.1) {
				if (!sel_g) sol.free();
				continue;
			}
			bool const update_pass(sel_g && !no_move && ((int(tfticks)+j)&31) == 0);
			if (!update_pass && sol.gen != 0 && !univ_sphere_vis(pos3, sol.radius)) continue;
			bool const sel_sun(sel_s && (clobj.type == UTYPE_STAR || clobj.type == UTYPE_SYSTEM));
			bool const skip_s(p_system && ((pass == 1 && sel_sun) || (pass == 2 && !sel_sun)));
			bool const has_sun(sol.sun.is_ok());
			current.system  = j;
			current.cluster = sol.cluster_id;

			if (has_sun && !skip_s && univ_sphere_vis(pos3, 2.0*sradius)) {
				if (!sol.sun.draw(pos3, cmvs, usg, 1.0)) continue;
			}
			bool const planets_visible(PLANET_MAX_SIZE*sizes >= 0.3*sradius || !sclip);
			if (planets_visible) sol.process();
			if (sol.planets.empty()) continue;

			if (planets_visible) { // calc_sphere_size(pos3, camera, PLANET_MAX_SIZE)
				if (!has_sun) { // sun is gone
					set_light_galaxy_ambient_only();
				}
				else {
					set_sun_loc_color(pos3, sol.sun.color, sradius, 0, 0, BASE_AMBIENT, BASE_DIFFUSE); // slow - can this be made faster?
				}
				set_ambient_color(sol.get_galaxy_color());
			}
			else { // we know all planets are too far away
				if (!sel_g) sol.free_planets(); // optional
				if (!update_pass) continue;
			}
			for (unsigned k = 0; k < sol.planets.size(); ++k) {
				bool const sel_p(sel_s && (clobj.type == UTYPE_PLANET || clobj.type == UTYPE_MOON) && (int)k == clobj.planet);
				if (p_system && (pass == 2 && !sel_p)) continue; // draw in another pass
				uplanet &planet(sol.planets[k]);
				point_d const pos4(pos + planet.pos);
				float const pradius(PLANET_ATM_RSCALE*planet.radius), sizep(calc_sphere_size(pos4, camera, pradius));
				bool skip_draw(!planets_visible);

				if (sclip && sizep < (planet.ring_data.empty() ? 0.6 : 0.3)) {
					if (!sel_g && sizep < 0.3) planet.free();
					if (update_pass) {skip_draw = 1;} else {continue;}
				}
				current.planet = k;
				bool const sel_planet(sel_p && clobj.type == UTYPE_PLANET);
				bool const skip_p(p_system && ((pass == 1 && sel_planet) || (pass == 2 && !sel_planet)));
				if (!skip_p && !no_move) planet.do_update(sol.pos, has_sun, has_sun); // don't really have to do this if planet is not visible
				if (skip_draw || (planet.gen != 0 && !univ_sphere_vis(pos4, planet.mosize))) continue;
				bool const planet_visible((sizep >= 0.6 || !sclip) && univ_sphere_vis(pos4, pradius));
				
				if (planet_visible && planet.is_ok() && !skip_p) {
					//uobject const *sobj(sel_s ? get_shadowing_object(planet, sol.sun) : NULL);
					planet.check_gen_texture(int(sizep));
					float const rscale((planet.atmos > 0.5 && planet.tsize <= PLANET_ATM_TEX_SZ) ? PLANET_ATM_RSCALE : 1.0);
					planet.draw(pos4, cmvs, usg, rscale); // ignore return value?
				}
				planet.process();
				bool const skip_moons(p_system && sel_planet && !skip_p);

				if (sizep >= 1.0 && !skip_moons) {
					for (unsigned l = 0; l < planet.moons.size(); ++l) {
						bool const sel_m(sel_p && clobj.type == UTYPE_MOON && (int)l == clobj.moon);
						if (p_system && ((pass == 1 && sel_m) || (pass == 2 && !sel_m))) continue; // draw in another pass
						umoon &moon(planet.moons[l]);
						if (!moon.is_ok()) continue;
						point_d const pos5(pos + moon.pos);
						float const sizem(calc_sphere_size(pos5, camera, moon.radius));
						if ((sizem < 0.2 && sclip) || !univ_sphere_vis(pos5, moon.radius)) continue;
						current.moon = l;
						//uobject const *sobj(sel_s ? get_shadowing_object(moon, sol.sun) : NULL);
						moon.check_gen_texture(int(sizem));
						moon.draw(pos5, cmvs, usg, 1.0); // draw moon
					}
				}
				if (planet.is_ok() && !skip_p) {
					if (sizep > 0.1) {
						planet.ensure_rings_texture();
						planet.draw_prings(usg, pos4, sizep, pos3, (has_sun ? sradius : 0.0));
					}
					if (planet_visible && planet.atmos > 0.01 && planet.tsize > PLANET_ATM_TEX_SZ) {
						planet.draw_atmosphere(usg, pos4, sizep);
					}
				}
			} // planet k
		} // system j
	} // galaxy i
}


// *** UPDATE CODE ***


void universe_t::init() {

	//RESET_TIME;
	assert(U_BLOCKS & 1); // U_BLOCKS is odd

	for (unsigned i = 0; i < U_BLOCKS; ++i) { // z
		for (unsigned j = 0; j < U_BLOCKS; ++j) { // y
			for (unsigned k = 0; k < U_BLOCKS; ++k) { // x
				int const ii[3] = {k, j, i};
				cells[i][j][k].gen_cell(ii);
			}
		}
	}
	//PRINT_TIME("Init Cells");
}


void universe_t::shift_cells(int dx, int dy, int dz) {

	//RESET_TIME;
	assert((abs(dx) + abs(dy) + abs(dz)) == 1);
	vector3d const vxyz((float)dx, (float)dy, (float)dz);
	cell_block temp;

	for (unsigned i = 0; i < U_BLOCKS; ++i) { // z
		for (unsigned j = 0; j < U_BLOCKS; ++j) { // y
			for (unsigned k = 0; k < U_BLOCKS; ++k) { // x
				int const i2(i + dz), j2(j + dy), k2(k + dx);
				bool const zout(i2 < 0 || i2 >= int(U_BLOCKS));
				bool const yout(j2 < 0 || j2 >= int(U_BLOCKS));
				bool const xout(k2 < 0 || k2 >= int(U_BLOCKS));

				if (xout || yout || zout) { // allocate new cell
					int const ii[3]     = {k, j, i};
					temp.cells[i][j][k].gen = 0;
					temp.cells[i][j][k].gen_cell(ii);
				}
				else {
					cells[i2][j2][k2].gen           = 1;
					temp.cells[i][j][k]             = cells[i2][j2][k2];
					temp.cells[i][j][k].rel_center -= vxyz*CELL_SIZE;
				}
			}
		}
	}
	for (unsigned i = 0; i < U_BLOCKS; ++i) { // z
		for (unsigned j = 0; j < U_BLOCKS; ++j) { // y
			for (unsigned k = 0; k < U_BLOCKS; ++k) { // x
				if (cells[i][j][k].gen == 0) cells[i][j][k].free();
				cells[i][j][k]     = temp.cells[i][j][k];
				cells[i][j][k].gen = 0;
			}
		}
	}
	//PRINT_TIME("Shift Cells");
}


// *** GENERATION CODE ***
inline int gen_rand_seed1(point const &center) {

	return 196613*(int(RS_SCALE*center.x+0.5)) +
		   393241*(int(RS_SCALE*center.y+0.5)) +
		   786433*(int(RS_SCALE*center.z+0.5)) + RAND_CONST*123;
}


inline int gen_rand_seed2(point const &center) {

	return 6291469*(int(RS_SCALE*center.x+0.5)) +
		   3145739*(int(RS_SCALE*center.y+0.5)) +
		   1572869*(int(RS_SCALE*center.z+0.5)) + RAND_CONST*456;
}


void ucell::gen_cell(int const ii[3]) {

	if (gen) return; // already generated
	UNROLL_3X(rel_center[i_] = CELL_SIZE*(float(ii[i_] - (int)U_BLOCKSo2));)
	pos    = rel_center + get_scaled_upt();
	radius = 0.5*CELL_SIZE;
	set_rand2_state(gen_rand_seed1(pos), gen_rand_seed2(pos));
	get_rseeds();
	gen      = 0;
	galaxies = new vector<ugalaxy>;
	galaxies->resize(rand2()%MAX_GALAXIES_PER_CELL + 1);

	for (unsigned l = 0; l < galaxies->size(); ++l) {
		if (!(*galaxies)[l].create(*this, l)) { // can't place the galaxy
			galaxies->resize(l); // so remove it
			break;
		}
	}
}


// fix fp resolution/shifting?
bool ugalaxy::create(ucell const &cell, int index) {

	current.type = UTYPE_GALAXY;
	gen_rseeds();
	clear_systems();
	gen      = 0;
	radius   = rand_uniform2(GALAXY_MIN_SIZE, GALAXY_MAX_SIZE);
	xy_angle = rand_uniform2(0.0, TWO_PI);
	axis     = signed_rand_vector2_norm();
	scale    = vector3d(1.0, rand_uniform2(0.6, 1.0), rand_uniform2(0.07, 0.2));
	lrq_rad  = 0.0;
	lrq_pos  = all_zeros;
	gen_name(current);
	float const dbase(SYSTEM_MIN_SPACING + radius);
	float d[3][2];
	point galaxy_ext(all_zeros), pts[8];

	for (unsigned j = 0; j < 3; ++j) {
		d[j][0] = -radius*scale[j];
		d[j][1] =  radius*scale[j];
	}
	get_cube_points(d, pts);

	for (unsigned p = 0; p < 8; ++p) {
		rotate_vector3d(axis, -xy_angle, pts[p]);

		for (unsigned j = 0; j < 3; ++j) {
			galaxy_ext[j] = max(galaxy_ext[j], fabs(pts[p][j]));
		}
	}
	for (unsigned j = 0; j < 3; ++j) {
		galaxy_ext[j] = (CELL_SIZEo2 - MAX_SYSTEM_EXTENT - min(GALAXY_OVERLAP*radius, galaxy_ext[j]));
		assert(galaxy_ext[j] >= 0.0);
	}
	for (unsigned i = 0; i < MAX_TRIES; ++i) {
		for (unsigned j = 0; j < 3; ++j) {
			pos[j] = galaxy_ext[j]*signed_rand_float2();
		}
		bool too_close(0);

		for (int j = 0; j < index && !too_close; ++j) {
			too_close = is_close_to((*cell.galaxies)[j], GALAXY_OVERLAP);
		}
		if (!too_close) return 1;
	}
	return 0;
}


void ugalaxy::apply_scale_transform(point &pos_) const {

	pos_[0] *= scale[0];
	pos_[1] *= scale[1];
	pos_[2] *= scale[2];
	rotate_vector3d(axis, xy_angle, pos_);
}


float ugalaxy::get_radius_at(point const &pos_) const {

	if (lrq_rad > 0.0 && p2p_dist_sq(pos_, lrq_pos) < 0.000001*min(radius*radius, p2p_dist_sq(pos_, pos))) {
		return 1.001*lrq_rad; // point is so close to last point that we can used the cached value (be conservative)
	}
	vector3d dir(pos_);
	rotate_vector3d(axis, -xy_angle, dir);
	dir[0] *= scale[0];
	dir[1] *= scale[1];
	dir[2] *= scale[2];
	float const rval(radius*dir.mag());
	lrq_rad = rval;
	lrq_pos = pos_;
	return rval;
}


bool ugalaxy::is_close_to(ugalaxy const &g, float overlap_amount) const {

	vector3d const delta(pos - g.pos);
	float const dist(delta.mag());
	return (dist < ((overlap_amount/dist)*(get_radius_at(delta*-1.0) + g.get_radius_at(delta)) + SYSTEM_MIN_SPACING));
}


void ugalaxy::calc_color() {

	color.assign(0.0, 0.0, 0.0, 1.0);

	for (unsigned j = 0; j < sols.size(); ++j) {
		color += sols[j].sun.get_ambient_color_val();
		//sols[j].calc_color(); // delay until later
	}
	color *= 1.5/MAX_SYSTEMS_PER_GALAXY;
	color.set_valid_color();
}


void ugalaxy::calc_bounding_sphere() {

	float dist_sq_max(0.0);

	for (unsigned i = 0; i < sols.size(); ++i) {
		dist_sq_max = max(p2p_dist_sq(pos, sols[i].pos), dist_sq_max);
	}
	radius = sqrt(dist_sq_max);
}


void ugalaxy::process(ucell const &cell) {

	if (gen) return;
	//RESET_TIME;
	current.type = UTYPE_GALAXY;
	set_rseeds();
	unsigned num_systems(rand2()%(MAX_SYSTEMS_PER_GALAXY+1));
	vector<point> placed;
	vector<ugalaxy> const &galaxies(*cell.galaxies);

	for (unsigned i = 0; i < cell.galaxies->size(); ++i) { // find galaxies that overlap this one
		ugalaxy const &g((*cell.galaxies)[i]);
		if (&g == this || !is_close_to(g, 1.0)) continue;

		for (unsigned j = 0; j < g.sols.size(); ++j) { // find systems in other galaxies that overlap this one
			point const spos(g.pos + g.sols[j].pos);
			vector3d const sdelta(spos - pos);
			float const sdist(sdelta.mag());
			
			if (sdist < TOLERANCE || (sdist < (radius/sdist + MAX_SYSTEM_EXTENT) &&
				sdist < (get_radius_at(sdelta)/sdist + MAX_SYSTEM_EXTENT)))
			{
				placed.push_back(spos);
			}
		}
	}
	for (unsigned i = 0; i < num_systems; ++i) {
		point pos2;
		if (!gen_system_loc(pos2, placed)) num_systems = i; // can't place it, give up
	}
	sols.resize(num_systems);
	unsigned tot_systems(0);

	for (unsigned c = 0, cur = 0; c < clusters.size(); ++c) { // calculate actual center and radius values for each cluster
		system_cluster &cl(clusters[c]);
		unsigned const nsystems((unsigned)cl.systems.size());
		cl.radius = 0.0;
		cl.bounds = 0.0;
		cl.center = all_zeros;
		cl.s1     = cur;
		current.cluster = c;
		
		for (unsigned i = 0; i < nsystems; ++i) {
			cl.center += cl.systems[i];
		}
		cl.center /= nsystems;

		for (unsigned i = 0; i < nsystems; ++i, ++cur) {
			cl.radius        = max(cl.radius, p2p_dist_sq(cl.center, cl.systems[i]));
			current.system   = cur;
			sols[cur].galaxy = this;
			sols[cur].cluster_id = c;
			sols[cur].create(cl.systems[i]);
		}
		cl.systems.clear();
		cl.radius    = sqrt(cl.radius);
		cl.bounds    = cl.radius + MAX_SYSTEM_EXTENT;
		tot_systems += nsystems;
		cl.s2        = cur;
	}
	assert(tot_systems == num_systems);
	calc_bounding_sphere();
	calc_color();
	lrq_rad = 0.0;
	gen     = 1;
	//if (num_systems > 480) PRINT_TIME("Galaxy Process");
}


bool ugalaxy::gen_system_loc(point &pos2, vector<point> const &placed) {

	pos2 = zero_vector;

	for (unsigned i = 0; i < MAX_TRIES; ++i) {
		float const rsize(radius*(1.0 - sqrt(rand2d())));
		pos2  = gen_rand_vector2(rsize);
		apply_scale_transform(pos2);
		pos2 += pos;
		bool bad_pos(0);
		
		for (unsigned j = 0; j < 3 && !bad_pos; ++j) {
			if (fabs(pos2[j]) > (CELL_SIZEo2 - MAX_SYSTEM_EXTENT)) bad_pos = 1; // shouldn't really get here
		}
		for (unsigned j = 0; j < placed.size() && !bad_pos; ++j) {
			bad_pos = dist_less_than(pos2, placed[j], SYSTEM_MIN_SPACING);
		}
		for (unsigned c = 0; c < clusters.size() && !bad_pos; ++c) {
			if (dist_less_than(pos2, clusters[c].center, clusters[c].bounds)) { // close to a cluster
				vector<point> const &cs(clusters[c].systems);

				for (unsigned s = 0; s < cs.size() && !bad_pos; ++s) {
					bad_pos = dist_less_than(pos2, cs[s], SYSTEM_MIN_SPACING);
				}
			}
		}
		if (bad_pos) continue;
		unsigned in_cluster((unsigned)clusters.size());
		float dmin(0.0);

		for (unsigned c = 0; c < clusters.size(); ++c) {
			float const test_dist((dmin == 0.0) ? clusters[c].radius : min(clusters[c].radius, dmin));

			if (dist_less_than(pos2, clusters[c].center, test_dist)) {
				in_cluster = c;
				dmin       = p2p_dist(pos2, clusters[c].center);
			}
		}
		if (in_cluster == clusters.size()) { // create initial clusters with fixed radius around a starting point
			float const cluster_size(0.1*radius + 0.3*p2p_dist(pos2, pos));
			clusters.push_back(system_cluster(cluster_size, pos2));
			clusters.back().systems.reserve(MAX_SYSTEMS_PER_GALAXY/10);
		}
		assert(in_cluster < clusters.size());
		system_cluster &cl(clusters[in_cluster]);
		cl.systems.push_back(pos2);

		if (cl.systems.size() == 2) { // use the center of the two points
			cl.center = (cl.systems[0] + cl.systems[1])*0.5;
			cl.bounds = 0.0;
		}
		cl.bounds = max(cl.bounds, (p2p_dist(pos2, cl.center) + SYSTEM_MIN_SPACING));
		return 1;
	}
	return 0;
}


void ussystem::create(point const &pos_) {

	current.type = UTYPE_SYSTEM;
	gen_rseeds();
	planets.clear();
	gen    = 0;
	radius = 0.0;
	pos    = pos_;
	galaxy_color.alpha = 0.0; // set to an invalid state
	sun.create(pos);
}


void ustar::create(point const &pos_) {

	current.type = UTYPE_STAR;
	set_defaults();
	gen_rseeds();
	rot_ang  = rot_ang0 = 0.0;
	pos      = pos_;
	temp     = rand_gaussian2(55.0, 10.0);
	radius   = 0.25*rand_uniform2(STAR_MIN_SIZE, STAR_MAX_SIZE) + (37.5*STAR_MAX_SIZE/temp)*rand_gaussian2(0.3, 0.1);
	radius   = max(radius, STAR_MIN_SIZE); // you get a lot of stars of size exactly STAR_MIN_SIZE
	gen_color();
	density  = rand_uniform2(3.0, 5.0);
	set_grav_mass();
	rot_axis = signed_rand_vector2_norm(); // for orbital plane
	gen      = 1; // get_name() is called later
	current.type = UTYPE_STAR;
	if (current.is_destroyed()) status = 1;
}


colorRGBA const &ussystem::get_galaxy_color() {

	if (galaxy_color.alpha == 0.0) calc_color(); // lazy update
	assert(galaxy_color.alpha > 0.0);
	return galaxy_color;
}


void ussystem::calc_color() { // delayed until the color is actually needed

	assert(galaxy);
	vector<ussystem> const &sols(galaxy->sols);
	assert(!sols.empty());
	galaxy_color.assign(0.0, 0.0, 0.0, 1.0);
	float const max_dist(5.0*SYSTEM_MIN_SPACING), max_dist_sq(max_dist*max_dist);
	float sum(0.0);

	for (unsigned j = 0; j < galaxy->clusters.size(); ++j) {
		ugalaxy::system_cluster const &cl(galaxy->clusters[j]);
		if (p2p_dist_sq(pos, cl.center) > (max_dist + cl.radius)*(max_dist + cl.radius)) continue; // too far away

		for (unsigned s = cl.s1; s < cl.s2; ++s) {
			if (&sols[s] == this)   continue; // skip ourselves
			float const dist(p2p_dist_sq(pos, sols[s].sun.pos));
			if (dist > max_dist_sq) continue; // too far away
			float const weight(1.0/dist); // 1/r^2 falloff
			sum          += weight;
			galaxy_color += sols[s].sun.get_ambient_color_val()*weight;
		}
	}
	if (sum == 0.0) return; // no other systems
	galaxy_color *= 1.5/((sum/sols.size())*MAX_SYSTEMS_PER_GALAXY); // remember to normalize
	galaxy_color.set_valid_color();
}


void ussystem::process() {

	if (gen) return;
	current.type = UTYPE_STAR;
	sun.set_rseeds();
	sun.gen_name(current);
	current.type = UTYPE_SYSTEM;
	set_rseeds();
	planets.resize((unsigned)sqrt(float((rand2()%(MAX_PLANETS_PER_SYSTEM+1))*(rand2()%(MAX_PLANETS_PER_SYSTEM+1)))));
	float const sradius(sun.radius);
	radius = sradius;

	for (unsigned i = 0; i < planets.size(); ++i) {
		current.planet    = i;
		planets[i].system = this;

		if (!planets[i].create_orbit(planets, i, pos, sun.rot_axis, sradius, PLANET_MAX_SIZE, PLANET_MIN_SIZE,
			INTER_PLANET_MIN_SPACING, PLANET_TO_SUN_MAX_SPACING, PLANET_TO_SUN_MIN_SPACING, 0.0))
		{ // failed to place planet
			planets.resize(i);
			remove_excess_cap(planets);
			break;
		}
		float const dmax(planets[i].orbit + planets[i].radius + MOON_TO_PLANET_MAX_SPACING + MOON_MAX_SIZE);
		radius = max(radius, dmax); // too bad we can't use p.mosize
	}
	radius = max(radius, 0.5f*(PLANET_TO_SUN_MIN_SPACING + PLANET_TO_SUN_MAX_SPACING)); // set min radius so that hyperspeed coll works
	gen    = 1;
}


bool urev_body::land_temp_ok() const {

	return (temp >= MIN_LAND_TEMP && temp <= MAX_LAND_TEMP);
}


bool urev_body::colonizable() const {

	return (is_ok() && temp >= MIN_COLONY_TEMP && temp <= MAX_COLONY_TEMP && colonizable_int());
}


bool urev_body::liveable() const { // only planets are liveable

	return (is_ok() && water > 0.15 && atmos > 0.25 && temp >= MIN_LIVE_TEMP && temp <= MAX_LIVE_TEMP);
}


void uplanet::calc_temperature() {

	assert(orbit > TOLERANCE);
	temp = system->sun.get_energy()/(orbit*orbit); // ~2-50
}


void uplanet::create(bool phase) {

	if (phase == 1) return; // no phase 1, only phase 0
	current.type = UTYPE_PLANET;
	gen_rotrev();
	mosize = radius;
	moons.clear();
	ring_data.clear();

	// atmosphere, water, temperature, gravity
	calc_temperature();
	density = rand_uniform2(0.8, 1.2);
	if (temp < CGAS_TEMP) density *= 0.5 + 0.5*(temp/CGAS_TEMP); // cold gas
	set_grav_mass();
	
	if (temp < FREEZE_TEMP) { // cold
		atmos = rand_uniform2(-0.2, 1.0);
		water = rand_uniform2(0.0, MAX_WATER); // ice
	}
	else if (temp > BOIL_TEMP) { // hot
		atmos = rand_uniform2(-0.9, 0.5);
		water = 0.0;
	}
	else if (temp > NO_AIR_TEMP) { // very hot
		atmos = 0.0;
		water = 0.0;
	}
	else { // average temp
		atmos = rand_uniform2(-0.3, 1.5);
		water = max(0.0f, min(MAX_WATER, 0.5f*(atmos + rand_uniform2(-MAX_WATER, 0.9*MAX_WATER))));
	}
	atmos     = CLIP_TO_01(atmos);
	float const rsc_scale(liveable() ? 2.0 : (colonizable() ? 1.0 : 0.5));
	resources = 750.0*radius*rsc_scale*(1.0 + 0.25*atmos - 0.25*fabs(0.5 - water))*(1.0 - fabs(1.0 - density));
	check_owner(current); // must be after setting of resources
	gen_color();
	gen_name(current);
	cloud_tex_repeat = (1 + (rand2()&1)); // 1 or 2
	current.type     = UTYPE_PLANET;
	if (current.is_destroyed()) status = 1;
}


float uplanet::get_vegetation() const {

	return ((has_vegetation() && temp > MIN_PLANT_TEMP && temp < MAX_PLANT_TEMP) ? sqrt(atmos*water) : 0.0);
}


void uplanet::process() {

	if (gen) return;
	current.type = UTYPE_PLANET;
	set_rseeds();
	if (temp < CGAS_TEMP && (rand2()&1)) {gen_prings();} // rings
	unsigned num_moons(0);

	if (rand2()&1) { // has moons
		num_moons = (unsigned)sqrt(float((rand2()%(MAX_MOONS_PER_PLANET+1))*(rand2()%(MAX_MOONS_PER_PLANET+1))));
	}
	moons.resize(num_moons);

	for (unsigned i = 0; i < moons.size(); ++i) {
		current.moon    = i;
		moons[i].planet = this;

		if (!moons[i].create_orbit(moons, i, pos, rot_axis, radius, MOON_MAX_SIZE, MOON_MIN_SIZE,
			INTER_MOON_MIN_SPACING, MOON_TO_PLANET_MAX_SPACING, MOON_TO_PLANET_MIN_SPACING, MOON_TO_PLANET_MIN_GAP))
		{ // failed to place moon
			moons.resize(i);
			remove_excess_cap(moons);
			break;
		}
		mosize = max(mosize, (radius + moons[i].orbit + moons[i].radius));
	}
	if (!moons.empty()) { // calculate rotation rate about rotation axis due to moons (Note: Stars can rotate as well.)
		// rk_term = r/(2*PI*a*k);
		// T^2 = k*(4*PI*PI*a*a*a/(G*(M + m)*cosf(i)*cosf(i)))*((m/M)*(r/R) + (M/m)*(D/d)*rk_term*rk_term);
		// a = average moon orbit, M = mass, m = sum of moon mass, R = radius, r = average moon radius
		// rot_rate = 360.0/t = 360.0/(3600*TICKS_PER_SECOND*T) = 1.0/(10.0*TICKS_PER_SECOND*sqrt(...))
		float rav(0.0), aav(0.0), dav(0.0), cav(0.0), mtot(0.0);

		for (unsigned i = 0; i < moons.size(); ++i) { // technically, this should be reclaculated if a moon is destroyed
			mtot += moons[i].mass;
			rav  += moons[i].radius*moons[i].mass;
			aav  += moons[i].orbit*moons[i].mass;
			dav  += moons[i].density*moons[i].mass;
			cav  += (1.0 - fabs(dot_product(rot_axis, moons[i].rev_axis)))*moons[i].mass; // probably not quite right
		}
		rav /= mtot;
		aav /= mtot;
		dav /= mtot;
		cav /= mtot;
		float const k(rand_uniform2(0.05, 0.5)), ci(cosf(cav)), rk_term(rav/(2*PI*aav*k));
		float const T_sq(k*(4*PI*PI*aav*aav*aav/(mass + mtot)*ci*ci)*((mtot/mass)*(rav/radius) + (mass/mtot)*(density/dav)*rk_term*rk_term));
		assert(T_sq > 0.0);
		rot_rate = ROT_RATE_CONST/(10.0*TICKS_PER_SECOND*sqrt(T_sq));
	}
	gen = 1;
}


point_d uplanet::do_update(point_d const &p0, bool update_rev, bool update_rot) {

	bool const has_sun(system->sun.is_ok());
	if (!has_sun) temp = 0.0;
	point_d const planet_pos(urev_body::do_update(p0, (has_sun && update_rev), (has_sun && update_rot)));
	bool const ok(is_ok());
	
	for (unsigned i = 0; i < moons.size(); ++i) {
		moons[i].do_update(planet_pos, ok, ok);
		if (!has_sun) moons[i].temp = 0.0;
	}
	return planet_pos;
}


void uplanet::update_shadows() {

	for (unsigned i = 0; i < moons.size(); ++i) {
		moons[i].update_shadows();
	}
}


struct upring {
	float radius1, radius2;
};


void uplanet::gen_prings() {

	unsigned const nr((rand2()%10)+1);
	float const sr(4.0/nr);
	float lastr(rand_uniform2(radius, 1.1*radius));
	vector<upring> rings(nr);

	for (unsigned i = 0; i < nr; ++i) {
		upring &ring(rings[i]);
		ring.radius1 = lastr        + sr*radius*rand_uniform2(0.0,  0.1);
		ring.radius2 = ring.radius1 + sr*radius*rand_uniform2(0.05, 0.3);
		lastr = ring.radius2;
	}
	ring_data.resize(RING_TEX_SZ);
	for (unsigned i = 0; i < RING_TEX_SZ; ++i) {ring_data[i].set_c4(ALPHA0);}
	ring_ri = rings.front().radius1;
	ring_ro = rings.back().radius2;
	float const rdiv((RING_TEX_SZ-3)/(ring_ro - ring_ri));
	colorRGBA rcolor(color);
	UNROLL_3X(rcolor[i_] += rand_uniform2(0.1, 0.6);)
	rcolor.alpha = rand_uniform2(0.4, 0.9);

	for (vector<upring>::const_iterator i = rings.begin(); i != rings.end(); ++i) {
		unsigned const tri(1+(i->radius1 - ring_ri)*rdiv), tro(1+(i->radius2 - ring_ri)*rdiv);
		assert(tri > 0 && tro+1 < RING_TEX_SZ && tri < tro);
		color_wrapper cw;
		UNROLL_3X(rcolor[i_] = CLIP_TO_01(rcolor[i_]*(1.0f + rand_uniform2(-0.15, 0.15)));)
		rcolor.A = CLIP_TO_01(rcolor.A*(1.0f + rand_uniform2(-0.1, 0.1)));
		cw.set_c4(rcolor);
		for (unsigned j = tri; j < tro; ++j) {ring_data[j] = cw;}
	}
	float max_rs(0.0);
	UNROLL_3X(rscale[i_] = rand_uniform2(0.9, 2.2);)
	UNROLL_3X(max_rs     = max(max_rs, rscale[i_]);)
	mosize = max(mosize, max_rs*lastr); // extend planet effective size
}


void uplanet::get_valid_orbit_r(float &orbit_r, float obj_r) const { // for satellites

	assert(is_ok());

	for (unsigned i = 0; i < moons.size(); ++i) {
		if (!moons[i].is_ok()) continue;
		float const orad(ORBIT_SPACE_MARGIN*(obj_r + moons[i].radius));

		if ((moons[i].orbit - orad) < orbit_r && (moons[i].orbit + orad) > orbit_r) { // orbits overlap (rarely happens)
			orbit_r = (moons[i].orbit - orad);
		}
	}
	orbit_r = max(orbit_r, (radius + ORBIT_SPACE_MARGIN*obj_r)); // apply min_orbit
}


void umoon::get_valid_orbit_r(float &orbit_r, float obj_r) const { // for satellites

	assert(is_ok());
	if (!planet || !planet->is_ok()) return;
	float const orad(ORBIT_SPACE_MARGIN*(obj_r + planet->radius));
	orbit_r = min(orbit_r, (orbit - orad));
	orbit_r = max(orbit_r, (radius + ORBIT_SPACE_MARGIN*obj_r)); // apply min_orbit
}


void umoon::calc_temperature() {

	temp = planet->system->sun.get_energy()/max(TOLERANCE, p2p_dist_sq(planet->system->sun.pos, pos));
	if (shadowed) temp *= 0.75; // cooler in shadow
}


void umoon::create(bool phase) { // no rotation due to satellites

	current.type = UTYPE_MOON;
	
	if (phase == 0) {
		gen_rotrev();
		gen = 2;
	}
	else {
		assert(gen == 2);
		density = rand_uniform2(0.8, 1.2);
		set_grav_mass();
		update_shadows();
		temp = planet->temp;
		gen_color();
		gen_name(current);
		resources = 750.0*radius*(colonizable() ? 2.0 : 1.0)*(1.0 - fabs(1.0 - density));
		check_owner(current); // must be after setting of resources
		calc_temperature(); // has to be after setting of resources - resources must be independent of moon position/temperature
		gen       = 1;
	}
	current.type = UTYPE_MOON;
	if (current.is_destroyed()) status = 1;
}


void urev_body::gen_rotrev() {

	set_defaults();
	gen_rseeds();
	tid = tsize = 0;
	rot_ang  = rev_ang = 0.0;
	rot_rate = 0.0;
	rev_rate = 0.0;
	rot_ang0 = 360.0*rand2d(); // degrees in OpenGL
	rev_ang0 = 360.0*rand2d(); // degrees in OpenGL
	rot_axis = signed_rand_vector2_norm();
	// inclination angle = angle between rot_axis and rev_axis

	// calculate revolution rate around parent
	// mu = G*M, a = orbit, R = radius, T = 2*PI*sqrt(a*a*a/mu)
	// T = 1.4*sqrt(a*a*a/R*R*R) hours (Standard planet), 1 hour = 3600*TICKS_PER_SECOND ticks
	// rev_rate = 360.0/t = 360.0/(3600*TICKS_PER_SECOND*T) = 1.0/(0.14*TICKS_PER_SECOND*sqrt(a*a*a/R*R*R))
	float const aoR(orbit/radius);
	rev_rate = REV_RATE_CONST/(0.14*TICKS_PER_SECOND*aoR*sqrt(aoR));
}


// Note: Update is only done when the objects solar system is visible to the player
point_d urev_body::do_update(point_d const &p0, bool update_rev, bool update_rot) {

	// what about fp errors when tfticks gets large?
	float const ra1(rev_ang);
	if ((animate2 || rot_ang == 0.0) && update_rot) rot_ang = rot_ang0 + double(tfticks)*double(rot_rate);
	if ((animate2 || rev_ang == 0.0) && update_rev) rev_ang = rev_ang0 + double(tfticks)*double(rev_rate);
	point_d new_pos(v_orbit);
	rotate_vector3d(vector3d_d(rev_axis), rev_ang/TO_DEG, new_pos); // more accurate (is this necessary?)
	new_pos *= double(orbit);
	new_pos += point_d(p0);
	pos      = new_pos;

	if (update_rev && int(10*ra1) != int(10*rev_ang)) {
		update_shadows(); // update every 0.1 degree
		calc_temperature();
	}
	return new_pos;
}


template<typename T>
bool urev_body::create_orbit(vector<T> const &objs, int i, point const &pos0, vector3d const &raxis, float radius0,
							 float max_size, float min_size, float rspacing, float ispacing, float minspacing, float min_gap)
{
	radius   = (min(0.4f*radius0, max_size) - min_size)*((float)rand2d()) + min_size;
	rev_axis = raxis;
	float const rad2(radius + rspacing), min_orbit(max((MIN_RAD_SPACE_FACTOR*(radius + radius0) + min_gap), minspacing));
	rev_axis += signed_rand_vector2_norm()*ORBIT_PLANE_DELTA;
	rev_axis.normalize();
	vector3d const start_vector(signed_rand_vector2_norm()); // doesn't matter, any will do
	cross_product(rev_axis, start_vector, v_orbit);
	v_orbit.normalize();
	bool too_close(1);
	unsigned counter;

	for (counter = 0; counter < MAX_TRIES && too_close; ++counter) {
		orbit     = rand_uniform2(min_orbit, ispacing);
		too_close = 0;

		for (int j = 0; j < i; ++j) { // slightly inefficient
			if (fabs(objs[j].orbit - orbit) < ORBIT_SPACE_MARGIN*(rad2 + objs[j].radius)) { // what about checking against objs[j].rad_spacing?
				too_close = 1; break;
			}
		}
	}
	if (too_close) return 0;
	create(0);
	do_update(pos0);
	create(1);
	return 1;
}


void urev_body::dec_orbiting_refs(s_object const &sobj) {

	assert(sobj.object == this);
	//assert(orbiting_refs > 0); // too strong - can fail if player leaves the galaxy and planets/moons are deleted (refs are reset)
	if (orbiting_refs >  0) --orbiting_refs;
	if (orbiting_refs == 0) set_owner(sobj, NO_OWNER);
}


// *** COLORS ***

void ustar::gen_color() {

	if (temp < 25.0) { // black: 0-25 (black hole)
		color.assign(0.0, 0.0, 0.0, 1.0);
	}
	else if (temp < 30.0) { // deep red: 25-30
		color.assign(0.2*(temp - 25.0), 0.0, 0.0, 1.0);
	}
	else if (temp < 40.0) { // red-orange-yellow: 30-40
		color.assign(1.0, 0.1*(temp - 30.0), 0.0, 1.0);
	}
	else if (temp < 65.0) { // yellow-white: 40-65
		color.assign(1.0, 1.0, 0.04*(temp - 40.0), 1.0);
	}
	else if (temp < 75.0) { // white-blue: 65-75
		color.assign((0.6 + 0.05*(75.0  - temp)), (0.8 + 0.025*(75.0 - temp)), 1.0, 1.0);
	}
	else { // blueish: > 75
		color.assign(0.6, 0.8, 1.0);
	}
	color.set_valid_color();
	gen_colorAB(0.8*MP_COLOR_VAR);
	if (temp < 30.0) colorA.G = colorA.B = colorB.G = colorB.B = 0.0; // make sure it's just red
}


inline colorRGBA ustar::get_ambient_color_val() const {

	return (is_ok() ? colorRGBA(color.R*STAR_MAX_SIZE_INV, color.G*STAR_MAX_SIZE_INV, color.B*STAR_MAX_SIZE_INV, color.A) : BLACK);
}


void uplanet::gen_color() {

	float const bright(rand_uniform2(0.5, 0.75));
	color.assign((0.75*bright + 0.40*rand2d()), (0.50*bright + 0.30*rand2d()), (0.25*bright + 0.15*rand2d()), 1.0);
	color.set_valid_color();
	
	if (has_vegetation()) { // override with Earth/terran colors (replace the above code?)
		colorA = colorRGBA(0.05, 0.35, 0.05, 1.0);
		colorB = colorRGBA(0.60, 0.45, 0.25, 1.0);
		adjust_colorAB(0.25*MP_COLOR_VAR);
		blend_color(color, colorA, colorB, 0.5, 0); // average the two colors
	}
	else {
		gen_colorAB(MP_COLOR_VAR);
	}
	if (water > 0.0) { // blend water with land for distant views?
		blend_color(color, ((temp < FREEZE_TEMP) ? P_ICE_C : P_WATER_C), color, water, 0);
	}
	if (atmos > 0.0) blend_color(color, CLOUD_C, color, 0.5*atmos, 0); // blend atmosphere with planet color for distant views
	color.set_valid_color();
}


void umoon::gen_color() {
	
	float const brightness(rand_uniform2(0.5, 0.75));
	float const mult(shadowed ? BASE_DIFFUSE/BASE_AMBIENT : 1.0);

	for (unsigned i = 0; i < 3; ++i) {
		color[i] = mult*(0.75*brightness + 0.25*rand2d());
	}
	color.alpha = 1.0;
	color.set_valid_color();
	gen_colorAB(1.4*MP_COLOR_VAR);
}


void uobj_solid::adjust_colorAB(float delta) {

	for (unsigned i = 0; i < 3; ++i) {
		float const d(delta*rand2d());
		colorA[i] += d;
		colorB[i] -= d;
	}
	colorA.set_valid_color();
	colorB.set_valid_color();
}


void uobj_solid::gen_colorAB(float delta) {

	colorA = colorB = color;
	adjust_colorAB(delta);
}


// *** TEXTURES ***


void urev_body::check_gen_texture(unsigned size) {

	if (size <= MIN_TEX_OBJ_SZ) return; // too small to be textured
	unsigned const tsize0(get_texture_size(size));
	bool do_gen(0);

	if (!glIsTexture(tid)) { // texture has not been generated
		gen_surface();
		do_gen = 1;
	}
	else if (tsize0 != tsize) { // new texture size
		::free_texture(tid); // delete old texture
		do_gen = 1;
	}
	if (do_gen) {
		tsize = tsize0; // new size
		create_texture(tsize0); // new texture
	}
}


void urev_body::create_texture(unsigned size) {

	assert(size <= MAX_TEXTURE_SIZE);
	vector<unsigned char> data(3*size*size);
	gen_texture_data(&data.front(), size, USE_HEIGHTMAP);
	setup_texture(tid, GL_MODULATE, 0, 0, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, size, size, 0, GL_RGB, GL_UNSIGNED_BYTE, &data.front());
}


void urev_body::get_surface_color(unsigned char *data, float val, float phi) const { // val in [0,1]

	bool const frozen(temp < FREEZE_TEMP);
	unsigned char const white[3] = {255, 255, 255};
	unsigned char const gray[3]  = {100, 100, 100};
	float const coldness(fabs(phi - PI_TWO)*2.0*PI_INV);

	if (val < water) { // underwater
		RGB_BLOCK_COPY(data, wic[frozen]);

		if (coldness > 0.8) { // ice
			float const blend_val(CLIP_TO_01(5.0f*(coldness - 0.8f) + 2.0f*(val - water)));
			BLEND_COLOR(data, white, data, blend_val);
		}
		return;
	}
	float const water_adj(0.07), val_ws((water > 0.0) ? wr_scale*(val - water) : val); // rescale for water

	if (water > 0.2 && atmos > 0.1) { // Earthlike planet
		if (val_ws < 0.1) { // low ground
			RGB_BLOCK_COPY(data, b);
		}
		else if (val_ws < 0.4) {
			BLEND_COLOR(data, a, b, 3.3333*(val_ws - 0.1));
		}
		else if (val_ws < 0.45) { // medium ground
			RGB_BLOCK_COPY(data, a);
		}
		else if (val_ws < 0.6) {
			BLEND_COLOR(data, gray, a, 6.6667*(val_ws - 0.45));
		}
		else { // high ground
			RGB_BLOCK_COPY(data, gray);
		}
	}
	else { // alien-like planet
		BLEND_COLOR(data, a, b, val_ws);
	}
	if (temp < BOIL_TEMP) {
		if (val < water + water_adj) { // close to water line (can have a little water even if water == 0)
			BLEND_COLOR(data, data, wic[frozen], (val - water)/water_adj);
			if (coldness > 0.9) {BLEND_COLOR(data, white, data, 10.0*(coldness - 0.9));} // ice
			return;
		}
		// phi=PI/2 => equator, phi=0.0 => north pole, phi=PI => south pole
		float const coldness(fabs(phi - PI_TWO)*2.0*PI_INV), st((1.0 - coldness*coldness)*snow_thresh);

		if (val > (st + 1.0E-6)) { // blend in some snow
			BLEND_COLOR(data, white, data, (val - st)/(1.0 - st));
			return;
		}
	}
}


void urev_body::calc_snow_thresh() {

	float const snow_temp(CLIP_TO_01(2.0f*((0.5f*FREEZE_TEMP + 0.5f*BOIL_TEMP) - temp))/(BOIL_TEMP - FREEZE_TEMP));
	float const snow_val(CLIP_TO_01(2.0f*(water - 0.1f))*snow_temp);
	snow_thresh = max(water, (1.0f - snow_val));
}


void uobj_solid::get_colors(unsigned char ca[3], unsigned char cb[3]) const {

	for (unsigned d = 0; d < 3; ++d) {
		assert(colorA[d] >= 0.0 && colorA[d] <= 1.0);
		assert(colorB[d] >= 0.0 && colorB[d] <= 1.0);
	}
	unpack_color(ca, colorA);
	unpack_color(cb, colorB);
}


unsigned get_texture_size(float psize) {

	unsigned i;
	unsigned const ps2(unsigned(2.0*psize));

	for (i = 8; i < MAX_TEXTURE_SIZE; i <<= 1) {
		if (ps2 < i) return i;
	}
	return i;
}


void free_universe_textures() {

	universe.free_textures();
}


void universe_t::free_textures() { // should be OK even if universe isn't setup

	for (unsigned z = 0; z < U_BLOCKS; ++z) { // z
		for (unsigned y = 0; y < U_BLOCKS; ++y) { // y
			for (unsigned x = 0; x < U_BLOCKS; ++x) { // x
				ucell &cell(cells[z][y][x]);
				if (cell.galaxies == NULL) continue;
				
				for (unsigned i = 0; i < cell.galaxies->size(); ++i) {
					ugalaxy &galaxy((*cell.galaxies)[i]);
					
					for (unsigned j = 0; j < galaxy.sols.size(); ++j) {
						ussystem &sol(galaxy.sols[j]);
						
						for (unsigned k = 0; k < sol.planets.size(); ++k) {
							uplanet &planet(sol.planets[k]);
							planet.free_texture();
							
							for (unsigned l = 0; l < planet.moons.size(); ++l) {
								planet.moons[l].free_texture();
							}
						}
					}
				}
			}
		}
	}
}


// *** MEMORY - FREE CODE ***

void ucell::free() {

	if (galaxies != NULL) {
		for (unsigned i = 0; i < galaxies->size(); ++i) (*galaxies)[i].free();
		delete galaxies;
		galaxies = NULL;
	}
	gen = 0;
}


void ugalaxy::clear_systems() {

	sols.clear();
	clusters.clear();
}


void ugalaxy::free() {

	for (unsigned i = 0; i < sols.size(); ++i) sols[i].free();
	clear_systems();
	gen = 0;
}


void ussystem::free_planets() {

	for (unsigned i = 0; i < planets.size(); ++i) planets[i].free();
}


void ussystem::free() {

	free_planets();
	planets.clear();
	sun.free();
	galaxy_color.alpha = 0.0; // set to an invalid state
	gen = 0;
}


void uplanet::free() {

	for (unsigned i = 0; i < moons.size(); ++i) moons[i].free();
	moons.clear();
	ring_data.clear();
	urev_body::free();
}


void urev_body::free() {

	if (gen) free_texture();
	delete surface; // OK if already NULL
	surface = NULL;
	//unset_owner();
	uobj_solid::free();
}


// *** DRAW CODE ***


void uobj_solid::apply_gl_rotate() const {

	rotate_from_v2v(rot_axis, plus_z);
	if (rot_ang != 0.0) rotate_about(rot_ang, plus_z); // in degrees
}


void uobj_solid::rotate_vector(vector3d &v) const {

	rotate_vector3d_by_vr(rot_axis, plus_z, v);
	if (rot_ang != 0.0) rotate_vector3d(plus_z, rot_ang/TO_DEG, v); // in radians
}


void uobj_solid::rotate_vector_inv(vector3d &v) const {

	if (rot_ang != 0.0) rotate_vector3d(plus_z, -rot_ang/TO_DEG, v); // in radians
	rotate_vector3d_by_vr(plus_z, rot_axis, v);
}


float get_pixel_size(float radius, float dist) {
	return min(10000.0f, ((float)window_width)*radius/max(dist, TOLERANCE)); // approx. in pixels
}


void move_in_front_of_far_clip(point_d &pos, point const &camera, float &size, float dist) {

	if (dist > 0.75*FAR_CLIP) { // behind far clipping plane - move closer and scale
		float const pscale(FAR_CLIP/(1.5*dist));
		size *= pscale;
		pos   = camera - (camera - pos)*pscale;
	}
}


bool ustar::draw(point_d pos_, camera_mv_speed const &cmvs, ushader_group &usg, float rscale) {

	vector3d const vcp(cmvs.camera, pos_);
	float const radius0(rscale*radius), dist(vcp.mag() - radius0);
	if (dist > U_VIEW_DIST) return 0; // too far away
	float size(get_pixel_size(radius0, dist)); // approx. in pixels
	// view volume has already been checked
	if (size < 0.02) return 0; // too small
	float const st_prod(STAR_BRIGHTNESS*size*temp);
	if (st_prod < 2.0 || st_prod*temp < 4.0) return 0; // too dim
	move_in_front_of_far_clip(pos_, cmvs.camera, size, (dist + radius0));
	colorRGBA ocolor(color);
		
	if (st_prod < 30.0) {
		float const cmult(0.00111*st_prod*st_prod); // small, attenuate (divide by 900) (what if behind far clipping plane?)
		blend_color(ocolor, ocolor, bkg_color, cmult, 1);
	}
	if (size < 2.5) { // both point cases below, normal is camera->object vector
		bool const small(size < 1.5), draw_as_line(cmvs.speed > 0.25); // lines of light - "warp speed"
		point const normal(cmvs.camera - pos_);
		pos_ = make_pt_global(pos_);

		if (small) {
			if (draw_as_line) {
				universe_pld.add_line(pos_, normal, ocolor, (pos_ - cmvs.motion_vector), normal, ocolor);
			}
			else {
				universe_pld.add_pt(pos_, normal, ocolor);
			}
		}
		else {
			ocolor.do_glColor();
		
			if (draw_as_line) {
				if (!small) glLineWidth(2.0); // 2 pixel width
				draw_line(pos_, (pos_ - cmvs.motion_vector)); // blend and smooth?
				if (!small) glLineWidth(1.0);
			}
			else {
				if (!small) glPointSize(2.0); // 2 pixel diameter
				draw_point(pos_);
				if (!small) glPointSize(1.0);
			}
		}
	}
	else { // sphere
		int ndiv(max(4, min(56, int(NDIV_SIZE_SCALE*sqrt(size)))));
		if (ndiv > 16) ndiv = ndiv & 0xFFFC;
		if (world_mode != WMODE_UNIVERSE) {ndiv = max(4, ndiv/2);} // lower res when in background
		assert(ndiv > 0);
		set_fill_mode();
		glPushMatrix();
		global_translate(pos_);
		usg.enable_star_shader(colorA, colorB);

		if (size > 6) {
			ocolor.do_glColor();
			draw_flare(pos_, all_zeros, 2.5*radius0, 2.5*radius0); // draw star's flare
		}
		apply_gl_rotate();
		glEnable(GL_TEXTURE_2D);
		select_texture(WHITE_TEX);
		WHITE.do_glColor();
		draw_sphere_dlist(all_zeros, radius0, ndiv, 1); // small sphere - use display list
		usg.disable_star_shader();
		glDisable(GL_TEXTURE_2D);
		if (size >= 64) {draw_flares(ndiv, 1);}
		glPopMatrix();
	} // end sphere draw
	return 1;
}


bool urev_body::draw(point_d pos_, camera_mv_speed const &cmvs, ushader_group &usg, float rscale) {

	vector3d const vcp(cmvs.camera, pos_);
	float const radius0(rscale*radius), dist(max(TOLERANCE, (vcp.mag() - radius0)));
	if (dist > U_VIEW_DIST) return 0; // too far away
	float size(get_pixel_size(radius0, dist)); // approx. in pixels
	if (size < 0.25 && !(display_mode & 0x01)) return 0; // too small
	if (!univ_sphere_vis(pos_, radius0))       return 1; // check if in the view volume
		
	if (world_mode == WMODE_UNIVERSE && !(display_mode & 0x01) && dist < FAR_CLIP && get_owner() != NO_OWNER) { // owner color
		set_owner_color(); // lighting is already disabled
		draw_sphere_dlist(make_pt_global(pos_), radius0*max(1.2, 3.0/size), 8, 0); // at least 3 pixels
		return 1;
	}
	move_in_front_of_far_clip(pos_, cmvs.camera, size, (dist + radius0));
	bool enable_l0(0);
	colorRGBA ocolor(color);

	if (world_mode == WMODE_UNIVERSE && !(display_mode & 2)) {
		show_colonizable_liveable(pos_, radius0); // show liveable/colonizable planets/moons
	}
	if (shadowed) {
		// *** use exact shadow? ***
		if (glIsEnabled(GL_LIGHT0)) enable_l0 = 1;
		glDisable(GL_LIGHT0);
	}
	glEnable(GL_LIGHTING);
	if (size < 1.0) {blend_color(ocolor, ocolor, bkg_color, size*size, 1);} // small, attenuate

	if (size < 2.5) { // both point cases below, normal is camera->object vector
		bool const small(size < 1.5);
		pos_ = make_pt_global(pos_);
		ocolor.do_glColor();
		(cmvs.camera - pos_).do_glNormal(); // orient towards camera
		if (!small) glPointSize(2.0); // 2 pixel diameter
		draw_point(pos_);
		if (!small) glPointSize(1.0);
	}
	else { // sphere
		bool const texture(size > MIN_TEX_OBJ_SZ && tid > 0);
		int ndiv(NDIV_SIZE_SCALE*sqrt(size));
		
		if (size < 64.0) {
			ndiv = max(4, min(48, ndiv));
			if (ndiv > 16) {ndiv = ndiv & 0xFFFC;}
		}
		else {
			int const pref_ndiv(min((int)min(MAX_TEXTURE_SIZE, SPHERE_MAX_ND), ndiv));
			for (ndiv = 1; ndiv < pref_ndiv; ndiv <<= 1) {} // make a power of 2
		}
		if (world_mode != WMODE_UNIVERSE) {ndiv = max(4, ndiv/2);} // lower res when in background
		assert(ndiv > 0);
		set_fill_mode();
		glPushMatrix();
		global_translate(pos_);
		if (texture) {usg.enable_planet_shader();}
		apply_gl_rotate();
		
		if (texture) { // texture map
			glEnable(GL_TEXTURE_2D);
			assert(tid > 0);
			bind_2d_texture(tid);
			WHITE.do_glColor();
		}
		else {
			ocolor.do_glColor();
		}
		if (ndiv >= N_SPHERE_DIV) {
			draw_surface(pos_, radius0, size, ndiv);
		}
		else {
			if (surface != NULL) {surface->clear_cache();} // only gets here when the object is visible
			draw_sphere_dlist(all_zeros, radius0, ndiv, texture); // small sphere - use display list
		}
		if (texture) {
			usg.disable_planet_shader();
			glDisable(GL_TEXTURE_2D);
		}
		glPopMatrix();
	} // end sphere draw
	if (enable_l0) glEnable(GL_LIGHT0);
	glDisable(GL_LIGHTING);
	return 1;
}


vector<float> &get_empty_perturb_map(int ndiv) {

	static vector<float> perturb_map;
	perturb_map.resize(ndiv*(ndiv+1));
	return perturb_map;
}


void urev_body::draw_surface(point_d const &pos_, float radius0, float size, int ndiv) {

	RESET_TIME;
	assert(ndiv > 0);
	point viewed_from(get_player_pos() - pos_);
	rotate_vector(viewed_from);
	float *pmap(NULL);
	bool using_hmap(0), gen_data(0);
	bool const SD_TIMETEST(SHOW_SPHERE_TIME && size >= 256.0);

	if (USE_HEIGHTMAP && surface != NULL && surface->has_heightmap()) { // sphere heightmap for planet or moon
		float const hmap_scale(get_hmap_scale());

		if (univ_planet_lod && ndiv == SPHERE_MAX_ND) {
			vector3d dir(get_player_dir()), upv(get_player_up());
			rotate_vector(dir);
			rotate_vector(upv);
			pos_dir_up const pdu(viewed_from, dir, upv, tan_term, sin_term, NEAR_CLIP_SCALED, FAR_CLIP);
			glDisable(GL_TEXTURE_2D);
			
			if (display_mode & 0x20) {
				surface->draw_view_clipped_sphere(pdu, radius0, hmap_scale, this);
			}
			else {
				surface->draw_cube_mapped_sphere(pdu, radius0, hmap_scale, this);
			}
			if (SHOW_SPHERE_TIME) PRINT_TIME("Draw VCS");
			return;
		}
		vector<float> &perturb_map(get_empty_perturb_map(ndiv));
		pmap = &perturb_map.front();

		if (!CACHE_SPHERE_DATA || !surface->sd.equal(all_zeros, radius0, ndiv)) {
			surface->free_dlist();
			gen_data = 1;
			float const cutoff(surface->min_cutoff), omcinv(1.0/(1.0 - cutoff)), rscale(hmap_scale*radius);
			vector<float> const &heightmap(surface->heightmap);
			unsigned const ssize(surface->ssize);
			assert(ssize > 0 && tsize > 0);

			for (int i = 0; i < ndiv; ++i) { // maybe could skip this if perturb_map is the same?
				unsigned const iy((i*ssize)/ndiv), offset((ndiv+1)*i);

				for (int j = 0; j <= ndiv; ++j) {
					//unsigned const ix((j*ssize)/(ndiv+1));
					unsigned const ix((j == ndiv) ? (ssize-1) : (j*ssize)/ndiv);
					float const val(heightmap[iy + ssize*ix]); // x and y are swapped
					perturb_map[j + offset] = rscale*(omcinv*(max(cutoff, val) - cutoff) - 0.5);
				}
			}
		}
		using_hmap = 1;
	} // USE_HEIGHTMAP
	//if (SD_TIMETEST) PRINT_TIME("Sphere Setup");
	bool const use_dlist(USE_SPHERE_DLISTS && surface != NULL);

	if (use_dlist && surface->exec_or_init_dlist()) {
		if (SD_TIMETEST) PRINT_TIME("Dlist Draw");
		return;
	}
	if (CACHE_SPHERE_DATA && using_hmap) {
		assert(surface != NULL);
		if (gen_data) surface->setup_draw_sphere(all_zeros, radius0, -0.5*get_hmap_scale()*radius, ndiv, pmap);
		surface->sd.draw_subdiv_sphere(viewed_from, 1, use_dlist);
		if (SD_TIMETEST) PRINT_TIME("Sphere Draw Fast"); // 39
	}
	else {
		draw_subdiv_sphere(all_zeros, radius0, ndiv, viewed_from, pmap, 1, use_dlist);
		if (SD_TIMETEST) PRINT_TIME("Sphere Draw v2"); // 25
	}
	if (use_dlist) glEndList();
}


void ustar::draw_flares(int ndiv, bool texture) {

	if (solar_flares.empty()) solar_flares.resize(4 + (rand()&3));
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	enable_blend();
	
	for (unsigned i = 0; i < solar_flares.size(); ++i) {
		solar_flares[i].update(color);
		solar_flares[i].draw(radius, max(3, ndiv/4), texture);
	}
	disable_blend();
	glDisable(GL_CULL_FACE);
}


void ustar::solar_flare::gen(colorRGBA const &color) {

	length   = rand_uniform(0.08, 0.16);
	radius   = rand_uniform(0.01, 0.04);
	angle    = rand_uniform(0.00, 360.0);
	time     = 0;
	lifetime = unsigned(rand_uniform(0.5, 2.0)*TICKS_PER_SECOND);
	dir      = signed_rand_vector().get_norm();
	color1   = colorRGBA((color.R + 0.25), color.G, (color.B - 0.1), 0.75*color.A);
	color1.set_valid_color();
	color2   = ALPHA0;
}


void ustar::solar_flare::update(colorRGBA const &color) {

	if (time >= lifetime) time = lifetime = 0;
	if (lifetime == 0) gen(color);
	if (lifetime == 0) return; // gen failed
	time    = min(lifetime, (time + iticks));
	length += 0.02*length*fticks;
	radius += 0.02*radius*fticks;
	angle  += 0.2*fticks;
}


void ustar::solar_flare::draw(float size, int ndiv, bool texture) const {

	if (lifetime == 0) return;
	assert(length > 0.0 && radius > 0.0 && size > 0.0);
	float const t(float(time)/float(lifetime));
	colorRGBA c;
	blend_color(c, color2, color1, t, 1);
	c.do_glColor();
	glPushMatrix();
	uniform_scale(size);
	rotate_about(angle, dir);
	draw_fast_cylinder(point(0.0, 0.0, 0.9), point(0.0, 0.0, (length + 0.9)), radius, 0.0, ndiv, texture);
	glPopMatrix();
}


void urev_body::show_colonizable_liveable(point const &pos_, float radius0) const {

	if (liveable()) {
		GREEN.do_glColor();
	}
	else if (colonizable()) {
		RED.do_glColor();
	}
	else {
		return;
	}
	glPushMatrix();
	global_translate(pos_);
	draw_sphere_dlist(all_zeros, 1.2*radius0, 12, 0);
	glPopMatrix();
}


void uplanet::ensure_rings_texture() {

	if (ring_data.empty() || ring_tid > 0) return; // no rings, or texture already created
	bool const mipmap = 1;
	setup_1d_texture(ring_tid, GL_MODULATE, mipmap, 0, 0, 0);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA8, RING_TEX_SZ, 0, GL_RGBA, GL_UNSIGNED_BYTE, &ring_data.front());
	if (mipmap) {gen_mipmaps(1);}
	check_gl_error(506);
}


void uplanet::draw_prings(ushader_group &usg, upos_point_type const &pos_, float size_, point const &sun_pos, float sun_radius) const {

	if (ring_data.empty()) return;
	assert(ring_tid > 0);
	bind_1d_texture(ring_tid);
	assert(ring_ri > 0.0 && ring_ri < ring_ro);
	usg.enable_ring_shader(make_pt_global(pos_), radius, ring_ri, ring_ro, make_pt_global(sun_pos), sun_radius); // even in !do_texture mode?
	enable_blend(); // must be drawn last
	glPushMatrix();
	global_translate(pos_);
	scale_by(rscale);
	rotate_into_plus_z(rot_axis); // rotate so that rot_axis is in +z
	WHITE.do_glColor();
	draw_tquad(ring_ro, ring_ro, 0.0, 1);
	usg.disable_ring_shader();
	glPopMatrix();
	disable_blend();
}


void uplanet::draw_atmosphere(ushader_group &usg, upos_point_type const &pos_, float size_) const { // *** must be drawn last ***

	if (size_ < 1.5 || atmos == 0) return;
	float const cloud_radius(PLANET_ATM_RSCALE*radius);
	unsigned const ndiv(max(4, min(48, int(4.8*size_))));
	enable_blend();
	glEnable(GL_LIGHTING);
	set_fill_mode();
	glPushMatrix();
	global_translate(pos_);
	apply_gl_rotate();
	glEnable(GL_CULL_FACE);

	if (usg.enable_atmospheric_shader(atmos)) {
		WHITE.do_glColor();
		glCullFace(GL_FRONT);
		draw_subdiv_sphere(all_zeros, cloud_radius, ndiv, 0, 1);
		glCullFace(GL_BACK);
		usg.disable_atmospheric_shader();
	}
	colorRGBA(1.0, 1.0, 1.0, atmos).do_glColor();
	usg.enable_cloud_shader(cloud_tex_repeat);
	draw_subdiv_sphere(all_zeros, cloud_radius, ndiv, 0, 1);
	usg.disable_cloud_shader();
	glDisable(GL_CULL_FACE);
	glPopMatrix();
	disable_blend();
	glDisable(GL_LIGHTING);
}


// *** PROCESSING/QUERY CODE ***


bool umoon::shadowed_by_planet() {

	assert(planet != NULL && planet->system != NULL);
	vector3d const v1(pos, planet->pos), v2(planet->pos, planet->system->sun.pos);
	float const dotp(dot_product(v1, v2));
	if (dotp < 0) return 0;
	float const dps(planet->orbit), rp(planet->radius), rs(planet->system->sun.radius);
	assert(orbit > TOLERANCE && dps > TOLERANCE);
	float const dx(orbit*sin(safe_acosf(dotp/(orbit*dps)))), rx(rp - (orbit/dps)*(rs - rp));
	return (dx < rx);
}


string umoon::get_info() const {

	ostringstream oss;
	oss << "Radius: " << radius << ", Temp: " << temp << ", Resources: " << resources
		<< ", Can Land: " << land_temp_ok() << ", Colonizable: " << colonizable();
	get_owner_info(oss);
	return oss.str();
}

string uplanet::get_info() const {

	ostringstream oss;
	oss << "Radius: " << radius << ", Temp: " << temp << ", Resources: " << resources
		<< ", Water: " << water << ", Atmos: " << atmos 
		<< endl << "Can Land: " << land_temp_ok() << ", Colonizable: " << colonizable() << ", Liveable: " << liveable();
	get_owner_info(oss);
	return oss.str();
}

string ustar::get_info() const {

	ostringstream oss;
	oss << "Radius: " << radius << ", Temp: " << 100.0*temp;
	return oss.str();
}


// if not find_largest then find closest
int universe_t::get_largest_closest_object(s_object &result, point pos, int find_largest, int max_level,
										   int offset, float expand, bool get_destroyed) const
{
	float max_gsize(0.0), min_gdist(CELL_SIZE);
	if (offset) offset_pos(pos);
	point posc(pos);
	UNROLL_3X(posc[i_] += CELL_SIZEo2;)
	result.init();

	// find the correct cell
	point const cell_origin(cells[0][0][0].pos), pos_trans(posc, cell_origin);
	
	for (unsigned d = 0; d < 3; ++d) {
		result.cellxyz[d] = int(pos_trans[d]/CELL_SIZE);
	}
	if (bad_cell_xyz(result.cellxyz)) {
		UNROLL_3X(result.cellxyz[i_] = -1;)
		return 0;
	}
	ucell const &cell(get_cell(result.cellxyz));

	if (max_level == UTYPE_CELL) { // cell
		result.val = 1;
		return 1;
	}
	if (cell.galaxies == NULL) return 0; // not yet allocated
	pos -= cell.pos;
	float const planet_thresh(expand*4.0*MAX_PLANET_EXTENT), moon_thresh(expand*2.0*MAX_PLANET_EXTENT);
	float const pt_sq(planet_thresh*planet_thresh), mt_sq(moon_thresh*moon_thresh);
	static int last_galaxy(-1), last_cluster(-1), last_system(-1);
	unsigned const ng((unsigned)cell.galaxies->size());
	unsigned const go((last_galaxy >= 0 && last_galaxy < int(ng)) ? last_galaxy : 0);
	bool found_system(0);

	for (unsigned gc_ = 0; gc_ < ng && !found_system; ++gc_) { // find galaxy
		unsigned gc(gc_);
		if (gc == 0) gc = go; else if (gc == go) gc = 0;
		ugalaxy const &galaxy((*cell.galaxies)[gc]);
		float const distg(p2p_dist(pos, galaxy.pos));
		if (distg > expand*(galaxy.radius + MAX_SYSTEM_EXTENT)) continue;
		float const galaxy_radius(galaxy.get_radius_at((pos - galaxy.pos).get_norm()));
		if (distg > expand*(galaxy_radius + MAX_SYSTEM_EXTENT)) continue;

		if (max_level == UTYPE_GALAXY || result.object == NULL) { // galaxy
			float const size(galaxy.radius/max(TOLERANCE, distg)); // galaxy_radius?

			if (max_level == UTYPE_GALAXY) {
				if ((find_largest && size > result.size) || (!find_largest && distg < result.dist)) {
					result.assign(gc, -1, -1, size, distg, UTYPE_GALAXY, NULL);
				}
				continue;
			}
			else if ((find_largest && size > max_gsize) || (!find_largest && distg < min_gdist)) {
				result.galaxy = gc;
				result.type   = UTYPE_GALAXY;
				max_gsize     = size;
				min_gdist     = distg;
			}
		}
		unsigned const num_clusters((unsigned)galaxy.clusters.size());
		unsigned const co((last_cluster >= 0 && last_cluster < int(num_clusters) && gc == go) ? last_cluster : 0);

		for (unsigned cl_ = 0; cl_ < num_clusters && !found_system; ++cl_) { // find cluster
			unsigned cl(cl_);
			if (cl == 0) cl = co; else if (cl == co) cl = 0;
			ugalaxy::system_cluster const &cluster(galaxy.clusters[cl]);
			float const testval(expand*cluster.bounds);
			if (p2p_dist_sq(pos, cluster.center) > testval*testval) continue;
			unsigned const cs1(cluster.s1), cs2(cluster.s2);
			unsigned const so((last_system >= int(cs1) && last_system < int(cs2) && cl == co) ? last_system : cs1);

			for (unsigned s_ = cs1; s_ < cs2 && !found_system; ++s_) {
				unsigned s(s_);
				if (s == cs1) s = so; else if (s == so) s = cs1;
				ussystem const &system(galaxy.sols[s]);
				assert(system.cluster_id == cl); // testing
				float const dists_sq(p2p_dist_sq(pos, system.pos)), testval(expand*(system.radius + MAX_PLANET_EXTENT));
				if (dists_sq > testval*testval) continue;
				float dists(sqrt(dists_sq));
				found_system = (expand <= 1.0 && dists < system.radius);
				
				if (system.sun.is_ok() || get_destroyed) {
					dists -= system.sun.radius;
					float const size(system.sun.radius/max(TOLERANCE, dists));

					if (dists <= 0.0) { // sun collision
						result.assign(gc, cl, s, size, dists, UTYPE_STAR, &system.sun);
						result.val = 2;
						return 2; // system
					}
					if ((find_largest && size > result.size) || (!find_largest && dists < result.dist)) {
						result.assign(gc, cl, s, size, dists, UTYPE_SYSTEM, &system.sun);
					}
				}
				if (max_level == UTYPE_SYSTEM || max_level == UTYPE_STAR) continue; // system/star
				unsigned const np((unsigned)system.planets.size());
				
				for (unsigned pc = 0; pc < np; ++pc) { // find planet
					uplanet const &planet(system.planets[pc]);
					float distp_sq(p2p_dist_sq(pos, planet.pos));
					if (distp_sq > pt_sq) continue;
					float distp(sqrt(distp_sq));
					
					if (planet.is_ok() || get_destroyed) {
						distp -= planet.radius;
						float const size(planet.radius/max(TOLERANCE, distp));
					
						if (distp <= 0.0) { // planet collision
							result.assign(gc, cl, s, size, distp, UTYPE_PLANET, &planet);
							result.planet = pc;
							result.val    = 2;
							return 2;
						}
						if ((find_largest && size > result.size) || (!find_largest && distp < result.dist)) {
							result.assign(gc, cl, s, size, distp, UTYPE_PLANET, &planet);
							result.planet = pc;
						}
					}
					if (max_level == UTYPE_PLANET) continue; // planet
					unsigned const nm((unsigned)planet.moons.size());
					
					for (unsigned mc = 0; mc < nm; ++mc) { // find moon
						umoon const &moon(planet.moons[mc]);
						if (!moon.is_ok() && !get_destroyed) continue;
						float const distm_sq(p2p_dist_sq(pos, moon.pos));
						if (distm_sq > mt_sq)                continue;
						float const distm(sqrt(distm_sq) - moon.radius), size(moon.radius/max(TOLERANCE, distm));

						if (distm <= 0.0) { // moon collision
							result.assign(gc, cl, s, size, distm, UTYPE_MOON, &moon);
							result.planet = pc;
							result.moon   = mc;
							result.val    = 1;
							return 2;
						}
						if ((find_largest && size > result.size) || (!find_largest && distm < result.dist)) {
							result.assign(gc, cl, s, size, distm, UTYPE_MOON, &moon);
							result.planet = pc;
							result.moon   = mc;
						}
					} // moon
				} // planet
			} // system
		} // cluster
	} // galaxy
	result.val   = ((result.size > 0.0) ? 1 : -1);
	last_galaxy  = result.galaxy;
	last_cluster = result.cluster;
	last_system  = result.system;
	//cout << last_galaxy << ":" << last_cluster << ":" << last_system << " "; // 0:4:151
	return (result.size > 0.0);
}


bool universe_t::get_trajectory_collisions(s_object &result, point &coll, vector3d dir, point start,
										   float dist, float line_radius) const
{
	int cell_status(1), cs[3];
	float rdist, ldist, t;
	vector3d c1, c2, val, tv, step;
	coll_test ctest;
	static vector<coll_test> gv, sv, pv;

	// check for simple point
	if (dist <= 0.0) {
		int const ival(get_largest_closest_object(result, start, 0, UTYPE_MOON, 1, 1.0));

		if (ival == 2 || (ival == 1 && result.dist <= line_radius)) { // collision at start
			coll = start; return 1;
		}
		coll.assign(0.0, 0.0, 0.0);
		return 0;
	}

	// offset start point
	offset_pos(start);

	// get end point
	dir.normalize();
	point end(start + dir*dist);

	// calculate cell block boundaries
	point const p_low(cells[0][0][0].pos), p_hi(cells[U_BLOCKS-1][U_BLOCKS-1][U_BLOCKS-1].pos);

	for (unsigned d = 0; d < 3; ++d) {
		c1[d] = p_low[d] - CELL_SIZEo2;
		c2[d] = p_hi[d]  + CELL_SIZEo2;
		if (start[d] <= c1[d] || start[d] >= c2[d]) return 0; // start is out of the cell block
	}

	// clip line by cell block boundaries (might be unneccessary due to the above check)
	float const cd[3][2] = {{c1.x, c2.x}, {c1.y, c2.y}, {c1.z, c2.z}};
	if (!do_line_clip(start, end, cd)) return 0; // out of cell block

	// make sure it's inside the cell block
	start += dir*RANGE_OFFSET;
	end   -= dir*RANGE_OFFSET;

	{ // check for start point collision (slow but important)
		bool const fast_test(1);
		int const ival(get_largest_closest_object(result, start, 0, (fast_test ? UTYPE_CELL : UTYPE_MOON), 0, 1.0));
		assert(ival != 0 || result.val != 0);

		if (!fast_test && (ival == 2 || (ival == 1 && result.dist <= line_radius))) { // collision at start
			coll = start;
			return 1;
		}
	}
	coll.assign(0.0, 0.0, 0.0);
	float T(0.0);

	// get points referenced to current cell
	ucell const *cell(&get_cell(result));
	point p0(start);
	start -= cell->pos;
	end   -= cell->pos;
	point curr(start);
	vector3d const dv(end, start);

	for (unsigned d = 0; d < 3; ++d) {
		if (dir[d] > 0.0) {cs[d] = 1; val[d] = c2[d];} else {cs[d] = -1; val[d] = c1[d];} // determine cell step direction
		step[d] = CELL_SIZE*((float)cs[d]); // determine line distance parameters
	}
	gv.resize(0);
	sv.resize(0);
	pv.resize(0);

	// main loop
	while (T <= 1.0) {
		if (cell_status == 0) { // search for next cell
			// find t values
			for (unsigned d = 0; d < 3; ++d) {
				if (dv[d] == 0) {tv[d] = 2.0;} else {tv[d] = (val[d] - start[d])/dv[d];}
			}

			// determine closest edge and update t and cell
			unsigned const dim((tv.x < tv.y) ? ((tv.x < tv.z) ? 0 : 2) : ((tv.y < tv.z) ? 1 : 2)); // smallest dim
			T = tv[dim];
			result.cellxyz[dim] += cs[dim];
			val[dim] += step[dim];
			if (T > 1.0)           return 0; // out of cells
			if (result.bad_cell()) return 0; // off the array
			cell = &get_cell(result);
			curr = p0 - cell->pos;
		}
		if (cell->galaxies == NULL) return 0; // galaxies not yet allocated
		cell_status   = 0;
		result.galaxy = result.cluster = result.system = result.planet = result.moon = -1;
		vector<ugalaxy> const &galaxies(*cell->galaxies);
		assert(galaxies.size() <= MAX_GALAXIES_PER_CELL); // just a consistency check

		for (unsigned i = 0; i < galaxies.size(); ++i) { // search for galaxies
			float const g_radius(galaxies[i].radius + MAX_SYSTEM_EXTENT);
			if (!dist_less_than(curr, galaxies[i].pos, (g_radius + dist))) continue;

			if (line_intersect_sphere(curr, dir, galaxies[i].pos, (g_radius+line_radius), rdist, ldist, t)) {
				ctest.index = i; // line passes through galaxy
				ctest.dist  = ldist;
				gv.push_back(ctest);
			}
		}
		if (gv.empty()) continue; // no intersecting galaxies
		std::sort(gv.begin(), gv.end());

		for (unsigned gc = 0; gc < gv.size(); ++gc) {
			ugalaxy const &galaxy(galaxies[gv[gc].index]);

			for (unsigned c = 0; c < galaxy.clusters.size(); ++c) {
				ugalaxy::system_cluster const &cl(galaxy.clusters[c]);
				if (!dist_less_than(curr, cl.center, (cl.bounds + dist))) continue;
				if (!line_intersect_sphere(curr, dir, cl.center, (cl.bounds+line_radius), rdist, ldist, t)) continue;
				
				for (unsigned i = cl.s1; i < cl.s2; ++i) { // search for systems
					float const s_radius(galaxy.sols[i].radius + MAX_PLANET_EXTENT);
					if (!dist_less_than(curr, galaxy.sols[i].pos, (s_radius + dist))) continue;

					if (line_intersect_sphere(curr, dir, galaxy.sols[i].pos, (s_radius+line_radius), rdist, ldist, t)) {
						ctest.index = i; // line passes through system
						ctest.dist  = ldist;
						ctest.rad   = rdist;
						ctest.t     = t;
						sv.push_back(ctest);
					}
				}
			}
			if (sv.empty()) continue; // no intersecting systems
			std::sort(sv.begin(), sv.end());
			result.galaxy = gv[gc].index;

			for (unsigned sc = 0; sc < sv.size(); ++sc) {
				ussystem const &system(galaxy.sols[sv[sc].index]);
				
				if (system.sun.is_ok() && sv[sc].rad <= system.sun.radius && sv[sc].t > 0.0 && sv[sc].dist <= dist) {
					ctest.index = -1;
					ctest.dist  = sv[sc].dist;
					ctest.rad   = sv[sc].rad;
					ctest.t     = t;
					pv.push_back(ctest); // line intersects sun
				}
				for (unsigned i = 0; i < system.planets.size(); ++i) { // search for planets
					float const p_radius(system.planets[i].mosize);
					if (!dist_less_than(curr, system.planets[i].pos, (p_radius + dist))) continue;

					// *** test against exact planet contour? ***
					if (line_intersect_sphere(curr, dir, system.planets[i].pos, (p_radius+line_radius), rdist, ldist, t)) {
						ctest.index = i;
						ctest.dist  = ldist;
						ctest.rad   = rdist;
						ctest.t     = t;
						pv.push_back(ctest); // line passes through planet orbit
					}
				}
				if (pv.empty()) continue; // no intersecting planets
				std::sort(pv.begin(), pv.end());
				result.system  = sv[sc].index;
				result.cluster = system.cluster_id;

				for (unsigned pc = 0; pc < pv.size(); ++pc) {
					if (pv[pc].index < 0) { // sun is closer
						result.dist   = pv[pc].dist;
						result.object = &(system.sun);
						result.size   = system.sun.radius/result.dist;
						result.type   = UTYPE_STAR;
						coll          = system.sun.pos; // center, not collision point, fix?
						return 1;
					}
					uplanet const &planet(system.planets[pv[pc].index]);

					if (planet.is_ok() && pv[pc].rad <= planet.radius && pv[pc].t > 0.0 && pv[pc].dist <= dist) {
						ctest.index = -1; // line intersects planet
						ctest.dist  = pv[pc].dist;
					}
					else {
						ctest.dist = 2.0*dist;
					}
					for (unsigned i = 0; i < planet.moons.size(); ++i) { // search for moons
						if (planet.moons[i].is_ok()) {
							float const m_radius(planet.moons[i].radius);
							if (!dist_less_than(curr, planet.moons[i].pos, (m_radius + dist))) continue;

							// *** test against exact moon contour? ***
							if (line_intersect_sphere(curr, dir, planet.moons[i].pos, (m_radius+line_radius), rdist, ldist, t)) {
								if (t > 0.0 && ldist <= dist && ldist < ctest.dist) {
									ctest.index = i; // line intersects moon
									ctest.dist  = ldist;
								}
							}
						}
					}
					if (ctest.dist > dist) continue;
					result.planet = pv[pc].index;
					result.dist   = ctest.dist;

					if (ctest.index < 0) { // planet is closer
						result.object = &planet;
						result.size   = planet.radius/result.dist;
						result.type   = UTYPE_PLANET;
					}
					else { // moon is closer
						result.moon   = ctest.index;
						result.object = &(planet.moons[ctest.index]);
						result.size   = planet.moons[ctest.index].radius/result.dist;
						result.type   = UTYPE_MOON;
					}
					coll = result.object->pos; // center, not collision point, fix?
					return 1;
				} // pc
				pv.resize(0);
			} // sc
			sv.resize(0);
		} // gc
		gv.resize(0);
	} // while t
	return 0;
}


float get_temp_in_system(s_object const &clobj, point const &pos) {

	assert(clobj.system >= 0);
	ussystem const &system(clobj.get_system());
	ustar const &sun(system.sun);
	if (!sun.is_ok()) return 0.0; // sun is dead
	point pos2(pos);
	offset_pos(pos2);
	float const sr(sun.radius), rdist(p2p_dist((pos2 - clobj.get_ucell().pos), system.pos));
	return ((rdist < sr) ? 4.0 : 1.0)*sun.get_energy()*min(1.0/(rdist*rdist), 10.0/(sr*sr));
}


float universe_t::get_point_temperature(s_object const &clobj, point const &pos) const {

	if (clobj.system >= 0) return get_temp_in_system(clobj, pos); // existing system is valid
	s_object result; // invalid system - expand the search radius and try again
	if (!get_largest_closest_object(result, pos, 0, UTYPE_SYSTEM, 1, 1.0) || result.system < 0) return 0.0;
	return get_temp_in_system(result, pos);
}


bool get_gravity(s_object &result, point pos, vector3d &gravity, int offset) {

	gravity.assign(0.0, 0.0, 0.0);
	if (result.val != 1) return 0; // only if in the cell block but not collided
	if (result.bad_cell() || result.galaxy < 0 || result.system < 0) return 0;
	if (offset) offset_pos(pos);
	ussystem const &system(result.get_system());
	pos -= result.get_ucell().pos;
	system.sun.add_gravity_vector(gravity, pos); // add sun's gravity
	
	for (unsigned i = 0; i < system.planets.size(); ++i) { // add planets' gravity
		system.planets[i].add_gravity_vector(gravity, pos);
	}
	if (result.planet >= 0) {
		uplanet const &planet(result.get_planet());

		for (unsigned i = 0; i < planet.moons.size(); ++i) { // add moons' gravity
			planet.moons[i].add_gravity_vector(gravity, pos);
		}
	}
	return 1;
}


inline void uobj_solid::set_grav_mass() {

	gravity = radius*density;
	mass    = MASS_SCALE*gravity*radius*radius; // used to be f(r^2) not f(r^3)
}


// *** PARAMS CODE ***


void set_sun_loc_color(point const &pos, colorRGBA const &color, float radius, bool shadowed, bool no_ambient, float a_scale, float d_scale) {

	float uambient[4], udiffuse[4];
	int const light(GL_LIGHT0);
	point const lpos(make_pt_global(pos));

	// set position - cache this?
	set_gl_light_pos(light, lpos, 1.0); // point light source

	// set color
	for (unsigned i = 0; i < 3; ++i) {
		float const ci(WHITE_COMP_D + OM_WCD*color[i]);
		uambient[i]  = (no_ambient ? 0.0 : a_scale*GLOBAL_AMBIENT*ATTEN_AMB_VAL*OM_WCA*color[i]);
		udiffuse[i]  = (shadowed   ? 0.0 : d_scale*ci);
		sun_color[i] = ci;
	}
	uambient[3] = (no_ambient ? 0.0 : 1.0);
	udiffuse[3] = (shadowed   ? 0.0 : 1.0);
	set_colors_and_enable_light(light, uambient, udiffuse);
	//cout << "UA: "; ((colorRGBA *)(&uambient))->print(); cout << endl;
	//cout << "UD: "; ((colorRGBA *)(&udiffuse))->print(); cout << endl;

	// set light attenuation - cache this?
	float const atten2(max(0.25f, min(1.0f, STAR_MAX_SIZE/radius)));
	glLightf(light, GL_CONSTANT_ATTENUATION,  STAR_CONST);
	glLightf(light, GL_LINEAR_ATTENUATION,    atten2*STAR_LINEAR_SCALED);
	glLightf(light, GL_QUADRATIC_ATTENUATION, atten2*STAR_QUAD_SCALED);
}


void set_light_galaxy_ambient_only() {

	glDisable(GL_LIGHT0);
}


void set_ambient_color(colorRGBA const &color) {

	float uambient[4];
	int const a_light(get_universe_ambient_light());
	glEnable(a_light);

	for (unsigned i = 0; i < 3; ++i) {
		uambient[i] = GLOBAL_AMBIENT*BASE_AMBIENT*(WHITE_COMP_A + OM_AAV*OM_WCA*color[i]);
	}
	uambient[3] = 1.0;
	glLightfv(a_light, GL_AMBIENT, uambient);
}


void set_lighting_params() {

	float const ambient[4] = {0.5, 0.5, 0.5, 1.0}, diffuse[4] = {1.0, 1.0, 1.0, 1.0}, zero4[4] = {0.0, 0.0, 0.0, 0.0};
	int const a_light(get_universe_ambient_light()), s_light(GL_LIGHT0);
	set_colors_and_enable_light(s_light, ambient, diffuse); // single star diffuse + ambient
	set_gl_light_pos(s_light, all_zeros, 0.0);
	set_colors_and_enable_light(a_light, ambient, zero4); // universe + galaxy ambient
	set_gl_light_pos(a_light, all_zeros, 0.0);
	set_light_atten(a_light, 1.0);
	glEnable(GL_NORMALIZE);
	//glEnable(GL_BLEND);
	//glEnable(GL_POINT_SMOOTH);
	//glEnable(GL_LINE_SMOOTH);
}


// ****************** CLASS HIERARCHY CODE **********************


void uobject::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass, int align, unsigned eflags, free_obj const *parent_) {

	if (!is_ok()) return;
	status = 1;
	assert(etype < NUM_ETYPES);
	if (etype == ETYPE_NONE) return;
	assert(bradius > 0.0);
	exp_type_params const &ep(et_params[etype]);
	int const etime((exp_time == 0) ? int(ep.duration*max(6, min(20, int(3.0*bradius/radius)))) : exp_time);
	add_blastr(pos, edir, bradius, damage, etime, align, ep.c1, ep.c2, etype, parent_);
	if (damage > 0.0) register_explosion(pos, bradius, damage, eflags, wclass, this, parent_); // delay one frame
}

void uobject::def_explode(float size, int etype, vector3d const &edir, int wclass, int align, unsigned eflags, free_obj const *parent_, float dscale) {
	
	explode(EXPLOSION_DAMAGE*dscale*size*radius, size*radius, etype, edir, 0, wclass, align, eflags, parent_);
	assert(!is_ok());
}


void uobject::add_gravity_vector_base(vector3d &vgravity, point const &mpos, float gfactor, float gmax) const {

	if (!is_ok()) return;
	vector3d const dir(pos, mpos);
	float const d_sq(dir.mag_sq());
	// this is the gravity acceleration - multiply by object mass to get force
	if (d_sq < 1000.0*radius*radius) vgravity += dir*(min(gfactor/d_sq, gmax)*InvSqrt(d_sq));
}


void urev_body::free_texture() { // and also free display list

	if (surface != NULL) surface->free_dlist();
	::free_texture(tid);
	tsize = 0;
}


void uplanet::free_texture() {

	urev_body::free_texture();
	::free_texture(ring_tid);
}


void urev_body::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass,
						int align, unsigned eflags, free_obj const *parent_)
{
	gen_fragments();
	uobject::explode(damage, bradius, etype, edir, exp_time, wclass, align, eflags, parent_);
}


void uobj_rgen::gen_rseeds() { // is this really OS/machine independent (even 32-bit vs. 64-bit)?

	urseed1 = rand2();
	urseed2 = rand2();
}

void uobj_rgen::get_rseeds() {

	urseed1 = global_rand_gen.rseed1;
	urseed2 = global_rand_gen.rseed2;
}

void uobj_rgen::set_rseeds() const {

	global_rand_gen.rseed1 = urseed1;
	global_rand_gen.rseed2 = urseed2;
}


// s_object


bool s_object::write(ostream &out) const {

	if (!out.good()) return 0;
	out << type << " " << cellxyz[0] << " " << cellxyz[1] << " " << cellxyz[2] << " "
		<< galaxy << " " << cluster << " " << system << " " << planet << " " << moon << " " << id;
	return out.good();
}


bool s_object::read(istream &in) {

	if (!in.good()) return 0;
	return ((in >> type >> cellxyz[0] >> cellxyz[1] >> cellxyz[2] >> galaxy >> cluster >> system >> planet >> moon >> id) && in.good());
}


void s_object::init() {

	dist   = CELL_SIZE;
	size   = 0.0;
	galaxy = cluster = system = planet = moon = -1;
	type   = UTYPE_NONE;
	val    = id = 0;
	object = NULL;
	cellxyz[0] = cellxyz[1] = cellxyz[2] = 0;
}


void s_object::assign(int gc, int cl, int sy, float sz, float di, int ty, uobj_solid const *obj) {

	galaxy  = gc;
	cluster = cl;
	system  = sy;
	size    = sz;
	dist    = di;
	type    = ty;
	object  = obj;
}


bool s_object::operator<(const s_object &I) const {

	if (type   < I.type)      return 1; if (type   > I.type)   return 0;
	if (type == UTYPE_NONE)   return 0;
	if (cellxyz[0] < I.cellxyz[0]) return 1; if (cellxyz[0] > I.cellxyz[0]) return 0;
	if (cellxyz[1] < I.cellxyz[1]) return 1; if (cellxyz[1] > I.cellxyz[1]) return 0;
	if (cellxyz[2] < I.cellxyz[2]) return 1; if (cellxyz[2] > I.cellxyz[2]) return 0;
	if (type == UTYPE_CELL)   return 0;
	if (galaxy < I.galaxy)    return 1; if (galaxy > I.galaxy)  return 0;
	if (type == UTYPE_GALAXY) return 0;
	if (cluster < I.cluster)  return 1; if (cluster > I.cluster)return 0; // maybe unnecessary
	if (system < I.system)    return 1; if (system > I.system)  return 0;
	if (type == UTYPE_SYSTEM || type == UTYPE_STAR)             return 0;
	if (planet < I.planet)    return 1; if (planet > I.planet)  return 0;
	if (type == UTYPE_PLANET) return 0;
	if (moon   < I.moon)      return 1; if (moon   > I.moon)    return 0;
	if (type == UTYPE_MOON)   return 0;
	return (id < I.id);
}


void s_object::print() const {

	cout << "type: " << type << ", cell: " << cellxyz[0] << ", " << cellxyz[1] << ", " << cellxyz[2] << ", galaxy: "
		 << galaxy << ", cluster: " << cluster << ", system: " << system << ", planet = " << planet << ", moon = "
		 << moon << ", id = " << id << endl;
}

bool s_object::bad_cell() const {
	return universe.bad_cell_xyz(cellxyz);
}

ucell &s_object::get_ucell() const {
	assert(type >= UTYPE_CELL);
	return universe.get_cell(*this);
}




