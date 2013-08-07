// 3D World - Universe control - code that interfaces with ships and fre_objs
// by Frank Gennari
// 11/11/05

#include "universe.h"
#include "ship.h"
#include "ship_util.h"
#include "asteroid.h"
#include "timetest.h"
#include "openal_wrap.h"


bool const TIMETEST           = (GLOBAL_TIMETEST || 0);
bool const ORBITAL_REGEN      = 0;
bool const PRINT_OWNERSHIP    = 0;
unsigned const GRAV_CHECK_MOD = 4; // must be a multiple of 2


float last_temp(-100.0);
s_object clobj0; // closest object to player
pos_dir_up player_pdu;
cobj_vector_t const empty_cobjs; // always empty
unsigned owner_counts[NUM_ALIGNMENT] = {0};
float resource_counts[NUM_ALIGNMENT] = {0.0};


extern bool univ_planet_lod; // smaller near_clip if true?
extern int uxyz[], window_width, window_height, do_run, fire_key, display_mode, DISABLE_WATER, frame_counter;
extern float zmax, zmin, fticks, univ_temp, temperature, atmosphere, vegetation, base_gravity, urm_static;
extern float def_water_level, water_plane_z, water_h_off_rel, tan_term, sin_term, init_temperature, camera_shake;
extern unsigned team_credits[];
extern string user_text;
extern vector<free_obj *> uobjs;
extern vector<cached_obj> stat_objs;
extern vector<temp_source> temp_sources;
extern universe_t universe;


void process_univ_objects();
void check_shift_universe();
void draw_universe_sun_flare();


void setup_universe_fog(s_object const &closest);
void set_current_system_light(s_object const &clobj, point const &pspos, float a_scale, float d_scale);
void destroy_sobj(s_object const &target);
bool get_gravity(s_object &result, point pos, vector3d &gravity, int offset);
void set_lighting_params();
point get_scaled_upt();



void init_universe_display() {

	setup_ships();
	set_perspective(PERSP_ANGLE, UNIV_NCLIP_SCALE); // that's all (closer near_clip for univ_planet_lod?)
	check_shift_universe();
}


void set_univ_pdu() {

	u_ship const &player(player_ship());
	player_pdu = pos_dir_up(player.get_pos(), player.get_dir(), player.get_up(), tan_term, sin_term, 0.0, U_VIEW_DIST);
}


void check_asserts() {

	assert(CELL_SIZE                  > 2*(GALAXY_MAX_SIZE + MAX_SYSTEM_EXTENT));
	assert(SYSTEM_MIN_SPACING         > 2*MAX_SYSTEM_EXTENT); // galaxies/systems don't overlap
	assert(INTER_PLANET_MIN_SPACING   > PLANET_MAX_SIZE/*2*MAX_PLANET_EXTENT*/); // planets don't overlap
	assert(INTER_MOON_MIN_SPACING     > MOON_MAX_SIZE); // moons don't overlap
	assert(MOON_TO_PLANET_MAX_SPACING > PLANET_MAX_SIZE + MOON_TO_PLANET_MIN_SPACING);
	assert(PLANET_TO_SUN_MAX_SPACING  > STAR_MAX_SIZE   + PLANET_TO_SUN_MIN_SPACING);
	assert(MOON_TO_PLANET_MIN_SPACING > MOON_MAX_SIZE);
	assert(INTER_PLANET_MIN_SPACING   > PLANET_MAX_SIZE);
	assert(PLANET_TO_SUN_MIN_SPACING  > STAR_MAX_SIZE + PLANET_MAX_SIZE);
}


void do_univ_init() {

	static bool univ_inited(0);
	if (univ_inited) return;
	setup_ships(); // just in case
	univ_inited     = 1;
	set_rand2_state(1,1);
	universe.init();
	check_asserts();
	import_default_modmap();
}


void setup_current_system() {
	
	bool regen_mesh(0);
	static int inited(0);
	static float last_water(-1.0);
	do_univ_init();

	if (!inited) {
		inited = 1;
		set_univ_pdu();
		universe.draw_all_cells(clobj0, 1, 1, 2); // required to gen galaxies/systems
	}
	point const &pos(get_player_pos2());
	universe.get_object_closest_to_pos(clobj0, pos, 0, 4.0);
	colorRGBA c1(ALPHA0), c2(ALPHA0);
	float water(0.0);
	atmosphere = 0.0;
	vegetation = 0.0;

	if (clobj0.type == UTYPE_PLANET || clobj0.type == UTYPE_MOON) { // planet or moon
		point const &opos(clobj0.object->get_pos());
		float const oradius(clobj0.object->get_radius());

		if (dist_less_than(pos, opos, 10.0*oradius)) { // fairly close to the planet
			if (!dist_less_than(pos, opos, 1.01*oradius)) { // not at the surface
				//cout << "distance to planet: " << p2p_dist(opos, pos) << ", planet radius: " << oradius << endl;
				point const p_int(opos + (oradius + get_player_radius())*(pos - opos).get_norm());
				player_ship().move_to(p_int); // move the player ship to the surface of the planet/moon
			}
		}
		else {
			// FIXME: what do we do here? still move the player ship? clear the closest planet? refuse to switch world mode? draw empty space?
		}
	}
	if (clobj0.type == UTYPE_PLANET) { // determine atmosphere and cloud cover
		uplanet const &planet(clobj0.get_planet());
		atmosphere = planet.atmos;
		water      = planet.water;
		vegetation = planet.get_vegetation();

		if (!planet.has_vegetation()) {
			c1 = planet.colorA;
			c2 = planet.colorB;
		}
	}
	else if (clobj0.type == UTYPE_MOON) {
		umoon const &moon(clobj0.get_moon());
		atmosphere = max(0.1f, moon.atmos); // 0.1
		water      = moon.water; // 0.0
		c1         = moon.colorA;
		c2         = moon.colorB;
	}
	if (water != last_water) {
		last_water = water;

		if (water == 0.0) {
			if (DISABLE_WATER == 0) DISABLE_WATER = 2; // temporary disable
		}
		else {
			if (DISABLE_WATER == 2) DISABLE_WATER = 0; // enable
			float const water_level(0.5 + 0.5*(water - 0.5)), rel_wpz(get_rel_wpz());

			if (fabs(water_level - rel_wpz) > 0.001) {
				cout << "water: " << water << ", water_level: " << water_level << ", rel_wpz: " << rel_wpz << ", water_h_off_rel: " << water_h_off_rel << endl;
				water_h_off_rel = 0.0; // so that get_rel_wpz() will return the base water level
				water_h_off_rel = water_level - get_rel_wpz();
				regen_mesh      = 1; // regen texture (is this needed?)
				def_water_level = water_plane_z = get_water_z_height(); // ???
				calc_watershed();
			}
		}
	}
	if (clobj0.type >= UTYPE_STAR) { // determine gravity
		assert(clobj0.object != NULL);
		base_gravity = 70.0*clobj0.object->gravity;
	}
	float const a_scale(0.5 + 1.0*atmosphere); // higher atmosphere = more clouds = more ambient lighting (FIXME: cloud color for ambient?)
	set_current_system_light(clobj0, pos, a_scale, 0.5);

	if (fabs(univ_temp - last_temp) > 0.1) { // update temperature
		cout << "temperature: " << univ_temp << endl;
		last_temp   = univ_temp;
		temperature = univ_temp;
		regen_mesh  = 1; // regen texture (is this needed?)
		//init_temperature = univ_temp; // ???
	}
	//if (regen_mesh) *** regen mesh ***
	setup_landscape_tex_colors(c1, c2);

	if (regen_mesh) {
		init_terrain_mesh();
		clear_tiled_terrain();
	}
}


void proc_uobjs_first_frame() {

	for (unsigned i = 0; i < uobjs.size(); ++i) {
		assert(uobjs[i]);
		uobjs[i]->first_frame_hook();
	}
}


void draw_universe(bool static_only, bool skip_closest, int no_distant, bool gen_only) { // should be process_universe()

	RESET_TIME;
	static int inited(0), first_frame_drawn(0);
	set_lighting_params();
	WHITE.do_glColor();
	do_univ_init();

	// clobj0 will not be set - need to draw cells before there are any sobjs
	if (!inited) {static_only = 0;} // force full universe init the first time
	
	if (!static_only) {
		u_ship &ps(player_ship());

		if (fire_key) {
			fire_key = 0;
			ps.try_fire_weapon(); // has to be before process_univ_objects()
		}
		process_univ_objects();
		if (TIMETEST) PRINT_TIME(" Proc Univ Objs");
		apply_explosions();
		if (TIMETEST) PRINT_TIME(" Explosion");
		//camera_origin = get_player_pos2();
		bool burning(ps.is_burning()), damaged(ps.was_damaged());
		ps.clear_damaged();

		if (!ps.is_ok()) {
			damaged = 1; // adds a red filter to the player's view while dying
			destroy_player_ship(0);
		}
		if (damaged) {
			if (camera_shake == 0.0) camera_shake = 1.0;
			add_camera_filter(colorRGBA(1.0, 0.0, 0.0, 0.04), 5, -1, CAM_FILT_DAMAGE);
		}
		if (burning) {
			float const tratio(ps.get_temp()/ps.specs().max_t);
			add_camera_filter(colorRGBA(1.0, 0.0, 0.0, min(0.5, max(0.1, 0.1*tratio))), 3, -1, CAM_FILT_BURN); // NOISE_TEX
		}
	}
	universe.get_object_closest_to_pos(clobj0, get_player_pos2(), 0, 4.0);
	if (!static_only) {setup_universe_fog(clobj0);}
	glEnable(GL_COLOR_MATERIAL);
	check_gl_error(120);
	universe.draw_all_cells(clobj0, skip_closest, skip_closest, no_distant, gen_only);
	check_gl_error(121);
	
	if (!gen_only && !first_frame_drawn) {
		proc_uobjs_first_frame();
		first_frame_drawn = 1;
	}
	if (!gen_only && !static_only) {
		if (TIMETEST) PRINT_TIME(" Universe Draw");
		check_gl_error(122);
		if (!gen_only) {draw_univ_objects(all_zeros);} // draw free objects
		check_gl_error(123);
		if (TIMETEST) PRINT_TIME(" Free Obj Draw");
	}
	check_shift_universe();
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(get_universe_ambient_light());
	glEnable(GL_LIGHT0);
	//draw_universe_sun_flare(); // doesn't look right
	inited = 1;
	if (TIMETEST) PRINT_TIME(" Final Universe");
}


void process_univ_objects() {

	vector<free_obj const*> stat_obj_query_res;

	for (unsigned i = 0; i < uobjs.size(); ++i) { // can we use cached_objs?
		free_obj *const uobj(uobjs[i]);
		bool const no_coll(uobj->no_coll()), particle(uobj->is_particle()), projectile(uobj->is_proj());
		if (no_coll && particle)   continue; // no collisions, gravity, or temperature on this object
		if (uobj->is_stationary()) continue;
		bool const is_ship(uobj->is_ship()), orbiting(uobj->is_orbiting());
		bool const calc_gravity(((uobj->get_time() + (unsigned(uobj)>>8)) & (GRAV_CHECK_MOD-1)) == 0);
		bool const lod_coll(univ_planet_lod && is_ship && uobj->is_player_ship());
		float const radius(uobj->get_c_radius()*(no_coll ? 0.5 : 1.0));
		upos_point_type const &obj_pos(uobj->get_pos());
		vector3d gravity(zero_vector); // sum of gravity from sun, planets, possibly some moons, and possibly asteroids
		point sun_pos(all_zeros);

		// skip orbiting objects (no collisions or gravity effects, temperature is mostly constant)
		s_object clobj; // closest object
		bool const include_asteroids(!particle); // disable particle-asteroid collisions because they're too slow
		int const found_close(orbiting ? 0 : universe.get_object_closest_to_pos(clobj, obj_pos, include_asteroids));
		bool temp_known(0);

		if (found_close) {
			if (clobj.type == UTYPE_ASTEROID) {
				uasteroid const &asteroid(clobj.get_asteroid());
				float const dist_to_cobj(clobj.dist - (asteroid.radius + radius));
				uobj->set_sobj_dist(dist_to_cobj);

				if (dist_to_cobj < 0.0) { // possible collision
					float const elastic((lod_coll ? 0.1 : 1.0)*SBODY_COLL_ELASTIC);
					vector3d const norm(obj_pos, asteroid.pos);
					vector3d const &ascale(asteroid.get_scale());
					double const nmag(norm.mag()), rsum(asteroid.radius*(norm*ascale).mag()/nmag + radius);
					
					if (nmag < rsum) {
						// FIXME: detailed collision?
						// FIXME: projectile explosions damage the asteroid (reduce its radius?)
						point const cpos(asteroid.pos + norm*(rsum/nmag)); // normalize, multiply, and add
						uobj->set_sobj_coll_tid(asteroid.get_fragment_tid(obj_pos));
						uobj->move_to(cpos); // setup correct position for explode?
						uobj->collision(asteroid.pos, asteroid.get_velocity(), S_BODY_DENSITY, asteroid.radius, NULL, elastic); // large mass
						uobj->move_to(cpos); // more accurate since this takes into account the terrain
					}
				}
			}
			else {
				assert(clobj.object != NULL);
				float const clobj_radius(clobj.object->get_radius());
				point const clobj_pos(clobj.object->get_pos());
				float const temperature(universe.get_point_temperature(clobj, obj_pos, sun_pos)*(FOBJ_TEMP_SCALE - uobj->get_shadow_val())); // shadow_val = 0-3
				uobj->set_temp(temperature, sun_pos);
				temp_known = 1;
				float hmap_scale(0.0);
				if (clobj.type == UTYPE_MOON  ) {hmap_scale = MOON_HMAP_SCALE;  }
				if (clobj.type == UTYPE_PLANET) {hmap_scale = PLANET_HMAP_SCALE;}
				float const dist_to_cobj(clobj.dist - (hmap_scale*clobj_radius + radius)); // (1.0 + HMAP_SCALE)*radius?
				uobj->set_sobj_dist(dist_to_cobj);

				if (clobj.type == UTYPE_PLANET || clobj.type == UTYPE_MOON) {
					int coll(0);

					if (dist_to_cobj < 0.0) { // collision (except for stars)
						float coll_r;
						point cpos;
						coll = 1;

						// player_ship and possibly other ships need the more stable but less accurate algorithm
						bool const simple_coll(!is_ship && !projectile);
						float const radius_coll(lod_coll ? 1.25*NEAR_CLIP_SCALED : radius);
						float const elastic((lod_coll ? 0.1 : 1.0)*SBODY_COLL_ELASTIC);

						if (clobj.object->collision(obj_pos, radius_coll, uobj->get_velocity(), cpos, coll_r, simple_coll)) {
							assert(clobj.object->mass > 0.0);
							uobj->set_sobj_coll_tid(clobj.object->get_fragment_tid(obj_pos));
							uobj->move_to(cpos); // setup correct position for explode?
							uobj->collision(clobj_pos, zero_vector, S_BODY_DENSITY*clobj.object->mass, coll_r, NULL, elastic);
							uobj->move_to(cpos); // more accurate since this takes into account the terrain
							coll = 2;
						}
					} // collision
					if (is_ship) uobj->near_sobj(clobj, coll);
				} // planet or moon
				if (calc_gravity) get_gravity(clobj, obj_pos, gravity, 1);
			}
		} // found_close
		if (!temp_known) {
			float temperature(0.0);
			if (!particle && !projectile) {temperature = universe.get_point_temperature(clobj, obj_pos, sun_pos)*FOBJ_TEMP_SCALE;}
			uobj->set_temp(temperature, sun_pos);
		}
		if (calc_gravity) {
			bool near_b_hole(0);
			vector3d swp_accel(zero_vector);

			if (!stat_objs.empty()) {
				all_query_data qdata(&stat_objs, obj_pos, 10.0, urm_static, uobj, stat_obj_query_res);
				get_all_close_objects(qdata);
				
				for (unsigned j = 0; j < stat_obj_query_res.size(); ++j) { // asteroid/black hole gravity
					near_b_hole |= (stat_obj_query_res[j]->get_gravity(gravity, obj_pos) == 2);
				}
			}
			if (clobj.has_valid_system()) {
				swp_accel = clobj.get_star().get_solar_wind_accel(obj_pos, uobj->get_mass(), uobj->get_surf_area());
			}
			uobj->add_gravity_swp(gravity, swp_accel, float(GRAV_CHECK_MOD), near_b_hole);
		}
		if (is_ship) {
			for (unsigned t = 0; t < temp_sources.size(); ++t) { // check for temperature of weapons - inefficient
				temp_source const &ts(temp_sources[t]);
				if (ts.source == uobj) continue; // no self damage
				float const dist_sq(p2p_dist_sq(obj_pos, ts.pos)), rval(ts.radius + radius);
				if (dist_sq > rval*rval)   continue;
				assert(ts.radius > TOLERANCE);
				float const temp(ts.temp*min(1.0f, (rval - sqrt(dist_sq))/ts.radius)*min(1.0, 0.5*max(1.0f, ts.radius/radius)));
				
				if (temp > uobj->get_temp()) {
					uobj->set_temp(temp, ts.pos, ts.source); // source should be valid (and should register as an attacker)
				}
			}
			if (!orbiting) {
				float speed_factor(uobj->get_max_sf()), speed_factor2(1.0);
				if (clobj.val > 0) speed_factor2 = max((lod_coll ? 0.002f : 0.01f), min(1.0f, 0.7f*clobj.dist)); // clip to [0.01, 1.0]
				uobj->set_speed_factor(min(speed_factor, speed_factor2));
			}
		}
	}
}


void reset_player_universe() {

	change_speed_mode(do_run);
	update_cpos();
	reset_player_ship();
	if (MOVE_PLAYER_RPOS) return; // ???
	uxyz[0] = uxyz[1] = uxyz[2] = 0;
	shift_univ_objs(get_scaled_upt(), 0);
	universe.init();
}


void check_shift_universe() {

	point camera(get_player_pos2());
	vector3d move(zero_vector);
	bool moved(0);

	for (unsigned d = 0; d < 3; ++d) { // max move distance is CELL_SIZE
		int sh[3] = {0, 0, 0};

		for (int sign = -1; sign <= 1; sign += 2) {
			while ((camera[d] < sign*CELL_SIZEo2) ^ (sign > 0)) {
				camera[d] -= sign*CELL_SIZE;
				move[d]   -= sign*CELL_SIZE;
				uxyz[d]   += sign;
				sh[d]      = sign;
				moved      = 1;
				universe.shift_cells(sh[0], sh[1], sh[2]);
			}
		}
	}
	if (moved) {shift_univ_objs(move, 1);} // advance all free objects by a cell
}


void fire_planet_killer(u_ship const *const ship, point const &ship_pos, vector3d const &fire_dir, float fire_range, int obj_types) {

	point coll; // unused
	s_object target;

	// test for collisions with sloid stellar objects
	bool const sobj_coll((obj_types & OBJ_TYPE_UOBJ) ? universe.get_trajectory_collisions(target, coll, fire_dir, ship_pos, fire_range, 0.0) : 0);
	bool sobjc(sobj_coll && target.is_solid());

	// test for other ship collisions
	int const itype(obj_types & ~OBJ_TYPE_UOBJ);
	free_obj *fobj = NULL;
	uobject const *obj = NULL;
	float f_dist(0.0);

	if (itype != 0) {
		line_int_data li_data(ship_pos, fire_dir, fire_range, ship, NULL, 0, 0);
		li_data.even_ncoll = 1;
		obj    = line_intersect_objects(li_data, fobj, itype);
		f_dist = li_data.dist;
	}
	bool fobjc(obj != NULL);

	if (sobjc && fobjc) { // two collisions - determine which collision was first
		if (f_dist < target.dist) sobjc = 0; else fobjc = 0;
	}
	if (fobjc) { // hit a ship or projectile before the stellar object
		assert(fobj != NULL);
		vector3d const hit_pos(ship_pos, fobj->get_pos());
		fobj->damage(1.0E6, DAMAGE_DESTROY, hit_pos.get_norm(), ship, UWEAP_DESTROY);
	}
	if (sobjc) { // destroy a stellar object
		if (!ship->invalid() && (target.type == UTYPE_PLANET || target.type == UTYPE_MOON)) {
			int const owner(target.get_world().owner);
			if (owner != NO_OWNER) register_attack_from(ship, owner); // team takes offense
		}
		destroy_sobj(target);
	}
}


// test for intersections with sloid stellar objects (unused)
bool universe_intersection_test(point const &pos, vector3d const &dir, float range) {

	point coll; // unused
	s_object target;
	return (universe.get_trajectory_collisions(target, coll, dir, pos, range, 0.0) && target.is_solid());
}


bool universe_ray_intersect(point const &start, point const &end, int obj_types, uobject const *cur=NULL, free_obj const *ignore=NULL) {

	vector3d const dir((end - start).get_norm());
	float const range(p2p_dist(start, end));
	//if ((obj_types & OBJ_TYPE_UOBJ) && universe_intersection_test(start, dir, range)) return 1;
	line_int_data li_data(start, dir, range, cur, ignore, 0, 0); // not sure how cur and ignore are used, or if it makes sense to pass them in here
	li_data.even_ncoll = 1;
	free_obj *fobj = NULL; // unused
	uobject const *const uobj(line_intersect_objects(li_data, fobj, obj_types));
	//if (uobj != NULL && uobj != cur && uobj != ignore) {cout << "intersects " << uobj->get_name() << endl;}
	return (uobj != NULL && uobj != cur && uobj != ignore);
}


void draw_universe_sun_flare() {

	if (clobj0.type < UTYPE_SYSTEM) return; // no star
	ustar const &sun(clobj0.get_star());
	if (!sun.is_ok() || !univ_sphere_vis((sun.pos + clobj0.get_ucell().rel_center), sun.radius)) return;
	float intensity(1.0);
	u_ship const &ps(player_ship());
	point const viewer(ps.get_pos());
	unsigned const npts = 16;
	static point pts[npts];
	static bool pts_valid(0);
	unsigned nvis(0);
	
	for (unsigned i = 0; i < npts; ++i) {
		if (!pts_valid) {pts[i] = signed_rand_vector_norm();}
		point const pos(sun.pos + pts[i]*sun.radius);
		if (!universe_ray_intersect(viewer, pos, OBJ_TYPE_ALL, &sun, &ps)) {++nvis;}
	}
	pts_valid = 1;
	if (nvis == 0) return;
	intensity = 0.1 + 0.9*float(nvis)/float(npts);
	point const gv(make_pt_global(viewer));
	DoFlares(gv, (gv + player_ship().get_dir()), make_pt_global(sun.pos), 0.01, 0.02, intensity, 4); // draw only some flares
}


void send_warning_message(string const &msg) {

	print_text_onscreen(msg.c_str(), RED, 1.0, 3*TICKS_PER_SECOND/2, 1);
	static int last_warning_frame(0);
	if ((frame_counter - last_warning_frame) > 5.0*TICKS_PER_SECOND) gen_sound(SOUND_ALERT, get_player_pos2());
	last_warning_frame = frame_counter;
}


void disable_player_ship() {

	if (player_ship().is_resetting()) return; // already destroyed
	print_text_onscreen("Ship Disabled", RED, 1.2, 2*TICKS_PER_SECOND, 2);
	add_camera_filter(colorRGBA(1.0, 0.5, 0.0, 0.1), 8, -1, CAM_FILT_DAMAGE);
}


void destroy_player_ship(bool captured) {

	if (player_ship().is_resetting()) return; // already destroyed
	do_run = 0;
	print_text_onscreen((captured ? "Ship Captured" : "Ship Destroyed"), RED, 1.2, 2*TICKS_PER_SECOND, 2);
	add_camera_filter(colorRGBA(1.0, 0.0, 0.0, 0.6), 8, -1, CAM_FILT_DAMAGE);
	if (!captured) gen_sound(SOUND_EXPLODE, get_player_pos2());
	player_ship().reset_after(TICKS_PER_SECOND);
}


void draw_universe_stats() {

	char text[128];
	float const aspect_ratio((float)window_width/(float)window_height);
	u_ship const &ps(player_ship());

	// draw position and orientation
	point const camera(ps.get_pos()), camera_scaled(camera/CELL_SIZE);
	vector3d const dir(ps.get_dir().get_norm());
	float const player_temp(player_ship().get_true_temp());
	YELLOW.do_glColor();
	//sprintf(text, "Loc: (%i: %3.3f, %i: %3.3f, %i: %3.3f)  Dir: (%1.3f, %1.3f, %1.3f)", uxyz[0], camera_scaled.x, uxyz[1], camera_scaled.y, uxyz[2], camera_scaled.z, dir.x, dir.y, dir.z);
	sprintf(text, "Loc: (%3.4f, %3.4f, %3.4f)  Dir: (%1.3f, %1.3f, %1.3f)  T: %3.1f",
		uxyz[0]+camera_scaled.x, uxyz[1]+camera_scaled.y, uxyz[2]+camera_scaled.z, dir.x, dir.y, dir.z, player_temp);
	draw_text(-0.009*aspect_ratio, -0.014, -0.028, text);

	// draw shields, armor, weapon status, etc.
	int const shields(int(100.0*ps.get_shields()/ps.get_max_shields()));
	int const armor(int(100.0*ps.get_armor()/ps.get_max_armor()));
	RED.do_glColor();

	if (ps.need_ammo()) {
		int const ammo(ps.get_ammo());
		sprintf(text, "Shields %d  Armor %d  Ammo %d %s", shields, armor, ammo, ((ps.get_wnum() > 0) ? "" : "[No Weapon]"));
	}
	else {
		sprintf(text, "Shields %d  Armor %d", shields, armor);
	}
	draw_text(-0.006*aspect_ratio, -0.011, -0.02, text);

	// draw credits, kills, total kills
	WHITE.do_glColor();
	sprintf(text, "Credits: %u", (ps.ncredits + team_credits[ALIGN_PLAYER]));
	draw_text(0.0086*aspect_ratio, 0.0135, -0.025, text);
	sprintf(text, "Kills: %u/%u", ps.get_num_kills(), ps.get_tot_kills());
	draw_text(0.0086*aspect_ratio, 0.0126, -0.025, text);
	sprintf(text, "Deaths: %u", ps.get_num_deaths());
	draw_text(0.0086*aspect_ratio, 0.0117, -0.025, text);

	// draw message text (user typed)
	WHITE.do_glColor();
	if (!user_text.empty()) draw_text(-0.010*aspect_ratio, 0.010, -0.02, user_text.c_str()); // x and z are scaled
}


void exec_universe_text(string const &text) {

	if (text.empty())   return;
	if (text[0] != '.') return; // not a command
	
	if (text == ".self-destruct") {
		destroy_player_ship(0);
		return;
	}
	// *** WRITE ***
}


bool rename_obj(uobject *obj, unsigned alignment) { // a little difficult to use, have to enter the text first then click on an object

	assert(obj);
	int const owner(obj->get_owner());
	if (/*owner != NO_OWNER &&*/ owner != alignment) return 0;
	if (user_text.empty()) return 0;

	if (obj->rename(user_text)) {
		string const str(string("Renaming ") + obj->get_name() + " to " + user_text);
		print_text_onscreen(str, PURPLE, 0.8, TICKS_PER_SECOND);
		return 1;
	}
	return 0;
}


bool sphere_intersect_uobject(point const &pos, float radius, bool include_asteroids) {

	s_object result;
	if (!universe.get_closest_object(result, pos, UTYPE_MOON, include_asteroids, 1, 1.0)) return 0;
	if (!result.object || !result.object->is_ok()) return 0; // can this happen?
	if (!result.object->sphere_intersection(pos, radius)) return 0;
	//cout << "intersects with " << result.object->get_name() << endl;
	return 1;
}


bool get_closest_object(point const &pos, s_object &result, int obj_type, bool include_asteroids, bool get_destroyed=0) {

	if (!universe.get_closest_object(result, pos, obj_type, include_asteroids, 1, 4.0, get_destroyed)) return 0;
	return (result.type == obj_type && result.object != NULL && (get_destroyed || result.object->is_ok()));
}


uobject const *get_closest_world_ptr(point const &pos, int type) {

	s_object result;
	return (get_closest_object(pos, result, type, 0) ? result.object : NULL);
}


float urev_body::get_land_value(unsigned align, point const &cur_pos, float sradius) const {

	float owner_val(0.5);
	
	if (owner == NO_OWNER) {
		unsigned const res_c(sclasses[USC_COLONY].cost + max(sclasses[USC_DEFSAT].cost, sclasses[USC_ANTI_MISS].cost));
		owner_val = ((team_credits[align] >= res_c) ? 2.5 : 0.75);
	}
	else if (align == owner || align == ALIGN_PIRATE) {
		owner_val = 1.0;
	}
	else if (TEAM_ALIGNED(align) && TEAM_ALIGNED(owner)) {
		owner_val = 2.0;
	}
	float value(2.0*rand_float());
	value += 1.0*liveable();
	value += (0.1*resources + 0.5)*owner_val;
	value += sradius/PLANET_TO_SUN_MAX_SPACING; // large system bonus
	value -= 0.5*p2p_dist(cur_pos, pos);
	return value;
}


bool line_intersect_sun(point const &p1, point const &p2, s_object const &result, float halo) {

	// avoid choosing a destination that requires flying through a star
	if (!result.has_valid_system()) return 0;
	ustar const &sun(result.get_system().sun);
	return (sun.is_ok() && line_sphere_intersect(p1, p2, sun.pos, halo*sun.radius));
}


uobject const *choose_dest_world(point const &pos, int exclude_id, unsigned align) {

	s_object result;
	float const g_expand(CELL_SIZE/GALAXY_MIN_SIZE); // Note: can be in more than one galaxy, but should be OK
	if (!universe.get_closest_object(result, pos, UTYPE_SYSTEM, 0, 1, 4.0, 0, g_expand) || result.type < UTYPE_GALAXY) return NULL;
	ugalaxy const &galaxy(result.get_galaxy());
	uobject const *dest = NULL;
	float const distval(2.5*galaxy.get_radius());
	float max_svalue(0.0);

	for (unsigned s = 0; s < galaxy.sols.size(); ++s) { // clusters don't help much here
		ussystem const &system(galaxy.sols[s]);
		if (system.planets.empty()) continue; // no planets to visit (could be that the planets aren't generated yet)
		float const svalue(distval*rand_float() - p2p_dist(pos, system.get_pos()));
		
		if (max_svalue == 0.0 || svalue > max_svalue) {
			float const sradius(system.get_radius());
			float max_pvalue(0.0);

			for (unsigned j = 0; j < system.planets.size(); ++j) {
				uplanet const &planet(system.planets[j]);
				bool const planet_acceptable(planet.colonizable() && planet.get_id() != exclude_id);
				if (!planet_acceptable && planet.moons.empty()) continue;
				float const pvalue(planet.get_land_value(align, pos, sradius));

				if (max_pvalue == 0.0 || pvalue > max_pvalue) {
					if (planet_acceptable && !line_intersect_sun(pos, planet.pos, result, 2.0)) {
						max_pvalue = pvalue;
						max_svalue = svalue;
						dest       = &planet;
					}
					for (unsigned k = 0; k < planet.moons.size(); ++k) {
						umoon const &moon(planet.moons[k]);
						if (!moon.colonizable() || moon.get_id() == exclude_id) continue;
						float const mvalue(moon.get_land_value(align, pos, sradius));

						if ((max_pvalue == 0.0 || mvalue > max_pvalue) && !line_intersect_sun(pos, moon.pos, result, 2.0)) {
							max_pvalue = mvalue;
							max_svalue = svalue;
							dest       = &moon;
						}
					}
				}
			}
		}
	}
	if (!dest) cout << "no dest" << endl; // testing
	return dest;
}


inline string get_owner_name(unsigned owner) {

	assert(owner < NUM_ALIGNMENT);
	return align_names[owner];
}


bool check_dest_ownership(int uobj_id, point const &pos, free_obj *own, bool check_for_land, bool homeworld) {
	
	assert(uobj_id >= 0 && own); // obj may be an invalid pointer though
	unsigned const owner(own->get_align());
	s_object result;
	int const world_types[2] = {UTYPE_PLANET, UTYPE_MOON};

	for (unsigned t = 0; t < 2; ++t) {
		if (!get_closest_object(pos, result, world_types[t], 0)) continue;
		if (result.object->get_id() != uobj_id)                  continue;
		urev_body &world(result.get_world());
		if (world.owner != NO_OWNER || !world.colonizable())     continue;
		if (check_for_land && !world.can_land_at(pos))           continue; // currently player only
		if (PRINT_OWNERSHIP) cout << world.get_name() << " is claimed by " << get_owner_name(owner) << "." << endl;
		float defend(world.resources);
		bool owned(0);

		while (defend > 8.0) {
			owned  |= (add_orbiting_ship(USC_DEFSAT, 0, 0, 0, own, &world) != NULL); // put a defense satellite in orbit
			defend -= 12.0;
		}
		if (defend > 0.0 || !owned) {
			owned  |= (add_orbiting_ship(USC_ANTI_MISS, 0, 0, 0, own, &world) != NULL); // put an anti-missile drone in orbit
		}
		float const rsc_val(world.resources - (homeworld ? 0.0 : 20.0));

		if (rsc_val > 0.0 && t == 0) { // planets only - colonies (first homeworld only?)
			unsigned const colony_types[5] = {USC_COLONY, USC_ARMED_COL, USC_HW_COL, USC_STARPORT, USC_HW_SPORT};
			unsigned start_val(0);
			for (start_val = 0; start_val < 4 && world.resources > 10.0*(start_val+1); ++start_val) {}
			if (start_val == 3 || (start_val == 2 && (rand()&1))) ++start_val;
			orbiting_ship const *oship(NULL);

			for (int i = start_val; i >= 0; --i) {
				oship = add_orbiting_ship(colony_types[i], 0, 1, 1, own, &world); // put a colony on world
				
				if (oship != NULL) {
					vector3d const dir((own->pos - oship->pos).get_norm());
					own->move_to(oship->get_pos() + dir*(1.1*(own->get_c_radius() + oship->get_c_radius()))); // move away from object
					break;
				}
			}
			owned |= (oship != NULL);
		}
		if (owned) {
			world.set_owner(result, owner); // only if inhabitable?
		}
		else {
			assert(world.get_owner() == NO_OWNER);
		}
		return 1;
	}
	return 0;
}


// ************ UOBJ_SOLID/UREV_BODY ************


bool uobj_solid::collision(point const &p, float rad, vector3d const &v, point &cpos, float &coll_r, bool simple) const { // maybe should be in universe.cpp

	coll_r = radius;
	if (!surface_test(rad, p, coll_r, simple)) return 0;
	vector3d const norm(p, pos);
	double const rsum(coll_r + rad), nmag(norm.mag());
	if (nmag > rsum) return 0;
	double const vmag(v.mag());

	if ((simple && nmag > TOLERANCE) || vmag < TOLERANCE || dot_product(v, norm) > -0.1*(nmag*vmag)) { // use for shallow coll angles
		cpos = pos + norm*(rsum/nmag); // normalize, multiply, and add
	}
	else {
		get_sphere_mov_sphere_int_pt(pos, p, v, rsum, cpos);
		if (nmag > TOLERANCE) cpos += norm*(0.05*rad/nmag); // slight adjustment to improve stability
	}
	return 1;
}


void urev_body::get_owner_info(ostringstream &oss) const {

	if (owner == NO_OWNER) {
		oss << endl << "Uninhabited";
	}
	else {
		assert(unsigned(owner) < NUM_ALIGNMENT);
		oss << endl << "Owned by " << get_owner_name(owner);
	}
}


void urev_body::set_owner(s_object const &sobj, unsigned owner_) {

	set_owner_int(owner_);
	sobj.set_owner(owner_);
}


void urev_body::set_owner_int(unsigned owner_) {

	if (owner == owner_) return; // already set
	unset_owner();
	
	if (owner_ != NO_OWNER) {
		assert(resources > 0.0);
		assert(owner_ < NUM_ALIGNMENT);
		++owner_counts[owner_];
		resource_counts[owner_] += resources;
	}
	owner = owner_; // may overwrite an old value
}


void urev_body::unset_owner() {

	if (owner != NO_OWNER) {
		assert(unsigned(owner) < NUM_ALIGNMENT);
		assert(owner_counts[owner] > 0); // testing
		--owner_counts[owner]; // shouldn't go negative even if read from a modmap
		assert(resources > 0.0);
		assert(resource_counts[owner] >= resources-1.0); // testing (account for fp error)
		resource_counts[owner] -= resources;
		resources = max(0.0f, resources); // in case it's slightly negative due to fp error
		owner = NO_OWNER;
	}
}


void urev_body::check_owner(s_object const &sobj) {

	// Note: nothing done here
	// if the orbiting object(s) have been destroyed since the planet was deleted, it is no longer owned
	// if the orbiting object(s) still exist, they will set ownership back in their update functions
	//int const new_owner(sobj.get_owner());
	//if (new_owner != NO_OWNER && owner == NO_OWNER) {++orbiting_refs;} // ???
	//set_owner_int(new_owner);
}


void urev_body::set_owner_color() const {

	if (owner == NO_OWNER) return;
	assert(unsigned(owner) < NUM_ALIGNMENT);
	alignment_colors[owner].do_glColor();
}


// ************ U_SHIP ************


bool have_resources_to_colonize(unsigned alignment) {

	unsigned s_types[2] = {USC_DEFSAT, USC_ANTI_MISS};

	for (unsigned i = 0; i < sizeof(s_types)/sizeof(unsigned); ++i) {
		if (get_ship_cost(s_types[i], alignment) >= 0.0) return 1;
	}
	return 0;
}


void u_ship::near_sobj(s_object &clobj, int coll) {

	if (invalid_priv() || !powered_priv()) return;
	assert(clobj.object != NULL);

	if (is_player_ship()) {
		if (!univ_planet_lod && coll == 2) {
			bool const homeworld(clobj.type == UTYPE_PLANET && !has_homeworld());

			if (check_dest_ownership(clobj.object->get_id(), pos, this, 1, homeworld)) {
				if (homeworld) claim_world(clobj.object);
				string const str(string("You have claimed ") + clobj.object->get_name() + (homeworld ? " as your homeworld" : ""));
				print_text_onscreen(str, CYAN, 1.2, 2*TICKS_PER_SECOND);
				gen_sound((homeworld ? SOUND_POWERUP : SOUND_ITEM), get_player_pos2());
			}
		}
	}
	else if (clobj.object->get_owner() == NO_OWNER && clobj.get_world().colonizable() &&
		!is_fighter() && !is_orbiting() && !is_rand_spawn() && can_move() && have_resources_to_colonize(alignment))
	{
		float const odist_sq(p2p_dist_sq(pos, clobj.object->get_pos()));

		if (coll || !dest_mgr.is_valid() || dest_mgr.is_cur_obj(clobj.object->get_id()) || odist_sq < 0.05*p2p_dist_sq(pos, dest_mgr.get_pos())) {
			if (coll || target_obj == NULL || odist_sq < 0.5*p2p_dist_sq(pos, target_obj->get_pos())) {
				dest_override = 1;
				dest_mgr.set_object(clobj.object); // the original destination is lost (if there was one)
			}
		}
	}
}


// ************ ORBITING_SHIP ************


orbiting_ship *add_orbiting_ship(unsigned sclass, bool guardian, bool on_surface, bool pos_from_parent, free_obj const *parent, urev_body *obj) {

	assert(obj && parent);
	//assert(obj->get_owner() == parent->get_align());
	if (!alloc_resources_for(sclass, parent->get_align(), 0)) return NULL;
	obj->inc_orbiting_refs();
	float angle(rand_uniform(0.0, TWO_PI));
	point const parent_pos(parent->get_pos());

	if (pos_from_parent) { // set initial location based on position of ship relative to planet
		vector3d dir((parent->get_pos() - obj->get_pos()).get_norm());
		obj->rotate_vector(dir);
		angle = atan2(dir.y, dir.x);
	}

	// not sure what to set orbit_radius to - could collide with other orbiting objects
	float const orbit_radius(on_surface ? 0.0 : 2.0*obj->get_radius()), rate(0.0);
	point const start_pos((pos_from_parent && on_surface) ? point(parent->get_pos() - obj->get_pos()) : all_zeros);
	orbiting_ship *const ship(new orbiting_ship(sclass, parent->get_align(), guardian, obj, zero_vector, start_pos,
		orbit_radius, angle, rate));
	ship->set_parent(parent);
	add_uobj_ship(ship);
	return ship;
}


// geostationary orbit (GEO) and geosynchronous orbit (GSO)
orbiting_ship::orbiting_ship(unsigned sclass_, unsigned align, bool guardian, urev_body const *obj,
							 vector3d const &axis_, point const &start_pos, float rad, float start_ang, float rate)
	: u_ship(sclass_, all_zeros, align, (AI_ATT_ENEMY | (guardian ? AI_GUARDIAN : 0)), TARGET_CLOSEST, 0), GSO(rate == 0.0),
	fixed_pos(start_pos != all_zeros), has_sobj(0), sobj_liveable(obj->liveable()), orbiting_type(obj->type),
	last_build_time(time), orbit_r(rad), start_angle(start_ang), rot_rate(rate), axis(axis_), rel_pos(start_pos)
{
	assert(obj);
	assert(orbiting_type == UTYPE_PLANET || orbiting_type == UTYPE_MOON);
	assert(!can_move());
	homeworld.set_object(obj);
	flags |= OBJ_FLAGS_ORBT;

	if (orbit_r == 0.0) {
		assert(specs().mass > 0.0);
		// *** determine radius of orbit orbit_r based on mass? ***
	}
	else {
		obj->get_valid_orbit_r(orbit_r, c_radius);
	}
	assert(radius < obj->radius); // too strict?
	assert(orbit_r == 0.0 || orbit_r >= obj->radius + c_radius); // too strict?
	if (orbit_r == 0.0) orbit_r = obj->radius + 0.5*radius; // + radius?
	if (GSO || axis == zero_vector) axis = obj->rot_axis; // GEO (orbits along the equator)

	if (fixed_pos) {
		obj->rotate_vector(rel_pos); // convert to local object space
		rel_pos.normalize();
	}
	set_pos_from_sobj(obj);
}


void orbiting_ship::update_state() {

	if (!is_ok()) return;
	bool const exploding_now(explode_now());
	has_sobj = 0;
	
	if (!(flags & OBJ_FLAGS_DIST) || exploding_now || (time&31) == 0) { // not too far away
		s_object result;

		if (get_closest_object(pos, result, orbiting_type, 0, 1)) {
			urev_body &world(result.get_world());

			if (!world.is_ok()) { // was destroyed
				if (!is_exploding()) destroy_ship(0.0);
				return;
			}
			if (homeworld.update_pos_if_close(&world)) { // the object we are orbiting is still there
				set_pos_from_sobj(&world);
				has_sobj = 1;
			}
			if (!ORBITAL_REGEN && exploding_now && world.get_owner() == alignment) { // is this correct?
				world.dec_orbiting_refs(result); // will die this frame
			}
			else if (!is_exploding() && world.get_owner() == NO_OWNER) {
				world.set_owner(result, alignment); // have to reset - world must have been regenerated
			}
			if (world.temp > get_temp()) set_temp(FOBJ_TEMP_SCALE*world.temp, world.get_pos(), NULL);
		}
	}
	if (!has_sobj) velocity = zero_vector; // stopped
}


void orbiting_ship::set_pos_from_sobj(urev_body const *const sobj) {

	vector3d delta(1.0, 0.0, 0.0); // should delta start out in x or y?

	if (fixed_pos) { // doesn't have to be at the equator (colony)
		delta = rel_pos;
		sobj->rotate_vector_inv(delta);
		delta.normalize();
	}
	else {
		angle = (GSO ? (start_angle + sobj->rot_ang/TO_DEG) : (angle + fticks*rot_rate));
		rotate_vector3d_by_vr(plus_z, axis, delta); // switch to local coordinate system
		rotate_vector3d_norm (axis, -angle, delta); // negate angle?
	}
	pos       = sobj->pos + delta*double(orbit_r);
	reset_pos = pos;
	if (specs().max_turn  == 0.0) dir = delta; // face outward (no turning?)
	if (specs().roll_rate == 0.0) orthogonalize_dir(axis, dir, upv, 1);
	// ignore the rotation and revolution of the object being orbited for now
	velocity = (GSO ? zero_vector : cross_product(delta, axis)*TWO_PI*rot_rate); // similar to calc_angular_vel()
	invalidate_rotv();
}


void orbiting_ship::apply_physics() {
	
	update_state();
	u_ship::apply_physics();
}


bool orbiting_ship::regen_enabled() const { // if has_sobj, shouldn't have to check homeworld
	
	return (ORBITAL_REGEN && has_sobj && has_homeworld());
}


cobj_vector_t const &uobject::get_cobjs() const {return empty_cobjs;}
bool uobject::sphere_intersection(point const &c, float r) const {return dist_less_than(c, pos, (r + radius));}


void uobject::gen_fragments(upos_point_type const &pos_offset, float rscale) const {

	unsigned const num_fragments((rand()&7) + 8);
	int const tex_id(get_fragment_tid(all_zeros));

	for (unsigned i = 0; i < num_fragments; ++i) {
		add_uobj(uobj_asteroid::create((pos_offset + pos + signed_rand_vector(1.2*rscale*radius)),
			0.2*rscale*radius*rand_uniform(0.5, 1.0), AS_MODEL_SPHERE, tex_id, (10*TICKS_PER_SECOND + rand()%TICKS_PER_SECOND))); // temporary
	}
	gen_moving_fragments((pos_offset + pos), (40 + (rand()%20)), tex_id, rscale);
}


void uobject::gen_moving_fragments(point const &hit_pos, unsigned num, int tid, float rscale, float vscale, vector3d const &vadd, colorRGBA const &pcolor) const {

	for (unsigned i = 0; i < num; ++i) {
		unsigned const ltime(5*TICKS_PER_SECOND + rand()%TICKS_PER_SECOND);
		float const fragment_radius(rscale*radius*rand_uniform(0.05, 0.1));

		for (unsigned n = 0; n < 4; ++n) { // make 4 attempts at generating the particle in a valid location
			point ppos(hit_pos + signed_rand_vector(0.5*rscale*radius));
			ppos += (ppos - pos).get_norm()*0.1*radius;

			if (!check_fragment_self_coll() || !sphere_intersection(ppos, fragment_radius)) {
				vector3d const vel(vadd + ((ppos - pos).get_norm() + signed_rand_vector(0.25))*radius*vscale*0.02);
				gen_particle(PTYPE_SPHERE, pcolor, pcolor, ltime, ppos, vel, fragment_radius, 0.0, ALIGN_NEUTRAL, 1, tid);
				break;
			}
		}
	}
}


void ustar::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass,
					int align, unsigned eflags, free_obj const *parent_)
{
	unsigned const num_parts(80 + rand()%40); // make size-dependent?
	colorRGBA color_a0(color);
	color_a0.alpha = 0.0;

	for (unsigned i = 0; i < num_parts; ++i) {
		unsigned const plifetime(unsigned(TICKS_PER_SECOND*rand_uniform(30.0, 45.0)));
		vector3d const vel(gen_rand_vector_uniform(0.01));
		float const sz(radius*rand_uniform(0.1, 0.2));
		//type, c1, c2, lt, pos, vel, size, damage, align, coll, texture_id=-1
		gen_particle(PTYPE_GLOW, color, color_a0, plifetime, pos, vel, sz, 1000.0, ALIGN_NEUTRAL, 1);
	}
	uobject::explode(damage, bradius, etype, edir, exp_time, wclass, align, eflags, parent_);
}



