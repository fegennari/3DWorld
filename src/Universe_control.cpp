// 3D World - Universe control - code that interfaces with ships and fre_objs
// by Frank Gennari
// 11/11/05

#include "universe.h"
#include "ship.h"
#include "ship_util.h"
#include "asteroid.h"
#include "timetest.h"
#include "openal_wrap.h"
#ifdef _OPENMP
#include <omp.h>
#endif


bool const TIMETEST           = (GLOBAL_TIMETEST || 0);
bool const ORBITAL_REGEN      = 0;
bool const PRINT_OWNERSHIP    = 0;
bool const PLAYER_SLOW_PLANET_APPROACH = 1;
unsigned const GRAV_CHECK_MOD = 4; // must be a multiple of 2


float last_temp(-100.0);
string player_killer;
s_object clobj0; // closest object to player
pos_dir_up player_pdu;
cobj_vector_t const empty_cobjs; // always empty
unsigned owner_counts[NUM_ALIGNMENT] = {0};
float resource_counts[NUM_ALIGNMENT] = {0.0};


extern bool claim_planet, water_is_lava, no_shift_universe;
extern int uxyz[], window_width, window_height, do_run, fire_key, display_mode, DISABLE_WATER, frame_counter;
extern unsigned NUM_THREADS;
extern float zmax, zmin, fticks, univ_temp, temperature, atmosphere, vegetation, base_gravity, urm_static;
extern float water_h_off_rel, init_temperature, camera_shake;
extern double tfticks;
extern unsigned char **water_enabled;
extern unsigned team_credits[];
extern point universe_origin;
extern colorRGBA base_cloud_color, base_sky_color;
extern string user_text;
extern water_params_t water_params;
extern vector<free_obj *> uobjs;
extern vector<cached_obj> stat_objs;
extern vector<temp_source> temp_sources;
extern vector<hyper_inhibit_t> hyper_inhibits;
extern universe_t universe;


void process_univ_objects();
void check_shift_universe();
void draw_universe_sun_flare();
void sort_uobjects();

void setup_universe_fog(s_object const &closest);
void set_current_system_light(s_object const &clobj, point const &pspos, float a_scale, float d_scale);
void destroy_sobj(s_object const &target);
bool get_gravity(s_object &result, point pos, vector3d &gravity, int offset);
void set_universe_lighting_params(bool for_universe_draw);
point get_scaled_upt();
void add_player_ship_engine_light();


#ifdef _OPENMP
int omp_get_thread_num_3dw() {return omp_get_thread_num();} // where does this belong?
#else
int omp_get_thread_num_3dw() {return 0;}
#endif

void init_universe_display() {

	setup_ships();
	set_perspective(PERSP_ANGLE, UNIV_NCLIP_SCALE); // that's all
	check_shift_universe();
}

void set_univ_pdu() {
	u_ship const &player(player_ship());
	player_pdu = pos_dir_up(player.get_pos(), player.get_dir(), player.get_up(), 0.0, 0.0, U_VIEW_DIST);
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


bool player_near_system() {return (clobj0.system >= 0);}

bool player_inside_system() {
	if (!player_near_system()) return 0;
	ussystem const &system(clobj0.get_system());
	return (dist_less_than(player_ship().pos, system.pos, system.radius));
}


void setup_current_system(float sun_intensity) { // called in ground mode
	
	int regen_mesh(0);
	static bool last_is_lava(0);
	static int inited(0);
	static float last_water(-1.0);
	do_univ_init();

	if (!inited) {
		inited = 1;
		set_univ_pdu();
		universe.draw_all_cells(clobj0, 1, 1, 2, 0, 1); // required to gen galaxies/systems
	}
	upos_point_type const &pos(get_player_pos2());
	universe.get_object_closest_to_pos(clobj0, pos, 0, 4.0);
	colorRGBA c1(ALPHA0), c2(ALPHA0);
	float water(0.0);
	atmosphere    = 0.0;
	vegetation    = 0.0;
	water_is_lava = 0;
	base_cloud_color = WHITE;
	base_sky_color   = BACKGROUND_DAY;

	if (clobj0.type == UTYPE_PLANET || clobj0.type == UTYPE_MOON) { // planet or moon
		upos_point_type const &opos(clobj0.object->get_pos());
		float const oradius(clobj0.object->get_radius());

		if (dist_less_than(pos, opos, 10.0*oradius)) { // fairly close to the planet
			if (!dist_less_than(pos, opos, 1.01*oradius)) { // not at the surface
				upos_point_type const p_int(opos + (oradius + get_player_radius())*(pos - opos).get_norm());
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
		base_cloud_color = planet.ai_color;
		base_sky_color   = planet.ao_color;

		if (planet.lava > 0.0) {
			assert(water == 0.0);
			water = planet.lava;
			water_is_lava = 1;
		}
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
	if (water == 0.0) {
		if (DISABLE_WATER == 0) DISABLE_WATER = 2; // temporary disable
	}
	else {
		if (DISABLE_WATER == 2) DISABLE_WATER = 0; // enable
	}
	if (water != last_water && water > 0.0 && water_enabled == NULL) { // only adjust water level if a custom water_enable mask hasn't been set
		last_water = water;
		float const water_level(0.5 + 0.5*(water - 0.5)), rel_wpz(get_rel_wpz());

		if (fabs(water_level - rel_wpz) > 0.001) {
			cout << TXT(water) << TXT(water_level) << TXT(rel_wpz) << TXT(water_h_off_rel) << TXT(water_is_lava) << endl;
			change_water_level(water_level);
			regen_mesh = 2; // regen texture (is this needed?)
		}
	}
	if (water_is_lava != last_is_lava) {
		last_is_lava = water_is_lava;
		if (water_is_lava) {water_params.set_def_lava();} else {water_params.set_def_water();}
	}
	if (clobj0.type >= UTYPE_STAR) { // determine gravity
		assert(clobj0.object != NULL);
		base_gravity = 70.0*clobj0.object->gravity;
	}
	float const a_scale(0.5 + 1.0*atmosphere); // higher atmosphere = more clouds = more ambient lighting
	set_current_system_light(clobj0, pos, a_scale*sun_intensity, 0.5*sun_intensity);

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
		clear_tiled_terrain(regen_mesh == 1);
	}
}


void proc_uobjs_first_frame() {

	for (unsigned i = 0; i < uobjs.size(); ++i) {
		assert(uobjs[i]);
		uobjs[i]->first_frame_hook();
	}
}


void process_ships(int timer1) {

	sort_uobjects();
	if (TIMETEST) PRINT_TIME(" Sort uobjs");
	add_player_ship_engine_light();
	update_blasts();
	if (TIMETEST) PRINT_TIME(" Process BRs");
	process_univ_objects();
	if (TIMETEST) PRINT_TIME(" Proc Univ Objs");
	apply_explosions();
	if (TIMETEST) PRINT_TIME(" Explosion");
	u_ship &ps(player_ship());
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
		if (tratio > 1.2) {add_camera_filter(colorRGBA(1.0, 0.0, 0.0, min(0.5, max(0.1, 0.1*tratio))), 3, -1, CAM_FILT_BURN);} // NOISE_TEX
	}
}


void draw_universe_all(bool static_only, bool skip_closest, bool no_move, int no_distant, bool gen_only, bool no_asteroid_dust) {

	set_universe_lighting_params(1); // for universe drawing
	universe.get_object_closest_to_pos(clobj0, get_player_pos2(), 0, 4.0);
	if (!static_only) {setup_universe_fog(clobj0);}
	check_gl_error(120);
	universe.draw_all_cells(clobj0, skip_closest, no_move, no_distant, gen_only, no_asteroid_dust);
	check_gl_error(121);
}


void draw_universe(bool static_only, bool skip_closest, bool no_move, int no_distant, bool gen_only, bool no_asteroid_dust) { // should be process_universe()

	RESET_TIME;
	static int inited(0), first_frame_drawn(0);
	do_univ_init();
	if (!inited) {static_only = 0;} // force full universe init the first time
	create_univ_cube_map();

	if (!static_only && fire_key) {
		fire_key = 0;
		player_ship().try_fire_weapon(); // must be before process_univ_objects(), on master thread, since this can destroy objects and free VBOs
	}
	// clobj0 will not be set - need to draw cells before there are any sobjs
#ifdef _OPENMP
	// disable multiple threads when the player is away from the starting galaxy center to avoid crashing when allocating/freeing galaxies, systems, and clusters
	bool const near_init_galaxy(dist_less_than(get_player_pos2(), universe_origin, GALAXY_MIN_SIZE));

	if (inited && !static_only && NUM_THREADS > 1 && !(display_mode & 0x40) && near_init_galaxy) {
		// is this legal when a query object that tries to access a planet/moon/star through clobj as the uobject is being deleted?
		#pragma omp parallel num_threads(2)
		{
			if (omp_get_thread_num_3dw() == 1) {process_ships(timer1);}
			else {draw_universe_all(static_only, skip_closest, no_move, no_distant, gen_only, no_asteroid_dust);} // *must* be done by master thread
		}
	}
	else
#endif
	{
		if (!static_only) {process_ships(timer1);}
		draw_universe_all(static_only, skip_closest, no_move, no_distant, gen_only, no_asteroid_dust);
	}
	if (!gen_only && !first_frame_drawn) {
		proc_uobjs_first_frame();
		first_frame_drawn = 1;
	}
	if (!gen_only && !static_only) {
		if (TIMETEST) PRINT_TIME(" Universe Draw");
		check_gl_error(122);
		draw_univ_objects(); // draw free objects
		check_gl_error(123);
		if (TIMETEST) PRINT_TIME(" Free Obj Draw");
	}
	check_shift_universe();
	disable_light(get_universe_ambient_light(1)); // for universe draw
	enable_light(0);
	draw_universe_sun_flare(); // looks a bit odd, maybe it should be changed
	inited = 1;
	if (TIMETEST) PRINT_TIME(" Final Universe");
}


void proc_collision(free_obj *const uobj, upos_point_type const &cpos, point const &coll_pos, float radius, vector3d const &velocity, float mass, float elastic, int coll_tid) {

	assert(mass > 0.0);
	uobj->set_sobj_coll_tid(coll_tid);
	uobj->move_to(cpos); // setup correct position for explode?
	uobj->collision(coll_pos, velocity, S_BODY_DENSITY*mass, radius, NULL, elastic); // large mass
	uobj->move_to(cpos); // more accurate since this takes into account the terrain
}


void process_univ_objects() {

	vector<free_obj const*> stat_obj_query_res;

	for (unsigned i = 0; i < uobjs.size(); ++i) { // can we use cached_objs?
		free_obj *const uobj(uobjs[i]);
		bool const no_coll(uobj->no_coll()), particle(uobj->is_particle()), projectile(uobj->is_proj());
		if (no_coll && particle)   continue; // no collisions, gravity, or temperature on this object
		if (uobj->is_stationary()) continue;
		bool const is_ship(uobj->is_ship()), orbiting(uobj->is_orbiting());
		bool const calc_gravity(((uobj->get_time() + unsigned(size_t(uobj)>>8)) & (GRAV_CHECK_MOD-1)) == 0);
		bool const lod_coll(PLAYER_SLOW_PLANET_APPROACH && is_ship && uobj->is_player_ship()); // enable if we want to do close planet flyby
		float const radius(uobj->get_c_radius()*(no_coll ? 0.5 : 1.0));
		upos_point_type const &obj_pos(uobj->get_pos());
		vector3d gravity(zero_vector); // sum of gravity from sun, planets, possibly some moons, and possibly asteroids
		point sun_pos(all_zeros);

		// skip orbiting objects (no collisions or gravity effects, temperature is mostly constant)
		s_object clobj; // closest object
		bool const include_asteroids(!particle); // disable particle-asteroid collisions because they're too slow
		int const found_close(orbiting ? 0 : universe.get_object_closest_to_pos(clobj, obj_pos, include_asteroids, 1.0, (no_coll ? 0.0 : radius)));
		bool temp_known(0), has_rings(0);
		float limit_speed_dist(clobj.dist);

		if (found_close) {
			if (clobj.type == UTYPE_ASTEROID) {
				uasteroid const &asteroid(clobj.get_asteroid());
				float const dist_to_cobj(clobj.dist - (asteroid.radius + radius));
				uobj->set_sobj_dist(dist_to_cobj);

				if (dist_to_cobj < 0.0) { // possible collision
					upos_point_type norm(obj_pos, asteroid.pos);
					vector3d const &ascale(asteroid.get_scale());
					double const dist(norm.mag());
					if (dist > TOLERANCE) {norm /= dist;} else {norm = plus_z;} // normalize
					double const a_radius(asteroid.radius*(norm*upos_point_type(ascale)).mag()), rsum(a_radius + radius);
					
					if (dist < rsum) {
						// TODO: detailed collision?
						if (projectile) {} // projectile explosions damage the asteroid (reduce its radius? what if it's instanced?)
						float const elastic((lod_coll ? 0.1 : 1.0)*SBODY_COLL_ELASTIC);
						upos_point_type const cpos(asteroid.pos + norm*min(rsum, 1.1*dist)); // move away from the asteroid, but limit the distance to smooth the response
						proc_collision(uobj, cpos, asteroid.pos, asteroid.radius, asteroid.get_velocity(), 1.0, elastic, asteroid.get_fragment_tid(obj_pos));

						if (is_ship && clobj.asteroid_field == AST_BELT_ID) { // ship collision with asteroid belt
							//clobj.get_asteroid_belt().detach_asteroid(clobj.asteroid); // incomplete
						}
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
				float dist_to_cobj(clobj.dist - (hmap_scale*clobj_radius + radius)); // (1.0 + hmap_scale)*radius?
				
				if (dist_to_cobj > 0.0 && is_ship && clobj.has_valid_system()) {
					ussystem const &system(clobj.get_system());

					if (system.asteroid_belt) {
						// check distance to system asteroid fields (planet asteroid fields should be close enough to the planet already)
						dist_to_cobj = min(dist_to_cobj, system.asteroid_belt->get_dist_to_boundary(obj_pos));
					}
				}
				uobj->set_sobj_dist(dist_to_cobj);

				if (clobj.type == UTYPE_PLANET || clobj.type == UTYPE_MOON) {
					int coll(0);

					if (dist_to_cobj < 0.0) { // collision (except for stars)
						float coll_r;
						upos_point_type cpos;
						coll = 1;

						// player_ship and possibly other ships need the more stable but less accurate algorithm
						bool const simple_coll(!is_ship && !projectile);
						float const radius_coll(lod_coll ? 1.25*NEAR_CLIP_SCALED : radius);
						float const elastic((lod_coll ? 0.1 : 1.0)*SBODY_COLL_ELASTIC);

						if (clobj.object->collision(obj_pos, radius_coll, uobj->get_velocity(), cpos, coll_r, simple_coll)) {
							proc_collision(uobj, cpos, clobj_pos, coll_r, zero_vector, clobj.object->mass, elastic, clobj.object->get_fragment_tid(obj_pos));
							coll = 2;
						}
					} // collision
					if (is_ship) {uobj->near_sobj(clobj, coll);}
				} // planet or moon
				if (calc_gravity) {get_gravity(clobj, obj_pos, gravity, 1);}

				if (clobj.type == UTYPE_PLANET) {
					// when near a planet with rings, use the dist to the outer rings to limit speed so that we don't fly through the rings too quickly
					uplanet const &planet(clobj.get_planet());
					has_rings = (planet.ring_ro > 0.0);
					if (has_rings) {limit_speed_dist = clobj.dist - (planet.ring_ro - planet.radius);} // can be negative
				}
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
				if (dist_sq > rval*rval) continue;
				assert(ts.radius > TOLERANCE);
				float const temp(ts.temp*min(1.0f, (rval - sqrt(dist_sq))/ts.radius)*min(1.0, 0.5*max(1.0f, ts.radius/radius)));
				
				if (temp > uobj->get_temp()) {
					uobj->set_temp(temp, ts.pos, ts.source); // source should be valid (and should register as an attacker)
				}
			} // for t
			if (!orbiting) {
				float const speed_factor(uobj->get_max_sf()); // SLOW_SPEED_FACTOR = 0.04, FAST_SPEED_FACTOR = 1.0
				float speed_factor2(1.0);
				
				if (clobj.val > 0) {
					float min_sf(0.25*SLOW_SPEED_FACTOR);
					
					if (lod_coll && (clobj.type == UTYPE_PLANET || clobj.type == UTYPE_MOON)) {
						assert(clobj.object != nullptr);
						if (dot_product_ptv(upos_point_type(uobj->get_velocity()), obj_pos, clobj.object->get_pos()) < 0.0) {min_sf = (has_rings ? 0.0025 : 0.001);} // only on approach
					}
					speed_factor2 = max(min_sf, min(1.0f, 0.7f*limit_speed_dist)); // clip to [0.01, 1.0]
				}
				if (min(speed_factor, speed_factor2) > SLOW_SPEED_FACTOR) { // faster than slow speed
					for (auto h = hyper_inhibits.begin(); h != hyper_inhibits.end(); ++h) {
						float const dist_sq(p2p_dist_sq(obj_pos, h->pos));
						if (dist_sq > h->radius*h->radius) continue; // too far away to take effect
						if (uobj == h->parent) continue; // don't inhibit self
						if (h->parent->is_related(uobj)) continue; // don't inhibit our own fighters or parent
						//if (h->parent->is_enemy(uobj)) continue; // should we only inhibit enemies?
						//uobj->register_attacker(h->parent); // no attacker registration (yet)
						float const val(sqrt(dist_sq)/h->radius), val2(val*val); // 0.0 - 1.0
						min_eq(speed_factor2, ((1.0f - val2)*SLOW_SPEED_FACTOR + val2*speed_factor));
						// WRITE
					} // for h
				}
				uobj->set_speed_factor(min(speed_factor, speed_factor2));
			}
		}
	} // for i
	claim_planet = 0; // unset the flag - should have been used by this point
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

	static bool had_init_shift(0);
	if (no_shift_universe && had_init_shift) return;
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
	had_init_shift = 1;
}


void fire_planet_killer(u_ship const *const ship, point const &ship_pos, vector3d const &fire_dir, float fire_range, int obj_types) {

	point coll; // unused
	s_object target;

	// test for collisions with solid stellar objects
	line_query_state lqs;
	bool const sobj_coll((obj_types & OBJ_TYPE_UOBJ) ? universe.get_trajectory_collisions(lqs, target, coll, fire_dir, ship_pos, fire_range, 0.0) : 0);
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


// test for intersections with solid stellar objects
bool universe_intersection_test(line_query_state &lqs, point const &pos, vector3d const &dir, float range, bool include_asteroids) {

	point coll; // unused
	s_object target;
	return (universe.get_trajectory_collisions(lqs, target, coll, dir, pos, range, 0.0, include_asteroids) && target.is_solid());
}


bool universe_ray_intersect(point const &start, point const &end, int obj_types, uobject const *cur=NULL, free_obj const *ignore=NULL) {

	vector3d const dir((end - start).get_norm());
	float const range(p2p_dist(start, end));
	//if ((obj_types & OBJ_TYPE_UOBJ) && universe_intersection_test(lqs, start, dir, range)) return 1;
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
	if (!sun.is_ok()) return;
	point const offset(clobj0.get_ucell().rel_center);
	if (!univ_sphere_vis((sun.pos + offset), sun.radius)) return;
	float intensity(1.0);
	u_ship const &ps(player_ship());
	point const viewer(ps.get_pos());

	if (universe_ray_intersect(viewer, sun.pos, OBJ_TYPE_ALL, &sun, &ps)) { // center not visible, do detailed ray query
		unsigned const npts = 32;
		static point pts[npts];
		static bool pts_valid(0);
		unsigned nvis(0);
	
		for (unsigned i = 0; i < npts; ++i) {
			if (!pts_valid) {pts[i] = signed_rand_vector_norm();}
			point const pos(sun.pos + pts[i]*(1.1*sun.radius)); // move slightly away from the sun so that it doesn't intersect
			if (univ_sphere_vis((pos + offset), 0.0) && !universe_ray_intersect(viewer, pos, OBJ_TYPE_ALL, &sun, &ps)) {++nvis;}
		}
		pts_valid = 1;
		if (nvis == 0) return;
		intensity = 0.1 + 0.9*min(1.0, 2.0*nvis/float(npts)); // intensity starts to fall off when > 50% occluded
	}
	point const gv(make_pt_global(viewer));
	DoFlares(gv, (gv + ps.get_dir()), make_pt_global(sun.pos), 0.01, 0.02, 0.5*intensity, 4); // draw only some flares with half intensity
}


void send_warning_message(string const &msg, bool no_duplicate) {

	static int last_warning_tfticks(0);
	static string last_msg;
	
	if (no_duplicate) {
		if (msg == last_msg) return;
		last_msg = msg;
	}
	if ((tfticks - last_warning_tfticks) > 5.0*TICKS_PER_SECOND) {
		print_text_onscreen(msg.c_str(), RED, 1.0, 1.0*TICKS_PER_SECOND, 1);
		gen_sound(SOUND_ALERT, get_player_pos2(), 0.75);
		last_warning_tfticks = tfticks;
	}
}


void disable_player_ship() {

	if (player_ship().is_resetting()) return; // already destroyed
	print_text_onscreen("Ship Disabled", RED, 1.2, 2*TICKS_PER_SECOND, 2);
	add_camera_filter(colorRGBA(1.0, 0.5, 0.0, 0.1), 8, -1, CAM_FILT_DAMAGE);
}


void destroy_player_ship(bool captured) {

	if (player_ship().is_resetting()) return; // already destroyed
	do_run = 0;
	string msg(captured ? "Ship Captured" : "Ship Destroyed");
	if (!player_killer.empty()) {msg += " by " + player_killer;}
	print_text_onscreen(msg, RED, 1.2, 2*TICKS_PER_SECOND, 2);
	add_camera_filter(colorRGBA(1.0, 0.0, 0.0, 0.6), 8, -1, CAM_FILT_DAMAGE);
	if (!captured) gen_sound(SOUND_EXPLODE, get_player_pos2());
	player_ship().reset_after(TICKS_PER_SECOND);
}


point get_universe_display_camera_pos() {

	point const camera(get_player_pos2()), camera_scaled(camera/CELL_SIZE);
	return point(uxyz[0]+camera_scaled.x, uxyz[1]+camera_scaled.y, uxyz[2]+camera_scaled.z);
}


void draw_universe_stats() {

	char text[128];
	float const aspect_ratio((float)window_width/(float)window_height);
	u_ship const &ps(player_ship());

	// draw position and orientation
	point const cpos(get_universe_display_camera_pos());
	vector3d const dir(ps.get_dir().get_norm());
	float const player_temp(player_ship().get_true_temp());
	//sprintf(text, "Loc: (%i: %3.3f, %i: %3.3f, %i: %3.3f)  Dir: (%1.3f, %1.3f, %1.3f)", uxyz[0], camera_scaled.x, uxyz[1], camera_scaled.y, uxyz[2], camera_scaled.z, dir.x, dir.y, dir.z);
	sprintf(text, "Loc: (%3.4f, %3.4f, %3.4f)  Dir: (%1.3f, %1.3f, %1.3f)  T: %3.1f", cpos.x, cpos.y, cpos.z, dir.x, dir.y, dir.z, player_temp);
	draw_text(YELLOW, -0.009*aspect_ratio, -0.014, -0.028, text);

	// draw shields, armor, weapon status, etc.
	int const shields(int(100.0*ps.get_shields()/ps.get_max_shields()));
	int const armor(int(100.0*ps.get_armor()/ps.get_max_armor()));

	if (ps.need_ammo()) {
		int const ammo(ps.get_ammo());
		sprintf(text, "Shields %d  Armor %d  Ammo %d %s", shields, armor, ammo, ((ps.get_wnum() > 0) ? "" : "[No Weapon]"));
	}
	else {
		sprintf(text, "Shields %d  Armor %d", shields, armor);
	}
	draw_text(RED, -0.006*aspect_ratio, -0.011, -0.02, text);
	draw_health_bar(armor, shields);

	// draw credits, kills, total kills
	sprintf(text, "Credits: %u", (ps.ncredits + team_credits[ALIGN_PLAYER]));
	draw_text(WHITE, 0.0086*aspect_ratio, 0.0135, -0.025, text);
	sprintf(text, "Kills: %u/%u", ps.get_num_kills(), ps.get_tot_kills());
	draw_text(WHITE, 0.0086*aspect_ratio, 0.0126, -0.025, text);
	sprintf(text, "Deaths: %u", ps.get_num_deaths());
	draw_text(WHITE, 0.0086*aspect_ratio, 0.0117, -0.025, text);

	// draw message text (user typed)
	if (!user_text.empty()) {draw_text(WHITE, -0.010*aspect_ratio, 0.010, -0.02, user_text);} // x and z are scaled
}


void exec_universe_text(string const &text) {

	if (text.empty())   return;
	if (text[0] != '.') return; // not a command
	
	if (text == ".self-destruct") {
		player_killer = "Self Destruct";
		destroy_player_ship(0);
		return;
	}
	// *** WRITE ***
}


bool rename_obj(uobject *obj, unsigned alignment) { // a little difficult to use, have to enter the text first then click on an object

	assert(obj);
	int const owner(obj->get_owner());
	if (/*owner != NO_OWNER &&*/ owner != (int)alignment) return 0;
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
	if (!universe.get_closest_object(result, pos, UTYPE_MOON, include_asteroids, 1, 1.0, 0, 1.0, radius)) return 0;
	if (!result.object || !result.object->is_ok()) return 0; // can this happen?
	if (!result.object->sphere_intersection(pos, radius)) return 0;
	//cout << "intersects with " << result.object->get_name() << endl;
	return 1;
}


bool get_closest_object(point const &pos, s_object &result, int obj_type, bool include_asteroids, bool get_destroyed=0) {

	if (!universe.get_closest_object(result, pos, obj_type, include_asteroids, 1, 4.0, get_destroyed)) return 0;
	if (result.type != obj_type)  return 0; // incorrect type
	if (obj_type <= UTYPE_SYSTEM) return 1; // no object to check
	return (result.object != NULL && (get_destroyed || result.object->is_ok()));
}

uobject const *get_closest_world_ptr(point const &pos, int type) {

	s_object result;
	return (get_closest_object(pos, result, type, 0) ? result.object : NULL);
}


float urev_body::get_land_value(unsigned align, point const &cur_pos, float sradius) const {

	float owner_val(0.5), value(2.0*rand_float());
	
	if (!is_owned()) { // not yet owned
		unsigned const res_c(sclasses[USC_COLONY].cost + max(sclasses[USC_DEFSAT].cost, sclasses[USC_ANTI_MISS].cost));
		owner_val = ((team_credits[align] >= res_c) ? 2.5 : 0.75);
	}
	else if ((int)align == owner || align == ALIGN_PIRATE) { // owned by self or pirate
		owner_val = 1.0;
	}
	else if (TEAM_ALIGNED(align) && TEAM_ALIGNED(owner)) { // owned by enemy team
		if (resource_counts[align] > 100) {owner_val = 2.0;} // only target enemy planets once we have an economy going
		if (have_excess_credits(align)) { // excess_credits case - more aggressively go after enemy colonies
			value += 0.5*get_wealthy_value(align)*GALAXY_MIN_SIZE;
		}
		else {
			float tot_resources(0.0);
			for (unsigned i = 0; i < NUM_ALIGNMENT; ++i) {tot_resources += resource_counts[i];}
			if (tot_resources > 100.0 && resource_counts[align] > 0.5*tot_resources) {value += 0.25*GALAXY_MIN_SIZE;} // dominant resources case - more aggressive
		}
	}
	value += 2.0*liveable(); // bonus for colonies that can build ships
	value += (0.1*resources + 0.5)*owner_val;
	value += sradius/PLANET_TO_SUN_MAX_SPACING; // large system bonus
	value -= 0.5*p2p_dist(cur_pos, pos);
	return value;
}


bool line_intersect_sun(point const &p1, point const &p2, ussystem const &system, float halo) {
	// avoid choosing a destination that requires flying through a star
	return (system.sun.is_ok() && line_sphere_intersect(p1, p2, (point)system.sun.pos, halo*system.sun.radius));
}


uobject const *choose_dest_world(point const &pos, int exclude_id, unsigned align, float tmax) {

	// use starting pos rather than ship pos for selecting the closest galaxy, which should return the home galaxy even when the ship is near the edge of an adjacent galaxy
	point const galaxy_query_pos(universe_origin);
	float const g_expand(CELL_SIZE/GALAXY_MIN_SIZE); // Note: can be in more than one galaxy, but should be OK
	s_object result;
	if (!universe.get_closest_object(result, galaxy_query_pos, UTYPE_GALAXY, 0, 1, 4.0, 0, g_expand) || result.type < UTYPE_GALAXY) return nullptr; // choose closest galaxy
	ugalaxy const &galaxy(result.get_galaxy());
	unsigned const nsystems(galaxy.sols.size());
	if (nsystems == 0) return nullptr; // no systems generated (can this happen?)
	uobject const *dest = NULL;
	unsigned const nqueries(8), max_queries(nsystems/2);
	float max_pvalue(0.0);
	static rand_gen_t rgen;

	for (unsigned n = 0, ngood = 0; (ngood < nqueries && n < max_queries); ++n) {
		unsigned const six(rgen.rand()%nsystems);
		ussystem const &system(galaxy.sols[six]); // chose a random system
		float const sradius(system.get_radius());
		bool sol_good(0);

		for (auto p = system.planets.begin(); p != system.planets.end(); ++p) {
			if (p->colonizable() && p->get_id() != exclude_id && p->temp < tmax) { // planet is acceptable
				float const pvalue(p->get_land_value(align, pos, sradius));
				sol_good = 1;

				if ((max_pvalue == 0.0 || pvalue > max_pvalue) && !line_intersect_sun(pos, p->pos, system, 2.0)) { // best planet/moon so far
					max_pvalue = pvalue;
					dest       = &(*p); // choose this planet as a candidate
				}
			}
			for (auto m = p->moons.begin(); m != p->moons.end(); ++m) {
				if (!m->colonizable() || m->get_id() == exclude_id || m->temp >= tmax) continue;
				float const mvalue(m->get_land_value(align, pos, sradius));
				sol_good = 1;

				if ((max_pvalue == 0.0 || mvalue > max_pvalue) && !line_intersect_sun(pos, m->pos, system, 2.0)) { // best planet/moon so far
					max_pvalue = mvalue;
					dest       = &(*m); // choose this moon as a candidate
				}
			} // for m
		} // for p
		if (sol_good) {++ngood;}
	} // for n
	//if (dest) {cout << "align: " << align << ", name: " << dest->get_name() << ", value: " << max_pvalue << endl;}
	if (!dest) {cout << "no dest" << endl;} // testing
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
		if (world.is_owned() || !world.colonizable())            continue;
		if (check_for_land && !world.can_land_at(pos))           continue; // currently player only
		if (PRINT_OWNERSHIP) {cout << world.get_name() << " is claimed by " << get_owner_name(owner) << "." << endl;}
		float defend(world.resources);
		bool owned(0);

		while (defend > 8.0) {
			if (add_orbiting_ship(USC_DEFSAT, 0, 0, 0, own, result) == NULL) break; // try to put a defense satellite in orbit
			owned   = 1;
			defend -= 12.0;
		}
		if (defend > 0.0 || !owned) {
			owned |= (add_orbiting_ship(USC_ANTI_MISS, 0, 0, 0, own, result) != NULL); // try to put an anti-missile drone in orbit
		}
		float const rsc_val(world.resources - (homeworld ? 0.0 : 20.0));

		if (rsc_val > 0.0 && t == 0) { // planets only - colonies (first homeworld only?)
			unsigned const colony_types[5] = {USC_COLONY, USC_ARMED_COL, USC_HW_COL, USC_STARPORT, USC_HW_SPORT};
			unsigned start_val(0);
			for (start_val = 0; start_val < 4 && world.resources > 10.0f*(start_val+1); ++start_val) {}
			if (start_val == 3 || (start_val == 2 && (rand()&1))) {++start_val;}
			orbiting_ship const *oship(NULL);

			for (int i = start_val; i >= 0; --i) {
				oship = add_orbiting_ship(colony_types[i], 0, 1, 1, own, result); // put a colony on world
				
				if (oship != NULL) {
					upos_point_type const dir((own->pos - oship->pos).get_norm());
					own->move_to(oship->get_pos() + dir*(1.1*((double)own->get_c_radius() + oship->get_c_radius()))); // move away from object
					break;
				}
			}
			owned |= (oship != NULL);
		}
		if (owned) {world.set_owner(result, owner);} // only if inhabitable?
		else {assert(!world.is_owned());}
		return 1;
	}
	return 0;
}


// ************ UOBJ_SOLID/UREV_BODY ************


bool uobj_solid::collision(upos_point_type const &p, float rad, vector3d const &v, upos_point_type &cpos, float &coll_r, bool simple) const { // maybe should be in universe.cpp

	coll_r = radius;
	if (!surface_test(rad, p, coll_r, simple)) return 0;
	upos_point_type const norm(p, pos);
	double const rsum((double)coll_r + (double)rad), nmag(norm.mag());
	if (nmag > rsum) return 0;
	double const vmag(v.mag());

	if (nmag > TOLERANCE && (simple || vmag < TOLERANCE || dot_product(v, norm) > -0.1*(nmag*vmag))) { // use for shallow coll angles
		cpos = pos + norm*(rsum/nmag); // normalize, multiply, and add
	}
	else {
		get_sphere_mov_sphere_int_pt(pos, p, v, rsum, cpos);
		if (nmag > TOLERANCE) {cpos += norm*(0.05*rad/nmag);} // slight adjustment to improve stability
	}
	return 1;
}


void urev_body::get_owner_info(ostringstream &oss, bool show_uninhabited) const {

	if (is_owned()) {
		assert(unsigned(owner) < NUM_ALIGNMENT);
		oss << endl << "Owned by " << get_owner_name(owner);
	}
	else if (show_uninhabited) {oss << endl << "Uninhabited";}
}


void urev_body::set_owner(s_object const &sobj, int owner_) {
	set_owner_int(owner_);
	sobj.set_owner(owner_);
}


void urev_body::set_owner_int(int owner_) {

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

	if (is_owned()) {
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


colorRGBA urev_body::get_owner_color() const {

	if (!is_owned()) return BLACK;
	assert(unsigned(owner) < NUM_ALIGNMENT);
	return alignment_colors[owner];
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
	urev_body const &world(clobj.get_world());

	if (is_player_ship()) {
		//if (coll == 2) {
		if (claim_planet && (coll || dist_less_than(pos, world.get_pos(), 2.0*radius))) {
			bool const homeworld(clobj.type == UTYPE_PLANET && !has_homeworld());

			if (check_dest_ownership(world.get_id(), pos, this, 1, homeworld)) {
				if (homeworld) {claim_world(&world);}
				string const str(string("You have claimed ") + world.get_name() + (homeworld ? " as your homeworld" : ""));
				print_text_onscreen(str, CYAN, 1.2, 2*TICKS_PER_SECOND);
				gen_sound((homeworld ? (int)SOUND_POWERUP : (int)SOUND_ITEM), get_player_pos2());
			}
		}
	}
	else if (world.get_owner() == NO_OWNER && world.colonizable() && can_colonize() && // or last orbiting ship is exploding?
		!is_fighter() && !is_orbiting() && !is_rand_spawn() && can_move() && have_resources_to_colonize(alignment))
	{
		float const odist_sq(p2p_dist_sq(pos, world.get_pos()));

		if (coll || !dest_mgr.is_valid() || dest_mgr.is_cur_obj(world.get_id()) || odist_sq < 0.05*p2p_dist_sq(pos, dest_mgr.get_pos())) {
			if (coll || target_obj == NULL || odist_sq < 0.5*p2p_dist_sq(pos, target_obj->get_pos())) {
				dest_override = 1;
				dest_mgr.set_object(&world); // the original destination is lost (if there was one)
			}
		}
	}
}


// ************ ORBITING_SHIP ************


orbiting_ship *add_orbiting_ship(unsigned sclass, bool guardian, bool on_surface, bool pos_from_parent, free_obj const *parent, s_object const &world_path) {

	assert(parent);
	urev_body &world(world_path.get_world());
	//assert(world.get_owner() == parent->get_align());
	assert(sclass < sclasses.size());
	if (parent->get_temp() > 1.0*TEMP_FACTOR*sclasses[sclass].max_t) return NULL; // planet/moon is too hot
	if (!alloc_resources_for(sclass, parent->get_align(), 0)) return NULL;
	world.inc_orbiting_refs();
	float angle(rand_uniform(0.0, TWO_PI));
	point const parent_pos(parent->get_pos());

	if (pos_from_parent) { // set initial location based on position of ship relative to planet
		vector3d dir((parent->get_pos() - world.get_pos()).get_norm());
		world.rotate_vector(dir);
		angle = atan2(dir.y, dir.x);
	}

	// not sure what to set orbit_radius to - could collide with other orbiting objects
	float const orbit_radius(on_surface ? 0.0 : 2.0*world.get_radius()), rate(0.0);
	point const start_pos((pos_from_parent && on_surface) ? point(parent->get_pos() - world.get_pos()) : all_zeros);
	orbiting_ship *const ship(new orbiting_ship(sclass, parent->get_align(), guardian, world_path, zero_vector, start_pos, orbit_radius, angle, rate));
	ship->set_parent(parent);
	add_uobj_ship(ship);
	return ship;
}


// geostationary orbit (GEO) and geosynchronous orbit (GSO)
orbiting_ship::orbiting_ship(unsigned sclass_, unsigned align, bool guardian, s_object const &world_path,
							 vector3d const &axis_, point const &start_pos, float rad, float start_ang, float rate)
	: u_ship(sclass_, all_zeros, align, (AI_ATT_ENEMY | (guardian ? AI_GUARDIAN : 0)), TARGET_CLOSEST, 0), GSO(rate == 0.0),
	  fixed_pos(start_pos != all_zeros), has_sobj(0), has_decremented_owner(0), last_build_time(time), system_ix(world_path.system),
	planet_ix(world_path.planet), moon_ix(world_path.moon), orbit_r(rad), rot_rate(rate), start_angle(start_ang), axis(axis_), rel_pos(start_pos)
{
	urev_body &world(world_path.get_world());
	assert(world.type == UTYPE_PLANET || world.type == UTYPE_MOON);
	assert(system_ix >= 0 && planet_ix >= 0);
	if (world.type == UTYPE_MOON) {assert(moon_ix >= 0);}
	ustar const &sun(world_path.get_system().sun);
	sobj_liveable = world.liveable();
	orbiting_type = world.type;
	assert(!can_move());
	homeworld.set_object(&world);
	flags |= OBJ_FLAGS_ORBT;

	if (orbit_r == 0.0) {
		assert(specs().mass > 0.0);
		// *** determine radius of orbit orbit_r based on mass? ***
	}
	else {
		world.get_valid_orbit_r(orbit_r, c_radius);
	}
	assert(radius < world.radius); // too strict?
	assert(orbit_r == 0.0 || orbit_r >= world.radius + c_radius); // too strict?
	if (orbit_r == 0.0) {orbit_r = world.radius + 0.65*radius;} // + radius?
	if (GSO || axis == zero_vector) {axis = world.rot_axis;} // GEO (orbits along the equator)

	if (fixed_pos) {
		world.rotate_vector(rel_pos); // convert to local object space
		rel_pos.normalize();
	}
	set_pos_from_sobj(&world);
	sobj_radius = world.radius;
	sun_pos     = sun.pos;
	sun_energy  = sun.get_energy();
}


void orbiting_ship::update_state() {

	if (!is_ok()) return;
	bool const exploding_now(explode_now());
	if ((flags & OBJ_FLAGS_DIST) && !exploding_now && (time&31) != 0) return; // don't update this frame when far away
	s_object result;
	if (!get_closest_object(pos, result, UTYPE_GALAXY, 0)) return; // galaxy not found (player is in a different cell?)
	ugalaxy &galaxy(result.get_galaxy());
	if (!galaxy.gen || galaxy.sols.empty()) return; // systems not generated for this galaxy
	assert((unsigned)system_ix < galaxy.sols.size());
	ussystem &system(galaxy.sols[system_ix]);
	if (system.planets.empty()) return; // planets not generated for this system
	assert((unsigned)planet_ix < system.planets.size());
	uplanet &planet(system.planets[planet_ix]);
	urev_body *world(nullptr);
	
	if (orbiting_type == UTYPE_PLANET) {world = &planet;}
	else {
		assert(orbiting_type == UTYPE_MOON);
		if (planet.moons.empty()) return; // moons not generated for this planet
		assert((unsigned)moon_ix < planet.moons.size());
		world = &planet.moons[moon_ix];
	}
	result.object = world;
	//velocity = zero_vector; // ???

	if (!world->is_ok()) { // was destroyed
		if (!is_exploding()) {destroy_ship(0.0);}
		return;
	}
	if (homeworld.update_pos_if_close(world)) { // the object we are orbiting is still there
		set_pos_from_sobj(world);
	}
	// release ownership on explosion
	//if (!ORBITAL_REGEN && exploding_now && world->get_owner() == (int)alignment) {world->dec_orbiting_refs(result);} // will die this frame
	if (!ORBITAL_REGEN && is_exploding() && !has_decremented_owner && world->get_owner() == (int)alignment) { // when exploding, only once
		world->dec_orbiting_refs(result);
		has_decremented_owner = 1;
	}
	else if (!is_exploding() && !world->is_owned()) {
		world->set_owner(result, alignment); // have to reset - world must have been regenerated
	}
	if (world->temp > get_temp()) {set_temp(FOBJ_TEMP_SCALE*world->temp, world->get_pos(), NULL);}
}


void orbiting_ship::set_pos_from_sobj(urev_body const *const sobj) {

	vector3d delta(1.0, 0.0, 0.0); // should delta start out in x or y?

	if (fixed_pos) { // doesn't have to be at the equator (colony)
		delta = rel_pos;
		sobj->rotate_vector_inv(delta);
		delta.normalize();
	}
	else {
		angle = (GSO ? (start_angle + sobj->rot_ang/TO_DEG) : ((double)angle + fticks*(double)rot_rate));
		rotate_norm_vector3d_into_plus_z(axis, delta, -1.0); // switch to local coordinate system (inverse rotate)
		rotate_vector3d_norm (axis, -angle, delta); // negate angle?
	}
	sobj_pos  = sobj->pos;
	pos       = sobj->pos + delta*double(orbit_r);
	reset_pos = pos;
	if (specs().max_turn  == 0.0) {dir = delta;} // face outward (no turning?)
	if (specs().roll_rate == 0.0) {orthogonalize_dir(axis, dir, upv, 1);}
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


float const *uobject::get_sphere_shadow_pmap(point const &sun_pos, point const &obj_pos, int ndiv) const {

	if (!has_custom_shadow_profile()) return NULL;
	assert(ndiv >= 3);
	float const dist_to_sun(p2p_dist(pos, sun_pos)), shadow_scale_val((dist_to_sun + p2p_dist(pos, obj_pos))/dist_to_sun);
	static vector<float> pmap_vector;
	pmap_vector.resize(ndiv);
	point const ce[2] = {pos, sun_pos};
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius, 0.0, ndiv));

	for (unsigned i = 0; i < (unsigned)ndiv; ++i) { // assumes the cylinder is more or less constant radius
		pmap_vector[i] = shadow_scale_val*(get_radius_at(vpn.p[i<<1]) - radius);
	}
	return &pmap_vector.front();
}


void ustar::explode(float damage, float bradius, int etype, vector3d const &edir, int exp_time, int wclass,
					int align, unsigned eflags, free_obj const *parent_)
{
	unsigned const num_parts(80 + rand()%40); // make size-dependent?
	colorRGBA color_a0(color, 0.0);

	for (unsigned i = 0; i < num_parts; ++i) {
		unsigned const plifetime(unsigned(TICKS_PER_SECOND*rand_uniform(30.0, 45.0)));
		vector3d const vel(gen_rand_vector_uniform(0.01));
		float const sz(radius*rand_uniform(0.1, 0.2));
		//type, c1, c2, lt, pos, vel, size, damage, align, coll, texture_id=-1
		gen_particle(PTYPE_GLOW, color, color_a0, plifetime, pos, vel, sz, 1000.0, ALIGN_NEUTRAL, 1);
	}
	uobject::explode(damage, bradius, etype, edir, exp_time, wclass, align, eflags, parent_);
	colorRGBA ci(colorA), co(colorB);
	blend_color(co, co, RED, 0.5); // red shift outer color
	add_uparticle_cloud(pos, 1.0*radius, 10.0*radius, ci, co, colorRGBA(ci, 0.0), colorRGBA(co, 0.0), 60*TICKS_PER_SECOND, 1000.0, 0.3, 0.4);
}


struct owner_stats_t {
	unsigned np, nm;
	float res, pop;
	owner_stats_t() : np(0), nm(0), res(0.0), pop(0.0) {}
};

void print_univ_owner_stats() {

	owner_stats_t stats[NUM_ALIGNMENT];
	unsigned np(0), nm(0);

	for (unsigned i = 0; i < U_BLOCKS; ++i) { // z
		for (unsigned j = 0; j < U_BLOCKS; ++j) { // y
			for (unsigned k = 0; k < U_BLOCKS; ++k) { // x
				int const ii[3] = {(int)k, (int)j, (int)i};
				ucell const &cell(universe.get_cell(ii));
				if (cell.galaxies == nullptr) continue;
				for (unsigned g = 0; g < cell.galaxies->size(); ++g) {
					ugalaxy const &galaxy((*cell.galaxies)[g]);
					for (unsigned s = 0; s < galaxy.sols.size(); ++s) {
						ussystem const &system(galaxy.sols[s]);
						for (unsigned p = 0; p < system.planets.size(); ++p) {
							uplanet const &planet(system.planets[p]);
							for (unsigned m = 0; m < planet.moons.size(); ++m) {
								umoon const &moon(planet.moons[m]);
								++nm;
								if (moon.owner < 0) continue;
								assert(moon.owner < NUM_ALIGNMENT);
								owner_stats_t &os(stats[moon.owner]);
								os.res += moon.resources;
								os.pop += moon.population;
								++os.nm;
							}
							++np;
							if (planet.owner < 0) continue;
							assert(planet.owner < NUM_ALIGNMENT);
							owner_stats_t &os(stats[planet.owner]);
							os.res += planet.resources;
							os.pop += planet.population;
							++os.np;
						} // for p
					} // for s
				} // for g
			} // for k
		} // for j
	} // for i
	unsigned maxvals[3] = {1000000, 1000000, 1000000}; // starting values set the min column width

	for (unsigned i = 0; i < NUM_ALIGNMENT; ++i) {
		max_eq(maxvals[0], stats[i].np);
		max_eq(maxvals[1], stats[i].nm);
		max_eq(maxvals[2], (unsigned)stats[i].res);
	}
	update_maxvals(maxvals, 3);
	cout << endl << "Total Planets    : " << np << "  Total Moons: " << nm << endl;
	cout << "Stats            : <planets> <moons> <resources> <population>" << endl;

	for (unsigned i = 0; i < NUM_ALIGNMENT; ++i) {
		owner_stats_t const &os(stats[i]);
		if (os.np == 0 && os.nm == 0) continue; // no resources
		cout << align_names[i];
		print_n_spaces(18 - (int)align_names[i].size());
		cout << ": ";
		write_uint_pad(os.np,  maxvals[0], "   ");
		write_uint_pad(os.nm,  maxvals[1]);
		write_uint_pad(os.res, maxvals[2], "     ");
		if (os.pop > 1000) {cout << (int)(os.pop)/1000.0 << " B";} // in billions
		else {cout << (int)os.pop << " M";} // in millions
		cout << endl;
	}
	cout << endl;
}



