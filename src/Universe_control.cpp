// 3D World - Universe control - code that interfaces with ships and fre_objs
// by Frank Gennari
// 11/11/05

#include "universe.h"
#include "ship.h"
#include "ship_util.h"
#include "timetest.h"


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
extern int uxyz[], window_width, window_height, do_run, fire_key, display_mode, DISABLE_WATER;
extern float zmax, zmin, fticks, univ_temp, temperature, atmosphere, vegetation, base_gravity, urm_static;
extern float def_water_level, water_plane_z, tan_term, sin_term, init_temperature, camera_shake;
extern unsigned team_credits[];
extern point last_camera;
extern string user_text;
extern vector<free_obj *> uobjs;
extern vector<cached_obj> stat_objs;
extern vector<temp_source> temp_sources;
extern universe_t universe;


void process_univ_objects();
void check_shift_universe();


void setup_universe_fog(s_object const &closest, bool damaged);
void set_current_system_light(s_object const &clobj, point const &pspos, float a_scale, float d_scale);
void destroy_sobj(s_object const &target);
bool get_gravity(s_object &result, point pos, vector3d &gravity, int offset);
void set_lighting_params();
point get_scaled_upt();



void init_universe_display() {

	setup_ships();
	set_perspective(PERSP_ANGLE, UNIV_NCLIP_SCALE); // that's all (closer near_clip for univ_planet_lod?)
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
	rseed1 = rseed2 = 1;
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
		universe.draw_all_cells(clobj0, 1, 1, 1); // required to gen galaxies/systems
	}
	point const &pos(get_player_pos2());
	universe.get_object_closest_to_pos(clobj0, pos);
	set_current_system_light(clobj0, pos, 0.5, 0.5);
	colorRGBA c1(ALPHA0), c2(ALPHA0);
	float water(0.0);
	atmosphere = 0.0;
	vegetation = 0.0;

	if (univ_temp != last_temp) { // determine temperature
		cout << "temperature: " << univ_temp << endl;
		last_temp   = univ_temp;
		temperature = univ_temp;
		//init_temperature = univ_temp; // ???
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
			float const water_rel(calc_glaciated_rel_value(water));
			
			if (water_rel != def_water_level || water_rel != water_plane_z) {
				//cout << "water: " << def_water_level << "/" << water_plane_z << " => " << water << "/" << water_rel << endl;
				//def_water_level = water_plane_z = water_rel; // *** FIX ***
				//regen_mesh      = 1; // regen texture?
				// *** recalculate water ***
			}
		}
	}
	if (clobj0.type >= UTYPE_STAR) { // determine gravity
		assert(clobj0.object != NULL);
		base_gravity = 70.0*clobj0.object->gravity;
	}
	//if (regen_mesh) *** regen mesh ***
	setup_landscape_tex_colors(c1, c2);
}


void draw_universe(bool static_only, bool skip_closest, bool no_distant) { // should be process_universe()

	RESET_TIME;
	static int inited(0);
	set_lighting_params();
	WHITE.do_glColor();
	do_univ_init();

	if (!inited) { // clobj0 will not be set - need to draw cells before there are any sobjs
		inited      = 1;
		static_only = 0; // force full universe init the first time
	}
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
			damaged = 1; // adds lots of red fog to the player's view while dying
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
		universe.get_object_closest_to_pos(clobj0, get_player_pos2());
		setup_universe_fog(clobj0, (burning || damaged));
	}
	else {
		universe.get_object_closest_to_pos(clobj0, get_player_pos2());
	}
	glEnable(GL_COLOR_MATERIAL);
	check_gl_error(120);
	universe.draw_all_cells(clobj0, skip_closest, skip_closest, no_distant);
	check_gl_error(121);

	if (!static_only) {
		if (TIMETEST) PRINT_TIME(" Universe Draw");
		check_gl_error(122);
		draw_univ_objects(all_zeros); // draw free objects
		check_gl_error(123);
		if (TIMETEST) PRINT_TIME(" Free Obj Draw");
		check_shift_universe();
	}
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(get_universe_ambient_light());
	glEnable(GL_LIGHT0);
	if (TIMETEST) PRINT_TIME(" Final Universe");
}


void process_univ_objects() {

	vector<free_obj const*> stat_obj_query_res;

	for (unsigned i = 0; i < uobjs.size(); ++i) { // can we use cached_objs?
		free_obj *const uobj(uobjs[i]);
		bool const no_coll(uobj->no_coll());
		if (no_coll && uobj->is_particle()) continue; // no collisions, gravity, or temperature on this object
		if (uobj->is_stationary())          continue;
		bool const is_ship(uobj->is_ship()), orbiting(uobj->is_orbiting());
		bool const calc_gravity(((uobj->get_time() + (unsigned(uobj)>>8)) & (GRAV_CHECK_MOD-1)) == 0);
		bool const lod_coll(univ_planet_lod && is_ship && uobj->is_player_ship());
		float const radius(uobj->get_c_radius()*(no_coll ? 0.5 : 1.0));
		upos_point_type const &obj_pos(uobj->get_pos());
		vector3d gravity(zero_vector); // sum of gravity from sun, planets, possibly some moons, and possibly asteroids

		// skip orbiting objects (no collisions or gravity effects, temperature is mostly constant)
		s_object clobj; // closest object
		int const found_close(orbiting ? 0 : universe.get_object_closest_to_pos(clobj, obj_pos));

		if (found_close) {
			float const clobj_radius(clobj.object->get_radius());
			point const clobj_pos(clobj.object->get_pos());
			float const temperature(universe.get_point_temperature(clobj, obj_pos)*(5 - uobj->get_shadow_val())); // shadow_val = 0-3
			uobj->set_temp(temperature, clobj_pos);
			float const dist_to_cobj(clobj.dist - (HMAP_SCALE*clobj_radius + radius)); // (1.0 + HMAP_SCALE)*radius?
			uobj->set_sobj_dist(dist_to_cobj);
			int coll(0);

			if (clobj.type == UTYPE_PLANET || clobj.type == UTYPE_MOON) {
				if (dist_to_cobj < 0.0) { // collision (except for stars)
					assert(clobj.object != NULL);
					float coll_r;
					point cpos;
					coll = 1;

					// player_ship and possibly other ships need the more stable but less accurate algorithm
					bool const simple_coll(!is_ship && !uobj->is_proj());
					float const radius_coll(lod_coll ? 1.25*NEAR_CLIP_SCALED : radius);
					float const elastic((lod_coll ? 0.1 : 1.0)*SBODY_COLL_ELASTIC);

					if (clobj.object->collision(obj_pos, radius_coll, uobj->get_velocity(), cpos, coll_r, simple_coll)) {
						assert(clobj.object->mass > 0.0);
						uobj->set_sobj_coll();
						uobj->move_to(cpos); // setup correct position for explode?
						uobj->collision(clobj_pos, zero_vector, S_BODY_DENSITY*clobj.object->mass, coll_r, NULL, elastic);
						uobj->move_to(cpos); // more accurate since this takes into account the terrain
						coll = 2;
					}
				} // collision
				if (is_ship) uobj->near_sobj(clobj, coll);
			} // planet or moon
			if (calc_gravity) get_gravity(clobj, obj_pos, gravity, 1);
		} // found_close
		else {
			uobj->set_temp(0.0, all_zeros);
		}
		if (calc_gravity) {
			bool near_b_hole(0);

			if (!stat_objs.empty()) {
				all_query_data qdata(&stat_objs, obj_pos, 10.0, urm_static, uobj, stat_obj_query_res);
				get_all_close_objects(qdata);
				
				for (unsigned j = 0; j < stat_obj_query_res.size(); ++j) { // asteroid gravity
					near_b_hole |= (stat_obj_query_res[j]->get_gravity(gravity, obj_pos) == 2);
				}
			}
			uobj->add_gravity(gravity, float(GRAV_CHECK_MOD), near_b_hole);
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
	last_camera = get_player_pos2();
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
			if ((camera[d] < sign*CELL_SIZEo2) ^ (sign > 0)) {
				move[d] -= sign*CELL_SIZE;
				uxyz[d] += sign;
				sh[d]    = sign;
				moved    = 1;
				universe.shift_cells(sh[0], sh[1], sh[2]);
			}
		}
	}
	if (moved) { // advance all free objects by a cell
		shift_univ_objs(move, 1);
		last_camera += move;
	}
}


void fire_planet_killer(u_ship const *const ship, point const &ship_pos, vector3d const &fire_dir, float fire_range, int obj_types) {

	point coll; // unused
	s_object target;

	// test for collisions with sloid stellar objects
	bool const sobj_coll((obj_types & OBJ_TYPE_UOBJ) ?
		universe.get_trajectory_collisions(target, coll, fire_dir, ship_pos, fire_range, 0.0) : 0);
	bool sobjc(sobj_coll && target.is_solid());

	// test for other ship collisions
	int const itype(obj_types & ~OBJ_TYPE_UOBJ);
	free_obj *fobj = NULL;
	uobject *obj = NULL;
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


void send_warning_message(string const &msg) {

	print_text_onscreen(msg.c_str(), RED, 1.0, 3*TICKS_PER_SECOND/2, 1);
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
	player_ship().reset_after(TICKS_PER_SECOND);
}


void draw_universe_stats() {

	char text[128];
	float const aspect_ratio((float)window_width/(float)window_height);
	u_ship const &ps(player_ship());

	// draw position and orientation
	point const camera(ps.get_pos()), camera_scaled(camera/CELL_SIZE);
	vector3d const dir(ps.get_dir().get_norm());
	YELLOW.do_glColor();
	sprintf(text, "Loc: (%i: %3.3f, %i: %3.3f, %i: %3.3f)  Dir: (%1.3f, %1.3f, %1.3f)",
		uxyz[0], camera_scaled.x, uxyz[1], camera_scaled.y, uxyz[2], camera_scaled.z, dir.x, dir.y, dir.z);
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
	string const str(string("Renaming ") + obj->get_name() + " to " + user_text);
	print_text_onscreen(str, PURPLE, 0.8, TICKS_PER_SECOND);
	obj->rename(user_text);
	return 1;
}


bool get_closest_object(point const &pos, s_object &result, int obj_type, bool get_destroyed=0) {

	if (!universe.get_largest_closest_object(result, pos, 0, obj_type, 1, 4.0, get_destroyed)) return 0;
	return (result.type == obj_type && (result.object->is_ok() || get_destroyed));
}


uobject const *get_closest_world_ptr(point const &pos, int type) {

	s_object result;
	return (get_closest_object(pos, result, type) ? result.object : NULL);
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


uobject const *choose_dest_world(point const &pos, int exclude_id, unsigned align) {

	s_object result;
	float const expand(CELL_SIZE/GALAXY_MIN_SIZE); // Note: can be in more than one galaxy, but should be OK
	if (!universe.get_largest_closest_object(result, pos, 0, UTYPE_GALAXY, 1, expand) || result.type != UTYPE_GALAXY) return NULL;
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
					if (planet_acceptable) {
						max_pvalue = pvalue;
						max_svalue = svalue;
						dest       = &planet;
					}
					for (unsigned k = 0; k < planet.moons.size(); ++k) {
						umoon const &moon(planet.moons[k]);
						if (!moon.colonizable() || moon.get_id() == exclude_id) continue;
						float const mvalue(moon.get_land_value(align, pos, sradius));

						if (max_pvalue == 0.0 || mvalue > max_pvalue) {
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


bool check_dest_ownership(int uobj_id, point const &pos, free_obj const *own, bool check_for_land, bool homeworld) {
	
	assert(uobj_id >= 0 && own); // obj may be an invalid pointer though
	unsigned const owner(own->get_align());
	s_object result;
	int const world_types[2] = {UTYPE_PLANET, UTYPE_MOON};

	for (unsigned t = 0; t < 2; ++t) {
		if (!get_closest_object(pos, result, world_types[t])) continue;
		if (result.object->get_id() != uobj_id)               continue;
		urev_body &world(result.get_world());
		if (world.owner != NO_OWNER || !world.colonizable())  continue;
		if (check_for_land && !world.can_land_at(pos))        continue; // currently player only
		if (PRINT_OWNERSHIP) cout << world.get_name() << " is claimed by " << get_owner_name(owner) << "." << endl;
		float defend(world.resources);
		bool owned(0);

		while (defend > 8.0) {
			owned  |= add_orbiting_ship(USC_DEFSAT, 0, 0, 0, own, &world); // put a defense satellite in orbit
			defend -= 12.0;
		}
		if (defend > 0.0 || !owned) {
			owned  |= add_orbiting_ship(USC_ANTI_MISS, 0, 0, 0, own, &world); // put an anti-missile drone in orbit
		}
		float const rsc_val(world.resources - (homeworld ? 0.0 : 20.0));

		if (rsc_val > 0.0 && t == 0) { // planets only - colonies (first homeworld only?)
			unsigned const colony_types[5] = {USC_COLONY, USC_ARMED_COL, USC_HW_COL, USC_STARPORT, USC_HW_SPORT};
			bool placed(0);
			unsigned start_val(0);
			for (start_val = 0; start_val < 4 && world.resources > 10.0*(start_val+1); ++start_val) {}
			if (start_val == 3 || (start_val == 2 && (rand()&1))) ++start_val;

			for (int i = start_val; i >= 0 && !placed; --i) {
				placed = add_orbiting_ship(colony_types[i], 0, 1, 1, own, &world); // put a colony on world if homeworld
			}
			owned |= placed;
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

	int const new_owner(sobj.get_owner());
	//if (new_owner != NO_OWNER && owner == NO_OWNER) ++orbiting_refs; // ???
	set_owner_int(new_owner);
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

	if (is_player_ship()) {
		if (!univ_planet_lod && coll == 2) {
			bool const homeworld(clobj.type == UTYPE_PLANET && !has_homeworld());

			if (check_dest_ownership(clobj.object->get_id(), pos, this, 1, homeworld)) {
				if (homeworld) claim_world(clobj.object);
				string const str(string("You have claimed ") + clobj.object->get_name() + (homeworld ? " as your homeworld" : ""));
				print_text_onscreen(str, CYAN, 1.2, 2*TICKS_PER_SECOND);
			}
		}
	}
	else if (clobj.object->get_owner() == NO_OWNER && clobj.get_world().colonizable() &&
		!is_fighter() && !is_orbiting() && can_move() && have_resources_to_colonize(alignment))
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


bool add_orbiting_ship(unsigned sclass, bool guardian, bool on_surface, bool pos_from_parent,
					   free_obj const *parent, urev_body *obj)
{
	assert(obj && parent);
	//assert(obj->get_owner() == parent->get_align());
	if (!alloc_resources_for(sclass, parent->get_align(), 0)) return 0;
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
	point const start_pos((pos_from_parent && on_surface) ? (parent->get_pos() - obj->get_pos()) : all_zeros);
	orbiting_ship *const ship(new orbiting_ship(sclass, parent->get_align(), guardian, obj, zero_vector, start_pos,
		orbit_radius, angle, rate));
	ship->set_parent(parent);
	add_uobj_ship(ship);
	return 1;
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

		if (get_closest_object(pos, result, orbiting_type, 1)) {
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
				world.set_owner(result, alignment); // have to reset - world must have been regenerated (untested)
			}
			if (world.temp > get_temp()) set_temp(world.temp, world.get_pos(), NULL);
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


cobj_vector_t const &uobject::get_cobjs() const {
	
	return empty_cobjs;
}


void uobject::gen_fragments() const {

	unsigned const num_fragments((rand()&7) + 8);

	for (unsigned i = 0; i < num_fragments; ++i) {
		add_uobj(new uobj_asteroid((pos + signed_rand_vector(1.2*radius)),
			0.2*radius*rand_uniform(0.5, 1.0), AS_MODEL_SPHERE, (10*TICKS_PER_SECOND + rand()%TICKS_PER_SECOND))); // temporary
	}
	gen_moving_fragments(pos, unsigned(rand_uniform(40, 60)), 1.0);
}


void uobject::gen_moving_fragments(point const &hit_pos, unsigned num, float rscale, float vscale) const {

	colorRGBA const &pcolor(/*texture_color(ROCK_SPHERE_TEX)*/WHITE);

	for (unsigned i = 0; i < num; ++i) {
		unsigned const ltime(5*TICKS_PER_SECOND + rand()%TICKS_PER_SECOND);
		point ppos(hit_pos + signed_rand_vector(0.5*rscale*radius));
		ppos += (ppos - pos).get_norm()*0.1*radius;
		vector3d const vel(((ppos - pos).get_norm() + signed_rand_vector(0.25))*radius*vscale*0.02);
		gen_particle(PTYPE_SPHERE, pcolor, pcolor, ltime, ppos, vel, rscale*radius*rand_uniform(0.05, 0.1), 0.0, ALIGN_NEUTRAL, 1, ROCK_SPHERE_TEX);
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



