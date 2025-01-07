// 3D World - Building Animals (rats, spiders, snakes, flies, cockroaches, etc.)
// by Frank Gennari 1/16/22

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h" // for object_model_loader_t
#include "openal_wrap.h"

bool  const SPIDERS_FOLLOW_PLAYER = 0;
float const RAT_FOV_DEG      = 60.0; // field of view in degrees
float const RAT_VIEW_FLOORS  = 4.0; // view distance in floors
float const RAT_FEAR_SPEED   = 1.3; // multiplier
float const RAT_ATTACK_SPEED = 1.2; // multiplier
float const RAT_FOV_DP(cos(0.5*RAT_FOV_DEG*TO_RADIANS));
float const SPIDER_VIEW_FLOORS = 4.0; // view distance in floors

extern bool player_attracts_flies, player_wait_respawn, camera_in_building, player_in_tunnel;
extern int animate2, camera_surf_collide, frame_counter, display_mode;
extern float fticks, NEAR_CLIP;
extern double tfticks, camera_zh;
extern building_dest_t cur_player_building_loc;
extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;
extern building_t const *player_building;

float get_closest_building_sound(point const &at_pos, point &sound_pos, float floor_spacing);
sphere_t get_cur_frame_loudest_sound();
void apply_building_gravity(float &vz, float dt_ticks);
void apply_fc_cube_max_merge_xy(vect_cube_t &cubes);
void register_fly_attract(bool no_msg);
bool phone_is_ringing();
unsigned get_fishtank_coll_cubes(room_object_t const &c, cube_t cubes[7]);
template<typename T> bool line_int_cubes_exp(point const &p1, point const &p2, vector<T> const &cubes, vector3d const &expand, cube_t const &line_bcube);


void building_animal_t::sleep_for(float time_secs_min, float time_secs_max) {
	wake_time = (float)tfticks + rand_uniform(time_secs_min, time_secs_max)*TICKS_PER_SECOND;
	dist_since_sleep = 0.0; // reset the counter
}
float building_animal_t::move(float timestep, bool can_move_forward, float anim_time_scale) { // returns movement distance
	// update animation time using position change; note that we can't just do the update in the rat movement code below because pos may be reset in case of collision
	anim_time += anim_time_scale*p2p_dist(pos, last_pos)/radius; // scale with size so that small rat/spider legs move faster
	last_pos   = pos;
	if (anim_time > 100000.0) {anim_time = 0.0;} // reset animation after awhile to avoid precision problems; should this be done when the player isn't looking?
	float move_dist(0.0);

	if (is_sleeping()) {
		if ((float)tfticks > wake_time) {wake_time = speed = 0.0;} // time to wake up
	}
	else if (speed == 0.0) {
		anim_time = 0.0; // reset animation to rest pos
	}
	else if (can_move_forward) { // apply movement
		move_dist = timestep*speed;
		pos = pos + move_dist*dir;
		dist_since_sleep += move_dist;
	}
	return move_dist;
}

void building_t::update_animals(point const &camera_bs, unsigned building_ix) {
	if (!animate2 || is_rotated() || !has_room_geom() || interior->rooms.empty()) return;
	update_rats         (camera_bs, building_ix);
	update_sewer_rats   (camera_bs, building_ix);
	update_pet_rats     (camera_bs, building_ix);
	update_spiders      (camera_bs, building_ix);
	update_sewer_spiders(camera_bs, building_ix);
	update_snakes       (camera_bs, building_ix);
	update_insects      (camera_bs, building_ix);
	interior->room_geom->last_animal_update_frame = frame_counter;
}

float select_animal_radius(float sz_min, float sz_max, float floor_spacing, rand_gen_t &rgen) {
	// there's no error check for min_sz <= max_sz, so just use min_sz in that case
	return 0.5f*floor_spacing*((sz_min >= sz_max) ? sz_min : rgen.rand_uniform(sz_min, sz_max));
}
template<typename T> void building_t::add_animals_on_floor(T &animals, unsigned building_ix, unsigned num_min, unsigned num_max, float sz_min, float sz_max) const {
	if (animals.placed) return; // already placed
	animals.placed = 1; // even if there were no animals placed
	if (num_max == 0) return; // none to place
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, mat_ix+1); // unique per building
	float const floor_spacing(get_window_vspace());
	unsigned const num(num_min + ((num_min == num_max) ? 0 : (rgen.rand() % (num_max - num_min + 1))));
	animals.reserve(num);

	for (unsigned n = 0; n < num; ++n) {
		float const radius(select_animal_radius(sz_min, sz_max, floor_spacing, rgen));
		point const pos(gen_animal_floor_pos(radius, T::value_type::allow_in_attic(), 0, 0, T::value_type::not_by_ext_door(), rgen)); // not_player_visible=0, pref_dark_room=0
		if (pos == all_zeros) continue; // bad pos? skip this animal
		animals.add(typename T::value_type(pos, radius, rgen.signed_rand_vector_spherical_xy_norm(), n));
	}
}

point building_t::gen_animal_floor_pos(float radius, bool place_in_attic, bool not_player_visible, bool pref_dark_room, bool not_by_ext_door, rand_gen_t &rgen) const {
	vector3d const xlate(get_tiled_terrain_model_xlate()); // too difficult to pass this in through the call stack, so just recalculate
	point const camera_bs(camera_pdu.pos - xlate);
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());

	for (unsigned n = 0; n < 100; ++n) { // make up to 100 tries
		cube_t place_area;
		int room_ix(-1); // attic will be -1

		if (place_in_attic && has_attic() && (rgen.rand()%10) == 0) { // 10% of the time in the attic; counts as a dark room even if the light is on
			place_area = get_attic_part();
			place_area.z1() = place_area.z2();
		}
		else {
			room_ix = (rgen.rand() % interior->rooms.size());
			room_t const &room(interior->rooms[room_ix]);
			if (room.z1() > ground_floor_z1) continue; // not on the ground floor or basement
			if (has_parking_garage && room.z1() == ground_floor_z1 && (rgen.rand()&3))  continue; // prefer parking garage and backrooms 4x, since each is a single room
			if (pref_dark_room && n < 50 && is_room_lit(room_ix, (room.z1() + radius))) continue; // check for dark room in first 50 iterations
			if (has_pool() && room_ix == interior->pool.room_ix)                        continue; // don't place in a room with a pool
			place_area = room; // will represent the usable floor area
			place_area.z1() += get_fc_thickness(); // on top of the floor

			if (has_parking_garage && room.part_id == basement_part_ix) { // allow placement on multiple floors of parking garage
				unsigned const num_floors(calc_num_floors_room(room, floor_spacing, get_floor_thickness()));
				if (num_floors > 1) {place_area.z1() += floor_spacing*(rgen.rand()%num_floors);} // select a random floor
			}
		}
		place_area.expand_by_xy(-(radius + wall_thickness));
		if (min(place_area.dx(), place_area.dy()) < 4.0*radius) continue; // room too small (can happen for has_complex_floorplan office buildings)
		point const pos(gen_xy_pos_in_area(place_area, radius, rgen, place_area.z1()));
		if (!is_valid_ai_placement(pos, radius, 0)) continue; // check room objects; start in the open, not under something; skip_nocoll=0
		
		if (not_by_ext_door && room_ix >= 0 && pos.z > ground_floor_z1 && pos.z < (ground_floor_z1 + floor_spacing)) { // on ground floor
			float const min_door_dist(1.5*floor_spacing);
			room_t const &room(interior->rooms[room_ix]);
			bool near_door(0);
			
			for (tquad_with_ix_t const &door : doors) { // exterior doors
				if (door.type == tquad_with_ix_t::TYPE_RDOOR) continue; // ignore roof access doors
				cube_t door_bcube(door.get_bcube());
				door_bcube.expand_by_xy(wall_thickness); // make sure it overlaps the room
				if (!door_bcube.intersects(room)) continue;
				door_bcube.expand_by_xy(min_door_dist - wall_thickness);
				if (door_bcube.contains_pt_xy(pos)) {near_door = 1; break;}
			}
			if (near_door) continue; // bad placement; don't spawn venemous spiders and snakes near the door that the player could run into immediately
		}
		if (not_player_visible && n < 50 && fabs(camera_bs.z - pos.z) < floor_spacing && camera_pdu.sphere_visible_test((pos + xlate), radius)) {
			// may be visible to the player; checked for the first 50 iterations
			bool const same_room(get_room_containing_pt(pos) == get_room_containing_camera(camera_bs));
			if (!line_intersect_walls(pos, camera_bs, same_room)) continue; // line of sight, skip
		}
		return pos; // success
	} // for n
	return all_zeros; // failed
}

bool building_t::is_pos_inside_building(point const &pos, float xy_pad, float hheight, bool inc_attic) const {
	float bcube_pad(xy_pad);
	if (inc_attic && has_attic() && pos.z >= interior->attic_access.z2()) {bcube_pad += get_attic_beam_depth();} // add extra spacing for attic beams (approximate)
	
	if (has_basement() && pos.z < ground_floor_z1) { // in the basement
		cube_t const &basement(get_basement());

		if (has_ext_basement()) { // check union of basement + extended basement
			cube_t sc; sc.set_from_sphere(pos, bcube_pad); // sphere bounding cube
			float cont_area(0.0);
			// Note: not guaranteed to be correct since ext basement bcube can overlap basement and double count the area, but good enough for initial rejection
			accumulate_shared_xy_area(basement, sc, cont_area);
			accumulate_shared_xy_area(interior->basement_ext_bcube, sc, cont_area);
			if (cont_area < 0.99*sc.get_area_xy()) return 0; // not contained
			if (!get_full_basement_bcube().contains_pt_xy(pos)) return 0; // check union cube as an additional test in case there was double counting
		}
		else if (!basement.contains_pt_xy_exp(pos, -bcube_pad)) return 0; // check the basement
	}
	else if (!bcube.contains_pt_xy_exp(pos, -bcube_pad)) return 0; // check for end point inside building bcube
	cube_t req_area(pos, pos);
	req_area.expand_by_xy(xy_pad);
	req_area.z2() += hheight;
	if (!is_cube_contained_in_parts(req_area) && !(inc_attic && cube_in_attic(req_area))) return 0;
	if (!check_cube_within_part_sides(req_area)) return 0; // handle non-cube walls; probably slow; room isn't known
	return 1;
}

void update_dir_incremental_no_zero_check(vector3d &dir, vector3d const &new_dir, float turn_rate, float timestep) {
	float const delta_dir(turn_rate*min(1.0f, 1.5f*(1.0f - pow(0.7f, timestep))));
	dir = (delta_dir*new_dir + (1.0 - delta_dir)*dir).get_norm();
}
void update_dir_incremental(vector3d &dir, vector3d const &new_dir, float turn_rate, float timestep, rand_gen_t &rgen) {
	if (new_dir != zero_vector) {update_dir_incremental_no_zero_check(dir, new_dir, turn_rate, timestep);} // update dir if new_dir is valid
	if (dir == zero_vector) {dir = rgen.signed_rand_vector_spherical_xy_norm();} // dir must always be valid
}
void try_resolve_coll(vector3d const &dir, vector3d const &upv, vector3d const &coll_dir, vector3d &new_dir, bool try_tangent) {
	if (dot_product(coll_dir, new_dir) > 0.0) { // must move away from the collision direction
		if (try_tangent) { // try to preserve direction by moving in a tangent to the collider
			new_dir = cross_product(coll_dir, upv).get_norm(); // should be orthogonal, but we can be safe and normalize
			if (dot_product(new_dir, dir) < 0.0) {new_dir.negate();} // there are two solutions; choose the one closer to our current dir
		}
		else {new_dir.negate();} // otherwise reverse direction
	}
}

bool play_attack_sound(point const &pos, float gain, float pitch, rand_gen_t &rgen) {
	static double last_sound_time(0.0);

	if (tfticks - last_sound_time > 0.4*TICKS_PER_SECOND) {
		gen_sound_thread_safe(SOUND_SQUISH, pos, gain, pitch);
		last_sound_time = tfticks + 0.2f*TICKS_PER_SECOND*rgen.rand_float(); // add some randomness
		return 1;
	}
	return 0;
}


// *** Rats ***

rat_t::rat_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_, bool dead_, unsigned tix) :
	building_animal_t(pos_, radius_, dir_, id_), dest(pos), tunnel_tank_ix(tix), dead(dead_) {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT)); // L=3878, W=861, H=801
	hwidth = radius*sz.y/sz.x; // scale radius by ratio of width to length
	height = 2.0*radius*sz.z/max(sz.x, sz.y); // use max of x/y size; the x/y size represents the bcube across rotations
}
cube_t get_obj_model_bcube_for_dir(point const &pos, vector3d const &dir, float radius, float height, unsigned model_id) {
	bool const pri_dim(fabs(dir.x) < fabs(dir.y));
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id));
	cube_t bcube(pos, pos);
	bcube.expand_in_dim( pri_dim, radius); // larger dim
	bcube.expand_in_dim(!pri_dim, radius*min(sz.x, sz.y)/max(sz.x, sz.y)); // smaller dim
	bcube.z2() += height;
	return bcube;
}
cube_t rat_t::get_bcube_with_dir() const {
	return get_obj_model_bcube_for_dir(pos, dir, radius, height, OBJ_MODEL_RAT);
}
bool rat_t::is_facing_dest() const {
	return (dot_product((dest - pos).get_norm(), dir) > 0.75); // only move if we're facing our dest, to avoid walking through an object
}

// Note: rat_obj is used for dead/broken flag logic
bool building_t::add_rat(point const &pos, float hlength, vector3d const &dir, point const &placed_from, bool &dead) {
	if (!rat_t::allow_in_attic() && point_in_attic(pos)) return 0;
	point rat_pos(pos);
	if (!get_zval_of_floor(pos, hlength, rat_pos.z)) return 0; // place on the floor, skip if there's no floor here; doesn't work with glass floors
	if (point_in_water_area(rat_pos, 0)) return 0; // don't place rat in water; full_room_height=0
	rat_t rat(rat_pos, hlength, vector3d(dir.x, dir.y, 0.0).get_norm(), interior->room_geom->rats.size(), dead); // dir in XY plane
	
	if (check_line_coll_expand(pos, rat_pos, hlength, rat.height)) { // something is in the way
		point const test_pos(rat_pos + vector3d(0.0, 0.0, rat.height));
		auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) { // check for a toilet that we can drop the rat into
			if (i->room_id != get_room_containing_pt(placed_from)) continue; // wrong room

			if (i->type == TYPE_TOILET && i->contains_pt(test_pos)) {
				point const sound_origin(i->get_cube_center());
				gen_sound_thread_safe(SOUND_FLUSH, local_to_camera_space(sound_origin));
				register_building_sound(sound_origin, 0.5);
				register_achievement("Sleep with the Fishes");
				return 1;
			}
			bool killed(0);

			if (i->type == TYPE_MWAVE && i->is_powered() && check_line_clip(pos, rat_pos, i->d)) { // Note: not on the floor, so can't check test_pos
				gen_sound_thread_safe(SOUND_BEEP, local_to_camera_space(i->get_cube_center()), 0.25);
				killed = 1;
			}
			else if (i->type == TYPE_STOVE && i->item_flags != 0 && i->contains_pt(test_pos)) { // at least one burner is on (gas, not electric - always powered)
				killed = 1;
			}
			if (killed) { // cook/kill
				bool const new_achievement(register_achievement("Tastes Like Chicken"));
				register_fly_attract(new_achievement); // don't print a message if the achievement message was printed
				dead = 1;

				if (i->type == TYPE_MWAVE && i->is_open() && !i->is_nonempty()) { // open and empty, microwave, put the rat inside
					// see building_room_geom_t::add_mwave()
					bool const open_dir(i->dim ^ i->dir ^ 1);
					cube_t const panel(get_mwave_panel_bcube(*i));
					cube_t body(*i);
					body.d[!i->dim][!open_dir] = panel.d[!i->dim][open_dir]; // the other half
					float const tray_height(0.25*panel.get_sz_dim(!i->dim));
					i->flags |= RO_FLAG_NONEMPTY;
					rat.pos   = cube_bot_center(body) + tray_height*plus_z; // centered on the microwave tray
					rat.dead  = 1;
					interior->room_geom->rats.add(rat);
					interior->room_geom->modified_by_player = 1;
					return 1;
				}
				return 0; // rat is not dropped
			}
		} // for i
		return 0; // can't place the rat here
	}
	rat.fear_pos = placed_from;
	rat.fear     = 1.0; // starts off with max fear
	interior->room_geom->rats.add(rat);
	interior->room_geom->modified_by_player = 1;
	return 1;
}

void building_t::update_rats(point const &camera_bs, unsigned building_ix) {
	vect_rat_t &rats(interior->room_geom->rats);
	if (rats.placed && rats.empty()) return; // no rats placed in this building
	if (global_building_params.num_rats_max == 0) return;
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) return; // no rat model
	//timer_t timer("Update Rats"); // multi-part: 1.1ms, open door 1.7ms; office building 1.7ms
	add_animals_on_floor(rats, building_ix, global_building_params.num_rats_min, global_building_params.num_rats_max,
		global_building_params.rat_size_min, global_building_params.rat_size_max);
	// update rats; rats always attack when the player is dead or a phone is ringing
	unsigned const min_attack_rats((player_wait_respawn || phone_is_ringing()) ? 1 : global_building_params.min_attack_rats);
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	unsigned num_near_player(0);
	point rat_alert_pos;

	for (rat_t &rat : rats) { // must be done before sorting
		rat.move(timestep, rat.is_facing_dest());
		num_near_player += rat.near_player;
		if (num_near_player == min_attack_rats) {rat_alert_pos = rat.pos;}
	}
	static bool prev_can_attack_player(0);
	bool const can_attack_player(num_near_player >= min_attack_rats);
	
	if (can_attack_player && !prev_can_attack_player && !player_wait_respawn) { // play sound on first attack if player is alive
		gen_sound_thread_safe(SOUND_RAT_SQUEAK, local_to_camera_space(rat_alert_pos), 1.0, 1.2); // high pitch
	}
	prev_can_attack_player = can_attack_player;
	rats.do_sort();
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame

	if (frame_counter & 1) { // reverse iteration, to avoid directional bias
		for (auto r = rats.rbegin(); r != rats.rend(); ++r) {update_rat(*r, camera_bs, timestep, rats.max_xmove, can_attack_player, rgen);}
	}
	else { // forward iteration; ~0.004ms per rat
		for (auto r = rats. begin(); r != rats. end(); ++r) {update_rat(*r, camera_bs, timestep, rats.max_xmove, can_attack_player, rgen);}
	}
}

void building_t::update_sewer_rats(point const &camera_bs, unsigned building_ix) {
	vect_rat_t &rats(interior->room_geom->sewer_rats);
	if (rats.placed && rats.empty()) return; // no sewer rats placed in this building
	if (global_building_params.num_rats_max == 0 || interior->tunnels.empty()) return;
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) return; // no rat model
	float const floor_spacing(get_window_vspace()), timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, building_ix+123); // unique per building
	
	if (!rats.placed) { // add initial rats in tunnels
		for (auto i = interior->tunnels.begin(); i != interior->tunnels.end(); ++i) {
			tunnel_seg_t const &t(*i);
			if (t.room_conn)      continue; // don't add to room connector tunnel segment
			if (!t.has_gate)      continue; // only add a rat if it has a gate to hide behind
			if (rgen.rand_bool()) continue; // add 50% of the time
			float const radius(select_animal_radius(global_building_params.rat_size_min, global_building_params.rat_size_max, floor_spacing, rgen));
			if (radius < t.water_level) continue; // skip if underwater
			float v1(t.p[0][t.dim]), v2(t.p[1][t.dim]); // placement range
			if      (t.closed_ends[0]) {v1 = t.gate_pos;}
			else if (t.closed_ends[1]) {v2 = t.gate_pos;}
			v1 += 2.0*radius; v2 -= 2.0*radius; // add padding
			if (v1 >= v2) continue; // not enough space for a rat
			point pos;
			vector3d dir;
			pos[ t.dim] = rgen.rand_uniform(v1, v2);
			pos[!t.dim] = t.p[0][!t.dim]; // either point should work
			pos.z       = t.p[0].z - t.radius; // on the bottom of the tunnel
			dir[t.dim]  = ((pos[t.dim] < t.gate_pos) ? 1.0 : -1.0); // face toward the gate
			rats.add(rat_t(pos, radius, dir, 0, 0, (i - interior->tunnels.begin()))); // id=0, dead=0; store tunnel index
		} // for i
		rats.placed = 1;
	}
	for (rat_t &rat : rats) { // update logic
		if (rat.speed > 0.0) { // moving
			assert(rat.tunnel_tank_ix < interior->tunnels.size());
			tunnel_seg_t const &t(interior->tunnels[rat.tunnel_tank_ix]);
			bool const dir(rat.dir[t.dim] > 0.0); // movement direction along the tunnel
			float const stop_pos(0.5*(t.p[dir][t.dim] + t.gate_pos));
			
			if ((rat.pos[t.dim] < stop_pos) != dir) { // stop
				rat.speed     = 0.0;
				rat.is_hiding = 1; // safe behind the bars
			}
			else {rat.move(timestep);}
		}
		else if (player_in_tunnel && !rat.is_hiding) { // stopped and not safe
			if (dist_xy_less_than(rat.pos, camera_bs, 2.0*floor_spacing)) { // check for player proximity
				gen_sound_thread_safe(SOUND_RAT_SQUEAK, local_to_camera_space(rat.pos), 1.0, 1.1); // medium pitch
				rat.speed = global_building_params.rat_speed*rgen.rand_uniform(1.6, 2.0); // fast
			}
		}
	} // for rat
}

void building_t::update_pet_rats(point const &camera_bs, unsigned building_ix) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) return; // no rat model
	vect_rat_t &rats(interior->room_geom->pet_rats);
	vect_room_object_t const &objs(interior->room_geom->objs);

	if (!rats.placed && has_mall()) { // add pet store rats on first update frame
		for (pet_tank_t const &t : interior->mall_info->pet_tanks) {
			if (t.animal_type != TYPE_RAT) continue;
			assert(t.obj_ix < objs.size());
			room_object_t const &obj(objs[t.obj_ix]);
			if (obj.type != TYPE_FISHTANK) continue; // taken by the player?
			assert(obj.item_flags == TYPE_RAT);
			rand_gen_t rgen;
			rgen.set_state(building_ix+1, t.obj_ix+1); // unique per building and per tank
			rgen.rand_mix();
			unsigned const num_rats((rgen.rand() % 3) + 2); // 2-4
			float const height(obj.dz()), zval(obj.z1() + 0.1*height); // around substrate height

			for (unsigned n = 0; n < num_rats; ++n) {
				float const radius(rgen.rand_uniform(0.5, 1.0)*0.15*height);
				point const pos(gen_xy_pos_in_area(obj, vector3d(radius, radius, radius), rgen, zval));
				rats.emplace_back(pos, radius, rgen.signed_rand_vector_spherical_xy_norm(), rats.size(), 0, t.obj_ix); // dead=0
				rats.back().dest = rats.back().dir;
			}
		} // for t
		rats.placed = 1;
	}
	if (rats.empty()) return; // no pet rats placed in this building
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	bool any_removed(0);
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame
	auto tank_start(rats.begin());

	for (auto i = rats.begin(); i != rats.end(); ++i) { // update logic
		rat_t &rat(*i);
		if (rat.tunnel_tank_ix != tank_start->tunnel_tank_ix) {tank_start = i;} // start a new tank
		// check for tank removed
		assert(rat.tunnel_tank_ix < objs.size());
		room_object_t const &obj(objs[rat.tunnel_tank_ix]);

		if (obj.type != TYPE_FISHTANK) { // taken by the player?
			any_removed = rat.dead = 1; // will be removed below
			continue;
		}
		if (rat.is_sleeping()) {
			if ((float)tfticks > rat.wake_time) {rat.wake_time = 0.0;} // time to wake up
		}
		else if (rat.dist_since_sleep > 0.75*(obj.dx() + obj.dy())) { // maybe stop and rest
			rat.sleep_for(1.0, 10.0); // 2-10s
			rat.speed = 0.0; // will reset anim_time in the next frame
		}
		else if (rat.speed > 0.0) { // moving
			rat.move(timestep, 1, 0.5); // can_move_forward=1, anim_time_scale=0.5
			// check for invalid pos; should this be predicted with lookahead?
			point const old_pos(rat.pos);
			cube_t valid_area(obj);
			valid_area.expand_by_xy(-(0.02*obj.dz() + rat.radius));
			valid_area.clamp_pt_xy(rat.pos);
			vector3d coll_normal;

			if (rat.dest != rat.dir) { // turning
				update_dir_incremental(rat.dir, rat.dest, 0.5, timestep, rgen); // turn_rate=0.5
				if (dot_product(rat.dir, rat.dest) > 0.99) {rat.dest = rat.dir;} // turn complete
			}
			else if (rat.pos != old_pos) { // moved due to coll
				coll_normal = (old_pos - rat.pos);
			}
			else { // check for collisions with other rats
				for (auto j = tank_start; j != rats.end() && j->tunnel_tank_ix == i->tunnel_tank_ix; ++j) {
					if (j == i || !dist_xy_less_than(i->pos, j->pos, (i->radius + j->radius))) continue;
					coll_normal = (j->pos - rat.pos);
					break; // only handles one coll
				}
			}
			if (coll_normal != zero_vector && dot_product(rat.dir, coll_normal) > 0.0) {
				// gradually turn; dest is the destination dir here, not the destination pos
				rat.dest = rgen.signed_rand_vector_spherical_xy_norm();
				if (dot_product(rat.dest, coll_normal) > 0.0) {rat.dest.negate();} // must point away from the collision

				if (rgen.rand_float() < 0.25) { // maybe stop and rest
					rat.sleep_for(1.0, 5.0); // 1-5s
					rat.speed = 0.0; // will reset anim_time in the next frame
				}
			}
		}
		else { // stopped
			rat.speed = global_building_params.rat_speed*rgen.rand_uniform(0.1, 0.15); // slower than free rats
		}
	} // for rat i
	if (any_removed) {rats.erase(remove_if(rats.begin(), rats.end(), [](rat_t const &r) {return r.dead;}), rats.end());}
}

bool can_hide_under(room_object_t const &c, cube_t &hide_area) {
	cube_t dishwasher; // used below

	if (c.is_open() && c.is_small_closet()) { // open small closet
		hide_area = c;
		hide_area.expand_by(-get_closet_wall_thickness(c)); // we want the inside of the closet, excluding the walls
		hide_area.z1() += 0.5*hide_area.dz(); // use the halfway point; somewhat arbitrary, but will affect the score
		return 1;
	}
	else if (c.type == TYPE_BED) {
		cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
		get_bed_cubes(c, cubes);
		hide_area = cubes[0]; // frame
		return 1;
	}
	// skip fallen over furniture and tables with a central support (round and mall)
	else if ((c.type == TYPE_DESK || (c.type == TYPE_TABLE && c.shape != SHAPE_CYLIN && c.item_flags == 0)) && !c.is_on_floor()) {
		cube_t cubes[5];
		get_table_cubes(c, cubes); // body and legs
		hide_area = cubes[0]; // body
		return 1;
	}
	else if (c.type == TYPE_CONF_TABLE) { // can't currently hide under
		//cube_t cubes[2]; // {top, base}
		//get_conf_table_cubes(c, cubes);
	}
	else if (c.type == TYPE_DRESSER || c.type == TYPE_NIGHTSTAND) {
		hide_area = get_dresser_middle(c);
		return 1;
	}
	else if (c.type == TYPE_CHAIR && !c.is_on_floor()) { // skip fallen over chairs
		cube_t cubes[3]; // seat, back, legs_bcube
		get_chair_cubes(c, cubes);
		hide_area = cubes[0]; // seat
		return 1;
	}
	else if (c.type == TYPE_BCASE && !c.is_on_floor()) { // skip fallen over bookcases
		cube_t top, middle, back, lr[2];
		get_bookcase_cubes(c, top, middle, back, lr);
		hide_area = middle;
		return 1;
	}
	else if (c.type == TYPE_BENCH) {
		cube_t cubes[4]; // seat, lo side, hi side, [back]
		get_bench_cubes(c, cubes); // Note: return value ignored
		hide_area = cubes[0]; // seat
		return 1;
	}
	else if (c.type == TYPE_KSINK && get_dishwasher_for_ksink(c, dishwasher)) {
		hide_area = dishwasher;
		hide_area.d[c.dim][!c.dir] = c.d[c.dim][!c.dir]; // use the back of the cabinet, not the back of the dishwasher door
		return 1;
	}
	else if (c.type == TYPE_COUCH) {
		hide_area = c;
		hide_area.z1() += 0.06*c.dz(); // there's space under the couch
		return 1;
	}
	else if (c.type == TYPE_POOL_TABLE) {
		hide_area = c;
		hide_area.z1() += 0.5*c.dz(); // there's space under the pool table; could call get_pool_table_cubes(), but this should be good enough
		return 1;
	}
	else if (c.type == TYPE_RCHAIR) {
		return 0; // not a hiding spot - yet
	}
	else if (c.type == TYPE_BRSINK) { // office building bathroom sink
		// not a very good hiding spot, but there aren't many in office buildings;
		// this is just a placeholder anyway, since sinks don't extend down to the floor and won't pass the rat zval test
		hide_area = c;
		return 1;
	}
	else if (c.type == TYPE_SHELVES) {
		if (c.get_num_shelves() < 4) return 0; // only 4 shelves provide a small enough space at the bottom to hide under
		hide_area = c;
		hide_area.z1() += 0.05*c.dz();
		return 1;
	}
	else if (c.type == TYPE_SHELFRACK) {
		cube_t back, top, sides[2], shelves[5];
		unsigned const num_shelves(get_shelf_rack_cubes(c, back, top, sides, shelves));

		if (num_shelves > 0) { // can possibly hide under the bottom shelf, though there's not much space
			hide_area = shelves[0]; // ignores the back
			set_cube_zvals(hide_area, c.z1(), shelves[0].z1());
			return 1;
		}
	}
	else if (c.is_parked_car()) { // parked car
		hide_area = c;
		hide_area.z1() += 0.12*c.dz(); // there's space under the car
		return 1;
	}
	//else if (c.type == TYPE_STAIR) {} // what about hiding under the stairs? only for SHAPE_WALLED_SIDES?
	return 0;
}

template<bool xy_only> class dir_gen_t {
	vector<vector3d> dirs;
	unsigned dir_ix;

	void gen_dirs() {
		rand_gen_t rgen;
		dirs.resize(1009); // make it a prime number
		for (auto &dir : dirs) {dir = (xy_only ? rgen.signed_rand_vector_spherical_xy() : rgen.signed_rand_vector_spherical()).get_norm();}
	}
public:
	dir_gen_t() : dir_ix(0) {}

	vector3d const &gen_dir() {
		if (dirs.empty()) {gen_dirs();}
		vector3d const &dir(dirs[dir_ix++]);
		if (dir_ix == dirs.size()) {dir_ix = 0;} // wrap around
		return dir;
	}
};
dir_gen_t<1> dir_gen_xy;
dir_gen_t<0> dir_gen_xyz;

void building_t::rat_bite_player(point const &pos, float damage, rand_gen_t &rgen) {
	play_attack_sound(local_to_camera_space(pos), 1.0, 1.0, rgen);
	bool const player_dead(player_take_damage(damage));
	if (player_dead) {register_achievement("Rat Food");} // damage over time; achievement if the player dies
	if (player_dead) {register_player_death(cur_player_building_loc.pos);}
}

// Note: non-const because biting the player can add blood and the player's inventory items to the building
void building_t::update_rat(rat_t &rat, point const &camera_bs, float timestep, float &max_xmove, bool can_attack_player, rand_gen_t &rgen) {

	if (rat.dead) return; // nothing to do
	float const floor_spacing(get_window_vspace()), trim_thickness(get_trim_thickness()), view_dist(RAT_VIEW_FLOORS*floor_spacing);
	float const hlength(rat.get_hlength()), hwidth(rat.hwidth), height(rat.height), hheight(0.5*height);
	float const squish_hheight(0.75*hheight); // rats can squish to get under low objects and walk onto small steps
	float const coll_radius(1.2f*hwidth); // slightly larger than half-width; maybe should use length so that the rat doesn't collide when turning?
	float const line_project_dist(max(1.1f*(hlength - coll_radius), 0.0f)); // extra space in front of the target destination
	// set dist_thresh based on the distance we can move this frame; if set too low, we may spin in circles trying to turn to stop on the right spot
	float const dist_thresh(2.0f*timestep*max(rat.speed, global_building_params.rat_speed));
	float const xy_pad(hlength + trim_thickness);
	vector3d const center_dz(0.0, 0.0, hheight); // or squish_hheight?
	bool update_path(0);
	vector3d coll_dir;
	point const prev_pos(rat.pos); // capture the pre-collision point
	static bool wrong_orient_warning(0);

	if (hwidth > hlength && !wrong_orient_warning) {
		cout << "*** Warning: Rat model is likely in the wrong orient: " << TXT(hlength) << TXT(hwidth) << endl;
		wrong_orient_warning = 1;
	}
	rgen.rand_mix(); // make sure it's different per rat

	if (rat.is_sleeping() && rat.fear == 0.0) {} // peacefully sleeping, no collision needed
	else if (check_and_handle_dynamic_obj_coll(rat.pos, rat.pos, rat.radius, rat.pos.z, (rat.pos.z + height), camera_bs, 0)) { // check for collisions; for_spider=0
		coll_dir = (prev_pos - rat.pos).get_norm(); // points toward the collider in the XY plane

		// check if new pos is valid, and has a path to dest
		if (!is_pos_inside_building(rat.pos, xy_pad, hheight)) {
			rat.pos = prev_pos; // restore previous pos before collision
			rat.sleep_for(0.1, 0.2); // wait 0.1-0.2s so that we don't immediately collide and get pushed out again
		}
		else if (check_line_coll_expand((rat.pos + center_dz), (rat.dest + center_dz), coll_radius, squish_hheight)) {
			rat.pos = prev_pos; // restore previous pos before collision
		}
		else {
			max_eq(max_xmove, fabs(rat.pos.x - prev_pos.x));
			// update the path about every 30 frames of colliding; this prevents the rat from being stuck while also avoiding jittering due to frequent dest updates;
			// using randomness rather than updating every 30 frames makes this process less regular and mechanical
			update_path = ((rgen.rand()%30) == 0);
		}
	}
	point const p1(rat.pos + center_dz);
	vector3d dir_to_fear;
	bool has_fear_dest(0);
	// apply scare logic
	bool const was_scared(rat.fear > 0.0);
	scare_rat(rat, camera_bs);
	rat.attacking = (rat.near_player && can_attack_player);
	if (rat.attacking) {rat.fear = 0.0;} // no fear when attacking
	bool const is_scared(rat.fear > 0.0), newly_scared(is_scared && !was_scared);

	// determine destination
	if (rat.attacking) {
		float const player_radius(get_scaled_player_radius()), min_dist(player_radius + hlength);
		point target(camera_bs.x, camera_bs.y, rat.pos.z);
		vector3d const vdir((target - rat.pos).get_norm());
		target -= vdir*(1.01*min_dist); // get within attacking range, but not at the center of the player
		
		if (is_pos_inside_building(target, xy_pad, hheight)) { // check if outside the valid area
			point const p1_ext(p1 + coll_radius*vdir); // move the line slightly toward the dest to prevent collisions at the initial location
			
			if (!check_line_coll_expand(p1_ext, target, coll_radius, squish_hheight)) {
				rat.dest      = target;
				rat.speed     = RAT_ATTACK_SPEED*global_building_params.rat_speed;
				rat.wake_time = 0.0; // wake up
				update_path   = 0;
				if (dist_xy_less_than(rat.pos, target, 0.05*min_dist)) {rat_bite_player(rat.pos, 0.004, rgen);} // do damage when nearly colliding with the player
			}
		}
	}
	if (is_scared) {
		// find hiding spot (pref in opposite direction from fear_pos);
		// we must check this each frame in case the player took or moved the object we were hiding under
		float const rat_z1(rat.pos.z), rat_z2(rat.pos.z + height), rat_squish_z2(p1.z + squish_hheight);
		point best_dest;
		float best_score(0.0);
		dir_to_fear   = (rat.fear_pos - rat.pos);
		dir_to_fear.z = 0.0; // XY plane only
		dir_to_fear.normalize();
		rat.wake_time = 0.0; // wake up
		vect_room_object_t::const_iterator b, e;
		get_begin_end_room_objs_on_ground_floor(rat_z2, 0, b, e); // for_spider=0

		for (auto c = b; c != e; ++c) {
			if (c->z1() > rat_z2 || c->z2() < rat_z1) continue; // wrong floor, or object not on the floor
			cube_t hide_area; // will be a subset of the object
			if (c->shape != SHAPE_CUBE || !can_hide_under(*c, hide_area)) continue; // only cubes for now
			float const top_gap(hide_area.z1() - rat_squish_z2); // space between top of rat and bottom of object
			if (top_gap < 0.0) continue; // rat can't fit under this object
			if (!dist_xy_less_than(hide_area.get_cube_center(), p1, view_dist)) continue; // too far away to see
			// select our destination under this hiding spot;
			// this must be unique per rat so that rats don't compete for the exact same spot, and must be the same across calls for stability;
			// also use the obj_id to mix things up between objects, and throw in the type in case obj_id is left at 0;
			// we can't use the obj vector position because it may change if the player takes or drops objects
			rand_gen_t my_rgen;
			my_rgen.set_state((rat.id + 1), (c->obj_id + (c->type << 16) + 1));
			cube_t safe_area(hide_area);
			point cand_dest(0.0, 0.0, p1.z); // x/y will be set below

			for (unsigned d = 0; d < 2; ++d) {
				// shrink by half length so that any point inside this cube will be covered; add an extra factor of 1.5 to avoid table/chair/desk legs, etc.;
				// make sure safe_area is strictly normalized
				safe_area.expand_in_dim(d, -min(1.5f*hlength, 0.49f*safe_area.get_sz_dim(d)));
				cand_dest[d] = my_rgen.rand_uniform(safe_area.d[d][0], safe_area.d[d][1]);
			}
			float const dist(p2p_dist(p1, cand_dest));
			
			if (dist < dist_thresh) { // already at this location
				if (check_line_coll_expand(p1, cand_dest, coll_radius, squish_hheight)) {update_path = 1; continue;} // location is invalid, need to update the path below
				has_fear_dest = 1; // it's valid, so stay there
				rat.speed     = 0.0;
				break;
			}
			float side_cov(0.5f*min(hide_area.dx(), hide_area.dy()) - hlength); // amount of overhang of the object around the rat's extents
			float const dist_to_fear(p2p_dist(rat.fear_pos, cand_dest));
			float score(side_cov - 0.5f*top_gap + 0.2f*dist_to_fear - 0.1f*max(dist, dist_thresh)); // can be positive or negative
			if (best_score != 0.0 && score <= best_score) continue; // check score before iterating over other rats; it can only decrease below
			if (check_line_coll_expand(p1, cand_dest, coll_radius, squish_hheight)) continue; // use center before checking other rats so that the entire path is valid
			float tot_mdist(0.0);
			bool skip(0);
			float const radius_scale = 0.8; // smaller dist (head can overlap tail)
			vect_rat_t const &rats(interior->room_geom->rats);
			float const rsum_max(radius_scale*(rat.radius + rats.max_radius) + rats.max_xmove), coll_x1(cand_dest.x - rsum_max), coll_x2(cand_dest.x + rsum_max);
			auto start(rats.get_first_with_xv_gt(coll_x1)); // use a binary search to speed up iteration

			for (auto r = start; r != rats.end(); ++r) {
				if (r->pos.x > coll_x2) break; // no rat after this can overlap - done
				if (&(*r) == &rat) continue; // skip ourself
				float const r_sum(radius_scale*(rat.radius + r->radius)); // smaller dist (head can overlap tail)
				if (!dist_xy_less_than(cand_dest, r->pos, r_sum)) continue; // no rat in this spot
				float const move_dist(1.01*r_sum - p2p_dist_xy(cand_dest, r->pos)); // slightly larger than r_sum to prevent collisions
				cand_dest += (p1 - cand_dest).get_norm()*move_dist; // move our target in front of this other rat
				side_cov  -= move_dist; // moving to this misaligned position loses side coverage
				score      = 4.0*side_cov - 0.5f*top_gap + 0.25f*dist_to_fear - 0.1*max(dist, dist_thresh); // update score
				score     -= 0.2*dist; // less desirable when occupied
				score     -= 2.0*move_dist; // use prev accumulated move dist; even less desirable if there are many rats in the way
				tot_mdist += move_dist;
				if (tot_mdist > 4.0*rat.radius || !hide_area.contains_pt_xy(cand_dest)) {skip = 1; break;} // moved too far, too many other rats at this dest, skip it
				if (best_score != 0.0 && score <= best_score) {skip = 1; break;} // score dropped too low
				r = start-1; // go back and test the other rats for collisions with this new position; maybe creates a bit more determinism/less chaos
			} // for r
			if (skip) continue;
			if (tot_mdist > 0.0 && !is_pos_inside_building(cand_dest, xy_pad, hheight)) continue; // check if outside the valid area if center was moved
			best_dest  = point(cand_dest.x, cand_dest.y, rat.pos.z); // keep zval on the floor
			best_score = score;
			if (cand_dest.x == rat.dest.x && cand_dest.y == rat.dest.y) break; // keep the same dest (optimization)
		} // for c
		if (!has_fear_dest && best_score != 0.0) { // found a valid hiding place; score can be positive or negative
			rat.dest = best_dest;
			if (dist_less_than(rat.pos, rat.dest, dist_thresh)) {rat.speed = 0.0;} // close enough - stop
			else {rat.speed = RAT_FEAR_SPEED*global_building_params.rat_speed;} // high speed if not at dest
			has_fear_dest = 1; // set to avoid triggering the code below if close to dest
			assert(rat.pos.z == rat.dest.z);
		}
		rat.fear = max(0.0f, (rat.fear - 0.2f*(timestep/TICKS_PER_SECOND))); // reduce fear over 5s
	}
	bool const is_at_dest(dist_less_than(rat.pos, rat.dest, dist_thresh));

	if (!is_scared && !rat.is_sleeping() && is_at_dest && rat.dist_since_sleep > 1.5*floor_spacing && rgen.rand_bool()) { // 50% chance of taking a rest
		rat.sleep_for(0.0, 4.0); // 0-4s
		rat.speed = 0.0; // will reset anim_time in the next frame
	}
	else if (!has_fear_dest && !rat.is_sleeping() &&
		(rat.speed == 0.0 || newly_scared || update_path || is_at_dest || check_line_coll_expand(rat.pos, rat.dest, coll_radius, hheight)))
	{
		// stopped, no dest, at dest, collided, or newly scared - choose a new dest
		float target_fov_dp(RAT_FOV_DP), target_max_dist(view_dist); // start at nominal/max values
		float const dist_upper_bound(0.12 + 0.88*(1.0 - rat.fear)); // shorten the distance based on the amount of fear to evade more easily
		float const min_step(min(dist_thresh, 0.05f*rat.radius));
		rat.speed = 0.0; // stop until we've found a valid destination

		for (unsigned n = 0; n < 200; ++n) { // make 200 tries
			if (n > 50) { // we've been at this for a while, maybe we need to relax our constraints, maybe follow the walls?
				target_fov_dp   -= 0.02; // allow for turns outside our field of view
				target_max_dist *= 0.96;  // decrease the max distance considered
			}
			vector3d vdir(dir_gen_xy.gen_dir()); // random XY direction

			if (coll_dir != zero_vector) { // resolve the collision; target_fov_dp is ignored in this case
				try_resolve_coll(rat.dir, plus_z, coll_dir, vdir, (n <= 10)); // if earlier in the iteration, try moving in a tangent
			}
			else { // not colliding; check if the new direction is close enough to our current direction
				float dp(dot_product(rat.dir, vdir));
				if (n < 180 && dp < 0.0) {vdir.negate(); dp = -dp;} // only allow switching directions in the last 20 iterations; can still fail the test below
				if (dp < target_fov_dp) continue; // not in field of view, use a new direction
			}
			if (is_scared && n <= 100 && dot_product(dir_to_fear, vdir) > 0.0) continue; // don't move toward danger; may make the rat move into a corner
			float dist(rgen.rand_uniform(0.1, dist_upper_bound)*target_max_dist); // random distance out to max view dist
			max_eq(dist, min_step); // make sure distance isn't too short
			point const cand(rat.pos + dist*vdir);
			if (!is_pos_inside_building(cand, xy_pad, hheight)) continue; // check if outside the valid area
			point const p2(cand + line_project_dist*vdir + center_dz); // extend in vdir so that the head doesn't collide
			point const p1_ext(p1 + coll_radius*vdir); // move the line slightly toward the dest to prevent collisions at the initial location
			if (check_line_coll_expand(p1_ext, p2, coll_radius, squish_hheight)) continue;
			rat.dest  = cand;
			rat.speed = global_building_params.rat_speed*rgen.rand_uniform(0.5, 1.0)*(is_scared ? 1.5 : 1.0); // random speed
			break; // success
		} // for n
		assert(rat.pos.z == rat.dest.z);
	}
	// update direction
	vector3d new_dir;

	if (!dist_less_than(rat.pos, rat.dest, dist_thresh)) {
		new_dir = (rat.dest - rat.pos).get_norm(); // point toward our destination
	}
	else if (has_fear_dest) { // stop, rest, and point toward what we fear
		max_eq(max_xmove, fabs(rat.pos.x - rat.dest.x));
		new_dir   = dir_to_fear;
		rat.speed = rat.dist_since_sleep = 0.0;
		rat.pos   = rat.dest; // just move it to the dest to prevent instability
	}
	// else dir is unchanged
	rat.is_hiding = (has_fear_dest && dist_less_than(rat.pos, rat.dest, 2.0*dist_thresh)); // close to fear_dest
	float const turn_rate(is_scared ? 1.1 : 1.0); // higher turning rate when scared
	update_dir_incremental(rat.dir, new_dir, turn_rate, timestep, rgen);
	assert(rat.dir.z == 0.0); // must be in XY plane
}

void building_t::scare_rat(rat_t &rat, point const &camera_bs) const {
	// Note: later calls to scare_rat_at_pos() have priority and will set rat.fear_pos, but all calls will accumulate fear
	assert(interior);
	point fire_pos;

	if (has_room_geom() && interior->room_geom->fire_manager.get_closest_fire(rat.pos, 2.0*rat.radius, rat.pos.z, rat.pos.z+rat.height, &fire_pos)) { // fire
		//scare_rat_at_pos(rat, fire_pos, 1.0, 1);
		rat.fear     = 1.0; // max fear
		rat.fear_pos = fire_pos;
		return; // done
	}
	float const sight_scare_amt = 0.5;
	rat.near_player = 0;

	for (person_t const &p : interior->people) { // other people in the building scare the rats
		if (p.is_waiting_or_stopped()) continue; // only scare if moving
		scare_rat_at_pos(rat, point(p.pos.x, p.pos.y, p.get_z1()), sight_scare_amt, 1);
	}
	if (camera_surf_collide && camera_in_building) {
		if (global_building_params.min_attack_rats > 0 && in_building_gameplay_mode()) { // rat attacks are enabled in gameplay mode
			// determine if the player is close and visible for attack strength; can't use a return value of scare_rat_at_pos() due to early termination
			if (fabs(rat.pos.z - camera_bs.z) < get_window_vspace()) { // same floor
				if (dist_less_than(rat.pos, camera_bs, RAT_VIEW_FLOORS*get_window_vspace())) { // close enough; doesn't have to be in the same room
					if (check_line_of_sight_large_objs(rat.pos, camera_bs)) {rat.near_player = 1;} // what about interior windows?
				}
			}
		}
		scare_rat_at_pos(rat, camera_bs, sight_scare_amt, 1); // the sight of the player walking in the building scares the rats
	}
	sphere_t const cur_sound(get_cur_frame_loudest_sound());
	if (cur_sound.radius > 0.0) {scare_rat_at_pos(rat, cur_sound.pos, 4.0*cur_sound.radius, 0);}
}

void building_t::scare_rat_at_pos(rat_t &rat, point const &scare_pos, float amount, bool by_sight) const {
	assert(amount > 0.0);
	if (fabs(rat.pos.z - scare_pos.z) > get_window_vspace()) return; // on a different floor, ignore
	if (rat.fear > 0.99 && dist_less_than(rat.fear_pos, scare_pos, rat.radius)) return; // already max fearful of this location (optimization)
	point const pos(rat.get_center()); // use center zval, not floor zval
	int const scare_room(get_room_containing_pt(scare_pos)), rat_room(get_room_containing_pt(pos));
	//assert(rat_room >= 0); // this generally doesn't fail, unless rats collide with each other and push each other outside of a room
	if (rat_room != scare_room) {amount *= 0.67;} // less fearful if in a different room
	float const max_scare_dist(RAT_VIEW_FLOORS*get_window_vspace()), scare_dist(max_scare_dist*min(amount, 1.0f));
	float const fear((scare_dist - p2p_dist(pos, scare_pos))/max_scare_dist);
	if (fear <= 0.0) return;
	if (by_sight && !check_line_of_sight_large_objs(pos, scare_pos)) return; // check line of sight
	rat.fear     = min(1.0f, (rat.fear + fear));
	rat.fear_pos = scare_pos;
}


// *** Spiders ***

// Note: radius is really used for height and body size; legs can extend out to ~2x radius, so we decrease radius to half the user-specified value to account for this
spider_t::spider_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_, unsigned tix, bool in_tank_) :
	building_animal_t(pos_, 0.5*radius_, dir_, id_), upv(plus_z), tunnel_tank_ix(tix), in_tank(in_tank_)
{
	pos.z += radius; // shift upward so that the center is off the ground
}
cube_t spider_t::get_bcube() const {
	cube_t bcube(pos);
	bcube.expand_by(vector3d(2.0*radius, 2.0*radius, radius)); // conservative?
	return bcube;
}
vector3d spider_t::get_size() const {
	float const xy_radius(get_xy_radius());
	vector3d sz(xy_radius, xy_radius, xy_radius);
	sz[get_max_dim(upv)] = radius; // use height for primary up dir, which should at least work for walking on cube surfaces
	return sz;
}
void spider_t::choose_new_dir(rand_gen_t &rgen) {
	if (upv == zero_vector) {upv = plus_z;} // can this ever happen?
	dir = cross_product(rgen.signed_rand_vector_spherical(), upv).get_norm(); // must be orthogonal to upv
	if (SPIDERS_FOLLOW_PLAYER && dot_product_ptv(dir, pos, get_camera_building_space()) > 0.0) {dir = -dir;} // invert dir if moving away from the player; approximate
}
void spider_t::jump(float vel) {
	assert(vel > 0.0);
	jump_vel_z = vel; // overwrite if already set
	speed      = 0.5*vel; // set forward speed lower for a steeper jump; dir must be set (and should probably be in the XY plane)
}
void spider_t::end_jump() {
	if (is_jumping()) {jump_vel_z = speed = 0.0;}
}

bool building_room_geom_t::maybe_spawn_spider_in_drawer(room_object_t const &c, cube_t const &drawer, unsigned drawer_id, float floor_spacing, bool is_door) {
	if (global_building_params.spider_drawer_prob == 0.0) return 0; // no spiders
	if (!spider_t::allow_in_attic() && c.in_attic())      return 0; // no spiders in the attic
	rand_gen_t rgen;
	rgen.set_state((unsigned(c.obj_id) << 8)+drawer_id+1, c.room_id+1);
	if (!rgen.rand_probability(global_building_params.spider_drawer_prob)) return 0; // no spider
	float radius(0.5f*floor_spacing*rgen.rand_uniform(global_building_params.spider_size_min, global_building_params.spider_size_max));
	min_eq(radius, rgen.rand_uniform(0.6, 0.9)*min(drawer.dz(), 0.5f*min(drawer.dx(), drawer.dy()))); // make sure it fits in the drawer
	vector3d dir;
	dir[ c.dim] = (c.dir ? 1.0 : -1.0); // face the outside of the drawer
	dir[!c.dim] = 0.25*rgen.signed_rand_float(); // not straight out
	dir.normalize();
	spiders.emplace_back(cube_bot_center(drawer), radius, dir, spiders.size());
	if (!is_door) {spiders.back().jump(0.002*rgen.rand_uniform(1.0, 1.4));} // jump out of drawers
	return 1;
}

void building_t::update_spiders(point const &camera_bs, unsigned building_ix) {
	vect_spider_t &spiders(interior->room_geom->spiders);
	bool const was_placed(spiders.placed);
	if (was_placed && spiders.empty()) return; // no spiders placed in this building
	//timer_t timer("Update Spiders");
	add_animals_on_floor(spiders, building_ix, global_building_params.num_spiders_min, global_building_params.num_spiders_max,
		global_building_params.spider_size_min, global_building_params.spider_size_max);

	if (!was_placed && has_mall()) { // add pet store spiders on first update frame
		for (pet_tank_t const &t : interior->mall_info->pet_tanks) {
			if (t.animal_type != TYPE_SPIDER) continue;
			vect_room_object_t const &objs(interior->room_geom->objs);
			assert(t.obj_ix < objs.size());
			room_object_t const &obj(objs[t.obj_ix]);
			if (obj.type != TYPE_FISHTANK) continue; // taken by the player?
			assert(obj.item_flags == TYPE_SPIDER);
			rand_gen_t rgen;
			rgen.set_state(building_ix+1, t.obj_ix+1); // unique per building and per tank
			rgen.rand_mix();
			bool const one_big_spider(rgen.rand_bool());
			unsigned const num_spiders((one_big_spider ? 0 : (rgen.rand() % 4)) + 1); // 1-5
			float const height(obj.dz()), zval(obj.z1() + 0.1*height); // around substrate height

			for (unsigned n = 0; n < num_spiders; ++n) {
				float const radius(rgen.rand_uniform(0.5, 1.0)*(one_big_spider ? 0.2 : 0.1)*height);
				point const pos(gen_xy_pos_in_area(obj, vector3d(radius, radius, radius), rgen, zval));
				spiders.emplace_back(pos, radius, rgen.signed_rand_vector_spherical_xy_norm(), spiders.size(), t.obj_ix, 1); // in_tank=1
			}
		} // for t
	}
	// update spiders
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	
	for (spider_t &spider : spiders) {
		if (!spider.squished) {spider.move(timestep);}
	}
	spiders.do_sort();
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame
	for (spider_t &spider : spiders) {update_spider(spider, camera_bs, timestep, spiders.max_xmove, rgen);}
}

void set_spider_speed(spider_t &spider, rand_gen_t &rgen, float floor_spacing) {
	spider.speed = global_building_params.spider_speed*rgen.rand_uniform(0.5, 1.0);
	if (spider.in_tank) {spider.speed *= min(1.0f, 2.0f*spider.radius/(floor_spacing*global_building_params.spider_size_min));} // scale pet store spider speed with radius
}
void building_t::update_sewer_spiders(point const &camera_bs, unsigned building_ix) {
	vect_spider_t &spiders(interior->room_geom->sewer_spiders);
	if (spiders.placed && spiders.empty()) return; // no spiders placed in this building
	if (global_building_params.num_spiders_max == 0 || interior->tunnels.empty()) return;
	float const floor_spacing(get_window_vspace()), timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, building_ix+123); // unique per building

	if (!spiders.placed) {
		for (auto i = interior->tunnels.begin(); i != interior->tunnels.end(); ++i) {
			tunnel_seg_t const &t(*i);
			if (t.room_conn)      continue; // don't add to room connector tunnel segment
			if (rgen.rand_bool()) continue; // add 50% of the time
			float const radius(select_animal_radius(global_building_params.spider_size_min, global_building_params.spider_size_max, floor_spacing, rgen));
			float const v1(t.p[0][t.dim] + 4.0*radius), v2(t.p[1][t.dim] - 4.0*radius); // placement range
			if (v1 >= v2) continue; // not enough space for a spider
			float const val(rgen.rand_uniform(v1, v2));
			// check to avoid placing spiders over vertical shafts
			bool bad_place(0);

			for (tunnel_conn_t const &c : t.conns) {
				if (c.dim == 2 && fabs(val - c.pos) < (c.radius + 2.0*radius)) {bad_place = 1; break;}
			}
			if (bad_place) continue; // skip this spider; could generate a new val
			point pos;
			vector3d dir;
			pos[ t.dim] = val;
			pos[!t.dim] = t.p[0][!t.dim]; // either point should work
			pos.z       = t.p[0].z + t.radius - radius; // on the top of the tunnel
			dir[t.dim]  = (rgen.rand_bool() ? 1.0 : -1.0); // face a random direction
			spider_t spider(pos, radius, dir, 0, (i - interior->tunnels.begin())); // id=0, store tunnel index
			spider.upv  = -plus_z; // upside down
			set_spider_speed(spider, rgen, floor_spacing);
			spiders.add(spider);
		} // for i
		spiders.placed = 1;
	}
	if (player_in_tunnel) {
		for (spider_t &spider : spiders) { // update logic
			assert(spider.tunnel_tank_ix < interior->tunnels.size());
			tunnel_seg_t const &t(interior->tunnels[spider.tunnel_tank_ix]);
			float &val(spider.pos[t.dim]);
			float v1(t.p[0][t.dim] + 2.0*spider.radius), v2(t.p[1][t.dim] - 2.0*spider.radius); // movement range

			// avoid vertical shafts
			for (tunnel_conn_t const &c : t.conns) {
				if (c.dim < 2) continue; // not vertical
				if (val < c.pos) {min_eq(v2, (c.pos - c.radius - 2.0f*spider.radius));} // spider to the left
				else             {max_eq(v1, (c.pos + c.radius + 2.0f*spider.radius));} // spider to the right
			}
			if (val < v1) {spider.dir[t.dim] =  1.0;} // off low  end - reverse
			if (val > v2) {spider.dir[t.dim] = -1.0;} // off high end - reverse
			if (!spider.squished) {spider.move(timestep);} // can't be squished?
		}
	}
}

class surface_orienter_t {
	vector3d size;
	cube_t check_cube;
	vect_cube_t colliders;
public:
	cube_t get_check_cube() const {return check_cube;}

	void init(point const &p1, point const &p2, vector3d const &size_) {
		size = size_;
		check_cube.set_from_point(p1); // include both cur and last points
		check_cube.union_with_pt (p2);
		check_cube.expand_by(size.get_max_val()); // we could expand by size, but expanding by the max dim is more conservative and possibly better
		colliders.clear();
	}
	void register_cube(cube_t const &c, bool check_z_merge=0) {
		if (!c.intersects(check_cube)) return;
		
		if (check_z_merge) { // attempt to merge ceilings into floors
			for (cube_t &c2 : colliders) {
				if (c2.z1() == c.z2()) {c2.z1() = c.z1(); return;}
			}
		}
		colliders.push_back(c); // record for later processing in align_to_surfaces()
	}
	void register_cubes(cube_t const *const cubes, unsigned num) {
		for (unsigned n = 0; n < num; ++n) {register_cube(cubes[n]);}
	}
	template<typename T> void register_cubes(vector<T> const &cubes, bool check_z_merge=0) { // cube, cube_with_ix_t, etc.
		for (cube_t const &c : cubes) {register_cube(c, check_z_merge);}
	}
	void clip_and_max_expand_cubes() {
		if (colliders.size() <= 1) return; // not needed
		for (cube_t &c : colliders) {c.intersect_with_cube_xy(check_cube);}
		apply_fc_cube_max_merge_xy(colliders);
		sort_and_unique(colliders);
	}
	bool align_to_surfaces(spider_t &s, float delta_dir, rand_gen_t &rgen) {
		if (colliders.empty()) return 0; // floating in midair
		// Note: assumes last_pos is valid and non-intersecting; may not hold for initial placement or when objects are moved
		float const dist(p2p_dist(s.pos, s.last_pos)), r_inner(size.get_min_val()), r_outer(size.get_max_val());
		vector3d dir, best_dir(s.dir), best_up(s.upv);
		point best_pos(s.pos);
		float best_score(0.0); // distance squared

		// generate a number of possible new vectors in the forward direction, perform sphere-cube intersection with each one, and choose the new location closest to the target
		for (unsigned n = 0; n <= 50; ++n) { // include forward on first iteration
			if (n == 0) {dir = s.dir;}
			else { // maybe this should be more consistent, such as all multiples of 45 degrees from "forward"
				dir = dir_gen_xyz.gen_dir(); // random XYZ direction
				if (dot_product(dir, s.dir) < 0.0) {dir.negate();} // face forward
			}
			point cand(s.last_pos + dist*dir);
			bool had_coll(0);
			for (auto const &c : colliders) {had_coll |= ellipse_cube_intersect(cand, 1.05*size, c);} // slightly larger radius to pick up nearby objects
			if (!had_coll) continue; // floating in space, not a valid dir
			had_coll = 0;
			vector3d coll_norm(s.upv);
			for (auto const &c : colliders) {had_coll |= sphere_cube_int_update_pos(cand, r_inner, c, s.last_pos, 0, &coll_norm);}

			if (n == 0 && !had_coll) { // forward vector, no coll
				// rerun with outer radius to ensure coll_norm is correct
				for (auto const &c : colliders) {point tmp(cand); sphere_cube_int_update_pos(tmp, r_outer, c, s.last_pos, 0, &coll_norm);}
			}
			// squared movement distance, with higher weight for moving in the correct direction, and a preference for moving orthogonal in upv
			float const dir_dp(dot_product(dir, s.dir)), up_dp(fabs(dot_product(dir, s.upv)));
			float const score(p2p_dist_sq(cand, s.last_pos)*(dir_dp + 0.25*up_dp));
			if (score > best_score && coll_norm != zero_vector) {best_dir = dir; best_up = coll_norm; best_pos = cand; best_score = score;}
			if (n == 0 && !had_coll) break; // forward direction is still valid, done (early termination optimization)
		} // for n
		if (best_score == 0.0) { // no valid directions, must be floating in space
			if (colliders.size() > 1) return 0; // multiple colliders is too complex to handle here
			point const contact_pt(colliders.back().closest_pt(s.pos));
			best_up  = (s.pos - contact_pt).get_norm(); // point away from the collider
			if (best_up == zero_vector) {best_up = plus_z;} // handle invalid upv (s.pos == contact_pt); can this happen?
			best_pos = contact_pt + r_inner*best_up;
			orthogonalize_dir(best_dir, best_up, best_dir, 1);
			delta_dir *= 0.25; // update more smoothly for stability
		}
		vector3d const delta(best_pos - s.pos);
		float const delta_mag(delta.mag());
		s.pos = ((delta_mag > s.radius) ? (s.pos + (s.radius/delta_mag)*delta) : best_pos); // limit movement distance to radius
		s.dir = ((delta_dir*best_dir + (1.0 - delta_dir)*s.dir)).get_norm(); // slowly reorient
		s.upv = ((delta_dir*best_up  + (1.0 - delta_dir)*s.upv)).get_norm(); // slowly reorient
		orthogonalize_dir(s.dir, s.upv, s.dir, 1);
		return 1;
	}
}; // surface_orienter_t

surface_orienter_t surface_orienter; // reused across spiders

class obj_avoid_t {
	point &pos;
	point const &p_last;
	float radius;
public:
	bool had_coll;
	obj_avoid_t(point &pos_, point const &p_last_, float radius_) : pos(pos_), p_last(p_last_), radius(radius_), had_coll(0) {}
	void register_avoid_cube(cube_t const &c) {had_coll |= sphere_cube_int_update_pos(pos, radius, c, p_last);}
};

bool building_t::update_spider_pos_orient(spider_t &spider, point const &camera_bs, float timestep, rand_gen_t &rgen) const {
	assert(interior);
	float const trim_thickness(get_trim_thickness()), coll_radius(spider.get_xy_radius());
	vector3d const size(spider.get_size());
	bool const in_attic(point_in_attic(spider.pos));
	obj_avoid_t obj_avoid(spider.pos, spider.last_pos, coll_radius); // use xy_radius for all dims
	surface_orienter.init(spider.pos, spider.last_pos, size);

	if (spider.in_tank) { // pet store tank spider
		assert(spider.tunnel_tank_ix < interior->room_geom->objs.size());
		room_object_t const &tank(interior->room_geom->objs[spider.tunnel_tank_ix]);

		if (tank.type != TYPE_FISHTANK) { // tank was taken, set the spider free
			spider.in_tank = 0;
			return 0;
		}
		bool const has_lid(tank.has_lid());

		if (!has_lid && spider.pos.z > tank.z2()) { // climb out of the open top
			spider.in_tank = 0;
			return 0;
		}
		cube_t cubes[7];
		unsigned const num_cubes(get_fishtank_coll_cubes(tank, cubes));
		surface_orienter.register_cubes(cubes, num_cubes);
		cube_t clamp_area(tank);
		clamp_area.expand_by(-spider.radius);
		if (!has_lid) {clamp_area.z2() = tank.z2() + spider.radius;} // don't clamp on the top if the tank is open
		clamp_area.clamp_pt(spider.pos); // restrict to the tank, just in case they can get out
	}
	else { // free building spider
		// Note: we can almost use fc_occluders, except this doesn't contain the very bottom floor because it's not an occluder, and maybe the overlaps would cause problems
		surface_orienter.register_cubes(interior->floors);
		surface_orienter.register_cubes(interior->ceilings, 1); // check_z_merge=1
		surface_orienter.clip_and_max_expand_cubes(); // required to remove false splits between ceilings and floors
		if (has_room_geom()) {surface_orienter.register_cubes(interior->room_geom->glass_floors);}
	
		if (!in_attic) {
			for (unsigned d = 0; d < 2; ++d) {surface_orienter.register_cubes(interior->walls[d]);} // XY walls
		}
		cube_t const tc(surface_orienter.get_check_cube());

		if (!in_attic) { // check interior doors
			for (auto const &ds : interior->door_stacks) {
				if (ds.z1() > tc.z2() || ds.z2() < tc.z1()) continue; // wrong floor for this stack/part
				// calculate door bounds for bcube test, assuming it's open
				cube_t door_bounds(ds);
				door_bounds.expand_by_xy(ds.get_width());
				if (!tc.intersects(door_bounds)) continue; // optimization
				assert(ds.first_door_ix < interior->doors.size());

				for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
					door_t const &door(interior->doors[dix]);
					if (!ds.is_same_stack(door)) break; // moved to a different stack, done
					if (door.z1() > tc.z2() || door.z2() < tc.z1()) continue; // wrong floor

					if (door.open) { // how to handle open doors? they're not cubes; avoid them entirely? use their bcubes?
						cube_t door_bcube(get_door_bounding_cube(door));
						bool const dir(door.get_check_dirs()); // side of the door frame the door opens to
						door_bcube.d[!door.dim][!dir] += (dir ? 1.0 : -1.0)*2.0*spider.radius; // shift edge away from door frame to allow spider to walk on inside of the frame
				
						if (tc.intersects(door_bcube)) {
							obj_avoid.register_avoid_cube(door_bcube);
							// while the extruded polygon collision check below is more accurate, it can result in spiders getting stuck behind open doors
							//if (sphere_ext_poly_intersect(door_tq.pts, 4, normal, spider.pos, coll_radius, door.get_thickness(), 0.0)) {obj_avoid.had_coll = 1;}
						}
					}
					else if (tc.intersects(door)) {surface_orienter.register_cube(door);} // closed door
				} // for dix
			} // for door_stacks
		}
		// check interior objects
		static vect_cube_t cubes, avoid;
		vect_room_object_t::const_iterator b, e;
		get_begin_end_room_objs_on_ground_floor(tc.z2(), 1, b, e); // for_spider=1

		for (auto i = b; i != e; ++i) {
			if (i->z1() > tc.z2()) continue; // object is too high
			if (spider.on_web && spider.web_dir == 0 && i->type == TYPE_LIGHT && i->shape == SHAPE_CUBE) continue; // ignore cube lights if descending on a web to avoid getting stuck
			if (i->z2() < tc.z1() && !(i->type == TYPE_DESK && i->shape == SHAPE_TALL))                  continue; // sigh, tall desks are special
			if (!tc.intersects_xy((i->type == TYPE_CLOSET) ? get_true_room_obj_bcube(*i) : (cube_t)*i))  continue; // no intersection with this object
			if (!i->is_spider_collidable()) continue;
			if (i->get_max_dim_sz() < spider.radius) continue; // too small, skip
			get_room_obj_cubes(*i, spider.pos, cubes, cubes, avoid); // climb on large and small objects and avoid non-cube objects
			surface_orienter.register_cubes(cubes);
			for (cube_t const &c : avoid) {obj_avoid.register_avoid_cube(c);}
			cubes.clear();
			avoid.clear();
		} // for i
		if (has_pool()) {obj_avoid.register_avoid_cube(interior->pool);}
		// check elevators and escalators
		surface_orienter.register_cubes(interior->elevators); // should we avoid open elevators?
	
		for (escalator_t const &e : interior->escalators) {
			if (!tc.intersects(e)) continue;
			cube_t cubes[7];
			e.get_all_cubes(cubes);
			surface_orienter.register_cubes(cubes, 7);
			obj_avoid.register_avoid_cube(e.get_ramp_bcube(0)); // avoid ramp since it's not a cube; exclude_sides=0
		}
		if (in_attic) {
			// the attic roof is not a cube we can walk on and the beams aren't real objects;
			// also, it's not really possible to move from the attic floor to the roof without getting stuck in/on the corner
		}
		else if (!is_cube()) {
			// not cube walls, skip exterior
		}
		else { // check exterior walls; the spider can walk on them, but only if they're closed
			float const wall_thickness(get_wall_thickness()), door_open_dist(get_door_open_dist());
			auto const parts_end(get_real_parts_end_inc_sec());
			point const query_pt(get_inv_rot_pos(camera_bs)); // for open exterior door tests

			for (auto i = parts.begin(); i != parts_end; ++i) {
				for (unsigned dim = 0; dim < 2; ++dim) {
					for (unsigned dir = 0; dir < 2; ++dir) {
						cube_t cube(*i);
						cube.d[dim][!dir] = i->d[dim][dir] + (dir ? -1.0 : 1.0)*trim_thickness; // shrink to trim thickness
						if (!cube.intersects(tc)) continue; // optimization
						cubes.clear();
						cubes.push_back(cube); // start with entire length

						for (auto j = parts.begin(); j != get_real_parts_end(); ++j) { // clip against other parts
							if (j == i) continue; // skip self
							subtract_cube_from_cubes(*j, cubes); // subtract this part from current cubes by clipping in XY
						}
						if (is_basement(i) && has_ext_basement()) { // check for extended basement door and exclude it if open
							door_t const &door(interior->get_ext_basement_door());
						
							if (door.open && door.z1() < tc.z2() && door.z2() > tc.z1()) { // open door that overlaps spider zval
								cube_t wall_cut(door);
								wall_cut.expand_in_dim(door.dim, wall_thickness);
								subtract_cube_from_cubes(wall_cut, cubes, nullptr, 1); // clip_in_z=1
							}
						}
						for (cube_t const &c : cubes) {
							if (!c.intersects(tc)) continue; // optimization

							if (has_complex_floorplan && spider.pos.z > ground_floor_z1) { // handle intersecting parts
								// it's not easy to clip this wall to the other parts, so instead we'll project the spider pos to the other side of the wall
								// and see if it's outside the building; if so, it's a true exterior wall
								if (c.z1() <= spider.pos.z && c.z2() >= spider.pos.z && c.d[!dim][0] <= spider.pos[!dim] && c.d[!dim][1] >= spider.pos[!dim]) {
									point proj_pos(spider.pos);
									proj_pos[dim] = i->d[dim][dir] + (dir ? 1.0 : -1.0)*wall_thickness;
									bool contained(0);

									for (auto j = parts.begin(); j != parts_end; ++j) {
										if (j != i && j->contains_pt(proj_pos)) {contained = 1; break;}
									}
									if (contained) continue; // not a true wall, skip
								}
							}
							// check for open exterior doors
							bool on_open_door(0);

							for (tquad_with_ix_t const &door : doors) {
								cube_t door_bc(door.get_bcube());
								if (!door_bc.contains_pt_exp(query_pt, door_open_dist)) continue; // not open
								door_bc.expand_in_dim(dim, spider.radius);
								if (door_bc.intersects(tc)) {on_open_door = 1; break;}
							}
							if (!on_open_door) {surface_orienter.register_cube(c);}
						} // for c
					} // for dir
				} // for dim
			} // for i
		}
		if (obj_avoid.had_coll) {
			if (spider.on_web && spider.web_dir == 0) { // collided with un unwalkable object while descending on a web
				// if coll is ignored, spider will clip through the object; if spider stops descending, it will get stuck; so instead climb back up the web
				spider.web_dir = 1;
				// if pointing nearly straight down, add some randomness to the XY component so that we can turn around
				if (spider.dir.z < -0.95) {spider.dir += rgen.signed_rand_vector_xy(0.25); spider.dir.normalize();}
			}
			spider.end_jump();
		}
	} // end free building spider case
	float const delta_dir(min(1.0f, 1.5f*(1.0f - pow(0.7f, timestep))));

	if (!surface_orienter.align_to_surfaces(spider, delta_dir, rgen)) { // not on a surface
		if (!spider.is_jumping()) { // if jumping, we continue the jump; otherwise, drop to a surface below using a web
			if (!spider.on_web) {
				spider.web_start_zval = max(spider.pos.z, spider.last_pos.z) + spider.get_xy_radius();
				spider.on_web         = 1;
				spider.web_dir        = 0; // descending
			}
			float const dz_sign(spider.web_dir ? 1.0 : -1.0);
			spider.pos    = spider.last_pos;
			spider.dir    = ((delta_dir*dz_sign)*plus_z + (1.0 - delta_dir)*spider.dir).get_norm(); // slowly reorient to +/- z
			spider.pos.z += 0.5*dz_sign*timestep*spider.speed; // drop at half speed
			orthogonalize_dir(spider.upv, spider.dir, spider.upv, 1);
		}
	}
	else { // on a surface
		spider.end_jump();
		spider.on_web = 0;
	}
	return obj_avoid.had_coll;
}

// Note: non-const because biting the player can add blood and the player's inventory items to the building
void building_t::maybe_bite_and_poison_player(point const &pos, point const &camera_bs, vector3d const &dir,
	float coll_radius, float damage, int poison_type, rand_gen_t &rgen)
{
	if (!in_building_gameplay_mode()) return;
	if (dot_product_ptv(dir, camera_bs, pos) < 0.0) return; // facing the wrong direction
	if (get_floor_for_zval(camera_bs.z) != get_floor_for_zval(pos.z)) return; // wrong floor

	if (dist_xy_less_than(pos, camera_bs, (get_scaled_player_radius() + coll_radius))) { // do damage when nearly colliding with the player
		bool const played_sound(play_attack_sound(local_to_camera_space(pos), 0.5, 1.5, rgen)); // quieter and higher pitch than rats
		
		if (played_sound) { // bite and maybe poisoned every so often
			bool const player_dead(player_take_damage(damage, (poison_type > 0), poison_type)); // scream if poisoned
			if (player_dead) {register_player_death(cur_player_building_loc.pos);}
		}
	}
}

// Note: non-const because biting the player can add blood and the player's inventory items to the building
void building_t::update_spider(spider_t &spider, point const &camera_bs, float timestep, float &max_xmove, rand_gen_t &rgen) {
	
	if (spider.squished || spider.is_sleeping()) return; // horribly squished or peacefully sleeping
	float const radius(spider.radius), height(2.0*radius), coll_radius(2.0f*radius);
	float const spider_z1(spider.pos.z - height), spider_z2(spider.pos.z + height);

	if (spider.is_jumping()) { // jumping
		spider.jump_dist += p2p_dist_xy(spider.pos, spider.last_pos);
		spider.pos.z     += timestep*spider.jump_vel_z;
		apply_building_gravity(spider.jump_vel_z, timestep);
		// if we're still in the initial phase of our jump, skip collision and movement logic below to avoid initial coll with our starting object
		if (spider.jump_dist < 2.0*spider.radius && is_pos_inside_building(spider.pos, radius, radius)) return;
	}
	rgen.rand_mix(); // make sure it's different per spider
	// Note: we use a small xy_pad to allow the spider to climb on exterior walls
	if (spider.speed == 0.0) {set_spider_speed(spider, rgen, get_window_vspace());} // random speed

	if (!is_pos_inside_building(spider.pos, radius, radius)) {
		spider.end_jump();
		spider.pos    = spider.last_pos; // restore previous pos before collision
		spider.on_web = 0;

		if (!is_pos_inside_building(spider.pos, radius, radius)) { // still not valid
			if (spider.last_valid_pos == all_zeros) { // bad spawn pos - retry
				spider.pos = gen_animal_floor_pos(radius, spider_t::allow_in_attic(), 1, 0, 1, rgen); // not_player_visible=1, pref_dark_room=0, not_by_ext_door=1
				return;
			}
			spider.pos = spider.last_valid_pos; // restore to prev frame pos
		}
		assert(is_pos_inside_building(spider.pos, radius, radius));
		spider.choose_new_dir(rgen);
		return; // or could continue below?
	}
	spider.last_valid_pos = spider.pos;
	bool const had_coll(update_spider_pos_orient(spider, camera_bs, timestep, rgen));
	if (had_coll) {spider.end_jump();}

	// regenerate dir if collided and not on a web, or if dir is somehow bad
	if ((had_coll && !spider.on_web) || spider.dir.mag_sq() < 0.25) {spider.choose_new_dir(rgen);}
	else if ((float)tfticks > spider.update_time) { // direction change or sleep
		if (spider.on_web) {
			spider.update_time = (float)tfticks + 1.0*TICKS_PER_SECOND; // wait another 1s before updating
		}
		else if (spider.dist_since_sleep > 2.0*get_window_vspace() && rgen.rand_bool()) { // 50% chance of taking a rest
			spider.sleep_for(0.1, 5.0); // 0.1-5s
			spider.speed = 0.0; // will reset anim_time in the next frame
		}
		else {
			spider.update_time = (float)tfticks + rand_uniform(5.0, 15.0)*TICKS_PER_SECOND; // 5-15s
			vector3d const prev_dir(spider.dir);
			spider.choose_new_dir(rgen);
			spider.dir = (spider.dir + prev_dir).get_norm(); // 50% mix of prev and new dir to avoid sharp turns
		}
	}
	if (spider.on_web) { // rotate around the vertical axis
		bool const rotate_dir(spider.id & 1);
		float const angle(0.04*timestep*(rotate_dir ? 1.0 : -1.0));
		rotate_vector3d(plus_z, angle, spider.upv);
		spider.dir.normalize();
	}
	vector3d coll_dir;
	point const prev_pos(spider.pos); // capture the pre-collision point

	if (!spider.is_jumping() && check_and_handle_dynamic_obj_coll(spider.pos, spider.pos, coll_radius, spider_z1, spider_z2, camera_bs, 1)) { // check for collisions; for_spider=1
		spider.end_jump();
		coll_dir = (prev_pos - spider.pos).get_norm(); // points toward the collider in the XY plane

		// check if new pos is valid, and has a path to dest
		if (!is_pos_inside_building(spider.pos, radius, radius)) {
			spider.pos = prev_pos; // restore previous pos before collision
			spider.sleep_for(0.1, 0.2); // wait 0.1-0.2s so that we don't immediately collide and get pushed out again
		}
		else {
			max_eq(max_xmove, fabs(spider.pos.x - prev_pos.x));
			// update the path about every 30 frames of colliding; this prevents the spider from being stuck while also avoiding jittering due to frequent dest updates;
			// using randomness rather than updating every 30 frames makes this process less regular and mechanical
			if ((rgen.rand()%30) == 0) {
				vector3d const old_dir(spider.dir);
				spider.choose_new_dir(rgen);
				try_resolve_coll(old_dir, spider.upv, coll_dir, spider.dir, rgen.rand_bool()); // try_tangent=random
			}
		}
	}
	if (spider.upv.z > 0.5) { // handle biting the player, but not if upside down on the ceiling
		maybe_bite_and_poison_player(spider.pos, camera_bs, spider.dir, coll_radius, 0.1, 1, rgen); // 0.1 damage with poison
	}
}

bool building_t::maybe_squish_animals(room_object_t const &obj, point const &player_pos) { // spiders and cockroaches
	assert(has_room_geom());
	bool any_squished(0);

	for (spider_t &spider : interior->room_geom->spiders) {
		if (spider.squished) continue; // already squished
		if (!obj.contains_pt_xy(spider.pos) || !obj.intersects(spider.get_bcube())) continue;
		if (obj.get_size().get_max_val() < spider.get_xy_radius()) continue; // object is too small to squish this spider
		add_blood_decal(spider.pos, 1.5*spider.get_xy_radius(), colorRGBA(0.4, 1.0, 0.2, 1.0)); // yellow-green
		spider.pos.z -= 0.4*spider.get_height(); // move it near the ground since it will be drawn flattened
		any_squished  = spider.squished = 1;
		register_achievement("Splat the Spider");
	} // for spider
	for (insect_t &insect : interior->room_geom->insects) { // Note: no size check, achievement, or height reduction
		if (insect.squished || insect.type != INSECT_TYPE_ROACH) continue; // already squished, or not a cockroach
		if (!obj.intersects(insect.get_bcube())) continue;

		if (!obj.contains_pt_xy(insect.pos)) { // partial intersection: push cockroach out of the way
			point const orig_pos(insect.pos);
			if (sphere_cube_int_update_pos(insect.pos, insect.radius, obj, insect.last_pos)) {insect.last_pos = orig_pos;}
			continue;
		}
		add_blood_decal(insect.pos, 1.6*insect.get_xy_radius(), colorRGBA(1.0, 1.0, 0.2, 1.0)); // yellow
		any_squished = insect.squished = 1;
	} // for insect
	if (any_squished) {gen_sound_thread_safe(SOUND_SQUISH, (obj.get_cube_center() + (get_camera_pos() - player_pos)));} // get xlate from player delta
	return any_squished;
}


// *** Snakes ***

snake_t::snake_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_) : building_animal_t(pos_, radius_, dir_, id_) {
	rand_gen_t rgen;
	rgen.set_state(id+1, 3*id+7);
	rgen.rand_mix();
	length  = 2.0*radius; // input radius is half length
	xy_radius = radius; // will be updated later
	radius *= 0.04;
	color   = WHITE*(1.0 - rgen.rand_float()*rgen.rand_float()); // random color shade, weighted toward lighter
	has_rattle = rgen.rand_bool();
	unsigned const NUM_SEGS = 20; // head + 18 segments + tail
	float const seg_length(length/NUM_SEGS);
	vector3d const seg_step(-seg_length*dir); // head -> tail
	segments.resize(NUM_SEGS, pos);
	for (unsigned n = 0; n < NUM_SEGS; ++n) {segments[n] += seg_step*(n - NUM_SEGS/2.0 + 0.5);} // set segment centers
}
void snake_t::calc_xy_radius() {
	float xy_radius_sq(0.0);
	for (point const &p : segments) {max_eq(xy_radius_sq, p2p_dist_xy_sq(pos, p));}
	xy_radius = sqrt(xy_radius_sq) + radius; // don't forget to add the body radius
}
cube_t snake_t::get_bcube() const {
	cube_t bcube;
	bcube.set_from_points(segments);
	bcube.expand_by_xy(radius);
	bcube.z2() += radius;
	return bcube;
}
// Note: seg_ix is a float so that we can specify fractional segments; actual segment index is offset by 0.5 for the segment center
float snake_t::get_seg_radius(float seg_ix) const {
	assert(!segments.empty());
	float const t(seg_ix/segments.size()); // varies from 0.0 at tip of nose to 1.0 at tip of tail
	
	if (t < 0.4) { // head
		float const v(t/0.4); // 0-1
		return 0.5*(1.0 + v*v)*radius; // 0.5 => 1.0
	}
	if (t > 0.6) { // tail
		float const v((t - 0.6)/0.4); // 0-1
		return sqrt(1.0 - v)*radius; // 1.0 => 0.0
	}
	return radius; // middle of body
}
void snake_t::move_segments(float dist) {
	if (dist == 0.0) return;
	assert(dist > 0.0);
	assert(segments.size() >= 2);
	float const seg_len(get_seg_length());
	segments[0] += dist*dir; // move the head forward incrementally

	// work backwards toward the tail and apply the new segment directions
	for (unsigned n = 1; n < segments.size(); ++n) {
		vector3d const seg_delta(segments[n-1] - segments[n]);
		segments[n] = segments[n-1] - seg_delta*(seg_len/seg_delta.mag());
	}
	// update pos as new bcube center; pos.z remains unchanged
	cube_t const bcube(get_bcube());
	pos.x = bcube.xc();
	pos.y = bcube.yc();
	calc_xy_radius();
}
vector2d v2_from_v3_xy(vector3d const &v) {return vector2d(v.x, v.y);}

// Note: ignores radius; seg_dir points toward the head
bool snake_t::check_line_int_xy(point const &p1, point const &p2, bool skip_head, vector3d *seg_dir) const {
	for (unsigned n = (skip_head ? 2 : 1); n < segments.size(); ++n) {
		point const &s1(segments[n-1]), &s2(segments[n]);
		if (!line_segs_intersect_2d(v2_from_v3_xy(p1), v2_from_v3_xy(p2), v2_from_v3_xy(s1), v2_from_v3_xy(s2))) continue;
		if (seg_dir) {*seg_dir = (s1 - s2).get_norm();}
		return 1;
	}
	return 0;
}
bool snake_t::check_sphere_int(point const &sc, float sr, bool skip_head, vector3d *seg_dir, point *closest_pos) const {
	float const r_sum(sr + radius), r_sum_sq(r_sum*r_sum), r_ext(r_sum + get_seg_length()), r_ext_sq(r_ext*r_ext);

	for (unsigned n = (skip_head ? 2 : 1); n < segments.size(); ++n) {
		point const &s1(segments[n-1]), &s2(segments[n]);
		if (p2p_dist_sq(s1, sc) > r_ext_sq)                 continue; // optimization
		if (!sphere_test_comp(s1, sc, (s1 - s2), r_sum_sq)) continue;
		if (seg_dir    ) {*seg_dir = (s1 - s2).get_norm();}
		if (closest_pos) {*closest_pos = get_closest_pt_on_line(sc, s1, s2);}
		return 1;
	}
	return 0;
}
bool snake_t::detailed_sphere_coll(point const &sc, float sr, point &coll_pos, float &coll_radius) const {
	if (!check_sphere_int(sc, sr, 0, nullptr, &coll_pos)) return 0; // skip_head=0
	coll_radius = radius;
	return 1;
}
float snake_t::get_curve_factor() const {
	float curve_factor(0.0);

	for (auto i = segments.begin()+1; i+1 != segments.end(); ++i) {
		curve_factor += cross_product((*i - *(i-1)), (*(i+1) - *i)).z; // Warning: doesn't handle acute angles
	}
	float const seg_len(get_seg_length());
	return curve_factor / (seg_len*seg_len); // normalize based on segment length
}

void building_t::update_snakes(point const &camera_bs, unsigned building_ix) {
	vect_snake_t &snakes(interior->room_geom->snakes);
	if (snakes.placed && snakes.empty()) return; // no snakes placed in this building
	//timer_t timer("Update Snakes");
	add_animals_on_floor(snakes, building_ix, global_building_params.num_snakes_min, global_building_params.num_snakes_max,
		global_building_params.snake_size_min, global_building_params.snake_size_max);
	// update snakes
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	snakes.do_sort(); // is this necessary?
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame
	for (snake_t &snake : snakes) {update_snake(snake, camera_bs, timestep, snakes.max_xmove, rgen);}
}

// Note: non-const because biting the player can add blood and the player's inventory items to the building
void building_t::update_snake(snake_t &snake, point const &camera_bs, float timestep, float &max_xmove, rand_gen_t &rgen) {
	if (snake.is_sleeping()) return; // peacefully sleeping, no collision or movement needed
	if (snake.speed == 0.0) {snake.speed = global_building_params.snake_speed*rgen.rand_uniform(0.5, 1.0);} // random speed
	float const lookahead_amt(1.0*snake.radius);
	float const dist(timestep*snake.speed);
	point const &head_pos(snake.get_head_pos()); // only the head moves
	vector3d const center_dz(0.0, 0.0, 0.5*snake.get_height());
	vector3d dir(snake.dir); // this will be the new snake direction
	bool dir_valid(0);

	for (unsigned n = 0; n < 100; ++n) { // 100 attempts to select a new direction to resolve a collision
		point const new_head_pos(head_pos + dir*dist); // move the head one timestep
		vector3d const lookahead(dir*lookahead_amt);
		point const query_pos(new_head_pos + lookahead); // apply lookahead
		vector3d coll_dir; // collision normal; points from the collider in the XY plane
		// coll_type: 0=no coll, 1=outside building, 2=static object, 3=dynamic object, 4=snake itself
		int const coll_type(check_for_snake_coll(snake, camera_bs, timestep, head_pos, query_pos, coll_dir));
		
		if (coll_type == 0 || (snake.stuck_counter > 60 && coll_type >= 3)) { // no collision, or stuck on a dynamic object or itself
			// dir is valid; check for coll ahead of us
			if (check_for_snake_coll(snake, camera_bs, timestep, head_pos, (query_pos + 2.0*lookahead), coll_dir)) {
				vector3d const side_dir(cross_product(dir, plus_z));
				vector3d cand_dir((dir + 0.2*rgen.rand_uniform(-1.0, 1.0)*side_dir).get_norm()); // minor rotation
				if (!check_for_snake_coll(snake, camera_bs, timestep, head_pos, (head_pos + cand_dir*(dist + 3.0*lookahead_amt)), coll_dir)) {dir = cand_dir;}
			}
			if (coll_type == 0) {snake.stuck_counter = 0;} // only reset stuck_counter if there was no coll
			dir_valid = 1;
			break; // done
		}
		// collision, try a new direction
		bool const use_coll_dir(coll_dir != zero_vector);
		float const min_allowed_dp(0.9 - 0.014*n); // 0.9 => -0.5; prefer a slight turn

		for (unsigned m = 0; m < 100; ++m) { // 100 attempts to chose a valid dir; should amost always be successful
			dir = rgen.signed_rand_vector_spherical_xy_norm();
			// if dir is too close to the X or Y axis, choose a new dir; this prevents head-on collisions with common axis aligned cubes such as walls;
			// maybe initial dir should follow this logic as well?
			if (max(fabs(dir.x), fabs(dir.y)) > 0.95) continue;
			if (use_coll_dir && dot_product(dir, coll_dir) > 0.0) {dir.negate();} // don't move toward the collider

			if (dot_product(dir, snake.dir) < min_allowed_dp) {
				if (use_coll_dir) continue; // dir must be opposite coll_dir, can't negate again
				dir.negate();
			}
			break; // done
		} // for m
	} // for n
	if (dir_valid) {
		// update direction
		snake.dir = dir;
		update_dir_incremental(snake.last_valid_dir, snake.dir, 1.0, timestep, rgen);
		// move snake forward
		float const prev_pos_x(snake.pos.x), move_dist(snake.move(timestep)); // snake moves here
		snake.move_segments(move_dist);
		max_eq(max_xmove, fabs(prev_pos_x - snake.pos.x));

		if (move_dist > 0.0) { // move head in a winding motion if moving
			vector3d const side_dir(cross_product(snake.dir, plus_z));
			float const speed_factor(snake.speed/global_building_params.snake_speed); // [0.5, 1.0]
			float const rot_amt(sin(0.1*snake.anim_time)); // rotation amount should be independent of speed
			rotate_vector3d(plus_z, 0.02*fticks*PI*rot_amt*speed_factor, snake.dir);
			snake.dir.normalize(); // is this needed?
		}
	}
	else { // stuck; can happen when two snakes collide with each other or when a snake forms an inner spiral with itself
		++snake.stuck_counter;
		//snake.pos = snake.last_pos; snake.get_head_pos() -= (snake.pos - snake.last_pos); // revert to prev valid pos - doesn't work
		// move away from coll dir?
		// ignore coll and proceed (if coll with dynamic object or self)?
		// allow hard turn?
	}
	maybe_bite_and_poison_player((head_pos + center_dz), camera_bs, snake.dir, 2.0*snake.radius, 0.5, (snake.has_rattle ? 2 : 1), rgen); // 0.5 damage, poison if has a rattle
}

void get_xy_dir_to_closest_cube_edge(point const &pos, cube_t const &c, vector3d &dir) { // for pos outside c
	float const dx1(fabs(c.x1() - pos.x)), dx2(fabs(c.x2() - pos.x)), dy1(fabs(c.y1() - pos.y)), dy2(fabs(c.y2() - pos.y));
	float dmin(0.0);
	if (1         ) {dmin = dx1; dir = -plus_x;}
	if (dx2 < dmin) {dmin = dx2; dir =  plus_x;}
	if (dy1 < dmin) {dmin = dy1; dir = -plus_y;}
	if (dy2 < dmin) {dmin = dy2; dir =  plus_y;}
}

// applies to snakes and flies
// return values: 0=no coll, 1=outside building, 2=static object, 3=dynamic object, 4=ourself (for snakes)
// coll_dir points in the direction of the collision, opposite the collision normal / the direction we want to move to avoid a collision (new pos - old pos)
int building_t::check_for_animal_coll(building_animal_t const &A, float hheight, float z_center_offset, bool on_floor_only, bool skip_player,
	point const &camera_bs, float timestep, point const &old_pos, point const &query_pos, vector3d &coll_dir) const
{
	if (!bcube.contains_pt_xy_exp(query_pos, -A.radius) && !interior->basement_ext_bcube.contains_pt_xy_exp(query_pos, -A.radius)) { // outside the building interior
		get_xy_dir_to_closest_cube_edge(query_pos, bcube, coll_dir); // find closest bcube edge to head_pos
		return 1; // outside building bcube
	}
	vector3d const center_dz(0.0, 0.0, z_center_offset);
	point const query_center_z(query_pos + center_dz);

	if (!is_pos_inside_building(query_pos, A.radius, hheight, 0)) { // inc_attic=0
		for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
			if (p->contains_pt(query_center_z)) {get_xy_dir_to_closest_cube_edge(query_pos, *p, coll_dir);}
		}
		return 1; // outside building
	}
	// coll_ret: 0=no coll, 1=d0 wall, 2=d1 wall, 3=closed door d0, 4=closed door d1, 5=open door, 6=stairs, 7=elevator, 8=exterior wall, 9=room object, 10=closet, 11=cabinet
	int const coll_ret(check_line_coll_expand((old_pos + center_dz), query_center_z, A.radius, hheight, !on_floor_only)); // for_spider=!on_floor_only

	if (coll_ret) {
		if (coll_ret == 1 || coll_ret == 3) {coll_dir.x = ((query_pos.x < old_pos.x) ? -1.0 : 1.0);} // dim=0 wall/door, separates in X
		if (coll_ret == 2 || coll_ret == 4) {coll_dir.y = ((query_pos.y < old_pos.y) ? -1.0 : 1.0);} // dim=1 wall/door, separates in Y
		return 2; // collision with static room object
	}
	if (query_pos.z != old_pos.z) { // non-horizontal (insect); check for ceiling and floor collisions
		vector3d const expand(A.radius, A.radius, hheight);
		cube_t line_bcube(old_pos, query_pos);
		line_bcube.expand_by(expand);
		if (line_int_cubes_exp(old_pos, query_pos, interior->ceilings, expand, line_bcube) ||
		    line_int_cubes_exp(old_pos, query_pos, interior->floors,   expand, line_bcube)) {coll_dir.z = ((query_pos.z < old_pos.z) ? -1.0 : 1.0); return 2;}
	}
	point query_pos_coll(query_pos); // may be updated below on collision

	if (check_and_handle_dynamic_obj_coll(query_pos_coll, A.pos, A.radius, query_pos.z, (query_pos.z + 2.0*hheight), camera_bs, 0, skip_player)) { // for_spider=0
		coll_dir = (query_pos - query_pos_coll).get_norm();
		return 3; // collision with dynamic object
	}
	return 0;
}
int building_t::check_for_snake_coll(snake_t const &snake, point const &camera_bs, float timestep, point const &old_pos, point const &query_pos, vector3d &coll_dir) const {
	float const hheight(0.5*snake.get_height());
	int const ret(check_for_animal_coll(snake, hheight, hheight, 1, 0, camera_bs, timestep, old_pos, query_pos, coll_dir)); // on_floor_only=1, skip_player=0
	if (ret) return ret;
	vector3d seg_dir;

	// check for self intersection in a line in front and a sphere for the side; skip_head=1
	if (snake.check_line_int_xy(old_pos, (query_pos + snake.dir*snake.radius), 1, &seg_dir) ||
		snake.check_sphere_int(query_pos, snake.radius, 1, &seg_dir))
	{
		// set coll_dir to the direction of the segment pointing toward the head - we don't want to go that way and get stuck in a spiral
		coll_dir = seg_dir;
		return 4; // collision with ourself
	}
	return 0;
}


// *** Insects ***

float get_cockroach_height_from_radius(float radius) {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_ROACH));
	return 2.0*radius*sz.z/max(sz.x, sz.y);
}
float insect_t::get_height() const {
	if (flies()) {return 2.0*radius;} // spherical
	if (type != INSECT_TYPE_ROACH || !building_obj_model_loader.is_model_valid(OBJ_MODEL_ROACH)) {return 0.4*radius;} // roach squished sphere
	return get_cockroach_height_from_radius(radius);
}
cube_t insect_t::get_bcube() const {
	if (type == INSECT_TYPE_ROACH) {return get_cube_height_radius(pos, radius, get_height());}
	cube_t bcube(pos);
	bcube.expand_by(radius); // default spherical (INSECT_TYPE_FLY)
	return bcube;
}
cube_t insect_t::get_bcube_with_dir() const {
	if (type == INSECT_TYPE_ROACH && building_obj_model_loader.is_model_valid(OBJ_MODEL_ROACH)) {
		return get_obj_model_bcube_for_dir(pos, dir, radius, get_height(), OBJ_MODEL_ROACH);
	}
	return get_bcube();
}

void building_t::update_insects(point const &camera_bs, unsigned building_ix) {
	vect_insect_t &insects(interior->room_geom->insects);
	if (insects.placed && insects.empty()) return; // no insects placed in this building
	//timer_t timer("Update Insects");
	// spawn insects on the floor and let them lift off and fly away if they're a flying type
	bool const was_placed(insects.placed);
	add_animals_on_floor(insects, building_ix, global_building_params.num_insects_min, global_building_params.num_insects_max,
		global_building_params.insect_size_min, global_building_params.insect_size_max);
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame

	if (!was_placed) { // newly placed (on the floor); shift by radius in Z
		for (insect_t &insect : insects) {
			// add flies to lit rooms and cockroaches to dark rooms
			int const room_id(get_room_containing_pt(insect.pos));
			if (is_room_lit(room_id, insect.get_z2())) {insect.type = INSECT_TYPE_FLY;} // make it a fly (which it already should be)
			else { // make it a cockroach
				insect.type    = INSECT_TYPE_ROACH;
				insect.radius *= 5.0; // larger than a fly
			}
			insect.pos.z += 0.5*insect.get_height(); // for insects, pos.z is the center rather than the bottom
		} // for insect
	}
	// update insects
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms

	for (insect_t &insect : insects) {
		if (!insect.squished) {insect.move(timestep);}
	}
	vector<pair<float, point>> targets; // used for flies
	for (insect_t &insect : insects) {update_insect(insect, camera_bs, timestep, targets, rgen);}
}

void building_t::update_insect(insect_t &insect, point const &camera_bs, float timestep, vector<pair<float, point>> &targets, rand_gen_t &rgen) const {
	// run logic that's common to all insects
	if (insect.squished) return;

	if (insect.speed == 0.0) { // generate initial speed and dir
		insect.speed = global_building_params.insect_speed*rgen.rand_uniform(0.5, 1.0);
		insect.dir   = rgen.signed_rand_vector_norm();
		if     (!insect.flies())     {insect.dir.z  =  0.0;} // moves in the XY plane
		else if (insect.dir.z < 0.0) {insect.dir.z *= -1.0;} // assume starting on the floor, and fly upward
	}
	if      (insect.type == INSECT_TYPE_FLY  ) {update_fly  (insect, camera_bs, timestep, targets, rgen);}
	else if (insect.type == INSECT_TYPE_ROACH) {update_roach(insect, camera_bs, timestep, rgen);}
	//else if (insect.type == INSECT_TYPE_CENTIPEDE) {update_centipede(insect, camera_bs, timestep, rgen);}
	else {assert(0);} // unsupported insect type
}

void building_t::update_fly(insect_t &fly, point const &camera_bs, float timestep, vector<pair<float, point>> &targets, rand_gen_t &rgen) const {
	fly.has_target = fly.target_player = 0;
	float const radius(fly.radius), hheight(0.5*fly.get_height());
	point &pos(fly.pos);

	// check for collision and bounce off the object for now
	if (!is_pos_inside_building(pos, radius, hheight)) {
		pos = fly.last_pos; // restore previous pos before collision

		if (!is_pos_inside_building(pos, radius, hheight)) { // still not valid, respawn; error?
			pos = gen_animal_floor_pos(radius, insect_t::allow_in_attic(), 1, 0, 0, rgen); // not_player_visible=1, pref_dark_room=0, not_by_ext_door=0
			return;
		}
		fly.dir *= -1.0; // reverse direction; maybe should use direction to closest building wall?
		fly.delta_dir = zero_vector; // reset delta_dir
		return; // continue below?
	}
	float const dist_to_player(p2p_dist(pos, camera_bs)), sight_range(2.0*get_window_vspace());
	bool const target_player((player_attracts_flies || player_wait_respawn) && dist_to_player < sight_range);
	unsigned const update_freq(1 + interior->room_geom->insects.size()/250); // reduced update rate for many insects

	if (((frame_counter + fly.id) % update_freq) == 0) { // run collision/update logic every few frames
		vector3d const lookahead(fly.dir*(2.0*radius));
		vector3d coll_dir; // opposite the collision normal; points toward the collider in the XY plane
		// coll_type: 0=no coll, 1=outside building, 2=static object, 3=dynamic object; z_center_offset=0.0, on_floor_only=0
		int const ret(check_for_animal_coll(fly, hheight, 0.0, 0, target_player, camera_bs, timestep, pos, (pos + lookahead), coll_dir));
	
		if (ret) { // collision
			// if coll_dir is not set, use the direction we last moved in, or our dir if that's zero
			if (coll_dir == zero_vector) {coll_dir = ((fly.last_pos == fly.pos) ? fly.dir : (fly.pos - fly.last_pos).get_norm());}

			if (ret == 2 && !target_player) { // static object collision
				// land on the object? we need to get the object index first (if a room object) and check that it's a cube; or only land on walls?
				// we would also need to change the fly's up vector from +z to align to a vertical surface, similar to how spiders are drawn
			}
			if (ret == 3) { // dynamic object, such as the player
				//fly.pos -= radius*fticks*coll_dir; // push away slowly, rather than spinning in place
				fly.dir = -coll_dir; // fly away
			}
			else { // static object
				fly.pos = fly.last_pos; // move back to a point where we didn't collide (assuming update is frequent enough)
				fly.dir = rgen.signed_rand_vector_norm();
				if (dot_product(fly.dir, coll_dir) > 0.0) {fly.dir.negate();} // don't fly toward the collider
			}
			fly.delta_dir = zero_vector; // reset delta_dir on coll
			return;
		}
		if (rgen.rand_float() < 0.25) { // every 4 updates send out a longer range collision query
			if (check_for_animal_coll(fly, hheight, 0.0, 0, target_player, camera_bs, timestep, pos, (pos + 8.0*lookahead), coll_dir)) { // Note: coll_dir is unused
				fly.delta_dir *= 0.9f; // reduce direction change
				fly.dir       += 0.25*rgen.signed_rand_vector_norm(); // adjust direction
				fly.dir.normalize();
				min_eq(fly.accel, 0.0f); // stop accelerating
			}
		}
	}
	// run player/zombie targeting logic
	targets.clear();
	if (target_player) {targets.emplace_back(dist_to_player, camera_bs);}
	bool follow_mode(0);
	
	if (in_building_gameplay_mode()) { // look for zombies to follow
		for (person_t const &p : interior->people) {
			point const eye_pos(p.get_eye_pos());
			float const dist_sq(p2p_dist_sq(pos, eye_pos));
			if (dist_sq > sight_range*sight_range) continue; // too far away to see
			targets.emplace_back(sqrt(dist_sq), eye_pos);
		}
	}
	sort(targets.begin(), targets.end()); // sort min to max distance

	for (auto const &target : targets) {
		bool const is_player(target.second == camera_bs);
		
		if (get_room_containing_pt(pos) != (is_player ? cur_player_building_loc.room_ix : get_room_containing_pt(target.second))) { // different rooms
			if (!is_pt_visible(pos, target.second))                 continue; // not visible
			if (check_line_int_interior_window(pos, target.second)) continue; // not reachable
		}
		fly.has_target = 1;
		if (is_player) {fly.target_player = 1;}
		follow_mode = (target.first > 1.2*get_scaled_player_radius()); // follow if not very close to the target
		vector3d const dir_to_target((target.second - pos).get_norm());
		update_dir_incremental_no_zero_check(fly.dir, dir_to_target, 0.5, timestep); // slow turn to target direction
	} // for target
	// apply a slow random dir change
	fly.delta_dir += (0.1f*timestep)*rgen.signed_rand_vector();
	fly.dir       += (0.1f*timestep)*fly.delta_dir;
	fly.dir.normalize();
	if (fabs(fly.dir.z) > 0.99) {fly.dir = fly.delta_dir;} // don't point straight up or straight down
	// apply a slow random acceleration
	float const max_speed(global_building_params.insect_speed);
	fly.accel += (0.04f*timestep)*rgen.signed_rand_float();
	fly.accel  = CLIP_TO_pm1(fly.accel);
	fly.speed  = (follow_mode ? 1.6 : 1.0)*min(max_speed, max(0.5f*max_speed, (fly.speed + (0.05f*timestep)*fly.accel))); // faster when following
	
	// play buzz sound if near player and also attracted to the player
	if ((player_attracts_flies || player_wait_respawn) && dist_to_player < 1.1*get_scaled_player_radius()) {
		gen_sound_thread_safe(SOUND_FLY_BUZZ, local_to_camera_space(pos), 1.0, 1.0, 1.0, 1); // skip_if_already_playing=1
	}
	// make sure we're not in front of the camera near clip plane
	float const camera_dmin(1.2*(radius + NEAR_CLIP));
	if (dist_to_player < camera_dmin) {pos = camera_bs + camera_dmin*(pos - camera_bs).get_norm();}
}

void building_t::update_roach(insect_t &roach, point const &camera_bs, float timestep, rand_gen_t &rgen) const {
	float const hheight(0.5*roach.get_height());
	float const radius(0.67*roach.radius); // use a smaller collision radius to allow roaches to partially enter a wall or object before disappearing
	point &pos(roach.pos);

	if (!roach.is_sleeping()) { // check for collisions
		bool spawn_new_pos(0);
		if (!is_pos_inside_building(pos, radius, hheight, 0)) {spawn_new_pos = 1;} // outside building, respawn, inc_attic=0

		if (!spawn_new_pos) {
			// coll_ret: 0=no coll, 1=d0 wall, 2=d1 wall, 3=closed door d0, 4=closed door d1, 5=open door, 6=stairs, 7=elevator, 8=exterior wall, 9=room object, 10=closet, 11=cabinet
			int const coll_ret(check_line_coll_expand(roach.last_pos, pos, radius, hheight, 0)); // for_spider=0
			bool maybe_stuck(0);

			if (coll_ret == 9) { // room object (except closet or kitchen cabinet); open door collisions are ignored (can go under doors)
				vector3d cnorm;
				float hardness(0.0);
				int obj_ix(-1);
				if (pos == roach.last_pos) {spawn_new_pos = 1;} // stuck?
				else if (interior->check_sphere_coll_room_objects(*this, roach.last_pos, roach.last_pos, radius, interior->room_geom->objs.end(), cnorm, hardness, obj_ix) &&
					dot_product(cnorm, roach.dir) >= 0.0)
				{
					// last_pos is also colliding, and dir is away from the collision - continue without direction change; applies to books dropped by the player
				}
				else {
					vector3d const new_dir(rgen.signed_rand_vector_spherical_xy_norm());
					roach.dir = (dot_product(roach.dir, new_dir) < 0.0) ? new_dir : -new_dir; // make sure it's at least 180 degrees different
					if (roach.last_pos != all_zeros) {pos = roach.last_pos;} // reset pos to prev frame
					roach.is_scared = 0; // distracted, no longer scared
					roach.no_scare  = 1;
					maybe_stuck     = 1;
				}
			}
			else if (coll_ret > 0 && coll_ret != 5) { // wall, closed door, stairs, elevator, or closet
				spawn_new_pos = 1; // disappear under it and respawn
			}
			if (!maybe_stuck) {roach.stuck_counter = 0;}
			else if (roach.stuck_counter++ > 60) {spawn_new_pos = 1;} // respawn if stuck colliding for 60 consecutive frames
		}
		if (spawn_new_pos) {
			pos = gen_animal_floor_pos(roach.radius, 0, 1, 1, 0, rgen); // place_in_attic=0, not_player_visible=1, pref_dark_room=1, not_by_ext_door=0
			roach.is_scared = roach.no_scare = 0; // no longer scared
			return;
		}
	}
	// Note: no need to check for dynamic object collisions since we can likely just go under them
	int room_id(-1);
	float const floor_spacing(get_window_vspace()), scare_dist(0.9*floor_spacing);
	point const camera_bot(camera_bs.x, camera_bs.y, camera_bs.z-get_bldg_player_height());
	point run_from;
	
	if (!roach.no_scare) { // run scare logic; similar to building_t::scare_rat()
		if (has_room_geom() && interior->room_geom->fire_manager.get_closest_fire(pos, 2.0*radius, pos.z-hheight, pos.z+hheight, &run_from)) { // fire
			roach.is_scared = 1;
		}
		else if (camera_surf_collide && dist_less_than(pos, camera_bot, scare_dist)) {
			run_from = camera_bot;
			roach.is_scared = 1;
		}
		else {
			for (person_t const &p : interior->people) { // other people in the building scare the rats
				if (p.is_waiting_or_stopped()) continue; // only scare if moving
				point const ppos(p.pos.x, p.pos.y, p.get_z1());
				if (dist_less_than(pos, ppos, scare_dist)) {run_from = ppos; roach.is_scared = 1; break;}
			}
			sphere_t const cur_sound(get_cur_frame_loudest_sound());
			if (cur_sound.radius > 0.0 && dist_less_than(pos, cur_sound.pos, 4.0*cur_sound.radius)) {run_from = cur_sound.pos; roach.is_scared = 1;}
			else if (!roach.is_scared && ((frame_counter + roach.id) & 3) == 0) { // check for light every 4 frames
				room_id = get_room_containing_pt(pos);
				roach.is_scared |= is_room_lit(room_id, roach.get_z2()); // run from the light
			}
		}
	}
	float const nom_speed((roach.is_scared ? 1.0 : 0.25)*global_building_params.insect_speed), max_speed(2.0*nom_speed);
	roach.speed += 0.01*max_speed*fticks*rgen.signed_rand_float(); // add a bit of random variation over time
	max_eq(roach.speed, nom_speed);
	min_eq(roach.speed, max_speed);

	if (roach.is_scared) {
		roach.wake_time = 0.0; // wake up
		if (run_from != all_zeros) {roach.dir = (pos - point(run_from.x, run_from.y, pos.z)).get_norm();} // run away
		else if (room_id >= 0) { // in a room, run toward the nearest wall
			cube_t const &room(get_room(room_id));
			roach.dir = room.closest_side_dir(pos, 4); // skip_dims=4 (z)
		}
		// else run quickly in the current direction
	}
	else { // slow random walk and stop
		if (!roach.is_sleeping() && roach.dist_since_sleep > roach.dist_to_sleep) { // sleep
			roach.sleep_for(0.0, 4.0); // 0-4s
			roach.no_scare  = 0; // allow scaring again
			roach.speed     = 0.0;
			roach.delta_dir = rgen.signed_rand_vector_spherical_xy_norm(); // choose a new random dir after sleep is over
			roach.dist_to_sleep = floor_spacing*rgen.rand_uniform(0.2, 1.0); // walk for a while before sleeping again
		}
		else if (roach.dir != roach.delta_dir) { // slowly turn
			update_dir_incremental(roach.dir, roach.delta_dir, 0.25, timestep, rgen);
		}
	}
}

void register_fly_attract(bool no_msg) {
	if (player_attracts_flies) return; // no change
	player_attracts_flies = 1;
	
	if (!no_msg && global_building_params.num_insects_max > 0 && player_building && !player_building->interior->room_geom->insects.empty()) {
		print_text_onscreen("You have attracted the flies", ORANGE, 1.0, 2.5*TICKS_PER_SECOND, 1);
	}
}

