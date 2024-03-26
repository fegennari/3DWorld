// 3D World - Building vs. Player/AI Interaction Logic
// by Frank Gennari 1/14/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "openal_wrap.h"

// physics constants, currently applied to balls
float const KICK_VELOCITY  = 0.0025;
float const MIN_VELOCITY   = 0.0001;
float const SMOKE_VELOCITY = 0.0006;
float const OBJ_DECELERATE = 0.008;
float const OBJ_GRAVITY    = 0.0003; // unsigned magnitude
float const TERM_VELOCITY  = 1.0;
float const OBJ_ELASTICITY = 0.8;

extern bool tt_fire_button_down, flashlight_on, use_last_pickup_object, city_action_key, player_wait_respawn;
extern int player_in_closet, camera_surf_collide, can_pickup_bldg_obj, animate2, frame_counter, player_in_elevator, player_in_attic;
extern float fticks, CAMERA_RADIUS, office_chair_rot_rate;
extern double tfticks;
extern building_dest_t cur_player_building_loc;
extern building_t const *player_building;


bool player_can_open_door(door_t const &door);
unsigned player_has_room_key();
bool player_has_pool_cue();
void register_broken_object(room_object_t const &obj);
void record_building_damage(float damage);
void refill_thirst();
colorRGBA get_glow_color(float stime, bool fade);
void play_hum_sound(point const &pos, float gain, float pitch);
bool ceiling_fan_is_on(room_object_t &obj, vect_room_object_t const &objs);
cube_t get_pool_table_top_surface(room_object_t const &c);

// Note: pos is in camera space
void gen_sound_thread_safe(unsigned id, point const &pos, float gain, float pitch, float gain_scale, bool skip_if_already_playing) {
	assert(gain > 0.0 && pitch > 0.0 && gain_scale > 0.0);
	float const dist(p2p_dist(get_camera_pos(), pos)), dscale(10.0*CAMERA_RADIUS*gain_scale); // distance at which volume is halved
	gain *= dscale/(dist + dscale);
	if (gain < 0.025) return; // too soft to hear
#pragma omp critical(gen_sound)
	gen_sound(id, pos, gain, pitch, 0, zero_vector, skip_if_already_playing);
}

// lights

float get_radius_for_room_light(room_object_t const &obj);

bool is_motion_detected(point const &activator, cube_t const &light, cube_t const &room, float floor_spacing) {
	return (room.contains_pt(activator) && activator.z < light.z1() && activator.z > (light.z2() - floor_spacing));
}

void building_t::run_light_motion_detect_logic(point const &camera_bs) {
	if (!animate2)             return;
	if (is_house || !interior) return; // office buildings only
	if (player_in_elevator)    return; // skip so that we don't have a lot of clicking when lights switch on due to AIs while passing floors
	float const floor_spacing(get_window_vspace());
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_LIGHT || !i->is_active() || !i->is_powered()) continue; // not a light, unpowered, or not motion activated
		assert(i->room_id < interior->rooms.size());
		room_t const &room(interior->rooms[i->room_id]);
		bool const is_player(is_motion_detected(camera_bs, *i, room, floor_spacing));
		bool activated(is_player);
		float &off_time(i->light_amt); // store auto off time in the light_amt field

		for (auto p = interior->people.begin(); p != interior->people.end() && !activated; ++p) {
			if (p->is_waiting_or_stopped()) continue; // skip if stopped/waiting
			if (fabs(p->pos.z - camera_bs.z) > 0.75*floor_spacing) continue; // player is on a different floor, skip (too many clicking sounds)
			activated |= is_motion_detected(p->pos, *i, room, floor_spacing);
		}
		if (activated) {
			off_time = tfticks + 10.0*TICKS_PER_SECOND; // automatically turn off 10s since last activation
			if (i->is_lit()) continue; // stays lit - no change
		}
		else {
			if (!animate2)          continue; // no auto off
			if (!i->is_lit())       continue; // already off, and stays off
			if (tfticks < off_time) continue; // already on, and not yet time to switch off
		}
		i->toggle_lit_state();
		set_obj_lit_state_to(i->room_id, i->z2(), i->is_lit());
		register_light_state_change(*i, i->get_cube_center());
		if (1 || is_player) {register_indir_lighting_state_change(i - interior->room_geom->objs.begin());} // only update for player activations?
	} // for i
}

// Note: called by the player; closest_to is in building space, not camera space
bool building_t::toggle_room_light(point const &closest_to, bool sound_from_closest_to, int room_id, bool inc_lamps, bool closet_light, bool known_in_attic) {
	if (!has_room_geom()) return 0; // error?
	// attic lights on posts are exactly between roof tquads and point_in_attic() may not return the correct value, so known_in_attic should be passed in
	bool const in_attic(known_in_attic || point_in_attic(closest_to));

	if (room_id < 0 && !in_attic) { // caller has not provided a valid room_id, so determine it now
		point query_pt(closest_to);
		if (is_rotated()) {do_xy_rotate_inv(bcube.get_cube_center(), query_pt);}
		room_id = get_room_containing_pt(query_pt);
		if (room_id < 0) return 0; // closest_to is not contained in a room of this building
	}
	bool const ignore_floor(in_attic || get_room(room_id).is_single_floor);
	vect_room_object_t &objs(interior->room_geom->objs);
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	bool const in_closet(closet_light || bool(player_in_closet)); // while in the closet, player can only toggle closet lights and not room lights
	float closest_dist_sq(0.0);
	int closest_light(-1);

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_type() || (!inc_lamps && i->type == TYPE_LAMP)) continue; // not a light
		if ( in_attic && !i->in_attic()) continue;
		if (!in_attic && i->room_id != room_id) continue; // wrong room
		if (i->in_closet() != in_closet) continue;
		if (i->in_elevator()) continue; // can't toggle elevator light
		if (!ignore_floor && get_floor_for_zval(i->z1()) != get_floor_for_zval(closest_to.z)) continue; // wrong floor (skip garages and sheds)
		point center(i->get_cube_center());
		if (is_rotated()) {do_xy_rotate(bcube.get_cube_center(), center);}
		float const dist_sq(p2p_dist_sq(closest_to, center));
		if (closest_dist_sq == 0.0 || dist_sq < closest_dist_sq) {closest_dist_sq = dist_sq; closest_light = int(i - objs.begin());}
	} // for i
	if (closest_light < 0) return 0;
	room_object_t const &light(objs[closest_light]);
	toggle_light_object(light, (sound_from_closest_to ? closest_to : light.get_cube_center()));
	return 1;
}

void building_t::toggle_light_object(room_object_t const &light, point const &sound_pos) { // called by the player
	assert(has_room_geom());
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	bool updated(0), is_lamp(1);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) { // toggle all lights on this floor of this room
		if (!i->is_light_type() || i->room_id != light.room_id || !i->is_powered()) continue;
		if (light.in_closet() != i->in_closet()) continue; // closet + room light are toggled independently
		if (i->z2() != light.z2()) continue; // Note: uses light z2 rather than z1 so that thin lights near doors are handled correctly
		i->toggle_lit_state(); // Note: doesn't update indir lighting
		i->flags &= ~RO_FLAG_IS_ACTIVE; // disable motion detection feature if the player manually switches lights off
		register_indir_lighting_state_change(i - interior->room_geom->objs.begin());
		updated = 1;
		if (i->type == TYPE_LAMP)  continue; // lamps don't affect room object ambient lighting, and don't require regenerating the vertex data, so skip the step below
		if (is_lamp) {set_obj_lit_state_to(light.room_id, light.z2(), i->is_lit());} // update object lighting flags as well, for first non-lamp light
		is_lamp = 0;
	} // for i
	if (!updated) return; // can we get here?
	register_light_state_change(light, sound_pos, is_lamp);
	//interior->room_geom->modified_by_player = 1; // should light state always be preserved?
}
void building_t::register_light_state_change(room_object_t const &light, point const &sound_pos, bool is_lamp) {
	if (!is_lamp) {interior->room_geom->invalidate_lights_geom();} // recreate light geom with correct emissive properties if not a lamp; deferred until next draw pass
	gen_sound_thread_safe(SOUND_CLICK, local_to_camera_space(sound_pos));
	register_building_sound(sound_pos, 0.1);
	float const fear_amt((light.is_light_on() ? 1.0 : 0.5)*(is_lamp ? 0.5 : 1.0)); // max fear from lights turning on; lamps are half as much fear

	for (rat_t &rat : interior->room_geom->rats) { // light change scares rats
		if (get_room_containing_pt(rat.pos) == light.room_id) {scare_rat_at_pos(rat, sound_pos, fear_amt, 0);} // scare if in the same room
	}
}

bool building_t::set_room_light_state_to(room_t const &room, float zval, bool make_on) { // called by AI people
	if (!has_room_geom()) return 0; // error?
	if (room.is_hallway)  return 0; // don't toggle lights for hallways, which can have more than one light
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	float const window_vspacing(get_window_vspace());
	bool updated(0);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_LIGHT) continue; // not a light (excludes lamps)
		if (i->flags & (RO_FLAG_IN_CLOSET | RO_FLAG_IN_ELEV)) continue; // skip lights in closets or elevators, these can only be toggled by the player
		if (i->z1() < zval || i->z1() > (zval + window_vspacing) || !room.contains_cube_xy(*i)) continue; // light is on the wrong floor or in the wrong room
		if (i->is_lit() != make_on) {i->toggle_lit_state(); updated = 1;} // Note: doesn't update indir lighting or room light value
	} // for i
	if (updated) {interior->room_geom->invalidate_lights_geom();} // recreate light geom with correct emissive properties; will flag for update next frame
	return updated;
}

void building_t::set_obj_lit_state_to(unsigned room_id, float light_z2, bool lit_state) {
	assert(has_room_geom());
	room_t const &room(get_room(room_id));
	float const light_intensity(room.light_intensity);
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	float const obj_zmin(room.is_single_floor ? room.z1() : (light_z2 - get_floor_ceil_gap()));

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->room_id != room_id || i->z1() < obj_zmin || i->z1() > light_z2) continue; // wrong room or floor
		if (i->is_obj_model_type()) continue; // light_amt currently does not apply to 3D models; should it?
		bool invalidate(0);

		if (i->type == TYPE_STAIR || i->type == TYPE_STAIR_WALL || i->type == TYPE_ELEVATOR || i->type == TYPE_LIGHT || i->type == TYPE_BLOCKER ||
			i->type == TYPE_COLLIDER || i->type == TYPE_SIGN || i->type == TYPE_WALL_TRIM || i->type == TYPE_RAILING || i->type == TYPE_BLINDS ||
			i->type == TYPE_SWITCH || i->type == TYPE_OUTLET || i->type == TYPE_PG_WALL || i->type == TYPE_PG_PILLAR || i->type == TYPE_PG_BEAM ||
			i->type == TYPE_PARK_SPACE || i->type == TYPE_RAMP || i->type == TYPE_VENT)
		{
			continue; // not a type that uses light_amt
		}
		else if (i->type == TYPE_WINDOW) {
			if (lit_state) {i->flags |= RO_FLAG_LIT;} else {i->flags &= ~RO_FLAG_LIT;}
			invalidate = 1;
		}
		else {
			float const prev_light_amt(i->light_amt);
			if (lit_state) {i->light_amt += light_intensity;} else {i->light_amt = max((i->light_amt - light_intensity), 0.0f);} // shouldn't be negative, but clamp to 0 just in case
			invalidate = (fabs(i->light_amt - prev_light_amt) > 0.1); // generally always true, but good to have this check/optimization in the future
		}
		if (invalidate) {interior->room_geom->invalidate_draw_data_for_obj(*i);} // Note: can't clear here if called from building AI (not in the draw thread)
	} // for i
}

bool building_room_geom_t::closet_light_is_on(cube_t const &closet) const {
	auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type == TYPE_LIGHT && i->in_closet() && closet.contains_cube(*i)) {return i->is_light_on();}
	}
	return 0;
}

breaker_zone_t building_interior_t::get_circuit_breaker_info(unsigned zone_id, unsigned num_zones, float floor_spacing) const {
	assert(zone_id < num_zones);

	// Note: if there are multiple panels, they will affect the same set of zones; it seems too difficult to assign rooms/zones across panels;
	// this means that zones will follow the state of the last breaker that was toggled to a different state
	if (!elevators.empty()) { // elevators are always zone 0 (lower left or right breaker)
		if (zone_id == 0) return breaker_zone_t(RTYPE_ELEVATOR, 0, 0);
		--zone_id; --num_zones; // exclude elevator
	}
	unsigned const num_rooms(rooms.size());
	if (num_zones == 0 || num_rooms == 0) return breaker_zone_t(); // no zones left, or no rooms
	// wrap around if there are more zones than rooms; this will assign the same room to multiple zones/breakers;
	// labels will use room names on upper floors for wrapped zones; maybe this should also split toggle behavior by floor groups?
	// but that would be extra complex because we would need per-floor room unpowered flags, and it's unclear if the player would even understand how this works
	unsigned const floor_block_ix(zone_id/num_rooms);
	zone_id = zone_id % num_rooms;
	// determine which rooms this breaker controls;
	// we really should have breakers control lights on separate floors rather than vertical rooms stacks, but this is much easier;
	// note that the first breaker/room (after the elevator) will be the primary hallway in office buildings and will also control all cameras,
	// and the last breaker will be for the extended basement/backrooms if there is one, otherwise the basement/parking garage
	float const rooms_per_zone(max(1.0f, float(num_rooms)/num_zones));
	unsigned const room_start(round_fp(zone_id*rooms_per_zone)), room_end(min((unsigned)round_fp((zone_id+1)*rooms_per_zone), num_rooms));
	if (room_start >= room_end) return breaker_zone_t(); // no rooms
	// pick a room with the highest priority for the label
	unsigned const room_priorities[NUM_RTYPES] = {0, 2, 1, 1, 2, 2, 3, 3, 3, 2, 1, 3, 2, 3, 3, 3, 2, 2, 2, 2, 3, 3, 0, 3, 3, 0, 4, 3, 4, 4, 4, 0};
	unsigned ret_rtype(0), highest_priority(0);

	for (unsigned r = room_start; r < room_end; ++r) {
		room_t const &room(rooms[r]);
		// use room on the first floor, since it's more likely to be special, unless the zones wrap around (more zones than rooms)
		unsigned const floor_ix(min(floor_block_ix, ((unsigned)round_fp(room.dz()/floor_spacing)-1)));
		unsigned const rtype(room.get_room_type(floor_ix));
		assert(rtype < NUM_RTYPES);
		unsigned const priority(room_priorities[rtype] + 1); // add one to be nonzero
		if (priority > highest_priority) {ret_rtype = rtype; highest_priority = priority;}
	}
	return breaker_zone_t(ret_rtype, room_start, room_end);
}
void building_t::toggle_circuit_breaker(bool is_on, unsigned zone_id, unsigned num_zones) {
	assert(has_room_geom());
	breaker_zone_t const zone(interior->get_circuit_breaker_info(zone_id, num_zones, get_window_vspace()));
	if (zone.invalid()) return; // no rooms for this zone

	if (zone.rtype == RTYPE_ELEVATOR) { // disable elevator; as long as we don't place breakers in elevators, the player can't get trapped in an elevator
		interior->elevators_disabled = !is_on;
		interior->room_geom->modified_by_player = 1;
		interior->room_geom->invalidate_lights_geom  ();
		interior->room_geom->update_dynamic_draw_data(); // needed for lit buttons
		auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
			if (i->type != TYPE_LIGHT || !i->in_elevator() || i->is_lit() == is_on) continue;
			i->flags ^= RO_FLAG_LIT; // toggle elevator light lit state
		}
		return;
	}
	for (unsigned r = zone.room_start; r < zone.room_end; ++r) {interior->rooms[r].unpowered = !is_on;} // update room unpowered flags; not yet used, but they may be later
	bool updated(0);
	auto objs_start(interior->room_geom->objs.begin());
	auto objs_end  (interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs_start; i != objs_end; ++i) {
		if ((i->is_powered() == is_on) || i->room_id < zone.room_start || i->room_id >= zone.room_end) continue; // no state change, or wrong zone

		if (i->is_light_type()) { // light
			if (i->in_elevator()) continue; // handled above
			bool const was_on(i->is_light_on());
			i->flags ^= RO_FLAG_NO_POWER; // need to disable this light and not allow the player/AI/motion detector to turn it back on
			bool const state_change(i->is_light_on() != was_on); // update if light on state changed
			updated  |= state_change;
			if (state_change) {register_indir_lighting_state_change(i - objs_start);} // Note: could be slow
		}
		else if (i->type == TYPE_MONITOR || i->type == TYPE_TV) { // interactive + drawn powered devices
			if (i->obj_id != 1) {interior->room_geom->invalidate_draw_data_for_obj(*i);}
			i->obj_id = 1; // turn it off
			i->flags ^= RO_FLAG_NO_POWER;
		}
		else if (i->type == TYPE_MWAVE || i->type == TYPE_CEIL_FAN || i->type == TYPE_CAMERA || i->type == TYPE_CLOCK || i->type == TYPE_LAVALAMP || i->type == TYPE_FISHTANK) {
			i->flags ^= RO_FLAG_NO_POWER; // interactive powered devices; stove is gas and not electric powered
		}
		// Note: stoves use gas rather than electricity and don't need power; lit exit signs are always on
	} // for i
	interior->room_geom->modified_by_player = 1; // I guess we need to set this, to be safe, as this breaker will likely have some effect
	if (!updated) return; // that's it, don't need to update geom
	// since the state of at least one light has changed, it's likely that other geom has been invalidated, so just update it all
	interior->room_geom->invalidate_lights_geom();
	interior->room_geom->invalidate_static_geom();
	interior->room_geom->invalidate_small_geom ();
	interior->room_geom->update_text_draw_data ();
}

// doors and other interactive objects

// used for drawing open doors
int building_t::find_ext_door_close_to_point(tquad_with_ix_t &door, point const &pos, float dist) const {
	if (doors.empty()) return -1;
	point query_pt(pos);
	if (is_rotated()) {do_xy_rotate_inv(bcube.get_cube_center(), query_pt);}
	int const room_id(get_room_containing_pt(query_pt));
	cube_t room_exp;

	if (room_id >= 0) { // pos is inside a room
		room_exp = get_room(room_id);
		room_exp.expand_by(get_wall_thickness()); // make sure it contains the door
	}
	// Note: returns the first exterior door found, assuming there can be at most one within dist of pos
	for (auto d = doors.begin(); d != doors.end(); ++d) {
		cube_t c(d->get_bcube());
		if (room_id >= 0 && !room_exp.contains_cube(c)) continue; // door not in the same room as pos - there is likely a wall between them
		c.expand_by_xy(dist);
		if (c.contains_pt(query_pt)) {door = *d; return (d - doors.begin());}
	} // for d
	return -1; // not found
}
bool building_t::point_near_ext_door(point const &pos, float dist) const { // simplified version of above function
	if (doors.empty()) return 0;
	point query_pt(pos);
	if (is_rotated()) {do_xy_rotate_inv(bcube.get_cube_center(), query_pt);}

	for (auto d = doors.begin(); d != doors.end(); ++d) {
		if (d->get_bcube().contains_pt_exp(query_pt, dist)) return 1;
	}
	return 0;
}

// used for pedestrians; pos should be outside the building
bool building_t::get_building_door_pos_closest_to(point const &target_pos, point &door_pos, bool inc_garage_door) const {
	float dmin_sq(0.0);

	for (auto d = doors.begin(); d != doors.end(); ++d) {
		if (d->type == tquad_with_ix_t::TYPE_GDOOR && !inc_garage_door) continue; // skip garage doors
		if (d->type == tquad_with_ix_t::TYPE_RDOOR) continue; // skip rooftop doors
		point const center(d->get_bcube().get_cube_center());
		float const dsq(p2p_dist_xy_sq(target_pos, center)); // ignore zval
		if (dmin_sq == 0.0 || dsq < dmin_sq) {door_pos = center; dmin_sq = dsq;}
	}
	if (dmin_sq == 0.0) return 0; // doors not added for some reason
	return 1;
}

void building_t::register_open_ext_door_state(int door_ix) {
	bool const is_open(door_ix >= 0), was_open(open_door_ix >= 0);
	bool const ring_doorbell(is_open && is_house && door_ix == 0 && city_action_key); // action key when front door of a house is open
	if (is_open == was_open && !ring_doorbell) return; // no state change
	unsigned const dix(is_open ? (unsigned)door_ix : (unsigned)open_door_ix);
	assert(dix < doors.size());
	auto const &door(doors[dix]);
	point const door_center(door.get_bcube().get_cube_center());

	if (ring_doorbell) {
		if (!camera_pdu.point_visible_test(door_center + get_camera_coord_space_xlate())) return; // not looking at the door
		gen_sound_thread_safe(SOUND_DOORBELL, local_to_camera_space(door_center)); // convert to camera space
		return;
	}
	static float last_sound_tfticks(0);

	if ((tfticks - last_sound_tfticks) > 0.25*TICKS_PER_SECOND) { // play at most once every 0.25 second
		last_sound_tfticks = tfticks;
		play_door_open_close_sound(door_center, is_open);
		vector3d const normal(door.get_norm());
		bool const dim(fabs(normal.x) < fabs(normal.y)), dir(normal[dim] < 0.0);
		point pos_interior(door_center);
		pos_interior[dim] += (dir ? 1.0 : -1.0)*CAMERA_RADIUS; // move point to the building interior so that it's a valid AI position
		register_building_sound(pos_interior, 0.4); // slightly quieter than interior doors because the user has no control over this
	}
	open_door_ix = door_ix;
}

bool check_obj_dir_dist(point const &closest_to, vector3d const &in_dir, cube_t const &c, point const &center, float dmax) { // Note: only zvals of door are used, no rotate required
	if (!c.closest_dist_less_than(closest_to, dmax)) return 0; // too far
	if (in_dir == zero_vector) return 1; // no direction filter specified
	point vis_pt(center.x, center.y, closest_to.z); // use query point zval
	max_eq(vis_pt.z, c.z1()); // clamp visibility test point to z-range of door to allow the player to open the door even looking at the top or bottom of it
	min_eq(vis_pt.z, c.z2());
	return (dot_product(in_dir, (vis_pt - closest_to).get_norm()) > 0.5); // door is not in the correct direction, skip
}

bool can_open_bathroom_stall_or_shower(room_object_t const &stall, point const &pos, vector3d const &from_dir) {
	if (stall.is_broken()) return 0; // broken, can't open
	// Note: since there currently aren't any objects the player can push/pull in bathrooms, we don't need to check if the door is blocked from opening; may need to revisit later
	bool dim(0), dir(0); // dim/dir that door is on
	if      (stall.type == TYPE_STALL ) {dim = stall.dim; dir = !stall.dir;} // bathroom stall
	else if (stall.type == TYPE_SHOWER) { // show stall
		dim = (stall.dx() < stall.dy());
		dir = !(dim ? stall.dir : stall.dim); // xdir=stall.dim, ydir=stall.dir
	}
	else {assert(0);}
	point door_center;
	door_center[ dim] = stall.d[dim][dir];
	door_center[!dim] = stall.get_center_dim(!dim);
	door_center.z = pos.z;
	return (dot_product_ptv(from_dir, door_center, pos) > 0.0); // facing the stall door
}

bool building_t::chair_can_be_rotated(room_object_t const &chair) const {
	if (chair.rotates()) return 1;
	// check if blocked by another object such as a desk
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_DESK && i->type != TYPE_TABLE) continue; // these should be the only objects a chair can be pushed under
		if (i->intersects(chair)) return 0;
	}
	return 1;
}

// called for the player; mode: 0=normal, 1=pull
bool building_t::apply_player_action_key(point const &closest_to_in, vector3d const &in_dir_in, int mode, bool check_only, bool no_check_conn_building) {
	if (!interior) return 0; // error?
	float const dmax(4.0*CAMERA_RADIUS), floor_spacing(get_window_vspace());
	float closest_dist_sq(0.0), t(0.0); // t is unused
	unsigned door_ix(0), obj_ix(0);
	bool found_item(0), is_obj(0);
	vector3d in_dir(in_dir_in);
	point closest_to(closest_to_in);
	maybe_inv_rotate_pos_dir(closest_to, in_dir);
	point const query_ray_end(closest_to + dmax*in_dir);

	if (mode == 0) { // if the player is in the closet, only the closet door can be opened
		for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
			if (player_in_closet && !i->is_closet_door())         continue; // only allow the player to open closet doors when in the closet
			float const door_z2(i->z2() + (i->on_stairs ? 0.25*floor_spacing : 0.0)); // increase height when on stairs; needed for steep basement stairs
			if (i->z1() > closest_to.z || door_z2 < closest_to.z) continue; // wrong floor, skip
			point const center(i->get_cube_center());
			float const dist_sq(p2p_dist_sq(closest_to, center));
			if (found_item && dist_sq >= closest_dist_sq) continue; // not the closest
			if (!check_obj_dir_dist(closest_to, in_dir, *i, center, (player_in_closet ? 0.5 : 1.0)*dmax)) continue; // door not in the correct direction or too far away
			cube_t const door_bcube(i->get_true_bcube()); // expand to nonzero area

			if (!door_bcube.line_intersects(closest_to, query_ray_end)) { // if camera ray doesn't intersect the door frame, check for ray intersection with opened door
				if (!i->open || (closest_to[i->dim] < i->d[i->dim][i->open_dir]) == i->open_dir) continue; // closed, or player not on the side the door opens to
				tquad_with_ix_t const door(set_interior_door_from_cube(*i));
				if (!line_poly_intersect(closest_to, query_ray_end, door.pts, door.npts, door.get_norm(), t)) continue; // test camera ray intersection with door plane
			}
			closest_dist_sq = dist_sq;
			door_ix    = (i - interior->doors.begin());
			found_item = 1;
		} // for i
	}
	if (has_room_geom()) { // check for closet doors in houses, bathroom stalls in office buildings, and other objects that can be interacted with
		vect_room_object_t &objs(interior->room_geom->objs), &expanded_objs(interior->room_geom->expanded_objs);
		auto objs_end(interior->room_geom->get_stairs_start());
		cube_t active_area;

		// make a first pass over all the large objects to determine if the player is inside one; in that case, the player can't reach out and interact with an object outside it
		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->type != TYPE_STALL && i->type != TYPE_SHOWER && i->type != TYPE_ELEVATOR && !(i->type == TYPE_CLOSET && i->is_open())) continue; // TYPE_CUBICLE?
			if (!i->contains_pt(closest_to)) continue;
			active_area = *i;
			break; // there can be only one - done
		}
		for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
			if (mode > 0) continue; // pull object only mode, skip this step
			auto const &obj_vect((vect_id == 1) ? expanded_objs : objs);
			unsigned const obj_id_offset((vect_id == 1) ? objs.size() : 0);
			auto obj_vect_end((vect_id == 1) ? expanded_objs.end() : objs.end()); // include all objects because blinds are added at the end

			for (auto i = obj_vect.begin(); i != obj_vect_end; ++i) {
				room_object const type(i->type);
				if (cur_player_building_loc.room_ix >= 0 && i->room_id != cur_player_building_loc.room_ix && type != TYPE_BUTTON) continue; // not in the same room as the player
				if (!active_area.is_all_zeros() && !i->intersects(active_area)) continue; // out of reach for the player
				// check for objects not in the attic when the player is in the attic and vice versa
				if (bool(player_in_attic) != i->in_attic() && type != TYPE_ATTIC_DOOR) continue;
				bool keep(0);
				if (type == TYPE_BOX && !i->is_open()) {keep = 1;} // box can only be opened once; check first so that selection works for boxes in closets
				else if (type == TYPE_CLOSET) {
					if (i->is_small_closet()) continue; // uses regular door now
					if (in_dir.z > 0.5)       continue; // not looking up at the light
					keep = 1; // closet door can be opened
				}
				else if (!player_in_closet) {
					if      ((type == TYPE_TOILET || type == TYPE_URINAL) && !i->is_broken()) {keep = 1;} // toilet/urinal can be flushed unless broken
					else if (type == TYPE_STALL && i->shape == SHAPE_CUBE && can_open_bathroom_stall_or_shower(*i, closest_to, in_dir)) {keep = 1;} // bathroom stall can be opened
					else if (type == TYPE_MIRROR && i->is_house())  {keep = 1;} // medicine cabinet
					else if (i->is_sink_type() || type == TYPE_TUB) {keep = 1;} // sink/tub
					else if (i->is_light_type() || type == TYPE_LAVALAMP) {keep = 1;} // room light or lamp
					else if (type == TYPE_FISHTANK && (i->flags & RO_FLAG_ADJ_TOP)) {keep = 1;} // fishtank with a lid and light
					else if (type == TYPE_PICTURE || type == TYPE_TPROLL || type == TYPE_MWAVE || type == TYPE_STOVE ||
						/*type == TYPE_FRIDGE ||*/ type == TYPE_TV || type == TYPE_MONITOR || type == TYPE_BLINDS || type == TYPE_SHOWER ||
						type == TYPE_SWITCH || type == TYPE_BOOK || type == TYPE_BRK_PANEL || type == TYPE_BREAKER || type == TYPE_ATTIC_DOOR ||
						type == TYPE_OFF_CHAIR || type == TYPE_WFOUNTAIN || type == TYPE_FALSE_DOOR || type == TYPE_LG_BALL) {keep = 1;}
					else if (type == TYPE_BUTTON && i->in_elevator() == bool(player_in_elevator)) {keep = 1;} // check for buttons inside/outside elevator
					else if (type == TYPE_PIZZA_BOX && !i->was_expanded()) {keep = 1;} // can't open if on a shelf
					else if (i->is_parked_car() && !i->is_broken()) {keep = 1;} // parked car with unbroken windows
					else if (!check_only && type == TYPE_SHELFRACK && !i->obj_expanded()) {keep = 1;} // expand shelfrack when action key is actually applied
					else if (type == TYPE_POOL_BALL && player_has_pool_cue()) {keep = 1;} // can only push pool ball if holding a pool cue
				}
				else if (type == TYPE_LIGHT) {keep = 1;} // closet light
				if (!keep) continue;
				cube_t obj_bc(*i), dishwasher;
				if (type == TYPE_KSINK && get_dishwasher_for_ksink(*i, dishwasher) && dishwasher.line_intersects(closest_to, query_ray_end)) {obj_bc = dishwasher;}
				else if (type == TYPE_KSINK || type == TYPE_BRSINK) {obj_bc = get_sink_cube(*i);} // the sink itself is actually smaller
				// shrink lamps in XY to a cube interior to their building cylinder to make drawers under lamps easier to select
				else if (type == TYPE_LAMP      ) {obj_bc.expand_by(vector3d(-i->dx(), -i->dy(), 0.0)*(0.5*(1.0 - 1.0/SQRT2)));}
				else if (type == TYPE_ATTIC_DOOR) {obj_bc = get_attic_access_door_cube(*i, 1);} // inc_ladder=1, to make it easier to select when in the attic
				point center;

				if (type == TYPE_CLOSET) {
					center = i->get_cube_center();
					center[i->dim] = i->d[i->dim][i->dir]; // use center of door, not center of closet
				}
				else {center = obj_bc.closest_pt(closest_to);}
				if (fabs(center.z - closest_to.z) > 0.7*floor_spacing) continue; // wrong floor
				// use dmax for closets and open breaker boxes to prioritize objects inside
				bool const low_priority(type == TYPE_CLOSET || (type == TYPE_BRK_PANEL && i->is_open()));
				float const dist_sq(low_priority ? dmax*dmax : p2p_dist_sq(closest_to, center));
				if (found_item && dist_sq >= closest_dist_sq)          continue; // not the closest
				if (!obj_bc.closest_dist_less_than(closest_to, dmax))  continue; // too far
				if (in_dir != zero_vector && !obj_bc.line_intersects(closest_to, query_ray_end)) continue; // player is not pointing at this object
				// checking for office chair rotation is expensive, so it's done last, just before updating closest
				if (type == TYPE_OFF_CHAIR && !chair_can_be_rotated(*i)) continue;
				closest_dist_sq = dist_sq;
				obj_ix = (i - obj_vect.begin()) + obj_id_offset;
				is_obj = found_item = 1;
			} // for i
		} // for vect_id
		if (!player_in_closet && mode == 0) {
			float const drawer_dist(found_item ? sqrt(closest_dist_sq) : 2.5*CAMERA_RADIUS);
			
			if (interior->room_geom->open_nearest_drawer(*this, closest_to, in_dir, drawer_dist, 0, check_only)) { // drawer is closer - open or close it
				return (check_only ? 1 : 0); // check_only returns 1 here because this counts as interactive
			}
		}
		if (!found_item && !check_only && !player_in_closet) {move_nearest_object(closest_to, in_dir, 3.0*CAMERA_RADIUS, mode);} // try to move an object instead
	}
	if (!found_item) { // no door or object found
		if (no_check_conn_building)   return 0; // avoid infinite recursion
		building_t *const conn_building(get_conn_bldg_for_pt(closest_to)); // check for door in connected building; should we check if room is ext conn first?
		if (conn_building == nullptr) return 0; // no other building
		return conn_building->apply_player_action_key(closest_to_in, in_dir_in, mode, check_only, 1); // chain call to the other building; no_check_conn_building=1
	}
	if (check_only) return 1;

	if (is_obj) { // interactive object
		if (!interact_with_object(obj_ix, closest_to, query_ray_end, in_dir)) return 0; // generate sound from the player height
	}
	else { // interior door
		door_t &door(interior->doors[door_ix]);
		if (!player_can_open_door(door)) return 0; // locked/blocked
		if (door.is_padlocked() && !door.open) {remove_padlock_from_door(door_ix, closest_to);}
		if (door.locked && !player_has_room_key()) {door.locked = 0;} // don't lock door when closing, to prevent the player from locking themselves in a room
		toggle_door_state(door_ix, 1, 1, closest_to); // toggle state if interior door; player_in_this_building=1, by_player=1, at player pos
		//interior->room_geom->modified_by_player = 1; // should door state always be preserved?
	}
	return 1;
}

bool building_t::interact_with_object(unsigned obj_ix, point const &int_pos, point const &query_ray_end, vector3d const &int_dir) {
	auto &obj(interior->room_geom->get_room_object_by_index(obj_ix));
	point const sound_origin(obj.xc(), obj.yc(), int_pos.z), local_center(local_to_camera_space(sound_origin)); // generate sound from the player height
	float sound_scale(0.5); // for building sound level
	bool update_draw_data(0);
	cube_t dishwasher;

	if (obj.type == TYPE_TOILET || obj.type == TYPE_URINAL) { // toilet/urinal can be flushed, but otherwise is not modified
		bool const is_urinal(obj.type == TYPE_URINAL); // urinal is quieter and higher pitch
		gen_sound_thread_safe(SOUND_FLUSH, local_center, (is_urinal ? 0.5 : 1.0), (is_urinal ? 1.25 : 1.0));
		sound_scale = 0.5;
		//refill_thirst(); // player can drink from toilet?
	}
	else if (obj.type == TYPE_KSINK && get_dishwasher_for_ksink(obj, dishwasher) && dishwasher.line_intersects(int_pos, query_ray_end)) { // dishwasher
		gen_sound_thread_safe_at_player(SOUND_METAL_DOOR, 0.2, 0.75);
		obj.flags       ^= RO_FLAG_OPEN; // toggle open/closed
		sound_scale      = 0.5;
		update_draw_data = 1;

		// since TYPE_KSINK already uses the RO_FLAG_EXPANDED flag for cabinet doors, we have to use the RO_FLAG_USED for dishwasher expansion
		if (obj.is_open() && !obj.is_used()) { // newly opened
			interior->room_geom->expand_dishwasher(obj, dishwasher);
			obj.flags |= RO_FLAG_USED; // can't expand again
		}
		else if (!obj.is_open() && obj.is_used()) { // closed
			interior->room_geom->unexpand_dishwasher(obj, dishwasher);
			obj.flags &= ~RO_FLAG_USED; // can now expand again
		}
	}
	else if (obj.is_sink_type() || obj.type == TYPE_TUB) { // sink or tub
		if (!obj.is_active() && obj.type == TYPE_TUB) {
			gen_sound_thread_safe(SOUND_SINK, local_center); // play sound when turning the tub on
			
			if (obj.item_flags < 4) { // water level is 0-4
				++obj.item_flags;
				interior->room_geom->invalidate_static_geom();
			}
			//refill_thirst(); // player can drink from tub?
		}
		if (obj.is_sink_type()) {
			obj.flags ^= RO_FLAG_IS_ACTIVE; // toggle active bit, only for sinks for now

			if (obj.is_active() && obj.item_flags == 0) { // no water yet
				obj.item_flags ^= 1; // mark as filled with water
				interior->room_geom->invalidate_static_geom();
			}
			refill_thirst(); // player can drink from sink
		}
		sound_scale = 0.4;
	}
	else if (obj.is_light_type()) {
		toggle_light_object(obj, obj.get_cube_center());
		sound_scale = 0.0; // sound has already been registered above
	}
	else if (obj.type == TYPE_TPROLL) {
		if (!obj.is_hanging() && !obj.was_expanded()) {
			gen_sound_thread_safe(SOUND_FOOTSTEP, local_center, 0.5, 1.5); // could be better
			obj.flags |= RO_FLAG_HANGING; // pull down the roll
			update_draw_data = 1;
		}
		sound_scale = 0.0; // no sound
	}
	else if (obj.type == TYPE_PICTURE) { // tilt the picture
		obj.flags |= RO_FLAG_RAND_ROT;
		++obj.item_flags; // choose a different random rotation
		gen_sound_thread_safe(SOUND_SLIDING, local_center, 0.25, 2.0); // higher pitch
		sound_scale      = 0.0; // no sound
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_OFF_CHAIR) { // handle rotate of office chair
		office_chair_rot_rate += 0.1;
		obj.flags |= RO_FLAG_ROTATING; // Note: this is a model, no need to regen vertex data
		gen_sound_thread_safe(SOUND_SQUEAK, local_center, 0.25, 0.5); // lower pitch
		sound_scale = 0.2;
	}
	else if (obj.type == TYPE_MWAVE) {
		cube_t const panel(get_mwave_panel_bcube(obj));
		float cur_tmin(0.0), cur_tmax(1.0);

		if (!get_line_clip(int_pos, query_ray_end, panel.d, cur_tmin, cur_tmax)) { // not pointing at the panel - open and close the door
			obj.flags       ^= RO_FLAG_OPEN; // toggle open/closed
			update_draw_data = 1;
			gen_sound_thread_safe((obj.is_open() ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_DOOR_CLOSE), local_center, 0.5, 1.6);
		}
		else if (obj.is_powered()) { // pointing at the panel - make it beep
			gen_sound_thread_safe(SOUND_BEEP, local_center, 0.25);
			sound_scale = 0.6;
		}
	}
	else if (obj.type == TYPE_STOVE) { // toggle burners; doesn't need power
		float const height(obj.dz());
		bool const dim(obj.dim), dir(obj.dir);
		unsigned burner_id(0);
		float tmin(1.0);
		
		for (unsigned w = 0; w < 2; ++w) { // width dim
			for (unsigned d = 0; d < 2; ++d) { // depth dim
				bool wv(bool(w) ^ dim ^ dir ^ 1), dv(bool(d) ^ dir);
				cube_t c(obj);
				set_cube_zvals(c, (c.z1() + 0.7*height), (c.z1() + 0.8*height)); // select the cook top area
				((dim ? wv : dv) ? c.x1() : c.x2()) = c.xc();
				((dim ? dv : wv) ? c.y1() : c.y2()) = c.yc();
				float cur_tmin(0.0), cur_tmax(1.0);
				if (!get_line_clip(int_pos, query_ray_end, c.d, cur_tmin, cur_tmax) || tmin < cur_tmin) continue;
				tmin      = cur_tmin;
				burner_id = (2U*w + d); // this point intersects earlier - select this burner
			} // for d
		} // for w
		if (tmin < 1.0) { // found an intersection
			unsigned const flag_mask(1U<<burner_id);
			bool const is_on(obj.item_flags & flag_mask);
		
			if (is_on) { // currently on, turn off
				obj.item_flags &= ~flag_mask;
				gen_sound_thread_safe(SOUND_CLICK, local_center, 0.75);
				sound_scale = 0.1;
			}
			else { // currently off, turn on
				obj.item_flags |= flag_mask;
				gen_sound_thread_safe(SOUND_HISS, local_center, 0.25, 0.6);
				sound_scale = 0.3;
			}
		}
	}
	else if (obj.type == TYPE_FRIDGE) {
		obj.flags       ^= RO_FLAG_OPEN; // toggle open/closed
		update_draw_data = 1;
		gen_sound_thread_safe((obj.is_open() ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_DOOR_CLOSE), local_center, 0.4, 0.67);
	}
	else if (obj.type == TYPE_TV || obj.type == TYPE_MONITOR) {
		if (obj.is_powered()) {
			if (!obj.is_broken()) { // no visual effect if broken, but still clicks
				if (obj.type == TYPE_MONITOR && (obj.obj_id & 1)) {--obj.obj_id;} // toggle on and off, but don't change the desktop
				else {++obj.obj_id;} // toggle on/off, and also change the picture
				update_draw_data = 1;
			}
			gen_sound_thread_safe(SOUND_CLICK, local_center, 0.4);
		}
	}
	else if (obj.type == TYPE_BUTTON) { // Note: currently, buttons are only used for elevators
		if (!obj.is_active() && !interior->elevators_disabled) { // if not already active
			register_button_event(obj);
			obj.flags |= RO_FLAG_IS_ACTIVE;
			interior->room_geom->invalidate_draw_data_for_obj(obj); // need to regen object data due to lit state change; don't have to set modified_by_player
		}
	}
	else if (obj.type == TYPE_SWITCH) {
		// should select the correct light(s) for the room containing the switch
		toggle_room_light(obj.get_cube_center(), 1, obj.room_id, 0, obj.in_closet(), obj.in_attic()); // exclude lamps; select closet lights if a closet light switch
		gen_sound_thread_safe_at_player(SOUND_CLICK, 0.5);
		obj.flags       ^= RO_FLAG_OPEN; // toggle on/off
		sound_scale      = 0.1;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_BREAKER) {
		gen_sound_thread_safe_at_player(SOUND_CLICK, 1.0);
		obj.flags       ^= RO_FLAG_OPEN; // toggle on/off
		toggle_circuit_breaker(obj.is_open(), obj.obj_id, obj.item_flags);
		sound_scale      = 0.25;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_BLINDS) { // see building_t::add_window_blinds()
		if (!adjust_blinds_state(obj_ix)) return 0;
		gen_sound_thread_safe_at_player(SOUND_SLIDING, 0.5);
		sound_scale      = 0.3;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_BOOK) {
		if (!obj.is_open()) { // check if there's space to open the book
			room_object_t open_area(obj); // the area that must be clear if we want to open this book
			open_area.z1() += 0.75*obj.dz(); // shift the bottom up to keep it from intersecting whatever this book is resting on
			open_area.translate_dim(obj.dim, 1.01*(obj.dir ? -1.0 : 1.0)*obj.get_sz_dim(obj.dim)); // translate more than width to keep it from overlapping obj
			if (!is_obj_pos_valid(open_area, 0, 1, 0)) return 0; // intersects some part of the building; keep_in_room=0, allow_block_door=1, check_stairs=0
			if (overlaps_any_placed_obj(open_area))    return 0;
		}
		if (!check_for_water_splash(sound_origin, 0.4, 1)) {gen_sound_thread_safe_at_player(SOUND_OBJ_FALL, 0.25);} // splash or drop; full_room_height=1
		obj.flags       ^= RO_FLAG_OPEN; // toggle open/closed
		sound_scale      = 0.1; // very little sound
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_SHOWER) { // shower
		// if (interior->room_geom->cube_intersects_moved_obj(c_test)) continue; // not yet needed
		if (can_open_bathroom_stall_or_shower(obj, int_pos, int_dir)) { // open/close shower door
			obj.flags       ^= RO_FLAG_OPEN; // toggle open/close
			sound_scale      = 0.35;
			update_draw_data = 1;
			play_open_close_sound(obj, sound_origin);
		}
		else { // turn on shower water
			gen_sound_thread_safe_at_player(SOUND_SINK);
			sound_scale = 0.5;
			if (!obj.is_open()) {register_achievement("Squeaky Clean");}

			if (!obj.item_flags) {
				obj.item_flags = 1; // mark as filled with water
				interior->room_geom->invalidate_static_geom();
			}
		}
	}
	else if (obj.type == TYPE_BOX) {
		if (!check_for_water_splash(sound_origin, 0.6)) {gen_sound_thread_safe_at_player(SOUND_OBJ_FALL, 0.5);}
		obj.flags       |= RO_FLAG_OPEN; // mark as open
		sound_scale      = 0.2;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_PIZZA_BOX) {
		gen_sound_thread_safe_at_player(SOUND_SLIDING, 0.1);
		obj.flags       ^= RO_FLAG_OPEN; // toggle open/close
		sound_scale      = 0.0; // no sound
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_CLOSET || obj.type == TYPE_STALL) {
		if (!obj.is_open()) { // not yet open
			// remove any spraypaint or marker that's on the door; would be better if we could move it with the door, or add it back when the door is closed
			cube_t door(get_open_closet_door(obj));
			door.expand_in_dim(obj.dim, get_wall_thickness());
			remove_paint_in_cube(door); // use the door before it's opened
		}
		obj.flags ^= RO_FLAG_OPEN; // toggle open/close

		if (obj.type == TYPE_CLOSET) {
			interior->room_geom->expand_object(obj, *this); // expand any boxes so that the player can pick them up
			sound_scale = 0.25; // closets are quieter, to allow players to more easily hide
		}
		play_open_close_sound(obj, sound_origin);
		register_indir_lighting_geom_change();
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_MIRROR && obj.is_house()) { // medicine cabinet
		obj.flags ^= RO_FLAG_OPEN; // toggle open/close
		interior->room_geom->expand_object(obj, *this);
		play_open_close_sound(obj, sound_origin);
		sound_scale      = 0.4;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_BRK_PANEL) { // breaker panel
		obj.flags ^= RO_FLAG_OPEN; // toggle open/close
		interior->room_geom->expand_object(obj, *this);
		play_open_close_sound(obj, sound_origin);
		sound_scale      = 0.6;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_ATTIC_DOOR) {
		gen_sound_thread_safe_at_player(SOUND_SLIDING, 1.0); // better sound?
		obj.flags       ^= RO_FLAG_OPEN; // open/close
		sound_scale      = 0.5;
		update_draw_data = 1;
		interior->attic_access_open ^= 1;
		// toggle the attic light as well
		auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
			if (i->type != TYPE_LIGHT || !i->in_attic()) continue; // not attic light
			if (i->is_lit() != obj.is_open()) {toggle_light_object(*i, sound_origin);}
		}
	}
	else if (obj.type == TYPE_LAVALAMP) {
		obj.flags       ^= RO_FLAG_LIT; // toggle lit
		update_draw_data = 1;
		gen_sound_thread_safe(SOUND_CLICK, local_center, 0.35);
	}
	else if (obj.type == TYPE_FISHTANK) {
		obj.flags       ^= RO_FLAG_LIT; // toggle the light on the lid; no draw data update
		gen_sound_thread_safe(SOUND_CLICK, local_center, 0.4);
	}
	else if (obj.type == TYPE_WFOUNTAIN) {
		refill_thirst();
		gen_sound_thread_safe(SOUND_GULP, local_center, 0.5);
		sound_scale = 0.1; // very little sound
	}
	else if (obj.type == TYPE_FALSE_DOOR) { // locked, can't open
		print_text_onscreen("Door is locked", RED, 1.0, 2.0*TICKS_PER_SECOND, 0);
		gen_sound_thread_safe_at_player(SOUND_CLICK, 1.0, 0.6);
		return 0;
	}
	else if (obj.type == TYPE_SHELFRACK) { // expand shelfrack
		interior->room_geom->expand_object(obj, *this);
	}
	else if (is_ball_type(obj.type)) { // push the ball
		room_obj_dstate_t &dstate(interior->room_geom->get_dstate(obj));
		dstate.velocity.x += 0.5*KICK_VELOCITY*int_dir.x;
		dstate.velocity.y += 0.5*KICK_VELOCITY*int_dir.y;

		if (!obj.is_dynamic()) { // make it dynamic
			obj.flags |= RO_FLAG_DYNAMIC;
			interior->room_geom->invalidate_small_geom();
		}
	}
	else if (obj.is_parked_car()) {
		gen_sound_thread_safe_at_player(SOUND_GLASS);
		register_broken_object(obj);
		add_broken_glass_to_floor(int_pos, 0.8*CAMERA_RADIUS);
		assert(!obj.is_broken());
		obj.flags  |= RO_FLAG_BROKEN;
		sound_scale = 1.0; // loud sound, but no update of draw data
	}
	else {assert(0);} // unhandled type
	if (update_draw_data) {interior->room_geom->update_draw_state_for_room_object(obj, *this, 0);}
	if (sound_scale > 0.0) {register_building_sound(sound_origin, sound_scale);}
	if (obj.type == TYPE_BOX) {add_box_contents(obj);} // must be done last to avoid reference invalidation
	return 1;
}

bool building_t::adjust_blinds_state(unsigned obj_ix) {
	auto &obj(interior->room_geom->get_room_object_by_index(obj_ix));

	if (obj.is_hanging()) { // hanging horizontal blinds
		float const floor_spacing(get_window_vspace()), window_v_border(0.94*get_window_v_border()); // border_mult=0.94 to account for the frame
		float const window_height(floor_spacing*(1.0 - 2.0*window_v_border)), blinds_height(obj.dz());
		assert(window_height > 0.0);
		bool const mostly_open(blinds_height < 0.5*window_height);
		if (mostly_open) {obj.z1() = obj.z2() - floor_spacing*(1.0 - window_v_border) + 0.05*floor_spacing;} // close the blinds fully
		else             {obj.z1() = obj.z2() - (1.8*get_wall_thickness() + 0.05*floor_spacing);} // open the blinds
		// set new thickness (matches building_t::add_window_blinds())
		obj.d[obj.dim][!obj.dir] = obj.d[obj.dim][obj.dir] + (mostly_open ? 0.0927 : 0.227)*get_wall_thickness()*(obj.dir ? -1.0 : 1.0);
	}
	else { // vertical blinds - in pairs
		assert(obj.flags & (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI)); // should have had one of these flags set
		bool const move_dir(obj.flags & RO_FLAG_ADJ_HI);
		unsigned other_blinds_ix(0);

		if (move_dir) { // this is the left side blind, the other side is to the right in the next slot
			other_blinds_ix = obj_ix + 1;
		}
		else { // this is the right side blind, the other side is to the left in the previous slot
			assert(obj_ix > 0);
			other_blinds_ix = obj_ix - 1;
		}
		auto &other_blinds(interior->room_geom->get_room_object_by_index(other_blinds_ix));
		if (other_blinds.type != TYPE_BLINDS) {assert(0); return 0;} // was taken, etc.
		float const fixed_end(obj.d[!obj.dim][!move_dir]), width(obj.get_width());
		float const window_center(0.5f*(fixed_end + other_blinds.d[!obj.dim][move_dir])); // center of the span of the pair of left/right blinds
		bool const mostly_open(width < 0.5*fabs(fixed_end - window_center));
		float &move_edge(obj.d[!obj.dim][move_dir]);
		if (mostly_open) {move_edge = window_center;} // close the blinds fully
		else             {move_edge = 0.25*window_center + 0.75*fixed_end;} // open the blinds
	}
	assert(obj.is_strictly_normalized());
	register_blinds_state_change();
	return 1;
}

void building_t::toggle_door_state(unsigned door_ix, bool player_in_this_building, bool by_player, point const &actor_pos) { // called by the player or AI
	assert(interior);
	door_t &door(interior->get_door(door_ix));
	door.toggle_open_state(/*by_player*/player_in_this_building); // allow partial open/animated door if player is in this building
	// we changed the door state, but navigation should adapt to this, except for doors on stairs (which are special)
	if ( door.on_stairs ) {invalidate_nav_graph();} // any in-progress paths may have people walking to and stopping at closed/locked doors
	if (!door.for_closet) {interior->door_state_updated = 1;} // required for AI navigation logic to adjust to this change; what about backrooms doors?
	if (has_room_geom() ) {interior->room_geom->invalidate_mats_mask |= (1 << MAT_TYPE_DOORS);} // need to recreate doors VBO
	check_for_water_splash(cube_bot_center(door), 2.0); // big splash

	if (player_in_this_building || by_player) { // is it really safe to call this from the AI thread?
		point door_center(door.xc(), door.yc(), actor_pos.z);
		// if door was opened or is fully closed play a sound; otherwise, play the close sound later when fully closed
		if (door.open || door.open_amt == 0.0) {play_door_open_close_sound(door_center, door.open);}
		
		if (by_player) {
			// bias the sound slightly toward the side of the door the player is on to force a zombie to go through the doorway to get there,
			// rather than targeting the exact center of the doorway and possibly clipping through the wall (for example, if in a hallway)
			door_center[door.dim] += get_wall_thickness()*((door_center[door.dim] < actor_pos[door.dim]) ? 1.0 : -1.0);
			register_building_sound(door_center, 0.5);
		}
		// update indir lighting state if needed; for now this is only for player actions to avoid thread safety issues and too many updates
		if (by_player && enable_building_indir_lighting()) {
			// Note: only have to register geom change if the light is on the same floor as the player, but this should always be true if the player just closed this door
			register_indir_lighting_geom_change();
			static vector<unsigned> light_ids;
			get_lights_near_door(door_ix, light_ids);
			for (unsigned light_ix : light_ids) {register_indir_lighting_state_change(light_ix, 1);} // is_door_change=1
		}
	}
	handle_items_intersecting_closed_door(door); // check if we need to move any objects out of the way
	
	if (door.open) { // was closed and now open
		notify_door_fully_closed_state(door);
		cube_t door_exp(door);
		door_exp.expand_in_dim(door.dim, 0.6*get_wall_thickness()); // make sure decals are included
		remove_paint_in_cube(door_exp); // remove any paint that was over the closed door
	}
}
void building_t::notify_door_fully_closed_state(door_t const &door) {
	if (door.obj_ix < 0)  return; // no associated object
	if (!has_room_geom()) return; // error?
	assert((unsigned)door.obj_ix < interior->room_geom->objs.size());
	room_object_t &obj(interior->room_geom->objs[door.obj_ix]);
	
	if (obj.type == TYPE_CLOSET) { // this was a closet door
		if (door.open) {obj.flags |= RO_FLAG_OPEN;} else {obj.flags &= ~RO_FLAG_OPEN;}
		if (door.open) {interior->room_geom->expand_object(obj, *this);} // expand any boxes so that the player can pick them up
	}
}
void building_t::handle_items_intersecting_closed_door(door_t const &door) {
	if (door.open || door.open_amt > 0.0) return; // not fully closed or about to open
	if (!has_room_geom()) return;
	cube_t const door_bcube(door.get_true_bcube());

	for (room_object_t &obj : interior->room_geom->expanded_objs) { // currently only expanded objects such as books
		if (!obj.intersects(door_bcube)) continue;
		float const door_center(door_bcube.get_center_dim(door.dim)), obj_center(obj.get_center_dim(door.dim));
		bool const move_dir(door_center < obj_center);
		float const move_dist(door_bcube.d[door.dim][move_dir] - obj.d[door.dim][!move_dir]);
		obj.translate_dim(door.dim, move_dist);
		interior->room_geom->invalidate_draw_data_for_obj(obj);
	} // for obj
}

bool door_t::check_key_mask_unlocks(unsigned key_mask) const {
	if (!is_closed_and_locked()) return 1;
	if (key_mask == 0)           return 0; // no key
	if (is_locked_unlockable())  return 0;
	if (!is_padlocked())         return 1; // any key color opens a non-padlocked door
	return (key_mask & (1 << get_padlock_color_ix())); // only a key matching the padlock color unlocks a padlocked door
}
void door_t::toggle_open_state(bool allow_partial_open) {
	open ^= 1;
	if (!allow_partial_open || on_stairs) {make_fully_open_or_closed();} // update open_amt immediately if needed
}
bool door_t::next_frame() { // returns true if state changed
	if (!is_partially_open()) return 0;
	open_amt += (open ? 1.0 : -(auto_close ? 0.1 : 1.0))*4.0*(fticks/TICKS_PER_SECOND); // 0.25s to open/close; auto closing doors close slowly
	open_amt  = CLIP_TO_01(open_amt);
	return 1;
}
void building_t::doors_next_frame(point const &player_pos) {
	if (!has_room_geom()) return;
	interior->last_active_door_ix = -1;

	for (auto d = interior->doors.begin(); d != interior->doors.end(); ++d) {
		if (d->auto_close && d->open_amt > 0.0) { // handle auto closing
			cube_t open_area(*d);
			open_area.expand_in_dim(d->dim, 0.5*d->get_width()); // expand halfway (for the non-opening side)
			open_area.union_with_cube(d->get_open_door_path_bcube()); // include the path of the door
			bool is_blocked(open_area.contains_pt_exp_xy_only(player_pos, get_scaled_player_radius())); // check if blocked by player
			
			if (!is_blocked) { // check building AI people blocking the door
				for (person_t const &person : interior->people) {
					if (open_area.contains_pt_exp_xy_only(person.pos, person.radius)) {is_blocked = 1; break;}
				}
			}
			if (is_blocked) {d->open = 1;} // push open
			else if (d->open_amt == 1.0) { // auto close
				toggle_door_state((d - interior->doors.begin()), 1, 1, d->get_cube_center());
			}
		}
		if (!d->next_frame()) continue;
		handle_items_intersecting_closed_door(*d);
		
		if (!d->open && d->open_amt == 0.0) { // door closes fully
			play_door_open_close_sound(point(d->xc(), d->yc(), camera_pos.z), 0); // play close sound at player z; open=0
			notify_door_fully_closed_state(*d);
		}
		interior->last_active_door_ix = (d - interior->doors.begin());
	} // for d
	if (interior->last_active_door_ix >= 0) {interior->room_geom->invalidate_mats_mask |= (1 << MAT_TYPE_DOORS);} // need to recreate doors VBO
}

point building_t::local_to_camera_space(point const &pos) const {
	point pos_rot(pos);
	if (is_rotated()) {do_xy_rotate(bcube.get_cube_center(), pos_rot);}
	return (pos_rot + get_camera_coord_space_xlate()); // convert to camera space
}

void building_t::play_door_open_close_sound(point const &pos, bool open, float gain, float pitch) const {
	gen_sound_thread_safe((open ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_DOOR_CLOSE), local_to_camera_space(pos), gain, pitch);
}
// called for both player and AI actions
void building_t::play_open_close_sound(room_object_t const &obj, point const &sound_origin) const {
	if (obj.is_small_closet() || obj.type == TYPE_STALL) {
		play_door_open_close_sound(sound_origin, obj.is_open(), 1.0, ((obj.type == TYPE_STALL) ? 2.0 : 1.0)); // higher pitch for stalls
	}
	else if (obj.type == TYPE_CLOSET   ) {gen_sound_thread_safe_at_player(SOUND_SLIDING);}
	else if (obj.type == TYPE_SHOWER   ) {gen_sound_thread_safe_at_player((obj.is_open() ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_METAL_DOOR));}
	else if (obj.type == TYPE_MIRROR   ) {play_door_open_close_sound(sound_origin, obj.is_open(), 0.4, 1.6);} // medicine cabinet
	else if (obj.type == TYPE_BRK_PANEL) {play_door_open_close_sound(sound_origin, obj.is_open(), 0.6, 1.8);} // breaker panel
	else {assert(0);} // not implemented
}

// dynamic objects: elevators and balls

bool apply_object_bounce(vector3d &velocity, vector3d const &cnorm, float hardness, bool on_floor) {
	float const vmag(velocity.mag()), elasticity(OBJ_ELASTICITY*hardness);
	if (vmag < TOLERANCE) return 0;
	vector3d v_ref;
	calc_reflection_angle(velocity/vmag, v_ref, cnorm);
	velocity = vmag*v_ref;
	if (on_floor && cnorm == plus_z) {velocity.z *= elasticity;} else {velocity *= elasticity;} // only attenuate velocity in Z for a floor collision
	return 1;
}
void apply_object_bounce_with_sound(building_t const &building, vector3d &velocity, vector3d const &cnorm, point const &pos, float hardness, bool on_floor) {
	if (!apply_object_bounce(velocity, cnorm, hardness, on_floor)) return;
	float const bounce_volume(min(1.0f, ((cnorm == plus_z) ? 0.5f : 1.0f)*velocity.mag()/KICK_VELOCITY)); // relative to kick velocity; lower if vertical/floor coll

	if (bounce_volume > 0.25) { // apply bounce sound
		if (bounce_volume > 0.5) {
			if (bounce_volume > 0.7 && building.check_for_water_splash(pos, 0.75*bounce_volume, 0, 1)) return; // handled as splash; full_room_height=0, draw_splash=1
			gen_sound_thread_safe(SOUND_KICK_BALL, building.local_to_camera_space(pos), 0.75*bounce_volume*bounce_volume);
		}
		register_building_sound(pos, 0.7*bounce_volume);
	}
}
void apply_floor_vel_thresh(vector3d &velocity, vector3d const &cnorm) {
	if (cnorm == plus_z) { // collision with the floor or the top surface of something
		if (fabs(velocity.z) < 0.25*OBJ_GRAVITY*fticks) {velocity.z = 0.0;} // zero velocity z component if near zero to reduce instability
	}
}

bool check_ball_kick(room_object_t &ball, vector3d &velocity, point &new_center, point const &p_pos, float pz1, float pz2, float pradius) {
	if (ball.z2() < pz1 || ball.z1() > pz2) return 0; // no collision
	point const center(ball.get_cube_center()), ppos(p_pos.x, p_pos.y, center.z); // use zval of object (Z range was checked earlier)
	float const r_sum(ball.get_radius() + pradius);
	if (!dist_xy_less_than(ppos, center, r_sum)) return 0; // no collision
	vector3d const dir((center - ppos).get_norm());
	new_center = (center + (1.05*r_sum - p2p_dist_xy(ppos, center))*dir); // move so that it no longer collides with a bit of tolerance
	velocity.x = KICK_VELOCITY*dir.x; velocity.y = KICK_VELOCITY*dir.y; // keep existing velocity.z
	return 1;
}

void apply_building_gravity(float &vz, float dt_ticks) {
	vz -= OBJ_GRAVITY*dt_ticks; // apply gravitational acceleration
	max_eq(vz, -TERM_VELOCITY);
}

void obj_dynamic_to_static(room_object_t &obj, building_interior_t &interior) {
	obj.flags &= ~RO_FLAG_DYNAMIC; // clear dynamic flag
	interior.update_dynamic_draw_data(); // remove from dynamic objects and schedule an update
	interior.room_geom->invalidate_draw_data_for_obj(obj); // add to small static objects
}

void building_t::run_ball_update(vector<room_object_t>::iterator ball_it, point const &player_pos, float player_z1, bool player_is_moving) {
	room_object_t &ball(*ball_it);

	if (ball.type == TYPE_POOL_BALL && ball.is_dynamic()) { // check of moving on the pool table
		assert(ball.state_flags < interior->room_geom->objs.size());
		room_object_t const &pool_table(interior->room_geom->objs[ball.state_flags]);

		if (pool_table.type == TYPE_POOL_TABLE && pool_table.contains_pt(cube_bot_center(ball))) { // ball on pool table
			update_pool_table(ball);
			return;
		}
	}
	assert(is_ball_type(ball.type)); // currently, only balls have has_dstate()
	float const player_radius(get_scaled_player_radius()), player_z2(player_pos.z), radius(ball.get_radius());
	float const fc_thick(get_fc_thickness()), fticks_stable(min(fticks, 1.0f)); // cap to 1/40s to improve stability
	room_obj_dstate_t &dstate(interior->room_geom->get_dstate(ball));
	vector3d &velocity(dstate.velocity);
	ball_type_t const &bt(ball.get_ball_type());
	point const center(ball.get_cube_center());
	bool const was_dynamic(ball.is_dynamic()), is_moving_fast(velocity.mag() > 0.5*KICK_VELOCITY), can_kick(bt.can_kick && !ball.no_coll());
	bool on_floor(0), kicked(0);
	point new_center(center);

	// check the player, but not if they're looking directly at the ball; assume in that case they intend to pick it up instead
	if (can_kick && camera_surf_collide && player_is_moving && !player_wait_respawn && dot_product_ptv(cview_dir, new_center, player_pos) < 0.9*p2p_dist(new_center, player_pos)) {
		kicked |= check_ball_kick(ball, velocity, new_center, player_pos, player_z1, player_z2, player_radius);
	}
	for (auto p = interior->people.begin(); p != interior->people.end(); ++p) { // check building AI people
		if (is_moving_fast) { // treat collision as a bounce
			vector3d cnorm;

			if (sphere_cube_int_update_pos(new_center, radius, p->get_bcube(), center, 0, &cnorm)) {
				register_person_hit((p - interior->people.begin()), ball, velocity);
				apply_object_bounce_with_sound(*this, velocity, cnorm, new_center, 0.75, on_floor); // hardness=0.75
			}
		}
		else if (can_kick) { // treat collision as a kick
			kicked |= check_ball_kick(ball, velocity, new_center, p->pos, p->get_z1(), p->get_z2(), 0.6*p->get_width());
		}
	} // for p
	if (kicked) {
		static float last_sound_tfticks(0);
		static point last_sound_pt(all_zeros);
		ball.flags |= RO_FLAG_DYNAMIC; // make it dynamic

		if ((tfticks - last_sound_tfticks) > 1.0*TICKS_PER_SECOND && !dist_less_than(new_center, last_sound_pt, radius)) { // play at most once per second
			check_for_water_splash(new_center, 1.5); // check for splash, but also play kick sound; will be extra loud
			gen_sound_thread_safe(SOUND_KICK_BALL, local_to_camera_space(new_center), 0.5);
			register_building_sound(new_center, 0.75);
			last_sound_tfticks = tfticks;
			last_sound_pt      = new_center;
		}
	}
	else if (was_dynamic) { // not colliding, but is moving
		float const max_timestep(0.1); // in ticks (40 per second)
		unsigned const num_steps(max(1, round_fp(fticks_stable/max_timestep)));
		float const step_sz(fticks_stable/num_steps);

		for (unsigned step = 0; step < num_steps; ++step) {
			point const test_pt(new_center.x, new_center.y, (new_center.z - radius - 0.1*fc_thick));
			on_floor = 0; // reset for this iteration

			for (auto f = interior->floors.begin(); f != interior->floors.end(); ++f) {
				if (f->contains_pt(test_pt)) {on_floor = 1; break;}
			}
			if (on_floor) { // moving on the floor, apply surface friction
				velocity *= (1.0f - min(1.0f, bt.friction*OBJ_DECELERATE*step_sz));
				if (velocity.mag_sq() < MIN_VELOCITY*MIN_VELOCITY) {velocity = zero_vector;} // zero velocity if stopped
			}
			else { // in the air - apply gravity
				apply_building_gravity(velocity.z, step_sz);
			}
			if (velocity == zero_vector) { // stopped
				obj_dynamic_to_static(ball, *interior);
				break; // done
			}
			new_center += velocity*step_sz; // move based on velocity
		} // for step
	}
	if (new_center == center) return; // not moving, done
	// check for collisions and move to new location
	vector3d const delta(new_center - center);
	unsigned const num_steps(max(1U, unsigned(delta.mag()/radius))); // don't move more than ball radius
	float const step_len(1.0/num_steps);
	vector3d cnorm;
	int obj_ix(-1);
	float hardness(0.0);
	bool had_coll(0);

	for (unsigned n = 0; n < num_steps && !had_coll; ++n) {
		new_center = center + ((n+1)*step_len)*delta;
		had_coll   = interior->check_sphere_coll(*this, new_center, center, radius, ball_it, cnorm, hardness, obj_ix, 1); // is_ball=1
	}
	if (had_coll) {
		hardness *= bt.elastic;
		apply_floor_vel_thresh(velocity, cnorm);
		apply_object_bounce_with_sound(*this, velocity, cnorm, new_center, hardness, on_floor);

		if (obj_ix >= 0) { // collided with a room object
			auto &obj(interior->room_geom->get_room_object_by_index(obj_ix));
			bool handled(0);

			// break the glass if not already broken; should windows get broken as well?
			if (bt.breaks_glass && (obj.type == TYPE_TV || obj.type == TYPE_MONITOR || obj.type == TYPE_DRESS_MIR || (obj.type == TYPE_MIRROR && !obj.is_open())) &&
				velocity.mag() > 2.0*MIN_VELOCITY && !obj.is_broken())
			{
				vector3d front_dir(all_zeros);
				front_dir[obj.dim] = (obj.dir ? 1.0 : -1.0);

				if (dot_product(cnorm, front_dir) > 0.9) { // hit the front side of the screen
					if (dist_less_than(new_center, obj.get_cube_center(), (radius + 0.5*obj.get_length() + 0.2*obj.dz()))) { // near the screen center
						// capture value before breaking; if the player then takes this object, damage will be higher, but we can attribute this to making a mess of broken glass
						register_broken_object(obj);
						obj.flags |= RO_FLAG_BROKEN;
						point const sound_origin(obj.xc(), obj.yc(), new_center.z); // generate sound from the player height
						gen_sound_thread_safe(SOUND_GLASS, local_to_camera_space(sound_origin), 0.7);
						register_building_sound(sound_origin, 0.7);
						interior->room_geom->update_draw_state_for_room_object(obj, *this, 0);
						if (obj.type == TYPE_DRESS_MIR || obj.type == TYPE_MIRROR) {register_achievement("7 Years of Bad Luck");}
						else if ((obj.type == TYPE_TV || obj.type == TYPE_MONITOR) && obj.is_powered()/*!(obj.obj_id & 1)*/) { // only if turned on?
							unsigned const obj_id(ball_it - interior->room_geom->objs.begin());
							interior->room_geom->particle_manager.add_for_obj(ball, 0.06*radius, front_dir, 1.0*KICK_VELOCITY, 50, 60, PART_EFFECT_SPARK, obj_id);
						}
						handled = 1;
					}
				}
			}
			if (obj.type == TYPE_PICTURE || obj.type == TYPE_TV || obj.type == TYPE_MONITOR || obj.type == TYPE_BUTTON || obj.type == TYPE_SWITCH ||
				obj.type == TYPE_BREAKER || (obj.type == TYPE_OFF_CHAIR && obj.rotates()))
			{
				if (!handled) {interact_with_object(obj_ix, new_center, new_center, velocity.get_norm());}
			}
		}
	}
	point const prev_new_center(new_center);

	if (move_sphere_to_valid_part(new_center, center, radius) && new_center != prev_new_center) { // collision with exterior wall
		apply_object_bounce_with_sound(*this, velocity, (new_center - prev_new_center).get_norm(), new_center, 1.0, on_floor); // hardness=1.0
		// add TYPE_CRACK if collides with a window?
	}
	if (new_center != center) { // is moving
		interior->update_dynamic_draw_data();

		if (has_water()) { // check for water collisions
			float const ceil_zval(get_bcube_z1_inc_ext_basement() + get_floor_for_zval(new_center.z)*get_window_vspace() + get_floor_ceil_gap());
			float const prev_zval(new_center.z);

			if (set_float_height(new_center, radius, ceil_zval, bt.density)) {
				float const target_zval(new_center.z);
				min_eq(new_center.z, (prev_zval + 1.0f*OBJ_GRAVITY*fticks_stable)); // limit max float velocity based on negative gravity
				velocity *= (1.0f - min(1.0f, 25.0f*OBJ_DECELERATE*fticks_stable)); // apply water dampening
				if (point_in_water_area(new_center)) {velocity.z = 0.0;} // remove vertical velocity component if center is underwater

				if (new_center.z == target_zval && velocity.mag_sq() < MIN_VELOCITY*MIN_VELOCITY) { // zero velocity if stopped and no longer rising
					velocity = zero_vector;
					obj_dynamic_to_static(ball, *interior);
				}
				static float last_splash_time(0.0);

				if ((tfticks - last_splash_time) > 0.5*TICKS_PER_SECOND) { // at most once every 0.5s
					bool const draw_splash(prev_zval-radius <= interior->water_zval && center.z-radius > interior->water_zval); // check if prev above the water line
					float const splash_size(0.75*sqrt(bt.weight));
					check_for_water_splash(new_center, splash_size, 1, draw_splash); // full_room_height=1
					last_splash_time = tfticks;
				}
				on_floor = 0; // not rolling on the floor
			}
			else if (point_in_water_area((new_center - radius*plus_z), 0)) { // check bottom point
				velocity *= (1.0f - min(1.0f, 10.0f*OBJ_DECELERATE*fticks_stable)); // apply (lower) water dampening
			}
		}
		apply_roll_to_matrix(dstate.rot_matrix, new_center, center, plus_z, radius, (on_floor ? 0.0 : 0.01), (on_floor ? 1.0 : 0.2));
		ball.translate(new_center - center);
		room_object_t squish_obj(ball);
		squish_obj.expand_by_xy(0.5*radius); // increase the radius to account for spiders and roaches being pushed out of the way of moving balls
		maybe_squish_animals(squish_obj, player_pos);
		int const new_room_id(get_room_containing_pt(new_center));
		if (new_room_id >= 0) {ball.room_id = new_room_id;} // needed for light_amt recompute when toggling lights; should we always recompute light_amt on room change?
		if (!was_dynamic) {interior->room_geom->invalidate_small_geom();} // static => dynamic transition, need to remove from static object vertex data
	}
	// check for collision with closed door separating the adjacent building at the end of the connecting room
	building_t *const cont_bldg(get_bldg_containing_pt(new_center));

	if (cont_bldg != nullptr && cont_bldg != this) { // switched buildings
		if (cont_bldg->interior->room_geom->add_room_object(ball, *this, 1, velocity)) {ball.remove();} // move ball from this to other_bldg
	}
}

void building_t::update_pool_table(room_object_t &ball) {
	// update velocity
	room_obj_dstate_t &dstate(interior->room_geom->get_dstate(ball));
	vector3d &velocity(dstate.velocity);
	ball_type_t const &bt(ball.get_ball_type());
	float const fticks_stable(min(fticks, 1.0f)); // cap to 1/40s to improve stability
	velocity *= (1.0f - min(1.0f, bt.friction*OBJ_DECELERATE*fticks_stable)); // apply dampening
	
	if (velocity.mag_sq() < 0.25*MIN_VELOCITY*MIN_VELOCITY) { // zero velocity if stopped
		velocity = zero_vector;
		obj_dynamic_to_static(ball, *interior);
		return;
	}
	// apply movement
	point const pos(ball.get_cube_center());
	point new_pos(pos + velocity*fticks_stable); // move based on velocity
	// check for collisions
	float const radius(ball.get_radius());
	room_object_t const &table(interior->room_geom->objs[ball.state_flags]);
	unsigned const balls_start(ball.state_flags + 1), balls_end(balls_start + 16); // 16 pool balls placed after the pool table
	auto &objs(interior->room_geom->objs);
	assert(balls_end <= objs.size());
	vector3d coll_dir;

	// collision detection with balls
	for (auto i = objs.begin()+balls_start; i != objs.begin()+balls_end; ++i) {
		if (i->obj_id == ball.obj_id)  continue; // skip ourself
		if (i->type != TYPE_POOL_BALL) continue; // ball was taken
		float const radius2(i->get_radius()), rsum(radius + radius2); // radius2 should be equal to radius
		point const pos2(i->get_cube_center());
		if (!dist_less_than(new_pos, pos2, rsum)) continue; // no collision
		// TODO
	}
	// collision detection with table edges
	cube_t play_area(get_pool_table_top_surface(table));
	play_area.expand_by_xy(-radius);

	if (!play_area.contains_pt_xy(new_pos)) {
		point const pre_clamp_pos(new_pos);
		play_area.clamp_pt_xy(new_pos);
		coll_dir = (new_pos - pre_clamp_pos).get_norm();
	}
	if (coll_dir != zero_vector) {
		// TODO: apply collision to velocity
		velocity = coll_dir*velocity.mag();
	}
	if (pos == new_pos) return; // no change in pos
	apply_roll_to_matrix(dstate.rot_matrix, new_pos, pos, plus_z, radius);
	ball.translate(new_pos - pos);
	interior->update_dynamic_draw_data();
}

void building_t::update_player_interact_objects(point const &player_pos) { // Note: player_pos is in building space
	assert(interior);
	interior->update_elevators(*this, player_pos);
	update_creepy_sounds(player_pos);
	if (!has_room_geom()) return; // nothing else to do
	float const floor_spacing(get_window_vspace());
	float const player_radius(get_scaled_player_radius()), player_z1(player_pos.z - get_bldg_player_height()), player_z2(player_pos.z);
	bool const player_in_this_building(this == player_building);
	static point last_player_pos(all_zeros);
	bool const player_is_moving(player_pos != last_player_pos);
	point camera_rot(player_pos);
	int player_room_ix(-1);
	float hum_amt(0.0), hum_freq(0.0);
	
	if (player_in_this_building) {
		last_player_pos = player_pos;
		maybe_inv_rotate_point(camera_rot); // rotate camera pos into building space; should we use camera_rot elsewhere below?
		player_room_ix = get_room_containing_pt(camera_rot);
		if (player_in_elevator >= 3) {hum_amt = 0.2; hum_freq = 100.0;} // moving elevator sound
	}
	// update dynamic objects; run for current and connected buildings
	auto &objs(interior->room_geom->objs);

	for (auto c = objs.begin(); c != objs.end(); ++c) {
		if (camera_surf_collide && c->type == TYPE_OFF_CHAIR && c->rotates()) { // player push office chair
			if (c->z1() > player_z2 || c->z2() < player_z1) continue; // no zval overlap
			float const chair_radius(0.25*(c->dx() + c->dy())); // treat as a cylinder
			vector3d const delta((c->xc() - player_pos.x), (c->yc() - player_pos.y), 0.0); // in XY plane
			float const min_dist(1.2*(player_radius + chair_radius)); // use a larger radius, since collision detection should prevent intersections
			if (delta.xy_mag_sq() > min_dist*min_dist) continue; // no intersection
			float const dist(delta.xy_mag());
			cube_t new_obj(*c);
			new_obj.translate(delta*(1.01*(min_dist - dist)/dist)); // move slightly more so that we don't intersect on the next frame
			room_t const &room(get_room(c->room_id));
			if (!get_walkable_room_bounds(room).contains_cube_xy(new_obj)) continue; // not contained in interior part of the room
			if (is_obj_placement_blocked(new_obj, room, 1))                continue; // inc_open_doors=1
			auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
			bool overlaps(0);

			for (auto i = objs.begin(); i != objs_end; ++i) {
				if (!i->no_coll() && i->intersects(new_obj) && i != c) {overlaps = 1; break;} // make sure to exclude ourself
			}
			if (overlaps) continue;
			c->copy_from(new_obj);
			c->flags |= RO_FLAG_MOVED; // not needed?
			// chairs block the AI, so must update backrooms path finding
			if (room.is_ext_basement()) {invalidate_nav_grid(room.get_floor_containing_zval(c->zc(), get_window_vspace()));}
			static float last_sound_tfticks(0);

			if ((tfticks - last_sound_tfticks) > 0.5*TICKS_PER_SECOND) { // play at most twice per second
				if (check_for_water_splash(cube_bot_center(new_obj), 0.9)) {} // check for splash
				else { // no splash, play squeak
					gen_sound_thread_safe_at_player(SOUND_SQUEAK, 0.2, 0.4); // lower pitch; should really use a rolling sound
					register_building_sound(player_pos, 0.2);
				}
				last_sound_tfticks = tfticks;
			}
			continue;
		}
		if (player_room_ix >= 0 && (int)c->room_id == player_room_ix) {
			if (c->type == TYPE_CEIL_FAN) {
				if (player_pos.z < c->z2() && player_pos.z > c->z1() - floor_spacing && ceiling_fan_is_on(*c, objs)) {
					float const dist(p2p_dist(camera_rot, c->get_cube_center())), sound_dist(2.0*floor_spacing);
					if (dist < sound_dist) {hum_amt = 0.15*(1.0 - dist/sound_dist); hum_freq = 90.0;}
				}
			}
			else if (c->type == TYPE_BRK_PANEL) {
				if (c->is_open()) {
					float const dist(p2p_dist(camera_rot, c->get_cube_center())), sound_dist(0.75*floor_spacing);
					if (dist < sound_dist) {hum_amt = 0.2*(1.0 - dist/sound_dist); hum_freq = 60.0;}
				}
			}
			//else if (c->type == TYPE_FURNACE) {} // or AC unit?
			//else if (c->type == TYPE_FRIDGE ) {}
		}
		if (!c->has_dstate()) continue; // Note: no test of player_coll flag
		run_ball_update(c, player_pos, player_z1, player_is_moving);
	} // for c
	if (player_in_this_building) { // interactions only run for player building
		if (player_in_closet) { // check for collisions with expanded objects in closets
			auto &expanded_objs(interior->room_geom->expanded_objs);

			for (auto c = expanded_objs.begin(); c != expanded_objs.end(); ++c) {
				if (c->type == TYPE_CLOTHES && sphere_cube_intersect(player_pos, player_radius, *c)) { // shirt in a closet with the player
					assert(c != objs.begin());
					room_object_t &hanger(*(c-1)); // hanger is the previous object
					assert(hanger.type == TYPE_HANGER);
					bool const rot_dir(dot_product(c->get_dir(), (player_pos - c->get_cube_center())) < 0);
					unsigned const add_flag(rot_dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI), rem_flag(rot_dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
					c->flags     |= ((c->type == TYPE_CLOTHES) ? RO_FLAG_ROTATING : 0) | add_flag;
					c->flags     &= ~rem_flag;
					hanger.flags |= RO_FLAG_ROTATING | add_flag; // play a sound as well?
					hanger.flags &= ~rem_flag;
				}
			} // for c
		}
		if (use_last_pickup_object || (tt_fire_button_down && !flashlight_on)) { // use object not active, and not using fire key without flashlight (space bar)
			maybe_use_last_pickup_room_object(player_pos);
			use_last_pickup_object = 0; // reset for next frame
		}
		maybe_update_tape(player_pos, 0); // end_of_tape=0
		if (interior->room_geom->fire_manager.get_closest_fire(player_pos, player_radius, player_z1, player_z2)) {player_take_damage(0.006);} // small amount of fire damage

		if (player_room_ix >= 0 /*&& !is_house*/) { // check for sounds; should this be for office buildings only?
			room_t const &room(get_room(player_room_ix));
			unsigned const camera_floor(room.get_floor_containing_zval(camera_rot.z, get_window_vspace()));
			unsigned const room_type(room.get_room_type(camera_floor));
			assert(room_type < NUM_RTYPES);
			if      (room_type == RTYPE_UTILITY ) {hum_amt = 0.1; hum_freq =  60.0;}
			else if (room_type == RTYPE_SERVER  ) {hum_amt = 0.2; hum_freq = 120.0;}
			//else if (room_type == RTYPE_SECURITY) {}
			//else if (room_type == RTYPE_SWIM    ) {}
			//else if (room_type == RTYPE_LAUNDRY   ) {}
		}
	} // end player_in_this_building
	if (hum_amt > 0.0) {play_hum_sound(player_pos, hum_amt, 0.01*hum_freq);}
	doors_next_frame(player_pos); // run for current and connected buildings
	interior->room_geom->particle_manager.next_frame(*this);
	interior->room_geom->fire_manager.next_frame(interior->room_geom->particle_manager);
}

// particle_manager_t

void particle_manager_t::add_for_obj(room_object_t &obj, float pradius, vector3d const &dir, float part_vel,
	unsigned min_parts, unsigned max_parts, unsigned effect, int parent_obj_id)
{
	assert(min_parts > 0 && min_parts <= max_parts);
	unsigned const num(min_parts + ((max_parts == min_parts) ? 0 : (rgen.rand()%(max_parts - min_parts + 1))));

	for (unsigned n = 0; n < num; ++n) {
		vector3d part_dir(rgen.signed_rand_vector_norm());
		if (dot_product(dir, dir) < 0.0) {part_dir.negate();}
		point pos;
		switch (obj.shape) {
		case SHAPE_CUBE  : pos = rgen.gen_rand_cube_point(obj); break;
		case SHAPE_SPHERE: pos = obj.get_cube_center() + sqrt(obj.get_radius()*obj.get_radius()*rgen.rand_float())*part_dir; break;
		case SHAPE_CYLIN : pos = obj.get_cube_center() + sqrt(obj.get_radius()*obj.get_radius()*rgen.rand_float())*part_dir; pos.z = rgen.rand_uniform(obj.z1(), obj.z2()); break;
		default: assert(0); // unsupported shape
		}
		vector3d const v(part_vel*rgen.rand_uniform(0.8, 1.25)*part_dir);
		particles.emplace_back(pos, v, WHITE, pradius*rgen.rand_uniform(0.8, 1.25), effect, parent_obj_id);
	}
}
cube_t particle_manager_t::get_bcube() const {
	cube_t bcube;
	for (particle_t const &p : particles) {bcube.assign_or_union_with_sphere(p.pos, p.radius);}
	return bcube;
}
void particle_manager_t::next_frame(building_t &building) {
	if (particles.empty()) return;
	float const fticks_stable(min(fticks, 4.0f)); // clamp to 0.1s
	auto const &objs(building.interior->room_geom->objs);
	float const lifetimes[NUM_PART_EFFECTS] = {0.0, 2.5, 3.0, 2.0, 0.25, 2.0}; // none, sparks, clouds, smoke, splash, bubble

	for (particle_t &p : particles) {
		point const p_last(p.pos);
		p.pos  += fticks_stable*p.vel;
		p.time += fticks_stable;
		float const lifetime(p.time/(lifetimes[p.effect]*TICKS_PER_SECOND));
		if (lifetime > 1.0) {p.effect = PART_EFFECT_NONE; continue;} // end of life
		
		if (p.effect == PART_EFFECT_SPARK) {
			p.color = get_glow_color(2.0*lifetime, 1); // fade=1
			apply_building_gravity(p.vel.z, 0.5*fticks_stable); // half gravity
		}
		else if (p.effect == PART_EFFECT_CLOUD) { // form fire extinguisher
			p.radius = p.init_radius*(1.0 + 4.0*lifetime); // radius increases over lifetime
			apply_building_gravity(p.vel.z, 0.05*fticks_stable); // very small gravity
		}
		else if (p.effect == PART_EFFECT_SMOKE) { // floats up
			p.radius = p.init_radius*(1.0 + 3.0*lifetime); // radius increases over lifetime
			p.color  = colorRGBA(WHITE*(0.25*(1.0 - lifetime)), p.alpha*(1.0 - lifetime)); // gray => transparent black
		}
		else if (p.effect == PART_EFFECT_SPLASH) {
			p.radius  = p.init_radius*(1.0 + 1.0*lifetime); // radius increases over lifetime
			p.color.A = 1.0 - lifetime; // transition to transparent
			continue; // no collision detection
		}
		else if (p.effect == PART_EFFECT_BUBBLE) { // floats up
			for (cube_t const &c : building.interior->ceilings) { // check if we hit a building ceiling
				if (sphere_cube_intersect(p.pos, p.radius, c)) {p.pos.z = c.z1() - p.radius;} // keep under the bottom of the ceiling
			}
			if (p.pos.z > building.interior->water_zval) {p.effect = PART_EFFECT_NONE;} // remove when it hits the water line
			continue; // done - no other collisions
		}
		else {assert(0);}
		// check for collisions and apply bounce, similar to balls
		float const bounce_scale = 0.5;
		vector3d cnorm;
		int obj_ix(-1);
		float hardness(0.0);
		auto self(objs.end());
		if (p.parent_obj_id >= 0) {assert((unsigned)p.parent_obj_id < objs.size()); self = objs.begin() + p.parent_obj_id;}

		if (building.interior->check_sphere_coll(building, p.pos, p_last, p.radius, self, cnorm, hardness, obj_ix)) {
			if (p.effect == PART_EFFECT_CLOUD || p.effect == PART_EFFECT_SMOKE) {p.effect = PART_EFFECT_NONE; continue;} // no bounce
			apply_floor_vel_thresh(p.vel, cnorm);
			bool const bounced(apply_object_bounce(p.vel, cnorm, bounce_scale*hardness, 0)); // on_floor=0
			
			// sparks create fires on floor coll on the first bounce when still hot, 10% of the time
			if (p.effect == PART_EFFECT_SPARK && bounced && p.bounce_count == 0 && lifetime < 0.5 && cnorm == plus_z && rgen.rand_float() < 0.1) {
				building.register_spark_floor_hit(p.pos);
			}
			if (bounced) {++p.bounce_count;}
		}
		point const prev_pos(p.pos);

		if (building.move_sphere_to_valid_part(p.pos, p_last, p.radius) && p.pos != prev_pos) { // collision with exterior wall
			apply_object_bounce(p.vel, (p.pos - prev_pos).get_norm(), bounce_scale, 0); // on_floor=0
		}
	} // for p
	// remove dead particles
	particles.erase(std::remove_if(particles.begin(), particles.end(), [](particle_t const &p) {return (p.effect == PART_EFFECT_NONE);}), particles.end());
}

// fire_manager_t

cube_t fire_manager_t::fire_t::get_bcube() const {
	cube_t bcube(pos);
	bcube.expand_by_xy(radius);
	bcube.z2() += get_height();
	return bcube;
}

void fire_manager_t::spawn_fire(point const &pos, float size) {
	size *= rgen.rand_uniform(0.8, 1.25); // randomize size a bit

	for (fire_t &f : fires) { // merge if near other fires
		if (!dist_less_than(pos, f.pos, (size + f.radius))) continue;
		f.max_radius = sqrt(f.max_radius*f.max_radius + size*size); // increase radius
		return;
	}
	fires.emplace_back(pos, size);
}
cube_t fire_manager_t::get_bcube() const {
	cube_t bcube;
	for (fire_t const &f : fires) {bcube.assign_or_union_with_cube(f.get_bcube());}
	return bcube;
}
bool fire_manager_t::get_closest_fire(point const &pos, float xy_radius, float z1, float z2, point *fire_pos) const {
	bool found(0);
	float dmin_sq(0.0);

	for (fire_t const &f : fires) {
		if (!dist_xy_less_than(pos, f.pos, (xy_radius + f.radius))) continue;
		if (f.pos.z > z2 || (f.pos.z + f.get_height()) < z1)        continue;
		float const dist_sq(p2p_dist_xy_sq(pos, f.pos));
		if (found && dist_sq >= dmin_sq) continue; // not the closest
		dmin_sq = dist_sq;
		found   = 1;
		if (fire_pos != nullptr) {*fire_pos = f.pos;}
	}
	return found;
}
void fire_manager_t::add_fire_bcubes_for_cube(cube_t const &sel_cube, vect_cube_t &fire_bcubes) const {
	for (fire_t const &f : fires) {
		cube_t const bcube(f.get_bcube());
		if (bcube.intersects(sel_cube)) {fire_bcubes.push_back(bcube);}
	}
}
void fire_manager_t::put_out_fires(point const &p1, point const &p2, float radius) { // area can be cylindrical
	for (fire_t &f : fires) {
		bool hit(dist_less_than(p1, f.pos, radius)); // sphere
		if (p1 != p2) {hit |= (dist_less_than(p2, f.pos, radius) || pt_line_seg_dist_less_than(f.pos, p1, p2, radius));} // extend to capsule shape
		if (hit) {f.max_radius = 0.0;}
	}
}
void fire_manager_t::next_frame(particle_manager_t &particle_manager) {
	float const fticks_stable(min(fticks, 4.0f)); // clamp to 0.1s

	for (fire_t &f : fires) {
		f.time += fticks_stable;
		float const lifetime(f.time/(4.0*TICKS_PER_SECOND)); // 4s
		if (lifetime < 0.3) {f.radius = (lifetime/0.3)*f.max_radius;} // grow at start of life
		if (lifetime > 0.7) {f.radius = (1.0 - (lifetime - 0.7)/0.3)*f.max_radius;} // shink at end of life
		// slow random drift? would need to check object collisions
		if (f.time < f.next_smoke_time) continue;
		float const radius(1.0*f.radius);
		particle_manager.add_particle((f.pos + 1.1*radius*plus_z), SMOKE_VELOCITY*plus_z, GRAY, radius, PART_EFFECT_SMOKE); // generate smoke
		f.next_smoke_time = f.time + rgen.rand_uniform(0.25, 0.5)*TICKS_PER_SECOND;
	} // for f
	// remove dead fires
	fires.erase(std::remove_if(fires.begin(), fires.end(), [](fire_t const &f) {return (f.max_radius == 0.0 || f.radius < 0.0);}), fires.end());
}

void building_t::register_spark_floor_hit(point const &pos) {
	if (!has_room_geom()) return;
	float const wall_thickness(get_wall_thickness()), z_range(0.1*wall_thickness), z1(pos.z - z_range), z2(pos.z + z_range), fire_size(1.2*wall_thickness);
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_RUG) continue; // currently, only rugs can be burned
		if (!i->contains_pt_xy(pos) || i->z1() > z2 || i->z2() < z1) continue; // no collision
		point const fpos(pos.x, pos.y, (i->z2() + 0.01*wall_thickness));
		interior->room_geom->fire_manager .spawn_fire   (fpos,     fire_size);
		interior->room_geom->decal_manager.add_burn_spot(fpos, 2.4*fire_size); // add black spot on the rug
		record_building_damage(10.0);
		break; // spawn only one fire
	} // for i
}

void building_decal_manager_t::add_burn_spot(point const &pos, float radius) {
	// if there are too many existing spots, remove the first (oldest) one; this can be slow, but shouldn't happen very often
	unsigned const max_spots = 100;
	unsigned const num_spots(burn_qbd.verts.size()/6); // 6 verts/2 triangles per quad
	if (num_spots >= max_spots) {burn_qbd.verts.erase(burn_qbd.verts.begin(), (burn_qbd.verts.begin() + 6*(num_spots - max_spots + 1)));}
	burn_qbd.add_quad_dirs(pos, -plus_x*radius, plus_y*radius, BLACK); // -x!
}

// elevators

bool building_interior_t::update_elevators(building_t const &building, point const &player_pos) { // Note: player_pos is in building space
	if (elevators_disabled) return 0;
	float const z_space(0.05*building.get_floor_thickness()); // to prevent z-fighting
	float const delta_open_amt(min(1.0f, 2.0f*fticks/TICKS_PER_SECOND)); // 0.5s for full open
	float const elevator_wait_time(5.0*TICKS_PER_SECOND); // 5s
	static int prev_move_dir(2); // starts at not-moving
	vect_room_object_t &objs(room_geom->objs);
	bool was_updated(0), update_ddd(0);

	// Note: the player can only be in one elevator at a time, but they can push the call button for one elevator and get into another, so we have to check all elevators
	for (auto e = elevators.begin(); e != elevators.end(); ++e) { // find containing elevator (optimization + need to know z-range of elevator shaft)
		e->is_moving = 0; // reset this frame

		if (e->at_dest || !e->was_called()) { // stopped on a floor (either in between stops or at the end of the call requests)
			bool const time_for_doors_to_close(e->at_dest_frame > 0 && frame_counter > (e->at_dest_frame + elevator_wait_time) && !e->hold_doors);

			if (!e->was_called() && e->open_amt > 0.0 && time_for_doors_to_close) { // inactive, close the doors
				e->open_amt = max((e->open_amt - delta_open_amt), 0.0f);
				if (e->open_amt == 0.0) {e->at_dest_frame = 0;} // done waiting/closing
				e->at_dest  = 0; // reset for next cycle
				update_ddd  = 1; // regen verts for door
			}
			else if (e->open_amt == 1.0) { // doors are fully open
				if (e->at_dest_frame == 0) { // not yet waiting
					e->at_dest_frame = frame_counter; // record the time when the doors are first fully open
				}
				else if (e->was_called() && time_for_doors_to_close) { // called to another floor, and we've waited long enough
					e->at_dest_frame = 0; // done waiting
					e->at_dest = 0; // should no longer get into this case on the next frame and should instead get into the movement logic below
				}
			}
			else if (e->open_amt > 0.0) { // doors are partially open - continue to open fully
				e->open_amt = min((e->open_amt + delta_open_amt), 1.0f);
				update_ddd  = 1; // regen verts for door
			}
			continue;
		}
		if (e->hold_doors) {
			if (e->open_amt < 1.0) { // open doors
				e->open_amt = min((e->open_amt + delta_open_amt), 1.0f);
				update_ddd  = 1; // regen verts for door
			}
			continue;
		}
		e->at_dest_frame = 0; // reset just in case, since we're not in a waiting state
		assert(e->car_obj_id < objs.size());
		room_object_t &obj(objs[e->car_obj_id]); // elevator car for this elevator
		assert(obj.type == TYPE_ELEVATOR && obj.room_id == (e - elevators.begin())); // sanity check
		float const target_zval(e->get_target_zval()), ez1(obj.z1());

		if (e->open_amt > 0.0 && target_zval != ez1) { // doors not yet closed, and not at target zval
			e->open_amt = max((e->open_amt - delta_open_amt), 0.0f); // close the doors
			update_ddd  = 1; // regen verts for door
			continue;
		}
		if (e->hold_movement) continue; // hold for this frame only
		bool const move_dir(target_zval > ez1); // 0=down, 1=up
		if (target_zval != ez1) {e->going_up = move_dir;} // only update if the elevator is moving/not stopped at the dest
		assert(e->button_id_start < e->button_id_end && e->button_id_end <= objs.size());
		float dist(min(0.5f*CAMERA_RADIUS, 0.04f*obj.dz()*fticks)*(move_dir ? 1.0 : -1.0)); // clamp to half camera radius to avoid falling through the floor for low framerates
		// check if any called floors are closer than the current one and in the correct direction; if so, move them to the front of the queue
		e->move_closest_in_dir_to_front(ez1, move_dir);
		if (e->was_called() && fabs(target_zval - ez1) < fabs(dist)) {dist = (target_zval - ez1);} // move to position
		else if (move_dir) {min_eq(dist, (e->z2() - obj.z2() - z_space));} // going up
		else               {max_eq(dist, (e->z1() - obj.z1() + z_space));} // going down
		update_ddd = 1;

		if (fabs(dist) < 0.001*z_space) { // no movement; maybe at target_zval or top/bottom of elevator shaft (check with a tolerance)
			max_eq(e->open_amt, delta_open_amt); // begin to open if not already open
			e->register_at_dest(); // will remove the current call request
			obj.flags |= RO_FLAG_OPEN;

			// disable all call buttons for this elevator that are no longer called
			for (auto j = objs.begin() + e->button_id_start; j != objs.begin() + e->button_id_end; ++j) {
				if (j->type == TYPE_BLOCKER) continue; // button was removed?
				assert(j->type == TYPE_BUTTON);
				unsigned const up_down_mask((j->flags & RO_FLAG_ADJ_TOP) ? 2 : ((j->flags & RO_FLAG_ADJ_BOT) ? 1 : 3)); // top=up, bot=down, neither=both
				if (!j->is_active() || e->was_floor_called(j->obj_id, up_down_mask)) continue; // already unlit, or this floor has also been called
				j->flags &= ~RO_FLAG_IS_ACTIVE; // clear active/lit state
				room_geom->invalidate_small_geom(); // need to regen object data due to lit state change
			}
			point const sound_pos(obj.get_cube_center());
			gen_sound_thread_safe(SOUND_BEEP, building.local_to_camera_space(sound_pos), 0.5, 0.75); // lower frequency beep
			register_building_sound(sound_pos, 0.5);
			continue;
		}
		obj.translate_dim(2, dist); // translate in Z
		obj.item_flags = uint16_t(floor((obj.zc() - e->z1())/building.get_window_vspace())); // set current floor
		e->is_moving   = 1;
		assert(e->light_obj_id < objs.size());
		room_object_t &light(objs[e->light_obj_id]); // light for this elevator

		if (light.type != TYPE_BLOCKER) { // translate light as well, if not taken
			assert(light.type == TYPE_LIGHT);
			light.translate_dim(2, dist); // translate in Z
			room_geom->invalidate_lights_geom(); // or make elevator lights part of the elevator instead?
			building.register_indir_lighting_state_change(e->light_obj_id);
		}
		for (auto j = objs.begin() + e->button_id_start; j != objs.begin() + e->button_id_end; ++j) {
			if (j->type == TYPE_BLOCKER) continue; // button was removed?
			assert(j->type == TYPE_BUTTON);
			if (j->in_elevator()) {j->translate_dim(2, dist);} // interior panel button, translate in Z
		}
		if ((int)move_dir != prev_move_dir && obj.contains_pt(player_pos)) { // moving, and player is in the elevator
			gen_sound_thread_safe_at_player(SOUND_SLIDING, 0.75); // play this sound when the elevator starts moving or changes direction
			register_building_sound(player_pos, 0.4);
		}
		prev_move_dir = move_dir;
		was_updated   = 1;
	} // for e
	for (auto &e : elevators) {e.hold_doors = e.hold_movement = 0;} // reset frame state
	if (update_ddd) {update_dynamic_draw_data();} // clear dynamic material vertex data (for all elevators) and recreate their VBOs
	if (was_updated) return 1;
	prev_move_dir = 2;
	return 0;
}

bool elevator_t::was_floor_called(unsigned floor_ix, unsigned up_down_mask) const {
	for (call_request_t const &cr : call_requests) {
		if (cr.floor_ix == floor_ix && (cr.req_dirs & up_down_mask) != 0) return 1;
	}
	return 0;
}
void elevator_t::call_elevator(unsigned floor_ix, float targ_z, unsigned req_dirs, bool inside_press) {
	// Note: the only case I'm not sure about is if someone on floor 2 pushes the down call button and someone on floor 3 pushes the down call button
	// before the other person gets into the elevator; I think the elevator should go to floor 2 and then down;
	// but what if the person getting in on floor 2 actually presses the button for floor 3? does it go up in that case even though they pressed the down call button?
	for (call_request_t &cr : call_requests) {
		if (cr.floor_ix != floor_ix) continue; // duplicate press for this floor
		bool const prev_ip(cr.inside_press);
		cr.req_dirs     |= req_dirs; // combine up/down flags
		cr.inside_press |= inside_press;
		if (inside_press && !prev_ip) {stable_sort(call_requests.begin(), call_requests.end());} // prioritize inside press if changed
		return;
	} // for cr
	call_requests.emplace_back(floor_ix, targ_z, req_dirs, inside_press); // place the request at the end to make this the last floor visited
	if (inside_press) {stable_sort(call_requests.begin(), call_requests.end());} // prioritize inside press
}
void elevator_t::register_at_dest() {
	if (at_dest) return; // duplicate call; can this happen?
	//assert(!call_requests.empty()); // too strong?

	if (!call_requests.empty()) {
		uint8_t &req_dirs(call_requests.front().req_dirs);
		//cout << TXTi(req_dirs) << TXT(stop_on_passing_floor) << TXT(going_up) << TXT(call_requests.size()) << TXT(call_requests[0].floor_ix) << endl; // TESTING

		if (stop_on_passing_floor && call_requests.size() > 1) { // only valid if there's another call request
			// if both up and down buttons were pressed (req_dirs == 3), we should unset one dir but leave the other;
			// if this is a stop along the way (stop_on_passing_floor==1), we should unset the one in the current elevator dir;
			// TODO: but if this is the end stop of the elevator, we don't know which dir because we're not tracking which of the up/down button was pressed first;
			// furthermore, even setting req_dirs in this case can cause the elevator to deadlock if the occupant presses a button that causes it to move in the wrong direction
			// Note: we can't rely on going_up because we may be at a stop where the elevator reverses direction, and going_up hasn't yet been updated
			bool const real_going_up(call_requests[0].floor_ix < call_requests[1].floor_ix); // call_requests[0].floor_ix should be our current floor
			req_dirs &= (real_going_up ? 1 : 2);
			// if we still have a call, move this CR to the back so that it's picked up after we reverse direction; will be popped below
			if (req_dirs != 0) {call_requests.push_back(call_requests.front());}
		}
		call_requests.pop_front();
	}
	at_dest = 1;
	stop_on_passing_floor = 0; // reset for next cycle
}
void elevator_t::move_closest_in_dir_to_front(float zval, bool dir) {
	assert(was_called());
	auto new_front(call_requests.end());
	float dzmin(dz()); // start at a large value

	for (auto cr = call_requests.begin(); cr != call_requests.end(); ++cr) {
		if ((cr->zval > zval) != dir) continue; // floor in the wrong direction compared to elevator movement
		if (cr != call_requests.begin() && !(cr->req_dirs & (1<<unsigned(dir)))) continue; // call in the wrong direction compared to elevator movement (okay if the front)
		float const dz(fabs(cr->zval - zval));
		if (dz < dzmin) {new_front = cr; dzmin = dz;}
	}
	assert(new_front != call_requests.end()); // a dest floor must have been found

	if (new_front != call_requests.begin()) { // closest is not the front; swap with the front
		call_request_t const v(*new_front); // deep copy
		call_requests.erase(new_front);
		call_requests.push_front(v);
		stop_on_passing_floor = 1; // flag so that we can set the up/down call button state correctly
	}
}

void building_t::register_button_event(room_object_t const &button) {
	// here room_id is elevator_id (buttons are only used with elevators)
	bool const is_up(button.flags & RO_FLAG_ADJ_TOP);
	call_elevator_to_floor(get_elevator(button.room_id), button.obj_id, button.in_elevator(), is_up); // floor_ix=button.obj_id
}
void building_t::call_elevator_to_floor(elevator_t &elevator, unsigned floor_ix, bool is_inside_elevator, bool is_up) {
	if (interior->elevators_disabled) return; // nope
	float const targ_z(elevator.z1() + max(get_window_vspace()*floor_ix, 0.05f*get_floor_thickness())); // bottom of elevator car for this floor
	assert(targ_z <= bcube.z2()); // sanity check
	// inside: check if button is above or below current elevator car pos, otherwise depends on the call button pressed
	bool const dest_is_up(targ_z > interior->room_geom->objs[elevator.car_obj_id].z1());
	unsigned const req_dirs(is_inside_elevator ? (dest_is_up ? 2 : 1) : (is_up ? 2 : 1));
	elevator.call_elevator(floor_ix, targ_z, req_dirs, is_inside_elevator);
}

void clamp_sphere_xy(point &pos, cube_t const &c, float radius) {
	cube_t bounds(c);
	bounds.expand_by_xy(-radius); // must fit entire sphere
	bounds.clamp_pt_xy(pos);
}
// applies to balls; Note: only moves in XY
bool building_t::move_sphere_to_valid_part(point &pos, point const &p_last, float radius) const {
	point const init_pos(pos);

	if (has_attic() && pos.z > interior->attic_access.z2()) { // special case handling for attic
		for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
			if (!i->contains_pt_xy(p_last)) continue; // not the part containing the previous pos (assumes parts are not stacked)
			clamp_sphere_xy(pos, *i, (4.0*radius + get_attic_beam_depth())); // include extra space for attics (approximate)
			if (pos != init_pos) return 1;

			if (pos.z > p_last.z) { // rising; check if above attic roof
				float const move_zmin(max(p_last.z, interior->attic_access.z2()+radius));

				while (pos.z > move_zmin && !point_under_attic_roof(point(pos.x, pos.y, pos.z+radius))) {
					pos.z = min(pos.z-0.1f*radius, move_zmin); // shift down by 10% of radius until we hit the prev zval or attic floor, or we no longer collide with the roof
				}
			}
			return (pos != init_pos);
		} // for i
	}
	else {
		float xy_area_contained(0.0);
		cube_t sphere_bcube; sphere_bcube.set_from_sphere(pos, radius);
		cube_t const clamp_cube((has_basement() && pos.z < ground_floor_z1) ? get_full_basement_bcube() : bcube);
		clamp_sphere_xy(pos, clamp_cube, radius); // keep pos within the valid building bcube

		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (pos.z >= i->z1() && pos.z < i->z2()) {accumulate_shared_xy_area(*i, sphere_bcube, xy_area_contained);}
		}
		if (has_ext_basement() && (pos.z - radius) < ground_floor_z1) { // use the ext basement hallway if pos is in the basement, otherwise use the entire ext basement
			cube_t const &basement_cube(point_in_extended_basement_not_basement(pos) ? interior->basement_ext_bcube : get_ext_basement_entrance());
			accumulate_shared_xy_area(basement_cube, sphere_bcube, xy_area_contained);
		}
		if (xy_area_contained > 0.99*sphere_bcube.dx()*sphere_bcube.dy()) { // sphere contained in union of parts (not outside the building)
			if (pos != init_pos) return 1;

			if (!check_cube_within_part_sides(sphere_bcube)) { // outside the building
				int const part_ix(get_part_ix_containing_pt(pos));
				if (part_ix >= 0 && do_sphere_coll_polygon_sides(pos, parts[part_ix], radius, 1, get_part_ext_verts(part_ix), nullptr)) return 1; // interior_coll=1
			}
			return 0;
		}
		// find part containing p_last and clamp to that part
		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (!i->contains_pt(p_last)) continue; // not the part containing the previous pos
			clamp_sphere_xy(pos, *i, radius);
			return 1;
		}
		if (point_in_extended_basement(p_last)) { // Note: no need to check for individual ext basement rooms, as interior wall collisions should handle this
			clamp_sphere_xy(pos, interior->basement_ext_bcube, radius);
			return 1;
		}
	}
	//assert(0); // should never get here (except for FP error?)
	pos = p_last; // restore original position
	return 1;
}

// ray queries

// center, -x, +x, -y, +y, -z, +z; assumes pts is at least size 7/5 for skip_z=0/1
void get_sphere_boundary_pts(point const &center, float radius, point *pts, bool skip_z=0) {
	pts[0] = center;

	for (unsigned dim = 0, ix = 1; dim < (skip_z ? 2U : 3U); ++dim) {
		vector3d dir(zero_vector);
		dir[dim]  = radius;
		pts[ix++] = center - dir;
		pts[ix++] = center + dir;
	}
}

// returns true if ray ends inside cubes; similar to code in building_t::refine_light_bcube(), but also clips in Z
bool trace_ray_through_cubes(vect_cube_t const &cubes, point const &p1, point const &p2, float tolerance) {
	point cur_pt(p1);
	auto prev_part(cubes.end());

	for (unsigned n = 0; n < cubes.size(); ++n) { // could be while(1), but this is safer in case we run into an infinite loop due to FP errors
		bool found(0);

		for (auto p = cubes.begin(); p != cubes.end(); ++p) {
			cube_t tc(*p);
			tc.expand_by(tolerance);
			if (p == prev_part || !tc.contains_pt(cur_pt)) continue; // ray does not continue into this new part
			point tmp_pt(cur_pt), new_pt(p2);
			if (do_line_clip(tmp_pt, new_pt, p->d)) {cur_pt = new_pt; prev_part = p; found = 1; break;} // ray continues into this part
		}
		if (!found)       return 0; // ray has exited the cubes, done
		if (cur_pt == p2) return 1; // we've reached the end point, done
	} // for n
	return 0; // ray has exited the cubes
}
bool building_t::check_line_intersect_doors(point const &p1, point const &p2, bool inc_open) const {
	// is it more efficient to iterate over interior->door_stacks? only if inc_open=0 (the common case)?
	float const zmin(min(p1.z, p2.z)), zmax(max(p1.z, p2.z));

	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
		if (!inc_open && i->open)             continue; // check only closed doors
		if (i->z1() > zmax || i->z2() < zmin) continue; // check Z-range

		if (i->open) {
			cube_t door_bounds(*i);
			door_bounds.expand_by_xy(i->get_width());
			if (!door_bounds.line_intersects(p1, p2)) continue; // check intersection with rough/conservative door bounds (optimization)
			tquad_with_ix_t const door(set_interior_door_from_cube(*i));
			//float const thickness(i->get_thickness()); // for now, we ignore the door thickness since it shouldn't really make a difference
			float t(0.0); // unused
			if (line_poly_intersect(p1, p2, door.pts, door.npts, door.get_norm(), t)) return 1;
		}
		else if (i->get_true_bcube().line_intersects(p1, p2)) return 1; // test line bcube?
	} // for i
	return 0;
}
bool building_t::is_pt_visible(point const &p1, point const &p2) const {
	if (!interior) return 1;
	if (is_light_occluded(p1, p2)) return 0; // okay to call for non-light point; checks walls, ceilings, and floors
	
	if (parts.size() == 1) {} // single part, no check necessary
	else if (p1.z < ground_floor_z1) { // in the basement
		if (has_ext_basement()) {
			vect_cube_t cubes; // static?
			cubes.push_back(get_basement());
			cubes.push_back(interior->basement_ext_bcube);
			if (!trace_ray_through_cubes(cubes, p1, p2, 0.01*get_wall_thickness())) return 0;
		} // else only basement, not need to trace exterior walls
	}
	else {
		if (!trace_ray_through_cubes(parts, p1, p2, 0.01*get_wall_thickness())) return 0; // view blocked by exterior wall (ignores windows)
	}
	if (check_line_intersect_doors(p1, p2)) return 0;
	return 1;
}
bool building_t::is_sphere_visible(point const &center, float radius, point const &pt) const {
	if (!interior) return 1;
	point pts[7];
	get_sphere_boundary_pts(center, radius, pts);
	for (unsigned n = 0; n < 7; ++n) {if (is_pt_visible(pts[n], pt)) return 1;}
	return 0;
}

bool building_t::is_pt_lit(point const &pt) const {
	if (!has_room_geom()) return 0; // no lights
	int const room_id(get_room_containing_pt(pt)); // call this only once on center in is_sphere_lit()?
	if (room_id < 0) return 0; // outside building?
	room_t const &room(get_room(room_id));
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_type() || !i->is_light_on() || i->in_closet()) continue; // not a light, or light not on; skip closet lights
		bool const same_room((int)i->room_id == room_id);
		//if (!same_room) continue; // different room (optimization); too strong?
		//bool const same_floor(fabs(i->z1() - pt.z) < floor_spacing); // doesn't work with lamps
		bool const same_floor(room.is_single_floor || get_floor_for_zval(pt.z) == get_floor_for_zval(i->z1()));
		if (!i->has_stairs() && !same_floor) continue; // different floors, and no stairs (optimization)
		if (same_floor && same_room) return 1; // same floor of same room, should be visible (optimization)
		point const center(i->get_cube_center());
		if (!dist_less_than(center, pt, 0.95*get_radius_for_room_light(*i))) continue;
		if (is_pt_visible(center, pt)) return 1; // likely returns true if same room
	} // for i
	return 0;
}
bool building_t::is_sphere_lit(point const &center, float radius) const {
	if (!has_room_geom()) return 0; // no lights (optimization)
	point pts[7];
	get_sphere_boundary_pts(center, radius, pts);
	for (unsigned n = 0; n < 7; ++n) {if (is_pt_lit(pts[n])) return 1;}
	return 0;
}
bool building_t::is_room_lit(int room_id, float zval) const {
	if (room_id < 0) return 0; // outside building?
	room_t const &room(get_room(room_id));
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_type() || !i->is_light_on() || i->light_is_out() || i->in_closet()) continue; // not a light, light not on, or broken; skip closet lights
		if ((int)i->room_id != room_id) continue; // different room
		if (room.is_single_floor || get_floor_for_zval(zval) == get_floor_for_zval(i->z1())) return 1; // same floor - lit
	}
	return 0;
}

// room objects

// used for pushing objects and opening books
bool building_t::is_obj_pos_valid(room_object_t const &obj, bool keep_in_room, bool allow_block_door, bool check_stairs) const {
	assert(interior);
	room_t const &room(get_room(obj.room_id));

	if (keep_in_room) { // check room bounds; should this apply to the attic as well?
		cube_t place_area(room);
		place_area.expand_by_xy(-0.99*get_trim_thickness()); // shrink to exclude wall trim
		if (!place_area.contains_cube(obj)) return 0; // outside the room
	}
	bool contained_in_part(0);

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		if (p->contains_cube(obj)) {contained_in_part = 1; break;}
	}
	if (!contained_in_part && !cube_in_attic(obj) && !interior->basement_ext_bcube.contains_cube(obj)) return 0;

	for (unsigned d = 0; d < 2; ++d) { // check for wall intersection
		if (has_bcube_int_no_adj(obj, interior->walls[d])) return 0;
	}
	for (auto const &ds : interior->door_stacks) {
		if (ds.z1() > obj.z2() || ds.z2() < obj.z1()) continue; // wrong floor for this stack/part
		// calculate door bounds for bcube test, assuming it's open
		cube_t door_bounds(ds);
		door_bounds.expand_by_xy(ds.get_width());
		if (!door_bounds.intersects(obj)) continue; // optimization
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t const &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done

			if (!door.open) { // closed
				if (!keep_in_room && door.get_true_bcube().intersects(obj)) return 0; // door blocking room boundary
				continue; // ignore
			}
			if (!door_opens_inward(door, room)) continue; // or opens into the adjacent room, ignore
			
			if (allow_block_door) {
				// check current position of open door; ignore closed position; if the door is closed later, it may intersect the object
				if (get_door_bounding_cube(door).intersects(obj)) return 0;
			}
			else if (is_cube_close_to_door(obj, 0.0, 1, door, door.get_check_dirs(), door.open_dir, allow_block_door)) return 0;
		} // for dix
	} // for door_stacks
	if (check_stairs) {
		cube_t obj_exp(obj);
		obj_exp.expand_by_xy(0.25*get_wall_thickness()); // add a bit of extra padding to account for the railing thickness extending outside the stairs bcube (approx)
		if (has_bcube_int(obj_exp, interior->stairwells)) return 0;
	}
	if (has_bcube_int(obj, interior->elevators)) return 0;
	return 1;
}

unsigned building_room_geom_t::allocate_dynamic_state() {
	unsigned const ix(obj_dstate.size());
	obj_dstate.push_back(room_obj_dstate_t());
	return ix;
}
room_obj_dstate_t &building_room_geom_t::get_dstate(room_object_t const &obj) {
	assert(obj.has_dstate());
	assert(obj.obj_id < obj_dstate.size());
	return obj_dstate[obj.obj_id];
}

bool is_blocking_obj_on_top_surface(room_object_t const &obj) { // objects on tables, counters, desks, etc.
	return (obj.type == TYPE_PLANT || obj.type == TYPE_KEYBOARD || obj.type == TYPE_BOTTLE || obj.type == TYPE_MWAVE || is_ball_type(obj.type) ||
		obj.type == TYPE_PLATE || obj.type == TYPE_LAPTOP || obj.type == TYPE_PAN || obj.type == TYPE_VASE || obj.type == TYPE_URN || obj.type == TYPE_MONITOR ||
		obj.type == TYPE_LAMP || obj.type == TYPE_CUP || obj.type == TYPE_TOASTER || obj.type == TYPE_SILVER || obj.type == TYPE_PIZZA_BOX || obj.type == TYPE_PLANT_MODEL);
}
bool can_place_on_object(room_object_t const &obj, point const &pos, float radius, float z_bias, float start_zval, float &zval) {
	if (!obj.can_place_onto() && obj.type != TYPE_RUG) return 0; // can't place on this object type
	if (!obj.contains_pt_xy(pos)) return 0; // center of mass not contained
	cube_t c(obj);

	if (obj.type == TYPE_BED) {
		cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
		get_bed_cubes(obj, cubes);
		if (!cubes[3].contains_pt_xy(pos)) return 0; // check again
		if (cubes[1].contains_pt_xy_exp(pos, radius) || cubes[2].contains_pt_xy_exp(pos, radius)) return 0; // intersects the head or foot, skip
		c = cubes[3]; // mattress
	}
	if (c.z2() < zval || c.z2() > start_zval) return 0; // below the floor or above the object's starting position; maybe can place under this object

	if (obj.shape == SHAPE_CYLIN) {
		if (!dist_xy_less_than(pos, obj.get_cube_center(), obj.get_radius())) return 0; // round table
	}
	else {assert(obj.shape != SHAPE_SPHERE);} // SHAPE_CUBE, SHAPE_TALL, SHAPE_SHORT are okay; others don't make sense
	zval = c.z2() + z_bias; // place on top of this object
	return 1;
}

bool building_t::get_zval_for_pool_bottom(point const &pos, float &zval) const {
	if (!has_pool     ()) return 0;
	if (!has_room_geom()) return 0; // if the room geom hasn't been generated yet, we can't determine pool depth, so just ignore the pool
	indoor_pool_t const &pool(interior->pool);
	if (!pool.contains_pt_xy(pos)) return 0;
	room_t const &pool_room(get_room(pool.room_ix));
	if (pos.z > pool_room.z2() || pos.z < pool.z1()) return 0;
	vect_room_object_t const &objs(interior->room_geom->objs);
	unsigned const pr_ix(interior->room_geom->pool_ramp_obj_ix), ps_six(interior->room_geom->pool_stairs_start_ix);
	zval = pool.z1(); // start on the bottom of the pool

	if (pr_ix > 0) { // pool has a ramp on the bottom
		assert(pr_ix+1 < objs.size()); // must have ramp and upper surface
		room_object_t const &ramp(objs[pr_ix]), &upper(objs[pr_ix+1]);
		assert(ramp .type == TYPE_POOL_TILE && ramp .shape == SHAPE_ANGLED);
		assert(upper.type == TYPE_POOL_TILE && upper.shape == SHAPE_CUBE);

		if (upper.contains_pt_xy(pos)) {zval = upper.z2();} // on upper surface
		else if (ramp.contains_pt_xy(pos)) { // on the ramp
			float const t(CLIP_TO_01((pos[ramp.dim] - ramp.d[ramp.dim][0])/ramp.get_length()));
			zval = ramp.z1() + ramp.dz()*(ramp.dir ? t : (1.0-t));
		}
	}
	if (ps_six > 0) { // check for stairs; should always be true
		assert(ps_six < objs.size());
		assert(objs[ps_six].type == TYPE_STAIR);

		for (auto i = objs.begin()+ps_six; i != objs.end(); ++i) {
			if (i->type != TYPE_STAIR) break; // done with stairs
			if (i->contains_pt_xy(pos)) {max_eq(zval, i->z2()); break;} // can only be on one stair
		}
	}
	return 1;
}
bool building_t::get_zval_of_floor(point const &pos, float radius, float &zval) const {
	if (!interior) return 0; // error?
	cube_t cur_bcube;
	cur_bcube.set_from_sphere(pos, radius);
	if (!bcube.contains_cube_xy(cur_bcube) && !interior->basement_ext_bcube.contains_cube(cur_bcube)) return 0; // not contained/too close to walls
	if (get_zval_for_pool_bottom(pos, zval)) return 1; // on the bottom of the pool
	float const floor_spacing(get_window_vspace());

	for (auto i = interior->floors.begin(); i != interior->floors.end(); ++i) { // blood can only be placed on floors
		if (pos.z < i->z2() || pos.z > (i->z2() + floor_spacing) || !i->contains_cube_xy(cur_bcube)) continue; // wrong floor, or not contained
		zval = (i->z2() + 0.0015*floor_spacing); // slightly above rugs (0.0015 vs. 0.001) and flooring (0.0015 vs. 0.0012)
		return 1;
	}
	return 0; // no suitable floor found
}
bool building_t::get_zval_for_obj_placement(point const &pos, float radius, float &zval, bool add_z_bias) const {
	if (!has_room_geom()) return 0; // error?
	float const start_zval(pos.z);
	if (!get_zval_of_floor(pos, radius, zval)) return 0; // if there's no floor, then there's probably no object to place on either

	for (unsigned d = 0; d < 2; ++d) { // check walls
		for (cube_t const &wall : interior->walls[d]) {
			if (pos.z > wall.z1() && pos.z < wall.z2() && sphere_cube_intersect(pos, radius, wall)) return 0; // intersects a wall, can't place here
		}
	}
	float const z_bias(add_z_bias ? 0.0005*get_window_vspace() : 0.0); // maybe add a tiny bias to prevent z-fighting
	auto objs_end(interior->room_geom->objs.end()); // must include stairs and elevators, which means we have to check all objects
	bool is_placed_on_obj(0);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (can_place_on_object(*i, pos, radius, z_bias, start_zval, zval)) { // can place on this object
			is_placed_on_obj = 1; // Note: any object placed on top of this object will be added/seen after it
			continue;
		}
		if (is_placed_on_obj && is_blocking_obj_on_top_surface(*i)) {} // check it; skips the code below
		else if (!i->is_floor_collidable())          continue; // ignore
		if (!sphere_cube_intersect(pos, radius, *i)) continue; // no intersection
		if (i->is_round() && !dist_xy_less_than(pos, i->get_cube_center(), (i->get_radius() + radius))) continue; // round object (approx)
		return 0; // object in the way, can't place here
	} // for i
	for (auto i = interior->room_geom->expanded_objs.begin(); i != interior->room_geom->expanded_objs.end(); ++i) { // check books, etc.
		if (!i->can_place_onto() || !i->contains_pt_xy(pos) || i->z2() < zval || i->z2() > start_zval) continue; // not a valid placement
		zval = i->z2() + z_bias; // place on top of this object
	}
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) { // check for door intersection
		if (!i->open && sphere_cube_intersect(pos, radius, i->get_true_bcube())) return 0; // blocked by closed door (open doors are more difficult)
	}
	return 1;
}

