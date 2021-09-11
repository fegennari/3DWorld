// 3D World - Building vs. Player/AI Interaction Logic
// by Frank Gennari 1/14/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "openal_wrap.h"

// physics constants, currently applied to balls
float const KICK_VELOCITY  = 0.0025;
float const THROW_VELOCITY = 0.0050;
float const MIN_VELOCITY   = 0.0001;
float const OBJ_DECELERATE = 0.008;
float const OBJ_GRAVITY    = 0.0003;
float const TERM_VELOCITY  = 1.0;
float const OBJ_ELASTICITY = 0.8;
float const ALERT_THRESH   = 0.08; // min sound alert level for AIs

bool do_room_obj_pickup(0), use_last_pickup_object(0), show_bldg_pickup_crosshair(0), player_near_toilet(0), player_in_elevator(0), city_action_key(0);
int can_pickup_bldg_obj(0);
float office_chair_rot_rate(0.0), cur_building_sound_level(0.0);
carried_item_t player_held_object;
bldg_obj_type_t bldg_obj_types[NUM_ROBJ_TYPES];
vector<sphere_t> cur_sounds; // radius = sound volume

extern bool building_action_key, camera_in_building, tt_fire_button_down, flashlight_on, player_is_hiding;
extern int window_width, window_height, display_framerate, player_in_closet, frame_counter, display_mode, game_mode, animate2, camera_surf_collide;
extern float fticks, CAMERA_RADIUS;
extern double tfticks, camera_zh;
extern building_params_t global_building_params;
extern building_dest_t cur_player_building_loc;


void place_player_at_xy(float xval, float yval);
room_object_t get_dresser_middle(room_object_t const &c);
room_object_t get_desk_drawers_part(room_object_t const &c);
cube_t get_sink_cube(room_object_t const &c);
bool player_can_unlock_door();
void show_key_icon();
bool player_has_room_key();
void register_broken_object(room_object_t const &obj);

bool in_building_gameplay_mode() {return (game_mode == 2);} // replaces dodgeball mode

// Note: pos is in camera space
void gen_sound_thread_safe(unsigned id, point const &pos, float gain=1.0, float pitch=1.0, float gain_scale=1.0, bool skip_if_already_playing=0) {
	float const dist(p2p_dist(get_camera_pos(), pos)), dscale(10.0*CAMERA_RADIUS*gain_scale); // distance at which volume is halved
	gain *= dscale/(dist + dscale);
	if (gain < 0.025) return; // too soft to hear
#pragma omp critical(gen_sound)
	gen_sound(id, pos, gain, pitch, 0, zero_vector, skip_if_already_playing);
}
void gen_sound_thread_safe_at_player(unsigned id, float gain=1.0, float pitch=1.0) {
	gen_sound_thread_safe(id, get_camera_pos(), gain, pitch);
}

// lights

float get_radius_for_room_light(room_object_t const &obj);
void register_building_sound(point const &pos, float volume);
void register_building_sound_at_player(float volume);

// Note: called by the player; closest_to is in building space, not camera space
bool building_t::toggle_room_light(point const &closest_to, bool sound_from_closest_to, int room_id) {
	if (!has_room_geom()) return 0; // error?

	if (room_id < 0) { // caller has not provided a valid room_id, so determine it now
		point query_pt(closest_to);
		if (is_rotated()) {do_xy_rotate_inv(bcube.get_cube_center(), query_pt);}
		room_id = get_room_containing_pt(query_pt);
		if (room_id < 0) return 0; // closest_to is not contained in a room of this building
	}
	room_t const &room(get_room(room_id));
	vector<room_object_t> &objs(interior->room_geom->objs);
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
	float closest_dist_sq(0.0);
	int closest_light(-1);

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_type() || i->room_id != room_id) continue; // not a light, or the wrong room
		if (bool(i->flags & RO_FLAG_IN_CLOSET) != bool(player_in_closet)) continue; // while in the closet, player can only toggle closet lights and not room lights
		if (i->flags & RO_FLAG_IN_ELEV) continue; // can't toggle elevator light
		if (!room.is_sec_bldg && get_floor_for_zval(i->z1()) != get_floor_for_zval(closest_to.z)) continue; // wrong floor (skip garages and sheds)
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

void building_t::toggle_light_object(room_object_t const &light, point const &sound_pos) {
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
	bool updated(0);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) { // toggle all lights on this floor of this room
		if (i->is_light_type() && i->room_id == light.room_id && i->z1() == light.z1()) { // Note: closet light should have a different z1 that should not match the room lights
			i->toggle_lit_state(); // Note: doesn't update indir lighting
			if (i->type == TYPE_LAMP) continue; // lamps don't affect room object ambient lighting, and don't require regenerating the vertex data, so skip the step below
			set_obj_lit_state_to(light.room_id, light.z2(), i->is_lit()); // update object lighting flags as well
			updated = 1;
		}
	} // for i
	if (updated) {interior->room_geom->clear_and_recreate_lights();} // recreate light geom with correct emissive properties
	gen_sound_thread_safe(SOUND_CLICK, local_to_camera_space(sound_pos));
	register_building_sound(sound_pos, 0.1);
	//interior->room_geom->modified_by_player = 1; // should light state always be preserved?
}

bool building_t::set_room_light_state_to(room_t const &room, float zval, bool make_on) { // called by AI people
	if (!has_room_geom()) return 0; // error?
	if (room.is_hallway)  return 0; // don't toggle lights for hallways, which can have more than one light
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
	float const window_vspacing(get_window_vspace());
	bool updated(0);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->type != TYPE_LIGHT) continue; // not a light (excludes lamps)
		if (i->flags & (RO_FLAG_IN_CLOSET | RO_FLAG_IN_ELEV)) continue; // skip lights in closets or elevators, these can only be toggled by the player
		if (i->z1() < zval || i->z1() > (zval + window_vspacing) || !room.contains_cube_xy(*i)) continue; // light is on the wrong floor or in the wrong room
		if (i->is_lit() != make_on) {i->toggle_lit_state(); updated = 1;} // Note: doesn't update indir lighting or room light value
	} // for i
	if (updated) {interior->room_geom->clear_and_recreate_lights();} // recreate light geom with correct emissive properties; will flag for update next frame
	return updated;
}

void building_t::set_obj_lit_state_to(unsigned room_id, float light_z2, bool lit_state) {
	assert(has_room_geom());
	float const light_intensity(get_room(room_id).light_intensity);
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
	float const obj_zmin(light_z2 - get_window_vspace()); // get_floor_thickness()?
	bool was_updated(0);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->room_id != room_id || i->z1() < obj_zmin || i->z1() > light_z2) continue; // wrong room or floor
		if (i->is_obj_model_type()) continue; // light_amt currently does not apply to 3D models; should it?

		if (i->type == TYPE_STAIR || i->type == TYPE_STAIR_WALL || i->type == TYPE_ELEVATOR || i->type == TYPE_LIGHT || i->type == TYPE_BLOCKER ||
			i->type == TYPE_COLLIDER || i->type == TYPE_SIGN || i->type == TYPE_WALL_TRIM || i->type == TYPE_RAILING || i->type == TYPE_BLINDS)
		{
			continue; // not a type that uses light_amt
		}
		else if (i->type == TYPE_WINDOW) {
			if (lit_state) {i->flags |= RO_FLAG_LIT;} else {i->flags &= ~RO_FLAG_LIT;}
		}
		else {
			if (lit_state) {i->light_amt += light_intensity;} else {i->light_amt = max((i->light_amt - light_intensity), 0.0f);} // shouldn't be negative, but clamp to 0 just in case
		}
		was_updated = 1;
	} // for i
	if (was_updated) {interior->room_geom->materials_invalid = 1;} // need to recreate them; can't clear here if called from building AI (not in the draw thread)
}

bool building_room_geom_t::closet_light_is_on(cube_t const &closet) const {
	auto objs_end(get_std_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type == TYPE_LIGHT && (i->flags & RO_FLAG_IN_CLOSET) && closet.contains_cube(*i)) {return i->is_lit();}
	}
	return 0;
}

// doors and other interactive objects

// used for drawing open doors
int building_t::find_ext_door_close_to_point(tquad_with_ix_t &door, point const &pos, float dist) const {
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

// used for pedestrians; pos should be outside the building
bool building_t::get_building_door_pos_closest_to(point const &target_pos, point &door_pos) const {
	float dmin_sq(0.0);

	for (auto d = doors.begin(); d != doors.end(); ++d) {
		if (d->type == tquad_with_ix_t::TYPE_GDOOR || d->type == tquad_with_ix_t::TYPE_RDOOR) continue; // skip garage and rooftop doors
		point const center(d->get_bcube().get_cube_center());
		float const dsq(p2p_dist_xy_sq(target_pos, center)); // ignore zval
		if (dmin_sq == 0.0 || dsq < dmin_sq) {door_pos = center; dmin_sq = dsq;}
	}
	if (dmin_sq == 0.0) return 0; // doors not added for some reason
	return 1;
}

void building_t::register_open_ext_door_state(int door_ix) {
	bool const is_open(door_ix >= 0), was_open(open_door_ix >= 0);
	if (is_open == was_open) return; // no state change
	unsigned const dix(is_open ? (unsigned)door_ix : (unsigned)open_door_ix);
	assert(dix < doors.size());
	auto const &door(doors[dix]);
	vector3d const normal(door.get_norm());
	bool const dim(fabs(normal.x) < fabs(normal.y)), dir(normal[dim] < 0.0);
	point const door_center(door.get_bcube().get_cube_center()), sound_pos(local_to_camera_space(door_center)); // convert to camera space
	gen_sound_thread_safe((is_open ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_DOOR_CLOSE), sound_pos);
	point pos_interior(door_center);
	pos_interior[dim] += (dir ? 1.0 : -1.0)*CAMERA_RADIUS; // move point to the building interior so that it's a valid AI position
	register_building_sound(pos_interior, 0.4); // slightly quieter than interior doors because the user has no control over this
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

bool can_open_bathroom_stall(room_object_t const &stall, point const &pos, vector3d const &dir) {
	point door_center;
	door_center[ stall.dim] = stall.d[stall.dim][!stall.dir];
	door_center[!stall.dim] = stall.get_center_dim(!stall.dim);
	door_center.z = pos.z;
	return (dot_product_ptv(dir, door_center, pos) > 0.0); // facing the stall door
}

bool building_t::apply_player_action_key(point const &closest_to_in, vector3d const &in_dir_in) { // called for the player
	if (!interior) return 0; // error?
	float const dmax(4.0*CAMERA_RADIUS), floor_spacing(get_window_vspace());
	float closest_dist_sq(0.0);
	unsigned door_ix(0), obj_ix(0);
	bool found_item(0), is_obj(0);
	vector3d in_dir(in_dir_in);
	point closest_to(closest_to_in);
	maybe_inv_rotate_pos_dir(closest_to, in_dir);

	if (!player_in_closet) { // if the player is in the closet, only the closet door can be opened
		for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
			if (i->z1() > closest_to.z || i->z2() < closest_to.z) continue; // wrong floor, skip
			point const center(i->get_cube_center());
			float const dist_sq(p2p_dist_sq(closest_to, center));
			if (found_item && dist_sq >= closest_dist_sq) continue; // not the closest
			if (!check_obj_dir_dist(closest_to, in_dir, *i, center, dmax)) continue; // door is not in the correct direction or too far away, skip
			closest_dist_sq = dist_sq;
			door_ix    = (i - interior->doors.begin());
			found_item = 1;
		} // for i
	}
	if (interior->room_geom) { // check for closet doors in houses, bathroom stalls in office buildings, and other objects that can be interacted with
		vector<room_object_t> &objs(interior->room_geom->objs), &expanded_objs(interior->room_geom->expanded_objs);
		auto objs_end(interior->room_geom->get_stairs_start());
		point const query_ray_end(closest_to + dmax*in_dir);
		cube_t active_area;

		// make a first pass over all the large objects to determinf if the player is inside one; in that case, the player can't reach out and interact with an object outside it
		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->type != TYPE_STALL && i->type != TYPE_SHOWER && i->type != TYPE_ELEVATOR && !(i->type == TYPE_CLOSET && i->is_open())) continue; // TYPE_CUBICLE?
			if (!i->contains_pt(closest_to)) continue;
			active_area = *i;
			break; // there can be only one - done
		}
		for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
			auto const &obj_vect((vect_id == 1) ? expanded_objs : objs);
			unsigned const obj_id_offset((vect_id == 1) ? objs.size() : 0);
			auto obj_vect_end((vect_id == 1) ? expanded_objs.end() : objs_end); // skip stairs and elevators

			for (auto i = obj_vect.begin(); i != obj_vect_end; ++i) {
				if (cur_player_building_loc.room_ix >= 0 && i->room_id != cur_player_building_loc.room_ix && i->type != TYPE_BUTTON) continue; // not in the same room as the player
				if (!active_area.is_all_zeros() && !i->intersects(active_area)) continue; // out of reach for the player
				bool keep(0);
				if (i->type == TYPE_BOX && !i->is_open()) {keep = 1;} // box can only be opened once; check first so that selection works for boxes in closets
				else if (i->type == TYPE_CLOSET && in_dir.z < 0.5) {keep = 1;} // closet door can be opened; not looking up at the light
				else if (!player_in_closet) {
					if      (i->type == TYPE_TOILET || i->type == TYPE_URINAL) {keep = 1;} // toilet/urinal can be flushed
					else if (i->type == TYPE_STALL && i->shape == SHAPE_CUBE && can_open_bathroom_stall(*i, closest_to, in_dir)) {keep = 1;} // cube bathroom stall can be opened
					else if (i->type == TYPE_OFF_CHAIR && (i->flags & RO_FLAG_RAND_ROT)) {keep = 1;} // office chair can be rotated
					else if (i->is_sink_type() || i->type == TYPE_TUB) {keep = 1;} // sink/tub
					else if (i->is_light_type()) {keep = 1;} // room light or lamp
					else if (i->type == TYPE_PICTURE || i->type == TYPE_TPROLL || i->type == TYPE_BUTTON || i->type == TYPE_MWAVE || i->type == TYPE_TV || i->type == TYPE_MONITOR) {keep = 1;}
					else if (i->type == TYPE_BLINDS || i->type == TYPE_SHOWER || i->type == TYPE_SWITCH /*|| i->type == TYPE_BOOK*/) {keep = 1;}
				}
				else if (i->type == TYPE_LIGHT) {keep = 1;} // closet light
				if (!keep) continue;
				cube_t obj_bc(*i);
				if (i->type == TYPE_KSINK || i->type == TYPE_BRSINK) {obj_bc = get_sink_cube(*i);} // the sink itself is actually smaller
				point center;

				if (i->type == TYPE_CLOSET) {
					center = i->get_cube_center();
					center[i->dim] = i->d[i->dim][i->dir]; // use center of door, not center of closet
				}
				else {center = obj_bc.closest_pt(closest_to);}
				if (fabs(center.z - closest_to.z) > 0.7*floor_spacing) continue; // wrong floor
				float const dist_sq((i->type == TYPE_CLOSET) ? dmax*dmax : p2p_dist_sq(closest_to, center)); // use dmax for closets to prioritize objects inside closets
				if (found_item && dist_sq >= closest_dist_sq)          continue; // not the closest
				if (!obj_bc.closest_dist_less_than(closest_to, dmax))  continue; // too far
				if (in_dir != zero_vector && !obj_bc.line_intersects(closest_to, query_ray_end)) continue; // player is not pointing at this object
				closest_dist_sq = dist_sq;
				obj_ix = (i - obj_vect.begin()) + obj_id_offset;
				is_obj = found_item = 1;
			} // for i
		} // for vect_id
		if (!player_in_closet) {
			float const drawer_dist(found_item ? sqrt(closest_dist_sq) : 2.5*CAMERA_RADIUS);
			if (interior->room_geom->open_nearest_drawer(*this, closest_to, in_dir, drawer_dist, 0)) return 0; // drawer is closer - open or close it
		}
		if (!found_item && !player_in_closet) {move_nearest_object(closest_to, in_dir, 3.0*CAMERA_RADIUS);} // try to move an object instead
		if (!found_item) return 0; // no door or object found
	}
	if (is_obj) { // interactive object
		if (!interact_with_object(obj_ix, closest_to, in_dir)) return 0; // generate sound from the player height
	}
	else { // interior door
		door_t &door(interior->doors[door_ix]);
		if (door.is_closed_and_locked() && !player_can_unlock_door()) return 0; // locked
		if (door.locked && !player_has_room_key()) {door.locked = 0;} // don't lock door when closing, to prevent the player from locking themselves in a room
		toggle_door_state(door_ix, 1, 1, closest_to.z); // toggle state if interior door; player_in_this_building=1, by_player=1, at player height
		//interior->room_geom->modified_by_player = 1; // should door state always be preserved?
	}
	return 1;
}

bool building_t::interact_with_object(unsigned obj_ix, point const &int_pos, vector3d const &int_dir) {
	auto &obj(interior->room_geom->get_room_object_by_index(obj_ix));
	float const pitch((obj.type == TYPE_STALL) ? 2.0 : 1.0); // higher pitch for stalls
	point const sound_origin(obj.xc(), obj.yc(), int_pos.z), local_center(local_to_camera_space(sound_origin)); // generate sound from the player height
	float sound_scale(0.5); // for building sound level
	bool update_draw_data(0);

	if (obj.type == TYPE_TOILET || obj.type == TYPE_URINAL) { // toilet/urinal can be flushed, but otherwise is not modified
		gen_sound_thread_safe(SOUND_FLUSH, local_center);
		sound_scale = 0.5;
	}
	else if (obj.is_sink_type() || obj.type == TYPE_TUB) { // sink or tub
		if (!(obj.flags & RO_FLAG_IS_ACTIVE) && obj.type == TYPE_TUB) {gen_sound_thread_safe(SOUND_SINK, local_center);} // play sound when turning the tub on
		if (obj.is_sink_type()) {obj.flags ^= RO_FLAG_IS_ACTIVE;} // toggle active bit, only for sinks for now
		sound_scale = 0.4;
	}
	else if (obj.is_light_type()) {
		toggle_light_object(obj, obj.get_cube_center());
		sound_scale = 0.0; // sound has already been registered above
	}
	else if (obj.type == TYPE_TPROLL) {
		if (!(obj.flags & (RO_FLAG_HANGING | RO_FLAG_WAS_EXP))) {
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
	else if (obj.type == TYPE_MWAVE) { // beeps
		gen_sound_thread_safe(SOUND_BEEP, local_center, 0.25);
		sound_scale = 0.6;
	}
	else if (obj.type == TYPE_TV || obj.type == TYPE_MONITOR) {
		if (!(obj.flags & RO_FLAG_BROKEN)) { // no visual effect if broken, but still clicks
			if (obj.type == TYPE_MONITOR && (obj.obj_id & 1)) {--obj.obj_id;} // toggle on and off, but don't change the desktop
			else {++obj.obj_id;} // toggle on/off, and also change the picture
			update_draw_data = 1;
		}
		gen_sound_thread_safe(SOUND_CLICK, local_center, 0.4);
	}
	else if (obj.type == TYPE_BUTTON) {
		if (!(obj.flags & RO_FLAG_IS_ACTIVE)) { // if not already active
			register_button_event(obj);
			interior->room_geom->clear_static_small_vbos(); // need to regen object data due to lit state change; don't have to set modified_by_player
			obj.flags |= RO_FLAG_IS_ACTIVE;
		}
	}
	else if (obj.type == TYPE_SWITCH) {
		gen_sound_thread_safe_at_player(SOUND_CLICK, 0.5);
		toggle_room_light(obj.get_cube_center(), 1, obj.room_id); // should select the correct light(s) for the room containing the switch
		obj.flags       ^= RO_FLAG_OPEN; // toggle on/off
		sound_scale      = 0.1;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_BLINDS) { // see building_t::add_window_blinds()
		if (!adjust_blinds_state(obj_ix)) return 0;
		gen_sound_thread_safe_at_player(SOUND_SLIDING, 0.5);
		sound_scale      = 0.3;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_BOOK) {
		obj.flags       ^= RO_FLAG_OPEN; // toggle open/closed
		sound_scale      = 0.0; // no sound
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_SHOWER) { // shower
		if (can_open_bathroom_stall(obj, int_pos, int_dir)) { // open/close shower door
			obj.flags ^= RO_FLAG_OPEN; // toggle open/close
			gen_sound_thread_safe_at_player((obj.is_open() ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_METAL_DOOR));
			sound_scale      = 0.35;
			update_draw_data = 1;
		}
		else { // turn on shower water
			gen_sound_thread_safe_at_player(SOUND_SINK);
			sound_scale = 0.5;
		}
	}
	else if (obj.type == TYPE_BOX) {
		gen_sound_thread_safe_at_player(SOUND_OBJ_FALL, 0.5);
		obj.flags       |= RO_FLAG_OPEN; // mark as open
		sound_scale      = 0.2;
		update_draw_data = 1;
	}
	else if (obj.type == TYPE_CLOSET || obj.type == TYPE_STALL || obj.type == TYPE_SHOWER) {
		obj.flags ^= RO_FLAG_OPEN; // toggle open/close

		if (obj.type == TYPE_CLOSET) {
			interior->room_geom->expand_object(obj); // expand any boxes so that the player can pick them up
			sound_scale = 0.25; // closets are quieter, to allow players to more easily hide
		}
		if (obj.is_small_closet()) {play_door_open_close_sound(sound_origin, obj.is_open(), 1.0, pitch);}
		else {gen_sound_thread_safe_at_player(SOUND_SLIDING);}
		update_draw_data = 1;
	}
	else {assert(0);} // unhandled type
	if (update_draw_data) {interior->room_geom->update_draw_state_for_room_object(obj, *this, 0);}
	if (sound_scale > 0.0) {register_building_sound(sound_origin, sound_scale);}
	if (obj.type == TYPE_BOX) {add_box_contents(obj);} // must be done last to avoid reference invalidation
	return 1;
}

bool building_t::adjust_blinds_state(unsigned obj_ix) {
	auto &obj(interior->room_geom->get_room_object_by_index(obj_ix));

	if (obj.flags & RO_FLAG_HANGING) { // hanging horizontal blinds
		float const floor_spacing(get_window_vspace()), window_v_border(0.94*get_window_v_border()); // border_mult=0.94 to account for the frame
		float const window_height(floor_spacing*(1.0 - 2.0*window_v_border)), blinds_height(obj.dz());
		assert(window_height > 0.0);
		bool const mostly_open(blinds_height < 0.5*window_height);
		if (mostly_open) {obj.z1() = obj.z2() - floor_spacing*(1.0 - window_v_border) + 0.05*floor_spacing;} // close the blinds fully
		else             {obj.z1() = obj.z2() - (1.8*get_wall_thickness() + 0.05*floor_spacing);} // open the blinds
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
		float const fixed_end(obj.d[!obj.dim][!move_dir]), width(obj.get_sz_dim(!obj.dim));
		float const window_center(0.5f*(fixed_end + other_blinds.d[!obj.dim][move_dir])); // center of the span of the pair of left/right blinds
		bool const mostly_open(width < 0.5*fabs(fixed_end - window_center));
		float &move_edge(obj.d[!obj.dim][move_dir]);
		if (mostly_open) {move_edge = window_center;} // close the blinds fully
		else             {move_edge = 0.25*window_center + 0.75*fixed_end;} // open the blinds
	}
	assert(obj.is_strictly_normalized());
	return 1;
}

void building_t::toggle_door_state(unsigned door_ix, bool player_in_this_building, bool by_player, float zval) { // called by the player or AI
	assert(interior && door_ix < interior->doors.size());
	door_t &door(interior->doors[door_ix]);
	door.open ^= 1; // toggle open state
	// we changed the door state, but navigation should adapt to this, except for doors on stairs (which are special)
	if (door.on_stairs) {invalidate_nav_graph();} // any in-progress paths may have people walking to and stopping at closed/locked doors
	interior->door_state_updated = 1; // required for AI navigation logic to adjust to this change
	if (has_room_geom()) {interior->room_geom->mats_doors.clear();} // need to recreate doors VBO

	if (player_in_this_building) { // is it really safe to call this from the AI thread?
		point const door_center(door.xc(), door.yc(), zval);
		play_door_open_close_sound(door_center, door.open);
		if (by_player) {register_building_sound(door_center, 0.5);}
	}
}

point building_t::local_to_camera_space(point const &pos) const {
	point pos_rot(pos);
	if (is_rotated()) {do_xy_rotate(bcube.get_cube_center(), pos_rot);}
	return (pos_rot + get_camera_coord_space_xlate()); // convert to camera space
}

void building_t::play_door_open_close_sound(point const &pos, bool open, float gain, float pitch) const {
	gen_sound_thread_safe((open ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_DOOR_CLOSE), local_to_camera_space(pos), gain, pitch);
}

// dynamic objects: elevators and balls

void apply_object_bounce(building_t const &building, vector3d &velocity, vector3d const &cnorm, point const &pos, float hardness, bool on_floor) {
	bool const vert_coll(cnorm == plus_z);
	float const vmag(velocity.mag()), elasticity(OBJ_ELASTICITY*hardness);
	if (vmag < TOLERANCE) return;
	vector3d v_ref;
	calc_reflection_angle(velocity/vmag, v_ref, cnorm);
	velocity = vmag*v_ref;
	if (on_floor && vert_coll) {velocity.z *= elasticity;} else {velocity *= elasticity;} // only attenuate velocity in Z for a floor collision
	float const bounce_volume(min(1.0f, (vert_coll ? 0.5f : 1.0f)*vmag/KICK_VELOCITY)); // relative to kick velocity

	if (bounce_volume > 0.25) { // apply bounce sound
		if (bounce_volume > 0.5) {gen_sound_thread_safe(SOUND_KICK_BALL, building.local_to_camera_space(pos), 0.75*bounce_volume*bounce_volume);}
		register_building_sound(pos, 0.7*bounce_volume);
	}
}

bool check_ball_kick(room_object_t &ball, vector3d &velocity, point &new_center, point const &p_pos, float pz1, float pz2, float pradius) {
	point const center(ball.get_cube_center()), ppos(p_pos.x, p_pos.y, center.z); // use zval of object (Z range was checked above)
	float const radius(ball.get_radius()), r_sum(radius + pradius);
	if (ball.z2() < pz1 || ball.z1() > pz2 || !dist_xy_less_than(ppos, center, r_sum)) return 0; // no collision
	vector3d const dir((center - ppos).get_norm());
	new_center = (center + (1.05*r_sum - p2p_dist_xy(ppos, center))*dir); // move so that it no longer collides with a bit of tolerance
	velocity.x = KICK_VELOCITY*dir.x; velocity.y = KICK_VELOCITY*dir.y; // keep existing velocity.z
	return 1;
}

void building_t::update_player_interact_objects(point const &player_pos, unsigned building_ix, int first_ped_ix) {
	assert(interior);
	interior->update_elevators(*this, player_pos);
	if (!has_room_geom()) return; // nothing else to do
	float const player_radius(get_scaled_player_radius()), player_z1(player_pos.z - camera_zh - player_radius), player_z2(player_pos.z);
	float const fc_thick(0.5*get_floor_thickness()), fticks_stable(min(fticks, 1.0f)); // cap to 1/40s to improve stability
	static float last_sound_tfticks(0);
	static point last_sound_pt(all_zeros), last_player_pos(all_zeros);
	vect_cube_t ped_bcubes;
	if (first_ped_ix >= 0) {get_ped_bcubes_for_building(first_ped_ix, building_ix, ped_bcubes);}
	bool const player_is_moving(player_pos != last_player_pos);
	last_player_pos = player_pos;

	for (auto c = interior->room_geom->objs.begin(); c != interior->room_geom->objs.end(); ++c) { // check for other objects to collide with (including stairs)
		if (c->no_coll() || !c->has_dstate()) continue; // Note: no test of player_coll flag
		assert(c->type == TYPE_LG_BALL); // currently, only large balls have has_dstate()
		assert(c->obj_id < interior->room_geom->obj_dstate.size());
		room_obj_dstate_t &dstate(interior->room_geom->get_dstate(*c));
		vector3d &velocity(dstate.velocity);
		point const center(c->get_cube_center());
		float const radius(c->get_radius());
		bool const was_dynamic(c->is_dynamic());
		bool on_floor(0);
		point new_center(center);
		bool kicked(camera_surf_collide && player_is_moving && check_ball_kick(*c, velocity, new_center, player_pos, player_z1, player_z2, player_radius)); // check the player
		bool const is_moving_fast(velocity.mag() > 0.5*KICK_VELOCITY);

		for (auto p = ped_bcubes.begin(); p != ped_bcubes.end(); ++p) { // check building AI people
			if (is_moving_fast) { // treat collision as a bounce
				vector3d cnorm;

				if (sphere_cube_int_update_pos(new_center, radius, *p, center, 1, 0, &cnorm)) {
					register_person_hit((first_ped_ix + (p - ped_bcubes.begin())), *c, velocity);
					apply_object_bounce(*this, velocity, cnorm, new_center, 0.75, on_floor); // hardness=0.75
				}
			}
			else { // treat collision as a kick
				kicked = check_ball_kick(*c, velocity, new_center, p->get_cube_center(), p->z1(), p->z2(), 0.6*p->dx()); // assume dx == dy == radius
			}
		}
		if (kicked) {
			c->flags |= RO_FLAG_DYNAMIC; // make it dynamic

			if ((tfticks - last_sound_tfticks) > 1.0*TICKS_PER_SECOND && !dist_less_than(new_center, last_sound_pt, radius)) { // play at most once per second
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
					velocity *= (1.0f - min(1.0f, OBJ_DECELERATE*step_sz));
					if (velocity.mag_sq() < MIN_VELOCITY*MIN_VELOCITY) {velocity = zero_vector;} // zero velocity if stopped
				}
				else { // in the air - apply gravity
					velocity.z -= OBJ_GRAVITY*step_sz; // apply gravitational acceleration
					max_eq(velocity.z, -TERM_VELOCITY);
				}
				if (velocity == zero_vector) { // stopped
					c->flags &= ~RO_FLAG_DYNAMIC; // clear dynamic flag
					interior->update_dynamic_draw_data(); // remove from dynamic objects and schedule an update
					interior->room_geom->clear_static_small_vbos(); // add to small static objects
					break; // done
				}
				new_center += velocity*step_sz; // move based on velocity
			} // for step
		}
		if (new_center != center) { // check for collisions and move to new location
			vector3d cnorm;
			int obj_ix(-1);
			float hardness(0.0);

			if (interior->check_sphere_coll(*this, new_center, center, radius, c, cnorm, hardness, obj_ix)) {
				if (cnorm == plus_z) { // collision with the floor or the top surface of something
					if (fabs(velocity.z) < 0.25*OBJ_GRAVITY*fticks) {velocity.z = 0.0;} // zero velocity z component if near zero to reduce instability
				}
				apply_object_bounce(*this, velocity, cnorm, new_center, hardness, on_floor);
				
				if (obj_ix >= 0) { // collided with a room object
					auto &obj(interior->room_geom->get_room_object_by_index(obj_ix));
					bool handled(0);

					// break the glass if not already broken
					if ((obj.type == TYPE_TV || obj.type == TYPE_MONITOR) && velocity.mag() > 2.0*MIN_VELOCITY && !(obj.flags & RO_FLAG_BROKEN)) {
						vector3d front_dir(all_zeros);
						front_dir[obj.dim] = (obj.dir ? 1.0 : -1.0);
						
						if (dot_product(cnorm, front_dir) > 0.9) { // hit the front side of the screen
							if (dist_less_than(center, obj.get_cube_center(), (radius + 0.5*obj.get_sz_dim(obj.dim) + 0.2*obj.dz()))) { // near the screen center
								// capture value before breaking; if the player then takes this object, damage will be higher, but we can attribute this to making a mess of broken glass
								register_broken_object(obj);
								obj.flags |= RO_FLAG_BROKEN;
								point const sound_origin(obj.xc(), obj.yc(), center.z); // generate sound from the player height
								gen_sound_thread_safe(SOUND_GLASS, local_to_camera_space(sound_origin), 0.7);
								register_building_sound(sound_origin, 0.7);
								interior->room_geom->update_draw_state_for_room_object(obj, *this, 0);
								handled = 1;
							}
						}
					}
					if (obj.type == TYPE_PICTURE || obj.type == TYPE_TV || obj.type == TYPE_MONITOR || obj.type == TYPE_BUTTON || obj.type == TYPE_SWITCH ||
						(obj.type == TYPE_OFF_CHAIR && (obj.flags & RO_FLAG_RAND_ROT)))
					{
						if (!handled) {interact_with_object(obj_ix, center, velocity.get_norm());}
					}
				}
			}
			point const prev_new_center(new_center);
			
			if (move_sphere_to_valid_part(new_center, center, radius) && new_center != prev_new_center) { // collision with exterior wall
				apply_object_bounce(*this, velocity, (new_center - prev_new_center).get_norm(), new_center, 1.0, on_floor); // hardness=1.0
				// add TYPE_CRACK if collides with a window?
			}
			if (new_center != center) {apply_roll_to_matrix(dstate.rot_matrix, new_center, center, plus_z, radius, (on_floor ? 0.0 : 0.01), (on_floor ? 1.0 : 0.2));}
			if (!was_dynamic) {interior->room_geom->clear_static_small_vbos();} // static => dynamic transition, need to remove from static object vertex data
			interior->update_dynamic_draw_data();
			c->translate(new_center - center);
		}
	} // for c
	if (use_last_pickup_object || (tt_fire_button_down && !flashlight_on)) { // use object not active, and not using fire key without flashlight (space bar)
		maybe_use_last_pickup_room_object(player_pos);
		use_last_pickup_object = 0; // reset for next frame
	}
}

bool building_interior_t::update_elevators(building_t const &building, point const &player_pos) { // Note: player_pos is in building space
	float const z_space(0.05*building.get_floor_thickness()); // to prevent z-fighting
	float const delta_open_amt(min(1.0f, 2.0f*fticks/TICKS_PER_SECOND)); // 0.5s for full open
	static int prev_move_dir(2); // starts at not-moving
	vector<room_object_t> &objs(room_geom->objs);
	bool was_updated(0), update_ddd(0);

	// Note: the player can only be in one elevator at a time, but they can push the call button for one elevator and get into another, so we have to check all elevators
	for (auto e = elevators.begin(); e != elevators.end(); ++e) { // find containing elevator (optimization + need to know z-range of elevator shaft)
		if (!e->was_called) { // stopped on a floor
			if (e->open_amt > 0.0 && e->open_amt < 1.0) { // partially open - continue to open fully
				e->open_amt = min((e->open_amt + delta_open_amt), 1.0f);
				update_ddd  = 1; // regen verts for open door
			}
			continue;
		}
		assert(e->car_obj_id < objs.size());
		room_object_t &obj(objs[e->car_obj_id]); // elevator car for this elevator
		assert(obj.type == TYPE_ELEVATOR && obj.room_id == (e - elevators.begin())); // sanity check

		if (e->open_amt > 0.0 && e->target_zval != obj.z1()) { // doors not yet closed, and not at target zval
			e->open_amt = max((e->open_amt - delta_open_amt), 0.0f);
			update_ddd  = 1; // regen verts for open door
			continue;
		}
		bool const move_dir(e->target_zval > obj.z1()); // 0=down, 1=up
		assert(e->button_id_start < e->button_id_end && e->button_id_end <= objs.size());
		float dist(min(0.5f*CAMERA_RADIUS, 0.04f*obj.dz()*fticks)*(move_dir ? 1.0 : -1.0)); // clamp to half camera radius to avoid falling through the floor for low framerates
		if (e->was_called && fabs(e->target_zval - obj.z1()) < fabs(dist)) {dist = (e->target_zval - obj.z1());} // move to position
		else if (move_dir) {min_eq(dist, (e->z2() - obj.z2() - z_space));} // going up
		else               {max_eq(dist, (e->z1() - obj.z1() + z_space));} // going down
		update_ddd = 1;

		if (fabs(dist) < 0.001*z_space) { // no movement, at target_zval or top/bottom of elevator shaft (check with a tolerance)
			max_eq(e->open_amt, delta_open_amt); // begin to open if not already open
			e->was_called = 0;
			obj.flags    |= RO_FLAG_OPEN;
			bool was_updated(0);

			for (auto j = objs.begin() + e->button_id_start; j != objs.begin() + e->button_id_end; ++j) { // disable all call buttons for this elevator
				if (j->type == TYPE_BLOCKER) continue; // button was removed?
				assert(j->type == TYPE_BUTTON);
				if (j->flags & RO_FLAG_IS_ACTIVE) {j->flags &= ~RO_FLAG_IS_ACTIVE; was_updated = 1;} // clear active/lit state
			}
			if (was_updated) {room_geom->clear_static_small_vbos();} // need to regen object data due to lit state change
			point const sound_pos(obj.get_cube_center());
			gen_sound_thread_safe(SOUND_BEEP, building.local_to_camera_space(sound_pos), 0.5, 0.75); // lower frequency beep
			register_building_sound(sound_pos, 0.5);
			break;
		}
		obj.translate_dim(2, dist); // translate in Z
		obj.item_flags = uint16_t(floor((obj.zc() - e->z1())/building.get_window_vspace())); // set current floor
		assert(e->light_obj_id < objs.size());
		room_object_t &light(objs[e->light_obj_id]); // light for this elevator

		if (light.type != TYPE_BLOCKER) { // translate light as well, if not taken
			assert(light.type == TYPE_LIGHT);
			light.translate_dim(2, dist); // translate in Z
			room_geom->clear_and_recreate_lights(); // or make elevator lights part of the elevator instead?
		}
		for (auto j = objs.begin() + e->button_id_start; j != objs.begin() + e->button_id_end; ++j) {
			if (j->type == TYPE_BLOCKER) continue; // button was removed?
			assert(j->type == TYPE_BUTTON);
			if (j->flags & RO_FLAG_IN_ELEV) {j->translate_dim(2, dist);} // interior panel button, translate in Z
		}
		if ((int)move_dir != prev_move_dir && obj.contains_pt(player_pos)) { // moving, and player is in the elevator
			gen_sound_thread_safe_at_player(SOUND_SLIDING, 0.75); // play this sound when the elevator starts moving or changes direction
			register_building_sound(player_pos, 0.4);
		}
		prev_move_dir = move_dir;
		was_updated   = 1;
	} // for e
	if (update_ddd) {update_dynamic_draw_data();} // clear dynamic material vertex data (for all elevators) and recreate their VBOs
	if (was_updated) return 1;
	prev_move_dir = 2;
	return 0;
}

void building_t::register_button_event(room_object_t const &button) {
	assert(interior);
	assert(button.room_id < interior->elevators.size()); // here room_id is elevator_id (buttons are only used with elevators)
	unsigned const floor_ix(button.obj_id);
	elevator_t &elevator(interior->elevators[button.room_id]);
	elevator.call_elevator(elevator.z1() + max(get_window_vspace()*floor_ix, 0.05f*get_floor_thickness())); // bottom of elevator car for this floor
}

bool building_t::move_sphere_to_valid_part(point &pos, point const &p_last, float radius) const { // Note: only moves in XY
	point const init_pos(pos);
	float xy_area_contained(0.0);
	cube_t sphere_bcube; sphere_bcube.set_from_sphere(pos, radius);
	cube_t valid_region(bcube);
	valid_region.expand_by(-radius);
	valid_region.clamp_pt(pos); // keep pos within the valid building bcube

	for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
		if (!i->intersects(sphere_bcube)) continue;
		cube_t overlap(sphere_bcube);
		overlap.intersect_with_cube(*i);
		xy_area_contained += overlap.dx()*overlap.dy();
	}
	if (xy_area_contained > 0.99*sphere_bcube.dx()*sphere_bcube.dy()) return (pos != init_pos); // sphere contained in union of parts (not outside the building)

	// find part containing p_last and clamp to that part
	for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
		if (!i->contains_pt(p_last)) continue; // not the part containing the previous pos
		cube_t bounds(*i);
		bounds.expand_by_xy(-radius);
		bounds.clamp_pt_xy(pos);
		return 1;
	}
	//assert(0); // should never get here (except for FP error?)
	pos = p_last; // restore original position
	return 1;
}

// ray queries

// center, -x, +x, -y, +y, -z, +z
void get_sphere_boundary_pts(point const &center, float radius, point pts[7]) {
	pts[0] = center;

	for (unsigned dim = 0, ix = 1; dim < 3; ++dim) {
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
bool building_t::check_line_intersect_doors(point const &p1, point const &p2) const {
	float const wall_thickness(get_wall_thickness());

	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
		if (i->open) continue; // check only closed doors
		cube_t door(*i);
		door.expand_in_dim(i->dim, 0.5*wall_thickness); // increase door thickness
		if (door.line_intersects(p1, p2)) return 1;
	}
	return 0;
}
bool building_t::is_pt_visible(point const &p1, point const &p2) const {
	if (!interior) return 1;
	if (is_light_occluded(p1, p2)) return 0; // okay to call for non-light point; checks walls, ceilings, and floors
	if (parts.size() > 1 && !trace_ray_through_cubes(parts, p1, p2, 0.01*get_wall_thickness())) return 0; // view blocked by exterior wall (ignores windows)
	return !check_line_intersect_doors(p1, p2);
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
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (!i->is_light_type() || !i->is_lit()) continue; // not a light, or light not on
		bool const same_room((int)i->room_id == room_id);
		//if (!same_room) continue; // different room (optimization); too strong?
		//bool const same_floor(fabs(i->z1() - pt.z) < floor_spacing); // doesn't work with lamps
		bool const same_floor(room.is_sec_bldg || get_floor_for_zval(pt.z) == get_floor_for_zval(i->z1()));
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

// object types/pickup

void setup_bldg_obj_types() {
	static bool was_setup(0);
	if (was_setup) return; // nothing to do
	was_setup = 1;
	// player_coll, ai_coll, pickup, attached, is_model, lg_sm, value, weight, name [capacity]
	//                                                pc ac pu at im ls value  weight  name
	bldg_obj_types[TYPE_TABLE     ] = bldg_obj_type_t(1, 1, 1, 0, 0, 1, 70.0,  40.0,  "table");
	bldg_obj_types[TYPE_CHAIR     ] = bldg_obj_type_t(0, 1, 1, 0, 0, 1, 50.0,  25.0,  "chair"); // skip player collisions because they can be in the way and block the path in some rooms
	bldg_obj_types[TYPE_STAIR     ] = bldg_obj_type_t(1, 0, 0, 1, 0, 1, 0.0,   0.0,   "stair");
	bldg_obj_types[TYPE_STAIR_WALL] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "stairs wall");
	bldg_obj_types[TYPE_ELEVATOR  ] = bldg_obj_type_t(1, 1, 0, 1, 0, 0, 0.0,   0.0,   "elevator");
	bldg_obj_types[TYPE_LIGHT     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 0, 40.0,  5.0,   "light");
	bldg_obj_types[TYPE_RUG       ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 50.0,  20.0,  "rug");
	bldg_obj_types[TYPE_PICTURE   ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 100.0, 1.0,   "picture"); // should be random value
	bldg_obj_types[TYPE_WBOARD    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 50.0,  25.0,  "whiteboard");
	bldg_obj_types[TYPE_BOOK      ] = bldg_obj_type_t(0, 0, 1, 0, 0, 3, 10.0,  1.0,   "book");
	bldg_obj_types[TYPE_BCASE     ] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 150.0, 100.0, "bookcase"); // Note: can't pick up until bookcase can be expanded and books taken off
	bldg_obj_types[TYPE_TCAN      ] = bldg_obj_type_t(0, 1, 1, 0, 0, 2, 12.0,  2.0,   "trashcan"); // skip player collisions because they can be in the way and block the path in some rooms
	bldg_obj_types[TYPE_DESK      ] = bldg_obj_type_t(1, 1, 0, 0, 0, 3, 100.0, 80.0,  "desk"); // drawers are small items
	bldg_obj_types[TYPE_BED       ] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 300.0, 200.0, "bed"); // pillows are small, and the rest is large
	bldg_obj_types[TYPE_WINDOW    ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 0.0,   0.0,   "window");
	bldg_obj_types[TYPE_BLOCKER   ] = bldg_obj_type_t(0, 0, 0, 0, 0, 0, 0.0,   0.0,   "<blocker>");  // not a drawn object; block other objects, but not the player or AI
	bldg_obj_types[TYPE_COLLIDER  ] = bldg_obj_type_t(1, 1, 0, 0, 0, 0, 0.0,   0.0,   "<collider>"); // not a drawn object; block the player and AI
	bldg_obj_types[TYPE_CUBICLE   ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 500.0, 250.0, "cubicle"); // skip collisions because they have their own colliders
	bldg_obj_types[TYPE_STALL     ] = bldg_obj_type_t(1, 1, 1, 1, 0, 1, 40.0,  20.0,  "bathroom divider"); // can pick up short sections of bathroom stalls (urinal dividers)
	bldg_obj_types[TYPE_SIGN      ] = bldg_obj_type_t(0, 0, 1, 0, 0, 3, 10.0,  1.0,   "sign");
	bldg_obj_types[TYPE_COUNTER   ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "kitchen counter");
	bldg_obj_types[TYPE_CABINET   ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 0.0,   0.0,   "kitchen cabinet");
	bldg_obj_types[TYPE_KSINK     ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "kitchen sink");
	bldg_obj_types[TYPE_BRSINK    ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "bathroom sink"); // for office building bathrooms
	bldg_obj_types[TYPE_PLANT     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 3, 18.0,  8.0,   "potted plant");
	bldg_obj_types[TYPE_DRESSER   ] = bldg_obj_type_t(1, 1, 0, 0, 0, 3, 120.0, 110.0, "dresser"); // Note: can't pick up until drawers can be opened and items removed from them
	bldg_obj_types[TYPE_NIGHTSTAND] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 60.0,  45.0,  "nightstand");
	bldg_obj_types[TYPE_FLOORING  ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 0.0,   0.0,   "flooring");
	// closets can't be picked up, but they can block a pickup; marked as large because small objects are not modified; marked as is_model because closets can contain lamps
	bldg_obj_types[TYPE_CLOSET    ] = bldg_obj_type_t(1, 1, 1, 1, 1, 1, 0.0,   0.0,   "closet");
	bldg_obj_types[TYPE_WALL_TRIM ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "wall trim");
	bldg_obj_types[TYPE_RAILING   ] = bldg_obj_type_t(1, 0, 0, 1, 0, 2, 0.0,   0.0,   "railing");
	bldg_obj_types[TYPE_CRATE     ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 10.0,  12.0,  "crate"); // should be random value
	bldg_obj_types[TYPE_BOX       ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 5.0,   8.0,   "box");   // should be random value
	bldg_obj_types[TYPE_MIRROR    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 40.0,  15.0,  "mirror");
	bldg_obj_types[TYPE_SHELVES   ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 0.0,   0.0,   "shelves");
	bldg_obj_types[TYPE_KEYBOARD  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 15.0,  2.0,   "keyboard");
	bldg_obj_types[TYPE_SHOWER    ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   0.0,   "shower");
	bldg_obj_types[TYPE_RDESK     ] = bldg_obj_type_t(1, 1, 0, 0, 0, 1, 800.0, 300.0, "reception desk");
	bldg_obj_types[TYPE_BOTTLE    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 1.0,   1.0,   "bottle");
	bldg_obj_types[TYPE_WINE_RACK ] = bldg_obj_type_t(1, 1, 1, 0, 0, 3, 75.0,  40.0,  "wine rack");
	bldg_obj_types[TYPE_COMPUTER  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 500.0, 20.0,  "computer");
	bldg_obj_types[TYPE_MWAVE     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 1, 100.0, 50.0,  "microwave oven");
	bldg_obj_types[TYPE_PAPER     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.0,   0.0,   "sheet of paper"); // will have a random value that's often 0
	bldg_obj_types[TYPE_BLINDS    ] = bldg_obj_type_t(0, 0, 0, 1, 0, 1, 50.0,  7.0,   "window blinds");
	bldg_obj_types[TYPE_PEN       ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.10,  0.02,  "pen");
	bldg_obj_types[TYPE_PENCIL    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.10,  0.02,  "pencil");
	bldg_obj_types[TYPE_PAINTCAN  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 12.0,  8.0,   "paint can");
	bldg_obj_types[TYPE_LG_BALL   ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 15.0,  1.2,   "ball");
	bldg_obj_types[TYPE_HANGER_ROD] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 10.0,  5.0,   "hanger rod");
	bldg_obj_types[TYPE_DRAIN     ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "drain pipe");
	bldg_obj_types[TYPE_MONEY     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 20.0,  0.0,   "pile of money"); // $20 bills
	bldg_obj_types[TYPE_PHONE     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 200.0, 0.1,   "cell phone");
	bldg_obj_types[TYPE_TPROLL    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.25,  0.1,   "TP roll", 200);
	bldg_obj_types[TYPE_SPRAYCAN  ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 2.0,   1.0,   "spray paint", 5000);
	bldg_obj_types[TYPE_MARKER    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.20,  0.05,  "marker",      10000);
	bldg_obj_types[TYPE_BUTTON    ] = bldg_obj_type_t(0, 0, 1, 1, 0, 2, 1.0,   0.05,  "button");
	bldg_obj_types[TYPE_CRACK     ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "crack");
	bldg_obj_types[TYPE_SWITCH    ] = bldg_obj_type_t(0, 0, 0, 1, 0, 2, 0.0,   0.0,   "switch");
	bldg_obj_types[TYPE_PLATE     ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 6.0,   0.25,  "plate");
	bldg_obj_types[TYPE_LAPTOP    ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 600.0, 8.0,   "laptop");
	bldg_obj_types[TYPE_FPLACE    ] = bldg_obj_type_t(1, 1, 0, 1, 0, 1, 0.0,   2000.0,"fireplace");
	bldg_obj_types[TYPE_LBASKET   ] = bldg_obj_type_t(1, 1, 1, 0, 0, 2, 12.0,  2.0,   "laundry basket");
	bldg_obj_types[TYPE_WHEATER   ] = bldg_obj_type_t(1, 1, 0, 1, 0, 2, 300.0, 500.0, "water heater");
	// player_coll, ai_coll, pickup, attached, is_model, lg_sm, value, weight, name [capacity]
	// 3D models
	bldg_obj_types[TYPE_TOILET    ] = bldg_obj_type_t(1, 1, 1, 1, 1, 0, 120.0, 88.0,  "toilet");
	bldg_obj_types[TYPE_SINK      ] = bldg_obj_type_t(1, 1, 1, 1, 1, 0, 80.0,  55.0,  "sink");
	bldg_obj_types[TYPE_TUB       ] = bldg_obj_type_t(1, 1, 0, 1, 1, 1, 250.0, 200.0, "bathtub");
	bldg_obj_types[TYPE_FRIDGE    ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 700.0, 300.0, "refrigerator"); // no pickup, too large and may want to keep it for future hunger bar
	bldg_obj_types[TYPE_STOVE     ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 400.0, 150.0, "stove");
	bldg_obj_types[TYPE_TV        ] = bldg_obj_type_t(1, 1, 1, 0, 1, 1, 400.0, 70.0,  "TV");
	bldg_obj_types[TYPE_MONITOR   ] = bldg_obj_type_t(1, 1, 1, 0, 1, 1, 250.0, 15.0,  "computer monitor");
	bldg_obj_types[TYPE_COUCH     ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 600.0, 300.0, "couch");
	bldg_obj_types[TYPE_OFF_CHAIR ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 150.0, 60.0,  "office chair");
	bldg_obj_types[TYPE_URINAL    ] = bldg_obj_type_t(1, 1, 1, 1, 1, 0, 100.0, 80.0,  "urinal");
	bldg_obj_types[TYPE_LAMP      ] = bldg_obj_type_t(0, 0, 1, 0, 1, 0, 25.0,  12.0,  "lamp");
	bldg_obj_types[TYPE_WASHER    ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 300.0, 150.0, "washer");
	bldg_obj_types[TYPE_DRYER     ] = bldg_obj_type_t(1, 1, 1, 0, 1, 0, 300.0, 160.0, "dryer");
	// keys are special because they're potentially either a small object or an object model (in a drawer)
	bldg_obj_types[TYPE_KEY       ] = bldg_obj_type_t(0, 0, 1, 0, 0, 2, 0.0,   0.05,  "room key"); // drawn as an object, not a model
	bldg_obj_types[TYPE_HANGER    ] = bldg_obj_type_t(0, 0, 1, 0, 1, 2, 0.25,  0.05,  "clothes hanger");
	//                                                pc ac pu at im ls value  weight  name [capacity]
}

float carried_item_t::get_remaining_capacity_ratio() const {
	assert(type >= 0 && type < NUM_ROBJ_TYPES);
	unsigned const capacity(bldg_obj_types[type].capacity);
	return ((capacity == 0) ? 1.0 : (1.0 - float(use_count)/float(capacity))); // Note: zero capacity is unlimited and ratio returned is always 1.0
}

bldg_obj_type_t const &get_room_obj_type(room_object_t const &obj) {
	assert(obj.type < NUM_ROBJ_TYPES);
	return bldg_obj_types[obj.type];
}
bldg_obj_type_t get_taken_obj_type(room_object_t const &obj) {
	if (obj.type == TYPE_PICTURE && (obj.flags & RO_FLAG_TAKEN1)) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 20.0, 6.0, "picture frame");} // second item to take from picture
	if (obj.type == TYPE_TPROLL  && (obj.flags & RO_FLAG_TAKEN1)) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 6.0,  0.5, "toilet paper holder");} // second item to take from tproll

	if (obj.type == TYPE_BED) { // player_coll, ai_coll, pickup, attached, is_model, lg_sm, value, weight, name
		if (obj.flags & RO_FLAG_TAKEN2) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 250.0, 80.0, "mattress"  );} // third item to take from bed
		if (obj.flags & RO_FLAG_TAKEN1) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 80.0,  4.0,  "bed sheets");} // second item to take from bed
		return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 20.0, 1.0, "pillow"); // first item to take from bed
	}
	if (obj.type == TYPE_PLANT && !(obj.flags & RO_FLAG_ADJ_BOT)) { // plant not on a table/desk
		if (obj.flags & RO_FLAG_TAKEN2) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 10.0, 10.0, "plant pot");} // third item to take
		if (obj.flags & RO_FLAG_TAKEN1) {return bldg_obj_type_t(0, 0, 1, 0, 0, 1, 1.0,  10.0, "dirt"     );} // second item to take
		return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 25.0, 5.0, "plant"); // first item to take
	}
	if (obj.type == TYPE_COMPUTER && (obj.flags & RO_FLAG_WAS_EXP)) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2, 100.0, 20.0, "old computer");}
	if (obj.type == TYPE_BOX      && (obj.flags & RO_FLAG_OPEN   )) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2,   0.0, 0.05, "opened box"  );}
	if (obj.type == TYPE_CRATE    && (obj.flags & RO_FLAG_OPEN   )) {return bldg_obj_type_t(0, 0, 1, 0, 0, 2,   2.0, 0.5,  "opened crate");}
	if (obj.type == TYPE_TV       && (obj.flags & RO_FLAG_BROKEN )) {return bldg_obj_type_t(1, 1, 1, 0, 1, 1,  20.0, 70.0, "broken TV"   );}
	if (obj.type == TYPE_MONITOR  && (obj.flags & RO_FLAG_BROKEN )) {return bldg_obj_type_t(1, 1, 1, 0, 1, 1,  10.0, 15.0, "broken computer monitor");}

	if (obj.type == TYPE_LG_BALL) {
		bldg_obj_type_t type(get_room_obj_type(obj));
		type.name = ((obj.item_flags & 1) ? "basketball" : "soccer ball"); // use a more specific type name; all other fields are shared across balls
		return type;
	}
	if (obj.type == TYPE_BOTTLE) {
		bottle_params_t const &bparams(bottle_params[obj.get_bottle_type()]);
		bldg_obj_type_t type(0, 0, 1, 0, 0, 2,  bparams.value, 1.0, bparams.name);

		if (obj.is_bottle_empty()) {
			type.name    = "empty " + type.name;
			type.weight *= 0.25;
			type.value   = 0.0;
		}
		return type;
	}
	return get_room_obj_type(obj); // default value
}
rand_gen_t rgen_from_obj(room_object_t const &obj) {
	rand_gen_t rgen;
	rgen.set_state(12345*abs(obj.x1()), 67890*abs(obj.y1()));
	return rgen;
}
float get_obj_value(room_object_t const &obj) {
	float value(get_taken_obj_type(obj).value);
	if (obj.type == TYPE_CRATE || obj.type == TYPE_BOX) {value *= (1 + (rgen_from_obj(obj).rand() % 20));}
	else if (obj.type == TYPE_PAPER) {
		rand_gen_t rgen(rgen_from_obj(obj));
		if (rgen.rand_float() < 0.25) { // 25% of papers have some value
			float const val_mult((rgen.rand_float() < 0.25) ? 10 : 1); // 25% of papers have higher value
			value = val_mult*(2 + (rgen.rand()%10))*(1 + (rgen.rand()%10));
		}
	}
	else if (obj.type == TYPE_MONEY) {
		unsigned const num_bills(round_fp(obj.dz()/(0.01*obj.get_sz_dim(obj.dim))));
		value *= num_bills;
	}
	if (obj.flags & RO_FLAG_USED) {value = 0.01*floor(50.0*value);} // used objects have half value, rounded down to the nearest cent
	return value;
}
float get_obj_weight(room_object_t const &obj) {
	return get_taken_obj_type(obj).weight; // constant per object type, for now, but really should depend on object size/volume
}
bool is_consumable(room_object_t const &obj) {
	return (in_building_gameplay_mode() && obj.type == TYPE_BOTTLE && !obj.is_bottle_empty() && !(obj.flags & RO_FLAG_NO_CONS));
}

void show_weight_limit_message() {
	std::ostringstream oss;
	oss << "Over weight limit of " << global_building_params.player_weight_limit << " lbs";
	print_text_onscreen(oss.str(), RED, 1.0, 1.5*TICKS_PER_SECOND, 0);
}

class phone_manager_t {
	bool is_enabled=0, is_ringing = 0, is_on=0;
	double stop_ring_time=0.0, next_ring_time=0.0, next_cycle_time=0.0, auto_off_time=0.0, next_button_time=0.0;
	rand_gen_t rgen;

	void schedule_next_ring() {next_ring_time = tfticks + (double)rgen.rand_uniform(30.0, 120.0)*TICKS_PER_SECOND;} // 30s to 2min
public:
	bool is_phone_ringing() const {return is_ringing;}
	bool is_phone_on     () const {return is_on     ;}

	void next_frame() {
		if (!is_enabled || !camera_in_building) {} // do nothing
		else if (is_ringing) {
			if (tfticks > stop_ring_time) {is_ringing = 0; schedule_next_ring();} // stop automatically
			else if (tfticks > next_cycle_time) { // start a new ring cycle
				gen_sound_thread_safe_at_player(SOUND_PHONE_RING, 1.0);
				register_building_sound_at_player(1.0);
				next_cycle_time += 4.2*TICKS_PER_SECOND; // 4.2s between rings
			}
		}
		else {
			if (tfticks > next_ring_time) { // start a new ring cycle
				is_ringing      = 1;
				stop_ring_time  = tfticks + (double)rgen.rand_uniform(12.0, 24.0)*TICKS_PER_SECOND; // 10-20s into the future
				next_cycle_time = tfticks; // cycle begins now
			}
			if (is_on && tfticks > auto_off_time) {is_on = 0;} // auto off
		}
	}
	void enable() {
		if (is_enabled) return; // already enabled
		is_enabled = 1;
		is_ringing = is_on = 0;
		schedule_next_ring();
	}
	void disable() {
		is_enabled = is_ringing = is_on = 0;
	}
	void player_action() {
		if (tfticks < next_button_time) return; // skip if pressed immediately after the last press (switch debouncer)
		next_button_time = tfticks + 0.25*TICKS_PER_SECOND;

		if (is_ringing) { // switch off
			is_ringing     = is_on = 0;
			stop_ring_time = tfticks; // now
			schedule_next_ring();
		}
		else {
			is_on ^= 1;
			if (is_on) {auto_off_time = tfticks + 4.0*TICKS_PER_SECOND;} // 4s auto off delay
		}
		gen_sound_thread_safe_at_player(SOUND_CLICK, 0.5);
	}
};

phone_manager_t phone_manager;

class player_inventory_t { // manages player inventory, health, and other stats
	vector<carried_item_t> carried; // interactive items the player is currently carrying
	float cur_value, cur_weight, tot_value, tot_weight, damage_done, best_value, player_health, drunkenness, bladder, bladder_time, prev_player_zval;
	bool prev_in_building, has_key;

	void register_player_death(unsigned sound_id, std::string const &why) {
		point const xlate(get_camera_coord_space_xlate());
		place_player_at_xy(xlate.x, xlate.y); // move back to the origin/spawn location
		gen_sound_thread_safe_at_player(sound_id);
		print_text_onscreen(("You Have Died" + why), RED, 2.0, 2*TICKS_PER_SECOND, 10);
		clear(); // respawn
	}
public:
	player_inventory_t() {clear_all();}

	void clear() { // called on player death
		max_eq(best_value, tot_value);
		cur_value     = cur_weight = tot_value = tot_weight = damage_done = 0.0;
		drunkenness   = bladder = bladder_time = prev_player_zval = 0.0;
		player_health = 1.0; // full health
		prev_in_building = has_key = 0;
		phone_manager.disable();
		carried.clear();
	}
	void clear_all() { // called on game mode init
		tot_value = best_value = 0.0;
		clear();
	}
	void take_damage(float amt) {player_health -= amt*(1.0f - 0.75f*min(drunkenness, 1.0f));} // up to 75% damage reduction when drunk
	void record_damage_done(float amt) {damage_done += amt;}
	void return_object_to_building(room_object_t const &obj) {damage_done -= get_obj_value(obj);}
	bool check_weight_limit(float weight) const {return ((cur_weight + weight) <= global_building_params.player_weight_limit);}
	bool can_pick_up_item(room_object_t const &obj) const {return check_weight_limit(get_obj_weight(obj));}
	float get_carry_weight_ratio() const {return min(1.0f, cur_weight/global_building_params.player_weight_limit);}
	float get_speed_mult () const {return (1.0f - 0.4f*get_carry_weight_ratio())*((bladder > 0.9) ? 0.6 : 1.0);} // 40% reduction for heavy load, 40% reduction for full bladder
	float get_drunkenness() const {return drunkenness;}
	bool  player_is_dead () const {return (player_health <= 0.0);}
	bool  player_has_key () const {return has_key;}

	bool can_unlock_door() const {
		if (has_key) return 1;
		print_text_onscreen("Door is locked", RED, 1.0, 2.0*TICKS_PER_SECOND, 0);
		gen_sound_thread_safe_at_player(SOUND_CLICK, 1.0, 0.6);
		return 0;
	}
	void switch_item(bool dir) { // Note: current item is always carried.back()
		if (carried.size() <= 1) return; // no other item to switch to
		if (dir) {std::rotate(carried.begin(), carried.begin()+1, carried.end());}
		else     {std::rotate(carried.begin(), carried.end  ()-1, carried.end());}
		gen_sound_thread_safe_at_player(SOUND_CLICK, 0.5);
	}
	void add_item(room_object_t const &obj) {
		float health(0.0), drunk(0.0); // add these fields to bldg_obj_type_t?
		bool const bladder_was_full(bladder >= 0.9);
		float const value(get_obj_value(obj));
		damage_done += value;
		colorRGBA text_color(GREEN);
		std::ostringstream oss;
		oss << get_taken_obj_type(obj).name;

		if (is_consumable(obj)) { // nonempty bottle, consumable
			switch (obj.get_bottle_type()) {
			case 0: health =  0.25; break; // water
			case 1: health =  0.50; break; // Coke
			case 2: drunk  =  0.25; break; // beer
			case 3: drunk  =  0.50; break; // wine (entire bottle)
			case 4: health = -0.50; break; // poison - take damage
			default: assert(0);
			}
		}
		if (health > 0.0) { // heal
			player_health = min(1.0f, (player_health + health));
			oss << ": +" << round_fp(100.0*health) << "% Health";
		}
		if (health < 0.0) { // take damage
			player_health += health;
			oss << ": " << round_fp(100.0*health) << "% Health";
			text_color = RED;
			add_camera_filter(colorRGBA(RED, 0.25), 4, -1, CAM_FILT_DAMAGE); // 4 ticks of red damage
		}
		if (obj.type == TYPE_KEY) {
			has_key = 1; // mark as having the key, but it doesn't go into the inventory or contribute to weight or value
		}
		else if (health == 0.0 && drunk == 0.0) { // print value and weight if item is not consumed
			float const weight(get_obj_weight(obj));
			cur_value  += value;
			cur_weight += weight;
			
			if (obj.is_interactive()) {
				carried.push_back(obj);
				room_object_t &co(carried.back());

				if (obj.type == TYPE_BOOK) { // clear dim and dir for books
					float const dx(co.dx()), dy(co.dy()), dz(co.dz());

					if (dz > min(dx, dy)) { // upright book from a bookcase, put it on its side facing the player
						co.x2() = co.x1() + dz;
						co.y2() = co.y1() + max(dx, dy);
						co.z2() = co.z1() + min(dx, dy);
					}
					else if (co.dim) { // swap aspect ratio to make dim=0
						co.x2() = co.x1() + dy;
						co.y2() = co.y1() + dx;
					}
					co.dim = co.dir = 0;
					co.flags &= ~RO_FLAG_RAND_ROT; // remove the rotate bit
				}
				else if (obj.type == TYPE_PHONE) {
					if (co.dim) { // swap aspect ratio to make dim=0
						float const dx(co.dx()), dy(co.dy());
						co.x2() = co.x1() + dy;
						co.y2() = co.y1() + dx;
					}
					co.dim = co.dir = 0; // clear dim and dir
					phone_manager.enable();
				}
			}
			oss << ": value $";
			if (value < 1.0 && value > 0.0) {oss << ((value < 0.1) ? "0.0" : "0.") << round_fp(100.0*value);} // make sure to print the leading/trailing zero for cents
			else {oss << value;}
			oss << " weight " << get_obj_weight(obj) << " lbs";
		}
		else { // add one drink to the bladder, 25% of capacity
			bladder = min(1.0f, (bladder + 0.25f));
		}
		if (drunk > 0.0) {
			drunkenness += drunk;
			oss << ": +" << round_fp(100.0*drunk) << "% Drunkenness";
			if (drunkenness > 0.99f && (drunkenness - drunk) <= 0.99f) {oss << "\nYou are drunk"; text_color = DK_GREEN;}
		}
		if (!bladder_was_full && bladder >= 0.9f) {oss << "\nYou need to use the bathroom"; text_color = YELLOW;}
		print_text_onscreen(oss.str(), text_color, 1.0, 3*TICKS_PER_SECOND, 0);
	}
	bool take_person(bool &person_has_key, float person_height) {
		if (drunkenness < 1.5) { // not drunk enough
			print_text_onscreen("Not drunk enough", RED, 1.0, 2.0*TICKS_PER_SECOND, 0);
			return 0;
		}
		float const value(1000), weight((person_height > 0.025) ? 180.0 : 80.0); // always worth $1000; use height to select man vs. girl
		if (!check_weight_limit(weight)) {show_weight_limit_message(); return 0;}
		has_key    |= person_has_key; person_has_key = 0; // steal their key
		cur_value  += value;
		cur_weight += weight;
		std::ostringstream oss;
		oss << "zombie: value $" << value << " weight " << weight << " lbs";
		print_text_onscreen(oss.str(), GREEN, 1.0, 4*TICKS_PER_SECOND, 0);
		return 1; // success
	}
	bool try_use_last_item(room_object_t &obj) {
		if (carried.empty()) return 0; // no interactive carried item
		obj = carried.back(); // deep copy
		if (!obj.has_dstate()) {return obj.can_use();} // not a droppable/throwable item(ball); should always return 1
		remove_last_item(); // drop the item - remove it from our inventory
		return 1;
	}
	void mark_last_item_used() {
		assert(!carried.empty());
		carried_item_t &obj(carried.back());
		obj.flags |= RO_FLAG_USED;
		unsigned const capacity(bldg_obj_types[obj.type].capacity);

		if (capacity > 0) {
			++obj.use_count;
			if (obj.use_count >= capacity) {remove_last_item();} // remove after too many uses
		}
	}
	void remove_last_item() {
		assert(!carried.empty());
		room_object_t const &obj(carried.back());
		cur_value  -= get_obj_value (obj);
		cur_weight -= get_obj_weight(obj);
		cur_value   = 0.01*round_fp(100.0*cur_value ); // round to nearest cent
		cur_weight  = 0.01*round_fp(100.0*cur_weight); // round to nearest 0.01 lb
		assert(cur_value > -0.01 && cur_weight > -0.01); // sanity check for math
		max_eq(cur_value,  0.0f); // handle FP rounding error
		max_eq(cur_weight, 0.0f);
		carried.pop_back(); // Note: invalidates obj
	}
	void collect_items(bool keep_interactive) {
		if (!keep_interactive) {has_key = 0;} // key only good for current building
		phone_manager.disable(); // phones won't ring when taken out of their building, since the player can't switch to them anyway
		if (carried.empty() && cur_weight == 0.0 && cur_value == 0.0) return; // nothing to add
		float keep_value(0.0), keep_weight(0.0);

		if (keep_interactive) { // carried items don't contribute to collected value and weight; their value and weight remain in the inventory after collection
			for (auto i = carried.begin(); i != carried.end(); ++i) {
				keep_value  += get_obj_value(*i);
				keep_weight += get_obj_weight(*i);
			}
			cur_value  -= keep_value;
			cur_weight -= keep_weight;
		}
		else {carried.clear();}
		if (cur_weight == 0.0 && cur_value == 0.0) return; // no update to print
		std::ostringstream oss;
		oss << "Added value $" << cur_value << " Added weight " << cur_weight << " lbs\n";
		tot_value  += cur_value;  cur_value  = keep_value;
		tot_weight += cur_weight; cur_weight = keep_weight;
		oss << "Total value $" << tot_value << " Total weight " << tot_weight << " lbs";
		print_text_onscreen(oss.str(), GREEN, 1.0, 4*TICKS_PER_SECOND, 0);
	}
	void show_stats() const {
		if (!carried.empty()) {
			player_held_object = carried.back(); // deep copy last pickup object if usable
			
			if (player_held_object.type == TYPE_PHONE) {
				if (phone_manager.is_phone_ringing()) {player_held_object.flags |= RO_FLAG_EMISSIVE;} // show ring screen
				else if (phone_manager.is_phone_on()) {player_held_object.flags |= RO_FLAG_OPEN    ;} // show lock screen
			}
		}
		if (display_framerate) { // controlled by framerate toggle
			float const aspect_ratio((float)window_width/(float)window_height);

			if (cur_weight > 0.0 || tot_weight > 0.0 || best_value > 0.0) { // don't show stats until the player has picked something up
				std::ostringstream oss;
				oss << "Cur $" << cur_value << " / " << cur_weight << " lbs  Total $" << tot_value << " / " << tot_weight
					<< " lbs  Best $" << best_value << "  Damage $" << damage_done;
				
				if (!carried.empty()) {
					unsigned const capacity(bldg_obj_types[player_held_object.type].capacity);
					oss << "  [" << get_taken_obj_type(player_held_object).name; // print the name of the throwable object
					if (capacity > 0) {oss << " " << (capacity - carried.back().use_count) << "/" << capacity;} // print use/capacity
					oss << "]";
				}
				draw_text(GREEN, -0.010*aspect_ratio, -0.011, -0.02, oss.str(), 0.8); // size=0.8
			}
			if (in_building_gameplay_mode()) {
				float const lvl(min(cur_building_sound_level, 1.0f));
				unsigned const num_bars(round_fp(20.0*lvl));

				if (num_bars > 0) { // display sound meter
					colorRGBA const color(lvl, (1.0 - lvl), 0.0, 1.0); // green => yellow => orange => red
					draw_text(color, -0.005*aspect_ratio, -0.010, -0.02, std::string(num_bars, '#'));
				}
				if (player_is_hiding) {draw_text(LT_BLUE, -0.001*aspect_ratio, -0.009, -0.02, "[Hiding]");}
			}
		}
		if (in_building_gameplay_mode()) {
			// Note: shields is used for drunkenness; values are scaled from 0-1 to 0-100; powerup values are for bladder fullness
			draw_health_bar(100.0*player_health, 100.0*drunkenness, bladder, YELLOW, get_carry_weight_ratio(), WHITE);
			if (has_key) {show_key_icon();}
		}
	}
	void next_frame() {
		show_stats();
		phone_manager.next_frame(); // even if not in gameplay mode?
		if (!in_building_gameplay_mode()) return;
		// handle player fall damage logic
		point const camera_pos(get_camera_pos());
		float const fall_damage_start(3.0*CAMERA_RADIUS); // should be a function of building floor spacing?
		float const player_zval(camera_pos.z), delta_z(prev_player_zval - player_zval);
		
		if (camera_in_building != prev_in_building) {prev_in_building = camera_in_building;}
		else if (prev_player_zval != 0.0 && delta_z > fall_damage_start && camera_in_building) {
			// only take fall damage when inside the building (no falling off the roof for now)
			player_health -= 1.0f*(delta_z - fall_damage_start)/fall_damage_start;
			if (player_is_dead()) {register_player_death(SOUND_SQUISH, " of a fall"); return;} // dead
			gen_sound_thread_safe_at_player(SOUND_SQUISH, 0.5);
			add_camera_filter(colorRGBA(RED, 0.25), 4, -1, CAM_FILT_DAMAGE); // 4 ticks of red damage
			register_building_sound_at_player(1.0);
		}
		prev_player_zval = player_zval;
		// handle death events
		if (player_is_dead() ) {register_player_death(SOUND_SCREAM3, ""); return;} // dead
		if (drunkenness > 2.0) {register_player_death(SOUND_DROWN,   " of alcohol poisoning"); return;}
		// update state for next frame
		drunkenness = max(0.0f, (drunkenness - 0.0001f*fticks)); // slowly decrease over time
		
		if (player_near_toilet) { // empty bladder
			if (bladder > 0.9) {gen_sound_thread_safe_at_player(SOUND_GASP);} // urinate
			if (bladder > 0.0) { // toilet flush
#pragma omp critical(gen_sound)
				gen_delayed_sound(1.0, SOUND_FLUSH, camera_pos); // delay by 1s
				register_building_sound_at_player(0.5);
			}
			bladder = 0.0;
		}
		else if (bladder > 0.9) {
			bladder_time += fticks;

			if (bladder_time > 5.0*TICKS_PER_SECOND) { // play the "I have to go" sound
				gen_sound_thread_safe_at_player(SOUND_HURT);
				bladder_time = 0.0;
			}
		}
		player_near_toilet = 0;
	}
};

player_inventory_t player_inventory;

float get_player_drunkenness() {return player_inventory.get_drunkenness();}
float get_player_building_speed_mult() {return player_inventory.get_speed_mult();}
bool player_can_unlock_door() {return player_inventory.can_unlock_door();}

void register_building_sound_for_obj(room_object_t const &obj, point const &pos) {
	float const weight(get_obj_weight(obj)), volume((weight <= 1.0) ? 0.0 : min(1.0f, 0.01f*weight)); // heavier objects make more sound
	register_building_sound(pos, volume);
}

void register_broken_object(room_object_t const &obj) {player_inventory.record_damage_done(get_obj_value(obj));}

bool register_player_object_pickup(room_object_t const &obj, point const &at_pos) {
	bool const can_pick_up(player_inventory.can_pick_up_item(obj));

	if (!do_room_obj_pickup) { // player has not used the pickup key, but we can still use this to notify the player that an object can be picked up
		can_pickup_bldg_obj = (can_pick_up ? 1 : 2);
		return 0;
	}
	if (!can_pick_up) {
		show_weight_limit_message();
		return 0;
	}
	if (is_consumable(obj)) {gen_sound_thread_safe_at_player(SOUND_GULP, 1.0 );}
	else                    {gen_sound_thread_safe_at_player(SOUND_ITEM, 0.25);}
	register_building_sound_for_obj(obj, at_pos);
	do_room_obj_pickup = 0; // no more object pickups
	return 1;
}

bool building_t::player_pickup_object(point const &at_pos, vector3d const &in_dir) {
	if (!has_room_geom()) return 0;
	return interior->room_geom->player_pickup_object(*this, at_pos, in_dir);
}
bool building_room_geom_t::player_pickup_object(building_t &building, point const &at_pos, vector3d const &in_dir) {
	point at_pos_rot(at_pos);
	vector3d in_dir_rot(in_dir);
	building.maybe_inv_rotate_pos_dir(at_pos_rot, in_dir_rot);
	float const range(3.0*CAMERA_RADIUS), drawer_range_max(2.5*CAMERA_RADIUS);
	float drawer_range(drawer_range_max), obj_dist(0.0);
	int const obj_id(find_nearest_pickup_object(building, at_pos_rot, in_dir_rot, range, obj_dist));
	if (obj_id >= 0) {min_eq(drawer_range, obj_dist);} // only include drawers that are closer than the pickup object
	if (open_nearest_drawer(building, at_pos_rot, in_dir_rot, drawer_range_max, 1)) return 1; // try objects in drawers; pickup_item=1
	if (obj_id < 0) return 0; // no object to pick up
	room_object_t &obj(get_room_object_by_index(obj_id));

	if (obj.type == TYPE_SHELVES || (obj.type == TYPE_WINE_RACK && !(obj.flags & RO_FLAG_EXPANDED))) { // shelves or unexpanded wine rack
		assert(!(obj.flags & RO_FLAG_EXPANDED)); // should not have been expanded
		expand_object(obj);
		bool const picked_up(player_pickup_object(building, at_pos, in_dir)); // call recursively on contents
		// if we picked up an object, assume the VBOs have already been updated; otherwise we need to update them to expand this object
		if (!picked_up) {create_small_static_vbos(building);} // assumes expanded objects are all "small"
		return picked_up;
	}
	if (obj.type == TYPE_BCASE) {
		static vector<room_object_t> books;
		books.clear();
		get_bookcase_books(obj, books);
		int closest_obj_id(-1);
		float dmin_sq(0.0);
		point const p2(at_pos + in_dir*range);

		for (auto i = books.begin(); i != books.end(); ++i) {
			point p1c(at_pos), p2c(p2);
			if (!do_line_clip(p1c, p2c, i->d))  continue; // test ray intersection vs. bcube
			float const dsq(p2p_dist(at_pos, p1c)); // use closest intersection point
			if (dmin_sq > 0.0 && dsq > dmin_sq) continue; // not the closest
			closest_obj_id = (i - books.begin()); // valid pickup object
			dmin_sq = dsq; // this object is the closest
		} // for i
		if (dmin_sq == 0.0) return 0; // no book to pick up
		room_object_t &book(books[closest_obj_id]);
		if (!register_player_object_pickup(book, at_pos)) return 0;
		obj.set_combined_flags(obj.get_combined_flags() | (1<<(book.item_flags&31))); // set flag bit to remove this book from the bookcase
		player_inventory.add_item(book);
		update_draw_state_for_room_object(book, building, 1);
		return 1;
	}
	if (!register_player_object_pickup(obj, at_pos)) return 0;
	remove_object(obj_id, building);
	return 1;
}

void building_t::register_player_enter_building() const {
	// nothing to do yet
}
void building_t::register_player_exit_building() const {
	// only collect items in gameplay mode where there's a risk the player can lose them; otherwise, let the player carry items between buildings
	player_inventory.collect_items(!in_building_gameplay_mode());
}

bool has_cube_line_coll(point const &p1, point const &p2, vect_cube_t const &cubes) {
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {
		if (i->line_intersects(p1, p2)) return 1;
	}
	return 0;
}
bool building_t::check_for_wall_ceil_floor_int(point const &p1, point const &p2) const { // and interior doors
	if (!interior) return 0;
	for (unsigned d = 0; d < 2; ++d) {if (has_cube_line_coll(p1, p2, interior->walls[d])) return 1;}
	if (has_cube_line_coll(p1, p2, interior->ceilings) || has_cube_line_coll(p1, p2, interior->floors)) return 1; // or is only checking one good enough?
	return check_line_intersect_doors(p1, p2);
}

bool object_has_something_on_it(room_object_t const &obj, vector<room_object_t> const &objs, vector<room_object_t>::const_iterator objs_end) {
	// only these types can have objects placed on them (what about TYPE_SHELF?)
	if (obj.type != TYPE_TABLE && obj.type != TYPE_DESK && obj.type != TYPE_COUNTER && obj.type != TYPE_DRESSER && obj.type != TYPE_NIGHTSTAND &&
		obj.type != TYPE_BOX && obj.type != TYPE_CRATE && obj.type != TYPE_WINE_RACK && obj.type != TYPE_BOOK) return 0;

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type == TYPE_BLOCKER) continue; // ignore blockers (from removed objects)
		if (*i == obj)               continue; // skip self (bcube check)
		if (obj.type == TYPE_WINE_RACK && obj.contains_pt(i->get_cube_center())) return 1; // check for wine bottles left in wine rack
		if (i->z1() == obj.z2() && i->intersects_xy(obj))                        return 1; // zval has to match exactly
		if (obj.type == TYPE_BOX && obj.is_open() && obj.contains_cube(*i))      return 1; // open box with an object inside
	}
	return 0;
}

void building_t::maybe_inv_rotate_pos_dir(point &pos, vector3d &dir) const {
	if (is_rotated()) {
		do_xy_rotate_inv(bcube.get_cube_center(), pos);
		do_xy_rotate_normal_inv(dir);
	}
}

cube_t get_true_obj_bcube(room_object_t const &obj) {
	if (obj.type == TYPE_PEN || obj.type == TYPE_PENCIL) {
		cube_t obj_bcube(obj);
		obj_bcube.expand_in_dim(!obj.dim, obj.get_sz_dim(!obj.dim));
		return obj_bcube;
	}
	if (obj.type == TYPE_BED) { // do more accurate check with various parts of the bed
		cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
		get_bed_cubes(obj, cubes);
		return cubes[3]; // check mattress only, since we can only take the mattress, sheets, and pillows
	}
	if (obj.is_obj_model_type()) {
		cube_t obj_bcube(obj);
		obj_bcube.expand_by(-0.1*obj.get_size()); // since models don't fill their bcubes, shrink them a bit when doing a ray query
		return obj_bcube;
	}
	return obj; // unmodified
}

bool obj_has_open_drawers(room_object_t const &obj) {return ((obj.type == TYPE_NIGHTSTAND || obj.type == TYPE_DRESSER || obj.type == TYPE_DESK) && obj.drawer_flags);}

int building_room_geom_t::find_nearest_pickup_object(building_t const &building, point const &at_pos, vector3d const &in_dir, float range, float &obj_dist) const {
	int closest_obj_id(-1);
	float dmin_sq(0.0);
	point const p2(at_pos + in_dir*range);

	for (unsigned vect_id = 0; vect_id < 2; ++vect_id) {
		auto const &obj_vect((vect_id == 1) ? expanded_objs : objs), &other_obj_vect((vect_id == 1) ? objs : expanded_objs);
		unsigned const obj_id_offset((vect_id == 1) ? objs.size() : 0);
		auto objs_end((vect_id == 1) ? expanded_objs.end() : get_stairs_start()); // skip stairs and elevators
		auto other_objs_end((vect_id == 1) ? get_stairs_start() : expanded_objs.end());

		for (auto i = obj_vect.begin(); i != objs_end; ++i) {
			assert(i->type < NUM_ROBJ_TYPES);
			if (!bldg_obj_types[i->type].pickup) continue; // this object type can't be picked up
			cube_t const obj_bcube(get_true_obj_bcube(*i));
			point p1c(at_pos), p2c(p2);
			if (!do_line_clip(p1c, p2c, obj_bcube.d)) continue; // test ray intersection vs. bcube
			float const dsq(p2p_dist(at_pos, p1c)); // use closest intersection point
			if (dmin_sq > 0.0 && dsq > dmin_sq)       continue; // not the closest
		
			if (i->type == TYPE_CLOSET || (i->type == TYPE_STALL && i->shape != SHAPE_SHORT)) { // can only take short stalls (separating urinals)
				if (!i->is_open() && !i->contains_pt(at_pos)) { // stalls/closets block the player from taking toilets/boxes unless open, or the player is inside
					closest_obj_id = -1;
					dmin_sq = dsq;
				}
				continue;
			}
			if (i->type == TYPE_CHAIR) { // separate back vs. seat vs. legs check for improved accuracy
				cube_t cubes[3]; // seat, back, legs_bcube
				get_chair_cubes(*i, cubes);
				bool intersects(0);
				for (unsigned n = 0; n < 3 && !intersects; ++n) {intersects |= cubes[n].line_intersects(p1c, p2c);}
				if (!intersects) continue;
			}
			if (i->type == TYPE_MIRROR  && !i->is_house())                continue; // can only pick up mirrors from houses, not office buildings
			if (i->type == TYPE_TABLE   && i->shape == SHAPE_CUBE)        continue; // can only pick up short (TV) tables and cylindrical tables
			if (i->type == TYPE_BED     && (i->flags & RO_FLAG_TAKEN3))   continue; // can only take pillow, sheets, and mattress - not the frame
			if (i->type == TYPE_SHELVES && (i->flags & RO_FLAG_EXPANDED)) continue; // shelves are already expanded, can no longer select this object
			if (obj_has_open_drawers(*i))                                 continue; // can't take if any drawers are open
			if (object_has_something_on_it(*i,       obj_vect, objs_end)) continue; // can't remove a table, etc. that has something on it
			if (object_has_something_on_it(*i, other_obj_vect, other_objs_end)) continue; // check the other one as well
			if (building.check_for_wall_ceil_floor_int(at_pos, p1c))      continue; // skip if it's on the other side of a wall, ceiling, or floor
			closest_obj_id = (i - obj_vect.begin()) + obj_id_offset; // valid pickup object
			dmin_sq = dsq; // this object is the closest, even if it can't be picked up
		} // for i
	} // for vect_id
	obj_dist = sqrt(dmin_sq);
	return closest_obj_id;
}

bool is_counter(room_object_t const &obj) {return (obj.type == TYPE_COUNTER || obj.type == TYPE_KSINK);}

bool building_room_geom_t::open_nearest_drawer(building_t &building, point const &at_pos, vector3d const &in_dir, float range, bool pickup_item) {
	int closest_obj_id(-1);
	float dmin_sq(0.0);
	point const p2(at_pos + in_dir*range);
	auto objs_end(get_std_objs_end()); // skip buttons/stairs/elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		bool const is_counter_type(is_counter(*i) || i->type == TYPE_CABINET);
		if (!(i->type == TYPE_DRESSER || i->type == TYPE_NIGHTSTAND || i->type == TYPE_DESK || // drawers that can be opened or items picked up from
			(!pickup_item && is_counter_type))) continue; // doors that can be opened (no item pickup)
		cube_t bcube(*i);
		if (!is_counter_type) {bcube.d[i->dim][i->dir] += 0.75*(i->dir ? 1.0 : -1.0)*i->get_sz_dim(i->dim);} // expand outward to include open drawers
		point p1c(at_pos), p2c(p2);
		if (!do_line_clip(p1c, p2c, bcube.d)) continue; // test ray intersection vs. bcube
		float const dsq(p2p_dist_sq(at_pos, p1c)); // use closest intersection point
		if (dmin_sq > 0.0 && dsq > dmin_sq) continue; // not the closest
		if (building.check_for_wall_ceil_floor_int(at_pos, p1c)) continue; // skip if it's on the other side of a wall, ceiling, or floor
		closest_obj_id = (i - objs.begin());
		dmin_sq = dsq; // this object is the closest, even if it can't be picked up
	} // for i
	if (closest_obj_id < 0) return 0; // no object
	room_object_t &obj(objs[closest_obj_id]);
	room_object_t drawers_part;
	vect_cube_t drawers; // or doors
	bool const has_doors(is_counter(obj) || obj.type == TYPE_CABINET);

	// Note: this is a messy solution and must match the drawing code, but it's unclear how else we can get the location of the drawers
	if (has_doors) {get_cabinet_or_counter_doors(obj, drawers);}
	else {
		if (obj.type == TYPE_DESK) {
			if (!(obj.room_id & 3)) return 0; // no drawers for this desk
			drawers_part = get_desk_drawers_part(obj);
			bool const side(obj.obj_id & 1);
			drawers_part.d[!obj.dim][side] -= (side ? 1.0 : -1.0)*0.85*get_tc_leg_width(obj, 0.06);
		}
		else {
			drawers_part = get_dresser_middle(obj);
			drawers_part.expand_in_dim(!obj.dim, -0.5*get_tc_leg_width(obj, 0.10));
		}
		get_drawer_cubes(drawers_part, drawers, 0); // front_only=0
	}
	dmin_sq        = 0.0;
	closest_obj_id = -1;

	for (auto i = drawers.begin(); i != drawers.end(); ++i) {
		point p1c(at_pos), p2c(p2);
		if (!do_line_clip(p1c, p2c, i->d)) continue; // test ray intersection vs. drawer
		float const dsq(p2p_dist_sq(at_pos, p1c)); // use closest intersection point
		if (dmin_sq == 0.0 || dsq < dmin_sq) {closest_obj_id = (i - drawers.begin()); dmin_sq = dsq;} // update if closest
	}
	if (closest_obj_id < 0) return 0; // no drawer
	
	if (pickup_item && !has_doors) { // pick up item in drawer rather than opening drawer; no pickup items behind doors yet
		if (!(obj.drawer_flags & (1U << closest_obj_id))) return 0; // drawer is not open
		// Note: drawer cube passed in is not shrunk to the interior part, but that should be okay because we're not doing a line test against that object
		room_object_t const item(get_item_in_drawer(obj, drawers[closest_obj_id], closest_obj_id));
		if (item.type == TYPE_NONE) return 0; // no item
		if (!register_player_object_pickup(item, at_pos)) return 0;
		obj.item_flags |= (1U << closest_obj_id); // flag item as taken
		player_inventory.add_item(item);
		update_draw_state_for_room_object(item, building, 1);
	}
	else { // open or close the drawer/door
		obj.drawer_flags ^= (1U << (unsigned)closest_obj_id); // toggle flag bit for selected drawer
		cube_t const &drawer(drawers[closest_obj_id]);
		point const drawer_center(drawer.get_cube_center());
		if (has_doors) {building.play_door_open_close_sound(drawer_center, obj.is_open(), 0.5, 1.5);}
		else {gen_sound_thread_safe(SOUND_SLIDING, building.local_to_camera_space(drawer_center), 0.5);}
		register_building_sound(drawer_center, 0.4);
		
		if (has_doors) { // expand any items in the cabinet so that the player can pick them up
			// find any cabinets adjacent to this one in the other dim (inside corner) and ensure the opposing door is closed so that they don't intersect
			cube_t c_test(drawer);
			c_test.d[obj.dim][obj.dir] += (obj.dir ? 1.0 : -1.0)*drawer.get_sz_dim(!obj.dim); // expand outward by the width of the door

			for (auto i = objs.begin(); i != objs_end; ++i) {
				if ((is_counter(obj) && !is_counter(*i)) || (obj.type == TYPE_CABINET && i->type != TYPE_CABINET)) continue; // wrong object type
				if (i->dim == obj.dim) continue; // not opposing dim (also skips obj itself)
				float const dir_sign(i->dir ? 1.0 : -1.0);
				cube_t i_exp(*i);
				i_exp.d[i->dim][i->dir] += dir_sign*i->get_sz_dim(i->dim); // expand other counter/cabinet to account for open doors
				if (!i_exp.intersects(c_test)) continue;
				get_cabinet_or_counter_doors(*i, drawers);

				for (auto j = drawers.begin(); j != drawers.end(); ++j) {
					cube_t drawer_exp(*j);
					drawer_exp.d[i->dim][i->dir] += dir_sign*j->get_sz_dim(!i->dim); // expand outward by the width of the door
					if (drawer_exp.intersects(c_test)) {i->drawer_flags &= ~(1U << (j - drawers.begin()));} // make sure any intersecting doors are closed
				}
			} // for i
			// Note: expanding cabinets by opening a single door will allow the player to take items from anywhere in the cabinet, even if behind a closed door
			expand_object(obj);
			update_draw_state_for_room_object(obj, building, 0);
		}
		else {
			create_small_static_vbos(building); // only need to update small objects for drawers
		}
	}
	return 1;
}

room_object_t &building_room_geom_t::get_room_object_by_index(unsigned obj_id) {
	if (obj_id < objs.size()) {return objs[obj_id];}
	unsigned const exp_obj_id(obj_id - objs.size());
	assert(exp_obj_id < expanded_objs.size());
	return expanded_objs[exp_obj_id];
}

void building_room_geom_t::remove_object(unsigned obj_id, building_t &building) {
	room_object_t &obj(get_room_object_by_index(obj_id));
	room_object_t const old_obj(obj); // deep copy
	assert(obj.type != TYPE_ELEVATOR); // elevators require special updates for drawing logic and cannot be removed at this time
	player_inventory.add_item(obj);
	bldg_obj_type_t const type(get_taken_obj_type(obj)); // capture type before updating obj
	bool const is_light(obj.type == TYPE_LIGHT);

	if (obj.type == TYPE_PICTURE && !(obj.flags & RO_FLAG_TAKEN1)) {obj.flags |= RO_FLAG_TAKEN1;} // take picture, leave frame
	else if (obj.type == TYPE_TPROLL && !(obj.flags & (RO_FLAG_TAKEN1 | RO_FLAG_WAS_EXP))) {obj.flags |= RO_FLAG_TAKEN1;} // take toilet paper roll, leave holder; not for expanded TP rolls
	else if (obj.type == TYPE_BED) {
		if      (obj.flags & RO_FLAG_TAKEN2) {obj.flags |= RO_FLAG_TAKEN3;} // take mattress
		else if (obj.flags & RO_FLAG_TAKEN1) {obj.flags |= RO_FLAG_TAKEN2;} // take sheets
		else {obj.flags |= RO_FLAG_TAKEN1;} // take pillow(s)
	}
	else if (obj.type == TYPE_PLANT && !(obj.flags & RO_FLAG_ADJ_BOT)) { // plant not on a table/desk
		if      (obj.flags & RO_FLAG_TAKEN2) {obj.type = TYPE_BLOCKER;} // take pot - gone
		else if (obj.flags & RO_FLAG_TAKEN1) {obj.flags |= RO_FLAG_TAKEN2;} // take dirt
		else {obj.flags |= RO_FLAG_TAKEN1;} // take plant
	}
	else if (obj.type == TYPE_TOILET || obj.type == TYPE_SINK) { // leave a drain in the floor
		cube_t drain;
		drain.set_from_point(point(obj.xc(), obj.yc(), obj.z1()));
		drain.expand_by_xy(0.065*obj.dz());
		drain.z2() += 0.02*obj.dz();
		obj = room_object_t(drain, TYPE_DRAIN, obj.room_id, 0, 0, RO_FLAG_NOCOLL, obj.light_amt, SHAPE_CYLIN, DK_GRAY);
		create_small_static_vbos(building);
	}
	else { // replace it with an invisible blocker that won't collide with anything
		obj.type  = TYPE_BLOCKER;
		obj.flags = (RO_FLAG_NOCOLL | RO_FLAG_INVIS);
	}
	if (is_light) {clear_and_recreate_lights();}
	update_draw_state_for_room_object(old_obj, building, 1);
}

// Note: called when adding, removing, or moving objects
void building_room_geom_t::update_draw_state_for_room_object(room_object_t const &obj, building_t &building, bool was_taken) {
	// reuild necessary VBOs and other data structures
	if (obj.is_dynamic()) {mats_dynamic.clear();} // dynamic object
	else if (obj.type == TYPE_BUTTON && (obj.flags & RO_FLAG_IN_ELEV)) {update_dynamic_draw_data();} // interior elevator buttons are drawn as dynamic objects
	else { // static object
		assert(obj.type < NUM_ROBJ_TYPES);
		bldg_obj_type_t const type(was_taken ? get_taken_obj_type(obj) : bldg_obj_types[obj.type]);
		if (type.lg_sm & 2) {create_small_static_vbos(building);} // small object
		if (type.lg_sm & 1) {create_static_vbos      (building);} // large object
		if (type.is_model ) {create_obj_model_insts  (building);} // 3D model
		//if (type.ai_coll  ) {building.invalidate_nav_graph();} // removing this object should not affect the AI navigation graph
	}
	modified_by_player = 1; // flag so that we avoid re-generating room geom if the player leaves and comes back
}

int building_room_geom_t::find_avail_obj_slot() const {
	auto objs_end(get_stairs_start()); // skip stairs and elevators

	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->type == TYPE_BLOCKER) {return int(i - objs.begin());} // blockers are used as temporaries for room object placement and to replace removed objects
	}
	return -1; // no slot found
}
void building_room_geom_t::add_expanded_object(room_object_t const &obj) {
	for (auto i = expanded_objs.begin(); i != expanded_objs.end(); ++i) {
		if (i->type == TYPE_BLOCKER) {*i = obj; return;} // found a slot - done
	}
	expanded_objs.push_back(obj); // not found - in this case we can add a new object
}
bool building_room_geom_t::add_room_object(room_object_t const &obj, building_t &building, bool set_obj_id, vector3d const &velocity) {
	assert(obj.type != TYPE_LIGHT && get_room_obj_type(obj).pickup); // currently must be a pickup object, and not a light

	if (!set_obj_id && (obj.flags & RO_FLAG_WAS_EXP)) { // if object was expanded, and it's not a dynamic object, use an expanded slot (books, etc.)
		assert(velocity == zero_vector);
		add_expanded_object(obj);
	}
	else {
		int const obj_id(find_avail_obj_slot());
		if (obj_id < 0) return 0; // no slot found
		room_object_t &added_obj(get_room_object_by_index(obj_id));
		added_obj = obj; // overwrite with new object
		if (set_obj_id) {added_obj.obj_id = (uint16_t)(obj.has_dstate() ? allocate_dynamic_state() : obj_id);}
		if (velocity != zero_vector) {get_dstate(added_obj).velocity = velocity;}
	}
	update_draw_state_for_room_object(obj, building, 0);
	return 1;
}

bool is_movable(room_object_t const &obj) {
	if (obj.no_coll() || obj.type == TYPE_BLOCKER) return 0; // no blockers
	assert(obj.type < NUM_ROBJ_TYPES);
	bldg_obj_type_t const &bot(bldg_obj_types[obj.type]);
	return (bot.weight >= 40.0 && !bot.attached); // heavy non-attached objects, including tables
}
bool building_t::move_nearest_object(point const &at_pos, vector3d const &in_dir, float range) {
	assert(has_room_geom());
	int closest_obj_id(-1);
	float dmin_sq(0.0);
	point const p2(at_pos + in_dir*range);
	vector<room_object_t> &objs(interior->room_geom->objs), &expanded_objs(interior->room_geom->expanded_objs);
	auto objs_end(interior->room_geom->get_stairs_start()); // skip stairs and elevators

	// determine which object the player may be choosing to move
	for (auto i = objs.begin(); i != objs_end; ++i) {
		if (i->no_coll() || i->type == TYPE_BLOCKER) continue; // not interactive
		//if (!is_movable(*i)) continue; // not movable; should this check be done later, after setting closest_obj_id?
		cube_t const obj_bcube(get_true_obj_bcube(*i));
		point p1c(at_pos), p2c(p2);
		if (!do_line_clip(p1c, p2c, obj_bcube.d)) continue; // test ray intersection vs. bcube
		float const dsq(p2p_dist(at_pos, p1c)); // use closest intersection point
		if (dmin_sq > 0.0 && dsq > dmin_sq)       continue; // not the closest
		if (obj_has_open_drawers(*i))             continue; // can't move if any drawers are open
		if (object_has_something_on_it(*i, objs,          objs_end           )) continue; // can't move a table, etc. that has something on it
		if (object_has_something_on_it(*i, expanded_objs, expanded_objs.end())) continue; // check the other one as well
		if (check_for_wall_ceil_floor_int(at_pos, p1c)) continue; // skip if it's on the other side of a wall, ceiling, or floor
		closest_obj_id = (i - objs.begin()); // valid pickup object
		dmin_sq = dsq; // this object is the closest
	} // for i
	cout << TXT(closest_obj_id) << TXT(dmin_sq) << endl;
	if (closest_obj_id < 0) return 0;

	// determine move direction and distance
	room_object_t &obj(objs[closest_obj_id]);
	if (!is_movable(obj))        {cout << "not movable"   << endl; return 0;} // closest object isn't movable
	if (obj.contains_pt(at_pos)) {cout << "player inside" << endl; return 0;} // player is inside this object?
	float const move_dist(rand_uniform(0.5, 1.0)*CAMERA_RADIUS*(100.0f/max(75.0f, bldg_obj_types[obj.type].weight))); // heavier objects move less; add some global randomness
	vector3d delta(obj.closest_pt(at_pos) - at_pos);
	delta.z = 0.0; // XY only
	delta.normalize();

	// attempt to move the object
	for (unsigned mdir = 0; mdir < 3; ++mdir) { // X+Y, closer dim, further dim
		vector3d move_vector(zero_vector);
		if (mdir == 0) {move_vector = delta*move_dist;} // move diag in XY
		else { // move in one dim
			if (delta.x == 0.0 || delta.y == 0.0) break; // no more dims to try (only one mdir iteration)
			bool const dim(fabs(delta.x) < fabs(delta.y));
			move_vector[dim] = delta[dim]*move_dist;
		}
		for (unsigned n = 0; n < 5; ++n, move_vector *= 0.5) { // move in several incrementally smaller steps
			cout << TXT(mdir) << TXT(n) << TXT(move_vector.str()) << endl;
			room_object_t moved_obj(obj);
			moved_obj += move_vector; // only the position changes
			if (!is_obj_pos_valid(moved_obj, 1)) {cout << "bad pos" << endl; continue;} // try a smaller movement; keep_in_room=1
			bool bad_placement(0);

			for (auto i = objs.begin(); i != objs_end; ++i) { // do we need to check expanded_objs?
				if (i == objs.begin() + closest_obj_id)      continue; // skip self
				if (i->no_coll() || i->type == TYPE_BLOCKER) continue; // skip non-colliding objects and blockers that add clearance between objects as these won't block this object
				//if (i->intersects(obj)) continue; // assume that if the current object intersects something, it's okay to continue intersecting it -- doesn't work for chairs
				if (i->intersects(moved_obj)) {bad_placement = 1; break;}
			}
			if (bad_placement) {cout << "intersects obj" << endl; continue;} // intersects another object, try a smaller movement
			obj = moved_obj; // keep this placement
			interior->room_geom->update_draw_state_for_room_object(obj, *this, 0);
			gen_sound_thread_safe_at_player(SOUND_SLIDING);
			register_building_sound_at_player(0.7);
			cout << "success!" << endl;
			return 1; // success
		} // for n
	} // for mdir
	return 0; // failed
}

bool building_t::is_obj_pos_valid(room_object_t const &obj, bool keep_in_room) const {
	assert(interior);
	room_t const &room(get_room(obj.room_id));
	if (keep_in_room && !room.contains_cube(obj)) {cout << "outside room " << TXT(obj.str()) << TXT(room.str()) << endl; return 0;} // outside the room
	bool contained_in_part(0);

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		if (p->contains_cube(obj)) {contained_in_part = 1; break;}
	}
	if (!contained_in_part) {cout << "outside part" << endl; return 0;}

	for (unsigned d = 0; d < 2; ++d) { // check for wall intersection
		if (has_bcube_int_no_adj(obj, interior->walls[d])) {cout << "wall intersect" << endl; return 0;}
	}
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) { // check for door intersection
		if (is_cube_close_to_door(obj, 0.0, (i->open && door_opens_inward(*i, room)), *i, (i->dim ^ i->open_dir ^ i->hinge_side ^ 1))) {cout << "door intersect" << endl; return 0;}
	}
	if (has_bcube_int(obj, interior->stairwells) || has_bcube_int(obj, interior->elevators)) {cout << "stair/elevator intersect" << endl; return 0;}
	return 1;
}

void play_obj_fall_sound(room_object_t const &obj, point const &player_pos) {
	gen_sound_thread_safe(SOUND_OBJ_FALL, (get_camera_pos() + (obj.get_cube_center() - player_pos)));
	register_building_sound_for_obj(obj, player_pos);
}

bool building_t::maybe_use_last_pickup_room_object(point const &player_pos) {
	if (player_in_elevator) return 0; // can't use items in elevators
	assert(has_room_geom());
	static bool delay_use(0);
	static double last_use_time(0.0);
	if (delay_use && (tfticks - last_use_time) < 0.5*TICKS_PER_SECOND) return 0; // half second delay on prev item use or switch
	delay_use = 0;
	room_object_t obj;
	if (!player_inventory.try_use_last_item(obj)) return 0;

	if (obj.has_dstate()) { // it's a dynamic object (ball), throw it; only activated with use_object/'E' key
		float const cradius(get_scaled_player_radius());
		point dest(player_pos + (1.2f*(cradius + obj.get_radius()))*cview_dir);
		dest.z -= 0.5*cradius; // slightly below the player's face
		obj.translate(dest - point(obj.xc(), obj.yc(), obj.z1()));
		obj.flags |= RO_FLAG_DYNAMIC; // make it dynamic, assuming it will be dropped/thrown
		if (!interior->room_geom->add_room_object(obj, *this, 1, THROW_VELOCITY*cview_dir)) return 0;
		player_inventory.return_object_to_building(obj); // re-add this object's value
		play_obj_fall_sound(obj, player_pos);
		delay_use = 1;
	}
	else if (obj.can_use()) { // active with either use_object or fire key
		if (obj.type == TYPE_TPROLL) {
			point const dest(player_pos + (1.5f*get_scaled_player_radius())*cview_dir);
			if (!apply_toilet_paper(dest, cview_dir, 0.5*obj.dz())) return 0;
			player_inventory.mark_last_item_used();
		}
		else if (obj.type == TYPE_SPRAYCAN || obj.type == TYPE_MARKER) { // spraypaint or marker
			if (!apply_paint(player_pos, cview_dir, obj.color, obj.type)) return 0;
			player_inventory.mark_last_item_used();
		}
		else if (obj.type == TYPE_BOOK) {
			float const half_width(0.5*max(max(obj.dx(), obj.dy()), obj.dz()));
			point dest(player_pos + (1.2f*(get_scaled_player_radius() + half_width))*cview_dir);
			if (!get_zval_for_obj_placement(dest, half_width, dest.z, 0)) return 0; // no suitable placement found; add_z_bias=0
			// orient based on the player's primary direction
			bool const place_dim(fabs(cview_dir.y) < fabs(cview_dir.x));

			if (obj.dim != place_dim) {
				float const dx(obj.dx()), dy(obj.dy());
				obj.x2() = obj.x1() + dy;
				obj.y2() = obj.y1() + dx;
			}
			obj.dim    = place_dim;
			obj.dir    = ((cview_dir[!place_dim] > 0) ^ place_dim);
			obj.flags |= (RO_FLAG_TAKEN1 | RO_FLAG_WAS_EXP);
			obj.translate(dest - point(obj.xc(), obj.yc(), obj.z1()));
			if (!interior->room_geom->add_room_object(obj, *this)) return 0;
			player_inventory.return_object_to_building(obj); // re-add this object's value
			player_inventory.remove_last_item(); // used
			play_obj_fall_sound(obj, player_pos);
			delay_use = 1;
		}
		else if (obj.type == TYPE_PHONE) {phone_manager.player_action();}
		else {assert(0);}
	}
	else {assert(0);}
	last_use_time = tfticks;
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

// spraypaint, markers, and decals

bool line_int_cube_get_t(point const &p1, point const &p2, cube_t const &cube, float &tmin) {
	float tmin0(0.0), tmax0(1.0);
	if (get_line_clip(p1, p2, cube.d, tmin0, tmax0) && tmin0 < tmin) {tmin = tmin0; return 1;}
	return 0;
}
bool line_int_cubes_get_t(point const &p1, point const &p2, vect_cube_t const &cubes, float &tmin, cube_t &target) {
	bool had_int(0);

	for (auto c = cubes.begin(); c != cubes.end(); ++c) {
		if (line_int_cube_get_t(p1, p2, *c, tmin)) {target = *c; had_int = 1;}
	}
	return had_int;
}
vector3d get_normal_for_ray_cube_int_xy(point const &p, cube_t const &c, float tolerance) {
	vector3d n(zero_vector);

	for (unsigned d = 0; d < 2; ++d) { // find the closest intersecting cube XY edge, which will determine the normal vector
		if (fabs(p[d] - c.d[d][0]) < tolerance) {n[d] = -1.0; break;} // test low  edge
		if (fabs(p[d] - c.d[d][1]) < tolerance) {n[d] =  1.0; break;} // test high edge
	}
	return n;
}

// decals

class paint_manager_t : public paint_draw_t { // for paint on exterior walls/windows, viewed from inside the building
	building_t const *paint_bldg = nullptr;
public:
	bool have_paint_for_building() const { // only true if the building contains the player
		return (paint_bldg && !(qbd[0].empty() && qbd[1].empty()) && paint_bldg->bcube.contains_pt(get_camera_building_space()));
	}
	quad_batch_draw &get_paint_qbd(building_t const *const building, bool is_marker) {
		if (building != paint_bldg) { // paint switches to this building
			for (unsigned d = 0; d < 2; ++d) {qbd[d].clear();}
			paint_bldg = building;
		}
		return qbd[is_marker];
	}
};

paint_manager_t ext_paint_manager;
bool have_buildings_ext_paint() {return ext_paint_manager.have_paint_for_building();}
void draw_buildings_ext_paint() {ext_paint_manager.draw_paint();}

bool building_t::apply_paint(point const &pos, vector3d const &dir, colorRGBA const &color, room_object const obj_type) const { // spraypaint or marker
	bool const is_spraypaint(obj_type == TYPE_SPRAYCAN), is_marker(obj_type == TYPE_MARKER);
	assert(is_spraypaint || is_marker); // only these two are supported
	// find intersection point and normal; assumes pos is inside the building
	assert(has_room_geom());
	float const max_dist((is_spraypaint ? 16.0 : 3.0)*CAMERA_RADIUS), tolerance(0.01*get_wall_thickness());
	point const pos2(pos + max_dist*dir);
	float tmin(1.0);
	vector3d normal;
	cube_t target;
	
	for (unsigned d = 0; d < 2; ++d) {
		if (line_int_cubes_get_t(pos, pos2, interior->walls[d], tmin, target)) {
			normal    = zero_vector;
			normal[d] = -SIGN(dir[d]); // normal is opposite of ray dir in this dim
		}
	}
	if (line_int_cubes_get_t(pos, pos2, interior->floors  , tmin, target)) {normal =  plus_z;}
	if (line_int_cubes_get_t(pos, pos2, interior->ceilings, tmin, target)) {normal = -plus_z;}
	
	// include exterior walls; okay to add spraypaint and markers over windows
	cube_t const part(get_part_containing_pt(pos));
	float tmin0(0.0), tmax0(1.0);
	bool exterior_wall(0);

	if (get_line_clip(pos, pos2, part.d, tmin0, tmax0) && tmax0 < tmin) { // part edge is the closest intersection point
		// check other parts to see if ray continues into them; if not, it exited the building; this implementation isn't perfect but should be close enough
		point const cand_p_int(pos + tmax0*(pos2 - pos));
		bool found(0);

		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (*p == part || !p->contains_pt_exp(cand_p_int, tolerance)) continue; // ray does not continue into this new part
			if (check_line_clip(cand_p_int, pos2, p->d)) {found = 1; break;} // ray continues into this part
		}
		if (!found) { // ray has exited the building
			vector3d const n(-get_normal_for_ray_cube_int_xy(cand_p_int, part, tolerance)); // negate the normal because we're looking for the exit point from the cube
			if (n != zero_vector) {tmin = tmax0; normal = n; target = part; exterior_wall = 1;}
		}
	}
	// check for rugs, pictures, and whiteboards, which can all be painted over; also check for walls from closets
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
	bool const is_wall(normal.x != 0.0 || normal.y != 0.0), is_floor(normal == plus_z);
	bool walls_blocked(0);

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
 		if ((is_wall && (i->type == TYPE_PICTURE || i->type == TYPE_WBOARD || i->type == TYPE_MIRROR)) || (is_floor && (i->type == TYPE_RUG || i->type == TYPE_FLOORING))) {
			if (line_int_cube_get_t(pos, pos2, *i, tmin)) {target = *i;} // Note: return value is ignored, we only need to update tmin and target; normal should be unchanged
		}
		else if (i->type == TYPE_CLOSET && line_int_cube_get_t(pos, pos2, *i, tmin)) {
			point const cand_p_int(pos + tmin*(pos2 - pos));
			normal = get_normal_for_ray_cube_int_xy(cand_p_int, *i, tolerance); // should always return a valid normal
			target = *i;
		}
		else if (i->type == TYPE_STALL || i->type == TYPE_CUBICLE) {
			cube_t c(*i);

			if (i->type == TYPE_STALL && i->shape != SHAPE_SHORT) { // toilet stall, clip cube to wall height
				float const dz(c.dz());
				c.z2() -= 0.35*dz; c.z1() += 0.15*dz;
			}
			float tmin0(tmin);
			if (!line_int_cube_get_t(pos, pos2, c, tmin0)) continue;
			if (i->contains_pt(pos)) continue; // inside stall/cubicle, can't paint the exterior
			vector3d const n(get_normal_for_ray_cube_int_xy((pos + tmin0*(pos2 - pos)), c, tolerance)); // should always return a valid normal
			if (n[i->dim] != 0) {walls_blocked = 1; continue;} // only the side walls count
			tmin = tmin0; normal = n; target = c;
		}
	} // for i
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		float tmin0(tmin);
		if (!line_int_cube_get_t(pos, pos2, *i, tmin0)) continue;
		if (i->contains_pt(pos)) {walls_blocked = 1; continue;} // can't spraypaint the outside of the elevator when standing inside it
		vector3d const n(get_normal_for_ray_cube_int_xy((pos + tmin0*(pos2 - pos)), *i, tolerance)); // should always return a valid normal
		if (n[i->dim] == (i->dir ? 1.0 : -1.0)) {walls_blocked = 1; continue;} // skip elevator opening, even if not currently open
		tmin = tmin0; normal = n; target = *i;
	}
	for (auto i = interior->stairwells.begin(); i != interior->stairwells.end(); ++i) {
		if (i->shape == SHAPE_STRAIGHT) continue; // no walls, skip
		// expand by wall half-width; see building_t::add_stairs_and_elevators()
		float const step_len_pos(i->get_sz_dim(i->dim)/i->get_num_stairs());
		cube_t c(*i);
		c.expand_by_xy(0.15*step_len_pos); // wall half width
		float tmin0(tmin);
		if (!line_int_cube_get_t(pos, pos2, c, tmin0)) continue;
		if (c.contains_pt(pos)) {walls_blocked = 1; continue;} // can't spraypaint the outside of the stairs when standing inside them
		vector3d const n(get_normal_for_ray_cube_int_xy((pos + tmin0*(pos2 - pos)), c, tolerance)); // should always return a valid normal

		if (i->shape == SHAPE_U) {
			if (n[i->dim] == (i->dir ? -1.0 : 1.0)) {walls_blocked = 1; continue;} // skip stairs opening
		}
		else if (i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES) {
			// Note: we skip the end for SHAPE_WALLED and only check the sides because it depends on the floor we're on
			if (n[i->dim] != 0) {walls_blocked = 1; continue;} // skip stairs opening, either side
		}
		else {assert(0);} // unsupported stairs type
		tmin = tmin0; normal = n; target = c;
	}
	if (normal == zero_vector) return 0; // no walls, ceilings, floors, etc. hit
	if (walls_blocked && normal.z == 0.0) return 0; // can't spraypaint walls through elevator, stairs, etc.
	point p_int(pos + tmin*(pos2 - pos));
	if (check_line_intersect_doors(pos, p_int)) return 0; // blocked by door, no spraypaint; can't add spraypaint over door in case door is opened
	float const max_radius((is_spraypaint ? 2.0 : 0.035)*CAMERA_RADIUS);
	float const dist(p2p_dist(pos, p_int)), radius(is_spraypaint ? min(max_radius, max(0.05f*max_radius, 0.1f*dist)) : max_radius); // modified version of get_spray_radius()
	float const alpha((is_spraypaint && radius > 0.5*max_radius) ? (1.0 - (radius - 0.5*max_radius)/max_radius) : 1.0); // 0.5 - 1.0
	p_int += 0.01*radius*normal; // move slightly away from the surface
	assert(bcube.contains_pt(p_int));
	unsigned const dim(get_max_dim(normal)), d1((dim+1)%3), d2((dim+2)%3);

	// check that entire circle is contained in the target
	for (unsigned e = 0; e < 2; ++e) {
		unsigned const d(e ? d2 : d1);
		if (p_int[d] - 0.9*radius < target.d[d][0] || p_int[d] + 0.9*radius > target.d[d][1]) return 0; // extends outside the target surface in this dim
	}
	static point last_p_int(all_zeros);
	if (dist_less_than(p_int, last_p_int, 0.25*radius)) return 1; // too close to previous point, skip (to avoid overlapping sprays at the same location); still return 1
	last_p_int = p_int;
	vector3d dir1, dir2; // unit vectors
	dir1[d1] = 1.0; dir2[d2] = 1.0;
	float const winding_order_sign(-SIGN(normal[dim])); // make sure to invert the winding order to match the normal sign
	// Note: interior spraypaint draw uses back face culling while exterior draw does not; invert the winding order for exterior quads so that they show through windows correctly
	vector3d const dx(radius*dir1*winding_order_sign*(exterior_wall ? -1.0 : 1.0));
	interior->room_geom->decal_manager.paint_draw[exterior_wall].qbd[is_marker].add_quad_dirs(p_int, dx, radius*dir2, colorRGBA(color, alpha), normal); // add interior/exterior paint
	if (exterior_wall) {ext_paint_manager.get_paint_qbd(this, is_marker).add_quad_dirs(p_int, dx, radius*dir2, colorRGBA(color, alpha), normal);} // add exterior paint only
	static double next_sound_time(0.0);

	if (tfticks > next_sound_time) { // play sound if sprayed/marked, but not too frequently; marker has no sound
		gen_sound_thread_safe_at_player((is_spraypaint ? (int)SOUND_SPRAY : (int)SOUND_SQUEAK), 0.25);
		if (is_spraypaint) {register_building_sound(pos, 0.1);}
		next_sound_time = tfticks + double(is_spraypaint ? 0.5 : 0.25)*TICKS_PER_SECOND;
	}
	player_inventory.record_damage_done(is_spraypaint ? 1.0 : 0.1); // spraypaint does more damage than markers
	return 1;
}

bool room_object_t::can_place_onto() const {
	return (type == TYPE_TABLE || type == TYPE_DESK || type == TYPE_DRESSER || type == TYPE_NIGHTSTAND || type == TYPE_COUNTER || type == TYPE_KSINK ||
		type == TYPE_BRSINK || type == TYPE_BED || type == TYPE_BOX || type == TYPE_CRATE || type == TYPE_KEYBOARD || type == TYPE_BOOK);
}

bool building_t::get_zval_of_floor(point const &pos, float radius, float &zval) const {
	if (!interior) return 0; // error?
	cube_t cur_bcube;
	cur_bcube.set_from_sphere(pos, radius);
	if (!bcube.contains_cube_xy(cur_bcube)) return 0; // not contained/too close to walls
	float const floor_spacing(get_window_vspace());

	for (auto i = interior->floors.begin(); i != interior->floors.end(); ++i) { // blood can only be placed on floors
		if (pos.z < i->z2() || pos.z > (i->z2() + floor_spacing) || !i->contains_cube_xy(cur_bcube)) continue; // wrong floor, or not contained
		zval = (i->z2() + 0.0015*floor_spacing); // slightly above rugs (0.0015 vs. 0.001) and flooring (0.0015 vs. 0.0012)
		return 1;
	}
	return 0; // no suitable floor found
}
bool building_t::get_zval_for_obj_placement(point const &pos, float radius, float &zval, bool add_z_bias) const {
	float const start_zval(pos.z);
	if (!get_zval_of_floor(pos, radius, zval)) return 0; // if there's no floor, then there's probably no object to place on either
	if (!has_room_geom()) return 1; // probably can't get here
	float const z_bias(add_z_bias ? 0.0005*get_window_vspace() : 0.0); // maybe add a tiny bias to prevent z-fighting
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (!i->can_place_onto())    continue; // can't place on this object type
		if (!i->contains_pt_xy(pos)) continue; // center of mass not contained
		cube_t c(*i);

		if (i->type == TYPE_BED) {
			cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
			get_bed_cubes(*i, cubes);
			if (!cubes[3].contains_pt_xy(pos)) continue; // check again
			if (cubes[1].contains_pt_xy_exp(pos, radius) || cubes[2].contains_pt_xy_exp(pos, radius)) continue; // intersects the head or foot, skip
			c = cubes[3]; // mattress
		}
		if (c.z2() < zval || c.z2() > start_zval) continue; // below the floor or above the object's starting position

		if (i->shape == SHAPE_CYLIN) {
			if (!dist_xy_less_than(pos, i->get_cube_center(), i->get_radius())) continue; // round table
		}
		else {assert(i->shape != SHAPE_SPHERE);} // SHAPE_CUBE, SHAPE_TALL, SHAPE_SHORT are okay; others don't make sense
		zval = c.z2() + z_bias; // place on top of this object
	} // for i
	for (auto i = interior->room_geom->expanded_objs.begin(); i != interior->room_geom->expanded_objs.end(); ++i) { // check books, etc.
		if (!i->can_place_onto() || !i->contains_pt_xy(pos) || i->z2() < zval || i->z2() > start_zval) continue; // not a valid placement
		zval = i->z2() + z_bias; // place on top of this object
	}
	return 1;
}

bool building_t::apply_toilet_paper(point const &pos, vector3d const &dir, float half_width) {
	// for now, just drop a square of TP on the floor; could do better; should the TP roll shrink in size as this is done?
	assert(has_room_geom());
	static point last_tp_pos;
	if (dist_xy_less_than(pos, last_tp_pos, 1.5*half_width)) return 0; // too close to prev pos
	last_tp_pos = pos;
	float zval(pos.z);
	if (!get_zval_for_obj_placement(pos, half_width, zval, 1)) return 0; // no suitable placement found; add_z_bias=1
	vector3d d1(dir.x, dir.y, 0.0);
	if (d1 == zero_vector) {d1 = plus_x;} else {d1.normalize();}
	vector3d d2(cross_product(d1, plus_z));
	if (d2 == zero_vector) {d2 = plus_y;} else {d2.normalize();}
	interior->room_geom->decal_manager.tp_qbd.add_quad_dirs(point(pos.x, pos.y, zval), d1*half_width, d2*half_width, WHITE, plus_z);
	interior->room_geom->modified_by_player = 1; // make sure TP stays in this building
	// Note: no damage done for TP
	return 1;
}

void building_t::add_blood_decal(point const &pos) {
	assert(has_room_geom());
	float const radius(get_scaled_player_radius());
	float zval(pos.z);
	if (!get_zval_of_floor(pos, radius, zval)) return; // no suitable floor found
	tex_range_t const tex_range(tex_range_t::from_atlas((rand()&1), (rand()&1), 2, 2)); // 2x2 texture atlas
	interior->room_geom->decal_manager.blood_qbd.add_quad_dirs(point(pos.x, pos.y, zval), -plus_x*radius, plus_y*radius, WHITE, plus_z, tex_range); // Note: never cleared
	interior->room_geom->modified_by_player = 1; // make sure blood stays in this building
	player_inventory.record_damage_done(100.0); // blood is a mess to clean up (though damage will be reset on player death anyway)
}

// sound/audio tracking

void register_building_sound(point const &pos, float volume) {
	if (volume == 0.0 || !(show_bldg_pickup_crosshair || in_building_gameplay_mode())) return; // only when in gameplay/item pickup mode
	assert(volume > 0.0); // can't be negative
#pragma omp critical(building_sounds_update)
	{ // since this can be called by both the draw thread and the AI update thread, it should be in a critical section
		if (volume > ALERT_THRESH && cur_sounds.size() < 100) { // cap at 100 sounds in case they're not being cleared
			float const max_merge_dist(0.5*CAMERA_RADIUS);
			bool merged(0);

			for (auto i = cur_sounds.begin(); i != cur_sounds.end(); ++i) { // attempt to merge with an existing nearby sound
				if (dist_less_than(pos, i->pos, max_merge_dist)) {i->radius += volume; merged = 1;}
			}
			if (!merged) {cur_sounds.emplace_back(pos, volume);} // Note: volume is stored in radius field of sphere_t
		}
		cur_building_sound_level += volume;
	}
}
void register_building_sound_at_player(float volume) {
	register_building_sound(get_camera_building_space(), 1.0);
}

bool get_closest_building_sound(point const &at_pos, point &sound_pos, float floor_spacing) {
	if (cur_sounds.empty()) return 0;
	float max_vol(0.0); // 1.0 at a sound=1.0 volume at a distance of floor_spacing

	for (auto i = cur_sounds.begin(); i != cur_sounds.end(); ++i) {
		float vol(i->radius/max(0.01f*floor_spacing, p2p_dist(i->pos, at_pos)));
		if (fabs(i->pos.z - at_pos.z) > 0.75f*floor_spacing) {vol *= 0.5;} // half the volume when the sound comes from another floor
		if (vol > max_vol) {max_vol = vol; sound_pos = i->pos;}
	} // for i
	//cout << TXT(cur_sounds.size()) << TXT(max_vol) << endl;
	return (max_vol*floor_spacing > 0.06f);
}

void maybe_play_zombie_sound(point const &sound_pos_bs, unsigned zombie_ix, bool alert_other_zombies, bool high_priority) {
	unsigned const NUM_ZSOUNDS = 5;
	static rand_gen_t rgen;
	static double next_time_all(0.0), next_times[NUM_ZSOUNDS] = {};
	if (!high_priority && tfticks < next_time_all) return; // don't play any sound too frequently
	if (!high_priority && (rgen.rand()&3) != 0)    return; // only generate a sound 25% of the time (each frame), to allow more than one zombie to get a chance
	unsigned const sound_id(zombie_ix%NUM_ZSOUNDS); // choose one of the zombie sounds, determined by the current zombie
	double &next_time(next_times[sound_id]);
	if (!high_priority && tfticks < next_time) return; // don't play this particular sound too frequently
	next_time_all = tfticks + double(rgen.rand_uniform(1.0, 2.0))*TICKS_PER_SECOND; // next sound of any  type can play between 0.8 and 2.0s in the future
	next_time     = tfticks + double(rgen.rand_uniform(2.5, 5.0))*TICKS_PER_SECOND; // next sound of this type can play between 2.5 and 5.0s in the future
	gen_sound_thread_safe((SOUND_ZOMBIE1 + sound_id), (sound_pos_bs + get_camera_coord_space_xlate()));
	if (alert_other_zombies) {register_building_sound(sound_pos_bs, 0.4);}
}

void water_sound_manager_t::register_running_water(room_object_t const &obj, building_t const &building) {
	if (!(obj.flags & RO_FLAG_IS_ACTIVE)) return; // not turned on
	if (fabs(obj.z2() - camera_bs.z) > building.get_window_vspace()) return; // on the wrong floor
	point const pos(obj.get_cube_center());
	float const dsq(p2p_dist_sq(pos, camera_bs));
	if (dmin_sq == 0.0 || dsq < dmin_sq) {closest = pos; dmin_sq = dsq;}
}
void water_sound_manager_t::finalize() {
	if (dmin_sq == 0.0) return; // no water found
	static point prev_closest;
	bool const skip_if_already_playing(closest == prev_closest); // don't reset sound loop unless it moves to a different sink
	prev_closest = closest;
	gen_sound_thread_safe(SOUND_SINK, closest, 1.0, 1.0, 0.06, skip_if_already_playing); // fast distance falloff; will loop at the end if needed
}

// gameplay logic

bool player_has_room_key() {return player_inventory.player_has_key();}

// return value: 0=no effect, 1=player is killed, 2=this person is killed
int register_ai_player_coll(bool &has_key, float height) {
	if (do_room_obj_pickup && player_inventory.take_person(has_key, height)) {
		gen_sound_thread_safe_at_player(SOUND_ITEM, 0.5);
		do_room_obj_pickup = 0; // no more object pickups
		return 2;
	}
	static double last_coll_time(0.0);
	
	if (tfticks - last_coll_time > 2.0*TICKS_PER_SECOND) {
		gen_sound_thread_safe_at_player(SOUND_SCREAM1);
		last_coll_time = tfticks;
	}
	add_camera_filter(colorRGBA(RED, 0.25), 1, -1, CAM_FILT_DAMAGE); // 1 tick of red damage
	player_inventory.take_damage(0.04*fticks); // take damage over time
	
	if (player_inventory.player_is_dead()) {
		if (player_has_room_key()) {has_key = 1;}
		return 1;
	}
	return 0;
}

void building_gameplay_action_key(int mode, bool mouse_wheel) {
	if (camera_in_building) { // building interior action
		if (mouse_wheel && mode <= 1) {player_inventory.switch_item(mode);}
		// show crosshair on first pickup because it's too difficult to pick up objects without it
		else if (mode == 1) {do_room_obj_pickup = show_bldg_pickup_crosshair = 1;} // 'e'
		else if (mode == 2) {use_last_pickup_object = 1;} // 'E'
		else                {building_action_key    = 1;} // 'q'
	}
	else { // building exterior/city/road/car action
		if      (mode != 0) {city_action_key = 1;} // 'e'/'E'
		else                {} // 'q'
	}
}

void building_gameplay_next_frame() {
	if (office_chair_rot_rate != 0.0) { // update office chair rotation
		office_chair_rot_rate *= exp(-0.05*fticks); // exponential slowdown
		if (office_chair_rot_rate < 0.001) {office_chair_rot_rate = 0.0;} // stop rotating
	}
	if (in_building_gameplay_mode()) { // run gameplay update logic
		show_bldg_pickup_crosshair = 1;
		// update sounds used by AI
		auto i(cur_sounds.begin()), o(i);

		for (; i != cur_sounds.end(); ++i) {
			i->radius *= exp(-0.04*fticks);
			if (i->radius > ALERT_THRESH) {*(o++) = *i;} // keep if above thresh
		}
		cur_sounds.erase(o, cur_sounds.end());
	}
	player_held_object = carried_item_t();
	player_inventory.next_frame();
	// reset state for next frame
	cur_building_sound_level = min(1.2f, max(0.0f, (cur_building_sound_level - 0.01f*fticks))); // gradual decrease
	can_pickup_bldg_obj = 0;
	do_room_obj_pickup  = city_action_key = 0;
}

void enter_building_gameplay_mode() {player_inventory.clear_all();}

