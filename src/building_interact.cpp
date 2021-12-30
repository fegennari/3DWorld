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

extern bool tt_fire_button_down, flashlight_on, player_is_hiding, use_last_pickup_object, city_action_key;
extern int player_in_closet, camera_surf_collide, building_action_key, can_pickup_bldg_obj;
extern float fticks, CAMERA_RADIUS, office_chair_rot_rate;
extern double tfticks, camera_zh;
extern building_dest_t cur_player_building_loc;


cube_t get_sink_cube(room_object_t const &c);
bool player_can_open_door(door_t const &door);
bool player_has_room_key();
void register_broken_object(room_object_t const &obj);

// Note: pos is in camera space
void gen_sound_thread_safe(unsigned id, point const &pos, float gain, float pitch, float gain_scale, bool skip_if_already_playing) {
	float const dist(p2p_dist(get_camera_pos(), pos)), dscale(10.0*CAMERA_RADIUS*gain_scale); // distance at which volume is halved
	gain *= dscale/(dist + dscale);
	if (gain < 0.025) return; // too soft to hear
#pragma omp critical(gen_sound)
	gen_sound(id, pos, gain, pitch, 0, zero_vector, skip_if_already_playing);
}

// lights

float get_radius_for_room_light(room_object_t const &obj);

// Note: called by the player; closest_to is in building space, not camera space
bool building_t::toggle_room_light(point const &closest_to, bool sound_from_closest_to, int room_id, bool inc_lamps) {
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
		if (!i->is_light_type() || (!inc_lamps && i->type == TYPE_LAMP) || i->room_id != room_id) continue; // not a light, or the wrong room
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
	bool const ring_doorbell(is_open && is_house && door_ix == 0 && city_action_key); // action key when front door of a house is open
	if (is_open == was_open && !ring_doorbell) return; // no state change
	unsigned const dix(is_open ? (unsigned)door_ix : (unsigned)open_door_ix);
	assert(dix < doors.size());
	auto const &door(doors[dix]);
	point const door_center(door.get_bcube().get_cube_center()), sound_pos(local_to_camera_space(door_center)); // convert to camera space

	if (ring_doorbell) {
		gen_sound_thread_safe(SOUND_DOORBELL, sound_pos);
		return;
	}
	gen_sound_thread_safe((is_open ? (unsigned)SOUND_DOOR_OPEN : (unsigned)SOUND_DOOR_CLOSE), sound_pos);
	vector3d const normal(door.get_norm());
	bool const dim(fabs(normal.x) < fabs(normal.y)), dir(normal[dim] < 0.0);
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

// called for the player; mode: 0=normal, 1=pull
bool building_t::apply_player_action_key(point const &closest_to_in, vector3d const &in_dir_in, int mode, bool check_only) {
	if (!interior) return 0; // error?
	float const dmax(4.0*CAMERA_RADIUS), floor_spacing(get_window_vspace());
	float closest_dist_sq(0.0), t(0.0); // t is unused
	unsigned door_ix(0), obj_ix(0);
	bool found_item(0), is_obj(0);
	vector3d in_dir(in_dir_in);
	point closest_to(closest_to_in);
	maybe_inv_rotate_pos_dir(closest_to, in_dir);
	point const query_ray_end(closest_to + dmax*in_dir);

	if (!player_in_closet && mode == 0) { // if the player is in the closet, only the closet door can be opened
		for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
			if (i->z1() > closest_to.z || i->z2() < closest_to.z) continue; // wrong floor, skip
			point const center(i->get_cube_center());
			float const dist_sq(p2p_dist_sq(closest_to, center));
			if (found_item && dist_sq >= closest_dist_sq) continue; // not the closest
			if (!check_obj_dir_dist(closest_to, in_dir, *i, center, dmax)) continue; // door is not in the correct direction or too far away, skip
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
		vector<room_object_t> &objs(interior->room_geom->objs), &expanded_objs(interior->room_geom->expanded_objs);
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
			auto obj_vect_end((vect_id == 1) ? expanded_objs.end() : objs_end); // skip stairs and elevators

			for (auto i = obj_vect.begin(); i != obj_vect_end; ++i) {
				if (cur_player_building_loc.room_ix >= 0 && i->room_id != cur_player_building_loc.room_ix && i->type != TYPE_BUTTON) continue; // not in the same room as the player
				if (!active_area.is_all_zeros() && !i->intersects(active_area)) continue; // out of reach for the player
				bool keep(0);
				if (i->type == TYPE_BOX && !i->is_open()) {keep = 1;} // box can only be opened once; check first so that selection works for boxes in closets
				else if (i->type == TYPE_CLOSET) {
					if (in_dir.z > 0.5) continue; // not looking up at the light

					if (/*!i->is_open() &&*/ i->is_small_closet() && !interior->room_geom->moved_obj_ids.empty()) { // only applies to small closets
						cube_t c_test(*i);
						float const width(i->get_sz_dim(!i->dim)), wall_width(0.5*(width - 0.5*i->dz())); // see get_closet_cubes()
						c_test.d[i->dim][ i->dir] += (i->dir ? 1.0f : -1.0f)*(width - 2.0f*wall_width); // extend outward
						c_test.d[i->dim][!i->dir]  = i->d[i->dim][i->dir]; // back is flush with front of closet
						c_test.expand_in_dim(!i->dim, -wall_width); // shrink to door width
						if (interior->room_geom->cube_intersects_moved_obj(c_test)) continue; // blocked, can't open
					}
					keep = 1; // closet door can be opened
				}
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
		if (!player_in_closet && mode == 0) {
			float const drawer_dist(found_item ? sqrt(closest_dist_sq) : 2.5*CAMERA_RADIUS);
			
			if (interior->room_geom->open_nearest_drawer(*this, closest_to, in_dir, drawer_dist, 0, check_only)) { // drawer is closer - open or close it
				return (check_only ? 1 : 0); // check_only returns 1 here because this counts as interactive
			}
		}
		if (!found_item && !check_only && !player_in_closet) {move_nearest_object(closest_to, in_dir, 3.0*CAMERA_RADIUS, mode);} // try to move an object instead
		if (!found_item) return 0; // no door or object found
	}
	if (check_only) return 1;

	if (is_obj) { // interactive object
		if (!interact_with_object(obj_ix, closest_to, in_dir)) return 0; // generate sound from the player height
	}
	else { // interior door
		door_t &door(interior->doors[door_ix]);
		if (!player_can_open_door(door)) return 0; // locked/blocked
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
		toggle_room_light(obj.get_cube_center(), 1, obj.room_id, 0); // should select the correct light(s) for the room containing the switch; exclude lamps
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
	auto &objs(interior->room_geom->objs), &expanded_objs(interior->room_geom->expanded_objs);

	for (auto c = objs.begin(); c != objs.end(); ++c) { // check for other objects to collide with (including stairs)
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
	if (player_in_closet) { // check for collisions with expanded objects in closets
		bool updated_rot(0);

		for (auto c = expanded_objs.begin(); c != expanded_objs.end(); ++c) {
			if (c->type == TYPE_CLOTHES && sphere_cube_intersect(player_pos, player_radius, *c)) { // shirt in a closet with the player
				assert(c != objs.begin());
				room_object_t &hanger(*(c-1)); // hanger is the previous object
				assert(hanger.type == TYPE_HANGER);
				bool const rot_dir(dot_product(c->get_dir(), (player_pos - c->get_cube_center())) < 0);
				unsigned const orig_flags(c->flags), add_flag(rot_dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI), rem_flag(rot_dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
				c->flags     |= ((c->type == TYPE_CLOTHES) ? RO_FLAG_ROTATING : 0) | add_flag;
				c->flags     &= ~rem_flag;
				hanger.flags |= RO_FLAG_ROTATING | add_flag;
				hanger.flags &= ~rem_flag;
				updated_rot  |= (c->flags != orig_flags);
			}
		} // for c
		if (updated_rot) {interior->room_geom->create_small_static_vbos(*this);} // play a sound as well?
	}
	if (use_last_pickup_object || (tt_fire_button_down && !flashlight_on)) { // use object not active, and not using fire key without flashlight (space bar)
		maybe_use_last_pickup_room_object(player_pos);
		use_last_pickup_object = 0; // reset for next frame
	}
	maybe_update_tape(player_pos, 0); // end_of_tape=0
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
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
		if (i->open) continue; // check only closed doors
		if (i->get_true_bcube().line_intersects(p1, p2)) return 1;
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

// room objects

bool building_t::is_obj_pos_valid(room_object_t const &obj, bool keep_in_room) const {
	assert(interior);
	room_t const &room(get_room(obj.room_id));

	if (keep_in_room) {
		cube_t place_area(room);
		place_area.expand_by_xy(-0.99*get_trim_thickness()); // shrink to exclude wall trim
		if (!place_area.contains_cube(obj)) return 0; // outside the room
	}
	bool contained_in_part(0);

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		if (p->contains_cube(obj)) {contained_in_part = 1; break;}
	}
	if (!contained_in_part) return 0;

	for (unsigned d = 0; d < 2; ++d) { // check for wall intersection
		if (has_bcube_int_no_adj(obj, interior->walls[d])) return 0;
	}
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) { // check for door intersection
		if (!i->open || !door_opens_inward(*i, room)) continue; // closed, or opens into the adjacent room, ignore
		if (is_cube_close_to_door(obj, 0.0, 1, *i, i->get_check_dirs())) return 0;
	}
	if (has_bcube_int(obj, interior->stairwells) || has_bcube_int(obj, interior->elevators)) return 0;
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

