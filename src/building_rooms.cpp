// 3D World - Building Interior Room Assignment, etc.
// by Frank Gennari 4/30/2020

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t

extern object_model_loader_t building_obj_model_loader;

void setup_bldg_obj_types();
bool get_wall_quad_window_area(vect_vnctcc_t const &wall_quad_verts, unsigned i, cube_t &c, float &tx1, float &tx2, float &tz1, float &tz2);
void get_balcony_pillars(room_object_t const &c, float ground_floor_z1, cube_t pillar[2]);
void expand_convex_polygon_xy(vect_point &points, point const &center, float expand);


unsigned light_ix_assign_t::get_ix_for_light(cube_t const &c) {
	point2d<float> const pos(c.x1(), c.y1());

	for (auto i = cur.begin(); i != cur.end(); ++i) {
		if (i->first == pos) return i->second; // existing light is part of the same stack and is valid to return
	}
	cur.emplace_back(pos, get_next_ix()); // allocate a new light
	return cur.back().second;
}

colorRGBA get_light_color_temp(float t) {
	// 0.0: 1.0 1.0 0.5
	// 0.5: 1.0 1.0 1.0
	// 1.0: 0.5 0.5 1.0
	if (t > 0.5) {return colorRGBA(1.5-t, 1.5-t, 1.0  );} // high temp blue   spectrum
	else         {return colorRGBA(1.0,   1.0,   t+0.5);} // low  temp yellow spectrum
}
colorRGBA get_light_color_temp_range(float tmin, float tmax, rand_gen_t &rgen) {
	return get_light_color_temp(((tmin == tmax) ? tmin : rgen.rand_uniform(tmin, tmax)));
}

// Note: applies to both houses and office buildings, but only houses will have bedrooms
bool building_t::can_be_bedroom_or_bathroom(room_t const &room, unsigned floor, bool skip_conn_check) const { // check room type and existence of exterior door
	if (room.has_stairs_on_floor(floor) || room.has_elevator || room.is_hallway || room.is_office || room.is_sec_bldg) return 0; // no bed/bath in these cases
	
	if (floor == 0 && room.z1() == ground_floor_z1) {
		// run special logic for bedrooms and bathrooms (private rooms) on the first floor of a house
		if (is_room_adjacent_to_ext_door(room)) return 0; // door to house does not open into a bedroom/bathroom
		if (skip_conn_check) return 1;

		if (room.dz() > 1.5*get_window_vspace()) { // more than one floor
			if (interior->stairwells.empty()) return 1; // failed to place stairs in this house, maybe because it was too small; I guess we just return 1 here
			// determine if this room is on the shortest path from an exterior door to the stairs; if so, it can't be a bedroom or bathroom;
			// okay, that's not easy/fast to do, so determine if there is any path from the exterior door to the stairs that doesn't go through this room;
			// this won't work when there are two paths from the door to the stairs and this room is only on one of the paths, so we could put a BR/BR on both paths
			int cur_room(-1);
			vector<unsigned> door_rooms, stairs_rooms;

			for (unsigned i = 0; i < interior->rooms.size(); ++i) {
				room_t const &r(interior->rooms[i]);
				if (r == room) {cur_room = i; continue;} // this room; we know it can't have stairs or an exterior door
				if (r.is_sec_bldg || r.z2() <= ground_floor_z1) continue; // skip basement rooms, garages, and sheds
				if (r.has_stairs_on_floor(floor))  {stairs_rooms.push_back(i);}
				if (is_room_adjacent_to_ext_door(r)) {door_rooms.push_back(i);}
			}
			if (cur_room < 0 || stairs_rooms.empty()) {cout << TXT(bcube.str());}
			assert(cur_room >= 0); // must be found
			assert(!stairs_rooms.empty());

			if (!is_rotated() && is_cube() && !has_complex_floorplan) { // too strong for rotated or non-cube buildings, where door placement can sometimes fail
				assert(!doors.empty());
				assert(!door_rooms.empty());
			}
			for (auto d = door_rooms.begin(); d != door_rooms.end(); ++d) {
				for (auto s = stairs_rooms.begin(); s != stairs_rooms.end(); ++s) {
					if (!are_rooms_connected_without_using_room(*d, *s, cur_room)) return 0;
				}
			}
		}
	}
	return 1;
}
bool building_t::can_be_bathroom(room_t const &room) const { // Note: assumes caller has checked can_be_bedroom_or_bathroom()
	float const vspace(get_window_vspace());
	return (min(room.dx(), room.dy()) < 2.4*vspace && max(room.dx(), room.dy()) < 3.2*vspace && count_num_int_doors(room) == 1);
}
unsigned building_t::count_num_int_doors(room_t const &room) const {
	cube_t room_exp(room);
	float const wall_thickness(get_wall_thickness());
	room_exp.expand_by(wall_thickness, wall_thickness, -wall_thickness); // expand in XY and shrink in Z
	unsigned num(0);
	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {num += i->intersects(room_exp);}
	return num;
}

void expand_to_nonzero_area(cube_t &c, float exp_amt, bool dim) {
	while (c.get_sz_dim(dim) == 0.0) {
		c.expand_in_dim(dim, exp_amt);
		exp_amt *= 2.0;
	}
}
void set_light_xy(cube_t &light, point const &center, float light_size, bool light_dim, room_obj_shape light_shape) {
	for (unsigned dim = 0; dim < 2; ++dim) {
		float const sz(((light_shape == SHAPE_CYLIN || light_shape == SHAPE_SPHERE) ? 1.6 : ((bool(dim) == light_dim) ? 2.2 : 1.0))*light_size);
		light.d[dim][0] = center[dim] - sz;
		light.d[dim][1] = center[dim] + sz;
	}
}
unsigned calc_num_floors_room(room_t const &r, float window_vspacing, float floor_thickness) {
	return (r.is_sec_bldg ? 1 : calc_num_floors(r, window_vspacing, floor_thickness));
}

void building_t::gen_room_details(rand_gen_t &rgen, unsigned building_ix) {

	assert(interior);
	if (interior->room_geom) return; // already generated?
	setup_bldg_obj_types(); // initialize object types if not already done
	//highres_timer_t timer("Gen Room Details");
	// Note: people move from room to room, so using their current positions for room object generation is both nondeterministic and unnecessary
	vect_cube_t blockers, valid_lights; // blockers are used for fireplaces
	interior->room_geom.reset(new building_room_geom_t(bcube.get_llc()));
	vect_room_object_t &objs(interior->room_geom->objs);
	vector<room_t> &rooms(interior->rooms);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const light_thick(0.025*window_vspacing), def_light_size(0.1*window_vspacing);
	interior->room_geom->obj_scale = window_vspacing; // used to scale room object textures
	unsigned tot_num_rooms(0), num_bathrooms(0), num_bedrooms(0);
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {tot_num_rooms += calc_num_floors_room(*r, window_vspacing, floor_thickness);}
	objs.reserve(tot_num_rooms); // placeholder - there will be more than this many
	float const extra_bathroom_prob((is_house ? 2.0 : 1.0)*0.02*min((int(tot_num_rooms) - 4), 20));
	room_obj_shape const light_shape(is_house ? SHAPE_CYLIN : SHAPE_CUBE);
	unsigned cand_bathroom(rooms.size()); // start at an invalid value
	unsigned added_kitchen_mask(0); // per-floor
	unsigned added_bathroom_objs_mask(0);
	bool added_bedroom(0), added_living(0), added_library(0), added_dining(0), added_laundry(0), added_basement_utility(0), added_fireplace(0);
	light_ix_assign_t light_ix_assign;
	interior->create_fc_occluders(); // not really part of room geom, but needed for generating and drawing room geom, so we create them here
	has_int_fplace = 0; // reset for this generation

	if (rooms.size() > 1) { // choose best room assignments for required rooms; if a single room, skip this step
		float min_score(0.0);

		// Note: assigning cand_bathroom when has_pri_hall() is not strictly necessary, but may help add a bathroom to an upper stacked part
		for (auto r = rooms.begin(); r != rooms.end(); ++r) {
			if (r->is_sec_bldg) continue; // garage/shed excluded - not a normal room
			if (has_basement() && r->part_id == (int)basement_part_ix) continue; // skip the basement
			unsigned const num_floors(calc_num_floors_room(*r, window_vspacing, floor_thickness));

			// find best bathroom with no hard size constraints;
			// use the top floor for the test since it's less restrictive than the ground floor; will be checked per-floor later
			if (can_be_bedroom_or_bathroom(*r, (num_floors-1), 0)) { // skip_conn_check=0
				if (has_chimney == 2 && num_floors == 1) { // can't be a bathroom if there's a fireplace
					cube_t test_cube(*r);
					test_cube.expand_by_xy(floor_thickness);
					if (test_cube.intersects(get_fireplace())) continue;
				}
				float score(r->dx() + r->dy()); // starts as half the perimeter
				score *= (1.0 + 10.0*(max(count_num_int_doors(*r), 1U) - 1U)); // multiply by a large value if there are mult doors so we only choose this if there are no alternatives
				if (min_score == 0.0 || score < min_score) {cand_bathroom = (r - rooms.begin()); min_score = score;}
			}
		} // for r
	}
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {
		bool const is_basement(has_basement() && r->part_id == (int)basement_part_ix); // includes extended basement
		float light_amt(is_basement ? 0.0f : window_vspacing*r->get_light_amt()); // exterior light: multiply perimeter/area by window spacing to make unitless; none for basement rooms
		if (!is_house && r->is_hallway) {light_amt *= 2.0;} // double the light in office building hallways because they often connect to other lit hallways
		float const floor_height(r->is_sec_bldg ? r->dz() : window_vspacing); // secondary buildings are always one floor
		unsigned const num_floors(calc_num_floors_room(*r, floor_height, floor_thickness)), room_id(r - rooms.begin());
		point room_center(r->get_cube_center());

		// determine light pos and size for this stack of rooms
		float const dx(r->dx()), dy(r->dy());
		bool const room_dim(dx < dy); // longer room dim
		bool const must_be_bathroom(room_id == cand_bathroom && num_bathrooms == 0); // cand bathroom, and bathroom not already placed
		bool const is_parking_garage(r->get_room_type(0) == RTYPE_PARKING   ); // all floors should be parking garage
		bool const is_unfinished    (r->get_room_type(0) == RTYPE_UNFINISHED); //  // unfinished room, for example in a non-cube shaped office building
		bool const is_ext_basement(r->is_ext_basement());
		float light_size(def_light_size); // default size for houses
		unsigned const room_objs_start(objs.size());
		unsigned nx(1), ny(1); // number of lights in X and Y for this room

		if (!is_cube()) { // somewhat more lights for non-cube shaped building pie slices
			nx = max(1U, unsigned(0.4*dx/window_vspacing));
			ny = max(1U, unsigned(0.4*dy/window_vspacing));
		}
		else if (r->is_office) { // more lights for large offices; light size varies by office size; parking garages are handled later
			nx = max(1U, unsigned(0.5*dx/window_vspacing));
			ny = max(1U, unsigned(0.5*dy/window_vspacing));
			float const room_size(dx + dy); // normalized to office size
			light_size = max(0.015f*room_size, 0.67f*def_light_size);
		}
		else if (r->is_hallway) { // light size varies by hallway size
			float const room_size(min(dx, dy)); // normalized to hallway width
			light_size = max(0.06f*room_size, 0.67f*def_light_size);
		}
		if (r->is_sec_bldg) {
			if    (has_garage) {r->assign_all_to(RTYPE_GARAGE);}
			else if (has_shed) {r->assign_all_to(RTYPE_SHED);}
		}
		float const light_val(22.0*light_size);
		r->light_intensity = light_val*light_val/r->get_area_xy(); // average for room, unitless; light surface area divided by room surface area with some fudge constant
		cube_t light;
		set_light_xy(light, room_center, light_size, room_dim, light_shape);
		bool added_bathroom(0);
		float z(r->z1());
		if (!r->interior) {r->interior = (is_basement || get_part_for_room(*r).contains_cube_xy_no_adj(*r));} // AKA windowless; calculate if not already set
		// make chair colors consistent for each part by using a few variables for a hash
		colorRGBA chair_colors[12] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN};
		colorRGBA chair_color(chair_colors[(13*r->part_id + 123*tot_num_rooms + 617*mat_ix + 1367*num_floors) % 12]);
		light_ix_assign.next_room();
		rand_gen_t room_rgen(rgen); // shared across all floors of this room
		// select light color for this room
		colorRGBA color;
		if (r->is_ext_basement_conn()) {color = RED;}
		else if (is_house)          {color = get_light_color_temp(0.4);} // house - yellowish
		else if (is_parking_garage) {color = get_light_color_temp_range(0.2, 0.5, rgen);} // parking garage - yellow-white
		else if (r->is_office)      {color = get_light_color_temp(0.6);} // office - blueish
		else if (r->is_hallway)     {color = get_light_color_temp(0.6);} // office building hallway - blueish
		else                        {color = get_light_color_temp(0.5);} // small office - white

		// place objects on each floor for this room
		for (unsigned f = 0; f < num_floors; ++f, z += floor_height) {
			float const floor_zval(z + fc_thick);
			room_center.z = floor_zval;
			// top floor may have stairs connecting to upper stack
			bool const top_floor(f+1 == num_floors);
			bool const has_stairs_this_floor(r->has_stairs_on_floor(f));
			bool is_lit(0), has_light(1), light_dim(room_dim), wall_light(0), has_stairs(has_stairs_this_floor), top_of_stairs(has_stairs && top_floor);
			float light_delta_z(0.0);

			if (is_parking_garage) { // parking garage; added first because this sets the number of lights
				r->interior = 1;
				add_parking_garage_objs(rgen, *r, room_center.z, room_id, f, num_floors, nx, ny, light_delta_z);
				for (auto i = objs.begin() + room_objs_start; i != objs.end(); ++i) {i->flags |= RO_FLAG_INTERIOR;}
			}
			if ((!has_stairs && (f == 0 || top_floor) && interior->stairwells.size() > 1) || top_of_stairs) { // should this be outside the loop?
				// check for stairwells connecting stacked parts (is this still needed?); check for roof access stairs and set top_of_stairs=0
				for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
					if (!r->contains_cube_xy(*s)) continue; // stairs not in this room
					// Note: here we adjust stairs zval by floor_thickness to include stairs in the floor but not in the room above
					if (s->z1() + floor_thickness > r->z2()) continue; // stairs above the room
					if (s->z2() + floor_thickness < r->z1()) continue; // stairs below the room
					if (s->roof_access) {top_of_stairs = 0;}
					has_stairs = 1;
				} // for s
			}
			int light_obj_ix(-1);
			unsigned num_lights(r->num_lights), flags(0);
			float const light_z2(z + floor_height - fc_thick + light_delta_z);

			// motion detection lights for large office building office; limit to interior rooms so that we still have some lit rooms viewed through windows
			if (!is_house && has_pri_hall() && r->is_office && r->interior) {
				flags |= RO_FLAG_IS_ACTIVE; // leave unlit and enable motion detection for lights
			}
			else {
				// 50% of lights are on, 75% for top of stairs, 100% for non-basement hallways, 100% for parking garages
				is_lit = ((r->is_hallway && !is_basement) || is_parking_garage || ((rgen.rand() & (top_of_stairs ? 3 : 1)) != 0));

				if (!is_lit) { // check people and set is_lit if anyone is in this floor of this room
					for (person_t const &p : interior->people) {
						cube_t const bc(p.get_bcube());
						if (!bc.intersects_xy(*r)) continue; // person not in this room
						if (bc.z2() < light_z2 && bc.z1() + floor_height > light_z2) {is_lit = 1; break;} // on this floor
					}
				}
			}
			if (is_lit)     {flags |= RO_FLAG_LIT | RO_FLAG_EMISSIVE;}
			if (has_stairs) {flags |= RO_FLAG_RSTAIRS;}
			if (r->is_ext_basement_conn()) {flags |= RO_FLAG_EXTERIOR;} // flag as exterior since this light may reach the connected building
			// add a light to the ceiling of this room if there's space (always for top of stairs);
			set_cube_zvals(light, (light_z2 - light_thick), light_z2);
			valid_lights.clear();

			if (num_lights > 1) { // r->is_hallway or ext basement
				// hallway: place a light on each side (of the stairs if they exist), and also between stairs and elevator if there are both
				if (r->has_elevator && r->has_stairs == 255) {max_eq(num_lights, (2U + r->has_elevator));} // main hallway with elevator + stairs on all floors; 3+ lights
				min_eq(num_lights, 6U);
				float const offsets[6] = {0.0, -0.2, -0.3, -0.36, -0.4, -0.43}, steps[6] = {0.0, 0.4, 0.3, 0.24, 0.2, 0.172}; // indexed by num_lights-1
				float const hallway_len(r->get_sz_dim(light_dim));

				for (unsigned d = 0; d < num_lights; ++d) {
					float const delta((offsets[num_lights-1] + d*steps[num_lights-1])*hallway_len);
					cube_t hall_light(light);
					hall_light.translate_dim(light_dim, delta);
					try_place_light_on_ceiling(hall_light, *r, room_dim, fc_thick, 0, 0, valid_lights, rgen); // allow_rot=0, allow_mult=0
				}
			}
			else if (nx > 1 || ny > 1) { // office or parking garage with multiple lights
				float const xstep(dx/nx), ystep(dy/ny);
				vector3d const shrink(0.5*light.dx()*sqrt((nx - 1)/nx), 0.5*light.dy()*sqrt((ny - 1)/ny), 0.0);

				for (unsigned y = 0; y < ny; ++y) {
					for (unsigned x = 0; x < nx; ++x) {
						cube_t cur_light(light);
						cur_light.expand_by_xy(-shrink);
						cur_light.translate(point((-0.5f*dx + (x + 0.5)*xstep), (-0.5f*dy + (y + 0.5)*ystep), 0.0));
						try_place_light_on_ceiling(cur_light, *r, room_dim, fc_thick, 0, 0, valid_lights, rgen); // allow_rot=0, allow_mult=0
					}
				} // for y
			}
			else { // normal room with a single light
				if (is_house && is_basement && !is_ext_basement) { // house basement cylindrical lights only
					if (max(dx, dy) < 2.5*window_vspacing && min(dx, dy) < 2.0*window_vspacing) { // small rooms only
						try_place_light_on_wall(light, *r, room_dim, floor_zval, valid_lights, rgen);
						wall_light = !valid_lights.empty(); // if this fails, fall back to a ceiling light
					}
				}
				if (!wall_light) {try_place_light_on_ceiling(light, *r, room_dim, fc_thick, 1, 1, valid_lights, rgen);} // allow_rot=1, allow_mult=1
				if (!valid_lights.empty()) {light_obj_ix = objs.size();} // this will be the index of the light to be added later
			}
			rand_gen_t rgen_lights(rgen); // copy state so that we don't modify rgen
			unsigned const objs_start_inc_lights(objs.size());

			for (cube_t const &l : valid_lights) {
				bool dim(l.dx() < l.dy()), dir(0); // dir is only used for wall lights
				unsigned l_flags(flags);
				
				if (wall_light) {
					l_flags |= RO_FLAG_ADJ_HI; // flag as on wall
					dim     ^= 1; // use shorter dim
					dir      = (l.get_center_dim(dim) < r->get_center_dim(dim)); // direction the light is pointing, opposite the wall
				}
				else {l_flags |= RO_FLAG_NOCOLL;} // no collision detection for ceiling lights
				if (check_skylight_intersection(l)) {l_flags |= RO_FLAG_ADJ_TOP; has_skylight_light = 1;} // if attached to a skylight, draw top surface
				room_object_t light_obj(l, TYPE_LIGHT, room_id, dim, dir, l_flags, light_amt, light_shape, color);
				light_obj.obj_id = light_ix_assign.get_ix_for_light(l);
				unsigned const flicker_mod(is_parking_garage ? 50 : (is_ext_basement ? 20 : 0)); // 2% chance for parking garage, 5% chance for ext basement
				
				if (flicker_mod > 0 && (((rgen_lights.rand() + 3*f)%flicker_mod) == 13)) {light_obj.flags |= RO_FLAG_BROKEN;} // maybe make this a flickering light
				else if (is_ext_basement && valid_lights.size() == 1 && (rgen_lights.rand() & 7) == 0) { // broken ext basement light; not for hallways with multiple lights
					light_obj.flags |= RO_FLAG_BROKEN2;
					light_obj.flags &= ~RO_FLAG_LIT; // off by default
				}
				objs.emplace_back(light_obj);
			} // for l
			if (is_lit) {r->lit_by_floor |= (1ULL << (f&63));} // flag this floor as being lit (for up to 64 floors)
			if (is_parking_garage) continue; // generated above, done; no outlets or light switches
			if (is_unfinished    ) continue; // no objects for now; if adding objects later, need to make sure they stay inside the building bounds
			float tot_light_amt(light_amt); // unitless, somewhere around 1.0
			if (is_lit) {tot_light_amt += r->light_intensity;}
			bool const is_ground_floor(f == 0 && !is_basement && r->z1() <= ground_floor_z1), is_garage_or_shed(r->is_garage_or_shed(f));
			unsigned const objs_start(wall_light ? objs_start_inc_lights : objs.size()); // wall light counts as an object since it must be avoided
			rgen.rand_mix();

			if (r->is_ext_basement_conn()) { // room connecting extended basements of two buildings
				bool const conn_dim(r->interior == 4);
				// add blockers at both room ends to avoid placing objects there, since there's a door to the other building blocking it
				for (unsigned d = 0; d < 2; ++d) {
					cube_t blocker(*r);
					blocker.d[conn_dim][!d] = blocker.d[conn_dim][d] + (d ? -1.0 : 1.0)*2.0*get_wall_thickness();
					objs.emplace_back(blocker, TYPE_BLOCKER, room_id, conn_dim, d, RO_FLAG_INVIS);
				}
			}
			if (r->no_geom || is_garage_or_shed) {
				if (is_garage_or_shed) {
					if (r->get_room_type(0) == RTYPE_GARAGE) {
						room_center.z = add_flooring(*r, room_center.z, room_id, tot_light_amt, FLOORING_CONCRETE);
						add_garage_objs(rgen, *r, room_center.z, room_id, tot_light_amt);
					}
					// is there enough clearance between shelves and a car parked in the garage? there seems to be in all the cases I've seen
					add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
				}
				add_outlets_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor, is_basement);
				if (has_light) {add_light_switches_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor, is_basement);} // shed, garage, or hallway

				if (is_house && r->is_hallway) { // allow pictures, rugs, and light switches in the hallways of houses
					hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement);
					if (rgen.rand_bool()) {add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);} // 50% of the time; not all rugs will be placed
				}
				if (!is_house && r->is_hallway && *r == pri_hall) { // office building primary hallway
					add_pri_hall_objs(rgen, room_rgen, *r, room_center.z, room_id, tot_light_amt, f);
					if (is_ground_floor) {r->assign_to(RTYPE_LOBBY, f);} // first floor primary hallway, make it the lobby
				}
				continue; // no other geometry for this room
			}
			//if (has_stairs && !pri_hall.is_all_zeros()) continue; // no other geometry in office building base part rooms that have stairs
			unsigned const floor_mask(1<<f);
			bool added_tc(0), added_desk(0), added_obj(0), can_place_onto(0), no_whiteboard(0);
			bool is_bathroom(0), is_bedroom(0), is_kitchen(0), is_living(0), is_dining(0), is_storage(0), is_utility(0), no_plants(0), is_play_art(0);
			unsigned num_chairs(0);
			// unset room type if not locked on this floor during floorplanning; required to generate determinstic room geom
			if (!r->is_rtype_locked(f)) {r->assign_to(RTYPE_NOTSET, f);}

			// place room objects
			bool const allow_br(!is_house || must_be_bathroom || f > 0 || num_floors == 1 || (rgen.rand_float() < 0.33f*(added_living + (added_kitchen_mask&1) + 1))); // bed/bath
			bool is_office_bathroom(is_room_office_bathroom(*r, room_center.z, f)), has_fireplace(0);
			blockers.clear(); // clear for this new room
			
			if (has_chimney == 2 && !is_basement && is_entry_floor && !added_fireplace) { // handle fireplaces on the first floor
				has_fireplace = added_fireplace = maybe_add_fireplace_to_room(rgen, *r, blockers, room_center.z, room_id, tot_light_amt);
			}
			if (is_office_bathroom) { // bathroom is already assigned
				added_obj = is_bathroom = added_bathroom = no_whiteboard =
					add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement, added_bathroom_objs_mask); // add bathroom
			}
			else if (!is_house && f == 0) { // office building special first floor rooms; can be in a stacked part
				if (r->get_room_type(f) == RTYPE_UTILITY) {
					added_obj = no_whiteboard = no_plants = is_utility = add_office_utility_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
				else if (r->get_room_type(f) == RTYPE_SERVER) {
					added_obj = no_whiteboard = no_plants = add_server_room_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
			}
			// bedroom or bathroom case; need to check first floor even if must_be_bathroom
			if (!added_obj && allow_br && can_be_bedroom_or_bathroom(*r, f)) {
				// Note: num_bedrooms is summed across all floors, while num_bathrooms is per-floor
				bool const pref_sec_bath(is_house && num_bathrooms == 1 && num_bedrooms > 1 && rooms.size() >= 6 && !must_be_bathroom && !has_fireplace && can_be_bathroom(*r));
				float const bedroom_prob(pref_sec_bath ? 0.25 : 0.75), bathroom_prob((pref_sec_bath ? 2.0 : 1.0)*extra_bathroom_prob);
				// place a bedroom 75% of the time unless this must be a bathroom; if we got to the second floor and haven't placed a bedroom, always place it;
				// houses only, and must have a window (exterior wall)
				if (is_house && !must_be_bathroom && !is_basement && ((f > 0 && !added_bedroom) || rgen.rand_float() < bedroom_prob)) {
					// if haven't added a bedroom, force if last floor of last room (excluding the extended basement)
					bool const force(!added_bedroom && f+1 == num_floors && (r+1 == rooms.end() || (r+1)->is_ext_basement()));
					added_obj = can_place_onto = added_bedroom = is_bedroom =
						add_bedroom_objs(rgen, *r, blockers, chair_color, room_center.z, room_id, f, tot_light_amt, objs_start, is_lit, is_basement, force, light_ix_assign);
					if (is_bedroom) {r->assign_to(RTYPE_BED, f);}
					num_bedrooms += is_bedroom;
				}
				if (!added_obj && !has_fireplace && (must_be_bathroom || (can_be_bathroom(*r) && (num_bathrooms == 0 || rgen.rand_float() < bathroom_prob)))) {
					// bathrooms can be in both houses and office buildings
					added_obj = is_bathroom = added_bathroom =
						add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement, added_bathroom_objs_mask); // add bathroom
					if (is_bathroom) {r->assign_to(RTYPE_BATH, f);}
				}
			}
			if (!added_obj && r->is_office) { // add cubicles if this is a large office
				added_obj = no_whiteboard = create_office_cubicles(rgen, *r, room_center.z, room_id, tot_light_amt);
			}
			if (!added_obj && rgen.rand_float() < (is_basement ? 0.4 : (r->is_office ? 0.6 : (is_house ? 0.95 : 0.5)))) {
				// place a table and maybe some chairs near the center of the room if it's not a hallway;
				// 60% of the time for offices, 95% of the time for houses, and 50% for other buildings
				unsigned const num_tcs(add_table_and_chairs(rgen, *r, blockers, room_id, room_center, chair_color, 0.1, tot_light_amt));
				if (num_tcs > 0) {added_tc = added_obj = can_place_onto = 1; num_chairs = num_tcs - 1;}
				// on ground floor, try to make this a kitchen; not all houses will have a kitchen with this logic - maybe we need fewer bedrooms?
				if (!(added_kitchen_mask & floor_mask) && (!is_house || is_ground_floor) && !is_basement) { // office buildings can also have kitchens, even on non-ground floors
					if (added_tc || (is_house && (r+1) == rooms.end())) { // make it a kitchen if it's the last room in a house, even if there's no table
						if (add_kitchen_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, added_living)) {
							r->assign_to(RTYPE_KITCHEN, f);
							added_kitchen_mask |= floor_mask;
							is_kitchen = added_obj = 1;
						}
					}
				}
			}
			if (!added_obj && (is_basement || (r->is_office && r->interior && f == 0 /*&& r->z1() == ground_floor_z1*/)) && rgen.rand_bool()) {
				// if we haven't added any objects yet, and this room is an interior office on the first floor or basement, make it a storage room 50% of the time
				added_obj = no_whiteboard = is_storage = add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
				if (added_obj) {r->assign_to(RTYPE_STORAGE, f);}
			}
			if (!added_obj && (!is_basement || rgen.rand_bool())) { // try to place a desk if there's no table, bed, etc.
				added_obj = can_place_onto = added_desk = (is_house ?
					add_desk_to_room(rgen, *r, blockers, chair_color, room_center.z, room_id, f, tot_light_amt, objs_start, is_basement) :
					add_office_objs (rgen, *r, blockers, chair_color, room_center.z, room_id, f, tot_light_amt, objs_start, is_basement));
				if (added_obj && !has_stairs_this_floor) {r->assign_to((is_house ? (room_type)RTYPE_STUDY : (room_type)RTYPE_OFFICE), f);} // or other room type - may overwrite below
			}
			if (is_house && (added_tc || added_desk) && !is_kitchen && is_ground_floor) { // don't add second living room unless we added a kitchen and have enough rooms
				if ((!added_living && !r->has_center_stairs && rooms.size() >= 8 && (added_kitchen_mask || rgen.rand_bool())) || is_room_adjacent_to_ext_door(*r, 1)) { // front_door_only=1
					// add a living room on the ground floor if it has a table or desk but isn't a kitchen
					added_living = is_living = add_livingroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
					if (is_living) {r->assign_to(RTYPE_LIVING, f);}
				}
			}
			if (is_house && added_tc && num_chairs > 0 && !is_living && !is_kitchen) { // room with table and chair that's not a kitchen
				if (is_ground_floor) { // dining room, must be on the ground floor
					if (light_obj_ix >= 0) { // handle dining room light (assume there is only one): extend downward and make it a sphere
						assert((unsigned)light_obj_ix < objs.size());
						room_object_t &light(objs[light_obj_ix]);
						light.shape = SHAPE_SPHERE;
						light.z2() += 0.5f*light.dz();
						light.z1() -= 0.22f*(light.dx() + light.dy());
					}
					if (!added_dining) {add_diningroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);} // only one room is the primary dining room
					r->assign_to(RTYPE_DINING, f);
					is_dining = added_dining = 1;
				}
				else if (!added_library && add_library_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement)) { // add library, at most one
					r->assign_to(RTYPE_LIBRARY, f);
					added_library = 1;
				}
			}
			if (!is_house && r->is_office && !no_whiteboard && (rgen.rand() % (pri_hall.is_all_zeros() ? 30U : max(50U, (unsigned)interior->rooms.size()))) == 0) {
				// office, no cubicles or bathroom - try to make it a library (in rare cases)
				if (add_library_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement)) {r->assign_to(RTYPE_LIBRARY, f);}
			}
			if (can_place_onto) { // an object was placed (table, desk, counter, etc.), maybe add a book or bottle on top of it
				place_objects_onto_surfaces(rgen, *r, room_id, tot_light_amt, objs_start, f, is_basement);
			}
			if (is_house) { // place house-specific items
				if (!is_bathroom && !is_kitchen && rgen.rand_float() < (is_basement ? 0.25 : 0.8)) { // place bookcase 80% of the time, but not in bathrooms or kitchens
					rand_gen_t rgen2(rgen); // copy so that rgen isn't updated in the call below
					add_bookcase_to_room(rgen2, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
				}
				if (!has_stairs && (rgen.rand()&3) <= (added_tc ? 0 : 2) && !is_kitchen) { // maybe add a rug, 25% of the time if there's a table and 75% of the time otherwise
					add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
			}
			else if (has_pri_hall() && r->part_id == 0 && f == 0 && added_desk) { // office building part with primary hallway, first floor of any part
				add_office_door_sign(rgen, *r, room_center.z, room_id, tot_light_amt);
			}
			bool const room_type_was_not_set(r->get_room_type(f) == RTYPE_NOTSET);

			if (room_type_was_not_set) { // attempt to assign it with an optional room type
				if (is_ground_floor && is_room_adjacent_to_ext_door(*r)) { // entryway/lobby if on ground floor, has exterior door, and unassigned
					r->assign_to((is_house ? (room_type)RTYPE_ENTRY : (room_type)RTYPE_LOBBY), f); // office building lobby can have a whiteboard - is that okay?
				}
				else if (!is_house) {r->assign_to(RTYPE_OFFICE, f);} // any unset room in an office building is an office
				// else house
				else if (has_stairs && !is_basement) {} // will be marked as RTYPE_STAIRS below
				else if ((!added_obj || is_basement) && f == 0 && !added_laundry && !has_fireplace &&
					add_laundry_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, added_bathroom_objs_mask))
				{
					r->assign_to(RTYPE_LAUNDRY, f);
					added_laundry = 1;
				}
				else if (!added_obj && !has_fireplace) { // make it a storage room until we add some other room type that it can be
					add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
					r->assign_to(RTYPE_STORAGE, f);
					is_storage = 1; // mark it as a storage room whether or not we've added anything to it
				}
				else if (is_basement) {r->assign_to(RTYPE_CARD, f);} // basement card room
				else { // unassigned room of house on upper floor with added object/table
					// this case is relatively rare, and we've already added a table, so it's too late to make this a bedroom/bathroom if can_be_bedroom_or_bathroom(*r, f)
					r->assign_to((rgen.rand_bool() ? (room_type)RTYPE_PLAY : (room_type)RTYPE_ART), f); // play room or art room
					is_play_art = 1;
				}
			}
			if (is_house && is_basement && !added_basement_utility && !has_stairs && (is_storage || room_type_was_not_set) && rgen.rand_bool()) {
				// basement laundry, storage, or card room; should this be placed before adding boxes to the floor of storage rooms?
				added_basement_utility = is_utility = no_plants = add_basement_utility_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				if (added_basement_utility) {r->assign_to(RTYPE_UTILITY, f);}
			}
			if (!is_bathroom && !is_bedroom && !is_kitchen && !is_storage && !no_plants && !is_basement) { // add potted plants to some room types
				// 0-2 for living/dining rooms, 50% chance for houses, 25% (first floor) / 10% (other floors) chance for offices
				unsigned const num(is_house ? (rgen.rand() % ((is_living || is_dining) ? 3 : 2)) : ((rgen.rand()%((f == 0) ? 4 : 10)) == 0));
				if (num > 0) {add_plants_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, num);}
			}
			add_outlets_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor, is_basement);
			if (has_light) {add_light_switches_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor, is_basement);} // add light switch if room has a light
			
			if (!r->is_hallway) { // no vents in hallways; vents use orig floor zval, not adjusted for bathroom tile floor
				if (is_house) {add_ceil_vent_to_room(rgen, *r, floor_zval, room_id, objs_start_inc_lights );} // house vents
				else          {add_wall_vent_to_room(rgen, *r, floor_zval, room_id, objs_start, is_utility);} // office building vents
			}
			// pictures and whiteboards must not be placed behind anything, excluding trashcans; so we add them here
			bool const can_hang((is_house || !(is_bathroom || is_kitchen || no_whiteboard)) && !is_storage); // no whiteboards in office bathrooms or kitchens
			bool const was_hung(can_hang && hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement));

			if (is_bathroom || is_kitchen || rgen.rand_float() < 0.8) { // 80% of the time, always in bathrooms and kitchens
				add_trashcan_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, (was_hung && !is_house)); // no trashcans on same wall as office whiteboard
			}
			if (is_bedroom || is_living || is_dining || is_play_art) {
				add_floor_clutter_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
			}
			if (is_house && !(is_bathroom || is_kitchen || is_storage) && rgen.rand_float() < ((f > 0) ? 0.15 : 0.25)) {
				unsigned const max_num(is_bedroom ? 1 : 2);
				add_boxes_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, max_num); // place boxes in this room
			}
			if (has_stairs_this_floor && r->get_room_type(f) == RTYPE_NOTSET) {r->assign_to(RTYPE_STAIRS, f);} // default to stairs if not set above
		} // for f (floor)
		if (added_bathroom) {++num_bathrooms;}

		if (r->interior) { // tag objects as interior if room is interior
			for (auto i = objs.begin() + room_objs_start; i != objs.end(); ++i) {i->flags |= RO_FLAG_INTERIOR;}
		}
	} // for r (room)
	if (is_house) {interior->assign_master_bedroom(window_vspacing, floor_thickness);}

	if (is_rotated() || !is_cube()) {} // skip for rotated and non-cube buildings, since toilets, etc. may not be placed
	else if (num_bathrooms == 0) { // can happen, but very rare
		cout << "no bathroom in building " << bcube.xc() << " " << bcube.yc() << endl;
		if (cand_bathroom < rooms.size()) {cout << "cand bathroom was at " << rooms[cand_bathroom].str() << endl;}
	}
	else {
		if (!(added_bathroom_objs_mask & PLACED_TOILET)) {cout << "no toilet in building " << bcube.xc() << " " << bcube.yc() << endl;}
		if (!(added_bathroom_objs_mask & PLACED_SINK  )) {cout << "no sink in building "   << bcube.xc() << " " << bcube.yc() << endl;}
		//if (is_house && !(added_bathroom_objs_mask & (PLACED_TUB | PLACED_SHOWER))) {cout << "no bathtub or shower in building " << bcube.xc() << " " << bcube.yc() << endl;} // common
	}
	// add trim + window coverings; must be done after room assignment; not implemented for rotated buildings
	if (!is_rotated()) {add_window_trim_and_coverings(0, 1, 1);} // add_trim=0, add_coverings=1, add_ext_sills=1
	if (is_house && has_basement()) {add_basement_electrical_house(rgen);}
	if (is_house && has_basement_pipes) {add_house_basement_pipes (rgen);}
	if (has_attic()) {add_attic_objects(rgen);}
	maybe_add_fire_escape  (rgen);
	add_balconies          (rgen);
	add_exterior_door_items(rgen);
	add_extra_obj_slots(); // needed to handle balls taken from one building and brought to another
	add_stairs_and_elevators(rgen); // the room objects - stairs and elevators have already been placed within a room
	objs.shrink_to_fit(); // Note: currently up to around 15K objs max for large office buildings
	interior->room_geom->light_bcubes.resize(light_ix_assign.get_next_ix()); // allocate but don't fill un until needed
	// randomly vary wood color for this building
	colorRGBA &wood_color(interior->room_geom->wood_color);
	float const luminance(rgen.rand_uniform(0.4, 1.6));
	for (unsigned i = 0; i < 3; ++i) {wood_color[i] = luminance*WOOD_COLOR[i]*rgen.rand_uniform(0.9, 1.1);}
	wood_color.set_valid_color();
	max_eq(wood_color.R, max(wood_color.G, wood_color.B)); // make sure wood isn't blue or green tinted
}

void building_interior_t::assign_master_bedroom(float window_vspacing, float floor_thickness) {
	float best_area(0.0);
	unsigned master_br(0), mbr_floor(0);

	for (auto r = rooms.begin(); r != rooms.end(); ++r) {
		unsigned const num_floors(calc_num_floors_room(*r, window_vspacing, floor_thickness));
		
		for (unsigned f = 0; f < num_floors; ++f) {
			if (r->get_room_type(f) != RTYPE_BED) continue;
			float const area(r->get_area_xy());
			if (area > best_area) {master_br = (r - rooms.begin()); mbr_floor = f; best_area = area;} // Note: prioritizes lower floors
		}
	} // for r
	if (best_area > 0) {rooms[master_br].assign_to(RTYPE_MASTER_BED, mbr_floor);}
}

// *** Exterior Objects - not really room objects ***

void building_t::maybe_add_fire_escape(rand_gen_t &rgen) {
	if (!is_house) return; // houses only for now
	// our hard-coded fire escape model is designed for a 5 story building; but the max number of floors for a 'house' is 5-6 anyway, which makes them relatively rare
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_FESCAPE)) return;
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fe_height(4.25*window_vspacing);

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		unsigned const num_floors(calc_num_floors(*p, window_vspacing, floor_thickness));
		if (num_floors != 5 && num_floors != 6) continue; // not 5-6 stories
		unsigned const pref_dim_dir(rgen.rand() & 3);
		// it's uncommon to get here, so we only check if the model size here
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FESCAPE)); // D, W, H
		float const fe_hwidth(0.5*fe_height*sz.y/sz.z), fe_depth(fe_height*sz.x/sz.z);

		for (unsigned d = 0; d < 4; ++d) {
			unsigned const dd((d + pref_dim_dir) & 3);
			bool const dim(dd >> 1), dir(dd & 1);
			if (p->d[dim][dir] != bcube.d[dim][dir]) continue; // not on the building bcube - could intersect another part, porch, etc.
			if (p->get_sz_dim(!dim) < 3.0*fe_hwidth) continue; // wall is too narrow
			cube_t fe_bc;
			set_cube_zvals(fe_bc, p->z1(), (p->z1() + fe_height));
			set_wall_width(fe_bc, rgen.rand_uniform((p->d[!dim][0] + 1.2*fe_hwidth), (p->d[!dim][1] - 1.2*fe_hwidth)), fe_hwidth, !dim);
			fe_bc.d[dim][0] = fe_bc.d[dim][1] = p->d[dim][dir];
			fe_bc.d[dim][dir] += (dir ? 1.0 : -1.0)*fe_depth;
			if (has_bcube_int_no_adj(fe_bc, parts))              continue; // check for intersection with other parts, in particular the chimney and fireplace
			if (has_driveway() && fe_bc.intersects_xy(driveway)) continue; // skip if intersects driveway or garage
			if (is_room_adjacent_to_ext_door(fe_bc, 0))          continue; // check exterior doors; front_door_only=0
			interior->room_geom->objs.emplace_back(fe_bc, TYPE_FESCAPE, 0, dim, dir, 0, 1.0, SHAPE_CUBE, BLACK); // room_id=0
			details.emplace_back(fe_bc, DETAIL_OBJ_COLLIDER);
			union_with_coll_bcube(fe_bc);
			return; // success/done
		} // for d
	} // for p
}

void building_t::add_balconies(rand_gen_t &rgen) {
	if (!is_house || !has_room_geom()) return; // houses only for now
	if (rgen.rand_bool()) return; // only add balconies to 50% of houses
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	float const balcony_depth(floor_spacing*rgen.rand_uniform(0.45, 0.6)); // constant per house
	float const room_min_z2(ground_floor_z1 + 1.5*floor_spacing); // > 1 floor
	unsigned const balcony_style(rgen.rand()); // shared across all balconies of this house
	unsigned const max_balconies(1 + rgen.rand_bool()); // 1-2 per house
	unsigned num_balconies(0);
	auto &objs(interior->room_geom->objs);
	cube_t avoid1, avoid2;
	if (!objs.empty() && objs.back().type == TYPE_FESCAPE) {avoid1 = objs.back();}

	// find suitable rooms for balconies; since room walls will never intersect windows, we can make the balcony the same width to avoid intersecting windows
	for (auto room = interior->rooms.begin(); room != interior->rooms.end(); ++room) {
		if (room->interior)           continue; // no windows
		if (room->is_sec_bldg)        continue; // no garage or shed, even if it's multiple stories tall
		if (room->z2() < room_min_z2) continue; // ground floor only
		if (rgen.rand_float() < 0.75) continue;
		float const balcony_z1(room->z2() - floor_spacing /*+ get_fc_thickness()*/); // floor level of top floor of room
		unsigned const floor_ix(room->get_floor_containing_zval((balcony_z1 + 0.1*floor_spacing), floor_spacing));
		unsigned const room_id(room - interior->rooms.begin());
		room_type const rtype(room->get_room_type(floor_ix));
		if (rtype == RTYPE_BATH) continue; // no bathroom balconies as that would be weird
		assert(room->part_id < parts.size());
		cube_t const &part(parts[room->part_id]);
		bool added(0);

		for (unsigned dim = 0; dim < 2 && !added; ++dim) {
			for (unsigned dir = 0; dir < 2 && !added; ++dir) {
				if (classify_room_wall(*room, balcony_z1, dim, dir, 1) != ROOM_WALL_EXT) continue; // not fully exterior wall, skip
				cube_t balcony(*room);
				balcony.z1() = balcony_z1;
				balcony.d[dim][!dir]  = balcony.d[dim][dir]; // abuts the exterior wall of the room
				balcony.d[dim][ dir] += (dir ? 1.0 : -1.0)*balcony_depth; // extend outward from the house
				if (!avoid1.is_all_zeros() && avoid1.intersects(balcony)) continue; // check for fire escape intersection
				if (!avoid2.is_all_zeros() && avoid2.intersects(balcony)) continue; // check for previous balcony intersection
				if (check_cube_intersect_non_main_part(balcony))          continue; // porch roof, porch support, and chimney, etc.
				cube_t balcony_ext_down(balcony);
				balcony_ext_down.z1() = ground_floor_z1; // extend down to the ground
				bool part_int(0);

				for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // check for any intersecting parts (likely below this one, for stacked parts)
					if ((p - parts.begin()) != room->part_id && p->intersects_no_adj(balcony_ext_down)) {part_int = 1; break;}
				}
				if (part_int) continue;
				balcony.z2() -= 0.6*floor_spacing; // reduce wall height to 40%
				balcony.expand_in_dim(!dim, wall_thickness); // expand slightly to include window frame and merge adj balcony shared walls
				max_eq(balcony.d[!dim][0], (part.d[!dim][0] + 0.25f*wall_thickness)); // clamp slightly smaller than the containing part in !dim
				min_eq(balcony.d[!dim][1], (part.d[!dim][1] - 0.25f*wall_thickness));
				// if the space below the balcony is blocked by something, flag as hanging so that we don't try to draw vertical supports later
				cube_t area_below(balcony);
				set_cube_zvals(area_below, ground_floor_z1, balcony.z1());
				bool const hanging((!avoid1.is_all_zeros() && avoid1.intersects(area_below)) ||
					(!avoid2.is_all_zeros() && avoid2.intersects(area_below)) ||
					check_cube_intersect_non_main_part(area_below) ||
					(has_driveway() && driveway.intersects_xy(area_below)));
				room_object_t balcony_obj(balcony, TYPE_BALCONY, room_id, dim, dir, (hanging ? RO_FLAG_HANGING : 0), 1.0, SHAPE_CUBE, WHITE);
				balcony_obj.obj_id = balcony_style; // set so that we can select from multiple balcony styles
				objs.push_back(balcony_obj);
				avoid2 = balcony; // shouldn't need to consider area_below
				avoid2.expand_by_xy(wall_thickness);
				// maybe add plants to balconies; note that they won't be properly lit since plants use indoor lighting,
				// and the plants won't be drawn when the player is outside the building; this should be okay because an empty pot works on a balcony as well
				unsigned const num_plants(rgen.rand() % 3); // 0-2

				if (num_plants > 0) { // place plants
					cube_t cubes[4], plant_avoid;
					get_balcony_cubes(objs.back(), cubes);
					float const wall_width(cubes[1].get_sz_dim(dim)); // use front wall
					cube_t floor_inner(cubes[0]);
					floor_inner.expand_by_xy(-wall_width);

					for (unsigned n = 0; n < num_plants; ++n) {
						if (place_plant_on_obj(rgen, floor_inner, room_id, 1.0, plant_avoid)) {plant_avoid = objs.back();} // tot_light_amt=1.0
					}
				}
				// else place table and chairs?
				if (!hanging) { // add colliders for vertical supports
					cube_t pillar[2];
					get_balcony_pillars(balcony_obj, ground_floor_z1, pillar);
					for (unsigned d = 0; d < 2; ++d) {details.emplace_back(pillar[d], DETAIL_OBJ_COLL_SHAD);} // collider + shadow caster
				}
				union_with_coll_bcube(balcony_ext_down);
				++num_balconies;
				added = 1;
			} // for dir
		} // for dim
		if (num_balconies == max_balconies) break; // done
	} // for room
}

void building_t::add_extra_obj_slots() {
	assert(has_room_geom());
	vect_room_object_t &objs(interior->room_geom->objs);
	if (objs.empty()) return; // if there are no objects (empty building), don't allocate any extra slots
	unsigned num_slots(0);
	for (auto i = objs.begin(); i != objs.end(); ++i) {num_slots += (i->type == TYPE_BLOCKER);}
	if (num_slots >= 10) return;
	// make sure there are at least 10 blockers that will create free slots when adding dynamic objects
	float const v(get_wall_thickness()); // arbitrary
	cube_t const c(all_zeros, vector3d(v, v, v)); // arbitrary, doesn't have to be inside building, only needs to be strictly normalized
	for (unsigned n = num_slots; n < 20; ++n) {objs.emplace_back(c, TYPE_BLOCKER, 0, 0, 0, (RO_FLAG_INVIS | RO_FLAG_NOCOLL));}
}

// *** Wall and Door Trim ***

void building_t::add_wall_and_door_trim_if_needed() {
	if (!interior || (interior->walls[0].empty() && interior->walls[1].empty())) return; // no interior or walls
	if (!interior->room_geom->trim_objs.empty()) return; // trim already generated
	add_wall_and_door_trim();
	interior->room_geom->trim_objs.shrink_to_fit();
}

void cut_trim_around_doors(vector<tquad_with_ix_t> const &doors, vect_cube_t &trim_cubes, float door_expand, bool dim) {
	for (auto d = doors.begin(); d != doors.end(); ++d) {
		cube_t door(d->get_bcube());
		bool const door_dim(door.dy() < door.dx());
		if (door_dim != bool(dim)) continue;
		door.expand_in_dim(door_dim, door_expand); // expand to nonzero area; use a larger expand to account for distance door is offset away from ext wall
		subtract_cube_from_cubes(door, trim_cubes); // subtract this door from current trim cubes by clipping in XY
	}
}

void building_t::add_wall_and_door_trim() { // and window trim
	//highres_timer_t timer("Add Wall And Door Trim");
	assert(has_room_geom());
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness), wall_thickness(get_wall_thickness());
	float const trim_height(get_trim_height()), trim_thickness(get_trim_thickness()), expand_val(2.0*trim_thickness), door_trim_offset(0.025*window_vspacing);
	float const door_trim_exp(2.0*trim_thickness + 0.5*wall_thickness), door_trim_width(0.5*wall_thickness), floor_to_ceil_height(window_vspacing - floor_thickness);
	float const trim_toler(0.1*trim_thickness); // required to handle wall intersections that were calculated with FP math and may misalign due to FP rounding error
	float const ext_wall_toler(0.01*trim_thickness); // required to prevent z-fighting when AA is disabled
	unsigned const flags(RO_FLAG_NOCOLL);
	// ceiling trim disabled for large office buildings with outside corners because there's a lot of trim to add, and outside corners don't join correctly;
	// ceiling trim also disabled for non-houses (all office buildings), because it doesn't really work with acoustic paneling
	bool const has_outside_corners(!is_house && !pri_hall.is_all_zeros()), has_ceil_trim(!has_outside_corners && is_house);
	colorRGBA const &trim_color(is_house ? WHITE : DK_GRAY);
	vect_room_object_t &objs(interior->room_geom->trim_objs);
	vect_cube_t trim_cubes;

	for (auto d = interior->door_stacks.begin(); d != interior->door_stacks.end(); ++d) { // vertical strips on each side + strip on top of interior doors
		if (d->on_stairs) continue; // no frame for stairs door, skip
		cube_t trim(*d);
		trim.expand_in_dim(d->dim, door_trim_exp);
		trim.z2() -= 0.1*trim_toler; // shift top down oh so slightly to prevent z-fighting with top of wall when drawn under a skylight

		for (unsigned side = 0; side < 2; ++side) { // left/right of door
			trim.d[!d->dim][0] = d->d[!d->dim][side] - (side ? trim_thickness : door_trim_width);
			trim.d[!d->dim][1] = d->d[!d->dim][side] + (side ? door_trim_width : trim_thickness);
			bool const draw_top(check_skylight_intersection(trim)); // draw top edge of trim for top floor?
			unsigned const flags2(flags | RO_FLAG_ADJ_BOT | (draw_top ? 0 : RO_FLAG_ADJ_TOP));
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, d->dim, side, flags2, 1.0, SHAPE_TALL, trim_color); // abuse tall flag
		}
		// add trim at top of door
		unsigned const num_floors(calc_num_floors(*d, window_vspacing, floor_thickness));
		float z(d->z1() + floor_to_ceil_height);
		trim.d[!d->dim][0] = d->d[!d->dim][0] + trim_thickness;
		trim.d[!d->dim][1] = d->d[!d->dim][1] - trim_thickness;

		for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
			set_cube_zvals(trim, z-trim_thickness, z); // z2=ceil height
			bool const draw_top(f+1 == num_floors && check_skylight_intersection(trim)); // draw top edge of trim for top floor?
			unsigned const flags2(flags | (draw_top ? 0 : RO_FLAG_ADJ_TOP));
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, d->dim, 0, flags2, 1.0, SHAPE_SHORT, trim_color);
		}
	} // for d
	for (auto d = doors.begin(); d != doors.end(); ++d) { // exterior doors
		if (d->type == tquad_with_ix_t::TYPE_RDOOR) continue; // roof access door - requires completely different approach to trim and has not been implemented
		cube_t door(d->get_bcube()), trim(door);
		bool const dim(door.dy() < door.dx()), garage_door(d->type == tquad_with_ix_t::TYPE_GDOOR);
		trim.expand_in_dim(dim, door_trim_exp);
		bool dir(0);
		unsigned ext_flags(flags);
		colorRGBA const &ext_trim_color(garage_door ? WHITE : door_color); // garage doors are always white
		float const trim_width(((garage_door && has_int_garage) ? 1.5 : 1.0)*door_trim_width); // interior garage doors have thicker trim

		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (!i->intersects_no_adj(trim)) continue;
			trim.intersect_with_cube_xy(*i); // clip to containing part
			dir = (i->get_center_dim(dim) < trim.get_center_dim(dim));

			if (!is_cube()) { // handle polygon sides
				vect_point const &points(get_part_ext_verts(i - parts.begin()));

				for (auto v = points.begin(); v != points.end(); ++v) {
					point const &p1(*v), &p2((v == points.begin()) ? points.back() : *(v-1));
					if (fabs(p1[dim] - p2[dim]) > 0.1*wall_thickness) continue; // not nearly parallel in dim
					cube_t const edge_bc(p1, p2);
					if (!edge_bc.intersects_xy(trim)) continue; // wrong edge
					trim.d[dim][dir] = 0.5*(p1[dim] + p2[dim]); // clip to edge pos
					break; // should be only one
				} // for i
			}
			trim.d[dim][dir] -= (dir ? -1.0 : 1.0)*door_trim_offset; // move to the same offset for door
			ext_flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
			break;
		}
		unsigned const side_trim_flags(ext_flags | RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP);

		// add side trim; we could split this into interior vs. exterior geom, but there's no static assignment that works for both open and closed doors,
		// and it's not very noticeable anyway due to how thin the edges of the door trim are
		for (unsigned side = 0; side < 2; ++side) { // left/right of door
			trim.d[!dim][0] = door.d[!dim][side] - (side ? trim_width : 0.0);
			trim.d[!dim][1] = door.d[!dim][side] + (side ? 0.0 : trim_width);
			assert(trim.is_strictly_normalized());
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, side, side_trim_flags, 1.0, SHAPE_TALL, ext_trim_color); // abuse tall flag
		}
		// add trim at bottom of door for threshold
		trim.d[!dim][0] = door.d[!dim][0];
		trim.d[!dim][1] = door.d[!dim][1];
		set_cube_zvals(trim, door.z1()+fc_thick, door.z1()+fc_thick+2.0*trim_thickness); // floor height
		objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, (ext_flags | RO_FLAG_ADJ_BOT), 1.0, SHAPE_SHORT, ext_trim_color);

		if (d->type == tquad_with_ix_t::TYPE_HDOOR || d->is_building_door() || garage_door) { // add trim at top of exterior door, houses and office buildings
			set_cube_zvals(trim, door.z2()-0.03*door.dz(), door.z2()); // ends at top of door texture; see logic in clip_door_to_interior()
		}
		if (d->is_building_door()) { // different logic for building doors
			ext_flags = flags; // unlike hdoors, need to draw the back face to hide the gap betweeen ceiling and floor above
			trim.d[dim][dir] += (dir ? -1.0 : 1.0)*0.2*door_trim_offset; // minor shift back toward building to prevent z-fighting
		}
		objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_SHORT, ext_trim_color); // top of door
	} // for d
	for (unsigned dim = 0; dim < 2; ++dim) { // add horizontal strips along each wall at each floor/ceiling
		for (auto w = interior->walls[dim].begin(); w != interior->walls[dim].end(); ++w) {
			cube_t trim(*w);
			trim.expand_in_dim(dim, trim_thickness);

			if (has_outside_corners) { // handle outside corners of office building hallway intersections
				for (auto W = interior->walls[!dim].begin(); W != interior->walls[!dim].end(); ++W) { // check walls in other dim for an outside corner
					for (unsigned d = 0; d < 2; ++d) {
						if (W->z1() > w->z2() || W->z2() < w->z1()) continue; // no z overlap, wrong stack
						if (W->d[!dim][0] > w->d[!dim][d]+trim_toler || W->d[!dim][1] < w->d[!dim][d]-trim_toler) continue; // not adjacent/overlapping
						if (W->d[ dim][0] < w->d[ dim][0]-trim_toler && W->d[ dim][1] > w->d[ dim][1]+trim_toler) continue; // skip T junctions
						trim.d[!dim][d] = W->d[!dim][d] + (d ? 1.0 : -1.0)*trim_thickness; // expand to cover gap at outside corners of hallway walls
					}
				} // for W
			}
			unsigned const num_floors(calc_num_floors(*w, window_vspacing, floor_thickness));
			// snap to the nearest floor to handle short walls due to cut out stairs
			float const ground_wall_z1(bcube.z1() + fc_thick);
			float z(ground_wall_z1 + window_vspacing*round_fp((w->z1() - ground_wall_z1)/window_vspacing));

			for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
				set_cube_zvals(trim, z, z+trim_height); // starts at floor height
				bool ext_dirs[2] = {0,0};

				if (z < ground_floor_z1 && has_ext_basement() && !get_basement().intersects(trim)) { // check for exterior wall of extended basement
					for (unsigned dir = 0; dir < 2; ++dir) { // for each side of wall
						cube_t test_cube(trim);
						test_cube.d[dim][!dir] = w->d[dim][dir];
						point const center(test_cube.get_cube_center());
						ext_dirs[dir] = !interior->point_in_ext_basement_room(center);
					}
				}
				unsigned const trim_flags(flags | (ext_dirs[0] ? RO_FLAG_ADJ_LO : 0) | (ext_dirs[1] ? RO_FLAG_ADJ_HI : 0)); // disable exterior faces
				objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, 0, trim_flags, 1.0, SHAPE_CUBE, trim_color); // floor trim
				if (!has_ceil_trim) continue;
				trim.z2() = z + floor_to_ceil_height; // ceil height
				trim.z1() = trim.z2() - trim_height;

				for (unsigned dir = 0; dir < 2; ++dir) { // for each side of wall
					if (ext_dirs[dir]) continue; // skip
					cube_t ceil_trim(trim);
					ceil_trim.d[dim][!dir] = w->d[dim][dir];
					objs.emplace_back(ceil_trim, TYPE_WALL_TRIM, 0, dim, dir, flags, 1.0, SHAPE_ANGLED, trim_color); // ceiling trim
				}
			} // for f
		} // for w
	} // for d
	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // add trim for exterior walls
		if (is_basement(i)) continue; // skip basement walls because they're bare concrete
		bool const is_sec_bldg(i == get_real_parts_end());
		unsigned const num_floors(is_sec_bldg ? 1 : calc_num_floors(*i, window_vspacing, floor_thickness));

		if (!is_cube()) { // add floor trim around the interior of the exterior walls for each floor
			float const toler(0.1*wall_thickness); // add a bit of tolerance to account for FP error
			vect_point points(get_part_ext_verts(i - parts.begin())); // deep copy so that we can modify it
			expand_convex_polygon_xy(points, i->get_cube_center(), -trim_thickness); // shrink by trim thickness
			float z(i->z1() + fc_thick);

			for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
				for (auto p = points.begin(); p != points.end(); ++p) {
					point const &p1(*p), &p2((p+1 == points.end()) ? points.front() : *(p+1));
					cube_t trim(p1, p2); // bcube of the edge
					set_cube_zvals(trim, z, z+trim_height); // starts at floor height
					bool added(0);

					for (unsigned d = 0; d < 2; ++d) { // check for horizontal or vertical dim; maybe need to clip to door
						if (fabs(p1[d] - p2[d]) > toler) continue; // not aligned in this axis
						float const side_pos(0.5*(p1[d] + p2[d]));
						bool const dir(side_pos > i->get_center_dim(d));
						trim.d[d][ dir] = side_pos;
						trim.d[d][!dir] = side_pos + (dir ? -1.0 : 1.0)*trim_thickness;
						trim_cubes.clear();
						trim_cubes.push_back(trim); // start with entire length
						if (f == 0) {cut_trim_around_doors(doors, trim_cubes, (expand_val + wall_thickness), d);} // first floor, cut out areas for ext doors
						for (cube_t const &c : trim_cubes) {objs.emplace_back(c, TYPE_WALL_TRIM, 0, d, dir, flags, 1.0, SHAPE_CUBE, trim_color);}
						added = 1;
						break;
					}
					if (added) continue;

					for (unsigned d = 0; d < 2; ++d) {
						if (trim.get_sz_dim(d) < trim_thickness) {trim.expand_in_dim(d, 0.1*trim_thickness);} // expand slightly to make nonzero area
					}
					bool const dim(p1.x < p2.x), dir(p1.y < p2.y); // encode edge X/Y signs in dim and dir
					objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, flags, 1.0, SHAPE_CYLIN, trim_color); // floor trim - cylin section
				} // for p
			} // for f
			continue;
		}
		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				cube_t trim(*i);
				trim.d[dim][!dir]  = i->d[dim][dir] + (dir ? -1.0 : 1.0)*trim_thickness;
				trim.d[dim][ dir] += (dir ? -1.0 : 1.0)*ext_wall_toler; // slight bias away from the exterior wall
				unsigned const ext_flags(flags | (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));
				float z(i->z1() + fc_thick);

				for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
					set_cube_zvals(trim, z, z+trim_height); // starts at floor height
					trim_cubes.clear();
					trim_cubes.push_back(trim); // start with entire length
					float const room_height(is_sec_bldg ? (i->dz() - floor_thickness) : floor_to_ceil_height);
					float const ceil_trim_z2(z + room_height), ceil_trim_z1(ceil_trim_z2 - trim_height); // ceil height

					for (auto j = parts.begin(); j != get_real_parts_end(); ++j) { // clip against other parts
						if (j == i) continue; // skip self
						cube_t clip_cube(*j);
						clip_cube.expand_in_dim(dim, expand_val); // expand to clip trim on the other side of the split wall
						subtract_cube_from_cubes(clip_cube, trim_cubes); // subtract this part from current trim cubes by clipping in XY
					}
					if (has_ceil_trim && is_house) { // houses have shorter doors and ceiling trim extends above the door, so draw full range
						for (auto c = trim_cubes.begin(); c != trim_cubes.end(); ++c) {
							cube_t trim(*c); // copy so that we can modify it
							set_cube_zvals(trim, ceil_trim_z1, ceil_trim_z2);
							objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, !dir, flags, 1.0, SHAPE_ANGLED, trim_color); // ceiling trim
						}
					}
					if (f == 0) {cut_trim_around_doors(doors, trim_cubes, (expand_val + wall_thickness), dim);} // first floor, cut out areas for ext doors

					for (cube_t &c : trim_cubes) {
						objs.emplace_back(c, TYPE_WALL_TRIM, 0, dim, 0, ext_flags, 1.0, SHAPE_CUBE, trim_color); // floor trim
						if (!has_ceil_trim || is_house) continue;
						set_cube_zvals(c, ceil_trim_z1, ceil_trim_z2); // okay to edit in-place here
						objs.emplace_back(c, TYPE_WALL_TRIM, 0, dim, !dir, flags, 1.0, SHAPE_ANGLED, trim_color); // ceiling trim
					}
				} // for f
			} // for dir
		} // for dim
	} // for i
	if (!is_cube() || is_rotated()) return; // window trim is not yet working for non-cube and rotated buildings
	add_window_trim_and_coverings(1, 0); // add_trim=1, add_coverings=0
}

void building_t::add_window_trim_and_coverings(bool add_trim, bool add_coverings, bool add_ext_sills) {
	assert(add_trim || add_coverings || add_ext_sills);
	if (!has_windows()) return; // no windows
	float const trim_thickness(get_trim_thickness()), ext_wall_toler(0.01*trim_thickness); // required to prevent z-fighting when AA is disabled
	float const window_h_border(WINDOW_BORDER_MULT*get_window_h_border()), window_v_border(WINDOW_BORDER_MULT*get_window_v_border()); // (0, 1) range
	// Note: depth must be small to avoid object intersections; this applies to the windowsill as well
	float const window_trim_width(0.75*get_wall_thickness()), window_trim_depth(1.0*trim_thickness), windowsill_depth(1.0*trim_thickness);
	float const floor_spacing(get_window_vspace()), window_offset(0.01*floor_spacing); // must match building_draw_t::add_section()
	colorRGBA const &trim_color(is_house ? WHITE : DK_GRAY);
	colorRGBA const &frame_color(get_material().window_color);
	vect_room_object_t &trim_objs(interior->room_geom->trim_objs), &objs(interior->room_geom->objs);
	static vect_vnctcc_t wall_quad_verts;
	wall_quad_verts.clear();
	get_all_drawn_window_verts_as_quads(wall_quad_verts);
	rand_gen_t rgen;
	rgen.set_state(wall_quad_verts.size(), interior->rooms.size());
	rgen.rand_mix();
	
	if (add_ext_sills && !rgen.rand_bool()) { // 50% of the time
		add_ext_sills = 0;
		if (!add_trim && !add_coverings) return; // nothing else to add
	}
	for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
		cube_t c;
		float tx1, tx2, tz1, tz2;
		if (!get_wall_quad_window_area(wall_quad_verts, i, c, tx1, tx2, tz1, tz2)) continue;
		bool const dim(c.dy() < c.dx()), dir(wall_quad_verts[i].get_norm()[dim] > 0.0);
		assert(c.get_sz_dim(dim) == 0.0); // must be zero size in one dim (X or Y oriented); could also use the vertex normal
		float const d_tx_inv(1.0f/(tx2 - tx1)), d_tz_inv(1.0f/(tz2 - tz1));
		float const window_width(c.get_sz_dim(!dim)*d_tx_inv), window_height(c.dz()*d_tz_inv); // window_height should be equal to window_vspacing
		float const border_xy(window_width*window_h_border), border_z(window_height*window_v_border), dscale(dir ? -1.0 : 1.0);
		cube_t window(c); // copy dim <dim>
		window.translate_dim(dim, dscale*window_offset);
		window.d[dim][!dir] += dscale*window_trim_depth; // add thickness on interior of building
		window.d[dim][ dir] += dscale*ext_wall_toler; // slight bias away from the exterior wall
		unsigned const ext_flags(RO_FLAG_NOCOLL | (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));

		for (float z = tz1; z < tz2; z += 1.0) { // each floor
			float const bot_edge(c.z1() + (z - tz1)*window_height);
			set_cube_zvals(window, bot_edge+border_z, bot_edge+window_height-border_z);
			bool const add_separators(is_house && rgen.rand_bool()); // 50% of walls/floors for houses; not for glass block windows?
			bool const one_dim_separators(add_separators && rgen.rand_bool()); // 1=vert/horiz separators only, 0=cross shaped separators
			bool const sep_dim((window_height - 2.0*border_z) < 0.9*(window_width - 2.0*border_xy)); // split in larger-ish dim

			for (float xy = tx1; xy < tx2; xy += 1.0) { // windows along each wall
				float const low_edge(c.d[!dim][0] + (xy - tx1)*window_width);
				window.d[!dim][0] = low_edge + border_xy;
				window.d[!dim][1] = low_edge + window_width - border_xy;
				if (add_coverings) {add_window_coverings(window, dim, dir);}

				if (add_ext_sills) {
					cube_t sill(window);
					sill.z1() -= 0.04*window_height;
					sill.z2()  = sill.z1() + 0.035*window_height;
					sill.expand_in_dim(!dim, 0.050*window_height);
					sill.d[dim][!dir] -= dscale*window_offset; // flush with exterior wall to avoid clipping through interior
					sill.d[dim][ dir] -= dscale*0.06*window_height; // extend out from the wall
					objs.emplace_back(sill, TYPE_WIND_SILL, 0, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR), 1.0, SHAPE_CUBE, frame_color);
				}
				if (!add_trim) continue;
				// add window trim
				float const window_ar(window.get_sz_dim(!dim)/window.dz());
				float const side_trim_width(window_trim_width*((window_ar > 1.5) ? (window_ar - 0.5) : 1.0)); // widen for very wide windows to cover holes at stretched edges
				cube_t top(window), bot(window), side(window);
				top.z1()  = window.z2();
				top.z2() += window_trim_width;
				bot.z2()  = window.z1();
				bot.z1() -= window_trim_width;
				bot.d[dim][!dir] += dscale*(windowsill_depth - window_trim_depth); // shift out further for windowsill
				top.expand_in_dim(!dim, side_trim_width);
				bot.expand_in_dim(!dim, side_trim_width);
				// interior trim
				unsigned const trim_start(trim_objs.size());
				trim_objs.emplace_back(top, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL, trim_color);
				trim_objs.emplace_back(bot, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL, trim_color);

				for (unsigned s = 0; s < 2; ++s) { // left/right sides
					side.d[!dim][ s] = window.d[!dim][s] - (s ? -1.0 : 1.0)*side_trim_width;
					side.d[!dim][!s] = window.d[!dim][s];
					trim_objs.emplace_back(side, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL, trim_color);
				}
				unsigned const trim_end(trim_objs.size());
				// exterior trim
				for (unsigned ix = trim_start; ix < trim_end; ++ix) {
					room_object_t trim(trim_objs[ix]); // start with interior trim
					trim.dir  ^= 1; // flip to exterior
					trim.flags = RO_FLAG_EXTERIOR | RO_FLAG_HAS_EXTRA | RO_FLAG_NOCOLL | (dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI); // swap lo/hi flags; flag as fully exterior
					trim.color = frame_color;
					trim.d[dim][!dir]  = trim.d[dim][dir]; // shift to other side of the wall
					trim.d[dim][ dir] -= dscale*2.5*window_trim_depth; // shift away from the wall
					trim_objs.push_back(trim);
				} // for i
				// add cross shaped window pane separators
				if (add_separators) {
					bool is_split(0), has_bathroom_block_window(0);
					int const room_id(get_room_id_for_window(window, dim, dir, is_split));
					unsigned floor_ix(0); // unused

					if (get_room_type_and_floor(room_id, window.zc(), floor_ix) == RTYPE_BATH) { // check for bathroom block windows
						cube_t window_exp(window);
						window_exp.expand_by_xy(2.0*trim_thickness); // window is shifted by trim_thickness, so need to expand by more than that

						for (room_object_t const &obj : objs) {
							if (obj.type == TYPE_WINDOW && obj.intersects(window_exp)) {has_bathroom_block_window = 1; break;}
						}
					}
					if (!has_bathroom_block_window) {
						unsigned const base_flags(RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR);
						float const sep_hwidth(0.2*side_trim_width*(one_dim_separators ? 2.0 : 1.0));

						if (!one_dim_separators || sep_dim == 1) { // vertical separator
							cube_t sep(window);
							set_wall_width(sep, window.get_center_dim(!dim), sep_hwidth, !dim);
							trim_objs.emplace_back(sep, TYPE_WALL_TRIM, 0, dim, dir, (base_flags | RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP), 1.0, SHAPE_TALL, trim_color);
						}
						if (!one_dim_separators || sep_dim == 0) { // horizontal separator
							cube_t sep(window);
							set_wall_width(sep, window.zc(), sep_hwidth, 2);
							trim_objs.emplace_back(sep, TYPE_WALL_TRIM, 0, dim, dir, base_flags, 1.0, SHAPE_SHORT, trim_color); // skip !dim ends
						}
					}
				} // end add_separators
			} // for xy
		} // for z
	} // for i
	if (is_house && add_trim && rgen.rand_float() < 0.2) { // add exterior house wall first floor trim
		float const zval(ground_floor_z1 + floor_spacing), width(4.0*window_trim_depth);
		vect_cube_t trims;

		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if ((p->z2() - ground_floor_z1) < 1.5*floor_spacing) continue; // single story, skip

			for (unsigned dim = 0; dim < 2; ++dim) {
				for (unsigned dir = 0; dir < 2; ++dir) {
					unsigned const flags(RO_FLAG_NOCOLL | (dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI) | RO_FLAG_EXTERIOR | RO_FLAG_HAS_EXTRA);
					cube_t trim(*p);
					set_cube_zvals(trim, zval, zval+width);
					trim.d[dim][!dir]  = trim.d[dim][dir]; // set to zero area
					trim.d[dim][ dir] += (dir ? 1.0 : -1.0)*width; // set width adjacent to wall
					trim.expand_in_dim(!dim, width); // overlap the ends so that corners join properly (with some overlap)
					trims.clear();
					trims.push_back(trim);

					for (auto p2 = parts.begin(); p2 != get_real_parts_end(); ++p2) {
						if (p2 == p) continue; // skip self
						cube_t sub_cube(*p2);
						sub_cube.z2() += 0.1*floor_spacing; // clip if part ends just below the trim because the roof will clip through it
						subtract_cube_from_cubes(sub_cube, trims);
					}
					for (cube_t const &trim : trims) {
						trim_objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, flags, 1.0, SHAPE_TALL, trim_color); // not tall, but we need to draw the bottom surface
					}
				} // for dir
			} // for dim
		} // for p
	}
}

// *** Windows ***

room_type building_t::get_room_type_and_floor(int room_id, float zval, unsigned &floor_ix) const {
	if (room_id < 0) return RTYPE_NOTSET; // error?
	room_t const &room(get_room(room_id));
	floor_ix = room.get_floor_containing_zval(zval, get_window_vspace());
	return room.get_room_type(floor_ix);
}
void building_t::add_window_coverings(cube_t const &window, bool dim, bool dir) {
	// add blinds to some windows based on the containing room type for this floor
	bool is_split(0);
	int const room_id(get_room_id_for_window(window, dim, dir, is_split));
	if (is_split) return; // window split across multiple rooms - how do we handle this? for now skip it
	unsigned floor_ix(0);

	switch (get_room_type_and_floor(room_id, window.zc(), floor_ix)) {
	case RTYPE_BED: case RTYPE_MASTER_BED: add_window_blinds(window, dim, dir, room_id, floor_ix); break; // bedroom
	case RTYPE_BATH: add_bathroom_window(window, dim, dir, room_id, floor_ix); break; // bathroom
	} // end switch
}

void building_t::add_window_blinds(cube_t const &window, bool dim, bool dir, unsigned room_id, unsigned floor) {
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), extend(0.9*wall_thickness); // extend is 15% larger than window trim width
	vect_room_object_t &objs(interior->room_geom->objs);
	bool vertical((mat_ix + interior->rooms.size() + parts.size()) & 1); // something that's per-building
	
	if (vertical) { // check for horizontal wall clearance
		float const min_space(2.0*wall_thickness);
		room_t const &room(get_room(room_id));
		if ((window.d[!dim][0] - min_space) < room.d[!dim][0] || (window.d[!dim][1] + min_space) > room.d[!dim][1]) {vertical = 0;} // not enough space for vertical
		else if (is_house) { // check for blocking bedroom closets
			cube_t window_exp(window);
			window_exp.expand_by_xy(min_space);

			for (auto i = objs.begin(); i != objs.end(); ++i) { // iterate over room objects placed so far
				if (i->type == TYPE_CLOSET && i->intersects(window_exp)) {vertical = 0; break;}
			}
		}
	}
	rand_gen_t rgen;
	rgen.set_state((123*room_id + 211*interior->rooms.size()), (777*floor + 1));
	// open_amt is a mix of 50% room-based and 50% window-based to get somewhat consistent levels per room
	bool const full_open(rgen.rand_float() < 0.75);
	float const open_amt(0.9*(full_open ? 1.0 : 0.5*(rgen.rand_float() + fract(1123.7*objs.size())))); // 0-90% 25% the time, 90% for the rest
	float const thickness(0.15*wall_thickness*(vertical ? 0.05 : (open_amt + 0.025))); // vertical blinds have no furniture clearance and can't bunch up
	float window_gap(0.02*wall_thickness);
	cube_t c(window);
	float &window_edge(c.d[dim][dir]);

	// add a slight gap for the wall trim; increase until it's different from the wall to avoid z-fighting when very far from the origin
	while (1) {
		window_edge += (dir ? -1.0 : 1.0)*window_gap;
		if (window_edge != window.d[dim][dir]) break;
		window_gap *= 2.0;
	}
	c.d[dim][!dir] += (dir ? -1.0 : 1.0)*thickness;
	expand_to_nonzero_area(c, thickness, dim);
	c.expand_in_dim(2, extend); // extend in Z to cover window trim

	if (vertical) { // vertical, moves horizontally
		c.expand_in_dim(!dim, 1.5*wall_thickness); // larger expand value (beyond the wall trim)
		float const center(c.get_center_dim(!dim)), half_width(0.5*c.get_sz_dim(!dim));
		float const shift_val(1.44*max(0.0f, (open_amt - 0.4f))*half_width); // more likely to be fully closed

		for (unsigned d = 0; d < 2; ++d) { // left, right
			cube_t c2(c);
			c2.d[!dim][!d] = center - (d ? -1.0 : 1.0)*shift_val;
			objs.emplace_back(c2, TYPE_BLINDS, room_id, dim, dir, (RO_FLAG_NOCOLL | (d ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI)), 1.0); // always fully lit
		}
	}
	else { // horizontal, moves vertically
		c.expand_in_dim(!dim, extend); // expand width to cover trim +15% WT
		c.z2() += extend + 0.05*floor_spacing; // expand height to allow space for it to bunch up at the top
		c.z1() += open_amt*window.dz(); // raise amount is random per-room
		objs.emplace_back(c, TYPE_BLINDS, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), 1.0); // always fully lit
	}
}

void building_t::add_bathroom_window(cube_t const &window, bool dim, bool dir, unsigned room_id, unsigned floor) { // frosted window blocks, for houses or office buildings
	room_t const &room(get_room(room_id));
	if (count_ext_walls_for_room(room, window.z1()) != 1) return; // looks odd to have window block walls at the corner of a building, so only enable this for single exterior walls
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c(window);
	c.translate_dim(dim, (dir ? 1.0 : -1.0)*get_trim_thickness());
	unsigned const flags(RO_FLAG_NOCOLL | (room.is_lit_on_floor(floor) ? RO_FLAG_LIT : 0));
	objs.emplace_back(c, TYPE_WINDOW, room_id, dim, dir, flags, 1.0, SHAPE_CUBE, WHITE); // always lit
}

int building_t::get_room_id_for_window(cube_t const &window, bool dim, bool dir, bool &is_split) const {
	assert(interior);
	float const wall_thickness(get_wall_thickness());
	point const center(window.get_cube_center());

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (center[!dim] < r->d[!dim][0] || center[!dim] >= r->d[!dim][1]) continue; // test center point for windows that straddle two rooms
		if (center.z < r->z1() || center.z > r->z2()) continue;
		if (fabs(center[dim] - r->d[dim][dir]) > wall_thickness) continue; // wrong wall
		is_split = (window.d[!dim][0] < r->d[!dim][0] || window.d[!dim][1] > r->d[!dim][1]);
		return (r - interior->rooms.begin()); // found
	} // for r
	return -1; // not found
}

// *** Stairs and Elevators ***

void add_elevator_button(point const &pos, float button_radius, bool dim, bool dir, unsigned elevator_id, unsigned floor_id, bool inside, bool is_up, vect_room_object_t &objs) {
	float const button_thickness(0.25*button_radius);
	cube_t c; c.set_from_point(pos);
	c.expand_in_dim(!dim, button_radius);
	c.expand_in_dim(2, button_radius); // Z
	c.d[dim][dir] += (dir ? 1.0 : -1.0)*button_thickness;
	unsigned flags(RO_FLAG_NOCOLL);
	if (inside) {flags |= RO_FLAG_IN_ELEV;}
	else        {flags |= (is_up ? RO_FLAG_ADJ_TOP : RO_FLAG_ADJ_BOT);}
	expand_to_nonzero_area(c, button_thickness, dim);
	objs.emplace_back(c, TYPE_BUTTON, elevator_id, dim, dir, flags, 1.0, SHAPE_CYLIN, colorRGBA(1.0, 0.9, 0.5)); // room_id=elevator_id
	objs.back().obj_id = floor_id; // encode floor index as obj_id
}
void add_floor_number(unsigned floor_ix, unsigned floor_offset, bool has_parking_garage, ostringstream &oss) { // Note: floor_ix=1 is ground floor
	oss.str("");
	int const adj_floor_ix(int(floor_ix) - int(floor_offset));
	if (adj_floor_ix <= 0) {oss << (has_parking_garage ? "P" : "B") << (1 - adj_floor_ix);} // basement floors
	else {oss << adj_floor_ix;} // above ground floors
}
void set_floor_text_for_sign(room_object_t &sign, unsigned floor_ix, unsigned floor_offset, bool has_parking_garage, ostringstream &oss) { // Note: floor_ix=1 is ground floor
	add_floor_number(floor_ix, floor_offset, has_parking_garage, oss);
	sign.obj_id = register_sign_text(oss.str());
	int const adj_floor_ix(int(floor_ix) - int(floor_offset));
	float width_adj(0.0);
	if      (adj_floor_ix <= 0 ) {width_adj =  ((adj_floor_ix <  0) ? 0.1 : 0.0);} // basement floors; widen if lower than B1
	else if (adj_floor_ix <  10) {width_adj = -((adj_floor_ix == 1) ? 0.2 : 0.1);} // 1-10: make narrow
	else if (adj_floor_ix >= 20) {width_adj =  0.1;} // 20+: widen
	if (width_adj != 0.0) {sign.expand_in_dim(!sign.dim, width_adj*sign.get_sz_dim(!sign.dim));}
}
unsigned building_t::calc_floor_offset(float zval) const { // for basements
	return ((zval < ground_floor_z1) ? round_fp((ground_floor_z1 - zval)/get_window_vspace()) : 0);
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), half_thick(0.5*floor_thickness);
	float const wall_thickness(get_wall_thickness()), elevator_car_z1_add(0.05*floor_thickness), fc_thick_scale(get_elevator_fc_thick_scale());
	float const stairs_sign_width(1.0*wall_thickness);
	vect_room_object_t &objs(interior->room_geom->objs);
	ostringstream oss; // reused across elevators/floors

	// add floor signs for U-shaped stairs
	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator || i->for_ramp || i->shape != SHAPE_U) continue; // not U-shaped stairs
		// stacked conn stairs start at floor 0 but are really the top floor of the part below; i->floor is not a global index and can't be used
		unsigned const floor_offset(calc_floor_offset(bcube.z1())); // use building z1 - should return number of underground levels
		unsigned const real_floor(round_fp((i->z1() - bcube.z1())/get_window_vspace()));
		unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING);
		point center;
		center[ i->dim] = i->d[i->dim][!i->dir]; // front of stairs
		center[!i->dim] = i->get_center_dim(!i->dim);
		if (has_parking_garage && i->z1() < ground_floor_z1) {center[!i->dim] += 0.25*i->get_sz_dim(!i->dim);} // shift to the side for parking garages to avoid center beams
		center.z = i->z1();
		cube_t sign(center);
		sign.d[i->dim][!i->dir] += (i->dir ? -1.0 : 1.0)*0.25*wall_thickness; // set sign thickness
		sign.expand_in_dim(!i->dim, stairs_sign_width); // set sign width
		sign.z1() -= 2.5*wall_thickness; // set sign height
		objs.emplace_back(sign, TYPE_SIGN, 0, i->dim, !i->dir, flags, 1.0, SHAPE_CUBE, DK_BLUE); // no room_id
		set_floor_text_for_sign(objs.back(), real_floor, floor_offset, has_parking_garage, oss);

		// if this is the top landing, we need to add a floor sign on the ceiling above it for the top floor
		if (i->is_at_top && !i->roof_access) {
			sign.translate_dim(2, window_vspacing); // move up one floor

			if (check_skylight_intersection(sign)) { // check for sklight and move sign to the side
				sign.translate_dim(!i->dim, 0.5*(i->get_sz_dim(!i->dim) - stairs_sign_width)); // translate to the positive side
				flags |= RO_FLAG_ADJ_TOP; // dra the top surface
			}
			objs.emplace_back(sign, TYPE_SIGN, 0, i->dim, !i->dir, flags, 1.0, SHAPE_CUBE, DK_BLUE); // no room_id
			set_floor_text_for_sign(objs.back(), (real_floor + 1), floor_offset, has_parking_garage, oss);
		}
	} // for i

	// add elevator lights, and signs on each floor; must be done before setting buttons_start
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		// add light
		i->light_obj_id = objs.size();
		float const light_zval(max(ground_floor_z1, i->z1()) + elevator_car_z1_add + (1.0 - fc_thick_scale)*window_vspacing); // starts on the ground floor
		cube_t light(point(i->xc(), i->yc(), light_zval));
		light.z1() -= 0.02*window_vspacing;
		light.expand_by_xy(0.06*window_vspacing);
		objs.emplace_back(light, TYPE_LIGHT, i->room_id, i->dim, i->dir, (RO_FLAG_NOCOLL | RO_FLAG_IN_ELEV | RO_FLAG_LIT), 0.0, SHAPE_CYLIN, WHITE);
		objs.back().obj_id = uint16_t(i - interior->elevators.begin()); // encode elevator index as obj_id
		// add floor signs
		unsigned const num_floors(calc_num_floors(*i, window_vspacing, floor_thickness));
		unsigned const floor_offset(calc_floor_offset(i->z1()));
		float const ewidth(i->get_width());
		cube_t sign;
		sign.d[i->dim][0] = sign.d[i->dim][1] = i->d[i->dim][i->dir];
		sign.d[i->dim][i->dir] += (i->dir ? 1.0 : -1.0)*0.1*wall_thickness; // front of sign
		set_wall_width(sign, (i->d[!i->dim][1] - 0.1*ewidth), 0.04*ewidth, !i->dim); // to the high side, opposite the call button

		for (unsigned f = 0; f < num_floors; ++f) { // Note: floor number starts at 1 even if the elevator doesn't extend to the ground floor
			sign.z1() = i->z1()   + (f + 0.5)*window_vspacing;
			sign.z2() = sign.z1() + 0.1*ewidth;
			objs.emplace_back(sign, TYPE_SIGN, i->room_id, i->dim, i->dir, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, DK_BLUE);
			set_floor_text_for_sign(objs.back(), f+1, floor_offset, has_parking_garage, oss);
		}
		// add concrete flooring at the base of the elevator, over the carpet
		cube_t efloor(*i);
		efloor.expand_by_xy(-0.25*half_thick); // less than the width of the elevator walls, to prevent z-fighting
		efloor.z1() += half_thick; // on top of the carpet
		efloor.z2()  = efloor.z1() + 0.01*half_thick; // set thickness (very thin)
		objs.emplace_back(efloor, TYPE_FLOORING, i->room_id, i->dim, i->dir, (RO_FLAG_NOCOLL | RO_FLAG_IN_ELEV), 1.0, SHAPE_CUBE, WHITE, FLOORING_CONCRETE);
	} // for e
	interior->room_geom->buttons_start = objs.size();

	// add elevator buttons for each floor; must be done before setting stairs_start
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		float const button_radius(0.3*wall_thickness), ewidth(i->get_width());
		unsigned const num_floors(calc_num_floors(*i, window_vspacing, floor_thickness)), elevator_id(i - interior->elevators.begin());
		assert(num_floors > 1); // otherwise, why have an elevator?
		i->button_id_start = objs.size();

		// call buttons on each floor outside the elevator
		for (unsigned f = 0; f < num_floors; ++f) {
			point pos;
			pos[ i->dim] = i->d[i->dim][i->dir]; // front of the elevator
			pos[!i->dim] = i->d[!i->dim][0] + 0.1*ewidth; // to the low side

			for (unsigned d = 0; d < 2; ++d) { // {down, up} call buttons
				if ((d == 0 && f == 0) || (d == 1 && f == num_floors-1)) continue; // no floor above/below
				pos.z = i->z1() + (f + 0.05*d + 0.45)*window_vspacing;
				add_elevator_button(pos, button_radius, i->dim, i->dir, elevator_id, f, 0, d, objs); // inside=0, is_up=d
			}
		} // for f
		// call buttons for each floor inside the elevator car; first find the panel location for the starting elevator car position
		cube_t elevator_car(*i);
		max_eq(elevator_car.z1(), ground_floor_z1); // always starts on the ground floor, not the bottom of the basement
		elevator_car.z1() += elevator_car_z1_add;
		elevator_car.z2()  = elevator_car.z1() + window_vspacing; // currently at the bottom floor
		cube_t const panel(get_elevator_car_panel(room_object_t(elevator_car, TYPE_ELEVATOR, elevator_id, i->dim, i->dir, 0), fc_thick_scale));
		float const dz(panel.dz()), button_spacing(dz/(num_floors + 1)); // add extra spacing on bottom and top of panel
		float const inner_button_radius(min(button_radius, min(0.35f*button_spacing, 0.25f*panel.get_sz_dim(!i->dim)))); // may need to be smaller
		point pos;
		pos[ i->dim] = panel.d[i->dim][!i->dir]; // front face of inside panel
		pos[!i->dim] = panel.get_center_dim(!i->dim) + 0.8*inner_button_radius; // a bit right of center to make room for floor number text
		
		for (unsigned f = 0; f < num_floors; ++f) {
			pos.z = panel.z1() + (f + 1)*button_spacing;
			add_elevator_button(pos, inner_button_radius, i->dim, !i->dir, elevator_id, f, 1, 0, objs); // inside=1, is_up=0, pointing in opposite dir
		}
		i->button_id_end = objs.size();
	} // for e
	interior->room_geom->stairs_start  = objs.size();
	interior->room_geom->has_elevators = (!interior->elevators.empty());
	colorRGBA const railing_colors[3]  = {GOLD, LT_GRAY, BLACK};
	colorRGBA const railing_color(railing_colors[rgen.rand()%3]); // set per-building

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator || i->for_ramp) continue; // for elevator or ramp, not stairs
		unsigned const num_stairs(i->get_num_stairs());
		float const stair_dz(window_vspacing/(num_stairs+1)), stair_height(stair_dz + floor_thickness);
		bool const dim(i->dim), dir(i->dir), has_side_walls(i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES || i->shape == SHAPE_U);
		bool const has_wall_both_sides(i->against_wall[0] && i->against_wall[1]); // ext basement stairs
		bool const side(dir); // for U-shaped stairs; for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
		// Note: stairs always start at floor_thickness above the landing z1, ignoring landing z2/height
		float const tot_len(i->get_sz_dim(dim)), floor_z(i->z1() + floor_thickness - window_vspacing), step_len_pos(tot_len/num_stairs);
		float const wall_hw(min(STAIRS_WALL_WIDTH_MULT*max(step_len_pos, stair_dz), 0.25f*stair_dz));
		float const stairs_zmin(i->in_ext_basement ? interior->basement_ext_bcube.z1() : bcube.z1());
		float step_len((dir ? 1.0 : -1.0)*step_len_pos), z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i);

		if (i->shape != SHAPE_U) { // straight stairs
			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				// tag basement bottom stair so that no width extension is added, since this may clip through the basement door
				unsigned const flags((n == 0 && !i->in_ext_basement && i->z1() < ground_floor_z1) ? RO_FLAG_RSTAIRS : 0);
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				set_cube_zvals(stair, max(stairs_zmin, (z + 0.4f*stair_height)), z+stair_height); // don't go below the floor (Note: z1 was (z + 0.5f*half_thick))
				assert(stair.z1() < stair.z2());
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir, flags); // Note: room_id=0, not tracked, unused
			} // for n
		}
		else { // U-shaped stairs
			float const mid(i->get_center_dim(!dim));
			stair.d[!dim][side] = mid;
			step_len *= 2.0;

			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				if (n == num_stairs/2) { // reverse direction and switch to other side
					step_len *= -1.0;
					stair.d[!dim][ side] = i->d[!dim][side];
					stair.d[!dim][!side] = mid;
				}
				assert(!(num_stairs & 1)); // require num_stairs to be an even number
				bool const is_rev(n >= num_stairs/2), stairs_dir(dir^is_rev);
				stair.d[dim][!stairs_dir] = pos; stair.d[dim][stairs_dir] = pos + step_len;
				set_cube_zvals(stair, max(floor_z, z), z+stair_height);
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, side^is_rev); // Note: room_id=0, not tracked, unused
				objs.back().shape = SHAPE_STAIRS_U;
			} // for n
		}
		// add walls and railings
		bool const extend_walls_up(i->is_at_top && !i->roof_access); // space above is open, add a wall so that people can't fall down the stairs
		float const railing_z2(i->z2() + (i->roof_access ? 0.025*i->dz() : 0.0)); // capture z2 before we change it; move roof access railing up a bit to offset the shrink resize
		float const wall_bottom(floor_z - half_thick), railing_side_dz(0.5*stair_dz); // for U-shaped stairs
		cube_t wall(*i);
		if (extend_walls_up) {wall.z2() += window_vspacing - floor_thickness;}
		else {wall.z2() -= 0.5*floor_thickness;} // prevent z-fighting on top floor
		wall.z1() = max((stairs_zmin + half_thick), wall_bottom); // full height
		set_wall_width(wall, i->d[dim][dir], wall_hw, dim);

		if ((i->shape == SHAPE_WALLED && !(i->against_wall[0] || i->against_wall[1]) && (!i->stack_conn || !i->is_at_top)) || i->shape == SHAPE_U) {
			objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir); // add wall at back/end of stairs
		}
		else if ((i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES) && extend_walls_up) { // add upper section only
			cube_t wall_upper(wall);
			set_wall_width(wall_upper, (i->d[dim][!dir] + (dir ? 1.0 : -1.0)*wall_hw), wall_hw, dim); // move to the other side
			wall_upper.z1() = railing_z2;

			for (unsigned d = 0; d < 2; ++d) {
				// if there's no wall, extend to cover the gap where the wall would be; slightly smaller to avoid z-fighting on exterior walls (happens to be the trim thickness);
				// we can't just skip this face for exterior walls because it may be visible through a window
				if (i->against_wall[d]) {wall_upper.d[!dim][d] += (d ? 1.0 : -1.0)*0.9*wall_thickness;}
			}
			objs.emplace_back(wall_upper, TYPE_STAIR_WALL, 0, dim, dir); // add wall at back/end of stairs
		}
		wall.d[dim][!dir] = i->d[dim][!dir];

		for (unsigned d = 0; d < 2; ++d) { // sides of stairs
			set_wall_width(wall, i->d[!dim][d], wall_hw, !dim);
			wall.expand_in_dim(dim, 0.01*wall_hw); // just enough to avoid z-fighting with stairs
			bool const add_wall(has_side_walls && !i->against_wall[d]); // don't add a wall if the stairs are already against a wall
			if (add_wall) {objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir);} // add walls around stairs for this floor

			if (i->has_railing) { // add railings
				bool railing_dir(dir);
				cube_t railing(wall);
				unsigned flags(add_wall ? RO_FLAG_NOCOLL : 0);
				if (!has_side_walls && !has_wall_both_sides/*!i->against_wall[d]*/) {flags |= RO_FLAG_OPEN;} // use this flag to indicate no walls, need balusters
				railing.z2() = railing_z2;

				if (add_wall || i->roof_access) {
					railing.translate_dim(!dim, (d ? -1.0 : 1.0)*2.0*wall_hw); // shift railing inside of walls
					railing.expand_in_dim( dim, -(i->roof_access ? 2.0 : 1.0)*wall_hw); // shrink slightly to avoid clipping through an end wall
				}
				if (i->shape == SHAPE_U) { // adjust railing height/angle to match stairs
					float const z_split(railing.zc());
					if (bool(d) == side) {railing.z1() = z_split + railing_side_dz; flags |= RO_FLAG_ADJ_HI; railing_dir ^= 1;}
					else                 {railing.z2() = z_split - railing_side_dz; flags |= RO_FLAG_ADJ_LO;}
				}
				objs.emplace_back(railing, TYPE_RAILING, 0, dim, railing_dir, flags, 1.0, SHAPE_CUBE, railing_color);
			}
		} // for d
		if (i->has_railing && i->shape == SHAPE_U) { // add a railing for the back wall of U-shaped stairs
			float const railing_zc(wall_bottom + 0.819*window_vspacing); // determined experimentally
			cube_t railing(*i);
			set_wall_width(railing, (i->d[dim][dir] + (dir ? -1.0 : 1.0)*2.0*wall_hw), wall_hw, dim);
			set_wall_width(railing, railing_zc, 1.4*railing_side_dz, 2); // set zvals
			objs.emplace_back(railing, TYPE_RAILING, 0, !dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_HI | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_BOT), 1.0, SHAPE_CUBE, railing_color); // no ends
		}
		else if (i->has_railing && !has_wall_both_sides && (i->stack_conn || (extend_walls_up && i->shape == SHAPE_STRAIGHT))) {
			// add railings around the top if: straight + top floor with no roof access, connector stairs, or basement stairs
			room_object_t railing(*i, TYPE_RAILING, 0, !dim, dir, (RO_FLAG_TOS | RO_FLAG_ADJ_BOT), 1.0, SHAPE_CUBE, railing_color); // flag to skip drawing ends
			railing.z1()  = railing.z2(); // starts at the floor
			railing.z2() += window_vspacing - floor_thickness;
			set_wall_width(railing, (i->d[dim][!dir] + (dir ? -1.0 : 1.0)*wall_hw), wall_hw, dim); // no overlap with stairs cutout

			for (unsigned d = 0; d < 2; ++d) {
				if (has_side_walls && !i->against_wall[d]) {railing.d[!dim][d] += (d ? -1.0 : 1.0)*2.0*wall_hw;} // shift railing inside of walls
			}
			// if the stairs extend to a wall, the back railing can be omitted; this is rare but happens in one house near the starting area
			cube_t railing_center(railing);
			railing_center.d[dim][dir] = railing_center.d[dim][!dir]; // shrink to zero area in this dim
			bool in_wall(0);

			for (auto const &w : interior->walls[dim]) {
				if (w.contains_cube(railing_center)) {in_wall = 1;}
			}
			if (!in_wall) {objs.emplace_back(railing);} // back railing
			railing.d[dim][dir] = i->d[dim][dir]; // extend to the front of the stairs
			railing.dim  ^= 1;
			railing.flags = RO_FLAG_TOS | RO_FLAG_ADJ_TOP; // flag so that no vertical pole is added

			for (unsigned d = 0; d < 2; ++d) { // sides of stairs
				railing.dir = bool(d);
				set_wall_width(railing, i->d[!dim][d], wall_hw, !dim);
				if (has_side_walls && !i->against_wall[d]) {railing.translate_dim(!dim, (d ? -1.0 : 1.0)*2.0*wall_hw);} // shift railing inside of walls
				objs.emplace_back(railing);
			}
		}
	} // for i (landings)
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		unsigned const elevator_id(i - interior->elevators.begin()); // used for room_object_t::room_id
		cube_t elevator_car(*i);
		max_eq(elevator_car.z1(), ground_floor_z1); // always starts on the ground floor, not the bottom of the basement
		elevator_car.z1() += elevator_car_z1_add; // to prevent z-fighting when looking at the building from the bottom
		elevator_car.z2()  = elevator_car.z1() + window_vspacing; // one floor of height
		i->car_obj_id = objs.size();
		objs.emplace_back(elevator_car, TYPE_ELEVATOR, elevator_id, i->dim, i->dir, RO_FLAG_DYNAMIC);
		objs.back().drawer_flags = (uint16_t)calc_num_floors(*i, window_vspacing, floor_thickness); // store the number of floors in drawer_flags; used for drawing
		objs.back().item_flags   = (uint16_t)calc_floor_offset(i->z1()); // use correct starting floor index
	} // for i
}

int building_t::get_ext_door_dir(cube_t const &door_bcube, bool dim) const { // erturn value of 2 means 'not found'
	float const width(door_bcube.get_sz_dim(!dim));

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find part containing this door so that we can get the correct dir
		if (is_basement(p)) continue; // skip the basement
		if (p->z1() != ground_floor_z1) continue; // not ground floor
		if (p->d[!dim][1] < door_bcube.d[!dim][1] || p->d[!dim][0] > door_bcube.d[!dim][0]) {continue;} // not contained in this dim
		if      (fabs(p->d[dim][0] - door_bcube.d[dim][0]) < 0.1*width) return 0;
		else if (fabs(p->d[dim][1] - door_bcube.d[dim][1]) < 0.1*width) return 1;
	} // for p
	if (is_cube()) { // some non-cube buildings have no exterior doors
		cout << "Warning: Failed to find building exterior door: " << TXT(bcube.str()) << TXT(door_bcube.str()) << TXT(is_house) << endl; // debug printout
		//assert(0); // never gets here (too strong?)
	}
	return 2; // not found
}

void building_t::add_doorbell_and_lamp(tquad_with_ix_t const &door, rand_gen_t &rgen) { // and porch packages
	cube_t const door_bcube(door.get_bcube());
	bool const dim(door_bcube.dy() < door_bcube.dx());
	int const dir_ret(get_ext_door_dir(door_bcube, dim));
	if (dir_ret > 1) return; // not found, skip doorbell and lamp
	bool dir(dir_ret != 0);
	bool const side(dir ^ dim); // currently always to the right, which matches the door handle side
	float const door_width(door_bcube.get_sz_dim(!dim)), half_width(0.016*door_width), half_height(1.8*half_width), button_thickness(0.1*half_width);
	float const zval(door_bcube.z1() + 0.55*door_bcube.dz());
	float const pos(door_bcube.d[!dim][side] + (side ? 1.0 : -1.0)*5.0*half_width);
	cube_t c;
	c.d[dim][0  ]  = c.d[dim][1] = door_bcube.d[dim][dir] - 0.02*(dir ? 1.0 : -1.0)*get_window_vspace(); // slightly in front of exterior wall
	c.d[dim][dir] += (dir ? 1.0 : -1.0)*button_thickness;
	set_cube_zvals(c, (zval - half_height), (zval + half_height));
	set_wall_width(c, pos, half_width, !dim);
	expand_to_nonzero_area(c, button_thickness, dim);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const room_id(0); // not valid, should be unused
	float const tot_light_amt(1.0);
	unsigned const base_flags(RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR);
	objs.emplace_back(c, TYPE_BUTTON, room_id, dim, dir, (base_flags | RO_FLAG_LIT), tot_light_amt, SHAPE_CYLIN); // always lit

	// add a wall lamp above the button if there's a porch, garage, or shed (L-shaped house)
	if ((has_porch() || has_garage || has_shed) && building_obj_model_loader.is_model_valid(OBJ_MODEL_WALL_LAMP)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_WALL_LAMP)); // D, W, H
		float const width(0.3*door_width), height(width*sz.z/sz.y), depth(width*sz.x/sz.y); // scale to the width of the wall lamp
		float const z1(door_bcube.z1() + 0.7*door_bcube.dz());
		float const lamp_pos(door_bcube.d[!dim][side] + (side ? 1.0 : -1.0)*0.6*width);
		cube_t lamp(c);
		lamp.d[dim][dir] = c.d[dim][!dir] + (dir ? 1.0 : -1.0)*depth;
		set_cube_zvals(lamp, z1, (z1 + height));
		set_wall_width(lamp, lamp_pos, 0.5*width, !dim);
		unsigned flags(base_flags);
		if (rgen.rand_bool()) {flags |= RO_FLAG_LIT;} // light is on 50% of the time
		objs.emplace_back(lamp, TYPE_WALL_LAMP, room_id, dim, dir, flags, tot_light_amt);
		if (objs.back().is_lit()) {ext_lights.emplace_back(lamp.get_cube_center(), 20.0*width, WALL_LAMP_COLOR);}
	}
	if (has_porch() && rgen.rand_bool()) { // maybe add a package on the porch
		// find the front door
		float const wall_thickness(get_wall_thickness());
		cube_t front_door;

		for (auto const &d : doors) {
			cube_t dbc(d.get_bcube());
			dbc.expand_by_xy(wall_thickness);
			if (porch.intersects(dbc)) {front_door = dbc;}
		}
		if (!front_door.is_all_zeros()) { // should never fail?
			vector3d sz; // half size relative to window_vspacing
			gen_crate_sz(sz, rgen, get_window_vspace());
			bool const door_dim(front_door.dy() < front_door.dx());
			cube_t place_area(front_door), avoid(front_door);
			place_area.expand_by(vector3d(4.0*sz.x, 4.0*sz.y, 0.0));
			cube_t porch_shrink(porch);
			porch_shrink.expand_by_xy(-1.0*wall_thickness); // add space for window sill (approximate)
			place_area.intersect_with_cube(porch_shrink);
			avoid.expand_in_dim(door_dim, 4.0*sz[door_dim]); // add clearance in front of the door; package must be left to one side

			if (place_area.dx() > 2.5*sz.x && place_area.dy() > 2.5*sz.y) { // check if large enough
				for (unsigned n = 0; n < 50; ++n) {
					point const pos(gen_xy_pos_in_area(place_area, sz, rgen, porch.z2()));
					cube_t const box(get_cube_height_radius(pos, sz, 2.0*sz.z)); // multiply by 2 since this is a size rather than half size/radius
					if (box.intersects(avoid)) continue;
					bool skip(0);

					for (cube_t const &part : parts) { // check for intersection with porch support
						if (box.intersects(part)) {skip = 1; break;}
					}
					if (skip) continue;
					objs.emplace_back(box, TYPE_BOX, room_id, rgen.rand_bool(), 0, base_flags, tot_light_amt, SHAPE_CUBE, gen_box_color(rgen));
					break; // done
				} // for n
			}
		}
	}
}

room_t::room_t(cube_t const &c, unsigned p, unsigned nl, bool is_hallway_, bool is_office_, bool is_sec_bldg_) :
	cube_t(c), no_geom(is_hallway_), is_hallway(is_hallway_), is_office(is_office_), is_sec_bldg(is_sec_bldg_), part_id(p), num_lights(nl) // no geom in hallways
{
	if      (is_sec_bldg) {assign_all_to(RTYPE_GARAGE);} // or RTYPE_SHED - will be set later
	else if (is_hallway)  {assign_all_to(RTYPE_HALL  );}
	else if (is_office)   {assign_all_to(RTYPE_OFFICE);}
	else if (has_stairs)  {assign_all_to(RTYPE_STAIRS);} // not really correct since has_stairs is now a per-floor bit flag, but this will likely be overwritten later anyway
	else                  {assign_all_to(RTYPE_NOTSET, 0);} // locked=0
}
void room_t::assign_all_to(room_type rt, bool locked) {
	for (unsigned n = 0; n < NUM_RTYPE_SLOTS; ++n) {rtype[n] = rt;}
	if (locked) {rtype_locked = 0xFF;} // room type is locked on all floors
}
void room_t::assign_to(room_type rt, unsigned floor, bool locked) {
	// room types are only tracked up to the 4th floor, and every floor above that has the same type as the 4th floor; good enough for houses at least
	floor = wrap_room_floor(floor);
	if (rtype[floor] == RTYPE_BATH) return; // assign unless already set to a bathroom, since we need that for has_room_of_type(RTYPE_BATH)
	rtype[floor] = rt;
	if (locked) {rtype_locked |= (1 << floor);} // lock this floor
}
bool room_t::has_room_of_type(room_type type) const {
	for (unsigned n = 0; n < NUM_RTYPE_SLOTS; ++n) {
		if (rtype[n] == type) return 1;
	}
	return 0;
}

