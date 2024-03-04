// 3D World - Building Interior Room Assignment, etc.
// by Frank Gennari 4/30/2020

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t

extern int world_mode;
extern object_model_loader_t building_obj_model_loader;

void setup_bldg_obj_types();
bool cube_int_tiled_terrain_trees(cube_t const &c);
bool get_wall_quad_window_area(vect_vnctcc_t const &wall_quad_verts, unsigned i, cube_t &c, float &tx1, float &tx2, float &tz1, float &tz2);
void get_balcony_pillars(room_object_t const &c, float ground_floor_z1, cube_t pillar[2]);
void expand_convex_polygon_xy(vect_point &points, point const &center, float expand);
bool is_pool_tile_floor(room_object_t const &obj);
void invalidate_tile_smap_in_region(cube_t const &region, bool repeat_next_frame=0);
void get_obj_drawers_or_doors(room_object_t const &obj, vect_cube_t &drawers, room_object_t &drawers_part, float &drawer_extend);


unsigned light_ix_assign_t::get_ix_for_light(cube_t const &c, bool walls_not_shared) {
	point2d<float> const pos(c.x1(), c.y1());

	if (!walls_not_shared) { // search for existing stack
		for (auto i = cur.begin(); i != cur.end(); ++i) {
			if (i->first == pos) return i->second; // existing light is part of the same stack and is valid to return
		}
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
bool building_t::can_be_bedroom_or_bathroom(room_t const &room, unsigned floor_ix, bool skip_conn_check) const { // check room type and existence of exterior door
	if (room.has_stairs_on_floor(floor_ix) || room.has_elevator || room.is_hallway || room.is_office || room.is_sec_bldg) return 0; // no bed/bath in these cases
	
	if (has_ext_door_this_floor(room.z1(), floor_ix)) {
		// run special logic for bedrooms and bathrooms (private rooms) on the first floor of a house
		if (is_room_adjacent_to_ext_door(room)) return 0; // door to house does not open into a bedroom/bathroom
		if (skip_conn_check) return 1;
		bool const is_multi_floor(room.dz() > 1.5*get_window_vspace());
		bool const has_stairs(is_multi_floor && !interior->stairwells.empty()); // more than one floor and stairs placement didn't fail

		// check paths if there are stairs or an interior garage; skip for single floor houses since there may be no valid bed/bath placement with these constraints
		if (is_multi_floor && (has_stairs || has_int_garage)) {
			// determine if this room is on the shortest path from an exterior door to the stairs or garage; if so, it can't be a bedroom or bathroom;
			// okay, that's not easy/fast to do, so determine if there is any path from the exterior door to the stairs that doesn't go through this room;
			// this won't work when there are two paths from the door to the stairs and this room is only on one of the paths, so we could put a BR/BR on both paths
			int cur_room(-1);
			static vector<unsigned> door_rooms, stairs_rooms;
			door_rooms.clear();
			stairs_rooms.clear();

			for (unsigned i = 0; i < interior->rooms.size(); ++i) {
				room_t const &r(interior->rooms[i]);
				if (r == room) {cur_room = i; continue;} // this room; we know it can't have stairs or an exterior door
				if (r.is_sec_bldg || r.z2() <= ground_floor_z1) continue; // skip basement rooms, garages, and sheds
				if (r.has_stairs_on_floor(floor_ix))  {stairs_rooms.push_back(i);}
				if (is_room_adjacent_to_ext_door(r))  {door_rooms  .push_back(i);}
			}
			if (is_multi_floor && stairs_rooms.empty()) {
				if (!has_missing_stairs) {cout << "Building with missing stairs: " << bcube.str() << endl;}
				has_missing_stairs = 1;
			}
			assert(cur_room >= 0); // must be found
			// Note: stairs_rooms can be empty if there are only basement stairs
			// Note: doors can be empty if door placement failed, which can happen for rotated and non-cube buildings
			for (unsigned d : door_rooms) {
				for (unsigned s : stairs_rooms) {
					if (!are_rooms_connected_without_using_room(d, s, cur_room)) return 0;
				}
				if (has_int_garage && !are_rooms_connected_without_using_room(d, interior->garage_room, cur_room)) return 0;
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
	return (r.is_single_floor ? 1 : calc_num_floors(r, window_vspacing, floor_thickness));
}

void building_t::setup_courtyard() {
	if (!has_courtyard) return;
	assert(has_room_geom());
	courtyard_t &courtyard(interior->room_geom->courtyard);
	courtyard.set_to_zeros();
	set_cube_zvals(courtyard, ground_floor_z1, parts.front().z2());
	
	// find couryard bounds formed from all non-exterior part edges; assume courtyard contains the building center
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		if (is_basement(p)) continue; // skip basement
		if (p->z1() < ground_floor_z1) continue; // not on ground floor; shouldn't happen?
		if (p->x2() < bcube.xc()) {courtyard.x1() = p->x2();}
		if (p->x1() > bcube.xc()) {courtyard.x2() = p->x1();}
		if (p->y2() < bcube.yc()) {courtyard.y1() = p->y2();}
		if (p->y1() > bcube.yc()) {courtyard.y2() = p->y1();}
	} // for p
	assert(courtyard.is_strictly_normalized());
	if (!has_courtyard_door) return; // return here with courtyard.door_ix and courtyard.room_ix left as -1
	// find courtyard door
	cube_t door_bc;
	unsigned num_doors(0), num_rooms(0);

	for (auto d = doors.begin(); d != doors.end(); ++d) { // check exterior doors; likely the last one, unless there's a roof exit door
		cube_t bc(d->get_bcube());
		bool const dim(bc.dy() < bc.dx());
		bc.expand_in_dim(dim, get_wall_thickness());
		if (bc.intersects(courtyard)) {courtyard.door_ix = int16_t(d - doors.begin()); door_bc = bc; ++num_doors;}
	}
	assert(num_doors == 1);

	// find courtyard room
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (r->intersects(door_bc)) {courtyard.room_ix = int16_t(r - interior->rooms.begin()); ++num_rooms;}
	}
	assert(num_rooms == 1);
}
bool building_t::point_in_courtyard(point const &pos_bs) const {
	return (has_courtyard && has_room_geom() && interior->room_geom->courtyard.contains_pt(pos_bs));
}

void building_t::gen_room_details(rand_gen_t &rgen, unsigned building_ix) {

	assert(has_room_geom());
	setup_bldg_obj_types(); // initialize object types if not already done
	setup_courtyard();
	//highres_timer_t timer("Gen Room Details");
	// Note: people move from room to room, so using their current positions for room object generation is both nondeterministic and unnecessary
	vect_cube_t blockers, valid_lights, segs; // blockers are used for fireplaces
	vect_room_object_t &objs(interior->room_geom->objs);
	vector<room_t> &rooms(interior->rooms);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness), wall_thickness(get_wall_thickness());
	float const light_thick(0.025*window_vspacing), def_light_size(0.1*window_vspacing);
	interior->room_geom->obj_scale = window_vspacing; // used to scale room object textures
	unsigned tot_num_rooms(0), num_bathrooms(0), num_bedrooms(0), num_storage_rooms(0);
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {tot_num_rooms += calc_num_floors_room(*r, window_vspacing, floor_thickness);}
	objs.reserve(tot_num_rooms); // placeholder - there will be more than this many
	float const extra_bathroom_prob((is_house ? 2.0 : 1.0)*0.02*min((int(tot_num_rooms) - 4), 20));
	room_obj_shape const light_shape(is_house ? SHAPE_CYLIN : SHAPE_CUBE);
	unsigned cand_bathroom(rooms.size()); // start at an invalid value
	unsigned added_kitchen_mask(0), added_living_mask(0), added_bath_mask(0); // per-floor
	unsigned added_bathroom_objs_mask(0);
	bool added_bedroom(0), added_library(0), added_dining(0), added_laundry(0), added_basement_utility(0), added_fireplace(0), added_pool_room(0);
	light_ix_assign_t light_ix_assign;
	interior->create_fc_occluders(); // not really part of room geom, but needed for generating and drawing room geom, so we create them here
	has_int_fplace = 0; // reset for this generation

	if (rooms.size() > 1) { // choose best room assignments for required rooms; if a single room, skip this step
		float min_score(0.0); // lower is better

		// Note: assigning cand_bathroom when has_pri_hall() is not strictly necessary, but may help add a bathroom to an upper stacked part
		for (auto r = rooms.begin(); r != rooms.end(); ++r) {
			if (r->is_sec_bldg) continue; // garage/shed excluded - not a normal room
			if (has_basement() && r->part_id == (int)basement_part_ix) continue; // skip the basement
			if (r->is_single_floor && r->dz() > 1.5*window_vspacing)   continue; // no tall ceiling rooms
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
				if (num_floors > 1 && r->has_stairs_on_floor(0)) {score *= 4.0;} // penalty for ground floors stairs connecting to the basement or a stacked part
				if (min_score == 0.0 || score < min_score) {cand_bathroom = (r - rooms.begin()); min_score = score;} // lower score is better
			}
		} // for r
	}
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {
		bool const is_basement(has_basement() && r->part_id == (int)basement_part_ix); // includes extended basement and parking garage
		float light_amt(is_basement ? 0.0f : window_vspacing*r->get_light_amt()); // exterior light: multiply perimeter/area by window spacing to make unitless; none for basement rooms
		if (!is_house && r->is_hallway) {light_amt *= 2.0;} // double the light in office building hallways because they often connect to other lit hallways
		float const floor_height(r->is_single_floor ? r->dz() : window_vspacing); // secondary buildings are always one floor
		unsigned const num_floors(calc_num_floors_room(*r, floor_height, floor_thickness)), room_id(r - rooms.begin());
		unsigned const min_br(multi_family ? num_floors : 1); // multi-family house requires one per floor; can apply to both bedrooms and bathrooms
		point room_center(r->get_cube_center());

		// determine light pos and size for this stack of rooms
		float const dx(r->dx()), dy(r->dy());
		bool const room_dim(dx < dy); // longer room dim
		room_type const init_rtype_f0(r->get_room_type(0));
		bool const is_parking_garage(init_rtype_f0 == RTYPE_PARKING   ); // all floors should be parking garage
		bool const is_unfinished    (init_rtype_f0 == RTYPE_UNFINISHED); //  // unfinished room, for example in a non-cube shaped office building
		bool const is_swim_pool_room(init_rtype_f0 == RTYPE_SWIM); // room with a swimming pool
		bool const is_retail_room   (init_rtype_f0 == RTYPE_RETAIL);
		bool const is_ext_basement(r->is_ext_basement()), is_backrooms(r->is_backrooms());
		float light_density(0.0), light_size(def_light_size); // default size for houses
		unsigned const room_objs_start(objs.size());
		unsigned nx(1), ny(1); // number of lights in X and Y for this room

		if (!is_cube()) { // somewhat more lights for non-cube shaped building pie slices
			light_density = 0.4;
		}
		else if (r->is_hallway) { // light size varies by hallway size
			float const room_size(min(dx, dy)); // normalized to hallway width
			light_size = max(0.06f*room_size, 0.67f*def_light_size);
		}
		else if (r->is_office || r->has_skylight) { // more lights for large offices; light size varies by office size; parking garages are handled later
			light_density = 0.5;
			float const room_size(dx + dy); // normalized to office size
			light_size = max(0.015f*room_size, 0.67f*def_light_size);
		}
		else if (is_backrooms) { // large office basement room
			light_density = 0.55;
			light_size   *= 0.75; // smaller
		}
		else if (is_retail_room) { // more lights in the shorter dim
			light_size *= 0.7*pow(double(retail_floor_levels), 0.4); // smaller; increase size for taller rooms
			nx = max(1U, unsigned((room_dim ? 0.7 : 0.4)*dx/window_vspacing));
			ny = max(1U, unsigned((room_dim ? 0.4 : 0.7)*dy/window_vspacing));
		}
		else if (is_swim_pool_room) {
			light_density = 0.4;
		}
		else if (r->is_single_floor) {
			light_size *= sqrt(r->dz()/window_vspacing); // larger lights for taller rooms
		}
		if (light_density > 0.0) { // uniform 2D grid of lights
			nx = max(1U, unsigned(light_density*dx/window_vspacing));
			ny = max(1U, unsigned(light_density*dy/window_vspacing));
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
		colorRGBA chair_color(chair_colors[(13*r->part_id + 123*tot_num_rooms + 617*mat_ix + 1367*num_floors) % NUM_CHAIR_COLORS]);
		light_ix_assign.next_room();
		rand_gen_t room_rgen(rgen); // shared across all floors of this room
		// select light color for this room
		colorRGBA color;
		if (r->is_ext_basement_conn()) {color = RED;}
		else if (is_house)          {color = get_light_color_temp(0.4);} // house - yellowish
		else if (is_parking_garage) {color = get_light_color_temp_range(0.2, 0.5, rgen);} // parking garage - yellow-white
		else if (is_backrooms)      {color = get_light_color_temp_range(0.2, 0.4, rgen);} // backrooms - yellow-white
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
			unsigned const floor_objs_start(objs.size()); // needed for backrooms lights
			bool is_lit(0), has_light(1), light_dim(room_dim), wall_light(0), has_stairs(has_stairs_this_floor), top_of_stairs(has_stairs && top_floor);
			float light_delta_z(0.0);
			vect_cube_t rooms_to_light;

			if (is_parking_garage) { // parking garage; added first because this sets the number of lights
				assert(r->interior == 1);
				add_parking_garage_objs(rgen, *r, room_center.z, room_id, f, num_floors, nx, ny, light_delta_z);
			}
			else if (is_backrooms) { // should be single floor only
				add_backrooms_objs(rgen, *r, room_center.z, room_id, f, rooms_to_light);
			}
			else if (is_retail_room) {
				add_retail_room_objs(rgen, *r, room_center.z, room_id, light_ix_assign);
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
			else if (r->is_sec_bldg) {is_lit = 0;} // garage and shed lights start off
			else {
				// 50% of lights are on, 75% for top of stairs, 100% for non-basement hallways, 100% for parking garages and backrooms
				is_lit  = ((r->is_hallway && !is_basement) || is_parking_garage || is_backrooms || is_retail_room || ((rgen.rand() & (top_of_stairs ? 3 : 1)) != 0));
				is_lit |= (r->is_ext_basement_conn() || (r->is_ext_basement() && r->intersects(get_basement()))); // ext basement conn or primary hallway

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
			// add a light to the ceiling of this room if there's space (always for top of stairs)
			unsigned const lcheck_start_ix(is_backrooms ? floor_objs_start : objs.size()); // must check lights vs. backrooms walls and pillars
			set_cube_zvals(light, (light_z2 - light_thick), light_z2);
			valid_lights.clear();

			if (num_lights > 1) { // r->is_hallway or ext basement
				// hallway: place a light on each side (of the stairs if they exist), and also between stairs and elevator if there are both
				if (r->is_hallway && r->has_elevator && r->has_stairs == 255) {max_eq(num_lights, (2U + r->has_elevator));} // hall with elevator + stairs on all floors; 3+ lights
				min_eq(num_lights, 6U);
				float const offsets[6] = {0.0, -0.2, -0.3, -0.36, -0.4, -0.43}, steps[6] = {0.0, 0.4, 0.3, 0.24, 0.2, 0.172}; // indexed by num_lights-1
				float const hallway_len(r->get_sz_dim(light_dim));
				unsigned num[2] = {1, 1};
				num[light_dim] = num_lights;

				for (unsigned d = 0; d < num_lights; ++d) {
					float const delta((offsets[num_lights-1] + d*steps[num_lights-1])*hallway_len);
					cube_t hall_light(light);
					hall_light.translate_dim(light_dim, delta);
					try_place_light_on_ceiling(hall_light, *r, room_dim, fc_thick, 0, 0, num[0], num[1], lcheck_start_ix, valid_lights, rgen); // allow_rot=0, allow_mult=0
				}
				if (r->is_hallway && has_pri_hall()) { // make sure to place lights between all stairs and elevators
					float const light_len(light.get_sz_dim(room_dim));
					cube_t centerline(*r);
					set_wall_width(centerline, r->get_center_dim(!room_dim), wall_thickness, !room_dim);
					set_cube_zvals(centerline, floor_zval, light_z2); // clip to the range of this floor
					segs.clear();
					segs.push_back(centerline);

					for (elevator_t  const &e : interior->elevators ) {
						if (centerline.intersects(e)) {subtract_cube_from_cubes(e, segs);}
					}
					for (stairwell_t const &s : interior->stairwells) {
						if (centerline.intersects(s)) {subtract_cube_from_cubes(s, segs);}
					}
					for (cube_t const &seg : segs) {
						if (seg.d[!room_dim][0] != centerline.d[!room_dim][0] || seg.d[!room_dim][1] != centerline.d[!room_dim][1]) continue; // not full width
						if (seg.get_sz_dim(room_dim) < 1.5*light_len) continue; // too short

						if (!has_bcube_int_xy(seg, valid_lights)) { // no light yet on this segment
							cube_t light2(light);
							set_wall_width(light2, seg.get_center_dim(room_dim), 0.5*light_len, room_dim); // place on segment center
							if (is_light_placement_valid(light2, *r, fc_thick) && !overlaps_other_room_obj(light2, lcheck_start_ix)) {valid_lights.push_back(light2);}
						}
					}
				}
			}
			else if (nx > 1 || ny > 1) { // office, parking garage, or backrooms with multiple lights
				vector3d const shrink(0.5*light.dx()*sqrt((nx - 1)/nx), 0.5*light.dy()*sqrt((ny - 1)/ny), 0.0);
				float xstep(dx/nx), ystep(dy/ny), xs(-0.5f*dx + 0.5*xstep), ys(-0.5f*dy + 0.5*ystep);

				if (is_retail_room) { // custom logic to align lights to aisles
					xs -= 0.25*xstep; ys -= 0.25*ystep;
					xstep = dx/(nx - 0.5); ystep = dy/(ny - 0.5);
				}
				for (unsigned y = 0; y < ny; ++y) {
					for (unsigned x = 0; x < nx; ++x) {
						cube_t cur_light(light);
						cur_light.expand_by_xy(-shrink);
						cur_light.translate(point((xs + x*xstep), (ys + y*ystep), 0.0));
						try_place_light_on_ceiling(cur_light, *r, room_dim, fc_thick, 0, 0, nx, ny, lcheck_start_ix, valid_lights, rgen); // allow_rot=0, allow_mult=0
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
				if (!wall_light) {try_place_light_on_ceiling(light, *r, room_dim, fc_thick, 1, 1, 1, 1, lcheck_start_ix, valid_lights, rgen);} // allow_rot=1, allow_mult=1
				if (!valid_lights.empty()) {light_obj_ix = objs.size();} // this will be the index of the light to be added later
			}
			rand_gen_t rgen_lights(rgen); // copy state so that we don't modify rgen
			unsigned const objs_start_inc_lights(objs.size());
			bool const walls_not_shared(is_backrooms); // multi-floor backrooms have different walls and can't share the light stack

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
				light_obj.obj_id = light_ix_assign.get_ix_for_light(l, walls_not_shared);
				unsigned const flicker_mod(is_parking_garage ? 50 : (is_ext_basement ? 20 : 0)); // 2% chance for parking garage, 5% chance for ext basement
				
				if (flicker_mod > 0 && (((rgen_lights.rand() + 3*f)%flicker_mod) == 13)) {light_obj.flags |= RO_FLAG_BROKEN;} // maybe make this a flickering light
				else if (is_ext_basement && valid_lights.size() == 1 && (rgen_lights.rand() & 7) == 0) { // broken ext basement light; not for hallways with multiple lights
					light_obj.flags |= RO_FLAG_BROKEN2;
					light_obj.flags &= ~RO_FLAG_LIT; // off by default
				}
				objs.emplace_back(light_obj);
			} // for l
			if (is_lit) {r->lit_by_floor |= (1ULL << (f&63));} // flag this floor as being lit (for up to 64 floors)

			if (is_backrooms) {
				room_object_t const ref_light(light, TYPE_LIGHT, room_id, room_dim, 0, (flags | RO_FLAG_NOCOLL), light_amt, light_shape, color);
				add_missing_backrooms_lights(rgen, room_center.z, room_id, floor_objs_start, objs_start_inc_lights, ref_light, rooms_to_light, light_ix_assign);
				continue; // nothing else to add
			}
			if (is_parking_garage || is_retail_room) continue; // generated above, done; no outlets or light switches
			if (is_unfinished) continue; // no objects for now; if adding objects later, need to make sure they stay inside the building bounds
			float tot_light_amt(light_amt); // unitless, somewhere around 1.0
			if (is_lit) {tot_light_amt += r->light_intensity;}
			bool const is_garage_or_shed(r->is_garage_or_shed(f));
			bool const is_ground_floor_part(!is_basement && r->z1() <= ground_floor_z1);
			bool const is_ground_floor(is_ground_floor_part && (1 << f) & floor_ext_door_mask); // really this is a check for doors on this floor
			bool const is_entry_floor(is_ground_floor || (is_ground_floor_part && multi_family)); // for placing entry level rooms/objs (fireplace, living, dining, kitchen, etc.)
			unsigned const objs_start(wall_light ? objs_start_inc_lights : objs.size()); // wall light counts as an object since it must be avoided
			rgen.rand_mix();
			blockers.clear(); // clear for this new room

			if (r->is_ext_basement_conn()) { // room connecting extended basements of two buildings
				bool const conn_dim(r->interior == 4);
				// add blockers at both room ends to avoid placing objects there, since there's a door to the other building blocking it
				for (unsigned d = 0; d < 2; ++d) {
					cube_t blocker(*r);
					blocker.d[conn_dim][!d] = blocker.d[conn_dim][d] + (d ? -1.0 : 1.0)*2.0*wall_thickness;
					objs.emplace_back(blocker, TYPE_BLOCKER, room_id, conn_dim, d, RO_FLAG_INVIS);
				}
			}
			if (num_sides > 4 && !is_basement) { // non-cube/triangle building pie slice - add room pillars
				// these could be added as walls or added for the bottom floor and extend through all levels, though that would require more changes
				add_office_pillars(rgen, *r, room_center.z, room_id, f, valid_lights, blockers);
			}
			if (r->no_geom || is_garage_or_shed || is_swim_pool_room) {
				if (is_garage_or_shed) {
					if (init_rtype_f0 == RTYPE_GARAGE) {
						room_center.z = add_flooring(*r, room_center.z, room_id, tot_light_amt, FLOORING_CONCRETE);
						add_garage_objs(rgen, *r, room_center.z, room_id, tot_light_amt);
					}
					// is there enough clearance between shelves and a car parked in the garage? there seems to be in all the cases I've seen
					add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
				}
				if (is_swim_pool_room) {
					assert(is_ext_basement); // for now, only in extended basements
					add_swimming_pool_room_objs(rgen, *r, room_center.z, room_id, tot_light_amt);
				}
				add_outlets_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor, is_basement);
				if (has_light) {add_light_switches_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor, is_basement);} // shed, garage, or hallway
				if (r->is_hallway && is_ext_basement) {add_false_door_to_extb_hallway_if_needed(*r, room_center.z, room_id);}

				if (is_house && r->is_hallway) { // allow pictures, rugs, and light switches in the hallways of houses
					hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement);
					if (rgen.rand_bool()) {add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);} // 50% of the time; not all rugs will be placed
				}
				if (!is_house && r->is_hallway && *r == pri_hall) { // office building primary hallway
					add_pri_hall_objs(rgen, room_rgen, *r, room_center.z, room_id, tot_light_amt, f);
					if (is_ground_floor) {r->assign_to(RTYPE_LOBBY, f);} // first floor primary hallway, make it the lobby
				}
				if (has_stairs_this_floor && r->get_room_type(f) == RTYPE_NOTSET) {r->assign_to(RTYPE_STAIRS, f);}
				continue; // no other geometry for this room
			}
			//if (has_stairs && !pri_hall.is_all_zeros()) continue; // no other geometry in office building base part rooms that have stairs
			unsigned const floor_mask(1<<f);
			// must be a BR if cand bathroom, and BR not already placed; applies to all floors of this room; if multi-family, we check for a BR prev placed on this floor
			bool const must_be_bathroom(room_id == cand_bathroom && (multi_family ? !(added_bath_mask & floor_mask) : (num_bathrooms == 0)));
			bool const is_tall_room(r->is_single_floor && r->dz() > 1.5*window_vspacing);
			bool added_tc(0), added_desk(0), added_obj(0), can_place_onto(0), no_whiteboard(0);
			bool is_bathroom(0), is_bedroom(0), is_kitchen(0), is_living(0), is_dining(0), is_storage(0), is_utility(0), no_plants(0), is_play_art(0);
			unsigned num_chairs(0);
			// unset room type if not locked on this floor during floorplanning; required to generate determinstic room geom
			if (!r->is_rtype_locked(f)) {r->assign_to(RTYPE_NOTSET, f);}

			// place room objects
			bool const added_living(added_living_mask & floor_mask);
			bool const allow_br(!is_house || must_be_bathroom || f > 0 || num_floors == 1 || (rgen.rand_float() < 0.33f*(added_living + (added_kitchen_mask&1) + 1))); // bed/bath
			bool is_office_bathroom(is_room_office_bathroom(*r, room_center.z, f)), has_fireplace(0);
			
			if (has_chimney == 2 && !is_basement && is_entry_floor && !added_fireplace) { // handle fireplaces on the first floor
				has_fireplace = added_fireplace = maybe_add_fireplace_to_room(rgen, *r, blockers, room_center.z, room_id, tot_light_amt);
			}
			if (is_office_bathroom) { // bathroom is already assigned
				added_obj = is_bathroom = added_bathroom = no_whiteboard =
					add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement, added_bathroom_objs_mask); // add bathroom
			}
			else if (!is_house && f == 0) { // office building special pre-assigned first floor rooms; can be in a stacked part
				if (init_rtype_f0 == RTYPE_UTILITY) {
					added_obj = no_whiteboard = no_plants = is_utility = add_office_utility_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
				else if (init_rtype_f0 == RTYPE_SERVER) {
					added_obj = no_whiteboard = no_plants = add_server_room_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
				else if (init_rtype_f0 == RTYPE_SECURITY) {
					added_obj = no_whiteboard = no_plants = add_security_room_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
			}
			// bedroom or bathroom case; need to check first floor even if must_be_bathroom
			if (!added_obj && allow_br && !is_tall_room && can_be_bedroom_or_bathroom(*r, f)) {
				// Note: num_bedrooms is summed across all floors, while num_bathrooms is per-floor
				// Note: min_br is applied to bedrooms, but could be applied to bathrooms in the same way
				bool const pref_sec_bath(is_house && num_bathrooms == 1 && num_bedrooms > min_br && rooms.size() >= 6 && !must_be_bathroom && !has_fireplace && can_be_bathroom(*r));
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
				// office buildings can also have kitchens, even on non-ground floors; no tall room kitchens because the cabinets and stove hood have no ceiling to connect to
				if (!(added_kitchen_mask & floor_mask) && (!is_house || is_entry_floor) && !is_basement && !is_tall_room) {
					if (added_tc || (is_house && (r+1) == rooms.end())) { // make it a kitchen if it's the last room in a house, even if there's no table
						if (add_kitchen_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, added_living)) {
							r->assign_to(RTYPE_KITCHEN, f);
							added_kitchen_mask |= floor_mask;
							is_kitchen = added_obj = 1;
						}
					}
				}
			}
			if (!added_obj && !added_pool_room && is_house && is_basement && add_pool_room_objs(rgen, *r, room_center.z, room_id, tot_light_amt)) { // pool room
				r->assign_to(RTYPE_POOL, f);
				added_pool_room = added_obj = 1;
			}
			// if we haven't added any objects yet, and this room is an interior office on the first floor or basement, make it a storage room 50% of the time; at most 4x
			if (!added_obj && num_storage_rooms <= 4 && (is_basement || (r->is_office && r->interior && f == 0 /*&& r->z1() == ground_floor_z1*/)) && rgen.rand_bool()) {
				added_obj = no_whiteboard = is_storage = add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
				if (added_obj) {r->assign_to(RTYPE_STORAGE, f); ++num_storage_rooms;}
			}
			if (!added_obj && (!is_basement || rgen.rand_bool())) { // try to place a desk if there's no table, bed, etc.
				added_obj = can_place_onto = added_desk = add_office_objs(rgen, *r, blockers, chair_color, room_center.z, room_id, f, tot_light_amt, objs_start, is_basement);
				if (added_obj && !has_stairs_this_floor) {r->assign_to((is_house ? (room_type)RTYPE_STUDY : (room_type)RTYPE_OFFICE), f);} // or other room type - may overwrite below
			}
			// Note: added_obj is not set to 1 below this point
			if (is_house && (added_tc || added_desk) && !is_kitchen && is_entry_floor) { // don't add second living room unless we added a kitchen and have enough rooms
				if ((!added_living && !r->has_center_stairs && rooms.size() >= 8 && (added_kitchen_mask || rgen.rand_bool())) || is_room_an_exit(*r, room_id, room_center.z)) {
					// add a living room on the ground floor if it has a table or desk but isn't a kitchen
					if (add_livingroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start)) {
						is_living = 1;
						added_living_mask |= floor_mask;
						r->assign_to(RTYPE_LIVING, f);
						// if living room is next to a dining room or kitchen, do we want to remove the door and wall between the rooms? would need to regen VBOs
					}
				}
			}
			if (is_house && added_tc && num_chairs > 0 && !is_living && !is_kitchen) { // room with table and chair that's not a kitchen
				if (is_entry_floor) { // dining room, must be on the ground floor or an entry floor
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
			if (!is_house && !added_obj && has_stairs_this_floor) { // office building "office" with stairs and no furniture
				r->assign_to(RTYPE_STAIRS, f);
				no_whiteboard = 1;
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
			// office building part with primary hallway, first floor of first non-retail part
			else if (has_pri_hall() && r->part_id == (has_retail() ? 1 : 0) && f == 0 && added_desk) {
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
				else if (!added_obj && !has_fireplace) { // unassigned empty room - make it a storage room
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
			if (is_bathroom) {added_bath_mask |= floor_mask;}
		} // for f (floor)
		if (added_bathroom) {++num_bathrooms;}

		if (r->interior) { // tag objects as interior if room is interior
			for (auto i = objs.begin() + room_objs_start; i != objs.end(); ++i) {i->flags |= RO_FLAG_INTERIOR;}
		}
	} // for r (room)
	if (is_house) {interior->assign_master_bedroom(window_vspacing, floor_thickness);}
	add_padlocks(rgen);

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
	add_exterior_ac_pipes(rgen);
	unsigned const ext_objs_start(objs.size());
	vect_cube_t balconies;
	ext_steps.clear(); // clear prev value in case this building's interior is recreated
	maybe_add_fire_escape  (rgen); // or ladder
	add_balconies          (rgen, balconies);
	add_gutter_downspouts  (rgen, balconies);
	add_exterior_door_items(rgen);
	if (has_chimney) {add_chimney_cap(rgen);}
	if (!is_rotated()) {add_ext_door_steps(ext_objs_start);} // must be after adding balconies and fire escape
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

void building_t::add_chimney_cap(rand_gen_t &rgen) {
	assert(interior);
	if (rgen.rand_float() < 0.25) return; // 25% chance of no chimney cap
	cube_t ccap(get_chimney()); // start with full chimney, then place the cap on top
	set_cube_zvals(ccap, ccap.z2(), (ccap.z2() + 0.2*get_window_vspace()));
	interior->room_geom->objs.emplace_back(ccap, TYPE_CHIM_CAP, 0); // room_id=0
	set_obj_id(interior->room_geom->objs); // used for the style
	bcube.union_with_cube(ccap); // extend bcube to include the chimnmey cap so that it's always drawn
}

void building_t::maybe_add_fire_escape(rand_gen_t &rgen) { // or ladder
	if (!is_house) return; // houses only for now
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), wall_thickness(get_wall_thickness()), fe_height(4.25*window_vspacing);

	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_FESCAPE)) {
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			unsigned const num_floors(calc_num_floors(*p, window_vspacing, floor_thickness));
			// our hard-coded fire escape model is designed for a 5 story building; but the max number of floors for a 'house' is 5-6 anyway, which makes them relatively rare
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
				fe_bc.d[dim][0]    = fe_bc.d[dim][1] = p->d[dim][dir];
				fe_bc.d[dim][dir] += (dir ? 1.0 : -1.0)*fe_depth;
				cube_t fe_bc_exp(fe_bc);
				fe_bc_exp.expand_by(wall_thickness); // needed for exterior door check
				if (has_bcube_int_no_adj(fe_bc, parts))              continue; // check for intersection with other parts, in particular the chimney and fireplace
				if (has_driveway() && fe_bc.intersects_xy(driveway)) continue; // skip if intersects driveway or garage
				if (cube_int_ext_door(fe_bc_exp))                    continue; // check exterior doors
				interior->room_geom->objs.emplace_back(fe_bc, TYPE_FESCAPE, 0, dim, dir, RO_FLAG_EXTERIOR, 1.0, SHAPE_CUBE, BLACK); // room_id=0
				details.emplace_back(fe_bc, DETAIL_OBJ_COLLIDER);
				union_with_coll_bcube(fe_bc);
				return; // success/done
			} // for d
		} // for p
	}
	// if no fire escape was added, maybe add a TYPE_LADDER to the roof between windows; only for hipped roofs, to avoid the gutter; skip houses with stacked parts
	if (roof_type == ROOF_TYPE_HIPPED && rgen.rand_bool()) {
		cube_t const &part(parts[0]); // add to first/primary part
		unsigned const pref_dim_dir(rgen.rand() & 3);
		float const window_h_border(get_window_h_border()), hwidth(0.3*get_doorway_width()), depth(0.3*hwidth);
		bool const check_stacked(real_num_parts > 1 && parts[1].z1() > parts[0].z1());

		for (unsigned d = 0; d < 4; ++d) {
			unsigned const dd((d + pref_dim_dir) & 3);
			bool const dim(dd >> 1), dir(dd & 1);
			if (part.d[dim][dir] != bcube.d[dim][dir]) continue; // not on the building bcube - could intersect another part, porch, etc.
			if (part.get_sz_dim(!dim) < 8.0*hwidth)    continue; // wall is too narrow

			for (unsigned n = 0; n < 10; ++n) { // make 10 tries
				cube_t bc(part); // copy zvals
				set_wall_width(bc, rgen.rand_uniform((part.d[!dim][0] + 2.0*hwidth), (part.d[!dim][1] - 2.0*hwidth)), hwidth, !dim);
				float const window_hspacing(get_hspacing_for_part(part, !dim));
				if (is_val_inside_window(part, !dim, bc.d[!dim][0]-wall_thickness, window_hspacing, window_h_border)) continue; // check for window intersection
				if (is_val_inside_window(part, !dim, bc.d[!dim][1]+wall_thickness, window_hspacing, window_h_border)) continue; // check for window intersection
				bc.d[dim][0]    = bc.d[dim][1] = part.d[dim][dir]; // at wall
				bc.d[dim][dir] += (dir ? 1.0 : -1.0)*depth; // extend outward

				if (check_stacked) {
					if (parts[1].d[dim][dir] != part.d[dim][dir]) continue; // walls not aligned - may clip through roof overhang, skip
					// if entire upper wall is not shared, then this is a bad placement because the upper windows are likely misaligned
					if (parts[1].d[!dim][0] != part.d[!dim][0] || parts[1].d[!dim][1] != part.d[!dim][1]) continue;
					bc.z2() = parts[1].z2(); // extend upward
				}
				cube_t bc_pad(bc);
				bc_pad.d[dim][dir] += (dir ? 1.0 : -1.0)*2.0*get_scaled_player_radius(); // add space for the player; may not be needed
				cube_t bc_exp(bc_pad);
				bc_exp.expand_by(wall_thickness); // needed for exterior door check
				if (has_bcube_int_no_adj(bc_pad, parts))          continue; // check for intersection with other parts, in particular the chimney and fireplace
				if (has_driveway() && bc.intersects_xy(driveway)) continue; // skip if intersects driveway or garage
				if (cube_int_ext_door(bc_exp))                    continue; // check exterior doors
				if (has_bcube_int(bc_exp, details))               continue; // check details; outdoor AC units can intersect
				interior->room_geom->objs.emplace_back(bc, TYPE_LADDER, 0, dim, dir, RO_FLAG_EXTERIOR, 1.0, SHAPE_CUBE, GRAY); // room_id=0
				union_with_coll_bcube(bc);
				ladder = bc;
				return; // success/done
			} // for n
		} // for d
	}
}

void building_t::add_balconies(rand_gen_t &rgen, vect_cube_t &balconies) {
	if (!is_house || !has_room_geom()) return; // houses only for now
	if (rgen.rand_bool()) return; // only add balconies to 50% of houses
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), door_width(get_doorway_width());
	float const min_depth(max(0.5f*floor_spacing, 2.5f*get_scaled_player_radius()/(1.0f - BALCONY_PILLAR_SCALE))); // make sure the player can fit around the pillars
	float const balcony_depth(min_depth*rgen.rand_uniform(1.0, 1.2)); // constant per house
	float const room_min_z2(ground_floor_z1 + 1.5*floor_spacing); // > 1 floor
	unsigned const balcony_style(rgen.rand()); // shared across all balconies of this house
	unsigned const max_balconies(1 + rgen.rand_bool()); // 1-2 per house
	unsigned num_balconies(0);
	auto &objs(interior->room_geom->objs);
	unsigned const init_objs_sz(objs.size());
	vect_cube_t avoid;
	if (!objs.empty() && (objs.back().type == TYPE_FESCAPE || objs.back().type == TYPE_LADDER)) {avoid.push_back(objs.back());} // avoid fire escape or ladder
	if (has_driveway()) {avoid.push_back(driveway);}
	if (!city_driveway.is_all_zeros()) {avoid.push_back(city_driveway);}

	// find suitable rooms for balconies; since room walls will never intersect windows, we can make the balcony the same width to avoid intersecting windows
	for (auto room = interior->rooms.begin(); room != interior->rooms.end(); ++room) {
		if (room->interior)           continue; // no windows
		if (room->is_sec_bldg)        continue; // no garage or shed, even if it's multiple stories tall
		if (room->is_single_floor)    continue; // no single floor tall rooms
		if (room->z2() < room_min_z2) continue; // ground floor only
		if (rgen.rand_float() < 0.75) continue;
		float const balcony_z1(room->z2() - floor_spacing /*+ get_fc_thickness()*/); // floor level of top floor of room
		unsigned const floor_ix(room->get_floor_containing_zval((balcony_z1 + 0.1*floor_spacing), floor_spacing));
		unsigned const room_id(room - interior->rooms.begin());
		room_type const rtype(room->get_room_type(floor_ix));
		if (is_bathroom(rtype)) continue; // no bathroom balconies as that would be weird
		assert(room->part_id < parts.size());
		cube_t const &part(parts[room->part_id]);
		bool added(0);

		for (unsigned dim = 0; dim < 2 && !added; ++dim) {
			for (unsigned dir = 0; dir < 2 && !added; ++dir) {
				if (classify_room_wall(*room, balcony_z1, dim, dir, 1) != ROOM_WALL_EXT) continue; // not fully exterior wall, skip
				float const ext_wall_pos(room->d[dim][dir]);
				cube_t balcony(*room);
				balcony.z1() = balcony_z1;
				balcony.d[dim][!dir]  = ext_wall_pos; // abuts the exterior wall of the room
				balcony.d[dim][ dir] += (dir ? 1.0 : -1.0)*balcony_depth; // extend outward from the house
				if (has_bcube_int(balcony, avoid)) continue; // blocked
				if (check_cube_intersect_non_main_part(balcony)) continue; // porch roof, porch support, and chimney, etc.
				bool bad_pos(0);

				if (room->has_stairs_on_floor(floor_ix)) { // check for stairs blocking balcony access
					cube_t test_cube(balcony);
					test_cube.d[dim][!dir] -= (dir ? 1.0 : -1.0)*door_width; // move into the interior
					set_wall_width(test_cube, balcony.get_center_dim(!dim), wall_thickness, !dim); // center of the balcony in long dim

					for (stairwell_t const &s : interior->stairwells) {
						if (s.intersects(test_cube)) {bad_pos = 1; break;}
					}
					if (bad_pos) continue;
				}
				if (real_num_parts > 1) { // check if adjacent to the second part (may block a window or clip through a windowsill)
					cube_t const &other_part(parts[(room->part_id == 0) ? 1 : 0]);
					cube_t test_cube(balcony);
					test_cube.z1() = ground_floor_z1; // extend down to the ground
					test_cube.expand_in_dim( dim,    -wall_thickness); // shrink to exclude adjacencies
					test_cube.expand_in_dim(!dim, 2.0*wall_thickness); // expand to include adjacencies
					if (test_cube.intersects(other_part)) continue;
				}
				if (!exterior_flag.is_all_zeros() && exterior_flag.intersects(balcony)) continue; // flag placed in the way, no balcony
				cube_t balcony_ext_down(balcony), balcony_ext_out(balcony);
				balcony_ext_down.z1() = ground_floor_z1; // extend down to the ground

				for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // check for any intersecting parts (likely below this one, for stacked parts)
					if ((p - parts.begin()) != room->part_id && p->intersects_no_adj(balcony_ext_down)) {bad_pos = 1; break;}
				}
				if (bad_pos) continue;
				balcony_ext_out.expand_by_xy(0.2*floor_spacing); // for testing against upper floor doors
				if (cube_int_ext_door(balcony_ext_out)) continue; // check exterior doors; we may want to create doors to balconies in the future
				if (world_mode == WMODE_INF_TERRAIN && cube_int_tiled_terrain_trees(balcony - get_tiled_terrain_model_xlate())) continue; // check trees (slow)
				balcony.z2() -= 0.6*floor_spacing; // reduce wall height to 40%
				balcony.expand_in_dim(!dim, wall_thickness); // expand slightly to include window frame and merge adj balcony shared walls
				max_eq(balcony.d[!dim][0], (part.d[!dim][0] + 0.25f*wall_thickness)); // clamp slightly smaller than the containing part in !dim
				min_eq(balcony.d[!dim][1], (part.d[!dim][1] - 0.25f*wall_thickness));
				// if the space below the balcony is blocked by something, flag as hanging so that we don't try to draw vertical supports later
				cube_t area_below(balcony);
				set_cube_zvals(area_below, ground_floor_z1, balcony.z1());
				bool const hanging(has_bcube_int(area_below, avoid) || check_cube_intersect_non_main_part(area_below));
				room_object_t balcony_obj(balcony, TYPE_BALCONY, room_id, dim, dir, (hanging ? RO_FLAG_HANGING : 0), 1.0, SHAPE_CUBE, WHITE);
				balcony_obj.obj_id = balcony_style; // set so that we can select from multiple balcony styles
				unsigned const draw_style(balcony_style & 3); // wooden walls, metal railing + bars, metal railing + 45 deg rotated bars, metal railing + wood sides
				unsigned const balcony_obj_ix(objs.size());
				objs.push_back(balcony_obj);
				balconies.push_back(balcony);
				avoid.push_back(balcony); // shouldn't need to consider area_below
				avoid.back().expand_by_xy(wall_thickness);
				// add exterior step for this balcony so that the player can stand on it and can't pass through the railings
				cube_t floor_slab(balcony);
				floor_slab.z2() = balcony.z1() + 0.12*balcony.dz(); // matches code in get_balcony_cubes()
				ext_steps.emplace_back(floor_slab, dim, 0, dir, 0, 0, 0, 1); // enclosed, no step dir
				details.emplace_back(floor_slab, DETAIL_OBJ_SHAD_ONLY);

				if (draw_style == 0 || draw_style == 3) { // add shadow casters for sides
					cube_t cubes[4]; // {bottom, front, left side, right side}
					get_balcony_cubes(balcony_obj, cubes);
					for (unsigned n = 1; n < 4; ++n) {details.emplace_back(cubes[n], DETAIL_OBJ_SHAD_ONLY);} // skip bottom, which was added above
				}
				// add door connecting to the house if possible
				if (has_windows()) { // find a space between two windows
					float const door_width(get_doorway_width()), door_hwidth(0.5*door_width), edge_pad(2.0*wall_thickness);
					float const window_hspacing(get_hspacing_for_part(part, !dim)), window_h_border(get_window_h_border());

					if (door_width > 0.0 && door_width < window_hspacing*2.0*window_h_border) { // door can fit between two windows
						unsigned const num_windows(round_fp(part.get_sz_dim(!dim)/window_hspacing));
						vect_cube_t cands;

						for (unsigned n = 0; n+1 < num_windows; ++n) { // iterate over each window pair (skip last window)
							float const pos(part.d[!dim][0] + n*window_hspacing); // midpoint between two windows
							float const lo(pos - door_hwidth), hi(pos + door_hwidth);
							if (lo < balcony.d[!dim][0]+edge_pad || hi > balcony.d[!dim][1]-edge_pad) continue; // not fully within the balcony
							float const door_z1(balcony_z1 + get_fc_thickness());
							cube_t door;
							set_cube_zvals(door, door_z1, (door_z1 + get_door_height()));
							set_wall_width(door, pos, door_hwidth, !dim);
							set_wall_width(door, part.d[dim][dir], 2.0*get_trim_thickness(), dim);
							cube_t door_exp(door);
							door_exp.d[dim][!dir] += (dir ? -1.0 : 1.0)*door_width; // add clearance for door to open
							if (interior->is_blocked_by_stairs_or_elevator(door_exp)) continue;
							bool blocked(0);

							for (auto i = objs.begin(); i != objs.begin()+init_objs_sz; ++i) { // check if blocked by furniture or something on the other side of the wall
								if (i->type == TYPE_BOX || i->type == TYPE_LG_BALL || i->type == TYPE_TCAN) continue; // easy to move, so can block the door
								if (i->type == TYPE_RUG || i->type == TYPE_FLOORING) continue; // not blocking
								
								if (i->intersects(door_exp)) {
									if (i->type == TYPE_OUTLET) {i->remove();} // if blocked by an outlet, just remove the outlet (assuming this cand is chosen)
									blocked = 1;
									break;
								}
							}
							if (!blocked) {cands.push_back(door);}
						} // for n
						if (!cands.empty()) {
							cube_t const &door(cands[rgen.rand() % cands.size()]);
							objs.emplace_back(door, TYPE_FALSE_DOOR, room_id, dim, dir, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, door_color);
						}
					}
				}
				// maybe add plants to balconies; note that they won't be properly lit since plants use indoor lighting,
				// and the plants won't be drawn when the player is outside the building; this should be okay because an empty pot works on a balcony as well
				unsigned const num_plants(rgen.rand() % 3); // 0-2

				if (num_plants > 0) { // place plants
					cube_t cubes[4];
					vect_cube_t plant_avoid;
					get_balcony_cubes(objs.back(), cubes);
					float const wall_width(cubes[1].get_sz_dim(dim)); // use front wall
					cube_t floor_inner(cubes[0]);
					floor_inner.expand_by_xy(-wall_width);

					for (unsigned n = 0; n < num_plants; ++n) {
						if (place_plant_on_obj(rgen, floor_inner, room_id, 1.0, 1.2, plant_avoid)) {plant_avoid.push_back(objs.back());} // tot_light_amt=1.0, sz_scale=1.2
					}
				}
				// else place table and chairs?
				if (!hanging) { // add colliders for vertical supports
					cube_t pillar[2];
					get_balcony_pillars(balcony_obj, ground_floor_z1, pillar);
					bool no_pillars(0);

					for (unsigned d = 0; d < 2; ++d) { // check for pillars blocking exterior door or clipping through the fence
						cube_t pillar_exp(pillar[d]);
						pillar_exp.d[dim][!dir] = ext_wall_pos; // extend to exterior wall
						pillar_exp.expand_by_xy(wall_thickness);
						if (cube_int_ext_door(pillar_exp) || has_bcube_int(pillar[d], fences)) {no_pillars = 1; break;}
					}
					if (no_pillars) {objs[balcony_obj_ix].flags |= RO_FLAG_HANGING;} // make it hanging instead
					else { // add pillars
						for (unsigned d = 0; d < 2; ++d) {details.emplace_back(pillar[d], DETAIL_OBJ_COLL_SHAD);} // collider + shadow caster
					}
				}
				union_with_coll_bcube(balcony_ext_down);
				++num_balconies;
				added = 1;
			} // for dir
		} // for dim
		if (num_balconies == max_balconies) break; // done
	} // for room
	if (num_balconies > 0) {invalidate_tile_smap_in_region(bcube + get_camera_coord_space_xlate());}
}

void building_t::add_gutter_downspouts(rand_gen_t &rgen, vect_cube_t const &balconies) {
	for (cube_with_ix_t const &g : gutters) {
		bool const dim(g.ix >> 1), dir(g.ix & 1);
		float const len(g.get_sz_dim(!dim)), width(g.get_sz_dim(dim)), ds_width(0.5f*width), edge_spacing(1.65*ds_width); // just enough to clear the fence post
		float const wall_pos(g.d[dim][!dir]), dir_sign(dir ? 1.0 : -1.0);
		assert(len > 0.0 && width > 0.0);
		assert(ground_floor_z1 < g.z1());
		cube_t ds;
		set_cube_zvals(ds, ground_floor_z1, g.z1());
		ds.d[dim][!dir] = wall_pos; // at the wall
		ds.d[dim][ dir] = wall_pos + dir_sign*0.6*ds_width; // extend out from the wall
		// find part associated with this gutter and clip interior gutter to this
		cube_t int_gutter(g), part;
		point query_pt;
		query_pt[ dim] = wall_pos - dir_sign*0.1*ds_width; // slightly inside the wall
		query_pt[!dim] = g.get_center_dim(!dim);
		query_pt.z = g.z1() - width; // slightly down
		bool skip_ends[2] = {0,0};

		for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
			if (!p->contains_pt(query_pt)) continue;
			// gutters usually extend beyond the building part, but they can end inside it for two-part house roofs meeting at an inside corner where the roof tquad is clipped;
			// in this case we skip the downspout because it may be too close to a lower window
			if (g.d[!dim][0] > p->d[!dim][0]) {skip_ends[0] = 1;}
			if (g.d[!dim][1] < p->d[!dim][1]) {skip_ends[1] = 1;}
			max_eq(int_gutter.d[!dim][0], p->d[!dim][0]); // clip !dim range to the building part
			min_eq(int_gutter.d[!dim][1], p->d[!dim][1]);
			assert(part.is_all_zeros());
			part  = *p; // part may be useful to have below
		}
		assert(!part.is_all_zeros()); // must be found

		for (unsigned e = 0; e < 2; ++e) { // add gutter at each end
			if (skip_ends[e]) continue;
			float const centerline(int_gutter.d[!dim][e] + (e ? -1.0 : 1.0)*edge_spacing);
			set_wall_width(ds, centerline, 0.5*ds_width, !dim);
			if (has_bcube_int(ds, balconies)) continue; // check balconies; is this necessary?
			// what about fire escapes? is it possible for the gutter to intersect them, and does that look wrong?
			if (has_chimney == 2 && (get_chimney().intersects(ds) || get_fireplace().intersects(ds))) continue; // check exterior chimney
			if (has_porch() && porch.intersects_xy(ds)) continue; // check porch (roof)
			cube_t ds_exp(ds);
			ds_exp.d[dim][!dir] -= dir_sign*get_wall_thickness(); // move toward the house so that the garage door is intersected
			
			if (has_sec_bldg() && get_sec_bldg().intersects(ds_exp)) { // gutter on garage or shed
				bool const inner_dir(centerline < bcube.get_center_dim(!dim));
				if (bool(e) == inner_dir) continue; // only add gutter on the ouyside facing direction
			}
			else if (real_num_parts > 1) { // check for intersections with lower parts when stacked
				cube_t test_cube(ds);
				test_cube.d[dim][!dir] += dir_sign*0.1*ds_width; // move away from the wall to prevent self-intersection with upper part
				if (cube_int_parts_no_sec(test_cube)) continue;
				// skip gutters from upper stacked parts down through lower parts as they may intersect a window
				if (part.z1() > ground_floor_z1 && g.d[!dim][e] > bcube.d[!dim][0] && g.d[!dim][e] < bcube.d[!dim][1]) continue;
			}
			// check for intersections with garage doors and corner front doors
			ds_exp.expand_in_dim(!dim, 4.0*ds_width); // add extra padding to avoid doorbells and lamps
			bool door_int(0);

			for (tquad_with_ix_t const &door : doors) {
				if (door.get_bcube().intersects(ds_exp)) {door_int = 1; break;}
			}
			if (door_int) continue;
			interior->room_geom->objs.emplace_back(ds, TYPE_DOWNSPOUT, 0, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR), 1.0, SHAPE_CUBE, WHITE);
		} // for e
	} // for g
}

void building_t::add_exterior_ac_pipes(rand_gen_t rgen) {
	if (!has_ac) return; // no AC
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR);
	colorRGBA const colors[3] = {BLACK, COPPER_C, GRAY }; // insulated, non-insulated, power
	float    const radius [3] = {0.06,  0.025,    0.035}; // relative to height
	float    const offsets[3] = {0.2,   0.32,     0.85 }; // from end
	vect_cube_t avoid_pipes;

	for (roof_obj_t const &ac : details) {
		if (ac.type != ROOF_OBJ_AC) continue;
		bool const min_dim(ac.dy() < ac.dx()); // adjacent to wall in short dim/along long dim (for ground AC units)
		float const depth(ac.get_sz_dim(min_dim)); // actual spacing is smaller than this
		float const height(ac.dz()), length(ac.get_sz_dim(!min_dim));
		
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (ac.z1() == ground_floor_z1) { // AC on ground (house); find closest wall of closest part
				cube_t bc(*p);
				bc.expand_in_dim(min_dim, 2.0*depth);
				if (!bc.contains_cube(ac)) continue; // wrong part
				bool const wall_dir(p->get_center_dim(min_dim) < ac.get_center_dim(min_dim));
				cube_t gap(ac);
				gap.d[min_dim][ wall_dir] = ac.d[min_dim][!wall_dir]; // edge of AC unit
				gap.d[min_dim][!wall_dir] = p->d[min_dim][ wall_dir]; // edge of part exterior wall
				cube_t pipe(gap);

				for (unsigned n = 0; n < 3; ++n) {
					float const r(radius[n]*height);
					set_wall_width(pipe, (ac.z1() + 0.12*height + r), r, 2);
					set_wall_width(pipe, (ac.d[!min_dim][0] + offsets[n]*length), r, !min_dim);
					interior->room_geom->objs.emplace_back(pipe, TYPE_PIPE, 0, min_dim, 0, (flags | RO_FLAG_HANGING), 1.0, SHAPE_CYLIN, colors[n]); // flat ends
				}
			}
			else { // AC on roof (office); find part whose roof contains the unit
				if (p->z2() != ac.z1() || !p->contains_cube_xy(ac)) continue; // wrong part
				bool const pref_dir(rgen.rand_bool());
				float const pipe_ext_z1(p->z1() - get_fc_thickness()), pipe_len(max(4.0f*radius[0]*height, rgen.rand_uniform(0.25, 1.0)*depth));
				cube_t valid_area(*p);
				valid_area.expand_by_xy(-get_wall_thickness()); // shrink to account for any roof wall

				for (unsigned d = 0; d < 2; ++d) { // try both dirs
					bool const dir(pref_dir ^ bool(d));
					float const edge_val(ac.d[min_dim][dir]), dsign(dir ? 1.0 : -1.0);
					cube_t pipe(ac);
					pipe.d[min_dim][!dir] = edge_val; // edge of AC unit
					pipe.d[min_dim][ dir] = edge_val + dsign*pipe_len; // edge of part exterior wall
					pipe.z1() = pipe_ext_z1; // extend downward to include helipads, etc.; zval will be overwritten below
					if (!valid_area.contains_cube_xy(pipe))            continue;
					if (has_bcube_int(pipe, skylights))                continue;
					if (d == 0 && has_bcube_int_no_adj(pipe, details)) continue; // if pref dir is blocked, try the other dir

					for (unsigned n = 0; n < 3; ++n) { // 3 pipes
						float const r(radius[n]*height), pipe_end(pipe.d[min_dim][dir]), vert_ext(pipe_end + dsign*r);
						if (length < 5.0*r) continue; // too small for pipe - shouldn't happen
						set_wall_width(pipe, (ac.z1() + 0.24*height + r), r, 2);
						set_wall_width(pipe, (ac.d[!min_dim][0] + offsets[n]*length), r, !min_dim);

						for (unsigned m = 0; m < 4; ++m) { // 4 placement tries
							if (m > 0) { // select a random position
								set_wall_width(pipe, rgen.rand_uniform(ac.d[!min_dim][0]+2.0*r, ac.d[!min_dim][1]-2.0*r), r, !min_dim);
							}
							cube_t pipe_ext(pipe);
							pipe_ext.z1() = pipe_ext_z1; // extend downward to include helipads, etc.
							pipe_ext.d[min_dim][dir] = vert_ext; // space for the bend + vertical section
							if (has_bcube_int_no_adj(pipe_ext, details)) continue; // no adj to exclude the AC unit the pipe is connected to
							if (has_bcube_int(pipe_ext, avoid_pipes))    continue;
							avoid_pipes.push_back(pipe_ext);
							interior->room_geom->objs.emplace_back(pipe, TYPE_PIPE, 0, min_dim, 0, (flags | RO_FLAG_HANGING), 1.0, SHAPE_CYLIN, colors[n]); // flat ends
							// add bend and vertical section
							cube_t vpipe(pipe);
							vpipe.d[min_dim][!dir] = pipe_end - dsign*r;
							vpipe.d[min_dim][ dir] = vert_ext;
							vpipe.z2() -= r;
							vpipe.z1()  = p->z2();
							interior->room_geom->objs.emplace_back(vpipe, TYPE_PIPE, 0, 0, 1, (flags | RO_FLAG_ADJ_HI), 1.0, SHAPE_CYLIN, colors[n]); // round top end
							break;
						} // for m
					} // for n
					break;
				} // for d
			}
			break; // done - there should only be one part
		} // for p
	} // for i
}

void building_t::add_padlocks(rand_gen_t rgen) {
	if (!is_house) return; // houses only for now, to avoid adding too many locks
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_PADLOCK)) return; // no model
	assert(has_room_geom());
	// determine which color of keys we have hidden in drawers, so that we don't place a padlock that can't be opened;
	// this doesn't guarantee we can open it because the key may be behind the door, but it's a good check anyway
	unsigned key_color_mask(0);
	room_object_t drawers_part;
	vect_cube_t drawers;
	float drawer_extend(0.0); // unused

	for (room_object_t const &obj : interior->room_geom->objs) { // iterate over all objects; buttons and stairs haven't been placed yet
		if (obj.type != TYPE_DRESSER && obj.type != TYPE_NIGHTSTAND && (obj.type != TYPE_DESK || !obj.desk_has_drawers())) continue; // item doesn't have drawers
		room_object_t obj_drawers_open(obj);
		obj_drawers_open.drawer_flags = ~uint16_t(0); // make all drawers open so that we get the correct drawers bounds
		get_obj_drawers_or_doors(obj_drawers_open, drawers, drawers_part, drawer_extend);

		for (unsigned dix = 0; dix < drawers.size(); ++dix) {
			float stack_z1(0.0);

			for (unsigned item_ix = 0; item_ix < 16; ++item_ix) { // take the *last* item in the drawer first, which will be the top item if stacked
				room_object_t const item(interior->room_geom->get_item_in_drawer(drawers_part, drawers[dix], dix, item_ix, stack_z1));
				if (item.type == TYPE_NONE) break; // no more items
				if (item.type != TYPE_KEY ) continue;
				assert(item.obj_id < NUM_LOCK_COLORS);
				key_color_mask |= (1 << item.obj_id);
			}
		} // for dix
	} // for i
	if (key_color_mask == 0) return; // no keys, so no padlocks

	for (auto d = interior->doors.begin(); d != interior->doors.end(); ++d) {
		if (d->open || !d->locked || rgen.rand_bool()) continue;
		add_padlock_to_door((d - interior->doors.begin()), key_color_mask, rgen);
	}
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
	if (!has_room_geom() || (interior->walls[0].empty() && interior->walls[1].empty())) return; // no interior or walls
	if (interior->room_geom->trim_was_added) return; // trim already generated
	interior->room_geom->trim_was_added = 1;
	add_wall_and_door_trim();
	interior->room_geom->trim_objs.shrink_to_fit();
}

void cut_trim_around_doors(vect_tquad_with_ix_t const &doors, vect_cube_t &trim_cubes, float door_expand, bool dim) {
	for (auto d = doors.begin(); d != doors.end(); ++d) {
		cube_t door(d->get_bcube());
		bool const door_dim(door.dy() < door.dx());
		if (door_dim != bool(dim)) continue;
		door.expand_in_dim(door_dim, door_expand); // expand to nonzero area; use a larger expand to account for distance door is offset away from ext wall
		subtract_cube_from_cubes(door, trim_cubes); // subtract this door from current trim cubes by clipping in XY
	}
}
void clip_trim_cube(cube_t const &trim, cube_t const &trim_exclude, vect_cube_t &trim_cubes) {
	trim_cubes.clear();
	if (!trim_exclude.is_all_zeros() && trim_exclude.intersects(trim)) {subtract_cube_from_cube(trim, trim_exclude, trim_cubes);}
	else {trim_cubes.push_back(trim);}
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
	vect_cube_t trim_cubes, trim_parts;
	cube_t trim_exclude;

	// exclude trim at intermediate floors of tall rooms
	for (room_t const &room : interior->rooms) {
		if (room.is_single_floor && room.dz() > 1.5*window_vspacing) {
			trim_exclude = room;
			trim_exclude.expand_by_xy(0.5*wall_thickness); // include half the wall
			trim_exclude.expand_in_dim(2, -0.5*window_vspacing); // allow trim at floor and ceiling, but not at floors in between
			break; // can only have one room of this type
		}
	}
	for (auto d = interior->door_stacks.begin(); d != interior->door_stacks.end(); ++d) { // vertical strips on each side + strip on top of interior doors
		if (d->on_stairs) continue; // no frame for stairs door, skip
		cube_t trim(*d);
		set_wall_width(trim, trim.get_center_dim(d->dim), door_trim_exp, d->dim);
		trim.z2() -= 0.1*trim_toler; // shift top down oh so slightly to prevent z-fighting with top of wall when drawn under a skylight

		for (unsigned side = 0; side < 2; ++side) { // left/right of door
			trim.d[!d->dim][0] = d->d[!d->dim][side] - (side ? trim_thickness : door_trim_width);
			trim.d[!d->dim][1] = d->d[!d->dim][side] + (side ? door_trim_width : trim_thickness);
			bool const draw_top(d->mult_floor_room || check_skylight_intersection(trim)); // draw top edge of trim for top floor?
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
			bool const draw_top(d->mult_floor_room || (f+1 == num_floors && check_skylight_intersection(trim))); // draw top edge of trim for top floor?
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
		float const trim_width(garage_door ? 0.016*door.get_sz_dim(!dim) : door_trim_width); // garage door trim is based on width

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
		bool const draw_bot(door.z1() > ground_floor_z1 + floor_thickness); // door is above the ground floor, draw the bottom edge
		trim.d[!dim][0] = door.d[!dim][0];
		trim.d[!dim][1] = door.d[!dim][1];
		set_cube_zvals(trim, door.z1(), (door.z1() + 0.1*fc_thick + trim_thickness)); // floor height + extend slightly above
		objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, (ext_flags | (draw_bot ? 0 : RO_FLAG_ADJ_BOT)), 1.0, SHAPE_SHORT, ext_trim_color);

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
			bool const in_basement(w->zc() < ground_floor_z1);
			if (in_basement && interior->has_backrooms && !get_basement().intersects_no_adj(*w)) continue; // no trim in basement backrooms
			cube_t trim(*w);
			trim.expand_in_dim(dim, trim_thickness);

			// handle outside corners of office building hallway intersections; not needed for basements
			if (has_outside_corners && !in_basement) {
				auto const end(interior->walls[!dim].begin()+interior->extb_walls_start[!dim]); // exclude extended basement walls

				for (auto W = interior->walls[!dim].begin(); W != end; ++W) { // check walls in other dim for an outside corner
					for (unsigned d = 0; d < 2; ++d) { // check both ends of the current wall/trim
						if (W->z1() > w->z2() || W->z2() < w->z1()) continue; // no z overlap, wrong stack
						if (W->d[!dim][0] > w->d[!dim][d]+trim_toler || W->d[!dim][1] < w->d[!dim][d]-trim_toler) continue; // not adjacent/overlapping
						if (W->d[ dim][0] > w->d[ dim][1]+trim_toler || W->d[ dim][1] < w->d[ dim][0]-trim_toler) continue; // not adjacent/overlapping in other dim
						if (W->d[ dim][0] < w->d[ dim][0]-trim_toler && W->d[ dim][1] > w->d[ dim][1]+trim_toler) continue; // skip T junctions
						trim.d[!dim][d] = W->d[!dim][d] + (d ? 1.0 : -1.0)*trim_thickness; // expand to cover gap at outside corners of hallway walls
					}
				} // for W
			}
			if (w->dz() < 0.5*window_vspacing) continue; // short wall segment from tall room extension, no trim
			unsigned const num_floors(calc_num_floors(*w, window_vspacing, floor_thickness));
			// snap to the nearest floor to handle short walls due to cut out stairs
			float const ground_wall_z1(bcube.z1() + fc_thick);
			float z(ground_wall_z1 + window_vspacing*round_fp((w->z1() - ground_wall_z1)/window_vspacing));

			for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
				set_cube_zvals(trim, z, z+trim_height); // starts at floor height
				bool ext_dirs[2] = {0,0};

				if (z < ground_floor_z1 && has_ext_basement() && !get_basement().intersects(trim)) { // check for exterior wall of extended basement
					for (unsigned dir = 0; dir < 2; ++dir) { // for each side of wall
						point test_pt(0.0, 0.0, trim.z2());
						test_pt[dim] = w->d[dim][dir];
						float const test_vals[3] = {w->d[!dim][0], w->get_center_dim(!dim), w->d[!dim][1]}; // try both ends and center point
						bool any_inside(0);
						
						for (unsigned d = 0; d < 3 && !any_inside; ++d) {
							test_pt[!dim] = test_vals[d];
							any_inside   |= interior->point_in_ext_basement_room(test_pt);
						}
						ext_dirs[dir] = !any_inside;
					}
				}
				unsigned const trim_flags(flags | (ext_dirs[0] ? RO_FLAG_ADJ_LO : 0) | (ext_dirs[1] ? RO_FLAG_ADJ_HI : 0)); // disable exterior faces
				clip_trim_cube(trim, trim_exclude, trim_parts);
				for (cube_t const &t : trim_parts) {objs.emplace_back(t, TYPE_WALL_TRIM, 0, dim, 0, trim_flags, 1.0, SHAPE_CUBE, trim_color);} // floor trim
				if (!has_ceil_trim) continue;
				trim.z2() = z + floor_to_ceil_height; // ceil height
				trim.z1() = trim.z2() - trim_height;

				for (unsigned dir = 0; dir < 2; ++dir) { // for each side of wall
					if (ext_dirs[dir]) continue; // skip
					cube_t ceil_trim(trim);
					ceil_trim.d[dim][!dir] = w->d[dim][dir];
					clip_trim_cube(ceil_trim, trim_exclude, trim_parts);
					for (cube_t const &t : trim_parts) {objs.emplace_back(t, TYPE_WALL_TRIM, 0, dim, dir, flags, 1.0, SHAPE_ANGLED, trim_color);} // ceiling trim
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
						if (has_ext_door_this_floor(i->z1(), f)) {cut_trim_around_doors(doors, trim_cubes, (expand_val + wall_thickness), d);} // cut out areas for ext doors
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
							cube_t trim2(*c); // copy so that we can modify it
							set_cube_zvals(trim2, ceil_trim_z1, ceil_trim_z2);
							clip_trim_cube(trim2, trim_exclude, trim_parts);
							for (cube_t const &t : trim_parts) {objs.emplace_back(t, TYPE_WALL_TRIM, 0, dim, !dir, flags, 1.0, SHAPE_ANGLED, trim_color);} // ceiling trim
						}
					}
					if (has_ext_door_this_floor(i->z1(), f)) { // cut out areas for ext doors; not for stacked parts
						cut_trim_around_doors(doors, trim_cubes, (expand_val + wall_thickness), dim);
					}
					for (cube_t &c : trim_cubes) {
						clip_trim_cube(c, trim_exclude, trim_parts);
						for (cube_t const &t : trim_parts) {objs.emplace_back(t, TYPE_WALL_TRIM, 0, dim, 0, ext_flags, 1.0, SHAPE_CUBE, trim_color);} // floor trim
						if (!has_ceil_trim || is_house) continue;
						set_cube_zvals(c, ceil_trim_z1, ceil_trim_z2); // okay to edit in-place here
						clip_trim_cube(c, trim_exclude, trim_parts);
						for (cube_t const &t : trim_parts) {objs.emplace_back(t, TYPE_WALL_TRIM, 0, dim, !dir, flags, 1.0, SHAPE_ANGLED, trim_color);} // ceiling trim
					}
				} // for f
			} // for dir
		} // for dim
	} // for i
	// add trim on door thresholds for bathroom flooring and pool tile
	auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

	for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
		if (i->type == TYPE_FLOORING && (i->item_flags == FLOORING_TILE || i->item_flags == FLOORING_MARBLE || i->item_flags == FLOORING_CONCRETE)) {} // bath/server room flooring
		else if (is_pool_tile_floor(*i)) {} // pool tile
		else {continue;} // no trim
		cube_t flooring_exp(*i);
		flooring_exp.expand_by_xy(0.5*wall_thickness);
		flooring_exp.z2() += 0.5*i->dz(); // slightly taller

		for (door_stack_t const &ds : interior->door_stacks) {
			if (ds.on_stairs || ds.for_closet) continue; // skip basement and closet doors
			cube_t door_bc(ds.get_true_bcube());
			if (!i->intersects(door_bc))       continue;
			door_bc.expand_in_dim(ds.dim, -0.05*door_bc.get_sz_dim(ds.dim)); // shrink slightly to prevent Z-fighting
			assert(i->intersects(door_bc));
			cube_t trim(flooring_exp);
			trim.intersect_with_cube(door_bc);
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, ds.dim, 0, flags, 1.0, SHAPE_CUBE, trim_color);
		} // for ds
	} // for i
	if (!is_cube() || is_rotated()) return; // window trim is not yet working for non-cube and rotated buildings
	add_window_trim_and_coverings(1, 0, 0); // add_trim=1, add_coverings=0, add_ext_sills=0
}

void building_t::add_window_trim_and_coverings(bool add_trim, bool add_coverings, bool add_ext_sills) {
	assert(add_trim || add_coverings || add_ext_sills);
	if (!has_windows()) return; // no windows
	float const trim_thickness(get_trim_thickness()), ext_wall_toler(0.01*trim_thickness); // required to prevent z-fighting when AA is disabled
	float const window_h_border(WINDOW_BORDER_MULT*get_window_h_border()), window_v_border(WINDOW_BORDER_MULT*get_window_v_border()); // (0, 1) range
	// Note: depth must be small to avoid object intersections; this applies to the windowsill as well
	float const window_trim_width(0.75*get_wall_thickness()), window_trim_depth(1.0*trim_thickness), windowsill_depth(1.0*trim_thickness);
	float const floor_spacing(get_window_vspace()), window_offset(get_door_shift_dist());
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
		bool const is_attic(has_attic() && c.z1() >= get_attic_part().z2());
		cube_t window(c); // copy dim <dim>
		window.translate_dim(dim, dscale*window_offset);
		window.d[dim][!dir] += dscale*(is_attic ? (get_attic_beam_depth() + 0.5*window_trim_depth) : window_trim_depth); // add thickness on interior (beam depth for attic)
		window.d[dim][ dir] += dscale*ext_wall_toler; // slight bias away from the exterior wall
		unsigned const ext_flags(RO_FLAG_NOCOLL | (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));

		for (float z = tz1; z < tz2; z += 1.0) { // each floor
			float const bot_edge(c.z1() + (z - tz1)*window_height);
			set_cube_zvals(window, bot_edge+border_z, bot_edge+window_height-border_z);
			// add separators for 50% of walls/floors for houses; not for attic windows; not for glass block windows?
			bool const add_separators(is_house && !is_attic && rgen.rand_bool());
			bool const one_dim_separators(add_separators && rgen.rand_bool()); // 1=vert/horiz separators only, 0=cross shaped separators
			bool const sep_dim((window_height - 2.0*border_z) < 0.9*(window_width - 2.0*border_xy)); // split in larger-ish dim

			for (float xy = tx1; xy < tx2; xy += 1.0) { // windows along each wall
				float const low_edge(c.d[!dim][0] + (xy - tx1)*window_width);
				window.d[!dim][0] = low_edge + border_xy;
				window.d[!dim][1] = low_edge + window_width - border_xy;
				if (add_coverings && !is_attic) {add_window_coverings(window, dim, dir);}

				if (add_ext_sills && !is_attic) {
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

					if (is_bathroom(get_room_type_and_floor(room_id, window.zc(), floor_ix))) { // check for bathroom block windows
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

cube_t building_t::get_step_for_ext_door(tquad_with_ix_t const &door) const {
	cube_t const c(door.get_bcube());
	bool const dim(c.dy() < c.dx()), dir(door.get_norm()[dim] > 0.0);
	float length(((door.type == tquad_with_ix_t::TYPE_GDOOR) ? 0.6 : 0.5)*c.dz());
	max_eq(length, 2.4f*get_scaled_player_radius()); // make sure step is wide enough for the player to walk on
	cube_t step(c);
	set_cube_zvals(step, (c.z1() - get_fc_thickness()), c.z1());
	step.d[dim][ dir] += (dir ? 1.0 : -1.0)*length; // extend outward
	return step;
}
void building_t::add_ext_door_steps(unsigned ext_objs_start) {
	float const floor_spacing(get_window_vspace()), fc_thickness(get_fc_thickness());
	float const door_shift_dist(2.5*get_door_shift_dist()); // 1x for door shift and 1.5x offset in add_door()
	colorRGBA const step_color(LT_GRAY);
	vect_room_object_t &objs(interior->room_geom->objs);
	vector<unsigned> to_add_stairs;

	// add step at the base of each exterior door
	for (auto const &d : doors) {
		if (d.type == tquad_with_ix_t::TYPE_RDOOR) continue; // skip roof access door
		cube_t const c(d.get_bcube());
		bool const above_ground(c.z1() > ground_floor_z1 + 2.0*fc_thickness);
		bool const dim(c.dy() < c.dx()), dir(d.get_norm()[dim] > 0.0);
		bool const is_garage(d.type == tquad_with_ix_t::TYPE_GDOOR);
		cube_t step(get_step_for_ext_door(d));
		step.d[dim][!dir] -= (dir ? 1.0 : -1.0)*door_shift_dist; // shift slightly away from the building to prevent Z-fighting with the exterior wall
		assert(step.is_strictly_normalized());
		unsigned const obj_ix(objs.size());
		// must draw the bottom surface in case it's on a hill, unless this is a city building
		unsigned const flags(RO_FLAG_EXTERIOR | (above_ground ? RO_FLAG_ADJ_BOT : 0) | ((is_in_city && !above_ground) ? RO_FLAG_ADJ_LO : 0));
		room_obj_shape const shape(is_garage ? SHAPE_ANGLED : SHAPE_CUBE); // garage door has a sloped ramp
		objs.emplace_back(step, TYPE_EXT_STEP, 0, dim, !dir, flags, 1.0, shape, step_color);
		if (above_ground && !is_garage) {to_add_stairs.push_back(obj_ix);} // add steps up to this door
	} // for d
	if (to_add_stairs.empty()) return; // done
	cube_t const &part(parts[0]); // assumes door is on parts[0] (single part)
	bool const add_step_gaps(objs.size() & 1); // something random-ish per building
	float const base_step_height(floor_spacing/NUM_STAIRS_PER_FLOOR), head_clearance(0.8*get_floor_ceil_gap()), railing_thickness(0.5*get_wall_thickness());
	cube_t stairs_bcube(part);
	stairs_bcube.expand_by_xy(0.25*floor_spacing); // allow stairs to go slightly outside the building bcube, but not enough to need a railing on the other side
	vect_cube_t cand_steps;
	vector<room_object_t> railings;

	// add stairs going to upper level doors
	for (unsigned ix : to_add_stairs) {
		room_object_t &s(objs[ix]);
		float const delta_z(s.z1() - ground_floor_z1);
		if (delta_z <= base_step_height) continue; // no stairs are needed; should always be false based on the above_ground check above
		bool const dim(s.dim), dir(s.dir);
		unsigned const flags(s.flags | RO_FLAG_HANGING); // draw the side facing the building because it may be visible through a window
		unsigned const num_steps(round_fp(delta_z/base_step_height));
		unsigned const num_floors(round_fp((s.z1() - ground_floor_z1)/floor_spacing));
		float const step_height(delta_z/num_steps), max_step_len(2.25*step_height), init_step_len(s.get_sz_dim(!dim)), step_overlap(1.0*step_height), dir_sign(dir ? 1.0 : -1.0);
		bool step_dir(s.get_center_dim(!dim) < part.get_center_dim(!dim)); // preferred steps go down toward longer wall segment
		s.z1() = s.z2() - step_height; // set correct step height for the first step
		cube_t const door_step(s);
		colorRGBA const railing_color(BLACK);
		float const railing_inside_edge(door_step.d[dim][!dir] + dir_sign*railing_thickness);
		cube_t top_railing(door_step);
		top_railing.z2() += floor_spacing;
		top_railing.z1() += step_height;
		bool success(0);
		// Note: s reference is invalidated beyond this point

		for (unsigned d = 0; d < 2; ++d) { // try both dirs
			float const sdir_sign(step_dir ? 1.0 : -1.0), init_translate(sdir_sign*(init_step_len - step_overlap));
			float step_len(init_step_len), max_step_len_dir(max_step_len);
			cube_t step(door_step);
			step.d[dim][dir] -= dir_sign*door_shift_dist; // move slightly away from the building to prevent Z-fighting with interior wall
			cand_steps.clear();
			// constrain steps to fit inside the building stairs_bcube by making them steeper if needed
			min_eq(max_step_len_dir, (step_overlap + fabs(door_step.d[!dim][step_dir] - stairs_bcube.d[!dim][step_dir])/num_steps));

			if (step_len > max_step_len_dir) { // shorten steps if they're too long
				step.d[!dim][step_dir] -= sdir_sign*(step_len - max_step_len_dir);
				step_len = max_step_len_dir;
			}
			vector3d translate(0.0, 0.0, -step_height);
			translate[!dim] = sdir_sign*(step_len - step_overlap); // overlap by step_height
			step.translate_dim(!dim, (init_translate - translate[!dim])); // first translate
			success = 1;

			for (unsigned n = 0; n <= num_steps; ++n) { // one extra iteration to check for collisions at the bottom
				step += translate;
				cube_t check_cube(step);
				check_cube.z2() += head_clearance;

				// check for collisions with previous steps, balconies, and fire escapes
				for (auto i = objs.begin()+ext_objs_start; i != objs.end(); ++i) {
					if ((i - objs.begin()) == ix) continue; // skip our starting step
					cube_t no_block(*i);
					no_block.z2() += head_clearance; // the other object can be walked on as well
					if (check_cube.intersects(no_block)) {success = 0; break;}
				}
				if (!success) break;

				for (roof_obj_t const &ro : details) {
					if (ro.type == ROOF_OBJ_AC && check_cube.intersects(ro)) {success = 0; break;} // check AC unit (may be on the ground rather than rooftop)
				}
				if (!success) break;

				if (has_chimney == 2) { // exterior chimney
					if (check_cube.intersects(get_chimney  ())) {success = 0; break;} // chimney
					if (check_cube.intersects(get_fireplace())) {success = 0; break;} // fireplace
				}
				// what about residential city objects such as fences and trashcans? currently, we don't have multi-family houses in cities, so maybe it's okay
				if (n < num_steps) {cand_steps.push_back(step);} // don't add the last step
			} // for n
			if (!success) {step_dir ^= 1; continue;} // try other dir
			assert(!cand_steps.empty());

			for (cube_t const &cand : cand_steps) { // cube, dim, step_dir, wall_dir, at_door
				cube_t step(cand);
				if (add_step_gaps) {step.z1() += 0.5*cand.dz();} // move the bottom halfway up
				objs.emplace_back(step, TYPE_EXT_STEP, 0, dim, dir, flags, 1.0, SHAPE_CUBE, step_color);
				ext_steps.emplace_back(cand, !dim, step_dir, dir); // Note: here dim is wall dim, but we store stairs dim
				
				if (step.z1() > ground_floor_z1+0.5*floor_spacing && step.z1() < ground_floor_z1+0.75*floor_spacing) {
					cube_t collider(step);
					set_cube_zvals(collider, ground_floor_z1, (ground_floor_z1 + 0.5*(step.z1() - ground_floor_z1)));
					details.emplace_back(collider, DETAIL_OBJ_COLLIDER);
				}
				details.emplace_back(step, DETAIL_OBJ_SHAD_ONLY);
			} // for cand
			ext_steps.back().is_base = ext_steps.back().at_ground = 1; // last step is at the bottom
			// add side railing, with balusters
			cube_t railing(cand_steps[min(size_t(1), cand_steps.size()-1U)]); // second from the top
			railing.d[!dim][!step_dir] -= sdir_sign*railing_thickness; // lengthen slightly to meet the top railing
			railing.union_with_cube(cand_steps.back());
			railing.d[dim][dir] = railing_inside_edge;
			railing.z1() += step_height;
			railing.z2()  = railing.z1() + num_floors*floor_spacing;
			railings.emplace_back(railing, TYPE_RAILING, 0, !dim, !step_dir, (RO_FLAG_OPEN | RO_FLAG_EXTERIOR), 1.0, SHAPE_CUBE, railing_color);
			railings.back().item_flags = max(num_floors, 1U) - 1; // store the number of floors-1 in item_flags
			// add end railing
			railing = top_railing;
			railing.d[!dim][step_dir]  = railing.d[!dim][!step_dir] + sdir_sign*railing_thickness;
			railing.d[ dim][     dir] -= dir_sign*railing_thickness; // move away from the building
			railings.emplace_back(railing, TYPE_RAILING, 0, dim, dir, (RO_FLAG_TOS | RO_FLAG_ADJ_BOT | RO_FLAG_EXTERIOR), 1.0, SHAPE_CUBE, railing_color);
			// add top railing
			railing.d[ dim][     dir] = railing_inside_edge;
			railing.d[!dim][step_dir] = door_step.d[!dim][step_dir];
			railings.emplace_back(railing, TYPE_RAILING, 0, !dim, !step_dir, (RO_FLAG_TOS | RO_FLAG_EXTERIOR), 1.0, SHAPE_CUBE, railing_color);
			break; // done
		} // for d
		if (!success) { // failed to connect stairs in either dim; add railings and turn into a tiny balcony instead
			for (unsigned d = 0; d < 2; ++d) { // add end railings
				cube_t railing(top_railing);
				railing.d[!dim][d  ]  = railing.d[!dim][!d] + (d ? 1.0 : -1.0)*railing_thickness;
				railing.d[ dim][dir] -= dir_sign*railing_thickness; // move away from the building
				railings.emplace_back(railing, TYPE_RAILING, 0, dim, dir, (RO_FLAG_TOS | RO_FLAG_ADJ_BOT | RO_FLAG_EXTERIOR), 1.0, SHAPE_CUBE, railing_color);
			}
			// add top railing
			cube_t railing(top_railing);
			railing.d[dim][dir] = railing_inside_edge;
			railings.emplace_back(railing, TYPE_RAILING, 0, !dim, !step_dir, (RO_FLAG_TOS | RO_FLAG_EXTERIOR), 1.0, SHAPE_CUBE, railing_color);
		}
		ext_steps.emplace_back(door_step, !dim, step_dir, dir, 1); // add the door step; here dim is wall dim, but we store stairs dim
		details.emplace_back(door_step, DETAIL_OBJ_SHAD_ONLY);
	} // for ix
	vector_add_to(railings, objs); // add railings at the end
	for (ext_step_t const &step : ext_steps) {union_with_coll_bcube(step);}
	if (!ext_steps.empty()) {invalidate_tile_smap_in_region(bcube + get_camera_coord_space_xlate());}
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
	case RTYPE_BED : case RTYPE_MASTER_BED: add_window_blinds(window, dim, dir, room_id, floor_ix); break; // bedroom
	case RTYPE_BATH: case RTYPE_MENS: case RTYPE_WOMENS: add_bathroom_window(window, dim, dir, room_id, floor_ix); break; // bathroom
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
	c.translate_dim(dim, (dir ? 1.0 : -1.0)*0.5*get_trim_thickness()); // half the previous translate to prevent Z-fighting in mirror reflections
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

unsigned get_L_stairs_first_flight_count(stairs_landing_base_t const &s, float landing_width) {
	float const length(s.get_length() - landing_width), width(s.get_width() - landing_width), tot_len(length + width);
	assert(length > 0.0 && width > 0.0);
	unsigned const num_stairs_add(s.get_num_stairs() - 1); // Note: the -1 accounts for the landing, which servers as a stair
	assert(num_stairs_add > 1); // must be at least one in each dim
	unsigned const len_ratio(round_fp(num_stairs_add*length/tot_len));
	return max(1U, min(num_stairs_add-1, len_ratio));
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), half_thick(0.5*floor_thickness);
	float const wall_thickness(get_wall_thickness()), elevator_car_z1_add(0.05*floor_thickness), fc_thick_scale(get_elevator_fc_thick_scale());
	float const stairs_sign_width(1.0*wall_thickness);
	vect_room_object_t &objs(interior->room_geom->objs);
	ostringstream oss; // reused across elevators/floors

	// add floor signs for U-shaped stairs
	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator || i->for_ramp || !i->is_u_shape() || i->not_an_exit) continue; // not U-shaped stairs, or no exit
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
				sign.translate_dim(!i->dim, 0.5*(i->get_width() - stairs_sign_width)); // translate to the positive side
				flags |= RO_FLAG_ADJ_TOP; // dra the top surface
			}
			objs.emplace_back(sign, TYPE_SIGN, 0, i->dim, !i->dir, flags, 1.0, SHAPE_CUBE, DK_BLUE); // no room_id
			set_floor_text_for_sign(objs.back(), (real_floor + 1), floor_offset, has_parking_garage, oss);
		}
	} // for i

	// add elevator lights, and signs on each floor; must be done before setting buttons_start
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		bool const dim(i->dim), dir(i->dir);
		// add light
		i->light_obj_id = objs.size();
		float const light_zval(max(ground_floor_z1, i->z1()) + elevator_car_z1_add + (1.0 - fc_thick_scale)*window_vspacing); // starts on the ground floor
		cube_t light(point(i->xc(), i->yc(), light_zval));
		light.z1() -= 0.02*window_vspacing;
		light.expand_by_xy(0.06*window_vspacing);
		objs.emplace_back(light, TYPE_LIGHT, i->room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_IN_ELEV | RO_FLAG_LIT), 0.0, SHAPE_CYLIN, WHITE);
		objs.back().obj_id = uint16_t(i - interior->elevators.begin()); // encode elevator index as obj_id
		// add floor signs
		unsigned const num_floors(calc_num_floors(*i, window_vspacing, floor_thickness));
		unsigned const floor_offset(calc_floor_offset(i->z1()));
		float const ewidth(i->get_width()), dsign(dir ? 1.0 : -1.0), front_wall(i->d[dim][dir]);
		cube_t sign;
		sign.d[dim][0] = sign.d[dim][1] = front_wall;
		sign.d[dim][dir] += dsign*0.1*wall_thickness; // front of sign
		set_wall_width(sign, (i->d[!dim][1] - 0.1*ewidth), 0.04*ewidth, !dim); // to the high side, opposite the call button

		for (unsigned f = 0; f < num_floors; ++f) { // Note: floor number starts at 1 even if the elevator doesn't extend to the ground floor
			if (i->skip_floor_ix(f)) continue;
			sign.z1() = i->z1()   + (f + 0.5)*window_vspacing;
			sign.z2() = sign.z1() + 0.1*ewidth;
			objs.emplace_back(sign, TYPE_SIGN, i->room_id, dim, dir, (RO_FLAG_NOCOLL /*| RO_FLAG_HAS_EXTRA*/), 1.0, SHAPE_CUBE, DK_BLUE); // no frame?
			set_floor_text_for_sign(objs.back(), f+1, floor_offset, has_parking_garage, oss);
		}
		// add concrete flooring at the base of the elevator, over the carpet
		cube_t efloor(*i);
		efloor.expand_by_xy(-0.25*half_thick); // less than the width of the elevator walls, to prevent z-fighting
		efloor.z1() += half_thick; // on top of the carpet
		efloor.z2()  = efloor.z1() + 0.01*half_thick; // set thickness (very thin)
		objs.emplace_back(efloor, TYPE_FLOORING, i->room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_IN_ELEV), 1.0, SHAPE_CUBE, WHITE, FLOORING_CONCRETE);

		if (i->skip_floors_mask > 0) { // block off any unreachable floors
			for (unsigned f = 0; f < 64; ++f) {
				if (!(i->skip_floors_mask & (1ULL << f))) continue;
				cube_t wall(*i);
				wall.d[dim][0] = wall.d[dim][1] = front_wall;
				wall.d[dim][dir] += 0.5*dsign*wall_thickness;
				wall.z1() = i->z1() + f*window_vspacing - half_thick;
				wall.z2() = wall.z1() + window_vspacing;
				objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, !dir, RO_FLAG_HANGING); // hanging so that the bottom surface is drawn
			} // for f
		}
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
			if (i->skip_floor_ix(f)) continue;
			point pos;
			pos[ i->dim] = i->d[i->dim][i->dir]; // front of the elevator
			pos[!i->dim] = i->d[!i->dim][0] + 0.1*ewidth; // to the low side

			for (unsigned d = 0; d < 2; ++d) { // {down, up} call buttons
				if ((d == 0 && f == 0) || (d == 1 && f == num_floors-1)) continue; // no floor above/below
				pos.z = i->z1() + (f + 0.05*d + 0.45)*window_vspacing;
				add_elevator_button(pos, button_radius, i->dim, i->dir, elevator_id, f, 0, d, objs); // inside=0, is_up=d
			}
		} // for f
		// call buttons for each floor inside the elevator car; first find the panel location for the starting elevator car position;
		// floor numbers are added in building_room_geom_t::add_elevator();
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
		float cur_z(panel.z1() + button_spacing);
		
		for (unsigned f = 0; f < num_floors; ++f) {
			if (i->skip_floor_ix(f)) continue; // also skips cur_z update to avoid a gap in the buttons, but there's still a gap in the floor numbers
			pos.z  = cur_z;
			cur_z += button_spacing;
			add_elevator_button(pos, inner_button_radius, i->dim, !i->dir, elevator_id, f, 1, 0, objs); // inside=1, is_up=0, pointing in opposite dir
		}
		i->button_id_end = objs.size();
	} // for e
	interior->room_geom->stairs_start  = objs.size();
	interior->room_geom->has_elevators = (!interior->elevators.empty());
	colorRGBA const railing_colors[3]  = {GOLD, DK_BROWN, BLACK};
	colorRGBA const railing_color(railing_colors[rgen.rand()%3]); // set per-building

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator || i->for_ramp) continue; // for elevator or ramp, not stairs
		unsigned const num_stairs(i->get_num_stairs());
		float const stair_dz(window_vspacing/(num_stairs+1)), stair_height(stair_dz + floor_thickness), stair_z1h(0.4f*stair_height);
		bool const dim(i->dim), dir(i->dir), has_side_walls(i->has_walled_sides() || i->is_u_shape());
		bool const has_wall_both_sides(i->against_wall[0] && i->against_wall[1]); // ext basement stairs
		bool const side(dir); // for U-shaped stairs; for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
		// Note: stairs always start at floor_thickness above the landing z1, ignoring landing z2/height
		float const floor_z(i->z1() + floor_thickness - window_vspacing), step_len_pos(i->get_step_length());
		float const wall_hw(min(STAIRS_WALL_WIDTH_MULT*max(step_len_pos, stair_dz), 0.25f*stair_dz));
		float const stairs_zmin(i->in_ext_basement ? interior->basement_ext_bcube.z1() : bcube.z1());
		float step_len((dir ? 1.0 : -1.0)*step_len_pos), z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i), landing; // Note: landing is for L-shaped stairs
		unsigned num_stairs1(0), num_stairs2(0); // used for L-shaped stairs

		if (i->is_straight()) { // straight stairs
			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				// tag basement bottom stair so that no width extension is added, since this may clip through the basement door
				unsigned const flags((n == 0 && !i->in_ext_basement && i->z1() < ground_floor_z1) ? RO_FLAG_RSTAIRS : 0);
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				set_cube_zvals(stair, max(stairs_zmin, z+stair_z1h), z+stair_height); // don't go below the floor (Note: z1 was (z + 0.5f*half_thick))
				assert(stair.z1() < stair.z2());
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir, flags); // Note: room_id=0, not tracked, unused
			} // for n
		}
		else if (i->is_u_shape()) { // U-shaped stairs
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
		else if (i->is_l_shape()) {
			bool const dir2(i->bend_dir);
			// determine the number of stairs in each dim by looking at stairs spans
			float const landing_width(get_landing_width());
			num_stairs1 = get_L_stairs_first_flight_count(*i, landing_width);
			num_stairs2 = num_stairs-1 - num_stairs1;
			// set width of first flight equal to the width of the landing
			float const inner_edge(stair.d[!dim][!dir2] + (dir2 ? 1.0 : -1.0)*landing_width);
			stair.d[!dim][dir2] = inner_edge;
			// add stairs in dim
			step_len = (dir ? 1.0 : -1.0)*(i->get_length() - landing_width)/num_stairs1;

			for (unsigned n = 0; n < num_stairs1; ++n, z += stair_dz, pos += step_len) {
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				set_cube_zvals(stair, max(stairs_zmin, z+stair_z1h), z+stair_height); // don't go below the floor
				assert(stair.z1() < stair.z2());
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir, 0);
			}
			// add landing
			stair.d[dim][!dir] = pos; stair.d[dim][dir] = i->d[dim][dir]; // extend to the end
			set_cube_zvals(stair, z+stair_z1h, z+stair_height);
			objs.emplace_back(stair, TYPE_STAIR, 0, !dim, !dir, 0); // place in !dim but use first dir, so that player coll is enabled since we don't have a wall here
			objs.back().shape = SHAPE_STAIRS_L; // only the landing needs a tag because it can be entered from right angle sides
			landing = stair;
			z += stair_dz;
			// add wall supporting landing
			cube_t wall(stair);
			set_cube_zvals(wall, floor_z, stair.z2()); // top of floor below to top of landing
			wall.d[dim][!dir] = i->d[dim][dir]; // against landing
			wall.d[dim][ dir] = i->d[dim][dir] + (dir ? 1.0 : -1.0)*wall_thickness; // extend out by thickness
			wall.expand_in_dim(!dim, 0.05*wall_thickness); // expand slightly to prevent Z-fighting with side of stair and end of railing
			objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir, 0);
			// add stairs in !dim
			pos = inner_edge;
			step_len = (dir2 ? 1.0 : -1.0)*(i->get_width() - landing_width)/num_stairs2;

			for (unsigned n = 0; n < num_stairs2; ++n, z += stair_dz, pos += step_len) {
				stair.d[!dim][!dir2] = pos; stair.d[!dim][dir2] = pos + step_len;
				set_cube_zvals(stair, z+stair_z1h, z+stair_height);
				objs.emplace_back(stair, TYPE_STAIR, 0, !dim, dir2, 0);
			}
		}
		else {assert(0);}
		// add walls and railings
		bool const extend_walls_up(i->is_at_top && !i->roof_access); // space above is open, add a wall so that people can't fall down the stairs
		float const railing_z2(i->z2() + (i->roof_access ? 0.025*i->dz() : 0.0)); // capture z2 before we change it; move roof access railing up a bit to offset the shrink resize
		float const wall_bottom(floor_z - half_thick), railing_side_dz(0.5*stair_dz); // for U-shaped stairs
		cube_t wall(*i);
		if (extend_walls_up) {wall.z2() += window_vspacing - floor_thickness;}
		else {wall.z2() -= 0.5*floor_thickness;} // prevent z-fighting on top floor
		wall.z1() = max((stairs_zmin + half_thick), wall_bottom); // full height
		set_wall_width(wall, i->d[dim][dir], wall_hw, dim);
		float walls_extend_to(0.0);

		if ((i->shape == SHAPE_WALLED && !(i->against_wall[0] || i->against_wall[1]) && (!i->stack_conn || !i->is_at_top)) || i->is_u_shape()) {
			objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir); // add wall at back/end of stairs

			if (i->not_an_exit && i->is_u_shape()) { // blocked U-shaped stairs
				cube_t front_wall(*i);
				front_wall.z2() -= floor_thickness;
				front_wall.z1()  = floor_z - 0.5*floor_thickness;
				front_wall.expand_in_dim(!dim, wall_hw); // widen slightly
				// create a box for the landing so that the player and AI can walk there when changing floors
				float const landing_width(i->get_retail_landing_width(window_vspacing)), front_pos(i->d[dim][!dir]);
				set_wall_width(front_wall, (front_pos + (dir ? -1.0 : 1.0)*(wall_hw + landing_width)), wall_hw, dim); // move to the front
				objs.emplace_back(front_wall, TYPE_STAIR_WALL, 0, dim, !dir, 0); // add wall in front of stairs
				walls_extend_to = front_wall.d[dim][dir];
				cube_t sub_floor(*i);
				sub_floor.d[dim][ dir] = front_pos - 0.25*(dir ? -1.0 : 1.0)*wall_hw; // move toward stairs slightly to prevent Z-fighting and fill the gap
				sub_floor.d[dim][!dir] = front_wall.d[dim][!dir]; // furthest extent of wall
				sub_floor.z1() = front_wall.z1() - half_thick;
				sub_floor.z2() = sub_floor.z1() + floor_thickness;
				sub_floor.expand_in_dim(!dim, wall_hw); // cover the bottoms of the side walls
				// add as both a stairs wall (for drawing) and a floor (for player collision)
				objs.emplace_back(sub_floor, TYPE_STAIR_WALL, 0, dim, !dir, RO_FLAG_HANGING); // hanging so that the bottom surface is drawn
				interior->floors.push_back(sub_floor);
			}
		}
		else if (i->has_walled_sides() && extend_walls_up) { // add upper section only
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
			if (i->is_l_shape()) continue; // nothing to add in this loop
			set_wall_width(wall, i->d[!dim][d], wall_hw, !dim);
			wall.expand_in_dim(dim, 0.01*wall_hw); // just enough to avoid z-fighting with stairs
			bool const add_wall(has_side_walls && !i->against_wall[d]); // don't add a wall if the stairs are already against a wall
			
			if (add_wall) { // add walls around stairs for this floor
				// clip basement stairs wall to the basement to avoid drawing artifacts at the bottom of a house exterior wall; or just skip drawing the wall?
				cube_t wall_clipped(wall);
				if (walls_extend_to != 0.0) {wall_clipped.d[dim][!dir] = walls_extend_to;} // extend outward to meet front wall
				if (i->z1() < ground_floor_z1 && !i->in_ext_basement) {wall_clipped.intersect_with_cube_xy(get_basement());}
				assert(wall_clipped.is_strictly_normalized());
				objs.emplace_back(wall_clipped, TYPE_STAIR_WALL, 0, dim, dir);
			}
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
				if (i->is_u_shape()) { // adjust railing height/angle to match stairs
					flags |= RO_FLAG_HAS_EXTRA; // make it taller
					float const z_split(railing.zc());
					if (bool(d) == side) {railing.z1() = z_split + railing_side_dz; flags |= RO_FLAG_ADJ_HI; railing_dir ^= 1;}
					else                 {railing.z2() = z_split - railing_side_dz; flags |= RO_FLAG_ADJ_LO;}
				}
				objs.emplace_back(railing, TYPE_RAILING, 0, dim, railing_dir, flags, 1.0, SHAPE_CUBE, railing_color);
			}
		} // for d
		if (i->has_railing && i->is_u_shape()) { // add a railing for the back wall of U-shaped stairs
			float const railing_zc(wall_bottom + 0.819*window_vspacing); // determined experimentally
			cube_t railing(*i);
			set_wall_width(railing, (i->d[dim][dir] + (dir ? -1.0 : 1.0)*2.0*wall_hw), wall_hw, dim);
			set_wall_width(railing, railing_zc, 1.4*railing_side_dz, 2); // set zvals
			unsigned const railing_flags(RO_FLAG_NOCOLL | RO_FLAG_ADJ_HI | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_BOT | RO_FLAG_HAS_EXTRA); // make taller
			objs.emplace_back(railing, TYPE_RAILING, 0, !dim, dir, railing_flags, 1.0, SHAPE_CUBE, railing_color); // no ends
		}
		else if (i->has_railing && i->is_l_shape()) { // add railings to the sides of each segment, along the hole at the top on 3 sides, and around the two landing sides
			bool const dir2(i->bend_dir);
			// add side railings
			cube_t segs[2] = {*i, *i}; // {dim/dir, !dim/dir2}
			set_cube_zvals(segs[0], floor_z, landing.z1()+stair_dz);
			set_cube_zvals(segs[1], landing.z2(), railing_z2);
			segs[0].d[!dim][ dir2] = segs[1].d[!dim][!dir2] = landing.d[!dim][ dir2]; // inside of lower/first  flight
			segs[1].d[ dim][!dir ] = segs[0].d[ dim][ dir ] = landing.d[ dim][!dir ]; // inside of upper/second flight
			float const railing_hw(0.75*wall_hw);
			bool const dirs[2] = {dir, dir2};

			for (unsigned d = 0; d < 2; ++d) { // sides of stairs
				for (unsigned s = 0; s < 2; ++s) { // segs
					bool const rdim(dim ^ bool(s)), rdir(dirs[s]);
					cube_t railing(segs[s]);
					set_wall_width(railing, segs[s].d[!rdim][d], railing_hw, !rdim);
					objs.emplace_back(railing, TYPE_RAILING, 0, rdim, rdir, RO_FLAG_OPEN, 1.0, SHAPE_CUBE, railing_color); // with balusters
					objs.back().state_flags = (s ? num_stairs2+1 : num_stairs1); // encode num_stairs in state_flags so that railing height can be drawn correctly
				}
			}
			// add landing railings
			unsigned const railing_flags(RO_FLAG_OPEN | RO_FLAG_HAS_EXTRA | RO_FLAG_TOS); // make taller, with balusters

			for (unsigned d = 0; d < 2; ++d) {
				bool const rdim(dim ^ bool(d)), rdir(dirs[d] ^ bool(d));
				cube_t railing(landing);
				set_cube_zvals(railing, landing.z2(), railing_z2);
				set_wall_width(railing, landing.d[rdim][rdir], railing_hw, rdim);
				objs.emplace_back(railing, TYPE_RAILING, 0, !rdim, rdir, railing_flags, 1.0, SHAPE_CUBE, railing_color);
			}
			// add top railings
			cube_t railing(*i);
			set_cube_zvals(railing, i->z2(), (i->z2() + window_vspacing - floor_thickness)); // starts at the floor
			if (!i->is_at_top) {railing.d[!dim][!dir2] = landing.d[!dim][dir2];} // inside edge of landing to make room for the flight above, if there is one
			
			// upper end next to landing and lower end; vertical poles are only needed for one end, except for the top floor
			for (unsigned d = 0; d < 2; ++d) {
				set_wall_width(railing, i->d[dim][d], railing_hw, dim);
				unsigned const flags(RO_FLAG_TOS | ((bool(d) ^ dir ^ 1) ? 0 : RO_FLAG_OPEN)); // balusters on one side
				objs.emplace_back(railing, TYPE_RAILING, 0, !dim, d, flags, 1.0, SHAPE_CUBE, railing_color);
			}
			// long edge
			railing.d[dim][ dir] = landing.d[dim][!dir]; // ends at the landing/upper stairs
			railing.d[dim][!dir] = i->d[dim][!dir]; // ends at stairs cutout to meet the other railing
			railing.expand_in_dim(dim, railing_hw); // extend to cover the vertical poles
			set_wall_width(railing, i->d[!dim][dir2], railing_hw, !dim);
			objs.emplace_back(railing, TYPE_RAILING, 0, dim, dir2, (RO_FLAG_TOS | RO_FLAG_OPEN | RO_FLAG_ADJ_TOP), 1.0, SHAPE_CUBE, railing_color); // no vertical pole
		}
		else if (i->has_railing && !has_wall_both_sides && (i->stack_conn || (extend_walls_up && i->shape == SHAPE_STRAIGHT))) {
			// add railings around the top if: straight + top floor with no roof access, connector stairs, or basement stairs
			room_object_t railing(*i, TYPE_RAILING, 0, !dim, dir, (RO_FLAG_TOS | RO_FLAG_ADJ_BOT), 1.0, SHAPE_CUBE, railing_color); // flag to skip drawing ends
			set_cube_zvals(railing, i->z2(), (i->z2() + window_vspacing - floor_thickness)); // starts at the floor
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
	if (has_pool()) { // add pool stairs
		interior->room_geom->pool_stairs_start_ix = objs.size();
		indoor_pool_t const &pool(interior->pool);
		cube_t pool_shallow(pool);
		pool_shallow.z1() = pool.shallow_zval;
		bool const dim(pool.dim), dir(pool.dir);
		float const stairs_height(window_vspacing/(NUM_STAIRS_PER_FLOOR+1)), pool_depth(pool_shallow.dz());
		unsigned const num_stairs(round_fp(pool_depth/stairs_height)); // same spacing is regular stairs
		float const step_height(pool_depth/(num_stairs+1)), step_stride((dir ? -1.0 : 1.0)*1.2*step_height); // last step up to the edge counts
		cube_t step(pool_shallow); // copy the correct width (spans to entire pool width)
		step.d[dim][!dir] = pool.d[dim][dir] + step_stride; // extend into pool

		for (unsigned n = 0; n < num_stairs; ++n) {
			step.z2() -= step_height; // shift down first, since the first step is below the pool edge
			objs.emplace_back(step, TYPE_STAIR, pool.room_ix, dim, !dir, RO_FLAG_IN_POOL);
			if (n+1 < num_stairs) {step.translate_dim(dim, step_stride);} // don't need to translate the last step
		}
		// add stairs railings
		cube_t railing(pool_shallow);
		railing.z2() += 0.5*step_height + get_trim_thickness(); // starts on pool deck
		railing.d[ dim][ dir] -= 0.5*step_stride; // on the pool deck
		railing.d[ dim][!dir]  = step.d[dim][!dir]; // to the end of the last step
		railing.expand_in_dim(!dim, -0.5*wall_thickness); // shrink slightly
		float const positions[3] = {railing.d[!dim][0], railing.d[!dim][1], railing.get_center_dim(!dim)}; // lo, hi, center
		unsigned const num_railings(2 + rgen.rand_bool()); // 2-3

		for (unsigned d = 0; d < num_railings; ++d) { // each side of the pool, and maybe the center
			set_wall_width(railing, positions[d], 0.375*wall_thickness, !dim);
			objs.emplace_back(railing, TYPE_RAILING, pool.room_ix, dim, dir, RO_FLAG_IN_POOL, 1.0, SHAPE_CUBE, GOLD); // no balusters
		}
	}
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

void building_t::add_doorbell_lamp_and_porch_items(tquad_with_ix_t const &door, rand_gen_t &rgen) { // and porch packages
	cube_t const door_bcube(door.get_bcube());
	bool const dim(door_bcube.dy() < door_bcube.dx());
	int const dir_ret(get_ext_door_dir(door_bcube, dim));
	if (dir_ret > 1) return; // not found, skip doorbell, lamp, and chair
	// add doorbell
	bool dir(dir_ret != 0);
	bool const side(dir ^ dim); // currently always to the right, which matches the door handle side
	float const door_width(door_bcube.get_sz_dim(!dim)), half_width(0.016*door_width), half_height(1.8*half_width), button_thickness(0.1*half_width);
	float const db_zval(door_bcube.z1() + 0.55*door_bcube.dz()), floor_spacing(get_window_vspace());
	float const pos(door_bcube.d[!dim][side] + (side ? 1.0 : -1.0)*5.0*half_width), dsign(dir ? 1.0 : -1.0);
	cube_t c;
	c.d[dim][0  ]  = c.d[dim][1] = door_bcube.d[dim][dir] - 0.02*dsign*floor_spacing; // slightly in front of exterior wall
	c.d[dim][dir] += dsign*button_thickness;
	set_cube_zvals(c, (db_zval - half_height), (db_zval + half_height));
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
		lamp.d[dim][dir] = c.d[dim][!dir] + dsign*depth;
		set_cube_zvals(lamp, z1, (z1 + height));
		set_wall_width(lamp, lamp_pos, 0.5*width, !dim);
		unsigned flags(base_flags);
		if (rgen.rand_bool()) {flags |= RO_FLAG_LIT;} // light is on 50% of the time
		objs.emplace_back(lamp, TYPE_WALL_LAMP, room_id, dim, dir, flags, tot_light_amt);
		if (objs.back().is_lit()) {ext_lights.emplace_back(lamp.get_cube_center(), 20.0*width, WALL_LAMP_COLOR);}
	}
	if (has_porch()) { // add porch items
		// find the front door
		float const wall_thickness(get_wall_thickness());
		cube_t front_door;

		for (auto const &d : doors) {
			cube_t dbc(d.get_bcube());
			dbc.expand_by(wall_thickness);
			if (porch.intersects(dbc)) {front_door = dbc;}
		}
		if (front_door.is_all_zeros()) return; // should never fail?
		float const rand_val(rgen.rand_float());

		if (rand_val > 0.6) { // maybe add a package on the porch
			vector3d sz; // half size relative to window_vspacing
			gen_crate_sz(sz, rgen, floor_spacing);
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
		else if (rand_val < 0.4 && building_obj_model_loader.is_model_valid(OBJ_MODEL_RCHAIR)) { // add a chair on the porch
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RCHAIR)); // D, W, H
			float const height(0.5*floor_spacing), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z), zbot(porch.z2());
			float const back_pos(porch.d[dim][!dir] + dsign*0.1*depth);
			bool const chair_side(!side); // opposite the doorbell
			float const center(door_bcube.d[!dim][chair_side] + (chair_side ? 1.0 : -1.0)*2.25*hwidth);
			vect_room_object_t &objs(interior->room_geom->objs);
			cube_t chair;
			set_cube_zvals(chair, zbot, zbot+height);
			chair.d[dim][!dir] = back_pos;
			chair.d[dim][ dir] = back_pos + dsign*depth;
			set_wall_width(chair, center, hwidth, !dim);
			
			if (porch.contains_cube_xy(chair)) {
				objs.emplace_back(chair, TYPE_RCHAIR, room_id, dim, dir, RO_FLAG_EXTERIOR, tot_light_amt);
				details.emplace_back(chair, DETAIL_OBJ_COLLIDER);
			}
		}
	}
}

room_t::room_t(cube_t const &c, unsigned p, unsigned nl, bool is_hallway_, bool is_office_, bool is_sec_bldg_) :
	cube_t(c), no_geom(is_hallway_), // no geom in hallways
	is_hallway(is_hallway_), is_office(is_office_), is_sec_bldg(is_sec_bldg_), // secondary buildings are always a single floor
	is_single_floor(is_sec_bldg), part_id(p), num_lights(nl)
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
	// assign unless already set to a bathroom, unless we're refining the bathroom type to men's or women's
	if (is_bathroom(rtype[floor]) && !is_bathroom(rt)) return;
	rtype[floor] = rt;
	if (locked) {rtype_locked |= (1 << floor);} // lock this floor
}
bool room_t::has_room_of_type(room_type type) const {
	for (unsigned n = 0; n < NUM_RTYPE_SLOTS; ++n) {
		if (rtype[n] == type) return 1;
	}
	return 0;
}

