// 3D World - Building Basement and Parking Garage Logic
// by Frank Gennari 03/11/2022

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for car_t


extern city_params_t city_params; // for num_cars

car_t car_from_parking_space(room_object_t const &o);
void subtract_cube_from_floor_ceil(cube_t const &c, vect_cube_t &fs);
colorRGBA get_light_color_temp_range(float tmin, float tmax, rand_gen_t &rgen);
void set_light_xy(cube_t &light, point const &center, float light_size, bool light_dim, room_obj_shape light_shape);


bool enable_parked_cars() {return (city_params.num_cars > 0 && !city_params.car_model_files.empty());}

void subtract_cubes_from_cube_split_in_dim(cube_t const &c, vect_cube_t const &sub, vect_cube_t &out, vect_cube_t &out2, unsigned dim) {
	out.clear();
	out.push_back(c);

	for (auto s = sub.begin(); s != sub.end(); ++s) {
		if (!c.intersects(c)) continue; // no overlap with orig cube (optimization)
		out2.clear();

		// clip all of out against *s, write results to out2, then swap with out
		for (auto i = out.begin(); i != out.end(); ++i) {
			if (!i->intersects(*s)) {out2.push_back(*i); continue;} // no overlap, keep entire cube

			if (i->d[dim][0] < s->d[dim][0]) { // lo side
				out2.push_back(*i);
				out2.back().d[dim][1] = s->d[dim][0];
			}
			if (i->d[dim][1] > s->d[dim][1]) { // hi side
				out2.push_back(*i);
				out2.back().d[dim][0] = s->d[dim][1];
			}
		} // for i
		out.swap(out2);
	} // for s
}

// *** Utilities ***

unsigned building_t::add_water_heaters(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool single_only) {
	float const floor_spacing(get_window_vspace()), height(get_floor_ceil_gap());
	float const radius((is_house ? 0.18 : 0.20)*height); // larger radius for office buildings
	cube_t const room_bounds(get_walkable_room_bounds(room));
	cube_t place_area(room_bounds);
	place_area.expand_by(-(1.05*radius + get_trim_thickness())); // account for the pan
	if (!place_area.is_strictly_normalized()) return 0; // too small to place water heater
	vect_room_object_t &objs(interior->room_geom->objs);
	point center(0.0, 0.0, zval);
	bool first_xdir(rgen.rand_bool()), first_ydir(rgen.rand_bool()), first_dim(0);

	if (is_house) { // random corner selection
		first_dim = rgen.rand_bool();
	}
	else { // select corners furthest from the door (back wall)
		first_dim = (room.dy() < room.dx()); // face the shorter dim for office buildings so that we can place longer rows

		for (door_stack_t const &ds : interior->door_stacks) {
			if (!ds.is_connected_to_room(room_id)) continue;
			bool const dir(room.get_center_dim(ds.dim) < ds.get_center_dim(ds.dim));
			(ds.dim ? first_ydir : first_xdir) = !dir; // opposite the door
			//first_dim = ds.dim; // place against back wall; too restrictive?
		}
	}
	// make 5 attempts to place a water heater - one in each corner and 1 along a random wall for variety; 20 attempts for office buildings/apartments/hotels
	for (unsigned n = 0; n < (is_house ? 5U : 20U); ++n) {
		bool dim(0), dir(0), xdir(0), ydir(0);

		if (n < 4) { // corner
			dim  = (first_dim  ^ (n >= 2));
			xdir = (first_xdir ^ bool(n & 1));
			ydir = (first_ydir ^ bool(n & 1));
			dir  = (dim ? ydir : xdir);
			center.x = place_area.d[0][xdir];
			center.y = place_area.d[1][ydir];
		}
		else { // wall
			dim = rgen.rand_bool();
			dir = rgen.rand_bool(); // choose a random wall
			center[ dim] = place_area.d[dim][dir]; // against this wall
			center[!dim] = rgen.rand_uniform(place_area.d[!dim][0], place_area.d[!dim][1]);
		}
		cube_t c(get_cube_height_radius(center, radius, height));
		if (is_obj_placement_blocked(c, room, 1)) continue;
		cube_t c_exp(c);
		c_exp.expand_by_xy(0.2*radius); // small keepout in XY
		c_exp.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.25*radius; // add more keepout in front where the controls are
		c_exp.intersect_with_cube(room); // don't pick up objects on the other side of the wall
		if (overlaps_other_room_obj(c_exp, objs_start)) continue; // check existing objects, in particular storage room boxes that will have already been placed
		if (zval > ground_floor_z1 && check_if_against_window(c, room, dim, dir)) continue; // can fail for apartment and hotel utility rooms
		unsigned const flags((is_house ? RO_FLAG_IS_HOUSE : 0) | RO_FLAG_INTERIOR);
		objs.emplace_back(c, TYPE_WHEATER, room_id, dim, !dir, flags, tot_light_amt, SHAPE_CYLIN);
		unsigned num_added(1);

		if (is_house && has_attic()) { // add rooftop vent above the water heater
			float const attic_floor_zval(get_attic_part().z2()), vent_radius(0.15*radius);
			point const vent_bot_center(center.x, center.y, attic_floor_zval);
			add_attic_roof_vent(vent_bot_center, vent_radius, room_id, 1.0); // light_amt=1.0; room_id is for the basement because there's no attic room
		}
		if (!single_only && !is_house && n < 4) { // office building, placed at corner; try to add additional water heaters along the wall
			vector3d step;
			step[!dim] = 2.2*radius*((dim ? xdir : ydir) ? -1.0 : 1.0); // step in the opposite dim
			// ideally we would iterate over all the plumbing fixtures to determine water usage, but they might not be placed yet, so instead count rooms;
			// but rooms can be of very different sizes, so maybe counting volume is the best?
			float tot_volume(0.0);
			for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {tot_volume += i->get_volume();}
			tot_volume /= floor_spacing*floor_spacing*floor_spacing;
			unsigned const max_wh(max(1, min(5, round_fp(tot_volume/1000.0))));

			for (unsigned m = 1; m < max_wh; ++m) { // one has already been placed
				c.translate(step);
				if (!room_bounds.contains_cube(c)) break; // went outside the room, done
				if (is_obj_placement_blocked(c, room, 1)) continue; // bad placement, skip
				c_exp.translate(step);
				if (overlaps_other_room_obj(c_exp, objs_start)) continue; // check existing objects
				objs.emplace_back(c, TYPE_WHEATER, room_id, dim, !dir, flags, tot_light_amt, SHAPE_CYLIN);
				++num_added;
			} // for m
		}
		return num_added; // done
	} // for n
	return 0;
}

bool gen_furnace_cand(cube_t const &place_area, float floor_spacing, bool near_wall, rand_gen_t &rgen, cube_t &furnace, bool &dim, bool &dir) {
	// my furnace is 45" tall, 17.5" wide, and 21" deep, with a 9" tall duct below (total height = 54"); assume floor spacing is 96"
	float const height(0.563*floor_spacing), hwidth(0.182*height), hdepth(0.219*height), room_min_sz(min(place_area.dx(), place_area.dy()));
	if (hdepth > 5.0*room_min_sz || 2.1*(hdepth + building_t::get_scaled_player_radius()) > room_min_sz) return 0; // place area is too small
	dim = rgen.rand_bool();
	point center;
	float const lo(place_area.d[dim][0] + hdepth), hi(place_area.d[dim][1] - hdepth);

	if (near_wall) {
		dir = rgen.rand_bool();
		center[dim] = (dir ? lo : hi); // dir is which way it's facing, not which wall it's against
	}
	else {
		center[dim] = rgen.rand_uniform(lo, hi);
		dir = (center[dim] < place_area.get_center_dim(dim)); // face the center of the room (attic or above hallway/stairs)
	}
	center[!dim] = rgen.rand_uniform(place_area.d[!dim][0]+hwidth, place_area.d[!dim][1]-hwidth);
	set_wall_width(furnace, center[ dim], hdepth,  dim);
	set_wall_width(furnace, center[!dim], hwidth, !dim);
	set_cube_zvals(furnace, place_area.z1(), place_area.z1()+height);
	return 1;
}

void building_t::add_breaker_panel(rand_gen_t &rgen, cube_t const &c, float ceil_zval, bool dim, bool dir, unsigned room_id, float tot_light_amt) {
	assert(has_room_geom());
	auto &objs(interior->room_geom->objs);
	objs.emplace_back(c, TYPE_BRK_PANEL, room_id, dim, dir, RO_FLAG_INTERIOR, tot_light_amt, SHAPE_CUBE, colorRGBA(0.5, 0.6, 0.7));
	set_obj_id(objs);

	if (c.z1() < ground_floor_z1) { // add conduit if in the basement
		assert(c.z2() < ceil_zval);
		cube_t conduit(cube_top_center(c));
		conduit.z2() = ceil_zval;
		float const bp_hdepth(0.5*c.get_sz_dim(dim)), conduit_radius(rgen.rand_uniform(0.38, 0.46)*bp_hdepth);
		conduit.expand_by_xy(conduit_radius);
		objs.emplace_back(conduit, TYPE_PIPE, room_id, 0, 1, (RO_FLAG_NOCOLL | RO_FLAG_INTERIOR), tot_light_amt, SHAPE_CYLIN, LT_GRAY); // vertical pipe
	}
}
bool connect_furnace_to_breaker_panel(room_object_t const &furnace, cube_t const &bp, bool dim, bool dir, float conn_height, cube_t &conn) {
	assert(furnace.type == TYPE_FURNACE);
	if (furnace.dim != dim || furnace.dir == dir) return 0; // wrong orient (note that dir is backwards)
	if (furnace.d[dim][dir] < bp.d[dim][0] || furnace.d[dim][dir] > bp.d[dim][1]) return 0; // back against the wrong wall (can't check wall_pos here)
	float const radius(0.67*0.2*bp.get_sz_dim(dim)); // 67% of average conduit radius
	set_wall_width(conn, conn_height, radius, 2); // set zvals
	if (conn.z1() < furnace.z1() || conn.z2() > furnace.z2()) return 0;
	bool const pipe_dir(furnace.get_center_dim(!dim) < bp.get_center_dim(!dim));
	set_wall_width(conn, conn_height, radius, 2);
	float const wall_pos(bp.d[dim][dir]);
	conn.d[dim][ dir] = wall_pos;
	conn.d[dim][!dir] = wall_pos + (dir ? -1.0 : 1.0)*2.0*radius; // extend outward from wall
	conn.d[!dim][!pipe_dir] = furnace.d[!dim][ pipe_dir];
	conn.d[!dim][ pipe_dir] = bp     .d[!dim][!pipe_dir];
	return 1;
}

// Note: for houses
bool building_t::add_basement_utility_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!has_room_geom()) return 0;
	// place water heater
	bool was_added(add_water_heaters(rgen, room, zval, room_id, tot_light_amt, objs_start));
	if (interior->furnace_type == FTYPE_BASEMENT) {was_added |= add_furnace_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);} // place furnace in basement
	return was_added;
}

// Note: this function is here rather than in building_rooms.cpp because utility rooms are connected to utilities in the basement, and it's similar to the code above
bool building_t::add_office_utility_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	zval       = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_CONCRETE); // add concreate and move the effective floor up
	objs_start = interior->room_geom->objs.size(); // exclude this from collision checks
	unsigned const num_water_heaters(add_water_heaters(rgen, room, zval, room_id, tot_light_amt, objs_start));
	if (num_water_heaters == 0 && !is_apt_or_hotel()) return 0; // apartments and hotels have utility rooms even if there's no water heater
	// add one furnace per water heater
	auto &objs(interior->room_geom->objs);
	unsigned const furnaces_start(objs.size());
	for (unsigned n = 0; n < num_water_heaters; ++n) {add_furnace_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);}
	unsigned const furnaces_end(objs.size());
	cube_t place_area(get_walkable_room_bounds(room));
	place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6); // place janitorial sink
	// add breaker panel
	float const floor_spacing(get_window_vspace()), floor_height(floor_spacing - 2.0*get_fc_thickness()), ceil_zval(zval + get_floor_ceil_gap());
	float const bp_hwidth(rgen.rand_uniform(0.15, 0.25)*(is_house ? 0.7 : 1.0)*floor_height), bp_hdepth(rgen.rand_uniform(0.05, 0.07)*(is_house ? 0.5 : 1.0)*floor_height);
	
	if (bp_hwidth < 0.25*min(room.dx(), room.dy())) { // if room is large enough
		cube_t c;
		set_cube_zvals(c, (ceil_zval - 0.75*floor_height), (ceil_zval - rgen.rand_uniform(0.2, 0.35)*floor_height));

		for (unsigned n = 0; n < 20; ++n) { // 20 tries
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
			float const dir_sign(dir ? -1.0 : 1.0), wall_pos(place_area.d[dim][dir]);
			float const bp_center(rgen.rand_uniform(place_area.d[!dim][0]+bp_hwidth, place_area.d[!dim][1]-bp_hwidth));
			set_wall_width(c, bp_center, bp_hwidth, !dim);
			c.d[dim][ dir] = wall_pos;
			c.d[dim][!dir] = wall_pos + dir_sign*2.0*bp_hdepth; // extend outward from wall
			assert(c.is_strictly_normalized());
			cube_t test_cube(c);
			test_cube.d[dim][!dir] += dir_sign*2.0*bp_hwidth; // add a width worth of clearance in the front so that the door can be opened
			test_cube.z2() = ceil_zval; // extend up to ceiling to ensure space for the conduit
			if (is_obj_placement_blocked(test_cube, room, 1) || overlaps_other_room_obj(test_cube, objs_start)) continue;
			if (zval > ground_floor_z1 && check_if_against_window(c, room, dim, dir)) continue; // can fail for apartment and hotel utility rooms
			add_breaker_panel(rgen, c, ceil_zval, dim, dir, room_id, tot_light_amt);
			// connect furnaces on the same wall to the breaker box
			float const conn_height(c.z1() + rgen.rand_uniform(0.25, 0.75)*c.dz());

			for (unsigned f = furnaces_start; f < furnaces_end; ++f) {
				if (objs[f].type == TYPE_BLOCKER) continue; // skip blockers
				cube_t conn;
				if (!connect_furnace_to_breaker_panel(objs[f], c, dim, dir, conn_height, conn)) continue;
				if (is_obj_placement_blocked(conn, room, 1) || overlaps_other_room_obj(conn, objs_start, 0, &furnaces_start)) continue; // check water heaters only
				objs.emplace_back(conn, TYPE_PIPE, room_id, !dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, LT_GRAY); // horizontal
			}
			break; // done
		} // for n
	}
	// don't add signs for interior utility rooms in apartments and hotels
	if (!room.is_apt_or_hotel_room() || room.get_is_entryway()) {add_door_sign("Utility", room, zval, room_id);}
	return 1;
}

bool building_t::add_furnace_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const floor_spacing(get_window_vspace());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-get_trim_thickness());
	place_area.z1() = zval;

	for (unsigned n = 0; n < 100; ++n) { // 100 tries
		cube_t furnace;
		bool dim(0), dir(0);
		if (!gen_furnace_cand(place_area, floor_spacing, 1, rgen, furnace, dim, dir)) break; // near_wall=1
		cube_t test_cube(furnace);
		test_cube.d[dim][dir] += (dir ? 1.0 : -1.0)*0.5*furnace.get_sz_dim(dim); // add clearance in front
		if (zval < ground_floor_z1) {test_cube.z2() = zval + get_floor_ceil_gap();} // basement furnace; extend to the ceiling to make sure there's space for the vent
		if (is_obj_placement_blocked(test_cube, room, 1) || overlaps_other_room_obj(test_cube, objs_start)) continue;
		if (zval > ground_floor_z1 && check_if_against_window(furnace, room, dim, !dir)) continue; // can fail for apartment and hotel utility rooms
		unsigned const flags((is_house ? RO_FLAG_IS_HOUSE : 0) | RO_FLAG_INTERIOR);
		interior->room_geom->objs.emplace_back(furnace, TYPE_FURNACE, room_id, dim, dir, flags, tot_light_amt);
		// add a blocker for clearance to avoid placing two furnaces together at a utility room corner
		cube_t blocker(test_cube);
		blocker.d[dim][!dir] = furnace.d[dim][dir];
		interior->room_geom->objs.emplace_back(blocker, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS);
		// no room for exhaust vent, so I guess it's inside the intake vent at the top
		return 1; // success/done
	} // for n
	return 0; // failed
}

// *** Parking Garages ***

vector3d building_t::get_parked_car_size() const {
	vector3d car_sz(get_nom_car_size());
	if (car_sz.x != 0.0) return car_sz; // valid car size, use this
	return get_window_vspace()*vector3d(1.67, 0.73, 0.45); // no cars, use size relative to building floor spacing
}

void add_pg_obstacles(vect_room_object_t const &objs, unsigned objs_start, unsigned objs_end, vect_cube_t &walls, vect_cube_t &beams, vect_cube_t &obstacles) {
	assert(objs_start <= objs_end);

	for (auto i = objs.begin()+objs_start; i != objs.begin()+objs_end; ++i) {
		if (i->type == TYPE_PG_WALL) {walls.push_back(*i);} // wall
		else if (i->type == TYPE_PG_PILLAR) { // pillar
			walls    .push_back(*i); // included in walls
			obstacles.push_back(*i); // pillars also count as obstacles
		}
		else if (i->type == TYPE_PG_BEAM) {beams    .push_back(*i);} // ceiling beam
		//else if (i->type == TYPE_RAMP   ) {obstacles.push_back(*i);} // ramps are obstacles for pipes, but they're already added
	} // for i
}

cube_t building_t::get_ext_basement_door_blocker() const {
	if (!has_ext_basement()) return cube_t();
	door_t const &eb_door(interior->get_ext_basement_door());
	cube_t avoid(eb_door.get_true_bcube());
	avoid.expand_in_dim( eb_door.dim, get_doorway_width ());
	avoid.expand_in_dim(!eb_door.dim, get_wall_thickness());
	return avoid;
}

void building_t::add_parking_garage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	unsigned num_floors, unsigned &nlights_x, unsigned &nlights_y, float &light_delta_z, light_ix_assign_t &light_ix_assign)
{
	assert(has_room_geom());
	rgen.rseed1 += 123*floor_ix; // make it unique per floor
	rgen.rseed2 += room_id;
	// rows are separated by walls and run in dim, with a road and parking spaces on either side of it;
	// spaces are arranged in !dim, with roads along the edges of the building that connect to the roads of each row
	bool const dim(room.dx() < room.dy()); // long/primary dim; cars are lined up along this dim, oriented along the other dim
	vector3d const car_sz(get_parked_car_size()), parking_sz(1.1*car_sz.x, 1.4*car_sz.y, 1.5*car_sz.z); // space is somewhat larger than a car; car length:width = 2.3
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness());
	float const int_wall_thick(get_wall_thickness()), wall_thickness(1.2*int_wall_thick), wall_hc(0.5*wall_thickness); // thicker
	float const ceiling_z(zval + window_vspacing - floor_thickness); // Note: zval is at floor level, not at the bottom of the room
	float const pillar_width(0.5*car_sz.y), pillar_hwidth(0.5*pillar_width), beam_hwidth(0.5*pillar_hwidth), road_width(2.3*car_sz.y); // road wide enough for two cars
	float const wid_sz(room.get_sz_dim(dim)), len_sz(room.get_sz_dim(!dim)), wid_sz_spaces(wid_sz - 2.0*road_width);
	float const min_strip_sz(2.0*parking_sz.x + road_width + max(wall_thickness, pillar_width)); // road + parking spaces on each side + wall/pillar
	assert(car_sz.z < (window_vspacing - floor_thickness)); // sanity check; may fail for some user parameters, but it's unclear what we do in that case
	unsigned const num_space_wid(wid_sz_spaces/parking_sz.y), num_full_strips(max(1U, unsigned(len_sz/min_strip_sz))); // take the floor
	bool const half_strip((num_full_strips*min_strip_sz + parking_sz.x + road_width + wall_thickness) < len_sz); // no space for a full row, add a half row
	bool const half_row_side(half_strip ? rgen.rand_bool() : 0); // pick a random side
	unsigned const num_rows(2*num_full_strips + half_strip), num_strips(num_full_strips + half_strip), num_walls(num_strips - 1);
	unsigned const capacity(num_rows*num_space_wid); // ignoring space blocked by stairs and elevators
	unsigned &nlights_len(dim ? nlights_x : nlights_y), &nlights_wid(dim ? nlights_y : nlights_x);
	nlights_len = num_rows; // lights over each row of parking spaces
	nlights_wid = round_fp(0.25*wid_sz/parking_sz.y); // 4 parking spaces per light on average, including roads
	//cout << TXT(nlights_len) << TXT(nlights_wid) << TXT(num_space_wid) << TXT(num_rows) << TXT(capacity) << endl; // TESTING
	assert(num_space_wid >= 4); // must fit at least 4 cars per row
	
	// add walls and pillars between strips
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());
	colorRGBA const wall_color(WHITE);
	cube_t room_floor_cube(room), virt_room_for_wall(room);
	set_cube_zvals(room_floor_cube, zval, ceiling_z);
	cube_t wall(room_floor_cube), pillar(room_floor_cube), beam(room_floor_cube);
	wall.expand_in_dim(dim, -road_width); // wall ends at roads that line the sides of the room; include pillar for better occluder and in case the pillar is skipped
	assert(wall.is_strictly_normalized());
	float wall_spacing(len_sz/(num_walls + 1));
	float const pillar_shift(0.01*pillar_width); // small value to avoid z-fighting
	float const wall_len(wall.get_sz_dim(dim) + 2.0f*pillar_shift), pillar_start(wall.d[dim][0] + pillar_hwidth - pillar_shift);
	float const row_width(wall_spacing - wall_thickness), space_length(0.5f*(row_width - road_width)), beam_spacing(len_sz/num_rows);
	unsigned const num_pillars(max(2U, unsigned(round_fp(0.25*wall_len/parking_sz.y)))); // every 4 spaces, at least 2 at the ends of the wall
	float const pillar_spacing((wall_len - pillar_width)/(num_pillars - 1)), beam_delta_z(0.95*wall.dz()), tot_light_amt(room.light_intensity);
	bool short_sides[2] = {0,0};

	if (half_strip) {
		short_sides[half_row_side] = 1;
		virt_room_for_wall.d[!dim][half_row_side] += (half_row_side ? 1.0 : -1.0)*space_length;
		wall_spacing = virt_room_for_wall.get_sz_dim(!dim)/(num_walls + 1); // recalculate wall spacing
	}
	light_delta_z = beam_delta_z - wall.dz(); // negative
	beam.z1()    += beam_delta_z; // shift the bottom up to the ceiling
	vect_cube_t obstacles, obstacles_exp, obstacles_ps, wall_parts, temp;
	// get obstacles for walls with and without clearance; maybe later add entrance/exit ramps, etc.
	interior->get_stairs_and_elevators_bcubes_intersecting_cube(room_floor_cube, obstacles,         wall_thickness ); // with small clearance in front, for beams and pipes
	interior->get_stairs_and_elevators_bcubes_intersecting_cube(room_floor_cube, obstacles_exp, 0.9*window_vspacing); // with more  clearance in front, for walls and pillars

	if (has_ext_basement()) { // everything should avoid the extended basement door
		cube_t const &avoid(get_ext_basement_door_blocker());
		obstacles    .push_back(avoid);
		obstacles_exp.push_back(avoid);
		obstacles_ps .push_back(avoid);
	}
	cube_with_ix_t const &ramp(interior->pg_ramp);
	bool const is_top_floor(floor_ix+1 == num_floors);
	
	// add ramp if one was placed during floorplanning, before adding parking spaces
	// Note: lights can be very close to ramps, but I haven't actually seen them touch; do we need to check for and handle that case?
	if (!ramp.is_all_zeros()) {
		bool const dim(ramp.ix >> 1), dir(ramp.ix & 1), is_blocked(is_top_floor && interior->ignore_ramp_placement);
		cube_t rc(ramp); // ramp clipped to this parking garage floor
		set_cube_zvals(rc, zval, (zval + window_vspacing));
		unsigned const flags(is_blocked ? 0 : RO_FLAG_OPEN); // ramp is open if the top exit is open
		objs.emplace_back(rc, TYPE_RAMP, room_id, dim, dir, flags, tot_light_amt, SHAPE_ANGLED, wall_color);
		obstacles    .push_back(rc); // don't place parking spaces next to the ramp
		obstacles_exp.push_back(rc); // clip beams to ramp
		obstacles_exp.back().expand_in_dim( dim,      road_width); // keep entrance and exit areas clear of parking spaces, even the ones against exterior walls
		obstacles_exp.back().expand_in_dim(!dim, 0.75*road_width); // keep walls and pillars away from the sides of ramps
		obstacles_ps .push_back(obstacles_exp.back()); // also keep parking spaces away from ramp
		// add ramp railings
		bool const side(ramp.get_center_dim(!dim) < room.get_center_dim(!dim)); // which side of the ramp the railing is on (opposite the wall the ramp is against)
		float const railing_thickness(0.4*wall_thickness), ramp_length(rc.get_sz_dim(dim)), dir_sign(dir ? 1.0 : -1.0), side_sign(side ? 1.0 : -1.0), shorten_factor(0.35);
		cube_t railing(rc);
		railing.d[!dim][!side] = railing.d[!dim][side] - side_sign*railing_thickness;
		railing.z1() += 0.5*railing_thickness; // place bottom of bar along ramp/floor
		cube_t ramp_railing(railing);
		ramp_railing.d[dim][dir] -= dir_sign*shorten_factor*ramp_length; // shorten length to only the lower part
		ramp_railing.z2() -= shorten_factor*railing.dz(); // shorten height by the same amount to preserve the slope
		colorRGBA const railing_color(LT_GRAY);
		objs.emplace_back(ramp_railing, TYPE_RAILING, room_id, dim, dir, RO_FLAG_OPEN, tot_light_amt, SHAPE_CUBE, railing_color); // lower railing
		set_cube_zvals(railing, rc.z2(), (rc.z2() + window_vspacing));
		railing.translate_dim(!dim, side_sign*railing_thickness); // shift off the ramp and onto the ajdacent floor

		if (!is_top_floor) { // add side railing for lower level
			railing.d[dim][!dir] += dir_sign*shorten_factor*ramp_length; // shorten length to only the upper part
			objs.emplace_back(railing, TYPE_RAILING, room_id, dim, 0, (RO_FLAG_OPEN | RO_FLAG_TOS), tot_light_amt, SHAPE_CUBE, railing_color);
		}
		else if (!is_blocked) { // add upper railings at the top for the full length
			railing.translate_dim( dim, -0.5*dir_sign*railing_thickness); // shift down the ramp a bit
			objs.emplace_back(railing, TYPE_RAILING, room_id, dim, 0, (RO_FLAG_OPEN | RO_FLAG_TOS), tot_light_amt, SHAPE_CUBE, railing_color);
			cube_t back_railing(rc);
			set_cube_zvals(back_railing, railing.z1(), railing.z2());
			back_railing.translate_dim( dim, -dir_sign*railing_thickness); // shift onto the ajdacent floor
			back_railing.translate_dim(!dim, 0.5*side_sign*railing_thickness); // shift away from the exterior wall
			back_railing.d[dim][dir] = back_railing.d[dim][!dir] + dir_sign*railing_thickness;
			objs.emplace_back(back_railing, TYPE_RAILING, room_id, !dim, 0, (RO_FLAG_OPEN | RO_FLAG_TOS), tot_light_amt, SHAPE_CUBE, railing_color);
		}
	}
	// add equipment room for first elevator extending to parking garage, only on the lowest floor
	if (floor_ix == 0) {
		for (elevator_t const &e : interior->elevators) {
			if (e.z1() > zval || !room.contains_cube_xy(e)) continue;
			bool const dim(!e.dim), dir(rgen.rand_bool());
			float const dir_sign(dir ? 1.0 : -1.0), door_width(get_doorway_width()), room_len(max(1.25f*e.get_sz_dim(dim), 1.0f*e.get_sz_dim(!dim)));
			cube_t sub_room(e);
			
			if (e.adj_elevator_ix >= 0) {
				elevator_t const &e2(get_elevator(e.adj_elevator_ix));
				if (e2.z1() == e.z1()) {sub_room.union_with_cube(e2);} // extend to cover both elevators
			}
			if (sub_room.get_sz_dim(!dim) < (1.5*door_width + 2.0*int_wall_thick)) continue; // too narrow
			float const ceil_zval(wall.z2());
			set_cube_zvals(sub_room, zval, ceil_zval);
			sub_room.d[dim][!dir] = e.d[dim][dir]; // adjacent to the elevator
			sub_room.d[dim][ dir] = e.d[dim][dir] + dir_sign*room_len; // extend outward
			cube_t avoid(sub_room);
			avoid.d[dim][dir] += dir_sign*door_width; // add space for the door
			if (!room.contains_cube(avoid)) continue; // outside the parking garage
			if (has_bcube_int_no_adj(avoid, obstacles)) continue; // check other elevators, stairs, and ext basement door
			// add door
			cube_t door(sub_room);
			door.d[dim][0] = door.d[dim][1] = sub_room.d[dim][dir] - 0.5*dir_sign*int_wall_thick; // shrink to zero area at wall centerline
			set_wall_width(door, sub_room.get_center_dim(!dim), 0.5*door_width, !dim);
			door_t Door(door, dim, !dir, rgen.rand_bool());
			Door.for_small_room = 1; // mark so that it only opens 90 degrees
			add_interior_door(Door, 0, 0, 1); // is_bathroom=0, make_unlocked=0, make_closed=1
			room_t equipment_room(room, sub_room); // keep flags, copy cube
			equipment_room.interior = 1; // treated as basement but not extended basement (no wall padding)
			equipment_room.expand_in_dim(!dim, -int_wall_thick); // subtract off walls
			equipment_room.d[dim][dir] -= dir_sign*int_wall_thick;
			// add a ceiling light; do we need a light switch as well?
			float const light_size(0.02*(sub_room.dx() + sub_room.dy()));
			cube_t light;
			set_light_xy(light, sub_room.get_cube_center(), light_size, !dim, SHAPE_CUBE);
			bool const recessed(fabs(light.d[dim][dir] - door.d[dim][!dir]) < door_width);
			float const light_dz((recessed ? 0.01 : 0.025)*window_vspacing);
			set_cube_zvals(light, (ceil_zval - light_dz), ceil_zval);
			colorRGBA const light_color(get_light_color_temp_range(0.2, 0.5, rgen));
			room_object_t light_obj(light, TYPE_LIGHT, room_id, !dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, light_color);
			add_sub_room_light(light_obj, equipment_room, !dim, objs.size(), light_ix_assign, rgen);
			// add machines
			add_machines_to_room(rgen, equipment_room, zval, room_id, tot_light_amt, objs.size(), 1); // objs_start at end; less_clearance=1
			// add a breaker panel
			//add_breaker_panel(rgen, bp, ceil_zval, dim, dir, room_id, tot_light_amt); // TODO

			// add walls, one on each side of elevator, and one on each side of the door; added last since these are occluders; what about back of room/side of elevator?
			if (interior->room_geom->wall_ps_start == 0) {interior->room_geom->wall_ps_start = objs.size();} // set if not set, on first level

			for (unsigned d = 0; d < 2; ++d) {
				cube_t side_wall(sub_room), door_wall(sub_room);
				door_wall.d[!dim][ d  ] = side_wall.d[!dim][!d] = side_wall.d[!dim][d] - (d ? 1.0 : -1.0)*int_wall_thick;
				door_wall.d[!dim][!d  ] = door.d[!dim][d];
				door_wall.d[ dim][!dir] = sub_room.d[dim][dir] - dir_sign*int_wall_thick;
				objs.emplace_back(side_wall, TYPE_PG_WALL, room_id, !dim, d,   0, tot_light_amt, SHAPE_CUBE, wall_color);
				objs.emplace_back(door_wall, TYPE_PG_WALL, room_id,  dim, dir, 0, tot_light_amt, SHAPE_CUBE, wall_color);
			} // d
			obstacles    .push_back(avoid);
			obstacles_exp.push_back(avoid);
			obstacles_ps .push_back(avoid);
			interior->elevator_equip_room = sub_room; // needed for parking garage light placement, and may be useful in other situations
			break; // done
		} // for e
	}
	// add walls and pillars
	bool const no_sep_wall(num_walls == 0 || (capacity < 100 && (room_id & 1))); // use room_id rather than rgen so that this agrees between floors
	bool const split_sep_wall(!no_sep_wall && (num_pillars >= 5 || (num_pillars >= 4 && rgen.rand_bool())));
	float sp_const(0.0);
	if      (no_sep_wall)           {sp_const = 0.25;} // no separator wall, minimal clearance around stairs
	else if (pri_hall_stairs_to_pg) {sp_const = 0.50;} // central stairs should open along the wall, not a tight space, need less clearance
	else if (split_sep_wall)        {sp_const = 0.75;} // gap should provide access, need slightly less clearance
	else                            {sp_const = 1.00;} // stairs may cut through/along full wall, need max clearance
	float const space_clearance(sp_const*max(0.5f*window_vspacing, parking_sz.y));
	// obstacles with clearance all sides, for parking spaces
	interior->get_stairs_and_elevators_bcubes_intersecting_cube(room_floor_cube, obstacles_ps, max(space_clearance, 0.9f*window_vspacing), space_clearance);
	if (interior->room_geom->wall_ps_start == 0) {interior->room_geom->wall_ps_start = objs.size();} // set if not set, on first level
	float center_pos(wall.get_center_dim(dim));
	// if there's an odd number of pillars, move the gap between two pillars on one side or the other
	if (split_sep_wall && (num_pillars & 1)) {center_pos += (rgen.rand_bool() ? -1.0 : 1.0)*0.5*pillar_spacing;}
	vect_cube_t pillars; // added after wall segments

	for (unsigned n = 0; n < num_walls+2; ++n) { // includes room far walls
		if (n < num_walls) { // interior wall
			float const pos(virt_room_for_wall.d[!dim][0] + (n + 1)*wall_spacing); // reference from the room far wall, assuming we can fit a full width double row strip
			set_wall_width(wall,   pos, wall_hc, !dim);
			set_wall_width(pillar, pos, pillar_hwidth, !dim);
			
			if (!no_sep_wall) {
				cube_t walls[2] = {wall, wall};

				if (split_sep_wall) { // add a gap between the walls for people to walk through
					walls[0].d[dim][1] = center_pos - 0.4*window_vspacing;
					walls[1].d[dim][0] = center_pos + 0.4*window_vspacing;
				}
				for (unsigned side = 0; side < (split_sep_wall ? 2U : 1U); ++side) {
					subtract_cubes_from_cube_split_in_dim(walls[side], obstacles_exp, wall_parts, temp, dim);
			
					for (auto const &w : wall_parts) {
						if (w.get_sz_dim(dim) < 2.0*window_vspacing) continue; // too short, skip
						objs.emplace_back(w, TYPE_PG_WALL, room_id, !dim, 0, 0, tot_light_amt, SHAPE_CUBE, wall_color);
						interior->room_geom->pgbr_walls[!dim].push_back(w); // save for occlusion culling
					}
				} // for side
			}
		}
		else { // room wall
			bool const side(n == num_walls+1);
			pillar.d[!dim][ side] = room.d[!dim][side];
			pillar.d[!dim][!side] = room.d[!dim][side] + (side ? -1.0 : 1.0)*pillar_hwidth; // half the width of an interior wall pillar
		}
		for (unsigned p = 0; p < num_pillars; ++p) { // add support pillars
			float const ppos(pillar_start + p*pillar_spacing);
			set_wall_width(pillar, ppos, pillar_hwidth, dim);
			if (has_bcube_int_xy(pillar, obstacles_exp)) continue; // skip entire pillar if it intersects stairs or an elevator
			pillars.push_back(pillar);
		} // for p
	} // for n
	for (auto const &p : pillars) {objs.emplace_back(p, TYPE_PG_PILLAR, room_id, !dim, 0, 0, tot_light_amt, SHAPE_CUBE, wall_color);}

	// add beams in !dim, at and between pillars
	unsigned const beam_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING);

	for (unsigned p = 0; p < (4*(num_pillars - 1) + 1); ++p) { // add beams, 4 per pillar
		float const ppos(pillar_start + 0.25*p*pillar_spacing);
		set_wall_width(beam, ppos, beam_hwidth, dim);
		subtract_cubes_from_cube_split_in_dim(beam, obstacles, wall_parts, temp, !dim);
		
		for (auto const &w : wall_parts) {
			if (min(w.dx(), w.dy()) > beam_hwidth) {objs.emplace_back(w, TYPE_PG_BEAM, room_id, !dim, 0, beam_flags, tot_light_amt, SHAPE_CUBE, wall_color);}
		}
	} // for p
	// add beams in dim for each row of lights
	for (unsigned n = 0; n < num_rows; ++n) {
		float const pos(room.d[!dim][0] + (n + 0.5)*beam_spacing);
		cube_t beam(room_floor_cube);
		beam.z1() += beam_delta_z; // shift the bottom up to the ceiling
		set_wall_width(beam, pos, beam_hwidth, !dim);
		subtract_cubes_from_cube_split_in_dim(beam, obstacles, wall_parts, temp, dim);

		for (auto const &w : wall_parts) {
			if (min(w.dx(), w.dy()) > beam_hwidth) {objs.emplace_back(w, TYPE_PG_BEAM, room_id, dim, 0, beam_flags, tot_light_amt, SHAPE_CUBE, wall_color);}
		}
	}

	// add parking spaces on both sides of each row (one side if half row)
	cube_t row(wall); // same length as the wall; includes the width of the pillars
	row.z2() = row.z1() + get_rug_thickness(); // slightly above the floor
	float const space_width(row.get_sz_dim(dim)/num_space_wid), strips_start(virt_room_for_wall.d[!dim][0]), wall_half_gap(2.0*wall_hc), space_shrink(row_width - space_length);
	bool const add_cars(enable_parked_cars() && !is_rotated()); // skip cars for rotated buildings
	unsigned const max_handicap_spots(capacity/20 + 1);
	unsigned num_handicap_spots(0);
	rand_gen_t rgen2(rgen); // make a copy to use with cars so that enabling cars doesn't change the parking garage layout

	for (unsigned n = 0; n < num_strips; ++n) {
		assert(space_length > 0.0);
		assert(space_width  > 0.0);
		assert(wall_spacing > 2.0*wall_half_gap);
		row.d[!dim][0] = strips_start + (n + 0)*wall_spacing + wall_half_gap;
		row.d[!dim][1] = strips_start + (n + 1)*wall_spacing - wall_half_gap;
		assert(row.is_strictly_normalized());
		if (row.get_sz_dim(!dim) < 1.2*space_shrink) continue; // space too small, likely due to incorrect building scale relative to car size; skip parking space placement

		for (unsigned d = 0; d < 2; ++d) { // for each side of the row
			bool const at_ext_wall[2] = {(n == 0 && d == 0), (n+1 == num_strips && d == 1)}, at_either_ext_wall(at_ext_wall[0] || at_ext_wall[1]);
			if ((short_sides[0] && at_ext_wall[0]) || (short_sides[1] && at_ext_wall[1])) continue; // skip this row
			float row_left_edge(row.d[dim][0]); // spaces start flush with the row, or flush with the room if this is the exterior wall
			unsigned num_spaces_per_row(num_space_wid);

			if (at_either_ext_wall) { // at either room exterior wall - can extend spaces up to the wall
				float row_right_edge(row.d[dim][1]); // opposite end of the row
				while ((row_left_edge  - space_width) > room.d[dim][0]) {row_left_edge  -= space_width; ++num_spaces_per_row;} // add rows to the left
				while ((row_right_edge + space_width) < room.d[dim][1]) {row_right_edge += space_width; ++num_spaces_per_row;} // add rows to the right
			}
			float const d_sign(d ? 1.0 : -1.0);
			cube_t space(row);
			space.d[!dim][!d] += d_sign*space_shrink; // shrink
			space.d[ dim][0]   = row_left_edge;
			bool last_was_space(0);

			for (unsigned s = 0; s < num_spaces_per_row; ++s) {
				space.d[dim][1] = space.d[dim][0] + space_width; // set width
				assert(space.is_strictly_normalized());
				
				if (has_bcube_int_xy(space, obstacles_ps)) { // skip entire space if it intersects stairs or an elevator
					if (last_was_space) {objs.back().flags &= ~RO_FLAG_ADJ_HI;} // no space to the right for the previous space
					last_was_space = 0;
				}
				else {
					unsigned flags(RO_FLAG_NOCOLL);
					if (last_was_space          ) {flags |= RO_FLAG_ADJ_LO;} // adjacent space to the left
					if (s+1 < num_spaces_per_row) {flags |= RO_FLAG_ADJ_HI;} // not the last space - assume there will be a space to the right
					bool const add_car(add_cars && rgen2.rand_float() < 0.5); // 50% populated with cars

					// make it a handicap spot if near an elevator and there aren't already too many
					if (num_handicap_spots < max_handicap_spots) {
						cube_t hc_area(space);
						hc_area.expand_by(1.5*space_width);
						if (!no_sep_wall) {hc_area.intersect_with_cube_xy(row);} // keep within the current row if there are walls in between rows

						for (elevator_t const &e : interior->elevators) {
							if (e.z1() > space.z2()) continue; // doesn't extend down to this level
							if (e.intersects_xy(hc_area)) {flags |= RO_FLAG_IS_ACTIVE; ++num_handicap_spots; break;}
						}
					}
					room_object_t pspace(space, TYPE_PARK_SPACE, room_id, !dim, d, flags, tot_light_amt, SHAPE_CUBE, wall_color); // floor_color?

					if (add_car) { // add a collider to block this area from the player, people, and rats; add first so that objs.back() is correct for the next iter
						car_t car(car_from_parking_space(pspace));
						objs.emplace_back(car.bcube, TYPE_COLLIDER, room_id, !dim, d, (RO_FLAG_INVIS | RO_FLAG_FOR_CAR));
						pspace.obj_id = (uint16_t)(objs.size() + rgen2.rand()); // will be used for the car model and color
						pspace.flags |= RO_FLAG_USED;
						objs.back().obj_id = pspace.obj_id; // will be used for loot collected from the car
					}
					if (no_sep_wall && !at_either_ext_wall) { // add small yellow curbs to block cars
						float const curb_height(0.04*window_vspacing), curb_width(1.5*curb_height);
						cube_t curb(space);
						curb.z2() += curb_height; // set height
						curb.d[!dim][!d] += d_sign*(space.get_sz_dim(!dim) - curb_width);   // shrink to the correct width
						curb.translate_dim(!dim, -d_sign*(0.5*pillar_hwidth + curb_width)); // move inward to avoid pillars and walls
						curb.expand_in_dim(dim, -0.2*space_width);
						objs.emplace_back(curb, TYPE_CURB, room_id, dim, 0, 0, 1.0, SHAPE_CUBE, colorRGBA(1.0, 0.8, 0.3)); // dir=0
					}
					objs.push_back(pspace);
					last_was_space = 1;
				}
				space.d[dim][0] = space.d[dim][1]; // shift to next space
			} // for s
		} // for d
	} // for n
	if (is_top_floor) { // place pipes on the top level parking garage ceiling (except for sprinkler pipes, which go on every floor)
		// avoid intersecting lights, pillars, walls, stairs, elevators, and ramps;
		// note that lights haven't been added yet though, but they're placed on beams, so we can avoid beams instead
		vect_cube_t walls, beams;
		add_pg_obstacles(objs, objs_start, objs.size(), walls, beams, obstacles);
		add_basement_electrical(obstacles, walls, beams, room_id, rgen);
		// get pipe ends (risers) coming in through the ceiling
		vect_riser_pos_t sewer, cold_water, hot_water, gas_pipes;
		get_pipe_basement_water_connections(sewer, cold_water, hot_water, rgen);
		vect_cube_t pipe_cubes;
		// hang sewer pipes under the ceiling beams; hang water pipes from the ceiling, above sewer pipes and through the beams
		float const ceil_zval(beam.z1()), water_ceil_zval(beam.z2());
		add_basement_pipes(obstacles, walls, beams, sewer,      pipe_cubes, room_id, num_floors, objs_start, ceil_zval,      rgen, PIPE_TYPE_SEWER, 0); // sewer
		add_to_and_clear(pipe_cubes, obstacles); // add sewer pipes to obstacles
		add_basement_pipes(obstacles, walls, beams, cold_water, pipe_cubes, room_id, num_floors, objs_start, water_ceil_zval, rgen, PIPE_TYPE_CW,   0); // cold water
		add_to_and_clear(pipe_cubes, obstacles); // add cold water pipes to obstacles
		add_basement_pipes(obstacles, walls, beams, hot_water,  pipe_cubes, room_id, num_floors, objs_start, water_ceil_zval, rgen, PIPE_TYPE_HW,   1); // hot water
		add_to_and_clear(pipe_cubes, obstacles); // add hot water pipes to obstacles
		get_pipe_basement_gas_connections(gas_pipes);
		add_basement_pipes(obstacles, walls, beams, gas_pipes,  pipe_cubes, room_id, num_floors, objs_start, water_ceil_zval, rgen, PIPE_TYPE_GAS,  1); // gas
		add_to_and_clear(pipe_cubes, obstacles); // add gas pipes to obstacles

		// if there are multiple parking garage floors, lights have already been added on the floor(s) below; add them as occluders for sprinkler pipes;
		// lights on the top floor will be added later and will check for pipe intersections; elevator equipment room lights are ignored as the entire room is an obstacle
		for (auto i = objs.begin()+interior->room_geom->wall_ps_start; i < objs.begin()+objs_start; ++i) {
			if (i->type == TYPE_LIGHT && room.contains_cube(*i)) {obstacles.push_back(*i);}
		}
		add_sprinkler_pipes(obstacles, walls, beams, pipe_cubes, room_id, num_floors, objs_start, rgen);
	}
}

// *** Electrical ***

bool is_good_conduit_placement(cube_t const &c, cube_t const &avoid, vect_room_object_t const &objs, unsigned end_ix, unsigned skip_ix) {
	if (!avoid.is_all_zeros() && c.intersects(avoid)) return 0;

	for (unsigned i = 0; i < end_ix; ++i) {
		if (i != skip_ix && objs[i].intersects(c)) return 0;
	}
	return 1;
}
void building_t::add_basement_electrical(vect_cube_t &obstacles, vect_cube_t const &walls, vect_cube_t const &beams, int room_id, rand_gen_t &rgen) {
	float const tot_light_amt = 1.0; // to offset the darkness of the basement
	cube_t const &basement(get_basement());
	float const floor_spacing(get_window_vspace()), fc_thickness(get_fc_thickness()), floor_height(floor_spacing - 2.0*fc_thickness), ceil_zval(basement.z2() - fc_thickness);
	unsigned const num_panels(is_house ? 1 : (2 + (rgen.rand()&3))); // 1 for houses, 2-4 for office buildings
	auto &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());

	for (unsigned n = 0; n < num_panels; ++n) {
		float const bp_hwidth(rgen.rand_uniform(0.15, 0.25)*(is_house ? 0.7 : 1.0)*floor_height), bp_hdepth(rgen.rand_uniform(0.05, 0.07)*(is_house ? 0.5 : 1.0)*floor_height);
		if (bp_hwidth > 0.25*min(basement.dx(), basement.dy())) continue; // basement too small
		cube_t c;
		set_cube_zvals(c, (ceil_zval - 0.75*floor_height), (ceil_zval - rgen.rand_uniform(0.2, 0.35)*floor_height));

		for (unsigned t = 0; t < ((n == 0) ? 100U : 1U); ++t) { // 100 tries for the first panel, one try after that
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
			float const dir_sign(dir ? -1.0 : 1.0), wall_pos(basement.d[dim][dir]), bp_center(rgen.rand_uniform(basement.d[!dim][0]+bp_hwidth, basement.d[!dim][1]-bp_hwidth));
			set_wall_width(c, bp_center, bp_hwidth, !dim);
			c.d[dim][ dir] = wall_pos;
			c.d[dim][!dir] = wall_pos + dir_sign*2.0*bp_hdepth; // extend outward from wall
			assert(c.is_strictly_normalized());
			cube_t test_cube(c);
			test_cube.d[dim][!dir] += dir_sign*2.0*bp_hwidth; // add a width worth of clearance in the front so that the door can be opened
			test_cube.z2() = ceil_zval; // extend up to ceiling to ensure space for the conduit
			if (has_bcube_int(test_cube, obstacles) || has_bcube_int(test_cube, walls) || has_bcube_int(test_cube, beams)) continue; // bad breaker box or conduit position
			if (is_cube_close_to_doorway(test_cube, basement, 0.0, 1)) continue; // needed for ext basement doorways; inc_open=1
			unsigned cur_room_id((room_id < 0) ? get_room_containing_pt(c.get_cube_center()) : (unsigned)room_id); // calculate room_id if needed
			add_breaker_panel(rgen, c, ceil_zval, dim, dir, cur_room_id, tot_light_amt);
			cube_t blocker(c);
			set_cube_zvals(blocker, ceil_zval-floor_height, ceil_zval); // expand to floor-to-ceiling
			obstacles.push_back(blocker); // block off from pipes

			if (is_house) { // try to reroute outlet conduits that were previously placed on the same wall horizontally to the breaker box
				float const conn_height(c.z1() + rgen.rand_uniform(0.25, 0.75)*c.dz());
				cube_t avoid;
				if (has_ext_basement()) {avoid = get_ext_basement_door_blocker();}

				// Note: only need to check basement objects, but there's no easy way to do this (index not recorded), so we check them all
				for (unsigned i = 0; i < objs_start; ++i) { // can't use an iterator as it may be invalidated
					room_object_t &obj(objs[i]);

					if (obj.type == TYPE_FURNACE) { // connect to furnace if on the same wall; straight segment with no connectors
						cube_t conn;
						if (!connect_furnace_to_breaker_panel(obj, c, dim, dir, conn_height, conn)) continue;
						if (!is_good_conduit_placement(conn, avoid, objs, objs_start, i)) continue;
						objs.emplace_back(conn, TYPE_PIPE, room_id, !dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, LT_GRAY); // horizontal
						continue;
					}
					// handle power outlets
					if (obj.type != TYPE_PIPE || obj.shape != SHAPE_CYLIN || obj.dim != 0 || obj.dir != 1) continue; // vertical pipes only
					if (obj.d[dim][dir] != wall_pos) continue; // wrong wall; though we could always bend the pipe around the corner?
					if (obj.intersects_xy(c))        continue; // too close to the panel; may intersect anyway; is this possibe?
					float const radius(obj.get_radius());
					if ((obj.z1() + 4.0*radius) > conn_height) continue; // too high (possibly light switch rather than outlet)
					float const pipe_center(obj.get_center_dim(!dim));
					bool const pipe_dir(pipe_center < bp_center);
					float const bp_near_edge(c.d[!dim][!pipe_dir]);
					room_object_t h_pipe(obj);
					set_wall_width(h_pipe, conn_height, radius, 2); // set zvals
					h_pipe.d[!dim][!pipe_dir] = pipe_center;
					h_pipe.d[!dim][ pipe_dir] = bp_near_edge;
					assert(h_pipe.is_strictly_normalized());
					// check for valid placement/blocked by ext basement door or other objects;
					// we may place two horizontal conduits that overlap, but they should be at the same positions and look correct enough
					if (!is_good_conduit_placement(h_pipe, avoid, objs, objs_start, i)) continue;
					h_pipe.dim = !dim;
					h_pipe.dir = 0; // horizontal
					obj.z2()   = conn_height; // shorten original conduit

					if (min(obj.dz(), h_pipe.get_sz_dim(!dim)) > 4.0*radius) { // not a short segment; add connector parts
						float const conn_exp(0.2*radius), xlate((dir ? -1.0 : 1.0)*conn_exp);
						float const conn_len(2.0*radius), signed_conn_len((pipe_dir ? 1.0 : -1.0)*conn_len);
						obj   .translate_dim(dim, xlate); // move away from the wall
						h_pipe.translate_dim(dim, xlate);
						room_object_t v_conn(obj), h_conn(h_pipe);
						h_conn.expand_in_dim(dim,  conn_exp);
						v_conn.expand_in_dim(dim,  conn_exp);
						h_conn.expand_in_dim(2,    conn_exp);
						v_conn.expand_in_dim(!dim, conn_exp);
						room_object_t v_conn2(v_conn), h_conn2(h_conn); // these will be the short connectors going to the outlet and breaker box, respectively
						v_conn .z1()  = conn_height -     conn_len;
						v_conn2.z2()  = obj.z1()    + 0.6*conn_len;
						h_conn .d[!dim][ pipe_dir] = pipe_center  +     signed_conn_len;
						h_conn2.d[!dim][!pipe_dir] = bp_near_edge - 0.6*signed_conn_len;
						v_conn .flags |=  RO_FLAG_ADJ_HI; // add joint connecting to horizontal pipe
						h_conn .flags |= (RO_FLAG_HANGING | (pipe_dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO)); // one flat end
						v_conn2.flags |= (RO_FLAG_HANGING | RO_FLAG_ADJ_HI); // flat end on top
						h_conn2.flags |= (RO_FLAG_HANGING | (pipe_dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI)); // one flat end
						v_conn .color  = h_conn.color = v_conn2.color  = h_conn2.color = GRAY; // darker
						objs.push_back(v_conn );
						objs.push_back(h_conn );
						objs.push_back(v_conn2);
						objs.push_back(h_conn2);
						// duplicate to get flat ends on the vertical part as well
						v_conn.flags |=  RO_FLAG_ADJ_LO | RO_FLAG_HANGING;
						v_conn.flags &= ~RO_FLAG_ADJ_HI;
						objs.push_back(v_conn);
					}
					else { // short segment
						obj.flags |= RO_FLAG_ADJ_HI; // add joint connecting to horizontal pipe
					}
					objs.push_back(h_pipe); // not added to obstacles since this should not affect pipes routed in the ceiling
				} // for i
			}
			break; // done
		} // for t
	} // for n
}

void building_t::add_basement_electrical_house(rand_gen_t &rgen) {
	assert(has_room_geom());
	float const obj_expand(0.5*get_wall_thickness());
	cube_t const &basement(get_basement());
	vect_cube_t obstacles, walls;

	for (unsigned d = 0; d < 2; ++d) { // add basement walls
		for (cube_t const &wall : interior->walls[d]) {
			if (wall.zc() < ground_floor_z1 && basement.intersects(wall)) {walls.push_back(wall);}
		}
	}
	// add basement objects; include them all, since it's not perf critical; we haven't added objects such as trim yet;
	// what about blocking the breaker box with something, is that possible?
	for (room_object_t const &c : interior->room_geom->objs) {
		if (c.z1() < ground_floor_z1 && basement.intersects(c)) {obstacles.push_back(c); obstacles.back().expand_by(obj_expand);} // with some clearance
	}
	vector_add_to(interior->stairwells, obstacles); // add stairs; what about open doors?
	vector_add_to(interior->elevators,  obstacles); // there probably are none, but it doesn't hurt to add them
	add_basement_electrical(obstacles, walls, vect_cube_t(), -1, rgen); // no beams, room_id=-1 (to be calculated)
}

void building_t::add_parking_garage_ramp(rand_gen_t &rgen) {
	assert(interior && !is_house && has_parking_garage);
	cube_with_ix_t &ramp(interior->pg_ramp);
	assert(ramp.is_all_zeros()); // must not have been set
	cube_t const &basement(get_basement());
	bool const dim(basement.dx() < basement.dy()); // long/primary dim
	// see building_t::add_parking_garage_objs(); make sure there's space for a ramp plus both exit dirs within the building width
	float const width(basement.get_sz_dim(!dim)), road_width(min(0.25f*width, 2.3f*get_parked_car_size().y));
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const z1(basement.z1() + fc_thick), z2(basement.z2() + fc_thick); // bottom level basement floor to first floor floor
	bool const ramp_pref_xdir(rgen.rand_bool()), ramp_pref_ydir(rgen.rand_bool());
	bool added_ramp(0), dir(0);

	for (unsigned pass = 0; pass < 2 && !added_ramp; ++pass) {
		for (unsigned xd = 0; xd < 2 && !added_ramp; ++xd) {
			for (unsigned yd = 0; yd < 2; ++yd) {
				bool const xdir(bool(xd) ^ ramp_pref_xdir), ydir(bool(yd) ^ ramp_pref_ydir);
				float const xsz((dim ? 2.0 : 1.0)*road_width), ysz((dim ? 1.0 : 2.0)*road_width); // longer in !dim
				unsigned const num_ext(unsigned(basement.d[0][xdir] == bcube.d[0][xdir]) + unsigned(basement.d[1][ydir] == bcube.d[1][ydir]));
				if (num_ext < 2-pass) continue; // must be on the exterior edge of the building in both dims for pass 0, and one dim for pass 1
				dir = (dim ? xdir : ydir);
				point corner(basement.d[0][xdir], basement.d[1][ydir], z1);
				corner[!dim] += (dir ? -1.0 : 1.0)*road_width; // shift away from the wall so that cars have space to turn onto the level floor
				point const c1((corner.x - 0.001*(xdir ? 1.0 : -1.0)*xsz), (corner.y - 0.001*(ydir ? 1.0 : -1.0)*ysz), z1); // slight inward shift to prevent z-fighting
				point const c2((corner.x + (xdir ? -1.0 : 1.0)*xsz), (corner.y + (ydir ? -1.0 : 1.0)*ysz), z2);
				cube_t const ramp_cand(c1, c2);
				assert(ramp_cand.is_strictly_normalized());
				cube_t test_cube(ramp_cand);
				test_cube.expand_in_dim(!dim, road_width); // extend outward for clearance to enter/exit the ramp (ramp dim is actually !dim)
				if (interior->is_blocked_by_stairs_or_elevator(test_cube)) continue;
				// check for backrooms door, in case it was placed already (but currently it's not)
				if (has_ext_basement() && get_ext_basement_door_blocker().intersects(test_cube)) continue; // blocked by extended basement door
				ramp = cube_with_ix_t(ramp_cand, (((!dim)<<1) + dir)); // encode dim and dir in ramp index field
				added_ramp = 1;
				break; // done
			} // for yd
		} // for xd
	} // for pass
	if (!added_ramp) return; // what if none of the 4 corners work for a ramp?
	// add landings, which are used to draw the vertical edges of the cutout
	unsigned num_floors(calc_num_floors(basement, window_vspacing, floor_thickness));
	float z(basement.z1() + window_vspacing); // start at upper floor rather than lower floor

	if (1) { // FIXME: rooms on the ground floor above ramps aren't yet handled, so clip ramps to avoid disrupting their floors until this is fixed
		ramp.z2() -= 2.0*floor_thickness;
		--num_floors;
		interior->ignore_ramp_placement = 1; // okay to place room objects over ramps because the floor has not been removed
	}
	for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) { // skip first floor - draw pairs of floors and ceilings
		landing_t landing(ramp, 0, f, !dim, dir, 0, SHAPE_RAMP, 0, (f+1 == num_floors), 0, 1); // for_ramp=1
		set_cube_zvals(landing, (z - fc_thick), (z + fc_thick));
		interior->landings.push_back(landing);
	}
	// cut out spaces from floors and ceilings
	subtract_cube_from_floor_ceil(ramp, interior->floors  );
	subtract_cube_from_floor_ceil(ramp, interior->ceilings);
	// make rooms over the ramp of type RTYPE_RAMP_EXIT
}

