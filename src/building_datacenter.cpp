// 3D World - Data Center Buildings
// by Frank Gennari 5/25/26

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

extern object_model_loader_t building_obj_model_loader;

int select_dc_battery_model();


// returns the hallway, if one was created
cube_t building_t::create_datacenter_floorplan(unsigned part_id, float window_hspacing[2], int num_windows_per_side[2], rand_gen_t &rgen) {
	// add hallway; two server rooms to either side, office(s) at one end,
	// power/AC utility room(s) at the other, plus bathroom, with stairs and elevator to the side of the hallway
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const doorway_width(get_nominal_doorway_width()), doorway_hwidth(0.5*doorway_width);
	float const wall_thick(get_wall_thickness()), wall_hthick(0.5*wall_thick), wall_edge_spacing(0.05*wall_thick);
	cube_t const &part(parts[part_id]);
	vector3d const psz(part.get_size());
	bool const min_dim(psz.y < psz.x), max_dim(!min_dim);
	vector<room_t> &rooms(interior->rooms);
	assert(rooms.empty()); // must call this first

	if (psz[min_dim] < 4.0*window_vspacing) { // too small, assign entire part to a single server room; shouldn't happen
		add_assigned_room(part, part_id, RTYPE_SERVER);
		return cube_t(); // no hallway
	}
	if (part_id == 0) {hallway_dim = max_dim;} // long dim
	float num_hall_windows, hall_width, room_width; // unused
	cube_t const hall(get_hallway_for_part(part, num_hall_windows, hall_width, room_width));
	auto &room_walls(interior->walls[max_dim]), &hall_walls(interior->walls[min_dim]); // room_walls: perpendicular to hallway; hall_walls: parallel to hallway
	bool se_end(rgen.rand_bool()); // also office end, opposite utility end, near front door

	// if there was a cooling tower placed on the roof, place stairs on the opposite end so that the utility room is closer to the cooling tower
	for (roof_obj_t const &d : details) {
		if (d.type != ROOF_OBJ_COOLING) continue;
		se_end = (d.get_center_dim(max_dim) < part.get_center_dim(max_dim));
		break; // only need one
	}
	// place stairs/elevator/offices on the side with more space (further from the part center), or random if centered
	float const hall_center(hall.get_center_dim(min_dim));
	bool const se_side(hall_center < part.get_center_dim(min_dim) + wall_thick*rgen.signed_rand_float());
	float const se_side_sign(se_side ? 1.0 : -1.0), se_end_sign(se_end ? 1.0 : -1.0), wind_space(window_hspacing[max_dim]);
	unsigned const hall_walls_start(hall_walls.size());
	unsigned const num_short_wind(num_windows_per_side[min_dim]); // 6-11
	unsigned const num_long_wind (num_windows_per_side[max_dim]); // 6-12: {S,O,U}: 321 421 422 522 532 632 732
	unsigned const num_office_wind((num_long_wind >= 10) ? 3 : 2), num_util_wind((num_long_wind >= 8) ? 2 : 1), num_bath_wind(1);
	assert(num_office_wind + num_util_wind < num_long_wind);
	float const office_pos(part.d[max_dim][ se_end] - se_end_sign*num_office_wind*wind_space);
	float const util_pos  (part.d[max_dim][!se_end] + se_end_sign*num_util_wind  *wind_space);
	float const bath_pos  (office_pos + se_end_sign*num_bath_wind*wind_space); // between server room and office on side opposite stairs/elevator
	bool const skip_sr_util_windows = 1;
	interior->dc_info.reset(new bldg_datacenter_info_t(se_end, office_pos, util_pos, bath_pos, skip_sr_util_windows)); // needed for window creation, etc.

	if (num_short_wind >= 9 && num_long_wind >= 9) { // larger building, at least 4 windows per side; add secondary hallway and two rows of side rooms
		// is this needed?
	}
	// add hallway room
	unsigned const num_lights(num_long_wind); // likely will be reduced
	unsigned const hall_room_id(add_room(hall, part_id, num_lights, 1, 0)); // add primary hallway
	//rooms.back().mark_open_wall_dim(min_dim); // flag primary hallway as open on sides if there are secondary hallways
	if (part_id == 0) {pri_hall = hall;}
	// place stairs and elevator along the hallway
	float const se_wall_pos(hall.d[min_dim][se_side] - se_side_sign*wall_hthick);
	float const ewidth(1.8*doorway_width), edepth(1.8*doorway_width), stairs_width(2.5*doorway_width);
	float stairs_depth(window_hspacing[min_dim]); // one window width
	if (stairs_depth < 2.4*doorway_width) {stairs_depth *= 2.0;} // increase to 2 windows if needed
	float const stairs_start(part.d[max_dim][se_end] - 0.8*se_end_sign*wall_hthick);
	float const stairs_end(stairs_start - se_end_sign*stairs_width), elevator_start(stairs_end - 0.25*se_end_sign*wall_thick);
	cube_t stairs(part), elevator(part);
	stairs  .d[min_dim][!se_side] = elevator.d[min_dim][!se_side] = se_wall_pos;
	stairs  .d[min_dim][ se_side] = se_wall_pos + se_side_sign*stairs_depth;
	elevator.d[min_dim][ se_side] = se_wall_pos + se_side_sign*edepth;
	stairs  .d[max_dim][ se_end ] = stairs_start;
	stairs  .d[max_dim][!se_end ] = stairs_end;
	elevator.d[max_dim][ se_end ] = elevator_start; // small gap between stairs and elevator
	elevator.d[max_dim][!se_end ] = elevator_start - se_end_sign*ewidth;
	interior->stairwells.emplace_back(stairs, 0, min_dim, se_side, SHAPE_U); // add temp stairs so that we can extract these variables later
	interior->stairwells.back().against_wall[se_end] = 1;
	have_hall_side_stairs = 1;
	get_room(hall_room_id).has_stairs = 255; // stairs on all floors
	elevator_t E(elevator, hall_room_id, min_dim, !se_side, 0, 1); // elevator shaft; at_edge=0, interior_room=1 (considered interior-enough)
	add_or_extend_elevator(E, 1);
	get_room(hall_room_id).has_elevator = 1;
	// setup server room area
	cube_t server_area(part), office_area(part), util_area(part);
	server_area.d[max_dim][ se_end] = office_area.d[max_dim][!se_end] = office_pos;
	server_area.d[max_dim][!se_end] = util_area  .d[max_dim][ se_end] = util_pos;
	assert(server_area.is_strictly_normalized());
	assert(office_area.is_strictly_normalized());
	assert(util_area  .is_strictly_normalized());
	float const server_area_len(server_area.get_sz_dim(max_dim)), door_end_space(max(doorway_width, 0.25f*server_area_len));
	float const server_door_pos(rgen.rand_uniform(server_area.d[max_dim][0]+door_end_space, server_area.d[max_dim][1]-door_end_space));
	bool const conn_office_to_server(rgen.rand_float() < 0.8), conn_util_to_server(rgen.rand_float() < 0.4);
	unsigned server_room_ids[2] = {};

	for (unsigned d = 0; d < 2; ++d) { // each side of hallway
		float const hall_side(hall.d[min_dim][d]);
		// add server rooms
		server_room_ids[d] = rooms.size();
		cube_t server_room(server_area);
		server_room.d[min_dim][!d] = hall_side;
		add_assigned_room(server_room, part_id, RTYPE_SERVER);
		rooms.back().set_is_large(); // for AI navigation
		// add office and bathroom
		bool const br_side(bool(d) != se_side);
		cube_t office(office_area);
		office.d[min_dim][!d] = hall_side;
		float br_door_pos(0.0);
		
		if (br_side) { // bathroom side is opposite stairs and elevator
			cube_t bathroom(office);
			bathroom.d[max_dim][se_end] = office.d[max_dim][!se_end] = bath_pos;
			// only the bottom two floors are bathrooms (Men/Women); upper floors are offices
			add_assigned_room(bathroom, part_id, RTYPE_OFFICE); // set the offices first
			for (unsigned f = 0; f < 2; ++f) {rooms.back().assign_to(RTYPE_BATH, f, 1, 1);} // locked=1, force=1
			br_door_pos = bathroom.get_center_dim(max_dim); // centered
			// small office next to bathroom is security room on the first floor
			float const rlen(office.get_sz_dim(min_dim)), rwidth(office.get_sz_dim(max_dim));
			unsigned const num_windows_this_side(round_fp(office.get_sz_dim(min_dim)/window_hspacing[min_dim]));

			if (num_windows_this_side >= 2 && rlen > 4.0*window_vspacing && (rlen + rwidth) > 6.0*window_vspacing && rwidth > 2.0*doorway_width) {
				// if large, split into office and security room
				// need to shift wall pos to handle odd number of windows; use half the border since we have plenty of space
				float split_pos(office.d[min_dim][0] + rgen.rand_uniform(0.4, 0.6)*rlen);
				split_pos = shift_val_to_not_intersect_window(part, split_pos, window_hspacing[min_dim], 0.5*get_window_h_border(), min_dim);
				cube_t outer_office(office);
				outer_office.d[min_dim][d] = office.d[min_dim][!d] = split_pos;
				add_assigned_room(outer_office, part_id, RTYPE_OFFICE);
				cube_t walls[2] = {office, office};
				create_wall(walls[0], min_dim, split_pos, fc_thick, wall_hthick, wall_edge_spacing);
				float const door_pos(rgen.rand_uniform(office.d[max_dim][0]+1.5*doorway_hwidth, office.d[max_dim][1]-1.5*doorway_hwidth));
				remove_section_from_cube_and_add_door(walls[0], walls[1], (door_pos - doorway_hwidth), (door_pos + doorway_hwidth), max_dim, d, 0, 0, 1); // closed=1
				for (unsigned e = 0; e < 2; ++e) {hall_walls.push_back(walls[e]);}
			}
			add_assigned_room(office, part_id, RTYPE_OFFICE);
			rooms.back().assign_to(RTYPE_SECURITY, 0, 1, 1); // floor=0, locked=1, force=1
		}
		else { // stairs/elevator side; add a large office with stairs/elevator cut out
			add_assigned_room(office, part_id, RTYPE_OP_CENTER);
		}
		// add utility rooms
		cube_t utility(util_area);
		utility.d[min_dim][!d] = hall_side;
		add_assigned_room(utility, part_id, RTYPE_UTILITY);
		// add walls and doors along hallway
		cube_t walls[2] = {part, part}; // {lo, hi} sides
		create_wall(walls[0], min_dim, hall_side, fc_thick, wall_hthick, wall_edge_spacing);
		remove_section_from_cube_and_add_door(walls[0], walls[1], (server_door_pos - doorway_hwidth), (server_door_pos + doorway_hwidth), max_dim, d); // for security room
		
		if (br_side) { // add bathroom door
			cube_t br_wall;
			remove_section_from_cube_and_add_door(walls[se_end], br_wall, (br_door_pos - doorway_hwidth), (br_door_pos + doorway_hwidth), max_dim, d);
			if (se_end) {swap(walls[se_end], br_wall);} // swap lo/hi order
			hall_walls.push_back(br_wall);
		}
		for (unsigned e = 0; e < 2; ++e) { // for each end of the hallway
			bool const is_office_side(bool(e) == se_end);
			cube_t &adj_room(is_office_side ? office : utility);
			cube_t wall(walls[e]);
			if (wall.intersects(stairs)) {wall.d[max_dim][se_end] = elevator.d[max_dim][!se_end];} // shorten wall in front of stairs to end at edge of elevator
			float const wall_lo(max(wall.d[max_dim][0], adj_room.d[max_dim][0])), wall_hi(min(wall.d[max_dim][1], adj_room.d[max_dim][1]));
			assert(wall_lo < wall_hi);

			if ((wall_hi - wall_lo) > 1.2*doorway_width) { // not too narrow for a door; should always get here
				float const rgen_lo(wall_lo + doorway_hwidth + doorway_hwidth), rgen_hi(wall_hi - doorway_hwidth - doorway_hwidth);
				float const door_pos((rgen_lo < rgen_hi) ? rgen.rand_uniform(rgen_lo, rgen_hi) : 0.5*(rgen_lo + rgen_hi)); // centered if wall is too narrow
				cube_t wall2;
				remove_section_from_cube_and_add_door(wall, wall2, (door_pos - doorway_hwidth), (door_pos + doorway_hwidth), max_dim, d);
				hall_walls.push_back(wall2);
			}
			hall_walls.push_back(wall);
		} // for d
		// add walls and doors between rooms
		float const room_split_wall_pos[3] = {office_pos, util_pos, bath_pos};

		for (unsigned e = 0; e < (br_side ? 3U : 2U); ++e) { // each room to split with a wall
			cube_t rwall(part);
			rwall.d[min_dim][!d] = hall_side + (d ? 1.0 : -1.0)*wall_hthick;
			create_wall(rwall, max_dim, room_split_wall_pos[e], fc_thick, wall_hthick, wall_edge_spacing);

			if ((conn_office_to_server && e == 0 && !br_side) || (conn_util_to_server && e == 1)) { // add a door between this room and the server room
				bool const open_dir(bool(e) != se_end);
				float const lo(rwall.d[min_dim][0] + doorway_width), hi(rwall.d[min_dim][1] - doorway_width);

				if (lo < hi) { // should be true
					float const door_pos(rgen.rand_uniform(lo, hi));
					cube_t rwall2;
					remove_section_from_cube_and_add_door(rwall, rwall2, (door_pos - doorway_hwidth), (door_pos + doorway_hwidth), min_dim, open_dir);
					room_walls.push_back(rwall2);
				}
			}
			room_walls.push_back(rwall);
		} // for e
	} // for d
	// add interior windows to walls separating server rooms from hallway
	unsigned const hall_walls_end(hall_walls.size());

	for (unsigned w = hall_walls_start; w < hall_walls_end; ++w) {
		cube_t &wall(hall_walls[w]);
		if (!wall.intersects(hall) || !wall.intersects(server_area)) continue;
		cube_t wind_area(wall);
		wind_area.intersect_with_cube(server_area);
		float const wall_len(wind_area.get_sz_dim(max_dim)), window_edge_space(0.25*wall_len);
		if (wall_len < 2.0*window_vspacing) continue; // too short
		bool const wdir(hall_center < wall.get_center_dim(min_dim));
		unsigned const room_id(server_room_ids[wdir]);
		float const wind_lo(wind_area.d[max_dim][0] + window_edge_space), wind_hi(wind_area.d[max_dim][1] - window_edge_space);
		cube_t wall2(wall), window(wall);
		wall .d[max_dim][1] = window.d[max_dim][0] = wind_lo; // low  edge of window
		wall2.d[max_dim][0] = window.d[max_dim][1] = wind_hi; // high edge of window
		hall_walls.push_back(wall2);
		get_room(     room_id).set_interior_window();
		get_room(hall_room_id).set_interior_window();
		interior->int_windows.emplace_back(window, room_id);
	} // for w
	return hall;
}

void building_t::add_datacenter_outdoor_objs(rand_gen_t &rgen) {
	// add outdoor AC units
	if ((real_num_parts - has_basement()) != 1) return; // only handles a single part
	bool const dim(!bool(hallway_dim));
	cube_t const &part(get_first_part());
	float const floor_spacing(get_window_vspace()), bldg_length(part.get_sz_dim(!dim)), edge_space(rgen.rand_uniform(0.2, 0.3)*bldg_length);
	float const depth(0.40*floor_spacing), width(1.5*depth), height(0.35*floor_spacing);
	unsigned const num(min(16.0, rgen.rand_uniform(0.4, 0.6)*(bldg_length - 2.0*edge_space)/width));
	if (num < 2) return; // shouldn't happen
	float const start(part.d[!dim][0] + edge_space), end(part.d[!dim][1] - edge_space), spacing((end - start)/num);
	roof_obj_t ac(ROOF_OBJ_AC);
	set_cube_zvals(ac, part.z1(), (part.z1() + height));

	for (unsigned dir = 0; dir < 2; ++dir) {
		ac.d[dim][!dir] = part.d[dim][ dir] + (dir ? 1.0 : -1.0)*0.07*floor_spacing; // place slightly away from the exterior wall to avoid window sills
		ac.d[dim][ dir] = ac  .d[dim][!dir] + (dir ? 1.0 : -1.0)*depth;

		for (unsigned n = 0; n < num; ++n) {
			set_wall_width(ac, (start + n*spacing), 0.5*width, !dim);
			details.push_back(ac);
			union_with_coll_bcube(ac);
		}
	} // for dir
	has_ac = 1;
}

void building_t::add_datacenter_rooftop_objs(rand_gen_t &rgen) { // Note: interior hasn't been setup yet
	cube_t const &part(get_first_part()); // main part
	bool const dim(part.dx() < part.dy()); // long dim
	cube_t roof(part);
	roof.expand_by_xy(-get_roof_wall_thick()); // shrink off the outer roof wall
	cube_t ct_area(roof);
	ct_area.expand_by_xy(-0.05*roof.get_sz_dim(!dim)); // add a border
	float const roof_zval(roof.z2()), roof_width(ct_area.get_sz_dim(!dim)), roof_length(ct_area.get_sz_dim(dim));
	float const ct_width(rgen.rand_uniform(0.16, 0.22)*roof_width), ct_height(rgen.rand_uniform(0.7, 1.0)*ct_width);
	float const ct_len(min(0.3f*roof_length, rgen.rand_uniform(1.5, 1.75)*ct_width)), ct_hlen(0.5*ct_len), gap(0.25*ct_len);
	bool const side(rgen.rand_bool()); // choose a random side; the utility room will be added to this side
	bool const add_two(0.65*roof_length > 2.0*(ct_len + gap));
	ct_area.d[dim][side] -= (side ? 1.0 : -1.0)*(add_two ? 0.35 : 0.5)*roof_length; // shrink to 65%/50% of the area
	cube_t ct_areas[2] = {ct_area, ct_area}, blockers[2];
	
	if (add_two) { // split in half in dim
		float const split_pos(ct_area.get_center_dim(dim));
		ct_areas[0].d[dim][1] = split_pos - gap;
		ct_areas[1].d[dim][0] = split_pos + gap;
	}
	for (unsigned n = 0; n < (add_two ? 2 : 1); ++n) {
		cube_t const &cta(ct_areas[bool(n) ^ side ^ 1]); // starts with furthest to the end
		cube_t ctower;
		set_cube_zvals(ctower, roof_zval, roof_zval+ct_height);
		set_wall_width(ctower, rgen.rand_uniform((cta.d[dim][0] + ct_hlen), (cta.d[dim][1] - ct_hlen)), ct_hlen, dim); // for both cooling towers

		for (unsigned dir = 0; dir < 2; ++dir) {
			float const outer_edge(cta.d[!dim][dir]);
			ctower.d[!dim][ dir] = outer_edge;
			ctower.d[!dim][!dir] = outer_edge - (dir ? 1.0 : -1.0)*ct_width;
			details.emplace_back(ctower, ROOF_OBJ_COOLING);
			blockers[dir].assign_or_union_with_cube(ctower);
		}
	} // for n
	// add colliders between pairs of cooling towers and to the ends to avoid placing smaller AC units there so we can run pipes
	for (unsigned dir = 0; dir < 2; ++dir) {
		blockers[dir].d[dim][!side] = roof.d[dim][!side]; // extend to the edge of the roof on this side to provide space for coolant pipe bends
		details.emplace_back(blockers[dir], DETAIL_OBJ_KEEPOUT);
	}
}

bool building_t::add_server_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned lights_start) {
	assert(has_room_geom());
	bool const long_dim(room.dx() < room.dy()), mult_rows(is_datacenter());
	bool const check_windows(!(interior->dc_info && interior->dc_info->skip_sr_util_windows));
	float const window_vspacing(get_window_vspace()), ceiling_zval(zval + get_floor_ceil_gap()), wall_thickness(get_wall_thickness());
	float const server_height(0.7*window_vspacing*(mult_rows ? rgen.rand_uniform(0.9, 1.0) : rgen.rand_uniform(0.9, 1.1))); // slightly shorter if multi-row to avoid blocking lights
	float       server_width (0.3*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hwidth(0.5*server_width);
	float const server_depth (0.4*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hdepth(0.5*server_depth);
	float const comp_height  (0.2*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const min_spacing  (0.1*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const comp_hwidth(0.5*0.44*comp_height), comp_hdepth(0.5*0.9*comp_height); // fixed AR=0.44 to match the texture
	float const server_period(server_width + min_spacing), conduit_radius(0.05*server_width);
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*wall_thickness); // server spacing from walls
	cube_t inner_area(place_area);
	zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_METAL); // add server room metal tile and move the effective floor up
	cube_t server, computer;
	set_cube_zvals(server,   zval, (zval + server_height));
	set_cube_zvals(computer, zval, (zval + comp_height  ));
	point center;
	unsigned num_servers(0), num_comps(0);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const servers_start(objs.size());

	// try to line servers up against each interior wall wherever they fit
	for (unsigned D = 0; D < 2; ++D) {
		bool const dim(bool(D) ^ long_dim); // place along walls in long dim first
		float const room_len(place_area.get_sz_dim(dim));
		unsigned const num(room_len/server_period); // take the floor
		if (num == 0) continue; // not enough space for a server in this dim
		float const server_spacing(room_len/num);
		center[dim] = place_area.d[dim][0] + 0.5*server_spacing; // first position at half spacing
		bool is_ext_wall[2]={};
		
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (check_windows && classify_room_wall(room, zval, !dim, dir, 0) == ROOM_WALL_EXT) {is_ext_wall[dir] = 1;} // will skip placement on this wall
			else {inner_area.d[!dim][dir] += (dir ? -1.0 : 1.0)*server_depth;} // shrink inner area by server depth
		}
		for (unsigned n = 0; n < num; ++n, center[dim] += server_spacing) {
			set_wall_width(server, center[dim], server_hwidth, dim); // position along the wall

			for (unsigned dir = 0; dir < 2; ++dir) {
				if (is_ext_wall[dir]) continue; // not against exterior walls
				float const dir_sign(dir ? -1.0 : 1.0);
				float const wall_pos(place_area.d[!dim][dir]);
				center[!dim] = wall_pos + dir_sign*server_hdepth;
				set_wall_width(server, center[!dim], server_hdepth, !dim); // position from the wall
				cube_t server_exp(server);
				server_exp.expand_in_dim(dim, server_hwidth); // check for more side/width spacing for doors

				// Note: overlaps_other_room_obj includes previously placed servers, so we don't have to check for intersections at the corners of rooms
				if (is_obj_placement_blocked(server_exp, room, 1) || overlaps_other_room_obj(server, objs_start) || overlaps_or_adj_int_window(server_exp)) {
					// no space for server; try computer instead
					set_wall_width(computer,  center[ dim], comp_hwidth, dim); // position along the wall
					set_wall_width(computer, (wall_pos + 1.2*dir_sign*comp_hdepth), comp_hdepth, !dim); // position from the wall
					if (is_obj_placement_blocked(computer, room, 1) || overlaps_other_room_obj(computer, objs_start) || overlaps_or_adj_int_window(computer)) continue;
					objs.emplace_back(computer, TYPE_COMPUTER, room_id, !dim, !dir, 0, tot_light_amt);
					++num_comps;
					continue;
				}
				objs.emplace_back(server, TYPE_SERVER, room_id, !dim, !dir, 0, tot_light_amt);
				cube_t blocker(server);
				blocker.d[!dim][ dir]  = server.d[!dim][!dir]; // front of server
				blocker.d[!dim][!dir] += dir_sign*server_width; // add space in the front for the door to open (don't block with another server)
				objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, 0, RO_FLAG_INVIS);
				// add wire conduit on the top
				cube_t conduit;
				set_wall_width(conduit, (wall_pos + 0.75*dir_sign*server_hdepth), conduit_radius, !dim); // further toward the back wall
				set_wall_width(conduit, center[dim], conduit_radius, dim);
				set_cube_zvals(conduit, server.z2(), ceiling_zval);
				objs.emplace_back(conduit, TYPE_PIPE, room_id, 0, 1, (RO_FLAG_NOCOLL | RO_FLAG_LIT), tot_light_amt, SHAPE_CYLIN, LT_GRAY); // vertical with shadows for DC
				++num_servers;
			} // for dir
		} // for n
	} // for dim
	if (num_servers > 0) { // add a keyboard to the master server
		unsigned const master_server(rgen.rand() % num_servers);

		for (unsigned i = servers_start, server_ix = 0; i < objs.size(); ++i) {
			room_object_t const &server(objs[i]);
			if (server.type != TYPE_SERVER  ) continue;
			if (server_ix++ != master_server) continue;
			float const kbd_hwidth(0.8*server_hwidth), kbd_depth(0.6*kbd_hwidth), kbd_height(0.04*kbd_hwidth); // slightly flatter than regular keyboards
			bool const dim(server.dim), dir(server.dir);
			float const kbd_z1(server.z1() + 0.57*server.dz()), server_front(server.d[dim][dir]);
			cube_t keyboard;
			set_cube_zvals(keyboard, kbd_z1, kbd_z1+kbd_height);
			keyboard.d[dim][!dir] = server_front; // at front of server
			keyboard.d[dim][ dir] = server_front + (dir ? 1.0 : -1.0)*kbd_depth; // sticks out of the front
			set_wall_width(keyboard, server.get_center_dim(!dim), kbd_hwidth, !dim);
			if (is_obj_placement_blocked(keyboard, room, 1)) break; // Note: not checking overlaps_other_room_obj() because it will overlap server blockers
			objs.emplace_back(keyboard, TYPE_KEYBOARD, room_id, dim, dir, RO_FLAG_HANGING, tot_light_amt); // add as white, will be drawn with gray/black texture
			break;
		} // for i
	}
	if (!mult_rows) {
		// maybe add laptops on top of some servers, to reward the player for finding this room
		for (unsigned i = servers_start; i < objs.size(); ++i) {
			room_object_t const &server(objs[i]);
			if (server.type != TYPE_SERVER) continue;
			if (rgen.rand_float() > 0.2)    continue; // place laptops 20% of the time
			bool const dim(server.dim), dir(server.dir);
			float const server_front(server.d[dim][dir]); // copy before reference is invalidated
			if (!place_laptop_on_obj(rgen, server, room_id, tot_light_amt)) continue; // no avoid, use_dim_dir=0
			// make the laptop hang over the edge of the front of the server so that the player can see and take it
			room_object_t &laptop(objs.back());
			float const xlate(server_front - laptop.d[dim][dir] + (dir ? 1.0 : -1.0)*rgen.rand_uniform(0.05, 0.35)*laptop.get_sz_dim(dim));
			laptop.translate_dim(dim, xlate);
			laptop.flags |= RO_FLAG_HANGING; // make sure bottom is drawn
		} // for i
	}
	if (mult_rows) {
		// add interior rows of servers for data centers grouped into blocks
		unsigned const num_per_block = 8;
		float const inner_len(inner_area.get_sz_dim(long_dim)), inner_width(inner_area.get_sz_dim(!long_dim));
		float const clearance(get_min_front_clearance_inc_people());
		float const front_clearance(1.2*clearance), back_clearance(1.0*clearance), edge_gap(max(server_depth, 1.6f*clearance)), side_gap(0.02*server_width);
		float block_width(num_per_block*server_width + (num_per_block-1)*side_gap);
		float aisle_gap(max(server_depth, 1.3f*clearance)), block_gap(max(2.0f*server_width, 1.25f*clearance));
		float row_spacing(server_depth + aisle_gap), block_spacing(block_width + block_gap);
		float const avail_depth(inner_width - 2*edge_gap + aisle_gap), avail_width(inner_len - 2*edge_gap + block_gap);
		float const num_blocks_fp(avail_width/block_spacing); // take floor
		unsigned const num_rows(avail_depth/row_spacing); // take floor
		unsigned num_blocks(num_blocks_fp);

		if (num_blocks <= 3 && num_blocks_fp - num_blocks > 0.75) { // close to having space for an extra block
			server_width *= 0.9; // shrink servers slightly and try again
			block_width   = num_per_block*server_width + (num_per_block-1)*side_gap;
			block_spacing = block_width + block_gap;
			num_blocks    = avail_width/block_spacing;
		}
		// use a consistent direction for center rows; face the room interior door if there is one
		vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors
		bool const sdim(!long_dim);
		bool const sdir(doorways.empty() ? rgen.rand_bool() : (room.get_center_dim(sdim) < doorways.front().get_center_dim(sdim)));
		float const dsign(sdir ? 1.0 : -1.0);
		vector<vector2d> floor_vent_pos; // added at isle intersections
		assert(place_area.contains_cube(inner_area));

		if (num_rows > 0 && num_blocks > 0) {
			row_spacing   = avail_depth/num_rows;
			block_spacing = avail_width/num_blocks;
			aisle_gap     = row_spacing   - server_depth;
			block_gap     = block_spacing - block_width;
			float const block_edge_gap(0.5*(inner_len   + block_gap - avail_width));
			float const row_edge_gap  (0.5*(inner_width + aisle_gap - avail_depth));
			float const block_start(inner_area.d[!sdim][0] + block_edge_gap); // include servers along the wall and their gap
			float const server_pos_step(server_width + side_gap), row_start(inner_area.d[sdim][0] + row_edge_gap);
			float row_pos(row_start);

			for (unsigned r = 0; r < num_rows; ++r) {
				float block_pos(block_start);
				server.d[sdim][0] = row_pos;
				server.d[sdim][1] = row_pos + server_depth;

				for (unsigned b = 0; b < num_blocks; ++b) { // should culling be per-block?
					float server_pos(block_pos);

					for (unsigned n = 0; n < num_per_block; ++n) {
						server.d[!sdim][0] = server_pos;
						server.d[!sdim][1] = server_pos + server_width;
						// check for front and back clearance, but no side clearance; only check for stairs, elevators, and doors, not other objects/servers
						cube_t server_exp(server);
						server_exp.d[sdim][ sdir] += dsign*front_clearance;
						server_exp.d[sdim][!sdir] -= dsign*back_clearance;

						if (!is_obj_placement_blocked(server_exp, room, 1)) {
							objs.emplace_back(server, TYPE_SERVER, room_id, sdim, sdir, RO_FLAG_ON_FLOOR, tot_light_amt); // flag so that back is drawn
							// no conduits, since there will be too many and they may block ceiling lights
							++num_servers;
						}
						server_pos += server_pos_step;
					} // for n
					block_pos += block_spacing;
				} // for b
				row_pos += row_spacing;
			} // for r
			vector2d vent_pos;
			vent_pos[sdim] = row_start + 0.5*(server_depth - row_spacing); // centered between rows

			for (unsigned r = 0; r <= num_rows; r += 2) { // every other row
				vent_pos[!sdim] = block_start + 0.5*(num_per_block*server_pos_step - block_spacing); // centered between blocks

				for (unsigned b = 0; b <= num_blocks; ++b) {
					floor_vent_pos.push_back(vent_pos);
					vent_pos[!sdim] += block_spacing;
				}
				vent_pos[sdim] += 2.0*row_spacing;
			} // for r
		}
		// add an array of ceiling vents between the lights
		vector3d light_edge, light_space;
		cube_t ref_light;

		for (unsigned i = lights_start; i < objs_start; ++i) {
			room_object_t const &light(objs[i]);
			if (light.type != TYPE_LIGHT) continue;
			
			if (ref_light.is_all_zeros()) {
				ref_light  = light;
				light_edge = light.get_cube_center() - room.get_llc();
			}
			if (light_space.x == 0.0 && light.x1() != ref_light.x1()) {light_space.x = light.x1() - ref_light.x1();}
			if (light_space.y == 0.0 && light.y1() != ref_light.y1()) {light_space.y = light.y1() - ref_light.y1();}
		} // for i
		assert(!ref_light.is_all_zeros());
		unsigned const vent_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING);
		float const thickness(0.1*wall_thickness), hlen(2.0*wall_thickness), hwid(2.0*wall_thickness);
		cube_t vent;

		if (light_space.x > 0.0 && light_space.y > 0.0) {
			set_cube_zvals(vent, ceiling_zval-thickness, ceiling_zval);

			for (float x = place_area.x1()+0.5*light_edge.x; x < place_area.x2(); x += light_space.x) {
				for (float y = place_area.y1()+0.5*light_edge.y; y < place_area.y2(); y += light_space.y) {
					vector2d const vc(x, y);
					set_wall_width(vent, vc[ long_dim], hlen,  long_dim);
					set_wall_width(vent, vc[!long_dim], hwid, !long_dim);
					objs.emplace_back(vent, TYPE_VENT, room_id, long_dim, 0, vent_flags, 1.0); // dir=0; fully lit
				}
			} // for x
		}
		// add an array of floor vents
		set_cube_zvals(vent, zval, zval+thickness);

		for (vector2d const &vc : floor_vent_pos) {
			set_wall_width(vent, vc[ long_dim], hlen,  long_dim);
			set_wall_width(vent, vc[!long_dim], hwid, !long_dim);
			if (!is_obj_placement_blocked(vent, room, 1)) {objs.emplace_back(vent, TYPE_VENT, room_id, long_dim, 1, vent_flags, 1.0);} // dir=1; fully lit
		}
	}
	if (is_datacenter() && !check_windows && !interior->int_windows.empty()) {
		// data center with windowless walls; add interior wall surfaces so that window blending works correctly (draw order is wrong for exterior walls)
		cube_t const &part(parts[room.part_id]);
		bool const dim(!bool(hallway_dim));
		float const wall_thick(0.5*get_trim_thickness());

		for (unsigned dir = 0; dir < 2; ++dir) {
			float const wall_pos(room.d[dim][dir]);
			if (wall_pos != part.d[dim][dir]) continue; // not an exterior wall
			cube_t wall(room);
			wall.d[dim][!dir] = wall_pos - (dir ? 1.0 : -1.0)*wall_thick;
			objs.emplace_back(wall, TYPE_POOL_TILE, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_HI)); // flag as plaster wall
		} // for dir
	}
	if (num_servers == 0 && num_comps == 0) return 0; // both servers and computers count
	add_door_sign("Server Room", room, zval, room_id);
	return 1;
}

void move_to_not_intersect(cube_t const &keepout, cube_t &obj, bool dim, unsigned pref_dir=2) { // Note: caller must check for intersection
	bool const move_dir((pref_dir < 2) ? pref_dir : (keepout.get_center_dim(dim) < obj.get_center_dim(dim))); // move away from AC
	obj.translate_dim(dim, (keepout.d[dim][move_dir] - obj.d[dim][!move_dir]));
}
void move_lights_to_not_intersect(vect_room_object_t &objs, cube_t const &keepout, bool dim, unsigned objs_start, unsigned lights_start, unsigned pref_dir=2) {
	for (unsigned j = lights_start; j < objs_start; ++j) {
		room_object_t& light(objs[j]);
		if (light.type != TYPE_LIGHT) continue;
		if (light.d[dim][0] > keepout.d[dim][1] || light.d[dim][1] < keepout.d[dim][0]) continue; // not overlapping
		move_to_not_intersect(keepout, light, dim, pref_dir);
	}
}
cube_t get_exhaust_duct_for_generator(room_object_t const &generator, float duct_z2) {
	assert(generator.type == TYPE_GENERATOR);
	assert(duct_z2 > generator.z2());
	bool const gdim(generator.dim), gdir(generator.dir);
	float const gen_width(generator.get_width()), gen_length(generator.get_length());
	cube_t duct;
	set_cube_zvals(duct, generator.z2(), duct_z2);
	set_wall_width(duct, generator.get_center_dim(!gdim), 0.2*gen_width, !gdim);
	set_wall_width(duct, generator.d[gdim][gdir] - (gdir ? 1.0 : -1.0)*0.68*gen_length, 0.03*gen_length, gdim);
	return duct;
}
void add_exterior_duct(cube_t const &h_duct, bool dim, bool dir, unsigned room_id, vect_room_object_t &objs) {
	float const wall_pos(h_duct.d[dim][dir]), duct_height(h_duct.dz());
	cube_t vent(h_duct);
	vent.d[dim][!dir] = wall_pos;
	vent.d[dim][ dir] = wall_pos + (dir ? 1.0 : -1.0)*0.2*duct_height; // set depth/extension
	vent.expand_in_dim(!dim, 0.2*duct_height);
	vent.expand_in_z(0.1*duct_height);
	objs.emplace_back(vent, TYPE_VENT, room_id, dim, !dir, (RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR), 1.0); // fully lit
}

void add_pipe_through_floors(float v_pipe_radius, float zval, float dz_mult, float front_mult, float out_mult, bool dim, bool dir, bool is_ground_floor, unsigned num_floors,
	unsigned room_id, float dsign, float light_amt, float floor_spacing, unsigned pipe_type, colorRGBA const &pipe_color, colorRGBA const &fit_color,
	cube_t const &obj, vect_room_object_t &objs, vector<pipe_conn_t> &pipe_conn)
{
	unsigned const pipe_flags(RO_FLAG_NOCOLL  | RO_FLAG_LIT    ); // shadow casting
	unsigned const fitting_flags(pipe_flags   | RO_FLAG_HANGING);
	float const h_pipe_radius(0.75*v_pipe_radius), pipe_conn_z(zval + dz_mult*obj.dz());
	point pipe_pos(0.0, 0.0, zval);
	pipe_pos[ dim] = obj.get_center_dim(dim) - dsign*front_mult*obj.get_sz_dim(dim); // toward the front
	pipe_pos[!dim] = obj.d[!dim][dir]  + (dir ? 1.0 : -1.0)*out_mult*v_pipe_radius;
	cube_t v_pipe(pipe_pos);
	v_pipe.expand_by_xy(v_pipe_radius);

	if (is_ground_floor) { // pipe extends through all floors; add only for ground floor room
		v_pipe.z2() = pipe_conn_z + (num_floors-1)*floor_spacing; // extend to top floor
		objs.emplace_back(v_pipe, TYPE_PIPE, room_id, 0, 1, pipe_flags, light_amt, SHAPE_CYLIN, pipe_color); // vertical
		pipe_conn.emplace_back(pipe_pos.x, pipe_pos.y, v_pipe_radius, pipe_type);
	}
	cube_t h_pipe(pipe_pos);
	h_pipe.expand_in_dim(dim, h_pipe_radius);
	set_wall_width(h_pipe, pipe_conn_z, h_pipe_radius, 2); // Z
	h_pipe.d[!dim][!dir] = obj.get_center_dim(!dim); // connects to generator
	objs.emplace_back(h_pipe, TYPE_PIPE, room_id, !dim, 0, pipe_flags, light_amt, SHAPE_CYLIN, pipe_color); // horizontal
	// add fittings
	float const fitting_exp(0.2*v_pipe_radius);
	cube_t v_fitting(v_pipe), h_fitting(h_pipe);
	v_fitting.expand_by_xy(fitting_exp);
	set_wall_width(v_fitting, pipe_conn_z, (h_pipe_radius + fitting_exp), 2); // Z
	objs.emplace_back(v_fitting, TYPE_PIPE, room_id, 0, 1, (fitting_flags | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI), light_amt, SHAPE_CYLIN, fit_color); // vertical
	h_fitting.expand_in_z(fitting_exp);
	h_fitting.expand_in_dim(dim, fitting_exp);
	h_fitting.d[!dim][ dir] = pipe_pos[!dim]; // centerline of v_pipe
	h_fitting.d[!dim][!dir] = pipe_pos[!dim] - (dir ? 1.0 : -1.0)*2.0*v_pipe_radius;
	objs.emplace_back(h_fitting, TYPE_PIPE, room_id, !dim, 0, (fitting_flags | (dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI)), light_amt, SHAPE_CYLIN, fit_color); // horizontal
}

void building_t::add_dc_utility_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id,
	unsigned floor_ix, unsigned num_floors, float tot_light_amt, unsigned objs_start, unsigned lights_start)
{
	assert(interior->dc_info);
	zval       = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_CONCRETE); // add concrete and move the effective floor up
	objs_start = interior->room_geom->objs.size(); // exclude this from collision checks
	auto &objs(interior->room_geom->objs);
	cube_t const place_area(get_walkable_room_bounds(room));
	bool const dim(hallway_dim), dir(interior->dc_info->se_dir); // server room is opposite the stairs and elevators
	bool const bldg_side(pri_hall.get_center_dim(!dim) < room.get_center_dim(!dim)), is_ground_floor(floor_ix == 0);
	bool const has_fan_model(building_obj_model_loader.is_model_valid(OBJ_MODEL_RAD_FAN));
	float const floor_spacing(get_window_vspace()), ceiling_zval(zval + get_floor_ceil_gap()), clearance(get_min_front_clearance_inc_people()), dsign(dir ? 1.0 : -1.0);
	float cur_place_pos(place_area.d[dim][dir]); // start at front wall adjacent to server room
	// place batteries, represented as kitchen fridges, against the wall shared with the server room
	int const bat_model(select_dc_battery_model());
	cube_t batteries_bcube;

	if (bat_model >= 0) {
		float const bat_height(0.55*floor_spacing), wall_pos(cur_place_pos);
		unsigned const bat_start(objs.size());
		add_row_of_models(place_area, zval, room_id, tot_light_amt, bat_height, 0.08, bat_model, TYPE_KITCH_APP, dim, dir, dim, !dir, get_sub_model_id(bat_model), cur_place_pos);
		unsigned const bat_end(objs.size());
		// add conduits to the top; maybe these should connect together?
		bool const conduit_side(rgen.rand_bool());
		float const conduit_radius(0.02*floor_spacing), csign(conduit_side ? 1.0 : -1.0);
		cube_t conduit;
		set_cube_zvals(conduit, zval+bat_height, ceiling_zval);

		for (unsigned i = bat_start; i != bat_end; ++i) {
			room_object_t const &bat(objs[i]);
			assert(bat.type == TYPE_KITCH_APP);
			set_wall_width(conduit, (wall_pos                  - 0.2*dsign*bat.get_depth()), conduit_radius,  dim); // further toward the back wall
			set_wall_width(conduit, (bat.d[!dim][conduit_side] - 0.1*csign*bat.get_width()), conduit_radius, !dim); // off to one side
			objs.emplace_back(conduit, TYPE_PIPE, room_id, 0, 1, (RO_FLAG_NOCOLL | RO_FLAG_LIT), tot_light_amt, SHAPE_CYLIN, GRAY_BLACK); // vertical, with shadows
			batteries_bcube.assign_or_union_with_cube(bat);
		}
		if (has_fan_model) { // add fans to the top
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAD_FAN)); // D, W, H

			for (unsigned i = bat_start; i != bat_end; ++i) {
				room_object_t const &bat(objs[i]);
				float const scale(0.67*min(bat.get_depth()/sz.z, bat.get_width()/sz.y)), height(scale*sz.x), width(scale*sz.y), depth(scale*sz.z);
				float const fan_z1(bat.z2() - 1.2*height); // shift slightly down into the top of the battery, more than we normally would due to something wrong in the model translate
				cube_t fan;
				set_cube_zvals(fan, fan_z1, fan_z1+height);
				set_wall_width(fan, bat.get_center_dim( dim), 0.5*depth,  dim);
				set_wall_width(fan, bat.get_center_dim(!dim), 0.5*width, !dim);
				objs.emplace_back(fan, TYPE_RAD_FAN, room_id, dim, dir, (RO_FLAG_ADJ_TOP | RO_FLAG_NOCOLL), tot_light_amt); // vertical
			} // for i
		}
	}
	cube_t inner_area(place_area);
	inner_area.expand_in_dim(!dim, -1.5*clearance);
	
	// place AC units
	float &ac_height(interior->dc_info->ac_height), &ac_width(interior->dc_info->ac_width), &ac_depth(interior->dc_info->ac_depth); // shared across all building utility rooms
	bool &cw_pipe_side(interior->dc_info->cw_pipe_side), &drain_pipe_side(interior->dc_info->drain_pipe_side), &ac_pipe_end(interior->dc_info->ac_pipe_end);
	bool &gen_pipe_side(interior->dc_info->gen_pipe_side), &gen_dir(interior->dc_info->gen_dir), &xg_side(interior->dc_info->xg_side);

	if (ac_height == 0.0) { // assign these values on the first floor of the first utility room
		ac_height       = rgen.rand_uniform(0.45, 0.54)*floor_spacing;
		ac_width        = rgen.rand_uniform(0.20, 0.24)*floor_spacing;
		ac_depth        = rgen.rand_uniform(0.40, 0.60)*floor_spacing;
		cw_pipe_side    = rgen.rand_bool();
		drain_pipe_side = rgen.rand_bool();
		ac_pipe_end     = rgen.rand_bool();
		gen_pipe_side   = rgen.rand_bool();
		gen_dir         = rgen.rand_bool(); // either facing toward or away from the door
		xg_side         = rgen.rand_bool();
	}
	colorRGBA const ac_color(0.5, 0.55, 0.6, 1.0); // blue-gray
	cube_t ac_area(inner_area);
	ac_area.expand_in_dim(!dim, -0.25*ac_width); // add extra padding for fans
	unsigned const ac_start(objs.size()), item_flags(0), obj_flags(RO_FLAG_IN_FACTORY); // flag as factor so that metal is painted and less reflective
	add_row_of_objects(ac_area, zval, room_id, tot_light_amt, ac_height, ac_width, ac_depth, 1.33, TYPE_METAL_BAR, dim, dir, dim, dir, item_flags, obj_flags, ac_color, cur_place_pos);
	unsigned const ac_end(objs.size());

	if (ac_start < ac_end && has_fan_model) { // add AC side fans
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAD_FAN)); // D, W, H
		float const scale(0.85*min(ac_depth/sz.y, min(0.5f*ac_width/sz.x, 0.5f*ac_height/sz.z)));
		float const height(scale*sz.z), width(scale*sz.y), depth(scale*sz.x);

		for (unsigned i = ac_start; i < ac_end; ++i) {
			cube_t const ac(objs[i]);
			cube_t fan;
			fan.z2() = ac .z2() - 0.1*height;
			fan.z1() = fan.z2() -     height;
			set_wall_width(fan, ac.get_center_dim(dim), 0.5*width, dim);

			for (unsigned d = 0; d < 2; ++d) { // each side
				float const edge_pos(ac.d[!dim][d]);
				fan.d[!dim][!d] = edge_pos;
				fan.d[!dim][ d] = edge_pos + (d ? 1.0 : -1.0)*depth;
				objs.emplace_back(fan, TYPE_RAD_FAN, room_id, !dim, d, 0, tot_light_amt, SHAPE_CUBE, WHITE);
			}
		} // for i
	}
	// add AC ducts and pipes
	float const h_duct_hwidth(0.5*ac_width), h_duct_height(0.5*ac_width), v_duct_radius(0.3*ac_width), h_duct_z1(ceiling_zval - h_duct_height);
	float const bs_sign(bldg_side ? 1.0 : -1.0), v_pipe_z2(0.5*(h_duct_z1 + ceiling_zval)); // Z midpoint of h duct
	float const ac_pipe_radius(0.07*ac_width), ac_hp_radius(1.33*ac_pipe_radius), ac_fit_thick(0.1*ac_hp_radius);
	float const far_wall_pos(room.d[!dim][bldg_side]), far_wall_inner(far_wall_pos - bs_sign*h_duct_height);
	unsigned const duct_flags(RO_FLAG_ADJ_TOP | RO_FLAG_ADJ_BOT); // skip top and bottom
	unsigned const pipe_flags(RO_FLAG_NOCOLL  | RO_FLAG_LIT    ); // shadow casting
	unsigned const fitting_flags(pipe_flags   | RO_FLAG_HANGING);
	unsigned const skip_end_flag(bldg_side ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO); // end against the far wall
	colorRGBA const ac_pipe_color(0.2, 0.5, 1.0), ac_fit_color(BRASS_C);
	cube_t h_pipe;

	for (unsigned i = ac_start; i < ac_end; ++i) { // add ducts connecting AC units
		bool const closest_to_door(i == (bldg_side ? ac_start : ac_end-1));
		cube_t const ac(objs[i]);
		vector2d const ac_center(ac.xc(), ac.yc());
		// add vertical duct connecting to the top of each AC unit
		cube_t duct;
		set_cube_zvals(duct, ac.z2(), h_duct_z1);
		for (unsigned d = 0; d < 2; ++d) {set_wall_width(duct, ac_center[d], v_duct_radius, d);}
		objs.emplace_back(duct, TYPE_DUCT, room_id, 0, 1, duct_flags, tot_light_amt, SHAPE_CYLIN); // vertical
		
		if (closest_to_door) {
			// add main duct for first AC unit
			duct.d[!dim][!bldg_side] -= bs_sign*0.5*v_duct_radius; // extend out slightly
			duct.d[!dim][ bldg_side]  = far_wall_pos; // overlaps vertical duct so that it connects even if misaligned
			set_cube_zvals(duct, h_duct_z1, ceiling_zval);
			set_wall_width(duct, ac_center[dim], h_duct_hwidth, dim);
			objs.emplace_back(duct, TYPE_DUCT, room_id, !dim, 0, (RO_FLAG_ADJ_TOP | skip_end_flag), tot_light_amt, SHAPE_CUBE); // horizontal; skip top and back
			// first in row; check if we need to move an entire row of lights
			cube_t keepout(ac); // use most of entire AC as keepout so that we have space to route pipes into the ceiling
			keepout.expand_in_dim(dim, -0.05*ac_depth);
			move_lights_to_not_intersect(objs, keepout, dim, objs_start, lights_start);
			// add vertical duct connecting to the air intake; this should at least partially overlap and connect to the horizontal duct
			duct.d[!dim][!bldg_side] = far_wall_inner;
			duct.d[!dim][ bldg_side] = far_wall_pos;
			duct.expand_in_dim(dim, 0.5*h_duct_hwidth); // double the width

			if (is_ground_floor) { // ground floor duct only extends upward unless there's a parking garage below to extend into
				if (has_parking_garage) {
					duct.z1() = zval; // ground_floor_z1?
					interior->dc_info->intake_ducts[bldg_side] = duct;
				}
				else {duct.z1() -= 0.25*h_duct_height;}
			}
			else {duct.z1() = zval;}
			unsigned const v_duct_flags((is_ground_floor ? RO_FLAG_ADJ_TOP : duct_flags) | skip_end_flag); // draw bottom for ground floor since it doesn't extend down to the floor
			objs.emplace_back(duct, TYPE_DUCT, room_id, !dim, 1, v_duct_flags, tot_light_amt, SHAPE_CUBE); // vertical; skip top, bottom, and back
			duct.z1() = h_duct_z1; // restore original z1 for the vent
			add_exterior_duct(duct, !dim, bldg_side, room_id, objs); // add exterior vent

			if (zval + floor_spacing >= room.z2()) { // add vertical roof duct on top floor
				float const outer_edge(far_wall_pos - bs_sign*get_roof_wall_thick()); // shift away from the roof wall
				cube_t top_vent(duct);
				top_vent.d[!dim][ bldg_side] = outer_edge;
				top_vent.d[!dim][!bldg_side] = outer_edge - 1.5*bs_sign*h_duct_height; // set depth
				top_vent.expand_in_dim(dim, 0.15*h_duct_height); // widen
				set_cube_zvals(top_vent, room.z2(), (room.z2() + 1.0*h_duct_height)); // set height
				objs.emplace_back(top_vent, TYPE_VENT, room_id, dim, 1, (RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_EXTERIOR), 1.0); // dir=1; fully lit
			}
		}
		// add vertical pipes on each side
		for (unsigned side = 0; side < 2; ++side) {
			// add intake and outflow coolant pipes
			bool const ac_pipe_side(bldg_side ^ bool(side)); // AC pipe dir; opposite corners
			point pipe_center(0.0, 0.0, ac.z2());
			pipe_center[ dim] = ac_center[dim] + (side ? 1.0 : -1.0)*(h_duct_hwidth + ac_hp_radius + ac_fit_thick); // side of the duct including pipe and fitting radius
			pipe_center[!dim] = ac.d[!dim][ac_pipe_side] - (ac_pipe_side ? 1.0 : -1.0)*2.0*ac_pipe_radius; // near the edge
			cube_t v_pipe(pipe_center);
			v_pipe.expand_by_xy(ac_pipe_radius);
			v_pipe.z2() = v_pipe_z2;
			objs.emplace_back(v_pipe, TYPE_PIPE, room_id, 0, 1, pipe_flags, tot_light_amt, SHAPE_CYLIN, ac_pipe_color); // vertical
			// add brass fittings on bottom (AC) and top (h_pipe)
			cube_t fitting(v_pipe);
			fitting.expand_by_xy(ac_fit_thick);
			fitting.z2() = v_pipe.z1() + 1.0*ac_pipe_radius;
			objs.emplace_back(fitting, TYPE_PIPE, room_id, 0, 1, (fitting_flags | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, ac_fit_color); // vertical, draw top end
			set_cube_zvals(fitting, (v_pipe.z2() - 2.0*ac_pipe_radius), v_pipe.z2());
			objs.emplace_back(fitting, TYPE_PIPE, room_id, 0, 1, (fitting_flags | RO_FLAG_ADJ_LO), tot_light_amt, SHAPE_CYLIN, ac_fit_color); // vertical, draw bot end

			if (closest_to_door) {
				// add horizontal pipe(s), then bend and go up into the ceiling
				float const pos_offset((side ? 1.0 : -1.0)*2.1*(floor_ix + 1)*ac_hp_radius); // adjacent with a small gap, and an extra gap at the edge for the vent
				float const bend_pos1(far_wall_inner - bs_sign*ac_hp_radius); // !dim
				float const bend_pos2(duct.d[dim][side] + pos_offset); // dim; unique per-floor
				float const bend_pos3(far_wall_pos - bs_sign*ac_hp_radius); // !dim
				set_wall_width(h_pipe, v_pipe_z2, ac_hp_radius, 2); // set zvals
				set_wall_width(h_pipe, pipe_center[dim], ac_hp_radius, dim);
				h_pipe.d[!dim][!bldg_side] = v_pipe.d[!dim][!bldg_side];
				h_pipe.d[!dim][ bldg_side] = bend_pos1;
				objs.emplace_back(h_pipe, TYPE_PIPE, room_id, !dim, 0, (pipe_flags | skip_end_flag), tot_light_amt, SHAPE_CYLIN, ac_pipe_color); // horizontal, bend at far wall
				cube_t h_pipe2(h_pipe); // copy zvals
				set_wall_width(h_pipe2, bend_pos1, ac_hp_radius, !dim);
				h_pipe2.d[dim][!side] = pipe_center[dim];
				h_pipe2.d[dim][ side] = bend_pos2;
				objs.emplace_back(h_pipe2, TYPE_PIPE, room_id, dim, 0, pipe_flags, tot_light_amt, SHAPE_CYLIN, ac_pipe_color); // horizontal
				cube_t h_pipe3(h_pipe2); // copy zvals
				set_wall_width(h_pipe3, bend_pos2, ac_hp_radius, dim);
				h_pipe3.d[!dim][!bldg_side] = bend_pos1;
				h_pipe3.d[!dim][ bldg_side] = bend_pos3;
				objs.emplace_back(h_pipe3, TYPE_PIPE, room_id, !dim, 0, (pipe_flags | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, ac_pipe_color); // horizontal, both bends
				cube_t v_pipe2(h_pipe3); // copy dim
				set_cube_zvals(v_pipe2, v_pipe.z2(), room.z2()); // all the way to the top of the building
				set_wall_width(v_pipe2, bend_pos3, ac_hp_radius, !dim);
				objs.emplace_back(v_pipe2, TYPE_PIPE, room_id, 0, 1, pipe_flags, tot_light_amt, SHAPE_CYLIN, ac_pipe_color); // vertical
				// route coolant pipes to the nearest cooling tower(s) on the roof
				vect_cube_t towers;
				cube_t towers_bcube;

				for (roof_obj_t const &d : details) {
					if (d.type != ROOF_OBJ_COOLING) continue;
					if ((bcube.get_center_dim(!dim) < d.get_center_dim(!dim)) != bldg_side) continue; // wrong side of building
					towers.push_back(d);
					towers_bcube.assign_or_union_with_cube(d);
				}
				unsigned const rp_flags(RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR);
				float const min_edge_space(ac_hp_radius), rise_pos(bend_pos3 - bs_sign*1.25*get_roof_wall_thick()); // along the roof wall
				float const conn_pt(towers_bcube.d[!dim][bldg_side] - bs_sign*((dir ? -1.0 : 1.0)*pos_offset + ((side ^ dir) ? 0.4 : 0.25)*towers_bcube.get_sz_dim(!dim)));
				cube_t rv_pipe(v_pipe2), rh_pipe(v_pipe2);
				set_cube_zvals(rv_pipe, room.z2(), room.z2()+    ac_hp_radius);
				set_cube_zvals(rh_pipe, room.z2(), room.z2()+2.0*ac_hp_radius); // on top of the roof
				set_wall_width(rv_pipe, rise_pos, ac_hp_radius, !dim);
				objs.emplace_back(rv_pipe, TYPE_PIPE, room_id, 0, 1, (rp_flags | RO_FLAG_ADJ_HI), 1.0, SHAPE_CYLIN, ac_pipe_color); // vertical, draw top end
				cube_t rc_pipe(rh_pipe); // connector pipe; copy zvals
				bool was_connected(0);
				rh_pipe.d[!dim][bldg_side] = rise_pos;

				for (cube_t const &t : towers) {
					if (bend_pos2 > t.d[dim][0]+min_edge_space && bend_pos2 < t.d[dim][1]-min_edge_space) { // can connect straight across
						rh_pipe.d[!dim][!bldg_side] = t.d[!dim][bldg_side]; // connects to the side of the cooling tower
						objs.emplace_back(rh_pipe, TYPE_PIPE, room_id, !dim, 0, rp_flags, 1.0, SHAPE_CYLIN, ac_pipe_color); // horizontal
						was_connected = 1;
						break;
					}
				} // for t
				if (towers.size() > 1 && floor_ix <= num_floors/2) { // add horizontal pipe connecting towers in this row; only for half the floors, rounded up
					cube_t conn_area(towers_bcube);
					conn_area.expand_in_dim(dim, -towers.front().get_sz_dim(dim)); // space between towers; assume all lengths are the same
					assert(conn_area.is_strictly_normalized());
					copy_dim(rc_pipe, conn_area, dim);
					set_wall_width(rc_pipe, conn_pt, ac_hp_radius, !dim);
					objs.emplace_back(rc_pipe, TYPE_PIPE, room_id, dim, 0, rp_flags, 1.0, SHAPE_CYLIN, ac_pipe_color); // horizontal

					if (!was_connected && bend_pos2 > towers_bcube.d[dim][0] && bend_pos2 < towers_bcube.d[dim][1]) { // add a T-junction
						// connect to rc_pipe; must route up and over; not implemented because this case may not be reachable
						was_connected = 1;
					}
				}
				if (!was_connected) { // need to add a bend
					bool const bdir(bend_pos2 < towers_bcube.get_center_dim(dim)); // bend toward the group of towers
					rh_pipe.d[!dim][!bldg_side] = conn_pt; // bend point
					objs.emplace_back(rh_pipe, TYPE_PIPE, room_id, !dim, 0, rp_flags, 1.0, SHAPE_CYLIN, ac_pipe_color); // horizontal
					cube_t rh_pipe2(rh_pipe); // copy zvals
					set_wall_width(rh_pipe2, conn_pt, ac_hp_radius, !dim);
					rh_pipe2.d[dim][!dir] = bend_pos2;
					rh_pipe2.d[dim][ dir] = towers_bcube.d[dim][!dir];
					if (rh_pipe2.get_sz_dim(dim) <= 0.0) {rh_pipe2.d[dim][dir] += (dir ? 1.0 : -1.0)*min_edge_space;} // extend into cooling tower if too short
					objs.emplace_back(rh_pipe2, TYPE_PIPE, room_id, dim, 0, (rp_flags | (dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI)), 1.0, SHAPE_CYLIN, ac_pipe_color); // horizontal
				}
			}
			// add fittings to h_pipe
			fitting = h_pipe;
			fitting.expand_in_z(ac_fit_thick);
			fitting.expand_in_dim(dim, ac_fit_thick);
			set_wall_width(fitting, pipe_center[!dim], 1.0*ac_hp_radius, !dim);
			objs.emplace_back(fitting, TYPE_PIPE, room_id, !dim, 0, (fitting_flags | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, ac_fit_color); // horizontal, draw ends
			// add water pipe and drain pipe
			add_pipe_through_floors(0.010*floor_spacing, zval, 0.9, (ac_pipe_end ? 1.0 : -1.0)*0.45, 4.0, dim, cw_pipe_side,    is_ground_floor, num_floors, room_id,
				dsign, tot_light_amt, floor_spacing, PIPE_TYPE_CW,    COPPER_C, BRASS_C, ac, objs, interior->dc_info->pipe_conn);
			add_pipe_through_floors(0.012*floor_spacing, zval, 0.1, (ac_pipe_end ? -1.0 : 1.0)*0.45, 4.0, dim, drain_pipe_side, is_ground_floor, num_floors, room_id,
				dsign, tot_light_amt, floor_spacing, PIPE_TYPE_SEWER, WHITE,    WHITE,   ac, objs, interior->dc_info->pipe_conn);
		} // for side
	} // for i
	cube_t xfmr_area(inner_area), gen_area(inner_area);
	unsigned gen_start(0);

	for (unsigned pass = 0; pass < 2; ++pass) {
		float const pre_place_pos(cur_place_pos);
		unsigned const xg_obj_size(objs.size());
		// place transformers
		float const xfmr_height(0.4*floor_spacing), xfmr_zval(zval - TRANSFORMER_Z_SHIFT*xfmr_height); // shift base to the floor
		add_row_of_models(xfmr_area, xfmr_zval, room_id, tot_light_amt, xfmr_height, 0.15, OBJ_MODEL_SUBSTATION, TYPE_XFORMER, dim, dir, dim, dir, 0, cur_place_pos);
		// connect with wire conduits? but they already have conduits into the floor
		if (pass == 1) {cur_place_pos = pre_place_pos;} // use the same row
		// place generators near the back wall
		gen_start = objs.size();
		float const gen_height(0.56*floor_spacing), gap((pass == 0) ? 1.0 : 0.25); // more gap/sparser on first pass
		if (add_row_of_models(gen_area, zval, room_id, tot_light_amt, gen_height, gap, OBJ_MODEL_GENERATOR, TYPE_GENERATOR, dim, dir,  dim,  dir, 0, cur_place_pos)) break;
		// can't place lengthwise; try sideways
		if (add_row_of_models(gen_area, zval, room_id, tot_light_amt, gen_height, 0.1, OBJ_MODEL_GENERATOR, TYPE_GENERATOR, dim, dir, !dim, (gen_dir ^ bldg_side), 0, cur_place_pos)) break;
		if (pass == 1) break; // failed, done
		// try again, but this time split the width in half and try to place each type
		objs.resize(xg_obj_size); // remove any models added above
		xfmr_area.d[!dim][xg_side ^ bldg_side] = gen_area.d[!dim][!xg_side ^ bldg_side] = inner_area.get_center_dim(!dim);
		cur_place_pos = pre_place_pos;
	} // for pass
	if (gen_start > 0) { // add ducts and fuel lines for generators
		unsigned const gen_end(objs.size());

		for (unsigned i = gen_start; i < gen_end; ++i) {
			room_object_t const generator(objs[i]); // deep copy to avoid invalidating the reference
			cube_t const v_duct(get_exhaust_duct_for_generator(generator, ceiling_zval));
			// add horizontal duct out to the back of the building
			vector2d const duct_sz(v_duct.get_size_xy());
			float const duct_height(duct_sz.get_min_val()), h_duct_z1(ceiling_zval - duct_height), back_wall_pos(room.d[dim][!dir]);
			unsigned h_duct_flags(RO_FLAG_ADJ_TOP | RO_FLAG_ADJ_HI | RO_FLAG_ADJ_LO); // skip top and both ends
			cube_t h_duct(v_duct);
			set_cube_zvals(h_duct, h_duct_z1, ceiling_zval);
			h_duct.d[dim][ dir] = v_duct.d[dim][!dir];
			h_duct.d[dim][!dir] = back_wall_pos; // extend to the back wall
			
			if (duct_sz[!dim] < duct_sz[dim]) { // sideways generator placement
				h_duct.expand_in_dim(!dim, 0.5*(duct_sz[dim] - duct_sz[!dim])); // expand to keep fixed cross section area
				h_duct.d[dim][dir] = v_duct.d[dim][dir] + dsign*0.1*duct_sz[dim]; // overlaps and covers v_duct with a bit of extension
				h_duct_flags &= ~(dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO); // back face is visible
			}
			objs.emplace_back(h_duct, TYPE_DUCT, room_id, dim, dir, h_duct_flags, tot_light_amt, SHAPE_CUBE); // horizontal
			objs.emplace_back(v_duct, TYPE_DUCT, room_id, 0,   1,     duct_flags, tot_light_amt, SHAPE_CUBE); // vertical
			add_exterior_duct(h_duct, dim, !dir, room_id, objs); // add exterior vent

			if (i == gen_start) { // first in row; check if we need to move an entire row of lights
				cube_t keepout(v_duct);
				keepout.union_with_cube(h_duct); // avoid this as well
				keepout.expand_in_dim(dim, 0.1*h_duct_hwidth);
				move_lights_to_not_intersect(objs, keepout, dim, objs_start, lights_start, dir); // pref_dir=dir so that lights move toward the interior and aren't blocked by h_duct
			}
			// add fuel pipe
			add_pipe_through_floors(0.008 * floor_spacing, zval, 0.36, 0.28, 2.0, generator.dim, gen_pipe_side, is_ground_floor, num_floors, room_id,
				(generator.dir ? 1.0 : -1.0), tot_light_amt, floor_spacing, PIPE_TYPE_GAS, YELLOW, LT_GRAY, generator, objs, interior->dc_info->pipe_conn);
		} // for i
	}
	// add breaker panel
	for (unsigned n = 0; n < 4; ++n) { // 4 attempts to avoid battery
		if (!add_breaker_panel_by_door(rgen, room, zval, room_id, tot_light_amt)) break;
		room_object_t breaker(objs.back());
		assert(breaker.type == TYPE_BRK_PANEL);
		float const bp_width(breaker.get_width());
		breaker.expand_in_dim( breaker.dim, 1.1*bp_width); // add clearance for the door to open
		breaker.expand_in_dim(!breaker.dim, 0.1*bp_width);
		if (cube_int_if_nonzero(batteries_bcube, breaker)) {objs.pop_back();} else {break;}
	}
	add_door_sign("Utility", room, zval, room_id);
}

bool building_t::add_row_of_objects(cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, float height, float width, float depth,
	float gap_mult, unsigned type, bool dim, bool dir, bool obj_dim, bool obj_dir, unsigned item_flags, unsigned obj_flags, colorRGBA const &color, float &cur_pos)
{
	float const clearance(get_min_front_clearance_inc_people()), min_spacing((1.0 + gap_mult)*width), row_space(1.2*clearance), avail_width(place_area.get_sz_dim(!dim));
	unsigned const num(avail_width/min_spacing);
	if (num == 0) return 0;
	float const spacing(avail_width/num), dsign(dir ? 1.0 : -1.0);
	float tpos(place_area.d[!dim][0] + 0.5*spacing);
	cube_t obj;
	set_cube_zvals(obj, zval, zval+height);
	obj.d[dim][ dir] = cur_pos;
	float new_pos(cur_pos - dsign*depth); // move away from the wall
	obj.d[dim][!dir] = new_pos;
	if (new_pos < place_area.d[dim][0] || new_pos > place_area.d[dim][1]) return 0; // can't fit in the space
	cur_pos = new_pos - dsign*row_space; // add a gap between the next row; max extend outside the room, but that's okay
	bool placed(0);

	for (unsigned n = 0; n < num; ++n, tpos += spacing) {
		set_wall_width(obj, tpos, 0.5*width, !dim);
		if (is_obj_placement_blocked(obj, place_area, 1, 1)) continue; // inc_open_doors=1, check_open_dir=1
		interior->room_geom->objs.emplace_back(obj, type, room_id, obj_dim, obj_dir, obj_flags, tot_light_amt, SHAPE_CUBE, color, item_flags);
		placed = 1;
	}
	return placed;
}
bool building_t::add_row_of_models(cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, float height, float gap_mult,
	unsigned model_id, unsigned type, bool dim, bool dir, bool obj_dim, bool obj_dir, unsigned item_flags, float &cur_pos)
{
	if (!building_obj_model_loader.is_model_valid(model_id)) return 0;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
	float width(height*sz.y/sz.z), depth(height*sz.x/sz.z);
	if (obj_dim != dim) {swap(width, depth);}
	return add_row_of_objects(place_area, zval, room_id, tot_light_amt, height, width, depth, gap_mult, type, dim, dir, obj_dim, obj_dir, item_flags, 0, WHITE, cur_pos);
}

bool building_t::add_breaker_panel_by_door(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors
	if (doorways.empty()) return 0; // should always have a door, unless we happen to have a room in a small part at the top of a skyscraper?
	door_stack_t const &door(doorways.front()); // choose the first door (there is likely only one)
	bool const side(!door.get_check_dirs()), dim(door.dim), dir(door.get_center_dim(dim) > room.get_center_dim(dim)); // the wall the door is on
	cube_t const room_bounds(get_room_wall_bounds(room));
	float const door_edge(door.d[!dim][side]), wall_edge(room_bounds.d[!dim][side]), ceil_zval(zval + get_floor_ceil_gap());
	float const wall_len(fabs(door_edge - wall_edge)), wall_center(0.5*(door_edge + wall_edge)), wall_pos(room_bounds.d[dim][dir]), floor_spacing(get_window_vspace());
	if (wall_len < 0.5*floor_spacing) return 0; // too narrow
	float const width(min(0.5f*wall_len, rgen.rand_uniform(0.25, 0.35)*floor_spacing)), depth(0.04*floor_spacing);
	cube_t breaker_panel;
	set_cube_zvals(breaker_panel, (ceil_zval - 0.75*floor_spacing), (ceil_zval - rgen.rand_uniform(0.25, 0.3)*floor_spacing));
	set_wall_width(breaker_panel, wall_center, 0.5*width, !dim);
	breaker_panel.d[dim][ dir] = wall_pos;
	breaker_panel.d[dim][!dir] = wall_pos + (dir ? -1.0 : 1.0)*depth;
	add_breaker_panel(rgen, breaker_panel, ceil_zval, door.dim, dir, room_id, tot_light_amt);
	breaker_panel.d[dim][!dir] += (dir ? -1.0 : 1.0)*width; // add padding for desk placement
	return 1;
}

void clip_to_exclude(cube_t &c, cube_t const &exclude, bool dim) {
	bool const dir(c.get_center_dim(dim) < exclude.get_center_dim(dim));
	c.d[dim][dir] = exclude.d[dim][!dir];
}
void building_t::add_op_center_objs(rand_gen_t rgen, room_t const &room, colorRGBA const &chair_color, float zval,
	unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start)
{
	vect_cube_t blockers;
	auto &objs(interior->room_geom->objs);
	bool const dim(hallway_dim), dir(!interior->dc_info->se_dir); // operation center room is in the direction of the stairs and elevators
	float const floor_spacing(get_window_vspace()), clearance(get_min_front_clearance()), width(0.8*floor_spacing*rgen.rand_uniform(1.0, 1.2));
	float const depth(0.4*floor_spacing*rgen.rand_uniform(1.0, 1.2)), height(0.25*floor_spacing*rgen.rand_uniform(1.08, 1.2)); // slightly larger than normal
	cube_t const place_area(get_walkable_room_bounds(room));
	float cur_place_pos(place_area.d[dim][dir]); // start at front wall adjacent to server room
	unsigned const desks_start(objs.size());
	add_row_of_objects(place_area, zval, room_id, tot_light_amt, height, width, depth, 0.05, TYPE_DESK, dim, dir, dim, !dir, 0, 0, WHITE, cur_place_pos);

	if (room.get_sz_dim(dim) > 2.0*depth + 0.5*floor_spacing + 3.0*clearance) { // if wide enough, add a second row of desks along the windows
		cur_place_pos = place_area.d[dim][!dir]; // opposite wall
		add_row_of_objects(place_area, zval, room_id, tot_light_amt, height, width, depth, 0.05, TYPE_DESK, dim, !dir, dim, dir, 0, 0, WHITE, cur_place_pos);
	}
	unsigned const desks_end(objs.size());

	for (unsigned i = desks_start; i < desks_end; ++i) {
		objs[i].obj_id = rgen.rand(); // randomize drawer contents
		objs[i].flags |= RO_FLAG_UNTEXTURED; // mark as plastic
		add_desk_objects(rgen, i, chair_color, room, blockers, 1, 1, 1.0); // add_computer=1, add_phone=1, comp_sz_scale=1.0
	}
	// add up to to 8 random desks along other walls
	for (unsigned n = 0; n < 8; ++n) {add_desk_to_room(rgen, room, blockers, chair_color, zval, room_id, tot_light_amt, objs_start, 0);} // is_basement=0
	// maybe add a conference table in the center of the room
	cube_t room_clipped(room);

	for (stairwell_t const &s : interior->stairwells) { // avoid stairs
		if (s.intersects(room_clipped)) {clip_to_exclude(room_clipped, s, !dim);}
	}
	for (elevator_t const &e : interior->elevators) { // avoid elevator
		if (e.intersects(room_clipped)) {clip_to_exclude(room_clipped, e, !dim);}
	}
	cube_t ct_area(place_area), room_inner(place_area);
	room_inner.expand_in_dim(!dim, -depth); // add space around edges of room for desks
	room_clipped.intersect_with_cube(room_inner);
	room_clipped.expand_in_dim(!dim, -clearance);
	ct_area.intersect_with_cube(room_clipped);
	bool added_conf_table(0);
	if (ct_area.is_strictly_normalized()) {added_conf_table = add_conference_table(rgen, room_clipped, zval, room_id, tot_light_amt, !dim, ct_area);}

	if (!added_conf_table) { // if no conference table, add table with up to 6 chairs; should fit in the center of the room
		add_table_and_chairs(rgen, room, blockers, room_id, get_cube_center_zval(room_clipped, zval), chair_color, 0.1, tot_light_amt, 6);
	}
	// add filing cabinets
	unsigned const num_filing_cabinets(2 + (rgen.rand()%4)); // 2-5

	for (unsigned n = 0; n < num_filing_cabinets; ++n) {
		if (add_filing_cabinet_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start)) {blockers.push_back(objs.back());}
	}
	add_clock_to_room_wall(rgen, room, zval, room_id, tot_light_amt, objs_start, 1); // digital
}

