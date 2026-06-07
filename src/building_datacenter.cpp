// 3D World - Data Center Buildings
// by Frank Gennari 5/25/26

#include "function_registry.h"
#include "buildings.h"

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
	bool const se_end(rgen.rand_bool()); // also office end, opposite utility end, near front door
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

	if (num_short_wind >= 9 && num_long_wind >= 9) { // larger building, at least 4 windows per side; add secondary hallway and two rows of side rooms
		// TODO
	}
	// add hallway room
	unsigned const hall_room_id(add_room(hall, part_id, 3, 1, 0)); // add primary hallway as room with 3+ lights
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
			add_assigned_room(bathroom, part_id, RTYPE_BATH);
			br_door_pos = bathroom.get_center_dim(max_dim); // centered
		}
		add_assigned_room(office, part_id, RTYPE_OFFICE);
		if (br_side) {rooms.back().assign_to(RTYPE_SECURITY, 0, 1);} // small office next to bathroom is security room on the first floor
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

		for (unsigned e = 0; e < (br_side ? 3 : 2); ++e) { // each room to split with a wall
			cube_t rwall(part);
			rwall.d[min_dim][!d] = hall_side + (d ? 1.0 : -1.0)*wall_hthick;
			create_wall(rwall, max_dim, room_split_wall_pos[e], fc_thick, wall_hthick, wall_edge_spacing);

			if ((conn_office_to_server && e == 0 && !br_side) || (conn_util_to_server && e == 1)) { // add a door between this room and the server room
				bool const open_dir(bool(e) != se_end);
				float const door_pos(rgen.rand_uniform(rwall.d[min_dim][0]+doorway_hwidth, rwall.d[min_dim][1]-doorway_hwidth));
				cube_t rwall2;
				remove_section_from_cube_and_add_door(rwall, rwall2, (door_pos - doorway_hwidth), (door_pos + doorway_hwidth), min_dim, open_dir);
				room_walls.push_back(rwall2);
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

bool building_t::add_server_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) { // for office buildings
	bool const long_dim(room.dx() < room.dy()), mult_rows(is_datacenter());
	float const window_vspacing(get_window_vspace()), ceiling_zval(zval + get_floor_ceil_gap());
	float const server_height(0.7*window_vspacing*rgen.rand_uniform(0.9, 1.1)*(mult_rows ? 0.9 : 1.0)); // slightly shorter if mulri-row to avoid blocking lights
	float const server_width (0.3*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hwidth(0.5*server_width);
	float const server_depth (0.4*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hdepth(0.5*server_depth);
	float const comp_height  (0.2*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const min_spacing  (0.1*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const comp_hwidth(0.5*0.44*comp_height), comp_hdepth(0.5*0.9*comp_height); // fixed AR=0.44 to match the texture
	float const server_period(server_width + min_spacing), conduit_radius(0.05*server_width);
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // server spacing from walls
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
			if (classify_room_wall(room, zval, !dim, dir, 0) == ROOM_WALL_EXT) {is_ext_wall[dir] = 1;} // will skip placement on this wall
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
				objs.emplace_back(conduit, TYPE_PIPE, room_id, 0, 1, RO_FLAG_NOCOLL, 1.0, SHAPE_CYLIN, LT_GRAY); // vertical
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
		float const front_clearance(1.2*clearance), back_clearance(1.0*clearance), edge_gap(max(server_depth, 1.6f*clearance));
		float const side_gap(0.02*server_width), block_width(num_per_block*server_width + (num_per_block-1)*side_gap);
		float aisle_gap(max(server_depth, 1.5f*clearance)), block_gap(max(2.0f*server_width, 1.25f*clearance));
		float row_spacing(server_depth + aisle_gap), block_spacing(block_width + block_gap);
		float const avail_depth(inner_width - 2*edge_gap + aisle_gap), avail_width(inner_len - 2*edge_gap + block_gap);
		unsigned const num_rows(avail_depth/row_spacing), num_blocks(avail_width/block_spacing); // take floor
		// use a consistent direction for center rows; face the room interior door if there is one
		vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors
		bool const sdim(!long_dim);
		bool const sdir(doorways.empty() ? rgen.rand_bool() : (room.get_center_dim(sdim) < doorways.front().get_center_dim(sdim)));
		float const dsign(sdir ? 1.0 : -1.0);
		assert(place_area.contains_cube(inner_area));

		if (num_rows > 0 && num_blocks > 0) {
			row_spacing   = avail_depth/num_rows;
			block_spacing = avail_width/num_blocks;
			aisle_gap     = row_spacing   - server_depth;
			block_gap     = block_spacing - block_width;
			float const block_start(inner_area.d[!sdim][0] + edge_gap); // include servers along the wall and their gap
			float row_pos(inner_area.d[sdim][0] + edge_gap);

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
							// TODO: add wire conduit? but what about ceiling lights?
							++num_servers;
						}
						server_pos += server_width + side_gap;
					} // for n
					block_pos += block_spacing;
				} // for b
				row_pos += row_spacing;
			} // for r
		}
	}
	if (num_servers == 0 && num_comps == 0) return 0; // both servers and computers count
	add_door_sign("Server Room", room, zval, room_id);
	return 1;
}

