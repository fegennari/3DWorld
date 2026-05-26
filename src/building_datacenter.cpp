// 3D World - Data Center Buildings
// by Frank Gennari 5/25/26

#include "function_registry.h"
#include "buildings.h"
//#include "city_model.h"

void building_t::create_datacenter_floorplan(unsigned part_id, rand_gen_t &rgen) {
	cube_t const &part(parts[part_id]);
	assert(interior->rooms.empty()); // must call this first
	cube_t room(part);
	add_assigned_room(room, part_id, RTYPE_SERVER);
}

bool building_t::add_server_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) { // for office buildings
	float const window_vspacing(get_window_vspace()), ceiling_zval(zval + get_floor_ceil_gap());
	float const server_height(0.7*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const server_width (0.3*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hwidth(0.5*server_width);
	float const server_depth (0.4*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hdepth(0.5*server_depth);
	float const comp_height  (0.2*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const min_spacing  (0.1*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const comp_hwidth(0.5*0.44*comp_height), comp_hdepth(0.5*0.9*comp_height); // fixed AR=0.44 to match the texture
	float const server_period(server_width + min_spacing), conduit_radius(0.05*server_width);
	bool const long_dim(room.dx() < room.dy()), mult_rows(is_datacenter());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // server spacing from walls
	zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_METAL); // add server room metal tile and move the effective floor up
	cube_t server, computer;
	set_cube_zvals(server,   zval, (zval + server_height));
	set_cube_zvals(computer, zval, (zval + comp_height  ));
	point center;
	unsigned num_servers(0), num_comps(0);
	vect_room_object_t &objs(interior->room_geom->objs);

	// try to line servers up against each wall wherever they fit
	for (unsigned D = 0; D < 2; ++D) {
		bool const dim(bool(D) ^ long_dim); // place along walls in long dim first
		float const room_len(place_area.get_sz_dim(dim));
		unsigned const num(room_len/server_period); // take the floor
		if (num == 0) continue; // not enough space for a server in this dim
		float const server_spacing(room_len/num);
		center[dim] = place_area.d[dim][0] + 0.5*server_spacing; // first position at half spacing

		for (unsigned n = 0; n < num; ++n, center[dim] += server_spacing) {
			set_wall_width(server, center[dim], server_hwidth, dim); // position along the wall

			for (unsigned dir = 0; dir < 2; ++dir) {
				float const dir_sign(dir ? -1.0 : 1.0);
				float const wall_pos(place_area.d[!dim][dir]);
				center[!dim] = wall_pos + dir_sign*server_hdepth;
				set_wall_width(server, center[!dim], server_hdepth, !dim); // position from the wall
				cube_t server_exp(server);
				server_exp.expand_in_dim(dim, server_hwidth); // check for more side/width spacing for doors

				// Note: overlaps_other_room_obj includes previously placed servers, so we don't have to check for intersections at the corners of rooms
				if (is_obj_placement_blocked(server_exp, room, 1) || overlaps_other_room_obj(server, objs_start)) { // no space for server; try computer instead
					set_wall_width(computer,  center[ dim], comp_hwidth, dim); // position along the wall
					set_wall_width(computer, (wall_pos + 1.2*dir_sign*comp_hdepth), comp_hdepth, !dim); // position from the wall
					if (is_obj_placement_blocked(computer, room, 1) || overlaps_other_room_obj(computer, objs_start)) continue;
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
	if (num_servers == 0 && num_comps == 0) return 0; // both servers and computers count

	if (num_servers > 0) { // add a keyboard to the master server
		unsigned const master_server(rgen.rand() % num_servers);

		for (unsigned i = objs_start, server_ix = 0; i < objs.size(); ++i) {
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
	// maybe add laptops on top of some servers, to reward the player for finding this room
	for (unsigned i = objs_start; i < objs.size(); ++i) {
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
		laptop.flags |= RO_FLAG_HANGING;
	} // for i
	// add interior rows of servers for data centers grouped into blocks
	if (mult_rows) {
		unsigned const num_per_block = 8;
		float const room_len(place_area.get_sz_dim(long_dim)), room_width(place_area.get_sz_dim(!long_dim));
		float const clearance(get_min_front_clearance_inc_people());
		float const front_clearance(1.2*clearance), back_clearance(1.0*clearance), edge_gap(max(server_depth, 1.6f*clearance));
		float const side_gap(0.02*server_width), block_width(num_per_block*server_width + (num_per_block-1)*side_gap);
		float aisle_gap(max(server_depth, 1.5f*clearance)), block_gap(max(2.0f*server_width, 1.25f*clearance));
		float row_spacing(server_depth + aisle_gap), block_spacing(block_width + block_gap);
		float const avail_depth(room_width - 2*(server_depth + edge_gap) + aisle_gap), avail_width(room_len - 2*(server_width + edge_gap) + block_gap);
		unsigned const num_rows(avail_depth/row_spacing), num_blocks(avail_width/block_spacing); // take floor
		bool const sdim(!long_dim), sdir(rgen.rand_bool()); // consistent direction for center rows
		float const dsign(sdir ? 1.0 : -1.0);

		if (num_rows > 0 && num_blocks > 0) {
			row_spacing   = avail_depth/num_rows;
			block_spacing = avail_width/num_blocks;
			aisle_gap     = row_spacing   - server_depth;
			block_gap     = block_spacing - block_width;
			float const block_start(place_area.d[!sdim][0] + server_width + edge_gap); // include servers along the wall and their gap
			float row_pos(place_area.d[sdim][0] + server_depth + edge_gap);

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
							objs.emplace_back(server, TYPE_SERVER, room_id, sdim, sdir, 0, tot_light_amt);
							// TODO: add wire conduit?
						}
						server_pos += server_width + side_gap;
					}
					block_pos += block_spacing;
				} // for b
				row_pos += row_spacing;
			} // for r
		}
	}
	add_door_sign("Server Room", room, zval, room_id);
	return 1;
}

