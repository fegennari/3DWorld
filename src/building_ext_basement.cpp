// 3D World - Building Extended Basements and Backrooms
// by Frank Gennari 03/11/2022

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
//#include "profiler.h"
#include <cfloat> // for FLT_MAX

extern int player_in_basement;
extern float DX_VAL_INV, DY_VAL_INV;
extern building_params_t global_building_params;
extern building_t const *player_building;

bool using_hmap_with_detail();


bool building_t::extend_underground_basement(rand_gen_t rgen) {
	if (!has_basement() || is_rotated() || !interior) return 0;
	//highres_timer_t timer("Extend Underground Basement"); // 540ms total
	float const height(get_window_vspace() - get_fc_thickness()); // full height of floor to avoid a gap at the top (not get_floor_ceil_gap())
	cube_t basement(get_basement());
	basement.z2() = basement.z1() + get_window_vspace(); // limit basement to the bottom floor if a parking garage
	bool dim(rgen.rand_bool()), dir(rgen.rand_bool());

	for (unsigned len = 4; len >= 2; --len) { // 100%, 75%, 50% of basement length
		for (unsigned d = 0; d < 2; ++d, dim ^= 1) { // try both dims
			for (unsigned e = 0; e < 2; ++e, dir ^= 1) { // try both dirs
				if (basement.d[dim][dir] != bcube.d[dim][dir]) continue; // wall not on the building bcube
				cube_t cand_door(place_door(basement, dim, dir, height, 0.0, 0.0, 0.25, DOOR_WIDTH_SCALE, 1, 0, rgen));
				if (cand_door.is_all_zeros()) continue; // can't place a door on this wall
				float const fc_thick(get_fc_thickness());
				set_cube_zvals(cand_door, basement.z1()+fc_thick, basement.z2()-fc_thick); // change z to span floor to ceiling for interior door
				cand_door.translate_dim(dim, (dir ? 1.0 : -1.0)*0.25*get_wall_thickness()); // zero width, centered on the door
				if (add_underground_exterior_rooms(rgen, cand_door, basement, dim, dir, 0.25*len)) return 1; // exit on success
			} // for e
		} // for d
		if (!is_house) return 0; // not large enough for office building
	} // for len
	return 0;
}

float query_min_height(cube_t const &c, float stop_at) {
	float hmin(FLT_MAX);

	if (using_tiled_terrain_hmap_tex() && !using_hmap_with_detail()) { // optimized flow when using heightmap texture
		float x1((c.x1() + X_SCENE_SIZE)*DX_VAL_INV + 0.5 + xoff2), x2((c.x2() + X_SCENE_SIZE)*DX_VAL_INV + 0.5 + xoff2);
		float y1((c.y1() + Y_SCENE_SIZE)*DY_VAL_INV + 0.5 + yoff2), y2((c.y2() + Y_SCENE_SIZE)*DY_VAL_INV + 0.5 + yoff2);

		for (float y = y1-0.5; y < y2+0.5; y += 0.5) {
			for (float x = x1-0.5; x < x2+0.5; x += 0.5) {
				min_eq(hmin, get_tiled_terrain_height_tex(x, y, 1)); // check every grid point with the X/Y range; nearest_texel=1
				if (hmin < stop_at) return hmin;
			}
		}
	}
	else { // we don't have the float heightmap here, so we have to do an expensive get_exact_zval() for each grid point
		float const x_step(0.5*DX_VAL), y_step(0.5*DY_VAL);

		for (float y = c.y1()-y_step; y < c.y2()+y_step; y += y_step) {
			for (float x = c.x1()-x_step; x < c.x2()+x_step; x += x_step) {
				min_eq(hmin, get_exact_zval(min(x, c.x2()), min(y, c.y2()))); // check every grid point with the X/Y range
				if (hmin < stop_at) return hmin;
			}
		}
	}
	return hmin;
}

struct ext_basement_room_params_t {
	vect_cube_t wall_exclude, wall_segs, temp_cubes;
	vect_extb_room_t rooms;
	vector<stairs_place_t> stairs;
};

bool building_t::is_basement_room_under_mesh_not_int_bldg(cube_t &room, building_t const *exclude) const {
	float const ceiling_zval(room.z2() - get_fc_thickness());
	if (query_min_height(room, ceiling_zval) < ceiling_zval)  return 0; // check for terrain clipping through ceiling
	// check for other buildings, including their extended basements;
	// Warning: not thread safe, since we can be adding basements to another building at the same time
	if (check_buildings_cube_coll(room, 0, 1, this, exclude)) return 0; // xy_only=0, inc_basement=1, exclude ourself
	cube_t const grid_bcube(get_grid_bcube_for_building(*this));
	assert(!grid_bcube.is_all_zeros()); // must be found
	assert(grid_bcube.contains_cube_xy(bcube)); // must contain our building
	if (!grid_bcube.contains_cube_xy(room)) return 0; // outside the grid (tile or city) bcube
	return 1;
}
bool building_t::is_basement_room_placement_valid(cube_t &room, ext_basement_room_params_t &P, bool dim, bool dir, bool *add_end_door, building_t const *exclude) const {
	float const wall_thickness(get_wall_thickness());
	cube_t test_cube(room);
	test_cube.expand_in_dim(2, -0.01*test_cube.dz()); // shrink slightly so that rooms on different floors can cross over each other
	
	if (!P.rooms.empty()) { // not the first hallway; check if too close to the basement such that the wall or trim will clip through the basement wall
		cube_t room_exp(test_cube);
		room_exp.expand_by_xy(wall_thickness + get_trim_thickness());
		if (room_exp.intersects(P.rooms.front())) return 0;
	}
	test_cube.expand_in_dim(dim, -0.1*wall_thickness); // shrink slightly to avoid intersections with our parent room
	float const room_len(room.get_sz_dim(dim)), room_width(room.get_sz_dim(!dim));
	extb_room_t *end_conn_room(nullptr);

	for (auto r = P.rooms.begin(); r != P.rooms.end(); ++r) {
		if (!r->intersects(test_cube)) continue;
		if (add_end_door == nullptr)   return 0; // no end door enabled
		if (r == P.rooms.begin())      return 0; // basement is first room - can't reconnect to it
		if (r->d[!dim][0] > room.d[!dim][0] || r->d[!dim][1] < room.d[!dim][1]) return 0; // doesn't span entire room/clips off corner - invalid
		if (r->z1() != room.z1() || r->z2() != room.z2()) return 0; // floor or ceiling zval not shared
		float const edge_pos(r->d[dim][!dir]); // intersection edge pos on other room
		float const clip_len((edge_pos - room.d[dim][!dir])*(dir ? 1.0 : -1.0));
		if (clip_len < max(room_width, 0.5f*room_len)) return 0; // clipped room length is less than room width or half unclipped room length; handles negative size
		room.d[dim][dir] = test_cube.d[dim][dir] = edge_pos; // clip room to this shorter length and add an end door; may be clipped smaller for another room
		assert(room.is_strictly_normalized());
		end_conn_room = &(*r);
		*add_end_door = 1;
	} // for r
	for (stairs_place_t const &s : P.stairs) {
		cube_t avoid(s);
		avoid.expand_in_dim(!s.dim, 0.25*s.get_sz_dim(!s.dim)); // expand to the sides to avoid placing a door too close to the stairs
		if (avoid.intersects(room)) return 0;
	}
	if (!is_basement_room_under_mesh_not_int_bldg(room, exclude)) return 0;
	if (end_conn_room) {end_conn_room->conn_bcube.assign_or_union_with_cube(room);} // include this room in our connected bcube
	return 1;
}

// add rooms to the basement that may extend outside the building's bcube
bool building_t::add_underground_exterior_rooms(rand_gen_t &rgen, cube_t const &door_bcube, cube_t const &basement, bool wall_dim, bool wall_dir, float length_mult) {
	// start by placing a hallway in ext_wall_dim/dir using interior walls;
	assert(interior);
	float const ext_wall_pos(basement.d[wall_dim][wall_dir]);
	float const hallway_len(length_mult*basement.get_sz_dim(wall_dim)), door_width(door_bcube.get_sz_dim(!wall_dim)), hallway_width(1.6*door_width);
	// Note: misnamed: hallway for houses, but backrooms for offices with parking garages
	extb_room_t hallway(basement, 0); // is_hallway=0; will likely be set below
	set_wall_width(hallway, door_bcube.get_center_dim(!wall_dim), 0.5*hallway_width, !wall_dim);
	hallway.d[wall_dim][!wall_dir] = ext_wall_pos; // flush with the exterior wall/door
	hallway.d[wall_dim][ wall_dir] = ext_wall_pos + (wall_dir ? 1.0 : -1.0)*hallway_len;
	assert(hallway.is_strictly_normalized());
	ext_basement_room_params_t P;
	if (!is_basement_room_placement_valid(hallway, P, wall_dim, wall_dir)) return 0; // try to place the hallway; add_end_door=nullptr
	// valid placement; now add the door, hallway, and connected rooms
	has_basement_door = 1;
	// Note: recording the door_stack index rather than the door index allows us to get either the first door or the first stack
	interior->ext_basement_door_stack_ix = interior->door_stacks.size();
	float const fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	P.wall_exclude.push_back(basement);
	P.wall_exclude.back().expand_in_dim(wall_dim, 1.1*get_trim_thickness()); // add slightly expanded basement to keep interior wall trim from intersecting exterior walls
	P.wall_exclude.push_back(door_bcube);
	P.wall_exclude.back().expand_in_dim(wall_dim, 2.0*wall_thickness); // make sure the doorway covers the entire wall thickness
	interior->ext_basement_hallway_room_id = interior->rooms.size();
	door_t Door(door_bcube, wall_dim, wall_dir, rgen.rand_bool());
	add_interior_door(Door, 0, 1); // open 50% of the time; is_bathroom=0, make_unlocked=1
	P.rooms.emplace_back(basement, 0);
	P.rooms.push_back(hallway);
	hallway.conn_bcube = basement; // make sure the basement is included

	if (!is_house && has_parking_garage && max_expand_underground_room(hallway, wall_dim, wall_dir, rgen)) { // office building with parking garage
		// currently, the extended basement can only be a network of connected hallways with leaf rooms, or a single large basement room (this case);
		// if we want to allow both (either a large room connected to a hallway or a large room with hallways coming off of it), we need per-room flags
		hallway.is_hallway = 0; // should already be set to 0, but this makes it more clear
		interior->has_backrooms = 1;
	}
	else { // recursively add rooms connected to this hallway in alternating dimensions
		// Note: if we get here for office buildings and global_building_params.max_ext_basement_room_depth == 0, this will only generate the hallway
		if (add_ext_basement_rooms_recur(hallway, P, door_width, !wall_dim, 1, rgen)) { // dept=1, since we already added a hallway
			end_ext_basement_hallway(hallway, P.rooms[1].conn_bcube, P, door_width, wall_dim, wall_dir, 0, rgen);
			hallway.is_hallway = 1;
		}
	}
	// place rooms, now that wall_exclude has been calculated, starting with the hallway
	cube_t wall_area(hallway);
	wall_area.d[wall_dim][!wall_dir] += (wall_dir ? 1.0 : -1.0)*0.5*wall_thickness; // move separator wall inside the hallway to avoid clipping exterior wall
	assert(P.rooms.size() >= 2); // must have at least basement and primary hallway
	interior->place_exterior_room(hallway, wall_area, fc_thick, wall_thickness, P, basement_part_ix, 0, hallway.is_hallway); // use basement part_ix; num_lights=0
	if (interior->has_backrooms) {interior->rooms.back().assign_all_to(RTYPE_BACKROOMS);} // make it backrooms

	for (auto r = P.rooms.begin()+2; r != P.rooms.end(); ++r) { // skip basement and primary hallway
		interior->place_exterior_room(*r, *r, fc_thick, wall_thickness, P, basement_part_ix, 0, r->is_hallway);
	}
	for (stairs_place_t const &stairs : P.stairs) {
		landing_t landing(stairs, 0, 0, stairs.dim, stairs.dir, stairs.add_railing, SHAPE_STRAIGHT, 0, 1, 1, 0, 1); // roof_access=0, is_at_top=1, stack_conn=1, for_ramp=0, ieb=1
		bool const against_wall[2] = {1, 1};
		landing.set_against_wall(against_wall);
		stairs_landing_base_t stairwell(landing);
		landing  .z1()  = stairs.z2() - 2.0*fc_thick;
		stairwell.z2() += 0.99*get_floor_ceil_gap(); // bottom of ceiling of upper part; must cover z-range of upper floor for AIs and room object collisions
		interior->landings.push_back(landing);
		interior->stairwells.emplace_back(stairwell, 1); // num_floors=1
	} // for stairs
	return 1;
}

bool building_t::add_ext_basement_rooms_recur(extb_room_t &parent_room, ext_basement_room_params_t &P, float door_width, bool dim, unsigned depth, rand_gen_t &rgen) {
	// add doors and other rooms along hallway; currently, all rooms are hallways
	float const parent_len(parent_room.get_sz_dim(!dim)), parent_width(parent_room.get_sz_dim(dim));
	float const end_spacing(0.75*door_width), min_length(max(4.0f*door_width, 0.5f*parent_len)), max_length(max(parent_len, 2.0f*min_length));
	float const pos_lo(parent_room.d[!dim][0] + end_spacing), pos_hi(parent_room.d[!dim][1] - end_spacing);
	if (pos_lo >= pos_hi) return 0; // not enough space to add a door
	bool const is_end_room(depth >= global_building_params.max_ext_basement_room_depth);
	float const min_width_scale(is_end_room ? 1.0 : 0.9), max_width_scale(is_end_room ? 3.0 : 1.5);
	bool was_added(0);

	for (unsigned n = 0; n < global_building_params.max_ext_basement_hall_branches; ++n) {
		for (unsigned N = 0; N < 2; ++N) { // make up to 2 tries to place this room
			bool const dir(rgen.rand_bool());
			float const conn_edge(parent_room.d[dim][dir]), room_pos(rgen.rand_uniform(pos_lo, pos_hi));
			float const room_length(rgen.rand_uniform(min_length, max_length));
			float const room_width(max(1.5f*door_width, rgen.rand_uniform(min_width_scale, max_width_scale)*parent_width));
			extb_room_t room(parent_room); // sets correct zvals
			room.d[dim][!dir] = conn_edge;
			room.d[dim][ dir] = conn_edge + (dir ? 1.0 : -1.0)*room_length;
			set_wall_width(room, room_pos, 0.5*room_width, !dim);
			bool add_end_door(0);
			if (!is_basement_room_placement_valid(room, P, dim, dir, &add_end_door)) continue; // can't place the room here
			bool const add_doors[2] = {(dir == 1 || add_end_door), (dir == 0 || add_end_door)}; // one or both
			cube_t const cur_room(add_and_connect_ext_basement_room(room, P, door_width, dim, dir, is_end_room, depth, add_doors, rgen));
			parent_room.conn_bcube.assign_or_union_with_cube(cur_room);
			was_added = 1;
			break; // done/success
		} // for N
	} // for n
	return was_added;
}

bool building_t::max_expand_underground_room(cube_t &room, bool dim, bool dir, rand_gen_t &rgen) const {
	float const floor_spacing(get_window_vspace()), step_len(1.0*floor_spacing), room_len(room.get_sz_dim(dim));
	float const room_width_min(0.5*room_len), room_sz_max(2.0*room_len);
	bool cant_expand[4] = {};
	cant_expand[2*dim + (!dir)] = 1; // can't expand back into the basement where we came from
	unsigned const start_ix(rgen.rand() & 3); // start with a random side
	cube_t exp_room(room);

	while (1) {
		bool any_valid(0);

		for (unsigned i = 0; i < 4; ++i) { // check all 4 dims
			unsigned const d((i + start_ix) & 3);
			if (cant_expand[d]) continue;
			bool const edim(d >> 1), edir(d & 1);
			cube_t exp_slice(exp_room);
			exp_slice.d[edim][ edir] += (edir ? 1.0 : -1.0)*step_len; // move the edge outward
			if (exp_slice.get_sz_dim(edim) > room_sz_max) {cant_expand[d] = 1; continue;} // too large
			exp_slice.d[edim][!edir]  = exp_room.d[edim][edir]; // shrink to zero area since we've already checked exp_room
			assert(exp_slice.is_strictly_normalized());
			if (!is_basement_room_under_mesh_not_int_bldg(exp_slice)) {cant_expand[d] = 1; continue;} // can't expand this edge any more
			exp_room.d[edim][edir] = exp_slice.d[edim][edir]; // keep edge movement
			any_valid = 1;
		} // for i
		if (!any_valid) break;
	} // end while
	assert(exp_room.contains_cube(room));
	if (exp_room.get_sz_dim(!dim) < room_width_min) return 0; // room is too narrow, make it a hallway instead
	room = exp_room;
	return 1;
}

cube_t building_t::add_ext_basement_door(cube_t const &room, float door_width, bool dim, bool dir, bool is_end_room, rand_gen_t &rgen) {
	float const fc_thick(get_fc_thickness());
	cube_t door;
	set_cube_zvals(door, room.z1()+fc_thick, room.z2()-fc_thick);
	set_wall_width(door, room.get_center_dim(!dim), 0.5*door_width, !dim);
	door.d[dim][0] = door.d[dim][1] = room.d[dim][dir]; // one end of the room
	door_t Door(door, dim, !dir, rgen.rand_bool());
	add_interior_door(Door, 0, !is_end_room); // open 50% of the time; is_bathroom=0, make_unlocked=!is_end_room
	door.expand_in_dim(dim, 2.0*get_wall_thickness());
	return door;
}
cube_t building_t::add_and_connect_ext_basement_room(extb_room_t &room, ext_basement_room_params_t &P,
	float door_width, bool dim, bool dir, bool is_end_room, unsigned depth, bool const add_doors[2], rand_gen_t &rgen)
{
	assert(room.is_strictly_normalized());
	unsigned const cur_room_ix(P.rooms.size());
	P.rooms.push_back(room);
	// add a connecting door at one or both ends
	for (unsigned d = 0; d < 2; ++d) {
		if (!add_doors[d]) continue; // no door at this end
		cube_t const door(add_ext_basement_door(room, door_width, dim, d, is_end_room, rgen));
		P.wall_exclude.push_back(door);
		room.conn_bcube.assign_or_union_with_cube(door);
	}
	// recursively add rooms connecting to this one
	bool is_hallway(0);
	if (!is_end_room) {is_hallway = add_ext_basement_rooms_recur(room, P, door_width, !dim, depth+1, rgen);}
	extb_room_t &cur_room(P.rooms[cur_room_ix]);
	cur_room.is_hallway = is_hallway; // all non-end rooms are hallways
	// for hallway with no door at the end: clip off the extra
	if (is_hallway) {end_ext_basement_hallway(cur_room, room.conn_bcube, P, door_width, dim, dir, depth+1, rgen);} // Note: may invalidate cur_room reference
	return P.rooms[cur_room_ix];
}

void building_t::end_ext_basement_hallway(extb_room_t &room, cube_t const &conn_bcube, ext_basement_room_params_t &P,
	float door_width, bool dim, bool dir, unsigned depth, rand_gen_t &rgen)
{
	room.conn_bcube.assign_or_union_with_cube(conn_bcube); // combine conns from child rooms and end connected rooms
	float const stairs_len(3.0*door_width), stairs_end(room.d[dim][dir]), dsign(dir ? 1.0 : -1.0);
	
	if (depth < global_building_params.max_ext_basement_room_depth && dsign*(room.d[dim][dir] - room.conn_bcube.d[dim][dir]) > stairs_len /*&& rgen.rand_bool()*/) {
		float const ceil_below(room.z1()), floor_below(ceil_below - room.dz());

		if (floor_below > get_max_sea_level()) {
			// connect downward with stairs; we know that there aren't any other rooms/doors coming off the end that could be blocked
			float const fc_thick(get_fc_thickness()), stairs_start(stairs_end - dsign*stairs_len);
			float const hall_below_len(max(4.0f*door_width, rgen.rand_uniform(0.5, 1.5)*room.get_sz_dim(dim)));
			extb_room_t hall_below(room, 0, 1); // copy !dim values; is_hallway will be set later; has_stairs=1
			set_cube_zvals(hall_below, floor_below, ceil_below);
			hall_below.d[dim][!dir] = stairs_start;
			hall_below.d[dim][ dir] = stairs_end + dsign*hall_below_len;
			assert(hall_below.is_strictly_normalized());
			bool add_end_door(0);

			if (is_basement_room_placement_valid(hall_below, P, dim, dir, &add_end_door)) {
				// create stairs
				cube_t stairs(room); // copy room.d[dim][dir] (far end/bottom of stairs)
				set_cube_zvals(stairs, (floor_below + fc_thick), (room.z1() + fc_thick));
				stairs.d[dim][!dir] = stairs_start; // near end/top of stairs
				bool const add_railing(rgen.rand_bool()); // 50% of the time
				float const wall_half_thick(0.5*get_wall_thickness());
				stairs.expand_in_dim(!dim, -(add_railing ? 2.0 : 1.0)*wall_half_thick); // shrink on the sides; more if there are railings
				stairs.expand_in_dim( dim, -(wall_half_thick + get_trim_thickness()));  // shrink on the ends
				assert(!stairs.intersects_xy(room.conn_bcube));
				P.stairs.emplace_back(stairs, dim, !dir, add_railing);
				hall_below.conn_bcube = stairs; // must include the stairs
				room.has_stairs = 1;
				// add the room
				bool const add_doors[2] = {(dir == 0 && add_end_door), (dir == 1 && add_end_door)}; // at most one at the end
				add_and_connect_ext_basement_room(hall_below, P, door_width, dim, dir, 0, depth, add_doors, rgen); // is_end_room=0
				return; // done
			}
		}
	}
	room.clip_hallway_to_conn_bcube(dim); // else clip
}

void extb_room_t::clip_hallway_to_conn_bcube(bool dim) { // clip off the unconnected length at the end of the hallway
	if (!is_hallway) return;
	assert(conn_bcube.is_strictly_normalized());
	assert(conn_bcube.intersects(*this));
	max_eq(d[dim][0], conn_bcube.d[dim][0]);
	min_eq(d[dim][1], conn_bcube.d[dim][1]);
	assert(is_strictly_normalized());
}

void building_interior_t::place_exterior_room(extb_room_t const &room, cube_t const &wall_area, float fc_thick, float wall_thick, ext_basement_room_params_t &P,
	unsigned part_id, unsigned num_lights, bool is_hallway, unsigned is_building_conn, unsigned wall_skip_dim, unsigned thin_wall_dir)
{
	assert(room.is_strictly_normalized());
	bool const long_dim(room.dx() < room.dy());
	if (num_lights == 0) {num_lights = max(1U, min(8U, (unsigned)round_fp(0.33*room.get_sz_dim(long_dim)/room.get_sz_dim(!long_dim))));} // auto calculate num_lights
	room_t Room(room, part_id, num_lights, is_hallway);
	Room.interior = is_building_conn + 2; // mark as extended basement, or possibly connecting room between two buildings if is_building_conn == 1|2
	if (room.has_stairs) {Room.has_stairs = 1;} // stairs on the first/only floor
	if (is_hallway) {Room.assign_all_to(RTYPE_HALL, 0);} // initially all hallways; locked=0
	rooms.push_back(Room);
	cube_t ceiling(room), floor(room);
	ceiling.z1() = room.z2() - fc_thick;
	floor  .z2() = room.z1() + fc_thick;
	subtract_cubes_from_cube(ceiling, P.stairs, P.wall_segs, P.temp_cubes, 2); // cut out stairs; zval_mode=2 (check for zval overlap)
	vector_add_to(P.wall_segs, ceilings);
	subtract_cubes_from_cube(floor,   P.stairs, P.wall_segs, P.temp_cubes, 2); // cut out stairs; zval_mode=2 (check for zval overlap)
	vector_add_to(P.wall_segs, floors);
	basement_ext_bcube.assign_or_union_with_cube(room);
	// add walls; Note: two adjoining rooms may share overlapping walls
	float const wall_half_thick(0.5*wall_thick);

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			float half_thick(wall_half_thick);
			
			if (dim == wall_skip_dim) {
				if (dir == thin_wall_dir) {half_thick *= 0.5;} // this wall needed for shadows, but make it thinner to avoid Z-fighting
				else continue; // no walls in this dim/dir
			}
			cube_t wall(wall_area);
			set_wall_width(wall, wall_area.d[dim][dir], half_thick, dim);
			set_cube_zvals(wall, floor.z2(), ceiling.z1());
			if (bool(dim) != long_dim) {wall.expand_in_dim(!dim, -half_thick);} // remove the overlaps at corners in the long dim (house exterior wall dim)
			// remove overlapping walls from shared rooms
			unsigned const wall_exclude_sz(P.wall_exclude.size());
			assert(extb_walls_start[dim] <= walls[dim].size());

			for (auto w = walls[dim].begin()+extb_walls_start[dim]; w != walls[dim].end(); ++w) { // check all prev added ext basement walls in this dim
				if (w->d[ dim][0] == wall.d[ dim][0] && w->d[ dim][1] == wall.d[dim][1] &&
					w->d[!dim][0] <  wall.d[!dim][1] && w->d[!dim][1] >  wall.d[dim][0] &&
					w->z1() == wall.z1() && w->z2() == wall.z2()) {P.wall_exclude.push_back(*w);} // same width and height, overlap in length
			}
			subtract_cubes_from_cube(wall, P.wall_exclude, P.wall_segs, P.temp_cubes, 2); // cut out doorways, etc.; zval_mode=2 (check for zval overlap)
			vector_add_to(P.wall_segs, walls[dim]);
			P.wall_exclude.resize(wall_exclude_sz); // remove the wall_exclude cubes we just added
		} // for dir
	} // for dim
}

bool try_place_wall(cube_t const &place_area, cube_t const &wall_area, cube_t const &cent_area, bool dim, float len_min, float len_max, float half_thick, float min_gap,
	vect_cube_t &walls, vect_cube_t &blockers, rand_gen_t &rgen)
{
	assert(len_min < len_max);
	float const wall_len   (rgen.rand_uniform(len_min, len_max));
	float const wall_pos   (rgen.rand_uniform(wall_area.d[ dim][0], wall_area.d[ dim][1])); // position of wall centerline
	float const wall_center(rgen.rand_uniform(cent_area.d[!dim][0], cent_area.d[!dim][1])); // center of wall
	cube_t wall(wall_area); // copy zvals
	set_wall_width(wall, wall_pos,    half_thick,    dim);
	set_wall_width(wall, wall_center, 0.5*wall_len, !dim);
	// wall can extend outside the room, and we generally want it to end at the room edge in some cases
	if (wall.d[!dim][0] < place_area.d[!dim][0]) {wall.d[!dim][0] = place_area.d[!dim][0];} // clip to room bounds
	else if (wall.d[!dim][0] - place_area.d[!dim][0] < min_gap) {wall.d[!dim][0] += min_gap;} // force min_gap from wall
	if (wall.d[!dim][1] > place_area.d[!dim][1]) {wall.d[!dim][1] = place_area.d[!dim][1];} // clip to room bounds
	else if (place_area.d[!dim][1] - wall.d[!dim][1] < min_gap) {wall.d[!dim][1] -= min_gap;} // force min_gap from wall
	if (wall.get_sz_dim(!dim) < len_min) return 0; // too short
	if (has_bcube_int(wall, blockers))   return 0; // too close to a previous wall in this dim
	walls.push_back(wall);
	wall.expand_by_xy(min_gap); // require a gap around this wall - but can still have intersections in the other dim
	blockers.push_back(wall);
	return 1;
}
void partition_cubes_into_conn_groups(vect_cube_t const &cubes, vector<vect_cube_t> &groups, float pad=0.0) {
	groups.clear();
	vect_cube_t ungrouped(cubes); // all cubes start ungrouped

	while (!ungrouped.empty()) {
		vect_cube_t group;
		group.push_back(ungrouped.back());
		ungrouped.pop_back();

		for (unsigned ix = 0; ix < group.size(); ++ix) { // process each cube added to the group
			cube_t cur(group[ix]);
			cur.expand_by_xy(pad);

			for (unsigned i = 0; i < ungrouped.size(); ++i) {
				cube_t &cand(ungrouped[i]);
				if (!cand.intersects(cur)) continue; // includes adjacency
				group.push_back(cand);
				cand = ungrouped.back(); // remove cand by replacing it with the last ungrouped cube
				ungrouped.pop_back();
				--i;
			} // for i
		} // for ix
		groups.push_back(group);
	} // while
}

struct cmp_cube_x1_y1 {
	bool operator()(cube_t const &a, cube_t const &b) const {return ((a.x1() == b.x1()) ? (a.y1() < b.y1()) : (a.x1() < b.x1()));}
};
void invert_walls(cube_t const &room, vect_cube_t const walls[2], vect_cube_t &space, float pad=0.0) {
	space.clear();
	space.push_back(room);
	
	for (unsigned d = 0; d < 2; ++d) {
		for (cube_t const &wall : walls[d]) {
			cube_t sub(wall);
			sub.expand_by_xy(pad);
			subtract_cube_from_cubes(sub, space);
		}
	}
	if (space.empty()) return;
	// max merge adjacent cubes in Y
	sort(space.begin(), space.end(), cmp_cube_x1_y1());
	auto i(space.begin()+1), o(i); // skip first

	for (; i != space.end(); ++i) {
		cube_t &c(*(o-1)); // last output cube
		if (c.x1() == i->x1() && c.x2() == i->x2() && c.y2() == i->y1()) {c.y2() = i->y2();} // extend in +Y
		else {*o++ = *i;} // keep this cube
	} // for i
	space.erase(o, space.end());
}
void resize_cubes_xy(vect_cube_t &cubes, float val) { // val can be positive or negative
	for (cube_t &c : cubes) {c.expand_by_xy(val);}
}
void transpose_cube_xy(cube_t &c) {
	swap(c.x1(), c.y1());
	swap(c.x2(), c.y2());
}
void transpose_cubes_xy(vect_cube_t &cubes) {
	for (cube_t &c : cubes) {transpose_cube_xy(c);}
}

struct group_range_t {
	unsigned gix;
	float lo, hi;
	group_range_t(unsigned gix_, float lo_, float hi_) : gix(gix_), lo(lo_), hi(hi_) {}
	void update(float lo_, float hi_) {min_eq(lo, lo_); max_eq(hi, hi_);}
	float len() const {return (hi - lo);}
};
struct cube_by_center_dim_descending {
	unsigned d;
	cube_by_center_dim_descending(unsigned dim) : d(dim) {}
	bool operator()(cube_t const &a, cube_t const &b) const {return (b.get_center_dim(d) < a.get_center_dim(d));}
};
struct cube_by_xy_dim {
	unsigned d;
	cube_by_xy_dim(unsigned dim) : d(dim) {}
	bool operator()(cube_t const &a, cube_t const &b) const {return ((a.d[d][0] == b.d[d][0]) ? (a.d[!d][0] < b.d[!d][0]) : (a.d[d][0] < b.d[d][0]));}
};

void building_t::add_backrooms_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, vect_cube_t &rooms_to_light) {
	//highres_timer_t timer("Add Backrooms Objs"); // up to ~2ms
	assert(has_room_geom());
	float const floor_spacing(get_window_vspace()), wall_thickness(1.2*get_wall_thickness()), wall_half_thick(0.5*wall_thickness); // slightly thicker than regular walls
	float const ceiling_z(zval + get_floor_ceil_gap()); // Note: zval is at floor level, not at the bottom of the room
	//float const tot_light_amt(room.light_intensity /*+ floor_spacing*room.get_light_amt()*/); // ???
	float const tot_light_amt(0.75);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());
	if (interior->room_geom->backrooms_start == 0) {interior->room_geom->backrooms_start = objs_start;}

	// find the shared wall with the basement/parking garage and calculate true room bounds
	bool sw_dim(0), sw_dir(0), adj_found(0);
	cube_t const &parking_garage(get_basement());

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (room.d[dim][dir] == parking_garage.d[dim][!dir]) {sw_dim = dim; sw_dir = dir; adj_found = 1;}
		}
	}
	assert(adj_found);
	float const shared_extend((sw_dir ? -1.0 : 1.0)*0.5*get_wall_thickness()); // account for extra shift applied to shared wall to keep it from clipping into the basement/PG
	room_t true_room(room);
	true_room.d[sw_dim][sw_dir] += shared_extend;
	vector3d const sz(true_room.get_size());
	vect_cube_t blockers_per_dim[2], walls_per_dim[2];

	// find entrance door and add it to wall blockers
	assert((unsigned)interior->ext_basement_door_stack_ix < interior->door_stacks.size());
	auto const &ent_door(interior->door_stacks[interior->ext_basement_door_stack_ix]);
	blockers_per_dim[!ent_door.dim].push_back(ent_door.get_clearance_bcube());

	// add random interior walls to create an initial maze
	float const doorway_width(get_doorway_width()), doorway_hwidth(0.5*doorway_width), min_gap(1.2*doorway_width);
	cube_t const place_area(get_walkable_room_bounds(true_room));
	cube_t wall_area(place_area);
	wall_area.expand_by_xy(-min_gap);
	if (min(wall_area.dx(), wall_area.dy()) < 2.0*floor_spacing) return; // room too small to place walls - shouldn't happen
	set_cube_zvals(wall_area, zval, ceiling_z);
	float const min_side(min(sz.x, sz.y)), area(sz.x*sz.y);
	float const wall_len_min(1.0*floor_spacing), wall_len_max(max(0.25f*min_side, 1.5f*wall_len_min)), wall_len_avg(0.5*(wall_len_min + wall_len_max)); // min_gap/wall_len_min = 0.6
	unsigned const num_walls(round_fp(rgen.rand_uniform(1.6, 2.0)*area/(wall_len_avg*wall_len_avg))); // add a bit of random density variation
	colorRGBA const wall_color(WHITE); // to match exterior walls, ceilings, and floors

	for (unsigned n = 0; n < num_walls; ++n) {
		for (unsigned m = 0; m < 10; ++m) { // 10 tries to place a wall
			bool const dim(rgen.rand_bool()); // should this be weighted by the room aspect ratio?
			if (try_place_wall(place_area, wall_area, place_area, dim, wall_len_min, wall_len_max, wall_half_thick, min_gap, walls_per_dim[dim], blockers_per_dim[dim], rgen)) break;
		}
	}

	// find long lines of sight and add extra walls to block them
	float const max_space_factor = 0.5; // relative to room size
	vect_cube_t space, extra_walls;

	for (unsigned dim = 0; dim < 2; ++dim) {
		// Note: wall ends may be moved slightly to snap to orthogonal walls, so this test isn't perfectly accurate
		if (dim) { // transpose cubes so that max merge is in Y rather than X
			vect_cube_t walls_transpose[2] = {walls_per_dim[1], walls_per_dim[0]};
			for (unsigned d = 0; d < 2; ++d) {transpose_cubes_xy(walls_transpose[d]);}
			cube_t place_area_transpose(place_area);
			transpose_cube_xy(place_area_transpose);
			invert_walls(place_area_transpose, walls_transpose, space); // no padding
			transpose_cubes_xy(space); // transpose back
		}
		else {invert_walls(place_area, walls_per_dim, space);} // no padding
		float const max_space(max_space_factor*sz[dim]);
		
		for (cube_t &s : space) {
			float const space_len(s.get_sz_dim(dim));
			if (space_len < max_space)         continue; // small enough
			if (has_bcube_int(s, extra_walls)) continue; // covered (at least partially) by a previously placed extra wall
			cube_t cent_area(s);
			cent_area.expand_in_dim(dim, -0.375*space_len); // restrict to central 25%
			set_cube_zvals(cent_area, wall_area.z1(), wall_area.z2()); // set zvals the same as walls
			float const space_width(s.get_sz_dim(!dim)), min_len(max(wall_len_min, space_width)), max_len(max(wall_len_max, 2.0f*space_width)); // must cross the space

			for (unsigned m = 0; m < 10; ++m) { // 10 tries to place a wall
				if (try_place_wall(place_area, cent_area, cent_area, dim, min_len, max_len, wall_half_thick, min_gap, walls_per_dim[dim], blockers_per_dim[dim], rgen)) {
					extra_walls.push_back(walls_per_dim[dim].back()); // record in case this blocks other long spaces in this loop
					break; // done
				}
			}
		} // for s
		extra_walls.clear();
	} // for dim
	// what if we instead calculate the shortest path from the entrance door to the far wall and add walls until it's some min length? is the runtime acceptable?

	// shift wall ends to remove small gaps and stubs
	float const wall_end_ext(doorway_width); // needed to handle two orthogonal walls with nearby corners but no projection

	for (unsigned dim = 0; dim < 2; ++dim) {
		vect_cube_t &walls(walls_per_dim[dim]);

		for (cube_t &wall : walls) {
			for (unsigned d = 0; d < 2; ++d) { // check each end
				float &val(wall.d[!dim][d]);
				if (val == place_area.d[!dim][d]) continue; // at exterior wall, skip

				for (cube_t &w : walls_per_dim[!dim]) {
					if (wall.d[dim][0] > w.d[dim][1] || wall.d[dim][1] < w.d[dim][0]) { // no projection
						if (wall.d[dim][0] > w.d[dim][1]+wall_end_ext || wall.d[dim][1] < w.d[dim][0]-wall_end_ext) continue; // no extended projection
						// handle nearby corners by moving wall away from the corner if needed; doesn't break because the edge may be moved again in the else case
						if      (d == 0 && val > w.d[!dim][1]) {max_eq(val, w.d[!dim][1]+min_gap);} // ensure left  space between walls is at least min_gap
						else if (d == 1 && val < w.d[!dim][0]) {min_eq(val, w.d[!dim][0]-min_gap);} // ensure right space between walls is at least min_gap
					}
					else {
						// if we already aligned <w> to <wall> then we now have a corner; extend val to the far edge of <w> to fill in the corner area
						bool const is_corner(wall.d[dim][0] == w.d[dim][1] || wall.d[dim][1] == w.d[dim][0]);
						float const edge_pos(w.d[!dim][bool(!d) ^ is_corner]);
						if (fabs(val - edge_pos) < min_gap) {val = edge_pos; break;}
					}
				} // for w
			} // for d
			if (wall.get_sz_dim(!dim) < wall_len_min) {wall = cube_t();} // too short, drop
		} // for wall
		walls.erase(std::remove_if(walls.begin(), walls.end(), [](cube_t const &c) {return c.is_all_zeros();}), walls.end());
	} // for dim

	// find areas of empty space
	float const pad(wall_half_thick), nav_pad(doorway_hwidth); // min radius for player and AI navigation
	vector<vect_cube_t> space_groups;
	invert_walls(place_area, walls_per_dim, space, nav_pad);
	resize_cubes_xy(space, nav_pad); // restore padding (under-over)
	partition_cubes_into_conn_groups(space, space_groups, pad);
	vect_cube_t small_rooms, door_keepout;

#if 0 // TESTING
	rand_gen_t rgen2(rgen);
	for (vect_cube_t const &g : space_groups) {
		colorRGBA const color(rgen2.rand_float(), rgen2.rand_float(), rgen2.rand_float());
		for (cube_t const &c : g) {objs.emplace_back(c, TYPE_DBG_SHAPE, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_BACKROOM), 1.0, SHAPE_CUBE, color);}
	}
#endif

	// Add doorways + doors to guarantee full connectivity using space
	if (space_groups.size() > 1) { // multiple disconnected sub-graphs
		float const door_min_spacing(0.5*doorway_width), min_shared_edge(doorway_width + 2*wall_thickness); // allow space for door frame
		vector<group_range_t> adj;
		vect_cube_t doors_to_add, walls_to_add, all_doors;
		set<pair<unsigned, unsigned>> connected;
		bool const first_dim(rgen.rand_bool()); // mix it up to avoid favoring one dim

		for (unsigned d = 0; d < 2; ++d) {
			unsigned const dim(bool(d) ^ first_dim);
			walls_to_add.clear();
			vect_cube_t &walls(walls_per_dim[dim]);

			for (cube_t &wall : walls) {
				cube_t query(wall);
				query.expand_in_dim(dim, pad); // increase wall width to overlap space cubes
				adj.clear();

				for (unsigned gix = 0; gix < space_groups.size(); ++gix) {
					bool hit(0);

					for (cube_t const &s : space_groups[gix]) {
						if (!s.intersects(query)) continue;
						float const lo(max(query.d[!dim][0], s.d[!dim][0])), hi(min(query.d[!dim][1], s.d[!dim][1]));
						if (!hit) {adj.emplace_back(gix, lo, hi);} else {adj.back().update(lo, hi);}
						hit = 1;
					}
					if (hit && adj.back().len() < min_shared_edge) {adj.pop_back();} // remove if too short
				} // for gix
				if (adj.size() <= 1) continue; // not adjacent to multiple space group
				//objs.emplace_back(query, TYPE_DBG_SHAPE, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_BACKROOM), 1.0, SHAPE_CUBE, RED); // TESTING
				doors_to_add.clear();

				// test every pair of groups
				for (unsigned i = 0; i < adj.size(); ++i) {
					for (unsigned j = i+1; j < adj.size(); ++j) {
						group_range_t const &g1(adj[i]), &g2(adj[j]);
						float const lo(max(g1.lo, g2.lo)), hi(min(g1.hi, g2.hi));
						if (hi - lo < min_shared_edge) continue; // not enough space to add a door
						pair<unsigned, unsigned> const ix_pair(min(g1.gix, g2.gix), max(g1.gix, g2.gix));
						if (connected.find(ix_pair) != connected.end()) continue; // these two groups were already connected
						float const vmin(lo + doorway_hwidth), vmax(hi - doorway_hwidth);
							
						// add a doorway here if possible
						for (unsigned n = 0; n < 10; ++n) { // make up to 10 tries
							float const door_pos(rgen.rand_uniform(vmin, vmax));
							cube_t door_cand(wall);
							set_wall_width(door_cand, door_pos, doorway_hwidth, !dim);
							// check for collisions with walls in other dims and prev placed doors; shouldn't collide with walls in this dim due to checks during wall placement
							cube_t wall_test_cube(door_cand);
							wall_test_cube.expand_in_dim(dim, doorway_width); // make room for the door to open
							if (has_bcube_int(wall_test_cube, walls_per_dim[!dim])) continue;
							cube_t door_test_cube(wall);
							door_test_cube.expand_in_dim(!dim, door_min_spacing);
							if (has_bcube_int(door_test_cube, all_doors)) continue; // check all doors, even ones on other walls in case there are two overlapping doors
							doors_to_add.push_back(door_cand);
							all_doors   .push_back(door_cand);
							connected.insert(ix_pair); // mark these two space groups as being connected
							break;
						} // for n
					} // for j
				} // for i
				sort(doors_to_add.begin(), doors_to_add.end(), cube_by_center_dim_descending(!dim)); // sort in dim !dim, high to low

				for (cube_t const &door : doors_to_add) {
					assert(door.is_strictly_normalized());
					// select an open direction that doesn't block another door
					bool open_dir(rgen.rand_bool());
					cube_t open_area(door);
					open_area.d[dim][open_dir] += (open_dir ? 1.0 : -1.0)*doorway_width;
					if (has_bcube_int(open_area, door_keepout)) {open_dir ^= 1;} // swap the open direction; hopefully we don't block a different door now
					// cut the doorway into the wall and add the door
					cube_t wall2;
					bool const make_unlocked = 1; // makes exploring easier
					bool const make_closed   = 1; // makes it easier to tell if the door has been used
					remove_section_from_cube_and_add_door(wall, wall2, door.d[!dim][0], door.d[!dim][1], !dim, open_dir, 0, make_unlocked, make_closed); // is_bathroom=0
					walls_to_add.push_back(wall2); // keep high side as it won't be used with any other doors
					// add a blocker so that no ceiling lights are placed in the path of this door
					cube_t blocker(door);
					blocker.d[dim][open_dir] += (open_dir ? 1.0 : -1.0)*doorway_width; // add clearance in front
					objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, 0, (RO_FLAG_INVIS | RO_FLAG_BACKROOM), tot_light_amt);
					// keep other doors from opening into this door's space
					cube_t keepout(door);
					keepout.expand_in_dim(dim, doorway_width); // don't block either end
					door_keepout.push_back(keepout);
				} // for door
			} // for wall
			vector_add_to(walls_to_add, walls);
		} // for d
		for (unsigned gix = 0; gix < space_groups.size(); ++gix) {
			vect_cube_t const &group(space_groups[gix]);
			if (group.size() != 1) continue;
			cube_t sub_room(group.front());
			sub_room.intersect_with_cube(true_room); // can't go outside the backrooms (under-over can move exterior walls)

			for (auto const &c : connected) { // add if it was connected with a door
				if (c.first == gix || c.second == gix) {small_rooms.push_back(sub_room); break;}
			}
		} // for gix
	}

	// Add some random pillars in large open spaces
	float const pillar_grid_step(2.5*min_gap), pillar_grid_exp(0.75*pillar_grid_step), pillar_hwidth(0.07*floor_spacing), merge_wall_min(1.0*pillar_grid_step);
	unsigned const xdiv(ceil(wall_area.dx()/pillar_grid_step)), ydiv(ceil(wall_area.dy()/pillar_grid_step));
	float const xstep(wall_area.dx()/(xdiv-1)), ystep(wall_area.dy()/(ydiv-1));
	vect_cube_t big_space;
	cube_t grid(wall_area); // copy zvals

	for (unsigned y = 0; y < ydiv; ++y) {
		for (unsigned x = 0; x < xdiv; ++x) {
			set_wall_width(grid, (wall_area.y1() + y*ystep), pillar_grid_exp, 1);
			set_wall_width(grid, (wall_area.x1() + x*xstep), pillar_grid_exp, 0);	
			if (!place_area.contains_cube_xy(grid) || has_bcube_int(grid, walls_per_dim[0]) || has_bcube_int(grid, walls_per_dim[1])) continue;
			big_space.push_back(grid);
		} // for x
	} // for y
	partition_cubes_into_conn_groups(big_space, space_groups, pad);
	cube_t pillar(wall_area); // copy zvals

	for (vect_cube_t &group : space_groups) {
		if (rgen.rand_float() < 0.4) continue; // no pillars 40% of the time
		cube_t prev_space, cur_row;
		sort(group.begin(), group.end(), cube_by_xy_dim(rgen.rand_bool())); // make adjacent space cubes adjacent in order, with dim chosen randomly

		for (auto s = group.begin(); s != group.end(); ++s) {
			for (unsigned d = 0; d < 2; ++d) {set_wall_width(pillar, s->get_center_dim(d), pillar_hwidth, d);}
			objs.emplace_back(pillar, TYPE_PG_PILLAR, room_id, 0, 0, RO_FLAG_BACKROOM, tot_light_amt, SHAPE_CUBE, wall_color); // dim=0, dir=0
			// maybe merge rows of adjacent pillars into a single wall
			cube_t merge_cand;

			if (!s->intersects(prev_space) || ((pillar.x1() != cur_row.x1() || pillar.x2() != cur_row.x2()) && (pillar.y1() != cur_row.y1() || pillar.y2() != cur_row.y2()))) {
				// different row, gap, or first pillar
				merge_cand = cur_row;
				cur_row    = pillar; // seed for next row
			}
			else { // extend row
				cur_row.union_with_cube(pillar);
			}
			if (s+1 == group.end()) {merge_cand = cur_row;} // last pillar
			prev_space = *s;

			if (!merge_cand.is_all_zeros() && rgen.rand_float() < 0.5) { // add wall 50% of the time
				for (unsigned dim = 0; dim < 2; ++dim) {
					if (merge_cand.get_sz_dim(!dim) < merge_wall_min) continue; // too short to create a wall; will get here in at least one dim
					cube_t wall(merge_cand);
					wall.expand_by_xy(wall_half_thick - pillar_hwidth); // shrink
					walls_per_dim[dim].push_back(wall);
				}
			}
		} // for d
	} // for g

	// add walls
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (cube_t &wall : walls_per_dim[dim]) {
			objs.emplace_back(wall, TYPE_PG_WALL, room_id, dim, 0, RO_FLAG_BACKROOM, tot_light_amt, SHAPE_CUBE, wall_color); // dir=0
		}
		vector_add_to(walls_per_dim[dim], interior->room_geom->pgbr_walls[dim]); // store walls for occlusion and door opening checks
	}
	
	// Add vents, light switch, and outlets (on all walls - or should it be only on the wall adjacent to building or next to the door?)
	add_light_switches_to_room(rgen, true_room, zval, room_id, objs_start, 0, 1); // is_ground_floor=0, is_basement=1
	add_outlets_to_room       (rgen, true_room, zval, room_id, objs_start, 0, 1); // is_ground_floor=0, is_basement=1
	add_wall_vent_to_room     (rgen, true_room, zval, room_id, objs_start, 0   ); // check_for_ducts=0
	rgen.rand_bool(); // mix up in between
	add_ceil_vent_to_room     (rgen, true_room, zval, room_id, objs_start); // add a ceiling vent as well

	// Add occasional random items/furniture
	// TODO: chairs, boxes, crates, balls, fire extinguisher, etc.
	
	// Make small rooms with doors bathrooms, etc.
	for (cube_t const &r : small_rooms) {
		unsigned num_doors(0);

		for (cube_t const &dk : door_keepout) {
			if (dk.intersects_xy(r)) {++num_doors;}
		}
		if (num_doors == 0) continue; // not connected with a door, not reachable, skip (and don't need a light either)
		room_t sub_room(room);
		sub_room.copy_from(r); // keep flags, copy cube
		sub_room.interior = 1; // treated as basement but not extended basement (no wall padding)
		unsigned const sub_objs_start(objs.size()); // no objects have been placed in this sub-room yet

		if (num_doors == 1) { // only make bathroom if there's a single door
			float floor_zval(zval); // may be modified below, but otherwise unused
			unsigned const floor_ix(0); // pass this in, or always zero?
			unsigned added_bathroom_objs_mask(0); // unused
			add_bathroom_objs(rgen, sub_room, floor_zval, room_id, tot_light_amt, sub_objs_start, floor_ix, 1, added_bathroom_objs_mask); // is_basement=1
			room.has_mirror |= sub_room.has_mirror;
		}
		else { // 2 or more rooms
			if (rgen.rand_bool()) {add_furnace_to_room(rgen, sub_room, zval, room_id, tot_light_amt, sub_objs_start);}
			else {add_water_heaters(rgen, sub_room, zval, room_id, tot_light_amt, sub_objs_start, 1);} // single_only=1
		}
		rooms_to_light.push_back(r);
	} // for r
}

void building_t::add_missing_backrooms_lights(rand_gen_t rgen, float zval, unsigned room_id, unsigned objs_start, unsigned lights_start,
	vect_cube_t const &rooms_to_light, light_ix_assign_t &light_ix_assign)
{
	if (rooms_to_light.empty()) return; // nothing to do
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const lights_end(objs.size());
	assert(lights_start < lights_end);
	if (lights_start == lights_end) return; // no lights; error?
	room_object_t const ref_light(objs[lights_start]); // use this as a reference for the light size/shape/color
	point const ref_light_center(ref_light.get_cube_center());
	vect_cube_t to_add;

	for (cube_t const &r : rooms_to_light) {
		bool has_light(0);

		for (auto i = objs.begin()+lights_start; i != objs.begin()+lights_end; ++i) {
			assert(i->type == TYPE_LIGHT);
			if (i->intersects(r)) {has_light = 1; break;}
		}
		if (has_light) continue;
		room_object_t light(ref_light);
		light += vector3d((r.xc() - ref_light_center.x), (r.yc() - ref_light_center.y), 0.0);
		bool const room_dim(r.dx() < r.dy()); // longer room dim
		to_add.clear();
		try_place_light_on_ceiling(light, r, room_dim, get_fc_thickness(), 1, 0, 1, 1, objs_start, to_add, rgen); // or wall light?

		for (cube_t const &L : to_add) { // should be size 1
			light.copy_from(L);
			light.obj_id = light_ix_assign.get_ix_for_light(light);
			objs.push_back(light);
		}
	} // for r
}

bool building_room_geom_t::cube_int_backrooms_walls(cube_t const &c) const { // used for door opening collision checks
	// no dim is passed in, so we check both dims; includes parking garage walls, which we can probably ignore, but they should be small in size
	for (unsigned d = 0; d < 2; ++d) {
		if (has_bcube_int(c, pgbr_walls[d])) return 1;
	}
	return 0;
}

cube_t building_t::get_bcube_inc_extensions() const {
	cube_t ret(bcube);
	if (has_ext_basement()) {ret.union_with_cube(interior->basement_ext_bcube);}
	return ret;
}
cube_t building_t::get_full_basement_bcube() const {
	assert(has_basement());
	cube_t ret(get_basement());
	if (has_ext_basement()) {ret.union_with_cube(interior->basement_ext_bcube);}
	return ret;
}
room_t const &building_t::get_ext_basement_hallway() const {
	assert(interior);
	assert(interior->ext_basement_hallway_room_id >= 0);
	return *interior->ext_basement_rooms_start();
}

vector<room_t>::const_iterator building_interior_t::ext_basement_rooms_start() const {
	if (ext_basement_hallway_room_id < 0) return rooms.end(); // no ext basement rooms
	assert((unsigned)ext_basement_hallway_room_id < rooms.size());
	return rooms.begin() + ext_basement_hallway_room_id;
}
bool building_interior_t::point_in_ext_basement_room(point const &pos, float expand) const {
	if (ext_basement_hallway_room_id < 0)                 return 0; // no ext basement rooms
	if (!basement_ext_bcube.contains_pt_exp(pos, expand)) return 0;

	for (auto r = ext_basement_rooms_start(); r != rooms.end(); ++r) {
		if (r->contains_pt_exp(pos, expand)) return 1;
	}
	return 0;
}
// returns true if cube is completely contained in any single room
bool building_interior_t::cube_in_ext_basement_room(cube_t const &c, bool xy_only) const {
	if (ext_basement_hallway_room_id < 0)        return 0; // no ext basement rooms
	if (!basement_ext_bcube.contains_cube_xy(c)) return 0;

	for (auto r = ext_basement_rooms_start(); r != rooms.end(); ++r) {
		if (xy_only ? r->contains_cube_xy(c) : r->contains_cube(c)) return 1;
	}
	return 0;
}
door_t const &building_interior_t::get_ext_basement_door() const {
	assert(ext_basement_door_stack_ix >= 0 && (unsigned)ext_basement_door_stack_ix < door_stacks.size());
	unsigned const door_ix(door_stacks[ext_basement_door_stack_ix].first_door_ix);
	return get_door(door_ix);
}


// *** Code to join exterior basements of two nearby buildings ***

void populate_params_from_building(building_interior_t const &bi, ext_basement_room_params_t &P) {
	if (!P.rooms.empty()) return; // already populated
	for (auto r = bi.rooms.begin()+bi.ext_basement_hallway_room_id+1; r != bi.rooms.end(); ++r) {P.rooms.emplace_back(*r, r->is_hallway, r->has_stairs);}
	for (auto s = bi.stairwells.begin(); s != bi.stairwells.end(); ++s) {P.stairs.emplace_back(*s, s->dim, s->dir, 0);} // add_railing=0 (unused)
}
void building_t::try_connect_ext_basement_to_building(building_t &b) {
	assert(has_ext_basement() && b.has_ext_basement());
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), z_toler(0.1*get_trim_thickness());
	float const doorway_width(get_doorway_width()), wall_hwidth(0.8*doorway_width), min_shared_wall_len(2.01*(wall_hwidth + 2.0*wall_thickness));
	float const max_connect_dist(EXT_BASEMENT_JOIN_DIST*floor_spacing), min_connect_dist(2.1*doorway_width); // need enough space to fit two open doors
	assert(b.get_window_vspace() == floor_spacing);
	cube_t const &other_eb_bc(b.interior->basement_ext_bcube);
	assert(interior->basement_ext_bcube.z2() == other_eb_bc.z2()); // must be at same elevation
	assert((unsigned)  interior->ext_basement_hallway_room_id <   interior->rooms.size());
	assert((unsigned)b.interior->ext_basement_hallway_room_id < b.interior->rooms.size());
	ext_basement_room_params_t P, Pb, Padd; // P=input rooms for *this, Pb=input rooms for b, Padd=new rooms output for *this
	rand_gen_t rgen;
	rgen.set_state(interior->rooms.size(), b.interior->rooms.size());

	// find nearby candidate rooms
	for (auto r1 = interior->rooms.begin()+interior->ext_basement_hallway_room_id; r1 != interior->rooms.end(); ++r1) {
		cube_t search_area(*r1);
		search_area.expand_by(max_connect_dist);
		if (!search_area.intersects(other_eb_bc)) continue; // too far

		for (auto r2 = b.interior->rooms.begin()+b.interior->ext_basement_hallway_room_id; r2 != b.interior->rooms.end(); ++r2) {
			if (!search_area.intersects(*r2))        continue; // too far
			if (fabs(r1->z1() - r2->z1()) > z_toler) continue; // different floors/levels; do we need to check toler?
			
			if (r1->intersects(*r2)) { // previously failed at -1.12, -15.7
				cout << "Error: Invalid intersection of rooms at " << r1->str() << " and " << r2->str() << endl;
				continue; // uuuuh, just leave the rooms be, I guess...
			}
			for (unsigned d = 0; d < 2; ++d) { // r1/r2 join dim {x, y}
				float const shared_lo(max(r1->d[!d][0], r2->d[!d][0])), shared_hi(min(r1->d[!d][1], r2->d[!d][1]));
				if (shared_hi - shared_lo < min_shared_wall_len) continue; // check for projection in dim !d long enough to place a hallway and door
				// add an extra wall_thickness padding on either side to avoid wall trim clipping through a nearby room
				float const door_center(rgen.rand_uniform(shared_lo+wall_hwidth+wall_thickness, shared_hi-wall_hwidth-wall_thickness));
				bool const dir(r1->d[d][0] < r2->d[d][0]); // dir sign from r1 => r2 in dim d
				cube_t cand_join(*r1);
				cand_join.d[d][ dir] = r2->d[d][!dir];
				cand_join.d[d][!dir] = r1->d[d][ dir];
				if (cand_join.get_sz_dim(d) < min_connect_dist) continue;
				set_wall_width(cand_join, door_center, wall_hwidth, !d);
				assert(cand_join.is_strictly_normalized());
				cube_t test_cube(cand_join);
				test_cube.expand_in_dim(d, -wall_thickness); // shrink ends to avoid false intersection with rooms at either end
				populate_params_from_building(*  interior, P );
				populate_params_from_building(*b.interior, Pb);
				if (!  is_basement_room_placement_valid(test_cube, P,  d,  dir, nullptr, &b  )) continue; // add_end_door=nullptr
				if (!b.is_basement_room_placement_valid(test_cube, Pb, d, !dir, nullptr, this)) continue; // add_end_door=nullptr
				Padd.rooms.emplace_back(cand_join, 1, 0, d, dir); // is_hallway=1, has_stairs=0
				Padd.rooms.back().conn_bcube = *r2; // store room in the other building that we're connecting to in conn_bcube
			} // for d
		} // for r2
	} // for r1
	if (Padd.rooms.empty()) return; // failed to connect
	building_t *const buildings[2] = {this, &b};

	for (unsigned bix = 0; bix < 2; ++bix) {
		if (!buildings[bix]->interior->conn_info) {buildings[bix]->interior->conn_info.reset(new building_conn_info_t);}
	}
	for (auto const &r : Padd.rooms) { // add any new rooms from above
		// examples: (-0.9, -8.8), (2.1, -5.7), (2.47, -6.2), (8.96, -9.6), (5.02, -8.5), (-3.9, -8.13), (-1.94, -6.1)
		//if (fabs(r.x1()) < 10.0 && fabs(r.y1()) < 10.0) {cout << r.str() << endl;}
		unsigned const is_building_conn(r.hallway_dim ? 2 : 1);
		// skip one end in hallway_dim and make the other end (bordering the other building) thinner to avoid Z-fighting but still cast shadows
		interior->place_exterior_room(r, r, get_fc_thickness(), wall_thickness, P, basement_part_ix, 0, r.is_hallway, is_building_conn, r.hallway_dim, r.connect_dir);
		unsigned const conn_door_ix(b.interior->doors.size()); // index of door that will be added to the other building, and separates the two buildings
		// place doors at each end
		for (unsigned dir = 0; dir < 2; ++dir) {
			building_t *door_dest(buildings[bool(dir) ^ r.connect_dir ^ 1]); // add door to the building whose room it connects to
			cube_t const door(door_dest->add_ext_basement_door(r, doorway_width, r.hallway_dim, dir, 0, rgen)); // is_end_room=0
			// subtract door from walls of each building
			for (unsigned bix = 0; bix < 2; ++bix) {subtract_cube_from_cubes(door, buildings[bix]->interior->walls[r.hallway_dim]);}
		} // for dir
		b.interior->doors[conn_door_ix].is_bldg_conn = 1;
		cube_t ext_bcube(r);
		ext_bcube.d[r.hallway_dim][r.connect_dir] = r.conn_bcube.d[r.hallway_dim][r.connect_dir]; // extend to cover the entire width of the adjacent hallway in the other building

		for (unsigned bix = 0; bix < 2; ++bix) { // connect both ways
			// door belongs to b, which is the first building passed in
			buildings[bix]->interior->conn_info->add_connection(buildings[!bix], r, conn_door_ix, r.hallway_dim, r.connect_dir, (bix == 0));
			buildings[bix]->interior->basement_ext_bcube.union_with_cube(ext_bcube);
		}
	} // for r
	for (unsigned bix = 0; bix < 2; ++bix) {buildings[bix]->interior->remove_excess_capacity();} // optional optimization
}

void try_join_house_ext_basements(vect_building_t &buildings) {
	//timer_t timer("Join House Basements"); // ~10ms
	vector<vector<unsigned>> houses_by_city;

	for (auto b = buildings.begin(); b != buildings.end(); ++b) {
		if (!b->is_in_city || !b->is_house || !b->has_ext_basement()) continue; // Note: houses only, for now
		if (b->city_ix >= houses_by_city.size()) {houses_by_city.resize(b->city_ix+1);}
		houses_by_city[b->city_ix].push_back(b - buildings.begin());
	}
	for (vector<unsigned> const &work : houses_by_city) { // could be run in parallel, but not needed
		// do a quadratic iteration to find nearby houses in this city that can potentially be connected
		for (unsigned i = 0; i < work.size(); ++i) {
			building_t &b1(buildings[work[i]]);
			cube_t search_area(b1.interior->basement_ext_bcube);
			search_area.expand_by_xy(EXT_BASEMENT_JOIN_DIST*b1.get_window_vspace());

			for (unsigned j = i+1; j < work.size(); ++j) {
				building_t &b2(buildings[work[j]]);
				if (!search_area.intersects(b2.interior->basement_ext_bcube)) continue; // too far
				b1.try_connect_ext_basement_to_building(b2);
			}
		} // for i
	} // for work
}

building_t *building_t::get_conn_bldg_for_pt(point const &p, float radius) const {
	if (!player_in_basement || !has_conn_info()) return nullptr; // only active when player is in the basement
	return interior->conn_info->get_conn_bldg_for_pt(p, radius);
}
building_t *building_t::get_bldg_containing_pt(point const &p) {
	if (!player_in_basement) return nullptr; // only active when player is in the basement
	if (!has_conn_info()) return nullptr; // not really meant to be called in this case; caller must check for null ret and run some default logic
	return interior->conn_info->get_bldg_containing_pt(*this, p);
}
bool building_t::is_visible_through_conn(building_t const &b, vector3d const &xlate, float view_dist, bool expand_for_light) const {
	if (!player_in_basement) return 0; // only active when player is in the basement
	return (has_conn_info() && interior->conn_info->is_visible_through_conn(*this, b, xlate, view_dist, expand_for_light));
}
bool building_t::interior_visible_from_other_building_ext_basement(vector3d const &xlate, bool expand_for_light) const {
	if (!player_in_basement) return 0; // player not in basement; note that it's possible for the player to be in the basement and still see the conn room in the ext basement
	if (player_building == nullptr || player_building == this || !has_conn_info()) return 0;
	float const view_dist(8.0*get_window_vspace()); // arbitrary constant, should reflect length of largest hallway
	return player_building->is_visible_through_conn(*this, xlate, view_dist, expand_for_light);
}
cube_t building_t::get_conn_room_closest_to(point const &pos_bs) const { // in reference to the player's current building
	if (!player_in_basement || player_building == nullptr || player_building == this || !has_conn_info()) return cube_t();
	return interior->conn_info->get_conn_room_closest_to(*this, *player_building, pos_bs);
}

void building_conn_info_t::add_connection(building_t *b, cube_t const &room, unsigned door_ix, bool dim, bool dir, bool door_is_b) {
	if (conn.empty() || conn.back().b != b) {conn.emplace_back(b);} // register a new building if needed
	conn.back().rooms.emplace_back(room, door_ix, dim, dir, door_is_b);
}
building_t *building_conn_info_t::get_conn_bldg_for_pt(point const &p, float radius) const {
	for (conn_pt_t const &c : conn) {
		for (conn_room_t const &room : c.rooms) {
			if ((radius == 0.0) ? room.contains_pt(p) : sphere_cube_intersect(p, radius, room)) return c.b;
		}
	}
	return nullptr;
}
building_t *building_conn_info_t::get_bldg_containing_pt(building_t &parent, point const &p) const {
	for (conn_pt_t const &c : conn) {
		for (conn_room_t const &room : c.rooms) {
			if (room.contains_pt(p)) return (room.door_is_b ? &parent : c.b); // room belongs to one building and door belongs to the other
			cube_t other_side_of_door(room);
			other_side_of_door.d[room.dim][!room.dir] = room.d[room.dim][room.dir]; // flush with the door
			other_side_of_door.d[room.dim][ room.dir] = room.d[room.dim][room.dir] + (room.dir ? 1.0 : -1.0)*parent.get_doorway_width(); // extend into adj room
			if (other_side_of_door.contains_pt(p)) return (room.door_is_b ? c.b : &parent);
		}
	} // for c
	return nullptr;
}
bool building_conn_info_t::is_visible_through_conn(building_t const &parent, building_t const &target, vector3d const &xlate, float view_dist, bool expand_for_light) const {
	for (conn_pt_t const &c : conn) {
		if (c.b != &target) continue; // skip wrong building

		for (conn_room_t const &room : c.rooms) {
			cube_t room_cs(room + xlate);
			if (!room_cs.closest_dist_less_than(camera_pdu.pos, view_dist)) continue; // too far away
			if (expand_for_light) {room_cs.expand_by(view_dist);} // increase the bounds in case room is behind the player but light cast from it is visible
			if (!camera_pdu.cube_visible(room_cs)) continue;
			if (!expand_for_light) return 1; // can't ignore closed doors for room objects because they we may not draw the door itself
			// if this is a light, check if the connecting door is open
			door_t const &door((room.door_is_b ? target : parent).get_door(room.door_ix));
			if (door.open || door.open_amt > 0.0) return 1; // true if either about to open or not fully closed
		} // for room
	} // for c
	return 0;
}
door_t const *building_conn_info_t::get_door_to_conn_part(building_t const &parent, point const &pos_bs) const {
	for (conn_pt_t const &c : conn) {
		for (conn_room_t const &room : c.rooms) {
			if (room.contains_pt(pos_bs)) return &(room.door_is_b ? *c.b : parent).get_door(room.door_ix);
		}
	}
	return nullptr;
}
cube_t building_conn_info_t::get_conn_room_closest_to(building_t const &parent, building_t const &target, point const &pos_bs) const {
	cube_t closest;
	float dmin_sq(0.0);

	for (conn_pt_t const &c : conn) {
		if (c.b != &target) continue; // wrong building
		for (conn_room_t const &room : c.rooms) {
			if (room.contains_pt(pos_bs)) return room; // contained case returns immediately
			float const dist_sq(p2p_dist_sq(room.closest_pt(pos_bs), pos_bs));
			if (dmin_sq == 0.0 || dist_sq < dmin_sq) {dmin_sq = dist_sq; closest = room;}
		}
	}
	return closest;
}

