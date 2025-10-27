// 3D World - Building Extended Basements
// by Frank Gennari 03/11/2022

#include "function_registry.h"
#include "buildings.h"
#include <cfloat> // for FLT_MAX

extern int player_in_basement;
extern float DX_VAL_INV, DY_VAL_INV;
extern building_params_t global_building_params;
extern building_t const *player_building;

bool using_hmap_with_detail();
float get_ped_coll_radius();
bool cube_int_underground_obj(cube_t const &c);
unsigned choose_backrooms_wall_tex(rand_gen_t &rgen);


bool building_t::extend_underground_basement(rand_gen_t rgen) {
	if (!has_basement() || is_rotated() || !interior) return 0;
	if (is_prison() && !has_parking_garage)           return 0; // no extended basements in prison "dungeons"
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
				
				if (has_pg_ramp()) { // check for ramp coll
					cube_t test_cube(cand_door);
					test_cube.expand_by_xy(get_wall_thickness());
					if (interior->pg_ramp.intersects(test_cube)) continue;
				}
				float const fc_thick(get_fc_thickness());
				set_cube_zvals(cand_door, basement.z1()+fc_thick, basement.z2()-fc_thick); // change z to span floor to ceiling for interior door
				cand_door.translate_dim(dim, (dir ? 1.0 : -1.0)*0.25*get_wall_thickness()); // zero width, centered on the door
				bool ret(0);
				ret = add_underground_exterior_rooms(rgen, cand_door, basement, dim, dir, 0.25*len);
				if (ret) return 1; // exit on success
			} // for e
		} // for d
		if (!is_house) return 0; // not large enough for office building
	} // for len
	return 0;
}

float query_min_height(cube_t const &c, float stop_at) { // c_in is in global building space
	float hmin(FLT_MAX);

	if (using_tiled_terrain_hmap_tex() && !using_hmap_with_detail()) { // optimized flow when using heightmap texture; not adding xoff2/yoff2
		float const pad(1.0), step(0.5); // set pad to 1.0 rather than 0.5 to handle steep terrain
		float x1((c.x1() + X_SCENE_SIZE)*DX_VAL_INV + 0.5), x2((c.x2() + X_SCENE_SIZE)*DX_VAL_INV + 0.5);
		float y1((c.y1() + Y_SCENE_SIZE)*DY_VAL_INV + 0.5), y2((c.y2() + Y_SCENE_SIZE)*DY_VAL_INV + 0.5);

		for (float y = y1-pad; y < y2+pad; y += step) {
			for (float x = x1-pad; x < x2+pad; x += step) {
				min_eq(hmin, get_tiled_terrain_height_tex(x, y, 1)); // check every grid point with the X/Y range; nearest_texel=1
				if (hmin < stop_at) return hmin;
			}
		}
	}
	else { // we don't have the float heightmap here, so we have to do an expensive get_exact_zval() for each grid point
		float const x_step(0.5*DX_VAL), y_step(0.5*DY_VAL);
		cube_t c2(c);
		c2 += vector3d(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0); // cancel out xoff2/yoff2 translate

		for (float y = c2.y1()-y_step; y < c2.y2()+y_step; y += y_step) {
			for (float x = c2.x1()-x_step; x < c2.x2()+x_step; x += x_step) {
				min_eq(hmin, get_exact_zval(min(x, c2.x2()), min(y, c2.y2()))); // check every grid point with the X/Y range
				if (hmin < stop_at) return hmin;
			}
		}
	}
	return hmin;
}

struct ext_basement_room_params_t {
	vect_cube_t wall_exclude, wall_segs, temp_cubes, stairs_exp;
	vect_extb_room_t rooms;
	vector<stairs_place_t> stairs;
	
	void subtract_stairs_from_floor_ceil(cube_t const &c, vect_cube_t &out) {
		subtract_cubes_from_cube(c, stairs, out, temp_cubes, 2); // cut out stairs; zval_mode=2 (check for zval overlap)
	}
};

bool building_t::is_basement_room_not_int_bldg(cube_t const &room, building_t const *exclude, bool allow_outside_grid) const {
	// check for other buildings, including their extended basements
	if (check_buildings_cube_coll(room, 0, 1, this, exclude)) return 0; // xy_only=0, inc_basement=1, exclude ourself

	if (!allow_outside_grid) {
		cube_t const grid_bcube(get_grid_bcube_for_building(*this));
		assert(!grid_bcube.is_all_zeros()); // must be found
		assert(grid_bcube.contains_cube_xy(bcube)); // must contain our building
		if (!grid_bcube.contains_cube_xy(room)) return 0; // outside the grid (tile or city) bcube
	}
	if (cube_int_underground_obj(room)) return 0; // check tunnels, in-ground pools, etc.
	return 1;
}
bool building_t::is_basement_room_under_mesh_not_int_bldg(cube_t const &room, building_t const *exclude, bool allow_outside_grid) const {
	float const ceiling_zval(room.z2() - get_fc_thickness());
	if (query_min_height(room, ceiling_zval) < ceiling_zval) return 0; // check for terrain clipping through ceiling
	return is_basement_room_not_int_bldg(room, exclude, allow_outside_grid);
}
bool building_t::is_basement_room_placement_valid(cube_t &room, ext_basement_room_params_t &P, bool dim, bool dir, bool *add_end_door, building_t const *exclude) const {
	float const wall_thickness(get_wall_thickness()), wall_expand_toler(0.1*wall_thickness);
	cube_t test_cube(room);
	test_cube.expand_in_z(-0.01*test_cube.dz()); // shrink slightly so that rooms on different floors can cross over each other
	
	if (!P.rooms.empty()) { // not the first hallway; check if too close to the basement such that the wall or trim will clip through the basement wall
		cube_t room_exp(test_cube);
		room_exp.expand_by_xy(wall_thickness + get_trim_thickness());
		if (room_exp.intersects(P.rooms.front())) return 0;
	}
	test_cube.d[dim][!dir] -= (dir ? -1.0 : 1.0)*wall_expand_toler; // shrink slightly to avoid intersections with our parent room; or pass in the parent room?
	test_cube.d[dim][ dir] -= (dir ? -1.0 : 1.0)*wall_expand_toler; // expand the end slightly
	test_cube.expand_in_dim(!dim, wall_expand_toler); // expand slightly on the sides to avoid adjacent rooms
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
	// Note: misnamed: hallway for houses, but may be backrooms for offices with parking garages
	extb_room_t hallway(basement, 0); // is_hallway=0; will likely be set below
	set_wall_width(hallway, door_bcube.get_center_dim(!wall_dim), 0.5*hallway_width, !wall_dim);
	hallway.d[wall_dim][!wall_dir] = ext_wall_pos; // flush with the exterior wall/door
	hallway.d[wall_dim][ wall_dir] = ext_wall_pos + (wall_dir ? 1.0 : -1.0)*hallway_len;
	assert(hallway.is_strictly_normalized());
	ext_basement_room_params_t P;
	if (!is_basement_room_placement_valid(hallway, P, wall_dim, wall_dir)) return 0; // try to place the hallway; add_end_door=nullptr
	// valid placement; now add the door, hallway, and connected rooms
	has_basement_door       = 1;
	interior->extb_wall_dim = wall_dim;
	interior->extb_wall_dir = wall_dir;
	// Note: recording the door_stack index rather than the door index allows us to get either the first door or the first stack
	interior->ext_basement_door_stack_ix = interior->door_stacks.size();
	float const fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	P.wall_exclude.push_back(basement);
	P.wall_exclude.back().expand_in_dim(wall_dim, 1.1*get_trim_thickness()); // add slightly expanded basement to keep interior wall trim from intersecting exterior walls
	P.wall_exclude.push_back(door_bcube);
	P.wall_exclude.back().expand_in_dim(wall_dim, 2.0*wall_thickness); // make sure the doorway covers the entire wall thickness
	vector<room_t> &rooms(interior->rooms);
	interior->ext_basement_hallway_room_id = rooms.size();
	door_t Door(door_bcube, wall_dim, wall_dir, rgen.rand_bool());
	add_interior_door(Door, 0, 1); // open 50% of the time; is_bathroom=0, make_unlocked=1
	P.rooms.emplace_back(basement, 0);
	P.rooms.push_back(hallway);
	hallway.conn_bcube = basement; // make sure the basement is included
	bool added_lg_room(0);

	if (!is_house && has_parking_garage) { // office building with parking garage
		bool const add_mall(global_building_params.max_mall_levels > 0 && global_building_params.mall_prob > 0.0 && (!global_building_params.no_retail_and_mall || !has_retail()));
		bool const try_mall_first(add_mall && rgen.rand_probability(global_building_params.mall_prob));

		for (unsigned n = 0; n < 2; ++n) {
			bool const is_mall((n == 0) == try_mall_first);
			if (is_mall && !add_mall)  continue;
			unsigned const num_floors_added(max_expand_underground_room(hallway, wall_dim, wall_dir, is_mall, rgen));
			if (num_floors_added == 0) continue;
			interior->num_extb_floors = num_floors_added;
			hallway.is_hallway = 0; // should already be set to 0, but this makes it more clear
			// currently, the extended basement can only be a network of connected hallways with leaf rooms, or a single large basement room (this case);
			// if we want to allow both (either a large room connected to a hallway or a large room with hallways coming off of it), we need per-room flags

			if (is_mall) {
				interior->mall_info.reset(new building_mall_info_t);
				water_damage = crack_damage = 0.0; // no damage for malls
				door_t &ent_door(interior->doors.back());
				door_stack_t &ent_ds(interior->door_stacks.back());
				ent_door.open_dir ^= 1; // door opens into the parking garage rather than the mall
				ent_ds  .open_dir ^= 1;
				
				if (ent_door.z2() < (hallway.z2() - 2.0*fc_thick)) { // counts as multi-floor (for drawing top edge)
					ent_door.set_mult_floor();
					ent_ds  .set_mult_floor();
				}
				setup_mall_concourse(hallway, wall_dim, wall_dir, rgen);
			}
			else { // backrooms, possibly flooded
				unsigned const num_floors(setup_multi_floor_room(hallway, Door, wall_dim, wall_dir, rgen));
				interior->has_backrooms = 1;
				interior->backrooms_tid = choose_backrooms_wall_tex(rgen); // may be 0/none

				if (num_floors > 1) { // lowest level of multilevel rooms has water; no water if there's a single floor
					float wmin(global_building_params.basement_water_level_min), wmax(global_building_params.basement_water_level_max);
					if (wmax < wmin) {swap(wmin, wmax);} // user specfied the values backwards? this isn't error checked in the option parsing, so swap the values

					if (wmax > 0.0) {
						float const ftv(get_floor_thick_val());
						float water_level((wmin == wmax) ? wmin : rgen.rand_uniform(wmin, wmax)); // can be a single value or a range
						min_eq(water_level, (num_floors - 1.0f)); // top floor can't have water
						// handle water near the level of an upper floor; offset to prevent Z-fighting
						if (water_level > 0.5 && fract(water_level + 0.5*ftv) < 0.6*ftv) {water_level -= 0.6*ftv;}
						if (water_level > 0.0) {interior->water_zval = hallway.z1() + fc_thick + water_level*get_window_vspace();}
					}
				}
			}
			added_lg_room = 1;
			break; // done/success
		} // for n
	}
	if (!added_lg_room) { // not a larger underground room; recursively add rooms connected to this hallway in alternating dimensions
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
	tid_nm_pair_t ceil_tex;
	if (!is_house && hallway.is_hallway) {get_ceil_tex_and_color(hallway, ceil_tex);}
	bool const remove_ceil_tiles(ceil_tex.tid == get_material().ceil_tex.tid);
	if (remove_ceil_tiles) {remove_ceiling_tiles(hallway, ceil_tex, P, rgen);}
	interior->place_exterior_room(hallway, wall_area, fc_thick, wall_thickness, P, basement_part_ix, 0, hallway.is_hallway); // use basement part_ix; num_lights=0
	if (interior->has_backrooms) {rooms.back().assign_all_to(RTYPE_BACKROOMS);} // make it backrooms
	else if (has_mall())         {rooms.back().assign_all_to(RTYPE_MALL     );} // make it a mall concourse
	if (has_mall()) {rooms.back().is_single_floor = 1;}
	reserve_extra(rooms, ((P.rooms.size()-2) + 1)); // allocate an extra room for a possible connector to an adjacent building

	for (auto r = P.rooms.begin()+2; r != P.rooms.end(); ++r) { // skip basement and primary hallway
		if (remove_ceil_tiles && r->is_hallway) {remove_ceiling_tiles(*r, ceil_tex, P, rgen);}
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
	if (!added_lg_room) {maybe_assign_extb_room_as_swimming(rgen);}
	if (has_mall()) {add_mall_stores(hallway, wall_dim, !wall_dir, rgen);}

	if (!has_backrooms_or_mall() && global_building_params.add_basement_tunnels) { // maybe add tunnel connections to hallways
		for (unsigned r = interior->ext_basement_hallway_room_id; r < rooms.size(); ++r) {try_place_tunnel_at_extb_hallway_end(rooms[r], r, rgen);}
	}
	return 1;
}

// num_rooms=5255 num_tiles=10453 num_verts=416764
void building_t::remove_ceiling_tiles(cube_t const &room_in, tid_nm_pair_t const &ceil_tex, ext_basement_room_params_t &P, rand_gen_t &rgen) {
	assert(interior);
	unsigned const num_missing(rgen.rand() % 4); // 0-3
	if (num_missing == 0) return;
	float const wall_thick(get_wall_thickness()), fc_thick(get_fc_thickness()), tile_thick(0.2*fc_thick);
	cube_t ceiling(room_in), room;
	ceiling.z1() = room_in.z2() - fc_thick;
	P.subtract_stairs_from_floor_ceil(ceiling, P.wall_segs);

	for (cube_t const &c : P.wall_segs) { // find the ceiling part with the largest area
		if (c.get_area_xy() > 0.5*room_in.get_area_xy()) {room = c; break;}
	}
	if (room.is_all_zeros()) return; // mostly stairs?
	float tscale[2] = {ceil_tex.get_drawn_tscale_x(), 2.0f*ceil_tex.get_drawn_tscale_y()}; // width is half the texture
	vector2d const room_sz(room.get_size_xy());
	bool const dim(room_sz.x < room_sz.y); // long dim
	if (dim) {swap(tscale[0], tscale[1]);}
	for (unsigned d = 0; d < 2; ++d) {tscale[d] = max(1, round_fp(tscale[d]*room_sz[d]))/room_sz[d];} // exact tiling
	vector2d const tile_sz(1.0/tscale[0], 1.0/tscale[1]);
	unsigned const num[2] = {unsigned(room_sz.x/tile_sz.x), unsigned(room_sz.y/tile_sz.y)};
	if (num[0] == 0 || num[1] == 0) return; // room narrower than one tile
	// add ceiling space
	unsigned const room_ix(interior->rooms.size()); // will be the next room added
	ceiling_space_t ceil_space(room, room_ix, num[0]+1, num[1]+1, tile_sz); // verts is tiles+1
	ceil_space.expand_by_xy(-0.52*wall_thick); // slightly more than half to avoid Z-fighting with pool room walls
	ceil_space.z1() = room.z2() - fc_thick;
	ceil_space.z2() += 0.01*fc_thick; // expand up slightly above the reach of the upward pointing ceiling lights, since they're not shadowed properly
	cube_t cs_ext(ceil_space);
	cs_ext.z2() += fc_thick;
	if (is_basement_room_under_mesh_not_int_bldg(cs_ext)) {ceil_space.z2() = cs_ext.z2();} // expand upward if there's space
	// add missing tiles
	unsigned const missing_tiles_start(interior->missing_ceil_tiles.size());
	cube_t cut, cut_clip_cube(ceil_space);
	cut_clip_cube.expand_by_xy(-0.02*wall_thick); // slight shrink to make sure it's in front of the inside wall
	set_cube_zvals(cut, ceiling.z1(), (ceiling.z1() + tile_thick)); // make it relatively thin

	for (unsigned n = 0; n < num_missing; ++n) {
		for (unsigned N = 0; N < 10; ++N) { // N tries
			unsigned ix[2]={};

			for (unsigned d = 0; d < 2; ++d) { // xy
				ix[d] = rgen.rand() % num[d];
				float const pos(room.d[d][0] + ix[d]*tile_sz[d]);
				cut.d[d][0] = pos + ((d ^ dim) ? 0.085 : 0.042)*tile_sz[d]; // clip off the lower edge to preserve the brown frame
				cut.d[d][1] = pos + tile_sz[d];
			}
			assert(cut.intersects_xy(cut_clip_cube));
			cut.intersect_with_cube_xy(cut_clip_cube);
			bool dup(0);
			for (auto i = interior->missing_ceil_tiles.begin()+missing_tiles_start; i != interior->missing_ceil_tiles.end(); ++i) {dup |= (cut == *i);}
			if (dup) continue; // we already removed this tile
			interior->missing_ceil_tiles.emplace_back(cut, room_ix); // record so that tiles can be added on the floor later
			
			for (unsigned x = ix[0]; x <= ix[0]+1; ++x) { // add light for this opening
				for (unsigned y = ix[1]; y <= ix[1]+1; ++y) {ceil_space.set_light_val(x, y, 255);}
			}
			break;
		} // for N
	} // for n
	interior->ceiling_spaces.push_back(ceil_space);
}

void building_t::add_ceiling_tile_objects(rand_gen_t &rgen) {
	assert(has_room_geom());
	float const fc_thick(get_fc_thickness()), light_amt(1.0);
	unsigned const pipe_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING), wire_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_IN_HALLWAY);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_end(objs.size());
	unsigned cur_obj_ix(0); // used for lights iteration; lights should be in the same order as ceiling_spaces
	auto tile_iter(interior->missing_ceil_tiles.begin());
	vect_cube_t miss_tiles, pipe_avoid, light_bcs, tile_block;

	for (ceiling_space_t const &cs : interior->ceiling_spaces) {
		miss_tiles.clear();
		pipe_avoid.clear();
		light_bcs .clear();
		tile_block.clear();

		for (; cur_obj_ix != objs_end && objs[cur_obj_ix].room_id <= cs.room_ix; ++cur_obj_ix) {
			room_object_t const &light(objs[cur_obj_ix]);
			if (light.room_id != cs.room_ix || light.type != TYPE_LIGHT) continue; // not a light in this room
			light_bcs.push_back(light);
			cube_t base(light);
			base.z1()  = light.z2(); // flush with top of light
			base.z2() += 1.0*light.dz();
			objs.emplace_back(base, TYPE_METAL_BAR, cs.room_ix, light.dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, LT_GRAY, EF_Z1); // skip bottom
			// add two support posts
			float const post_radius(0.04*light.get_width());
			cube_t post;
			set_cube_zvals(post, base.z2(), cs.z2());
			set_wall_width(post, light.get_center_dim(!light.dim), post_radius, !light.dim);
			if (pipe_avoid.empty()) {pipe_avoid.push_back(post);} // only need to add the first light, since they should be in a line

			for (unsigned d = 0; d < 2; ++d) {
				set_wall_width(post, (light.d[light.dim][d] + (d ? -1.0 : 1.0)*4.0*post_radius), post_radius, light.dim);
				objs.emplace_back(post, TYPE_METAL_BAR, cs.room_ix, 0, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, LT_GRAY, EF_Z12); // skip top and bottom
			}
		} // for cur_obj
		// place ceiling tiles on the floor
		room_t const &room(get_room(cs.room_ix));
		cube_t place_area(get_walkable_room_bounds(room));
		place_area.expand_by_xy(-get_trim_thickness());
		bool const dim(cs.dx() < cs.dy()); // long dim
		float const floor_zval(room.z1() + fc_thick), fc_gap(get_floor_ceil_gap());

		for (; tile_iter != interior->missing_ceil_tiles.end() && tile_iter->ix <= cs.room_ix; ++tile_iter) {
			assert(tile_iter->ix == cs.room_ix); // can't have tile without a ceiling space
			// place tile on the floor
			miss_tiles.push_back(*tile_iter); // needed for visible pipe placement
			cube_t const tile(*tile_iter);

			for (unsigned n = 0; n < 40; ++n) { // 40 placement attempts
				cube_t ftile(*tile_iter); // tile on floor
				ftile.translate_dim(2, (floor_zval - tile.z1()));
				for (unsigned d = 0; d < 2; ++d) {ftile.translate_dim(d, 2.0*rgen.signed_rand_float()*tile.get_sz_dim(d));} // move +/- 2 tiles away in both dims
				float const rot_angle(PI*rgen.rand_float()*(1.0 - n/40.0)); // reduce angle with later iterations
				point const rot_pt(ftile.get_cube_center());
				cube_t tile_rot(ftile - rot_pt); // rotate about the tile center
				tile_rot = rotate_cube(tile_rot, plus_z, rot_angle) + rot_pt;
				if (!place_area.contains_cube_xy(tile_rot))        continue; // clipping through walls
				if (has_bcube_int(tile_rot, tile_block))           continue; // check previously placed tiles
				if (has_bcube_int(tile_rot, interior->stairwells)) continue; // don't place on stairs
				if (interior->is_cube_close_to_doorway(tile_rot, room, 0.0, 1, 1)) continue; // no clipping through doors
				objs.emplace_back(ftile, TYPE_CEIL_TILE, cs.room_ix, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, colorRGBA(WHITE, rot_angle)); // store rot angle in alpha
				tile_block.push_back(tile_rot);
				break; // success
			} // for n
			if (rgen.rand_float() < 0.5) { // add hanging wire
				bool const cx(rgen.rand_bool()), cy(rgen.rand_bool()); // choose a random corner
				float const radius(rgen.rand_uniform(0.04, 0.08)*fc_thick);
				point const corner((tile.d[0][cx] + (cx ? -1.0 : 1.0)*radius), (tile.d[1][cy] + (cy ? -1.0 : 1.0)*radius), cs.z2()); // starts at ceiling
				cube_t wire(corner, corner);
				wire.expand_by_xy(radius);

				if (!has_bcube_int_xy(wire, light_bcs)) {
					wire.z1() = cs.z1() - rgen.rand_uniform(0.1, 0.3)*fc_gap; // set length
					objs.emplace_back(wire, TYPE_ELEC_WIRE, cs.room_ix, 0, 1, wire_flags, light_amt, SHAPE_CYLIN, BLACK); // vertical
				}
			}
			if (rgen.rand_float() < 0.5) { // add spider web
				bool const tdim(rgen.rand_bool()), tdir(rgen.rand_bool()); // choose a random edge
				cube_t web(tile); // copy tile width
				set_cube_zvals(web, tile.z2(), cs.z2()); // top of tile to top of ceiling
				set_wall_width(web, tile.d[tdim][tdir], 0.02*fc_thick, tdim); // close to zero area
				float const web_hwidth(0.6*web.dz());

				if (web.get_sz_dim(!tdim) > 1.2*web_hwidth) { // too long, shrink to random segment
					set_wall_width(web, rgen.rand_uniform((tile.d[!tdim][0] + web_hwidth), (tile.d[!tdim][1] - web_hwidth)), web_hwidth, !tdim);
				}
				objs.emplace_back(web, TYPE_SPIWEB, cs.room_ix, tdim, tdir, RO_FLAG_NOCOLL, light_amt);
			}
		} // for tile_iter
		// add pipes in the ceiling
		//float const pipe_z1(cs.z1() + 0.2*fc_thick); // resting on ceiling tile
		float const pipe_zc(cs.zc()); // center of space
		unsigned const num_pipes(1 + (rgen.rand() % 5)); // 1-5
		cube_t pipe(cs); // full length of ceiling space
		unsigned const NCOLORS = 7;
		colorRGBA const colors[NCOLORS] = {DARK_BRASS_C, COPPER_C, COPPER_C, WHITE, GRAY, LT_GRAY, DK_GRAY};

		for (unsigned n = 0; n < num_pipes; ++n) {
			float const radius(rgen.rand_uniform(0.12, 0.25)*fc_thick), edge_spacing(2.0*radius);
			colorRGBA const &pipe_color(colors[rgen.rand() % NCOLORS]);

			for (unsigned N = 0; N < 10; ++N) { // 10 attempts to place a pipe
				float const pipe_pos(rgen.rand_uniform((cs.d[!dim][0] + edge_spacing), (cs.d[!dim][1] - edge_spacing)));
				set_wall_width(pipe, pipe_zc,  radius, 2);
				set_wall_width(pipe, pipe_pos, radius, !dim);
				if (!has_bcube_int_xy(pipe, miss_tiles)) break; // not visible through a missing tile, skip (but counts as a pipe)
				if ( has_bcube_int_xy(pipe, pipe_avoid)) continue; // blocked by light post or another pipe
				interior->room_geom->objs.emplace_back(pipe, TYPE_PIPE, cs.room_ix, dim, 0, pipe_flags, light_amt, SHAPE_CYLIN, pipe_color); // horizontal
				pipe_avoid.push_back(pipe);
				// TODO: hanging brackets over miss_tiles?
				break; // success
			} // for N
		} // for n
	} // for cs
}

void extend_adj_cubes(cube_t const &oldc, cube_t const &newc, vect_cube_t &cubes, float wall_thickness, unsigned wall_dim=2) {
	cube_t query(oldc);
	query.expand_by(0.5*wall_thickness);

	for (cube_t &c : cubes) {
		if (!query.intersects(c)) continue;

		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				float const old_pos(oldc.d[dim][dir]), new_pos(newc.d[dim][dir]);
				if (new_pos == old_pos) continue; // no edge movement
				if (fabs(c.d[dim][dir] - old_pos) > wall_thickness) continue; // not along this edge
				if (dim == wall_dim) {c.translate_dim(dim, (new_pos - old_pos));} // translate the wall
				else {c.d[dim][dir] = new_pos;} // translate edge by the same amount
			}
		}
		float const z2_ext(newc.z2() - oldc.z2());
		if (z2_ext == 0.0) {} // no z extend
		else if (wall_dim < 2) { // wall
			if (c.z2() <= oldc.z2()) {c.z2() += z2_ext;} // extend wall upward (not wall above the door)
		}
		else if (c.zc() > oldc.zc() && c.contains_pt_xy(oldc.get_cube_center())) {c.translate_dim(2, z2_ext);} // extend ceiling upward
	} // for c
}
bool building_t::check_pool_room_slice_valid(cube_t const &slice, int skip_room_ix) const {
	if (slice.intersects(get_basement())) return 0;

	for (auto r = interior->ext_basement_rooms_start(); r != interior->rooms.end(); ++r) {
		if (int(r - interior->rooms.begin()) != skip_room_ix && r->intersects_no_adj(slice)) return 0;
	}
	return is_basement_room_under_mesh_not_int_bldg(slice);
}
void building_t::maybe_assign_extb_room_as_swimming(rand_gen_t &rgen) {
	// swimming pools are in the basement so that we don't need to cut out the terrain, and in the extended basement so that we can create a custom lower floor area
	assert(has_ext_basement());
	vector<room_t> &rooms(interior->rooms);
	int const start_room_ix(interior->ext_basement_hallway_room_id);
	assert(start_room_ix >= 0 && unsigned(start_room_ix) < rooms.size());
	float const floor_spacing(get_window_vspace()), min_dmin(3.0*floor_spacing), pool_max_depth(1.0*floor_spacing), sea_level(get_max_sea_level());
	auto const first_room(rooms.begin() + start_room_ix + 1); // skip ext basement connector hallway
	int largest_valid_room(-1);
	float best_dmin(min_dmin);

	for (auto r = first_room; r != rooms.end(); ++r) {
		if (r->is_hallway || r->has_stairs)       continue;
		if (r->z1() - pool_max_depth < sea_level) continue; // too deep
		float const dmin(min(r->dx(), r->dy()));
		if (dmin < best_dmin) continue; // too small, or smaller than best room
		bool invalid(0);

		for (auto r2 = first_room; r2 != rooms.end(); ++r2) { // check for any rooms below this one
			if (r2 != r && r2->z1() < r->z1() && r2->intersects_xy(*r)) {invalid = 1; break;}
		}
		if (invalid) continue;
		largest_valid_room = (r - rooms.begin());
		best_dmin = dmin;
	} // for r
	if (largest_valid_room < 0) return; // no valid room found
	room_t &room(rooms[largest_valid_room]);
	// make sure there's space below the room for the pool
	cube_t room_ext_down(room);
	room_ext_down.z1() -= pool_max_depth; // account for the bottom of the pool
	if (!is_basement_room_not_int_bldg(room_ext_down)) return; // only fails 1-2% of the time
	// choose the first door as the main entrance; there's likely only one
	vect_door_stack_t &doorways(get_doorways_for_room(room, room.z1()));
	if (doorways.empty()) {std::cerr << "Error: No doorway found for pool in room " << room.str() << endl;} // can this ever fail?

	// attempt to expand the room so that we can fit a larger pool; assume connected to a hallway with <door>, and expand away from door
	float const wall_thickness(get_wall_thickness()), wall_pad(2.0*wall_thickness + get_trim_thickness()), exp_step(0.25*floor_spacing), fc_thickness(get_fc_thickness());
	room_t const orig_room(room);
	bool long_dim(room.dx() < room.dy()); // likely door.dim, unless there are multiple doors
	bool const exp_dim(doorways.empty() ? long_dim : doorways.front().dim);
	bool const exp_dir(doorways.empty() ? rgen.rand_bool() : (doorways.front().get_center_dim(exp_dim) < room.get_center_dim(exp_dim)));
	cube_t const &basement(get_basement());
	// try to expand away from the first door, then expand to either side of the door, starting with a random side
	bool const pref_exp_other_dir(rgen.rand_bool());
	bool const exp_dims[3] = {exp_dim, !exp_dim, !exp_dim}, exp_dirs[3] = {exp_dir, pref_exp_other_dir, !pref_exp_other_dir};

	for (unsigned d = 0; d < 3; ++d) { // {away from door, to one side, to the other side}
		bool const edim(exp_dims[d]), edir(exp_dirs[d]);
		float &end_wall(room.d[edim][edir]);
		float const step_dist((edir ? 1.0 : -1.0)*exp_step), exp_max(((d == 0) ? 2.0 : 1.0)*room.get_sz_dim(edim)); // 2-3x longer
		unsigned const num_steps(round_fp(exp_max/exp_step));

		for (unsigned n = 0; n < num_steps; ++n) {
			cube_t slice(room);
			slice.d[edim][!edir]  = end_wall;
			slice.d[edim][ edir] += step_dist;
			slice.expand_by_xy(wall_pad); // add some space around it for the walls
			slice.z1() -= pool_max_depth; // account for the bottom of the pool
			if (!check_pool_room_slice_valid(slice, largest_valid_room)) break;
			end_wall += step_dist;
		} // for n
	} // for d
	if (room.z2() < basement.z2()) { // try to expand upward for a higher ceiling if on a lower level (and ceiling doesn't go above ground_floor_z1)
		unsigned const z_exp_num = 10;
		float const z_exp_step(0.1*floor_spacing);

		for (unsigned n = 0; n < z_exp_num; ++n) {
			cube_t slice(room);
			set_cube_zvals(slice, room.z2(), (room.z2() + z_exp_step));
			if (!check_pool_room_slice_valid(slice, largest_valid_room)) break;
			room.z2() = slice.z2(); // extend upward
			if (room.z2() > ground_floor_z1) {room.z2() = ground_floor_z1; break;} // stop and clamp if too tall
		}
		if (room.z2() > orig_room.z2()) { // was extended vertically; add missing wall sections above doors
			cube_t room_exp(room);
			room_exp.expand_by_xy(wall_thickness);

			for (door_stack_t &ds : interior->door_stacks) {
				if (!ds.intersects(room_exp)) continue; // door not connected to this room
				add_wall_section_above_pool_room_door(ds, room);
			}
		}
	}
	assert(room.is_strictly_normalized());
	long_dim = (room.dx() < room.dy()); // recalculate, in case the aspect ratio of the room changed when expanding
	bool const was_extended(room != orig_room);
	room.is_single_floor = 1; // even if it was extended upward

	if (was_extended) { // room was extended; move or extend any connected walls, ceilings, and floors
		extend_adj_cubes(orig_room, room, interior->floors,   wall_thickness);
		extend_adj_cubes(orig_room, room, interior->ceilings, wall_thickness);
		for (unsigned d = 0; d < 2; ++d) {extend_adj_cubes(orig_room, room, interior->walls[d], wall_thickness, d);}
		interior->basement_ext_bcube.union_with_cube(room);
	}
	float const doorway_width(get_doorway_width()), min_spacing(1.5*doorway_width), floor_zval(room.z1() + fc_thickness);
	indoor_pool_t &pool(interior->pool);
	pool.copy_from(get_walkable_room_bounds(room));
	pool.expand_by_xy(-min_spacing);
	assert(pool.is_strictly_normalized());
	float const pool_length(pool.get_sz_dim(long_dim)), pool_depth(min(pool_max_depth, 0.3f*pool_length));
	set_cube_zvals(pool, (floor_zval - pool_depth), floor_zval);
	pool.orig_z1    = pool.z1();
	pool.bottomless = ((pool.z1() > sea_level + 4.0*floor_spacing) && rgen.rand_float() < 0.2); // bottomless 20% of the time, if high enough above sea level

	if (pool.bottomless) {
		pool.z1() = max((pool.z1() - 10.0f*pool_depth), (sea_level + fc_thickness)); // not quite bottomless, but very deep
		
		if (!is_basement_room_not_int_bldg(pool)) { // some other extended basement is below the pool, can't make bottomless
			pool.z1() = pool.orig_z1;
			pool.bottomless = 0;
		}
	}
	float const pool_width(pool.get_sz_dim(!long_dim)), extra_width(pool_width - 2.0*floor_spacing), pool_bottom(pool.z1() - fc_thickness);
	if (pool_width < floor_spacing) return; // too small; shouldn't happen unless the door width to floor spacing values are wrong
	if (extra_width > 0.0) {pool.expand_in_dim(!long_dim, -0.25*extra_width);} // can make narrower if there's extra space
	pool.dim = long_dim; // or door.dim?
	bool const dir(doorways.empty() ? rgen.rand_bool() : (room.get_center_dim(pool.dim) < doorways.front().get_center_dim(pool.dim)));
	float const door_shift((dir ? -1.0 : 1.0)*0.5*doorway_width);
	if (was_extended) {pool.d[pool.dim][dir] += door_shift;} // shift edge away from door
	else {pool.translate_dim(pool.dim, door_shift);} // translate away from door
	if (pool_length > 2.0*floor_spacing) {pool.d[pool.dim][!dir] -= door_shift;} // shrink pool on end opposite the door if long to make room for a diving board
	pool.dir     = dir;
	pool.room_ix = largest_valid_room;
	pool.valid   = 1;
	pool.shallow_zval = pool.orig_z1; // default is all deep
	assert(pool.is_strictly_normalized());
	
	// cut out a space in the floor for the pool
	for (cube_t &f : interior->floors) {
		if (!room.contains_cube(f)) continue;
		vect_cube_t floor_parts;
		subtract_cube_from_cube(f, pool, floor_parts); // should only get here for one floor
		assert(floor_parts.size() == 4);
		f = pool;
		set_cube_zvals(f, pool_bottom, pool.z1()); // move to the bottom of the pool
		vector_add_to(floor_parts, interior->floors);
	} // for f
	room.assign_to(RTYPE_SWIM, 0);
	if (pool.bottomless || rgen.rand_float() < 0.8) {interior->water_zval = pool.z2() - 0.05*pool_depth;} // add water to the pool 80% of the time, always if bottomless
	min_eq(interior->basement_ext_bcube.z1(), pool_bottom); // is this a good idea? it certainly makes other logic easier
}

void building_t::add_wall_section_above_pool_room_door(door_stack_t &ds, room_t const &room) {
	float const ceil_zval(room.z2() - get_fc_thickness());
	if (ds.z2() >= ceil_zval) return; // no gap above door
	ds.set_mult_floor(); // counts as multi-floor (for drawing top edge)
	interior->get_door(ds.first_door_ix).set_mult_floor();
	cube_t wall(ds);
	set_wall_width(wall, ds.get_center_dim(ds.dim), 0.5*get_wall_thickness(), ds.dim);
	set_cube_zvals(wall, ds.z2(), ceil_zval);
	interior->walls[ds.dim].push_back(wall);
}

bool building_t::add_ext_basement_rooms_recur(extb_room_t &parent_room, ext_basement_room_params_t &P, float door_width, bool dim, unsigned depth, rand_gen_t &rgen) {
	// add doors and other rooms along hallway; currently, all rooms are hallways
	float const parent_len(parent_room.get_sz_dim(!dim)), parent_width(parent_room.get_sz_dim(dim));
	float const min_length(max(4.0f*door_width, 0.5f*parent_len)), max_length(max(parent_len, 2.0f*min_length));
	// add at least a doorway's worth of spacing to the connecting hallway so that the doors don't block each other when open
	float const pos_lo(parent_room.d[!dim][0] + (parent_room.connect_dir ? 1.5 : 0.7)*door_width);
	float const pos_hi(parent_room.d[!dim][1] - (parent_room.connect_dir ? 0.7 : 1.5)*door_width);
	if (pos_lo >= pos_hi) return 0; // not enough space to add a door
	bool const is_end_room(depth >= global_building_params.max_ext_basement_room_depth);
	float const min_width_scale(is_end_room ? 1.0 : 0.9), max_width_scale(is_end_room ? 3.0 : 1.5);
	bool was_added(0);

	for (unsigned n = 0; n < global_building_params.max_ext_basement_hall_branches; ++n) {
		if (interior->rooms.size() + P.rooms.size() >= 255) break; // cap the number of rooms at 255 so that we can store room_ix in a uint8_t

		for (unsigned N = 0; N < 2; ++N) { // make up to 2 tries to place this room
			bool const dir(rgen.rand_bool());
			float const conn_edge(parent_room.d[dim][dir]), room_pos(rgen.rand_uniform(pos_lo, pos_hi));
			float const room_length(rgen.rand_uniform(min_length, max_length));
			float const room_width(max(1.5f*door_width, rgen.rand_uniform(min_width_scale, max_width_scale)*parent_width));
			extb_room_t room((cube_t)parent_room); // copy bcube but not flags; sets correct zvals
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

// is_mall=0: backrooms, is_mall=1: mall concourse
unsigned building_t::max_expand_underground_room(cube_t &room, bool dim, bool dir, bool is_mall, rand_gen_t &rgen) const {
	float const floor_spacing(get_window_vspace()), step_len(1.0*floor_spacing), room_len(room.get_sz_dim(dim));
	float const room_width_min(is_mall ? 4.0*floor_spacing : 0.5*room_len), room_width_max(is_mall ? max(8.0*floor_spacing, 0.25*room_len) : 2.0*room_len);
	float const room_len_min((is_mall ? 2.0 : 1.0)*room_len), room_len_max((is_mall ? 4.0 : 2.0)*room_len);
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
			if (exp_slice.get_sz_dim(edim) > ((edim == dim) ? room_len_max : room_width_max)) {cant_expand[d] = 1; continue;} // too large
			exp_slice.d[edim][!edir]  = exp_room.d[edim][edir]; // shrink to zero area since we've already checked exp_room
			assert(exp_slice.is_strictly_normalized());
			if (!is_basement_room_under_mesh_not_int_bldg(exp_slice)) {cant_expand[d] = 1; continue;} // can't expand this edge any more
			exp_room.d[edim][edir] = exp_slice.d[edim][edir]; // keep edge movement
			any_valid = 1;
		} // for i
		if (!any_valid) break;
	} // end while
	assert(exp_room.contains_cube(room));
	if (exp_room.get_sz_dim(!dim) < room_width_min) return 0; // room is too narrow, make it a hallway or mall concourse instead
	if (exp_room.get_sz_dim( dim) < room_len_min  ) return 0; // room is too short, make it a hallway or backrooms instead
	cube_t const orig_room(room);
	room = exp_room;
	unsigned num_floors_add(0);
	float const room_floor_spacing((is_mall ? MALL_FLOOR_HEIGHT : 1.0)*floor_spacing); // mall has larger (2x) floor spacing

	if (is_mall) {
		unsigned const max_levels(global_building_params.max_mall_levels);
		assert(max_levels > 0);
		if (max_levels <= 2) {num_floors_add = max_levels - 1;} // 1-2 levels
		else {num_floors_add = rgen.rand_uniform_uint(2, max_levels) - 1;} // 2-max_levels levels

		if (room_floor_spacing > floor_spacing) { // check if we can lower the floor to increase room height
			float mall_floor_z(room.z1() - (room_floor_spacing - floor_spacing)); // of upper level
			cube_t cand(room);
			set_cube_zvals(cand, mall_floor_z, room.z1()); // one floor below
			// shift upward if this is a lower basement floor; shift is limited by the min of the gap above the mall celing and the gap below the basement floor (stairs height)
			float const max_shift_up(min((ground_floor_z1 - room.z2()), (get_basement().z1() - mall_floor_z)));
			bool valid(0);

			if (max_shift_up > get_floor_thickness()) {
				cube_t cand_shift(room);
				set_cube_zvals(cand_shift, room.z2(), (room.z2() + max_shift_up)); // old to new ceiling
				float const ceiling_zval(cand_shift.z2() - get_fc_thickness());

				if (!check_buildings_cube_coll(cand_shift, 0, 1, this) && query_min_height(room, ceiling_zval) > ceiling_zval) { // can shift up
					room.z2()    += max_shift_up;
					mall_floor_z += max_shift_up;
					valid = 1;
				}
			}
			if (!valid && check_buildings_cube_coll(cand, 0, 1, this)) {room = orig_room; return 0;} // not enough space below
			room.z1() = mall_floor_z;
		}
	}
	else {
		float const max_depth(room.z2() - get_max_sea_level());
		unsigned const max_num_floors(max(1U, min(global_building_params.max_ext_basement_room_depth, unsigned(floor(max_depth/floor_spacing)))));
		if (max_num_floors > 1) {num_floors_add = rgen.rand() % max_num_floors;} // maybe expand downward for additional floors
	}
	for (unsigned n = 0; n < num_floors_add; ++n) {
		cube_t cand(room);
		set_cube_zvals(cand, room.z1()-room_floor_spacing, room.z1()); // one floor below
		if (check_buildings_cube_coll(cand, 0, 1, this)) {num_floors_add = n; break;} // check for ext basement and tunnels below; xy_only=0, inc_basement=1, exclude ourself
		room.z1() = cand.z1();
	}
	return (num_floors_add + 1); // return the total number of floors
}

cube_t building_t::add_ext_basement_door(cube_t const &room, float door_width, bool dim, bool dir, bool is_end_room, bool is_tall_room, rand_gen_t &rgen, bool opens_other_side) {
	float const fc_thick(get_fc_thickness());
	cube_t door;
	set_cube_zvals(door, room.z1()+fc_thick, room.z2()-fc_thick);
	set_wall_width(door, room.get_center_dim(!dim), 0.5*door_width, !dim);
	door.d[dim][0] = door.d[dim][1] = room.d[dim][dir]; // one end of the room
	door_t Door(door, dim, (!dir ^ opens_other_side), rgen.rand_bool()); // open 50% of the time
	if (is_tall_room) {Door.set_mult_floor();}
	add_interior_door(Door, 0, !is_end_room); // is_bathroom=0, make_unlocked=!is_end_room
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
		cube_t const door(add_ext_basement_door(room, door_width, dim, d, is_end_room, 0, rgen)); // is_tall_room=0
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
				stairs.expand_in_dim(!dim, - wall_half_thick); // shrink on the sides
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
	if (room.has_stairs) {Room.has_stairs = 255;} // stairs on the first/only floor, or all floors if backroom
	if (is_hallway) {Room.assign_all_to(RTYPE_HALL, 0);} // initially all hallways; locked=0
	bool const is_first_extb_room((int)rooms.size() == ext_basement_hallway_room_id);
	rooms.push_back(Room);
	cube_t ceiling(room), floor(room);
	ceiling.z1() = room.z2() - fc_thick;
	floor  .z2() = room.z1() + fc_thick;
	P.subtract_stairs_from_floor_ceil(ceiling, P.wall_segs);
	
	if (is_first_extb_room && has_mall()) { // subtract mall elevator shaft and skylights from mall concourse ceiling
		if (mall_info->city_elevator_ix >= 0) {subtract_cube_from_cubes(elevators[mall_info->city_elevator_ix], P.wall_segs);}
		for (cube_t const &skylight : mall_info->skylights) {subtract_cube_from_cubes(skylight, P.wall_segs);}
	}
	vector_add_to(P.wall_segs, ceilings);
	P.subtract_stairs_from_floor_ceil(floor, P.wall_segs);
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

unsigned building_t::get_ext_basement_floor_ix(float zval) const {
	assert(has_ext_basement());
	return unsigned(max(0.0f, (zval - interior->basement_ext_bcube.z1())/get_window_vspace()));
}
void building_t::get_pgbr_wall_ix_for_pos(point const &pos, index_pair_t &start, index_pair_t &end) const { // pos is in building space
	if (!has_room_geom() || !is_pos_in_pg_or_backrooms(pos)) return;
	auto const &pgbr_wall_ixs(interior->room_geom->pgbr_wall_ixs);

	if (get_basement().contains_pt(pos)) { // inside parking garage
		if (pgbr_wall_ixs.empty()) {end = index_pair_t(interior->room_geom->pgbr_walls);} // not using indices, so use full range
		else {end = pgbr_wall_ixs.front();} // ends at first index (backrooms)
	}
	else if (has_ext_basement() && interior->basement_ext_bcube.contains_pt(pos)) { // inside backrooms
		unsigned const floor_ix(get_ext_basement_floor_ix(pos.z)); // floor containing pos.z

		if (floor_ix+1 < pgbr_wall_ixs.size()) { // if outside the valid floor range, start==end, the range will be empty, and we skip all walls
			start = pgbr_wall_ixs[floor_ix];
			end   = pgbr_wall_ixs[floor_ix+1];
		}
	}
}
bool building_t::point_in_extended_basement(point const &pos) const {
	if (!has_basement() || pos.z > ground_floor_z1 || !interior) return 0;
	if (interior->basement_ext_bcube.contains_pt(pos))           return 1;
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
cube_t building_t::get_ext_basement_entrance() const {
	assert(interior);
	assert(interior->ext_basement_hallway_room_id >= 0);
	door_t const &door(interior->get_ext_basement_door());
	cube_t room(*interior->ext_basement_rooms_start()); // hallway or backrooms
	float const expand(max(get_wall_thickness(), get_scaled_player_radius()));
	// clamp to padded door bounds to prevent objects from passing through the walls adjacent to the door
	max_eq(room.d[!door.dim][0], door.d[!door.dim][0]-expand);
	min_eq(room.d[!door.dim][1], door.d[!door.dim][1]+expand);
	assert(room.is_strictly_normalized());
	return room;
}

vector<room_t>::const_iterator building_interior_t::ext_basement_rooms_start() const {
	if (ext_basement_hallway_room_id < 0) return rooms.end(); // no ext basement rooms
	assert((unsigned)ext_basement_hallway_room_id < rooms.size());
	return rooms.begin() + ext_basement_hallway_room_id;
}
bool building_interior_t::point_in_ext_basement_room(point const &pos, float floor_spacing, float expand) const {
	if (ext_basement_hallway_room_id < 0)                 return 0; // no ext basement rooms
	if (!basement_ext_bcube.contains_pt_exp(pos, expand)) return 0;

	for (auto r = ext_basement_rooms_start(); r != rooms.end(); ++r) {
		if (r->contains_pt_exp(pos, expand)) return 1;
	}
	if (point_in_tunnel(pos, expand)) return 1;
	if (pool.valid && pool.contains_pt_exp(pos, expand)) return 1;

	if (has_mall()) { // check player in mall elevator or U-shaped stairs, which may be outside building rooms
		for (cube_t const &c : mall_info->ext_stairs_elevators) {
			if (c.contains_pt_exp(pos, expand)) return 1;
		}
		if (point_in_U_stairwell(pos, floor_spacing, 1)) return 1; // in_mall=1; handles landing, but doesn't handle expand
	}
	return 0;
}
// returns true if cube is completely contained in any single room; tunnels are ignored
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
	assert(ground_floor_z1 == b.ground_floor_z1); // must be at same elevation
	float const floor_spacing(get_window_vspace());
	if (b.get_window_vspace() != floor_spacing) return; // can't connect if floor spacing differs; can happen with city office buildings
	float const wall_thickness(get_wall_thickness()), z_toler(0.1*get_trim_thickness()), doorway_width(get_doorway_width());
	float const wall_hwidth(0.8*doorway_width), min_shared_wall_len(2.01*(wall_hwidth + 2.0*wall_thickness));
	float const max_connect_dist(EXT_BASEMENT_JOIN_DIST*floor_spacing), min_connect_dist(2.1*doorway_width); // need enough space to fit two open doors
	cube_t const &other_eb_bc(b.interior->basement_ext_bcube);
	vector<room_t> const &rooms1(interior->rooms), &rooms2(b.interior->rooms);
	unsigned const rstart1(interior->ext_basement_hallway_room_id), rstart2(b.interior->ext_basement_hallway_room_id);
	assert(rstart1 < rooms1.size() && rstart2 < rooms2.size());
	auto r1_begin(rooms1.begin() + rstart1), r2_begin(rooms2.begin() + rstart2);
	ext_basement_room_params_t P, Pb, Padd; // P=input rooms for *this, Pb=input rooms for b, Padd=new rooms output for *this
	rand_gen_t rgen;
	rgen.set_state(rooms1.size(), rooms2.size());

	// find nearby candidate rooms
	for (auto r1 = r1_begin; r1 != rooms1.end(); ++r1) {
		cube_t search_area(*r1);
		search_area.expand_by(max_connect_dist);
		if (!search_area.intersects(other_eb_bc)) continue; // too far
		if (interior->has_backrooms && has_water() && r1->z1() < interior->water_zval) continue; // don't connect if underwater

		for (auto r2 = r2_begin; r2 != rooms2.end(); ++r2) {
			if (!search_area.intersects(*r2))        continue; // too far
			if (fabs(r1->z1() - r2->z1()) > z_toler) continue; // different floors/levels; do we need to check toler?
			if (b.interior->has_backrooms && b.has_water() && r1->z1() < b.interior->water_zval) continue; // don't connect if underwater
			
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
				cand_join.z2() = cand_join.z1() + floor_spacing; // make it exactly one floor, in case this room connects to a tall pool room
				cand_join.d[d][ dir] = r2->d[d][!dir];
				cand_join.d[d][!dir] = r1->d[d][ dir];
				if (cand_join.get_sz_dim(d) < min_connect_dist) continue;
				set_wall_width(cand_join, door_center, wall_hwidth, !d);
				assert(cand_join.is_strictly_normalized());
				cube_t test_cube(cand_join);
				test_cube.expand_in_dim(d, -wall_thickness); // shrink ends to avoid false intersection with rooms at either end
				// check for intersections with starting hallways, since these aren't valid to connect to and aren't included in populate_params_from_building()
				if (test_cube.intersects(*r1_begin) || test_cube.intersects(*r2_begin)) continue;
				populate_params_from_building(*  interior, P );
				populate_params_from_building(*b.interior, Pb);
				if (!  is_basement_room_placement_valid(test_cube, P,  d,  dir, nullptr, &b  )) continue; // add_end_door=nullptr
				if (!b.is_basement_room_placement_valid(test_cube, Pb, d, !dir, nullptr, this)) continue; // add_end_door=nullptr
				Padd.rooms.emplace_back(cand_join, 1, 0, d, dir); // is_hallway=1, has_stairs=0
				Padd.rooms.back().conn_bcube = *r2; // store room in the other building that we're connecting to in conn_bcube
				if (r2 == r2_begin) {b.interior->conn_room_in_extb_hallway = 1;} // flag if connected to ext basement starting room
			} // for d
		} // for r2
	} // for r1
	if (Padd.rooms.empty()) return; // failed to connect
	unsigned ds_start[2] = {};
	building_t *const buildings[2] = {this, &b};

	for (unsigned bix = 0; bix < 2; ++bix) {
		ds_start[bix] = buildings[bix]->interior->door_stacks.size();
		if (!buildings[bix]->interior->conn_info) {buildings[bix]->interior->conn_info.reset(new building_conn_info_t);}
	}
	for (auto const &r : Padd.rooms) { // add any new rooms from above
		unsigned const is_building_conn(r.hallway_dim ? 2 : 1);
		// skip one end in hallway_dim and make the other end (bordering the other building) thinner to avoid Z-fighting but still cast shadows
		interior->place_exterior_room(r, r, get_fc_thickness(), wall_thickness, P, basement_part_ix, 0, r.is_hallway, is_building_conn, r.hallway_dim, r.connect_dir);
		unsigned const conn_door_ix(b.interior->doors.size()); // index of door that will be added to the other building, and separates the two buildings
		
		// place doors at each end
		for (unsigned dir = 0; dir < 2; ++dir) {
			building_t *door_dest(buildings[bool(dir) ^ r.connect_dir ^ 1]); // add door to the building whose room it connects to
			cube_t const door(door_dest->add_ext_basement_door(r, doorway_width, r.hallway_dim, dir, 0, 0, rgen)); // is_end_room=0, is_tall_room=0
			// subtract door from walls of each building
			for (unsigned bix = 0; bix < 2; ++bix) {subtract_cube_from_cubes(door, buildings[bix]->interior->walls[r.hallway_dim], nullptr, 1);} // no holes, clip_in_z=1
		} // for dir
		b.interior->doors      .back().set_bldg_conn(); // door added to the other building, and separates the two buildings
		b.interior->door_stacks.back().set_bldg_conn();
		cube_t ext_bcube(r);
		ext_bcube.d[r.hallway_dim][r.connect_dir] = r.conn_bcube.d[r.hallway_dim][r.connect_dir]; // extend to cover the entire width of the adjacent hallway in the other building

		for (unsigned bix = 0; bix < 2; ++bix) { // connect both ways
			// door belongs to b, which is the first building passed in
			buildings[bix]->interior->conn_info->add_connection(buildings[!bix], r, conn_door_ix, r.hallway_dim, r.connect_dir, (bix == 0));
			buildings[bix]->interior->basement_ext_bcube.union_with_cube(ext_bcube);
		}
	} // for r
	for (unsigned bix = 0; bix < 2; ++bix) {buildings[bix]->finalize_extb_conn_rooms(ds_start[bix]);}
}

void building_t::finalize_extb_conn_rooms(unsigned ds_start) {
	assert(interior);
	interior->assign_door_conn_rooms(ds_start); // assign room connections to any doors that were added

	if (has_pool()) { // check for tall pool rooms and add extra wall segments above the door; maybe just clip walls in Z instead?
		int const room_ix(interior->pool.room_ix);
		room_t const &pool_room(get_room(room_ix));

		if (pool_room.dz() > get_window_vspace()) { // tall pool room
			for (auto d = interior->door_stacks.begin()+ds_start; d != interior->door_stacks.end(); ++d) {
				for (unsigned s = 0; s < 2; ++s) { // for each side of the door
					if ((int)d->conn_room[s] != room_ix) continue; // not the pool room
					add_wall_section_above_pool_room_door(*d, pool_room);
				}
			}
		}
	}
	interior->remove_excess_capacity(); // optional optimization
}

void try_join_city_building_ext_basements(vect_building_t &buildings) {
	//timer_t timer("Join Building Basements"); // ~10ms
	vector<vector<unsigned>> bldgs_by_city;

	for (auto b = buildings.begin(); b != buildings.end(); ++b) {
		if (!b->is_in_city || !b->has_ext_basement()) continue;
		if (b->city_ix >= bldgs_by_city.size()) {bldgs_by_city.resize(b->city_ix+1);}
		bldgs_by_city[b->city_ix].push_back(b - buildings.begin());
	}
	for (vector<unsigned> const &work : bldgs_by_city) { // could be run in parallel, but not needed
		// do a quadratic iteration to find nearby buildings in this city that can potentially be connected
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

	if (player_in_basement < 3) { // not in extended basement - we can do some other checks
		if (!player_building->interior || !player_building->interior->conn_room_in_extb_hallway) return 0; // shouldn't be visible
		door_t const &door(player_building->interior->get_ext_basement_door());
		if (!door.can_see_through()) return 0; // fully closed opaque door
		cube_t dbc(door.get_true_bcube());
		dbc.expand_in_dim(door.dim, get_wall_thickness()); // expand a bit to handle player in the doorway
		if (!camera_pdu.cube_visible(dbc + xlate)) return 0; // check ext basement entrance visible
	}
	float const view_dist(8.0*get_window_vspace()); // arbitrary constant, should reflect length of largest hallway
	return player_building->is_visible_through_conn(*this, xlate, view_dist, expand_for_light);
}
cube_t building_t::get_conn_room_closest_to(point const &pos_bs) const { // in reference to the player's current building
	if (!player_in_basement || player_building == nullptr || player_building == this || !has_conn_info()) return cube_t();
	return interior->conn_info->get_conn_room_closest_to(*this, *player_building, pos_bs);
}
bool building_t::point_in_extb_conn_room(point const &pos_bs) const {
	return (interior->conn_info && interior->conn_info->point_in_conn_room(pos_bs));
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
	float const min_dist(target.get_wall_thickness());

	for (conn_pt_t const &c : conn) {
		if (c.b != &target) continue; // skip wrong building

		for (conn_room_t const &room : c.rooms) {
			cube_t room_cs(room + xlate);
			if (!room_cs.closest_dist_less_than(camera_pdu.pos, view_dist)) continue; // too far away
			if ( room_cs.closest_dist_less_than(camera_pdu.pos, min_dist )) return 1; // in doorway
			if (expand_for_light) {room_cs.expand_by(view_dist);} // increase the bounds in case room is behind the player but light cast from it is visible
			if (!camera_pdu.cube_visible(room_cs)) continue;
			if (!expand_for_light) return 1; // can't ignore closed doors for room objects because they we may not draw the door itself
			// if this is a light, check if the connecting door is open
			door_t const &door((room.door_is_b ? target : parent).get_door(room.door_ix));
			if (door.open || door.can_see_through()) return 1; // true if either about to open or not fully closed
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
bool building_conn_info_t::point_in_conn_room(point const &pos_bs) const {
	for (conn_pt_t const &c : conn) {
		if (check_vect_cube_contains_pt(c.rooms, pos_bs)) return 1;
	}
	return 0;
}


