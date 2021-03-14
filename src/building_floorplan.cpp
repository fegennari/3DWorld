// 3D World - Building Interior Floorplan Generation (walls, ceilings, floors, rooms, and doors)
// by Frank Gennari 2/20/20

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "profiler.h"

extern building_params_t global_building_params;


void building_t::add_interior_door(door_t &door) {
	assert(interior);
	interior->door_stacks.push_back(door);
	if (!SPLIT_DOOR_PER_FLOOR || door.on_stairs) {add_interior_door_for_floor(door); return;} // add a single door across all floors
	float const floor_spacing(get_window_vspace()), door_height(floor_spacing - get_floor_thickness());

	// Note: door.dz() should be an exact multiple of floor_spacing except for an extra floor thickness at the bottom
	for (float zval = door.z1(); zval + 0.5f*floor_spacing < door.z2(); zval += floor_spacing) { // continue until we don't have enough space left to add a door
		door_t door_seg(door);
		set_cube_zvals(door_seg, zval, zval+door_height); // clip to ceiling
		add_interior_door_for_floor(door_seg);
	}
}
void building_t::add_interior_door_for_floor(door_t &door) {
	if (!door.on_stairs) { // don't set open/locked state for stairs doors
		door.open   = (              fract(interior->doors.size()*1.61803) < global_building_params.open_door_prob  ); // use the golden ratio
		door.locked = (!door.open && fract(interior->doors.size()*3.14159) < global_building_params.locked_door_prob); // use pi
	}
	interior->doors.push_back(door);
}

void building_t::remove_section_from_cube_and_add_door(cube_t &c, cube_t &c2, float v1, float v2, bool xy, bool open_dir) {
	// remove a section from this cube; c is input+output cube, c2 is other output cube
	assert(v1 > c.d[xy][0] && v1 < v2 && v2 < c.d[xy][1]); // v1/v2 must be interior values for cube
	c2 = c; // clone first cube
	c.d[xy][1] = v1; c2.d[xy][0] = v2; // c=low side, c2=high side
	// add a door
	door_t door(c, !xy, open_dir);
	door.d[!xy][0] = door.d[!xy][1] = c.get_center_dim(!xy); // zero area at wall centerline
	door.d[ xy][0] = v1; door.d[ xy][1] = v2;
	add_interior_door(door);
}

void building_t::insert_door_in_wall_and_add_seg(cube_t &wall, float v1, float v2, bool dim, bool open_dir, bool keep_high_side) {
	cube_t wall2;
	remove_section_from_cube_and_add_door(wall, wall2, v1, v2, dim, open_dir);
	if (keep_high_side) {swap(wall, wall2);} // swap left and right
	interior->walls[!dim].push_back(wall2);
}

float cube_rand_side_pos(cube_t const &c, bool dim, float min_dist_param, float min_dist_abs, rand_gen_t &rgen) {
	assert(min_dist_param < 0.5f); // aplies to both ends
	float const lo(c.d[dim][0]), hi(c.d[dim][1]), delta(hi - lo), gap(max(min_dist_abs, min_dist_param*delta)), v1(lo + gap), v2(hi - gap);
	//if (v2 <= v1) {cout << TXT(dim) << TXT(lo) << TXT(hi) << TXT(min_dist_abs) << TXT(delta) << TXT(gap) << endl;}
	//assert(v1 <= v2); // too strong?
	if (v1 >= v2) {return 0.5f*(v1 + v2);} // if range is denormalized, use the center
	return rgen.rand_uniform(v1, v2);
}

// see global_building_params.window_xspace/window_width
int building_t::get_num_windows_on_side(float xy1, float xy2) const {
	assert(xy1 < xy2);
	float const tscale(get_material().get_floorplan_window_xscale());
	float t0(tscale*xy1), t1(tscale*xy2);
	clip_low_high(t0, t1);
	return round_fp(t1 - t0);
}
float building_t::get_window_h_border() const {return 0.5*(1.0 - global_building_params.get_window_width_fract ());} // (0.0, 1.0)
float building_t::get_window_v_border() const {return 0.5*(1.0 - global_building_params.get_window_height_fract());} // (0.0, 1.0)
float building_t::get_hspacing_for_part(cube_t const &part, bool dim) const {return part.get_sz_dim(dim)/get_num_windows_on_side(part.d[dim][0], part.d[dim][1]);}

void set_wall_width(cube_t &wall, float pos, float half_thick, unsigned dim) {
	wall.d[dim][0] = pos - half_thick;
	wall.d[dim][1] = pos + half_thick;
}
void clip_wall_to_ceil_floor(cube_t &wall, float fc_thick) {
	wall.z1() += fc_thick; // start at the floor
	wall.z2() -= fc_thick; // start at the ceiling
}
// Note: wall should start out equal to the room bcube
void create_wall(cube_t &wall, bool dim, float wall_pos, float fc_thick, float wall_half_thick, float wall_edge_spacing) {
	clip_wall_to_ceil_floor(wall, fc_thick);
	set_wall_width(wall, wall_pos, wall_half_thick, dim);
	// move a bit away from the exterior wall to prevent z-fighting; we might want to add walls around the building exterior and cut window holes
	wall.expand_in_dim(!dim, -wall_edge_spacing);
}

// Note: assumes edge is not clipped and doesn't work when clipped
bool is_val_inside_window(cube_t const &c, bool dim, float val, float window_spacing, float window_border) {
	window_border *= 0.9; // adjust based on window frame so that wall doesn't end right at the edge
	float const uv(fract((val - c.d[dim][0])/window_spacing));
	return (uv > window_border && uv < 1.0f-window_border);
}
float shift_val_to_not_intersect_window(cube_t const &c, float val, float hspace, float window_border, bool dim) {
	window_border *= 0.9; // adjust based on window frame so that wall doesn't end right at the edge
	float const uv(fract((val - c.d[dim][0])/hspace));
	if (!(uv > window_border && uv < 1.0f-window_border)) return val; // okay as is
	float const uv_target((uv < 0.5) ? window_border : 1.0f-window_border);
	return val + (uv_target - uv)*hspace; // shift val exactly to the window border in the closer direction
}

struct split_cube_t : public cube_t {
	float door_lo[2][2], door_hi[2][2]; // per {dim x dir}
	
	split_cube_t(cube_t const &c) : cube_t(c) {
		door_lo[0][0] = door_lo[0][1] = door_lo[1][0] = door_lo[1][1] = door_hi[0][0] = door_hi[0][1] = door_hi[1][0] = door_hi[1][1] = 0.0f;
	}
	bool bad_pos(float val, bool dim) const {
		for (unsigned d = 0; d < 2; ++d) { // check both dirs (wall end points)
			if (door_lo[dim][d] < door_hi[dim][d] && val > door_lo[dim][d] && val < door_hi[dim][d]) return 1;
		}
		return 0;
	}
};

unsigned calc_num_floors(cube_t const &c, float window_vspacing, float floor_thickness) {
	float const z_span(c.dz() - floor_thickness);
	assert(z_span > 0.0);
	unsigned const num_floors(round_fp(z_span/window_vspacing)); // round - no partial floors
	if (num_floors > 1000) {std::cerr << "Error: building with more than 1000 floors (" << num_floors << ")" << endl; exit(0);} // sanity check
	return num_floors;
}

void subtract_cube_xy(cube_t const &c, cube_t const &r, cube_t *out) { // subtract r from c; ignores zvals
	assert(c.contains_cube_xy(r));
	for (unsigned i = 0; i < 4; ++i) {out[i] = c;}
	out[0].y2() = r.y1(); // bottom -y
	out[1].y1() = r.y2(); // top    +y
	out[2].y1() = r.y1(); out[2].y2() = r.y2(); out[2].x2() = r.x1(); // left  center -x
	out[3].y1() = r.y1(); out[3].y2() = r.y2(); out[3].x1() = r.x2(); // right center +x
}

bool building_t::interior_enabled() const {
	if (!ADD_BUILDING_INTERIORS) return 0; // disabled
	if (world_mode != WMODE_INF_TERRAIN) return 0; // tiled terrain mode only
	if (!global_building_params.windows_enabled()) return 0; // no windows, can't assign floors and generate interior
	//if (has_overlapping_cubes) return; // overlapping cubes buildings are more difficult to handle
	if (!is_cube()) return 0; // only generate interiors for cube buildings for now
	if (!global_building_params.add_city_interiors && !get_material().add_windows) return 0; // not a building type that has generated windows (skip buildings with windows baked into textures)
	return 1;
}

int building_t::classify_room_wall(room_t const &room, float zval, bool dim, bool dir, bool ret_sep_if_part_int_part_ext) const { // Note: zval is for the floor
	if (room.d[dim][dir] == bcube.d[dim][dir]) return ROOM_WALL_EXT; // at bcube border
	if (!ret_sep_if_part_int_part_ext && (room.ext_sides & (1 << (2*dim + dir)))) return ROOM_WALL_EXT; // use ext_sides (optimization, but may conservatively return ext)
	cube_t const &part(get_part_for_room(room));
	float const wall_thickness(get_wall_thickness()), max_gap(1.5*wall_thickness), part_edge(part.d[dim][dir]);
	if (dir ? (room.d[dim][1] + max_gap < part_edge) : (room.d[dim][0] - max_gap > part_edge)) return ROOM_WALL_INT; // interior to part, allowing for wall_thickness of gap
	if (real_num_parts == 1) return ROOM_WALL_EXT; // optimization

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		//if (is_basement(p)) {} // special window-less wall type for basements?
		if (p->d[dim][!dir] != part_edge) continue; // not opposite wall
		if (p->z1() >= room.z2() || p->z2() <= room.z1()) continue; // no z overlap (wrong stack)

		if (ret_sep_if_part_int_part_ext) { // return sep if any part of the wall separates parts
			if (p->d[!dim][0] >= room.d[!dim][1] || p->d[!dim][1] <= room.d[!dim][0]) continue; // wall not overlapping
			return ROOM_WALL_SEP;
		}
		else { // only return sep if the entire wall separates parts (none is exterior)
			if (p->d[!dim][0] > room.d[!dim][0] || p->d[!dim][1] < room.d[!dim][1]) continue; // wall not contained
			if (zval + get_floor_thickness() < p->z2()) return ROOM_WALL_SEP; // this part covers the wall in z (assuming no overhangs), so wall is interior split between parts
		}
	} // for p
	return ROOM_WALL_EXT; // exterior
}

void building_t::gen_interior(rand_gen_t &rgen, bool has_overlapping_cubes) { // Note: contained in building bcube, so no bcube update is needed
	if (!interior_enabled()) return;

	// make up to 20 attempts to generate a connected interior; the first attempt almost always succeeds; currently it only fails if a house basement is unconnected
	// 64 houses are unconnected, this loop fixes all but 6 of them after 3 iterations; 12 iterations is required for 100% success
	for (unsigned n = 0; n < 20; ++n) {
		gen_interior_int(rgen, has_overlapping_cubes);
		if (!interior->is_unconnected) break; // done
	}
	interior->finalize();
}

void building_t::gen_interior_int(rand_gen_t &rgen, bool has_overlapping_cubes) { // Note: contained in building bcube, so no bcube update is needed

	// defer this until the building is close to the player?
	interior.reset(new building_interior_t);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const doorway_width(0.5*window_vspacing), doorway_hwidth(0.5*doorway_width);
	float const wall_thick(get_wall_thickness()), wall_half_thick(0.5*wall_thick), wall_edge_spacing(0.05*wall_thick), min_wall_len(4.0*doorway_width);
	float const window_border(get_window_h_border());
	point bldg_door_open_dir_tp(bcube.get_cube_center()); // used to determine in which direction doors open; updated base on central hallway
	// houses have at most two parts; exclude garage, shed, porch, porch support, etc.
	auto parts_end(get_real_parts_end());
	vector<split_cube_t> to_split;
	uint64_t must_split[2] = {0,0};
	unsigned first_wall_to_split[2] = {0,0};
	// allocate space for all floors
	unsigned tot_num_floors(0), tot_num_stairwells(0), tot_num_landings(0); // num floor/ceiling cubes, not number of stories; used only for reserving vectors

	for (auto p = parts.begin(); p != parts_end; ++p) {
		bool const has_stairs(!is_house || p == parts.begin()); // assumes one set of stairs or elevator per part
		unsigned const cubes_per_floor(has_stairs ? 4 : 1); // account for stairwell cutouts
		unsigned const num_floors(calc_num_floors(*p, window_vspacing, floor_thickness));
		tot_num_floors     += cubes_per_floor*(num_floors - 1) + 1; // first floor has no cutout
		tot_num_stairwells += (has_stairs && num_floors > 1);
		tot_num_landings   += (has_stairs ? (num_floors - 1) : 0);
	}
	if (has_sec_bldg()) {++tot_num_floors;}
	interior->ceilings.reserve(tot_num_floors);
	interior->floors  .reserve(tot_num_floors);
	interior->landings.reserve(tot_num_landings);
	interior->stairwells.reserve(tot_num_stairwells);
	
	// generate walls and floors for each part;
	// this will need to be modified to handle buildings that have overlapping parts, or skip those building types completely
	for (auto p = parts.begin(); p != parts_end; ++p) {
		unsigned const num_floors(calc_num_floors(*p, window_vspacing, floor_thickness));
		if (num_floors == 0) continue; // not enough space to add a floor (can this happen?)
		// for now, assume each part has the same XY bounds and can use the same floorplan; this means walls can span all floors and don't need to be duplicated for each floor
		vector3d const psz(p->get_size());
		bool const min_dim(psz.y < psz.x); // hall dim
		float const cube_width(psz[min_dim]);
		bool const first_part(p == parts.begin()), first_part_this_stack(first_part || (p-1)->z1() < p->z1());
		bool const use_hallway(!is_house && !has_complex_floorplan && first_part_this_stack && (p+1 == parts.end() || (p+1)->z1() > p->z1()) && cube_width > 4.0*min_wall_len);
		unsigned const rooms_start(interior->rooms.size()), part_id(p - parts.begin());
		cube_t hall, place_area(*p);
		place_area.expand_by_xy(-wall_edge_spacing); // shrink slightly to avoid z-fighting with walls
		float window_hspacing[2] = {0.0};
		int num_windows_per_side[2] = {0};

		for (unsigned d = 0; d < 2; ++d) {
			num_windows_per_side[d] = get_num_windows_on_side(p->d[d][0], p->d[d][1]);
			window_hspacing[d] = psz[d]/num_windows_per_side[d];
		}
		if (use_hallway) {
			// building with rectangular slice (no adjacent exterior walls at this level), generate rows of offices
			// Note: we could probably make these unsigned, but I want to avoid unepected negative numbers in the math
			int const num_windows   (num_windows_per_side[!min_dim]);
			int const num_windows_od(num_windows_per_side[min_dim]); // other dim, for use in hallway width calculation
			int windows_per_room((num_windows >= 7 && num_windows_od >= 7) ? 2 : 1); // 1-2 windows per room (only assign 2 windows if we can get into the secondary hallway case below)
			float const cube_len(psz[!min_dim]), wind_hspacing(cube_len/num_windows), min_hall_width(3.6*doorway_width);
			float room_len(wind_hspacing*windows_per_room);

			while (room_len < 0.9*min_wall_len) { // add more windows to increase room size if too small
				++windows_per_room;
				room_len = wind_hspacing*windows_per_room;
			}
			int const num_rooms((num_windows+windows_per_room-1)/windows_per_room); // round up
			bool const partial_room((num_windows % windows_per_room) != 0); // an odd number of windows leaves a small room at the end
			assert(num_rooms >= 0 && num_rooms < 1000); // sanity check
			float num_hall_windows((num_windows_od & 1) ? 1.4 : 1.8); // hall either contains 1 (odd) or 2 (even) windows, wider for single window case to make room for stairs
			max_eq(num_hall_windows, min_hall_width*num_windows_od/cube_width); // enforce min_hall_width (may split a window, but this limit is only hit for non-window city office buildings)
			float const hall_width(num_hall_windows*cube_width/num_windows_od);
			float const room_width(0.5f*(cube_width - hall_width)); // rooms are the same size on each side of the hallway
			float const hwall_extend(0.5f*(room_len - doorway_width - wall_thick));
			float const wall_pos(p->d[!min_dim][0] + room_len); // pos of first wall separating first from second rooms
			float const hall_wall_pos[2] = {(p->d[min_dim][0] + room_width), (p->d[min_dim][1] - room_width)};
			if (hallway_dim == 2) {hallway_dim = !min_dim;} // cache in building for later use, only for first part (ground floor)
			vect_cube_t &room_walls(interior->walls[!min_dim]), &hall_walls(interior->walls[min_dim]);
			hall = *p;
			for (unsigned e = 0; e < 2; ++e) {hall.d[min_dim][e] = hall_wall_pos[e];}
			
			if (num_windows_od >= 7 && num_rooms >= 4) { // at least 7 windows (3 on each side of hallway)
				float const min_hall_width(1.5f*doorway_width), max_hall_width(2.5f*doorway_width);
				float const sh_width(max(min(0.4f*hall_width, max_hall_width), min_hall_width)), hspace(window_hspacing[!min_dim]);

				if (rgen.rand_bool()) { // ring hallway
					float const room_depth(0.5f*(room_width - sh_width)); // for inner and outer rows of rooms
					assert(room_depth > 2.0f*doorway_width); // I'm not sure if this can fail or what we should do in that case - no secondary hallways?
					float const hall_offset(room_len); // outer edge of hallway to building exterior on each end
					bool const add_doors_to_main_wall(rgen.rand_bool());
					unsigned const num_cent_rooms(num_rooms - 2); // skip the rooms on each side
					unsigned const num_doors_inner_rooms(add_doors_to_main_wall ? 3U : 2U);
					unsigned const num_offices(4*(num_cent_rooms + 4)); // X/Y mirror symmetry
					assert(num_cent_rooms > 0);
					hall_walls.reserve(2*(11 + num_doors_inner_rooms*num_cent_rooms)); // long dim (along hall dir)
					room_walls.reserve(2*(10 + 2*num_cent_rooms)); // short dim
					interior->rooms.reserve(num_offices + 7); // num_offices + pri hall + 2 sec hall + 4 conn hall
					interior->doors.reserve(num_offices + 2*add_doors_to_main_wall*num_cent_rooms + 4); // at least one per office
					interior->exclusion.reserve(6); // 2 sec hallways + 4 conn hallways
					unsigned const bathroom_ix(rgen.rand_bool() ? 0 : num_cent_rooms-1); // place bathrooms on one of the end center rooms

					for (unsigned d = 0; d < 2; ++d) { // for each side of main hallway
						float const dsign(d ? -1.0 : 1.0);
						float const wall_edge(p->d[min_dim][d]), targ_hall_outer(wall_edge + dsign*room_depth);
						float const hall_outer(shift_val_to_not_intersect_window(*p, targ_hall_outer, window_hspacing[min_dim], window_border, min_dim));
						float const hall_inner(hall_outer + dsign*sh_width), conn_hall_len(dsign*(hall_wall_pos[d] - hall_inner));
						float const targ_side_room_split(hall_outer + 0.5f*dsign*(conn_hall_len + sh_width)); // split room halfway
						float const side_room_split(shift_val_to_not_intersect_window(*p, targ_side_room_split, window_hspacing[min_dim], window_border, min_dim));
						cube_t s_hall(*p), c_hall(*p);
						s_hall.d[ min_dim][ d] = hall_outer;
						s_hall.d[ min_dim][!d] = hall_inner;
						s_hall.d[!min_dim][ 0] = p->d[!min_dim][0] + hall_offset;
						s_hall.d[!min_dim][ 1] = p->d[!min_dim][1] - hall_offset;
						c_hall.d[ min_dim][ d] = hall_inner;
						c_hall.d[ min_dim][!d] = hall_wall_pos[d];
						add_room(s_hall, part_id, 3, 1, 0); // add sec hallway as room with 3 lights
						interior->exclusion.push_back(s_hall); // excluded from placing stairs and elevators

						// walls along sec hallway
						cube_t long_swall(s_hall), short_swall(s_hall); // sec outer, sec inner; both have rows of doors in them
						clip_wall_to_ceil_floor(long_swall,  fc_thick);
						clip_wall_to_ceil_floor(short_swall, fc_thick);
						short_swall.expand_in_dim(!min_dim, -(sh_width - wall_half_thick)); // expand to fill the corner + shrink wall length
						long_swall.d[!min_dim][0] = place_area.d[!min_dim][0]; // crosses the entire building (creates sec hall + room on each end)
						long_swall.d[!min_dim][1] = place_area.d[!min_dim][1];
						cube_t main_wall(short_swall); // also short
						set_wall_width(main_wall,   hall.d[min_dim][d], wall_half_thick, min_dim);
						set_wall_width(long_swall,  hall_outer, wall_half_thick, min_dim);
						set_wall_width(short_swall, hall_inner, wall_half_thick, min_dim);

						// add inner (windowless) and outer rooms
						float const rooms_start(short_swall.d[!min_dim][0]), rooms_end(short_swall.d[!min_dim][1]);
						float const span(rooms_end - rooms_start), cent_room_width(span/num_cent_rooms);
						float room_pos(rooms_start);
						cube_t room_outer(*p), room_inner(*p);
						room_outer.d[min_dim][!d] = hall_outer; // other dim stays at wall_edge
						room_inner.d[min_dim][ d] = hall_inner;
						room_inner.d[min_dim][!d] = hall_wall_pos[d];

						for (unsigned e = 0; e < 2; ++e) { // for each end of building in long dim
							float const esign(e ? -1.0 : 1.0); // Note: offset_inner == ((e ? rooms_end : rooms_start) + esign*wall_half_thick)
							float const other_edge(p->d[!min_dim][e]), offset_outer(s_hall.d[!min_dim][e]), offset_inner(offset_outer + esign*sh_width);
							c_hall.d[!min_dim][ e] = offset_outer;
							c_hall.d[!min_dim][!e] = offset_inner;
							unsigned const num_lights((c_hall.get_sz_dim(min_dim) > 0.25*s_hall.get_sz_dim(!min_dim)) ? 2 : 1); // 2 lights if it's long enough
							add_room(c_hall, part_id, num_lights, 1, 0); // add conn hallway as room
							cube_t exclude(c_hall);
							exclude.d[min_dim][!d] += dsign*doorway_width; // expand out a bit into the main hallway to ensure there's space to enter this hallway
							interior->exclusion.push_back(exclude); // excluded from placing stairs and elevators

							for (unsigned side = 0; side < 2; ++side) { // add walls along connector hallway
								cube_t conn_wall(c_hall);
								clip_wall_to_ceil_floor(conn_wall, fc_thick);
								set_wall_width(conn_wall, c_hall.d[!min_dim][side], wall_half_thick, !min_dim);

								if (side == e) { // long wall, this is the one we want to extend add add doors to
									conn_wall.d[min_dim][d] = place_area.d[min_dim][d]; // extend long wall to edge of building to create corner room
									float const door_pos[3] = {0.5f*(wall_edge + hall_outer), 0.5f*(hall_outer + hall_inner), 0.5f*(hall_inner + hall_wall_pos[d])};

									for (unsigned n = 0; n < 3; ++n) { // add doors for corner room, end room, and side room
										float const lo_pos(door_pos[n] - doorway_hwidth), hi_pos(door_pos[n] + doorway_hwidth);
										insert_door_in_wall_and_add_seg(conn_wall, lo_pos, hi_pos, min_dim, e, !d);
									}
								}
								room_walls.push_back(conn_wall);
							} // for side
							cube_t end_wall(main_wall);
							end_wall.d[!min_dim][ e] = place_area.d[!min_dim][e]; // end of main hall/building
							end_wall.d[!min_dim][!e] = offset_outer + esign*wall_half_thick; // end of conn hall
							hall_walls.push_back(end_wall); // end of main hall
							set_wall_width(end_wall, side_room_split, wall_half_thick, min_dim);
							hall_walls.push_back(end_wall); // extension of sec hall short wall that separates corner office

							// add 2x end of hall, corner, and side rooms
							// Note: this should be the same as the shifted room_pos below so that the room bounds align with the wall positions
							float const shifted_offset_inner(shift_val_to_not_intersect_window(*p, (e ? rooms_end : rooms_start), hspace, window_border, !min_dim));
							cube_t end_room(room_outer), cor_room(room_outer); // z and min_dim vals are the same as the outer room
							cor_room.d[!min_dim][ e] = other_edge; // corner of building
							cor_room.d[!min_dim][!e] = offset_outer;
							add_room(cor_room, part_id, 1, 0, 1); // corner office
							end_room.d[!min_dim][ e] = offset_outer; // aligned with end of conn hall
							end_room.d[!min_dim][!e] = shifted_offset_inner;
							cube_t const conn_hall_end_room(end_room);
							add_room(end_room, part_id, 1, 0, 1); // office at end of conn hall
							end_room = cor_room; // next to corner room, copy !min_dim values
							end_room.d[min_dim][ d] = hall_outer;
							end_room.d[min_dim][!d] = side_room_split;
							add_room(end_room, part_id, 1, 0, 1); // office at end of sec hall
							cube_t side_room(end_room); // next to end_room
							side_room.d[min_dim][ d] = side_room_split;
							side_room.d[min_dim][!d] = hall_wall_pos[d];
							add_room(side_room, part_id, 1, 0, 1); // split this into multiple rooms if there's space? seems like there usually isn't
							float const door_pos[2] = {0.5f*(other_edge + offset_outer), 0.5f*(offset_outer + shifted_offset_inner)}; // {corner, conn hall end}

							for (unsigned n = 0; n < 2; ++n) { // add doors for corner room and side room on each side of long_swall
								float const lo_pos(door_pos[n] - doorway_hwidth), hi_pos(door_pos[n] + doorway_hwidth);
								// omit the door if there isn't enough space in the wall for it (wall was shifted too close to offset_outer)
								if (n == 1 && !(lo_pos-wall_half_thick > conn_hall_end_room.d[!min_dim][0] && hi_pos+wall_half_thick < conn_hall_end_room.d[!min_dim][1])) continue;
								insert_door_in_wall_and_add_seg(long_swall, lo_pos, hi_pos, !min_dim, d, !e);
							}
						} // for e
						for (unsigned n = 0; n <= num_cent_rooms; ++n) {
							float const start_pos(shift_val_to_not_intersect_window(*p, room_pos, hspace, window_border, !min_dim));

							if (n < num_cent_rooms) { // add a room
								float const next_pos(shift_val_to_not_intersect_window(*p, (room_pos + cent_room_width), hspace, window_border, !min_dim));
								room_pos += cent_room_width; // move to next row
								room_outer.d[!min_dim][0] = start_pos;
								room_inner.d[!min_dim][0] = ((n == 0) ? (rooms_start + wall_half_thick) : start_pos); // first inner room is relative to the sec hallway
								room_outer.d[!min_dim][1] = next_pos;
								room_inner.d[!min_dim][1] = ((n+1 == num_cent_rooms) ? (rooms_end - wall_half_thick) : next_pos); // last inner room is relative to sec hallway
								add_room(room_outer, part_id, 1, 0, 1); // office
								add_room(room_inner, part_id, 1, 0, 1); // office or bathroom
								if (n == bathroom_ix) {interior->rooms.back().assign_all_to(RTYPE_BATH);} // make it a bathroom
								// add doors to 2-3 walls
								float const door_pos(0.5f*(start_pos + next_pos)), lo_pos(door_pos - doorway_hwidth), hi_pos(door_pos + doorway_hwidth);
								cube_t *to_split[3] = {&long_swall, &short_swall, &main_wall};

								for (unsigned n = 0; n < num_doors_inner_rooms; ++n) {
									insert_door_in_wall_and_add_seg(*to_split[n], lo_pos, hi_pos, !min_dim, (d^(n&1)), 1);
								}
							}
							// add walls separating rooms
							cube_t wall_outer(room_outer);
							clip_wall_to_ceil_floor(wall_outer, fc_thick);
							wall_outer.d[min_dim][d] += dsign*wall_edge_spacing;
							set_wall_width(wall_outer, start_pos, wall_half_thick, !min_dim);
							room_walls.push_back(wall_outer);

							if (n > 0 && n < num_cent_rooms) { // can skip end walls adjacent to conn hallways
								cube_t wall_inner(room_inner);
								clip_wall_to_ceil_floor(wall_inner, fc_thick);
								set_wall_width(wall_inner, start_pos, wall_half_thick, !min_dim);
								room_walls.push_back(wall_inner);
							}
						} // for n
						hall_walls.push_back(main_wall  );
						hall_walls.push_back(long_swall ); // add remaining wall segment
						hall_walls.push_back(short_swall); // add remaining wall segment
					} // for d
				} // end ring hallways

				else { // secondary hallways with rooms on each side
					float const sh_len(room_width);
					int const num_sec_halls(num_rooms/2); // round down if odd
					int const windows_per_side_od((num_windows_od - round_fp(num_hall_windows))/2); // of hallway
					int windows_per_room_od((windows_per_side_od & 1) ? 1 : 2), rooms_per_side(0); // rooms along each sec hall
					float room_sub_width(0.0);

					while (1) { // increase the number of windows per room until room is large enough
						rooms_per_side = max(windows_per_side_od/windows_per_room_od, 1);
						room_sub_width = sh_len/rooms_per_side;
						if (room_sub_width > 1.5*doorway_width) break; // done
						++windows_per_room_od;
					}
					float const sh_spacing(cube_len/num_sec_halls - sh_width), end_spacing(0.5*sh_spacing); // half spacing at both ends
					assert(sh_spacing > min_wall_len); // I'm not sure if this can fail or what we should do in that case - use fewer secondary hallways?
					float room_start(p->d[!min_dim][0]), wall_pos(room_start + end_spacing); // first sec hall wall pos
					int const num_offices(4*rooms_per_side*num_sec_halls);
					vect_cube_t &split_walls(hall_walls);
					room_walls.reserve(4*(rooms_per_side+1)*num_sec_halls + 2*(num_sec_halls-1)); // walls with doors + room dividers
					split_walls.reserve(2*rooms_per_side*(num_sec_halls+1));
					interior->rooms.reserve(num_offices + 2*num_sec_halls + 1); // offices + sec hallways + pri hallway
					interior->doors.reserve(num_offices); // one per office
					interior->exclusion.reserve(2*num_sec_halls);
					unsigned const bathroom_ix((num_sec_halls <= 2) ? 0 : (rgen.rand()%(num_sec_halls-1))); // place bathrooms in rooms near the central hallway

					for (int i = 0; i <= num_sec_halls; ++i) { // actually iterates over the number of room blocks between halls (num halls + 1)
						// shift the position of the hallway walls to avoid intersecting windows;
						// this may make halls either larger (contain a window) or smaller (fit between two windows);
						// as long as the space between windows is large enough (> doorway_width) this should be okay
						float hall_start_pos(shift_val_to_not_intersect_window(*p, wall_pos,            hspace, window_border, !min_dim));
						float const target_hall_end_pos(wall_pos + sh_width); // use post-shifted wall for better fit to orig hall width
						float hall_end_pos  (shift_val_to_not_intersect_window(*p, target_hall_end_pos, hspace, window_border, !min_dim));
						float const sec_hall_width(hall_end_pos - hall_start_pos);
						bool const div_room(i > 0 && i < num_sec_halls), add_sec_hall(i < num_sec_halls);

						if (add_sec_hall && sec_hall_width < min_hall_width) { // hall is too narrow, try shifting start pos
							float const target_hall_start_pos(hall_start_pos - (min_hall_width - sec_hall_width));
							hall_start_pos = shift_val_to_not_intersect_window(*p, target_hall_start_pos, hspace, window_border, !min_dim);
						}
						if (add_sec_hall && sec_hall_width > max_hall_width) { // hall is too wide, try shifting start pos
							float const target_hall_start_pos(hall_start_pos - (max_hall_width - sec_hall_width));
							hall_start_pos = shift_val_to_not_intersect_window(*p, target_hall_start_pos, hspace, window_border, !min_dim);
						}
						for (unsigned d = 0; d < 2; ++d) { // left, right of main hall
							float const dsign(d ? -1.0 : 1.0);
							cube_t split_wall(*p);
							split_wall.d[!min_dim][0] = room_start     + wall_edge_spacing; // start of building or end of prev sec hall
							split_wall.d[!min_dim][1] = hall_start_pos - wall_edge_spacing; // start of this hall or end of building
							float room_split_pos(p->d[min_dim][d]), div_pos(0); // start at edge of building (outermost room with windows)
							cube_t div_wall(*p); // sets z1/z2
							clip_wall_to_ceil_floor(div_wall, fc_thick);

							if (div_room) { // divide this block of rooms into two rows, one facing the sec hallway on each side
								div_pos = shift_val_to_not_intersect_window(*p, split_wall.get_center_dim(!min_dim), hspace, window_border, !min_dim);
								set_wall_width(div_wall, div_pos, wall_half_thick, !min_dim);
							}
							div_wall.d[min_dim][ d] = place_area.d[min_dim][d];
							div_wall.d[min_dim][!d] = hall.d[min_dim][d];
							cube_t sep_walls[2];

							if (add_sec_hall) { // add sec hall
								cube_t s_hall(*p);
								s_hall.d[ min_dim][!d] = hall.d[min_dim][d]; // ends at main hall
								s_hall.d[!min_dim][ 0] = hall_start_pos;
								s_hall.d[!min_dim][ 1] = hall_end_pos;
								add_room(s_hall, part_id, 2, 1, 0); // add sec hallway as room with 2 lights
								cube_t exclude(s_hall);
								exclude.d[min_dim][!d] += dsign*doorway_width; // expand out a bit into the main hallway to ensure there's space to enter this hallway
								interior->exclusion.push_back(exclude); // excluded from placing stairs and elevators
								
								for (unsigned dir = 0; dir < 2; ++dir) { // add walls between hall and rooms on each side
									sep_walls[dir] = div_wall; // copy z and min_dim from div_wall
									sep_walls[dir].d[min_dim][!d] += dsign*wall_half_thick; // extend to meet the wall in the other dim
									set_wall_width(sep_walls[dir], s_hall.d[!min_dim][dir], wall_half_thick, !min_dim);
								}
							}
							for (int r = 0; r < rooms_per_side; ++r) {
								float const next_split_pos(room_split_pos + dsign*room_sub_width);
								cube_t room(split_wall);
								room.d[!min_dim][0] = room_start; // end exactly at part bcube
								room.d[!min_dim][1] = hall_start_pos; // end exactly at part bcube
								room.d[min_dim][ d] = room_split_pos;
								room.d[min_dim][!d] = next_split_pos;

								if (div_room) { // interior rooms, split into 2 rooms with a wall in between
									room.d[!min_dim][1] = div_pos; // left room
									add_room(room, part_id, 1, 0, 1); // office along sec hallway
									room.d[!min_dim][0] = div_pos; // right room
									room.d[!min_dim][1] = split_wall.d[!min_dim][1]; // restore orig value
								}
								add_room(room, part_id, 1, 0, 1); // office or bathroom along sec hallway
								if (r+1 == rooms_per_side && (unsigned)i == (bathroom_ix+1)) {interior->rooms.back().assign_all_to(RTYPE_BATH);} // bathroom must be an interior/windowless room

								if (add_sec_hall) { // add doorways + doors
									float const doorway_pos(0.5f*(room_split_pos + next_split_pos)); // room center
									float const lo_pos(doorway_pos - doorway_hwidth), hi_pos(doorway_pos + doorway_hwidth);

									for (unsigned dir = 0; dir < 2; ++dir) {
										insert_door_in_wall_and_add_seg(sep_walls[dir], lo_pos, hi_pos, min_dim, dir, !d);
									}
								}
								room_split_pos = next_split_pos;
								set_wall_width(split_wall, room_split_pos, wall_half_thick, min_dim);
								split_walls.push_back(split_wall);
								clip_wall_to_ceil_floor(split_walls.back(), fc_thick); // clip post-add so that the room continues to use part z1/z2
							} // for r
							if (add_sec_hall) { // add sec hall remaining walls after doorway insertion
								for (unsigned dir = 0; dir < 2; ++dir) {room_walls.push_back(sep_walls[dir]);}
							}
							if (div_room) {room_walls.push_back(div_wall);} // add a divider wall
						} // for d
						room_start = hall_end_pos;
						wall_pos  += sh_spacing + sh_width;
						min_eq(wall_pos, p->d[!min_dim][1]); // clamp to far side of building to handle the final row of rooms
					} // for i
				} // end secondary hallways
			} // end multiple hallways case

			else { // single main hallway
				room_walls.reserve(2*(num_rooms-1));
				hall_walls.reserve(2*(num_rooms+1));
				cube_t rwall(*p); // copy from part; shared zvals, but X/Y will be overwritten per wall
				create_wall(rwall, !min_dim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing); // room walls, create first wall
				vector<float> doorway_vals(2*num_rooms, 0.0);

				for (int i = 0; i+1 < num_rooms; ++i) { // num_rooms-1 walls
					for (unsigned d = 0; d < 2; ++d) {
						room_walls.push_back(rwall);
						room_walls.back().d[min_dim][!d] = hall_wall_pos[d];
						cube_t hwall(room_walls.back());
						set_wall_width(hwall, hall_wall_pos[d], wall_half_thick, min_dim);
						for (unsigned e = 0; e < 2; ++e) {hwall.d[!min_dim][e] += (e ? 1.0f : -1.0f)*hwall_extend;}
						if (partial_room && i+2 == num_rooms) {hwall.d[!min_dim][1] -= 1.5*doorway_width;} // pull back a bit to make room for a doorway at the end of the hall
						hall_walls.push_back(hwall); // longer sections that form T-junctions with room walls
						doorway_vals[2*i+d+1] = hwall.d[!min_dim][d];
					}
					rwall.translate_dim(!min_dim, room_len);
				} // for i
				for (unsigned s = 0; s < 2; ++s) { // add half length hall walls at each end of the hallway
					cube_t hwall(rwall); // copy to get correct zvals
					float const hwall_len((partial_room && s == 1) ? doorway_width : hwall_extend); // hwall for partial room at end is only length doorway_width
					hwall.d[!min_dim][ s] = place_area.d[!min_dim][s]; // end at the wall
					hwall.d[!min_dim][!s] = hwall.d[!min_dim][s] + (s ? -1.0f : 1.0f)*hwall_len; // end at first doorway
					doorway_vals[s*(2*num_rooms-1)] = hwall.d[!min_dim][!s];

					for (unsigned d = 0; d < 2; ++d) {
						set_wall_width(hwall, hall_wall_pos[d], wall_half_thick, min_dim);
						hall_walls.push_back(hwall);
					}
				} // for s
				// add rooms and doors
				interior->rooms.reserve(2*num_rooms + 1); // two rows of rooms + hallway
				interior->doors.reserve(2*num_rooms);
				float const wall_end(p->d[!min_dim][1]);
				float pos(p->d[!min_dim][0]);

				for (int i = 0; i < num_rooms; ++i) {
					// Note: it would probably be simpler to create one wall and cut doorways into it like what's done in the secondary hallways case above
					float const next_pos(min(wall_end, (pos + room_len))); // clamp to end of building to last row handle partial room)

					for (unsigned d = 0; d < 2; ++d) { // left, right sides of hallway
						cube_t c(*p); // copy zvals and exterior wall pos
						c.d[ min_dim][!d] = hall_wall_pos[d];
						c.d[!min_dim][ 0] = pos;
						c.d[!min_dim][ 1] = next_pos;
						add_room(c, part_id, 1, 0, 1); // office
						if (i == num_rooms/2) {interior->rooms.back().assign_all_to(RTYPE_BATH);} // assign the middle room to be a bathroom
						door_t door(c, min_dim, d); // copy zvals and wall pos
						clip_wall_to_ceil_floor(door, fc_thick);
						door.d[ min_dim][d] = hall_wall_pos[d]; // set to zero area at hallway
						for (unsigned e = 0; e < 2; ++e) {door.d[!min_dim][e] = doorway_vals[2*i+e];}
						add_interior_door(door);
					}
					pos = next_pos;
				} // for i
			} // end single main hallway case
			add_room(hall, part_id, 3, 1, 0); // add hallway as room with 3 lights
			if (p->z1() == ground_floor_z1 || pri_hall.is_all_zeros()) {pri_hall = hall;} // assign to primary hallway if on first floor of hasn't yet been assigned
			for (unsigned d = 0; d < 2; ++d) {first_wall_to_split[d] = interior->walls[d].size();} // don't split any walls added up to this point
		} // end use_hallway

		else { // generate random walls using recursive 2D slices
			float const min_wall_len2(0.85*min_wall_len); // a somewhat shorter value that applies to some tests (but not wall_split_thresh)
			bool const no_walls(min(p->dx(), p->dy()) < min_wall_len2); // not enough space to add a room (chimney, porch support, garage, shed, etc.)
			float const min_split_len(max(global_building_params.wall_split_thresh, 1.0f)*min_wall_len);
			assert(to_split.empty());
			if (no_walls) {add_room(*p, part_id, 1, 0, 0);} // add entire part as a room
			else {to_split.emplace_back(*p);} // seed room is entire part, no door
			bool is_first_split(1);
			point part_door_open_dir_tp(p->get_cube_center()); // used to determine in which direction doors open; updated base on central hallway
			
			if (first_part) { // reserve walls/rooms/doors - take a guess at the correct size
				for (unsigned d = 0; d < 2; ++d) {interior->walls[d].reserve(8*parts.size());}
				interior->rooms.reserve(8*parts.size()); // two rows of rooms + optional hallway
				interior->doors.reserve(4*parts.size());
			}
			while (!to_split.empty()) {
				split_cube_t const c(to_split.back()); // Note: non-const because door_lo/door_hi is modified during T-junction insert
				to_split.pop_back();
				vector3d const csz(c.get_size());
				bool wall_dim(0); // which dim the room is split by
				if      (csz.y > min_wall_len2 && csz.x > 1.25*csz.y) {wall_dim = 0;} // split long room in x
				else if (csz.x > min_wall_len2 && csz.y > 1.25*csz.x) {wall_dim = 1;} // split long room in y
				else {wall_dim = rgen.rand_bool();} // choose a random split dim for nearly square rooms
				
				//if (csz[!wall_dim] < min_wall_len || csz[wall_dim] < 1.5*doorway_width) {
				if (min(csz.x, csz.y) < min_wall_len2) {
					add_room(c, part_id, 1, 0, 0);
					continue; // not enough space to add a wall
				}
				float const min_dist_abs(1.5*doorway_width + wall_thick);
				bool const on_edge(c.d[!wall_dim][0] == p->d[!wall_dim][0] || c.d[!wall_dim][1] == p->d[!wall_dim][1]); // at edge of the building - make sure walls don't intersect windows
				float wall_pos(0.0);
				bool pos_valid(0);

				for (unsigned num = 0; num < 20; ++num) { // 20 tries to choose a wall pos that's not inside a window
					wall_pos = cube_rand_side_pos(c, wall_dim, 0.25, min_dist_abs, rgen);
					if (on_edge && is_val_inside_window(*p, wall_dim, wall_pos, window_hspacing[wall_dim], window_border)) continue; // try a new wall_pos
					if (c.bad_pos(wall_pos, wall_dim)) continue; // intersects doorway from prev wall, try a new wall_pos
					pos_valid = 1; break; // done, keep wall_pos
				}
				if (!pos_valid) { // no valid pos, skip this split
					add_room(c, part_id, 1, 0, 0);
					continue;
				}
				cube_t wall(c), wall2; // copy from cube; shared zvals, but X/Y will be overwritten per wall
				create_wall(wall, wall_dim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing);
				float const doorway_pos(cube_rand_side_pos(c, !wall_dim, 0.25, doorway_width, rgen));
				float const lo_pos(doorway_pos - doorway_hwidth), hi_pos(doorway_pos + doorway_hwidth);
				bool const open_dir(wall_pos > part_door_open_dir_tp[wall_dim]); // doors open away from the building center
				remove_section_from_cube_and_add_door(wall, wall2, lo_pos, hi_pos, !wall_dim, open_dir);
				interior->walls[wall_dim].push_back(wall);
				interior->walls[wall_dim].push_back(wall2);
				float door_lo[2] = {lo_pos, lo_pos}, door_hi[2] = {hi_pos, hi_pos}; // passed to next split step to avoid placing a wall that intersects this doorway
				bool const do_split(csz[wall_dim] > min_split_len); // split into two smaller rooms

				if (is_house && is_first_split && csz[!wall_dim] > 1.2*min_split_len) { // wall/hall len is at least enough to place 2-3 rooms
					// maybe create a hallway: create another split parallel to this one offset a bit and make the room in between a hallway
					bool const dir((wall_pos - c.d[wall_dim][0]) < (c.d[wall_dim][1] - wall_pos)); // further part edge
					float other_wall_pos(0.0);
					bool other_pos_valid(0);

					for (unsigned num = 0; num < 10; ++num) { // 10 tries to choose a wall pos that's not inside a window
						other_wall_pos = wall_pos + ((dir ? 1.0 : -1.0)*rgen.rand_uniform(1.6, 2.4)*doorway_width); // opposite edge of hallway
						if (on_edge && is_val_inside_window(*p, wall_dim, other_wall_pos, window_hspacing[wall_dim], window_border)) continue; // try a new wall_pos
						other_pos_valid = 1; break; // done, keep wall_pos
					}
					float const first_side_room_len(dir ? (wall_pos - c.d[wall_dim][0]) : (c.d[wall_dim][1] - wall_pos));
					float const other_side_room_len(dir ? (c.d[wall_dim][1] - other_wall_pos) : (other_wall_pos - c.d[wall_dim][0]));

					if (other_pos_valid && min(first_side_room_len, other_side_room_len) > 1.2*min_wall_len) { // rooms on both sides of the hallway are large enough
						// create hallway room
						cube_t hall(c);
						hall.d[wall_dim][ dir] = other_wall_pos; // Note: building navigation code wants hallways to end at wall centerlines
						hall.d[wall_dim][!dir] = wall_pos;
						assert(hall.is_strictly_normalized());
						unsigned const num_lights(2.0*csz[!wall_dim]/min_split_len); // more lights for longer hallways
						add_room(hall, part_id, num_lights, 1, 0);
						// add other wall parts and doorway, with a different random doorway pos
						cube_t o_wall1(c), o_wall2;
						create_wall(o_wall1, wall_dim, other_wall_pos, fc_thick, wall_half_thick, wall_edge_spacing);
						float const door_pos(cube_rand_side_pos(c, !wall_dim, 0.25, doorway_width, rgen));
						float const o_lo_pos(door_pos - doorway_hwidth), o_hi_pos(door_pos + doorway_hwidth);
						remove_section_from_cube_and_add_door(o_wall1, o_wall2, o_lo_pos, o_hi_pos, !wall_dim, !open_dir); // opens in other dir
						interior->walls[wall_dim].push_back(o_wall1);
						interior->walls[wall_dim].push_back(o_wall2);
						wall.d[wall_dim][dir] = o_wall1.d[wall_dim][dir]; // make original wall the entire width of the hall so that the room on the other side doesn't overlap the hall
						door_lo[dir] = o_lo_pos;
						door_hi[dir] = o_hi_pos;
						part_door_open_dir_tp[wall_dim] = 0.5f*(wall_pos + other_wall_pos); // use the center of the hallway so that doors open into the rooms
						if (get_real_num_parts() == 1) {bldg_door_open_dir_tp[wall_dim] = part_door_open_dir_tp[wall_dim];} // this is also the test point for the entire building
					}
				}
				for (unsigned d = 0; d < 2; ++d) { // still have space to split in other dim, add the two parts to the stack
					split_cube_t c_sub(c);
					c_sub.d[wall_dim][d] = wall.d[wall_dim][!d]; // clip to wall pos
					c_sub.door_lo[!wall_dim][d] = door_lo[!d] - wall_half_thick; // set new door pos in this dim (keep door pos in other dim, if set)
					c_sub.door_hi[!wall_dim][d] = door_hi[!d] + wall_half_thick;
					if (do_split) {to_split.push_back(c_sub);} else {add_room(c_sub, part_id, 1, 0, 0);} // leaf case (unsplit), add a new room
				}
				is_first_split = 0;
			} // end while()
			// insert walls to split up parts into rectangular rooms
			float const min_split_wall_len(0.5*min_wall_len); // allow a shorter than normal wall because these walls have higher priority
			
			if (min(p->dx(), p->dy()) > min_split_wall_len) { // if not too small
				for (auto p2 = parts.begin(); p2 != get_real_parts_end(); ++p2) {
					if (p2 == p) continue; // skip self
					if (min(p2->dx(), p2->dy()) < min_split_wall_len) continue; // too small, skip

					for (unsigned dim = 0; dim < 2; ++dim) {
						for (unsigned dir = 0; dir < 2; ++dir) {
							float const val(p->d[!dim][dir]);
							if (p2->d[!dim][!dir] != val) continue; // not adjacent
							if (p2->z1() >= p->z2() || p2->z2() <= p->z1()) continue; // no overlap in Z
							if (p2->d[dim][0] > p->d[dim][0] || p2->d[dim][1] < p->d[dim][1]) continue; // not contained in dim (don't have to worry about Z-shaped case)

							if (p2->d[dim][0] == p->d[dim][0] && p2->d[dim][1] == p->d[dim][1]) { // same xy values, must only vary in z
								if (p2->z2() < p->z2()) continue; // add wall only on one side (arbitrary)
							}
							cube_t wall;
							wall.z1() = max(p->z1(), p2->z1()) + fc_thick; // shared Z range
							wall.z2() = min(p->z2(), p2->z2()) - fc_thick;

							for (unsigned d = 0; d < 2; ++d) {
								wall.d[dim][d] = place_area.d[dim][d]; // shorter part side with slight offset
								if (p->d[dim][d] != p2->d[dim][d]) {wall.d[dim][d] += (d ? 1.0 : -1.0)*0.8*wall_edge_spacing;} // reduce the gap at the corner between the two parts
							}
							if (wall.get_sz_dim(dim) < min_split_wall_len) continue; // wall is too short to add (can this happen?)
							wall.d[!dim][ dir] = val;
							wall.d[!dim][!dir] = val + (dir ? -1.0 : 1.0)*wall_thick;
							must_split[!dim] |= (1ULL << (interior->walls[!dim].size() & 63)); // flag this wall for extra splitting
							interior->walls[!dim].push_back(wall);
							assert(rooms_start < interior->rooms.size()); // must have added at least one room

							for (auto r = (interior->rooms.begin() + rooms_start); r != interior->rooms.end(); ++r) {
								assert(p->contains_cube(*r));
								if (r->d[!dim][dir] == val) {r->d[!dim][dir] = wall.d[!dim][!dir];} // clip room to exclude this newly added wall
							}
						} // for dir
					} // for dim
				} // for p2
			} // end !too_small
		} // end wall placement
		add_ceilings_floors_stairs(rgen, *p, hall, (p - parts.begin()), num_floors, rooms_start, use_hallway, first_part_this_stack, window_hspacing, window_border);
	} // for p (parts)

	if (has_sec_bldg()) { // add garage/shed floor and ceiling
		assert(parts_end < parts.end());
		cube_t const &garage(*parts_end);
		cube_t C(garage);
		C.z2() = C.z1() + fc_thick;
		interior->floors.push_back(C);
		C.z2() = garage.z2();
		C.z1() = C.z2() - fc_thick;
		interior->ceilings.push_back(C);
		add_room(garage, (parts_end - parts.begin()), 1, 0, 0, 1); // is_sec_bldg=1
		interior->rooms.back().no_geom = 1;
	}
	// attempt to cut extra doorways into long walls if there's space to produce a more connected floorplan
	for (unsigned d = 0; d < 2; ++d) { // x,y: dim in which the wall partitions the room (wall runs in dim !d)
		vect_cube_t &walls(interior->walls[d]);
		vect_cube_t const &perp_walls(interior->walls[!d]);

		// Note: iteration will include newly added all segments to recursively split long walls
		for (unsigned w = first_wall_to_split[d]; w < walls.size(); ++w) { // skip hallway walls
			bool pref_split(must_split[d] & (1ULL << (w & 63)));

			for (unsigned nsplits = 0; nsplits < 4; ++nsplits) { // at most 4 splits
				cube_t &wall(walls[w]); // take a reference here because a prev iteration push_back() may have invalidated it
				float const len(wall.get_sz_dim(!d)), min_split_len((pref_split ? 0.5 : 1.5)*min_wall_len); // = 2.0/6.0 * doorway_width
				if (len < min_split_len) break; // not long enough to split - done
				float const min_dist_abs(min(1.5f*doorway_width, 0.5f*min_split_len));
				// walls currently don't run along the inside of exterior building walls, so we don't need to handle that case yet
				bool was_split(0);

				for (unsigned ntries = 0; ntries < (pref_split ? 40U : 10U); ++ntries) { // choose random doorway positions and check against perp_walls for occlusion
					float const doorway_pos(cube_rand_side_pos(wall, !d, 0.2, min_dist_abs, rgen));
					float const lo_pos(doorway_pos - doorway_hwidth), hi_pos(doorway_pos + doorway_hwidth);
					bool valid(1);

					for (auto p = (perp_walls.begin() + first_wall_to_split[!d]); p != perp_walls.end(); ++p) { // skip hallway walls
						if (p->z1() != wall.z1()) continue; // not the same zval/story
						if (p->d[!d][1] < lo_pos-wall_thick || p->d[!d][0] > hi_pos+wall_thick) continue; // no overlap with wall
						if (p->d[ d][1] > wall.d[d][0]-wall_thick && p->d[d][0] < wall.d[d][1]+wall_thick) {valid = 0; break;} // has perp intersection
					}
					if (valid && !pref_split) { // don't split walls into small segments that border the same two rooms on both sides (two doorways between the same pair of rooms)
						float const lo[2] = {wall.d[!d][0]-wall_thick, lo_pos}, hi[2] = {hi_pos, wall.d[!d][1]+wall_thick}; // ranges of the two split wall segments, grown a bit

						for (unsigned s = 0; s < 2; ++s) { // check both wall segments
							bool contained[2] = {0,0};

							for (unsigned e = 0; e < 2; ++e) { // check both directions from the wall
								for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
									if (r->z1() >= wall.z2() || r->z2() <= wall.z1()) continue; // no overlap in Z
									// skip wall edges co-incident with room edges where the entire wall is contained in the room; this can happen at part boundaries
									if (wall.d[d][e] == r->d[d][e] && wall.d[d][!e] < r->d[d][1] && wall.d[d][!e] > r->d[d][0]) continue;
									if (wall.d[d][e] < r->d[d][0] || wall.d[d][e] > r->d[d][1]) continue; // wall not inside room in dim d/dir e
									if (lo[s] > r->d[!d][0] && hi[s] < r->d[!d][1]) {contained[e] = 1; break;} // entire wall contained in span of room
								}
							}
							if (contained[0] && contained[1]) {valid = 0; break;} // wall seg contained in rooms on both sides => two doors in same wall between rooms => drop
						} // for s
					}
					if (!valid) continue;
					cube_t cand(wall); // sub-section of wall that will become a doorway
					cand.d[!d][0] = lo_pos; cand.d[!d][1] = hi_pos;
					bool const elevators_only(pref_split && ntries > 20); // allow blocking stairs if there's no other way to insert a door
					if (interior->is_blocked_by_stairs_or_elevator(cand, doorway_width, elevators_only)) continue; // stairs in the way, skip; should we assert !pref_split?
					bool const open_dir(wall.get_center_dim(d) > bldg_door_open_dir_tp[d]); // doors open away from the building center
					insert_door_in_wall_and_add_seg(wall, lo_pos, hi_pos, !d, open_dir, 0); // Note: modifies wall
					was_split = 1;
					break;
				} // for ntries
				if (!was_split) break; // no more splits
				pref_split = 0; // already split, no longer preferred
			} // for nsplits
		} // for w
	} // for d

	// add stairs to connect together stacked parts for office buildings; must be done last after all walls/ceilings/floors have been assigned
	for (auto p = parts.begin(); p != parts_end; ++p) {connect_stacked_parts_with_stairs(rgen, *p);}
}

void building_t::add_ceilings_floors_stairs(rand_gen_t &rgen, cube_t const &part, cube_t const &hall, unsigned part_ix, unsigned num_floors,
	unsigned rooms_start, bool use_hallway, bool first_part_this_stack, float window_hspacing[2], float window_border)
{
	// increase floor thickness if !is_house? but then we would probably have to increase the space between floors as well, which involves changing the texture scale
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const doorway_width(0.5*window_vspacing), wall_thickness(get_wall_thickness());
	float ewidth(1.5*doorway_width); // for elevators
	float z(part.z1());
	cube_t stairs_cut, elevator_cut;
	bool stairs_dim(0), add_elevator(0), stairs_have_railing(1), stairs_against_wall[2] = {0, 0};
	stairs_shape sshape(SHAPE_STRAIGHT); // straight by default
	bool const must_add_stairs(first_part_this_stack || (has_complex_floorplan && part == parts.back())); // first part in stack, or tallest/last part of complex building

	// add stairwells and elevator shafts
	if (num_floors == 1) {} // no need for stairs or elevator
	else if (use_hallway) { // part is the hallway cube
		add_elevator = 1;
		if (interior->landings.empty()) {interior->landings.reserve(add_elevator ? 1 : (num_floors-1));} // lower bound
		assert(!interior->rooms.empty());
		room_t &room(interior->rooms.back()); // hallway is always the last room to be added
		bool const long_dim(hall.dx() < hall.dy());
		// U-shape if there's enough room
		if (room.get_sz_dim(!long_dim) > 6.0*doorway_width) {sshape = SHAPE_U; ewidth *= 1.6;} // increase the width of both the stairs and elevator
		else {sshape = SHAPE_WALLED_SIDES;} // walled sides to meet fire codes
		cube_t stairs(hall); // start as hallway

		if (add_elevator) {
			point center(room.get_cube_center());
			float const center_shift(0.125*room.get_sz_dim(long_dim)*(rgen.rand_bool() ? -1.0 : 1.0));
			center[long_dim] += center_shift; // make elevator off-center
			elevator_t elevator(room, long_dim, rgen.rand_bool(), rgen.rand_bool(), 0); // elevator shaft
			elevator.x1() = center.x - 0.5*ewidth; elevator.x2() = center.x + 0.5*ewidth;
			elevator.y1() = center.y - 0.5*ewidth; elevator.y2() = center.y + 0.5*ewidth;
			add_or_extend_elevator(elevator, 1);
			room.has_elevator = 1;
			elevator_cut      = elevator;
			stairs.translate_dim(long_dim, -center_shift); // shift stairs in the opposite direction
		}
		// always add stairs
		for (unsigned dim = 0; dim < 2; ++dim) { // shrink in XY
			bool const is_step_dim(bool(dim) == long_dim); // same orientation as the hallway
			float shrink(stairs.get_sz_dim(dim) - (is_step_dim ? 4.0*doorway_width : 0.9*ewidth)); // set max size of stairs opening
			stairs.expand_in_dim(dim, -0.5*shrink); // centered in the hallway
		}
		room.has_stairs = 1;
		stairs_cut      = stairs;
		stairs_dim      = long_dim;
	}
	else if (!is_house || interior->stairwells.empty()) { // only add stairs to first part of a house unless we haven't added stairs yet
		// sometimes add an elevator to building parts, but not the first part in a stack (to guarantee we have at least one set of stairs)
		// it might not be possible to place an elevator a part with no interior rooms, but that should be okay, because some other part will still have stairs
		// do we need support for multiple floor cutouts stairs + elevator in this case as well?
		if (!is_house && !must_add_stairs) {
			float const elevator_prob[4] = {0.9, 0.5, 0.3, 0.1}; // higher chance of adding an elevator if there are more existing elevators
			add_elevator = (rgen.rand_float() < elevator_prob[min(interior->elevators.size(), size_t(3))]);
		}
		unsigned const rooms_end(interior->rooms.size()), num_avail_rooms(rooms_end - rooms_start);
		assert(num_avail_rooms > 0); // must have added at least one room
		float stairs_scale(1.0);

		for (unsigned N = 0; N < 4; ++N) {
			unsigned const rand_ix(rgen.rand()); // choose a random starting room to make into a stairwell

			for (unsigned n = 0; n < num_avail_rooms; ++n) { // try all available rooms starting with the selected one to see if we can fit a stairwell/elevator in any of them
				unsigned const stairs_room(rooms_start + (rand_ix + n)%num_avail_rooms);
				room_t &room(interior->rooms[stairs_room]);
				assert(room.part_id == part_ix); // sanity check

				if (add_elevator) {
					if (min(room.dx(), room.dy()) < 2.0*ewidth) continue; // room is too small to place an elevator
					bool placed(0);

					for (unsigned y = 0; y < 2 && !placed; ++y) { // try all 4 corners
						for (unsigned x = 0; x < 2 && !placed; ++x) {
							int const wtype_x(classify_room_wall(room, room.z1(), 0, x, 1)), wtype_y(classify_room_wall(room, room.z1(), 1, y, 1)); // include partial sep walls
							if (wtype_x == ROOM_WALL_SEP || wtype_y == ROOM_WALL_SEP) continue; // don't place elevators between parts where they could block doorways
							float const xval(room.d[0][x] + (x ? -ewidth : ewidth)), yval(room.d[1][y] + (y ? -ewidth : ewidth)), shrink(0.01*ewidth);
							// check room interior edge for intersection with windows
							if (wtype_x == ROOM_WALL_EXT && is_val_inside_window(part, 0, xval, window_hspacing[0], window_border)) continue;
							if (wtype_y == ROOM_WALL_EXT && is_val_inside_window(part, 1, yval, window_hspacing[1], window_border)) continue;
							bool const dim(rgen.rand_bool()), is_open(rgen.rand_bool());
							elevator_t elevator(room, dim, !(dim ? y : x), is_open, (wtype_x == ROOM_WALL_EXT || wtype_y == ROOM_WALL_EXT)); // elevator shaft
							elevator.d[0][!x] = xval;
							elevator.d[1][!y] = yval;
							// shrink to leave a small gap between the outer wall to prevent z-fighting
							if (wtype_x == ROOM_WALL_EXT) {elevator.d[0][x] += (x ? -shrink : shrink);}
							if (wtype_y == ROOM_WALL_EXT) {elevator.d[1][y] += (y ? -shrink : shrink);}
							if (has_bcube_int(elevator, interior->exclusion)) continue; // try again
							if (is_cube_close_to_doorway(elevator, room))     continue; // try again
							add_or_extend_elevator(elevator, 1);
							elevator_cut = elevator;
							placed       = 1; // successfully placed
						} // for x
					} // for y
					if (!placed) continue; // try another room
					room.has_elevator = 1;
					room.no_geom      = 1;
				}
				else { // stairs
					cube_t cutout(room);
					cutout.expand_by_xy(-floor_thickness); // padding around walls
					float const dx(cutout.dx()), dy(cutout.dy()); // choose longer dim if high aspect ratio
					float const stairs_sz(stairs_scale*doorway_width);
					if      (dx > 1.2*dy) {stairs_dim = 0;}
					else if (dy > 1.2*dx) {stairs_dim = 1;}
					else {stairs_dim = rgen.rand_bool();} // close to square
					if (cutout.get_sz_dim(stairs_dim) < 4.0*stairs_sz || cutout.get_sz_dim(!stairs_dim) < 3.0*stairs_sz) continue; // not enough space for stairs

					for (unsigned dim = 0; dim < 2; ++dim) { // shrink in XY
						bool const is_step_dim(bool(dim) == stairs_dim);
						float shrink(cutout.get_sz_dim(dim) - (is_step_dim ? 4.0 : 1.2)*stairs_sz); // set max size of stairs opening
						max_eq(shrink, 2.0f*stairs_sz); // allow space for doors to open and player to enter/exit
						cutout.expand_in_dim(dim, -0.5*shrink); // centered in the room

						if (!is_step_dim) { // see if we can push the stairs to the wall on one of the sides without blocking a doorway
							bool const first_dir(rgen.rand_bool());

							for (unsigned d = 0; d < 2; ++d) {
								bool const dir(bool(d) ^ first_dir);
								// if the room is on the edge of the part that's not on the building exterior, then this room connects two parts and we need to place a door here later
								if (classify_room_wall(room, room.z1(), dim, dir, 1) == ROOM_WALL_SEP) continue; // include partial sep walls
								cube_t cand(cutout);
								// add small gap to prevent z-fighting and FP accuracy asserts
								float const shift((cand.d[dim][dir] - room.d[dim][dir]) - (dir ? -1.0 : 1.0)*wall_thickness); // negative if dir==1
								cand.d[dim][0] -= shift; cand.d[dim][1] -= shift; // close the gap - flush with the wall
								if (!is_cube_close_to_doorway(cand, room)) {cutout = cand; stairs_against_wall[dir] = 1; break;} // keep if it's good
							} // for d
						}
					} // for dim
					if (!is_house && (stairs_against_wall[0] || stairs_against_wall[1])) {sshape = SHAPE_WALLED_SIDES;} // add wall between room and office stairs if against a room wall
					if (interior->landings.empty()) {interior->landings.reserve(num_floors-1);}
					assert(cutout.is_strictly_normalized());
					stairs_cut      = cutout;
					room.has_stairs = 1;
					if (use_hallway || !pri_hall.is_all_zeros()) {room.no_geom = 1;} // no geom in an office with stairs for buildings with hallways
				}
				break; // success - done
			} // for n
			if (!must_add_stairs || add_elevator || !stairs_cut.is_all_zeros()) break; // successfully placed stairs, or not required to place stairs
			stairs_scale -= 0.1; // shrink stairs a bit and try again
		} // for N
	} // end stairs/elevator placement

	// add ceilings and floors; we have num_floors+1 separators; the first is only a floor, and the last is only a ceiling
	cube_t C(part);
	C.z1() = z; C.z2() = z + fc_thick;
	unsigned const floors_start(interior->floors.size());
	interior->floors.push_back(C); // ground floor, full area
	z += window_vspacing; // move to next floor
	bool const has_stairs(!stairs_cut.is_all_zeros()), has_elevator(!elevator_cut.is_all_zeros());
	bool const stairs_dir(has_stairs ? rgen.rand_bool() : 0); // same for every floor, could maybe alternate for stairwells
	cube_t &first_cut(has_elevator ? elevator_cut : stairs_cut); // elevator is larger
	unsigned last_landing_ix(0);

	for (unsigned f = 1; f < num_floors; ++f, z += window_vspacing) { // skip first floor - draw pairs of floors and ceilings
		cube_t to_add[8]; // up to 2 cuts for stairs + elevator
		float const zc(z - fc_thick), zf(z + fc_thick);

		if (!has_stairs && !has_elevator) {to_add[0] = part;} // neither - add single cube
		else {
			bool const is_at_top(f+1 == num_floors);
			subtract_cube_xy(part, first_cut, to_add);

			if (has_stairs && has_elevator) { // both
				bool found(0);

				for (unsigned n = 0; n < 4; ++n) { // find the cube where the stairs are placed
					if (!to_add[n].intersects_xy_no_adj(stairs_cut)) continue;
					subtract_cube_xy(to_add[n], stairs_cut, to_add+4); // append up to 4 more cubes
					to_add[n].set_to_zeros(); // this cube was replaced
					found = 1;
					break; // assume there is only one; it's up to the placer step to ensure this
				}
				assert(found);
			}
			if (has_stairs) { // add landings and stairwells
				// make sure to enable back wall for the first flight of stairs
 				landing_t landing(stairs_cut, 0, f, stairs_dim, stairs_dir, stairs_have_railing, ((f == 1 && sshape == SHAPE_WALLED_SIDES) ? (stairs_shape)SHAPE_WALLED : sshape), 0, is_at_top);
				landing.z1() = zc; landing.z2() = zf;
				landing.set_against_wall(stairs_against_wall);
				last_landing_ix = interior->landings.size();
				interior->landings.push_back(landing);

				if (f == 1) { // only add for first floor
					interior->stairwells.emplace_back(stairs_cut, num_floors, stairs_dim, stairs_dir, sshape);
					interior->stairwells.back().set_against_wall(stairs_against_wall);
				}
			}
			if (has_elevator) {
				assert(!interior->elevators.empty());
				landing_t landing(elevator_cut, 1, f, interior->elevators.back().dim, interior->elevators.back().dir);
				landing.z1() = zc; landing.z2() = zf;
				interior->landings.push_back(landing);
			}
		}
		for (unsigned i = 0; i < 8; ++i) { // skip zero area cubes from stairs/elevator shafts along an exterior wall
			cube_t &c(to_add[i]);
			if (c.is_zero_area()) continue;
			c.z1() = zc; c.z2() = z;  interior->ceilings.push_back(c);
			c.z1() = z;  c.z2() = zf; interior->floors  .push_back(c);
			//c.z1() = zf; c.z2() = zc + window_vspacing; // add per-floor walls, door cutouts, etc. here
		}
	} // for f
	bool has_roof_access(0);

	if (must_add_stairs && has_stairs && !is_house && roof_type == ROOF_TYPE_FLAT && !has_helipad) { // add roof access for stairs
		bool const is_sloped(sshape != SHAPE_U);
		cube_t box(stairs_cut);
		if (!is_sloped) {box.expand_by_xy(fc_thick);}
		box.z1() = z + floor_thickness; box.z2() = z + window_vspacing;
		box.z2() -= (is_sloped ? 0.15 : 0.2)*window_vspacing; // slightly lower than a normal floor
		cube_t check_box(box);
		check_box.d[stairs_dim][stairs_dir] += (stairs_dir ? 1.0 : -1.0)*doorway_width; // expand at stairs exit to ensure clearance

		if (!has_bcube_int_no_adj(check_box, parts)) { // no overlap with other parts (should we check in front?)
			float const zc(z - fc_thick);
			cube_t to_add[4]; // only one cut / 4 cubes (-y, +y, -x, +x)
			subtract_cube_xy(part, stairs_cut, to_add);
			interior->landings[last_landing_ix].is_at_top = 0; // previous landing is no longer at the top
			landing_t landing(stairs_cut, 0, num_floors, stairs_dim, stairs_dir, 1, sshape, 1, 1); // stairs_have_railing=1, roof_access=1, is_at_top=1
			landing.z1() = zc; landing.z2() = z; // no floor above
			landing.set_against_wall(stairs_against_wall);
			interior->landings.push_back(landing);
			interior->stairwells.back().z2() += fc_thick; // extend upward
			interior->stairwells.back().z1() += fc_thick; // requiured to trick roof clipping into treating this as a stack connector stairwell

			for (unsigned i = 0; i < 4; ++i) { // skip zero area cubes from stairs/elevator shafts along an exterior wall
				cube_t &c(to_add[i]);
				if (c.is_zero_area()) continue;
				c.z1() = zc; c.z2() = z; interior->ceilings.push_back(c);
				c.set_to_zeros();
			}
			bool const dir(stairs_dir^(sshape == SHAPE_U));
			box.z1() = z;
			// add a door to the roof
			float const door_shift(0.2*(dir ? -1.0 : 1.0)*fc_thick);
			cube_t door(box);
			door.d[stairs_dim][ dir] += door_shift; // shift slightly to fill the gap
			door.d[stairs_dim][!dir]  = door.d[stairs_dim][dir];
			if (!is_sloped) {door.d[!stairs_dim][1] = door.get_center_dim(!stairs_dim);} // only the open half
			add_door(door, part_ix, stairs_dim, dir, 1, 1); // roof_access=1
			// clear any roof objects that are in the way
			cube_t clear_cube(box);
			clear_cube.d[stairs_dim][dir] += (dir ? 1.0 : -1.0)*window_vspacing; // clear out space in front of the door
			remove_intersecting_roof_cubes(clear_cube);
			// add a small 3-sided box around the stairs using roof blocks
			unsigned const opening_ix(2*(1 - stairs_dim) + dir);

			if (is_sloped) { // sloped roof
				float const z1(box.z1()), z2(box.z2());
				point const pts[4] = {point(box.x1(), box.y1(), z1), point(box.x2(), box.y1(), z1), point(box.x2(), box.y2(), z1), point(box.x1(), box.y2(), z1)};
				unsigned const hi_side1[4] = {0, 2, 3, 1}, hi_side2[4] = {1, 3, 0, 2}, lo_side1[4] = {3, 1, 2, 0}, lo_side2[4] = {2, 0, 1, 3};
				unsigned const top_vix[2] = {hi_side1[opening_ix], hi_side2[opening_ix]}, bot_vix[2] = {lo_side1[opening_ix], lo_side2[opening_ix]};
				tquad_t top(4); // quad
				for (unsigned n = 0; n < 4; ++n) {top.pts[n] = pts[n];}
				top.pts[top_vix[0]].z = top.pts[top_vix[1]].z = z2; // these are the higher points
				roof_tquads.emplace_back(top, (unsigned)tquad_with_ix_t::TYPE_ROOF_ACC);

				for (unsigned s = 0; s < 2; ++s) {
					tquad_t side(3); // triangle
					side.pts[0] = pts[bot_vix[s]]; // corner
					side.pts[1] = pts[top_vix[s]]; // bottom
					side.pts[2] = top.pts[top_vix[s]]; // top
					if (s) {swap(side.pts[0], side.pts[1]);} // make normal point outward for correct BFC
					roof_tquads.emplace_back(side, (unsigned)tquad_with_ix_t::TYPE_ROOF_ACC);

					tquad_t frame(4); // quads on sides of door frame
					frame.pts[0] = frame.pts[3] = pts[top_vix[s]]; // bottom
					frame.pts[1] = frame.pts[2] = top.pts[top_vix[s]] - vector3d(0.0, 0.0, 0.05*fc_thick); // top, shift slightly down to match the slope
					float const frame_width(0.08*((bool(s)^dir^stairs_dim^1) ? -1.0 : 1.0)*box.get_sz_dim(!stairs_dim));
					for (unsigned d = 0; d < 2; ++d) {frame.pts[(s<<1)+d][!stairs_dim] +=    frame_width;} // move toward center of door
					for (unsigned n = 0; n < 4; ++n) {frame.pts[n]       [ stairs_dim] += 0.6*door_shift;} // shift back behind the door
					roof_tquads.emplace_back(frame, (unsigned)tquad_with_ix_t::TYPE_ROOF_ACC);
				} // for s
			}
			else { // box roof
				cube_t hole(stairs_cut), front(box);
				hole.expand_by_xy(0.1*fc_thick); // to prevent z-fighting
				front.d[ stairs_dim][!dir] = hole.d[stairs_dim][dir];
				hole .d[ stairs_dim][ dir] = box. d[stairs_dim][dir]; // move edge flush with box to remove this wall and create an opening
				front.d[!stairs_dim][   0] = front.get_center_dim(!stairs_dim); // block off the non-opening half
				subtract_cube_xy(box, hole, to_add);

				for (unsigned i = 0; i < 4; ++i) {
					cube_t &c(to_add[i]);
					if (!c.is_zero_area()) {details.emplace_back(c, (uint8_t)ROOF_OBJ_SCAP);} // skip open side
				}
				box.z1() = front.z2() = box.z2() - fc_thick;
				details.emplace_back(box,   ROOF_OBJ_SCAP); // top
				details.emplace_back(front, ROOF_OBJ_SCAP); // front half
				max_eq(bcube.z2(), box.z2());
			}
			interior->stairwells.back().roof_access = 1;
			has_roof_access = 1;
		}
	}
	if (!has_roof_access) { // roof ceiling, full area
		interior->top_ceilings_mask |= (uint64_t(1) << (interior->ceilings.size() & 63)); // mark this as a top ceiling so that it can be drawn; okay if wraps around
		C.z1() = z - fc_thick; C.z2() = z;
		interior->ceilings.push_back(C);
	}
	std::reverse(interior->floors.begin()+floors_start, interior->floors.end()); // order floors top to bottom to reduce overdraw when viewed from above
}

bool building_t::check_cube_intersect_walls(cube_t const &c) const {
	for (unsigned d = 0; d < 2; ++d) {
		if (has_bcube_int(c, interior->walls[d])) return 1;
	}
	return 0;
}

bool building_t::is_valid_stairs_elevator_placement(cube_t const &c, float pad, bool check_walls) const {
	// check if any previously placed walls intersect this cand stairs/elevator; we really only need to check the walls from <part> and *p though
	if (interior->is_blocked_by_stairs_or_elevator(c, pad)) return 0;
	if (check_walls && check_cube_intersect_walls(c)) return 0;
	// if we're not checking walls, then at least check for open doors to avoid having the stairs intersect an open door
	if (is_cube_close_to_doorway(c, cube_t(), pad, !check_walls)) return 0;
	return 1;
}

void subtract_cube_from_cube(cube_t const &c, cube_t const &s, vect_cube_t &out) { // XY only
	cube_t C;
	if (c.y1() < s.y1()) {C = c; C.y2() = s.y1(); out.push_back(C);} // bottom
	if (c.y2() > s.y2()) {C = c; C.y1() = s.y2(); out.push_back(C);} // top
	if (c.x1() < s.x1()) {C = c; max_eq(C.y1(), s.y1()); min_eq(C.y2(), s.y2()); C.x2() = s.x1(); out.push_back(C);} // left center
	if (c.x2() > s.x2()) {C = c; max_eq(C.y1(), s.y1()); min_eq(C.y2(), s.y2()); C.x1() = s.x2(); out.push_back(C);} // right center
}
void subtract_cube_from_cube_inplace(cube_t const &s, vect_cube_t &cubes, unsigned &ix, unsigned &iter_end) { // Note: ix is an index to cubes
	unsigned const prev_sz(cubes.size());
	assert(ix < prev_sz);
	cube_t const c(cubes[ix]); // deep copy - reference will become invalid
	subtract_cube_from_cube(c, s, cubes);
	cubes[ix] = cubes.back(); cubes.pop_back(); // reuse this slot for one of the output cubes (or move the last cube here if there are no output cubes)
	if (cubes.size() <= prev_sz) {--ix; --iter_end;} // no cubes added, last cube was swapped into this slot and needs to be reprocessed
}
template<typename T> void subtract_cubes_from_cube(cube_t const &c, T const &sub, vect_cube_t &out, vect_cube_t &out2) { // XY only
	out.clear();
	out.push_back(c);

	for (auto s = sub.begin(); s != sub.end(); ++s) {
		if (s->z1() <= c.z1() || s->z1() >= c.z2() || s->z2() <= c.z2()) continue; // not correct floor
		if (!c.intersects_xy(c)) continue; // no overlap with orig cube (optimization)
		out2.clear();

		// clip all of out against *s, write results to out2, then swap with out
		for (auto i = out.begin(); i != out.end(); ++i) {
			if (!i->intersects_xy(*s)) {out2.push_back(*i); continue;} // no overlap, keep entire cube
			subtract_cube_from_cube(*i, *s, out2);
		}
		out.swap(out2);
	} // for s
}
template void subtract_cubes_from_cube(cube_t const &c, vect_cube_t const &sub, vect_cube_t &out, vect_cube_t &out2); // explicit instantiation

bool subtract_cube_from_cubes(cube_t const &s, vect_cube_t &cubes, vect_cube_t *holes, bool clip_in_z, bool include_adj) {
	unsigned iter_end(cubes.size()); // capture size before splitting
	bool was_clipped(0);

	for (unsigned i = 0; i < iter_end; ++i) {
		cube_t const &c(cubes[i]);
		if (!(include_adj ? c.intersects(s) : c.intersects_no_adj(s))) continue; // keep it
		
		if (holes) {
			cube_t hole(c);
			hole.intersect_with_cube(s);
			bool merged(0);

			// see if we can merge this hole with another hole to remove an interior face
			for (auto h = holes->begin(); h != holes->end() && !merged; ++h) {
				if (h->z1() != hole.z1() || h->z2() != hole.z2()) continue;

				for (unsigned d = 0; d < 2 && !merged; ++d) {
					if (h->d[!d][0] != hole.d[!d][0] || h->d[!d][1] != hole.d[!d][1]) continue;

					for (unsigned e = 0; e < 2; ++e) {
						if (h->d[d][e] == hole.d[d][!e]) {h->d[d][e] = hole.d[d][e]; merged = 1; break;}
					}
				}
			} // for h
			if (!merged) {holes->push_back(hole);}
		} // end holes handling
		if (clip_in_z) { // Note: remember that c reference will be invalidated below
			cube_t shared(c);
			shared.intersect_with_cube_xy(s);
			if (shared.z1() < s.z1()) {cube_t bot(shared); bot.z2() = s.z1(); cubes.push_back(bot);} // bottom part
			if (shared.z2() > s.z2()) {cube_t top(shared); top.z1() = s.z2(); cubes.push_back(top);} // top part
		}
		subtract_cube_from_cube_inplace(s, cubes, i, iter_end); // Note: invalidates c reference
		was_clipped = 1;
	} // for i
	return was_clipped;
}
template<typename T> void subtract_cubes_from_cubes(T const &sub, vect_cube_t &cubes) {
	for (auto i = sub.begin(); i != sub.end(); ++i) {subtract_cube_from_cubes(*i, cubes);}
}

void subtract_cube_from_floor_ceil(cube_t const &c, vect_cube_t &fs) {
	unsigned iter_end(fs.size()); // capture orig size

	for (unsigned i = 0; i < iter_end; ++i) {
		cube_t const &cur(fs[i]);
		if (cur.z1() > c.z2() || cur.z2() < c.z1()) continue; // no z overlap
		if (cur.intersects_no_adj(c)) {subtract_cube_from_cube_inplace(c, fs, i, iter_end);} // Note: invalidates cur reference
	}
}

void building_t::connect_stacked_parts_with_stairs(rand_gen_t &rgen, cube_t const &part) { // and extend elevators vertically; part is on the bottom

	//highres_timer_t timer("Connect Stairs"); // 72ms (serial)
	float const window_vspacing(get_window_vspace()), fc_thick(0.5*get_floor_thickness()), wall_thickness(get_wall_thickness());
	float const doorway_width(0.5*window_vspacing), stairs_len(4.0*doorway_width);
	bool const is_basement(has_basement() && part == parts[basement_part_ix]);

	if (part.z2() < bcube.z2()) { // if this is the top floor, there is nothing above it (but roof geom may get us into this case anyway)
		bool connected(0);

		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (*p == part) continue; // skip self
			if (p->z1() != part.z2()) continue; // *p not on top of part
			if (!part.intersects_xy(*p)) continue; // no XY overlap
			cube_t shared(part);
			shared.intersect_with_cube(*p); // dz() == 0
			cube_t pref_shared(shared);
			// Note: parts are sorted top to bottom, so any part above <part> should be before it in parts - but we don't want to rely on that here;
			// however, this does mean that the part above this one has already been processed
			float stairs_width(1.2*doorway_width); // relatively small
			float stairs_pad(doorway_width), len_with_pad(stairs_len + 2.0*stairs_pad); // pad both ends of stairs to make sure player has space to enter/exit
			if (max(shared.dx(), shared.dy()) < 1.0*len_with_pad || min(shared.dx(), shared.dy()) < 1.2*stairs_width) continue; // too small to add stairs between these parts

			if (!pri_hall.is_all_zeros() && part.contains_cube(pri_hall) && pri_hall.intersects_xy(shared)) { // have a primary hallway in this part
				pref_shared.intersect_with_cube(pri_hall);
				if (max(pref_shared.dx(), pref_shared.dy()) < 1.2*len_with_pad || min(pref_shared.dx(), pref_shared.dy()) < 1.5*stairs_width) {pref_shared = shared;} // too small
			}
			// place stairs in shared area if there's space and no walls are in the way for either the room above or below
			cube_t cand;
			cand.z1() = part.z2() - window_vspacing + fc_thick; // top of top floor for this part
			cand.z2() = part.z2() + fc_thick; // top of bottom floor of upper part *p

			// is it better to extend the existing stairs in *p, or the stairs we're creating here (stairs_cut) if they line up?
			// iterations: 0-19: place in pri hallway, 20-39: place anywhere, 40-159: shrink size, 150-179: compact stairs, 180-199: allow cut walls
			for (unsigned n = 0; n < 200; ++n) { // make 200 tries to add stairs
				cube_t place_region((n < 20) ? pref_shared : shared); // use preferred shared area from primary hallway for first 20 iterations

				if (n >= 40 && n < 160 && (n%10) == 0) { // decrease stairs size slightly every 10 iterations, 12 times
					stairs_width -= 0.025*doorway_width; // 1.2*DW => 0.9*DW
					stairs_pad   -= 0.030*doorway_width; // 1.0*WD => 0.64*DW
					len_with_pad -= 0.230*doorway_width*(is_basement ? 1.2 : 0.0); // 6.0*DW => 3.24*DW / 2.689*DW ; basement can have steeper stairs
					max_eq(stairs_width, get_min_front_clearance()); // ensure the player can fit
					max_eq(stairs_pad,   get_min_front_clearance()); // ensure the player can fit
				}
				bool dim(0), too_small(0);
				if (min(place_region.dx(), place_region.dy()) < 1.5*len_with_pad) {dim = (place_region.dx() < place_region.dy());} // use larger dim
				else {dim = rgen.rand_bool();}
				bool const stairs_dir(rgen.rand_bool());

				for (unsigned d = 0; d < 2; ++d) {
					float const stairs_sz((bool(d) == dim) ? len_with_pad : stairs_width);
					float const v1(place_region.d[d][0]), v2(place_region.d[d][1] - stairs_sz);
					if (v2 <= v1) {too_small = 1; break;}
					cand.d[d][0] = rgen.rand_uniform(v1, v2); // LLC
					cand.d[d][1] = cand.d[d][0] + stairs_sz; // URC
				}
				if (too_small) continue;
				cube_t cand_test[2] = {cand, cand}; // {lower, upper} parts, starts on lower floor
				cand_test[0].z1() += 0.1*window_vspacing; cand_test[0].z2() -= 0.1*window_vspacing; // shrink to lower part
				cand_test[1].z1() += 1.1*window_vspacing; cand_test[1].z2() += 0.9*window_vspacing; // move to upper part
				cand_test[ stairs_dir].d[dim][0] += stairs_pad; // subtract off padding on one side
				cand_test[!stairs_dir].d[dim][1] -= stairs_pad; // subtract off padding on one side
				// clipped walls don't look right in some cases and may block hallways and rooms, use as a last resort; disable for houses since basement is optional anyway
				bool const allow_clip_walls(n > 180 && !is_house);
				bool bad_place(0), wall_clipped(0);

				for (unsigned d = 0; d < 2; ++d) {
					if (has_bcube_int(cand_test[d], interior->exclusion)) {bad_place = 1; break;} // bad placement
					if (!is_valid_stairs_elevator_placement(cand_test[d], stairs_pad, !allow_clip_walls)) {bad_place = 1; break;} // bad placement
				}
				if (bad_place) continue;

				if (allow_clip_walls) { // clip out walls around stairs
					cand_test[0].z1() = cand.z1(); cand_test[0].z2() = part.z2(); // lower
					cand_test[1].z1() = part.z2(); cand_test[1].z2() = part.z2() + window_vspacing - fc_thick; // upper (bot of ceiling for first floor on upper part)

					for (unsigned e = 0; e < 2; ++e) {
						for (unsigned d = 0; d < 2; ++d) {wall_clipped |= subtract_cube_from_cubes(cand_test[e], interior->walls[d], nullptr, 1);} // clip_in_z=1
					}
				}
				assert(cand.is_strictly_normalized());
				cand.expand_in_dim(dim, -stairs_pad); // subtract off padding
				if (!cand.is_strictly_normalized()) continue; // not enough space, likely because the player radius/front clearance is too large
				// add walls around stairs if room walls were clipped or this is the basement; otherwise, make stairs straight with railings;
				// basement stairs only have walls on the bottom floor, so we set is_at_top=0
				stairs_shape const sshape((is_basement || wall_clipped) ? (stairs_shape)SHAPE_WALLED : (stairs_shape)SHAPE_STRAIGHT);
				landing_t landing(cand, 0, 0, dim, stairs_dir, !wall_clipped, sshape, 0, !is_basement, 1); // roof_access=0, is_at_top=!is_basement, stacked_conn=1
				landing.z1() = part.z2() - fc_thick; // only include the ceiling of this part and the floor of *p
				cube_t stairwell(cand);
				stairwell.z2() = part.z2() + window_vspacing - fc_thick; // bottom of ceiling of upper part; must cover z-range of upper floor for AIs and room object collisions
				interior->landings.push_back(landing);
				interior->stairwells.emplace_back(stairwell, 1, dim, stairs_dir, sshape, 0, 1); // roof_access=0, stack_conn=1

				if (is_basement) { // add a basement door at the bottom of the stairs
					float const pos_shift((stairs_dir ? 1.0 : -1.0)*0.8*wall_thickness);
					door_t door(cand, dim, !stairs_dir, 0, 1); // open=0, on_stairs=1
					door.z2() -= fc_thick; // bottom of basement ceiling, not the floor above
					door.d[dim][stairs_dir] = door.d[dim][!stairs_dir] + pos_shift;
					if (!stairs_dir) {door.translate_dim(dim, -pos_shift);} // why the asymmetry?
					door.translate_dim( dim, -0.2*pos_shift); // shift so that the door doesn't intersect the railing, covers the stairs overhang, and the top edge can't be seen
					door.expand_in_dim(!dim, -0.15*cand.get_sz_dim(dim)/NUM_STAIRS_PER_FLOOR); // shrink by stairs wall half width
					assert(door.is_strictly_normalized());
					interior->stairwells.back().stairs_door_ix = (int16_t)interior->doors.size(); // record door index gating stairs for AI navigation
					add_interior_door(door);
				}
				// attempt to cut holes in ceiling of this part and floor of above part
				cube_t cut_cube(cand);
				cut_cube.z1() += fc_thick; // shrink to avoid clipping floors exactly at the base of the stairs
				subtract_cube_from_floor_ceil(cut_cube, interior->floors);
				subtract_cube_from_floor_ceil(cut_cube, interior->ceilings);

				for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
					if (r->intersects(stairwell) && r->contains_cube_xy(stairwell)) {r->has_stairs = 1;} // Note: may be approximate
				}
				connected = 1;
				break; // success
			} // for n
			if (connected && is_basement) break; // only need to connect one part for the basement
		} // for p
		if (!connected && is_basement) {interior->is_unconnected = 1;} // failed to connect basement with stairs - flag as unconnected
	}

	// now attempt to extend elevators into floors above/below
	vect_cube_t holes;

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find parts above/below this elevator
			cube_t needed_space(*e);
			needed_space.expand_by_xy(fc_thick);
			if (!p->contains_cube_xy(needed_space)) continue; // elevator doesn't extend inside this cube
			bool const is_above(p->z1() == e->z2()), is_below(p->z2() == e->z1());
			if (!is_above && !is_below) continue;
			cube_t extension(*e);
			extension.z1() = p->z1(); extension.z2() = p->z2();
			cube_t cand_test(extension);
			cand_test.expand_in_dim(2, -fc_thick); // shrink slightly in Z so that we don't intersect the original elevator *e
			cand_test.d[e->dim][e->dir] += doorway_width*(e->dir ? 1.0 : -1.0); // add extra space in front of the elevator
			if (!p->contains_cube_xy(cand_test)) continue; // not enough space at elevator entrance
			bool const allow_clip_walls = 1; // optional
			if (!is_valid_stairs_elevator_placement(cand_test, doorway_width, !allow_clip_walls)) continue; // bad placement
			if (has_bcube_int(cand_test, interior->exclusion)) continue; // bad placement

			if (allow_clip_walls) { // clip out walls around extended elevator
				for (unsigned d = 0; d < 2; ++d) {
					holes.clear();
					subtract_cube_from_cubes(cand_test, interior->walls[d], (((bool)d == e->dim) ? &holes : nullptr));

					// see if we need to add any walls to close gaps created between the front of the elevator and a removed wall section
					for (auto h = holes.begin(); h != holes.end(); ++h) {
						float const elevator_front(e->d[e->dim][e->dir]), inner_face(h->d[e->dim][!e->dir]), outer_face(h->d[e->dim][e->dir]);
						if ((inner_face > elevator_front) ^ e->dir) continue; // wall not in front of the elevator
						
						for (unsigned g = 0; g < 2; ++g) { // left and right of the elevator
							float const elevator_edge(e->d[!e->dim][g]);
							if (h->d[!e->dim][g] != elevator_edge) continue; // not clipped to this edge
							cube_t wall(*h); // copy zvals from hole
							wall.d[ e->dim][ e->dir] = outer_face;
							wall.d[ e->dim][!e->dir] = elevator_front;
							wall.d[!e->dim][ g] = elevator_edge;
							wall.d[!e->dim][!g] = elevator_edge + (g ? -1.0 : 1.0)*wall_thickness;
							interior->walls[!d].push_back(wall);
						} // for g
					} // for h
				} // for d
			}
			min_eq(e->z1(), extension.z1()); max_eq(e->z2(), extension.z2()); // perform extension in Z
			float const shift((is_above ? -1.1 : 1.1)*fc_thick);
			extension.z1() += shift; extension.z2() += shift; // also cut a hole in the lower ceiling/upper floor
			holes.clear();
			subtract_cube_from_cubes(extension, interior->ceilings);
			subtract_cube_from_cubes(extension, interior->floors, &holes); // capture holes from floors
			if (is_above) {add_or_extend_elevator(*e, 0);} // extend only
			
			for (auto h = holes.begin(); h != holes.end(); ++h) {
				landing_t landing(*h, 1, 0, e->dim, e->dir);
				landing.z1() -= fc_thick; // since we only captured floor cutouts, extend them downward to include the ceiling below
				interior->landings.push_back(landing);
			}
			for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
				if (r->intersects(cand_test)) {r->has_elevator = 1;} // find all intersecting rooms and set has_elevator flag
			}
		} // for p
	} // for e
}

bool building_t::are_parts_stacked(cube_t const &p1, cube_t const &p2) const {
	if (is_house) return 0; // houses are never stacked
	if (p1.z2() == p2.z1() && p1.contains_cube_xy(p2)) return 1; // p2 stacked on p1
	if (p2.z2() == p1.z1() && p2.contains_cube_xy(p1)) return 1; // p1 stacked on p2
	return 0;
}

bool building_t::clip_part_ceiling_for_stairs(cube_t const &c, vect_cube_t &out, vect_cube_t &temp) const { // and elevators
	if (!interior || interior->stairwells.empty()) return 0;
	subtract_cubes_from_cube(c, interior->stairwells, out, temp);

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) { // handle elevators that span multiple parts
		if (e->z2() <= c.z1()) continue; // elevator shaft doesn't extend this high
		if (e->z2() == c.z2()) continue; // don't cut a hole in the roof of the building where the top of the elevator shaft ends
		subtract_cube_from_cubes(*e, out);
	}
	return 1;
}

void building_t::add_room(cube_t const &room, unsigned part_id, unsigned num_lights, bool is_hallway, bool is_office, bool is_sec_bldg) {
	assert(interior);
	assert(room.is_strictly_normalized());
	room_t r(room, part_id, num_lights, is_hallway, is_office, is_sec_bldg);
	cube_t const &part(parts[part_id]);

	for (unsigned d = 0; d < 4; ++d) { // find exterior sides
		bool const dim(d>>1), dir(d&1);
		float const part_edge(part.d[dim][dir]);
		if (room.d[dim][dir] != part_edge) continue; // interior to this part
		bool is_exterior(1);

		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // simplified version of classify_room_wall()
			if (p->d[dim][!dir] != part_edge) continue; // not opposite wall
			if (p->z1() >= room.z2() || p->z2() <= room.z1()) continue; // no z overlap (wrong stack)
			if (p->d[!dim][0] > room.d[!dim][0] || p->d[!dim][1] < room.d[!dim][1]) continue; // wall not contained
			if (room.z2() <= p->z2()) {is_exterior = 0; break;} // this part covers the wall in z (assuming no overhangs), so wall is interior split between parts (not exterior)
		} // for p
		if (is_exterior) {r.ext_sides |= (1 << d);}
	} // for d
	interior->rooms.push_back(r);
}

void building_t::add_or_extend_elevator(elevator_t const &elevator, bool add) {
	if (add) {interior->elevators.push_back(elevator);}
	if (is_house || roof_type != ROOF_TYPE_FLAT || has_helipad) return; // sloped roof, not flat, can't add elevator cap
	float const window_vspacing(get_material().get_floor_spacing());
	cube_t ecap(elevator);
	ecap.z1()  = elevator.z2();
	ecap.z2() += 0.25*window_vspacing; // set height
	if (!elevator.at_edge) {ecap.expand_by_xy(0.025*window_vspacing);}
	
	// check to see if the elevator is at the top of the building
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		if (p->z1() != elevator.z2()) continue; // not on top of the elevator
		if (p->intersects(ecap)) return; // part over elevator - should we add some sort of cap in this case, or block off the first floor, or add something to the interior?
	}
	remove_intersecting_roof_cubes(ecap);
	details.emplace_back(ecap, ROOF_OBJ_ECAP);
	max_eq(bcube.z2(), ecap.z2()); // extend bcube z2 to contain ecap
}

void building_t::remove_intersecting_roof_cubes(cube_t const &c) {
	for (unsigned i = 0; i < details.size(); ++i) { // remove any existing objects that overlap ecap
		auto &obj(details[i]);
		if (obj.type != ROOF_OBJ_BLOCK && obj.type != ROOF_OBJ_AC && obj.type != ROOF_OBJ_ANT) continue; // only remove blocks, AC units, and antennas
		if (!obj.intersects(c)) continue;
		swap(obj, details.back());
		details.pop_back();
		--i; // wraparound okay
	}
}
