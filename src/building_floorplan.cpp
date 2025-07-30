// 3D World - Building Interior Floorplan Generation (walls, ceilings, floors, rooms, and doors)
// by Frank Gennari 2/20/20

#include "function_registry.h"
#include "buildings.h"
#include "profiler.h"

unsigned const MAX_OFFICE_UTILITY_ROOMS = 1;

extern building_params_t global_building_params;


void building_t::add_interior_door(door_t &door, bool is_bathroom, bool make_unlocked, bool make_closed) { // add door stack + doors
	assert(interior);
	interior->door_stacks.emplace_back(door, interior->doors.size());
	if (door.on_stairs) {add_interior_door_for_floor(door, is_bathroom, make_unlocked); return;} // add a single door across all floors
	float const floor_spacing(get_window_vspace()), door_height(get_floor_ceil_gap());

	// Note: door.dz() should be an exact multiple of floor_spacing except for an extra floor thickness at the bottom
	for (float zval = door.z1(); zval + 0.5f*floor_spacing < door.z2(); zval += floor_spacing) { // continue until we don't have enough space left to add a door
		door_t door_seg(door);
		set_cube_zvals(door_seg, zval, zval+door_height); // clip to ceiling
		add_interior_door_for_floor(door_seg, is_bathroom, make_unlocked, make_closed);
	}
	interior->door_stacks.back().num_doors = (interior->doors.size() - interior->door_stacks.back().first_door_ix);
}
void building_t::add_interior_door_for_floor(door_t &door, bool is_bathroom, bool make_unlocked, bool make_closed) {
	if (is_bathroom) {door.open = 0; door.locked = 0;} // bathroom doors are always closed but unlocked
	else if (!door.on_stairs) { // don't set open/locked state for stairs doors
		door.open   = (!make_closed   &&               fract(interior->doors.size()*1.61803) < global_building_params.open_door_prob  ); // use the golden ratio
		door.locked = (!make_unlocked && !door.open && fract(interior->doors.size()*3.14159) < global_building_params.locked_door_prob); // use pi
	}
	door.make_fully_open_or_closed();
	interior->doors.push_back(door);
}

void building_t::remove_section_from_cube_and_add_door(cube_t &c, cube_t &c2, float v1, float v2, bool xy,
	bool open_dir, bool is_bathroom, bool make_unlocked, bool make_closed, bool jail_door)
{
	// remove a section from this cube; c is input+output cube, c2 is other output cube
	assert(v1 < v2);
	assert(v1 > c.d[xy][0] && v2 < c.d[xy][1]); // v1/v2 must be interior values for cube
	bool const hinge_side(xy ^ open_dir ^ (v1-c.d[xy][0] < c.d[xy][1]-v1) ^ 1); // put the hinge on the side closer to the end of the wall
	c2 = c; // clone first cube
	c.d[xy][1] = v1; c2.d[xy][0] = v2; // c=low side, c2=high side
	// add a door stack and doors
	door_t door(c, !xy, open_dir, 1, 0, hinge_side); // open=1, on_stairs=0
	door.d[!xy][0] = door.d[!xy][1] = c.get_center_dim(!xy); // zero area at wall centerline
	door.d[ xy][0] = v1;
	door.d[ xy][1] = v2;
	if (jail_door) {door.for_jail = 2;} // prison bar doors; opaque with a barred window and a frame
	add_interior_door(door, is_bathroom, make_unlocked, make_closed);
}

void building_t::insert_door_in_wall_and_add_seg(cube_t &wall, float v1, float v2, bool dim, bool open_dir,
	bool keep_high_side, bool is_bathroom, bool make_unlocked, bool make_closed, bool jail_door)
{
	cube_t wall2;
	remove_section_from_cube_and_add_door(wall, wall2, v1, v2, dim, open_dir, is_bathroom, make_unlocked, make_closed, jail_door);
	if (keep_high_side) {swap(wall, wall2);} // swap left and right
	interior->walls[!dim].push_back(wall2);
}

cube_t door_base_t::get_open_door_path_bcube() const { // independent of room
	cube_t bcube(get_true_bcube());
	float const width(get_width());
	bool const dir(get_check_dirs());
	bcube.d[!dim][dir     ] += (dir      ? 1.0 : -1.0)*width; // include door fully open position
	bcube.d[ dim][open_dir] += (open_dir ? 1.0 : -1.0)*width;
	return bcube;
}
void building_t::reverse_door_hinges_if_needed() { // quadratic in the number of door stacks, which should be relatively small
	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
		if (i->on_stairs) continue;
		cube_t const bc1(i->get_open_door_path_bcube());
		bool reverse_hinges(0);

		for (auto j = interior->door_stacks.begin(); j != interior->door_stacks.end(); ++j) {
			if (i == j || j->on_stairs) continue; // skip self and stairs doors
			if (j->z1() >= i->z2() || j->z2() <= i->z1()) continue; // no Z overlap
			cube_t const bc2(j->get_open_door_path_bcube());
			if (bc1.intersects_no_adj(bc2)) {reverse_hinges = 1; break;}
		}
		if (reverse_hinges) {
			i->hinge_side ^= 1; // reverse the door stack's hinges
			assert(i->first_door_ix < interior->doors.size());
			
			for (unsigned dix = i->first_door_ix; dix < interior->doors.size(); ++dix) { // reverse the doors themselves
				door_t &door(interior->doors[dix]);
				if (!i->is_same_stack(door)) break; // moved to a different stack, done
				door.hinge_side ^= 1;
			}
		}
	} // for i
}

void building_t::ensure_doors_to_room_are_closed(room_t const &room, unsigned doors_start, bool ensure_locked) {
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), wall_thick(get_wall_thickness());
	cube_t test_cube(room);
	test_cube.expand_by_xy(wall_thick);
	set_cube_zvals(test_cube, room.z1()+floor_thickness, room.z1()+window_vspacing-floor_thickness); // shrink to first floor

	for (auto i = interior->doors.begin()+doors_start; i != interior->doors.end(); ++i) {
		if (i->open && i->get_true_bcube().intersects(test_cube)) {i->open = 0;} // make sure door is closed
		i->make_fully_open_or_closed();
		if (ensure_locked) {i->locked = 1;}
	}
}

float cube_rand_side_pos(cube_t const &c, bool dim, float min_dist_param, float min_dist_abs, rand_gen_t &rgen, bool for_door) {
	assert(min_dist_param < 0.5f); // aplies to both ends
	float const lo(c.d[dim][0]), hi(c.d[dim][1]);
	
	if (for_door && global_building_params.put_doors_in_corners) { // place near the wall to keep doors to the edges of rooms; not using min_dist_param
		float const gap(0.6*min_dist_abs), v1(lo + gap), v2(hi - gap), shift(0.25*min_dist_param*rgen.rand_float());
		bool const side(rgen.rand_bool());
		return max(v1, min(v2, (side ? (v2 - shift) : (v1 + shift))));
	}
	float const delta(hi - lo), gap(max(min_dist_abs, min_dist_param*delta)), v1(lo + gap), v2(hi - gap);
	//if (v2 <= v1) {cout << TXT(dim) << TXT(lo) << TXT(hi) << TXT(min_dist_abs) << TXT(delta) << TXT(gap) << endl;}
	if (v1 >= v2) {return 0.5f*(v1 + v2);} // if range is denormalized, use the center
	return rgen.rand_uniform(v1, v2);
}

// see global_building_params.window_xspace/window_width
int building_t::get_num_windows_on_side(float xy1, float xy2) const {
	assert(xy1 < xy2);
	float const tscale(get_material().get_floorplan_window_xscale());
	float t0(tscale*xy1), t1(tscale*xy2);
	clip_low_high_tc(t0, t1);
	return round_fp(t1 - t0);
}
float building_t::get_window_h_border() const {return 0.5*(1.0 - global_building_params.get_window_width_fract ());} // (0.0, 1.0)
float building_t::get_window_v_border() const {return 0.5*(1.0 - global_building_params.get_window_height_fract());} // (0.0, 1.0)
float building_t::get_hspacing_for_part(cube_t const &part, bool dim) const {return part.get_sz_dim(dim)/get_num_windows_on_side(part, dim);}

void set_wall_width(cube_t &wall, float pos, float half_thick, unsigned dim) {
	wall.d[dim][0] = pos - half_thick;
	wall.d[dim][1] = pos + half_thick;
}
void resize_around_center_xy(cube_t &c, float radius) {
	for (unsigned d = 0; d < 2; ++d) {set_wall_width(c, c.get_center_dim(d), radius, d);}
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

// Note: assumes edge is clipped to a whole window
bool is_val_inside_window(cube_t const &c, bool dim, float val, float window_spacing, float window_border) {
	window_border *= 0.9; // adjust based on window frame so that wall doesn't end right at the edge
	float const uv(fract((val - c.d[dim][0])/window_spacing));
	return (uv > window_border && uv < 1.0f-window_border);
}
// shift_edges_mode: 0=normal, 1=opposite dir, 2=none
float shift_val_to_not_intersect_window(cube_t const &c, float val, float hspace, float window_border, bool dim, int shift_edges_mode) {
	if (shift_edges_mode == 2) return val; // no shift
	window_border *= 0.9; // adjust based on window frame so that wall doesn't end right at the edge
	float const uv(fract((val - c.d[dim][0])/hspace));
	if (!(uv > window_border && uv < 1.0f-window_border)) return val; // okay as is
	float const uv_target(((uv < 0.5) ^ bool(shift_edges_mode)) ? window_border : 1.0f-window_border); // move to the closer edge, unless shift_other_dir=1
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

void subtract_cube_xy(cube_t const &c, cube_t const &ri, cube_t *out) { // subtract r from c; ignores zvals
	//assert(c.contains_cube_xy(ri));
	cube_t r(ri);
	r.intersect_with_cube_xy(c);
	for (unsigned i = 0; i < 4; ++i) {out[i] = c;}
	out[0].y2() = r.y1(); // bottom -y
	out[1].y1() = r.y2(); // top    +y
	out[2].y1() = r.y1(); out[2].y2() = r.y2(); out[2].x2() = r.x1(); // left  center -x
	out[3].y1() = r.y1(); out[3].y2() = r.y2(); out[3].x1() = r.x2(); // right center +x
}

bool building_t::interior_enabled() const {
	if (world_mode != WMODE_INF_TERRAIN)                return 0; // tiled terrain mode only
	if (!global_building_params.gen_building_interiors) return 0; // disabled
	if (!global_building_params.windows_enabled())      return 0; // no windows, can't assign floors and generate interior
	if (!is_cube() && has_complex_floorplan)            return 0; // not handling non-cube buildings with complex floorplans here
	// skip buildings with windows baked into textures
	if (!global_building_params.add_city_interiors && !has_windows()) return 0; // not a building type that has generated windows
	return 1;
}

int building_t::classify_room_wall(room_t const &room, float zval, bool dim, bool dir, bool ret_sep_if_part_int_part_ext) const { // Note: zval is for the floor
	if (room.is_ext_basement()) return ROOM_WALL_INT; // treated as interior wall (optimization)
	//if (has_basement() && zval < ground_floor_z1) return ROOM_WALL_BASEMENT; // callers don't expect this, but we may want to add/handle this case in the future
	if (room.d[dim][dir] == bcube.d[dim][dir]) return ROOM_WALL_EXT; // at bcube border
	if (!ret_sep_if_part_int_part_ext && (room.ext_sides & (1 << (2*dim + dir)))) return ROOM_WALL_EXT; // use ext_sides (optimization, but may conservatively return ext)
	cube_t const &part(get_part_for_room(room));
	float const wall_thickness(get_wall_thickness()), max_gap(1.5*wall_thickness), part_edge(part.d[dim][dir]);
	if (dir ? (room.d[dim][1] + max_gap < part_edge) : (room.d[dim][0] - max_gap > part_edge)) return ROOM_WALL_INT; // interior to part, allowing for wall_thickness of gap
	if (real_num_parts == 1) return ROOM_WALL_EXT; // optimization

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
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
unsigned building_t::count_ext_walls_for_room(room_t const &room, float zval) const {
	unsigned num_ext_walls(0);

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {num_ext_walls += (classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT);}
	}
	return num_ext_walls;
}

void building_t::gen_interior(rand_gen_t &rgen, bool has_overlapping_cubes) { // Note: contained in building bcube, so no bcube update is needed
	if (!interior_enabled()) return;
	unsigned const details_size(details.size()), roof_tquads_size(roof_tquads.size()), doors_size(doors.size());

	// make up to 16 attempts to generate a connected interior; the first attempt almost always succeeds; currently it only fails if a basement is unconnected
	// 64 house basements are unconnected, this loop fixes all but 6 of them after 3 iterations; 12 iterations is required for 100% success
	for (unsigned n = 0; n < 16; ++n) {
		// remove any objects added in the previous (failed) iteration
		details.resize(details_size);
		roof_tquads.resize(roof_tquads_size);
		doors.resize(doors_size);
		ext_lights.clear(); // generated as part of the interior
		gen_interior_int(rgen, has_overlapping_cubes);
		if (!interior->is_unconnected) break; // done
	} // for n
	for (unsigned d = 0; d < 2; ++d) {interior->extb_walls_start[d] = interior->walls[d].size();}
	if (is_parking()) {add_parking_roof_lights();}
	// calculate and cache interior_z2
	interior_z2 = ground_floor_z1;
	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {max_eq(interior_z2, i->z2());}
}

// Note: these are used in gen_interior_int() and maybe_add_skylight()
float building_t::get_min_hallway_width() const {
	// need wider hallway for U-shaped stairs, but not quite the 6.0 factor used in add_ceilings_floors_stairs()
	return ((has_tall_retail() ? 5.4 : 3.6)*get_nominal_doorway_width());
}
bool building_t::can_use_hallway_for_part(unsigned part_id) const {
	if (is_house || has_complex_floorplan || !is_cube() || is_parking() || (int)part_id == basement_part_ix) return 0;
	assert(part_id < parts.size());
	cube_t const &p(parts[part_id]);
	bool const first_part_this_stack(part_id == 0 || parts[part_id-1].z1() < p.z1());
	if (!first_part_this_stack) return 0;
	bool const next_diff_stack(part_id+1 == parts.size() || parts[part_id+1].z1() != p.z1());
	if (!next_diff_stack)       return 0;
	float const min_wall_len(get_min_wall_len());
	return (min(p.dx(), p.dy()) > max(4.0f*min_wall_len, (get_min_hallway_width() + 2.0f*min_wall_len)));
}
cube_t building_t::get_hallway_for_part(cube_t const &part, float &num_hall_windows, float &hall_width, float &room_width) {
	bool const min_dim(part.dy() < part.dx());
	int const num_windows_od(get_num_windows_on_side(part, min_dim)); // in short dim
	float const part_width(part.get_sz_dim(min_dim)), min_hall_width(get_min_hallway_width()); // need wider hallway for U-shaped stairs
	bool const is_odd(num_windows_od & 1);
	num_hall_windows = (is_odd ? 1.4 : 1.8); // hall either contains 1 (odd) or 2 (even) windows, wider for single window case to make room for stairs
	max_eq(num_hall_windows, min_hall_width*num_windows_od/part_width); // enforce min_hall_width (may split a window, but this limit is only hit for non-window city office buildings)
	//if (is_odd) {min_eq(num_hall_windows, 1.5f);} // hard limit for single window case to avoid hall walls clipping through windows
	hall_width = num_hall_windows*part_width/num_windows_od;
	room_width = 0.5f*(part_width - hall_width); // rooms are the same size on each side of the hallway
	cube_t hall(part);
	hall.expand_in_dim(min_dim, -room_width); // shink rooms off of each end
	return hall;
}

void building_t::gen_interior_int(rand_gen_t &rgen, bool has_overlapping_cubes) { // Note: contained in building bcube, so no bcube update is needed

	// defer this until the building is close to the player?
	interior.reset(new building_interior_t);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const doorway_width(get_nominal_doorway_width()), doorway_hwidth(0.5*doorway_width);
	float const wall_thick(get_wall_thickness()), wall_half_thick(0.5*wall_thick), wall_edge_spacing(0.05*wall_thick), min_wall_len(get_min_wall_len());
	float const window_border(get_window_h_border()), door_to_wall_dist(1.5*doorway_hwidth + wall_thick);
	bool const put_doors_in_corners(global_building_params.put_doors_in_corners);
	vector3d const car_sz(get_nom_car_size());
	point bldg_door_open_dir_tp(bcube.get_cube_center()); // used to determine in which direction doors open; updated base on central hallway
	// houses have at most two parts; exclude garage, shed, porch, porch support, etc.
	auto parts_end(get_real_parts_end());
	vector<split_cube_t> to_split;
	uint64_t must_split[2] = {0,0};
	unsigned first_wall_to_split[2] = {0,0};
	cube_t pref_conn_to; // house hallway, etc.
	// allocate space for all floors; this is now likely a large undercount of the actual number of objects needed
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
	interior->int_door_width = get_doorway_width(); // set the default value until doors are added
	interior->ceilings  .reserve(tot_num_floors);
	interior->floors    .reserve(tot_num_floors);
	interior->landings  .reserve(tot_num_landings);
	interior->stairwells.reserve(tot_num_stairwells);
	vector<room_t> &rooms(interior->rooms);
	
	// generate walls and floors for each part;
	// this will need to be modified to handle buildings that have overlapping parts, or skip those building types completely
	for (auto p = parts.begin(); p != parts_end; ++p) {
		unsigned num_floors(calc_num_floors(*p, window_vspacing, floor_thickness)), part_id(p - parts.begin());
		if (num_floors == 0) continue; // not enough space to add a floor (can this happen?)
		if (is_retail_part(*p)) {num_floors = 1;} // retail area is always one floor
		// for now, assume each part has the same XY bounds and can use the same floorplan; this means walls can span all floors and don't need to be duplicated for each floor
		vector3d const psz(p->get_size());
		bool const is_basement_part(is_basement(p)), first_part(part_id == 0), first_part_this_stack(first_part || is_basement_part || (p-1)->z1() < p->z1());
		// office building hallways only; house hallways are added later
		bool const is_industrial_part(is_industrial() && first_part);
		if (is_industrial_part) {num_floors = 1;} // industrial buildings are a single floor
		bool const use_hallway(!is_industrial_part && can_use_hallway_for_part(part_id)), min_dim(psz.y < psz.x);
		unsigned const rooms_start(rooms.size()), doors_start(interior->doors.size()), num_doors_per_stack(num_floors);
		cube_t hall, place_area(*p);
		place_area.expand_by_xy(-wall_edge_spacing); // shrink slightly to avoid z-fighting with walls
		bool no_split_walls_this_part(0);
		float window_hspacing   [2] = {0.0};
		int num_windows_per_side[2] = {0};
		vector<unsigned> utility_room_cands, special_room_cands; // used for parts with a hallway + assign_special_room_types()

		for (unsigned d = 0; d < 2; ++d) {
			num_windows_per_side[d] = get_num_windows_on_side(*p, d);
			window_hspacing     [d] = psz[d]/num_windows_per_side[d];
		}
		if (!is_cube()) { // cylinder, etc.
			float const min_dim_sz(min(p->dx(), p->dy()));
			bool const is_office(!is_house);

			if (min_dim_sz > 2.0*min_wall_len) { // large cylinder or N-gon
				// create a pie slice split for cylindrical parts; since we can only add X or Y walls, place one of each that crosses the entire part
				point const center(p->xc(), p->yc(), p->z1());
				bool const clip_to_ext_walls(!use_cylinder_coll()); // if this building isn't a full cylinder, we may need to clip the walls shorter
				vector<point> const &points(get_part_ext_verts(part_id));
				// if the part is on the small side, only add a single wall; choose randomly; maybe should add the shorter wall if asymmetric?
				bool const only_one_wall(min_dim_sz < 3.0*min_wall_len);
				unsigned const N((flat_side_amt == 0.0) ? num_sides : 4); // assume worst case of flat side making a right angle
				float const corner_angle((N-2)*PI/N);
				float const end_pullback(1.5*(wall_half_thick + get_trim_thickness())/tan(0.5*corner_angle));
				unsigned skip_dim(only_one_wall ? rgen.rand_bool() : 2);

				for (unsigned d = 0; d < 2; ++d) {
					if (d == skip_dim) continue;
					cube_t wall(*p);
					clip_wall_to_ceil_floor(wall, fc_thick);
					set_wall_width(wall, center[d], wall_half_thick, d);
					
					if (clip_to_ext_walls) {
						for (unsigned e = 0; e < 2; ++e) { // for each end
							for (auto p = points.begin(); p != points.end(); ++p) { // check each edge - only one can intersect
								point const &p1(*p), &p2((p+1 == points.end()) ? points.front() : *(p+1));
								if (max(p1[d], p2[d]) <= center[d] || min(p1[d], p2[d]) >= center[d]) continue; // doesn't cross the wall in dim d
								float const t((p1[d] - center[d])/(p1[d] - p2[d])), int_pos(p1[!d] + t*(p2[!d] - p1[!d]));
								if      ( e && int_pos > center[!d]) {min_eq(wall.d[!d][1], int_pos); break;}
								else if (!e && int_pos < center[!d]) {max_eq(wall.d[!d][0], int_pos); break;}
							} // for p
						} // for e
					}
					wall.expand_in_dim(!d, -end_pullback); // shrink slightly to avoid clipping through the exterior wall
					assert(wall.is_strictly_normalized());

					// cut two doorways in each wall if there's space
					for (unsigned e = 0; e < 2; ++e) {
						if (fabs(wall.d[!d][e] - center[!d]) < 2.0*doorway_width) continue; // wall too short to add a door
						float const door_pos(0.5*(wall.d[!d][e] + center[!d])); // midpoint of the half-wall
						insert_door_in_wall_and_add_seg(wall, (door_pos - doorway_hwidth), (door_pos + doorway_hwidth), !d, 0, 1); // keep_high_side=1
					}
					interior->walls[d].push_back(wall); // add remainder
				} // for d
				// now add 2-4 rooms
				unsigned br_floors_used(0); // bit mask for bathrooms

				for (unsigned r = 0; r < 4; ++r) {
					bool const xside(r & 1), yside(r >> 1);
					if ((skip_dim == 0 && xside) || (skip_dim == 1 && yside)) continue; // part not split in this dim
					cube_t room(*p);
					float const shrink_val(is_office ? 0.0 : wall_half_thick); // wall not included in room bounds unless this is an office
					if (skip_dim != 0) {room.d[0][xside] = center.x - (xside ? 1.0 : -1.0)*shrink_val;}
					if (skip_dim != 1) {room.d[1][yside] = center.y - (yside ? 1.0 : -1.0)*shrink_val;}
					add_room(room, part_id, 1, 0, is_office);
					
					// assign one floor of this room as a bathroom unless there are only a few floors; should guarantee each building has at least one bathroom
					if (r <= num_floors) {
						unsigned const floors_end(min(num_floors, NUM_RTYPE_SLOTS-1)); // not too high a floor so that slot clamp occurs
						unsigned const floors_start((floors_end < 3) ? 0 : 1); // skip first floor to avoid ext doors (which haven't been placed yet), unless 1-2 floors
						unsigned floor_ix(floors_start + (rgen.rand() % (floors_end - floors_start)));
						if (!(br_floors_used & (1<<floor_ix))) {rooms.back().assign_to(RTYPE_BATH, floor_ix);}
						br_floors_used |= (1<<floor_ix); // at most one bathroom per floor
					}
				} // for r
			}
			else {
				add_room(*p, part_id, 1, 0, is_office); // add entire part as a room; num_lights will be calculated later
			}
		} // end non-cube room
		else if (is_industrial_part) { // first part is industrial space: factory, warehouse, power plant, etc.
			create_industrial_floorplan(part_id, window_hspacing, window_border, rgen);
		}
		else if (!is_house && is_basement_part && min(psz.x, psz.y) > 5.0*car_sz.x && max(psz.x, psz.y) > 12.0*car_sz.y) { // make this a parking garage
			add_room(*p, part_id); // add entire part as a room; num_lights will be calculated later
			rooms.back().assign_all_to(RTYPE_PARKING); // make it a parking garage
			has_parking_garage = 1;
		}
		else if (has_retail() && part_id == 0) {
			add_room(*p, part_id); // add entire part as a room; num_lights will be calculated later
			rooms.back().assign_all_to(RTYPE_RETAIL);
			rooms.back().is_single_floor = 1;
		}
		else if (is_parking() && part_id == 0) { // parking structure
			cube_t room(*p); // start with full part
			room.expand_by_xy(-get_park_struct_wall_thick()); // shrink off exterior walls
			add_room(room, part_id); // add entire part as a room; num_lights will be calculated later
			rooms.back().assign_all_to(RTYPE_PARKING);
		}
		else if (use_hallway) {
			// building with rectangular slice (no adjacent exterior walls at this level), generate rows of offices
			// Note: we could probably make these unsigned, but I want to avoid unepected negative numbers in the math;
			bool const apt_or_hotel(is_apt_or_hotel()), is_a_school(is_school()); // use a different floorplan for apartments, hotels, and schools
			int const num_windows   (num_windows_per_side[!min_dim]);
			int const num_windows_od(num_windows_per_side[ min_dim]); // other dim, for use in hallway width calculation
			int windows_per_room((num_windows >= 7 && num_windows_od >= 7) ? 2 : 1); // 1-2 windows per room (only assign 2 windows if we can get into the secondary hallway case below)
			float const hall_len(psz[!min_dim]), wind_hspacing(hall_len/num_windows);
			float room_len(wind_hspacing*windows_per_room);

			while (room_len < 0.9*min_wall_len) { // add more windows to increase room size if too small
				++windows_per_room;
				room_len = wind_hspacing*windows_per_room;
			}
			bool const partial_room((num_windows % windows_per_room) != 0 && !apt_or_hotel && !is_a_school); // an odd number of windows leaves a small room at the end
			bool const is_ground_floor(is_ground_floor_excluding_retail(p->z1()));
			float const hwall_extend(0.5f*(room_len - doorway_width - wall_thick));
			float const wall_pos(p->d[!min_dim][0] + room_len); // pos of first wall separating first from second rooms
			float num_hall_windows, hall_width, room_width;
			hall = get_hallway_for_part(*p, num_hall_windows, hall_width, room_width);
			float const *hall_wall_pos(hall.d[min_dim]);
			// round down for apts/hotels/schools, otherwise round up
			int const num_rooms(((apt_or_hotel || is_a_school) ? num_windows : (num_windows+windows_per_room-1))/windows_per_room);
			int const windows_per_side_od((num_windows_od - round_fp(num_hall_windows))/2); // of hallway
			assert(num_rooms >= 0 && num_rooms < 1000); // sanity check
			auto& room_walls(interior->walls[!min_dim]), & hall_walls(interior->walls[min_dim]); // room_walls: perpendicular to hallway; hall_walls: parallel to hallway
			unsigned const hall_walls_start(hall_walls.size());
			if (hallway_dim == 2) {hallway_dim = !min_dim;} // cache in building for later use, only for first part (ground floor)
			
			// add secondary or ring hallways if there are at least 7 windows (3 on each side of hallway); not for apartments and hotels; schools must have an even number of windows
			bool const has_sec_hallways(!apt_or_hotel && num_windows_od >= 7 && num_rooms >= 4 && (!is_a_school || !(windows_per_side_od & 1)));

			if (has_sec_hallways) {
				float const min_hall_width(1.5f*doorway_width), max_hall_width(2.5f*doorway_width);
				float const sh_width(max(min(0.4f*hall_width, max_hall_width), min_hall_width)), hspace(window_hspacing[!min_dim]);
				float const ring_hall_room_depth(0.5f*(room_width - sh_width)); // for inner and outer rows of rooms
				interior->has_sec_hallways = 1;

				// use a ring hallway, if large enough and not a school (because it creates more small and windowless rooms)
				// Note: the ring_hall_room_depth check can fail for two level retail areas that force wider hallways for U-shaped stairs
				if (ring_hall_room_depth > 2.0f*doorway_width && hall_len > 12.5*window_vspacing && !is_a_school && rgen.rand_bool()) {
					float const hall_offset(room_len); // outer edge of hallway to building exterior on each end
					bool const add_doors_to_main_wall(rgen.rand_bool());
					unsigned const num_cent_rooms(num_rooms - 2); // skip the rooms on each side
					unsigned const num_doors_inner_rooms(add_doors_to_main_wall ? 3U : 2U);
					unsigned const num_offices(4*(num_cent_rooms + 4)); // X/Y mirror symmetry
					unsigned const min_num_doors(num_offices + 2*add_doors_to_main_wall*num_cent_rooms + 4); // at least one per office
					assert(num_cent_rooms > 0);
					hall_walls.reserve(2*(11 + num_doors_inner_rooms*num_cent_rooms)); // long dim (along hall dir)
					room_walls.reserve(2*(10 + 2*num_cent_rooms)); // short dim
					rooms.reserve(num_offices + 7); // num_offices + pri hall + 2 sec hall + 4 conn hall
					interior->door_stacks.reserve(min_num_doors);
					interior->doors.reserve(num_doors_per_stack*min_num_doors);
					interior->exclusion.reserve(6); // 2 sec hallways + 4 conn hallways
					unsigned const bathroom_ix(rgen.rand_bool() ? 0 : num_cent_rooms-1); // place bathrooms on one of the end center rooms

					for (unsigned d = 0; d < 2; ++d) { // for each side of main hallway
						float const dsign(d ? -1.0 : 1.0);
						float const wall_edge(p->d[min_dim][d]), targ_hall_outer(wall_edge + dsign*ring_hall_room_depth);
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
						unsigned const num_lights(min(6U, max(2U, unsigned(0.5*s_hall.get_sz_dim(!min_dim)/s_hall.get_sz_dim(min_dim)))));
						add_room(s_hall, part_id, num_lights, 1, 0); // add sec hallway as room with several lights
						rooms.back().mark_open_wall(min_dim, !d); // adjacent to connector hallways
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
						float const rooms_range[2] = {short_swall.d[!min_dim][0], short_swall.d[!min_dim][1]};
						float const span(rooms_range[1] - rooms_range[0]), cent_room_width(span/num_cent_rooms);
						cube_t room_outer(*p), room_inner(*p);
						room_outer.d[min_dim][!d] = hall_outer; // other dim stays at wall_edge
						room_inner.d[min_dim][ d] = hall_inner;
						room_inner.d[min_dim][!d] = hall_wall_pos[d];
						int shift_edges_mode[2] = {0,0}; // shift: 0=normal, 1=opposite, 2=none

						for (unsigned e = 0; e < 2; ++e) { // for each end of building in long dim
							float const esign(e ? -1.0 : 1.0); // Note: offset_inner == (rooms_range[e] + esign*wall_half_thick)
							float const other_edge(p->d[!min_dim][e]), offset_outer(s_hall.d[!min_dim][e]), offset_inner(offset_outer + esign*sh_width);
							c_hall.d[!min_dim][ e] = offset_outer;
							c_hall.d[!min_dim][!e] = offset_inner;
							unsigned const num_lights2((c_hall.get_sz_dim(min_dim) > 0.25*s_hall.get_sz_dim(!min_dim)) ? 2 : 1); // 2 lights if it's long enough
							add_room(c_hall, part_id, num_lights2, 1, 0); // add conn hallway as room
							rooms.back().mark_open_wall_dim(min_dim); // adjacent to primary and secondary hallway on each side
							cube_t exclude(c_hall);
							exclude.d[min_dim][!d] += 1.25*dsign*doorway_width; // expand out into the main hallway to ensure there's space to enter this hallway
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
							float const nom_wall_edge(rooms_range[e]), min_room_width(1.6*doorway_width);
							float shifted_offset_inner(shift_val_to_not_intersect_window(*p, nom_wall_edge, hspace, window_border, !min_dim));

							if (fabs(shifted_offset_inner - offset_outer) < min_room_width) { // end room is too narrow, shift in other dir
								shifted_offset_inner = shift_val_to_not_intersect_window(*p, nom_wall_edge, hspace, window_border, !min_dim, 1);
								shift_edges_mode[e]  = 1; // opposite shift

								if (cent_room_width - fabs(shifted_offset_inner - nom_wall_edge) - 0.5*window_hspacing[!min_dim] < min_room_width) {
									// center room is now too narrow (assuming half window width shift on opposite wall), split the window
									cube_t unshifted_wall(room_outer);
									set_wall_width(unshifted_wall, nom_wall_edge, wall_thick, !min_dim); // side of building
									unshifted_wall.expand_in_dim(min_dim, wall_thick); // extend outside of building to include windows
									interior->split_window_walls.push_back(unshifted_wall);
									shifted_offset_inner = nom_wall_edge;
									shift_edges_mode[e]  = 2; // no shift
								}
							}
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
						bool const doors_in_corners(put_doors_in_corners && cent_room_width > 2.0*doorway_width);
						bool const door_side(doors_in_corners ? rgen.rand_bool() : 0); // side to push the door when in a corner of the room; consistent per wall
						float room_pos(rooms_range[0]);

						for (unsigned n = 0; n <= num_cent_rooms; ++n) { // center rooms, inside and outside
							// this iteration includes both end walls of outer rooms, so it must check for shift_edges_other_dir at both ends
							bool const is_first(n == 0), is_last(n+1 == num_cent_rooms);
							int const sem_lo(is_first ? shift_edges_mode[0] : ((n == num_cent_rooms) ? shift_edges_mode[1] : 0));
							float const start_pos(shift_val_to_not_intersect_window(*p, room_pos, hspace, window_border, !min_dim, sem_lo));

							if (n < num_cent_rooms) { // add a room
								room_pos += cent_room_width; // move to next row
								int const sem_hi(is_last ? shift_edges_mode[1] : 0);
								float const next_pos(shift_val_to_not_intersect_window(*p, room_pos, hspace, window_border, !min_dim, sem_hi)); // start_pos of next n
								room_outer.d[!min_dim][0] = start_pos;
								room_outer.d[!min_dim][1] = next_pos;
								room_inner.d[!min_dim][0] = (is_first ? (rooms_range[0] + wall_half_thick) : start_pos); // first inner room is relative to the sec hallway
								room_inner.d[!min_dim][1] = (is_last  ? (rooms_range[1] - wall_half_thick) : next_pos ); // last inner room is relative to sec hallway
								add_room(room_outer, part_id, 1, 0, 1); // office
								add_room(room_inner, part_id, 1, 0, 1); // office or bathroom
								bool const is_bathroom(n == bathroom_ix);

								if (is_bathroom) {rooms.back().assign_all_to(RTYPE_BATH);} // make it a bathroom (windowless)
								else if (is_ground_floor) { // previous inner room room (windowless) for the ground floor
									bool const next_to_br(n+1 == bathroom_ix || n == bathroom_ix+1);
									(next_to_br ? utility_room_cands : special_room_cands).push_back(rooms.size() - 1); // utility rooms must be next to the bathroom
								}
								// add doors to 2-3 walls
								float door_pos(0.5f*(start_pos + next_pos)); // centered by default
								cube_t *to_split[3] = {&long_swall, &short_swall, &main_wall}; // outer, inner, inner

								for (unsigned k = 0; k < num_doors_inner_rooms; ++k) {
									if (doors_in_corners) {
										cube_t const &room(k ? room_inner : room_outer);
										door_pos = (door_side ? (room.d[!min_dim][1] - door_to_wall_dist) : (room.d[!min_dim][0] + door_to_wall_dist));
									}
									float const lo_pos(door_pos - doorway_hwidth), hi_pos(door_pos + doorway_hwidth);
									insert_door_in_wall_and_add_seg(*to_split[k], lo_pos, hi_pos, !min_dim, (d^(k&1)), 1, is_bathroom);
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
					float const sh_len(room_width), min_room_width(1.5*doorway_width);
					int const num_sec_halls(num_rooms/2); // round down if odd
					int windows_per_room_od((windows_per_side_od & 1) ? 1 : 2), rooms_per_side(0); // rooms along each sec hall
					float room_sub_width(0.0);

					while (1) { // increase the number of windows per room until room is large enough
						rooms_per_side = max(windows_per_side_od/windows_per_room_od, 1);
						room_sub_width = sh_len/rooms_per_side;
						if (room_sub_width > min_room_width) break; // done
						++windows_per_room_od;
					}
					float const sh_spacing(hall_len/num_sec_halls - sh_width), end_spacing(0.5*sh_spacing); // half spacing at both ends
					assert(sh_spacing > min_wall_len); // I'm not sure if this can fail or what we should do in that case - use fewer secondary hallways?
					float room_start(p->d[!min_dim][0]), wall_pos(room_start + end_spacing); // first sec hall wall pos
					int const num_offices(4*rooms_per_side*num_sec_halls), rand_val(rgen.rand());
					bool const add_div_wall_door(is_office_bldg() && rooms_per_side > 2 && (rand_val & 255)); // 50% of the time when conditions are met
					int const min_num_door_stacks(num_offices + add_div_wall_door*2*(num_sec_halls-1)); // one per office + extras
					auto &split_walls(hall_walls);
					room_walls .reserve(4*(rooms_per_side+1)*num_sec_halls + 2*(num_sec_halls-1)); // walls with doors + room dividers
					split_walls.reserve(2*rooms_per_side*(num_sec_halls+1));
					rooms      .reserve(num_offices + 2*num_sec_halls + 1); // offices + sec hallways + pri hallway
					interior->door_stacks.reserve(min_num_door_stacks);
					interior->doors      .reserve(num_doors_per_stack*min_num_door_stacks);
					interior->exclusion  .reserve(2*num_sec_halls);
					unsigned const bathroom_ix((num_sec_halls <= 2) ? 0 : (rand_val%(num_sec_halls-1))); // place bathrooms in rooms near the central hallway

					for (int i = 0; i <= num_sec_halls; ++i) { // actually iterates over the number of room blocks between halls (num halls + 1)
						// shift the position of the hallway walls to avoid intersecting windows;
						// this may make halls either larger (contain a window) or smaller (fit between two windows);
						// as long as the space between windows is large enough (> doorway_width) this should be okay
						float target_hall_end_pos(wall_pos + sh_width); // use post-shifted wall for better fit to orig hall width
						float hall_start_pos(shift_val_to_not_intersect_window(*p, wall_pos,            hspace, window_border, !min_dim));
						float hall_end_pos  (shift_val_to_not_intersect_window(*p, target_hall_end_pos, hspace, window_border, !min_dim));
						float sec_hall_width(hall_end_pos - hall_start_pos);
						bool const div_room(i > 0 && i < num_sec_halls), add_sec_hall(i < num_sec_halls);

						// handle hallways that are too wide or too narrow; do narrow check last because it's more important
						if (add_sec_hall && sec_hall_width > max_hall_width) { // hall is too wide, try shifting start pos
							float const target_hall_start_pos(hall_start_pos - (max_hall_width - sec_hall_width));
							hall_start_pos = shift_val_to_not_intersect_window(*p, target_hall_start_pos, hspace, window_border, !min_dim);
							sec_hall_width = hall_end_pos - hall_start_pos;
							
							if (sec_hall_width > 1.01*max_hall_width) { // still too large; try shifting to the other side of the window
								target_hall_end_pos = hall_end_pos + (max_hall_width - sec_hall_width);
								hall_end_pos   = shift_val_to_not_intersect_window(*p, target_hall_end_pos, hspace, window_border, !min_dim);
								sec_hall_width = hall_end_pos - hall_start_pos;
							}
						}
						if (add_sec_hall && sec_hall_width < min_hall_width) { // hall is too narrow, try shifting start pos
							float const target_hall_start_pos(hall_start_pos - (min_hall_width - sec_hall_width));
							hall_start_pos = shift_val_to_not_intersect_window(*p, target_hall_start_pos, hspace, window_border, !min_dim);
							sec_hall_width = hall_end_pos - hall_start_pos;

							if (sec_hall_width < 0.99*min_hall_width) { // still too small; try shifting to the other side of the window
								hall_start_pos = shift_val_to_not_intersect_window(*p, target_hall_start_pos, hspace, window_border, !min_dim, 1); // shift_edges_mode=1
								sec_hall_width = hall_end_pos - hall_start_pos;
							}
						}
						for (unsigned d = 0; d < 2; ++d) { // left, right of main hall
							float const dsign(d ? -1.0 : 1.0);
							bool const doors_in_corners(put_doors_in_corners && room_sub_width > 2.0*doorway_width);
							// should the door always be close to the hallway (shorter walking path), or close to the exterior wall (no dead end hallway)?
							bool const door_side(d); // closer to exterior wall
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
							cube_t sep_walls[2], first_room;

							if (add_sec_hall) { // add sec hall
								cube_t s_hall(*p);
								s_hall.d[ min_dim][!d] = hall.d[min_dim][d]; // ends at main hall
								s_hall.d[!min_dim][ 0] = hall_start_pos;
								s_hall.d[!min_dim][ 1] = hall_end_pos;
								unsigned const num_lights(min(4U, max(2U, unsigned(0.35*s_hall.get_sz_dim(min_dim)/window_vspacing))));
								add_room(s_hall, part_id, num_lights, 1, 0); // add sec hallway as room
								rooms.back().mark_open_wall(min_dim, !d); // adjacent to main hallway on this side
								cube_t exclude(s_hall);
								exclude.d[min_dim][!d] += 1.25*dsign*doorway_width; // expand out into the main hallway to ensure there's space to enter this hallway
								interior->exclusion.push_back(exclude); // excluded from placing stairs and elevators
								
								for (unsigned dir = 0; dir < 2; ++dir) { // add walls between hall and rooms on each side
									sep_walls[dir] = div_wall; // copy z and min_dim from div_wall
									sep_walls[dir].d[min_dim][!d] += dsign*wall_half_thick; // extend to meet the wall in the other dim
									set_wall_width(sep_walls[dir], s_hall.d[!min_dim][dir], wall_half_thick, !min_dim);
								}
							}
							for (int r = 0; r < rooms_per_side; ++r) { // add office rooms, working from the windows inside
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
								bool const is_br_aisle((unsigned)i == (bathroom_ix+1));
								bool const is_bathroom(r+1 == rooms_per_side && is_br_aisle); // bathroom must be an interior/windowless room
								if (r == 0 && !is_bathroom) {first_room = room;} // used for splitting divider wall

								if (is_bathroom) {rooms.back().assign_all_to(RTYPE_BATH);}
								else if (is_ground_floor && div_room && r > 0) { // windowless room on ground floor, not at ext wall
									bool const next_to_br((rooms_per_side <= 2) ? ((unsigned)i == (bathroom_ix) || (unsigned)i == (bathroom_ix+2)) : is_br_aisle);
									(next_to_br ? utility_room_cands : special_room_cands).push_back(rooms.size() - 1); // utility rooms must be next to the bathroom
								}
								if (add_sec_hall) { // add doorways + doors
									float door_pos(0.0);
									if (!doors_in_corners) {door_pos = 0.5f*(room_split_pos + next_split_pos);} // centered	
									else {door_pos = (door_side ? (room.d[min_dim][1] - door_to_wall_dist) : (room.d[min_dim][0] + door_to_wall_dist));}
									float const lo_pos(door_pos - doorway_hwidth), hi_pos(door_pos + doorway_hwidth);

									for (unsigned dir = 0; dir < 2; ++dir) {
										insert_door_in_wall_and_add_seg(sep_walls[dir], lo_pos, hi_pos, min_dim, dir, !d, is_bathroom); // bathroom doors are always open
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
							if (div_room) { // add a divider wall
								if (add_div_wall_door && !first_room.is_all_zeros()) { // maybe split wall and add a door connecting opposite rooms
									bool const open_dir(rooms.size() & 1); // pseudo random
									float const door_pos(first_room.get_center_dim(min_dim)), lo_pos(door_pos - doorway_hwidth), hi_pos(door_pos + doorway_hwidth);
									insert_door_in_wall_and_add_seg(div_wall, lo_pos, hi_pos, min_dim, open_dir);
								}
								room_walls.push_back(div_wall);
							}
						} // for d
						room_start = hall_end_pos;
						wall_pos  += sh_spacing + sh_width;
						min_eq(wall_pos, p->d[!min_dim][1]); // clamp to far side of building to handle the final row of rooms
					} // for i
				} // end secondary hallways
			} // end multiple hallways case

			else { // single main hallway
				// Note: these reserves may be unnecessary now that buildings commonly have stacked parts and basements
				room_walls.reserve(2*(num_rooms-1));
				hall_walls.reserve(2*(num_rooms+1));
				cube_t rwall(*p); // copy from part; shared zvals, but X/Y will be overwritten per wall
				create_wall(rwall, !min_dim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing); // room walls, create first wall
				//bool const doors_in_corners(put_doors_in_corners && !apt_or_hotel && room_len > 2.0*doorway_width); // this case is more complex and has not been implemented
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
					float const adj_door_val(doorway_vals[s ? (2*num_rooms-2) : 1]);
					hwall.d[!min_dim][ s] = place_area.d[!min_dim][s]; // end at the wall
					hwall.d[!min_dim][!s] = doorway_vals[s ? (2*num_rooms-1) : 0] = adj_door_val + (s ? 1.0f : -1.0f)*doorway_width; // end at first doorway

					for (unsigned d = 0; d < 2; ++d) {
						set_wall_width(hwall, hall_wall_pos[d], wall_half_thick, min_dim);
						hall_walls.push_back(hwall);
					}
				} // for s
				// add rooms and doors
				rooms.reserve(2*num_rooms + 1); // two rows of rooms + hallway
				interior->door_stacks.reserve(2*num_rooms);
				interior->doors      .reserve(num_doors_per_stack*2*num_rooms);
				float const wall_end(p->d[!min_dim][1]);
				float pos(p->d[!min_dim][0]);

				for (int i = 0; i < num_rooms; ++i) {
					// Note: it would probably be simpler to create one wall and cut doorways into it like what's done in the secondary hallways case above
					float const next_pos((i == num_rooms-1) ? wall_end : (pos + room_len)); // handle different size of last room; must end exactly at the building edge

					for (unsigned d = 0; d < 2; ++d) { // left, right sides of hallway
						cube_t c(*p); // copy zvals and exterior wall pos
						c.d[ min_dim][!d] = hall_wall_pos[d];
						c.d[!min_dim][ 0] = pos;
						c.d[!min_dim][ 1] = next_pos;
						// apartments and hotels have utility rooms on the first floor, but not on corner/end units
						if (apt_or_hotel && i > 0 && i+1 < num_rooms) {utility_room_cands.push_back(rooms.size());}
						add_room(c, part_id, 1, 0, 1); // office or bathroom; no utility rooms for now (except for apartments or hotels) since these rooms tend to be large
						bool const is_bathroom(i == num_rooms/2 && !apt_or_hotel);
						if (is_bathroom) {rooms.back().assign_all_to(RTYPE_BATH);} // assign the middle room to be a bathroom
						door_t door(c, min_dim, d); // copy zvals and wall pos
						clip_wall_to_ceil_floor(door, fc_thick);
						door.d[min_dim][d] = hall_wall_pos[d]; // set to zero area at hallway
						for (unsigned e = 0; e < 2; ++e) {door.d[!min_dim][e] = doorway_vals[2*i+e];} // set doorway left and right edges
						add_interior_door(door, is_bathroom);
						if (apt_or_hotel) {divide_last_room_into_apt_or_hotel(i, num_rooms, num_windows, windows_per_room, windows_per_side_od, !min_dim, d, rgen);}
					} // for d
					pos = next_pos;
				} // for i
			} // end single main hallway case

			if ((num_windows_od & 1) && num_hall_windows > 1.5) {
				// walls may intersect windows, so flag them so that we can cover the wall ends with extra window trim
				for (auto w = hall_walls.begin()+hall_walls_start; w != hall_walls.end(); ++w) {
					for (unsigned d = 0; d < 2; ++d) {
						if (fabs(w->d[!min_dim][d] - p->d[!min_dim][d]) > wall_thick) continue; // doesn't meet the exterior wall
						cube_t sw_wall(*w);
						set_wall_width(sw_wall, w->get_center_dim(min_dim), wall_thick, min_dim); // side of building
						sw_wall.expand_in_dim(!min_dim, wall_thick); // extend outside of building to include windows
						interior->split_window_walls.push_back(sw_wall);
					}
				} // for w
			}
			if (is_hospital()) {add_hospital_bathrooms(rooms_start, rgen);}
			add_room(hall, part_id, 3, 1, 0); // add primary hallway as room with 3+ lights
			if (has_sec_hallways) {rooms.back().mark_open_wall_dim(min_dim);} // flag primary hallway as open on sides if there are secondary hallways
			if (is_ground_floor || pri_hall.is_all_zeros()) {pri_hall = hall;} // assign to primary hallway if on first floor and hasn't yet been assigned
			no_split_walls_this_part = 1;
		} // end use_hallway

		else if (is_prison()) {
			divide_part_into_jail_cells(*p, part_id, rgen);
			add_part_sep_walls(p, place_area, rooms_start, must_split);
		}
		else { // generate random walls using recursive 2D slices
			float const min_wall_len2(0.85*min_wall_len); // a somewhat shorter value that applies to some tests (but not wall_split_thresh)
			bool const no_walls(min(p->dx(), p->dy()) < min_wall_len2); // not enough space to add a room (chimney, porch support, garage, shed, etc.)
			float const min_split_len(max(global_building_params.wall_split_thresh, 1.0f)*min_wall_len);
			assert(to_split.empty());
			if (no_walls) {add_room(*p, part_id, 1);} // add entire part as a room
			else {to_split.emplace_back(*p);} // seed room is entire part, no door
			bool is_first_split(1);
			point part_door_open_dir_tp(p->get_cube_center()); // used to determine in which direction doors open; updated base on central hallway
			
			if (first_part) { // reserve walls/rooms/doors - take a guess at the correct size
				unsigned const num_doors_est(4*real_num_parts + has_basement());
				for (unsigned d = 0; d < 2; ++d) {interior->walls[d].reserve(8*real_num_parts);}
				rooms.reserve(8*real_num_parts + has_sec_bldg()); // two rows of rooms + optional hallway
				interior->door_stacks.reserve(num_doors_est);
				interior->doors.reserve(num_doors_per_stack*num_doors_est);
			}
			unsigned const ds_start(interior->door_stacks.size());
			// see if we have a skylight to work with
			vect_cube_t part_skylights; // should be at most size 1 with currently skylight addition code
			bool added_whirl(0);

			for (cube_t const &skylight : skylights) {
				if (skylight.intersects(*p)) {part_skylights.push_back(skylight);}
			}
			while (!to_split.empty()) { // recursively split rooms until all rooms are too small to split or there are no more valid splits
				split_cube_t const c(to_split.back()); // Note: non-const because door_lo/door_hi is modified during T-junction insert
				to_split.pop_back();
				vector3d const csz(c.get_size());
				bool wall_dim(0); // which dim the room is split by
				if      (csz.y > min_wall_len2 && csz.x > 1.25*csz.y) {wall_dim = 0;} // split long room in x
				else if (csz.x > min_wall_len2 && csz.y > 1.25*csz.x) {wall_dim = 1;} // split long room in y
				else {wall_dim = rgen.rand_bool();} // choose a random split dim for nearly square rooms

				if (min(csz.x, csz.y) < min_wall_len2) {
					add_room(c, part_id, 1);
					continue; // not enough space to add a wall
				}
				// add central hallway if wall/hall len is at least enough to place 2-3 rooms; 50% chance if part is a basement; houses and small office buildings
				bool const is_basement_room(has_basement() && part_id == (unsigned)basement_part_ix);
				bool const add_hallway(is_first_split && csz[!wall_dim] > 1.5*min_split_len && (!is_basement_room || rgen.rand_bool()));

				if (!add_hallway && min(csz.x, csz.y) > 3.0*min_wall_len && rgen.rand_bool()) { // create a whirl pattern
					cube_t cr(c), sr[4] = {c, c, c, c}; // center room and 4 side rooms (in CCW order)
					bool is_valid(1);

					for (unsigned d = 0; d < 2 && is_valid; ++d) { // {x, y}
						bool const on_edge(c.d[!d][0] == p->d[!d][0] || c.d[!d][1] == p->d[!d][1]); // at edge of the building
						float const cc(c.get_center_dim(d));

						for (unsigned e = 0; e < 2 && is_valid; ++e) { // {lo, hi}
							float const lo(e ? (cc+0.5*min_wall_len) : (c.d[d][0]+min_wall_len2)), hi(e ? (c.d[d][1]-min_wall_len2) : (cc-0.5*min_wall_len));
							is_valid = 0;

							for (unsigned n = 0; n < 10; ++n) {
								cr.d[d][e] = rgen.rand_uniform(lo, hi);
								if (on_edge && is_val_inside_window(*p, d, cr.d[d][e], window_hspacing[d], window_border)) continue; // check for windows
								is_valid = 1; break;
							}
						} // for e
					} // for d
					if (is_valid) {
						for (cube_t const &sl : part_skylights) { // center room must contain all skylights
							if (!cr.contains_cube_xy(sl)) {is_valid = 0; break;}
						}
					}
					if (is_valid) { // walls added successfully
						bool const ccw(rgen.rand_bool()); // CCW whirl
						sr[0].x2() = sr[1].x1() = cr.d[0][ ccw];
						sr[2].x1() = sr[3].x2() = cr.d[0][!ccw];
						sr[0].y2() = sr[3].y1() = cr.d[1][!ccw];
						sr[1].y2() = sr[2].y1() = cr.d[1][ ccw];
						float const edges[4] = {cr.y1(), cr.x2(), cr.y2(), cr.x1()};
					
						for (unsigned d = 0; d < 4; ++d) {
							bool const wdim(bool(d&1) ^ ccw);
							float const wall_pos(edges[(d + (ccw ? 0 : 3)) & 3]);
							cube_t wall(sr[d]), wall2;
							create_wall(wall, wdim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing);
							float const doorway_pos(cube_rand_side_pos(cr, !wdim, 0.25, doorway_width, rgen, 1)); // for_door=1
							bool const open_dir(wall_pos > part_door_open_dir_tp[wdim]); // doors open away from the building center
							remove_section_from_cube_and_add_door(wall, wall2, (doorway_pos - doorway_hwidth), (doorway_pos + doorway_hwidth), !wdim, open_dir);
							interior->walls[wdim].push_back(wall );
							interior->walls[wdim].push_back(wall2);
						} // for d
						for (unsigned r = 0; r < 5; ++r) {
							cube_t &room(r ? sr[r-1] : cr);

							for (unsigned d = 0; d < 2; ++d) { // subtract interior wall thickness from rooms
								if (room.d[d][0] > c.d[d][0]) {room.d[d][0] += wall_half_thick;}
								if (room.d[d][1] < c.d[d][1]) {room.d[d][1] -= wall_half_thick;}
							}
							if (max(room.dx(), room.dy()) > min_split_len) {to_split.push_back(room);} // split sub-room
							else {add_room(room, part_id);} // add whole sub-room
						}
						added_whirl = 1;
						continue;
					}
				} // end whirl pattern
				float const min_dist_abs(1.5*doorway_width + wall_thick);
				bool const on_edge(c.d[!wall_dim][0] == p->d[!wall_dim][0] || c.d[!wall_dim][1] == p->d[!wall_dim][1]); // at edge of the building - walls don't intersect windows
				// create a wall to split up this room
				// what about allowing adjacent rooms not separated by a wall?
				// this would work for connected public rooms (living, dining, kitchen), but would limit our ability to later assign these rooms as bed, bath, etc.
				float wall_pos(0.0), min_dist_param(0.25);
				bool pos_valid(0);
				cube_t wall, wall2; // copy from cube; shared zvals, but X/Y will be overwritten per wall

				for (unsigned num = 0; num < 20; ++num) { // 20 tries to choose a wall pos that's not inside a window or skylight or intersecting a wall
					wall_pos = cube_rand_side_pos(c, wall_dim, min_dist_param, min_dist_abs, rgen, 0); // for_door=0
					if (on_edge && is_val_inside_window(*p, wall_dim, wall_pos, window_hspacing[wall_dim], window_border)) continue; // try a new wall_pos
					if (c.bad_pos(wall_pos, wall_dim)) continue; // intersects doorway from prev wall, try a new wall_pos
					wall = c;
					create_wall(wall, wall_dim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing);
					
					if (has_bcube_int_xy(wall, part_skylights)) { // check for skylight intersection
						if (num >= 10 && min_dist_param > 0.1) {min_dist_param -= 0.05;} // allow walls further from the center to avoid the skylight
						continue;
					}
					pos_valid = 1;

					if (added_whirl) { // check for wall intersection with center room door of whirl
						cube_t test_cube(wall);
						test_cube.expand_by_xy(wall_thick);

						for (auto ds = interior->door_stacks.begin()+ds_start; ds != interior->door_stacks.end(); ++ds) {
							if (ds->get_true_bcube().intersects(test_cube)) {pos_valid = 0; break;}
						}
					}
					if (pos_valid) break; // done, keep wall_pos
				} // for num
				if (!pos_valid) { // no valid pos, skip this split
					add_room(c, part_id, 1);
					continue;
				}
				// insert a doorway into the wall
				float const doorway_pos(cube_rand_side_pos(c, !wall_dim, 0.25, doorway_width, rgen, 1)); // for_door=1
				float const lo_pos(doorway_pos - doorway_hwidth), hi_pos(doorway_pos + doorway_hwidth);
				bool const open_dir(wall_pos > part_door_open_dir_tp[wall_dim]); // doors open away from the building center
				remove_section_from_cube_and_add_door(wall, wall2, lo_pos, hi_pos, !wall_dim, open_dir);
				auto &walls(interior->walls[wall_dim]);
				walls.push_back(wall );
				walls.push_back(wall2);
				float door_lo[2] = {lo_pos, lo_pos}, door_hi[2] = {hi_pos, hi_pos}; // passed to next split step to avoid placing a wall that intersects this doorway

				if (add_hallway) { // add central hallway
					// maybe create a hallway: create another split parallel to this one offset a bit and make the room in between a hallway
					bool const dir((wall_pos - c.d[wall_dim][0]) < (c.d[wall_dim][1] - wall_pos)); // further part edge
					float other_wall_pos(0.0);
					bool other_pos_valid(0);

					for (unsigned num = 0; num < 40; ++num) { // 40 tries to choose a wall pos that's not inside a window
						float const upper_sz((num < 20) ? 2.4 : 3.0); // allow wider hallways if the first 20 placements fail
						other_wall_pos = wall_pos + ((dir ? 1.0 : -1.0)*rgen.rand_uniform(1.6, upper_sz)*doorway_width); // opposite edge of hallway
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
						pref_conn_to = hall; // prefer to connect doors to this room
						pref_conn_to.expand_by_xy(2.0*wall_thick);
						// add other wall parts and doorway, with a different random doorway pos
						cube_t o_wall1(c), o_wall2;
						create_wall(o_wall1, wall_dim, other_wall_pos, fc_thick, wall_half_thick, wall_edge_spacing);
						float const door_pos(cube_rand_side_pos(c, !wall_dim, 0.25, doorway_width, rgen, 1)); // for_door=1
						float const o_lo_pos(door_pos - doorway_hwidth), o_hi_pos(door_pos + doorway_hwidth);
						remove_section_from_cube_and_add_door(o_wall1, o_wall2, o_lo_pos, o_hi_pos, !wall_dim, !open_dir); // opens in other dir
						walls.push_back(o_wall1);
						walls.push_back(o_wall2);
						wall.d[wall_dim][dir] = o_wall1.d[wall_dim][dir]; // make original wall the entire width of the hall so that the room on the other side doesn't overlap the hall
						door_lo[dir] = o_lo_pos;
						door_hi[dir] = o_hi_pos;
						part_door_open_dir_tp[wall_dim] = 0.5f*(wall_pos + other_wall_pos); // use the center of the hallway so that doors open into the rooms
						if (real_num_parts == 1) {bldg_door_open_dir_tp[wall_dim] = part_door_open_dir_tp[wall_dim];} // this is also the test point for the entire building
					}
				}
				// split into two smaller rooms; ensure we can split at least once per dim
				bool const do_split(csz[wall_dim] > min(min_split_len, max(min_wall_len, 0.9f*bcube.get_sz_dim(wall_dim))));

				for (unsigned d = 0; d < 2; ++d) { // still have space to split in other dim, add the two parts to the stack
					split_cube_t c_sub(c);
					c_sub.d[wall_dim][d] = wall.d[wall_dim][!d]; // clip to wall pos
					c_sub.door_lo[!wall_dim][d] = door_lo[!d] - wall_half_thick; // set new door pos in this dim (keep door pos in other dim, if set)
					c_sub.door_hi[!wall_dim][d] = door_hi[!d] + wall_half_thick;
					if (do_split) {to_split.push_back(c_sub);} else {add_room(c_sub, part_id, 1);} // leaf case (unsplit), add a new room
				}
				is_first_split = 0;
			} // end while()
			add_part_sep_walls(p, place_area, rooms_start, must_split);
		} // end wall placement
		add_ceilings_floors_stairs(rgen, *p, hall, part_id, num_floors, rooms_start, use_hallway, first_part_this_stack, window_hspacing, window_border, is_industrial_part);
		assign_special_room_types(utility_room_cands, special_room_cands, doors_start, rgen); // for rooms with hallways
		
		if (no_split_walls_this_part) { // don't split any walls added up to this point (for fixed floorplans with primary hallways)
			for (unsigned d = 0; d < 2; ++d) {first_wall_to_split[d] = interior->walls[d].size();}
		}
	} // for p (parts)
	if (has_sec_bldg()) { // add garage/shed floor and ceiling
		assert(parts_end < parts.end());
		cube_t const &garage(*parts_end);
		cube_t C(garage);
		C.z2() = garage.z1() + fc_thick;
		interior->floors.push_back(C);
		C.z2() = garage.z2();
		C.z1() = C.z2() - fc_thick;
		interior->ceilings.push_back(C);
		add_room(garage, (parts_end - parts.begin()), 1, 0, 0, 1); // is_sec_bldg=1
		rooms.back().set_no_geom();
	}
	// attempt to cut extra doorways into long walls if there's space to produce a more connected floorplan
	for (unsigned d = 0; d < 2; ++d) { // x,y: dim in which the wall partitions the room (wall runs in dim !d)
		auto &walls(interior->walls[d]);
		auto const &perp_walls(interior->walls[!d]);

		// Note: iteration will include newly added all segments to recursively split long walls
		for (unsigned w = first_wall_to_split[d]; w < walls.size(); ++w) { // skip hallway walls
			bool pref_split(must_split[d] & (1ULL << (w & 63)));
			float stairs_elev_pad(doorway_width);

			for (unsigned nsplits = 0; nsplits < 4; ++nsplits) { // at most 4 splits
				cube_t &wall(walls[w]); // take a reference here because a prev iteration push_back() may have invalidated it
				if (is_industrial() && wall.z1() >= ground_floor_z1) break; // don't split warehouse/factory walls as they already have doors
				float const len(wall.get_sz_dim(!d)), min_split_len((pref_split ? 0.5 : 1.5)*min_wall_len); // = 2.0/6.0 * doorway_width
				if (len < min_split_len) break; // not long enough to split - done
				float const min_dist_abs(min(1.5f*doorway_width, max(0.5f*doorway_width, 0.5f*min_split_len)));
				bool const pref_conn(pref_conn_to.contains_cube(wall)); // house hallway, etc.
				// walls currently don't run along the inside of exterior building walls, so we don't need to handle that case yet
				bool was_split(0);

				for (unsigned ntries = 0; ntries < (pref_split ? 40U : 10U); ++ntries) { // choose random doorway positions and check against perp_walls for occlusion
					float const doorway_pos(cube_rand_side_pos(wall, !d, 0.2, min_dist_abs, rgen, 1)); // for_door=1
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

							for (room_t const &r : rooms) {
								if (r.z1() >= wall.z2() || r.z2() <= wall.z1()) continue; // no overlap in Z

								for (unsigned e = 0; e < 2; ++e) { // check both directions from the wall (front and back)
									// skip wall edges co-incident with room edges where the entire wall is contained in the room; this can happen at part boundaries
									if (wall.d[d][e] == r.d[d][e] && wall.d[d][!e] < r.d[d][1] && wall.d[d][!e] > r.d[d][0]) continue;
									if (wall.d[d][e] <  r.d[d][0] || wall.d[d][ e] > r.d[d][1]) continue; // wall not inside room in dim d/dir e
									if (lo[s] > r.d[!d][0] && hi[s] < r.d[!d][1]) {contained[e] = 1;} // entire wall contained in span of room
								}
							} // for r
							if (contained[0] && contained[1]) {valid = 0; break;} // wall seg contained in rooms on both sides => two doors in same wall between rooms => drop
						} // for s
					}
					if (!valid) continue;
					cube_t cand(wall); // sub-section of wall that will become a doorway
					cand.d[!d][0] = lo_pos; cand.d[!d][1] = hi_pos;
					// if the wall must be split, and we've failed after 20 tries, reduce the stairs and elevator padding; this happens to a building near the player start
					if (pref_split && ntries >= 20 && (ntries%10) == 0) {stairs_elev_pad *= 0.5;} // 1.0 for n=0-19, 0.5 for n=20-29, 0.25 for n=30-39
					bool const elevators_only(pref_split && ntries > 20); // allow blocking stairs if there's no other way to insert a door
					int  const no_check_enter_exit((ntries > 5) ? (pref_conn ? 2 : 1) : 0); // no stairs enter/exit pad after first 5 tries, no enter/exit at all if pref_conn
					if (interior->is_blocked_by_stairs_or_elevator(cand, stairs_elev_pad, elevators_only, no_check_enter_exit)) continue; // should we assert !pref_split?
					// check for other nearby doorways
					bool pos_valid(1);

					for (door_stack_t const &ds : interior->door_stacks) {
						if (ds.z1() < cand.z2() && ds.z2() > cand.z1() && ds.get_open_door_path_bcube().intersects(cand)) {pos_valid = 0; break;}
					}
					if (!pos_valid) continue;
					bool open_dir(wall.get_center_dim(d) > bldg_door_open_dir_tp[d]); // doors open away from the building center
					if (is_prison() && !is_prison_door_valid(cand, d, open_dir)) continue; // check if too close to jail cells
					was_split = 1;

					// Note: this code doesn't work for multiple reasons but is left in for reference in case I figure this out later
					if (0 && is_house && wall.z1() >= ground_floor_z1) { // not the basement
						// find rooms on each side, check if one is adjacent to an exterior door (likely the living room or entryway),
						// then remove the entire wall on the first floor rather than adding a door; but this doesn't work because doors haven't been placed yet
						// but this will make the floorplans different, is that legal? only for single floor rooms? and it will require updates to AI navigation logic
						cube_t door_area(wall);
						door_area.d[!d][0] = lo_pos; // clip to door range
						door_area.d[!d][1] = hi_pos;
						door_area.expand_in_dim(d, wall_thick); // make sure the wall overlaps the rooms on both sides
						room_t rl, rh; // room to the {low, high} sides

						for (room_t const &r : rooms) {
							if (r.intersects_no_adj(door_area)) {((r.get_center_dim(d) < wall.get_center_dim(d)) ? rl : rh) = r;}
						}
						if (!rl.is_hallway && !rh.is_hallway && !rl.is_all_zeros() && !rh.is_all_zeros()) { // don't connect a hallway
							if (is_room_adjacent_to_ext_door(rl) || is_room_adjacent_to_ext_door(rh)) { // one room must be adjacent to an exterior door
								float const cut_lo(max(rl.d[!d][0], rh.d[!d][0])), cut_hi(min(rl.d[!d][1], rh.d[!d][1]));

								if (cut_lo <= lo_pos && cut_hi >= hi_pos) { // check for valid range; can fail occasionally when the same wall is split multiple times
									cube_t wall2(wall);
									wall .d[!d][1] = cut_lo; // lo side
									wall2.d[!d][0] = cut_hi; // hi side
								
									if (wall2.is_strictly_normalized()) {
										if (!wall.is_strictly_normalized()) {wall = wall2;}
										else {interior->walls[d].push_back(wall2);}
									}
									else if (!wall.is_strictly_normalized()) {
										wall = walls.back(); // remove this wall
										walls.pop_back();
										nsplits = 0; // reset outer loop iteration (well, will start at 1 out of 4 anyway; could set to -1)
									}
									break; // done
								}
							}
						}
					}
					bool const jail_door(is_prison());
					insert_door_in_wall_and_add_seg(wall, lo_pos, hi_pos, !d, open_dir, 0, 0, 0, 0, jail_door); // high_side=is_br=unlocked=closed=0; modifies wall
					break;
				} // for ntries
				if (!was_split) break; // no more splits
				pref_split = 0; // already split, no longer preferred
			} // for nsplits
		} // for w
	} // for d
	reverse_door_hinges_if_needed();

	// add stairs to connect together stacked parts for office buildings; must be done last after all walls/ceilings/floors have been assigned
	for (unsigned p = 0; p < real_num_parts; ++p) {connect_stacked_parts_with_stairs(rgen, parts[p], p);}
	if (has_parking_garage || is_parking()) {add_parking_garage_ramp(rgen);}

	if (is_parking()) {
		add_parking_structure_entrance(rgen);
		add_parking_structure_bathroom(rgen);
	}
	// furnace logic
	if (has_attic()) {
		if (has_basement()) {interior->furnace_type = (rgen.rand_bool() ? FTYPE_BASEMENT : FTYPE_ATTIC);} // both attic and basement: place furnace in one at random
		else {interior->furnace_type = FTYPE_ATTIC;} // attic only: place furnace in attic
	}
	else if (has_basement()) {interior->furnace_type = FTYPE_BASEMENT;} // basement only: place furnace in basement (at least if it's a house)
	// else no furnace
} // end gen_interior_int()

void building_t::add_part_sep_walls(vect_cube_t::const_iterator p, cube_t const &place_area, unsigned rooms_start, uint64_t must_split[2]) {
	// insert walls to split up parts into rectangular rooms
	float const wall_thick(get_wall_thickness()), fc_thick(get_fc_thickness());
	auto parts_end(get_real_parts_end());
	float const min_split_wall_len(0.5*get_min_wall_len()); // allow a shorter than normal wall because these walls have higher priority

	if (min(p->dx(), p->dy()) < min_split_wall_len) { // too small, skip
		has_small_part = 1;
		return;
	}
	for (auto p2 = parts.begin(); p2 != parts_end; ++p2) {
		if (p2 == p) continue; // skip self
		if (min(p2->dx(), p2->dy()) < min_split_wall_len) {has_small_part = 1; continue;} // too small, skip

		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				float const val(p->d[!dim][dir]);
				if (p2->d[!dim][!dir] != val) continue; // not adjacent
				if (p2->z1() >= p->z2() || p2->z2() <= p->z1()) continue; // no overlap in Z
				if (p2->d[dim][0] > p->d[dim][0] || p2->d[dim][1] < p->d[dim][1]) continue; // not contained in dim (don't have to worry about Z-shaped case)

				if (p2->d[dim][0] == p->d[dim][0] && p2->d[dim][1] == p->d[dim][1]) { // same xy values, must only vary in z
					if (p2->z2() < p->z2()) continue; // add wall only on one side (arbitrary)
					if (has_complex_floorplan && p2->z2() == p->z2() && dir == 0) continue; // two abutting walls in complex floorplan: skip the one with dir==0
				}
				cube_t wall;
				wall.z1() = max(p->z1(), p2->z1()) + fc_thick; // shared Z range
				wall.z2() = min(p->z2(), p2->z2()) - fc_thick;

				for (unsigned d = 0; d < 2; ++d) {
					if (p->d[dim][d] != p2->d[dim][d]) { // corner between the two parts
						wall.d[dim][d] = p->d[dim][d]; // shorter part side with no gap
						//wall.d[dim][d] -= (d ? 1.0 : -1.0)*0.2*wall_edge_spacing; // reduced gap (but still causes a noticeable split in the interior wall and trim)
					}
					else {
						wall.d[dim][d] = place_area.d[dim][d]; // shorter part side with slight offset
					}
				} // for d
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
}

void building_t::assign_special_room_types(vector<unsigned> &utility_room_cands, vector<unsigned> &special_room_cands, unsigned doors_start, rand_gen_t &rgen) {
	vector<room_t> &rooms(interior->rooms);
	vector<unsigned> conf_rooms;
	
	if (!is_residential() && !is_industrial()) { // office building, hospital, school, etc.
		// try to assign conference room(s); this may be overwritten with a special room on the ground floor;
		// conference rooms are currently interior and windowless so that we can have glass walls without blending problems with windows
		float max_area(0.0);
		vector<unsigned> conf_room_cands;

		for (unsigned d = 0; d < 2; ++d) { // conference rooms can optionally be next to a bathroom
			for (unsigned room_ix : (d ? special_room_cands : utility_room_cands)) {
				room_t &room(get_room(room_ix));
				if (room.has_stairs || room.has_elevator || room.is_hallway || interior->is_blocked_by_stairs_or_elevator(room)) continue; // not a valid conference room
				if (get_conf_room_table_length_width(get_walkable_room_bounds(room)) == vector2d()) continue; // can't fit a conference room table
				conf_room_cands.push_back(room_ix);
				max_eq(max_area, room.get_area_xy());
			}
		} // for d
		for (unsigned room_ix : conf_room_cands) {
			if (rooms[room_ix].get_area_xy() < 0.99*max_area) continue; // skip non-max area rooms
			if (!rgen.rand_bool()) continue; // 50% chance of conference room
			rooms[room_ix].assign_all_to(RTYPE_CONF, 1); // locked=1
			conf_rooms.push_back(room_ix);
		}
	}
	// add special ground floor room types
	vector<unsigned> special_room_types; // placed in this priority order
	special_room_types.push_back(RTYPE_UTILITY);
	if (is_hotel()       || is_hospital()) {special_room_types.push_back(RTYPE_LAUNDRY );}
	if (is_office_bldg() || is_hospital()) {special_room_types.push_back(RTYPE_SECURITY);} // add for some schools as well?
	if (is_office_bldg() || is_hospital() || is_school()) {special_room_types.push_back(RTYPE_SERVER);}

	for (unsigned rtype : special_room_types) {
		// allow special rooms to use utility room cands when all utitility rooms have been assigned and there are no special rooms (no secondary hallways)
		bool const is_utility(rtype == RTYPE_UTILITY);
		vector<unsigned> &room_cands((is_utility || special_room_cands.empty()) ? utility_room_cands : special_room_cands);

		for (unsigned n = 0; n < (is_utility ? MAX_OFFICE_UTILITY_ROOMS : 1); ++n) {
			if (room_cands.empty()) break; // no more rooms to assign
			unsigned &room_ix(room_cands[rgen.rand() % room_cands.size()]);
			room_t &room(get_room(room_ix));
			room.assign_to(rtype, 0, 1); // assign this room on floor 0; locked=1
			bool const ensure_locked(0); // probably should be locked, but unlocked makes these rooms easier to explore
			ensure_doors_to_room_are_closed(room, doors_start, ensure_locked);

			if (room.is_apt_or_hotel_room()) { // make other rooms in this unit the same type; should be RTYPE_UTILITY or RTYPE_LAUNDRY
				for (auto r = rooms.begin()+room_ix+1; r != rooms.end() && r->unit_id == room.unit_id; ++r) {r->assign_to(rtype, 0, 1);} // locked=1
			}
			room_ix = room_cands.back();
			room_cands.pop_back(); // remove this room from consideration
		} // for n
	} // for rtype
	for (unsigned room_ix : conf_rooms) {add_conference_room_window(room_ix);}
}

void building_t::add_conference_room_window(unsigned room_ix) {
	room_t &room(get_room(room_ix));
	assert(!room.is_hallway);
	// find the longest wall segment shared between the room and an adjacent hallway
	float const wall_thickness(get_wall_thickness()), wall_hthick(0.5*wall_thickness), fc_thick(get_fc_thickness()), door_width(get_doorway_width());
	float const window_edge_space(0.25*door_width), min_window_width(0.5*door_width), min_wall_len(min_window_width + 2.0*window_edge_space);
	bool best_dim(0);
	float best_lo(0.0), best_hi(0.0);
	cube_t *best_wall(nullptr);
	room_t *conn_hall(nullptr);

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			cube_t wall_area(room);
			set_wall_width(wall_area, room.d[dim][dir], wall_thickness, dim);
			wall_area.expand_in_dim(!dim, -wall_hthick); // shrink to room interior
			wall_area.expand_in_z(-fc_thick); // shrink off floor and ceiling to avoid picking up adjacent rooms

			for (room_t &r : interior->rooms) {
				if (!r.is_hallway || !r.intersects(wall_area)) continue;
				cube_t shared_area(wall_area); // shared wall between conference room and hallway
				max_eq(shared_area.d[!dim][0], r.d[!dim][0]);
				min_eq(shared_area.d[!dim][1], r.d[!dim][1]);

				for (cube_t &w : interior->walls[dim]) {
					if (!w.intersects(shared_area)) continue;
					float const lo(max(shared_area.d[!dim][0], w.d[!dim][0])), hi(min(shared_area.d[!dim][1], w.d[!dim][1]));
					if ((hi - lo) <= max((best_hi - best_lo), min_wall_len)) continue; // too short
					best_dim = dim; best_lo = lo; best_hi = hi; best_wall = &w; conn_hall = &r;
				}
			} // for r
		} // for dir
	} // for dim
	if (best_wall != nullptr) {
		assert(conn_hall != nullptr);
		cube_t &wall(*best_wall);
		float const wind_lo(best_lo + window_edge_space), wind_hi(best_hi - window_edge_space);
		cube_t wall2(wall), window(wall);
		wall .d[!best_dim][1] = window.d[!best_dim][0] = wind_lo; // low  edge of window
		wall2.d[!best_dim][0] = window.d[!best_dim][1] = wind_hi; // high edge of window
		interior->walls[best_dim].push_back(wall2);
		room.set_interior_window();
		conn_hall->set_interior_window();

		if (room.get_room_type(0) != RTYPE_CONF) { // re-add ground floor wall if re-assigned as a utility or special room
			cube_t short_wall(window);
			short_wall.z2() = window.z1() = short_wall.z1() + get_window_vspace(); // only floor tall
			interior->walls[best_dim].push_back(short_wall);
		}
		interior->int_windows.emplace_back(window, room_ix);
	}
}

void building_t::divide_last_room_into_apt_or_hotel(unsigned room_row_ix, unsigned hall_num_rooms, unsigned tot_num_windows,
	unsigned windows_per_room, unsigned windows_per_room_side, bool hall_dim, bool hall_dir, rand_gen_t &rgen) // hall_dim is also window dim
{
	assert(interior);
	assert(hall_num_rooms  > 0 && room_row_ix < hall_num_rooms);
	assert(tot_num_windows > 0 && windows_per_room > 0 && windows_per_room <= tot_num_windows);
	assert(is_apt_or_hotel());
	vector<room_t> &rooms(interior->rooms);
	// the current room and door stack are the last ones that were added; note that adding new rooms or doors below will invalidate these references
	assert(!rooms.empty());
	assert(!interior->door_stacks.empty());
	room_t &room(rooms.back());
	unsigned const new_rooms_start(rooms.size()-1); // including current room
	//if (is_room_adjacent_to_ext_door(room)) {} // can't check, walkway hasn't been added yet
	door_stack_t const &ds(interior->door_stacks.back());
	bool at_lo_end(room_row_ix == 0), at_hi_end(room_row_ix == hall_num_rooms-1);
	unsigned num_windows(windows_per_room);
	if (is_apartment() && num_windows == 1) {btype = BTYPE_HOTEL;} // too small for apartment; need at least two windows for an apartment
	if (at_hi_end) {num_windows += (tot_num_windows % windows_per_room);} // last room is larger with more windows; hopefully only one more
	unsigned const part_id(room.part_id);
	float const door_width(ds.get_width()), door_hwidth(0.5*door_width); // newly added doors will be this width as well
	float const wall_thickness(get_wall_thickness()), wall_half_thick(0.5*wall_thickness), door_to_wall_min_space(2.0*wall_thickness), fc_thick(get_fc_thickness());
	float const door_center(ds.get_center_dim(hall_dim)), room_center(room.get_center_dim(hall_dim));
	bool const lg_door_side(door_center < room_center); // more space on this side of door
	vect_cube_t &hall_para_walls(interior->walls[!hall_dim]), &hall_perp_walls(interior->walls[hall_dim]);
	bool make_small_apt(is_apartment() && num_windows < 3), make_three_room(is_hotel());
	room.is_office = 0; // but room.office_floorplan remains 1
	assert(ds.first_door_ix < interior->doors.size());

	// make sure hotel room doors auto close and apartment doors start closed
	for (auto d = interior->doors.begin()+ds.first_door_ix; d != interior->doors.end(); ++d) {
		if (is_hotel()) {d->make_auto_close();} else {d->open = 0; d->open_amt = 0.0;}
	}
	if (make_three_room || make_small_apt) { // 3-4 rooms (hotel or small apartment)
		// divide into entryway/living room by door, bedroom by window, and bathroom
		// divide room in short dim; side by window becomes bedroom, divide other side, part connected to door is living room/entryway, other part is bathroom
		float bed_lb_split_pos(room.get_center_dim(!hall_dim));
		float liv_bath_split_pos(ds.d[hall_dim][lg_door_side] + (lg_door_side ? 1.0 : -1.0)*door_to_wall_min_space); // door goes to living room, which is larger
		if (lg_door_side) {max_eq(liv_bath_split_pos, room_center);} else {min_eq(liv_bath_split_pos, room_center);} // large side must be at least half the room width

		if (has_int_windows() && (at_lo_end || at_hi_end) && windows_per_room_side > 1 && (windows_per_room_side & 1)) { // shift bed_lb_split_pos to prevent window intersection
			float const window_h_space(room.get_sz_dim(!hall_dim)/windows_per_room_side);
			bed_lb_split_pos += (hall_dir ? -1.0 : 1.0)*(0.5*(1.0 - get_window_h_border())*window_h_space + get_wind_trim_thick()); // shift near edge of window frame
		}
		// add new rooms; rooms tile exactly and have no space for walls
		cube_t bed(room), lb(room);
		bed.d[!hall_dim][!hall_dir] = lb.d[!hall_dim][hall_dir] = bed_lb_split_pos;
		cube_t living(lb), bath(lb);
		bath.d[hall_dim][!lg_door_side] = living.d[hall_dim][lg_door_side] = liv_bath_split_pos;
		room.copy_from(living); // ext_sides doesn't change
		calc_room_ext_sides(room); // update since ext_sides may have changed
		room.assign_all_to(is_hotel() ? RTYPE_COMMON : RTYPE_LIVING); // public first; common room is similar to living room but without the table
		room.set_is_entryway();
		unsigned const living_rid(new_rooms_start), bed_rid(add_room(bed, part_id)), bath_rid(add_room(bath, part_id));
		get_room(bed_rid ).assign_all_to(RTYPE_BED);
		get_room(bath_rid).assign_all_to(make_small_apt ? RTYPE_KITCHEN : RTYPE_BATH); // small apartment uses a kitchen rather than a bathroom for this room
		// add interior walls and doors; all doors are unlocked
		bool const no_div_wall_or_door(is_hotel() && rgen.rand_float() < 0.75); // open wall 75% of the time; should this be consistent per building?
		float const living_center(living.get_center_dim(hall_dim));
		float bath_center(bath.get_center_dim(!hall_dim));
		cube_t bed_wall(bed), bath_wall(bath); // if no_div_wall_or_door=1, bedroom wall only borders the bathroom
		clip_wall_to_ceil_floor(bed_wall,  fc_thick);
		clip_wall_to_ceil_floor(bath_wall, fc_thick);
		set_wall_width(bed_wall,  bed_lb_split_pos,   wall_half_thick, !hall_dim);
		set_wall_width(bath_wall, liv_bath_split_pos, wall_half_thick,  hall_dim);

		if (no_div_wall_or_door) { // make bedroom/living room split an open wall
			cube_t open_wall(bed_wall);
			// split the bedroom wall at the living room/bathroom boundary, with half a wall width of overlap to fill the gap
			bed_wall .d[ hall_dim][!lg_door_side] = open_wall.d[hall_dim][lg_door_side] = liv_bath_split_pos + (lg_door_side ? -1.0 : 1.0)*0.5*wall_thickness;
			bath_wall.d[!hall_dim][hall_dir] += (hall_dir ? 1.0 : -1.0)*0.5*wall_thickness; // half a wall width of overlap to fill the gap
			interior->open_walls.push_back(open_wall);
			get_room(bed_rid   ).mark_open_wall(!hall_dim, !hall_dir);
			get_room(living_rid).mark_open_wall(!hall_dim,  hall_dir);
		}
		else { // add living room to bedroom door
			if (bed_wall.get_sz_dim( hall_dim) <= door_width) {cout << "bedroom too small: " << bed .str() << endl;} // error?
			else {insert_door_in_wall_and_add_seg(bed_wall, living_center-door_hwidth, living_center+door_hwidth, hall_dim, hall_dir, 0, 0, 1);} // opens into bedroom
		}
		if (make_small_apt && rgen.rand_float() < 0.75) { // make living room/kitchen split an open wall
			interior->open_walls.push_back(bath_wall); // kitchen wall
			get_room(bath_rid  ).mark_open_wall(hall_dim, !lg_door_side);
			get_room(living_rid).mark_open_wall(hall_dim,  lg_door_side);
		}
		else {
			if (bath_wall.get_sz_dim(!hall_dim) <= door_width) {cout << "bathroom too small: " << bath.str() << endl;} // error?
			else {insert_door_in_wall_and_add_seg(bath_wall, bath_center-door_hwidth, bath_center+door_hwidth, !hall_dim, lg_door_side, 0, 1, 1);} // opens into bathroom
			hall_perp_walls.push_back(bath_wall);
		}
		hall_para_walls.push_back(bed_wall );

		if (make_small_apt) { // 4 rooms (apartment); bathroom becomes kitchen, and bathroom is added to bedroom
			float bed_bath_split_pos(liv_bath_split_pos);

			if (has_int_windows()) { // prevent walls from intersecting windows
				cube_t const &part(parts[part_id]);
				float const window_hspacing(part.get_sz_dim(hall_dim)/get_num_windows_on_side(part, hall_dim));
				bed_bath_split_pos = shift_val_to_not_intersect_window(part, bed_bath_split_pos, window_hspacing, get_window_h_border(), hall_dim);
				set_wall_width(bath_wall, bed_bath_split_pos, wall_half_thick, hall_dim); // update wall pos in case it was moved
			}
			bath = bed;
			bath.d[hall_dim][!lg_door_side] = bed.d[hall_dim][lg_door_side] = bed_bath_split_pos;
			get_room(bed_rid).copy_from(bed); // update with smaller bedroom
			calc_room_ext_sides(get_room(bed_rid)); // update since ext_sides may have changed
			get_room(add_room(bath, part_id)).assign_all_to(RTYPE_BATH);
			// add wall and door
			bath_center = bath.get_center_dim(!hall_dim);
			for (unsigned d = 0; d < 2; ++d) {bath_wall.d[!hall_dim][d] = bath.d[!hall_dim][d];} // move to new bathroom pos
			if (bath_wall.get_sz_dim(!hall_dim) <= door_width) {cout << "bathroom 2 too small: " << bath.str() << endl;} // error?
			else {insert_door_in_wall_and_add_seg(bath_wall, bath_center-door_hwidth, bath_center+door_hwidth, !hall_dim, lg_door_side, 0, 1, 1);} // opens into bathroom
			hall_perp_walls.push_back(bath_wall);
		}
	}
	else if (is_apartment()) { // 5 rooms (apartment)
		// entryway connected to front door, or could be living room if much longer in hallway dim
		// bedroom; must have at least one window
		// living room; must have at least one window
		// bathroom
		// kitchen
		float front_back_split_pos(room.get_center_dim(!hall_dim) + (hall_dir ? -1.0 : 1.0)*0.1*room.get_sz_dim(!hall_dim)); // bed+living are slightly larger
		float living_bed_split_pos(room_center);
		float const entry_door_space_signed((lg_door_side ? 1.0 : -1.0)*(0.5*door_width + door_to_wall_min_space));
		float kitchen_split_pos(ds.d[hall_dim][ lg_door_side] + entry_door_space_signed);
		float bath_split_pos   (ds.d[hall_dim][!lg_door_side] - entry_door_space_signed);
		float const shift(entry_door_space_signed * min(0.5f, 0.5f*fabs(door_center - room_center)/door_width)); // move split closer to the center of the room
		kitchen_split_pos += shift;
		bath_split_pos    += shift;

		if (has_int_windows()) { // prevent walls from intersecting windows
			if (num_windows & 1) { // odd number of windows; shift living_bed_split_pos to not intersect a window
				float const window_h_space(room.get_sz_dim(hall_dim)/num_windows);
				living_bed_split_pos += (lg_door_side ? -1.0 : 1.0)*(0.5*(1.0 - get_window_h_border())*window_h_space + get_wind_trim_thick()); // shift near edge of bedroom window frame
			}
			if (at_lo_end || at_hi_end) { // check side windows
				cube_t const &part(parts[part_id]);
				float const window_hspacing(part.get_sz_dim(!hall_dim)/get_num_windows_on_side(part, !hall_dim));
				front_back_split_pos = shift_val_to_not_intersect_window(part, front_back_split_pos, window_hspacing, get_window_h_border(), !hall_dim);
			}
		}
		// add new rooms; rooms tile exactly and have no space for walls
		cube_t const room_area(room);
		cube_t bed(room), living(room), bath(room), kitchen(room), entry(room);
		bed.d[!hall_dim][!hall_dir] = living.d[!hall_dim][!hall_dir] = /* window side */
			bath.d[!hall_dim][hall_dir] = kitchen.d[!hall_dim][hall_dir] = entry.d[!hall_dim][hall_dir] = front_back_split_pos; // entry side
		living .d[hall_dim][!lg_door_side] = bed  .d[hall_dim][ lg_door_side] = living_bed_split_pos; // living room should be larger, or equal size
		kitchen.d[hall_dim][!lg_door_side] = entry.d[hall_dim][ lg_door_side] = kitchen_split_pos;
		bath   .d[hall_dim][ lg_door_side] = entry.d[hall_dim][!lg_door_side] = bath_split_pos;
		room.copy_from(entry);
		calc_room_ext_sides(room); // update since ext_sides may have changed
		room.assign_all_to(RTYPE_ENTRY);
		room.set_is_entryway();
		unsigned const bed_rid(add_room(bed, part_id)), bath_rid(add_room(bath, part_id)), kitchen_rid(add_room(kitchen, part_id)), living_rid(add_room(living, part_id));
		get_room(bed_rid    ).assign_all_to(RTYPE_BED    );
		get_room(bath_rid   ).assign_all_to(RTYPE_BATH   );
		get_room(kitchen_rid).assign_all_to(RTYPE_KITCHEN);
		get_room(living_rid ).assign_all_to(RTYPE_LIVING );
		// add interior walls
		cube_t fb_wall(room_area), lb_wall(bed), ke_wall(kitchen), be_wall(bath); // front-back, living room-bedroom, kitchen-entryway, bathroom-entryway
		clip_wall_to_ceil_floor(fb_wall, fc_thick);
		clip_wall_to_ceil_floor(lb_wall, fc_thick);
		clip_wall_to_ceil_floor(ke_wall, fc_thick);
		clip_wall_to_ceil_floor(be_wall, fc_thick);
		set_wall_width(fb_wall, front_back_split_pos, wall_half_thick, !hall_dim);
		set_wall_width(lb_wall, living_bed_split_pos, wall_half_thick,  hall_dim);
		set_wall_width(ke_wall, kitchen_split_pos,    wall_half_thick,  hall_dim);
		set_wall_width(be_wall, bath_split_pos,       wall_half_thick,  hall_dim);
		// add doors; all are unlocked
		float const front_center  (ke_wall.get_center_dim(!hall_dim)), back_center(lb_wall.get_center_dim(!hall_dim));
		float const kitchen_center(kitchen.get_center_dim( hall_dim)), bath_center(bath   .get_center_dim( hall_dim));
		float const le_lo(min(kitchen_split_pos, living_bed_split_pos)), le_hi(max(kitchen_split_pos, living_bed_split_pos));
		float const be_lo(min(bath_split_pos,    living_bed_split_pos)), be_hi(max(bath_split_pos,    living_bed_split_pos));
		bool const keep_high_side(!lg_door_side), conn_bed_to_bath(rgen.rand_bool()); // randomly add one of two bathroom doors
		insert_door_in_wall_and_add_seg(lb_wall, back_center -door_hwidth, back_center +door_hwidth, !hall_dim, !lg_door_side, 0, 0, 1); // opens into bedroom
		insert_door_in_wall_and_add_seg(ke_wall, front_center-door_hwidth, front_center+door_hwidth, !hall_dim,  lg_door_side, 0, 0, 1); // opens into kitchen
		// entry-bathroom door is optional; opens into bathroom
		if (!conn_bed_to_bath) {insert_door_in_wall_and_add_seg(be_wall, front_center-door_hwidth, front_center+door_hwidth, !hall_dim, !lg_door_side, 0, 1, 1);}
		insert_door_in_wall_and_add_seg(fb_wall, kitchen_center-door_hwidth, kitchen_center+door_hwidth, hall_dim, hall_dir,  keep_high_side, 0, 1); // opens into living room

		if (le_hi - le_lo > 1.2*door_width) { // have space for entryway-living room door
			float const le_center(0.5*(le_lo + le_hi));
			insert_door_in_wall_and_add_seg(fb_wall, le_center-door_hwidth, le_center+door_hwidth, hall_dim, hall_dir,  keep_high_side, 0, 1); // opens into living room
		}
		if (be_hi - be_lo > 1.2*door_width) { // have space for entryway-bedroom room door
			float const be_center(0.5*(be_lo + be_hi));
			insert_door_in_wall_and_add_seg(fb_wall, be_center-door_hwidth, be_center+door_hwidth, hall_dim, hall_dir,  keep_high_side, 0, 1); // opens into bedroom
		}
		// bedroom-bathroom door is optional; opens into bathroom
		if (conn_bed_to_bath) {insert_door_in_wall_and_add_seg(fb_wall, bath_center-door_hwidth, bath_center+door_hwidth, hall_dim, hall_dir, !keep_high_side, 0, 1);}
		hall_para_walls.push_back(fb_wall);
		hall_perp_walls.push_back(lb_wall);
		hall_perp_walls.push_back(ke_wall);
		hall_perp_walls.push_back(be_wall);
	}
	else {assert(0);} // unsupported building type

	for (auto r = rooms.begin()+new_rooms_start; r != rooms.end(); ++r) {
		r->set_office_floorplan();
		r->unit_id = next_unit_id;
	}
	++next_unit_id;
}

bool building_t::maybe_assign_interior_garage(bool &gdim, bool &gdir) {
	if (interior == nullptr) return 0;
	vector<room_t> &rooms(interior->rooms);
	unsigned const num_rooms(rooms.size());
	if (!is_house || has_sec_bldg() || num_rooms < 8U) return 0; // no garage for this building (including the case where we have a shed)

	if (has_basement()) { // count non-basement rooms
		unsigned non_basement_rooms(0);
		for (auto r = rooms.begin(); r != rooms.end(); ++r) {non_basement_rooms += (!(r->z1() < ground_floor_z1));}
		if (non_basement_rooms < 8) return 0; // not enough non-basement rooms
	}
	rand_gen_t rgen;
	rgen.set_state(mat_ix+1, num_rooms+1);
	unsigned const room_start(rgen.rand() % num_rooms);
	float const wall_thickness(get_wall_thickness());
	bool const pref_street_dim(street_dir ? ((street_dir-1) >> 1) : 0), pref_street_dir(street_dir ? ((street_dir-1)&1) : 0);
	int best_room(-1);
	float best_score(0.0);

	for (unsigned rix = 0; rix < num_rooms; ++rix) {
		unsigned const cur_room((rix + room_start) % num_rooms);
		room_t &r(rooms[cur_room]);
		if (r.has_stairs_on_floor(0) || r.is_hallway || r.z1() != ground_floor_z1) continue;
		if (has_basement() && r.part_id == (int)basement_part_ix) continue; // skip basement rooms
		if (get_part_for_room(r).contains_cube_xy_no_adj(r)) continue; // skip interior rooms
		if (r.get_room_type(0) != RTYPE_NOTSET) continue; // already assigned
		if (r.is_single_floor) continue; // no tall ceiling rooms; those should be used for living rooms; likely not needed because these aren't assigned yet
		bool const dim(street_dir ? pref_street_dim : (r.dx() < r.dy())); // use larger dim unless there's a preference
		if (r.d[!dim][0] != bcube.d[!dim][0] && r.d[!dim][1] != bcube.d[!dim][1]) continue; // require other dim at either side of the building
		cube_t room_interior(r);
		room_interior.expand_by_xy(-wall_thickness); // shrink slightly to allow a bit of extra padding at walls
		if (!car_can_fit(room_interior)) continue; // too small
		float const length(r.get_sz_dim(dim)), width(r.get_sz_dim(!dim)), ar(length/width);
		if (ar > 2.6 || ar < 1.1) continue; // aspect ratio too high or too low
		bool const pref_dir(street_dir ? pref_street_dir : rgen.rand_bool());

		for (unsigned d = 0; d < 2; ++d) {
			bool const dir(bool(d) ^ pref_dir);
			if (r.d[dim][dir] != bcube.d[dim][dir])   continue; // exterior walls on the building bcube only
			if (street_dir && dir != pref_street_dir) continue; // wrong dir
			float const score(-(float)count_num_int_doors(r) - fabs(ar - 1.75f)); // prefer rooms with fewer interior doors, and aspect ratios around 1.75
			if (best_room == -1 || score > best_score) {best_room = cur_room; best_score = score; gdim = dim; gdir = dir;}
			break; // no need to try the other dir
		} // for d
	} // for r
	if (best_room < 0) return 0; // failed
	assert((unsigned)best_room < rooms.size());
	room_t &room(rooms[best_room]);
	room.assign_to(RTYPE_GARAGE, 0, 1); // first floor, locked=1
	cube_t room_exp(room);
	room_exp.expand_by_xy(wall_thickness);

	// ensure all doors connected to the garage open outward from the garage so that they won't clip through cars (though it's generally the opposite for real houses)
	for (door_stack_t &ds : interior->door_stacks) {
		if (!ds.intersects(room_exp)) continue;

		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				cube_t wall(room);
				wall.d[dim][!dir] = wall.d[dim][dir]; // shrink to zero width
				wall.expand_in_dim(dim, wall_thickness); // expand wall outward
				if (!ds.intersects(wall)) continue;
				ds.open_dir = dir; // update the door stack
				assert(ds.first_door_ix < interior->doors.size());

				for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
					door_t &door(interior->doors[dix]);
					if (!ds.is_same_stack(door)) break; // moved to a different stack, done
					door.open_dir = dir; // update the doors themselves
				}
			} // for dir
		} // for dim
	} // for d
	interior->garage_room = best_room;
	return 1;
}

void find_and_merge_with_landing(vector<landing_t> &landings, cube_t const &stairs, stairs_shape &sshape, unsigned num_floors, bool is_above) {
	for (landing_t &landing : landings) { // process and maybe update all landings for this stairwell
		if (!landing.intersects(stairs)) continue; // wrong landing; only one landing stack should intersect the stairs
		if (landing.shape == SHAPE_WALLED) {sshape = SHAPE_WALLED; landing.shape = SHAPE_WALLED_SIDES;} // shift bottom wall down a floor
		// if the parking garage is multiple levels, exclude the back wall so that we can connect the lower level(s) with this same stairwell
		if (sshape == SHAPE_WALLED && num_floors > 1) {sshape = SHAPE_WALLED_SIDES;}
		if (is_above) {landing.is_at_top = 0;} // landing below is no longer at the top
		return;
	}
	assert(0); // should never get here
}

void building_t::add_ceilings_floors_stairs(rand_gen_t &rgen, cube_t const &part, cube_t const &hall, unsigned part_ix, unsigned num_floors,
	unsigned rooms_start, bool use_hallway, bool first_part_this_stack, float window_hspacing[2], float window_border, bool is_single_floor)
{
	// increase floor thickness if !is_house? but then we would probably have to increase the space between floors as well, which involves changing the texture scale
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const doorway_width(get_nominal_doorway_width()), wall_thickness(get_wall_thickness());
	float min_ewidth(1.5*doorway_width), ewidth(min_ewidth); // square; for elevators
	float z(part.z1());
	cube_t stairs_cut, elevator_cut;
	bool stairs_dim(0), bend_dir(0), stairs_have_railing(1), extended_from_above(0);
	bool stairs_against_wall[2] = {0, 0};
	stairs_shape sshape(SHAPE_STRAIGHT); // straight by default
	bool must_add_stairs(first_part_this_stack || (is_prison() && interior->stairwells.empty()));
	bool const is_basement_room((int)part_ix == basement_part_ix);
	int force_stairs_dir(2); // 2=unset
	unsigned stairs_ix(0);

	if (has_complex_floorplan) { // the tallest part of a complex floorplan building contains the stairs; should be the last part, unless there's a basement
		must_add_stairs = 1;

		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (p->z2() > part.z2()) {must_add_stairs = 0; break;} // don't need stairs if another part is taller
		}
	}
	// add stairwells and elevator shafts
	if (num_floors == 1) {} // no need for stairs or elevator
	else if (use_hallway) { // part is the hallway cube
		assert(!interior->rooms.empty());
		room_t &room(interior->rooms.back()); // hallway is always the last room to be added
		bool const long_dim(hall.dx() < hall.dy());
		// U-shape if there's enough room in width; also required for two floor retail; increase the width of both the stairs and elevator
		if (room.get_sz_dim(!long_dim) > 6.0*doorway_width || has_tall_retail()) {sshape = SHAPE_U; ewidth *= 1.6;}
		else {sshape = SHAPE_WALLED_SIDES;} // walled sides to meet fire codes
		cube_t stairs(hall); // start as hallway
		// add elevator(s)
		bool const add_side_elevator(!is_apt_or_hotel() && rgen.rand_bool()); // as opposed to central elevators; not valid for fixed floorplans of apartments and hotels
		float const hall_len(room.get_sz_dim(long_dim)), stairs_len(4.0*doorway_width), ehwidth(0.5*ewidth);
		unsigned const num_elevators((hall_len > 10.0*ewidth) ? 2 : 1); // two elevators if there's space; will only add 1 if there's a side elevator
		unsigned const room_ix(interior->rooms.size() - 1);
		bool elevator_added(0);
		if (interior->landings.empty()) {interior->landings.reserve(num_elevators + (num_floors-1));} // lower bound

		if (add_side_elevator) {
			bool const dir(rgen.rand_bool());
			float const wall_pos(room.d[!long_dim][dir] + (dir ? -1.0 : 1.0)*0.5*wall_thickness); // side of hallway
			elevator_t elevator(room, room_ix, !long_dim, !dir, 0, 1); // elevator shaft; at_edge=0, interior_room=1 (considered interior-enough)
			elevator.d[!long_dim][!dir] = wall_pos;
			elevator.d[!long_dim][ dir] = wall_pos + (dir ? 1.0 : -1.0)*ewidth; // cuts into rooms
			float const stairs_spacing(0.5*stairs_len), end_spacing(max((stairs_len + stairs_spacing), 0.2f*hall_len));
			float const place_lo(room.d[long_dim][0] + end_spacing), place_hi(room.d[long_dim][1] - end_spacing);
			cube_t const room_bc(interior->has_sec_hallways ? cube_t() : room); // for door coll test; use an empty room because the door may open to a side hallway
			
			if (place_lo < place_hi) { // hallway should be long enough
				vect_cube_t avoid; // avoid bathrooms; what about utility rooms?

				for (room_t const &r : interior->rooms) {
					if (!is_bathroom(r.get_room_type(0))) continue; // not a bathroom; assumes bathrooms are assigned to all floors, and checking the first floor is correct
					if (r.intersects(part)) {avoid.push_back(r);}
				}
				for (unsigned n = 0; n < 100; ++n) { // try 100 random positions along the wall
					set_wall_width(elevator, rgen.rand_uniform(place_lo, place_hi), ehwidth, long_dim);
					assert(part.contains_cube(elevator));
					if (has_bcube_int(elevator, avoid))               continue;
					if (is_cube_close_to_doorway(elevator, room_bc))  continue; // try again
					if (has_bcube_int(elevator, interior->exclusion)) continue; // try again
					if (check_skylight_intersection(elevator))        continue; // check skylights; is this necessary?
					add_or_extend_elevator(elevator, 1);
					cube_t const sub_cube(elevator.get_bcube_padded(window_vspacing)); // trim is excluded
					for (unsigned d = 0; d < 2; ++d) {subtract_cube_from_cubes(sub_cube, interior->walls[d], nullptr, 1);}
					float const room_center(room.get_center_dim(long_dim));
					set_wall_width(stairs, room_center, 0.5*stairs_len, long_dim); // start with stairs centered
					cube_t avoid(elevator);
					avoid.expand_in_dim(long_dim, stairs_spacing);

					if (stairs.d[long_dim][0] <= avoid.d[long_dim][1] && avoid.d[long_dim][0] <= stairs.d[long_dim][1]) { // overlapping or adjacent in long_dim
						bool const stairs_dir(avoid.get_center_dim(long_dim) < room_center); // set stairs to the side with more space
						stairs.translate_dim(long_dim, (avoid.d[long_dim][stairs_dir] - stairs.d[long_dim][!stairs_dir])); // shift stairs away from elevator (req for subtract)
					}
					room.has_elevator = 1;
					elevator_cut      = elevator;
					elevator_added    = 1;
					break; // done
				} // for N
			}
		}
		if (num_elevators > 0 && !elevator_added) { // add central elevator(s) if no side elevator was added
			point center(room.get_cube_center());
			// choose elevator dir on first elevator, and reuse this dir for later parts;
			// increases symmetry and possibly improves the chances of connecting elevators vertically through multiple stacked parts
			if (interior->elevators.empty()) {interior->elevator_dir = rgen.rand_bool();}
			float const center_shift(0.125*hall_len*(interior->elevator_dir ? -1.0 : 1.0));
			center[long_dim] += center_shift; // make elevator off-center
			elevator_t elevator(room, room_ix, long_dim, rgen.rand_bool(), 0, 1); // elevator shaft; at_edge=0, interior_room=1 (considered interior-enough)
			for (unsigned d = 0; d < 2; ++d) {set_wall_width(elevator, center[d], ehwidth, d);}

			if (num_elevators == 1) {add_or_extend_elevator(elevator, 1);} // single elevator
			else { // double back-to-back elevators
				assert(num_elevators == 2);
				elevator.expand_in_dim(long_dim, min(0.5f*min_ewidth, ehwidth)); // increase the depth up to 2x

				for (unsigned e = 0; e < 2; ++e) {
					elevator_t E(elevator);
					E.d[long_dim][bool(e) ^ E.dir ^ 1] = center[long_dim]; // back-to-back
					E.dir        ^=     bool(e); // facing opposite directions
					E.adj_pair_ix = 1 + bool(e); // flag so that we don't include this in our hallway elevator count
					unsigned const elevator_ix(interior->elevators.size());
					if (e) {assert(elevator_ix > 0); E.adj_elevator_ix = elevator_ix-1;} // adjacent to previous elevator
					else {E.adj_elevator_ix = elevator_ix+1;} // adjacent to next elevator
					add_or_extend_elevator(E, 1);
				}
			}
			room.has_elevator = 1;
			elevator_cut      = elevator;
			stairs.translate_dim(long_dim, -center_shift); // shift stairs in the opposite direction
		}
		// always add stairs
		for (unsigned dim = 0; dim < 2; ++dim) { // shrink in XY
			bool const is_step_dim(bool(dim) == long_dim); // same orientation as the hallway
			float shrink(stairs.get_sz_dim(dim) - (is_step_dim ? stairs_len : 0.9*ewidth)); // set max size of stairs opening
			stairs.expand_in_dim(dim, -0.5*shrink); // centered in the hallway
		}
		room.has_stairs = 255; // stairs on all floors
		room.set_has_center_stairs();
		stairs_cut = stairs;
		stairs_dim = long_dim;
	}
	else if (is_basement_room && part.contains_cube_xy(pri_hall) && can_extend_stairs_to_pg(stairs_ix)) { // multi-floor parking garage case
		stairwell_t &s(interior->stairwells[stairs_ix]);
		s.extends_below = 1;
		// copy fields from these stairs and extend down
		stairs_cut = s;
		stairs_dim = s.dim;
		force_stairs_dir    = s.dir;
		extended_from_above = 1;
		sshape       = s.shape;
		// assume we can extend the existing hallway stairs downward
		copy_zvals(stairs_cut, part);
		room_t &room(interior->rooms.back()); // should be the last room
		room.has_stairs = 255; // stairs on all floors
		cube_t stairs_bot(s);
		stairs_bot.z2() = stairs_bot.z1() + window_vspacing; // limit to bottom landing
		find_and_merge_with_landing(interior->landings, stairs_bot, sshape, num_floors, 0); // merge with bottom landing; is_above=0
	}
	// only add stairs to first part of a house unless we haven't added stairs yet, or if it's the top floor of a stacked part
	else if (!is_house || interior->stairwells.empty() || (first_part_this_stack && part.z1() > ground_floor_z1)) {
		bool const is_parking_str(is_parking() && !is_basement_room);
		bool add_elevator(0);
		// sometimes add an elevator to building parts, but not the first part in a stack (to guarantee we have at least one set of stairs)
		// it might not be possible to place an elevator a part with no interior rooms, but that should be okay, because some other part will still have stairs
		// do we need support for multiple floor cutouts stairs + elevator in this case as well?
		if (!is_house && !is_parking_str && !must_add_stairs) {
			float const elevator_prob[4] = {0.9, 0.5, 0.3, 0.1}; // higher chance of adding an elevator if there are more existing elevators
			add_elevator = (rgen.rand_float() < elevator_prob[min(interior->elevators.size(), size_t(3))]);
		}
		unsigned const rooms_end(interior->rooms.size()), num_avail_rooms(rooms_end - rooms_start);
		assert(num_avail_rooms > 0); // must have added at least one room
		float stairs_scale(1.0);

		for (unsigned N = 0; N < 4; ++N) {
			unsigned const rand_ix(rgen.rand()); // choose a random starting room to make into a stairwell or add an elevator to

			// try all available rooms starting with the selected one to see if we can fit a stairwell/elevator in any of them
			for (unsigned n = 0; n < num_avail_rooms; ++n) {
				unsigned const stairs_room(rooms_start + (rand_ix + n)%num_avail_rooms);
				room_t &room(interior->rooms[stairs_room]);
				assert(room.part_id == part_ix); // sanity check
				if (room.is_nested()) continue; // nested rooms typically shouldn't have stairs or elevators

				if (is_parking_str) {
					// always add an elevator to parking structures; place against a wall;
					// it should be by the main exterior door, but that hasn't been placed yet, and it may be far from the building center
					bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
					float const dsign(dir ? -1.0 : 1.0), wall_pos(room.d[dim][dir]), end_spacing(0.35*room.get_sz_dim(!dim)); // center 50% of wall to avoid ramps
					elevator_t elevator(room, stairs_room, dim, !dir, 1, 1); // at_edge=1, interior=1 (needed for terrain culling)
					elevator.d[dim][!dir] = wall_pos + dsign*ewidth;
					float const center_pos(rgen.rand_uniform((room.d[!dim][0] + end_spacing), (room.d[!dim][1] - end_spacing)));
					set_wall_width(elevator, center_pos, 0.5*ewidth, !dim);
					// can't be blocked by anything?
					add_or_extend_elevator(elevator, 1); // add=1
					elevator_cut = elevator;
					// maybe add U-shaped stairs to the side of the elevator
					int const stairs_mode(rgen.rand_bool() ? 1 : 2); // 0=straight side stairs, 1=U-shaped stairs parallel to elevator, 2=U shaped stairs perpendicular

					if (stairs_mode > 0) {
						bool const perp(stairs_mode == 2), sdir(center_pos < room.get_center_dim(!dim)); // closer to the center
						float const elevator_side(elevator.d[!dim][sdir] + (sdir ? 1.0 : -1.0)*0.25*wall_thickness); // account for ecap width
						float const width(2.0*doorway_width), depth(2.2*doorway_width);
						cube_t stairs(room);
						stairs.d[!dim][!sdir] = elevator_side; // abuts elevator
						stairs.d[!dim][ sdir] = elevator_side + (sdir ? 1.0 : -1.0)*(perp ? depth : width);
						stairs.d[ dim][  dir] = wall_pos + dsign*0.3*wall_thickness; // shift slightly away from wall so that exterior of stairs wall matches elevator
						stairs.d[ dim][! dir] = wall_pos + dsign*(perp ? width : depth);
						assert(stairs.is_strictly_normalized());
						
						if (room.contains_cube_xy(stairs)) { // can this ever fail?
							sshape     = SHAPE_U;
							stairs_cut = stairs;
							stairs_dim = (dim ^ perp);
							force_stairs_dir = (perp ? !sdir : dir);
							room.has_stairs  = 255; // stairs on all floors
							break; // done - no need to add stairs below
						}
					}
				}
				if (add_elevator) {
					if (min(room.dx(), room.dy()) < 2.0*ewidth) continue; // room is too small to place an elevator
					bool const no_ext_walls(!(is_prison() && interior->elevators.empty()));
					cube_t clip_bounds(part);
					clip_bounds.expand_by_xy(-get_trim_thickness());
					bool placed(0);

					for (unsigned y = 0; y < 2 && !placed; ++y) { // try all 4 corners
						for (unsigned x = 0; x < 2 && !placed; ++x) {
							int const wtype_x(classify_room_wall(room, room.z1(), 0, x, 1)), wtype_y(classify_room_wall(room, room.z1(), 1, y, 1)); // include partial sep walls
							// don't place elevators between parts where they could block doorways, except for prisons, which have fewer placement options
							if (no_ext_walls && (wtype_x == ROOM_WALL_SEP || wtype_y == ROOM_WALL_SEP)) continue;
							if (!is_cube  () && (wtype_x == ROOM_WALL_EXT || wtype_y == ROOM_WALL_EXT)) continue; // no exterior walls of non-cube buildings
							float const xval(room.d[0][x] + (x ? -ewidth : ewidth)), yval(room.d[1][y] + (y ? -ewidth : ewidth)), shrink(0.01*ewidth);
							// check room interior edge for intersection with windows if not a prison
							if (no_ext_walls) {
								if (wtype_x == ROOM_WALL_EXT && is_val_inside_window(part, 0, xval, window_hspacing[0], window_border)) continue;
								if (wtype_y == ROOM_WALL_EXT && is_val_inside_window(part, 1, yval, window_hspacing[1], window_border)) continue;
							}
							bool const dim(rgen.rand_bool()), at_edge(wtype_x == ROOM_WALL_EXT || wtype_y == ROOM_WALL_EXT);
							elevator_t elevator(room, stairs_room, dim, !(dim ? y : x), at_edge, room.interior); // elevator shaft
							elevator.d[0][!x] = xval;
							elevator.d[1][!y] = yval;
							// shrink to leave a small gap between the outer wall to prevent z-fighting
							if (wtype_x == ROOM_WALL_EXT) {elevator.d[0][x] += (x ? -shrink : shrink);}
							if (wtype_y == ROOM_WALL_EXT) {elevator.d[1][y] += (y ? -shrink : shrink);}
							elevator.intersect_with_cube_xy(clip_bounds); // last attempt to clip to shrunk part in case elevator was on a sep wall corner (for prisons)
							if (has_bcube_int(elevator, interior->exclusion)) continue; // try again
							if (is_cube_close_to_doorway(elevator, room))     continue; // try again
							if (check_skylight_intersection(elevator))        continue; // check skylights; is this necessary?
							if (!check_cube_within_part_sides(elevator))      continue; // outside building; do we need to check for clearance?
							if (stairs_or_elevator_blocked_by_nested_room(elevator, stairs_room)) continue;
							add_or_extend_elevator(elevator, 1); // add=1
							elevator_cut = elevator;
							placed       = 1; // successfully placed
						} // for x
					} // for y
					if (!placed) continue; // try another room
					room.has_elevator = 1;
					add_elevator      = 0;
					//room.set_no_geom();
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
					}
					// see if we can push the stairs to the wall on one of the sides without blocking a doorway
					bool const first_dir(rgen.rand_bool());

					for (unsigned d = 0; d < 2; ++d) {
						bool const dim(!stairs_dim), dir(bool(d) ^ first_dir);
						int const wall_type(classify_room_wall(room, room.z1(), dim, dir, 1));
						// if the room is on the edge of the part that's not on the building exterior, then this room connects two parts and we need to place a door here later
						if (              wall_type == ROOM_WALL_SEP) continue; // include partial sep walls
						if (!is_cube() && wall_type == ROOM_WALL_EXT) continue; // no exterior walls of non-cube buildings
						cube_t cand(cutout);
						// add small gap to prevent z-fighting and FP accuracy asserts
						float const shift((cand.d[dim][dir] - room.d[dim][dir]) - (dir ? -1.0 : 1.0)*wall_thickness); // negative if dir==1
						cand.translate_dim(dim, -shift); // close the gap - flush with the wall
						if (is_cube_close_to_doorway(cand, room)) continue;
						if (stairs_or_elevator_blocked_by_nested_room(cand, stairs_room))  continue;
						if (!elevator_cut.is_all_zeros() && cand.intersects(elevator_cut)) continue; // too close to parking structure elevator

						if (!is_cube()) { // check for stairs outside or too close to building walls
							cube_t stairs_ext(cand);
							stairs_ext.expand_in_dim(stairs_dim, doorway_width);
							if (!check_cube_within_part_sides(stairs_ext)) continue; // outside building
						}
						cutout = cand;
						stairs_against_wall[dir] = 1;
						break; // keep if it's good
					} // for d
					bool const against_wall(stairs_against_wall[0] || stairs_against_wall[1]);
					float const room_width(room.get_sz_dim(!stairs_dim));
					// skip if we can't push against a wall and the room is too narrow for space around the stairs to allow doors to open and people to walk
					if (!against_wall && room_width < (1.1*2.0*doorway_width + stairs_sz)) continue;

					if (!against_wall && stairs_or_elevator_blocked_by_nested_room(cutout, stairs_room)) { // check if left in room center
						// if this is a prison, consider areas next to interior walls where there are no cells
						bool wall_dir(0);
						if (!is_prison() || !place_stairs_in_prison_room(cutout, stairs_room, stairs_dim, wall_dir, rgen)) continue;
						stairs_against_wall[wall_dir] = 1;
					}
					if (is_house && against_wall) { // try to convert to L-shaped stairs
						float const stairs_len(cutout.get_sz_dim(stairs_dim)), stairs_width(cutout.get_sz_dim(!stairs_dim)); // primary/lower segment

						if (stairs_len > 2.0*stairs_width) { // not too short and wide
							bool const bdir(stairs_against_wall[0]); // bend away from the wall
							float const room_edge_min_space(max(2.0f*stairs_sz, 0.5f*room_width)), dscale(bdir ? 1.0 : -1.0);
							float const max_expand(min(2.0f*doorway_width, (stairs_len - stairs_width))); // primary/lower segment is always longer
							float expand_l(max_expand);

							for (unsigned M = 0; M < 5; ++M) { // 5 expansion attempts (50% reduction max)
								cube_t cutout_l(cutout);
								cutout_l.d[!stairs_dim][bdir] += dscale*expand_l; // extend out for the L
								cube_t cutout_l_pad(cutout_l), cutout_with_wall(cutout_l);
								// expand on ends to include space for landing support wall to not intersect doors; expand in both dirs since dir hasn't been determined yet;
								// only needed for stairs spanning more than 2 floors since the bottom landing support wall is under the landing rather than beside it
								if (num_floors > 2) {cutout_with_wall.expand_in_dim(stairs_dim, wall_thickness);}
								cutout_l_pad.d[!stairs_dim][bdir] += dscale*room_edge_min_space; // extend out for the stairs exit

								if (room.contains_cube_xy(cutout_l_pad) && !is_cube_close_to_doorway(cutout_with_wall, room)) {
									sshape   = SHAPE_L;
									cutout   = cutout_l;
									bend_dir = bdir;
									break; // done
								}
								expand_l -= 0.1*max_expand; // shorten it
							} // for M
						}
					}
					if (!is_cube()) { // check for stairs outside or too close to building walls
						cube_t stairs_ext(cutout);
						stairs_ext.expand_in_dim(stairs_dim, doorway_width);
						
						if (!check_cube_within_part_sides(stairs_ext)) { // outside building
							// try to move toward building center/away from exterior
							bool const move_dir(cutout.get_center_dim(!stairs_dim) < bcube.get_center_dim(!stairs_dim));
							float const move_amt((move_dir ? 1.0 : -1.0)*0.2*room.get_sz_dim(!stairs_dim));
							cutout    .translate_dim(!stairs_dim, move_amt);
							stairs_ext.translate_dim(!stairs_dim, move_amt);
							if (is_cube_close_to_doorway(cutout, room) || !check_cube_within_part_sides(stairs_ext)) continue; // still bad, skip this room
						}
					}
					if (!is_house && against_wall) {sshape = SHAPE_WALLED_SIDES;} // add wall between room and office stairs if against a room wall
					if (interior->landings.empty()) {interior->landings.reserve(num_floors-1);}
					assert(cutout.is_strictly_normalized());
					stairs_cut      = cutout;
					room.has_stairs = 255; // stairs on all floors
					if (!against_wall) {room.set_has_center_stairs();}
					if (use_hallway || !pri_hall.is_all_zeros()) {room.set_no_geom();} // no geom in an office with stairs for buildings with hallways
				}
				break; // success - done
			} // for n
			if (!must_add_stairs || add_elevator || !stairs_cut.is_all_zeros()) break; // successfully placed stairs, or not required to place stairs
			stairs_scale -= 0.1; // shrink stairs a bit and try again
		} // for N
		//if (is_house && interior->stairwells.empty()) {interior->is_unconnected = 1;} // fails too often/too slow
	} // end stairs/elevator placement

	// add ceilings and floors; we have num_floors+1 separators; the first is only a floor, and the last is only a ceiling
	cube_t C(part);
	C.z1() = z; C.z2() = z + fc_thick;
	unsigned const floors_start(interior->floors.size());
	interior->floors.push_back(C); // ground floor, full area
	bool const has_stairs(!stairs_cut.is_all_zeros()), has_elevator(!elevator_cut.is_all_zeros());
	bool const stairs_dir((force_stairs_dir < 2) ? force_stairs_dir : (has_stairs ? rgen.rand_bool() : 0)); // same for every floor, could maybe alternate for stairwells
	float const floor_vert_spacing(is_single_floor ? part.dz() : (is_retail_part(part) ? retail_floor_levels : 1)*window_vspacing);
	cube_t &first_cut(has_elevator ? elevator_cut : stairs_cut); // elevator is larger
	//unsigned const first_ceiling_ix(max(1U, unsigned(retail_floor_levels))); // skip first floor and any ground floor retail space
	unsigned last_landing_ix(0);
	z += floor_vert_spacing; // move to next floor

	for (unsigned f = 1; f < num_floors; ++f, z += floor_vert_spacing) { // skip first floor; draw pairs of floors and ceilings
		cube_t to_add[8]; // up to 2 cuts for stairs + elevator
		float const zc(z - fc_thick), zf(z + fc_thick);

		if (!has_stairs && !has_elevator) {to_add[0] = part;} // neither - add single cube
		else {
			bool const is_at_top(f+1 == num_floors && !extended_from_above);
			assert(part.contains_cube_xy(first_cut));
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
 				landing_t landing(stairs_cut, 0, f, stairs_dim, stairs_dir, stairs_have_railing,
					((f == 1 && sshape == SHAPE_WALLED_SIDES) ? (stairs_shape)SHAPE_WALLED : sshape), 0, is_at_top);
				set_cube_zvals(landing, zc, zf);
				landing.set_against_wall(stairs_against_wall);
				landing.bend_dir = bend_dir;
				last_landing_ix  = interior->landings.size();
				interior->landings.push_back(landing);

				if (f == 1) { // only add for first floor
					interior->stairwells.emplace_back(stairs_cut, num_floors, stairs_dim, stairs_dir, sshape);
					interior->stairwells.back().set_against_wall(stairs_against_wall);
					interior->stairwells.back().bend_dir = bend_dir;
				}
			}
			if (has_elevator) {
				assert(!interior->elevators.empty());
				landing_t landing(elevator_cut, 1, f, interior->elevators.back().dim, interior->elevators.back().dir);
				set_cube_zvals(landing, zc, zf);
				interior->landings.push_back(landing);
			}
		}
		for (unsigned i = 0; i < 8; ++i) { // skip zero area cubes from stairs/elevator shafts along an exterior wall
			cube_t &cf(to_add[i]);
			if (!cf.is_zero_area()) {interior->add_ceil_floor_pair(cf, zc, z, zf);}
		}
	} // for f
	bool has_roof_access(0);

	if (must_add_stairs && has_stairs && !is_house && roof_type == ROOF_TYPE_FLAT) { // add roof access for stairs
		bool const is_sloped(sshape != SHAPE_U);
		cube_t box(stairs_cut);
		if (!is_sloped) {box.expand_by_xy(fc_thick);}
		box.z1()  = z + floor_thickness;
		box.z2()  = z + window_vspacing - (is_sloped ? 0.15 : 0.2)*window_vspacing; // slightly lower than a normal floor
		cube_t check_box(box);
		bool const opening_dir(stairs_dir ^ (sshape == SHAPE_U)); // U-shaped stairs have exit on same side as entrance
		check_box.d[stairs_dim][opening_dir] += (opening_dir ? 1.0 : -1.0)*doorway_width; // expand at stairs exit to ensure clearance

		// check for overlap with other parts or skylights (should we check in front?)
		if (!has_bcube_int_no_adj(check_box, parts) && !check_skylight_intersection(check_box) && (!has_helipad || !get_helipad_bcube().intersects_xy(check_box))) {
			float const zc(z - fc_thick);
			cube_t to_add[4]; // only one cut / 4 cubes (-y, +y, -x, +x)
			subtract_cube_xy(part, stairs_cut, to_add);
			interior->landings[last_landing_ix].is_at_top = 0; // previous landing is no longer at the top
			landing_t landing(stairs_cut, 0, num_floors, stairs_dim, stairs_dir, 1, sshape, 1, 1); // stairs_have_railing=1, roof_access=1, is_at_top=1
			set_cube_zvals(landing, zc, z); // no floor above
			landing.set_against_wall(stairs_against_wall);
			interior->landings.push_back(landing);
			interior->stairwells.back().z2() += fc_thick; // extend upward
			interior->stairwells.back().z1() += fc_thick; // required to trick roof clipping into treating this as a stack connector stairwell

			for (unsigned i = 0; i < 4; ++i) { // skip zero area cubes from stairs/elevator shafts along an exterior wall
				cube_t &c(to_add[i]);
				if (c.is_zero_area()) continue;
				set_cube_zvals(c, zc, z);
				add_ceiling_cube_no_skylights(c);
				c.set_to_zeros();
			}
			bool const dir(stairs_dir ^ (sshape == SHAPE_U)), u_side(dir); // see logic in building_t::add_stairs_and_elevators() for side
			box.z1() = z;
			// add a door to the roof
			float const door_shift(0.2*(dir ? -1.0 : 1.0)*fc_thick);
			cube_t door(box);
			door.d[stairs_dim][ dir] += door_shift; // shift slightly to fill the gap
			door.d[stairs_dim][!dir]  = door.d[stairs_dim][dir];
			if (!is_sloped) {door.d[!stairs_dim][u_side] = door.get_center_dim(!stairs_dim);} // U-shaped stairs; only the open half
			add_door(door, part_ix, stairs_dim, dir, 1, 1); // roof_access=1
			// clear any roof objects that are in the way
			cube_t clear_cube(box);
			clear_cube.d[stairs_dim][dir] += (dir ? 1.0 : -1.0)*window_vspacing; // clear out space in front of the door
			clear_cube.expand_in_dim(!stairs_dim, doorway_width); // add clearance to the sides to make sure the player can reach other areas of the roof
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
				unsigned const tquad_type(tquad_with_ix_t::TYPE_ROOF_ACC);
				roof_tquads.emplace_back(top, tquad_type);
				tquad_t frame_top(4); // top of door frame

				for (unsigned s = 0; s < 2; ++s) {
					tquad_t side(3); // triangle
					side.pts[0] = pts[bot_vix[s]]; // corner
					side.pts[1] = pts[top_vix[s]]; // bottom
					side.pts[2] = top.pts[top_vix[s]]; // top
					if (s) {swap(side.pts[0], side.pts[1]);} // make normal point outward for correct BFC
					roof_tquads.emplace_back(side, tquad_type);

					tquad_t frame(4); // quads on sides of door frame
					frame.pts[0] = frame.pts[3] =     pts[top_vix[s]]; // bottom
					frame.pts[1] = frame.pts[2] = top.pts[top_vix[s]] - vector3d(0.0, 0.0, 0.05*fc_thick); // top, shift slightly down to match the slope
					float const frame_width(0.08*((bool(s)^dir^stairs_dim^1) ? -1.0 : 1.0)*box.get_sz_dim(!stairs_dim));
					for (unsigned d = 0; d < 2; ++d) {frame.pts[(s<<1)+d][!stairs_dim] +=    frame_width;} // move toward center of door
					for (unsigned n = 0; n < 4; ++n) {frame.pts[n]       [ stairs_dim] += 0.6*door_shift;} // shift back behind the door
					roof_tquads.emplace_back(frame, tquad_type);
					frame_top.pts[2*s +   s] = frame.pts[s+1];
					frame_top.pts[2*s + 1-s] = frame.pts[s+1] - vector3d(0.0, 0.0, 0.5*fc_thick);
				} // for s
				roof_tquads.emplace_back(frame_top, tquad_type);
			}
			else { // U-shaped stairs; box roof
				cube_t hole(stairs_cut), front(box);
				hole.expand_by_xy(0.1*fc_thick); // to prevent z-fighting
				front.d[ stairs_dim][!dir   ] = hole.d[stairs_dim][dir];
				hole .d[ stairs_dim][ dir   ] = box. d[stairs_dim][dir]; // move edge flush with box to remove this wall and create an opening
				front.d[!stairs_dim][!u_side] = front.get_center_dim(!stairs_dim); // block off the non-opening half
				subtract_cube_xy(box, hole, to_add);

				for (unsigned i = 0; i < 4; ++i) {
					cube_t &c(to_add[i]);
					if (!c.is_zero_area()) {details.emplace_back(c, (uint8_t)ROOF_OBJ_SCAP);} // skip open side
				}
				box.z1() = front.z2() = box.z2() - fc_thick;
				details.emplace_back(box,   (uint8_t)ROOF_OBJ_SCAP); // top
				details.emplace_back(front, (uint8_t)ROOF_OBJ_SCAP); // front half
				max_eq(bcube.z2(), box.z2());
			}
			interior->stairwells.back().roof_access = 1;
			has_roof_access = 1;
		}
	}
	if (!has_roof_access) { // roof ceiling, full area
		set_cube_zvals(C, (z - fc_thick), z);
		
		if (is_house && part_ix == get_attic_part_ix() && add_attic_access_door(C, part_ix, num_floors, rooms_start, rgen)) { // primary/upper part only
			cube_t ceiling_parts[4];
			subtract_cube_xy(C, interior->attic_access, ceiling_parts);
			float const fc_mid_z(C.zc()); // split between the ceiling and floor parts

			for (unsigned i = 0; i < 4; ++i) {
				cube_t &c(ceiling_parts[i]);
				if (c.is_zero_area()) continue;
				c.z2() = fc_mid_z;
				interior->ceilings.push_back(c);
				c.z1() = fc_mid_z;
				c.z2() = C.z2();
				interior->floors.push_back(c);
			} // for i
		}
		else {
			add_ceiling_cube_no_skylights(C);
		}
	}
	std::reverse(interior->floors.begin()+floors_start, interior->floors.end()); // order floors top to bottom to reduce overdraw when viewed from above
}

bool building_t::can_extend_stairs_to_pg(unsigned &stairs_ix) const {
	if (!has_parking_garage || (!has_pri_hall() && !is_parking())) return 0; // no parking garage or primary hall, except for parking structure (but it doesn't help?)
	float stairs_zmax(ground_floor_z1 + get_floor_thickness());

	// check ground floor stairs, then possibly stairs above the retail floor; prefer stairs on the ground floor when there are both
	for (unsigned pass = 0; pass < 2; ++pass) {
		if (pass == 1) {
			if (!has_retail()) break; // only one pass
			stairs_zmax += retail_floor_levels*get_window_vspace(); // assume stairs can be extended down to retail ground floor
		}
		for (unsigned i = 0; i < interior->stairwells.size(); ++i) {
			stairwell_t const &s(interior->stairwells[i]);
			if (s.z1() < ground_floor_z1 || s.z1() > stairs_zmax) continue; // not ground floor stairs (or just above ground floor if retail)
			if (has_pri_hall() && !pri_hall.contains_cube_xy(s))  continue; // not in primary hall or the retail room below
			stairs_ix = i; // can be primary hall stairs or stairs extended into retail area below
			return 1;
		}
	} // for pass
	return 0;
}

void building_t::add_ceiling_cube_no_skylights(cube_t const &c) {
	assert(interior);

	if (!check_skylight_intersection(c)) {
		interior->ceilings.push_back(c);
		return;
	}
	vect_cube_t cubes;
	cubes.push_back(c);
	for (cube_t const &skylight : skylights) {subtract_cube_from_cubes(skylight, cubes);}
	for (cube_t const &cube : cubes) {interior->ceilings.push_back(cube);}
}

bool building_t::check_cube_intersect_walls(cube_t const &c) const {
	assert(interior);
	return (has_bcube_int(c, interior->walls[0]) || has_bcube_int(c, interior->walls[1]));
}

// Note: when querying for elevators and stairs, c should be padded to include the entrance/exit; c_nopad is the unpadded cube
bool building_t::is_valid_stairs_elevator_placement(cube_t const &c, cube_t const &c_nopad, float pad, int dim, bool check_walls, bool check_private_rooms) const {
	assert(interior);
	// check if any previously placed walls intersect this cand stairs/elevator; we really only need to check the walls from <part> and *p though
	if (interior->is_blocked_by_stairs_or_elevator(c_nopad, pad)) return 0;
	if (check_walls && check_cube_intersect_walls(c))             return 0; // Note: extb_walls_start may not yet be set because extended basement has not been added
	if (is_cube_close_to_doorway(c, cube_t(), pad, 1))            return 0; // check for open doors to avoid having the stairs intersect an open door
	if (check_skylight_intersection(c))                           return 0;

	if (dim <= 1 && pad > 0.0) { // dim was specified (not AI placement); check if containing rooms have space on either side for the player to walk
		for (room_t const &r : interior->rooms) {
			if (r.is_ext_basement()) break; // no need to check extended basement rooms
			if (!r.contains_cube_xy_overlaps_z(c)) continue; // stairs/elevator not contained in this room
			if (max((c.d[!dim][0] - r.d[!dim][0]), (r.d[!dim][1] - c.d[!dim][1])) < pad) return 0;
		}
	}
	if (!is_house && has_pri_hall() && c.z2() > ground_floor_z1) { // office building with primary hallway; not in basement
		if (c.z1() < ground_floor_z1 + get_window_vspace()) { // overlaps the ground floor
			// add extra padding around exterior doors to avoid blocking them with stairs/elevator
			float const door_width(DOOR_WIDTH_SCALE_OFFICE*get_office_bldg_door_height());
			point end_pt;
			end_pt[!hallway_dim] = pri_hall.get_center_dim(!hallway_dim); // assumes door is centered in the hallway

			for (unsigned d = 0; d < 2; ++d) { // check both hallway ends
				end_pt[hallway_dim] = pri_hall.d[hallway_dim][d];
				cube_t blocked(end_pt);
				blocked.expand_by_xy(door_width); // ensure at least one door width of spacing around the door
				blocked.intersect_with_cube_xy(pri_hall);
				if (blocked.intersects_xy(c)) return 0;
			}
		}
	}
	if (check_private_rooms && is_apt_or_hotel()) { // Note: may not yet be assigned
		for (room_t const &r : interior->rooms) {
			if (r.is_ext_basement()) break; // no need to check extended basement rooms
			if (r.is_apt_or_hotel_room() && r.intersects(c)) return 0;
		}
	}
	if (is_industrial() && c.z1() > ground_floor_z1) { // check entrance door
		bool const dim(interior->ind_info->entrance_dim), dir(interior->ind_info->entrance_dir);
		float const doorway_width(get_doorway_width());
		cube_t door;
		set_cube_zvals(door, ground_floor_z1+get_fc_thickness(), ground_floor_z1+get_floor_ceil_gap());
		set_wall_width(door, interior->ind_info->entrance_pos,  doorway_width, !dim); // 2x total width for double office door
		set_wall_width(door, get_industrial_area().d[dim][dir], doorway_width,  dim);
		if (c.intersects(door)) return 0;
	}
	return 1;
}

template<typename T> void subtract_cube_from_cube(T const &c, cube_t const &s, vector<T> &out, bool clear_out) { // XY only
	if (clear_out) {out.clear();}
	if (c.y1() < s.y1()) {T C(c); C.y2() = s.y1(); out.push_back(C);} // bottom
	if (c.y2() > s.y2()) {T C(c); C.y1() = s.y2(); out.push_back(C);} // top
	if (c.x1() < s.x1()) {T C(c); max_eq(C.y1(), s.y1()); min_eq(C.y2(), s.y2()); C.x2() = s.x1(); out.push_back(C);} // left center
	if (c.x2() > s.x2()) {T C(c); max_eq(C.y1(), s.y1()); min_eq(C.y2(), s.y2()); C.x1() = s.x2(); out.push_back(C);} // right center
}
template<typename T> void subtract_cube_from_cube_inplace(cube_t const &s, vector<T> &cubes, unsigned &ix, unsigned &iter_end) { // Note: ix is an index to cubes
	unsigned const prev_sz(cubes.size());
	assert(ix < prev_sz);
	T const c(cubes[ix]); // deep copy - reference will become invalid
	subtract_cube_from_cube(c, s, cubes);
	bool const none_added(cubes.size() == prev_sz);
	cubes[ix] = cubes.back(); cubes.pop_back(); // reuse this slot for one of the output cubes (or move the last cube here if there are no output cubes)
	if (none_added) {--ix; --iter_end;} // no cubes added, last cube was swapped into this slot and needs to be reprocessed
}
// subtracts in X and Y only; zval_mode: 0=check floor/building ext walls, 1=ignore zvals, 2=check zval overlap
template<typename T> void subtract_cubes_from_cube(cube_t const &c, vector<T> const &sub, vect_cube_t &out, vect_cube_t &out2, int zval_mode) {
	out.clear();
	out.push_back(c);

	for (auto s = sub.begin(); s != sub.end(); ++s) {
		if (zval_mode == 0 && (s->z1() <= c.z1() || s->z1() >= c.z2() || s->z2() <= c.z2())) continue; // not correct floor
		if (zval_mode == 2 && (s->z1() >= c.z2() || s->z2() <= c.z1())) continue; // no zval overlap
		if (!c.intersects_xy_no_adj(c)) continue; // no overlap with orig cube (optimization)
		out2.clear();

		// clip all of out against *s, write results to out2, then swap with out
		for (auto i = out.begin(); i != out.end(); ++i) {
			if (!i->intersects_xy_no_adj(*s)) {out2.push_back(*i); continue;} // no overlap, keep entire cube
			subtract_cube_from_cube(*i, *s, out2);
		}
		out.swap(out2);
	} // for s
}
template void subtract_cubes_from_cube(cube_t const &c, vector<cube_t>         const &sub, vect_cube_t &out, vect_cube_t &out2, int zval_mode); // explicit instantiation
template void subtract_cubes_from_cube(cube_t const &c, vector<stairs_place_t> const &sub, vect_cube_t &out, vect_cube_t &out2, int zval_mode); // explicit instantiation
template void subtract_cubes_from_cube(cube_t const &c, vector<elevator_t>     const &sub, vect_cube_t &out, vect_cube_t &out2, int zval_mode); // explicit instantiation

template<typename T> bool subtract_cube_from_cubes(cube_t const &s, vector<T> &cubes, vect_cube_t *holes, bool clip_in_z, bool include_adj, bool no_z_test) {
	unsigned iter_end(cubes.size()); // capture size before splitting
	bool was_clipped(0);
	vect_cube_t top_bot_cubes;

	for (unsigned i = 0; i < iter_end; ++i) {
		T const &c(cubes[i]);

		if (no_z_test) {
			if (!(include_adj ? c.intersects_xy(s) : c.intersects_xy_no_adj(s))) continue; // keep it
		}
		else {
			if (!(include_adj ? c.intersects(s) : c.intersects_no_adj(s))) continue; // keep it
		}
		if (holes) {
			cube_t hole(c); // always a cube
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
			T shared(c);
			shared.intersect_with_cube_xy(s);
			if (shared.z1() < s.z1()) {T bot(shared); bot.z2() = s.z1(); top_bot_cubes.push_back(bot);} // bottom part
			if (shared.z2() > s.z2()) {T top(shared); top.z1() = s.z2(); top_bot_cubes.push_back(top);} // top part
		}
		subtract_cube_from_cube_inplace(s, cubes, i, iter_end); // Note: invalidates c reference
		was_clipped = 1;
	} // for i
	vector_add_to(top_bot_cubes, cubes);
	return was_clipped;
}
template bool subtract_cube_from_cubes(cube_t const &s, vector<cube_t>         &cubes, vect_cube_t *holes, bool clip_in_z, bool include_adj, bool no_z_test);
template bool subtract_cube_from_cubes(cube_t const &s, vector<cube_with_ix_t> &cubes, vect_cube_t *holes, bool clip_in_z, bool include_adj, bool no_z_test);

void subtract_cube_from_floor_ceil(cube_t const &c, vect_cube_t &fs) {
	unsigned iter_end(fs.size()); // capture orig size

	for (unsigned i = 0; i < iter_end; ++i) {
		cube_t const &cur(fs[i]);
		if (cur.z1() > c.z2() || cur.z2() < c.z1()) continue; // no z overlap
		if (cur.intersects_no_adj(c)) {subtract_cube_from_cube_inplace(c, fs, i, iter_end);} // Note: invalidates cur reference
	}
}

void get_room_cands(vector<room_t> const &rooms, unsigned part_ix, cube_t const &place_region, vector2d const &min_sz, bool check_private_rooms, vect_cube_with_ix_t &room_cands) {
	room_cands.clear();

	for (unsigned room_id = 0; room_id < rooms.size(); ++room_id) {
		room_t const &r(rooms[room_id]);
		if (r.part_id != part_ix || !r.intersects(place_region) || r.is_bathroom_rtype() || r.is_jail_cell()) continue; // Note: r.is_nested() allowed for some room types
		if (check_private_rooms && r.is_apt_or_hotel_room()) continue;
		cube_t cand(r);
		cand.intersect_with_cube_xy(place_region);
		vector2d const room_sz(cand.get_size_xy());
		if (room_sz[0] < min_sz[0] || room_sz[1] < min_sz[1]) continue; // room too small
		room_cands.emplace_back(cand, room_id);
	} // for r
}

void building_t::connect_stacked_parts_with_stairs(rand_gen_t &rgen, cube_t const &part, unsigned lower_part_ix) { // and extend elevators vertically; part is lower

	//highres_timer_t timer("Connect Stairs", 1, 1, 1); // track_not_print=1; 256ms total
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness), wall_thickness(get_wall_thickness());
	float const doorway_width(get_nominal_doorway_width()), stairs_len(4.0*doorway_width), railing_pad(0.5*wall_thickness);
	bool const is_basement(has_basement() && part == get_basement()), use_basement_stairs(is_basement && is_house); // office basement has regular stairs
	bool const is_retail(is_retail_part(part));
	bool const check_private_rooms = 0; // this could go either way; which is worse - an unconnected stacked part, or public stairs in a hotel room or apartment?
	// use fewer iterations on tiled buildings to reduce the frame spikes when new tiles are generated
	unsigned const iter_mult_factor(global_building_params.gen_inf_buildings() ? 5 : 10), num_iters(20*iter_mult_factor);
	unsigned const num_floors(is_retail ? 1 : calc_num_floors(part, window_vspacing, floor_thickness)); // retail area is always one floor
	assert(num_floors > 0);

	if (part.z2() < bcube.z2()) { // if this is the top floor, there is nothing above it (but roof geom may get us into this case anyway)
		vect_door_stack_t doorways;
		vect_cube_with_ix_t room_cands[2]; // {upper, lower}
		bool connected(0);

		for (unsigned upper_part_ix = 0; upper_part_ix < real_num_parts; ++upper_part_ix) { // find the part on the top
			cube_t const &p(parts[upper_part_ix]);
			if (upper_part_ix == lower_part_ix) continue; // skip self
			if (p.z1() != part.z2())            continue; // *p not on top of part
			if (!part.intersects_xy(p))         continue; // no XY overlap
			cube_t shared(part);
			shared.intersect_with_cube(p); // dz() == 0
			cube_t pref_shared(shared);
			// Note: office building parts are sorted top to bottom, so any part above <part> should be before it in parts - but we don't want to rely on that here;
			// however, this does mean that the part above this one has already been processed; except for stacked houses, which are ordered {bottom, top}
			float stairs_width(1.2*doorway_width); // relatively small
			float stairs_pad(doorway_width), len_with_pad(stairs_len + 2.0*stairs_pad); // pad both ends of stairs to make sure player has space to enter/exit
			if (max(shared.dx(), shared.dy()) < 1.0*len_with_pad || min(shared.dx(), shared.dy()) < 1.2*stairs_width) continue; // too small to add stairs between these parts

			if (has_pri_hall() && part.contains_cube(pri_hall) && pri_hall.intersects_xy(shared)) { // have a primary hallway in this part
				pref_shared.intersect_with_cube(pri_hall);
				if (max(pref_shared.dx(), pref_shared.dy()) < 1.2*len_with_pad || min(pref_shared.dx(), pref_shared.dy()) < 1.5*stairs_width) {pref_shared = shared;} // too small
			}
			// place stairs in shared area if there's space and no walls are in the way for either the room above or below
			cube_t cand;
			stairs_shape sshape(SHAPE_STRAIGHT);
			bool cand_is_valid(0), dim(0), stairs_dir(0), bend_dir(0), add_railing(1), stack_conn(1), is_at_top(0);
			float const stairs_height((is_retail ? retail_floor_levels : 1)*window_vspacing);
			float const cand_z2(part.z2() + fc_thick); // top of bottom floor of upper part *p
			float const cand_z1(cand_z2 - stairs_height); // top of top floor for this/lower part
			unsigned stairs_ix(0);
			
			// try to extend primary hallway or parking structure stairs down to parking garage below; should this apply to all ground floor stairwells?
			// single floor parking garage, or upper floor connect only
			if (is_basement && (is_parking() || part.contains_cube_xy(pri_hall)) && can_extend_stairs_to_pg(stairs_ix)) {
				stairwell_t &s(interior->stairwells[stairs_ix]);
				s.extends_below       = 1;
				pri_hall_stairs_to_pg = 1;
				cand = s; dim = s.dim; stairs_dir = s.dir; sshape = s.shape; // copy fields from these stairs and extend down
				stack_conn    = 0; // not stacked - extended main stairs
				cand_is_valid = 1;
				// Note: we likey extended the stairs down to the lower parking garage level in add_ceilings_floors_stairs() for the basement
				cube_t stairs_bot(s);
				stairs_bot.z2() = s.z1() + window_vspacing; // limit to bottom landing
				find_and_merge_with_landing(interior->landings, stairs_bot, sshape, 1, 0); // merge with bottom landing; num_floors=1; is_above=0
			}
			else if (!(is_basement && (is_house || is_parking())) && is_cube()) { // not for house basements, parking garages, or non-cube buildings
				// try to extend an existing stairwell on the part above or below upward/downward
				for (unsigned ab = 0; ab < 2 && !cand_is_valid; ++ab) { // extend {below, above}
					cube_t const &targ_part(ab ? part : p);

					for (stairwell_t &s : interior->stairwells) {
						if (s.in_ext_basement) continue; // not stacked/stackable, skip (also should skip basement stairs themselves?)
						if (!shared.contains_cube_xy(s) || s.z2() < targ_part.zc() || s.z1() > targ_part.zc()) continue; // stairs not contained in both and crossing one part
						// check for clearance on the other part
						cube_t s_bc(s);
						set_cube_zvals(s_bc, cand_z1, cand_z2);
						if (ab == 0) {s_bc.z1() += floor_thickness; s_bc.z2() -= floor_thickness;} // shrink to lower part
						else         {s_bc.z1() += window_vspacing + floor_thickness; s_bc.z2() += window_vspacing - floor_thickness;} // move to upper part
						cube_t ext_cube(s_bc);
						if (bool(ab) ^ s.dir) {ext_cube.d[s.dim][0] -= stairs_pad;} // add padding on exit side
						else                  {ext_cube.d[s.dim][1] += stairs_pad;} // add padding on exit side
						if (!shared.contains_cube_xy(ext_cube))           continue; // test for space to enter and exit
						if (has_bcube_int(ext_cube, interior->exclusion)) continue; // bad placement
						bool const allow_clip_walls(0); // clipping walls rarely helps and tends to create some strange stairs
						if (!is_valid_stairs_elevator_placement(ext_cube, s_bc, stairs_pad, s.dim, !allow_clip_walls, check_private_rooms)) continue; // bad placement
						(ab ? s.extends_above : s.extends_below) = 1;
						cand = s; dim = s.dim; stairs_dir = s.dir; bend_dir = s.bend_dir; sshape = s.shape; // copy fields from these stairs and extend down/up
						stack_conn    = 0; // not stacked - extended main stairs
						cand_is_valid = 1;
						if (ab) {is_at_top = 1;} // adding to the top, so this landing is the new top
						// clip stairs cube to the correct top or bottom floor to find the adjacent landing
						float const landing_z1(part.z2() - fc_thick + (ab ? -1.0 : 1.0)*window_vspacing);
						set_cube_zvals(s_bc, landing_z1, landing_z1+floor_thickness);
						find_and_merge_with_landing(interior->landings, s_bc, sshape, 1, ab); // num_floors=1
						
						if (allow_clip_walls) { // Note: conservative and not well tested, but this case is likely to be disabled anyway
							cube_t clip_cube(cand);
							if (ab) {clip_cube.z2() += window_vspacing;} else {clip_cube.z1() -= window_vspacing;}
							for (unsigned d = 0; d < 2; ++d) {clip_cube.d[s.dim][d] = ext_cube.d[s.dim][d];} // include padding
							for (unsigned d = 0; d < 2; ++d) {has_clipped_wall |= (uint8_t)subtract_cube_from_cubes(clip_cube, interior->walls[d], nullptr, 1);}
						}
						break; // done
					} // for s
				} // for ab
			}
			set_cube_zvals(cand, cand_z1, cand_z2);
			// check if any previously placed stairs span these zvals; since rooms are generally horizontally connected on each floor, there's a valid existing path
			bool have_spanning_stairs(0);

			if (!cand_is_valid) {
				for (stairwell_t const &s: interior->stairwells) {
					float const stairs_z2(s.z2() - (s.roof_access ? window_vspacing : 0.0)); // roof access flight doesn't count
					if (s.z1() <= cand_z1 && stairs_z2 >= cand_z2) {have_spanning_stairs = 1; break;}
				}
			}
			unsigned const cur_num_iters(have_spanning_stairs ? num_iters/2 : num_iters); // fewer iterations and no compact/cut stairs if we have existing spanning stairs
			std::set<pair<cube_t, unsigned>> edges_used;

			// iterations: 0-19: place in pri hallway, 20-39: place anywhere, 40-159: shrink size, 150-179: compact stairs, 180-199: allow cut walls
			for (unsigned n = 0; n < cur_num_iters && !cand_is_valid; ++n) { // make up to 200 tries to add stairs
				// clipped walls don't look right in some cases and may block hallways and rooms, use as a last resort; disable for houses since basement is optional anyway
				bool const allow_clip_walls(n >= 180 && !is_house && !is_prison()), use_pref_shared(n <= 2*iter_mult_factor);
				cube_t place_region(use_pref_shared ? pref_shared : shared); // use preferred shared area from primary hallway for bottom part for first 20 iterations

				if (n >= 4*iter_mult_factor && n < 16*iter_mult_factor && (n%iter_mult_factor) == 0) { // decrease stairs size slightly every 10 iterations, 12 times
					stairs_width -= 0.025*doorway_width; // 1.2*DW => 0.90*DW
					stairs_pad   -= 0.030*doorway_width; // 1.0*WD => 0.64*DW
					len_with_pad -= 0.230*doorway_width*(use_basement_stairs ? 1.2 : 1.0); // 6.0*DW => 3.24*DW / 2.689*DW ; basement can have steeper stairs
					// should this be get_min_front_clearance_inc_people()? that makes more sense, but using this value creates very steep stairs
					float const min_clearance(get_min_front_clearance()); // ensure the player can fit
					max_eq(stairs_width, min_clearance);
					max_eq(stairs_pad,   min_clearance);
				}
				if (min(place_region.dx(), place_region.dy()) < 1.5*len_with_pad) {dim = (place_region.dx() < place_region.dy());} // use larger dim
				else {dim = rgen.rand_bool();}
				// shrink place_region sides slightly to allow for the railing in office buildings and avoid z-fighting with the walls in houses
				place_region.expand_in_dim(!dim, -railing_pad);
				vector2d min_sz;
				min_sz[ dim] = len_with_pad;
				min_sz[!dim] = stairs_width;
				unsigned sel_room[2] = {0,0}; // {upper, lower}
				bool too_small(0);

				if (!allow_clip_walls) { // select rooms above and below to constrain the place region
					get_room_cands(interior->rooms, upper_part_ix, place_region, min_sz, check_private_rooms, room_cands[0]);
					if (room_cands[0].empty()) continue; // no valid upper room
					if (room_cands[0].size() > 1) {std::shuffle(room_cands[0].begin(), room_cands[0].end(), rand_gen_wrap_t(rgen));}
					bool success(0);

					for (cube_with_ix_t &rc1 : room_cands[0]) {
						get_room_cands(interior->rooms, lower_part_ix, rc1, min_sz, check_private_rooms, room_cands[1]); // set place_region=rc for use_pref_shared?
						if (room_cands[1].empty()) continue; // not a valid lower room
						cube_with_ix_t const &rc2(room_cands[1][rgen.rand()%room_cands[1].size()]);
						sel_room[0]  = rc1.ix;
						sel_room[1]  = rc2.ix;
						place_region = rc2;
						success      = 1;
						if (is_parking()) {place_region.expand_by_xy(-0.5*wall_thickness);} // shrink to prevent Z-fighting with inside surface of ext wall
						break; // success/done
					}
					if (!success) continue;
				}
				for (unsigned d = 0; d < 2; ++d) { // {x, y}
					float const stairs_sz(min_sz[d]), v1(place_region.d[d][0]), v2(place_region.d[d][1] - stairs_sz);
					if (v2 <= v1) {too_small = 1; break;}
					float &lo_edge(cand.d[d][0]);
					bool val_set(0);
					if (is_basement && bool(d) != dim && !(n&1)) {lo_edge = (rgen.rand_bool() ? v2 : v1); val_set = 1;} // choose edge of region, against a wall
					else if (!allow_clip_walls && bool(d) != dim) { // width dim
						bool edge(rgen.rand_bool());

						for (unsigned e = 0; e < 2; ++e) { // choose a side edge; if already tried, try the opposite edge
							if (edges_used.insert(make_pair(place_region, (2*dim + edge))).second) {lo_edge = (edge ? v2 : v1); val_set = 1; break;}
							edge ^= 1;
						}
					}
					if (!val_set) {lo_edge = rgen.rand_uniform(v1, v2);} // LLC
					cand.d[d][1] = min((lo_edge + stairs_sz), place_region.d[d][1]); // URC; clamp to avoid assert due to FP error
				} // for d
				if (too_small) continue;
				assert(place_region.contains_cube_xy(cand));
				bool const pri_hall_stairs(is_basement && !is_house && has_pri_hall() && pri_hall.z1() == ground_floor_z1 && dim == (hallway_dim == 1));
				// basement stairs placed in a first floor office building primary hallway should face the door
				if (pri_hall_stairs && pri_hall.contains_cube_xy(cand)) {stairs_dir = (pri_hall.get_center_dim(dim) < cand.get_center_dim(dim));}
				else {stairs_dir = rgen.rand_bool();} // the direction we move in when going up the stairs
				cube_t cand_test[2] = {cand, cand}; // {lower, upper} parts, starts on lower floor
				cand_test[0].z1() += 0.1*window_vspacing; cand_test[0].z2() -= 0.1*window_vspacing; // shrink to lower part
				cand_test[1].z1() += 1.1*window_vspacing; cand_test[1].z2() += 0.9*window_vspacing; // move to upper part
				bool bad_place(0), wall_clipped(0);

				if (!pri_hall_stairs && !allow_clip_walls) { // prefer stairs_dir that puts entrances and exits near doors
					point const stairs_center(cand.get_cube_center());
					point stair_ends[2] = {stairs_center, stairs_center};
					for (unsigned d = 0; d < 2; ++d) {stair_ends[d][dim] = cand.d[dim][d];}
					unsigned pref_dir(0); // two-bit mask

					for (room_t const &r : interior->rooms) {
						if (!r.intersects_xy(cand)) continue; // stairs not in this room

						for (unsigned d = 0; d < 2; ++d) {
							cube_t const &c(cand_test[d]);
							if (!r.intersects(c)) continue; // stairs not in this room
							if (!r.contains_cube_xy(c)) {bad_place = 1; break;} // stairs overlapping but not contained in this room
							get_doorways_for_room(r, r.z1(), doorways); // get interior doors on first floor of this room using thread safe function
							
							for (door_stack_t const &ds : doorways) {
								point const door_center(ds.get_cube_center());
								bool const dir(p2p_dist_xy_sq(door_center, stair_ends[0]) < p2p_dist_xy_sq(door_center, stair_ends[1]));
								pref_dir |= (1 << (dir ^ bool(d))); // record closer door dir, inverted for lower vs. upper entrances
							}
						} // for d
					} // for r
					if (bad_place) continue;
					if (pref_dir == 1) {stairs_dir = 0;} else if (pref_dir == 2) {stairs_dir = 1;} // pref agrees for all doors of both rooms
				}
				assert(cand.is_strictly_normalized());
				cand.expand_in_dim(dim, -stairs_pad); // subtract off padding
				if (!cand.is_strictly_normalized()) continue; // not enough space, likely because the player radius/front clearance is too large
				cand_test[ stairs_dir].d[dim][0] += stairs_pad; // subtract off padding on one side
				cand_test[!stairs_dir].d[dim][1] -= stairs_pad; // subtract off padding on one side

				for (unsigned d = 0; d < 2; ++d) { // {lower, upper}
					if (has_bcube_int(cand_test[d], interior->exclusion)) {bad_place = 1; break;} // bad placement
					cube_t cand_nopad(cand_test[d]), cand_pad(cand_nopad);
					cand_nopad.expand_in_dim( dim, -stairs_pad); // subtract off padding
					cand_pad  .expand_in_dim(!dim, railing_pad); // extra space at sides for railing
					// Note: check_walls is still needed to handle nested subrooms such as hospital room bathrooms
					if (!is_valid_stairs_elevator_placement(cand_pad, cand_nopad, stairs_pad, dim, !allow_clip_walls, check_private_rooms)) {bad_place = 1; break;} // bad place
					if (stairs_or_elevator_blocked_by_nested_room(cand_nopad, sel_room[!d])) {bad_place = 1; break;}
					// what about stairs intersecting bathrooms when allow_clip_walls=1? I've seen that happen once
					if (is_cube()) continue;
					// handle non-cube building; need to check both parts above and below, so clip our test cube to each part
					cube_t tb[2] = {cand_test[d], cand_test[d]};
					copy_zvals(tb[0], part);
					copy_zvals(tb[1], p   );
					if (!check_cube_within_part_sides(tb[0]) || !check_cube_within_part_sides(tb[1])) {bad_place = 1; break;} // check both top/bot parts
				} // for d
				if (bad_place) continue;

				if (allow_clip_walls) { // clip out walls around stairs
					cand_test[0].z1() = cand.z1(); cand_test[0].z2() = part.z2(); // lower
					cand_test[1].z1() = part.z2(); cand_test[1].z2() = part.z2() + window_vspacing - fc_thick; // upper (bot of ceiling for first floor on upper part)
					cand_test[ stairs_dir].d[dim][1] += 0.5*stairs_pad; // increase padding
					cand_test[!stairs_dir].d[dim][0] -= 0.5*stairs_pad; // increase padding

					for (unsigned e = 0; e < 2; ++e) {
						for (unsigned d = 0; d < 2; ++d) {wall_clipped |= subtract_cube_from_cubes(cand_test[e], interior->walls[d], nullptr, 1);} // clip_in_z=1
					}
					has_clipped_wall |= (uint8_t)wall_clipped;
				}
				// add walls around stairs if room walls were clipped or this is the basement; otherwise, make stairs straight with railings;
				// basement stairs only have walls on the bottom floor, so we set is_at_top=0; skip basement back stairs wall to prevent the player from getting stuck
				sshape        = (use_basement_stairs ? (stairs_shape)SHAPE_WALLED_SIDES : (wall_clipped ? (stairs_shape)SHAPE_WALLED : (stairs_shape)SHAPE_STRAIGHT));
				add_railing   = !wall_clipped; // or awlays? wall clipped stairs can sometimes have problems with railings, but sometimes they're needed
				is_at_top     = !is_basement; // no walls around the top of basement stairs; top floor signs are handled by the hallway case above
				cand_is_valid = 1;
				break;
			} // for n
			if (!cand_is_valid) { // no valid candidate found
				if (is_house && !is_basement) {interior->is_unconnected = 1;} // failed to connect house stacked parts
				continue;
			}
			if (!is_house) { // houses should only have one set of stairs on each floor
				// extend stairs to other floors of the part above and/or below if the parts don't already have stairs?
			}
			landing_t landing(cand, 0, 0, dim, stairs_dir, add_railing, sshape, 0, is_at_top, stack_conn); // roof_access=0
			landing.bend_dir = bend_dir;
			stairs_landing_base_t stairwell(landing);
			landing  .z1() = part.z2() - fc_thick; // only include the ceiling of this part and the floor of *p
			stairwell.z2() = part.z2() + window_vspacing - fc_thick; // bottom of ceiling of upper part; must cover z-range of upper floor for AIs and room object collisions
			interior->landings.push_back(landing);
			interior->stairwells.emplace_back(stairwell, 1); // num_floors=1

			if (is_retail) { // add intermediate landings/flights of stairs for multi-story retail areas
				for (unsigned n = 1; n < retail_floor_levels; ++n) { // skip first
					interior->landings.back().not_an_exit = 1; // all but the last (ground floor) landing are not exits
					if (n < 16) {interior->stairwells.back().not_an_exit_mask |= (1 << n);} // flag stairwells as well (used by building AI)
					landing.translate_dim(2, -window_vspacing); // shift down by a floor
					interior->landings.push_back(landing);
				}
			}
			// attempt to cut holes in ceiling of this part and floor of above part
			cube_t cut_cube(cand);
			cut_cube.z1() += fc_thick; // shrink to avoid clipping floors exactly at the base of the stairs
			subtract_cube_from_floor_ceil(cut_cube, interior->floors);
			subtract_cube_from_floor_ceil(cut_cube, interior->ceilings);

			// set has_stairs flags for containing rooms
			for (room_t &r : interior->rooms) {
				if (!r.intersects_no_adj(stairwell)) continue; // no stairs in this room
				// set the test point to the stairs entrance on the correct level using stairs_dir; the room contains stairs if it contains the stairs entrance
				point test_pt(stairwell.get_cube_center());

				if (part.contains_cube(r)) { // bottom of stairs
					test_pt[dim] = stairwell.d[dim][!stairs_dir];
					if (r.contains_pt_xy(test_pt)) {r.has_stairs |= (1U << min(num_floors-1U, 7U));} // top floor (clamp to 7 to avoid 8-bit overflow)
				}
				else if (p.contains_cube(r)) { // top of stairs
					test_pt[dim] = stairwell.d[dim][stairs_dir];
					if (r.contains_pt_xy(test_pt)) {r.has_stairs |= 1;} // bottom floor
				}
				else {cout << TXT(stairwell.str()) << TXT(r.str()) << TXT(part.str()) << TXT(p.str()) << endl; assert(0);} // something bad happened
			} // for r
			if (use_basement_stairs) { // add a basement door at the bottom of the stairs
				float const pos_shift((stairs_dir ? 1.0 : -1.0)*0.8*wall_thickness);
				door_t door(cand, dim, !stairs_dir, 0, 1); // open=0, on_stairs=1
				door.z2() -= fc_thick; // bottom of basement ceiling, not the floor above
				door.d[dim][stairs_dir] = door.d[dim][!stairs_dir] + pos_shift; // shift from the edge slightly into the stairwell, inserting the door as an end cap
				if (!stairs_dir) {door.translate_dim(dim, -pos_shift);} // why the asymmetry?
				door.translate_dim( dim, -0.2*pos_shift); // shift so that the door doesn't intersect the railing, covers the stairs overhang, and the top edge can't be seen
				door.expand_in_dim(!dim, -STAIRS_WALL_WIDTH_MULT*cand.get_sz_dim(dim)/NUM_STAIRS_PER_FLOOR); // shrink by stairs wall half width
				assert(door.is_strictly_normalized());
				interior->stairwells.back().stairs_door_ix = (int16_t)interior->doors.size(); // record door index gating stairs for AI navigation
				add_interior_door(door);
			}
			connected = 1;
			if (use_basement_stairs) break; // only need to connect one part for the basement
		} // for p
		if (!connected && is_basement) {interior->is_unconnected = 1;} // flag as unconnected if failed to connect basement with stairs
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
			cand_test.expand_in_z(-fc_thick); // shrink slightly in Z so that we don't intersect the original elevator *e
			cube_t const cand_test_nopad(cand_test);
			cand_test.d[e->dim][e->dir] += doorway_width*(e->dir ? 1.0 : -1.0); // add extra space in front of the elevator
			if (!p->contains_cube_xy(cand_test)) continue; // not enough space at elevator entrance
			bool const allow_clip_walls = 1; // optional
			// this check prevents us from extending both elevators in a back-to-back pair up or down at the same time, since they'll be too close to each other;
			// to work around this, we temporarily remove the adjacent elevator by mapping it to the building LLC
			cube_t orig_adj_elevator;

			if (e->adj_elevator_ix >= 0) {
				assert((unsigned)e->adj_elevator_ix < interior->elevators.size());
				elevator_t &adj(interior->elevators[e->adj_elevator_ix]);
				orig_adj_elevator = adj;
				adj.set_from_point(bcube.get_llc());
			}
			// Note: we've already added a doorway width of padding to the front of the elevator, so we don't need to add full padding in the call below
			bool const is_valid(is_valid_stairs_elevator_placement(cand_test, cand_test_nopad, doorway_width, e->dim, !allow_clip_walls, 1)); // check_private_rooms=1
			if (e->adj_elevator_ix >= 0) {interior->elevators[e->adj_elevator_ix].copy_from(orig_adj_elevator);} // restore original pos
			if (!is_valid)                                     continue; // bad placement
			if (!check_cube_within_part_sides(cand_test))      continue; // bad placement; do we need to check for clearance?
			if (has_bcube_int(cand_test, interior->exclusion)) continue; // bad placement

			if (allow_clip_walls) { // clip out walls around extended elevator
				cube_t clip_cube(cand_test);
				clip_cube.d[e->dim][e->dir] += 0.5*doorway_width*(e->dir ? 1.0 : -1.0); // add some extra wall clearance

				for (unsigned d = 0; d < 2; ++d) {
					holes.clear();
					has_clipped_wall |= (subtract_cube_from_cubes(clip_cube, interior->walls[d], (((bool)d == e->dim) ? &holes : nullptr)) << 1);

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
			if (is_below) { // extend below case
				if (is_retail_part(*p)) { // disable this floor on any elevators in the retail area
					for (unsigned n = 1; n < retail_floor_levels; ++n) {e->set_skip_floor(n);} // skip first
				}
				else { // shift skip_floors_mask up by the number of floors added to the bottom
					unsigned const num_floors_add(round_fp((e->z1() - extension.z1())/window_vspacing));
					e->skip_floors_mask <<= num_floors_add;
				}
			}
			float const shift((is_above ? -1.1 : 1.1)*fc_thick);
			min_eq(e->z1(), extension.z1()); max_eq(e->z2(), extension.z2()); // perform extension in Z
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
			// find all intersecting rooms and set has_elevator flag
			for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
				if (r->intersects(cand_test)) {++r->has_elevator;} // add another elevator
			}
		} // for p
	} // for e
}

bool building_t::are_parts_stacked(cube_t const &p1, cube_t const &p2) const {
	if (p1.z2() == p2.z1() && p1.contains_cube_xy(p2)) return 1; // p2 stacked on p1
	if (p2.z2() == p1.z1() && p2.contains_cube_xy(p1)) return 1; // p1 stacked on p2
	return 0;
}

// clip exterior ceiling for stairs, elevators, ramps, and skylights
bool building_t::clip_part_ceiling_for_stairs(cube_t const &c, vect_cube_t &out, vect_cube_t &temp) const {
	if (!interior) return 0;
	subtract_cubes_from_cube(c, interior->stairwells, out, temp);

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) { // handle elevators that span multiple parts
		if (e->z2() <= c.z1()) continue; // elevator shaft doesn't extend this high
		if (e->z2() == c.z2()) continue; // don't cut a hole in the roof of the building where the top of the elevator shaft ends
		subtract_cube_from_cubes(*e, out);
	}
	if (has_pg_ramp()) {subtract_cube_from_cubes(interior->pg_ramp, out);} // is this needed?
	for (cube_t const &skylight : skylights) {subtract_cube_from_cubes(skylight, out);} // check skylights
	return 1;
}

void building_interior_t::assign_door_conn_rooms(unsigned start_ds_ix) {
	assert(start_ds_ix <= door_stacks.size());

	for (auto d = door_stacks.begin()+start_ds_ix; d != door_stacks.end(); ++d) {
		if (d->get_for_closet() || d->get_backrooms()) continue; // excluded
		if (d->for_jail == 1) continue; // jail cell door should have already been assigned
		unsigned const dsix(d - door_stacks.begin());
		unsigned rooms_start(0), rooms_end(rooms.size());

		if (ext_basement_door_stack_ix >= 0) { // if we have an extended basement, determine whether or not this door is part of it and split the room range (optimization)
			if      ((int)dsix < ext_basement_door_stack_ix) {rooms_end   = ext_basement_hallway_room_id;} // main building door can skip extended basement rooms
			else if ((int)dsix > ext_basement_door_stack_ix) {rooms_start = ext_basement_hallway_room_id;} // extended basement door can skip main building rooms
			assert(rooms_start < rooms_end && rooms_end <= rooms.size());
		}
		point const door_center(d->get_cube_center());
		float const test_pt_shift(0.5*d->get_width());
		assert(d->first_door_ix < doors.size());
		point test_pts[2] = {door_center, door_center};

		for (unsigned s = 0; s < 2; ++s) { // for each side of the door
			if (d->on_stairs) {test_pts[s].z += (s ? d->dz() : 0.0);} // stairs door, test below and above
			else {test_pts[s][d->dim] += (s ? 1.0 : -1.0)*test_pt_shift;} // normal door, test front and back
		}
		for (unsigned s = 0; s < 2; ++s) { // for each side of the door
			point test_pt(test_pts[s]);
			int ds_room_ix(-1);

			for (unsigned r = rooms_start; r < rooms_end; ++r) {
				if (!rooms[r].contains_pt(test_pt)) continue;
				ds_room_ix = r;
				// only break if room does not contain the other side of the door; if it does, then this is likely a larger containing room and should be skipped
				if (!rooms[r].contains_pt(test_pts[!s])) break;
			}
			if (ds_room_ix == -1) { // adj room not found
				if (d->get_bldg_conn()) { // door connecting adjacent building with no room for this building on the other side
					ds_room_ix = 0; // set to 0 and hope it's unused; this can't be the first room, so the assert below won't fail
				}
				else { // can only happen with complex floorplan buildings where a wall ends exactly at a doorway so that neither room contains the point
					test_pt[!d->dim] = d->d[!d->dim][0]; // choose edge of door rather than center

					for (unsigned r = rooms_start; r < rooms_end; ++r) {
						if (rooms[r].contains_pt(test_pt)) {ds_room_ix = r; break;}
					}
					if (ds_room_ix < 0) {cout << "bad door " << door_center.str() << endl;} // TESTING
					assert(ds_room_ix >= 0);
					if (ds_room_ix < 0) {ds_room_ix = 0;} // don't fail if the assert above is commented out
				}
			}
			d->conn_room[s] = ds_room_ix;
		} // for s
		assert(d->conn_room[0] != d->conn_room[1]); // can't be connected to the same room on both sides

		for (unsigned dix = d->first_door_ix; dix < doors.size(); ++dix) {
			door_t &door(doors[dix]);
			if (!d->is_same_stack(door)) break; // moved to a different stack, done
			for (unsigned n = 0; n < 2; ++n) {door.conn_room[n] = d->conn_room[n];} // copy rooms to doors
		}
	} // for i
}

void building_t::create_two_story_tall_rooms(rand_gen_t &rgen) {
	if (!is_house || !interior)     return; // houses only, for now
	if (interior->rooms.size() < 6) return; // not enough rooms
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness), wall_thickness(get_wall_thickness());
	float const min_tall_room_sz(1.6*floor_spacing);

	// Note: wall trim top/bottom aren't drawn, but they generally aren't visible by the player standing on the bottom floor
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		room_t &room(*r);
		if (room.z1() != ground_floor_z1) continue; // ground floor only
		if (room.has_stairs || room.has_elevator || room.is_hallway || room.is_sec_bldg || room.is_single_floor) continue;
		unsigned const room_ix(r - interior->rooms.begin());
		if (has_int_garage && (int)room_ix == interior->garage_room)      continue; // no interior garages
		if (calc_num_floors(room, floor_spacing, floor_thickness) != 2)   continue; // two story rooms only
		if (has_attic() && room.contains_cube_xy(interior->attic_access)) continue; // don't make the attic access unreachable
		if (min(room.dx(), room.dy()) < min_tall_room_sz)                 continue; // room too small
		if (!is_room_adjacent_to_ext_door(room))                          continue; // only consider entrance rooms that may become living rooms
		
		// gather list of all connected doors on the upper floor
		float const upper_floor_zval_thresh(room.z1() + 1.5*floor_spacing); // anything above this is definitely on the second floor
		vect_door_stack_t &door_stacks(interior->door_stacks);
		vect_door_t &idoors(interior->doors);
		vector<unsigned> stack_ixs;

		for (unsigned i = 0; i < door_stacks.size(); ++i) {
			door_stack_t &ds(door_stacks[i]);
			if ( ds.not_a_room_separator())        continue; // skip basement and closet doors
			if (!ds.is_connected_to_room(room_ix)) continue;
			ds.set_mult_floor(); // counts as multi-floor (for drawing top edge), even if not extending to upper floor
			assert(ds.first_door_ix < idoors.size());
			idoors[ds.first_door_ix].set_mult_floor();
			if (ds.z2() > upper_floor_zval_thresh) {stack_ixs.push_back(i);} // add if extends to second floor
		} // for i
		if (stack_ixs.size() > 1) { // only need to check if there are multiple connecting doors, since a single door must connect as this room has no stairs
			// check for connected rooms on the upper floor that would become unreachable if their door was removed
			int first_room_ix(-1);
			bool is_disconnected(0);

			for (unsigned six : stack_ixs) { // check connectivity to other rooms connecting to the other side of this door stack
				unsigned const ds_room_ix(door_stacks[six].get_conn_room(room_ix));
				if (first_room_ix < 0) {first_room_ix = ds_room_ix;} // use the first room as a reference
				else if (!are_rooms_connected_without_using_room_or_door(first_room_ix, ds_room_ix, room_ix)) {is_disconnected = 1; break;}
			}
			if (is_disconnected) continue;
		}
		// replace doors with walls on upper floors
		float const ceil_zval(room.z1() + floor_spacing - fc_thick);
		vector<unsigned> dixs_to_remove;

		for (unsigned i : stack_ixs) {
			door_stack_t &ds(door_stacks[i]);
			unsigned const dix_to_remove(ds.first_door_ix+1); // second/upper door in the stack
			assert(dix_to_remove < idoors.size()); // first *and* second door must be valid
			unsigned const door_ix_end((i+1 == door_stacks.size()) ? idoors.size() : door_stacks[i+1].first_door_ix);
			assert(door_ix_end == ds.first_door_ix+2); // must be a stack of exactly two doors
			ds.z2() = ceil_zval; // clip to a single door height
			bool const dim(idoors[dix_to_remove].dim);
			cube_t wall(idoors[dix_to_remove]);
			wall.z1() = ceil_zval; // starts at the top of the lower door
			set_wall_width(wall, wall.get_center_dim(dim), 0.5*wall_thickness, dim);
			interior->walls[dim].push_back(wall);
			dixs_to_remove.push_back(dix_to_remove);

			for (auto &s : interior->stairwells) { // fix up basement stairs door index
				if (s.stairs_door_ix >= 0 && s.stairs_door_ix > (int)dix_to_remove) {--s.stairs_door_ix;}
			}
		} // for i
		// remove doors in reverse order since later indices will change; ext basement conn should be added later, so those door_ix values don't need to be updated
		for (auto i = dixs_to_remove.rbegin(); i != dixs_to_remove.rend(); ++i) {idoors.erase(idoors.begin() + *i);}
		// update door stack first_door_ix values
		for (unsigned dsi=0, six=0; dsi < door_stacks.size(); ++dsi) {
			unsigned &dix(door_stacks[dsi].first_door_ix);
			assert(dix >= six);
			dix -= six; // this number of doors was removed from earlier door stacks
			assert(dix < idoors.size());
			if (six < stack_ixs.size() && dsi == stack_ixs[six]) {++six;} // move to next update stack
		}
		// extend any wall adjacent to a shorter (single story) part upward by the floor thickness to fill the gap
		if (real_num_parts > 1) {
			cube_t const &part(get_part_for_room(room));

			for (unsigned dim = 0; dim < 2; ++dim) {
				for (unsigned dir = 0; dir < 2; ++dir) {
					float const wall_pos(room.d[dim][dir]);
					if (part.d[dim][dir] != wall_pos) continue; // not a part exterior edge

					for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
						if (*p == part) continue; // skip self
						if (p->d[dim][!dir] != wall_pos) continue; // not the adjacent part
						if (p->z2() >= part.z2()) continue; // not shorter
						// we found a shorter adjacent part; now create a new wall segment that covers the gap
						cube_t wall(room);
						set_cube_zvals(wall, ceil_zval, ceil_zval+fc_thick);
						wall.d[dim][!dir] = wall_pos;
						wall.d[dim][ dir] = wall_pos + (dir ? 1.0 : -1.0)*wall_thickness;
						interior->walls[dim].push_back(wall);
					} // for p
				} // for dir
			} // for dim
		}
		// remove the floor and ceiling between the two levels
		cube_t to_remove(room);
		to_remove.z1() += (floor_spacing - fc_thick); // first  floor ceiling
		to_remove.z2() -= (floor_spacing - fc_thick); // second floor floor
		subtract_cube_from_cubes(to_remove, interior->ceilings);
		subtract_cube_from_cubes(to_remove, interior->floors  );
		room.is_single_floor = 1;
		break; // at most one per house
	} // for room
}

void building_t::calc_room_ext_sides(room_t &room) const {
	cube_t const &part(parts[room.part_id]);
	room.ext_sides = 0;

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
		}
		if (is_exterior) {room.ext_sides |= (1 << d);}
	} // for d
}
unsigned building_t::add_room(cube_t const &room, unsigned part_id, unsigned num_lights, bool is_hallway, bool is_office, bool is_sec_bldg) {
	assert(interior);
	assert(room.is_strictly_normalized());
	room_t r(room, part_id, num_lights, is_hallway, is_office, is_sec_bldg);
	calc_room_ext_sides(r);
	if (check_skylight_intersection(room)) {r.set_has_skylight();}
	unsigned const room_id(interior->rooms.size());
	interior->rooms.push_back(r);
	return room_id;
}

void building_t::add_or_extend_elevator(elevator_t const &elevator, bool add) {
	assert(interior);
	if (add) {interior->elevators.push_back(elevator);}
	if (is_house || roof_type != ROOF_TYPE_FLAT) return; // sloped roof, not flat, can't add elevator cap
	float const window_vspacing(get_window_vspace());
	cube_t ecap(elevator);
	ecap.z1()  = elevator.z2();
	ecap.z2() += 0.25*window_vspacing; // set height
	if (!elevator.at_edge) {ecap.expand_by_xy(0.025*window_vspacing);}
	
	// check to see if the elevator is at the top of the building
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		if (p->z1() != elevator.z2()) continue; // not on top of the elevator
		if (p->intersects(ecap)) return; // part over elevator - should we add some sort of cap in this case, or block off the first floor, or add something to the interior?
	}
	if (check_skylight_intersection(ecap)) { // can this happen? probably not if elevator is placed correctly with this check
		if (add) {interior->elevators.back().under_skylight = 1;}
		return;
	}
	if (has_helipad && get_helipad_bcube().intersects_xy(ecap)) return; // check for helipad intersection
	remove_intersecting_roof_cubes(ecap);
	details.emplace_back(ecap, ROOF_OBJ_ECAP);
	max_eq(bcube.z2(), ecap.z2()); // extend bcube z2 to contain ecap
}

void building_t::remove_intersecting_roof_cubes(cube_t const &c) {
	vect_cube_t ac_to_remove;

	for (unsigned i = 0; i < details.size(); ++i) { // remove any existing objects that overlap ecap or stairs roof access
		auto &obj(details[i]);
		uint8_t const type(obj.type);
		// only remove some object types; may cause ducts to become disconnected from AC units
		if (type != ROOF_OBJ_BLOCK && type != ROOF_OBJ_AC && type != ROOF_OBJ_DUCT && type != ROOF_OBJ_ANT && type != ROOF_OBJ_WTOWER && type != ROOF_OBJ_SMOKESTACK &&
			type != ROOF_OBJ_SAT_DISH && type != ROOF_OBJ_TV_ANT) continue;
		if (!obj.intersects(c)) continue;
		if (type == ROOF_OBJ_AC) {ac_to_remove.push_back(obj);} // need to remove ducts connected to this AC unit

		if (type == ROOF_OBJ_BLOCK) { // see if there's a door associated with this block
			cube_t test_cube(obj);
			test_cube.expand_by_xy(get_wall_thickness());

			for (auto j = roof_tquads.begin(); j != roof_tquads.end(); ++j) {
				if (j->type == tquad_with_ix_t::TYPE_RDOOR2 && j->get_bcube().intersects(test_cube)) {roof_tquads.erase(j); break;} // there can be only one
			}
		}
		swap(obj, details.back());
		details.pop_back();
		--i; // wraparound okay
	} // for i
	for (cube_t const &ac : ac_to_remove) {remove_intersecting_roof_cubes(ac);} // remove connecting ducts as well
}
