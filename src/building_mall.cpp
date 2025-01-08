// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

bool const EXTEND_MALL_ELEVATOR_TO_CITY = 1;

extern bool camera_in_building;
using std::string;

string choose_store_name(unsigned store_type, unsigned item_category, rand_gen_t &rgen);
colorRGBA choose_pot_color(rand_gen_t &rgen);
colorRGBA choose_sign_color(rand_gen_t &rgen, bool emissive=0);
bool is_invalid_city_placement_for_cube(cube_t const &c);
void add_city_plot_cut(cube_t const &cut);
void add_city_ug_elevator_entrance(ug_elev_info_t const &uge);
void offset_hanging_tv(room_object_t &obj);
void rotate_obj_cube(cube_t &c, cube_t const &bc, bool in_dim, bool dir);
colorRGBA get_couch_color(rand_gen_t &rgen);
bool try_add_lamp(cube_t const &place_area, float floor_spacing, unsigned room_id, unsigned flags, float light_amt,
	vect_cube_t &cubes, vect_room_object_t &objects, rand_gen_t &rgen);

extern object_model_loader_t building_obj_model_loader;


float building_t::get_mall_floor_spacing(cube_t const &room) const { // special function that allows for larger than normal floor spacing
	assert(has_mall());
	return room.dz()/interior->num_extb_floors;
}
cube_t building_t::get_mall_center(cube_t const &room) const {
	bool const dim(interior->extb_wall_dim);
	float const ww_width(0.25*room.get_sz_dim(!dim));
	cube_t center(room); // center open area
	center.expand_in_dim(!dim, -1.0*ww_width); // short dim
	center.expand_in_dim( dim, -2.0*ww_width); // long dim
	return center;
}
void building_t::get_mall_open_areas(cube_t const &room, vect_cube_t &openings) const { // Note: applies to single floor malls as well
	bool const dim(interior->extb_wall_dim);
	cube_t const center(get_mall_center(room));
	float const length(center.get_sz_dim(dim)), width(center.get_sz_dim(!dim)), gap(0.75*width), half_gap(0.5*gap);
	unsigned const num_splits(max(1U, unsigned(round_fp(0.25*length/width))));
	float const step((length + gap)/num_splits);
	float pos(center.d[dim][0] - half_gap);
	openings.reserve(num_splits);

	for (unsigned n = 0; n < num_splits; ++n) {
		cube_t c(center);
		c.d[dim][0] = pos + half_gap;
		pos += step;
		c.d[dim][1] = pos - half_gap;
		openings.push_back(c);
	}
}

bool building_t::inside_mall_hallway(point const &pos) const {
	if (!has_mall()) return 0;

	for (auto r = interior->rooms.begin()+interior->ext_basement_hallway_room_id; r != interior->rooms.end(); ++r) {
		if (r->is_hallway && r->contains_pt(pos)) return 1;
	}
	return 0;
}
bool building_t::is_inside_mall_stores(point const &pos) const {
	if (!has_mall() || pos.z > ground_floor_z1) return 0;
	if (get_basement().contains_pt(pos))        return 0; // in basement, not mall
	// Note: inside_mall_hallway() check is needed to handle hallways on floors with no end stores; if this is too slow, the hallways can be stored inside interior
	return (interior->mall_info->store_bounds.contains_pt(pos) && !inside_mall_hallway(pos));
}

unsigned choose_one_center(unsigned num, rand_gen_t &rgen) {
	unsigned ix(num/2); // center value
	if (!(num & 1) && rgen.rand_bool()) {--ix;} // tie breaker if even
	return ix;
}

void building_t::add_mall_se_landing(cube_t const &c, bool is_escalator, bool se_dim, bool se_dir, bool ww_dir) {
	assert(c.is_strictly_normalized());
	interior->add_ceil_floor_pair(c, c.z1(), c.zc(), c.z2());
	interior->mall_info->landings.emplace_back(c, (8*is_escalator + 4*se_dim + 2*se_dir + ww_dir));
}

void building_t::setup_mall_concourse(cube_t const &room, bool dim, bool dir, rand_gen_t &rgen) {
	assert(interior);
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), door_width(get_doorway_width());
	float const floor_thickness(get_floor_thickness()), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness()), trim_thickness(get_trim_thickness());
	
	// add wall section below the door on lower floors, since the entrace door is on the top level
	float const room_end(room.d[dim][!dir]); // entrance end
	door_t const &door(interior->get_ext_basement_door()); // entrance door to mall
	assert(door.dim == dim);
	cube_t wall(door);
	set_cube_zvals(wall, room.z1()+fc_thick, door.z1()); // below the door
	wall.d[dim][!dir] = room_end; // doesn't exactly match the regular wall width, but not visible anyway
	wall.d[dim][ dir] = room_end + (dir ? 1.0 : -1.0)*wall_thickness; // extend into the room
	if (wall.dz() > 0.0) {interior->walls[dim].push_back(wall);}

	// handle upper floors
	unsigned const num_floors(interior->num_extb_floors);
	if (num_floors == 1) return; // single floor; nothing else to do
	vect_cube_t openings;
	get_mall_open_areas(room, openings);
	unsigned const num_openings(openings.size()), num_se_added((num_openings + 1)*(num_floors - 1));
	assert(num_openings > 0);
	reserve_extra(interior->stairwells, num_se_added);
	reserve_extra(interior->escalators, num_se_added);

	for (unsigned f = 1; f < num_floors; ++f) { // skip first floor
		float const z(room.z1() + f*floor_spacing), zc(z - fc_thick), zf(z + fc_thick), floor_below_z(zf - floor_spacing), floor_above_z(zf + floor_spacing);
		rand_gen_t se_rgen(rgen); // copy current rgen and create a temp one so that all floors have the same stacked stairs and escalator placement
		
		// add side walkways
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			cube_t ww(room);
			ww.d[!dim][!d] = openings.front().d[!dim][d];
			interior->add_ceil_floor_pair(ww, zc, z, zf);
		}
		for (unsigned n = 0; n <= num_openings; ++n) { // ends and gaps between openings
			bool const first(n == 0), last(n == num_openings);
			cube_t ww(openings.front());
			ww.d[dim][0] = (first ? room.d[dim][0] : openings[n-1].d[dim][1]);
			ww.d[dim][1] = (last  ? room.d[dim][1] : openings[n  ].d[dim][0]);
			interior->add_ceil_floor_pair(ww, zc, z, zf);
			// add stairs
			bool run_dir(first ? 1 : (last ? 0 : se_rgen.rand_bool())), side_dir(se_rgen.rand_bool());
			float ww_edge(ww.d[dim][run_dir]), ww_side(ww.d[!dim][side_dir]);
			float const stairs_len(1.5*floor_spacing), stairs_width(0.75*window_vspace);
			float rdir_sign(run_dir ? 1.0 : -1.0), side_sign(side_dir ? -1.0 : 1.0);
			cube_t stairs_bc;
			set_cube_zvals(stairs_bc, floor_below_z, zf); // top of floor below to top of current floor
			stairs_bc.d[ dim][ !run_dir] = ww_edge; // at walkway
			stairs_bc.d[ dim][  run_dir] = ww_edge + rdir_sign*stairs_len; // extend away from walkway
			stairs_bc.d[!dim][ side_dir] = ww_side + side_sign*wall_thickness;
			stairs_bc.d[!dim][!side_dir] = ww_side + side_sign*stairs_width;
			assert(stairs_bc.is_strictly_normalized());
			landing_t landing(stairs_bc, 0, 0, dim, !run_dir, 1, SHAPE_STRAIGHT, 0, 1, 0, 0, 1); // add_railing=1, roof_access=0, is_at_top=1, stack_conn=0, for_ramp=0, in_extb=1
			landing.num_stairs = round_fp(NUM_STAIRS_PER_FLOOR*floor_spacing/window_vspace);
			landing.in_mall    = 1;
			stairs_landing_base_t stairwell(landing);
			landing  .z1() = zc;
			stairwell.z2() = floor_above_z;
			interior->landings.push_back(landing);
			interior->stairwells.emplace_back(stairwell, 1); // num_floors=1

			if (f > 1) { // middle floor; add landing at bottom of stairs
				float const stairs_end(landing.d[dim][run_dir]);
				cube_t stairs_lc(landing);
				stairs_lc.translate_dim(2, -floor_spacing); // move to the floor below
				stairs_lc.d[!dim][ side_dir]  = ww_side;
				stairs_lc.d[!dim][!side_dir] += side_sign*wall_thickness; // extend out a bit further
				stairs_lc.d[ dim][!run_dir ]  = stairs_end;
				stairs_lc.d[ dim][ run_dir ]  = stairs_end + 1.2*rdir_sign*stairs_width;
				add_mall_se_landing(stairs_lc, 0, dim, run_dir, side_dir); // is_escalator=0
			}
			// add escalators
			if (!first && !last) {run_dir ^= 1; rdir_sign *= -1.0;} // opposite end except at concourse ends
			side_dir ^= 1; side_sign *= -1.0; // opposite side
			ww_edge = ww.d[ dim][run_dir ];
			ww_side = ww.d[!dim][side_dir];
			float const delta_z(floor_spacing + 2.0*trim_thickness), e_width(0.8*door_width); // make slightly taller so that top surface is above the railing bottom trim
			cube_t epair;
			set_cube_zvals(epair, floor_below_z, floor_above_z);
			epair.d[ dim][ !run_dir] = ww_edge - rdir_sign*door_width; // overlapping walkway
			epair.d[ dim][  run_dir] = ww_edge + rdir_sign*(1.0*delta_z + door_width); // extend away from walkway, 45 degree slope
			epair.d[!dim][ side_dir] = ww_side + side_sign*wall_thickness;
			epair.d[!dim][!side_dir] = ww_side + side_sign*2.2*e_width; // set width of the pair
			assert(epair.is_strictly_normalized());
			
			for (unsigned s = 0; s < 2; ++s) { // place two side-by-side escalators with opposite directions
				// extend 90% of floor thickness below; enough to hide building people animated feet, but not enough to clip through the ceiling below
				escalator_t e(epair, dim, !run_dir, (bool(s) ^ dim ^ run_dir ^ 1), door_width, delta_z, 0.9*floor_thickness, 1); // in_mall=1
				e.d[!dim][!s] = epair.d[!dim][s] + (s ? -1.0 : 1.0)*e_width;
				interior->escalators.push_back(e);
			}
			if (f > 1) { // middle floor; add landing at bottom escalators
				float const esc_end(epair.d[dim][run_dir]);
				cube_t esc_lc(epair);
				set_cube_zvals(esc_lc, zc-floor_spacing, zf-floor_spacing); // set ceil-floor zvals on the floor below
				esc_lc.d[!dim][ side_dir]  = ww_side;
				esc_lc.d[!dim][!side_dir] += side_sign*wall_thickness; // extend out a bit further
				esc_lc.d[ dim][!run_dir ]  = esc_end - rdir_sign*(door_width + wall_thickness); // under bottom escalator entrance, shifted a bit further back to avoid clipping
				esc_lc.d[ dim][ run_dir ]  = esc_end + 1.2*rdir_sign*stairs_width;
				add_mall_se_landing(esc_lc, 1, dim, run_dir, side_dir); // is_escalator=1
			}
		} // for n
	} // for f
	if (!openings.empty()) {
		// add elevator, which may extend up to ground level for city buildings;
		// note that this won't be added for single floor malls even though it would be valid if extending to ground level, since control flow doesn't get here
		unsigned const opening_ix(choose_one_center(num_openings, rgen));
		cube_t const opening(openings[opening_ix]);
		bool edir((num_openings & 1) ? rgen.rand_bool() : (opening_ix == num_openings/2)); // prefer closer to center; random if tied
		float const width(1.6*door_width), depth(width);
		elevator_t elevator(room, interior->ext_basement_hallway_room_id, dim, !edir, 0, 1, 1); // at_edge=0, interior_room=1, in_mall=1
		set_wall_width(elevator, opening.get_center_dim(!dim), 0.5*width, !dim); // set width

		for (unsigned d = 0; d < 3; ++d) { // try both sides, then wrap around to the first side again if both can't extend
			float const ww_edge(opening.d[dim][!edir] + (edir ? 1.0 : -1.0)*0.2*wall_thickness);
			// Note: elevator extends half a floor width below and above the room; is this okay, or can it clip through other objects?
			elevator.d[dim][!edir] = ww_edge; // front is adjacent to walkway
			elevator.d[dim][ edir] = ww_edge + (edir ? 1.0 : -1.0)*depth; // extend back away from walkway by depth

			if (EXTEND_MALL_ELEVATOR_TO_CITY && is_in_city && d < 2) { // extend elevator up to street level if there's space?
				cube_t entrance(elevator);
				entrance.expand_by_xy(0.25*wall_thickness + 0.1*depth); // account for city exterior entrance wall thickness == 0.1*depth
				cube_t test_cube(entrance);
				set_cube_zvals(test_cube, elevator.z2(), elevator.z2()+floor_spacing);
				test_cube.d[dim][!edir] -= (edir ? 1.0 : -1.0)*depth; // extend by depth in front of elevator

				if (!is_cube_city_placement_invalid(test_cube)) { // valid
					// Note: should extend window_vspace above, but elevator won't have a stop at this pos, unless parking garage is two floors; so we extend to floor_spacing
					elevator.z2() += floor_spacing;
					entrance.z2()  = elevator.z2();
					entrance.z1()  = ground_floor_z1; // at city level
					add_city_plot_cut(elevator);
					add_city_ug_elevator_entrance(ug_elev_info_t(entrance, (ground_floor_z1 + window_vspace), dim, !edir));
					interior->mall_info->city_elevator_ix = interior->elevators.size();
				}
				else { // try the other dir
					edir         ^= 1;
					elevator.dir ^= 1;
					continue;
				}
			}
			break; // success
		} // for d
		interior->elevators.push_back(elevator);
	}
	if (is_in_city) { // add skylight(s)?
		// must call is_cube_city_placement_invalid() and add_city_plot_cut()
	}
}

bool building_t::top_of_mall_elevator_visible(point const &camera_bs, vector3d const &xlate) const {
	if (!has_mall() || interior->mall_info->city_elevator_ix < 0)     return 0;
	if (camera_bs.z < ground_floor_z1 || camera_bs.z > bcube.z2())    return 0; // player under the ground or far above
	if (bcube.contains_pt(camera_bs))                                 return 0; // player in main building, not outside by elevator
	if (!interior->mall_info->store_bounds.contains_pt_xy(camera_bs)) return 0; // not above the mall / too far away (should end at the front of the building)
	elevator_t const &e(get_elevator(interior->mall_info->city_elevator_ix));
	if (e.contains_pt(camera_bs)) return 1; // player in elevator
	if ((e.d[e.dim][e.dir] < camera_bs[e.dim]) ^ e.dir)  return 0;
	cube_t elevator(e);
	elevator.z1() = ground_floor_z1; // at city level
	return camera_pdu.cube_visible(elevator + xlate); // VFC
}

bool building_t::is_cube_city_placement_invalid(cube_t const &c) const { // for mall skylights, elevator, etc.
	if (!is_basement_room_not_int_bldg(c, nullptr, 1)) return 1; // check for buildings above; no exclude, allow_outside_grid=1
	return is_invalid_city_placement_for_cube(c); // Note: city objects may not have been placed yet
}
bool building_t::is_store_placement_invalid(cube_t const &store) const {
	// here we don't need to check rooms of our own building, but we need to check other building basements;
	// we normally check for rooms outside the building's tile, but this doesn't seem to be a problem for malls; we still check city bounds for city malls
	return !is_basement_room_under_mesh_not_int_bldg(store, nullptr, !is_in_city); // allow_outside_grid=!is_in_city
}
void building_t::add_extb_room_floor_and_ceil(cube_t const &room) {
	float const fc_thick(get_fc_thickness());
	cube_t ceiling(room), floor(room);
	ceiling.z1() = room.z2() - fc_thick;
	floor  .z2() = room.z1() + fc_thick;
	interior->ceilings.push_back(ceiling);
	interior->floors  .push_back(floor);
	interior->basement_ext_bcube.assign_or_union_with_cube(room);
}
void add_mall_room_walls(cube_t const &room, float wall_thickness, bool dim, bool dir, bool at_mall_end, bool &has_adj_store, vect_cube_t walls[2]) {
	for (unsigned wdim = 0; wdim < 2; ++wdim) {
		for (unsigned wdir = 0; wdir < 2; ++wdir) {
			if (bool(wdim) != dim && bool(wdir) != dir && !at_mall_end) continue; // already have walls on this side
			if (bool(wdim) == dim && bool(wdir) == 0 && has_adj_store ) continue; // wall shared with adjacent store
			cube_t wall(room);
			set_wall_width(wall, room.d[wdim][wdir], 0.5*wall_thickness, wdim);
			walls[wdim].push_back(wall);
		}
	} // for dim
	has_adj_store = 1;
}

void building_t::add_mall_store(cube_t const &store, cube_t const &window_area, bool dim, bool dir, bool &has_adj_store, rand_gen_t &rgen) {
	bool const at_mall_end(window_area != store);
	float const floor_spacing(get_mall_floor_spacing()), window_vspace(get_window_vspace()), wall_thickness(get_wall_thickness()), fc_thick(get_fc_thickness());
	unsigned const room_ix(interior->rooms.size());
	room_t Room(store, basement_part_ix);
	Room.assign_all_to(RTYPE_STORE);
	Room.interior        = 2; // mark as extended basement
	Room.is_single_floor = 1;
	Room.set_interior_window();
	interior->get_extb_start_room().set_interior_window(); // mall concourse also has an interior window
	interior->rooms.push_back(Room);
	add_extb_room_floor_and_ceil(store);
	add_mall_room_walls(store, wall_thickness, dim, dir, at_mall_end, has_adj_store, interior->walls);
	// add door/window openings
	float const ceil_gap(get_floor_thick_val()*floor_spacing - fc_thick + get_mall_top_window_gap(floor_spacing, window_vspace));
	float const doorway_width(get_doorway_width()), wall_pos(store.d[!dim][!dir]);
	cube_t walls_cut(window_area);
	walls_cut.z1() += fc_thick;
	walls_cut.z2() -= ceil_gap;
	walls_cut.expand_in_dim(dim, -(0.25*window_vspace + 0.15*window_area.get_sz_dim(dim))*(at_mall_end ? 0.25 : 1.0)); // shrink, less at mall end stores
	set_wall_width(walls_cut, wall_pos, 2.0*wall_thickness, !dim); // add extra thickness
	subtract_cube_from_cubes(walls_cut, interior->walls[!dim], nullptr, 1); // no holes; clip_in_z=1
	// cut an opening in the center
	cube_t opening(walls_cut);
	set_wall_width(opening, walls_cut.get_center_dim(dim), 1.0*doorway_width, dim); // twice door width
	cube_t doorway(opening);
	set_wall_width(doorway, wall_pos, 0.5*wall_thickness, !dim);
	bool const closed(rgen.rand_float() < 0.1); // 10% of the time
	interior->mall_info->store_doorways.emplace_back(doorway, room_ix, !dim, dir, closed);
	
	// add window on each side of the doorway
	for (unsigned side = 0; side < 2; ++side) {
		cube_t window(walls_cut);
		window.d[dim][!side] = opening.d[dim][side];
		set_wall_width(window, wall_pos, 0.25*wall_thickness, !dim);
		assert(window.is_strictly_normalized());
		interior->int_windows.emplace_back(window, room_ix);
	}
}

void building_t::add_mall_stores(cube_t const &room, bool dim, bool entrance_dir, rand_gen_t &rgen) { // and bathrooms
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), wall_thickness(get_wall_thickness());
	float const min_width(4.0*window_vspace), max_width(9.0*window_vspace), min_depth(6.0*window_vspace), max_depth(8.0*window_vspace);
	float const fc_thick(get_fc_thickness()), doorway_width(get_doorway_width()), bathroom_width(4.0*window_vspace);
	unsigned const num_floors(interior->num_extb_floors);
	unsigned const max_bathrooms(max(1, round_fp(0.1*num_floors*room.get_sz_dim(dim)/room.get_sz_dim(!dim)))); // scale with number of stores
	unsigned num_bathrooms(0);
	bool at_stack_end[2] = {0, 0};
	vect_cube_t &side_walls(interior->walls[!dim]);
	vector<room_t> &rooms(interior->rooms);
	vector<unsigned> hall_stacks[2]; // one per end hallway
	interior->mall_info->store_bounds = room; // start with the mall concourse
	// pre-calculate depths of stores in each direction so that they vertically align correctly across stores;
	// use a consistent depth for each side + floor so that walls can be shared
	float depths[4] = {}; // for two sides and two ends
	for (unsigned n = 0; n < 4; ++n) {depths[n] = rgen.rand_uniform(min_depth, max_depth);}

	// pre-split walls into horizontal strips for each floor
	for (unsigned d = 0; d < 2; ++d) {
		vect_cube_t &walls(interior->walls[d]);
		vect_cube_t walls_to_add;
		assert(interior->extb_walls_start[d] < walls.size());

		for (auto w = walls.begin()+interior->extb_walls_start[d]; w != walls.end(); ++w) {
			for (unsigned f = 1; f < num_floors; ++f) { // skip first floor
				float const z(room.z1() + f*floor_spacing);
				if (w->z2() <= z || w->z1() >= z) continue; // doesn't span z
				walls_to_add.push_back(*w);
				w->z1() = walls_to_add.back().z2() = z;
			}
		} // for w
		vector_add_to(walls_to_add, walls);
	} // for d
	for (unsigned f = 0; f < num_floors; ++f) {
		bool const is_top_floor(f+1 == num_floors);
		unsigned const rooms_start(rooms.size()), store_walls_start(side_walls.size());
		bool added_bathrooms(0); // at most one pair per floor
		cube_t store, floor_bcube(room);
		store.z1() = floor_bcube.z1() = room .z1() + f*floor_spacing;
		store.z2() = floor_bcube.z2() = store.z1() +   floor_spacing;

		// place stores on each side of concourse
		for (unsigned d = 0; d < 2; ++d) { // sides of mall
			float const wall_pos(room.d[!dim][d]), depth(depths[d]);
			float pos(room.d[dim][0]), pos_end(room.d[dim][1]);
			float const middle(0.5*(pos + pos_end));
			// prevent exterior wall of store from clipping through parking garage wall
			if (entrance_dir) {pos_end -= wall_thickness;} else {pos += wall_thickness;}
			store.d[!dim][!d] = wall_pos;
			store.d[!dim][ d] = wall_pos + (d ? 1.0 : -1.0)*depth;
			bool has_adj_store(0);
			
			while (pos + min_width < pos_end) { // continue until we can't fit a min width room
				float const store_width(rgen.rand_uniform(min_width, max_width));
				float next_pos(pos);
				bool is_bathroom(0);

				if (!added_bathrooms && num_bathrooms < max_bathrooms && pos < middle && pos+store_width > middle) { // crosses the middle; make a bathroom
					next_pos   += bathroom_width;
					is_bathroom = 1;
				}
				else { // make a store
					next_pos += store_width;
				}
				if (next_pos + min_width > pos_end) {next_pos = pos_end;} // clamp to far end of mall, and prevent a narrow store
				store.d[dim][0] = pos;
				store.d[dim][1] = next_pos;
				assert(store.is_strictly_normalized());

				if (is_store_placement_invalid(store)) {
					is_bathroom = 0;
					store.d[dim][1] = next_pos = pos + min_width; // try min width store
					if (is_store_placement_invalid(store)) {pos = next_pos; has_adj_store = 0; continue;} // invalid, skip this store
				}
				if (is_bathroom) {
					bool const use_low_ceiling(1); // looks better since stalls don't extend very high
					bool const wm_first(rgen.rand_bool());
					float const sep_pos(0.5*(pos + next_pos)); // separator between men's and women's rooms
					cube_t door, bathrooms(store);
					set_cube_zvals(door, store.z1()+fc_thick, store.z1()+window_vspace-fc_thick);
					door.d[!dim][0] = door.d[!dim][1] = wall_pos;
					if (use_low_ceiling) {bathrooms.z2() = bathrooms.z1() + window_vspace;} // set normal ceiling height
					interior->mall_info->bathrooms.push_back(bathrooms);
					
					for (unsigned e = 0; e < 2; ++e) {
						room_type const rtype((bool(e) ^ wm_first) ? RTYPE_WOMENS : RTYPE_MENS);
						cube_t bathroom(bathrooms);
						bathroom.d[dim][!e] = sep_pos;
						room_t Room(bathroom, basement_part_ix);
						Room.assign_all_to(rtype);
						Room.interior        = 2; // mark as extended basement
						Room.is_single_floor = 1;
						Room.is_office       = 1; // required for creating office building bathroom
						rooms.push_back(Room);
						add_mall_room_walls(bathroom, wall_thickness, dim, d, 0, has_adj_store, interior->walls); // at_mall_end=0
						// add door
						set_wall_width(door, bathroom.get_center_dim(dim), 0.5*doorway_width, dim);
						door_t Door(door, !dim, d, 1); // open=1
						Door.make_auto_close();
						Door.set_mult_floor(); // counts as multi-floor (for drawing top edge)
						Door.rtype = rtype; // flag so that it has the correct sign
						add_interior_door(Door, 1); // is_bathroom=1
						cube_t walls_cut(door);
						walls_cut.expand_in_dim(!dim, wall_thickness);
						subtract_cube_from_cubes(walls_cut, interior->walls[!dim], nullptr, 1); // no holes; clip_in_z=1
					} // for e
					add_extb_room_floor_and_ceil(bathrooms); // add for both bathrooms together
					if (use_low_ceiling) {has_adj_store = 0;} // bathroom wall does not fully cover the height needed by an adjacent store wall
					added_bathrooms = 1; // on this floor
					++num_bathrooms;
				}
				else { // regular store
					add_mall_store(store, store, dim, d, has_adj_store, rgen); // window_area is full storefront
				}
				floor_bcube.union_with_cube(store);
				pos = next_pos;
			} // while
		} // for d
		// place a store at each end
		for (unsigned d = 0; d < 2; ++d) { // ends of mall
			bool const entrance_side(bool(d) == entrance_dir);
			if (entrance_side && is_top_floor) continue; // blocked by top floor entrance
			for (unsigned e = 0; e < 2; ++e) {store.d[!dim][e] = room.d[!dim][e];} // width of mall concourse
			float const wall_pos(room.d[dim][d]), depth(depths[d+2]);
			store.d[dim][!d] = wall_pos + (d ? 1.0 : -1.0)*(entrance_side ? -0.5*wall_thickness : 0.0); // shift slightly on entrance side
			store.d[dim][ d] = wall_pos + (d ? 1.0 : -1.0)*depth;
			if (is_store_placement_invalid(store)) continue;
			cube_t const window_area(store); // window area is clipped to the mall concourse

			// try to extend to the left and right to increase width
			for (unsigned e = 0; e < 2; ++e) {
				cube_t store_exp(store);
				store_exp.d[!dim][e] += (e ? 1.0 : -1.0)*depths[e];
				if (!is_store_placement_invalid(store_exp)) {store = store_exp;}
			}
			bool has_adj_store(0); // unused
			add_mall_store(store, window_area, !dim, d, has_adj_store, rgen);
			floor_bcube.union_with_cube(store);
		} // for d
		interior->mall_info->store_bounds.union_with_cube(floor_bcube);
		// add a narrower non-public hallway behind each row of stores
		bool const is_tall_room(floor_spacing > 1.1*window_vspace);
		unsigned const rooms_end(rooms.size());
		float const hall_width(2.0*doorway_width);
		cube_t hall(floor_bcube); // extends full length of the mall, including the end stores
		hall.expand_in_dim(dim, wall_thickness); // extend a bit further to allow connecting a parallel hallway between the two side hallways at each end
		hall.z2() = floor_bcube.z1() + window_vspace; // normal height
		cube_t side_halls[2];

		for (unsigned d = 0; d < 2; ++d) { // sides of mall
			float const wall_pos(floor_bcube.d[!dim][d]);
			if (room.d[!dim][d] == wall_pos) continue; // no stores added to this side, skip; excludes bathrooms
			hall.d[!dim][0]  = hall.d[!dim][1] = wall_pos; // move to outside of rooms
			hall.d[!dim][d] += (d ? 1.0 : -1.0)*hall_width;
			if (is_store_placement_invalid(hall)) continue;
			rooms.emplace_back(hall, basement_part_ix, 0, 1, 0, 0, 2); // num_lights=0, is_hallway=1, is_office=0, is_sec_floor=0, interior=1 (extb)
			add_extb_room_floor_and_ceil(hall);
			unsigned const hall_walls_start(side_walls.size());
			bool has_adj_store(0); // unused
			add_mall_room_walls(hall, wall_thickness, !dim, d, 1, has_adj_store, interior->walls); // at_mall_end=1 (to cover missing stores)
			side_halls[d] = hall;

			// split walls between stores and hallways into two half slices so that they can be textured differently; is there a way to skip drawing interior faces?
			for (auto w = side_walls.begin()+store_walls_start; w != side_walls.end(); ++w) {
				if (w->d[!dim][0] >= wall_pos || w->d[!dim][1] <= wall_pos) continue; // wrong wall
				w->d[!dim][bool(d) ^ (w < side_walls.begin()+hall_walls_start) ^ 1] = wall_pos; // set to half width
			}
			// connect to rooms with doors
			for (auto r = rooms.begin()+rooms_start; r != rooms.begin()+rooms_end; ++r) {
				if (!r->intersects(hall)) continue; // wrong side
				room_type const rtype(r->get_room_type(0));
				bool const is_bath(is_bathroom(rtype));
				//if (is_bath) continue; // no door to bathrooms?
				cube_t conn_room(*r);
				copy_zvals(conn_room, hall);
				cube_t const door_cut(add_ext_basement_door(conn_room, doorway_width, !dim, d, 1, (is_tall_room && !is_bath), rgen)); // is_end_room=1 (may be locked)
				subtract_cube_from_cubes(door_cut, interior->walls[!dim], nullptr, 1); // no holes; clip_in_z=1
				interior->doors.back().rtype = rtype; // flag door (in particular if bathroom) so that the correct sign is added
			}
		} // for d
		if (!side_halls[0].is_all_zeros() && !side_halls[1].is_all_zeros()) {
			// both side hallways are valid, try to add perpendicular hallways connecting them along the ends of the mall
			cube_t hall_union(side_halls[0]);
			hall_union.union_with_cube(side_halls[1]);

			for (unsigned d = 0; d < 2; ++d) { // each end
				float const dsign(d ? 1.0 : -1.0), wall_pos(hall_union.d[dim][d]);
				cube_t hall(hall_union);
				hall.d[dim][!d] = wall_pos;
				hall.d[dim][ d] = wall_pos + dsign*hall_width;
				if (is_store_placement_invalid(hall)) continue;
				unsigned const hall_room_ix(rooms.size());
				rooms.emplace_back(hall, basement_part_ix, 0, 1, 0, 0, 2); // num_lights=0, is_hallway=1, is_office=0, is_sec_floor=0, interior=1 (extb)
				add_extb_room_floor_and_ceil(hall);
				bool has_adj_store(0); // unused
				add_mall_room_walls(hall, wall_thickness, !dim, d, 1, has_adj_store, interior->walls); // at_mall_end=1 (to cover missing stores)
				
				for (unsigned e = 0; e < 2; ++e) { // add doors to each end connecting to the other hallways
					cube_t const door_cut(add_ext_basement_door(side_halls[e], doorway_width, dim, d, 0, 0, rgen)); // is_end_room=0 (unlocked), is_tall_room=0
					subtract_cube_from_cubes(door_cut, interior->walls[dim], nullptr, 1); // no holes; clip_in_z=1
				}
				if (at_stack_end[d]) continue; // done with this stack
				vector<unsigned> &stacks(hall_stacks[d]);
				if (stacks.empty() || get_room(stacks.back()).d[dim][!d] == wall_pos) {stacks.push_back(hall_room_ix);} // add if first or aligned to prev
				else {at_stack_end[d] = 1;} // flag so that the rare case of a missing middle hallway won't cause problems
			} // for d
		}
	} // for f
	if (MALL_FLOOR_HEIGHT == 2.0) { // connect back hallways with stairs and/or elevator; only valid for exactly N*floor (N=2) spacing
		for (unsigned d = 0; d < 2; ++d) { // each end
			vector<unsigned> const &stack(hall_stacks[d]);
			unsigned const stack_height(stack.size());
			if (stack_height < 2) continue; // need at least two floors to connect with elevator or stairs
			cube_t hall_span(get_room(stack.front())); // lowest level
			hall_span.union_with_cube(get_room(stack.back())); // highest level
			float const hall_center(hall_span.get_center_dim(!dim));
			// add an elevator connected to hall_span
			bool const elevator_side(rgen.rand_bool());
			float const hall_offset((elevator_side ? 1.0 : -1.0)*1.4*doorway_width), dsign(d ? 1.0 : -1.0), wall_pos(hall_span.d[dim][d]);
			cube_t elevator_bc;
			set_cube_zvals(elevator_bc, hall_span.z1(), (hall_span.z2() + floor_spacing - window_vspace)); // extend a full mall floor on the top
			set_wall_width(elevator_bc, (hall_center + hall_offset), 0.8*doorway_width, !dim);
			elevator_bc.d[dim][!d] = wall_pos;
			elevator_bc.d[dim][ d] = wall_pos + dsign*1.5*doorway_width; // depth
			assert(elevator_bc.is_strictly_normalized());

			if (!is_store_placement_invalid(elevator_bc)) {
				interior->elevators.emplace_back(elevator_bc, stack.front(), dim, !d, 1, 1, 2); // at_edge=1, interior_room=1, in_mall=2 (back hallway)
				cube_t elevator_cut(elevator_bc);
				elevator_cut.d[dim][!d] -= dsign*2.0*wall_thickness; // extend a bit into the hallway
				subtract_cube_from_cubes(elevator_cut, interior->walls[dim], nullptr, 0); // no holes; clip_in_z=0
				interior->basement_ext_bcube.assign_or_union_with_cube(elevator_bc); // required for elevator lights logic
				for (unsigned i : stack) {get_room(i).has_elevator = 1;} // elevator connects to all rooms in this stack
			}
			// add U-shaped stairs connected to hall_span
			cube_t stairs(hall_span);
			set_wall_width(stairs, (hall_center - hall_offset), 1.0*doorway_width, !dim); // shift opposite the elevator dir
			stairs.d[dim][!d] = wall_pos;
			stairs.d[dim][ d] = wall_pos + dsign*3.2*doorway_width; // depth
			assert(stairs.is_strictly_normalized());

			if (!is_store_placement_invalid(elevator_bc)) {
				unsigned const stairs_num_floors(2*stack_height - 1); // 2 floors per hallway except for the top
				stairwell_t stairwell(stairs, stairs_num_floors, dim, d, SHAPE_U);
				stairwell.in_mall = 2;

				for (unsigned f = 1; f < stairs_num_floors; ++f) { // skip bottom floor
					float const z(hall_span.z1() + f*window_vspace), zc(z - fc_thick), zf(z + fc_thick);
					landing_t landing(stairs, 0, f, dim, d, 1, SHAPE_U, 0, (f+1 == stairs_num_floors), 0, 0, 1);
					landing.in_mall = 2;
					if (!(f & 1)) {landing.not_an_exit = 1;}
					set_cube_zvals(landing, zc, zf);
					interior->landings.push_back(landing);
					if (f < 16 && (f & 1)) {stairwell.not_an_exit_mask |= (1 << f);} // odd floors are not exits
				} // for f
				interior->stairwells.emplace_back(stairwell);
				cube_t stairs_cut(stairs);
				stairs_cut.d[dim][!d] -= dsign*2.0*wall_thickness; // extend a bit into the hallway
				subtract_cube_from_cubes(stairs_cut, interior->walls[dim], nullptr, 0); // no holes; clip_in_z=0
				interior->basement_ext_bcube.assign_or_union_with_cube(stairs);
				for (unsigned i : stack) {get_room(i).has_stairs = 1;} // all rooms in this stack have stairs (on their first floors)
				cube_t bot_floor(stairs), top_ceil(stairs);
				bot_floor.z2() = hall_span.z1() + fc_thick;
				top_ceil .z1() = hall_span.z2() - fc_thick;
				interior->floors  .push_back(bot_floor);
				interior->ceilings.push_back(top_ceil );
			}
		} // for d
	}
}

void building_t::add_mall_stairs() { // connecting to the entrance door
	if (!has_mall()) return;
	room_t const &room(get_mall_concourse());
	door_t const &door(interior->get_ext_basement_door());
	bool const dim(interior->extb_wall_dim), dir(interior->extb_wall_dir);
	assert(door.dim == dim);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const room_id(interior->ext_basement_hallway_room_id);
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	// add stairs under the door if needed; maybe we should add to stairwells with SHAPE_FAN for building AI path finding to work?
	float const upper_floor_zval(room.z2() - floor_spacing + fc_thick), delta_z(door.z1() - upper_floor_zval);
	unsigned const num_steps(max(0, (int)ceil(NUM_STAIRS_PER_FLOOR*delta_z/get_floor_ceil_gap())));
	if (num_steps == 0) return; // no stairs needed
	float const step_height(delta_z/num_steps), step_len(2.0*step_height);
	float const wall_edge(room.d[dim][!dir]), dsign(dir ? 1.0 : -1.0), front_step_dist(dsign*step_len);
	cube_t stair(door);
	set_cube_zvals(stair, door.z1()-step_height, door.z1());
	stair.d[dim][!dir] = wall_edge; // starts at wall
	stair.d[dim][ dir] = wall_edge + 2.0*front_step_dist;

	for (unsigned n = 0; n < num_steps; ++n) { // top down
		stair.expand_in_dim(!dim, step_len);  // widen sides
		stair.d[dim][dir] += front_step_dist; // add length
		stair.intersect_with_cube_xy(room); // don't let stair extend outside mall in case door is close to a wall
		objs.emplace_back(stair, TYPE_STAIR, room_id, dim, dir, 0, 1.0, SHAPE_STAIRS_FAN);
		stair.translate_dim(2, -step_height);
	}
	// add railings against the wall
	cube_t railing;
	set_cube_zvals(railing, upper_floor_zval, door.z1());
	railing.d[dim][!dir] = wall_edge + dsign*1.0*wall_thickness;
	railing.d[dim][ dir] = wall_edge + dsign*1.6*wall_thickness;

	for (unsigned d = 0; d < 2; ++d) {
		railing.d[!dim][!d] = door .d[!dim][d] + (d ? 1.0 : -1.0)*1.5*wall_thickness;
		railing.d[!dim][ d] = stair.d[!dim][d] + (d ? 1.0 : -1.0)*1.5*wall_thickness;
		if (!room.contains_cube_xy(railing)) continue; // skip if outside mall in case door is close to a wall
		objs.emplace_back(railing, TYPE_RAILING, room_id, !dim, !d, 0, 1.0, SHAPE_CUBE, GOLD); // no balusters
	}
	stair.z2() = door.z2(); // set height
	interior->mall_info->ent_stairs = stairwell_t(stair, 1, dim, !dir, SHAPE_FAN);
}

// and back hallway stairs lights
void building_t::add_mall_lower_floor_lights(room_t const &room, unsigned room_id, unsigned lights_start, light_ix_assign_t &light_ix_assign) {
	// Note: lights should be in 2-3 rows, and no lights should be partially overlapping an opening in !dim
	bool const dim(interior->extb_wall_dim);
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness()), opening_pad(2.0*get_wall_thickness());
	vect_cube_t openings;
	get_mall_open_areas(room, openings);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_end(objs.size());
	assert(lights_start <= objs_end);

	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const zc(room.z1() + f*floor_spacing - fc_thick); // bottom of ceiling

		for (unsigned i = lights_start; i < objs_end; ++i) {
			room_object_t const &obj(objs[i]);
			if (obj.type != TYPE_LIGHT) continue; // should this ever fail?
			room_object_t light(obj);
			cube_t light_pad(light);
			light_pad.expand_in_dim(dim, opening_pad);
			bool skip(0);

			for (cube_t const &c : openings) {
				if (!c.intersects_xy(light_pad)) continue;
				if (c.contains_cube_xy(obj)) {skip = 1; break;} // fully over opening, skip
				bool const dir(c.get_center_dim(dim) < obj.get_center_dim(dim));
				light.translate_dim(dim, (c.d[dim][dir] - light_pad.d[dim][!dir])); // shift to not overlap opening
				break; // can only intersect one opening
			}
			if (skip) continue;
			light.translate_dim(2, (zc - obj.z2()));
			light.obj_id = light_ix_assign.get_ix_for_light(light);
			objs.push_back(light);
		} // for i
	} // for f
	for (stairwell_t const &s : interior->stairwells) {
		if (!s.in_mall || !s.is_u_shape() || s.not_an_exit_mask == 0) continue;
		unsigned const light_ix(light_ix_assign.get_next_ix()); // can this be reused?
		
		for (unsigned n = 0; n < s.num_floors; ++n) {
			if (!(s.not_an_exit_mask & (1 << n))) continue;
			add_U_stair_landing_lights(s, room_id, light_ix, (s.z1() + n*get_window_vspace()));
		}
	}
}

struct plant_loc_t : public sphere_t {
	bool upper_floor;
	colorRGBA color;
	plant_loc_t(point const &p, float r, colorRGBA const &c, bool uf=0) : sphere_t(p, r), upper_floor(uf), color(c) {}
};

// this is for the central mall concourse; store objects are added in add_mall_store_objs() below; treated as a single floor
unsigned building_t::add_mall_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, vect_cube_t &rooms_to_light) {
	bool const mall_dim(interior->extb_wall_dim);
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), fc_thick(get_fc_thickness()), doorway_width(get_doorway_width());
	float const wall_thickness(get_wall_thickness()), trim_thick(get_trim_thickness()), room_centerline(room.get_center_dim(!mall_dim)), pillar_hwidth(2.0*wall_thickness);
	float const railing_height(0.42*window_vspace), railing_top_bar_thick(0.04*window_vspace), vbar_hwidth(0.35*wall_thickness);
	float const plate_thickness(0.1*wall_thickness), bot_bar_thickness(0.4*wall_thickness), top_bar_thickness(0.5*wall_thickness);
	float const light_amt = 1.0; // fully lit, for now
	unsigned const num_floors(interior->num_extb_floors);
	cube_t const mall_center(get_mall_center(room));
	vect_cube_t openings, railing_cuts, railing_segs, temp, pillars, blockers; // blockers are on the first/ground floor only
	vector<plant_loc_t> plant_locs;
	point plant_loc(0.0, 0.0, zval);
	get_mall_open_areas(room, openings);
	unsigned const num_openings(openings.size());
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned num_elevators(0);
	colorRGBA const bar_color(LT_GRAY);
	colorRGBA pot_colors[3]; // plants at {stairs, escalators, pillars}
	for (unsigned n = 0; n < 3; ++n) {pot_colors[n] = choose_pot_color(rgen);}
	cube_t entrance_stairs_bcube(interior->mall_info->ent_stairs);
	entrance_stairs_bcube.expand_by_xy(doorway_width); // add plenty of clearance around the door for the stairs
	if (num_floors == 1) {blockers.push_back(entrance_stairs_bcube);} // entrance is on the first/ground floor
	cube_t pillar(room); // copy room zvals

	// gather railing cuts and plant locs
	for (stairwell_t const &s : interior->stairwells) {
		if (s.in_mall != 1) continue; // not in mall concourse
		cube_t cut(s);
		set_cube_zvals(cut, s.zc(), s.z2()-0.5*floor_spacing); // extends from floor at top of stairs half a floor up
		cut.expand_in_dim( s.dim, doorway_width); // add padding on both ends
		cut.expand_in_dim(!s.dim, 0.8*wall_thickness); // add space for railings
		railing_cuts.push_back(cut);
		room.has_stairs = 255; // should this be set earlier?
		// place plants to the sides of stairs at the bottom
		float const plant_radius(0.15*window_vspace);

		if (s.z1() < room.z1() + 0.5*floor_spacing) { // only on ground floor
			plant_loc[s.dim] = 0.9*s.d[s.dim][!s.dir] + 0.1*s.d[s.dim][s.dir]; // bottom end of stairs

			for (unsigned d = 0; d < 2; ++d) {
				plant_loc[!s.dim] = s.d[!s.dim][d] + (d ? 1.0 : -1.0)*1.25*plant_radius;
				plant_locs.emplace_back(plant_loc, plant_radius, pot_colors[0]);
			}
		}
		if (s.z1() <= zval) { // ground floor, add a blocker
			blockers.push_back(s);
			blockers.back().d[s.dim][!s.dir] += 1.25*(s.dir ? -1.0 : 1.0)*s.get_width(); // extend from bottom of stairs
			blockers.back().expand_in_dim(!s.dim, 2.25f*plant_radius); // add side clearance for plants
		}
	} // for s
	for (escalator_t const &e : interior->escalators) {
		if (!e.in_mall) continue;
		cube_t cut(e);
		set_cube_zvals(cut, e.zc(), e.z2()-0.5*floor_spacing); // extends from floor at top of escalator half a floor up
		cut.expand_in_dim(!e.dim, 0.5*wall_thickness); // add space for railings
		railing_cuts.push_back(cut);
		// add plants; escalators are placed in pairs, so we add one plant to each side
		bool const side(e.dir ^ e.move_dir ^ e.dim);
		float const plant_radius(0.16*window_vspace);
		plant_loc[!e.dim] = e.d[!e.dim][side] + (side ? 1.0 : -1.0)*1.25*plant_radius;

		if (e.z1() < room.z1() + 0.5*floor_spacing) { // only on ground floor
			plant_loc[ e.dim] = 0.9*e.d[e.dim][!e.dir] + 0.1*e.d[e.dim][e.dir]; // bottom end of escalator
			plant_locs.emplace_back(plant_loc, plant_radius, pot_colors[1]);
		}
		// add to upper walkway if not too close to railing
		if (plant_loc[!mall_dim] - plant_radius > mall_center.d[!mall_dim][0] && plant_loc[!mall_dim] + plant_radius < mall_center.d[!mall_dim][1]) {
			plant_loc[e.dim] = 0.1*e.d[e.dim][!e.dir] + 0.9*e.d[e.dim][e.dir]; // top end of escalator
			point loc(plant_loc);
			loc.z = e.z1() + floor_spacing; // one floor above
			plant_locs.emplace_back(loc, plant_radius, pot_colors[1], 1); // upper_floor=1
		}
		if (e.z1() <= zval) { // ground floor, add a blocker
			blockers.push_back(e);
			blockers.back().d[e.dim][!e.dir] += 2.0*(e.dir ? -1.0 : 1.0)*e.get_width(); // extend from bottom of escalator
			blockers.back().expand_in_dim(!e.dim, 2.25f*plant_radius); // add side clearance for plants
		}
	} // for e
	for (elevator_t const &e : interior->elevators) {
		if (e.in_mall != 1) continue; // not in mall concourse
		railing_cuts.push_back(e.get_bcube_padded(doorway_width));
		++num_elevators;
		// place clocks on back of elevator for each floor, and front if there's space; we know there's vertical space since elevators are only added if there are multiple floors
		bool const digital(rgen.rand_bool());
		float const floor_extra_spacing(floor_spacing - window_vspace);

		for (unsigned f = 0; f < num_floors; ++f) {
			float const z(zval + f*floor_spacing + 0.5*floor_extra_spacing);
			add_clock_to_cube(e, z, room_id, light_amt, e.dim, !e.dir, digital); // back

			if (floor_extra_spacing > 0.4*window_vspace) {
				cube_t place_cube(e);
				place_cube.d[e.dim][e.dir] += (e.dir ? 1.0 : -1.0)*0.5*wall_thickness; // account for front blocker
				add_clock_to_cube(place_cube, z, room_id, light_amt, e.dim, e.dir, digital); // front
			}
		} // for f
		blockers.push_back(e);
		blockers.back().d[e.dim][e.dir] += (e.dir ? 1.0 : -1.0)*e.get_length(); // add space in front
	} // for e
	room.has_elevator = num_elevators; // should this be set earlier?

	// add railings for mall landings
	for (cube_with_ix_t const &c : interior->mall_info->landings) {
		bool const is_escalator(c.ix & 8), se_dim(c.ix & 4), se_dir(c.ix & 2), ww_dir(c.ix & 1);
		float const zb1(c.z1() - trim_thick), zb2(c.z1() + 2.0*fc_thick + trim_thick), zt1(c.z1() + fc_thick + railing_height), zt2(zt1 + railing_top_bar_thick);
		float const inner_edge(c.d[!se_dim][!ww_dir]); // bordering the opening
		cube_t railing_area(c), vbar;
		if (is_escalator) {railing_area.d[se_dim][!se_dir] += (se_dir ? 1.0 : -1.0)*(doorway_width + 2.0*wall_thickness);} // clip off the area under the elevator entrance
		set_cube_zvals(vbar, zb2, zt1);
		set_wall_width(vbar, inner_edge, vbar_hwidth, !se_dim);

		for (unsigned d = 0; d < 2; ++d) { // add vertical bars
			set_wall_width(vbar, railing_area.d[se_dim][d], vbar_hwidth, se_dim);
			objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, EF_Z12); // skip top and bottom
		}
		for (unsigned dim = 0; dim < 2; ++dim) { // add railings and glass plates; somewhat duplicated with the openings code below, but difficult to factor out
			unsigned const skip_mask(get_skip_mask_for_xy(!dim));
			bool const dir((bool(dim) == se_dim) ? se_dir : !ww_dir);
			cube_t railing(railing_area);
			set_wall_width(railing, railing_area.d[dim][dir], plate_thickness, dim); // start with panel width
			// add bottom bar(s)
			cube_t bot_bar(railing);
			if (is_escalator && bool(dim) != se_dim) {bot_bar.d[se_dim][!se_dir] = c.d[se_dim][!se_dir];} // use full area for bottom bar
			bot_bar.expand_by_xy (bot_bar_thickness - plate_thickness); // additional expand from plate
			bot_bar.expand_in_dim(!dim, plate_thickness); // fill the corner notch
			set_cube_zvals(bot_bar, zb1, zb2);
			objs.emplace_back(bot_bar, TYPE_METAL_BAR, room_id, !dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, skip_mask);

			if (bool(dim) == se_dim) { // add second bottom bar to the other end, under the stairs or escalator, to close the ceiling-floor gap
				bot_bar.translate_dim(se_dim, (se_dir ? -1.0 : 1.0)*c.get_sz_dim(se_dim));
				objs.emplace_back(bot_bar, TYPE_METAL_BAR, room_id, !dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, skip_mask);
			}
			// add top bar
			cube_t bar(railing);
			bar.expand_by_xy (top_bar_thickness - plate_thickness); // additional expand from plate
			bar.expand_in_dim(!dim, plate_thickness); // fill the corner notch
			set_cube_zvals(bar, zt1, zt2);
			objs.emplace_back(bar, TYPE_METAL_BAR, room_id, !dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color); // can't skip drawing of ends if clipped by stairs
			// add glass
			cube_t panel(railing);
			set_cube_zvals(panel, zb2, zt1);
			objs.emplace_back(panel, TYPE_INT_WINDOW, room_id, dim, 0, 0, light_amt);
		} // for dim
		// add a vertical support next to the landing
		float const attach_pos(inner_edge + (ww_dir ? -1.0 : 1.0)*bot_bar_thickness);
		pillar.d[!se_dim][ ww_dir] = attach_pos;
		pillar.d[!se_dim][!ww_dir] = attach_pos + (ww_dir ? -1.0 : 1.0)*2.0*pillar_hwidth; // extend into the opening
		set_wall_width(pillar, (railing_area.d[se_dim][0] + (se_dir ? 0.1 : 0.9)*railing_area.get_sz_dim(se_dim)), pillar_hwidth, se_dim); // close to stairs/escalator
		pillars .push_back(pillar);
		blockers.push_back(pillar);
		// add railing cut where landing connects to walkway
		cube_t cut(c);
		cut.z2() += 0.5*floor_spacing; // extend above top of railing
		cut.d[!se_dim][ ww_dir] += (ww_dir ? 1.0 : -1.0)*1.00*doorway_width; // extend toward walkway
		cut.d[ se_dim][!se_dir] -= (se_dir ? 1.0 : -1.0)*0.75*doorway_width; // extend a bit under the stairs and escalator to avoid blocking the player with the railing
		railing_cuts.push_back(cut);
	} // for c

	// add vertical support pillars and fire extinguishers
	float fe_height(0.0), fe_radius(0.0);
	bool const add_fire_extinguishers(get_fire_ext_height_and_radius(window_vspace, fe_height, fe_radius));
	bool add_fe(rgen.rand_bool()); // random start on even vs. odd pillars
	vect_cube_t main_pillars; // used for food court trashcans

	for (cube_t const &opening : openings) {
		for (unsigned dir = 0; dir < 2; ++dir) { // each side of opening
			float const pillar_pos(opening.d[!mall_dim][dir] + (dir ? -1.0 : 1.0)*0.7*pillar_hwidth);
			set_wall_width(pillar, pillar_pos, pillar_hwidth, !mall_dim);
			set_wall_width(pillar, opening.get_center_dim(mall_dim), pillar_hwidth, mall_dim);
			main_pillars.push_back(pillar);
			railing_cuts.push_back(pillar);
			// place plants to the sides of pillars
			float const plant_radius(0.2*window_vspace), plant_clearance(1.25*plant_radius);
			plant_loc[!mall_dim] = pillar_pos;

			for (unsigned d = 0; d < 2; ++d) {
				plant_loc[mall_dim] = pillar.d[mall_dim][d] + (d ? 1.0 : -1.0)*plant_clearance;
				plant_locs.emplace_back(plant_loc, plant_radius, pot_colors[2]);
			}
			if (add_fire_extinguishers) { // maybe add fire extinguisher on pillar
				if (add_fe ^ bool(dir)) {
					bool const dim(!mall_dim), dir(pillar.get_center_dim(dim) < room_centerline);
					add_fire_ext(fe_height, fe_radius, zval, pillar.d[dim][!dir], pillar.get_center_dim(!dim), room_id, light_amt, dim, dir);
				}
			}
			blockers.push_back(pillar);
			blockers.back().expand_in_dim( mall_dim, 2.0*plant_clearance); // add clearance for plants
			blockers.back().expand_in_dim(!mall_dim, 2.0*fe_radius      ); // add clearance for fire extinguisherts
		} // for dir
		add_fe ^= 1; // swap every pair of pillars, alternating sides
	} // for opening
	vector_add_to(main_pillars, pillars);

	// add walkway railings
	for (unsigned f = 1; f < num_floors; ++f) { // skip first floor
		float const z(room.z1() + f*floor_spacing), zb1(z - fc_thick - trim_thick), zb2(z + fc_thick + trim_thick), zt1(z + railing_height), zt2(zt1 + railing_top_bar_thick);
		cube_t vbar;
		set_cube_zvals(vbar, zb2, zt1);

		for (cube_t const &opening : openings) {
			for (unsigned dim = 0; dim < 2; ++dim) {
				for (unsigned dir = 0; dir < 2; ++dir) {
					unsigned const skip_mask(get_skip_mask_for_xy(!dim));
					float const centerline(opening.d[dim][dir]);
					cube_t railing(opening);
					set_wall_width(railing, centerline, plate_thickness, dim); // start with panel width
					// add bottom bar - not clipped
					cube_t bot_bar(railing);
					bot_bar.expand_by_xy (bot_bar_thickness - plate_thickness); // additional expand from plate
					bot_bar.expand_in_dim(!dim, plate_thickness); // fill the corner notch
					set_cube_zvals(bot_bar, zb1, zb2);
					objs.emplace_back(bot_bar, TYPE_METAL_BAR, room_id, !dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, skip_mask);
					set_cube_zvals(railing, zb1, zt2); // set to full Z-range of railing for proper zval test in clipping
					subtract_cubes_from_cube(railing, railing_cuts, railing_segs, temp, 2); // check zval overlap

					for (cube_t const &r : railing_segs) {
						// add top bar
						cube_t bar(r);
						bar.expand_by_xy (top_bar_thickness - plate_thickness); // additional expand from plate
						bar.expand_in_dim(!dim, plate_thickness); // fill the corner notch
						set_cube_zvals(bar, zt1, zt2);
						objs.emplace_back(bar, TYPE_METAL_BAR, room_id, !dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color); // can't skip drawing of ends if clipped by stairs
						// add glass
						cube_t panel(r);
						set_cube_zvals(panel, zb2, zt1);
						assert(panel.is_strictly_normalized());
						objs.emplace_back(panel, TYPE_INT_WINDOW, room_id, dim, 0, 0, light_amt);

						// add vertical bars where railing was clipped
						for (unsigned d = 0; d < 2; ++d) {
							if (r.d[!dim][d] == railing.d[!dim][d]) continue; // not clipped
							set_wall_width(vbar, r.d[!dim][d], vbar_hwidth, !dim);
							set_wall_width(vbar, centerline,   vbar_hwidth,  dim);
							objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, EF_Z12); // skip top and bottom
						}
						// add vertical bars at intervals
						float const length(r.get_sz_dim(!dim));
						unsigned const num_vbars(round_fp(1.0*length/window_vspace));
						float const vbar_spacing(length/(num_vbars+1)); // bars at ends are already added
						float vbar_pos(r.d[!dim][0]);

						for (unsigned n = 0; n < num_vbars; ++n) {
							vbar_pos += vbar_spacing; // first bar is offset
							set_wall_width(vbar, vbar_pos,   vbar_hwidth, !dim);
							set_wall_width(vbar, centerline, vbar_hwidth,  dim);
							objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, EF_Z12); // skip top and bottom
						}
					} // for r
				} // for dir
			} // for dim
			// add corner bars; will need to add to both ends and both dims once there are gaps for stairs and escalators
			for (unsigned n = 0; n < 4; ++n) {
				bool const dirs[2] = {bool(n&1), bool(n&2)}; // x, y
				for (unsigned d = 0; d < 2; ++d) {set_wall_width(vbar, opening.d[d][dirs[d]], vbar_hwidth, d);}
				cube_t test_cube(vbar);
				test_cube.expand_by_xy(-0.75*vbar_hwidth); // shrink to allow a bit of railing overlap
				if (has_bcube_int(test_cube, railing_cuts)) continue; // skip if intersecting; shouldn't happen?
				objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, EF_Z12); // skip top and bottom
			}
		}
	} // for f

	// add objects for store doors, which are always between pairs of interior windows
	for (store_doorway_t const &d : interior->mall_info->store_doorways) {
		//room_t const &store(get_room(d.room_id));
		bool const dim(d.dim);
		// can either use room_id or d.room_id for these objects
		// add ceiling box where the gate would come down from
		cube_t cbox(d);
		cbox.z1() = d.z2() - 0.1*window_vspace;
		cbox.expand_in_dim( dim,  1.0*wall_thickness); // grow
		cbox.expand_in_dim(!dim, -0.5*wall_thickness); // shrink
		objs.emplace_back(cbox, TYPE_METAL_BAR, room_id, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, GRAY);
		
		if (d.closed) { // add closed gate
			cube_t gate(d);
			gate.z2()  = cbox.z1();
			gate.z1() += trim_thick;
			gate.expand_in_dim( dim, -0.4*d.get_sz_dim(dim)); // shrink
			gate.expand_in_dim(!dim, -0.5*wall_thickness); // shrink (same as cbox)
			objs.emplace_back(gate, TYPE_STORE_GATE, room_id, dim, 0, 0, light_amt, SHAPE_CUBE, LT_GRAY);
		}
	} // for d

	// add a fountain in the center of an opening with benches around it
	unsigned const NUM_MALL_BENCH_COLORS = 6;
	colorRGBA const mall_bench_colors[NUM_MALL_BENCH_COLORS] = {WHITE, LT_BROWN, DK_GRAY, LT_BROWN, colorRGBA(0.1, 0.4, 0.8), LT_BROWN}; // LT_BROWN becomes wood texture
	float const bench_height(0.42*window_vspace), bench_depth(0.23*window_vspace);
	int fountain_opening_ix(-1), fc_opening_ix(-1);

	if (!openings.empty() && building_obj_model_loader.is_model_valid(OBJ_MODEL_FOUNTAIN)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FOUNTAIN)); // W, D, H
		fountain_opening_ix = choose_one_center(num_openings, rgen);
		cube_t const &opening(openings[fountain_opening_ix]);
		float const max_radius(0.25*min(opening.dx(), opening.dy()));
		float height(0.75*floor_spacing*rgen.rand_uniform(0.9, 1.1)), radius(0.5*height*0.5*(sz.x + sz.y)/sz.z); // use average of width and depth for radius
		if (radius > max_radius) {height *= max_radius/radius; radius = max_radius;} // reduce size if radius is too large
		cube_t fbc(opening);
		set_cube_zvals(fbc, zval, zval+height);
		resize_around_center_xy(fbc, radius);
		unsigned const item_flags(rgen.rand()); // select a random sub_model_id
		objs.emplace_back(fbc, TYPE_BLDG_FOUNT, room_id, rgen.rand_bool(), rgen.rand_bool(), 0, light_amt, SHAPE_CYLIN, WHITE, item_flags); // random dim/dir
		blockers.push_back(fbc);
		// maybe add benches to either side of the fountain
		cube_t bench;
		set_cube_zvals(bench, zval, zval+bench_height); // set height
		colorRGBA const &bench_color(mall_bench_colors[rgen.rand() % NUM_MALL_BENCH_COLORS]);

		for (unsigned d = 0; d < 2; ++d) {
			if (rgen.rand_float() < 0.25) continue; // skip 25% of the time
			float const back_pos(fbc.d[mall_dim][d] + (d ? 1.0 : -1.0)*1.5*wall_thickness);
			bench.d[mall_dim][!d] = back_pos; // against the fountain
			bench.d[mall_dim][ d] = back_pos + (d ? 1.0 : -1.0)*bench_depth;
			set_wall_width(bench, fbc.get_center_dim(!mall_dim), 0.5*radius, !mall_dim);
			objs.emplace_back(bench, TYPE_BENCH, room_id, mall_dim, d, RO_FLAG_IN_MALL, light_amt, SHAPE_CUBE, bench_color);
			blockers.push_back(bench);
		} // for d
	}
	// trashcan setup for trashcans along walls and food court pillars
	bool const is_cylin_tcan(rgen.rand_bool());
	float const tcan_height((is_cylin_tcan ? 0.28 : 0.4)*window_vspace), tcan_radius(0.12*window_vspace);
	room_obj_shape const tcan_shape(is_cylin_tcan ? SHAPE_CYLIN : SHAPE_CUBE);
	colorRGBA const tcan_color(is_cylin_tcan ? LT_GRAY : colorRGBA(0.8, 0.6, 0.4)); // tan for cube trashcans
	// track available openings
	vector<unsigned> avail_openings;

	for (unsigned i = 0; i < num_openings; ++i) {
		if ((int)i != fountain_opening_ix) {avail_openings.push_back(i);}
	}
	// add a food court to one of the openings
	if (!avail_openings.empty()) {
		unsigned const aoix(rgen.rand() % avail_openings.size());
		fc_opening_ix = avail_openings[aoix];
		avail_openings.erase(avail_openings.begin() + aoix);
		cube_t place_area(openings[fc_opening_ix]);
		place_area.expand_by_xy(0.06*room.get_sz_dim(!mall_dim)); // allow a bit of overlap with the walkway bounds if there are multiple floors
		if (num_floors*num_openings < 4) {place_area.d[mall_dim][rgen.rand_bool()] = place_area.get_center_dim(mall_dim);} // half sized food court for small malls

		// add trashcans in food court next to main (not landing support) pillars on the ground floor
		for (cube_t const &p : main_pillars) {
			if (!p.intersects_xy(place_area)) continue;
			bool const pdir(room_centerline < p.get_center_dim(!mall_dim));
			cube_t tcan;
			set_cube_zvals(tcan, zval, zval+tcan_height); // set height
			float const wall_pos(p.d[!mall_dim][!pdir] + (pdir ? -1.0 : 1.0)*0.5*wall_thickness); // with a bit of padding
			tcan.d[!mall_dim][ pdir] = wall_pos; // against the wall
			tcan.d[!mall_dim][!pdir] = wall_pos + (pdir ? -1.0 : 1.0)*(is_cylin_tcan ? 2.0 : 1.6)*tcan_radius;
			set_wall_width(tcan, p.get_center_dim(mall_dim), tcan_radius, mall_dim);
			room_object_t const tcan_obj(tcan, TYPE_TCAN, room_id, !mall_dim, !pdir, RO_FLAG_IN_MALL, light_amt, tcan_shape, tcan_color);
			objs.push_back(tcan_obj);
			add_large_trashcan_contents(rgen, tcan_obj, room_id, light_amt);
			blockers.push_back(tcan);
		} // for p
		add_food_court_objs(rgen, place_area, zval, room_id, light_amt, blockers);
		interior->mall_info->food_court_bounds = place_area;
	}
	// add objects to remaining openings
	unsigned tree_opening_ix(0);
	if (!avail_openings.empty() && num_floors > 1) {tree_opening_ix = (rgen.rand() % avail_openings.size());} // select a random available opening for a tree

	for (unsigned i = 0; i < avail_openings.size(); ++i) {
		unsigned const oix(avail_openings[i]);
		assert(oix < openings.size());
		cube_t const &opening(openings[oix]);
		cube_t tree_avoid;
		
		// add palm or pine tree if more than one floor tall, 50% of the time and at least one available opening
		if (num_floors > 1 && (i == tree_opening_ix || rgen.rand_bool())) {
			bool const is_pine(0); // 0=palm, 1=pine; palm trees look better up close
			// make pine tree taller; palm tree has to be shorter since fronds are added to the top to increase the effective height
			float const height((is_pine ? 1.1 : 1.0)*rgen.rand_uniform(0.35, 0.4)*room.dz());
			cube_t tree_bc(point(opening.xc(), opening.yc(), zval));
			tree_bc.z2() += height;
			tree_bc.expand_by_xy(rgen.rand_uniform(0.20, 0.24)*height); // set radius
			objs.emplace_back(tree_bc, TYPE_TREE, room_id, 0, 0, RO_FLAG_IN_MALL, light_amt, SHAPE_CYLIN, choose_pot_color(rgen), (is_pine ? 1 : 0));
			tree_avoid = tree_bc;
			tree_avoid.expand_by_xy(0.1*height); // add extra padding for vases below
			//continue; // allow both trees and vases
		}
		// add vases/sculptures
		unsigned const num_vases((rgen.rand() & 3) + 1); // 1-4
		if (num_vases == 0) continue; // no vases
		float const opening_len(opening.get_sz_dim(mall_dim)), opening_width(opening.get_sz_dim(!mall_dim));
		float const place_area_hlen(min(0.2*opening_len, 0.5*opening_width)), place_area_hwidth(min(0.2*opening_len, 0.2*opening_width));
		point const opening_center(opening.xc(), opening.yc(), zval);
		cube_t place_area(opening);
		set_wall_width(place_area, opening_center[ mall_dim], place_area_hlen,    mall_dim);
		set_wall_width(place_area, opening_center[!mall_dim], place_area_hwidth, !mall_dim);
		colorRGBA const vase_color(gen_vase_color(rgen)); // consistent per opening
		unsigned const rand_ix(rgen.rand()); // consistent per opening
		point center;

		for (unsigned n = 0; n < num_vases; ++n) {
			float const height(window_vspace*rgen.rand_uniform(0.5, 1.0)), radius(window_vspace*rgen.rand_uniform(0.1, 0.25));
			
			for (unsigned N = 0; N < ((num_vases >= 3) ? 10U : 1U); ++N) { // 10 placement tries if at least 3 vases
				if (num_vases == 1) {
					center = opening_center;
				}
				else if (num_vases == 2) {
					if (n == 0) {center  = gen_xy_pos_in_area(place_area, radius, rgen);} // first vase
					else        {center += (opening_center - center);} // second vase diagonally opposite
				}
				else if (num_vases >= 3) {
					center[!mall_dim] = opening_center[!mall_dim];
					center[ mall_dim] = rgen.rand_uniform(place_area.d[mall_dim][0], place_area.d[mall_dim][1]);
				}
				cube_t vase;
				vase.set_from_sphere(center, radius);
				set_cube_zvals(vase, zval, zval+height);
				if (has_bcube_int(vase, blockers)) continue;
				if (!tree_avoid.is_all_zeros() && tree_avoid.intersects(vase)) continue;
				objs.emplace_back(vase, TYPE_VASE, room_id, 0, 0, RO_FLAG_IN_MALL, light_amt, SHAPE_CYLIN, vase_color);
				objs.back().obj_id = rand_ix;
				blockers.push_back(vase);
				break; // success
			} // for N
		} // for n
	} // for i
	// if there are bathrooms, add a water fountain between each pair
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_WFOUNTAIN)) {
		for (cube_t const &bathrooms : interior->mall_info->bathrooms) { // each bathroom pair
			bool const wf_dir(bathrooms.get_center_dim(!mall_dim) > room_centerline);
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_WFOUNTAIN)); // D, W, H
			float const height(0.25*window_vspace), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z);
			float const z1(bathrooms.z1() + fc_thick + 0.18*window_vspace), wall_pos(room.d[!mall_dim][wf_dir] + (wf_dir ? -1.0 : 1.0)*0.5*wall_thickness);
			cube_t wf;
			set_wall_width(wf, bathrooms.get_center_dim(mall_dim), hwidth, mall_dim);
			set_cube_zvals(wf, z1, z1+height);
			wf.d[!mall_dim][ wf_dir] = wall_pos;
			wf.d[!mall_dim][!wf_dir] = wall_pos - (wf_dir ? 1.0 : -1.0)*depth;
			objs.emplace_back(wf, TYPE_WFOUNTAIN, room_id, !mall_dim, wf_dir, 0, light_amt, SHAPE_CUBE);
		} // for bathrooms
	}
	// add potted plants
	float const plant_shift(1.1*get_flooring_thick());

	for (plant_loc_t &p : plant_locs) {
		p.pos.z += plant_shift; // move up slightly to avoid z-fighting of bottom when the dirt is taken
		cube_t const plant(get_cube_height_radius(p.pos, p.radius, 4.0*p.radius));
		unsigned const flags(p.upper_floor ? RO_FLAG_HAS_EXTRA : 0); // flag upper floor plants so that the bottom sides of leaves are drawn
		objs.emplace_back(plant, TYPE_PLANT, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN, p.color);
		set_obj_id(objs);
	}
	// place objects along mall-store walls on each floor, but not in front of glass
	cube_t place_area(room);
	place_area.expand_by_xy(-wall_thickness);
	vect_cube_t wall_blockers;
	unsigned const num_benches  (num_openings/2 + 1); // per floor per side
	unsigned const num_trashcans(num_openings/4 + 1); // per floor per side
	unsigned const num_pictures (num_openings/1 + 1); // per floor per side
	float const clearance(get_min_front_clearance_inc_people()), bench_hlen(0.5*window_vspace*rgen.rand_uniform(0.6, 0.8));
	float const tcan_end_pad(2.0*tcan_radius), bench_end_pad(2.0*bench_hlen);
	colorRGBA const &bench_color(mall_bench_colors[rgen.rand() % NUM_MALL_BENCH_COLORS]);

	for (unsigned f = 0; f < num_floors; ++f) {
		float const z(zval + f*floor_spacing);
		cube_t tcan, bench, picture;
		set_cube_zvals(tcan,  z, z+tcan_height ); // set height
		set_cube_zvals(bench, z, z+bench_height); // set height

		for (unsigned d = 0; d < 2; ++d) { // side
			float const dsign(d ? -1.0 : 1.0), wall_pos(place_area.d[!mall_dim][d]); // with a bit of padding
			// add benches
			bench.d[!mall_dim][ d] = wall_pos; // against the wall
			bench.d[!mall_dim][!d] = wall_pos + dsign*bench_depth;

			for (unsigned n = 0; n < num_benches; ++n) {
				for (unsigned N = 0; N < 10; ++N) { // 10 attempts to place bench
					float const pos(rgen.rand_uniform(place_area.d[mall_dim][0]+bench_end_pad, place_area.d[mall_dim][1]-bench_end_pad));
					set_wall_width(bench, pos, bench_hlen, mall_dim);
					cube_t test_cube(bench);
					test_cube.expand_in_dim( mall_dim, 0.5*bench_hlen);
					test_cube.expand_in_dim(!mall_dim, clearance);
					if (!is_valid_placement_at_mall_wall(test_cube, entrance_stairs_bcube, wall_blockers)) continue;
					wall_blockers.push_back(bench);
					wall_blockers.back().expand_by_xy(0.1*window_vspace); // don't place too close to other benches or trashcans
					objs.emplace_back(bench, TYPE_BENCH, room_id, !mall_dim, !d, RO_FLAG_IN_MALL, light_amt, SHAPE_CUBE, bench_color);
					break;
				} // for N
			} // for n
			// add trashcans
			tcan.d[!mall_dim][ d] = wall_pos; // against the wall
			tcan.d[!mall_dim][!d] = wall_pos + dsign*(is_cylin_tcan ? 2.0 : 1.6)*tcan_radius;

			for (unsigned n = 0; n < num_trashcans; ++n) {
				for (unsigned N = 0; N < 10; ++N) { // 10 attempts to place trashcan
					float const pos(rgen.rand_uniform(place_area.d[mall_dim][0]+tcan_end_pad, place_area.d[mall_dim][1]-tcan_end_pad));
					set_wall_width(tcan, pos, tcan_radius, mall_dim);
					cube_t test_cube(tcan);
					test_cube.expand_by_xy(tcan_radius + wall_thickness);
					if (!is_valid_placement_at_mall_wall(test_cube, entrance_stairs_bcube, wall_blockers)) continue;
					wall_blockers.push_back(tcan);
					wall_blockers.back().expand_by_xy(window_vspace); // don't place two nearby trashcans
					room_object_t const tcan_obj(tcan, TYPE_TCAN, room_id, !mall_dim, !d, RO_FLAG_IN_MALL, light_amt, tcan_shape, tcan_color);
					objs.push_back(tcan_obj);
					add_large_trashcan_contents(rgen, tcan_obj, room_id, light_amt);
					break;
				} // for N
			} // for n
			// add pictures
			picture.d[!mall_dim][ d] = wall_pos; // against the wall
			picture.d[!mall_dim][!d] = wall_pos + dsign*0.1*wall_thickness;

			for (unsigned n = 0; n < num_pictures; ++n) {
				float const hheight(0.5*rgen.rand_uniform(0.3, 0.4)*floor_spacing), hwidth(hheight*rgen.rand_uniform(1.5, 1.8)), pic_end_pad(1.0*hwidth);
				set_wall_width(picture, max((z + 0.5f*floor_spacing), (bench.z1() + hheight + 0.05f*floor_spacing)), hheight, 2); // set zvals; above benches

				for (unsigned N = 0; N < 20; ++N) { // 20 attempts to place picture
					float const pos(rgen.rand_uniform(place_area.d[mall_dim][0]+pic_end_pad, place_area.d[mall_dim][1]-pic_end_pad));
					set_wall_width(picture, pos, hwidth, mall_dim);
					cube_t test_cube(picture);
					test_cube.expand_in_dim(!mall_dim, 1.0*wall_thickness);
					test_cube.expand_in_dim( mall_dim, 2.0*wall_thickness);
					if (!is_valid_placement_at_mall_wall(test_cube, entrance_stairs_bcube, wall_blockers)) continue;
					objs.emplace_back(picture, TYPE_PICTURE, room_id, !mall_dim, !d, (RO_FLAG_NOCOLL | RO_FLAG_IN_MALL), light_amt); // picture faces dir opposite the wall
					objs.back().obj_id = rgen.rand(); // determines picture texture
					wall_blockers.push_back(test_cube);
				} // for N
			} // for n
			wall_blockers.clear();
		} // for d
	} // for f
	// add a reception desk and TV in front of end walls that have no store
	float const desk_width(1.5*window_vspace), desk_depth(0.5*desk_width), desk_height(0.32*window_vspace), centerline(room.get_center_dim(!mall_dim));
	cube_t desk;
	set_wall_width(desk, centerline, 0.5*desk_width, !mall_dim);

	for (unsigned f = 0; f < num_floors; ++f) {
		float const z(zval + f*floor_spacing);
		set_cube_zvals(desk, z, z+desk_height);

		for (unsigned d = 0; d < 2; ++d) { // ends
			if (f+1 == num_floors && bool(d) == !interior->extb_wall_dir) continue; // likely blocked by entrance stairs
			float const wall_pos(room.d[mall_dim][d]), d_sign(d ? -1.0 : 1.0), back_pos(wall_pos + d_sign*0.5*desk_depth);
			desk.d[mall_dim][ d] = back_pos - d_sign*desk_depth; // back, extended back more - for end store window coll test
			desk.d[mall_dim][!d] = back_pos + d_sign*desk_depth; // front
			if (has_bcube_int(desk, interior->int_windows) || has_bcube_int(desk, interior->mall_info->store_doorways)) continue; // too close to store door or window
			desk.d[mall_dim][ d] = back_pos; // back
			add_reception_desk(rgen, desk, mall_dim, !d, room_id, light_amt); // ignore return value

			if (building_obj_model_loader.is_model_valid(OBJ_MODEL_TV)) {
				vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
				float const tv_height(0.35*floor_spacing*rgen.rand_uniform(1.0, 1.2)), tv_hwidth(0.5*tv_height*sz.y/sz.z), tv_depth(tv_height*sz.x/sz.z);
				cube_t tv;
				tv.z1() = z + 0.3*floor_spacing;
				tv.z2() = tv.z1() + tv_height;
				set_wall_width(tv, centerline, tv_hwidth, !mall_dim);
				tv.d[mall_dim][ d] = wall_pos + d_sign*0.5*wall_thickness; // on the wall
				tv.d[mall_dim][!d] = tv.d[mall_dim][d] + (d ? -1.0 : 1.0)*tv_depth;
				objs.emplace_back(tv, TYPE_TV, room_id, mall_dim, !d, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), light_amt, SHAPE_CUBE, BLACK);
				offset_hanging_tv(objs.back());
				set_obj_id(objs);
			}
		} // for d
	} // for f

	// add pillars last so that we can check lights against them
	unsigned const pillars_start(objs.size());
	for (cube_t const &pillar : pillars) {objs.emplace_back(pillar, TYPE_OFF_PILLAR, room_id, !mall_dim, 0, 0, light_amt, SHAPE_CUBE, WHITE, EF_Z12);}
	return pillars_start;
}

bool building_t::is_valid_placement_at_mall_wall(cube_t const &c, cube_t const &entrance_stairs_bcube, vect_cube_t const &wall_blockers) const {
	if (entrance_stairs_bcube.intersects(c))              return 0; // too close to entrance stairs
	if (has_bcube_int(c, interior->mall_info->bathrooms)) return 0; // too close to bathroom doors and water fountain
	if (has_bcube_int(c, interior->int_windows) || has_bcube_int(c, interior->mall_info->store_doorways)) return 0; // too close to store door or window
	if (has_bcube_int(c, wall_blockers)) return 0;
	return 1;
}

void building_t::add_large_trashcan_contents(rand_gen_t &rgen, room_object_t const &tcan, unsigned room_id, float tot_light_amt) { // with trash
	vect_room_object_t &objs(interior->room_geom->objs);
	// place trash in trashcans
	bool const is_cylin(tcan.shape == SHAPE_CYLIN);
	float const tc_radius(0.5*tcan.get_sz_dim(!tcan.dim)); // cylinder radius or cube width
	unsigned const num_objs(rgen.rand() & 3); // 0-3

	for (unsigned n = 0; n < num_objs; ++n) {
		float trash_radius(tc_radius*rgen.rand_uniform(0.18, 0.3));
		cube_t trash;
		gen_xy_pos_for_round_obj(trash, tcan, trash_radius, 2.0*trash_radius, (trash_radius + (is_cylin ? 0.2 : 0.1)*tc_radius), rgen, 1); // place_at_z1=1
		trash.translate_dim(2, (trash_radius + (is_cylin ? 0.75 : 0.1)*tc_radius)); // shift up, more for cylinders
		colorRGBA const color(trash_colors[rgen.rand() % NUM_TRASH_COLORS]);
		objs.emplace_back(trash, TYPE_TRASH, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_SPHERE, color);
		set_obj_id(objs);
	} // for n
}

bool building_t::add_mall_table_with_chairs(rand_gen_t &rgen, cube_t const &table, cube_t const &place_area, colorRGBA const &chair_color,
	unsigned room_id, float tot_light_amt, bool dim, unsigned tid_tag, vect_cube_t &blockers)
{
	if (!place_area.contains_cube_xy(table) || has_bcube_int(table, blockers)) return 0;
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const table_obj_id(objs.size());
	room_object_t table_obj(table, TYPE_TABLE, room_id, 0, 0, RO_FLAG_IN_MALL, tot_light_amt, SHAPE_CUBE, WHITE, tid_tag); // tid_tag sets table texture
	objs.push_back(table_obj);
	set_obj_id(objs);
	// place chairs around the table
	float const table_len(table.get_sz_dim(dim)), table_width(table.get_sz_dim(!dim)), table_center(table.get_center_dim(!dim));
	unsigned const chairs_per_side(round_fp(table_len/table_width));
	float const chair_spacing(table_len/chairs_per_side);
	float const window_vspacing(get_window_vspace()), chair_height(0.45*window_vspacing), chair_hwidth(0.1*window_vspacing);
	bool const first_dir(rgen.rand_bool());

	for (unsigned D = 0; D < 2; ++D) {
		bool const dir(bool(D) ^ first_dir); // prevent bias

		for (unsigned n = 0; n < chairs_per_side; ++n) {
			if (D == 0 && rgen.rand_float() < 0.33) continue; // skip 33% of the time, but only on one side
			point chair_pos(0.0, 0.0, table.z1());
			chair_pos[!dim] = table_center + (dir ? -1.0f : 1.0f)*0.5*table_width;
			chair_pos[ dim] = table.d[dim][0] + (n + 0.5)*chair_spacing;
			// customized version of add_chair()
			chair_pos[!dim] += (dir ? -1.0 : 1.0)*rgen.rand_uniform(-0.5, 1.2)*chair_hwidth;
			cube_t chair(get_cube_height_radius(chair_pos, chair_hwidth, chair_height));
			if (!place_area.contains_cube_xy(chair) || has_bcube_int(chair, blockers)) continue;
			unsigned flags(RO_FLAG_IN_MALL);

			if (rgen.rand_float() < 0.05) { // fallen chair 5% of the time
				// rotate 90 degrees about back legs bottom, tilting backwards
				cube_t new_chair(chair);
				rotate_obj_cube(new_chair, chair, !dim, dir);

				if (!has_bcube_int(chair, blockers)) { // has space to fall
					chair  = new_chair;
					flags |= RO_FLAG_ON_FLOOR;
				}
			}
			objs.emplace_back(chair, TYPE_CHAIR, room_id, !dim, dir, flags, tot_light_amt, SHAPE_CUBE, chair_color, 1); // item_flags=1 to specify a plastic chair
			blockers.push_back(objs.back());
		} // for n
	} // for D
	blockers.push_back(table); // add the table last, so that it doesn't block its own chairs
	unsigned const place_obj_id(rgen.rand() & 7);

	switch (place_obj_id) { // TYPE_PHONE, TYPE_FOOD_BOX, TYPE_TRASH?
	case 0: place_bottle_on_obj(rgen, table_obj, room_id, tot_light_amt); break;
	case 1: place_dcan_on_obj  (rgen, table_obj, room_id, tot_light_amt); break;
	case 2: place_cup_on_obj   (rgen, table_obj, room_id, tot_light_amt); break;
	case 3: place_plate_on_obj (rgen, table_obj, room_id, tot_light_amt); break;
	case 4: if (rgen.rand_float() < 0.5) {place_pizza_on_obj (rgen, table_obj, room_id, tot_light_amt);} break; // less common
	case 5: if (rgen.rand_float() < 0.5) {place_banana_on_obj(rgen, table_obj, room_id, tot_light_amt);} break; // less common
	case 6: if (rgen.rand_float() < 0.1) {place_laptop_on_obj(rgen, table_obj, room_id, tot_light_amt);} break; // very rare
	case 7: if (rgen.rand_float() < 0.4) {place_eating_items_on_table(rgen, table_obj_id); break;} // less common
	} // default = place nothing
	return 1;
}

bool building_t::add_food_court_objs(rand_gen_t &rgen, cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, vect_cube_t const &blockers) {
	bool const dim(!interior->extb_wall_dim); // row dim runs perpendicular to mall dim for more shorter rows
	float const window_vspacing(get_window_vspace()), clearance(get_min_front_clearance_inc_people()), row_span(place_area.get_sz_dim(!dim));
	float const table_height(0.3*window_vspacing), table_width(0.4*window_vspacing), table_len_min(table_width), table_len_max(3.0*table_len_min), table_gap(0.25*table_width);
	float row_spacing(2.5*table_width);
	unsigned const num_rows(row_span/row_spacing); // floor
	if (num_rows == 0) return 0; // shouldn't happen?
	row_spacing = row_span/num_rows;
	float const col_start(place_area.d[dim][0]), col_end(place_area.d[dim][1]), max_cont_row_len(8.0*table_len_min);
	// find blockers intersecting the place area
	vect_cube_t fc_blockers;

	for (cube_t const &c : blockers) {
		if (!place_area.intersects_xy(c)) continue; // only consider blockers in the food court area
		fc_blockers.push_back(c);
		fc_blockers.back().expand_by_xy(clearance); // add space for walking around tables and chairs
	}
	unsigned const NUM_MALL_CHAIR_COLORS = 5;
	colorRGBA const mall_chair_colors[NUM_MALL_CHAIR_COLORS] = {WHITE, LT_GRAY, GRAY, ORANGE, LT_BROWN};
	colorRGBA const &chair_color(mall_chair_colors[rgen.rand() % NUM_MALL_CHAIR_COLORS]);
	unsigned const tid_tag(rgen.rand() + 1); // sets table texture; make nonzero to flag as a textured surface table
	cube_t table;
	set_cube_zvals(table, zval, zval+table_height);

	for (unsigned r = 0; r < num_rows; ++r) {
		float const row_pos(place_area.d[!dim][0] + (r + 0.5)*row_spacing);
		set_wall_width(table, row_pos, 0.5*table_width, !dim);
		float cur_row_len(0.0), pos(col_start + table_len_min*rgen.rand_float()); // starting point; use a random offset

		while (pos + 1.05*table_len_min < col_end) { // while we can place a table; slight bias to account for FP error
			float const table_len(rgen.rand_uniform(table_len_min, min(table_len_max, (col_end - pos))));
			table.d[dim][0] = pos;
			table.d[dim][1] = pos + table_len;
			float next_table_start(table_len + table_gap);
			cur_row_len += next_table_start;
			if (cur_row_len > max_cont_row_len) {next_table_start += 1.5*clearance; cur_row_len = 0.0;} // add a gap if run length gets too large
			pos += next_table_start;
			// what about tall tables with TYPE_BAR_STOOL, such as in coffe shops?
			add_mall_table_with_chairs(rgen, table, place_area, chair_color, room_id, tot_light_amt, dim, tid_tag, fc_blockers);
		} // while
	} // for r
	return 1;
}

// Note: room is non-const because the has_mirror flag may get set
void building_t::add_mall_store_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned &type_mask, light_ix_assign_t &light_ix_assign) {
	float const door_width(get_doorway_width()), floor_spacing(room.dz()), window_vspace(get_window_vspace());
	float const wall_thickness(get_wall_thickness()), wall_hthick(0.5*wall_thickness), fc_thick(get_fc_thickness()), pillar_width(2.0*wall_thickness);
	float const light_amt = 1.0; // fully lit, for now
	// get doorway bcube
	cube_t doorway;

	for (store_doorway_t const &d : interior->mall_info->store_doorways) {
		if (d.room_id != room_id) continue;
		assert(doorway.is_all_zeros()); // must be exactly one
		doorway = d;
	}
	assert(!doorway.is_all_zeros()); // must be found
	// select store type
	bool const dim(doorway.dy() < doorway.dx()), dir(room.get_center_dim(dim) < doorway.get_center_dim(dim)); // points from room center toward doorway; doorway wall
	bool const mall_dim(interior->extb_wall_dim), is_end_store(dim == mall_dim), tall_retail(floor_spacing > 1.5*window_vspace);
	float const room_len(room.get_sz_dim(dim)), room_width(room.get_sz_dim(!dim)), pillar_z2(room.z2() - fc_thick);
	room_t const &mall_room(get_mall_concourse());
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(room_id > 0); // can't be the first room
	room_t const &prev_room(interior->rooms[room_id-1]); // mall or store
	bool const on_new_floor(room.z1() != prev_room.z1()), on_new_row(room.d[dim][0] != prev_room.d[dim][0]);
	if (on_new_floor || on_new_row) {type_mask = 0;}
	unsigned const NUM_STORE_SELECT = 9;
	unsigned const store_selects[NUM_STORE_SELECT] = {STORE_CLOTHING, STORE_CLOTHING, STORE_BOOK, STORE_FURNITURE, STORE_PETS, STORE_APPLIANCE, STORE_RETAIL, STORE_RETAIL, STORE_RETAIL};
	unsigned const objs_start(objs.size());
	unsigned store_type(store_selects[rgen.rand() % NUM_STORE_SELECT]);
	cube_t const &fc_area(interior->mall_info->food_court_bounds);

	if (0 && !is_end_store && !fc_area.is_zero_area() && room.z1() == fc_area.z1()) { // food stores should be placed on ground floors near the food court
		// always place if contained in food court
		if (room.d[mall_dim][0] > fc_area.d[mall_dim][0] && room.d[mall_dim][1] < fc_area.d[mall_dim][1]) {store_type = STORE_FOOD;}
		// place 50% of the time if partially overlaps food court
		else if (room.d[mall_dim][0] < fc_area.d[mall_dim][1] && room.d[mall_dim][1] > fc_area.d[mall_dim][0] && rgen.rand_bool()) {store_type = STORE_FOOD;}
	}
	// bookstore and clothing stores are too expensive for the larger end stores, and pet stores/appliance stores should be small, so make them retail or furniture stores instead
	if (is_end_store && (store_type == STORE_BOOK || store_type == STORE_CLOTHING || store_type == STORE_PETS || store_type == STORE_APPLIANCE))
		{store_type = (rgen.rand_bool() ? STORE_FURNITURE : STORE_RETAIL);}
	// furniture stores should be larger, so make them book or clothing stores if small
	else if (store_type == STORE_FURNITURE && room_width < 0.8*room_len) {store_type = (rgen.rand_bool() ? STORE_BOOK : STORE_CLOTHING);}
	// mixed: 225FPS, 2991MB, 249ms
	// bookstores: 295FPS, 3463MB, 838ms
	// clothing stores: 163FPS, 2737MB, 19ms
	// retail stores: 300FPS, 2897MB, 177ms
	// don't place two pet stores in the same row; assign as retail or clothing instead
	if (store_type == STORE_PETS && (type_mask & (1U << store_type))) {store_type = (rgen.rand_bool() ? STORE_RETAIL : STORE_CLOTHING);}
	bool const is_retail(store_type == STORE_RETAIL);
	unsigned item_category(is_retail ? (rgen.rand() % NUM_RETAIL_CAT) : 0); // same category for each rack with equal probability
	if (is_end_store && is_retail && item_category == RETAIL_BOXED) {item_category = RETAIL_FOOD;} // make end retail stores food rather than boxes
	unsigned type_ix(is_retail ? (NUM_STORE_TYPES + item_category) : store_type);
	type_mask |= (1U << type_ix);
	// generate or select a store name
	string store_name;
	
	for (unsigned n = 0; n < 10; ++n) { // 10 attempts to generate a store name unique to this mall
		store_name = choose_store_name(store_type, item_category, rgen);
		bool is_duplicate(0);

		for (store_info_t const &i : interior->mall_info->stores) {
			if (store_name == i.name) {is_duplicate = 1; break;}
		}
		if (!is_duplicate) break; // unique name, done
	}
	//cout << store_name << endl; // TESTING
	interior->mall_info->stores.emplace_back(dim, dir, room_id, store_type, item_category, store_name);

	// place items
	if (1) { // add store name on sign above the entrance
		float const sign_z1(room.z1() + 0.7*floor_spacing), sign_height(0.3*window_vspace), sign_thick(wall_hthick);
		float const ext_wall_pos(mall_room.d[dim][!dir] + (dir ? 1.0 : -1.0)*(is_end_store ? 1.0 : 0.5)*wall_thickness);
		// stores to the sides of mall concourses have their doors centered, but stores on the end may not be centered on the door, but the mall concourse is
		float const door_center((is_end_store ? mall_room : room).get_center_dim(!dim));
		float sign_hwidth(0.25*sign_height*(store_name.size() + 2));
		min_eq(sign_hwidth, 0.5f*room_width); // can't be wider than the store
		cube_t sign;
		set_cube_zvals(sign, sign_z1, sign_z1+sign_height);
		sign.d[dim][!dir] = ext_wall_pos;
		sign.d[dim][ dir] = ext_wall_pos + (dir ? 1.0 : -1.0)*sign_thick;
		set_wall_width(sign, door_center, sign_hwidth, !dim);
		bool const emissive(0);
		unsigned const flags(RO_FLAG_LIT | RO_FLAG_NOCOLL | (emissive ? RO_FLAG_EMISSIVE : 0) | RO_FLAG_HANGING);
		colorRGBA const sign_color(choose_sign_color(rgen, emissive));
		objs.emplace_back(sign, TYPE_SIGN, interior->ext_basement_hallway_room_id, dim, dir, flags, light_amt, SHAPE_CUBE, sign_color); // always lit
		objs.back().obj_id = register_sign_text(store_name);
	}
	// add theft sensors to either side of the doorway for retail and clothing stores
	if (store_type == STORE_RETAIL || store_type == STORE_CLOTHING) {
		cube_t ts_area(doorway);
		ts_area.z2() = doorway.z1() + 0.6*window_vspace; // set height
		ts_area.d[dim][ dir] = doorway.d[dim][!dir] + (dir ? -1.0 : 1.0)*0.75*wall_thickness; // move slightly away from the doorway
		ts_area.d[dim][!dir] = ts_area.d[dim][ dir] + (dir ? -1.0 : 1.0)*0.22*window_vspace ; // extend into the store
		ts_area.expand_in_dim(!dim, 0.07*window_vspace);

		for (unsigned e = 0; e < 2; ++e) {
			cube_t ts(ts_area);
			ts.d[!dim][!e] = ts_area.d[!dim][e] + (e ? -1.0 : 1.0)*0.05*window_vspace; // set thickness
			objs.emplace_back(ts, TYPE_THEFT_SENS, room_id, !dim, e, 0, light_amt, SHAPE_CUBE, WHITE);
		}
	}
	cube_t blocked, room_area(room);
	room_area.expand_by_xy(-wall_hthick);

	// add checkout counter(s)/cash register(s) to the side of the door
	if (store_type == STORE_CLOTHING || store_type == STORE_RETAIL || store_type == STORE_BOOK) {
		bool const side(rgen.rand_bool());
		float const checkout_width(0.6*window_vspace), door_edge(doorway.d[!dim][side]);
		cube_t checkout_area(room);
		checkout_area.d[!dim][!side] = door_edge; // to the side of the doorway
		if (checkout_area.get_sz_dim(!dim) > checkout_width) {checkout_area.d[!dim][side] = door_edge + (side ? 1.0 : -1.0)*checkout_width;} // reduce width if needed
		checkout_area.expand_in_dim(!dim, -1.0*door_width); // add padding
		checkout_area.d[dim][!dir] = room.get_center_dim(dim); // only add to the front half of the store
		add_checkout_objs(checkout_area, zval, room_id, light_amt, objs_start, dim, side, (side ^ dim ^ 1));
		blocked = checkout_area;
		blocked.expand_in_dim(!dim, 2.0*door_width); // add extra space to both sides
	}
	if (store_type == STORE_FOOD || (store_type == STORE_RETAIL && item_category == RETAIL_FOOD)) { // maybe add wine racks
		float const width(0.4*window_vspace*rgen.rand_uniform(1.0, 1.5)), depth(0.16*window_vspace), height(0.5*window_vspace*rgen.rand_uniform(1.0, 1.5));

		for (unsigned yc = 0; yc < 2; ++yc) {
			for (unsigned xc = 0; xc < 2; ++xc) {
				if (rgen.rand_bool()) continue;
				bool const dim(rgen.rand_bool()), dir(dim ? yc : xc), dir2(dim ? xc : yc);
				cube_t c(point(room_area.d[0][xc], room_area.d[1][yc], zval)); // start at the room corner
				c.d[ dim][!dir ] += (dir  ? -1.0 : 1.0)*depth;
				c.d[!dim][!dir2] += (dir2 ? -1.0 : 1.0)*width;
				c.z2() += height;
				if (!blocked.is_all_zeros() && c.intersects_xy(blocked)) continue; // blocked; shouldn't happen?
				objs.emplace_back(c, TYPE_WINE_RACK, room_id, dim, !dir, 0, light_amt); // Note: dir faces into the room, not the wall
				set_obj_id(objs);
			} // for xc
		} // for yc
	}
	// add rows of shelves for retail and book stores, and rows of clothing to clothing stores, and rows of shelves to pet stores
	if (store_type == STORE_RETAIL || store_type == STORE_BOOK || store_type == STORE_CLOTHING || store_type == STORE_PETS) {
		cube_t place_area(room_area);
		place_area.expand_in_dim(dim, -0.5*door_width); // add extra padding in front and back for doors
		// simplified version of building_t::add_retail_room_objs() with no escalators, checkout counters, wall light, or short racks
		float const dx(place_area.dx()), dy(place_area.dy()), spacing(0.8), nom_aisle_width(1.5*door_width);
		unsigned const nx(max(1U, unsigned(spacing*dx/window_vspace))), ny(max(1U, unsigned(spacing*dy/window_vspace)));
		float const length(dim ? dy : dx), width(dim ? dx : dy), max_rack_width(0.5*window_vspace);
		unsigned const nrows((dim ? nx : ny)-1), nracks(max(2U, (dim ? ny : nx)/4));
		
		if (width > 4.0*nom_aisle_width && nrows > 0) { // can fit shelf racks
			float row_aisle_width(nom_aisle_width), aisle_spacing((width - row_aisle_width)/nrows), rack_width(aisle_spacing - row_aisle_width);
			assert(rack_width > 0.0);

			if (rack_width > max_rack_width) { // rack is too wide; widen the aisle instead
				rack_width      = max_rack_width;
				row_aisle_width = aisle_spacing - rack_width;
				aisle_spacing   = (width - row_aisle_width)/nrows;
			}
			float const rack_spacing((length - nom_aisle_width)/nracks), rack_length(rack_spacing - nom_aisle_width);
			assert(rack_length > 0.0);
			cube_t rack, pillar, pillar_area(room);
			set_cube_zvals(rack,   zval, (zval + SHELF_RACK_HEIGHT_FS*window_vspace));
			set_cube_zvals(pillar, zval, pillar_z2); // up to the ceiling
			pillar_area.expand_in_dim(dim, -0.25*room_len); // center 50% of room
			unsigned const style_id(rgen.rand()); // same style for each rack
			unsigned rack_id(0);

			for (unsigned n = 0; n < nrows; ++n) { // n+1 aisles
				float const rack_lo(place_area.d[!dim][0] + row_aisle_width + n*aisle_spacing), rack_center(rack_lo + 0.5*rack_width);
				rack.d[!dim][0] = rack_lo;
				rack.d[!dim][1] = rack_lo + rack_width;
				set_wall_width(pillar, rack_center, 0.5*pillar_width, !dim); // centered on the rack

				for (unsigned r = 0; r < nracks; ++r) {
					float const start(place_area.d[dim][0] + nom_aisle_width + r*rack_spacing);
					rack.d[dim][0] = start;
					rack.d[dim][1] = start + rack_length;
					if (!blocked.is_all_zeros() && rack.intersects_xy(blocked)) continue; // blocked
					cube_t test_cube(rack);
					test_cube.expand_in_dim( dim, 0.5*door_width); // add extra padding at ends
					test_cube.expand_in_dim(!dim, 1.0*door_width); // add extra padding at sides
					if (overlaps_other_room_obj(test_cube, objs_start)) continue; // blocked; not needed?

					if (store_type == STORE_BOOK) { // add bookcases
						add_row_of_bookcases(rack, zval, room_id, light_amt, !dim, 0); // place_inside=0
					}
					else if (store_type == STORE_CLOTHING) { // add clothes racks
						// Note: grouping into a rack object doesn't really help here because these are 3D models and must be drawn individually anyway
						float const hr_radius(0.007*window_vspace), frame_hwidth(1.2*hr_radius), frame_width(2.0*frame_hwidth);
						unsigned const flags(RO_FLAG_INTERIOR | RO_FLAG_IN_MALL | RO_FLAG_NOCOLL);
						// add an invisible collider around the clothes rack for player and AI collisions
						cube_t collider(rack);
						collider.z2() = zval + 0.5*window_vspace;
						set_wall_width(collider, rack_center, 0.2*rack_width, !dim);
						objs.emplace_back(collider, TYPE_COLLIDER, room_id, !dim, 0, (RO_FLAG_IN_MALL | RO_FLAG_INVIS), light_amt);
						unsigned const num_segs(2 + (room_id & 1)); // 2-3
						int const hanger_model_id(rgen.rand()), clothing_model_id(rgen.rand()); // consistent for each rack
						float const seg_len(rack_length/num_segs);

						for (unsigned n = 0; n < num_segs; ++n) { // split into segments
							// add hanger rod
							float const lo_val(start + n*seg_len);
							room_object_t hanger_rod(rack, TYPE_HANGER_ROD, room_id, !dim, 0, flags); // SHAPE_CUBE, even though it's a horizontal cylinder
							hanger_rod.d[dim][0] = lo_val;
							hanger_rod.d[dim][1] = lo_val + seg_len;
							hanger_rod.z1() = collider  .z2();
							hanger_rod.z2() = hanger_rod.z1() + 2.0*hr_radius;
							set_wall_width(hanger_rod, rack_center, hr_radius, !dim);
							// add a frame that holds the hanger rod; can construct from metal bars
							colorRGBA const frame_color(DK_GRAY);

							for (unsigned d = 0; d < 2; ++d) { // each end
								if (d == 1 && n+1 < num_segs) continue; // not the end, use frame from the next segment
								float const end_pos(hanger_rod.d[dim][d]);
								cube_t hbar(collider); // copy width from collider
								set_wall_width(hbar, hanger_rod.zc(), frame_hwidth, 2); // Z
								hbar.d[dim][ d] = end_pos;
								hbar.d[dim][!d] = end_pos + (d ? -1.0 : 1.0)*frame_width;
								objs.emplace_back(hbar, TYPE_METAL_BAR, room_id, 0, 0, flags, light_amt, SHAPE_CUBE, frame_color, 0); // draw all sides

								for (unsigned e = 0; e < 2; ++e) { // each side leg
									cube_t vbar(hbar);
									set_cube_zvals(vbar, rack.z1(), hbar.z1());
									vbar.d[!dim][!e] = hbar.d[!dim][e] + (e ? -1.0 : 1.0)*frame_width;
									objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, flags, light_amt, SHAPE_CUBE, frame_color, EF_Z12); // skip top and bottom
								}
							} // for d
							objs.push_back(hanger_rod); // must be added just before clothes
							// add hangers and hanging clothes
							unsigned const num_hangers(round_fp(15.0*rgen.rand_uniform(1.0, 1.5)*rack_length/(num_segs*window_vspace)));
							building_room_geom_t::add_hangers_and_clothing(window_vspace, num_hangers, flags, hanger_model_id, clothing_model_id, objs, rgen);
						} // for n
					}
					else if (store_type == STORE_PETS) {
						float const shelf_height(0.85*window_vspace), shelf_depth(0.22*window_vspace);
						cube_t center_wall(rack);
						center_wall.z2() = zval + shelf_height + fc_thick; // increase height above top of shelf wall anchors
						set_wall_width(center_wall, rack_center, 0.38*wall_thickness, !dim);
						add_shelves_along_walls(center_wall, zval, room_id, light_amt, !dim, store_type, shelf_height, shelf_depth, 0, rgen); // place_inside=0
						objs.emplace_back(center_wall, TYPE_PG_WALL, room_id, !dim, 0, (RO_FLAG_IN_MALL | RO_FLAG_ADJ_TOP), light_amt, SHAPE_CUBE); // draw top; or TYPE_STAIR_WALL?
					}
					else { // add retail shelf racks
						add_shelf_rack(rack, dim, style_id, rack_id, room_id, RO_FLAG_IN_MALL, item_category+1, 0, rgen); // add_occluders=0
					}
					if (n > 0 && n+1 < nrows && (n & 1)) { // use alternating/odd rows, skipping ends
						pillar.d[dim][0] = pillar.d[dim][1] = rack.d[dim][dir]; // end of back rack facing front of store - less likely to be blocked by checkout counter
						pillar.d[dim][dir] += (dir ? 1.0 : -1.0)*pillar_width;
						if (pillar_area.contains_cube_xy(pillar)) {add_retail_pillar(pillar, zval, room_id, tall_retail);}
					}
				} // for r
			} // for n
		}
	}
	if (store_type == STORE_BOOK) { // add bookcases along side walls
		cube_t bc_area(room_area);
		bc_area.expand_in_dim(dim, -0.1*room_len); // shrink ends
		add_row_of_bookcases(bc_area, zval, room_id, light_amt, !dim, 1); // place_inside=1
	}
	else if (store_type == STORE_CLOTHING) { // add shelves of TYPE_FOLD_SHIRT along walls; clothes racks are added above
		float const shelf_height(0.63*window_vspace), shelf_depth(0.25*window_vspace); // set height so that the top shelf is below the camera height
		add_shelves_along_walls(room_area, zval, room_id, light_amt, !dim, store_type, shelf_height, shelf_depth, 1, rgen); // place_inside=1
	}
	else if (store_type == STORE_FURNITURE || store_type == STORE_APPLIANCE) {
		// divide the store up into a 2D square grid of "rooms": bedrooms, dining rooms, living rooms, etc. and populate these
		float const trim_thick(get_trim_thickness()), room_pad(1.0*door_width), room_size(1.5*window_vspace), room_spacing(room_size + room_pad);
		float const size_delta(room_pad - wall_thickness - 2.0*trim_thick), pad_len(room_len + size_delta), pad_width(room_width + size_delta);
		unsigned const rooms_long(max(1U, unsigned(pad_len/room_spacing))), rooms_wide(max(1U, unsigned(pad_width/room_spacing)));
		bool const pdirs[2] = {rgen.rand_bool(), rgen.rand_bool()};
		float const rlen(pad_len/rooms_long), rwidth(pad_width/rooms_wide);
		cube_t door_blocker(doorway), div_area(room);
		door_blocker.expand_by_xy(0.2*window_vspace);
		div_area.expand_by_xy(0.5*room_pad - wall_hthick - trim_thick); // offset for for the shrink of rooms
		div_area.d[dim][dir] -= (dir ? 1.0 : -1.0)*wall_thickness; // shink for front window/wall clearance
		vect_cube_t blockers; // may be empty
		cube_t r; // sub-room bounds
		set_cube_zvals(r, room.z1(), room.z1()+window_vspace);
		unsigned types_used(0); // bit mask for appliance/plumbing stores

		for (unsigned row = 0; row < rooms_long; ++row) { // front to back
			for (unsigned col = 0; col < rooms_wide; ++col) { // side to side
				r.d[ dim][0] = div_area.d[ dim][0] + row*rlen;
				r.d[ dim][1] = r       .d[ dim][0] + rlen;
				r.d[!dim][0] = div_area.d[!dim][0] + col*rwidth;
				r.d[!dim][1] = r       .d[!dim][0] + rwidth;
				r.expand_by_xy(-0.5*room_pad); // shrink to add padding between rooms for people to walk
				if (r.intersects(door_blocker)) continue; // leave the room(s) in front of the entrance empty
				unsigned const start_ix(objs.size());
				blockers.clear();

				if (row > 0 && col > 0 && row+1 < rooms_long && col+1 < rooms_wide) { // skip edge sub-rooms
					// place pillars at alternating rows and columns, from the side closest to the center if the number of rows/columns is even
					if (bool(row & 1) == (bool(rooms_long & 1) || pdirs[dim]) && bool(col & 1) == (bool(rooms_wide & 1) || pdirs[!dim])) {
						cube_t pillar(r);
						pillar.z2() = pillar_z2;
						for (unsigned d = 0; d < 2; ++d) {pillar.d[d][!pdirs[d]] = r.d[d][pdirs[d]] + (pdirs[d] ? -1.0 : 1.0)*pillar_width;}
						add_retail_pillar(pillar, zval, room_id, tall_retail);
						blockers.push_back(pillar); // needed for bed placement
					}
				}
				if (store_type == STORE_FURNITURE) {
					room_t sub_room(r, room.part_id);
					sub_room.assign_all_to(RTYPE_STORE);
					// set open_wall flags; required for drawing backs of bookcases and avoiding tall desks
					if (row   > 0         ) {sub_room.mark_open_wall( dim, 0);}
					if (row+1 < rooms_long) {sub_room.mark_open_wall( dim, 1);}
					if (col   > 0         ) {sub_room.mark_open_wall(!dim, 0);}
					if (col+1 < rooms_wide) {sub_room.mark_open_wall(!dim, 1);}
					sub_room.mark_open_wall(dim, dir); // front wall is glass and always considered open
					point const room_center(r.xc(), r.yc(), zval);
					colorRGBA const &chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
					unsigned const NUM_RTYPES = 6;
					unsigned const rtypes[NUM_RTYPES] = {RTYPE_BED, RTYPE_BED, RTYPE_LIVING, RTYPE_LIVING, RTYPE_OFFICE, RTYPE_DINING};
					unsigned const rtype(rtypes[rgen.rand() % NUM_RTYPES]);

					if (rtype == RTYPE_BED) { // Note: light_ix_assign is needed for closet lights, but these aren't added here
						unsigned const bed_ix(objs.size());

						if (add_bedroom_objs(rgen, sub_room, blockers, chair_color, zval, room_id, 0, light_amt, start_ix, 1, 0, 1, light_ix_assign)) {
							for (auto i = objs.begin()+bed_ix+1; i != objs.end(); ++i) { // set has_mirror flag if a dresser mirror was placed
								if (i->type == TYPE_DRESS_MIR) {room.set_has_mirror();}
							}
							// add wall at head of bed
							room_object_t const &bed(objs[bed_ix]); // should be the first object placed
							assert(bed.type == TYPE_BED);
							bool const wdim(bed.dim), wdir(bed.dir);
							cube_t wall(sub_room);
							float const wall_pos(sub_room.d[wdim][wdir]);
							wall.d[wdim][!wdir] = wall_pos;
							wall.d[wdim][ wdir] = wall_pos + (wdir ? 1.0 : -1.0)*0.5*wall_thickness; // shift outward from room
							cube_t wall_area(room);
							wall_area.d[dim][dir] -= (dir ? 1.0 : -1.0)*0.5*rlen; // shink to avoid blocking front window

							// only if interior; skip if overlaps the room wall; check for back hallway doorway
							if (wall_area.contains_cube_xy_exp(wall, wall_thickness) && !is_cube_close_to_doorway(wall, sub_room, 0.0, 1, 1)) {
								objs.emplace_back(wall, TYPE_PG_WALL, room_id, wdim, wdir, (RO_FLAG_IN_MALL | RO_FLAG_ADJ_TOP), light_amt); // draw top
							
								if (rgen.rand_float() < 0.67) { // add a picture on the wall 67% of the time
									float const height(window_vspace*rgen.rand_uniform(0.28, 0.42));
									float const width(min(height, 0.4f*wall.get_sz_dim(!wdim))*rgen.rand_uniform(1.5, 1.8)); // width > height
									cube_t c;
									set_wall_width(c, (zval + rgen.rand_uniform(0.54, 0.62)*window_vspace), 0.5*height, 2); // Z
									set_wall_width(c, wall.get_center_dim(!wdim), 0.5*width, !wdim);
									c.d[wdim][ wdir] = wall_pos;
									c.d[wdim][!wdir] = wall_pos + (wdir ? -1.0 : 1.0)*0.04*wall_thickness; // shift into room
									objs.emplace_back(c, TYPE_PICTURE, room_id, wdim, !wdir, RO_FLAG_NOCOLL, light_amt); // picture faces dir opposite the wall
									set_obj_id(objs);
								}
							}
						}
					}
					else if (rtype == RTYPE_LIVING) {
						colorRGBA const color(get_couch_color(rgen));
						unsigned const place_start(objs.size());
						unsigned const num_couches(1 + rgen.rand_bool()); // 1-2

						for (unsigned n = 0; n < num_couches; ++n) { // place up to 2 couches
							if (place_model_along_wall(OBJ_MODEL_COUCH, TYPE_COUCH, sub_room, 0.40, rgen, zval, room_id, light_amt, sub_room, start_ix, 0.0, 4, 1, color)) {
								blockers.push_back(objs.back()); // add couch to blockers for table
								blockers.back().expand_by_xy(0.5*door_width); // add extra padding
							}
						}
						if (rgen.rand_bool() && place_model_along_wall(OBJ_MODEL_RCHAIR, TYPE_RCHAIR, sub_room, 0.5, rgen, zval, room_id, light_amt, sub_room, start_ix, 0.0)) {
							blockers.push_back(objs.back());
						}
						for (unsigned n = 0; n < 10; ++n) { // 10 attempts to place a valid table; wood, no chairs, tall=0
							if (add_table_and_chairs(rgen, sub_room, blockers, room_id, room_center, chair_color, 0.25, light_amt, 0, 0, 0)) {
								cube_t const table(objs.back());

								if (building_obj_model_loader.is_model_valid(OBJ_MODEL_LAMP) && rgen.rand_bool()) {
									// maybe add a lamp on the table; lamps aren't added to regular living rooms, but we can add them here
									float const height(0.25*window_vspace), radius(0.5*height*get_lamp_width_scale());

									if (radius < 3.0*min(table.dx(), table.dy())) { // add if it fits
										point center(gen_xy_pos_in_area(table, 1.2*radius, rgen, table.z2()));
										cube_t lamp(get_cube_height_radius(center, radius, height));
										unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_IN_MALL); // no collisions
										if (rgen.rand_bool()) {flags |= RO_FLAG_LIT;} // 50% chance of being lit
										objs.back().flags |= RO_FLAG_ADJ_TOP; // flag table as having something on it
										objs.emplace_back(lamp, TYPE_LAMP, room_id, rgen.rand_bool(), 0, flags, light_amt, SHAPE_CYLIN, lamp_colors[rgen.rand()%NUM_LAMP_COLORS]);
									}
								}
								break;
							}
							rgen.rand_mix(); // needed to get different rand values
						} // for n
						unsigned const num_plants(rgen.rand() % 3); // 0-2
						add_plants_to_room(rgen, sub_room, zval, room_id, light_amt, place_start, num_plants);
						// try to place a lamp 50% of the time; these go on the floor but aren't really meant to be floor lamps; they can be placeholders for future tall lamps
						//if (rgen.rand_bool()) {try_add_lamp(sub_room, window_vspace, room_id, RO_FLAG_IN_MALL, light_amt, blockers, objs, rgen);}
					}
					else if (rtype == RTYPE_KITCHEN || rtype == RTYPE_DINING) { // tables only, for now
						sub_room.expand_by_xy(room_pad); // make it larger to allow chairs to stick out a bit, since they don't block the player
						bool const use_tall_table(rgen.rand_float() < 0.25);
						int const wooden_or_plastic(rgen.rand() % 3); // randomly select between {wooden, plastic, wooden table with plastic chairs}
						add_table_and_chairs(rgen, sub_room, blockers, room_id, room_center, chair_color, 0.0, light_amt, 4, use_tall_table, wooden_or_plastic);
					}
					else if (rtype == RTYPE_OFFICE) {
						add_desk_to_room(rgen, sub_room, blockers, chair_color, zval, room_id, light_amt, start_ix, 1, 0, 1); // is_basement=1, desk_ix=0, no_computer=1
						add_bookcase_to_room(rgen, sub_room, zval, room_id, light_amt, start_ix, 0);
					}
					if (rtype == RTYPE_BED || rtype == RTYPE_LIVING || rtype == RTYPE_DINING) { // maybe add a rug
						if (rgen.rand_float() < 0.75) {add_rug_to_room(rgen, sub_room, zval, room_id, light_amt, start_ix);} // 75% of the time
					}
				} // end STORE_FURNITURE
				else if (store_type == STORE_APPLIANCE) {
					for (unsigned N = 0; N < 20; ++N) { // 20 attempts to place a new and valid model type
						float const val(rgen.rand_float());
						room_object obj_type(TYPE_NONE);
						float hscale(1.0);

						if (val < 0.4) { // 40% kitchen appliances
							if (val < 0.2) {obj_type = TYPE_FRIDGE; hscale = 0.75;}
							else           {obj_type = TYPE_STOVE ; hscale = 0.46;}
							// TODO: dishwasher (TYPE_KSINK with special case)?
						}
						else if (val < 0.6) { // 20% laundry appliances
							if (val < 0.5) {obj_type = TYPE_WASHER; hscale = 0.42;}
							else           {obj_type = TYPE_DRYER ; hscale = 0.38;}
						}
						else { // 40% bathroom plumbing fixtures
							if      (val < 0.70) {obj_type = TYPE_TUB   ; hscale = 0.20;}
							else if (val < 0.85) {obj_type = TYPE_TOILET; hscale = 0.35;}
							else                 {obj_type = TYPE_SINK  ; hscale = 0.45;}
							// TYPE_SHOWERTUB, WHEATER, TYPE_FURNACE?
						}
						// TYPE_CEIL_FAN?
						if (obj_type == TYPE_NONE) continue;
						assert(obj_type >= TYPE_TOILET);
						unsigned const model_type_id(obj_type - TYPE_TOILET);
						if (N < 10 && (types_used & (1 << model_type_id))) continue; // type already placed
						unsigned const model_id(model_type_id + OBJ_MODEL_TOILET);
						if (!building_obj_model_loader.is_model_valid(model_id)) continue; // no model
						vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
						float const height(hscale*window_vspace);
						bool const pri_dir(rgen.rand_bool());
						cube_t app;
						set_cube_zvals(app, zval, zval+height);
						bool any_placed(0);
					
						// place appliances around the perimeter of the sub-room facing outward
						for (unsigned D = 0; D < 2; ++D) {
							bool const sdim(bool(D) ^ pri_dir);
							bool const add_hoods(D == 0 && obj_type == TYPE_STOVE && building_obj_model_loader.is_model_valid(OBJ_MODEL_HOOD));
							float const width(height*sz.y/sz.z), depth(height*sz.x/sz.z), min_spacing(1.1*width + wall_thickness), sub_room_sz(r.get_sz_dim(sdim));
							unsigned const num(sub_room_sz/min_spacing); // take the floor
							if (num == 0) continue; // too large for this sub-room; shouldn't happen
							float const spacing(sub_room_sz/num);
							cube_t place_area(room_area);
							place_area.expand_in_dim(!sdim, -depth); // shink sides to avoid placing appliance facing a wall
							cube_t walls[2]; // for stove hoods

							for (unsigned n = 0; n < num; ++n) {
								set_wall_width(app, (r.d[sdim][0] + (n + 0.5)*spacing), 0.5*width, sdim);

								for (unsigned sdir = 0; sdir < 2; ++sdir) {
									float const edge_pos(r.d[!sdim][sdir]), dsign(sdir ? 1.0 : -1.0);
									app.d[!sdim][ sdir] = edge_pos; // front
									app.d[!sdim][!sdir] = edge_pos - dsign*depth; // back
									if (!place_area.contains_cube_xy(app))           continue; // against a wall or store window
									if (is_cube_close_to_doorway(app, r, 0.0, 1, 1)) continue; // blocking back hallway door
									if (door_blocker.intersects_xy(app))             continue; // blocking front door
									if (has_bcube_int(app, blockers))                continue; // blocked
									objs.emplace_back(app, obj_type, room_id, !sdim, sdir, RO_FLAG_IN_MALL, light_amt);
									blockers.push_back(app);
									any_placed = 1;

									if (add_hoods) { // add hood above the stove
										vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_HOOD)); // D, W, H
										float const height(width*sz.z/sz.y), depth(width*sz.x/sz.y); // scale to the width of the stove
										float const ceiling_z(zval + get_floor_ceil_gap()), z_top(ceiling_z + get_fc_thickness()); // shift up a bit because it's too low
										float const back_pos(app.d[!sdim][!sdir]);
										cube_t hood(app);
										set_cube_zvals(hood, z_top-height, z_top);
										hood.d[!sdim][sdir] = back_pos + dsign*depth; // extend front outward from wall
										objs.emplace_back(hood, TYPE_HOOD, room_id, !sdim, sdir, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, LT_GRAY);
										// add wall segment connecting this hood
										cube_t wall_seg(app);
										wall_seg.z2() = hood.z2();
										wall_seg.d[!sdim][ sdir] = back_pos; // front
										wall_seg.d[!sdim][!sdir] = back_pos - dsign*0.8*wall_thickness; // back
										walls[sdir].assign_or_union_with_cube(wall_seg);
									}
								} // for sdir
							} // for n
							for (unsigned d = 0; d < 2; ++d) {
								if (walls[d].is_all_zeros()) continue;
								objs.emplace_back(walls[d], TYPE_PG_WALL, room_id, !sdim, d, (RO_FLAG_IN_MALL | RO_FLAG_ADJ_TOP), light_amt); // draw top
								blockers.push_back(walls[d]);
							}
						} // for sdim
						if (any_placed) {
							types_used |= (1 << model_type_id);
							break; // done/success
						}
					} // for N
				} // end STORE_APPLIANCE
				else {assert(0);} // invalid store type
			} // for col
		} // for row
	}
	else if (store_type == STORE_FOOD) { // restaurant, coffee shop, etc.
		// TODO
	}
	else if (store_type == STORE_PETS) { // rats, snakes, birds, spiders, fish, etc.
		// add fish tanks along walls
		float const shelf_height(0.85*window_vspace), shelf_depth(0.25*window_vspace);
		add_shelves_along_walls(room_area, zval, room_id, light_amt, !dim, store_type, shelf_height, shelf_depth, 1, rgen); // place_inside=1
		// add door blocker to avoid placing a fishtank in front of the entry doorway
		cube_t door_blocker(doorway);
		door_blocker.expand_by_xy(0.25*window_vspace);
		objs.emplace_back(door_blocker, TYPE_BLOCKER, room_id, 0, 0, 0, light_amt, SHAPE_CUBE);
		// add a few more fishtanks if there's any extra space along front and back walls
		unsigned const num_fishtanks(rgen.rand() % 5); // 0-4
		cube_t place_area(room_area);
		place_area.d[dim][dir] -= (dir ? 1.0 : -1.0)*0.5*wall_thickness; // shink for front window frame clearance
		for (unsigned n = 0; n < num_fishtanks; ++n) {add_fishtank_to_room(rgen, room, zval, room_id, light_amt, objs_start, place_area);}
	}
	// add ducts and vents in the ceiling
	unsigned const num_vents(max(2U, (unsigned)round_fp(0.5*room_len/window_vspace))); // per side
	float const vent_spacing(room_len/num_vents), duct_height(0.2*window_vspace);
	unsigned const skip_dir((room_width < 0.8*room_len) ? rgen.rand_bool() : 2); // skip one side if room is narrow
	unsigned const duct_flags(RO_FLAG_IN_MALL | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI); // skip both ends
	cube_t duct(room);
	duct.z2() -= fc_thick;
	duct.z1()  = duct.z2() - duct_height;
	duct.expand_in_dim(dim, -wall_hthick); // shrink

	for (unsigned d = 0; d < 2; ++d) { // each side
		if (unsigned(d) == skip_dir) continue;
		float const dscale(d ? 1.0 : -1.0);
		duct.d[!dim][ d] = room.d[!dim][d] - dscale*wall_hthick;
		duct.d[!dim][!d] = duct.d[!dim][d] - dscale*0.3*window_vspace;
		objs.emplace_back(duct, TYPE_DUCT, room_id, dim, 0, (duct_flags | RO_FLAG_ADJ_TOP), light_amt, SHAPE_CUBE); // skip top as well; dir=0 (XY)
		// add vents along the duct
		cube_t duct_ext(duct); // duct => vent extension
		duct_ext.expand_in_z(-0.15*duct_height); // shrink
		duct_ext.d[!dim][ d] = duct.d[!dim][!d]; // flush with duct
		duct_ext.d[!dim][!d] = duct.d[!dim][!d] - dscale*0.1*window_vspace; // extend a bit further out

		for (unsigned n = 0; n < num_vents; ++n) {
			set_wall_width(duct_ext, (duct.d[dim][0] + (0.5 + n)*vent_spacing), 0.12*window_vspace, dim);
			objs.emplace_back(duct_ext, TYPE_DUCT, room_id, !dim, 0, duct_flags, light_amt, SHAPE_CUBE); // dir=0 (XY)
			cube_t vent(duct_ext);
			vent.d[!dim][d] = duct_ext.d[!dim][!d]; // shrink to end of duct
			vent.expand_by(0.005*window_vspace);
			objs.emplace_back(vent, TYPE_VENT, room_id, !dim, d, RO_FLAG_IN_MALL, light_amt, SHAPE_CUBE);
		} // for n
	} // for d
}

void building_t::add_row_of_bookcases(cube_t const &row, float zval, unsigned room_id, float light_amt, bool dim, bool place_inside) {
	float const window_vspace(get_window_vspace());
	float const bc_width(0.45*window_vspace), bc_depth(0.13*window_vspace), bc_height(0.75*window_vspace);
	float const wall_length(row.get_sz_dim(!dim));
	unsigned const num_bc(max(1, round_fp(wall_length/bc_width)));
	float const bcase_width(wall_length/num_bc);
	vect_room_object_t &objs(interior->room_geom->objs);

	for (unsigned d = 0; d < 2; ++d) { // for each side
		float const back_pos(place_inside ? row.d[dim][!d] : row.get_center_dim(dim));
		cube_t c;
		set_cube_zvals(c, zval, zval+bc_height);
		c.d[dim][!d] = back_pos;
		c.d[dim][ d] = back_pos + (d ? 1.0 : -1.0)*bc_depth;

		for (unsigned m = 0; m < num_bc; ++m) {
			c.d[!dim][0] = row.d[!dim][0] + m*bcase_width;
			c.d[!dim][1] = c  .d[!dim][0] + bcase_width;
			objs.emplace_back(c, TYPE_BCASE, room_id, dim, d, RO_FLAG_IN_MALL, light_amt);
			set_obj_id(objs);
		}
	} // for d
}

void building_t::add_shelves_along_walls(cube_t const &room_area, float zval, unsigned room_id, float light_amt, bool dim,
	unsigned store_type, float height, float depth, bool place_inside, rand_gen_t &rgen)
{
	cube_t c(room_area);
	set_cube_zvals(c, zval, zval+height);

	for (unsigned d = 0; d < 2; ++d) { // for each side
		float const back_pos(room_area.d[dim][bool(d) ^ place_inside]);
		c.d[dim][!d] = back_pos;
		c.d[dim][ d] = back_pos + (d ? 1.0 : -1.0)*depth;
		unsigned shelf_flags(RO_FLAG_INTERIOR | RO_FLAG_IN_MALL);
		if (!place_inside) {shelf_flags |= RO_FLAG_OPEN | RO_FLAG_NONEMPTY;} // shelf in middle of rooms; flag as open so that back is drawn; no empty shelves
		else if (rgen.rand_float() < 0.9) {shelf_flags |= RO_FLAG_NONEMPTY;} // shelf against room wall; 10% chance of being empty
		add_shelves(c, dim, !d, room_id, light_amt, shelf_flags, store_type, rgen);
	} // for d
}

int building_interior_t::get_store_id_for_room(unsigned room_id) const { // uses a slow linear iteration over stores, but not called much
	if (!get_room(room_id).is_store()) return -1; // not a store
	assert(has_mall());

	for (auto s = mall_info->stores.begin(); s != mall_info->stores.end(); ++s) {
		if (s->room_id == room_id) return (s - mall_info->stores.begin()); // found
	}
	return -1; // not found
}

