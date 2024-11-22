// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

using std::string;
string choose_store_name(rand_gen_t rgen);
colorRGBA choose_pot_color(rand_gen_t &rgen);
bool is_invalid_city_placement_for_cube(cube_t const &c);
void add_city_plot_cut(cube_t const &cut);

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

unsigned choose_one_center(unsigned num, rand_gen_t &rgen) {
	unsigned ix(num/2); // center value
	if (!(num & 1) && rgen.rand_bool()) {--ix;} // tie breaker if even
	return ix;
}

void building_t::add_mall_se_landing(cube_t const &c, bool is_escalator, bool se_dim, bool se_dir, bool ww_dir) {
	assert(c.is_strictly_normalized());
	interior->add_ceil_floor_pair(c, c.z1(), c.zc(), c.z2());
	interior->mall_landings.emplace_back(c, (8*is_escalator + 4*se_dim + 2*se_dir + ww_dir));
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
	if (!openings.empty()) { // add elevator
		unsigned const opening_ix(choose_one_center(num_openings, rgen));
		cube_t const opening(openings[opening_ix]);
		bool const edir((num_openings & 1) ? rgen.rand_bool() : (opening_ix == num_openings/2)); // closer to center; random if tied
		float const ww_edge(opening.d[dim][!edir]), width(1.6*door_width), depth(width);
		// Note: elevator extends half a floor width below and above the room; is this okay, or can it clip through other objects?
		elevator_t elevator(room, interior->ext_basement_hallway_room_id, dim, !edir, 1, 1, 1); // at_edge=1, interior_room=1, in_mall=1
		elevator.d[dim][!edir] = ww_edge; // front is adjacent to walkway
		elevator.d[dim][ edir] = ww_edge + (edir ? 1.0 : -1.0)*depth; // extend back away from walkway by depth
		set_wall_width(elevator, opening.get_center_dim(!dim), 0.5*width, !dim); // set width

		if (0 && is_in_city) { // extend elevator up to street level if there's space?
			cube_t test_cube(elevator);
			set_cube_zvals(test_cube, elevator.z2(), elevator.z2()+floor_spacing);
			test_cube.d[dim][!edir] -= (edir ? 1.0 : -1.0)*depth; // extend by depth in front of elevator

			if (!is_cube_city_placement_invalid(test_cube)) {
				// TODO: need exterior building surrounding elevator, and a way for the player to interact with it on street level
				// TODO: needs to be a city blocker for pedestrians, etc.
				// TODO: need to cut a hole in the plot for the elevator shaft
				// TODO: should extend window_vspace above, but elevator won't have a stop at this pos, unless parking garage is two floors
				elevator.z2() += floor_spacing;
				add_city_plot_cut(elevator);
			}
		}
		interior->elevators.push_back(elevator);
	}
	if (is_in_city) { // add skylight(s)?
		// must call is_cube_city_placement_invalid() and add_city_plot_cut()
	}
}

bool building_t::is_cube_city_placement_invalid(cube_t const &c) const { // for mall skylights, elevator, etc.
	if (!is_basement_room_not_int_bldg(c, nullptr, 1)) return 1; // check for buildings above; no exclude, allow_outside_grid=1
	return is_invalid_city_placement_for_cube(c);
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

void building_t::add_mall_store(cube_t const &store, cube_t const &window_area, bool dim, bool dir, bool &has_adj_store) {
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
	interior->store_doorways.emplace_back(doorway, room_ix);
	
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
	float const bathroom_width(4.0*window_vspace);
	unsigned const num_floors(interior->num_extb_floors);
	bool added_bathrooms(0);
	vect_cube_t &side_walls(interior->walls[!dim]);
	interior->mall_store_bounds = room; // start with the mall concourse
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
		unsigned const rooms_start(interior->rooms.size()), store_walls_start(side_walls.size());
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

				if (!added_bathrooms && pos < middle && pos+store_width > middle) { // crosses the middle; make a bathroom
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
					float const fc_thick(get_fc_thickness()), doorway_width(get_doorway_width());
					float const sep_pos(0.5*(pos + next_pos)); // separator between men's and women's rooms
					cube_t door, bathrooms(store);
					set_cube_zvals(door, store.z1()+fc_thick, store.z1()+window_vspace-fc_thick);
					door.d[!dim][0] = door.d[!dim][1] = wall_pos;
					if (use_low_ceiling) {bathrooms.z2() = bathrooms.z1() + window_vspace;} // set normal ceiling height
					interior->mall_bathrooms = bathrooms;
					
					for (unsigned e = 0; e < 2; ++e) {
						room_type const rtype((bool(e) ^ wm_first) ? RTYPE_WOMENS : RTYPE_MENS);
						cube_t bathroom(bathrooms);
						bathroom.d[dim][!e] = sep_pos;
						room_t Room(bathroom, basement_part_ix);
						Room.assign_all_to(rtype);
						Room.interior        = 2; // mark as extended basement
						Room.is_single_floor = 1;
						Room.is_office       = 1; // required for creating office building bathroom
						interior->rooms.push_back(Room);
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
					added_bathrooms = 1;
				}
				else { // regular store
					add_mall_store(store, store, dim, d, has_adj_store); // window_area is full storefront
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
			add_mall_store(store, window_area, !dim, d, has_adj_store);
			floor_bcube.union_with_cube(store);
		} // for d
		interior->mall_store_bounds.union_with_cube(floor_bcube);
		// add a narrower non-public hallway behind each row of stores
		bool const is_tall_room(floor_spacing > 1.1*window_vspace);
		unsigned const rooms_end(interior->rooms.size());
		float const doorway_width(get_doorway_width()), hall_width(2.0*doorway_width);
		cube_t hall(floor_bcube); // extends full length of the mall, including the end stores
		hall.z2() = floor_bcube.z1() + window_vspace; // normal height

		for (unsigned d = 0; d < 2; ++d) { // sides of mall
			float const wall_pos(floor_bcube.d[!dim][d]);
			if (room.d[!dim][d] == wall_pos) continue; // no stores added to this side, skip; excludes bathrooms
			hall.d[!dim][0]  = hall.d[!dim][1] = wall_pos; // move to outside of rooms
			hall.d[!dim][d] += (d ? 1.0 : -1.0)*hall_width;
			if (is_store_placement_invalid(hall)) continue;
			room_t Room(hall, basement_part_ix, 0, 1); // num_lights=0, is_hallway=1
			Room.interior = 2; // mark as extended basement
			interior->rooms.push_back(Room);
			add_extb_room_floor_and_ceil(hall);
			unsigned const hall_walls_start(side_walls.size());
			bool has_adj_store(0); // unused
			add_mall_room_walls(hall, wall_thickness, !dim, d, 1, has_adj_store, interior->walls); // at_mall_end=1 (to cover missing stores)

			// split walls between stores and hallways into two half slices so that they can be textured differently; is there a way to skip drawing interior faces?
			for (auto w = side_walls.begin()+store_walls_start; w != side_walls.end(); ++w) {
				if (w->d[!dim][0] >= wall_pos || w->d[!dim][1] <= wall_pos) continue; // wrong wall
				w->d[!dim][bool(d) ^ (w < side_walls.begin()+hall_walls_start) ^ 1] = wall_pos; // set to half width
			}
			// connect to rooms with doors
			for (auto r = interior->rooms.begin()+rooms_start; r != interior->rooms.begin()+rooms_end; ++r) {
				if (!r->intersects(hall)) continue; // wrong side
				room_type const rtype(r->get_room_type(0));
				bool const is_bath(is_bathroom(rtype));
				//if (is_bath) continue; // no door to bathrooms?
				cube_t conn_room(*r);
				set_cube_zvals(conn_room, hall.z1(), hall.z2());
				cube_t const door_cut(add_ext_basement_door(conn_room, doorway_width, !dim, d, 1, (is_tall_room && !is_bath), rgen)); // is_end_room=1
				subtract_cube_from_cubes(door_cut, interior->walls[!dim], nullptr, 1); // no holes; clip_in_z=1
				interior->doors.back().rtype = rtype; // flag door (in particular if bathroom) so that the correct sign is added
			}
		} // for d
	} // for f
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

	if (num_steps > 0) {
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
	}
}

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
	cube_t entrance_stairs_bcube(interior->get_ext_basement_door()); // approximate
	entrance_stairs_bcube.expand_by_xy(2.0*(floor_spacing - window_vspace) + doorway_width); // add plenty of clearance around the door for the stairs
	if (num_floors == 1) {blockers.push_back(entrance_stairs_bcube);} // entrance is on the first/ground floor
	cube_t pillar(room); // copy room zvals

	// gather railing cuts and plant locs
	for (stairwell_t const &s : interior->stairwells) {
		if (!s.in_mall) continue;
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
		if (!e.in_mall) continue;
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
	for (cube_with_ix_t const &c : interior->mall_landings) {
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
	for (cube_t const &d : interior->store_doorways) {
		bool const dim(d.dy() < d.dx());
		// add ceiling box where the gate would come down from
		cube_t cbox(d);
		cbox.z1() = d.z2() - 0.1*window_vspace;
		cbox.expand_in_dim( dim,  1.0*wall_thickness); // grow
		cbox.expand_in_dim(!dim, -0.5*wall_thickness); // shrink
		objs.emplace_back(cbox, TYPE_METAL_BAR, room_id, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, GRAY);
		// add closed gate or other objects?
	} // for d

	// add a fountain in the center of an opening
	int fountain_opening_ix(-1), fc_opening_ix(-1);

	if (!openings.empty() && building_obj_model_loader.is_model_valid(OBJ_MODEL_FOUNTAIN)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FLAG)); // W, D, H
		fountain_opening_ix = choose_one_center(num_openings, rgen);
		cube_t const &opening(openings[fountain_opening_ix]);
		float const max_radius(0.25*min(opening.dx(), opening.dy()));
		float height(0.9*floor_spacing*rgen.rand_uniform(0.9, 1.1)), radius(0.5*height*0.5*(sz.x + sz.y)/sz.z); // use average of width and depth for radius
		if (radius > max_radius) {height *= max_radius/radius; radius = max_radius;} // reduce size if radius is too large
		cube_t fbc;
		set_cube_zvals(fbc, zval, zval+height);
		for (unsigned d = 0; d < 2; ++d) {set_wall_width(fbc, opening.get_center_dim(d), radius, d);}
		objs.emplace_back(fbc, TYPE_BLDG_FOUNT, room_id, rgen.rand_bool(), rgen.rand_bool(), 0, light_amt, SHAPE_CYLIN); // random dim/dir
		objs.back().item_flags = rgen.rand(); // select a random sub_model_id
		blockers.push_back(fbc);
	}
	// trashcan setup for trashcans along walls and food court pillars
	bool const is_cylin_tcan(rgen.rand_bool());
	float const tcan_height((is_cylin_tcan ? 0.28 : 0.4)*window_vspace), tcan_radius(0.12*window_vspace);
	room_obj_shape const tcan_shape(is_cylin_tcan ? SHAPE_CYLIN : SHAPE_CUBE);
	colorRGBA const tcan_color(is_cylin_tcan ? LT_GRAY : colorRGBA(0.8, 0.6, 0.4)); // tan for cube trashcans

	// add a food court to one of the openings
	if (num_openings >= ((fountain_opening_ix >= 0) ? 2 : 1)) {
		while (1) {
			fc_opening_ix = rgen.rand() % num_openings;
			if (fc_opening_ix != fountain_opening_ix) break;
		}
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
			add_mall_trashcan_contents(rgen, tcan_obj, room_id, light_amt);
			blockers.push_back(tcan);
		} // for p
		add_food_court_objs(rgen, place_area, zval, room_id, light_amt, blockers);
	}
	// add objects to remaining openings
	for (unsigned i = 0; i < num_openings; ++i) {
		if ((int)i == fountain_opening_ix || (int)i == fc_opening_ix) continue; // already occupied
		// TODO: set of TYPE_VASE
	} // for i
	// if there are bathrooms, add a water fountain between them
	if (!interior->mall_bathrooms.is_all_zeros() && building_obj_model_loader.is_model_valid(OBJ_MODEL_WFOUNTAIN)) {
		bool const wf_dir(interior->mall_bathrooms.get_center_dim(!mall_dim) > room_centerline);
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_WFOUNTAIN)); // D, W, H
		float const height(0.25*window_vspace), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z);
		float const z1(zval + 0.18*window_vspace), wall_pos(room.d[!mall_dim][wf_dir] + (wf_dir ? -1.0 : 1.0)*0.5*wall_thickness);
		cube_t wf;
		set_wall_width(wf, interior->mall_bathrooms.get_center_dim(mall_dim), hwidth, mall_dim);
		set_cube_zvals(wf, z1, z1+height);
		wf.d[!mall_dim][ wf_dir] = wall_pos;
		wf.d[!mall_dim][!wf_dir] = wall_pos - (wf_dir ? 1.0 : -1.0)*depth;
		objs.emplace_back(wf, TYPE_WFOUNTAIN, room_id, !mall_dim, wf_dir, 0, light_amt, SHAPE_CUBE);
	}
	// add potted plants
	for (plant_loc_t const &p : plant_locs) {
		cube_t const plant(get_cube_height_radius(p.pos, p.radius, 4.0*p.radius));
		unsigned const flags(RO_FLAG_ADJ_BOT | (p.upper_floor ? RO_FLAG_HAS_EXTRA : 0)); // flag upper floor plants so that the bottom sides of leaves are drawn
		objs.emplace_back(plant, TYPE_PLANT, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN, p.color);
		set_obj_id(objs);
	}
	// place objects along mall-store walls on each floor, but not in front of glass
	cube_t place_area(room);
	place_area.expand_by_xy(-wall_thickness);
	vect_cube_t wall_blockers;
	unsigned const num_benches  (num_openings/2 + 1); // per floor per side
	unsigned const num_trashcans(num_openings/4 + 1); // per floor per side
	float const clearance(get_min_front_clearance_inc_people());
	float const bench_height(0.42*window_vspace), bench_depth(0.23*window_vspace), bench_hlen(0.5*window_vspace*rgen.rand_uniform(0.6, 0.8));
	float const tcan_end_pad(2.0*tcan_radius), bench_end_pad(2.0*bench_hlen);
	unsigned const NUM_MALL_BENCH_COLORS = 6;
	colorRGBA const mall_bench_colors[NUM_MALL_BENCH_COLORS] = {WHITE, LT_BROWN, DK_GRAY, LT_BROWN, colorRGBA(0.1, 0.4, 0.8), LT_BROWN}; // LT_BROWN becomes wood texture
	colorRGBA const &bench_color(mall_bench_colors[rgen.rand() % NUM_MALL_BENCH_COLORS]);

	for (unsigned f = 0; f < num_floors; ++f) {
		float const z(zval + f*floor_spacing);
		cube_t tcan, bench;
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
					if (entrance_stairs_bcube.intersects(test_cube))    continue; // too close to entrance stairs
					if (interior->mall_bathrooms.intersects(test_cube)) continue; // too close to bathroom doors and water fountain
					if (has_bcube_int(test_cube, interior->int_windows) || has_bcube_int(test_cube, interior->store_doorways)) continue; // too close to store door or window
					if (has_bcube_int(test_cube, wall_blockers))        continue;
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
					if (entrance_stairs_bcube.intersects(test_cube))    continue; // too close to entrance stairs
					if (interior->mall_bathrooms.intersects(test_cube)) continue; // too close to bathroom doors and water fountain
					if (has_bcube_int(test_cube, interior->int_windows) || has_bcube_int(test_cube, interior->store_doorways)) continue; // too close to store door or window
					if (has_bcube_int(test_cube, wall_blockers))        continue;
					wall_blockers.push_back(tcan);
					wall_blockers.back().expand_by_xy(window_vspace); // don't place two nearby trashcans
					room_object_t const tcan_obj(tcan, TYPE_TCAN, room_id, !mall_dim, !d, RO_FLAG_IN_MALL, light_amt, tcan_shape, tcan_color);
					objs.push_back(tcan_obj);
					add_mall_trashcan_contents(rgen, tcan_obj, room_id, light_amt);
					break;
				} // for N
			} // for n
			wall_blockers.clear();
		} // for d
	} // for f
	// add a reception desk in front of end walls that have no store
	for (unsigned f = 0; f < num_floors; ++f) {
		for (unsigned d = 0; d < 2; ++d) { // ends
			if (f+1 == num_floors && d == interior->extb_wall_dir) continue; // likely blocked by entrance stairs
			//float const z(zval + f*floor_spacing);
		}
	} // for f
	// TODO: palm trees, TYPE_PICTURE?, TYPE_BAR_STOOL?

	// add pillars last so that we can check lights against them
	unsigned const pillars_start(objs.size());
	for (cube_t const &pillar : pillars) {objs.emplace_back(pillar, TYPE_OFF_PILLAR, room_id, !mall_dim, 0, 0, light_amt, SHAPE_CUBE, WHITE, EF_Z12);}
	return pillars_start;
}

void building_t::add_mall_trashcan_contents(rand_gen_t &rgen, room_object_t const &tcan, unsigned room_id, float tot_light_amt) { // with trash
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
	room_object_t table_obj(table, TYPE_TABLE, room_id, 0, 0, RO_FLAG_IN_MALL, tot_light_amt, SHAPE_CUBE, WHITE);
	table_obj.item_flags = tid_tag; // sets table texture
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
			objs.emplace_back(chair, TYPE_CHAIR, room_id, !dim, dir, RO_FLAG_IN_MALL, tot_light_amt, SHAPE_CUBE, chair_color);
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
	unsigned const tid_tag(rgen.rand()); // sets table texture
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
			add_mall_table_with_chairs(rgen, table, place_area, chair_color, room_id, tot_light_amt, dim, tid_tag, fc_blockers);
		} // while
	} // for r
	// TODO: trashcans
	return 1;
}

void building_t::add_mall_store_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id) {
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), doorway_width(get_doorway_width()), wall_thickness(get_wall_thickness());
	float const light_amt = 1.0; // fully lit, for now
	// get doorway bcube
	cube_t doorway;

	for (cube_with_ix_t const &d : interior->store_doorways) {
		if (d.ix != room_id) continue;
		assert(doorway.is_all_zeros()); // must be exactly one
		doorway = d;
	}
	assert(!doorway.is_all_zeros()); // must be found
	bool const dim(doorway.dy() < doorway.dx()), dir(room.get_center_dim(dim) < doorway.get_center_dim(dim));
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());

	// add store name on sign above the entrance
	// TODO: TYPE_SIGN

	// add checkout counter(s)/cash register(s) to the side of the door
	if (rgen.rand_float() < 0.67) {
		bool const side(rgen.rand_bool());
		cube_t checkout_area(room);
		checkout_area.d[!dim][!side] = doorway.d[!dim][side]; // to the side of the doorway
		checkout_area.expand_in_dim(!dim, -1.0*doorway_width); // add padding
		checkout_area.d[dim][!dir] = room.get_center_dim(dim); // only add to the front half of the store
		add_checkout_objs(checkout_area, zval, room_id, light_amt, objs_start, dim, side, (side ^ dim ^ 1));
	}
	// TODO: TYPE_SHELFRACK, TYPE_DUCT, etc.
}

