// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

using std::string;
string choose_store_name(rand_gen_t rgen);
colorRGBA choose_pot_color(rand_gen_t &rgen);

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
void building_t::get_mall_open_areas(cube_t const &room, vect_cube_t &openings) const {
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

void building_t::setup_mall_concourse(cube_t const &room, bool dim, bool dir, rand_gen_t &rgen) {
	assert(interior);
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), door_width(get_doorway_width());
	float const floor_thickness(get_floor_thickness()), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness()), trim_thickness(get_trim_thickness());
	
	// add wall section below the door on lower floors, since the entrace door is on the top level
	float const room_end(room.d[dim][!dir]); // entrance end
	door_t const &door(interior->get_ext_basement_door());
	assert(door.dim == dim);
	cube_t wall(door);
	set_cube_zvals(wall, room.z1()+fc_thick, door.z1()); // below the door
	wall.d[dim][!dir] = room_end; // doesn't exactly match the regular wall width, but not visible anyway
	wall.d[dim][ dir] = room_end + (dir ? 1.0 : -1.0)*wall_thickness; // extend into the room
	if (wall.dz() > 0.0) {interior->walls[dim].push_back(wall);}

	// handle upper floors
	if (interior->num_extb_floors == 1) return; // single floor; nothing else to do
	vect_cube_t openings;
	get_mall_open_areas(room, openings);
	assert(!openings.empty());
	reserve_extra(interior->stairwells, (openings.size() + 1));
	reserve_extra(interior->escalators, (openings.size() + 1));

	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const z(room.z1() + f*floor_spacing), zc(z - fc_thick), zf(z + fc_thick), floor_below_z(zf - floor_spacing), floor_above_z(zf + floor_spacing);
		
		// add side walkways
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			cube_t ww(room);
			ww.d[!dim][!d] = openings.front().d[!dim][d];
			interior->add_ceil_floor_pair(ww, zc, z, zf);
		}
		for (unsigned n = 0; n <= openings.size(); ++n) { // ends and gaps between openings
			bool const first(n == 0), last(n == openings.size());
			cube_t ww(openings.front());
			ww.d[dim][0] = (first ? room.d[dim][0] : openings[n-1].d[dim][1]);
			ww.d[dim][1] = (last  ? room.d[dim][1] : openings[n  ].d[dim][0]);
			interior->add_ceil_floor_pair(ww, zc, z, zf);
			// add stairs
			bool run_dir(first ? 1 : (last ? 0 : rgen.rand_bool())), side_dir(rgen.rand_bool());
			float ww_edge(ww.d[dim][run_dir]), ww_side(ww.d[!dim][side_dir]);
			float const stairs_len(1.5*floor_spacing), stairs_width(0.75*window_vspace);
			cube_t stairs_bc;
			set_cube_zvals(stairs_bc, floor_below_z, zf); // top of floor below to top of current floor
			stairs_bc.d[ dim][ !run_dir] = ww_edge; // at walkway
			stairs_bc.d[ dim][  run_dir] = ww_edge + (run_dir  ?  1.0 : -1.0)*stairs_len; // extend away from walkway
			stairs_bc.d[!dim][ side_dir] = ww_side + (side_dir ? -1.0 :  1.0)*wall_thickness;
			stairs_bc.d[!dim][!side_dir] = ww_side + (side_dir ? -1.0 :  1.0)*stairs_width;
			landing_t landing(stairs_bc, 0, 0, dim, !run_dir, 1, SHAPE_STRAIGHT, 0, 1, 0, 0, 1); // add_railing=1, roof_access=0, is_at_top=1, stack_conn=0, for_ramp=0, in_extb=1
			landing.num_stairs = round_fp(NUM_STAIRS_PER_FLOOR*floor_spacing/window_vspace);
			landing.in_mall    = 1;
			stairs_landing_base_t stairwell(landing);
			landing  .z1() = zc;
			stairwell.z2() = floor_above_z;
			interior->landings.push_back(landing);
			interior->stairwells.emplace_back(stairwell, 1); // num_floors=1
			// add escalators
			if (!first && !last) {run_dir ^= 1;} // opposite end except at concourse ends
			side_dir ^= 1; // opposite side
			ww_edge = ww.d[ dim][run_dir ];
			ww_side = ww.d[!dim][side_dir];
			float const delta_z(floor_spacing + 2.0*trim_thickness), e_width(0.8*door_width); // make slightly taller so that top surface is above the railing bottom trim
			cube_t epair;
			set_cube_zvals(epair, floor_below_z, floor_above_z);
			epair.d[ dim][ !run_dir] = ww_edge - (run_dir  ?  1.0 : -1.0)*door_width; // overlapping walkway
			epair.d[ dim][  run_dir] = ww_edge + (run_dir  ?  1.0 : -1.0)*(1.0*delta_z + door_width); // extend away from walkway, 45 degree slope
			epair.d[!dim][ side_dir] = ww_side + (side_dir ? -1.0 :  1.0)*wall_thickness;
			epair.d[!dim][!side_dir] = ww_side + (side_dir ? -1.0 :  1.0)*2.2*e_width; // set width of the pair
			
			for (unsigned s = 0; s < 2; ++s) { // place two side-by-side escalators with opposite directions
				// extend 90% of floor thickness below; enough to hide building people animated feet, but not enough to clip through the ceiling below
				escalator_t e(epair, dim, !run_dir, (bool(s) ^ dim ^ run_dir ^ 1), door_width, delta_z, 0.9*floor_thickness, 1); // in_mall=1
				e.d[!dim][!s] = epair.d[!dim][s] + (s ? -1.0 : 1.0)*e_width;
				interior->escalators.push_back(e);
			}
		} // for n
	} // for f
	if (!openings.empty()) { // add elevator
		unsigned const opening_ix(choose_one_center(openings.size(), rgen));
		cube_t const opening(openings[opening_ix]);
		bool const edir((openings.size() & 1) ? rgen.rand_bool() : (opening_ix == openings.size()/2)); // closer to center; random if tied
		float const ww_edge(opening.d[dim][!edir]);
		// Note: elevator extends half a floor width below and above the room; is this okay, or can it clip through other objects?
		elevator_t elevator(room, interior->ext_basement_hallway_room_id, dim, !edir, 1, 1, 1); // at_edge=1, interior_room=1, in_mall=1
		elevator.d[dim][!edir] = ww_edge; // adjacent to walkway
		elevator.d[dim][ edir] = ww_edge + (edir ? 1.0 : -1.0)*1.6*door_width; // extend away from walkway by depth
		set_wall_width(elevator, opening.get_center_dim(!dim), 0.8*door_width, !dim); // set width

		if (is_in_city) {
			// extend elevator up to street level if there's space?
		}
		interior->elevators.push_back(elevator);
	}
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
	interior->store_doorways.push_back(doorway);
	
	// add window on each side of the doorway
	for (unsigned side = 0; side < 2; ++side) {
		cube_t window(walls_cut);
		window.d[dim][!side] = opening.d[dim][side];
		set_wall_width(window, wall_pos, 0.25*wall_thickness, !dim);
		assert(window.is_strictly_normalized());
		interior->int_windows.emplace_back(window, room_ix);
	}
}

void building_t::add_mall_stores(cube_t const &room, bool dim, bool entrance_dir, rand_gen_t &rgen) {
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), wall_thickness(get_wall_thickness());
	float const min_width(4.0*window_vspace), max_width(9.0*window_vspace), min_depth(6.0*window_vspace), max_depth(8.0*window_vspace);
	float const bathroom_width(4.0*window_vspace);
	unsigned const num_floors(interior->num_extb_floors);
	bool added_bathrooms(0);
	vect_cube_t &side_walls(interior->walls[!dim]);
	interior->store_bounds_by_floor.resize(num_floors);

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
		float depths[2] = {};
		cube_t store, floor_bcube(room);
		store.z1() = floor_bcube.z1() = room .z1() + f*floor_spacing;
		store.z2() = floor_bcube.z2() = store.z1() +   floor_spacing;

		// place stores on each side of concourse
		for (unsigned d = 0; d < 2; ++d) { // sides of mall
			float const wall_pos(room.d[!dim][d]), depth(rgen.rand_uniform(min_depth, max_depth)); // consistent depth for each side + floor so that walls can be shared
			float pos(room.d[dim][0]), pos_end(room.d[dim][1]);
			float const middle(0.5*(pos + pos_end));
			// prevent exterior wall of store from clipping through parking garage wall
			if (entrance_dir) {pos_end -= wall_thickness;} else {pos += wall_thickness;}
			depths[d] = depth;
			store.d[!dim][!d] = wall_pos;
			store.d[!dim][ d] = wall_pos + (d ? 1.0 : -1.0)*depth;
			bool has_adj_store(0);
			
			while (pos + min_width < pos_end) { // continue until we can't fit a min width room
				float next_pos(pos);
				bool is_bathroom(0);

				if (!added_bathrooms && pos < middle && pos+bathroom_width > middle) {
					next_pos += bathroom_width;
					is_bathroom = 1;
				}
				else {
					next_pos += rgen.rand_uniform(min_width, max_width);
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
						cube_t bathroom(bathrooms);
						bathroom.d[dim][!e] = sep_pos;
						room_t Room(bathroom, basement_part_ix);
						Room.assign_all_to((bool(e) ^ wm_first) ? RTYPE_WOMENS : RTYPE_MENS);
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
			float const wall_pos(room.d[dim][d]), depth(rgen.rand_uniform(min_depth, max_depth));
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
		interior->store_bounds_by_floor[f] = floor_bcube;
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
				bool const is_bathroom(is_bathroom(r->get_room_type(0)));
				//if (is_bathroom) continue; // no door to bathrooms?
				cube_t conn_room(*r);
				set_cube_zvals(conn_room, hall.z1(), hall.z2());
				cube_t const door_cut(add_ext_basement_door(conn_room, doorway_width, !dim, d, 1, (is_tall_room && !is_bathroom), rgen)); // is_end_room=1
				subtract_cube_from_cubes(door_cut, interior->walls[!dim], nullptr, 1); // no holes; clip_in_z=1
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
	colorRGBA color;
	plant_loc_t(point const &p, float r, colorRGBA const &c) : sphere_t(p, r), color(c) {}
};

// this is for the central mall concourse; store objects are added in add_mall_store_objs() below; treated as a single floor
unsigned building_t::add_mall_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, vect_cube_t &rooms_to_light) {
	bool const mall_dim(interior->extb_wall_dim);
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), fc_thick(get_fc_thickness()), doorway_width(get_doorway_width());
	float const wall_thickness(get_wall_thickness()), trim_thick(get_trim_thickness()), room_centerline(room.get_center_dim(!mall_dim));
	float const light_amt = 1.0; // fully lit, for now
	unsigned const num_floors(interior->num_extb_floors);
	cube_t const mall_center(get_mall_center(room));
	vect_cube_t openings, railing_cuts, railing_segs, temp, pillars;
	vector<plant_loc_t> plant_locs;
	point plant_loc(0.0, 0.0, zval);
	get_mall_open_areas(room, openings);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned num_elevators(0);
	colorRGBA pot_colors[3]; // {stairs, escalators, pillars}
	for (unsigned n = 0; n < 3; ++n) {pot_colors[n] = choose_pot_color(rgen);}

	// gather railing cuts and plant locs
	for (stairwell_t const &s : interior->stairwells) {
		if (!s.in_mall) continue;
		railing_cuts.push_back(s);
		railing_cuts.back().expand_in_dim( s.dim, doorway_width); // add padding on both ends
		railing_cuts.back().expand_in_dim(!s.dim, 0.8*wall_thickness); // add space for railings
		room.has_stairs = 255; // should this be set earlier?
		// place plants to the sides of stairs
		float const plant_radius(0.15*window_vspace);
		plant_loc[s.dim] = 0.9*s.d[s.dim][!s.dir] + 0.1*s.d[s.dim][s.dir]; // bottom end of stairs

		for (unsigned d = 0; d < 2; ++d) {
			plant_loc[!s.dim] = s.d[!s.dim][d] + (d ? 1.0 : -1.0)*1.25*plant_radius;
			plant_locs.emplace_back(plant_loc, plant_radius, pot_colors[0]);
		}
	} // for s
	for (escalator_t const &e : interior->escalators) {
		if (!e.in_mall) continue;
		railing_cuts.push_back(e);
		railing_cuts.back().expand_in_dim(!e.dim, 0.5*wall_thickness); // add space for railings
		// escalators are placed in pairs, so we add one plant to each side
		float const plant_radius(0.16*window_vspace);
		bool const side(e.dir ^ e.move_dir);
		plant_loc[ e.dim] = 0.9*e.d[e.dim][!e.dir] + 0.1*e.d[e.dim][e.dir]; // bottom end of escalator
		plant_loc[!e.dim] = e.d[!e.dim][side] + (side ? 1.0 : -1.0)*1.25*plant_radius;
		plant_locs.emplace_back(plant_loc, plant_radius, pot_colors[1]);

		// add to upper walkway(s) if not too close to railing
		if (plant_loc[!mall_dim] - plant_radius > mall_center.d[!mall_dim][0] && plant_loc[!mall_dim] + plant_radius < mall_center.d[!mall_dim][1]) {
			plant_loc[e.dim] = 0.1*e.d[e.dim][!e.dir] + 0.9*e.d[e.dim][e.dir]; // top end of escalator

			for (unsigned f = 1; f < num_floors; ++f) { // skip first floor
				plant_loc.z += floor_spacing;
				plant_locs.emplace_back(plant_loc, plant_radius, pot_colors[1]);
			}
			plant_loc.z = zval; // restore ground floor zval
		}
	} // for e
	for (elevator_t const &e : interior->elevators) {
		if (!e.in_mall) continue;
		railing_cuts.push_back(e.get_bcube_padded(doorway_width));
		++num_elevators;
		// place clock on back of elevator and front if there's space
		bool const digital(rgen.rand_bool());
		add_clock_to_cube(e, (zval + floor_spacing) , room_id, light_amt, e.dim, !e.dir, digital); // back

		if (floor_spacing >= 1.5*window_vspace) {
			cube_t place_cube(e);
			place_cube.d[e.dim][e.dir] += (e.dir ? 1.0 : -1.0)*0.5*wall_thickness; // account for front blocker
			add_clock_to_cube(place_cube, (zval + floor_spacing + 0.5*window_vspace) , room_id, light_amt, e.dim, e.dir, digital); // front
		}
	} // for e
	room.has_elevator = num_elevators; // should this be set earlier?
	// add vertical support pillars
	float const pillar_hwidth(2.0*wall_thickness);
	cube_t pillar(room); // copy room zvals

	for (cube_t const &opening : openings) {
		for (unsigned dir = 0; dir < 2; ++dir) { // each side of opening
			float const pillar_pos(opening.d[!mall_dim][dir] + (dir ? -1.0 : 1.0)*0.7*pillar_hwidth);
			set_wall_width(pillar, pillar_pos, pillar_hwidth, !mall_dim);
			set_wall_width(pillar, opening.get_center_dim(mall_dim), pillar_hwidth, mall_dim);
			pillars     .push_back(pillar);
			railing_cuts.push_back(pillar);
			// place plants to the sides of pillars
			float const plant_radius(0.2*window_vspace);
			plant_loc[!mall_dim] = pillar_pos;

			for (unsigned d = 0; d < 2; ++d) {
				plant_loc[mall_dim] = pillar.d[mall_dim][d] + (d ? 1.0 : -1.0)*1.25*plant_radius;
				plant_locs.emplace_back(plant_loc, plant_radius, pot_colors[2]);
			}
		} // for dir
	} // for opening
	// add walkway railings
	for (unsigned f = 1; f < num_floors; ++f) { // skip first floor
		float const z(room.z1() + f*floor_spacing), plate_thickness(0.1*wall_thickness), vbar_hwidth(0.35*wall_thickness);
		float const zb1(z - fc_thick - trim_thick), zb2(z + fc_thick + trim_thick), zt1(z + 0.42*window_vspace), zt2(zt1 + 0.04*window_vspace);
		colorRGBA const bar_color(LT_GRAY);
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
					bot_bar.expand_by_xy (0.3*wall_thickness); // additional expand from plate
					bot_bar.expand_in_dim(!dim, plate_thickness); // fill the corner notch
					set_cube_zvals(bot_bar, zb1, zb2);
					objs.emplace_back(bot_bar, TYPE_METAL_BAR, room_id, !dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, bar_color, skip_mask);
					subtract_cubes_from_cube(railing, railing_cuts, railing_segs, temp, 2); // check zval overlap

					for (cube_t const &r : railing_segs) {
						// add top bar
						cube_t bar(r);
						bar.expand_by_xy (0.4*wall_thickness); // additional expand from plate
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
	} // for d

	if (is_in_city) {
		// add skylight(s)?
	}
	// add a fountain in the center of an opening
	if (!openings.empty() && building_obj_model_loader.is_model_valid(OBJ_MODEL_FOUNTAIN)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FLAG)); // W, D, H
		cube_t const opening(openings[choose_one_center(openings.size(), rgen)]);
		float const max_radius(0.25*min(opening.dx(), opening.dy()));
		float height(0.9*floor_spacing*rgen.rand_uniform(0.9, 1.1)), radius(0.5*height*0.5*(sz.x + sz.y)/sz.z); // use average of width and depth for radius
		if (radius > max_radius) {height *= max_radius/radius; radius = max_radius;} // reduce size if radius is too large
		cube_t fbc;
		set_cube_zvals(fbc, zval, zval+height);
		for (unsigned d = 0; d < 2; ++d) {set_wall_width(fbc, opening.get_center_dim(d), radius, d);}
		objs.emplace_back(fbc, TYPE_BLDG_FOUNT, room_id, rgen.rand_bool(), rgen.rand_bool(), 0, light_amt, SHAPE_CYLIN); // random dim/dir
		objs.back().item_flags = rgen.rand(); // select a random sub_model_id
	}
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
		objs.emplace_back(plant, TYPE_PLANT, room_id, 0, 0, RO_FLAG_ADJ_BOT, light_amt, SHAPE_CYLIN, p.color);
		set_obj_id(objs);
	}

	// TODO: palm trees, TYPE_TABLE, TYPE_CHAIR, TYPE_PICTURE, TYPE_TCAN, TYPE_SIGN, TYPE_RDESK, TYPE_DUCT, TYPE_VASE, TYPE_BENCH, TYPE_TV, TYPE_BAR_STOOL
	//cube_t walk_area(room);
	//walk_area.expand_by_xy(-wall_thickness);

	// add pillars last so that we can check lights against them; must be added last
	unsigned const pillars_start(objs.size());
	float fe_height(0.0), fe_radius(0.0);
	bool const add_fire_extinguishers(get_fire_ext_height_and_radius(window_vspace, fe_height, fe_radius));
	unsigned pillar_ix(0);
	bool add_fe(rgen.rand_bool()); // random start on even vs. odd pillars

	for (cube_t const &pillar : pillars) {
		objs.emplace_back(pillar, TYPE_OFF_PILLAR, room_id, !mall_dim, 0, 0, light_amt, SHAPE_CUBE, WHITE, EF_Z12);
		
		if (add_fire_extinguishers) { // maybe add fire extinguisher on pillar
			if (add_fe) {
				bool const dim(!mall_dim), dir(pillar.get_center_dim(dim) < room_centerline);
				add_fire_ext(fe_height, fe_radius, zval, pillar.d[dim][!dir], pillar.get_center_dim(!dim), room_id, light_amt, dim, dir);
			}
			if (!(pillar_ix++ & 1)) {add_fe ^= 1;} // swap every pair of pillars, alternating sides
		}
	} // for pillar
	return pillars_start;
}

void building_t::add_mall_store_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id) {
	// TODO: TYPE_SHELFRACK, TYPE_CASHREG, etc.
}

