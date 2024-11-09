// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"

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
		unsigned opening_ix(openings.size()/2); // center opening
		if (!(openings.size() & 1) && rgen.rand_bool()) {--opening_ix;} // tie breaker if even
		cube_t const opening(openings[opening_ix]);
		bool const edir(rgen.rand_bool());
		float const ww_edge(opening.d[dim][!edir]);
		// Note: elevator extends half a floor width below and above the room; is this okay, or can it clip through other objects?
		elevator_t elevator(room, interior->ext_basement_hallway_room_id, dim, !edir, 1, 1, 1); // at_edge=1, interior_room=1, in_mall=1
		elevator.d[dim][!edir] = ww_edge; // adjacent to walkway
		elevator.d[dim][ edir] = ww_edge + (edir ? 1.0 : -1.0)*1.6*door_width; // extend away from walkway by depth
		set_wall_width(elevator, opening.get_center_dim(!dim), 0.8*door_width, !dim); // set width
		interior->elevators.push_back(elevator);
	}
}

void building_t::add_mall_stores(cube_t const &room, bool dim, rand_gen_t &rgen) {
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace());
	float const wall_thickness(get_wall_thickness()), fc_thick(get_fc_thickness());
	float const wall_half_thick(0.5*wall_thickness), ceil_gap(get_floor_thick_val()*floor_spacing - fc_thick);
	float const min_width(2.0*window_vspace), max_width(8.0*window_vspace), min_depth(4.0*window_vspace), max_depth(6.0*window_vspace);

	for (unsigned f = 0; f < interior->num_extb_floors; ++f) {
		cube_t store;
		store.z1() = room .z1() + f*floor_spacing;
		store.z2() = store.z1() +   floor_spacing;

		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			float const pos_end(room.d[dim][1] - wall_thickness), wall_pos(room.d[!dim][d]);
			float pos(room.d[dim][0] + wall_thickness);
			store.d[!dim][!d] = wall_pos;
			
			while (pos < pos_end) {
				float next_pos(pos + rgen.rand_uniform(min_width, max_width));
				if (next_pos + min_width > pos_end) {next_pos = pos_end;} // clamp to far end of mall, and prevent a narrow store
				store.d[ dim][0] = pos;
				store.d[ dim][1] = next_pos;
				store.d[!dim][d] = wall_pos + (d ? 1.0 : -1.0)*rgen.rand_uniform(min_depth, max_depth); // set depth
				pos = next_pos;

				if (!is_basement_room_under_mesh_not_int_bldg(store)) { // invalid, try min depth
					store.d[!dim][d] = wall_pos + (d ? 1.0 : -1.0)*min_depth;
					if (!is_basement_room_under_mesh_not_int_bldg(store)) continue; // still invalid, skip this store
				}
				assert(store.is_strictly_normalized());
				unsigned const room_ix(interior->rooms.size());
				room_t Room(store, basement_part_ix);
				Room.assign_all_to(RTYPE_STORE);
				Room.interior        = 2; // mark as extended basement
				Room.is_single_floor = 1;
				interior->rooms.push_back(Room);
				cube_t ceiling(store), floor(store);
				ceiling.z1() = store.z2() - fc_thick;
				floor  .z2() = store.z1() + fc_thick;
				interior->ceilings.push_back(ceiling);
				interior->floors  .push_back(floor);
				interior->basement_ext_bcube.assign_or_union_with_cube(store);

				// add walls
				for (unsigned wdim = 0; wdim < 2; ++wdim) {
					for (unsigned wdir = 0; wdir < 2; ++wdir) {
						if (wdim == !dim && wdir != d) continue; // already have walls on this side
						float half_thick(wall_half_thick);
						cube_t wall(store);
						set_wall_width(wall, store.d[wdim][wdir], half_thick, wdim);
						interior->walls[wdim].push_back(wall);
					}
				} // for dim
				// add door/window openings
				cube_t walls_cut(store);
				walls_cut.z1() += fc_thick;
				walls_cut.z2() -= ceil_gap;
				walls_cut.expand_in_dim(dim, -(0.25*window_vspace + 0.15*store.get_sz_dim(dim))); // shrink
				set_wall_width(walls_cut, wall_pos, 2.0*wall_thickness, !dim);
				subtract_cube_from_cubes(walls_cut, interior->walls[!dim], nullptr, 1); // no holes; clip_in_z=1
				// add window
				// TODO: to cut a doorway
				cube_t window(walls_cut);
				set_wall_width(window, wall_pos, 0.25*wall_thickness, !dim);
				interior->int_windows.emplace_back(window, room_ix);
			} // while
		} // for d
	} // for n
}

void building_t::add_mall_stairs() { // connecting to the entrance door
	if (!has_mall()) return;
	room_t const &room(interior->get_extb_start_room());
	door_t const &door(interior->get_ext_basement_door());
	bool const dim(interior->extb_wall_dim), dir(interior->extb_wall_dir);
	assert(door.dim == dim);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const room_id(interior->ext_basement_hallway_room_id);
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	// add stairs under the door if needed
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
			objs.emplace_back(railing, TYPE_RAILING, room_id, !dim, !d, RO_FLAG_HAS_EXTRA, 1.0, SHAPE_CUBE, GOLD); // no balusters
		}
	}
}

void building_t::add_mall_lower_floor_lights(room_t const &room, unsigned room_id, unsigned lights_start, light_ix_assign_t &light_ix_assign) {
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness());
	vect_cube_t openings;
	get_mall_open_areas(room, openings);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_end(objs.size());
	assert(lights_start <= objs_end);

	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const zc(room.z1() + f*floor_spacing - fc_thick); // bottom of ceiling

		for (unsigned i = lights_start; i < objs_end; ++i) {
			room_object_t const &obj(objs[i]);
			if (obj.type != TYPE_LIGHT)       continue; // should this ever fail?
			if (has_bcube_int(obj, openings)) continue; // skip lights over the openings
			room_object_t light(obj);
			light.translate_dim(2, (zc - obj.z2()));
			light.obj_id = light_ix_assign.get_ix_for_light(light);
			objs.push_back(light);
		} // for i
	} // for f
}

unsigned building_t::add_mall_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t &rooms_to_light) {
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), fc_thick(get_fc_thickness()), doorway_width(get_doorway_width());
	float const wall_thickness(get_wall_thickness()), trim_thick(get_trim_thickness());
	bool const mall_dim(interior->extb_wall_dim);
	vect_cube_t openings, railing_cuts, railing_segs, temp, pillars;
	get_mall_open_areas(room, openings);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned num_elevators(0);

	// gather railing cuts
	for (stairwell_t const &s : interior->stairwells) {
		if (!s.in_mall) continue;
		railing_cuts.push_back(s);
		railing_cuts.back().expand_in_dim( s.dim, doorway_width); // add padding on both ends
		railing_cuts.back().expand_in_dim(!s.dim, 0.8*wall_thickness); // add space for railings
		room.has_stairs = 255; // should this be set earlier?
	}
	for (escalator_t const &e : interior->escalators) {
		if (!e.in_mall) continue;
		railing_cuts.push_back(e);
		railing_cuts.back().expand_in_dim(!e.dim, 0.5*wall_thickness); // add space for railings
	}
	for (elevator_t const &e : interior->elevators) {
		if (!e.in_mall) continue;
		railing_cuts.push_back(e.get_bcube_padded(doorway_width));
		++num_elevators;
	}
	room.has_elevator = num_elevators; // should this be set earlier?
	// add vertical support pillars
	float const pillar_hwidth(2.0*wall_thickness);
	cube_t pillar(room); // copy room zvals

	for (cube_t const &opening : openings) {
		for (unsigned dir = 0; dir < 2; ++dir) { // each side of opening
			set_wall_width(pillar, (opening.d[!mall_dim][dir] + (dir ? -1.0 : 1.0)*0.7*pillar_hwidth), pillar_hwidth, !mall_dim);
			set_wall_width(pillar, opening.get_center_dim(mall_dim), pillar_hwidth, mall_dim);
			pillars.push_back(pillar);
			railing_cuts.push_back(pillar);
		}
	} // for opening
	// add walkway railings
	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
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
					objs.emplace_back(bot_bar, TYPE_METAL_BAR, room_id, !dim, 0, 0, 1.0, SHAPE_CUBE, bar_color, skip_mask);
					subtract_cubes_from_cube(railing, railing_cuts, railing_segs, temp, 2); // check zval overlap

					for (cube_t const &r : railing_segs) {
						// add top bar
						cube_t bar(r);
						bar.expand_by_xy (0.4*wall_thickness); // additional expand from plate
						bar.expand_in_dim(!dim, plate_thickness); // fill the corner notch
						set_cube_zvals(bar, zt1, zt2);
						objs.emplace_back(bar, TYPE_METAL_BAR, room_id, !dim, 0, 0, 1.0, SHAPE_CUBE, bar_color); // can't skip drawing of ends if clipped by stairs
						// add glass
						cube_t panel(r);
						set_cube_zvals(panel, zb2, zt1);
						objs.emplace_back(panel, TYPE_INT_WINDOW, room_id, dim, 0, 0, 1.0);

						// add vertical bars where railing was clipped
						for (unsigned d = 0; d < 2; ++d) {
							if (r.d[!dim][d] == railing.d[!dim][d]) continue; // not clipped
							set_wall_width(vbar, r.d[!dim][d], vbar_hwidth, !dim);
							set_wall_width(vbar, centerline,   vbar_hwidth,  dim);
							objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, 0, 1.0, SHAPE_CUBE, bar_color, EF_Z12); // skip top and bottom
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
							objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, 0, 1.0, SHAPE_CUBE, bar_color, EF_Z12); // skip top and bottom
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
				objs.emplace_back(vbar, TYPE_METAL_BAR, room_id, 0, 0, 0, 1.0, SHAPE_CUBE, bar_color, EF_Z12); // skip top and bottom
			}
		}
	} // for f
	
	// TODO: potted plants, palm trees, TYPE_TABLE, TYPE_CHAIR, TYPE_PICTURE, TYPE_WBOARD, TYPE_TCAN, TYPE_SIGN, TYPE_PLANT, TYPE_RDESK,
	// TYPE_VENT, TYPE_DUCT, TYPE_VASE, TYPE_BENCH, TYPE_CLOCK, TYPE_TV, TYPE_FIRE_EXT, TYPE_BAR_STOOL, TYPE_WFOUNTAIN, TYPE_BLDG_FOUNT
	//cube_t walk_area(room);
	//walk_area.expand_by_xy(-wall_thickness);

	// add pillars last so that we can check lights against them
	unsigned const pillars_start(objs.size());
	for (cube_t const &pillar : pillars) {objs.emplace_back(pillar, TYPE_OFF_PILLAR, room_id, 0, 0, 0, 1.0, SHAPE_CUBE, WHITE, EF_Z12);}
	return pillars_start;
}

void building_t::add_mall_store_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id) {
	// TODO
}

