// 3D World - Building Bathrooms
// by Frank Gennari 12/13/2025

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t


extern object_model_loader_t building_obj_model_loader;


bool building_t::is_room_office_bathroom(room_t &room, float zval, unsigned floor) const { // Note: may also update room flags
	if (!room.is_office || !is_bathroom(room.get_room_type(floor))) return 0;
	if (!room_has_stairs_or_elevator(room, zval, floor))            return 1;
	room.clear_room_type(floor); // not a bathroom; can't call assign_to() because it skips bathrooms
	return 0;
}

bool building_t::add_bathroom_objs(rand_gen_t rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt,
	unsigned lights_start, unsigned objs_start, unsigned floor, bool is_basement, bool add_shower_tub, unsigned &added_bathroom_objs_mask)
{
	// Note: zval passed by reference due to add_flooring()
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());

	if (!is_house && !skylights.empty()) { // check for skylights if it's an office building; houses can have bathroom skylights (my old house did)
		// should we allow bathroom stalls (with ceilings?) even if there's a skylight?
		cube_t test_cube(room);
		set_cube_zvals(test_cube, zval, zval+floor_spacing);
		if (check_skylight_intersection(test_cube)) return 0;
	}
	cube_t const room_bounds(get_walkable_room_bounds(room));
	cube_t place_area(room_bounds);
	place_area.expand_by(-0.5*wall_thickness);
	vector2d const place_area_sz(place_area.dx(), place_area.dy());
	float const min_room_dim(min(place_area_sz.x, place_area_sz.y));
	if (min_room_dim < 0.7*floor_spacing) return 0; // room is too small (should be rare)
	bool const have_toilet(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOILET)), have_sink(building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK));
	vect_room_object_t &objs(interior->room_geom->objs);

	if ((have_toilet || have_sink) && is_cube() && !is_industrial() && !is_restroom()) {
		// bathroom with at least a toilet or sink; cube shaped parts only; no industrial or park restroom; add flooring
		int const flooring_type(is_residential() ? (is_basement ? (int)FLOORING_CONCRETE : (int)FLOORING_TILE) : (int)FLOORING_MARBLE);
		if (flooring_type == FLOORING_CONCRETE && get_material().basement_floor_tex.tid == get_concrete_tid()) {} // already concrete
		else { // replace carpet/wood with marble/tile/concrete
			zval = add_flooring(room, zval, room_id, tot_light_amt, flooring_type); // move the effective floor up
		}
	}
	if (have_toilet && (room.is_office || (!is_basement && (is_prison() || is_restaurant() || is_restroom())))) {
		// office, above ground prison bathroom, restaurant, and park restroom have stalls
		if (min_room_dim > (has_house_floorplan() ? 1.4 : 1.5)*floor_spacing && max(place_area_sz.x, place_area_sz.y) > 2.0*floor_spacing) {
			if (divide_bathroom_into_stalls(rgen, room, zval, room_id, tot_light_amt, floor, lights_start, objs_start)) { // large enough, divide into bathroom stalls
				added_bathroom_objs_mask |= (PLACED_TOILET | PLACED_SINK);
				return 1;
			}
		}
	}
	float const tub_height_factor(0.2), sink_height_factor(0.45); // in units of floor spacing
	unsigned const vanity_obj_ix(objs.size());
	bool placed_obj(0), placed_toilet(0);
	bool added_vanity(is_house && !is_basement && rgen.rand_float() < 0.75 && add_vanity_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start)); // maybe add vanity

	// place toilet first because it's in the corner out of the way and higher priority
	if (have_toilet) { // have a toilet model
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
		float const height(0.35*floor_spacing), width(height*sz.y/sz.z), length(height*sz.x/sz.z); // for toilet
		unsigned const first_corner(rgen.rand() & 3);
		bool const first_dim(rgen.rand_bool());

		for (unsigned pass = 0; pass < (added_vanity ? 2U : 1U) && !placed_toilet; ++pass) { // {without, with} vanity
			if (pass == 1) { // vanity pass; remove the vanity; not a bathroom without a toilet
				objs.resize(vanity_obj_ix);
				room.clear_has_mirror(); // no more mirror
				added_vanity = 0;
			}
			for (unsigned n = 0; n < 4 && !placed_toilet; ++n) { // try 4 room corners
				unsigned const corner_ix((first_corner + n)&3);
				bool const xdir(corner_ix&1), ydir(corner_ix>>1);
				point const corner(place_area.d[0][xdir], place_area.d[1][ydir], zval);
				if (!check_pt_within_part_sides(corner)) continue; // invalid corner

				for (unsigned d = 0; d < 2 && !placed_toilet; ++d) { // try both dims
					bool const dim(bool(d) ^ first_dim), dir(dim ? ydir : xdir);
					cube_t c(corner, corner);
					c.d[0][!xdir] += (xdir ? -1.0 : 1.0)*(dim ? width : length);
					c.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(dim ? length : width);
					for (unsigned e = 0; e < 2; ++e) {c.d[!dim][e] += ((dim ? xdir : ydir) ? -1.5 : 1.5)*wall_thickness;} // extra padding on left and right sides
					c.z2() += height;
					cube_t c2(c); // used for placement tests
					c2.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.8*length; // extra padding in front of toilet, to avoid placing other objects there (sink and tub)
					c2.expand_in_dim(!dim, 0.4*width); // more padding on the sides
					if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, room, 0.0, 1)) continue; // bad placement
					objs.emplace_back(c,  TYPE_TOILET,  room_id, dim, !dir, 0, tot_light_amt);
					add_bathroom_plumbing(objs.back());
					objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker cube to ensure no other object overlaps this space
					add_poi_dim_dir(c, room_id, dim, !dir, 0.5); // dscale=0.5
					placed_obj = placed_toilet = 1; // done
					added_bathroom_objs_mask  |= PLACED_TOILET;

					// try to place a roll of toilet paper on the adjacent wall
					bool const tp_dir(dim ? xdir : ydir);
					float const length(0.18*height), wall_pos(c.get_center_dim(dim)), far_edge_pos(wall_pos + (dir ? -1.0 : 1.0)*0.5*length);
					cube_t const &part(get_part_for_room(room));

					// if this wall has windows and bathroom has multiple exterior walls (which means it has non-glass block windows), don't place a TP roll
					if (is_basement || !has_int_windows() || classify_room_wall(room, zval, !dim, tp_dir, 0) != ROOM_WALL_EXT ||
						!is_val_inside_window(part, dim, far_edge_pos, get_hspacing_for_part(part, dim), get_window_h_border()) || count_ext_walls_for_room(room, zval) <= 1)
					{
						add_tp_roll(room_bounds, room_id, tot_light_amt, !dim, tp_dir, length, (c.z1() + 0.7*height), wall_pos);
					}
				} // for d
			} // for n
		} // for pass
		if (!placed_toilet) { // if the toilet can't be placed in a corner, allow it to be placed anywhere; needed for small offices
			placed_toilet = place_model_along_wall(OBJ_MODEL_TOILET, TYPE_TOILET, room, 0.35, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8);
			placed_obj   |= placed_toilet;
			added_bathroom_objs_mask |= PLACED_TOILET;

			if (placed_toilet) { // if toilet was placed, try to place a roll of toilet paper on the same wall as the toilet
				room_object_t const &toilet(objs.back()); // okay if this is the blocker
				add_poi_dim_dir(toilet, room_id, toilet.dim, toilet.dir, 0.5); // dscale=0.5
				
				// Note: not calling is_val_inside_window() here because I don't have a test case for that and it may not even be possible to get here when the toilet is next to a window
				if (is_basement || !has_int_windows() || classify_room_wall(room, zval, toilet.dim, !toilet.dir, 0) != ROOM_WALL_EXT) { // check for possible windows
					bool place_dir(rgen.rand_bool()); // pick a random starting side

					for (unsigned d = 0; d < 2; ++d) {
						float const length(0.18*height), wall_pos(toilet.d[!toilet.dim][place_dir] + (place_dir ? 1.0 : -1.0)*0.5*width);
						if (add_tp_roll(room_bounds, room_id, tot_light_amt, toilet.dim, !toilet.dir, length, (toilet.z1() + 0.7*height), wall_pos, 1)) break; // check_valid=1
						place_dir ^= 1; // try the other dir
					} // for d
				}
			}
		}
	} // end have_toilet
	unsigned const pre_shower_tub_sink_ix(objs.size());
	float tub_clearance(0.0);
	bool no_tub(0);

	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_TUB) && building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK)) {
		// add clearance around the tub if a sink placed on the opposite wall blocks the player
		vector3d const tub_sz (building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TUB )); // D, W, H
		vector3d const sink_sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SINK)); // D, W, H
		float const tub_depth (tub_height_factor *floor_spacing*tub_sz .x/tub_sz .z);
		float const sink_depth(sink_height_factor*floor_spacing*sink_sz.x/sink_sz.z);
		float const clearance(get_min_front_clearance_inc_people());
		if (min_room_dim < (clearance + tub_depth + sink_depth)) {tub_clearance = clearance;}
	}
	for (unsigned n = 0; n < 20; ++n) { // 20 tries to add a tub/shower and sink
		unsigned bathroom_objs_mask(0);

		// try to add a shower; 50% chance if on first floor of a house; not in basements (due to drawing artifacts)
		if (add_shower_tub && (is_apt_or_hotel() || floor > 0 || rgen.rand_bool())) {
			float shower_dx(0.0), shower_dy(0.0), wall_thick(0.0);
			bool hdim(0); // handle dim/long dim
			// use a shower + tub combo for basements due to glass drawing artifacts, and for arpartments and hotels with small bathrooms; requires tub model
			bool const place_shower_tub((is_basement || is_apt_or_hotel() || rgen.rand_float() < 0.25) && building_obj_model_loader.is_model_valid(OBJ_MODEL_TUB));
			float const shower_height(place_shower_tub ? get_floor_ceil_gap() : 0.8*floor_spacing), tub_height(tub_height_factor*floor_spacing);

			if (place_shower_tub) {
				vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TUB)); // D, W, H
				hdim       = rgen.rand_bool();
				shower_dx  = tub_height*(sz.y/sz.z); // width
				shower_dy  = tub_height*(sz.x/sz.z); // depth
				wall_thick = 0.05*shower_dx;
				shower_dx += wall_thick;
				if (hdim) {swap(shower_dx, shower_dy);}
			}
			else {
				for (unsigned d = 0; d < 2; ++d) {(d ? shower_dy : shower_dx) = rgen.rand_uniform(0.4, 0.5)*floor_spacing;}
				hdim = (shower_dx < shower_dy); // larger dim, must match handle/door drawing code
			}
			unsigned const first_corner(rgen.rand() & 3);
			bool placed_shower(0), is_ext_wall[2][2] = {0};
		
			if (!is_basement && has_int_windows()) { // precompute which walls are exterior, {dim}x{dir}; basement walls are not considered exterior because there are no windows
				for (unsigned d = 0; d < 4; ++d) {is_ext_wall[d>>1][d&1] = (classify_room_wall(room, zval, (d>>1), (d&1), 0) == ROOM_WALL_EXT);}
			}
			for (unsigned ar = 0; ar < 2; ++ar) { // try both aspect ratios/door sides
				for (unsigned n = 0; n < 4; ++n) { // try 4 room corners
					unsigned const corner_ix((first_corner + n)&3);
					bool const xdir(corner_ix&1), ydir(corner_ix>>1), dirs[2] = {xdir, ydir};
					point const corner(room_bounds.d[0][xdir], room_bounds.d[1][ydir], zval); // flush against the wall
					cube_t c(corner, corner);
					c.d[0][!xdir] += (xdir ? -1.0 : 1.0)*shower_dx;
					c.d[1][!ydir] += (ydir ? -1.0 : 1.0)*shower_dy;
					c.z2() += shower_height; // set height

					if (!place_shower_tub || has_int_windows()) {
						// exterior walls aren't drawn in the correct order for shower glass alpha blend, so skip any exterior walls;
						// also don't place shower tubs near windows; assume they're long enough to always overlap a window when placed on an exterior wall
						bool is_bad(0);

						for (unsigned d = 0; d < 2; ++d) { // check for window intersection
							if (is_ext_wall[!d][dirs[!d]]) {is_bad = 1; break;}
						}
						if (is_bad) continue;
					}
					cube_t c2(c); // used for placement tests

					if (!place_shower_tub) { // shower: extend out by door width on the side that opens, and a small amount on the other side
						c2.d[0][!xdir] += (xdir ? -1.0 : 1.0)*((!hdim) ? 1.1*shower_dy : 0.2*shower_dx);
						c2.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(  hdim  ? 1.1*shower_dx : 0.2*shower_dy);
					}
					if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, room, 0.0, 1)) continue; // bad placement
					bool const dim(place_shower_tub ? !hdim : xdir), dir(place_shower_tub ? !(hdim ? xdir : ydir) : ydir); // different encodings
					objs.emplace_back(c, (place_shower_tub ? TYPE_SHOWERTUB : TYPE_SHOWER), room_id, dim, dir, 0, tot_light_amt);
					set_obj_id(objs); // selects tile texture/color

					if (place_shower_tub) { // add the tub part as well
						bool const wall_dir(hdim ? ydir : xdir);
						// set flag to indicate which side is the wall for adding the shower head, and make open by default
						unsigned const lo_hi_flag(wall_dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
						objs.back().flags |= lo_hi_flag | RO_FLAG_OPEN;
						objs.back().color  = (is_basement ? WHITE : wall_color); // color of the end wall
						cube_t tub(c);
						tub.z2() = c.z1() + tub_height;
						tub.d[!dim][!wall_dir] -= (wall_dir ? -1.0 : 1.0)*wall_thick; // shrink off the wall
						tub.translate_dim(dim, (dir ? 1.0 : -1.0)*0.1*get_trim_thickness()); // shift away from the wall slightly to prevent Z-fighting
						objs.emplace_back(tub, TYPE_TUB, room_id, dim, dir, lo_hi_flag, tot_light_amt);
						bathroom_objs_mask |= PLACED_TUB;
						no_tub = 1;
					}
					objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker cube to ensure no other object overlaps this space
					placed_obj = placed_shower = 1;
					bathroom_objs_mask  |= PLACED_SHOWER;
					break; // done
				} // for n
				if (placed_shower) break; // done
				swap(shower_dx, shower_dy); // try the other aspect ratio
				hdim ^= 1;
			} // for ar
		}
		// place a tub, but not in office buildings; placed before the sink because it's the largest and the most limited in valid locations;
		// skip for the last 5 iterations so that we can add a sink instead
		if (!no_tub && add_shower_tub && n < 15 && (!is_basement || rgen.rand_bool())) { // 50% of the time if in the basement
			cube_t place_area_tub(room_bounds);
			place_area_tub.expand_by(-get_trim_thickness()); // just enough to prevent z-fighting and intersecting the wall trim
			unsigned const tub_obj_ix(objs.size());
		
			if (place_model_along_wall(OBJ_MODEL_TUB, TYPE_TUB, room, tub_height_factor, rgen, zval, room_id, tot_light_amt, place_area_tub, objs_start, 0.4)) {
				room_object_t const &tub(objs[tub_obj_ix]);
				placed_obj = 1;
				bathroom_objs_mask |= PLACED_TUB;

				if (tub_clearance > 0.0 && objs.back().type == TYPE_BLOCKER) {
					cube_t min_blocker(tub);
					min_blocker.expand_by_xy(tub_clearance);
					min_blocker.intersect_with_cube_xy(place_area_tub);
					objs.back().union_with_cube(min_blocker);
				}
				if (rgen.rand_bool()) { // add a bar of soap on the edge of the tub
					float const one_inch(get_one_inch()), soap_hlen(2.0*one_inch), soap_hwidth(1.25*one_inch), soap_height(1.0*one_inch); // 4x2.5x1
					colorRGBA const soap_color(soap_colors[rgen.rand() % NUM_SOAP_COLORS]);
					cube_t soap;
					set_cube_zvals(soap, tub.z2(), tub.z2()+soap_height);
					set_wall_width(soap, (tub.d[tub.dim][!tub.dir] + (tub.dir ? 1.0 : -1.0)*1.75*soap_hwidth), soap_hwidth, tub.dim);
					set_wall_width(soap, (rgen.rand_uniform((tub.d[!tub.dim][0] + 2.0*soap_hlen), (tub.d[!tub.dim][1] - 2.0*soap_hlen))), soap_hlen, !tub.dim);
					objs.emplace_back(soap, TYPE_BAR_SOAP, room_id, !tub.dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_ROUNDED_CUBE, soap_color);
				}
			}
		}
		unsigned const sink_obj_ix(objs.size());

		if (added_vanity) { // added vanity, no need to place a sink
			bathroom_objs_mask |= PLACED_SINK; // vanity includes a sink
		}
		else if (place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, sink_height_factor, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6)) {
			placed_obj = 1;
			bathroom_objs_mask |= PLACED_SINK;
			assert(sink_obj_ix < objs.size());
			room_object_t const &sink(objs[sink_obj_ix]); // sink, not blocker
			add_poi_dim_dir(sink, room_id, sink.dim, sink.dir, 0.5); // dscale=0.5
		
			if (point_in_water_area(sink.get_cube_center())) {} // no medicine cabinet, because the reflection system doesn't support both a mirror and water reflection
			else if (is_parking() && !is_basement) {} // no mirror in parking structure since reflections don't work there
			else if (is_basement || classify_room_wall(room, zval, sink.dim, !sink.dir, 0) != ROOM_WALL_EXT) { // interior wall only
				// add a mirror/medicine cabinet above the sink
				float const mirror_expand(0.1*sink.get_sz_dim(!sink.dim));
				cube_t mirror(sink); // start with the sink left and right position
				mirror.expand_in_dim(!sink.dim, mirror_expand); // make slightly wider
				set_cube_zvals(mirror, sink.z2(), sink.z2()+0.3*floor_spacing);
				mirror.d[sink.dim][!sink.dir] = room_bounds.d[sink.dim][!sink.dir];
				mirror.d[sink.dim][ sink.dir] = mirror.d[sink.dim][!sink.dir] + (sink.dir ? 1.0 : -1.0)*1.0*wall_thickness; // thickness

				if (!overlaps_other_room_obj(mirror, objs_start, 0, &sink_obj_ix)) { // check_all=0; skip sink + blocker
					if (interior->is_cube_close_to_doorway(mirror, room, 0.0, 1, 1)) { // check doorways as well, since mirror is wider thank sink
						mirror.expand_in_dim(!sink.dim, -mirror_expand); // undo expand and try for a narrow mirror
					}
					// this mirror is actually 3D, so we enable collision detection; treat as a house even if it's in an office building
					unsigned flags(RO_FLAG_IS_HOUSE); // Note: not necessarily a house
					if (count_ext_walls_for_room(room, mirror.z1()) == 1) {flags |= RO_FLAG_INTERIOR;} // flag as interior if windows are opaque glass blocks
					objs.emplace_back(mirror, TYPE_MED_CAB, room_id, sink.dim, sink.dir, flags, tot_light_amt); // Note: invalidates sink reference
					set_obj_id(objs); // for crack texture selection/orient
					room.set_has_mirror();
				}
			}
		}
		else if (n+1 < 20) { // not the last try, rip up and retry
			objs.resize(pre_shower_tub_sink_ix);
			continue;
		}
		added_bathroom_objs_mask |= bathroom_objs_mask;
		break; // success or last try
	} // for n
	if (is_house) {maybe_add_radiator_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);}
	else if (floor <= NUM_RTYPE_SLOTS && room.get_room_type(floor) == RTYPE_MENS && building_obj_model_loader.is_model_valid(OBJ_MODEL_URINAL)) {
		// assigned to men's room; add a urinal at the wall, not against an exterior wall
		float const urinal_zval(zval + 0.1*floor_spacing); // shift it up
		cube_t bounds(room_bounds);
		bounds.expand_by(-0.01*wall_thickness); // shrink slightly to prevent Z-figthing
		place_model_along_wall(OBJ_MODEL_URINAL, TYPE_URINAL, room, 0.4, rgen, urinal_zval, room_id, tot_light_amt, bounds, objs_start, 2.0, 4, 0, WHITE, 0, 0, 0, 0, 0.0, 1);
	}
	if ((room.is_office || is_conv_store()) && !room.is_nested()) {add_door_sign("Restroom", room, zval, room_id);} // add office/store bathroom sign; not for hospital rooms
	return placed_obj;
}

bool building_t::add_vanity_to_room(rand_gen_t &rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t const room_bounds(get_walkable_room_bounds(room));
	vector2d const room_sz(room_bounds.dx(), room_bounds.dy());
	vect_room_object_t &objs(interior->room_geom->objs);
	// select an interior wall that doesn't have a door
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors; there should be no exterior doors
	vector<pair<float, unsigned>> avail_walls;
	unsigned num_ext_walls(0);

	for (unsigned d = 0; d < 4; ++d) {
		bool const dim(d >> 1), dir(d & 1);
		if (classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) {++num_ext_walls; continue;} // skip exterior wall
		cube_t wall(room);
		set_wall_width(wall, room.d[dim][dir], wall_thickness, dim);
		if (has_bcube_int(wall, doorways)) continue;
		avail_walls.emplace_back(rgen.rand_uniform(1.0, 1.1)*room_sz[dim], d); // prefer long wall, with some randomness
	} // for d
	sort(avail_walls.begin(), avail_walls.end());

	for (auto const &cand : avail_walls) {
		bool const dim(cand.second >> 1), dir(cand.second & 1);
		float const wall_len(room_sz[!dim]), wall_edge(room_bounds.d[dim][dir]); // no trim padding
		float const height(0.34*floor_spacing), depth(0.74*height), length(rgen.rand_uniform(0.8, 1.0)*floor_spacing), dsign(dir ? -1.0 : 1.0);
		unsigned flags(is_house ? RO_FLAG_IS_HOUSE : 0);
		cube_t vanity(room_bounds);
		vanity.expand_in_dim(!dim, -get_trim_thickness()); // avoid Z-fighting with exterior wall
		set_cube_zvals(vanity, zval, zval+height);
		vanity.d[dim][!dir] = wall_edge + dsign*depth;

		if (wall_len <= length) { // short wall: vanity is full width
			flags |= (RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP); // flag to not draw either end
		}
		else if (wall_len < 1.5*length) { // middle size: vanity against 2 walls/corner
			bool const shift_dir(rgen.rand_bool());
			vanity.d[!dim][shift_dir] += (shift_dir ? -1.0 : 1.0)*(wall_len - length);
			flags |= (shift_dir ? RO_FLAG_ADJ_BOT : RO_FLAG_ADJ_TOP); // flag to not draw one end
		}
		else { // long wall: vanity centered
			set_wall_width(vanity, (room_bounds.get_center_dim(!dim) + 0.1*rgen.signed_rand_float()*length), 0.5*length, !dim);
		}
		cube_t blocker(vanity);
		blocker.d[dim][!dir] += 0.9*dsign*depth;
		blocker.expand_in_dim(!dim, 0.05*depth); // include overhang on sides
		if (overlaps_obj_or_placement_blocked(blocker, room, objs_start)) continue; // bad placement; need to check elevators for hospitals
		// Note: has doors but no drawers
		unsigned const vanity_obj_ix(objs.size());
		objs.emplace_back(vanity,  TYPE_VANITY,  room_id, dim, !dir, flags, tot_light_amt);
		objs.emplace_back(blocker, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker in front
		add_poi_dim_dir(vanity, room_id, dim, !dir, 0.75); // dscale=0.75
		// add a mirror/medicine cabinet above the sink; shouldn't need to check for overlaps since vanity is placed first
		cube_t mirror(vanity); // start with the sink left and right position
		mirror.z1() = vanity.z2() + 0.14*floor_spacing;
		mirror.z2() = mirror.z1() + 0.30*floor_spacing;
		set_wall_width(mirror, vanity.get_center_dim(!dim), 0.13*floor_spacing, !dim);
		mirror.d[dim][!dir] = wall_edge + dsign*wall_thickness; // thickness
		unsigned const mirror_flags(RO_FLAG_IS_HOUSE | ((num_ext_walls == 1) ? RO_FLAG_INTERIOR : 0)); // flag as interior if windows are opaque glass blocks
		objs.emplace_back(mirror, TYPE_MED_CAB, room_id, dim, !dir, mirror_flags, tot_light_amt);
		set_obj_id(objs); // for crack texture selection/orient
		room.set_has_mirror();

		if (rgen.rand_bool()) { // add a bar of soap next to the vanity sink
			float const one_inch(get_one_inch()), soap_hlen(min(0.5*depth, 2.0*one_inch)), soap_hwidth(1.25*one_inch), soap_height(1.0*one_inch); // 4x2.5x1
			colorRGBA const soap_color(soap_colors[rgen.rand() % NUM_SOAP_COLORS]);
			vector3d const soap_sz((dim ? soap_hwidth : soap_hlen), (dim ? soap_hlen : soap_hwidth), soap_height);
			cube_t const sink(get_sink_cube(objs[vanity_obj_ix]));
			bool const side(rgen.rand_bool());
			float const sink_edge(sink.d[!dim][side]);
			cube_t soap, soap_area(vanity);
			soap_area.d[ dim][!dir ] = vanity.get_center_dim(dim); // front half space
			soap_area.d[!dim][!side] = sink_edge + (side ? 1.0 : -1.0)*soap_hwidth;
			soap_area.d[!dim][ side] = sink_edge + (side ? 1.0 : -1.0)*4.0*one_inch; // move away from sink edge
			gen_xy_pos_for_cube_obj(soap, soap_area, soap_sz, soap_height, rgen, 0); // place_at_z1=0
			objs.emplace_back(soap, TYPE_BAR_SOAP, room_id, dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_ROUNDED_CUBE, soap_color);
		}
		return 1; // success/done
	} // for cand
	return 0;
}

void building_t::add_bathroom_plumbing(room_object_t const &obj) { // only water pipes; drains are added when the plumbing fixture is removed
	assert(has_room_geom());
	bool const dim(obj.dim), dir(obj.dir);
	float const wall_thickness(get_wall_thickness());
	float pipe_radius(0.12*wall_thickness);
	unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_LIT); // RO_FLAG_LIT means the pipe casts shadows
	flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO); // cap the exposed end so that the wall and interior aren't visible
	cube_t pipes[2];
	point pipe_p1;
	pipe_p1[!dim] = obj.get_center_dim(!dim); // starts centered
	pipe_p1[ dim] = obj.d[dim][!dir]; // on the back facing the wall

	if (obj.type == TYPE_TOILET) {
		pipe_p1.z = obj.z1() + 0.52*obj.dz(); // directly runs to the tank, not on the bottom like most toilets
		pipes[0].set_from_point(pipe_p1);
	}
	else if (obj.type == TYPE_SINK) {
		float const spacing(0.2*obj.get_width());
		pipe_p1.z = obj.z1() + 0.78*obj.dz(); // near the top
		pipe_p1[!dim] -= 0.5*spacing;
		pipes[0].set_from_point(pipe_p1);
		pipe_p1[!dim] += 1.0*spacing;
		pipes[1].set_from_point(pipe_p1);
	}
	else if (obj.type == TYPE_URINAL) {
		pipe_p1.z = obj.z1() + 0.75*obj.dz(); // near the top
		pipes[0].set_from_point(pipe_p1);
		pipe_radius *= 1.25; // wider
	}
	else if (obj.type == TYPE_TUB || obj.type == TYPE_SHOWERTUB) {
		// not yet handled because a tub/shower+tub can't be taken
	}
	else {assert(0);} // not a plumbing fixture

	for (unsigned d = 0; d < 2; ++d) { // this loop invalidates obj
		cube_t &pipe(pipes[d]);
		if (pipe.is_all_zeros()) continue;
		pipe.expand_in_dim( dim, 0.5*wall_thickness); // set length
		pipe.expand_in_dim(!dim, pipe_radius);
		pipe.expand_in_dim(2,    pipe_radius);
		if (obj.z1() > ground_floor_z1 && !bcube.contains_cube(pipe)) continue; // outside building; needed for park restrooms
		interior->room_geom->objs.emplace_back(pipe, TYPE_PIPE, obj.room_id, dim, 0, flags, obj.light_amt, SHAPE_CYLIN, COPPER_C); // horizontal
	} // for d
}

bool building_t::add_tp_roll(cube_t const &room, unsigned room_id, float tot_light_amt, bool dim, bool dir, float length, float zval, float wall_pos, bool check_valid_pos) {
	float const diameter(length);
	cube_t tp;
	set_cube_zvals(tp, zval, (zval + diameter));
	set_wall_width(tp, wall_pos, 0.5*length, !dim); // set length
	tp.d[dim][ dir] = room.d[dim][dir]; // against the wall
	tp.d[dim][!dir] = tp  .d[dim][dir] + (dir ? -1.0 : 1.0)*1.1*diameter; // set the diameter with a bit of extra space for clearance
	// Note: not checked against other bathroom objects because the toilet is placed first
	if (check_valid_pos && (!room.contains_cube(tp) || is_obj_placement_blocked(tp, room, 1))) return 0;
	if (has_small_part && !check_if_placed_on_wall(tp, get_room(room_id), dim, dir))           return 0; // need to check for missing walls
	interior->room_geom->objs.emplace_back(tp, TYPE_TPROLL, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, WHITE);
	set_obj_id(interior->room_geom->objs);
	return 1;
}

bool building_t::divide_bathroom_into_stalls(rand_gen_t &rgen, room_t &room, float zval, unsigned room_id,
	float tot_light_amt, unsigned floor, unsigned lights_start, unsigned lights_end)
{
	// Note: assumes no prior placed objects other than ceiling lights
	bool const use_sink_model(0), is_park_restroom(is_restroom());
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	vector3d const tsz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
	float const theight(0.35*floor_spacing), twidth(theight*tsz.y/tsz.z), tlength(theight*tsz.x/tsz.z), stall_depth(2.2*tlength);
	float sheight(0), swidth(0), slength(0), uheight(0), uwidth(0), ulength(0);

	if (use_sink_model) {
		vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SINK)); // L, W, H
		sheight = 0.45*floor_spacing; swidth = sheight*ssz.y/ssz.z; slength = sheight*ssz.x/ssz.z;
	}
	else {
		sheight = 0.36*floor_spacing; swidth = 0.3*floor_spacing; slength = 0.32*floor_spacing;
		//slength = (has_parking_garage ? (tlength + 2.0*wall_thickness) : 0.32*floor_spacing); // align sink drain to toilets for parking garage pipes?
	}
	float stall_width(2.0*twidth), sink_spacing(1.75*swidth);
	// depth_val is around 2*clearance; allow slightly tighter spaces in restaurants and restrooms
	float const depth_val(max(stall_depth, slength)), req_depth((has_house_floorplan() ? 1.6 : 2.0)*depth_val);
	bool br_dim(room.dy() < room.dx()), sink_side(0), sink_side_set(0), mens_room(0); // br_dim is the smaller dim
	cube_t place_area(room), avoid;

	if ((is_park_restroom || is_restaurant()) && zval >= ground_floor_z1) { // handle exterior walls if above ground
		cube_t room_border(parts[room.part_id]);
		room_border.expand_by_xy(-get_trim_thickness()); // add padding around the walls
		place_area = get_walkable_room_bounds(room);
		place_area.intersect_with_cube_xy(room_border);
	}
	else {place_area.expand_by(-0.5*wall_thickness);}
	float const room_width(place_area.get_sz_dim(br_dim));
	if (room_width < req_depth) return 0; // too narrow
	cube_t const place_area_uc(place_area); // capture before clipping

	// check for any stairs or elevators that may partiall overlap a bathroom wall from the adjacent room and avoid them; there should be at most one
	for (stairwell_t const &s : interior->stairwells) {
		if (s.intersects(room) && zval >= s.z1() && zval < s.z2()) {avoid = s; break;}
	}
	for (elevator_t const &e : interior->elevators) {
		if (e.intersects(room) && zval >= e.z1() && zval < e.z2()) {avoid = e; break;}
	}
	// determine men's room vs. women's room
	room_type const init_rtype(room.get_room_type(floor));

	if (room.is_rtype_locked(floor) && (init_rtype == RTYPE_MENS || init_rtype == RTYPE_WOMENS)) { // pre-assigned
		mens_room = (init_rtype == RTYPE_MENS);
	}
	else { // determine if this is the men's or women's room
		if (is_prison() && !interior->has_room_type(RTYPE_BATH)) {mens_room = 1;} // prisons more commonly have men, so select men's room as the first bathroom
		else {
			point const part_center(get_part_for_room(room).get_cube_center()), room_center(room.get_cube_center());
			mens_room = ((part_center.x < room_center.x) ^ (part_center.y < room_center.y));
		}
		if (!is_cube()) { // alternate between men's and women's restrooms
			if      (interior->room_geom->mens_count < interior->room_geom->womens_count) {mens_room = 1;}
			else if (interior->room_geom->mens_count > interior->room_geom->womens_count) {mens_room = 0;}
		}
		else {
			bool has_second_bathroom(0);

			// if there are two bathrooms (one on each side of the building), assign a gender to each side; if only one, alternate gender per floor
			for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
				if (r->part_id != room.part_id || &(*r) == &room) continue; // different part or same room
				if (is_room_office_bathroom(*r, zval, floor)) {has_second_bathroom = 1; break;}
			}
			if (!has_second_bathroom) {mens_room ^= (floor & 1);}
		}
	}
	bool const add_urinals(mens_room && building_obj_model_loader.is_model_valid(OBJ_MODEL_URINAL));
	bool no_end_urinal_divider(0);

	if (add_urinals) { // use urinal model
		vector3d const usz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_URINAL)); // L, W, H
		uheight = 0.4*floor_spacing; uwidth = uheight*usz.y/usz.z; ulength = uheight*usz.x/usz.z;
	}
	if (is_park_restroom) { // restroom only, check exterior doors
		for (tquad_t const &door : doors) {
			cube_t const door_bc(door.get_bcube());
			cube_t door_exp(door_bc);
			door_exp.expand_by_xy(wall_thickness);
			if (!room.intersects(door_exp)) continue;
			bool const door_dim(door_bc.dy() < door_bc.dx());
			
			if (door_dim == br_dim) { // door opens to short end of bathroom, choose a random dir
				sink_side = rgen.rand_bool();
			}
			else {
				sink_side = (door_bc.get_center_dim(!br_dim) > room.get_center_dim(!br_dim)); // side the door is on
				place_area.d[!br_dim][sink_side] += (sink_side ? -1.0 : 1.0)*(door_bc.get_sz_dim(br_dim) - 0.25*swidth); // add sink clearance for the door to close
			}
			sink_side_set = 1;
			break; // should only have one door
		} // for door
	}
	else { // check interior doors
		for (unsigned d = 0; d < 2 && !sink_side_set; ++d) { // try long then short dim
			bool first_sink_side(0);

			if (d == 0 && point_in_mall(room.get_cube_center())) { // place sink on the side closer to the mall concourse
				first_sink_side = (room.get_center_dim(!br_dim) < get_mall_concourse().get_center_dim(!br_dim));
			}
			for (unsigned e = 0; e < 2 && !sink_side_set; ++e) {
				bool const side(bool(e) ^ first_sink_side);
				cube_t c(room);
				set_cube_zvals(c, zval, zval+wall_thickness); // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
				c.d[!br_dim][!side] = c.d[!br_dim][side] + (side ? -1.0 : 1.0)*wall_thickness; // shrink to near zero area in this dim

				for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
					if (i->dim == br_dim) continue; // door in wrong dim
					if (!i->is_connected_to_room(room_id) || !is_cube_close_to_door(c, 0.0, 0, *i, 2)) continue; // check both dirs
					sink_side = side; sink_side_set = 1;
					float const door_width(i->get_sz_dim(br_dim));
					float const door_wall_space(max((i->d[br_dim][0] - place_area.d[br_dim][0]), (place_area.d[br_dim][1] - i->d[br_dim][1])));
					assert(door_wall_space > 0.0); // too strong?

					if (door_wall_space < door_width + slength + get_min_front_clearance()) { // add sink clearance for the door to close if needed
						cube_t const door_bc(i->get_open_door_path_bcube()); // add full clearance if sink in door path
						bool const in_door_path(door_bc.d[br_dim][0] < place_area.d[br_dim][0]+slength || door_bc.d[br_dim][1] > place_area.d[br_dim][1]-slength);
						place_area.d[!br_dim][sink_side] += (sink_side ? -1.0 : 1.0)*(door_width - (in_door_path ? 0.0 : 0.25*swidth)); // allow closer if not blocked
					}
					else {no_end_urinal_divider = 1;}
					break; // sinks are on the side closest to the first door
				}
			} // for e
			if (d == 0 && !sink_side_set) {br_dim ^= 1;} // door not found on long dim - R90 and try short dim
		} // for d
	}
	//if (!sink_side_set) return 0;
	assert(sink_side_set);
	float const room_len(place_area.get_sz_dim(!br_dim));
	float sinks_len((/*is_park_restroom ? 0.4 :*/ 0.4)*room_len), stalls_len(room_len - sinks_len);
	if (stalls_len <= 0.0) return 0; // shouldn't happen
	if (sinks_len < 2.0*sink_spacing) {sink_spacing *= 0.8;} // reduce sink spacing a bit to try and fit at least two

	if (sinks_len < 2.0*sink_spacing) { // still too small; allocate more space to sinks/urinals and less to stalls; target 2 of each
		float const targ_sinks_len(2.01*sink_spacing), new_stalls_len(room_len - targ_sinks_len);
		if (new_stalls_len > 2.0*stall_width) {sinks_len = targ_sinks_len; stalls_len = new_stalls_len;}
	}
	unsigned const num_stalls(std::floor(stalls_len/stall_width)), num_sinks(std::floor(sinks_len/sink_spacing));
	if (num_stalls < 2 || num_sinks < 1) return 0; // not enough space for 2 stalls and a sink
	stall_width  = stalls_len/num_stalls; // reclaculate to fill the gaps
	sink_spacing = sinks_len/num_sinks;
	bool const two_rows(room_width > 3.0*depth_val); // could be as small as 2.6x
	// put stalls on back wall rather than sides if short and wide? but most likely the end opposite the door is an exterior wall, and we can't put stalls there
	//if (!two_rows && room_width > stalls_len && stalls_len > stall_depth && zval > ground_floor_z1 && !is_prison()) {}
	bool skip_stalls_side(room_id & 1); // put stalls on a side consistent across floors
	if (classify_room_wall(room, zval, br_dim, !skip_stalls_side, 0) == ROOM_WALL_EXT) {skip_stalls_side ^= 1;} // no stalls/sinks along exterior walls
	unsigned const showers_dir(is_prison() ? skip_stalls_side : 2); // one side of prison bathroom has showers rather than stalls
	float const sink_side_sign(sink_side ? 1.0 : -1.0), stall_step(sink_side_sign*stall_width), sink_step(-sink_side_sign*sink_spacing);
	float const floor_thickness(get_floor_thickness());
	unsigned const NUM_STALL_COLORS = 4;
	colorRGBA const stall_colors[NUM_STALL_COLORS] = {colorRGBA(0.75, 1.0, 0.9, 1.0), colorRGBA(0.7, 0.8, 1.0), WHITE, DK_GRAY}; // blue-green, light blue
	colorRGBA const stall_color(stall_colors[(3*interior->doors.size() + 7*mat_ix) % NUM_STALL_COLORS]); // random, but constant for each building
	room_obj_shape const stall_shape(room.has_tall_ceil(floor_spacing) ? SHAPE_TALL : SHAPE_CUBE); // tall for tall ceiling rooms like restaurants
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());
	room_object_t mirrors[2]; // candidate mirrors for each dir
	++(mens_room ? interior->room_geom->mens_count : interior->room_geom->womens_count);

	for (unsigned dir = 0; dir < 2; ++dir) { // each side of the wall
		if (!two_rows && dir == (unsigned)skip_stalls_side) continue; // no stalls/sinks on this side
		// add stalls or showers
		float const dir_sign(dir ? -1.0 : 1.0), wall_pos(place_area.d[br_dim][dir]), stall_from_wall(wall_pos + dir_sign*(0.5*tlength + wall_thickness));
		float stall_pos(place_area.d[!br_dim][!sink_side] + 0.5*stall_step);
		bool last_ooo(0);

		for (unsigned n = 0; n < num_stalls; ++n, stall_pos += stall_step) {
			//if (n == 0 && num_stalls > 2) {} // what about larger handicap stall?
			point center(stall_from_wall, stall_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t stall(center);
			stall.z2() = stall.z1() + floor_spacing - floor_thickness; // set stall height to room height
			stall.expand_in_dim(!br_dim, 0.5*stall_width);
			stall.d[br_dim][ dir] = wall_pos; // + wall_thickness?
			stall.d[br_dim][!dir] = wall_pos + dir_sign*stall_depth; // extend the front outward from the wall
			if (interior->is_cube_close_to_doorway(stall, room, 0.0, 1)) continue; // skip if close to a door (for rooms with doors at both ends); inc_open=1
			if (!avoid.is_all_zeros() && stall.intersects(avoid))        continue;

			if (!is_cube()) {
				cube_t stall_exp(stall);
				stall_exp.d[br_dim][!dir] += dir_sign*1.0*get_scaled_player_radius(); // make sure the player can fit
				if (!check_cube_within_part_sides(stall_exp)) continue; // outside the building
			}
			bool const is_open(rgen.rand_bool()); // 50% chance of stall door being open
			bool const out_of_order(!last_ooo && !is_open && !room.is_ext_basement() && rgen.rand_float() < 0.2); // not for mall bathrooms
			unsigned flags(out_of_order ? RO_FLAG_BROKEN : 0); // toilet can't be flushed and door can't be opened if out of order
			last_ooo = out_of_order; // no two OOO in a row

			if (dir == showers_dir) { // add shower rather than stall
				flags |= RO_FLAG_IN_JAIL;
				
				if (rgen.rand_bool()) { // maybe add a bar of soap
					float const one_inch(get_one_inch()), soap_hlen(2.0*one_inch), soap_hwidth(1.25*one_inch), soap_height(1.0*one_inch); // 4x2.5x1
					colorRGBA const soap_color(soap_colors[rgen.rand() % NUM_SOAP_COLORS]);
					vector3d const soap_sz((!br_dim ? soap_hwidth : soap_hlen), (!br_dim ? soap_hlen : soap_hwidth), soap_height);
					float const signed_depth(dir_sign*stall_depth), inner_wall_pos(stall.d[br_dim][dir] + 0.02*signed_depth);
					cube_t shelf(stall); // from building_room_geom_t::add_br_stall()
					shelf.d[br_dim][ dir] = inner_wall_pos + 0.02*signed_depth; // don't clip through tile
					shelf.d[br_dim][!dir] = inner_wall_pos + 0.10*signed_depth;
					shelf.z2() = stall.z1() + 0.28*stall.dz();
					cube_t soap;
					gen_xy_pos_for_cube_obj(soap, shelf, soap_sz, soap_height, rgen, 0); // place_at_z1=0
					objs.emplace_back(soap, TYPE_BAR_SOAP, room_id, !br_dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_ROUNDED_CUBE, soap_color);
				}
			}
			else { // toilet stall
				cube_t toilet(center);
				toilet.expand_in_dim( br_dim, 0.5*tlength);
				toilet.expand_in_dim(!br_dim, 0.5*twidth);
				toilet.z2() += theight;
				objs.emplace_back(toilet, TYPE_TOILET, room_id, br_dim, !dir, flags, tot_light_amt);
				add_bathroom_plumbing(objs.back());
				float const tp_length(0.18*theight), wall_pos(toilet.get_center_dim(br_dim));
				cube_t stall_inner(stall);
				stall_inner.expand_in_dim(!br_dim, -0.0125*stall.dz()); // subtract off stall wall thickness
				add_tp_roll(stall_inner, room_id, tot_light_amt, !br_dim, dir, tp_length, (zval + 0.7*theight), wall_pos);
			}
			assert(stall.is_strictly_normalized());
			objs.emplace_back(stall, TYPE_STALL, room_id, br_dim, dir, (flags | (is_open ? RO_FLAG_OPEN : 0)), tot_light_amt, stall_shape, stall_color);
			if (dir == showers_dir) {objs.back().obj_id = room_id;} // sets shower style (handles); matches across rooms

			if (out_of_order) { // add out-of-order sign
				cube_t stall_clipped(stall);
				stall_clipped.z2() -= 0.33*stall.dz(); // make it shorter; really only need the stall door, but the dim dimension size isn't used anyway
				add_out_or_order_sign(stall_clipped, br_dim, !dir, room_id);
			}
			// move any lights that intersect stalls, since stalls extend to the ceiling
			for (auto i = objs.begin()+lights_start; i < objs.begin()+lights_end; ++i) {
				if (i->type != TYPE_LIGHT || !i->intersects(stall)) continue;
				float const move_amt(stall.d[br_dim][!dir] - i->d[br_dim][dir] + dir_sign*wall_thickness);
				i->translate_dim(br_dim, move_amt); // shouldn't intersect anything else
			}
		} // for n
		if (add_urinals && dir == (unsigned)skip_stalls_side) continue; // no urinals and sinks are each on one side
		// add sinks
		float const sink_start(place_area.d[!br_dim][sink_side] + 0.5f*sink_step);
		float const sink_from_wall(wall_pos + dir_sign*(0.5f*slength + (use_sink_model ? wall_thickness : 0.0f)));
		float sink_pos(sink_start);
		bool hit_mirror_end(0);
		unsigned last_sink_ix(0);
		cube_t sinks_bcube;

		for (unsigned n = 0; n < num_sinks; ++n, sink_pos += sink_step) {
			point center(sink_from_wall, sink_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t sink(center);
			sink.expand_in_dim(br_dim, 0.5*slength);
			sink.z2() += sheight;
			if (interior->is_cube_close_to_doorway(sink, room, 0.0, 1)) continue; // skip if close to a door, inc_open=1, pre expand
			sink.expand_in_dim(!br_dim, 0.5*(use_sink_model ? swidth : fabs(sink_step))); // tile exactly with the adjacent sink
			if (interior->is_cube_close_to_doorway(sink, room, 0.0, 0)) continue; // skip if close to a door
			if (!check_cube_within_part_sides(sink))                    continue; // outside the building
			if (!avoid.is_all_zeros() && sink.intersects(avoid))        continue;
			if (use_sink_model) {objs.emplace_back(sink, TYPE_SINK,   room_id, br_dim, !dir, 0, tot_light_amt);} // sink 3D model
			else                {objs.emplace_back(sink, TYPE_BRSINK, room_id, br_dim, !dir, 0, tot_light_amt);} // flat basin sink
			if (use_sink_model) {add_bathroom_plumbing(objs.back());}
			// if we started the mirror, but we have a gap with no sink (blocked by a door, etc.), then end the mirror
			hit_mirror_end |= (n > last_sink_ix+1 && !sinks_bcube.is_all_zeros());
			if (!hit_mirror_end) {sinks_bcube.assign_or_union_with_cube(sink);}
			last_sink_ix = n;
		} // for n
		if (add_urinals) {
			// add urinals opposite the sinks, using same spacing as sinks;
			// should urinals be on the same wall as sinks between sinks and stalls when long and thin such as mall restrooms?
			bool const is_tight(room_width < 1.3*req_depth);
			float const div_wall_len(0.25*floor_spacing);
			bool const extra_urinals(!two_rows && room_width > 2.0*depth_val + max(ulength, div_wall_len)); // add urinals opposite stalls if wide enough
			unsigned num_urinals((extra_urinals ? 2 : 1)*num_sinks); // should be able to fit twice as many
			float const u_wall(place_area.d[br_dim][!dir]), u_from_wall(u_wall - dir_sign*(0.5*ulength + 0.01*wall_thickness));
			float u_pos(sink_start);
			cube_t sep_wall;
			set_cube_zvals(sep_wall, zval+0.15*uheight, zval+1.25*uheight);
			sep_wall.d[br_dim][!dir] = u_wall;
			sep_wall.d[br_dim][ dir] = u_wall - dir_sign*div_wall_len;
			// school bathrooms have a short urinal at one end if there is more than one urinal
			unsigned const short_urinal_ix((is_school() && num_urinals > 1) ? ((room_id & 1) ? 0 : num_urinals-1) : num_urinals);
			bool last_added_urinal(0);

			for (unsigned n = 0; n < num_urinals; ++n, u_pos += sink_step) {
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				point center(u_from_wall, u_pos, (zval + ((n == short_urinal_ix) ? 0.125 : 0.25)*uheight));
				if (br_dim) {swap(center.x, center.y);} // R90 about z
				cube_t urinal(center);
				urinal.expand_in_dim( br_dim, 0.5*ulength);
				urinal.expand_in_dim(!br_dim, 0.5*uwidth);
				urinal.z2() += uheight;
				last_added_urinal = 0; // in case we continue below
				if (interior->is_cube_close_to_doorway(urinal, room, 0.0, 1)) continue; // skip if close to a door
				if (!check_cube_within_part_sides(urinal))                    continue; // outside the building
				if (!avoid.is_all_zeros() && urinal.intersects(avoid))        continue;
				if (is_park_restroom && is_cube_close_to_exterior_doorway(urinal, 0.0, 0)) continue; // inc_open=0
				
				if (no_end_urinal_divider && n == 0) {} // no starting divider if against a wall
				else if (!interior->is_cube_close_to_doorway(sep_wall, room, 0.0, 1)) { // check for doors, when the bathroom door is not centered on the room
					if (!is_park_restroom || !is_cube_close_to_exterior_doorway(sep_wall, 0.0, 0)) {
						if (room.contains_cube_xy(sep_wall)) { // needed for park restroom
							objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, RO_FLAG_HANGING, tot_light_amt, SHAPE_SHORT, stall_color);
						}
					}
				}
				objs.emplace_back(urinal, TYPE_URINAL, room_id, br_dim, dir, 0, tot_light_amt);
				add_poi_dim_dir(urinal, room_id, br_dim, dir, 0.5); // dscale=0.5
				add_bathroom_plumbing(objs.back());
				last_added_urinal = 1;
			} // for n
			if (!two_rows && !is_tight && last_added_urinal) { // skip first wall if adjacent to a stall, or opposite a stall if space is tight
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				if (!avoid.is_all_zeros() && sep_wall.intersects(avoid)) {} // skip
				else if (is_park_restroom && is_cube_close_to_exterior_doorway(sep_wall, 0.0, 0)) {} // skip
				else {objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, RO_FLAG_HANGING, tot_light_amt, SHAPE_SHORT, stall_color);}
			}
		}
		if (!sinks_bcube.is_all_zeros()) { // add a long mirror above the sink
			cube_t mirror(sinks_bcube);
			mirror.expand_in_dim(!br_dim, -0.25*wall_thickness); // slightly smaller
			mirror.d[br_dim][ dir] = wall_pos;
			mirror.d[br_dim][!dir] = wall_pos + dir_sign*0.1*wall_thickness;
			mirror.z1() = sinks_bcube.z2() + 0.25*floor_thickness;
			mirror.z2() = zval + 0.9*floor_spacing - floor_thickness;
			if (mirror.is_strictly_normalized()) {mirrors[dir] = room_object_t(mirror, TYPE_MIRROR, room_id, br_dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);}
			add_poi_dim_dir(sinks_bcube, room_id, br_dim, !dir, 0.75); // dscale=0.75
		}
		// add wall items on each side that has sinks
		float const side_wall_pos(place_area_uc.d[!br_dim][sink_side]), sink_edge(place_area.d[!br_dim][sink_side]);
		float const z1(zval + max((sheight + 0.15*floor_spacing), 0.5*floor_spacing)); // for paper towel and soap dispensers
		cube_t gap_area(place_area_uc);
		gap_area.d[!br_dim][!sink_side] = sink_edge; // ends at edge of sink
		float const gap_width(gap_area.get_sz_dim(!br_dim));
		float hd_dist_from_wall(((gap_width > 0.25*floor_spacing) ? 0.3 : 1.0)*slength); // start at sink length (depth) if no gap, less if there's a gap

		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_TOWEL_DISP)) { // paper towel dispenser
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOWEL_DISP)); // D, W, H
			float const height(0.21*floor_spacing), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z);
			bool const place_on_sink_wall(gap_width > 3.2*hwidth);
			cube_t td;
			set_cube_zvals(td, z1, z1+height);
			bool td_dim(0), td_dir(0);

			if (place_on_sink_wall) { // have space to place between sink and door wall
				td_dim = br_dim; td_dir = !dir;
				td.d[td_dim][!td_dir] = wall_pos;
				td.d[td_dim][ td_dir] = wall_pos + dir_sign*depth;
				set_wall_width(td, (side_wall_pos - sink_side_sign*2.0*hwidth), hwidth, !td_dim); // close to side wall
			}
			else { // not enough space; try side wall next to door
				td_dim = !br_dim; td_dir = !sink_side;
				td.d[td_dim][!td_dir] = side_wall_pos;
				td.d[td_dim][ td_dir] = side_wall_pos - sink_side_sign*depth;
				set_wall_width(td, (wall_pos + dir_sign*3.0*hwidth), hwidth, !td_dim); // width from back wall
			}
			if (!is_cube_close_to_doorway(td, room, 0.0, 2, 1)) { // inc_open=2 (interior only), check_open_dir=1
				objs.emplace_back(td, TYPE_TOWEL_DISP, room_id, td_dim, td_dir, 0, tot_light_amt);
				if ( place_on_sink_wall) {gap_area.d[!br_dim][!sink_side] = td.d[!br_dim][!sink_side] - sink_side_sign*hwidth;} // shrink gap for soap dispenser placement
				if (!place_on_sink_wall) {hd_dist_from_wall = 4.0f*hwidth;} // update if placed on side wall (no/small gap case)
			}
		}
		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_SOAP_DISP)) { // soap dispenser
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SOAP_DISP)); // D, W, H
			float const height(0.12*floor_spacing), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z);

			if (gap_area.get_sz_dim(!br_dim) > 3.2*hwidth) { // have space to place between sink and [paper towel dispenser or door wall]
				bool const sd_dim(br_dim), sd_dir(!dir);
				cube_t sd;
				set_cube_zvals(sd, z1, z1+height);
				sd.d[sd_dim][!sd_dir] = wall_pos;
				sd.d[sd_dim][ sd_dir] = wall_pos + dir_sign*depth;
				set_wall_width(sd, (sink_edge + sink_side_sign*2.0*hwidth), hwidth, !sd_dim); // close to sink
				// it really shouldn't be possible for the soap dispenser to intersect a door or another object; no check needed?
				objs.emplace_back(sd, TYPE_SOAP_DISP, room_id, sd_dim, sd_dir, 0, tot_light_amt); // nocoll?
			}
		}
		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_HAND_DRYER)) { // hand dryer
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_HAND_DRYER)); // D, W, H
			float const height(0.16*floor_spacing), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z), hd_z1(zval + 0.36*floor_spacing);
			bool const hd_dim(!br_dim), hd_dir(!sink_side);
			cube_t hd;
			set_cube_zvals(hd, hd_z1, hd_z1+height); // below the park restroom window
			hd.d[hd_dim][!hd_dir] = side_wall_pos;
			hd.d[hd_dim][ hd_dir] = side_wall_pos - sink_side_sign*depth;
			set_wall_width(hd, (wall_pos + dir_sign*(hd_dist_from_wall + 3.0*hwidth)), hwidth, !hd_dim);

			if (!is_cube_close_to_doorway(hd, room, 0.0, 2, 1)) { // inc_open=2 (interior only), check_open_dir=1
				// try to push closer to the door
				float const stop_pos(room.get_center_dim(br_dim) - dir_sign*hwidth), step_val(dir_sign*hwidth);

				for (unsigned n = 0; n < 10; ++n) { // up to 10 steps
					cube_t hd2(hd);
					hd2.translate_dim(br_dim, step_val);
					if (hd2.d[br_dim][1] > stop_pos && hd2.d[br_dim][0] < stop_pos) break; // reached stop_pos; done
					cube_t hd_tc(hd2);
					hd_tc.expand_in_dim(br_dim, 0.25*hwidth); // extra padding for door frame
					if (is_cube_close_to_doorway(hd_tc, room, 0.0, 2, 1)) break; // can't step this far
					hd = hd2; // valid step
				} // for n
				objs.emplace_back(hd, TYPE_HAND_DRYER, room_id, hd_dim, hd_dir, 0, tot_light_amt);
				cube_t hd_blocker(hd); // add a blocker under the hand dryer to avoid placing a trashcan there
				set_cube_zvals(hd_blocker, hd_z1-height, hd_z1);
				objs.emplace_back(hd_blocker, TYPE_BLOCKER, room_id, hd_dim, hd_dir, 0, tot_light_amt);
			}
		}
	} // for dir
	for (unsigned d = 0; d < 2; ++d) { // each candidate mirror
		if (mirrors[d].is_all_zeros()) continue;
		if (ENABLE_MIRROR_REFLECTIONS && d == (unsigned)skip_stalls_side && !mirrors[!d].is_all_zeros()) continue; // select a single side if reflections are enabled
		objs.push_back(mirrors[d]);
		set_obj_id(objs); // for crack texture selection/orient
		room.set_has_mirror();
	}
	// add a large round trashcan; should we not add a smaller trashcan if this is placed?
	add_corner_trashcans(rgen, room, zval, room_id, tot_light_amt, objs_start, !br_dim, 0); // both_ends=0
	room_type const rtype(mens_room ? RTYPE_MENS : RTYPE_WOMENS);
	room.assign_to(rtype, floor);
	room.set_has_br_stalls();
	
	// make sure doors start closed and unlocked, and flag them as auto_close;
	// if (!is_cube()) we also want to make sure the door opens inward, but we can't change it for only one door in the stack
	for (door_stack_t const &ds : interior->door_stacks) {
		if (!ds.is_connected_to_room(room_id)) continue;
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() > zval || door.z2() < zval) continue; // wrong floor
			door.rtype  = rtype; // tag type so that the correct sign is added
			door.locked = 0; // only needed for non-cube buildings
			if (!door.is_bars()) {door.make_auto_close();} // what about DOOR_TYPE_JAIL?
			if ((!is_cube() || is_prison()) && !door_opens_inward(ds, room)) {door.opens_out_of_br = 1;}
		} // for ds
	} // for ds
	// add a sign outside the bathroom door; should schools use boys/girls? but it doesn't match the text on the signs placed on the doors
	string const sign_text(/*is_school() ? (mens_room ? "Boys" : "Girls") :*/ (mens_room ? "Men" : "Women"));
	add_door_sign(sign_text, room, zval, room_id, 1); // no_check_adj_walls=1

	// make this door/room out of order 10% of the time; only for cube buildings (others need the connectivity), and not for mall or restaurant bathrooms
	if (is_cube() && !(has_mall() && room.is_ext_basement()) && !is_restaurant() && !is_restroom() && rgen.rand_float() < 0.1) {
		// there may be multiple doors, and they must all have the OOO sign
		vect_door_stack_t &doorways(get_doorways_for_room(room, zval));
		for (door_stack_t &ds : doorways) {make_door_out_of_order(room, zval, room_id, ds);}
		room.set_has_out_of_order(); // flag if any floor is out of order
	}
	return 1;
}

void building_t::add_bathroom_window(cube_t const &window, bool dim, bool dir, unsigned room_id, unsigned floor) { // frosted window blocks, for houses or office buildings
	if (!has_int_windows()) return; // no interior (or exterior) drawn windows
	room_t const &room(get_room(room_id));
	// exterior looks odd to have window block walls at the corner of a building,
	// so only enable this for single exterior walls, when there are no exterior windows, when there are stalls, or for industrial/rest/conv bathrooms (which look odd without it)
	if (!is_industrial() && !is_restaurant() && !is_conv_store() && !is_restroom() && has_windows() && !room.has_br_stalls() && count_ext_walls_for_room(room, window.z1()) != 1) return;
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c(window);
	c.translate_dim(dim, (dir ? 1.0 : -1.0)*0.5*get_trim_thickness()); // half the previous translate to prevent Z-fighting in mirror reflections
	unsigned flags(RO_FLAG_NOCOLL | (room.is_lit_on_floor(floor) ? RO_FLAG_LIT : 0));
	if (has_windows()) {flags |= RO_FLAG_HAS_EXTRA;} // draw exterior face if there are exterior windows
	objs.emplace_back(c, TYPE_WINDOW, room_id, dim, dir, flags, 1.0, SHAPE_CUBE, WHITE); // always lit
}

void building_t::create_restroom_floorplan(unsigned part_id, rand_gen_t &rgen) {
	assert(part_id < parts.size());
	cube_t const &part(parts[part_id]);
	assert(interior->rooms.empty()); // must call this first
	// divide into men's and women's rooms
	mw_restroom_side = rgen.rand_bool();
	bool const dim(!get_street_dim()); // split dim
	float const wall_thick(get_wall_thickness()), split_pos(part.get_center_dim(dim));
	cube_t wall(part);
	create_wall(wall, dim, split_pos, get_fc_thickness(), 0.5*wall_thick, 0.05*wall_thick);
	interior->walls[dim].push_back(wall);

	for (unsigned d = 0; d < 2; ++d) {
		cube_t room(part);
		room.d[dim][!d] = split_pos;
		add_assigned_room(room, part_id, ((bool(d) ^ mw_restroom_side) ? RTYPE_MENS : RTYPE_WOMENS));
		interior->rooms.back().is_single_floor = 1; // probably not needed
	}
	if (street_side) { // doors at front case; add a walkway in front in the form of a driveway
		bool const dir(get_street_side());
		float const front_wall(bcube.d[!dim][dir]);
		set_cube_zvals(driveway, ground_floor_z1, ground_floor_z1+0.5*get_fc_thickness());
		driveway.d[!dim][!dir] = front_wall;
		driveway.d[!dim][ dir] = front_wall + (dir ? 1.0 : -1.0)*0.75*get_window_vspace(); // extend outward in front of the building
		for (unsigned d = 0; d < 2; ++d) {driveway.d[dim][d] = bcube.d[dim][d];}
		driveway.expand_in_dim(dim, -0.5*wall_thick); // shrink inward slightly inside the exterior fences/walls
	}
}

