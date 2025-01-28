// 3D World - Factory Buildings
// by Frank Gennari 1/26/25

#include "function_registry.h"
#include "buildings.h"

extern float CAMERA_RADIUS;


float shift_val_to_not_intersect_window(cube_t const &c, float val, float hspace, float window_border, bool dim);

void building_t::create_factory_floorplan(unsigned part_id, float window_hspacing[2], float window_border, rand_gen_t &rgen) {
	assert(part_id < parts.size());
	cube_t const &part(parts[part_id]);
	vector<room_t> &rooms(interior->rooms);
	bool const dim(part.dx() < part.dy()), dir(rgen.rand_bool()); // long dim
	float const door_width(get_doorway_width()), door_hwidth(0.5*door_width), floor_spacing(get_window_vspace());
	float const wall_thick(get_wall_thickness()), wall_hthick(0.5*wall_thick), floor_thick(get_floor_thickness()), fc_thick(get_fc_thickness());
	float const door_ent_pad(2.0*door_width), room_len(part.get_sz_dim(dim)), room_width(part.get_sz_dim(!dim)), dsign(dir ? -1.0 : 1.0);
	float const sub_room_len(max(1.5*floor_spacing, min(3.0*floor_spacing, 0.2*room_len))*rgen.rand_uniform(0.9, 1.0));
	float const centerline(part.get_center_dim(!dim) + (rgen.rand_bool() ? 1.0 : -1.0)*rgen.rand_uniform(0.15, 0.25)*room_width); // biased to a random side
	unsigned const num_floors(calc_num_floors(part, floor_spacing, floor_thick));
	assert(num_floors >= 2); // main factory must be at least 2 floors tall
	// add bathroom and office to either side of a potential placement of the front entrance door
	float split_pos(part.d[dim][dir] + dsign*sub_room_len);
	split_pos = shift_val_to_not_intersect_window(part, split_pos, window_hspacing[dim], window_border, dim);
	cube_t sub_rooms(part), floor_space(part);
	sub_rooms.z2() = part.z1() + floor_spacing;
	sub_rooms.d[dim][!dir] = floor_space.d[dim][dir] = split_pos;
	cube_t sub_room[2] = {sub_rooms, sub_rooms}; // to each side of the entrance
	float wall_edge[2] = {(centerline - door_ent_pad), (centerline + door_ent_pad)};
	cube_t entrance_area(sub_rooms);

	for (unsigned d = 0; d < 2; ++d) { // determine wall positions to avoid intersecting windows
		wall_edge[d] = shift_val_to_not_intersect_window(part, wall_edge[d], window_hspacing[!dim], window_border, !dim);
		sub_room[d].d[!dim][!d] = entrance_area.d[!dim][d] = wall_edge[d];
	}
	floor_space.expand_in_z(-fc_thick); // shrink to remove ceiling and floor
	float const entrance_pos(0.5*(wall_edge[0] + wall_edge[1]));
	interior->factory_info.reset(new bldg_factory_info_t(dim, dir, entrance_pos, floor_space, entrance_area));
	bool const larger_room(sub_room[0].get_sz_dim(!dim) < sub_room[1].get_sz_dim(!dim));

	for (unsigned d = 0; d < 2; ++d) { // add walls, doors, and ceilings/floors
		// TODO: concrete interior walls?
		bool const is_larger(bool(d) == larger_room), is_bathroom(!is_larger);
		cube_t const &r(sub_room[d]);
		cube_t wall_z_range(r);
		wall_z_range.translate_dim(2, fc_thick); // extend up to cover the floor above the room
		cube_t lwall(wall_z_range), swall(wall_z_range), lwall2, swall2;
		swall.d[ dim][!dir] += dsign*wall_hthick; // extend out
		lwall.d[!dim][!d  ] += (d ? -1.0 : 1.0)*wall_hthick; // extend out
		set_wall_width(lwall, split_pos,     wall_hthick,  dim);
		set_wall_width(swall, r.d[!dim][!d], wall_hthick, !dim);

		// add door(s) and walls
		for (unsigned e = 0; e < 2; ++e) {
			bool const wdim(dim ^ bool(e)), ddir(e ? d : dir);
			cube_t wall(e ? swall : lwall); // copy so that it can be clipped
			
			if (is_larger || e == 1) { // no door in factory floor side of bathroom
				float const wall_center(wall.get_center_dim(!wdim)), door_lo(wall_center - door_hwidth), door_hi(wall_center + door_hwidth);
				insert_door_in_wall_and_add_seg(wall, door_lo, door_hi, !wdim, ddir, 0, is_bathroom);
				interior->door_stacks.back().set_mult_floor(); // counts as multi-floor (for drawing top edge)
				interior->walls[wdim].push_back(wall);
				// add section of wall above the door
				wall.d[!wdim][0] = door_lo;
				wall.d[!wdim][1] = door_hi;
				wall.z1() = interior->doors.back().z2();
			}
			interior->walls[wdim].push_back(wall);
		} // for e
		cube_t room_inner(r);
		room_inner.d[ dim][!dir] -= dsign*wall_hthick; // extend in
		room_inner.d[!dim][!d  ] -= (d ? -1.0 : 1.0)*wall_hthick; // extend in
		// add ceiling and floor
		cube_t room_ceil(room_inner), room_floor(r); // floor above the room
		room_ceil .z1() = r.z2() - fc_thick;
		room_floor.z1() = r.z2();
		room_floor.z2() = r.z2() + fc_thick;
		room_floor.d[ dim][!dir] = swall.d[ dim][!dir];
		room_floor.d[!dim][!d  ] = lwall.d[!dim][!d  ];
		interior->ceilings.push_back(room_ceil );
		interior->floors  .push_back(room_floor);
		// add room itself; will overlap main factory room
		add_room(room_inner, part_id, 2); // 2 lights; not a typical office building office
		rooms.back().assign_all_to(is_larger ? RTYPE_OFFICE : RTYPE_BATH); // office or bathroom
	} // for d
	// should there be an entryway room, then the factory doesn't overlap the sub-rooms? but then there will be empty space above them
	// add entire part as a room (factory floor); must be done last so that smaller contained rooms are picked up in early exit queries (model occlusion, light toggle, door conn)
	add_room(part, part_id); // num_lights will be calculated later
	rooms.back().assign_all_to(RTYPE_FACTORY);
	rooms.back().is_single_floor = 1;
}

void building_t::add_factory_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id) {
	assert(interior->factory_info);
	float const light_amt(1.0); // always lit?
	bool const edim(interior->factory_info->entrance_dim), edir(interior->factory_info->entrance_dir), beam_dim(!edim); // edim is the long dim
	float const window_vspace(get_window_vspace()), wall_thick(get_wall_thickness()), fc_thick(get_fc_thickness());
	vect_room_object_t &objs(interior->room_geom->objs);
	// add support pillars around the exterior, between windows; add ceiling beams
	float const support_width(FACTORY_BEAM_THICK*wall_thick), support_hwidth(0.5*support_width);
	float const ceil_zval(room.z2() - fc_thick), beams_z1(ceil_zval - support_width);
	cube_t support_bounds(room);
	support_bounds.expand_by_xy(-support_hwidth);
	cube_t support, beam;
	set_cube_zvals(support, zval,     beams_z1 );
	set_cube_zvals(beam,    beams_z1, ceil_zval);
	vect_cube_t support_parts, nested_rooms;
	float const shift_vals[6] = {-0.1, 0.2, -0.3, 0.4, -0.5, 0.6}; // cumulative version of {-0.1, 0.1, -0.2, 0.2, -0.3, 0.3}; not enough shift to overlap a window

	for (room_t const &r : interior->rooms) {
		if (!r.intersects_no_adj(room)) break; // done with above ground factory rooms
		if (r != room) {nested_rooms.push_back(r);} // skip self (factory)
	}
	for (unsigned dim = 0; dim < 2; ++dim) {
		unsigned const num_windows(get_num_windows_on_side(room, !dim));
		float const spacing(room.get_sz_dim(!dim)/num_windows);

		for (unsigned dir = 0; dir < 2; ++dir) { // walls
			float const wall_pos(room.d[dim][dir]);
			support.d[dim][ dir] = wall_pos;
			support.d[dim][!dir] = wall_pos + (dir ? -1.0 : 1.0)*support_width;

			for (unsigned n = 0; n <= num_windows; ++n) {
				if (dim == 1 && (n == 0 || n == num_windows)) continue; // don't duplicate add corners
				float const centerline(max(support_bounds.d[!dim][0], min(support_bounds.d[!dim][1], room.d[!dim][0] + n*spacing))); // clamp to support_bounds
				set_wall_width(support, centerline, support_hwidth, !dim);
				cube_t test_cube(support);
				test_cube.expand_by_xy(wall_thick);
				bool valid(0);

				// check and either move or skip if blocked by exterior door or basement stairs
				for (unsigned m = 0; m < 6; ++m) {
					if (!cube_int_ext_door(test_cube) && !interior->is_blocked_by_stairs_or_elevator(test_cube)) {valid = 1; break;}
					float const shift(shift_vals[m]*spacing);
					support  .translate_dim(!dim, shift);
					test_cube.translate_dim(!dim, shift);
				}
				if (!valid) continue; // failed, skip
				unsigned skip_faces(~get_face_mask(dim, dir));
				if (n == 0          ) {skip_faces |= ~get_face_mask(!dim, 0);}
				if (n == num_windows) {skip_faces |= ~get_face_mask(!dim, 1);}
				objs.emplace_back(support, TYPE_IBEAM, room_id, dim, 1, 0, light_amt, SHAPE_CUBE, WHITE, skip_faces); // vertical

				for (cube_t const &r : nested_rooms) { // clip in Z if intersects a room
					if (!r.intersects(support)) continue;
					cube_t bot(support);
					objs.back().z1() = r.z2() + fc_thick; // top of room cube
					bot.z2() = objs.back().z1() + 0.1*fc_thick; // shift up slightly
					cube_t sub(r);
					sub.expand_by(0.5*wall_thick);
					subtract_cube_from_cube(bot, sub, support_parts, 1); // clear_out=1
					for (cube_t const &c : support_parts) {objs.emplace_back(c, TYPE_IBEAM, room_id, dim, 1, RO_FLAG_ADJ_TOP, light_amt, SHAPE_CUBE, WHITE, skip_faces);}
					break;
				} // for r
			} // for n
		} // for dir
		// add beams
		unsigned const num_hdiv(2*num_windows); // add intermediate beams and hang lights on them
		for (unsigned d = 0; d < 2; ++d) {beam.d[dim][d] = room.d[dim][d];}
		if (bool(dim) == beam_dim) {beam.expand_in_dim(beam_dim, -support_hwidth);} // half overlap of vert supports

		for (unsigned n = 0; n <= num_hdiv; ++n) {
			if (bool(dim) != beam_dim && n > 0 && n < num_hdiv) continue; // only add edge beams in this dim
			float const centerline(max(support_bounds.d[!dim][0], min(support_bounds.d[!dim][1], room.d[!dim][0] + n*0.5f*spacing))); // clamp to support_bounds
			set_wall_width(beam, centerline, support_hwidth, !dim);
			unsigned skip_faces(EF_Z2);
			if (n == 0       ) {skip_faces |= ~get_face_mask(!dim, 0);}
			if (n == num_hdiv) {skip_faces |= ~get_face_mask(!dim, 1);}
			objs.emplace_back(beam, TYPE_IBEAM, room_id, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, WHITE, skip_faces);
		} // for n
	} // for dim
#if 0 // debug visualization
	cube_t dbg(room);
	set_cube_zvals(dbg, bcube.z2(), bcube.z2()+bcube.dz());
	objs.emplace_back(dbg, TYPE_DBG_SHAPE, room_id, 0, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, RED);
#endif
	float const clearance(get_min_front_clearance_inc_people()), player_height(get_player_height());
	cube_t place_area(room);
	place_area.expand_by_xy(-support_width); // inside the supports
	place_area.intersect_with_cube(interior->factory_info->floor_space); // clip off side rooms and floor/ceiling
	unsigned const objs_start(objs.size());

	// add ladders to walls
	for (cube_t const &r : nested_rooms) {
		float const wall_pos(r.d[edim][!edir] + 0.5*(edir ? -1.0 : 1.0)*wall_thick); // partially inside the wall
		float const lo(max(r.d[!edim][0], place_area.d[!edim][0])), hi(min(r.d[!edim][1], place_area.d[!edim][1]));
		float const ladder_hwidth(rgen.rand_uniform(0.1, 0.11)*window_vspace);
		if ((hi - lo) < 4.0*ladder_hwidth) continue; // too narrow
		cube_t ladder;
		set_cube_zvals(ladder, zval, (r.z2() + fc_thick + player_height + CAMERA_RADIUS));
		ladder.d[edim][ edir] = wall_pos;
		ladder.d[edim][!edir] = wall_pos + (edir ? -1.0 : 1.0)*0.06*window_vspace; // set depth

		for (unsigned n = 0; n < 10; ++n) { // 10 attempts to place a ladder that doesn't block the door
			set_wall_width(ladder, rgen.rand_uniform((lo + ladder_hwidth), (hi - ladder_hwidth)), ladder_hwidth, !edim);
			if (cube_int_ext_door(ladder) || interior->is_blocked_by_stairs_or_elevator(ladder)) continue; // blocked
			objs.emplace_back(ladder, TYPE_INT_LADDER, room_id, edim, !edir, RO_FLAG_IN_FACTORY, light_amt);
			cube_t blocker(ladder);
			blocker.d[edim][!edir] += (edir ? -1.0 : 1.0)*clearance;
			objs.emplace_back(blocker, TYPE_BLOCKER, room_id, edim, !edir, RO_FLAG_INVIS);
			break;
		} // for n
	} // for r
	// add machines
	add_machines_to_factory(rgen, room, place_area, zval, room_id, light_amt, objs_start);
	// TODO: catwalks
	// TODO: large fans in the ceiling
	// TODO: fire extinguishers, stacks of boxes and crates, paint cans, buckets, breaker panels, fire sprinklers, clock, transformer, water fountain?
}