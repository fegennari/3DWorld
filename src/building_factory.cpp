// 3D World - Factory Buildings
// by Frank Gennari 1/26/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

float const SMOKE_VEL_FACTORY = 0.0008;

extern float CAMERA_RADIUS, fticks;
extern object_model_loader_t building_obj_model_loader;


float shift_val_to_not_intersect_window(cube_t const &c, float val, float hspace, float window_border, bool dim);
cube_t place_cylin_object(rand_gen_t rgen, cube_t const &place_on, float radius, float height, float dist_from_edge, bool place_at_z1=0);

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
		add_room(room_inner, part_id, (is_larger ? 2 : 1)); // 2 lights in larger room; not a typical office building office
		rooms.back().assign_all_to(is_larger ? RTYPE_OFFICE : RTYPE_BATH); // office or bathroom
		interior->factory_info->sub_rooms.push_back(room_inner);
	} // for d
	// should there be an entryway room, then the factory doesn't overlap the sub-rooms? but then there will be empty space above them
	// add entire part as a room (factory floor); must be done last so that smaller contained rooms are picked up in early exit queries (model occlusion, light toggle, door conn)
	add_room(part, part_id); // num_lights will be calculated later
	rooms.back().assign_all_to(RTYPE_FACTORY);
	rooms.back().is_single_floor = 1;
}

void add_ladder_with_blocker(cube_t const &ladder, unsigned room_id, bool dim, bool dir, float light_amt, float clearance, vect_room_object_t &objs, unsigned flags=0) {
	objs.emplace_back(ladder, TYPE_INT_LADDER, room_id, dim, dir, (flags | RO_FLAG_IN_FACTORY), light_amt);
	cube_t blocker(ladder);
	blocker.d[dim][dir] += (dir ? 1.0 : -1.0)*clearance;
	objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, dir, RO_FLAG_INVIS);
}
void add_hanging_ladder(cube_t const &ladder, unsigned room_id, bool dim, bool dir, float light_amt, float clearance, float collider_z2, vect_room_object_t &objs) {
	add_ladder_with_blocker(ladder, room_id, dim, dir, light_amt, clearance, objs, RO_FLAG_HANGING);
	// add an invisible collider to the back of the ladder to block the player
	float const edge_pos(ladder.d[dim][!dir]);
	cube_t collider(ladder);
	collider.z2() = collider_z2;
	collider.d[dim][ dir] = edge_pos;
	collider.d[dim][!dir] = edge_pos - (dir ? 1.0 : -1.0)*0.1*ladder.get_sz_dim(dim); // set depth
	objs.emplace_back(collider, TYPE_COLLIDER, room_id, dim, !dir, (RO_FLAG_INVIS | RO_FLAG_ADJ_TOP), light_amt);
}

void building_t::add_factory_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start_inc_lights) {
	assert(interior->factory_info);
	bool const edim(interior->factory_info->entrance_dim), edir(interior->factory_info->entrance_dir), beam_dim(!edim); // edim is the long dim; beam_dim is short dim
	float const window_vspace(get_window_vspace()), wall_thick(get_wall_thickness()), fc_thick(get_fc_thickness()), edir_sign(edir ? 1.0 : -1.0);
	vect_room_object_t &objs(interior->room_geom->objs);
	// add support pillars around the exterior, between windows; add ceiling beams
	float const support_width(FACTORY_BEAM_THICK*wall_thick), support_hwidth(0.5*support_width);
	float const ceil_zval(room.z2() - fc_thick), beams_z1(ceil_zval - support_width), room_center_short(room.get_center_dim(!edim));
	vector2d const room_sz(room.get_size_xy());
	cube_t support_bounds(room), support, beam;
	support_bounds.expand_by_xy(-support_hwidth);
	set_cube_zvals(support, zval,     beams_z1 );
	set_cube_zvals(beam,    beams_z1, ceil_zval);
	vect_cube_t support_parts, beams, supports, lights;
	vector<float> beam_pos; // in short dim; needed for hanging catwalks
	vect_cube_t const &nested_rooms(interior->factory_info->sub_rooms);
	unsigned const objs_start_inc_beams(objs.size());
	float const shift_vals[6] = {-0.1, 0.2, -0.3, 0.4, -0.5, 0.6}; // cumulative version of {-0.1, 0.1, -0.2, 0.2, -0.3, 0.3}; not enough shift to overlap a window
	unsigned support_count(0);
	// gather all lights and calculate light amount
	float light_amt(0.0);

	for (auto i = objs.begin()+objs_start_inc_lights; i != objs.end(); ++i) {
		if (i->type != TYPE_LIGHT) continue;
		lights.push_back(*i);
		if (i->is_lit()) {light_amt += 1.0;}
	}
	assert(!lights.empty());
	light_amt /= lights.size(); // average

	for (unsigned dim = 0; dim < 2; ++dim) {
		bool const short_dim(bool(dim) == beam_dim);
		unsigned const num_windows(get_num_windows_on_side(room, !dim));
		float const spacing(room_sz[!dim]/num_windows);

		for (unsigned dir = 0; dir < 2; ++dir) { // walls
			float const wall_pos(room.d[dim][dir]);
			support.d[dim][ dir] = wall_pos;
			support.d[dim][!dir] = wall_pos + (dir ? -1.0 : 1.0)*support_width;

			for (unsigned n = 0; n <= num_windows; ++n) {
				if (!short_dim && (n == 0 || n == num_windows)) continue; // don't duplicate add corners
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
				if ((bool)dim == edim) {supports.push_back(support);} // currently only need to track in one dim for catwalk + ladder placement
				bool was_clipped(0);

				for (cube_t const &r : nested_rooms) { // clip in Z if intersects a room
					if (!r.intersects(support)) continue;
					cube_t bot(support);
					objs.back().z1() = r.z2() + fc_thick; // top of room cube
					bot.z2() = objs.back().z1() + 0.1*fc_thick; // shift up slightly
					cube_t sub(r);
					sub.expand_by(0.5*wall_thick);
					subtract_cube_from_cube(bot, sub, support_parts, 1); // clear_out=1
					for (cube_t const &c : support_parts) {objs.emplace_back(c, TYPE_IBEAM, room_id, dim, 1, RO_FLAG_ADJ_TOP, light_amt, SHAPE_CUBE, WHITE, skip_faces);}
					was_clipped = 1;
					break;
				} // for r
				if (!was_clipped && (support_count % 4) == 0) { // add a roof drain pipe every 4th support, if unclipped
					bool const side(centerline < room.get_center_dim(!dim)); // side of beam to place the pipe; closer to building center
					float const pipe_radius(0.6*support_hwidth);
					point pipe_center(0.0, 0.0, zval);
					pipe_center[ dim] = wall_pos + (dir ? -1.0 : 1.0)*pipe_radius; // against the wall
					pipe_center[!dim] = support.d[!dim][side] + (side ? 1.0 : -1.0)*1.25*pipe_radius;
					cube_t pipe(pipe_center);
					pipe.expand_by_xy(pipe_radius);
					pipe.z2() = ceil_zval; // all the way up to the ceiling (beams_z1?)

					if (!has_bcube_int(pipe, nested_rooms) && !cube_int_ext_door(pipe) && !interior->is_blocked_by_stairs_or_elevator(pipe)) {
						objs.emplace_back(pipe, TYPE_PIPE, room_id, 0, 1, (RO_FLAG_NOCOLL | RO_FLAG_LIT), light_amt, SHAPE_CYLIN, DK_BROWN); // vertical, shadowed
					}
				}
				++support_count;
			} // for n
		} // for dir
		// add horizontal beams in the roof
		unsigned const num_hdiv(2*num_windows); // add intermediate beams and hang lights on them
		for (unsigned d = 0; d < 2; ++d) {beam.d[dim][d] = room.d[dim][d];}
		if (short_dim) {beam.expand_in_dim(beam_dim, -support_hwidth);} // half overlap of vert supports

		for (unsigned n = 0; n <= num_hdiv; ++n) {
			if (!short_dim && n > 0 && n < num_hdiv) continue; // only add edge beams in this dim
			float const centerline(max(support_bounds.d[!dim][0], min(support_bounds.d[!dim][1], room.d[!dim][0] + n*0.5f*spacing))); // clamp to support_bounds
			set_wall_width(beam, centerline, support_hwidth, !dim);
			unsigned skip_faces(EF_Z2);
			if (n == 0       ) {skip_faces |= ~get_face_mask(!dim, 0);}
			if (n == num_hdiv) {skip_faces |= ~get_face_mask(!dim, 1);}
			objs.emplace_back(beam, TYPE_IBEAM, room_id, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, WHITE, skip_faces);
			beams.push_back(beam); // needed for sprinkler pipe hanging brackets
			if (short_dim) {beam_pos.push_back(centerline);}
		} // for n
		if (!short_dim) { // add center beam; will intersect/overlap beams in other dim
			set_wall_width(beam, room_center_short, support_hwidth, !dim);
			objs.emplace_back(beam, TYPE_IBEAM, room_id, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, WHITE, EF_Z2);
			beams.push_back(beam); // needed for sprinkler pipe hanging brackets
		}
	} // for dim
#if 0 // debug visualization
	cube_t dbg(room);
	set_cube_zvals(dbg, bcube.z2(), bcube.z2()+bcube.dz());
	objs.emplace_back(dbg, TYPE_DBG_SHAPE, room_id, 0, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, RED);
#endif
	float const clearance(/*get_min_front_clearance_inc_people()*/get_min_front_clearance());
	float const player_height(get_player_height()), ladder_extend_up(player_height + CAMERA_RADIUS);
	cube_t place_area(room);
	place_area.expand_by_xy(-support_width); // inside the supports
	cube_t const place_area_upper(place_area); // includes space above office and bathroom
	place_area.intersect_with_cube(interior->factory_info->floor_space); // clip off side rooms and floor/ceiling
	cube_t const &entry(interior->factory_info->entrance_area);
	unsigned const objs_start(objs.size());

	// add ladders to walls
	for (cube_t const &r : nested_rooms) {
		float const wall_pos(r.d[edim][!edir] - 0.5*edir_sign*wall_thick); // partially inside the wall
		float const lo(max(r.d[!edim][0], place_area.d[!edim][0])), hi(min(r.d[!edim][1], place_area.d[!edim][1]));
		float const ladder_hwidth(rgen.rand_uniform(0.1, 0.11)*window_vspace), edge_spacing(2.0*ladder_hwidth);
		if ((hi - lo) < 4.0*edge_spacing) continue; // too narrow
		cube_t ladder;
		set_cube_zvals(ladder, zval, (r.z2() + fc_thick + ladder_extend_up));
		ladder.d[edim][ edir] = wall_pos;
		ladder.d[edim][!edir] = wall_pos - edir_sign*0.06*window_vspace; // set depth

		for (unsigned n = 0; n < 10; ++n) { // 10 attempts to place a ladder that doesn't block the door
			set_wall_width(ladder, rgen.rand_uniform((lo + edge_spacing), (hi - edge_spacing)), ladder_hwidth, !edim);
			if (cube_int_ext_door(ladder) || interior->is_blocked_by_stairs_or_elevator(ladder)) continue; // blocked
			add_ladder_with_blocker(ladder, room_id, edim, !edir, light_amt, clearance, objs);
			break; // success
		}
	} // for r
	// add catwalk above the entryway connecting the roofs of the two sub-rooms
	float const doorway_width(get_doorway_width()), catwalk_width(1.1*doorway_width), catwalk_hwidth(0.5*catwalk_width), catwalk_height(0.5*window_vspace);
	float const cw_lo(entry.d[edim][0] + 1.5*catwalk_hwidth), cw_hi(entry.d[edim][1] - 1.2*catwalk_hwidth);
	float const trim_thickness(get_trim_thickness());

	if (cw_lo < cw_hi) { // should always be true
		cube_t catwalk(entry);
		catwalk.z1() = entry  .z2() + fc_thick;
		catwalk.z2() = catwalk.z1() + catwalk_height;
		set_wall_width(catwalk, rgen.rand_uniform(cw_lo, cw_hi), catwalk_hwidth, edim);
		objs.emplace_back(catwalk, TYPE_CATWALK, room_id, !edim, rgen.rand_bool(), RO_FLAG_IN_FACTORY, light_amt); // random mesh texture
	}
	if (1) { // add central upper catwalk in long dim
		// find a location near the center not obstructed by a door, stairs, or beam
		float const max_shift(0.3*room_sz[!edim]), wall_pos(room.d[edim][!edir]); // wall opposite the main entrance; don't shift so much that we intersect ceiling fans
		float catwalk_center(room_center_short), cur_shift(1.0*support_width*(rgen.rand_bool() ? 1.0 : -1.0));
		bool add_catwalk(0);
		cube_t cand(room); // copy zvals
		//cand.translate_dim(edim, -edir_sign*doorway_width);
		set_wall_width(cand, wall_pos, 2.0*clearance, edim); // make sure it overlaps exterior doors only on one side

		while (fabs(cur_shift) < max_shift) {
			set_wall_width(cand, catwalk_center, catwalk_hwidth, !edim);
			
			if (!has_bcube_int(cand, supports) && !cube_int_ext_door(cand) && !interior->is_blocked_by_stairs_or_elevator(cand)) { // valid placement
				bool bad_place(0); // don't let supports intersect lights
				
				for (unsigned d = 0; d < 2; ++d) { // each side of catwalk
					cube_t side(cand);
					side.d[!edim][!d] = side.d[!edim][d]; // shrink to zero area
					if (has_bcube_int(side, lights)) {bad_place = 1; break;}
				}
				if (!bad_place) {add_catwalk = 1; break;} // success
			}
			catwalk_center = room_center_short + cur_shift;
			cur_shift      = -1.2*cur_shift; // bad placement; increase step and swap sides
		} // while
		if (add_catwalk) {
			// determine catwalk placement
			float const ladder_depth(0.1*window_vspace), end_pad(clearance - ladder_depth - support_width);
			cube_t catwalk(place_area_upper); // set the length
			catwalk.z1() = room   .z2() - window_vspace + fc_thick - support_width; // would be upper floor zval - support_width
			catwalk.z2() = catwalk.z1() + catwalk_height;
			set_wall_width(catwalk, catwalk_center, catwalk_hwidth, !edim);
			// connect to floor with ladder at one end
			float const ladder_hwidth(0.5*catwalk_hwidth);
			cube_t ladder;
			set_cube_zvals(ladder, zval, (catwalk.z1() + ladder_extend_up)); // extend up to catwalk
			set_wall_width(ladder, catwalk_center, ladder_hwidth, !edim);
			ladder.d[edim][!edir] = wall_pos + edir_sign*trim_thickness; // prevent Z-fighting
			ladder.d[edim][ edir] = wall_pos + edir_sign*ladder_depth; // set depth
			add_ladder_with_blocker(ladder, room_id, edim, edir, light_amt, clearance, objs);
			objs.back().expand_in_dim(!edim, (catwalk_hwidth - ladder_hwidth)); // expand blocker to the width of the ladder; required to avoid nearby pipes, etc.
			// try to add a ladder on the other end; it may be on top of a sub-room
			cube_t ladder2(ladder);
			ladder2.translate_dim(edim, edir_sign*(room_sz[edim] - ladder_depth - 2.0*trim_thickness)); // translate to the other wall
			bool add_sec_ladder(1);

			for (cube_t const &r : nested_rooms) {
				if (r.contains_cube_xy(ladder2)) {ladder2.z1()   = r.z2() + fc_thick; break;} // ladd starts on the roof of this room
				if (r.intersects_xy   (ladder2)) {add_sec_ladder = 0; break;} // skip if partially overlapping
			}
			if (add_sec_ladder && has_bcube_int(ladder2, supports) || cube_int_ext_door(ladder2) || interior->is_blocked_by_stairs_or_elevator(ladder2)) {add_sec_ladder = 0;}
			if (add_sec_ladder) {add_ladder_with_blocker(ladder2, room_id, edim, !edir, light_amt, clearance, objs);}
			// create catwalk
			catwalk.expand_in_dim(edim, -end_pad); // shrink
			unsigned const flags(RO_FLAG_IN_FACTORY | RO_FLAG_HANGING);
			room_object_t catwalk_obj(catwalk, TYPE_CATWALK, room_id, edim, rgen.rand_bool(), flags, light_amt); // random mesh texture
			// add bars to hang the catwalk from beams
			float const rod_radius(0.008*window_vspace);
			bool last_added(0);

			for (float bpos : beam_pos) {
				if (bpos < catwalk.d[edim][0] || bpos > catwalk.d[edim][1]) continue; // not over catwalk
				if (last_added) {last_added = 0; continue;} // alternate every other beam
				last_added = 1;

				for (unsigned d = 0; d < 2; ++d) { // each side of catwalk
					point center;
					center[ edim] = bpos;
					center[!edim] = catwalk.d[!edim][d] + (d ? 1.0 : -1.0)*0.9*rod_radius; // just outside the edge and slightly overlapping
					cube_t rod(center);
					rod.expand_by_xy(rod_radius);
					set_cube_zvals(rod, (catwalk.z1() + 1.0*rod_radius), beams_z1);
					unsigned const pflags(RO_FLAG_IN_FACTORY | RO_FLAG_HANGING | RO_FLAG_ADJ_LO | RO_FLAG_NOCOLL | RO_FLAG_LIT); // draw bottom surface
					objs.emplace_back(rod, TYPE_PIPE, room_id, 0, 1, pflags, light_amt, SHAPE_CYLIN); // vertical
				} // for d
			} // for bpos
			// try adding a ladder to the middle of the catwalk
			bool const ldir(rgen.rand_bool());
			float const edge_pos(catwalk.d[!edim][ldir]), ldsign(ldir ? 1.0 : -1.0);
			float const llo(catwalk.d[edim][0] + catwalk_hwidth), lhi(catwalk.d[edim][1] - catwalk_hwidth);
			set_wall_width(ladder, edge_pos, 0.5*ladder_depth, !edim); // centered on the catwalk edge; if adjacent, the player will fall

			for (unsigned n = 0; n < 20; ++n) { // some random attempts
				if (llo >= lhi) continue; // shouldn't happen
				set_wall_width(ladder, rgen.rand_uniform(llo, lhi), ladder_hwidth, edim);
				cube_t tc(ladder);
				tc.d[!edim][ldir] += ldsign*clearance; // add space in front
				if (overlaps_obj_or_placement_blocked(tc, room, objs_start) || has_bcube_int(tc, nested_rooms)) continue;
				// ladder pos is valid; split the catwalk that was added into left, right, and center parts so that we can insert a cut into the railing
				cube_t seg_split(ladder);
				seg_split.expand_in_dim(edim, 0.25*doorway_width); // increase width beyond the ladder
				room_object_t &far_end(catwalk_obj); // opposite the end ladder
				room_object_t near_end(catwalk_obj), middle(catwalk_obj);
				middle.flags |= (ldir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO); // disable side bars on ladder side for middle catwalk segment
				far_end .d[edim][!edir] = middle.d[edim][ edir] = seg_split.d[edim][ edir];
				near_end.d[edim][ edir] = middle.d[edim][!edir] = seg_split.d[edim][!edir];
				objs.push_back(near_end);
				objs.push_back(middle  );
				add_hanging_ladder(ladder, room_id, !edim, ldir, light_amt, clearance, catwalk.z1(), objs);
				break; // success/done
			} // for n
			if (!add_sec_ladder) {catwalk_obj.flags |= (edir ? RO_FLAG_ADJ_TOP : RO_FLAG_ADJ_BOT);} // add end bars to non-ladder side unless there's a second ladder
			objs.push_back(catwalk_obj);
			// TODO: connect with stairs along a wall?
		}
	} // end upper catwalk
	// add transformer
	float const tzval(zval - 0.02*window_vspace); // transformer is slightly below floor level
	cube_t xfmr_area(place_area);
	xfmr_area.d[edim][edir] -= edir_sign*1.2*doorway_width; // don't place too close to sub-rooms or entrance
	place_model_along_wall(OBJ_MODEL_SUBSTATION, TYPE_XFORMER, room, 0.6, rgen, tzval, room_id, light_amt, xfmr_area, objs_start, 0.0, 4, 0, WHITE, 0, 0, 0, 1); // sideways
	// add machines
	unsigned const machines_start(objs.size());
	add_machines_to_factory(rgen, room, place_area, zval, room_id, light_amt, objs_start, objs_start_inc_beams);

	for (auto i = objs.begin()+machines_start; i != objs.end(); ++i) { // add some smoke emitters
		if (i->type != TYPE_MACHINE)  continue;
		if (rgen.rand_float() > 0.15) continue; // 15% chance
		float const smoke_radius(0.12*(i->dx() + i->dy()));
		point const pos(i->xc(), i->yc(), (i->z1() + 0.8*i->dz())); // smoke starts near machine top and rises up to the ceiling
		interior->factory_info->smoke_emitters.emplace_back(pos, smoke_radius, (i - objs.begin())); // store machine obj id
	}
	// add fire extinguisher
	add_fire_ext_along_wall(entry, zval, room_id, light_amt, !edim, rgen.rand_bool(), rgen); // choose a random side
	
	if (!nested_rooms.empty()) { // add breaker panel
		float const bp_ceil_zval(zval + get_floor_ceil_gap());

		for (unsigned n = 0; n < 10; ++n) { // 10 attempts
			unsigned const rix(rgen.rand() % nested_rooms.size());
			cube_t const &bpr(nested_rooms[rix]);
			bool const bpdim(rgen.rand_bool()), bpdir((bpdim == edim) ? !edir : (bpr.get_center_dim(bpdim) < entry.get_center_dim(bpdim)));
			float const hwidth(0.5*rgen.rand_uniform(0.25, 0.35)*window_vspace), depth(0.04*window_vspace), edge_space(hwidth + support_width);
			float const lo(bpr.d[!bpdim][0] + edge_space), hi(bpr.d[!bpdim][1] - edge_space);
			if (lo >= hi) continue; // wall too short
			float const dsign(bpdir ? 1.0 : -1.0), wall_pos(bpr.d[bpdim][bpdir] + dsign*wall_thick);
			cube_t breaker_panel;
			set_cube_zvals(breaker_panel, (bp_ceil_zval - 0.7*window_vspace), (bp_ceil_zval - rgen.rand_uniform(0.25, 0.3)*window_vspace));
			set_wall_width(breaker_panel, rgen.rand_uniform(lo, hi), hwidth, !bpdim);
			breaker_panel.d[bpdim][!bpdir] = wall_pos;
			breaker_panel.d[bpdim][ bpdir] = wall_pos + dsign*depth;
			cube_t tc(breaker_panel);
			tc.d[bpdim][bpdir] += dsign*2.0*hwidth; // add clearance so that it can open
			if (is_cube_close_to_doorway(tc, room, 0.0, 1) || overlaps_other_room_obj(tc, objs_start)) continue; // avoid doors and machines
			add_breaker_panel(rgen, breaker_panel, bp_ceil_zval, bpdim, !bpdir, room_id, light_amt);
			objs.emplace_back(tc, TYPE_BLOCKER, room_id, bpdim, !bpdir, RO_FLAG_INVIS); // don't place objects such as boxes in front of the panel
			break; // success
		} // for n
	}
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_WFOUNTAIN)) { // add water fountain
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_WFOUNTAIN)); // D, W, H
		float const height(0.25*window_vspace), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z), edge_space(hwidth + support_width), z1(zval + 0.18*window_vspace);

		for (unsigned n = 0; n < 10; ++n) { // 10 attempts
			unsigned const rix(rgen.rand() % nested_rooms.size());
			cube_t const &wfr(nested_rooms[rix]);
			bool const wfdim(rgen.rand_bool()), wfdir((wfdim == edim) ? !edir : (wfr.get_center_dim(wfdim) < entry.get_center_dim(wfdim)));
			float const lo(wfr.d[!wfdim][0] + edge_space), hi(wfr.d[!wfdim][1] - edge_space);
			if (lo >= hi) continue; // wall too short
			float const wall_pos(wfr.d[wfdim][wfdir] + (wfdir ? 1.0 : -1.0)*wall_thick);
			cube_t wf;
			set_wall_width(wf, rgen.rand_uniform(lo, hi), hwidth, !wfdim);
			set_cube_zvals(wf, z1, z1+height);
			wf.d[wfdim][!wfdir] = wall_pos;
			wf.d[wfdim][ wfdir] = wall_pos + (wfdir ? 1.0 : -1.0)*depth;
			if (is_cube_close_to_doorway(wf, room, 0.0, 1) || overlaps_other_room_obj(wf, objs_start)) continue; // avoid doors and machines
			objs.emplace_back(wf, TYPE_WFOUNTAIN, room_id, wfdim, !wfdir, 0, light_amt, SHAPE_CUBE);
			break; // success
		} // for n
	}
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_CEIL_FAN)) { // add ceiling fans
		unsigned const fans_per_beam = 2;
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_CEIL_FAN)); // D, W, H
		float const diameter(0.8*window_vspace), height(diameter*sz.z/sz.y); // assumes width = depth = diameter
		float const length(room_sz[!edim]), spacing(length/fans_per_beam);
		point top_center(0.0, 0.0, beams_z1);
		bool last_added(0);

		for (float bpos : beam_pos) {
			if (last_added) {last_added = 0; continue;} // alternate every other beam
			last_added = 1;
			top_center[edim] = bpos;

			for (unsigned n = 0; n < fans_per_beam; ++n) {
				top_center[!edim] = room.d[!edim][0] + (n + 0.5)*spacing;
				cube_t fan(top_center);
				fan.expand_by_xy(0.5*diameter);
				fan.z1() -= height;
				if (!support_bounds.contains_cube(fan))          continue; // skip for edge beams
				if (overlaps_other_room_obj(fan, objs_start, 1)) continue; // including other fans; check_all=1
				unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_IN_FACTORY);
				if (rgen.rand_float() < 0.75) {flags |= RO_FLAG_ROTATING;} // make fan rotate when turned on 75% of the time
				objs.emplace_back(fan, TYPE_CEIL_FAN, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN, WHITE);
			} // for n
		} // for bpos
	}
	// add fire sprinkler pipes
	float const custom_floor_spacing(room.dz() - support_width); // place sprinklers under ceiling beams
	float const wall_pad(1.05*support_width); // add a gap between the wall supports
	vect_cube_t obstacles(nested_rooms), walls, pipe_cubes; // avoid nested rooms; walls and pipe_cubes start empty and are unused
	cube_t entry_area(entry);
	entry_area.z2() = room.z1() + window_vspace;
	entry_area.expand_in_dim(edim, doorway_width);
	obstacles.push_back(entry_area); // don't block the entryway with the vertical pipe

	for (tquad_with_ix_t const &door : doors) {
		obstacles.push_back(door.get_bcube());
		obstacles.back().expand_by_xy(doorway_width); // don't block an exterior door
	}
	for (stairwell_t const &s : interior->stairwells) { // is this needed?
		if (s.intersects(room)) {obstacles.push_back(s);}
	}
	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (!i->no_coll()) {obstacles.push_back(*i);}
	}
	// main pipe not under lights? but lights are placed later
	unsigned const main_pipe_ix(objs.size());
	bool const added_sprinklers(add_sprinkler_pipes(obstacles, walls, beams, pipe_cubes, room_id, 1, objs_start, rgen, custom_floor_spacing, wall_pad)); // num_floors=1

	// add ceiling ducts and vents (similar to malls)
	float const ducts_z2(room.z2() - 0.8*window_vspace); // under the first window, below beams, lights, and pipes
	cube_t duct_bounds(room);
	duct_bounds.expand_by_xy(-(support_width - 0.5*wall_thick));

	if (added_sprinklers) { // shorten ends so as to not overlap sprinkler pipe
		room_object_t const &main_pipe(objs[main_pipe_ix]);
		assert(main_pipe.type == TYPE_PIPE);
		duct_bounds.expand_in_dim(edim, -1.5*main_pipe.get_sz_dim(edim));
	}
	bool const cylin_ducts(rgen.rand_bool());
	unsigned const ducts_start(objs.size());
	add_ceiling_ducts(duct_bounds, ducts_z2, room_id, edim, 2, light_amt, cylin_ducts, 0, 0, rgen); // draw ends and top
	unsigned const ducts_end(objs.size());

	// add vertical section of duct up to the ceiling; maybe should add HVAC unit along the vent somewhere?
	for (unsigned i = ducts_start; i < ducts_end; ++i) {
		room_object_t const &obj(objs[i]);
		if (obj.type != TYPE_DUCT || obj.dim != edim) continue; // not the long duct
		// if there are an odd number of vertical supports, offset to a random side to avoid intersecting the center support
		float const hwidth(0.5*obj.dz()), center(obj.get_center_dim(edim) + (rgen.rand_bool() ? 1.0 : -1.0)*(support_hwidth + hwidth));
		cube_t duct(obj);
		set_cube_zvals(duct, obj.zc(), ceil_zval);
		set_wall_width(duct, center, hwidth, edim);
		unsigned const duct_flags(RO_FLAG_IN_MALL | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI);
		objs.emplace_back(duct, TYPE_DUCT, room_id, 0, 1, duct_flags, light_amt, obj.shape); // vertical, same shape
	} // for i
	// add boxes and crates in piles
	unsigned const num_piles(4 + (rgen.rand() % 5)); // 4-8
	vect_cube_t exclude; // empty

	for (unsigned n = 0; n < num_piles; ++n) {
		unsigned const side(rgen.rand() % 3); // all sides but the entrance
		bool const pdim((side == 2) ? edim : !edim), pdir((pdim == edim) ? !edir : bool(side));
		point center;
		center[ pdim] = place_area.d[pdim][pdir]; // against the wall/pillars
		center[!pdim] = rgen.rand_uniform(place_area.d[!pdim][0], place_area.d[!pdim][1]); // random pos along wall
		cube_t crate_bounds(center);
		crate_bounds.expand_in_dim( pdim, window_vspace*rgen.rand_uniform(0.5, 1.0)); // random depth; half the space is outside the bounds and wasted
		crate_bounds.expand_in_dim(!pdim, window_vspace*rgen.rand_uniform(0.5, 1.0)); // random width; may be clipped at the corner of the building
		crate_bounds.intersect_with_cube(place_area);
		if (interior->is_blocked_by_stairs_or_elevator(crate_bounds)) continue; // don't block stairs
		add_boxes_and_crates(rgen, room, zval, room_id, light_amt, objs_start, 13, 0, interior->factory_info->floor_space, crate_bounds, exclude); // 4-16; is_basement=0
	} // for n
	// add buckets
	unsigned const num_buckets((rgen.rand() % 4) + 1); // 1-4
	add_buckets_to_room(rgen, place_area, zval, room_id, light_amt, objs_start, num_buckets);
	// add paint cans
	unsigned const num_pcans((rgen.rand() % 4) + 1); // 1-4

	for (unsigned n = 0; n < num_pcans; ++n) {
		float const height(rgen.rand_uniform(0.08, 0.1)*window_vspace), radius(0.44*height);
		cube_t const pcan(place_cylin_object(rgen, place_area, radius, height, (radius + trim_thickness), 1)); // place_at_z1=1
		if (overlaps_obj_or_placement_blocked(pcan, place_area, objs_start)) continue; // bad placement
		objs.emplace_back(pcan, TYPE_PAINTCAN, room_id, 0, 0, 0, light_amt, SHAPE_CYLIN, WHITE);
	}
	// add floor clutter and stains
	bool const add_bottles(1), add_papers(0), add_glass(1), add_trash(rgen.rand_float() < 0.65); // 65% of rooms
	place_area.z1() = zval; // is this needed/correct?
	add_floor_clutter_objs(rgen, room, place_area, zval, room_id, light_amt, objs_start, add_bottles, add_trash, add_papers, add_glass);
	unsigned const num_floor_stains(rgen.rand() % 9); // 0-8
	float const stain_rmax(0.25*min(window_vspace, min(room.dx(), room.dy())));
	add_floor_stains(rgen, place_area, zval, room_id, light_amt, objs_start, num_floor_stains, stain_rmax);
}

void building_t::add_factory_office_objs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, unsigned floor, float tot_light_amt, unsigned objs_start) {
	// add a clock on the wall
	bool const digital(rgen.rand_bool());
	float const floor_spacing(get_window_vspace()), clock_height((digital ? 0.08 : 0.16)*floor_spacing), clock_z1(zval + get_floor_ceil_gap() - 1.4*clock_height);
	float const clock_width((digital ? 4.0 : 1.0)*clock_height), clock_depth((digital ? 0.05 : 0.08)*clock_width);
	cube_t const place_area(get_walkable_room_bounds(room));
	cube_t const &factory_area(get_factory_area());
	cube_t clock;
	set_cube_zvals(clock, clock_z1, clock_z1+clock_height);

	for (unsigned n = 0; n < 10; ++n) { // 10 attempts
		bool const dim(rgen.rand_bool()), dir(room.d[dim][0] == factory_area.d[dim][0]); // interior wall
		float const edge_space(max(1.5*clock_width, 0.25*room.get_sz_dim(!dim))); // somewhat centered
		float const lo(place_area.d[!dim][0] + edge_space), hi(place_area.d[!dim][1] - edge_space);
		if (lo >= hi) continue; // wall too short
		float const dsign(dir ? -1.0 : 1.0), wall_pos(place_area.d[dim][dir]);
		set_wall_width(clock, rgen.rand_uniform(lo, hi), 0.5*clock_width, !dim);
		clock.d[dim][ dir] = wall_pos;
		clock.d[dim][!dir] = wall_pos + dsign*clock_depth;
		cube_t tc(clock);
		tc.d[dim][!dir] += dsign*0.25*floor_spacing; // add clearance
		if (is_cube_close_to_doorway(tc, room, 0.0, 1) || overlaps_other_room_obj(tc, objs_start)) continue;
		add_clock(clock, room_id, tot_light_amt, dim, !dir, digital);
		break; // success
	} // for n
	// add boxes and crates
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_pre(objs.size());
	add_boxes_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, 8); // up to 8 boxes along walls

	if (objs_pre == objs.size() || rgen.rand_bool()) { // no boxes were added; try adding to the middle of the room; also try half the time
		cube_t const room_bounds(get_walkable_room_bounds(room)), crate_bounds(room_bounds);
		vect_cube_t exclude;
	
		for (door_stack_t const &ds : interior->door_stacks) {
			if (ds.is_connected_to_room(room_id)) {exclude.push_back(ds.get_open_door_bcube_for_room(room));}
		}
		add_boxes_and_crates(rgen, room, zval, room_id, tot_light_amt, objs_start, 3, 0, room_bounds, crate_bounds, exclude); // 4-6; is_basement=0
	}
}

void building_t::add_factory_smokestack(rand_gen_t &rgen) {
	assert(is_factory());
	assert(!parts.empty());
	cube_t const &base(parts[0]);
	float const ss_radius(rgen.rand_uniform(0.01, 0.02)*(base.dx() + base.dy()));
	unsigned const num(1 + (rgen.rand()&1)); // 1-2

	for (unsigned n = 0; n < num; ++n) {
		for (unsigned m = 0; m < 10; ++m) { // 10 attempts to place smokestack
			point const ss_center(gen_xy_pos_in_area(base, 2.5*ss_radius, rgen, base.z2()));
			cube_t smokestack(ss_center);
			smokestack.expand_by_xy(ss_radius);
			smokestack.z2() += rgen.rand_uniform(0.75, 1.0)*base.dz(); // set height; should be above roof peak
			if (!has_bcube_int(smokestack, details)) {details.emplace_back(smokestack, (uint8_t)ROOF_OBJ_SMOKESTACK);}
			break; // success/done
		} // for m
	} // for n
}

void bldg_factory_info_t::next_frame(particle_manager_t &particle_manager) {
	for (smoke_source_t &s : smoke_emitters) { // generate smoke
		s.time += fticks;
		if (s.time < s.next_smoke_time) continue;
		particle_manager.add_particle(s.pos, SMOKE_VEL_FACTORY*plus_z, GRAY, s.radius, PART_EFFECT_SMOKE, s.pid, 0.5); // coll_radius=0.5
		s.next_smoke_time = s.time + rgen.rand_uniform(0.5, 0.8)*TICKS_PER_SECOND;
	}
}