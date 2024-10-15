// 3D World - Building Underground/Basement Tunnels
// by Frank Gennari 10/10/2024

#include "function_registry.h"
#include "buildings.h"

extern double tfticks;

float query_min_height(cube_t const &c, float stop_at);
colorRGBA choose_pipe_color(rand_gen_t &rgen);

// *** tunnel_seg_t ***

tunnel_seg_t::tunnel_seg_t(point const &p1, point const &p2, float radius_) : radius(radius_) {
	if      (p1.x == p2.x) {dim = 1;}
	else if (p1.y == p2.y) {dim = 0;}
	else {assert(0);} // no vertical tunnels for now
	p[0] = p1; p[1] = p2;
	assert(p[0][dim] != p[1][dim]);
	if (p[1][dim] < p[0][dim]) {swap(p[0], p[1]);}
	bcube = cube_t(p[0], p[1]);
	bcube.expand_in_z(radius);
	bcube.expand_in_dim(!dim, radius);
	bcube_ext = bcube;
}
void tunnel_seg_t::set_as_room_conn(bool rdir, float wall_gap) {
	room_conn = 1;
	room_dir  = rdir;
	bcube_ext.d[!dim][room_dir] += (room_dir ? 1.0 : -1.0)*wall_gap;
}
cube_t tunnel_seg_t::get_player_walk_area(point const &player_pos, float player_radius) const {
	cube_t walk_area(bcube_ext);
	float const walk_width(0.1*radius), blocked_width(radius - walk_width); // area inside the tunnel near the center where the player can walk
	if (room_conn) {walk_area.d[!dim][!room_dir] += (room_dir ? 1.0 : -1.0)*blocked_width;} // shrink on non-room side
	else {walk_area.expand_in_dim(!dim, -blocked_width);} // shrink on both sides
	// prevent the player from walking off a closed end or through a gate
	if (closed_ends[0]) {walk_area.d[dim][0] += player_radius;}
	if (closed_ends[1]) {walk_area.d[dim][1] -= player_radius;}

	if (has_gate) {
		if (player_pos[dim] < gate_pos) {min_eq(walk_area.d[dim][1], gate_pos-player_radius);} // player at low end, clamp high end
		else                            {max_eq(walk_area.d[dim][0], gate_pos+player_radius);} // player at high end, clamp low end
	}
	return walk_area;
}
cube_t tunnel_seg_t::get_room_conn_block() const {
	cube_t block(bcube_ext);
	block.d[!dim][ room_dir] += (room_dir ? -1.0 : 1.0)*0.05*radius;
	block.d[!dim][!room_dir]  = block.d[!dim][room_dir] + (room_dir ? -1.0 : 1.0)*0.5*radius;
	block.z2() = block.z1() + 0.22*radius; // set the height
	return block;
}

// *** Placement ***

void building_t::get_valid_extb_room_end_doors(room_t const &room, float zval, unsigned room_id, float end_pad_ext, cube_with_ix_t doors[2]) const {
	assert(interior);
	if (has_pool() && (int)room_id == interior->pool.room_ix) return; // no false doors or tunnels in swimming pool room

	if (!room.is_hallway) {
		float const dx(room.dx()), dy(room.dy());
		if (min(dx, dy) > 0.25*max(dx, dy)) return; // not high aspect ratio
	}
	if (room.is_ext_basement_conn()) return; // don't add a door to the end that connects to the adjacent building's extended basement
	bool const dim(room.dx() < room.dy()); // long dim
	float const door_width(get_doorway_width()), wall_thickness(get_wall_thickness()), expand_val(2.0*wall_thickness);
	assert(interior->ext_basement_hallway_room_id >= 0 && (unsigned)interior->ext_basement_hallway_room_id < interior->rooms.size());

	for (unsigned dir = 0; dir < 2; ++dir) { // check both ends
		cube_t query_region(room);
		query_region.d[dim][!dir]  = room.d[dim][dir]; // shrink to the end of the hallway
		query_region.d[dim][ dir] += (dir ? 1.0 : -1.0)*end_pad_ext; // extend away from the room
		query_region.expand_in_dim(dim, 1.5*door_width); // distance from end to nearest door/room in either direction
		query_region.expand_by_xy(expand_val);
		query_region.expand_in_z(-get_fc_thickness()); // subtract floors and ceilings
		if (interior->is_cube_close_to_doorway(query_region, room))   continue; // bad placement/not needed
		if (interior->is_blocked_by_stairs_or_elevator(query_region)) continue; // check stairs
		bool near_other_room(0);

		for (unsigned r = interior->ext_basement_hallway_room_id; r < interior->rooms.size(); ++r) {
			if (r != room_id && interior->rooms[r].intersects(query_region)) {near_other_room = 1; break;}
		}
		if (near_other_room) continue; // too close to another room to the side of or opposite the door
		cube_t door;
		set_cube_zvals(door, zval, (zval + get_floor_ceil_gap())); // extend to the ceiling
		set_wall_width(door, room.get_center_dim(!dim), 0.5*door_width, !dim);
		set_wall_width(door, (room.d[dim][dir] + (dir ? -1.0 : 1.0)*0.5*wall_thickness), 2.0*get_trim_thickness(), dim);
		doors[dir] = cube_with_ix_t(door, (2*dim + dir));
	} // for dir
}

void building_t::add_false_door_to_extb_room_if_needed(room_t const &room, float zval, unsigned room_id) {
	assert(has_room_geom());
	cube_with_ix_t doors[2];
	get_valid_extb_room_end_doors(room, zval, room_id, 0.0, doors); // end_pad_ext=0.0

	for (unsigned d = 0; d < 2; ++d) { // check both door locations
		cube_with_ix_t const &door(doors[d]);
		if (door.is_all_zeros()) continue;
		bool const dim(door.ix >> 1), dir(door.ix & 1);
		unsigned flags(RO_FLAG_NOCOLL);

		if (room.has_tunnel_conn()) { // add an open door; should this have collisions enabled?
			flags |= (RO_FLAG_HAS_EXTRA | RO_FLAG_OPEN); // make this an open vault/blast door
			//cube_t c(door); c.z1() = c.z2(); c.z2() += 1.0; interior->room_geom->objs.emplace_back(c, TYPE_DBG_SHAPE, room_id, dim, dir, 0, 1.0, SHAPE_CUBE, RED); // TESTING
		}
		else { // add a closed door
			if      (is_house)           {flags |= RO_FLAG_IS_HOUSE ;}
			else if ((room_id & 3) == 0) {flags |= RO_FLAG_HAS_EXTRA;} // make this a vault/blast door 25% of the time
		}
		interior->room_geom->objs.emplace_back(door, TYPE_FALSE_DOOR, room_id, dim, dir, flags, 1.0, SHAPE_CUBE, WHITE);
		// add door blocker to avoid placing objects in front of the door; needed even for hallways to avoid placing an overlapping outlet
		cube_t blocker(door);
		blocker.d[dim][!dir] = room.d[dim][dir] + (dir ? -1.0 : 1.0)*get_doorway_width(); // add clearance
		interior->room_geom->objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, dir, RO_FLAG_INVIS);
	} // for d
}

bool building_t::is_tunnel_placement_valid(point const &p1, point const &p2, float radius) const {
	cube_t const tunnel_bc(tunnel_seg_t(p1, p2, radius).bcube);
	if (cube_intersects_basement_or_extb_room(tunnel_bc))             return 0;
	if (query_min_height(tunnel_bc, tunnel_bc.z2()) < tunnel_bc.z2()) return 0; // check for terrain clipping through ceiling
	if (!is_basement_room_not_int_bldg(tunnel_bc))                    return 0;
	return 1;
}
bool building_t::try_place_tunnel_at_extb_hallway_end(room_t &room, unsigned room_id, rand_gen_t &rgen) {
	if (!room.is_hallway && rgen.rand_bool()) return 0; // 50% chance for non-hallways
	float const zval(room.z1() + get_fc_thickness()), end_pad_ext(2.0*get_doorway_width());
	cube_with_ix_t doors[2];
	get_valid_extb_room_end_doors(room, zval, room_id, end_pad_ext, doors);

	for (unsigned d = 0; d < 2; ++d) { // check both door locations
		cube_with_ix_t const &door(doors[d]);
		if (door.is_all_zeros()) continue;
		bool const dim(door.ix >> 1), dir(door.ix & 1);
		float const wall_thickness(get_wall_thickness()), floor_spacing(get_window_vspace()), sm_shift_val(0.5*get_rug_thickness());
		float const min_len(8.0*floor_spacing), max_len(20.0*floor_spacing); // in each direction
		float const radius(0.5*door.dz()), wall_gap(2.0*wall_thickness), check_radius(radius + wall_thickness), dist_from_door(radius + wall_gap);
		cube_t wall_clip(door);
		wall_clip.expand_in_dim(dim, 2.0*wall_thickness); // make sure it contains the wall
		subtract_cube_from_cubes(wall_clip, interior->walls[dim]); // remove door from wall
		point middle(door.get_cube_center());
		middle[dim] += (dir ? 1.0 : -1.0)*dist_from_door;
		middle.z -= sm_shift_val; // shift down slightly to prevent Z-fighting with concrete on ceiling and floor
		point p1(middle), p2(middle);
		p1[!dim] -= min_len; p2[!dim] += min_len; // start at min length in each dim
		if (!is_tunnel_placement_valid(p1, p2, check_radius)) continue; // can't place a tunnel of min length
		unsigned const num_steps = 10;
		float const step_len((max_len - min_len)/num_steps);

		// try extending tunnel in both directions
		for (unsigned d = 0; d < 2; ++d) {
			for (unsigned n = 0; n < num_steps; ++n) {
				if (d) {
					point p2e(p2);
					p2e[!dim] += step_len;
					if (!is_tunnel_placement_valid(p2, p2e, check_radius)) break; // can't extend further in this dir
					p2 = p2e; // accept the new length
				}
				else {
					point p1e(p1);
					p1e[!dim] -= step_len;
					if (!is_tunnel_placement_valid(p1e, p1, check_radius)) break; // can't extend further in this dir
					p1 = p1e; // accept the new length
				}
			}
		} // for dir
		// split into three segments (center, left, right)
		// TODO: add bends at the ends if there's space
		unsigned const tseg_ix(interior->tunnels.size());
		float const gate_dist_from_end(5.0*floor_spacing);
		point pa(p1), pb(p2);
		pa[!dim] = door.d[!dim][0] + sm_shift_val; // left  end of room connection
		pb[!dim] = door.d[!dim][1] - sm_shift_val; // right end of room connection
		tunnel_seg_t tseg_c(pa, pb, radius), tseg_l(p1, pa, radius), tseg_r(pb, p2, radius);
		tseg_c.set_as_room_conn(!dir, wall_gap);
		tseg_l.closed_ends[0] = tseg_r.closed_ends[1] = tseg_l.has_gate = tseg_r.has_gate = 1;
		tseg_l.gate_pos = p1[!dim] + gate_dist_from_end;
		tseg_r.gate_pos = p2[!dim] - gate_dist_from_end;
		tseg_l.conn_ix[1] = tseg_ix;
		tseg_r.conn_ix[0] = tseg_ix;
		for (unsigned d = 0; d < 2; ++d) {tseg_c.conn_ix[d] = tseg_ix + d + 1;}
		tunnel_seg_t *to_add[3] = {&tseg_c, &tseg_l, &tseg_r};
		float const water_level(rgen.rand_uniform(0.0, 1.0)*0.2*radius), water_flow(rgen.signed_rand_float());

		for (unsigned n = 0; n < 3; ++n) {
			tunnel_seg_t &tseg(*to_add[n]);
			tseg.water_level  = water_level;
			tseg.water_flow   = water_flow;
			tseg.conn_room_ix = room_id;
			interior->tunnels.emplace_back(tseg);
		}
		interior->basement_ext_bcube.assign_or_union_with_cube(tunnel_seg_t(p1, p2, radius).bcube); // add full tunnel to bcube
		room.set_has_tunnel_conn();
		return 1; // done (only add tunnel conn to one end)
	} // for d
	return 0;
}

void building_t::add_tunnel_objects(rand_gen_t rgen) {
	assert(has_room_geom());
	vect_room_object_t &objs(interior->room_geom->objs);

	for (tunnel_seg_t const &t : interior->tunnels) {
		if (t.room_conn) continue; // nothing added to these tunnels
		// add smaller pipes
		unsigned const num_pipes(rgen.rand() % 3); // 0-2

		for (unsigned n = 0; n < num_pipes; ++n) {
			float const radius(0.05*t.radius*rgen.rand_uniform(0.5, 1.0));
			float const v1(t.p[0][t.dim] + 2.0*radius), v2(t.p[1][t.dim] - 2.0*radius);
			if (v1 >= v2) continue; // too short a tunnel; shouldn't happen
			float const pos(rgen.rand_uniform(v1, v2)), height(t.radius*rgen.rand_uniform(0.7, 0.9));
			if (t.has_gate && fabs(pos - t.gate_pos) < 2.0*radius) continue; // too close to the gate
			float const zval(t.p[0].z + height), hlen(sqrt(t.radius*t.radius - height*height) + 2.0*radius);
			cube_t pipe;
			set_wall_width(pipe, t.p[0][!t.dim], hlen, !t.dim);
			set_wall_width(pipe, pos,  radius, t.dim);
			set_wall_width(pipe, zval, radius, 2);
			objs.emplace_back(pipe, TYPE_PIPE, 0, !t.dim, 0, RO_FLAG_NOCOLL, 1.0, SHAPE_CYLIN, choose_pipe_color(rgen)); // horizontal, room_id=0
		} // for n
	} // for t
}

// *** Drawing ***

// tunnels really should be drawn as building interior verts rather than small static objects, but the code here is much more flexible and more efficient
void building_room_geom_t::add_tunnel(tunnel_seg_t const &t) {
	bool const shadowed(0); // not shadowed, since there's no light
	bool const dim(t.dim);
	unsigned const ndiv(48);
	float const length(t.get_length()), circumference(PI*t.radius), centerline(t.bcube.get_center_dim(!dim));
	float const side_tscale(2.0), len_tscale(side_tscale*length/circumference), end_tscale(len_tscale*2.0*t.radius/length);
	colorRGBA const wall_color(WHITE);
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_concrete_tid(), 16.0, shadowed), shadowed, 0, 1));
	unsigned const verts_start_ix(mat.itri_verts.size()), ixs_start_ix(mat.indices.size());
	// only the tunnel interior needs to be drawn, since it can't be viewed from the exterior; but drawing the exterior helps with debugging
	mat.add_ortho_cylin_to_verts(t.bcube, wall_color, dim, 0, 0, 1, 0, 1.0, 1.0, side_tscale, end_tscale, 0, ndiv, 0.0, 0, len_tscale, 0.0, t.room_conn);

	if (t.room_conn) { // rotate half cylinder into the proper orient
		float angle(0.0);
		if (!dim      ) {angle += PI_TWO;} // 90  degrees
		if (t.room_dir) {angle += PI    ;} // 180 degrees
		if (angle != 0.0) {rotate_verts(mat.itri_verts, (dim ? plus_y : plus_x), angle, t.p[0], verts_start_ix);}
		// draw wall sections connecting to the door
		cube_t conn_area(t.bcube_ext);
		conn_area.d[!dim][!t.room_dir] = centerline; // cut in half

		for (unsigned invert = 0; invert < 2; ++invert) { // draw ceiling and floor, both inner and outer
			mat.add_cube_to_verts(conn_area, wall_color, all_zeros, ~EF_Z12, 0, 0, 0, invert);
		}
		// draw another cube at the inside of the door to block the water
		mat.add_cube_to_verts(t.get_room_conn_block(), wall_color, all_zeros, (EF_Z1 | get_skip_mask_for_xy(dim))); // draw top, front, and back
		// draw the two C-shaped sides
		size_t const num_verts(mat.itri_verts.size() - verts_start_ix), num_ixs(mat.indices.size() - ixs_start_ix);
		float const door_pos(t.bcube_ext.d[!dim][t.room_dir]), center_dim(t.bcube.get_center_dim(dim)), tscale(1.0/PI);

		for (unsigned d = 0; d < 2; ++d) {
			size_t const vert_ix_off(mat.itri_verts.size() - verts_start_ix);

			// copy cylin side verts and transform them into the correct place
			for (unsigned i = 0; i < num_verts; ++i) {
				auto v(mat.itri_verts[verts_start_ix+i]);

				if ((v.v[dim] < center_dim) ^ bool(d)) { // keep this vertex on the cylinder, but move to the other hemisphere
					v.v[!dim] = 2.0*centerline - v.v[!dim]; // reflect about center_dim
				}
				else { // move this vertex to the wall next to the door
					v.v[!dim] = door_pos;
					v.v[ dim] = t.bcube.d[dim][d]; // move to the other end
				}
				v.set_ortho_norm(dim, !d); // recalculate normal
				v.t[0] *= 4.0*tscale; // hack to make the texture coordinates look less distorted
				v.t[1] *= tscale;
				mat.itri_verts.push_back(v);
			} // for i
			// indices have the same quad topology, so copy and offset them
			for (unsigned i = 0; i < num_ixs; ++i) {mat.indices.push_back(mat.indices[ixs_start_ix+i] + vert_ix_off);}
		} // for d
	}
	// draw closed ends in black so that they appear to extend into darkness
	if (t.closed_ends[0] || t.closed_ends[1]) {
		assert(!t.room_conn); // not supported
		unsigned const start_ix(mat.itri_verts.size());
		mat.add_ortho_cylin_to_verts(t.bcube, BLACK, t.dim, t.closed_ends[0], t.closed_ends[1], 1, 0, 1.0, 1.0, 1.0, 1.0, 1);
		for (auto v = mat.itri_verts.begin()+start_ix; v != mat.itri_verts.end(); ++v) {v->set_norm_to_zero();} // zero normal to prevent wet effect
	}
	// draw gate if present
	if (t.has_gate) {
		assert(!t.room_conn); // not supported
		unsigned const num_bars(8), bar_ndiv(16);
		float const bar_spacing(2*t.radius/(num_bars + 1)), bar_radius(0.025*t.radius);
		float const zc(t.bcube.zc()), rsq(t.radius*t.radius);
		float bar_pos(t.bcube.d[!dim][0] + bar_spacing);
		rgeom_mat_t &bar_mat(get_material(tid_nm_pair_t(get_texture_by_name("metals/67_rusty_dirty_metal.jpg")), shadowed, 0, 1)); // inc_shadows=0 (no light), small=1
		colorRGBA const bar_color(DK_GRAY);

		for (unsigned n = 0; n < num_bars; ++n, bar_pos += bar_spacing) {
			cube_t bar;
			float const dist(bar_pos - centerline), hheight(sqrt(rsq - dist*dist) + bar_radius), len_tc(0.25*hheight/bar_radius);
			set_wall_width(bar, zc,         hheight,     2  );
			set_wall_width(bar, t.gate_pos, bar_radius,  dim);
			set_wall_width(bar, bar_pos,    bar_radius, !dim);
			bar_mat.add_vcylin_to_verts(bar, bar_color, 0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, bar_ndiv, 0.0, 0, len_tc); // draw sides but not ends
		}
	}
}

void building_room_geom_t::add_tunnel_water(tunnel_seg_t const &t) {
	if (t.water_level <= 0.0) return;
	// draw water surface
	bool const dim(t.dim);
	float const tscale(1.0/t.radius);
	float tscale_xy[2] = {tscale, tscale};
	tscale_xy[!dim] *= 1.0/(1.0 + 10.0*fabs(t.water_flow)); // stretch along tunnel length more with higher water flow
	rgeom_mat_t &water_mat(get_material(tid_nm_pair_t(FOAM_TEX, -1, tscale_xy[0], tscale_xy[1]), 0, 1)); // unshadowed, dynamic
	cube_t water(t.bcube);
	water.z2() = t.bcube.z1() + t.water_level;
	float const dist(t.radius - t.water_level), water_hwidth(sqrt(t.radius*t.radius - dist*dist));
	set_wall_width(water, t.bcube.get_center_dim(!dim), water_hwidth, !dim);
	if (t.room_conn) {water.d[!dim][t.room_dir] = t.get_room_conn_block().d[!dim][!t.room_dir];} // ends flush with conn block
	(dim ? water_mat.tex.txoff : water_mat.tex.tyoff) = fract(0.15*(tfticks/TICKS_PER_SECOND)*t.water_flow); // animate the water texture
	water_mat.add_cube_to_verts(water, DK_BROWN, all_zeros, ~EF_Z2); // draw top surface only
}


