// 3D World - Building Underground/Basement Tunnels
// by Frank Gennari 10/10/2024

#include "function_registry.h"
#include "buildings.h"

extern double tfticks;

float query_min_height(cube_t const &c, float stop_at);
colorRGBA choose_pipe_color(rand_gen_t &rgen);
room_object_t get_open_false_door(room_object_t const &c);

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
	bcube_ext = bcube_draw = bcube;
}
void tunnel_seg_t::set_as_room_conn(bool rdir, float wall_gap) {
	room_conn = 1;
	room_dir  = rdir;
	bcube_ext.d[!dim][room_dir] += (room_dir ? 1.0 : -1.0)*wall_gap;
}
cube_t tunnel_seg_t::get_walk_area(point const &pos, float user_radius) const {
	assert(radius > user_radius); // otherwise we can't walk in this tunnel
	cube_t walk_area(bcube_ext);
	float const walk_width(0.1*radius), blocked_width(radius - walk_width); // area inside the tunnel near the center where the player can walk
	if (room_conn) {walk_area.d[!dim][!room_dir] += (room_dir ? 1.0 : -1.0)*blocked_width;} // shrink on non-room side
	else {walk_area.expand_in_dim(!dim, -blocked_width);} // shrink on both sides
	// prevent the player from walking off a closed end, through a gate, or through the connecting room door frame
	if (closed_ends[0]) {walk_area.d[dim][0] += user_radius;}
	if (closed_ends[1]) {walk_area.d[dim][1] -= user_radius;}

	if (room_conn) {
		float const shrink_amt(min(user_radius, (fabs(pos[!dim] - bcube.get_center_dim(!dim)) - radius + user_radius)));
		if (shrink_amt > 0.0) {walk_area.expand_in_dim(dim, -shrink_amt);} // 45 degree XY edges at the entrance door
	}
	if (has_gate) {
		if (pos[dim] < gate_pos) {min_eq(walk_area.d[dim][1], gate_pos-user_radius);} // player at low end, clamp high end
		else                     {max_eq(walk_area.d[dim][0], gate_pos+user_radius);} // player at high end, clamp low end
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
point tunnel_seg_t::get_room_conn_pt(float zval) const {
	assert(room_conn);
	point pt(0.0, 0.0, zval);
	pt[ dim] = bcube_ext.get_center_dim(dim);
	pt[!dim] = bcube_ext.d[!dim][room_dir];
	return pt;
}
bool tunnel_seg_t::is_blocked_by_gate(point const &p1, point const &p2) const {
	return (has_gate && gate_pos > min(p1[dim], p2[dim]) && gate_pos < max(p1[dim], p2[dim]));
}
tunnel_seg_t tunnel_seg_t::connect_segment_to(point const &pos, bool conn_dir, unsigned new_tseg_ix, bool is_end) {
	tunnel_seg_t tseg(p[conn_dir], pos, radius); // same radius
	bool const tseg_conn_dir(tseg.p[0] == pos);
	assert(conn_ix[conn_dir] < 0); // can't connect two tunnels to the same point
	conn_ix[conn_dir] = new_tseg_ix; // parent => child
	tseg.conn_ix[tseg_conn_dir] = tseg_ix; // child => parent
	tseg.conn_room_ix = conn_room_ix;
	tseg.tseg_ix      = new_tseg_ix;
	tseg.water_level  = water_level;
	tseg.water_flow   = water_flow;
	if ((tseg.dim != dim) && (conn_dir ^ dim ^ tseg_conn_dir)) {tseg.water_flow *= -1.0;}
	if (is_end) {tseg.closed_ends[!tseg_conn_dir] = 1;} // closed at unconnected end
	closed_ends[conn_dir] = 0; // parent is no longer closed
	return tseg;
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
	vect_room_object_t &objs(interior->room_geom->objs);

	for (unsigned d = 0; d < 2; ++d) { // check both door locations
		cube_with_ix_t const &door(doors[d]);
		if (door.is_all_zeros()) continue;
		bool const dim(door.ix >> 1), dir(door.ix & 1), is_open(room.has_tunnel_conn());
		unsigned flags(RO_FLAG_NOCOLL);

		if (is_open) { // add an open door; should this have collisions enabled?
			flags |= (RO_FLAG_HAS_EXTRA | RO_FLAG_OPEN); // make this an open vault/blast door
			//cube_t c(door); c.z1() = c.z2(); c.z2() += 1.0; objs.emplace_back(c, TYPE_DBG_SHAPE, room_id, dim, dir, 0, 1.0, SHAPE_CUBE, RED); // TESTING
		}
		else { // add a closed door
			if      (is_house)           {flags |= RO_FLAG_IS_HOUSE ;}
			else if ((room_id & 3) == 0) {flags |= RO_FLAG_HAS_EXTRA;} // make this a vault/blast door 25% of the time
		}
		objs.emplace_back(door, TYPE_FALSE_DOOR, room_id, dim, dir, flags, 1.0, SHAPE_CUBE, WHITE);
		if (is_open) {objs.emplace_back(get_open_false_door(objs.back()), TYPE_COLLIDER, room_id, dim, dir, RO_FLAG_INVIS);} // open door is a collider for player
		// add door blocker to avoid placing objects in front of the door; needed even for hallways to avoid placing an overlapping outlet
		cube_t blocker(door);
		blocker.d[dim][!dir] = room.d[dim][dir] + (dir ? -1.0 : 1.0)*get_doorway_width(); // add clearance
		objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, dir, RO_FLAG_INVIS);
	} // for d
}

bool building_t::is_tunnel_bcube_placement_valid(cube_t const &tunnel_bc) const {
	if (cube_intersects_basement_or_extb_room(tunnel_bc))             return 0;
	if (query_min_height(tunnel_bc, tunnel_bc.z2()) < tunnel_bc.z2()) return 0; // check for terrain clipping through ceiling
	if (!is_basement_room_not_int_bldg(tunnel_bc))                    return 0;
	return 1;
}
bool building_t::is_tunnel_placement_valid(point const &p1, point const &p2, float radius) const {
	return is_tunnel_bcube_placement_valid(tunnel_seg_t(p1, p2, radius).bcube);
}
void building_t::try_extend_tunnel(point &p1, point &p2, float max_extend, float check_radius, bool dim, bool dir) const {
	unsigned const num_steps = 10;
	float const step_len(max_extend/num_steps);

	for (unsigned n = 0; n < num_steps; ++n) {
		if (dir) {
			point p2e(p2);
			p2e[dim] += step_len;
			if (!is_tunnel_placement_valid(p2, p2e, check_radius)) break; // can't extend further in this dir
			p2 = p2e; // accept the new length
		}
		else {
			point p1e(p1);
			p1e[dim] -= step_len;
			if (!is_tunnel_placement_valid(p1e, p1, check_radius)) break; // can't extend further in this dir
			p1 = p1e; // accept the new length
		}
	}
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
		float const max_extend(max_len - min_len);
		// try extending tunnel in both directions
		for (unsigned e = 0; e < 2; ++e) {try_extend_tunnel(p1, p2, max_extend, check_radius, !dim, e);}
		// split into three segments (center, left, right)
		vect_tunnel_seg_t &tunnels(interior->tunnels);
		unsigned const tseg_ix(tunnels.size());
		float const gate_dist_from_end(5.0*floor_spacing);
		float const water_level(rgen.rand_uniform(0.0, 1.0)*0.2*radius), water_flow(rgen.signed_rand_float());
		point pa(p1), pb(p2);
		pa[!dim] = door.d[!dim][0] + sm_shift_val; // left  end of room connection
		pb[!dim] = door.d[!dim][1] - sm_shift_val; // right end of room connection
		tunnel_seg_t tseg_c(pa, pb, radius);
		tseg_c.set_as_room_conn(!dir, wall_gap);
		tseg_c.water_level  = water_level;
		tseg_c.water_flow   = water_flow;
		tseg_c.conn_room_ix = room_id;
		tseg_c.tseg_ix      = tseg_ix;
		interior->add_tunnel_seg(tseg_c);

		for (unsigned e = 0; e < 2; ++e) { // {left, right}
			unsigned const cur_seg_ix(tunnels.size());
			point const &end_pt(e ? p2 : p1);
			// check if we can extend this tunnel before placing it, so that the query doesn't pick up a self intersection
			cube_t end_conn_bcube; end_conn_bcube.set_from_sphere(end_pt, check_radius);
			bool const can_extend(is_tunnel_bcube_placement_valid(end_conn_bcube));
			interior->connect_and_add_tunnel_seg(tunnels[tseg_ix], end_pt, e, !e, gate_dist_from_end);
			if (!can_extend) continue;
			// attempt to add a 90 degree bend
			bool const first_dir(rgen.rand_bool());

			for (unsigned n = 0; n < 2; ++n) {
				bool const bend_dir(bool(n) ^ first_dir);
				point p1b(end_pt), p2b(end_pt);
				point &new_end_pt(bend_dir ? p2b : p1b);
				new_end_pt[dim] += (bend_dir ? 1.0 : -1.0)*min_len;
				point p1bt(p1b), p2bt(p2b); // pull back the connecting end so that it doesn't intersect the parent tunnel segment
				(bend_dir ? p1bt : p2bt)[dim] += (bend_dir ? 1.0 : -1.0)*(check_radius + wall_thickness);
				if (!is_tunnel_placement_valid(p1bt, p2bt, check_radius)) continue; // can't place a tunnel of min length
				try_extend_tunnel(p1b, p2b, max_extend, check_radius, dim, bend_dir);
				interior->connect_and_add_tunnel_seg(tunnels[cur_seg_ix], new_end_pt, e, !bend_dir, gate_dist_from_end);
				tunnel_seg_t &parent(tunnels[cur_seg_ix]), &child(tunnels.back());
				// pull back cylinder draw bcubes to remove overlaps
				float const child_ext((bend_dir ? -1.0 : 1.0)*radius), parent_ext((e ? 1.0 : -1.0)*radius);
				child .bcube_draw.d[child .dim][!bend_dir] -= child_ext;
				parent.bcube_draw.d[parent.dim][ e       ] -= parent_ext;
				// extend bcubes to cover the gap at the corner so that the player stays in the tunnel
				child .bcube_ext .d[child .dim][!bend_dir] += child_ext;
				parent.bcube_ext .d[parent.dim][ e       ] += parent_ext;
				parent.add_bend_dir[e] = int(bend_dir); // bend is drawn by the parent
				parent.has_gate = 0; // remove gate from parent
				break; // don't need to check the other dir
			} // for n
		} // for e
		room.set_has_tunnel_conn();
		return 1; // done (only add tunnel conn to one end)
	} // for d
	return 0;
}

void building_interior_t::connect_and_add_tunnel_seg(tunnel_seg_t &parent, point const &conn_pt, bool parent_conn_dir, bool child_conn_dir, float gate_dist_from_end) {
	bool is_end(1); // starts as an end, but may be extended
	tunnel_seg_t tseg(parent.connect_segment_to(conn_pt, parent_conn_dir, tunnels.size(), is_end));
	if (is_end) {tseg.set_gate(conn_pt[tseg.dim] + gate_dist_from_end*(child_conn_dir ? 1.0 : -1.0));}
	add_tunnel_seg(tseg);
}
void building_interior_t::add_tunnel_seg(tunnel_seg_t const &tseg) {
	tunnels.emplace_back(tseg);
	basement_ext_bcube.assign_or_union_with_cube(tseg.bcube); // add to bcube
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

// *** queries and pathing ***

bool building_interior_t::point_in_tunnel(point const &pos, float expand) const {
	for (tunnel_seg_t const &t : tunnels) {
		if (t.bcube_ext.contains_pt_exp(pos, expand)) return 1;
	}
	return 0;
}
bool building_interior_t::point_near_tunnel_entrance(point const &pos) const {
	for (tunnel_seg_t const &t : tunnels) {
		if (!t.room_conn) continue;
		cube_t c(t.bcube_ext);
		c.union_with_cube(get_room(t.conn_room_ix)); // include the connected room
		if (c.contains_pt(pos)) return 1;
	}
	return 0;
}
int building_interior_t::get_tunnel_ix_for_point(point const &pos) const {
	for (auto t = tunnels.begin(); t != tunnels.end(); ++t) { // bcube containment first
		if (t->bcube.contains_pt(pos)) return (t - tunnels.begin());
	}
	for (auto t = tunnels.begin(); t != tunnels.end(); ++t) { // now include extended bcubes
		if (t->bcube_ext.contains_pt(pos)) return (t - tunnels.begin());
	}
	return -1; // not found
}
int building_interior_t::get_tunnel_ix_for_room(unsigned room_ix) const {
	for (auto t = tunnels.begin(); t != tunnels.end(); ++t) {
		if (t->room_conn && t->conn_room_ix == room_ix) return (t - tunnels.begin());
	}
	return -1; // not found
}

// Note: paths include end_pt and start_pt
void add_clamped_pt(tunnel_seg_t const &tseg, point pos, float radius, ai_path_t &path) {
	tseg.get_walk_area(pos, radius).clamp_pt_xy(pos);
	path.add(pos, 1); // fixed=1
}
bool building_interior_t::get_tunnel_path_from_room(point const &end_pt, unsigned room_ix, float radius, ai_path_t &path) const { // room => tunnel
	path.clear();
	room_t const &room(get_room(room_ix));
	if (!room.has_tunnel_conn()) return 0;
	int const pt_tix(get_tunnel_ix_for_point(end_pt));
	if (pt_tix < 0) return 0; // end_pt not in a tunnel
	if (tunnels[pt_tix].conn_room_ix != room_ix) return 0; // end_pt tunnel not connected to this room
	int const room_tix(get_tunnel_ix_for_room(room_ix));
	assert(room_tix >= 0); // should be found
	tunnel_seg_t const &tseg_r(tunnels[room_tix]);
	path.add(tseg_r.get_room_conn_pt(end_pt.z)); // room to tunnel transition
	
	if (room_tix != pt_tix) { // not in the room connector segment
		path.add(tseg_r.bcube.xc(), tseg_r.bcube.yc(), end_pt.z, 1); // center of tunnel entrance segment; fixed=1
		if (!get_tunnel_path(room_tix, pt_tix, -1, path)) {path.clear(); return 0;}
		reverse(path.begin()+2, path.end());
	}
	if (tunnels[pt_tix].is_blocked_by_gate(end_pt, path.back())) {path.clear(); return 0;} // blocked by gate
	add_clamped_pt(tunnels[pt_tix], end_pt, radius, path);
	return 1;
}
bool building_interior_t::get_tunnel_path_to_room(point const &start_pt, unsigned &room_ix, ai_path_t &path) const { // tunnel => room
	path.clear();
	int const pt_tix(get_tunnel_ix_for_point(start_pt));
	if (pt_tix < 0) return 0; // start_pt not in a tunnel
	room_ix = tunnels[pt_tix].conn_room_ix;
	int const room_tix(get_tunnel_ix_for_room(room_ix));
	assert(room_tix >= 0); // should be found
	tunnel_seg_t const &tseg(tunnels[room_tix]);
	path.add(start_pt);

	if (room_tix != pt_tix) { // not in the room connector segment
		if (!get_tunnel_path(pt_tix, room_tix, -1, path)) {path.clear(); return 0;}
		reverse(path.begin()+1, path.end());
		path.add(tseg.bcube.xc(), tseg.bcube.yc(), start_pt.z, 1); // center of tunnel entrance segment; fixed=1
	}
	path.add(tseg.get_room_conn_pt(start_pt.z), 1); // tunnel to room transition; fixed=1
	return 1;
}
bool building_interior_t::get_tunnel_path_two_pts(point const &start_pt, point const &end_pt, float radius, ai_path_t &path) const {
	path.clear();
	int const s_tix(get_tunnel_ix_for_point(start_pt));
	if (s_tix < 0) return 0;
	int const e_tix(get_tunnel_ix_for_point(end_pt));
	if (s_tix < 0) return 0;
	if (tunnels[s_tix].conn_room_ix != tunnels[e_tix].conn_room_ix) return 0; // different tunnel networks
	path.add(start_pt);

	if (s_tix != e_tix) { // not in the same segment
		if (!get_tunnel_path(s_tix, e_tix, -1, path)) {path.clear(); return 0;}
		reverse(path.begin()+1, path.end());
	}
	if (tunnels[e_tix].is_blocked_by_gate(end_pt, path.back())) {path.clear(); return 0;} // blocked by gate
	add_clamped_pt(tunnels[e_tix], end_pt, radius, path);
	return 1;
}
bool building_interior_t::get_tunnel_path(unsigned tix1, unsigned tix2, int prev_tix, ai_path_t &path) const {
	if (tix1 == tix2) return 1; // paths connected
	tunnel_seg_t const &tseg(tunnels[tix1]);

	for (unsigned d = 0; d < 2; ++d) { // try both ends
		int const new_tix(tseg.conn_ix[d]);
		if (new_tix < 0 || new_tix == prev_tix)          continue;
		if (!get_tunnel_path(new_tix, tix2, tix1, path)) continue;
		assert(!path.empty());
		point const &tseg_pt(tseg.p[d]);
		if (tseg.is_blocked_by_gate(tseg_pt, path.back())) continue; // blocked by gate
		path.add(tseg_pt.x, tseg_pt.y, path.back().z, 1); // fixed=1
		return 1;
	} // for d
	return 0;
}

// *** Drawing ***

// tunnels really should be drawn as building interior verts rather than small static objects, but the code here is much more flexible and more efficient
void building_room_geom_t::add_tunnel(tunnel_seg_t const &t) {
	bool const shadowed(0); // not shadowed, since there's no light
	bool const dim(t.dim);
	unsigned const ndiv(48);
	float const length(t.get_length()), circumference(PI*t.radius), centerline(t.bcube.get_center_dim(!dim)), bar_radius(0.025*t.radius);
	float const side_tscale(2.0), len_tscale(side_tscale*length/circumference), end_tscale(len_tscale*2.0*t.radius/length);
	colorRGBA const wall_color(WHITE);
	rgeom_mat_t &mat(get_material(tid_nm_pair_t(get_concrete_tid(), 16.0, shadowed), shadowed, 0, 1));

	// add smaller side pipes for non-short, non-connector segments; must use the same material and be drawn before the main tunnel walls
	if (!t.room_conn && length > 4.0*t.radius) {
		rand_gen_t rgen;
		rgen.set_state(t.conn_room_ix+1, t.tseg_ix+1);
		rgen.rand_mix();
		unsigned const max_pipes(3), mod_val(max_pipes + 1), num_pipes(rgen.rand() % mod_val); // 0-3
		float avoid_vals[mod_val]={}, avoid_radii[mod_val]={};
		unsigned num_avoid(0);
		if (t.has_gate) {avoid_vals[num_avoid] = t.gate_pos; avoid_radii[num_avoid++] = bar_radius;}

		for (unsigned n = 0; n < num_pipes; ++n) {
			float const radius(rgen.rand_uniform(0.1, 0.3)*t.radius), end_pad(2.0*radius);
			float const pos(rgen.rand_uniform(t.bcube_draw.d[dim][0]+end_pad, t.bcube_draw.d[dim][1]-end_pad));
			bool blocked(0);

			for (unsigned N = 0; N < num_avoid; ++N) {
				if (fabs(pos - avoid_vals[N]) < (end_pad + avoid_radii[N])) {blocked = 1; break;} // too close to gate or another pipe
			}
			if (blocked) continue;
			bool const dir(rgen.rand_bool());
			point center(0.0, 0.0, t.bcube.zc());
			center[ dim] = pos; // pos along tunnel
			center[!dim] = t.bcube.d[!dim][dir]; // wall of tunnel
			cube_t pipe(center);
			pipe.expand_in_dim(   2,     radius);
			pipe.expand_in_dim( dim,     radius);
			pipe.expand_in_dim(!dim, 0.5*radius); // length
			mat.add_ortho_cylin_to_verts(pipe, wall_color, !dim, 0, 0, 1); // skip ends; two_sided=1
			// draw the open/interior end transparent to create a hole in the outer pipe
			mat.add_ortho_cylin_to_verts(pipe, ALPHA0, !dim,  dir, !dir, 0, 0, 1.0, 1.0, 1.0, 1.0, 1); // transparent cut; skip_sides=1
			mat.add_ortho_cylin_to_verts(pipe, BLACK,  !dim, !dir,  dir, 0, 1, 1.0, 1.0, 1.0, 1.0, 1); // black interior end cap; skip_sides=1
			avoid_vals[num_avoid] = pos; avoid_radii[num_avoid++] = radius; // should we check dir or track it separately? may not matter much
		} // for n
	}
	unsigned const verts_start_ix(mat.itri_verts.size()), ixs_start_ix(mat.indices.size());
	// only the tunnel interior needs to be drawn, since it can't be viewed from the exterior; but drawing the exterior helps with debugging
	mat.add_ortho_cylin_to_verts(t.bcube_draw, wall_color, dim, 0, 0, 1, 0, 1.0, 1.0, side_tscale, end_tscale, 0, ndiv, 0.0, 0, len_tscale, 0.0, t.room_conn);

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
	// draw walls for curved bends
	for (unsigned d = 0; d < 2; ++d) {
		if (t.add_bend_dir[d] < 0) continue; // no bend on this end
		bool const conn_dir(t.add_bend_dir[d]);
		int const off_ixs[8] = {1, 2, 0, 3, 3, 2, 0, 1};
		unsigned const off_ix(off_ixs[4*dim + 2*d + conn_dir]);
		float const s_offset(0.25*off_ix);
		point pos(t.p[d]);
		pos.x += ((off_ix == 1 || off_ix == 2) ? 1.0 : -1.0)*t.radius;
		pos.y += ((off_ix == 2 || off_ix == 3) ? 1.0 : -1.0)*t.radius;
		unsigned itris_start(mat.itri_verts.size()), ixs_start(mat.indices.size());
		mat.add_vert_torus_to_verts(pos, t.radius, t.radius, wall_color, side_tscale, 0, 2, s_offset, ndiv); // low_detail=0, half_or_quarter=2 (quarter)
		add_inverted_triangles(mat.itri_verts, mat.indices, itris_start, ixs_start); // draw the back side
	} // for d
	// draw closed ends in black so that they appear to extend into darkness
	if (t.closed_ends[0] || t.closed_ends[1]) {
		assert(!t.room_conn); // not supported
		unsigned const start_ix(mat.itri_verts.size());
		mat.add_ortho_cylin_to_verts(t.bcube_draw, BLACK, t.dim, t.closed_ends[0], t.closed_ends[1], 1, 0, 1.0, 1.0, 1.0, 1.0, 1);
		for (auto v = mat.itri_verts.begin()+start_ix; v != mat.itri_verts.end(); ++v) {v->set_norm_to_zero();} // zero normal to prevent wet effect
	}
	// draw gate if present
	if (t.has_gate) {
		assert(!t.room_conn); // not supported
		unsigned const num_bars(8), bar_ndiv(16);
		float const bar_spacing(2*t.radius/(num_bars + 1)), zc(t.bcube.zc()), rsq(t.radius*t.radius);
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
	cube_t water(t.bcube_draw);
	water.z2() = t.bcube.z1() + t.water_level;
	float const dist(t.radius - t.water_level), water_hwidth(sqrt(t.radius*t.radius - dist*dist));
	set_wall_width(water, t.bcube.get_center_dim(!dim), water_hwidth, !dim);
	if (t.room_conn) {water.d[!dim][t.room_dir] = t.get_room_conn_block().d[!dim][!t.room_dir];} // ends flush with conn block
	float const flow_val(0.15*(tfticks/TICKS_PER_SECOND)*t.water_flow);
	float tex_add[2] = {0.0, 0.0};
	tex_add[!dim] = fract(flow_val); // animate the water texture
	water_mat.add_cube_to_verts(water, DK_BROWN, all_zeros, ~EF_Z2, 0, 0, 0, 0, 0, tex_add[0], tex_add[1]); // draw top surface only

	// draw water at bends
	for (unsigned d = 0; d < 2; ++d) {
		if (t.add_bend_dir[d] < 0) continue; // no bend on this end
		bool const conn_dir(t.add_bend_dir[d]);
		cube_t bend_water(t.bcube_ext);
		bend_water.d[dim][!d] = t.bcube_draw.d[dim][d]; // constrain to end only
		bend_water.z2() = water.z2();

		if (water_mat.tex.tscale_x != water_mat.tex.tscale_y) { // average out the texture scale
			water_mat.tex.tscale_x = water_mat.tex.tscale_y = sqrt(water_mat.tex.tscale_x*water_mat.tex.tscale_y);
		}
		// water moves diagonally; can we rotate the quad to create seamless movement?
		tex_add[!dim] = fract(SQRTOFTWOINV*flow_val);
		tex_add[ dim] = tex_add[!dim] * ((d ^ t.dim ^ conn_dir) ? 1.0 : -1.0);
		water_mat.add_cube_to_verts(bend_water, DK_BROWN, all_zeros, ~EF_Z2, 0, 0, 0, 0, 0, tex_add[0], tex_add[1]); // draw top surface only
	} // for d
}


