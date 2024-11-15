// 3D World - Building Basement Backrooms
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"


bool try_place_wall(cube_t const &place_area, cube_t const &wall_area, cube_t const &cent_area, bool dim, float len_min, float len_max, float half_thick, float min_gap,
	vect_cube_t &walls, vect_cube_t &blockers, rand_gen_t &rgen)
{
	assert(len_min < len_max);
	float const wall_len   (rgen.rand_uniform(len_min, len_max));
	float const wall_pos   (rgen.rand_uniform(wall_area.d[ dim][0], wall_area.d[ dim][1])); // position of wall centerline
	float const wall_center(rgen.rand_uniform(cent_area.d[!dim][0], cent_area.d[!dim][1])); // center of wall
	cube_t wall(wall_area); // copy zvals
	set_wall_width(wall, wall_pos,    half_thick,    dim);
	set_wall_width(wall, wall_center, 0.5*wall_len, !dim);
	// wall can extend outside the room, and we generally want it to end at the room edge in some cases
	if (wall.d[!dim][0] < place_area.d[!dim][0]) {wall.d[!dim][0] = place_area.d[!dim][0];} // clip to room bounds
	else if (wall.d[!dim][0] - place_area.d[!dim][0] < min_gap) {wall.d[!dim][0] += min_gap;} // force min_gap from wall
	if (wall.d[!dim][1] > place_area.d[!dim][1]) {wall.d[!dim][1] = place_area.d[!dim][1];} // clip to room bounds
	else if (place_area.d[!dim][1] - wall.d[!dim][1] < min_gap) {wall.d[!dim][1] -= min_gap;} // force min_gap from wall
	if (wall.get_sz_dim(!dim) < len_min) return 0; // too short
	if (has_bcube_int(wall, blockers))   return 0; // too close to a previous wall in this dim
	walls.push_back(wall);
	wall.expand_by_xy(min_gap); // require a gap around this wall - but can still have intersections in the other dim
	blockers.push_back(wall);
	return 1;
}

void partition_cubes_into_conn_groups(vect_cube_t const &cubes, vector<vect_cube_t> &groups, float pad=0.0) {
	groups.clear();
	vect_cube_t ungrouped(cubes); // all cubes start ungrouped

	while (!ungrouped.empty()) {
		vect_cube_t group;
		group.push_back(ungrouped.back());
		ungrouped.pop_back();

		for (unsigned ix = 0; ix < group.size(); ++ix) { // process each cube added to the group
			cube_t cur(group[ix]);
			cur.expand_by_xy(pad);

			for (unsigned i = 0; i < ungrouped.size(); ++i) {
				cube_t &cand(ungrouped[i]);
				if (!cand.intersects(cur)) continue; // includes adjacency
				group.push_back(cand);
				cand = ungrouped.back(); // remove cand by replacing it with the last ungrouped cube
				ungrouped.pop_back();
				--i;
			} // for i
		} // for ix
		groups.push_back(group);
	} // while
}

struct cmp_cube_x1_y1 {
	bool operator()(cube_t const &a, cube_t const &b) const {return ((a.x1() == b.x1()) ? (a.y1() < b.y1()) : (a.x1() < b.x1()));}
};
void invert_walls(cube_t const &room, vect_cube_t const walls[2], vect_cube_t &space, float pad=0.0) {
	space.clear();
	space.push_back(room);
	
	for (unsigned d = 0; d < 2; ++d) {
		for (cube_t const &wall : walls[d]) {
			cube_t sub(wall);
			sub.expand_by_xy(pad);
			subtract_cube_from_cubes(sub, space);
		}
	}
	if (space.empty()) return;
	// max merge adjacent cubes in Y
	sort(space.begin(), space.end(), cmp_cube_x1_y1());
	auto i(space.begin()+1), o(i); // skip first

	for (; i != space.end(); ++i) {
		cube_t &c(*(o-1)); // last output cube
		if (c.x1() == i->x1() && c.x2() == i->x2() && c.y2() == i->y1()) {c.y2() = i->y2();} // extend in +Y
		else {*o++ = *i;} // keep this cube
	} // for i
	space.erase(o, space.end());
}
void resize_cubes_xy(vect_cube_t &cubes, float val) { // val can be positive or negative
	for (cube_t &c : cubes) {c.expand_by_xy(val);}
}
void transpose_cube_xy(cube_t &c) {
	swap(c.x1(), c.y1());
	swap(c.x2(), c.y2());
}
void transpose_cubes_xy(vect_cube_t &cubes) {
	for (cube_t &c : cubes) {transpose_cube_xy(c);}
}

struct group_range_t {
	unsigned gix;
	float lo, hi;
	group_range_t(unsigned gix_, float lo_, float hi_) : gix(gix_), lo(lo_), hi(hi_) {}
	void update(float lo_, float hi_) {min_eq(lo, lo_); max_eq(hi, hi_);}
	float len() const {return (hi - lo);}
};
struct cube_by_center_dim_descending {
	unsigned d;
	cube_by_center_dim_descending(unsigned dim) : d(dim) {}
	bool operator()(cube_t const &a, cube_t const &b) const {return (b.get_center_dim(d) < a.get_center_dim(d));}
};
struct cube_by_xy_dim {
	unsigned d;
	cube_by_xy_dim(unsigned dim) : d(dim) {}
	bool operator()(cube_t const &a, cube_t const &b) const {return ((a.d[d][0] == b.d[d][0]) ? (a.d[!d][0] < b.d[!d][0]) : (a.d[d][0] < b.d[d][0]));}
};

void building_t::add_backrooms_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t &rooms_to_light) {
	//highres_timer_t timer("Add Backrooms Objs"); // up to ~2ms
	assert(has_room_geom());
	float const floor_spacing(get_window_vspace()), wall_thickness(1.2*get_wall_thickness()), wall_half_thick(0.5*wall_thickness); // slightly thicker than regular walls
	float const ceiling_z(zval + get_floor_ceil_gap()); // Note: zval is at floor level, not at the bottom of the room
	//float const tot_light_amt(room.light_intensity /*+ floor_spacing*room.get_light_amt()*/); // ???
	float const tot_light_amt(0.75);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());
	if (interior->room_geom->backrooms_start == 0) {interior->room_geom->backrooms_start = objs_start;}
	rgen.rseed1 += 123*floor_ix; // make it unique per floor

	// find the shared wall with the basement/parking garage and calculate true room bounds
	bool const sw_dim(interior->extb_wall_dim), sw_dir(!interior->extb_wall_dir);
	float const shared_extend((sw_dir ? -1.0 : 1.0)*0.5*get_wall_thickness()); // account for extra shift applied to shared wall to keep it from clipping into the basement/PG
	room_t true_room(room);
	true_room.d[sw_dim][sw_dir] += shared_extend;
	vector3d const sz(true_room.get_size());
	vect_cube_t blockers_per_dim[2], walls_per_dim[2], pillar_avoid;
	float const doorway_width(get_doorway_width()), doorway_hwidth(0.5*doorway_width), min_gap(1.2*doorway_width), min_edge_gap(min_gap + wall_thickness);

	// find entrance door and add it to wall blockers
	door_t const &ent_door(interior->get_ext_basement_door());
	pillar_avoid.push_back(ent_door.get_clearance_bcube());
	blockers_per_dim[!ent_door.dim].push_back(pillar_avoid.back());

	// add ext basement connector doors as well
	for (door_stack_t const &ds : interior->door_stacks) {
		if (!ds.get_bldg_conn()) continue;
		pillar_avoid.push_back(ds.get_clearance_bcube());
		blockers_per_dim[!ds.dim].push_back(pillar_avoid.back());
	}

	// find any stairs in this room and add to both wall blockers
	for (stairwell_t const &s : interior->stairwells) {
		if (s.z1() >= ceiling_z || s.z2() <= zval) continue; // wrong floor
		cube_t stairs_bcube(get_stairs_bcube_expanded(s, doorway_width, wall_thickness, doorway_width));
		if (!room.intersects(stairs_bcube)) continue;
		bool const dir(rgen.rand_bool());
		stairs_bcube.d[s.dim][dir] += (dir ? 1.0 : -1.0)*(doorway_width - wall_thickness); // add an extra doorway width gap on a random side to avoid blocking a hallway
		pillar_avoid.push_back(stairs_bcube);
		for (unsigned d = 0; d < 2; ++d) {blockers_per_dim[d].push_back(stairs_bcube);}
	}

	// add random interior walls to create an initial maze
	cube_t const place_area(get_walkable_room_bounds(true_room));
	cube_t wall_area(place_area);
	wall_area.expand_by_xy(-min_edge_gap);
	if (min(wall_area.dx(), wall_area.dy()) < 2.0*floor_spacing) return; // room too small to place walls - shouldn't happen
	set_cube_zvals(wall_area, zval, ceiling_z);
	float const min_side(min(sz.x, sz.y)), area(sz.x*sz.y);
	float const wall_len_min(1.0*floor_spacing), wall_len_max(max(0.25f*min_side, 1.5f*wall_len_min)), wall_len_avg(0.5*(wall_len_min + wall_len_max)); // min_gap/wall_len_min = 0.6
	unsigned const num_walls(round_fp(rgen.rand_uniform(1.6, 2.0)*area/(wall_len_avg*wall_len_avg))); // add a bit of random density variation
	colorRGBA const wall_color(WHITE); // to match exterior walls, ceilings, and floors

	for (unsigned n = 0; n < num_walls; ++n) {
		for (unsigned m = 0; m < 10; ++m) { // 10 tries to place a wall
			bool const dim(rgen.rand_bool()); // should this be weighted by the room aspect ratio?
			if (try_place_wall(place_area, wall_area, place_area, dim, wall_len_min, wall_len_max, wall_half_thick, min_gap, walls_per_dim[dim], blockers_per_dim[dim], rgen)) break;
		}
	}

	// find long lines of sight and add extra walls to block them
	float const max_space_factor = 0.5; // relative to room size
	vect_cube_t space, extra_walls;

	for (unsigned dim = 0; dim < 2; ++dim) {
		// Note: wall ends may be moved slightly to snap to orthogonal walls, so this test isn't perfectly accurate
		if (dim) { // transpose cubes so that max merge is in Y rather than X
			vect_cube_t walls_transpose[2] = {walls_per_dim[1], walls_per_dim[0]};
			for (unsigned d = 0; d < 2; ++d) {transpose_cubes_xy(walls_transpose[d]);}
			cube_t place_area_transpose(place_area);
			transpose_cube_xy(place_area_transpose);
			invert_walls(place_area_transpose, walls_transpose, space); // no padding
			transpose_cubes_xy(space); // transpose back
		}
		else {invert_walls(place_area, walls_per_dim, space);} // no padding
		float const max_space(max_space_factor*sz[dim]);
		
		for (cube_t &s : space) {
			float const space_len(s.get_sz_dim(dim));
			if (space_len < max_space)         continue; // small enough
			if (has_bcube_int(s, extra_walls)) continue; // covered (at least partially) by a previously placed extra wall
			cube_t cent_area(s);
			cent_area.expand_in_dim(dim, -0.375*space_len); // restrict to central 25%
			set_cube_zvals(cent_area, wall_area.z1(), wall_area.z2()); // set zvals the same as walls
			float const space_width(s.get_sz_dim(!dim)), min_len(max(wall_len_min, space_width)), max_len(max(wall_len_max, 2.0f*space_width)); // must cross the space

			for (unsigned m = 0; m < 10; ++m) { // 10 tries to place a wall
				if (try_place_wall(place_area, cent_area, cent_area, dim, min_len, max_len, wall_half_thick, min_gap, walls_per_dim[dim], blockers_per_dim[dim], rgen)) {
					extra_walls.push_back(walls_per_dim[dim].back()); // record in case this blocks other long spaces in this loop
					break; // done
				}
			}
		} // for s
		extra_walls.clear();
	} // for dim
	// what if we instead calculate the shortest path from the entrance door to the far wall and add walls until it's some min length? is the runtime acceptable?

	// shift wall ends to remove small gaps and stubs
	float const wall_end_ext(doorway_width); // needed to handle two orthogonal walls with nearby corners but no projection

	for (unsigned dim = 0; dim < 2; ++dim) {
		vect_cube_t &walls(walls_per_dim[dim]);

		for (cube_t &wall : walls) {
			for (unsigned d = 0; d < 2; ++d) { // check each end
				float &val(wall.d[!dim][d]);
				if (val == place_area.d[!dim][d]) continue; // at exterior wall, skip

				for (cube_t &w : walls_per_dim[!dim]) {
					if (wall.d[dim][0] > w.d[dim][1] || wall.d[dim][1] < w.d[dim][0]) { // no projection
						if (wall.d[dim][0] > w.d[dim][1]+wall_end_ext || wall.d[dim][1] < w.d[dim][0]-wall_end_ext) continue; // no extended projection
						// handle nearby corners by moving wall away from the corner if needed; doesn't break because the edge may be moved again in the else case
						if      (d == 0 && val > w.d[!dim][1]) {max_eq(val, w.d[!dim][1]+min_gap);} // ensure left  space between walls is at least min_gap
						else if (d == 1 && val < w.d[!dim][0]) {min_eq(val, w.d[!dim][0]-min_gap);} // ensure right space between walls is at least min_gap
					}
					else {
						// if we already aligned <w> to <wall> then we now have a corner; extend val to the far edge of <w> to fill in the corner area
						bool const is_corner(wall.d[dim][0] == w.d[dim][1] || wall.d[dim][1] == w.d[dim][0]);
						float const edge_pos(w.d[!dim][bool(!d) ^ is_corner]);
						if (fabs(val - edge_pos) < min_gap) {val = edge_pos; break;}
					}
				} // for w
			} // for d
			if (wall.get_sz_dim(!dim) < wall_len_min) {wall = cube_t();} // too short, drop
		} // for wall
		walls.erase(std::remove_if(walls.begin(), walls.end(), [](cube_t const &c) {return c.is_all_zeros();}), walls.end());
	} // for dim

	// find areas of empty space
	float const pad(wall_half_thick), nav_pad(doorway_hwidth); // min radius for player and AI navigation
	vector<vect_cube_t> space_groups;
	invert_walls(place_area, walls_per_dim, space, nav_pad);
	resize_cubes_xy(space, nav_pad); // restore padding (under-over)
	partition_cubes_into_conn_groups(space, space_groups, pad);
	vect_cube_t small_rooms, door_keepout, unreachable_rooms;

#if 0 // TESTING
	rand_gen_t rgen2(rgen);
	for (vect_cube_t const &g : space_groups) {
		colorRGBA const color(rgen2.rand_float(), rgen2.rand_float(), rgen2.rand_float());
		for (cube_t c : g) {
			set_cube_zvals(c, zval, zval+0.5*floor_spacing);
			objs.emplace_back(c, TYPE_DBG_SHAPE, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_BACKROOM), 1.0, SHAPE_CUBE, color);
		}
	}
#endif
	// Add doorways + doors to guarantee full connectivity using space
	if (space_groups.size() > 1) { // multiple disconnected sub-graphs
		float const door_min_spacing(0.5*doorway_width), min_shared_edge(doorway_width + 2*wall_thickness); // allow space for door frame
		float const min_wall_spacing(get_trim_thickness());
		vector<group_range_t> adj;
		vect_cube_t doors_to_add, walls_to_add, all_doors;
		set<pair<unsigned, unsigned>> connected;
		vector<uint8_t> reachable(space_groups.size(), 0); // start as all unreachable
		bool const first_dim(rgen.rand_bool()); // mix it up to avoid favoring one dim

		for (unsigned d = 0; d < 2; ++d) {
			unsigned const dim(bool(d) ^ first_dim);
			walls_to_add.clear();
			vect_cube_t &walls(walls_per_dim[dim]);

			for (cube_t &wall : walls) {
				cube_t query(wall);
				query.expand_in_dim(dim, pad); // increase wall width to overlap space cubes
				adj.clear();

				for (unsigned gix = 0; gix < space_groups.size(); ++gix) {
					bool hit(0);

					for (cube_t const &s : space_groups[gix]) {
						if (!s.intersects(query)) continue;
						float const lo(max(query.d[!dim][0], s.d[!dim][0])), hi(min(query.d[!dim][1], s.d[!dim][1]));
						if (!hit) {adj.emplace_back(gix, lo, hi);} else {adj.back().update(lo, hi);}
						hit = 1;
					}
					if (hit && adj.back().len() < min_shared_edge) {adj.pop_back();} // remove if too short
				} // for gix
				if (adj.size() <= 1) continue; // not adjacent to multiple space group
				//objs.emplace_back(query, TYPE_DBG_SHAPE, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_BACKROOM), 1.0, SHAPE_CUBE, RED); // TESTING
				doors_to_add.clear();

				// test every pair of groups
				for (unsigned i = 0; i < adj.size(); ++i) {
					for (unsigned j = i+1; j < adj.size(); ++j) {
						group_range_t const &g1(adj[i]), &g2(adj[j]);
						float const lo(max(g1.lo, g2.lo)), hi(min(g1.hi, g2.hi));
						if (hi - lo <= min_shared_edge) continue; // not enough space to add a door
						pair<unsigned, unsigned> const ix_pair(min(g1.gix, g2.gix), max(g1.gix, g2.gix));
						if (connected.find(ix_pair) != connected.end()) continue; // these two groups were already connected
						float const vmin(lo + doorway_hwidth + min_wall_spacing), vmax(hi - doorway_hwidth - min_wall_spacing);
						if (vmin >= vmax) continue; // not enough space to add a door (should be rare)
							
						// add a doorway here if possible
						for (unsigned n = 0; n < 10; ++n) { // make up to 10 tries
							float const door_pos(rgen.rand_uniform(vmin, vmax));
							cube_t door_cand(wall);
							set_wall_width(door_cand, door_pos, doorway_hwidth, !dim);
							// check for collisions with walls in other dims and prev placed doors; shouldn't collide with walls in this dim due to checks during wall placement
							cube_t wall_test_cube(door_cand);
							wall_test_cube.expand_in_dim(dim, doorway_width); // make room for the door to open
							if (has_bcube_int(wall_test_cube, walls_per_dim[!dim])) continue;
							cube_t door_test_cube(wall);
							door_test_cube.expand_in_dim(!dim, door_min_spacing);
							if (has_bcube_int(door_test_cube, all_doors)) continue; // check all doors, even ones on other walls in case there are two overlapping doors
							if (room.has_stairs && interior->is_blocked_by_stairs_or_elevator(door_test_cube)) continue; // check for stairs; doesn't check open path
							doors_to_add.push_back(door_cand);
							all_doors   .push_back(door_cand);
							connected.insert(ix_pair); // mark these two space groups as being connected
							reachable[g1.gix] = reachable[g2.gix] = 1; // assume both of these are now reachable; likely, but not guaranteed
							break;
						} // for n
					} // for j
				} // for i
				sort(doors_to_add.begin(), doors_to_add.end(), cube_by_center_dim_descending(!dim)); // sort in dim !dim, high to low

				for (cube_t const &door : doors_to_add) {
					assert(door.is_strictly_normalized());
					// select an open direction that doesn't block another door
					bool open_dir(rgen.rand_bool());
					cube_t open_area(door);
					open_area.d[dim][open_dir] += (open_dir ? 1.0 : -1.0)*doorway_width;
					// swap the open direction if blocked by another door or stairs; hopefully we don't block a different door now
					if (has_bcube_int(open_area, door_keepout) || (room.has_stairs && interior->is_blocked_by_stairs_or_elevator(open_area))) {open_dir ^= 1;}
					// cut the doorway into the wall and add the door
					cube_t wall2;
					bool const make_unlocked = 1; // makes exploring easier
					bool const make_closed   = 1; // makes it easier to tell if the door has been used
					// this should add one door and one door stack
					assert(door.d[!dim][0] > wall.d[!dim][0] && door.d[!dim][1] < wall.d[!dim][1]);
					remove_section_from_cube_and_add_door(wall, wall2, door.d[!dim][0], door.d[!dim][1], !dim, open_dir, 0, make_unlocked, make_closed); // is_bathroom=0
					interior->door_stacks.back().set_backrooms();
					interior->doors      .back().set_backrooms();
					walls_to_add.push_back(wall2); // keep high side as it won't be used with any other doors
					// add a blocker so that no ceiling lights are placed in the path of this door
					cube_t blocker(door);
					blocker.d[dim][open_dir] += (open_dir ? 1.0 : -1.0)*doorway_width; // add clearance in front
					objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, 0, (RO_FLAG_INVIS | RO_FLAG_BACKROOM), tot_light_amt);
					// keep other doors from opening into this door's space
					cube_t keepout(door);
					keepout.expand_in_dim(dim, doorway_width); // don't block either end
					door_keepout.push_back(keepout);
				} // for door
			} // for wall
			vector_add_to(walls_to_add, walls);
		} // for d
		for (unsigned gix = 0; gix < space_groups.size(); ++gix) {
			vect_cube_t const &group(space_groups[gix]);

			if (group.size() == 1 && !reachable[gix]) {
				unreachable_rooms.push_back(group.front()); // track unreachable single cube rooms
				continue; // not counted as small rooms
			}
			if (group.size() != 1) {
				if (group.size() < space.size()/8) { // small area, but not a single rectangle; still may need to add room light
					float tot_area(0.0);
					cube_t group_bounds;

					for (cube_t const &c : group) {
						group_bounds.assign_or_union_with_cube(c);
						tot_area += c.get_area_xy();
					}
					if (tot_area > 0.5*group_bounds.get_area_xy() && tot_area < 0.1*room.get_area_xy()) {rooms_to_light.push_back(group_bounds);}
				}
				continue;
			}
			cube_t sub_room(group.front());
			sub_room.intersect_with_cube(true_room); // can't go outside the backrooms (under-over can move exterior walls)

			for (auto const &c : connected) { // add if it was connected with a door
				if (c.first == gix || c.second == gix) {small_rooms.push_back(sub_room); break;}
			}
		} // for gix
	}

	// Add some random pillars in large open spaces
	float const pillar_grid_step(2.5*min_gap), pillar_grid_exp(0.75*pillar_grid_step), pillar_hwidth(0.07*floor_spacing), merge_wall_min(1.0*pillar_grid_step);
	unsigned const xdiv(ceil(wall_area.dx()/pillar_grid_step)), ydiv(ceil(wall_area.dy()/pillar_grid_step));
	float const xstep(wall_area.dx()/(xdiv-1)), ystep(wall_area.dy()/(ydiv-1));
	vect_cube_t big_space;
	cube_t grid(wall_area); // copy zvals

	for (unsigned y = 0; y < ydiv; ++y) {
		for (unsigned x = 0; x < xdiv; ++x) {
			set_wall_width(grid, (wall_area.y1() + y*ystep), pillar_grid_exp, 1);
			set_wall_width(grid, (wall_area.x1() + x*xstep), pillar_grid_exp, 0);	
			if (!place_area.contains_cube_xy(grid)) continue;
			if (has_bcube_int(grid, walls_per_dim[0]) || has_bcube_int(grid, walls_per_dim[1]) || has_bcube_int(grid, pillar_avoid)) continue;
			big_space.push_back(grid);
		} // for x
	} // for y
	partition_cubes_into_conn_groups(big_space, space_groups, pad);
	cube_t pillar(wall_area); // copy zvals
	room_obj_shape const pillar_shape(rgen.rand_bool() ? SHAPE_CYLIN : SHAPE_CUBE); // randomly cube or cylinder per-floor; should this be per building?

	for (vect_cube_t &group : space_groups) {
		if (rgen.rand_float() < 0.4) continue; // no pillars 40% of the time
		cube_t prev_space, cur_row;
		sort(group.begin(), group.end(), cube_by_xy_dim(rgen.rand_bool())); // make adjacent space cubes adjacent in order, with dim chosen randomly

		for (auto s = group.begin(); s != group.end(); ++s) {
			for (unsigned d = 0; d < 2; ++d) {set_wall_width(pillar, s->get_center_dim(d), pillar_hwidth, d);}
			objs.emplace_back(pillar, TYPE_PG_PILLAR, room_id, 0, 0, RO_FLAG_BACKROOM, tot_light_amt, pillar_shape, wall_color); // dim=0, dir=0
			// maybe merge rows of adjacent pillars into a single wall
			cube_t merge_cand;

			if (!s->intersects(prev_space) || ((pillar.x1() != cur_row.x1() || pillar.x2() != cur_row.x2()) && (pillar.y1() != cur_row.y1() || pillar.y2() != cur_row.y2()))) {
				// different row, gap, or first pillar
				merge_cand = cur_row;
				cur_row    = pillar; // seed for next row
			}
			else { // extend row
				cur_row.union_with_cube(pillar);
			}
			if (s+1 == group.end()) {merge_cand = cur_row;} // last pillar
			prev_space = *s;

			if (!merge_cand.is_all_zeros() && rgen.rand_float() < 0.5) { // add wall 50% of the time
				for (unsigned dim = 0; dim < 2; ++dim) {
					if (merge_cand.get_sz_dim(!dim) < merge_wall_min) continue; // too short to create a wall; will get here in at least one dim
					if (has_bcube_int(merge_cand, pillar_avoid))      continue;
					cube_t wall(merge_cand);
					wall.expand_by_xy(wall_half_thick - pillar_hwidth); // shrink
					walls_per_dim[dim].push_back(wall);
				}
			}
		} // for d
	} // for g

	// add walls
	auto pbgr_walls(interior->room_geom->pgbr_walls);
	auto &pgbr_wall_ixs(interior->room_geom->pgbr_wall_ixs);
	if (pgbr_wall_ixs.empty()) {pgbr_wall_ixs.emplace_back(pbgr_walls);} // add first index

	for (unsigned d = 0; d < 2; ++d) {
		sort(walls_per_dim[d].begin(), walls_per_dim[d].begin(), cube_by_sz(!d)); // sort walls longest to shortest to improve occlusion culling time
		
		for (cube_t &wall : walls_per_dim[d]) {
			objs.emplace_back(wall, TYPE_PG_WALL, room_id, d, 0, RO_FLAG_BACKROOM, tot_light_amt, SHAPE_CUBE, wall_color); // dir=0
			pbgr_walls[d].push_back(wall); // store walls for occlusion and door opening checks
#if 0 // enable for pedestrian navigation debugging
			cube_t c(wall);
			c.z2() -= 0.5*wall.dz(); // half height
			c.expand_by_xy(get_ped_coll_radius());
			objs.emplace_back(c, TYPE_DBG_SHAPE, room_id, d, 0, (RO_FLAG_NOCOLL | RO_FLAG_BACKROOM), 1.0, SHAPE_CUBE, RED);
#endif
		}
	} // for d
	pgbr_wall_ixs.emplace_back(pbgr_walls); // end of range
	
	// Add vents, light switch, and outlets (on all walls - or should it be only on the wall adjacent to building or next to the door?)
	add_light_switches_to_room(rgen, true_room, zval, room_id, objs_start, 0, 1); // is_ground_floor=0, is_basement=1
	add_outlets_to_room       (rgen, true_room, zval, room_id, objs_start, 0, 1); // is_ground_floor=0, is_basement=1
	add_wall_vent_to_room     (rgen, true_room, zval, room_id, objs_start, 0   ); // check_for_ducts=0
	rgen.rand_bool(); // mix up in between
	add_ceil_vent_to_room     (rgen, true_room, zval, room_id, objs_start); // add a ceiling vent as well
	
	// Make small rooms with doors bathrooms, etc.
	for (cube_t const &r : small_rooms) {
		unsigned num_doors(0);

		for (cube_t const &dk : door_keepout) {
			if (dk.intersects_xy(r)) {++num_doors;}
		}
		if (interior->get_ext_basement_door().get_true_bcube().intersects(r)) {++num_doors;} // backrooms entrance door counts
		if (num_doors == 0) continue; // not connected with a door, not reachable, skip (and don't need a light either)
		rooms_to_light.push_back(r);
		room_t sub_room(room, r); // keep flags, copy cube
		sub_room.interior = 1; // treated as basement but not extended basement (no wall padding)
		sub_room.intersect_with_cube(place_area); // clip off areas adjacent to the exterior wall
		unsigned const sub_objs_start(objs.size()); // no objects have been placed in this sub-room yet
		bool objs_added(0), no_boxes(0);

		if (num_doors == 1) { // only make bathroom if there's a single door
			float floor_zval(zval); // may be modified below, but otherwise unused
			unsigned const floor_ix(0); // pass this in, or always zero?
			unsigned added_bathroom_objs_mask(0); // unused
			objs_added = no_boxes = add_bathroom_objs(rgen, sub_room, floor_zval, room_id, tot_light_amt, sub_objs_start, floor_ix, 1, 0, added_bathroom_objs_mask); // is_basement=1
			if (sub_room.get_has_mirror()) {room.set_has_mirror();}
		}
		// 2 or more rooms
		if (!objs_added) { // maybe make a machine room
			objs_added = add_machines_to_room(rgen, sub_room, zval, room_id, tot_light_amt, objs_start);
		}
		if (!objs_added && rgen.rand_bool()) { // add furnace
			objs_added = add_furnace_to_room(rgen, sub_room, zval, room_id, tot_light_amt, sub_objs_start);
		}
		if (!objs_added) { // add water heater
			objs_added = add_water_heaters(rgen, sub_room, zval, room_id, tot_light_amt, sub_objs_start, 1); // single_only=1
		}
		if (!objs_added) {} // try other objects or room types?
		if (!/*objs_added*/no_boxes) { // add some random boxes if we haven't added anything else/if not a bathroom
			unsigned const max_num_boxes(rgen.rand() % 4); // 0-3
			add_boxes_to_room(rgen, sub_room, zval, room_id, tot_light_amt, objs_start, max_num_boxes);
		}
	} // for r

	// maybe add a fire extinguisher to a random wall
	float fe_height(0.0), fe_radius(0.0);

	if (get_fire_ext_height_and_radius(floor_spacing, fe_height, fe_radius)) { // add a fire extinguisher on the wall
		unsigned const num_fire_ext(1 + (rgen.rand() % 3)); // 1-3

		for (unsigned N = 0; N < num_fire_ext; ++N) {
			for (unsigned n = 0; n < 10; ++n) { // make 10 tries
				bool const dim(rgen.rand_bool()); // choose a random wall dim
				vect_cube_t const &walls(walls_per_dim[dim]);
				if (walls.empty()) continue;
				cube_t const &wall(walls[rgen.rand() % walls.size()]);
				float const min_clearance(2.0*fe_radius), wall_pos_lo(wall.d[!dim][0] + min_clearance), wall_pos_hi(wall.d[!dim][1] - min_clearance);
				if (wall_pos_lo >= wall_pos_hi) continue; // wall too short
				bool const dir(rgen.rand_bool()); // random, but the same across all floors
				float const wall_edge(wall.d[dim][!dir]), val(rgen.rand_uniform(wall_pos_lo, wall_pos_hi));
				cube_t bounds(wall);
				set_wall_width(bounds, val, fe_radius, !dim);
				bounds.d[dim][!dir] = wall_edge + (dir ? -1.0 : 1.0)*0.1*wall_thickness; // move slightly away from wall to avoid colliding with it
				bounds.d[dim][ dir] = wall_edge + (dir ? -1.0 : 1.0)*2.0*fe_radius; // outer edge
				if (overlaps_other_room_obj(bounds, objs_start) || is_obj_placement_blocked(bounds, true_room, 1, 0)) continue; // inc_open_doors=1, check_open_dir=0
				if (has_bcube_int(bounds, unreachable_rooms)) continue; // don't place in an unreachable room
				add_fire_ext(fe_height, fe_radius, zval, wall_edge, val, room_id, tot_light_amt, dim, dir);
				break; // done
			} // for n
		} // for N
	}

	// Add occasional random items/furniture: chairs, office chairs, boxes, crates, balls, etc.
	unsigned const num_chairs(rgen.rand() % 5); // 0-4

	if (num_chairs > 0) { // add chairs
		// note that chairs don't block the player (and can be taken/moved), but they may block the AI
		float const max_radius(0.25*floor_spacing); // should be at least as large as the max chair size
		colorRGBA chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
		vect_cube_t blockers;

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (!i->no_coll()) {blockers.push_back(*i);} // use the simple no_coll() check because it works with the type of objects placed in this room
		}
		for (unsigned n = 0; n< num_chairs; ++n) {
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()), office_chair(rgen.rand_bool()); // 50% chace of office chair with random_rotation=1

			for (unsigned N = 0; N < 4; ++N) { // make up to 4 attempts to place a chair
				point const chair_pos(gen_xy_pos_in_area(place_area, max_radius, rgen, zval));
				if (check_vect_cube_contains_pt_xy(unreachable_rooms, chair_pos)) continue; // don't place in an unreachable room
				
				if (add_chair(rgen, true_room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, office_chair, 1)) {
					objs.back().flags |= RO_FLAG_BACKROOM; // flag so that player collisions are enabled
					blockers.push_back(objs.back());
					break;
				}
			}
		} // for n
	}
	// note that boxes/crates are only placed along the exterior walls; this helps prevent the AI from getting stuck on them, since it can't take boxes like the player
	unsigned const num_boxes(rgen.rand() % 21); // 0-20
	add_boxes_to_room(rgen, true_room, zval, room_id, tot_light_amt, objs_start, num_boxes);
	unsigned const num_balls(rgen.rand() % 4); // 0-3
	for (unsigned n = 0; n < num_balls; ++n) {add_ball_to_room(rgen, true_room, place_area, zval, room_id, tot_light_amt, objs_start);}
}

void building_t::add_missing_backrooms_lights(rand_gen_t rgen, float zval, unsigned room_id, unsigned objs_start, unsigned lights_start,
	room_object_t const &ref_light, vect_cube_t const &rooms_to_light, light_ix_assign_t &light_ix_assign)
{
	vect_room_object_t &objs(interior->room_geom->objs);
	room_t const &room(get_room(room_id));
	// apply custom colors to some lights; at the moment this only includes grid lights and not small room lights added below
	unsigned const NUM_LIGHT_COLORS = 5;
	colorRGBA const light_colors[NUM_LIGHT_COLORS] = {RED, GREEN, BLUE, YELLOW, ORANGE};
	float const zone_spacing(8.0*get_window_vspace());
	unsigned const zone_x_add(rgen.rand()), zone_y_add(rgen.rand());
	rand_gen_t zone_rgen;

	for (auto i = objs.begin()+lights_start; i != objs.end(); ++i) {
		assert(i->type == TYPE_LIGHT);
		unsigned const zone_x(zone_x_add + round_fp(i->xc()/zone_spacing)), zone_y(zone_y_add + round_fp(i->yc()/zone_spacing));
		zone_rgen.set_state(zone_x, zone_y);
		zone_rgen.rand_mix();
		if (zone_rgen.rand_float() > 0.2) continue; // 20% chance of colored light
		i->color = light_colors[zone_rgen.rand() % NUM_LIGHT_COLORS];
	} // for i
	unsigned const lights_end(objs.size());
	point const ref_light_center(ref_light.get_cube_center());

	for (cube_t const &r : rooms_to_light) {
		bool has_light(0);

		for (auto i = objs.begin()+lights_start; i != objs.begin()+lights_end; ++i) {
			assert(i->type == TYPE_LIGHT);
			if (i->intersects(r)) {has_light = 1; break;}
		}
		if (has_light) continue;
		room_object_t light(ref_light); // what if this is a short light that was blocked by a doorway? use a different light?
		light += vector3d((r.xc() - ref_light_center.x), (r.yc() - ref_light_center.y), 0.0);
		bool const room_dim(r.dx() < r.dy()); // longer room dim
		add_sub_room_light(light, room_t(room, r), room_dim, objs_start, light_ix_assign, rgen);
	} // for r
}
void building_t::add_sub_room_light(room_object_t light, room_t const &room, bool dim, unsigned objs_start, light_ix_assign_t &light_ix_assign, rand_gen_t &rgen) {
	vect_cube_t to_add;
	try_place_light_on_ceiling(light, room, dim, get_fc_thickness(), 1, 0, 1, 1, objs_start, to_add, rgen); // or wall light?

	for (cube_t const &L : to_add) { // should be size 1
		light.copy_from(L);
		light.obj_id = light_ix_assign.get_ix_for_light(light);
		interior->room_geom->objs.push_back(light);
	}
}

unsigned building_t::setup_multi_floor_room(extb_room_t &room, door_t const &door, bool wall_dim, bool wall_dir, rand_gen_t &rgen) { // for backrooms
	if (!interior) return 1; // shouldn't call?
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness), wall_thickness(get_wall_thickness());
	unsigned const num_floors(calc_num_floors(room, floor_spacing, floor_thickness));
	assert(num_floors > 0);
	if (num_floors == 1) return 1;
	assert(interior->fc_occluders.empty()); // must be called before setup_fc_occluders() because it adds to floors and ceilings
	// add wall segment under the door on lower floors
	float const room_entrance_edge(room.d[wall_dim][!wall_dir]), dir_sign(wall_dir ? 1.0 : -1.0);
	cube_t wall(door);
	set_cube_zvals(wall, room.z1()+fc_thick, door.z1());
	wall.d[wall_dim][!wall_dir] = room_entrance_edge + 0.1*dir_sign*wall_thickness; // shift out slightly
	wall.d[wall_dim][ wall_dir] = room_entrance_edge + 1.0*dir_sign*wall_thickness;
	assert(wall.is_strictly_normalized());
	interior->walls[wall_dim].push_back(wall);

	cube_t door_avoid(door.get_clearance_bcube());
	door_avoid.union_with_cube(door.get_open_door_path_bcube()); // make sure it's path is clear as well
	vect_cube_t avoid, cf_to_add;
	avoid.push_back(door_avoid);
	//for (door_stack_t const &ds : interior->door_stacks) {if (ds.intersects(room)) {avoid.push_back(ds.get_open_door_path_bcube());}} // doors not yet placed
	stairs_shape const sshape(SHAPE_STRAIGHT);
	float const door_width(door.get_width()), stairs_hwidth(0.6*door_width), stairs_hlen(rgen.rand_uniform(2.0, 3.0)*stairs_hwidth), wall_spacing(1.0*door_width);
	float z(room.z1() + floor_spacing); // move to next floor

	for (unsigned f = 1; f < num_floors; ++f, z += floor_spacing) { // skip first floor - draw pairs of floors and ceilings
		float const zc(z - fc_thick), zf(z + fc_thick);
		cf_to_add.clear();
		cf_to_add.push_back(room); // seed with entire room

		// add stairs
		if (4.0*(stairs_hlen + wall_spacing) < min(room.dx(), room.dy())) { // condition should generally be true
			unsigned const num_stairs(1 + (rgen.rand()&1)); // 1-2 stairs per floor

			for (unsigned n = 0; n < num_stairs; ++n) {
				bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
				vector3d size;
				size[ dim] = stairs_hlen;
				size[!dim] = stairs_hwidth;
				cube_t place_area(room);
				for (unsigned d = 0; d < 2; ++d) {place_area.expand_in_dim(d, -(size[d] + wall_spacing));}

				for (unsigned n = 0; n < 100; ++n) { // 100 iterations to place stairs
					cube_t cand;
					for (unsigned d = 0; d < 2; ++d) {set_wall_width(cand, rgen.rand_uniform(place_area.d[d][0], place_area.d[d][1]), size[d], d);}
					set_cube_zvals(cand, zf-floor_spacing, zc+floor_spacing); // top of floor below to bottom of ceiling on this floor
					if (has_bcube_int(cand, avoid)) continue; // bad placement
					assert(cand.is_strictly_normalized());
					subtract_cube_from_cubes(cand, cf_to_add);
					landing_t landing(cand, 0, 0, dim, dir, 1, sshape, 0, 1, 1, 0, 1); // elevator=0, floor=0, railing=1, roof=0, at_top=1, stacked=1, ramp=0, stacked=1, extb=1
					set_cube_zvals(landing, zc, zf);
					interior->landings.push_back(landing);
					stairwell_t const S(cand, 1, dim, dir, sshape, 0, 1, 1); // floors=1, roof=0, stacked=1, ext_basement=1
					interior->stairwells.push_back(S);
					avoid.push_back(get_stairs_bcube_expanded(S, door_width, wall_thickness, door_width));
					room.has_stairs = 1;
					break; // success
				}
			} // for n
		}
		for (cube_t &cf : cf_to_add) {interior->add_ceil_floor_pair(cf, zc, z, zf);}
	} // for f
	return num_floors;
}

bool building_t::is_pos_in_pg_or_backrooms(point const &pos) const { // Note: typically used to guard interior->room_geom->pgbr_walls logic
	return (has_parking_garage && pos.z < ground_floor_z1 && ((interior && interior->has_backrooms) || get_basement().contains_pt(pos)));
}

bool building_room_geom_t::cube_int_backrooms_walls(cube_t const &c) const { // used for door opening collision checks
	// no dim is passed in, so we check both dims; includes parking garage walls, which we can probably ignore, but they should be small in size
	// could accelerate with interior->room_geom->pgbr_wall_ixs, but the extra complexity may no be necessary
	for (unsigned d = 0; d < 2; ++d) {
		if (has_bcube_int(c, pgbr_walls[d])) return 1;
	}
	return 0;
}


