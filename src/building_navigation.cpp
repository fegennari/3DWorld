// 3D World - Building Navigation for AIs
// by Frank Gennari 3/14/20

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for pedestrian_t
#include <queue>
#pragma warning(disable : 26812) // prefer enum class over enum


bool const STAY_ON_ONE_FLOOR = 0; // Note: multi-floor movement not yet fully supported

extern float fticks;

point get_cube_center_zval(cube_t const &c, float zval) {return point(c.xc(), c.yc(), zval);}

// Note: this should go into building_t/buildings.h at some point, but is temporarily here
class building_nav_graph_t {
	struct pt_with_ix_t { // size=12
		unsigned ix;
		vector2d pt[2]; // stairs store {up, down}
		pt_with_ix_t(unsigned ix_, vector2d const &ptu, vector2d const &ptd) : ix(ix_) {pt[0] = ptu; pt[1] = ptd;}
	};
	struct node_t { // represents one room or one stairwell
		bool has_exit, is_hallway, is_stairs; // has_exit and is_stairs are not yet used
		cube_t bcube;
		vector<pt_with_ix_t> conn_rooms;
		node_t() : has_exit(0), is_hallway(0), is_stairs(0) {}
		point get_center(float zval) const {return get_cube_center_zval(bcube, zval);}

		void add_conn_room(unsigned room, cube_t const &cu, cube_t const &cd) {
			for (auto i = conn_rooms.begin(); i != conn_rooms.end(); ++i) {if (i->ix == room) return;} // ignore duplicates
			conn_rooms.emplace_back(room, vector2d(cu.xc(), cu.yc()), vector2d(cd.xc(), cd.yc()));
		}
	};
	struct a_star_node_state_t {
		int came_from_ix;
		point path_pt;
		float g_score, h_score, f_score;
		a_star_node_state_t() : came_from_ix(-1), g_score(0), h_score(0), f_score(0) {}
	};

	unsigned num_rooms, num_stairs;
	float stairs_extend;
	vector<node_t> nodes;
	node_t       &get_node(unsigned room)       {assert(room < nodes.size()); return nodes[room];}
	node_t const &get_node(unsigned room) const {assert(room < nodes.size()); return nodes[room];}

	void remove_connection(unsigned from, unsigned to) {
		auto &conn(get_node(from).conn_rooms);

		for (auto i = conn.begin(); i != conn.end(); ++i) {
			if (i->ix == to) {swap(*i, conn.back()); conn.pop_back(); return;}
		}
		assert(0); // must be found - should not get here
	}
public:
	building_nav_graph_t(float stairs_extend_) : num_rooms(0), num_stairs(0), stairs_extend(stairs_extend_) {}

	void set_num_rooms(unsigned num_rooms_, unsigned num_stairs_) {
		num_rooms  = num_rooms_;
		num_stairs = num_stairs_;
		nodes.resize(num_rooms + num_stairs);
		for (unsigned n = num_rooms; n < nodes.size(); ++n) {nodes[n].is_stairs = 1;}
	}
	void set_room_bcube  (unsigned room,   cube_t const &c) {get_node(room).bcube = c;}
	void set_stairs_bcube(unsigned stairs, cube_t const &c) {get_node(stairs + num_rooms).bcube = c;}
	void mark_hallway(unsigned room) {get_node(room).is_hallway = 1;}
	void mark_exit   (unsigned room) {get_node(room).has_exit   = 1;}

	void connect_stairs(unsigned room, unsigned stairs, bool dim, bool dir) {
		assert(room < num_rooms && stairs < num_stairs);
		unsigned const node_ix2(num_rooms + stairs);
		node_t &n2(get_node(node_ix2));
		cube_t entry_u(n2.bcube), entry_d(n2.bcube);
		float const extend((dir ? -1.0 : 1.0)*stairs_extend); // extend away from stairs for entrance/exit area; will be denormalized in this dim
		entry_u.d[dim][ dir] = entry_u.d[dim][!dir] + extend; // shrink to zero area at the entrance to the stairs when going up
		entry_d.d[dim][!dir] = entry_d.d[dim][ dir] - extend; // shrink to zero area at the entrance to the stairs when going down
		get_node(room).add_conn_room(node_ix2, entry_u, entry_d);
		n2.add_conn_room(room, entry_u, entry_d);
	}
	void connect_rooms(unsigned room1, unsigned room2, cube_t const &conn_bcube) { // graph is bidirectional
		assert(room1 < num_rooms && room2 < num_rooms);
		get_node(room1).add_conn_room(room2, conn_bcube, conn_bcube);
		get_node(room2).add_conn_room(room1, conn_bcube, conn_bcube);
	}
	void disconnect_room_pair(unsigned room1, unsigned room2) { // remove connections in both directions
		assert(room1 != room2 && room1 < num_rooms && room2 < num_rooms);
		remove_connection(room1, room2);
		remove_connection(room2, room1);
	}
	bool is_room_connected_to(unsigned room1, unsigned room2) const { // Note: likely faster than running full A* algorithm
		assert(room1 < num_rooms && room2 < num_rooms);
		if (room1 == room2) return 1;
		vector<uint8_t> seen(nodes.size(), 0);
		vector<unsigned> pend;
		pend.push_back(room1);
		seen[room1] = 1;

		while (!pend.empty()) {
			node_t const &node(get_node(pend.back()));
			pend.pop_back();

			for (auto i = node.conn_rooms.begin(); i != node.conn_rooms.end(); ++i) {
				if (i->ix == room2) return 1; // found, done
				if (!seen[i->ix]) {pend.push_back(i->ix); seen[i->ix] = 1;}
			}
		} // end while()
		return 0; // not found
	}
	unsigned count_connected_components() const {
		if (nodes.empty()) return 0;
		vector<uint8_t> seen(nodes.size(), 0);
		vector<unsigned> pend;
		unsigned ncomp(0);
		
		for (unsigned n = 0; n < nodes.size(); ++n) {
			if (seen[n]) continue; // node already processed
			++ncomp; // start a new component
			pend.push_back(n);
			seen[n] = 1;

			while (!pend.empty()) {
				node_t const &node(get_node(pend.back()));
				pend.pop_back();

				for (auto i = node.conn_rooms.begin(); i != node.conn_rooms.end(); ++i) {
					if (!seen[i->ix]) {pend.push_back(i->ix); seen[i->ix] = 1;}
				}
			} // end while()
		} // for n
		return ncomp;
	}
	bool is_fully_connected() const {return (count_connected_components() == 1);}

	bool is_valid_pos(vect_cube_t const &avoid, point const &pos, float radius, float height) const { // Note: assumes zvals are already checked
		cube_t c(pos, pos);
		c.expand_by_xy(radius);

		for (auto i = avoid.begin(); i != avoid.end(); ++i) {
			if (i->intersects_xy(c) && sphere_cube_intersect_xy(pos, radius, *i)) return 0;
		}
		return 1;
	}
	static void choose_pt_xy_in_room(point &pos, cube_t const &c, rand_gen_t &rgen) {
		for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(c.d[d][0], c.d[d][1]);}
	}
	point get_stairs_entrance_pt(float zval, unsigned node_ix, bool up_or_down) const {
		node_t const &node(get_node(node_ix));
		assert(!node.conn_rooms.empty());
		vector2d const &pt(node.conn_rooms.front().pt[up_or_down]); // Note: all conn_rooms should be the same value
		return point(pt.x, pt.y, zval);
	}
	point find_valid_room_dest(vect_cube_t const &avoid, float radius, float height, float zval, unsigned node_ix, bool up_or_down, bool &not_room_center, rand_gen_t &rgen) const {
		node_t const &node(get_node(node_ix));
		if (node.is_stairs) {return get_stairs_entrance_pt(zval, node_ix, up_or_down);}
		cube_t place_area(node.bcube);
		place_area.expand_by_xy(-2.0*radius); // shrink by twice the radius
		point const center(get_cube_center_zval(node.bcube, zval));
		if (!place_area.is_strictly_normalized()) return center; // should generally not be true
		not_room_center = 1;
		point pos(center); // first candidate is the center of the room

		for (unsigned n = 0; n < 100; ++n) { // 100 random tries to find a valid dest_pos
			if (is_valid_pos(avoid, pos, radius, height)) return pos; // success
			choose_pt_xy_in_room(pos, place_area, rgen); // choose a random new point in the room
		}
		return center; // failed, return room center
	}

	static bool check_line_int_xy(vect_cube_t const &c, point const &p1, point const &p2) {
		for (auto i = c.begin(); i != c.end(); ++i) {
			if (check_line_clip_xy(p1, p2, i->d)) return 1;
		}
		return 0;
	}
	static bool check_pt_contained_xy(vect_cube_t const &c, point const &p) {
		for (auto i = c.begin(); i != c.end(); ++i) {
			if (i->contains_pt_xy(p)) return 1;
		}
		return 0;
	}
	bool connect_room_endpoints(vect_cube_t const &avoid, cube_t const &walk_area, point const &p1, point const &p2, float radius,
		vector<point> &path, vect_cube_t &keepout, rand_gen_t &rgen, bool ignore_p1_coll=0, bool ignore_p2_coll=0) const
	{
		assert(p1.z == p2.z);
		bool is_path_valid(1);
		keepout.clear();

		for (auto i = avoid.begin(); i != avoid.end(); ++i) {
			if (!i->intersects_xy(walk_area)) continue;
			cube_t c(*i);
			c.expand_by_xy(radius);

			if (c.contains_pt_xy(p1)) {
				if (ignore_p1_coll) continue;
				return 0;
			}
			if (c.contains_pt_xy(p2)) {
				if (ignore_p2_coll) continue;
				return 0;
			}
			if (check_line_clip_xy(p1, p2, c.d)) {is_path_valid = 0;}
			keepout.push_back(c);
		}
		if (is_path_valid) return 1; // done

		// straight line not valid, need to run path finding
		assert(!keepout.empty());
		// do something simple and brute force:
		// since rooms tend to only have objects in the middle, try to find the best point that creates the shortest non-intersecting 2-part or 3-part path
		for (unsigned npts = 1; npts <= 2; ++npts) { // try to add 1 point; if that fails, try to add 2 points
			unsigned const num_tries((npts == 1) ? 200 : 100);
			float dmin(0.0);
			bool use_pos2(0);
			point best_pt, best_pt2, pos, pos2;
			pos.z = pos2.z = p1.z;

			for (unsigned n = 0; n < num_tries; ++n) { // make num_tries attempts
				choose_pt_xy_in_room(pos, walk_area, rgen); // choose a rand point in the room
				if (check_pt_contained_xy(keepout, pos)) continue; // bad point
				bool const seg1_bad(check_line_int_xy(keepout, p1, pos)), seg2_bad(check_line_int_xy(keepout, pos, p2));

				if (npts == 1 || !(seg1_bad || seg2_bad)) { // single point or both segments valid
					if (seg1_bad || seg2_bad) continue; // bad point (either line intersects)
					float const dist(p2p_dist_xy(p1, pos) + p2p_dist_xy(pos, p2));
					if (dmin == 0.0 || dist < dmin) {best_pt = pos; dmin = dist; use_pos2 = 0;}
				}
				else { // npts == 2
					if (seg1_bad && seg2_bad) continue; // bad point (both lines intersect)
					bool success(0);

					for (unsigned m = 0; m < 20; ++m) { // make 20 random attempts
						choose_pt_xy_in_room(pos2, walk_area, rgen); // choose a rand point in the room
						if (!check_pt_contained_xy(keepout, pos2) && !check_line_int_xy(keepout, pos, pos2)) {success = 1; break;} // good point
					}
					if (!success) continue; // no good point

					for (unsigned ordering = 0; ordering < 2; ++ordering) { // try {p1-pos-pos2-p2} and {p1-pos2-pos-p2}
						if (check_line_int_xy(keepout, pos, p1) || check_line_int_xy(keepout, pos2, p2)) {swap(pos, pos2); continue;} // bad point (could optimize this if needed)
						float const dist(p2p_dist_xy(p1, pos) + p2p_dist_xy(pos, pos2) + p2p_dist_xy(pos2, p2));
						if (dmin == 0.0 || dist < dmin) {best_pt = pos; best_pt2 = pos2; dmin = dist; use_pos2 = 1;}
						break;
					} // for ordering
				}
			} // for n
			if (dmin == 0.0) continue;
			path.push_back(best_pt);
			if (use_pos2) {path.push_back(best_pt2);}
			return 1; // success
		} // for npts
		return 0; // failed
	}
	static point closest_room_pt(cube_t const &c, point const &pos) {
		return point(max(c.x1(), min(c.x2(), pos.x)), max(c.y1(), min(c.y2(), pos.y)), pos.z);
	}
	bool reconstruct_path(vector<a_star_node_state_t> const &state, vect_cube_t const &avoid, point const &cur_pt,
		float radius, float height, unsigned start_ix, unsigned end_ix, bool is_first_path, bool up_or_down, vector<point> &path) const
	{
		unsigned n(start_ix);
		rand_gen_t rgen;
		static unsigned call_ix(1);
		rgen.set_state(start_ix, call_ix++);
		vect_cube_t keepout;

		while (1) {
			node_t const &node(get_node(n));
			point const &next(state[n].path_pt);
			int const came_from(state[n].came_from_ix);
			bool const is_first_pt(path.empty());
			cube_t walk_area(node.bcube);
			if (!node.is_stairs) {walk_area.expand_by_xy(-radius);} // shrink by radius
			
			if (node.is_hallway) {
				bool const min_dim(walk_area.dy() < walk_area.dx());
				float const shrink(min(0.5f*radius, max(0.0f, 0.9f*0.5f*walk_area.get_sz_dim(min_dim)))); // shrink by an extra half radius if hallway is wide enough
				walk_area.d[min_dim][0] += shrink; walk_area.d[min_dim][1] -= shrink;
			}
			if (is_first_pt) { // last point in path (first point in reverse path)
				assert(came_from >= 0);
				bool success(0);

				for (unsigned n = 0; n < 10; ++n) { // keep retrying until we find a point that is reachable from the doorway
					bool not_room_center(0);
					point const end_point(find_valid_room_dest(avoid, radius, height, cur_pt.z, start_ix, up_or_down, not_room_center, rgen));
					path.push_back(end_point);
					if (node.is_stairs) {success = 1; break;} // done, don't need to run code below
					point const room_exit(closest_room_pt(walk_area, next)); // first doorway
					if (connect_room_endpoints(avoid, walk_area, end_point, room_exit, radius, path, keepout, rgen)) {path.push_back(room_exit); success = 1; break;}
					path.clear(); // failed, reset for next iteration
					if (!not_room_center) break; // if we did choose the room center, and there is no path to it, we've failed
				} // for n
				if (!success) {assert(path.empty()); return 0;} // failed to connect to a point in dest room
			}
			else if (came_from < 0) { // done (next is not valid here)
				assert(n == end_ix);
				if (node.is_stairs) return 1; // success
				point const final_pt(closest_room_pt(walk_area, path.back())); // walk from room into last doorway
				path.push_back(final_pt);

				// find path to first doorway; ignore collisions with p2 (cur_pt) in case this person was pushed into an object by another person
				if (!connect_room_endpoints(avoid, walk_area, final_pt, cur_pt, radius, path, keepout, rgen, 0, 1)) {
					// allow a failure if this is the first path taken by this AI so that it's not stuck behind an object due to bad initial placement
					if (!is_first_path) {path.clear(); return 0;}
				}
				return 1; // success
			}
			else if (!node.is_stairs) { // adjust the path through a room
				assert(!path.empty());
				point const &prev(path.back());
				assert(prev.z == next.z);
				point const p1(closest_room_pt(walk_area, prev)), p2(closest_room_pt(walk_area, next));
				path.push_back(p1); // walk out of doorway and into room
				
				if (!connect_room_endpoints(avoid, walk_area, p1, p2, radius, path, keepout, rgen)) { // unreachable
					path.clear();
					// TODO: try another path?
					//disconnect_room_pair(n, came_from); // ???
					return 0;
				}
				path.push_back(p2); // walk from room into doorway
			}
			path.push_back(next); // doorway
			n = came_from;
		} // end while()
		return 0; // never gets here
	}
	
	// A* algorithm; Note: path is stored backwards
	bool find_path_points(unsigned room1, unsigned room2, float radius, float height, bool use_stairs,
		bool is_first_path, bool up_or_down, vect_cube_t const &avoid, point const &cur_pt, vector<point> &path) const
	{
		assert(room1 < nodes.size() && room2 < nodes.size());
		assert(room1 != room2); // or just return an empty path?
		path.clear();
		vector<a_star_node_state_t> state(nodes.size());
		vector<uint8_t> open(nodes.size(), 0), closed(nodes.size(), 0); // tentative/already evaluated nodes
		std::priority_queue<pair<float, unsigned> > open_queue;
		point const dest_pos(get_node(room2).get_center(cur_pt.z)); // Note: approximate, actual dest may be different
		a_star_node_state_t &start(state[room1]);
		start.g_score = 0.0;
		start.h_score = start.f_score = p2p_dist_xy(get_node(room1).get_center(cur_pt.z), dest_pos); // estimated total cost from start to goal through current
		open[room1]   = 1;
		open_queue.push(make_pair(-start.f_score, room1));

		while (!open_queue.empty()) {
			unsigned const cur(open_queue.top().second);
			open_queue.pop();
			assert(!closed[cur]);
			a_star_node_state_t const &cs(state[cur]);
			node_t const &cur_node(get_node(cur));
			point const center(cur_node.get_center(cur_pt.z));
			assert(!closed[cur]);
			closed[cur] = 1;
			open[cur]   = 0;

			for (auto i = cur_node.conn_rooms.begin(); i != cur_node.conn_rooms.end(); ++i) {
				assert(i->ix < nodes.size());
				if (closed[i->ix]) continue; // already closed (duplicate)
				node_t const &conn_node(get_node(i->ix));
				if (conn_node.is_stairs && !use_stairs && i->ix != room2) continue; // skip stairs in this mode
				point const conn_center(conn_node.get_center(cur_pt.z));
				a_star_node_state_t &sn(state[i->ix]);
				vector2d const &pt(i->pt[up_or_down]);
				float const new_g_score(sn.g_score + p2p_dist_xy(center, pt) + p2p_dist_xy(pt, conn_center));
				if (!open[i->ix]) {open[i->ix] = 1;}
				else if (new_g_score >= sn.g_score) continue; // not better
				sn.came_from_ix = cur;
				sn.path_pt.assign(pt.x, pt.y, cur_pt.z);
				if (i->ix == room2) {return reconstruct_path(state, avoid, cur_pt, radius, height, i->ix, room1, is_first_path, up_or_down, path);} // done, reconstruct path (in reverse)
				sn.g_score = new_g_score;
				sn.h_score = p2p_dist_xy(conn_center, dest_pos);
				sn.f_score = sn.g_score + sn.h_score;
				open_queue.push(make_pair(-sn.f_score, i->ix));
			} // for i
		} // end while()
		return 0; // failed - no path from room1 to room2
	}
}; // end building_nav_graph_t

void building_t::build_nav_graph() const {

	assert(interior);
	if (interior->nav_graph) return; // already built
	interior->nav_graph.reset(new building_nav_graph_t(0.5*get_window_vspace())); // set stairs_extend == doorway width
	building_nav_graph_t &ng(*interior->nav_graph);
	float const wall_width(0.5*get_floor_thickness());
	unsigned const num_rooms(interior->rooms.size()), num_stairs(interior->stairwells.size());
	ng.set_num_rooms(num_rooms, num_stairs);
	for (unsigned s = 0; s < num_stairs; ++s) {ng.set_stairs_bcube(s, interior->stairwells[s]);}

	for (unsigned r = 0; r < num_rooms; ++r) {
		room_t const &room(interior->rooms[r]);
		if (room.is_hallway) {ng.mark_hallway(r);}
		ng.set_room_bcube(r, room);
		cube_t c(room);
		c.expand_by_xy(wall_width); // to include adjacent doors

		for (auto d = doors.begin(); d != doors.end(); ++d) { // exterior doors
			if (!d->is_exterior_door() || d->type == tquad_with_ix_t::TYPE_RDOOR) continue;
			if (c.intersects(d->get_bcube())) {ng.mark_exit(r); break;}
		}
		for (auto d = interior->doors.begin(); d != interior->doors.end(); ++d) {
			if (!c.intersects_no_adj(*d)) continue; // door not adjacent to this room
			cube_t dc(*d);
			dc.expand_by_xy(wall_width); // to include adjacent rooms

			for (unsigned r2 = r+1; r2 < num_rooms; ++r2) { // check rooms with higher index (since graph is bidirectional)
				if (dc.intersects_no_adj(interior->rooms[r2])) {ng.connect_rooms(r, r2, *d); break;}
			}
		} // for d
		for (unsigned s = 0; s < num_stairs; ++s) { // stairs
			stairwell_t const &stairwell(interior->stairwells[s]);
			//if (!stairwell.stack_conn) continue;
			if (room.intersects_no_adj(stairwell)) {ng.connect_stairs(r, s, stairwell.dim, stairwell.dir);}
		}
		if (room.is_hallway) { // check for connected hallways
			for (unsigned r2 = r+1; r2 < num_rooms; ++r2) { // check rooms with higher index (since graph is bidirectional)
				room_t const &room2(interior->rooms[r2]);
				if (!room2.is_hallway || room2.z1() != room.z1() || !room2.intersects(c)) continue; // not a connected hallway
				cube_t conn_cube(c);
				conn_cube.intersect_with_cube(room2);
				ng.connect_rooms(r, r2, conn_cube);
			} // for r2
		}
		//for (unsigned e = 0; e < interior->elevators.size(); ++e) {} // elevators are not yet used by AIs so are ignored here
	} // for r
}

unsigned building_t::count_connected_room_components() const {
	if (!interior) return 0;
	build_nav_graph();
	unsigned const num(interior->nav_graph->count_connected_components());
	interior->nav_graph.reset(); // no longer needed
	return num;
}

point building_t::get_center_of_room(unsigned room_ix) const {
	assert(interior);
	assert(room_ix < interior->rooms.size());
	return interior->rooms[room_ix].get_cube_center();
}

bool building_t::choose_dest_room(building_ai_state_t &state, pedestrian_t &person, rand_gen_t &rgen, bool same_floor) const {

	assert(interior && interior->nav_graph);
	building_loc_t const loc(get_building_loc_for_pt(person.pos));
	if (loc.room_ix < 0) return 0; // not in a room
	state.cur_room   = loc.room_ix;
	state.dest_room  = state.cur_room; // set to the same room and pos just in case
	person.target_pos = person.pos; // set but not yet used
	if (interior->rooms.size() == 1) return 0; // no other room to move to

	for (unsigned n = 0; n < 100; ++n) { // make 100 attempts at finding a valid room
		unsigned const cand_room(rgen.rand() % interior->rooms.size());
		if (cand_room == loc.room_ix) continue;
		room_t const room(interior->rooms[cand_room]);
		if (room.is_hallway) continue; // don't select a hallway

		if (1 || same_floor) { // for now, always do this so that we don't have to handle walking between stacked parts
			if (person.pos.z < room.z1() || person.pos.z > room.z2()) continue; // room above or below the current pos
		}
		if (!interior->nav_graph->is_room_connected_to(loc.room_ix, cand_room)) continue;
		state.dest_room     = cand_room; // set but not yet used
		person.target_pos   = get_center_of_room(cand_room);
		person.target_pos.z = person.pos.z; // keep orig zval to stay on the same floor

		if (!same_floor) { // allow moving to a different floor, currently only one floor at a time
			assert(room.part_id < parts.size());
			cube_t const &part(parts[room.part_id]);
			unsigned const rand_val(rgen.rand() & 3); // 0-3

			if (rand_val == 0) { // try one floor below
				float const new_z(person.target_pos.z - get_window_vspace());
				if (new_z > part.z1()) {person.target_pos.z = new_z;} // change if there is a floor below
			}
			else if (rand_val == 1) { // try one floor above
				float const new_z(person.target_pos.z + get_window_vspace());
				if (new_z < part.z2()) {person.target_pos.z = new_z;} // change if there is a floor above
			}
		}
		return 1;
	} // for n
	return 0; // failed
}

template<typename T> void add_bcube_if_overlaps_zval(vector<T> const &cubes, vect_cube_t &out, float z1, float z2) {
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {
		if (i->z1() < z2 && i->z2() > z1) {out.push_back(*i);}
	}
}

void building_interior_t::get_avoid_cubes(vect_cube_t &avoid, float z1, float z2) const { // for AI
	add_bcube_if_overlaps_zval(stairwells, avoid, z1, z2);
	add_bcube_if_overlaps_zval(elevators,  avoid, z1, z2);
	if (!room_geom) return; // no room objects

	for (auto c = room_geom->objs.begin(); c != (room_geom->objs.begin() + room_geom->stairs_start); ++c) {
		if (c->no_coll() || c->type == TYPE_ELEVATOR || c->type == TYPE_STAIR || c->type == TYPE_LIGHT) continue; // the object types are not collided with
		if (c->z1() < z2 && c->z2() > z1) {avoid.push_back(*c);}
	}
}

int building_t::find_nearest_stairs(point const &p1, point const &p2, bool straight_only, int part_ix) const { // returns -1 on failure
	assert(interior);
	if (interior->stairwells.empty()) return -1; // no stairs
	assert(part_ix < 0 || (unsigned)part_ix < parts.size());
	float const zmin(min(p1.z, p2.z)), zmax(max(p1.z, p2.z));
	int ret(-1);
	float dmin(0.0);

	for (unsigned s = 0; s < interior->stairwells.size(); ++s) {
		stairwell_t const &stairs(interior->stairwells[s]);
		if (straight_only && stairs.shape == SHAPE_U) continue; // skip U-shaped stairs
		if (zmin < stairs.z1() || zmax > stairs.z2()) continue; // stairs don't span the correct floors
		if (part_ix >= 0 && !parts[part_ix].contains_cube(stairs)) continue; // stairs don't belong to this part
		point const center(stairs.get_cube_center());
		float const dist(p2p_dist(p1, center) + p2p_dist(center, p2));
		if (dmin == 0.0 || dist < dmin) {ret = s; dmin = dist;}
	} // for s
	return ret;
}

bool building_t::find_route_to_point(point const &from, point const &to, float radius, bool is_first_path, vector<point> &path) const {

	assert(interior && interior->nav_graph);
	path.clear();
	if (from == to) return 1; // already there???
	building_loc_t const loc1(get_building_loc_for_pt(from)), loc2(get_building_loc_for_pt(to));
	if (loc1.part_ix < 0 || loc2.part_ix < 0 || loc1.room_ix < 0 || loc2.room_ix < 0) return 0; // not in a room
	assert((unsigned)loc1.part_ix < parts.size() && (unsigned)loc2.part_ix < parts.size());
	assert((unsigned)loc1.room_ix < interior->rooms.size() && (unsigned)loc2.room_ix < interior->rooms.size());
	if (loc1 == loc2) {path.push_back(to); return 1;} // same part/room/floor, move to <to> and we're done
	bool use_stairs(0);
	if (loc1.floor != loc2.floor) {use_stairs = 1;} // different floors
	else if (loc1.part_ix != loc2.part_ix) {
		if (parts[loc1.part_ix].z1() != parts[loc2.part_ix].z1()) {use_stairs = 1;} // stacked parts
	}
	float const floor_spacing(get_window_vspace()), height(0.7*floor_spacing); // approximate, since we're not tracking actual heights
	static vect_cube_t avoid; // reuse across frames/people
	avoid.clear();
	interior->get_avoid_cubes(avoid, (from.z - height), (from.z + height)); // if using stairs, don't avoid stairs

	if (use_stairs) { // find path from <from> to nearest stairs, then find path from stairs to <to>
		int const stairs_ix(find_nearest_stairs(from, to, 1)); // straight_only=1; pass in loc1.part_ix if both loc part_ix values are equal?
		if (stairs_ix < 0) return 0; // no stairs can take us from <from> to <to>
		assert((unsigned)stairs_ix < interior->stairwells.size());
		stairwell_t const &stairs(interior->stairwells[stairs_ix]);
		bool const up_or_down(loc1.floor > loc2.floor); // 0=up, 1=down
		unsigned stairs_room_ix(stairs_ix + interior->rooms.size()); // map to graph space
		vector<point> from_path;
		// Note: passing use_stairs=0 here because it's unclear if we want to go through stairs nodes in our A* algorithm
		if (!interior->nav_graph->find_path_points(loc1.room_ix, stairs_room_ix, radius, height, 0, is_first_path, up_or_down, avoid, from, from_path)) return 0; // from => stairs
		point const seg2_start(interior->nav_graph->get_stairs_entrance_pt(to.z, (stairs_ix + interior->rooms.size()), !up_or_down)); // other end
		if (!interior->nav_graph->find_path_points(stairs_room_ix, loc2.room_ix, radius, height, 0, is_first_path, !up_or_down, avoid, seg2_start, path)) return 0; // stairs => to
		assert(!path.empty() && !from_path.empty());
		path.push_back(seg2_start); // other end of the stairs
		// add two more points to straighten the entrance and exit paths; this segment doesn't check for intersection with stairs
		cube_t stairs_ext(stairs);
		stairs_ext.d[stairs.dim][!stairs.dir] += (stairs.dir ? -1.0 : 1.0)*stairs.get_sz_dim(stairs.dim)/NUM_STAIRS_PER_FLOOR; // location of step up
		path.push_back(stairs_ext.closest_pt(path.back()));
		path.push_back(stairs_ext.closest_pt(from_path.front()));
		vector_add_to(from_path, path); // concatenate the two path segments in reverse order
	}
	else {
		if (!interior->nav_graph->find_path_points(loc1.room_ix, loc2.room_ix, radius, height, use_stairs, is_first_path, 0, avoid, from, path)) return 0; // failed to find a path
	}
	assert(!path.empty());
	return 1;
}

void building_ai_state_t::next_path_pt(pedestrian_t &person, bool same_floor) {
	assert(!path.empty());
	person.target_pos = path.back();
	if (same_floor) {person.target_pos.z = person.pos.z;} // make sure zvals are equal (ignore zval of room center)
	path.pop_back();
}

bool building_t::place_person(point &ppos, float radius, rand_gen_t &rgen) const {

	if (!interior || interior->rooms.empty()) return 0; // should be error case
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);

	for (unsigned n = 0; n < 100; ++n) { // make 100 attempts
		room_t const &room(interior->rooms[rgen.rand() % interior->rooms.size()]); // select a random room
		if (room.is_sec_bldg) continue; // don't place people in garages and sheds
		if (min(room.dx(), room.dy()) < 4.0*radius) continue; // room to small to place a person
		unsigned const num_floors(calc_num_floors(room, window_vspacing, floor_thickness));
		assert(num_floors > 0);
		unsigned const floor_ix(rgen.rand() % num_floors); // place person on a random floor
		// Note: people are placed before lights are assigned to rooms, so this may not work and must be handled during light placement
		if (room.lit_by_floor && !(room.lit_by_floor & (1ULL << (floor_ix&63)))) continue; // don't place person in an unlit room
		point pos;
		pos.z = room.z1() + fc_thick + window_vspacing*floor_ix;
		for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(room.d[d][0]+radius, room.d[d][1]-radius);} // random XY point inside this room
		cube_t bcube(pos);
		bcube.expand_by(radius); // expand more in Z?
		if (!is_valid_stairs_elevator_placement(bcube, radius, radius)) continue;
		bool bad_place(0);

		// Note: people are placed before room geom is generated for all buildings, so this may not work and will have to be handled during room geom placement
		if (interior->room_geom) { // check placement against room geom objects
			vector<room_object_t> const &objs(interior->room_geom->objs);
			auto objs_end(objs.begin() + interior->room_geom->stairs_start); // skip stairs and elevators

			for (auto i = objs.begin(); i != objs_end; ++i) {
				if (i->intersects(bcube)) {bad_place = 1; break;}
			}
		}
		if (bad_place) continue;
		ppos = pos;
		return 1;
	} // for n
	return 0;
}

int building_t::ai_room_update(building_ai_state_t &state, rand_gen_t &rgen, vector<pedestrian_t> &people, float delta_dir, unsigned person_ix, bool stay_on_one_floor) const {

	assert(person_ix < people.size());
	pedestrian_t &person(people[person_ix]);
	if (person.speed == 0.0) {person.anim_time = 0.0; return AI_STOP;} // stopped
	bool choose_dest(person.target_pos == all_zeros);
	float const radius_scale = 0.75; // somewhat smaller than radius
	float const coll_dist(radius_scale*person.radius);
	float &wait_time(person.waiting_start); // reuse this field

	if (wait_time > 0) {
		if (wait_time > fticks) { // waiting
			for (auto p = people.begin()+person_ix+1; p < people.end(); ++p) { // check for other people colliding with this person and handle it
				if (p->dest_bldg != person.dest_bldg) break; // done with this building
				if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
				float const rsum(coll_dist + radius_scale*p->radius);
				if (!dist_xy_less_than(person.pos, p->pos, rsum)) continue; // not intersecting
				move_person_to_not_collide(person, *p, person.pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
			} // for p
			wait_time -= fticks;
			return AI_WAITING;
		}
		wait_time = 0.0;
		choose_dest = 1;
	}
	assert(interior);
	assert(bcube.contains_pt(person.pos)); // person must be inside the building

	if (choose_dest) { // no current destination - choose a new one
		person.anim_time = 0.0; // reset animation
		if (!interior->nav_graph) {build_nav_graph();}
		// if there's no valid room or valid path, set the speed to 0 so that we don't check this every frame; movement will be stopped from now on
		if (!choose_dest_room(state, person, rgen, stay_on_one_floor)) {person.speed = 0.0; return AI_STOP;}

		if (!find_route_to_point(person.pos, person.target_pos, coll_dist, state.is_first_path, state.path)) {
			person.anim_time = 0.0;
			wait_time = 1.0*TICKS_PER_SECOND; // stop for 1 second then try again
			return AI_WAITING;
		}
		state.is_first_path = 0;
		state.next_path_pt(person, stay_on_one_floor);
		return AI_BEGIN_PATH;
	}
	float const max_dist(person.speed*fticks);

	if (dist_less_than(person.pos, person.target_pos, 1.1f*max_dist)) { // at dest
		assert(bcube.contains_pt(person.target_pos));
		person.pos = person.target_pos;
		if (!state.path.empty()) {state.next_path_pt(person, stay_on_one_floor); return AI_NEXT_PT;} // move to next path point
		person.anim_time = 0.0; // reset animation
		wait_time = TICKS_PER_SECOND*rgen.rand_uniform(1.0, 10.0); // stop for 1-10 seconds
		return AI_AT_DEST;
	}
	vector3d const new_dir(person.target_pos - person.pos);
	float const new_dir_mag(new_dir.mag());
	point new_pos;
	
	if (dot_product(new_dir, person.dir) < 0.999*new_dir_mag) { // dir not perfectly aligned
		//if (person.is_close_to_player()) {cout << TXT(new_dir.str()) << TXT(person.dir.str()) << TXT(new_dir_mag) << TXT(delta_dir) << TXT(max_dist) << TXT(person.radius) << endl;}
		assert(new_dir_mag > TOLERANCE); // should be guaranteed by dist_less_than() test, assuming zvals are equal (which they should be)
		float const step_scale(max(0.1f, dot_product(person.dir, new_dir)/new_dir_mag)); // move more slowly when direction misaligns to avoid overshooting target_pos
		person.dir = (delta_dir/new_dir_mag)*new_dir + (1.0 - delta_dir)*person.dir; // merge new_dir into dir gradually for smooth turning
		person.dir.normalize();
		new_pos = person.pos + (max_dist*step_scale)*person.dir;
		cube_t clip_cube(bcube);
		clip_cube.expand_by_xy(-coll_dist); // shrink
		clip_cube.clamp_pt_xy(new_pos); // make sure person stays within building bcube; can't clip to room because person may be exiting it
	}
	else { // optimization for aligned dir
		new_pos = person.pos + max_dist*person.dir;
	}
	// don't do collision detection while on stairs because it doesn't work properly; just let people walk through each other
	if (fabs(person.pos.z - person.target_pos.z) < 0.01*person.radius) {
		for (auto p = people.begin()+person_ix+1; p < people.end(); ++p) { // check all other people in the same building after this one and attempt to avoid them
			if (p->dest_bldg != person.dest_bldg) break; // done with this building
			if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
			float const rsum(coll_dist + radius_scale*p->radius);
			if (!dist_xy_less_than(new_pos, p->pos, rsum)) continue; // new pos not close
			if (!dist_xy_less_than(person.pos, p->pos, rsum)) return AI_STOP; // old pos not intersecting, stop
			person.anim_time = 0.0; // pause animation in case this person is mid-step
			move_person_to_not_collide(person, *p, new_pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
			return AI_MOVING; // return here, but don't update animation or dir; only handles a single collision
		} // for p
	}
	person.pos        = new_pos;
	person.anim_time += max_dist;
	return AI_MOVING;
}

void building_t::move_person_to_not_collide(pedestrian_t &person, pedestrian_t const &other, point const &new_pos, float rsum, float coll_dist) const {
	point other_pos(other.pos.x, other.pos.y, person.pos.z); // use same zval to ignore height differences
	float const sep_dist(p2p_dist_xy(person.pos, other_pos)), move_dist(rsum - sep_dist); // distance we have to move
	// move away from the other person, hopefully not through a wall
	point const orig_pos(person.pos);
	if (sep_dist > 0.01*coll_dist) {person.pos += (move_dist/sep_dist)*(person.pos - other_pos);}
	else {person.pos.x += rsum;} // avoid divide-by-zero, choose +X direction arbitrarily
	int const room_ix(get_building_loc_for_pt(orig_pos).room_ix);
	cube_t clip_bounds((room_ix >= 0 && (unsigned)room_ix < interior->rooms.size()) ? interior->rooms[room_ix] : bcube);
	clip_bounds.expand_by_xy(-coll_dist); // shrink
	clip_bounds.union_with_pt(orig_pos); // we know this point was valid
	clip_bounds.union_with_pt(new_pos);  // we know this point is valid
	clip_bounds.clamp_pt_xy(person.pos); // force player into the room

	if (!bcube.contains_pt(person.pos)) { // this can happen on rare occasions, due to fp inaccuracy or multiple collisions
		//cout << TXT(rsum) << TXT(sep_dist) << TXT(move_dist) << TXT(room_ix) << TXT(other.pos.str()) << TXT(person.pos.str()) << TXT(bcube.str()) << endl;
		bcube.clamp_pt_xy(person.pos); // just clamp pos so that it doesn't assert later
	}
}

void vect_building_t::ai_room_update(vector<building_ai_state_t> &ai_state, vector<pedestrian_t> &people, float delta_dir, rand_gen_t &rgen) const {
	//timer_t timer("Building People Update"); // ~3.7ms for 50K people, 0.55ms with distance check
	point const camera_bs(get_camera_pos() - get_tiled_terrain_model_xlate());
	float const dmax(1.5f*(X_SCENE_SIZE + Y_SCENE_SIZE));
	unsigned const num_people(people.size());
	ai_state.resize(num_people);

	for (unsigned i = 0; i < num_people; ++i) {
		if (!dist_less_than(people[i].pos, camera_bs, dmax)) continue; // too far away, no updates
		unsigned const bix(people[i].dest_bldg);
		assert(bix < size());
		operator[](bix).ai_room_update(ai_state[i], rgen, people, delta_dir, i, STAY_ON_ONE_FLOOR); // dispatch to the correct building
	}
}

building_loc_t building_t::get_building_loc_for_pt(point const &pt) const {
	building_loc_t loc;

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		if (p->contains_pt(pt)) {loc.part_ix = (p - parts.begin()); break;}
	}
	if (interior) { // rooms and stairwells, no elevators yet
		for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
			if (!r->contains_pt(pt)) continue;
			loc.room_ix = (r - interior->rooms.begin());
			if (r->is_sec_bldg) {loc.floor = 0;} // only one floor
			else {loc.floor = unsigned((pt.z - r->z1())/get_window_vspace());}
			break;
		}
		for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
			if (s->contains_pt(pt)) {loc.stairs_ix = (s - interior->stairwells.begin()); break;}
		}
	}
	return loc;
}

// these must be here to handle deletion of building_nav_graph_t, which is only defined in this file
building_interior_t::building_interior_t() {}
building_interior_t::~building_interior_t() {}
