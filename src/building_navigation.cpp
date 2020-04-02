// 3D World - Building Navigation for AIs
// by Frank Gennari 3/14/20

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for pedestrian_t
#include <queue>
#pragma warning(disable : 26812) // prefer enum class over enum


extern float fticks;

point get_cube_center_z1(cube_t const &c) {return point(c.get_center_dim(0), c.get_center_dim(1), c.z1());}

// Note: this should go into building_t/buildings.h at some point, but is temporarily here
class building_nav_graph_t {
	struct pt_with_ix_t {
		unsigned ix;
		point pt;
		pt_with_ix_t(unsigned ix_, point const &pt_) : ix(ix_), pt(pt_) {}
	};
	struct node_t { // represents one room or one stairwell
		bool has_exit, is_hallway, is_stairs;
		cube_t bcube;
		vector<pt_with_ix_t> conn_rooms;
		node_t() : has_exit(0), is_hallway(0), is_stairs(0) {}
		point get_center() const {return get_cube_center_z1(bcube);}

		void add_conn_room(unsigned room, cube_t const &c) {
			for (auto i = conn_rooms.begin(); i != conn_rooms.end(); ++i) {if (i->ix == room) return;} // ignore duplicates
			conn_rooms.emplace_back(room, get_cube_center_z1(c));
		}
		point get_conn_pt(unsigned room) const {
			for (auto i = conn_rooms.begin(); i != conn_rooms.end(); ++i) {if (i->ix == room) return i->pt;}
			assert(0);
			return all_zeros; // never gets here
		}
	};
	struct a_star_node_state_t {
		int came_from_ix;
		point path_pt;
		float g_score, h_score, f_score;
		a_star_node_state_t() : came_from_ix(-1), g_score(0), h_score(0), f_score(0) {}
	};

	unsigned num_rooms, num_stairs;
	vector<node_t> nodes;
	node_t       &get_node(unsigned room)       {assert(room < nodes.size()); return nodes[room];}
	node_t const &get_node(unsigned room) const {assert(room < nodes.size()); return nodes[room];}

	void add_room_to_path(unsigned room1, unsigned room2, vector<point> &path) const {
		path.push_back(get_node(room1).get_conn_pt(room2)); // door pos
		path.push_back(get_node(room2).get_center()); // next room center
	}
	void remove_connection(unsigned from, unsigned to) {
		auto &conn(get_node(from).conn_rooms);

		for (auto i = conn.begin(); i != conn.end(); ++i) {
			if (i->ix == to) {swap(*i, conn.back()); conn.pop_back(); return;}
		}
		assert(0); // must be found - should not get here
	}
public:
	building_nav_graph_t() : num_rooms(0), num_stairs(0) {}

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
		cube_t entry(n2.bcube);
		entry.d[dim][dir] = entry.d[dim][!dir]; // shrink to zero area at the entrance to the stairs when going up (need to reverse when going down)
		get_node(room).add_conn_room(node_ix2, entry);
		n2.add_conn_room(room, entry);
	}
	void connect_rooms(unsigned room1, unsigned room2, cube_t const &conn_bcube) { // graph is bidirectional
		assert(room1 < num_rooms && room2 < num_rooms);
		get_node(room1).add_conn_room(room2, conn_bcube);
		get_node(room2).add_conn_room(room1, conn_bcube);
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
	point find_valid_room_dest(vect_cube_t const &avoid, float radius, float height, unsigned node_ix, bool &not_room_center) const {
		node_t const &node(get_node(node_ix));
		cube_t place_area(node.bcube);
		place_area.expand_by_xy(-2.0*radius); // shrink by twice the radius
		point const center(get_cube_center_z1(node.bcube));
		if (!place_area.is_strictly_normalized()) return center; // should generally not be true
		not_room_center = 1;
		static unsigned call_ix(1); // unique value for every call
		rand_gen_t rgen;
		rgen.set_state(node_ix, call_ix++);
		point pos(center); // first candidate is the center of the room

		for (unsigned n = 0; n < 100; ++n) { // 100 random tries to find a valid dest_pos
			if (is_valid_pos(avoid, pos, radius, height)) return pos; // success
			for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(place_area.d[d][0], place_area.d[d][1]);} // choose a random new point in the room
		}
		return center; // failed, return room center
	}

	bool connect_room_endpoints(vect_cube_t const &avoid, cube_t const &walk_area, point const &p1, point const &p2, float radius, vector<point> &path) const {
		bool is_path_valid(1);

		for (auto i = avoid.begin(); i != avoid.end(); ++i) {
			if (!i->intersects_xy(walk_area)) continue;
			cube_t c(*i);
			c.expand_by_xy(radius);
			if (check_line_clip_xy(p1, p2, c.d)) {is_path_valid = 0; break;}
		}
		if (is_path_valid) return 1; // done
		// straight line not valid, need to run path finding
		vect_cube_t keepout, nav_mesh, temp;

		for (auto i = avoid.begin(); i != avoid.end(); ++i) {
			if (i->intersects_xy(walk_area)) {keepout.push_back(*i); keepout.back().expand_by_xy(radius);}
		}
		assert(!keepout.empty());
		subtract_cubes_from_cube(walk_area, keepout, nav_mesh, temp);
		// TODO: find route from p1 to p2
		return 1;
	}
	bool add_path_for_room(vect_cube_t const &avoid, cube_t const &room, point const &next, float radius, vector<point> &path) const {
		cube_t walk_area(room);
		walk_area.expand_by_xy(-radius); // shrink by radius
		point const &prev(path.back());
		point const p1(walk_area.closest_pt(prev)), p2(walk_area.closest_pt(next));
		if (p1 != prev) {path.push_back(p1);} // walk out of doorway and into room
		if (!connect_room_endpoints(avoid, walk_area, p1, p2, radius, path)) return 0;
		if (p2 != next) {path.push_back(p2);} // walk from room into doorway
		return 1;
	}
	bool reconstruct_path(vector<a_star_node_state_t> const &state, vect_cube_t const &avoid,
		float radius, float height, unsigned start_ix, unsigned end_ix, vector<point> &path) const
	{
		unsigned n(start_ix);
		rand_gen_t rgen;
		rgen.set_state(start_ix, nodes.size());

		while (1) {
			assert(n < nodes.size());
			node_t const &node(get_node(n));
			point const &next(state[n].path_pt);
			int const came_from(state[n].came_from_ix);
			bool const is_first_pt(path.empty());

			if (is_first_pt) { // last point in path (first point in reverse path)
				assert(came_from >= 0); // too strong?
				cube_t walk_area(node.bcube);
				walk_area.expand_by_xy(-radius); // shrink by radius
				bool success(0);

				for (unsigned n = 0; n < 10; ++n) { // keep retrying until we find a point that is reachable from the doorway
					bool not_room_center(0);
					point const end_point(find_valid_room_dest(avoid, radius, height, start_ix, not_room_center));
					path.push_back(end_point);
					if (connect_room_endpoints(avoid, walk_area, end_point, next, radius, path)) {success = 1; break;}
					path.clear(); // failed, reset for next iteration
					if (!not_room_center) break; // if we did choose the room center, and there is no path to it, we've failed
				} // for n
				if (!success) {assert(path.empty()); return 0;} // failed to connect to a point in dest room
			}
			else { // use room center for intermediate points; person is unlikely to go through these points anyway
				path.push_back(node.get_center());
			}
			if (came_from < 0) {assert(n == end_ix); break;} // done

			if (!is_first_pt && min(node.bcube.dx(), node.bcube.dy()) > 3.0*radius) { // adjust the path through a room
				assert(!path.empty());
				path.pop_back(); // remove room center point

				if (!add_path_for_room(avoid, node.bcube, next, radius, path)) {
					path.clear();
					// TODO: try another path?
					//disconnect_room_pair(n, came_from); // ???
					return 0;
				}
			}
			path.push_back(next); // doorway
			n = came_from;
		} // end while()
		return 1;
	}
	
	// A* algorithm; Note: path is stored backwards
	bool find_path_points(unsigned room1, unsigned room2, float radius, float height, bool use_stairs, vect_cube_t const &avoid, vector<point> &path) const {
		assert(room1 < nodes.size() && room2 < nodes.size());
		assert(room1 != room2); // or just return an empty path?
		path.clear();
		// TODO: if use_stairs==1, or room1 and room2 are in different stacks, the path must include stairs
		vector<a_star_node_state_t> state(nodes.size());
		vector<uint8_t> open(nodes.size(), 0), closed(nodes.size(), 0); // tentative/already evaluated nodes
		std::priority_queue<pair<float, unsigned> > open_queue;
		point const dest_pos(get_node(room2).get_center());
		a_star_node_state_t &start(state[room1]);
		start.g_score = 0.0;
		start.h_score = start.f_score = p2p_dist(get_node(room1).get_center(), dest_pos); // estimated total cost from start to goal through current
		open[room1]   = 1;
		open_queue.push(make_pair(-start.f_score, room1));

		while (!open_queue.empty()) {
			unsigned const cur(open_queue.top().second);
			open_queue.pop();
			if (closed[cur]) continue; // already closed (duplicate)
			if (cur == room2) {return reconstruct_path(state, avoid, radius, height, cur, room1, path);} // done, reconstruct path (in reverse)
			a_star_node_state_t const &cs(state[cur]);
			node_t const &cur_node(get_node(cur));
			point const center(cur_node.get_center());
			assert(!closed[cur]);
			closed[cur] = 1;
			open[cur]   = 0;

			for (auto i = cur_node.conn_rooms.begin(); i != cur_node.conn_rooms.end(); ++i) {
				assert(i->ix < nodes.size());
				if (closed[i->ix]) continue; // already closed (duplicate)
				node_t const &conn_node(get_node(i->ix));
				if (conn_node.is_stairs && !use_stairs) continue; // skip stairs in this mode
				point const conn_center(conn_node.get_center());
				a_star_node_state_t &sn(state[i->ix]);
				float const new_g_score(sn.g_score + p2p_dist(center, i->pt) + p2p_dist(i->pt, conn_center));
				if (!open[i->ix]) {open[i->ix] = 1;}
				else if (new_g_score >= sn.g_score) continue; // not better
				sn.came_from_ix = cur;
				sn.path_pt = i->pt;
				sn.g_score = new_g_score;
				sn.h_score = p2p_dist(conn_center, dest_pos);
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
	interior->nav_graph.reset(new building_nav_graph_t);
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

		if (same_floor) {
			cube_t const room(interior->rooms[cand_room]);
			if (person.pos.z < room.z1() || person.pos.z > room.z2()) continue; // room above or below the current pos
		}
		if (!interior->nav_graph->is_room_connected_to(loc.room_ix, cand_room)) continue;

		if (!same_floor) {
			// do some stuff with stairs
		}
		state.dest_room     = cand_room; // set but not yet used
		person.target_pos   = get_center_of_room(cand_room);
		person.target_pos.z = person.pos.z; // keep orig zval to stay on the same floor
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

bool building_t::find_route_to_point(point const &from, point const &to, float radius, vector<point> &path) const {

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
	float const height(0.7*get_window_vspace()); // approximate, since we're not tracking actual heights
	static vect_cube_t avoid; // reuse across frames/people
	avoid.clear();
	interior->get_avoid_cubes(avoid, (from.z - height), (from.z + height));
	if (!interior->nav_graph->find_path_points(loc1.room_ix, loc2.room_ix, radius, height, use_stairs, avoid, path)) return 0; // failed to find a path
	assert(!path.empty());
	if (path.size() > 1) {path.pop_back();} // don't need to go through the center of the first room
	return 1;
}

void building_ai_state_t::next_path_pt(pedestrian_t &person, bool same_floor) {
	assert(!path.empty());
	person.target_pos = path.back();
	if (same_floor) {person.target_pos.z = person.pos.z;} // make sure zvals are equal (ignore zval of room center)
	path.pop_back();
}

int building_t::ai_room_update(building_ai_state_t &state, rand_gen_t &rgen, vector<pedestrian_t> &people, float delta_dir, unsigned person_ix, bool stay_on_one_floor) const {

	assert(person_ix < people.size());
	pedestrian_t &person(people[person_ix]);
	if (person.speed == 0.0) {person.anim_time = 0.0; return AI_STOP;} // stopped
	bool choose_dest(person.target_pos == all_zeros);
	float const coll_dist(0.7*person.radius); // somewhat smaller than radius

	if (state.wait_time > 0) {
		if (state.wait_time > fticks) { // waiting
			for (auto p = people.begin()+person_ix+1; p < people.end(); ++p) { // check for other people colliding with this person and handle it
				if (p->dest_bldg != person.dest_bldg) break; // done with this building
				if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
				float const rsum(coll_dist + 0.7*p->radius);
				if (!dist_xy_less_than(person.pos, p->pos, rsum)) continue; // not intersecting
				move_person_to_not_collide(person, *p, person.pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
			} // for p
			state.wait_time -= fticks;
			return AI_WAITING;
		}
		state.wait_time = 0.0;
		choose_dest = 1;
	}
	assert(interior);
	assert(bcube.contains_pt(person.pos)); // person must be inside the building

	if (choose_dest) { // no current destination - choose a new one
		person.anim_time = 0.0; // reset animation
		if (!interior->nav_graph) {build_nav_graph();}
		// if there's no valid room or valid path, set the speed to 0 so that we don't check this every frame; movement will be stopped from now on
		if (!choose_dest_room(state, person, rgen, stay_on_one_floor)) {person.speed = 0.0; return AI_STOP;}
		if (!find_route_to_point(person.pos, person.target_pos, coll_dist, state.path)) {cout << "*** Bad Path ***"; person.speed = 0.0; return AI_STOP;}
		state.next_path_pt(person, stay_on_one_floor);
		return AI_AT_DEST;
	}
	float const max_dist(person.speed*fticks);

	if (dist_less_than(person.pos, person.target_pos, max_dist)) { // at dest
		if (!stay_on_one_floor && person.pos.z != person.target_pos.z) {} // handle this case
		assert(bcube.contains_pt(person.target_pos));
		person.anim_time = 0.0; // reset animation
		person.pos = person.target_pos;
		
		if (!state.path.empty()) { // move to next path point
			state.next_path_pt(person, stay_on_one_floor);
			return AI_NEXT_PT;
		}
		state.wait_time = TICKS_PER_SECOND*rgen.rand_uniform(1.0, 10.0); // stop for 1-10 seconds
		return AI_AT_DEST;
	}
	vector3d new_dir(person.target_pos - person.pos);
	new_dir.z = 0.0; // XY only, even if on stairs
	float const new_dir_mag(new_dir.mag());
	point new_pos;
	
	if (dot_product(new_dir, person.dir) < 0.999*new_dir_mag) { // dir not perfectly aligned
		assert(new_dir_mag > TOLERANCE); // should be guaranteed by dist_less_than() test, assuming zvals are equal (which they should be)
		person.dir = (delta_dir/new_dir_mag)*new_dir + (1.0 - delta_dir)*person.dir; // merge new_dir into dir gradually for smooth turning
		person.dir.normalize();
		new_pos = person.pos + max_dist*person.dir;
		cube_t clip_cube(bcube);
		clip_cube.expand_by_xy(-coll_dist); // shrink
		clip_cube.clamp_pt_xy(new_pos); // make sure person stays within building bcube; can't clip to room because person may be exiting it
	}
	else { // optimization for aligned dir
		new_pos = person.pos + max_dist*person.dir;
	}
	for (auto p = people.begin()+person_ix+1; p < people.end(); ++p) { // check all other people in the same building after this one and attempt to avoid them
		if (p->dest_bldg != person.dest_bldg) break; // done with this building
		if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
		float const rsum(coll_dist + 0.7*p->radius);
		if (!dist_xy_less_than(new_pos, p->pos, rsum)) continue; // new pos not close
		if (!dist_xy_less_than(person.pos, p->pos, rsum)) return AI_STOP; // old pos not intersecting, stop
		person.anim_time = 0.0; // pause animation in case this person is mid-step
		move_person_to_not_collide(person, *p, new_pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
		return AI_MOVING; // return here, but don't update animation or dir; only handles a single collision
	} // for p
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
	//timer_t timer("Building People Update"); // ~2.7ms
	bool const stay_on_one_floor = 1; // multi-floor movement not yet supported
	unsigned const num_people(people.size());
	ai_state.resize(num_people);

	for (unsigned i = 0; i < num_people; ++i) {
		unsigned const bix(people[i].dest_bldg);
		assert(bix < size());
		operator[](bix).ai_room_update(ai_state[i], rgen, people, delta_dir, i, stay_on_one_floor); // dispatch to the correct building
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
