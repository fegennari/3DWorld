// 3D World - Building Navigation for AIs
// by Frank Gennari 3/14/20

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for pedestrian_t
#include <queue>


bool const STAY_ON_ONE_FLOOR  = 0;
float const COLL_RADIUS_SCALE = 0.75; // somewhat smaller than radius, but larger than PED_WIDTH_SCALE

int cpbl_update_frame(0);
building_dest_t cur_player_building_loc, prev_player_building_loc;

extern int frame_counter, display_mode, player_in_closet, animate2;
extern float fticks;
extern double camera_zh;
extern building_params_t global_building_params;
extern bldg_obj_type_t bldg_obj_types[];

bool in_building_gameplay_mode();
bool ai_follow_player() {return (global_building_params.ai_follow_player || in_building_gameplay_mode());}
bool can_ai_follow_player(pedestrian_t const &person);
bool get_closest_building_sound(point const &at_pos, point &sound_pos, float floor_spacing);
void maybe_play_zombie_sound(point const &sound_pos_bs, unsigned zombie_ix, bool alert_other_zombies);
bool register_ai_player_coll(bool &has_key);

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
	bool invalid;
	building_nav_graph_t(float stairs_extend_) : num_rooms(0), num_stairs(0), stairs_extend(stairs_extend_), invalid(0) {}

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

	void connect_stairs(unsigned room, unsigned stairs, bool dim, bool dir, bool is_u) {
		assert(room < num_rooms && stairs < num_stairs);
		unsigned const node_ix2(num_rooms + stairs);
		node_t &n2(get_node(node_ix2));
		cube_t entry_u(n2.bcube), entry_d(n2.bcube);
		float const extend((dir ? -1.0 : 1.0)*stairs_extend); // extend away from stairs for entrance/exit area; will be denormalized in this dim

		if (is_u) { // U-shaped stairs: entrances are on the same side
			bool const side(dir); // Note: see code in add_stairs_and_elevators()
			entry_u.d[dim][dir] = entry_d.d[dim][dir] = entry_u.d[dim][!dir] + extend; // shrink to extend length at the entrance to the stairs
			float const mid(n2.bcube.get_center_dim(!dim));
			entry_u.d[!dim][ side] = mid; // bottom
			entry_d.d[!dim][!side] = mid; // top
		}
		else { // straight stairs: entrances are on opposite ends
			entry_u.d[dim][ dir] = entry_u.d[dim][!dir] + extend; // shrink to extend length at the entrance to the stairs when going up
			entry_d.d[dim][!dir] = entry_d.d[dim][ dir] - extend; // shrink to extend length at the entrance to the stairs when going down
		}
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
		cube_t walk_area_exp(walk_area);
		walk_area_exp.expand_by_xy(radius); // to capture objects that we could intersect when our center is in walk_area
		keepout.clear();

		for (auto i = avoid.begin(); i != avoid.end(); ++i) {
			if (!i->intersects_xy(walk_area_exp)) continue;
			cube_t c(*i);
			c.expand_by_xy(radius);

			if (c.contains_pt_xy(p1)) {
				if (ignore_p1_coll) {
					if (i->contains_pt_xy(p1)) continue; // even the unexpanded cube intersects, skip it
					c = *i; // use unexpanded cube instead
				}
				else {return 0;} // fail
			}
			if (c.contains_pt_xy(p2)) {
				if (ignore_p2_coll) {
					if (i->contains_pt_xy(p2)) continue; // even the unexpanded cube intersects, skip it
					c = *i; // use unexpanded cube instead
				}
				else {return 0;} // fail
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
		return point(max(c.x1(), min(c.x2(), pos.x)), max(c.y1(), min(c.y2(), pos.y)), pos.z); // c.clamp_pt_xy(), but returning a new point
	}
	static cube_t calc_walkable_room_area(node_t const &node, float radius) {
		cube_t walk_area(node.bcube);
		if (!node.is_stairs) {walk_area.expand_by_xy(-radius);} // shrink by radius

		if (node.is_hallway) {
			bool const min_dim(walk_area.dy() < walk_area.dx());
			float const shrink(min(0.5f*radius, max(0.0f, 0.9f*0.5f*walk_area.get_sz_dim(min_dim)))); // shrink by an extra half radius if hallway is wide enough
			walk_area.expand_in_dim(min_dim, -shrink);
		}
		return walk_area;
	}
	bool reconstruct_path(vector<a_star_node_state_t> const &state, vect_cube_t const &avoid, point const &cur_pt, float radius,
		float height, unsigned start_ix, unsigned end_ix, unsigned ped_ix, bool is_first_path, bool up_or_down, unsigned ped_rseed, vector<point> &path) const
	{
		unsigned n(start_ix);
		rand_gen_t rgen;
		rgen.set_state((ped_ix + 17*start_ix + 1), (ped_rseed + 1));
		vect_cube_t keepout;

		while (1) {
			node_t const &node(get_node(n));
			point const &next(state[n].path_pt);
			int const came_from(state[n].came_from_ix);
			bool const is_first_pt(path.empty());
			cube_t const walk_area(calc_walkable_room_area(node, radius));
			
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
					// try another path? this case is rare; on failure, the person will wait a second then choose a different destination room
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
	bool complete_path_within_room(point const &p1, point const &p2, unsigned room, unsigned ped_ix, float radius, unsigned ped_rseed, vect_cube_t const &avoid, vector<point> &path) const {
		// used for reaching a goal such as the player within the same room
		cube_t const walk_area(calc_walkable_room_area(get_node(room), radius));
		rand_gen_t rgen;
		rgen.set_state((ped_ix + 13*room + 1), (ped_rseed + 1));
		vect_cube_t keepout;
		path.push_back(p2); // Note: path is constructed backwards, so p2 is added first and connect_room_endpoints takes swapped arguments
		// ignore starting collisions, for example collisions with stairwell when exiting stairs?
		if (!connect_room_endpoints(avoid, walk_area, p2, p1, radius, path, keepout, rgen, 0, 1)) {path.clear(); return 0;} // ignore initial coll with p1
		// maybe add an extra path point to prevent clipping through walls when walking throug a doorway
		point const p1_extend_pt(closest_room_pt(walk_area, p1));
		if (p1_extend_pt != p1) {path.push_back(p1_extend_pt);} // add if p1 was clamped
		return 1;
	}
	
	// A* algorithm; Note: path is stored backwards
	bool find_path_points(unsigned room1, unsigned room2, unsigned ped_ix, float radius, float height, bool use_stairs, bool is_first_path,
		bool up_or_down, unsigned ped_rseed, vect_cube_t const &avoid, point const &cur_pt, vector<point> &path) const
	{
		// Note: opening and closing doors updates the nav graph; an AI encountering a closed door after choosing a path can either open it or stop and wait
		assert(room1 < nodes.size() && room2 < nodes.size());
		assert(room1 != room2);
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
				
				if (i->ix == room2) { // done, reconstruct path (in reverse)
					return reconstruct_path(state, avoid, cur_pt, radius, height, i->ix, room1, ped_ix, is_first_path, up_or_down, ped_rseed, path);
				}
				sn.g_score = new_g_score;
				sn.h_score = p2p_dist_xy(conn_center, dest_pos);
				sn.f_score = sn.g_score + sn.h_score;
				open_queue.push(make_pair(-sn.f_score, i->ix));
			} // for i
		} // end while()
		return 0; // failed - no path from room1 to room2
	}
}; // end building_nav_graph_t

cube_t building_t::get_walkable_room_bounds(room_t const &room) const {
	cube_t c(room);
	// Note: regular house rooms start and end at the walls; offices and hallways tile exactly and include half the walls, so we have to subtract those back off
	if (room.is_hallway || room.is_office) {c.expand_by_xy(-0.5*get_wall_thickness());}
	return c;
}

void building_t::build_nav_graph() const {

	assert(interior);
	if (interior->nav_graph && !interior->nav_graph->invalid) return; // already built
	interior->nav_graph.reset(new building_nav_graph_t(0.5*get_window_vspace())); // set stairs_extend == doorway width
	building_nav_graph_t &ng(*interior->nav_graph);
	float const wall_width(get_wall_thickness());
	unsigned const num_rooms(interior->rooms.size()), num_stairs(interior->stairwells.size());
	ng.set_num_rooms(num_rooms, num_stairs);
	for (unsigned s = 0; s < num_stairs; ++s) {ng.set_stairs_bcube(s, interior->stairwells[s]);}

	for (unsigned r = 0; r < num_rooms; ++r) {
		room_t const &room(interior->rooms[r]);
		cube_t c(get_walkable_room_bounds(room));
		if (room.is_hallway) {ng.mark_hallway(r);}
		ng.set_room_bcube(r, c);
		c.expand_by_xy(wall_width); // to include adjacent doors
		if (is_room_adjacent_to_ext_door(c)) {ng.mark_exit(r);}

		for (auto d = interior->doors.begin(); d != interior->doors.end(); ++d) {
			// door starts off closed/locked, treat it as a barrier for now and don't connect the rooms; we could check person.has_key, but this graph is shared across all people
			if (d->is_closed_and_locked() || (!d->open && global_building_params.ai_opens_doors < 2)) continue;
			if (!c.intersects_no_adj(*d)) continue; // door not adjacent to this room
			cube_t dc(*d);
			dc.expand_by_xy(wall_width); // to include adjacent rooms

			for (unsigned r2 = r+1; r2 < num_rooms; ++r2) { // check rooms with higher index (since graph is bidirectional)
				if (dc.intersects_no_adj(interior->rooms[r2])) {ng.connect_rooms(r, r2, *d); break;}
			}
		} // for d
		for (unsigned s = 0; s < num_stairs; ++s) { // stairs
			stairwell_t const &stairwell(interior->stairwells[s]);

			if (stairwell.stairs_door_ix >= 0 && global_building_params.ai_opens_doors < 2) { // check for open doors; doors on stairs can't be locked
				assert((unsigned)stairwell.stairs_door_ix < interior->doors.size());
				if (!interior->doors[stairwell.stairs_door_ix].open) continue; // stairs blocked by closed door, don't connect
			}
			if (room.intersects_no_adj(stairwell)) {ng.connect_stairs(r, s, stairwell.dim, stairwell.dir, (stairwell.shape == SHAPE_U));}
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
	//if (is_house && has_basement()) {cout << TXT(ng.is_fully_connected()) << endl;}
}

void building_t::invalidate_nav_graph() { // Note: this is safe to call in one thread while using in another
	if (interior && interior->nav_graph) {interior->nav_graph->invalid = 1;}
}

unsigned building_t::count_connected_room_components() {
	if (!interior) return 0;
	build_nav_graph();
	unsigned const num(interior->nav_graph->count_connected_components());
	interior->nav_graph.reset(); // no longer needed
	return num;
}

bool building_t::is_room_adjacent_to_ext_door(cube_t const &room, bool front_door_only) const {
	cube_t room_exp(room);
	room_exp.expand_by_xy(get_wall_thickness());

	for (auto d = doors.begin(); d != doors.end(); ++d) { // exterior doors
		if (!d->is_exterior_door() || d->type == tquad_with_ix_t::TYPE_RDOOR) continue;
		if (room_exp.intersects(d->get_bcube())) return 1;
		if (front_door_only) return 0; // assumes the first door is the front door
	}
	return 0;
}

// Warning: this may be called from a different thread from the one that uses it for AI updates
void building_t::register_player_in_building(point const &camera_bs, unsigned building_id) const {
	if (animate2) {prev_player_building_loc = cur_player_building_loc;} // only update previous pos when AI is running so that it doesn's miss a floor or room change
	cur_player_building_loc = building_dest_t(get_building_loc_for_pt(camera_bs), camera_bs, building_id);
	cpbl_update_frame       = frame_counter;
}
void end_register_player_in_building() {
	if (cpbl_update_frame != frame_counter) {prev_player_building_loc = cur_player_building_loc = building_dest_t();} // player not in building, reset
}

bool building_t::choose_dest_goal(building_ai_state_t &state, pedestrian_t &person, rand_gen_t &rgen, bool same_floor) const {

	assert(interior && interior->nav_graph);
	building_loc_t const loc(get_building_loc_for_pt(person.pos));
	if (loc.room_ix < 0) return 0; // not in a room
	float const floor_spacing(get_window_vspace());
	building_dest_t goal;
	point sound_pos;

	if ((global_building_params.ai_target_player || cur_player_building_loc.same_room_floor(loc)) && can_target_player(state, person)) {
		goal = cur_player_building_loc; // player is in a different room of our building, or we're following the player's position
		state.goal_type = GOAL_TYPE_PLAYER;
	}
	else if (can_ai_follow_player(person) && get_closest_building_sound(person.pos, sound_pos, floor_spacing)) { // target the loudest sound
		goal = building_dest_t(get_building_loc_for_pt(sound_pos), sound_pos, cur_player_building_loc.building_ix); // same building as player (current building)
		state.goal_type = GOAL_TYPE_SOUND;
	}
	if (!goal.is_valid()) return 0; // player or sound
	unsigned const cand_room(goal.room_ix);
	room_t const &room(get_room(cand_room)); // target room
	if (!interior->nav_graph->is_room_connected_to(loc.room_ix, cand_room)) return 0; // unreachable
	state.cur_room      = loc.room_ix;
	state.dest_room     = cand_room; // set but not yet used
	person.target_pos   = (global_building_params.ai_target_player ? goal.pos: get_center_of_room(cand_room));
	person.target_pos.z = person.pos.z; // keep orig zval to stay on the same floor
	person.is_on_stairs = 0;

	if (!same_floor) { // allow moving to a different floor, currently only one floor at a time
		//room_t const &room(get_room(loc.room_ix));
		if (goal.floor_ix < loc.floor_ix) { // try one floor below
			float const new_z(person.target_pos.z - floor_spacing);
			if (new_z > room.z1()) {person.target_pos.z = new_z;} // change if there is a floor below
		}
		else if (goal.floor_ix > loc.floor_ix) { // try one floor above
			float const new_z(person.target_pos.z + floor_spacing);
			if (new_z < room.z2()) {person.target_pos.z = new_z;} // change if there is a floor above
		}
		// else if floors differ by more than 1, we'll end up visiting the room on the wrong floor
	}
	if (global_building_params.ai_target_player) { // ensure target is a valid location in this building; this must be done *after* adjacent floor zval adjustment
		float const z2_add(person.get_height() - person.radius), coll_dist(COLL_RADIUS_SCALE*person.radius);
		cube_t legal_area(bcube);
		legal_area.expand_by_xy(-coll_dist);
		legal_area.z1() += person.radius;
		legal_area.z2() -= z2_add;
		assert(legal_area.is_strictly_normalized());
		legal_area.clamp_pt(person.target_pos); // clamp to building interior
		float dmin_sq(bcube.get_max_extent()); // start at a large value
		cube_t closest_part;

		for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) { // there shouldn't be any people in secondary buildings, but include them anyway
			if (p->contains_pt(person.target_pos)) {dmin_sq = 0; break;} // done
			float const dsq(p2p_dist_sq(person.target_pos, p->closest_pt(person.target_pos)));
			if (dsq < dmin_sq) {closest_part = *p;}
		}
		if (dmin_sq > 0.0 && !closest_part.is_all_zeros()) {closest_part.clamp_pt(person.target_pos);} // clamp to closest part
		static vect_cube_t avoid; // reuse across frames/people
		interior->get_avoid_cubes(avoid, (person.target_pos.z - person.radius), (person.target_pos.z + z2_add), get_floor_thickness(), 1); // same_as_player=1

		for (unsigned n = 0; n < 4; ++n) { // iterate a few times in case a collision moves pos into another object
			for (auto i = avoid.begin(); i != avoid.end(); ++i) { // move target_pos to avoid room objects
				cube_t c(*i);
				c.expand_by_xy(coll_dist);
				sphere_cube_int_update_pos(person.target_pos, 1.01*coll_dist, c, person.pos, 1, 1); // check_int=1, skip_z=1, ignore return value
			}
		}
	}
	return 1;
}

int building_t::choose_dest_room(building_ai_state_t &state, pedestrian_t &person, rand_gen_t &rgen, bool same_floor) const {

	assert(interior && interior->nav_graph);
	building_loc_t const loc(get_building_loc_for_pt(person.pos));
	if (loc.room_ix < 0) return 0; // not in a room
	state.dest_room     = state.cur_room = loc.room_ix; // set to the same room and pos just in case
	person.target_pos   = person.pos; // set but not yet used
	person.is_on_stairs = 0;
	if (interior->rooms.size() == 1) return 0; // no other room to move to
	float const floor_spacing(get_window_vspace());

	for (unsigned n = 0; n < 100; ++n) { // make 100 attempts at finding a valid room
		unsigned const cand_room(rgen.rand() % interior->rooms.size());
		if (cand_room == (unsigned)loc.room_ix) continue;
		room_t const &room(interior->rooms[cand_room]);
		if (room.is_hallway) continue; // don't select a hallway

		if (1 || same_floor) { // for now, always do this so that we don't have to handle walking between stacked parts
			if (person.pos.z < room.z1() || person.pos.z > room.z2()) continue; // room above or below the current pos (won't walk between stacks)
		}
		if (!interior->nav_graph->is_room_connected_to(loc.room_ix, cand_room)) continue;
		state.dest_room     = cand_room; // set but not yet used
		person.target_pos   = get_center_of_room(cand_room);
		person.target_pos.z = person.pos.z; // keep orig zval to stay on the same floor

		if (!same_floor) { // allow moving to a different floor, currently only one floor at a time
			cube_t const &part(get_part_for_room(room));
			unsigned const rand_val(rgen.rand() & 3); // 0-3

			if (rand_val == 0) { // try one floor below
				float const new_z(person.target_pos.z - floor_spacing);
				if (new_z > part.z1()) {person.target_pos.z = new_z;} // change if there is a floor below
			}
			else if (rand_val == 1) { // try one floor above
				float const new_z(person.target_pos.z + floor_spacing);
				if (new_z < part.z2()) {person.target_pos.z = new_z;} // change if there is a floor above
			}
		}
		state.goal_type = GOAL_TYPE_ROOM;
		return 1;
	} // for n
	return 2; // failed, but can retry
}

template<typename T> void add_bcube_if_overlaps_zval(vector<T> const &cubes, vect_cube_t &out, float z1, float z2) {
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {
		if (i->z1() < z2 && i->z2() > z1) {out.push_back(*i);}
	}
}

void building_interior_t::get_avoid_cubes(vect_cube_t &avoid, float z1, float z2, float floor_thickness, bool same_as_player) const { // for AI
	avoid.clear();
	add_bcube_if_overlaps_zval(stairwells, avoid, z1, z2); // clearance not required
	add_bcube_if_overlaps_zval(elevators,  avoid, z1, z2); // clearance not required
	if (!room_geom) return; // no room objects
	auto objs_end(room_geom->objs.begin() + room_geom->stairs_start); // skip stairs and elevators

	for (auto c = room_geom->objs.begin(); c != objs_end; ++c) {
		// these object types are not collided with by people and can be skipped
		if (c->no_coll() || !(same_as_player ? bldg_obj_types[c->type].player_coll : bldg_obj_types[c->type].ai_coll)) continue;
		if (c->z1() < z2 && c->z2() > z1) {avoid.push_back(*c);}
	}
}

void building_t::find_nearest_stairs(point const &p1, point const &p2, vector<unsigned> &nearest_stairs, bool straight_only, int part_ix) const {
	nearest_stairs.clear();
	assert(interior);
	if (interior->stairwells.empty()) return; // no stairs
	assert(part_ix < 0 || (unsigned)part_ix < parts.size());
	float const zmin(min(p1.z, p2.z)), zmax(max(p1.z, p2.z));
	vector<pair<float, unsigned>> sorted;

	for (unsigned s = 0; s < interior->stairwells.size(); ++s) {
		stairwell_t const &stairs(interior->stairwells[s]);
		if (straight_only && stairs.shape == SHAPE_U) continue; // skip U-shaped stairs
		if (zmin < stairs.z1() || zmax > stairs.z2()) continue; // stairs don't span the correct floors
		if (part_ix >= 0 && !parts[part_ix].contains_cube(stairs)) continue; // stairs don't belong to this part
		point const center(stairs.get_cube_center());
		float const dist(p2p_dist(p1, center) + p2p_dist(center, p2));
		sorted.emplace_back(dist, s);
	} // for s
	sort(sorted.begin(), sorted.end()); // sort by distance, min first
	for (auto s = sorted.begin(); s != sorted.end(); ++s) {nearest_stairs.push_back(s->second);}
}

cube_t get_stairs_plus_step_up(stairwell_t const &stairs) {
	cube_t stairs_ext(stairs);
	stairs_ext.d[stairs.dim][!stairs.dir] += (stairs.dir ? -1.0 : 1.0)*stairs.get_sz_dim(stairs.dim)/NUM_STAIRS_PER_FLOOR; // location of step up
	return stairs_ext;
}

bool building_t::find_route_to_point(pedestrian_t const &person, float radius, bool is_first_path, bool is_moving_target, bool following_player, vector<point> &path) const {

	assert(interior && interior->nav_graph);
	point const &from(person.pos), &to(person.target_pos);
	path.clear();
	building_loc_t const loc1(get_building_loc_for_pt(from)), loc2(get_building_loc_for_pt(to));
	if (loc1.part_ix < 0 || loc2.part_ix < 0 || loc1.room_ix < 0 || loc2.room_ix < 0) return 0; // not in a room
	assert((unsigned)loc1.part_ix < parts.size() && (unsigned)loc2.part_ix < parts.size());
	assert((unsigned)loc1.room_ix < interior->rooms.size() && (unsigned)loc2.room_ix < interior->rooms.size());
	float const floor_spacing(get_window_vspace()), height(0.7*floor_spacing), z2_add(height - radius); // approximate, since we're not tracking actual heights
	static vect_cube_t avoid; // reuse across frames/people
	interior->get_avoid_cubes(avoid, (from.z - radius), (from.z + z2_add), get_floor_thickness(), following_player);

	if (loc1.same_room_floor(loc2)) { // same room/floor (not checking stairs_ix)
		assert(from.z == to.z);
		return interior->nav_graph->complete_path_within_room(from, to, loc1.room_ix, person.ssn, radius, person.cur_rseed, avoid, path);
	}
	if (loc1.floor_ix != loc2.floor_ix) { // different floors: find path from <from> to nearest stairs, then find path from stairs to <to>
		bool const straight_only = 0;
		vector<unsigned> nearest_stairs;
		find_nearest_stairs(from, to, nearest_stairs, straight_only); // pass in loc1.part_ix if both loc part_ix values are equal?
		bool const up_or_down(loc1.floor_ix > loc2.floor_ix); // 0=up, 1=down FIXME: handle part_ix

		for (auto s = nearest_stairs.begin(); s != nearest_stairs.end(); ++s) { // try using stairs, closest to furthest
			assert(*s < interior->stairwells.size());
			stairwell_t const &stairs(interior->stairwells[*s]);
			unsigned const stairs_room_ix(*s + interior->rooms.size()); // map to graph space
			path.clear();
			vector<point> from_path;
			// Note: passing use_stairs=0 here because it's unclear if we want to go through stairs nodes in our A* algorithm
			// from => stairs
			if (!interior->nav_graph->find_path_points(loc1.room_ix, stairs_room_ix, person.ssn, radius, height, 0, is_first_path,  up_or_down, person.cur_rseed, avoid, from, from_path)) continue;
			point const seg2_start(interior->nav_graph->get_stairs_entrance_pt(to.z, stairs_room_ix, !up_or_down)); // other end
			interior->get_avoid_cubes(avoid, (seg2_start.z - radius), (seg2_start.z + z2_add), get_floor_thickness(), following_player); // new floor, new zval, new avoid cubes
			// stairs => to
			if (!interior->nav_graph->find_path_points(stairs_room_ix, loc2.room_ix, person.ssn, radius, height, 0, is_first_path, !up_or_down, person.cur_rseed, avoid, seg2_start, path)) continue;
			assert(!path.empty() && !from_path.empty());
			path.push_back(seg2_start); // other end of the stairs
			// add two more points to straighten the entrance and exit paths; this segment doesn't check for intersection with stairs
			cube_t const stairs_ext(get_stairs_plus_step_up(stairs));
			point const stairs_enter(stairs_ext.closest_pt(from_path.front())), stairs_exit(stairs_ext.closest_pt(path.back()));
			path.push_back(stairs_exit);

			if (stairs.shape == SHAPE_U) { // add 2 extra points on mid-level landing; entrance and exit will be on the same side
				bool const dim(stairs.dim), dir(stairs.dir); // Note: see code in add_stairs_and_elevators()
				float const turn_pt(stairs.d[dim][dir] - 0.1*(dir ? 1.0 : -1.0)*stairs.get_sz_dim(dim)), seg_delta_z(0.45f*(to.z - from.z));
				point exit_turn(stairs_exit.x, stairs_exit.y, (to.z - seg_delta_z));
				exit_turn[dim] = turn_pt;
				path.push_back(exit_turn); // turning point for exit side of stairs
				point enter_turn(stairs_enter.x, stairs_enter.y, (from.z + seg_delta_z));
				enter_turn[dim] = turn_pt;
				path.push_back(enter_turn); // turning point for entrance side of stairs
			}
			path.push_back(stairs_enter);
			vector_add_to(from_path, path); // concatenate the two path segments in reverse order
			assert(!path.empty());
			return 1; // done/success
		} // for s
		path.clear(); // not necessary?
		return 0; // failed
	}
	assert(loc1.room_ix != loc2.room_ix);
	if (!interior->nav_graph->find_path_points(loc1.room_ix, loc2.room_ix, person.ssn, radius, height, 0, is_first_path, 0, person.cur_rseed, avoid, from, path)) return 0;
	assert(!path.empty());
	return 1;
}

void building_ai_state_t::next_path_pt(pedestrian_t &person, bool same_floor, bool starting_path) {
	assert(!path.empty());
	person.is_on_stairs = (!same_floor && !starting_path && person.target_pos.z != path.back().z);
	person.target_pos   = path.back();
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
		if (room.lit_by_floor && !room.is_lit_on_floor(floor_ix)) continue; // don't place person in an unlit room
		point pos;
		pos.z = room.z1() + fc_thick + window_vspacing*floor_ix;
		for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(room.d[d][0]+radius, room.d[d][1]-radius);} // random XY point inside this room
		cube_t bcube(pos);
		bcube.expand_by(radius); // expand more in Z?
		if (!is_valid_stairs_elevator_placement(bcube, radius)) continue;
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

bool can_ai_follow_player(pedestrian_t const &person) {
	if (!ai_follow_player()) return 0; // disabled
	if (!cur_player_building_loc.is_valid()) return 0; // no target
	if (cur_player_building_loc.building_ix != (int)person.dest_bldg) return 0; // wrong building
	if (player_in_closet >= 2) return 0; // ignore player if in the closet with the door closed
	return 1;
}
bool has_nearby_sound(pedestrian_t const &person, float floor_spacing) {
	if (!can_ai_follow_player(person)) return 0; // no need to track sounds
	point sound_pos; // unused
	return get_closest_building_sound(person.pos, sound_pos, floor_spacing);
}

bool building_t::same_room_and_floor_as_player(building_ai_state_t const &state, pedestrian_t const &person) const {
	return (cur_player_building_loc.room_ix == (int)state.cur_room && cur_player_building_loc.floor_ix == get_floor_for_zval(person.pos.z));
}
bool building_t::is_player_visible(building_ai_state_t const &state, pedestrian_t const &person, unsigned vis_test) const {
	if (vis_test == 0) return 1; // no visibility test
	building_dest_t const &target(cur_player_building_loc);
	float const player_radius(get_scaled_player_radius());
	point const pp2(target.pos - vector3d(0.0, 0.0, camera_zh)); // player's bottom sphere
	bool const same_room_and_floor(same_room_and_floor_as_player(state, person));

	if (!same_room_and_floor) { // check visibility; assume LOS if in the same room
		point const eye_pos(person.pos + vector3d(0.0, 0.0, (0.9*person.get_height() - person.radius))); // for person
		if (!is_sphere_visible(target.pos, player_radius, eye_pos) && !is_sphere_visible(pp2, player_radius, eye_pos)) return 0; // check both the bottom and top of player
	}
	if (vis_test >= 2) { // check person FOV
		//if (dot_product((pp2 - eye_pos).get_norm(), person.dir) < 0.5) return 0; // 60 degree FOV => dp < 0.5
		float const view_dist(10.0*get_window_vspace()); // 10 stories (~100 feet) should be good enough
		pos_dir_up pdu(person.pos, person.dir, plus_z, 0.0, 0.1*person.radius, view_dist, 0.0, 1); // auto perspective angle, no_zoom=1
		if (!pdu.sphere_visible_test(target.pos, player_radius) && !pdu.sphere_visible_test(pp2, player_radius)) return 0;
	}
	if (vis_test >= 3 && !same_room_and_floor) { // check lit state if not in the same room and floor
		if (!is_sphere_lit(target.pos, player_radius) && !is_sphere_lit(pp2, player_radius)) return 0; // check both the bottom and top of player
	}
	return 1;
}
bool building_t::can_target_player(building_ai_state_t const &state, pedestrian_t const &person) const {
	if (!can_ai_follow_player(person)) return 0;
	return is_player_visible(state, person, global_building_params.ai_player_vis_test); // 0=no test, 1=LOS, 2=LOS+FOV, 3=LOS+FOV+lit
}

bool building_t::need_to_update_ai_path(building_ai_state_t const &state, pedestrian_t const &person) const {
	if (!global_building_params.ai_target_player || !can_ai_follow_player(person) || !interior) return 0; // disabled
	building_dest_t const &target(cur_player_building_loc);
	bool const same_room(same_room_and_floor_as_player(state, person)); // check room and floor

	if (global_building_params.ai_player_vis_test == 0) { // the below logic only applies if there's no line of sight test and is an optimization
		// if the player's room has not changed, and the person is not yet in this room, continue on the same path to the dest room;
		// however, if the path is empty, continue to choose a new path (needed for AI more than one floor away from target to take multiple flights of stairs)
		if (target.same_room_floor(prev_player_building_loc) && !same_room && !state.on_new_path_seg && !state.path.empty()) return 0;
	}
	//if (same_room && state.path.size() > 1) return 0; // same room but path has a jog, continue on existing path (faster, but slower to adapt to player position change)
	if (dist_less_than(person.pos, target.pos, person.radius)) return 0; // already close enough
	if (person.on_stairs()) return 0; // don't change paths when on the stairs
	float const floor_spacing(get_window_vspace());

	if (int(target.pos.z/floor_spacing) != int(prev_player_building_loc.pos.z/floor_spacing)) { // if player did not change floors
		if (fabs(person.pos.z - target.pos.z) > 2.0f*floor_spacing) return 0; // person and player are > 2 floors apart, continue toward stairs (or should it be one floor apart?) (optimization)
	}
	if (can_target_player(state, person)) { // have player visibility
		if (state.goal_type == GOAL_TYPE_PLAYER && target.same_room_floor(prev_player_building_loc) && !same_room && !state.on_new_path_seg && !state.path.empty()) return 0;
		return 1;
	}
	if (has_nearby_sound(person, floor_spacing)) return 1; // new sound source
	return 0; // continue on the previously chosen path
}

// Note: non-const because this updates room light and door state
int building_t::ai_room_update(building_ai_state_t &state, rand_gen_t &rgen, vector<pedestrian_t> &people, float delta_dir, unsigned person_ix, bool stay_on_one_floor) {

	assert(person_ix < people.size());
	pedestrian_t &person(people[person_ix]);
	if (person.speed == 0.0) {person.anim_time = 0.0; return AI_STOP;} // stopped
	if (!interior->room_geom && frame_counter < 60) {person.anim_time = 0.0; return AI_WAITING;} // wait until room geom is generated for this building
	float const coll_dist(COLL_RADIUS_SCALE*person.radius);
	float &wait_time(person.waiting_start); // reuse this field
	float speed_mult(1.0);
	person.following_player = 0; // reset for this frame

	if (wait_time > 0) {
		if (wait_time > fticks && !can_ai_follow_player(person)) { // waiting; don't wait if there's a player to follow
			for (auto p = people.begin()+person_ix+1; p < people.end(); ++p) { // check for other people colliding with this person and handle it
				if (p->dest_bldg != person.dest_bldg) break; // done with this building
				if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
				float const rsum(coll_dist + COLL_RADIUS_SCALE*p->radius);
				if (!dist_xy_less_than(person.pos, p->pos, rsum)) continue; // not intersecting
				move_person_to_not_collide(person, *p, person.pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
			} // for p
			wait_time -= fticks;
			person.anim_time = 0.0; // reset just in case (though should already be at 0.0)
			state.goal_type  = GOAL_TYPE_NONE;
			return AI_WAITING;
		}
		wait_time = 0.0;
		person.target_pos = all_zeros; // force choose_dest=1 below
	}
	assert(interior);
	assert(bcube.contains_pt(person.pos)); // person must be inside the building
	build_nav_graph();

	if (can_ai_follow_player(person) && dist_less_than(person.pos, cur_player_building_loc.pos, 1.2f*(person.radius + get_scaled_player_radius()))) {
		register_ai_player_coll(person.has_key); // Note: returns is_dead, so we could track kills here
	}
	bool choose_dest(!person.target_valid());
	bool const update_path(need_to_update_ai_path(state, person));

	if (update_path) { // need to update based on player movement; higher priority than choose_dest
		if (choose_dest_goal(state, person, rgen, stay_on_one_floor) != 1 || // check if person can reach the target
			!find_route_to_point(person, coll_dist, 0, 1, 1, state.path)) // is_first_path=0, is_moving_target=1, following_player=1
		{
			choose_dest = 1; // or increment person.cur_rseed and return AI_WAITING? or restore person and state to prev values?
		}
		else { // success
			state.next_path_pt(person, stay_on_one_floor, 1);
			person.following_player = 1;
			choose_dest = 0;
			speed_mult  = 1.5; // faster when the player is in the same room
			// run logic to play zombie sounds
			bool const same_room_and_floor(same_room_and_floor_as_player(state, person));
			bool play_sound(same_room_and_floor); // always play sound if in the same room and floor
			if (!play_sound && (person_ix & 1)) {play_sound |= is_player_visible(state, person, 1);} // 50% of zombies use line of sight test
			if (!play_sound && (person_ix & 2)) {play_sound |= has_nearby_sound(person, get_window_vspace());} // 50% of zombies use sound test
			if (play_sound) {maybe_play_zombie_sound(person.pos, person_ix, same_room_and_floor);} // alert other zombies if in the same room and floor as the player
		}
	}
	if (choose_dest) { // no current destination, choose a new destination
		++person.cur_rseed; // update rseed so that this person will choose a different path if the previous path finding call failed
		if (!update_path) {person.anim_time = 0.0;} // reset animation if not a failed update path frame
		int const ret(choose_dest_room(state, person, rgen, stay_on_one_floor)); // 0=failed, 1=success, 2=failed but can retry

		if (ret == 2 && interior->door_state_updated) { // wait rather than stopping in case the player trapped this person in a room by closing the door
			person.wait_for(5.0); // stop for 5 seconds, then try again
			return AI_WAITING;
		}
		else if (ret != 1) { // if there's no valid room or valid path, set the speed to 0 so that we don't check this every frame
			if (!ai_follow_player()) { // if not following the player
				//person.speed = 0.0; // movement will be stopped from now on - not valid if we get into this state before gameplay is enabled
				person.wait_for(2.0); // stop for 2 seconds, then try again
			}
			return AI_STOP;
		}
		if (!find_route_to_point(person, coll_dist, state.is_first_path, 0, 0, state.path)) { // following_player=0, is_moving_target=0
			person.wait_for(1.0); // stop for 1 second, then try again
			return AI_WAITING;
		}
		state.is_first_path = 0;
		state.next_path_pt(person, stay_on_one_floor, 1);
		return AI_BEGIN_PATH;
	}
	float const max_dist(person.speed*speed_mult*fticks);
	state.on_new_path_seg = 0; // clear flag for this frame
	//person.following_player = can_target_player(state, person); // for debugging visualization

	if (dist_less_than(person.pos, person.target_pos, 1.1f*max_dist)) { // at dest
		assert(bcube.contains_pt(person.target_pos));
		person.pos = person.target_pos;
		if (!state.path.empty()) {state.next_path_pt(person, stay_on_one_floor, 0); return AI_NEXT_PT;} // move to next path point
		// don't wait if we can follow the player
		bool const no_wait(global_building_params.ai_target_player && (can_target_player(state, person) || has_nearby_sound(person, get_window_vspace())));
		if (!no_wait) {person.wait_for(rgen.rand_uniform(1.0, 10.0));} // stop for 1-10 seconds
		state.on_new_path_seg = 1; // allow player following AI update logic to rerun this frame
		return AI_AT_DEST;
	}
	// this person is walking to a destination point
	vector3d const new_dir(person.target_pos - person.pos);
	float const new_dir_mag(new_dir.mag());
	point new_pos;

	if (dot_product(new_dir, person.dir) < 0.999*new_dir_mag) { // dir not perfectly aligned
		//if (person.is_close_to_player()) {cout << TXT(new_dir.str()) << TXT(person.dir.str()) << TXT(new_dir_mag) << TXT(delta_dir) << TXT(max_dist) << TXT(person.radius) << endl;}
		assert(new_dir_mag > TOLERANCE); // should be guaranteed by dist_less_than() test, assuming zvals are equal (which they should be)
		float const step_scale(max(0.1f, dot_product(person.dir, new_dir)/new_dir_mag)); // move more slowly when direction misaligns to avoid overshooting target_pos
		person.dir = (delta_dir/new_dir_mag)*new_dir + (1.0 - delta_dir)*person.dir; // merge new_dir into dir gradually for smooth turning
		if (person.on_stairs()) {person.dir.z = new_dir.z/new_dir_mag;} // dir.z tracks exactly
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
	if (!person.on_stairs()) {
		for (auto p = people.begin()+person_ix+1; p < people.end(); ++p) { // check all other people in the same building after this one and attempt to avoid them
			if (p->dest_bldg != person.dest_bldg) break; // done with this building
			if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
			float const rsum(coll_dist + COLL_RADIUS_SCALE*p->radius);
			if (!dist_xy_less_than(new_pos, p->pos, rsum)) continue; // new pos not close
			if (!dist_xy_less_than(person.pos, p->pos, rsum)) return AI_STOP; // old pos not intersecting, stop
			person.anim_time = 0.0; // pause animation in case this person is mid-step
			move_person_to_not_collide(person, *p, new_pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
			return AI_MOVING; // return here, but don't update animation or dir; only handles a single collision
		} // for p
	}
	bool const player_in_this_building(cur_player_building_loc.building_ix == (int)person.dest_bldg); // basement door only counts if the player is in this building
	bool const might_have_closed_door(global_building_params.open_door_prob < 1.0 || (player_in_this_building && has_basement()));

	if (interior->door_state_updated || (global_building_params.ai_opens_doors == 2 && might_have_closed_door)) {
		// check for any doors that started open or the player has closed;
		// this can be slow, so we only enable it for buildings where the player changed the door state, or when the AI can open all doors
		for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
			if (i->open) continue; // doors tend to block the player, don't collide with them unless they're closed
			if (new_pos.z < i->z1() || new_pos.z > i->z2())         continue; // wrong part/floor
			if (!sphere_cube_intersect(new_pos, person.radius, *i)) continue; // no intersection with door
			cube_t door(*i);
			door.expand_in_dim(i->dim, 0.5*get_wall_thickness()); // increase door thickness to a nonzero value
			if (!door.line_intersects(person.pos, person.target_pos)) continue; // check if our path goes through the door, to allow for "glancing blows" when pushed or turning

			if (global_building_params.ai_opens_doors && (!i->is_closed_and_locked() || person.has_key)) { // can open the door
				toggle_door_state((i - interior->doors.begin()), player_in_this_building, 0, person.pos.z); // by_player=0
			}
			else { // can't open the door
				person.wait_for(5.0); // wait for 5s and then choose a new desination
				return AI_WAITING; // cut the path short at this closed door
			}
		} // for i
	}
	// logic to clip this person to correct room Z-bounds in case something went wrong; remove if/when this is fixed
	// Note: we probably can't use the room Z bounds here becase the person may be on the stairs connecting two stacked parts
	float const min_valid_zval(bcube.z1() + 0.5f*get_floor_thickness() + person.radius), max_valid_zval(bcube.z2() - person.radius);

	if (player_in_this_building && !person.on_stairs()) { // movement in XY, not on stairs, snap to nearest floor
		// this is optional and is done just in case something went wrong
		assert(state.cur_room < interior->rooms.size());
		room_t const &room(interior->rooms[state.cur_room]);
		float const floor_spacing(get_window_vspace());
		int cur_floor(max(0, round_fp((new_pos.z - min_valid_zval)/floor_spacing)));
		min_eq(cur_floor, (round_fp(room.dz()/floor_spacing) - 1)); // clip to the valid floors for this room
		float const adj_zval(cur_floor*floor_spacing + min_valid_zval);
		if (fabs(adj_zval - new_pos.z) > 0.1*person.radius) {person.target_pos = all_zeros; state.path.clear();} // if we snap to the floor, reset the target and path
		new_pos.z = adj_zval;
	}
	max_eq(new_pos.z, min_valid_zval); // don't let the person go below the ground floor
	min_eq(new_pos.z, max_valid_zval); // don't let the person go above the room ceiling
	// update state
	person.pos        = new_pos; // Note: new_pos.z should equal person.poz.z unless on stairs, which is difficult to accurately check for in this function
	person.anim_time += max_dist;
	if (player_in_this_building) {state.cur_room = get_room_containing_pt(person.pos);} // update cur_room after moving
	ai_room_lights_update(state, person, people, person_ix); // non-const part
	return AI_MOVING;
}

void building_t::ai_room_lights_update(building_ai_state_t const &state, pedestrian_t const &person, vector<pedestrian_t> const &people, unsigned person_ix) {
	if (!(display_mode & 0x20)) return; // disabled by default, enable with key '6'
	if (ai_follow_player() && global_building_params.ai_player_vis_test >= 3) return; // if AI tests that the player is lit, then we shouldn't be turning on and off the lights
	if ((frame_counter + person_ix) & 7) return; // update room info only every 8 frames
	int const room_ix(get_room_containing_pt(person.pos));
	if (room_ix < 0) return; // room is not valid (between rooms, etc.)
	room_t const &cur_room(get_room(room_ix));
	set_room_light_state_to(cur_room, person.pos.z, 1); // make sure current room light is on
	if ((unsigned)room_ix == state.cur_room) return; // same room as last time - done
	assert(state.cur_room < interior->rooms.size());
	bool other_person_in_room(0);

	// check for other people in the room before turning the lights off on them
	for (unsigned i = person_ix; i-- ;) { // look behind
		pedestrian_t const &p(people[i]);
		if (p.dest_bldg != person.dest_bldg) break; // done with this building
		if (get_room_containing_pt(p.pos) == (int)state.cur_room && fabs(person.pos.z - p.pos.z) < get_window_vspace()) {other_person_in_room = 1; break;}
	}
	for (unsigned i = person_ix+1; i < people.size() && !other_person_in_room; ++i) { // look ahead
		pedestrian_t const &p(people[i]);
		if (p.dest_bldg != person.dest_bldg) break; // done with this building
		if (get_room_containing_pt(p.pos) == (int)state.cur_room && fabs(person.pos.z - p.pos.z) < get_window_vspace()) {other_person_in_room = 1; break;}
	}
	if (!other_person_in_room) {set_room_light_state_to(get_room(state.cur_room), person.pos.z, 0);} // make sure old room light is off
}

void building_t::move_person_to_not_collide(pedestrian_t &person, pedestrian_t const &other, point const &new_pos, float rsum, float coll_dist) const {
	point const other_pos(other.pos.x, other.pos.y, person.pos.z); // use same zval to ignore height differences
	float const sep_dist(p2p_dist_xy(person.pos, other_pos)), move_dist(rsum - sep_dist); // distance we have to move
	// move away from the other person, hopefully not through a wall
	point const orig_pos(person.pos);
	if (sep_dist > 0.01*coll_dist) {person.pos += (move_dist/sep_dist)*(person.pos - other_pos);}
	else {person.pos.x += rsum;} // avoid divide-by-zero, choose +X direction arbitrarily
	int const room_ix(get_building_loc_for_pt(orig_pos).room_ix);
	cube_t clip_bounds((room_ix >= 0 && (unsigned)room_ix < interior->rooms.size()) ? get_room(room_ix) : bcube);
	clip_bounds.expand_by_xy(-coll_dist); // shrink
	clip_bounds.union_with_pt(orig_pos); // we know this point was valid
	clip_bounds.union_with_pt(new_pos);  // we know this point is valid
	clip_bounds.clamp_pt_xy(person.pos); // force player into the room

	if (!bcube.contains_pt(person.pos)) { // this can happen on rare occasions, due to fp inaccuracy or multiple collisions
		//cout << TXT(rsum) << TXT(sep_dist) << TXT(move_dist) << TXT(room_ix) << TXT(other.pos.str()) << TXT(person.pos.str()) << TXT(bcube.str()) << endl;
		bcube.clamp_pt_xy(person.pos); // just clamp pos so that it doesn't assert later
	}
}

// Note: non-const because this updates room lights
void vect_building_t::ai_room_update(vector<building_ai_state_t> &ai_state, vector<pedestrian_t> &people, float delta_dir, rand_gen_t &rgen) {
	//timer_t timer("Building People Update"); // ~3.7ms for 50K people, 0.55ms with distance check
	point const camera_bs(get_camera_building_space());
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

int building_t::get_room_containing_pt(point const &pt) const {
	assert(interior);
	float const wall_thickness(get_wall_thickness());

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		cube_t tc(*r);
		tc.expand_by_xy(wall_thickness); // to include point in doorway
		if (tc.contains_pt(pt)) {return (r - interior->rooms.begin());}
	}
	return -1; // room not found
}
building_loc_t building_t::get_building_loc_for_pt(point const &pt) const {
	building_loc_t loc;

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		if (p->contains_pt(pt)) {loc.part_ix = (p - parts.begin()); break;}
	}
	if (interior) { // rooms and stairwells, no elevators yet
		loc.room_ix = get_room_containing_pt(pt);
		if (loc.room_ix >= 0) {loc.floor_ix = get_floor_for_zval(pt.z);}

		for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
			if (s->contains_pt(pt)) {loc.stairs_ix = (s - interior->stairwells.begin()); break;} // Note: stairs_ix is not currently used, except for >= 0 test
		}
	}
	return loc;
}
bool building_t::room_containing_pt_has_stairs(point const &pt) const {
	int const room_ix(get_room_containing_pt(pt));
	if (room_ix < 0) return 0; // no room contains this point
	return get_room(room_ix).has_stairs;
}

// these must be here to handle deletion of building_nav_graph_t, which is only defined in this file
building_interior_t::building_interior_t() : top_ceilings_mask(0), door_state_updated(0), is_unconnected(0) {}
building_interior_t::~building_interior_t() {}
