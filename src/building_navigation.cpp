// 3D World - Building Navigation for AIs
// by Frank Gennari 3/14/20

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for person_t
#include <queue>


bool  const ENABLE_AI_ELEVATORS  = 0;
unsigned const ELEVATOR_CAPACITY = 1; // number of people that can use the elevator at once; nonzero
float const COLL_RADIUS_SCALE    = 0.75; // somewhat smaller than radius, but larger than PED_WIDTH_SCALE
float const RETREAT_TIME         = 4.0f*TICKS_PER_SECOND; // 4s

int cpbl_update_frame(0);
building_dest_t cur_player_building_loc, prev_player_building_loc;

extern bool player_is_hiding;
extern int frame_counter, display_mode, animate2;
extern float fticks;
extern building_params_t global_building_params;
extern bldg_obj_type_t bldg_obj_types[];

bool in_building_gameplay_mode();
bool ai_follow_player() {return (global_building_params.ai_follow_player || in_building_gameplay_mode());}
bool can_ai_follow_player(person_t const &person);
float get_closest_building_sound(point const &at_pos, point &sound_pos, float floor_spacing);
void maybe_play_zombie_sound(point const &sound_pos_bs, unsigned zombie_ix, bool alert_other_zombies, bool high_priority=0);
int register_ai_player_coll(bool &has_key, float height);

point get_cube_center_zval(cube_t const &c, float zval) {return point(c.xc(), c.yc(), zval);}

// Note: this should go into building_t/buildings.h at some point, but is temporarily here
class building_nav_graph_t {
	struct conn_room_t { // size=16
		unsigned ix; // node index
		int door_ix; // index of connecting door on first floor; -1 for connections with out a door
		vector2d pt[2]; // stairs entry/exit point {up, down}
		conn_room_t(unsigned ix_, int dix, vector2d const &ptu, vector2d const &ptd) : ix(ix_), door_ix(dix) {pt[0] = ptu; pt[1] = ptd;}
	};
	struct node_t { // represents one room or one stairwell
		bool has_exit=0, is_hallway=0, is_stairs=0, is_ramp=0; // has_exit is not yet used
		cube_t bcube;
		vector<conn_room_t> conn_rooms;

		point get_center(float zval) const {return get_cube_center_zval(bcube, zval);}
		bool is_vert_conn() const {return (is_stairs || is_ramp);} // or is_elevator if added later

		void add_conn_room(unsigned room, int door_ix, cube_t const &cu, cube_t const &cd) {
			for (auto i = conn_rooms.begin(); i != conn_rooms.end(); ++i) {if (i->ix == room) return;} // ignore duplicates
			conn_rooms.emplace_back(room, door_ix, vector2d(cu.xc(), cu.yc()), vector2d(cd.xc(), cd.yc()));
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
	bool has_pg_ramp;
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

	void set_num_rooms(unsigned num_rooms_, unsigned num_stairs_, bool has_pg_ramp_) {
		num_rooms   = num_rooms_;
		num_stairs  = num_stairs_;
		has_pg_ramp = has_pg_ramp_;
		nodes.resize(num_rooms + num_stairs + has_pg_ramp);
		for (unsigned n = num_rooms; n < (num_rooms + num_stairs); ++n) {nodes[n].is_stairs = 1;}
		if (has_pg_ramp) {nodes.back().is_ramp = 1;}
	}
	void set_room_bcube  (unsigned room,   cube_t const &c) {get_node(room).bcube = c;}
	void set_stairs_bcube(unsigned stairs, cube_t const &c) {get_node(stairs + num_rooms).bcube = c;}
	void set_ramp_bcube                   (cube_t const &c) {get_node(num_stairs + num_rooms).bcube = c;}
	void mark_hallway(unsigned room) {get_node(room).is_hallway = 1;}
	void mark_exit   (unsigned room) {get_node(room).has_exit   = 1;}

	void connect_stairs(unsigned room, unsigned stairs, bool dim, bool dir, bool is_u, bool is_ramp) {
		assert(room < num_rooms && (is_ramp || stairs < num_stairs));
		unsigned const node_ix2(num_rooms + (is_ramp ? num_stairs : stairs)); // pg_ramp comes after all stairs
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
		get_node(room).add_conn_room(node_ix2, -1, entry_u, entry_d);
		n2.add_conn_room(room, -1, entry_u, entry_d);
	}
	void connect_rooms(unsigned room1, unsigned room2, int door_ix, cube_t const &conn_bcube) { // graph is bidirectional
		assert(room1 < num_rooms && room2 < num_rooms);
		get_node(room1).add_conn_room(room2, door_ix, conn_bcube, conn_bcube);
		get_node(room2).add_conn_room(room1, door_ix, conn_bcube, conn_bcube);
	}
	void disconnect_room_pair(unsigned room1, unsigned room2) { // remove connections in both directions
		assert(room1 != room2 && room1 < num_rooms && room2 < num_rooms);
		remove_connection(room1, room2);
		remove_connection(room2, room1);
	}
	bool is_room_connected_to(unsigned room1, unsigned room2, vect_door_t const &doors, float zval, bool has_key) const {
		// Note: likely faster than running full A* algorithm
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
				if (seen[i->ix]) continue;
				if (!can_use_conn(*i, doors, zval, has_key)) continue; // blocked by closed or locked door
				if (i->ix == room2) return 1; // found, done
				pend.push_back(i->ix);
				seen[i->ix] = 1;
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

	static bool is_valid_pos(vect_cube_t const &avoid, point const &pos, float radius, float height) { // Note: assumes zvals are already checked
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
	static bool find_valid_pt_in_room(vect_cube_t const &avoid, float radius, float height, float zval, cube_t const &room, rand_gen_t &rgen, point &pos, bool no_use_init=0) {
		if (!no_use_init && avoid.empty()) return 1; // no colliders, any point is valid, choose the initial point (which should be the room center)
		cube_t place_area(room);
		place_area.expand_by_xy(-2.0*radius); // shrink by twice the radius
		if (!place_area.is_strictly_normalized()) return 0; // should generally not be true

		if (no_use_init) { // chose a new initial point
			choose_pt_xy_in_room(pos, place_area, rgen);
			if (avoid.empty()) return 1;
		}
		point orig_pos(pos);

		for (unsigned n = 0; n < 100; ++n) { // 100 random tries to find a valid dest_pos
			if (is_valid_pos(avoid, pos, radius, height)) return 1; // success
			choose_pt_xy_in_room(pos, place_area, rgen); // choose a random new point in the room
		}
		pos = orig_pos; // use orig value as failed point
		return 0;
	}
	point find_valid_room_dest(vect_cube_t const &avoid, float radius, float height, float zval, unsigned node_ix,
		bool up_or_down, bool &not_room_center, rand_gen_t &rgen, bool no_use_init, point const *const custom_dest) const
	{
		node_t const &node(get_node(node_ix));
		if (node.is_vert_conn()) {return get_stairs_entrance_pt(zval, node_ix, up_or_down);}
		point pos((custom_dest != nullptr) ? *custom_dest : get_cube_center_zval(node.bcube, zval)); // first candidate is the center of the room, if not custom
		if (find_valid_pt_in_room(avoid, radius, height, zval, node.bcube, rgen, pos, no_use_init)) {not_room_center = 1;} // success
		return pos;
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
		} // for i
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
		if (!node.is_vert_conn()) {walk_area.expand_by_xy(-radius);} // shrink by radius if not stairs or a ramp

		if (node.is_hallway) {
			bool const min_dim(walk_area.dy() < walk_area.dx());
			float const shrink(min(0.5f*radius, max(0.0f, 0.9f*0.5f*walk_area.get_sz_dim(min_dim)))); // shrink by an extra half radius if hallway is wide enough
			walk_area.expand_in_dim(min_dim, -shrink);
		}
		return walk_area;
	}
	bool reconstruct_path(vector<a_star_node_state_t> const &state, vect_cube_t const &avoid, point const &cur_pt, float radius,
		float height, unsigned start_ix, unsigned end_ix, unsigned ped_ix, bool is_first_path, bool up_or_down, unsigned ped_rseed,
		point const *const custom_dest, vector<point> &path) const
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
					bool const no_use_init(n > 0); // choose a random new point on iterations after the first one
					bool not_room_center(0);
					point const end_point(find_valid_room_dest(avoid, radius, height, cur_pt.z, start_ix, up_or_down, not_room_center, rgen, no_use_init, custom_dest));
					path.push_back(end_point);
					if (node.is_vert_conn()) {success = 1; break;} // done, don't need to run code below
					point const room_exit(closest_room_pt(walk_area, next)); // first doorway
					if (connect_room_endpoints(avoid, walk_area, end_point, room_exit, radius, path, keepout, rgen)) {path.push_back(room_exit); success = 1; break;}
					path.clear(); // failed, reset for next iteration
					if (!not_room_center) break; // if we did choose the room center, and there is no path to it, we've failed
				} // for n
				if (!success) {assert(path.empty()); return 0;} // failed to connect to a point in dest room
			}
			else if (came_from < 0) { // done (next is not valid here)
				assert(n == end_ix);
				if (node.is_vert_conn()) return 1; // success
				point const final_pt(closest_room_pt(walk_area, path.back())); // walk from room into last doorway
				path.push_back(final_pt);

				// find path to first doorway; ignore collisions with p2 (cur_pt) in case this person was pushed into an object by another person
				if (!connect_room_endpoints(avoid, walk_area, final_pt, cur_pt, radius, path, keepout, rgen, 0, 1)) {
					// allow a failure if this is the first path taken by this AI so that it's not stuck behind an object due to bad initial placement
					if (!is_first_path) {path.clear(); return 0;}
				}
				return 1; // success
			}
			else if (!node.is_vert_conn()) { // adjust the path through a room
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
	bool complete_path_within_room(point const &from, point const &to, unsigned room, unsigned ped_ix, float radius,
		unsigned ped_rseed, bool is_first_path, bool following_player, vect_cube_t const &avoid, vector<point> &path) const
	{
		// used for reaching a goal such as the player within the same room
		cube_t const walk_area(calc_walkable_room_area(get_node(room), radius));
		rand_gen_t rgen;
		rgen.set_state((ped_ix + 13*room + 1), (ped_rseed + 1));
		vect_cube_t keepout;
		path.push_back(to); // Note: path is constructed backwards, so "to" is added first and connect_room_endpoints takes swapped arguments
		
		// ignore starting collisions, for example collisions with stairwell when exiting stairs?
		// ignore initial coll with "from", and coll with "to" when following the player
		if (!connect_room_endpoints(avoid, walk_area, to, from, radius, path, keepout, rgen, 1, following_player)) { // ignore_p1_coll (to) = 1
			if (!is_first_path) { // ignore failure on first path to allow person to get out from an object they spawn in
				bool success(0);

				if (following_player) {
					// try to find a partial path, starting at "to" and working toward "from"; since rooms are rectangular, all points on the line will be contained
					for (unsigned n = 1; n <= 9; ++n) { // {10% ... 90%}
						point const new_to(to + (float(n)/10)*(from - to));
						assert(walk_area.contains_pt(new_to));
						if (connect_room_endpoints(avoid, walk_area, new_to, from, radius, path, keepout, rgen, 0, following_player)) {success = 1; break;} // ignore_p1_coll = 0
					} // for n
				}
				if (!success) {path.clear(); return 0;}
			}
		}
		// maybe add an extra path point to prevent clipping through walls when walking through a doorway
		point const from_extend_pt(closest_room_pt(walk_area, from));
		if (from_extend_pt != from) {path.push_back(from_extend_pt);} // add if p1 was clamped
		return 1;
	}

	bool can_use_conn(conn_room_t const &conn, vect_door_t const &doors, float zval, bool has_key) const {
		if (conn.door_ix < 0) return 1; // no door
		if (global_building_params.ai_opens_doors && has_key) return 1; // locked door won't stop us
		unsigned const door_ix(conn.door_ix);
		assert(door_ix < doors.size());
		door_t const &first_door(doors[door_ix]);

		for (auto i = (doors.begin() + door_ix); i != doors.end(); ++i) {
			if (!i->is_same_stack(first_door)) break; // we've reached the end of the vertical stack of doors
			if (zval < i->z1() || zval > i->z2()) continue; // not the correct floor
			return (i->open || (global_building_params.ai_opens_doors && (!i->locked || has_key)));
		}
		return 1; // Note: we can get here for complex floorplan office buildings with bad interior walls (-4.18, 4.28, -3.46)
	}
	
	// A* algorithm; Note: path is stored backwards
	bool find_path_points(unsigned room1, unsigned room2, unsigned ped_ix, float radius, float height, bool use_stairs, bool is_first_path,
		bool up_or_down, unsigned ped_rseed, vect_cube_t const &avoid, point const &cur_pt, vect_door_t const &doors, bool has_key,
		point const *const custom_dest, vector<point> &path) const
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
				if (conn_node.is_vert_conn() && !use_stairs && i->ix != room2) continue; // skip stairs/ramp in this mode
				point const conn_center(conn_node.get_center(cur_pt.z));
				a_star_node_state_t &sn(state[i->ix]);
				vector2d const &pt(i->pt[up_or_down]);
				float const new_g_score(sn.g_score + p2p_dist_xy(center, pt) + p2p_dist_xy(pt, conn_center));
				if (!open[i->ix]) {open[i->ix] = 1;}
				else if (new_g_score >= sn.g_score) continue; // not better
				if (!can_use_conn(*i, doors, cur_pt.z, has_key)) continue; // blocked by closed or locked door
				sn.came_from_ix = cur;
				sn.path_pt.assign(pt.x, pt.y, cur_pt.z);
				
				if (i->ix == room2) { // done, reconstruct path (in reverse)
					return reconstruct_path(state, avoid, cur_pt, radius, height, i->ix, room1, ped_ix, is_first_path, up_or_down, ped_rseed, custom_dest, path);
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
	// Note: regular house rooms start and end at the walls;
	// offices, hallways, and extended basement rooms tile exactly and include half the walls, so we have to subtract those back off
	if (room.inc_half_walls()) {c.expand_by_xy(-0.5*get_wall_thickness());}
	return c;
}

void building_t::build_nav_graph() const { // Note: does not depend on room geom

	assert(interior);
	if (interior->nav_graph && !interior->nav_graph->invalid) return; // already built
	interior->nav_graph.reset(new building_nav_graph_t(0.5*get_window_vspace())); // set stairs_extend == doorway width
	building_nav_graph_t &ng(*interior->nav_graph);
	float const wall_width(get_wall_thickness());
	unsigned const num_rooms(interior->rooms.size()), num_stairs(interior->stairwells.size());
	bool const has_ramp(has_pg_ramp());
	ng.set_num_rooms(num_rooms, num_stairs, has_ramp);
	for (unsigned s = 0; s < num_stairs; ++s) {ng.set_stairs_bcube(s, interior->stairwells[s]);}
	if (has_ramp) {ng.set_ramp_bcube(interior->pg_ramp);}

	for (unsigned r = 0; r < num_rooms; ++r) {
		room_t const &room(interior->rooms[r]);
		cube_t c(get_walkable_room_bounds(room));
		if (room.is_hallway) {ng.mark_hallway(r);}
		ng.set_room_bcube(r, c);
		c.expand_by_xy(wall_width); // to include adjacent doors
		if (is_room_adjacent_to_ext_door(c)) {ng.mark_exit(r);}

		for (auto d = interior->door_stacks.begin(); d != interior->door_stacks.end(); ++d) {
			// if SPLIT_DOOR_PER_FLOOR, we should only add this door if it's on the same floor as our graph but that doesn't work because the graph is shared across all floors,
			// so instead we'll have to record the door index and check the correct door during path finding; it's not valid to test door open/locked state here
			if (!c.intersects_no_adj(*d)) continue; // door not adjacent to this room
			cube_t dc(*d);
			dc.expand_by_xy(wall_width); // to include adjacent rooms
			assert(d->first_door_ix < interior->doors.size());

			for (unsigned r2 = r+1; r2 < num_rooms; ++r2) { // check rooms with higher index (since graph is bidirectional)
				if (dc.intersects_no_adj(interior->rooms[r2])) {ng.connect_rooms(r, r2, d->first_door_ix, *d); break;}
			}
		} // for d
		for (unsigned s = 0; s < num_stairs; ++s) { // stairs
			stairwell_t const &stairwell(interior->stairwells[s]);

			if (stairwell.stairs_door_ix >= 0 && global_building_params.ai_opens_doors < 2) { // check for open doors; doors on stairs can't be locked
				if (!get_door(stairwell.stairs_door_ix).open) continue; // stairs blocked by closed door, don't connect (even if unlocked)
			}
			if (room.intersects_no_adj(stairwell)) {ng.connect_stairs(r, s, stairwell.dim, stairwell.dir, (stairwell.shape == SHAPE_U), 0);} // is_ramp=0
		}
		if (room.is_hallway) { // check for connected hallways
			for (unsigned r2 = r+1; r2 < num_rooms; ++r2) { // check rooms with higher index (since graph is bidirectional)
				room_t const &room2(interior->rooms[r2]);
				if (!room2.is_hallway || room2.z1() != room.z1() || !room2.intersects(c)) continue; // not a connected hallway
				cube_t conn_cube(c);
				conn_cube.intersect_with_cube(room2);
				ng.connect_rooms(r, r2, -1, conn_cube);
			} // for r2
		}
		if (room.get_room_type(0) == RTYPE_PARKING && has_ramp && room.intersects_no_adj(interior->pg_ramp)) { // include parking garage ramp
			bool const dim(interior->pg_ramp.ix >> 1), dir(interior->pg_ramp.ix & 1);
			ng.connect_stairs(r, 0, dim, dir, 0, 1); // stairs_ix=0, is_u=0, is_ramp=1
		}
		//for (unsigned e = 0; e < interior->elevators.size(); ++e) {} // elevators are ignored here by the AI; should check interior->elevators_disabled
	} // for r
	//if (is_house && !has_sec_bldg() && has_basement() && !ng.is_fully_connected()) {cout << "bcube " << bcube.str() << endl;}
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
		if (room_exp.contains_pt(d->get_bcube().get_cube_center())) return 1;
		if (front_door_only) return 0; // assumes the first door is the front door
	}
	return 0;
}

// Note: this is somewhat slow, should we build and use the nav graph? maybe not since this is only called for room assiggnment on the first floor of houses
bool building_t::are_rooms_connected_without_using_room(unsigned room1, unsigned room2, unsigned room_exclude) const {
	unsigned const num_rooms(interior->rooms.size());
	assert(room1 < num_rooms && room2 < num_rooms);
	assert(room_exclude != room1 && room_exclude != room2 && room_exclude < num_rooms);
	if (room1 == room2) return 1;
	float const wall_width(get_wall_thickness());
	bool const use_bit_mask(num_rooms <= 64); // almost always true
	static vector<unsigned> pend; // reused across calls
	pend.clear();
	pend.push_back(room1);
	static vector<uint8_t> seen; // reused across calls
	uint64_t seen_mask(0);

	if (use_bit_mask) {seen_mask |= ((1ULL << room1) | (1ULL << room_exclude));}
	else { // must use seen vector
		seen.clear();
		seen.resize(num_rooms, 0);
		seen[room1       ] = 1;
		seen[room_exclude] = 1; // mark as seen so that we won't visit it
	}
	while (!pend.empty()) { // flood fill - maybe A* is better, but that's a lot of work
		unsigned const cur(pend.back());
		pend.pop_back();
		cube_t room(interior->rooms[cur]);
		room.expand_by_xy(wall_width); // to include adjacent doors

		for (auto d = interior->door_stacks.begin(); d != interior->door_stacks.end(); ++d) {
			if (!room.intersects_no_adj(*d)) continue; // door not adjacent to this room
			cube_t dc(*d);
			dc.expand_by_xy(wall_width); // to include adjacent rooms

			for (unsigned r = 0; r < num_rooms; ++r) {
				if ((use_bit_mask ? (seen_mask & (1ULL << r)) : seen[r]) || !dc.intersects_no_adj(interior->rooms[r])) continue;
				if (r == room2) return 1; // found, done
				pend.push_back(r);
				if (use_bit_mask) {seen_mask |= (1ULL << r);} else {seen[r] = 1;}
			}
		} // for d
	} // end while()
	return 0; // room2 not found
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

bool building_t::choose_dest_goal(person_t &person, rand_gen_t &rgen) const { // used for following the player in gameplay mode

	assert(interior && interior->nav_graph);
	building_loc_t const loc(get_building_loc_for_pt(person.pos));
	if (loc.room_ix < 0) return 0; // not in a room
	float const floor_spacing(get_window_vspace());
	building_dest_t goal;
	point sound_pos;

	if ((global_building_params.ai_target_player || cur_player_building_loc.same_room_floor(loc)) && can_target_player(person)) {
		goal = cur_player_building_loc; // player is in a different room of our building, or we're following the player's position
		person.goal_type = GOAL_TYPE_PLAYER;
	}
	else if (can_ai_follow_player(person) && get_closest_building_sound(person.pos, sound_pos, floor_spacing)) { // target the loudest sound
		goal = building_dest_t(get_building_loc_for_pt(sound_pos), sound_pos, cur_player_building_loc.building_ix); // same building as player (current building)
		person.goal_type = GOAL_TYPE_SOUND;
	}
	if (!goal.is_valid()) return 0; // player or sound
	unsigned const cand_room(goal.room_ix);
	room_t const &room(get_room(cand_room)); // target room
	if (!interior->nav_graph->is_room_connected_to(loc.room_ix, cand_room, interior->doors, person.pos.z, person.has_key)) return 0; // unreachable
	person.cur_room     = loc.room_ix;
	person.dest_room    = cand_room; // set but not yet used
	person.target_pos   = (global_building_params.ai_target_player ? goal.pos: get_center_of_room(cand_room));
	person.target_pos.z = person.pos.z; // keep orig zval to stay on the same floor
	person.is_on_stairs = 0;

	// allow moving to a different floor, currently only one floor at a time
	if (goal.floor_ix < loc.floor_ix) { // try one floor below
		float const new_z(person.target_pos.z - floor_spacing);
		if (new_z > room.z1()) {person.target_pos.z = new_z;} // change if there is a floor below
	}
	else if (goal.floor_ix > loc.floor_ix) { // try one floor above
		float const new_z(person.target_pos.z + floor_spacing);
		if (new_z < room.z2()) {person.target_pos.z = new_z;} // change if there is a floor above
	}
	// else if floors differ by more than 1, we'll end up visiting the room on the wrong floor
	
	if (global_building_params.ai_target_player) { // ensure target is a valid location in this building; this must be done *after* adjacent floor zval adjustment
		// handle the case where the player is standing on the stairs on the same floor by moving zval to a different floor to force this person to use the stairs; what about pg_ramp?
		if (person.goal_type == GOAL_TYPE_PLAYER && loc.floor_ix == goal.floor_ix && goal.stairs_ix >= 0) {
			float const person_z1(person.get_z1()), player_z1(cur_player_building_loc.pos.z - CAMERA_RADIUS - get_player_height()), fc_thick(get_fc_thickness());
			// make destination exactly one floor above or below of where we currently are; some hysteresis is required to handle the case where the player is at the same zval
			if      (player_z1 + fc_thick < person_z1) {person.target_pos.z -= floor_spacing;} // move down one floor
			else if (player_z1 - fc_thick > person_z1) {person.target_pos.z += floor_spacing;} // move up   one floor
		}
		float const z2_add(person.get_height() - person.radius), coll_dist(COLL_RADIUS_SCALE*person.radius);
		cube_t ai_bcube;
		if (has_basement() && person.target_pos.z < ground_floor_z1) {ai_bcube = get_full_basement_bcube();} // in the basement
		else {ai_bcube = bcube;} // above ground
		cube_t legal_area(ai_bcube);
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
		if (has_ext_basement()) {
			float const dsq(p2p_dist_sq(person.target_pos, interior->basement_ext_bcube.closest_pt(person.target_pos)));
			if (dsq < dmin_sq) {closest_part = interior->basement_ext_bcube;}
		}
		if (dmin_sq > 0.0 && !closest_part.is_all_zeros()) {closest_part.clamp_pt(person.target_pos);} // clamp to closest part
		static vect_cube_t avoid; // reuse across frames/people
		// same_as_player=1, skip_stairs=1
		interior->get_avoid_cubes(avoid, (person.target_pos.z - person.radius), (person.target_pos.z + z2_add), 0.5*person.radius, get_floor_thickness(), 1, 1);

		// check for initial collisions at the player's location, but exclude stairs in case the player is standing on them;
		// this may no longer be required since complete_path_within_room() now ignores initial collisions with the dest, but is likely still a good idea
		for (unsigned n = 0; n < 4; ++n) { // iterate a few times in case a collision moves pos into another object
			bool any_updated(0);

			for (auto i = avoid.begin(); i != avoid.end(); ++i) { // move target_pos to avoid room objects
				cube_t c(*i);
				c.expand_by_xy(coll_dist);
				any_updated |= sphere_cube_int_update_pos(person.target_pos, 1.01*coll_dist, c, person.pos, 1, 1); // check_int=1, skip_z=1, ignore return value
			}
			if (!any_updated) break; // done
		} // for n
	}
	return 1;
}

bool building_t::select_person_dest_in_room(person_t &person, rand_gen_t &rgen, room_t const &room) const {
	float const height(0.7*get_window_vspace()), radius(COLL_RADIUS_SCALE*person.radius);
	point dest_pos(room.get_cube_center());
	static vect_cube_t avoid; // reuse across frames/people
	get_avoid_cubes(person.target_pos.z, height, radius, avoid, 0); // following_player=0
	bool const no_use_init(room.get_room_type(0) == RTYPE_PARKING); // don't use the room center for a parking garage
	if (!interior->nav_graph->find_valid_pt_in_room(avoid, radius, height, person.target_pos.z, room, rgen, dest_pos, no_use_init)) return 0;
	person.target_pos.x = dest_pos.x;
	person.target_pos.y = dest_pos.y;
	person.goal_type    = GOAL_TYPE_ROOM;
	return 1;
}

point get_pos_to_stand_for_elevator(person_t const &person, elevator_t const &e, float floor_spacing) {
	point pos(person.pos);
	pos[!e.dim] = e.get_center_dim(!e.dim); // centered on the elevator
	pos[ e.dim] = e.d[e.dim][e.dir] + (e.dir ? 1.0 : -1.0)*0.4*floor_spacing; // stand in front of the elevator door
	return pos;
}

int building_t::choose_dest_room(person_t &person, rand_gen_t &rgen) const { // used for randomly walking around the building

	assert(interior && interior->nav_graph);
	building_loc_t const loc(get_building_loc_for_pt(person.pos));
	if (loc.room_ix < 0) return 0; // not in a room
	person.dest_room    = person.cur_room = loc.room_ix; // set to the same room and pos just in case
	person.target_pos   = person.pos; // set but not yet used
	person.is_on_stairs = 0;
	if (interior->rooms.size() == 1) return 0; // no other room to move to
	float const floor_spacing(get_window_vspace());

	// maybe choose a destination elevator (for office buildings) if our prev dest wasn't an elevator
	if (ENABLE_AI_ELEVATORS && !is_house && !interior->elevators_disabled && !person.last_used_elevator && person.goal_type != GOAL_TYPE_ELEVATOR /*&& rgen.rand_float() < 0.25*/) {
		int const nearest_elevator(find_nearest_elevator_this_floor(person.pos));

		if (nearest_elevator >= 0) {
			elevator_t const &e(get_elevator(nearest_elevator));
			int const elevator_room(get_room_containing_pt(point(e.xc(), e.yc(), person.pos.z))); // room containing elevator center at person zval
			assert(elevator_room >= 0); // elevator must be in a valid room

			// Note: to simplify multiple AI logic, we could check e.in_use here and not even go to the elevator if someone else is using it
			if (interior->nav_graph->is_room_connected_to(loc.room_ix, elevator_room, interior->doors, person.pos.z, person.has_key)) { // check if reachable
				person.target_pos   = get_pos_to_stand_for_elevator(person, e, floor_spacing);
				person.dest_room    = elevator_room;
				person.goal_type    = GOAL_TYPE_ELEVATOR;
				person.cur_elevator = (uint8_t)nearest_elevator;
				person.last_used_elevator = 1;
				return 1;
			}
		}
	}

	// make 100 attempts at finding a valid room
	for (unsigned n = 0; n < 100; ++n) {
		unsigned const cand_room(rgen.rand() % interior->rooms.size());
		if (cand_room == (unsigned)loc.room_ix) continue;
		room_t const &room(interior->rooms[cand_room]);
		if (room.is_hallway) continue; // don't select a hallway
		if ((person.pos.z + floor_spacing) < room.z1() || (person.pos.z - floor_spacing) > room.z2()) continue; // room more than one floor above/below current pos
		// allow move to a different stacked part 25% of the time; 100% of the time for parking garages, since they're more rare
		if ((person.pos.z < room.z1() || person.pos.z > room.z2()) && !(has_parking_garage && room.z1() < ground_floor_z1) && (rgen.rand()&3) != 0) continue;
		if (!interior->nav_graph->is_room_connected_to(loc.room_ix, cand_room, interior->doors, person.pos.z, person.has_key)) continue;
		person.dest_room    = cand_room; // set but not yet used
		person.target_pos   = get_center_of_room(cand_room);
		person.target_pos.z = person.pos.z; // keep orig zval to stay on the same floor

		if (person.pos.z > room.z2()) { // room is below the person
			float const new_z(person.target_pos.z - floor_spacing);
			assert(new_z > room.z1() && new_z < room.z2());
			person.target_pos.z = new_z; // target the floor below
		}
		else if (person.pos.z < room.z1()) { // room is above the person
			float const new_z(person.target_pos.z + floor_spacing);
			assert(new_z > room.z1() && new_z < room.z2());
			person.target_pos.z = new_z; // target the floor above
		}
		else if (interior->ext_basement_hallway_room_id < 0 || cand_room < (unsigned)interior->ext_basement_hallway_room_id) { // skip for ext basement rooms
			// room covers floor this person is on; allow moving to a different floor, currently only one floor at a time
			cube_t const &part(get_part_for_room(room)); // or just use the room?
			unsigned const rand_val(rgen.rand() & 3); // 0-3

			if (rand_val == 0) { // try one floor below
				float const new_z(person.target_pos.z - floor_spacing);
				if (new_z > part.z1()) {person.target_pos.z = new_z;} // change if there is a floor below
			}
			else if (rand_val == 1) { // try one floor above
				float const new_z(person.target_pos.z + floor_spacing);
				if (new_z < part.z2()) {person.target_pos.z = new_z;} // change if there is a floor above
			}
			assert(person.target_pos.z > room.z1() && person.target_pos.z < room.z2());
		}
		person.goal_type = GOAL_TYPE_ROOM;
		return 1;
	} // for n
	if (loc.room_ix < 0) return 2; // failed - room not vali (error?)
	room_t const &room(get_room(loc.room_ix));
	bool const is_parking_garage(room.get_room_type(0) == RTYPE_PARKING);

	// how about a different floor of the same room? only check this 50% of the time for parking garages to allow movement within a level
	if (room.has_stairs == 255 && (!is_parking_garage || rgen.rand_bool())) {
		float const new_z(person.target_pos.z + (rgen.rand_bool() ? -1.0 : 1.0)*floor_spacing); // one floor above or below

		if (new_z > room.z1() && new_z < room.z2()) { // valid if this floor is inside the room
			person.target_pos.z = new_z;
			if (select_person_dest_in_room(person, rgen, room)) return 1;
		}
	}
	// how about a different location in the same room? this will at least get the person unstuck from an object and moving inside a parking garage
	if (person.is_first_path || is_parking_garage) {
		if (select_person_dest_in_room(person, rgen, room)) return 1;
	}
	return 2; // failed, but can retry
}

template<typename T> void add_bcube_if_overlaps_zval(vector<T> const &cubes, vect_cube_t &out, float z1, float z2) {
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {
		if (i->z1() < z2 && i->z2() > z1) {out.push_back(*i);}
	}
}

// for AI collision detection
void building_interior_t::get_avoid_cubes(vect_cube_t &avoid, float z1, float z2, float r_shrink_if_low, float floor_thickness, bool same_as_player, bool skip_stairs) const {
	avoid.clear();
	if (!skip_stairs) {add_bcube_if_overlaps_zval(stairwells, avoid, z1, z2);} // clearance not required
	add_bcube_if_overlaps_zval(elevators, avoid, z1, z2); // clearance not required
	if (!room_geom) return; // no room objects
	auto objs_end(room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators
	float const z_thresh(z1 + 0.35*(z2 - z1)); // used with r_shrink_if_low

	for (auto c = room_geom->objs.begin(); c != objs_end; ++c) {
		// these object types are not collided with by people and can be skipped
		// what about TYPE_CURB? should people step on/over these, or generally avoid them?
		if (c->no_coll() || c->is_dynamic() || c->type == TYPE_LG_BALL) continue; // skip dynamic objects (balls, etc.)
		if (!(same_as_player ? c->is_player_collidable() : bldg_obj_types[c->type].ai_coll)) continue;
		// skip_stairs also skips ramps? no, this seems to result in zombies walking in midair as if the ramp was a solid floor when approaching from above;
		// maybe zombies should walk down the ramp instead? but ramps are wider than stairs, if the player stands to the side then the zombie may walk right by;
		// so that means we need navigation to the side while on a ramp? this seems quite difficult for the current system to support
		//if (skip_stairs && c->type == TYPE_RAMP) continue;
		//if (c->type == TYPE_ATTIC_DOOR && (c->flags & RO_FLAG_IN_HALLWAY)) continue; // skip open attic doors in hallways because they block the path too much
		cube_t bc(get_true_room_obj_bcube(*c)); // needed for open attic door
		if (bc.z1() > z2 || bc.z2() < z1) continue;

		if (r_shrink_if_low > 0.0 && c->z2() < z_thresh && c->shape == SHAPE_CUBE) { // shrink cube if it's low; applies to boxes and crates on the floor
			bc.expand_by_xy(-min(0.95f*0.5f*min(c->dx(), c->dy()), r_shrink_if_low)); // make sure it doesn't shrink to zero area
		}
		avoid.push_back(bc);
		
		if (same_as_player && c->type == TYPE_TABLE && c->shape == SHAPE_CYLIN) {
			// special handling of round table so that the player can't hide in the area unreachable from the bounding cube
			float const shrink_val(c->get_radius()/SQRT2); // small shrink
			avoid.back().expand_by(-shrink_val);

			for (unsigned d = 0; d < 2; ++d) { // create a cross shape
				avoid.push_back(bc);
				avoid.back().expand_in_dim(d, -2.0*shrink_val);
			}
		}
	} // for c
}

bool building_t::stairs_contained_in_part(stairwell_t const &s, cube_t const &p) const {
	cube_t sc(s);
	if (s.roof_access) {sc.z2() -= get_window_vspace();} // clip off top floor roof access
	return p.contains_cube(sc);
}
void building_t::find_nearest_stairs_or_ramp(point const &p1, point const &p2, vector<unsigned> &nearest_stairs, int part_ix) const {
	nearest_stairs.clear();
	assert(interior);
	assert(part_ix < 0 || (unsigned)part_ix < parts.size());
	float const zmin(min(p1.z, p2.z)), zmax(max(p1.z, p2.z));
	vector<pair<float, unsigned>> sorted;

	for (unsigned s = 0; s < interior->stairwells.size(); ++s) {
		stairwell_t const &stairs(interior->stairwells[s]);
		if (zmin < stairs.z1() || zmax > stairs.z2()) continue; // stairs don't span the correct floors
		if (part_ix >= 0 && !stairs_contained_in_part(stairs, parts[part_ix])) continue; // stairs don't belong to this part (Note: this option is currently unused)
		point const center(stairs.get_cube_center());
		sorted.emplace_back((p2p_dist(p1, center) + p2p_dist(center, p2)), s);
	} // for s
	if ((part_ix < 0 || part_ix == basement_part_ix) && has_pg_ramp() && zmax < ground_floor_z1) { // parking garage ramp
		point const center(interior->pg_ramp.get_cube_center());
		sorted.emplace_back((p2p_dist(p1, center) + p2p_dist(center, p2)), interior->stairwells.size());
	}
	sort(sorted.begin(), sorted.end()); // sort by distance, min first
	for (auto s = sorted.begin(); s != sorted.end(); ++s) {nearest_stairs.push_back(s->second);}
}
int building_t::find_nearest_elevator_this_floor(point const &pos) const {
	int nearest(-1); // -1 => none found
	float dmin_sq(0.0);

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		if (e->z1() > pos.z || e->z2() < pos.z) continue; // doesn't span the correct floor
		// skip if elevator is currently in use?
		float const dsq(p2p_dist_xy_sq(pos, e->get_cube_center()));
		if (nearest < 0 || dsq < dmin_sq) {nearest = (e - interior->elevators.begin()); dmin_sq = dsq;}
	}
	return nearest;
}

cube_t get_stairs_plus_step_up(stairwell_t const &stairs) {
	cube_t stairs_ext(stairs);
	stairs_ext.d[stairs.dim][!stairs.dir] += (stairs.dir ? -1.0 : 1.0)*stairs.get_sz_dim(stairs.dim)/NUM_STAIRS_PER_FLOOR; // location of step up
	return stairs_ext;
}

void building_t::get_avoid_cubes(float zval, float height, float radius, vect_cube_t &avoid, bool following_player) const {
	assert(interior);
	interior->get_avoid_cubes(avoid, (zval - radius), (zval + (height - radius)), 0.5*radius, get_floor_thickness(), following_player);
}
bool building_t::find_route_to_point(person_t const &person, float radius, bool is_first_path, bool following_player, vector<point> &path) const {

	assert(interior && interior->nav_graph);
	point const &from(person.pos), &to(person.target_pos);
	path.clear();
	building_loc_t const loc1(get_building_loc_for_pt(from)), loc2(get_building_loc_for_pt(to));
	if (loc1.part_ix < 0 || loc2.part_ix < 0 || loc1.room_ix < 0 || loc2.room_ix < 0) return 0; // not in a room
	assert((unsigned)loc1.part_ix < parts.size() && (unsigned)loc2.part_ix < parts.size());
	assert((unsigned)loc1.room_ix < interior->rooms.size() && (unsigned)loc2.room_ix < interior->rooms.size());
	float const floor_spacing(get_window_vspace()), height(0.7*floor_spacing), z2_add(height - radius); // approximate, since we're not tracking actual heights
	static vect_cube_t avoid; // reuse across frames/people
	get_avoid_cubes(from.z, height, radius, avoid, following_player);

	if (loc1.same_room_floor(loc2)) { // same room/floor (not checking stairs_ix)
		assert(from.z == to.z);
		return interior->nav_graph->complete_path_within_room(from, to, loc1.room_ix, person.ssn, radius, person.cur_rseed, is_first_path, following_player, avoid, path);
	}
	if (loc1.floor_ix != loc2.floor_ix) { // different floors: find path from <from> to nearest stairs, then find path from stairs to <to>
		vector<unsigned> nearest_stairs_or_ramp;
		find_nearest_stairs_or_ramp(from, to, nearest_stairs_or_ramp); // pass in loc1.part_ix if both loc part_ix values are equal?
		bool const up_or_down(loc1.floor_ix > loc2.floor_ix); // 0=up, 1=down; Note: floor_ix is relative to bcube.z1() (or ext_basement z1), so is consistent across parts

		for (unsigned s : nearest_stairs_or_ramp) { // try using stairs or ramp, closest to furthest
			bool const is_ramp(s == interior->stairwells.size());
			assert(is_ramp || s < interior->stairwells.size());
			unsigned const stairs_room_ix(s + interior->rooms.size()); // map to graph space (should work for stairs or ramp)
			path.clear();
			vector<point> from_path;
			// Note: passing use_stairs=0 here because it's unclear if we want to go through stairs nodes in our A* algorithm
			// from => stairs/ramp
			if (!interior->nav_graph->find_path_points(loc1.room_ix, stairs_room_ix, person.ssn, radius, height, 0, is_first_path,
				up_or_down, person.cur_rseed, avoid, from, interior->doors, person.has_key, nullptr, from_path)) continue; // no custom_dest
			point const seg2_start(interior->nav_graph->get_stairs_entrance_pt(to.z, stairs_room_ix, !up_or_down)); // other end
			// new floor, new zval, new avoid cubes
			interior->get_avoid_cubes(avoid, (seg2_start.z - radius), (seg2_start.z + z2_add), 0.5*radius, get_floor_thickness(), following_player);
			// stairs/ramp => to
			if (!interior->nav_graph->find_path_points(stairs_room_ix, loc2.room_ix, person.ssn, radius, height, 0, is_first_path,
				!up_or_down, person.cur_rseed, avoid, seg2_start, interior->doors, person.has_key, nullptr, path)) continue; // no custom_dest
			assert(!path.empty() && !from_path.empty());
			path.push_back(seg2_start); // other end of the stairs
			// add two more points to straighten the entrance and exit paths; this segment doesn't check for intersection with stairs
			point enter_pt;

			if (is_ramp) {
				cube_t const walk_area(interior->pg_ramp);
				enter_pt = walk_area.closest_pt(from_path.front());
				path.push_back(walk_area.closest_pt(path.back())); // exit point
			}
			else { // stairs
				stairwell_t const &stairs(interior->stairwells[s]);
				cube_t const stairs_ext(get_stairs_plus_step_up(stairs));
				enter_pt = stairs_ext.closest_pt(from_path.front());
				point const exit_pt(stairs_ext.closest_pt(path.back()));
				path.push_back(exit_pt);

				if (stairs.shape == SHAPE_U) { // add 2 extra points on mid-level landing; entrance and exit will be on the same side
					bool const dim(stairs.dim), dir(stairs.dir); // Note: see code in add_stairs_and_elevators()
					float const turn_pt(stairs.d[dim][dir] - 0.1*(dir ? 1.0 : -1.0)*stairs.get_sz_dim(dim)), seg_delta_z(0.45f*(to.z - from.z));
					point exit_turn(exit_pt.x, exit_pt.y, (to.z - seg_delta_z));
					exit_turn[dim] = turn_pt;
					path.push_back(exit_turn); // turning point for exit side of stairs
					point enter_turn(enter_pt.x, enter_pt.y, (from.z + seg_delta_z));
					enter_turn[dim] = turn_pt;
					path.push_back(enter_turn); // turning point for entrance side of stairs
				}
			}
			path.push_back(enter_pt);
			vector_add_to(from_path, path); // concatenate the two path segments in reverse order
			assert(!path.empty());
			return 1; // done/success
		} // for s
		path.clear(); // not necessary?
		return 0; // failed
	}
	assert(loc1.room_ix != loc2.room_ix);
	// if the target is an elevator, use that as the preferred destination rather than the center of the room
	point const *const custom_dest((person.goal_type == GOAL_TYPE_ELEVATOR) ? &person.target_pos : nullptr);
	if (!interior->nav_graph->find_path_points(loc1.room_ix, loc2.room_ix, person.ssn, radius, height, 0, is_first_path,
		0, person.cur_rseed, avoid, from, interior->doors, person.has_key, custom_dest, path)) return 0;
	assert(!path.empty());
	return 1;
}

void person_t::next_path_pt(bool starting_path) {
	assert(!path.empty());
	is_on_stairs = (!starting_path && target_pos.z != path.back().z); // or ramp
	target_pos   = path.back();
	path.pop_back();
}

bool building_t::is_valid_ai_placement(point const &pos, float radius, bool skip_nocoll) const { // for people and animals
	if (!is_pos_inside_building(pos, radius, radius)) return 0; // required for attic
	cube_t ai_bcube(pos);
	ai_bcube.expand_by(radius); // expand more in Z?
	if (!is_valid_stairs_elevator_placement(ai_bcube, radius)) return 0;

	// Note: people are placed before room geom is generated for all buildings, so this may not work and will have to be handled during room geom placement
	if (interior->room_geom) { // check placement against room geom objects
		auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) {
			if (skip_nocoll && i->no_coll()) continue;
			if (i->type == TYPE_FLOORING || i->type == TYPE_BLOCKER) continue; // okay to place on flooring; ignore blockers (used for placement clearance)
			if (i->intersects(ai_bcube)) return 0;
		}
	}
	return 1;
}

struct room_cand_t {
	unsigned room_ix, floor_ix;
	room_cand_t(unsigned r, unsigned f) : room_ix(r), floor_ix(f) {}
};

bool building_t::place_people_if_needed(unsigned building_ix, float radius, vector<point> &locs) const {
	if (!interior || interior->rooms.empty() || is_rotated()) return 0; // no people in these cases
	if (interior->placed_people) return 0; // already placed
	interior->placed_people = 1; // set, even if no people are placed below
	if (!global_building_params.gen_building_interiors) return 0; // no interiors, no people
	unsigned num_min(is_house ? global_building_params.people_per_house_min : global_building_params.people_per_office_min);
	unsigned num_max(is_house ? global_building_params.people_per_house_max : global_building_params.people_per_office_max);
	if (num_max < num_min) {swap(num_min, num_max);} // or error? if so, it should be checked during option processing
	if (num_max == 0) return 0;
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, mat_ix); // should be canonical per building
	unsigned const num_people(num_min + (rgen.rand()%(num_max - num_min + 1)));
	if (num_people == 0) return 0;
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	// weight rooms by number of floors to avoid placing too many people in basements
	static vector<room_cand_t> room_cands;
	room_cands.clear();
	unsigned first_basement_room(0);

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) { // add room_cands
		if (r->is_sec_bldg) continue; // don't place people in garages and sheds
		if (min(r->dx(), r->dy()) < 3.0*radius) continue; // room to small to place a person
		unsigned const num_floors(calc_num_floors(*r, window_vspacing, floor_thickness));
		assert(num_floors > 0);
		if (first_basement_room == 0 && r->z1() < ground_floor_z1) {first_basement_room = room_cands.size();}
		unsigned const room_ix(r - interior->rooms.begin());

		for (unsigned f = 0; f < num_floors; ++f) {
			if (r->lit_by_floor && !r->is_lit_on_floor(f)) continue; // don't place person in an unlit room; only applies if room lighting has been calculated
			room_cands.emplace_back(room_ix, f);
		}
	} // for r
	for (unsigned N = 0; N < num_people && !room_cands.empty(); ++N) {
		unsigned const max_cand_ix((N == 0 && first_basement_room > 0) ? first_basement_room : room_cands.size()); // place the first person in a non-basement room

		for (unsigned n = 0; n < 10; ++n) { // make 10 attempts at choosing a room
			unsigned const cand_ix(rgen.rand() % max_cand_ix);
			room_cand_t const &cand(room_cands[cand_ix]);
			room_t const &room(get_room(cand.room_ix));
			float const zval(room.z1() + fc_thick + window_vspacing*cand.floor_ix);
			bool success(0);

			for (unsigned m = 0; m < 100; ++m) { // make 100 attempts at choosing a position in this room
				point const pos(gen_xy_pos_in_area(room, radius, rgen, zval)); // random XY point inside this room
				if (!is_valid_ai_placement(pos, radius, 1)) continue; // skip_nocoll=1
				locs.push_back(pos);
				success = 1;
				break; // done/success
			} // for n
			if (success) {
				// don't place another person on this floor of this room unless there aren't enough rooms for everyone
				if (room_cands.size() >= (num_people - N)) {room_cands.erase(room_cands.begin() + cand_ix);}
				break;
			}
		} // for n
	} // for N
	return 1;
}

bool can_ai_follow_player(person_t const &person) {
	if (!ai_follow_player()) return 0; // disabled
	if (!cur_player_building_loc.is_valid()) return 0; // no target
	if (cur_player_building_loc.building_ix != person.cur_bldg) return 0; // wrong building
	if (player_is_hiding) return 0; // ignore player if in the closet, bathroom stall, or shower with the door closed
	if (person.retreat_time > 0.0) return 0; // ignore the player if retreating
	return 1;
}
bool has_nearby_sound(person_t const &person, float floor_spacing) {
	if (!can_ai_follow_player(person)) return 0; // no need to track sounds
	point sound_pos; // unused
	return get_closest_building_sound(person.pos, sound_pos, floor_spacing);
}

bool building_t::same_room_and_floor_as_player(person_t const &person) const {
	return (cur_player_building_loc.room_ix == person.cur_room && cur_player_building_loc.floor_ix == get_floor_for_zval(person.pos.z) && cur_player_building_loc.stairs_ix < 0);
}
bool building_t::is_player_visible(person_t const &person, unsigned vis_test) const {
	if (vis_test == 0) return 1; // no visibility test
	building_dest_t const &target(cur_player_building_loc);
	float const player_radius(get_scaled_player_radius());
	point const pp2(target.pos - vector3d(0.0, 0.0, get_player_height())); // player's bottom sphere
	bool const same_room(person.cur_room >= 0 && cur_player_building_loc.room_ix == person.cur_room);
	unsigned const person_floor_ix(get_floor_for_zval(person.pos.z));
	unsigned const floor_delta(abs((int)cur_player_building_loc.floor_ix - (int)person_floor_ix));
	bool const same_room_and_floor(same_room && floor_delta == 0); // Note: doesn't check cur_player_building_loc.stairs_ix
	bool has_los(same_room_and_floor);

	if (!has_los && same_room && floor_delta == 1 && person.pos.z > ground_floor_z1) {
		// if the person and the player are on adjacent floors of the same room connected by stairs (and not in a parking garage), cheat and say they have a line of sight
		room_t const &room(get_room(person.cur_room));
		unsigned const room_floor_start(get_floor_for_zval(room.z1()));
		assert(room_floor_start <= cur_player_building_loc.floor_ix && room_floor_start <= person_floor_ix);
		unsigned const f1(cur_player_building_loc.floor_ix - room_floor_start), f2(person_floor_ix - room_floor_start);
		has_los = (room.has_stairs_on_floor(f1) && room.has_stairs_on_floor(f2));
	}
	if (!has_los) { // check visibility; assume LOS if in the same room
		point const eye_pos(person.get_eye_pos());
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
bool building_t::can_target_player(person_t const &person) const {
	if (!can_ai_follow_player(person)) return 0;
	return is_player_visible(person, global_building_params.ai_player_vis_test); // 0=no test, 1=LOS, 2=LOS+FOV, 3=LOS+FOV+lit
}

bool building_t::need_to_update_ai_path(person_t const &person) const {
	if (!global_building_params.ai_target_player || !can_ai_follow_player(person) || !interior) return 0; // disabled
	building_dest_t const &target(cur_player_building_loc);
	bool const same_room(same_room_and_floor_as_player(person)); // check room and floor

	if (global_building_params.ai_player_vis_test == 0) { // the below logic only applies if there's no line of sight test and is an optimization
		// if the player's room has not changed, and the person is not yet in this room, continue on the same path to the dest room;
		// however, if the path is empty, continue to choose a new path (needed for AI more than one floor away from target to take multiple flights of stairs)
		if (target.same_room_floor(prev_player_building_loc) && !same_room && !person.on_new_path_seg && !person.path.empty()) return 0;
	}
	//if (same_room && person.path.size() > 1) return 0; // same room but path has a jog, continue on existing path (faster, but slower to adapt to player position change)
	if (dist_less_than(person.pos, target.pos, person.radius)) return 0; // already close enough
	if (person.on_stairs()) return 0; // don't change paths when on the stairs
	float const floor_spacing(get_window_vspace());

	if (int(target.pos.z/floor_spacing) != int(prev_player_building_loc.pos.z/floor_spacing)) { // if player did not change floors
		if (fabs(person.pos.z - target.pos.z) > 2.0f*floor_spacing) return 0; // person and player are > 2 floors apart, continue toward stairs (or should it be one floor apart?) (optimization)
	}
	if (can_target_player(person)) { // have player visibility
		if (person.goal_type == GOAL_TYPE_PLAYER && target.same_room_floor(prev_player_building_loc) && !same_room && !person.on_new_path_seg && !person.path.empty()) return 0;
		return 1;
	}
	if (has_nearby_sound(person, floor_spacing)) return 1; // new sound source
	return 0; // continue on the previously chosen path
}

void building_t::all_ai_room_update(rand_gen_t &rgen, float delta_dir) {
	assert(interior);

	for (unsigned i = 0; i < interior->people.size(); ) { // Note: no increment
		person_t &person(interior->people[i]);
		person.ai_state = ai_room_update(person, delta_dir, i, rgen);
		if (person.ai_state == AI_TO_REMOVE) {interior->people.erase(interior->people.begin() + i);} // remove this person
		else {++i;}
	}
}

unsigned get_elevator_floor(float zval, elevator_t const &e, float floor_spacing) { // floor index relative to this elevator
	return max(0.0f, (zval - e.z1()))/floor_spacing;
}
float get_person_max_move_dist(person_t const &person, float speed_mult=1.0) {
	return person.speed*speed_mult*min(fticks, 4.0f); // clamp fticks to 100ms
}
float move_person_forward_to_target(person_t &person) { // Note: expected that person.target_pos.z == person.pos.z
	float const move_dist(get_person_max_move_dist(person));
	vector3d const delta(person.target_pos - person.pos);
	person.pos       += (move_dist/delta.mag())*delta;
	person.anim_time += move_dist;
	return move_dist;
}
bool person_slow_turn(person_t &person, point const &target, float delta_dir) {
	vector3d dir_to_target((target.x - person.pos.x), (target.y - person.pos.y), 0.0); // zval is always 0
	dir_to_target.normalize();
	if (dot_product(person.dir, dir_to_target) > 0.99) return 0; // turn is complete
	// if directions are nearly opposite, pick a side to turn using the cross product to get an orthogonal vector
	if (dot_product(dir_to_target, person.dir) < -0.9) {dir_to_target = cross_product(person.dir, plus_z);}
	person.dir = delta_dir*dir_to_target + (1.0 - delta_dir)*person.dir; // merge new_dir into dir gradually for smooth turning
	person.dir.normalize();
	return 1;
}

void building_t::call_elevator_to_floor_and_light_nearest_button(elevator_t &elevator, unsigned floor_ix, bool is_inside_elevator, bool is_up) {
	call_elevator_to_floor(elevator, floor_ix, is_inside_elevator, is_up);
	// light the button that was pressed
	if (!has_room_geom()) return; // can't light the button
	assert(elevator.button_id_start <= elevator.button_id_end && elevator.button_id_end <= interior->room_geom->objs.size());

	for (auto i = interior->room_geom->objs.begin() + elevator.button_id_start; i != interior->room_geom->objs.begin() + elevator.button_id_end; ++i) {
		if (i->type == TYPE_BLOCKER) continue; // button was removed?
		assert(i->type == TYPE_BUTTON);
		if (i->obj_id != floor_ix) continue; // wrong floor
		if (i->in_elevator() != is_inside_elevator) continue; // wrong inside/outside
		if (!is_inside_elevator && (i->flags & (is_up ? RO_FLAG_ADJ_BOT : RO_FLAG_ADJ_TOP))) continue; // wrong up/down button
		i->flags |= RO_FLAG_IS_ACTIVE; // set active/lit state
		interior->room_geom->invalidate_small_geom(); // need to regen object data due to lit state change; should be thread safe
		break; // only one button
	} // for i
}

int building_t::run_ai_elevator_logic(person_t &person, float delta_dir, rand_gen_t &rgen) {
	assert(person.goal_type == GOAL_TYPE_ELEVATOR);
	if (!has_room_geom()) return person.ai_state; // if player has moved away and room_geom was deleted, remain in the current state
	float const floor_spacing(get_window_vspace());
	elevator_t &e(get_elevator(person.cur_elevator));
	room_object_t const &ecar(interior->get_elevator_car(e));

	// Note: when AI_LOCKS_ELEVATOR==1, only one person can be in the elevator at once, including walking into and out of the elevator
	// Note: if interior->elevators_disabled==1, the AI will either give up waiting or get stuck in the elevator if already inside
	if (person.ai_state == AI_WAIT_ELEVATOR) {
		// waiting for the elevator to arrive; we allow the AI to be pushed and to give up if waiting too long;
		// when the elevator opens, we could check e.going_up to see if it's headed in the correct direction;
		// however, since we haven't decided on a dest floor yet, it makes sense to just go with it and select a floor in the direction the elevator is headed
		point const elevator_center(e.xc(), e.yc(), person.pos.z);

		if (e.open_amt == 1.0 && e.num_occupants < ELEVATOR_CAPACITY) { // doors are fully open, we can fit
			if (get_elevator_floor(ecar.zc(), e, floor_spacing) == get_elevator_floor(person.pos.z, e, floor_spacing)) { // wait for elevator to reach our current floor
				person.dir = (elevator_center - person.pos).get_norm(); // snap our direction to forward, in the rare case the elevator arrives before we've completed our turn
				person.waiting_start = 0.0; // no longer waiting for elevator
				++e.num_occupants; // make space for ourselves in the elevator
				return AI_ENTER_ELEVATOR;
			}
		}
		person_slow_turn(person, elevator_center, 0.5*delta_dir); // slow turn to face the elevator
	}
	else if (person.ai_state == AI_ENTER_ELEVATOR) {
		person.target_pos.assign(e.xc(), e.yc(), person.pos.z);
		float const move_dist(move_person_forward_to_target(person)); // walk into elevator

		if (dist_xy_less_than(person.pos, person.target_pos, move_dist)) {
			person.anim_time = 0.0; // stop and turn
			return AI_ACTIVATE_ELEVATOR;
		}
	}
	else if (person.ai_state == AI_ACTIVATE_ELEVATOR) {
		vector3d const prev_dir(person.dir);
		bool const is_turning(person_slow_turn(person, get_pos_to_stand_for_elevator(person, e, floor_spacing), 0.5*delta_dir)); // slow turn to face the elevator door

		if (!is_turning) { // turn completed
			call_elevator_to_floor_and_light_nearest_button(e, person.dest_elevator_floor, 1, 0); // is_inside_elevator=1, is_up=0
			person.anim_time = 0.0; // stop and wait
			return AI_RIDE_ELEVATOR;
		}
	}
	else if (person.ai_state == AI_RIDE_ELEVATOR) {
		person.pos.z = ecar.z1() + person.radius + get_fc_thickness(); // move with the elevator

		if (e.open_amt == 1.0) { // doors are fully open
			// Note: we could probably query e.get_target_floor() for this
			if (get_elevator_floor(ecar.zc(), e, floor_spacing) == person.dest_elevator_floor) { // at destination floor
				return AI_EXIT_ELEVATOR;
			}
		}
	}
	else if (person.ai_state == AI_EXIT_ELEVATOR) {
		person.target_pos = get_pos_to_stand_for_elevator(person, e, floor_spacing);
		float const move_dist(move_person_forward_to_target(person)); // exit elevator to a point in front, then select a new non-elevator destination

		if (dist_xy_less_than(person.pos, person.target_pos, move_dist)) {
			person.target_pos = all_zeros;
			//assert(e.num_occupants > 0); // invalid if people are despawned/respawned when the player is far away?
			if (e.num_occupants > 0) {--e.num_occupants;} // we're no longer in this elevator
			return AI_AT_DEST;
		}
	}
	else {assert(0);} // invalid state
	return person.ai_state; // fallthrough case, no change to ai_state
}

// Note: non-const because this updates room light and door state
int building_t::ai_room_update(person_t &person, float delta_dir, unsigned person_ix, rand_gen_t &rgen) {

	if (person.speed == 0.0) {person.anim_time = 0.0; return AI_STOP;} // stopped
	assert(interior);
	if (!interior->room_geom && frame_counter < 60) {person.anim_time = 0.0; return AI_WAITING;} // wait until room geom is generated for this building
	float const coll_dist(COLL_RADIUS_SCALE*person.radius);
	float &wait_time(person.waiting_start); // reuse this field
	float speed_mult(1.0);
	person.following_player = 0; // reset for this frame

	if (person.retreat_time > 0.0) {
		if (person.retreat_time == RETREAT_TIME) { // first retreating frame - clear path
			person.path.clear();
			person.target_pos = all_zeros;
			//person.is_first_path = 1; // probably not needed
		}
		wait_time = 0.0; // no waiting while retreating
		person.retreat_time -= fticks;
		max_eq(person.retreat_time, 0.0f);
	}
	if (wait_time > 0) {
		if (wait_time > fticks && !can_ai_follow_player(person)) { // waiting; don't wait if there's a player to follow
			// check for other people colliding with this person and handle it; should this apply when waiting for the elevator?
			for (auto p = interior->people.begin()+person_ix+1; p < interior->people.end(); ++p) {
				if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
				float const rsum(coll_dist + COLL_RADIUS_SCALE*p->radius);
				if (!dist_xy_less_than(person.pos, p->pos, rsum)) continue; // not intersecting
				move_person_to_not_collide(person, *p, person.pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
			} // for p
			wait_time -= fticks;
			person.anim_time = 0.0; // reset just in case (though should already be at 0.0)

			if (person.ai_state != AI_WAIT_ELEVATOR) { // don't reset goal and return here if waiting at an elevator
				person.goal_type = GOAL_TYPE_NONE; // our orig goal is invalid
				return AI_WAITING;
			}
		}
		else { // wait is up
			wait_time = 0.0;
			person.target_pos = all_zeros; // force choose_dest=1 below
			person.ai_state   = AI_WAITING; // waiting, but no longer for an elevator
		}
	}
	if (!point_in_building_or_basement_bcube(person.pos)) { // person must be inside the building
		cout << TXT(person.pos.str()) << TXT(bcube.str()) << TXT(interior->basement_ext_bcube.str()) << endl;
		assert(0);
	}
	if (person.ai_state >= AI_WAIT_ELEVATOR) {return run_ai_elevator_logic(person, delta_dir, rgen);} // handle elevator case
	build_nav_graph();

	if (can_ai_follow_player(person)) {
		// use zval of the feet to handle cases where the person and the player are different heights
		point const feet_pos(person.pos.x, person.pos.y, person.get_z1()), player_feet_pos(cur_player_building_loc.pos - vector3d(0.0, 0.0, CAMERA_RADIUS+get_player_height()));

		if (dist_less_than(feet_pos, player_feet_pos, 1.2f*(person.radius + get_scaled_player_radius()))) {
			if (!check_for_wall_ceil_floor_int(person.pos, cur_player_building_loc.pos)) {
				int const ret_status(register_ai_player_coll(person.has_key, person.get_height())); // return value: 0=no effect, 1=player is killed, 2=this person is killed
				// player is killed, we could track kills here
				if (ret_status == 1) {add_blood_decal(cur_player_building_loc.pos, get_scaled_player_radius());}
				else if (ret_status == 2) {return AI_TO_REMOVE;} // player defeats zombie, remove it
			}
		}
	}
	bool choose_dest(!person.target_valid());
	bool const update_path(need_to_update_ai_path(person));

	if (update_path) { // need to update based on player movement; higher priority than choose_dest
		if (choose_dest_goal(person, rgen) != 1 || // check if person can reach the target
			!find_route_to_point(person, coll_dist, 0, 1, person.path)) // is_first_path=0, following_player=1
		{
			choose_dest = 1; // or increment person.cur_rseed and return AI_WAITING? or restore person to prev value?
		}
		else { // success
			person.next_path_pt(1);
			person.following_player = 1;
			choose_dest = 0;
			speed_mult  = 1.5; // faster when the player is in the same room
			// run logic to play zombie sounds
			bool const same_room_and_floor(same_room_and_floor_as_player(person));
			bool play_sound(same_room_and_floor); // always play sound if in the same room and floor
			if (!play_sound && (person_ix & 1)) {play_sound |= is_player_visible(person, 1);} // 50% of zombies use line of sight test
			if (!play_sound && (person_ix & 2)) {play_sound |= has_nearby_sound(person, get_window_vspace());} // 50% of zombies use sound test
			if (play_sound) {maybe_play_zombie_sound(person.pos, person_ix, same_room_and_floor);} // alert other zombies if in the same room and floor as the player
		}
	}
	if (choose_dest) { // no current destination, choose a new destination
		++person.cur_rseed; // update rseed so that this person will choose a different path if the previous path finding call failed
		if (!update_path) {person.anim_time = 0.0;} // reset animation if not a failed update path frame
		int const ret(choose_dest_room(person, rgen)); // 0=failed, 1=success, 2=failed but can retry

		if (ret == 2 && interior->door_state_updated) { // wait rather than stopping in case the player trapped this person in a room by closing the door
			person.wait_for(2.0); // stop for 2 seconds, then try again
			return AI_WAITING;
		}
		else if (ret != 1) { // if there's no valid room or valid path, set the speed to 0 so that we don't check this every frame
			if (!ai_follow_player()) { // if not following the player
				//person.speed = 0.0; // movement will be stopped from now on - not valid if we get into this state before gameplay is enabled
				person.wait_for(2.0); // stop for 2 seconds, then try again
			}
			return AI_STOP;
		}
		if (!find_route_to_point(person, coll_dist, person.is_first_path, 0, person.path)) { // following_player=0
			person.wait_for(1.0); // stop for 1 second, then try again
			return AI_WAITING;
		}
		if (has_room_geom()) {person.is_first_path = 0;} // treat the path as the first path until room geom is generated
		person.next_path_pt(1);
		return AI_BEGIN_PATH;
	}
	float const max_dist(get_person_max_move_dist(person, speed_mult));
	person.on_new_path_seg = 0; // clear flag for this frame
	//person.following_player = can_target_player(person); // for debugging visualization

	if (dist_less_than(person.pos, person.target_pos, 1.1f*max_dist)) { // at dest
		assert(point_in_building_or_basement_bcube(person.target_pos));
		person.pos = person.target_pos;
		if (!person.path.empty()) {person.next_path_pt(0); return AI_NEXT_PT;} // move to next path point
		float const floor_spacing(get_window_vspace());

		if (person.goal_type == GOAL_TYPE_ELEVATOR) {
			elevator_t &e(get_elevator(person.cur_elevator));
			// floor index relative to this elevator, not the room or building
			unsigned const cur_floor(get_elevator_floor(person.pos.z, e, floor_spacing)), num_floors(round_fp(e.dz()/floor_spacing));
			assert(num_floors > 1);
			assert(cur_floor  < num_floors);

			// select the destination floor, different from the current floor; this must be done before calling the elevator so that we know if we're going up or down
			while (1) {
				person.dest_elevator_floor = rgen.rand() % num_floors;
				if (person.dest_elevator_floor != cur_floor) break; // floor is valid
			}
			bool const is_up(person.dest_elevator_floor > cur_floor);
			call_elevator_to_floor_and_light_nearest_button(e, cur_floor, 0, is_up); // is_inside_elevator=0
			vector3d const dir_to_elevator(e.get_cube_center() - person.pos);
			// snap to face the elevator; would be better to gradually turn, but it's not clear how to do that; at least we should be facing in that general direction
			person.dir.assign(dir_to_elevator.x, dir_to_elevator.y, 0.0);
			person.dir.normalize();
			person.wait_for(60.0); // wait there for the elevator to arrive; if we're waiting for more than 60s, give up and choose another dest
			return AI_WAIT_ELEVATOR;
		}
		// don't wait if we can follow the player
		bool const no_wait(global_building_params.ai_target_player && (can_target_player(person) || has_nearby_sound(person, floor_spacing)));
		if (!no_wait) {person.wait_for(rgen.rand_uniform(1.0, (can_ai_follow_player(person) ? 2.0 : 10.0)));} // stop for 1-10s, 1-2s if player is in this building in gameplay mode
		person.on_new_path_seg    = 1; // allow player following AI update logic to rerun this frame
		person.last_used_elevator = 0;
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
	}
	else { // optimization for aligned dir
		new_pos = person.pos + max_dist*person.dir;
	}
	// make sure the person is inside the building, in case they were pushed by another person
	cube_t clip_cube;
		
	if (has_basement() && new_pos.z < ground_floor_z1) { // in the basement
		cube_t const &basement(get_basement());

		if (has_ext_basement()) {
			cube_t sc; sc.set_from_sphere(new_pos, coll_dist); // sphere bounding cube
			float cont_area(0.0);
			accumulate_shared_xy_area(basement, sc, cont_area);
			clip_cube = get_full_basement_bcube(); // start with full union bcube for basement

			if (basement.contains_pt(new_pos)) { // primarily in the basement
				accumulate_shared_xy_area(interior->basement_ext_bcube, sc, cont_area);
				if (cont_area < 0.99*sc.get_area_xy()) {clip_cube = basement;} // not contained - force into the basement
			}
			else { // primarily in the extended basement
				cube_t cur_room_bcube(basement); // start at the basement in case we're not even in a room, so we can at least snap here as a fallback

				for (auto r = interior->ext_basement_rooms_start(); r != interior->rooms.end(); ++r) {
					if (new_pos.z < r->z1() || new_pos.z > r->z2()) continue; // wrong level
					if (r->contains_pt(new_pos)) {cur_room_bcube = *r;} // this is the room we'll be pushed into if needed
					accumulate_shared_xy_area(*r, sc, cont_area);
				}
				if (cont_area < 0.99*sc.get_area_xy()) {clip_cube = cur_room_bcube;} // not contained - force into containing room
			}
		}
		else {clip_cube = basement;} // basement only
	}
	else {clip_cube = bcube;} // above ground
	clip_cube.expand_by_xy(-coll_dist); // shrink
	clip_cube.clamp_pt_xy(new_pos); // make sure person stays within building bcube; can't clip to room because person may be exiting it
	
	if (!point_in_building_or_basement_bcube(new_pos)) { // person must be inside the building
		cout << TXT(new_pos.str()) << TXT(bcube.str()) << TXT(interior->basement_ext_bcube.str()) << endl;
		assert(0);
	}
	// don't do collision detection while on stairs because it doesn't work properly; just let people walk through each other
	if (!person.on_stairs()) {
		for (auto p = interior->people.begin()+person_ix+1; p < interior->people.end(); ++p) { // check all other people in the same building after this one and attempt to avoid them
			if (fabs(person.pos.z - p->pos.z) > coll_dist) continue; // different floors
			float const rsum(coll_dist + COLL_RADIUS_SCALE*p->radius);
			if (!dist_xy_less_than(new_pos, p->pos, rsum)) continue; // new pos not close
			if (!dist_xy_less_than(person.pos, p->pos, rsum)) return AI_STOP; // old pos not intersecting, stop
			person.anim_time = 0.0; // pause animation in case this person is mid-step
			move_person_to_not_collide(person, *p, new_pos, rsum, coll_dist); // if we get here, we have to actively move out of the way
			return AI_MOVING; // return here, but don't update animation or dir; only handles a single collision
		} // for p
	}
	bool const player_in_this_building(cur_player_building_loc.building_ix == person.cur_bldg); // basement door only counts if the player is in this building
	bool const might_have_closed_door(global_building_params.open_door_prob < 1.0 || (player_in_this_building && is_house && has_basement()));

	if (interior->door_state_updated || (global_building_params.ai_opens_doors == 2 && might_have_closed_door)) {
		cube_t sc; sc.set_from_sphere(new_pos, person.radius); // sphere bounding cube

		for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) { // can be slow, but not as slow as iterating over doors
			if (new_pos.z < i->z1() || new_pos.z > i->z2())         continue; // wrong part/floor
			if (!i->intersects(sc)) continue; // no intersection with door
			if (!i->get_true_bcube().line_intersects(person.pos, person.target_pos)) continue; // check if path goes through door, to allow for "glancing blows" when pushed or turning
			assert(i->first_door_ix < interior->doors.size());

			for (unsigned dix = i->first_door_ix; dix < interior->doors.size(); ++dix) {
				door_t const &door(interior->doors[dix]);
				if (!i->is_same_stack(door)) break; // moved to a different stack, done
				if (door.z1() > person.pos.z || door.z2() < person.pos.z) continue; // wrong floor
				if (door.open) continue; // doors tend to block the AI, don't collide with them unless they're closed

				if (global_building_params.ai_opens_doors && !door.is_locked_or_blocked(person.has_key)) { // can open the door
					toggle_door_state(dix, player_in_this_building, 0, person.pos.z); // by_player=0
				}
				else { // can't open the door
					person.wait_for(5.0); // wait for 5s and then choose a new desination
					return AI_WAITING; // cut the path short at this closed door
				}
			} // for dix
		} // for i
	}
	handle_vert_cylin_tape_collision(new_pos, person.pos, person.get_z1(), person.get_z2(), person.radius, 0); // should be okay to use zvals from old pos; is_player=0
	// logic to clip this person to correct room Z-bounds in case something went wrong; remove if/when this is fixed
	// Note: we probably can't use the room Z bounds here becase the person may be on the stairs connecting two stacked parts
	float const true_z1(point_in_extended_basement_not_basement(new_pos) ? get_bcube_z1_inc_ext_basement() : bcube.z1());
	float const min_valid_zval(true_z1 + get_fc_thickness() + person.radius), max_valid_zval(bcube.z2() - person.radius); // Note: max should include the attic

	if (/*player_in_this_building &&*/ !person.on_stairs() && person.cur_room >= 0) { // movement in XY, not on stairs, room is valid: snap to nearest floor
		// this is optional and is done just in case something went wrong
		room_t room(get_room(person.cur_room));

		if (!room.contains_pt(new_pos)) { // maybe we just exited the stairs into a different part and cur_room hasn't been updated yet
			int const new_room_ix(get_room_containing_pt(new_pos));
			if (new_room_ix >= 0) {room = get_room(new_room_ix);}
		}
		float const floor_spacing(get_window_vspace());
		int cur_floor(max(0, round_fp((new_pos.z - min_valid_zval)/floor_spacing)));
		int const max_floor(round_fp((room.z2() - true_z1)/floor_spacing) - 1);
		min_eq(cur_floor, max_floor); // clip to the valid floors for this room relative to lowest building floor
		float const adj_zval(cur_floor*floor_spacing + min_valid_zval);
		if (fabs(adj_zval - new_pos.z) > 0.1*person.radius) {person.target_pos = all_zeros; person.path.clear();} // if we snap to the floor, reset the target and path
		new_pos.z = adj_zval;
	}
	max_eq(new_pos.z, min_valid_zval); // don't let the person go below the ground floor
	min_eq(new_pos.z, max_valid_zval); // don't let the person go above the room ceiling
	// update state
	person.pos        = new_pos; // Note: new_pos.z should equal person.poz.z unless on stairs, which is difficult to accurately check for in this function
	person.anim_time += max_dist;
	bool enable_bl_update(display_mode & 0x20); // disabled by default, enable with key '6'
	if (ai_follow_player() && global_building_params.ai_player_vis_test >= 3) {enable_bl_update = 0;} // if AI tests that player is lit, don't turn on/off lights
	if (enable_bl_update) {ai_room_lights_update(person);} // non-const part
	if (player_in_this_building) {person.cur_room = get_room_containing_pt(person.pos);} // update cur_room after moving and lights update
	return AI_MOVING;
}

void building_t::ai_room_lights_update(person_t const &person) {
	int const room_ix(get_room_containing_pt(person.pos));
	if (room_ix < 0) return; // room is not valid (between rooms, etc.)
	set_room_light_state_to(get_room(room_ix), person.pos.z, 1); // make sure current room light is on
	if (room_ix == person.cur_room) return; // same room as last time - done
	bool other_person_in_room(0);

	// check for other people in the room before turning the lights off on them
	for (person_t const &p : interior->people) {
		if (p.pos == person.pos) continue; // skip ourself
		if (get_room_containing_pt(p.pos) == person.cur_room && fabs(person.pos.z - p.pos.z) < get_window_vspace()) {other_person_in_room = 1; break;}
	}
	if (!other_person_in_room) {set_room_light_state_to(get_room(person.cur_room), person.pos.z, 0);} // make sure old room light is off
}

void building_t::move_person_to_not_collide(person_t &person, person_t const &other, point const &new_pos, float rsum, float coll_dist) const {
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

	if (!point_in_building_or_basement_bcube(person.pos)) { // this can happen on rare occasions, due to fp inaccuracy or multiple collisions
		//cout << TXT(rsum) << TXT(sep_dist) << TXT(move_dist) << TXT(room_ix) << TXT(other.pos.str()) << TXT(person.pos.str()) << TXT(bcube.str()) << endl;
		bcube.clamp_pt_xy(person.pos); // just clamp pos so that it doesn't assert later
	}
}

// Note: non-const because this updates room lights
void vect_building_t::ai_room_update(float delta_dir, float dmax, point const &camera_bs, rand_gen_t &rgen) {
	//timer_t timer("Building People Update"); // 0.25ms, mostly iteration overhead, for sparse update with 2-6 people per building (avg for 2 calls city + secondary)

	for (iterator b = begin(); b != end(); ++b) {
		if (!b->interior || b->interior->people.empty() || !b->bcube.closest_dist_less_than(camera_bs, dmax)) continue; // no people or too far away, no updates
		b->all_ai_room_update(rgen, delta_dir);
	}
}

int building_t::get_room_containing_pt(point const &pt) const {
	assert(interior);
	float const wall_thickness(get_wall_thickness());
	bool const in_ext_basement(point_in_extended_basement_not_basement(pt));
	auto rooms_start(in_ext_basement ? interior->ext_basement_rooms_start() : interior->rooms.begin());
	auto rooms_end  (in_ext_basement ? interior->rooms.end() : interior->ext_basement_rooms_start());

	for (auto r = rooms_start; r != rooms_end; ++r) {
		if (r->contains_pt_exp_xy_only(pt, wall_thickness)) {return (r - interior->rooms.begin());} // expand to include point in doorway
	}
	return -1; // room not found
}
building_loc_t building_t::get_building_loc_for_pt(point const &pt) const {
	building_loc_t loc;

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		if (p->contains_pt(pt)) {loc.part_ix = (p - parts.begin()); break;}
	}
	if (loc.part_ix < 0 && point_in_extended_basement(pt)) {loc.part_ix = basement_part_ix;} // use basement part if in extended basement, even though point is outside the cube

	if (interior) { // rooms and stairwells, no elevators yet
		loc.room_ix = get_room_containing_pt(pt);
		if (loc.room_ix >= 0) {loc.floor_ix = get_floor_for_zval(pt.z);}

		for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
			if (s->contains_pt(pt)) {loc.stairs_ix = (s - interior->stairwells.begin()); break;} // Note: stairs_ix is not currently used, except for >= 0 test
		}
		// what about interior->pg_ramp? should that set stairs_ix to stairwells.size()-1?
	}
	return loc;
}
bool building_t::room_containing_pt_has_stairs(point const &pt) const { // Note: only used in building_t::get_all_drawn_window_verts()
	int const room_ix(get_room_containing_pt(pt));
	if (room_ix < 0) return 0; // no room contains this point
	return get_room(room_ix).has_stairs; // Note: used for drawing, can be conservative and return true if any floor of this room has stairs
}

void building_t::register_person_hit(unsigned person_ix, room_object_t const &obj, vector3d const &velocity) {
	if (velocity == zero_vector) return; // stationary object, ignore it
	if (!ai_follow_player())     return; // not in gameplay mode, ignore it
	assert(interior && person_ix < interior->people.size());
	person_t &person(interior->people[person_ix]);

	if (obj.type == TYPE_LG_BALL) { // currently this is the only throwable/dynamic object
		if (obj.zc() < (person.get_z1() + 0.25*person.get_height())) return; // less than 25% up, coll with legs, assume this is kicking a ball that's on the floor
		if (person.retreat_time == 0.0) {maybe_play_zombie_sound(person.pos, person_ix, 1, 1);} // player sound on first retreat: alert_other_zombies=1, high_priority=1
		// Note: this isn't really thread safe, but it should be okay to modify this state while the AI thread is running
		person.retreat_time = RETREAT_TIME; // retreat
		register_achievement("Zombie Bashing");
	}
}

// these must be here to handle deletion of building_nav_graph_t, which is only defined in this file
building_interior_t::building_interior_t() :
	garage_room(-1), ext_basement_hallway_room_id(-1), ext_basement_door_stack_ix(-1), furnace_type(FTYPE_NONE), attic_type(ATTIC_TYPE_RAFTERS),
	door_state_updated(0), is_unconnected(0), ignore_ramp_placement(0), placed_people(0), elevators_disabled(0), attic_access_open(0) {}
building_interior_t::~building_interior_t() {}
