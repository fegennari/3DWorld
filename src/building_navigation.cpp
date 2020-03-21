// 3D World - Building Navigation for AIs
// by Frank Gennari 3/14/20

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#pragma warning(disable : 26812) // prefer enum class over enum


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
	unsigned num_rooms, num_stairs;
	vector<node_t> nodes;
	node_t       &get_node(unsigned room)       {assert(room < nodes.size()); return nodes[room];}
	node_t const &get_node(unsigned room) const {assert(room < nodes.size()); return nodes[room];}

	void add_room_to_path(unsigned room1, unsigned room2, vector<point> &path) const {
		path.push_back(get_node(room1).get_conn_pt(room2)); // door pos
		path.push_back(get_node(room2).get_center()); // next room center
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
		
	bool find_path_points(unsigned room1, unsigned room2, bool use_stairs, vector<point> &path) const {
		assert(room1 != room2); // or just return an empty path?
		path.clear();
		//path.push_back(get_node(room1).get_center()); // center of current room
		// if use_stairs==1, or room1 and room2 are in different stacks, the path must include stairs
		//unsigned cur(room1);
		// TODO: run A* algorithm
		//path[i] = get_node(ix).get_center();
		return 0;
	}
};

void building_t::build_nav_graph(building_nav_graph_t &ng) const {

	assert(interior); // too strong?
	if (!interior) return;
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
	building_nav_graph_t ng;
	build_nav_graph(ng);
	return ng.count_connected_components();
}

bool building_t::find_route_to_point(point const &from, point const &to, vector<point> &path) const {

	path.clear();
	if (from == to) return 1; // already there???
	if (!interior)  return 0; // error?
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
	building_nav_graph_t ng; // TODO: cache this in the building (by unique_ptr) or somewhere else?
	build_nav_graph(ng);
	if (!ng.find_path_points(loc1.room_ix, loc2.room_ix, use_stairs, path)) return 0; // failed to find a path
	assert(!path.empty());
	if (path.back() != to) {path.push_back(to);} // add dest pos if not the center of the final room
	return 1;
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
