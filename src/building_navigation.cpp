// 3D World - Building Navigation for AIs
// by Frank Gennari 3/14/20

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#pragma warning(disable : 26812) // prefer enum class over enum


// Note: this should go into building_t/buildings.h at some point, but is temporarily here
class building_nav_graph_t {
	struct node_t { // represents one room or one stairwell
		bool has_exit, is_hallway, is_stairs;
		vector<unsigned> conn_rooms;
		node_t() : has_exit(0), is_hallway(0), is_stairs(0) {}

		void add_conn_room(unsigned room) {
			for (auto i = conn_rooms.begin(); i != conn_rooms.end(); ++i) {if (*i == room) return;} // ignore duplicates
			conn_rooms.push_back(room);
		}
	};
	unsigned num_rooms, num_stairs;
	vector<node_t> nodes;
	node_t &get_node(unsigned room) {assert(room < nodes.size()); return nodes[room];}
public:
	building_nav_graph_t() : num_rooms(0), num_stairs(0) {}

	void set_num_rooms(unsigned num_rooms_, unsigned num_stairs_) {
		num_rooms  = num_rooms_;
		num_stairs = num_stairs_;
		nodes.resize(num_rooms + num_stairs);
		for (unsigned n = num_rooms; n < nodes.size(); ++n) {nodes[n].is_stairs = 1;}
	}
	void mark_hallway(unsigned room) {get_node(room).is_hallway = 1;}
	void mark_exit   (unsigned room) {get_node(room).has_exit   = 1;}

	void connect_stairs(unsigned room, unsigned stairs) {
		assert(room < num_rooms && stairs < num_stairs);
		connect_rooms(room, (num_rooms + stairs));
	}

	void connect_rooms(unsigned room1, unsigned room2) { // graph is bidirectional
		assert(room1 < num_rooms && room2 < num_rooms);
		get_node(room1).add_conn_room(room2);
		get_node(room2).add_conn_room(room1);
	}
	bool is_fully_connected() const {
		return 0; // TODO
	}
	bool find_path(unsigned room1, unsigned room2, bool different_floor, vector<unsigned> &path) const {
		assert(room1 != room2); // or just return an empty path?
		path.clear();
		// if different_floor==1, or room1 and room2 are in different stacks, the path must include stairs
		unsigned cur(room1);
		// TODO
		return 0;
	}
};

void building_t::build_nav_graph() {

	assert(interior); // too strong?
	if (!interior) return;
	float const wall_width(0.5*get_floor_thickness());
	unsigned const num_rooms(interior->rooms.size()), num_stairs(interior->stairwells.size());
	building_nav_graph_t ng; // TODO: take this by reference or cache it inside interior by unique_ptr
	ng.set_num_rooms(num_rooms, num_stairs);

	for (unsigned r = 0; r < num_rooms; ++r) {
		room_t const &room(interior->rooms[r]);
		if (room.is_hallway) {ng.mark_hallway(r);}
		cube_t c(room);
		c.expand_by_xy(wall_width); // to include adjacent doors

		for (auto d = doors.begin(); d != doors.end(); ++d) { // exterior doors
			if (!d->is_exterior_door() || d->type == tquad_with_ix_t::TYPE_RDOOR) continue;
			if (c.intersects(d->get_bcube())) {ng.mark_exit(r); break;}
		}
		for (auto d = interior->doors.begin(); d != interior->doors.end(); ++d) {
			cube_t dc(*d);
			dc.expand_by_xy(wall_width); // to include adjacent rooms

			for (unsigned r2 = 0; r2 < num_rooms; ++r2) {
				if (r2 == r) continue; // skip self
				if (dc.intersects(interior->rooms[r2])) {ng.connect_rooms(r, r2); break;}
			}
		} // for d
		for (unsigned s = 0; s < num_stairs; ++s) {
			if (c.intersects(interior->stairwells[s])) {ng.connect_stairs(r, s);}
		}
		//for (unsigned e = 0; e < interior->elevators.size(); ++e) {} // elevators are not yet used by AIs so are ignored here
	} // for r
}
