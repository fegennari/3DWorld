// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"

void building_t::setup_mall_concourse(cube_t &room, bool dim, bool dir, rand_gen_t &rgen) {
	assert(interior);
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	unsigned const num_floors(calc_num_floors(room, floor_spacing, floor_thickness));

	if (num_floors > 1) { // add wall section below the door on lower floors, since the entrace door is on the top level
		float const room_end(room.d[dim][!dir]); // entrance end
		door_t const &door(interior->get_ext_basement_door());
		assert(door.dim == dim);
		cube_t wall(door);
		set_cube_zvals(wall, room.z1()+fc_thick, door.z1()); // below the door
		//set_wall_width(wall, room_end, 0.5*wall_thickness, dim);
		wall.d[dim][!dir] = room_end;
		wall.d[dim][ dir] = room_end + (dir ? 1.0 : -1.0)*wall_thickness;
		interior->walls[dim].push_back(wall);
	}
	// handle upper floors
	for (unsigned f = 1; f < num_floors; ++f) {
		// add side walkways

		// add railings

		// add stairs and escalators
	} // for n
}

void building_t::add_mall_stores(cube_t &room, bool dim, bool dir, rand_gen_t &rgen) {
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness());
	unsigned const num_floors(calc_num_floors(room, floor_spacing, floor_thickness));

	for (unsigned f = 0; f < num_floors; ++f) {
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			// TODO
		} // for d
	} // for n
}

void building_t::add_mall_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t &rooms_to_light) {
	// TODO: railings, potted plants, fountain, benches, tables, chairs, etc.
}

