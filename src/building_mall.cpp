// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"

void building_t::setup_mall_concourse(cube_t &room, bool dim, bool dir, rand_gen_t &rgen) {
	assert(interior);
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	unsigned const num_floors(calc_num_floors(room, floor_spacing, floor_thickness));
	if (num_floors == 1) return; // single floor; nothing else to do
	//cube_t walk_area(room);
	//walk_area.expand_by_xy(-wall_thickness);
	// add wall section below the door on lower floors, since the entrace door is on the top level
	float const room_end(room.d[dim][!dir]); // entrance end
	door_t const &door(interior->get_ext_basement_door());
	assert(door.dim == dim);
	cube_t wall(door);
	set_cube_zvals(wall, room.z1()+fc_thick, door.z1()); // below the door
	wall.d[dim][!dir] = room_end;
	wall.d[dim][ dir] = room_end + (dir ? 1.0 : -1.0)*wall_thickness;
	interior->walls[dim].push_back(wall);
	// handle upper floors
	float const ww_width(0.25*room.get_sz_dim(!dim));
	cube_t center(room); // center open area
	center.expand_in_dim(!dim, -1.0*ww_width); // short dim
	center.expand_in_dim( dim, -1.5*ww_width); // long dim

	for (unsigned f = 1; f < num_floors; ++f) {
		float const zc(room.z1() + f*floor_spacing), z1(zc - fc_thick), z2(zc + fc_thick);
		
		// add side walkways
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			cube_t ww(room);
			ww.d[!dim][!d] = center.d[!dim][d];
			interior->add_ceil_floor_pair(ww, z1, zc, z2);
		}
		for (unsigned d = 0; d < 2; ++d) { // each end of concourse
			cube_t ww(center);
			ww.d[dim][ d] = room  .d[dim][d];
			ww.d[dim][!d] = center.d[dim][d];
			interior->add_ceil_floor_pair(ww, z1, zc, z2);
		}
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
	// add lights to the underside of upper floor walkways
}

