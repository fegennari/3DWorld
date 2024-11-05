// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"

float building_t::get_mall_floor_spacing(cube_t const &room) const { // special function that allows for larger than normal floor spacing
	assert(has_mall());
	return room.dz()/interior->num_extb_floors;
}
cube_t building_t::get_mall_center(cube_t const &room) const {
	bool const dim(room.dx() < room.dy());
	float const ww_width(0.25*room.get_sz_dim(!dim));
	cube_t center(room); // center open area
	center.expand_in_dim(!dim, -1.0*ww_width); // short dim
	center.expand_in_dim( dim, -1.5*ww_width); // long dim
	return center;
}

void building_t::setup_mall_concourse(cube_t &room, bool dim, bool dir, rand_gen_t &rgen) {
	assert(interior);
	float const floor_spacing(get_mall_floor_spacing(room)), floor_thickness(get_floor_thickness()), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	if (interior->num_extb_floors == 1) return; // single floor; nothing else to do
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
	cube_t const center(get_mall_center(room));

	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const z(room.z1() + f*floor_spacing), zc(z - fc_thick), zf(z + fc_thick);
		
		// add side walkways
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			cube_t ww(room);
			ww.d[!dim][!d] = center.d[!dim][d];
			interior->add_ceil_floor_pair(ww, zc, z, zf);
		}
		for (unsigned d = 0; d < 2; ++d) { // each end of concourse
			cube_t ww(center);
			ww.d[dim][ d] = room  .d[dim][d];
			ww.d[dim][!d] = center.d[dim][d];
			interior->add_ceil_floor_pair(ww, zc, z, zf);
		}
		// add railings

		// add stairs and escalators
	} // for n
}

void building_t::add_mall_stores(cube_t &room, bool dim, bool dir, rand_gen_t &rgen) {
	float const floor_spacing(get_mall_floor_spacing(room)), floor_thickness(get_floor_thickness());

	for (unsigned f = 0; f < interior->num_extb_floors; ++f) {
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			// TODO
		} // for d
	} // for n
}

void building_t::add_mall_lower_floor_lights(room_t const &room, unsigned room_id, unsigned lights_start, light_ix_assign_t &light_ix_assign) {
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness());
	cube_t const center(get_mall_center(room));
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_end(objs.size());
	assert(lights_start <= objs_end);

	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const zc(room.z1() + f*floor_spacing - fc_thick); // bottom of ceiling

		for (unsigned i = lights_start; i < objs_end; ++i) {
			room_object_t const &obj(objs[i]);
			if (obj.type != TYPE_LIGHT) continue; // should this ever fail?
			if (center.intersects(obj)) continue; // skip lights over the open center
			room_object_t light(obj);
			light.translate_dim(2, (zc - obj.z2()));
			light.obj_id = light_ix_assign.get_ix_for_light(light);
			objs.push_back(light);
		} // for i
	} // for f
}

void building_t::add_mall_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t &rooms_to_light) {
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness());
	// TODO: railings, potted plants, fountain, benches, tables, chairs, etc.
}

