// 3D World - Building Underground Shopping Malls
// by Frank Gennari 11/03/2024

#include "function_registry.h"
#include "buildings.h"

float building_t::get_mall_floor_spacing(cube_t const &room) const { // special function that allows for larger than normal floor spacing
	assert(has_mall());
	return room.dz()/interior->num_extb_floors;
}
cube_t building_t::get_mall_center(cube_t const &room) const {
	bool const dim(interior->extb_wall_dim);
	float const ww_width(0.25*room.get_sz_dim(!dim));
	cube_t center(room); // center open area
	center.expand_in_dim(!dim, -1.0*ww_width); // short dim
	center.expand_in_dim( dim, -2.0*ww_width); // long dim
	return center;
}
void building_t::get_mall_open_areas(cube_t const &room, vect_cube_t &openings) const {
	bool const dim(interior->extb_wall_dim);
	cube_t const center(get_mall_center(room));
	float const length(center.get_sz_dim(dim)), width(center.get_sz_dim(!dim)), gap(0.75*width), half_gap(0.5*gap);
	unsigned const num_splits(max(1U, unsigned(round_fp(0.25*length/width))));
	float const step((length + gap)/num_splits);
	float pos(center.d[dim][0] - half_gap);
	openings.reserve(num_splits);

	for (unsigned n = 0; n < num_splits; ++n) {
		cube_t c(center);
		c.d[dim][0] = pos + half_gap;
		pos += step;
		c.d[dim][1] = pos - half_gap;
		openings.push_back(c);
	}
}

void building_t::setup_mall_concourse(cube_t const &room, bool dim, bool dir, rand_gen_t &rgen) {
	assert(interior);
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	
	// add wall section below the door on lower floors, since the entrace door is on the top level
	float const room_end(room.d[dim][!dir]); // entrance end
	door_t const &door(interior->get_ext_basement_door());
	assert(door.dim == dim);
	cube_t wall(door);
	set_cube_zvals(wall, room.z1()+fc_thick, door.z1()); // below the door
	wall.d[dim][!dir] = room_end; // doesn't exactly match the regular wall width, but not visible anyway
	wall.d[dim][ dir] = room_end + (dir ? 1.0 : -1.0)*wall_thickness; // extend into the room
	if (wall.dz() > 0.0) {interior->walls[dim].push_back(wall);}

	// handle upper floors
	if (interior->num_extb_floors == 1) return; // single floor; nothing else to do
	vect_cube_t openings;
	get_mall_open_areas(room, openings);
	assert(!openings.empty());

	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const z(room.z1() + f*floor_spacing), zc(z - fc_thick), zf(z + fc_thick);
		
		// add side walkways
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			cube_t ww(room);
			ww.d[!dim][!d] = openings.front().d[!dim][d];
			interior->add_ceil_floor_pair(ww, zc, z, zf);
		}
		for (unsigned n = 0; n <= openings.size(); ++n) { // ends and gaps between openings
			cube_t ww(openings.front());
			ww.d[dim][0] = ((n == 0              ) ? room.d[dim][0] : openings[n-1].d[dim][1]);
			ww.d[dim][1] = ((n == openings.size()) ? room.d[dim][1] : openings[n  ].d[dim][0]);
			interior->add_ceil_floor_pair(ww, zc, z, zf);
		}
		// add stairs and escalators
		// TODO
	} // for f
}

void building_t::add_mall_stores(cube_t &room, bool dim, bool dir, rand_gen_t &rgen) {
	float const floor_spacing(get_mall_floor_spacing(room)), floor_thickness(get_floor_thickness());

	for (unsigned f = 0; f < interior->num_extb_floors; ++f) {
		for (unsigned d = 0; d < 2; ++d) { // each side of concourse
			// TODO
		} // for d
	} // for n
}

void building_t::add_mall_stairs() {
	if (!has_mall()) return;
	room_t const &room(interior->get_extb_start_room());
	door_t const &door(interior->get_ext_basement_door());
	bool const dim(interior->extb_wall_dim), dir(interior->extb_wall_dir);
	assert(door.dim == dim);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const room_id(interior->ext_basement_hallway_room_id);
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness()), wall_thickness(get_wall_thickness());
	// add stairs under the door if needed
	float const upper_floor_zval(room.z2() - floor_spacing + fc_thick), delta_z(door.z1() - upper_floor_zval);
	unsigned const num_steps(max(0, (int)ceil(NUM_STAIRS_PER_FLOOR*delta_z/get_floor_ceil_gap())));

	if (num_steps > 0) {
		float const step_height(delta_z/num_steps), step_len(2.0*step_height);
		float const wall_edge(room.d[dim][!dir]), dsign(dir ? 1.0 : -1.0), front_step_dist(dsign*step_len);
		cube_t stair(door);
		set_cube_zvals(stair, door.z1()-step_height, door.z1());
		stair.d[dim][!dir] = wall_edge; // starts at wall
		stair.d[dim][ dir] = wall_edge + 2.0*front_step_dist;

		for (unsigned n = 0; n < num_steps; ++n) { // top down
			stair.expand_in_dim(!dim, step_len);  // widen sides
			stair.d[dim][dir] += front_step_dist; // add length
			objs.emplace_back(stair, TYPE_STAIR, room_id, dim, dir, 0, 1.0, SHAPE_STAIRS_FAN);
			stair.translate_dim(2, -step_height);
		}
		// add railings against the wall
		cube_t railing;
		set_cube_zvals(railing, upper_floor_zval, door.z1());
		railing.d[dim][!dir] = wall_edge + dsign*1.0*wall_thickness;
		railing.d[dim][ dir] = wall_edge + dsign*1.6*wall_thickness;

		for (unsigned d = 0; d < 2; ++d) {
			railing.d[!dim][!d] = door .d[!dim][d] + (d ? 1.0 : -1.0)*1.5*wall_thickness;
			railing.d[!dim][ d] = stair.d[!dim][d] + (d ? 1.0 : -1.0)*1.5*wall_thickness;
			objs.emplace_back(railing, TYPE_RAILING, room_id, !dim, !d, RO_FLAG_HAS_EXTRA, 1.0, SHAPE_CUBE, GOLD); // no balusters
		}
	}
}

void building_t::add_mall_lower_floor_lights(room_t const &room, unsigned room_id, unsigned lights_start, light_ix_assign_t &light_ix_assign) {
	float const floor_spacing(get_mall_floor_spacing(room)), fc_thick(get_fc_thickness());
	vect_cube_t openings;
	get_mall_open_areas(room, openings);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_end(objs.size());
	assert(lights_start <= objs_end);

	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const zc(room.z1() + f*floor_spacing - fc_thick); // bottom of ceiling

		for (unsigned i = lights_start; i < objs_end; ++i) {
			room_object_t const &obj(objs[i]);
			if (obj.type != TYPE_LIGHT)       continue; // should this ever fail?
			if (has_bcube_int(obj, openings)) continue; // skip lights over the openings
			room_object_t light(obj);
			light.translate_dim(2, (zc - obj.z2()));
			light.obj_id = light_ix_assign.get_ix_for_light(light);
			objs.push_back(light);
		} // for i
	} // for f
}

void building_t::add_mall_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t &rooms_to_light) {
	float const floor_spacing(get_mall_floor_spacing(room)), window_vspace(get_window_vspace()), fc_thick(get_fc_thickness());
	float const wall_thickness(get_wall_thickness()), trim_thick(get_trim_thickness());
	vect_cube_t openings;
	get_mall_open_areas(room, openings);
	vect_room_object_t &objs(interior->room_geom->objs);

	// add walkway railings
	for (unsigned f = 1; f < interior->num_extb_floors; ++f) { // skip first floor
		float const z(room.z1() + f*floor_spacing);
		float const zb1(z - fc_thick - trim_thick), zb2(z + fc_thick + 2.0*trim_thick), zt1(z + 0.45*window_vspace), zt2(zt1 + 0.04*window_vspace);
		colorRGBA const bar_color(LT_GRAY);

		for (cube_t const &opening : openings) {
			for (unsigned dim = 0; dim < 2; ++dim) {
				for (unsigned dir = 0; dir < 2; ++dir) {
					cube_t railing(opening);
					set_wall_width(railing, opening.d[dim][dir], 0.0, dim);
					cube_t bar(railing);
					// bottom bar
					bar.expand_by_xy(0.4*wall_thickness);
					set_cube_zvals(bar, zb1, zb2);
					objs.emplace_back(bar, TYPE_METAL_BAR, room_id, !dim, 0, 0, 1.0, SHAPE_CUBE, bar_color);
					// top bar
					bar.expand_by_xy(0.1*wall_thickness); // additional expand
					set_cube_zvals(bar, zt1, zt2);
					objs.emplace_back(bar, TYPE_METAL_BAR, room_id, !dim, 0, 0, 1.0, SHAPE_CUBE, bar_color);
					// glass
					cube_t panel(railing);
					set_cube_zvals(panel, zb2, zt1);
					panel.expand_in_dim(dim, 0.1*wall_thickness);
					objs.emplace_back(panel, TYPE_INT_WINDOW, room_id, dim, 0, 0, 1.0);
				} // for dir
			} // for dim
			// add corner bars; will need to add to both ends and both dims once there are gaps for stairs and escalators
			float const vbar_hwidth(0.35*wall_thickness);

			for (unsigned n = 0; n < 4; ++n) {
				bool const dirs[2] = {(n&1), (n&2)}; // x, y
				cube_t bar;
				set_cube_zvals(bar, zb2, zt1);
				for (unsigned d = 0; d < 2; ++d) {set_wall_width(bar, opening.d[d][dirs[d]], vbar_hwidth, d);}
				objs.emplace_back(bar, TYPE_METAL_BAR, room_id, 0, 0, 0, 1.0, SHAPE_CUBE, bar_color);
			}
		}
	} // for f
	// TODO: potted plants, fountain, benches, tables, chairs, etc.
	//cube_t walk_area(room);
	//walk_area.expand_by_xy(-wall_thickness);
}

