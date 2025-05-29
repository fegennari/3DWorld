// 3D World - Parking Garages/Structures
// by Frank Gennari 5/27/25

#include "function_registry.h"
#include "buildings.h"
//#include "city_model.h"


// index bits: enable dims: 1=x, 2=y, 4=z | disable cube faces: 8=x1, 16=x2, 32=y1, 64=y2, 128=z1, 256=z2
void building_t::get_parking_struct_ext_walls(vect_cube_with_ix_t &walls, bool exterior_surfaces) const {
	assert(is_parking());
	assert(real_num_parts == (1 + has_basement()));
	cube_t const &part(parts.front()); // above ground part
	float const floor_spacing(get_window_vspace()), floor_thick(get_floor_thickness()), wall_thick(get_wall_thickness());
	float const lower_wall_height(0.35*floor_spacing), upper_wall_height(0.15*floor_spacing);
	unsigned num_floors(calc_num_floors(part, floor_spacing, floor_thick));
	assert(num_floors > 0);
	vect_cube_t wall_parts, door_cuts, temp;

	for (tquad_with_ix_t const &d : doors) { // find all doors on the ground floor
		if (d.type == tquad_with_ix_t::TYPE_RDOOR) continue; // roof access door - skip
		cube_t dbc(d.get_bcube());
		bool const dim(dbc.dy() < dbc.dx()), dir(part.get_center_dim(dim) < dbc.get_center_dim(dim));
		dbc.expand_in_dim(dim, 2.0*wall_thick);
		door_cuts.push_back(dbc);
		// add extra vertical wall segments to either side of the door
		cube_t wall(part);
		set_cube_zvals(wall, (part.z1() + lower_wall_height), (part.z1() + floor_spacing - upper_wall_height));
		wall.d[dim][!dir] = wall.d[dim][dir] + (dir ? -1.0 : 1.0)*wall_thick; // set wall thickness

		for (unsigned side = 0; side < 2; ++side) {
			float const door_edge(dbc.d[!dim][side]);
			cube_t w(wall);
			w.d[!dim][!side] = door_edge;
			w.d[!dim][ side] = door_edge + (side ? 1.0 : -1.0)*0.75*floor_spacing;
			w.intersect_with_cube(part); // don't extend outside the part
			unsigned const sf[2][2] = {{8, 16}, {32, 64}}; // {dim x dir}
			unsigned face_mask(0);
			if (exterior_surfaces) {face_mask = 3 + sf[dim][!dir] + sf[!dim][!side];} // outside and exposed end
			else                   {face_mask = (1 << dim) + sf[dim][dir];} // inside only
			walls.emplace_back(w, face_mask);
		} // for side
	}
	// one extra floor; each wall slice goes from the ceiling of the level below to the top of the wall of the level above, except for the two ends
	for (unsigned f = 0; f <= num_floors; ++f) {
		float const zval(part.z1() + f*floor_spacing), wall_bot(zval - upper_wall_height), wall_top(zval + lower_wall_height);
		cube_t slice(part);
		set_cube_zvals(slice, max(part.z1(), wall_bot), min(wall_top, part.z2()));
		if (exterior_surfaces) {walls.emplace_back(slice, 3);} // XY exterior walls around entire perimeter; don't need to handle doors here
		cube_t hole(slice);
		hole.expand_by_xy(-wall_thick);
		cube_t sides[4]; // {-y, +y, center -x, center +x}
		subtract_cube_xy(slice, hole, sides);
		unsigned const face_masks[4] = {127-64, 127-32, 127-16, 127-8}; // enable XYZ but skip all XY but {+y, -y, +x, -x} in XY

		for (unsigned n = 0; n < 4; ++n) { // exterior: top and bottom; interior: inside faces
			if (f < 2) { // cut out slots for doors on the ground floor lower wall and second floor upper wall
				for (unsigned lu = 0; lu < 2; ++lu) { // split into lower and upper sections
					if (f == (lu ? num_floors : 0)) continue; // no lower/upper section for this floor
					cube_t wall(sides[n]);
					if (lu) {wall.z1() = zval + 0.5*floor_thick;} // upper wall
					else    {wall.z2() = zval - 0.5*floor_thick;} // lower wall
					wall_parts.clear();
					subtract_cubes_from_cube(wall, door_cuts, wall_parts, temp, 2); // check zval overlap
					unsigned const face_mask(exterior_surfaces ? 4 : face_masks[n]);
					for (cube_t const &w : wall_parts) {walls.emplace_back(w, face_mask);}
				} // for lu
			}
			else {walls.emplace_back(sides[n], (exterior_surfaces ? 4 : face_masks[n]));} // add entire wall
		} // for n
	} // for f
}

void building_t::add_parking_struct_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	unsigned num_floors, unsigned &nlights_x, unsigned &nlights_y, float &light_delta_z, light_ix_assign_t &light_ix_assign)
{
	assert(has_room_geom());
	cube_t const &part(parts.front()); // above ground part
	vect_room_object_t &objs(interior->room_geom->objs);
	// add vertical support pillars spanning all floors
	// TODO
}
