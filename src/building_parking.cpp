// 3D World - Parking Garages/Structures
// by Frank Gennari 5/27/25

#include "function_registry.h"
#include "buildings.h"
//#include "city_model.h"


// index bits: enable dims: 1=x, 2=y, 4=z | disable cube faces: 8=x1, 16=x2, 32=y1, 64=y2, 128=z1, 256=z2
void building_t::get_parking_garage_ext_walls(vect_cube_with_ix_t &walls, bool exterior_surfaces) {
	assert(is_parking());
	assert(real_num_parts == (1 + has_basement()));
	cube_t const &part(parts.front()); // above ground part
	float const floor_spacing(get_window_vspace()), floor_thick(get_floor_thickness()), wall_thick(get_wall_thickness());
	unsigned num_floors(calc_num_floors(part, floor_spacing, floor_thick));
	assert(num_floors > 0);

	// one extra floor; each wall slice goes from the ceiling of the level below to the top of the wall of the level above, except for the two ends
	for (unsigned f = 0; f <= num_floors; ++f) {
		float const zval(part.z1() + f*floor_spacing), wall_bot(zval - 0.15*floor_spacing), wall_top(zval + 0.35*floor_spacing);
		cube_t slice(part);
		set_cube_zvals(slice, max(part.z1(), wall_bot), min(wall_top, part.z2()));
		if (exterior_surfaces) {walls.emplace_back(slice, 3);} // XY exterior walls around entire perimeter
		cube_t hole(slice);
		hole.expand_by_xy(-wall_thick);
		cube_t sides[4]; // {-y, +y, center -x, center +x}
		subtract_cube_xy(slice, hole, sides);
		unsigned const face_masks[4] = {127-64, 127-32, 127-16, 127-8}; // enable XYZ but skip all XY but {+y, -y, +x, -x} in XY

		for (unsigned n = 0; n < 4; ++n) {
			walls.emplace_back(sides[n], (exterior_surfaces ? 4 : face_masks[n])); // exterior: top and bottom; interior: inside faces
		}
	} // for f
}
