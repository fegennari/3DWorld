// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"

extern building_params_t global_building_params;


void remove_section_from_cube(cube_t &c, cube_t &c2, float v1, float v2, bool xy) { // c is input+output cube, c2 is other output cube
	//if (!(v1 > c.d[xy][0] && v1 < v2 && v2 < c.d[xy][1])) {cout << TXT(v1) << TXT(v2) << TXT(c.d[xy][0]) << TXT(c.d[xy][1]) << TXT(xy) << endl;}
	assert(v1 > c.d[xy][0] && v1 < v2 && v2 < c.d[xy][1]); // v1/v2 must be interior values for cube
	c2 = c; // clone first cube
	c.d[xy][1] = v1; c2.d[xy][0] = v2; // c=low side, c2=high side
}
float cube_rand_side_pos(cube_t const &c, int dim, float min_dist_param, float min_dist_abs, rand_gen_t &rgen) {
	assert(dim < 3);
	assert(min_dist_param < 0.5f); // aplies to both ends
	float const lo(c.d[dim][0]), hi(c.d[dim][1]), delta(hi - lo), gap(max(min_dist_abs, min_dist_param*delta));
	//if ((hi-gap) <= (lo+gap)) {cout << TXT(dim) << TXT(lo) << TXT(hi) << TXT(min_dist_abs) << TXT(delta) << TXT(gap) << endl;}
	return rgen.rand_uniform((lo + gap), (hi - gap));
}

// see global_building_params.window_xspace/window_width
int building_t::get_num_windows_on_side(float xy1, float xy2) const {
	assert(xy1 < xy2);
	building_mat_t const &mat(get_material());
	float tscale(2.0f*mat.get_window_tx()), t0(tscale*xy1), t1(tscale*xy2);
	clip_low_high(t0, t1);
	return round_fp(t1 - t0);
}

// Note: wall should start out equal to the room bcube
void create_wall(cube_t &wall, bool dim, float wall_pos, float fc_thick, float wall_half_thick, float wall_edge_spacing) {
	wall.z1() += fc_thick; // start at the floor
	wall.z2() -= fc_thick; // start at the ceiling
	wall.d[ dim][0] = wall_pos - wall_half_thick;
	wall.d[ dim][1] = wall_pos + wall_half_thick;
	// move a bit away from the exterior wall to prevent z-fighting; we might want to add walls around the building exterior and cut window holes
	wall.d[!dim][0] += wall_edge_spacing;
	wall.d[!dim][1] -= wall_edge_spacing;
}

// Note: assumes edge is not clipped and doesn't work when clipped
bool is_val_inside_window(cube_t const &c, bool dim, float val, float window_spacing, float window_border) {
	float const uv(fract((val - c.d[dim][0])/window_spacing));
	return (uv > window_border && uv < 1.0-window_border);
}

void building_t::gen_interior(rand_gen_t &rgen, bool has_overlapping_cubes) { // Note: contained in building bcube, so no bcube update is needed

	if (!ADD_BUILDING_INTERIORS) return; // disabled
	if (world_mode != WMODE_INF_TERRAIN) return; // tiled terrain mode only
	if (!global_building_params.windows_enabled()) return; // no windows, can't assign floors and generate interior
	//if (has_overlapping_cubes) return; // overlapping cubes buildings are more difficult to handle
	if (!is_cube()) return; // only generate interiors for cube buildings for now
	building_mat_t const &mat(get_material());
	if (!mat.add_windows) return; // not a building type that has generated windows (skip office buildings with windows baked into textures)
	// defer this until the building is close to the player?
	interior.reset(new building_interior_t);
	float const window_vspacing(mat.get_floor_spacing());
	float const floor_thickness(0.1*window_vspacing), fc_thick(0.5*floor_thickness);
	float const doorway_width(0.5*window_vspacing), doorway_hwidth(0.5*doorway_width);
	float const wall_thick(0.5*floor_thickness), wall_half_thick(0.5*wall_thick), wall_edge_spacing(0.05*wall_thick), min_wall_len(4.0*doorway_width);
	float const wwf(global_building_params.get_window_width_fract()), window_border(0.5*(1.0 - wwf)); // (0.0, 1.0)
	unsigned wall_seps_placed[2][2] = {0}; // bit masks for which wall separators have been placed per part, one per {dim x dir}; scales to 32 parts, which should be enough
	vect_cube_t to_split;
	
	// generate walls and floors for each part;
	// this will need to be modified to handle buildings that have overlapping parts, or skip those building types completely
	for (auto p = parts.begin(); p != (parts.end() - has_chimney); ++p) {
		if (is_house && (p - parts.begin()) > 1) break; // houses have at most two parts; exclude garage, shed, porch, porch support, etc.
		float const z_span(p->dz() - floor_thickness);
		assert(z_span > 0.0);
		unsigned const num_floors(round_fp(z_span/window_vspacing)); // round down - no partial floors; add a slight ajustment to account for fp error
		assert(num_floors <= 100); // sanity check
		if (num_floors == 0) continue; // not enough space to add a floor (can this happen?)
		// for now, assume each part has the same XY bounds and can use the same floorplan; this means walls can span all floors and don't need to be duplicated for each floor
		vector3d const psz(p->get_size());
		bool const min_dim(psz.y < psz.x); // hall dim
		float const cube_width(psz[min_dim]);

		if (!is_house && (p+1 == parts.end() || (p+1)->z1() > p->z1()) && cube_width > 4.0*min_wall_len) {
			// building with rectangular slice (no adjacent exterior walls at this level), generate rows of offices
			int const num_windows   (get_num_windows_on_side(p->d[!min_dim][0], p->d[!min_dim][1]));
			int const num_windows_od(get_num_windows_on_side(p->d[ min_dim][0], p->d[ min_dim][1])); // other dim, for use in hallway width calculation
			int const windows_per_room((num_windows > 5) ? 2 : 1); // 1-2 windows per room
			int const num_rooms((num_windows+windows_per_room-1)/windows_per_room); // round up
			bool const partial_room((num_windows % windows_per_room) != 0); // an odd number of windows leaves a small room at the end
			assert(num_rooms >= 0 && num_rooms < 1000); // sanity check
			float const window_hspacing(psz[!min_dim]/num_windows), room_len(window_hspacing*windows_per_room);
			float const hall_width(((num_windows_od & 1) ? 1 : 2)*psz[min_dim]/num_windows_od); // hall either contains 1 (odd) or 2 (even) windows
			float const room_width(0.5f*(cube_width - hall_width)); // rooms are the same size on each side of the hallway
			float const hwall_extend(0.5f*(room_len - doorway_width - wall_thick));
			float const hall_wall_pos[2] = {(p->d[min_dim][0] + room_width), (p->d[min_dim][1] - room_width)};
			vect_cube_t &room_walls(interior->walls[!min_dim]), &hall_walls(interior->walls[min_dim]);
			cube_t rwall(*p); // copy from part; shared zvals, but X/Y will be overwritten per wall
			float const wall_pos(p->d[!min_dim][0] + room_len); // pos of first wall separating first from second rooms
			create_wall(rwall, !min_dim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing); // room walls

			for (int i = 0; i+1 < num_rooms; ++i) { // num_rooms-1 walls
				for (unsigned d = 0; d < 2; ++d) {
					room_walls.push_back(rwall);
					room_walls.back().d[min_dim][!d] = hall_wall_pos[d];
					cube_t hwall(room_walls.back());
					for (unsigned e = 0; e < 2; ++e) {hwall.d[ min_dim][e]  = hall_wall_pos[d] + (e ? 1.0f : -1.0f)*wall_half_thick;}
					for (unsigned e = 0; e < 2; ++e) {hwall.d[!min_dim][e] += (e ? 1.0f : -1.0f)*hwall_extend;}
					if (partial_room && i+2 == num_rooms) {hwall.d[!min_dim][1] -= 1.5*doorway_width;} // pull back a bit to make room for a doorway at the end of the hall
					hall_walls.push_back(hwall); // longer sections that form T-junctions with room walls
				}
				for (unsigned e = 0; e < 2; ++e) {rwall.d[!min_dim][e] += room_len;}
			} // for i
			for (unsigned s = 0; s < 2; ++s) { // add half length hall walls at each end of the hallway
				cube_t hwall(rwall); // copy to get correct zvals
				float const hwall_len((partial_room && s == 1) ? doorway_width : hwall_extend); // hwall for partial room at end is only length doorway_width
				hwall.d[!min_dim][ s] = p->d   [!min_dim][s] + (s ? -1.0f : 1.0f)*wall_edge_spacing; // end at the wall
				hwall.d[!min_dim][!s] = hwall.d[!min_dim][s] + (s ? -1.0f : 1.0f)*hwall_len; // end at first doorway

				for (unsigned d = 0; d < 2; ++d) {
					for (unsigned e = 0; e < 2; ++e) {hwall.d[ min_dim][e] = hall_wall_pos[d] + (e ? 1.0f : -1.0f)*wall_half_thick;}
					hall_walls.push_back(hwall);
				}
			} // for s
		}
		else { // generate random walls using recursive 2D slices
			assert(to_split.empty());
			to_split.push_back(*p);
			float window_hspacing[2] = {0.0};
			
			for (unsigned d = 0; d < 2; ++d) {
				int const num_windows(get_num_windows_on_side(p->d[d][0], p->d[d][1]));
				window_hspacing[d] = psz[d]/num_windows;
				interior->walls[d].reserve(parts.size()); // likely at least this many
			}
			while (!to_split.empty()) {
				cube_t const c(to_split.back());
				to_split.pop_back();
				vector3d const csz(c.get_size());
				bool wall_dim(0); // which dim the room is split by
				if      (csz.y > min_wall_len && csz.x > 1.25*csz.y) {wall_dim = 0;} // split long room in x
				else if (csz.x > min_wall_len && csz.y > 1.25*csz.x) {wall_dim = 1;} // split long room in y
				else {wall_dim = rgen.rand_bool();} // choose a random split dim for nearly square rooms
				if (csz[!wall_dim] < min_wall_len) continue; // not enough space to add a wall (chimney, porch support, etc.)
				float wall_pos(0.0);
				bool const on_edge(c.d[wall_dim][0] == p->d[wall_dim][0] || c.d[wall_dim][1] == p->d[wall_dim][1]); // at edge of the building - make sure walls don't intersect windows
				
				for (unsigned num = 0; num < 20; ++num) { // 20 tries to choose a wall pos that's not inside a window
					wall_pos = cube_rand_side_pos(c, wall_dim, 0.25, wall_thick, rgen);
					if (!on_edge || !is_val_inside_window(*p, wall_dim, wall_pos, window_hspacing[wall_dim], window_border)) break; // okay, keep it
				}
				cube_t wall(c), wall2, wall3; // copy from cube; shared zvals, but X/Y will be overwritten per wall
				create_wall(wall, wall_dim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing);

				// determine if either end of the wall ends at an adjacent part and insert an extra wall there to form a T junction
				for (auto p2 = parts.begin(); p2 != parts.end(); ++p2) {
					for (unsigned dir = 0; dir < 2; ++dir) {
						float const val(c.d[!wall_dim][dir]);
						if (p2 == p) continue; // skip self
						if (p2->d[!wall_dim][!dir] != val) continue; // not adjacent
						if (p2->z1() >= c.z2() || p2->z2() <= c.z1()) continue; // no overlap in Z
						if (p2->d[wall_dim][0] >= wall_pos || p2->d[wall_dim][1] <= wall_pos) continue; // no overlap in wall_dim
						// TODO_INT: what if we try to cut a door into the area where the other part placed a wall?
						if (wall_seps_placed[wall_dim][!dir] & (1 << (p2 - parts.begin()))) continue; // already placed a separator for this part, don't add a duplicate
						wall3.z1() = max(c.z1(), p2->z1()) + fc_thick; // shared Z range
						wall3.z2() = min(c.z2(), p2->z2()) - fc_thick;
						wall3.d[ wall_dim][0] = max(c.d[wall_dim][0], p2->d[wall_dim][0]) + wall_edge_spacing; // shared wall_dim range with slight offset
						wall3.d[ wall_dim][1] = min(c.d[wall_dim][1], p2->d[wall_dim][1]) - wall_edge_spacing;
						wall3.d[!wall_dim][ dir] = val;
						wall3.d[!wall_dim][!dir] = val + (dir ? -1.0 : 1.0)*wall_thick;

						for (unsigned s = 0; s < 2; ++s) { // add doorways to both sides of wall_pos if there's space, starting with the high side
							if (fabs(wall3.d[wall_dim][!s] - wall_pos) > 2.0f*doorway_width) {
								float const doorway_pos(0.5f*(wall_pos + wall3.d[wall_dim][!s])); // centered, for now
								remove_section_from_cube(wall3, wall2, doorway_pos-doorway_hwidth, doorway_pos+doorway_hwidth, wall_dim);
								interior->walls[!wall_dim].push_back(wall2);
							}
						} // for s
						interior->walls[!wall_dim].push_back(wall3);
						wall_seps_placed[wall_dim][dir] |= (1 << (p - parts.begin())); // mark this wall as placed
					} // for dir
				} // for p2
				float const doorway_pos(cube_rand_side_pos(c, !wall_dim, 0.25, doorway_width, rgen));
				remove_section_from_cube(wall, wall2, doorway_pos-doorway_hwidth, doorway_pos+doorway_hwidth, !wall_dim);
				interior->walls[wall_dim].push_back(wall);
				interior->walls[wall_dim].push_back(wall2);

				if (csz[wall_dim] > max(global_building_params.wall_split_thresh, 1.0f)*min_wall_len) {
					for (unsigned d = 0; d < 2; ++d) { // still have space to split in other dim, add the two parts to the stack
						cube_t c_sub(c);
						c_sub.d[wall_dim][d] = wall.d[wall_dim][!d]; // clip to wall pos
						to_split.push_back(c_sub);
					}
				}
			} // end while()
		} // end wall placement
		// add ceilings and floors; we have num_floors+1 separators; the first is only a floor, and the last is only a ceiling
		interior->ceilings.reserve(num_floors);
		interior->floors  .reserve(num_floors);
		float z(p->z1());

		for (unsigned f = 0; f <= num_floors; ++f, z += window_vspacing) {
			cube_t c(*p);
			if (f > 0         ) {c.z1() = z - fc_thick; c.z2() = z; interior->ceilings.push_back(c);}
			if (f < num_floors) {c.z1() = z; c.z2() = z + fc_thick; interior->floors  .push_back(c);}
			c.z1() = z + fc_thick; c.z2() = z + window_vspacing - fc_thick;
			if (f == 0 && p->z1() == bcube.z1() && !doors.empty()) {} // TODO_INT: doors were placed in the prev step; use them to create doorway cutouts on the first floor
			// TODO_INT: add per-floor walls, door cutouts, etc.
		} // for f
	} // for p
}

