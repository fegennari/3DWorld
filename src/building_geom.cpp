// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "function_registry.h"
#include "buildings.h"

using std::string;


cube_t grass_exclude1, grass_exclude2;

extern bool begin_motion;
extern int animate2;
extern float grass_width, CAMERA_RADIUS, fticks;
extern building_params_t global_building_params;


/*static*/ float building_t::get_scaled_player_radius() {return CAMERA_RADIUS*global_building_params.player_coll_radius_scale;}
float get_bldg_player_height() {return (CAMERA_RADIUS + get_player_height());} // feet to eyes

void building_t::set_z_range(float z1, float z2) {
	bcube.z1() = z1; bcube.z2() = z2;
	adjust_part_zvals_for_floor_spacing(bcube);
	if (!parts.empty()) {parts[0].z1() = z1; parts[0].z2() = z2;}
}
building_mat_t const &building_t::get_material() const {return global_building_params.get_material(mat_ix);}

void building_t::gen_rotation(rand_gen_t &rgen) {
	float const max_rot_angle(get_material().max_rot_angle);
	if (max_rot_angle == 0.0) return;
	float const rot_angle(rgen.rand_uniform(0.0, TO_RADIANS*max_rot_angle)); // max_rot_angle is specified in degrees
	rot_sin = sin(rot_angle);
	rot_cos = cos(rot_angle);
	parts.clear();
	parts.push_back(bcube); // this is the actual building base
	set_bcube_from_rotated_cube(parts.back());
}
point building_t::get_inv_rot_pos(point const &pos) const {
	if (!is_rotated()) return pos;
	point pos_rot(pos);
	do_xy_rotate_inv(bcube.get_cube_center(), pos_rot);
	return pos_rot;
}

void building_t::set_bcube_from_rotated_cube(cube_t const &bc) {
	point const center(bc.get_cube_center());

	for (unsigned i = 0; i < 4; ++i) {
		point corner(bc.d[0][i&1], bc.d[1][i>>1], bc.d[2][i&1]);
		do_xy_rotate(center, corner);
		if (i == 0) {bcube.set_from_point(corner);} else {bcube.union_with_pt(corner);} // Note: detail cubes are excluded
	}
}

cube_t building_t::calc_parts_bcube() const {
	assert(!parts.empty());
	cube_t bc(parts[0]);
	for (auto i = parts.begin()+1; i != parts.end(); ++i) {bc.union_with_cube(*i);} // update bcube
	return bc;
}
void building_t::calc_bcube_from_parts() {
	cube_t const bc(calc_parts_bcube());
	if (!is_rotated()) {bcube = bc;} else {set_bcube_from_rotated_cube(bc);} // apply rotation if needed and set bcube
	coll_bcube = bcube; // initial value
}

void building_t::move_door_to_other_side_of_wall(tquad_with_ix_t &door, float dist_mult, bool invert_normal) const {
	cube_t const c(door.get_bcube());
	bool const dim(c.dy() < c.dx()), dir(door.get_norm()[dim] > 0.0); // closest cube side dir
	float door_shift(0.0);
	if (invert_normal) {swap(door.pts[0], door.pts[1]); swap(door.pts[2], door.pts[3]);} // swap vertex order to invert normal
	float const door_edge(door.pts[0][dim]);

	if (door.type == tquad_with_ix_t::TYPE_RDOOR) { // not on a wall, use shift relative to floor/wall thickness
		door_shift = (dir ? 1.0 : -1.0)*get_door_shift_dist();
	}
	else {
		float const large_val(bcube.dz());
		door_shift = large_val; // start with a large value

		for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) { // find the part that this door was added to (inc garages and sheds)
			if (p->z1() > c.z1() || p->z2() < c.z2()) continue; // Z-range not contained
			if (is_basement(p)) continue; // skip the basement; probably not needed with the above check
			float edge_pos(p->d[dim][dir]); // start with bcube edge

			if (!is_cube()) { // find polygon edge
				float const toler(0.1*get_wall_thickness()); // add a bit of tolerance to account for FP error
				vect_point const &points(get_part_ext_verts(p - parts.begin()));

				for (auto i = points.begin(); i != points.end(); ++i) {
					point const &p1(*i), &p2((i == points.begin()) ? points.back() : *(i-1));
					if (fabs(p1[dim] - p2[dim]) > toler) continue; // edge not in this dim
					edge_pos = 0.5*(p1[dim] + p2[dim]);
					break; // should be only one
				}
			}
			float const dist(door_edge - edge_pos); // signed
			if (fabs(dist) < fabs(door_shift)) {door_shift = dist;}
		} // for p
		if (door_shift == large_val) { // not found
			std::cerr << TXT(c.str()) << TXT(dim) << TXT(dir) << TXTi(real_num_parts) << TXT(is_cube()) << endl;
			assert(0);
		}
	}
	door_shift *= -(1.0 + dist_mult); // reflect on other side
	// double door_shift until it changes the value; needed for buildings far from the origin; use tolerance of 1E-7 since max precision is around 6E-8
	while (fabs(door_shift/door_edge) < 1.0E-7) {door_shift *= 2.0;}
	for (unsigned n = 0; n < door.npts; ++n) {door.pts[n][dim] += door_shift;} // move to opposite side of wall
}

void building_t::clip_door_to_interior(tquad_with_ix_t &door) const {

	cube_t clip_cube(door.get_bcube());
	float const dz(clip_cube.dz());
	float xy_border(0.0), z_border(0.0);
	if      (door.type == tquad_with_ix_t::TYPE_GDOOR) {xy_border = 0.016; z_border = 0.02;} // garage door
	else if (door.is_building_door())                  {xy_border = 0.04;  z_border = 0.03;} // building door
	else {xy_border = 0.06; z_border = 0.03;} // house door
	// clip off bottom for floor if clip_to_floor==1 and not a roof door; somewhat arbitrary, should we use interior->floors.back().z2() instead?
	if (door.type != tquad_with_ix_t::TYPE_RDOOR) {clip_cube.z1() += 0.1*get_floor_thickness();}
	clip_cube.z2() -= z_border*dz;
	bool const dim(clip_cube.dx() < clip_cube.dy()); // border dim
	clip_cube.expand_in_dim(dim, -xy_border*clip_cube.get_sz_dim(dim)); // shrink by border
	for (unsigned n = 0; n < door.npts; ++n) {clip_cube.clamp_pt(door.pts[n]);}
}

void building_t::get_garage_dim_dir(cube_t const &garage, bool &dim, bool &dir) const { // works with interior or exterior garages
	dim = (garage.dx() < garage.dy()); // long dim
	if (street_dir > 0 && bool((street_dir-1)>>1) == dim) {dir = !((street_dir-1)&1);} // use street_dir if it's set and dims agree
	else {dir = (garage.get_center_dim(dim) < bcube.get_center_dim(dim));} // assumes the garage is at an exterior wall and doesn't occupy the entire house width
}

int building_t::get_part_ix_containing_pt(point const &pt) const {
	auto parts_end(get_real_parts_end_inc_sec());

	for (auto i = parts.begin(); i != parts_end; ++i) { // includes garage/shed
		if (i->contains_pt(pt)) {return (i - parts.begin());}
	}
	return -1; // not found
}
cube_t building_t::get_part_containing_pt(point const &pt) const {
	auto parts_end(get_real_parts_end_inc_sec());

	for (auto i = parts.begin(); i != parts_end; ++i) { // includes garage/shed
		if (i->contains_pt(pt)) {return *i;}
	}
	// we can get here in rare cases due to FP precision problems, or we're on the roof access stairs; find the closest cube to pt, which should be very close
	float dmin_sq(0.0);
	cube_t closest;

	for (auto i = parts.begin(); i != parts_end; ++i) {
		float const dist_sq(p2p_dist(pt, i->closest_pt(pt)));
		if (dmin_sq == 0.0 || dist_sq < dmin_sq) {dmin_sq = dist_sq; closest = *i;}
	}
	assert(!closest.is_all_zeros()); // must be found
	return closest;
}
int building_t::get_part_ix_containing_cube(cube_t const &c) const {
	auto parts_end(get_real_parts_end_inc_sec());

	for (auto i = parts.begin(); i != parts_end; ++i) { // could call b.get_part_for_room() if we had a room_t
		if (i->contains_cube(c)) return (i - parts.begin());
	}
	return -1; // not found
}
cube_t building_t::get_part_containing_cube(cube_t const &c) const {
	int const part_ix(get_part_ix_containing_cube(c));
	return ((part_ix < 0) ? cube_t() : parts[part_ix]);
}
bool building_t::check_cube_within_part_sides(cube_t const &c) const {
	if (is_cube())   return 1; // assume caller has already checked parts/bcube
	int const part_ix(get_part_ix_containing_cube(c));
	if (part_ix < 0) return 0;
	vect_point const &points(get_part_ext_verts(part_ix));

	for (unsigned n = 0; n < 4; ++n) { // test all 4 top corners; if any is outside, placement is invalid; correct as long as points forms a convex shape
		if (!point_in_polygon_2d(c.d[0][n&1], c.d[1][n>>1], points.data(), points.size())) return 0;
	}
	return 1; // contained
}
bool building_t::check_pt_within_part_sides(point const &p) const {
	if (is_cube())   return 1; // assume caller has already checked parts/bcube
	int const part_ix(get_part_ix_containing_pt(p));
	if (part_ix < 0) return 0;
	vect_point const &points(get_part_ext_verts(part_ix));
	return point_in_polygon_2d(p.x, p.y, points.data(), points.size());
}

void building_t::adjust_part_zvals_for_floor_spacing(cube_t &c) const {
	if (!EXACT_MULT_FLOOR_HEIGHT) return;
	float const floor_spacing(get_window_vspace()), dz(c.dz());
	assert(dz > 0.0);
	assert(floor_spacing > 0.0);
	float const num_floors(dz/floor_spacing);
	int const targ_num_floors(max(1, round_fp(num_floors)));
	c.z2() += floor_spacing*(targ_num_floors - num_floors); // ensure c.dz() is an exact multiple of num_floors
}

// split c against existing cubes and add its non-overlapping pieces
void split_cubes_recur(cube_t c, vect_cube_t &cubes, unsigned search_start, unsigned search_end) {

	for (unsigned i = search_start; i < search_end; ++i) {
		cube_t &sc(cubes[i]);
		// while this generally holds, it can fail in rare cases (or used to?), I assume due to floating-point error or some such thing, so this check has been disabled
		//assert(sc.z2() >= c.z2()); // assumes cubes are ordered descending by ztop
		if (!sc.intersects_no_adj(c)) continue;
		if (sc.contains_cube(c)) return; // contained, done (remove all of c)

		if (sc.z1() == c.z1() && sc.z2() > c.z2() && c.contains_cube_xy(sc)) { // sc covers the base of c and spans c's z-range
			assert(sc.z1() < c.z2());
			sc.z1() = c.z2(); // create a stack instead
			continue;
		}
		// find a split plane
		for (unsigned d = 0; d < 2; ++d) { // dim
			for (unsigned e = 0; e < 2; ++e) { // dir
				float const split_pos(cubes[i].d[d][e]); // Note: can't use sc reference as it may have been invalidated by a push_back()
				
				if (split_pos > c.d[d][0] && split_pos < c.d[d][1]) { // this plane splits c
					cube_t hi_c(c);
					hi_c.d[d][0] = split_pos; // hi part
					c.   d[d][1] = split_pos; // lo part
					// recursively split the hi part starting at this cube if it's not contained; this split plane will no longer be active
					if (!cubes[i].contains_cube(hi_c)) {split_cubes_recur(hi_c, cubes, i, search_end);}
					if ( cubes[i].contains_cube(c)) return; // done (optimization)
				}
			} // for e
		} // for d
	} // for i
	cubes.push_back(c);
}

void building_t::gen_geometry(int rseed1, int rseed2) {

	if (!is_valid()) return; // invalid building
	if (!parts.empty()) {adjust_part_zvals_for_floor_spacing(parts.front());}
	cube_t const base(parts.empty() ? bcube : parts.back());
	assert(base.is_strictly_normalized());
	parts.clear();
	details.clear();
	roof_tquads.clear();
	doors.clear();
	interior.reset();
	building_mat_t const &mat(get_material());
	rand_gen_t rgen;
	rgen.set_state(123+rseed1, 345*rseed2);
	ao_bcz2         = bcube.z2(); // capture z2 before union with roof and detail geometry (which increases building height)
	ground_floor_z1 = bcube.z1(); // record before adding basement
	wall_color      = mat.wall_color; // start with default wall color

	if (btype == BTYPE_UNSET) { // building type not customized
		if (is_house)                                       {btype = BTYPE_HOUSE   ;} // may be flatted as BTYPE_MULT_FAM in gen_house()
		else if (rgen.rand_probability(mat.apartment_prob)) {btype = ((rseed1 & 1) ? BTYPE_HOTEL : BTYPE_APARTMENT);}
		else if (is_cube() && (rseed1&15) == 0)             {btype = BTYPE_HOSPITAL;} // 1/16 the time
		else if (is_cube() && (rseed1&15) == 1)             {btype = BTYPE_SCHOOL  ;} // 1/16 the time
		else                                                {btype = BTYPE_OFFICE  ;} // office is the default for non-residential buildings
	}
	assign_name(rgen);
	
	if (is_house) {
		gen_house(base, rgen);
		return;
	}
	// determine building shape (cube, cylinder, other)
	if (rgen.rand_probability(mat.round_prob)) {num_sides = MAX_CYLIN_SIDES;} // max number of sides for drawing rounded (cylinder) buildings
	else if (rgen.rand_probability(mat.cube_prob)) {num_sides = 4;} // cube
	else { // N-gon
		num_sides = mat.min_sides;
		if (mat.min_sides != mat.max_sides) {num_sides += (rgen.rand() % (1 + abs((int)mat.max_sides - (int)mat.min_sides)));}
	}
	bool const was_cube(is_cube()); // before num_sides increase due to ASF

	if (num_sides >= 6 && mat.max_fsa > 0.0) { // at least 6 sides
		flat_side_amt = max(0.0f, min(0.45f, rgen.rand_uniform(mat.min_fsa, mat.max_fsa)));
		if (flat_side_amt > 0.0 && rot_sin == 0.0) {start_angle = rgen.rand_uniform(0.0, TWO_PI);} // flat side, not rotated: add random start angle to break up uniformity
	}
	if ((num_sides == 3 || num_sides == 4 || num_sides == 6) && mat.max_asf > 0.0 && rgen.rand_probability(mat.asf_prob)) { // triangles/cubes/hexagons
		alt_step_factor = max(0.0f, min(0.99f, rgen.rand_uniform(mat.min_asf, mat.max_asf)));
		if (alt_step_factor > 0.0 && !(num_sides&1)) {half_offset = 1;} // chamfered cube/hexagon
		if (alt_step_factor > 0.0) {num_sides *= 2;}
	}

	// determine the number of levels and splits
	unsigned num_levels(mat.min_levels);

	if (mat.min_levels < mat.max_levels) { // have a range of levels
		if (was_cube || rgen.rand_bool()) {num_levels += rgen.rand() % (mat.max_levels - mat.min_levels + 1);} // only half of non-cubes are multilevel (unless min_level > 1)
	}
	if (mat.min_level_height > 0.0) {num_levels = max(mat.min_levels, min(num_levels, unsigned(bcube.dz()/mat.min_level_height)));}
	num_levels = max(num_levels, 1U); // min_levels can be zero to apply more weight to 1 level buildings
	bool const do_split(num_levels < 4 && is_cube() && rgen.rand_probability(mat.split_prob)); // don't split buildings with 4 or more levels, or non-cubes
	float const height(base.dz()), floor_spacing(get_window_vspace());

	if (num_levels == 1) { // single level
		if (do_split) { // generate L, T, or U shape
			split_in_xy(base, rgen);

			if (has_city_trees()) { // see if we can place a tree in the courtyard
				point const center(bcube.get_cube_center());
				cube_t place_area(center);
				place_area.expand_by_xy(0.05f*(bcube.dx() + bcube.dy()));
				if (!has_bcube_int(place_area, parts)) {tree_pos = place_area.get_cube_center(); tree_pos.z = ground_floor_z1;}
			}
			if (0 && btype == BTYPE_OFFICE && is_cube()) { // single level cube
				btype = BTYPE_PRISON;
				assign_name(rgen); // re-assign a name
			}
			gen_details(rgen, 0);
		}
		else { // single part, entire cube/cylinder
			parts.push_back(base);
			unsigned const rand_val(rgen.rand()), num_floors(round_fp(height/floor_spacing));
			if ((rand_val&3) != 0) {maybe_add_special_roof(rgen);} // 75% chance
			
			// consider a possible vertical split of the floorplan into two parts
			if (!interior_enabled() || num_floors < 2) {} // no interior, or single floor, can't split vertically
			else if (rgen.rand_probability(global_building_params.split_stack_floorplan_prob)) {
				// while this works, it doesn't seem to add much value, it only creates odd geometry and makes connecting stairs/elevators difficult
				// two stacked parts of the same x/y dimensions but different interior floorplans
				parts.push_back(base);
				parts[0].z2() = base.z1() + rgen.rand_uniform(0.4, 0.6)*height; // split in Z: parts[0] is the bottom, parts[1] is the top
				adjust_part_zvals_for_floor_spacing(parts[0]);
				parts[1].z1() = parts[0].z2();
			}
			else if (btype == BTYPE_OFFICE && is_cube() && num_floors <= 4) { // <= 4 floors
				btype = (rgen.rand_bool() ? BTYPE_WAREHOUSE : BTYPE_FACTORY); // make this a factory or warehouse
				assign_name(rgen); // re-assign a name
			}
			// parking garages are <= 8 floors and large footprint
			else if (btype == BTYPE_OFFICE && is_cube() && num_floors <= 8 && roof_type == ROOF_TYPE_FLAT && min(bcube.dx(), bcube.dy()) > 12.0*floor_spacing) {
				btype = BTYPE_PARKING; // make a parking garage
				assign_name(rgen); // re-assign a name
			}
			else if (is_cube() && num_floors >= 3 && !is_hospital() && rgen.rand_probability(global_building_params.retail_floorplan_prob)) { // 3+ floors, consider retail
				rand_gen_t rgen2(rgen); // create a new rgen to avoid affecting the other building parameters when this option is changed
				retail_floor_levels = 1;
				// only create a tall retail area if there are at least 4 floors (2 below and 2 above), otherwise the top part won't have central stairs to extend below
				if (num_floors >= 4 && rgen2.rand_probability(global_building_params.two_floor_retail_prob)) {retail_floor_levels = 2;}
				parts.push_back(base);
				parts[0].z2() = parts[1].z1() = base.z1() + retail_floor_levels*floor_spacing; // split in Z: parts[0] is the bottom, parts[1] is the top
			}
			gen_details(rgen, 1);
		}
		finish_gen_geometry(rgen, 0);
		return; // for now the bounding cube
	}
	// generate building levels and splits
	float const dz(height/num_levels), abs_min_edge_move(0.5*floor_spacing), align_dist(0.1*floor_spacing); // dz is same as door width
	bool const not_too_small(min(bcube.dx(), bcube.dy()) > 4.0*abs_min_edge_move);
	assert(height > 0.0);

	if (!do_split && not_too_small && (rgen.rand()&3) < (was_cube ? 2 : 3) && !has_windows()) {
		// oddly shaped multi-sided overlapping sections (50% chance for cube buildings and 75% chance for others)
		vector2d const sz(base.get_size_xy());
		parts.reserve(num_levels); // at least this many

		for (unsigned i = 0; i < num_levels; ++i) { // generate overlapping cube levels, tallest to shortest
			cube_t bc(base); // copy from base to start, keep z1
			bc.z2() = base.z1() + (num_levels - i)*dz; // z2
			if (i > 0) {bc.z2() += dz*rgen.rand_uniform(-0.45, 0.45); bc.z2() = min(bc.z2(), base.z2());}
			if (i > 0) {assert(bc.z2() <= parts.back().z2());}
			assert(bc.is_strictly_normalized());
			adjust_part_zvals_for_floor_spacing(bc);
			bool valid(0);

			// make 200 attempts to generate a cube that isn't contained in any existing cubes; most of the time should pass the first time, so it should rarely fail
			for (unsigned n = 0; n < 200; ++n) {
				for (unsigned d = 0; d < 2; ++d) { // x,y
					float const mv_lo(rgen.rand_uniform(-0.2, 0.45)), mv_hi(rgen.rand_uniform(-0.2, 0.45));
					if (mv_lo > 0.0) {bc.d[d][0] = base.d[d][0] + max(abs_min_edge_move, mv_lo*sz[d]);}
					if (mv_hi > 0.0) {bc.d[d][1] = base.d[d][1] - max(abs_min_edge_move, mv_hi*sz[d]);}
				}
				if (!bc.is_strictly_normalized()) continue; // can fail when building is very small compared to floor spacing/abs_min_edge_move
				bool contained(0);
				for (auto p = parts.begin(); p != parts.end(); ++p) {contained |= p->contains_cube(bc);}
				if (!contained) {valid = 1; break;} // success
			} // for n
			if (!valid) break; // remove this part and end the building here

			// align the edge with the nearby edges of any previous parts to avoid small jogs that can't fit a room
			for (unsigned d = 0; d < 2; ++d) {
				for (unsigned e = 0; e < 2; ++e) {
					float &edge(bc.d[d][e]);

					for (cube_t const &p : parts) {
						if (fabs(edge - p.d[d][e]) < align_dist) {edge = p.d[d][e]; break;}
					}
				} // for e
			} // for d
			if (i == 0 || !is_cube()) {parts.push_back(bc); continue;} // no splitting
			split_cubes_recur(bc, parts, 0, parts.size()); // split this cube against all previously added cubes and remove overlapping areas
		} // for i
		if (!parts.empty()) { // at least one part was placed; should (almost) always be true?
			// this tends to create some strange interior rooms, so skip this building type for secondary buildings with windows that the player can see through
			if (parts.size() > 1) {has_complex_floorplan = 1;}
			parts.shrink_to_fit(); // optional
			std::reverse(parts.begin(), parts.end()); // highest part should be last so that it gets the roof details (remember to update basement_part_ix if needed)
			calc_bcube_from_parts(); // update bcube
			gen_details(rgen, 1);
			finish_gen_geometry(rgen, 1);
			return;
		}
	}
	parts.resize(num_levels);

	for (unsigned i = 0; i < num_levels; ++i) {
		cube_t &bc(parts[i]);
		if (i == 0) {bc = base;} // use full building footprint
		else {
			cube_t const &prev(parts[i-1]);
			float const shift_mult(was_cube ? 1.0 : 0.5); // half the shift for non-cube buildings

			for (unsigned d = 0; d < 2; ++d) {
				float const len(prev.get_sz_dim(d)), min_edge_len((0.2f/shift_mult)*(bcube.get_sz_dim(d)));
				bool const inv(rgen.rand_bool());

				for (unsigned e = 0; e < 2; ++e) {
					float delta(0.0);
					if (rgen.rand()&3) {delta = shift_mult*rgen.rand_uniform(0.1, 0.4);} // 25% chance of no shift, 75% chance of 20-40% shift
					bc.d[d][e] = prev.d[d][e] + (e ? -delta : delta)*len;
				}
				for (unsigned E = 0; E < 2; ++E) {
					bool const e((E != 0) ^ inv); // no dir favoritism for 20% check
					if (bc.get_sz_dim(d) < min_edge_len) {bc.d[d][e] = prev.d[d][e];} // if smaller than 20% base width, revert the change
				}
			}
			bc.z1() = prev.z2(); // z1
		}
		bc.z2() = bc.z1() + dz; // z2
		bc.normalize(); // handle XY inversion due to shift
	} // for i
	for (unsigned i = 1; i < num_levels; ++i) {
		float const ddz(rgen.rand_uniform(-0.35*dz, 0.35*dz)); // random shift in z height
		parts[i-1].z2() += ddz;
		adjust_part_zvals_for_floor_spacing(parts[i-1]);
		parts[i].z1() = parts[i-1].z2(); // make top and bottom parts align
	}
	adjust_part_zvals_for_floor_spacing(parts[num_levels-1]); // last one
	max_eq(bcube.z2(), parts[num_levels-1].z2()); // adjust bcube if needed
	cube_t const split_cube(parts.back());

	if (do_split && split_cube.dx() > 0.4*bcube.dx() && split_cube.dy() > 0.4*bcube.dy()) { // generate L, T, or U shape for top level if not too small
		parts.pop_back();
		split_in_xy(split_cube, rgen);
		if (num_levels <= 3) {gen_details(rgen, 0);}
	}
	else {
		if ((rgen.rand()&3) != 0) {maybe_add_special_roof(rgen);} // 67% chance
		if (num_levels <= 3) {gen_details(rgen, 1);}
	}
	finish_gen_geometry(rgen, 0);
}

void building_t::setup_damage_vals() {
	if (has_mall()) return; // no damage for malls
	point const center(bcube.get_cube_center());
	float const floor_spacing(get_window_vspace());
	float const water_rand_val(fract((center.x + center.y + center.z)/floor_spacing)), crack_rand_val(fract(1.5*(center.x + center.y + center.z)/floor_spacing));
	water_damage = max(0.0f, (water_rand_val - 0.5f)); // 50% of buildings have up to 50% water damage
	crack_damage = ((crack_rand_val < 0.4) ? 2.5*crack_rand_val : 0.0); // random amount of cracks, 40% of the time
}

void building_t::create_per_part_ext_verts() {
	if (per_part_ext_verts.empty() && !is_cube()) { // generate exterior verts for each part if not yet created
		per_part_ext_verts.resize(parts.size());
		for (unsigned p = 0; p < parts.size(); ++p) {building_draw_utils::calc_poly_pts(*this, bcube, parts[p], per_part_ext_verts[p]);}
	}
}
void building_t::finish_gen_geometry(rand_gen_t &rgen, bool has_overlapping_cubes) { // for office buildings
	if (coll_bcube.is_all_zeros()) {coll_bcube = bcube;} // calculate if it hasn't been calculated yet
	if (global_building_params.add_office_basements) {maybe_add_basement(rgen);}
	assert(parts.size() > 0 && parts.size() < 256);
	real_num_parts = uint8_t(parts.size()); // no parts can be added after this point
	create_per_part_ext_verts();
	parts_generated = 1;

	// apartments, hotels, hospitals, and schools must have a primary hallway; if we can't add a hallway to the first/bottom part, then make it an office instead
	if ((is_apt_or_hotel() || is_hospital() || is_school()) && !can_use_hallway_for_part(0)) {
		btype = BTYPE_OFFICE;
		assign_name(rgen); // re-assign a name
	}
	gen_interior(rgen, has_overlapping_cubes);

	if (global_building_params.windows_enabled() && global_building_params.add_office_br_basements && !has_complex_floorplan) { // skip complex floorplan buildings
		extend_underground_basement(rgen); // maybe add door inside basement and connected backrooms area; rgen is copied, not modified
	}
	if (interior) {interior->assign_door_conn_rooms();} // must be after adding extended basement and before adding room geom
	if (interior) {interior->finalize();}
	gen_building_doors_if_needed(rgen);
	setup_damage_vals();
}

void building_t::split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen) {

	// generate L, T, U, H, +, O shape
	point const llc(seed_cube.get_llc());
	vector2d const sz(seed_cube.get_size_xy());
	bool const allow_courtyard(seed_cube.dx() < 1.6*seed_cube.dy() && seed_cube.dy() < 1.6*seed_cube.dx()); // AR < 1.6:1
	int const shape(rgen.rand()%(allow_courtyard ? 10 : 9)); // 0-9
	bool const has_hole(shape == 9);
	bool const is_hpo(shape >= 7);
	bool const dim(rgen.rand_bool()); // {x,y}
	bool const dir(is_hpo ? 1 : rgen.rand_bool()); // {neg,pos} - H/+/O shapes are symmetric and always pos
	float const smin(0.2), smax(has_hole ? 0.35 : 0.4); // outer and inner split sizes
	float const div(is_hpo ? rgen.rand_uniform(smin, smax) : rgen.rand_uniform(0.3, 0.7)); // split pos in 0-1 range
	float const s1(rgen.rand_uniform(smin, smax)), s2(rgen.rand_uniform(1.0-smax, 1.0-smin)); // split pos in 0-1 range
	float const dpos(llc[dim] + div*sz[dim]), spos1(llc[!dim] + s1*sz[!dim]), spos2(llc[!dim] + s2*sz[!dim]); // split pos in cube space
	float const dpos2(llc[dim] + (1.0 - div)*sz[dim]); // other end - used for H, +, and O
	unsigned const start(parts.size()), num((shape >= 9) ? 4 : ((shape >= 6) ? 3 : 2));
	parts.resize(start+num, seed_cube);
	parts[start+0].d[dim][ dir] = dpos; // full width part (except +)
	parts[start+1].d[dim][!dir] = dpos; // partial width part (except +)
	has_courtyard = (has_hole && seed_cube.z1() == ground_floor_z1); // courtyard is only on the ground floor

	switch (shape) {
	case 0: case 1: case 2: case 3: // L
		parts[start+1].d[!dim][shape>>1] = ((shape&1) ? spos2 : spos1);
		break;
	case 4: case 5: // T
		parts[start+1].d[!dim][0] = spos1;
		parts[start+1].d[!dim][1] = spos2;
		break;
	case 6: // U
		parts[start+2].d[ dim][!dir] = dpos; // partial width part
		parts[start+1].d[!dim][1   ] = spos1;
		parts[start+2].d[!dim][0   ] = spos2;
		break;
	case 7: { // H
		parts[start+1].d[ dim][ dir] = dpos2;
		parts[start+1].d[!dim][ 0  ] = spos1; // middle part
		parts[start+1].d[!dim][ 1  ] = spos2;
		parts[start+2].d[ dim][!dir] = dpos2; // full width part
		break;
	}
	case 8: { // +
		parts[start+0].d[!dim][ 0  ] = spos1;
		parts[start+0].d[!dim][ 1  ] = spos2;
		parts[start+2].d[!dim][ 0  ] = spos1;
		parts[start+2].d[!dim][ 1  ] = spos2;
		parts[start+1].d[ dim][ dir] = dpos2; // middle part
		parts[start+2].d[ dim][!dir] = dpos2; // partial width part
		break;
	}
	case 9: { // O (courtyard)
		parts[start+2].d[ dim][!dir] = dpos; // partial width part (same as U)
		parts[start+1].d[ dim][ dir] = dpos2;
		parts[start+2].d[ dim][ dir] = dpos2;
		parts[start+1].d[!dim][1   ] = spos1;
		parts[start+2].d[!dim][0   ] = spos2;
		parts[start+3].d[ dim][!dir] = dpos2; // other full width part
		break;
	}
	default: assert(0);
	}
}

bool get_largest_xy_dim(cube_t const &c) {return (c.dy() > c.dx());}

bool building_t::is_valid_door_pos(cube_t const &door, float door_width, bool dim) const {
	if (interior) {
		if (interior->is_blocked_by_stairs_or_elevator(door, door_width)) return 0;
		cube_t test_cube(door);
		test_cube.expand_in_dim(dim, door_width);
		if (interior->is_cube_close_to_doorway(test_cube, cube_t(), 0.0, 1, 1))   return 0; // check interior doors: null room, inc_open=1, check_open_dir=1
		if (door.zc() < ground_floor_z1 && check_cube_intersect_walls(test_cube)) return 0; // required for basement doors
	}
	cube_t test_cube(door);
	test_cube.expand_by_xy(2.0*door_width); // must have at least two door widths of space

	for (auto const &door : doors) { // check other exterior doors
		if (test_cube.intersects(door.get_bcube())) return 0;
	}
	if (has_chimney == 2 && test_cube.intersects(get_fireplace())) return 0; // too close to fireplace (Note: door is actually placed first, likely has no effect)
	return 1;
}

cube_t building_t::place_door(cube_t const &base, bool dim, bool dir, float door_height, float door_center, float door_pos,
	float door_center_shift, float width_scale, bool can_fail, bool opens_up, rand_gen_t &rgen, unsigned floor_ix) const
{
	float const door_width(width_scale*door_height), door_half_width(0.5*door_width);
	if (can_fail && base.get_sz_dim(!dim) < 2.0*door_width) return cube_t(); // part is too small to place a door
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), fc_thickness(get_fc_thickness());
	float const door_shift(get_door_shift_dist()), base_lo(base.d[!dim][0]), base_hi(base.d[!dim][1]);
	bool const calc_center(door_center == 0.0); // door not yet calculated
	bool const centered(door_center_shift == 0.0 || hallway_dim == (uint8_t)dim); // center doors connected to primary hallways
	// ideally we want the front (first) door to connect to the stairs in a multi-family house, but the stairs may be in the back, so we allow the back door as well
	bool const is_basement(base.zc() < ground_floor_z1), is_front_door(/*doors.empty() &&*/ !is_basement);
	bool const pref_near_stairs(interior && is_front_door && multi_family);
	unsigned const base_num_tries(10), num_tries((pref_near_stairs ? 2 : 1)*base_num_tries);
	cube_t door;
	door.z1() = base.z1() + floor_ix*floor_spacing + fc_thickness; // floor level relative to base part
	door.z2() = min((door.z1() + door_height), (base.z2() - fc_thickness));

	for (unsigned n = 0; n < num_tries; ++n) { // make up to 10 tries to place a valid door
		if (calc_center) { // add door to first part of house/building
			float const offset(centered ? 0.5 : rgen.rand_uniform(0.5-door_center_shift, 0.5+door_center_shift));
			door_center = offset*base_lo + (1.0 - offset)*base_hi;
		}
		if (calc_center || door_pos == 0.0) {door_pos = base.d[dim][dir];}

		if (pref_near_stairs && n < 15) { // check if room has stairs or is a hallway; basement stairs stairs don't count
			point center;
			center.z     = door.zc();
			center[ dim] = door_pos - wall_thickness*(dir ? 1.0 : -1.0); // move into the house interior
			center[!dim] = door_center;
			int const room_ix(get_room_containing_pt(center));
			//assert(room_ix >= 0);
			if (room_ix < 0) continue; // should never fail, but if it does, then the door placement must be bad
			room_t const &room(interior->get_room(room_ix));
			
			if (!room.has_stairs_on_floor(1)) { // check for stairs on second floor to exclude basement stairs; we know there must be at least two floors
				if (n < base_num_tries || !room.is_hallway) continue; // allow hallways as well for n=[10, 14]
			}
		}
		if (interior && (!has_pri_hall() || is_basement) && !opens_up) { // not on a hallway - check distance to interior walls to make sure the door has space to open
			auto const &walls(interior->walls[!dim]); // perpendicular to door
			float const door_lo(door_center - 1.2*door_half_width), door_hi(door_center + 1.2*door_half_width); // pos along wall with a small expand
			float const dpos_lo(door_pos    -     door_half_width), dpos_hi(door_pos    +     door_half_width); // expand width of the door

			for (auto w = walls.begin(); w != walls.end(); ++w) {
				if (w->z1() > door.z2() || w->z2() < door.z1())         continue; // wrong part/floor
				if (w->d[ dim][0] > dpos_hi || w->d[ dim][1] < dpos_lo) continue; // not ending at same wall as door
				if (w->d[!dim][0] > door_hi || w->d[!dim][1] < door_lo) continue; // not intersecting door
				// Note: since we know that all rooms are wider than the door width, we know that we have space for a door on either side of the wall
				// move the door so that it doesn't open into the end of the wall, but clamp to the base bounds in case the condition above doesn't hold
				float const lo_dist(w->d[!dim][0] - door_lo), hi_dist(door_hi - w->d[!dim][1]);
				if (lo_dist < hi_dist) {door_center = min((door_center + lo_dist), (base_hi - door_half_width));}
				else                   {door_center = max((door_center - hi_dist), (base_lo + door_half_width));}
				break;
			} // for w
		}
		door.d[dim][!dir] = door_pos + door_shift*(dir ? 1.0 : -1.0); // move slightly away from the building to prevent z-fighting
		door.d[dim][ dir] = door.d[dim][!dir]; // make zero size in this dim
		set_wall_width(door, door_center, door_half_width, !dim);

		// we're free to choose the door pos, and have the interior, so we can check if the door is in a good location;
		// if we get here on the last iteration, just keep the door even if it's an invalid location
		if (calc_center && interior && !centered && !is_valid_door_pos(door, door_width, dim)) {
			if (can_fail && n+1 == num_tries) return cube_t(); // last iteration, fail
			continue; // try a new location
		}
		break; // done
	} // for n
	return door;
}

bool building_t::check_walkway_door_clearance(cube_t const &c, bool dim) const {
	if (interior->is_blocked_by_stairs_or_elevator(c, 0.0, 0, 2)) return 0; // dmin=0.0, elevators_only=0, no_check_enter_exit=2

	// check wall ends; shouldn't need to check the other dim because rooms should be at least a doorway width wide
	for (cube_t const &w : interior->walls[!dim]) {
		if (w.intersects_no_adj(c)) return 0;
	}
	return 1;
}
bool building_t::add_walkway_door(building_walkway_geom_t &walkway, bool dir, unsigned part_ix) {
	float const door_width(get_office_ext_doorway_width()), door_shift(get_door_shift_dist()), floor_spacing(get_window_vspace());
	float const door_height(get_floor_ceil_gap()); // not using get_door_height() because we want to span the entire height, since there's no interior wall above
	bool const dim(walkway.dim);
	cube_t const &wbc(walkway.bcube);
	unsigned const num_floors(round_fp(wbc.dz()/floor_spacing));
	assert(num_floors > 0);
	float const center(wbc.get_center_dim(!dim));
	cube_t door;
	set_wall_width(door, center, 0.5*door_width, !dim);
	door.d[dim][!dir] = wbc .d[dim][!dir] + (dir ? 1.0 : -1.0)*door_shift; // move away from the building to prevent z-fighting
	door.d[dim][ dir] = door.d[dim][!dir]; // make zero size in this dim
	float zval(wbc.z1() + get_fc_thickness()); // bottom of lowest level door

	if (interior) { // check for clearance; should this be done per-floor?
		// walkway doors don't have a full door width of front clearance since they open outward and are also split in half
		float const front_clearance(max(0.5f*door_width, get_min_front_clearance_inc_people()));
		set_cube_zvals(door, zval, wbc.z2()); // full walkway height
		cube_t door_exp(door);
		door_exp.d[dim][!dir] -= (dir ? 1.0 : -1.0)*front_clearance;
		
		if (!check_walkway_door_clearance(door_exp, dim)) { // blocked, try points to the sides
			bool sdir(interior->rooms.size() & 1); // start with a pseudo random offset direction
			float const shift_amt(0.1*door_width);
			unsigned const num_shifts(0.5*max(0.0f, (wbc.get_sz_dim(!dim) - door_width))/shift_amt); // for each side, rounded down
			bool success(0);

			for (unsigned n = 0; n < num_shifts && !success; ++n) {
				for (unsigned d = 0; d < 2; ++d, sdir ^= 1) {
					float const shift((d ? 1.0 : -1.0)*(n+1)*shift_amt);
					cube_t cand(door_exp);
					cand.translate_dim(!dim, shift);
					if (!check_walkway_door_clearance(cand, dim)) continue; // still bad
					door.translate_dim(!dim, shift);
					success = 1;
					break;
				} // for d
			} // for n
			if (!success) return 0; // no valid door placement
		}
	}
	for (unsigned f = 0; f < num_floors; ++f, zval += floor_spacing) {
		set_cube_zvals(door, zval, zval+door_height);
		add_door(door, part_ix, dim, dir, 1, 0, 0, 1); // for_office_building=1, roof_access=0, courtyard=0, for_walkway=1
	}
	for (unsigned d = 0; d < 2; ++d) {walkway.door_bounds[!dir][d] = door.d[!dim][d];} // dir is relative to building and opposite walkway, so must be swapped
	have_walkway_ext_door = 1;
	return 1; // success
}

bool building_t::clip_cube_to_parts(cube_t &c, vect_cube_t &cubes) const { // use for fences
	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) {
		if (!p->intersects(c)) continue;
		subtract_cube_from_cube(c, *p, cubes, 1); // clear_out=1
		// each part should either not intersect the fence, or clip off one side of it; if this doesn't hold for some reason (bad house size, etc.), then return failured
		assert(!cubes.empty());
		if (cubes.size() > 1) return 0; // failed
		c = cubes.back();
	} // for p
	return 1; // success
}

void add_cube_top(cube_t const &c, vect_tquad_with_ix_t &tquads, unsigned type) {
	tquad_t tquad(4); // quad
	tquad.pts[0].assign(c.x1(), c.y1(), c.z2());
	tquad.pts[1].assign(c.x2(), c.y1(), c.z2());
	tquad.pts[2].assign(c.x2(), c.y2(), c.z2());
	tquad.pts[3].assign(c.x1(), c.y2(), c.z2());
	tquads.emplace_back(tquad, type);
}

bool building_t::add_chimney(bool two_parts, bool stacked_parts, bool hipped_roof[4], float roof_dz[4], unsigned force_dim[2], rand_gen_t &rgen) {
	unsigned part_ix(0); // default for single part and stacked_parts=1

	if (two_parts) { // start by selecting a part (if two parts)
		float const v0(parts[0].get_volume()), v1(parts[1].get_volume());
		if      (v0 > 2.0*v1) {part_ix = 0;} // choose larger part 0
		else if (v1 > 2.0*v0) {part_ix = 1;} // choose larger part 1
		else {part_ix = rgen.rand_bool();} // close in area - choose a random part
	} // else if stacked_parts use the first/bottom part
	assert(roof_dz[part_ix] > 0.0);
	unsigned const fdim(force_dim[part_ix]);
	cube_t const &part(parts[part_ix]);
	bool dir(rgen.rand_bool()), dim((fdim < 2) ? fdim : get_largest_xy_dim(part)); // use longest side if not forced

	if (stacked_parts) {
		if (parts[0].d[dim][dir] != parts[1].d[dim][dir]) {dir ^= 1;} // stacks not aligned, try other dir
		
		if (parts[0].d[dim][dir] != parts[1].d[dim][dir]) { // stacks still not aligned
			if (fdim < 2) return 0; // dim was forced, fail
			dim ^= 1; // try the other dim
			if (parts[0].d[dim][dir] != parts[1].d[dim][dir]) {dir ^= 1;} // try other dim, original dir
			if (parts[0].d[dim][dir] != parts[1].d[dim][dir]) return 0; // tried all four dim/dir values - fail
		}
	}
	else if (two_parts) {
		cube_t parts_union(parts[0]); // use union of parts to account for larger bcubes due to other objects or rotated house
		parts_union.union_with_cube(parts[1]);
		if (part.d[dim][dir] != parts_union.d[dim][dir]) {dir ^= 1;} // force dir to be on the edge of the house bcube (not at a point interior to the house)
	}
	float chimney_dz((hipped_roof[part_ix] ? 0.5 : 1.0)*roof_dz[part_ix]); // lower for hipped roof
	cube_t c(part);
	float const sz1(c.get_sz_dim(!dim)), sz2(c.get_sz_dim(dim)), center(c.get_center_dim(!dim));
	float const window_vspace(get_window_vspace()), sz1_sane(max(sz1, 2.0f*window_vspace)), chimney_depth(0.03f*(sz1_sane + sz2));
	float shift(0.0);

	if (stacked_parts && part_ix == 0 && parts[1].d[dim][dir] == part.d[dim][dir]) { // shared wall for stacked parts
		c.z2() = parts[1].z2(); // extend upward to top part
		max_eq(chimney_dz, (hipped_roof[1] ? 0.5f : 1.0f)*roof_dz[1]); // must extend above the upper part's roof
	}
	if ((rgen.rand()%3) != 0) { // make the chimney non-centered 67% of the time
		shift = sz1*rgen.rand_uniform(0.1, 0.25); // select a shift in +/- (0.1, 0.25) - no small offset from center
		if (rgen.rand_bool()) {shift = -shift;}
	}
	float chimney_height(rgen.rand_uniform(1.25, 1.5)*chimney_dz);

	if (rgen.rand()&3) { // chimney outside the bounds of the house, 75% of the time
		float const hwidth(0.04f*sz1_sane);
		c.d[dim][!dir]  = c.d[dim][dir] + 0.001*(dir ? 1.0 : -1.0)*chimney_depth; // slight shift to avoid Z-fighting
		c.d[dim][ dir] += (dir ? 1.0 : -1.0)*chimney_depth;

		for (unsigned n = 0; n < 10; ++n) { // make up to 10 tries to place the chimney
			set_wall_width(c, (center + shift), hwidth, !dim); // set chimney width
			// add bottom section, which will be the outside of the fireplace
			cube_t fplace(c);
			fplace.z1() = part.z1();
			fplace.z2() = c.z1() = part.z1() + 0.85*window_vspace; // 85% of floor height
			fplace.expand_in_dim(!dim, max(1.0f*hwidth, 0.6f*chimney_depth)); // widen for fireplace
			// check if blocked by a door, wall, stairs, or garage, and skip if it is
			cube_t test_cube(fplace);
			test_cube.expand_by_xy(0.5*hwidth);
			max_eq(test_cube.z2(), part.z2()); // extend upward to include chimney and avoid upper floor doors
			bool bad_pos(0);
			for (auto d = doors.begin(); d != doors.end() && !bad_pos; ++d) {bad_pos = test_cube.intersects(d->get_bcube());} // check exterior doors

			if (interior) { // check walls
				for (auto w = interior->walls[!dim].begin(); w != interior->walls[!dim].end() && !bad_pos; ++w) {
					if (w->z1() > fplace.z2() || w->z2() < fplace.z2()) continue; // not on the ground floor
					bad_pos = test_cube.intersects(*w);
				}
				if (!bad_pos) { // check stairs
					cube_t fireplace_ext(fplace);
					fireplace_ext.d[dim][!dir] += (dir ? -1.0 : 1.0)*1.5*chimney_depth;
					bad_pos = interior->is_blocked_by_stairs_or_elevator(fireplace_ext);
				}
				if (!bad_pos && has_int_garage) { // check garage, which shouldn't have a fireplace
					bad_pos = interior->get_garage_room().intersects(test_cube);
				}
			}
			if (bad_pos) { // failed to place chimney
				shift = sz1*rgen.rand_uniform(-0.3, 0.3); // try a new random position, with a bit more room to move
				continue;
			}
			parts.push_back(fplace); // Note: invalidates part
			union_with_coll_bcube(fplace);
			if (rgen.rand_bool()) {c.d[!dim][0] = fplace.d[!dim][0]; c.d[!dim][1] = fplace.d[!dim][1];} // widen chimney to include entire fireplace (for more modern houses)
			has_chimney = 2; // flag as exterior chimney
			break; // done
		} // for n
		if (!has_chimney) return 0; // failed to place
	}
	else { // chimney inside the bounds of the house; placement can't fail
		set_wall_width(c, (center + shift), 0.05*sz1_sane, !dim); // set chimney width
		c.d[dim][!dir]  = c.d[dim][dir] + (dir ? -1.0 : 1.0)*chimney_depth;
		c.d[dim][ dir] += (dir ? -1.0 : 1.0)*0.01*sz2; // slight shift from edge of house to avoid z-fighting
		c.z1() = c.z2();
		has_chimney = 1; // flag as interior chimney
		// no windows if there's an interior chimney that crossed the roof centerline and spans two roof tquads;
		// that case isn't handled because the clipping of the chimney to the roof is too conservative,
		// and we can't put an attic window centered on this side anyway
		float const centerline(part.get_center_dim(!dim));
		if (c.d[!dim][0] < centerline && c.d[!dim][1] > centerline) {has_attic_window = 0;}
	}
	chimney_height -= 0.4f*abs(shift); // lower it if it's not at the peak of the roof
	c.z2() += max(chimney_height, 0.75f*window_vspace); // make it at least 3/4 a story in height
	parts.push_back(c);
	max_eq(bcube.z2(), c.z2()); // update bcube to contain chimney
	return 1;
}

void building_t::maybe_gen_chimney_smoke() const {
	if (!animate2 || !begin_motion) return; // activate with 'b' key
	cube_t chimney;

	if (has_chimney == 2 && has_int_fplace) { // chimney with interior fireplace
		if (int(24534*bcube.x1()) & 3) return; // only 25% of houses have chimney smoke; use position as random seed
		chimney = get_chimney();
	}
	else if (has_smokestack) { // look for smokestack
		for (auto const &c : details) {
			if (c.type == ROOF_OBJ_SMOKESTACK) {chimney = c; break;}
		}
	}
	if (chimney.is_all_zeros()) return;
	static rand_gen_t smoke_rgen;
	if (smoke_rgen.rand_float() > 4.0f*fticks/TICKS_PER_SECOND) return; // randomly spawn every so often
	gen_arb_smoke(cube_top_center(chimney), GRAY, vector3d(0.0, 0.0, 0.75), 0.25*smoke_rgen.rand_uniform(0.015, 0.025),
		smoke_rgen.rand_uniform(0.5, 0.7), smoke_rgen.rand_uniform(0.5, 0.75), 0.0, NO_SOURCE, SMOKE, 1, 1.0, 1);
}

void building_t::gen_house(cube_t const &base, rand_gen_t &rgen) {

	assert(parts.empty());
	int const type(rgen.rand()%3); // 0=single cube, 1=L-shape, 2=two-part
	bool const two_parts(type != 0);
	unsigned force_dim[2] = {2,2}; // force roof dim to this value, per part; 2 = unforced/auto
	num_sides = 4;
	parts.reserve(two_parts ? 5 : 2); // two house sections + porch roof + porch support + chimney (upper bound)
	parts.push_back(base);
	bool const gen_door(global_building_params.windows_enabled());
	unsigned const rand_num(rgen.rand()); // bits used for random bools: {1=door_dim, 2=fence1, 4=fence2, 8=garage1, 16=garage2, 32=stacked_parts}
	float door_height(get_door_height()), door_center(0.0), door_pos(0.0), dist1(0.0), dist2(0.0);
	float const floor_spacing(get_window_vspace()), driveway_dz(0.004*door_height);
	bool const stacked_parts(!two_parts && (rand_num & 32) && bcube.dz() > 1.8*floor_spacing); // single part and at least two floors
	bool const pref_street_dim(street_dir ? ((street_dir-1) >> 1) : 0), pref_street_dir(street_dir ? ((street_dir-1)&1) : 0);
	bool door_dim(street_dir ? pref_street_dim : (rand_num & 1)), door_dir(0), dim(0), dir(0), dir2(0), skip_last_roof(0);
	unsigned door_part(0), detail_type(0), num_floors(0); // num_floors is only calculated for single cube houses and only used for multi-family houses
	real_num_parts = (two_parts ? 2 : 1); // only walkable parts: excludes shed, garage, porch roof, and chimney
	cube_t door_cube;

	if (two_parts) { // multi-part house; parts[1] is the lower height part
		dir = rgen.rand_bool(); // in dim; may be reassigned in street_dir case below
		float const split(rgen.rand_uniform(0.4, 0.6));
		float delta_height(0.0), shrink[2] = {0.0};
		parts.push_back(base); // add second part

		if (type == 1) { // L-shape
			if (street_dir) { // make sure L-type house faces the street; we can't change the garage dim if the aspect ratio is wrong, but we can at least get the dir right
				dim = rgen.rand_bool(); // need to know dim first
				((dim == pref_street_dim) ? dir : dir2) = pref_street_dir;
				((dim != pref_street_dim) ? dir : dir2) = rgen.rand_bool();
			}
			else {
				dir2 = rgen.rand_bool(); // in !dim
				dim  = rgen.rand_bool();
			}
			shrink[dir2] = rgen.rand_uniform(0.4, 0.6)*(dir2 ? -1.0 : 1.0);
			delta_height = max(0.0f, rgen.rand_uniform(-0.1, 0.5));
		}
		else if (type == 2) { // two-part
			dim          = get_largest_xy_dim(base); // choose longest dim
			delta_height = rgen.rand_uniform(0.1, 0.5);

			for (unsigned d = 0; d < 2; ++d) {
				if (rgen.rand_bool()) {shrink[d] = rgen.rand_uniform(0.2, 0.35)*(d ? -1.0 : 1.0);}
			}
			float const tot_shrink(shrink[0] - shrink[1]), un_shrink_each(0.5*(tot_shrink - 0.6)); // max shrink summed across each side is 0.6 to prevent too-small parts
			if (un_shrink_each > 0.0) {shrink[0] -= un_shrink_each; shrink[1] += un_shrink_each;}
		}
		vector2d const sz(base.get_size_xy());
		parts[0].d[ dim][ dir] += (dir ? -1.0 : 1.0)*split*sz[dim]; // split in dim
		parts[1].d[ dim][!dir]  = parts[0].d[dim][dir];
		cube_t const pre_shrunk_p1(parts[1]); // save for use in details below
		parts[1].z2() -= delta_height*parts[1].dz(); // lower height
		//ao_bcz2 = parts[1].z2(); // set this lower zval as the base AO so that AO values line up with both parts and the roof triangles
		if (global_building_params.gen_building_interiors) {adjust_part_zvals_for_floor_spacing(parts[1]);}

		if (shrink[0] == 0.0 && shrink[1] == 0.0 && parts[0].z2() == parts[1].z2()) { // same width and height - not a valid split
			bool const side(rgen.rand_bool()); // which side do we move
			shrink[side] = rgen.rand_uniform(0.2, 0.35)*(side ? -1.0 : 1.0)*sz[!dim];
		}
		for (unsigned d = 0; d < 2; ++d) {parts[1].d[!dim][d] += shrink[d]*sz[!dim];} // shrink this part in the other dim
		if (type == 1 && rgen.rand_bool()) {force_dim[0] = !dim; force_dim[1] = dim; roof_dims = 1;} // L-shape/perpendicular, half the time
		else if (type == 2) {force_dim[0] = force_dim[1] = dim; roof_dims = 2;} // two-part/parallel - force both parts to have roof along split dim
		detail_type = ((type == 1) ? (rgen.rand()%3) : 0); // 0=none, 1=porch, 2=detatched garage/shed
		door_dir    = ((door_dim == dim) ? dir : dir2); // if we have a porch/shed/garage, put the door on that side
		if (door_dim == dim && detail_type == 0) {door_dir ^= 1;} // put it on the opposite side so that the second part isn't in the way
		else if (street_dir && detail_type == 0) {door_dir = pref_street_dir;} // no porch/garage/shed
		maybe_add_basement(rgen);
		vect_cube_t cubes; // reused for tree and fence

		if (detail_type != 0) { // add details to L-shaped house
			cube_t c(pre_shrunk_p1);
			c.d[!dim][!dir2] = parts[1].d[!dim][dir2]; // other half of the shrunk part1
			float const spacing_scale((street_dir && detail_type == 2) ? 0.8 : 1.0); // larger room/smaller spacing for connecting to streets, to increase the chance of getting a garage
			dist1 = (c.d[!dim][!dir2] - base.d[!dim][dir2])*rgen.rand_uniform(0.4, 0.6)*spacing_scale;
			dist2 = (c.d[ dim][!dir ] - base.d[ dim][dir ])*rgen.rand_uniform(0.4, 0.6)*spacing_scale;
			float const base_dz(parts[1].dz()), height(min(base_dz, max(door_height/0.95f, rgen.rand_uniform(0.55, 0.7)*base_dz)));

			if (gen_door) { // add door in interior of L, centered under porch roof (if it exists, otherwise where it would be)
				door_cube   = c;
				door_center = door_cube.get_center_dim(!door_dim) + 0.5f*((door_dim == dim) ? dist1 : dist2);
				door_pos    = door_cube.d[door_dim][!door_dir];
				door_part   = ((door_dim == dim) ? 0 : 1); // which part the door is connected to
				min_eq(door_height, 0.95f*height);
			}
			if (detail_type == 1) { // porch; parts added will be {porch roof, porch support pillar}
				float const width(0.05f*(fabs(dist1) + fabs(dist2))); // width of support pillar
				c.d[!dim][dir2 ] += dist1; // move away from bcube edge
				c.d[ dim][ dir ] += dist2; // move away from bcube edge
				c.d[!dim][!dir2] -= 0.001*dist1; // adjust slightly so it's not exactly adjacent to the house and won't be considered internal face removal logic
				c.d[ dim][ !dir] -= 0.001*dist2;
				bool const add_porch_slab = 1;

				if (add_porch_slab) {
					porch = c;
					porch.z2() = porch.z1() + driveway_dz;
				}
				//c.z1() += height; // move up
				c.z1() += floor_spacing; // always one floor high
				c.z2()  = c.z1() + 0.05*parts[1].dz();
				parts.push_back(c); // porch roof
				c.z2() = c.z1();
				c.z1() = (add_porch_slab ? porch.z2() : pre_shrunk_p1.z1()); // support pillar
				c.d[!dim][!dir2] = c.d[!dim][dir2] + (dir2 ? -1.0 : 1.0)*width;
				c.d[ dim][!dir ] = c.d[ dim][ dir] + (dir  ? -1.0 : 1.0)*width;
				skip_last_roof = 1;
			}
			else if (detail_type == 2) { // detached garage/shed
				c.d[!dim][dir2 ]  = base.d[!dim][dir2]; // shove it into the opposite corner of the bcube
				c.d[ dim][dir  ]  = base.d[ dim][dir ]; // shove it into the opposite corner of the bcube
				c.d[!dim][!dir2] -= dist1; // move away from bcube edge
				c.d[ dim][!dir ] -= dist2; // move away from bcube edge
				c.z2() = c.z1() + max(floor_spacing, min(min(c.dx(), c.dy()), height)); // no taller than x or y size, but at least one floor high; Note: z1 same as part1
				bool is_garage(car_can_fit(c)), pri_dim(c.dx() < c.dy()); // garage must be able to fit a car
				// maybe we could extend the garage in the correct dim to fit a car, but that can lead to problems with the fence and driveway placement
				//if (street_dir && pri_dim != pref_street_dim) {is_garage = 0;} // garage is in wrong dim, make it a shed instead (we're allowing driveway bends, so this is okay now)
				(is_garage ? has_garage : has_shed) = 1;
				// add a door
				bool gdoor_dir(c.d[pri_dim][0] == bcube.d[pri_dim][0]); // interior face
				float wscale(DOOR_WIDTH_SCALE);
				float const gs_door_height(min(door_height, 0.9f*c.dz())); // clamp to 90% of garage/shed height to make sure there's room for the ceiling trim
				
				if (is_garage) {
					float const shelf_depth(0.15*floor_spacing), garage_width(c.get_sz_dim(!pri_dim));
					wscale     = min(0.9f*garage_width, (garage_width - 2.0f*shelf_depth))/gs_door_height; // avoid clipping through shelves
					gdoor_dir ^= 1; // facing away from the house, not toward it
				}
				add_door(place_door(c, pri_dim, gdoor_dir, gs_door_height, 0.0, 0.0, 0.0, wscale, 0, is_garage, rgen), parts.size(), pri_dim, gdoor_dir, 0);
				if (is_garage) {doors.back().type = tquad_with_ix_t::TYPE_GDOOR;} // make it a garage door rather than a house door

				if (is_garage) { // add driveway
					bool const ddir((pri_dim == dim) ? dir : dir2), has_r90_turn(street_dir && (2*pri_dim + ddir + 1) != street_dir);
					float const length(c.get_sz_dim(pri_dim)), width(c.get_sz_dim(!pri_dim));
					driveway = c;
					driveway.z2() = driveway.z1() + driveway_dz;
					driveway.expand_in_dim(!pri_dim, -0.05*width); // shrink slightly to a bit narrower than the garage
					driveway.d[pri_dim][!ddir] = driveway.d[pri_dim][ddir];
					driveway.d[pri_dim][ ddir] += (ddir ? 1.0 : -1.0)*max(0.4*length, (has_r90_turn ? 1.2 : 0.9)*width); // set length extension, wider for 90 degree turn
					if (!assigned_plot.is_all_zeros()) {driveway.intersect_with_cube_xy(assigned_plot);} // clip driveway to the assigned plot, if one was specified
					//garage_dim = pri_dim;
				}
			}
			parts.push_back(c); // support column or shed/garage
		} // end house details
		else if (type == 1 && has_city_trees() && !is_rotated()) { // place a tree for L-shaped house with no garage or shed
			cube_t empty_space(bcube);

			for (unsigned e = 0; e < 2; ++e) {
				subtract_cube_from_cube(empty_space, parts[e], cubes, 1); // clear_out=1
				assert(cubes.size() == 1); // Note: may fail for rotated buildings
				empty_space = cubes[0];
			}
			tree_pos = empty_space.get_cube_center(); // centered on where the garage/shed would have been
			vector3d const tdir(tree_pos - bcube.get_cube_center());
			tree_pos  += (0.05f*(bcube.dx() + bcube.dy())/tdir.mag())*tdir; // shift slightly away from house center so that it's less likely to intersect the house
			tree_pos.z = ground_floor_z1;
		}
		if (type == 1 && (rand_num & 6) != 0 && !is_rotated()) { // L-shaped house, add a fence 75% of the time (skip rotated due to assert in clip_cube_to_parts())
			// if stree_dim is set, make primary fence perpendicular; otherwise use same dim as shrunk house part for full fence section; or should we use garage_dim?
			bool const fdim(street_dir ? pref_street_dim : dim), fdir((fdim == dim) ? dir2 : dir);
			float const fence_thickness(0.08*door_height);
			cube_t fence(bcube), fence2(bcube);
			fence.z1() = fence2.z1() = ground_floor_z1;
			fence.z2() = fence2.z2() = ground_floor_z1 + 0.65*door_height; // set fence height
			fence.d[!fdim][!fdir] = fence.d[!fdim][fdir] + (fdir ? -1.0 : 1.0)*fence_thickness;
			// we can calculate the exact location of the fence, but it depends on detail_type, garage/shed position, etc.,
			// so it's easier to start with the entire edge and clip it by the house parts and optional garage/shed
			
			if (clip_cube_to_parts(fence, cubes)) {
				fence.expand_in_dim(fdim, -0.01*door_height); // shrink slightly to avoid clipping through the exterior wall
				fences.push_back(fence);

				if ((rand_num & 6) != 2) { // 67% of the time
					bool const fdir2((fdim == dim) ? dir : dir2);
					fence2.d[fdim][!fdir2] = fence2.d[fdim][fdir2] + (fdir2 ? -1.0 : 1.0)*fence_thickness;
					
					if (clip_cube_to_parts(fence2, cubes)) {
						float const center_pt(fence2.get_center_dim(!fdim)), gate_width(1.0*door_height);

						for (unsigned d = 0; d < 2; ++d) { // split fence, add gap for gate, and add each segment if it's large enough
							cube_t seg(fence2);
							seg.d[!fdim][d] = center_pt + (d ? -0.5f : 0.5f)*gate_width;
							if (seg.get_sz_dim(!fdim) < gate_width) continue; // too small
							seg.expand_in_dim(!fdim, -0.01*door_height); // shrink slightly to avoid clipping through the exterior wall
							fences.push_back(seg);
						}
					}
				}
			}
		}
		if (force_dim[1] == 2 && parts[0].z2() == parts[1].z2()) { // L-shaped house with both parts at the same height
			bool const pref_dim((force_dim[0] < 2) ? (1-force_dim[0]) : (!get_largest_xy_dim(parts[0])));

			if (parts[1].get_sz_dim(!pref_dim) < 1.2*parts[0].get_sz_dim(pref_dim)) { // not too high an aspect ratio between the widths of the two parts
				force_dim[1] = pref_dim; // force roof orients to be perpendicular to avoid an odd crease in the middle
				roof_dims    = 1; // mark as perpendicular
			}
		}
	} // end two_parts (multi-part house)
	else if (stacked_parts) {
		cube_t top_part(base); // parts[0] is the bottom
		parts[0].z2() = rgen.rand_uniform((base.z1() + 0.8*floor_spacing), (base.z2() - 0.8*floor_spacing)); // split zval
		if (global_building_params.gen_building_interiors) {adjust_part_zvals_for_floor_spacing(parts[0]);}
		top_part.z1() = parts[0].z2();
		if (global_building_params.gen_building_interiors) {adjust_part_zvals_for_floor_spacing(top_part);}
		float const sz[2] = {base.dx(), base.dy()};

		for (unsigned dim = 0; dim < 2; ++dim) { // x/y
			for (unsigned dir = 0; dir < 2; ++dir) {
				if (!rgen.rand_bool()) continue; // no shrink on this side
				top_part.d[dim][dir] += (dir ? -1.0 : 1.0)*sz[dim]*rgen.rand_uniform(0.1, 0.25); // shrink edge
			}
		} // for dim
		parts.push_back(top_part);
		++real_num_parts;
		maybe_add_basement(rgen);
	}
	else { // single cube house
		num_floors = calc_num_floors(parts[0], floor_spacing, get_floor_thickness());
		// make it a multi-family house if it's a single large part with at least three floors
		multi_family = (num_floors > 2 && parts[0].dx()*parts[0].dy() > 50.0*floor_spacing*floor_spacing);
		if (multi_family) {btype = BTYPE_MULT_FAM;}
		maybe_add_basement(rgen);

		if (gen_door) { // have exterior doors and windows
			door_dir  = (street_dir ? pref_street_dir : rgen.rand_bool()); // select a random dir if street_dir is not set
			door_part = 0; // only one part
		}
	}
	calc_bcube_from_parts(); // maybe calculate a tighter bounding cube
	gen_interior(rgen, 0); // before adding door

	if (gen_door) { // add exterior doors and possibly a garage + driveway and extended basement
		// attempt to add an interior garage when legal, always when along a street, else 75% of the time; not for multi-family, since we can't make them one per resident
		if (!has_garage && !multi_family && (street_dir || (rand_num & 24))) {
			bool gdim(0), gdir(0);

			if (maybe_assign_interior_garage(gdim, gdir)) { // assigned a garage
				has_int_garage = 1;
				room_t const &garage(interior->get_garage_room());
				float const wscale(0.9*garage.get_sz_dim(!gdim)/door_height);
				add_door(place_door(garage, gdim, gdir, door_height, 0.0, 0.0, 0.0, wscale, 0, 1, rgen), garage.part_id, gdim, gdir, 0); // centered in garage
				doors.back().type = tquad_with_ix_t::TYPE_GDOOR; // make it a garage door
				driveway = garage;
				driveway.z2() = garage.z1() + driveway_dz;
				driveway.d[gdim][!gdir] = garage.d[gdim][gdir]; // start pos
				driveway.d[gdim][ gdir] = garage.d[gdim][gdir] + (gdir ? 1.0 : -1.0)*0.4*garage.get_sz_dim(gdim); // extend outward
				if (!assigned_plot.is_all_zeros()) {driveway.intersect_with_cube_xy(assigned_plot);} // clip driveway to the assigned plot, if one was specified
			}
		}
		cube_t door;
		bool door_valid(0);
		unsigned const num_tries = 10;
		
		for (unsigned n = 0; n < num_tries; ++n) { // make several attempts at generating a valid interior where this door can be placed; this still fails a few times
			for (unsigned d = 0; d < 2; ++d) { // try both dims
				door = place_door(parts[door_part], door_dim, door_dir, door_height, door_center, door_pos, 0.25, DOOR_WIDTH_SCALE, 0, 0, rgen); // keep door_height
				if (is_valid_door_pos(door, 0.5*door_height, door_dim)) {door_valid = 1; break;} // done
				if (d == 1 && n+1 == num_tries) break; // no more tries available
				// swap door_dim and calculate new door_dir; this is duplicated from the code above, but not easier to factor out
				door_dim ^= 1;
				door_dir  = (two_parts ? ((door_dim == dim) ? dir : dir2) : (door_dir^1)); // if we have a porch/shed/garage, put the door on that side
				if (door_dim == dim && two_parts && detail_type == 0) {door_dir ^= 1;} // put it on the opposite side so that the second part isn't in the way

				if (door_center != 0.0) { // reclaculate for L-shaped house
					door_center = door_cube.get_center_dim(!door_dim) + 0.5f*((door_dim == dim) ? dist1 : dist2);
					door_pos    = door_cube.d[door_dim][!door_dir];
					door_part   = ((door_dim == dim) ? 0 : 1); // which part the door is connected to
				}
			} // for d
			if (door_valid) break;

			if (n+1 < num_tries) { // still invalid, regenerate interior
				if (has_int_garage) { // must also remove garage, garage door, and driveway
					assert(doors.size() == 1);
					driveway.set_to_zeros();
					doors.pop_back();
					has_int_garage = 0;
					interior->garage_room = -1;
				}
				gen_interior(rgen, 0);
			}
		} // for n
		if (/*door_valid &&*/ add_door(door, door_part, door_dim, door_dir, 0)) {floor_ext_door_mask |= 1;} // should we still add a door if it's invalid?
		if (doors.size() == 2) {swap(doors[0], doors[1]);} // make sure the house door comes before the garage/shed door
		float const tot_area(parts[0].get_area_xy() + (two_parts ? parts[1].get_area_xy() : 0.0f));

		if (tot_area > 25.0f*floor_spacing*floor_spacing) { // if house is large enough, add a back door
			bool added_door(0);

			for (unsigned p = 0; p < real_num_parts && !added_door; ++p) {
				unsigned door2_part(two_parts ? 1-door_part : 0); // try the other part first
				cube_t const &part(parts[door2_part]);

				for (unsigned d = 0; d < 2 && !added_door; ++d) {
					bool const door2_dim(door_dim ^ bool(d)); // try same dim first

					for (unsigned door2_dir = 0; door2_dir < 2 && !added_door; ++door2_dir) {
						if (door2_dim == door_dim && door_dir == bool(door2_dir) /*&& door2_part == door_part*/) continue; // don't place second door on the same side
						if (part.d[door2_dim][door2_dir] != bcube.d[door2_dim][door2_dir]) continue; // door on building bcube is always exterior/never interior face between two parts
						cube_t const door2(place_door(part, door2_dim, door2_dir, door_height, 0.0, 0.0, 0.25, DOOR_WIDTH_SCALE, 1, 0, rgen));
						if (!is_valid_door_pos(door2, 0.5*door_height, door2_dim)) continue; // bad placement
						added_door |= add_door(door2, door2_part, door2_dim, door2_dir, 0);
					} // for door_dir2
				} // for d
			} // for p
			if (added_door) {floor_ext_door_mask |= 1;}
		} // end back door
		if (multi_family) { // maybe place upper floor door(s); house should be a single cube
			assert(num_floors > 1);

			for (unsigned f = 1; f < num_floors; ++f) { // every floor above ground level
				bool const dim0(rgen.rand_bool()), dir0(rgen.rand_bool());
				bool added_door(0);

				for (unsigned d = 0; d < 2 && !added_door; ++d) {
					for (unsigned e = 0; e < 2 && !added_door; ++e) {
						unsigned const dim(dim0 ^ bool(d)), dir(dir0 ^ bool(e));
						cube_t const door2(place_door(parts[0], dim, dir, door_height, 0.0, 0.0, 0.25, DOOR_WIDTH_SCALE, 1, 0, rgen, f));
						if (is_valid_door_pos(door2, 0.5*door_height, dim)) {added_door |= add_door(door2, 0, dim, dir, 0);}
					}
				}
				if (added_door) {floor_ext_door_mask |= (1 << f);}
			} // for f
		}
		if (global_building_params.max_ext_basement_room_depth > 0) {
			extend_underground_basement(rgen); // maybe add door inside basement and connected extended basement; rgen is copied, not modified
		}
	} // end gen_door
	if (interior) {interior->assign_door_conn_rooms();} // must be after adding extended basement and before creating two story tall rooms
	create_two_story_tall_rooms(rgen); // must be after generating interior and adding exterior doors
	// add roof tquads
	float peak_height(rgen.rand_uniform(0.15, 0.5)); // same for all parts
	if (has_attic()) {max_eq(peak_height, 0.3f);} // set larger min size if there's an attic
	float roof_dz[4] = {};
	bool any_hipped(0), hipped_roof[4] = {};
	int last_hipped(2); // starts at <unset>

	for (auto i = parts.begin(); (i + skip_last_roof) != parts.end(); ++i) {
		if (is_basement(i)) continue; // skip the basement
		unsigned const ix(i - parts.begin());

		if (stacked_parts && ix == 0) { // bottom part of stacked parts
			roof_dz[ix] = gen_sloped_roof_for_stacked_parts(parts[0], parts[1]);
			continue;
		}
		bool const main_part(ix < real_num_parts);
		unsigned const fdim(main_part ? force_dim[ix] : 2);
		cube_t const &other((two_parts && main_part) ? parts[1-ix] : *i); // == self for single part houses
		bool const dim((fdim < 2) ? fdim : get_largest_xy_dim(*i)); // use longest side if not forced
		bool const other_dim(two_parts ? ((main_part && force_dim[1-ix] < 2) ? force_dim[1-ix] : get_largest_xy_dim(other)) : 0);
		float extend_to(0.0), max_dz(i->dz());

		if (has_sec_bldg() && ix == real_num_parts) { // garage or shed
			if (has_shed) {} // use ROOF_TYPE_SHED? need to add triangle wall sections for this; currently we can't have a per-part roof type either
		}
		if (type == 1 && ix == 1 && dim != other_dim && parts[0].z2() == parts[1].z2()) { // L-shape, same z2, opposite dim T-junction
			max_dz    = peak_height*parts[0].get_sz_dim(!other_dim); // clamp roof zval to other roof's peak
			extend_to = parts[0].get_center_dim(!other_dim); // extend lower part roof to center of upper part roof
		}
		bool can_be_hipped(main_part && extend_to == 0.0 && i->get_sz_dim(dim) > i->get_sz_dim(!dim)); // must be longer dim
		
		if (can_be_hipped && two_parts) {
			if (dim == other_dim && (i->d[!dim][0] == other.d[!dim][1] || i->d[!dim][1] == other.d[!dim][0])) {} // two parallel adjacent parts - hipped roof is allowed
			else {
				float const part_roof_z (i->z2()    + min(peak_height*i->get_sz_dim(!dim), i->dz()));
				float const other_roof_z(other.z2() + min(peak_height*other.get_sz_dim(!other_dim), other.dz()));
				can_be_hipped = (part_roof_z >= other_roof_z); // no hipped for lower part
			}
		}
		bool const hipped(can_be_hipped && ((last_hipped == 2) ? rgen.rand_bool() : last_hipped)); // hipped roof 50% of the time, with preference to share across parts
		last_hipped = hipped;
		if (hipped) {roof_dz[ix] = gen_hipped_roof(*i, peak_height, extend_to);}
		else {
			unsigned skip_side_tri(2); // default = skip neither
			
			if (two_parts && main_part && (dim == other_dim || hipped_roof[1-ix]) && i->d[!dim][0] >= other.d[!dim][0] && i->d[!dim][1] <= other.d[!dim][1] && i->z2() <= other.z2()) {
				// side of i contained in other
				if (ix == 1 && i->z2() == other.z2()) { // same height - check if we need to ajust the slope of this roof to match
					float const roof_slope(roof_dz[0]/other.get_sz_dim(!dim));
					min_eq(max_dz, roof_slope*i->get_sz_dim(!dim)); // peak_height cancels out
				}
				for (unsigned d = 0; d < 2; ++d) {
					if (dim == other_dim && i->d[dim][d] == other.d[dim][!d]) {skip_side_tri = d;} // remove smaller of two opposing/overlapping triangles to prevent z-fighting
				}
			}
			assert(max_dz > 0.0);
			roof_dz[ix] = gen_peaked_roof(*i, peak_height, dim, extend_to, max_dz, skip_side_tri);
		}
		assert(roof_dz[ix] > 0.0);
		hipped_roof[ix] = hipped;
		any_hipped     |= hipped;
	} // for i
	if ((rgen.rand()%3) != 0) { // add a chimney 67% of the time
		add_chimney(two_parts, stacked_parts, hipped_roof, roof_dz, force_dim, rgen); // Note: return value is ignored
	}
	// Note: driveway collisions are handled through check_road_seg_sphere_coll()
	parts_generated = 1; // must be after adding chimney
	roof_type = (any_hipped ? ROOF_TYPE_HIPPED : ROOF_TYPE_PEAK);
	add_roof_to_bcube();
	gen_grayscale_detail_color(rgen, 0.4, 0.8); // for roof
	door_color = (rgen.rand_bool() ? LT_BROWN : WHITE);
	// white, white, white, white, pink, peach, lt green, lt blue
	colorRGBA const wall_colors[8] = {WHITE, WHITE, WHITE, WHITE, colorRGBA(1.0, 0.85, 0.85), colorRGBA(1.0, 0.85, 0.75), colorRGBA(0.85, 1.0, 0.85), colorRGBA(0.85, 0.85, 1.0)};
	wall_color = wall_color.modulate_with(wall_colors[rgen.rand()%8]);
	unsigned const roof_obj_type((rgen.rand() & 7));
	if (roof_obj_type &  1) {add_solar_panels(rgen);} // maybe add solar panels - 50%
	if (roof_obj_type == 0) {add_sat_dish    (rgen);} // add a satellite dish - 1/8
	if (roof_obj_type == 2) {add_tv_antenna  (rgen);} // add a TV antenna - 1/8
	if (rgen.rand_bool()) {add_outdoor_ac_unit(rgen);} // place an outdoor AC unit against an exterior wall 50% of the time, not actually on the roof
	if (has_basement()) {has_basement_pipes = rgen.rand_bool();}
	if (interior) {interior->finalize();}
	setup_damage_vals();
}

bool building_t::add_outdoor_ac_unit(rand_gen_t &rgen) { // for houses
	float const door_height(get_door_height());
	float const depth(door_height*rgen.rand_uniform(0.26, 0.35)), width(1.5*depth), height(door_height*rgen.rand_uniform(0.32, 0.36)); // constant width to keep the texture square
	unsigned const ac_part_ix((real_num_parts == 2 && parts[0].z1() == parts[1].z1()) ? rgen.rand_bool() : 0); // select a random part if two are on the ground floor
	cube_t const &ac_part(parts[ac_part_ix]);
	bool const ac_dim(rgen.rand_bool()), ac_dir(rgen.rand_bool());
	if (ac_part.get_sz_dim(!ac_dim) < 4.0*width) return 0; // not enough space
	float const place_pos(rgen.rand_uniform(ac_part.d[!ac_dim][0]+width, ac_part.d[!ac_dim][1]-width));
	roof_obj_t ac(ROOF_OBJ_AC);
	ac.d[ac_dim][!ac_dir] = ac_part.d[ac_dim][ ac_dir] + (ac_dir ? 1.0 : -1.0)*0.07*get_window_vspace(); // place slightly away from the exterior wall to avoid window sills
	ac.d[ac_dim][ ac_dir] = ac     .d[ac_dim][!ac_dir] + (ac_dir ? 1.0 : -1.0)*depth;
	set_wall_width(ac, place_pos, 0.5*width, !ac_dim);
	set_cube_zvals(ac, ac_part.z1(), (ac_part.z1() + height));
	if (has_driveway() && ac.intersects_xy(driveway)) return 0; // driveway in the way
	if (has_porch   () && ac.intersects_xy(porch   )) return 0; // porch    in the way

	for (unsigned i = 0; i < parts.size(); ++i) {
		if (i != ac_part_ix && ac.intersects(parts[i])) return 0; // intersects another part, or the fireplace
	}
	if (is_cube_close_to_exterior_doorway(ac, width, 1)) return 0;
	details.push_back(ac);
	union_with_coll_bcube(ac);
	has_ac = 1;
	return 1;
}

bool is_valid_driveway_pos(cube_t const &driveway, cube_t const &bcube, vect_cube_t const &bcubes) {
	for (auto i = bcubes.begin(); i != bcubes.end(); ++i) {
		if (*i != bcube && i->intersects_xy(driveway)) return 0;
	}
	return 1;
}
bool add_driveway_if_legal(cube_t &dw, cube_t const &target, cube_t const &bcube, cube_t const &sub_plot, vect_cube_t const &avoid,
	vect_cube_t const &blockers, vect_cube_t const &bcubes, cube_t &ret, rand_gen_t &rgen, float hwidth, bool dim, bool dir)
{
	if (target.get_sz_dim(!dim) < 3.0*hwidth) return 0; // not wide enough for driveway
	dw.d[dim][!dir] = target.d[dim][dir];
	if (dw.get_sz_dim(dim) < 3.0*hwidth) return 0; // too short compared to the width, likely can't fit a car

	for (unsigned n = 0; n < 10; ++n) { // make up to 10 attempts
		float const pos(rgen.rand_uniform((target.d[!dim][0] + hwidth), (target.d[!dim][1] - hwidth)));
		set_wall_width(dw, pos, hwidth, !dim);
		if (has_bcube_int(dw, avoid))                  continue; // blocked, including adjacency
		if (has_bcube_int_xy_no_adj(dw, blockers))     continue; // blocked
		if (!is_valid_driveway_pos(dw, bcube, bcubes)) continue; // blocked
		if (!sub_plot.contains_cube_xy(dw))            continue; // extends outside the plot
		assert(dw.dx() > 0 && dw.dy() > 0); // strictly normalized in XY
		ret = dw;
		return 1;
	}
	return 0; // failed
}
bool extend_existing_driveway_in_dim(cube_t const &driveway, cube_t const &plot, cube_t const &bcube, vect_cube_t const &bcubes, cube_t &ret, bool dim, bool dir) {
	cube_t dw(driveway);
	dw.d[dim][ dir] = plot    .d[dim][dir];
	dw.d[dim][!dir] = driveway.d[dim][dir]; // shorten to connect to existing house driveway
	float const len(dw.get_sz_dim(dim));
	if (len < 0.01*dw.get_sz_dim(!dim)) return 1; // length is tiny, extension is not needed
	assert(dw.is_strictly_normalized());
	if (len > 0.3*plot.get_sz_dim(dim)) return 0; // don't let the driveway cross through the entire block on the wrong side of the house
	if (!is_valid_driveway_pos(dw, bcube, bcubes)) return 0;
	dw.z1() = dw.z2(); // shrink to zero height
	ret = dw;
	return 1; // success
}
bool extend_existing_driveway(cube_t const &driveway, cube_t const &plot, cube_t const &bcube, vect_cube_t const &bcubes, cube_t &ret) {
	if (!plot.contains_cube_xy(driveway)) return 1; // driveway already extends to plot, don't need to add an extension, mark as done
	bool dim(0), dir(0);
	if      (driveway.x1() < bcube.x1()) {dim = 0; dir = 0;}
	else if (driveway.x2() > bcube.x2()) {dim = 0; dir = 1;}
	else if (driveway.y1() < bcube.y1()) {dim = 1; dir = 0;}
	else if (driveway.y2() > bcube.y2()) {dim = 1; dir = 1;}
	else {cout << TXT(driveway.str()) << TXT(bcube.str()) << endl; assert(0); return 0;} // not adjacent?
	if (extend_existing_driveway_in_dim(driveway, plot, bcube, bcubes, ret, dim, dir)) return 1;
	// what about making a right angle turn?
	dim ^= 1; // try the other dim
	dir  = ((plot.d[dim][1] - driveway.d[dim][1]) < (driveway.d[dim][0] - plot.d[dim][0])); // closer dir
	return extend_existing_driveway_in_dim(driveway, plot, bcube, bcubes, ret, dim, dir);
}

void get_closest_dim_dir_xy(cube_t const &inner, cube_t const &outer, bool &dim, bool &dir) {
	float dmin(0.0);
	bool dmin_set(0);

	for (unsigned d = 0; d < 2; ++d) { // find closest edge of *i to edge of bcube
		for (unsigned e = 0; e < 2; ++e) {
			float const dist(fabs(inner.d[d][e] - outer.d[d][e]));
			if (!dmin_set || dist < dmin) {dim = d; dir = e; dmin = dist; dmin_set = 1;}
		}
	}
}
unsigned get_street_dir(cube_t const &inner, cube_t const &outer) {
	bool dim(0), dir(0);
	get_closest_dim_dir_xy(inner, outer, dim, dir);
	return 2*dim + dir + 1;
}

// Note: applies to city residential plots, called by city_obj_placer_t;
// shouldn't be const because this modifies city_driveway as it's needed for balcony placement logic later; but this must be const for city API, so city_driveway is mutable
bool building_t::maybe_add_house_driveway(cube_t const &plot, unsigned building_ix) const {
	if (!is_house) return 0;
	assert(plot.contains_cube_xy(bcube));
	cube_t sub_plot(plot);

	if (!assigned_plot.is_all_zeros()) { // check for residential sub-plot (could also be the same as plot)
		assert(plot.contains_cube_xy(assigned_plot));
		assert(assigned_plot.contains_cube_xy(bcube));
		sub_plot = assigned_plot;
	}
	static vect_cube_t bcubes, avoid; // reused across calls
	bcubes.clear();
	get_building_bcubes(plot, bcubes); // Note: could be passed in by the caller, but requires many changes
	// garages should be placed so that their driveways exit toward the edge of the plot
	if (has_driveway() && extend_existing_driveway(driveway, plot, bcube, bcubes, city_driveway)) return 1;
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, building_ix+111);
	rgen.rand_mix();
	float const hwidth(0.8*get_window_vspace()*rgen.rand_uniform(0.9, 1.1));
	avoid = fences;
	if (has_porch()) {avoid.push_back(porch);}
	if (has_chimney == 2) {avoid.push_back(get_fireplace());} // avoid placing the driveway across from the chimney/fireplace
	if (tree_pos != all_zeros) {avoid.push_back(cube_t()); avoid.back().set_from_sphere(tree_pos, hwidth);} // assume tree diameter is similar to driveway width
	cube_t ac_unit;
	
	if (has_ac) { // include outdoor AC unit; would be better if we can move the AC instead, but it's too late for that
		assert(!details.empty() && details.back().type == ROOF_OBJ_AC);
		ac_unit = details.back();
		avoid.push_back(ac_unit);
	}
	bool dim(0), dir(0);
	get_closest_dim_dir_xy(bcube, plot, dim, dir); // must use larger plot, not sub_plot

	for (unsigned n = 0; n < 2; ++n) { // check the closest dim, then the second closest dim
		cube_t dw(plot); // copy zvals and dim/dir from plot

		for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
			if (add_driveway_if_legal(dw, *i, bcube, sub_plot, avoid, parts, bcubes, city_driveway, rgen, hwidth, dim, dir)) return 1;
		}
		if (add_driveway_if_legal(dw, bcube, bcube, sub_plot, avoid, vect_cube_t(), bcubes, city_driveway, rgen, hwidth, dim, dir)) return 1;
		// maybe it was too short? try to place to one side of the house or the other
		dw.d[dim][!dir] = bcube.d[dim][!dir]; // extend up to far side of house
		float const length(dw.get_sz_dim(dim)), min_length(3.0*hwidth);
		assert(length > 0.0);
		bool const first_side(rgen.rand_bool()); // not biased

		if (length > min_length) { // shorten a random amount if it's long enough
			dw.d[dim][!dir] -= (dir ? -1.0 : 1.0)*min((length - min_length), hwidth*rgen.rand_uniform(0.0, 4.0));
		}
		for (unsigned S = 0; S < 2; ++S) {
			bool const s(bool(S) ^ first_side);
			set_wall_width(dw, (bcube.d[!dim][s] + (s ? 1.0 : -1.0)*hwidth), hwidth, !dim);
			if (!sub_plot.contains_cube_xy(dw)) continue; // extends outside the plot
			if (!is_valid_driveway_pos(dw, bcube, bcubes)) continue; // blocked (don't need to check parts or chimney/fireplace here)
			float const min_plot_edge_spacing(0.5*get_nom_car_size().x); // half a car length
			// check if too close to edge of plot/too close to intersection
			if (dw.d[!dim][0] < plot.d[!dim][0]+min_plot_edge_spacing || dw.d[!dim][1] > plot.d[!dim][1]-min_plot_edge_spacing) continue;
			if (has_ac && dw.intersects_xy(ac_unit)) continue;
			if (s) {dw.translate_dim(2, 0.001*hwidth);} // hack to prevent z-fighting when driveways overlap on the left and right of adjacent houses
			assert(dw.dx() > 0 && dw.dy() > 0); // strictly normalized in XY
			city_driveway = dw;
			return 1;
		} // for s
		if (n == 0) {
			dim ^= 1; // try the other dim
			dir  = (bcube.get_center_dim(dim) > plot.get_center_dim(dim)); // choose the closer dir
		}
	} // for n
	return 0; // failed to add driveway; this is rare, but can happen with a house that's close to the road, close to the plot edge on one side, and has an AC unit on the other side
}

void try_expand_into_xy(cube_t &c1, cube_t const &c2) {
	for (unsigned d = 0; d < 2; ++d) { // attempt to merge in this dim
		if (c2.d[d][0] > c1.d[d][0] || c2.d[d][1] < c1.d[d][1]) continue; // edge not contained

		if (c2.d[!d][0] <= c1.d[!d][1] && c1.d[!d][0] <= c2.d[!d][1]) { // overlap - take the union in this dim
			min_eq(c1.d[!d][0], c2.d[!d][0]);
			max_eq(c1.d[!d][1], c2.d[!d][1]);
		}
	} // for d
}

bool building_t::can_have_basement() const {
	if (!is_cube())          return 0; // simple cube shaped buildings only
	if (!interior_enabled()) return 0; // if there's no interior, there's no point in adding a basement
	return 1;
}

// for houses or office buildings
void building_t::maybe_add_basement(rand_gen_t rgen) { // rgen passed by value so that the original isn't modified
	if (!can_have_basement()) return;
	float basement_prob(is_house ? global_building_params.basement_prob_house : global_building_params.basement_prob_office);
	if (is_in_city && !is_house) {basement_prob *= 2.0;} // double the basement/parking garage probability for city office buildings
	if (!rgen.rand_probability(basement_prob)) return; // no basement
	float const floor_spacing(get_window_vspace()), max_sea_level(get_max_sea_level());
	float basement_z1(ground_floor_z1 - floor_spacing);
	if (basement_z1 < max_sea_level) return; // no basement below the water line
	basement_part_ix = (int8_t)parts.size(); // index of part that will be added below
	cube_t basement;
	
	if (is_house) {
		basement = parts[0]; // start at first part

		if (real_num_parts == 2) {
			unsigned const ix(parts[1].get_area_xy() > parts[0].get_area_xy());
			basement = parts[ix]; // start with the larger part
			// attempt to expand into the smaller part as long as it fits within the footprint of the upper floors
			try_expand_into_xy(basement, parts[!ix]);
		}
	}
	else { // office building
		if (global_building_params.max_office_basement_floors == 0) return; // disabled; but user should set basement_prob=0 rather than floor=0
		unsigned num_gf_parts(0);
		cube_t largest_gf_part;

		for (auto const &p : parts) {
			if (p.z1() != ground_floor_z1) continue; // only count ground floor parts
			++num_gf_parts;
			if (p.get_area_xy() > largest_gf_part.get_area_xy()) {largest_gf_part = p;}
		}
		assert(num_gf_parts > 0);
		// Note: player collision is still based on ground floor building shape, so having the basement extend outside of this footprint still doesn't work
		//if (num_gf_parts < 4) {basement = bcube;} // not an O-shaped building, don't need to handle tree in the middle, use full bcube
		// use largest ground floor part, then try to expand
		basement = largest_gf_part;
		real_num_parts = (uint8_t)parts.size(); // set now because it's needed in the call below
		expand_ground_floor_cube(basement);
		// maybe extend the basement downward with extra floors
		unsigned const num_basement_floors(1 + (rgen.rand() % global_building_params.max_office_basement_floors)); // 1+
		for (unsigned n = 1; n < num_basement_floors && (basement_z1 - floor_spacing) > max_sea_level; ++n) {basement_z1 -= floor_spacing;}
	}
	set_cube_zvals(basement, basement_z1, ground_floor_z1);
	parts.push_back(basement);
	min_eq(bcube.z1(), basement_z1); // not really necessary, will be updated later anyway, but good to have here for reference; orig bcube.z1() is saved in ground_floor_z1
	coll_bcube.union_with_cube(basement);
	++real_num_parts;
}

int building_t::find_main_roof_tquad_ix(rand_gen_t &rgen, bool skip_if_has_other_obj) const {
	assert(!roof_tquads.empty());
	int best_tquad(-1);
	float best_zmax(0.0), best_area(0.0);

	for (auto r = roof_tquads.begin(); r != roof_tquads.end(); ++r) { // find highest roof, otherwise largest area if tied
		if (!r->is_roof()) continue; // only include roofs
		if (r->npts != 4)  continue; // only quads, skip triangles
		float const zmax(max(max(r->pts[0].z, r->pts[1].z), max(r->pts[2].z, r->pts[3].z))), area(polygon_area(r->pts, r->npts));

		// determine the best candidate so far; use a random value to break ties
		if (best_area == 0.0 || zmax > best_zmax || (zmax == best_zmax && area > best_area) || (zmax == best_zmax && area == best_area && rgen.rand_bool())) {
			if (skip_if_has_other_obj) {
				cube_t bcube(r->get_bcube());
				bcube.expand_by(-0.1*bcube.get_size()); // shrink by 10% (upper bound of approx solar panel area)
				bool valid(1);

				for (auto r2 = roof_tquads.begin(); r2 != roof_tquads.end(); ++r2) { // check other objects, including anything on the roof that's not part of the roof
					if (r2 == r) continue; // skip self
					if (r2->type == tquad_with_ix_t::TYPE_ROOF_HIP && r2->npts == 3) continue; // skip hipped roof triangles as they won't intersect the adjacent quad
					if (r2->get_bcube().intersects(bcube)) {valid = 0; break;}
				}
				if (!valid) continue; // intersects another roof section, skip
			}
			best_zmax  = zmax;
			best_area  = area;
			best_tquad = int(r - roof_tquads.begin());
		}
	} // for r
	return best_tquad;
}

// Note: occasionally the chosen point will generate a wire that intersects some other part of the house for every nearby power pole and may be skipped
bool building_t::get_power_point(vector<point> &ppts) const {
	if (!is_house) return 0; // houses only for now
	static rand_gen_t rgen; // used for tie breaker when both sides of the roof are symmetric
	int const roof_tquad_ix(find_main_roof_tquad_ix(rgen, 0)); // skip_if_has_other_obj=0
	if (roof_tquad_ix < 0) return 0; // not found
	assert((unsigned)roof_tquad_ix < roof_tquads.size());
	tquad_with_ix_t const &roof(roof_tquads[roof_tquad_ix]);
	assert(roof.npts == 4);
	float zmax(0.0);
	unsigned ridge_ix(0);

	for (unsigned n = 0; n < 4; ++n) {
		if (n == 0 || roof.pts[n].z > zmax) {zmax = roof.pts[n].z; ridge_ix = n;}
	}
	assert(ridge_ix < 3);
	unsigned ridge_ix2(ridge_ix + 1); // assume next point
	if (ridge_ix == 0 && roof.pts[ridge_ix2].z != zmax) {ridge_ix2 = 3;} // wraps around from 3 => 0
	point const &p1(roof.pts[ridge_ix]), &p2(roof.pts[ridge_ix2]);
	assert(p2.z == zmax); // there must be a ridge of two adjacent points at max height
	float ridge_pos(0.0);

	if (street_dir) {
		bool const street_dim((street_dir-1) >> 1), pref_dir((street_dir-1)&1);
		
		if (p1[street_dim] != p2[street_dim]) { // choose the point closest to the street as it will likely also be closest to the power lines
			ridge_pos = (((p1[street_dim] < p2[street_dim]) ^ pref_dir) ? 0.01 : 0.99); // near the end of the roof, but leave some space for the connecting pole
		}
	}
	if (ridge_pos == 0.0) {ridge_pos = rgen.rand_uniform(0.1, 0.9);} // select a rand point along ridge if not set, not too close to the edge to avoid the chimney
	ppts.push_back((1.0 - ridge_pos)*p1 + ridge_pos*p2); // interpolate along the ridgeline
	return 1;
}

void building_t::add_solar_panels(rand_gen_t &rgen) { // for houses
	int const roof_tquad_ix(find_main_roof_tquad_ix(rgen, 1)); // skip_if_has_other_obj=1
	if (roof_tquad_ix < 0) return; // not found, possibly not a house or has other objects
	assert((unsigned)roof_tquad_ix < roof_tquads.size());
	tquad_with_ix_t const &roof(roof_tquads[roof_tquad_ix]);
	assert(roof.npts == 4);
	cube_t const roof_bcube(roof.get_bcube());
	float const thickness(0.075*get_window_vspace());
	vector3d const normal(roof.get_norm()), bias(thickness*normal); // slightly above the roof
	tquad_with_ix_t panel(4, tquad_with_ix_t::TYPE_SOLAR);
	point const *const pts(roof.pts);
	// shrink and make rectangular
	point const center(0.25*(pts[0] + pts[1] + pts[2] + pts[3]) + bias);
	vector3d const v10(pts[1] - pts[0]), v32(pts[2] - pts[3]), v30(pts[3] - pts[0]), v21(pts[2] - pts[1]); // 4 sides of roof quad
	float const mu1(v10.mag()), mu2(v32.mag()), mv1(v30.mag()), mv2(v21.mag()); // 4 side lengths
	// min side lengths shrunk, then halved; include vector sum because we want the orthogonal projection
	float const mu(rgen.rand_uniform(0.75, 0.85)*0.5*min(min(mu1, mu2), 0.5f*(v10 + v32).mag())), mv(rgen.rand_uniform(0.75, 0.85)*0.5*min(min(mv1, mv2), 0.5f*(v30 + v21).mag()));
	if (4.0f*mu*mv < 0.25f*roof_bcube.dx()*roof_bcube.dy()) return; // skip if less than 25% of rectangularized roof area (or try another root quad?)
	vector3d const u((v10/mu1 + v32/mu2).get_norm()), v((v30/mv1 + v21/mv2).get_norm());
	panel.pts[0] = center - mu*u - mv*v;
	panel.pts[1] = center + mu*u - mv*v;
	panel.pts[2] = center + mu*u + mv*v;
	panel.pts[3] = center - mu*u + mv*v;
	if (has_chimney && get_chimney().intersects(panel.get_bcube())) return; // chimney intersection
	roof_tquads.push_back(panel);
	// add sides, similar to thick_poly_to_sides()
	tquad_t bottom(panel);
	for (unsigned n = 0; n < 4; ++n) {bottom[n] -= bias;}
	tquad_with_ix_t side(4, tquad_with_ix_t::TYPE_TRIM);

	for (unsigned i = 0; i < 4; ++i) {
		unsigned const inext((i+1)&3);
		side[0] = panel [i];
		side[1] = bottom[i];
		side[2] = bottom[inext];
		side[3] = panel [inext];
		roof_tquads.push_back(side);
	}
}
void building_t::add_rooftop_sat_dish(cube_t const &top, float radius, rand_gen_t &rgen) {
	point dish_center;

	for (unsigned d = 0; d < 2; ++d) {
		bool const ddir(rgen.rand_bool());
		dish_center[d] = top.d[d][ddir] + (ddir ? -1.0 : 1.0)*1.0*radius;
	}
	roof_obj_t dish(ROOF_OBJ_SAT_DISH);
	dish.set_from_point(dish_center);
	dish.expand_by_xy(radius);
	set_cube_zvals(dish, top.z2(), (top.z2() + 2.5*radius));
	details.push_back(dish);
}
cube_t building_t::get_top_building_part() const {
	cube_t top; // uppermost part
	auto parts_end(get_real_parts_end_inc_sec());

	for (auto p = parts.begin(); p != parts_end; ++p) {
		if (top.is_all_zeros() || p->z2() > top.z2()) {top = *p;}
	}
	assert(top.is_strictly_normalized());
	return top;
}
void building_t::add_sat_dish(rand_gen_t &rgen) { // for houses
	float const radius(rgen.rand_uniform(0.2, 0.3)*get_window_vspace());
	add_rooftop_sat_dish(get_top_building_part(), radius, rgen);
}
void building_t::add_tv_antenna(rand_gen_t &rgen) { // for houses
	cube_t const top(get_top_building_part());
	bool const long_dim(rgen.rand_bool());
	float const hlen(0.6*get_window_vspace()), hwidth(0.6*hlen), height(1.8*hlen);
	roof_obj_t ant(ROOF_OBJ_TV_ANT);
	set_cube_zvals(ant, top.z2(), (top.z2() + height));

	for (unsigned d = 0; d < 2; ++d) {
		bool const ddir(rgen.rand_bool());
		set_wall_width(ant, (top.d[d][ddir] + (ddir ? -1.0 : 1.0)*hlen), ((d == long_dim) ? hlen : hwidth), d);
	}
	details.push_back(ant);
}

void building_t::maybe_inv_rotate_pos_dir(point &pos, vector3d &dir) const {
	if (is_rotated()) {
		do_xy_rotate_inv(bcube.get_cube_center(), pos);
		do_xy_rotate_normal_inv(dir);
	}
}

bool building_t::check_cube_contained_in_part(cube_t const &c) const {
	auto parts_end(get_real_parts_end_inc_sec());

	for (auto p = parts.begin(); p != parts_end; ++p) {
		if (p->contains_cube(c)) return 1;
	}
	return 0;
}

bool building_room_geom_t::cube_intersects_moved_obj(cube_t const &c, int ignore_obj_id) const {
	for (auto i = moved_obj_ids.begin(); i != moved_obj_ids.end(); ++i) {
		assert(*i < objs.size());
		if (int(*i) != ignore_obj_id && objs[*i].intersects(c)) return 1; // intersects a moved object; assumes object is collidable
	}
	return 0;
}

void rotate_and_shift_door(tquad_with_ix_t &door, float angle, float shift, bool dim, bool swap_sides) {
	float const rot_angle(-float(angle)*TO_RADIANS*(swap_sides ? -1.0 : 1.0)), sin_term(sin(rot_angle)), cos_term(cos(rot_angle));
	for (unsigned i = 1; i < 3; ++i) {do_xy_rotate(sin_term, cos_term, door.pts[0], door.pts[i]);}
	UNROLL_4X(door.pts[i_][dim] += shift;);
}
tquad_with_ix_t building_t::set_door_from_cube(cube_t const &c, bool dim, bool dir, unsigned type, float pos_adj, bool exterior,
	float open_amt, bool opens_out, bool opens_up, bool swap_sides, bool open_min_amt, door_rotation_t &drot) const
{
	tquad_with_ix_t door(4, type); // quad
	bool const opened(open_amt > 0.0);
	float const pos(c.d[dim][0] + (opened ? 0.0 : pos_adj*(dir ? 1.0 : -1.0))); // move away from wall slightly (not needed if opened)
	door.pts[0].z = door.pts[1].z = c.z1(); // bottom
	door.pts[2].z = door.pts[3].z = c.z2(); // top
	door.pts[0][!dim] = door.pts[3][!dim] = c.d[!dim][ dir]; //  dir side
	door.pts[1][!dim] = door.pts[2][!dim] = c.d[!dim][!dir]; // !dir side
	door.pts[0][ dim] = door.pts[1][ dim] = door.pts[2][dim] = door.pts[3][dim] = pos;
	if (dim == 0) {swap(door.pts[0], door.pts[1]); swap(door.pts[2], door.pts[3]);} // swap two corner points to flip winding dir and invert normal for doors oriented in X

	if (opened) { // rotate >= 90 degrees about pts[0]/pts[3] or pts[1]/pts[2]
		if (opens_up) { // rotates inward and upward 90 degrees
			door.pts[0].z    = door.pts[1].z    = c.z2();
			door.pts[0][dim] = door.pts[1][dim] = pos + c.dz()*((dir ^ opens_out) ? 1.0 : -1.0);
		}
		else { // rotates to the side
			bool const check_backroom_walls(interior && interior->has_backrooms && c.z1() < ground_floor_z1 && !get_basement().contains_cube(c));
			float const wall_thickness(get_wall_thickness()), floor_thickness(get_floor_thickness());
			float const width(c.get_sz_dim(!dim)), signed_width(width*((dir ^ opens_out) ? 1.0 : -1.0));
			float const offset(0.005*width*((dir ^ dim) ? 1.0 : -1.0)*open_amt); // move slightly away from the wall to prevent z-fighting
			for (unsigned i = 0; i < 4; ++i) {door.pts[i][!dim] += offset;}
			// open exactly 90 degrees to start
			door.pts[1-swap_sides][!dim] = door.pts[2+swap_sides][!dim] = door.pts[swap_sides][!dim]; // 1,2=0 / 0,3=1
			door.pts[1][dim] = door.pts[2][dim] = pos + signed_width;
			
			if (!exterior) { // interior
				bool const has_moved_objs(has_room_geom() && !interior->room_geom->moved_obj_ids.empty());
				bool const in_ext_basement(point_in_extended_basement_not_basement(c.get_cube_center()));
				float const doorway_width(get_doorway_width());
				unsigned max_angle(75); // in degrees
				drot.shift = 0.07*signed_width*open_amt;

				if (in_ext_basement && has_room_geom()) { // only open 105 degrees if door opens to a jail
					cube_t tc(c);
					tc.translate_dim(dim, signed_width); // move into the room on the side it opens to
					if (has_bcube_int(tc, interior->room_geom->jails)) {max_angle = 15.0;}
				}
				for (; max_angle > 0; max_angle -= 15) { // try to open door as much as 75 degrees in steps of 15 degrees
					if (open_min_amt) continue;
					tquad_with_ix_t door_rot(door); // cache orig 90 degree open door in case we need to revert it
					rotate_and_shift_door(door_rot, max_angle, drot.shift, dim, swap_sides);
					cube_t test_bcube(door_rot.get_bcube());
					test_bcube.expand_in_dim(dim, -wall_thickness ); // shrink in other dim to avoid intersecting with other part/walls when this door separates two parts
					test_bcube.expand_in_dim(2,   -floor_thickness); // shrink a bit in z to avoid picking up objects from stacks above or below
					if (!in_ext_basement && !check_cube_contained_in_part(test_bcube)) continue; // extends outside part
					cube_t walls_test_bcube(test_bcube);
					walls_test_bcube.expand_in_dim(!dim, wall_thickness); // expand slightly to leave a bit of a gap between walls, and space for whiteboards
					if (has_bcube_int(walls_test_bcube, interior->walls[!dim]))            continue; // hits perp wall
					// open doors don't really block the player from entering or exiting stairs since they can be walked through or closed
					if (interior->is_blocked_by_stairs_or_elevator(test_bcube, 0.0, 0, 1)) continue; // hits stairs or elevator; dmin=0, elevators_only=0, no_check_enter_exit=1
					if (check_backroom_walls && interior->room_geom->cube_int_backrooms_walls(walls_test_bcube)) continue;
					unsigned const first_door_ix(in_ext_basement ? interior->ext_basement_door_stack_ix : 0);
					unsigned const last_door_ix((in_ext_basement || !has_ext_basement()) ? interior->door_stacks.size() : interior->ext_basement_door_stack_ix);
					bool int_other_door(0);

					for (unsigned dsix = first_door_ix; dsix < last_door_ix; ++dsix) {
						door_stack_t const &ds(interior->door_stacks[dsix]);
						if (ds.z1() >= c.z2() || ds.z2() <= c.z1()) continue; // wrong floor
						if (ds.get_true_bcube().intersects(c))      continue; // skip self
						if (ds.get_open_door_path_bcube().intersects(test_bcube)) {int_other_door = 1; break;}
					}
					if (int_other_door) continue;

					if (have_walkway_ext_door) { // check if too close to a walkway exit door or its exit sign
						bool close_to_exit(0);

						for (tquad_with_ix_t const &d : doors) {
							cube_t c(d.get_bcube());
							c.expand_in_dim((c.dy() < c.dx()), doorway_width);
							if (c.intersects(test_bcube)) {close_to_exit = 1; break;}
						}
						if (close_to_exit) continue;
					}
					// check if the player moved an object that would block this door
					if (has_moved_objs) {
						cube_t union_cube(test_bcube);
						union_cube.union_with_cube(door.get_bcube()); // include the full path from 90 degrees to ensure the door can be swung open
						if (interior->room_geom->cube_intersects_moved_obj(union_cube)) continue;
					}
					break; // success - done
				} // for max_angle
				drot.angle = (max_angle + 90.0)*open_amt;
				float const target_angle(drot.angle - 90.0); // includes an XY swap; can be positive or negative
				rotate_and_shift_door(door, target_angle, drot.shift, dim, swap_sides);
			}
		}
	}
	return door;
}
tquad_with_ix_t building_t::set_interior_door_from_cube(door_t const &door) const {
	unsigned const type(is_residential() ? (unsigned)tquad_with_ix_t::TYPE_IDOOR : (unsigned)tquad_with_ix_t::TYPE_ODOOR); // house/hotel/apartment or office door
	door_rotation_t drot; // returned values are unused
	return set_door_from_cube(door, door.dim, door.open_dir, type, 0.0, 0, door.open_amt, 0, 0, door.hinge_side, door.use_min_open_amt(), drot);
}
cube_t building_t::get_door_bounding_cube(door_t const &door) const {
	tquad_with_ix_t const door_tq(set_interior_door_from_cube(door));
	vector3d const normal(door_tq.get_norm());
	cube_t door_bcube(door_tq.get_bcube());
	door_bcube.expand_by(door.get_thickness()*vector3d(fabs(normal.x), fabs(normal.y), 0.0));
	return door_bcube;
}

bool building_t::add_door(cube_t const &c, unsigned part_ix, bool dim, bool dir, bool for_office_building, bool roof_access, bool courtyard, bool for_walkway) { // exterior door

	if (c.is_all_zeros()) return 0;
	vector3d const sz(c.get_size());
	assert(sz[dim] == 0.0 && sz[!dim] > 0.0 && sz.z > 0.0);
	// if it's an office building with two doors already added, make this third door a back metal door
	unsigned const type(for_office_building ? ((doors.size() == 2 && !courtyard && !for_walkway) ?
		(unsigned)tquad_with_ix_t::TYPE_BDOOR2 : (unsigned)tquad_with_ix_t::TYPE_BDOOR) : (unsigned)tquad_with_ix_t::TYPE_HDOOR);
	door_rotation_t drot; // return value is unused
	// exterior=1, open_amt=0.0, opens_out=opens_up=swap_sides=open_min_amt=0
	doors.push_back(set_door_from_cube(c, dim, dir, type, 1.5*get_door_shift_dist(), 1, 0.0, 0, 0, 0, 0, drot));
	if (!roof_access && part_ix < 4) {door_sides[part_ix] |= 1 << (2*dim + dir);}
	if (roof_access) {doors.back().type = tquad_with_ix_t::TYPE_RDOOR;}
	return 1;
}

unsigned extend_roof(cube_t &top, float extend_to, bool dim) {
	if (extend_to == 0.0) return 2; // extend in neither dim
	if (extend_to < top.d[dim][0]) {top.d[dim][0] = extend_to; return 0;} // lo side extend
	if (extend_to > top.d[dim][1]) {top.d[dim][1] = extend_to; return 1;} // hi side extend
	return 0; // extend in neither dim
}

// roof made from two sloped quads and two triangles
float building_t::gen_peaked_roof(cube_t const &top_, float peak_height, bool dim, float extend_to, float max_dz, unsigned skip_side_tri) {

	cube_t top(top_); // deep copy
	unsigned const extend_dir(extend_roof(top, extend_to, dim)), prev_num_roof_tquads(roof_tquads.size());
	float const width(top.get_sz_dim(!dim)), roof_dz(min(max_dz, min(peak_height*width, top.dz())));
	float const z1(top.z2()), z2(z1 + roof_dz), x1(top.x1()), y1(top.y1()), x2(top.x2()), y2(top.y2());
	assert(roof_dz > 0.0);
	point pts[6] = {point(x1, y1, z1), point(x1, y2, z1), point(x2, y2, z1), point(x2, y1, z1), point(x1, y1, z2), point(x2, y2, z2)};
	if (dim == 0) {pts[4].y = pts[5].y = 0.5f*(y1 + y2);} // yc
	else          {pts[4].x = pts[5].x = 0.5f*(x1 + x2);} // xc
	unsigned const qixs[2][2][4] = {{{0,3,5,4}, {4,5,2,1}}, {{0,4,5,1}, {4,3,2,5}}}; // 2 quads
	roof_tquads.reserve(prev_num_roof_tquads + (3 + (extend_dir == 2))); // 2 roof quads + 1-2 side triangles
	unsigned const tixs[2][2][3] = {{{1,0,4}, {3,2,5}}, {{0,3,4}, {2,1,5}}}; // 2 triangles

	for (unsigned n = 0; n < 2; ++n) { // triangle section/wall from z1 up to roof
		if (n == extend_dir || n == skip_side_tri) continue; // skip this side
		bool skip(0);

		// exclude tquads contained in/adjacent to other parts, considering only the cube parts;
		// yes, a triangle side can be occluded by a cube + another opposing triangle side from a higher wall of the house, but it's uncommon, complex, and currently ignored
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (is_basement(p)) continue; // skip the basement
			if (p->d[dim][!n] != top.d[dim][n]) continue; // not the opposing face
			if (p->z1() <= z1 && p->z2() >= z2 && p->d[!dim][0] <= top.d[!dim][0] && p->d[!dim][1] >= top.d[!dim][1]) {skip = 1; break;}
		}
		if (skip) continue;
		tquad_t tquad(3); // triangle
		UNROLL_3X(tquad.pts[i_] = pts[tixs[dim][n][i_]];);
		if (z1 < bcube.z2() && tquad.pts[2][dim] != bcube.d[dim][n]) {tquad.pts[2][dim] += (n ? -1.0 : 1.0)*0.001*width;} // shift peak in slightly to prevent z-fighting with int walls
		roof_tquads.emplace_back(tquad, (unsigned)tquad_with_ix_t::TYPE_WALL); // tag as wall
	} // for n
	bool const extend_roof = 0; // disabled for now, but somewhat works

	if (extend_roof) { // extend the roof outside the wall a small amount
		// may require updating bcube for drawing; causes problems with L/T shaped houses with roof intersecting other walls and interiors; requires two sided drawing
		float const extend(0.05*width), extend_dz(2.0*extend*roof_dz/width);
		for (unsigned n = 0; n < 4; ++n) {pts[n].z -= extend_dz;} // so that slope is preserved
		if (dim == 0) {pts[0].y -= extend; pts[1].y += extend; pts[2].y += extend; pts[3].y -= extend;}
		else          {pts[0].x -= extend; pts[1].x -= extend; pts[2].x += extend; pts[3].x += extend;}
	}
	for (unsigned n = 0; n < 2; ++n) { // roof
		tquad_t tquad(4); // quad
		UNROLL_4X(tquad.pts[i_] = pts[qixs[dim][n][i_]];);

		if (prev_num_roof_tquads > 0 && has_L_shaped_roof_area()) {
			// clip to existing roof tquad from a previous part;
			// the top and bottom edges may meet the other tquad at exactly the same zval, which may not count as intersection, so use the centerline instead
			point const p1(0.5*(tquad.pts[0] + tquad.pts[dim ? 1 : 3])), p2(0.5*(tquad.pts[dim ? 3 : 1] + tquad.pts[2]));
			
			for (auto i = roof_tquads.begin(); i != roof_tquads.begin()+prev_num_roof_tquads; ++i) {
				if (!i->is_roof()) continue;
				vector3d const normal(i->get_norm());
				point p_int; // used in the m loop below
				if (!line_poly_intersect(p1, p2, i->pts, i->npts, normal, p_int)) continue;

				// if any edge intersects this tquad, clip all horizontal edges to it; this handles the case where edges meet exactly
				for (unsigned m = 0; m < 4; ++m) {
					point &A(tquad.pts[m]), &B(tquad.pts[(m+1)&3]);
					if (A.z != B.z) continue; // not a horizontal edge, skip
					float t(0.0); // unused
					if (!line_int_plane(A, B, i->pts[0], normal, p_int, t, 0)) continue; // should never fail
					bool const clip_side(dot_product_ptv(normal, A, B) < 0.0);
					(clip_side ? A : B) = p_int;
				} // for m
				break; // should only intersect once
			} // for i
		}
		roof_tquads.emplace_back(tquad, (unsigned)tquad_with_ix_t::TYPE_ROOF_PEAK); // tag as peaked roof
	} // for n
	has_attic_window = 1; // if there's an attic, it may have windows with this roof type
	return roof_dz;
}

// roof made from two sloped quads + two sloped triangles
float building_t::gen_hipped_roof(cube_t const &top_, float peak_height, float extend_to) {

	bool const dim(get_largest_xy_dim(top_)); // always the largest dim
	cube_t top(top_); // deep copy
	extend_roof(top, extend_to, dim); // Note: return value is unused
	float const width(top.get_sz_dim(!dim)), length(top.get_sz_dim(dim)), offset(0.5f*(length - width)), roof_dz(min(peak_height*width, top.dz()));
	float const z1(top.z2()), z2(z1 + roof_dz), x1(top.x1()), y1(top.y1()), x2(top.x2()), y2(top.y2());
	assert(roof_dz > 0.0);
	point const center(0.5f*(x1 + x2), 0.5f*(y1 + y2), z2);
	point pts[6] = {point(x1, y1, z1), point(x1, y2, z1), point(x2, y2, z1), point(x2, y1, z1), center, center};
	if (dim) {UNROLL_3X(swap(pts[0], pts[i_+1]);)} // rotate 1 vertex CCW
	pts[4][dim] -= offset; pts[5][dim] += offset; // move points away from center to create ridgeline
	unsigned const qixs[2][4] = {{0,3,5,4}, {2,1,4,5}};
	unsigned const tixs[2][3] = {{1,0,4}, {3,2,5}};
	unsigned const start_ix(roof_tquads.size());
	roof_tquads.reserve(start_ix + 4);

	for (unsigned n = 0; n < 2; ++n) {
		roof_tquads.emplace_back(4, tquad_with_ix_t::TYPE_ROOF_HIP); // quad
		UNROLL_4X(roof_tquads.back().pts[i_] = pts[qixs[n][i_]];)
		roof_tquads.emplace_back(3, tquad_with_ix_t::TYPE_ROOF_HIP); // triangle
		UNROLL_3X(roof_tquads.back().pts[i_] = pts[tixs[n][i_]];)
	}
	return roof_dz;
}

float building_t::gen_sloped_roof_for_stacked_parts(cube_t const &bot, cube_t const &top) {

	assert(bot.contains_cube_xy(top));
	float const window_vspace(get_window_vspace());
	float roof_dz(window_vspace*WINDOW_BORDER_MULT*get_window_v_border() - 0.04f*window_vspace); // bottom of window sill
	float const roof_z1(bot.z2()), roof_z2(roof_z1 + roof_dz);
	max_eq(roof_dz, 0.05f*window_vspace); // make sure it's not near zero
	
	for (unsigned dim = 0; dim < 2; ++dim) { // x/y
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (bot.d[dim][dir] == top.d[dim][dir]) continue; // no roof on this side
			// add sloped roof tquad; if either end is inset, treat it as a hipped roof with no gutters
			bool inset[2] = {};
			for (unsigned d = 0; d < 2; ++d) {inset[d] = (top.d[!dim][d] != bot.d[!dim][d]);}
			tquad_with_ix_t roof(4, ((inset[0] || inset[1]) ? tquad_with_ix_t::TYPE_ROOF_HIP : tquad_with_ix_t::TYPE_ROOF_PEAK)/*tquad_with_ix_t::TYPE_ROOF_SLOPE*/);
			roof.pts[0][ dim] = roof.pts[1][ dim] = bot.d[ dim][dir]; // lower edge
			roof.pts[2][ dim] = roof.pts[3][ dim] = top.d[ dim][dir]; // upper edge
			roof.pts[0].z     = roof.pts[1].z     = roof_z1; // lower edge
			roof.pts[2].z     = roof.pts[3].z     = roof_z2; // upper edge
			roof.pts[0][!dim] = bot.d[!dim][1];
			roof.pts[1][!dim] = bot.d[!dim][0];
			roof.pts[2][!dim] = top.d[!dim][0];
			roof.pts[3][!dim] = top.d[!dim][1];
			if (dim ^ dir) {swap(roof.pts[1], roof.pts[3]);} // use correct winding order
			roof_tquads.push_back(roof);
			
			for (unsigned d = 0; d < 2; ++d) { // add required wall triangles
				if (inset[d]) continue; // inset on this edge, not flush; no wall needed
				tquad_with_ix_t wall(3, tquad_with_ix_t::TYPE_WALL); // triangle
				wall.pts[0].z = wall.pts[1].z = roof_z1;
				wall.pts[2].z = roof_z2;
				wall.pts[0][dim] = bot.d[dim][dir];
				wall.pts[1][dim] = wall.pts[2][dim] = top.d[dim][dir];
				UNROLL_3X(wall.pts[i_][!dim] = bot.d[!dim][d];);
				if (d ^ dir ^ dim) {swap(wall.pts[1], wall.pts[2]);} // use correct winding order
				roof_tquads.push_back(wall);
			} // for d
		} // for dir
	} // for dim
	return roof_dz;
}

void building_t::gen_building_doors_if_needed(rand_gen_t &rgen) { // for office buildings

	if (num_sides > 8)                       return; // can't place doors on curved building sides
	if (!is_cube() && has_complex_floorplan) return; // this case isn't handled either
	assert(!parts.empty());
	float const door_height(get_office_bldg_door_height()), wscale(DOOR_WIDTH_SCALE_OFFICE); // a bit taller and a lot wider than house doors

	if (has_pri_hall()) { // building has primary hallway, place doors at both ends of first part
		for (unsigned d = 0; d < 2; ++d) {
			cube_t const door(place_door(parts.front(), bool(hallway_dim), d, door_height, 0.0, 0.0, 0.0, wscale, 0, 0, rgen));
			if (add_door(door, 0, bool(hallway_dim), d, 1)) {floor_ext_door_mask |= 1;}
		}
		return;
	}
	bool const pref_dim(rgen.rand_bool()), pref_dir(rgen.rand_bool()), bldg_has_windows(has_windows());
	bool used[4] = {0,0,0,0}; // per-side, not per-base cube
	// at least 2 doors unless it's a small rectangle (large rectangle will have a central hallway with doors at each end)
	unsigned const min_doors((parts.size() > 1 || is_industrial()) ? 2 : 1);
	unsigned const max_doors(is_parking() ? 2 : (bldg_has_windows ? 3 : 4)); // buildings with windows have at most 3 doors since they're smaller
	unsigned const num_doors(min_doors + (rgen.rand() % (max_doors - min_doors + 1)));
	cube_t const parts_bcube(get_unrotated_parts_bcube());

	for (unsigned num = 0; num < num_doors; ++num) {
		bool placed(0);

		for (auto b = parts.begin(); b != get_real_parts_end() && !placed; ++b) { // try all different ground floor parts
			if (is_basement(b)) continue; // skip the basement
			unsigned const part_ix(b - parts.begin());
			if (bldg_has_windows && part_ix >= 4) break; // only first 4 parts can have doors - must match first floor window removal logic
			if (b->z1() > ground_floor_z1)        break; // moved off the ground floor

			if (!is_cube()) { // N-gon building; See if there's an X or Y oriented side to add a door to; should get here for exactly one part
				vect_point const &points(get_part_ext_verts(part_ix));
				float const toler(0.1*get_wall_thickness()); // add a bit of tolerance to account for FP error
				unsigned match_side_ix(0);
				bool found(0);

				for (auto i = points.begin(); i != points.end() && !found; ++i) {
					point const &p1(*i), &p2((i == points.begin()) ? points.back() : *(i-1));
					
					for (unsigned d = 0; d < 2; ++d) { // check for horizontal or vertical dim
						if (fabs(p1[d] - p2[d]) > toler) continue; // not aligned in this axis
						if (match_side_ix != num)        continue; // not the selected side; should this be randomized? most of the time there's only one valid side
						++match_side_ix;
						float const side_pos(0.5*(p1[d] + p2[d])), side_center(0.5*(p1[!d] + p2[!d]));
						bool const dir(side_pos > b->get_center_dim(d)), allow_fail(!doors.empty());
						if (!add_door(place_door(*b, d, dir, door_height, side_center, side_pos, 0.0, wscale, allow_fail, 0, rgen), part_ix, d, dir, 1)) continue;
						placed = found = 1;
						break;
					} // for d
				} // for i
				continue;
			}
			if (interior && interior->ind_info) { // industrial building custom doors
				bool const dim(interior->ind_info->entrance_dim), dir(interior->ind_info->entrance_dir);

				if (num == 0) { // main entrance is in or between the smaller rooms
					if (add_door(place_door(*b, dim, dir, door_height, interior->ind_info->entrance_pos, 0.0, 0.0, wscale, 0, 0, rgen), part_ix, dim, dir, 1)) {
						used[2*dim + dir] = 1; // mark used
						placed = 1;
						continue;
					}
				}
				else if (num == 1 && is_warehouse()) { // second door is warehouse loading garage door, opposite the main entrance door
					bool added(0);

					for (unsigned n = 0; n < 2; ++n) { // {centered, non-centered}
						float const door_center_shift(n ? 0.4 : 0.05); // first pass is (nearly) centered, second pass is non-centered
						cube_t const door(place_door(*b, dim, !dir, door_height, 0.0, 0.0, door_center_shift, 2.0*wscale, 1, 1, rgen));

						if (add_door(door, part_ix, dim, !dir, 1)) {
							doors.back().type = tquad_with_ix_t::TYPE_GDOOR;
							used[2*dim + (!dir)] = 1; // mark used
							driveway = door;
							driveway.z1() = ground_floor_z1;
							driveway.z2() = ground_floor_z1 + 0.005*door_height;
							driveway.d[dim][ dir] = door.d[dim][!dir];
							driveway.d[dim][!dir] += (dir ? -1.0 : 1.0)*1.5*door_height;
							added = 1;
							break;
						}
					} // for n
					if (added) {continue;} // not counted as placed
				}
			}
			for (unsigned n = 0; n < 4; ++n) {
				bool const dim(pref_dim ^ bool(n>>1)), dir(pref_dir ^ bool(n&1));
				if (used[2*dim + dir]) continue; // door already placed on this side
				// find a side on the exterior to ensure door isn't obstructed by a building cube or in the courtyard
				if (b->d[dim][dir] != parts_bcube.d[dim][dir]) continue;
				bool const allow_fail(!doors.empty() || b+1 != get_real_parts_end()); // allow door placement to fail if we've already placed at least one door of not last part
				if (!add_door(place_door(*b, dim, dir, door_height, 0.0, 0.0, 0.1, wscale, allow_fail, 0, rgen), part_ix, dim, dir, 1)) continue;
				used[2*dim + dir] = 1; // mark used
				placed = 1;
				break;
			} // for n
		} // for b
	} // for num
	//assert(!doors.empty()); // must have placed at least one door - too strong for rotated buildings?

	if (has_courtyard) { // add a door opening into the courtyard
		assert(parts.size() >= 4);
		unsigned const pref_side(rgen.rand()&3);
		bool placed(0);

		for (unsigned p = 0; p < 4 && !placed; ++p) {
			unsigned const part_ix((p + pref_side) & 3);
			cube_t const &part(parts[part_ix]);
			if (part.z1() != ground_floor_z1) continue; // not on the ground floor

			for (unsigned n = 0; n < 4 && !placed; ++n) {
				bool const dim(n>>1), dir(n&1);
				if (part.d[dim][dir] == parts_bcube.d[dim][dir] || part.d[dim][!dir] != parts_bcube.d[dim][!dir]) continue; // find a side on the interior
				unsigned const door_ix(doors.size());
				placed = add_door(place_door(part, dim, dir, door_height, 0.0, 0.0, 0.0, wscale, 1, 0, rgen), part_ix, dim, dir, 1, 0, 1); // centered, courtyard
				if (placed) {courtyard_door_ix = int8_t(door_ix);}
			}
		} // for b
	}
	if (!doors.empty()) {floor_ext_door_mask |= 1;} // I suppose courtyard doors count here
}

cube_t building_t::register_deck_and_get_part_bounds(cube_t const &deck) {
	assert(deck_bounds.is_all_zeros()); // can only have one deck
	deck_bounds = deck;
	cube_t deck_exp(deck), part_bounds;
	deck_exp.expand_by(get_wall_thickness()); // find nearby parts that don't quite overlap; maybe checking for adjacency is close enough?

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // should garages and sheds count?
		if (i->z1() != ground_floor_z1) continue; // not ground floor
		if (i->intersects_xy(deck)) {part_bounds.assign_or_union_with_cube(*i);}
	}
	return part_bounds;
}

void building_t::place_roof_ac_units(unsigned num, float sz_scale, cube_t const &bounds, vect_cube_t const &avoid, rand_gen_t &rgen) {

	if (num == 0) return;
	vector3d cube_sz(sz_scale, 1.5*sz_scale, 1.2*sz_scale); // consistent across all AC units (note that x and y will be doubled)
	if (rgen.rand_bool()) {swap(cube_sz.x, cube_sz.y);} // other orientation
	unsigned const num_floors(round_fp(bounds.dz()/get_window_vspace()));
	float const duct_height_scale(rgen.rand_uniform(0.7, 0.9)), duct_width_scale(rgen.rand_uniform(0.4, 0.7)); // consistent per part
	bool add_array(num_floors > 1);

	for (unsigned i = 0; i < num; ++i) {
		unsigned const acs_start(details.size());
		roof_obj_t c(ROOF_OBJ_AC);
		bool placed(0);

		for (unsigned n = 0; n < 100 && !placed; ++n) { // limited to 100 attempts to prevent infinite loop
			c.set_from_point(point(rgen.rand_uniform(bounds.x1(), bounds.x2()), rgen.rand_uniform(bounds.y1(), bounds.y2()), bounds.z2()));
			c.expand_by_xy(cube_sz);
			c.z2() += cube_sz.z; // z2
			if (!bounds.contains_cube_xy(c)) continue; // not contained
			if (has_bcube_int_no_adj(c, avoid) || has_bcube_int_no_adj(c, details)) continue; // intersects avoid cubes or other detail objects (inc_adj=0)
			placed = has_ac = 1;
		} // for n
		if (!placed) break; // failed, exit loop
		details.push_back(c);
		if ((rgen.rand() & 7) == 0) {swap(cube_sz.x, cube_sz.y);} // swap occasionally

		if (add_array && rgen.rand_bool()) { // try to add a row or 2D grid of AC units
			unsigned const nx(max(1U, (unsigned)round_fp(rgen.rand_uniform(0.5, 1.5)*sqrt(num_floors))));
			unsigned const ny(max(1U, (unsigned)round_fp(rgen.rand_uniform(0.5, 1.5)*num_floors/nx)));
			float const spacing(5.0*sz_scale);
			unsigned num_added(1); // one was added above
			
			for (unsigned y = 0; y < ny; ++y) {
				for (unsigned x = 0; x < nx; ++x) {
					if (x == 0 && y == 0) continue; // already placed
					roof_obj_t ac(c);
					ac.translate(vector3d(x*spacing, y*spacing, 0.0));
					if (!bounds.contains_cube_xy(ac) || has_bcube_int_no_adj(ac, avoid) || has_bcube_int_no_adj(ac, details)) continue;
					details.push_back(ac);
					++num_added;
					++i; // this counts as an added AC unit
				} // for y
			} // for y
			if (num_added >= nx*ny/2) {add_array = 0;} // we placed our array - done
		}
		// add ducts
		unsigned const acs_end(details.size());
		bool const dim(c.dy() < c.dx()); // short dim
		bool dir(rgen.rand_bool());
		float const length(c.get_sz_dim(!dim)), duct_width((1.0 - duct_width_scale)*length);
		float const min_len(min(0.25f*length, 1.1f*duct_width)), init_len(max(min_len, rgen.rand_uniform(1.0, 3.0)*length));
		roof_obj_t duct(ROOF_OBJ_DUCT);
		set_cube_zvals(duct, c.z1(), (c.z1() + duct_height_scale*c.dz()));

		for (unsigned i = acs_start; i != acs_end; ++i) {
			roof_obj_t const &ac(details[i]);
			assert(ac.type == ROOF_OBJ_AC);
			set_wall_width(duct, ac.get_center_dim(!dim), 0.5*duct_width, !dim);

			for (unsigned d = 0; d < 2; ++d) { // dir
				float const edge(ac.d[dim][dir]);
				duct.d[dim][!dir] = edge;
				bool placed(0);

				for (float duct_len = init_len; duct_len > min_len; duct_len *= 0.75) { // try to shorten the length if placement fails
					duct.d[dim][dir] = edge + (dir ? 1.0 : -1.0)*duct_len;
					cube_t test_cube(duct);
					test_cube.d[dim][dir] = duct.d[dim][dir] + (dir ? 1.0 : -1.0)*max(0.6f*duct_width, duct.dz()); // add extra spacing at the end
					if (!bounds.contains_cube_xy(test_cube) || has_bcube_int_no_adj(test_cube, avoid) || has_bcube_int_no_adj(test_cube, details)) continue; // skip if bad
					details.push_back(duct);
					placed = 1;
					break;
				}
				if (placed) break;
				dir ^= 1;
			} // for d
		} // for i
	} // for n
}

void building_t::add_roof_walls(cube_t const &c, float wall_width, bool overlap_corners, cube_t out[4]) {
	// add walls around the roof; we don't need to draw all sides of all cubes, but it's probably not worth the trouble to sort it out
	float const height(6.0f*wall_width), width(wall_width), shorten(overlap_corners ? 0.0f : wall_width), z1(c.z2()), z2(c.z2()+height);
	out[0] = cube_t(c.x1(), c.x2(), c.y1(), c.y1()+width, z1, z2);
	out[1] = cube_t(c.x1(), c.x2(), c.y2()-width, c.y2(), z1, z2);
	out[2] = cube_t(c.x1(), c.x1()+width, c.y1()+shorten, c.y2()-shorten, z1, z2);
	out[3] = cube_t(c.x2()-width, c.x2(), c.y1()+shorten, c.y2()-shorten, z1, z2);
}

void subtract_part_from_walls(cube_t const &s, vect_cube_t &cubes, float wall_width) {
	unsigned iter_end(cubes.size()); // capture size before splitting

	for (unsigned i = 0; i < iter_end; ++i) {
		cube_t const &c(cubes[i]);
		if (!c.intersects_no_adj(s)) continue; // keep it
		cube_t shared(c);
		shared.intersect_with_cube(s);
		if (shared.dx() < 1.1*wall_width && shared.dy() < 1.1*wall_width) continue; // clipped off only the end of the wall (inside corner), keep entire wall
		subtract_cube_from_cube_inplace(s, cubes, i, iter_end); // Note: invalidates c reference
	} // for i
}

void merge_cubes_xy(vect_cube_t &cubes) {
	for (unsigned i = 0; i < cubes.size(); ++i) {
		for (unsigned j = 0; j < cubes.size(); ++j) {
			if (i == j) continue; // skip self
			cube_t &a(cubes[i]), &b(cubes[j]);
			if (a.is_all_zeros() || b.is_all_zeros()) continue; // one of the cubes has already been removed
			if (a.z1() != b.z1() || a.z2() != b.z2()) continue; // zvals not shared (doesn't actually happen for building walls)

			for (unsigned d = 0; d < 2; ++d) { // try to merge b into a
				if (a.d[!d][0] != b.d[!d][0] || a.d[!d][1] != b.d[!d][1]) continue; // range not shared in this dim
				if (a.d[d][1] >= b.d[d][0] && a.d[d][0] <= b.d[d][1]) {min_eq(a.d[d][0], b.d[d][0]); max_eq(a.d[d][1], b.d[d][1]); b = cube_t(); break;}
			}
		} // for j
	} // for i
}

void building_t::gen_details(rand_gen_t &rgen, bool is_rectangle) { // for the roof; office buildings, not houses

	bool const flat_roof(roof_type == ROOF_TYPE_FLAT), add_walls(is_cube() && flat_roof); // simple cube buildings with flat roofs
	cube_t flat_roof_area(get_flat_roof_section_bcube());
	bool const has_flat_upper_roof(!flat_roof_area.is_all_zeros()), has_flat_top(flat_roof || has_flat_upper_roof);
	unsigned num_ac_units(0); // cube buildings only for now, no parking structures; 0-8 for industrial, otherwise 0-6
	if (has_flat_top && is_cube() && !is_parking()) {num_ac_units = (rgen.rand() % (is_industrial() ? 9 : 7));}
	float const window_vspacing(get_window_vspace()), wall_width(0.049*window_vspacing); // slightly narrower than interior wall width to avoid z-fighting with roof access
	assert(!parts.empty());
	create_per_part_ext_verts(); // needed for roof containment queries
	add_company_sign(rgen);

	if (!is_rectangle) { // polygon roof, can only add walls and AC units
		if (add_walls && rgen.rand_bool()) { // add walls 50% of the time
			cube_t temp[4];
			vect_cube_t cubes, to_add;

			for (auto p = parts.begin(); p != parts.end(); ++p) {
				if (p->z2() < bcube.z2()) continue; // not the top floor
				add_roof_walls(*p, wall_width, 1, temp); // overlap_corners=1 so that clipping works correctly
				cubes.resize(4);
				for (unsigned i = 0; i < 4; ++i) {cubes[i] = temp[i];}
				float const z2(cubes[0].z2());

				for (auto p2 = parts.begin(); p2 != parts.end(); ++p2) { // remove any interior wall sections that overlap another top part
					if (p2 == p || p2->z2() < p->z2()) continue; // skip self and non-top floors
					
					for (unsigned d = 0; d < 2; ++d) {
						cube_t clip_cube(*p2);
						clip_cube.z2() = z2;
						float const exp(d ? wall_width : -wall_width);
						clip_cube.expand_by(exp, -exp, 0.0);
						subtract_part_from_walls(clip_cube, cubes, wall_width);
					}
				} // for p2
				vector_add_to(cubes, to_add);
			} // for p
			merge_cubes_xy(to_add); // optimization

			for (auto i = to_add.begin(); i != to_add.end(); ++i) {
				if (!i->is_all_zeros()) {details.emplace_back(*i, (uint8_t)ROOF_OBJ_WALL);}
			}
		}
		if (num_ac_units > 0) {
			float const xy_sz(0.75*bcube.get_size_xy().mag()*rgen.rand_uniform(0.012, 0.02)); // size based on bcube

			for (auto p = parts.begin(); p != parts.end(); ++p) {
				unsigned const num_this_part(min(num_ac_units, unsigned(num_ac_units*p->dx()*p->dy()/(bcube.dx()*bcube.dy()) + 1))); // distribute based on area
				place_roof_ac_units(num_this_part, xy_sz, *p, parts, rgen);
			}
		}
		details.shrink_to_fit();
		return;
	}
	cube_t const top_part(parts.back()); // top/last part
	cube_t top(has_flat_upper_roof ? flat_roof_area : top_part); // top/last part
	top.z1() = top_part.z1();
	vector3d const tsz(top.get_size()), tpsz(top_part.get_size());
	point top_center(top.get_cube_center());
	if (is_rotated() && !is_cube()) {do_xy_rotate_inv(bcube.get_cube_center(), top_center);} // put top_center in building bcube coordinate space
	float const helipad_radius(2.0*window_vspacing);
	bool const can_have_hp_or_sl(flat_roof && num_sides >= 4 && flat_side_amt == 0.0 && !is_house && !is_parking());
	has_helipad = (can_have_hp_or_sl && min(tsz.x, tsz.y) > (is_cube() ? 3.2 : 4.0)*helipad_radius && bcube.dz() > 8.0*window_vspacing && (rgen.rand() % 12) == 0);
	cube_t avoid_bcube, door_blocker;

	if (has_helipad) { // add helipad
		tquad_t helipad(4); // quad
		float const z(top.z2() + 0.02*window_vspacing); // slightly above the roof to avoid Z-fighting
		float const x1(top_center.x - helipad_radius), x2(top_center.x + helipad_radius), y1(top_center.y - helipad_radius), y2(top_center.y + helipad_radius);
		bool const dir(rgen.rand_bool()); // R90 50% of the time
		helipad.pts[ dir+0   ].assign(x1, y1, z);
		helipad.pts[ dir+1   ].assign(x2, y1, z);
		helipad.pts[ dir+2   ].assign(x2, y2, z);
		helipad.pts[(dir+3)&3].assign(x1, y2, z);
		roof_tquads.emplace_back(helipad, (uint8_t)tquad_with_ix_t::TYPE_HELIPAD);
		avoid_bcube = helipad.get_bcube();
	}
	else if (can_have_hp_or_sl && is_cube() && interior_enabled()) {
		maybe_add_skylight(rgen);
	}
	bool const add_rooftop_door(!is_cube() && has_complex_floorplan /*&& has_helipad*/); // add if there's no interior/stairs to the roof
	unsigned num_blocks(0);
	
	if (flat_roof && !is_parking() && skylights.empty()) { // no roof blocks if there are roof quads (houses, etc.), skylights, or this is a parking garage
		num_blocks = (rgen.rand() % 9); // 0-8
		if (add_rooftop_door) {max_eq(num_blocks, 1U);} // at least one for the door
	}
	bool const add_antenna((flat_roof || roof_type == ROOF_TYPE_SLOPE) && !has_helipad && !is_parking() && skylights.empty() && rgen.rand_bool());
	bool const can_add_special_obj(has_flat_top && !has_helipad && !add_antenna && !is_parking() && skylights.empty()); // open roof space with no blocker
	bool const add_water_tower(can_add_special_obj && (tsz.x < 2.0*tsz.y && tsz.y < 2.0*tsz.x) && (tsz.x > 0.5*tpsz.x && tsz.y > 0.5*tpsz.y) && rgen.rand_bool());
	bool const add_sat_dish(can_add_special_obj && !add_water_tower && is_cube() && !is_industrial());
	unsigned const num_details(num_blocks + num_ac_units + 4*add_walls + add_antenna + add_water_tower + add_sat_dish);
	if (num_details == 0) return; // nothing to do
	if (add_walls && min(tsz.x, tsz.y) < 4.0*wall_width) return; // too small
	float const xy_sz(tpsz.xy_mag()); // better to use bcube for size?
	cube_t bounds(top);
	if (add_walls) {bounds.expand_by_xy(-wall_width);}
	reserve_extra(details, num_details);

	if (add_water_tower) {
		float const radius(rgen.rand_uniform(0.04, 0.06)*(tpsz.x + tpsz.y)), height(rgen.rand_uniform(3.0, 5.0)*radius);
		roof_obj_t wtower(ROOF_OBJ_WTOWER);
		point wt_center(top_center);
		for (unsigned d = 0; d < 2; ++d) {wt_center[d] += rgen.rand_uniform(-1.0, 1.0)*0.25*tsz[d];} // apply a random shift
		wtower.set_from_point(wt_center);
		wtower.expand_by_xy(radius);

		if (check_part_contains_cube_xy(top, parts.size()-1, wtower)) {
			set_cube_zvals(wtower, top.z2(), (bcube.z2() + height)); // z2 uses bcube to include sloped roof
			details.push_back(wtower);
			avoid_bcube = wtower;
		}
	}
	if (add_sat_dish) {
		float const radius(rgen.rand_uniform(0.015, 0.025)*(tpsz.x + tpsz.y));
		add_rooftop_sat_dish(top, radius, rgen);
		avoid_bcube = details.back();
		// Note: orient is chosen to point away from the roof later later since there's no dim/dir in roof objects
	}
	for (unsigned i = 0; i < num_blocks; ++i) {
		roof_obj_t c(ROOF_OBJ_BLOCK); // generic block
		float const height_scale(0.0035f*(tsz.z + bcube.dz())); // based on avg height of current section and entire building
		float height(height_scale*rgen.rand_uniform(1.0, 4.0));
		bool placed(0);

		for (unsigned n = 0; n < 100; ++n) { // limited to 100 attempts to prevent infinite loop
			c.set_from_point(point(rgen.rand_uniform(bounds.x1(), bounds.x2()), rgen.rand_uniform(bounds.y1(), bounds.y2()), top.z2()));
			c.expand_by_xy(vector3d(xy_sz*rgen.rand_uniform(0.01, 0.07), xy_sz*rgen.rand_uniform(0.01, 0.07), 0.0));
			if (!bounds.contains_cube_xy(c)) continue; // not contained
			if (!avoid_bcube .is_all_zeros() && c.intersects_xy(avoid_bcube ))      continue; // bad placement
			if (!door_blocker.is_all_zeros() && c.intersects_xy(door_blocker))      continue; // bad placement
			if (!is_cube() && !check_part_contains_cube_xy(top, parts.size()-1, c)) continue; // not contained in roof
			placed = 1;
			break;
		} // for n
		if (!placed) break; // failed, exit loop

		if (add_rooftop_door && i == 0) { // make the first block tall and add a door (that the player can't open)
			float const door_width(get_doorway_width());

			if (min(c.dx(), c.dy()) > 1.5*door_width) { // block is large enough to place a door
				bool const door_dim(rgen.rand_bool());
				bool const door_dir((c.d[door_dim][0] - top.d[door_dim][0]) < (top.d[door_dim][1] - c.d[door_dim][1])); // closer to building center/further from roof edge
				float const door_height(get_door_height()), center(c.get_center_dim(!door_dim));
				float const wall_pos(c.d[door_dim][door_dir] + (door_dir ? 1.0 : -1.0)*0.02*door_width); // move slightly away from the wall
				float const lr_offset(((door_dim ^ door_dir) ? 1.0 : -1.0)*0.5*door_width);
				max_eq(height, 1.1f*door_height);
				tquad_t door(4);
				for (unsigned n = 0; n < 4; ++n) {door.pts[n][door_dim] = wall_pos;}
				door.pts[0].z = door.pts[1].z = c.z1(); // bottom
				door.pts[2].z = door.pts[3].z = c.z1() + door_height; // top
				door.pts[0][!door_dim] = door.pts[3][!door_dim] = center - lr_offset; // left  side
				door.pts[1][!door_dim] = door.pts[2][!door_dim] = center + lr_offset; // right side
				roof_tquads.emplace_back(door, (unsigned)tquad_with_ix_t::TYPE_RDOOR2);
				door_blocker = door.get_bcube();
				door_blocker.d[door_dim][door_dir] += (door_dir ? 1.0 : -1.0)*1.2*get_doorway_width(); // add door clearance
				has_fake_roof_door = 1;
			}
		}
		c.z2() += height; // z2
		details.push_back(c);
	} // for i
	if (add_antenna) { // add antenna
		float const height(max(2.0f*window_vspacing, rgen.rand_uniform(0.25, 0.5)*tsz.z));
		float const radius(min(0.025f*height, 0.003f*rgen.rand_uniform(1.0, 2.0)*(tsz.x + tsz.y)));
		roof_obj_t antenna(ROOF_OBJ_ANT);
		antenna.set_from_point(top_center); // always in the center of the roof
		antenna.expand_by_xy(radius);
		set_cube_zvals(antenna, top.z2(), (bcube.z2() + height)); // z2 uses bcube to include sloped roof
		details.push_back(antenna);
		has_antenna = 1;
	}
	if (num_ac_units > 0) {
		vect_cube_t ac_avoid;
		for (cube_t const &s : skylights) {ac_avoid.push_back(s); ac_avoid.back().expand_in_z(0.25*window_vspacing);} // avoid skylights
		if (!avoid_bcube .is_all_zeros()) {ac_avoid.push_back(avoid_bcube );}
		if (!door_blocker.is_all_zeros()) {ac_avoid.push_back(door_blocker);}
		place_roof_ac_units(num_ac_units, xy_sz*rgen.rand_uniform(0.012, 0.02), bounds, ac_avoid, rgen);
	}
	if (add_walls) {
		cube_t cubes[4];
		add_roof_walls(top, wall_width, 0, cubes); // overlap_corners=0
		for (unsigned i = 0; i < 4; ++i) {details.emplace_back(cubes[i], (uint8_t)ROOF_OBJ_WALL);}
	}
	if (is_factory() || is_powerplant()) {add_smokestack(rgen);}
	for (roof_obj_t const &o : details) {assert(o.is_strictly_normalized()); max_eq(bcube.z2(), o.z2());} // extend bcube z2 to contain details
	if (has_flat_top) {gen_grayscale_detail_color(rgen, 0.2, 0.6);} // for antenna and roof
}

void building_t::maybe_add_skylight(rand_gen_t &rgen) {
	if (is_industrial()) return; // no industrial building skylights; they work, but factories, etc. have enough windows already
	if (is_parking   ()) return; // not needed for parking structures
	// maybe add skylights; cube roofs only for now, since we can't cut holes in other shapes; only for office buildings with interiors;
	// could add skylights to house rooms such as bathrooms, but they haven't been assigned and we would need to cut holes in the roof and possibly attic
	// note that at this point there has been no floorplanning, so we don't know where primary hallways, etc. will be
	float part_zmax(bcube.z1());
	for (cube_t const &part : parts) {max_eq(part_zmax, part.z2());}

	// add skylights to top floor roofs
	for (auto p = parts.begin(); p != parts.end(); ++p) { // at this point, all parts should be main building parts
		if (p->z2() != part_zmax) continue; // not top floor
		cube_t roof_ceiling(*p), skylight;
		roof_ceiling.z1() = p->z2() - get_fc_thickness();

		if (can_use_hallway_for_part(p - parts.begin())) { // if we have a hallway, align the skylight to it
			float num_hall_windows, hall_width, room_width;
			cube_t const hall(get_hallway_for_part(*p, num_hall_windows, hall_width, room_width));
			bool const hall_dim(hall.dx() < hall.dy());
			skylight = hall;
			skylight.expand_in_dim( hall_dim, -rgen.rand_uniform(0.1, 0.3)*hall.get_sz_dim(hall_dim)); // shrink
			skylight.expand_in_dim(!hall_dim, -0.5*get_wall_thickness()); // shrink slightly
		}
		else {
			for (unsigned d = 0; d < 2; ++d) {
				float const hwidth(rgen.rand_uniform(0.1, 0.2)*roof_ceiling.get_sz_dim(d));
				if (rgen.rand_bool()) {set_wall_width(skylight, roof_ceiling.get_center_dim(d), hwidth, d);} // centered
				else {set_wall_width(skylight, rgen.rand_uniform((roof_ceiling.d[d][0] + 1.5*hwidth), (roof_ceiling.d[d][1] - 1.5*hwidth)), hwidth, d);} // misaligned
			}
		}
		copy_zvals(skylight, roof_ceiling);
		skylights.push_back(skylight);
	} // for p
}

cube_t building_t::get_helipad_bcube() const {
	assert(has_helipad);

	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
		if (i->type == tquad_with_ix_t::TYPE_HELIPAD) return i->get_bcube();
	}
	assert(0); // if has_helipad was set, a helipad must be found
	return cube_t(); // never gets here
}

void building_t::maybe_add_special_roof(rand_gen_t &rgen) {
	assert(!parts.empty());
	cube_t const &top(parts.back());
	vector3d const sz(top.get_size()); // top/last part

	if (global_building_params.onion_roof && num_sides >= 16 && flat_side_amt == 0.0) { // cylinder building
		if (sz.x < 1.2*sz.y && sz.y < 1.2*sz.x && sz.z > max(sz.x, sz.y)) {roof_type = ROOF_TYPE_ONION;}
	}
	else if (is_cube()) { // only simple cubes are handled
		if (global_building_params.dome_roof && sz.x < 1.2*sz.y && sz.y < 1.2*sz.x && sz.z > max(sz.x, sz.y)) {roof_type = ROOF_TYPE_DOME;} // roughly square, not too short
		else if (parts.size() == 1 && bcube.dz() < 4.5*get_window_vspace() && rgen.rand_bool()) { // shorter single cube building; likely to become factory or warehouse
			cube_t const &top(parts[0]);
			roof_type = ROOF_TYPE_CURVED;
			max_eq(bcube.z2(), (top.z2() + 0.125f*min(top.dx(), top.dy()))); // should curve in the short dim by 25% of radius
		}
		else {gen_sloped_roof(rgen, top);} // sloped roof
	}
	if      (roof_type == ROOF_TYPE_DOME ) {max_eq(bcube.z2(), (top.z2() + 0.5f*max(sz.x, sz.y)));}
	else if (roof_type == ROOF_TYPE_ONION) {max_eq(bcube.z2(), (top.z2() + 1.0f*max(sz.x, sz.y)));}
}
void building_t::gen_sloped_roof(rand_gen_t &rgen, cube_t const &top) { // Note: currently not supported for rotated buildings

	float const peak_height(rgen.rand_uniform(0.2, 0.5));
	float const wmin(min(top.dx(), top.dy())), z1(top.z2()), z2(z1 + peak_height*wmin), x1(top.x1()), y1(top.y1()), x2(top.x2()), y2(top.y2());
	point const pts[5] = {point(x1, y1, z1), point(x1, y2, z1), point(x2, y2, z1), point(x2, y1, z1), point(0.5f*(x1 + x2), 0.5f*(y1 + y2), z2)};
	float const d1(rgen.rand_uniform(0.0, 0.8));

	if (d1 < 0.2) { // pointed roof with 4 sloped triangles (hipped)
		unsigned const ixs[4][3] = {{1,0,4}, {3,2,4}, {0,3,4}, {2,1,4}};
		roof_tquads.reserve(4);

		for (unsigned n = 0; n < 4; ++n) {
			roof_tquads.emplace_back(3, tquad_with_ix_t::TYPE_ROOF_OFFICE); // triangle
			UNROLL_3X(roof_tquads.back().pts[i_] = pts[ixs[n][i_]];)
		}
	}
	else { // flat roof with center quad and 4 surrounding sloped quads
		point const center((1.0 - d1)*pts[4]);
		point pts2[8];
		for (unsigned n = 0; n < 4; ++n) {pts2[n] = pts[n]; pts2[n+4] = d1*pts[n] + center;}
		unsigned const ixs[5][4] = {{4,7,6,5}, {0,4,5,1}, {3,2,6,7}, {0,3,7,4}, {2,1,5,6}}; // add the flat quad first, which works better for sphere intersections
		roof_tquads.reserve(5);

		for (unsigned n = 0; n < 5; ++n) {
			roof_tquads.emplace_back(4, tquad_with_ix_t::TYPE_ROOF_OFFICE); // quad
			UNROLL_4X(roof_tquads.back().pts[i_] = pts2[ixs[n][i_]];)
		}
	}
	roof_type = ROOF_TYPE_SLOPE;
	add_roof_to_bcube(); // is max_eq(bcube.z2(), z2) good enough?
	gen_grayscale_detail_color(rgen, 0.4, 0.8); // for antenna and roof
}

void building_t::add_roof_to_bcube() {
	for (auto const &tq : roof_tquads) {
		tq.update_bcube(bcube); // technically should only need to update z2
		
		if (has_attic()) { // use roof tquads to include the attic space
			if (!tq.is_roof()) continue;
			for (unsigned n = 0; n < tq.npts; ++n) {max_eq(interior_z2, tq.pts[n].z);}
		}
	} // for tq
}
void building_t::gen_grayscale_detail_color(rand_gen_t &rgen, float imin, float imax) {
	float cscale(rgen.rand_uniform(imin, imax));
	if (is_parking() || is_industrial()) {cscale = 0.5 + 0.5*cscale;} // lighten for parking structures and industrial since they use metal and concrete
	detail_color = colorRGBA(cscale, cscale, cscale, 1.0);
}


// *** Interiors ***

void building_t::expand_ground_floor_cube(cube_t &cube, cube_t const &skip) const {
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // try to expand the cube to cover more parts (pri parts only)
		if (p->z1() != ground_floor_z1) continue; // only count ground floor parts
		if (skip.contains_cube_xy(*p))  continue; // already contained, skip
		cube_t cand_ge(cube);
		cand_ge.union_with_cube(*p);
		if (cand_ge.get_area_xy() < 1.05f*(cube.get_area_xy() + p->get_area_xy())) {cube = cand_ge;} // union mostly includes the two parts
		else if (real_num_parts && is_cube_contained_in_parts(cand_ge)) {cube = cand_ge;} // useful for cross-shaped building
		else {try_expand_into_xy(cube, *p);}
	}
}

void building_t::get_exclude_cube(point const &pos, cube_t const &skip, cube_t &exclude, bool camera_in_building) const {
	float const cube_pad(4.0*grass_width*(camera_in_building ? 2.0 : 1.0)), extent(bcube.get_max_dim_sz());
	float dmin_sq(extent*extent); // start with a large value, squared
	exclude.set_to_zeros();

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) { // find closest part, including garages/sheds
		if (p->z1() != ground_floor_z1) continue; // only count ground floor parts
		if (skip.contains_cube_xy(*p))  continue; // already contained, skip
		float const dist_sq(p2p_dist_sq(pos, p->closest_pt(pos)));
		if (dist_sq < dmin_sq) {exclude = *p; dmin_sq = dist_sq;} // keep if closest part to pos
	}
	if (exclude.is_all_zeros()) return; // not found (only ground floor part was skipped)
	expand_ground_floor_cube(exclude, skip);
	exclude.expand_by_xy(cube_pad); // exclude grass blades that partially intersect the building interior
}

void building_t::update_grass_exclude_at_pos(point const &pos, vector3d const &xlate, bool camera_in_building) const {
	get_exclude_cube(pos, cube_t(),       grass_exclude1, camera_in_building); // first/closest  cube
	get_exclude_cube(pos, grass_exclude1, grass_exclude2, camera_in_building); // second closest cube
	if (!grass_exclude1.is_all_zeros()) {grass_exclude1 = get_rotated_bcube(grass_exclude1) + xlate;}
	if (!grass_exclude2.is_all_zeros()) {grass_exclude2 = get_rotated_bcube(grass_exclude2) + xlate;}
}

void building_t::update_stats(building_stats_t &s) const { // calculate all of the counts that are easy to get
	++s.nbuildings;
	s.nparts   += parts.size();
	s.ndetails += details.size();
	s.ntquads  += roof_tquads.size();
	s.ndoors   += doors.size();
	if (!interior) return;
	++s.ninterior;
	s.nrooms  += interior->rooms.size();
	s.nceils  += interior->ceilings.size();
	s.nfloors += interior->floors.size();
	s.nwalls  += interior->walls[0].size() + interior->walls[1].size();
	s.ndoors  += interior->doors.size(); // I guess these also count as doors?
	if (!interior->room_geom) return;
	++s.nrgeom;
	s.nobjs  += interior->room_geom->objs.size();
	s.nverts += interior->room_geom->get_num_verts();
}

bool door_opens_inward(door_base_t const &door, cube_t const &room) {
	return (room.is_all_zeros() || (door.d[door.dim][0] < room.get_center_dim(door.dim)) == door.open_dir); // null room always returns 1 (conservative)
}
// check_dirs/open_dirs: 0=lo only, 1=hi only, 2+=both
bool is_cube_close_to_door(cube_t const &c, float dmin, bool inc_open, cube_t const &door, unsigned check_dirs, unsigned open_dirs, bool allow_block_door) {
	if (c.z2() < door.z1() || c.z1() > door.z2()) return 0;
	bool const dim(door.dy() < door.dx());
	float const width(door.get_sz_dim(!dim)), height(width/DOOR_WIDTH_SCALE); // door may be a stack, can't use dz(), calculate height from width
	float const trim_width(0.5*WALL_THICK_VAL*height); // width of door trim in both dims (half wall width)
	// we don't know how much the door is open and in which direction, so expand by door width in both dirs to be conservative
	float const keepout(inc_open ? width : trim_width);
	float const lo_edge(door.d[!dim][0] - ((check_dirs == 1) ? trim_width : keepout)), hi_edge(door.d[!dim][1] + ((check_dirs == 0) ? trim_width : keepout));
	if (c.d[!dim][0] > hi_edge || c.d[!dim][1] < lo_edge) return 0; // no overlap in !dim
	float dmin_lo(dmin), dmin_hi(dmin);
	
	if (!allow_block_door) { // max with door width so that door has space to open, in dirs where the door can open
		if (open_dirs != 1) {max_eq(dmin_lo, width+trim_width);}
		if (open_dirs != 0) {max_eq(dmin_hi, width+trim_width);}
	}
	return (c.d[dim][0] < door.d[dim][1]+dmin_hi && c.d[dim][1] > door.d[dim][0]-dmin_lo); // within min_dist
}
bool building_t::is_cube_close_to_exterior_doorway(cube_t const &c, float dmin, bool inc_open) const {
	for (auto i = doors.begin(); i != doors.end(); ++i) { // test exterior doors
		bool const is_garage_door(i->type == tquad_with_ix_t::TYPE_GDOOR), check_open(inc_open && !is_garage_door); // garage doors open up, not sideways
		if (is_cube_close_to_door(c, dmin, check_open, i->get_bcube(), 2)) return 1; // check both dirs
	}
	return 0;
}
// Note: inc_open only applies to interior doors
bool building_t::is_cube_close_to_doorway(cube_t const &c, cube_t const &room, float dmin, bool inc_open, bool check_open_dir) const {
	// Note: we want to test this for things like stairs, but exterior doors likely haven't been allocated at this point, so we have to check for that during door placement
	if (is_cube_close_to_exterior_doorway(c, dmin, inc_open)) return 1;
	return (interior ? interior->is_cube_close_to_doorway(c, room, dmin, inc_open, check_open_dir) : 0); // test interior doors
}
bool building_interior_t::is_cube_close_to_doorway(cube_t const &c, cube_t const &room, float dmin, bool inc_open, bool check_open_dir) const { // ignores zvals
	if (door_stacks.empty()) return 0; // single room house with no door?
	cube_t test_cube(c);
	// assume all doors are the same size and use the last for reference, but pad by 1.5x anyway; upper bound on the door bcube when open any amount in any direction
	float const door_width(max(dmin, max(door_stacks.back().dx(), door_stacks.back().dy())));
	test_cube.expand_by_xy(1.5*door_width);
	
	if (!room.is_all_zeros()) { // don't go outside the original room
		if (!room.intersects_xy(c)) {cout << "Bad room cube: " << TXT(room.str()) << TXT(c.str()) << endl;}
		assert(room.intersects_xy(c)); // can't test zval in case of ball in pool
		cube_t room_exp(room);
		room_exp.expand_by_xy(0.1*door_width); // small amount to include wall width
		test_cube.intersect_with_cube_xy(room_exp);
	}
	for (door_stack_t const &ds : door_stacks) { // interior doors
		if (!test_cube.intersects(ds)) continue; // optimization
		if (is_cube_close_to_door(c, dmin, (inc_open && door_opens_inward(ds, room)), ds, ds.get_check_dirs(), (check_open_dir ? ds.open_dir : 2))) return 1;
	}
	for (cube_t const &w : open_walls) { // open walls count as doorways, even though there's no door
		cube_t wall_exp(w);
		bool const dim(w.dy() < w.dx());
		wall_exp.expand_in_dim( dim,      door_width); // expand outward
		wall_exp.expand_in_dim(!dim, 0.05*door_width); // extend a bit to avoid placing objects right on the edge
		if (wall_exp.intersects(c)) return 1;
	}
	return 0;
}

cube_t get_stairs_bcube_expanded(stairwell_t const &s, float ends_clearance, float sides_clearance, float doorway_width) {
	cube_t tc(s);
	tc.expand_in_dim(s.dim, ends_clearance); // add extra space at both ends of stairs; may only need to add on open ends, but this is difficult to check for
	float const floor_spacing(doorway_width/DOOR_WIDTH_SCALE), wall_hw(s.get_wall_hwidth(floor_spacing));
	tc.expand_in_dim(!s.dim, (sides_clearance + wall_hw)); // add extra space to account for walls and railings on stairs
	return tc;
}
void get_L_stairs_entrances(stairs_landing_base_t const &s, float doorway_width, bool for_placement, cube_t entrances[2]) { // {lower, upper}
	bool const dirs[2] = {s.dir, s.bend_dir};
	float const landing_width(doorway_width), extend_amt((for_placement ? 1.0 : 1.1)*doorway_width); // must be >= doorway_width for correct AI path finding

	for (unsigned d = 0; d < 2; ++d) { // {lower, upper}
		bool const dim(s.dim ^ bool(d)), dir(dirs[d] ^ bool(d)), dir2(dirs[!d] ^ bool(d));
		cube_t &se(entrances[d]);
		se = s; // copy stairs; will return full zval range
		// set width, but include the full edge *plus* width for lower entrance placement queries, since we don't want to block the path between wall and railing
		bool const block_full_side(for_placement && d == 0);
		se.d[!dim][dir2] = s.d[!dim][dir2 ^ block_full_side ^ 1] + (dir2 ? 1.0 : -1.0)*landing_width; // set width
		se.d[dim][0   ]  = se.d[ dim][1] = s.d[dim][!dir]; // edge of stairs entrance
		se.d[dim][!dir] += (dir ? -1.0 : 1.0)*extend_amt; // extend away from stairs
	} // for d
	// extend the top of the stairs exit toward the wall to avoid blocking the space between the railing and wall on this side as well
	if (for_placement) {entrances[1].d[s.dim][s.dir] += (s.dir ? 1.0 : -1.0)*doorway_width;}
}
// no_check_enter_exit: 0=check enter and exit with expand, 1=check enter and exit with no expand, 2=don't check enter and exit at all
bool has_stairs_bcube_int(cube_t const &bcube, vect_stairwell_t const &stairs, float doorway_width, float pad_added_to_bcube, int no_check_enter_exit) {
	cube_t pre_test(bcube);
	if (no_check_enter_exit < 2) {pre_test.expand_by_xy(doorway_width);}
	float const approx_floor_spacing(doorway_width/DOOR_WIDTH_SCALE); // used for L-shaped stairs

	for (stairwell_t const &s : stairs) {
		if (!s.intersects(pre_test)) continue; // early termination test optimization
		cube_t const tc(get_stairs_bcube_expanded(s, ((no_check_enter_exit == 2) ? 0.0 : doorway_width), 0.0, doorway_width)); // sides_clearance=0.0
		
		if (tc.intersects(bcube)) { // intersects core stairs
			if (!s.is_l_shape()) return 1; // L-shaped stairs are special
			if (bcube.z2() >= s.z1() + approx_floor_spacing) return 1; // not on the bottom floor - counts as intersection
			float const landing_width(doorway_width);
			bool const dirs[2] = {s.dir, s.bend_dir};

			for (unsigned d = 0; d < 2; ++d) { // {lower, upper}
				bool const dim(s.dim ^ bool(d)), dir2(dirs[!d] ^ bool(d));
				cube_t seg(tc);
				seg.d[!dim][dir2] = s.d[!dim][!dir2] + (dir2 ? 1.0 : -1.0)*landing_width; // set width
				if (seg.intersects(bcube)) return 1;
			}
		}
		// extra check for objects blocking the entrance/exit to the side; this is really only needed for open ends, but helps to avoid squeezing objects behind stairs as well
		if (no_check_enter_exit || s.is_u_shape()) continue; // U-shaped stairs are only open on one side and generally placed in hallways, so ignore
		
		if (s.is_l_shape()) { // L-shaped stairs are special
			cube_t entrances[2]; // {lower, upper}
			get_L_stairs_entrances(s, doorway_width, 1, entrances); // for_placement=1

			if (bcube.z2() < s.z1() + approx_floor_spacing) { // bottom floor
				// lower entrance can be smaller - subtract off the width of the second/bend segment
				entrances[0].d[!s.dim][s.bend_dir] -= (s.bend_dir ? 1.0 : -1.0)*(s.get_width() - doorway_width);
			}
			//entrances[0].z2() -= approx_floor_spacing; // lower entrance is not on upper floor; but we don't want to block the path between railing and wall either
			entrances[1].z1() += approx_floor_spacing; // upper entrance is not on lower floor
			if (entrances[0].intersects(bcube) || entrances[1].intersects(bcube)) return 1;
			continue;
		}
		// expand end width, but don't need to expand as much if bcube has already been expanded from the original object
		float const end_width_expand(0.75*doorway_width - pad_added_to_bcube);
		if (end_width_expand <= 0.0) continue;

		for (unsigned e = 0; e < 2; ++e) { // for each end (entrance/exit)
			cube_t end(tc);
			end.d[s.dim][!e] = s.d[s.dim][e]; // shrink to the gap between *s and tc
			if (!s.against_wall[0]) {end.d[!s.dim][0] -= end_width_expand;}
			if (!s.against_wall[1]) {end.d[!s.dim][1] += end_width_expand;}
			if (end.intersects(bcube)) return 1; // Note: ignores zval test, which can affect bottom exit and top entrances of non-walled stairs
		}
	} // for s
	return 0;
}
bool has_elevator_bcube_int(cube_t const &bcube, vector<elevator_t> const &elevators, float front_pad, float trim_thickness) {
	for (elevator_t const &e : elevators) {
		if (e.get_bcube_padded(front_pad, trim_thickness).intersects(bcube)) return 1; // front pad matches wall clip logic
	}
	return 0;
}
bool has_escalator_bcube_int(cube_t const &bcube, vector<escalator_t> const &escalators, float doorway_width) {
	for (escalator_t const &e : escalators) {
		cube_t tc(e);
		tc.expand_in_dim(e.dim, doorway_width); // add extra space to both ends of the escalator (conservative)
		if (tc.intersects(bcube)) return 1;
	}
	return 0;
}
float building_interior_t::get_doorway_width() const {
	return (doors.empty() ? int_door_width : max(doors.front().dx(), doors.front().dy())); // calculate doorway width from first door
}
float building_t::get_doorway_width() const {
	float const width(interior ? interior->get_doorway_width() : 0.0);
	return (width ? width : DOOR_WIDTH_SCALE*get_door_height()); // calculate from window spacing/door height if there's no interior or no interior doors
}
// and ramps, and extb conn rooms/doors, and escalators
bool building_interior_t::is_blocked_by_stairs_or_elevator(cube_t const &c, float dmin, bool elevators_only, int no_check_enter_exit) const {
	cube_t tc(c);
	tc.expand_by_xy(dmin); // no pad in z
	float const doorway_width(get_doorway_width());
	float const front_pad(max(0.0, (1.5*doorway_width - dmin))); // don't need to double pad
	float const trim_thickness(0.1*WALL_THICK_VAL*doorway_width/(DOOR_WIDTH_SCALE*.95*(1.0 - FLOOR_THICK_VAL_OFFICE)));
	if (has_elevator_bcube_int (tc, elevators,  front_pad, trim_thickness)) return 1;
	if (has_escalator_bcube_int(tc, escalators, doorway_width)) return 1;
	if (elevators_only) return 0;
	tc.z1() -= 0.001*tc.dz(); // expand slightly to avoid placing an object exactly at the top of the stairs
	// must check zval to exclude stairs and elevators in parts with other z-ranges
	if (has_stairs_bcube_int(tc, stairwells, doorway_width, dmin, no_check_enter_exit)) return 1;
	
	if (!pg_ramp.is_all_zeros()) { // Note: ramp is placed during basement room geom generation, which may not have been run yet
		bool const dim(pg_ramp.ix >> 1), dir(pg_ramp.ix & 1);
		cube_t ramp_ext(pg_ramp);
		// if ramp extends upward, extend the top upward a bit so that it intersects anything placed on/near the floor
		if (!ignore_ramp_placement) {ramp_ext.z2() += 0.5*doorway_width;}
		// otherwise we still need to avoid the lower part of the ramp; extend the top down to avoid intersecting objects on the floor above
		else                        {ramp_ext.z2() -= 0.5*doorway_width;}
		ramp_ext.d[dim][dir] += (dir ? 1.0 : -1.0)*pg_ramp.get_sz_dim(!dim); // clear space in front of the ramp equal to its width
		if (ramp_ext.intersects(tc)) return 1;
	}
	if (conn_info != nullptr) { // include extended basement connector doors here because they aren't accounted for in the regular doors check
		for (auto const &c : conn_info->conn) {
			for (auto const &room : c.rooms) {
				if (room.intersects(tc)) return 1;
			}
		}
	}
	return 0;
}
// similar to above (without stairs pretest), but returns bounding cubes rather than checking for intersections; used for parking garage queries
void building_interior_t::get_stairs_and_elevators_bcubes_intersecting_cube(cube_t const &c, vect_cube_t &bcubes, float ends_clearance, float sides_clearance) const {
	float const doorway_width(get_doorway_width());

	for (auto const &s : stairwells) {
		cube_t const tc(get_stairs_bcube_expanded(s, ends_clearance, sides_clearance, doorway_width));
		if (tc.intersects(c)) {bcubes.push_back(tc);}
	}
	for (auto const &e : elevators) {
		cube_t tc(e);
		tc.expand_by_xy(sides_clearance);
		tc.d[e.dim][e.dir] += (ends_clearance - sides_clearance)*(e.dir ? 1.0 : -1.0); // add extra space in front of the elevator
		if (tc.intersects(c)) {bcubes.push_back(tc);}
	}
	// escalators are currently ignored here because they're not in basements/parking garages
}

struct cmp_mall_wall {
	building_interior_t const &interior;
	cmp_mall_wall(building_interior_t const &i) : interior(i) {}
	bool operator()(cube_t const &c) const {return interior.is_inside_mall_stores(c.get_cube_center());}
};

void building_interior_t::sort_for_optimal_culling() {
	for (unsigned d = 0; d < 2; ++d) { // sort walls longest to shortest to improve occlusion culling time
		vect_cube_t &v(walls[d]);
		assert(extb_walls_start[d] <= v.size());

		if (has_mall()) { // move mall back hallway walls after mall interior walls for improved culling and ease of drawing
			mall_hall_walls_start[d] = std::partition(v.begin()+extb_walls_start[d], v.end(), cmp_mall_wall(*this)) - v.begin();
		}
		else {mall_hall_walls_start[d] = v.size();} // no mall, set to end of range
		sort(v.begin(), v.begin()+extb_walls_start[d], cube_by_sz(!d)); // skip exterior basement walls
		sort(v.begin()+extb_walls_start[d], v.begin()+mall_hall_walls_start[d], cube_by_sz(!d)); // only exterior basement walls; excludes mall hall walls (not occluders)
	} // for d
	sort(floors  .begin(), floors  .end(), [](cube_t const &a, cube_t const &b) {return (a.z1() > b.z1());}); // top down,  for early z culling and improved occluder fusion
	sort(ceilings.begin(), ceilings.end(), [](cube_t const &a, cube_t const &b) {return (a.z2() < b.z2());}); // bottom up, for early z culling and improved occluder fusion
}
void building_interior_t::remove_excess_capacity() {
	remove_excess_cap(floors);
	remove_excess_cap(ceilings);
	remove_excess_cap(rooms);
	remove_excess_cap(doors);
	remove_excess_cap(door_stacks);
	remove_excess_cap(landings);
	remove_excess_cap(stairwells);
	remove_excess_cap(elevators);
	for (unsigned d = 0; d < 2; ++d) {remove_excess_cap(walls[d]);}
}
void building_interior_t::finalize() {
	if (rooms.size() > (1U << 16)) {
		std::cerr << "Error: Too many rooms: " << rooms.size() << endl;
		exit(0);
	}
	sort_for_optimal_culling();
	remove_excess_capacity();
}

void apply_fc_cube_max_merge_xy(vect_cube_t &cubes) {
	for (auto c = cubes.begin(); c != cubes.end();) { // Note: no increment
		auto floor_end(c);
		for (floor_end++; floor_end != cubes.end() && floor_end->z1() == c->z1(); ++floor_end) {} // find all ceilings on this floor

		for (auto c1 = c; c1 != floor_end; ++c1) { // source
			for (auto c2 = c; c2 != floor_end; ++c2) { // merge cand
				if (c1 != c2) {try_expand_into_xy(*c1, *c2);}
			}
		}
		c = floor_end;
	} // for c
}
void building_interior_t::create_fc_occluders() {
	if (!fc_occluders.empty() || ceilings.empty()) return; // already done, or no ceilings
	fc_occluders = ceilings; // start by copying ceilings
	cube_t bcube(ceilings[0]);
	
	for (auto &c : fc_occluders) {
		assert(c.is_strictly_normalized());
		c.z2() += c.dz(); // double the thickness to include the floor as well
		bcube.union_with_cube(c);
	}
	// max merge across the ceilings of each floor (occluder fusion)
	apply_fc_cube_max_merge_xy(fc_occluders);
	// remove small slivers
	auto i(fc_occluders.begin()), o(i);

	for (; i != fc_occluders.end(); ++i) {
		assert(i->is_strictly_normalized());
		if (i->dx() > 0.01*bcube.dx() && i->dy() > 0.01*bcube.dy()) {*(o++) = *i;}
	}
	fc_occluders.erase(o, fc_occluders.end());
}

void building_interior_t::add_ceil_floor_pair(cube_t cf, float zc, float z, float zf) {
	set_cube_zvals(cf, zc, z);
	ceilings.push_back(cf);
	set_cube_zvals(cf, z, zf);
	floors.push_back(cf);
}

