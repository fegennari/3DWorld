// 3D World - Building/House Attic Logic
// by Frank Gennari 06/10/2022

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"

colorRGBA get_light_color_temp(float t);
colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c);
void add_boxes_to_space(room_object_t const &c, vect_room_object_t &objects, cube_t const &bounds, vect_cube_t &cubes, rand_gen_t &rgen,
	unsigned num_boxes, float xy_scale, float hmin, float hmax, bool allow_crates, unsigned flags); // from building_room_obj_expand
bool try_add_lamp(cube_t const &place_area, float floor_spacing, unsigned room_id, unsigned flags, float light_amt,
	vect_cube_t &cubes, vect_room_object_t &objects, rand_gen_t &rgen);
bool gen_furnace_cand(cube_t const &place_area, float floor_spacing, bool near_wall, rand_gen_t &rgen, cube_t &furnace, bool &dim, bool &dir);
bool add_obj_to_closet(room_object_t const &c, cube_t const &interior, vect_room_object_t &objects, vect_cube_t &cubes, rand_gen_t &rgen,
	vector3d const &sz, unsigned obj_type, unsigned flags, room_obj_shape shape=SHAPE_CUBE, colorRGBA const &color=WHITE, bool against_back=0);
void narrow_furnace_intake(cube_t &duct, room_object_t const &c);


unsigned building_t::get_attic_part_ix() const { // Note: can be called before adding the attic, so can't check has_attic()
	assert(!parts.empty());
	return ((real_num_parts >= 2 && parts[1].z2() > parts[0].z2()) ? 1 : 0); // if there are at least two parts, use the taller one
}

// Note: may incorrectly return 0 for points exactly between roof tquads, such as those on the centerline of the roof
bool building_t::point_under_attic_roof(point const &pos, vector3d *const cnorm) const {
	if (!get_attic_part().contains_pt_xy(pos)) return 0;

	for (auto const &tq : roof_tquads) {
		if (!is_attic_roof(tq, 1)) continue; // type_roof_only=1
		if (!point_in_polygon_2d(pos.x, pos.y, tq.pts, tq.npts)) continue; // check 2D XY point containment
		vector3d const normal(tq.get_norm());
		if (normal.z == 0.0) continue; // skip vertical sides
		if (cnorm) {*cnorm = -normal;} // we're looking at the underside of the roof, so reverse the normal; set whether or not we're inside the attic
		if (dot_product_ptv(normal, pos, tq.pts[0]) < 0.0) return 1;
	}
	return 0;
}
bool building_t::point_in_attic(point const &pos, vector3d *const cnorm) const {
	if (!has_attic() || pos.z < interior->attic_access.z2() || pos.z > interior_z2) return 0; // test attic floor zval
	return point_under_attic_roof(pos, cnorm);
}
bool building_t::cube_in_attic(cube_t const &c) const {
	if (!has_attic() || c.z2() < interior->attic_access.z2() || c.z1() > interior_z2) return 0; // test attic floor zval
	// test the 4 top corners of the cube
	float const z2(c.z2() + 2.5*get_attic_beam_depth()); // account for attic beam depth, which reduces the ceiling height / increases our effective cube height (approximate)
	return (point_under_attic_roof(point(c.x1(), c.y1(), z2)) || point_under_attic_roof(point(c.x1(), c.y2(), z2)) ||
		    point_under_attic_roof(point(c.x2(), c.y2(), z2)) || point_under_attic_roof(point(c.x2(), c.y1(), z2)));
}

void building_t::get_attic_roof_tquads(vect_tquad_with_ix_t &tquads) const {
	tquads.clear();
	if (!has_attic()) return;

	for (auto const &tq : roof_tquads) {
		if (is_attic_roof(tq, 1)) {tquads.push_back(tq);} // type_roof_only=1
	}
}

void building_t::get_attic_windows(vect_tquad_with_ix_t &tquads, float offset_scale) const {
	if (!is_house || !has_attic() || !has_attic_window) return;
	float const floor_spacing(get_window_vspace()), shift_dist(offset_scale*get_door_shift_dist()), wall_thickness(get_wall_thickness());
	cube_t const &part(get_attic_part()); // attic is always in the tallest/first part
	cube_t avoid;

	if (has_chimney) {
		avoid = get_chimney();
		avoid.expand_by_xy(wall_thickness);
	}
	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
		if (i->type != tquad_with_ix_t::TYPE_WALL || i->npts != 3) continue; // not a wall triangle
		cube_t const bc(i->get_bcube());
		float const height(bc.dz());
		if (height < 0.5*floor_spacing)             continue; // too small (maybe porch roof)
		if (has_chimney && avoid.intersects_xy(bc)) continue; // blocked by chimney
		bool const dim(bc.dy() < bc.dx()); // should be zero size in X or Y
		float const wall_pos(i->pts[0][dim]); // all points should have the same pos
		bool const dir(part.get_center_dim(dim) < wall_pos);
		if (part.d[dim][dir] != bcube.d[dim][dir]) continue; // only exterior facing walls (avoid window clipping through adjacent part)
		cube_t part_test(bc);
		part_test.expand_in_dim( dim,  wall_thickness); // expand to pick up nearby walls
		part_test.expand_in_dim(!dim,- wall_thickness); // shrink to exclude adjacencies
		if (!part_test.intersects(part))           continue; // not on the main part with the attic
		// create square inscribed in triangle
		float const center(bc.get_center_dim(!dim)), base(bc.get_sz_dim(!dim));
		float const window_width(min(0.75f*floor_spacing, height*base/(height + base))); // window width == window height
		tquad_with_ix_t window(4, tquad_with_ix_t::TYPE_HDOOR); // mark as a door so that it's textured from [0.0, 1.0]
		float const window_pos(wall_pos + (dir ? 1.0 : -1.0)*shift_dist); // shift away from the wall
		for (unsigned n = 0; n < 4; ++n) {window.pts[n][dim] = window_pos;} // set window position
		float const bot_edge(bc.z1() + min(0.4f*window_width, 0.1f*bc.dz())); // raise off the floor a bit
		window.pts[0].z = window.pts[1].z = bot_edge; // bottom edge
		window.pts[2].z = window.pts[3].z = bot_edge + window_width; // top edge
		window.pts[0][!dim] = window.pts[3][!dim] = center - 0.5*window_width;
		window.pts[1][!dim] = window.pts[2][!dim] = center + 0.5*window_width;
		if (dim ^ dir ^ 1) {std::reverse(window.pts, window.pts+4);} // reverse the winding order
		tquads.push_back(window);
	} // for i
}

bool building_t::has_L_shaped_roof_area() const {
	if (real_num_parts == 1) return 0; // not L-shaped
	cube_t const &A(parts[0]), &B(parts[1]);
	if (A.z2() != B.z2())    return 0; // not at same level
	if (roof_dims == 2)      return 0; // parallel roof
	if (roof_dims == 1)      return 1; // perpendicular roof
	// secondary part's roof is oriented in long dim; if this is the dim adjacent to the primary part,
	// then the two attic areas are connected forming in L-shape; otherwise, there will be two parallel roof peaks with a valley in between
	bool const adj_x(A.x1() == B.x2() || A.x2() == B.x1()), adj_y(A.y1() == B.y2() || A.y2() == B.y1());
	assert(adj_x != adj_y); // must be adjacent in exactly one dim
	return (B.get_sz_dim(!adj_y) < B.get_sz_dim(adj_y));
}

bool building_t::add_attic_access_door(cube_t const &ceiling, unsigned part_ix, unsigned num_floors, unsigned rooms_start, rand_gen_t &rgen) {
	// roof tquads don't intersect correct on the interior for L-shaped house attics, so skip the attic in this case, for now
	//if (has_L_shaped_roof_area()) return 0;
	float const floor_spacing(get_window_vspace());
	cube_t const &part(parts[part_ix]);
	if (min(part.dx(), part.dy()) < 2.75*floor_spacing) return 0; // must be large enough
	// add a ceiling cutout for attic access
	float const half_len(0.24*floor_spacing), half_wid(0.16*floor_spacing);
	room_t best_room;
	float best_area(0.0);
	bool in_hallway(0);

	for (unsigned r = rooms_start; r < interior->rooms.size(); ++r) {
		room_t const &room(interior->rooms[r]);
		if (room.part_id != part_ix) continue;
		if (room.has_stairs_on_floor(num_floors-1)) continue; // skip room with stairs
		if (max(room.dx(), room.dy()) < 2.5*half_len || min(room.dx(), room.dy()) < 2.5*half_wid) continue; // too small
		if (room.is_hallway) {best_room = room; in_hallway = 1; break;} // hallway is always preferred
		// should we reject this room if there's not enough head clearance above it in the attic?
		float const area(room.dx()*room.dy());
		if (area > best_area) {best_room = room; best_area = area;} // choose room with the largest area
	}
	if (best_room.is_all_zeros()) return 0;
	bool const long_dim(best_room.dx() < best_room.dy());
	cube_t valid_area(best_room);
	valid_area.expand_in_dim( long_dim, -1.2*half_len); // add sufficient clearance
	valid_area.expand_in_dim(!long_dim, -1.2*half_wid); // add sufficient clearance
	if (!valid_area.is_strictly_normalized()) return 0; // not enough space for the door (shouldn't be the case)
	rand_gen_t rgen2(rgen); // deep copy to avoid disrupting rgen state
	point access_pos;

	if (in_hallway) {
		access_pos = best_room.get_cube_center();
		// place off center to avoid blocking center light
		access_pos[ long_dim] += (rgen2.rand_bool() ? -1.0 : 1.0)*0.1*best_room.get_sz_dim(long_dim);
		// place on the side of the hallway to avoid blocking the player and people from walking by
		access_pos[!long_dim]  = (rgen2.rand_bool() ? (best_room.d[!long_dim][0] + 1.1*half_wid) : (best_room.d[!long_dim][1] - 1.1*half_wid));
	}
	else {
		cube_t const &part(get_part_for_room(best_room));
		// if the room spans the entire part, make the attic access in the center so that the stairs have proper clearance
		bool const span_x(best_room.x1() == part.x1() && best_room.x2() == part.x2()), span_y(best_room.y1() == part.y1() && best_room.y2() == part.y2());
		bool const xd(best_room.xc() < part.xc()), yd(best_room.yc() < part.yc()); // closer to the center of the part to maximize head space
		access_pos.x = (span_x ? best_room.xc() : (0.7*best_room.d[0][xd] + 0.3*best_room.d[0][!xd]));
		access_pos.y = (span_y ? best_room.yc() : (0.7*best_room.d[1][yd] + 0.3*best_room.d[1][!yd]));
	}
	valid_area.clamp_pt_xy(access_pos);
	interior->attic_access.set_from_point(access_pos);
	interior->attic_access.expand_in_dim( long_dim, half_len); // long dim
	interior->attic_access.expand_in_dim(!long_dim, half_wid); // short dim
	copy_zvals(interior->attic_access, ceiling); // same zvals as ceiling
	bool const dir(best_room.get_center_dim(long_dim) < interior->attic_access.get_center_dim(long_dim));
	interior->attic_access.ix = 2*long_dim + dir;
	return 1;
}

cube_t building_t::get_attic_access_door_avoid() const {
	assert(has_attic());
	float const floor_spacing(get_window_vspace());
	cube_with_ix_t avoid(interior->attic_access);
	bool const dim(avoid.ix >> 1), dir(avoid.ix & 1);
	avoid.expand_by_xy(0.25*floor_spacing);
	avoid.d[dim][dir] += (dir ? 1.0 : -1.0)*0.5*floor_spacing; // more spacing in front where the ladder is
	avoid.z2() += 0.5*floor_spacing; // make it taller
	return avoid;
}

void building_t::assign_attic_type(rand_gen_t rgen) {
	if (rgen.rand_bool()) {interior->attic_type = ATTIC_TYPE_RAFTERS; return;} // rafters is the most common case
	// wood and plaster don't look good with hipped roofs because they have vertical beams
	bool const allow_no_rafters(roof_type != ROOF_TYPE_HIPPED && !rgen.rand_bool()); // make these cases rare
	if (allow_no_rafters) {interior->attic_type = rgen.rand()%NUM_ATTIC_TYPES;} // ATTIC_TYPE_RAFTERS, ATTIC_TYPE_FIBERGLASS, ATTIC_TYPE_WOOD, ATTIC_TYPE_PLASTER
	else                  {interior->attic_type = rgen.rand()%2;} // ATTIC_TYPE_RAFTERS, ATTIC_TYPE_FIBERGLASS
}

void find_roofline_beam_span(cube_t &beam, float roof_z2, point const pts[4], bool dim) {
	swap(beam.d[!dim][0], beam.d[!dim][1]); // start denormalized

	for (unsigned n = 0; n < 4; ++n) { // find the span of the top of the roofline
		if (pts[n].z != roof_z2) continue; // point not at peak of roof
		min_eq(beam.d[!dim][0], pts[n][!dim]);
		max_eq(beam.d[!dim][1], pts[n][!dim]);
	}
}
void create_attic_posts(building_t const &b, cube_t const &beam, bool dim, cube_t posts[2]) {
	assert(beam.is_strictly_normalized());
	cube_t const avoid(b.get_attic_access_door_avoid());

	for (unsigned d = 0; d < 2; ++d) {
		cube_t post(beam);
		set_cube_zvals(post, b.interior->attic_access.z2(), beam.z1()); // extends from attic floor to bottom of beam
		post.d[!dim][!d] = post.d[!dim][d] + (d ? -1.0 : 1.0)*beam.dz();
		assert(post.is_strictly_normalized());
		if (!post.intersects_xy(avoid)) {posts[d] = post;} // skip if too close to attic access door
	} // for d
}

unsigned count_num_poly_intersects(float x, float y, float radius, point const *const points, int npts) {
	float const dx[5] = {0.0, 1.0, -1.0, 0.0, 0.0}, dy[5] = {0.0, 0.0, 0.0, 1.0, -1.0};
	unsigned count(0);
	for (unsigned n = 0; n < 5; ++n) {count += point_in_polygon_2d((x + radius*dx[n]), (y + radius*dy[n]), points, npts);} // check 2D XY point contain
	return count;
}

unsigned building_t::get_attic_room_id() const {
	assert(has_attic());
	int const room_id(get_room_containing_pt(cube_bot_center(interior->attic_access) - get_floor_thickness()*plus_z)); // should we cache this during floorplanning?
	assert(room_id >= 0); // must be found
	return room_id;
}
void building_t::add_attic_objects(rand_gen_t rgen) {
	assign_attic_type(rgen); // must be done after roof is added, not in add_attic_access_door()
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size()), obj_flags(RO_FLAG_INTERIOR | RO_FLAG_IN_ATTIC);
	// add attic access door
	cube_with_ix_t adoor(interior->attic_access);
	assert(adoor.is_strictly_normalized());
	adoor.expand_in_z(-0.2*adoor.dz()); // shrink in z
	unsigned const room_id(get_attic_room_id());
	room_t const &room(get_room(room_id));
	bool const ddim(adoor.ix >> 1), ddir(adoor.ix & 1);
	unsigned const acc_flags(room.is_hallway ? RO_FLAG_IN_HALLWAY : 0), attic_door_ix(objs.size());
	float const light_amt(1.0); // always set to 1.0 here, since indir is special cased for attics
	// is light_amount=1.0 correct? since this door can be viewed from both inside and outside the attic, a single number doesn't really work anyway
	objs.emplace_back(adoor, TYPE_ATTIC_DOOR, room_id, ddim, ddir, acc_flags, light_amt, SHAPE_CUBE); // Note: player collides with open attic door
	cube_t const avoid(get_attic_access_door_avoid());
	vect_cube_t avoid_cubes;
	avoid_cubes.push_back(avoid);

	// add light(s)
	cube_t const &part(get_part_for_room(room)); // Note: assumes attic is a single part
	bool const long_dim(part.dx() < part.dy());
	float const floor_spacing(get_window_vspace()), beam_depth(get_attic_beam_depth()), z_floor(interior->attic_access.z2());
	float const sep_dist(part.get_sz_dim(long_dim) - part.get_sz_dim(!long_dim)), attic_height(interior_z2 - z_floor), light_radius(0.03*attic_height);
	point const light_center(part.xc(), part.yc(), (interior_z2 - 1.2*light_radius - beam_depth)); // center of the part near the ceiling
	cube_t light;
	point light_pos[2] = {light_center, light_center}; // start centered
	unsigned num_lights(1);

	if (sep_dist > 0.25*attic_height) { // consider adding two lights
		float const move_dist(0.5*sep_dist - light_radius - beam_depth); // allow extra space for vertical beams
		bool valid(1);

		for (unsigned d = 0; d < 2; ++d) {
			point test_pt(light_center); // spread apart/up an extra radius so that light doesn't partially clip through roof
			test_pt.z += light_radius;
			test_pt[ long_dim] += (d ? -1.0 : 1.0)*(move_dist + light_radius);
			test_pt[!long_dim] += 0.01*sep_dist; // move a tiny bit to the side to avoid incorrect results for queries lying exactly between two roof tquads
			if (!point_in_attic(test_pt)) {valid = 0; break;} // light is outside attic; must be due to hipped roof
		}
		if (valid) {
			light_pos[0][long_dim] -= move_dist;
			light_pos[1][long_dim] += move_dist;
			num_lights = 2;
		}
	}
	for (unsigned n = 0; n < num_lights; ++n) {
		light.set_from_sphere(light_pos[n], light_radius);
		unsigned const light_flags(RO_FLAG_EMISSIVE | RO_FLAG_NOCOLL | obj_flags); // start off and auto turn on when the player enters the attic
		objs.emplace_back(light, TYPE_LIGHT, room_id, 0, 0, light_flags, light_amt, SHAPE_SPHERE, get_light_color_temp(0.45)); // yellow-shite
	}
	if (has_chimney == 1) { // interior chimney; not drawn when player is in the attic because it's part of the exterior geometry
		cube_t chimney(get_chimney());

		if (part.intersects(chimney)) { // in the correct part
			max_eq(chimney.z1(), z_floor);
			min_eq(chimney.z2(), interior_z2); // clip to attic interior range
		
			if (chimney.z1() < chimney.z2()) { // skip if too short; shouldn't happen?
				chimney.expand_by_xy(-0.05*min(chimney.dx(), chimney.dy())); // shrink to make it inside the exterior chimney so that it doesn't show through when outside the attic

				if (!chimney.intersects(avoid)) { // don't block attic access door (probably won't/can't happen)
					objs.emplace_back(chimney, TYPE_CHIMNEY, room_id, 0, 0, obj_flags, light_amt, SHAPE_CUBE, side_color);
					avoid_cubes.push_back(chimney);
				}
			}
		}
	}
	// add vertical vent pipes as avoid cubes
	for (auto i = objs.begin(); i != objs.begin()+objs_start; ++i) {
		if (i->type == TYPE_PIPE && !i->is_exterior() && i->z2() > z_floor && i->intersects_xy(part)) {avoid_cubes.push_back(*i);}
	}
	// add posts as colliders; somewhat of a duplicate of the code in building_room_geom_t::add_attic_rafters()
	for (tquad_with_ix_t const &tq : roof_tquads) {
		if (tq.npts == 3 || !is_attic_roof(tq, 1)) continue; // not a roof tquad; type_roof_only=1
		vector3d const normal(tq.get_norm()); // points outside of the attic
		bool const dim(fabs(normal.x) < fabs(normal.y)), dir(normal[dim] > 0.0); // dim this tquad is facing; beams run in the other dim
		if (dir == 1) continue; // only need to add for one side due to symmetry
		float const beam_width(0.5*beam_depth);
		cube_t const bcube(tq.get_bcube());
		cube_t beam(bcube); // set the z1 base and exterior edge d[dim][dir]
		beam.z1() = beam.z2() - beam_depth; // approximate
		set_wall_width(beam, bcube.d[dim][!dir], 0.5*beam_width, dim); // inside/middle edge
		find_roofline_beam_span(beam, bcube.z2(), tq.pts, dim);
		if (beam.d[!dim][0] == bcube.d[!dim][0]) continue; // not a hipped roof
		if (beam.get_sz_dim(!dim) <= beam_depth) continue; // if it's long enough
		cube_t posts[2];
		create_attic_posts(*this, beam, dim, posts);
		point const access_center(adoor.get_cube_center());
		bool const ls_post(p2p_dist_xy_sq(posts[1].get_cube_center(), access_center) < p2p_dist_xy_sq(posts[0].get_cube_center(), access_center));

		for (unsigned d = 0; d < 2; ++d) {
			cube_t const &post(posts[d]);
			if (post.is_all_zeros()) continue;
			objs.emplace_back(post, TYPE_COLLIDER, room_id, dim, 0, (RO_FLAG_INVIS | obj_flags), 1.0);
			avoid_cubes.push_back(post);
			avoid_cubes.back().expand_by_xy(beam_width); // add extra spacing

			if (bool(d) == ls_post) { // add light switch on the side facing the attic access door for the closer post
				point const post_center(post.get_cube_center());
				vector3d const post_dir(post_center - access_center);
				bool const ls_dim(fabs(post_dir.x) < fabs(post_dir.y)), ls_dir(post_dir[ls_dim] > 0.0);
				cube_t c(get_light_switch_bounds(post.z1(), post.d[ls_dim][!ls_dir], post.get_center_dim(!ls_dim), ls_dim, ls_dir));
				objs.emplace_back(c, TYPE_SWITCH, room_id, ls_dim, ls_dir, (RO_FLAG_NOCOLL | obj_flags), 1.0); // fully lit
			}
		} // for d
	} // for i
	cube_t place_area(part);
	place_area.z1() = place_area.z2() = z_floor; // bottom of attic floor
	place_area.expand_by_xy(-0.75*floor_spacing); // keep away from corners; just a guess; applies to boxes and furnace
	bool has_furnace(0);

	if (interior->furnace_type == FTYPE_ATTIC) { // add furnace in the attic
		// place the furnace above the room where the air intake should be (hallway or stairs room), if there is one
		int const furnace_room(choose_air_intake_room());
		cube_t furnace_place_area(place_area);

		if (furnace_room >= 0) {
			furnace_place_area = get_room(furnace_room);
			furnace_place_area.z1() = place_area.z1();
			furnace_place_area.intersect_with_cube_xy(place_area); // must also be contained within the full place area
		}
		for (unsigned n = 0; n < 100; ++n) { // 100 tries
			if (n == 50) {furnace_place_area = place_area;} // revert to full place area if we haven't found a candidate after 50 iterations
			cube_t furnace;
			bool dim(0), dir(0);
			if (!gen_furnace_cand(furnace_place_area, floor_spacing, 0, rgen, furnace, dim, dir)) break; // near_wall=0
			float const width(furnace.get_sz_dim(!dim)), depth(furnace.get_sz_dim(dim));
			cube_t test_cube(furnace);
			test_cube.d[dim][dir] += (dir ? 1.0 : -1.0)*0.5*depth; // add clearance in front
			if (has_bcube_int(test_cube, avoid_cubes) || !cube_in_attic(furnace)) continue;
			unsigned const flags((is_house ? RO_FLAG_IS_HOUSE : 0) | obj_flags);
			room_object_t const furnace_obj(furnace, TYPE_FURNACE, room_id, dim, dir, flags, light_amt);
			objs.push_back(furnace_obj);
			add_attic_ductwork(rgen, furnace_obj, avoid_cubes);
			avoid_cubes.push_back(test_cube);

			if (furnace_room >= 0) { // place an intake vent in the ceiling of this room under the furnace
				cube_t vent(furnace);
				vent.expand_by_xy(-0.05*width); // shrink slightly
				vent.z2() = furnace.z1() - get_fc_thickness(); // ceiling of the room below
				vent.z1() = vent.z2() - 0.1*get_wall_thickness();
				bool place_ok(1);
				auto use_this_slot(objs.end());

				// vent should be guaranteed to fit inside the room, but may intersect the light; skip it in this case
				for (auto i = objs.begin(); i != objs.end(); ++i) {
					if (i->type == TYPE_LIGHT && i->intersects(vent)) {place_ok = 0; break;}
					if (i->type == TYPE_VENT  && i->intersects(vent)) {use_this_slot = i;} // intersects another room vent
				}
				if (place_ok) {
					if (use_this_slot != objs.end()) {use_this_slot->copy_from(vent); use_this_slot->dim = dim;} // replace ceiling vent with air return vent
					else {objs.emplace_back(vent, TYPE_VENT, furnace_room, dim, 0, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), 1.0);} // dir=0/horizontal; fully lit
				}
			}
			// add an exhaust vent up through the roof; it may clip through the rafters, but this is difficult to check for and avoid
			float const vent_radius(0.1*width), dir_offset(0.26);
			point vent_bot_center(cube_top_center(furnace));
			vent_bot_center[dim] = (dir_offset*furnace.d[dim][!dir] + (1.0 - dir_offset)*furnace.d[dim][dir]);
			add_attic_roof_vent(vent_bot_center, vent_radius, room_id, light_amt);
			has_furnace = 1;
			break; // success/done
		} // for n
	}
	unsigned const rug_avoid_cubes_end(avoid_cubes.size());
	vector3d sz;

	// add lamp(s)
	unsigned const num_lamps(rgen.rand() % (has_furnace ? 5 : 3)); // 0-4/2
	
	for (unsigned n = 0; n < num_lamps; ++n) {
		if (!try_add_lamp(place_area, floor_spacing, room_id, obj_flags, light_amt, avoid_cubes, objs, rgen)) continue;
		if (!cube_in_attic(objs.back())) {objs.pop_back();} // too tall, skip
	}
	// add chair(s)
	unsigned const num_chairs(rgen.rand() % (has_furnace ? 4 : 3)); // 0-3/2
	
	if (num_chairs > 0) {
		float const height(0.4*floor_spacing), hwidth(0.1*floor_spacing);
		sz.assign(hwidth, hwidth, height);
		colorRGBA chair_color(WHITE); // defaults to white

		for (auto i = objs.begin(); i != objs.begin()+objs_start; ++i) {
			if (i->type == TYPE_CHAIR) {chair_color = i->color; break;} // use the color of the first chair added to this building
		}
		for (unsigned n = 0; n < num_chairs; ++n) {
			if (!add_obj_to_closet(objs[attic_door_ix], place_area, objs, avoid_cubes, rgen, sz, TYPE_CHAIR, obj_flags)) continue;
			if (!cube_in_attic(objs.back())) {objs.pop_back(); continue;} // too tall, skip
			objs.back().dim   = rgen.rand_bool(); // random orient
			objs.back().dir   = rgen.rand_bool();
			objs.back().color = chair_color;
		}
	}
	// add nightstand(s)
	unsigned const num_nightstands(rgen.rand() % (has_furnace ? 4 : 3)); // 0-3/2
	
	for (unsigned n = 0; n < num_nightstands; ++n) {
		bool const dim(rgen.rand_bool());
		float const height(rgen.rand_uniform(0.24, 0.26)*floor_spacing), depth(rgen.rand_uniform(0.15, 0.2)*floor_spacing), width(rgen.rand_uniform(1.0, 2.0)*depth);
		sz[ dim] = 0.5*depth;
		sz[!dim] = 0.5*width;
		sz.z     = height;
		if (!add_obj_to_closet(objs[attic_door_ix], place_area, objs, avoid_cubes, rgen, sz, TYPE_NIGHTSTAND, obj_flags)) continue;
		if (!cube_in_attic(objs.back())) {objs.pop_back(); continue;} // too tall, skip
		objs.back().dim = dim;
		objs.back().dir = (objs.back().get_center_dim(dim) < place_area.get_center_dim(dim)); // face the center of the attic so that drawers can be opened
	} // for n
	// add paintcan(s)
	unsigned const num_paintcans(rgen.rand() % (has_furnace ? 6 : 4)); // 0-5/3
	
	if (num_paintcans > 0) {
		float const height(0.64*0.2*floor_spacing), radius(0.28*0.2*floor_spacing);
		sz.assign(radius, radius, height);

		for (unsigned n = 0; n < num_paintcans; ++n) {
			add_obj_to_closet(objs[attic_door_ix], place_area, objs, avoid_cubes, rgen, sz, TYPE_PAINTCAN, obj_flags, SHAPE_CYLIN);
		}
	}
	// add boxes; currently not stacked - should they be?
	unsigned const num_boxes(rgen.rand() % (has_furnace ? 100 : 60)); // 0-99/59
	float const box_sz(0.18*floor_spacing);
	add_boxes_to_space(objs[attic_door_ix], objs, place_area, avoid_cubes, rgen, num_boxes, box_sz, 0.5*box_sz, 1.5*box_sz, 1, obj_flags); // allow_crates=1

	// TYPE_BOOK, TYPE_BOTTLE, TYPE_PAPER?

	// add rug last, under any previous movable items
	point rug_center;
	vector3d rug_hsz; // half length/width
	rug_center.z = z_floor;

	for (unsigned n = 0; n < 20; ++n) { // 20 tries
		for (unsigned d = 0; d < 2; ++d) {rug_hsz   [d] = rgen.rand_uniform(0.2, 0.4)*place_area.get_sz_dim(d);}
		for (unsigned d = 0; d < 2; ++d) {min_eq(rug_hsz[d], 2.0f*rug_hsz[!d]);} // limit aspect ratio to 2.0
		for (unsigned d = 0; d < 2; ++d) {rug_center[d] = rgen.rand_uniform(place_area.d[d][0]+rug_hsz[d], place_area.d[d][1]-rug_hsz[d]);}
		cube_t rug(rug_center);
		rug.expand_by_xy(rug_hsz);
		rug.z2() += get_rug_thickness(); // set thickness/height
		bool bad_place(0);

		for (auto i = avoid_cubes.begin(); i != avoid_cubes.begin()+rug_avoid_cubes_end; ++i) { // skip objects that can be placed on rugs
			if (i->intersects_xy(rug)) {bad_place = 1; break;}
		}
		if (bad_place) continue;
		objs.emplace_back(rug, TYPE_RUG, room_id, 0, 0, (obj_flags | RO_FLAG_NOCOLL), light_amt);
		objs.back().obj_id = uint16_t(objs.size()); // determines rug texture
		break; // done
	} // for n
}

bool building_t::add_attic_roof_vent(point const &bot_center, float radius, unsigned room_id, float light_amt) {
	cube_t pipe(bot_center);
	pipe.expand_by_xy(radius);
	float pipe_z2(bot_center.z); // starting value zero height; should set based on roof height below

	for (auto const &tq : roof_tquads) {
		bool const is_solar_panel(tq.type == tquad_with_ix_t::TYPE_SOLAR);
		if (!is_solar_panel && !is_attic_roof(tq, 1)) continue; // type_roof_only=1
		unsigned const int_count(count_num_poly_intersects(bot_center.x, bot_center.y, radius, tq.pts, tq.npts));
		if (int_count == 0) continue; // no intersections, skip
		if (int_count < 5 || is_solar_panel) return 0; // don't allow vent to clip through solar panel or partially clip through a roof
		vector3d const normal(tq.get_norm());
		if (normal.z == 0.0) continue; // skip vertical sides
		float const denom(dot_product(normal, plus_z));
		if (fabs(denom) < TOLERANCE) break; // should never fail?
		pipe_z2  = bot_center.z + dot_product_ptv(normal, tq.pts[0], bot_center)/denom;
	} // for tq
	if (pipe_z2 <= bot_center.z) return 0; // no roof tquad found; can this happen?
	float const floor_spacing(get_window_vspace());
	pipe_z2 += 0.05*floor_spacing + radius; // move up a bit to get it above the roof
	set_cube_zvals(pipe, bot_center.z, pipe_z2);
	cube_t end_cap(pipe);
	set_cube_zvals(end_cap, pipe_z2, (pipe_z2 + 0.04*floor_spacing + radius));
	end_cap.expand_by_xy(0.8*radius);
	unsigned const end_cap_flags(RO_FLAG_LIT | RO_FLAG_EXTERIOR | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI | RO_FLAG_HANGING | RO_FLAG_NOCOLL); // draw top and bottom flat ends
	interior->room_geom->objs.emplace_back(pipe,    TYPE_PIPE, room_id, 0, 1, RO_FLAG_LIT,   light_amt, SHAPE_CYLIN, LT_GRAY); // dir=1 for vertical; casts shadows
	interior->room_geom->objs.emplace_back(end_cap, TYPE_PIPE, room_id, 0, 1, end_cap_flags, 1.0,       SHAPE_CYLIN, GRAY);
	return 1;
}

cube_t get_attic_access_door_cube(room_object_t const &c, bool inc_ladder) {
	if (!c.is_open()) return c;
	float const len(c.get_length()), thickness(c.dz()), delta(len - thickness);
	cube_t door(c);
	door.z1() -= delta; // open downward
	door.d[c.dim][!c.dir] -= (c.dir ? -1.0 : 1.0)*delta; // shorten to expose the opening
	if (inc_ladder) {door.union_with_cube(get_ladder_bcube_from_open_attic_door(c, door));}
	return door;
}
cube_t get_ladder_bcube_from_open_attic_door(room_object_t const &c, cube_t const &door) {
	float const door_len(c.get_length()), door_width(c.get_width()), door_inside_edge(door.d[c.dim][!c.dir]);
	cube_t ladder(door); // sets ladder step depth
	ladder.expand_in_dim(!c.dim, -0.05*door_width); // a bit narrower
	ladder.d[c.dim][ c.dir] = door_inside_edge; // flush with open side of door
	ladder.d[c.dim][!c.dir] = door_inside_edge + (c.dir ? -1.0 : 1.0)*2.0*c.dz();
	ladder.z1() = door.z2() - 0.95*(door_len/0.44); // matches door length calculation used in floorplanning step
	return ladder;
}

void building_room_geom_t::add_attic_door(room_object_t const &c, float tscale) {
	rgeom_mat_t &wood_mat(get_wood_material(tscale, 1, 0, 1)); // shadows + small
	colorRGBA const color(apply_light_color(c, c.color));
	float cord_z2(c.z1()), cord_len_pos(0.0);

	if (c.is_open()) {
		unsigned const qv_start1(wood_mat.quad_verts.size());
		cube_t const door(get_attic_access_door_cube(c));
		wood_mat.add_cube_to_verts(door, color, door.get_llc(), 0); // all sides
		// rotate 10 degrees
		point rot_pt;
		rot_pt[ c.dim] = door.d[c.dim][!c.dir]; // door inside edge
		rot_pt[!c.dim] = c.get_center_dim(!c.dim); // doesn't matter?
		rot_pt.z       = door.z2(); // top of door
		vector3d const rot_axis(c.dim ? -plus_x : plus_y);
		float const rot_angle((c.dir ? -1.0 : 1.0)*10.0*TO_RADIANS);
		rotate_verts(wood_mat.quad_verts, rot_axis, rot_angle, rot_pt, qv_start1);
		// draw the ladder
		colorRGBA const ladder_color(apply_light_color(c, LT_BROWN)); // slightly darker
		rgeom_mat_t &ladder_mat(get_wood_material(2.0*tscale, 1, 0, 1)); // shadows + small; larger tscale
		unsigned const qv_start2(ladder_mat.quad_verts.size());
		cube_t const ladder(get_ladder_bcube_from_open_attic_door(c, door));
		float const ladder_width(ladder.get_sz_dim(!c.dim));
		float const side_width_factor = 0.05; // relative to door_width

		for (unsigned n = 0; n < 2; ++n) { // sides
			cube_t side(ladder);
			side.d[!c.dim][!n] -= (n ? -1.0 : 1.0)*(1.0 - side_width_factor)*ladder_width;
			ladder_mat.add_cube_to_verts(side, ladder_color, side.get_llc(), EF_Z1, 1); // skip bottom, swap_tex_st=1
		}
		// draw the steps
		unsigned const num_steps = 10;
		float const length(ladder.dz()), step_spacing(length/(num_steps+1)), step_thickness(0.1*step_spacing);
		cube_t step(ladder);
		step.expand_in_dim(!c.dim, -side_width_factor*ladder_width);

		for (unsigned n = 0; n < num_steps; ++n) { // steps
			step.z1() = ladder.z1() + (n+1)*step_spacing;
			step.z2() = step  .z1() + step_thickness;
			ladder_mat.add_cube_to_verts(step, ladder_color, step.get_llc(), get_skip_mask_for_xy(!c.dim), 1); // skip sides, swap_tex_st=1
		}
		rotate_verts(ladder_mat.quad_verts, rot_axis, rot_angle, rot_pt, qv_start2);

		// draw the springs
		for (unsigned n = 0; n < 2; ++n) {
			float const spring_len(0.42*length), radius(0.03*ladder_width), r_wire(0.1*radius), coil_gap(6.0*r_wire);
			point pos;
			pos[ c.dim] = door.d[c.dim][c.dir] + (c.dir ? -1.0 : 1.0)*0.23*c.get_sz_dim(c.dim);
			pos[!c.dim] = door.d[!c.dim][n];
			pos.z = c.z1() - spring_len; // set bottom point
			add_spring(pos, radius, r_wire, spring_len, coil_gap, 2, LT_GRAY); // vertical
		} // for n
		cord_len_pos = -0.35; // closer to the back
		cord_z2     -= 0.36*ladder.dz(); // shift down
	}
	else { // draw only the top and bottom faces of the door
		wood_mat.add_cube_to_verts(c, color, c.get_llc(), ~EF_Z12); // shadows + small, top and bottom only
		cord_len_pos = 0.35; // closer to the front
	}
	// draw the cord and handle
	float const thickness(c.dz()), cord_len(3.0*thickness), cord_radius(0.07*thickness), ball_radius(0.3*thickness);
	cube_t cord;
	set_cube_zvals(cord, cord_z2-cord_len, cord_z2);
	set_wall_width(cord,  c.get_center_dim(!c.dim), cord_radius, !c.dim);
	set_wall_width(cord, (c.get_center_dim( c.dim) + (c.dir ? -1.0 : 1.0)*cord_len_pos*c.get_sz_dim(c.dim)), cord_radius, c.dim);
	rgeom_mat_t &cord_mat(get_untextured_material(0, 0, 1)); // unshadowed, small
	cord_mat.add_vcylin_to_verts(cord, WHITE, 0, 0); // draw sides only
	cord_mat.add_sphere_to_verts(cube_bot_center(cord), ball_radius, WHITE, 1); // low_detail=1
}

bool building_t::is_attic_roof(tquad_with_ix_t const &tq, bool type_roof_only) const {
	if (!has_attic()) return 0;
	if (!tq.is_roof() && (type_roof_only || tq.type != tquad_with_ix_t::TYPE_WALL)) return 0;
	cube_t const tq_bcube(tq.get_bcube());
	if (tq_bcube.z1() < interior->attic_access.z1()) return 0; // not the top section that has the attic (porch roof, lower floor roof)
	return get_attic_part().contains_pt_xy_inclusive(tq_bcube.get_cube_center()); // check for correct part
}

struct edge_t {
	point p[2];
	edge_t() {}
	edge_t(point const &A, point const &B, bool cmp_dim) {
		p[0] = A; p[1] = B;
		if (B[cmp_dim] < A[cmp_dim]) {swap(p[0], p[1]);} // make A less in cmp_dim
	}
};

// from https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
vector3d barycentric_interpolate(vector2d const &p0, vector2d const &p1, vector2d const &p2, vector2d const &P) {
	vector3d ret;
	float const inverse_det(1.0/((p1.y-p2.y) * (p0.x-p2.x) + (p2.x-p1.x) * (p0.y-p2.y)));
	ret.x = ((p1.y-p2.y) * (P.x-p2.x) + (p2.x-p1.x) * (P.y-p2.y)) * inverse_det;
	ret.y = ((p2.y-p0.y) * (P.x-p2.x) + (p0.x-p2.x) * (P.y-p2.y)) * inverse_det;
	ret.z = 1.0 - ret.x - ret.y;
	return ret;
}
void interpolate_tcs(vert_norm_comp_tc_color const &p0, vert_norm_comp_tc_color const &p1, vert_norm_comp_tc_color const &p2,
	vert_norm_comp_tc_color &P, unsigned d0, unsigned d1)
{
	vector3d const weights(barycentric_interpolate(vector2d(p0.v[d0], p0.v[d1]), vector2d(p1.v[d0], p1.v[d1]), vector2d(p2.v[d0], p2.v[d1]), vector2d(P.v[d0], P.v[d1])));
	for (unsigned d = 0; d < 2; ++d) {P.t[d] = weights[0]*p0.t[d] + weights[1]*p1.t[d] + weights[2]*p2.t[d];}
}
void add_attic_roof_geom(rgeom_mat_t &mat, colorRGBA const &color, float thickness, float tscale, bool swap_st, vect_cube_t const &window_holes, building_t const &b) {
	thickness *= b.get_attic_beam_depth();

	for (tquad_with_ix_t const &i : b.roof_tquads) {
		if (!b.is_attic_roof(i, 0)) continue; // type_roof_only=0
		tquad_with_ix_t tq(i);
		std::reverse(tq.pts, tq.pts+tq.npts); // reverse the normal and winding order
		vector3d const normal(tq.get_norm());
		bool const is_roof(tq.is_roof());
		
		for (unsigned n = 0; n < tq.npts; ++n) {
			if (is_roof) {assert(normal.z < 0.0); tq.pts[n].z += thickness/normal.z;} // roof: shift downward
			else {tq.pts[n] += thickness*normal;} // wall: shift inward
		}
		vert_norm_comp_tc_color vert;
		float const denom(0.5f*(b.bcube.dx() + b.bcube.dy())), tsx(tscale/denom), tsy(tscale/denom);
		vert.set_c4(color);
		vert.set_norm(normal);
		unsigned const verts_start(mat.itri_verts.size());

		for (unsigned i = 0; i < tq.npts; ++i) {
			vert.v = tq.pts[i];

			if (is_roof) { // roof
				// shrink bottom edge by thickness to avoid clipping through the floor below;
				// since the quad is a loop, we can check both the previous and next points to see if they're above us, and use the vector delta for the shift
				vector3d shift_dir;
				point const &prev(tq.pts[(i+tq.npts-1)%tq.npts]), &next(tq.pts[(i+1)%tq.npts]);
				if      (vert.v.z < prev.z) {shift_dir = prev - vert.v;}
				else if (vert.v.z < next.z) {shift_dir = next - vert.v;}
				if (shift_dir.z > 0.0) {vert.v += shift_dir*(thickness/shift_dir.z);}
				bool const swap_tc_xy((fabs(normal.x) < fabs(normal.y)) ^ swap_st);
				vert.t[!swap_tc_xy] = vert.v.x*tsx;
				vert.t[ swap_tc_xy] = vert.v.y*tsy;
			}
			else { // side wall
				bool const dim(tq.pts[0].x == tq.pts[1].x); // use nonzero width dim
				vert.t[ swap_st] = vert.v[dim]*tsx;
				vert.t[!swap_st] = vert.v.z*tsy;
			}
			mat.itri_verts.push_back(vert);
		} // for i
		if (!is_roof && !window_holes.empty()) { // cut holes for windows; verts are ordered {top, left, right} viewed from inside
			cube_t const bc(tq.get_bcube());
			bool const dim(bc.dy() < bc.dx()), invert((normal[dim] < 0.0) ^ dim); // should be zero size in X or Y
			float const wall_pos(bc.d[dim][0]);
			cube_t holes_test(bc);
			holes_test.expand_in_dim(dim, thickness);
			bool added_window(0);

			for (cube_t const &w : window_holes) {
				if (!w.intersects(holes_test)) continue;
				point pts[4];
				pts[0].z = pts[1].z = w.z1(); // bottom edge
				pts[2].z = pts[3].z = w.z2(); // top edge
				pts[0][!dim] = pts[3][!dim] = w.d[!dim][ invert];
				pts[1][!dim] = pts[2][!dim] = w.d[!dim][!invert];

				for (unsigned n = 0; n < 4; ++n) { // add 4 quad verts
					pts[n][dim] = wall_pos; // set wall position
					vert.v = pts[n];
					interpolate_tcs(mat.itri_verts[verts_start+0], mat.itri_verts[verts_start+1], mat.itri_verts[verts_start+2], vert, !dim, 2);
					mat.itri_verts.push_back(vert);
				}
				unsigned const ixs[7][3] = {{0,6,5}, {0,5,2}, {5,4,2}, {4,3,2}, {3,1,2}, {6,1,3}, {0,1,6}}; // indices for 7 triangles
				for (unsigned t = 0; t < 7; ++t) {UNROLL_3X(mat.indices.push_back(verts_start + ixs[t][i_]);)}
				added_window = 1;
				break; // there should only be one
			} // for w
			if (added_window) continue; // indices already added, don't need to add below
		}
		for (unsigned i = 0; i < ((tq.npts == 4) ? 6U : 3U); ++i) { // 3 indices for triangles, 6 indices (2 triangles) for quads
			mat.indices.push_back(verts_start + quad_to_tris_ixs[i]);
		}
	} // for i
}

void building_room_geom_t::add_attic_interior_and_rafters(building_t const &b, float tscale, bool detail_pass) {
	if (!b.has_attic()) return;
	if (!detail_pass && !b.has_attic_window) return; // nothing to do
	// determine attic window hole locations
	vect_tquad_with_ix_t window_tquads;
	b.get_attic_windows(window_tquads, 0.0); // offset_scale=0.0
	unsigned const small(window_tquads.empty() ? 2 : 1); // small if there are attic windows, otherwise detail
	if (!detail_pass && small == 2) return; // nothing to do
	unsigned const attic_type(b.interior->attic_type);
	float const window_expand(b.get_wall_thickness());
	float const window_h_border(1.05*b.get_window_h_border()), window_v_border(1.05*b.get_window_v_border()); // (0, 1) range, slightly expanded
	vect_cube_t window_holes;

	for (tquad_with_ix_t const &w : window_tquads) {
		cube_t bc(w.get_bcube());
		bool const dim(bc.dy() < bc.dx()); // should be zero size in X or Y
		bc.expand_in_dim(!dim, -window_h_border*bc.get_sz_dim(!dim));
		bc.expand_in_dim(2,    -window_v_border*bc.dz());
		bc.expand_by(window_expand);
		window_holes.push_back(bc);
	}
	if (detail_pass == (small == 2)) { // draw the attic interior
		if (attic_type == ATTIC_TYPE_WOOD) {
			add_attic_roof_geom(get_material(tid_nm_pair_t(get_plywood_tid()), 0, 0, small), WHITE, 1.0, 16.0, 0, window_holes, b); // no shadows
		}
		else if (attic_type == ATTIC_TYPE_PLASTER) { // or gypsum?
			add_attic_roof_geom(get_material(b.get_material().wall_tex, 0, 0, small), WHITE, 1.0, 16.0, 0, window_holes, b); // no shadows
		}
		else if (attic_type == ATTIC_TYPE_FIBERGLASS) {
			add_attic_roof_geom(get_material(tid_nm_pair_t(get_insulation_tid()), 0, 0, small), colorRGBA(1.0, 0.7, 0.6), 0.5, 16.0, 0, window_holes, b); // no shadows
		}
		else if (attic_type == ATTIC_TYPE_RAFTERS) {
			add_attic_roof_geom(get_material(b.get_attic_texture(), 0, 0, small), WHITE, 0.1, 8.0, 1, window_holes, b); // no shadows, swap_st=1
		}
		else {assert(0);} // unsupported type
	}
	if (!detail_pass || attic_type == ATTIC_TYPE_WOOD || attic_type == ATTIC_TYPE_PLASTER) return; // no rafters to draw

	// build the roof trusses
	get_wood_material(tscale, 0, 0, 2); // ensure unshadowed material
	rgeom_mat_t &wood_mat   (get_wood_material(tscale, 1, 0, 2)); // shadows + detail
	rgeom_mat_t &wood_mat_us(get_wood_material(tscale, 0, 0, 2)); // no shadows + detail
	float const floor_spacing(b.get_window_vspace()), delta_z(0.1*b.get_floor_thickness()); // matches value in get_all_drawn_verts()

	// Note: there may be a chimney in the attic, but for now we ignore it
	for (auto i = b.roof_tquads.begin(); i != b.roof_tquads.end(); ++i) {
		if (!b.is_attic_roof(*i, 0)) continue; // type_roof_only=0
		bool const is_roof(i->is_roof()); // roof tquad; not wall triangle
		// draw beams along inside of roof; start with a vertical cube and rotate to match roof angle
		tquad_with_ix_t tq(*i);
		for (unsigned n = 0; n < tq.npts; ++n) {tq.pts[n].z -= delta_z;} // shift down slightly
		cube_t const bcube(tq.get_bcube());
		vector3d const normal(tq.get_norm()); // points outside of the attic
		bool const dim(fabs(normal.x) < fabs(normal.y)), dir(normal[dim] > 0.0); // dim this tquad is facing; beams run in the other dim
		float const base_width(bcube.get_sz_dim(!dim)), run_len(bcube.get_sz_dim(dim)), height(bcube.dz()), height_scale(1.0/fabs(normal[dim]));
		float const beam_width(0.04*floor_spacing), beam_hwidth(0.5*beam_width), beam_depth(2.0*beam_width);
		float const epsilon(0.02*beam_hwidth), beam_edge_gap(beam_hwidth + epsilon), dir_sign(dir ? -1.0 : 1.0);
		unsigned const num_beams(max(2, round_fp(3.0f*base_width/floor_spacing)));
		float const beam_spacing((base_width - 2.0f*beam_edge_gap)/(num_beams - 1));
		// shift slightly for opposing roof sides to prevent Z-fighting on center beam
		float const beam_pos_start(bcube.d[!dim][0] + beam_edge_gap + dir_sign*0.5*epsilon);
		unsigned const qv_start(wood_mat.quad_verts.size());
		cube_t beam(bcube); // set the z1 base and exterior edge d[dim][dir]
		if (is_roof) {beam.z1() += beam_depth*run_len/height;} // shift up to avoid clipping through the ceiling of the room below
		// determine segments for our non-base edges
		edge_t edges[3]; // non-base edge segments: start plus: 1 for rectangle, 2 for triangle, 3 for trapezoid
		unsigned num_edges(0);

		for (unsigned n = 0; n < tq.npts; ++n) {
			point const &A(tq.pts[n]), &B(tq.pts[(n+1)%tq.npts]);
			if (A.z == bcube.z1() && B.z == bcube.z1()) continue; // base edge, skip
			if (A[!dim] == B[!dim]) continue; // non-angled edge, skip
			edges[num_edges++] = edge_t(A, B, !dim);
		}
		assert(num_edges > 0 && num_edges <= 3);
		bool const is_hipped_side(num_edges == 3 || (num_edges == 2 && tq.npts == 4)); // inc hipped roof for square house with two degenerate quads that should be triangles
		if (is_hipped_side) {assert(i->type == tquad_with_ix_t::TYPE_ROOF_HIP);}
		float const beam_shorten((is_roof ? 2.0 : 1.0)*beam_hwidth*height/(0.5*base_width)); // large for sloped roof to account for width of beams between tquads

		// add vertical beams, which will be rotated to follow the slope of the roof to form rafters
		for (unsigned n = 0; n < num_beams; ++n) {
			float const roof_pos(beam_pos_start + n*beam_spacing);
			set_wall_width(beam, roof_pos, beam_hwidth, !dim);
			beam.d[dim][!dir] = beam.d[dim][dir] + dir_sign*beam_depth;
			bool found(0);

			for (unsigned e = 0; e < num_edges; ++e) {
				edge_t const &E(edges[e]);
				if (roof_pos < E.p[0][!dim] || roof_pos >= E.p[1][!dim]) continue; // beam not contained in this edge
				if (E.p[0].z == E.p[1].z) {beam.z2() = E.p[0].z;} // horizontal edge
				else {beam.z2() = E.p[0].z + ((roof_pos - E.p[0][!dim])/(E.p[1][!dim] - E.p[0][!dim]))*(E.p[1].z - E.p[0].z);} // interpolate zval
				beam.z2() += (height_scale - 1.0)*(beam.z2() - bcube.z1()); // rescale to account for length post-rotate
				beam.z2() -= beam_shorten; // shorten to avoid clipping through the roof at the top
				assert(!found); // break instead?
				found = 1;
			} // for e
			assert(found);
			if (beam.dz() < 4.0f*beam_depth) continue; // too short, skip
			assert(beam.is_strictly_normalized());
			cube_t beam_seg(beam);
			// skip top, bottom and face against the roof (top may be partially visible when rotated)
			unsigned skip_faces(~get_face_mask(dim, dir) | EF_Z12);
			
			if (tq.type == tquad_with_ix_t::TYPE_ROOF_PEAK) { // peaked roof has end rafters against the triangle walls; those faces can be skipped
				if (n == 0          ) {skip_faces |= ~get_face_mask(!dim, 0);}
				if (n == num_beams-1) {skip_faces |= ~get_face_mask(!dim, 1);}
			}
			if (!is_roof) { // triangular wall section; cut out for window holes
				for (cube_t const &w : window_holes) {
					if (!beam_seg.intersects(w)) continue;
					cube_t beam_bot(beam_seg);
					beam_bot.z2() = w.z1();
					beam_seg.z1() = w.z2(); // this becomes the bottom segment
					wood_mat.add_cube_to_verts(beam_bot, WHITE, beam.get_llc(), (skip_faces & ~EF_Z2));
					skip_faces &= ~EF_Z1;
					break; // there should only be one
				}
			}
			wood_mat.add_cube_to_verts(beam_seg, WHITE, beam.get_llc(), skip_faces);
		} // for n
		if (!is_roof) continue; // below is for sloped roof tquads only
		// rotate to match slope of roof
		point rot_pt; // point where roof meets attic floor
		rot_pt[ dim] = bcube.d[dim][dir];
		rot_pt[!dim] = bcube.get_center_dim(dim); // doesn't matter?
		rot_pt.z     = bcube.z1(); // floor
		vector3d const rot_axis(dim ? -plus_x : plus_y);
		float const rot_angle((dir ? 1.0 : -1.0)*atan2(run_len, height));
		rotate_verts(wood_mat.quad_verts, rot_axis, rot_angle, rot_pt, qv_start);

		if (is_hipped_side) { // trapezoid case: add diag beam along both angled edges (hip truss); dim is long dim
			for (unsigned e = 0; e < num_edges; ++e) {
				edge_t const &E(edges[e]);
				if (E.p[0].z == E.p[1].z) continue; // not an angled edge
				bool const low_ix(E.p[1].z == bcube.z1());
				point const &lo(E.p[low_ix]), &hi(E.p[!low_ix]);
				vector3d const edge_delta(hi - lo);
				float const edge_len(edge_delta.mag());
				vector3d const edge_dir(edge_delta/edge_len);
				beam.set_from_point(lo);
				beam.z1() += beam_depth*run_len/height; // avoid clipping through the floor below
				beam.z2() += edge_len; // will be correct after rotation
				beam.expand_in_dim(!dim, beam_hwidth);
				beam.d[dim][!dir] = beam.d[dim][dir] + dir_sign*beam_depth;
				unsigned const qv_start_angled(wood_mat.quad_verts.size());
				wood_mat.add_cube_to_verts(beam, WHITE, beam.get_llc(), (~get_face_mask(dim, dir) | EF_Z12));
				// rotate into place
				vector3d const axis(cross_product(edge_dir, plus_z));
				float const angle(get_angle(plus_z, edge_dir));
				rotate_verts(wood_mat.quad_verts, axis, angle, lo, qv_start_angled);
				// rotate around edge_dir so that bottom surface is aligned with the average normal of the two meeting roof tquads; always 45 degrees
				rotate_verts(wood_mat.quad_verts, edge_dir*((e == num_edges-1) ? 1.0 : -1.0), 0.25*PI, lo, qv_start_angled);
				float const shift_down_val(beam_hwidth*height/run_len);
				for (auto v = wood_mat.quad_verts.begin() + qv_start_angled; v != wood_mat.quad_verts.end(); ++v) {v->v.z -= shift_down_val;}
			} // for e
		}
		if (tq.npts == 4 && dir == 0) {
			// add rafter along the roofline for this quad; dim is long dim
			float const centerline(bcube.d[dim][!dir]); // inside/middle edge
			beam = bcube;
			beam.z2() -= beam_hwidth*height/run_len; // shift to just touching the roof at the top
			beam.z1()  = beam.z2() - beam_depth;
			set_wall_width(beam, centerline, beam_hwidth, dim);
			assert(beam.is_strictly_normalized());
			// determine span of beam along roofline for the trapezoid case/hipped roof; can be zero length for square houses with four triangle roof sections
			if (is_hipped_side) {find_roofline_beam_span(beam, bcube.z2(), tq.pts, dim);}
			beam.expand_in_dim(!dim, -epsilon); // prevent Z-fighting
			
			if (beam.get_sz_dim(!dim) > beam_depth) { // if it's long enough
				unsigned skip_faces(EF_Z2); // skip top
				if (tq.type == tquad_with_ix_t::TYPE_ROOF_PEAK) {skip_faces |= get_skip_mask_for_xy(!dim);} // don't draw ends for peaked roofs as they abut triangle walls
				wood_mat_us.add_cube_to_verts(beam, WHITE, beam.get_llc(), skip_faces); // shadows not needed
				
				if (is_hipped_side) { // trapezoid: add vertical posts (king posts) at each end if there's space
					cube_t posts[2];
					create_attic_posts(b, beam, dim, posts);
				
					for (unsigned d = 0; d < 2; ++d) {
						if (!posts[d].is_all_zeros()) {wood_mat.add_cube_to_verts(posts[d], WHITE, posts[d].get_llc(), EF_Z12);} // skip top and bottom
					}
				}
			}
			if (num_edges == 1) { // tilted rectangle (not trapezoid)
				// add horizontal straining beams connecting each vertical beam to form an A-frame; make them unshadowed because shadows look bad when too close to the light
				beam.z2() -= 3.0*beam_depth; // below roofline beam
				beam.z1()  = beam.z2() - 0.8*beam_depth; // slightly smaller
				float const beam_hlen(((bcube.z2() - beam.z2())/bcube.dz())*run_len); // width of roof tquad at top of beam
				set_wall_width(beam, centerline, beam_hlen, dim);

				for (unsigned n = 1; n+1 < num_beams; ++n) { // same loop as above, but skip the ends
					float const roof_pos(beam_pos_start + n*beam_spacing);
					set_wall_width(beam, roof_pos, 0.9*beam_hwidth, !dim); // slightly thinner to avoid Z-fighting
					wood_mat_us.add_cube_to_verts(beam, WHITE, beam.get_llc(), get_skip_mask_for_xy(dim));
				}
			}
		}
	} // for i
}

// 0=not in attic, 1=in attic with clearance, 2=in attic without clearance
int building_t::vent_in_attic_test(cube_t const &vent, bool dim) const {
	point test_pt(cube_top_center(vent) + 1.1*get_fc_thickness()*plus_z);
	if (!point_in_attic(test_pt)) return 0; // not in attic
	test_pt.z += vent.get_sz_dim(dim) + get_attic_beam_depth(); // check duct height + beam clearance
	return (point_in_attic(test_pt) ? 1 : 2);
}

bool duct_merges_to_xy(cube_t const &from, cube_t const &to) { // Note: assumes <from> and <to> have the same z-ranges
	if (!from.intersects_xy(to)) return 0;
	return ((from.x1() >= to.x1() && from.x2() <= to.x2()) || (from.y1() >= to.y1() && from.y2() <= to.y2())); // check either X and Y dims for containment
}

bool maybe_clip_overlapping_duct(room_object_t &duct, vect_room_object_t const &objs, unsigned objs_start, bool is_cylin, vect_cube_t &sub_cubes) {
	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (!i->intersects(duct)) continue;
		subtract_cube_from_cube((cube_t)duct, *i, sub_cubes, 1); // clear_out=1
		if (sub_cubes.empty())    return 0; // contained? shouldn't happen
		if (sub_cubes.size() > 1) continue; // multiple parts, skip because this case is too complex; should be rare
		duct.copy_from(sub_cubes[0]); // single part - clip to shorter length to remove the overlap

		if (is_cylin && duct.dim != i->dim) { // right angle cylinder intersection, need to extend to merge cylinders
			float const extend(0.5*duct.get_width());
			if      (duct.d[duct.dim][0] == i->d[duct.dim][1]) {duct.d[duct.dim][0] -= extend;} // lo side
			else if (duct.d[duct.dim][1] == i->d[duct.dim][0]) {duct.d[duct.dim][1] += extend;} // hi side
		}
	} // for i
	return 1;
}

bool try_route_duct_with_jog(room_object_t &duct, cube_t const &dest, bool first_dim, vect_room_object_t &objs, unsigned objs_start,
	vect_cube_t const &avoid_cubes, vect_cube_t &sub_cubes, float extend_len, bool &use_extend, bool is_cylin)
{
	bool const dirs[2] = {(duct.xc() < dest.xc()), (duct.yc() < dest.yc())};

	for (unsigned n = 0; n < 2; ++n) { // try both routing directions
		bool const dim(first_dim ^ bool(n)), dir1(dirs[dim]), dir2(!dirs[!dim]);
		room_object_t cand1(duct), cand2(duct);
		cand1.dim  =  dim;
		cand2.dim  = !dim;
		cand2.x1() = dest.x1(); cand2.y1() = dest.y1(); cand2.x2() = dest.x2(); cand2.y2() = dest.y2();

		if (is_cylin) { // segments meet at their centerlines with a bend in between
			cand1.d[ dim][ dir1] = dest.get_center_dim( dim); // extend to the destination centerline
			cand2.d[!dim][ dir2] = duct.get_center_dim(!dim); // extend to the duct centerline
		}
		else { // cand2 ends flush with cand1; cand2 may overlap with dest, but it will be clipped later
			cand1.d[ dim][ dir1] = dest.d[ dim][!dir1]; // extend to the destination near side; may be denormalized and skipped below
			cand2.d[!dim][ dir2] = duct.d[!dim][ dir2]; // extend to the duct
		}
		assert(cand2.is_strictly_normalized());
		bool use_cand1(cand1.is_strictly_normalized()), use_cand2(1);
		if ((use_cand1 && has_bcube_int(cand1, avoid_cubes)) || has_bcube_int(cand2, avoid_cubes)) continue; // bad routing

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) { // merge to previously placed ducts
			// Note: if our cand1 intersects an earlier and shorter cand1 duct, we want to either clip our cand1 to the existing cand1,
			// or extend the existing cand1 to meet our duct; however, sorting from furthest to nearest above should prevent this from happening
			// in most cases, and that's the point of the sort, to avoid having to handle this case
			if (i->contains_cube(cand2)) {use_cand2 = 0; break;} // cand2 can be skipped
		}
		if (use_cand1) {use_cand1 = maybe_clip_overlapping_duct(cand1, objs, objs_start, is_cylin, sub_cubes);}
		if (use_cand2) {use_cand2 = maybe_clip_overlapping_duct(cand2, objs, objs_start, is_cylin, sub_cubes);}

		if (use_extend) { // check if our route is shorter or longer than the extend length and use whichever is shorter
			float tot_len(0.0);
			if (use_cand1) {tot_len += cand1.get_sz_dim( dim);}
			if (use_cand2) {tot_len += cand2.get_sz_dim(!dim);}
			if (tot_len < extend_len) {use_extend = 0;} else {break;}
		}
		if (is_cylin && use_cand1 && use_cand2) {cand1.flags |= (dir1 ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);} // place a bend at this location
		if (use_cand1) {objs.push_back(cand1);} // maybe add first  segment from vent
		if (use_cand2) {objs.push_back(cand2);} // maybe add second segment to furnace
		return 1; // success
	} // for n
	return 0; // failed
}

struct cmp_by_dist_decending { // furthest to closest
	cube_t const ref;
	cmp_by_dist_decending(cube_t const &ref_) : ref(ref_) {}
	
	float get_dist(cube_t const &c) const {
		point const closest_pt(ref.closest_pt(c.get_cube_center()));
		return (fabs(c.xc() - closest_pt.x) + fabs(c.yc() - closest_pt.y)); // Manhattan XY distance
	}
	bool operator()(cube_t const &a, cube_t const &b) const {return (get_dist(b) < get_dist(a));}
};

// should there also be an add_basement_ductwork() for office buildings?
void building_t::add_attic_ductwork(rand_gen_t rgen, room_object_t const &furnace, vect_cube_t &avoid_cubes) {
	assert(has_room_geom());
	bool const is_cylin((rgen.rand()%3) == 0); // 33% of the time
	room_obj_shape const duct_shape(is_cylin ? SHAPE_CYLIN : SHAPE_CUBE);
	float const fc_thick(get_fc_thickness());
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const duct_flags(RO_FLAG_INTERIOR | RO_FLAG_IN_ATTIC | RO_FLAG_ADJ_BOT);

	// attempt to add the intake vent to the top of the furnace if there's room under the roof
	float const intake_height(0.2*furnace.dz()); // slightly wider than vent ducts (0.167*dz)
	cube_t duct_top(furnace);
	narrow_furnace_intake(duct_top, furnace);
	duct_top.z2() += intake_height; // make it taller
	cube_t intake(duct_top);
	duct_top.z1() = furnace.z2();
	intake.d[furnace.dim][ furnace.dir]  = furnace.d[furnace.dim][!furnace.dir]; // back of furnace
	intake.d[furnace.dim][!furnace.dir] -= (furnace.dir ? 1.0 : -1.0)*intake_height; // back of intake

	if (cube_in_attic(duct_top) && cube_in_attic(intake)) {
		cube_t const cubes[2] = {duct_top, intake};
		
		for (unsigned d = 0; d < 2; ++d) {
			objs.emplace_back(cubes[d], TYPE_DUCT, furnace.room_id, furnace.dim, 1, duct_flags, 1.0, SHAPE_CUBE, DUCT_COLOR); // dir=1; vertical, always cubes
		}
	}
	// add ducts on the attic floor
	unsigned const ducts_start(objs.size()); // includes cubes above vents in cylindrical ducts case
	vect_room_object_t ducts, objs_to_add;

	// find all vents in the ceilings of upper floors just below the attic
	for (auto i = objs.begin(); i != objs.end(); ++i) { // Note: can't use get_placed_objs_end() because buttons_start hasn't been set yet
		if (i->type != TYPE_VENT) continue;
		if (vent_in_attic_test(*i, i->dim) != 1) continue; // check for attic floor/top ceiling and for roof clearance
		float const duct_width(i->get_width()), duct_height((is_cylin ? 0.75 : 0.72)*duct_width); // smaller than opening; allows player to walk over duct more easily
		cube_t duct(*i);
		duct.z1() = i->z2() + fc_thick; // attic floor
		duct.z2() = duct.z1() + duct_height;
		// copy room_id from the vent, even though it's not in the same room as the vent; add to ducts rather than objs to avoid iterator invalidation
		room_object_t duct_obj(duct, TYPE_DUCT, i->room_id, 0, 0, duct_flags, 1.0, duct_shape, DUCT_COLOR);

		if (is_cylin) { // add a rectangular box that the cylinder connects to
			objs_to_add.push_back(duct_obj);
			objs_to_add.back().shape = SHAPE_CUBE;

			for (unsigned d = 0; d < 2; ++d) { // shrink the width to equal the height to make it square
				duct_obj.expand_in_dim(d, 0.5*(duct_height - i->get_sz_dim(d)));
			}
		}
		ducts.push_back(duct_obj);
	} // for i
	vector_add_to(objs_to_add, objs); // add cylindrical duct boxes
	// sort ducts furthest to closest to the furnace so that shorter runs can be connected to existing longer runs
	sort(ducts.begin(), ducts.end(), cmp_by_dist_decending(furnace));
	vect_cube_t sub_cubes; // reused
	vect_room_object_t ducts_to_reroute;
	unsigned const horiz_ducts_start(objs.size());

	for (room_object_t &duct : ducts) {
		bool added(0);

		for (auto i = objs.begin()+horiz_ducts_start; i != objs.end(); ++i) {
			if (!duct_merges_to_xy(duct, *i)) continue;
			subtract_cube_from_cube((cube_t)duct, *i, sub_cubes, 1); // clear_out=1
			if (sub_cubes.size() == 1) {duct.copy_from(sub_cubes[0]);} // clip to premove overlaps; should generally (always?) be true
			objs.push_back(duct); // add vertical duct with no extension
			added = 1;
			break;
		}
		if (added) continue; // done with this duct
		// find shortest straight line route to nearest existing duct, record length, and use if below path fails or has a longer length
		room_object_t extend_duct;
		float extend_len(0.0);
		bool use_extend(0), was_connected(0);

		for (auto i = objs.begin()+horiz_ducts_start; i != objs.end(); ++i) {
			for (unsigned d = 0; d < 2; ++d) { // x,y
				if (i->d[!d][0] > duct.d[!d][0] || i->d[!d][1] < duct.d[!d][1]) continue; // duct not contained in the other dim
				bool const dir(duct.get_center_dim(d) < i->get_center_dim(d));
				room_object_t cand(duct);
				cand.d[d][dir] = i->d[d][!dir]; // extend to meet existing duct
				if (is_cylin) {cand.d[d][dir] += (dir ? 1.0 : -1.0)*0.5*i->get_sz_dim(d);} // extend to centerline
				assert(cand.is_strictly_normalized()); // is this too strong? continue in this case? or does the merge check above always set added=1 in this case?
				float const len(cand.get_sz_dim(d));
				if (extend_len > 0.0 && len > extend_len) continue; // longer than prev best connection
				if (has_bcube_int(cand, avoid_cubes))     continue;
				extend_duct = cand; extend_len = len; extend_duct.dim = d; use_extend = 1;
			} // for d
		} // for i
		float const duct_width(duct.get_length()), seg2_width_x(min(duct_width, 0.95f*furnace.dx())), seg2_width_y(min(duct_width, 0.95f*furnace.dy()));
		cube_t port(furnace);
		port.expand_by(vector3d(-0.5*(furnace.dx() - seg2_width_x), -0.5*(furnace.dy() - seg2_width_y), 0.0)); // shrink in X and Y
		bool const first_dim(furnace.dim); // make it consistent across vents to maximize sharing
		was_connected = try_route_duct_with_jog(duct, port, first_dim, objs, horiz_ducts_start, avoid_cubes, sub_cubes, extend_len, use_extend, is_cylin);
		if (use_extend) {objs.push_back(extend_duct);} // use straight extension
		else if (!was_connected) {ducts_to_reroute.push_back(duct);} // likely blocked by the attic door between the vent and the furnace
	} // for ducts
	for (room_object_t &duct : ducts_to_reroute) {
		// try one jog to an existing duct
		vect_room_object_t conns;

		for (auto i = objs.begin()+horiz_ducts_start; i != objs.end(); ++i) {
			// split into a number of candidate connection points along the length of this duct segment
			float const length(i->get_length()), width(i->get_width());
			unsigned const num_steps(round_fp(0.25f*length/width));
			if (num_steps == 0) continue; // too short, skip
			float const step((length - width)/num_steps);

			for (unsigned n = 0; n < num_steps; ++n) {
				room_object_t conn(*i);
				conn.d[i->dim][0] += n*step; // left edge
				conn.d[i->dim][1]  = conn.d[i->dim][0] + width; // right edge
				conns.push_back(conn);
			}
		} // for i
		sort(conns.begin(), conns.end(), cmp_by_dist_decending(duct)); // reuse our convenient sort function
		reverse(conns.begin(), conns.end()); // here we want closest first

		for (room_object_t const &conn : conns) {
			bool use_extend(0); // will remain at 0
			if (try_route_duct_with_jog(duct, conn, conn.dim, objs, horiz_ducts_start, avoid_cubes, sub_cubes, 0.0, use_extend, is_cylin)) {break;}
		}
	} // for ducts_to_reroute
	for (auto i = objs.begin()+ducts_start; i != objs.end(); ++i) {avoid_cubes.push_back(*i);} // add *all* ducts to avoid_cubes
}

int building_t::choose_air_intake_room() const { // for the air return
	int best_room(-1); // start at not found

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (has_attic() && !get_attic_part().contains_cube(*r)) continue; // if there's an attic, the room must be in the same part as it
		if (r->is_hallway) {return (r - interior->rooms.begin());} // hallway is always the best room, return it now
		if (r->has_stairs == 255) {best_room = (r - interior->rooms.begin());} // check for stairs on all floors; return this room unless a hallway is found later
	}
	return best_room;
}

void building_room_geom_t::add_chimney(room_object_t const &c, tid_nm_pair_t const &tex) { // inside attic
	tid_nm_pair_t tex2(tex);
	tex2.shadowed = 1;
	tex2.tscale_x *= 4.0; tex2.tscale_y *= 4.0;
	rgeom_mat_t &mat(get_material(tex2, 1, 0, 2));
	mat.add_cube_to_verts(c, c.color, c.get_llc(), (EF_Z12 | EF_Y12), 1); // X sides, swap_tex_st=1
	mat.add_cube_to_verts(c, c.color, c.get_llc(), (EF_Z12 | EF_X12), 0); // Y sides, swap_tex_st=0
}

// I guess skylights are similar to attics, so this code goes here?
void building_room_geom_t::add_skylights_details(building_t const &b) {
	for (cube_t const &skylight : b.skylights) {add_skylight_details(skylight, b.has_skylight_light);}
}
void add_skylight_bar_grid(cube_t const &skylight, float bar_hwidth, float spacing, unsigned num, bool dim, vect_cube_t &bars) {
	for (unsigned i = 0; i < num+1; ++i) { // includes bar at the end
		cube_t bar(skylight); // uses skylight zvals
		set_wall_width(bar, (skylight.d[dim][0] + bar_hwidth + i*spacing), bar_hwidth, dim);

		if (i == 0 || i == num) { // make the bars at the edges a bit wider and taller
			bar.expand_by_xy(0.25*bar_hwidth);
			bar.z2() += 0.25*bar_hwidth;
			bar.z1() -= 0.02*bar_hwidth; // move down slightly, enough to not Z-fight with the ceiling but not enough to clip through the elevator ceiling
		}
		bars.push_back(bar);
	} // for i
}
void building_room_geom_t::add_skylight_details(cube_t const &skylight, bool has_skylight_light) {
	vector3d const sz(skylight.get_size());
	bool const dmax(sz.x < sz.y);
	float const thickness(skylight.dz()), length(sz[dmax]), width(sz[!dmax]), bar_width(0.67*sz.z), bar_hwidth(0.5*bar_width);
	if (width < 10.0*thickness) return; // no bars needed
	unsigned const num_rows((unsigned)ceil(0.1*length/thickness)), num_cols((unsigned)ceil(0.04*width/thickness));
	float const row_spacing((length - bar_width)/num_rows), col_spacing((width - bar_width)/num_cols);
	vect_cube_t bars[2]; // per dim
	add_skylight_bar_grid(skylight, bar_hwidth, col_spacing, num_cols, !dmax, bars[!dmax]);
	add_skylight_bar_grid(skylight, bar_hwidth, row_spacing, num_rows,  dmax, bars[ dmax]);

	if (has_skylight_light) {
		cube_t test_cube(skylight);
		test_cube.z1() -= skylight.dz(); // shift down to include light z-val range
		auto objs_end(get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = objs.begin(); i != objs_end; ++i) {
			if (i->type != TYPE_LIGHT || !i->intersects(test_cube)) continue;
			cube_t bar_area(*i);
			bar_area.expand_by_xy(-0.25*i->get_size()); // shrink to half length and width
			bool connected(0);

			for (unsigned d = 0; d < 2 && !connected; ++d) {
				for (cube_t const &bar : bars[d]) {
					if (bar.intersects_xy(bar_area)) {connected = 1; break;}
				}
			}
			if (connected) continue; // already connected - done
			cube_t conn;
			set_cube_zvals(conn, skylight.z1(), skylight.zc()); // half height
			set_wall_width(conn, i->get_center_dim(!dmax), 0.8*bar_hwidth, !dmax); // centered on the light

			// find bars to either side in dmax; since they don't intersect the light, and the light is ocntained in the skylight, they must exist
			for (cube_t const &bar : bars[dmax]) {
				float const low_edge(bar.d[dmax][0]);
				if (low_edge <  bar_area.d[dmax][0]) continue; // not yet
				conn.d[dmax][1] = low_edge; // inner bar ends here
				conn.d[dmax][0] = low_edge - (row_spacing - bar_width); // inner bar starts here
				break; // done
			}
			assert(conn.is_strictly_normalized()); // must have been found in the loop above
			bars[!dmax].push_back(conn);
		} // for i
	}
	// Note: if drawn as exterior, lighting may look better, but it won't cast shadows on the interior;
	// can draw only the top as exterior, but that's the face oriented toward the light that casts the shadow
	bool const draw_as_exterior = 0;
	rgeom_mat_t &bar_mat(get_untextured_material(!draw_as_exterior, 0, 0, 0, draw_as_exterior)); // inc_shadows=1

	for (unsigned d = 0; d < 2; ++d) {
		for (cube_t const &bar : bars[d]) {bar_mat.add_cube_to_verts_untextured(bar, WHITE, get_skip_mask_for_xy(!d));} // skip ends
	}
}

