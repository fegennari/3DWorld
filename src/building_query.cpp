// 3D World - Building Query Functions
// by Frank Gennari 10/09/21

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"


extern bool draw_building_interiors, player_near_toilet, player_is_hiding, player_in_elevator;
extern int player_in_closet;
extern double camera_zh;
extern building_params_t global_building_params;
extern bldg_obj_type_t bldg_obj_types[];


float get_railing_height(room_object_t const &c);
cylinder_3dw get_railing_cylinder(room_object_t const &c);
bool sphere_vert_cylin_intersect_with_ends(point &center, float radius, cylinder_3dw const &c, vector3d *cnorm);


bool building_t::check_bcube_overlap_xy(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const {

	if (expand_rel == 0.0 && expand_abs == 0.0 && !bcube.intersects(b.bcube)) return 0;
	if (!is_rotated() && !b.is_rotated()) return 1; // above check is exact, top-level bcube check up to the caller
	if (b.bcube.contains_pt_xy(bcube.get_cube_center()) || bcube.contains_pt_xy(b.bcube.get_cube_center())) return 1; // slightly faster to include this check
	return (check_bcube_overlap_xy_one_dir(b, expand_rel, expand_abs, points) || b.check_bcube_overlap_xy_one_dir(*this, expand_rel, expand_abs, points));
}

// Note: only checks for point (x,y) value contained in one cube/N-gon/cylinder; assumes pt has already been rotated into local coordinate frame
bool building_t::check_part_contains_pt_xy(cube_t const &part, point const &pt, vector<point> &points) const {

	if (!part.contains_pt_xy(pt)) return 0; // check bounding cube
	if (is_simple_cube()) return 1; // that's it
	building_draw_utils::calc_poly_pts(*this, part, points);
	return point_in_polygon_2d(pt.x, pt.y, points.data(), points.size(), 0, 1); // 2D x/y containment
}

bool building_t::check_bcube_overlap_xy_one_dir(building_t const &b, float expand_rel, float expand_abs, vector<point> &points) const { // can be called before levels/splits are created

	// Note: easy cases are handled by check_bcube_overlap_xy() above
	point const center1(b.bcube.get_cube_center()), center2(bcube.get_cube_center());

	for (auto p1 = b.parts.begin(); p1 != b.parts.end(); ++p1) {
		point pts[9]; // {center, 00, 10, 01, 11, x0, x1, y0, y1}

		if (b.parts.size() == 1) {pts[0] = center1;} // single cube: we know we're rotating about its center
		else {
			pts[0] = p1->get_cube_center();
			b.do_xy_rotate(center1, pts[0]); // rotate into global space
		}
		cube_t c_exp(*p1);
		c_exp.expand_by_xy(expand_rel*p1->get_size() + vector3d(expand_abs, expand_abs, expand_abs));

		for (unsigned i = 0; i < 4; ++i) { // {00, 10, 01, 11}
			pts[i+1].assign(c_exp.d[0][i&1], c_exp.d[1][i>>1], 0.0); // XY only
			b.do_xy_rotate(center1, pts[i+1]); // rotate into global space
		}
		for (unsigned i = 0; i < 5; ++i) {do_xy_rotate_inv(center2, pts[i]);} // inverse rotate into local coord space - negate the sine term
		cube_t c_exp_rot(pts+1, 4); // use points 1-4
		pts[5] = 0.5*(pts[1] + pts[3]); // x0 edge center
		pts[6] = 0.5*(pts[2] + pts[4]); // x1 edge center
		pts[7] = 0.5*(pts[1] + pts[2]); // y0 edge center
		pts[8] = 0.5*(pts[3] + pts[4]); // y1 edge center

		for (auto p2 = parts.begin(); p2 != parts.end(); ++p2) {
			if (c_exp_rot.contains_pt_xy(p2->get_cube_center())) return 1; // quick and easy test for heavy overlap

			for (unsigned i = 0; i < 9; ++i) {
				if (check_part_contains_pt_xy(*p2, pts[i], points)) return 1; // Note: building geometry is likely not yet generated, this check should be sufficient
				//if (p2->contains_pt_xy(pts[i])) return 1;
			}
		}
	} // for p1
	return 0;
}

// called for players and butterfiles
bool building_t::test_coll_with_sides(point &pos, point const &p_last, float radius, cube_t const &part, vector<point> &points, vector3d *cnorm) const {

	building_draw_utils::calc_poly_pts(*this, part, points); // without the expand
	point quad_pts[4]; // quads
	unsigned const num_steps(max(1U, min(100U, (unsigned)ceil(2.0*p2p_dist_xy(pos, p_last)/radius))));
	vector3d const step_delta((pos - p_last)/num_steps);
	pos = p_last;

	// if the player is moving too quickly, the intersection with a side polygon may be missed, which allows the player to travel through the building;
	// so we split the test into multiple smaller sphere collision steps
	for (unsigned step = 0; step < num_steps; ++step) {
		pos += step_delta;

		for (unsigned S = 0; S < num_sides; ++S) { // generate vertex data quads
			for (unsigned d = 0, ix = 0; d < 2; ++d) {
				point const &p(points[(S+d)%num_sides]);
				for (unsigned e = 0; e < 2; ++e) {quad_pts[ix++].assign(p.x, p.y, part.d[2][d^e]);}
			}
			vector3d const normal(get_poly_norm(quad_pts));
			float const rdist(dot_product_ptv(normal, pos, quad_pts[0]));
			if (rdist < 0.0 || rdist >= radius) continue; // too far or wrong side
			if (!sphere_poly_intersect(quad_pts, 4, pos, normal, rdist, radius)) continue;
			pos += normal*(radius - rdist);
			if (cnorm) {*cnorm = normal;}
			return 1;
		} // for S
	} // for step
	if (max(pos.z, p_last.z) > part.z2() && point_in_polygon_2d(pos.x, pos.y, points.data(), num_sides, 0, 1)) { // test top plane (sphere on top of polygon?)
		pos.z = part.z2() + radius; // make sure it doesn't intersect the roof
		if (cnorm) {*cnorm = plus_z;}
		return 1;
	}
	return 0;
}

cube_t building_t::get_coll_bcube() const {
	if (!is_house || (!has_ac && !has_driveway())) return bcube;
	cube_t cc(bcube);
	if (has_ac) {cc.expand_by_xy(0.35*get_window_vspace());} // conservative
	if (has_driveway()) {cc.union_with_cube(driveway);}
	return cc;
}

// Note: used for the player
bool building_t::check_sphere_coll(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, vector3d const &xlate,
	float radius, bool xy_only, vector<point> &points, vector3d *cnorm_ptr, bool check_interior) const
{
	if (!is_valid()) return 0; // invalid building
	
	if (radius > 0.0) {
		cube_t const cc(get_coll_bcube());
		if (!(xy_only ? sphere_cube_intersect_xy((pos - xlate), radius, cc) : sphere_cube_intersect((pos - xlate), radius, cc))) return 0;
	}
	float const xy_radius(radius*global_building_params.player_coll_radius_scale);
	point pos2(pos), p_last2(p_last), center;
	bool had_coll(0), is_interior(0);
	float part_z2(bcube.z2());

	if (is_rotated()) {
		center = bcube.get_cube_center() + xlate;
		do_xy_rotate_inv(center, pos2); // inverse rotate - negate the sine term
		do_xy_rotate_inv(center, p_last2);
	}
	if (check_interior && draw_building_interiors && interior != nullptr) { // check for interior case first
		float const zval(max(pos2.z, p_last2.z) - xlate.z); // this is the real zval for use in collsion detection, in building space
		point const pos2_bs(pos2 - xlate);
		cube_t sc; sc.set_from_sphere(pos2_bs, radius); // sphere bounding cube

		// Note: first check uses min of the two zvals to reject the basement, which is actually under the mesh
		if ((min(pos2.z, p_last2.z) + radius) > ground_floor_z1 && zval < (ground_floor_z1 + get_door_height())) { // on the ground floor
			for (auto d = doors.begin(); d != doors.end(); ++d) {
				if (d->type == tquad_with_ix_t::TYPE_RDOOR) continue; // doesn't apply to roof door
				cube_t bc(d->get_bcube());
				bool const door_dim(bc.dy() < bc.dx());
				bc.expand_in_dim( door_dim, 1.1*radius); // expand by radius plus some tolerance in door dim
				bc.expand_in_dim(!door_dim, -0.5*xy_radius); // shrink slightly in the other dim to prevent the player from clipping through the wall next to the door
				bc.z1() -= max(radius, (float)camera_zh); // account for player on steep slope up to door - require player head above doorframe bottom
				if (bc.contains_pt(pos2_bs)) return 0; // check if we can use a door - disable collsion detection to allow the player to walk through
			}
		}
		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // include garages and sheds
			cube_t c(*i + xlate);
			if (!c.contains_pt(point(pos2.x, pos2.y, zval))) continue; // not interior to this part
			float cont_area(0.0);

			for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
				if (zval < p->z1() || zval > p->z2()) continue; // wrong floor/part in stack
				if (p->intersects_xy(sc)) {cont_area += (min(p->x2(), sc.x2()) - max(p->x1(), sc.x1()))*(min(p->y2(), sc.y2()) - max(p->y1(), sc.y1()));} // accumulate shared XY area
			}
			if (cont_area < 0.99*sc.get_area_xy()) { // sphere bounding cube not contained in union of parts - sphere is partially outside the building
				c.expand_by_xy(-radius); // shrink part by sphere radius
				c.clamp_pt_xy(pos2); // force pos2 into interior of the cube to prevent the sphere from intersecting the part
				had_coll = 1;
			}
			is_interior = 1;
			break; // flag for interior collision detection
		} // for i
		if (!is_interior) { // not interior to a part - check roof access
			float const floor_thickness(get_floor_thickness());

			for (auto i = interior->stairwells.begin(); i != interior->stairwells.end(); ++i) {
				if (!i->roof_access) continue;
				cube_t test_cube(*i);
				test_cube.expand_by_xy(-0.5f*radius);
				if (!test_cube.contains_pt_xy(pos2 - xlate)) continue; // pos not over stairs
				if (zval < i->z2() + radius + floor_thickness) {is_interior = 1; break;} // Note: don't have to check zval > i->z2() because we know that !is_interior
			}
		}
	}
	if (is_interior) {
		point pos2_bs(pos2 - xlate);
		if (check_sphere_coll_interior(pos2_bs, (p_last2 - xlate), ped_bcubes, radius, xy_only, cnorm_ptr)) {pos2 = pos2_bs + xlate; had_coll = 1;}
	}
	else {
		for (auto i = parts.begin(); i != parts.end(); ++i) {
			if (xy_only && i->z1() > ground_floor_z1) { // only need to check first level in this mode
				if (has_complex_floorplan) {continue;} else {break;}
			}
			if (!xy_only && ((pos2.z + radius < i->z1() + xlate.z) || (pos2.z - radius > i->z2() + xlate.z))) continue; // test z overlap
			if (radius == 0.0 && !(xy_only ? i->contains_pt_xy(pos2) : i->contains_pt(pos2))) continue; // no intersection; ignores p_last
			cube_t const part_bc(*i + xlate);
			bool part_coll(0);

			if (use_cylinder_coll()) {
				point const cc(part_bc.get_cube_center());
				float const crx(0.5*i->dx()), cry(0.5*i->dy()), r_sum(radius + max(crx, cry));
				if (!dist_xy_less_than(pos2, cc, r_sum)) continue; // no intersection

				if (fabs(crx - cry) < radius) { // close to a circle
					if (p_last2.z > part_bc.z2() && dist_xy_less_than(pos2, cc, max(crx, cry))) {
						pos2.z = part_bc.z2() + radius; // make sure it doesn't intersect the roof
						if (cnorm_ptr) {*cnorm_ptr = plus_z;}
					}
					else { // side coll
						vector2d const d((pos2.x - cc.x), (pos2.y - cc.y));
						float const mult(r_sum/d.mag());
						pos2.x = cc.x + mult*d.x;
						pos2.y = cc.y + mult*d.y;
						if (cnorm_ptr) {*cnorm_ptr = vector3d(d.x, d.y, 0.0).get_norm();} // no z-component
					}
					part_coll = 1;
				}
				else {
					part_coll |= test_coll_with_sides(pos2, p_last2, radius, part_bc, points, cnorm_ptr); // use polygon collision test
				}
			}
			else if (num_sides != 4) { // triangle, hexagon, octagon, etc.
				part_coll |= test_coll_with_sides(pos2, p_last2, radius, part_bc, points, cnorm_ptr);
			}
			else if (!xy_only && part_bc.contains_pt_xy_exp(pos2, radius) && p_last2.z > (i->z2() + xlate.z)) { // on top of building
				pos2.z = i->z2() + xlate.z + radius;
				if (cnorm_ptr) {*cnorm_ptr = plus_z;}
				part_coll = 1;
			}
			else if (sphere_cube_int_update_pos(pos2, radius, part_bc, p_last2, 1, xy_only, cnorm_ptr)) { // cube
				part_coll = 1; // flag as colliding, continue to look for more collisions (inside corners)
			}
			if (part_coll && pos2.z < part_bc.z1()) {pos2.z = part_bc.z2() + radius;} // can't be under a building - make it on top of the building instead
			if (part_coll) {part_z2 = i->z2();}
			had_coll |= part_coll;
		} // for i
		for (auto i = fences.begin(); i != fences.end(); ++i) {
			had_coll |= sphere_cube_int_update_pos(pos2, radius, (*i + xlate), p_last2, 1, xy_only, cnorm_ptr);
		}
		// Note: driveways are handled elsewhere in the control flow
		if (!xy_only) { // don't need to check details and roof in xy_only mode because they're contained in the XY footprint of the parts
			for (auto i = details.begin(); i != details.end(); ++i) {
				if (sphere_cube_int_update_pos(pos2, radius, (*i + xlate), p_last2, 1, xy_only, cnorm_ptr)) {had_coll = 1;} // cube, flag as colliding
			}
			for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
				point const pos_xlate(pos2 - xlate);

				if (check_interior && had_coll && pos2.z - xlate.z > part_z2) { // player standing on top of a building with a sloped roof or roof access cover
					if (point_in_polygon_2d(pos_xlate.x, pos_xlate.y, i->pts, i->npts, 0, 1)) {
						vector3d const normal(i->get_norm());
						if (normal.z == 0.0) continue; // skip vertical sides as the player can't stand on them
						float const rdist(dot_product_ptv(normal, pos_xlate, i->pts[0]));

						if (draw_building_interiors && i->type == tquad_with_ix_t::TYPE_ROOF_ACC) { // don't allow walking on roof access tquads
							if (rdist < -radius*normal.z) continue; // player is below this tquad
							else {pos2.x = p_last2.x; pos2.y = p_last2.y; break;} // block the player from walking here (can only walk through raised opening)
						}
						pos2.z += (radius - rdist)/normal.z; // determine the distance we need to move vertically to achieve this diag separation
					}
				}
				else { // normal case for bouncing object, etc.
					vector3d const normal(i->get_norm());
					float const rdist(dot_product_ptv(normal, pos_xlate, i->pts[0]));

					if (fabs(rdist) < radius && sphere_poly_intersect(i->pts, i->npts, pos_xlate, normal, rdist, radius)) {
						pos2 += normal*(radius - rdist); // update current pos
						if (cnorm_ptr) {*cnorm_ptr = ((normal.z < 0.0) ? -1.0 : 1.0)*normal;} // make sure normal points up
						had_coll = 1; // flag as colliding
						break; // only use first colliding tquad
					}
				}
			} // for i
		}
	} // end !is_interior case
	if (!had_coll) return 0; // Note: no collisions with windows or doors, since they're colinear with walls

	if (is_rotated()) {
		do_xy_rotate(center, pos2); // rotate back around center
		if (cnorm_ptr) {do_xy_rotate_normal(*cnorm_ptr);} // rotate normal back
	}
	pos = pos2;
	return 1;
}

float room_object_t::get_radius() const {
	if (shape == SHAPE_CYLIN ) {return 0.25f*(dx() + dy());} // vertical cylinder: return average of x/y diameter
	if (shape == SHAPE_SPHERE) {return 0.5*dx();} // sphere, should be the same dx()/dy()/dz() value (but can't assert due to FP precision errors)
	assert(0); // cubes don't have a radius
	return 0.0; // never gets here
}
cylinder_3dw room_object_t::get_cylinder() const {
	float const radius(get_radius());
	point const center(get_cube_center());
	return cylinder_3dw(point(center.x, center.y, z1()), point(center.x, center.y, z2()), radius, radius);
}

// Note: returns bit vector for each cube that collides; supports up to 32 cubes
unsigned check_cubes_collision(cube_t const *const cubes, unsigned num_cubes, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	unsigned coll_ret(0);

	for (unsigned n = 0; n < num_cubes; ++n) {
		if (!cubes[n].is_all_zeros() && sphere_cube_int_update_pos(pos, radius, cubes[n], p_last, 1, 0, cnorm)) {coll_ret |= (1<<n);} // skip_z=0
	}
	return coll_ret;
}
unsigned check_closet_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[5];
	get_closet_cubes(c, cubes, 1); // get cubes for walls and door; required to handle collision with closet interior; for_collision=1
	// skip collision check of open doors for large closets since this case is more complex
	return check_cubes_collision(cubes, ((c.is_open() && !c.is_small_closet()) ? 4U : 5U), pos, p_last, radius, cnorm);
}
unsigned check_bed_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
	get_bed_cubes(c, cubes);
	unsigned num_to_check(5); // skip legs_bcube
	if (c.flags & RO_FLAG_TAKEN1) {--num_to_check;} // skip pillows
	if (c.flags & RO_FLAG_TAKEN3) {--num_to_check;} // skip mattress
	unsigned coll_ret(check_cubes_collision(cubes, num_to_check, pos, p_last, radius, cnorm));
	get_tc_leg_cubes(cubes[5], 0.04, cubes); // head_width=0.04
	coll_ret |= (check_cubes_collision(cubes, 4, pos, p_last, radius, cnorm) << 5); // check legs
	return coll_ret;
}
// actually applies to tables, desks, dressers, and nightstands
unsigned check_table_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[5];
	get_table_cubes(c, cubes); // body and legs
	return check_cubes_collision(cubes, 5, pos, p_last, radius, cnorm);
}
unsigned check_chair_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	cube_t cubes[3]; // seat, back, legs_bcube
	get_chair_cubes(c, cubes);
	return check_cubes_collision(cubes, 3, pos, p_last, radius, cnorm);
}
// Note: these next two are intended to be called when maybe_inside_room_object() returns true
bool check_stall_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	float const width(c.get_sz_dim(!c.dim));
	cube_t sides[3] = {c, c, c};
	sides[0].d[!c.dim][1] -= 0.95*width;
	sides[1].d[!c.dim][0] += 0.95*width;
	if (!c.is_open()) {sides[2].d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*0.975*c.get_sz_dim(c.dim);} // check collision with closed door
	bool had_coll(0);
	for (unsigned d = 0; d < (c.is_open() ? 2U : 3U); ++d) {had_coll |= sphere_cube_int_update_pos(pos, radius, sides[d], p_last, 1, 0, cnorm);}
	return had_coll;
}
bool check_shower_collision(room_object_t const &c, point &pos, point const &p_last, float radius, vector3d *cnorm) {
	bool const door_dim(c.dx() < c.dy()), door_dir(door_dim ? c.dim : c.dir), side_dir(door_dim ? c.dir : c.dim); // {c.dim, c.dir} => {dir_x, dir_y}
	cube_t sides[2] = {c, c};
	sides[0].d[!door_dim][!side_dir] -= (side_dir ? -1.0 : 1.0)*0.95*c.get_sz_dim(!door_dim); // shrink to just the outer glass wall of the shower
	if (!c.is_open()) {sides[1].d[door_dim][!door_dir] += (door_dir ? 1.0 : -1.0)*0.95*c.get_sz_dim(door_dim);} // check collision with closed door
	bool had_coll(0);
	for (unsigned d = 0; d < (c.is_open() ? 1U : 2U); ++d) {had_coll |= sphere_cube_int_update_pos(pos, radius, sides[d], p_last, 1, 0, cnorm);}
	return had_coll;
}
bool maybe_inside_room_object(room_object_t const &obj, point const &pos, float radius) {
	return ((obj.is_open() && sphere_cube_intersect(pos, radius, obj)) || obj.contains_pt(pos));
}

cube_t get_closet_bcube_including_door(room_object_t const &c) {
	if (c.type != TYPE_CLOSET || !c.is_open() || !c.is_small_closet()) return c;
	cube_t bcube(c); // only applies to small closets with open doors
	float const width(c.get_sz_dim(!c.dim)), wall_width(0.5*(width - 0.5*c.dz())); // see get_closet_cubes()
	bcube.d[c.dim][c.dir] += (c.dir ? 1.0f : -1.0f)*(width - 2.0f*wall_width); // extend outward
	return bcube;
}

// Note: used for the player; pos and p_last are already in rotated coordinate space
// default player is actually too large to fit through doors and too tall to fit between the floor and celing,
// so player size/height must be reduced in the config file
bool building_t::check_sphere_coll_interior(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, float radius, bool xy_only, vector3d *cnorm) const {
	pos.z = bcube.z1(); // start at building z1 rather than the terrain height in case we're at the foot of a steep hill
	assert(interior);
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness());
	float const xy_radius(radius*global_building_params.player_coll_radius_scale); // XY radius can be smaller to allow player to fit between furniture
	bool had_coll(0), on_stairs(0);
	float obj_z(max(pos.z, p_last.z)); // use p_last to get orig zval
	
	if (!xy_only && 2.2f*radius < (floor_spacing - floor_thickness)) { // diameter is smaller than space between floor and ceiling
		// check Z collision with floors; no need to check ceilings; this will set pos.z correctly so that we can set skip_z=0 in later tests
		float const floor_test_zval(obj_z + floor_thickness); // move up by floor thickness to better handle steep stairs

		for (auto i = interior->floors.begin(); i != interior->floors.end(); ++i) {
			if (!i->contains_pt_xy(pos)) continue; // sphere not in this part/cube
			float const z1(i->z1());
			if (floor_test_zval < z1 || floor_test_zval > z1 + floor_spacing) continue; // this is not the floor the sphere is on
			if (pos.z < z1 + radius) {pos.z = z1 + radius; obj_z = max(pos.z, p_last.z); had_coll = 1;} // move up
			break; // only change zval once
		}
	}
	// Note: this check must be after pos.z is set from interior->floors
	// pass in radius as wall_test_extra_z as a hack to allow player to step over a wall that's below the stairs connecting stacked parts
	had_coll |= interior->check_sphere_coll_walls_elevators_doors(*this, pos, p_last, xy_radius, radius, 0, cnorm); // check_open_doors=0 (to avoid getting the player stuck)

	if (interior->room_geom) { // collision with room geometry
		vector<room_object_t> const &objs(interior->room_geom->objs);

		for (auto c = interior->room_geom->get_stairs_start(); c != objs.end(); ++c) { // check for and handle stairs first
			if (c->no_coll() || c->type != TYPE_STAIR) continue;
			if (!c->contains_pt_xy(pos))  continue; // sphere not on this stair
			if (obj_z < c->z1())          continue; // below the stair
			if (pos.z - radius > c->z1()) continue; // above the stair
			pos.z = c->z1() + radius; // stand on the stair - this can happen for multiple stairs
			obj_z = max(pos.z, p_last.z);
			bool const is_u(c->shape == SHAPE_STAIRS_U);
			if (!is_u || c->dir == 1) {max_eq(pos[!c->dim], (c->d[!c->dim][0] + xy_radius));} // force the sphere onto the stairs
			if (!is_u || c->dir == 0) {min_eq(pos[!c->dim], (c->d[!c->dim][1] - xy_radius));}
			had_coll = on_stairs = 1;
		} // for c
		for (auto c = objs.begin(); c != objs.end(); ++c) { // check for other objects to collide with (including stairs)
			if (c->no_coll() || !bldg_obj_types[c->type].player_coll) continue;

			if (c->type == TYPE_ELEVATOR) { // special handling for elevators
				if (!c->contains_pt_xy(pos)) continue;
				if      (obj_z >= c->z2()) {max_eq(pos.z, (c->z2() + radius)); had_coll = player_in_elevator = 1;} // standing on the roof of the elevator
				else if (obj_z >= c->z1()) {max_eq(pos.z, (c->z1() + radius)); had_coll = player_in_elevator = 1;} // inside the elevator
				continue;
			}
			if ((c->type == TYPE_STAIR || on_stairs) && (obj_z + radius) > c->z2()) continue; // above the stair - allow it to be walked on
			cube_t c_extended(*c);
			if      (c->type == TYPE_STAIR ) {c_extended.z1() -= camera_zh;} // handle the player's head (for stairs), or is pos already at this height?
			else if (c->type == TYPE_CLOSET) {c_extended = get_closet_bcube_including_door(*c);}
			if (!sphere_cube_intersect(pos, xy_radius, c_extended)) continue; // optimization

			if (c->type == TYPE_RAILING) { // only collide with railing at top of stairs, not when walking under stairs
				cylinder_3dw const railing(get_railing_cylinder(*c));
				float const t((pos[c->dim] - railing.p1[c->dim])/(railing.p2[c->dim] - railing.p1[c->dim]));
				float const railing_zval(railing.p1.z + CLIP_TO_01(t)*(railing.p2.z - railing.p1.z));
				if ((railing_zval - get_railing_height(*c)) > float(pos.z + camera_zh) || railing_zval < (pos.z - radius)) continue; // no Z collision
			}
			if (c->shape == SHAPE_CYLIN) { // vertical cylinder
				cylinder_3dw cylin(c->get_cylinder());
				cylin.p2.z += radius; // extend upward by radius
				had_coll |= sphere_vert_cylin_intersect_with_ends(pos, xy_radius, cylin, cnorm);
			}
			else if (c->shape == SHAPE_SPHERE) { // sphere
				point const center(c->get_cube_center());
				if (!dist_less_than(pos, center, (c->get_radius() + radius))) continue; // xy_radius?
				if (cnorm) {*cnorm = (pos - center).get_norm();}
				had_coll = 1;
			}
			else if (c->type == TYPE_CLOSET) { // special case to handle closet interiors
				had_coll |= (bool)check_closet_collision(*c, pos, p_last, xy_radius, cnorm);
				
				if (c->contains_pt(pos)) {
					player_in_closet |= RO_FLAG_IN_CLOSET;
					if (interior->room_geom->closet_light_is_on(*c)) {player_in_closet |= RO_FLAG_LIT;}
					if (c->is_open()) {player_in_closet |= RO_FLAG_OPEN;} else {player_is_hiding = 1;} // player is hiding if the closet door is closed
				}
			}
			else if (c->type == TYPE_STALL && maybe_inside_room_object(*c, pos, xy_radius)) {
				// stall is open and intersecting player, or player is inside stall; perform collision test with sides only
				had_coll |= check_stall_collision(*c, pos, p_last, xy_radius, cnorm);
			}
			else if (c->type == TYPE_SHOWER && maybe_inside_room_object(*c, pos, xy_radius)) {
				// shower is open and intersecting player, or player is inside shower; perform collision test with side only
				had_coll |= check_shower_collision(*c, pos, p_last, xy_radius, cnorm);
			}
			else if (sphere_cube_int_update_pos(pos, xy_radius, c_extended, p_last, 1, 0, cnorm)) { // assume it's a cube; skip_z=0
				if (c->type == TYPE_TOILET || c->type == TYPE_URINAL) {player_near_toilet = 1;}
				had_coll = 1;
			}
			if ((c->type == TYPE_STALL || c->type == TYPE_SHOWER) && !c->is_open() && c->contains_pt(pos)) {player_is_hiding = 1;} // player is hiding in the stall/shower
		} // for c
	}
	for (auto i = ped_bcubes.begin(); i != ped_bcubes.end(); ++i) {
		float const ped_radius(0.5*max(i->dx(), i->dy())); // determine radius from bcube X/Y
		point const center(i->get_cube_center());
		float const dist(p2p_dist(pos, center)), r_sum(xy_radius + 0.5*ped_radius); // ped_radius is a bit too large
		if (dist >= r_sum) continue; // no intersection
		vector3d const normal(vector3d(pos.x-center.x, pos.y-center.y, 0.0).get_norm()); // XY direction
		if (cnorm) {*cnorm = normal;}
		pos += normal*(r_sum - dist);
		had_coll = 1;
	} // for i
	handle_vert_cylin_tape_collision(pos, p_last, pos.z-radius, pos.z+camera_zh, xy_radius, 1); // is_player=1
	return had_coll; // will generally always be true due to floors
}

// Note: called on basketballs and soccer balls
bool building_interior_t::check_sphere_coll(building_t const &building, point &pos, point const &p_last, float radius,
	vector<room_object_t>::const_iterator self, vector3d &cnorm, float &hardness, int &obj_ix) const
{
	bool had_coll(check_sphere_coll_walls_elevators_doors(building, pos, p_last, radius, 0.0, 1, &cnorm)); // check_open_doors=1
	if (had_coll) {hardness = 1.0;}
	cube_t scube; scube.set_from_sphere(pos, radius);

	if (pos.z < p_last.z) { // check floor collision if falling
		max_eq(scube.z2(), (p_last.z - radius)); // extend up to touch the bottom of the last position to prevent it from going through the floor in one frame

		for (auto i = floors.begin(); i != floors.end(); ++i) {
			if (!i->intersects(scube)) continue; // overlap
			pos.z = i->z2() + radius; // move to just touch the top of the floor
			cnorm = plus_z; // collision with top surface of floor
			had_coll = 1; hardness = 1.0;
		}
	}
	else { // check ceiling collision if rising
		min_eq(scube.z1(), (p_last.z + radius)); // extend down

		for (auto i = ceilings.begin(); i != ceilings.end(); ++i) {
			if (!i->intersects(scube)) continue; // overlap
			pos.z = i->z1() - radius; // move to just touch the top of the ceiling
			cnorm = -plus_z; // collision with top surface of ceiling
			had_coll = 1; hardness = 1.0;
		}
	}
	if (!room_geom) {return had_coll;} // no room geometry

	// Note: no collision check with expanded_objs
	for (auto c = room_geom->objs.begin(); c != room_geom->objs.end(); ++c) { // check for other objects to collide with
		// ignore blockers and railings, but allow more than c->no_coll()
		if (c == self || c->type == TYPE_BLOCKER || c->type == TYPE_RAILING || c->type == TYPE_PAPER || c->type == TYPE_PEN || c->type == TYPE_PENCIL ||
			c->type == TYPE_BOTTLE || c->type == TYPE_FLOORING || c->type == TYPE_SIGN || c->type == TYPE_WBOARD || c->type == TYPE_WALL_TRIM ||
			c->type == TYPE_DRAIN || c->type == TYPE_CRACK || c->type == TYPE_SWITCH) continue;
		if (!sphere_cube_intersect(pos, radius, ((c->type == TYPE_CLOSET) ? get_closet_bcube_including_door(*c) : *c))) continue; // no intersection (optimization)
		unsigned coll_ret(0);
		// add special handling for things like elevators and cubicles? right now these are only in office buildings, where there are no dynamic objects

		if (c->shape == SHAPE_CYLIN) { // vertical cylinder (including table)
			cylinder_3dw const cylin(c->get_cylinder());

			if (c->type == TYPE_TABLE) {
				cylinder_3dw top(cylin), base(cylin);
				top.p1.z  = base.p2.z = c->z2() - 0.12*c->dz(); // top shifted down by 0.12
				base.r1  *= 0.4; base.r2 *= 0.4; // vertical support has radius 0.08, legs have radius of 0.6, so use something in between
				coll_ret |= unsigned(sphere_vert_cylin_intersect_with_ends(pos, radius, top, &cnorm) || sphere_vert_cylin_intersect_with_ends(pos, radius, base, &cnorm));
			}
			else {coll_ret |= (unsigned)sphere_vert_cylin_intersect_with_ends(pos, radius, cylin, &cnorm);}
		}
		else if (c->shape == SHAPE_SPHERE) { // sphere
			point const center(c->get_cube_center());
			if (!dist_less_than(pos, center, (c->get_radius() + radius))) continue;
			cnorm = (pos - center).get_norm();
			coll_ret |= 1;
		}
		else { // assume it's a cube
			// some object types are special because they're common collision objects and they're not filled cubes
			if      (c->type == TYPE_CLOSET) {coll_ret |= check_closet_collision(*c, pos, p_last, radius, &cnorm);} // special case to handle closet interiors
			else if (c->type == TYPE_BED  )  {coll_ret |= check_bed_collision   (*c, pos, p_last, radius, &cnorm);}
			else if (c->type == TYPE_TABLE)  {coll_ret |= check_table_collision (*c, pos, p_last, radius, &cnorm);}
			else if (c->type == TYPE_DESK )  {coll_ret |= check_table_collision (*c, pos, p_last, radius, &cnorm);}
			else if (c->type == TYPE_CHAIR)  {coll_ret |= check_chair_collision (*c, pos, p_last, radius, &cnorm);}
			else if (c->type == TYPE_STALL  && maybe_inside_room_object(*c, pos, radius)) {coll_ret |= (unsigned)check_stall_collision (*c, pos, p_last, radius, &cnorm);}
			else if (c->type == TYPE_SHOWER && maybe_inside_room_object(*c, pos, radius)) {coll_ret |= (unsigned)check_shower_collision(*c, pos, p_last, radius, &cnorm);}
			else {coll_ret |= (unsigned)sphere_cube_int_update_pos(pos, radius, *c, p_last, 1, 0, &cnorm);} // skip_z=0
		}
		if (coll_ret) { // collision with this object - set hardness
			if      (c->type == TYPE_COUCH ) {hardness = 0.6;} // couches are soft
			else if (c->type == TYPE_RUG   ) {hardness = 0.8;} // rug is somewhat soft
			else if (c->type == TYPE_BLINDS) {hardness = 0.6;} // blinds are soft
			else if (c->type == TYPE_BED && (coll_ret & 24)) {hardness = 0.5;} // pillow/mattress collision is very soft
			else {hardness = 1.0;}
			obj_ix   = (c - room_geom->objs.begin()); // may be overwritten, will be the last collided object
			had_coll = 1;
		}
	} // for c
	return had_coll;
}

// Note: should be valid for players and other spherical objects
bool building_interior_t::check_sphere_coll_walls_elevators_doors(building_t const &building, point &pos, point const &p_last, float radius,
	float wall_test_extra_z, bool check_open_doors, vector3d *cnorm) const
{
	float const obj_z(max(pos.z, p_last.z)), wall_test_z(obj_z + wall_test_extra_z); // use p_last to get orig zval
	bool had_coll(0);

	// Note: pos.z may be too small here and we should really use obj_z, so skip_z must be set to 1 in cube tests and obj_z tested explicitly instead
	for (unsigned d = 0; d < 2; ++d) { // check XY collision with walls
		for (auto i = walls[d].begin(); i != walls[d].end(); ++i) {
			if (wall_test_z < i->z1() || wall_test_z > i->z2()) continue; // wrong part/floor
			had_coll |= sphere_cube_int_update_pos(pos, radius, *i, p_last, 1, 1, cnorm); // skip_z=1 (handled by zval test above)
		}
	}
	for (auto e = elevators.begin(); e != elevators.end(); ++e) {
		if (obj_z < e->z1() || obj_z > e->z2()) continue; // wrong part/floor

		if (room_geom && (e->open_amt > 0.75 || e->contains_pt(point(pos.x, pos.y, obj_z)))) { // elevator is mostly open, can enter || already inside elevator
			assert(e->car_obj_id < room_geom->objs.size());
			room_object_t &obj(room_geom->objs[e->car_obj_id]); // elevator car for this elevator

			if (obj_z > obj.z1() && obj_z < obj.z2()) { // same floor as elevator car - can enter it; otherwise can't enter elevator shaft
				cube_t cubes[5];
				unsigned const num_cubes(e->get_coll_cubes(cubes));
				for (unsigned n = 0; n < num_cubes; ++n) {had_coll |= sphere_cube_int_update_pos(pos, radius, cubes[n], p_last, 1, 1, cnorm);} // skip_z=1
				continue; // done with elevator
			}
		}
		had_coll |= sphere_cube_int_update_pos(pos, radius, *e, p_last, 1, 1, cnorm); // handle as a single blocking cube; skip_z=1
	} // for e
	for (auto i = doors.begin(); i != doors.end(); ++i) {
		if (i->open) {
			if (!check_open_doors) continue; // doors tend to block the player and other objects, don't collide with them unless they're closed
			cube_t door_bounds(*i);
			door_bounds.expand_by_xy(i->get_width());
			if (!sphere_cube_intersect(pos, radius, door_bounds)) continue; // check intersection with rough/conservative door bounds (optimization)
			tquad_with_ix_t const door(building.set_interior_door_from_cube(*i));
			vector3d normal(door.get_norm());
			if (dot_product_ptv(normal, pos, door.pts[0]) < 0.0) {normal.negate();} // use correct normal sign
			float rdist, thick;
			
			if (sphere_ext_poly_int_base(door.pts[0], normal, pos, radius, i->get_thickness(), thick, rdist)) {
				if (sphere_poly_intersect(door.pts, door.npts, pos, normal, rdist, thick)) {
					pos += normal*(thick - rdist);
					if (cnorm) {*cnorm = normal;}
					had_coll = 1;
				}
			}
			continue;
		}
		if (obj_z < i->z1() || obj_z > i->z2()) continue; // wrong part/floor
		had_coll |= sphere_cube_int_update_pos(pos, radius, *i, p_last, 1, 0, cnorm); // skip_z=0
	} // for i
	return had_coll;
}

bool get_line_clip_update_t(point const &p1, point const &p2, cube_t const &c, float &t) {
	float tmin(0.0), tmax(1.0);
	if (get_line_clip(p1, p2, c.d, tmin, tmax) && tmin < t) {t = tmin; return 1;}
	return 0;
}
bool get_line_clip_update_t(point const &p1, point const &p2, vect_cube_t const &cubes, float &t) {
	bool had_coll(0);
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {had_coll |= get_line_clip_update_t(p1, p2, *i, t);}
	return had_coll;
}
bool building_interior_t::line_coll(building_t const &building, point const &p1, point const &p2, point &p_int) const {
	bool had_coll(0);
	float t(1.0), tmin(1.0);
	had_coll |= get_line_clip_update_t(p1, p2, floors,   t);
	had_coll |= get_line_clip_update_t(p1, p2, ceilings, t);
	for (unsigned d = 0; d < 2; ++d) {had_coll |= get_line_clip_update_t(p1, p2, walls[d], t);}

	for (auto e = elevators.begin(); e != elevators.end(); ++e) {
		cube_t cubes[5];
		unsigned const num_cubes(e->get_coll_cubes(cubes));
		for (unsigned n = 0; n < num_cubes; ++n) {had_coll |= get_line_clip_update_t(p1, p2, cubes[n], t);}
	}
	for (auto i = doors.begin(); i != doors.end(); ++i) {
		if (i->open) {
			cube_t door_bounds(*i);
			door_bounds.expand_by_xy(i->get_width());
			if (!check_line_clip(p1, p2, door_bounds.d)) continue; // check intersection with rough/conservative door bounds (optimization)
			tquad_with_ix_t const door(building.set_interior_door_from_cube(*i));
			vector3d const normal(door.get_norm()), v1(p2 - p1);
			point pts[2][4];
			gen_poly_planes(door.pts, door.npts, normal, i->get_thickness(), pts);
			float const dp(dot_product(v1, normal));
			bool const test_side(dp > 0.0);
			
			if (thick_poly_intersect(v1, p1, normal, pts, test_side, door.npts)) {
				tmin = dot_product_ptv(normal, pts[!test_side][0], p1)/dp; // calculate t, since thick_poly_intersect() doesn't return it; use the far edge of the door
				if (tmin > 0.0 && tmin < t) {t = tmin; had_coll = 1;}
			}
		}
		else {had_coll |= get_line_clip_update_t(p1, p2, *i, t);} // closed
	}
	if (room_geom) { // check room geometry
		for (auto c = room_geom->objs.begin(); c != room_geom->objs.end(); ++c) { // check for other objects to collide with (including stairs)
			if ((c->no_coll() && c->type != TYPE_WALL_TRIM) || c->type == TYPE_BLOCKER || c->type == TYPE_ELEVATOR) continue; // keep wall trim but skip blockers and elevators

			if (c->type == TYPE_CLOSET) { // special case to handle closet interiors
				cube_t cubes[5];
				get_closet_cubes(*c, cubes, 1); // get cubes for walls and door; for_collision=1
				unsigned const n_end((c->is_open() && !c->is_small_closet()) ? 4U : 5U); // skip open doors for large closets since this case is more complex
				for (unsigned n = 0; n < n_end; ++n) {had_coll |= get_line_clip_update_t(p1, p2, cubes[n], t);}
			}
			else if (c->shape == SHAPE_CYLIN) { // vertical cylinder
				if (line_intersect_cylinder_with_t(p1, p2, c->get_cylinder(), 1, tmin) && tmin < t) {t = tmin; had_coll = 1;}
			}
			else if (c->shape == SHAPE_SPHERE) { // sphere
				float const radius(c->get_radius());
				if (sphere_test_comp(p1, c->get_cube_center(), (p1 - p2), radius*radius, tmin) && tmin < t) {t = tmin; had_coll = 1;}
			}
			//else if (c->type == TYPE_STALL || c->type == TYPE_SHOWER) {}
			else {had_coll |= get_line_clip_update_t(p1, p2, *c, t);} // assume it's a cube
		} // for c
	}
	if (had_coll) {p_int = p1 + t*(p2 - p1);}
	return had_coll;
}

void update_closest_pt(cube_t const &cube, point const &pos, point &closest, float pad_dist, float &dmin_sq) {
	cube_t cube_pad(cube);
	cube_pad.expand_by(pad_dist);
	point const cand(cube_pad.closest_pt(pos));
	float const dsq(p2p_dist_sq(pos, cand));
	if (dmin_sq < 0.0 || dsq < dmin_sq) {closest = cand; dmin_sq = dsq;}
}
void update_closest_pt(vect_cube_t const &cubes, point const &pos, point &closest, float pad_dist, float &dmin_sq) {
	for (auto c = cubes.begin(); c != cubes.end(); ++c) {update_closest_pt(*c, pos, closest, pad_dist, dmin_sq);}
}
point building_interior_t::find_closest_pt_on_obj_to_pos(building_t const &building, point const &pos, float pad_dist, bool no_ceil_floor) const {
	float dmin_sq(-1.0); // start at an invalid value
	point closest(pos); // start at pt - will keep this value if there are no objects

	if (!no_ceil_floor) {
		update_closest_pt(floors,   pos, closest, pad_dist, dmin_sq);
		update_closest_pt(ceilings, pos, closest, pad_dist, dmin_sq);
	}
	for (unsigned d = 0; d < 2; ++d) {update_closest_pt(walls[d], pos, closest, pad_dist, dmin_sq);}
	for (auto e = elevators.begin(); e != elevators.end(); ++e) {update_closest_pt(*e, pos, closest, pad_dist, dmin_sq);} // ignores open elevator doors

	for (auto i = doors.begin(); i != doors.end(); ++i) {
		if (i->open) {} // handle open doors? - closest point on extruded polygon
		else {update_closest_pt(*i, pos, closest, pad_dist, dmin_sq);} // closed
	}
	if (room_geom) { // check room geometry
		for (auto c = room_geom->objs.begin(); c != room_geom->objs.end(); ++c) { // check for other objects to collide with (including stairs)
			if ((c->no_coll() && c->type != TYPE_WALL_TRIM) || c->type == TYPE_BLOCKER || c->type == TYPE_ELEVATOR) continue; // keep wall trim but skip blockers and elevators

			if (c->type == TYPE_CLOSET) { // special case to handle closet interiors
				cube_t cubes[5];
				get_closet_cubes(*c, cubes, 1); // get cubes for walls and door; for_collision=1
				unsigned const n_end((c->is_open() && !c->is_small_closet()) ? 4U : 5U); // skip open doors for large closets since this case is more complex
				for (unsigned n = 0; n < n_end; ++n) {update_closest_pt(cubes[n], pos, closest, pad_dist, dmin_sq);}
			}
			if (c->shape == SHAPE_CYLIN) { // vertical cylinder
				point const center(c->get_cube_center());
				float const radius(c->get_radius() + pad_dist), dsq_xy(max(0.0f, (p2p_dist_xy_sq(pos, center) - radius*radius)));
				float const zval(max(c->z1()-pad_dist, min(c->z2()+pad_dist, pos.z))), dz(pos.z - zval), dsq(dsq_xy + dz*dz);
				if (dmin_sq < 0.0 || dsq < dmin_sq) {closest = point(center.x, center.y, zval) + radius*(point(pos.x, pos.y, 0.0) - point(center.x, center.y, 0.0)).get_norm(); dmin_sq = dsq;}
			}
			else if (c->shape == SHAPE_SPHERE) { // sphere
				point const center(c->get_cube_center());
				float const radius(c->get_radius() + pad_dist), dsq(max(0.0f, (p2p_dist_sq(pos, center) - radius*radius)));
				if (dmin_sq < 0.0 || dsq < dmin_sq) {closest = center + radius*(pos - center).get_norm(); dmin_sq = dsq;}
			}
			else {update_closest_pt(*c, pos, closest, pad_dist, dmin_sq);} // assume it's a cube
		} // for c
	}
	// what about exterior building walls?
	return closest;
}

// Note: p1/p2 are in building space; return value: 0=none, 1=side, 2=roof, 3=details
unsigned building_t::check_line_coll(point const &p1, point const &p2, float &t, vector<point> &points, bool occlusion_only, bool ret_any_pt, bool no_coll_pt) const {
	point p1r(p1), p2r(p2); // copy before clipping
	if (!check_line_clip(p1r, p2r, bcube.d)) return 0; // no intersection

	if (is_rotated()) {
		point const center(bcube.get_cube_center());
		do_xy_rotate_inv(center, p1r); // inverse rotate - negate the sine term
		do_xy_rotate_inv(center, p2r);
	}
	float const pzmin(min(p1r.z, p2r.z)), pzmax(max(p1r.z, p2r.z));
	bool const vert(p1r.x == p2r.x && p1r.y == p2r.y);
	float tmin(0.0);
	unsigned coll(0); // 0=none, 1=side, 2=roof, 3=details

	for (auto i = parts.begin(); i != parts.end(); ++i) {
		if (pzmin > i->z2() || pzmax < i->z1()) continue; // no overlap in z
		bool hit(0);

		if (use_cylinder_coll()) { // vertical cylinder
			// Note: we know the line intersects the cylinder's bcube, and there's a good chance it intersects the cylinder, so we don't need any expensive early termination cases here
			point const cc(i->get_cube_center());
			vector3d const csz(i->get_size());

			if (vert) { // vertical line + vertical cylinder optimization + handling of ellipsoids
				if (!point_in_ellipse(p1r, cc, 0.5*csz.x, 0.5*csz.y)) continue; // no intersection (below test should return true as well)
				tmin = (i->z2() - p1r.z)/(p2r.z - p1r.z);
				if (tmin >= 0.0 && tmin < t) {t = tmin; hit = 1;}
			}
			else {
				float const radius(0.5*(occlusion_only ? min(csz.x, csz.y) : max(csz.x, csz.y))); // use conservative radius unless this is an occlusion query
				point const cp1(cc.x, cc.y, i->z1()), cp2(cc.x, cc.y, i->z2());
				if (!line_int_cylinder(p1r, p2r, cp1, cp2, radius, radius, 1, tmin) || tmin > t) continue; // conservative for non-occlusion rays

				if (!occlusion_only && csz.x != csz.y) { // ellipse
					vector3d const delta(p2r - p1r);
					float const rx_inv_sq(1.0/(0.25*csz.x*csz.x)), ry_inv_sq(1.0/(0.25*csz.y*csz.y));
					float t_step(0.1*max(csz.x, csz.y)/delta.mag());

					for (unsigned n = 0; n < 10; ++n) { // use an interative approach
						if (point_in_ellipse_risq((p1r + tmin*delta), cc, rx_inv_sq, ry_inv_sq)) {hit = 1; tmin -= t_step;} else {tmin += t_step;}
						if (hit) {t_step *= 0.5;} // converge on hit point
					}
					if (!hit) continue; // not actually a hit
				} // end ellipse case
				t = tmin; hit = 1;
			}
		}
		else if (num_sides != 4) {
			building_draw_utils::calc_poly_pts(*this, *i, points);
			float const tz((i->z2() - p1r.z)/(p2r.z - p1r.z)); // t value at zval = top of cube

			if (tz >= 0.0 && tz < t) {
				float const xval(p1r.x + tz*(p2r.x - p1r.x)), yval(p1r.y + tz*(p2r.y - p1r.y));
				if (point_in_polygon_2d(xval, yval, points.data(), points.size(), 0, 1)) {t = tz; hit = 1;} // XY plane test for vertical lines and top surface
			}
			if (!vert) { // test building sides
				point quad_pts[4]; // quads

				for (unsigned S = 0; S < num_sides; ++S) { // generate vertex data quads
					for (unsigned d = 0, ix = 0; d < 2; ++d) {
						point const &p(points[(S+d)%num_sides]);
						for (unsigned e = 0; e < 2; ++e) {quad_pts[ix++].assign(p.x, p.y, i->d[2][d^e]);}
					}
					if (line_poly_intersect(p1r, p2r, quad_pts, 4, get_poly_norm(quad_pts), tmin) && tmin < t) {t = tmin; hit = 1;} // Note: untested
				} // for S
			}
		}
		else {hit |= get_line_clip_update_t(p1r, p2r, *i, t);} // cube

		if (hit) {
			if (occlusion_only) return 1; // early exit
			if (vert) {coll = 2;} // roof
			else {
				float const zval(p1r.z + t*(p2.z - p1r.z));
				coll = ((fabs(zval - i->z2()) < 0.0001*i->dz()) ? 2 : 1); // test if clipped zval is close to the roof zval
			}
			if (ret_any_pt) return coll;
		}
	} // for i
	if (occlusion_only) return 0;

	for (auto i = details.begin(); i != details.end(); ++i) {
		if (get_line_clip_update_t(p1r, p2r, *i, t)) {coll = 3;} // details cube
	}
	if (!no_coll_pt || !vert) { // vert line already tested building cylins/cubes, and marked coll roof, no need to test again unless we need correct coll_pt t-val
		for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
			if (line_poly_intersect(p1r, p2r, i->pts, i->npts, i->get_norm(), tmin) && tmin < t) {t = tmin; coll = 2;} // roof quad
		}
	}
	if (!vert) { // don't need to check fences for vertical collisions since they're horizontal blockers
		for (auto i = fences.begin(); i != fences.end(); ++i) {
			if (get_line_clip_update_t(p1r, p2r, *i, t)) {coll = 3;} // counts as details cube
		}
	}
	return coll; // Note: no collisions with windows or doors, since they're colinear with walls; no collision with interior for now
}

// Note: if xy_radius == 0.0, this is a point test; otherwise, it's an approximate vertical cylinder test
bool building_t::check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points) const {

	if (xy_radius == 0.0 && !bcube.contains_pt(pos)) return 0; // no intersection (bcube does not need to be rotated)
	point pr(pos);
	maybe_inv_rotate_point(pr);

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
		if (pr.z > i->z2() || pr.z < i->z1()) continue; // no overlap in z

		if (use_cylinder_coll()) { // vertical cylinder
			point const cc(i->get_cube_center());
			vector3d const csz(i->get_size());
			float const dx(cc.x - pr.x), dy(cc.y - pr.y), rx(0.5*csz.x + xy_radius), ry(0.5*csz.y + xy_radius);
			if (dx*dx/(rx*rx) + dy*dy/(ry*ry) > 1.0f) continue; // no intersection (below test should return true as well)
			return 1;
		}
		else if (num_sides != 4) {
			building_draw_utils::calc_poly_pts(*this, *i, points);

			if (xy_radius > 0.0) { // cylinder case: expand polygon by xy_radius; assumes a convex polygon
				point const center(i->get_cube_center());

				for (auto p = points.begin(); p != points.end(); ++p) {
					vector3d dir(*p - center);
					dir.z = 0.0; // only want XY component
					*p += dir*(xy_radius/dir.mag());
				}
			}
			if (point_in_polygon_2d(pr.x, pr.y, &points.front(), points.size(), 0, 1)) return 1; // XY plane test for top surface
		}
		else { // cube
			if (xy_radius > 0.0) {
				cube_t cube(*i);
				cube.expand_by(xy_radius);
				if (cube.contains_pt(pr)) return 1;
			}
			else if (i->contains_pt(pr)) return 1;
		}
	} // for i
	return 0;
}

bool building_t::check_point_xy_in_part(point const &pos) const { // simpler/faster version of check_point_or_cylin_contained() with no z check
	if (!bcube.contains_pt_xy(pos)) return 0; // no intersection (bcube does not need to be rotated)
	point pr(pos);
	maybe_inv_rotate_point(pr);

	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
		if (i->contains_pt_xy(pr)) return 1;
	}
	return 0;
}

bool has_cube_line_coll(point const &p1, point const &p2, vect_cube_t const &cubes) {
	for (auto i = cubes.begin(); i != cubes.end(); ++i) {
		if (i->line_intersects(p1, p2)) return 1;
	}
	return 0;
}
bool building_t::check_for_wall_ceil_floor_int(point const &p1, point const &p2) const { // and interior doors
	if (!interior) return 0;
	for (unsigned d = 0; d < 2; ++d) {if (has_cube_line_coll(p1, p2, interior->walls[d])) return 1;}
	if (has_cube_line_coll(p1, p2, interior->ceilings) || has_cube_line_coll(p1, p2, interior->floors)) return 1; // or is only checking one good enough?
	return check_line_intersect_doors(p1, p2);
}

