// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#pragma warning(disable : 26812) // prefer enum class over enum

cube_t grass_exclude1, grass_exclude2;

extern bool draw_building_interiors;
extern int display_mode;
extern float grass_width;
extern building_params_t global_building_params;


void building_t::set_z_range(float z1, float z2) {
	bcube.z1() = z1; bcube.z2() = z2;
	adjust_part_zvals_for_floor_spacing(bcube);
	if (!parts.empty()) {parts[0].z1() = z1; parts[0].z2() = z2;}
}
building_mat_t const &building_t::get_material() const {return global_building_params.get_material(mat_ix);}

void building_t::gen_rotation(rand_gen_t &rgen) {

	float const max_rot_angle(get_material().max_rot_angle);
	if (max_rot_angle == 0.0) return;
	float const rot_angle(rgen.rand_uniform(0.0, max_rot_angle));
	rot_sin = sin(rot_angle);
	rot_cos = cos(rot_angle);
	parts.clear();
	parts.push_back(bcube); // this is the actual building base
	cube_t const &bc(parts.back());
	point const center(bc.get_cube_center());

	for (unsigned i = 0; i < 4; ++i) {
		point corner(bc.d[0][i&1], bc.d[1][i>>1], bc.d[2][i&1]);
		do_xy_rotate(rot_sin, rot_cos, center, corner);
		if (i == 0) {bcube.set_from_point(corner);} else {bcube.union_with_pt(corner);} // Note: detail cubes are excluded
	}
}

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
			do_xy_rotate(b.rot_sin, b.rot_cos, center1, pts[0]); // rotate into global space
		}
		cube_t c_exp(*p1);
		c_exp.expand_by_xy(expand_rel*p1->get_size() + vector3d(expand_abs, expand_abs, expand_abs));

		for (unsigned i = 0; i < 4; ++i) { // {00, 10, 01, 11}
			pts[i+1].assign(c_exp.d[0][i&1], c_exp.d[1][i>>1], 0.0); // XY only
			do_xy_rotate(b.rot_sin, b.rot_cos, center1, pts[i+1]); // rotate into global space
		}
		for (unsigned i = 0; i < 5; ++i) {do_xy_rotate(-rot_sin, rot_cos, center2, pts[i]);} // inverse rotate into local coord space - negate the sine term
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

bool building_t::test_coll_with_sides(point &pos, point const &p_last, float radius, cube_t const &part, vector<point> &points, vector3d *cnorm) const {

	building_draw_utils::calc_poly_pts(*this, part, points); // without the expand
	point quad_pts[4]; // quads
	bool updated(0);

	// FIXME: if the player is moving too quickly, the intersection with a side polygon may be missed,
	// which allows the player to travel through the building, but using a line intersection test from p_past2 to pos has other problems
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
		updated = 1;
	} // for S
	if (updated) return 1;

	if (max(pos.z, p_last.z) > part.z2() && point_in_polygon_2d(pos.x, pos.y, points.data(), num_sides, 0, 1)) { // test top plane (sphere on top of polygon?)
		pos.z = part.z2() + radius; // make sure it doesn't intersect the roof
		if (cnorm) {*cnorm = plus_z;}
		return 1;
	}
	return 0;
}

bool building_t::check_sphere_coll(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, vector3d const &xlate,
	float radius, bool xy_only, vector<point> &points, vector3d *cnorm_ptr, bool check_interior) const
{
	if (!is_valid()) return 0; // invalid building
	point p_int;
	vector3d cnorm; // unused
	unsigned cdir(0); // unused
	if (radius > 0.0 && !sphere_cube_intersect(pos, radius, (bcube + xlate), p_last, p_int, cnorm, cdir, 1, xy_only)) return 0;
	point pos2(pos), p_last2(p_last), center;
	bool had_coll(0), is_interior(0);
	float part_z2(bcube.z2());

	if (is_rotated()) {
		center = bcube.get_cube_center() + xlate;
		do_xy_rotate(-rot_sin, rot_cos, center, pos2); // inverse rotate - negate the sine term
		do_xy_rotate(-rot_sin, rot_cos, center, p_last2);
	}
	if (check_interior && draw_building_interiors && interior != nullptr) { // check for interior case first
		float const zval(max(pos2.z, p_last2.z) - xlate.z); // this is the real zval for use in collsion detection, in building space
		cube_t sc; sc.set_from_sphere((pos2 - xlate), radius); // sphere bounding cube

		if (zval > bcube.z1() && zval < (bcube.z1() + get_door_height())) { // on the ground floor
			for (auto d = doors.begin(); d != doors.end(); ++d) {
				if (d->get_bcube().intersects_xy(sc)) return 0; // check if we can use a door - disable collsion detection to allow the player to walk through
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
			for (auto i = interior->stairwells.begin(); i != interior->stairwells.end(); ++i) {
				if (!i->roof_access) continue;
				cube_t test_cube(*i);
				test_cube.expand_by_xy(-0.5f*radius);
				if (!test_cube.contains_pt_xy(pos2 - xlate)) continue; // pos not over stairs
				if (zval < i->z2() + 1.5f*radius) {is_interior = 1; break;} // Note: don't have to check zval > i->z2() because we know that !is_interior
			}
		}
	}
	if (is_interior) {
		point pos2_bs(pos2 - xlate);
		if (check_sphere_coll_interior(pos2_bs, (p_last - xlate), ped_bcubes, radius, xy_only, cnorm_ptr)) {pos2 = pos2_bs + xlate; had_coll = 1;}
	}
	else {
		for (auto i = parts.begin(); i != parts.end(); ++i) {
			if (xy_only && i->d[2][0] > bcube.d[2][0]) break; // only need to check first level in this mode
			if (!xy_only && ((pos2.z + radius < i->d[2][0] + xlate.z) || (pos2.z - radius > i->d[2][1] + xlate.z))) continue; // test z overlap
			if (radius == 0.0 && !(xy_only ? i->contains_pt_xy(pos2) : i->contains_pt(pos2))) continue; // no intersection; ignores p_last
			cube_t const part_bc(*i + xlate);
			bool part_coll(0);

			if (use_cylinder_coll()) {
				point const cc(part_bc.get_cube_center());
				float const crx(0.5*i->dx()), cry(0.5*i->dy()), r_sum(radius + max(crx, cry));
				if (!dist_xy_less_than(pos2, cc, r_sum)) continue; // no intersection

				if (fabs(crx - cry) < radius) { // close to a circle
					if (p_last2.z > part_bc.d[2][1] && dist_xy_less_than(pos2, cc, max(crx, cry))) {
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
			else if (!xy_only && i->contains_pt_xy_exp(pos2, radius) && (p_last2.z - 0.9f*radius) > (i->d[2][1] + xlate.z)) { // on top of building
				pos2.z = i->d[2][1] + xlate.z + radius;
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

						if (i->type == tquad_with_ix_t::TYPE_ROOF_ACC) { // don't allow walking on roof access tquads
							if (rdist < -radius*normal.z) continue; // player is below this tquad
							else {pos2.x = p_last.x; pos2.y = p_last.y; break;} // block the player from walking here (can only walk through raised opening)
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
		do_xy_rotate(rot_sin, rot_cos, center, pos2); // rotate back around center
		if (cnorm_ptr) {do_xy_rotate(rot_sin, rot_cos, all_zeros, *cnorm_ptr);} // rotate back (pure rotation)
	}
	pos = pos2;
	return had_coll;
}

// Note: pos and p_last are already in rotated coordinate space
// default player is actually too large to fit through doors and too tall to fit between the floor and celing, so player size/height must be reduced in the config file
bool building_t::check_sphere_coll_interior(point &pos, point const &p_last, vect_cube_t const &ped_bcubes, float radius, bool xy_only, vector3d *cnorm) const {
	assert(interior);
	float const floor_spacing(get_window_vspace());
	bool had_coll(0), on_stairs(0);
	float obj_z(max(pos.z, p_last.z)); // use p_last to get orig zval

	for (unsigned d = 0; d < 2; ++d) { // check XY collision with walls
		for (auto i = interior->walls[d].begin(); i != interior->walls[d].end(); ++i) {
			if (obj_z < i->z1() || obj_z > i->z2()) continue; // wrong part/floor
			had_coll |= sphere_cube_int_update_pos(pos, radius, *i, p_last, 1, 0, cnorm); // skip_z=0 (required for stacked parts that have diff walls)
		}
	}
	// for now, players aren't allowed in elevators
	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		if (e->open) {} // maybe they should be allowed to fall down elevator shafts if the elevator door is open?
		if (obj_z < e->z1() || obj_z > e->z2()) continue; // wrong part/floor
		had_coll |= sphere_cube_int_update_pos(pos, radius, *e, p_last, 1, 0, cnorm); // skip_z=0
	}
	/*for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) { // doors tend to block the player, don't collide with them
		if (obj_z < i->z1() || obj_z > i->z2()) continue; // wrong part/floor
		had_coll |= sphere_cube_int_update_pos(pos, radius, *i, p_last, 1, 0, cnorm); // skip_z=0
	}*/
	if (!xy_only && 2.2*radius < floor_spacing*(1.0 - FLOOR_THICK_VAL)) { // diameter is smaller than space between floor and ceiling
		// check Z collision with floors; no need to check ceilings
		obj_z = max(pos.z, p_last.z);

		for (auto i = interior->floors.begin(); i != interior->floors.end(); ++i) {
			if (!i->contains_pt_xy(pos)) continue; // sphere not in this part/cube
			float const z1(i->z1());
			if (obj_z < z1 || obj_z > z1 + floor_spacing) continue; // this is not the floor the sphere is on
			if (pos.z < z1 + radius) {pos.z = z1 + radius; had_coll = 1;} // move up
			break; // only change zval once
		}
	}
	if (interior->room_geom) { // collision with room cubes
		vector<room_object_t> const &objs(interior->room_geom->objs);
		obj_z = max(pos.z, p_last.z);

		for (auto c = (objs.begin() + interior->room_geom->stairs_start); c != objs.end(); ++c) { // check for and handle stairs first
			if (c->no_coll() || c->type != TYPE_STAIR) continue;
			if (!c->contains_pt_xy(pos))  continue; // sphere not on this stair
			if (obj_z < c->z1())          continue; // below the stair
			if (pos.z - radius > c->z1()) continue; // above the stair
			pos.z = c->z1() + radius; // stand on the stair - this can happen for multiple stairs
			obj_z = max(pos.z, p_last.z);
			bool const is_u(c->shape == SHAPE_STAIRS_U);
			if (!is_u || c->dir == 1) {max_eq(pos[!c->dim], (c->d[!c->dim][0] + radius));} // force the sphere onto the stairs
			if (!is_u || c->dir == 0) {min_eq(pos[!c->dim], (c->d[!c->dim][1] - radius));}
			had_coll = on_stairs = 1;
		} // for c
		for (auto c = objs.begin(); c != objs.end(); ++c) { // check for other objects to collide with
			if (c->no_coll()) continue;
			if ((c->type == TYPE_STAIR || on_stairs) && (obj_z + radius) > c->z2()) continue; // above the stair - allow it to be walked on
			had_coll |= sphere_cube_int_update_pos(pos, radius, *c, p_last, 1, 0, cnorm); // skip_z=0
		}
	}
	for (auto i = ped_bcubes.begin(); i != ped_bcubes.end(); ++i) {
		float const ped_radius(0.5*max(i->dx(), i->dy())); // determine radius from bcube X/Y
		point const center(i->get_cube_center());
		float const dist(p2p_dist(pos, center)), r_sum(radius + 0.5*ped_radius); // ped_radius is a bit too large
		if (dist >= r_sum) continue; // no intersection
		vector3d const normal(vector3d(pos.x-center.x, pos.y-center.y, 0.0).get_norm()); // XY direction
		if (cnorm) {*cnorm = normal;}
		pos += normal*(r_sum - dist);
		had_coll = 1;
	} // for i
	return had_coll; // will generally always be true due to floors
}

unsigned building_t::check_line_coll(point const &p1, point const &p2, vector3d const &xlate, float &t, vector<point> &points,
	bool occlusion_only, bool ret_any_pt, bool no_coll_pt) const
{
	if (!check_line_clip(p1-xlate, p2-xlate, bcube.d)) return 0; // no intersection
	point p1r(p1), p2r(p2);
	float tmin(0.0), tmax(1.0);
	unsigned coll(0); // 0=none, 1=side, 2=roof, 3=details

	if (is_rotated()) {
		point const center(bcube.get_cube_center() + xlate);
		do_xy_rotate(-rot_sin, rot_cos, center, p1r); // inverse rotate - negate the sine term
		do_xy_rotate(-rot_sin, rot_cos, center, p2r);
	}
	p1r -= xlate; p2r -= xlate;
	float const pzmin(min(p1r.z, p2r.z)), pzmax(max(p1r.z, p2r.z));
	bool const vert(p1r.x == p2r.x && p1r.y == p2r.y);

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
		else if (get_line_clip(p1r, p2r, i->d, tmin, tmax) && tmin < t) {t = tmin; hit = 1;} // cube

		if (hit) {
			if (occlusion_only) return 1; // early exit
			if (vert) {coll = 2;} // roof
			else {
				float const zval(p1.z + t*(p2.z - p1.z));
				coll = ((fabs(zval - i->d[2][1]) < 0.0001*i->dz()) ? 2 : 1); // test if clipped zval is close to the roof zval
			}
			if (ret_any_pt) return coll;
		}
	} // for i
	if (occlusion_only) return 0;

	for (auto i = details.begin(); i != details.end(); ++i) {
		if (get_line_clip(p1r, p2r, i->d, tmin, tmax) && tmin < t) {t = tmin; coll = 3;} // details cube
	}
	if (!no_coll_pt || !vert) { // vert line already tested building cylins/cubes, and marked coll roof, no need to test again unless we need correct coll_pt t-val
		for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {
			if (line_poly_intersect(p1r, p2r, i->pts, i->npts, i->get_norm(), tmin) && tmin < t) {t = tmin; coll = 2;} // roof quad
		}
	}
	return coll; // Note: no collisions with windows or doors, since they're colinear with walls; no collision with interior for now
}

// Note: if xy_radius == 0.0, this is a point test; otherwise, it's an approximate vertical cylinder test
bool building_t::check_point_or_cylin_contained(point const &pos, float xy_radius, vector<point> &points) const {

	if (xy_radius == 0.0 && !bcube.contains_pt(pos)) return 0; // no intersection
	point pr(pos);
	if (is_rotated()) {do_xy_rotate(-rot_sin, rot_cos, bcube.get_cube_center(), pr);} // inverse rotate - negate the sine term

	for (auto i = parts.begin(); i != parts.end(); ++i) {
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

void building_t::calc_bcube_from_parts() {
	assert(!parts.empty());
	bcube = parts[0];
	for (auto i = parts.begin()+1; i != parts.end(); ++i) {bcube.union_with_cube(*i);} // update bcube
}

void building_t::move_door_to_other_side_of_wall(tquad_with_ix_t &door, float dist_mult, bool invert_normal) const {
	cube_t const c(door.get_bcube());
	bool const dim(c.dy() < c.dx()), dir(door.get_norm()[dim] > 0.0); // closest cube side dir
	float door_shift(0.0);

	if (door.type == tquad_with_ix_t::TYPE_RDOOR) { // not on a wall, can't move it (should caller check this?)
		door_shift = 0.001*(dir ? 1.0 : -1.0)*get_material().get_floor_spacing();
	}
	else {
		door_shift = bcube.dz(); // start with a large value
		if (invert_normal) {swap(door.pts[0], door.pts[1]); swap(door.pts[2], door.pts[3]);} // swap vertex order to invert normal

		for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) { // find the part that this door was added to (inc garages and sheds)
			float const dist(door.pts[0][dim] - p->d[dim][dir]); // signed
			if (fabs(dist) < fabs(door_shift)) {door_shift = dist;}
		}
		assert(fabs(door_shift) < bcube.dz());
	}
	door_shift *= -(1.0 + dist_mult); // reflect on other side
	for (unsigned n = 0; n < door.npts; ++n) {door.pts[n][dim] += door_shift;} // move to opposite side of wall
}

void building_t::clip_door_to_interior(tquad_with_ix_t &door, bool clip_to_floor) const {

	cube_t clip_cube(door.get_bcube());
	float const dz(clip_cube.dz());
	float const border((door.type == tquad_with_ix_t::TYPE_BDOOR) ? 0.04 : 0.08);
	// clip off bottom for floor if clip_to_floor==1; somewhat arbitrary, should we use interior->floors.back().z2() instead?
	float const clip_bot(clip_to_floor ? 0.7*FLOOR_THICK_VAL*get_material().get_floor_spacing() : 0.04*dz);
	clip_cube.z1() += clip_bot;
	clip_cube.z2() -= 0.5*border*dz;
	bool const dim(clip_cube.dx() < clip_cube.dy()); // border dim
	float const xy_border(border*clip_cube.get_sz_dim(dim));
	clip_cube.d[dim][0] += xy_border; clip_cube.d[dim][1] -= xy_border; // shrink by border
	for (unsigned n = 0; n < door.npts; ++n) {clip_cube.clamp_pt(door.pts[n]);}
}

cube_t building_t::get_part_containing_pt(point const &pt) const {
	for (auto i = parts.begin(); i != (parts.end() - has_chimney); ++i) { // includes garage/shed
		if (i->contains_pt(pt)) {return *i;}
	}
	// we can get here in rare cases due to FP precision problems; find the closest cube to pt, which should be very close
	float dmin_sq(0.0);
	cube_t closest;

	for (auto i = parts.begin(); i != (parts.end() - has_chimney); ++i) {
		float const dist_sq(p2p_dist(pt, i->closest_pt(pt)));
		if (dmin_sq == 0.0 || dist_sq < dmin_sq) {dmin_sq = dist_sq; closest = *i;}
	}
	assert(!closest.is_all_zeros()); // must be found
	return closest;
}

void building_t::adjust_part_zvals_for_floor_spacing(cube_t &c) const {

	if (!EXACT_MULT_FLOOR_HEIGHT) return;
	float const floor_spacing(get_window_vspace()), dz(c.dz());
	assert(dz > 0.0 && floor_spacing > 0.0);
	float const num_floors(dz/floor_spacing);
	int const targ_num_floors(max(1, round_fp(num_floors)));
	c.z2() += floor_spacing*(targ_num_floors - num_floors); // ensure c.dz() is an exact multiple of num_floors
}

void split_cubes_recur(cube_t c, vect_cube_t &cubes, unsigned search_start, unsigned search_end) {

	for (unsigned i = search_start; i < search_end; ++i) {
		cube_t const &sc(cubes[i]);
		// while this generally holds, it can fail in rare cases, I assume due to floating-point error or some such thing, so this check has been disabled
		//assert(sc.z2() >= c.z2()); // assumes cubes are ordered descending by ztop
		if (!sc.intersects_no_adj(c)) continue;
		if (sc.contains_cube(c)) return; // contained, done (remove all of c)
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
	ao_bcz2 = bcube.z2(); // capture z2 before union with roof and detail geometry (which increases building height)
	if (is_house) {gen_house(base, rgen); return;}

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
	if (mat.min_level_height > 0.0) {num_levels = max(mat.min_levels, min(num_levels, unsigned(bcube.get_size().z/mat.min_level_height)));}
	num_levels = max(num_levels, 1U); // min_levels can be zero to apply more weight to 1 level buildings
	bool const do_split(num_levels < 4 && is_cube() && rgen.rand_probability(mat.split_prob)); // don't split buildings with 4 or more levels, or non-cubes

	if (num_levels == 1) { // single level
		if (do_split) { // generate L, T, or U shape
			split_in_xy(base, rgen);

			if (has_city_trees()) { // see if we can place a tree in the courtyard
				point const center(bcube.get_cube_center());
				cube_t place_area(center, center);
				place_area.expand_by_xy(0.05f*(bcube.dx() + bcube.dy()));
				if (!has_bcube_int(place_area, parts)) {tree_pos = place_area.get_cube_center(); tree_pos.z = bcube.z1();}
			}
			gen_details(rgen, 0);
		}
		else { // single part, entire cube/cylinder
			parts.push_back(base);
			if ((rgen.rand()&3) != 0) {maybe_add_special_roof(rgen);} // 75% chance
			
			if (0 && interior_enabled()) { // while this works, it doesn't seem to add much value, it only creates odd geometry and connecting stairs/elevators
				// two stacked parts of the same x/y dimensions but different interior floorplans
				float const dz(base.dz()), split_zval(rgen.rand_uniform(0.4, 0.6));
				parts.push_back(base);
				parts[0].z2() = parts[0].z1() + split_zval*dz; // split in Z: parts[0] is the bottom, parts[1] is the top
				adjust_part_zvals_for_floor_spacing(parts[0]);
				parts[1].z1() = parts[0].z2();
			}
			gen_details(rgen, 1);
		}
		end_add_parts();
		gen_interior(rgen, 0);
		gen_building_doors_if_needed(rgen);
		return; // for now the bounding cube
	}
	// generate building levels and splits
	float const height(base.dz()), dz(height/num_levels);
	float const abs_min_edge_move(0.5*mat.get_floor_spacing()); // same as door width
	bool const not_too_small(min(bcube.dx(), bcube.dy()) > 4.0*abs_min_edge_move);
	assert(height > 0.0);

	if (!do_split && not_too_small && (rgen.rand()&3) < (was_cube ? 2 : 3)) {
		// oddly shaped multi-sided overlapping sections (50% chance for cube buildings and 75% chance for others)
		point const llc(base.get_llc()), sz(base.get_size());
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
			if (i == 0 || !is_simple_cube()) {parts.push_back(bc); continue;} // no splitting
			split_cubes_recur(bc, parts, 0, parts.size()); // split this cube against all previously added cubes and remove overlapping areas
		} // for i
		if (!parts.empty()) { // at least one part was placed; should (almost) always be true?
			parts.shrink_to_fit(); // optional
			std::reverse(parts.begin(), parts.end()); // highest part should be last so that it gets the roof details
			calc_bcube_from_parts(); // update bcube
			gen_details(rgen, 1);
			end_add_parts();
			gen_interior(rgen, 1);
			gen_building_doors_if_needed(rgen);
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

	if (do_split) { // generate L, T, or U shape
		cube_t const split_cube(parts.back());
		parts.pop_back();
		split_in_xy(split_cube, rgen);
		if (num_levels <= 3) {gen_details(rgen, 0);}
	}
	else {
		if ((rgen.rand()&3) != 0) {maybe_add_special_roof(rgen);} // 67% chance
		if (num_levels <= 3) {gen_details(rgen, 1);}
	}
	end_add_parts();
	gen_interior(rgen, 0);
	gen_building_doors_if_needed(rgen);
}

void building_t::split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen) {

	// generate L, T, U, H, +, O shape
	point const llc(seed_cube.get_llc()), sz(seed_cube.get_size());
	bool const allow_courtyard(seed_cube.dx() < 1.6*seed_cube.dy() && seed_cube.dy() < 1.6*seed_cube.dx()); // AR < 1.6:1
	int const shape(rand()%(allow_courtyard ? 10 : 9)); // 0-9
	has_courtyard = (shape == 9);
	bool const is_hpo(shape >= 7);
	bool const dim(rgen.rand_bool()); // {x,y}
	bool const dir(is_hpo ? 1 : rgen.rand_bool()); // {neg,pos} - H/+/O shapes are symmetric and always pos
	float const smin(0.2), smax(has_courtyard ? 0.35 : 0.4); // outer and inner split sizes
	float const div(is_hpo ? rgen.rand_uniform(smin, smax) : rgen.rand_uniform(0.3, 0.7)); // split pos in 0-1 range
	float const s1(rgen.rand_uniform(smin, smax)), s2(rgen.rand_uniform(1.0-smax, 1.0-smin)); // split pos in 0-1 range
	float const dpos(llc[dim] + div*sz[dim]), spos1(llc[!dim] + s1*sz[!dim]), spos2(llc[!dim] + s2*sz[!dim]); // split pos in cube space
	float const dpos2(llc[dim] + (1.0 - div)*sz[dim]); // other end - used for H, +, and O
	unsigned const start(parts.size()), num((shape >= 9) ? 4 : ((shape >= 6) ? 3 : 2));
	parts.resize(start+num, seed_cube);
	parts[start+0].d[dim][ dir] = dpos; // full width part (except +)
	parts[start+1].d[dim][!dir] = dpos; // partial width part (except +)

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

cube_t building_t::place_door(cube_t const &base, bool dim, bool dir, float door_height, float door_center,
	float door_pos, float door_center_shift, float width_scale, bool can_fail, rand_gen_t &rgen)
{
	float const door_width(width_scale*door_height), door_half_width(0.5*door_width);
	float const door_shift(0.01*get_material().get_floor_spacing());
	bool const calc_center(door_center == 0.0); // door not yet calculated
	bool const centered(door_center_shift == 0.0 || hallway_dim == (uint8_t)dim); // center doors connected to primary hallways
	cube_t door;
	door.z1() = base.z1(); // same bottom as house
	door.z2() = door.z1() + door_height;

	for (unsigned n = 0; n < 10; ++n) { // make up to 10 tries to place a valid door
		if (calc_center) { // add door to first part of house/building
			float const offset(centered ? 0.5 : rgen.rand_uniform(0.5-door_center_shift, 0.5+door_center_shift));
			door_center = offset*base.d[!dim][0] + (1.0 - offset)*base.d[!dim][1];
			door_pos    = base.d[dim][dir];
		}
		if (interior && hallway_dim == 2) { // not on a hallway
			vect_cube_t const &walls(interior->walls[!dim]); // perpendicular to door
			float const door_lo(door_center - 1.2*door_half_width), door_hi(door_center + 1.2*door_half_width); // pos along wall with a small expand
			float const dpos_lo(door_pos    -     door_half_width), dpos_hi(door_pos    +     door_half_width); // expand width of the door

			for (auto w = walls.begin(); w != walls.end(); ++w) {
				if (w->d[ dim][0] > dpos_hi || w->d[ dim][1] < dpos_lo) continue; // not ending at same wall as door
				if (w->d[!dim][0] > door_hi || w->d[!dim][1] < door_lo) continue; // not intersecting door
				// Note: since we know that all rooms are wider than the door width, we know that we have space for a door on either side of the wall
				float const lo_dist(w->d[!dim][0] - door_lo), hi_dist(door_hi - w->d[!dim][1]);
				if (lo_dist < hi_dist) {door_center += lo_dist;} else {door_center -= hi_dist;} // move the door so that it doesn't open into the end of the wall
				break;
			}
		}
		door.d[ dim][!dir] = door_pos + door_shift*(dir ? 1.0 : -1.0); // move slightly away from the house to prevent z-fighting
		door.d[ dim][ dir] = door.d[dim][!dir]; // make zero size in this dim
		door.d[!dim][0] = door_center - door_half_width; // left
		door.d[!dim][1] = door_center + door_half_width; // right

		// we're free to choose the door pos, and have the interior, so we can check if the door is in a good location;
		// if we get here on the last iteration, just keep the door even if it's an invalid location
		if (calc_center && interior && !centered) {
			if (interior->is_blocked_by_stairs_or_elevator(door, door_width)) {
				if (can_fail && n == 9) return cube_t(); // last iteration, fail
				continue; // try a new location
			}
		}
		break; // done
	} // for n
	return door;
}

void building_t::gen_house(cube_t const &base, rand_gen_t &rgen) {

	assert(parts.empty());
	int const type(rgen.rand()%3); // 0=single cube, 1=L-shape, 2=two-part
	bool const two_parts(type != 0);
	unsigned force_dim[2] = {2}; // force roof dim to this value, per part; 2 = unforced/auto
	bool skip_last_roof(0);
	num_sides = 4;
	parts.reserve(two_parts ? 5 : 2); // two house sections + porch roof + porch support + chimney (upper bound)
	parts.push_back(base);
	// add a door
	bool const gen_door(global_building_params.windows_enabled());
	float const door_width_scale(0.5);
	float door_height(get_door_height()), door_center(0.0), door_pos(0.0), dist1(0.0), dist2(0.0);;
	bool door_dim(rgen.rand_bool()), door_dir(0), dim(0), dir(0), dir2(0);
	unsigned door_part(0), detail_type(0);
	real_num_parts = (two_parts ? 2 : 1); // only walkable parts: excludes shed, garage, porch roof, and chimney
	cube_t door_cube;

	if (two_parts) { // multi-part house; parts[1] is the lower height part
		parts.push_back(base); // add second part
		dir = rgen.rand_bool(); // in dim
		float const split(rgen.rand_uniform(0.4, 0.6)*(dir  ? -1.0 : 1.0));
		float delta_height(0.0), shrink[2] = {0.0};

		if (type == 1) { // L-shape
			dir2         = rgen.rand_bool(); // in !dim
			dim          = rgen.rand_bool();
			shrink[dir2] = rgen.rand_uniform(0.4, 0.6)*(dir2 ? -1.0 : 1.0);
			delta_height = max(0.0f, rand_uniform(-0.1, 0.5));
		}
		else if (type == 2) { // two-part
			dim          = get_largest_xy_dim(base); // choose longest dim
			delta_height = rand_uniform(0.1, 0.5);

			for (unsigned d = 0; d < 2; ++d) {
				if (rgen.rand_bool()) {shrink[d] = rgen.rand_uniform(0.2, 0.35)*(d ? -1.0 : 1.0);}
			}
		}
		vector3d const sz(base.get_size());
		parts[0].d[ dim][ dir] += split*sz[dim]; // split in dim
		parts[1].d[ dim][!dir]  = parts[0].d[dim][dir];
		cube_t const pre_shrunk_p1(parts[1]); // save for use in details below
		for (unsigned d = 0; d < 2; ++d) {parts[1].d[!dim][d] += shrink[d]*sz[!dim];} // shrink this part in the other dim
		parts[1].z2() -= delta_height*parts[1].dz(); // lower height
		//ao_bcz2 = parts[1].z2(); // set this lower zval as the base AO so that AO values line up with both parts and the roof triangles
		if (ADD_BUILDING_INTERIORS) {adjust_part_zvals_for_floor_spacing(parts[1]);}
		if (type == 1 && rgen.rand_bool()) {force_dim[0] = !dim; force_dim[1] = dim;} // L-shape, half the time
		else if (type == 2) {force_dim[0] = force_dim[1] = dim;} // two-part - force both parts to have roof along split dim
		detail_type = ((type == 1) ? (rgen.rand()%3) : 0); // 0=none, 1=porch, 2=detatched garage/shed
		door_dir = ((door_dim == dim) ? dir : dir2); // if we have a porch/shed/garage, put the door on that side
		if (door_dim == dim && detail_type == 0) {door_dir ^= 1;} // put it on the opposite side so that the second part isn't in the way

		if (detail_type != 0) { // add details to L-shaped house
			cube_t c(pre_shrunk_p1);
			c.d[!dim][!dir2] = parts[1].d[!dim][dir2]; // other half of the shrunk part1
			dist1 = (c.d[!dim][!dir2] - base.d[!dim][dir2])*rgen.rand_uniform(0.4, 0.6);
			dist2 = (c.d[ dim][!dir ] - base.d[ dim][dir ])*rgen.rand_uniform(0.4, 0.6);
			float const base_dz(parts[1].dz()), height(min(base_dz, max(door_height/0.95f, rgen.rand_uniform(0.55, 0.7)*base_dz)));

			if (gen_door) { // add door in interior of L, centered under porch roof (if it exists, otherwise where it would be)
				door_cube   = c;
				door_center = door_cube.get_center_dim(!door_dim) + 0.5f*((door_dim == dim) ? dist1 : dist2);
				door_pos    = door_cube.d[door_dim][!door_dir];
				door_part   = ((door_dim == dim) ? 0 : 1); // which part the door is connected to
				min_eq(door_height, 0.95f*height);
			}
			if (detail_type == 1) { // porch
				float const width(0.05f*(fabs(dist1) + fabs(dist2))); // width of support pillar
				c.d[!dim][dir2 ] += dist1; // move away from bcube edge
				c.d[ dim][ dir ] += dist2; // move away from bcube edge
				c.d[!dim][!dir2] -= 0.001*dist1; // adjust slightly so it's not exactly adjacent to the house and won't be considered internal face removal logic
				c.d[ dim][ !dir] -= 0.001*dist2;
				c.z1() += height; // move up
				c.z2()  = c.z1() + 0.05*parts[1].dz();
				parts.push_back(c); // porch roof
				c.z2() = c.z1();
				c.z1() = pre_shrunk_p1.z1(); // support pillar
				c.d[!dim][!dir2] = c.d[!dim][dir2] + (dir2 ? -1.0 : 1.0)*width;
				c.d[ dim][!dir ] = c.d[ dim][ dir] + (dir  ? -1.0 : 1.0)*width;
				skip_last_roof = 1;
			}
			else if (detail_type == 2) { // detatched garage/shed
				c.d[!dim][dir2 ]  = base.d[!dim][dir2]; // shove it into the opposite corner of the bcube
				c.d[ dim][dir  ]  = base.d[ dim][dir ]; // shove it into the opposite corner of the bcube
				c.d[!dim][!dir2] -= dist1; // move away from bcube edge
				c.d[ dim][!dir ] -= dist2; // move away from bcube edge
				c.z2() = c.z1() + min(min(c.dx(), c.dy()), height); // no taller than x or y size; Note: z1 same as part1
				vector3d const car_sz(get_nom_car_size());
				bool const is_garage(max(c.dx(), c.dy()) > 1.2f*car_sz.x && min(c.dx(), c.dy()) > 1.2f*car_sz.y && c.dz() > 1.2f*car_sz.z); // must be able to fit a car
				(is_garage ? has_garage : has_shed) = 1;
				// add a door
				bool const long_dim(c.dx() < c.dy());
				bool door_dir(c.d[long_dim][0] == bcube.d[long_dim][0]); // interior face
				if (is_garage) {door_dir ^= 1;} // facing away from the house, not toward it
				float const wscale(is_garage ? 0.9*c.get_sz_dim(!long_dim)/door_height : door_width_scale);
				add_door(place_door(c, long_dim, door_dir, door_height, 0.0, 0.0, 0.0, wscale, 0, rgen), parts.size(), long_dim, door_dir, 0);
				if (is_garage) {doors.back().type = tquad_with_ix_t::TYPE_GDOOR;} // make it a garage door rather than a house door
			}
			parts.push_back(c); // support column or shed/garage
		} // end house details
		else if (type == 1 && has_city_trees()) { // place a tree for L-shaped house with no garage or shed
			vect_cube_t cubes;
			cube_t empty_space(bcube);

			for (unsigned e = 0; e < 2; ++e) {
				cubes.clear();
				subtract_cube_from_cube(empty_space, parts[e], cubes);
				assert(cubes.size() == 1);
				empty_space = cubes[0];
			}
			tree_pos = empty_space.get_cube_center(); // centered on where the garage/shed would have been
			vector3d const dir(tree_pos - bcube.get_cube_center());
			tree_pos  += (0.05f*(bcube.dx() + bcube.dy())/dir.mag())*dir; // shift slightly away from house center so that it's less likely to intersect the house
			tree_pos.z = bcube.z1();
		}
		calc_bcube_from_parts(); // maybe calculate a tighter bounding cube
	} // end type != 0  (multi-part house)
	else if (gen_door) { // single cube house
		door_dir  = rgen.rand_bool(); // select a random dir
		door_part = 0; // only one part
	}
	gen_interior(rgen, 0); // before adding door

	if (gen_door) {
		cube_t door(place_door(parts[door_part], door_dim, door_dir, door_height, door_center, door_pos, 0.25, door_width_scale, 0, rgen));

		if (interior && interior->is_blocked_by_stairs_or_elevator(door, 0.5*door_height)) { // blocked by stairs - switch door to other side/dim
			door_dim ^= 1;
			door_dir  = (two_parts ? ((door_dim == dim) ? dir : dir2) : (door_dir^1)); // if we have a porch/shed/garage, put the door on that side
			if (door_dim == dim && two_parts && detail_type == 0) {door_dir ^= 1;} // put it on the opposite side so that the second part isn't in the way

			if (door_center != 0.0) { // reclaculate for L-shaped house
				door_center = door_cube.get_center_dim(!door_dim) + 0.5f*((door_dim == dim) ? dist1 : dist2);
				door_pos    = door_cube.d[door_dim][!door_dir];
				door_part   = ((door_dim == dim) ? 0 : 1); // which part the door is connected to
			}
			door = place_door(parts[door_part], door_dim, door_dir, door_height, door_center, door_pos, 0.25, door_width_scale, 0, rgen); // keep door_height
		}
		add_door(door, door_part, door_dim, door_dir, 0);
		if (doors.size() == 2) {swap(doors[0], doors[1]);} // make sure the house door comes before the garage/shed door
	}
	float const peak_height(rgen.rand_uniform(0.15, 0.5)); // same for all parts
	float roof_dz[3] = {0.0f};
	bool hipped_roof[3] = {0};

	for (auto i = parts.begin(); (i + skip_last_roof) != parts.end(); ++i) {
		unsigned const ix(i - parts.begin()), fdim(force_dim[ix]);
		cube_t const &other(two_parts ? parts[1-ix] : *i); // == self for single part houses
		bool const dim((fdim < 2) ? fdim : get_largest_xy_dim(*i)); // use longest side if not forced
		bool const other_dim(two_parts ? ((force_dim[1-ix] < 2) ? force_dim[1-ix] : get_largest_xy_dim(other)) : 0);
		float extend_to(0.0), max_dz(i->dz());

		if (type == 1 && ix == 1 && dim != other_dim && parts[0].z2() == parts[1].z2()) { // same z2, opposite dim T-junction
			max_dz    = peak_height*parts[0].get_sz_dim(!other_dim); // clamp roof zval to other roof's peak
			extend_to = parts[0].get_center_dim(!other_dim); // extend lower part roof to center of upper part roof
		}
		bool can_be_hipped(ix < real_num_parts && extend_to == 0.0 && i->get_sz_dim(dim) > i->get_sz_dim(!dim)); // must be longer dim
		
		if (can_be_hipped && two_parts) {
			float const part_roof_z (i->z2()    + min(peak_height*i->get_sz_dim(!dim), i->dz()));
			float const other_roof_z(other.z2() + min(peak_height*other.get_sz_dim(!other_dim), other.dz()));
			can_be_hipped = (part_roof_z >= other_roof_z); // no hipped for lower part
		}
		bool const hipped(can_be_hipped && rgen.rand_bool()); // hipped roof 50% of the time
		if (hipped) {roof_dz[ix] = gen_hipped_roof(*i, peak_height, extend_to);}
		else {
			unsigned skip_side_tri(2); // default = skip neither
			
			if (two_parts && (dim == other_dim || hipped_roof[1-ix]) && i->d[!dim][0] >= other.d[!dim][0] && i->d[!dim][1] <= other.d[!dim][1] && i->z2() <= other.z2()) {
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
	} // for i
	if ((rgen.rand()%3) != 0) { // add a chimney 67% of the time
		unsigned part_ix(0);

		if (two_parts) { // start by selecting a part (if two parts)
			float const v0(parts[0].get_volume()), v1(parts[1].get_volume());
			if      (v0 > 2.0*v1) {part_ix = 0;} // choose larger part 0
			else if (v1 > 2.0*v0) {part_ix = 1;} // choose larger part 1
			else {part_ix = rgen.rand_bool();} // close in area - choose a random part
		}
		assert(roof_dz[part_ix] > 0.0);
		unsigned const fdim(force_dim[part_ix]);
		cube_t const &part(parts[part_ix]);
		bool const dim((fdim < 2) ? fdim : get_largest_xy_dim(part)); // use longest side if not forced
		bool dir(rgen.rand_bool());
		if (two_parts && part.d[dim][dir] != bcube.d[dim][dir]) {dir ^= 1;} // force dir to be on the edge of the house bcube (not at a point interior to the house)
		cube_t c(part);
		float const sz1(c.get_sz_dim(!dim)), sz2(c.get_sz_dim(dim));
		float shift(0.0);

		if ((rgen.rand()%3) != 0) { // make the chimney non-centered 67% of the time
			shift = sz1*rgen.rand_uniform(0.1, 0.25); // select a shift in +/- (0.1, 0.25) - no small offset from center
			if (rgen.rand_bool()) {shift = -shift;}
		}
		float const center(c.get_center_dim(!dim) + shift);
		float const chimney_dz((hipped_roof[part_ix] ? 0.5 : 1.0)*roof_dz[part_ix]); // lower for hipped roof
		c.d[dim][!dir]  = c.d[dim][ dir] + (dir ? -0.03f : 0.03f)*(sz1 + sz2); // chimney depth
		c.d[dim][ dir] += (dir ? -0.01 : 0.01)*sz2; // slight shift from edge of house to avoid z-fighting
		c.d[!dim][0] = center - 0.05*sz1;
		c.d[!dim][1] = center + 0.05*sz1;
		c.z1()  = c.z2();
		c.z2() += rgen.rand_uniform(1.25, 1.5)*chimney_dz - 0.4f*abs(shift);
		parts.push_back(c);
		// add top quad to cap chimney (will also update bcube to contain chimney)
		tquad_t tquad(4); // quad
		tquad.pts[0].assign(c.x1(), c.y1(), c.z2());
		tquad.pts[1].assign(c.x2(), c.y1(), c.z2());
		tquad.pts[2].assign(c.x2(), c.y2(), c.z2());
		tquad.pts[3].assign(c.x1(), c.y2(), c.z2());
		roof_tquads.emplace_back(tquad, (unsigned)tquad_with_ix_t::TYPE_CCAP); // tag as chimney cap
		has_chimney = 1;
	}
	roof_type = ROOF_TYPE_PEAK; // peaked and hipped roofs are both this type
	add_roof_to_bcube();
	gen_grayscale_detail_color(rgen, 0.4, 0.8); // for roof
}

tquad_with_ix_t set_door_from_cube(cube_t const &c, bool dim, bool dir, unsigned type, float pos_adj, bool opened, bool opens_out, bool swap_sides) {

	tquad_with_ix_t door(4, type); // quad
	float const pos(c.d[dim][0] + (opened ? 0.0 : pos_adj*(dir ? 1.0 : -1.0))); // move away from wall slightly (not needed if opened)
	door.pts[0].z = door.pts[1].z = c.z1(); // bottom
	door.pts[2].z = door.pts[3].z = c.z2(); // top
	door.pts[0][!dim] = door.pts[3][!dim] = c.d[!dim][ dir]; //  dir side
	door.pts[1][!dim] = door.pts[2][!dim] = c.d[!dim][!dir]; // !dir side
	door.pts[0][ dim] = door.pts[1][ dim] = door.pts[2][dim] = door.pts[3][dim] = pos;
	if (dim == 0) {swap(door.pts[0], door.pts[1]); swap(door.pts[2], door.pts[3]);} // swap two corner points to flip windowing dir and invert normal for doors oriented in X

	if (opened) { // rotate 90 degrees about pts[0]/pts[3] - change pts[1]/pts[2]; this is just a placeholder for now
		float const width(c.get_sz_dim(!dim)), offset(0.01*width*((dir ^ dim) ? 1.0 : -1.0)); // move slightly away from the wall to prevent z-fighting
		if (swap_sides) {door.pts[0][!dim] = door.pts[3][!dim] = door.pts[1][!dim] + offset;}
		else            {door.pts[1][!dim] = door.pts[2][!dim] = door.pts[0][!dim] + offset;}
		door.pts[1][ dim] = door.pts[2][ dim] = pos + width*((dir ^ opens_out) ? 1.0 : -1.0);
	}
	return door;
}

bool building_t::add_door(cube_t const &c, unsigned part_ix, bool dim, bool dir, bool for_building) { // exterior doors

	if (c.is_all_zeros()) return 0;
	vector3d const sz(c.get_size());
	assert(sz[dim] == 0.0 && sz[!dim] > 0.0 && sz.z > 0.0);
	unsigned const type(for_building ? (unsigned)tquad_with_ix_t::TYPE_BDOOR : (unsigned)tquad_with_ix_t::TYPE_HDOOR);
	doors.push_back(set_door_from_cube(c, dim, dir, type, 0.01*sz[!dim], 0, 0, 0)); // opened=0, opens_out=0, swap_sides=0
	if (part_ix < 4) {door_sides[part_ix] |= 1 << (2*dim + dir);}
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
	unsigned const extend_dir(extend_roof(top, extend_to, dim));
	float const width(top.get_sz_dim(!dim)), roof_dz(min(max_dz, min(peak_height*width, top.dz())));
	float const z1(top.z2()), z2(z1 + roof_dz), x1(top.x1()), y1(top.y1()), x2(top.x2()), y2(top.y2());
	assert(roof_dz > 0.0);
	point pts[6] = {point(x1, y1, z1), point(x1, y2, z1), point(x2, y2, z1), point(x2, y1, z1), point(x1, y1, z2), point(x2, y2, z2)};
	if (dim == 0) {pts[4].y = pts[5].y = 0.5f*(y1 + y2);} // yc
	else          {pts[4].x = pts[5].x = 0.5f*(x1 + x2);} // xc
	unsigned const qixs[2][2][4] = {{{0,3,5,4}, {4,5,2,1}}, {{0,4,5,1}, {4,3,2,5}}}; // 2 quads
	roof_tquads.reserve(roof_tquads.size() + (3 + (extend_dir == 2))); // 2 roof quads + 1-2 side triangles
	unsigned const tixs[2][2][3] = {{{1,0,4}, {3,2,5}}, {{0,3,4}, {2,1,5}}}; // 2 triangles

	for (unsigned n = 0; n < 2; ++n) { // triangle section/wall from z1 up to roof
		if (n == extend_dir || n == skip_side_tri) continue; // skip this side
		bool skip(0);

		// exclude tquads contained in/adjacent to other parts, considering only the cube parts;
		// yes, a triangle side can be occluded by a cube + another opposing triangle side from a higher wall of the house, but it's uncommon, complex, and currently ignored
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
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
		roof_tquads.emplace_back(tquad, (unsigned)tquad_with_ix_t::TYPE_ROOF); // tag as roof
	}
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
	roof_tquads.resize(start_ix + 4); // defaults to TYPE_ROOF

	for (unsigned n = 0; n < 2; ++n) {
		unsigned const ix(start_ix + n);
		roof_tquads[ix+0].npts = 4; // quads
		UNROLL_4X(roof_tquads[ix+0].pts[i_] = pts[qixs[n][i_]];)
		roof_tquads[ix+2].npts = 3; // triangles
		UNROLL_3X(roof_tquads[ix+2].pts[i_] = pts[tixs[n][i_]];)
	}
	return roof_dz;
}

void building_t::gen_building_doors_if_needed(rand_gen_t &rgen) {

	if (!is_cube()) return; // for now, only cube buildings can have doors; doors can be added to N-gon (non cylinder) buildings later
	assert(!parts.empty());
	float const door_height(1.1*get_door_height()), wscale(0.7); // a bit taller and a lot wider than house doors

	if (hallway_dim < 2) { // building has primary hallway, place doors at both ends of first part
		for (unsigned d = 0; d < 2; ++d) {
			add_door(place_door(parts.front(), bool(hallway_dim), d, door_height, 0.0, 0.0, 0.0, wscale, 0, rgen), 0, bool(hallway_dim), d, 1);
		}
		return;
	}
	bool const pref_dim(rgen.rand_bool()), pref_dir(rgen.rand_bool());
	bool const has_windows(get_material().add_windows);
	bool used[4] = {0,0,0,0}; // per-side, not per-base cube
	unsigned const num_doors(1 + (rgen.rand() % (has_windows ? 3 : 4))); // 1-4; buildings with windows have at most 3 doors since they're smaller

	for (unsigned num = 0; num < num_doors; ++num) {
		bool placed(0);

		for (auto b = parts.begin(); b != get_real_parts_end() && !placed; ++b) { // try all different ground floor parts
			unsigned const part_ix(b - parts.begin());
			if (has_windows && part_ix >= 4) break; // only first 4 parts can have doors - must match first floor window removal logic
			if (b->z1() > bcube.z1()) break; // moved off the ground floor - done

			for (unsigned n = 0; n < 4; ++n) {
				bool const dim(pref_dim ^ bool(n>>1)), dir(pref_dir ^ bool(n&1));
				if (b->d[dim][dir] != bcube.d[dim][dir]) continue; // find a side on the exterior to ensure door isn't obstructed by a building cube
				if (used[2*dim + dir]) continue; // door already placed on this side
				used[2*dim + dir] = 1; // mark used
				bool const allow_fail(!doors.empty() || b+1 != get_real_parts_end()); // allow door placement to fail if we've already placed at least one door of not last part
				placed = add_door(place_door(*b, dim, dir, door_height, 0.0, 0.0, 0.1, wscale, allow_fail, rgen), part_ix, dim, dir, 1);
				break;
			} // for n
		} // for b
	} // for num
	if (has_courtyard) { // add a door opening into the courtyard
		assert(parts.size() >= 4);
		unsigned const pref_side(rgen.rand()&3);
		bool placed(0);

		for (unsigned p = 0; p < 4 && !placed; ++p) {
			unsigned const part_ix((p + pref_side) & 3);
			cube_t const &part(parts[part_ix]);
			if (part.z1() > bcube.z1()) continue; // not on the ground floor

			for (unsigned n = 0; n < 4; ++n) {
				bool const dim(n>>1), dir(n&1);
				if (part.d[dim][dir] == bcube.d[dim][dir] || part.d[dim][!dir] != bcube.d[dim][!dir]) continue; // find a side on the interior
				placed = add_door(place_door(part, dim, dir, door_height, 0.0, 0.0, 0.0, wscale, 1, rgen), part_ix, dim, dir, 1); // centered
			}
		} // for b
	}
}

void building_t::place_roof_ac_units(unsigned num, float sz_scale, cube_t const &bounds, vect_cube_t const &avoid, bool avoid_center, rand_gen_t &rgen) {

	if (num == 0) return;
	vector3d cube_sz(sz_scale, 1.5*sz_scale, 1.2*sz_scale); // consistent across all AC units (note that x and y will be doubled)
	if (rgen.rand_bool()) {swap(cube_sz.x, cube_sz.y);} // other orientation

	for (unsigned i = 0; i < num; ++i) {
		roof_obj_t c(ROOF_OBJ_AC);
		bool placed(0);

		for (unsigned n = 0; n < 100 && !placed; ++n) { // limited to 100 attempts to prevent infinite loop
			c.set_from_point(point(rgen.rand_uniform(bounds.x1(), bounds.x2()), rgen.rand_uniform(bounds.y1(), bounds.y2()), bounds.z2()));
			c.expand_by_xy(cube_sz);
			if (!bounds.contains_cube_xy(c)) continue; // not contained
			c.z2() += cube_sz.z; // z2
			if (has_bcube_int(c, avoid, 0) || has_bcube_int(c, details, 0)) continue; // intersects avoid cubes or other detail objects (inc_adj=0)
			if (avoid_center && c.contains_pt_xy(bounds.get_cube_center())) continue; // intersects antenna
			placed = 1;
		} // for n
		if (!placed) break; // failed, exit loop
		if ((rgen.rand() & 7) == 0) {swap(cube_sz.x, cube_sz.y);} // swap occasionally
		details.push_back(c);
	} // for n
}

void building_t::add_roof_walls(cube_t const &c, float wall_width, bool overlap_corners, cube_t out[4]) {
	// add walls around the roof; we don't need to draw all sides of all cubes, but it's probably not worth the trouble to sort it out
	float const height(2.0f*wall_width), width(wall_width), shorten(overlap_corners ? 0.0f : wall_width), z1(c.z2()), z2(c.z2()+height);
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

void building_t::gen_details(rand_gen_t &rgen, bool is_rectangle) { // for the roof

	bool const flat_roof(roof_type == ROOF_TYPE_FLAT), add_walls(is_simple_cube() && flat_roof); // simple cube buildings with flat roofs
	unsigned const num_ac_units((flat_roof && is_cube() && !is_rotated()) ? (rgen.rand() % 7) : 0); // cube buildings only for now
	float const window_vspacing(get_material().get_floor_spacing()), wall_width(0.2*window_vspacing);
	assert(!parts.empty());

	if (!is_rectangle) { // polygon roof, can only add AC units
		if (add_walls) {
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
				if (!i->is_all_zeros()) {details.emplace_back(*i, ROOF_OBJ_WALL);}
			}
		}
		if (num_ac_units > 0) {
			float const xy_sz(0.75*bcube.get_size().xy_mag()*rgen.rand_uniform(0.012, 0.02)); // size based on bcube

			for (auto p = parts.begin(); p != parts.end(); ++p) {
				unsigned const num_this_part(min(num_ac_units, unsigned(num_ac_units*p->dx()*p->dy()/(bcube.dx()*bcube.dy()) + 1))); // distribute based on area
				place_roof_ac_units(num_this_part, xy_sz, *p, parts, 0, rgen);
			}
		}
		details.shrink_to_fit();
		return;
	}
	unsigned const num_blocks(flat_roof ? (rgen.rand() % 9) : 0); // 0-8; 0 if there are roof quads (houses, etc.)
	bool const add_antenna((flat_roof || roof_type == ROOF_TYPE_SLOPE) && (rgen.rand() & 1));
	unsigned const num_details(num_blocks + num_ac_units + 4*add_walls + add_antenna);
	if (num_details == 0) return; // nothing to do
	cube_t const &top(parts.back()); // top/last part
	if (add_walls && min(top.dx(), top.dy()) < 4.0*wall_width) return; // too small
	float const xy_sz(top.get_size().xy_mag()); // better to use bcube for size?
	cube_t bounds(top);
	if (add_walls) {bounds.expand_by_xy(-wall_width);}
	details.reserve(details.size() + num_details);
	vector<point> points; // reused across calls

	for (unsigned i = 0; i < num_blocks; ++i) {
		roof_obj_t c(ROOF_OBJ_BLOCK); // generic block
		float const height_scale(0.0035f*(top.dz() + bcube.dz())); // based on avg height of current section and entire building
		float const height(height_scale*rgen.rand_uniform(1.0, 4.0));
		bool placed(0);

		for (unsigned n = 0; n < 100 && !placed; ++n) { // limited to 100 attempts to prevent infinite loop
			c.set_from_point(point(rgen.rand_uniform(bounds.x1(), bounds.x2()), rgen.rand_uniform(bounds.y1(), bounds.y2()), top.z2()));
			c.expand_by(vector3d(xy_sz*rgen.rand_uniform(0.01, 0.07), xy_sz*rgen.rand_uniform(0.01, 0.07), 0.0));
			if (!bounds.contains_cube_xy(c)) continue; // not contained
			placed = 1;
			if (is_simple_cube()) break; // success/done

			for (unsigned j = 0; j < 4; ++j) { // check cylinder/ellipse
				point const pt(c.d[0][j&1], c.d[1][j>>1], 0.0); // XY only
				if (!check_part_contains_pt_xy(top, pt, points)) {placed = 0; break;}
			}
		} // for n
		if (!placed) break; // failed, exit loop
		c.z2() += height; // z2
		details.push_back(c);
	} // for i
	if (num_ac_units > 0) {place_roof_ac_units(num_ac_units, xy_sz*rgen.rand_uniform(0.012, 0.02), bounds, vect_cube_t(), add_antenna, rgen);}

	if (add_walls) {
		cube_t cubes[4];
		add_roof_walls(top, wall_width, 0, cubes); // overlap_corners=0
		for (unsigned i = 0; i < 4; ++i) {details.emplace_back(cubes[i], ROOF_OBJ_WALL);}
	}
	if (add_antenna) { // add antenna
		float const radius(0.003f*rgen.rand_uniform(1.0, 2.0)*(top.dx() + top.dy()));
		float const height(rgen.rand_uniform(0.25, 0.5)*top.dz());
		roof_obj_t antenna(ROOF_OBJ_ANT);
		antenna.set_from_point(top.get_cube_center());
		antenna.expand_by(vector3d(radius, radius, 0.0));
		antenna.z1() = top.z2(); // z1
		antenna.z2() = bcube.z2() + height; // z2 (use bcube to include sloped roof)
		details.push_back(antenna);
	}
	for (auto i = details.begin(); i != details.end(); ++i) {assert(i->is_strictly_normalized()); max_eq(bcube.z2(), i->z2());} // extend bcube z2 to contain details
	if (roof_type == ROOF_TYPE_FLAT) {gen_grayscale_detail_color(rgen, 0.2, 0.6);} // for antenna and roof
}

void building_t::maybe_add_special_roof(rand_gen_t &rgen) {
	assert(!parts.empty());
	cube_t const &top(parts.back());
	vector3d const sz(top.get_size()); // top/last part

	if (global_building_params.onion_roof && num_sides >= 16 && flat_side_amt == 0.0) { // cylinder building
		if (sz.x < 1.2*sz.y && sz.y < 1.2*sz.x && sz.z > max(sz.x, sz.y)) {roof_type = ROOF_TYPE_ONION;}
	}
	else if (is_simple_cube()) { // only simple cubes are handled
		if (global_building_params.dome_roof && sz.x < 1.2*sz.y && sz.y < 1.2*sz.x && sz.z > max(sz.x, sz.y)) {roof_type = ROOF_TYPE_DOME;} // roughly square, not too short
		else {gen_sloped_roof(rgen, top);} // sloped roof
	}
	if      (roof_type == ROOF_TYPE_DOME ) {max_eq(bcube.z2(), (top.z2() + 0.5f*(sz.x, sz.y)));}
	else if (roof_type == ROOF_TYPE_ONION) {max_eq(bcube.z2(), (top.z2() + 1.0f*(sz.x, sz.y)));}
}
void building_t::gen_sloped_roof(rand_gen_t &rgen, cube_t const &top) { // Note: currently not supported for rotated buildings

	float const peak_height(rgen.rand_uniform(0.2, 0.5));
	float const wmin(min(top.dx(), top.dy())), z1(top.z2()), z2(z1 + peak_height*wmin), x1(top.x1()), y1(top.y1()), x2(top.x2()), y2(top.y2());
	point const pts[5] = {point(x1, y1, z1), point(x1, y2, z1), point(x2, y2, z1), point(x2, y1, z1), point(0.5f*(x1 + x2), 0.5f*(y1 + y2), z2)};
	float const d1(rgen.rand_uniform(0.0, 0.8));

	if (d1 < 0.2) { // pointed roof with 4 sloped triangles
		unsigned const ixs[4][3] = {{1,0,4}, {3,2,4}, {0,3,4}, {2,1,4}};
		roof_tquads.resize(4); // defaults to TYPE_ROOF

		for (unsigned n = 0; n < 4; ++n) {
			roof_tquads[n].npts = 3; // triangles
			UNROLL_3X(roof_tquads[n].pts[i_] = pts[ixs[n][i_]];)
		}
	}
	else { // flat roof with center quad and 4 surrounding sloped quads
		point const center((1.0 - d1)*pts[4]);
		point pts2[8];
		for (unsigned n = 0; n < 4; ++n) {pts2[n] = pts[n]; pts2[n+4] = d1*pts[n] + center;}
		unsigned const ixs[5][4] = {{4,7,6,5}, {0,4,5,1}, {3,2,6,7}, {0,3,7,4}, {2,1,5,6}}; // add the flat quad first, which works better for sphere intersections
		roof_tquads.resize(5); // defaults to TYPE_ROOF

		for (unsigned n = 0; n < 5; ++n) {
			roof_tquads[n].npts = 4; // quads
			UNROLL_4X(roof_tquads[n].pts[i_] = pts2[ixs[n][i_]];)
		}
	}
	roof_type = ROOF_TYPE_SLOPE;
	add_roof_to_bcube(); // is max_eq(bcube.z2(), z2) good enough?
	gen_grayscale_detail_color(rgen, 0.4, 0.8); // for antenna and roof
}

void building_t::add_roof_to_bcube() {
	for (auto i = roof_tquads.begin(); i != roof_tquads.end(); ++i) {i->update_bcube(bcube);} // technically should only need to update z2
}
void building_t::gen_grayscale_detail_color(rand_gen_t &rgen, float imin, float imax) {
	float const cscale(rgen.rand_uniform(imin, imax));
	detail_color = colorRGBA(cscale, cscale, cscale, 1.0);
}


// *** Interiors ***

bool building_t::add_table_and_chairs(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id,
	point const &place_pos, colorRGBA const &chair_color, float rand_place_off, float tot_light_amt, bool is_lit)
{
	float const window_vspacing(get_window_vspace());
	uint8_t const obj_flags(is_lit ? RO_FLAG_LIT : 0);
	vector3d const room_sz(room.get_size());
	vector<room_object_t> &objs(interior->room_geom->objs);
	point table_pos(place_pos);
	vector3d table_sz;
	for (unsigned d = 0; d < 2; ++d) {table_sz [d]  = 0.18*window_vspacing*(1.0 + rgen.rand_float());} // half size relative to window_vspacing
	for (unsigned d = 0; d < 2; ++d) {table_pos[d] += rand_place_off*room_sz[d]*rgen.rand_uniform(-1.0, 1.0);} // near the center of the room
	point llc(table_pos - table_sz), urc(table_pos + table_sz);
	llc.z = table_pos.z; // bottom
	urc.z = table_pos.z + 0.2*window_vspacing;
	cube_t table(llc, urc);
	if (!is_valid_placement_for_room(table, room, blockers)) return 0; // check proximity to doors
	objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, obj_flags, tot_light_amt);
	float const chair_sz(0.1*window_vspacing); // half size

	// place some chairs around the table
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (rgen.rand_bool()) continue; // 50% of the time
			point chair_pos(table_pos); // same starting center and z1
			chair_pos[dim] += (dir ? -1.0f : 1.0f)*(rgen.rand_uniform(-0.5, 1.2)*chair_sz + table_sz[dim]);
			cube_t chair(chair_pos, chair_pos);
			chair.z2() += 0.4*window_vspacing; // chair height
			chair.expand_by(vector3d(chair_sz, chair_sz, 0.0));
			if (!is_valid_placement_for_room(chair, room, blockers)) continue; // check proximity to doors
			objs.emplace_back(chair, TYPE_CHAIR, room_id, dim, dir, obj_flags, tot_light_amt);
			objs.back().color = chair_color;
		} // for dir
	} // for dim
	return 1;
}

void set_light_xy(cube_t &light, point const &center, float light_size, bool light_dim, room_obj_shape light_shape) {
	for (unsigned dim = 0; dim < 2; ++dim) {
		float const sz(((light_shape == SHAPE_CYLIN) ? 1.6 : ((bool(dim) == light_dim) ? 2.2 : 1.0))*light_size);
		light.d[dim][0] = center[dim] - sz;
		light.d[dim][1] = center[dim] + sz;
	}
	light.z1() = light.z2() = center.z; // set so that valid pos can be checked
}

template<typename T> bool has_bcube_int_exp(cube_t const &bcube, vector<T> const &bcubes, float expand) {
	cube_t bcube_exp(bcube);
	bcube_exp.expand_by(expand); // expand in all dirs, including z
	return has_bcube_int(bcube_exp, bcubes);
}

// Note: these three floats can be calculated from mat.get_floor_spacing(), but it's easier to change the constants if we just pass them in
void building_t::gen_room_details(rand_gen_t &rgen, vect_cube_t const &ped_bcubes) {

	assert(interior);
	if (interior->room_geom) return; // already generated?
	//timer_t timer("Gen Room Details");
	interior->room_geom.reset(new building_room_geom_t);
	vector<room_object_t> &objs(interior->room_geom->objs);
	float const window_vspacing(get_window_vspace()), floor_thickness(FLOOR_THICK_VAL*window_vspacing), fc_thick(0.5*floor_thickness);
	interior->room_geom->obj_scale = window_vspacing; // used to scale room object textures
	unsigned tot_num_rooms(0), num_light_stacks(0);
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {tot_num_rooms += calc_num_floors(*r, window_vspacing, floor_thickness);}
	objs.reserve(tot_num_rooms); // placeholder - there will be more than this many
	room_obj_shape const light_shape(is_house ? SHAPE_CYLIN : SHAPE_CUBE);

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		float const light_amt(window_vspacing*r->get_light_amt()); // multiply perimeter/area by window spacing to make unitless
		unsigned const num_floors(calc_num_floors(*r, window_vspacing, floor_thickness));
		unsigned const room_id(r - interior->rooms.begin());
		point room_center(r->get_cube_center());
		// determine light pos and size for this stack of rooms
		bool const room_dim(r->dx() < r->dy()); // longer room dim
		float light_size(floor_thickness); // default size for houses

		if (r->is_office) { // light size varies by office size
			float const room_size(r->dx() + r->dy()); // normalized to office size
			light_size = max(0.015f*room_size, 0.67f*floor_thickness);
		}
		if (r->is_hallway) { // light size varies by hallway size
			float const room_size(min(r->dx(), r->dy())); // normalized to hallway width
			light_size = max(0.06f*room_size, 0.67f*floor_thickness);
		}
		float const light_val(22.0*light_size), room_light_intensity(light_val*light_val/r->get_area_xy()); // average for room, unitless
		cube_t pri_light, sec_light;
		set_light_xy(pri_light, room_center, light_size, room_dim, light_shape);
		if (!r->contains_cube_xy(pri_light)) {pri_light.set_to_zeros();} // disable light if it doesn't fit (small room)
		bool const blocked_by_stairs(!r->is_hallway && interior->is_blocked_by_stairs_or_elevator(pri_light, fc_thick));
		bool use_sec_light(0);
		float z(r->z1());

		if (blocked_by_stairs) { // blocked by stairs - see if we can add a light off to the side in the other orient
			bool const first_dir(rgen.rand_bool());

			for (unsigned d = 0; d < 2; ++d) { // see if we can place it by moving on one direction
				point new_center(room_center);
				new_center[room_dim] += ((bool(d) ^ first_dir) ? -1.0 : 1.0)*0.33*r->get_sz_dim(room_dim);
				set_light_xy(sec_light, new_center, light_size, !room_dim, light_shape); // flip the light dim
				if (!interior->is_blocked_by_stairs_or_elevator(sec_light, fc_thick)) {use_sec_light = 1; break;} // add if not blocked
			}
		}
		// make chair colors consistent for each part by using a few variables for a hash
		colorRGBA chair_colors[12] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN};
		colorRGBA chair_color(chair_colors[(13*r->part_id + 123*tot_num_rooms + 617*mat_ix + 1337*num_floors) % 12]);
		unsigned num_lights_added(0);

		// place objects on each floor for this room
		for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
			room_center.z = z + fc_thick; // floor height
			bool const top_floor(f+1 == num_floors), check_stairs(!is_house && parts.size() > 1 && top_floor); // top floor of building that may have stairs connecting to upper stack
			bool is_lit(0), light_dim(room_dim), has_stairs(r->has_stairs);

			if (!has_stairs && (f == 0 || top_floor) && interior->stairwells.size() > 1) { // check for stairwells connecting stacked parts
				for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
					if (!r->contains_cube_xy(*s)) continue; // stairs not in this room
					// Note: here we adjust stairs zval by floor_thickness to include stairs in the floor but not in the room above
					if (s->z1() + floor_thickness > r->z2()) continue; // stairs above the room
					if (s->z2() + floor_thickness < r->z1()) continue; // stairs below the room
					has_stairs = 1;
				}
			}
			bool const top_of_stairs(has_stairs && top_floor);
			if (top_of_stairs || check_stairs) {num_light_stacks += num_lights_added;} // these cases may shift the light, so we allocate a new stack
			cube_t light;
			if (!blocked_by_stairs || top_of_stairs) {light = pri_light;}
			else if (use_sec_light) {light = sec_light; light_dim ^= 1;}

			if (!light.is_all_zeros()) { // add a light to the center of the ceiling of this room if there's space (always for top of stairs)
				light.z2() = z + window_vspacing - fc_thick;
				light.z1() = light.z2() - 0.5*fc_thick;
				is_lit = (r->is_hallway || ((rgen.rand() & (top_of_stairs ? 3 : 1)) != 0)); // 50% of lights are on, 75% for top of stairs, 100% for hallways

				// check ped_bcubes and set is_lit if any are people are in this floor of this room
				for (auto p = ped_bcubes.begin(); p != ped_bcubes.end() && !is_lit; ++p) {
					if (!p->intersects_xy(*r)) continue; // person not in this room
					if (p->z2() < light.z1() && p->z1() + window_vspacing > light.z2()) {is_lit = 1;} // on this floor
				}
				uint8_t flags(RO_FLAG_NOCOLL); // no collision detection with lights
				if (is_lit)        {flags |= RO_FLAG_LIT;}
				if (top_of_stairs) {flags |= RO_FLAG_TOS;}
				if (has_stairs)    {flags |= RO_FLAG_RSTAIRS;}
				colorRGBA color;
				if (is_house) {color = colorRGBA(1.0, 1.0, 0.85);} // house - yellowish
				else if (r->is_hallway || r->is_office) {color = colorRGBA(0.85, 0.85, 1.0);} // office building - blueish
				else {color = colorRGBA(1.0, 1.0, 1.0);} // white - small office
				unsigned num_lights(r->num_lights);
				
				if (r->is_hallway && num_lights > 1) { // place a light on each side (of the stairs if they exist), and also between stairs and elevator if there are both
					if (r->has_elevator && r->has_stairs) {num_lights = 3;} // we really should have 3 lights in this case
					float const offset(((num_lights == 3) ? 0.3 : 0.2)*r->get_sz_dim(light_dim)); // closer to the ends in the 3 lights case
					cube_t valid_bounds(*r);
					valid_bounds.expand_by_xy(-0.1*window_vspacing); // add some padding

					for (unsigned d = 0; d < num_lights; ++d) {
						float const delta((d == 2) ? 0.0 : (d ? -1.0 : 1.0)*offset); // last light is in the center
						cube_t hall_light(light);
						hall_light.translate_dim(delta, light_dim);

						if (check_stairs && has_bcube_int_exp(hall_light, interior->stairwells, fc_thick)) { // keep moving until not blocked by stairs
							cube_t const hall_light_start(hall_light);
							bool is_valid(0);

							for (unsigned shift_dir = 0; shift_dir < 2 && !is_valid; ++shift_dir) {
								hall_light = hall_light_start;

								for (unsigned n = 0; n < 40; ++n) {
									if (!has_bcube_int_exp(hall_light, interior->stairwells, fc_thick)) {is_valid = 1; break;}
									hall_light.translate_dim(0.04*delta*(shift_dir ? -1.0 : 1.0), light_dim);
									if (!valid_bounds.contains_cube_xy(hall_light)) break; // translated outside the hall, give up
								}
							} // for shift_dir
							if (!is_valid) continue; // skip adding this light
						} // end check_stairs
						objs.emplace_back(hall_light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape, color); // dir=0 (unused)
						objs.back().obj_id = num_light_stacks + d;
					} // for d
					num_lights_added = num_lights;
				}
				else { // normal room
					if (check_stairs && has_bcube_int_exp(light, interior->stairwells, fc_thick)) {is_lit = 0;} // disable if blocked by stairs
					else {
						objs.emplace_back(light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape, color); // dir=0 (unused)
						objs.back().obj_id = num_light_stacks;
					}
					num_lights_added = 1;
				}
				if (is_lit) {r->lit_by_floor |= (1ULL << (f&63));} // flag this floor as being lit (for up to 64 floors)
			} // end light placement
			if (r->no_geom) continue; // no other geometry for this room
			if (has_stairs && !pri_hall.is_all_zeros()) continue; // no other geometry in office building rooms that have stairs
			float tot_light_amt(light_amt); // unitless, somewhere around 1.0
			if (is_lit) {tot_light_amt += room_light_intensity;} // light surface area divided by room surface area with some fudge constant

			// place a table and maybe some chairs near the center of the room 95% of the time if it's not a hallway
			if (rgen.rand_float() < 0.95) {add_table_and_chairs(rgen, *r, ped_bcubes, room_id, room_center, chair_color, 0.1, tot_light_amt, is_lit);}
			//if (z == bcube.z1()) {} // any special logic that goes on the first floor is here
		} // for f
		num_light_stacks += num_lights_added;
	} // for r
	add_stairs_and_elevators(rgen);
	objs.shrink_to_fit();
	interior->room_geom->light_bcubes.resize(num_light_stacks); // allocate but don't fill un until needed
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	unsigned const num_stairs = 12;
	float const window_vspacing(get_window_vspace()), floor_thickness(FLOOR_THICK_VAL*window_vspacing);
	float const stair_dz(window_vspacing/(num_stairs+1)), stair_height(stair_dz + floor_thickness);
	vector<room_object_t> &objs(interior->room_geom->objs);
	interior->room_geom->stairs_start = objs.size();

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator) continue; // for elevator, not stairs
		bool const dim(i->dim), dir(i->dir);
		// Note: stairs always start at floor_thickness above the landing z1, ignoring landing z2/height
		float const tot_len(i->get_sz_dim(dim)), floor_z(i->z1() + floor_thickness - window_vspacing), step_len_pos(tot_len/num_stairs);
		float step_len((dir ? 1.0 : -1.0)*step_len_pos), z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i);

		if (i->shape == SHAPE_STRAIGHT || i->shape == SHAPE_WALLED) { // straight stairs
			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // Note: room_id=0, not tracked, unused
			}
		}
		else if (i->shape == SHAPE_U) { // U-shaped stairs
			bool const side(0); // for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
			stair.d[!dim][side] = i->get_center_dim(!dim);
			step_len *= 2.0;

			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				if (n == num_stairs/2) { // reverse direction and switch to other side
					step_len *= -1.0;
					stair.d[!dim][ side] = i->d[!dim][side];
					stair.d[!dim][!side] = i->get_center_dim(!dim);
				}
				bool const is_rev(n >= num_stairs/2);
				stair.d[dim][dir^is_rev^1] = pos; stair.d[dim][dir^is_rev] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, side^is_rev); // Note: room_id=0, not tracked, unused
				objs.back().shape = SHAPE_STAIRS_U;
			} // for n
		}
		else {assert(0);} // unknown stairs shape

		if (i->shape == SHAPE_U || i->shape == SHAPE_WALLED) { // add walls around stairs for this floor
			float const wall_hw(0.15*step_len_pos), half_thick(0.5*floor_thickness);
			stair = *i;
			stair.z2() -= 0.5*floor_thickness; // prevent z-fighting on top floor
			stair.z1()  = max(bcube.z1()+half_thick, floor_z-half_thick); // full height
			set_wall_width(stair, i->d[dim][dir], wall_hw, dim);
			objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // back/end of stairs
			stair.d[dim][!dir] = i->d[dim][!dir];

			for (unsigned d = 0; d < 2; ++d) { // sides of stairs
				set_wall_width(stair, i->d[!dim][d], wall_hw, !dim);
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir);
			}
		}
	} // for i
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		// add any internal elevator parts with type=TYPE_ELEVATOR here
	}
}

void building_t::gen_and_draw_room_geom(shader_t &s, vect_cube_t &ped_bcubes, unsigned building_ix, int ped_ix, bool shadow_only) {
	if (!interior) return;
	if (is_rotated()) return; // no room geom for rotated buildings

	if (!has_room_geom()) {
		rand_gen_t rgen;
		rgen.set_state(building_ix, parts.size()); // set to something canonical per building
		ped_bcubes.clear();
		if (ped_ix >= 0) {get_ped_bcubes_for_building(ped_ix, building_ix, ped_bcubes);}
		gen_room_details(rgen, ped_bcubes); // generate so that we can draw it
		assert(has_room_geom());
	}
	interior->room_geom->draw(s, shadow_only);
}

void building_t::clear_room_geom() {
	if (!has_room_geom()) return;
	interior->room_geom->clear(); // free VBO data before deleting the room_geom object
	interior->room_geom.reset();
}

bool building_t::place_person(point &ppos, float radius, rand_gen_t &rgen) const {
	if (!interior || interior->rooms.empty()) return 0; // should be error case
	float const window_vspacing(get_window_vspace()), floor_thickness(FLOOR_THICK_VAL*window_vspacing), fc_thick(0.5*floor_thickness);

	for (unsigned n = 0; n < 100; ++n) { // make 100 attempts
		room_t const &room(interior->rooms[rgen.rand() % interior->rooms.size()]); // select a random room
		if (min(room.dx(), room.dy()) < 4.0*radius) continue; // room to small to place a person
		unsigned const num_floors(calc_num_floors(room, window_vspacing, floor_thickness));
		assert(num_floors > 0);
		unsigned const floor_ix(rgen.rand() % num_floors); // place person on a random floor
		// Note: people are placed before lights are assigned to rooms, so this may not work and must be handled during light placement
		if (room.lit_by_floor && !(room.lit_by_floor & (1ULL << (floor_ix&63)))) continue; // don't place person in an unlit room
		point pos;
		pos.z = room.z1() + fc_thick + window_vspacing*floor_ix;
		for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(room.d[d][0]+radius, room.d[d][1]-radius);} // random XY point inside this room
		cube_t bcube(pos);
		bcube.expand_by(radius); // expand more in Z?
		if (!is_valid_stairs_elevator_placement(bcube, radius, radius)) continue;
		bool bad_place(0);

		// Note: people are placed before room geom is generated for all buildings, so this may not work and will have to be handled during room geom placement
		if (interior->room_geom) { // check placement against room geom objects
			vector<room_object_t> const &objs(interior->room_geom->objs);
			auto objs_end(objs.begin() + interior->room_geom->stairs_start); // skip stairs

			for (auto i = objs.begin(); i != objs_end; ++i) {
				if (i->intersects(bcube)) {bad_place = 1; break;}
			}
		}
		if (bad_place) continue;
		ppos = pos;
		return 1;
	} // for n
	return 0;
}

void building_t::get_exclude_cube(point const &pos, cube_t const &skip, cube_t &exclude) const {
	float const cube_pad(4.0*grass_width), extent(bcube.get_max_extent());
	float dmin_sq(extent*extent); // start with a large value, squared

	for (auto p = parts.begin(); p != get_real_parts_end_inc_sec(); ++p) { // find closest part, including garages/sheds
		if (skip.contains_cube_xy(*p)) continue; // already contained, skip
		float const dist_sq(p2p_dist_sq(pos, p->closest_pt(pos)));
		if (dist_sq < dmin_sq) {exclude = *p; dmin_sq = dist_sq;} // keep if closest part to pos
	}
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // try to expand the cube to cover more parts (pri parts only)
		if (skip.contains_cube_xy(*p)) continue; // already contained, skip
		cube_t cand_ge(exclude);
		cand_ge.union_with_cube(*p);
		if (cand_ge.get_area_xy() < 1.05f*(exclude.get_area_xy() + p->get_area_xy())) {exclude = cand_ge;} // union mostly includes the two parts
	}
	if (!exclude.is_all_zeros()) {exclude.expand_by_xy(cube_pad);} // exclude grass blades that partially intersect the building interior
}

void building_t::update_grass_exclude_at_pos(point const &pos, vector3d const &xlate) const {
	get_exclude_cube(pos, cube_t(),       grass_exclude1); // first/closest cube
	get_exclude_cube(pos, grass_exclude1, grass_exclude2); // second closest cube
	if (!grass_exclude1.is_all_zeros()) {grass_exclude1 += xlate;}
	if (!grass_exclude2.is_all_zeros()) {grass_exclude2 += xlate;}
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

bool is_cube_close_to_door(cube_t const &c, float dmin, cube_t const &door) {
	bool const dim(door.dy() < door.dx());
	if (c.d[!dim][0] > door.d[!dim][1] || c.d[!dim][1] < door.d[!dim][0]) return 0; // no overlap in !dim
	float const min_dist((dmin == 0.0f) ? door.get_sz_dim(!dim) : dmin); // if dmin==0, use door width (so that door has space to open)
	return (c.d[dim][0] < door.d[dim][1]+min_dist && c.d[dim][1] > door.d[dim][0]-min_dist); // within min_dist
}
bool building_t::is_cube_close_to_doorway(cube_t const &c, float dmin) const {
	// Note: we want to test this for things like stairs, but exterior doors likely haven't been allocated at this point, so we have to check for that during door placement
	for (auto i = doors.begin(); i != doors.end(); ++i) { // test exterior doors
		cube_t const door(i->get_bcube());
		if (c.z2() < door.z1() || c.z1() > door.z2()) continue;
		if (is_cube_close_to_door(c, dmin, door)) return 1;
	}
	return (interior ? interior->is_cube_close_to_doorway(c, dmin) : 0); // test interior doors
}
bool building_interior_t::is_cube_close_to_doorway(cube_t const &c, float dmin) const { // ignores zvals
	for (auto i = doors.begin(); i != doors.end(); ++i) { // interior doors
		if (is_cube_close_to_door(c, dmin, *i)) return 1;
	}
	return 0;
}
bool building_interior_t::is_blocked_by_stairs_or_elevator(cube_t const &c, float dmin) const {
	cube_t tc(c);
	tc.expand_by_xy(dmin); // no pad in z
	if (has_bcube_int(tc, elevators)) return 1;
	tc.z1() -= 0.001*tc.dz(); // expand slightly to avoid placing an object exactly at the top of the stairs
	return has_bcube_int(tc, stairwells); // must check zval to exclude stairs and elevators in parts with other z-ranges
}
bool building_t::is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, float dmin) const {
	cube_t place_area(room);
	if (dmin != 0.0f) {place_area.expand_by_xy(-dmin);} // shrink by dmin
	if (!place_area.contains_cube_xy(c))   return 0; // not contained in interior part of the room
	if (is_cube_close_to_doorway(c, dmin)) return 0; // too close to a doorway
	if (interior && interior->is_blocked_by_stairs_or_elevator(c, dmin)) return 0; // faster to check only one per stairwell, but then we need to store another vector?
	if (has_bcube_int(c, blockers)) return 0; // Note: ignores dmin
	return 1;
}

void building_interior_t::finalize() {
	remove_excess_cap(floors);
	remove_excess_cap(ceilings);
	remove_excess_cap(rooms);
	remove_excess_cap(doors);
	remove_excess_cap(landings);
	remove_excess_cap(stairwells);
	remove_excess_cap(elevators);
	for (unsigned d = 0; d < 2; ++d) {remove_excess_cap(walls[d]);}
}

colorRGBA const WOOD_COLOR(0.9, 0.7, 0.5); // light brown, multiplies wood texture color

// skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2 to match CSG cube flags
void rgeom_mat_t::add_cube_to_verts(cube_t const &c, colorRGBA const &color, unsigned skip_faces) {
	vertex_t v;
	v.set_c4(color);

	// Note: stolen from draw_cube() with tex coord logic, back face culling, etc. removed
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (skip_faces & (1 << (2*(2-n) + j))) continue; // skip this face
			v.set_ortho_norm(n, j);
			v.v[n] = c.d[n][j];

			for (unsigned s1 = 0; s1 < 2; ++s1) {
				v.v[d[1]] = c.d[d[1]][s1];
				v.t[0] = tex.tscale_x*v.v[d[1]];

				for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
					v.v[d[0]] = c.d[d[0]][k^j^s1^1]; // need to orient the vertices differently for each side
					v.t[1] = tex.tscale_y*v.v[d[0]];
					quad_verts.push_back(v);
				}
			}
		} // for j
	} // for i
}

void rgeom_mat_t::add_vcylin_to_verts(cube_t const &c, colorRGBA const &color) {
	float const radius(0.5*min(c.dx(), c.dy())); // should be equal/square
	point const center(c.get_cube_center());
	point const ce[2] = {point(center.x, center.y, c.z1()), point(center.x, center.y, c.z2())};
	unsigned const ndiv(N_CYL_SIDES);
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius, radius, ndiv, v12));
	float const ndiv_inv(1.0/ndiv);
	unsigned qix(quad_verts.size()), tix(tri_verts.size());
	quad_verts.resize(qix + 4*ndiv);
	color_wrapper cw(color);

	for (unsigned i = 0; i < ndiv; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			unsigned const S(i + j), s(S%ndiv);
			float const ts(1.0f - S*ndiv_inv);
			vector3d const normal(vpn.n[s] + vpn.n[(S+ndiv-1)%ndiv]); // normalize?
			quad_verts[qix++].assign(vpn.p[(s<<1)+ j], normal, ts, 1.0*( j), cw.c);
			quad_verts[qix++].assign(vpn.p[(s<<1)+!j], normal, ts, 1.0*(!j), cw.c);
		}
	} // for i
	assert(qix == quad_verts.size());
	// add bottom end cap using triangles, currently using all TCs=0.0
	tri_verts.resize(tix + 3*ndiv);

	for (unsigned i = 0; i < ndiv; ++i) {
		for (unsigned j = 0; j < 2; ++j) {tri_verts[tix++].assign(vpn.p[((i + j)%ndiv)<<1], -plus_z, 0.0, 0.0, cw.c);}
		tri_verts[tix++].assign(ce[0], -plus_z, 0.0, 0.0, cw.c); // center
	}
	assert(tix == tri_verts.size());
}

void rgeom_mat_t::create_vbo() {
	num_tverts = tri_verts.size();
	num_qverts = quad_verts.size();
	vector_add_to(tri_verts, quad_verts);
	vbo.create_and_upload(quad_verts);
	clear_container(tri_verts);  // no longer needed
	clear_container(quad_verts); // no longer needed
}

void rgeom_mat_t::draw(shader_t &s, bool shadow_only) {
	if (shadow_only && tex.emissive) return; // assume this is a light source and shouldn't produce shadows
	assert(vbo.vbo_valid());
	assert(num_tverts > 0 || num_qverts > 0);
	if (!shadow_only) {tex.set_gl(s);} // ignores texture scale for now
	vbo.pre_render();
	vertex_t::set_vbo_arrays();
	if (num_qverts > 0) {draw_quads_as_tris(num_qverts);}
	if (num_tverts > 0) {glDrawArrays(GL_TRIANGLES, num_qverts, num_tverts);}
	tex.unset_gl(s);
}

void get_tc_leg_cubes(cube_t const &c, float width, cube_t cubes[4]) {
	for (unsigned y = 0; y < 2; ++y) {
		for (unsigned x = 0; x < 2; ++x) {
			cube_t leg(c);
			leg.d[0][x] += (1.0f - width)*(x ? -1.0f : 1.0f)*c.dx();
			leg.d[1][y] += (1.0f - width)*(y ? -1.0f : 1.0f)*c.dy();
			cubes[2*y+x] = leg;
		}
	}
}
void building_room_geom_t::add_tc_legs(cube_t const &c, colorRGBA const &color, float width, float tscale) {
	rgeom_mat_t &mat(get_wood_material(tscale));
	cube_t cubes[4];
	get_tc_leg_cubes(c, width, cubes);
	for (unsigned i = 0; i < 4; ++i) {mat.add_cube_to_verts(cubes[i], color, (EF_Z1 | EF_Z2));} // skip top and bottom faces
}

colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c) {
	if (display_mode & 0x10) return c; // disable this when using indir lighting
	return c * (0.5f + 0.5f*min(sqrt(o.light_amt), 1.5f)); // use c.light_amt as an approximation for ambient lighting due to sun/moon
}

void building_room_geom_t::add_table(room_object_t const &c, float tscale) { // 6 quads for top + 4 quads per leg = 22 quads = 88 verts
	cube_t top(c), legs_bcube(c);
	top.z1() += 0.85*c.dz(); // 15% of height
	legs_bcube.z2() = top.z1();
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(top, color); // all faces drawn
	add_tc_legs(legs_bcube, color, 0.08, tscale);
}

void building_room_geom_t::add_chair(room_object_t const &c, float tscale) { // 6 quads for seat + 5 quads for back + 4 quads per leg = 27 quads = 108 verts
	float const height(c.dz());
	cube_t seat(c), back(c), legs_bcube(c);
	seat.z1() += 0.32*height;
	seat.z2()  = back.z1() = seat.z1() + 0.07*height;
	legs_bcube.z2() = seat.z1();
	back.d[c.dim][c.dir] += 0.88f*(c.dir ? -1.0f : 1.0f)*c.get_sz_dim(c.dim);
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.2*tscale)).add_cube_to_verts(seat, apply_light_color(c, c.color)); // all faces drawn
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(back, color, EF_Z1); // skip bottom face
	add_tc_legs(legs_bcube, color, 0.15, tscale);
}

void building_room_geom_t::add_stair(room_object_t const &c, float tscale) {
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.5*tscale)).add_cube_to_verts(c, colorRGBA(0.85, 0.85, 0.85)); // all faces drawn
}

void building_room_geom_t::add_elevator(room_object_t const &c, float tscale) {
	assert(0); // not yet implemented
}

void building_room_geom_t::add_light(room_object_t const &c, float tscale) {
	// Note: need to use a different texture (or -1) for is_on because emissive flag alone does not cause a material change
	bool const is_on(c.is_lit());
	tid_nm_pair_t tp((is_on ? (int)WHITE_TEX : (int)PLASTER_TEX), tscale);
	tp.emissive = is_on;
	rgeom_mat_t &mat(get_material(tp));
	if      (c.shape == SHAPE_CUBE ) {mat.add_cube_to_verts  (c, c.color, EF_Z2);} //untextured, skip top face
	else if (c.shape == SHAPE_CYLIN) {mat.add_vcylin_to_verts(c, c.color);}
	else {assert(0);}
}

void building_room_geom_t::clear() {
	for (auto m = materials.begin(); m != materials.end(); ++m) {m->clear();}
	materials.clear();
	objs.clear();
	light_bcubes.clear();
}

unsigned building_room_geom_t::get_num_verts() const {
	unsigned num_verts(0);
	for (auto m = materials.begin(); m != materials.end(); ++m) {num_verts += m->num_qverts + m->num_tverts;}
	return num_verts;
}

rgeom_mat_t &building_room_geom_t::get_material(tid_nm_pair_t const &tex) {
	// for now we do a simple linear search because there shouldn't be too many unique materials
	for (auto m = materials.begin(); m != materials.end(); ++m) {
		if (m->tex == tex) {return *m;}
	}
	materials.emplace_back(tex); // not found, add a new material
	return materials.back();
}
rgeom_mat_t &building_room_geom_t::get_wood_material(float tscale) {
	return get_material(tid_nm_pair_t(WOOD2_TEX, tscale)); // hard-coded for common material
}

colorRGBA room_object_t::get_color() const {
	switch (type) {
	case TYPE_TABLE: return WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX));
	case TYPE_CHAIR: return (color + WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX)))*0.5; // 50% seat color / 50% wood legs color
	case TYPE_STAIR: return LT_GRAY; // close enough
	case TYPE_ELEVATOR: return LT_GRAY; // ???
	case TYPE_BCASE: return WOOD_COLOR;
	case TYPE_DESK:  return WOOD_COLOR;
	case TYPE_TCAN:  return BLACK;
	}
	return color; // default case - probably should always set color so that we can return it here
}

void building_room_geom_t::create_vbos() {
	if (empty()) return; // no geom
	float const tscale(2.0/obj_scale);

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible()) continue;
		assert(i->is_strictly_normalized());
		switch (i->type) {
		case TYPE_NONE:  assert(0); // not supported
		case TYPE_TABLE: add_table(*i, tscale); break;
		case TYPE_CHAIR: add_chair(*i, tscale); break;
		case TYPE_STAIR: add_stair(*i, tscale); break;
		case TYPE_ELEVATOR: add_elevator(*i, tscale); break;
		case TYPE_LIGHT: add_light(*i, tscale); break; // light fixture
		case TYPE_BOOK:  assert(0); break; // book - WRITE
		case TYPE_BCASE: assert(0); break; // bookcase - WRITE
		case TYPE_DESK:  assert(0); break; // desk - WRITE
		case TYPE_TCAN:  assert(0); break; // trashcan - WRITE
		default: assert(0); // undefined type
		}
	} // for i
	// Note: verts are temporary, but cubes are likely needed for things such as collision detection with the player (if it ever gets implemented)
	for (auto m = materials.begin(); m != materials.end(); ++m) {m->create_vbo();}
}

void building_room_geom_t::draw(shader_t &s, bool shadow_only) { // non-const because it creates the VBO
	if (empty()) return; // no geom
	if (materials.empty()) {create_vbos();} // create materials if needed
	for (auto m = materials.begin(); m != materials.end(); ++m) {m->draw(s, shadow_only);}
	vbo_wrap_t::post_render();
}

