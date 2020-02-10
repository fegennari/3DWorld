// 3D World - Building Interior Generation
// by Frank Gennari 11/15/19

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"

float const FLOOR_THICK_VAL = 0.1; // 10% of floor spacing

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
		for (auto i = parts.begin(); i != get_real_parts_end(); ++i) { // garages and sheds are excluded since they have no doors
			cube_t c(*i + xlate);
			if (!c.contains_pt(pos2)) continue; // not interior to this part
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

				if (check_interior && had_coll && pos2.z - xlate.z > part_z2) { // player standing on top of a building with a sloped roof
					if (point_in_polygon_2d(pos_xlate.x, pos_xlate.y, i->pts, i->npts, 0, 1)) {
						vector3d const normal(i->get_norm());
						float const rdist(dot_product_ptv(normal, pos_xlate, i->pts[0]));
						pos2.z += i->get_norm().z*(radius - rdist);
					}
				}
				else { // normal case for bouncing object, etc.
					vector3d const normal(i->get_norm());
					float const rdist(dot_product_ptv(normal, pos_xlate, i->pts[0]));

					if (fabs(rdist) < radius && sphere_poly_intersect(i->pts, i->npts, pos_xlate, normal, rdist, radius)) {
						pos2 += normal*(radius - rdist); // update current pos
						had_coll = 1; // flag as colliding
						if (cnorm_ptr) {*cnorm_ptr = ((normal.z < 0.0) ? -1.0 : 1.0)*normal;} // make sure normal points up
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
			bool const dim(c->dx() < c->dy());
			max_eq(pos[dim], (c->d[dim][0] + radius)); // force the sphere onto the stairs
			min_eq(pos[dim], (c->d[dim][1] - radius));
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
	float door_shift(bcube.dz()); // start with a large value
	if (invert_normal) {swap(door.pts[0], door.pts[1]); swap(door.pts[2], door.pts[3]);} // swap vertex order to invert normal

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find the part that this door was added to
		float const dist(door.pts[0][dim] - p->d[dim][dir]); // signed
		if (fabs(dist) < fabs(door_shift)) {door_shift = dist;}
	}
	assert(fabs(door_shift) < bcube.dz());
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
		if (do_split) {split_in_xy(base, rgen);} // generate L, T, or U shape
		else { // single part, entire cube/cylinder
			parts.push_back(base);
			if ((rgen.rand()&3) != 0) {gen_sloped_roof(rgen);} // 75% chance
			gen_details(rgen);
		}
		end_add_parts();
		gen_interior(rgen, 0);
		gen_building_doors_if_needed(rgen);
		return; // for now the bounding cube
	}
	// generate building levels and splits
	float const height(base.dz()), dz(height/num_levels);
	assert(height > 0.0);

	if (!do_split && (rgen.rand()&3) < (was_cube ? 2 : 3)) { // oddly shaped multi-sided overlapping sections (50% chance for cube buildings and 75% chance for others)
		point const llc(base.get_llc()), sz(base.get_size());
		float const abs_min_edge_move(0.5*mat.get_floor_spacing()); // same as door width
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
				assert(bc.is_strictly_normalized());
				bool contained(0);
				for (auto p = parts.begin(); p != parts.end(); ++p) {contained |= p->contains_cube(bc);}
				if (!contained) {valid = 1; break;} // success
			} // for n
			if (!valid) break; // remove this part and end the building here
			if (i == 0 || !is_simple_cube()) {parts.push_back(bc); continue;} // no splitting
			split_cubes_recur(bc, parts, 0, parts.size()); // split this cube against all previously added cubes and remove overlapping areas
		} // for i
		parts.shrink_to_fit(); // optional
		std::reverse(parts.begin(), parts.end()); // highest part should be last so that it gets the roof details
		calc_bcube_from_parts(); // update bcube
		gen_details(rgen);
		end_add_parts();
		gen_interior(rgen, 1);
		gen_building_doors_if_needed(rgen);
		return;
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
	}
	else {
		if ((rgen.rand()&3) != 0) {gen_sloped_roof(rgen);} // 67% chance
		if (num_levels <= 3) {gen_details(rgen);}
	}
	end_add_parts();
	gen_interior(rgen, 0);
	gen_building_doors_if_needed(rgen);
}

void building_t::split_in_xy(cube_t const &seed_cube, rand_gen_t &rgen) {

	// generate L, T, U, H, + shape
	point const llc(seed_cube.get_llc()), sz(seed_cube.get_size());
	int const shape(rand()%9); // 0-8
	bool const is_hp(shape >= 7);
	bool const dim(rgen.rand_bool()); // {x,y}
	bool const dir(is_hp ? 1 : rgen.rand_bool()); // {neg,pos} - H/+ shapes are always pos
	float const div(is_hp ? rgen.rand_uniform(0.2, 0.4) : rgen.rand_uniform(0.3, 0.7)), s1(rgen.rand_uniform(0.2, 0.4)), s2(rgen.rand_uniform(0.6, 0.8)); // split pos in 0-1 range
	float const dpos(llc[dim] + div*sz[dim]), spos1(llc[!dim] + s1*sz[!dim]), spos2(llc[!dim] + s2*sz[!dim]); // split pos in cube space
	unsigned const start(parts.size()), num((shape >= 6) ? 3 : 2);
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
		float const dpos2(llc[dim] + (1.0 - div)*sz[dim]); // other end
		parts[start+1].d[ dim][ dir] = dpos2;
		parts[start+1].d[!dim][ 0  ] = spos1; // middle part
		parts[start+1].d[!dim][ 1  ] = spos2;
		parts[start+2].d[ dim][!dir] = dpos2; // full width part
		break;
	}
	case 8: { // +
		float const dpos2(llc[dim] + (1.0 - div)*sz[dim]); // other end
		parts[start+0].d[!dim][ 0  ] = spos1;
		parts[start+0].d[!dim][ 1  ] = spos2;
		parts[start+2].d[!dim][ 0  ] = spos1;
		parts[start+2].d[!dim][ 1  ] = spos2;
		parts[start+1].d[ dim][ dir] = dpos2; // middle part
		parts[start+2].d[ dim][!dir] = dpos2; // partial width part
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
	bool const centered(door_center_shift == 0.0 || hallway_dim == (unsigned char)dim); // center doors connected to primary hallways
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
				has_garage = 1;
				c.d[!dim][dir2 ]  = base.d[!dim][dir2]; // shove it into the opposite corner of the bcube
				c.d[ dim][dir  ]  = base.d[ dim][dir ]; // shove it into the opposite corner of the bcube
				c.d[!dim][!dir2] -= dist1; // move away from bcube edge
				c.d[ dim][!dir ] -= dist2; // move away from bcube edge
				c.z2() = c.z1() + min(min(c.dx(), c.dy()), height); // no taller than x or y size; Note: z1 same as part1
			}
			parts.push_back(c); // support column or shed/garage
		} // end house details
		calc_bcube_from_parts(); // maybe calculate a tighter bounding cube
	} // end type != 0  (multi-part house)
	else if (gen_door) { // single cube house
		door_dir  = rgen.rand_bool(); // select a random dir
		door_part = 0; // only one part
	}
	gen_interior(rgen, 0); // before adding door

	if (gen_door) {
		cube_t door(place_door(parts[door_part], door_dim, door_dir, door_height, door_center, door_pos, 0.25, 0.5, 0, rgen));

		if (interior && interior->is_blocked_by_stairs_or_elevator(door, 0.5*door_height)) { // blocked by stairs - switch door to other side/dim
			door_dim ^= 1;
			door_dir  = (two_parts ? ((door_dim == dim) ? dir : dir2) : (door_dir^1)); // if we have a porch/shed/garage, put the door on that side
			if (door_dim == dim && two_parts && detail_type == 0) {door_dir ^= 1;} // put it on the opposite side so that the second part isn't in the way

			if (door_center != 0.0) { // reclaculate for L-shaped house
				door_center = door_cube.get_center_dim(!door_dim) + 0.5f*((door_dim == dim) ? dist1 : dist2);
				door_pos    = door_cube.d[door_dim][!door_dir];
				door_part   = ((door_dim == dim) ? 0 : 1); // which part the door is connected to
			}
			door = place_door(parts[door_part], door_dim, door_dir, door_height, door_center, door_pos, 0.25, 0.5, 0, rgen); // keep door_height
		}
		add_door(door, door_part, door_dim, door_dir, 0);
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
	add_roof_to_bcube();
	gen_grayscale_detail_color(rgen, 0.4, 0.8); // for roof
}

tquad_with_ix_t set_door_from_cube(cube_t const &c, bool dim, bool dir, unsigned type, float pos_adj, bool opened, bool swap_sides) {

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
		door.pts[1][ dim] = door.pts[2][ dim] = pos + width*(dir ? 1.0 : -1.0);
	}
	return door;
}

bool building_t::add_door(cube_t const &c, unsigned part_ix, bool dim, bool dir, bool for_building) {

	if (c.is_all_zeros()) return 0;
	vector3d const sz(c.get_size());
	assert(sz[dim] == 0.0 && sz[!dim] > 0.0 && sz.z > 0.0);
	unsigned const type(for_building ? (unsigned)tquad_with_ix_t::TYPE_BDOOR : (unsigned)tquad_with_ix_t::TYPE_HDOOR);
	doors.push_back(set_door_from_cube(c, dim, dir, type, 0.01*sz[!dim], 0, 0)); // opened=0, swap_sides=0
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
}

void building_t::gen_details(rand_gen_t &rgen) { // for the roof

	unsigned const num_blocks(roof_tquads.empty() ? (rgen.rand() % 9) : 0); // 0-8; 0 if there are roof quads (houses, etc.)
	bool const add_walls(is_simple_cube() && roof_tquads.empty()); // simple cube buildings with flat roofs
	has_antenna = (rgen.rand() & 1);
	details.resize(num_blocks + 4*add_walls + has_antenna);
	assert(!parts.empty());
	if (details.empty()) return; // nothing to do
	cube_t const &top(parts.back()); // top/last part

	if (num_blocks > 0) {
		float const xy_sz(top.get_size().xy_mag());
		vector<point> points; // reused across calls

		for (unsigned i = 0; i < num_blocks; ++i) {
			cube_t &c(details[i]);
			float const height_scale(0.0035f*(top.dz() + bcube.dz())); // based on avg height of current section and entire building
			float const height(height_scale*rgen.rand_uniform(1.0, 4.0));

			while (1) {
				c.set_from_point(point(rgen.rand_uniform(top.x1(), top.x2()), rgen.rand_uniform(top.y1(), top.y2()), 0.0));
				c.expand_by(vector3d(xy_sz*rgen.rand_uniform(0.01, 0.08), xy_sz*rgen.rand_uniform(0.01, 0.06), 0.0));
				if (!top.contains_cube_xy(c)) continue; // not contained
				if (is_simple_cube()) break; // success/done
				bool contained(1);

				for (unsigned j = 0; j < 4; ++j) { // check cylinder/ellipse
					point const pt(c.d[0][j&1], c.d[1][j>>1], 0.0); // XY only
					if (!check_part_contains_pt_xy(top, pt, points)) {contained = 0; break;}
				}
				if (contained) break; // success/done
			} // end while
			c.z1() = top.z2(); // z1
			c.z2() = top.z2() + height; // z2
		} // for i
	}
	if (add_walls) { // add walls around the roof; we don't need to draw all sides of all cubes, but it's probably not worth the trouble to sort it out
		float const window_vspacing(get_material().get_floor_spacing());
		float const height(0.4*window_vspacing), width(0.2*window_vspacing), z1(top.z2()), z2(top.z2()+height);
		details[num_blocks+0] = cube_t(top.x1(), top.x2(), top.y1(), top.y1()+width, z1, z2);
		details[num_blocks+1] = cube_t(top.x1(), top.x2(), top.y2()-width, top.y2(), z1, z2);
		details[num_blocks+2] = cube_t(top.x1(), top.x1()+width, top.y1()+width, top.y2()-width, z1, z2);
		details[num_blocks+3] = cube_t(top.x2()-width, top.x2(), top.y1()+width, top.y2()-width, z1, z2);
	}
	if (has_antenna) { // add antenna
		float const radius(0.003f*rgen.rand_uniform(1.0, 2.0)*(top.dx() + top.dy()));
		float const height(rgen.rand_uniform(0.25, 0.5)*top.dz());
		cube_t &antenna(details.back());
		antenna.set_from_point(top.get_cube_center());
		antenna.expand_by(vector3d(radius, radius, 0.0));
		antenna.z1() = top.z2(); // z1
		antenna.z2() = bcube.z2() + height; // z2 (use bcube to include sloped roof)
	}
	for (auto i = details.begin(); i != details.end(); ++i) {max_eq(bcube.z2(), i->z2());} // extend bcube z2 to contain details
	if (roof_tquads.empty()) {gen_grayscale_detail_color(rgen, 0.2, 0.6);} // for antenna and roof
}

void building_t::gen_sloped_roof(rand_gen_t &rgen) { // Note: currently not supported for rotated buildings

	assert(!parts.empty());
	if (!is_simple_cube()) return; // only simple cubes are handled
	cube_t const &top(parts.back()); // top/last part
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
	add_roof_to_bcube();
	//max_eq(bcube.z2(), z2);
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

void remove_section_from_cube(cube_t &c, cube_t &c2, float v1, float v2, bool xy) { // c is input+output cube, c2 is other output cube
	//if (!(v1 > c.d[xy][0] && v1 < v2 && v2 < c.d[xy][1])) {cout << TXT(v1) << TXT(v2) << TXT(c.d[xy][0]) << TXT(c.d[xy][1]) << TXT(xy) << endl;}
	assert(v1 > c.d[xy][0] && v1 < v2 && v2 < c.d[xy][1]); // v1/v2 must be interior values for cube
	c2 = c; // clone first cube
	c.d[xy][1] = v1; c2.d[xy][0] = v2; // c=low side, c2=high side
}
void remove_section_from_cube_and_add_door(cube_t &c, cube_t &c2, float v1, float v2, bool xy, vect_cube_t &doors) {
	remove_section_from_cube(c, c2, v1, v2, xy);
	cube_t door(c);
	door.d[!xy][0] = door.d[!xy][1] = c.get_center_dim(!xy); // zero area at wall centerline
	door.d[ xy][0] = v1; door.d[ xy][1] = v2;
	doors.push_back(door);
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
	window_border *= 0.95; // adjust based on window frame so that wall doesn't end right at the edge
	float const uv(fract((val - c.d[dim][0])/window_spacing));
	return (uv > window_border && uv < 1.0f-window_border);
}

struct split_cube_t : public cube_t {
	float door_lo[2][2], door_hi[2][2]; // per {dim x dir}
	
	split_cube_t(cube_t const &c) : cube_t(c) {
		door_lo[0][0] = door_lo[0][1] = door_lo[1][0] = door_lo[1][1] = door_hi[0][0] = door_hi[0][1] = door_hi[1][0] = door_hi[1][1] = 0.0f;
	}
	bool bad_pos(float val, bool dim) const {
		for (unsigned d = 0; d < 2; ++d) { // check both dirs (wall end points)
			if (door_lo[dim][d] < door_hi[dim][d] && val > door_lo[dim][d] && val < door_hi[dim][d]) return 1;
		}
		return 0;
	}
};

unsigned calc_num_floors(cube_t const &c, float window_vspacing, float floor_thickness) {
	float const z_span(c.dz() - floor_thickness);
	assert(z_span > 0.0);
	unsigned const num_floors(round_fp(z_span/window_vspacing)); // round - no partial floors
	assert(num_floors <= 100); // sanity check
	return num_floors;
}

void subtract_cube_xy(cube_t const &c, cube_t const &r, cube_t *out) { // subtract r from c; ignores zvals
	assert(c.contains_cube_xy(r));
	for (unsigned i = 0; i < 4; ++i) {out[i] = c;}
	out[0].y2() = r.y1(); // bottom
	out[1].y1() = r.y2(); // top
	out[2].y1() = r.y1(); out[2].y2() = r.y2(); out[2].x2() = r.x1(); // left center
	out[3].y1() = r.y1(); out[3].y2() = r.y2(); out[3].x1() = r.x2(); // right center
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
	float const floor_thickness(FLOOR_THICK_VAL*window_vspacing), fc_thick(0.5*floor_thickness);
	float const doorway_width(0.5*window_vspacing), doorway_hwidth(0.5*doorway_width);
	float const wall_thick(0.5*floor_thickness), wall_half_thick(0.5*wall_thick), wall_edge_spacing(0.05*wall_thick), min_wall_len(4.0*doorway_width);
	float const wwf(global_building_params.get_window_width_fract()), window_border(0.5*(1.0 - wwf)); // (0.0, 1.0)
	// houses have at most two parts; exclude garage, shed, porch, porch support, etc.
	unsigned const num_parts(get_real_num_parts());
	vector<split_cube_t> to_split;
	uint64_t must_split[2] = {0,0};
	unsigned first_wall_to_split[2] = {0,0};
	// allocate space for all floors
	unsigned tot_num_floors(0), tot_num_stairwells(0), tot_num_landings(0); // num floor/ceiling cubes, not number of stories; used only for reserving vectors

	for (auto p = parts.begin(); p != (parts.begin() + num_parts); ++p) {
		bool const has_stairs(!is_house || p == parts.begin()); // assumes one set of stairs or elevator per part
		unsigned const cubes_per_floor(has_stairs ? 4 : 1); // account for stairwell cutouts
		unsigned const num_floors(calc_num_floors(*p, window_vspacing, floor_thickness));
		tot_num_floors     += cubes_per_floor*(num_floors - 1) + 1; // first floor has no cutout
		tot_num_stairwells += (has_stairs && num_floors > 1);
		tot_num_landings   += (has_stairs ? (num_floors - 1) : 0);
	}
	if (has_garage) {++tot_num_floors;}
	interior->ceilings.reserve(tot_num_floors);
	interior->floors  .reserve(tot_num_floors);
	interior->landings.reserve(tot_num_landings);
	interior->stairwells.reserve(tot_num_stairwells);
	
	// generate walls and floors for each part;
	// this will need to be modified to handle buildings that have overlapping parts, or skip those building types completely
	for (auto p = parts.begin(); p != (parts.begin() + num_parts); ++p) {
		unsigned const num_floors(calc_num_floors(*p, window_vspacing, floor_thickness));
		if (num_floors == 0) continue; // not enough space to add a floor (can this happen?)
		// for now, assume each part has the same XY bounds and can use the same floorplan; this means walls can span all floors and don't need to be duplicated for each floor
		vector3d const psz(p->get_size());
		bool const min_dim(psz.y < psz.x); // hall dim
		float const cube_width(psz[min_dim]);
		bool const first_part(p == parts.begin());
		bool const use_hallway(!is_house && (first_part || (p-1)->z1() < p->z1()) && (p+1 == parts.end() || (p+1)->z1() > p->z1()) && cube_width > 4.0*min_wall_len);
		unsigned const rooms_start(interior->rooms.size()), part_id(p - parts.begin());
		cube_t hall;

		if (use_hallway) {
			// building with rectangular slice (no adjacent exterior walls at this level), generate rows of offices
			int const num_windows   (get_num_windows_on_side(p->d[!min_dim][0], p->d[!min_dim][1]));
			int const num_windows_od(get_num_windows_on_side(p->d[ min_dim][0], p->d[ min_dim][1])); // other dim, for use in hallway width calculation
			int const windows_per_room((num_windows > 5) ? 2 : 1); // 1-2 windows per room
			int const num_rooms((num_windows+windows_per_room-1)/windows_per_room); // round up
			bool const partial_room((num_windows % windows_per_room) != 0); // an odd number of windows leaves a small room at the end
			assert(num_rooms >= 0 && num_rooms < 1000); // sanity check
			float const window_hspacing(psz[!min_dim]/num_windows), room_len(window_hspacing*windows_per_room);
			float const num_hall_windows((num_windows_od & 1) ? 1.4 : 1.8); // hall either contains 1 (odd) or 2 (even) windows, wider for single window case to make room for stairs
			float const hall_width(num_hall_windows*psz[min_dim]/num_windows_od);
			float const room_width(0.5f*(cube_width - hall_width)); // rooms are the same size on each side of the hallway
			float const hwall_extend(0.5f*(room_len - doorway_width - wall_thick));
			float const hall_wall_pos[2] = {(p->d[min_dim][0] + room_width), (p->d[min_dim][1] - room_width)};
			hallway_dim = !min_dim; // cache in building for later use
			vect_cube_t &room_walls(interior->walls[!min_dim]), &hall_walls(interior->walls[min_dim]);
			room_walls.reserve(2*(num_rooms-1));
			hall_walls.reserve(2*(num_rooms+1));
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
				rwall.translate_dim(room_len, !min_dim);
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
			if (num_windows_od >= 7) { // at least 7 windows (3 on each side of hallway)
				// TODO_INT: add secondary halls to each side of this one and create two more rows of rooms
			}
			// add rooms and doors
			interior->rooms.reserve(2*num_rooms + 1); // two rows of rooms + hallway
			interior->doors.reserve(2*num_rooms);
			float pos(p->d[!min_dim][0]);

			for (int i = 0; i < num_rooms; ++i) {
				float const wall_end(p->d[!min_dim][1]);
				bool const is_small_room((pos + room_len) > wall_end);
				float const next_pos(is_small_room ? wall_end : (pos + room_len)); // clamp to end of building to last row handle partial room)

				for (unsigned d = 0; d < 2; ++d) { // lo, hi
					cube_t c(*p); // copy zvals and exterior wall pos
					c.d[ min_dim][!d] = hall_wall_pos[d];
					c.d[!min_dim][ 0] = pos;
					c.d[!min_dim][ 1] = next_pos;
					add_room(c, part_id);
					interior->rooms.back().is_office = 1;
					cube_t door(c); // copy zvals and wall pos
					door.d[ min_dim][d] = hall_wall_pos[d]; // set to zero area at hallway
					door.d[!min_dim][0] += hwall_extend; // shrink to doorway width
					door.d[!min_dim][1] -= hwall_extend;
					interior->doors.push_back(door);
				}
				pos = next_pos;
			} // for i
			hall = *p;
			for (unsigned e = 0; e < 2; ++e) {hall.d[min_dim][e] = hall_wall_pos[e];}
			add_room(hall, part_id); // add hallway as room
			if (p->z1() == bcube.z1() || pri_hall.is_all_zeros()) {pri_hall = hall;} // assign to primary hallway if on first floor of hasn't yet been assigned
			interior->rooms.back().is_hallway = interior->rooms.back().no_geom = 1;
			for (unsigned d = 0; d < 2; ++d) {first_wall_to_split[d] = interior->walls[d].size();} // don't split any walls added up to this point
		}
		else { // generate random walls using recursive 2D slices
			bool const no_walls(min(p->dx(), p->dy()) < min_wall_len); // not enough space to add a room (chimney, porch support, garage, shed, etc.)
			assert(to_split.empty());
			if (no_walls) {add_room(*p, part_id);} // add entire part as a room
			else {to_split.emplace_back(*p);} // seed room is entire part, no door
			float window_hspacing[2] = {0.0};
			
			if (first_part) { // reserve walls/rooms/doors - take a guess at the correct size
				for (unsigned d = 0; d < 2; ++d) {interior->walls[d].reserve(8*parts.size());}
				interior->rooms.reserve(8*parts.size()); // two rows of rooms + optional hallway
				interior->doors.reserve(4*parts.size());
			}
			for (unsigned d = 0; d < 2; ++d) {
				int const num_windows(get_num_windows_on_side(p->d[d][0], p->d[d][1]));
				window_hspacing[d] = psz[d]/num_windows;
			}
			while (!to_split.empty()) {
				split_cube_t c(to_split.back()); // Note: non-const because door_lo/door_hi is modified during T-junction insert
				to_split.pop_back();
				vector3d const csz(c.get_size());
				bool wall_dim(0); // which dim the room is split by
				if      (csz.y > min_wall_len && csz.x > 1.25*csz.y) {wall_dim = 0;} // split long room in x
				else if (csz.x > min_wall_len && csz.y > 1.25*csz.x) {wall_dim = 1;} // split long room in y
				else {wall_dim = rgen.rand_bool();} // choose a random split dim for nearly square rooms
				
				if (min(csz.x, csz.y) < min_wall_len) {
					add_room(c, part_id);
					continue; // not enough space to add a wall
				}
				float wall_pos(0.0);
				bool const on_edge(c.d[wall_dim][0] == p->d[wall_dim][0] || c.d[wall_dim][1] == p->d[wall_dim][1]); // at edge of the building - make sure walls don't intersect windows
				bool pos_valid(0);
				
				for (unsigned num = 0; num < 20; ++num) { // 20 tries to choose a wall pos that's not inside a window
					wall_pos = cube_rand_side_pos(c, wall_dim, 0.25, (1.5*doorway_width + wall_thick), rgen);
					if (on_edge && is_val_inside_window(*p, wall_dim, wall_pos, window_hspacing[wall_dim], window_border)) continue; // try a new wall_pos
					if (c.bad_pos(wall_pos, wall_dim)) continue; // intersects doorway from prev wall, try a new wall_pos
					pos_valid = 1; break; // done, keep wall_pos
				}
				if (!pos_valid) { // no valid pos, skip this split
					add_room(c, part_id);
					continue;
				}
				cube_t wall(c), wall2; // copy from cube; shared zvals, but X/Y will be overwritten per wall
				create_wall(wall, wall_dim, wall_pos, fc_thick, wall_half_thick, wall_edge_spacing);
				float const doorway_pos(cube_rand_side_pos(c, !wall_dim, 0.25, doorway_width, rgen));
				float const lo_pos(doorway_pos - doorway_hwidth), hi_pos(doorway_pos + doorway_hwidth);
				remove_section_from_cube_and_add_door(wall, wall2, lo_pos, hi_pos, !wall_dim, interior->doors);
				interior->walls[wall_dim].push_back(wall);
				interior->walls[wall_dim].push_back(wall2);
				bool const do_split(csz[wall_dim] > max(global_building_params.wall_split_thresh, 1.0f)*min_wall_len); // split into two smaller rooms

				for (unsigned d = 0; d < 2; ++d) { // still have space to split in other dim, add the two parts to the stack
					split_cube_t c_sub(c);
					c_sub.d[wall_dim][d] = wall.d[wall_dim][!d]; // clip to wall pos
					c_sub.door_lo[!wall_dim][d] = lo_pos - wall_half_thick; // set new door pos in this dim (keep door pos in other dim, if set)
					c_sub.door_hi[!wall_dim][d] = hi_pos + wall_half_thick;
					if (do_split) {to_split.push_back(c_sub);} else {add_room(c_sub, part_id);} // leaf case (unsplit), add a new room
				}
			} // end while()
			// insert walls to split up parts into rectangular rooms
			float const min_split_wall_len(0.75*min_wall_len); // allow a shorter than normal wall because these walls have higher priority
			bool const too_small(min(p->dx(), p->dy()) < min_split_wall_len);

			for (auto p2 = parts.begin(); p2 != get_real_parts_end() && !too_small; ++p2) {
				if (p2 == p) continue; // skip self
				if (min(p2->dx(), p2->dy()) < min_split_wall_len) continue; // too small, skip

				for (unsigned dim = 0; dim < 2; ++dim) {
					for (unsigned dir = 0; dir < 2; ++dir) {
						float const val(p->d[!dim][dir]);
						if (p2->d[!dim][!dir] != val) continue; // not adjacent
						if (p2->z1() >= p->z2() || p2->z2() <= p->z1()) continue; // no overlap in Z
						if (p2->d[dim][0] > p->d[dim][0] || p2->d[dim][1] < p->d[dim][1]) continue; // not contained in dim (don't have to worry about Z-shaped case)

						if (p2->d[dim][0] == p->d[dim][0] && p2->d[dim][1] == p->d[dim][1]) { // same xy values, must only vary in z
							if (p2->z2() < p->z2()) continue; // add wall only on one side (arbitrary)
						}
						cube_t wall;
						wall.z1() = max(p->z1(), p2->z1()) + fc_thick; // shared Z range
						wall.z2() = min(p->z2(), p2->z2()) - fc_thick;
						wall.d[ dim][0] = p->d[dim][0] + wall_edge_spacing; // shorter part side with slight offset
						wall.d[ dim][1] = p->d[dim][1] - wall_edge_spacing;
						if (wall.get_sz_dim(dim) < min_split_wall_len) continue; // wall is too short to add (can this happen?)
						wall.d[!dim][ dir] = val;
						wall.d[!dim][!dir] = val + (dir ? -1.0 : 1.0)*wall_thick;
						must_split[!dim] |= (1ULL << (interior->walls[!dim].size() & 63)); // flag this wall for extra splitting
						interior->walls[!dim].push_back(wall);
					} // for dir
				} // for dim
			} // for p2
		} // end wall placement
		add_ceilings_floors_stairs(rgen, *p, hall, num_floors, rooms_start, use_hallway, first_part);
	} // for p (parts)

	if (has_garage) { // add garage/shed floor and ceiling
		assert(num_parts < parts.size());
		cube_t const &garage(parts[num_parts]);
		cube_t C(garage);
		C.z2() = C.z1() + fc_thick;
		interior->floors.push_back(C);
		C.z2() = garage.z2();
		C.z1() = C.z2() - fc_thick;
		interior->ceilings.push_back(C);
	}
	// attempt to cut extra doorways into long walls if there's space to produce a more connected floorplan
	for (unsigned d = 0; d < 2; ++d) { // x,y: dim in which the wall partitions the room (wall runs in dim !d)
		vect_cube_t &walls(interior->walls[d]);
		vect_cube_t const &perp_walls(interior->walls[!d]);

		// Note: iteration will include newly added all segments to recursively split long walls
		for (unsigned w = first_wall_to_split[d]; w < walls.size(); ++w) { // skip hallway walls
			bool pref_split(must_split[d] & (1ULL << (w & 63)));

			for (unsigned nsplits = 0; nsplits < 4; ++nsplits) { // at most 4 splits
				cube_t &wall(walls[w]); // take a reference here because a prev iteration push_back() may have invalidated it
				float const len(wall.get_sz_dim(!d)), min_split_len((pref_split ? 0.75 : 1.5)*min_wall_len);
				if (len < min_split_len) break; // not long enough to split - done
				// walls currently don't run along the inside of exterior building walls, so we don't need to handle that case yet
				bool was_split(0);

				for (unsigned ntries = 0; ntries < (pref_split ? 10U : 4U); ++ntries) { // choose random doorway positions and check against perp_walls for occlusion
					float const doorway_pos(cube_rand_side_pos(wall, !d, 0.2, 1.5*doorway_width, rgen));
					float const lo_pos(doorway_pos - doorway_hwidth), hi_pos(doorway_pos + doorway_hwidth);
					bool valid(1);

					for (auto p = (perp_walls.begin() + first_wall_to_split[!d]); p != perp_walls.end(); ++p) { // skip hallway walls
						if (p->z1() != wall.z1()) continue; // not the same zval/story
						if (p->d[!d][1] < lo_pos-wall_thick || p->d[!d][0] > hi_pos+wall_thick) continue; // no overlap with wall
						if (p->d[ d][1] > wall.d[d][0]-wall_thick && p->d[d][0] < wall.d[d][1]+wall_thick) {valid = 0; break;} // has perp intersection
					}
					if (valid && !pref_split) { // don't split walls into small segments that border the same two rooms on both sides (two doorways between the same pair of rooms)
						float const lo[2] = {wall.d[!d][0]-wall_thick, lo_pos}, hi[2] = {hi_pos, wall.d[!d][1]+wall_thick}; // ranges of the two split wall segments, grown a bit

						for (unsigned s = 0; s < 2; ++s) { // check both wall segments
							bool contained[2] = {0,0};

							for (unsigned e = 0; e < 2; ++e) { // check both directions from the wall
								for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
									if (wall.d[d][e] < r->d[d][0] || wall.d[d][e] > r->d[d][1]) continue; // wall not inside room in dim d/dir e
									if (lo[s] > r->d[!d][0] && hi[s] < r->d[!d][1]) {contained[e] = 1; break;} // entire wall contained in span of room
								}
							}
							if (contained[0] && contained[1]) {valid = 0; break;} // wall seg contained in rooms on both sides => two doors in same wall between rooms => drop
						} // for s
					}
					if (!valid) continue;
					cube_t cand(wall);
					cand.d[!d][0] = lo_pos; cand.d[!d][1] = hi_pos;
					if (interior->is_blocked_by_stairs_or_elevator(wall, doorway_width)) continue; // stairs in the way, skip; should we assert !pref_split?
					cube_t wall2;
					remove_section_from_cube_and_add_door(wall, wall2, lo_pos, hi_pos, !d, interior->doors);
					walls.push_back(wall2); // Note: invalidates wall reference
					was_split = 1;
					break;
				} // for ntries
				if (!was_split) break; // no more splits
				pref_split = 0; // already split, no longer preferred
			} // for nsplits
		} // for w
	} // for d

	// add stairs to connect together stacked parts for office buildings; must be done last after all walls/ceilings/floors have been assigned
	for (auto p = parts.begin(); p != (parts.begin() + num_parts); ++p) {connect_stacked_parts_with_stairs(rgen, *p);}
	interior->finalize();
}

void building_t::add_ceilings_floors_stairs(rand_gen_t &rgen, cube_t const &part, cube_t const &hall, unsigned num_floors, unsigned rooms_start, bool use_hallway, bool first_part) {

	float const window_vspacing(get_material().get_floor_spacing());
	float const floor_thickness(FLOOR_THICK_VAL*window_vspacing), fc_thick(0.5*floor_thickness), doorway_width(0.5*window_vspacing);
	float const ewidth(1.5*doorway_width); // for elevators
	float z(part.z1());
	cube_t stairs_cut, elevator_cut;
	bool stairs_dim(0), add_elevator(0);
	unsigned stairs_shape(0); // straight by default

	// add stairwells and elevator shafts
	if (num_floors == 1) {} // no need for stairs or elevator
	else if (use_hallway) { // part is the hallway cube
		add_elevator = 1; //rgen.rand_bool();
		if (interior->landings.empty()) {interior->landings.reserve(add_elevator ? 1 : (num_floors-1));}
		assert(!interior->rooms.empty());
		room_t &room(interior->rooms.back()); // hallway is always the last room to be added
		bool const long_dim(hall.dx() < hall.dy());
		cube_t stairs(hall); // start as hallway

		if (add_elevator) {
			point center(room.get_cube_center());
			float const center_shift(0.125*room.get_sz_dim(long_dim)*(rgen.rand_bool() ? -1.0 : 1.0));
			center[long_dim] += center_shift; // make elevator off-center
			elevator_t elevator(room, long_dim, rgen.rand_bool(), rgen.rand_bool()); // elevator shaft
			elevator.x1() = center.x - 0.5*ewidth; elevator.x2() = center.x + 0.5*ewidth;
			elevator.y1() = center.y - 0.5*ewidth; elevator.y2() = center.y + 0.5*ewidth;
			interior->elevators.push_back(elevator);
			room.has_elevator = 1;
			elevator_cut      = elevator;
			stairs.translate_dim(-center_shift, long_dim); // shift stairs in the opposite direction
		}
		// always add stairs
		for (unsigned dim = 0; dim < 2; ++dim) { // shrink in XY
			bool const is_step_dim(bool(dim) == long_dim); // same orientation as the hallway
			float shrink(stairs.get_sz_dim(dim) - (is_step_dim ? 4.0*doorway_width : 0.9*ewidth)); // set max size of stairs opening, slightly narrower than elevator
			stairs.d[dim][0] += 0.5*shrink; stairs.d[dim][1] -= 0.5*shrink; // centered in the hallway
		}
		room.has_stairs = 1;
		stairs_cut      = stairs;
	}
	else if (!is_house || interior->stairwells.empty()) { // only add stairs to first part of a house unless we haven't added stairs yet
		// add elevator half of the time to building parts, but not the first part (to guarantee we have at least one set of stairs)
		// it might not be possible to place an elevator a part with no interior rooms, but that should be okay, because some other part will still have stairs
		// do we need support for multiple floor cutouts stairs + elevator in this case as well?
		add_elevator = (!is_house && !first_part && rgen.rand_bool());
		unsigned const rooms_end(interior->rooms.size()), num_avail_rooms(rooms_end - rooms_start);
		assert(num_avail_rooms > 0); // must have added at least one room
		float stairs_scale(1.0);

		for (unsigned N = 0; N < 4; ++N) {
			unsigned const rand_ix(rgen.rand()); // choose a random starting room to make into a stairwell

			for (unsigned n = 0; n < num_avail_rooms; ++n) { // try all available rooms starting with the selected one to see if we can fit a stairwell/elevator in any of them
				unsigned const stairs_room(rooms_start + (rand_ix + n)%num_avail_rooms);
				room_t &room(interior->rooms[stairs_room]);

				if (add_elevator) {
					if (min(room.dx(), room.dy()) < 2.0*ewidth) continue; // room is too small to place an elevator
					bool placed(0);

					for (unsigned y = 0; y < 2 && !placed; ++y) { // try all 4 corners
						for (unsigned x = 0; x < 2 && !placed; ++x) {
							// don't place elevators on building exteriors blocking windows or between parts where they would block doorways
							if (room.d[0][x] == part.d[0][x] || room.d[1][y] == part.d[1][y]) continue;
							bool const dim(rgen.rand_bool()), is_open(rgen.rand_bool());
							elevator_t elevator(room, dim, !(dim ? y : x), is_open); // elevator shaft
							elevator.d[0][!x] = elevator.d[0][x] + (x ? -ewidth : ewidth);
							elevator.d[1][!y] = elevator.d[1][y] + (y ? -ewidth : ewidth);
							elevator.expand_by_xy(-0.01*ewidth); // shrink to leave a small gap between the outer wall to prevent z-fighting
							if (is_cube_close_to_doorway(elevator)) continue; // try again
							interior->elevators.push_back(elevator);
							elevator_cut = elevator;
							placed       = 1; // successfully placed
						} // for x
					} // for y
					if (!placed) continue; // try another room
					room.has_elevator = 1;
					room.no_geom      = 1;
				}
				else { // stairs
					cube_t cutout(room);
					cutout.expand_by_xy(-floor_thickness); // padding around walls
					float const dx(cutout.dx()), dy(cutout.dy()); // choose longer dim if high aspect ratio
					float const stairs_sz(stairs_scale*doorway_width);
					if      (dx > 1.2*dy) {stairs_dim = 0;}
					else if (dy > 1.2*dx) {stairs_dim = 1;}
					else {stairs_dim = rgen.rand_bool();} // close to square
					if (cutout.get_sz_dim(stairs_dim) < 4.0*stairs_sz || cutout.get_sz_dim(!stairs_dim) < 3.0*stairs_sz) continue; // not enough space for stairs

					for (unsigned dim = 0; dim < 2; ++dim) { // shrink in XY
						bool const is_step_dim(bool(dim) == stairs_dim);
						float shrink(cutout.get_sz_dim(dim) - (is_step_dim ? 4.0 : 1.2)*stairs_sz); // set max size of stairs opening
						max_eq(shrink, 2.0f*stairs_sz); // allow space for doors to open and player to enter/exit
						cutout.d[dim][0] += 0.5*shrink; cutout.d[dim][1] -= 0.5*shrink; // centered in the room

						if (!is_step_dim) { // see if we can push the stairs to the wall on one of the sides without blocking a doorway
							bool const first_dir(rgen.rand_bool());

							for (unsigned d = 0; d < 2; ++d) {
								bool const dir(bool(d) ^ first_dir);
								// if the room is on the edge of the part that's not on the building bcube exterior, then this room connects two parts and we need to place a door here later;
								// technically, it's more accurate to check for an adjacent part, but that's somewhat difficult to do here
								if (room.d[dim][dir] == part.d[dim][dir] && part.d[dim][dir] != bcube.d[dim][dir]) continue;
								cube_t cand(cutout);
								float const shift(0.95f*(cand.d[dim][dir] - room.d[dim][dir])); // negative if dir==1, add small gap to prevent z-fighting and FP accuracy asserts
								cand.d[dim][0] -= shift; cand.d[dim][1] -= shift; // close the gap - flush with the wall
								if (!is_cube_close_to_doorway(cand)) {cutout = cand; break;} // keep if it's good
							} // for d
						}
					} // for dim
					if (interior->landings.empty()) {interior->landings.reserve(num_floors-1);}
					assert(cutout.is_strictly_normalized());
					stairs_cut      = cutout;
					room.has_stairs = 1;
					//room.no_geom    = 1;
				}
				break; // success - done
			} // for n
			if (!first_part || add_elevator || !stairs_cut.is_all_zeros()) break; // successfully placed stairs, or not required to place stairs
			stairs_scale -= 0.1; // shrink stairs a bit and try again
		} // for N
	} // end stairs/elevator placement

	// add ceilings and floors; we have num_floors+1 separators; the first is only a floor, and the last is only a ceiling
	cube_t C(part);
	C.z1() = z; C.z2() = z + fc_thick;
	unsigned const floors_start(interior->floors.size());
	interior->floors.push_back(C); // ground floor, full area
	z += window_vspacing; // move to next floor
	bool const has_stairs(!stairs_cut.is_all_zeros()), has_elevator(!elevator_cut.is_all_zeros());
	cube_t &first_cut(has_elevator ? elevator_cut : stairs_cut); // elevator is larger

	for (unsigned f = 1; f < num_floors; ++f, z += window_vspacing) { // skip first floor - draw pairs of floors and ceilings
		cube_t to_add[8];
		float const zc(z - fc_thick), zf(z + fc_thick);

		if (!has_stairs && !has_elevator) {to_add[0] = part;} // neither - add single cube
		else {
			subtract_cube_xy(part, first_cut, to_add);

			if (has_stairs && has_elevator) { // both
				bool found(0);

				for (unsigned n = 0; n < 4; ++n) { // find the cube where the stairs are placed
					if (!to_add[n].intersects_xy_no_adj(stairs_cut)) continue;
					subtract_cube_xy(to_add[n], stairs_cut, to_add+4); // append up to 4 more cubes
					to_add[n].set_to_zeros(); // this cube was replaced
					found = 1;
					break; // assume there is only one; it's up to the placer step to ensure this
				}
				assert(found);
			}
			if (has_stairs) { // add landings and stairwells
				landing_t landing(stairs_cut, 0, stairs_dim, 0, stairs_shape); // dir is unused and has been set to 0
				landing.z1() = zc; landing.z2() = zf;
				interior->landings.push_back(landing);
				if (f == 1) {interior->stairwells.push_back(stairs_cut);} // only add for first floor
			}
			if (has_elevator) {
				assert(!interior->elevators.empty());
				landing_t landing(elevator_cut, 1, interior->elevators.back().dim, interior->elevators.back().dir);
				landing.z1() = zc; landing.z2() = zf;
				interior->landings.push_back(landing);
			}
		}
		for (unsigned i = 0; i < 8; ++i) { // skip zero area cubes from stairs/elevator shafts along an exterior wall
			cube_t &c(to_add[i]);
			if (c.is_zero_area()) continue;
			c.z1() = zc; c.z2() = z;  interior->ceilings.push_back(c);
			c.z1() = z;  c.z2() = zf; interior->floors  .push_back(c);
			//c.z1() = zf; c.z2() = zc + window_vspacing; // add per-floor walls, door cutouts, etc. here
		}
	} // for f
	C.z1() = z - fc_thick; C.z2() = z;
	interior->ceilings.push_back(C); // roof ceiling, full area
	std::reverse(interior->floors.begin()+floors_start, interior->floors.end()); // order floors top to bottom to reduce overdraw when viewed from above
}

bool building_t::is_valid_stairs_elevator_placement(cube_t const &c, float door_pad, float stairs_pad) const {
	if (is_cube_close_to_doorway(c, stairs_pad)) return 0; // bad
	if (interior->is_blocked_by_stairs_or_elevator(c, door_pad)) return 0; // bad

	// check if any previously placed walls intersect this cand stairs/elevator; we really only need to check the walls from <part> and *p though
	for (unsigned d = 0; d < 2; ++d) {
		if (has_bcube_int(c, interior->walls[d])) return 0;
	}
	return 1;
}

void subtract_cube_from_cube(cube_t const &c, cube_t const &s, vect_cube_t &out) {
	cube_t C;
	if (c.y1() < s.y1()) {C = c; C.y2() = s.y1(); out.push_back(C);} // bottom
	if (c.y2() > s.y2()) {C = c; C.y1() = s.y2(); out.push_back(C);} // top
	if (c.x1() < s.x1()) {C = c; max_eq(C.y1(), s.y1()); min_eq(C.y2(), s.y2()); C.x2() = s.x1(); out.push_back(C);} // left center
	if (c.x2() > s.x2()) {C = c; max_eq(C.y1(), s.y1()); min_eq(C.y2(), s.y2()); C.x1() = s.x2(); out.push_back(C);} // right center
}
void subtract_cube_from_cube_inplace(cube_t const &s, vect_cube_t &cubes, unsigned ix) { // Note: ix is an index to cubes
	unsigned const prev_sz(cubes.size());
	assert(ix < prev_sz);
	cube_t const c(cubes[ix]); // deep copy - reference will become invalid
	subtract_cube_from_cube(c, s, cubes);
	assert(cubes.size() > prev_sz); // must have added at least one cube
	cubes[ix] = cubes.back(); cubes.pop_back(); // reuse this slot for one of the output cubes
}
template<typename T> void subtract_cubes_from_cube(cube_t const &c, T const &sub, vect_cube_t &out, vect_cube_t &out2) { // XY only
	out.clear();
	out.push_back(c);

	for (auto s = sub.begin(); s != sub.end(); ++s) {
		if (s->z1() <= c.z1() || s->z1() >= c.z2() || s->z2() <= c.z2()) continue; // not correct floor
		if (!c.intersects_xy(c)) continue; // no overlap with orig cube (optimization)
		out2.clear();

		// clip all of out against *s, write results to out2, then swap with out
		for (auto i = out.begin(); i != out.end(); ++i) {
			if (!i->intersects_xy(*s)) {out2.push_back(*i); continue;} // no overlap, keep entire cube
			subtract_cube_from_cube(*i, *s, out2);
		}
		out.swap(out2);
	} // for s
}
void subtract_cube_from_cubes(cube_t const &s, vect_cube_t &cubes, vect_cube_t *holes=nullptr) {
	unsigned const num_cubes(cubes.size()); // capture size before splitting

	for (unsigned i = 0; i < num_cubes; ++i) {
		cube_t const &c(cubes[i]);
		if (!c.intersects(s)) continue; // keep it
		if (holes) {holes->push_back(c); holes->back().intersect_with_cube(s);}
		subtract_cube_from_cube_inplace(s, cubes, i); // Note: invalidates c reference
	}
}
template<typename T> void subtract_cubes_from_cubes(T const &sub, vect_cube_t &cubes) {
	for (auto i = sub.begin(); i != sub.end(); ++i) {subtract_cube_from_cubes(*i, cubes);}
}

void subtract_cube_from_floor_ceil(cube_t const &c, vect_cube_t &fs) {
	unsigned const fsz(fs.size()); // capture orig size

	for (unsigned i = 0; i < fsz; ++i) {
		cube_t const &cur(fs[i]);
		if (cur.z1() > c.z2() || cur.z2() < c.z1()) continue; // no z overlap
		if (cur.intersects(c)) {subtract_cube_from_cube_inplace(c, fs, i);} // Note: invalidates cur reference
	}
}

void building_t::connect_stacked_parts_with_stairs(rand_gen_t &rgen, cube_t const &part) {

	float const window_vspacing(get_material().get_floor_spacing()), fc_thick(0.5*FLOOR_THICK_VAL*window_vspacing);
	float const doorway_width(0.5*window_vspacing), stairs_len(4.0*doorway_width);

	if (part.z2() < bcube.z2()) { // if this is the top floor, there is nothing above it
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (*p == part) continue; // skip self
			if (p->z1() != part.z2()) continue; // *p not on top of part
			if (!part.intersects_xy(*p)) continue; // no XY overlap
			cube_t shared(part);
			shared.intersect_with_cube(*p); // dz() == 0
			cube_t pref_shared(shared);
			// Note: parts are sorted top to bottom, so any part above <part> should be before it in parts - but we don't want to rely on that here;
			// however, this does mean that the part above this one has already been processed
			float stairs_width(1.2*doorway_width); // relatively small
			float stairs_pad(1.0*doorway_width), len_with_pad(stairs_len + 2.0*stairs_pad); // pad both ends of stairs to make sure player has space to enter/exit
			if (max(shared.dx(), shared.dy()) < 1.5*len_with_pad || min(shared.dx(), shared.dy()) < 2.0*stairs_width) continue; // too small to add stairs between these parts

			if (!pri_hall.is_all_zeros() && part.contains_cube(pri_hall) && pri_hall.intersects_xy(shared)) { // have a primary hallway in this part
				pref_shared.intersect_with_cube(pri_hall);
				if (max(pref_shared.dx(), pref_shared.dy()) < 1.5*len_with_pad || min(pref_shared.dx(), pref_shared.dy()) < 2.0*stairs_width) {pref_shared = shared;} // too small
			}
			// place stairs in shared area if there's space and no walls are in the way for either the room or above
			cube_t cand;
			cand.z1() = part.z2() - window_vspacing + fc_thick; // top of top floor for this part
			cand.z2() = part.z2() + fc_thick; // top of bottom floor of upper part *p
			bool stairs_added(0);

			// is it better to extend the existing stairs in *p, or the stairs we're creating here (stairs_cut) if they line up?
			for (unsigned n = 0; n < 100; ++n) { // make 100 tries to add stairs
				cube_t place_region((n < 10) ? pref_shared : shared); // use preferred shared area from primary hallway for first 10 iterations

				if (n > 0 && n%10 == 0) { // decrease stairs size slightly every 10 iterations
					stairs_width -= 0.025*doorway_width;
					stairs_pad   -= 0.040*doorway_width;
					len_with_pad -= 0.230*doorway_width;
				}
				bool dim(rgen.rand_bool());
				if (place_region.get_sz_dim(dim) < 1.5*len_with_pad) {dim ^= 1;} // too narrow in this dim, try other dim

				for (unsigned d = 0; d < 2; ++d) {
					float const stairs_sz((bool(d) == dim) ? len_with_pad : stairs_width);
					cand.d[d][0] = rgen.rand_uniform(place_region.d[d][0], (place_region.d[d][1] - stairs_sz)); // LLC
					cand.d[d][1] = cand.d[d][0] + stairs_sz; // URC
				}
				cube_t cand_test(cand);
				cand_test.z1() += 0.5*window_vspacing; cand_test.z2() += 0.5*window_vspacing; // move up a bit so that it intersects exactly the floor below and the floor above
				if (!is_valid_stairs_elevator_placement(cand_test, doorway_width, stairs_pad)) continue; // bad placement
				cand.d[dim][0] += stairs_pad; cand.d[dim][1] -= stairs_pad; // subtract off padding
				unsigned const stairs_shape(0); // straight, for now
				landing_t landing(cand, 0, dim, 0, stairs_shape); // dir is unused and set to 0
				landing.z1() = part.z2() - fc_thick; // only include the ceiling of this part and the floor of *p
				interior->landings.push_back(landing);
				interior->stairwells.push_back(cand);
				// attempt to cut holes in ceiling of this part and floor of above part
				subtract_cube_from_floor_ceil(cand, interior->floors);
				subtract_cube_from_floor_ceil(cand, interior->ceilings);
				break; // success
			} // for n
		} // for p
	}

	// now attempt to extend elevators into floors above/below
	vect_cube_t holes;

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find parts above/below this elevator
			if (!p->contains_cube_xy(*e)) continue; // elevator doesn't extend inside this cube
			bool const is_above(p->z1() == e->z2()), is_below(p->z2() == e->z1());
			if (!is_above && !is_below) continue;
			cube_t extension(*e);
			extension.z1() = p->z1(); extension.z2() = p->z2();
			cube_t cand_test(extension);
			cand_test.z1() += fc_thick; cand_test.z2() -= fc_thick; // shrink slightly in Z so that we don't intersect the original elevator *e
			cand_test.d[e->dim][e->dir] += doorway_width*(e->dir ? 1.0 : -1.0); // add extra space in front of the elevator
			if (!is_valid_stairs_elevator_placement(cand_test, doorway_width, doorway_width)) continue; // bad placement
			min_eq(e->z1(), extension.z1()); max_eq(e->z2(), extension.z2()); // perform extension in Z
			extension.z1() -= fc_thick; extension.z2() += fc_thick; // also cut a hole in the lower ceiling/upper floor
			holes.clear();
			subtract_cube_from_cubes(extension, interior->ceilings);
			subtract_cube_from_cubes(extension, interior->floors, &holes); // capture holes from floors
			
			for (auto h = holes.begin(); h != holes.end(); ++h) {
				landing_t landing(*h, 1, e->dim, e->dir);
				landing.z1() -= fc_thick; // since we only captured floor cutouts, extend them downward to include the ceiling below
				interior->landings.push_back(landing);
			}
		} // for p
	} // for e
}

bool building_t::clip_part_ceiling_for_stairs(cube_t const &c, vect_cube_t &out, vect_cube_t &temp) const { // and elevators
	if (!interior || interior->stairwells.empty()) return 0;
	subtract_cubes_from_cube(c, interior->stairwells, out, temp);

	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) { // handle elevators that span multiple parts
		if (e->z2() <= c.z1()) continue; // elevator shaft doesn't extend this high
		if (e->z2() == c.z2()) continue; // don't cut a hole in the roof of the building where the top of the elevator shaft ends
		subtract_cube_from_cubes(*e, out);
	}
	return 1;
}

void building_t::add_room(cube_t const &room, unsigned part_id) {
	assert(interior);
	room_t r(room, part_id);
	cube_t const &part(parts[part_id]);
	for (unsigned d = 0; d < 4; ++d) {r.ext_sides |= (unsigned(room.d[d>>1][d&1] == part.d[d>>1][d&1]) << d);} // find exterior sides
	interior->rooms.push_back(r);
}

bool building_t::add_table_and_chairs(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers,
	unsigned room_id, point const &place_pos, float rand_place_off, float tot_light_amt, bool is_lit)
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

bool has_bcube_int_exp(cube_t const &bcube, vect_cube_t const &bcubes, float expand) {
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
	unsigned tot_num_rooms(0);
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
		float const light_size((r->is_hallway ? 2.0 : (r->is_office ? 1.5 : 1.0))*floor_thickness); // use larger light for offices and hallways
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
		// place objects on each floor for this room
		for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
			room_center.z = z + fc_thick; // floor height
			bool const top_of_stairs(r->has_stairs && f+1 == num_floors);
			bool is_lit(0), light_dim(room_dim);
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
				unsigned char flags(RO_FLAG_NOCOLL); // no collision detection with lights
				if (is_lit)        {flags |= RO_FLAG_LIT;}
				if (top_of_stairs) {flags |= RO_FLAG_TOS;}
				if (r->has_stairs) {flags |= RO_FLAG_RSTAIRS;}
				bool const check_stairs(!is_house && parts.size() > 1 && f+1 == num_floors); // top floor of building that may have stairs connecting to upper stack
				
				if (r->is_hallway) { // place a light on each side of the stairs, and also between stairs and elevator if there are both
					unsigned const num_lights((r->has_elevator && r->has_stairs) ? 3 : 2);
					float const offset(((num_lights == 3) ? 0.3 : 0.2)*r->get_sz_dim(light_dim)); // closer to the ends in the 3 lights case
					cube_t valid_bounds(*r);
					valid_bounds.expand_by_xy(-0.1*window_vspacing); // add some padding

					for (unsigned d = 0; d < num_lights; ++d) {
						float const delta((d == 2) ? 0.0 : (d ? -1.0 : 1.0)*offset); // last light is in the center
						cube_t hall_light(light);
						hall_light.translate_dim(delta, light_dim);

						if (check_stairs) { // keep moving until not blocked by stairs
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
						}
						objs.emplace_back(hall_light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape); // dir=0 (unused)
					} // for d
				}
				else { // normal room
					if (check_stairs && has_bcube_int_exp(light, interior->stairwells, fc_thick)) {is_lit = 0;} // disable if blocked by stairs
					else {objs.emplace_back(light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape);} // dir=0 (unused)
				}
				if (is_lit) {r->lit_by_floor |= (1ULL << (f&63));} // flag this floor as being lit (for up to 64 floors)
			} // end light placement
			if (r->no_geom) continue; // no other geometry for this room
			float tot_light_amt(light_amt); // unitless, somewhere around 1.0
			if (is_lit) {tot_light_amt += room_light_intensity;} // light surface area divided by room surface area with some fudge constant

			// place a table and maybe some chairs near the center of the room 95% of the time if it's not a hallway
			if (rgen.rand_float() < 0.95) {add_table_and_chairs(rgen, *r, ped_bcubes, room_id, room_center, 0.1, tot_light_amt, is_lit);}
			//if (z == bcube.z1()) {} // any special logic that goes on the first floor is here
		} // for f
	} // for r
	add_stairs_and_elevators(rgen);
	objs.shrink_to_fit();
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	unsigned const num_stairs = 12;
	float const window_vspacing(get_window_vspace()), floor_thickness(FLOOR_THICK_VAL*window_vspacing);
	float const stair_dz(window_vspacing/(num_stairs+1)), stair_height(stair_dz + floor_thickness);
	bool const dir(rgen.rand_bool()); // same for every floor, could alternate for stairwells if we were tracking it
	vector<room_object_t> &objs(interior->room_geom->objs);
	interior->room_geom->stairs_start = objs.size();

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator) continue; // for elevator, not stairs
		bool const dim(i->dx() < i->dy()); // longer dim
		float const tot_len(i->get_sz_dim(dim)), floor_z(i->z2() - window_vspacing);
		float step_len((dir ? 1.0 : -1.0)*tot_len/num_stairs);
		float z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i);

		if (i->shape == 0) { // straight stairs
			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // Note: room_id=0, not tracked, unused
			}
		}
		else { // U-shaped stairs
			bool const side(0); // for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
			stair.d[!dim][side] = i->get_center_dim(!dim);
			step_len *= 2.0;

			// TODO: need larger landing, wider but shorter cutout, fix for walking/collision, railings, etc.
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
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // Note: room_id=0, not tracked, unused
			} // for n
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
			for (auto i = interior->room_geom->objs.begin(); i != interior->room_geom->objs.end(); ++i) {
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

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find closest part
		if (skip.contains_cube_xy(*p)) continue; // already contained, skip
		float const dist_sq(p2p_dist_sq(pos, p->closest_pt(pos)));
		if (dist_sq < dmin_sq) {exclude = *p; dmin_sq = dist_sq;} // keep if closest part to pos
	}
	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // try to expand the cube to cover more parts
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

void building_room_geom_t::add_tc_legs(cube_t const &c, colorRGBA const &color, float width, float tscale) {
	rgeom_mat_t &mat(get_wood_material(tscale));

	for (unsigned y = 0; y < 2; ++y) {
		for (unsigned x = 0; x < 2; ++x) {
			cube_t leg(c);
			leg.d[0][x] += (1.0f - width)*(x ? -1.0f : 1.0f)*c.dx();
			leg.d[1][y] += (1.0f - width)*(y ? -1.0f : 1.0f)*c.dy();
			mat.add_cube_to_verts(leg, color, (EF_Z1 | EF_Z2)); // skip top and bottom faces
		}
	}
}

colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c) {
	// TODO_INT: use c.light_amt as an approximation for ambient lighting due to sun/moon? Do we need per-object material colors?
	return c * (0.5f + 0.5f*min(sqrt(o.light_amt), 1.5f));
	//return c;
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
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.2*tscale)).add_cube_to_verts(seat, apply_light_color(c, colorRGBA(0.2, 0.2, 1.0))); // light blue; all faces drawn
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
	if      (c.shape == SHAPE_CUBE ) {mat.add_cube_to_verts  (c, WHITE, EF_Z2);} // white, untextured, skip top face
	else if (c.shape == SHAPE_CYLIN) {mat.add_vcylin_to_verts(c, WHITE);}
	else {assert(0);}
}

void building_room_geom_t::clear() {
	for (auto m = materials.begin(); m != materials.end(); ++m) {m->clear();}
	materials.clear();
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

colorRGBA const &room_object_t::get_color() const {
	switch (type) {
	case TYPE_NONE:  assert(0); // not supported
	case TYPE_TABLE: return WOOD_COLOR;
	case TYPE_CHAIR: return WOOD_COLOR;
	case TYPE_STAIR: return LT_GRAY; // close enough
	case TYPE_ELEVATOR: return LT_GRAY; // ???
	case TYPE_LIGHT: return WHITE;
	case TYPE_BOOK:  return BLUE; // ???
	case TYPE_BCASE: return WOOD_COLOR;
	case TYPE_DESK:  return WOOD_COLOR;
	case TYPE_TCAN:  return BLACK;
	default: assert(0); // undefined type
	}
	return BLACK; // never gets here
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

