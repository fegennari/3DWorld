// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 5/1/02

#include "3DWorld.h"
#include "physics_objects.h"
#include "mesh.h"
#include "player_state.h"
#include "csg.h"
#include "openal_wrap.h"

set<unsigned> moving_cobjs;

extern bool scene_smap_vbo_invalid;
extern int num_groups;
extern float base_gravity, tstep, temperature;
extern coll_obj_group coll_objects;
extern player_state *sstates;
extern platform_cont platforms;
extern obj_type object_types[NUM_TOT_OBJS];
extern obj_group obj_groups[NUM_TOT_OBJS];


bool push_cobj(unsigned index, vector3d &delta, set<unsigned> &seen);


int cube_polygon_intersect(coll_obj const &c, coll_obj const &p) {

	for (int i = 0; i < p.npoints; ++i) { // check points (fast)
		if (c.contains_pt(p.points[i])) return 1; // definite intersection
	}
	for (int i = 0; i < p.npoints; ++i) { // check edges
		if (check_line_clip(p.points[i], p.points[(i+1)%p.npoints], c.d)) return 1; // definite intersection
	}
	static vector<point> clipped_pts;
	clip_polygon_to_cube(c, p.points, p.npoints, p, clipped_pts); // clip the polygon to the cube
	if (!clipped_pts.empty()) return 1;

	if (p.thickness > MIN_POLY_THICK) { // test extruded (3D) polygon
		point pts[2][4];
		gen_poly_planes(p.points, p.npoints, p.norm, p.thickness, pts);
				
		for (unsigned j = 0; j < 2; ++j) {
			for (int i = 0; i < p.npoints; ++i) {
				if (check_line_clip(pts[j][i], pts[j][(i+1)%p.npoints], c.d)) return 1; // definite intersection
			}
		}
		// need to handle cube completely insde of a thick polygon
		if (sphere_ext_poly_intersect(p.points, p.npoints, p.norm, c.get_cube_center(), 0.0, p.thickness, MIN_POLY_THICK)) return 1;
		static vector<tquad_t> side_pts;
		thick_poly_to_sides(p.points, p.npoints, p.norm, p.thickness, side_pts);

		for (auto i = side_pts.begin(); i != side_pts.end(); ++i) { // clip each face to the cube
			clip_polygon_to_cube(c, i->pts, i->npts, p, clipped_pts);
			if (!clipped_pts.empty()) return 1;
		}
		return 0;
	}
	return 0;
}

int cylin_cube_int_aa_via_circle_rect(coll_obj const &cube, coll_obj const &cylin) {

	bool deq[3];
	UNROLL_3X(deq[i_] = (cylin.points[0][i_] == cylin.points[1][i_]);)

	for (unsigned d = 0; d < 3; ++d) { // approximate projected circle test for x/y/z oriented cylinders
		if (!deq[(d+1)%3] || !deq[(d+2)%3]) continue; // not oriented in direction d
		if ( circle_rect_intersect(cylin.points[0], min(cylin.radius, cylin.radius2), cube, d)) return 1;
		if (!circle_rect_intersect(cylin.points[0], max(cylin.radius, cylin.radius2), cube, d)) return 0; // no intersection
	}
	return 2; // FIXME: finish
}

int cylin_cylin_int(coll_obj const &c1, coll_obj const &c2) {

	float const line_dist(line_line_dist(c2.points[0], c2.points[1], c1.points[0], c1.points[1]));
	if (line_dist > (max(c2.radius, c2.radius2) + max(c1.radius, c1.radius2))) return 0;
	if (c1.line_intersect(c2.points[0], c2.points[1])) return 1;
	if (c2.line_intersect(c1.points[0], c1.points[1])) return 1;
	if (!cylin_cube_int_aa_via_circle_rect(c1, c2))    return 0;
	if (!cylin_cube_int_aa_via_circle_rect(c2, c1))    return 0;
	return 2; // FIXME: finish
}

int poly_cylin_int(coll_obj const &p, coll_obj const &c) {

	if (p.line_intersect(c.points[0], c.points[1])) return 1;

	for (int i = 0; i < p.npoints; ++i) {
		if (c.line_intersect(p.points[i], p.points[(i+1)%p.npoints])) return 1; // definite intersection
	}
	if (!cube_polygon_intersect(c, p)) return 0;
	if (!cylin_cube_int_aa_via_circle_rect(p, c)) return 0;

	if (p.thickness > MIN_POLY_THICK) { // test extruded (3D) polygon
		// WRITE
	}
	return 2; // FIXME: finish
}

int poly_poly_int_test(coll_obj const &p1, coll_obj const &p2) {

	for (int i = 0; i < p1.npoints; ++i) {
		if (p2.line_intersect(p1.points[i], p1.points[(i+1)%p1.npoints])) return 1;
	}
	if (p1.thickness <= MIN_POLY_THICK) return 0;
	// Note: this is inefficient since all edges belong to two faces and are checked twice, but probably doesn't matter
	// maybe we should call gen_poly_planes() instead and check top/bot faces, but this will miss 4 of the edges parallel to the normal
	static vector<tquad_t> side_pts;
	thick_poly_to_sides(p1.points, p1.npoints, p1.norm, p1.thickness, side_pts);

	for (auto i = side_pts.begin(); i != side_pts.end(); ++i) { // intersect the edges from each face
		for (unsigned j = 0; j < i->npts; ++j) {
			if (p2.line_intersect(i->pts[j], i->pts[(j+1)%i->npts])) return 1;
		}
	}
	// need to handle one polygon completely insde of a thick polygon (Note: calls thick_poly_to_sides() again)
	if (sphere_ext_poly_intersect(p1.points, p1.npoints, p1.norm, p2.points[0], 0.0, p1.thickness, MIN_POLY_THICK)) return 1;
	return 0;
}

bool coll_sphere_cylin_int(point const &sc, float sr, coll_obj const &c) {
	if (!sphere_cube_intersect(sc, sr, c)) return 0; // test bcube
	return sphere_intersect_cylinder(sc, sr, c.points[0], c.points[1], c.radius, c.radius2);
}

bool sphere_def_coll_vert_cylin(point const &sc, float sr, point const &cp1, point const &cp2, float cr) {
	return (sc.z >= min(cp1.z, cp2.z) && sc.z <= max(cp1.z, cp2.z) && dist_xy_less_than(sc, cp1, (sr + cr)));
}


// 0: no intersection, 1: intersection, 2: maybe intersection (incomplete)
// 21 total: 15 complete, 5 partial (all cylinder cases), 1 incomplete (capsule-capsule)
// Note: pos toler => adjacency doesn't count; neg toler => adjacency counts
int coll_obj::intersects_cobj(coll_obj const &c, float toler) const {

	if (c.type < type) {return c.intersects_cobj(*this, toler);} // swap arguments
	if (!intersects(c, toler)) return 0; // cube-cube intersection
	float const r1(radius-toler), r2(radius2-toler), cr1(c.radius-toler), cr2(c.radius2-toler);

	// c.type >= type
	switch (type) {
	case COLL_CUBE:
		switch (c.type) {
		case COLL_CUBE:     return 1; // as simple as that
		case COLL_CYLINDER: return circle_rect_intersect(c.points[0], cr1, *this, 2); // in z
		case COLL_SPHERE:   return sphere_cube_intersect(c.points[0], cr1, *this);
		case COLL_CAPSULE:
			if (sphere_cube_intersect(c.points[0], cr1, *this) || sphere_cube_intersect(c.points[1], cr2, *this)) return 1;
			// fallthrough
		case COLL_CYLINDER_ROT: {
			if (check_line_clip(c.points[0], c.points[1], d)) return 1; // definite intersection
			return cylin_cube_int_aa_via_circle_rect(*this, c);
		}
		case COLL_POLYGON: return cube_polygon_intersect(*this, c);
		default: assert(0);
		}

	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		switch (c.type) {
		case COLL_CYLINDER:
			assert(type == COLL_CYLINDER);
			return dist_xy_less_than(points[0], c.points[0], (cr1+radius)); // bcube test guarantees overlap in z
		case COLL_SPHERE:
			if (type == COLL_CYLINDER && sphere_def_coll_vert_cylin(c.points[0], cr1, points[0], points[1], r1)) return 1;
			return coll_sphere_cylin_int(c.points[0], cr1, *this);
		case COLL_CAPSULE:
			if (type == COLL_CYLINDER && (sphere_def_coll_vert_cylin(c.points[0], cr1, points[0], points[1], radius) ||
				                          sphere_def_coll_vert_cylin(c.points[1], cr2, points[0], points[1], radius))) return 1;
			if (coll_sphere_cylin_int(c.points[0], cr1, *this) || coll_sphere_cylin_int(c.points[1], cr2, *this)) return 1;
			// fallthrough
		case COLL_CYLINDER_ROT: return cylin_cylin_int(c, *this);
		case COLL_POLYGON:      return poly_cylin_int (c, *this);
		default: assert(0);
		}

	case COLL_SPHERE:
		switch (c.type) {
		case COLL_SPHERE: return dist_less_than(points[0], c.points[0], (cr1+radius));
		case COLL_CAPSULE: if (dist_less_than(points[0], c.points[0], (cr1+radius)) || dist_less_than(points[0], c.points[1], (cr2+radius))) return 1;
			// fallthrough
		case COLL_CYLINDER_ROT: return coll_sphere_cylin_int(points[0], r1, c);
		case COLL_POLYGON: return sphere_ext_poly_intersect(c.points, c.npoints, c.norm, points[0], r1, c.thickness, MIN_POLY_THICK);
		default: assert(0);
		}

	case COLL_POLYGON: {
		switch (c.type) {
		case COLL_CAPSULE:
			if (sphere_ext_poly_intersect(points, npoints, norm, c.points[0], cr1, thickness, MIN_POLY_THICK)) return 1;
			if (sphere_ext_poly_intersect(points, npoints, norm, c.points[1], cr2, thickness, MIN_POLY_THICK)) return 1;
			return poly_cylin_int(*this, c);
		case COLL_POLYGON: {
			if (thickness   > MIN_POLY_THICK && !cube_polygon_intersect(c, *this)) return 0;
			if (c.thickness > MIN_POLY_THICK && !cube_polygon_intersect(*this, c)) return 0;
			float const poly_toler(max(toler, (thickness + c.thickness)*(1.0f - fabs(dot_product(norm, c.norm)))));

			if (poly_toler > 0.0) { // use toler for edge adjacency tests (for adjacent roof polygons, sponza polygons, etc.)
				for (int i = 0; i < c.npoints; ++i) { // test point adjacency
					for (int j = 0; j < npoints; ++j) {
						if (dist_less_than(points[j], c.points[i], poly_toler)) return 1;
					}
				}
				for (int i = 0; i < c.npoints; ++i) { // test edges
					for (int j = 0; j < npoints; ++j) {
						if (pt_line_seg_dist_less_than(c.points[i],   points[j],   points[(j+1)%  npoints], poly_toler)) return 1;
						if (pt_line_seg_dist_less_than(  points[j], c.points[i], c.points[(i+1)%c.npoints], poly_toler)) return 1;
					}
				}
			}
			return (poly_poly_int_test(c, *this) || poly_poly_int_test(*this, c));
		}
		default: assert(0);
		}

	case COLL_CAPSULE: {
			assert(c.type == COLL_CAPSULE);
			sphere_t const sa[2] = {sphere_t(  points[0],  r1), sphere_t(  points[1],  r2)};
			sphere_t const sb[2] = {sphere_t(c.points[0], cr1), sphere_t(c.points[1], cr2)};

			for (unsigned i = 0; i < 4; ++i) { // 4 sphere-sphere intersections
				if (dist_less_than(sa[i&1].pos, sb[i>>1].pos, (sa[i&1].radius + sb[i>>1].radius))) return 1;
			}
			for (unsigned i = 0; i < 2; ++i) { // 4 sphere-cylinder intersections
				if (coll_sphere_cylin_int(sa[i].pos, sa[i].radius, c    )) return 1;
				if (coll_sphere_cylin_int(sb[i].pos, sb[i].radius, *this)) return 1;
			}
			return cylin_cylin_int(c, *this); // cylinder-cylinder intersection
		}
	}
	default: assert(0);
	}
	return 0;
}


float get_max_cobj_move_delta(coll_obj const &c1, coll_obj const &c2, vector3d const &delta, float step_thresh, float tolerance) {

	float valid_t(0.0);

	// since there is no cobj-cobj closest distance function, binary split the delta range until cobj and c no longer intersect
	for (float t = 0.5, step = 0.25; step > step_thresh; step *= 0.5) { // initial guess is the midpoint
		coll_obj test_cobj(c1); // deep copy
		test_cobj.shift_by(t*delta);
		if (test_cobj.intersects_cobj(c2, tolerance)) {t -= step;} else {valid_t = t; t += step;}
	}
	return valid_t;
}

bool binary_step_moving_cobj_delta(coll_obj const &cobj, vector<unsigned> const &cobjs, vector3d &delta, float tolerance) {

	float step_thresh(0.001);

	for (auto i = cobjs.begin(); i != cobjs.end(); ++i) {
		coll_obj const &c(coll_objects.get_cobj(*i));
		if (cobj.intersects_cobj(c, tolerance)) return 0; // intersects at the starting location, don't allow it to move (stuck)
		float const valid_t(get_max_cobj_move_delta(cobj, c, delta, step_thresh, tolerance));
		if (valid_t == 0.0) return 0; // can't move
		step_thresh /= valid_t; // adjust thresh to avoid tiny steps for large number of cobjs
		delta       *= valid_t;
	} // for i
	return 1;
}

void check_moving_cobj_int_with_dynamic_objs(unsigned index) {

	coll_obj &cobj(coll_objects.get_cobj(index));
	vector<unsigned> cobjs;
	get_intersecting_cobjs_tree(cobj, cobjs, -1, 0.0, 1, 0, -1); // duplicates are okay
	if (cobjs.empty()) return;

	// some dynamic object collided, but we can't tell which one, so iterate over the groups and test them all
	for (int g = 0; g < num_groups; ++g) {
		obj_group &objg(obj_groups[g]);
		if (!objg.enabled || !objg.large_radius()) continue;
		float const robj(object_types[objg.type].radius);
		
		for (unsigned i = 0; i < objg.end_id; ++i) {
			dwobject &obj(objg.get_obj(i));
			if (obj.status != 4) continue; // not stopped
			if (!cobj.sphere_intersects(obj.pos, robj)) continue;
			obj.flags |= WAS_PUSHED;
			obj.flags &= ~ALL_COLL_STOPPED;
			obj.status = 1;
		}
	} // for g
}

bool is_point_supported(coll_obj const &cobj, point const &pos) {

	switch (cobj.type) {
	case COLL_CUBE:     return cobj.contains_pt_xy(pos);
	case COLL_CYLINDER: return dist_xy_less_than(pos, cobj.points[0], cobj.radius);
	case COLL_SPHERE:   return 0; // not flat
	case COLL_CYLINDER_ROT:
		if (cobj.points[0].x != cobj.points[1].x || cobj.points[0].y != cobj.points[1].y) return 0; // non-vertical/not flat
		return dist_xy_less_than(pos, cobj.points[0], ((cobj.points[0].z < cobj.points[1].z) ? cobj.radius2 : cobj.radius)); // use radius at the top
	case COLL_CAPSULE:  return 0; // not flat
	case COLL_POLYGON:
		if (fabs(cobj.norm.z) < 0.9) {return point_in_polygon_2d(pos.x, pos.y, cobj.points, cobj.npoints);}
		
		if (cobj.thickness > MIN_POLY_THICK) { // non-vertical thick polygon
			vector<tquad_t> pts;
			thick_poly_to_sides(cobj.points, cobj.npoints, cobj.norm, cobj.thickness, pts);
			
			for (auto i = pts.begin(); i != pts.end(); ++i) {
				if (fabs(i->get_norm().z) < 0.9) {return point_in_polygon_2d(pos.x, pos.y, i->pts, i->npts);}
			}
		}
		return 0; // not flat
	default: assert(0);
	} // end switch
	return 0; // never gets here
}

bool is_rolling_cobj(coll_obj const &cobj) {

	if (cobj.type == COLL_SPHERE) return 1;
#if 0 // the other cases don't work
	if ((cobj.type == COLL_CYLINDER || cobj.type == COLL_CAPSULE) && cobj.radius == cobj.radius2) { // cylinder/capsule with uniform radius
		vector3d const dir(cobj.points[1] - cobj.points[0]);
		return (dir.z == 0.0 && (dir.x == 0.0 || dir.y == 0.0)); // oriented in either X or Y
	}
#endif
	return 0;
}

float get_cobj_step_height() {return 0.4*C_STEP_HEIGHT*CAMERA_RADIUS;} // cobj can be lifted by 40% of the player step height

void try_drop_movable_cobj(unsigned index) {

	float const tolerance(1.0E-6), cobj_zmin(min(czmin, zbottom));
	float const accel(-0.5*base_gravity*GRAVITY*tstep); // half gravity
	coll_obj &cobj(coll_objects.get_cobj(index));
	float const cobj_height(cobj.d[2][1] - cobj.d[2][0]), prev_v_fall(cobj.v_fall), cur_v_fall(prev_v_fall + accel);
	cobj.v_fall = 0.0; // assume the cobj stops falling; compute the correct v_fall if we reach the end without returning
	float gravity_dz(-tstep*cur_v_fall), max_dz(min(gravity_dz, cobj.d[2][0]-cobj_zmin)); // usually positive
	if (max_dz < tolerance) return; // can't drop further
	float const test_dz(max(max_dz, 0.25f*cobj_height)); // use a min value to ensure the z-slice isn't too small, to prevent instability (for example when riding down on an elevator)
	cube_t bcube(cobj); // start at the current cobj xy
	//bcube.d[2][1]  = cobj.d[2][0]; // top = cobj bottom (Note: more efficient, but doesn't work correctly for elevators)
	bcube.d[2][0] -= test_dz; // bottom (height = test_dz)
	vector<unsigned> cobjs;
	get_intersecting_cobjs_tree(bcube, cobjs, index, tolerance, 0, 0, -1);

	// see if this cobj's bottom edge is colliding with a platform that's moving up (elevator)
	// also, if the cobj is currently intersecting another movable cobj, try to resolve the intersection so that stacking works by moving the cobj up
	for (auto i = cobjs.begin(); i != cobjs.end(); ++i) {
		coll_obj const &c(coll_objects.get_cobj(*i));
		float const dz(c.d[2][1] - cobj.d[2][0]);
		if (dz <= 0 || c.d[2][1] > cobj.d[2][1]) continue; // bottom cobj/platform edge not intersecting

		if (c.cp.flags & COBJ_MOVABLE) { // both cobjs are movable - is this a stack?
			if (c.type == COLL_CUBE || c.type == COLL_CYLINDER || (c.type == COLL_POLYGON && c.norm.x == 0.0 && c.norm.y == 0.0)) {} // flat cobjs can always be stacked
			else if (dz < get_cobj_step_height()) {} // c_top - cobj_bot < step_height
			else if (c.v_fall <= 0.0) continue; // not rising (stopped or falling)
			// assume it's a stack and treat it like a moving platform
		}
		else {
			if (c.platform_id < 0) continue; // not a platform
			if (!platforms.get_cobj_platform(c).is_active()) continue; // platform is not moving (is_moving() is faster but off by one frame on platform stop/change dir)
		}
		if (!cobj.intersects_cobj(c, tolerance)) continue; // no intersection
		cobj.shift_by(vector3d(0.0, 0.0, dz)); // move cobj up
		scene_smap_vbo_invalid = 1;
		return; // or test other cobjs?
	} // for i
	
	// check other cobjs and the mesh to see if this cobj can be dropped
	vector3d delta(0.0, 0.0, -test_dz);
	point const center(cobj.get_center_pt()); // Note: not an accurate center of mass for polygons and truncated cones
	
	if (!binary_step_moving_cobj_delta(cobj, cobjs, delta, tolerance)) { // stuck
		if (prev_v_fall < 4.0*accel) {gen_sound(SOUND_OBJ_FALL, center, 0.5, 1.2);}
		// check for rolling cobjs that can roll downhill
		if (cobjs.size() != 1) return; // can only handle a single supporting cobj
		if (!is_rolling_cobj(cobj)) return; // not rolling
		coll_obj const &c(coll_objects.get_cobj(cobjs.front()));
		if (is_point_supported(c, center)) return; // center of gravity is resting stably, it's stuck
		vector3d move_dir((center - c.get_center_pt()).get_norm()); // FIXME: incorrect for cube/polygon
		
		if (c.type == COLL_CUBE) { // chose +/- x|y for cubes
			if (fabs(move_dir.x) < fabs(move_dir.y)) { // y-edge
				move_dir = ((move_dir.y < 0.0) ? -plus_y : plus_y);
			}
			else { // x-edge
				move_dir = ((move_dir.x < 0.0) ? -plus_x : plus_x);
			}
		}
		delta = 0.05*cobj_height*move_dir; // move 5% of cobj height
		set<unsigned> seen;
		push_cobj(index, delta, seen); // return value is ignored
		return; // done
	}
	float mesh_zval(interpolate_mesh_zval(center.x, center.y, 0.0, 1, 1, 1)); // Note: uses center point, not max mesh height under the cobj; clamped xy

	// check for ice
	int xpos(get_xpos(center.x)), ypos(get_ypos(center.y));
	if (is_in_ice(xpos, ypos)) { // on ice
		float const water_zval(water_matrix[ypos][xpos]);
		if (center.z > water_zval) {mesh_zval = max(mesh_zval, water_zval);} // use water zval if cobj center is above the ice
	}
	float const mesh_dz(mesh_zval - cobj.d[2][0]); // Note: can be positive if cobj is below the mesh
	if (fabs(mesh_dz) < tolerance) {return;} // resting on the mesh (never happens?)

	if (max(delta.z, -max_dz) < mesh_dz) { // under the mesh
		if (prev_v_fall < 10.0*accel) {gen_sound(SOUND_OBJ_FALL, center, 0.2, 0.8);}
		delta.z = mesh_dz; // don't let it go below the mesh
		bool const dim(abs(delta.y) < abs(delta.x));
		float radius(0.5*(cobj.d[dim][1] - cobj.d[dim][0])); // perpendicular to direction of movement
		if      (cobj.type == COLL_CUBE    ) {radius *= 1.4;}
		else if (cobj.type == COLL_SPHERE  ) {radius *= 0.2;} // smaller since bottom surface area is small (maybe also not-vert cylinder?)
		else if (cobj.type == COLL_CYLINDER) {radius  = min(radius, cobj.radius);}
		modify_grass_at(center, radius, 1); // crush grass
	}
	else if (delta.z <= -max_dz) { // cobj falls the entire max distance without colliding, accelerate it
		// set terminal velocity to one cobj_height per timestep to avoid falling completely through another moving cobj (such as an elevator)
		cobj.v_fall = max(cur_v_fall, -cobj_height/tstep);
	} // else |delta.z| may be > |max_dz|, but it was only falling for one frame so should not be noticeable (and should be more stable for small tstep/accel)
	else if (delta.z < 0.0) { // falling the partial distance
		cobj.v_fall = prev_v_fall; // maintain previous falling velocity (no acceleration)
	}
	// handle water
	float depth(0.0);
	point const bot_cent(center.x, center.y, cobj.d[2][0]);

	if (is_underwater(bot_cent, 0, &depth)) {
		if (temperature > W_FREEZE_POINT) { // water
			float const cobj_mass(cobj.get_mass());
			float density(cobj.cp.density);

			if (cobj.type == COLL_CUBE || cobj.type == COLL_CYLINDER) { // has a flat top surface that other cobjs can rest on (what about extruded polygons?)
				cobjs.clear();
				cube_t bcube(cobj); // start at the current cobj xy
				bcube.d[2][1] += 0.1*cobj_height; // z-slice slightly above cobj
				bcube.d[2][0]  = cobj.d[2][0]; // top of cobj
				get_intersecting_cobjs_tree(bcube, cobjs, index, tolerance, 0, 0, -1);
				float tot_mass(cobj_mass);

				// since the z-slice is so thin/small, we simply assume that any cobj intersecting it is resting on the current cobj
				for (auto i = cobjs.begin(); i != cobjs.end(); ++i) {
					coll_obj const &c(coll_objects.get_cobj(*i));

					if (is_underwater(c.get_center_pt(), 0)) { // top cobj is also underwater
						if (c.cp.density > 1.0) {tot_mass += c.volume * (c.cp.density - 1.0);} // if c sinks, add it's net mass when submerged in water
					}
					else {tot_mass += c.get_mass();} // above the water - add the full mass
				}
				density *= tot_mass / cobj_mass; // density is effective total mass / volume
			}
			if (density <= 1.0) { // floats
				if (delta.z <= 0.0) { // falling
					float const target_z(depth - density*cobj_height); // adjust z for constant depth - cobj moves with the water surface
					delta.z     = 0.1*target_z; // smoothly transition from current z to target z
					cobj.v_fall = 0.0;
				}
			}
			else { // sinks
				if (depth > 0.5*cobj_height && mesh_dz < -0.1*cobj_height) { // sinking below the water
					for (unsigned n = 0; n < 4; ++n) {gen_bubble(global_rand_gen.gen_rand_cube_point(cobj));}
				}
				float const fall_rate((density - 1.0)/density); // falling velocity (approaches to 1.0 as density increases)
				max_dz      = min(gravity_dz*fall_rate, max_dz);
				cobj.v_fall = max(cobj.v_fall, -0.001f/tstep); // use a lower terminal velocity in water
			}
			if (prev_v_fall < 8.0*accel) {add_splash(center, xpos, ypos, min(100.0f, -5000*prev_v_fall*cobj_mass), min(2.0f*CAMERA_RADIUS, cobj.get_bsphere_radius()), 1);}
		}
		else if (depth > 0.5*cobj_height) return; // stuck in ice
	}
	if (cobj.v_fall == 0.0 && cobj.type == COLL_CUBE && 0) { // Note: disabled until rotation is handled
		cobj.convert_cube_to_ext_polygon();
		// FIXME: apply rotation if not stable at rest
	}
	delta.z = max(delta.z, -max_dz); // clamp to the real max value if in freefall
	cobj.shift_by(delta); // move cobj down
	scene_smap_vbo_invalid = 1;
	check_moving_cobj_int_with_dynamic_objs(index);
}

bool push_cobj(unsigned index, vector3d &delta, set<unsigned> &seen) {

	coll_obj &cobj(coll_objects.get_cobj(index));
	if (!(cobj.cp.flags & COBJ_MOVABLE)) return 0; // not movable
	delta.z = 0.0; // for now, objects can only be pushed in xy
	float const tolerance(1.0E-6), toler_sq(tolerance*tolerance);
	if (delta.mag_sq() < toler_sq) return 0;

	// make sure the cobj center stays within the scene bounds
	for (unsigned d = 0; d < 2; ++d) { // x/y
		if      (delta[d] > 0.0) {delta[d] = min(delta[d], ( SCENE_SIZE[d]-HALF_DXY - cobj.d[d][1]));} // pos scene edge
		else if (delta[d] < 0.0) {delta[d] = max(delta[d], (-SCENE_SIZE[d] - cobj.d[d][0]));} // neg scene edge
	}
	if (delta.mag_sq() < toler_sq) return 0;

	// determine if this cobj can be moved by checking other static objects
	// dynamic objects should move out of the way of this cobj by themselves
	cube_t bcube(cobj); // orig pos
	bcube += delta; // move to new pos
	bcube.union_with_cube(cobj); // union of original and new pos
	bcube.expand_by(-tolerance); // shrink slightly to avoid false collisions (in prticular with extruded polygons)
	vector<unsigned> cobjs;
	get_intersecting_cobjs_tree(bcube, cobjs, index, tolerance, 0, 0, -1); // duplicates should be okay
	vector3d const start_delta(delta);
	
	if (cobj.type == COLL_CUBE) { // check for horizontally stackable movable cobjs - limited to cubes for now
		seen.insert(index);

		for (unsigned i = 0; i < cobjs.size(); ++i) {
			unsigned const cid(cobjs[i]);
			coll_obj const &c(coll_objects.get_cobj(cid));
			if (!(c.cp.flags & COBJ_MOVABLE) || c.type != COLL_CUBE) continue;
			if (!cobj.intersects_cobj(c, -tolerance)) continue; // no initial intersection/adjacency
			if (seen.find(cid) != seen.end()) continue; // prevent infinite recursion
			seen.insert(cid);
			vector3d delta2(delta);
			if (!push_cobj(cid, delta2, seen)) return 0; // can't push (recursive call)
			delta = delta2; // update with maybe reduced delta
			cobjs[i] = cobjs.back(); cobjs.pop_back(); --i; // remove this cobj
		}
	}
	if (!binary_step_moving_cobj_delta(cobj, cobjs, delta, tolerance)) { // failed to move
		// if there is a ledge (cobj z top) slightly above the bottom of the cobj, maybe we can lift it up;
		// meant to work with ramps and small steps, but not stairs or tree trunks (obviously)
		float step_height(get_cobj_step_height());
		bool has_ledge(0), has_ramp(0), success(0);
		
		for (auto i = cobjs.begin(); i != cobjs.end(); ++i) {
			coll_obj const &c(coll_objects.get_cobj(*i));
			if (c.type == COLL_POLYGON) {has_ramp = 1;} // Note: can only push cobjs up polygon ramps, not cylinders
			else if (c.d[2][1] - cobj.d[2][0] < step_height) {has_ledge = 1; break;} // c_top - cobj_bot < step_height
		}
		if (!has_ledge) {
			if (!has_ramp) return 0; // stuck
			step_height = min(step_height, delta.mag()); // if there is no ledge, limit step height to delta length (45 degree inclination)
		}
		for (unsigned n = 0; n < 10; ++n) { // try raising the cobj in 10 steps
			delta = start_delta; // reset to try again
			cobj.shift_by(vector3d(0.0, 0.0, 0.1*step_height)); // shift up slightly to see if it can be moved now
			if (binary_step_moving_cobj_delta(cobj, cobjs, delta, tolerance)) {success = 1; break;}
		}
		if (!success) {
			cobj.shift_by(vector3d(0.0, 0.0, -step_height)); // shift down to original position
			return 0; // can't move
		}
	}
	// cobj moved, see if it gets teleported
	vector3d cobj_delta(delta);
	float const radius(cobj.get_bsphere_radius());
	point const center(cobj.get_center_pt());
	point cobj_pos(center + cobj_delta); // post-move pos
	
	if (maybe_teleport_object(cobj_pos, 0.5*radius, NO_SOURCE)) { // was teleported
		cobj_delta = cobj_pos - center + delta.get_norm()*radius; // move radius past the teleport end position so that the player doesn't get stuck on the cobj when going through
	}
	// move the cobj and update various state
	cobj.move_cobj(cobj_delta, 1); // move the cobj instead of the player and re-add to coll structure
	moving_cobjs.insert(index); // may already be there
	scene_smap_vbo_invalid = 1;
	check_moving_cobj_int_with_dynamic_objs(index);
	if ((rand2()%1000) == 0) {gen_sound(SOUND_SLIDING, center, 0.1, 1.0);}
	return 1; // moved
}

bool proc_movable_cobj(point const &orig_pos, point &player_pos, unsigned index, int type) {

	if (type == CAMERA && sstates != nullptr && sstates[CAMERA_ID].jump_time > 0) return 0; // can't push while jumping (what about smileys?)
	vector3d delta(orig_pos - player_pos);
	set<unsigned> seen;
	if (!push_cobj(index, delta, seen)) return 0;
	player_pos += delta; // restore player pos, at least partially
	return 1; // moved
}

void proc_moving_cobjs() {

	vector<pair<float, unsigned>> by_z1;

	for (auto i = moving_cobjs.begin(); i != moving_cobjs.end();) {
		if (coll_objects.get_cobj(*i).status != COLL_STATIC) {moving_cobjs.erase(i++);} // remove if destroyed
		else {by_z1.push_back(make_pair(coll_objects.get_cobj(*i).d[2][0], *i)); ++i;} // otherwise try to drop it
	}
	sort(by_z1.begin(), by_z1.end()); // sort by z1 so that stacked cobjs work correctly (processed bottom to top)
	for (auto i = by_z1.begin(); i != by_z1.end(); ++i) {try_drop_movable_cobj(i->second);}
}

