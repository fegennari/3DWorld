// 3D World - Collision Object Destruction Code
// by Frank Gennari
// 5/31/11
#include "3DWorld.h"
#include "mesh.h"
#include "csg.h"
#include "physics_objects.h"


bool const LET_COBJS_FALL    = 0;
bool const REMOVE_UNANCHORED = 1;

int destroy_thresh(0);
vector<unsigned> falling_cobjs;

extern bool scene_dlist_invalid;
extern float tstep, zmin, base_gravity;
extern int cobj_counter, coll_id[];
extern obj_type object_types[];
extern obj_group obj_groups[];
extern vector<coll_obj> coll_objects;
extern vector<portal> portals;
extern coll_cell_opt_batcher cco_batcher;


unsigned subtract_cube(vector<color_tid_vol> &cts, vector3d &cdir, csg_cube const &cube, int destroy_thresh);


// **************** Cobj Destroy Code ****************


void destroy_coll_objs(point const &pos, float damage, int shooter, bool big) {

	//RESET_TIME;
	assert(damage >= 0.0);
	if (damage < 100.0) return;
	float const radius((big ? 4.0 : 1.0)*sqrt(damage)/650.0);
	vector3d cdir;
	vector<color_tid_vol> cts;
	int const dmin((damage > 800.0) ? DESTROYABLE : ((damage > 200.0) ? SHATTERABLE : EXPLODEABLE));
	csg_cube cube(pos.x, pos.x, pos.y, pos.y, pos.z, pos.z);
	cube.expand_by(radius);
	unsigned nrem(subtract_cube(cts, cdir, cube, dmin));
	if (nrem == 0 || cts.empty()) return; // nothing removed
	int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));

	if (!point_outside_mesh(xpos, ypos)) {
		// check for waypoints that can be added near this cube (at the center only)
		vector<int> const &cvals(v_collision_matrix[ypos][xpos].cvals);

		for (vector<int>::const_iterator i = cvals.begin(); i != cvals.end(); ++i) {
			if (*i < 0) continue;
			assert((unsigned)*i < coll_objects.size());
			if (coll_objects[*i].waypt_id < 0) coll_objects[*i].add_connect_waypoint(); // slow
		}
	}

	// update voxel pflow map for removal
	update_flow_for_voxels(cube);

	for (unsigned i = 0; i < cts.size(); ++i) {
		if (cts[i].destroy >= SHATTERABLE || cts[i].unanchored) update_flow_for_voxels(cts[i]);
	}

	// create fragments
	float const cdir_mag(cdir.mag());

	for (unsigned i = 0; i < cts.size(); ++i) {
		if (cts[i].destroy == EXPLODEABLE) {
			float const val(float(pow(double(cts[i].volume), 1.0/3.0))), exp_damage(25000.0*val + 0.25*damage + 500.0);
			create_explosion(pos, shooter, 0, exp_damage, 10.0*val, BLAST_RADIUS, 0);
			gen_fire(pos, min(4.0, 12.0*val), shooter);
		}
		if (!cts[i].draw) continue;
		float size_scale(1.0), num_parts(0.0);
		float const thickness(cts[i].thickness);
		bool const shattered(cts[i].destroy >= SHATTERABLE);
		bool const tri_fragments(shattered || cts[i].is_2d); // thin polys shatter to triangle fragments
		float const frag_radius(object_types[FRAGMENT].radius), avg_frag_dia(2.0*frag_radius), max_frag_dia(3.0*frag_radius);
		if (cts[i].volume < (tri_fragments ? MIN_POLY_THICK : frag_radius)*frag_radius*frag_radius) continue;

		if (tri_fragments) {
			float const sll(cts[i].second_largest_len());
			if (sll < 1.2*max_frag_dia) size_scale *= sll/max_frag_dia;
			float const dia(size_scale*avg_frag_dia);
			num_parts = cts[i].volume/(thickness*dia*dia);
		}
		else {
			if (thickness < 1.2*max_frag_dia) size_scale *= thickness/max_frag_dia;
			float const dia(size_scale*avg_frag_dia);
			num_parts = cts[i].volume/(dia*dia*dia);
		}
		if (size_scale < 0.1) continue;
		unsigned const num(min(100U, max(((tri_fragments && !cts[i].is_2d) ? 6U : 1U), unsigned(num_parts)))); // no more than 200
		//cout << "shattered: " << shattered << ", tri: " << tri_fragments << ", volume: " << cts[i].volume << ", num_parts: " << num_parts << ", num: " << num << ", ss: " << size_scale << endl;
		bool const non_csg(shattered || cts[i].unanchored);
		csg_cube frag_cube(cts[i]);
		if (!non_csg && !cube.cube_intersection(frag_cube, frag_cube)) frag_cube = cts[i]; // intersect frag_cube with cube (should pass)

		for (unsigned o = 0; o < num; ++o) {
			vector3d velocity(cdir);
			point fpos(frag_cube.gen_rand_pt_in_cube()); // only accurate for COLL_CUBE

			if (non_csg) {
				vector3d const vadd(fpos - pos); // average cdir and direction from collision point to fragment location

				if (vadd.mag() > TOLERANCE) {
					velocity += vadd.get_norm()*(cdir_mag/vadd.mag());
					velocity *= 0.5;
				}
			}
			gen_fragment(fpos, velocity, size_scale, 0.5*rand_float(), cts[i].color, cts[i].tid, cts[i].tscale, shooter, tri_fragments);
		}
	} // for i
	//PRINT_TIME("Destroy Cobjs");
}


void coll_obj::create_portal() const {

	switch (type) {
	case COLL_POLYGON:
		{
			assert(npoints == 3 || npoints == 4);
			portal p;

			for (int i = 0; i < npoints; ++i) {
				p.pts[i] = points[i]; // ignore thickness - use base polygon only
			}
			if (npoints == 3) p.pts[3] = p.pts[2]; // duplicate the last point
			portals.push_back(p);
		}
	case COLL_CUBE:
		{
			portal p;
			float max_area(0.0);
		
			for (unsigned i = 0; i < 6; ++i) { // choose enabled side with max area
				unsigned const dim(i>>1), dir(i&1), d0((dim+1)%3), d1((dim+2)%3);
				if (cp.surfs & EFLAGS[dim][dir]) continue; // disabled side
				float const area(fabs(d[d0][1] - d[d0][0])*fabs(d[d1][1] - d[d1][0]));

				if (area > max_area) {
					max_area = area;
					point pos;
					pos[dim] = d[dim][dir];

					for (unsigned n = 0; n < 4; ++n) {
						pos[d0] = d[d0][n<2];
						pos[d1] = d[d1][(n&1)^(n<2)];
						p.pts[n] = pos;
					}
				}
			}
			if (max_area > 0.0) portals.push_back(p);
		}
	default:
		assert(0); // other types are not supported yet (cylinder, sphere)
	}
}


void get_all_connected(unsigned cobj, vector<unsigned> &out) {

	assert(cobj < coll_objects.size());
	get_intersecting_cobjs_tree(coll_objects[cobj], out, cobj, TOLERANCE, 0, 1, cobj);
}


void check_cobjs_anchored(vector<unsigned> to_check, set<unsigned> anchored[2]) {

	vector<unsigned> out;

	for (vector<unsigned>::const_iterator j = to_check.begin(); j != to_check.end(); ++j) {
		if (anchored[0].find(*j) != anchored[0].end()) continue; // already known to be unanchored
		if (anchored[1].find(*j) != anchored[1].end()) continue; // already known to be anchored

		if (coll_objects[*j].is_anchored()) {
			anchored[1].insert(*j);
			continue;
		}

		// perform a graph search until we find an anchored cobj or we run out of cobjs
		bool is_anchored(0);
		vector<unsigned> open, closed;
		open.push_back(*j);
		++cobj_counter;
		assert(coll_objects[*j].counter != cobj_counter);
		coll_objects[*j].counter = cobj_counter;

		while (!open.empty()) {
			unsigned const cur(open.back());
			open.pop_back();
			closed.push_back(cur);
			//assert(anchored[0].find(cur) == anchored[0].end()); // requires that intersects_cobj() be symmetric
			out.resize(0);
			get_all_connected(cur, out);

			for (vector<unsigned>::const_iterator i = out.begin(); i != out.end(); ++i) {
				assert(*i >= 0 && *i != cur);
				assert(coll_objects[*i].counter != cobj_counter); // may be too strong - we might want to allow duplicates and just continue here
				open.push_back(*i); // need to do this first

				if (anchored[1].find(*i) != anchored[1].end() || coll_objects[*i].is_anchored()) {
					is_anchored = 1;
					break;
				}
				coll_objects[*i].counter = cobj_counter;
			}
			if (is_anchored) break;
		}
		// everything in the closed set has the same is_anchored state and can be cached
		copy(closed.begin(), closed.end(), inserter(anchored[is_anchored], anchored[is_anchored].begin()));
		
		if (is_anchored) { // all open is anchored as well
			copy(open.begin(), open.end(), inserter(anchored[is_anchored], anchored[is_anchored].begin()));
		}
		else {
			assert(open.empty());
		}
	} // for j
}


void add_to_falling_cobjs(set<unsigned> const &ids) {

	for (set<unsigned>::const_iterator i = ids.begin(); i != ids.end(); ++i) {
		assert(*i < coll_objects.size());
		coll_objects[*i].falling = 1;
		falling_cobjs.push_back(*i);
	}
}


void invalidate_static_cobjs() {

	build_cobj_tree(0, 0);
	scene_dlist_invalid = 1;
}


unsigned subtract_cube(vector<color_tid_vol> &cts, vector3d &cdir, csg_cube const &cube, int min_destroy) {

	if (destroy_thresh >= EXPLODEABLE) return 0;
	if (cube.is_zero_area())           return 0;
	//RESET_TIME;
	vector<coll_obj> &cobjs(coll_objects); // so we don't have to rename everything and can keep the shorter code
	point center(cube.get_cube_center());
	float const clip_cube_colume(cube.get_volume());
	vector<int> just_added, to_remove;
	vector<coll_obj> new_cobjs;
	cdir = zero_vector;
	vector<cube_t> mod_cubes;
	mod_cubes.push_back(cube);
	build_moving_cobj_tree();
	vector<unsigned> int_cobjs;
	get_intersecting_cobjs_tree(cube, int_cobjs, -1, 0.0, 0, 0, -1);
	cco_batcher.begin_batch();

	// determine affected cobjs
	for (unsigned k = 0; k < int_cobjs.size(); ++k) {
		unsigned const i(int_cobjs[k]);
		assert(i < cobjs.size());
		assert(cobjs[i].status == COLL_STATIC);
		// require fixed cobjs? platforms work now
		int const D(cobjs[i].destroy);
		if (D <= max(destroy_thresh, (min_destroy-1))) continue;
		bool const is_cylinder(cobjs[i].is_cylinder()), is_cube(cobjs[i].type == COLL_CUBE), is_polygon(cobjs[i].type == COLL_POLYGON);
		bool const csg_obj(is_cube || is_cylinder || is_polygon), shatter(D >= SHATTERABLE);
		if (!shatter && !csg_obj)            continue;
		if (!cobjs[i].intersects(cube, 0.0)) continue; // no intersection
		csg_cube const cube2(cobjs[i], 1);
		//if (is_cube && !cube2.contains_pt(cube.get_cube_center())) {} // check for non-destroyable cobj between center and cube2?
		float volume(cobjs[i].volume);
		float const min_volume(0.01*min(volume, clip_cube_colume));
		float const int_volume(cube2.get_overlap_volume(cube));
		bool no_new_cobjs(shatter || volume < TOLERANCE);

		if (is_cube && !shatter && int_volume < min_volume) { // don't remove tiny bits from cobjs
			cube.unset_intersecting_edge_flags(cobjs[i]);
			continue;
		}
		if (shatter || subtract_cobj(new_cobjs, cube, cobjs[i], 1)) {
			if (no_new_cobjs) new_cobjs.clear(); // completely destroyed
			if (is_cube)      cdir += cube2.closest_side_dir(center); // inexact
			if (D == SHATTER_TO_PORTAL) cobjs[i].create_portal();

			for (unsigned j = 0; j < new_cobjs.size(); ++j) { // new objects
				int const index(new_cobjs[j].add_coll_cobj()); // not sorted by alpha
				assert(index >= 0 && (size_t)index < cobjs.size());
				just_added.push_back(index);
				volume -= cobjs[index].volume;
			}
			if (is_polygon) volume = max(0.0f, volume); // FIXME: remove this when polygon splitting is correct
			assert(volume >= -TOLERANCE); // usually > 0.0
			cts.push_back(color_tid_vol(cobjs[i], volume, cobjs[i].calc_min_dim(), 0));
			cobjs[i].clear_internal_data();
			to_remove.push_back(i);
			if (shatter) mod_cubes.push_back(cobjs[i]);
		}
		new_cobjs.clear();
	} // for k
	cco_batcher.end_batch();

	// remove destroyed cobjs
	for (unsigned i = 0; i < to_remove.size(); ++i) {
		cobjs[to_remove[i]].remove_waypoint();
		remove_coll_object(to_remove[i]); // remove old collision object
	}
	if (!to_remove.empty()) {
		invalidate_static_cobjs(); // after destroyed cobj removal
		build_moving_cobj_tree();
	}

	// add new waypoints (after build_cobj_tree and end_batch)
	for (unsigned i = 0; i < just_added.size(); ++i) {
		cobjs[just_added[i]].add_connect_waypoint(); // slow
	}

	// process unanchored cobjs
	if (LET_COBJS_FALL || REMOVE_UNANCHORED) {
		//RESET_TIME;
		set<unsigned> anchored[2]; // {unanchored, anchored}

		for (unsigned i = 0; i < to_remove.size(); ++i) { // cobjs in to_remove are freed but still valid
			vector<unsigned> start;
			++cobj_counter;
			assert(coll_objects[to_remove[i]].counter != cobj_counter);
			coll_objects[to_remove[i]].counter = cobj_counter;
			get_all_connected(to_remove[i], start);
			check_cobjs_anchored(start, anchored);
		}
#if 0
		// additional optional error check that no cobj is both anchored and unanchored - can fail for polygons due to inexact intersection test
		for (set<unsigned>::const_iterator i = anchored[0].begin(); i != anchored[0].end(); ++i) {
			assert(anchored[1].find(*i) == anchored[1].end());
		}
#endif
		if (REMOVE_UNANCHORED) {
			for (set<unsigned>::const_iterator i = anchored[0].begin(); i != anchored[0].end(); ++i) {
				if (cobjs[*i].destroy <= max(destroy_thresh, (min_destroy-1))) continue; // can't destroy (can't get here?)
				cts.push_back(color_tid_vol(cobjs[*i], cobjs[*i].volume, cobjs[*i].calc_min_dim(), 1));
				cobjs[*i].clear_internal_data();
				mod_cubes.push_back(cobjs[*i]);
				cobjs[*i].remove_waypoint();
				remove_coll_object(*i);
				to_remove.push_back(*i);
			}
		}
		else if (LET_COBJS_FALL) {
			add_to_falling_cobjs(anchored[0]);
		}
		//PRINT_TIME("Check Anchored");
	}
	if (!to_remove.empty()) {
		//calc_visibility(SUN_SHADOW | MOON_SHADOW); // FIXME: what about updating (removing) mesh shadows?
		cdir.normalize();
	}
	//PRINT_TIME("Subtract Cube");
	return to_remove.size();
}


void check_falling_cobjs() {

	if (falling_cobjs.empty()) return; // nothing to do
	//RESET_TIME;
	float const accel(-0.5*base_gravity*GRAVITY*tstep); // half gravity
	set<unsigned> anchored[2]; // {unanchored, anchored}

	for (unsigned i = 0; i < falling_cobjs.size(); ++i) {
		unsigned const ix(falling_cobjs[i]);
		assert(ix < coll_objects.size());
	
		if (coll_objects[ix].status != COLL_STATIC) { // disable
			falling_cobjs[i] = falling_cobjs.back();
			falling_cobjs.pop_back();
			--i; // wraparound is ok
			continue;
		}
		// translate, add the new, then remove the old
		coll_objects[ix].clear_internal_data();
		coll_obj cobj(coll_objects[ix]); // make a copy
		cobj.v_fall += accel; // terminal velocity?
		cobj.shift_by(point(0.0, 0.0, tstep*cobj.v_fall), 1); // translate down
		int const index(cobj.add_coll_cobj());
		remove_coll_object(ix);
		assert(ix != index);
		falling_cobjs[i] = index;
	}
	vector<unsigned> last_falling(falling_cobjs);
	sort(last_falling.begin(), last_falling.end());
	build_moving_cobj_tree();
	check_cobjs_anchored(falling_cobjs, anchored);
	falling_cobjs.resize(0);
	add_to_falling_cobjs(anchored[0]);
	if (falling_cobjs != last_falling) invalidate_static_cobjs();
	//PRINT_TIME("Check Falling Cobjs");
}


// **************** Cobj Connectivity Code ****************


bool is_pt_under_mesh(point const &p) {

	//return is_under_mesh(p); // too slow?
	int const xpos(get_xpos(p.x)), ypos(get_ypos(p.y));
	if (point_outside_mesh(xpos, ypos)) return 0;
	return (p.z < mesh_height[ypos][xpos]);
}


int coll_obj::is_anchored() const {

	if (platform_id >= 0 || status != COLL_STATIC) return 0; // platforms and dynamic objects are never connecting
	if (destroy <= destroy_thresh)                 return 2; // can't be destroyed, so it never moves
	if (d[2][0] <= min(zbottom, czmin))            return 1; // below the scene
	if (d[2][0] > ztop)                            return 0; // above the mesh

	switch (type) {
	case COLL_CUBE:
		{
			int const x1(get_xpos(d[0][0])), x2(get_xpos(d[0][1]));
			int const y1(get_ypos(d[1][0])), y2(get_ypos(d[1][1]));
			
			for (int y = max(0, y1); y <= min(MESH_Y_SIZE-1, y2); ++y) {
				for (int x = max(0, x1); x <= min(MESH_X_SIZE-1, x2); ++x) {
					if (d[2][0] < mesh_height[y][x]) return 1;
				}
			}
			return 0;
		}
	case COLL_SPHERE:
		return is_pt_under_mesh((points[0] - vector3d(0.0, 0.0, radius)));
	case COLL_CYLINDER: // should really test the entire top/bottom surface
	case COLL_CYLINDER_ROT:
		return is_pt_under_mesh(points[0]) || is_pt_under_mesh(points[1]);
	case COLL_POLYGON: // should really test the entire surface(s)
		for (int i = 0; i < npoints; ++i) {
			if (is_pt_under_mesh(points[i])) return 1;

			if (thickness > MIN_POLY_THICK2) {
				if (is_pt_under_mesh(points[i] + norm*(0.5*thickness))) return 1;
				if (is_pt_under_mesh(points[i] - norm*(0.5*thickness))) return 1;
			}
		}
		return 0;
	default:
		assert(0);
	}
	return 0;
}


int cylin_cylin_int(coll_obj const &c1, coll_obj const &c2) {

	if (line_line_dist(c2.points[0], c2.points[1], c1.points[0], c1.points[1]) > (max(c2.radius, c2.radius2) + max(c1.radius, c1.radius2))) return 0;
	if (c1.line_intersect(c2.points[0], c2.points[1])) return 1;
	if (c2.line_intersect(c1.points[0], c1.points[1])) return 1;
	return 2; // FIXME: finish
}


int poly_cylin_int(coll_obj const &p, coll_obj const &c) {

	if (p.line_intersect(c.points[0], c.points[1])) return 1;

	for (int i = 0; i < p.npoints; ++i) {
		if (c.line_intersect(p.points[i], p.points[(i+1)%p.npoints])) return 1; // definite intersection
	}
	return 2; // FIXME: finish
}


// 0: no intersection, 1: intersection, 2: maybe intersection (incomplete)
// 15 total: 7 complete, 8 partial
int coll_obj::intersects_cobj(coll_obj const &c, float toler) const {

	if (c.type < type) return c.intersects_cobj(*this, toler); // swap arguments
	if (!intersects(c, toler)) return 0; // cube-cube intersection

	// c.type >= type
	switch (type) {
	case COLL_CUBE:
		switch (c.type) {
		case COLL_CUBE:
			return 1; // as simple as that
		case COLL_CYLINDER:
			return circle_rect_intersect(c.points[0], c.radius, *this);
		case COLL_SPHERE:
			return sphere_cube_intersect(c.points[0], c.radius, *this);
		case COLL_CYLINDER_ROT:
			if (check_line_clip(c.points[0], c.points[1], d)) return 1; // definite intersection
			return 2; // FIXME
		case COLL_POLYGON:
			for (int i = 0; i < c.npoints; ++i) {
				if (check_line_clip(c.points[i], c.points[(i+1)%c.npoints], d)) return 1; // definite intersection
			}
			if (c.thickness > MIN_POLY_THICK2) { // test extruded (3D) polygon
				static vector<point> pts[2];
				gen_poly_planes(c.points, c.npoints, c.norm, c.thickness, pts);
				
				for (unsigned j = 0; j < 2; ++j) {
					for (unsigned i = 0; i < pts[j].size(); ++i) {
						if (check_line_clip(pts[j][i], pts[j][(i+1)%pts[j].size()], d)) return 1; // definite intersection
					}
				}
				// call sphere_ext_poly_intersect?
				// call something like csg_cube::subtract_from_polygon()?
				return 0; // FIXME - close, but need to handle cube completely insde of a thick polygon
			}
			return 0;
		default: assert(0);
		}

	case COLL_CYLINDER:
		switch (c.type) {
		case COLL_CYLINDER:
			return dist_xy_less_than(points[0], c.points[0], (c.radius+radius));
		case COLL_SPHERE:
			return dist_xy_less_than(points[0], c.points[0], (c.radius+radius)); // FIXME: inexact (return 2?)
		case COLL_CYLINDER_ROT:
			return cylin_cylin_int(c, *this);
		case COLL_POLYGON:
			return poly_cylin_int(c, *this);
		default: assert(0);
		}

	case COLL_SPHERE:
		switch (c.type) {
		case COLL_SPHERE:
			return dist_less_than(points[0], c.points[0], (c.radius+radius));
		case COLL_CYLINDER_ROT:
			return sphere_intersect_cylinder(points[0], radius, c.points[0], c.points[1], c.radius, c.radius2);
		case COLL_POLYGON:
			return sphere_ext_poly_intersect(c.points, c.npoints, c.norm, points[0], radius, c.thickness, MIN_POLY_THICK2);
		default: assert(0);
		}

	case COLL_CYLINDER_ROT:
		switch (c.type) {
		case COLL_CYLINDER_ROT:
			return cylin_cylin_int(c, *this);
		case COLL_POLYGON:
			return poly_cylin_int (c, *this);
		default: assert(0);
		}

	case COLL_POLYGON:
		assert(c.type == COLL_POLYGON);
		{
			float const poly_toler(max(toler, (thickness + c.thickness)*(1.0f - fabs(dot_product(norm, c.norm)))));

			if (poly_toler > 0.0) { // use toler for edge adjacency tests (for adjacent roof polygons, sponza polygons, etc.)
				for (int i = 0; i < c.npoints; ++i) {
					for (int j = 0; j < npoints; ++j) {
						if (dist_less_than(points[j], c.points[i], poly_toler)) return 1;
					}
				}
				for (int i = 0; i < c.npoints; ++i) {
					for (int j = 0; j < npoints; ++j) {
						if (pt_line_seg_dist_less_than(c.points[i],   points[j],   points[(j+1)%  npoints], poly_toler)) return 1;
						if (pt_line_seg_dist_less_than(  points[j], c.points[i], c.points[(i+1)%c.npoints], poly_toler)) return 1;
					}
				}
			}
		}
		for (int i = 0; i < c.npoints; ++i) {
			if (line_intersect(c.points[i], c.points[(i+1)%c.npoints])) return 1;
		}
		for (int i = 0; i <   npoints; ++i) {
			if (line_intersect(  points[i],   points[(i+1)%  npoints])) return 1;
		}
		// call sphere_ext_poly_intersect?
		return 0; // FIXME - close, but need to handle one polygon completely insde of a thick polygon

	default:
		assert(0);
	}
	return 0;
}



