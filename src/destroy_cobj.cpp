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
vector<int> falling_cobjs;

extern float fticks, zmin;
extern int coll_id[];
extern obj_type object_types[];
extern obj_group obj_groups[];
extern vector<coll_obj> coll_objects;
extern vector<portal> portals;


void add_connect_waypoint_for_cobj(coll_obj &c);
void remove_waypoint_for_cobj(coll_obj &c);



void destroy_coll_objs(point const &pos, float damage, int shooter, bool big) {

	assert(damage >= 0.0);
	if (damage < 100.0) return;
	float const r((big ? 4.0 : 1.0)*sqrt(damage)/(rand_uniform(600.0, 750.0)));
	vector3d cdir;
	vector<color_tid_vol> cts;
	int const dmin((damage > 800.0) ? DESTROYABLE : ((damage > 200.0) ? SHATTERABLE : EXPLODEABLE));
	unsigned nrem(subtract_cube(coll_objects, cts, cdir, (pos.x-r),(pos.x+r),(pos.y-r),(pos.y+r),(pos.z-r),(pos.z+r), dmin));
	if (nrem == 0 || cts.empty()) return;
	float const cdir_mag(cdir.mag());

	for (unsigned i = 0; i < cts.size(); ++i) {
		if (cts[i].destroy == EXPLODEABLE) {
			float const val(float(pow(double(cts[i].volume), 1.0/3.0))), exp_damage(25000.0*val + 0.25*damage + 500.0);
			create_explosion(pos, shooter, 0, exp_damage, 10.0*val, BLAST_RADIUS, 0);
			gen_fire(pos, min(4.0, 12.0*val), shooter);
		}
		if (!cts[i].draw) continue;
		bool const shattered(cts[i].destroy >= SHATTERABLE);
		float size_scale(1.0), num_parts(0.0);
		float const thickness(cts[i].thickness);
		float const frag_radius(object_types[FRAGMENT].radius), avg_frag_dia(2.0*frag_radius), max_frag_dia(3.0*frag_radius);
		if (cts[i].volume < frag_radius*frag_radius*frag_radius) continue;

		if (shattered) {
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
		unsigned const num(min(100U, max((shattered ? 6U : 1U), unsigned(num_parts)))); // no more than 200
		//cout << "shattered: " << shattered << ", volume: " << cts[i].volume << ", num_parts: " << num_parts << ", num: " << num << ", ss: " << size_scale << endl;

		for (unsigned o = 0; o < num; ++o) {
			vector3d velocity(cdir);
			point fpos(pos);

			if (shattered || cts[i].unanchored) {
				for (unsigned j = 0; j < 3; ++j) { // only accurate for COLL_CUBE
					fpos[j] = rand_uniform(cts[i].d[j][0], cts[i].d[j][1]); // generate inside of the shattered cobj's volume
				}
				vector3d const vadd(fpos - pos); // average cdir and direction from collision point to fragment location

				if (vadd.mag() > TOLERANCE) {
					velocity += vadd.get_norm()*(cdir_mag/vadd.mag());
					velocity *= 0.5;
				}
			}
			gen_fragment(fpos, velocity, size_scale, 0.5*rand_float(), cts[i].color, cts[i].tid, cts[i].tscale, shooter, shattered);
		}
	} // for i
}


void add_portal(coll_obj const &c) {

	if (c.type == COLL_POLYGON) {
		assert(c.npoints == 3 || c.npoints == 4);
		portal p;

		for (int i = 0; i < c.npoints; ++i) {
			p.pts[i] = c.points[i]; // ignore thickness - use base polygon only
		}
		if (c.npoints == 3) p.pts[3] = p.pts[2]; // duplicate the last point
		portals.push_back(p);
	}
	else if (c.type == COLL_CUBE) {
		portal p;
		float max_area(0.0);
		
		for (unsigned i = 0; i < 6; ++i) { // choose enabled side with max area
			unsigned const dim(i>>1), dir(i&1), d0((dim+1)%3), d1((dim+2)%3);
			if (c.cp.surfs & EFLAGS[dim][dir]) continue; // disabled side
			float const area(fabs(c.d[d0][1] - c.d[d0][0])*fabs(c.d[d1][1] - c.d[d1][0]));

			if (area > max_area) {
				max_area = area;
				point pos;
				pos[dim] = c.d[dim][dir];

				for (unsigned n = 0; n < 4; ++n) {
					pos[d0] = c.d[d0][n<2];
					pos[d1] = c.d[d1][(n&1)^(n<2)];
					p.pts[n] = pos;
				}
			}
		}
		if (max_area > 0.0) portals.push_back(p);
	}
	// else other types are not supported yet (cylinder, sphere)
}


void get_all_connected(int cobj, set<int> &connected, set<int> const &closed) {

	assert((unsigned)cobj < coll_objects.size());
	coll_obj const &c(coll_objects[cobj]);
	int const x1(get_xpos(c.d[0][0])), x2(get_xpos(c.d[0][1]));
	int const y1(get_ypos(c.d[1][0])), y2(get_ypos(c.d[1][1]));
			
	for (int y = max(0, y1); y <= min(MESH_Y_SIZE-1, y2); ++y) {
		for (int x = max(0, x1); x <= min(MESH_X_SIZE-1, x2); ++x) {
			vector<int> const &cvals(v_collision_matrix[y][x].cvals);

			for (vector<int>::const_iterator i = cvals.begin(); i != cvals.end(); ++i) {
				if (*i < 0 || *i == cobj)                   continue;
				assert((unsigned)*i < coll_objects.size());
				if (closed.find(*i)    != closed.end()   )  continue; // closed  - skip
				if (connected.find(*i) != connected.end())  continue; // already processed - skip
				if (coll_objects[*i].status != COLL_STATIC) continue; // not static
				if (c.intersects_cobj(coll_objects[*i], TOLERANCE) == 1) connected.insert(*i);
			}
		}
	}
}


void check_cobjs_anchored(vector<int> to_check, set<int> anchored[2]) {

	for (vector<int>::const_iterator j = to_check.begin(); j != to_check.end(); ++j) {
		if ((*j) < 0) continue; // skip
		if (anchored[0].find(*j) != anchored[0].end()) continue; // already known to be unanchored
		if (anchored[1].find(*j) != anchored[1].end()) continue; // already known to be anchored

		// perform a graph search until we find an anchored cobj or we run out of cobjs
		bool is_anchored(0);
		set<int> open, closed;
		open.insert(*j);

		while (!open.empty()) {
			int const cur(*open.begin());
			closed.insert(cur);
			assert(anchored[0].find(cur) == anchored[0].end()); // requires that intersects_cobj() be symmetric

			if (anchored[1].find(cur) != anchored[1].end() || coll_objects[cur].is_anchored()) {
				is_anchored = 1;
				break;
			}
			open.erase(cur);
			get_all_connected(cur, open, closed);
		}
		// everything in the closed set has the same is_anchored state and can be cached
		copy(closed.begin(), closed.end(), inserter(anchored[is_anchored], anchored[is_anchored].begin()));
	} // for j
}


unsigned subtract_cube(vector<coll_obj> &cobjs, vector<color_tid_vol> &cts, vector3d &cdir,
					   float x1, float x2, float y1, float y2, float z1, float z2, int min_destroy)
{
	//RESET_TIME;
	if (destroy_thresh >= EXPLODEABLE) return 0;
	csg_cube const cube(x1, x2, y1, y2, z1, z2);
	point center(cube.get_cube_center());
	if (cube.is_zero_area()) return 0;
	vector<int> indices, to_remove;
	vector<coll_obj> new_cobjs;
	vector<int> cvals;
	cdir = zero_vector;
	unsigned const cobjs_size(cobjs.size());
	float const maxlen(cube.max_len());
	bool const is_small(maxlen < HALF_DXY);
	unsigned ncobjs, last_cobj(0), extra_removed(0);

	if (is_small) { // not much faster
		int const xpos(get_xpos(center.x)), ypos(get_ypos(center.y));
		if (point_outside_mesh(xpos, ypos)) return 0;
		cvals  = v_collision_matrix[ypos][xpos].cvals; // make a copy because cvals can change during iteration
		ncobjs = cvals.size(); // so as not to retest newly created subcubes
	}
	else {
		ncobjs = cobjs_size; // so as not to retest newly created subcubes
	}
	for (unsigned k = 0; k < ncobjs; ++k) {
		unsigned const i(is_small ? cvals[k] : k);
		assert((size_t)i < cobjs_size);
		if (cobjs[i].status != COLL_STATIC /*|| !cobjs[i].fixed*/) continue; // require fixed cobjs? exclude platforms (but they seem to work)?
		bool const is_cylinder(cobjs[i].is_cylinder()), is_cube(cobjs[i].type == COLL_CUBE), csg_obj(is_cube || is_cylinder);
		int const D(cobjs[i].destroy);
		if (D <= max(destroy_thresh, (min_destroy-1))) continue;
		bool const shatter(D >= SHATTERABLE);
		if (!shatter && !csg_obj)         continue;
		csg_cube const cube2(cobjs[i], !csg_obj);
		if (!cube2.intersects(cube, 0.0)) continue; // no intersection
		//if (is_cube && !cube2.contains_pt(cube.get_cube_center())) {} // check for non-destroyable cobj between center and cube2?
		float volume(cobjs[i].volume);
		bool no_new_cubes(shatter || volume < TOLERANCE);

		if (!csg_obj || subtract_cobj(new_cobjs, cube, cobjs[i]) || (shatter && is_cylinder)) {
			if (no_new_cubes) new_cobjs.clear(); // completely destroyed
			if (is_cube)      cdir += cube2.closest_side_dir(center); // inexact
			if (D == SHATTER_TO_PORTAL) add_portal(cobjs[i]);
			remove_waypoint_for_cobj(cobjs[i]);
			indices.clear();

			for (unsigned j = 0; j < new_cobjs.size(); ++j) { // new objects
				int const index(new_cobjs[j].add_coll_cobj()); // not sorted by alpha
				assert((size_t)index < cobjs.size());
				indices.push_back(index);
				volume -= cobjs[index].volume;
				add_connect_waypoint_for_cobj(cobjs[index]);
			}
			assert(volume >= -TOLERANCE); // usually > 0.0
			cts.push_back(color_tid_vol(cobjs[i], volume, cobjs[i].calc_min_dim(), 0));
			cobjs[i].clear_internal_data(cobjs, indices, i);
			to_remove.push_back(i);
			// FIXME: adjust lightmap pflow value so that smoke can flow through the hole?
		}
		new_cobjs.clear();
	} // for k
	if (!to_remove.empty()) {
		//calc_visibility(SUN_SHADOW | MOON_SHADOW); // *** FIXME: what about updating (removing) mesh shadows? ***

		for (unsigned i = 0; i < to_remove.size(); ++i) {
			remove_coll_object(to_remove[i]); // remove old collision object
		}
		if (LET_COBJS_FALL || REMOVE_UNANCHORED) {
			set<int> anchored[2]; // {unanchored, anchored}

			for (unsigned i = 0; i < to_remove.size(); ++i) { // cobjs in to_remove are freed but still valid
				set<int> start_set;
				get_all_connected(to_remove[i], start_set, set<int>());
				vector<int> const start(start_set.begin(), start_set.end());
				check_cobjs_anchored(start, anchored);
			}

			// check that each sub-cobj index is in exactly one of anchored[{0,1}]
			for (unsigned i = 0; i < indices.size(); ++i) {
				assert((anchored[0].find(indices[i]) == anchored[0].end()) != (anchored[1].find(indices[i]) == anchored[1].end()));
			}
			if (REMOVE_UNANCHORED) {
				indices.clear();

				for (set<int>::const_iterator i = anchored[0].begin(); i != anchored[0].end(); ++i) {
					if (cobjs[*i].destroy <= max(destroy_thresh, (min_destroy-1))) continue; // can't destroy (can't get here?)
					cts.push_back(color_tid_vol(cobjs[*i], cobjs[*i].volume, cobjs[*i].calc_min_dim(), 1));
					cobjs[*i].clear_internal_data(cobjs, indices, *i);
					remove_coll_object(*i);
					++extra_removed;
				}
			}
			else if (LET_COBJS_FALL) {
				copy(anchored[0].begin(), anchored[0].end(), back_inserter(falling_cobjs));
			}
		}
		cdir.normalize();
	}
	//PRINT_TIME("Subtract Cube");
	return (to_remove.size() + extra_removed);
}


void check_falling_cobjs() {

	// FIXME: add shadows when a cobj begins to shadow a new target cobj
	// FIXME: add velocity/acceleration
	// FIXME: fix texture offset
	if (falling_cobjs.empty()) return; // nothing to do
	float const dz(-0.001*fticks); // FIXME: add velocity/acceleration due to gravity
	set<int> anchored[2]; // {unanchored, anchored}

	for (vector<int>::iterator i = falling_cobjs.begin(); i != falling_cobjs.end(); ++i) {
		if (*i < 0) continue; // skip
		assert((unsigned)(*i) < coll_objects.size());
	
		if (coll_objects[*i].status != COLL_STATIC) {
			*i = -1; // disable
			continue;
		}
		// translate, add the new, then remove the old
		vector<int> indices;
		//indices.push_back(index); // if only we could do this first...
		coll_objects[*i].clear_internal_data(coll_objects, indices, *i);
		coll_objects[*i].clear_lightmap(0); // need to do this first, before the copy
		coll_obj cobj(coll_objects[*i]); // make a copy
		cobj.shift_by(point(0.0, 0.0, dz), 1); // translate down
		int const index(cobj.add_coll_cobj());
		remove_coll_object(*i);
		assert(*i != index);
		*i = index;
	}
	check_cobjs_anchored(falling_cobjs, anchored);
	falling_cobjs.resize(0);
	copy(anchored[0].begin(), anchored[0].end(), back_inserter(falling_cobjs));
}


bool is_pt_under_mesh(point const &p) {

	int const xpos(get_xpos(p.x)), ypos(get_ypos(p.y));
	if (point_outside_mesh(xpos, ypos)) return 0;
	return (p.z < mesh_height[ypos][xpos]);
}


int coll_obj::is_anchored() const {

	if (platform_id >= 0 || status != COLL_STATIC) return 0; // platforms and dynamic objects are never connecting
	if (fixed && destroy <= destroy_thresh)        return 2; // can't be destroyed, so it never moves
	if (d[2][0] <= min(zmin, czmin))               return 1; // below the scene

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
// 15 total: 7 complete, 5 partial, 3 unwritten
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
			return poly_cylin_int(c, *this);
		default: assert(0);
		}

	case COLL_POLYGON:
		assert(c.type == COLL_POLYGON);
		for (int i = 0; i < c.npoints; ++i) { // FIXME: use toler (for adjacent roof polygons)
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



