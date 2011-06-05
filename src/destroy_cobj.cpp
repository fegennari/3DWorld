// 3D World - Collision Object Destruction Code
// by Frank Gennari
// 5/31/11
#include "3DWorld.h"
#include "mesh.h"
#include "csg.h"
#include "physics_objects.h"


int destroy_thresh(0);

extern int coll_id[];
extern obj_type object_types[];
extern obj_group obj_groups[];
extern vector<coll_obj> coll_objects;
extern vector<portal> portals;



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
		int const num(min(100, int((rand()%20 + 20)*(cts[i].volume/0.0007))));
		bool const shattered(cts[i].destroy >= SHATTERABLE);
		point fpos(pos);

		for (int o = 0; o < num; ++o) {
			vector3d velocity(cdir);

			if (shattered) {
				for (unsigned j = 0; j < 3; ++j) { // only accurate for COLL_CUBE
					fpos[j] = rand_uniform(cts[i].d[j][0], cts[i].d[j][1]); // generate inside of the shattered cobj's volume
				}
				vector3d const vadd(fpos - pos); // average cdir and direction from collision point to fragment location

				if (vadd.mag() > TOLERANCE) {
					velocity += vadd.get_norm()*(cdir_mag/vadd.mag());
					velocity *= 0.5;
				}
			}
			gen_fragment(fpos, velocity, 1.0, 0.5*rand_float(), cts[i].color, cts[i].tid, cts[i].tscale, shooter, shattered);
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


void get_all_connected(int cobj, set<int> &connected) {

	assert((unsigned)cobj < coll_objects.size());
	coll_obj const &c(coll_objects[cobj]);
	int const x1(get_xpos(c.d[0][0])), x2(get_xpos(c.d[0][1]));
	int const y1(get_ypos(c.d[1][0])), y2(get_ypos(c.d[1][1]));
			
	for (int y = max(0, y1); y <= min(MESH_Y_SIZE-1, y2); ++y) {
		for (int x = max(0, x1); x <= min(MESH_X_SIZE-1, x2); ++x) {
			vector<int> const &cvals(v_collision_matrix[y][x].cvals);

			for (vector<int>::const_iterator i = cvals.begin(); i != cvals.end(); ++i) {
				if (*i < 0 || *i == cobj)                  continue;
				if (connected.find(*i) != connected.end()) continue; // already processed
				assert((unsigned)*i < coll_objects.size());
				if (c.intersects_cobj(coll_objects[*i], TOLERANCE)) connected.insert(*i);
			}
		}
	}
}


unsigned subtract_cube(vector<coll_obj> &cobjs, vector<color_tid_vol> &cts, vector3d &cdir,
					   float x1, float x2, float y1, float y2, float z1, float z2, int min_destroy)
{
	//RESET_TIME;
	if (destroy_thresh >= EXPLODEABLE) return 0;
	csg_cube const cube(x1, x2, y1, y2, z1, z2);
	point center(cube.get_cube_center());
	if (cube.is_zero_area()) return 0;
	float const sub_volume(cube.get_volume());
	vector<int> indices, to_remove;
	vector<coll_obj> new_cobjs;
	vector<int> cvals;
	cdir = zero_vector;
	unsigned const cobjs_size(cobjs.size());
	float const maxlen(cube.max_len());
	bool const is_small(maxlen < HALF_DXY);
	unsigned ncobjs, last_cobj(0);

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
			indices.clear();

			for (unsigned j = 0; j < new_cobjs.size(); ++j) { // new objects
				//test_for_falling_cobj(new_cobjs[j], i);
				int const index(new_cobjs[j].add_coll_cobj()); // not sorted by alpha
				assert((size_t)index < cobjs.size());
				indices.push_back(index);
				volume -= cobjs[index].volume;
			}
			assert(volume >= -TOLERANCE); // usually > 0.0
			cts.push_back(color_tid_vol(cobjs[i], volume));
			cobjs[i].clear_internal_data(cobjs, indices, i);
			to_remove.push_back(i);
		}
		new_cobjs.clear();
	} // for k
	if (!to_remove.empty()) {
		//calc_visibility(SUN_SHADOW | MOON_SHADOW); // *** FIXME: what about updating (removing) mesh shadows? ***

		if (1) {
			// ***** BEGIN TESTING *****
			// FIXME: update cobj connectivity and make unconnected cobjs fall
			// * Use intersects_cobj() to find connected groups of shapes starting at all cobjs to_remove in a queue (open set)
			//     - Assert that all cobjs in indices are reached
			// * If we reach the following it is connected using is_anchored():
			// * Otherwise, add all objects in the group (closed set) to a falling objects group
			// * Retest the objects in each falling objects group every frame:
			//     If still unconnected, let them fall by time*grav_acc, each frame:
			//       - Remove using remove_coll_object(), maybe call clear_internal_data()
			//       - Re-add using add_coll_cobj()
			//     If connected, remove from the group (or remove the group)

			set<int> open, closed;
			copy(to_remove.begin(), to_remove.end(), inserter(open, open.begin()));
			bool anchored(0);

			while (!open.empty()) {
				int const cur(*open.begin());
				assert((unsigned)cur < coll_objects.size());

				if (coll_objects[cur].is_anchored()) {
					anchored = 1;
					break;
				}
				open.erase(cur);
				closed.insert(cur);
				set<int> connected;
				get_all_connected(cur, connected);

				for (set<int>::const_iterator i = connected.begin(); i != connected.end(); ++i) {
					if (closed.find(*i) != closed.end()) continue; // already processed
					open.insert(*i); // may already be there
				}
			}
			if (!anchored) {
				for (unsigned i = 0; i < indices.size(); ++i) {
					assert(closed.find(indices[i]) != closed.end()); // must have gotten to them all
				}
			}
			cout << "anchored: " << anchored << endl;
		}
		// ***** END TESTING *****

		for (unsigned i = 0; i < to_remove.size(); ++i) {
			remove_coll_object(to_remove[i]); // remove old collision object
		}
		cdir.normalize();
	}
	//PRINT_TIME("Subtract Cube");
	return to_remove.size();
}


bool is_pt_under_mesh(point const &p) {

	int const xpos(get_xpos(p.x)), ypos(get_ypos(p.y));
	if (point_outside_mesh(xpos, ypos)) return 0;
	return (p.z < mesh_height[ypos][xpos]);
}


int coll_obj::is_anchored() const {

	if (platform_id >= 0 || status != COLL_STATIC) return 0; // platforms and dynamic objects are never connecting
	if (fixed && destroy <= destroy_thresh)        return 2; // can't be destroyed, so it never moves

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
			// check cube edges for intersection with polygon
			return 2; // FIXME
		default: assert(0);
		}
		break;

	case COLL_CYLINDER:
		switch (c.type) {
		case COLL_CYLINDER:
			return dist_xy_less_than(points[0], c.points[0], (c.radius+radius));
		case COLL_SPHERE:
			return dist_xy_less_than(points[0], c.points[0], (c.radius+radius)); // FIXME: inexact (return 2?)
		case COLL_CYLINDER_ROT:
			if (line_line_dist(points[0], points[1], c.points[0], c.points[1]) > (radius + max(c.radius, c.radius2))) return 0;
			return 2; // FIXME
		case COLL_POLYGON:
			// could use line_intersect_cylinder() for each polygon edge
			return 2; // FIXME
		default: assert(0);
		}
		break;

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
		break;

	case COLL_CYLINDER_ROT:
		switch (c.type) {
		case COLL_CYLINDER_ROT:
			if (line_line_dist(points[0], points[1], c.points[0], c.points[1]) > (max(radius, radius2) + max(c.radius, c.radius2))) return 0;
			return 2; // FIXME
		case COLL_POLYGON:
			// could use line_intersect_cylinder() for each polygon edge
			return 2; // FIXME
		default: assert(0);
		}
		break;

	case COLL_POLYGON:
		switch (c.type) {
		case COLL_POLYGON:
			return 2; // FIXME - need to deal with thickness as well
		default: assert(0);
		}
		break;

	default:
		assert(0);
	}
	return 0;
}



