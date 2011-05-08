// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/13/02

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"


extern int DISABLE_WATER, camera_change;
extern float temperature, zmin;
extern obj_type object_types[];
extern dwobject def_objects[];
extern vector<coll_obj> coll_objects;


class waypoint_builder {

	vector<point> &waypoints;
	float radius;

	bool is_waypoint_valid(point const &pos, int coll_id) const {
		if (!is_over_mesh(pos) || pos.z < zmin) return 0;
		float const mesh_zval(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 0));
		if (pos.z - radius < mesh_zval) return 0; // bottom of smiley is under the mesh - should use a mesh waypoint here
		dwobject obj(def_objects[SMILEY]); // create a fake temporary smiley object
		obj.pos     = pos;
		obj.pos.z  += SMALL_NUMBER; // move up slightly
		obj.coll_id = coll_id; // ignore collisions with the current object
		return !obj.check_vert_collision(0, 0, 0); // return true if no collision
	}

	void add_if_valid(point const &pos, int coll_id) {
		if (is_waypoint_valid(pos, coll_id)) waypoints.push_back(pos);
	}

	void add_waypoint_rect(float x1, float y1, float x2, float y2, float z, int coll_id) {
		if (min(x2-x1, y2-y1) < radius) return; // too small to stand on
		point const center(0.5*(x1+x2), 0.5*(y1+y2), z+radius);
		add_if_valid(center, coll_id);
		// *** WRITE ***
	}

	void add_waypoint_circle(point const &p, float r, int coll_id) {
		if (r < radius) return; // too small to stand on
		add_if_valid(p + point(0.0, 0.0, radius), coll_id);
		// *** WRITE ***
	}

	void add_waypoint_triangle(point const &p1, point const &p2, point const &p3, int coll_id) {
		if (dist_less_than(p1, p2, 2*radius) || dist_less_than(p2, p3, 2*radius) || dist_less_than(p3, p1, 2*radius)) return; // too small to stand on
		point const center(((p1 + p2 + p3) / 3.0) + point(0.0, 0.0, radius));
		add_if_valid(center, coll_id);
		// *** WRITE ***
	}

public:
	waypoint_builder(vector<point> &waypoints_)
		: waypoints(waypoints_), radius(object_types[SMILEY].radius) {}
	
	void add_cobj_waypoints(vector<coll_obj> const &cobjs) {
		int const cc(camera_change);
		camera_change = 0; // messes up collision detection code

		for (vector<coll_obj>::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
			if (i->status != COLL_STATIC) continue; // only looking for static objects

			switch (i->type) {
			case COLL_CUBE: // can stand on the top
				add_waypoint_rect(i->d[0][0], i->d[1][0], i->d[0][1], i->d[1][1], i->d[2][1], i->id); // top rect
				break;

			case COLL_CYLINDER: // can stand on the top
				if (i->points[0].z > i->points[1].z) {
					add_waypoint_circle(i->points[0], i->radius,  i->id);
				}
				else {
					add_waypoint_circle(i->points[1], i->radius2, i->id);
				}
				break;

			case COLL_POLYGON:
				assert(i->npoints == 3 || i->npoints == 4); // triangle or quad
				if (i->thickness > MIN_POLY_THICK2) break; // extruded polygons not yet handled
				if (fabs(i->norm.z) < 0.5)          break; // need a mostly vertical polygon to stand on
				add_waypoint_triangle(i->points[0], i->points[1], i->points[2], i->id);
				if (i->npoints == 4) add_waypoint_triangle(i->points[0], i->points[2], i->points[3], i->id); // quad only
				break;

			case COLL_SPHERE:       break; // not supported (can't stand on)
			case COLL_CYLINDER_ROT: break; // not supported (can't stand on)
			default: assert(0);
			}
		}
		camera_change = cc;
	}

	void add_mesh_waypoints() {
		for (int y = 0; y < MESH_Y_SIZE; ++y) {
			for (int x = 0; x < MESH_X_SIZE; ++x) {
				if (is_mesh_disabled(x, y)) continue; // mesh disabled
				float zval(mesh_height[y][x]);
			
				if (!DISABLE_WATER && has_water(x, y) && zval < water_matrix[y][x]) { // underwater
					if (temperature <= W_FREEZE_POINT) { // ice
						zval = water_matrix[y][x]; // walk on ice
					}
					else { // water
						continue; // don't go underwater
					}
				}
				// *** WRITE ***
			}
		}
	}
};


void create_waypoints(vector<point> &waypoints) {

	RESET_TIME;
	waypoints.clear();
	waypoint_builder wb(waypoints);
	wb.add_cobj_waypoints(coll_objects);
	wb.add_mesh_waypoints();
	PRINT_TIME("  Waypoint Create");
	cout << "Waypoints: " << waypoints.size() << ", cobjs: " << coll_objects.size() << endl;
}





