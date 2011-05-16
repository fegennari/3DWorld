// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/13/02

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"
#include "player_state.h"


bool const SHOW_WAYPOINTS      = 1;
bool const SHOW_WAYPOINT_EDGES = 0;
int const WP_RESET_FRAMES      = 100; // Note: in frames, not ticks, fix?
int const WP_RECENT_FRAMES     = 200;

vector<waypoint_t> waypoints;

extern bool use_waypoints;
extern int DISABLE_WATER, camera_change, frame_counter, num_smileys;
extern float temperature, zmin, tfticks, water_plane_z;
extern obj_type object_types[];
extern dwobject def_objects[];
extern vector<coll_obj> coll_objects;


// ********** waypt_used_set **********


void waypt_used_set::clear() {

	used.clear();
	last_wp    = 0;
	last_frame = 0;
}

void waypt_used_set::insert(unsigned wp) { // called when a waypoint has been reached

	last_wp    = wp;
	last_frame = frame_counter;
	used[wp]   = frame_counter;
}

bool waypt_used_set::is_valid(unsigned wp) { // called to determine whether or not a waypoint is valid based on when it was last used

	if (last_frame > 0) {
		if ((frame_counter - last_frame) >= WP_RECENT_FRAMES) {
			clear(); // last_wp has expired, so all of the others must have expired as well
			return 1;
		}
		if (wp == last_wp) return 0; // too recent (lasts used)
	}
	map<unsigned, int>::iterator it(used.find(wp));
	if (it == used.end()) return 1; // new waypoint
	if ((frame_counter - it->second) < WP_RESET_FRAMES) return 0; // too recent
	used.erase(it); // lazy update - remove the waypoint when found to be expired
	return 1;
}


// ********** waypoint_t **********


waypoint_t::waypoint_t(point const &p=all_zeros, bool const up=0) : pos(p), user_placed(up), visited(0) {

	smiley_times.resize(num_smileys, tfticks);
}


void waypoint_t::mark_visited_by_smiley(unsigned const smiley_id) {

	assert(smiley_id < smiley_times.size());
	smiley_times[smiley_id] = tfticks;
	visited = 1;
}


float waypoint_t::get_time_since_last_visited(unsigned const smiley_id) const {

	assert(smiley_id < smiley_times.size());
	float const delta_t(tfticks - smiley_times[smiley_id]);
	assert(delta_t >= 0.0);
	return delta_t;
}


// ********** waypoint_builder **********


class waypoint_builder {

	float radius;

	bool check_cobj_placement(point &pos, int coll_id) const {
		dwobject obj(def_objects[WAYPOINT]); // create a fake temporary smiley object
		obj.pos     = pos;
		obj.coll_id = coll_id; // ignore collisions with the current object
		bool const ret(!obj.check_vert_collision(0, 0, 0)); // return true if no collision
		pos = obj.pos;
		return ret;
	}

	bool is_waypoint_valid(point const &pos, int coll_id) const {
		if (!is_over_mesh(pos) || pos.z < zmin || !point_interior_to_mesh(get_xpos(pos.x), get_ypos(pos.y))) return 0;
		float const mesh_zval(interpolate_mesh_zval(pos.x, pos.y, 0.0, 0, 0));
		if (pos.z - radius < mesh_zval) return 0; // bottom of smiley is under the mesh - should use a mesh waypoint here
		return check_cobj_placement(point(pos), coll_id);
	}

	void add_if_valid(point const &pos, int coll_id) {
		if (is_waypoint_valid(pos, coll_id)) waypoints.push_back(waypoint_t(pos, 0));
	}

	void add_waypoint_rect(float x1, float y1, float x2, float y2, float z, int coll_id) {
		if (min(x2-x1, y2-y1) < radius) return; // too small to stand on
		point const center(0.5*(x1+x2), 0.5*(y1+y2), z+radius);
		add_if_valid(center, coll_id);
		// FIXME: try more points?
	}

	void add_waypoint_circle(point const &p, float r, int coll_id) {
		if (r < radius) return; // too small to stand on
		add_if_valid(p + point(0.0, 0.0, radius), coll_id);
	}

	void add_waypoint_triangle(point const &p1, point const &p2, point const &p3, int coll_id) {
		if (dist_less_than(p1, p2, 2*radius) || dist_less_than(p2, p3, 2*radius) || dist_less_than(p3, p1, 2*radius)) return; // too small to stand on
		point const center(((p1 + p2 + p3) / 3.0) + point(0.0, 0.0, radius));
		add_if_valid(center, coll_id);
		// could try more points
	}

public:
	waypoint_builder(void) : radius(object_types[WAYPOINT].radius) {}
	
	void add_cobj_waypoints(vector<coll_obj> const &cobjs) {
		int const cc(camera_change);
		camera_change = 0; // messes up collision detection code
		unsigned const num_waypoints(waypoints.size());

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
		cout << "Added " << (waypoints.size() - num_waypoints) << " cobj waypoints" << endl;
		camera_change = cc;
	}

	void add_mesh_waypoints() {
		unsigned const mesh_skip_dist(16);
		unsigned const num_waypoints(waypoints.size());

		for (int yy = 0; yy <= MESH_Y_SIZE; yy += mesh_skip_dist) {
			for (int xx = 0; xx <= MESH_X_SIZE; xx += mesh_skip_dist) {
				int const x(max(1, min(MESH_X_SIZE-2, xx))), y(max(1, min(MESH_Y_SIZE-2, yy)));
				if (is_mesh_disabled(x, y)) continue; // mesh disabled
				float zval(mesh_height[y][x]);
			
				if (!DISABLE_WATER && has_water(x, y) && zval < water_plane_z) { // underwater
					if (temperature <= W_FREEZE_POINT) { // ice
						zval = water_matrix[y][x]; // walk on ice
					}
					else { // water
						continue; // don't go underwater
					}
				}
				// *** WRITE - more filtering ***
				point const pos(get_xval(x), get_yval(y), (zval + radius));
				add_if_valid(pos, -1);
			}
		}
		cout << "Added " << (waypoints.size() - num_waypoints) << " terrain waypoints" << endl;
	}

	void connect_waypoints() {
		unsigned cand_edges(0), num_edges(0), tot_steps(0);
		unsigned const num_waypoints(waypoints.size());
		float const step_size(0.25*radius), step_height(C_STEP_HEIGHT*radius);
		int cindex(-1);
		vector<pair<float, unsigned> > cands;

		for (unsigned i = 0; i < num_waypoints; ++i) {
			point const start(waypoints[i].pos);
			cands.resize(0);

			for (unsigned j = 0; j < num_waypoints; ++j) {
				point const end(waypoints[j].pos);

				if (cindex >= 0) {
					assert((unsigned)cindex < coll_objects.size());
					if (coll_objects[cindex].line_intersect(start, end)) continue; // hit last cobj
				}
				if (i == j || check_coll_line(start, end, cindex, -1, 1, 0)) continue;
				cands.push_back(make_pair(p2p_dist_sq(start, end), j));
			}
			sort(cands.begin(), cands.end()); // closest to furthest

			for (unsigned j = 0; j < cands.size(); ++j) {
				unsigned const k(cands[j].second);
				point const end(waypoints[k].pos);
				vector3d const dir(end - start), dir_xy(vector3d(dir.x, dir.y, 0.0).get_norm());
				float const dmag_inv(1.0/dir.xy_mag());
				vector<unsigned> const &next(waypoints[i].next_wpts);
				bool colinear(0);

				for (unsigned l = 0; l < next.size() && !colinear; ++l) {
					assert(next[l] < num_waypoints);
					vector3d const dir2(waypoints[next[l]].pos - start), dir_xy2(vector3d(dir2.x, dir2.y, 0.0).get_norm());
					colinear = (dot_product(dir_xy, dir_xy2) > 0.999);
				}
				if (colinear) continue;
				point cur(start); // first and last points are guaranteed to be valid
				float zvel(0.0);

				while (1) {
					if (dist_less_than(cur, end, 1.1*step_size)) {
						waypoints[i].next_wpts.push_back(k);
						++num_edges;
						break;
					}
					vector3d const delta((end - cur).get_norm());
					point lpos(cur);
					cur += delta*step_size;
					int const ret(set_true_obj_height(cur, lpos, C_STEP_HEIGHT, zvel, WAYPOINT, -2, 0, 0));
					if (ret == 3)                                           break; // stuck
					if ((cur.z - lpos.z) > C_STEP_HEIGHT*radius)            break; // too high of a step
					check_cobj_placement(cur, -1);
					if (dot_product_ptv(delta, cur, lpos) < 0.01*step_size) break; // not making progress
					float const d(fabs((end.x - start.x)*(start.y - cur.y) - (end.y - start.y)*(start.x - cur.x))*dmag_inv);
					if (d > 2.0*radius)                                     break; // path deviation too long
					// FIXME: allow more deviation from path to get around small obstacles
					++tot_steps;
					//waypoints.push_back(waypoint_t(cur, 1)); // testing
				}
				++cand_edges;
			}
		}
		cout << "cand edges: " << cand_edges << ", true edges: " << num_edges << ", tot steps: " << tot_steps << endl;
	}
};


// ********** waypoint top level code **********



void create_waypoints(vector<point> const &user_waypoints) {

	RESET_TIME;
	waypoints.clear();
	
	for (unsigned i = 0; i < user_waypoints.size(); ++i) {
		waypoints.push_back(waypoint_t(user_waypoints[i], 1));
	}
	waypoint_builder wb;

	if (use_waypoints) {
		wb.add_cobj_waypoints(coll_objects);
		wb.add_mesh_waypoints();
		PRINT_TIME("  Waypoint Generation");
	}
	wb.connect_waypoints();
	PRINT_TIME("  Waypoint Connectivity");
	cout << "Waypoints: " << waypoints.size() << endl;
}


void shift_waypoints(vector3d const &vd) {

	for (unsigned i = 0; i < waypoints.size(); ++i) {
		waypoints[i].pos += vd;
	}
}


void draw_waypoints() {

	if (!SHOW_WAYPOINTS) return;
	pt_line_drawer pld;

	for (vector<waypoint_t>::const_iterator i = waypoints.begin(); i != waypoints.end(); ++i) {
		unsigned const wix(i - waypoints.begin());
		set_color(i->visited ? ORANGE : (i->user_placed ? YELLOW : WHITE));
		draw_sphere_at(i->pos, 0.25*object_types[WAYPOINT].radius, N_SPHERE_DIV/2);
		if (!SHOW_WAYPOINT_EDGES) continue;

		for (vector<unsigned>::const_iterator j = i->next_wpts.begin(); j != i->next_wpts.end(); ++j) {
			assert(*j < waypoints.size());
			assert(*j != wix);
			waypoint_t const &w(waypoints[*j]);
			bool bidir(0);

			for (vector<unsigned>::const_iterator k = w.next_wpts.begin(); k != w.next_wpts.end() && !bidir; ++k) {
				bidir = (*k == wix);
			}
			if (!bidir || *j < wix) {
				pld.add_line(i->pos, plus_z, (bidir ? YELLOW : WHITE), w.pos, plus_z, (bidir ? YELLOW : ORANGE));
			}
		}
	}
	pld.draw();
}





