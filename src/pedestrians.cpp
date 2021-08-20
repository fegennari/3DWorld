// 3D World - Pedestrians for Procedural Cities
// by Frank Gennari
// 12/6/18
#include "city.h"
#include "shaders.h"

float const CROSS_SPEED_MULT = 1.8; // extra speed multiplier when crossing the road
float const CROSS_WAIT_TIME  = 60.0; // in seconds
bool const FORCE_USE_CROSSWALKS = 0; // more realistic and safe, but causes problems with pedestian collisions

extern bool tt_fire_button_down;
extern int display_mode, game_mode, animate2, frame_counter;
extern float FAR_CLIP;
extern double camera_zh;
extern point pre_smap_player_pos;
extern city_params_t city_params;


string gen_random_name(rand_gen_t &rgen); // from Universe_name.cpp
void get_closest_dim_dir_xy(cube_t const &inner, cube_t const &outer, bool &dim, bool &dir);

string pedestrian_t::get_name() const {
	rand_gen_t rgen;
	rgen.set_state(ssn, 123); // use ssn as name rand gen seed
	return gen_random_name(rgen); // for now, borrow the universe name generator to assign silly names
}
string pedestrian_t::str() const { // Note: no label_str()
	std::ostringstream oss;
	oss << get_name() << ": " << TXTn(ssn) << TXT(speed) << TXTn(radius) << TXT(city) << TXT(plot) << TXT(next_plot) << TXT(dest_plot) << TXTn(dest_bldg)
		<< TXTi(stuck_count) << TXT(collided) << TXTn(in_the_road) << TXT(is_stopped) << TXT(at_dest) << TXTn(has_dest_car) << TXT(target_valid())
		<< "wait_time=" << get_wait_time_secs(); // Note: pos, vel, dir not printed
	return oss.str();
}

float pedestrian_t::get_speed_mult() const {return (in_the_road ? CROSS_SPEED_MULT : 1.0);}

void pedestrian_t::stop() {
	//dir = vel.get_norm(); // ???
	vel = zero_vector;
	anim_time  = 0.0; // reset animation so that ped is standing normally and not mid-stride - should really transition this gradually somehow
	is_stopped = 1;
}
void pedestrian_t::go() {
	vel = dir * speed; // assumes dir is correct
	is_stopped = 0;
}
void pedestrian_t::wait_for(float seconds) {
	anim_time     = 0.0; // reset animation
	waiting_start = seconds*TICKS_PER_SECOND; // stop for N seconds
	target_pos    = all_zeros; // clear any previous target
}
cube_t pedestrian_t::get_bcube() const {
	cube_t c;
	c.set_from_sphere(pos, get_width());
	set_cube_zvals(c, get_z1(), get_z2());
	return c;
}

float get_sidewalk_width        () {return SIDEWALK_WIDTH*city_params.road_width;} // approx sidewalk width in the texture
float get_sidewalk_walkable_area() {return 0.65*get_sidewalk_width();} // walkable area of sidewalk on the street side; 65%, to avoid streetlights and traffic lights
float get_inner_sidewalk_width  () {return 1.00*get_sidewalk_width();} // walkable area of sidewalk on the plot side (not a sidewalk texture in residential neighborhoods)

bool pedestrian_t::check_inside_plot(ped_manager_t &ped_mgr, point const &prev_pos, cube_t const &plot_bcube, cube_t const &next_plot_bcube) {
	if (in_building) return 0; // not implemented yet
	//if (ssn == 2516) {cout << "in_the_road: " << in_the_road << ", pos: " << pos.str() << ", plot_bcube: " << plot_bcube.str() << ", npbc: " << next_plot_bcube.str() << endl;}
	if (plot_bcube.contains_pt_xy(pos)) {return 1;} // inside the plot
	stuck_count = 0; // no longer stuck
	if (next_plot == plot) return 0; // no next plot - clip to this plot
	
	if (next_plot_bcube.contains_pt_xy(pos)) {
		ped_mgr.move_ped_to_next_plot(*this);
		next_plot = ped_mgr.get_next_plot(*this); // FIXME: update next_plot_bcube?
		return 1;
	}
	cube_t union_plot_bcube(plot_bcube);
	union_plot_bcube.union_with_cube(next_plot_bcube);
	if (!union_plot_bcube.contains_pt_xy(pos) && union_plot_bcube.contains_pt_xy(prev_pos)) {return 0;} // went outside the valid area
	float const dx(min(fabs(pos.x - plot_bcube.x1()), fabs(pos.x - plot_bcube.x2()))), dy(min(fabs(pos.y - plot_bcube.y1()), fabs(pos.y - plot_bcube.y2())));

	if (max(dx, dy) < 0.75*city_params.road_width) { // near an intersection - near the road in both dims
		at_crosswalk = 1; // Note: should only be at crosswalks; but if we actually are corssing the road, this is the correct thing to do
	}
	in_the_road = 1;
	return 1; // allow peds to cross the road; don't need to check for building or other object collisions
}

bool pedestrian_t::check_road_coll(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube) const {
	if (!in_the_road) return 0;
	float const expand((get_sidewalk_width() - get_sidewalk_walkable_area()) + radius); // max dist from plot edge where a collision can occur
	cube_t pbce(plot_bcube), npbce(next_plot_bcube);
	pbce.expand_by_xy(expand);
	npbce.expand_by_xy(expand);
	if ((!pbce.contains_pt_xy(pos)) && (!npbce.contains_pt_xy(pos))) return 0; // ped is too far from the edge of the road to collide with streetlights or stoplights
	if (ped_mgr.check_isec_sphere_coll(*this)) return 1;
	if (ped_mgr.check_streetlight_sphere_coll(*this)) return 1;
	return 0;
}

bool pedestrian_t::is_valid_pos(vect_cube_t const &colliders, bool &ped_at_dest, ped_manager_t const *const ped_mgr) const {
	if (in_the_road || in_building) return 1; // not in a plot, no collision detection needed
	unsigned building_id(0);

	if (check_buildings_ped_coll(pos, radius, plot, building_id)) {
		if (!has_dest_bldg || building_id != dest_bldg) return 0;
		bool const ret(!at_dest);
		ped_at_dest = 1;
		return ret; // only valid if we just reached our dest
	}
	float const xmin(pos.x - radius), xmax(pos.x + radius);

	for (auto i = colliders.begin(); i != colliders.end(); ++i) {
		if (i->x2() < xmin) continue; // to the left
		if (i->x1() > xmax) break; // to the right - sorted from left to right, so no more colliders can intersect - done
		if (!sphere_cube_intersect(pos, radius, *i)) continue;
		if (!has_dest_car || !ped_mgr || plot != dest_plot) return 0; // not looking for car intersection
		
		if (i->intersects_xy(dest_car_center)) { // check if collider is a parking lot car group that contains the dest car
			// Note: here we consider a collision with any car in this block as at destination, even if it's not the dest car;
			// it's possible that the dest car is walled in and surrounded by cars (which may be poorly parked) such that the ped can't reach it without a collision
			if (!ped_mgr->has_car_at_pt(pos, city, 1)) continue; // no car at this location, continue into parking lot (thread safe and faster version)
			//if (ped_mgr->get_car_manager().get_car_at_pt(pos, 1) == nullptr) continue; // no car at this location, continue into parking lot (slow, but not called very often)
			bool const ret(!at_dest);
			ped_at_dest = 1;
			return ret; // only valid if we just reached our dest
		}
		return 0; // bad collision
	} // for i
	return 1;
}

void register_ped_coll(pedestrian_t &p1, pedestrian_t &p2, unsigned pid1, unsigned pid2) {
	p1.collided = p1.ped_coll = 1; p1.colliding_ped = pid2;
	p2.collided = p2.ped_coll = 1; p2.colliding_ped = pid1;
}

bool pedestrian_t::check_ped_ped_coll_range(vector<pedestrian_t> &peds, unsigned pid, unsigned ped_start, unsigned target_plot, float prox_radius, vector3d &force) {
	float const prox_radius_sq(prox_radius*prox_radius);

	for (auto i = peds.begin()+ped_start; i != peds.end(); ++i) { // check every ped until we exit target_plot
		if (i->plot != target_plot) break; // moved to a new plot, no collision, done; since plots are globally unique across cities, we don't need to check cities
		float const dist_sq(p2p_dist_xy_sq(pos, i->pos));
		if (dist_sq > prox_radius_sq) continue; // proximity test
		if (i->destroyed) continue; // dead
		float const r_sum(0.6f*(radius + i->radius)); // using a smaller radius to allow peds to get close to each other
		if (dist_sq < r_sum*r_sum) {register_ped_coll(*this, *i, pid, (i - peds.begin())); return 1;} // collision
		if (speed < TOLERANCE) continue;
		vector3d const delta_v(vel - i->vel), delta_p((pos.x - i->pos.x), (pos.y - i->pos.y), 0.0);
		float const dp(-dot_product_xy(delta_v, delta_p));
		if (dp <= 0.0) continue; // diverging, no avoidance needed
		float const dv_mag(delta_v.mag()), dist(sqrt(dist_sq)), fmag(dist/(dist - 0.9*r_sum));
		if (dv_mag < TOLERANCE) continue;
		vector3d const rejection(delta_p - (dp/(dv_mag*dv_mag))*delta_v); // component of velocity perpendicular to delta_p (avoid dir)
		float const rmag(rejection.mag()), rel_vel(max(dv_mag/speed, 0.5f)); // higher when peds are converging
		if (rmag < TOLERANCE) continue;
		float const force_mult(dp/(dv_mag*dist)); // stronger with head-on collisions
		force += rejection*(rel_vel*force_mult*fmag/rmag);
		//cout << TXT(r_sum) << TXT(dist) << TXT(fmag) << ", dv: " << delta_v.str() << ", dp: " << delta_p.str() << ", rej: " << rejection.str() << ", force: " << force.str() << endl;
	} // for i
	return 0;
}

bool pedestrian_t::check_ped_ped_coll(ped_manager_t const &ped_mgr, vector<pedestrian_t> &peds, unsigned pid, float delta_dir) {
	if (in_building) return 0; // no ped-ped collisions in buildings (yet)
	assert(pid < peds.size());
	float const timestep(2.0*TICKS_PER_SECOND), lookahead_dist(timestep*speed); // how far we can travel in 2s
	float const prox_radius(1.2*radius + lookahead_dist); // assume other ped has a similar radius
	vector3d force(zero_vector);
	if (check_ped_ped_coll_range(peds, pid, pid+1, plot, prox_radius, force)) return 1;

	if (in_the_road && next_plot != plot) {
		// need to check for coll between two peds crossing the street from different sides, since they won't be in the same plot while in the street
		unsigned const ped_ix(ped_mgr.get_first_ped_at_plot(next_plot));
		assert(ped_ix <= peds.size()); // could be at the end
		if (check_ped_ped_coll_range(peds, pid, ped_ix, next_plot, prox_radius, force)) return 1;
	}
	if (force != zero_vector) {set_velocity((0.1*delta_dir)*force + ((1.0 - delta_dir)/speed)*vel);} // apply ped repulsive force
	return 0;
}

bool pedestrian_t::check_ped_ped_coll_stopped(vector<pedestrian_t> &peds, unsigned pid) {
	if (in_building) return 0; // no ped-ped collisions in buildings (yet)
	assert(pid < peds.size());

	// Note: shouldn't have to check peds in the next plot, assuming that if we're stopped, they likely are as well, and won't be walking toward us
	for (auto i = peds.begin()+pid+1; i != peds.end(); ++i) { // check every ped until we exit target_plot
		if (i->plot != plot) break; // moved to a new plot, no collision, done; since plots are globally unique across cities, we don't need to check cities
		if (!dist_xy_less_than(pos, i->pos, 0.6f*(radius + i->radius))) continue; // no collision
		if (i->destroyed) continue; // dead
		i->collided = i->ped_coll = 1; i->colliding_ped = pid;
		return 1; // Note: could omit this return and continue processing peds
	} // for i
	return 0;
}

point rand_xy_pt_on_cube_edge(cube_t const &c, float radius, rand_gen_t &rgen) {
	bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
	point pt;
	pt[ dim] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*radius;
	pt[!dim] = rgen.rand_uniform(c.d[!dim][0]+radius, c.d[!dim][1]-radius);
	pt.z     = c.z1();
	return pt;
}
bool pedestrian_t::try_place_in_plot(cube_t const &plot_cube, vect_cube_t const &colliders, unsigned plot_id, rand_gen_t &rgen) {
	pos    = rand_xy_pt_on_cube_edge(plot_cube, radius, rgen);
	pos.z += radius; // place on top of the plot
	plot   = next_plot = dest_plot = plot_id; // set next_plot and dest_plot as well so that they're valid for the first frame
	bool temp_at_dest(0); // we don't want to set at_dest from this call
	if (!is_valid_pos(colliders, temp_at_dest, nullptr)) return 0; // plot == next_plot; return if failed
	return 1; // success
}

float path_finder_t::path_t::calc_length_up_to(const_iterator i) const {
	assert(i <= end());
	float len(0.0);
	for (auto p = begin(); p+1 != i; ++p) {len += p2p_dist(*p, *(p+1));}
	return len;
}
cube_t path_finder_t::path_t::calc_bcube() const { // Note: could probably get away with only x/y bounds
	assert(!empty());
	cube_t bcube(front(), front());
	for (auto i = begin()+1; i != end(); ++i) {bcube.union_with_pt(*i);}
	return bcube;
}

// path_finder_t
bool path_finder_t::add_pt_to_path(point const &p, path_t &path) const {
	if (!plot_bcube.contains_pt_xy(p)) return 0; // outside the plot - invalid
	path.push_back(p);
	return 1;
}

bool path_finder_t::add_pts_around_cube_xy(path_t &path, path_t const &cur_path, path_t::const_iterator p, cube_t const &c, bool dir) {
	point const &n(*(p+1));
	if (c.contains_pt_xy(*p) || c.contains_pt_xy(n)) return 0; // point contained in cube - likely two overlapping/adjacent building cubes - fail
	path.clear();
	path.insert(path.end(), cur_path.begin(), p+1); // add the prefix, including p
	cube_t ec(c);
	ec.expand_by_xy(gap); // expand cubes in x and y
	point const corners [4] = {point( c.x1(),  c.y1(), p->z), point( c.x1(),  c.y2(), p->z), point( c.x2(),  c.y2(), p->z), point( c.x2(),  c.y1(), p->z)}; // CW
	point const ecorners[4] = {point(ec.x1(), ec.y1(), p->z), point(ec.x1(), ec.y2(), p->z), point(ec.x2(), ec.y2(), p->z), point(ec.x2(), ec.y1(), p->z)}; // CW; expanded cube
	vector3d const delta(n - *p); // not normalized
	// find the two closest corners to the left and right
	float min_dp1(0.0), min_dp2(0.0);
	unsigned cix1(4), cix2(4); // start at invalid values
	//cout << TXT(p->x) << TXT(p->y) << TXT(n.x) << TXT(n.y) << TXT(delta.x) << TXT(delta.y) << TXT(c.x1()) << TXT(c.x2()) << TXT(c.y1()) << TXT(c.y2()) << endl;

	for (unsigned i = 0; i < 4; ++i) {
		vector3d const delta2(corners[i] - *p); // unexpanded corners
		float const dp(dot_product_xy(delta, delta2)/delta2.mag());
		bool const turn_dir(cross_product_xy(delta, delta2) < 0.0);
		//cout << TXT(delta2.x) << TXT(delta2.y) << TXT(i) << TXT(dp) << TXT(turn_dir) << endl;
		if ((turn_dir ? cix1 : cix2) == 4 || dp < (turn_dir ? min_dp1 : min_dp2)) {(turn_dir ? min_dp1 : min_dp2) = dp; (turn_dir ? cix1 : cix2) = i;}
	}
	//cout << TXT(dir) << TXT(min_dp1) << TXT(min_dp2) << TXT(cix1) << TXT(cix2) << endl;
	if (cix1 == 4 || cix2 == 4) return 0; // something bad happened (floating-point error?), fail
	unsigned dest_cix(dir ? cix2 : cix1);
	bool const move_dir((((dest_cix+1)&3) == (dir ? cix1 : cix2)) ? 0 : 1); // CCW/CW based on which dir moves around the other side of the cube
	if (check_line_clip_xy(*p, ecorners[dest_cix], c.d)) return 0; // something bad happened (floating-point error?), fail
	//if (!line_int_cubes_xy(*p, ecorners[dest_cix], avoid)) return 0; // TODO: should we test other avoid cubes that may be blocking the path?
	if (!add_pt_to_path(ecorners[dest_cix], path)) return 0; // expanded corner

	if (check_line_clip_xy(n, ecorners[dest_cix], c.d)) { // no path to dest, add another point
		if (move_dir) {dest_cix = (dest_cix+1)&3;} else {dest_cix = (dest_cix+3)&3;}
		assert(dest_cix != (dir ? cix1 : cix2));
		if (!add_pt_to_path(ecorners[dest_cix], path)) return 0; // expanded corner

		if (check_line_clip_xy(n, ecorners[dest_cix], c.d)) { // no path to dest, add another point
			if (move_dir) {dest_cix = (dest_cix+1)&3;} else {dest_cix = (dest_cix+3)&3;}
			if (!add_pt_to_path(ecorners[dest_cix], path)) return 0; // expanded corner
			//assert(!check_line_clip_xy(n, ecorners[dest_cix], c.d)); // must have a path now
		}
	}
	path.insert(path.end(), p+1, cur_path.end()); // add the suffix, including p+1
	path.calc_length();
	return 1;
}

void path_finder_t::find_best_path_recur(path_t const &cur_path, unsigned depth) {
	if (depth >= MAX_PATH_DEPTH) return; // depth is too high, fail (stack not allocated)
	if (cur_path.length >= best_path.length) return; // bound (best path length is set to an upper bound even when not valid)
	cube_t const bcube(cur_path.calc_bcube());
	path_t::const_iterator first_int_p(cur_path.end());
	float tmin(1.0);
	unsigned cix(0);

	for (unsigned ix = 0; ix < avoid.size(); ++ix) {
		if (used[ix]) continue; // done with this cube
		cube_t const &c(avoid[ix]);
		if (!c.intersects_xy(bcube)) continue;

		for (auto p = cur_path.begin(); p+1 != cur_path.end() && p <= first_int_p; ++p) { // iterate over line segments up to/including first_int_p, skip last point
			float c_tmin, c_tmax;
			if (!get_line_clip_xy(*p, *(p+1), c.d, c_tmin, c_tmax)) continue;
			if (p < first_int_p || c_tmin < tmin) {first_int_p = p; tmin = c_tmin; cix = ix;} // intersection
		}
	} // for ix
	if (first_int_p != cur_path.end()) {
		assert(tmin < 1.0);
		path_t &next_path(path_stack[depth]);
		assert(!used[cix]);
		used[cix] = 1; // mark this cube used so that we don't try to intersect it again (and to avoid floating-point errors with line adjacency)

		for (unsigned d = 0; d < 2; ++d) {
			if (add_pts_around_cube_xy(next_path, cur_path, first_int_p, avoid[cix], d)) {find_best_path_recur(next_path, depth+1);} // recursive call
		}
		assert(used[cix]);
		used[cix] = 0; // mark cube as unused
		if (first_int_p == cur_path.begin() || found_complete_path()) return; // path is no good, terminate
		// calculate the length of the partial path; add twice the distance we're short (to the destination) as a penalty
		float const partial_len(cur_path.calc_length_up_to(first_int_p+1) + 2.0*p2p_dist(*first_int_p, dest));
		if (partial_len >= partial_path.length) return; // not better
		partial_path.clear();
		partial_path.insert(partial_path.end(), cur_path.begin(), first_int_p+1); // record best partial path seen
		partial_path.length = partial_len;
		return;
	}
	if (cur_path.length < best_path.length) { // this test almost always succeeds
		best_path = cur_path; // if we got here without returning above, this is the best path seen so far
		partial_path.clear(); // not using partial_path after this point
		partial_path.length = 0.0;
	}
}

bool path_finder_t::shorten_path(path_t &path) const {
	if (path.size() <= 2) return 0; // nothing to do
	path_t::iterator i(path.begin()+1), o(i); // skip first point, which we must keep

	for (; i != path.end(); ++i) {
		if (i+1 == path.end()) {*(o++) = *i; continue;} // always keep last point
		point const &a(*(i-1)), &b(*(i+1));
		if (line_int_cubes_xy(a, b, avoid)) {*(o++) = *i;} // keep if removing this point generates an intersection
	}
	if (o == path.end()) return 0; // no update
	path.erase(o, path.end());
	path.calc_length();
	return 1; // shortened
}

bool path_finder_t::find_best_path() {
	used.clear();
	used.resize(avoid.size(), 0);
	best_path.clear();
	partial_path.clear();
	cur_path.clear();
	cur_path.init(pos, dest);
	best_path.length = partial_path.length = 5.0*cur_path.length; // add an upper bound of 4x length to avoid too much recursion
	find_best_path_recur(cur_path, 0); // depth=0
	shorten_path(best_path); // see if we can remove any path points; this rarely has a big effect on path length, so it's okay to save time by doing this after the length test
	shorten_path(partial_path);
	//cout << TXT(avoid.size()) << TXT(cur_path.length) << TXT(best_path.length) << found_path() << endl;
	return found_path();
}

// Note: avoid must be non-overlapping and should be non-adjacent; even better if cubes are separated enough that peds can pass between them (> 2*ped radius)
// return values: 0=failed, 1=valid path, 2=init contained, 3=straight path (no collisions)
unsigned path_finder_t::run(point const &pos_, point const &dest_, cube_t const &plot_bcube_, float gap_, point &new_dest) {
	if (!line_int_cubes_xy(pos_, dest_, avoid)) return 3; // no work to be done, leave dest as it is
	pos = pos_; dest = dest_; plot_bcube = plot_bcube_; gap = gap_;
	//if (any_cube_contains_pt_xy(avoid, dest)) return 0; // invalid dest pos - ignore for now and let path finding deal with it when we get to that pos
	unsigned next_pt_ix(1); // default: point after pos

	if (!plot_bcube.contains_pt_xy(pos)) { // keep pos inside the plot
		plot_bcube.clamp_pt_xy(pos); // clamp point to plot bounds
		next_pt_ix = 0; // start at the new pos
	}
	else { // if there are any other cubes containing pos, move away from this cube; could be initial ped positions, ped pushed by a collision, or some other problem
		for (auto i = avoid.begin(); i != avoid.end(); ++i) {
			if (!i->contains_pt_xy(pos)) continue;
			int const building_id(get_building_bcube_contains_pos(pos));
			point new_pos;
			float dmin(0.0);

			for (unsigned dim = 0; dim < 2; ++dim) {
				for (unsigned dir = 0; dir < 2; ++dir) {
					float const edge(i->d[dim][dir]), dist(abs(pos[dim] - edge));
					if (dmin > 0.0 && dist >= dmin) continue; // not a better direction
					point cand_pos(pos);
					cand_pos[dim] = i->d[dim][dir];
					if (building_id >= 0 && check_line_coll_building(pos, cand_pos, building_id)) continue; // this dir intersects the building
					new_pos = cand_pos;
					dmin    = dist; // valid direction
				} // for dir
			} // for dim
			if (dmin > 0.0) { // new_pos is valid
				pos = new_pos;
				next_pt_ix = 0; // start at the new pos
				break; // at most one cube should contain pos
			}
		} // for i
	}
	if (!find_best_path()) return 0; // if we fail to find a path, leave new_dest unchanged
	vector<point> const &path(get_best_path());
	assert(next_pt_ix < path.size());
	new_dest = path[next_pt_ix]; // set dest to next point on the best path
	return (next_pt_ix ? 1 : 2); // return 2 for the init contained case
}

// pedestrian_t
point pedestrian_t::get_dest_pos(cube_t const &plot_bcube, cube_t const &next_plot_bcube, ped_manager_t const &ped_mgr) const {
	if (is_stopped && target_valid()) {return target_pos;} // stay the course (this case only needed for debug drawing)

	if (plot == dest_plot) { // this plot contains our dest building/car
		if (!at_dest && has_dest_bldg) { // not there yet
			cube_t const dest_bcube(get_building_bcube(dest_bldg));
			//if (dest_bcube.contains_pt_xy(pos)) {at_dest = 1;} // could set this here, but requiring a collision also works
			point const dest_pos(dest_bcube.get_cube_center()); // slowly adjust dir to move toward dest_bldg
			return point(dest_pos.x, dest_pos.y, pos.z); // same zval
		}
		else if (!at_dest && has_dest_car) {return point(dest_car_center.x, dest_car_center.y, pos.z);} // same zval
	}
	else if (next_plot != plot) { // move toward next plot
		if (!next_plot_bcube.contains_pt_xy(pos)) { // not yet crossed into the next plot
			bool const in_cur_plot(plot_bcube.contains_pt_xy(pos));
			point dest_pos(pos);
			// should cross at intersection, not in the middle of the street; find closest crosswalk (corner of plot_bcube) to dest_pos;
			// while the code below is correct, it tends to make peds collide with the stoplight on the corner and each other and never actually reach their destinations,
			// so we allow peds to cross the street wherever they want until at the very least the stoplights can be moved
			if (FORCE_USE_CROSSWALKS) { // closest corner (crosswalk)
				cube_t const &cube(in_cur_plot ? plot_bcube : next_plot_bcube); // target the corner of the current plot, then the corner of the next plot
				float const val((in_cur_plot ? 0.01 : -0.01)*city_params.road_width); // slightly outside the cur plot / inside the next plot, to ensure a proper transition

				for (unsigned d = 0; d < 2; ++d) { // x,y
					dest_pos[d] = (((pos[d] - cube.d[d][0]) < (cube.d[d][1] - pos[d])) ? (cube.d[d][0] - val) : (cube.d[d][1] + val));
				}
			}
			else { // closest point
				dest_pos = next_plot_bcube.closest_pt(pos);
			}
			if (!in_cur_plot) { // went outside the current plot
				cube_t union_plot_bcube(plot_bcube);
				union_plot_bcube.union_with_cube(next_plot_bcube);
				if (!union_plot_bcube.contains_pt_xy(pos)) {dest_pos = plot_bcube.closest_pt(pos);} // went outside on the wrong side, go back inside the current plot
			}
			dest_pos.z = pos.z; // same zval
			return dest_pos;
		}
	}
	return pos; // no dest
}

bool pedestrian_t::choose_alt_next_plot(ped_manager_t const &ped_mgr) {
	reset_waiting(); // reset waiting state regardless of outcome; we don't want to get here every frame if we fail to find another plot
	if (plot == next_plot) return 0; // no next plot (error?)
	//if (next_plot == dest_plot) return 0; // the next plot is our desination, should we still choose another plot?
	unsigned const cand_next_plot(ped_mgr.get_next_plot(*this, next_plot));
	if (cand_next_plot == next_plot || cand_next_plot == plot) return 0; // failed
	next_plot = cand_next_plot;
	return 1; // success
}

void add_and_expand_ped_avoid_cube(cube_t const &c, vect_cube_t &avoid, float expand, float height) {
	avoid.push_back(c);
	avoid.back().expand_by_xy(expand*((c.dz() < 0.67*height) ? 0.5 : 1.0)); // reduce expand value for short objects that will only collide with our legs
}

void pedestrian_t::get_avoid_cubes(ped_manager_t const &ped_mgr, vect_cube_t const &colliders,
	cube_t const &plot_bcube, cube_t const &next_plot_bcube, point &dest_pos, vect_cube_t &avoid) const
{
	avoid.clear();
	if (in_building) return; // not yet implemented, but if it was we would get the nearby building walls, objects, etc.
	float const height(get_height()), expand(1.1*radius); // slightly larger than radius to leave some room for floating-point error
	road_plot_t const &cur_plot(ped_mgr.get_city_plot_for_peds(city, plot));

	if (cur_plot.is_residential && !cur_plot.is_park) { // apply special restrictions when walking through a residential block
		cube_t avoid_area(plot_bcube);
		avoid_area.expand_by_xy(0.5*radius - (get_inner_sidewalk_width() + get_sidewalk_walkable_area())); // shrink to plot interior, and undo the expand applied to the plot
		bool avoid_entire_plot(0);

		if (plot == dest_plot && plot_bcube == next_plot_bcube) { // plot contains our destination, and the plot bcube has been updated
			if (city_params.assign_house_plots && (has_dest_bldg || has_dest_car)) { // we can only walk through our own sub-plot
				cube_t dest_cube;
				if      (has_dest_bldg) {dest_cube = get_building_bcube(dest_bldg);}
				else if (has_dest_car ) {dest_cube.set_from_sphere(dest_car_center, city_params.get_nom_car_size().x);} // somewhat approximate/conservative
				assert(dest_cube.intersects_xy(plot_bcube)); // or contains, or is that too strong?
				bool dim(0), dir(0);
				get_closest_dim_dir_xy(dest_cube, plot_bcube, dim, dir);
				cube_t approach_area(dest_cube);
				approach_area.d[dim][dir] = plot_bcube.d[dim][dir]; // expand out to the plot
				approach_area.expand_by_xy(radius); // add a small fudge factor

				if (!approach_area.contains_pt(pos)) { // not in approach area - walk around the plot
					// update dest_pos to use the proxy point along the sidewalk across from our destination as the next path point
					dest_pos[ dim]    = plot_bcube.d[dim][dir];
					dest_pos[!dim]    = dest_cube.get_center_dim(!dim);
					avoid_entire_plot = 1;
				}
			} // else we can walk through this plot
		}
		else {avoid_entire_plot = 1;} // not our destination plot, we can't walk through any residential properties

		if (avoid_entire_plot) {
			avoid.push_back(avoid_area); // this is the highest priority

			for (auto i = colliders.begin(); i != colliders.end(); ++i) { // remove any cubes contained in the plot, since they're redundant
				if (!avoid_area.contains_cube_xy(*i)) {add_and_expand_ped_avoid_cube(*i, avoid, expand, height);}
			}
			return; // done
		}
	} // else we can walk through this plot
	get_building_bcubes(cur_plot, avoid);
	expand_cubes_by_xy(avoid, expand); // expand building cubes in x and y to approximate a cylinder collision (conservative)
	//remove_cube_if_contains_pt_xy(avoid, pos); // init coll cases (for example from previous dest_bldg) are handled by path_finder_t
	if (plot == dest_plot && has_dest_bldg) {remove_cube_if_contains_pt_xy(avoid, dest_pos);} // exclude our dest building, we do want to collide with it

	for (auto i = colliders.begin(); i != colliders.end(); ++i) { // remove any cubes contained in the plot, since they're redundant
		if (plot == dest_plot && has_dest_car && i->contains_pt_xy(dest_pos)) continue; // exclude our dest car, we do want to collide with it
		add_and_expand_ped_avoid_cube(*i, avoid, expand, height);
	}
}

bool pedestrian_t::check_for_safe_road_crossing(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, vect_cube_t *dbg_cubes) const {
	if (!in_the_road || speed < TOLERANCE) return 1;
	float const sw_width(get_sidewalk_width());
	if (!plot_bcube.closest_dist_xy_less_than(pos, sw_width)) return 1; // too far into the road to turn back (should this use get_sidewalk_walkable_area()?)
	cube_t union_plot_bcube(plot_bcube);
	union_plot_bcube.union_with_cube(next_plot_bcube); // this is the area the ped is constrained to (both plots + road in between)
	if (!union_plot_bcube.contains_pt_xy(pos)) return 1; // not crossing between plots - must be in the road, go back to the sidewalk
	// just exited the plot and about the cross the road - check for cars; use speed rather than vel in case we're already stopped and vel==zero_vector
	float const dx(min((pos.x - plot_bcube.x1()), (plot_bcube.x2() - pos.x))), dy(min((pos.y - plot_bcube.y1()), (plot_bcube.y2() - pos.y)));
	bool const road_dim(dx < dy); // if at crosswalk, need to know which direction/road the ped is crossing
	float const time_to_cross((city_params.road_width - 2.0f*sw_width)/(speed*get_speed_mult())); // road area where cars can drive excluding sidewalks on each side
	//cout << "plot_bcube: " << plot_bcube.str() << " " << TXT(dx) << TXT(dy) << TXT(road_dim) << TXT(time_to_cross) << endl;
	return !ped_mgr.has_nearby_car(*this, road_dim, time_to_cross, dbg_cubes);
}

void pedestrian_t::move(ped_manager_t const &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, float &delta_dir) {
	if (!in_building) { // in the city
		if (!check_for_safe_road_crossing(ped_mgr, plot_bcube, next_plot_bcube)) {stop(); return;}
	}
	reset_waiting();
	if (is_stopped) {go();}

	if (target_valid()) { // if facing away from the target, rotate in place rather than moving in a circle
		vector3d const delta(target_pos - pos);
		float const dist(delta.mag());
		if (dist > radius && dot_product_xy(vel, delta) < 0.01*speed*dist) {delta_dir = min(1.0f, 4.0f*delta_dir); return;} // rotate faster
	}
	float const timestep(fticks*get_speed_mult());
	pos       += timestep*vel;
	anim_time += timestep*speed;
}

void pedestrian_t::run_path_finding(ped_manager_t &ped_mgr, cube_t const &plot_bcube, cube_t const &next_plot_bcube, vect_cube_t const &colliders, vector3d &dest_pos) {
	vect_cube_t &avoid(ped_mgr.path_finder.get_avoid_vector());
	get_avoid_cubes(ped_mgr, colliders, plot_bcube, next_plot_bcube, dest_pos, avoid);
	target_pos = all_zeros;
	cube_t union_plot_bcube(plot_bcube);
	union_plot_bcube.union_with_cube(next_plot_bcube); // this is the area the ped is constrained to (both plots + road in between)
	// run path finding between pos and dest_pos using avoid cubes
	if (ped_mgr.path_finder.run(pos, dest_pos, union_plot_bcube, 0.1*radius, dest_pos)) {target_pos = dest_pos;}
}

void pedestrian_t::get_plot_bcubes_inc_sidewalks(ped_manager_t const &ped_mgr, cube_t &plot_bcube, cube_t &next_plot_bcube) const {
	// this approach is more visually pleasing because pedestrians will actually walk on the edges of the roads on what appears to be the sidewalks;
	// unfortunately, they also run into streetlights, traffic lights, and each other in this narrow area
	plot_bcube      = ped_mgr.get_city_plot_for_peds(city, plot);
	next_plot_bcube = ped_mgr.get_city_plot_for_peds(city, next_plot);
	float const sidewalk_width(get_sidewalk_walkable_area());
	plot_bcube.expand_by_xy(sidewalk_width);
	next_plot_bcube.expand_by_xy(sidewalk_width);
}

void pedestrian_t::next_frame(ped_manager_t &ped_mgr, vector<pedestrian_t> &peds, unsigned pid, rand_gen_t &rgen, float delta_dir) {
	if (destroyed)    return; // destroyed
	if (speed == 0.0) return; // not moving, no update needed
	if (in_building)  return; // building update/movement logic handled elsewhere

	// navigation with destination
	if (at_dest) {
		register_at_dest();
		ped_mgr.choose_new_ped_plot_pos(*this);
	}
	else if (!has_dest_bldg && !has_dest_car) {ped_mgr.choose_dest_building_or_parked_car(*this);}
	if (at_crosswalk) {ped_mgr.mark_crosswalk_in_use(*this);}
	cube_t plot_bcube, next_plot_bcube;
	get_plot_bcubes_inc_sidewalks(ped_mgr, plot_bcube, next_plot_bcube);
	// movement logic
	point const prev_pos(pos); // assume this ped starts out not colliding
	move(ped_mgr, plot_bcube, next_plot_bcube, delta_dir);

	if (is_stopped) { // ignore any collisions and just stand there, keeping the same target_pos; will go when path is clear
		if (get_wait_time_secs() > CROSS_WAIT_TIME && choose_alt_next_plot(ped_mgr)) { // give up and choose another destination if waiting for too long
			target_pos = all_zeros;
			go(); // back up or turn so that we don't walk forward into the street? move() should attempt to rotate in place
		}
		else {
			check_ped_ped_coll_stopped(peds, pid); // still need to check for other peds colliding with us; this doesn't always work
			collided = ped_coll = 0;
			return;
		}
	}
	at_crosswalk = in_the_road = 0; // reset state for next frame; these may be set back to 1 below
	vect_cube_t const &colliders(ped_mgr.get_colliders_for_plot(city, plot));
	bool outside_plot(0);

	if (collided) {} // already collided with a previous ped this frame, handled below
	else if (!check_inside_plot(ped_mgr, prev_pos, plot_bcube, next_plot_bcube)) {collided = outside_plot = 1;} // outside the plot, treat as a collision with the plot bounds
	else if (!is_valid_pos(colliders, at_dest, &ped_mgr)) {collided = 1;} // collided with a static collider
	else if (check_road_coll(ped_mgr, plot_bcube, next_plot_bcube)) {collided = 1;} // collided with something in the road (stoplight, streetlight, etc.)
	else if (check_ped_ped_coll(ped_mgr, peds, pid, delta_dir)) {collided = 1;} // collided with another pedestrian
	else { // no collisions
		//cout << TXT(pid) << TXT(plot) << TXT(dest_plot) << TXT(next_plot) << TXT(at_dest) << TXT(delta_dir) << TXT((unsigned)stuck_count) << TXT(collided) << endl;
		vector3d dest_pos(get_dest_pos(plot_bcube, next_plot_bcube, ped_mgr));

		if (dest_pos != pos) {
			bool update_path(0);
			if (dist_less_than(pos, get_camera_pos(), 1000.0*radius)) { // nearby pedestrian - higher update rate
				update_path = (((frame_counter + ssn) & 15) == 0 || (target_valid() && dist_xy_less_than(pos, target_pos, radius)));
			}
			else { // distant pedestrian - lower update rate
				update_path = (((frame_counter + ssn) & 63) == 0);
			}
			// run only every several frames to reduce runtime; also run when at dest and when close to the current target pos or at the destination
			if (at_dest || update_path) {run_path_finding(ped_mgr, plot_bcube, next_plot_bcube, colliders, dest_pos);}
			else if (target_valid()) {dest_pos = target_pos;} // use previous frame's dest if valid
			vector3d dest_dir((dest_pos.x - pos.x), (dest_pos.y - pos.y), 0.0); // zval=0, not normalized
			float const dmag(dest_dir.xy_mag());

			if (speed > TOLERANCE && dmag > TOLERANCE) { // avoid divide-by-zero
				dest_dir /= dmag; // normalize
				// if destination is in exactly the opposite dir, pick an orthogonal direction using our SSN to decide which way deterministically
				if (dot_product_xy(dest_dir, vel) < -0.99*speed) {dest_dir = cross_product(vel, plus_z*((ssn&1) ? 1.0 : -1.0)).get_norm();}
				set_velocity((0.1*delta_dir)*dest_dir + ((1.0 - delta_dir)/speed)*vel); // slowly blend in destination dir (to avoid sharp direction changes)
			}
		}
		stuck_count = 0;
	}
	if (collided) { // collision
		if (!outside_plot) {
			point const cur_pos(pos);
			pos = prev_pos; // restore to previous valid pos unless we're outside the plot
			// if prev pos is also invalid, undo the restore to avoid getting this ped stuck in a collision object
			if (!is_valid_pos(colliders, at_dest, &ped_mgr) || check_road_coll(ped_mgr, plot_bcube, next_plot_bcube)) {pos = cur_pos;}
		}
		vector3d new_dir;

		if (++stuck_count > 8) {
			if (target_valid()) {pos += (0.1*radius)*(target_pos - pos).get_norm();} // move toward target_pos if it's valid since this should be a good direction
			else if (stuck_count > 100) {pos += (0.1*radius)*(get_dest_pos(plot_bcube, next_plot_bcube, ped_mgr) - pos).get_norm();} // move toward dest if stuck count is high
			else {pos += rgen.signed_rand_vector_spherical_xy()*(0.1*radius); } // shift randomly by 10% radius to get unstuck
		}
		if (ped_coll) {
			assert(colliding_ped < peds.size());
			vector3d const coll_dir(peds[colliding_ped].pos - pos);
			new_dir = cross_product(vel, plus_z);
			if (dot_product_xy(new_dir, coll_dir) > 0.0) {new_dir = -new_dir;} // orient away from the other ped
		}
		else { // static object collision (should be rare if path_finder does a good job)
			new_dir = rgen.signed_rand_vector_spherical_xy(); // try a random new direction
			if (dot_product_xy(vel, new_dir) > 0.0) {new_dir *= -1.0;} // negate if pointing in the same dir
		}
		set_velocity(new_dir);
		target_pos = all_zeros; // reset and force path finding to re-route from this new direction/pos
	}
	if (vel != zero_vector) { // if stopped, don't update dir
		if (!collided && target_valid()) {delta_dir = min(1.0f, 4.0f*delta_dir);} // use a tighter turning radius when there's an unobstructed target_pos
		dir = (delta_dir/speed)*vel + (1.0 - delta_dir)*dir; // merge velocity into dir gradually for smooth turning
		dir.z = 0.0; // should be zero, but set just in case
		dir.normalize();
	}
	collided = ped_coll = 0; // reset for next frame
}

void pedestrian_t::register_at_dest() {
	assert(plot == dest_plot);
	//cout << get_name() << " at destination " << (has_dest_car ? "car " : (has_dest_bldg ? "building " : "")) << dest_bldg << " in plot " << dest_plot << endl; // placeholder
}

bool pedestrian_t::is_close_to_player() const { // for debug printouts, etc.
	return dist_less_than((pos + get_tiled_terrain_model_xlate()), get_camera_pos(), 4.0*radius);
}


unsigned ped_model_loader_t::num_models() const {return city_params.ped_model_files.size();}

city_model_t const &ped_model_loader_t::get_model(unsigned id) const {
	assert(id < num_models());
	return city_params.ped_model_files[id];
}
city_model_t &ped_model_loader_t::get_model(unsigned id) {
	assert(id < num_models());
	return city_params.ped_model_files[id];
}

void ped_city_vect_t::add_ped(pedestrian_t const &ped, unsigned road_ix) {
	if (ped.city >= peds.size()) {peds.resize(ped.city+1);} // allocate city if needed
	auto &city(peds[ped.city]);
	if (road_ix  >= city.size()) {city.resize(road_ix+1 );} // allocate road if needed
	city[road_ix].emplace_back(ped.pos, ped.radius);
}
void ped_city_vect_t::clear() {
	for (auto i = peds.begin(); i != peds.end(); ++i) {
		for (auto j = i->begin(); j != i->end(); ++j) {j->clear();}
	}
}


// ped_manager_t
/*static*/ float ped_manager_t::get_ped_radius() {return 0.05*city_params.road_width;} // or should this be relative to player/camera radius?

void ped_manager_t::expand_cube_for_ped(cube_t &cube) const {
	float const radius(get_ped_radius());
	cube.expand_by_xy(radius); // PED_WIDTH_SCALE*radius for models?
	cube.z2() += PED_HEIGHT_SCALE*radius;
}

void ped_manager_t::init(unsigned num_city, unsigned num_building) {
	if (num_city == 0 && num_building == 0) return;
	timer_t timer("Gen Peds");
	peds.reserve(num_city);
	float const radius(get_ped_radius()); // base radius

	// place pedestrians in cities
	for (unsigned n = 0; n < num_city; ++n) {
		pedestrian_t ped(radius); // start with a constant radius
		assign_ped_model(ped);

		if (gen_ped_pos(ped)) {
			if (city_params.ped_speed > 0.0) {
				ped.speed = city_params.ped_speed*rgen.rand_uniform(0.5, 1.0);
				ped.vel   = rgen.signed_rand_vector_spherical_xy().get_norm()*ped.speed;
			}
			ped.ssn = (unsigned short)peds.size(); // assign init peds index so that all are unique; won't change if peds are reordered
			peds.push_back(ped);
		}
	} // for n

	// place people in buildings
	vect_building_place_t locs;
	place_building_people(locs, radius, city_params.ped_speed, num_building);
	peds_b.reserve(locs.size());
	bool const enable_bp_ai(enable_building_people_ai());

	for (auto i = locs.begin(); i != locs.end(); ++i) {
		pedestrian_t ped(radius);
		assign_ped_model(ped);
		float const angle(rgen.rand_uniform(0.0, TWO_PI));
		ped.pos   = i->p + vector3d(0.0, 0.0, ped.radius);
		ped.dir   = vector3d(cosf(angle), sinf(angle), 0.0);
		ped.speed = (enable_bp_ai ? city_params.ped_speed*rgen.rand_uniform(0.5, 0.75) : 0.0f); // small range, slower than outdoor city pedestrians
		ped.ssn   = (unsigned short)(peds.size() + peds_b.size()); // may wrap
		ped.dest_bldg   = i->bix; // store building index in dest_bldg field
		ped.in_building = 1;
		peds_b.push_back(ped);
	} // for i
	cout << "City Pedestrians: " << peds.size() << ", Building Residents: " << peds_b.size() << endl; // testing
	sort_by_city_and_plot();
}

void ped_manager_t::assign_ped_model(pedestrian_t &ped) { // Note: non-const, modifies rgen
	unsigned const num_models(ped_model_loader.num_models());
	if (num_models == 0) {ped.model_id = 0; return;} // will be unused
	ped.model_id = rgen.rand()%num_models;
	ped.radius  *= ped_model_loader.get_model(ped.model_id).scale;
	assert(ped.radius > 0.0); // no zero/negative model scales
}

struct ped_by_plot {
	bool operator()(pedestrian_t const &a, pedestrian_t const &b) const {return (a.plot < b.plot);}
};

void ped_manager_t::sort_by_city_and_plot() {
	//timer_t timer("Ped Sort"); // 0.12ms
	if (peds.empty()) return;
	bool const first_sort(by_city.empty()); // since peds can't yet move between cities, we only need to sorty by city the first time

	if (first_sort) { // construct by_city
		sort(peds.begin(), peds.end());
		unsigned const max_city(peds.back().city), max_plot(peds.back().plot);
		by_city.resize(max_city + 2); // one per city + terminator
		need_to_sort_city.resize(max_city+1, 0);

		for (unsigned city = 0, pix = 0; city <= max_city; ++city) {
			while (pix < peds.size() && peds[pix].city == city) {++pix;}
			unsigned const cur_plot((pix < peds.size()) ? peds[pix].plot : max_plot+1);
			by_city[city+1].assign(pix, cur_plot); // next city begins here
		}
	}
	else { // sort by plot within each city
		for (unsigned city = 0; city+1 < by_city.size(); ++city) {
			if (!need_to_sort_city[city]) continue;
			need_to_sort_city[city] = 0;
			sort((peds.begin() + by_plot[by_city[city].plot_ix]), (peds.begin() + by_plot[by_city[city+1].plot_ix]), ped_by_plot());
		}
	}
	// construct by_plot
	unsigned const max_plot(peds.back().plot);
	by_plot.resize((max_plot + 2), 0); // one per by_plot + terminator

	for (unsigned plot = 0, pix = 0; plot <= max_plot; ++plot) {
		while (pix < peds.size() && peds[pix].plot == plot) {++pix;}
		by_plot[plot+1] = pix; // next plot begins here
	}
	need_to_sort_peds = 0; // peds are now sorted
}

bool ped_manager_t::proc_sphere_coll(point &pos, float radius, vector3d *cnorm) const { // Note: no p_last; for potential use with ped/ped collisions
	float const rsum(get_ped_radius() + radius);

	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		cube_t const city_bcube(get_expanded_city_bcube_for_peds(city));
		if (pos.z > city_bcube.z2() + rsum) continue; // above the peds
		if (!sphere_cube_intersect_xy(pos, radius, city_bcube)) continue;

		for (unsigned plot = by_city[city].plot_ix; plot < by_city[city+1].plot_ix; ++plot) {
			cube_t const plot_bcube(get_expanded_city_plot_bcube_for_peds(city, plot));
			if (!sphere_cube_intersect_xy(pos, radius, plot_bcube)) continue;
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				if (!dist_less_than(pos, peds[i].pos, rsum)) continue;
				if (cnorm) {*cnorm = (pos - peds[i].pos).get_norm();}
				return 1; // return on first coll
			}
		} // for plot
	} // for city
	return 0;
}

bool ped_manager_t::line_intersect_peds(point const &p1, point const &p2, float &t) const {
	bool ret(0);

	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		if (!get_expanded_city_bcube_for_peds(city).line_intersects(p1, p2)) continue;

		for (unsigned plot = by_city[city].plot_ix; plot < by_city[city+1].plot_ix; ++plot) {
			if (!get_expanded_city_plot_bcube_for_peds(city, plot).line_intersects(p1, p2)) continue;
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				float tmin(0.0);
				if (line_sphere_int_closest_pt_t(p1, p2, peds[i].pos, peds[i].radius, tmin) && tmin < t) {t = tmin; ret = 1;}
			}
		} // for plot
	} // for city
	return ret;
}

void ped_manager_t::destroy_peds_in_radius(point const &pos_in, float radius) {
	point const pos(pos_in - get_camera_coord_space_xlate());
	bool const is_pt(radius == 0.0);
	float const rsum(get_ped_radius() + radius);

	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		cube_t const city_bcube(get_expanded_city_bcube_for_peds(city));
		if (pos.z > city_bcube.z2() + rsum) continue; // above the peds
		if (is_pt ? !city_bcube.contains_pt_xy(pos) : !sphere_cube_intersect_xy(pos, radius, city_bcube)) continue;

		for (unsigned plot = by_city[city].plot_ix; plot < by_city[city+1].plot_ix; ++plot) {
			cube_t const plot_bcube(get_expanded_city_plot_bcube_for_peds(city, plot));
			if (is_pt ? !plot_bcube.contains_pt_xy(pos) : !sphere_cube_intersect_xy(pos, radius, plot_bcube)) continue;
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				if (!dist_less_than(pos, peds[i].pos, rsum)) continue;
				peds[i].destroy();
				ped_destroyed = 1;
			}
		} // for plot
	} // for city
}

void ped_manager_t::remove_destroyed_peds() {
	//remove_destroyed(peds); // invalidates indexing, can't do this yet
	ped_destroyed = 0;
}

void ped_manager_t::register_ped_new_plot(pedestrian_t const &ped) {
	if (!need_to_sort_city.empty()) {need_to_sort_city[ped.city] = 1;}
	need_to_sort_peds = 1;
}
void ped_manager_t::move_ped_to_next_plot(pedestrian_t &ped) {
	if (ped.next_plot == ped.plot) return; // already there (error?)
	ped.plot = ped.next_plot; // assumes plot is adjacent; doesn't actually do any moving, only registers the move
	register_ped_new_plot(ped);
}

void ped_manager_t::next_frame() {
	if (!animate2) return; // nothing to do (only applies to moving peds)
	float const delta_dir(1.2*(1.0 - pow(0.7f, fticks))); // controls pedestrian turning rate

	// Note: peds and peds_b can be processed in parallel, but that doesn't seem to make a significant difference in framerate
	if (!peds.empty()) {
		//timer_t timer("Ped Update"); // ~4.2ms for 10K peds
		// Note: should make sure this is after sorting cars, so that road_ix values are actually in order; however, that makes things slower, and is unlikely to make a difference
#pragma omp critical(modify_car_data)
		car_manager.extract_car_data(cars_by_city);

		if (ped_destroyed) {remove_destroyed_peds();} // at least one ped was destroyed in the previous frame - remove it/them
		static bool first_frame(1);

		if (first_frame) { // choose initial ped destinations (must be after building setup, etc.)
			for (auto i = peds.begin(); i != peds.end(); ++i) {choose_dest_building_or_parked_car(*i);}
		}
		for (auto i = peds.begin(); i != peds.end(); ++i) {i->next_frame(*this, peds, (i - peds.begin()), rgen, delta_dir);}
		if (need_to_sort_peds) {sort_by_city_and_plot();}
		first_frame = 0;
	}
	if (!peds_b.empty()) {update_building_ai_state(peds_b, delta_dir);} // update people in buildings
}

pedestrian_t const *ped_manager_t::get_ped_at(point const &p1, point const &p2) const { // Note: p1/p2 in local TT space
	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		if (!get_expanded_city_bcube_for_peds(city).line_intersects(p1, p2)) continue; // skip

		for (unsigned plot = by_city[city].plot_ix; plot < by_city[city+1].plot_ix; ++plot) {
			if (!get_expanded_city_plot_bcube_for_peds(city, plot).line_intersects(p1, p2)) continue; // skip
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				if (line_sphere_intersect(p1, p2, peds[i].pos, peds[i].radius)) {return &peds[i];}
			}
		} // for plot
	} // for city
	return nullptr; // no ped found
}

void ped_manager_t::get_peds_crossing_roads(ped_city_vect_t &pcv) const {
	//timer_t timer("Get Peds Corssing Roads");
	pcv.clear();

	for (auto i = peds.begin(); i != peds.end(); ++i) {
		if (!i->in_the_road || i->is_stopped || i->destroyed) continue; // not actively crossing the road
		bool const road_dim(fabs(i->vel.y) < fabs(i->vel.x)); // ped should be moving across the road, so velocity should give us the road dim (opposite of velocity dim)
		int const road_ix(get_road_ix_for_ped_crossing(*i, road_dim));
		if (road_ix >= 0) {pcv.add_ped(*i, road_ix);}
	} // for i
}

bool ped_manager_t::has_nearby_car(pedestrian_t const &ped, bool road_dim, float delta_time, vect_cube_t *dbg_cubes) const {
	int const road_ix(get_road_ix_for_ped_crossing(ped, road_dim));
	if (road_ix < 0) return 0; // failed for some reason, assume the answer is no
	// Note: we only use road_ix, not seg_ix, because we need to find cars that are in adjacent segments to the ped (and it's difficult to get seg_ix)
	return has_nearby_car_on_road(ped, road_dim, (unsigned)road_ix, delta_time, dbg_cubes);
}

bool ped_manager_t::has_nearby_car_on_road(pedestrian_t const &ped, bool dim, unsigned road_ix, float delta_time, vect_cube_t *dbg_cubes) const {
	if (ped.city >= cars_by_city.size()) return 0; // no cars in this city? should be rare, unless cars aren't enabled
	car_city_vect_t const &cv(cars_by_city[ped.city]);
	point const &pos(ped.pos);

	for (unsigned dir = 0; dir < 2; ++dir) { // look both ways before crossing
		auto const &cars(cv.cars[dim][dir]);
		car_base_t ref_car; ref_car.cur_city = ped.city; ref_car.cur_road = road_ix;
		auto range_start(std::lower_bound(cars.begin(), cars.end(), ref_car, comp_car_road())); // binary search acceleration
		float const speed_mult(CAR_SPEED_SCALE*city_params.car_speed), pos_min(pos[dim] - ped.radius), pos_max(pos[dim] + ped.radius);
		auto closest_car(cars.end());

		for (auto it = range_start; it != cars.end(); ++it) {
			car_base_t const &c(*it);
			assert(c.cur_city == ped.city && c.dim == dim && c.dir == (dir != 0));
			if (c.cur_road != road_ix) break; // different road, done
			float const val(c.bcube.d[dim][!dir]); // back end of the car
			if (dir) {if (val > pos_max) break;   } // already passed the ped, not a threat - done (cars are sorted in this dim)
			else     {if (val < pos_min) continue;} // already passed the ped, not a threat - skip to next car
			if (closest_car == cars.end()) {closest_car = it;} // first threatening car
			else {
				float const val2(closest_car->bcube.d[dim][!dir]);
				if (dir ? (val > val2) : (val < val2)) {closest_car = it;} // this car is closer
			}
			if (!dir) break; // no cars can be closer than this (cars are sorted in this dim)
		} // for it
		if (closest_car == cars.end()) continue; // no car found
		car_base_t const &c(*closest_car);
		float lo(c.bcube.d[dim][0]), hi(c.bcube.d[dim][1]), travel_dist(0.0);

		if (lo >= pos_max || hi <= pos_min) { // current car doesn't already overlap, do more work to determine if it will overlap pos sometime in the near future
			if (c.turn_dir != TURN_NONE) {} // car is turning; since it's already on this road, it must be turning off this road; don't update its future position
			else if (c.stopped_at_light) { // check if the light will change in time for this car to reach pos; would it be easier to check the crosswalk signal?
				float const dist_gap(dir ? (pos_min - hi) : (lo - pos_max)), time_to_close(dist_gap/(speed_mult*c.max_speed));
				if (dist_gap > 0.0 && get_car_isec(c).will_be_green_light_in(c, time_to_close/TICKS_PER_SECOND)) {travel_dist = 1.01*dist_gap;} // move it just enough to cover the gap
			}
			else if (!c.is_stopped()) { // moving and not turning; assume it may be accelerating, and could reach max_speed by the time it passes pos
				//travel_dist = delta_time*speed_mult*c.cur_speed; // Note: inaccurate if car is accelerating
				travel_dist = delta_time*speed_mult*c.max_speed; // conservative - safer
			}
			//max_eq(travel_dist, 0.5f*c.get_length()); // extend by half a car length to avoid letting pedestrians cross in between cars stopped at a light (no longer needed?)
			if (dir) {hi += travel_dist;} else {lo -= travel_dist;}
			assert(lo < hi);
		}
		if (dbg_cubes) {
			cube_t cube(c.bcube);
			cube.d[dim][0] = lo; cube.d[dim][1] = hi;
			dbg_cubes->push_back(cube);
		}
		if (lo < pos_max && hi > pos_min) return 1; // overlaps current or future car in dim
	} // for dir
	return 0;
}

bool ped_manager_t::has_car_at_pt(point const &pos, unsigned city, bool is_parked) const {
	
	assert(city < cars_by_city.size());
	car_city_vect_t const &cv(cars_by_city[city]);

	if (is_parked) { // handle parked cars case
		for (auto c = cv.parked_car_bcubes.begin(); c != cv.parked_car_bcubes.end(); ++c) {
			if (c->contains_pt_xy(pos)) return 1;
		}
	}
	else { // check all dims/dirs of non-parked cars
		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				auto const &cars(cv.cars[dim][dir]);

				for (auto c = cars.begin(); c != cars.end(); ++c) {
					if (c->bcube.contains_pt_xy(pos)) return 1;
				}
			}
		}
	}
	return 0;
}

bool ped_manager_t::choose_dest_parked_car(unsigned city_id, unsigned &plot_id, unsigned &car_ix, point &car_center) {
	assert(city_id < cars_by_city.size());
	car_city_vect_t const &cv(cars_by_city[city_id]);
	if (cv.parked_car_bcubes.empty()) return 0;
	car_ix     = rgen.rand() % cv.parked_car_bcubes.size(); // Note: car_ix is stored in ped dest_bldg and doesn't get used after that
	plot_id    = cv.parked_car_bcubes[car_ix].ix;
	car_center = cv.parked_car_bcubes[car_ix].get_cube_center();
	return 1;
}

// drawing
void begin_ped_sphere_draw(shader_t &s, colorRGBA const &color, bool &in_sphere_draw, bool textured) {
	if (in_sphere_draw) return;
	if (!textured) {select_texture(WHITE_TEX);} // currently not textured
	s.set_cur_color(color);
	begin_sphere_draw(textured);
	in_sphere_draw = 1;
}
void end_sphere_draw(bool &in_sphere_draw) {
	if (!in_sphere_draw) return;
	end_sphere_draw();
	in_sphere_draw = 0;
}

void draw_colored_cube(cube_t const &c, colorRGBA const &color, shader_t &s) {
	s.set_cur_color(PURPLE);
	draw_simple_cube(c);
}

void pedestrian_t::debug_draw(ped_manager_t &ped_mgr) const {
	cube_t plot_bcube, next_plot_bcube;
	get_plot_bcubes_inc_sidewalks(ped_mgr, plot_bcube, next_plot_bcube);
	point const orig_dest_pos(get_dest_pos(plot_bcube, next_plot_bcube, ped_mgr));
	point dest_pos(orig_dest_pos);
	if (dest_pos == pos) return; // no path, nothing to draw
	vect_cube_t dbg_cubes;
	bool const safe_to_cross(check_for_safe_road_crossing(ped_mgr, plot_bcube, next_plot_bcube, &dbg_cubes));
	if (!safe_to_cross) {assert(!dbg_cubes.empty());} // must find a blocking car
	path_finder_t path_finder(1); // debug=1
	vect_cube_t const &colliders(ped_mgr.get_colliders_for_plot(city, plot));
	vect_cube_t &avoid(path_finder.get_avoid_vector());
	get_avoid_cubes(ped_mgr, colliders, plot_bcube, next_plot_bcube, dest_pos, avoid);
	cube_t union_plot_bcube(plot_bcube);
	union_plot_bcube.union_with_cube(next_plot_bcube);
	vector<point> path;
	unsigned const ret(path_finder.run(pos, dest_pos, union_plot_bcube, 0.05*radius, dest_pos)); // 0=failed, 1=valid path, 2=init contained, 3=straight path (no collisions)
	if (ret == 0) return; // no path found
	bool const at_dest_plot(plot == dest_plot), complete(path_finder.found_complete_path());
	colorRGBA line_color(at_dest_plot ? RED : YELLOW); // paths
	colorRGBA node_color(complete ? YELLOW : ORANGE);

	if (ret == 3) { // straight line
		path.push_back(pos);
		path.push_back(dest_pos);
		if (!at_dest_plot) {line_color = ORANGE; node_color = (safe_to_cross ? GREEN : RED);} // straight line
	}
	else {path = path_finder.get_best_path();} // found a path
	vector<vert_color> line_pts;
	shader_t s;
	s.begin_color_only_shader();
	bool in_sphere_draw(0);
	begin_ped_sphere_draw(s, node_color, in_sphere_draw, 0);

	if (ret == 2) { // show segment from current pos to edge of building/car
		assert(!path.empty());
		draw_sphere_vbo(path[0], radius, 16, 0);
		line_pts.emplace_back(pos, BLUE);
		line_pts.emplace_back(path[0], BLUE);
	}
	for (auto p = path.begin(); p+1 != path.end(); ++p) { // iterate over line segments, skip last point
		point const &n(*(p+1));
		draw_sphere_vbo(n, radius, 16, 0);
		line_pts.emplace_back(*p, line_color);
		line_pts.emplace_back(n,  line_color);
	}
	if (at_dest_plot && (has_dest_car || !complete)) { // show destination when in dest plot with incomplete path or car
		s.set_cur_color(PURPLE);
		draw_sphere_vbo(orig_dest_pos, 1.5*radius, 16, 0);
	}
	end_sphere_draw(in_sphere_draw);

	for (auto i = dbg_cubes.begin(); i != dbg_cubes.end(); ++i) {
		s.set_cur_color((!safe_to_cross && i+1 == dbg_cubes.end()) ? RED : GREEN);
		draw_simple_cube(*i);
	}
	ensure_outlined_polygons();
	s.set_cur_color(CYAN);

	for (auto i = avoid.begin(); i != avoid.end(); ++i) { // draw avoid cubes
		cube_t c(*i);
		max_eq(c.z2(), (c.z1() + radius)); // make sure it's nonzero area
		draw_simple_cube(c);
	}
	if (has_dest_bldg   ) {draw_colored_cube(get_building_bcube(dest_bldg), PURPLE, s);} // draw dest building bcube
	if (collided        ) {draw_colored_cube(get_bcube(), RED, s);} // show marker if collided this frame
	else if (in_the_road) {draw_colored_cube(get_bcube(), GREEN, s);}
	set_fill_mode(); // reset
	draw_verts(line_pts, GL_LINES);
	s.end_shader();
}

void ped_manager_t::next_animation() {
	unsigned const NUM_ANIMATIONS = 7; // including null animation
	string const animation_names[NUM_ANIMATIONS] = {"The Slide", "Walking", "The Bunny Hop", "The Flip", "The Twirl", "Marching", "Walk Like an Alien"};
	animation_id = (animation_id + 1) % NUM_ANIMATIONS;
	print_text_onscreen(animation_names[animation_id], WHITE, 1.5, 2*TICKS_PER_SECOND, 1);
}

void ped_manager_t::draw(vector3d const &xlate, bool use_dlights, bool shadow_only, bool is_dlight_shadows) {
	if (empty()) return;
	if (is_dlight_shadows && !city_params.car_shadows) return; // use car_shadows as ped_shadows
	if (shadow_only && !is_dlight_shadows) return; // don't add to precomputed shadow map
	//timer_t timer("Ped Draw"); // ~1ms
	bool const use_models(ped_model_loader.num_models() > 0), enable_animations(use_models);
	float const def_draw_dist((use_models ? 500.0 : 2000.0)*get_ped_radius());
	float const draw_dist(is_dlight_shadows ? 0.8*camera_pdu.far_ : def_draw_dist), draw_dist_sq(draw_dist*draw_dist); // smaller view dist for models
	pos_dir_up pdu(camera_pdu); // decrease the far clipping plane for pedestrians
	pdu.far_ = draw_dist;
	pdu.pos -= xlate; // adjust for local translate
	dstate.xlate = xlate;
	dstate.set_enable_normal_map(use_models && use_model3d_bump_maps());
	fgPushMatrix();
	translate_to(xlate);
	if (enable_animations) {enable_animations_for_shader(dstate.s);}
	dstate.pre_draw(xlate, use_dlights, shadow_only);
	if (enable_animations) {dstate.s.add_uniform_int("animation_id", animation_id);}
	if (!shadow_only) {dstate.s.add_uniform_float("hemi_lighting_normal_scale", 0.0);} // disable hemispherical lighting normal because the transforms make it incorrect
	bool in_sphere_draw(0);

	for (unsigned city = 0; city+1 < by_city.size(); ++city) {
		if (!pdu.cube_visible(get_expanded_city_bcube_for_peds(city))) continue; // city not visible - skip
		unsigned const plot_start(by_city[city].plot_ix), plot_end(by_city[city+1].plot_ix);
		assert(plot_start <= plot_end);

		for (unsigned plot = plot_start; plot < plot_end; ++plot) {
			assert(plot < by_plot.size());
			cube_t const plot_bcube(get_expanded_city_plot_bcube_for_peds(city, plot));
			if (is_dlight_shadows && !plot_bcube.closest_dist_less_than(pdu.pos, draw_dist)) continue; // plot is too far away
			if (!pdu.cube_visible(plot_bcube)) continue; // plot not visible - skip
			unsigned const ped_start(by_plot[plot]), ped_end(by_plot[plot+1]);
			assert(ped_start <= ped_end);
			if (ped_start == ped_end) continue; // no peds on this plot
			dstate.ensure_shader_active(); // needed for use_smap=0 case
			if (!shadow_only) {dstate.begin_tile(plot_bcube.get_cube_center(), 1);} // use the plot's tile's shadow map

			for (unsigned i = ped_start; i < ped_end; ++i) { // peds iteration
				assert(i < peds.size());
				pedestrian_t const &ped(peds[i]);
				assert(ped.city == city && ped.plot == plot);
				if (!draw_ped(ped, dstate.s, pdu, xlate, def_draw_dist, draw_dist_sq, in_sphere_draw, shadow_only, is_dlight_shadows, enable_animations)) continue;

				if (dist_less_than(pdu.pos, ped.pos, 0.5*draw_dist)) { // fake AO shadow at below half draw distance
					float const ao_radius(0.6*ped.radius);
					float const zval(get_city_plot_for_peds(ped.city, ped.plot).z2() + 0.02*ped.radius); // at the feet
					point pao[4];
					
					for (unsigned n = 0; n < 4; ++n) {
						point &v(pao[n]);
						v.x = ped.pos.x + (((n&1)^(n>>1)) ? -ao_radius : ao_radius);
						v.y = ped.pos.y + ((n>>1)         ? -ao_radius : ao_radius);
						v.z = zval;
					}
					dstate.ao_qbd.add_quad_pts(pao, colorRGBA(0, 0, 0, 0.4), plus_z);
				}
			} // for i
		} // for plot
	} // for city
	end_sphere_draw(in_sphere_draw);
	if (!shadow_only) {dstate.s.add_uniform_float("hemi_lighting_normal_scale", 1.0);} // restore
	pedestrian_t const *selected_ped(nullptr);

	if (tt_fire_button_down && game_mode != 1) {
		point const p1(get_camera_pos() - xlate), p2(p1 + cview_dir*FAR_CLIP);
		pedestrian_t const *ped(get_ped_at(p1, p2));
		selected_ped_ssn = (ped ? ped->ssn : -1); // update and cache for use in later frames
		if (ped) {dstate.set_label_text(ped->str(), (ped->pos + xlate));} // found
	}
	else if (selected_ped_ssn >= 0) { // a ped was selected, iterate and find it by SSN
		for (auto i = peds.begin(); i != peds.end(); ++i) {
			if (i->ssn == selected_ped_ssn) {selected_ped = &(*i); break;}
		}
		assert(selected_ped); // must be found
	}
	dstate.post_draw();
	if (selected_ped) {selected_ped->debug_draw(*this);}
	fgPopMatrix();
	dstate.show_label_text();
}

void ped_manager_t::draw_peds_in_building(int first_ped_ix, ped_draw_vars_t const &pdv) {
	if (first_ped_ix < 0) return; // no peds - error?
	assert((unsigned)first_ped_ix < peds_b.size());
	assert(peds_b[first_ped_ix].dest_bldg == pdv.bix); // consistency check
	float const def_draw_dist(120.0*get_ped_radius()); // smaller than city peds
	float const draw_dist(pdv.shadow_only ? camera_pdu.far_ : def_draw_dist), draw_dist_sq(draw_dist*draw_dist);
	pos_dir_up pdu(camera_pdu); // decrease the far clipping plane for pedestrians
	pdu.pos -= pdv.xlate; // adjust for local translate
	bool const enable_animations(enable_building_people_ai());
	bool in_sphere_draw(0);
	if (enable_animations) {pdv.s.add_uniform_int("animation_id", animation_id);}

	// Note: no far clip adjustment or draw dist scale
	for (auto p = peds_b.begin()+first_ped_ix; p != peds_b.end(); ++p) {
		if (p->dest_bldg != pdv.bix) break; // done with this building
		
		if ((display_mode & 0x08) && !city_params.ped_model_files.empty()) { // occlusion culling, if using models
			if (pdv.building.check_obj_occluded(p->get_bcube(), pdu.pos, pdv.oc, pdv.reflection_pass)) continue;
		}
		draw_ped(*p, pdv.s, pdu, pdv.xlate, def_draw_dist, draw_dist_sq, in_sphere_draw, pdv.shadow_only, pdv.shadow_only, enable_animations);
	} // for p
	end_sphere_draw(in_sphere_draw);
	pdv.s.upload_mvm(); // seems to be needed after applying model transforms, not sure why
	if (enable_animations) {pdv.s.add_uniform_int("animation_id", 0);} // make sure to leave animations disabled so that they don't apply to buildings
}

void ped_manager_t::get_ped_bcubes_for_building(int first_ped_ix, unsigned bix, vect_cube_t &bcubes, bool moving_only) const {
	if (first_ped_ix < 0) return; // no peds
	assert((unsigned)first_ped_ix < peds_b.size());
	assert(peds_b[first_ped_ix].dest_bldg == bix); // consistency check

	for (auto p = peds_b.begin()+first_ped_ix; p != peds_b.end(); ++p) {
		if (p->dest_bldg != bix) break; // done with this building
		if (moving_only && p->is_waiting_or_stopped()) continue;
		if (!p->destroyed) {bcubes.push_back(p->get_bcube());}
	}
}

bool ped_manager_t::draw_ped(pedestrian_t const &ped, shader_t &s, pos_dir_up const &pdu, vector3d const &xlate, float def_draw_dist, float draw_dist_sq,
	bool &in_sphere_draw, bool shadow_only, bool is_dlight_shadows, bool enable_animations)
{
	if (ped.destroyed) return 0; // skip
	float const dist_sq(p2p_dist_sq(pdu.pos, ped.pos));
	if (dist_sq > draw_dist_sq) return 0; // too far - skip
	if (is_dlight_shadows && !dist_less_than(pre_smap_player_pos, ped.pos, 0.4*def_draw_dist)) return 0; // too far from the player
	if (is_dlight_shadows && !sphere_in_light_cone_approx(pdu, ped.pos, 0.5*ped.get_height())) return 0;

	if (ped_model_loader.num_models() == 0 || !ped_model_loader.is_model_valid(ped.model_id)) {
		if (!pdu.sphere_visible_test(ped.pos, ped.radius)) return 0; // not visible - skip
		if (enable_animations) {s.add_uniform_float("animation_time", 0.0);}
		begin_ped_sphere_draw(s, YELLOW, in_sphere_draw, 0);
		int const ndiv = 16; // currently hard-coded
		draw_sphere_vbo(ped.pos, ped.radius, ndiv, 0);
	}
	else {
		cube_t const bcube(ped.get_bcube());
		if (!pdu.sphere_visible_test(bcube.get_cube_center(), 0.5*ped.get_height())) return 0; // not visible - skip
		if (!ped.in_building && dstate.is_occluded(bcube)) return 0; // only check occlusion for expensive ped models, and for peds outside buildings
		end_sphere_draw(in_sphere_draw);
		bool const low_detail(!shadow_only && dist_sq > 0.25*draw_dist_sq); // low detail for non-shadow pass at half draw dist
		if (enable_animations) {s.add_uniform_float("animation_time", ped.anim_time);}
		vector3d dir_horiz(ped.dir);
		dir_horiz.z = 0.0; // always face a horizontal direction, even if walking on a slope
		dir_horiz.normalize();
		//colorRGBA const &color(ped.following_player ? RED : WHITE); // force red when following player, for debugging purposes
		//colorRGBA const &color(ped.on_stairs() ? RED : ALPHA0);
		//colorRGBA const &color((ped.retreat_time > 0.0) ? RED : ALPHA0);
		colorRGBA const &color(ALPHA0); // A=0.0, leave unchanged
		ped_model_loader.draw_model(s, ped.pos, bcube, dir_horiz, color, xlate, ped.model_id, shadow_only, low_detail, enable_animations);
	}
	return 1;
}

void ped_manager_t::draw_player_model(shader_t &s, vector3d const &xlate, bool shadow_only) {
	if (ped_model_loader.num_models() == 0) { // no model - draw as sphere
		if (!shadow_only) return; // sphere is only used for shadows
		point const player_pos(pre_smap_player_pos - vector3d(0.0, 0.0, 0.5f*camera_zh)); // shift to center of player height
		draw_sphere_vbo(player_pos, 0.5f*CAMERA_RADIUS, N_SPHERE_DIV, 0); // use a smaller radius
		return;
	}
	bool const enable_animations(enable_building_people_ai());
	static float player_anim_time(0.0);
	static point prev_player_pos;
	
	if (enable_animations && p2p_dist_xy(pre_smap_player_pos, prev_player_pos) > 0.01*CAMERA_RADIUS) { // don't include minor differences related to turning in place
		prev_player_pos   = pre_smap_player_pos;
		player_anim_time += fticks*city_params.ped_speed;
	}
	pos_dir_up pdu(camera_pdu);
	pdu.pos -= xlate; // adjust for local translate
	unsigned const model_id = 0; // player is always the first model specified/loaded
	if (city_params.num_peds == 0 && city_params.num_building_peds == 0) {ped_model_loader.load_model_id(model_id);} // only need to load this particular model
	if (enable_animations) {s.add_uniform_int("animation_id", animation_id);}
	float const player_eye_height(CAMERA_RADIUS + camera_zh), player_height(1.1*player_eye_height), player_radius(player_height/PED_HEIGHT_SCALE);
	point const pos(pre_smap_player_pos + vector3d(0.0, 0.0, (player_radius - player_eye_height)));
	vector3d const dir_horiz(vector3d(cview_dir.x, cview_dir.y, 0.0).get_norm()); // always face a horizontal direction, even if walking on a slope
	cube_t bcube;
	bcube.set_from_sphere(pos, PED_WIDTH_SCALE*player_radius);
	bcube.z1() = pos.z - player_radius;
	bcube.z2() = bcube.z1() + player_height;
	if (enable_animations) {s.add_uniform_float("animation_time", player_anim_time);}
	ped_model_loader.draw_model(s, pos, bcube, dir_horiz, ALPHA0, xlate, model_id, shadow_only, 0, enable_animations);
	s.upload_mvm(); // not sure if this is needed
	if (enable_animations) {s.add_uniform_int("animation_id", 0);} // make sure to leave animations disabled so that they don't apply to buildings
}

