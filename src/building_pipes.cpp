// 3D World - Building Basement Pipes
// by Frank Gennari 10/01/2024

#include "function_registry.h"
#include "buildings.h"
#include <cfloat> // for FLT_MAX


bool line_int_cubes_exp(point const &p1, point const &p2, vect_cube_t const &cubes, vector3d const &expand);
void add_pg_obstacles(vect_room_object_t const &objs, unsigned objs_start, unsigned objs_end, vect_cube_t &walls, vect_cube_t &beams, vect_cube_t &obstacles);
void subtract_cubes_from_cube_split_in_dim(cube_t const &c, vect_cube_t const &sub, vect_cube_t &out, vect_cube_t &out2, unsigned dim);


// find the closest wall (including room wall) to this location, avoiding obstacles, and shift outward by radius; routes in X or Y only, for now
point get_closest_wall_pos(point const &pos, float radius, cube_t const &room, vect_cube_t const &walls, vect_cube_t const &obstacles, bool vertical) {
	if (!room.contains_pt_xy_exp(pos, radius)) {return pos;} // error?
	// what about checking pos intersecting walls or obstacles? is that up to the caller to handle?
	vector3d const expand(radius, radius, radius);
	point best(pos);
	float dmin(room.dx() + room.dy()); // use an initial distance larger than what we can return

	for (unsigned dim = 0; dim < 2; ++dim) { // check room/part exterior walls first
		for (unsigned dir = 0; dir < 2; ++dir) {
			float const val(room.d[dim][dir] + (dir ? -radius : radius)), dist(fabs(val - pos[dim])); // shift val inward
			if (dist >= dmin) continue;
			point cand(pos);
			cand[dim] = val;
			// check walls as well, even though any wall hit should be replaced with a closer point below
			if (line_int_cubes_exp(pos, cand, obstacles, expand) || line_int_cubes_exp(pos, cand, walls, expand)) continue;
			if (vertical && line_int_cubes_exp(cand, point(cand.x, cand.y, room.z1()), obstacles, expand)) continue; // check for vertical obstacles
			best = cand; dmin = dist; // success
		} // for dir
	} // for dim
	for (cube_t const &wall : walls) { // check all interior walls
		for (unsigned dim = 0; dim < 2; ++dim) {
			if (pos[!dim] < wall.d[!dim][0]+radius || pos[!dim] > wall.d[!dim][1]-radius) continue; // doesn't project in this dim
			bool const dir(wall.get_center_dim(dim) < pos[dim]);
			float const val(wall.d[dim][dir] - (dir ? -radius : radius)), dist(fabs(val - pos[dim])); // shift val outward
			if (dist >= dmin) continue;
			point cand(pos);
			cand[dim] = val;
			if (line_int_cubes_exp(pos, cand, obstacles, expand)) continue; // check obstacles only
			if (vertical && line_int_cubes_exp(cand, point(cand.x, cand.y, room.z1()), obstacles, expand)) continue; // check for vertical obstacles
			best = cand; dmin = dist; // success
		} // for dim
	}
	return best;
}

float get_merged_pipe_radius(float r1, float r2, float exponent) {return pow((pow(r1, exponent) + pow(r2, exponent)), 1/exponent);}

float get_merged_risers_radius(vect_riser_pos_t const &risers, int exclude_flow_dir=2) {
	float radius(0.0);
	
	for (riser_pos_t const &p : risers) {
		if ((int)p.flow_dir == exclude_flow_dir) continue;
		radius = get_merged_pipe_radius(radius, + p.radius, 4.0); // higher exponent to avoid pipes that are too large
	}
	return radius;
}

enum {PIPE_RISER=0, PIPE_CONN, PIPE_MAIN, PIPE_MEC, PIPE_EXIT, PIPE_EXTB, PIPE_FITTING};

void expand_cube_except_in_dim(cube_t &c, float expand, unsigned not_dim) {
	c.expand_by(expand);
	c.expand_in_dim(not_dim, -expand); // oops, we shouldn't have expanded in this dim
}

struct pipe_t {
	point p1, p2;
	float radius;
	unsigned dim, type, end_flags; // end_flags: 1 bit is low end, 2 bit is high end; 4 bit for round lo end, 8 bit for round hi end
	bool connected=0, conn_dir=0, outside_pg=0, for_extb=0;

	pipe_t(point const &p1_, point const &p2_, float radius_, unsigned dim_, unsigned type_, unsigned end_flags_, bool for_extb_=0) :
		p1(p1_), p2(p2_), radius(radius_), dim(dim_), type(type_), end_flags(end_flags_), connected(type != PIPE_RISER), for_extb(for_extb_) {}
	float get_length() const {return fabs(p2[dim] - p1[dim]);}

	cube_t get_bcube() const {
		cube_t bcube(p1, p2);
		expand_cube_except_in_dim(bcube, radius, dim);
		return bcube;
	}
};

void add_insul_exclude(pipe_t const &pipe, cube_t const &fitting, vect_cube_t &insul_exclude) {
	cube_t ie(fitting);
	ie.expand_by(2.0*pipe.radius); // expand in all dims; needed for right angle and T junctions because they're in multiple dims
	ie.expand_in_dim(pipe.dim, 1.0*pipe.radius); // expand a bit more in pipe dim
	insul_exclude.push_back(ie);
}
bool has_int_obstacle_or_parallel_wall(cube_t const &c, vect_cube_t const &obstacles, vect_cube_t const &walls) {
	if (has_bcube_int(c, obstacles)) return 1;
	bool const pipe_dim(c.dx() < c.dy());

	for (cube_t const &wall : walls) {
		bool const wall_dim(wall.dy() < wall.dx());
		if (wall_dim != pipe_dim && wall.intersects(c)) return 1; // check if wall and pipe are parallel
	}
	return 0;
}
void add_pass_through_fittings(room_object_t const &pipe, vect_cube_t const &walls, float fitting_len,
	float fitting_expand, unsigned dim, colorRGBA const &color, vect_room_object_t &objs)
{
	float const min_ext(2.0*fitting_len);

	for (cube_t const &wall : walls) {
		if (!wall.intersects(pipe)) continue;
		if (wall.d[dim][0] < pipe.d[dim][0] - min_ext || wall.d[dim][1] > pipe.d[dim][1] + min_ext) continue; // no space for fitting
		room_object_t pf(pipe);
		pf.flags = (RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI); // draw both ends as flat
		pf.color = color;
		pf.d[dim][0] = wall.d[dim][0] - fitting_len;
		pf.d[dim][1] = wall.d[dim][1] + fitting_len;
		expand_cube_except_in_dim(pf, 0.9*fitting_expand, dim); // expand slightly, less than regular fittings to prevent Z-fighting when overlapped near ends
		objs.push_back(pf);
	} // for wall
}
void add_hanging_pipe_bracket(cube_t const &pipe, float len_pos, float ceiling_zval, float radius, float length_factor, bool dim, unsigned room_id,
	float tot_light_amt, unsigned pipe_conn_start, vect_room_object_t &objs, vect_cube_t const &obstacles, vect_cube_t const &walls)
{
	unsigned const pipe_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING);
	float const radius_expand(1.12*radius - 0.5*pipe.dz()); // larger than outer radius; radius can be larger than pipe radius (pipe.dz()) for insulated pipes
	cube_t bracket(pipe);
	set_wall_width(bracket, len_pos, length_factor*radius, dim);
	bracket.expand_in_dim(!dim, radius_expand);
	bracket.expand_in_dim(2,    radius_expand);
	cube_t bc(bracket);
	bc.z2() = ceiling_zval; // extend up to ceiling
	if (has_bcube_int(bc, obstacles) || has_bcube_int(bc, walls)) return;
	if (has_bcube_int(bc, objs, pipe_conn_start)) return; // is this possible given the existing constraints on pipe placement? maybe only for vertically stacked pipes
	objs.emplace_back(bracket, TYPE_PIPE, room_id, dim, 0, (pipe_flags | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, LT_GRAY);

	if (bracket.z2() < ceiling_zval) { // add a vertical bolt into the ceiling
		cube_t bolt;
		set_cube_zvals(bolt, pipe.z2(), ceiling_zval);
		set_wall_width(bolt, bracket.get_center_dim( dim), 0.2*radius,  dim);
		set_wall_width(bolt, bracket.get_center_dim(!dim), 0.2*radius, !dim);
		objs.emplace_back(bolt, TYPE_PIPE, room_id, 0, 1, pipe_flags, tot_light_amt, SHAPE_CYLIN, GRAY); // vertical, no ends
	}
}

float const FITTING_LEN(1.25), FITTING_RADIUS(1.1); // relative to radius

float get_pipe_dist_to_wall(float radius, float trim_thickness) {return max(FITTING_RADIUS*radius, (radius + 2.0f*trim_thickness));}

bool building_t::cube_intersects_basement_or_extb_room(cube_t const &c, bool check_tunnel_pipes) const {
	return (c.intersects(get_basement()) || cube_intersects_extb_room(c, check_tunnel_pipes));
}
bool building_t::cube_intersects_extb_room(cube_t const &c, bool check_tunnel_pipes) const {
	if (!has_ext_basement() || !c.intersects(interior->basement_ext_bcube)) return 0;

	for (auto r = interior->ext_basement_rooms_start(); r != interior->rooms.end(); ++r) {
		if (c.intersects(*r)) return 1;
	}
	for (tunnel_seg_t const &t : interior->tunnels) {
		if (c.intersects(t.bcube_ext)) return 1;
		if (!check_tunnel_pipes) continue;

		// check horizontal and vertical pipe segments connected to tunnels; these only exist for buildings where room geom has been generated
		for (tunnel_conn_t const &conn : t.conns) {
			if (c.intersects(t.get_conn_bcube(conn))) return 1;
		}
	} // for r
	if (has_pool() && c.intersects(interior->pool)) return 1;
	if (has_mall() && has_bcube_int(c, interior->mall_info->ext_stairs_elevators)) return 1; // check mall back hallway stairs and elevator
	return 0;
}

void make_pipes_dirty(vect_room_object_t &objs, unsigned pipes_start) {
	assert(pipes_start <= objs.size());

	for (auto i = objs.begin()+pipes_start; i != objs.end(); ++i) {
		if (i->type == TYPE_PIPE) {i->flags |= RO_FLAG_BROKEN;}
	}
}

bool building_t::add_basement_pipes(vect_cube_t const &obstacles, vect_cube_t const &walls, vect_cube_t const &beams, vect_riser_pos_t const &risers,
	vect_cube_t &pipe_cubes, unsigned room_id, unsigned num_floors, unsigned objs_start, float ceil_zval, rand_gen_t &rgen, unsigned pipe_type, bool allow_place_fail)
{
	//highres_timer_t timer("add_basement_pipes");
	assert(pipe_type < NUM_PIPE_TYPES);
	if (risers.empty()) return 0; // can happen for hot water pipes when there are no hot water fixtures
	float const max_pipe_radius_mult[NUM_PIPE_TYPES] = {0.75, 0.25, 0.2, 0.15}; // sewer, cold water, hot water, gas; in muptiples of wall thickness
	bool const is_hot_water(pipe_type == PIPE_TYPE_HW), is_closed_loop(is_hot_water), add_insul(is_hot_water);
	bool const extb_wall_dim(interior->extb_wall_dim), extb_wall_dir(interior->extb_wall_dir);
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(objs_start <= objs.size());
	cube_t const &basement(get_basement());
	float r_main(get_merged_risers_radius(risers, (is_closed_loop ? 0 : 2))); // exclude incoming water from hot water heaters for hot water pipes
	if (r_main == 0.0) return 0; // hot water heater but no hot water pipes?
	unsigned const pipes_start(objs.size());
	float const insul_thickness(0.4), min_insum_len(4.0); // both relative to pipe radius
	float const window_vspacing(get_window_vspace()), fc_thickness(get_fc_thickness()), wall_thickness(get_wall_thickness());
	float const radius_factor(add_insul ? 1.0+insul_thickness : 1.0), max_pipe_radius(max_pipe_radius_mult[pipe_type]*wall_thickness);
	min_eq(r_main, max_pipe_radius); // limit pipe radius; even with these limits, we can still have hot water and cold water clipping through each other with enough flow
	float const pipe_zval(ceil_zval   - FITTING_RADIUS*r_main); // includes clearance for fittings vs. beams (and lights - mostly)
	float const pipe_min_z1(pipe_zval - FITTING_RADIUS*r_main);
	float const align_dist(2.0*wall_thickness); // align pipes within this range (in particular sinks and stall toilets)
	float const r_main_spacing(radius_factor*r_main); // include insulation thickness for hot water pipes
	assert(pipe_zval > bcube.z1());
	vector<pipe_t> pipes, fittings, extb_conn, vert_conn_pipes;
	cube_t pipe_end_bcube;
	unsigned num_conn_segs(0);
	// build random shifts table; make consistent per pipe to preserve X/Y alignments
	unsigned const NUM_SHIFTS = 21; // {0,0} + 20 random shifts
	vector3d rshifts[NUM_SHIFTS] = {};
	for (unsigned n = 1; n < NUM_SHIFTS; ++n) {rshifts[n][rgen.rand_bool()] = 0.25*window_vspacing*rgen.signed_rand_float();} // random shift in a random dir
	// determine if we need to add an entrance/exit pipe
	bool add_exit_pipe(1);

	if (is_closed_loop) { // hot water pipes; sewer pipes always have an exit and cold water always has an entrance
		for (riser_pos_t const &p : risers) {
			if (p.flow_dir == 0) {add_exit_pipe = 0; break;} // hot water flowing out of a water heater
		}
		// if we got here and add_exit_pipe=1, just create the exit pipe and assume the hot water heater is somewhere else
	}

	// seed the pipe graph with valid vertical segments and build a graph of X/Y values
	for (riser_pos_t const &p : risers) {
		assert(p.radius > 0.0);
		assert(p.pos.z > pipe_zval);
		point pos(p.pos);
		
		if (p.in_extb) { // extended basement connection has special handling
			cube_t valid_bounds(basement);
			valid_bounds.expand_by_xy(-FITTING_RADIUS*p.radius);
			valid_bounds.clamp_pt(pos);
			extb_conn.emplace_back(p.pos, pos, p.radius, extb_wall_dim, PIPE_RISER, 0); // store for later
		}
		else {
			bool valid(0);

			for (unsigned n = 0; n < NUM_SHIFTS; ++n) { // try zero + random shifts
				cube_t c(pos);
				c.expand_by_xy(p.radius);
				//c.z1() = bcube.z1(); // extend all the way down to the floor of the lowest basement
				c.z1() = pipe_min_z1; // extend down to the lowest pipe z1

				// can't place outside building bcube, or over stairs/elevators/ramps/pillars/walls/beams;
				// here beams are included because lights are attached to the underside of them, so avoiding beams should hopefully also avoid lights
				if (!bcube.contains_cube_xy(c) || has_bcube_int(c, obstacles) || has_bcube_int(c, walls) || has_bcube_int(c, beams)) {
					pos = p.pos + rshifts[n]; // apply shift
					continue;
				}
				valid = 1;
				break;
			} // for n
			if (!valid) continue; // no valid shift, skip this connection
		}
		pipes.emplace_back(point(pos.x, pos.y, pipe_zval), pos, p.radius, 2, PIPE_RISER, 0, p.in_extb); // neither end capped
		pipe_end_bcube.assign_or_union_with_cube(pipes.back().get_bcube());
	} // for risers
	if (pipes.empty()) return 0; // no valid pipes

	// calculate unique positions of pipes along the main pipe
	bool const dim(pipe_end_bcube.dx() < pipe_end_bcube.dy()); // main pipe dim
	map<float, vector<unsigned>> xy_map;

	for (auto p = pipes.begin(); p != pipes.end(); ++p) {
		unsigned const pipe_ix(p - pipes.begin());
		float &v(p->p1[dim]);
		if (p->for_extb) {xy_map[v].push_back(pipe_ix); continue;} // add without merging
		auto it(xy_map.find(v));
		if (it != xy_map.end()) {it->second.push_back(pipe_ix); continue;} // found
		bool found(0);

		// try to find an existing map value that's within align_dist of this value; messy and inefficient, but I'm not sure how else to do this
		for (auto &i : xy_map) {
			if (fabs(i.first - v) > align_dist) continue; // too far
			i.second.push_back(pipe_ix);
			v = p->p2[dim] = i.first;
			found = 1;
			break;
		}
		if (!found) {xy_map[v].push_back(pipe_ix);}
	} // for pipes

	// create main pipe that runs in the longer dim (based on drain pipe XY bounds)
	pipe_end_bcube.expand_in_dim(dim, r_main);
	// use the center of the pipes bcube to minimize run length, but clamp to the interior of the basement;
	// in the case where all risers are outside of the basement perimeter, this will make pipes run against the exterior basement wall
	float const pipes_bcube_center(max(basement.d[!dim][0]+r_main, min(basement.d[!dim][1]-r_main, pipe_end_bcube.get_center_dim(!dim))));
	float centerline(pipes_bcube_center);
	point mp[2]; // {lo, hi} ends
	bool exit_dir(0);
	point exit_pos;

	for (unsigned d = 0; d < 2; ++d) { // dim
		mp[d][ dim] = pipe_end_bcube.d[dim][d];
		mp[d][!dim] = centerline;
		mp[d].z     = pipe_zval;
	}
	// shift pipe until it clears all obstacles
	float step_dist(2.0*r_main);
	float const step_area(bcube.get_sz_dim(!dim)); // step by pipe radius
	unsigned const max_steps(max(1U, min(unsigned(step_area/step_dist), 100U))); // limit to 100 steps
	step_dist = step_area/max_steps;
	bool success(0);

	for (unsigned n = 0; n < max_steps; ++n) {
		cube_t const c(pipe_t(mp[0], mp[1], r_main_spacing, dim, PIPE_MAIN, 3).get_bcube());
		if (!has_int_obstacle_or_parallel_wall(c, obstacles, walls)) {
			success = 1;

			// check for overlap with beam running parallel to the main pipe, and reject it; mostly there to avoid blocking lights that may be on the beam
			for (cube_t const &beam : beams) {
				if (beam.get_sz_dim(dim) < beam.get_sz_dim(!dim)) continue; // beam not parallel to pipe, ignore
				if (c.intersects_xy(beam)) {success = 0; break;}
			}
			if (success) break; // success/done
		}
		float const xlate(((n>>1)+1)*((n&1) ? -1.0 : 1.0)*step_dist);
		UNROLL_2X(mp[i_][!dim] += xlate;); // try the new position
		
		if (!basement.contains_cube_xy(pipe_t(mp[0], mp[1], r_main_spacing, dim, PIPE_MAIN, 3).get_bcube())) { // outside the basement
			UNROLL_2X(mp[i_][!dim] -= 2.0*xlate;); // try shifting in the other direction
			if (!basement.contains_cube_xy(pipe_t(mp[0], mp[1], r_main_spacing, dim, PIPE_MAIN, 3).get_bcube())) break; // outside the basement in both dirs, fail
		}
	} // for n
	if (success) {
		centerline = mp[0][!dim]; // update centerline based on translate
	}
	else if (allow_place_fail) { // fail case
		// try stripping off some risers from the end and connecting the remaining ones
		vect_riser_pos_t risers_sub;
		float riser_min(FLT_MAX), riser_max(-FLT_MAX);

		for (riser_pos_t const &r : risers) {
			min_eq(riser_min, r.pos[dim]);
			max_eq(riser_max, r.pos[dim]);
		}
		if (risers.size() >= 10) { // optimization for many risers case to avoid long quadratic runtime: clip off an additional 5% from each end
			float const extra_clip(0.05*(riser_max - riser_min));
			riser_min += extra_clip;
			riser_max -= extra_clip;
		}
		for (unsigned dir = 0; dir < 2; ++dir) {
			risers_sub.clear();
			// try to remove risers at each end recursively until successful
			for (riser_pos_t const &r : risers) {
				if (dir ? (r.pos[dim] < riser_max) : (r.pos[dim] > riser_min)) {risers_sub.push_back(r);}
			}
			if (risers_sub.empty()) continue;
			assert(risers_sub.size() < risers.size());
			if (add_basement_pipes(obstacles, walls, beams, risers_sub, pipe_cubes, room_id, num_floors, objs_start, ceil_zval, rgen, pipe_type, 1)) return 1;
		} // for dir
		return 0;
	}
	else { // somewhat rare failure case; can lead to intersections between hot water and cold water pipes
		float cur_pos(centerline);
		UNROLL_2X(mp[i_][!dim] = centerline;); // else use the centerline, even though it's invalid
		// must still avoid stairs and elevators
		cube_t c(pipe_t(mp[0], mp[1], r_main_spacing, dim, PIPE_MAIN, 3).get_bcube());
		vect_cube_t avoid;
		vector_add_to(interior->stairwells, avoid);
		vector_add_to(interior->elevators,  avoid); // add elevators, but escalators won't be in the basement

		for (cube_t &a : avoid) {
			a.expand_by_xy(0.25*wall_thickness); // account for stairs walls
			if (!c.intersects(a)) continue;
			bool const move_dir(a.get_center_dim(!dim) < c.get_center_dim(!dim));
			float const new_pos(a.d[!dim][move_dir] + (move_dir ? 1.0 : -1.0)*r_main_spacing), move_amt(new_pos - cur_pos);
			c.translate_dim(!dim, move_amt);
			UNROLL_2X(mp[i_][!dim] += move_amt;);
			cur_pos = new_pos;
		} // for a
		centerline = mp[0][!dim]; // update centerline
	}
	mp[0][dim] = bcube.d[dim][1]; mp[1][dim] = bcube.d[dim][0]; // make dim range denormalized; will recalculate below with correct range
	bool const d(!dim);
	float const conn_pipe_merge_exp = 3.0; // cubic
	unsigned const conn_pipes_start(pipes.size());
	vector<float> conn_pipe_pos;
	float extra_extb_fitting_extend(0.0);
	bool had_extb_conn(0);

	// connect drains/feeders to main pipe in !dim
	for (auto const &v : xy_map) { // for each unique position along the main pipe
		float radius(0.0), range_min(centerline), range_max(centerline), unconn_radius(0.0); // range of connector perpendicular to main pipe
		float mp_pos(v.first);
		point const &ref_p1(pipes[v.second.front()].p1);
		unsigned num_keep(0);

		for (unsigned ix : v.second) {
			assert(ix < pipes.size());
			pipe_t &pipe(pipes[ix]);
			if (had_extb_conn && pipe.for_extb) continue; // already have an extended basement connector pipe
			float const val(pipe.p1[d]);

			if (!pipe.for_extb && fabs(val - centerline) < r_main) {pipe.p1[d] = pipe.p2[d] = centerline;} // shift to connect directly to main pipe since it's close enough
			else {
				float const range_ext(pipe.for_extb ? 0.0 : pipe.radius);
				float lo(val - range_ext), hi(val + range_ext);
				point p1(ref_p1), p2(p1);
				if      (lo < range_min) {p1[d] = lo; p2[d] = range_min;} // on the lo side
				else if (hi > range_max) {p1[d] = range_max; p2[d] = hi;} // on the hi side
				float const r_test(radius_factor*pipe.radius);
				bool skip(has_int_obstacle_or_parallel_wall(pipe_t(p1, p2, r_test, d, PIPE_CONN, 3).get_bcube(), obstacles, walls));

				if (skip) { // blocked, can't connect
					if (pipe.for_extb) { // extended basement pipe
						if (d != extb_wall_dim) { // try moving away from the wall in case there was an outlet conduit there
							float const offset_dist(0.35*wall_thickness); // should be about the diameter of a conduit
							vector3d offset;
							offset[!d] += (extb_wall_dir ? -1.0 : 1.0)*offset_dist;
							p1 += offset; p2 += offset;
							
							if (!has_int_obstacle_or_parallel_wall(pipe_t(p1, p2, r_test, d, PIPE_CONN, 3).get_bcube(), obstacles, walls)) {
								for (pipe_t &ep : extb_conn) { // success, move pipes and extb conn so that it matches
									if (pipe.p2 == ep.p2) {ep.p2 += offset;}
								}
								pipe.p1 += offset; pipe.p2 += offset;
								extra_extb_fitting_extend = offset_dist;
								mp_pos = p1[!d]; // update position
								skip   = 0;
							}
						}
					}
					else if ((p2[d] - p1[d]) > 8.0*pipe.radius) { // long segment: try extending halfway
						if      (lo < range_min) {p1[d] = lo = 0.5*(lo + range_min); pipe.p1[d] = pipe.p2[d] = lo + pipe.radius;} // on the lo side
						else if (hi > range_max) {p2[d] = hi = 0.5*(hi + range_max); pipe.p1[d] = pipe.p2[d] = hi - pipe.radius;} // on the hi side
						skip = has_int_obstacle_or_parallel_wall(pipe_t(p1, p2, r_test, d, PIPE_CONN, 3).get_bcube(), obstacles, walls);

						if (!skip) { // okay so far; now check that the new shortened pipe has a riser pos that's clear of walls and beams
							point const new_riser_pos((lo < range_min) ? p1 : p2); // the end that was clipped
							cube_t test_cube; test_cube.set_from_sphere(new_riser_pos, pipe.radius);
							skip = (has_bcube_int(test_cube, walls) || has_bcube_int(test_cube, beams));
						}
					}
				}
				if (skip) {
					unconn_radius = get_merged_pipe_radius(unconn_radius, pipe.radius, conn_pipe_merge_exp); // add this capacity for use in another riser
					continue;
				}
				min_eq(range_min, lo); // update ranges
				max_eq(range_max, hi);
			}
			pipe.connected = 1;
			had_extb_conn |= pipe.for_extb;
			if (unconn_radius > 0.0) {pipe.radius = get_merged_pipe_radius(pipe.radius, unconn_radius, conn_pipe_merge_exp); unconn_radius = 0.0;} // add extra capacity
			radius = get_merged_pipe_radius(radius, pipe.radius, conn_pipe_merge_exp);
			++num_keep;
		} // for ix
		if (num_keep == 0) continue; // no valid connections for this row
		min_eq(radius, 0.8f*max_pipe_radius); // a bit smaller than r_main

		// we can skip adding a connector if short and under the main pipe
		if (range_max - range_min > r_main) {
			point p1, p2;
			p1[!d] = p2[!d] = mp_pos;
			p1[ d] = range_min; p2[d] = range_max;
			p1.z   = p2.z = ref_p1.z;
			pipes.emplace_back(p1, p2, radius, d, PIPE_CONN, 3); // cap both ends

			for (unsigned ix : v.second) { // add fittings
				pipe_t const &pipe(pipes[ix]);
				if (!pipe.connected || pipe.for_extb) continue; // pipe was not connected, don't add the fitting
				float const val(pipe.p1[d]), fitting_len(FITTING_LEN*radius);
				p1[d] = val - fitting_len; p2[d] = val + fitting_len;
				fittings.emplace_back(p1, p2, FITTING_RADIUS*radius, d, PIPE_FITTING, 3);
			}
		} // end connector
		// add fitting to the main pipe
		point p1(mp[0]), p2(p1);
		float const fitting_len(FITTING_LEN*r_main);
		p1[!d] = mp_pos - fitting_len; p2[!d] = mp_pos + fitting_len;
		fittings.emplace_back(p1, p2, FITTING_RADIUS*r_main, !d, PIPE_FITTING, 3);
		// update main pipe endpoints to include this connector pipe range
		min_eq(mp[0][dim], mp_pos-radius);
		max_eq(mp[1][dim], mp_pos+radius);
		conn_pipe_pos.push_back(mp_pos);
		++num_conn_segs;
	} // for v
	if (mp[0][dim] >= mp[1][dim]) return 0; // no pipes connected to main? I guess there's nothing to do here
	unsigned const conn_pipe_fittings_end(fittings.size()), conn_pipes_end(pipes.size());

	// reroute extended basement pipe(s) from next to door to the end of the hallway; only one of them should actually be connected
	for (pipe_t const &ep : extb_conn) {
		for (pipe_t &p : pipes) {
			if (!p.connected) continue; // unconnected drain, skip
			if (p.type != PIPE_RISER || p.p2 != ep.p2) continue; // wrong pipe, or it was moved
			float const pipe_zmax(get_room(interior->ext_basement_hallway_room_id).z2() - fc_thickness - FITTING_RADIUS*p.radius);

			if (p.p1.z > pipe_zmax) { // pipe is too high, likely on the floor above; run a vertical segment to connect it
				// Note: unconnected pipes have a dead end
				if (pipe_type != PIPE_TYPE_SEWER) {p.connected = 0; continue;} // only allow for sewer pipes since water pipes may intersect sewer pipes when set to the same zval
				pipe_t vpipe(p);
				vpipe.p1.assign(ep.p2.x, ep.p2.y, pipe_zmax); // low end
				vpipe.p2.assign(ep.p2.x, ep.p2.y, p.p1.z   ); // high end
				vpipe.dim  = 2;
				vpipe.type = PIPE_MEC; // similar to vertical exit pipe with fittings on both ends
				vpipe.end_flags = 3; // round end at bottom; top end could be round, but it needs to be pulled back on the connector pipe
				if (has_bcube_int(vpipe.get_bcube(), obstacles)) {p.connected = 0; continue;} // unclear if this can happen, but should be a good check
				vert_conn_pipes.push_back(vpipe);
				p.p1.z = pipe_zmax;
			}
			p.p2.assign(ep.p1.x, ep.p1.y, p.p1.z);
			p.dim      = ep.dim;
			p.conn_dir = !extb_wall_dir;
			p.type     = PIPE_EXTB;
			if (bool(p.dim) == dim) {p.end_flags |= (1 << (p.conn_dir+2));} // flag connection point as round end if it's a right angle
			pipe_t const parent(p); // make a copy to avoid invalidating the p reference
			add_ext_basement_hallway_pipes_recur(interior->ext_basement_hallway_room_id, extb_wall_dim, pipe_type, radius_factor, parent, pipes, fittings, rgen);
			break;
		} // for p
	} // for ep
	vector_add_to(vert_conn_pipes, pipes);
	unsigned main_pipe_end_flags(0); // start with both ends unconnected
	bool no_main_pipe(0);

	if (add_exit_pipe) {
		bool tried_horizontal(0);

		for (unsigned attempt = 0; attempt < 2; ++attempt) { // make up to 2 attempts with a horizontal or vertical connection
			float const short_main_pipe((mp[1][dim] - mp[0][dim]) < 4.0*r_main); // likely a single connector pipe; common for gas pipes; vert exit pippe doesn't look good

			if (is_closed_loop || num_floors > 1 || attempt == 1 || num_conn_segs == 1 || short_main_pipe || rgen.rand_bool()) { // exit into the wall of the building
				if (tried_horizontal) break; // we got here the previous iteration; give up, exit can't be found
				tried_horizontal = 1;
				// Note: if roads are added for secondary buildings, we should have the exit on the side of the building closest to the road
				bool const first_dir((basement.d[dim][1] - mp[1][dim]) < (mp[0][dim] - basement.d[dim][0])); // closer basement exterior wall
				bool added_exit(0);

				for (unsigned d = 0; d < 2; ++d) { // dir
					bool const dir(bool(d) ^ first_dir);
					point ext[2]  = {mp[dir], mp[dir]};
					ext[dir][dim] = basement.d[dim][dir]; // shift this end to the basement wall
					if (has_int_obstacle_or_parallel_wall(pipe_t(ext[0], ext[1], r_main_spacing, dim, PIPE_MAIN, 0).get_bcube(), obstacles, walls)) continue;
					if (dir) {max_eq(mp[1][dim], ext[1][dim]);} else {min_eq(mp[0][dim], ext[0][dim]);} // extend main pipe to exit point on basement wall
					main_pipe_end_flags = (dir ? 2 : 1); // connect the end going to the exit
					added_exit = 1;
					break;
				} // for d
				if (added_exit) break;
				// no straight segment? how about a right angle?
				bool first_side(0);
				if (centerline == pipes_bcube_center) {first_side = rgen.rand_bool();} // centered, choose a random side
				else {first_side = ((basement.d[!dim][1] - mp[0][!dim]) < (mp[0][!dim] - basement.d[!dim][0]));} // off-center, choose closer basement exterior wall

				for (unsigned d = 0; d < 2 && !added_exit; ++d) { // dir
					for (unsigned e = 0; e < 2; ++e) { // side
						bool const dir(bool(d) ^ first_dir), side(bool(e) ^ first_side);
						point ext[2] = {mp[dir], mp[dir]};
						ext[side][!dim] = basement.d[!dim][side]; // shift this end to the basement wall
						pipe_t const exit_pipe(ext[0], ext[1], r_main, !dim, PIPE_MEC, (side ? 1 : 2)); // add a bend in the side connecting to the main pipe
						cube_t const pipe_bcube(exit_pipe.get_bcube());
						if (has_int_obstacle_or_parallel_wall(pipe_bcube, obstacles, walls)) continue; // can't extend to the ext wall in this dim
						bool bad_place(0);

						// check if the pipe is too close to an existing conn pipe; allow it to contain the other pipe in dim
						for (auto p = pipes.begin()+conn_pipes_start; p != pipes.begin()+conn_pipes_end; ++p) {
							cube_t const other_bcube(p->get_bcube());
							if (!pipe_bcube.intersects(other_bcube)) continue;
							if (pipe_bcube.d[dim][0] > other_bcube.d[dim][0] || pipe_bcube.d[dim][1] < other_bcube.d[dim][1]) {bad_place = 1; break;}
						}
						if (bad_place) continue; // seems to usually fail
						pipes.push_back(exit_pipe);
						main_pipe_end_flags = (dir ? 2 : 1); // connect the end going to the exit connector pipe
						added_exit = 1;
						break;
					} // for e
				} // for d
				if (added_exit) break;
			}
			// create exit segment and vertical pipe into the floor
			float const exit_floor_zval(basement.z1() + fc_thickness); // on the bottom level floor
			unsigned exit_pipe_end_flags(2); // bend at the top only

			if (short_main_pipe && conn_pipe_pos.size() == 1) { // align to the single connector pipe
				point exit_conn(mp[0]); // can use either mp
				exit_conn[dim] = conn_pipe_pos.front(); // connect directly to the pipe
				exit_pos = get_closest_wall_pos(exit_conn, r_main_spacing, basement, walls, obstacles, 1);
				point exit_fc_pos(exit_pos);
				
				if (exit_pos == exit_conn) { // no valid wall found
					exit_fc_pos.z = ceil_zval; // exit back into the ceiling rather than the floor to avoid a vertical pipe in the center of the room
					exit_pipe_end_flags = 1; // bend at the bottom only
				}
				else {
					if (exit_pos[!dim] != exit_conn[!dim]) { // exit point not along the main pipe; create a right angle bend
						// this doesn't look quite right when the exit point is on the same side of the main pipe as the connector pipe, but this is relatively rare
						pipes.emplace_back(exit_conn, exit_pos, r_main, !dim, PIPE_MEC, 3); // main exit connector, bends at both ends
						exit_pipe_end_flags = 0; // the above pipe will provide the bend, so it's not needed at the top of the exit pipe
					}
					exit_fc_pos.z = exit_floor_zval; // exit down into the floor
				}
				pipes.emplace_back(exit_fc_pos, exit_pos, r_main, 2, PIPE_EXIT, exit_pipe_end_flags);
				if (conn_pipe_fittings_end > 0) {fittings.erase(fittings.begin()+conn_pipe_fittings_end-1);} // remove the conn pipe fitting as it's not needed
				no_main_pipe = 1;
				break; // done
			}
			// find the end closest to the wall
			float exit_dmin(0.0);

			for (unsigned d = 0; d < 2; ++d) { // choose the exit side
				point const cand_exit_pos(get_closest_wall_pos(mp[d], r_main_spacing, basement, walls, obstacles, 1)); // vertical=1
				float dist(p2p_dist(mp[d], cand_exit_pos));
				if (dist == 0.0) {dist = FLT_MAX;} // dist will be 0 if we fail to find a wall, so don't prefer it in that case
				if (exit_dmin == 0.0 || dist < exit_dmin) {exit_pos = cand_exit_pos; exit_dir = d; exit_dmin = dist;}
			}
			point &exit_conn(mp[exit_dir]);

			if (exit_pos[!dim] == exit_conn[!dim]) { // exit point is along the main pipe
				if (exit_conn[dim] == exit_pos[dim]) { // exit is exactly at the pipe end
					// maybe no valid wall was found above, and this location is also invalid; try with a vertical connection if this was the first attempt;
					// this is the only situation where we can get to attempt=1
					if (attempt == 0 && has_bcube_int(pipe_t(point(exit_pos.x, exit_pos.y, exit_floor_zval), exit_pos, r_main, 2, PIPE_EXIT, 0).get_bcube(), obstacles)) continue;
						
					if (has_parking_garage) { // check that no parking spaces are blocked by this pipe
						auto start(objs.begin() + objs_start);

						for (auto i = start; i != objs.end(); ++i) {
							if (i->type != TYPE_PARK_SPACE || !i->contains_pt_xy(exit_pos)) continue;
							i->remove(); // remove this parking space
							// now remove any car collider and curb associated with this parking space (placed before it, optional collider then optional curb)
							auto cur(i);
							if (cur > start && (cur-1)->type == TYPE_CURB    ) {(cur-1)->remove(); --cur;}
							if (cur > start && (cur-1)->type == TYPE_COLLIDER) {(cur-1)->remove();}
							//break; // maybe can't break here because there may be a parking space to remove on another level?
						} // for i
					}
				}
				else if ((exit_conn[dim] < exit_pos[dim]) == exit_dir) { // extend main pipe to exit point
					exit_conn[dim] = exit_pos[dim];
					main_pipe_end_flags = (exit_dir ? 2 : 1); // connect the end going to the exit
				}
				else { // exit is in the middle of the pipe; add fitting to the main pipe
					point p1(exit_pos), p2(p1);
					float const fitting_len(FITTING_LEN*r_main);
					p1[dim] -= fitting_len; p2[dim] += fitting_len;
					fittings.emplace_back(p1, p2, FITTING_RADIUS*r_main, !d, PIPE_FITTING, 3);
					exit_pipe_end_flags = 0; // no bend needed
				}
			}
			else { // create a right angle bend
				// extend by 2*radius to avoid overlapping a connector pipe or it's fittings, but keep it inside the building
				if (exit_dir) {exit_conn[dim] = exit_pos[dim] = min(exit_conn[dim]+2.0f*r_main, bcube.d[dim][1]-r_main);}
				else          {exit_conn[dim] = exit_pos[dim] = max(exit_conn[dim]-2.0f*r_main, bcube.d[dim][0]+r_main);}
				pipes.emplace_back(exit_conn, exit_pos, r_main, !dim, PIPE_MEC, 3); // main exit connector, bends at both ends
				exit_pipe_end_flags = 0; // the above pipe will provide the bend, so it's not needed at the top of the exit pipe
				main_pipe_end_flags = (exit_dir ? 2 : 1); // connect the end going to the exit connector pipe
			}
			point exit_floor_pos(exit_pos);
			exit_floor_pos.z = exit_floor_zval;
			pipes.emplace_back(exit_floor_pos, exit_pos, r_main, 2, PIPE_EXIT, exit_pipe_end_flags);
			break; // success
		} // for attempt
	} // end add_exit_pipe
	// add main pipe
	if (!no_main_pipe) {
		assert(mp[0] != mp[1]);
		pipe_t main_pipe(mp[0], mp[1], r_main, dim, PIPE_MAIN, main_pipe_end_flags);
	
		if (!main_pipe.get_bcube().is_strictly_normalized()) {
			cout << "Error: Invalid main pipe: " << TXT(r_main) << TXT(mp[0].str()) << TXT(mp[1].str()) << TXT(main_pipe.get_bcube().str()) << endl;
			assert(0);
		}
		pipes.push_back(main_pipe);
	}
	// add pipe objects: sewer: dark gray pipes / gray-brown fittings; water: copper pipes / brass fittings; hot water: white insulation
	colorRGBA const pcolors[4] = {DK_GRAY, COPPER_C, COPPER_C, GRAY}; // sewer, cw, hw, gas
	colorRGBA const fcolors[4] = {colorRGBA(0.7, 0.6, 0.5, 1.0), BRASS_C, BRASS_C, DK_GRAY}; // sewer, cw, hw, gas
	colorRGBA const &pipes_color(pcolors[pipe_type]), &fittings_color(fcolors[pipe_type]);
	float const tot_light_amt = 1.0; // to offset the darkness of the basement
	vect_cube_t insul_exclude;
	reserve_extra(pipe_cubes, pipes.size());

	for (pipe_t &p : pipes) {
		if (!p.connected) continue; // unconnected drain, skip
		cube_t const pbc(p.get_bcube());
		if (p.type != PIPE_EXTB && !cube_intersects_basement_or_extb_room(pbc)) {p.outside_pg = 1; continue;} // outside the basement, don't need to draw
		pipe_cubes.push_back(pbc);
		pipe_cubes.back().expand_in_dim(p.dim, 0.5*p.radius); // add a bit of extra space for the end cap
		bool const pdim(p.dim & 1), pdir(p.dim >> 1); // encoded as: X:dim=0,dir=0 Y:dim=1,dir=0, Z:dim=x,dir=1
		unsigned flags(0);
		if (p.type != PIPE_EXIT) {flags |= RO_FLAG_NOCOLL;} // only exit pipe has collisions enabled
		if (p.type == PIPE_CONN || p.type == PIPE_MAIN || p.type == PIPE_EXTB) {flags |= RO_FLAG_HANGING;} // hanging connector/main pipe with flat ends
		room_object_t const pipe(pbc, TYPE_PIPE, room_id, pdim, pdir, flags, tot_light_amt, SHAPE_CYLIN, pipes_color);
		objs.push_back(pipe);
		if (p.type == PIPE_EXIT && p.dim == 2) {objs.back().flags |= RO_FLAG_LIT;} // vertical exit pipes are shadow casting; applies to pipe but not fittings

		// add pipe fittings around ends and joins; only fittings have flat and round ends because raw pipe ends should never be exposed;
		// note that we may not need fittings at T-junctions for hot water pipes, but then we would need to cap the ends
		if (p.type == PIPE_RISER) continue; // not for vertical drain pipes, since they're so short and mostly hidden above the connector pipes
		float const fitting_len(FITTING_LEN*r_main), fitting_expand((FITTING_RADIUS - 1.0)*p.radius); // fitting len based on main pipe radius

		for (unsigned d = 0; d < 2; ++d) {
			if ((p.type == PIPE_CONN || p.type == PIPE_MAIN) && !(p.end_flags & (1<<d))) continue; // already have fittings added from connecting pipes
			room_object_t pf(pipe);
			pf.flags |= RO_FLAG_NOCOLL;
			if (d == 1 || p.type != PIPE_MAIN) {pf.flags |= RO_FLAG_ADJ_LO;} // skip end cap for main pipe exiting through the wall
			if (d == 0 || p.type != PIPE_MAIN) {pf.flags |= RO_FLAG_ADJ_HI;} // skip end cap for main pipe exiting through the wall
			if (p.end_flags & (1<<(d+2))) {pf.flags |= (d ? RO_FLAG_ADJ_TOP : RO_FLAG_ADJ_BOT);} // handle round end flags
			pf.color  = fittings_color;
			expand_cube_except_in_dim(pf, fitting_expand, p.dim); // expand slightly
			float fitting_len_ext(fitting_len);
			if (p.type == PIPE_EXTB && p.conn_dir == bool(d)) {fitting_len_ext += p.radius + wall_thickness + extra_extb_fitting_extend;} // extb pipe passes through wall at start pt
			pf.d[p.dim][!d] = pf.d[p.dim][d] + (d ? -1.0 : 1.0)*fitting_len_ext;
			if (p.type != PIPE_EXTB && !basement.intersects_xy(pf)) continue; // skip fittings outside the basement, except for extended basement pipe
			objs.push_back(pf);
			if (add_insul) {add_insul_exclude(p, pf, insul_exclude);}

			if (p.type == PIPE_MEC || p.type == PIPE_EXIT) {
				if (p.end_flags & (1<<d)) { // connector or exit pipe with a round bend needs special handling
					objs.back().flags &= ~(d ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI); // unset end flag on the end that was cut to length, since that's not a bend
					// create a second fitting segment for the flat end; the sides will overlap with the previous fitting
					pf.flags &= ~(d ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
					pf.flags |= RO_FLAG_HANGING; // flat ends
					objs.push_back(pf);
					if (add_insul) {add_insul_exclude(p, pf, insul_exclude);}
				}
				else { // connector or exit pipe entering the wall or floor
					objs.back().flags |= RO_FLAG_HANGING; // flat ends
				}
			}
		} // for d
		if (p.type == PIPE_CONN && pipe.d[p.dim][0] < centerline-fitting_len && pipe.d[p.dim][1] > centerline+fitting_len) { // crosses through main pipe, add fitting
			room_object_t pf(pipe);
			pf.flags |= RO_FLAG_NOCOLL;
			pf.color  = fittings_color;
			set_wall_width(pf, centerline, fitting_len, p.dim);
			expand_cube_except_in_dim(pf, fitting_expand, p.dim); // expand slightly
			objs.push_back(pf);
		}
		if (p.dim < 2 && !add_insul) { // horizontal pipes only; no fittings on insulated pipes as insulation is flush with the walls
			// basement walls, parking garage walls, parking garage pillars, and beams
			add_pass_through_fittings(pipe, walls, fitting_len, fitting_expand, p.dim, LT_GRAY, objs); // different from fittings color
			add_pass_through_fittings(pipe, beams, fitting_len, fitting_expand, p.dim, LT_GRAY, objs);
		}
		if (p.type == PIPE_EXTB) { // add pipe hanging brackets
			float const ceiling_zval(interior->get_extb_start_room().z2() - fc_thickness);
			float const pipe_len(pipe.get_sz_dim(p.dim)), bracket_radius(radius_factor*p.radius); // includes insulation
			unsigned const num_brackets(0.5*pipe_len/window_vspacing);
			float const bracket_spacing(pipe_len/(num_brackets+1));

			for (unsigned n = 0; n < num_brackets; ++n) {
				float const len_pos(pipe.d[p.dim][0] + (n+1)*bracket_spacing);
				// skip checking of walls and obstacles since they shouldn't intersect; pass in objs end as start index to avoid checking pipes;
				// may intersect stacked pipes, but that should be okay; increase length to cover this up better
				add_hanging_pipe_bracket(pipe, len_pos, ceiling_zval, bracket_radius, 1.2, p.dim, room_id, tot_light_amt, objs.size(), objs, vect_cube_t(), vect_cube_t());
			}
		}
	} // for p
	for (pipe_t const &p : fittings) {
		cube_t const pbc(p.get_bcube());
		if (p.type != PIPE_EXTB && !cube_intersects_basement_or_extb_room(pbc)) continue; // outside the basement, don't need to draw
		bool const pdim(p.dim & 1), pdir(p.dim >> 1);
		unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI); // non-colliding, flat ends on both sides
		objs.emplace_back(pbc, TYPE_PIPE, room_id, pdim, pdir, flags, tot_light_amt, SHAPE_CYLIN, fittings_color);
		if (add_insul) {add_insul_exclude(p, pbc, insul_exclude);}
	} // for p
	if (add_insul) { // add hot water pipe insulation
		colorRGBA const insul_color(WHITE);
		vect_cube_t insulation, temp;

		for (pipe_t const &p : pipes) {
			if (!p.connected || p.outside_pg || p.type == PIPE_RISER) continue; // unconnected, outside, or drain, skip
			float const min_len(min_insum_len*p.radius), radius_exp(min(insul_thickness*p.radius, 0.2f*max_pipe_radius));
			if (p.get_length() < min_len) continue; // length is too short
			cube_t const pbc(p.get_bcube());
			insulation.clear();
			subtract_cubes_from_cube_split_in_dim(pbc, insul_exclude, insulation, temp, p.dim);
			unsigned const d1((p.dim+1)%3), d2((p.dim+2)%3);
			bool const pdim(p.dim & 1), pdir(p.dim >> 1); // encoded as: X:dim=0,dir=0 Y:dim=1,dir=0, Z:dim=x,dir=1
			unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_ADJ_HI | RO_FLAG_ADJ_LO); // not collidable, always flat ends
			if (p.type == PIPE_EXIT && p.dim == 2) {objs.back().flags |= RO_FLAG_LIT;} // vertical exit pipes are shadow casting
			
			for (cube_t &i : insulation) {
				if (i.d[d1][0] != pbc.d[d1][0] || i.d[d1][1] != pbc.d[d1][1] || i.d[d2][0] != pbc.d[d2][0] || i.d[d2][1] != pbc.d[d2][1]) continue; // clipped on any side
				if (i.get_sz_dim(p.dim) < min_len) continue; // too short
				i.expand_in_dim(d1, radius_exp);
				i.expand_in_dim(d2, radius_exp);
				objs.emplace_back(i, TYPE_PIPE, room_id, pdim, pdir, flags, tot_light_amt, SHAPE_CYLIN, insul_color);
			}
		} // for p
	}
	if ((pipe_type == PIPE_TYPE_SEWER || pipe_type == PIPE_TYPE_GAS) && water_damage > 0.25) {make_pipes_dirty(objs, pipes_start);}
	return 1;
}

void building_t::add_ext_basement_hallway_pipes_recur(unsigned room_id, bool hall_dim, unsigned pipe_type, float radius_factor,
	pipe_t const &parent, vector<pipe_t> &pipes, vector<pipe_t> &fittings, rand_gen_t &rgen) const
{
	// Note: pipes added here may loop around and meet each other or intersect other basement pipes if the extended basement hallway runs under the building;
	// this should be okay because pipes are at the same zval, so they connect at reasonable looking junctions (though without fittings)
	assert(has_ext_basement());
	room_t const &room(interior->get_room(room_id));
	float const room_centerline(room.get_center_dim(!hall_dim));
	bool const hall_side(room_centerline < parent.p1[!hall_dim]);
	float const rmin(((pipe_type == PIPE_TYPE_SEWER) ? 1.0 : 0.5)*0.01*get_window_vspace());

	// find connecting doors; start with first door after extended basement entrance door
	for (auto ds = interior->door_stacks.begin()+interior->ext_basement_door_stack_ix+1; ds != interior->door_stacks.end(); ++ds) {
		assert(ds->num_doors == 1); // must be a single door stack
		door_t const &door(get_door(ds->first_door_ix));
		if (door.dim == hall_dim || !door.is_connected_to_room(room_id)) continue; // not connected to the side of this hallway
		unsigned const conn_room_id(door.get_conn_room(room_id));
		room_t const &conn_room(interior->get_room(conn_room_id));
		if (!conn_room.is_hallway) continue;
		bool const conn_dir(room_centerline < door.get_center_dim(!hall_dim));
		if (conn_dir != hall_side) continue; // wrong side of the hallway - might intersect a ceiling light
		bool const door_side(rgen.rand_bool()); // add a branch to a random side of the door
		cube_t const room_bounds(get_walkable_room_bounds(conn_room));
		pipe_t conn_pipe(parent);
		// 0.8x = ~0.5^3; limit to between parent pipe radius and base pipe radius
		conn_pipe.radius = min(parent.radius, max(0.8f*parent.radius, rmin));
		float const wall_dist(get_pipe_dist_to_wall(radius_factor*conn_pipe.radius, get_trim_thickness()));
		point conn_pt(parent.p1);
		conn_pt[hall_dim] = room_bounds.d[hall_dim][door_side] + (door_side ? -1.0 : 1.0)*wall_dist; // adjacent to wall on side of door closest to source
		conn_pipe.p1 = conn_pipe.p2 = conn_pt;
		conn_pipe.p2[!hall_dim] = room_bounds.d[!hall_dim][conn_dir]; // far end of connecting hall
		conn_pipe.dim      = !hall_dim;
		conn_pipe.conn_dir = !conn_dir;
		if (conn_room.has_tunnel_conn()) {} // extend pipe through tunnel? code appears to be unreachable
		pipes.push_back(conn_pipe);
		// add fitting at the junction
		float const fitting_len(FITTING_LEN*conn_pipe.radius);
		point fitting_p1(conn_pt), fitting_p2(conn_pt);
		fitting_p1[hall_dim] -= fitting_len; fitting_p2[hall_dim] += fitting_len;
		fittings.emplace_back(fitting_p1, fitting_p2, FITTING_RADIUS*parent.radius, hall_dim, PIPE_FITTING, 3);
		add_ext_basement_hallway_pipes_recur(conn_room_id, !hall_dim, pipe_type, radius_factor, conn_pipe, pipes, fittings, rgen);
	} // for ds
}

// return value: 0=failed to place, 1=placed full length, 2=placed partial length
// Note: beams are for the top floor, even though we may be routing sprinkler pipes on a floor below
int add_sprinkler_pipe(building_t const &b, point const &p1, float end_val, float radius, bool dim, bool dir, vect_cube_t const &obstacles,
	vect_cube_t const &walls, vect_cube_t const &beams, vect_cube_t const &pipe_cubes, cube_t &ramp, float ceiling_zval, unsigned room_id,
	float tot_light_amt, vect_room_object_t &objs, colorRGBA const &pcolor, colorRGBA const &ccolor, unsigned max_shorten_steps, int add_sprinklers)
{
	point p2(p1);
	p2[dim] = end_val;
	cube_t h_pipe(p1, p2);
	h_pipe.expand_in_dim(!dim, radius);
	h_pipe.expand_in_dim(2,    radius); // set zvals
	bool const in_basement(ceiling_zval <= b.ground_floor_z1);
	unsigned const num_steps(8), max_n_iter(min(num_steps-1, max_shorten_steps+1));
	unsigned const objs_start(objs.size()), pipe_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING);
	unsigned conn_flags(pipe_flags);
	float &pipe_end(h_pipe.d[dim][!dir]);
	float const step_val((p1[dim] - pipe_end)/num_steps); // signed
	int ret(0);

	// keep cutting off the ends until we can place a pipe
	for (unsigned n = 0; n < max_n_iter; ++n, pipe_end += step_val) {
		// Note: pipe should touch the bottom of beams (or lower), so we don't need to check intersections with beams (and it may fail due to FP error)
		if (has_bcube_int(h_pipe, obstacles) || has_bcube_int(h_pipe, pipe_cubes)) continue; // no need to check walls here
		if (!ramp.is_all_zeros() && h_pipe.intersects(ramp)) continue; // check ramps as well, since they won't be included for lower floors
		if (n > 0) {conn_flags |= (dir ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI);} // shortened pipe, draw the connector end caps
		ret = ((n > 0) ? 2 : 1);
		break;
	} // for n
	if (ret == 0) return 0; // failed to place a pipe at any length
	bool const is_partial(ret == 2);
	// encoded as: X:dim=0,dir=0 Y:dim=1,dir=0, Z:dim=x,dir=1
	room_object_t const pipe(h_pipe, TYPE_PIPE, room_id, dim, 0, pipe_flags, tot_light_amt, SHAPE_CYLIN, pcolor); // add to pipe_cubes?
	objs.push_back(pipe);
	float const conn_thickness(0.2*radius), conn_max_length(3.2*radius);

	for (unsigned d = 0; d < 2; ++d) { // add connector segments
		float const conn_length(conn_max_length*((bool(d) == dir) ? 1.0 : (is_partial ? 0.4 : 0.6))); // long connected to pipe, medium at the wall, and short partial
		cube_t conn(h_pipe);
		conn.d[dim][!d] = conn.d[dim][d] + (d ? -1.0 : 1.0)*conn_length; // set length
		conn.expand_in_dim(!dim, conn_thickness);
		conn.expand_in_dim(2,    conn_thickness);
		unsigned flags(conn_flags);
		if (in_basement) {flags |= (d ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI);} // inside ends are visible
		else {flags |= (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI);} // above ground; both sides may be visible through a window
		objs.emplace_back(conn, TYPE_PIPE, room_id, dim, 0, flags, tot_light_amt, SHAPE_CYLIN, ccolor);
	} // for d
	if (is_partial) { // partial pipe, add a bracket to suspend the end
		float const len_pos(pipe_end - (dir ? -1.0 : 1.0)*2.6*conn_max_length);
		add_hanging_pipe_bracket(h_pipe, len_pos, ceiling_zval, radius, 0.8, dim, room_id, tot_light_amt, objs_start+1, objs, obstacles, walls);
	}
	if (add_sprinklers) {
		bool const inverted(add_sprinklers & 1); // use LSB
		float const length(h_pipe.get_sz_dim(dim)), sprinkler_radius(0.9*radius);
		unsigned const num_sprinklers(max(1U, unsigned(length/(60.0*radius))));
		float const end_spacing(is_partial ? 1.5*conn_max_length : length/(num_sprinklers - 0.5)); // near the connector if truncated, otherwise one spacing from the wall
		float const end(pipe_end - (dir ? -1.0 : 1.0)*end_spacing);
		float const spacing((num_sprinklers == 1) ? (pipe_end - p1[dim]) : (end - p1[dim])/(num_sprinklers - 0.5)); // signed; always halfway if single sprinkler
		float const sprinkler_height(3.2*radius), sprinkler_dz((inverted ? -1.0 : 1.0)*sprinkler_height);
		point center(p1);
		center[dim] = p1[dim] + 0.5*spacing; // half spacing from the end so that sprinklers are centered on the main pipe
		bool any_added(0);

		for (unsigned n = 0; n < num_sprinklers; ++n, center[dim] += spacing) {
			cube_t s(center);
			s.d[2][!inverted] += sprinkler_dz;
			s.expand_by_xy(sprinkler_radius);
			cube_t s_ext(s); // includes some extra clearance
			s_ext.d[2][!inverted] += sprinkler_dz;
			s.expand_by_xy(0.5*sprinkler_radius);
			if (has_bcube_int(s_ext, obstacles) || has_bcube_int(s_ext, walls) || has_bcube_int_xy(s_ext, beams) || has_bcube_int(s_ext, pipe_cubes)) continue;
			if (b.interior->is_blocked_by_stairs_or_elevator(s_ext)) continue; // check extra padding in front of stairs and elevators
			objs.emplace_back(s, TYPE_SPRINKLER, room_id, 0, inverted, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, pcolor);
			// add the connector ring around the sprinkler, same color as the pipe
			cube_t conn(h_pipe);
			set_wall_width(conn, center[dim], 0.35*conn_max_length, dim);
			conn.expand_in_dim(!dim, 0.8*conn_thickness);
			conn.expand_in_dim(2,    0.8*conn_thickness);
			objs.emplace_back(conn, TYPE_PIPE, room_id, dim, 0, (pipe_flags | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, pcolor);
			any_added = 1;
		} // for n
		// if no sprinklers were added (maybe because the pipe was along a beam), remove the entire pipe
		if (!any_added) {objs.resize(objs_start); return 0;}
	}
	if (h_pipe.z2() < ceiling_zval) { // add hangers for each beam the pipe passes under; ignores beam zval (assumes beams stack on each floor)
		for (cube_t const &beam : beams) {
			if (beam.get_sz_dim(!dim) < beam.get_sz_dim(dim)) continue; // beam runs in the wrong dim (parallel, not perpendicular)
			if (beam.d[!dim][0] > h_pipe.d[!dim][0] || beam.d[!dim][1] < h_pipe.d[!dim][1]) continue; // beam length doesn't contain pipe
			if (beam.d[ dim][0] < h_pipe.d[ dim][0] || beam.d[ dim][1] > h_pipe.d[ dim][1]) continue; // pipe length doesn't contain beam
			float const len_pos(beam.get_center_dim(dim));
			add_hanging_pipe_bracket(h_pipe, len_pos, ceiling_zval, radius, 0.8, dim, room_id, tot_light_amt, objs_start+1, objs, obstacles, walls);
		}
	}
	// add fittings at walls and pillars
	float const fitting_len(FITTING_LEN*radius), fitting_expand((FITTING_RADIUS - 1.0)*radius);
	add_pass_through_fittings(pipe, walls, fitting_len, fitting_expand, dim, ccolor, objs);
	return ret;
}
bool building_t::add_sprinkler_pipes(vect_cube_t const &obstacles, vect_cube_t const &walls, vect_cube_t const &beams, vect_cube_t const &pipe_cubes, unsigned room_id,
	unsigned num_floors, unsigned objs_start, rand_gen_t &rgen, float custom_floor_spacing, float wall_pad, unsigned pref_dim, vect_cube_t const &vpipe_avoid)
{
	// add vertical red (possibly rusted brown) sprinkler system pipe
	cube_t room(get_room(room_id));
	bool const in_basement(room.z1() < ground_floor_z1); // else in factory
	float const tot_light_amt = 1.0; // to offset the darkness of the basement
	float const window_vspace(get_window_vspace()), floor_spacing((custom_floor_spacing > 0.0) ? custom_floor_spacing : window_vspace);
	float const fc_thickness(get_fc_thickness()), wall_thickness(get_wall_thickness());
	float const sp_radius((in_basement ? 1.2 : 0.9)*wall_thickness), spacing(2.0*sp_radius), flange_expand(0.3*sp_radius); // larger for basement, smaller for factory
	float const bolt_dist(sp_radius + 0.5*flange_expand), bolt_radius(0.32*flange_expand), bolt_height(0.1*fc_thickness);
	bool const inverted_sprinklers(room_id & 1); // random-ish
	// pipe color; fade to rusty brown for basement pipes when there's water damage
	colorRGBA const pcolor(in_basement ? blend_color(DK_BROWN, RED, min(1.0f, 2.0f*water_damage), 0) : RED);
	colorRGBA const ccolor(BRASS_C); // connector color
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const sprinker_pipes_start(objs.size());
	room.expand_by_xy(-wall_pad);
	assert(room.is_strictly_normalized());
	cube_t c;
	set_cube_zvals(c, (room.z1() + fc_thickness), (room.z2() - fc_thickness));

	for (unsigned n = 0; n < 100; ++n) { // 100 random tries
		bool const dim((pref_dim < 2 && n < 50) ? bool(pref_dim) : rgen.rand_bool()), dir(rgen.rand_bool()); // wall dim/dir
		set_wall_width(c, rgen.rand_uniform(room.d[!dim][0]+spacing, room.d[!dim][1]-spacing), sp_radius, !dim);
		c.d[dim][ dir] = room.d[dim][dir] + (dir ? -1.0 : 1.0)*flange_expand; // against the wall (with space for the flange)
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*2.0*sp_radius;
		cube_t c2(c);
		c2.expand_in_dim(!dim, 0.5*sp_radius); // add a bit of extra space to the sides for the flanges and valves
		if (has_bcube_int(c2, obstacles) || has_bcube_int(c2, walls) || has_bcube_int(c2, beams) || has_bcube_int(c2, pipe_cubes) || has_bcube_int(c2, vpipe_avoid)) continue;
		// skip if the pipe aligns with a pillar because we won't be able to place a horizontal sprinkler pipe
		bool is_blocked(0);
		
		for (cube_t const &wall : walls) { // includes both walls and pillars
			if (wall.dx() > 2.0*wall.dy() || wall.dy() > 2.0*wall.dx()) continue; // skip walls since they have pillars intersecting them as well
			if ((wall.x2() > c.x1() && wall.x1() < c.x2()) || (wall.y2() > c.y1() && wall.y1() < c.y2())) {is_blocked = 1; break;} // X or Y projection
		}
		if (is_blocked) continue;
		objs.emplace_back(c, TYPE_PIPE, room_id, 0, 1, RO_FLAG_LIT, tot_light_amt, SHAPE_CYLIN, pcolor); // dir=1 for vertical; casts shadows; add to pipe_cubes?
		cube_t flange(c);
		flange.expand_by_xy(flange_expand);
		bool add_bot_flange(0);

		if (!in_basement && has_parking_garage) { // extend into basement parking garage?
			cube_t const &basement(get_basement());

			if (basement.contains_cube_xy(c)) { // always true?
				cube_t c_ext(c);
				set_cube_zvals(c_ext, (basement.z1() + fc_thickness), (basement.z2() - fc_thickness));

				if (!is_obj_placement_blocked(c_ext, basement, 1)) {
					objs.emplace_back(objs.back()); // duplicate the vertical pipe and adjust it's zvals
					objs.back().copy_from(c_ext);
					// add flange (no bolts on bottom)
					flange.z1() = basement.z1() + 1.00*fc_thickness;
					flange.z2() = basement.z1() + 1.15*fc_thickness;
					objs.emplace_back(flange, TYPE_PIPE, room_id, 0, 1, (RO_FLAG_HANGING | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, pcolor);
					add_bot_flange = 1;
				}
			}
		}
		// add flanges at top and bottom of each floor
		point const center(c.get_cube_center());

		for (unsigned f = 0; f <= num_floors; ++f) { // flanges for each ceiling/floor
			bool const at_bot(!add_bot_flange && f == 0), at_top(f == num_floors);
			unsigned flags(RO_FLAG_HANGING | (!at_bot)*RO_FLAG_ADJ_LO | (!at_top || !in_basement)*RO_FLAG_ADJ_HI);
			float const z(room.z1() + f*floor_spacing);
			flange.z1() = z - (at_bot ? -1.0 : 1.15)*fc_thickness;
			flange.z2() = z + (at_top ? -1.0 : 1.15)*fc_thickness;
			objs.emplace_back(flange, TYPE_PIPE, room_id, 0, 1, flags, tot_light_amt, SHAPE_CYLIN, pcolor);
			unsigned const NUM_BOLTS = 8;
			float const angle_step(TWO_PI/NUM_BOLTS);

			for (unsigned m = 0; m < NUM_BOLTS; ++m) { // add bolts
				float const angle(m*angle_step), dx(bolt_dist*sin(angle)), dy(bolt_dist*cos(angle));
				cube_t bolt;
				bolt.set_from_sphere(point(center.x+dx, center.y+dy, 0.0), bolt_radius);
				set_cube_zvals(bolt, flange.z1()-bolt_height, flange.z2()+bolt_height);
				objs.emplace_back(bolt, TYPE_PIPE, room_id, 0, 1, (flags | RO_FLAG_RAND_ROT), tot_light_amt, SHAPE_CYLIN, pcolor);
			} // for m
		} // for f
		// add valves on each floor
		for (unsigned f = 0; f < num_floors; ++f) {
			float const z(room.z1() + (f + 0.55)*window_vspace), valve_radius(1.5*sp_radius);
			point center(c.xc(), c.yc(), z);
			center[dim] += (dir ? -1.0 : 1.0)*1.5*sp_radius;
			cube_t valve(center);
			valve.expand_in_dim(2,        valve_radius);
			valve.expand_in_dim(!dim,     valve_radius);
			valve.expand_in_dim( dim, 0.4*valve_radius);
			if (has_bcube_int(valve, obstacles) || has_bcube_int(valve, walls)) continue; // check for pillars, etc.; should be rare
			objs.emplace_back(valve, TYPE_VALVE, room_id, dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, pcolor);
			// add a brass band around the pipe where the valve connects
			cube_t band(c);
			band.expand_by_xy(0.1*sp_radius);
			set_cube_zvals(band, z-0.5*valve_radius, z+0.5*valve_radius);
			objs.emplace_back(band, TYPE_PIPE, room_id, 0, 1, (RO_FLAG_HANGING | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI), tot_light_amt, SHAPE_CYLIN, BRASS_C);
		} // for f
		// attempt to run horizontal pipes across the room ceiling
		float const h_pipe_radius(0.6*wall_thickness), conn_thickness(0.2*h_pipe_radius);
		float const ceil_gap(max(0.25f*fc_thickness, 0.05f*get_floor_ceil_gap())); // make enough room for both flange + bolts and ceiling beams
		unsigned const max_shorten_steps(n/5); // start off with no horizontal pipe shortening to get better coverage, but relax this constraint as we go
		bool place_failed(0);

		for (unsigned f = 0; f < num_floors; ++f) {
			bool const lf(f+1 < num_floors);
			vect_cube_t walls2, beams2, obstacles2;

			if (lf) { // query walls/beams/obstacles added to lower floors in previous PG placement steps
				obstacles2 = obstacles; // copy existing obstacles, since the stairs, elevators, ramps, etc. that aren't objects or per-floor are needed
				add_pg_obstacles(objs, interior->room_geom->wall_ps_start, objs_start, walls2, beams2, obstacles2);
			}
			vect_cube_t const &walls_(lf ? walls2 : walls), &beams_(lf ? beams2 : beams), &obstacles_(lf ? obstacles2 : obstacles);
			float const ceiling_zval(room.z1() + (f+1)*floor_spacing - fc_thickness);
			point p1(center); // vertical pipe center
			// place below the ceiling so that fittings just touch beams rather than clipping through them; however, pipes may still clip through lights placed on beams
			p1.z = ceiling_zval - ceil_gap - h_pipe_radius - conn_thickness;
			unsigned const pipe_obj_ix(objs.size());
			float const wpos(room.d[dim][!dir]); // extend to the opposite wall
			int const ret(add_sprinkler_pipe(*this, p1, wpos, h_pipe_radius, dim, dir, obstacles_, walls_, beams_, pipe_cubes,
				interior->pg_ramp, ceiling_zval, room_id, tot_light_amt, objs, pcolor, ccolor, max_shorten_steps, 0)); // sprinklers=0
			
			if (ret == 0) { // failed to place
				if (n < 50) {place_failed = 1; break;} // failed; retry with a different vertical pipe placement
				// try to run horizontal pipe in the opposite dim, but don't connect branch lines because the pipe may be too close to the wall and the code would be messy
				bool const other_dir(room.get_center_dim(!dim) < p1[!dim]);
				add_sprinkler_pipe(*this, p1, room.d[!dim][!other_dir], h_pipe_radius, !dim, other_dir, obstacles_, walls_, beams_, pipe_cubes,
					interior->pg_ramp, ceiling_zval, room_id, tot_light_amt, objs, pcolor, ccolor, max_shorten_steps, 0); // sprinklers=0
				continue;
			}
			// run smaller branch lines off this pipe in the other dim; we could add the actual sprinklers to these
			cube_t const h_pipe(objs[pipe_obj_ix]);
			float const h_pipe_len(h_pipe.get_sz_dim(dim)), conn_radius(0.5*h_pipe_radius);
			unsigned const pri_pipe_end_ix(objs.size()), num_conn(max(1U, unsigned(h_pipe_len/(1.5*window_vspace))));
			bool const end_flush(ret == 2); // if a partial pipe was added, place the last connector at the end of the pipe so that there's no unused length
			float const conn_step((dir ? -1.0 : 1.0)*h_pipe_len/(num_conn - (end_flush ? 0.5 : 0.0))); // signed step (pipe runs in -dir)
			p1[dim] += 0.5*conn_step; // first half step
			int added(0);
			float last_added_conn_pipe_pos(center[dim]); // start at center of vertical pipe

			for (unsigned n = 0; n < num_conn; ++n) { // add secondary connector pipes, with sprinklers
				added = 0; // reset for this conn

				for (unsigned d = 0; d < 2; ++d) { // extend to either side of the pipe
					float const wpos2(room.d[!dim][!d]); // extend to the opposite wall
					added |= add_sprinkler_pipe(*this, p1, wpos2, conn_radius, !dim, d, obstacles_, walls_, beams_, pipe_cubes,
						interior->pg_ramp, ceiling_zval, room_id, tot_light_amt, objs, pcolor, ccolor, 8, (inverted_sprinklers ? 1 : 2)); // add sprinklers
				}
				if (added) { // if conn was added in either dir, add a connector segment
					unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI); // flat ends on both sides
					cube_t conn(h_pipe);
					set_wall_width(conn, p1[dim], 1.2*h_pipe_radius, dim); // set length
					conn.expand_in_dim(!dim, conn_thickness);
					conn.expand_in_dim(2,    conn_thickness);
					objs.emplace_back(conn, TYPE_PIPE, room_id, dim, 0, flags, tot_light_amt, SHAPE_CYLIN, ccolor);
					last_added_conn_pipe_pos = p1[dim];
				}
				p1[dim] += conn_step;
			} // for n
			if (end_flush && !added && last_added_conn_pipe_pos != center[dim]) { // partial pipe with final connector not added, but some other connector added
				cube_t &pipe_bc(objs[pipe_obj_ix]);
				float &end_to_move(pipe_bc.d[dim][!dir]);
				float const pipe_end(end_to_move), xlate(last_added_conn_pipe_pos - pipe_end);
				end_to_move = last_added_conn_pipe_pos;

				for (unsigned i = pipe_obj_ix+1; i < pipe_obj_ix+3; ++i) { // check both end caps
					if (objs[i].d[dim][!dir] == pipe_end) {objs[i].translate_dim(dim, xlate);}
				}
				for (unsigned i = pipe_obj_ix+3; i < pri_pipe_end_ix; ++i) { // remove end caps no longer covering shortened pipe
					room_object_t &obj(objs[i]);
					if (obj.d[dim][0] < pipe_bc.d[dim][0] || obj.d[dim][1] > pipe_bc.d[dim][1]) {obj.remove();}
				}
			}
		} // for f
		if (place_failed) { // placement of secondary pipes failed; remove all sprinkler pipes and retry
			objs.resize(sprinker_pipes_start);
			continue;
		}
		if (in_basement && water_damage > 0.25) {make_pipes_dirty(objs, sprinker_pipes_start);} // make all pipes rusty
		return 1; // done
	} // for n
	return 0; // failed
}

// here each sphere represents the entry point of a pipe with this radius into the basement ceiling
// find all plumbing fixtures such as toilets, urinals, sinks, and showers; these should have all been placed in rooms by now
void building_t::get_pipe_basement_water_connections(vect_riser_pos_t &sewer, vect_riser_pos_t &cold_water, vect_riser_pos_t &hot_water, rand_gen_t &rgen) const {
	assert(has_room_geom());
	cube_t const &basement(get_basement());
	float const merge_dist = 4.0; // merge two pipes if their combined radius is within this distance
	// use reduced pipe radius for apartments and hotels since they have so many plumbing fixtures
	float const floor_spacing(get_window_vspace()), floor_thickness(get_floor_thickness()), second_floor_zval(ground_floor_z1 + floor_spacing);
	float const base_pipe_radius((is_apt_or_hotel() ? 0.008 : 0.01)*floor_spacing), base_pipe_area(base_pipe_radius*base_pipe_radius);
	float const merge_dist_sq(merge_dist*merge_dist), max_radius(0.3*get_wall_thickness()), ceil_zval(basement.z2() - get_fc_thickness());
	bool const inc_extb_conns(has_ext_basement() && !has_backrooms_or_mall());
	bool const inc_basement_wheaters = 1;
	float extb_pipe_radius(has_pool() ? base_pipe_radius : 0.0); // pools use cold water
	bool extb_pipe_has_hot(0);
	cube_t extb_area;
	vect_room_object_t water_heaters;

	// start with sewer pipes and water heaters
	for (room_object_t const &i : interior->room_geom->objs) { // check all objects placed so far
		if (i.in_mall()) continue; // skip appliance/plumbing store objects

		if (i.type == TYPE_WHEATER) { // water heaters are special because they take cold water and return hot water
			// maybe skip if in the basement, since this must connect directly to pipes rather than through a riser;
			// this can happen for houses (which don't have parking garages or pipes), but currently not for office buildings;
			// also skip water heaters in the extended basement (not inside basement)
			if (i.z1() < ground_floor_z1 && (!inc_basement_wheaters || !basement.intersects(i))) continue;
			water_heaters.push_back(i);
			continue;
		}
		// Note: the dishwasher is always next to the kitchen sink and uses the same water connections
		bool const hot_cold_obj (i.type == TYPE_SINK || i.type == TYPE_BRSINK || i.type == TYPE_KSINK || i.type == TYPE_TUB ||
			i.type == TYPE_SHOWER || i.type == TYPE_SHOWERTUB || i.type == TYPE_WASHER || i.type == TYPE_DWASHER || i.type == TYPE_VANITY);
		bool const cold_only_obj(i.type == TYPE_TOILET || i.type == TYPE_URINAL || i.type == TYPE_DRAIN);
		if (!hot_cold_obj && !cold_only_obj) continue;

		if (i.z1() < ground_floor_z1) { // object in the basement
			if (inc_extb_conns && point_in_extended_basement_not_basement(i.get_cube_center())) {
				extb_pipe_radius   = get_merged_pipe_radius(extb_pipe_radius, base_pipe_radius, 3.0);
				extb_pipe_has_hot |= hot_cold_obj;
			}
			continue; // unclear how to handle it here - should be a direct connection or go through the wall/floor
		}
		bool const upper_floor(i.z1() > second_floor_zval);
		point pos(i.xc(), i.yc(), ceil_zval);
		//if (!basement.contains_pt_xy(pos)) continue; // riser/pipe doesn't pass through the basement, but this is now allowed
		bool merged(0);

		// see if we can merge this sewer riser into an existing nearby riser
		for (auto &r : sewer) {
			float const p_area(r.radius*r.radius), sum_area(p_area + base_pipe_area);
			if (p2p_dist_xy_sq(r.pos, pos) > merge_dist_sq*sum_area) continue;
			r.pos      = (p_area*r.pos + base_pipe_area*pos)/sum_area; // merged position is weighted average area
			r.radius   = get_merged_pipe_radius(r.radius, base_pipe_radius, 3.0); // cubic
			r.has_hot |= hot_cold_obj;
			r.upper_floor |= upper_floor;
			merged     = 1;
			break;
		} // for p
		if (!merged) {sewer.emplace_back(pos, base_pipe_radius, hot_cold_obj, 0, upper_floor);} // create a new riser; flow_dir=0/out
	} // for i
	if (inc_extb_conns /*&& extb_pipe_radius > 0.0*/) { // we have extended basement water consumers (or maybe enable even if not), add a connection through the hallway
		extb_pipe_radius = max(min(0.5f*extb_pipe_radius, 2.0f*base_pipe_radius), base_pipe_radius); // not too large and not too small
		room_t const &hallway(interior->get_extb_start_room());
		extb_area = get_walkable_room_bounds(hallway);
		bool const wdim(interior->extb_wall_dim), wdir(interior->extb_wall_dir), first_side(rgen.rand_bool());
		float const wall_dist(get_pipe_dist_to_wall(extb_pipe_radius, get_trim_thickness())); // calculated for sewer pipe, but should be large enough for water/HW pipes
		point conn_pt(0.0, 0.0, ceil_zval);
		conn_pt[wdim] = extb_area.d[wdim][wdir]; // at the far end of the hallway

		for (unsigned side = 0; side < 2; ++side) { // try both sides, in the hope that we can connect to one
			bool const d(bool(side) ^ first_side);
			conn_pt[!wdim] = extb_area.d[!wdim][d] - (d ? 1.0 : -1.0)*wall_dist; // in the corner, but avoiding the door trim
			sewer.emplace_back(conn_pt, extb_pipe_radius, extb_pipe_has_hot, 0, 0, 1); // flow_dir=0/out, upper_floor=0, in_extb=1
		}
	}
	vect_riser_pos_t water_risers;
	water_risers.reserve(sewer.size());

	for (riser_pos_t &s : sewer) {
		min_eq(s.radius, max_radius); // clamp radius to a reasonable value after all merges
		water_risers.push_back(s);
		// try to find a nearby interior wall on the ground floor to route the riser to
		if (!basement.contains_pt_xy(s.pos)) continue; // outside the basement; riser not drawn anyway
		float dmin(0.5*floor_spacing); // max movement distance
		// if we're close to the basement wall, then leave pos unchanged; we don't want to move to the basement wall and have it half inside the basement
		for (unsigned d = 0; d < 2; ++d) {min_eq(dmin, min((basement.d[d][1] - s.pos[d]), (s.pos[d] - basement.d[d][0])));}
		vector3d best_delta;

		for (unsigned d = 0; d < 2; ++d) { // find nearby interior walls; d is the wall separating dim
			float const lo(s.pos[!d] - s.radius), hi(s.pos[!d] + s.radius);

			for (cube_t const &w : interior->walls[d]) {
				if (w.z1() > ground_floor_z1 + floor_thickness || w.z2() < ground_floor_z1 - floor_thickness) continue; // not on the ground floor
				if (w.d[!d][0] > lo || w.d[!d][1] < hi) continue; // riser not contained in wall length
				float const center(w.get_center_dim(d)), delta(center - s.pos[d]), dist(fabs(delta));
				if (dist > dmin) continue;
				dmin = dist;
				best_delta[!d] = 0.0;
				best_delta[ d] = delta;
			}
		} // for d
		if (best_delta.x == 0.0 && best_delta.y == 0.0) continue; // no valid wall to move to
		// always try to move hot and cold water pipes to a wall since they usually connect through the wall anyway;
		// sewer/drain pipes on the first floor typically run directly through the floor (and may not fit in the wall), but must be routed through the wall on upper floors
		if (s.upper_floor) {s.pos += best_delta;}
		water_risers.back()  .pos += best_delta;
		assert(basement.contains_pt_xy(water_risers.back().pos));
	} // for s
	// generate hot and cold water pipes
	// choose a shift direction 45 degrees diagonally in the XY plane to avoid collisions with sewer pipes placed in either dim
	vector3d const shift_dir(vector3d((rgen.rand_bool() ? -1.0 : 1.0), (rgen.rand_bool() ? -1.0 : 1.0), 0.0).get_norm());
	cold_water.reserve(water_risers.size()); // one cold water pipe per sewer pipe

	for (riser_pos_t &s : water_risers) {
		float const wp_radius(0.5*s.radius), pipe_spacing(2.0*(s.radius + wp_radius));
		cube_t place_area(s.in_extb ? extb_area : (basement.contains_pt_xy(s.pos) ? basement : bcube)); // force water riser inside basement if sewer riser is there
		place_area.expand_by_xy(-wp_radius);
		vector3d delta(shift_dir*(s.in_extb ? 0.0 : pipe_spacing));
		if (!place_area.contains_pt_xy(s.pos + delta)) {delta.negate();} // if shift takes pos outside placement area, shift in the other direction
		cold_water.emplace_back((s.pos + delta), wp_radius, 0, 1, s.upper_floor, s.in_extb); // has_hot=0, flow_dir=1/in
		if (!s.has_hot) continue; // no hot water, done
		float const hot_radius(0.75*wp_radius); // smaller radius than cold water pipe
		place_area.expand_by_xy(wp_radius - hot_radius); // resize for new radius
		if (!place_area.contains_pt_xy(s.pos - delta)) {delta *= 2.0;} // if shift takes pos outside placement area, shift further from the drain
		hot_water.emplace_back((s.pos - delta), hot_radius, 1, 1, s.upper_floor, s.in_extb); // shift in opposite dir; has_hot=1, flow_dir=1/in
	} // for sewer
	if (!water_heaters.empty() && !hot_water.empty()) { // add connections for water heaters if there are hot water pipes
		float const radius_hot(get_merged_risers_radius(hot_water)); // this is the radius of the main hot water supply; not limited to max
		float const per_wh_radius(radius_hot*pow(1.0/water_heaters.size(), 1/4.0)); // distribute evenly among the water heaters using the same merge exponent

		for (room_object_t const &wh : water_heaters) {
			float const radius(wh.get_radius()), shift_val(WHEATER_PIPE_SPACING*radius);
			point center(wh.xc(), wh.yc(), ceil_zval);
			vector3d wh_shift_dir(shift_dir);

			if (wh.z1() < ground_floor_z1) { // house basement water heater
				// should this connect directly? what if vent or other pipe is in the way? how to match pipe radius?
			}
			else { // office building first floor water heater: use exact location where bent over pipes meet the floor
				center[wh.dim] += (wh.dir ? 1.0 : -1.0)*radius*WHEATER_PIPE_H_DIST;
				wh_shift_dir[ wh.dim] = 0.0;
				wh_shift_dir[!wh.dim] = ((shift_dir[!wh.dim] < 0.0) ? -1.0 : 1.0);
			}
			cold_water.emplace_back((center + shift_val*wh_shift_dir), per_wh_radius, 0, 1); // has_hot=0, flow_dir=1/in
			hot_water .emplace_back((center - shift_val*wh_shift_dir), per_wh_radius, 1, 0); // has_hot=1, flow_dir=0/out
		} // for wh
	}
	if (!is_house) { // add pipes for office building rooftop water tower, if present; unclear if we ever get here
		for (auto i = details.begin(); i != details.end(); ++i) {
			if (i->type != ROOF_OBJ_WTOWER) continue;
			// assume water enters water tower through this pipe and exits water tower to provide water to upper floors, which we don't need to consider in the basement
			cold_water.emplace_back(point(i->xc(), i->yc(), ceil_zval), 2.0*base_pipe_radius, 0, 1); // 2x base radius, has_hot=0, flow_dir=1/in
		}
	}
}

void building_t::get_pipe_basement_gas_connections(vect_riser_pos_t &pipes) const {
	float const pipe_radius(0.005*get_window_vspace()), ceil_zval(get_basement().z2() - get_fc_thickness());

	for (room_object_t const &i : interior->room_geom->objs) { // check all objects placed so far
		if (i.z1() < ground_floor_z1) continue; // object in the house basement; unclear how to handle it here
		bool const is_gas_dryer(i.type == TYPE_DRYER && (i.obj_id & 3)); // gas dryer 75% of the time, since it makes the pipes more interesting
		if (i.type != TYPE_WHEATER && i.type != TYPE_FURNACE && i.type != TYPE_STOVE && i.type != TYPE_FPLACE && i.type != TYPE_HVAC_UNIT && !is_gas_dryer) continue;
		pipes.emplace_back(point(i.xc(), i.yc(), ceil_zval), pipe_radius, 0, 1); // flows in
	}
}

void get_water_heater_cubes(room_object_t const &wh, cube_t cubes[2]) { // {tank, pipes}
	cubes[0] = cubes[1] = wh;
	// shrink to include the pipes
	float const radius(wh.get_radius());
	set_wall_width(cubes[1], wh.get_center_dim( wh.dim), 0.2*radius,  wh.dim); // width of vent
	set_wall_width(cubes[1], wh.get_center_dim(!wh.dim), 0.7*radius, !wh.dim); // width of top pipes extent
	cubes[0].z2() -= 0.2*wh.dz(); // shorten the height for the tank; needed for vertical exit pipes
}

void building_t::add_house_basement_pipes(rand_gen_t &rgen) {
	if (!has_room_geom()) return; // error?
	float const fc_thick(get_fc_thickness());
	cube_t const &basement(get_basement());
	// hang sewer pipes under the ceiling beams; hang water pipes from the ceiling, above sewer pipes and through the beams
	unsigned const num_floors(1); // basement is always a single floor
	// houses have smaller radius pipes, so we should have enough space to stack sewer below hot water below cold water
	float const ceil_z(basement.z2()), sewer_zval(ceil_z - 1.8*fc_thick), cw_zval(ceil_z - 1.0*fc_thick), hw_zval(ceil_z - 1.4*fc_thick), gas_zval(ceil_z - 2.8*fc_thick);
	float const trim_thickness(get_trim_thickness()), wall_thickness(get_wall_thickness());
	vect_cube_t pipe_cubes, obstacles, walls, beams; // beams remains empty
	unsigned room_id(0);
	if (has_ext_basement()) {obstacles.push_back(get_ext_basement_door_blocker());}

	// we can't pass in a single valid room_id because house basements/pipes span multiple rooms, but we can at least use the ID of a room in the basement
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (basement.contains_cube(*r)) {room_id = (r - interior->rooms.begin()); break;}
	}
	for (unsigned d = 0; d < 2; ++d) { // add all basement walls
		for (cube_t &wall : interior->walls[d]) {
			if (wall.z1() >= ground_floor_z1 || !basement.intersects(wall)) continue; // not in the basement
			walls.push_back(wall);
			walls.back().expand_by_xy(trim_thickness); // include the trim
		}
	}
	// Note: elevators/buttons/stairs haven't been placed at this point, so iterate over all objects
	for (room_object_t const &i : interior->room_geom->objs) {
		if (i.z1() >= ground_floor_z1 || !basement.intersects(i)) continue; // not in the basement
		bool const no_blocking(i.type == TYPE_PICTURE || i.type == TYPE_WBOARD || i.type == TYPE_OUTLET || i.type == TYPE_SWITCH);
		// Note: bottles and trash on the floor counts, though it would be better to add pipes first and floor clutter later
		// Note: TYPE_PIPE (vertical electrical conduits from outlets) may block pipes from running horizontally along walls
		if (i.no_coll() && !no_blocking && i.type != TYPE_LIGHT && i.type != TYPE_PIPE && i.type != TYPE_VENT && !i.is_floor_clutter()) continue; // no collisions
		// Note: we could maybe skip if i.z2() < sewer_zval-pipe_radius, but we still need to handle collisions with vertical exit pipe segments

		if (i.type == TYPE_VENT) {
			obstacles.push_back(i);
			obstacles.back().z1() -= fc_thick; // add some clearance under the vent
		}
		else if (i.type == TYPE_WHEATER) {
			cube_t cubes[2];
			get_water_heater_cubes(i, cubes);
			UNROLL_2X(obstacles.push_back(cubes[i_]);)
		}
		else if (no_blocking) {
			cube_t c(i);
			c.d[i.dim][i.dir] += (i.dir ? 1.0 : -1.0)*wall_thickness; // add a wall thickness of clearance
			if (i.type == TYPE_PICTURE) {c.expand_in_dim(!i.dim, 0.05*i.get_sz_dim(!i.dim));} // expand slightly to include the frame
			obstacles.push_back(c);
		}
		else {obstacles.push_back(i);}
	} // for i
	// add doors as obstacles; maybe should move the ceiling up or move the tops of the doors down to avoid door collisions?
	// but this would need to be handled inside add_basement_pipes(), which uses obstacles in 8+ places;
	// so we would probably need to reduce the height of all doors and add extra wall segments above them to make this work,
	// and this is too late in the interior creation process to add walls without regenerating the VBO (unless we use stairs walls?)
	float const door_trim_exp(2.0*trim_thickness + 0.5*wall_thickness);

	for (door_t const &d : interior->doors) {
		if (d.z1() >= ground_floor_z1) continue; // not in the basement
		door_t door(d);
		door.open = 0; door.open_amt = 0.0; // start closed
		cube_t door_bcube(get_door_bounding_cube(door));
		if (has_ext_basement() && !basement.intersects(door_bcube)) continue; // skip extended basement doors
		door_bcube.d[d.dim][ d.open_dir] += (d.open_dir ? 1.0 : -1.0)*d.get_width(); // include space for door to swing open
		door_bcube.d[d.dim][!d.open_dir] -= (d.open_dir ? 1.0 : -1.0)*door_trim_exp; // include door trim width on the other side

		if (!d.on_stairs) { // [basement] stairs doors don't really open, so we only need clearance in front
			door.open = 1; door.open_amt = 1.0; // now try the open door to avoid blocking it when open
			door_bcube.union_with_cube(get_door_bounding_cube(door));
		}
		obstacles.push_back(door_bcube);
	} // for doors
	for (stairwell_t const &s : interior->stairwells) { // add stairwells (basement stairs); there should be no elevators
		if (s.z1() >= ground_floor_z1 || !basement.intersects(s)) continue; // not in the basement
		obstacles.push_back(s);
	}
	unsigned const objs_start(interior->room_geom->objs.size()); // no other basement objects added here that would interfere with pipes
	vect_riser_pos_t sewer, cold_water, hot_water, gas_pipes;
	get_pipe_basement_water_connections(sewer, cold_water, hot_water, rgen);
	add_basement_pipes(obstacles, walls, beams, sewer,      pipe_cubes, room_id, num_floors, objs_start, sewer_zval, rgen, PIPE_TYPE_SEWER, 1); // sewer
	add_to_and_clear(pipe_cubes, obstacles); // add sewer pipes to obstacles
	add_basement_pipes(obstacles, walls, beams, cold_water, pipe_cubes, room_id, num_floors, objs_start, cw_zval,    rgen, PIPE_TYPE_CW,    1); // cold water
	add_to_and_clear(pipe_cubes, obstacles); // add cold water pipes to obstacles
	add_basement_pipes(obstacles, walls, beams, hot_water,  pipe_cubes, room_id, num_floors, objs_start, hw_zval,    rgen, PIPE_TYPE_HW,    1); // hot water
	add_to_and_clear(pipe_cubes, obstacles); // add hot water pipes to obstacles
	get_pipe_basement_gas_connections(gas_pipes);
	add_basement_pipes(obstacles, walls, beams, gas_pipes,  pipe_cubes, room_id, num_floors, objs_start, gas_zval,   rgen, PIPE_TYPE_GAS,   1); // gas
}


