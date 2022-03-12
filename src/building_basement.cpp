// 3D World - Building Basement and Parking Garage Logic
// by Frank Gennari 03/11/2022

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for car_t
#include "profiler.h"

extern city_params_t city_params; // for num_cars

car_t car_from_parking_space(room_object_t const &o);
void subtract_cube_from_floor_ceil(cube_t const &c, vect_cube_t &fs);
bool line_int_cubes_exp(point const &p1, point const &p2, vect_cube_t const &cubes, vector3d const &expand);


bool building_t::add_basement_utility_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const height(get_window_vspace() - get_floor_thickness()), radius(0.18*height);
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-(1.05*radius + get_trim_thickness())); // account for the pan
	vect_room_object_t &objs(interior->room_geom->objs);
	point center(0.0, 0.0, zval);
	bool was_placed(0);

	for (unsigned n = 0; n < 5; ++n) { // make 14 attempts to place a water heater - one in each corner and 1 along a random wall for variety
		bool const dim(rgen.rand_bool());
		bool dir(0);

		if (n < 4) { // corner
			bool const xdir(rgen.rand_bool()), ydir(rgen.rand_bool());
			dir = (dim ? ydir : xdir);
			center.x = place_area.d[0][xdir];
			center.y = place_area.d[1][ydir];
		}
		else { // wall
			dir = rgen.rand_bool(); // choose a random wall
			center[ dim] = place_area.d[dim][dir]; // against this wall
			center[!dim] = rgen.rand_uniform(place_area.d[!dim][0], place_area.d[!dim][1]);
		}
		cube_t const c(get_cube_height_radius(center, radius, height));
		if (is_cube_close_to_doorway(c, room, 0.0, !room.is_hallway) || interior->is_blocked_by_stairs_or_elevator(c)) continue;
		cube_t c_exp(c);
		c_exp.expand_by_xy(0.2*radius); // small keepout in XY
		c_exp.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.25*radius; // add more keepout in front where the controls are
		c_exp.intersect_with_cube(room); // don't pick up objects on the other side of the wall
		if (overlaps_other_room_obj(c_exp, objs_start)) continue; // check existing objects, in particular storage room boxes that will have already been placed
		objs.emplace_back(c, TYPE_WHEATER, room_id, dim, !dir, 0, tot_light_amt, SHAPE_CYLIN);
		was_placed = 1;
		break; // done
	} // for n
	return was_placed;
}

void building_t::add_parking_garage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	unsigned num_floors, unsigned &nlights_x, unsigned &nlights_y, float &light_delta_z)
{
	assert(has_room_geom());
	rgen.rseed1 += 123*floor_ix; // make it unique per floor
	rgen.rseed2 += room_id;
	interior->room_geom->has_parking_garage = 1;
	// rows are separated by walls and run in dim, with a road and parking spaces on either side of it;
	// spaces are arranged in !dim, with roads along the edges of the building that connect to the roads of each row
	bool const dim(room.dx() < room.dy()); // long/primary dim; cars are lined up along this dim, oriented along the other dim
	vector3d const car_sz(get_nom_car_size()), parking_sz(1.1*car_sz.x, 1.4*car_sz.y, 1.5*car_sz.z); // space is somewhat larger than a car; car length:width = 2.3
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), wall_thickness(1.2*get_wall_thickness()), wall_hc(0.5*wall_thickness); // thicker
	float const ceiling_z(zval + window_vspacing - floor_thickness); // Note: zval is at floor level, not at the bottom of the room
	float const pillar_width(0.5*car_sz.y), pillar_hwidth(0.5*pillar_width), beam_hwidth(0.5*pillar_hwidth), road_width(2.3*car_sz.y); // road wide enough for two cars
	float const wid_sz(room.get_sz_dim(dim)), len_sz(room.get_sz_dim(!dim)), wid_sz_spaces(wid_sz - 2.0*road_width);
	float const min_strip_sz(2.0*parking_sz.x + road_width + max(wall_thickness, pillar_width)); // road + parking spaces on each side + wall/pillar
	assert(car_sz.z < (window_vspacing - floor_thickness)); // sanity check; may fail for some user parameters, but it's unclear what we do in that case
	unsigned const num_space_wid(wid_sz_spaces/parking_sz.y), num_full_strips(len_sz/min_strip_sz); // take the floor
	bool const half_strip((num_full_strips*min_strip_sz + parking_sz.x + road_width + wall_thickness) < len_sz); // no space for a full row, add a half row
	bool const half_row_side(half_strip ? rgen.rand_bool() : 0); // pick a random side
	unsigned const num_rows(2*num_full_strips + half_strip), num_strips(num_full_strips + half_strip), num_walls(num_strips - 1);
	unsigned const capacity(num_rows*num_space_wid); // ignoring space blocked by stairs and elevators
	unsigned &nlights_len(dim ? nlights_x : nlights_y), &nlights_wid(dim ? nlights_y : nlights_x);
	nlights_len = num_rows; // lights over each row of parking spaces
	nlights_wid = round_fp(0.25*wid_sz/parking_sz.y); // 4 parking spaces per light on average, including roads
	//cout << TXT(nlights_len) << TXT(nlights_wid) << TXT(num_space_wid) << TXT(num_rows) << TXT(capacity) << endl; // TESTING
	assert(num_space_wid   >= 4); // must fit at least 4 cars per row
	assert(num_full_strips >= 1);
	
	// add walls and pillars between strips
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());
	colorRGBA const wall_color(WHITE);
	cube_t room_floor_cube(room), virt_room_for_wall(room);
	set_cube_zvals(room_floor_cube, zval, ceiling_z);
	cube_t wall(room_floor_cube), pillar(room_floor_cube), beam(room_floor_cube);
	wall.expand_in_dim(dim, -road_width); // wall ends at roads that line the sides of the room; include pillar for better occluder and in case the pillar is skipped
	float wall_spacing(len_sz/(num_walls + 1));
	float const pillar_shift(0.01*pillar_width); // small value to avoid z-fighting
	float const wall_len(wall.get_sz_dim(dim) + 2.0f*pillar_shift), pillar_start(wall.d[dim][0] + pillar_hwidth - pillar_shift);
	float const row_width(wall_spacing - wall_thickness), space_length(0.5f*(row_width - road_width)), beam_spacing(len_sz/num_rows);
	unsigned const num_pillars(max(2U, unsigned(round_fp(0.25*wall_len/parking_sz.y)))); // every 4 spaces, at least 2 at the ends of the wall
	float const pillar_spacing((wall_len - pillar_width)/(num_pillars - 1)), beam_delta_z(0.95*wall.dz()), tot_light_amt(room.light_intensity);
	bool short_sides[2] = {0,0};

	if (half_strip) {
		short_sides[half_row_side] = 1;
		virt_room_for_wall.d[!dim][half_row_side] += (half_row_side ? 1.0 : -1.0)*space_length;
		wall_spacing = virt_room_for_wall.get_sz_dim(!dim)/(num_walls + 1); // recalculate wall spacing
	}
	light_delta_z = beam_delta_z - wall.dz(); // negative
	beam.z1()    += beam_delta_z; // shift the bottom up to the ceiling
	float const space_clearance(max(0.5f*window_vspacing, parking_sz.y)); // clearance between stairs/elevators and parking spaces so that cars and people can pass
	vect_cube_t obstacles, obstacles_exp, wall_parts, temp;
	// get obstacles for walls with and without clearance; maybe later add entrance/exit ramps, etc.
	interior->get_stairs_and_elevators_bcubes_intersecting_cube(room_floor_cube, obstacles, 0.0); // without clearance
	interior->get_stairs_and_elevators_bcubes_intersecting_cube(room_floor_cube, obstacles_exp, 0.9*window_vspacing); // with clearance
	cube_with_ix_t const &ramp(interior->pg_ramp);
	bool const is_top_floor(floor_ix+1 == num_floors);
	
	// add ramp if one was placed during floorplanning, before adding parking spaces
	// Note: lights can be very close to ramps, but I haven't actually seen them touch; do we need to check for and handle that case?
	if (!ramp.is_all_zeros()) {
		cube_t rc(ramp); // ramp clipped to this parking garage floor
		set_cube_zvals(rc, zval, (zval + window_vspacing));
		unsigned const flags((is_top_floor && interior->ignore_ramp_placement) ? 0 : RO_FLAG_OPEN); // ramp is open if the top exit is open
		objs.emplace_back(rc, TYPE_RAMP, room_id, (ramp.ix >> 1), (ramp.ix & 1), flags, tot_light_amt, SHAPE_ANGLED, wall_color);
		obstacles    .push_back(rc); // don't place parking spaces next to the ramp
		obstacles_exp.push_back(rc); // clip beams to ramp
	}
	// add walls and pillars
	if (interior->room_geom->pg_wall_start == 0) {interior->room_geom->pg_wall_start = objs.size();} // set if not set, on first level
	vect_cube_t pillars; // added after wall segments

	for (unsigned n = 0; n < num_walls+2; ++n) { // includes room far walls
		if (n < num_walls) { // interior wall
			float const pos(virt_room_for_wall.d[!dim][0] + (n + 1)*wall_spacing); // reference from the room far wall, assuming we can fit a full width double row strip
			set_wall_width(wall,   pos, wall_hc, !dim);
			set_wall_width(pillar, pos, pillar_hwidth, !dim);
			subtract_cubes_from_cube(wall, obstacles_exp, wall_parts, temp, 1); // ignore_zval=1
			
			for (auto const &w : wall_parts) {
				if (w.get_sz_dim(dim) < 2.0*window_vspacing) continue; // too short, skip
				objs.emplace_back(w, TYPE_PG_WALL, room_id, !dim, 0, 0, tot_light_amt, SHAPE_CUBE, wall_color, 0);
			}
		}
		else { // room wall
			bool const side(n == num_walls+1);
			pillar.d[!dim][ side] = room.d[!dim][side];
			pillar.d[!dim][!side] = room.d[!dim][side] + (side ? -1.0 : 1.0)*pillar_hwidth; // half the width of an interior wall pillar
		}
		for (unsigned p = 0; p < num_pillars; ++p) { // add support pillars
			float const ppos(pillar_start + p*pillar_spacing);
			set_wall_width(pillar, ppos, pillar_hwidth, dim);
			if (has_bcube_int_xy(pillar, obstacles_exp)) continue; // skip entire pillar if it intersects stairs or an elevator
			pillars.push_back(pillar);
		} // for p
	} // for n
	for (auto const &p : pillars) {objs.emplace_back(p, TYPE_PG_WALL, room_id, !dim, 0, 0, tot_light_amt, SHAPE_CUBE, wall_color, 1);}

	// add beams in !dim, at and between pillars
	unsigned const beam_flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING);

	for (unsigned p = 0; p < (4*(num_pillars - 1) + 1); ++p) { // add beams, 4 per pillar
		float const ppos(pillar_start + 0.25*p*pillar_spacing);
		set_wall_width(beam, ppos, beam_hwidth, dim);
		subtract_cubes_from_cube(beam, obstacles, wall_parts, temp, 1); // ignore_zval=1
		
		for (auto const &w : wall_parts) {
			if (min(w.dx(), w.dy()) > beam_hwidth) {objs.emplace_back(w, TYPE_PG_WALL, room_id, !dim, 0, beam_flags, tot_light_amt, SHAPE_CUBE, wall_color, 2);}
		}
	}
	// add beams in dim for each row of lights
	for (unsigned n = 0; n < num_rows; ++n) {
		float const pos(room.d[!dim][0] + (n + 0.5)*beam_spacing);
		cube_t beam(room_floor_cube);
		beam.z1() += beam_delta_z; // shift the bottom up to the ceiling
		set_wall_width(beam, pos, beam_hwidth, !dim);
		subtract_cubes_from_cube(beam, obstacles, wall_parts, temp, 1); // ignore_zval=1

		for (auto const &w : wall_parts) {
			if (min(w.dx(), w.dy()) > beam_hwidth) {objs.emplace_back(w, TYPE_PG_WALL, room_id, !dim, 0, beam_flags, tot_light_amt, SHAPE_CUBE, wall_color, 2);}
		}
	}

	// add parking spaces on both sides of each row (one side if half row)
	cube_t row(wall); // same length as the wall; includes the width of the pillars
	row.z2() = row.z1() + 0.001*window_vspacing; // slightly above the floor
	float const space_width(row.get_sz_dim(dim)/num_space_wid), strips_start(virt_room_for_wall.d[!dim][0]);
	bool const add_cars(city_params.num_cars > 0 && !city_params.car_model_files.empty() && !is_rotated()); // skip cars for rotated buildings

	for (unsigned n = 0; n < num_strips; ++n) {
		row.d[!dim][0] = strips_start + (n + 0)*wall_spacing + wall_hc;
		row.d[!dim][1] = strips_start + (n + 1)*wall_spacing - wall_hc;
		assert(space_length > 0.0);

		for (unsigned d = 0; d < 2; ++d) { // for each side of the row
			bool const at_ext_wall[2] = {(n == 0 && d == 0), (n+1 == num_strips && d == 1)};
			if ((short_sides[0] && at_ext_wall[0]) || (short_sides[1] && at_ext_wall[1])) continue; // skip this row
			float row_left_edge(row.d[dim][0]); // spaces start flush with the row, or flush with the room if this is the exterior wall
			unsigned num_spaces_per_row(num_space_wid);

			if (at_ext_wall[0] || at_ext_wall[1]) { // at either room exterior wall - can extend spaces up to the wall
				float row_right_edge(row.d[dim][1]); // opposite end of the row
				while ((row_left_edge  - space_width) > room.d[dim][0]) {row_left_edge  -= space_width; ++num_spaces_per_row;} // add rows to the left
				while ((row_right_edge + space_width) < room.d[dim][1]) {row_right_edge += space_width; ++num_spaces_per_row;} // add rows to the right
			}
			cube_t space(row);
			space.d[!dim][!d] += (d ? 1.0 : -1.0)*(row_width - space_length); // shrink
			space.d[ dim][0]   = row_left_edge;
			bool last_was_space(0);
			
			for (unsigned s = 0; s < num_spaces_per_row; ++s) {
				space.d[dim][1] = space.d[dim][0] + space_width; // set width
				assert(space.is_strictly_normalized());
				
				if (has_bcube_int_xy(space, obstacles, space_clearance)) { // skip entire space if it intersects stairs or an elevator, with padding
					if (last_was_space) {objs.back().flags &= ~RO_FLAG_ADJ_HI;} // no space to the right for the previous space
					last_was_space = 0;
				}
				else {
					unsigned flags(RO_FLAG_NOCOLL);
					if (last_was_space          ) {flags |= RO_FLAG_ADJ_LO;} // adjacent space to the left
					if (s+1 < num_spaces_per_row) {flags |= RO_FLAG_ADJ_HI;} // not the last space - assume there will be a space to the right
					bool const add_car(add_cars && rgen.rand_float() < 0.5); // 50% populated with cars
					room_object_t pspace(space, TYPE_PARK_SPACE, room_id, !dim, d, flags, tot_light_amt, SHAPE_CUBE, wall_color); // floor_color?

					if (add_car) { // add a collider to block this area from the player, people, and rats; add first so that objs.back() is correct for the next iter
						car_t car(car_from_parking_space(pspace));
						objs.emplace_back(car.bcube, TYPE_COLLIDER, room_id, !dim, d, RO_FLAG_INVIS);
						pspace.obj_id = (uint16_t)(objs.size() + rgen.rand()); // will be used for the car model and color
						pspace.flags |= RO_FLAG_USED;
					}
					objs.push_back(pspace);
					last_was_space = 1;
				}
				space.d[dim][0] = space.d[dim][1]; // shift to next space
			} // for s
		} // for d
	} // for n
	if (is_top_floor) {
		//highres_timer_t timer("Get Pipe Basement Connections");
		// move or remove pipes intersecting lights, pillars, walls, stairs, elevators, and ramps;
		// note that lights haven't been added yet though, so maybe pipes need to be added later?
		vect_cube_t walls;

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if      (i->type == TYPE_PG_WALL) {walls    .push_back(*i);} // walls and pillars
			else if (i->type == TYPE_RAMP   ) {obstacles.push_back(*i);} // ramps are obstacles for pipes
		}
		add_basement_pipes(obstacles, walls, room_id, tot_light_amt);
	}
}

// find the closest wall (including room wall) to this location, avoiding obstacles, and shift outward by radius; routes in X or Y only, for now
point get_closest_wall_pos(point const &pos, float radius, cube_t const &room, vect_cube_t const &walls, vect_cube_t const &obstacles) {
	if (!room.contains_pt_xy_exp(pos, radius)) {return pos;} // error?
	// what about checking pos intersecting walls or obstacles? is that up to the caller to handle?
	vector3d const expand(radius, radius, radius);
	point best(pos);
	float dmin(room.dx() + room.dy()); // use an initial distance larger than what we can return

	if (!room.is_all_zeros()) { // check room exterior walls first
		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				float const val(room.d[dim][dir] + (dir ? -radius : radius)), dist(fabs(val - pos[dim])); // shift val inward
				if (dist >= dmin) continue;
				point cand(pos);
				cand[dim] = val;
				// check walls as well, even though any wall hit should be replaced with a closer point below
				if (!line_int_cubes_exp(pos, cand, obstacles, expand) && !line_int_cubes_exp(pos, cand, walls, expand)) {best = cand; dmin = dist;}
			} // for dir
		} // for dim
	}
	for (cube_t const &wall : walls) { // check all interior walls
		for (unsigned dim = 0; dim < 2; ++dim) {
			bool const dir(wall.get_center_dim(dim) < pos[dim]);
			float const val(wall.d[dim][dir] - (dir ? -radius : radius)), dist(fabs(val - pos[dim])); // shift val outward
			if (dist >= dmin) continue;
			point cand(pos);
			cand[dim] = val;
			if (!line_int_cubes_exp(pos, cand, obstacles, expand)) {best = cand; dmin = dist;} // check obstacles only
		} // for dim
	}
	return best;
}

float get_merged_pipe_radius(float r1, float r2) {return pow((r1*r1*r1 + r2*r2*r2), 1/3.0);} // scales as cubic

struct pipe_t {
	point p1, p2;
	float radius;
	unsigned dim;
	bool connected, is_exit;

	pipe_t(point const &p1_, point const &p2_, float radius_, unsigned dim_, bool connected_=0) :
		p1(p1_), p2(p2_), radius(radius_), dim(dim_), connected(connected_), is_exit(0) {}

	cube_t get_bcube() const {
		cube_t bcube(p1, p2);
		bcube.expand_by(radius);
		bcube.expand_in_dim(dim, -radius); // oops, we shouldn't have expanded in this dim
		return bcube;
	}
};

void building_t::add_basement_pipes(vect_cube_t const &obstacles, vect_cube_t const &walls, unsigned room_id, float tot_light_amt) {
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t const &basement(get_basement());

	// get pipe ends coming in through the ceiling
	vector<sphere_t> pipe_ends;
	get_pipe_basement_connections(pipe_ends);
	if (pipe_ends.empty()) return; // can this happen?
	float rmax(0.0);
	for (sphere_t const &p : pipe_ends) {max_eq(rmax, p.radius);} // calculate max radius across all pipes
	float const ceil_zval(pipe_ends.front().pos.z); // should be the same for all pipe ends
	float const pipe_zval(ceil_zval - 1.5f*max(rmax, get_fc_thickness()));
	vector<pipe_t> pipes;
	map<float, vector<unsigned>> xy_map[2];
	unsigned num_connected(0);

	// seed the pipe graph with valid vertical segments and build a graph of X/Y values
	for (sphere_t const &p : pipe_ends) {
		assert(p.radius > 0.0);
		assert(p.pos.z > bcube.z1());
		cube_t c(p.pos);
		c.expand_by_xy(p.radius);
		c.z1() = bcube.z1(); // extend all the way down to the floor of the lowest basement

		if (has_bcube_int(c, obstacles) || has_bcube_int(c, walls)) { // can't place over stairs/elevators/ramps/pillars/walls
			continue; // skip for now; TODO: maybe should move the pipe so as not to collide, or merge with a nearby pipe?
		}
		for (unsigned d = 0; d < 2; ++d) {
			// TODO: align if very close together
			xy_map[d][p.pos[d]].push_back(pipes.size());
		}
		pipes.emplace_back(point(p.pos.x, p.pos.y, pipe_zval), p.pos, p.radius, 2);
		++num_connected;
	} // for pipe_ends

	// connect pipes in the same row or column
	for (unsigned d = 0; d < 2; ++d) {
		for (auto const &v : xy_map[!d]) {
			pipe_t exit_pipe(pipes[v.second.front()]); // deep copy
			if (v.second.size() == 1) continue; // skip singleton
			float radius(0.0), range_min(basement.d[d][1]), range_max(basement.d[d][0]); // start denormalized

			for (unsigned ix : v.second) {
				assert(ix < pipes.size());
				pipe_t &pipe(pipes[ix]);
				float const val(pipe.p1[d]);
				min_eq(range_min, val-pipe.radius);
				max_eq(range_max, val+pipe.radius);
				radius = get_merged_pipe_radius(radius, + pipe.radius);
				pipe.connected = 1;
			} // for ix
			assert(range_min < range_max);
			float const val(exit_pipe.p1[d]);
			min_eq(range_min, val-radius); // extend range to include the exit pipe
			max_eq(range_max, val+radius);
			point p1(exit_pipe.p1), p2(p1);
			p1[d] = range_min; p2[d] = range_max;
			pipe_t const pipe(p1, p2, radius, d, 1); // connected=1
			cube_t const c(pipe.get_bcube());
			// TODO: check c against obstacles and split/clip/drop
			// TODO: track connected and avoid loops
			pipes.push_back(pipe);
			// add exit pipe; TODO: only one per connected graph; must call get_closest_wall_pos()
			exit_pipe.is_exit = 1;
			exit_pipe.radius  = radius;
			exit_pipe.p2.z    = exit_pipe.p1.z;
			exit_pipe.p1.z    = basement.z1(); // exits to the basement floor
			pipes.push_back(exit_pipe);
		} // for v
	} // for d

	// handle any remaining unconnected pipes
	for (pipe_t &p : pipes) {
		if (p.connected) continue; // already connected, has an exit
		p.p1.z    = basement.z1(); // exits to the basement floor
		p.is_exit = 1;
	}

	// add pipe objects
	for (pipe_t const &p : pipes) {
		bool const pdim(p.dim & 1), pdir(p.dim >> 1); // encoded as: X:dim=0,dir=0 Y:dim=1,dir=0, Z:dim=x,dir=1
		unsigned flags(0);
		if      (p.is_exit ) {flags = RO_FLAG_ADJ_HI;} // collidable, round end at top
		else if (p.dim == 2) {flags = RO_FLAG_NOCOLL;}
		else                 {flags = (RO_FLAG_NOCOLL | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_HI | RO_FLAG_HANGING);} // add both flat ends
		cube_t const c(p.get_bcube());
		objs.emplace_back(c, TYPE_PIPE, room_id, pdim, pdir, flags, tot_light_amt, SHAPE_CYLIN, DK_GRAY);
	} // for p
	cout << TXT(obstacles.size()) << TXT(pipe_ends.size()) << TXT(pipes.size()) << TXT(num_connected) << endl;
}

// here each sphere represents the entry point of a pipe with this radius into the basement ceiling
// find all plumbing fixtures such as toilets, urinals, sinks, and showers; these should have all been placed in rooms by now
void building_t::get_pipe_basement_connections(vector<sphere_t> &pipes) const {
	float const merge_dist = 4.0; // merge two pipes if their combined radius is within this distance
	float const floor_spacing(get_window_vspace()), base_pipe_radius(0.01*floor_spacing), base_pipe_area(base_pipe_radius*base_pipe_radius);
	float const pipe_z2(ground_floor_z1 - get_fc_thickness()), merge_dist_sq(merge_dist*merge_dist), max_radius(0.4*get_wall_thickness());
	auto const &objs(interior->room_geom->objs);
	unsigned num_drains(0);
	cube_t const &basement(get_basement());

	for (auto i = objs.begin(); i != objs.end(); ++i) { // check all objects placed so far
		if (i->type != TYPE_TOILET && i->type != TYPE_SINK && i->type != TYPE_URINAL && i->type != TYPE_TUB && i->type != TYPE_SHOWER &&
			i->type != TYPE_BRSINK && i->type != TYPE_KSINK && i->type != TYPE_WASHER && i->type != TYPE_DRAIN) continue;
		point const pos(i->xc(), i->yc(), pipe_z2);
		if (!basement.contains_pt_xy(pos)) continue; // pipe doesn't pass through the basement, skip
		bool merged(0);

		// see if we can merge this pipe into an existing nearby pipe
		for (auto p = pipes.begin(); p != pipes.end(); ++p) {
			float const p_area(p->radius*p->radius), sum_area(p_area + base_pipe_area);
			if (!dist_xy_less_than(p->pos, pos, merge_dist_sq*sum_area)) continue;
			p->pos    = (p_area*p->pos + base_pipe_area*pos)/sum_area; // merged position is weighted average area
			p->radius = get_merged_pipe_radius(p->radius, base_pipe_radius);
			merged    = 1;
			break;
		} // for p
		if (!merged) {pipes.emplace_back(pos, base_pipe_radius);} // create a new pipe
		++num_drains;
	} // for i
	for (sphere_t &p : pipes) {min_eq(p.radius, max_radius);} // clamp radius to a reasonable value after all merges
	cout << TXT(objs.size()) << TXT(num_drains) << TXT(pipes.size()) << endl;
}

void building_t::add_parking_garage_ramp(rand_gen_t &rgen) {
	assert(interior && !is_house && has_parking_garage);
	cube_with_ix_t &ramp(interior->pg_ramp);
	assert(ramp.is_all_zeros()); // must not have been set
	cube_t const &basement(get_basement());
	bool const dim(basement.dx() < basement.dy()); // long/primary dim
	// see building_t::add_parking_garage_objs(); make sure there's space for a ramp plus both exit dirs within the building width
	float const width(basement.get_sz_dim(!dim)), road_width(min(0.25f*width, 2.3f*get_nom_car_size().y));
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const z1(basement.z1() + fc_thick), z2(basement.z2() + fc_thick); // bottom level basement floor to first floor floor
	bool const ramp_pref_xdir(rgen.rand_bool()), ramp_pref_ydir(rgen.rand_bool());
	bool added_ramp(0), dir(0);

	for (unsigned pass = 0; pass < 2 && !added_ramp; ++pass) {
		for (unsigned xd = 0; xd < 2 && !added_ramp; ++xd) {
			for (unsigned yd = 0; yd < 2; ++yd) {
				bool const xdir(bool(xd) ^ ramp_pref_xdir), ydir(bool(yd) ^ ramp_pref_ydir);
				float const xsz((dim ? 2.0 : 1.0)*road_width), ysz((dim ? 1.0 : 2.0)*road_width); // longer in !dim
				unsigned const num_ext(unsigned(basement.d[0][xdir] == bcube.d[0][xdir]) + unsigned(basement.d[1][ydir] == bcube.d[1][ydir]));
				if (num_ext < 2-pass) continue; // must be on the exterior edge of the building in both dims for pass 0, and one dim for pass 1
				dir = (dim ? xdir : ydir);
				point corner(basement.d[0][xdir], basement.d[1][ydir], z1);
				corner[!dim] += (dir ? -1.0 : 1.0)*road_width; // shift away from the wall so that cars have space to turn onto the level floor
				point const c1((corner.x - 0.001*(xdir ? 1.0 : -1.0)*xsz), (corner.y - 0.001*(ydir ? 1.0 : -1.0)*ysz), z1); // slight inward shift to prevent z-fighting
				point const c2((corner.x + (xdir ? -1.0 : 1.0)*xsz), (corner.y + (ydir ? -1.0 : 1.0)*ysz), z2);
				ramp = cube_with_ix_t(cube_t(c1, c2), (((!dim)<<1) + dir)); // encode dim and dir in ramp index field
				added_ramp = 1;
				break; // done
			} // for yd
		} // for xd
	} // for pass
	if (!added_ramp) return; // what if none of the 4 corners work for a ramp?
	// add landings, which are used to draw the vertical edges of the cutout
	unsigned num_floors(calc_num_floors(basement, window_vspacing, floor_thickness));
	float z(basement.z1() + window_vspacing); // start at upper floor rather than lower floor

	if (1) { // FIXME: rooms on the ground floor above ramps aren't yet handled, so clip ramps to avoid disrupting their floors until this is fixed
		ramp.z2() -= 2.0*floor_thickness;
		--num_floors;
		interior->ignore_ramp_placement = 1; // okay to place room objects over ramps because the floor has not been removed
	}
	for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) { // skip first floor - draw pairs of floors and ceilings
		landing_t landing(ramp, 0, f, !dim, dir, 0, SHAPE_RAMP, 0, (f+1 == num_floors), 0, 1); // for_ramp=1
		set_cube_zvals(landing, (z - fc_thick), (z + fc_thick));
		interior->landings.push_back(landing);
	}
	// cut out spaces from floors and ceilings
	subtract_cube_from_floor_ceil(ramp, interior->floors  );
	subtract_cube_from_floor_ceil(ramp, interior->ceilings);
}

