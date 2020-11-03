// 3D World - Building Interior Room Geometry Placement
// by Frank Gennari 4/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#include "profiler.h"
#pragma warning(disable : 26812) // prefer enum class over enum

extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;


bool building_t::overlaps_other_room_obj(cube_t const &c, unsigned objs_start) const {
	assert(has_room_geom());
	vector<room_object_t> &objs(interior->room_geom->objs);
	assert(objs_start <= objs.size());

	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (!(i->flags & RO_FLAG_NOCOLL) && i->intersects(c)) return 1;
	}
	return 0;
}

bool building_t::is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, bool inc_open_doors, float room_pad) const {
	cube_t place_area(room);
	if (room_pad != 0.0f) {place_area.expand_by_xy(-room_pad);} // shrink by dmin
	if (!place_area.contains_cube_xy(c)) return 0; // not contained in interior part of the room
	if (is_cube_close_to_doorway(c, room, 0.0, inc_open_doors)) return 0; // too close to a doorway
	if (interior && interior->is_blocked_by_stairs_or_elevator(c)) return 0; // faster to check only one per stairwell, but then we need to store another vector?
	if (has_bcube_int(c, blockers)) return 0; // Note: ignores dmin
	return 1;
}

float get_radius_for_square_model(unsigned model_id) {
	vector3d const chair_sz(building_obj_model_loader.get_model_world_space_size(model_id));
	return 0.5f*(chair_sz.x + chair_sz.y)/chair_sz.z; // assume square and take average of xsize and ysize
}
template<typename T> cube_t get_cube_height_radius(point const &center, T radius, float height) { // T can be float or vector3d
	cube_t c(center);
	c.expand_by_xy(radius);
	c.z2() += height;
	return c;
}

bool building_t::add_chair(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id, point const &place_pos,
	colorRGBA const &chair_color, bool dim, bool dir, float tot_light_amt, bool office_chair_model)
{
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR)) {office_chair_model = 0;}
	float const window_vspacing(get_window_vspace()), room_pad(4.0f*get_wall_thickness()), chair_height(0.4*window_vspacing);
	float chair_hwidth(0.0), push_out(0.0);
	point chair_pos(place_pos); // same starting center and z1

	if (office_chair_model) {
		chair_hwidth = 0.5f*chair_height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR);
		push_out     = 0.5; // pushed out a bit so that the arms don't intersect the table top
	}
	else {
		chair_hwidth = 0.1*window_vspacing; // half width
		push_out     = rgen.rand_uniform(-0.5, 1.2); // varible amount of pushed in/out
	}
	chair_pos[dim] += (dir ? -1.0f : 1.0f)*push_out*chair_hwidth;
	cube_t const chair(get_cube_height_radius(chair_pos, chair_hwidth, chair_height));
	if (!is_valid_placement_for_room(chair, room, blockers, 0, room_pad)) return 0; // check proximity to doors
	vector<room_object_t> &objs(interior->room_geom->objs);

	if (office_chair_model) {
		float const lum(0.4*chair_color.R + 0.2*chair_color.B + 0.1*chair_color.G); // calculate grayscale luminance
		objs.emplace_back(chair, TYPE_OFFICE_CHAIR, room_id, dim, dir, 0, tot_light_amt, room_obj_shape::SHAPE_CUBE, colorRGBA(lum, lum, lum));
	}
	else {
		objs.emplace_back(chair, TYPE_CHAIR, room_id, dim, dir, 0, tot_light_amt, room_obj_shape::SHAPE_CUBE, chair_color);
	}
	return 1;
}

// Note: must be first placed objects; returns the number of total objects added (table + optional chairs)
unsigned building_t::add_table_and_chairs(rand_gen_t rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id,
	point const &place_pos, colorRGBA const &chair_color, float rand_place_off, float tot_light_amt)
{
	float const window_vspacing(get_window_vspace()), room_pad(max(4.0f*get_wall_thickness(), get_min_front_clearance()));
	vector3d const room_sz(room.get_size());
	vector<room_object_t> &objs(interior->room_geom->objs);
	point table_pos(place_pos);
	vector3d table_sz;
	for (unsigned d = 0; d < 2; ++d) {table_sz [d]  = 0.18*window_vspacing*(1.0 + rgen.rand_float());} // half size relative to window_vspacing
	for (unsigned d = 0; d < 2; ++d) {table_pos[d] += rand_place_off*room_sz[d]*rgen.rand_uniform(-1.0, 1.0);} // near the center of the room
	bool const is_round((rgen.rand()&3) == 0); // 25% of the time
	if (is_round) {table_sz.x = table_sz.y = 0.6f*(table_sz.x + table_sz.y);} // round tables must have square bcubes for now (no oval tables yet); make radius slightly larger
	point llc(table_pos - table_sz), urc(table_pos + table_sz);
	llc.z = table_pos.z; // bottom
	urc.z = table_pos.z + rgen.rand_uniform(0.20, 0.22)*window_vspacing; // top
	cube_t table(llc, urc);
	if (!is_valid_placement_for_room(table, room, blockers, 0, room_pad)) return 0; // check proximity to doors and collision with blockers
	objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, 0, tot_light_amt, (is_round ? SHAPE_CYLIN : SHAPE_CUBE));
	unsigned num_added(1); // start with the table

	// place some chairs around the table
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (rgen.rand_bool()) continue; // 50% of the time
			point chair_pos(table_pos); // same starting center and z1
			chair_pos[dim] += (dir ? -1.0f : 1.0f)*table_sz[dim];
			num_added += add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, 0); // office_chair_model=0
		}
	}
	return num_added;
}
void building_t::shorten_chairs_in_region(cube_t const &region, unsigned objs_start) {
	for (auto i = interior->room_geom->objs.begin() + objs_start; i != interior->room_geom->objs.end(); ++i) {
		if (i->type != TYPE_CHAIR || !i->intersects(region)) continue;
		i->z2() -= 0.25*i->dz();
		i->shape = SHAPE_SHORT;
	}
}

void building_t::add_trashcan_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool check_last_obj) {
	unsigned const NUM_COLORS = 6;
	colorRGBA const colors[NUM_COLORS] = {BLUE, DK_GRAY, LT_GRAY, GRAY, BLUE, WHITE};
	int const rr(rgen.rand()%3), rar(rgen.rand()%3); // three sizes/ARs
	float const radius(0.02f*(3 + rr)*get_window_vspace()), height(0.55f*(3 + rar)*radius); // radius={0.06, 0.08, 0.10} x AR={1.65, 2.2, 2.75}
	cube_t room_bounds(get_walkable_room_bounds(room)), room_exp(room);
	room_bounds.expand_by_xy(-1.1*radius); // leave a slight gap between trashcan and wall
	if (!room_bounds.is_strictly_normalized()) return; // no space for trashcan (likely can't happen)
	int const floor_ix(int((zval - room.z1())/get_window_vspace()));
	bool const cylin(((mat_ix + 13*real_num_parts + 5*hallway_dim + 131*floor_ix) % 7) < 4); // varies per-building, per-floor
	point center;
	center.z = zval + 0.001*get_window_vspace(); // slightly above the floor to avoid z-fighting
	unsigned skip_wall(4); // start at an invalid value
	// find interior doorways connected to this room
	float const wall_thickness(get_wall_thickness());
	room_exp.expand_by(wall_thickness, wall_thickness, -wall_thickness); // expand in XY and shrink in Z
	vector<room_object_t> &objs(interior->room_geom->objs);
	cube_t avoid;
	vect_cube_t doorways;

	if (!objs.empty() && objs[objs_start].type == TYPE_TABLE) { // make sure there's enough space for the player to walk around the table
		avoid = objs[objs_start];
		avoid.expand_by_xy(get_min_front_clearance());
	}
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
		if (i->intersects(room_exp)) {doorways.push_back(*i);}
	}
	if (check_last_obj) {
		assert(!objs.empty());
		skip_wall = 2*objs.back().dim + (!objs.back().dir); // don't place trashcan on same wall as whiteboard (dir is opposite)
	}
	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a trashcan
		bool dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		if ((2U*dim + dir) == skip_wall) {dir ^= 1;} // don't place a trashcan on this wall, try opposite wall
		center[dim] = room_bounds.d[dim][dir]; // against this wall
		bool is_good(0);

		for (unsigned m = 0; m < 80; ++m) { // try to find a point near a doorway
			center[!dim] = rgen.rand_uniform(room_bounds.d[!dim][0], room_bounds.d[!dim][1]);
			if (doorways.empty()) break; // no doorways, keep this point
				
			for (auto i = doorways.begin(); i != doorways.end(); ++i) {
				float const dmin(radius + i->dx() + i->dy());
				if (!i->closest_dist_less_than(center, 2.0*dmin)) continue; // too far
				if ( i->closest_dist_less_than(center, dmin)) {is_good = 0; break;} // too close, reject this point
				is_good = 1; // close enough, keep this point
			}
			if (is_good) break; // done; may never get here if no points are good, but the code below will handle that
		} // for m
		cube_t const c(get_cube_height_radius(center, radius, height));
		if (!avoid.is_all_zeros() && c.intersects_xy(avoid)) continue; // bad placement
		if (is_cube_close_to_doorway(c, room, 0.0, !room.is_hallway) || interior->is_blocked_by_stairs_or_elevator(c) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_TCAN, room_id, dim, dir, 0, tot_light_amt, (cylin ? room_obj_shape::SHAPE_CYLIN : room_obj_shape::SHAPE_CUBE), colors[rgen.rand()%NUM_COLORS]);
		return; // done
	} // for n
}

// Note: no blockers for people
bool building_t::add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	cube_t room_bounds(get_walkable_room_bounds(room));
	room_bounds.expand_by_xy(-0.1*get_wall_thickness());
	float const vspace(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return 0; // room is too small
	rand_gen_t rgen2;
	rgen2.set_state((room_id + 1), (13*mat_ix + interior->rooms.size() + 1)); // local rgen that's per-building/room; ensures bookcases are all the same size in a library
	float const width(0.4*vspace*rgen2.rand_uniform(1.0, 1.2)), depth(0.12*vspace*rgen2.rand_uniform(1.0, 1.2)), height(0.7*vspace*rgen2.rand_uniform(1.0, 1.2));
	float const clearance(max(0.2f*vspace, get_min_front_clearance()));
	vector<room_object_t> &objs(interior->room_geom->objs);
	cube_t c;
	c.z1() = zval;
	c.z2() = zval + height;

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a bookcase
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		if (classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) continue; // don't place against an exterior wall/window, inc. partial ext walls
		c.d[dim][ dir] = room_bounds.d[dim][dir]; // against this wall
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		c.d[!dim][0] = pos - 0.5*width;
		c.d[!dim][1] = pos + 0.5*width;
		cube_t tc(c);
		tc.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance; // increase space to add clearance
		if (is_cube_close_to_doorway(tc, room, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(tc) || overlaps_other_room_obj(tc, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_BCASE, room_id, dim, !dir, 0, tot_light_amt); // Note: dir faces into the room, not the wall
		objs.back().obj_id = (uint16_t)objs.size();
		return 1; // done/success
	} // for n
	return 0; // not placed
}

bool building_t::room_has_stairs_or_elevator(room_t const &room, float zval) const {
	if (room.has_elevator) return 1; // elevator shafts extend through all rooms in a stack, don't need to check zval
	if (!room.has_stairs ) return 0; // no stairs
	assert(interior);
	cube_t c(room);
	c.z1() = zval;
	c.z2() = zval + 0.9*get_window_vspace();

	for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
		if (s->intersects(c)) return 1;
	}
	return 0;
}

// Note: must be first placed object
bool building_t::add_desk_to_room(rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers,
	colorRGBA const &chair_color, float zval, unsigned room_id, float tot_light_amt)
{
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return 0; // room is too small
	float const width(0.8*vspace*rgen.rand_uniform(1.0, 1.2)), depth(0.38*vspace*rgen.rand_uniform(1.0, 1.2)), height(0.21*vspace*rgen.rand_uniform(1.0, 1.2));
	vector<room_object_t> &objs(interior->room_geom->objs);
	unsigned num_placed(0);
	cube_t c, placed_desk;
	c.z1() = zval;
	c.z2() = zval + height;

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a desk
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		float const dsign(dir ? -1.0 : 1.0);
		c.d[dim][ dir] = room_bounds.d[dim][dir] + rgen.rand_uniform(0.1, 1.0)*dsign*get_wall_thickness(); // almost against this wall
		c.d[dim][!dir] = c.d[dim][dir] + dsign*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		c.d[!dim][0] = pos - 0.5*width;
		c.d[!dim][1] = pos + 0.5*width;
		if (num_placed > 0 && c.intersects(placed_desk)) continue; // intersects previously placed desk
		if (!is_valid_placement_for_room(c, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
		bool const is_tall(!room.is_office && rgen.rand_float() < 0.5 && classify_room_wall(room, zval, dim, dir, 0) != ROOM_WALL_EXT); // make short if against an exterior wall or in an office
		objs.emplace_back(c, TYPE_DESK, room_id, dim, !dir, 0, tot_light_amt, (is_tall ? SHAPE_TALL : SHAPE_CUBE));
		objs.back().obj_id = (uint16_t)objs.size();
		cube_t bc(c);
		bool const add_monitor(building_obj_model_loader.is_model_valid(OBJ_MODEL_TV) && rgen.rand_bool());

		if (add_monitor) { // add a computer monitor using the TV model
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
			float const tv_height(1.1*height), tv_hwidth(0.5*tv_height*sz.y/sz.z), tv_depth(tv_height*sz.x/sz.z), center(c.get_center_dim(!dim));
			cube_t tv;
			tv.z1() = c.z2();
			tv.z2() = c.z2() + tv_height;
			tv.d[dim][ dir] = c. d[dim][dir] + dsign*0.25*depth; // 25% of the way from the wall
			tv.d[dim][!dir] = tv.d[dim][dir] + dsign*tv_depth;
			set_wall_width(tv, center, tv_hwidth, !dim);
			objs.emplace_back(tv, TYPE_TV, room_id, dim, !dir, 0, tot_light_amt, room_obj_shape::SHAPE_SHORT, BLACK); // monitors are shorter than TVs
			objs.back().obj_id = (uint16_t)objs.size();
			// add a keyboard as well
			float const kbd_hwidth(0.7*tv_hwidth);
			cube_t keyboard;
			keyboard.z1() = c.z2();
			keyboard.z2() = c.z2() + 0.06*kbd_hwidth;
			keyboard.d[dim][!dir] = c.d[dim][!dir] - dsign*0.06*depth; // close to front edge
			keyboard.d[dim][ dir] = keyboard.d[dim][!dir] - dsign*0.6*kbd_hwidth;
			set_wall_width(keyboard, center, kbd_hwidth, !dim);
			objs.emplace_back(keyboard, TYPE_KEYBOARD, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // add as white, will be drawn with gray/black texture
		}
		if (rgen.rand_float() > 0.05) { // 5% chance of no chair
			point chair_pos;
			chair_pos.z = zval;
			chair_pos[dim]  = c.d[dim][!dir];
			chair_pos[!dim] = pos + rgen.rand_uniform(-0.1, 0.1)*width; // slightly misaligned
			// there are too many desks in office buildings, and they have office chairs in cubicles anyway, so only use chair models for desks in houses with computer monitors
			bool const office_chair_model(add_monitor && is_house);
			
			if (add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, office_chair_model)) {
				cube_t const &chair(objs.back());
				if (num_placed > 0 && chair.intersects(placed_desk)) {objs.pop_back();} // intersects previously placed desk, remove it
				else {bc.union_with_cube(chair);} // include the chair
			}
		}
		++num_placed;
		if (room.is_office && num_placed == 1 && rgen.rand_float() < 0.5 && !room_has_stairs_or_elevator(room, zval)) {placed_desk = bc; continue;} // allow two desks in one office
		break; // done/success
	} // for n
	return (num_placed > 0);
}

bool building_t::create_office_cubicles(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) { // assumes no prior placed objects
	if (!room.is_office) return 0; // offices only
	if (!room.interior && (rgen.rand()%3) == 0) return 0; // 66.7% chance for non-interior rooms
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const floor_spacing(get_window_vspace());
	// Note: we could choose the primary dim based on door placement like in office building bathrooms, but it seems easier to not place cubes by doors
	bool const long_dim(room.dx() < room.dy());
	float const rlength(room_bounds.get_sz_dim(long_dim)), rwidth(room_bounds.get_sz_dim(!long_dim)), midpoint(room_bounds.get_center_dim(!long_dim));
	if (rwidth < 2.5*floor_spacing || rlength < 3.5*floor_spacing) return 0; // not large enough
	unsigned const num_cubes(round_fp(rlength/(rgen.rand_uniform(0.75, 0.9)*floor_spacing))); // >= 4
	float const cube_width(rlength/num_cubes), cube_depth(cube_width*rgen.rand_uniform(0.8, 1.2)); // not quite square
	bool const add_middle_col(rwidth > 4.0*cube_depth + 2.0*interior->get_doorway_width()); // enough to fit 4 rows of cubes and 2 hallways in between
	uint16_t const bldg_id(uint16_t(mat_ix + interior->rooms.size())); // some value that's per-building
	cube_t const &part(get_part_for_room(room));
	vector<room_object_t> &objs(interior->room_geom->objs);
	bool const has_office_chair(building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR));
	float lo_pos(room_bounds.d[long_dim][0]), chair_height(0.0), chair_radius(0.0);
	cube_t c;
	c.z1() = zval;
	c.z2() = zval + 0.425*floor_spacing; // set height
	bool added_cube(0);

	if (has_office_chair) {
		chair_height = 0.425*floor_spacing;
		chair_radius = 0.5f*chair_height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR);
	}
	for (unsigned n = 0; n < num_cubes; ++n) {
		float const hi_pos(lo_pos + cube_width);
		c.d[long_dim][0] = lo_pos;
		c.d[long_dim][1] = hi_pos;

		for (unsigned is_middle = 0; is_middle < (add_middle_col ? 2U : 1U); ++is_middle) {
			if (is_middle && (n == 0 || n+1 == num_cubes)) continue; // skip end rows for middle section

			for (unsigned dir = 0; dir < 2; ++dir) {
				float const wall_pos(is_middle ? midpoint : room_bounds.d[!long_dim][dir]), dir_sign(dir ? -1.0 : 1.0);
				c.d[!long_dim][ dir] = wall_pos;
				c.d[!long_dim][!dir] = wall_pos + dir_sign*cube_depth;
				cube_t test_cube(c);
				test_cube.d[!long_dim][!dir] += dir_sign*0.5*cube_depth; // allow space for people to enter the cubicle
				if (interior->is_cube_close_to_doorway(test_cube, room, 0.0, 1, 1)) continue; // too close to a doorway; inc_open=1, check_zval=1
				if (interior->is_blocked_by_stairs_or_elevator(test_cube)) continue;
				bool const against_window(room.d[!long_dim][dir] == part.d[!long_dim][dir]);
				objs.emplace_back(c, TYPE_CUBICLE, room_id, !long_dim, dir, 0, tot_light_amt, ((against_window && !is_middle) ? SHAPE_SHORT : SHAPE_CUBE));
				objs.back().obj_id = bldg_id;
				added_cube = 1;
				// add colliders to allow the player to enter the cubicle but not cross the side walls
				cube_t c2(c), c3(c), c4(c);
				c2.d[long_dim][0] = hi_pos - 0.06*cube_width;
				c3.d[long_dim][1] = lo_pos + 0.06*cube_width;
				c4.d[!long_dim][!dir] = wall_pos + dir_sign*0.12*cube_depth;
				objs.emplace_back(c2, TYPE_COLLIDER, room_id, !long_dim, dir, RO_FLAG_INVIS, tot_light_amt); // side1
				objs.emplace_back(c3, TYPE_COLLIDER, room_id, !long_dim, dir, RO_FLAG_INVIS, tot_light_amt); // side2
				objs.emplace_back(c4, TYPE_COLLIDER, room_id, !long_dim, dir, RO_FLAG_INVIS, tot_light_amt); // back (against wall)

				if (has_office_chair && (rgen.rand()&3)) { // add office chair 75% of the time
					point center(c.get_cube_center());
					center[!long_dim] += dir_sign*0.2*cube_depth;
					for (unsigned d = 0; d < 2; ++d) {center[d] += 0.15*chair_radius*rgen.signed_rand_float();} // slightly random XY position
					center.z = zval;
					cube_t const chair(get_cube_height_radius(center, chair_radius, chair_height));
					objs.emplace_back(chair, TYPE_OFFICE_CHAIR, room_id, !long_dim, dir, RO_FLAG_RAND_ROT, tot_light_amt, room_obj_shape::SHAPE_CUBE, GRAY_BLACK);
				}
			} // for d
		} // for col
		lo_pos = hi_pos;
	} // for n
	return added_cube;
}

bool building_t::can_be_bedroom_or_bathroom(room_t const &room, bool on_first_floor) const { // check room type and existence of exterior door
	if (room.has_stairs || room.has_elevator || room.is_hallway || room.is_office) return 0; // no bed/bath in these cases (assumes a house)
	if (on_first_floor && is_room_adjacent_to_ext_door(room)) return 0; // door to house does not open into a bedroom/bathroom
	return 1;
}
bool building_t::can_be_bathroom(room_t const &room) const { // Note: assumes caller has checked can_be_bedroom_or_bathroom()
	float const vspace(get_window_vspace());
	return (min(room.dx(), room.dy()) < 2.0*vspace && max(room.dx(), room.dy()) < 3.0*vspace && count_num_int_doors(room) == 1);
}

unsigned building_t::count_num_int_doors(room_t const &room) const {
	cube_t room_exp(room);
	float const wall_thickness(get_wall_thickness());
	room_exp.expand_by(wall_thickness, wall_thickness, -wall_thickness); // expand in XY and shrink in Z
	unsigned num(0);
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {num += i->intersects(room_exp);}
	return num;
}

bool building_t::check_valid_closet_placement(cube_t const &c, room_t const &room, unsigned objs_start, float min_bed_space) const {
	if (min_bed_space > 0.0) {
		vector<room_object_t> &objs(interior->room_geom->objs);
		assert(objs_start < objs.size());
		room_object_t const &bed(objs[objs_start]);
		assert(bed.type == TYPE_BED);
		cube_t bed_exp(bed);
		bed_exp.expand_by_xy(min_bed_space);
		if (c.intersects_xy(bed_exp)) return 0; // too close to bed
	}
	return (!overlaps_other_room_obj(c, objs_start) && !is_cube_close_to_doorway(c, room, 0.0, 1));
}

bool building_t::add_bedroom_objs(rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers,
	float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool room_is_lit)
{
	vector<room_object_t> &objs(interior->room_geom->objs);
	unsigned const bed_obj_ix(objs.size()); // if placed, it will be this index
	if (!add_bed_to_room(rgen, room, blockers, zval, room_id, tot_light_amt)) return 0; // it's only a bedroom if there's bed
	assert(bed_obj_ix < objs.size());
	room_object_t const bed(objs[bed_obj_ix]); // deep copy so that we don't need to worry about invalidating the reference below
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.1*wall_thickness); // shrink to leave a small gap
	// closet
	float const doorway_width(interior->get_doorway_width()), floor_thickness(get_floor_thickness()), front_clearance(max(0.6f*doorway_width, get_min_front_clearance()));
	float const closet_min_depth(0.65*doorway_width), closet_min_width(1.5*doorway_width), min_dist_to_wall(1.0*doorway_width), min_bed_space(front_clearance);
	unsigned const first_corner(rgen.rand() & 3);
	bool const first_dim(rgen.rand_bool());
	cube_t const part(get_part_for_room(room));
	bool is_ext_wall[2][2] = {0}; // precompute which walls are exterior, {dim}x{dir}
	for (unsigned d = 0; d < 4; ++d) {is_ext_wall[d>>1][d&1] = (classify_room_wall(room, zval, (d>>1), (d&1), 0) == ROOM_WALL_EXT);}
	bool placed_closet(0);

	for (unsigned n = 0; n < 4 && !placed_closet; ++n) { // try 4 room corners
		unsigned const corner_ix((first_corner + n)&3);
		bool const xdir(corner_ix&1), ydir(corner_ix>>1);
		point const corner(room_bounds.d[0][xdir], room_bounds.d[1][ydir], zval);

		for (unsigned d = 0; d < 2 && !placed_closet; ++d) { // try both dims
			bool const dim(bool(d) ^ first_dim), dir(dim ? ydir : xdir), other_dir(dim ? xdir : ydir);
			if (room_bounds.get_sz_dim(!dim) < closet_min_width + min_dist_to_wall) continue; // room is too narrow to add a closet here
			if (is_ext_wall[dim][dir]) continue; // don't place closets against exterior walls where they would block a window
			float const dir_sign(dir ? -1.0 : 1.0), signed_front_clearance(dir_sign*front_clearance);
			float const window_hspacing(get_hspacing_for_part(part, dim));
			cube_t c(corner, corner);
			c.d[0][!xdir] += (xdir ? -1.0 : 1.0)*(dim ? closet_min_width : closet_min_depth);
			c.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(dim ? closet_min_depth : closet_min_width);
			if (is_ext_wall[!dim][other_dir] && is_val_inside_window(part, dim, c.d[dim][!dir], window_hspacing, get_window_h_border())) continue; // check for window intersection
			c.z2() += window_vspacing - floor_thickness;
			c.d[dim][!dir] += signed_front_clearance; // extra padding in front, to avoid placing too close to bed
			if (!check_valid_closet_placement(c, room, objs_start, min_bed_space)) continue; // bad placement
			// good placement, see if we can make the closet larger
			unsigned const num_steps = 10;
			float const req_dist(is_ext_wall[!dim][!other_dir] ? (other_dir ? -1.0 : 1.0)*min_dist_to_wall : 0.0); // signed; at least min dist from the opposite wall if it's exterior
			float const max_grow((room_bounds.d[!dim][!other_dir] - req_dist) - c.d[!dim][!other_dir]);
			float const len_step(max_grow/num_steps), depth_step(dir_sign*0.35*doorway_width/num_steps); // signed

			for (unsigned s1 = 0; s1 < num_steps; ++s1) { // try increasing width
				cube_t c2(c);
				c2.d[!dim][!other_dir] += len_step;
				if (!check_valid_closet_placement(c2, room, objs_start, min_bed_space)) break; // bad placement
				c = c2; // valid placement, update with larger cube
			}
			for (unsigned s2 = 0; s2 < num_steps; ++s2) { // now try increasing depth
				cube_t c2(c);
				c2.d[dim][!dir] += depth_step;
				if (is_ext_wall[!dim][other_dir] && is_val_inside_window(part, dim, (c2.d[dim][!dir] - signed_front_clearance), window_hspacing, get_window_h_border())) break; // bad placement
				if (!check_valid_closet_placement(c2, room, objs_start, min_bed_space)) break; // bad placement
				c = c2; // valid placement, update with larger cube
			}
			c.d[ dim][!dir] -= signed_front_clearance; // subtract off front clearance
			unsigned flags(0);
			if (c.d[!dim][0] == room_bounds.d[!dim][0]) {flags |= RO_FLAG_ADJ_LO;}
			if (c.d[!dim][1] == room_bounds.d[!dim][1]) {flags |= RO_FLAG_ADJ_HI;}
			objs.emplace_back(c, TYPE_CLOSET, room_id, dim, !dir, flags, tot_light_amt);
			placed_closet = 1; // done
		} // for d
	} // for n
	// dresser
	float const ds_height(rgen.rand_uniform(0.26, 0.32)*window_vspacing), ds_depth(rgen.rand_uniform(0.20, 0.25)*window_vspacing), ds_width(rgen.rand_uniform(0.6, 0.9)*window_vspacing);
	vector3d const ds_sz_scale(ds_depth/ds_height, ds_width/ds_height, 1.0);
	place_obj_along_wall(TYPE_DRESSER, room, ds_height, ds_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0);
	// nightstand
	unsigned const pref_orient(2*bed.dim + (!bed.dir)); // prefer the same orient as the bed so that it's placed on the same wall next to the bed
	float const ns_height(rgen.rand_uniform(0.24, 0.26)*window_vspacing), ns_depth(rgen.rand_uniform(0.15, 0.2)*window_vspacing), ns_width(rgen.rand_uniform(1.0, 2.0)*ns_depth);
	vector3d const ns_sz_scale(ns_depth/ns_height, ns_width/ns_height, 1.0);
	place_obj_along_wall(TYPE_DRESSER, room, ns_height, ns_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, pref_orient);

	// try to place a lamp on a dresser or nightstand that was added to this room
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_LAMP) && (rgen.rand()&3) != 0) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_LAMP)); // L, W, H
		float const height(0.25*window_vspacing), width(height*0.5f*(sz.x + sz.y)/sz.z);
		point pillow_center(bed.get_cube_center());
		pillow_center[bed.dim] += (bed.dir ? 1.0 : -1.0)*0.5*bed.get_sz_dim(bed.dim); // adjust from bed center to near the pillow(s)
		int obj_id(-1);
		float dmin_sq(0.0);

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) { // choose the dresser or nightstand closest to be bed
			if (i->type != TYPE_DRESSER) continue; // not a dresser or nightstand
			float const dist_sq(p2p_dist_xy_sq(i->get_cube_center(), pillow_center));
			if (dmin_sq == 0.0 || dist_sq < dmin_sq) {obj_id = (i - objs.begin()); dmin_sq = dist_sq;}
		}
		if (obj_id >= 0) { // found a valid object to place this on
			room_object_t const &obj(objs[obj_id]);
			point center(obj.get_cube_center());
			center.z = obj.z2();
			cube_t lamp(get_cube_height_radius(center, 0.5*width, height));
			lamp.translate_dim((obj.dir ? 1.0 : -1.0)*0.1*width, obj.dim); // move slightly toward the front to avoid clipping through the wall
			float const shift_range(0.4f*(obj.get_sz_dim(!obj.dim) - width)), obj_center(obj.get_center_dim(!obj.dim)), targ_pos(pillow_center[!obj.dim]);
			float shift_val(0.0), dmin(0.0);

			for (unsigned n = 0; n < 4; ++n) { // generate several random positions on the top of the object and choose the one closest to the bed
				float const cand_shift(rgen.rand_uniform(-1.0, 1.0)*shift_range), dist(fabs((obj_center + cand_shift) - targ_pos));
				if (dmin == 0.0 || dist < dmin) {shift_val = cand_shift; dmin = dist;}
			}
			lamp.translate_dim(shift_val, !obj.dim);
			unsigned const NUM_COLORS = 6;
			colorRGBA const colors[NUM_COLORS] = {WHITE, GRAY_BLACK, BROWN, LT_BROWN, DK_BROWN, OLIVE};
			unsigned flags(RO_FLAG_NOCOLL); // no collisions, as an optimization since the player and AI can't get onto the dresser/nightstand anyway
			if (!room_is_lit && rgen.rand_bool()) {flags |= RO_FLAG_LIT;} // 50% chance of being lit if the room is dark
			objs.emplace_back(lamp, TYPE_LAMP, room_id, obj.dim, obj.dir, flags, tot_light_amt, room_obj_shape::SHAPE_CYLIN, colors[rgen.rand()%NUM_COLORS]); // Note: invalidates obj ref
		}
	}
	return 1; // success
}

// Note: must be first placed object
bool building_t::add_bed_to_room(rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, float zval, unsigned room_id, float tot_light_amt) {
	unsigned const NUM_COLORS = 8;
	colorRGBA const colors[NUM_COLORS] = {WHITE, WHITE, WHITE, LT_BLUE, LT_BLUE, PINK, PINK, LT_GREEN}; // color of the sheets
	cube_t room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace()), wall_thick(get_wall_thickness());
	bool const dim(room_bounds.dx() < room_bounds.dy()); // longer dim
	vector3d expand, bed_sz;
	expand[ dim] = -wall_thick; // small amount of space
	expand[!dim] = -0.3f*vspace; // leave at least some space between the bed and the wall
	room_bounds.expand_by_xy(expand);
	float const room_len(room_bounds.get_sz_dim(dim)), room_width(room_bounds.get_sz_dim(!dim));
	if (room_len < 1.3*vspace || room_width < 0.7*vspace) return 0; // room is too small to fit a bed
	if (room_len > 4.0*vspace || room_width > 2.5*vspace) return 0; // room is too large to be a bedroom
	bool const first_head_dir(rgen.rand_bool()), first_wall_dir(rgen.rand_bool());
	vector<room_object_t> &objs(interior->room_geom->objs);
	cube_t c;
	c.z1() = zval;

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a bed
		float const sizes[6][2] = {{38, 75}, {38, 80}, {53, 75}, {60, 80}, {76, 80}, {72, 84}}; // twin, twin XL, full, queen, king, cal king
		unsigned const size_ix((room_width < 0.9*vspace) ? (rgen.rand() % 6) : (2 + (rgen.rand() % 4))); // only add twin beds to narrow rooms
		bed_sz[ dim] = 0.01f*vspace*(sizes[size_ix][1] + 8.0f); // length (mattress + headboard + footboard)
		bed_sz[!dim] = 0.01f*vspace*(sizes[size_ix][0] + 4.0f); // width  (mattress + small gaps)
		if (room_bounds.dx() < 1.5*bed_sz.x || room_bounds.dy() < 1.5*bed_sz.y) continue; // room is too small for a bed of this size
		bed_sz.z = 0.3*vspace*rgen.rand_uniform(1.0, 1.2); // height
		c.z2()   = zval + bed_sz.z;

		for (unsigned d = 0; d < 2; ++d) {
			float const min_val(room_bounds.d[d][0]), max_val(room_bounds.d[d][1] - bed_sz[d]);

			if (bool(d) == dim && n < 5) { // in the first few iterations, try to place the head of the bed against the wall (maybe not for exterior wall facing window?)
				c.d[d][0] = ((first_head_dir ^ bool(n&1)) ? min_val : max_val);
			}
			else if (bool(d) != dim && rgen.rand_bool()) { // try to place the bed against the wall sometimes
				c.d[d][0] = ((first_wall_dir ^ bool(n&1)) ? (min_val - 0.25*vspace) : (max_val + 0.25*vspace));
			}
			else {
				c.d[d][0] = rgen.rand_uniform(min_val, max_val);
			}
			c.d[d][1] = c.d[d][0] + bed_sz[d];
		} // for d
		if (!is_valid_placement_for_room(c, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
		bool const dir((room_bounds.d[dim][1] - c.d[dim][1]) < (c.d[dim][0] - room_bounds.d[dim][0])); // head of the bed is closer to the wall
		objs.emplace_back(c, TYPE_BED, room_id, dim, dir, 0, tot_light_amt);
		room_object_t &bed(objs.back());
		bed.obj_id = (uint16_t)objs.size();
		// use white color if a texture is assigned that's not close to white
		int const sheet_tid(bed.get_sheet_tid());
		if (sheet_tid < 0 || sheet_tid == WHITE_TEX || texture_color(sheet_tid).get_luminance() > 0.5) {bed.color = colors[rgen.rand()%NUM_COLORS];}
		return 1; // done/success
	} // for n
	return 0;
}

bool building_t::place_obj_along_wall(room_object type, room_t const &room, float height, vector3d const &sz_scale, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
	cube_t const &place_area, unsigned objs_start, float front_clearance, unsigned pref_orient, bool pref_centered, colorRGBA const &color, bool not_at_window)
{
	float const hwidth(0.5*height*sz_scale.y/sz_scale.z), depth(height*sz_scale.x/sz_scale.z), min_space(2.8*hwidth);
	vector3d const place_area_sz(place_area.get_size());
	if (max(place_area_sz.x, place_area_sz.y) <= min_space) return 0; // can't fit in either dim
	unsigned const force_dim((place_area_sz.x <= min_space) ? 0 : ((place_area_sz.y <= min_space) ? 1 : 2)); // *other* dim; 2=neither
	float const obj_clearance(depth*front_clearance), clearance(max(obj_clearance, get_min_front_clearance()));
	vector<room_object_t> &objs(interior->room_geom->objs);
	cube_t c;
	c.z1() = zval;
	c.z2() = zval + height;
	bool center_tried[4] = {};

	for (unsigned n = 0; n < 25; ++n) { // make 25 attempts to place the object
		bool const use_pref(pref_orient < 4 && n < 10); // use pref orient for first 10 tries
		bool const dim((force_dim < 2) ? force_dim : (use_pref ? (pref_orient >> 1) : rgen.rand_bool())); // choose a random wall unless forced
		bool const dir(use_pref ? !(pref_orient & 1) : rgen.rand_bool()); // dir is inverted for the model, so we invert pref dir as well
		unsigned const orient(2*dim + dir);
		float center(0.0);
		if (pref_centered && !center_tried[orient]) {center = place_area.get_center_dim(!dim); center_tried[orient] = 1;} // try centered
		else {center = rgen.rand_uniform(place_area.d[!dim][0]+hwidth, place_area.d[!dim][1]-hwidth);} // random position
		c.d[ dim][ dir] = place_area.d[dim][dir];
		c.d[ dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		c.d[!dim][   0] = center - hwidth;
		c.d[!dim][   1] = center + hwidth;

		if (not_at_window && classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) {
			cube_t const part(get_part_for_room(room));
			float const hspacing(get_hspacing_for_part(part, !dim)), border(get_window_h_border());
			// assume object is no larger than 2x window size and check left, right, and center positions
			if (is_val_inside_window(part, !dim, c.d[!dim][0], hspacing, border) ||
				is_val_inside_window(part, !dim, c.d[!dim][1], hspacing, border) ||
				is_val_inside_window(part, !dim, c.get_center_dim(!dim), hspacing, border)) continue;
		}
		cube_t c2(c); // used for collision tests
		c2.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance;
		if (overlaps_other_room_obj(c2, objs_start) || interior->is_blocked_by_stairs_or_elevator(c2)) continue; // bad placement
		// we don't need clearance for both the door and the object; test the object itself against the open door and the object with clearance against the closed door
		if (is_cube_close_to_doorway(c, room, 0.0, 1)) continue; // bad placement
		cube_t c3(c); // used for collision tests
		c3.d[dim][!dir] += (dir ? -1.0 : 1.0)*obj_clearance; // smaller clearance value (without player diameter)
		if (is_cube_close_to_doorway(c3, room, 0.0, 0)) continue; // bad placement
		objs.emplace_back(c, type, room_id, dim, !dir, 0, tot_light_amt, room_obj_shape::SHAPE_CUBE, color);
		objs.back().obj_id = (uint16_t)objs.size();
		if (front_clearance > 0.0) {objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS);} // add blocker cube to ensure no other object overlaps this space
		return 1; // done
	} // for n
	return 0; // failed
}
bool building_t::place_model_along_wall(unsigned model_id, room_object type, room_t const &room, float height, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
	cube_t const &place_area, unsigned objs_start, float front_clearance, unsigned pref_orient, bool pref_centered, colorRGBA const &color, bool not_at_window)
{
	if (!building_obj_model_loader.is_model_valid(model_id)) return 0; // don't have a model of this type
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
	return place_obj_along_wall(type, room, height*get_window_vspace(), sz, rgen, zval, room_id, tot_light_amt,
		place_area, objs_start, front_clearance, pref_orient, pref_centered, color, not_at_window);
}

bool building_t::add_bathroom_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned floor) {
	// Note: zval passed by reference
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.5*wall_thickness);
	if (min(place_area.dx(), place_area.dy()) < 0.7*floor_spacing) return 0; // room is too small (should be rare)
	bool const have_toilet(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOILET)), have_sink(building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK));
	vector<room_object_t> &objs(interior->room_geom->objs);

	if (!is_house && (have_toilet || have_sink)) { // office with at least a toilet or sink - replace carpet with tile
		float const new_zval(zval + 0.05*wall_thickness);
		cube_t floor(room_bounds);
		floor.z1() = zval;
		floor.z2() = new_zval;
		objs.emplace_back(floor, TYPE_FLOORING, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt);
		objs_start = objs.size(); // exclude this from collision checks
		zval = new_zval; // move the effective floor up
	}
	if (have_toilet && room.is_office && min(place_area.dx(), place_area.dy()) > 1.5*floor_spacing && max(place_area.dx(), place_area.dy()) > 2.0*floor_spacing) {
		if (divide_bathroom_into_stalls(rgen, room, zval, room_id, tot_light_amt, floor)) return 1; // large enough, try to divide into bathroom stalls
	}
	bool placed_obj(0), placed_toilet(0);
	
	// place toilet first because it's in the corner out of the way and higher priority
	if (have_toilet) { // have a toilet model
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
		float const height(0.35*floor_spacing), width(height*sz.y/sz.z), length(height*sz.x/sz.z);
		unsigned const first_corner(rgen.rand() & 3);
		bool const first_dim(rgen.rand_bool());

		for (unsigned n = 0; n < 4 && !placed_toilet; ++n) { // try 4 room corners
			unsigned const corner_ix((first_corner + n)&3);
			bool const xdir(corner_ix&1), ydir(corner_ix>>1);
			point const corner(place_area.d[0][xdir], place_area.d[1][ydir], zval);

			for (unsigned d = 0; d < 2 && !placed_toilet; ++d) { // try both dims
				bool const dim(bool(d) ^ first_dim), dir(dim ? ydir : xdir);
				cube_t c(corner, corner);
				c.d[0][!xdir] += (xdir ? -1.0 : 1.0)*(dim ? width : length);
				c.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(dim ? length : width);
				for (unsigned e = 0; e < 2; ++e) {c.d[!dim][e] += ((dim ? xdir : ydir) ? -1.5 : 1.5)*wall_thickness;} // extra padding on left and right sides
				c.z2() += height;
				cube_t c2(c); // used for placement tests
				c2.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.8*length; // extra padding in front of toilet, to avoid placing other objects there (sink and tub)
				c2.expand_in_dim(!dim, 0.4*width); // more padding on the sides
				if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, room, 0.0, 1)) continue; // bad placement
				objs.emplace_back(c, TYPE_TOILET, room_id, dim, !dir, 0, tot_light_amt);
				objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker cube to ensure no other object overlaps this space
				placed_obj = placed_toilet = 1; // done
			} // for d
		} // for n
		if (!placed_toilet) { // if the toilet can't be placed in a corner, allow it to be placed anywhere; needed for small offices
			placed_obj |= place_model_along_wall(OBJ_MODEL_TOILET, TYPE_TOILET, room, 0.35, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8);
		}
	}
	if (is_house) { // place a tub, but not in office buildings; placed before the sink because it's the largest and the most limited in valid locations
		cube_t place_area_tub(room_bounds);
		place_area_tub.expand_by(-0.05*wall_thickness); // just enough to prevent z-fighting
		placed_obj |= place_model_along_wall(OBJ_MODEL_TUB, TYPE_TUB, room, 0.2, rgen, zval, room_id, tot_light_amt, place_area_tub, objs_start, 0.4);
	}
	if (place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6)) {
		placed_obj = 1;
		room_object_t const &sink((objs.back().type == TYPE_SINK) ? objs.back() : objs[objs.size()-2]); // find sink, skip blocker
		
		if (classify_room_wall(room, zval, sink.dim, !sink.dir, 0) != ROOM_WALL_EXT) { // interior wall only
			// add a mirror above the sink; could later make into medicine cabinet
			cube_t mirror(sink); // start with the sink left and right position
			mirror.expand_in_dim(!sink.dim, 0.1*mirror.get_sz_dim(!sink.dim)); // make slightly wider
			mirror.z1() = sink.z2();
			mirror.z2() = mirror.z1() + 0.3*floor_spacing;
			mirror.d[sink.dim][!sink.dir] = room_bounds.d[sink.dim][!sink.dir];
			mirror.d[sink.dim][ sink.dir] = mirror.d[sink.dim][!sink.dir] + (sink.dir ? 1.0 : -1.0)*1.0*wall_thickness; // thickness
			// this mirror is actually 3D, so we enable collision detection; treat as a house even if it's in an office building
			objs.emplace_back(mirror, TYPE_MIRROR, room_id, sink.dim, sink.dir, RO_FLAG_IS_HOUSE, tot_light_amt);
		}
	}
	return placed_obj;
}

bool building_t::divide_bathroom_into_stalls(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned floor) {
	// Note: assumes no prior placed objects
	bool const use_sink_model(0 && building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK)); // not using sink models
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	vector3d const tsz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
	float const theight(0.35*floor_spacing), twidth(theight*tsz.y/tsz.z), tlength(theight*tsz.x/tsz.z), stall_depth(2.2*tlength);
	float sheight(0), swidth(0), slength(0), uheight(0), uwidth(0), ulength(0);

	if (use_sink_model) {
		vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SINK)); // L, W, H
		sheight = 0.45*floor_spacing; swidth = sheight*ssz.y/ssz.z; slength = sheight*ssz.x/ssz.z;
	}
	else {
		sheight = 0.36*floor_spacing; swidth = 0.3*floor_spacing; slength = 0.32*floor_spacing;
	}
	float stall_width(2.0*twidth), sink_spacing(1.75*swidth);
	bool br_dim(room.dy() < room.dx()), sink_side(0), sink_side_set(0);
	cube_t place_area(room), br_door;
	place_area.expand_by(-0.5*wall_thickness);

	// determine men's room vs. women's room
	point const part_center(get_part_for_room(room).get_cube_center()), room_center(room.get_cube_center());
	bool mens_room((part_center.x < room_center.x) ^ (part_center.y < room_center.y)), has_second_bathroom(0);

	// if there are two bathrooms (one on each side of the building), assign a gender to each side; if only one, alternate gender per floor
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (r->part_id != room.part_id || &(*r) == &room) continue; // different part or same room
		if (is_room_office_bathroom(*r, zval)) {has_second_bathroom = 1; break;}
	}
	if (!has_second_bathroom) {mens_room ^= (floor & 1);}
	bool const add_urinals(mens_room && building_obj_model_loader.is_model_valid(OBJ_MODEL_URINAL));

	if (add_urinals) { // use urinal model
		vector3d const usz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_URINAL)); // L, W, H
		uheight = 0.4*floor_spacing; uwidth = uheight*usz.y/usz.z; ulength = uheight*usz.x/usz.z;
	}
	for (unsigned d = 0; d < 2 && !sink_side_set; ++d) {
		for (unsigned side = 0; side < 2 && !sink_side_set; ++side) {
			cube_t c(room);
			c.z1() = zval; // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
			c.z2() = zval + wall_thickness;
			c.d[!br_dim][!side] = c.d[!br_dim][side] + (side ? -1.0 : 1.0)*wall_thickness; // shrink to near zero area in this dim

			for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
				if ((i->dy() < i->dx()) == br_dim) continue; // door in wrong dim
				if (!is_cube_close_to_door(c, 0.0, 0, *i, 2, 1)) continue; // check both dirs, check_zval=1
				sink_side = side; sink_side_set = 1;
				place_area.d[!br_dim][side] += (sink_side ? -1.0 : 1.0)*(i->get_sz_dim(br_dim) - 0.25*swidth); // add sink clearance for the door to close
				br_door = *i;
				break; // sinks are on the side closest to the door
			}
		} // for side
		if (d == 0 && !sink_side_set) {br_dim ^= 1;} // door not found on long dim - R90 and try short dim
	} // for d
	assert(sink_side_set);
	float const room_len(place_area.get_sz_dim(!br_dim)), room_width(place_area.get_sz_dim(br_dim));
	float const sinks_len(0.4*room_len), stalls_len(room_len - sinks_len), req_depth(2.0f*max(stall_depth, slength));
	if (room_width < req_depth) return 0;
	unsigned const num_stalls(std::floor(stalls_len/stall_width)), num_sinks(std::floor(sinks_len/sink_spacing));
	if (num_stalls < 2 || num_sinks < 2) return 0; // not enough space for 2 stalls and 2 sinks
	stall_width  = stalls_len/num_stalls; // reclaculate to fill the gaps
	sink_spacing = sinks_len/num_sinks;
	bool const two_rows(room_width > 1.5*req_depth), skip_stalls_side(room_id & 1); // put stalls on a side consistent across floors
	float const sink_side_sign(sink_side ? 1.0 : -1.0), stall_step(sink_side_sign*stall_width), sink_step(-sink_side_sign*sink_spacing);
	float const floor_thickness(get_floor_thickness());
	colorRGBA const stall_color(0.75, 1.0, 0.9, 1.0); // blue-green
	vector<room_object_t> &objs(interior->room_geom->objs);

	for (unsigned dir = 0; dir < 2; ++dir) { // each side of the wall
		if (!two_rows && dir == (unsigned)skip_stalls_side) continue; // no stalls/sinks on this side
		// add stalls
		float const dir_sign(dir ? -1.0 : 1.0), wall_pos(place_area.d[br_dim][dir]), stall_from_wall(wall_pos + dir_sign*(0.5*tlength + wall_thickness));
		float stall_pos(place_area.d[!br_dim][!sink_side] + 0.5*stall_step);

		for (unsigned n = 0; n < num_stalls; ++n) {
			point center(stall_from_wall, stall_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t toilet(center, center), stall(toilet);
			toilet.expand_in_dim( br_dim, 0.5*tlength);
			toilet.expand_in_dim(!br_dim, 0.5*twidth);
			toilet.z2() += theight;
			stall.z2() = stall.z1() + floor_spacing - floor_thickness; // set stall height to room height
			stall.expand_in_dim(!br_dim, 0.5*stall_width);
			stall.d[br_dim][ dir] = wall_pos; // + wall_thickness?
			stall.d[br_dim][!dir] = wall_pos + dir_sign*stall_depth;
			
			if (!interior->is_cube_close_to_doorway(stall, room, 0.0, 1, 1)) { // skip if close to a door (for rooms with doors at both ends); inc_open=1
				objs.emplace_back(toilet, TYPE_TOILET, room_id, br_dim, !dir, 0, tot_light_amt);
				objs.emplace_back(stall,  TYPE_STALL,  room_id, br_dim,  dir, 0, tot_light_amt, SHAPE_CUBE, stall_color);
			}
			stall_pos += stall_step;
		} // for n
		if (add_urinals && dir == (unsigned)skip_stalls_side) continue; // no urinals and sinks are each on one side
		// add sinks
		float const sink_start(place_area.d[!br_dim][sink_side] + 0.5f*sink_step);
		float const sink_from_wall(wall_pos + dir_sign*(0.5f*slength + (use_sink_model ? wall_thickness : 0.0f)));
		float sink_pos(sink_start);
		cube_t sinks_bcube;

		for (unsigned n = 0; n < num_sinks; ++n) {
			point center(sink_from_wall, sink_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t sink(center, center);
			sink.expand_in_dim( br_dim, 0.5*slength);
			sink.z2() += sheight;

			if (use_sink_model) { // sink 3D model
				sink.expand_in_dim(!br_dim, 0.5*swidth);
				objs.emplace_back(sink, TYPE_SINK, room_id, br_dim, !dir, 0, tot_light_amt);
			}
			else { // flat basin sink
				sink.expand_in_dim(!br_dim, 0.5*fabs(sink_step)); // tile exactly with the adjacent sink
				objs.emplace_back(sink, TYPE_BRSINK, room_id, br_dim, !dir, 0, tot_light_amt);
			}
			sinks_bcube.assign_or_union_with_cube(sink);
			sink_pos += sink_step;
		} // for n
		if (add_urinals) { // add urinals opposite the sinks, using same spacing as sinks
			float const u_wall(place_area.d[br_dim][!dir]), u_from_wall(u_wall - dir_sign*(0.5*ulength + 0.01*wall_thickness));
			float u_pos(sink_start);
			cube_t sep_wall;
			sep_wall.z1() = zval + 0.15*uheight;
			sep_wall.z2() = zval + 1.25*uheight;
			sep_wall.d[br_dim][!dir] = u_wall;
			sep_wall.d[br_dim][ dir] = u_wall - dir_sign*0.25*floor_spacing;

			for (unsigned n = 0; n < num_sinks; ++n) {
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, 0, tot_light_amt, SHAPE_SHORT, stall_color);
				point center(u_from_wall, u_pos, (zval + 0.2*uheight));
				if (br_dim) {swap(center.x, center.y);} // R90 about z
				cube_t urinal(center, center);
				urinal.expand_in_dim( br_dim, 0.5*ulength);
				urinal.expand_in_dim(!br_dim, 0.5*uwidth);
				urinal.z2() += uheight;
				objs.emplace_back(urinal, TYPE_URINAL, room_id, br_dim, dir, 0, tot_light_amt);
				u_pos += sink_step;
			} // for n
			if (!two_rows) { // skip first wall if adjacent to a stall
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, 0, tot_light_amt, SHAPE_SHORT, stall_color);
			}
		}
		if (!sinks_bcube.is_all_zeros()) { // add a long mirror above the sink
			if (!ENABLE_MIRROR_REFLECTIONS || dir != (unsigned)skip_stalls_side) { // don't add mirrors to both sides if reflections are enabled
				cube_t mirror(sinks_bcube);
				mirror.expand_in_dim(!br_dim, -0.25*wall_thickness); // slightly smaller
				mirror.d[br_dim][ dir] = wall_pos;
				mirror.d[br_dim][!dir] = wall_pos +  + dir_sign*0.1*wall_thickness;
				mirror.z1() = sinks_bcube.z2() + 0.25*floor_thickness;
				mirror.z2() = zval + 0.9*floor_spacing - floor_thickness;
				if (mirror.is_strictly_normalized()) {objs.emplace_back(mirror, TYPE_MIRROR, room_id, br_dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);}
			}
		}
	} // for dir
	// add a sign outside the bathroom door
	bool const shift_dir(room_center[br_dim] < part_center[br_dim]); // put the sign toward the outside of the building because there's more space and more light
	float const door_width(br_door.get_sz_dim(br_dim));
	cube_t sign(br_door);
	sign.z1() = zval + 0.50*floor_spacing;
	sign.z2() = zval + 0.55*floor_spacing;
	sign.translate_dim((shift_dir ? -1.0 : 1.0)*0.8*door_width, br_dim);
	sign.expand_in_dim(br_dim, -(mens_room ? 0.35 : 0.25)*door_width); // shrink a bit
	sign.translate_dim(sink_side_sign*0.5*wall_thickness, !br_dim); // move to outside wall
	sign.d[!br_dim][sink_side] += sink_side_sign*0.1*wall_thickness; // make nonzero area
	objs.emplace_back(sign, TYPE_SIGN, room_id, !br_dim, sink_side, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, DK_BLUE); // technically should use hallway room_id
	string const sign_text(mens_room ? "Men" : "Women");
	objs.back().obj_id = register_sign_text(sign_text);
	return 1;
}

void add_door_if_blocker(cube_t const &door, cube_t const &room, bool inc_open, bool dir, vect_cube_t &blockers) {
	bool const dim(door.dy() < door.dx()), edir(dim ^ dir ^ 1);
	float const width(door.get_sz_dim(!dim));
	cube_t door_exp(door);
	door_exp.expand_in_dim(dim, width);
	if (!door_exp.intersects(room)) return; // check against room before expanding along wall to exclude doors in adjacent rooms
	door_exp.expand_in_dim(!dim, width*0.25); // min expand value
	if (inc_open) {door_exp.d[!dim][edir] += (edir ? 1.0 : -1.0)*0.75*width;} // expand the remainder of the door width in this dir
	blockers.push_back(door_exp);
}
void building_t::gather_room_placement_blockers(cube_t const &room, unsigned objs_start, vect_cube_t &blockers, bool inc_open_doors, bool ignore_chairs) const {
	assert(has_room_geom());
	vector<room_object_t> &objs(interior->room_geom->objs);
	assert(objs_start <= objs.size());
	bool const first_floor(room.z1() < bcube.z1() + get_floor_thickness());
	blockers.clear();

	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (ignore_chairs && i->type == TYPE_CHAIR) continue;
		if (!(i->flags & RO_FLAG_NOCOLL) && i->intersects(room)) {blockers.push_back(*i);}
	}
	for (auto i = doors.begin(); i != doors.end(); ++i) {add_door_if_blocker(i->get_bcube(), room, 0, 0, blockers);} // exterior doors, inc_open=0
	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {add_door_if_blocker(*i, room, door_opens_inward(*i, room), i->open_dir, blockers);} // interior doors
	float const doorway_width(interior->get_doorway_width());

	for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
		cube_t tc(*s);
		if (first_floor) {tc.d[s->dim][!s->dir] += (s->dir ? -1.0 : 1.0);} // first floor, expand only in stairs entrance dim (could do the opposite for top floor)
		else {tc.expand_in_dim(s->dim, doorway_width);} // add extra space at both ends of stairs
		if (tc.intersects(bcube)) {blockers.push_back(tc);}
	}
	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		cube_t tc(*e);
		tc.d[e->dim][e->dir] += doorway_width*(e->dir ? 1.0 : -1.0); // add extra space in front of the elevator
		if (tc.intersects(bcube)) {blockers.push_back(tc);}
	}
}

bool building_t::add_kitchen_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool allow_adj_ext_door) {
	// Note: table and chairs have already been placed
	if (room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be kitchens
	if (!is_house && rgen.rand_bool()) return 0; // some office buildings have kitchens, allow it half the time
	// if it has an external door then reject the room half the time; most houses don't have a front door to the kitchen
	if (is_room_adjacent_to_ext_door(room, 1) && (!allow_adj_ext_door || rgen.rand_bool())) return 0; // front_door_only=1
	float const wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall for appliances
	bool placed_obj(0);
	placed_obj |= place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, room, 0.75, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.2, 4, 0, WHITE, 1); // not at window
	if (is_house) {placed_obj |= place_model_along_wall(OBJ_MODEL_STOVE, TYPE_STOVE, room, 0.46, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0);}
		
	if (is_house && placed_obj) { // if we have at least a fridge or stove, try to add countertops
		float const vspace(get_window_vspace()), height(0.345*vspace), depth(0.74*height), min_hwidth(0.6*height);
		float const front_clearance(max(0.6f*height, get_min_front_clearance()));
		cube_t cabinet_area(room_bounds);
		cabinet_area.expand_by(-0.05*wall_thickness); // smaller gap than place_area
		vector<room_object_t> &objs(interior->room_geom->objs);
		unsigned const counters_start(objs.size());
		cube_t c;
		c.z1() = zval;
		c.z2() = zval + height;
		cabinet_area.z1() = zval;
		cabinet_area.z2() = zval + vspace - get_floor_thickness();
		static vect_cube_t blockers;
		gather_room_placement_blockers(cabinet_area, objs_start, blockers, 1, 1); // inc_open_doors=1, ignore_chairs=1
		bool is_sink(1);

		for (unsigned n = 0; n < 50; ++n) { // 50 attempts
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
			float const center(rgen.rand_uniform(cabinet_area.d[!dim][0]+min_hwidth, cabinet_area.d[!dim][1]-min_hwidth)); // random position
			float const wall_pos(cabinet_area.d[dim][dir]), front_pos(wall_pos + (dir ? -1.0 : 1.0)*depth);
			c.d[ dim][ dir] = wall_pos;
			c.d[ dim][!dir] = front_pos + (dir ? -1.0 : 1.0)*front_clearance;
			c.d[!dim][   0] = center - min_hwidth;
			c.d[!dim][   1] = center + min_hwidth;
			cube_t c_min(c); // min runlength - used for collision tests
			for (unsigned e = 0; e < 2; ++e) {c.d[!dim][e] = cabinet_area.d[!dim][e];} // start at full room width
			bool bad_place(0);

			for (auto i = blockers.begin(); i != blockers.end(); ++i) {
				if (!i->intersects(c)) continue; // optimization - no cube interaction
				if (i->intersects(c_min)) {bad_place = 1; break;}
				if (i->d[!dim][1] < c_min.d[!dim][0]) {max_eq(c.d[!dim][0], i->d[!dim][1]);} // clip on lo side
				if (i->d[!dim][0] > c_min.d[!dim][1]) {min_eq(c.d[!dim][1], i->d[!dim][0]);} // clip on hi side
			}
			if (bad_place) continue;
			assert(c.contains_cube(c_min));
			c.d[dim][!dir] = front_pos; // remove front clearance

			for (auto i = objs.begin()+counters_start; i != objs.end(); ++i) { // find adjacencies to previously placed counters and flag to avoid placing doors
				if (i->dim == dim) continue; // not perpendicular
				if (i->d[!i->dim][dir] != wall_pos) continue; // not against the wall on this side
				if (i->d[i->dim][i->dir] != c.d[!dim][0] && i->d[i->dim][i->dir] != c.d[!dim][1]) continue; // not adjacent
				i->flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
			}
			objs.emplace_back(c, (is_sink ? TYPE_KSINK : TYPE_COUNTER), room_id, dim, !dir, 0, tot_light_amt);

			if (1) { // add upper cabinets, always (for now); should we remove cabinets in front of windows?
				cube_t c2(c);
				c2.z1() = zval + 0.66*vspace;
				c2.z2() = cabinet_area.z2(); // up to the ceiling
				if (!has_bcube_int_no_adj(c2, blockers)) {objs.emplace_back(c2, TYPE_CABINET, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);} // no collision detection
			}
			blockers.push_back(c); // add to blockers so that later counters don't intersect this one
			is_sink = 0; // sink is in first placed counter only
		} // for n
	}
	return placed_obj;
}

bool building_t::add_livingroom_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!is_house || room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be living rooms
	float const wall_thickness(get_wall_thickness());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall for appliances
	vector<room_object_t> &objs(interior->room_geom->objs);
	bool placed_couch(0), placed_tv(0);
	// place couches with a variety of colors
	unsigned const NUM_COLORS = 8;
	colorRGBA const colors[NUM_COLORS] = {GRAY_BLACK, WHITE, LT_GRAY, GRAY, DK_GRAY, LT_BROWN, BROWN, DK_BROWN};
	colorRGBA const &couch_color(colors[rgen.rand()%NUM_COLORS]);
	unsigned tv_pref_orient(4), couch_ix(objs.size()), tv_ix(0);
	
	if (place_model_along_wall(OBJ_MODEL_COUCH, TYPE_COUCH, room, 0.40, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.67, 4, 1, couch_color)) { // pref centered
		placed_couch   = 1;
		tv_pref_orient = (2*objs[couch_ix].dim + !objs[couch_ix].dir); // TV should be across from couch
	}
	tv_ix = objs.size();

	// place TV: pref centered; maybe should set not_at_window=1, but that seems too restrictive
	if (place_model_along_wall(OBJ_MODEL_TV, TYPE_TV, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 4.0, tv_pref_orient, 1, BKGRAY, 0)) {
		placed_tv = 1;
		// add a small table to place the TV on so that it's off the floor and not blocked as much by tables and chairs
		room_object_t &tv(objs[tv_ix]);
		float const height(0.4*tv.dz());
		cube_t table(tv); // same XY bounds as the TV
		tv.translate_dim(height, 2); // move TV up
		table.z2() = tv.z1();
		objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, 0, tot_light_amt);
	}
	if (placed_couch && placed_tv) {
		room_object_t const &couch(objs[couch_ix]), &tv(objs[tv_ix]);

		if (couch.dim == tv.dim && couch.dir != tv.dir) { // placed against opposite walls facing each other
			cube_t region(couch);
			region.union_with_cube(tv);
			shorten_chairs_in_region(region, objs_start); // region represents that space between the couch and the TV
		}
	}
	return (placed_couch || placed_tv);
}

bool building_t::add_library_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be libraries
	unsigned num_added(0);

	for (unsigned n = 0; n < 8; ++n) { // place up to 8 bookcases
		bool const added(add_bookcase_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start));
		if (added) {++num_added;} else {break;}
	}
	return (num_added > 0);
}

bool building_t::add_storage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), floor_thickness(get_floor_thickness());
	float const ceil_zval(zval + window_vspacing - floor_thickness), shelf_width(0.2*window_vspacing);
	cube_t room_bounds(get_walkable_room_bounds(room)), crate_bounds(room_bounds);
	vector<room_object_t> &objs(interior->room_geom->objs);
	unsigned const num_crates(4 + (rgen.rand() % 30)); // 4-33
	vect_cube_t exclude;
	cube_t test_cube(room);
	test_cube.z1() = zval; // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
	test_cube.z2() = zval + wall_thickness;
	unsigned num_placed(0);

	for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
		if (!is_cube_close_to_door(test_cube, 0.0, 0, *i, 2, 1)) continue; // check both dirs, check_zval=1
		exclude.push_back(*i);
		exclude.back().expand_in_dim( i->dim, 0.6*room.get_sz_dim(i->dim));
		exclude.back().expand_in_dim(!i->dim, 0.1*i->get_sz_dim(!i->dim));
	}
	// add shelves on walls (avoiding the door), and have crates avoid them
	for (unsigned dim = 0; dim < 2; ++dim) {
		if (room_bounds.get_sz_dim(dim) < 6.0*shelf_width) continue; // too narrow to add shelves in this dim

		for (unsigned dir = 0; dir < 2; ++dir) {
			if (rgen.rand_bool()) continue; // only add shelves to 50% of the walls
			cube_t shelves(room_bounds);
			shelves.z1() = zval;
			shelves.z2() = ceil_zval - floor_thickness;
			crate_bounds.d[dim][dir] = shelves.d[dim][!dir] = shelves.d[dim][dir] + (dir ? -1.0 : 1.0)*shelf_width; // outer edge of shelves, which is also the crate bounds
			shelves.expand_in_dim(!dim, -(shelf_width + 1.0f*wall_thickness)); // shorten shelves
			if (has_bcube_int(shelves, exclude)) continue; // too close to a doorway
			objs.emplace_back(shelves, TYPE_SHELVES, room_id, dim, dir, 0, tot_light_amt);
			objs.back().obj_id = (uint16_t)objs.size();
		} // for dir
	} // for dim
	// add a random office chair if there's space
	if (min(crate_bounds.dx(), crate_bounds.dy()) > 1.2*window_vspacing && building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR)) {
		float const chair_height(0.425*window_vspacing), chair_radius(0.5f*chair_height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR));
		cube_t place_area(crate_bounds);
		point pos(0.0, 0.0, zval);
		place_area.expand_by_xy(-chair_radius);
		for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(place_area.d[d][0], place_area.d[d][1]);}
		cube_t chair(get_cube_height_radius(pos, chair_radius, chair_height));
		
		// for now, just make one random attempt; if it fails then there's no chair in this room
		if (!has_bcube_int(chair, exclude) && !is_cube_close_to_doorway(chair, room, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(chair)) {
			objs.emplace_back(chair, TYPE_OFFICE_CHAIR, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_RAND_ROT, tot_light_amt);
		}
	}
	for (unsigned n = 0; n < 4*num_crates; ++n) { // make up to 4 attempts for every crate
		point pos(0.0, 0.0, zval);
		vector3d sz; // half size relative to window_vspacing
		for (unsigned d = 0; d < 3; ++d) {sz[d] = 0.06*window_vspacing*(1.0 + ((d == 2) ? 1.2 : 2.0)*rgen.rand_float());} // slightly more variation in XY
		cube_t place_area(crate_bounds);
		place_area.expand_by_xy(-sz);
		if (!place_area.is_strictly_normalized()) continue; // too large for this room
		for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(place_area.d[d][0], place_area.d[d][1]);}
		cube_t crate(get_cube_height_radius(pos, sz, 2.0*sz.z)); // multiply by 2 since this is a size rather than half size/radius
		if (has_bcube_int(crate, exclude)) continue; // don't place crates between the door and the center of the room
		bool bad_placement(0);

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (!i->intersects(crate)) continue;
			// only handle stacking of crates on other crates
			if (i->type == TYPE_CRATE && i->z1() == zval && (i->z2() + crate.dz() < ceil_zval) && i->contains_pt_xy(pos)) {crate.translate_dim(i->dz(), 2);}
			else {bad_placement = 1; break;}
		}
		if (bad_placement) continue;
		if (is_cube_close_to_doorway(crate, room, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(crate)) continue;
		objs.emplace_back(crate, TYPE_CRATE, room_id, 0, 0, 0, tot_light_amt);
		objs.back().obj_id = (uint16_t)objs.size(); // used to select texture
		if (++num_placed == num_crates) break; // we're done
	} // for n
	return 1; // it's always a storage room, even if it's empty
}

void building_t::place_book_on_obj(rand_gen_t rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, bool use_dim_dir) {
	point center(place_on.get_cube_center());
	for (unsigned d = 0; d < 2; ++d) {center[d] += 0.1*place_on.get_sz_dim(d)*rgen.rand_uniform(-1.0, 1.0);} // add a slight random shift
	float const book_sz(0.07*get_window_vspace());
	// book is randomly oriented for tables and rotated 90 degrees from desk orient
	bool const dim(use_dim_dir ? !place_on.dim : rgen.rand_bool()), dir(use_dim_dir ? (place_on.dir^place_on.dim) : rgen.rand_bool());
	cube_t book;
	vector3d book_scale(book_sz*rgen.rand_uniform(0.8, 1.2), book_sz*rgen.rand_uniform(0.8, 1.2), 0.0);
	book_scale[dim] *= 0.8; // slightly smaller in this dim
	book.set_from_point(point(center.x, center.y, place_on.z2()));
	book.expand_by(book_scale);
	book.z2() += book_sz*rgen.rand_uniform(0.1, 0.3);
	colorRGBA const color(book_colors[rgen.rand() % NUM_BOOK_COLORS]);
	vector<room_object_t> &objs(interior->room_geom->objs);
	objs.emplace_back(book, TYPE_BOOK, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, room_obj_shape::SHAPE_CUBE, color); // Note: invalidates place_on reference
	objs.back().obj_id = (uint16_t)objs.size();
}

void building_t::add_rug_to_room(rand_gen_t rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	if (!room_object_t::enable_rugs()) return; // disabled
	vector3d const room_sz(room.get_size());
	vector3d center(room.get_cube_center()); // Note: zvals ignored
	bool const min_dim(room_sz.y < room_sz.x);
	float const ar(rgen.rand_uniform(0.65, 0.85)), length(min(0.7f*room_sz[min_dim]/ar, room_sz[!min_dim]*rgen.rand_uniform(0.4, 0.7))), width(length*ar);
	cube_t rug;
	rug.z1() = zval;
	rug.z2() = rug.z1() + 0.001*get_window_vspace(); // almost flat

	for (unsigned d = 0; d < 2; ++d) {
		float const radius(0.5*((bool(d) == min_dim) ? width : length));
		center[d] += 0.05*room_sz[d]*rgen.rand_uniform(-1.0, 1.0); // slight random misalignment
		rug.d[d][0] = center[d] - radius;
		rug.d[d][1] = center[d] + radius;
	}
	vector<room_object_t> &objs(interior->room_geom->objs);
	objs.emplace_back(rug, TYPE_RUG, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt);
	objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix); // determines rug texture
}

// return value: 0=invalid, 1=valid and good, 2=valid but could be better
int building_t::check_valid_picture_placement(room_t const &room, cube_t const &c, float width, float zval, bool dim, bool dir, unsigned objs_start) const {
	float const wall_thickness(get_wall_thickness()), clearance(4.0*wall_thickness), side_clearance(1.0*wall_thickness);
	cube_t tc(c), keepout(c);
	tc.expand_in_dim(!dim, 0.1*width); // expand slightly to account for frame
	keepout.z1() = zval; // extend to the floor
	keepout.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance;
	keepout.expand_in_dim(!dim, side_clearance); // make sure there's space for the frame
	if (overlaps_other_room_obj(keepout, objs_start)) return 0;
	bool const inc_open(!is_house && !room.is_office);
	if (is_cube_close_to_doorway(tc, room, 0.0, inc_open)) return 0; // bad placement
	if ((room.has_stairs || room.has_elevator) && interior->is_blocked_by_stairs_or_elevator(tc, 4.0*wall_thickness)) return 0; // check stairs and elevators
	if (!inc_open && !room.is_hallway && is_cube_close_to_doorway(tc, room, 0.0, 1)) return 2; // success, but could be better (doors never open into hallway)
	return 1; // success
}

bool building_t::hang_pictures_in_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!room_object_t::enable_pictures()) return 0; // disabled
	
	if (!is_house && !room.is_office) {
		if (room.is_hallway) return 0; // no pictures or whiteboards in office building hallways (what about rooms with stairs?)
		// room in a commercial building - add whiteboard when there is a full wall to use
	}
	if (room.is_sec_bldg) return 0; // no pictures in secondary buildings
	if (room.rtype == RTYPE_STORAGE) return 0; // no pictures or whiteboards in storage rooms
	cube_t const &part(get_part_for_room(room));
	float const floor_height(get_window_vspace()), wall_thickness(get_wall_thickness());
	vector<room_object_t> &objs(interior->room_geom->objs);
	bool was_hung(0);

	if (!is_house || room.is_office) { // add whiteboards
		if (rgen.rand_float() < 0.2) return 0; // skip 20% of the time
		bool const pref_dim(rgen.rand_bool()), pref_dir(rgen.rand_bool());
		float const floor_thick(get_floor_thickness());

		for (unsigned dim2 = 0; dim2 < 2; ++dim2) {
			for (unsigned dir2 = 0; dir2 < 2; ++dir2) {
				bool const dim(bool(dim2) ^ pref_dim), dir(bool(dir2) ^ pref_dir);
				if (fabs(room.d[dim][dir] - part.d[dim][dir]) < 1.1*wall_thickness) continue; // on part boundary, likely exterior wall where there may be windows, skip
				cube_t c(room);
				c.z1() = zval + 0.25*floor_height;
				c.z2() = zval + 0.9*floor_height - floor_thick;
				c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*0.6*wall_thickness; // Note: offset by an additional half wall thickness
				c.expand_in_dim(!dim, -0.2*room.get_sz_dim(!dim)); // xy_space
				if (!check_valid_picture_placement(room, c, 0.6*room.get_sz_dim(!dim), zval, dim, dir, objs_start)) continue;
				objs.emplace_back(c, TYPE_WBOARD, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // whiteboard faces dir opposite the wall
				return 1; // done, only need to add one
			} // for dir
		} // for dim
		return 0;
	}
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			float const wall_pos(room.d[dim][dir]);
			if (fabs(room.d[dim][dir] - part.d[dim][dir]) < 1.1*wall_thickness) continue; // on part boundary, likely exterior wall where there may be windows, skip
			if (!room.is_hallway && rgen.rand_float() < 0.2) continue; // skip 20% of the time unless it's a hallway
			float const height(floor_height*rgen.rand_uniform(0.3, 0.6)), width(height*rgen.rand_uniform(1.5, 2.0)); // width > height
			if (width > 0.8*room.get_sz_dim(!dim)) continue; // not enough space
			point center;
			center[ dim] = wall_pos;
			center[!dim] = room.get_center_dim(!dim);
			center.z     = zval + rgen.rand_uniform(0.45, 0.55)*floor_height; // move up
			float const lo(room.d[!dim][0] + 0.7*width), hi(room.d[!dim][1] - 0.7*width);
			cube_t best_pos;

			for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to choose a position along the wall; first iteration is the center
				if (n > 0) { // try centered first, then non-centered
					if (hi - lo < width) break; // not enough space to shift, can't place this picture
					center[!dim] = rgen.rand_uniform(lo, hi);
				}
				cube_t c(center, center);
				c.expand_in_dim(2, 0.5*height);
				c.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.1*wall_thickness;
				if (room.is_hallway) {c.translate_dim((dir ? -1.0 : 1.0)*0.5*wall_thickness, dim);} // add an additional half wall thickness for hallways
				c.expand_in_dim(!dim, 0.5*width);
				int const ret(check_valid_picture_placement(room, c, width, zval, dim, dir, objs_start));
				if (ret == 0) continue; // invalid, retry
				best_pos = c;
				if (ret == 1) break; // valid and good - keep this pos
			} // for n
			if (best_pos.is_all_zeros()) continue; // failed placement
			objs.emplace_back(best_pos, TYPE_PICTURE, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // picture faces dir opposite the wall
			objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix + 61*dim + 123*dir); // determines picture texture
			was_hung = 1;
		} // for dir
	} // for dim
	return was_hung;
}

void building_t::add_plants_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num) {
	float const window_vspacing(get_window_vspace());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.1*get_wall_thickness()); // shrink to leave a small gap
	unsigned const num_colors = 8;
	colorRGBA const pot_colors[num_colors] = {LT_GRAY, GRAY, DK_GRAY, BKGRAY, WHITE, LT_BROWN, RED, colorRGBA(1.0, 0.35, 0.18)};
	
	for (unsigned n = 0; n < num; ++n) {
		float const height(rgen.rand_uniform(0.6, 0.9)*window_vspacing), width(rgen.rand_uniform(0.15, 0.35)*window_vspacing);
		vector3d const sz_scale(width/height, width/height, 1.0);
		place_obj_along_wall(TYPE_PLANT, room, height, sz_scale, rgen, zval, room_id, tot_light_amt,
			place_area, objs_start, 0.0, 4, 0, pot_colors[rgen.rand() % num_colors]); // no clearanc, pref_orient, or color
	}
}

void building_t::add_bathroom_windows(room_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	unsigned num_ext_walls(0);

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {num_ext_walls += (classify_room_wall(room, zval, dim, dir, 1) == ROOM_WALL_EXT);}
	}
	if (num_ext_walls != 1) return; // it looks odd to have window block walls at the corner of a building, so only enable this for single exterior walls
	float const floor_thickness(get_floor_thickness()), wall_thickness(get_wall_thickness()), window_thickness(0.05*wall_thickness);
	float const z2(zval + get_window_vspace() - floor_thickness);
	vector<room_object_t> &objs(interior->room_geom->objs);

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (classify_room_wall(room, zval, dim, dir, 1) != ROOM_WALL_EXT) continue; // exterior walls only
			cube_t c(room);
			c.z1() = zval; c.z2() = z2;
			c.d[dim][!dir]  = c.d[dim][dir] + (dir ? -1.0 : 1.0)*window_thickness;
			c.d[dim][ dir] -= (dir ? -1.0 : 1.0)*window_thickness;
			// shrink by wall thickness to avoid problems at the corners of buildings
			if (c.d[!dim][0] == bcube.d[!dim][0]) {c.d[!dim][0] += wall_thickness;}
			if (c.d[!dim][1] == bcube.d[!dim][1]) {c.d[!dim][1] -= wall_thickness;}
			objs.emplace_back(c, TYPE_WINDOW, room_id, dim, dir, RO_FLAG_NOCOLL, max(tot_light_amt, 1.0f), room_obj_shape::SHAPE_CUBE, WHITE); // always lit
		} // for dir
	} // for dim
}

void set_light_xy(cube_t &light, point const &center, float light_size, bool light_dim, room_obj_shape light_shape) {
	for (unsigned dim = 0; dim < 2; ++dim) {
		float const sz(((light_shape == SHAPE_CYLIN || light_shape == SHAPE_SPHERE) ? 1.6 : ((bool(dim) == light_dim) ? 2.2 : 1.0))*light_size);
		light.d[dim][0] = center[dim] - sz;
		light.d[dim][1] = center[dim] + sz;
	}
	light.z1() = light.z2() = center.z; // set so that valid pos can be checked
}

template<typename T> bool has_bcube_int_exp(cube_t const &bcube, vector<T> const &bcubes, float expand) {
	cube_t bcube_exp(bcube);
	bcube_exp.expand_by(expand); // expand in all dirs, including z
	return has_bcube_int(bcube_exp, bcubes);
}

// Note: these three floats can be calculated from mat.get_floor_spacing(), but it's easier to change the constants if we just pass them in
void building_t::gen_room_details(rand_gen_t &rgen, vect_cube_t const &ped_bcubes, unsigned building_ix) {

	assert(interior);
	if (interior->room_geom) return; // already generated?
	//highres_timer_t timer("Gen Room Details");
	interior->room_geom.reset(new building_room_geom_t(bcube.get_llc()));
	vector<room_object_t> &objs(interior->room_geom->objs);
	vector<room_t> &rooms(interior->rooms);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const light_thick(0.025*window_vspacing), def_light_size(0.1*window_vspacing);
	interior->room_geom->obj_scale = window_vspacing; // used to scale room object textures
	unsigned tot_num_rooms(0), num_light_stacks(0), num_bathrooms(0);
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {tot_num_rooms += calc_num_floors(*r, window_vspacing, floor_thickness);}
	objs.reserve(tot_num_rooms); // placeholder - there will be more than this many
	float const extra_bathroom_prob((is_house ? 2.0 : 1.0)*0.02*min((int(tot_num_rooms) - 4), 20));
	room_obj_shape const light_shape(is_house ? SHAPE_CYLIN : SHAPE_CUBE);
	unsigned cand_bathroom(rooms.size()); // start at an invalid value
	unsigned added_kitchen_mask(0); // per-floor
	bool added_bedroom(0), added_living(0), added_library(0);

	if (rooms.size() > 1) { // choose best room assignments for required rooms; if a single room, skip this step
		float min_score(0.0);

		for (auto r = rooms.begin(); r != rooms.end(); ++r) {
			if (r->is_sec_bldg) continue; // garage/shed excluded - not a normal room
			bool const on_first_floor(calc_num_floors(*r, window_vspacing, floor_thickness) == 1); // if a single floor, this must be on the first floor

			if (can_be_bedroom_or_bathroom(*r, on_first_floor)) { // find best bathroom with no hard size constraints
				float score(r->dx() + r->dy()); // starts as half the perimeter
				score *= (1.0 + 10.0*(max(count_num_int_doors(*r), 1U) - 1U)); // multiply by a large value if there are mult doors so we only choose this if there are no alternatives
				if (min_score == 0.0 || score < min_score) {cand_bathroom = (r - rooms.begin()); min_score = score;}
			}
		} // for r
	}
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {
		float const light_amt(window_vspacing*r->get_light_amt()); // multiply perimeter/area by window spacing to make unitless
		float const floor_height(r->is_sec_bldg ? r->dz() : window_vspacing); // secondary buildings are always one floor
		unsigned const num_floors(calc_num_floors(*r, floor_height, floor_thickness)), room_id(r - rooms.begin());
		point room_center(r->get_cube_center());
		// determine light pos and size for this stack of rooms
		bool const room_dim(r->dx() < r->dy()); // longer room dim
		bool const must_be_bathroom(room_id == cand_bathroom && num_bathrooms == 0); // cand bathroom, and bathroom not already placed
		float light_size(def_light_size); // default size for houses
		unsigned const room_objs_start(objs.size());

		if (r->is_sec_bldg) {
			if    (has_garage) {r->assign_to(RTYPE_GARAGE);}
			else if (has_shed) {r->assign_to(RTYPE_SHED);}
		}
		if (r->is_office) { // light size varies by office size
			float const room_size(r->dx() + r->dy()); // normalized to office size
			light_size = max(0.015f*room_size, 0.67f*def_light_size);
		}
		if (r->is_hallway) { // light size varies by hallway size
			float const room_size(min(r->dx(), r->dy())); // normalized to hallway width
			light_size = max(0.06f*room_size, 0.67f*def_light_size);
		}
		if (r->has_stairs && r->rtype == RTYPE_NOTSET) {r->assign_to(RTYPE_STAIRS);} // default to stairs, may be re-assigned below
		float const light_val(22.0*light_size);
		r->light_intensity = light_val*light_val/r->get_area_xy(); // average for room, unitless; light surface area divided by room surface area with some fudge constant
		cube_t pri_light, sec_light;
		set_light_xy(pri_light, room_center, light_size, room_dim, light_shape);
		if (!r->contains_cube_xy(pri_light)) {pri_light.set_to_zeros();} // disable light if it doesn't fit (small room)
		bool const blocked_by_stairs(!r->is_hallway && interior->is_blocked_by_stairs_or_elevator_no_expand(pri_light, fc_thick));
		bool use_sec_light(0), added_bathroom(0);
		float z(r->z1());

		if (blocked_by_stairs) { // blocked by stairs - see if we can add a light off to the side in the other orient
			bool const first_dir(rgen.rand_bool());

			for (unsigned d = 0; d < 2 && !use_sec_light; ++d) { // see if we can place it by moving on one direction
				for (unsigned n = 0; n < 5; ++n) { // try 5 different shift values: 0.2, 0.25, 0.3, 0.35, 0.4
					point new_center(room_center);
					new_center[room_dim] += ((bool(d) ^ first_dir) ? -1.0 : 1.0)*(0.2 + 0.05*n)*r->get_sz_dim(room_dim);
					set_light_xy(sec_light, new_center, light_size, !room_dim, light_shape); // flip the light dim
					if (!interior->is_blocked_by_stairs_or_elevator_no_expand(sec_light, fc_thick)) {use_sec_light = 1; break;} // add if not blocked
				}
			}
		}
		r->interior = get_part_for_room(*r).contains_cube_xy_no_adj(*r);
		// make chair colors consistent for each part by using a few variables for a hash
		colorRGBA chair_colors[12] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN};
		colorRGBA chair_color(chair_colors[(13*r->part_id + 123*tot_num_rooms + 617*mat_ix + 1367*num_floors) % 12]);
		unsigned num_lights_added(0);

		// place objects on each floor for this room
		for (unsigned f = 0; f < num_floors; ++f, z += floor_height) {
			room_center.z = z + fc_thick; // floor height
			bool const top_floor(f+1 == num_floors), check_stairs(!is_house && parts.size() > 1 && top_floor); // top floor of building that may have stairs connecting to upper stack
			bool is_lit(0), light_dim(room_dim), has_stairs(r->has_stairs);

			if (!has_stairs && (f == 0 || top_floor) && interior->stairwells.size() > 1) { // check for stairwells connecting stacked parts (is this still needed?)
				for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
					if (!r->contains_cube_xy(*s)) continue; // stairs not in this room
					// Note: here we adjust stairs zval by floor_thickness to include stairs in the floor but not in the room above
					if (s->z1() + floor_thickness > r->z2()) continue; // stairs above the room
					if (s->z2() + floor_thickness < r->z1()) continue; // stairs below the room
					has_stairs = 1;
				}
			}
			bool const top_of_stairs(has_stairs && top_floor);
			if (top_of_stairs || check_stairs) {num_light_stacks += num_lights_added;} // these cases may shift the light, so we allocate a new stack
			cube_t light;
			if (!blocked_by_stairs || top_of_stairs) {light = pri_light;}
			else if (use_sec_light) {light = sec_light; light_dim ^= 1;}
			int light_obj_ix(-1);

			if (!light.is_all_zeros()) { // add a light to the center of the ceiling of this room if there's space (always for top of stairs)
				light.z2() = z + floor_height - fc_thick;
				light.z1() = light.z2() - light_thick;
				is_lit = (r->is_hallway || ((rgen.rand() & (top_of_stairs ? 3 : 1)) != 0)); // 50% of lights are on, 75% for top of stairs, 100% for hallways

				// check ped_bcubes and set is_lit if any are people are in this floor of this room
				for (auto p = ped_bcubes.begin(); p != ped_bcubes.end() && !is_lit; ++p) {
					if (!p->intersects_xy(*r)) continue; // person not in this room
					if (p->z2() < light.z1() && p->z1() + floor_height > light.z2()) {is_lit = 1;} // on this floor
				}
				uint8_t flags(RO_FLAG_NOCOLL); // no collision detection with lights
				if (is_lit)        {flags |= RO_FLAG_LIT | RO_FLAG_EMISSIVE;}
				if (top_of_stairs) {flags |= RO_FLAG_TOS;}
				if (has_stairs)    {flags |= RO_FLAG_RSTAIRS;}
				colorRGBA color;
				if (is_house) {color = colorRGBA(1.0, 1.0, 0.9);} // house - yellowish
				else if (r->is_hallway || r->is_office) {color = colorRGBA(0.9, 0.9, 1.0);} // office building - blueish
				else {color = colorRGBA(1.0, 1.0, 1.0);} // white - small office
				unsigned num_lights(r->num_lights);

				if (r->is_hallway && num_lights > 1) { // hallway: place a light on each side (of the stairs if they exist), and also between stairs and elevator if there are both
					if (r->has_elevator && r->has_stairs) {num_lights = 3;} // we really should have 3 lights in this case
					float const offset(((num_lights == 3) ? 0.3 : 0.2)*r->get_sz_dim(light_dim)); // closer to the ends in the 3 lights case
					cube_t valid_bounds(*r);
					valid_bounds.expand_by_xy(-0.1*floor_height); // add some padding

					for (unsigned d = 0; d < num_lights; ++d) {
						float const delta((d == 2) ? 0.0 : (d ? -1.0 : 1.0)*offset); // last light is in the center
						cube_t hall_light(light);
						hall_light.translate_dim(delta, light_dim);

						if (check_stairs && has_bcube_int_exp(hall_light, interior->stairwells, fc_thick)) { // keep moving until not blocked by stairs
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
						} // end check_stairs
						objs.emplace_back(hall_light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape, color); // dir=0 (unused)
						objs.back().obj_id = num_light_stacks + d;
					} // for d
					num_lights_added = num_lights;
				}
				else if (r->is_office) { // office with possibly multiple lights
					float const dx(r->dx()), dy(r->dy()), ldx(light.dx()), ldy(light.dy());
					unsigned const nx(max(1U, unsigned(0.5*dx/window_vspacing))), ny(max(1U, unsigned(0.5*dy/window_vspacing))); // more lights for large offices
					float const xstep(dx/nx), ystep(dy/ny);
					vector3d const shrink(0.5*ldx*sqrt((nx - 1)/nx), 0.5*ldy*sqrt((ny - 1)/ny), 0.0);
					unsigned cur_light_ix(num_light_stacks);

					for (unsigned y = 0; y < ny; ++y) {
						for (unsigned x = 0; x < nx; ++x) {
							cube_t cur_light(light);
							cur_light.expand_by_xy(-shrink);
							cur_light.translate(point((-0.5f*dx + (x + 0.5)*xstep), (-0.5f*dy + (y + 0.5)*ystep), 0.0));
							if (check_stairs && has_bcube_int_exp(cur_light, interior->stairwells, fc_thick)) continue; // what about blocked_by_stairs flag?
							objs.emplace_back(cur_light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape, color); // dir=0 (unused)
							objs.back().obj_id = cur_light_ix++;
						} // for x
					} // for y
					num_lights_added = nx*ny;
				}
				else { // normal room with a single light
					if (check_stairs && has_bcube_int_exp(light, interior->stairwells, fc_thick)) {is_lit = 0;} // disable if blocked by stairs
					else {
						light_obj_ix = objs.size();
						objs.emplace_back(light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape, color); // dir=0 (unused)
						objs.back().obj_id = num_light_stacks;
					}
					num_lights_added = 1;
				}
				if (is_lit) {r->lit_by_floor |= (1ULL << (f&63));} // flag this floor as being lit (for up to 64 floors)
			} // end light placement
			float tot_light_amt(light_amt); // unitless, somewhere around 1.0
			if (is_lit) {tot_light_amt += r->light_intensity;}
			rgen.rand_mix();

			if (r->no_geom) {
				if (is_house && r->is_hallway) { // allow pictures and rugs in the hallways of houses
					hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs.size());
					if (rgen.rand_bool()) {add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt);}
				}
				continue; // no other geometry for this room
			}
			//if (has_stairs && !pri_hall.is_all_zeros()) continue; // no other geometry in office building base part rooms that have stairs
			unsigned const objs_start(objs.size()), floor_mask(1<<f);
			bool added_tc(0), added_desk(0), added_obj(0), can_place_book(0), is_bathroom(0), is_bedroom(0), is_kitchen(0), is_living(0), is_dining(0), is_storage(0), no_whiteboard(0);
			unsigned num_chairs(0);

			// place room objects
			bool const allow_br(!is_house || must_be_bathroom || f > 0 || num_floors == 1 || (rgen.rand_float() < 0.33f*(added_living + (added_kitchen_mask&1) + 1))); // bed/bath
			bool is_office_bathroom(is_room_office_bathroom(*r, room_center.z));

			if (is_office_bathroom) { // bathroom is already assigned
				added_obj = is_bathroom = added_bathroom = add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f); // add bathroom
			}
			if (!added_obj && allow_br && can_be_bedroom_or_bathroom(*r, (f == 0))) { // bedroom or bathroom case; need to check first floor even if is_cand_bathroom
				// place a bedroom 75% of the time unless this must be a bathroom; if we got to the second floor and haven't placed a bedroom, always place it; houses only
				if (is_house && !must_be_bathroom && ((f > 0 && !added_bedroom) || rgen.rand_float() < 0.75)) {
					added_obj = added_bedroom = is_bedroom = add_bedroom_objs(rgen, *r, ped_bcubes, room_center.z, room_id, tot_light_amt, objs_start, is_lit);
					if (is_bedroom) {r->assign_to(RTYPE_BED, f);}
					// Note: can't really mark room type as bedroom because it varies per floor; for example, there may be a bedroom over a living room connected to an exterior door
				}
				if (!added_obj && (must_be_bathroom || (can_be_bathroom(*r) && (num_bathrooms == 0 || rgen.rand_float() < extra_bathroom_prob)))) {
					// bathrooms can be in both houses and office buildings
					added_obj = is_bathroom = added_bathroom = add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f); // add bathroom
					if (is_bathroom) {r->assign_to(RTYPE_BATH, f);}
				}
			}
			if (!added_obj && r->is_office) { // add cubicles if this is a large office
				added_obj = no_whiteboard = create_office_cubicles(rgen, *r, room_center.z, room_id, tot_light_amt);
			}
			if (!added_obj && rgen.rand_float() < (r->is_office ? 0.6 : (is_house ? 0.95 : 0.5))) {
				// place a table and maybe some chairs near the center of the room if it's not a hallway;
				// 60% of the time for offices, 95% of the time for houses, and 50% for other buildings
				unsigned const num_tcs(add_table_and_chairs(rgen, *r, ped_bcubes, room_id, room_center, chair_color, 0.1, tot_light_amt));
				if (num_tcs > 0) {added_tc = added_obj = can_place_book = 1; num_chairs = num_tcs - 1;}
				// on ground floor, try to make this a kitchen; not all houses will have a kitchen with this logic - maybe we need fewer bedrooms?
				if (added_tc && !(added_kitchen_mask & floor_mask) && (!is_house || f == 0)) { // office buildings can also have kitchens, even on non-ground floors
					is_kitchen = add_kitchen_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, added_living);
					if (is_kitchen) {r->assign_to(RTYPE_KITCHEN, f); added_kitchen_mask |= floor_mask;}
				}
			}
			if (!added_obj && r->is_office && r->interior && f == 0 /*&& r->z1() == bcube.z1()*/ && rgen.rand_bool()) {
				// if we haven't added any objects yet, and this room is an interior office on the first floor, make it a storage room 50% of the time
				added_obj = no_whiteboard = is_storage = add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				if (added_obj) {r->assign_to(RTYPE_STORAGE, f);}
			}
			if (!added_obj) { // try to place a desk if there's no table, bed, etc.
				added_obj = can_place_book = added_desk = add_desk_to_room(rgen, *r, ped_bcubes, chair_color, room_center.z, room_id, tot_light_amt);
				if (added_obj && !r->has_stairs) {r->assign_to((is_house ? RTYPE_STUDY : RTYPE_OFFICE), f);} // or other room type - may be overwritten below
			}
			if (is_house && (added_tc || added_desk) && !is_kitchen && f == 0) { // don't add second living room unless we added a kitchen and have enough rooms
				if ((!added_living && rooms.size() >= 8 && (added_kitchen_mask || rgen.rand_bool())) || is_room_adjacent_to_ext_door(*r, 1)) { // front_door_only=1
					// add a living room on the ground floor if it has a table or desk but isn't a kitchen
					added_living = is_living = add_livingroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
					if (is_living) {r->assign_to(RTYPE_LIVING, f);}
				}
			}
			if (can_place_book) { // an object was placed (table, desk, counter, etc.), maybe add a book on top of it
				assert(objs.size() > objs_start);
				float const place_prob((is_house ? 1.0 : 0.5)*(r->is_office ? 0.8 : 1.0));
				unsigned const objs_end(objs.size());
				bool placed_on_counter(0);

				for (unsigned i = objs_start; i < objs_end; ++i) { // see if we can place books on any room objects
					room_object_t const &obj(objs[i]);
					bool place_book(0);
					if      (obj.type == TYPE_TABLE && i == objs_start) {place_book = (rgen.rand_float() < 0.40*place_prob);} // only first table (not TV table)
					else if (obj.type == TYPE_DESK) {place_book = ((i+1 == objs_end || objs[i+1].type != TYPE_TV) && rgen.rand_float() < 0.8*place_prob);} // only if no TV
					else if (obj.type == TYPE_COUNTER && !placed_on_counter) {place_book = placed_on_counter = (rgen.rand_float() < 0.5);}
					if (place_book) {place_book_on_obj(rgen, obj, room_id, tot_light_amt, (obj.type != TYPE_TABLE));}
				}
			}
			if (is_house) { // place house-specific items
				if (!is_bathroom && !is_kitchen && rgen.rand_float() < 0.8) { // place bookcase 80% of the time, but not in bathrooms or kitchens
					rand_gen_t rgen2(rgen); // copy so that rgen isn't updated in the call below
					add_bookcase_to_room(rgen2, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
				if (!has_stairs && (rgen.rand()&3) <= (added_tc ? 0 : 2)) { // maybe add a rug, 25% of the time if there's a table and 75% of the time otherwise
					add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt);
				}
			}
			if (is_house && added_tc && num_chairs > 0 && !is_living && !is_kitchen) { // room with table and chair that's not a kitchen
				if (f == 0) { // dining room, must be on the first floor
					if (light_obj_ix >= 0) { // handle dining room light: extend downward and make it a sphere
						assert((unsigned)light_obj_ix < objs.size());
						room_object_t &light(objs[light_obj_ix]);
						light.shape = SHAPE_SPHERE;
						light.z2() += 0.5f*light.dz();
						light.z1() -= 0.22f*(light.dx() + light.dy());
					}
					r->assign_to(RTYPE_DINING, f);
					is_dining = 1;
				}
				else if (!added_library && add_library_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start)) { // add library, at most one
					r->assign_to(RTYPE_LIBRARY, f);
					added_library = 1;
				}
			}
			if (r->rtype == RTYPE_NOTSET) {
				if (is_room_adjacent_to_ext_door(*r)) {r->assign_to(RTYPE_ENTRY, f);} // entryway if has exterior door and is unassigned
				else if (!is_house) {r->assign_to(RTYPE_OFFICE, f);} // any unset room in an office building is an office
			}
			if (!is_bathroom && !is_bedroom && !is_kitchen && !is_storage) { // add potted plants to some room types
				// 0-2 for living/dining rooms, 50% chance for houses, 25% (first floor) / 10% (other floors) chance for offices
				unsigned const num(is_house ? (rgen.rand() % ((is_living || is_dining) ? 3 : 2)) : ((rgen.rand()%((f == 0) ? 4 : 10)) == 0));
				if (num > 0) {add_plants_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, num);}
			}
			// pictures and whiteboards must not be placed behind anything, excluding trashcans; so we add them here
			bool const can_hang(is_house || !(is_bathroom || is_kitchen || no_whiteboard)); // no whiteboards in office bathrooms or kitchens
			bool const was_hung(can_hang && hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start));

			if (is_bathroom || is_kitchen || rgen.rand_float() < 0.8) { // 80% of the time, always in bathrooms and kitchens
				add_trashcan_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, (was_hung && !is_house)); // no trashcans on same wall as office whiteboard
			}
			if (is_bathroom) {add_bathroom_windows(*r, room_center.z, room_id, tot_light_amt);} // find all windows and add frosted windows
		} // for f (floor)
		num_light_stacks += num_lights_added;
		if (added_bathroom) {++num_bathrooms;}

		if (r->interior) { // tag objects as interior if room is interior
			for (auto i = objs.begin() + room_objs_start; i != objs.end(); ++i) {i->flags |= RO_FLAG_INTERIOR;}
		}
	} // for r (room)
	add_wall_and_door_trim();
	add_stairs_and_elevators(rgen); // the room objects - stairs and elevators have already been placed within a room
	add_exterior_door_signs(rgen);
	objs.shrink_to_fit();
	interior->room_geom->light_bcubes.resize(num_light_stacks); // allocate but don't fill un until needed
}

void building_t::add_wall_and_door_trim() { // and window trim

	//highres_timer_t timer("Add Wall And Door Trim");
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness), wall_thickness(get_wall_thickness());
	float const trim_height(0.04*window_vspacing), trim_thickness(0.1*wall_thickness), expand_val(2.0*trim_thickness);
	float const door_trim_exp(2.0*trim_thickness + 0.5*wall_thickness), door_trim_width(0.5*wall_thickness);
	float const door_height(get_door_height()), floor_to_ceil_height(window_vspacing - floor_thickness);
	unsigned const flags(RO_FLAG_NOCOLL);
	bool const has_ceil_trim(is_house); // ceiling trim on houses only, for now (though the code also handles office buildings)
	vector<room_object_t> &objs(interior->room_geom->objs);
	vect_cube_t trim_cubes;

	for (auto d = interior->doors.begin(); d != interior->doors.end(); ++d) { // vertical strips on each side + strip on top of interior doors
		cube_t trim(*d);
		trim.expand_in_dim(d->dim, door_trim_exp);

		for (unsigned side = 0; side < 2; ++side) { // left/right of door
			trim.d[!d->dim][0] = d->d[!d->dim][side] - (side ? door_trim_width : trim_thickness);
			trim.d[!d->dim][1] = d->d[!d->dim][side] + (side ? trim_thickness : door_trim_width);
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, d->dim, side, (flags | RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP), 1.0, SHAPE_TALL); // abuse tall flag
		}
		// add trim at top of door
		unsigned const num_floors(calc_num_floors(*d, window_vspacing, floor_thickness));
		float z(d->z1() + floor_to_ceil_height);
		trim.d[!d->dim][0] = d->d[!d->dim][0] + door_trim_width;
		trim.d[!d->dim][1] = d->d[!d->dim][1] - door_trim_width;

		for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
			trim.z1() = z - trim_thickness;
			trim.z2() = z; // ceil height
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, d->dim, 0, (flags | RO_FLAG_ADJ_TOP), 1.0, SHAPE_SHORT);
		}
	} // for d
	for (auto d = doors.begin(); d != doors.end(); ++d) { // exterior doors
		if (d->type == tquad_with_ix_t::TYPE_RDOOR) { // roof access door
			continue; // this requires a completely different approach to trim and has not yet been implemented
		}
		bool const garage_door(d->type == tquad_with_ix_t::TYPE_GDOOR);
		if (garage_door) continue; // no trim on garage door, needs different logic
		cube_t door(d->get_bcube());
		bool const dim(door.dy() < door.dx());
		//door.expand_in_dim(!dim, -0.1*door_trim_width); // shrink slightly so that the edge of the wall is contained in the trim
		cube_t trim(door);
		trim.expand_in_dim(dim, door_trim_exp);
		bool dir(0);
		unsigned ext_flags(flags);
		colorRGBA const &trim_color(garage_door ? WHITE : door_color); // garage doors are always white

		for (auto i = parts.begin(); i != get_real_parts_end(); ++i) {
			if (!i->intersects_no_adj(trim)) continue;
			trim.intersect_with_cube_xy(*i); // clip to containing part
			dir = (i->get_center_dim(dim) < trim.get_center_dim(dim));
			trim.d[dim][dir] -= (dir ? -1.0 : 1.0)*0.025*window_vspacing; // move to to same offset for door
			ext_flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
			break;
		}
		for (unsigned side = 0; side < 2; ++side) { // left/right of door
			trim.d[!dim][0] = door.d[!dim][side] - (side ? door_trim_width : 0.0);
			trim.d[!dim][1] = door.d[!dim][side] + (side ? 0.0 : door_trim_width);
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, side, (ext_flags | RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP), 1.0, SHAPE_TALL, trim_color); // abuse tall flag
		}
		// add trim at bottom of door for threshold
		trim.d[!dim][0] = door.d[!dim][0];
		trim.d[!dim][1] = door.d[!dim][1];
		trim.z1() = door.z1() + fc_thick; // floor height
		trim.z2() = trim.z1() + 2.0*trim_thickness;
		objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, (ext_flags | RO_FLAG_ADJ_BOT), 1.0, SHAPE_SHORT, trim_color);

		if (d->type == tquad_with_ix_t::TYPE_HDOOR || d->type == tquad_with_ix_t::TYPE_BDOOR) { // add trim at top of door, houses and office buildings
			trim.z1() = door.z2() - 0.03*door.dz(); // see logic in clip_door_to_interior()
			trim.z2() = door.z2(); // top of door texture
		}
		if (d->type == tquad_with_ix_t::TYPE_BDOOR) { // different logic for building doors
			ext_flags = flags; // unlike hdoors, need to draw the back face to hide the gap betweeen ceiling and floor above
			trim.d[dim][dir] += (dir ? -1.0 : 1.0)*0.005*window_vspacing; // minor shift back toward building to prevent z-fighting
		}
		objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_SHORT, trim_color); // top of door
	} // for d
	for (unsigned dim = 0; dim < 2; ++dim) { // add horizontal strips along each wall at each floor/ceiling
		for (auto w = interior->walls[dim].begin(); w != interior->walls[dim].end(); ++w) {
			cube_t trim(*w);
			trim.expand_in_dim(dim, trim_thickness);

			if (!is_house && !pri_hall.is_all_zeros()) { // handle outside corners of office building hallway intersections
				for (auto W = interior->walls[!dim].begin(); W != interior->walls[!dim].end(); ++W) { // check walls in other dim for an outside corner
					for (unsigned d = 0; d < 2; ++d) {
						if (W->z1() > w->z2() || W->z2() < w->z1()) continue; // no z overlap, wrong stack
						if (W->d[!dim][0] > w->d[!dim][d] || W->d[!dim][1] < w->d[!dim][d]) continue; // not adjacent/overlapping
						if (W->d[ dim][0] < w->d[ dim][0] && W->d[ dim][1] > w->d[ dim][1]) continue; // skip T junctions
						trim.d[!dim][d] = W->d[!dim][d] + (d ? 1.0 : -1.0)*trim_thickness; // expand to cover gap at outside corners of hallway walls
					}
				} // for W
			}
			unsigned const num_floors(calc_num_floors(*w, window_vspacing, floor_thickness));
			float z(w->z1());

			for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
				trim.z1() = z; // floor height
				trim.z2() = z + trim_height;
				objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, 0, flags, 1.0, SHAPE_CUBE); // floor trim
				if (!has_ceil_trim) continue;
				trim.z2() = z + floor_to_ceil_height; // ceil height
				trim.z1() = trim.z2() - trim_height;

				for (unsigned dir = 0; dir < 2; ++dir) { // for each side of wall
					cube_t ceil_trim(trim);
					ceil_trim.d[dim][!dir] = w->d[dim][dir];
					objs.emplace_back(ceil_trim, TYPE_WALL_TRIM, 0, dim, dir, flags, 1.0, SHAPE_ANGLED); // ceiling trim
				}
			} // for f
		} // for w
	} // for d
	for (auto i = parts.begin(); i != get_real_parts_end(); ++i) { // add trim for exterior walls
		unsigned const num_floors(calc_num_floors(*i, window_vspacing, floor_thickness));

		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				cube_t trim(*i);
				trim.d[dim][!dir] = i->d[dim][dir] + (dir ? -1.0 : 1.0)*trim_thickness;
				unsigned const ext_flags(flags | (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));
				float z(i->z1() + fc_thick);

				for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
					trim.z1() = z; // floor height
					trim.z2() = z + trim_height;
					trim_cubes.clear();
					trim_cubes.push_back(trim); // start with entire length
					float const ceil_trim_z2(z + floor_to_ceil_height), ceil_trim_z1(ceil_trim_z2 - trim_height); // ceil height

					for (auto j = parts.begin(); j != get_real_parts_end(); ++j) { // clip against other parts
						if (j == i) continue; // skip self
						cube_t clip_cube(*j);
						clip_cube.expand_in_dim(dim, expand_val); // expand to clip trim on the other side of the split wall
						subtract_cube_from_cubes(clip_cube, trim_cubes); // subtract this part from current trim cubes by clipping in XY
					}
					if (f == 0) { // first floor, cut out areas for exterior doors
						for (auto d = doors.begin(); d != doors.end(); ++d) {
							cube_t door(d->get_bcube());
							bool const door_dim(door.dy() < door.dx());
							if (door_dim != bool(dim)) continue;
							door.expand_in_dim(door_dim, (expand_val + wall_thickness)); // expand to nonzero area; use a larger expand to account for distance door is offset away from ext wall
							subtract_cube_from_cubes(door, trim_cubes); // subtract this door from current trim cubes by clipping in XY
						}
					}
					for (auto c = trim_cubes.begin(); c != trim_cubes.end(); ++c) {
						objs.emplace_back(*c, TYPE_WALL_TRIM, 0, dim, 0, ext_flags, 1.0, SHAPE_CUBE); // floor trim
						if (!has_ceil_trim) continue;
						if (is_house && c != trim_cubes.begin()) continue;
						cube_t ceil_trim(is_house ? trim : *c); // houses have shorter doors and ceiling trim extends above the door, so draw full range
						ceil_trim.z2() = ceil_trim_z2;
						ceil_trim.z1() = ceil_trim_z1;
						objs.emplace_back(ceil_trim, TYPE_WALL_TRIM, 0, dim, !dir, flags, 1.0, SHAPE_ANGLED); // ceiling trim
					}
				} // for f
			} // for dir
		} // for dim
	} // for i
	// add window trim
	if (is_rotated()) return; // not yet working for rotated buildings
	float const border_mult(0.94); // account for the frame part of the window texture, which is included in the interior cutout of the window
	float const window_h_border(border_mult*get_window_h_border()), window_v_border(border_mult*get_window_v_border()); // (0, 1) range
	// Note: depth must be small to avoid object intersections; this applies to the windowsill as well
	float const window_trim_width(0.75*wall_thickness), window_trim_depth(0.1*wall_thickness), windowsill_depth(0.1*wall_thickness);
	float const window_offset(0.01*window_vspacing); // must match building_draw_t::add_section()
	static vect_vnctcc_t wall_quad_verts;
	wall_quad_verts.clear();
	get_all_drawn_window_verts_as_quads(wall_quad_verts);
	assert((wall_quad_verts.size() & 3) == 0); // must be a multiple of 4

	for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
		auto const &v0(wall_quad_verts[i]);
		cube_t c(v0.v);
		float tx1(v0.t[0]), tx2(tx1), tz1(v0.t[1]), tz2(tz1); // tex coord ranges (xy, z); should generally be whole integers

		for (unsigned j = 1; j < 4; ++j) {
			auto const &vj(wall_quad_verts[i + j]);
			c.union_with_pt(vj.v);
			min_eq(tx1, vj.t[0]);
			max_eq(tx2, vj.t[0]);
			min_eq(tz1, vj.t[1]);
			max_eq(tz2, vj.t[1]);
		}
		if (tx1 == tx2 || tz1 == tz2) continue; // wall is too small to contain a window
		assert(tx2 - tx1 < 1000.0f && tz2 - tz1 < 1000.0f); // sanity check - less than 1000 windows in each dim
		assert(c.dz() > 0.0);
		bool const dim(c.dy() < c.dx()), dir(v0.get_norm()[dim] > 0.0);
		assert(c.get_sz_dim(dim) == 0.0); // must be zero size in one dim (X or Y oriented); could also use the vertex normal
		float const d_tx_inv(1.0f/(tx2 - tx1)), d_tz_inv(1.0f/(tz2 - tz1));
		float const window_width(c.get_sz_dim(!dim)*d_tx_inv), window_height(c.dz()*d_tz_inv); // window_height should be equal to window_vspacing
		float const border_xy(window_width*window_h_border), border_z(window_height*window_v_border);
		cube_t window(c); // copy dim <dim>
		window.translate_dim((dir ? -1.0 : 1.0)*window_offset, dim);
		window.d[dim][!dir] += (dir ? -1.0 : 1.0)*window_trim_depth; // add thickness on interior of building
		unsigned ext_flags(flags | (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));

		for (float z = tz1; z < tz2; z += 1.0) { // each floor
			float const bot_edge(c.z1() + (z - tz1)*window_height);
			window.z1() = bot_edge + border_z;
			window.z2() = bot_edge + window_height - border_z;

			for (float xy = tx1; xy < tx2; xy += 1.0) { // windows along each wall
				float const low_edge(c.d[!dim][0] + (xy - tx1)*window_width);
				window.d[!dim][0] = low_edge + border_xy;
				window.d[!dim][1] = low_edge + window_width - border_xy;
				cube_t top(window), bot(window), side(window);
				top.z1()  = window.z2();
				top.z2() += window_trim_width;
				bot.z2()  = window.z1();
				bot.z1() -= window_trim_width;
				bot.d[dim][!dir] += (dir ? -1.0f : 1.0f)*(windowsill_depth - window_trim_depth); // shift out further for windowsill
				top.expand_in_dim(!dim, window_trim_width);
				bot.expand_in_dim(!dim, window_trim_width);
				objs.emplace_back(top, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL);
				objs.emplace_back(bot, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL);

				for (unsigned s = 0; s < 2; ++s) { // left/right sides
					side.d[!dim][ s] = window.d[!dim][s] - (s ? -1.0 : 1.0)*window_trim_width;
					side.d[!dim][!s] = window.d[!dim][s];
					objs.emplace_back(side, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL);
				}
			} // for xy
		} // for z
	} // for i
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), half_thick(0.5*floor_thickness);
	vector<room_object_t> &objs(interior->room_geom->objs);
	interior->room_geom->stairs_start = objs.size();
	interior->room_geom->has_elevators = (!interior->elevators.empty());
	colorRGBA const railing_colors[3] = {GOLD, LT_GRAY, BLACK};
	colorRGBA const railing_color(railing_colors[rgen.rand()%3]); // set per-building

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator) continue; // for elevator, not stairs
		unsigned const num_stairs((i->shape == SHAPE_U) ? NUM_STAIRS_PER_FLOOR_U : NUM_STAIRS_PER_FLOOR);
		float const stair_dz(window_vspacing/(num_stairs+1)), stair_height(stair_dz + floor_thickness);
		bool const dim(i->dim), dir(i->dir), has_side_walls(i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES || i->shape == SHAPE_U);
		bool const side(dir); // for U-shaped stairs; for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
		// Note: stairs always start at floor_thickness above the landing z1, ignoring landing z2/height
		float const tot_len(i->get_sz_dim(dim)), floor_z(i->z1() + floor_thickness - window_vspacing), step_len_pos(tot_len/num_stairs);
		float step_len((dir ? 1.0 : -1.0)*step_len_pos), wall_hw(0.15*step_len_pos), z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i);

		if (i->shape != SHAPE_U) { // straight stairs
			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				stair.z1() = max(bcube.z1(), (z + 0.5f*half_thick)); // don't go below the floor
				stair.z2() = z + stair_height;
				assert(stair.z1() < stair.z2());
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // Note: room_id=0, not tracked, unused
			}
		}
		else { // U-shaped stairs
			stair.d[!dim][side] = i->get_center_dim(!dim);
			step_len *= 2.0;

			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				if (n == num_stairs/2) { // reverse direction and switch to other side
					step_len *= -1.0;
					stair.d[!dim][ side] = i->d[!dim][side];
					stair.d[!dim][!side] = i->get_center_dim(!dim);
				}
				assert(!(num_stairs & 1)); // require num_stairs to be an even number
				bool const is_rev(n >= num_stairs/2), stairs_dir(dir^is_rev);
				stair.d[dim][!stairs_dir] = pos; stair.d[dim][stairs_dir] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, side^is_rev); // Note: room_id=0, not tracked, unused
				objs.back().shape = SHAPE_STAIRS_U;
			} // for n
		}
		// add walls and railings
		bool const extend_walls_up(i->is_at_top && !i->roof_access); // space above is open, add a wall so that people can't fall down the stairs
		float const railing_z2(i->z2() + (i->roof_access ? 0.025*i->dz() : 0.0)); // capture z2 before we change it; move roof access railing up a bit to offset the shrink resize
		cube_t wall(*i);
		if (extend_walls_up) {wall.z2() += window_vspacing - floor_thickness;}
		else {wall.z2() -= 0.5*floor_thickness;} // prevent z-fighting on top floor
		wall.z1() = max(bcube.z1()+half_thick, floor_z-half_thick); // full height
		set_wall_width(wall, i->d[dim][dir], wall_hw, dim);

		if ((i->shape == SHAPE_WALLED && !(i->against_wall[0] || i->against_wall[1])) || i->shape == SHAPE_U) {
			objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir); // add wall at back/end of stairs
		}
		else if ((i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES) && extend_walls_up) { // add upper section only
			cube_t wall_upper(wall);
			set_wall_width(wall_upper, (i->d[dim][!dir] + (dir ? 1.0 : -1.0)*wall_hw), wall_hw, dim); // move to the other side
			wall_upper.z1() = railing_z2;
			objs.emplace_back(wall_upper, TYPE_STAIR_WALL, 0, dim, dir); // add wall at back/end of stairs
		}
		wall.d[dim][!dir] = i->d[dim][!dir];

		for (unsigned d = 0; d < 2; ++d) { // sides of stairs
			set_wall_width(wall, i->d[!dim][d], wall_hw, !dim);
			wall.expand_in_dim(dim, 0.01*wall_hw); // just enough to avoid z-fighting with stairs
			bool const add_wall(has_side_walls && !i->against_wall[d]); // don't add a wall if the stairs are already against a wall
			if (add_wall) {objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir);} // add walls around stairs for this floor

			if (i->has_railing) { // add railings
				bool railing_dir(dir);
				cube_t railing(wall);
				uint16_t flags(add_wall ? RO_FLAG_NOCOLL : 0);
				railing.z2() = railing_z2;

				if (add_wall || i->roof_access) {
					railing.translate_dim((d ? -1.0 : 1.0)*2.0*wall_hw, !dim); // shift railing inside of walls
					railing.expand_in_dim(dim, -(i->roof_access ? 2.0 : 1.0)*wall_hw); // shrink slightly to avoid clipping through an end wall
				}
				if (i->shape == SHAPE_U) { // adjust railing height/angle to match stairs
					float const z_split(railing.get_center_dim(2));
					if (bool(d) == side) {railing.z1() = z_split; flags |= RO_FLAG_ADJ_HI; railing_dir ^= 1;}
					else                 {railing.z2() = z_split; flags |= RO_FLAG_ADJ_LO;}
				}
				objs.emplace_back(railing, TYPE_RAILING, 0, dim, railing_dir, flags, 1.0, SHAPE_CUBE, railing_color); // collision works like a cube
			}
		} // for d
	} // for i
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		unsigned const elevator_id(i - interior->elevators.begin()); // used for room_object_t::room_id
		cube_t elevator_car(*i);
		elevator_car.z1() += 0.05*get_floor_thickness(); // to prevent z-fighting when looking at the building from the bottom
		elevator_car.z2() = elevator_car.z1() + window_vspacing; // currently at the bottom floor
		objs.emplace_back(elevator_car, TYPE_ELEVATOR, elevator_id, i->dim, i->dir, (i->open ? RO_FLAG_OPEN : 0));
	}
}

void building_t::add_sign_by_door(tquad_with_ix_t const &door, bool outside, std::string const &text, colorRGBA const &color, bool emissive) {
	cube_t const door_bcube(door.get_bcube());
	bool const dim(door_bcube.dy() < door_bcube.dx());
	float const width(door_bcube.get_sz_dim(!dim)), height(door_bcube.dz());
	cube_t c(door_bcube);

	if (outside) { // outside, place above the door
		c.z2() = door_bcube.z2() + 0.1*height;
	}
	else { // inside, place hanging near the top of the door
		c.z2() = door_bcube.z1() + get_window_vspace() - 0.5*get_floor_thickness(); // right against the ceiling
	}
	c.z1() = c.z2() - 0.05*height;
	float const sign_width(0.8*text.size()*c.dz()), shrink(0.5f*(width - sign_width));
	c.expand_in_dim(!dim, -shrink);
	vector<room_object_t> &objs(interior->room_geom->objs);

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find part containing this door so that we can get the correct dir
		if (p->z1() != bcube.z1()) continue; // not ground floor
		if (p->d[!dim][1] < door_bcube.d[!dim][1] || p->d[!dim][0] > door_bcube.d[!dim][0]) {continue;} // not contained in this dim
		bool dir(0);
		if      (fabs(p->d[dim][0] - door_bcube.d[dim][0]) < 0.1*width) {dir = 0;}
		else if (fabs(p->d[dim][1] - door_bcube.d[dim][1]) < 0.1*width) {dir = 1;}
		else {continue;} // wrong part
		if (!outside) {dir ^= 1; c.translate_dim((dir ? 1.0 : -1.0)*0.1*height, dim);} // move inside the building
		c.d[dim][dir] += (dir ? 1.0 : -1.0)*0.01*height;

		if (outside) {
			bool skip(0);
			for (auto p2 = get_real_parts_end_inc_sec(); p2 != parts.end(); ++p2) {skip |= p2->intersects(c);}
			if (skip) return; // sign intersects porch roof, skip this building
		}
		unsigned flags(RO_FLAG_LIT | RO_FLAG_NOCOLL | (emissive ? RO_FLAG_EMISSIVE : 0) | (outside ? 0 : RO_FLAG_HANGING));
		objs.emplace_back(c, TYPE_SIGN, 0, dim, dir, flags, 1.0, SHAPE_CUBE, color); // always lit; room_id is not valid
		objs.back().obj_id = register_sign_text(text);
		return; // done
	} // for p
	cout << TXT(bcube.str()) << TXT(door_bcube.str()) << TXT(is_house) << endl; // debug printout
	//assert(0); // never gets here (too strong?)
}

void building_t::add_exterior_door_signs(rand_gen_t &rgen) {
	if (is_house) { // maybe add welcome sign
		if (rgen.rand() % 5) return; // only 20% of houses have a welcome sign
		assert(!doors.empty());
		add_sign_by_door(doors.front(), 1, "Welcome", DK_BROWN, 0); // front door only, outside
	}
	else { // add exit signs
		if (pri_hall.is_all_zeros() && rgen.rand_bool()) return; // place exit signs on buildings with primary hallways and 50% of other buildings
		colorRGBA const exit_color(rgen.rand_bool() ? RED : GREEN);
		
		for (auto d = doors.begin(); d != doors.end(); ++d) {
			if (has_courtyard && (d+1) == doors.end()) break; // courtyard door is not an exit
			if (d->type == tquad_with_ix_t::TYPE_BDOOR) {add_sign_by_door(*d, 0, "Exit", exit_color, 1);} // inside, emissive
		}
	}
}

void building_t::draw_room_geom(shader_t &s, occlusion_checker_t &oc, vector3d const &xlate, unsigned building_ix,
	bool shadow_only, bool reflection_pass, bool inc_small, bool player_in_building)
{
	if (!interior || !interior->room_geom) return;
	if (ENABLE_MIRROR_REFLECTIONS && !shadow_only && !reflection_pass && player_in_building) {find_mirror_needing_reflection(xlate);}
	interior->room_geom->draw(s, *this, oc, xlate, get_material().wall_tex, building_ix, shadow_only, reflection_pass, inc_small, player_in_building);
}
void building_t::gen_and_draw_room_geom(shader_t &s, occlusion_checker_t &oc, vector3d const &xlate, vect_cube_t &ped_bcubes,
	unsigned building_ix, int ped_ix, bool shadow_only, bool reflection_pass, bool inc_small, bool player_in_building)
{
	if (!interior) return;
	if (!global_building_params.enable_rotated_room_geom && is_rotated()) return; // rotated buildings: need to fix texture coords, room object collision detection, mirrors, etc.

	if (!has_room_geom()) {
		rand_gen_t rgen;
		rgen.set_state(building_ix, parts.size()); // set to something canonical per building
		ped_bcubes.clear();
		if (ped_ix >= 0) {get_ped_bcubes_for_building(ped_ix, building_ix, ped_bcubes);}
		gen_room_details(rgen, ped_bcubes, building_ix); // generate so that we can draw it
		assert(has_room_geom());
	}
	draw_room_geom(s, oc, xlate, building_ix, shadow_only, reflection_pass, inc_small, player_in_building);
}

void building_t::clear_room_geom() {
	if (!has_room_geom()) return;
	interior->room_geom->clear(); // free VBO data before deleting the room_geom object
	interior->room_geom.reset();
}

room_t::room_t(cube_t const &c, unsigned p, unsigned nl, bool is_hallway_, bool is_office_, bool is_sec_bldg_) :
	cube_t(c), has_stairs(0), has_elevator(0), no_geom(is_hallway_), is_hallway(is_hallway_), is_office(is_office_), // no geom in hallways
	is_sec_bldg(is_sec_bldg_), interior(0), has_bathroom(0), ext_sides(0), part_id(p), num_lights(nl), lit_by_floor(0), light_intensity(0.0)
{
	if      (is_sec_bldg) {rtype = RTYPE_GARAGE;} // or RTYPE_SHED - will be set later
	else if (is_hallway)  {rtype = RTYPE_HALL;}
	else if (is_office)   {rtype = RTYPE_OFFICE;}
	else if (has_stairs)  {rtype = RTYPE_STAIRS;}
	else                  {rtype = RTYPE_NOTSET;}
}
void room_t::assign_to(room_type rt, unsigned floor) {
	if (rtype == RTYPE_NOTSET || floor == 0) {rtype = rt;} // rtype is not per floor: assign if not yet assigned or this is the first floor (first floor has priority)
	if (rt == RTYPE_BATH) {has_bathroom = 1;} // tag as has_bathroom even if bathroom is not on the ground floor
}

