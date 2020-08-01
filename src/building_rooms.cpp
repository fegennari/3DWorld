// 3D World - Building Interior Room Geometry Placement
// by Frank Gennari 4/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#pragma warning(disable : 26812) // prefer enum class over enum

extern object_model_loader_t building_obj_model_loader;


bool building_t::overlaps_other_room_obj(cube_t const &c, unsigned objs_start) const {
	assert(interior && interior->room_geom);
	vector<room_object_t> &objs(interior->room_geom->objs);
	assert(objs_start <= objs.size());

	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (i->intersects(c)) return 1;
	}
	return 0;
}

bool building_t::is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, bool inc_open_doors, float room_pad) const {
	cube_t place_area(room);
	if (room_pad != 0.0f) {place_area.expand_by_xy(-room_pad);} // shrink by dmin
	if (!place_area.contains_cube_xy(c)) return 0; // not contained in interior part of the room
	if (is_cube_close_to_doorway(c, 0.0, inc_open_doors)) return 0; // too close to a doorway
	if (interior && interior->is_blocked_by_stairs_or_elevator(c)) return 0; // faster to check only one per stairwell, but then we need to store another vector?
	if (has_bcube_int(c, blockers)) return 0; // Note: ignores dmin
	return 1;
}

bool building_t::add_chair(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id, point const &place_pos,
	colorRGBA const &chair_color, bool dim, bool dir, float tot_light_amt, bool is_lit)
{
	float const window_vspacing(get_window_vspace()), room_pad(4.0f*get_wall_thickness()), chair_sz(0.1*window_vspacing); // half size
	point chair_pos(place_pos); // same starting center and z1
	chair_pos[dim] += (dir ? -1.0f : 1.0f)*rgen.rand_uniform(-0.5, 1.2)*chair_sz;
	cube_t chair(chair_pos, chair_pos);
	chair.z2() += 0.4*get_window_vspace(); // chair height
	chair.expand_by_xy(chair_sz);
	if (!is_valid_placement_for_room(chair, room, blockers, 0, room_pad)) return 0; // check proximity to doors
	vector<room_object_t> &objs(interior->room_geom->objs);
	objs.emplace_back(chair, TYPE_CHAIR, room_id, dim, dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt);
	objs.back().color = chair_color;
	return 1;
}

// Note: must be first placed objects; returns the number of total objects added (table + optional chairs)
unsigned building_t::add_table_and_chairs(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id,
	point const &place_pos, colorRGBA const &chair_color, float rand_place_off, float tot_light_amt, bool is_lit)
{
	float const window_vspacing(get_window_vspace()), room_pad(4.0f*get_wall_thickness());
	uint8_t const obj_flags(is_lit ? RO_FLAG_LIT : 0);
	vector3d const room_sz(room.get_size());
	assert(interior && interior->room_geom);
	vector<room_object_t> &objs(interior->room_geom->objs);
	point table_pos(place_pos);
	vector3d table_sz;
	for (unsigned d = 0; d < 2; ++d) {table_sz [d]  = 0.18*window_vspacing*(1.0 + rgen.rand_float());} // half size relative to window_vspacing
	for (unsigned d = 0; d < 2; ++d) {table_pos[d] += rand_place_off*room_sz[d]*rgen.rand_uniform(-1.0, 1.0);} // near the center of the room
	point llc(table_pos - table_sz), urc(table_pos + table_sz);
	llc.z = table_pos.z; // bottom
	urc.z = table_pos.z + 0.2*window_vspacing;
	cube_t table(llc, urc);
	if (!is_valid_placement_for_room(table, room, blockers, 0, room_pad)) return 0; // check proximity to doors and collision with blockers
	objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, obj_flags, tot_light_amt);
	unsigned num_added(1); // start with the table

	// place some chairs around the table
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (rgen.rand_bool()) continue; // 50% of the time
			point chair_pos(table_pos); // same starting center and z1
			chair_pos[dim] += (dir ? -1.0f : 1.0f)*table_sz[dim];
			num_added += add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, is_lit);
		}
	}
	return num_added;
}
void building_t::shorten_chairs_in_region(cube_t const &region, unsigned objs_start) {
	for (auto i = interior->room_geom->objs.begin() + objs_start; i != interior->room_geom->objs.end(); ++i) {
		if (i->type != TYPE_CHAIR || !i->intersects(region)) continue;
		i->z2() -= 0.25*i->dz();
		i->type = TYPE_SM_CHAIR;
	}
}

void building_t::add_trashcan_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start, bool check_last_obj) {
	unsigned const NUM_COLORS = 6;
	colorRGBA const colors[NUM_COLORS] = {BLUE, DK_GRAY, LT_GRAY, GRAY, BLUE, WHITE};
	int const rr(rgen.rand()%3), rar(rgen.rand()%3); // three sizes/ARs
	float const radius(0.02f*(3 + rr)*get_window_vspace()), height(0.54f*(3 + rar)*radius);
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
	vect_cube_t doorways;

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
		cube_t c(center, center);
		c.expand_by_xy(radius);
		c.z2() += height;
		if (is_cube_close_to_doorway(c, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(c) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_TCAN, room_id, dim, dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt, (cylin ? room_obj_shape::SHAPE_CYLIN : room_obj_shape::SHAPE_CUBE));
		objs.back().color = colors[rgen.rand()%NUM_COLORS];
		return; // done
	} // for n
}

// Note: no blockers for people
bool building_t::add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return 0; // room is too small
	float const width(0.4*vspace*rgen.rand_uniform(1.0, 1.2)), depth(0.12*vspace*rgen.rand_uniform(1.0, 1.2)), height(0.7*vspace*rgen.rand_uniform(1.0, 1.2));
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
		tc.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.2*vspace; // increase space to add clearance
		if (is_cube_close_to_doorway(tc, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(tc) || overlaps_other_room_obj(tc, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_BCASE, room_id, dim, !dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt); // Note: dir faces into the room, not the wall
		objs.back().obj_id = (uint16_t)objs.size();
		return 1; // done/success
	} // for n
	return 0; // not placed
}

// Note: must be first placed object
bool building_t::add_desk_to_room(rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers,
	colorRGBA const &chair_color, float zval, unsigned room_id, float tot_light_amt, bool is_lit)
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
		c.d[dim][ dir] = room_bounds.d[dim][dir] + rgen.rand_uniform(0.1, 1.0)*(dir ? -1.0 : 1.0)*get_wall_thickness(); // almost against this wall
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		c.d[!dim][0] = pos - 0.5*width;
		c.d[!dim][1] = pos + 0.5*width;
		if (num_placed > 0 && c.intersects(placed_desk)) continue; // intersects previously placed desk
		if (!is_valid_placement_for_room(c, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
		bool const is_tall(!room.is_office && rgen.rand_float() < 0.5 && classify_room_wall(room, zval, dim, dir, 0) != ROOM_WALL_EXT); // make short if against an exterior wall or in an office
		objs.emplace_back(c, TYPE_DESK, room_id, dim, !dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt, (is_tall ? SHAPE_TALL : SHAPE_CUBE));
		objs.back().obj_id = (uint16_t)objs.size();
		cube_t bc(c);

		if (rgen.rand_float() > 0.05) { // 5% chance of no chair
			point chair_pos;
			chair_pos.z = zval;
			chair_pos[dim]  = c.d[dim][!dir];
			chair_pos[!dim] = pos + rgen.rand_uniform(-0.1, 0.1)*width; // slightly misaligned
			
			if (add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, is_lit)) {
				cube_t const &chair(objs.back());
				if (num_placed > 0 && chair.intersects(placed_desk)) {objs.pop_back();} // intersects previously placed desk, remove it
				else {bc.union_with_cube(chair);} // include the chair
			}
		}
		++num_placed;
		if (room.is_office && num_placed == 1 && rgen.rand_float() < 0.5) {placed_desk = bc; continue;} // allow two desks in one office
		break; // done/success
	} // for n
	return (num_placed > 0);
}

bool building_t::create_office_cubicles(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
	if (!room.is_office) return 0;
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const floor_spacing(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 2.5*floor_spacing || max(room_bounds.dx(), room_bounds.dy()) < 3.5*floor_spacing) return 0; // not large enough
	// TODO - WRITE
	return 0;
}

bool building_t::can_be_bedroom_or_bathroom(room_t const &room, bool on_first_floor) const { // check room type and existence of exterior door
	if (room.has_stairs || room.has_elevator || room.is_hallway || room.is_office) return 0; // no bed/bath in these cases
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

// Note: must be first placed object
bool building_t::add_bed_to_room(rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, float zval, unsigned room_id, float tot_light_amt, bool is_lit) {
	unsigned const NUM_COLORS = 8;
	colorRGBA const colors[NUM_COLORS] = {WHITE, WHITE, WHITE, LT_BLUE, LT_BLUE, PINK, PINK, LT_GREEN}; // color of the sheets
	cube_t room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace()), wall_thick(get_wall_thickness());
	bool const dim(room_bounds.dx() < room_bounds.dy()); // longer dim
	vector3d expand, bed_sz;
	expand[ dim] = -wall_thick; // small amount of space
	expand[!dim] = -0.3*vspace; // leave at least some space between the bed and the wall
	room_bounds.expand_by_xy(expand);
	if (room_bounds.get_sz_dim(dim) < 1.3*vspace || room_bounds.get_sz_dim(!dim) < 0.8*vspace) return 0; // room is too small to fit a bed
	if (room_bounds.get_sz_dim(dim) > 4.0*vspace || room_bounds.get_sz_dim(!dim) > 2.5*vspace) return 0; // room is too large to be a bedroom
	bool const first_head_dir(rgen.rand_bool()), first_wall_dir(rgen.rand_bool());
	vector<room_object_t> &objs(interior->room_geom->objs);
	cube_t c;
	c.z1() = zval;

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a bed
		bed_sz[ dim] = 0.8*vspace*rgen.rand_uniform(1.0, 1.2); // length
		bed_sz[!dim] = 0.5*vspace*rgen.rand_uniform(1.0, 1.5); // width
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
		objs.emplace_back(c, TYPE_BED, room_id, dim, dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt);
		room_object_t &bed(objs.back());
		bed.obj_id = (uint16_t)objs.size();
		// use white color if a texture is assigned that's not close to white
		int const sheet_tid(bed.get_sheet_tid());
		if (sheet_tid < 0 || sheet_tid == WHITE_TEX || texture_color(sheet_tid).get_luminance() > 0.5) {bed.color = colors[rgen.rand()%NUM_COLORS];}
		return 1; // done/success
	} // for n
	return 0;
}

bool building_t::place_obj_along_wall(room_object type, float height, vector3d const &sz_scale, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
	bool is_lit, cube_t const &place_area, unsigned objs_start, float front_clearance, unsigned pref_orient, bool pref_centered, colorRGBA const &color)
{
	float const hwidth(0.5*height*sz_scale.y/sz_scale.z), depth(height*sz_scale.x/sz_scale.z), min_space(2.8*hwidth);
	vector3d const place_area_sz(place_area.get_size());
	if (max(place_area_sz.x, place_area_sz.y) <= min_space) return 0; // can't fit in either dim
	unsigned const force_dim((place_area_sz.x <= min_space) ? 0 : ((place_area_sz.y <= min_space) ? 1 : 2)); // *other* dim; 2=neither
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
		cube_t c2(c); // used for collision tests
		c2.d[dim][!dir] += (dir ? -1.0 : 1.0)*depth*front_clearance;
		if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(c2)) continue; // bad placement
		objs.emplace_back(c, type, room_id, dim, !dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt, room_obj_shape::SHAPE_CUBE, color);
		objs.back().obj_id = (uint16_t)objs.size();
		objs.emplace_back(c2, TYPE_BLOCKER, room_id); // add blocker cube to ensure no other object overlaps this space
		return 1; // done
	} // for n
	return 0; // failed
}
bool building_t::place_model_along_wall(unsigned model_id, room_object type, float height, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
	bool is_lit, cube_t const &place_area, unsigned objs_start, float front_clearance, unsigned pref_orient, bool pref_centered, colorRGBA const &color)
{
	if (!building_obj_model_loader.is_model_valid(model_id)) return 0; // don't have a model of this type
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
	return place_obj_along_wall(type, height*get_window_vspace(), sz, rgen, zval, room_id, tot_light_amt,
		is_lit, place_area, objs_start, front_clearance, pref_orient, pref_centered, color);
}

bool building_t::add_bathroom_objs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.5*wall_thickness);
	if (min(place_area.dx(), place_area.dy()) < 0.7*floor_spacing) return 0; // room is too small (should be rare)
	bool const have_toilet(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOILET)), have_sink(building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK));

	if (have_toilet && have_sink && room.is_office && min(place_area.dx(), place_area.dy()) > 1.8*floor_spacing && max(place_area.dx(), place_area.dy()) > 2.6*floor_spacing) {
		if (divide_bathroom_into_stalls(rgen, room, zval, room_id, tot_light_amt, is_lit, objs_start)) return 1; // large enough, try to divide into bathroom stalls
	}
	vector<room_object_t> &objs(interior->room_geom->objs);
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
				if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, 0.0, 1)) continue; // bad placement
				objs.emplace_back(c, TYPE_TOILET, room_id, dim, !dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt);
				objs.emplace_back(c2, TYPE_BLOCKER, room_id); // add blocker cube to ensure no other object overlaps this space
				placed_obj = placed_toilet = 1; // done
			} // for d
		} // for n
		if (!placed_toilet) { // if the toilet can't be placed in a corner, allow it to be placed anywhere; needed for small offices
			placed_obj |= place_model_along_wall(OBJ_MODEL_TOILET, TYPE_TOILET, 0.35, rgen, zval, room_id, tot_light_amt, is_lit, place_area, objs_start, 0.8);
		}
	}
	if (is_house) { // place a tub, but not in office buildings; placed before the sink because it's the largest and the most limited in valid locations
		cube_t place_area_tub(room_bounds);
		place_area_tub.expand_by(-0.05*wall_thickness); // just enough to prevent z-fighting
		placed_obj |= place_model_along_wall(OBJ_MODEL_TUB, TYPE_TUB, 0.2, rgen, zval, room_id, tot_light_amt, is_lit, place_area_tub, objs_start, 0.4);
	}
	placed_obj |= place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, 0.45, rgen, zval, room_id, tot_light_amt, is_lit, place_area, objs_start, 0.6);
	return placed_obj;
}

bool building_t::divide_bathroom_into_stalls(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	vector3d const tsz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
	vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SINK  )); // L, W, H
	float const theight(0.35*floor_spacing), twidth(theight*tsz.y/tsz.z), tlength(theight*tsz.x/tsz.z), stall_depth(2.2*tlength);
	float const sheight(0.45*floor_spacing), swidth(sheight*ssz.y/ssz.z), slength(sheight*ssz.x/ssz.z);
	float stall_width(2.0*twidth), sink_spacing(1.75*swidth);
	bool br_dim(room.dy() < room.dx()), sink_side(0), sink_side_set(0);
	cube_t place_area(room);
	place_area.expand_by(-0.5*wall_thickness);

	// determine men's room vs. women's room (however, they are currently the same because there is no urinal model)
	assert(room.part_id < parts.size());
	point const part_center(parts[room.part_id].get_cube_center()), room_center(room.get_cube_center());
	bool const mens_room((part_center.x < room_center.x) ^ (part_center.y < room_center.y));

	for (unsigned d = 0; d < 2 && !sink_side_set; ++d) {
		for (unsigned side = 0; side < 2 && !sink_side_set; ++side) {
			cube_t c(room);
			c.z1() = zval; // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
			c.z2() = zval + wall_thickness;
			c.d[!br_dim][!side] = c.d[!br_dim][side] + (side ? -1.0 : 1.0)*wall_thickness; // shrink to near zero area in this dim

			for (auto i = interior->doors.begin(); i != interior->doors.end(); ++i) {
				if ((i->dy() < i->dx()) == br_dim) continue; // door in wrong dim
				if (!is_cube_close_to_door(c, 0.0, 0, *i, 1)) continue; // check_zval=1
				sink_side = side; sink_side_set = 1;
				place_area.d[!br_dim][side] += (sink_side ? -1.0 : 1.0)*(i->get_sz_dim(br_dim) - 0.25*swidth); // add sink clearance for the door to close
				break; // sinks are on the side closest to the door
			}
		} // for side
		if (d == 0 && !sink_side_set) {br_dim ^= 1;} // door not found on long dim - R90 and try short dim
	} // for d
	assert(sink_side_set);
	float const room_len(place_area.get_sz_dim(!br_dim)), room_width(place_area.get_sz_dim(br_dim));
	float const sinks_len(0.4*room_len), stalls_len(room_len - sinks_len), req_depth(2.0f*max(stall_depth, slength));
	if (room_width < req_depth) return 0;
	unsigned const num_stalls(floor(stalls_len/stall_width)), num_sinks(floor(sinks_len/sink_spacing));
	//cout << TXT(two_rows) << TXT(num_stalls) << TXT(num_sinks) << endl;
	if (num_stalls < 2 || num_sinks < 2) return 0; // not enough space for 2 stalls and 2 sinks
	stall_width  = stalls_len/num_stalls; // reclaculate to fill the gaps
	sink_spacing = sinks_len/num_sinks;
	bool const two_rows(room_width > 1.5*req_depth), skip_stalls_side(two_rows ? 0 : (room_id & 1)); // put stalls on a side consistent across floors
	float const stall_step((sink_side ? 1.0 : -1.0)*stall_width), sink_step((sink_side ? -1.0 : 1.0)*sink_spacing);
	unsigned const flags(is_lit ? RO_FLAG_LIT : 0);
	vector<room_object_t> &objs(interior->room_geom->objs);

	for (unsigned dir = 0; dir < 2; ++dir) { // each side of the wall
		if (!two_rows && dir == (unsigned)skip_stalls_side) continue; // no stalls/sinks on this side
		// add stalls
		float const dir_sign(dir ? -1.0 : 1.0), wall_pos(place_area.d[br_dim][dir]), stall_from_wall(wall_pos + dir_sign*(0.5*tlength + wall_thickness));
		float stall_pos(place_area.d[!br_dim][!sink_side] + 0.5*stall_step);

		for (unsigned n = 0; n < num_stalls; ++n) {
			point center(stall_from_wall, stall_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t toilet(center, center);
			toilet.expand_in_dim( br_dim, 0.5*tlength);
			toilet.expand_in_dim(!br_dim, 0.5*twidth);
			toilet.z2() += theight;
			objs.emplace_back(toilet, TYPE_TOILET, room_id, br_dim, !dir, flags, tot_light_amt);
			cube_t stall(center, center);
			stall.z2() = stall.z1() + floor_spacing; // set stall height to room height
			stall.expand_in_dim(!br_dim, 0.5*stall_width);
			stall.d[br_dim][ dir] = wall_pos; // + wall_thickness?
			stall.d[br_dim][!dir] = wall_pos + dir_sign*stall_depth;
			objs.emplace_back(stall, TYPE_STALL, room_id, br_dim, dir, flags, tot_light_amt, SHAPE_CUBE, colorRGBA(0.75, 1.0, 0.9, 1.0)); // blue-green
			stall_pos += stall_step;
		} // for n
		// add sinks
		float sink_pos(place_area.d[!br_dim][sink_side] + 0.5*sink_step), sink_from_wall(wall_pos + dir_sign*(0.5*slength + wall_thickness));

		for (unsigned n = 0; n < num_sinks; ++n) {
			point center(sink_from_wall, sink_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t sink(center, center);
			sink.expand_in_dim( br_dim, 0.5*slength);
			sink.expand_in_dim(!br_dim, 0.5*swidth);
			sink.z2() += sheight;
			objs.emplace_back(sink, TYPE_SINK, room_id, br_dim, !dir, flags, tot_light_amt);
			sink_pos += sink_step;
		} // for n
	} // for dir
	return 1;
}

bool building_t::add_kitchen_objs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start, bool allow_adj_ext_door) {
	// Note: table and chairs have already been placed
	if (room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be kitchens
	if (!is_house && rgen.rand_bool()) return 0; // some office buildings have kitchens, allow it half the time
	// if it has an external door then reject the room half the time; most houses don't have a front door to the kitchen
	if (is_room_adjacent_to_ext_door(room, 1) && (!allow_adj_ext_door || rgen.rand_bool())) return 0; // front_door_only=1
	float const wall_thickness(get_wall_thickness());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall for appliances
	bool placed_obj(0);
	placed_obj |= place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, 0.72, rgen, zval, room_id, tot_light_amt, is_lit, place_area, objs_start, 1.0);
	if (is_house) {placed_obj |= place_model_along_wall(OBJ_MODEL_STOVE, TYPE_STOVE, 0.50, rgen, zval, room_id, tot_light_amt, is_lit, place_area, objs_start, 0.8);}
	return placed_obj;
}

bool building_t::add_livingroom_objs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
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
	
	if (place_model_along_wall(OBJ_MODEL_COUCH, TYPE_COUCH, 0.40, rgen, zval, room_id, tot_light_amt, is_lit, place_area, objs_start, 0.67, 4, 1, couch_color)) { // pref centered
		placed_couch   = 1;
		tv_pref_orient = (2*objs[couch_ix].dim + !objs[couch_ix].dir); // TV should be across from couch
	}
	tv_ix = objs.size();

	if (place_model_along_wall(OBJ_MODEL_TV, TYPE_TV, 0.45, rgen, zval, room_id, tot_light_amt, is_lit, place_area, objs_start, 4.0, tv_pref_orient, 1, BKGRAY)) { // pref centered
		placed_tv = 1;
		// add a small table to place the TV on so that it's off the floor and not blocked as much by tables and chairs
		room_object_t &tv(objs[tv_ix]);
		float const height(0.4*tv.dz());
		cube_t table(tv); // same XY bounds as the TV
		tv.translate_dim(height, 2); // move TV up
		table.z2() = tv.z1();
		objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt);
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

bool building_t::add_library_objs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
	if (room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be libraries
	unsigned num_added(0);

	for (unsigned n = 0; n < 8; ++n) { // place up to 8 bookcases
		bool const added(add_bookcase_to_room(rgen, room, zval, room_id, tot_light_amt, is_lit, objs_start));
		if (added) {++num_added;} else {break;}
	}
	return (num_added > 0);
}

void building_t::place_book_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, bool is_lit, bool use_dim_dir) {
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
	uint8_t const flags((is_lit ? RO_FLAG_LIT : 0) | RO_FLAG_NOCOLL);
	vector<room_object_t> &objs(interior->room_geom->objs);
	objs.emplace_back(book, TYPE_BOOK, room_id, dim, dir, flags, tot_light_amt, room_obj_shape::SHAPE_CUBE, color);
	objs.back().obj_id = (uint16_t)objs.size();
}

void building_t::add_rug_to_room(rand_gen_t &rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit) {
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
	uint8_t const obj_flags((is_lit ? RO_FLAG_LIT : 0) | RO_FLAG_NOCOLL);
	vector<room_object_t> &objs(interior->room_geom->objs);
	objs.emplace_back(rug, TYPE_RUG, room_id, 0, 0, obj_flags, tot_light_amt);
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
	if (is_cube_close_to_doorway(tc, 0.0, inc_open)) return 0; // bad placement
	if ((room.has_stairs || room.has_elevator) && interior->is_blocked_by_stairs_or_elevator_no_expand(tc, 4.0*wall_thickness)) return 0; // check stairs and elevators
	if (!inc_open && is_cube_close_to_doorway(tc, 0.0, 1)) return 2; // success, but could be better
	return 1; // success
}

bool building_t::hang_pictures_in_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
	if (!room_object_t::enable_pictures()) return 0; // disabled
	
	if (!is_house && !room.is_office) {
		if (room.is_hallway) return 0; // no pictures or whiteboards in office building hallways (what about rooms with stairs?)
		// room in a commercial building - add whiteboard when there is a full wall to use
	}
	if (room.is_sec_bldg) return 0; // no pictures in secondary buildings
	assert(room.part_id < parts.size());
	cube_t const &part(parts[room.part_id]);
	float const floor_height(get_window_vspace()), wall_thickness(get_wall_thickness());
	uint8_t const obj_flags((is_lit ? RO_FLAG_LIT : 0) | RO_FLAG_NOCOLL);
	vector<room_object_t> &objs(interior->room_geom->objs);
	bool was_hung(0);

	if (!is_house || room.is_office) { // add whiteboards
		if (rgen.rand_float() < 0.2) return 0; // skip 20% of the time
		bool const pref_dim(rgen.rand_bool()), pref_dir(rgen.rand_bool());

		for (unsigned dim2 = 0; dim2 < 2; ++dim2) {
			for (unsigned dir2 = 0; dir2 < 2; ++dir2) {
				bool const dim(bool(dim2) ^ pref_dim), dir(bool(dir2) ^ pref_dir);
				if (fabs(room.d[dim][dir] - part.d[dim][dir]) < 1.1*wall_thickness) continue; // on part boundary, likely exterior wall where there may be windows, skip
				cube_t c(room);
				c.z1() = zval + 0.25*floor_height; c.z2() = zval + 0.8*floor_height;
				c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*0.6*wall_thickness; // Note: offset by an additional half wall thickness
				c.expand_in_dim(!dim, -0.2*room.get_sz_dim(!dim)); // xy_space
				if (!check_valid_picture_placement(room, c, 0.6*room.get_sz_dim(!dim), zval, dim, dir, objs_start)) continue;
				objs.emplace_back(c, TYPE_WBOARD, room_id, dim, !dir, obj_flags, tot_light_amt); // whiteboard faces dir opposite the wall
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
			objs.emplace_back(best_pos, TYPE_PICTURE, room_id, dim, !dir, obj_flags, tot_light_amt); // picture faces dir opposite the wall
			objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix + 61*dim + 123*dir); // determines picture texture
			was_hung = 1;
		} // for dir
	} // for dim
	return was_hung;
}

void building_t::add_bathroom_windows(room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit) {
	unsigned num_ext_walls(0);

	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {num_ext_walls += (classify_room_wall(room, zval, dim, dir, 1) == ROOM_WALL_EXT);}
	}
	if (num_ext_walls != 1) return; // it looks odd to have window block walls at the corner of a building, so only enable this for single exterior walls
	float const floor_thickness(get_floor_thickness()), wall_thickness(get_wall_thickness()), window_thickness(0.05*wall_thickness);
	float const z2(zval + get_window_vspace() - floor_thickness);
	vector<room_object_t> &objs(interior->room_geom->objs);
	uint8_t const obj_flags((is_lit ? RO_FLAG_LIT : 0) | RO_FLAG_NOCOLL);

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
			objs.emplace_back(c, TYPE_WINDOW, room_id, dim, dir, obj_flags, tot_light_amt, room_obj_shape::SHAPE_CUBE, WHITE); // always lit
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
void building_t::gen_room_details(rand_gen_t &rgen, vect_cube_t const &ped_bcubes) {

	assert(interior);
	if (interior->room_geom) return; // already generated?
	//timer_t timer("Gen Room Details");
	interior->room_geom.reset(new building_room_geom_t(bcube.get_llc()));
	vector<room_object_t> &objs(interior->room_geom->objs);
	vector<room_t> &rooms(interior->rooms);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const extra_bathroom_prob(is_house ? 0.5 : 0.25);
	interior->room_geom->obj_scale = window_vspacing; // used to scale room object textures
	unsigned tot_num_rooms(0), num_light_stacks(0), num_bathrooms(0);
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {tot_num_rooms += calc_num_floors(*r, window_vspacing, floor_thickness);}
	objs.reserve(tot_num_rooms); // placeholder - there will be more than this many
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
		bool is_office_bathroom(r->is_office && r->rtype == RTYPE_BATH && !(r->has_stairs || r->has_elevator));
		if (is_office_bathroom) {cout << (r->has_stairs || r->has_elevator);}
		float light_size(floor_thickness); // default size for houses
		unsigned const room_objs_start(objs.size());

		if (r->is_sec_bldg) {
			if    (has_garage) {r->assign_to(RTYPE_GARAGE);}
			else if (has_shed) {r->assign_to(RTYPE_SHED);}
		}
		if (r->is_office) { // light size varies by office size
			float const room_size(r->dx() + r->dy()); // normalized to office size
			light_size = max(0.015f*room_size, 0.67f*floor_thickness);
		}
		if (r->is_hallway) { // light size varies by hallway size
			float const room_size(min(r->dx(), r->dy())); // normalized to hallway width
			light_size = max(0.06f*room_size, 0.67f*floor_thickness);
		}
		if (r->has_stairs) {r->assign_to(RTYPE_STAIRS);}
		float const light_val(22.0*light_size), room_light_intensity(light_val*light_val/r->get_area_xy()); // average for room, unitless
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
		// make chair colors consistent for each part by using a few variables for a hash
		colorRGBA chair_colors[12] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN};
		colorRGBA chair_color(chair_colors[(13*r->part_id + 123*tot_num_rooms + 617*mat_ix + 1367*num_floors) % 12]);
		unsigned num_lights_added(0);

		// place objects on each floor for this room
		for (unsigned f = 0; f < num_floors; ++f, z += floor_height) {
			room_center.z = z + fc_thick; // floor height
			bool const top_floor(f+1 == num_floors), check_stairs(!is_house && parts.size() > 1 && top_floor); // top floor of building that may have stairs connecting to upper stack
			bool is_lit(0), light_dim(room_dim), has_stairs(r->has_stairs);

			if (!has_stairs && (f == 0 || top_floor) && interior->stairwells.size() > 1) { // check for stairwells connecting stacked parts
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
				light.z1() = light.z2() - 0.5*fc_thick;
				is_lit = (r->is_hallway || ((rgen.rand() & (top_of_stairs ? 3 : 1)) != 0)); // 50% of lights are on, 75% for top of stairs, 100% for hallways
				//is_lit |= r->is_sec_bldg; // lit for garages and sheds

				// check ped_bcubes and set is_lit if any are people are in this floor of this room
				for (auto p = ped_bcubes.begin(); p != ped_bcubes.end() && !is_lit; ++p) {
					if (!p->intersects_xy(*r)) continue; // person not in this room
					if (p->z2() < light.z1() && p->z1() + floor_height > light.z2()) {is_lit = 1;} // on this floor
				}
				uint8_t flags(RO_FLAG_NOCOLL); // no collision detection with lights
				if (is_lit)        {flags |= RO_FLAG_LIT;}
				if (top_of_stairs) {flags |= RO_FLAG_TOS;}
				if (has_stairs)    {flags |= RO_FLAG_RSTAIRS;}
				colorRGBA color;
				if (is_house) {color = colorRGBA(1.0, 1.0, 0.9);} // house - yellowish
				else if (r->is_hallway || r->is_office) {color = colorRGBA(0.9, 0.9, 1.0);} // office building - blueish
				else {color = colorRGBA(1.0, 1.0, 1.0);} // white - small office
				unsigned num_lights(r->num_lights);

				if (r->is_hallway && num_lights > 1) { // place a light on each side (of the stairs if they exist), and also between stairs and elevator if there are both
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
				else { // normal room
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
			if (is_lit) {tot_light_amt += room_light_intensity;} // light surface area divided by room surface area with some fudge constant

			if (r->no_geom) {
				if (is_house && r->is_hallway) { // allow pictures and rugs in the hallways of houses
					hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs.size());
					if (rgen.rand_bool()) {add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit);}
				}
				continue; // no other geometry for this room
			}
			if (has_stairs && !pri_hall.is_all_zeros()) continue; // no other geometry in office building rooms that have stairs
			unsigned const objs_start(objs.size()), floor_mask(1<<f);
			bool added_tc(0), added_obj(0), can_place_book(0), is_bathroom(0), is_bedroom(0), is_kitchen(0), is_living(0);
			unsigned num_chairs(0);

			// place room objects
			bool const allow_br(!is_house || must_be_bathroom || f > 0 || num_floors == 1 || (rgen.rand_float() < 0.33f*(added_living + (added_kitchen_mask&1) + 1))); // bed/bath

			if (is_office_bathroom) { // bathroom is already assigned (should this be a row of toilets and sinks if the room is large?)
				added_obj = is_bathroom = added_bathroom = add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start); // add bathroom
			}
			if (!added_obj && allow_br && can_be_bedroom_or_bathroom(*r, (f == 0))) { // bedroom or bathroom case; need to check first floor even if is_cand_bathroom
				// place a bedroom 75% of the time unless this must be a bathroom; if we got to the second floor and haven't placed a bedroom, always place it; houses only
				if (is_house && !must_be_bathroom && ((f > 0 && !added_bedroom) || rgen.rand_float() < 0.75)) {
					added_obj = added_bedroom = is_bedroom = add_bed_to_room(rgen, *r, ped_bcubes, room_center.z, room_id, tot_light_amt, is_lit);
					if (is_bedroom) {r->assign_to(RTYPE_BED, f);}
					// Note: can't really mark room type as bedroom because it varies per floor; for example, there may be a bedroom over a living room connected to an exterior door
				}
				if (!added_obj && (must_be_bathroom || (can_be_bathroom(*r) && (num_bathrooms == 0 || rgen.rand_float() < extra_bathroom_prob)))) {
					// bathrooms can be in both houses and office buildings
					added_obj = is_bathroom = added_bathroom = add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start); // add bathroom
					if (is_bathroom) {r->assign_to(RTYPE_BATH, f);}
				}
			}
			if (!added_obj && rgen.rand_float() < (r->is_office ? 0.6 : (is_house ? 0.95 : 0.5))) {
				// place a table and maybe some chairs near the center of the room if it's not a hallway;
				// 60% of the time for offices, 95% of the time for houses, and 50% for other buildings
				unsigned const num_tcs(add_table_and_chairs(rgen, *r, ped_bcubes, room_id, room_center, chair_color, 0.1, tot_light_amt, is_lit));
				if (num_tcs > 0) {added_tc = added_obj = can_place_book = 1; num_chairs = num_tcs - 1;}
				// on ground floor, try to make this a kitchen; not all houses will have a kitchen with this logic - maybe we need fewer bedrooms?
				if (added_tc && !(added_kitchen_mask & floor_mask) && (!is_house || f == 0)) { // office buildings can also have kitchens, even on non-ground floors
					is_kitchen = add_kitchen_objs(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start, added_living);
					if (is_kitchen) {r->assign_to(RTYPE_KITCHEN, f); added_kitchen_mask |= floor_mask;}
				}
			}
			if (!added_obj && r->is_office) {added_obj = create_office_cubicles(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start);} // handle large offices

			if (!added_obj) { // try to place a desk if there's no table or bed
				added_obj = can_place_book = add_desk_to_room(rgen, *r, ped_bcubes, chair_color, room_center.z, room_id, tot_light_amt, is_lit);
				if (added_obj && !r->has_stairs) {r->assign_to((is_house ? RTYPE_STUDY : RTYPE_OFFICE), f);} // or other room type - may be overwritten below
			}
			if (is_house && can_place_book && !is_kitchen && f == 0) { // don't add second living room unless we added a kitchen
				if (((!added_living && (added_kitchen_mask || rgen.rand_bool())) || is_room_adjacent_to_ext_door(*r, 1))) { // front_door_only=1
					// add a living room on the ground floor if it has a table or desk but isn't a kitchen
					added_living = is_living = add_livingroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start);
					if (is_living) {r->assign_to(RTYPE_LIVING, f);}
				}
			}
			if (can_place_book) { // an object was placed (table or desk), maybe add a book on top of it
				if (rgen.rand_float() < (added_tc ? 0.4 : 0.75)*(is_house ? 1.0 : 0.5)*(r->is_office ? 0.75 : 1.0)) {
					assert(objs.size() > objs_start);
					place_book_on_obj(rgen, objs[objs_start], room_id, tot_light_amt, is_lit, !added_tc);
				}
			}
			if (is_house) { // place house-specific items
				if (!is_bathroom && !is_kitchen && rgen.rand_float() < 0.8) { // place bookcase 80% of the time, but not in bathrooms or kitchens
					add_bookcase_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start);
				}
				if (!has_stairs && (rgen.rand()&3) <= (added_tc ? 0 : 2)) { // maybe add a rug, 25% of the time if there's a table and 75% of the time otherwise
					add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit);
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
				}
				else if (!added_library && add_library_objs(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start)) { // add library, at most one
					r->assign_to(RTYPE_LIBRARY, f);
					added_library = 1;
				}
			}
			if (r->rtype == RTYPE_NOTSET) {
				if (is_room_adjacent_to_ext_door(*r)) {r->assign_to(RTYPE_ENTRY, f);} // entryway if has exterior door and is unassigned
				else if (!is_house) {r->assign_to(RTYPE_OFFICE, f);} // any unset room in an office building is an office
			}
			bool const can_hang(is_house || !(is_bathroom || is_kitchen)); // no whiteboards in office bathrooms or kitchens
			bool const was_hung(can_hang && hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start));

			if (is_bathroom || is_kitchen || rgen.rand_float() < 0.8) { // 80% of the time, always in bathrooms and kitchens
				add_trashcan_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start, (was_hung && r->is_office)); // no trashcans on same wall as office whiteboard
			}
			if (is_bathroom) {add_bathroom_windows(*r, room_center.z, room_id, tot_light_amt, is_lit);} // find all windows and add frosted windows
			//if (z == bcube.z1()) {} // any special logic that goes on the first floor is here
		} // for f (floor)
		num_light_stacks += num_lights_added;
		if (added_bathroom) {++num_bathrooms;}

		// determine if room is interior and tag objects
		assert(r->part_id < parts.size());
		r->interior = parts[r->part_id].contains_cube_xy_no_adj(*r);

		if (r->interior) {
			for (auto i = objs.begin() + room_objs_start; i != objs.end(); ++i) {i->flags |= RO_FLAG_INTERIOR;}
		}
	} // for r (room)
	add_stairs_and_elevators(rgen);
	objs.shrink_to_fit();
	interior->room_geom->light_bcubes.resize(num_light_stacks); // allocate but don't fill un until needed
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness());
	float const stair_dz(window_vspacing/(NUM_STAIRS_PER_FLOOR+1)), stair_height(stair_dz + floor_thickness);
	vector<room_object_t> &objs(interior->room_geom->objs);
	interior->room_geom->stairs_start = objs.size();
	interior->room_geom->has_elevators = (!interior->elevators.empty());

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator) continue; // for elevator, not stairs
		bool const dim(i->dim), dir(i->dir);
		// Note: stairs always start at floor_thickness above the landing z1, ignoring landing z2/height
		float const tot_len(i->get_sz_dim(dim)), floor_z(i->z1() + floor_thickness - window_vspacing), step_len_pos(tot_len/NUM_STAIRS_PER_FLOOR);
		float step_len((dir ? 1.0 : -1.0)*step_len_pos), z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i);

		if (i->shape == SHAPE_STRAIGHT || i->shape == SHAPE_WALLED) { // straight stairs
			for (unsigned n = 0; n < NUM_STAIRS_PER_FLOOR; ++n, z += stair_dz, pos += step_len) {
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // Note: room_id=0, not tracked, unused
			}
		}
		else if (i->shape == SHAPE_U) { // U-shaped stairs
			bool const side(0); // for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
			stair.d[!dim][side] = i->get_center_dim(!dim);
			step_len *= 2.0;

			for (unsigned n = 0; n < NUM_STAIRS_PER_FLOOR; ++n, z += stair_dz, pos += step_len) {
				if (n == NUM_STAIRS_PER_FLOOR/2) { // reverse direction and switch to other side
					step_len *= -1.0;
					stair.d[!dim][ side] = i->d[!dim][side];
					stair.d[!dim][!side] = i->get_center_dim(!dim);
				}
				bool const is_rev(n >= NUM_STAIRS_PER_FLOOR/2);
				stair.d[dim][dir^is_rev^1] = pos; stair.d[dim][dir^is_rev] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, side^is_rev); // Note: room_id=0, not tracked, unused
				objs.back().shape = SHAPE_STAIRS_U;
			} // for n
		}
		else {assert(0);} // unknown stairs shape

		if (i->shape == SHAPE_U || i->shape == SHAPE_WALLED) { // add walls around stairs for this floor
			float const wall_hw(0.15*step_len_pos), half_thick(0.5*floor_thickness);
			stair = *i;
			stair.z2() -= 0.5*floor_thickness; // prevent z-fighting on top floor
			stair.z1()  = max(bcube.z1()+half_thick, floor_z-half_thick); // full height
			set_wall_width(stair, i->d[dim][dir], wall_hw, dim);
			objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // back/end of stairs
			stair.d[dim][!dir] = i->d[dim][!dir];

			for (unsigned d = 0; d < 2; ++d) { // sides of stairs
				set_wall_width(stair, i->d[!dim][d], wall_hw, !dim);
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir);
			}
		}
	} // for i
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		unsigned const elevator_id(i - interior->elevators.begin()); // used for room_object_t::room_id
		cube_t elevator_car(*i);
		elevator_car.z1() += 0.05*get_floor_thickness(); // to prevent z-fighting when looking at the building from the bottom
		elevator_car.z2() = elevator_car.z1() + window_vspacing; // currently at the bottom floor
		objs.emplace_back(elevator_car, TYPE_ELEVATOR, elevator_id, i->dim, i->dir, (i->open ? RO_FLAG_OPEN : 0));
	}
}

void building_t::draw_room_geom(shader_t &s, vector3d const &xlate, bool shadow_only, bool inc_small, bool player_in_building) {
	if (interior && interior->room_geom) {interior->room_geom->draw(s, xlate, shadow_only, inc_small, player_in_building);}
}
void building_t::gen_and_draw_room_geom(shader_t &s, vector3d const &xlate, vect_cube_t &ped_bcubes, unsigned building_ix,
	int ped_ix, bool shadow_only, bool inc_small, bool player_in_building)
{
	if (!interior) return;
	if (is_rotated()) return; // no room geom for rotated buildings

	if (!has_room_geom()) {
		rand_gen_t rgen;
		rgen.set_state(building_ix, parts.size()); // set to something canonical per building
		ped_bcubes.clear();
		if (ped_ix >= 0) {get_ped_bcubes_for_building(ped_ix, building_ix, ped_bcubes);}
		gen_room_details(rgen, ped_bcubes); // generate so that we can draw it
		assert(has_room_geom());
	}
	draw_room_geom(s, xlate, shadow_only, inc_small, player_in_building);
}

void building_t::clear_room_geom() {
	if (!has_room_geom()) return;
	interior->room_geom->clear(); // free VBO data before deleting the room_geom object
	interior->room_geom.reset();
}

room_t::room_t(cube_t const &c, unsigned p, unsigned nl, bool is_hallway_, bool is_office_, bool is_sec_bldg_) :
	cube_t(c), has_stairs(0), has_elevator(0), no_geom(is_hallway_), is_hallway(is_hallway_), is_office(is_office_),
	is_sec_bldg(is_sec_bldg_), interior(0), ext_sides(0), part_id(p), num_lights(nl), lit_by_floor(0) // no geom in hallways
{
	if      (is_sec_bldg) {rtype = RTYPE_GARAGE;} // or RTYPE_SHED - will be set later
	else if (is_hallway)  {rtype = RTYPE_HALL;}
	else if (is_office)   {rtype = RTYPE_OFFICE;}
	else if (has_stairs)  {rtype = RTYPE_STAIRS;}
	else                  {rtype = RTYPE_NOTSET;}
}
void room_t::assign_to(room_type rt, unsigned floor) {
	if (rtype == RTYPE_NOTSET || floor == 0) {rtype = rt;} // rtype is not per floor: assign if not yet assigned or this is the first floor (first floor has priority)
}

