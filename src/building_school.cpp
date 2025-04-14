// 3D World - Hospital Buildings
// by Frank Gennari 2/26/25

#include "function_registry.h"
#include "buildings.h"
//#include "city_model.h"


bool building_t::add_classroom_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	float tot_light_amt, unsigned objs_start, colorRGBA const &chair_color, unsigned &td_orient)
{
	if (room.has_elevator) return 0; // no classroom in an elevator
	float const vspace(get_window_vspace()), wall_thickness(get_wall_thickness());
	vector2d const room_sz(room.get_size_xy());
	if (room_sz.get_max_val() < 3.0*vspace || room_sz.get_min_val() < 1.8*vspace) return 0; // room is too small
	bool const dim(room_sz.x < room_sz.y); // long dim
	// front_dir is the side of the long dim not by the door
	bool valid_dirs[2] = {0,0}, door_sides[2] = {0,0};
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors

	for (unsigned dir = 0; dir < 2; ++dir) {
		cube_t fb_wall(room), side_wall(room);
		set_wall_width(fb_wall,   room.d[ dim][dir], wall_thickness,  dim);
		set_wall_width(side_wall, room.d[!dim][dir], wall_thickness, !dim);
		valid_dirs[dir] = !has_bcube_int(fb_wall,   doorways); // skip if there's a door on this wall
		door_sides[dir] =  has_bcube_int(side_wall, doorways);
	}
	bool dir(0); // front dir
	if (!valid_dirs[0] && !valid_dirs[1]) return 0; // no valid door-free walls
	if ( valid_dirs[0] &&  valid_dirs[1]) { // both ends are valid
		bool const ext[2] = {(classify_room_wall(room, zval, dim, 0, 0) == ROOM_WALL_EXT), (classify_room_wall(room, zval, dim, 1, 0) == ROOM_WALL_EXT)};
		if (ext[0] != ext[1]) {dir = ext[0];} // choose interior wall so that we can place a chalkboard behind the desk
		else {dir = rgen.rand_bool();} // choose a random end
	}
	else {dir = valid_dirs[1];} // only one valid wall

	float const clearance(get_min_front_clearance_inc_people());
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const td_width(0.8*vspace*rgen.rand_uniform(1.1, 1.2)), td_depth(0.38*vspace*rgen.rand_uniform(1.1, 1.2)), td_height(0.23*vspace*rgen.rand_uniform(1.12, 1.2));
	float const front_wall_pos(room_bounds.d[dim][dir]), room_center(room.get_center_dim(!dim)), dsign(dir ? -1.0 : 1.0);
	float const desk_front_pos(front_wall_pos + dsign*(max(1.2*clearance, 0.25*vspace) + wall_thickness)); // near wall, with space for chair
	float const desk_back_pos(desk_front_pos + dsign*td_depth);
	float const desk_width(0.48*vspace), desk_depth(0.34*vspace), desk_height(0.25*vspace);
	cube_t student_area(room_bounds);
	student_area.d[dim][ dir]  = desk_back_pos + dsign*1.3*clearance; // front side near teacher
	student_area.d[dim][!dir] -= dsign*0.5*clearance; // back wall

	for (unsigned dir = 0; dir < 2; ++dir) {
		student_area.d[!dim][dir] -= (dir ? 1.0 : -1.0)*(door_sides[dir] ? 1.0 : -0.45)*clearance; // leave space for doors but not side walls
	}
	float desk_wspacing(desk_width + clearance), desk_dspacing(desk_depth + clearance);
	float const avail_width(student_area.get_sz_dim(!dim)), avail_len(student_area.get_sz_dim(dim));
	unsigned const ncols(avail_width/desk_wspacing), nrows(avail_len/desk_dspacing);
	if (nrows < 1 || ncols < 1) return 0; // not enough space for student desks; shouldn't happen

	// place teacher desk at front
	cube_t tdesk;
	set_cube_zvals(tdesk, zval, zval+td_height);
	tdesk.d[dim][ dir] = desk_front_pos;
	tdesk.d[dim][!dir] = desk_back_pos ;
	set_wall_width(tdesk, room_center, 0.5*td_width, !dim);
	if (!add_classroom_desk(rgen, room, tdesk, room_id, tot_light_amt, chair_color, dim, dir, 0)) return 0; // desk_ix=0
	td_orient = 2*dim + dir; // place chalkboard behind the teacher desk
	// place rows and columns of student desks
	float const first_row_pos(student_area.d[dim][dir]);
	desk_wspacing = avail_width/ncols;
	desk_dspacing = avail_len  /nrows;

	for (unsigned row = 0; row < nrows; ++row) {
		for (unsigned col = 0; col < ncols; ++col) {
			cube_t desk;
			set_cube_zvals(desk, zval, zval+desk_height);
			desk.d[dim][ dir] = first_row_pos + dsign*desk_dspacing*row;
			desk.d[dim][!dir] = desk.d[dim][dir] + dsign*desk_depth;
			set_wall_width(desk, (student_area.d[!dim][0] + desk_wspacing*(col + 0.5)), 0.5*desk_width, !dim);
			if (!add_classroom_desk(rgen, room, desk, room_id, tot_light_amt, chair_color, dim, !dir, (1 + col + ncols*row))) continue;
		} // for col
	} // for row
	add_numbered_door_sign("Classroom ", room, zval, room_id, floor_ix);
	return 1;
}

bool building_t::add_classroom_desk(rand_gen_t &rgen, room_t const &room, cube_t const &desk, unsigned room_id, float tot_light_amt,
	colorRGBA const &chair_color, bool dim, bool dir, unsigned desk_ix)
{
	vect_room_object_t &objs(interior->room_geom->objs);
	if (is_obj_placement_blocked(desk, room, 1)) return 0; // check proximity to doors, etc.
	unsigned const pre_add_obj_id(objs.size());
	unsigned const flags((desk_ix == 0) ? RO_FLAG_HAS_EXTRA : 0); // teacher's desk always has drawers
	objs.emplace_back(desk, TYPE_DESK, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE); // no tall desks
	set_obj_id(objs);
	objs.back().obj_id += 123*desk_ix; // more random variation
	// add paper, pens, and pencils
	add_papers_to_surface      (desk, dim,  dir, 7, rgen, room_id, tot_light_amt); // add 0-7 sheet(s) of paper
	add_pens_pencils_to_surface(desk, dim, !dir, 4, rgen, room_id, tot_light_amt); // 0-4 pens/pencils
	// add chair
	point chair_pos;
	chair_pos.z     = desk.z1();
	chair_pos[ dim] = desk.d[dim][dir]; // front of desk
	chair_pos[!dim] = desk.get_center_dim(!dim) + 0.1*rgen.signed_rand_float()*desk.get_sz_dim(dim); // slightly misaligned
	
	if (!add_chair(rgen, room, vect_cube_t(), room_id, chair_pos, chair_color, dim, !dir, tot_light_amt, 0, 0, 0, 0, 2, 1)) { // no blockers, reduced_clearance=1
		objs.resize(pre_add_obj_id); // no chair; remove the desk as well, since it may be too close to a door
	}
	return 1;
}

void building_t::add_objects_next_to_classroom_chalkboard(rand_gen_t &rgen, room_object_t const &cb, room_t const &room, float zval, unsigned objs_start) {
	// add US flag; flags are two sided, so lighting doesn't look correct on the unlit side
	bool const dim(cb.dim), dir(cb.dir), side(dir ^ dim); // side of clock; always place on the right because the left side is lit with correct normals
	float const flag_pos(0.5*(cb.d[!dim][side] + room.d[!dim][side])), wall_pos(cb.d[dim][dir]); // halfway between edge of chalkboard and edge of room
	add_wall_us_flag(wall_pos, flag_pos, zval, dim, dir, cb.room_id, cb.light_amt);
	if (is_obj_placement_blocked(interior->room_geom->objs.back(), room, 1)) {interior->room_geom->objs.pop_back();} // remove if invalid placement; inc_open_doors=1
	// add clock to the left side
	float const floor_spacing(get_window_vspace()), clock_sz(0.16*floor_spacing), clock_z1(zval + get_floor_ceil_gap() - 1.4*clock_sz), clock_depth(0.08*clock_sz);
	cube_t clock;
	set_cube_zvals(clock, clock_z1, clock_z1+clock_sz);
	set_wall_width(clock, 0.5*(cb.d[!dim][!side] + room.d[!dim][!side]), 0.5*clock_sz, !dim);
	clock.d[dim][!dir] = wall_pos;
	clock.d[dim][ dir] = wall_pos + (dir ? 1.0 : -1.0)*clock_depth;
	cube_t tc(clock);
	tc.d[dim][dir] += (dir ? 1.0 : -1.0)*0.25*floor_spacing; // add clearance
	if (is_obj_placement_blocked(tc, room, 1)) return; // inc_open_doors=1
	add_clock(clock, cb.room_id, cb.light_amt, dim, dir, 0); // digital=0
	// add plants after chalkboard to avoid blocking it
	unsigned const num_plants(rgen.rand() % 30); // 0-2
	add_plants_to_room(rgen, room, zval, cb.room_id, cb.light_amt, objs_start, num_plants);
}

