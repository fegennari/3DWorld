// 3D World - Hospital Buildings
// by Frank Gennari 2/26/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

extern object_model_loader_t building_obj_model_loader;


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
	bool const add_bottles(0), add_trash(rgen.rand_float() < 0.4), add_papers(rgen.rand_float() < 0.4), add_glass(0);
	add_floor_clutter_objs(rgen, room, room_bounds, zval, room_id, tot_light_amt, objs_start, add_bottles, add_trash, add_papers, add_glass);
	add_numbered_door_sign("Classroom ", room, zval, room_id, floor_ix);
	return 1;
}

bool building_t::add_classroom_desk(rand_gen_t &rgen, room_t const &room, cube_t const &desk, unsigned room_id, float tot_light_amt,
	colorRGBA const &chair_color, bool dim, bool dir, unsigned desk_ix)
{
	vect_room_object_t &objs(interior->room_geom->objs);
	if (is_obj_placement_blocked(desk, room, 1)) return 0; // check proximity to doors, etc.
	bool const teacher_desk(desk_ix == 0);
	unsigned const desk_obj_ix(objs.size()), flags(teacher_desk ? RO_FLAG_HAS_EXTRA : 0); // teacher's desk always has drawers
	objs.emplace_back(desk, TYPE_DESK, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE); // no tall desks
	set_obj_id(objs);
	objs.back().obj_id += 123*desk_ix; // more random variation
	// add paper, pens, and pencils
	unsigned const objs_start(objs.size()); // excludes the desk
	if (rgen.rand_float() < 0.7) {add_papers_to_surface      (desk, dim,  dir, 7, rgen, room_id, tot_light_amt);} // add 0-7 sheet(s) of paper
	if (rgen.rand_float() < 0.7) {add_pens_pencils_to_surface(desk, dim, !dir, 4, rgen, room_id, tot_light_amt);} // 0-4 pens/pencils

	if (teacher_desk && rgen.rand_float() < 0.33) { // add a cup on the desk 33% of the time
		vect_cube_t const avoid(objs.begin()+objs_start, objs.end()); // add all papers, pens, and pencils
		place_cup_on_obj(rgen, desk, room_id, tot_light_amt, avoid);
	}
	if (teacher_desk && rgen.rand_float() < 0.5) { // place an apple on the desk
		vect_cube_t const avoid(objs.begin()+objs_start, objs.end()); // add all papers, pens, pencils, and cups
		place_apple_on_obj(rgen, desk, room_id, tot_light_amt, avoid);
	}
	if (is_school() && rgen.rand_float() < 0.67) { // maybe add a book on the desk; often skipped due to overlaps; only for schools (not hospitals)
		place_book_on_obj(rgen, objs[desk_obj_ix], room_id, tot_light_amt, objs_start, 1, RO_FLAG_USED, 1); // use_dim_dir=1; skip_if_overlaps=1
	}
	// add chair
	point chair_pos(0.0, 0.0, desk.z1());
	chair_pos[ dim] = desk.d[dim][dir]; // front of desk
	chair_pos[!dim] = desk.get_center_dim(!dim) + 0.1*rgen.signed_rand_float()*desk.get_sz_dim(dim); // slightly misaligned
	
	if (!add_chair(rgen, room, vect_cube_t(), room_id, chair_pos, chair_color, dim, !dir, tot_light_amt, 0, 0, 0, 0, 2, 1)) { // no blockers, reduced_clearance=1
		objs.resize(desk_obj_ix); // no chair; remove the desk as well, and any objects placed on it, since it may be too close to a door
	}
	return 1;
}

void building_t::add_objects_next_to_classroom_chalkboard(rand_gen_t &rgen, room_object_t const &cb, room_t const &room, float zval, unsigned objs_start) {
	// add US flag; flags are two sided, so lighting doesn't look correct on the unlit side
	bool const dim(cb.dim), dir(cb.dir), side(dir ^ dim); // side of clock; always place on the right because the left side is lit with correct normals
	float const flag_pos(0.5*(cb.d[!dim][side] + room.d[!dim][side])), wall_pos(cb.d[dim][dir]); // halfway between edge of chalkboard and edge of room
	add_wall_us_flag(wall_pos, flag_pos, zval, dim, dir, cb.room_id, cb.light_amt);
	vect_room_object_t &objs(interior->room_geom->objs);
	if (is_obj_placement_blocked(objs.back(), room, 1)) {objs.pop_back();} // remove if invalid placement; inc_open_doors=1
	// add clock to the left side
	float const floor_spacing(get_window_vspace()), clock_sz(0.16*floor_spacing), clock_z1(zval + get_floor_ceil_gap() - 1.4*clock_sz), clock_depth(0.08*clock_sz);
	cube_t clock;
	set_cube_zvals(clock, clock_z1, clock_z1+clock_sz);
	set_wall_width(clock, 0.5*(cb.d[!dim][!side] + room.d[!dim][!side]), 0.5*clock_sz, !dim);
	clock.d[dim][!dir] = wall_pos;
	clock.d[dim][ dir] = wall_pos + (dir ? 1.0 : -1.0)*clock_depth;
	cube_t tc(clock);
	tc.d[dim][dir] += (dir ? 1.0 : -1.0)*0.5*floor_spacing; // add clearance

	if (overlaps_obj_or_placement_blocked(tc, room, objs_start)) { // bad placement, try shifting down (below vent, etc.)
		clock.translate_dim(2, -0.1*floor_spacing);
		set_cube_zvals(tc, clock.z1(), clock.z2());
	}
	if (!overlaps_obj_or_placement_blocked(tc, room, objs_start)) {add_clock(clock, cb.room_id, cb.light_amt, dim, dir, 0);} // digital=0
	// add plants after chalkboard to avoid blocking it
	unsigned const num_plants(rgen.rand() % 3); // 0-2
	add_plants_to_room(rgen, room, zval, cb.room_id, cb.light_amt, objs_start, num_plants);
}

void building_t::add_hallway_lockers(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start) {
	bool const dim(room.dx() < room.dy()); // hallway dim
	cube_t const room_bounds(get_walkable_room_bounds(room));
	cube_t place_area(room_bounds);
	place_area.expand_in_dim(dim, -0.75*get_window_vspace()); // leave space at the ends for windows, etc.
	bool const add_padlocks(floor_ix == 0); // first floor only, to avoid having too many model objects (but may be added to stacked parts)
	add_room_lockers(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area, RTYPE_HALL, dim, 0, add_padlocks); // dir_skip_mask=0
	bool const add_bottles(rgen.rand_float() < 0.35), add_trash(rgen.rand_float() < 0.75), add_papers(rgen.rand_float() < 0.5), add_glass(0);
	add_floor_clutter_objs(rgen, room, room_bounds, zval, room_id, tot_light_amt, objs_start, add_bottles, add_trash, add_papers, add_glass);
}

bool building_t::add_room_lockers(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
	cube_t const &place_area, room_type rtype, bool dim, int dir_skip_mask, bool add_padlocks)
{
	float const floor_spacing(get_window_vspace()), locker_height(0.75*floor_spacing), locker_depth(0.25*locker_height);
	float locker_width(0.22*locker_height);
	vect_room_object_t &objs(interior->room_geom->objs);
	float const room_len(place_area.get_sz_dim(dim)); // long dim
	unsigned const lockers_start(objs.size()), num_lockers(room_len/locker_width); // floor
	// if there are enough lockers, increase their width slightly so that lockers tile to fill the exact wall length
	if (num_lockers >= 10) {locker_width = room_len/num_lockers;}
	bool const add_blockers(rtype != RTYPE_HALL); // add blockers in front of rows of lockers, except for school hallways (which have splits for secondary hallways, etc.)
	bool const single_side (rtype == RTYPE_GYM ), first_dir(single_side ? rgen.rand_bool() : 0);
	// add expanded blockers for stairs, elevators, etc. to ensure there's space for the player and people to walk on the sides
	float const clearance(get_min_front_clearance_inc_people()), se_clearance(2.0*clearance);
	add_padlocks &= building_obj_model_loader.is_model_valid(OBJ_MODEL_PADLOCK);
	vector3d const sz(add_padlocks ? building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_PADLOCK) : zero_vector); // D, W, H
	colorRGBA const lock_color(0.4, 0.4, 0.4); // gray
	unsigned const num_locker_colors = 6;
	colorRGBA const locker_colors[num_locker_colors] =
	{colorRGBA(0.4, 0.6, 0.7), colorRGBA(0.4, 0.7, 0.6), colorRGBA(0.2, 0.5, 0.8), colorRGBA(0.7, 0.05, 0.05), colorRGBA(0.6, 0.45, 0.25), GRAY};
	colorRGBA const locker_color(locker_colors[(7*mat_ix + 11*room_id + 13*interior->rooms.size())%num_locker_colors]); // random per part/room
	vect_cube_t blockers;

	for (stairwell_t const &s : interior->stairwells) {
		if (room.contains_cube_xy_overlaps_z(s)) {blockers.push_back(s);}
	}
	for (elevator_t const &e : interior->elevators) {
		if (room.contains_cube_xy_overlaps_z(e)) {blockers.push_back(e);}
	}
	for (cube_t &c : blockers) {c.expand_by_xy(se_clearance);}
	cube_t locker;
	set_cube_zvals(locker, zval, zval+locker_height);
	unsigned lix(0);

	for (unsigned D = 0; D < 2; ++D) { // for each side of the room
		bool const d(bool(D) ^ first_dir);
		if (dir_skip_mask & (1 << (unsigned)d)) continue;
		if (zval >= ground_floor_z1 && classify_room_wall(room, zval, !dim, d, 0) == ROOM_WALL_EXT) continue; // skip exterior walls with windows
		float const dsign(d ? -1.0 : 1.0), wall_edge(place_area.d[!dim][d]), locker_front(wall_edge + dsign*locker_depth);
		locker.d[!dim][ d] = wall_edge;
		locker.d[!dim][!d] = locker_front;
		cube_t row_bc;

		for (unsigned n = 0; n < num_lockers; ++n) {
			float const pos(place_area.d[dim][0] + n*locker_width);
			locker.d[dim][0] = pos;
			locker.d[dim][1] = pos + locker_width;
			assert(locker.intersects(room));
			if (has_bcube_int(locker, blockers)) continue;
			cube_t test_cube(locker);
			test_cube.expand_in_dim(dim, 2.0*locker_width); // add some padding to the sides
			test_cube.intersect_with_cube(room); // not blocked by objects in an adjacent room
			test_cube.d[!dim][!d] += dsign*max(locker_width, clearance); // add space in front for the door to open and the player to walk by
			bool invalid(0);

			for (auto i = objs.begin()+objs_start; i != objs.begin()+lockers_start; ++i) { // can skip other lockers
				if (i->type != TYPE_FLOORING && i->intersects(test_cube)) {invalid = 1; break;}
			}
			if (invalid || is_obj_placement_blocked(test_cube, room, 1,    0)) continue;
			if (!check_if_placed_on_wall           (test_cube, room, !dim, d)) continue; // ensure the locker is against a wall
			unsigned flags(is_industrial() ? RO_FLAG_IN_FACTORY : 0);

			if (add_padlocks && rgen.rand_float() < 0.25) { // add padlocks to some lockers
				float const height(0.04*floor_spacing), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z);
				cube_t lock;
				lock.z1() = zval + 0.365*floor_spacing;
				lock.z2() = lock.z1() + height;
				set_wall_width(lock, (locker.d[dim][0] + ((dim ^ bool(d)) ? 0.175 : 0.825)*locker_width), hwidth, dim);
				float const pos(locker.d[!dim][!d]);
				lock.d[!dim][ d] = pos;
				lock.d[!dim][!d] = pos + dsign*depth;
				objs.emplace_back(lock, TYPE_PADLOCK, room_id, !dim, d, (RO_FLAG_NOCOLL | RO_FLAG_IS_ACTIVE), 1.0, SHAPE_CUBE, lock_color); // attached
				flags |= RO_FLAG_NONEMPTY; // flag as locked
			}
			objs.emplace_back(locker, TYPE_LOCKER, room_id, !dim, !d, flags, tot_light_amt, SHAPE_CUBE, locker_color, lix++);
			set_obj_id(objs); // for random contents
			objs.back().state_flags = rtype; // store room type for correct object type addtion
			if (add_blockers) {row_bc.assign_or_union_with_cube(locker);}
			interior->room_geom->has_locker = 1;
		} // for n
		if (row_bc.is_all_zeros()) continue; // no lockers
		row_bc.d[!dim][ d] = locker_front;
		row_bc.d[!dim][!d] = locker_front + dsign*max(locker_width, clearance);
		objs.emplace_back(row_bc, TYPE_BLOCKER, room_id, !dim, !d, RO_FLAG_INVIS);
		if (single_side) break;
	} // for d
	return (objs.size() > lockers_start); // true if at least one locker was added
}

bool building_t::add_locker_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	bool const dim(room.dx() < room.dy()); // long dim
	cube_t const room_bounds(get_walkable_room_bounds(room));
	cube_t place_area(room_bounds);
	place_area.expand_by(-get_trim_thickness()); // shrink to leave a small gap
	if (!add_room_lockers(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area, RTYPE_LOCKER, dim, 0, 1)) return 0; // dir_skip_mask=0, add_padlocks=1
	float floor_spacing(get_window_vspace());
	vect_room_object_t &objs(interior->room_geom->objs);
	
	if (place_area.get_sz_dim(!dim) > 1.0*floor_spacing) { // add benches in the center of the room if room is wide enough
		float const bench_height(0.22*floor_spacing), bench_width(rgen.rand_uniform(0.7, 0.9)*bench_height), bench_len(rgen.rand_uniform(4.0, 5.0)*bench_height);
		float const clearance(1.1*get_min_front_clearance_inc_people());
		float const lo_end(place_area.d[dim][0] + clearance), hi_end(place_area.d[dim][1] - clearance), gap(hi_end - lo_end);
		unsigned const num_benches(gap/(bench_len + 0.5*bench_width)); // floor, with some min space in between (but maybe not enough for the player to walk)

		if (num_benches > 0) { // always true?
			float const bench_spacing(gap/num_benches);
			cube_t bench;
			set_cube_zvals(bench, zval, zval + bench_height);
			set_wall_width(bench, place_area.get_center_dim(!dim), 0.5*bench_width, !dim); // short dim

			for (unsigned n = 0; n < num_benches; ++n) {
				set_wall_width(bench, (lo_end + (n + 0.5)*bench_spacing), 0.5*bench_len, dim); // long  dim
				cube_t test_cube(bench);
				test_cube.expand_by_xy(clearance);
				if (is_obj_placement_blocked(test_cube, room, 1)) continue;
				unsigned const item_flags(1); // flag as metal mesh
				objs.emplace_back(bench, TYPE_BENCH, room_id, !dim, 0, 0, tot_light_amt, SHAPE_CUBE, WHITE, item_flags); // dir=0, flags=0

				if (rgen.rand_float() < 0.65) { // add teeshirt or pants on the bench
					unsigned const type(rgen.rand_bool() ? TYPE_PANTS : TYPE_TEESHIRT);
					float const length(((type == TYPE_TEESHIRT) ? 1.0 : 0.8)*bench_width), width(0.98*length), height(0.01*length);
					vector3d size(0.5*length, 0.5*width, height);
					bool const dim2(rgen.rand_bool()), dir2(rgen.rand_bool()); // choose a random orientation
					if (dim2) {std::swap(size.x, size.y);}
					cube_t c(gen_xy_pos_in_area(bench, size, rgen, bench.z2()));
					c.expand_by_xy(size);
					c.z2() += size.z;
					unsigned const sp_flags(RO_FLAG_NOCOLL | RO_FLAG_HAS_EXTRA); // flag as extra to indicate alpha mask material
					objs.emplace_back(c, type, room_id, dim2, dir2, sp_flags, tot_light_amt, SHAPE_CUBE, gen_shirt_pants_color(type, rgen));
				}
			} // for n
		}
	}
	// maybe add a water fountain
	add_wall_water_fountain(rgen, room, zval, room_id, tot_light_amt, objs_start);
	// add shoes on the floor
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_SHOE)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SHOE)); // L, W, H
		float const length(0.2*floor_spacing*rgen.rand_uniform(0.75, 1.0)), width(length*sz.y/sz.x), height(length*sz.z/sz.x), hlen(0.5*length), pair_hw(width);
		unsigned const num_shoes(1 + (rgen.rand() % 3)); // 1-3
		cube_t shoes;
		set_cube_zvals(shoes, zval, zval+height);

		for (unsigned n = 0, num_added = 0; n < 20 && num_added < num_shoes; ++n) { // try to place shoes
			bool const rdir(rgen.rand_bool()); // wall dir
			cube_t shoe_area(place_area);
			shoe_area.d[dim][!rdir] = place_area.d[dim][rdir] + (rdir ? -1.0 : 1.0)*1.5*length; // shrink area to near the short wall
			set_wall_width(shoes, rgen.rand_uniform((shoe_area.d[ dim][0] + hlen   ), (shoe_area.d[ dim][1] - hlen   )), hlen,     dim);
			set_wall_width(shoes, rgen.rand_uniform((shoe_area.d[!dim][0] + pair_hw), (shoe_area.d[!dim][1] - pair_hw)), pair_hw, !dim);
			if (is_obj_placement_blocked(shoes, room, 1) || overlaps_other_room_obj(shoes, objs_start)) continue;
			bool const sdir(rgen.rand_bool()); // shoe dir
			unsigned const item_flags(rgen.rand()); // random shoe sub-model
			add_obj_pair(room_object_t(shoes, TYPE_SHOE, room_id, dim, sdir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, WHITE, item_flags), objs);
			++num_added;
		} // for n
	}
	if (is_school() && rgen.rand_float() < 0.75) { // add shirt on pants on the floor of school locker rooms
		unsigned const type(rgen.rand_bool() ? TYPE_PANTS : TYPE_TEESHIRT);
		place_shirt_pants_on_floor(rgen, room, zval, room_id, tot_light_amt, place_area, objs_start, type);
	}
	string sign_text("Locker Room");

	if (is_school()) { // girls vs. boys
		bool boys(0);
		if      (interior->room_geom->mens_count < interior->room_geom->womens_count) {boys = 1;}
		else if (interior->room_geom->mens_count > interior->room_geom->womens_count) {boys = 0;}
		else {boys = rgen.rand_bool();} // tied
		++(boys ? interior->room_geom->mens_count : interior->room_geom->womens_count);
		sign_text = (boys ? "Boys Locker" : "Girls Locker");
	}
	add_door_sign(sign_text, room, zval, room_id);
	return 1;
}

void building_t::add_wall_water_fountain(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	unsigned const wf_obj_ix(interior->room_geom->objs.size());

	if (place_model_along_wall(OBJ_MODEL_WFOUNTAIN, TYPE_WFOUNTAIN, room, 0.25, rgen, (zval+0.18*get_window_vspace()),
		room_id, tot_light_amt, get_room_wall_bounds(room), objs_start, 1.0, 4, 0, WHITE, 1)) // not_at_window=1
	{
		interior->room_geom->objs[wf_obj_ix].dir ^= 1; // placed dir was backwards
	}
}

bool building_t::add_cafeteria_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start) {
	if (room_has_stairs_or_elevator(room, zval, floor_ix)) return 0;
	float const floor_spacing(get_window_vspace());
	cube_t const room_bounds(get_walkable_room_bounds(room));
	vector2d const room_sz(room_bounds.get_size_xy());
	if (min(room_sz.x, room_sz.y) < 3.0*floor_spacing || max(room_sz.x, room_sz.y) < 3.5*floor_spacing) return 0; // too small to be a cafeteria
	if (!has_tile_floor()) {zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_LGTILE);}
	bool const dim(room_sz.y < room_sz.x); // short dim
	unsigned const num_cols((room_sz[dim] > 3.5*floor_spacing) ? 2 : 1);
	float const clearance(get_min_front_clearance_inc_people()), col_space(1.25*clearance), trim_thick(get_trim_thickness());
	cube_t place_area(room_bounds);
	place_area.expand_by_xy(-1.1*clearance);
	float const row_span(place_area.get_sz_dim(!dim)), col_span(place_area.get_sz_dim(dim)), col_start(place_area.d[dim][0]);
	float const table_height(0.3*floor_spacing), table_width(0.4*floor_spacing), table_len((col_span + col_space)/num_cols - col_space);
	float row_spacing(2.5*table_width);
	unsigned const num_rows(max(1U, unsigned(row_span/row_spacing))); // floor
	row_spacing = row_span/num_rows;
	colorRGBA const &chair_color(mall_chair_colors[rgen.rand() % NUM_MALL_CHAIR_COLORS]);
	unsigned const tid_tag(rgen.rand() + 1); // sets table texture; make nonzero to flag as a textured surface table
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval));
	cube_t table;
	set_cube_zvals(table, zval, zval+table_height);
	vect_cube_t blockers;
	blockers.reserve(doorways.size());

	for (door_stack_t const &d : doorways) {
		float const width(d.get_width());
		blockers.push_back(d.get_true_bcube());
		blockers.back().expand_in_dim( d.dim, 1.5*width); // more space in front of door
		blockers.back().expand_in_dim(!d.dim, 1.0*width); // less space to the side
	}
	for (unsigned r = 0; r < num_rows; ++r) {
		float const row_pos(place_area.d[!dim][0] + (r + 0.5)*row_spacing);
		set_wall_width(table, row_pos, 0.5*table_width, !dim);
		float pos(col_start); // starting point

		for (unsigned c = 0; c < num_cols; ++c) {
			table.d[dim][0] = pos;
			table.d[dim][1] = pos + table_len;
			
			for (cube_t const &b : blockers) {
				if (!b.intersects(table)) continue;
				if      (b.d[dim][0] <= table.d[dim][0]) {table.d[dim][0] = b.d[dim][1] + trim_thick;} // clip lo
				else if (b.d[dim][1] >= table.d[dim][1]) {table.d[dim][1] = b.d[dim][0] - trim_thick;} // clip hi
			}
			if (table.get_sz_dim(dim) > table_width) {add_mall_table_with_chairs(rgen, table, place_area, chair_color, room_id, tot_light_amt, dim, tid_tag, blockers);}
			pos += table_len + col_space;
		} // for c
	} // for r
	add_corner_trashcans  (rgen, room, zval, room_id, tot_light_amt, rgen.rand_bool(), 1); // random dim, both_ends=1
	add_clock_to_room_wall(rgen, room, zval, room_id, tot_light_amt, objs_start);
	add_door_sign("Cafeteria", room, zval, room_id);
	return 1;
}

bool building_t::add_library_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement, bool add_tables) {
	if (room.is_hallway || room.is_sec_bldg) return 0; // these can't be libraries
	bool const is_lg((is_prison() || is_school()));
	unsigned const num_bookcases(is_lg ? 16 : 8);

	for (unsigned n = 0; n < num_bookcases; ++n) { // place bookcases
		if (add_bookcase_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, is_basement)) {
			if (n == 0 && is_school()) { // add a row of computers on a long table after the first bookcase
				float const window_vspacing(get_window_vspace()), table_height(rgen.rand_uniform(0.4, 0.42)*window_vspacing);
				vector3d const sz_scale(rgen.rand_uniform(0.75, 0.8), rgen.rand_uniform(2.8, 3.2), 1.0); // depth, width, height
				cube_t const place_area(get_room_bounds_inside_trim(room));
				vect_room_object_t &objs(interior->room_geom->objs);
				unsigned const table_obj_ix(objs.size());
				if (!place_obj_along_wall(TYPE_TABLE, room, table_height, sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 1)) return 0;
				room_object_t const table(objs[table_obj_ix]); // deep copy to avoid reference invalidation
				unsigned const num_computers(3);
				float const comp_spacing(table.get_sz_dim(!table.dim)/num_computers);

				for (unsigned n = 0; n < num_computers; ++n) {
					cube_t sub_table(table);
					sub_table.d[!table.dim][0] += n*comp_spacing;
					sub_table.d[!table.dim][1]  = sub_table.d[!table.dim][0] + comp_spacing;
					add_computer_to_desk(sub_table, table_obj_ix, table.dim, !table.dir, rgen, room_id, tot_light_amt, 0.5); // sz_scale=0.5
				}
			}
		}
		else { // failed to add
			if (n == 0) return 0; // can't add a single bookcase
			break;
		}
	} // for n
	if (add_tables) {fill_room_with_tables_and_chairs(rgen, room, zval, room_id, tot_light_amt, objs_start, 0, 4);} // add tables; plastic_tc=0, max_books=4
	if (!is_house ) {add_door_sign_remove_existing("Library", room, zval, room_id, objs_start);} // add office building library sign
	return 1;
}

bool building_t::fill_room_with_tables_and_chairs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt,
	unsigned objs_start, bool plastic_tc, unsigned max_books, unsigned max_num_xy)
{
	float const vspace(get_window_vspace()), clearance(get_min_front_clearance_inc_people());
	float const table_spacing(0.5*vspace + 2.0*clearance); // placed tables will be rectangular, but we use square spacing
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by_xy(-(0.1*table_spacing + 0.05*min(room.dx(), room.dy()))); // add extra space at room walls
	vector2d const place_sz(place_area.get_size_xy());
	unsigned nx(place_sz.x/table_spacing), ny(place_sz.y/table_spacing);
	if (nx < 1 || ny < 1) return 0; // not enough space for any tables
	if (max_num_xy > 0) {min_eq(nx, max_num_xy); min_eq(ny, max_num_xy);}
	float const xspace(place_sz.x/nx), yspace(place_sz.y/ny);
	colorRGBA const &chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const tc_start(objs.size());
	unsigned num_added(0);
	vect_cube_t blockers;

	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (i->no_coll() || !i->intersects(place_area)) continue;
		blockers.push_back(*i);
		if (i->type == TYPE_BCASE) {blockers.back().expand_in_dim(i->dim, 0.5*clearance);} // add clearance in front of bookcases
	}
	for (unsigned y = 0; y < ny; ++y) {
		for (unsigned x = 0; x < nx; ++x) {
			point const center((place_area.x1() + (x + 0.5)*xspace), (place_area.y1() + (y + 0.5)*yspace), zval);
			num_added += add_table_and_chairs(rgen, room, blockers, room_id, center, chair_color, 0.0, tot_light_amt, 4, 0, plastic_tc); // no offset, 4 chairs, short table
		}
	} // for y
	if (max_books > 0) { // add books on tables
		vector<room_object_t> tables;

		for (auto i = objs.begin()+tc_start; i != objs.end(); ++i) {
			if (i->type == TYPE_TABLE) {tables.push_back(*i);}
		}
		for (room_object_t const &table : tables) {
			unsigned const pp_start(objs.size()), num_books(rgen.rand() % (max_books+1));
			unsigned const flags((is_school() && rgen.rand_bool()) ? RO_FLAG_USED : 0); // flag as school book half the time
			float const shift_amt((table.shape == SHAPE_CYLIN) ? 0.25 : 0.35); // less shift for cylindrical tables to avoid overlapping the edges
			for (unsigned n = 0; n < num_books; ++n) {place_book_on_obj(rgen, table, room_id, tot_light_amt, pp_start, 0, flags, 1, shift_amt);} // skip_if_overlaps=1
		}
	}
	// other objects? deck of cards, pack of cigarettes, or ashtray for prison lounge tables?
	return (num_added > 0);
}

