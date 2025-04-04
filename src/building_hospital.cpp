// 3D World - Hospital Buildings
// by Frank Gennari 2/26/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

extern object_model_loader_t building_obj_model_loader;


void clip_wall_to_ceil_floor(cube_t &wall, float fc_thick);
colorRGBA get_pastic_chair_color(colorRGBA const &color);

bool can_create_hospital_room() {return building_obj_model_loader.is_model_valid(OBJ_MODEL_HOSP_BED);}

void building_t::add_hospital_bathrooms(unsigned rooms_start, rand_gen_t &rgen) {
	assert(interior);
	if (!can_create_hospital_room()) return;
	vector<room_t> &rooms(interior->rooms);
	unsigned const rooms_end(rooms.size());
	assert(rooms_start <= rooms_end);
	rand_gen_t local_rgen(rgen); // copy so as not to disrupt rgen
	for (unsigned r = rooms_start; r < rooms_end; ++r) {maybe_create_nested_bathroom(rooms[r], local_rgen);} // can't iterate because rooms are added to
}

bool building_t::maybe_create_nested_bathroom(room_t &room, rand_gen_t &rgen) { // for hospital rooms
	if (is_room_windowless         (room)) return 0; // windowless room can't be a hospital bedroom, but it can be a bathroom, storage, office, etc.
	if (check_skylight_intersection(room)) return 0; // unlikely, but not handled here
	float const floor_spacing(get_window_vspace()), rdx(room.dx()), rdy(room.dy());
	if (min(rdx, rdy) < 2.0*floor_spacing || max(rdx, rdy) < 2.5*floor_spacing) return 0; // too small
	float const door_width(get_doorway_width()), door_hwidth(0.5*door_width), wall_thick(get_wall_thickness()), wall_hthick(0.5*wall_thick);
	float const min_sz(2.0*door_width), max_sz(4.0*door_width), rzc(room.zc());
	bool const pref_dx(rgen.rand_bool()), pref_dy(rgen.rand_bool());

	// choose a valid dir in each dim, not along an exterior wall
	for (unsigned DX = 0; DX < 2; ++DX) {
		bool const dx(bool(DX) ^ pref_dx);
		if (classify_room_wall(room, rzc, 0, dx, 0) == ROOM_WALL_EXT) continue; // skip if exterior wall

		for (unsigned DY = 0; DY < 2; ++DY) {
			bool const dy(bool(DY) ^ pref_dy);
			if (classify_room_wall(room, rzc, 1, dy, 0) == ROOM_WALL_EXT) continue; // skip if exterior wall
			cube_t bathroom(room);

			for (unsigned d = 0; d < 2; ++d) {
				bool const dir(d ? dy : dx);
				bathroom.d[d][!dir] = room.d[d][dir] + (dir ? -1.0 : 1.0)*max(min_sz, min(max_sz, 0.35f*(d ? rdy : rdx)));
			}
			assert(bathroom.is_strictly_normalized());
			if (is_cube_close_to_doorway(bathroom, room, 0.0, 1)) continue; // inc_open=1
			// Note: stairs have not yet been placed, so we don't need to avoid them
			room_t orig_room(room); // copy as room reference may be invalidated below
			room.copy_from(bathroom); // first (contained) room must be the bathroom
			calc_room_ext_sides(room);
			room.assign_all_to(RTYPE_BATH);
			room.set_is_nested();
			orig_room.set_has_subroom();
			interior->rooms.push_back(orig_room); // re-add full room; invalidates room reference
			// add wall sections and door
			float const fc_thick(get_fc_thickness());
			bool door_dim(0);
			cube_t walls[2] = {bathroom, bathroom};
			vect_door_stack_t doorways;
			get_doorways_for_room(orig_room, orig_room.zc(), doorways, 1); // get interior doors; all_floors=1
			
			for (unsigned d = 0; d < 2; ++d) { // wall dim
				cube_t &wall(walls[d]);
				clip_wall_to_ceil_floor(wall, fc_thick);
				set_wall_width(wall, bathroom.d[d][!(d ? dy : dx)], wall_hthick, d);
				if (d == 0) {wall.d[1][!dy] += (dy ? -1.0 : 1.0)*wall_hthick;} // expand to overlap the corner in X
				else        {wall.d[0][!dx] += (dx ? -1.0 : 1.0)*wall_hthick;} // shrink to remove overlap in Y
			}
			if (doorways.empty()) {door_dim = rgen.rand_bool(); assert(0);} // shouldn't happen
			else { // select wall side closer to room doorway
				point const door_center(doorways.front().get_cube_center());
				door_dim = (p2p_dist_xy_sq(door_center, walls[1].get_cube_center()) < p2p_dist_xy_sq(door_center, walls[0].get_cube_center()));
			}
			for (unsigned d = 0; d < 2; ++d) { // wall dim
				cube_t &wall(walls[d]);

				if (bool(d) == door_dim) { // add door in this dim
					float const door_center(bathroom.get_center_dim(!d));
					// opens into bathroom; keep_high_side=0, is_bathroom=0 (not always closed), make_unlocked=1, make_closed=0
					insert_door_in_wall_and_add_seg(wall, (door_center - door_hwidth), (door_center + door_hwidth), !d, (d ? dy : dx), 0, 0, 1, 0);
				}
				interior->walls[d].push_back(wall);
			} // for d
			return 1; // success/done
		} // for DY
	} // for DX
	return 0; // no valid location found
}

bool building_t::add_hospital_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	float tot_light_amt, unsigned objs_start, int &nested_room_ix)
{
	if (!can_create_hospital_room()) return 0; // no hospital bed
	float const floor_spacing(get_window_vspace()), clerance(1.25*get_min_front_clearance_inc_people()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds), bathroom;
	place_area.expand_by(-1.0*wall_thickness); // add extra padding, since bed models are slightly different sizes
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const beds_start(objs.size());
	unsigned const max_beds(max(1U, unsigned(0.25*room_bounds.get_area_xy()/(floor_spacing*floor_spacing))));
	unsigned num_beds(0), pref_orient(4);

	if (room.has_subroom()) {
		// first, determine if this room has a nested bathroom, and add a blocker over the entire room;
		// bathroom objects should have been placed already, so the blocker won't cause problems, and we can stop iterating at the current room
		if (nested_room_ix < 0) { // not yet found; find for the first floor
			for (unsigned r = 0; r < room_id; ++r) {
				if (room.contains_cube(get_room(r))) {nested_room_ix = r; break;} // there should only be one
			}
		}
		if (nested_room_ix >= 0) { // found
			bathroom = get_room(nested_room_ix);
			cube_t blocker(bathroom);
			blocker.expand_by_xy(0.5*wall_thickness); // include room walls
			objs.emplace_back(blocker, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS);
		}
	}
	for (unsigned n = 0; n < max_beds; ++n) {
		unsigned const bed_ix(objs.size());
		if (!place_model_along_wall(OBJ_MODEL_HOSP_BED, TYPE_HOSP_BED, room, 0.42, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.5, pref_orient)) continue;
		assert(bed_ix < objs.size());
		room_object_t &bed(objs[bed_ix]);
		assert(bed.type == TYPE_HOSP_BED);
		bed.item_flags = (7U*mat_ix + 11U*room.part_id + interior->rooms.size()); // select a random sub-model per building part
		room_object_t &blocker(objs.back());
		assert(blocker.type == TYPE_BLOCKER);
		blocker.expand_in_dim(!blocker.dim, clerance); // add extra clearance to the sides of the bed
		blocker.intersect_with_cube(room); // but don't go outside the room if near a wall
		pref_orient = 2*bed.dim + bed.dir; // try to place later beds along the same wall
		++num_beds;
	} // for n
	if (num_beds == 0) return 0;
	bool const add_tall_table(rgen.rand_bool());
	bool const add_curtains(num_beds > 1 && building_obj_model_loader.is_model_valid(OBJ_MODEL_HOSP_CURT));
	unsigned const beds_end(objs.size());
	point const table_pos(room.xc(), room.yc(), zval); // approximate
	vect_cube_t blockers;
	vect_cube_with_ix_t curtains;

	for (auto i = objs.begin()+beds_start; i != objs.end(); ++i) {
		if (i->type == TYPE_BLOCKER) {blockers.push_back(*i);}

		if (add_curtains && i->type == TYPE_HOSP_BED) { // place curtains between adjacent beds
			for (auto j = i+1; j != objs.end(); ++j) {
				if (j->type != TYPE_HOSP_BED || j->dim != i->dim || j->dir != i->dir) continue; // not a bed along the same wall
				vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_HOSP_CURT)); // D, W, H
				float const fc_gap(get_floor_ceil_gap()), height(0.9*fc_gap);
				cube_t merged(*i), curtain;
				merged.union_with_cube(*j);
				bool is_blocked(0);

				for (auto k = objs.begin()+beds_start; k != objs.end(); ++k) {
					if (k != i && k != j && k->type == TYPE_HOSP_BED && k->intersects(merged)) {is_blocked = 1; break;}
				}
				if (is_blocked) continue; // another bed is between i and j
				curtain.z2() = zval + fc_gap; // ceiling
				curtain.z1() = curtain.z2() - height; // hanging down near the floor
				set_wall_width(curtain, merged.get_center_dim(!i->dim), 0.5*height*sz.x/sz.z, !i->dim); // between the beds
				set_wall_width(curtain, merged.get_center_dim( i->dim), 0.5*height*sz.y/sz.z,  i->dim);
				curtains.emplace_back(curtain, !i->dim); // don't invalidate references
			} // for j
		}
	} // for i
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_TV)) { // add TVs on walls opposite beds
		for (unsigned i = beds_start; i != beds_end; ++i) {
			room_object_t const& bed(objs[i]);
			if (bed.type != TYPE_HOSP_BED) continue;
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
			float const tv_height(0.35*floor_spacing*rgen.rand_uniform(1.0, 1.2)), tv_hwidth(0.5*tv_height*sz.y/sz.z), tv_depth(tv_height*sz.x/sz.z);
			bool const dim(bed.dim), dir(bed.dir);
			cube_t tv;
			tv.z1() = zval    + 0.45*floor_spacing;
			tv.z2() = tv.z1() + tv_height;
			set_wall_width(tv, bed.get_center_dim(!dim), tv_hwidth, !dim);
			float wall_pos(room_bounds.d[dim][dir]); // start at opposite room wall

			if (!bathroom.is_all_zeros()) { // check bathoom wall
				if (tv.d[!dim][0] > bathroom.d[!dim][0] && tv.d[!dim][1] < bathroom.d[!dim][1]) {
					wall_pos = bathroom.d[dim][!dir] + (dir ? -1.0 : 1.0)*0.5*wall_thickness;
				}
			}
			tv.d[dim][ dir] = wall_pos;
			tv.d[dim][!dir] = tv.d[dim][dir] + (dir ? -1.0 : 1.0)*tv_depth;
			cube_t test_cube(tv);
			test_cube.d[dim][!dir] = bed.d[dim][!dir]; // extend to the bed; should this ignore open doors?
			
			if (!overlaps_obj_or_placement_blocked(test_cube, room, objs_start) && !check_if_against_window(tv, room, dim, dir)) { // valid placement
				add_tv_to_wall(tv, room_id, tot_light_amt, dim, !dir, 0, 2); // use_monitor_image=0, on_off=2 (random); invalidates bed reference
			}
		} // for i
	}
	for (cube_with_ix_t const &c : curtains) {objs.emplace_back(c, TYPE_HOSP_CURT, room_id, c.ix, rgen.rand_bool(), 0);}
	unsigned const table_ix(objs.size());

	for (unsigned n = 0; n < 4; ++n) { // 4 attempts to place a table and chair
		if (add_table_and_chairs(rgen, room, blockers, room_id, table_pos, WHITE, 0.25, tot_light_amt, 1, add_tall_table) == 0) continue; // 1 chair

		if (rgen.rand_bool()) { // maybe place a phone on the table using a random dim/dir; should it face the chair?
			assert(table_ix < objs.size());
			room_object_t &table(objs[table_ix]);
			assert(table.type == TYPE_TABLE);
			if (place_phone_on_obj(rgen, table, room_id, tot_light_amt, rgen.rand_bool(), rgen.rand_bool())) {table.flags |= RO_FLAG_ADJ_TOP;}
		}
		break; // done/success
	} // for n
	add_numbered_door_sign("Room ", room, zval, room_id, floor_ix);
	return 1;
}

bool building_t::add_waiting_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	float tot_light_amt, unsigned objs_start, colorRGBA const &chair_color)
{
	unsigned const counts[4] = {1, 1, 2, 2}; // 1-2
	add_couches_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, counts);
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t const room_bounds(get_walkable_room_bounds(room));

	if (rgen.rand_bool()) { // add some reading material
		add_bookcase_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, 0); // is_basement=0
	}
	else { // add some snacks
		cube_t vm_area(room_bounds);
		vm_area.expand_by_xy(-0.25*wall_thickness);
		add_vending_machine(rgen, room, zval, room_id, tot_light_amt, objs_start, vm_area);
	}
	// can this block a door sign in a room with ring hallways?
	// there's not much we can do about that here, since the sign is part of another room that may either be before or after this one
	add_wall_tv(rgen, room, zval, room_id, tot_light_amt, objs_start);
	unsigned const num_chairs = 20; // up to this many, whatever we can fit
	bool const is_plastic(rgen.rand_bool());
	float const chair_hscale(is_plastic ? 0.44 : 0.4), chair_height(chair_hscale*floor_spacing), chair_xy_scale(0.2/chair_hscale);
	vector3d const chair_sz(chair_xy_scale, chair_xy_scale, 1.0);
	colorRGBA const ccolor(is_plastic ? get_pastic_chair_color(chair_color) : chair_color);
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t chair_place_area(room_bounds);
	chair_place_area.expand_by_xy(-wall_thickness);

	for (unsigned n = 0; n < num_chairs; ++n) {
		unsigned const chair_obj_ix(objs.size());
		if (!place_obj_along_wall(TYPE_CHAIR, room, chair_height, chair_sz, rgen, zval, room_id, tot_light_amt, chair_place_area, objs_start,
			1.0, 0, 4, 0, ccolor, 0, SHAPE_CUBE, 0.35*chair_height)) break; // end when failed to place
		assert(chair_obj_ix < objs.size());
		if (is_plastic) {objs[chair_obj_ix].item_flags = 1;} // flag as plastic
	}
	// add some more chairs to the center of the room
	// TODO
	add_clock_to_room_wall(rgen, room, zval, room_id, tot_light_amt, objs_start);
	unsigned const num_plants(1 + rgen.rand_bool()); // add 1-2 plants
	add_plants_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, num_plants);
	add_door_sign("Waiting Area", room, zval, room_id);

	// make sure doors start open and unlocked
	for (door_stack_t const &ds : interior->door_stacks) {
		if (!ds.is_connected_to_room(room_id)) continue;
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() > zval || door.z2() < zval) continue; // wrong floor
			door.locked   = 0;
			door.open     = 1;
			door.open_amt = 1.0;
		}
	} // for ds
	return 1;
}

bool building_t::add_exam_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	float tot_light_amt, unsigned objs_start, colorRGBA const &chair_color)
{
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-1.0*get_wall_thickness()); // add extra padding, since bed models are slightly different sizes
	if (!place_model_along_wall(OBJ_MODEL_HOSP_BED, TYPE_HOSP_BED, room, 0.42, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.5)) return 0;
	// TODO: should be a small desk, equipment, etc.
	add_desk_to_room(rgen, room, vect_cube_t(), chair_color, zval, room_id, tot_light_amt, objs_start, 0, 0, 0, 1); // force_computer=1
	//place_phone_on_obj(rgen, place_on, room_id, tot_light_amt, dim, dir); // place on the desk?
	place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6);
	// should be a short rotating stool?
	place_model_along_wall(OBJ_MODEL_BAR_STOOL, TYPE_BAR_STOOL, room, 0.4, rgen, zval, room_id, tot_light_amt, place_area, objs_start);
	add_filing_cabinet_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);
	// the room should already have a trashcan, but we can add a second one for medical waste
	add_trashcan_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, 0); // check_last_obj=0
	add_numbered_door_sign("Exam ", room, zval, room_id, floor_ix);
	return 1;
}

bool building_t::add_operating_room_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id,
	unsigned floor_ix, float &tot_light_amt, unsigned objs_start, unsigned lights_start)
{
	// TODO: operating table, equipment (with TYPE_VALVE/TYPE_GAUGE)
	assert(lights_start <= objs_start);
	bool const long_dim(room.dx() < room.dy());
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_HOSP_BED)) return 0; // don't have a model of this type
	// brighter lights in OR
	room.light_intensity *= 2.0;
	tot_light_amt        *= 2.0;
	float const floor_spacing(get_window_vspace()), height(0.42*floor_spacing);
	vect_room_object_t &objs(interior->room_geom->objs);
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_HOSP_BED)); // D, W, H
	point table_center(room.xc(), room.yc(), zval); // start at the center of the room
	// determine table_dir based on door pos if in same dim
	bool table_dir(rgen.rand_bool());
	vect_door_stack_t doorways;
	get_doorways_for_room(room, zval, doorways); // get interior doors

	if (!doorways.empty()) { // should always be true
		door_stack_t const door(doorways.front());
		bool const door_dir(door.get_center_dim(door.dim) < table_center[door.dim]);
		if (door.dim == long_dim) {table_dir = door_dir;} // feet facing the door
		table_center[door.dim] += (door_dir ? 1.0 : -1.0)*0.5*door.get_width(); // move table back from the door
	}
	// place table if valid
	cube_t op_table(table_center);
	op_table.z2() += height;
	op_table.expand_in_dim( long_dim, 0.5*height*(sz.x/sz.z));
	op_table.expand_in_dim(!long_dim, 0.5*height*(sz.y/sz.z));
	if (is_obj_placement_blocked(op_table, room, 1)) return 0; // inc_open_doors=1; no need to check for other objects since this is first
	objs.emplace_back(op_table, TYPE_HOSP_BED, room_id, long_dim, table_dir, 0, tot_light_amt);
	// add 1-2 machines
	cube_t place_area(get_walkable_room_bounds(room));
	unsigned const num_machines(1 + rgen.rand_bool());

	for (unsigned n = 0; n < num_machines; ++n) {
		float const machine_height(rgen.rand_uniform(0.6, 0.8)*floor_spacing), xy_scale(rgen.rand_uniform(0.4, 0.7));
		place_obj_along_wall(TYPE_MACHINE, room, machine_height, vector3d(xy_scale, xy_scale, 1.0), rgen, zval, room_id, tot_light_amt,
			place_area, objs_start, 0.0, 0, 4, 0, WHITE, 1); // not_at_window=1
	}
	add_clock_to_room_wall(rgen, room, zval, room_id, tot_light_amt, objs_start);
	add_numbered_door_sign("OR ", room, zval, room_id, floor_ix);

	// make ceiling lights larger and round
	for (auto i = objs.begin()+lights_start; i != objs.begin()+objs_start; ++i) {
		if (i->type != TYPE_LIGHT) continue;
		float const dx(i->dx()), dy(i->dy()); // make square
		if      (dx < dy) {i->expand_in_x(0.5*(dy - dx));}
		else if (dy < dx) {i->expand_in_y(0.5*(dx - dy));}
		i->shape = SHAPE_CYLIN;
	} // for i
	return 1;
}

