// 3D World - Hospital Buildings
// by Frank Gennari 2/26/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

extern object_model_loader_t building_obj_model_loader;


void clip_wall_to_ceil_floor(cube_t &wall, float fc_thick);
colorRGBA get_pastic_chair_color(colorRGBA const &color);
colorRGBA get_couch_color(rand_gen_t &rgen);

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

bool building_t::get_hospital_room_bathroom(room_t const &room, unsigned room_id, int &nested_room_ix, cube_t &bathroom) const {
	if (!room.has_subroom()) return 0;
	assert(has_room_geom());

	// first, determine if this room has a nested bathroom, and add a blocker over the entire room;
	// bathroom objects should have been placed already, so the blocker won't cause problems, and we can stop iterating at the current room
	if (nested_room_ix < 0) { // not yet found; find for the first floor
		for (unsigned r = 0; r < room_id; ++r) {
			if (room.contains_cube(get_room(r))) {nested_room_ix = r; break;} // there should only be one
		}
	}
	if (nested_room_ix < 0) return 0; // not found
	bathroom = get_room(nested_room_ix);
	return 1;
}

bool building_t::add_hospital_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	float tot_light_amt, unsigned objs_start, int &nested_room_ix)
{
	if (!can_create_hospital_room()) return 0; // no hospital bed
	float const floor_spacing(get_window_vspace()), clerance(1.25*get_min_front_clearance_inc_people()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds), bathroom;
	place_area.expand_by(-1.0*wall_thickness); // add extra padding, since bed models are slightly different sizes
	vect_room_object_t &objs(interior->room_geom->objs);
	bool const has_bathroom(get_hospital_room_bathroom(room, room_id, nested_room_ix, bathroom));
	unsigned const beds_start(objs.size());
	unsigned const max_beds(max(1U, unsigned(0.25*room_bounds.get_area_xy()/(floor_spacing*floor_spacing))));
	unsigned num_beds(0), pref_orient(4);

	if (has_bathroom) { // add blocker for bathroom
		objs.emplace_back(bathroom, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS);
		objs.back().expand_by_xy(0.5*wall_thickness); // include room walls
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
	// add curtains between beds
	vect_cube_t blockers, tv_blockers;
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
				if (is_obj_placement_blocked(curtain, room, 1, 0)) continue; // in particular, check for doors
				curtains.emplace_back(curtain, !i->dim); // don't invalidate references
				curtain.d[i->dim][!i->dir] = room.d[i->dim][!i->dir]; // extend to the wall
				tv_blockers.push_back(curtain);
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
			float const center(max((place_area.d[!dim][0] + tv_hwidth), min((place_area.d[!dim][1] - tv_hwidth), bed.get_center_dim(!dim)))); // clamp to place area
			cube_t tv;
			tv.z1() = zval    + 0.45*floor_spacing;
			tv.z2() = tv.z1() + tv_height;
			set_wall_width(tv, center, tv_hwidth, !dim);
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
			if (has_bcube_int(tv, tv_blockers) || overlaps_obj_or_placement_blocked(test_cube, room, objs_start) || check_if_against_window(tv, room, dim, dir)) continue;
			add_tv_to_wall(tv, room_id, tot_light_amt, dim, !dir, 0, 2); // use_monitor_image=0, on_off=2 (random); invalidates bed reference
			// add a blocker in front of the TV to avoid placing objects such as plants there
			cube_t blocker(tv);
			blocker.d[dim][!dir] += (dir ? -1.0 : 1.0)*floor_spacing;
			objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, !dir, RO_FLAG_INVIS);
		} // for i
	}
	for (cube_with_ix_t const &c : curtains) {objs.emplace_back(c, TYPE_HOSP_CURT, room_id, c.ix, rgen.rand_bool(), 0);}
	unsigned const table_ix(objs.size());

	for (unsigned n = 0; n < 8; ++n) { // 8 attempts to place a table and chair
		if (add_table_and_chairs(rgen, room, blockers, room_id, table_pos, WHITE, 0.25, tot_light_amt, 1, add_tall_table) == 0) { // 1 chair
			rgen.rand_mix(); // needed to get different rand values
			continue;
		}
		vect_cube_t avoid;

		if (rgen.rand_bool()) { // maybe place a phone on the table using a random dim/dir; should it face the chair?
			assert(table_ix < objs.size());
			room_object_t &table(objs[table_ix]);
			assert(table.type == TYPE_TABLE);

			if (place_phone_on_obj(rgen, table, room_id, tot_light_amt, rgen.rand_bool(), rgen.rand_bool())) {
				objs[table_ix].flags |= RO_FLAG_ADJ_TOP;
				avoid.push_back(objs.back());
			}
		}
		place_plant_on_obj(rgen, objs[table_ix], room_id, tot_light_amt, 0.7, avoid);
		break; // done/success
	} // for n
	if (1) { // maybe add chairs
		unsigned const num_chairs(rgen.rand() % (num_beds+1)); // up to 1 per bed
		bool const is_plastic(rgen.rand_bool());
		colorRGBA const &chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
		place_chairs_along_walls(rgen, room, zval, room_id, tot_light_amt, objs_start, chair_color, is_plastic, num_chairs);
	}
	bool placed_couch(0);

	if (num_beds > 1 && rgen.rand_bool()) { // maybe add a couch
		colorRGBA const color(get_couch_color(rgen));
		placed_couch = place_model_along_wall(OBJ_MODEL_COUCH, TYPE_COUCH, room, 0.40, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 4, 0, color);
	}
	if (!placed_couch && rgen.rand_bool()) {
		place_model_along_wall(OBJ_MODEL_RCHAIR, TYPE_RCHAIR, room, 0.5, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0);
	}
	add_plants_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, 1); // add 1 plant
	add_numbered_door_sign("Room ", room, zval, room_id, floor_ix);
	return 1;
}

bool building_t::add_waiting_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	float tot_light_amt, unsigned objs_start, int &nested_room_ix)
{
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t const room_bounds(get_walkable_room_bounds(room));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t bathroom;

	// add blocker for bathroom, in case this room has stairs and was assigned to a waiting room
	if (get_hospital_room_bathroom(room, room_id, nested_room_ix, bathroom)) {
		objs.emplace_back(bathroom, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS);
		objs.back().expand_by_xy(0.5*wall_thickness); // include room walls
	}
	unsigned const counts[4] = {1, 1, 2, 2}; // 1-2
	add_couches_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, counts);

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
	// add chairs along walls
	unsigned const num_chairs = 20; // up to this many, whatever we can fit
	bool const is_plastic(1); // always plastic
	colorRGBA const &chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
	place_chairs_along_walls(rgen, room, zval, room_id, tot_light_amt, objs_start, chair_color, is_plastic, num_chairs);
	// maybe add some more chairs to the center of the room
	cube_t chair_place_area(room_bounds);
	chair_place_area.expand_by_xy(-wall_thickness);
	vector2d const room_sz(chair_place_area.get_size_xy());
	bool const long_dim(room_sz.x < room_sz.y);
	float const length(room_sz[long_dim]), width(room_sz[!long_dim]);
	float const chair_hscale(is_plastic ? 0.44 : 0.4), chair_height(chair_hscale*floor_spacing), chair_gap(0.35*chair_height);
	float const chair_width(0.2*chair_height/chair_hscale), chair_hwidth(0.5*chair_width), min_chair_spacing(chair_width + chair_gap);
	float const clearance(get_min_front_clearance_inc_people()), chair_back_space(3.0*wall_thickness);
	float const path_clearance(2.0*clearance + 2.0*chair_width), room_min_sz(2.0*path_clearance), chair_place_len(length - room_min_sz);

	if (chair_place_len > 2.0*min_chair_spacing && width > (room_min_sz + chair_width + wall_thickness)) { // space for at least two chairs
		float const centerline(room.get_center_dim(!long_dim));
		unsigned const num_chairs(chair_place_len/min_chair_spacing); // take the floor
		float const chair_spacing(chair_place_len/num_chairs), chair_start(chair_place_area.d[long_dim][0] + path_clearance + 0.5*chair_spacing);
		// add a wall between the rows of chairs
		cube_t wall;
		set_cube_zvals(wall, zval, zval+0.8*chair_height);
		set_wall_width(wall, centerline, 0.45*wall_thickness, !long_dim); // set width
		set_wall_width(wall, room.get_center_dim(long_dim), 0.5*chair_place_len, long_dim); // set length
		objs.emplace_back(wall, TYPE_STAIR_WALL, room_id, !long_dim, 0, RO_FLAG_ADJ_TOP, tot_light_amt, SHAPE_CUBE); // draw top
		// add chairs
		colorRGBA const ccolor(is_plastic ? get_pastic_chair_color(chair_color) : chair_color);
		cube_t chair0;
		set_cube_zvals(chair0, zval, zval+chair_height);

		for (unsigned d = 0; d < 2; ++d) { // each side
			float const dsign(d ? 1.0 : -1.0), back_pos(centerline + dsign*0.5*chair_back_space);
			chair0.d[!long_dim][!d] = back_pos;
			chair0.d[!long_dim][ d] = back_pos + dsign*chair_width;

			for (unsigned n = 0; n < num_chairs; ++n) {
				set_wall_width(chair0, (chair_start + n*chair_spacing), chair_hwidth, long_dim);
				cube_t chair(chair0);
				for (unsigned e = 0; e < 2; ++e) {chair.translate_dim(e, 0.5*wall_thickness*rgen.signed_rand_float());} // add a bit of misalignment
				cube_t chair_pad(chair);
				chair_pad.d[!long_dim][d] += dsign*clearance; // extend further outward
				if (overlaps_obj_or_placement_blocked(chair_pad, chair_place_area, objs_start)) continue; // bad placement
				objs.emplace_back(chair, TYPE_CHAIR, room_id, !long_dim, d, RO_FLAG_PLCOLL, tot_light_amt, SHAPE_CUBE, ccolor);
				if (is_plastic) {objs.back().item_flags = 1;} // flag as plastic
			} // for n
		} // for d
	}
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

void building_t::place_chairs_along_walls(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt,
	unsigned objs_start, colorRGBA const &chair_color, bool is_plastic, unsigned num)
{
	if (num == 0) return; // error?
	float const chair_hscale(is_plastic ? 0.44 : 0.4), chair_height(chair_hscale*get_window_vspace()), chair_xy_scale(0.2/chair_hscale), chair_gap(0.35*chair_height);
	vector3d const chair_sz(chair_xy_scale, chair_xy_scale, 1.0);
	colorRGBA const ccolor(is_plastic ? get_pastic_chair_color(chair_color) : chair_color);
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t chair_place_area(get_walkable_room_bounds(room));
	chair_place_area.expand_by_xy(-get_wall_thickness());

	for (unsigned n = 0; n < num; ++n) {
		unsigned const chair_obj_ix(objs.size());
		if (!place_obj_along_wall(TYPE_CHAIR, room, chair_height, chair_sz, rgen, zval, room_id, tot_light_amt, chair_place_area, objs_start,
			1.0, 0, 4, 0, ccolor, 0, SHAPE_CUBE, chair_gap, RO_FLAG_PLCOLL)) break; // end when failed to place
		assert(chair_obj_ix < objs.size());
		if (is_plastic) {objs[chair_obj_ix].item_flags = 1;} // flag as plastic
	}
}

bool building_t::add_exam_room_objs(rand_gen_t rgen, room_t &room, float zval, unsigned room_id, unsigned floor_ix, float tot_light_amt, unsigned objs_start) {
	cube_t const room_area(get_walkable_room_bounds(room));
	cube_t place_area(room_area);
	place_area.expand_by(-1.0*get_wall_thickness()); // add extra padding, since bed models are slightly different sizes
	if (!place_model_along_wall(OBJ_MODEL_HOSP_BED, TYPE_HOSP_BED, room, 0.42, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.5)) return 0;
	vect_room_object_t &objs(interior->room_geom->objs);
	// TODO: should be a small desk, equipment, etc.
	colorRGBA const &chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
	add_desk_to_room(rgen, room, vect_cube_t(), chair_color, zval, room_id, tot_light_amt, objs_start, 0, 0, 0, 1); // force_computer=1
	//place_phone_on_obj(rgen, place_on, room_id, tot_light_amt, dim, dir); // place on the desk?

	if (rgen.rand_bool()) { // add a simple sink
		place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6);
	}
	else { // add a vanity
		add_vanity_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);
	}
	place_chairs_along_walls(rgen, room, zval, room_id, tot_light_amt, objs_start, chair_color, 1, 1); // is_plastic=1, num_chairs=1
	// should be a short rotating stool?
	place_model_along_wall(OBJ_MODEL_BAR_STOOL, TYPE_BAR_STOOL, room, 0.4, rgen, zval, room_id, tot_light_amt, place_area, objs_start);
	add_filing_cabinet_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);

	if (1) { // add a small wall mounted computer monitor
		unsigned const tv_obj_ix(objs.size());
		float const z1(zval + 0.5*get_window_vspace());
		
		if (place_model_along_wall(OBJ_MODEL_TV, TYPE_MONITOR, room, 0.25, rgen, z1, room_id, tot_light_amt, room_area, objs_start, 4.0, 4, 1, BKGRAY, 1, RO_FLAG_HANGING)) {
			offset_hanging_tv(objs[tv_obj_ix]);
			objs[tv_obj_ix].obj_id |= 1; // make it turned off
		}
	}
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

