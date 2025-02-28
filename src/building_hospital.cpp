// 3D World - Hospital Buildings
// by Frank Gennari 2/26/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

extern object_model_loader_t building_obj_model_loader;


void clip_wall_to_ceil_floor(cube_t &wall, float fc_thick);

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

				if (d == door_dim) { // add door in this dim
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

bool building_t::add_hospital_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, int &nested_room_ix) {
	if (!can_create_hospital_room()) return 0; // no hospital bed
	float const floor_spacing(get_window_vspace()), clerance(1.25*get_min_front_clearance_inc_people()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
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
			cube_t blocker(get_room(nested_room_ix));
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
	return 1;
}

