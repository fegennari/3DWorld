// 3D World - Building Kitchens
// by Frank Gennari 12/13/2025

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t

extern object_model_loader_t building_obj_model_loader;

void get_stove_burner_locs(room_object_t const &stove, point locs[4]);
float get_cockroach_height_from_radius(float radius);
float get_plate_radius(rand_gen_t &rgen, cube_t const &place_on, float window_vspacing);
bool cube_map_reflect_active();
void add_stack_of_plates(cube_t const &place_area, float radius, unsigned room_id, float light_amt, unsigned flags,
	rand_gen_t &rgen, vect_cube_t &blockers, vect_room_object_t &objects);


void add_door_if_blocker(cube_t const &door, cube_t const &room, bool inc_open, bool dir, bool hinge_side, bool exterior, vect_cube_t &blockers) {
	bool const dim(door.dy() < door.dx()), edir(dim ^ dir ^ hinge_side ^ 1);
	float const width(door.get_sz_dim(!dim));
	cube_t door_exp(door);
	door_exp.expand_in_dim(dim, width); // width
	if (!door_exp.intersects(room)) return; // check against room before expanding along wall to exclude doors in adjacent rooms
	if (exterior) {door_exp.expand_in_dim(dim, width);} // add extra clearance in front
	door_exp.expand_in_dim(!dim, width*0.25); // min expand value
	if (inc_open) {door_exp.d[!dim][edir] += (edir ? 1.0 : -1.0)*0.75*width;} // expand the remainder of the door width in this dir
	blockers.push_back(door_exp);
}
// used for kitchen cabinets
int building_t::gather_room_placement_blockers(cube_t const &room, unsigned objs_start, vect_cube_t &blockers, bool inc_open_doors, bool ignore_chairs) const {
	assert(has_room_geom());
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(objs_start <= objs.size());
	blockers.clear();
	int table_blocker_ix(-1);

	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (ignore_chairs && i->type == TYPE_CHAIR) continue;
		
		if (!i->no_coll() && i->intersects(room)) {
			if (i->type == TYPE_TABLE) {table_blocker_ix = int(blockers.size());} // track which blocker is the table, for use with kitchen counters
			blockers.push_back(*i);
		}
	}
	for (auto const &door : doors) { // exterior doors
		add_door_if_blocker(door.get_bcube(), room, 0, 0, 0, 1, blockers);
	}
	for (door_stack_t const &ds : interior->door_stacks) { // interior doors
		add_door_if_blocker(ds, room, door_opens_inward(ds, room), ds.open_dir, ds.hinge_side, 0, blockers);
	}
	// Note: caller must call is_blocked_by_stairs_or_elevator() to handle stairs and elevators
	return table_blocker_ix;
}

void building_t::add_fridge_sticky_notes(rand_gen_t rgen, unsigned fridge_obj_ix, float zval, unsigned room_id, float tot_light_amt) {
	if (rgen.rand_float() < 0.25) return; // no sticky notes
	vect_room_object_t &objs(interior->room_geom->objs);
	room_object_t const fridge(objs[fridge_obj_ix]); // deep copy to avoid invaliding the reference
	room_object_t note(fridge);
	float const one_inch(get_one_inch()), fridge_width(fridge.get_width()), fridge_height(fridge.get_height());
	if (fridge_width < 4.0*4.0*one_inch || fridge_height < 6.0*6.0*one_inch) return; // shouldn't fail
	bool const dim(fridge.dim), dir(fridge.dir);
	unsigned const NUM_SIZES(7), NUM_COLORS(8);
	vector2d  const note_sizes [NUM_SIZES ] = {{1.5, 2}, {2, 2}, {3, 3}, {4, 2}, {4, 4}, {3, 5}, {4, 6}}; // standard sizes, in inches
	colorRGBA const note_colors[NUM_COLORS] = {LT_YELLOW, LT_YELLOW, LT_YELLOW, PINK, ORANGE, CYAN, MAGENTA, LT_GREEN};
	float const front_face(fridge.d[dim][dir] - (dir ? 1.0 : -1.0)*0.08*fridge.get_depth()); // exclude the handles depth
	note.d[dim][!dir] = front_face; // back on front of fridge
	note.d[dim][ dir] = front_face + (dir ? 1.0 : -1.0)*0.1*one_inch;
	note.type = TYPE_STICK_NOTE;
	unsigned const notes_start(objs.size()), num_notes((rgen.rand() % 16) + 1); // 1-16
	vector2d note_size;

	// what about fridge sides (not against walls)?
	for (unsigned n = 0; n < num_notes; ++n) {
		if (n == 0 || (rgen.rand()&3) == 0) { // change color or size
			note.color = note_colors[rgen.rand() % NUM_COLORS];
			note_size  = note_sizes [rgen.rand() % NUM_SIZES ];
		}
		float const note_hwidth(0.5*note_size.x*one_inch), note_hheight(0.5*note_size.y*one_inch);
		float left(fridge.d[!dim][0] + 1.5*note_hwidth), right(fridge.d[!dim][1] - 1.5*note_hwidth); // of placement ranges
		float const bot(fridge.z1() + 0.4*fridge_height + 1.5*note_hheight), top(fridge.z2() - 0.1*fridge_height - 1.5*note_hheight); // of placement ranges
		if (dim ^ dir) {left += 0.6*fridge_width;} else {right -= 0.6*fridge_width;} // on side of door opposite water dispenser, avoiding door handle

		for (unsigned N = 0; N < 4; ++N) { // 4 attempts to find a non-overlapping note
			set_wall_width(note, rgen.rand_uniform(left, right), note_hwidth, !dim);
			set_wall_width(note, rgen.rand_uniform(bot,  top  ), note_hheight, 2);
			if (has_bcube_int(note, objs, notes_start)) continue;
			objs.push_back(note);
			break;
		}
	} // for n
}

colorRGBA get_toaster_color(rand_gen_t &rgen) {
	if (cube_map_reflect_active()) return WHITE; // reflective metal
	return toaster_colors[rgen.rand()%NUM_TOASTER_COLORS]; // matte colored
}
bool building_t::add_kitchen_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt,
	unsigned objs_start, bool allow_adj_ext_door, light_ix_assign_t &light_ix_assign)
{
	// Note: table and chairs have already been placed
	bool const residential(is_residential());
	if (room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be kitchens
	if (!residential && rgen.rand_bool()) return 0; // some office buildings have kitchens, allow it half the time
	// if it has an external door then reject the room half the time; most houses don't have a front door to the kitchen
	if (is_room_adjacent_to_ext_door(room, zval, 1) && (!allow_adj_ext_door || rgen.rand_bool())) return 0; // front_door_only=1
	float const vspace(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall for appliances
	vect_room_object_t &objs(interior->room_geom->objs);
	
	if (residential && min(room.dx(), room.dy()) > 2.0*vspace) { // add a pantry if large, which is of type TYPE_CLOSET
		float const clearance(get_min_front_clearance_inc_people());
		unsigned closet_obj_id(0); // unused
		add_closet_to_room(rgen, room, zval, room_id, objs_start, RTYPE_KITCHEN, 0, clearance, closet_obj_id, light_ix_assign); // bed_obj_ix=0 (not set)
	}
	// place a fridge
	unsigned const fridge_obj_ix(objs.size());
	bool placed_obj(0);

	if (place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, room, 0.75, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.2, 4, 0, WHITE, 1)) { // not at window
		if (is_house) {add_fridge_sticky_notes(rgen, fridge_obj_ix, zval, room_id, tot_light_amt);}
		placed_obj = 1;
	}
	if (residential) { // try to place a stove
		unsigned const stove_ix(objs.size()); // can't use objs.back() because there's a blocker
		
		if (place_model_along_wall(OBJ_MODEL_STOVE, TYPE_STOVE, room, 0.46, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0)) {
			assert(stove_ix < objs.size());

			if (building_obj_model_loader.is_model_valid(OBJ_MODEL_HOOD)) { // add hood above the stove
				room_object_t const &stove(objs[stove_ix]);
				vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_HOOD)); // D, W, H
				float const width(stove.get_sz_dim(!stove.dim)), height(width*sz.z/sz.y), depth(width*sz.x/sz.y); // scale to the width of the stove
				float const ceiling_z(zval + get_floor_ceil_gap()), z_top(ceiling_z + get_fc_thickness()); // shift up a bit because it's too low
				cube_t hood(stove);
				set_cube_zvals(hood, z_top-height, z_top);
				hood.d[stove.dim][stove.dir] = stove.d[stove.dim][!stove.dir] + (stove.dir ? 1.0 : -1.0)*depth;

				if (!check_if_against_window(hood, room, stove.dim, !stove.dir)) { // don't place hoods against windows
					objs.emplace_back(hood, TYPE_HOOD, room_id, stove.dim, stove.dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, LT_GRAY);
				
					if (has_attic()) { // add rooftop vent above the hood; should be close to the edge of the house in many cases
						float const attic_floor_zval(get_attic_part().z2()), vent_radius(0.075*(width + depth)); // kitchen can be in a part below the one with the attic
						point const vent_bot_center(hood.xc(), hood.yc(), attic_floor_zval);
						add_attic_roof_vent(vent_bot_center, vent_radius, room_id, 1.0); // light_amt=1.0; room_id is for the kitchen because there's no attic room
					}
				}
			}
			if (!rgen.rand_bool()) { // maybe add a pan on one of the stove burners
				room_object_t const &stove(objs[stove_ix]);
				float const stove_height(stove.dz()), delta_z(0.018*stove_height);
				float const pan_radius(rgen.rand_uniform(0.075, 0.09)*stove_height), pan_height(rgen.rand_uniform(0.035, 0.045)*stove_height);
				point locs[4];
				get_stove_burner_locs(stove, locs);
				unsigned const burner_ix(rgen.rand() & 3); // 0-3
				point &loc(locs[burner_ix]);
				loc.z += delta_z;
				cube_t const pan(get_cube_height_radius(loc, pan_radius, pan_height));
				objs.emplace_back(pan, TYPE_PAN, room_id, stove.dim, stove.dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, DK_GRAY);
			}
			placed_obj = 1;
		}
	}	
	if (residential && placed_obj) { // if we have at least a fridge or stove, try to add countertops
		float const height(0.345*vspace), depth(0.74*height), floor_thickness(get_floor_thickness()), min_clearance(get_min_front_clearance_inc_people());
		float min_hwidth(0.6*height), front_clearance(max(0.6f*height, min_clearance));
		cube_t cabinet_area(room_bounds);
		cabinet_area.expand_by(-0.05*wall_thickness); // smaller gap than place_area; this is needed to prevent z-fighting with exterior walls
		if (min(cabinet_area.dx(), cabinet_area.dy()) < 4.0*min_hwidth) return placed_obj; // no space for cabinets, room is too small
		unsigned const counters_start(objs.size());
		cube_t c;
		set_cube_zvals(c, zval, zval+height);
		set_cube_zvals(cabinet_area, zval, (zval + vspace - floor_thickness));
		static vect_cube_t blockers;
		int const table_blocker_ix(gather_room_placement_blockers(cabinet_area, objs_start, blockers, 1, 1)); // inc_open_doors=1, ignore_chairs=1
		bool const have_toaster(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOASTER));
		bool const have_milk   (building_obj_model_loader.is_model_valid(OBJ_MODEL_MILK   ));
		vector3d const toaster_sz(have_toaster ? building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOASTER) : zero_vector); // L, D, H
		bool is_sink(1), placed_mwave(0), placed_toaster(0), placed_milk(0), had_counter(0);
		unsigned num_paper_towels(0);
		cube_t mwave, toaster;

		for (unsigned n = 0; n < (had_counter ? 50U : 60U); ++n) { // 50 attempts, plus an extra 10 if no counters were placed
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
			if (room.has_open_wall(dim, dir)) continue; // don't place against open walls
			bool const is_ext_wall(classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT); // assumes not in basement
			// only consider exterior walls in the first 20 attempts to prioritize these so that we don't have splits visible through windows; also places kitchen sinks near windows
			if (n < 20 && !is_ext_wall) continue;

			if (n > 30 && !had_counter) { // reduce min width and clearance for later iterations if no counters were placed
				min_hwidth      -= 0.01*height;
				front_clearance -= 0.01*height;
			}
			float const center(rgen.rand_uniform(cabinet_area.d[!dim][0]+min_hwidth, cabinet_area.d[!dim][1]-min_hwidth)); // random position
			float const dir_sign(dir ? -1.0 : 1.0), wall_pos(cabinet_area.d[dim][dir]), front_pos(wall_pos + dir_sign*depth);
			c.d[ dim][ dir] = wall_pos;
			c.d[ dim][!dir] = front_pos + dir_sign*front_clearance;
			c.d[!dim][   0] = center - min_hwidth;
			c.d[!dim][   1] = center + min_hwidth;
			cube_t c_min(c); // min runlength - used for collision tests
			for (unsigned e = 0; e < 2; ++e) {c.d[!dim][e] = cabinet_area.d[!dim][e];} // start at full room width
			bool bad_place(0);

			for (auto i = blockers.begin(); i != blockers.end(); ++i) {
				cube_t b(*i); // expand tables by an extra clearance to allow the player to fit in the diagonal gap between the table and the counter
				if (int(i - blockers.begin()) == table_blocker_ix) {b.expand_in_dim(!dim, min_clearance);}
				if (!b.intersects(c)) continue; // optimization - no cube interaction
				if (b.intersects(c_min)) {bad_place = 1; break;}
				if (b.d[!dim][1] < c_min.d[!dim][0]) {max_eq(c.d[!dim][0], b.d[!dim][1]);} // clip on lo side
				if (b.d[!dim][0] > c_min.d[!dim][1]) {min_eq(c.d[!dim][1], b.d[!dim][0]);} // clip on hi side
			} // for i
			if (bad_place) continue;
			if (interior->is_blocked_by_stairs_or_elevator(c)) continue; // check stairs
			assert(c.contains_cube(c_min));

			for (auto i = objs.begin()+objs_start; i != objs.begin()+counters_start; ++i) { // remove any chairs blocking the counter doors/drawers
				if (i->type == TYPE_CHAIR && i->intersects(c)) {i->remove();}
			}
			c.d[dim][!dir] = front_pos; // remove front clearance
			bool const add_backsplash(!is_ext_wall); // only add to interior walls to avoid windows; assuming not in basement

			for (auto i = objs.begin()+counters_start; i != objs.end(); ++i) { // find adjacencies to previously placed counters and flag to avoid placing doors
				if (i->dim == dim) continue; // not perpendicular
				if (i->d[!i->dim][dir] != wall_pos) continue; // not against the wall on this side
				if (i->d[i->dim][i->dir] != c.d[!dim][0] && i->d[i->dim][i->dir] != c.d[!dim][1]) continue; // not adjacent
				i->flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
				if (add_backsplash) {i->flags |= RO_FLAG_HAS_EXTRA;}
			}
			unsigned const cabinet_id(objs.size());
			objs.emplace_back(c, (is_sink ? TYPE_KSINK : TYPE_COUNTER), room_id, dim, !dir, 0, tot_light_amt);
			set_obj_id(objs);
			had_counter = 1;
			
			if (add_backsplash) {
				objs.back().flags |= (RO_FLAG_ADJ_BOT | RO_FLAG_HAS_EXTRA); // flag back as having a back backsplash
				cube_t bs(c);
				bs.z1()  = c.z2();
				bs.z2() += BACKSPLASH_HEIGHT*c.dz();
				bs.d[dim][!dir] -= (dir ? -1.0 : 1.0)*0.99*depth; // matches building_room_geom_t::add_counter()
				objs.emplace_back(bs, TYPE_BLOCKER, room_id, dim, !dir, RO_FLAG_INVIS); // add blocker to avoid placing light switches here
			}
			// add upper cabinets
			cube_t c2(c);
			set_cube_zvals(c2, (zval + 0.65*vspace), cabinet_area.z2()); // up to the ceiling

			if (is_ext_wall) { // possibly against a window
				max_eq(c2.z1(), (c2.z2() - vspace*get_window_v_border() + 0.5f*floor_thickness)); // increase bottom of cabinet to the top of the window
			}
			if (c2.dz() > 0.1*vspace && !has_bcube_int_no_adj(c2, blockers)) { // add if it's not too short and not blocked
				objs.emplace_back(c2, TYPE_CABINET, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // no collision detection
				set_obj_id(objs);
			}
			blockers.push_back(c); // add to blockers so that later counters don't intersect this one
			bool added_obj(0);
			cube_t ptroll;

			// place a microwave on a counter 50% of the time
			if (!is_sink && !placed_mwave && c.get_sz_dim(!dim) > 0.5*vspace && rgen.rand_bool()) {
				float const mheight(rgen.rand_uniform(1.0, 1.2)*0.14*vspace), mwidth(1.7*mheight), mdepth(1.2*mheight); // fixed AR=1.7 to match the texture
				float const pos(rgen.rand_uniform((c.d[!dim][0] + 0.6*mwidth), (c.d[!dim][1] - 0.6*mwidth)));
				set_cube_zvals(mwave, c.z2(), c.z2()+mheight);
				set_wall_width(mwave, pos, 0.5*mwidth, !dim);
				mwave.d[dim][ dir] = wall_pos + dir_sign*0.05*mdepth;
				mwave.d[dim][!dir] = mwave.d[dim][dir] + dir_sign*mdepth;
				objs.emplace_back(mwave, TYPE_MWAVE, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);
				objs[cabinet_id].flags |= RO_FLAG_ADJ_TOP; // flag as having a microwave so that we don't add a book or bottle that could overlap it
				placed_mwave = added_obj = 1;
			}
			// place a toaster on a counter 90% of the time
			if (!is_sink && !placed_toaster && have_toaster && rgen.rand_float() < 0.9) {
				float const theight(0.09*vspace), twidth(theight*toaster_sz.x/toaster_sz.z), tdepth(theight*toaster_sz.y/toaster_sz.z);

				if (c.get_sz_dim(!dim) > 1.25*twidth && c.get_sz_dim(dim) > 1.25*tdepth) { // add if it fits
					float const pos_w(rgen.rand_uniform((c.d[!dim][0] + 0.6*twidth), (c.d[!dim][1] - 0.6*twidth)));
					float const pos_d(rgen.rand_uniform((c.d[ dim][0] + 0.6*tdepth), (c.d[ dim][1] - 0.6*tdepth)));
					set_cube_zvals(toaster, c.z2(), c.z2()+theight);
					set_wall_width(toaster, pos_w, 0.5*twidth, !dim);
					set_wall_width(toaster, pos_d, 0.5*tdepth,  dim);

					if (!placed_mwave || !mwave.intersects(toaster)) { // don't overlap the microwave
						objs.emplace_back(toaster, TYPE_TOASTER, room_id, !dim, rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt); // random dir
						objs.back().color = get_toaster_color(rgen);
						objs[cabinet_id].flags |= RO_FLAG_ADJ_TOP; // flag as having an object so that we don't add a book or bottle that could overlap it
						placed_toaster = added_obj = 1;
					}
				}
			}
			// place a paper towel roll on the counter 75% of the time if there's no toaster, microwave, or sink; max 2 per kitchen
			if (!is_sink && !added_obj && num_paper_towels < 2 && rgen.rand_float() < 0.75) {
				float const ptheight(0.115*vspace), ptradius(0.25*ptheight);

				if (min(c.dx(), c.dy()) > 2.5*ptradius) { // add if it fits
					gen_xy_pos_for_round_obj(ptroll, c, ptradius, ptheight, 1.2*ptradius, rgen);
					objs.emplace_back(ptroll, TYPE_TPROLL, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_HAS_EXTRA), tot_light_amt);
					objs[cabinet_id].flags |= RO_FLAG_ADJ_TOP; // flag as having an object so that we don't add a book or bottle that could overlap it
					++num_paper_towels;
				}
			}
			// place milk carton, at most one per kitchen
			if (have_milk && !is_sink && !placed_milk && rgen.rand_float() < 0.5) {
				static vect_cube_t avoid;
				avoid.clear();
				if (placed_mwave  ) {avoid.push_back(mwave  );}
				if (placed_toaster) {avoid.push_back(toaster);}
				if (!ptroll.is_all_zeros()) {avoid.push_back(ptroll);}
				placed_milk = place_milk_on_obj(rgen, c, room_id, tot_light_amt, avoid);
			}
			if (is_sink) { // kitchen sink; add cups, plates, and cockroaches
				cube_t sink(get_sink_cube(objs[cabinet_id]));
				sink.z2() = sink.z1(); // shrink to zero area at the bottom
				unsigned const objs_start(objs.size()), num_objs(1 + rgen.rand_bool()); // 1-2 objects

				for (unsigned n = 0; n < num_objs; ++n) {
					unsigned const obj_type(rgen.rand()%3);
					static vect_cube_t avoid;
					avoid.clear();
					if (objs.size() > objs_start) {avoid.push_back(objs.back());} // avoid the last object that was placed, if there was one

					if      (obj_type == 0) {place_plate_on_obj(rgen, sink, room_id, tot_light_amt, avoid);} // add a plate
					else if (obj_type == 1) {place_cup_on_obj  (rgen, sink, room_id, tot_light_amt, avoid, 1);} // add a cup; make_empty=1
					else if (obj_type == 2 && building_obj_model_loader.is_model_valid(OBJ_MODEL_ROACH)) { // add a cockroach (upside down?)
						sink.d[dim][!dir] = sink.get_center_dim(dim); // use the half area near the back wall to make sure the roach is visible to the player
						cube_t roach;
						float const radius(sink.get_sz_dim(dim)*rgen.rand_uniform(0.08, 0.12)), height(get_cockroach_height_from_radius(radius));
						gen_xy_pos_for_round_obj(roach, sink, radius, height, 1.1*radius, rgen);
						objs.emplace_back(roach, TYPE_ROACH, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt);
						if (rgen.rand_bool()) {objs.back().flags |= RO_FLAG_BROKEN;} // 50% chance it's dead
					}
				} // for n
			}
			is_sink = 0; // sink is in first placed counter only
		} // for n
		if (had_counter && rgen.rand_float() < 0.75) { // maybe place an island near the room center
			room_object_t const &sink(objs[counters_start]); // first counter is the kitchen sink
			bool const dim(sink.dim);

			if (cabinet_area.get_sz_dim(dim) > 6.0*depth) { // kitchen is wide enough
				float const wall_spacing(2.6*depth);
				cube_t island(sink);
				island.translate_dim(dim, (sink.dir ? 1.0 : -1.0)*3.0*depth);
				max_eq(island.d[!dim][0], cabinet_area.d[!dim][0]+wall_spacing);
				min_eq(island.d[!dim][1], cabinet_area.d[!dim][1]-wall_spacing);
				
				if (island.get_sz_dim(!dim) > depth && cabinet_area.contains_cube(island) &&
					!has_bcube_int(island, blockers) && !interior->is_blocked_by_stairs_or_elevator(island))
				{
					for (auto i = objs.begin()+objs_start; i != objs.begin()+counters_start; ++i) { // remove any chairs intersecting the island
						if (i->type == TYPE_CHAIR && i->intersects(island)) {i->remove();}
					}
					objs.emplace_back(island, TYPE_COUNTER, room_id, dim, !sink.dir, 0, tot_light_amt); // faces opposite the sink
					set_obj_id(objs);
				}
			}
		}
	}
	if (placed_obj && building_obj_model_loader.is_model_valid(OBJ_MODEL_BAN_PEEL) && rgen.rand_bool()) { // maybe place a banana peel on the floor
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_BAN_PEEL));
		float length(0.083*vspace), width(length*sz.y/sz.x), height(length*sz.z/sz.x);
		cube_t valid_area(place_area);
		valid_area.expand_by_xy(-0.25*vspace); // not too close to a wall
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random orientation
		vector3d size(0.5*length, 0.5*width, height);
		if (dim) {swap(size.x, size.y);}

		if (valid_area.dx() > 2.0*size.x && valid_area.dy() > 2.0*size.y) { // should always be true
			for (unsigned n = 0; n < 4; ++n) { // make 4 attempts to place the object
				cube_t c(gen_xy_pos_in_area(valid_area, size, rgen, zval));
				c.expand_by_xy(size);
				c.z2() += height;
				if (overlaps_other_room_obj(c, objs_start) || is_obj_placement_blocked(c, room, 1)) continue; // bad placement
				objs.emplace_back(c, TYPE_BAN_PEEL, room_id, dim, dir, RO_FLAG_RAND_ROT, tot_light_amt);
				break; // done
			} // for n
		}
	}
	return placed_obj;
}

// commercial kitchen appliances, read from model config file; should agree with the set of models specified in the config file
enum {KCA_BI_OVEN=0, KCA_DEEP_FRYER, KCA_SINK, KCA_SINK2, KCA_OVEN, KCA_OMWC, KCA_LP_FRYER, KCA_LP_DISHWASHER, KCA_LP_FRIDGE, KCA_LP_GRILL, KCA_LP_OVEN, KCA_LP_STOVE, NUM_KC_APP};
unsigned const cka_classes[NUM_KC_APP] = {1, 0, 2, 2, 1, 1, 0, 2, 2, 0, 1, 0}; // 0=needs hood, 1=cooks, 2=doesn't cook

cube_t get_com_kitchen_app_coll_bcube(room_object_t const &app) {
	cube_t c(app);
	assert(app.type == TYPE_KITCH_APP);
	assert(app.item_flags < NUM_KC_APP);
	if      (app.item_flags == KCA_DEEP_FRYER) {c.z2() -= 0.17*app.dz();}
	else if (app.item_flags == KCA_SINK      ) {c.z2() -= 0.51*app.dz();}
	else if (app.item_flags == KCA_LP_FRYER  ) {c.z2() -= 0.22*app.dz();}
	else if (app.item_flags == KCA_LP_GRILL  ) {c.z2() -= 0.18*app.dz();}
	else if (app.item_flags == KCA_LP_STOVE  ) {c.z2() -= 0.13*app.dz();}
	return c;
}
unsigned get_com_kitchen_app_coll_cubes(room_object_t const &app, cube_t cubes[3]) {
	assert(app.type == TYPE_KITCH_APP);

	if (app.item_flags == KCA_LP_DISHWASHER) { // 3 parts, excluding plates: {bottom, back middle, top}
		cubes[0] = cubes[1] = cubes[2] = app;
		cubes[0].z2() = cubes[1].z1() = app.z1() + 0.35*app.dz();
		cubes[1].z2() = cubes[2].z1() = app.z1() + 0.65*app.dz();
		cubes[1].d[app.dim][app.dir] -= (app.dir ? 1.0 : -1.0)*0.95*app.get_depth(); // shrink to the back
		return 3;
	}
	cubes[0] = get_com_kitchen_app_coll_bcube(app);

	if (app.item_flags == KCA_DEEP_FRYER || app.item_flags == KCA_LP_FRYER || app.item_flags == KCA_LP_GRILL || app.item_flags == KCA_LP_STOVE) {
		float val(0.0);
		if      (app.item_flags == KCA_DEEP_FRYER) {val = 0.96;}
		else if (app.item_flags == KCA_LP_FRYER  ) {val = 0.91;}
		else if (app.item_flags == KCA_LP_GRILL  ) {val = 0.93;}
		else if (app.item_flags == KCA_LP_STOVE  ) {val = 0.96;}
		cubes[1] = app;
		cubes[1].d[app.dim][app.dir] -= (app.dir ? 1.0 : -1.0)*val*app.get_depth(); // shrink to the back
		return 2; // 2 parts
	}
	return 1; // 1 part
}

struct ck_app_model_t {
	unsigned app_type=0, model_id=0, obj_counts[NUM_KC_APP]={};
	vector<unsigned> cands;

	void assign(rand_gen_t &rgen, unsigned cka_class) {
		unsigned min_count(1000);
		cands.clear();

		for (unsigned n = 0; n < NUM_KC_APP; ++n) {
			if (cka_classes[n] == cka_class) {min_eq(min_count, obj_counts[n]);}
		}
		for (unsigned n = 0; n < NUM_KC_APP; ++n) {
			if (cka_classes[n] == cka_class && obj_counts[n] == min_count) {cands.push_back(n);}
		}
		assert(!cands.empty()); // must have at least one cand for each class
		app_type = cands[rgen.rand() % cands.size()];
		model_id = combine_model_submodel_id(OBJ_MODEL_CK_APP, app_type);
	}
	float    get_height  () const {return 0.5*building_obj_model_loader.get_model(model_id).scale;} // in units of floor spacing
	vector3d get_sz_scale() const {return     building_obj_model_loader.get_model_world_space_size(model_id);}
	void post_add() {++obj_counts[app_type];}
};

void set_specular_for_low_poly_kitchen_models() {
	static bool specular_was_set=0;
	if (specular_was_set) return;
	specular_was_set = 1;

	for (unsigned i = KCA_LP_FRYER; i < NUM_KC_APP; ++i) {
		unsigned const model_id(combine_model_submodel_id(OBJ_MODEL_CK_APP, i));
		building_obj_model_loader.set_all_material_specular(model_id, WHITE, 60.0, 1.0); // fully specular metal
	}
}
void move_lights_to_not_intersect(vect_room_object_t &objs, unsigned objs_start, unsigned objs_end, cube_t const &avoid) {
	assert(objs_start <= objs_end && objs_end <= objs.size());

	for (auto i = objs.begin()+objs_start; i != objs.begin()+objs_end; ++i) { // move light away from avoid; this is rare
		if (i->type != TYPE_LIGHT || !i->intersects(avoid)) continue; // or intersects_xy()?
		bool const dim(i->dy() < i->dx()), dir(avoid.get_center_dim(dim) < i->get_center_dim(dim)); // move in narrow dim
		i->translate_dim(dim, (avoid.d[dim][dir] - i->d[dim][!dir])); // hopefully this pos is inside the room and not intersecting anything
	}
}

void building_t::add_commercial_kitchen_app_post(unsigned obj_ix, unsigned app_type, cube_t &hood, unsigned cclass_counts[3], rand_gen_t &rgen) {
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(obj_ix < objs.size());
	room_object_t const app(objs[obj_ix]); // deep copy
	bool const dim(app.dim), dir(app.dir);
	float const height(app.dz());
	cube_t frame(app);
	float z_shift(0.0);
	unsigned skip_faces(0); // 0=don't draw
	bool add_platform(0);

	if (app_type == KCA_BI_OVEN || app_type == KCA_OVEN || app_type == KCA_OMWC) { // front face only; add sides, top, and back
		unsigned const front_face_mask(get_face_mask(dim, dir));
		float const fb_shift((app_type == KCA_OMWC || app_type == KCA_OVEN) ? 0.06 : 0.09);
		skip_faces = ~front_face_mask | EF_Z1;
		frame.d[dim][dir] -= (dir ? 1.0 : -1.0)*fb_shift*app.get_depth(); // shift front edge behind door

		if (app_type != KCA_OMWC) {
			z_shift      = 0.6*height;
			add_platform = 1;
			frame.translate_dim(2, z_shift);
		}
		if (app_type == KCA_OVEN) { // need an extra plate near the front
			cube_t c(frame);
			c.z1() = frame.z2() - 0.25*height;
			objs.emplace_back(c, TYPE_METAL_BAR, app.room_id, dim, dir, 0, app.light_amt, SHAPE_CUBE, WHITE, front_face_mask);
		}
	}
	else if (app_type == KCA_SINK) { // top face only; draw front, sides, and back
		z_shift = 0.3*height;
		frame.z2() = app.z1() + z_shift + 0.49*height;
		skip_faces = EF_Z12;
	}
	else if (app_type == KCA_SINK2) { // top face only; draw front, sides, and back
		z_shift    = 1.1*height;
		frame.z2() = app.z1() + z_shift + 0.96*height;
		skip_faces = EF_Z12;
	}
	else if (app_type == KCA_DEEP_FRYER) {
		z_shift = 0.3*height;
		add_platform = 1;
	}
	else if (app_type == KCA_LP_GRILL) {add_pan_on_grill(app, rgen);}
	else if (app_type == KCA_LP_STOVE) {add_pan_on_stove(app, rgen);}

	if (add_platform) { // place this object on a platform/table; should this be merged across adjacent objects at the same height?
		assert(z_shift > 0.0);
		cube_t platform(app);
		platform.z2() = app.z1() + z_shift;
		objs.emplace_back(platform, TYPE_METAL_BAR, app.room_id, dim, dir, 0, app.light_amt, SHAPE_CUBE, WHITE, EF_Z1);
	}
	if (z_shift > 0.0) {objs[obj_ix].translate_dim(2, z_shift);} // shift up
	// add missing sides of this model
	if (skip_faces) {objs.emplace_back(frame, TYPE_METAL_BAR, app.room_id, dim, dir, RO_FLAG_NOCOLL, app.light_amt, SHAPE_CUBE, WHITE, skip_faces);}
	unsigned const app_cclass(cka_classes[app_type]);
	++cclass_counts[app_cclass];
	if (app_cclass == 0) {hood.assign_or_union_with_cube(app);}
}

bool building_t::add_commercial_kitchen_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, unsigned floor_ix,
	float light_amt, unsigned objs_start, unsigned lights_start, light_ix_assign_t &light_ix_assign)
{
	if (room_has_stairs_or_elevator(room, zval, floor_ix)) return 0; // works, but stairs may be blocked by a trolley
	float const floor_spacing(get_window_vspace()), wall_thick(get_wall_thickness()), trim_thick(get_trim_thickness());
	vector2d const room_sz(room.get_size_xy());
	if (room_sz.get_min_val() < 2.0*floor_spacing || room_sz.get_max_val() < 3.0*floor_spacing) return 0; // too small
	bool const in_mall(room.is_ext_basement() && has_mall()), dim(room_sz.x < room_sz.y); // long dim
	bool const add_island(room.get_sz_dim(!dim) > 3.0*floor_spacing); // if room is wide enough, add a center island
	float const ceil_zval((in_mall ? room.z2() : (zval + floor_spacing)) - get_fc_thickness()), clearance(get_min_front_clearance_inc_people());
	cube_t place_area(get_room_wall_bounds(room));
	place_area.expand_by_xy(-trim_thick);
	if (!has_tile_floor() && !in_mall) {zval = add_flooring(room, zval, room_id, light_amt, FLOORING_LGTILE);} // prison and maybe school
	vect_room_object_t &objs(interior->room_geom->objs);

	if (add_island) { // add a large table in the center of the room
		float const room_len(room.get_sz_dim(!dim));
		float const table_len(min(0.5f*room_len, (room_len - 2.0f*floor_spacing))), table_width(min(0.4f*table_len, 0.6f*floor_spacing));
		unsigned const item_flags(get_table_item_flags(rgen, 0, 1)); // is_plastic=0, is_metal=1
		cube_t table;
		set_cube_zvals(table, zval, (zval + 0.3*floor_spacing + wall_thick));
		set_wall_width(table, room.get_center_dim(!dim), 0.5*table_width, !dim);
		set_wall_width(table, room.get_center_dim( dim), 0.5*table_len,    dim);
		objs.emplace_back(table, TYPE_TABLE, room_id, dim, 0, RO_FLAG_ADJ_TOP, light_amt, SHAPE_CUBE, WHITE, item_flags); // metal; assumes table has something on it
		unsigned const avoid_start(objs.size());
		vect_cube_t avoid;
		
		// add microwave to table
		bool const mdim(!dim), mdir(rgen.rand_bool());
		float const mheight(rgen.rand_uniform(1.0, 1.2)*0.16*floor_spacing), mwidth(1.7*mheight), mdepth(1.2*mheight); // fixed AR=1.7 to match the texture
		float const pos(rgen.rand_uniform((table.d[!mdim][0] + 0.6*mwidth), (table.d[!mdim][1] - 0.6*mwidth)));
		cube_t mwave;
		set_cube_zvals(mwave, table.z2(), table.z2()+mheight);
		set_wall_width(mwave, pos, 0.5*mwidth, !mdim);
		mwave.d[mdim][!mdir] = table.d[mdim][!mdir] + (mdir ? 1.0 : -1.0)*rgen.rand_uniform(0.0, 0.1)*mdepth;
		mwave.d[mdim][ mdir] = mwave.d[mdim][!mdir] + (mdir ? 1.0 : -1.0)*mdepth;
		objs.emplace_back(mwave, TYPE_MWAVE, room_id, mdim, !mdir, RO_FLAG_NOCOLL, light_amt);
		avoid.push_back(mwave);
		// add deep fryer to table
		unsigned const model_id(combine_model_submodel_id(OBJ_MODEL_CK_APP, KCA_DEEP_FRYER));
		float const dp_height(0.7*0.5*building_obj_model_loader.get_model(model_id).scale); // smaller version
		
		if (place_model_along_wall(model_id, TYPE_KITCH_APP, room, dp_height, rgen, table.z2(), room_id, light_amt, table, avoid_start)) {
			objs.back().dir ^= 1; // backwards
			avoid.push_back(objs.back());
		}
		// add toaster to table
		if (place_model_along_wall(OBJ_MODEL_TOASTER, TYPE_TOASTER, room, 0.09, rgen, table.z2(), room_id, light_amt, table, avoid_start)) {
			avoid.push_back(objs.back());
		}
		// place objects on the table; similar to building_t::add_restaurant_counter()
		unsigned const num_cont (5 + (rgen.rand() % 5)); // 5-9
		unsigned const num_pizza(0 + (rgen.rand() % 3)); // 0-2

		for (unsigned n = 0; n < num_cont; ++n) { // containers
			switch (rgen.rand() % 5) {
			case 0: // cup
				if (place_cup_on_obj  (rgen, table, room_id, light_amt, avoid)) {avoid.push_back(objs.back());}
				break;
			case 1: // plate or bowl
				if (place_plate_on_obj(rgen, table, room_id, light_amt, avoid, (rgen.rand_float() < 0.35))) {avoid.push_back(objs.back());}
				break;
			case 2: // pan
				if (place_pan_on_obj  (rgen, table, room_id, light_amt, avoid)) {avoid.push_back(get_pan_bcube_inc_handle(objs.back()));}
				break;
			case 3: // milk
				if (place_milk_on_obj (rgen, table, room_id, light_amt, avoid)) {avoid.push_back(objs.back());}
			case 4: // tray
				add_cafeteria_tray_to_surface(table, !dim, room_id, light_amt, avoid, rgen);
				break;
			} // end switch
		} // for n
		for (unsigned n = 0; n < num_pizza; ++n) { // pizza
			if (place_pizza_on_obj(rgen, table, room_id, light_amt, avoid)) {avoid.push_back(objs.back());} // random dim/dir
		}
	}
	// add walk-in freezer as a "closet" type
	unsigned closet_obj_id(0);
	vect_cube_t duct_avoid;
	
	if (add_closet_to_room(rgen, room, zval, room_id, objs_start, RTYPE_KITCHEN, 0, clearance, closet_obj_id, light_ix_assign)) { // bed_obj_ix=0 (not set)
		room_object_t &freezer(objs[closet_obj_id]);
		cube_t avoid(freezer);
		avoid.expand_in_dim(freezer.dim, 1.1*get_doorway_width()); // allow space for door to open
		duct_avoid.push_back(avoid);
		if (in_mall) {freezer.flags |= RO_FLAG_IN_MALL;}
		move_lights_to_not_intersect(objs, lights_start, objs_start, avoid); // or freezer?
	}
	// place commerial kitchen appliances (grills, deep fryers, ovens, sinks, etc.); only legal if all models have been loaded
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_CK_APP) && building_obj_model_loader.get_num_sub_models(OBJ_MODEL_CK_APP) == NUM_KC_APP) {
		set_specular_for_low_poly_kitchen_models();
		float const app_gap(trim_thick); // must be large enough to prevent cube intersection
		unsigned fail_count(0), num_fridge(0), cclass_counts[3] = {0, 1, 1}; // prefer hood items
		cclass_counts[rgen.rand_bool() + 1] = 2; // randomize second two classes
		ck_app_model_t app_model; // reused across calls; tracks model stats
		cube_t hood;
		bool add_fridge_pass(1); // add fridge on first and every alternating non-cooking class pass

		for (unsigned n = 0; n < 20; ++n) { // place up to 20 kitchen appliance models groups; stop after enough failures
			unsigned const cka_class(min_element(cclass_counts, cclass_counts+3) - cclass_counts);
			bool const non_cook_class(cka_class == 2), add_fridge(non_cook_class && add_fridge_pass);
			unsigned obj_ix(objs.size());
			bool success(0);
			hood.set_to_zeros();

			if (add_fridge) {
				success = place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, room, 0.75, rgen, zval, room_id, light_amt, place_area, objs_start, 1.2, 4, 0, WHITE, 1);
				num_fridge += success;
			}
			else {
				app_model.assign(rgen, cka_class);
				float const height(app_model.get_height());
				success = place_model_along_wall(app_model.model_id, TYPE_KITCH_APP, room, height, rgen, zval, room_id, light_amt, place_area, objs_start, 1.0, 4, 0, WHITE, 1);
			}
			if (!success) {
				if (++fail_count >= 10) break; // stop after 10 failures
				continue;
			}
			if (non_cook_class && num_fridge > 1) {add_fridge_pass ^= 1;}
			room_object_t const app(objs[obj_ix]); // deep copy to avoid reference invalidation
			bool const adim(app.dim), adir(app.dir);
			
			if (!add_fridge) {
				assert(app.item_flags == app_model.app_type);
				app_model.post_add();
				add_commercial_kitchen_app_post(obj_ix, app_model.app_type, hood, cclass_counts, rgen);
			}
			for (unsigned d = 0; d < 2; ++d) { // for each side of this appliance
				cube_t prev(app);

				for (unsigned m = 0; m < 4; ++m) { // up to 4 neighbors
					app_model.assign(rgen, cka_class); // same class
					vector3d const sz_scale(app_model.get_sz_scale()); // D, W, H
					float const height(floor_spacing*app_model.get_height()), depth(height*sz_scale.x/sz_scale.z), width(height*sz_scale.y/sz_scale.z);
					room_object_t app2(app);
					app2.d[!adim][!d  ] = prev.d[!adim][ d   ] + (d    ? 1.0 : -1.0)*app_gap; // adj edge
					app2.d[!adim][ d  ] = app2.d[!adim][!d   ] + (d    ? 1.0 : -1.0)*width  ; // far edge
					app2.d[ adim][adir] = app .d[ adim][!adir] + (adir ? 1.0 : -1.0)*depth  ; // front
					set_cube_zvals(app2, zval, (zval + height));
					float const front(app2.d[adim][adir]);
					cube_t tc(app2);
					tc.d[adim][adir] = front + (adir ? 1.0 : -1.0)*clearance; // extend in front
					if (!place_area.contains_cube_xy(tc)        || is_obj_placement_blocked(tc, room, 1)) break;
					if (overlaps_other_room_obj(tc, objs_start) || check_if_against_window (tc, room, adim, !adir)) break;
					app2.type = TYPE_KITCH_APP; // in case app was a fridge
					app2.item_flags = app_model.app_type;
					objs.push_back(app2);
					app_model.post_add();
					add_commercial_kitchen_app_post(objs.size()-1, app_model.app_type, hood, cclass_counts, rgen);
					prev = app2;
					cube_t blocker(tc); // add a blocker for front clearance
					blocker.d[adim][!adir] = front;
					objs.emplace_back(blocker, TYPE_BLOCKER, room_id, adim, adir, RO_FLAG_INVIS, light_amt);
				} // for m
			} // for d
			if (!hood.is_all_zeros() && hood.get_sz_dim(!adim) > 0.35*floor_spacing) { // add hoods/vents where needed, if not too small
				float const hood_z1(max((ceil_zval - 1.0f*floor_spacing), (zval + 0.6f*floor_spacing)));
				set_cube_zvals(hood, hood_z1, ceil_zval);
				objs.emplace_back(hood, TYPE_VENT_HOOD, room_id, adim, adir, 0, light_amt, SHAPE_CUBE, WHITE);
				move_lights_to_not_intersect(objs, lights_start, objs_start, hood);
				duct_avoid.push_back(hood);
			}
		} // for n
	}
	else { // can't load all appliance models; place only fridge instead
		unsigned num_fridges((rgen.rand() % 3) + 2); // 2-4

		for (unsigned i = 0; i < num_fridges; ++i) {
			place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, room, 0.75, rgen, zval, room_id, light_amt, place_area, objs_start, 1.2, 4, 0, WHITE, 1);
		}
	}
	if (!in_mall) { // add ceiling ducts; mall already has them added
		unsigned const skip_dir(2), ducts_start(objs.size()); // allow both dirs, depending on what placement is legal
		add_ceiling_ducts(room, ceil_zval, room_id, dim, skip_dir, light_amt, 0, 1, 1, rgen, 0.5, duct_avoid); // cylin_ducts=0, skip_ends=1, skip_top=1, sz_scale=0.5

		for (auto i = objs.begin()+ducts_start; i != objs.end(); ++i) {
			if (i->type == TYPE_DUCT) {move_lights_to_not_intersect(objs, lights_start, objs_start, *i);}
		}
	}
	// what about placing appliances on the sides of the center table?
	add_mwave_on_table  (rgen, room, zval, room_id, light_amt, objs_start, place_area, 0, 1); // plastic=0, metal=1
	add_corner_trashcans(rgen, room, zval, room_id, light_amt, objs_start, dim, 1); // both_ends=1
	// add trolleys with plates; seems like this can work in a kitchen
	unsigned num_trolleys((rgen.rand() % 2) + 2); // 2-3
	vect_cube_t blockers;
	
	for (unsigned i = 0; i < num_trolleys; ++i) {
		unsigned const obj_ix(objs.size());
		if (!add_trolley(rgen, place_area, cube_t(), zval, room_id, light_amt, objs_start)) continue;
		if (rgen.rand_float() < 0.35) continue; // no plates
		room_object_t const trolley(objs[obj_ix]);
		cube_t plate_area(trolley);
		plate_area.expand_in_dim( trolley.dim, -0.2*trolley.get_sz_dim( trolley.dim)); // shorten length
		plate_area.expand_in_dim(!trolley.dim, -0.1*trolley.get_sz_dim(!trolley.dim)); // shorten width
		float const height(trolley.dz()), plate_radius(rgen.rand_uniform(0.25, 0.3)*min(plate_area.dx(), plate_area.dy()));
		set_cube_zvals(plate_area, plate_area.z2()-0.13*height, plate_area.z2()+0.4*height);
		add_stack_of_plates(plate_area, plate_radius, room_id, light_amt, RO_FLAG_NOCOLL, rgen, blockers, objs);
	} // for i
	if (rgen.rand_bool()) {add_buckets_to_room(rgen, place_area, zval, room_id, light_amt, objs_start, 1);} // maybe add a bucket
	if (!in_mall) {add_door_sign("Kitchen", room, zval, room_id);}
	// add floor stains
	unsigned const num_stains(rgen.rand() % 5); // 0-4
	float const rmax(min(0.2f*floor_spacing, 0.1f*min(place_area.dx(), place_area.dy())));
	add_floor_stains(rgen, place_area, zval, room_id, light_amt, objs_start, num_stains, rmax, 1); // is_food=1
	return 1;
}

bool building_t::place_pan_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid) {
	// somewhat larger than pans placed on stoves
	float const floor_spacing(get_window_vspace()), height(rgen.rand_uniform(0.02, 0.025)*floor_spacing);
	float const radius(rgen.rand_uniform(0.05, 0.06)*floor_spacing), edge_dist(2.0*radius);
	if (min(place_on.dx(), place_on.dy()) < 2.1*edge_dist) return 0; // too small
	cube_t pan_bc(place_cylin_object(rgen, place_on, radius, height, edge_dist)); // add space for the handle
	pan_bc.translate_dim(2, 0.01*height); // fix for Z-fighting
	room_object_t const pan(pan_bc, TYPE_PAN, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, DK_GRAY);
	pan_bc = get_pan_bcube_inc_handle(pan); // include the handle
	if (!place_on.contains_cube_xy(pan_bc) || has_bcube_int(pan_bc, avoid)) return 0; // only make one attempt
	interior->room_geom->objs.push_back(pan);
	return 1;
}
void building_t::add_pan_on_grill(room_object_t const &grill, rand_gen_t &rgen) {
	if (rgen.rand_bool()) return; // does a pan belong on a grill? only add half the time
	cube_t top(grill);
	top.z2() -= 0.17*grill.dz();
	top.d[grill.dim][!grill.dir] += (grill.dir ? 1.0 : -1.0)*0.1*grill.get_depth(); // shrink off back part
	place_pan_on_obj(rgen, top, grill.room_id, grill.light_amt); // no avoid
}
void building_t::add_pan_on_stove(room_object_t const &stove, rand_gen_t &rgen) { // commercial stove
	float const depth(stove.get_depth());
	cube_t top(stove);
	top.z2() -= 0.07*stove.dz();
	top.expand_by_xy(-0.07*depth);
	cube_t burner(top);
	for (unsigned d = 0; d < 2; ++d) {burner.d[d][rgen.rand_bool()] = top.get_center_dim(d);} // pick a random quadrant
	float const floor_spacing(get_window_vspace()), height(rgen.rand_uniform(0.02, 0.025)*floor_spacing), radius(rgen.rand_uniform(0.05, 0.06)*floor_spacing);
	cube_t pan;
	pan.set_from_point(cube_top_center(burner));
	pan.expand_by_xy(radius);
	pan.z2() += height;
	interior->room_geom->objs.emplace_back(pan, TYPE_PAN, stove.room_id, stove.dim, rgen.rand_bool(), RO_FLAG_NOCOLL, stove.light_amt, SHAPE_CYLIN, DK_GRAY);
}
bool building_t::place_milk_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_MILK)) return 0;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_MILK));
	float const height(0.06*get_window_vspace()), radius(height*(sz.x + sz.y)/sz.z); // assumes square
	if (min(place_on.dx(), place_on.dy()) < 2.5*radius) return 0; // surface is too small to place this milk
	cube_t milk;
	gen_xy_pos_for_round_obj(milk, place_on, radius, height, 1.1*radius, rgen);
	if (has_bcube_int(milk, avoid)) return 0; // only make one attempt
	interior->room_geom->objs.emplace_back(milk, TYPE_MILK, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt);
	return 1;
}

bool building_t::add_mwave_on_table(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id,
	float tot_light_amt, unsigned objs_start, cube_t const &place_area, bool is_plastic, bool is_metal)
{
	float const window_vspacing(get_window_vspace()), table_height(rgen.rand_uniform(0.33, 0.36)*window_vspacing);
	float const mheight(rgen.rand_uniform(1.0, 1.2)*0.14*window_vspacing), mwidth(1.7*mheight), mdepth(1.2*mheight); // fixed AR=1.7 to match the texture
	vector3d const sz_scale(rgen.rand_uniform(1.5, 1.8)*mdepth, rgen.rand_uniform(1.5, 1.8)*mwidth, table_height); // depth, width, height
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const table_obj_ix(objs.size());
	if (!place_obj_along_wall(TYPE_TABLE, room, table_height, sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 1)) return 0;
	room_object_t &table(objs[table_obj_ix]);
	bool const dim(table.dim), dir(table.dir);
	float const pos(rgen.rand_uniform((table.d[!dim][0] + 0.6*mwidth), (table.d[!dim][1] - 0.6*mwidth)));
	table.flags     |= RO_FLAG_ADJ_TOP; // flag as having a microwave so that we don't add a book or bottle that could overlap it
	table.item_flags = get_table_item_flags(rgen, is_plastic, is_metal); // sets table texture
	cube_t mwave;
	set_cube_zvals(mwave, table.z2(), table.z2()+mheight);
	set_wall_width(mwave, pos, 0.5*mwidth, !dim);
	mwave.d[dim][!dir] = table.d[dim][!dir] + (dir ? 1.0 : -1.0)*0.05*mdepth;
	mwave.d[dim][ dir] = mwave.d[dim][!dir] + (dir ? 1.0 : -1.0)*mdepth;
	objs.emplace_back(mwave, TYPE_MWAVE, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt);
	return 1;
}

bool building_t::place_eating_items_on_table(rand_gen_t &rgen, unsigned table_obj_id) {
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(table_obj_id < objs.size());
	room_object_t const table(objs[table_obj_id]); // deep copy to avoid invalidating the reference
	float const floor_spacing(get_window_vspace()), plate_radius(get_plate_radius(rgen, table, floor_spacing)), plate_height(0.1*plate_radius), spacing(1.33*plate_radius);
	unsigned const objs_size(objs.size());
	bool added_obj(0);

	for (unsigned i = (table_obj_id + 1); i < objs_size; ++i) { // iterate over chairs (by index, since we're adding to objs)
		if (objs[i].type != TYPE_CHAIR) break; // done with chairs for this table
		point const chair_center(objs[i].get_cube_center()), table_center(table.get_cube_center());
		point pos;

		if (table.shape == SHAPE_CYLIN) { // circular
			float const dist(table.get_radius() - spacing);
			pos = table_center + dist*(chair_center - table_center).get_norm();
		}
		else { // rectangular
			cube_t place_bounds(table);
			place_bounds.expand_by_xy(-spacing);
			pos = place_bounds.closest_pt(chair_center);
		}
		unsigned const plate_flags(RO_FLAG_NOCOLL | (table.is_glass_table() ? RO_FLAG_ADJ_BOT : 0)); // flag for bottom to be drawn if on a glass table
		cube_t plate;
		plate.set_from_sphere(pos, plate_radius);
		set_cube_zvals(plate, table.z2(), table.z2()+plate_height); // place on the table
		objs.emplace_back(plate, TYPE_PLATE, table.room_id, 0, 0, plate_flags, table.light_amt, SHAPE_CYLIN); // or bowl?
		set_obj_id(objs);
		place_food_on_plate(rgen, plate, table.room_id, table.light_amt);

		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_SILVER)) {
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SILVER)); // D, W, H
			float const sw_height(0.0075*floor_spacing), sw_hwidth(0.5*sw_height*sz.x/sz.z), sw_hlen(0.5*sw_height*sz.y/sz.z); // Note: x/y swapped
			vector3d const offset(pos - table_center);
			bool const dim(objs[i].dim), dir(!objs[i].dir); // use chair dim/dir
			cube_t sw_bc;
			set_cube_zvals(sw_bc, table.z2()+0.1*sw_height, table.z2()+sw_height);
			set_wall_width(sw_bc, pos[!dim] + ((dim ^ dir) ? 1.0 : -1.0)*1.2*(plate_radius + sw_hlen), sw_hlen, !dim);
			set_wall_width(sw_bc, pos[ dim], sw_hwidth, dim);
			objs.emplace_back(sw_bc, TYPE_SILVER, table.room_id, dim, dir, RO_FLAG_NOCOLL, table.light_amt, SHAPE_CUBE, GRAY);
		}
		added_obj = 1;
	} // for i
	if (added_obj && !is_school()) { // place a vase in the center of the table, but not for schools
		float const vase_radius(rgen.rand_uniform(0.35, 0.6)*plate_radius), vase_height(rgen.rand_uniform(2.0, 6.0)*vase_radius);
		cube_t vase;
		vase.set_from_sphere(table.get_cube_center(), vase_radius);
		set_cube_zvals(vase, table.z2(), table.z2()+vase_height); // place on the table
		objs.emplace_back(vase, TYPE_VASE, table.room_id, 0, 0, RO_FLAG_NOCOLL, table.light_amt, SHAPE_CYLIN, gen_vase_color(rgen));
		set_obj_id(objs);
	}
	return added_obj;
}

