// 3D World - Building Interior Room Geometry Drawing
// by Frank Gennari 7/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"

float const TAPE_HEIGHT_TO_RADIUS = 0.6;

float get_lamp_width_scale();
vect_cube_t &get_temp_cubes();
bool get_dishwasher_for_ksink(room_object_t const &c, cube_t &dishwasher);

bool add_if_not_intersecting(room_object_t const &obj, vector<room_object_t> &objects, vect_cube_t &cubes) {
	if (has_bcube_int(obj, cubes)) return 0;
	objects.push_back(obj);
	cubes.push_back(obj);
	return 1;
}
void gen_xy_pos_for_cube_obj(cube_t &C, cube_t const &S, vector3d const &sz, float height, rand_gen_t &rgen) {
	point center;
	for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(S.d[d][0]+sz[d], S.d[d][1]-sz[d]);} // randomly placed within the bounds of the shelf
	C.set_from_point(center);
	set_cube_zvals(C, S.z2(), S.z2()+height);
	C.expand_by_xy(sz);
}
void gen_xy_pos_for_round_obj(cube_t &C, cube_t const &S, float radius, float height, float spacing, rand_gen_t &rgen, bool place_at_z1=0) {
	point center;
	for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform((S.d[d][0] + spacing), (S.d[d][1] - spacing));} // place at least spacing from edge
	C.set_from_sphere(center, radius);
	float const place_z(place_at_z1 ? S.z1() : S.z2());
	set_cube_zvals(C, place_z, place_z+height);
}

void add_boxes_to_space(room_object_t const &c, vector<room_object_t> &objects, cube_t const &bounds, vect_cube_t &cubes, rand_gen_t &rgen,
	unsigned num_boxes, float xy_scale, float hmin, float hmax, bool allow_crates, unsigned flags)
{
	float const bounds_sz[2] = {bounds.dx(), bounds.dy()};
	room_object_t C(c);
	C.flags = flags; // Note: also clears open flag
	vector3d sz;
	point center;

	for (unsigned n = 0; n < num_boxes; ++n) {
		for (unsigned d = 0; d < 2; ++d) {
			sz    [d] = min(xy_scale*rgen.rand_uniform(0.5, 1.0), 0.99f*0.5f*bounds_sz[d]); // x,y half width; clamp to slightly smaller than bounds to avoid an assert
			center[d] = rgen.rand_uniform(bounds.d[d][0]+sz[d], bounds.d[d][1]-sz[d]); // randomly placed within the bounds of the closet
		}
		C.set_from_point(center);
		set_cube_zvals(C, bounds.z1(), (bounds.z1() + rgen.rand_uniform(hmin, hmax)));
		C.expand_by_xy(sz);
		if (has_bcube_int(C, cubes)) continue; // intersects - just skip it, don't try another placement
		C.color  = gen_box_color(rgen);
		C.dim    = c.dim ^ bool(rgen.rand()&3) ^ 1; // make the box label face outside 75% of the time
		C.obj_id = rgen.rand(); // used to select crate texture and box contents
		C.type   = ((!allow_crates || rgen.rand_bool()) ? (room_object)TYPE_BOX : (room_object)TYPE_CRATE);
		objects.push_back(C);
		cubes.push_back(C);
	} // for n
}

cube_t get_closet_interior_space(room_object_t const &c, cube_t const cubes[5]) {
	cube_t interior(c);
	if (!cubes[1].is_all_zeros()) {interior.d[!c.dim][0] = cubes[1].d[!c.dim][1];} // left  side (if wall exists)
	if (!cubes[3].is_all_zeros()) {interior.d[!c.dim][1] = cubes[3].d[!c.dim][0];} // right side (if wall exists)
	float const extra_space(c.is_small_closet() ? 0.0 : (c.dir ? -1.0 : 1.0)*1.2*cubes[4].get_sz_dim(c.dim)); // add some extra space for sliding/folding doors
	interior.d[c.dim][c.dir] = cubes[2].d[c.dim][!c.dir] + extra_space; // set to inside of front wall
	interior.expand_by_xy(-0.02*interior.min_len()); // shrink slightly to clip off the area taken by the wall trim where objects shouldn't be placed
	assert(interior.is_strictly_normalized());
	return interior;
}

void add_obj_to_closet(room_object_t const &c, cube_t const &interior, vector<room_object_t> &objects, vect_cube_t &cubes,
	rand_gen_t &rgen, vector3d const &sz, unsigned obj_type, unsigned flags, room_obj_shape shape=SHAPE_CUBE)
{
	for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts
		point center;
		for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(interior.d[d][0]+sz[d], interior.d[d][1]-sz[d]);}
		cube_t obj;
		obj.set_from_point(center);
		set_cube_zvals(obj, interior.z1(), (interior.z1() + sz.z));
		obj.expand_by_xy(sz);

		if (!has_bcube_int(obj, cubes)) { // check for intersection with boxes
			objects.emplace_back(obj, obj_type, c.room_id, c.dim, c.dir, flags, c.light_amt, shape);
			cubes.push_back(obj);
			break;
		}
	} // for n
}

void building_room_geom_t::add_closet_objects(room_object_t const &c, vector<room_object_t> &objects) {
	cube_t ccubes[5]; // only used to get interior space
	get_closet_cubes(c, ccubes);
	cube_t const interior(get_closet_interior_space(c, ccubes));
	float const depth(interior.get_sz_dim(c.dim)), box_sz(0.25*depth), window_vspacing(c.dz()*(1.0 + FLOOR_THICK_VAL_HOUSE));
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	unsigned const num_boxes((rgen.rand()%3) + (rgen.rand()%4)); // 0-5
	vect_cube_t &cubes(get_temp_cubes());
	add_boxes_to_space(c, objects, interior, cubes, rgen, num_boxes, box_sz, 0.8*box_sz, 1.5*box_sz, 0, flags); // allow_crates=0
	vector3d sz;
	point center;

	if (!c.is_small_closet()) { // larger closets have more random items
		if (rgen.rand_bool()) { // maybe add a lamp in the closet
			float const height(0.25*window_vspacing), width(height*get_lamp_width_scale()), radius(0.5*width);

			if (width > 0.0 && width < 0.9*min(interior.dx(), interior.dy())) { // check if lamp model is valid and lamp fits in closet
				center.assign(0.0, 0.0, interior.z1());

				for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts
					for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(interior.d[d][0]+radius, interior.d[d][1]-radius);}
					cube_t lamp(get_cube_height_radius(center, radius, height));

					if (!has_bcube_int(lamp, cubes)) { // check for intersection with boxes
						objects.emplace_back(lamp, TYPE_LAMP, c.room_id, 0, 0, flags, c.light_amt, SHAPE_CYLIN, lamp_colors[rgen.rand()%NUM_LAMP_COLORS]);
						cubes.push_back(lamp);
						break;
					}
				} // for n
			}
		}
		if (rgen.rand_bool()) { // maybe add a computer in the closet
			float const height(0.21*window_vspacing*rgen.rand_uniform(1.0, 1.2)), cheight(0.75*height), cwidth(0.44*cheight), cdepth(0.9*cheight);
			sz[c.dim] = 0.5*cdepth; sz[!c.dim] = 0.5*cwidth; sz.z = cheight;
			add_obj_to_closet(c, interior, objects, cubes, rgen, sz, TYPE_COMPUTER, flags);
		}
		if (rgen.rand_bool()) { // maybe add a keyboard in the closet
			float const kbd_hwidth(0.12*window_vspacing), kbd_depth(0.6*kbd_hwidth), kbd_height(0.06*kbd_hwidth);
			sz[c.dim] = 0.5*kbd_depth; sz[!c.dim] = kbd_hwidth; sz.z = kbd_height;
			add_obj_to_closet(c, interior, objects, cubes, rgen, sz, TYPE_KEYBOARD, flags);
		}
		if (rgen.rand_bool()) { // maybe add a paint can in the closet
			float const height(0.64*0.2*window_vspacing), radius(0.28*0.2*window_vspacing);
			sz.assign(radius, radius, height);
			add_obj_to_closet(c, interior, objects, cubes, rgen, sz, TYPE_PAINTCAN, flags, SHAPE_CYLIN);
		}
	}
	// add hanger rod
	float const hr_radius(0.007*window_vspacing);
	room_object_t hanger_rod(interior, TYPE_HANGER_ROD, c.room_id, c.dim, c.dir, (RO_FLAG_NOCOLL | RO_FLAG_INTERIOR));
	hanger_rod.z1() = c.z1() + 0.8*window_vspacing;
	hanger_rod.z2() = hanger_rod.z1() + 2.0*hr_radius;
	set_wall_width(hanger_rod, (0.45*c.d[c.dim][c.dir] + 0.55*c.d[c.dim][!c.dir]), hr_radius, c.dim); // move slightly toward the back
	unsigned const hanger_rod_ix(objects.size());
	objects.push_back(hanger_rod);

	// add hangers
	unsigned const num_hangers(rgen.rand() % (c.is_small_closet() ? 6 : 11)); // 0-5/10
	float const wire_radius(0.25*hr_radius);

	if (num_hangers > 0 && hanger_rod.get_sz_dim(!c.dim) > 10.0*wire_radius && has_hanger_model()) {
		room_object_t hanger(hanger_rod);
		hanger.type = TYPE_HANGER;
		set_cube_zvals(hanger, (hanger_rod.z1() - 0.09*window_vspacing), (hanger_rod.z2() + 2.0*wire_radius));
		hanger.expand_in_dim(c.dim, 0.09*window_vspacing); // set width
		float const pos_min(hanger_rod.d[!c.dim][0] + wire_radius), pos_max(hanger_rod.d[!c.dim][1] - wire_radius);
		float const pos_delta(pos_max - pos_min), slot_spacing(pos_delta/63.0);
		uint64_t slots_used(0); // divide the space into 64 slots, initially all empty
		bool const use_model(has_clothes_model());

		for (unsigned i = 0; i < num_hangers; ++i) { // since hangers are so narrow, we probably don't need to check for intersections
			unsigned slot_ix(0);
			bool found_slot(0);
			
			for (unsigned n = 0; n < 10; ++n) { // 10 attempts to find an unused slot
				slot_ix = rgen.rand()&63;
				uint64_t const slot_mask(uint64_t(1) << slot_ix);
				if (slots_used & slot_mask) continue;
				slots_used |= slot_mask;
				found_slot  = 1;
				break; // success
			}
			if (!found_slot) continue; // skip this hanger
			set_wall_width(hanger, (pos_min + slot_ix*slot_spacing), wire_radius, !c.dim);
			objects.push_back(hanger);
			
			if (rgen.rand_float() < 0.67) { // maybe add a shirt to the hanger
				objects.back().flags |= RO_FLAG_HANGING; // flag the hanger has having the shirt hanging on it
				room_object_t shirt(hanger); // or pants
				shirt.expand_by_xy(0.01*hanger.dz()); // expand slightly to avoid z-fighting with hanger
				shirt.expand_in_dim(c.dim, 0.045*c.dz()); // slightly wider
				shirt.z2() -= 0.55*hanger.dz(); // top
				shirt.z1() -= 0.3*c.dz(); // bottom
				unsigned const sflags(flags | RO_FLAG_HANGING);

				if (use_model) {
					objects.emplace_back(shirt, TYPE_CLOTHES, c.room_id, c.dim, c.dir, sflags, c.light_amt, SHAPE_CUBE, WHITE);
				}
				else {
					objects.emplace_back(shirt, TYPE_SHIRT, c.room_id, c.dim, c.dir, sflags, c.light_amt, SHAPE_CUBE, shirt_colors[rgen.rand()%NUM_SHIRT_COLORS]);
				}
			}
		} // for i
		objects[hanger_rod_ix].item_flags = uint16_t(objects.size() - hanger_rod_ix); // number of objects hanging on the hanger rod, including hangers and shirts
	}
}

void building_room_geom_t::expand_cabinet(room_object_t const &c) { // called on cabinets, counters, and kitchen sinks
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	rgen.rand_mix();
	vect_cube_t &cubes(get_temp_cubes());
	float const wall_thickness(0.04*c.dz()), light_amt(0.25*c.light_amt);
	cube_t interior(c), dishwasher;
	interior.expand_by(-wall_thickness);
	vector3d const c_sz(interior.get_size());
	
	if (c.type == TYPE_KSINK && get_dishwasher_for_ksink(c, dishwasher)) { // avoid placing objects that overlap the dishwasher
		dishwasher.d[c.dim][!c.dir] = c.d[c.dim][!c.dir]; // extend to the back of the cabinet
		dishwasher.expand_by_xy(wall_thickness);
		cubes.push_back(dishwasher);
	}
	unsigned const start_num_cubes(cubes.size()), flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);

	if ((c.type == TYPE_COUNTER || c.type == TYPE_KSINK) && rgen.rand_bool()) {
		float const tcan_height(c_sz.z*rgen.rand_uniform(0.35, 0.55)), tcan_radius(min(tcan_height/rgen.rand_uniform(1.6, 2.8), 0.4f*min(c_sz.x, c_sz.y)));
		cube_t tcan;
		gen_xy_pos_for_round_obj(tcan, interior, tcan_radius, tcan_height, 1.1*tcan_radius, rgen, 1); // place_at_z1=1
		room_object_t obj(tcan, TYPE_TCAN, c.room_id, c.dim, c.dir, flags, light_amt, (rgen.rand_bool() ? SHAPE_CYLIN : SHAPE_CUBE), tcan_colors[rgen.rand()%NUM_TCAN_COLORS]);
		add_if_not_intersecting(obj, expanded_objs, cubes);
	}
	// add boxes
	unsigned const num_boxes(rgen.rand()%4); // 0-3
	float const box_sz(0.3*c.get_sz_dim(c.dim));
	room_object_t cb(c);
	cb.light_amt = light_amt;
	add_boxes_to_space(cb, expanded_objs, interior, cubes, rgen, num_boxes, box_sz, 0.8*box_sz, 1.5*box_sz, 0, flags); // allow_crates=0
	// add paint cans (slightly smaller than normal)
	float const sz_scale(0.7*c_sz.z), pc_height(0.6*sz_scale), pc_radius(0.24*sz_scale);

	if (3*pc_radius < min(c_sz.x, c_sz.y)) { // have enough space for for paint cans
		unsigned const num_pcans(rgen.rand()%3); // 0-2

		for (unsigned n = 0; n < num_pcans; ++n) {
			cube_t pcan;
			gen_xy_pos_for_round_obj(pcan, interior, pc_radius, pc_height, 1.2*pc_radius, rgen, 1); // place_at_z1=1
			room_object_t obj(pcan, TYPE_PAINTCAN, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
			add_if_not_intersecting(obj, expanded_objs, cubes);
		}
	}
	// add plates
	unsigned const sz_ratio(round_fp(c_sz[!c.dim]/c_sz.z));
	unsigned const max_plates(0 + 1*sz_ratio), num_plates(rgen.rand() % max_plates); // wider cabinet has more plates
	float const plate_radius(min(sz_scale*rgen.rand_uniform(0.30, 0.35), 0.45f*min(c_sz.x, c_sz.y))), plate_height(0.1*plate_radius);

	for (unsigned n = 0; n < num_plates; ++n) {
		cube_t plate;
		gen_xy_pos_for_round_obj(plate, interior, plate_radius, plate_height, 1.1*plate_radius, rgen, 1); // place_at_z1=1
		room_object_t obj(plate, TYPE_PLATE, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
		if (!add_if_not_intersecting(obj, expanded_objs, cubes)) continue; // can't place the bottom plate
		unsigned const stack_height(1 + (rgen.rand()%5)); // 1-6

		for (unsigned s = 1; s < stack_height; ++s) {
			obj.translate_dim(2, plate_height); // shift up in z
			if (obj.z2() + plate_height > interior.z2()) break; // stack is too high, end it here
			expanded_objs.push_back(obj);
		}
	} // for n
	// add bottles
	unsigned const max_bottles(3 + 2*sz_ratio), num_bottles(rgen.rand() % max_bottles); // wider cabinet has more bottles

	for (unsigned n = 0; n < num_bottles; ++n) {
		float const bottle_height(sz_scale*rgen.rand_uniform(0.4, 0.65)), bottle_radius(sz_scale*rgen.rand_uniform(0.07, 0.1));
		if (min(c_sz.x, c_sz.y) < 3.0*bottle_radius) continue; // cabinet not wide/deep enough to add this bottle
		cube_t bottle;
		gen_xy_pos_for_round_obj(bottle, interior, bottle_radius, bottle_height, 1.5*bottle_radius, rgen, 1); // place_at_z1=1
		room_object_t obj(bottle, TYPE_BOTTLE, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN); // vertical
		obj.set_as_bottle(rgen.rand(), NUM_BOTTLE_TYPES-1, 1); // all bottle types, no_empty=1
		add_if_not_intersecting(obj, expanded_objs, cubes);
	}
	if (cubes.size() > start_num_cubes) {clear_static_small_vbos();} // some object was added
}

unsigned building_room_geom_t::get_shelves_for_object(room_object_t const &c, cube_t shelves[4]) {
	unsigned const num_shelves(2 + (c.room_id % 3)); // 2-4 shelves
	float const thickness(0.02*c.dz()), bracket_thickness(0.8*thickness), z_step(c.dz()/(num_shelves + 1)); // include a space at the bottom
	cube_t shelf(c);
	shelf.z2() = shelf.z1() + thickness; // set shelf thickness
	shelf.d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*bracket_thickness; // leave space behind the shelf for brackets

	for (unsigned s = 0; s < num_shelves; ++s) {
		shelf.translate_dim(2, z_step); // move up one step
		shelves[s] = shelf; // record for later use
	}
	return num_shelves;
}

void building_room_geom_t::get_shelf_objects(room_object_t const &c_in, cube_t const shelves[4], unsigned num_shelves, vector<room_object_t> &objects) {
	room_object_t c(c_in);
	c.flags |= RO_FLAG_WAS_EXP;
	bool const is_house(c.is_house());
	vector3d const c_sz(c.get_size());
	float const dz(c_sz.z), width(c_sz[c.dim]), thickness(0.02*dz), bracket_thickness(0.75*thickness);
	float const z_step(dz/(num_shelves + 1)), shelf_clearance(z_step - thickness - bracket_thickness), sz_scale(is_house ? 0.7 : 1.0), box_zscale(shelf_clearance*sz_scale);
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);

	for (unsigned s = 0; s < num_shelves; ++s) {
		cube_t const &S(shelves[s]);
		vect_cube_t &cubes(get_temp_cubes());
		room_object_t C(c);
		vector3d sz;
		// add crates/boxes
		unsigned const num_boxes(rgen.rand() % (is_house ? 8 : 13)); // 0-12
		cube_t bounds(S);
		bounds.z1() = S.z2(); // place on top of shelf
		add_boxes_to_space(c, objects, bounds, cubes, rgen, num_boxes, 0.42*width*sz_scale*(is_house ? 1.5 : 1.0), 0.4*box_zscale, 0.98*box_zscale, 1, c.flags); // allow_crates=1
		// add computers; what about monitors?
		bool const top_shelf(s+1 == num_shelves);
		float const h_val(0.21*1.1*dz), cheight(0.75*h_val), cwidth(0.44*cheight), cdepth(0.9*cheight); // fixed AR=0.44 to match the texture
		sz[ c.dim] = 0.5*cdepth;
		sz[!c.dim] = 0.5*cwidth;

		if (2.1*sz.x < c_sz.x && 2.1*sz.y < c_sz.y && (top_shelf || cheight < shelf_clearance)) { // if it fits in all dims
			unsigned const num_comps(rgen.rand() % (is_house ? 3 : 6)); // 0-5
			C.dim  = c.dim; C.dir = !c.dir; // reset dim, flip dir
			C.type = TYPE_COMPUTER;

			for (unsigned n = 0; n < num_comps; ++n) {
				gen_xy_pos_for_cube_obj(C, S, sz, cheight, rgen);
				add_if_not_intersecting(C, objects, cubes);
			}
		}
		// add keyboards
		float const kbd_hwidth(0.7*0.5*1.1*2.0*h_val), kbd_depth(0.6*kbd_hwidth), kbd_height(0.06*kbd_hwidth);
		sz[ c.dim] = 0.5*kbd_depth;
		sz[!c.dim] = kbd_hwidth;

		if (2.1*sz.x < c_sz.x && 2.1*sz.y < c_sz.y && (top_shelf || kbd_height < shelf_clearance)) { // if it fits in all dims
			unsigned const num_kbds(rgen.rand() % 5); // 0-4
			C.type = TYPE_KEYBOARD;

			for (unsigned n = 0; n < num_kbds; ++n) {
				gen_xy_pos_for_cube_obj(C, S, sz, kbd_height, rgen);
				C.dir = rgen.rand_bool();
				add_if_not_intersecting(C, objects, cubes);
			}
		}
		// add bottles
		unsigned const num_bottles(((rgen.rand()&3) == 0) ? 0 : (rgen.rand() % 9)); // 0-8, 75% chance
		C.dir  = C.dim = 0;
		C.type = TYPE_BOTTLE;

		for (unsigned n = 0; n < num_bottles; ++n) {
			// same as building_t::place_bottle_on_obj()
			float const bottle_height(z_step*rgen.rand_uniform(0.4, 0.7)), bottle_radius(z_step*rgen.rand_uniform(0.07, 0.11));
			if (min(c_sz.x, c_sz.y) < 5.0*bottle_radius) continue; // shelf not wide/deep enough to add this bottle
			gen_xy_pos_for_round_obj(C, S, bottle_radius, bottle_height, 2.0*bottle_radius, rgen);
			C.set_as_bottle(rgen.rand(), NUM_BOTTLE_TYPES-1, 1); // all bottle types, no_empty=1
			add_if_not_intersecting(C, objects, cubes);
		}
		// add paint cans
		float const pc_height(0.64*z_step), pc_radius(0.28*z_step), pc_edge_spacing(1.1*pc_radius);

		if (2.1*pc_edge_spacing < min(c_sz.x, c_sz.y)) { // shelf is wide/deep enough for paint cans
			unsigned const num_pcans(((rgen.rand()&3) == 0) ? 0 : (rgen.rand() % 7)); // 0-6, 75% chance
			C.color = WHITE;
			C.type  = TYPE_PAINTCAN;
			C.shape = SHAPE_CYLIN;

			for (unsigned n = 0; n < num_pcans; ++n) {
				gen_xy_pos_for_round_obj(C, S, pc_radius, pc_height, pc_edge_spacing, rgen);
				add_if_not_intersecting(C, objects, cubes);
			}
			C.shape = SHAPE_CUBE; // reset for next object type
		}
		// add spraypaint cans
		float const spcan_height(0.55*z_step), spcan_radius(0.17*spcan_height); // fixed size

		if (min(c_sz.x, c_sz.y) > 5.0*spcan_radius) { // add if shelf wide/deep enough
			unsigned const num_spcans(((rgen.rand()%5) < 3) ? (rgen.rand() % 4) : 0); // 0-3, 60% chance
			C.dir  = C.dim = 0;
			C.type = TYPE_SPRAYCAN;

			for (unsigned n = 0; n < num_spcans; ++n) {
				gen_xy_pos_for_round_obj(C, S, spcan_radius, spcan_height, 1.5*spcan_radius, rgen);
				C.color = spcan_colors[rgen.rand() % NUM_SPCAN_COLORS]; // random color
				add_if_not_intersecting(C, objects, cubes);
			}
		}
		// add large balls to houses
		float const ball_radius(0.048*1.1*dz); // 4.7 inches

		if (is_house && 2.1*ball_radius < min(c_sz.x, c_sz.y)) { // shelf is wide/deep enough for paint cans
			unsigned const num_balls(rgen.rand() % 3); // 0-2
			C.color = WHITE;
			C.type  = TYPE_LG_BALL;
			C.shape = SHAPE_SPHERE;

			for (unsigned n = 0; n < num_balls; ++n) {
				gen_xy_pos_for_round_obj(C, S, ball_radius, 2.0*ball_radius, ball_radius, rgen);
				C.item_flags = rgen.rand_bool(); // random type
				add_if_not_intersecting(C, objects, cubes);
			}
			C.shape      = SHAPE_CUBE; // reset for next object type
			C.item_flags = 0; // reset for next object type
		}
		// add tape rolls
		float const tape_radius(0.22*z_step), tape_height(TAPE_HEIGHT_TO_RADIUS*tape_radius); // fixed size

		if (min(c_sz.x, c_sz.y) > 3.0*tape_radius) { // add if shelf wide/deep enough
			unsigned const num_tapes(((rgen.rand()%4) < 3) ? (rgen.rand() % 4) : 0); // 0-3, 75% chance
			C.dir  = C.dim = 0;
			C.type = TYPE_TAPE;

			for (unsigned n = 0; n < num_tapes; ++n) {
				gen_xy_pos_for_round_obj(C, S, tape_radius, tape_height, 1.25*tape_radius, rgen);
				C.color = tape_colors[rgen.rand() % NUM_TAPE_COLORS]; // random color
				add_if_not_intersecting(C, objects, cubes);
			}
		}
	} // for s
}

void building_room_geom_t::expand_shelves(room_object_t const &c) {
	cube_t shelves[4]; // max number of shelves
	unsigned const num_shelves(get_shelves_for_object(c, shelves));
	get_shelf_objects(c, shelves, num_shelves, expanded_objs);
}

void building_room_geom_t::add_wine_rack_bottles(room_object_t const &c, vector<room_object_t> &objects) {
	float const height(c.dz()), width(c.get_sz_dim(!c.dim)), depth(c.get_sz_dim(c.dim)), shelf_thick(0.1*depth);
	unsigned const num_rows(max(1, round_fp(2.0*height/depth))), num_cols(max(1, round_fp(2.0*width/depth)));
	float const row_step((height - shelf_thick)/num_rows), col_step((width - shelf_thick)/num_cols);
	float const space_w(col_step - shelf_thick), space_h(row_step - shelf_thick), diameter(min(space_w, space_h) - shelf_thick);
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	room_object_t bottle(c);
	bottle.type   = TYPE_BOTTLE;
	bottle.dir   ^= 1;
	bottle.flags |= (RO_FLAG_WAS_EXP | RO_FLAG_NO_CONS); // wine bottles on wine racks are not consumable
	set_wall_width(bottle, (c.d[!c.dim][0] + 0.5f*(col_step + shelf_thick)), 0.5*diameter, !c.dim); // center in this dim

	for (unsigned i = 0; i < num_cols; ++i) { // columns/vertical
		bottle.z1() = c.z1() + shelf_thick; // rest on the top of the shelf
		bottle.z2() = bottle.z1() + diameter;

		for (unsigned j = 0; j < num_rows; ++j) { // rows/horizontal
			if (rgen.rand()%3) { // add a bottle 67% of the time; add_bottom=1
				bottle.set_as_bottle(3 + 64*rgen.rand_bool()); // always wine, but mixed cap color
				objects.push_back(bottle);
			}
			bottle.translate_dim(2, row_step); // translate in Z
		}
		bottle.translate_dim(!c.dim, col_step); // translate in !dim
	} // for i
}

void set_rand_pos_for_sz(cube_t &c, bool dim, float length, float width, rand_gen_t &rgen) {
	assert(c.get_sz_dim( dim) >= length);
	assert(c.get_sz_dim(!dim) >= width );
	c.d[ dim][0] = rgen.rand_uniform(c.d[ dim][0], (c.d[ dim][1] - length));
	c.d[ dim][1] = c.d[ dim][0] + length;
	c.d[!dim][0] = rgen.rand_uniform(c.d[!dim][0], (c.d[!dim][1] - width));
	c.d[!dim][1] = c.d[!dim][0] + width;
}

/*static*/ room_object_t building_room_geom_t::get_item_in_drawer(room_object_t const &c, cube_t const &drawer, unsigned drawer_ix) {
	assert(drawer_ix < 16);
	if (c.item_flags & (1U << drawer_ix)) {return room_object_t();} // item has been taken
	vector3d const sz(drawer.get_size()); // Note: drawer is the interior area
	rand_gen_t rgen;
	rgen.set_state((123*drawer_ix + 1), (456*c.room_id + 777*c.obj_id + 1));
	unsigned const type_ix(rgen.rand() % 11); // 0-10
	room_object_t obj; // starts as no item

	switch (type_ix) {
	case 0: // box
	{
		float const length(rgen.rand_uniform(0.4, 0.9)*sz[c.dim]), width(rgen.rand_uniform(0.4, 0.9)*sz[!c.dim]);
		obj = room_object_t(drawer, TYPE_BOX, c.room_id, rgen.rand_bool(), rgen.rand_bool());
		obj.color = gen_box_color(rgen);
		obj.z2()  = (obj.z1() + rgen.rand_uniform(0.4, 0.8)*sz.z);
		set_rand_pos_for_sz(obj, c.dim, length, width, rgen);
		break;
	}
	case 1: // paper
	{
		float const length(2.0*sz.z), width(0.77*length);

		if (length < 0.9*sz[c.dim] && width < 0.9*sz[!c.dim]) { // if it can fit
			obj = room_object_t(drawer, TYPE_PAPER, c.room_id, c.dim, c.dir);
			obj.obj_id = rgen.rand();
			obj.color  = paper_colors[rgen.rand()%NUM_PAPER_COLORS];
			obj.z2()   = (obj.z1() + 0.01*sz.z);
			set_rand_pos_for_sz(obj, c.dim, length, width, rgen);
		}
		break;
	}
	case 2: case 3: // pen/pencil/marker
	{
		unsigned const types[3] = {TYPE_PEN, TYPE_PENCIL, TYPE_MARKER}, type(types[rgen.rand()%3]);
		bool const dim(!c.dim); // always opposite orient of the drawer
		float const length(min(1.7f*sz.z, 0.9f*sz[dim])), diameter(((type == TYPE_MARKER) ? 0.08 : 0.036)*length);
		obj = room_object_t(drawer, type, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN);
		obj.color = ((obj.type == TYPE_MARKER) ? marker_colors[rgen.rand()&7] : ((obj.type == TYPE_PEN) ? pen_colors[rgen.rand()&3] : pencil_colors[rgen.rand()&1]));
		obj.z2()  = (obj.z1() + diameter);
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case 4: // book
	{
		float const length(rgen.rand_uniform(0.6, 0.9)*min(sz.x, sz.y)), width(rgen.rand_uniform(0.6, 1.0)*length);
		obj = room_object_t(drawer, TYPE_BOOK, c.room_id, !c.dim, (c.dir ^ c.dim));
		obj.obj_id = rgen.rand();
		obj.color  = book_colors[rgen.rand() % NUM_BOOK_COLORS];
		obj.z2()   = (obj.z1() + min(0.3f*width, rgen.rand_uniform(0.1, 0.35)*sz.z));
		set_rand_pos_for_sz(obj, c.dim, length, width, rgen);
		break;
	}
	case 5: // key (Note: aspect ratio of key depends on aspect ratio of door, but key model is always a constant aspect ratio)
	{
		obj = room_object_t(drawer, TYPE_KEY, c.room_id, rgen.rand_bool(), rgen.rand_bool());
		obj.expand_in_dim( obj.dim, -0.40*sz[ obj.dim]); // long  dim
		obj.expand_in_dim(!obj.dim, -0.46*sz[!obj.dim]); // short dim
		for (unsigned d = 0; d < 2; ++d) {obj.translate_dim(d, 0.35*sz[d]*rgen.rand_uniform(-1.0, 1.0));}
		obj.z2() = obj.z1() + 0.05*sz.z;
		break;
	}
	case 6: // bottle
	{
		bool const dim(c.dim ^ rgen.rand_bool() ^ 1); // random orient
		float const length(rgen.rand_uniform(0.7, 0.9)*min(1.8f*sz.z, min(sz.x, sz.y))), diameter(length*rgen.rand_uniform(0.26, 0.34));
		obj = room_object_t(drawer, TYPE_BOTTLE, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN);
		obj.set_as_bottle(rgen.rand()); // can be empty
		obj.z2() = (obj.z1() + diameter);
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case 7: // money
	{
		float const length(0.135*c.dz()), width(2.35*length);

		if (length < 0.9*sz[c.dim] && width < 0.9*sz[!c.dim]) { // if it can fit
			obj = room_object_t(drawer, TYPE_MONEY, c.room_id, c.dim, c.dir);
			obj.z2() = (obj.z1() + 0.01f*length*((rgen.rand()%20) + 1)); // 1-20 bills
			set_rand_pos_for_sz(obj, c.dim, length, width, rgen);
		}
		break;
	}
	case 8: // phone
	{
		bool const dim(c.dim ^ rgen.rand_bool()); // random orient
		float const length(0.3*c.dz()), width(0.45*length);

		if (length < 0.9*sz[dim] && width < 0.9*sz[!dim]) { // if it can fit
			unsigned const NUM_PHONE_COLORS = 7; // for the case
			colorRGBA const phone_colors[NUM_PHONE_COLORS] = {WHITE, GRAY, DK_GRAY, GRAY_BLACK, BLUE, RED, PINK};
			obj = room_object_t(drawer, TYPE_PHONE, c.room_id, dim, c.dir);
			obj.color = phone_colors[rgen.rand() % NUM_PHONE_COLORS];
			obj.z2()  = (obj.z1() + 0.045*length);
			set_rand_pos_for_sz(obj, dim, length, width, rgen);
		}
		break;
	}
	case 9: // spray paint can
	{
		bool const dim(c.dim ^ rgen.rand_bool() ^ 1); // random orient
		float const length(rgen.rand_uniform(0.8, 0.9)*min(1.8f*sz.z, min(sz.x, sz.y))), diameter(0.34*length);
		obj = room_object_t(drawer, TYPE_SPRAYCAN, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN, spcan_colors[rgen.rand() % NUM_SPCAN_COLORS]);
		obj.z2() = (obj.z1() + diameter);
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case 10: // tape roll
	{
		float const diameter(0.3*c.dz());

		if (diameter < 0.9*min(sz.x, sz.y)) { // if it can fit
			obj = room_object_t(drawer, TYPE_TAPE, c.room_id, 0, 0, 0, 1.0, SHAPE_CYLIN, tape_colors[rgen.rand() % NUM_TAPE_COLORS]); // dim/dir don't matter, so use 0
			obj.z2() = (obj.z1() + 0.5f*TAPE_HEIGHT_TO_RADIUS*diameter); // set height
			set_rand_pos_for_sz(obj, 0, diameter, diameter, rgen);
		}
		break;
	}
	} // end switch
	obj.flags    |= (RO_FLAG_WAS_EXP | RO_FLAG_NOCOLL);
	obj.light_amt = c.light_amt;
	return obj;
}

bool place_objects_in_box(cube_t const &box, vect_cube_t &obj_bcubes, float radius, float height) {
	if (height > box.dz()) return 0; // too tall to place in this box
	obj_bcubes.clear();
	unsigned const nx(floor(box.dx()/(2.0*radius))), ny(floor(box.dy()/(2.0*radius))); // truncate
	if (nx == 0 || ny == 0) return 0; // box is too small to place this object
	float const xspace(box.dx()/nx), yspace(box.dy()/ny); // >= width

	for (unsigned y = 0; y < ny; ++y) {
		float const ypos(box.y1() + (y+0.5)*yspace);

		for (unsigned x = 0; x < nx; ++x) {
			float const xpos(box.x1() + (x+0.5)*xspace);
			obj_bcubes.push_back(get_cube_height_radius(point(xpos, ypos, box.z1()), radius, height));
		}
	}
	return 1; // success
}

void building_t::add_box_contents(room_object_t const &box) {
	rand_gen_t rgen;
	box.set_rand_gen_state(rgen);
	rgen.rand_mix();
	cube_t c(box);
	c.expand_by(-0.01*box.get_size()); // shrink to interior area
	vector3d const sz(c.get_size());
	bool const dim(sz.x < sz.y); // long dim
	float const light_amt(box.light_amt), floor_spacing(get_window_vspace()), base_height(0.2*floor_spacing); // avg shelf height
	unsigned const flags(RO_FLAG_WAS_EXP | RO_FLAG_NOCOLL);
	uint8_t const room_id(box.room_id);
	auto &objs(interior->room_geom->expanded_objs);
	vect_cube_t obj_bcubes;

	// Note: the code below may invalidate the reference to box, so we can't use it after this point
	for (unsigned n = 0; n < 10; ++n) { // make up to 10 attempts at placing valid item(s) in this box
		unsigned const obj_type((n == 9) ? 0 : (rgen.rand()%7)); // {book, bottles, ball, paint can, spraypaint, toilet paper, tape, [box]}; place book on last iteration

		if (obj_type == 0) { // books; can always be placed
			unsigned const num_books(1 + (rgen.rand()&3)); // 1-4 books
			float cur_zval(c.z1());

			for (unsigned n = 0; n < num_books; ++n) {
				float const length(rgen.rand_uniform(0.7, 0.95)*min(sz[dim], 2.0f*sz[!dim])), width(min(rgen.rand_uniform(0.6, 1.0)*length, 0.95f*sz[!dim]));
				room_object_t obj(c, TYPE_BOOK, room_id, !dim, rgen.rand_bool(), flags, light_amt);
				obj.obj_id = rgen.rand();
				obj.color  = book_colors[rgen.rand() % NUM_BOOK_COLORS];
				set_cube_zvals(obj, cur_zval, (cur_zval + min(0.3f*width, rgen.rand_uniform(0.1, 0.2)*sz.z)));
				if (obj.z2() > c.z2()) break; // book doesn't fit - the stack is too tall; can't fail on the first book
				set_rand_pos_for_sz(obj, dim, length, width, rgen);
				objs.push_back(obj);
				cur_zval = obj.z2();
			} // for n
		}
		else if (obj_type == 1) { // bottles; not comsumable, as this would make things too easy for the player
			float const height(base_height*rgen.rand_uniform(0.4, 0.7)), radius(base_height*rgen.rand_uniform(0.07, 0.11));
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			unsigned const bottle_id(rgen.rand()); // same type for all bottles

			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_BOTTLE, room_id, 0, 0, (flags | RO_FLAG_NO_CONS), light_amt, SHAPE_CYLIN);
				objs.back().set_as_bottle(bottle_id, 3, 1); // 0-3; excludes poison; no_empty=1
			}
		}
		else if (obj_type == 2) { // ball - only 1
			if (!is_house) continue; // balls in houses only, not office buildings
			float const radius(0.048*floor_spacing);
			if (c.min_len() < 2.0*radius) continue; // can't fit any of this item
			point const center(c.xc(), c.yc(), (c.z1() + radius));
			cube_t ball; ball.set_from_sphere(center, radius);
			objs.emplace_back(ball, TYPE_LG_BALL, room_id, 0, 0, flags, light_amt, SHAPE_SPHERE);
			objs.back().item_flags = rgen.rand_bool(); // selects ball type
		}
		else if (obj_type == 3) { // paint cans
			float const height(0.64*base_height), radius(0.28*base_height);
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			
			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_PAINTCAN, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
			}
		}
		else if (obj_type == 4) { // spraypaint cans
			float const height(0.55*base_height), radius(0.17*height);
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			unsigned color_ix(rgen.rand()); // random starting color
			
			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_SPRAYCAN, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
				objs.back().color = spcan_colors[(color_ix++) % NUM_SPCAN_COLORS];
			}
		}
		else if (obj_type == 5) { // toilet paper rolls
			float const height(0.35*0.18*floor_spacing), radius(0.4*height);
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			
			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_TPROLL, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
			}
		}
		else if (obj_type == 6) { // rolls of tape
			float height(0.032*floor_spacing), radius(height/TAPE_HEIGHT_TO_RADIUS);

			for (unsigned m = 0; m < 2; ++m) {
				if (2.0*radius < 0.95*min(c.dx(), c.dy())) break; // size is okay
				height *= 0.9; radius *= 0.9; // can't fit any, try making it smaller
			}
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			unsigned color_ix(rgen.rand()); // random color

			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				// dim/dir don't matter, so use 0; flag as RO_FLAG_IN_CLOSET so that we know that we don't need to draw the bottom (not in a drawer)
				objs.emplace_back(*i, TYPE_TAPE, room_id, 0, 0, (flags | RO_FLAG_IN_CLOSET), light_amt, SHAPE_CYLIN);
				objs.back().color = tape_colors[color_ix % NUM_TAPE_COLORS];
			}
		}
		else if (obj_type == 7) { // nested box (not currently enabled)
			cube_t box2(box);
			box2.expand_by(-0.05*box.get_size()); // shrink by 10% in all dims
			objs.emplace_back(box2, TYPE_BOX, room_id, box.dim, box.dir, (flags | (box.flags & ~RO_FLAG_OPEN)), light_amt);
			objs.back().obj_id += rgen.rand();
		}
		else {continue;} // empty box?
		interior->room_geom->clear_static_small_vbos();
		break; // if we got here, something was placed in the box
	} // for n
}

void building_room_geom_t::expand_object(room_object_t &c) {
	if (c.flags & RO_FLAG_EXPANDED) return; // already expanded
	switch (c.type) {
	case TYPE_CLOSET:    expand_closet   (c); break;
	case TYPE_SHELVES:   expand_shelves  (c); break;
	case TYPE_WINE_RACK: expand_wine_rack(c); break;
	case TYPE_CABINET: case TYPE_COUNTER: case TYPE_KSINK: expand_cabinet(c); break;
	default: assert(0); // not a supported expand type
	}
	c.flags |= RO_FLAG_EXPANDED; // flag as expanded
}

