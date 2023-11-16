// 3D World - Building Interior Room Geometry Drawing
// by Frank Gennari 7/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

float const TAPE_HEIGHT_TO_RADIUS = 0.6;

extern object_model_loader_t building_obj_model_loader;

float get_lamp_width_scale();
vect_cube_t &get_temp_cubes();
bool get_dishwasher_for_ksink(room_object_t const &c, cube_t &dishwasher);
void set_wall_width(cube_t &wall, float pos, float half_thick, unsigned dim);
float get_med_cab_wall_thickness(room_object_t const &c);

void resize_model_cube_xy(cube_t &cube, float dim_pos, float not_dim_pos, unsigned id, bool dim) {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(id)); // L, W, H
	set_wall_width(cube, dim_pos,     0.5*cube.dz()*sz.x/sz.z,  dim); // length
	set_wall_width(cube, not_dim_pos, 0.5*cube.dz()*sz.y/sz.z, !dim); // width
}

bool is_shirt_model     (room_object_t const &obj) {return building_obj_model_loader.model_filename_contains(obj.get_model_id(), "shirt", "Shirt");}
bool is_pants_model     (room_object_t const &obj) {return building_obj_model_loader.model_filename_contains(obj.get_model_id(), "pants", "Pants");}
bool is_bar_hanger_model(room_object_t const &obj) {return building_obj_model_loader.model_filename_contains(obj.get_model_id(), "bar hanger", "Bar Hanger");}

bool add_if_not_intersecting(room_object_t const &obj, vect_room_object_t &objects, vect_cube_t &cubes) {
	if (has_bcube_int(obj, cubes)) return 0;
	objects.push_back(obj);
	cubes.push_back(obj);
	return 1;
}
point gen_xy_pos_in_area(cube_t const &S, vector3d const &sz, rand_gen_t &rgen, float zval) {
	point center(0.0, 0.0, zval);
	
	for (unsigned d = 0; d < 2; ++d) {
		float const lo(S.d[d][0]+sz[d]), hi(S.d[d][1]-sz[d]);
		assert(lo <= hi);
		center[d] = rgen.rand_uniform(lo, hi); // randomly placed within the specified bounds
	}
	return center;
}
point gen_xy_pos_in_area(cube_t const &S, float radius, rand_gen_t &rgen, float zval) {
	return gen_xy_pos_in_area(S, vector3d(radius, radius, 0.0), rgen, zval);
}
void gen_xy_pos_for_cube_obj(cube_t &C, cube_t const &S, vector3d const &sz, float height, rand_gen_t &rgen) {
	C.set_from_point(gen_xy_pos_in_area(S, sz, rgen));
	set_cube_zvals(C, S.z2(), S.z2()+height);
	C.expand_by_xy(sz);
}
void gen_xy_pos_for_round_obj(cube_t &C, cube_t const &S, float radius, float height, float spacing, rand_gen_t &rgen, bool place_at_z1) {
	C.set_from_sphere(gen_xy_pos_in_area(S, spacing, rgen), radius); // place at least spacing from edge
	float const place_z(place_at_z1 ? S.z1() : S.z2());
	set_cube_zvals(C, place_z, place_z+height);
}

void add_boxes_to_space(room_object_t const &c, vect_room_object_t &objects, cube_t const &bounds, vect_cube_t &cubes, rand_gen_t &rgen,
	unsigned num_boxes, float xy_scale, float hmin, float hmax, bool allow_crates, unsigned flags)
{
	float const bounds_sz[2] = {bounds.dx(), bounds.dy()};
	room_object_t C(c);
	C.flags = flags; // Note: also clears open flag
	vector3d sz;

	for (unsigned n = 0; n < num_boxes; ++n) {
		for (unsigned d = 0; d < 2; ++d) {
			sz[d] = min(xy_scale*rgen.rand_uniform(0.5, 1.0), 0.99f*0.5f*bounds_sz[d]); // x,y half width; clamp to slightly smaller than bounds to avoid an assert
		}
		C.set_from_point(gen_xy_pos_in_area(bounds, sz, rgen)); // randomly placed within the bounds of the closet
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

bool add_obj_to_closet(room_object_t const &c, cube_t const &interior, vect_room_object_t &objects, vect_cube_t &cubes,
	rand_gen_t &rgen, vector3d const &sz, unsigned obj_type, unsigned flags, room_obj_shape shape=SHAPE_CUBE)
{
	for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts
		cube_t obj;
		obj.set_from_point(gen_xy_pos_in_area(interior, sz, rgen));
		set_cube_zvals(obj, interior.z1(), (interior.z1() + sz.z));
		obj.expand_by_xy(sz);

		if (!has_bcube_int(obj, cubes)) { // check for intersection with boxes
			objects.emplace_back(obj, obj_type, c.room_id, c.dim, c.dir, flags, c.light_amt, shape);
			cubes.push_back(obj);
			return 1; // success
		}
	} // for n
	return 0; // failed
}

bool try_add_lamp(cube_t const &place_area, float floor_spacing, unsigned room_id, unsigned flags, float light_amt,
	vect_cube_t &cubes, vect_room_object_t &objects, rand_gen_t &rgen)
{
	float const height(0.25*floor_spacing), width(height*get_lamp_width_scale()), radius(0.5*width);
	if (width == 0.0 || width > 0.9*min(place_area.dx(), place_area.dy())) return 0; // check if lamp model is valid and lamp fits in closet
		
	for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts
		point const center(gen_xy_pos_in_area(place_area, radius, rgen, place_area.z1()));
		cube_t lamp(get_cube_height_radius(center, radius, height));

		if (!has_bcube_int(lamp, cubes)) { // check for intersection with other objects
			objects.emplace_back(lamp, TYPE_LAMP, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN, lamp_colors[rgen.rand()%NUM_LAMP_COLORS]);
			cubes.push_back(lamp);
			return 1;
		}
	} // for n
	return 0;
}

rand_gen_t room_object_t::create_rgen() const {
	rand_gen_t rgen;
	rgen.set_state(obj_id+1, room_id+1);
	rgen.rand_mix(); // optional, but improves randomness
	return rgen;
}

void building_room_geom_t::add_closet_objects(room_object_t const &c, vect_room_object_t &objects) {
	rand_gen_t rgen(c.create_rgen());
	cube_t ccubes[5]; // only used to get interior space
	get_closet_cubes(c, ccubes);
	cube_t const interior(get_closet_interior_space(c, ccubes));
	float const depth(interior.get_sz_dim(c.dim)), box_sz(0.25*depth), window_vspacing(c.dz()*(1.0 + FLOOR_THICK_VAL_HOUSE));
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);
	unsigned const num_boxes((rgen.rand()%3) + (rgen.rand()%4)); // 0-5
	vect_cube_t &cubes(get_temp_cubes());
	add_boxes_to_space(c, objects, interior, cubes, rgen, num_boxes, box_sz, 0.8*box_sz, 1.5*box_sz, 0, flags); // allow_crates=0

	if (!c.is_small_closet()) { // larger closets have more random items
		vector3d sz;

		if (rgen.rand_bool()) { // maybe add a lamp in the closet
			try_add_lamp(interior, window_vspacing, c.room_id, flags, c.light_amt, cubes, objects, rgen);
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

	// add hangers and hanging clothes
	unsigned const num_hangers(c.is_small_closet() ? ((rgen.rand()%7) + 2) : ((rgen.rand()%13) + 4)); // 2-8 / 4-16
	float const wire_radius(0.25*hr_radius);

	if (num_hangers > 0 && hanger_rod.get_sz_dim(!c.dim) > 10.0*wire_radius && building_obj_model_loader.is_model_valid(OBJ_MODEL_HANGER)) {
		room_object_t hanger(hanger_rod);
		hanger.type       = TYPE_HANGER;
		hanger.item_flags = rgen.rand(); // choose a random hanger sub_model_id per-closet
		set_cube_zvals(hanger, (hanger_rod.z1() - 0.07*window_vspacing), (hanger_rod.z2() + 1.0*wire_radius));
		unsigned const hid(hanger.get_model_id());
		float const edge_spacing(max(4.0f*hr_radius, 0.1f*depth));
		float const pos_min(hanger_rod.d[!c.dim][0] + edge_spacing), pos_max(hanger_rod.d[!c.dim][1] - edge_spacing);
		float const pos_delta(pos_max - pos_min), slot_spacing(pos_delta/63.0);
		uint64_t slots_used(0); // divide the space into 64 slots, initially all empty
		bool const use_model(building_obj_model_loader.is_model_valid(OBJ_MODEL_CLOTHES));
		bool const mix_hangers(rgen.rand_bool());

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
			resize_model_cube_xy(hanger, hanger.get_center_dim(c.dim), (pos_min + slot_ix*slot_spacing), hid, c.dim);
			expand_to_nonzero_area(hanger, 0.1*hr_radius, !c.dim);
			objects.push_back(hanger);
			
			if (use_model && rgen.rand_float() < 0.8) { // maybe add clothing to the hanger
				objects.back().flags |= RO_FLAG_HANGING; // flag the hanger has having the shirt hanging on it
				room_object_t clothes(hanger, TYPE_CLOTHES, c.room_id, c.dim, c.dir, (flags | RO_FLAG_HANGING), c.light_amt, SHAPE_CUBE, WHITE);
				clothes.z2() -= 0.55*hanger.dz(); // top
				clothes.z1() -= 0.3*c.dz(); // bottom
				clothes.item_flags = rgen.rand(); // choose a random clothing sub_model_id
				if (is_pants_model(clothes) && is_bar_hanger_model(hanger)) {++clothes.item_flags;} // hack to avoid placing pants on bar hangers
				unsigned const cid(clothes.get_model_id());
				float const scale(building_obj_model_loader.get_model_scale(cid));
				if (scale != 0.0) {set_wall_width(clothes, clothes.zc(), 0.5*scale*clothes.dz(), 2);} // rescale zvals around the center
				resize_model_cube_xy(clothes, clothes.get_center_dim(c.dim), clothes.get_center_dim(!c.dim), cid, c.dim);
				bool skip(0);

				for (auto m = objects.begin()+hanger_rod_ix+1; m != objects.end(); ++m) { // skip if this object intersects a previous hanging clothing item
					if (m->type == TYPE_CLOTHES && m->intersects(clothes)) {skip = 1; break;}
				}
				if (skip) continue;

				if (is_shirt_model(clothes) && rgen.rand_float() < 0.67) { // 67% of shirts and randomly colored rather than colored + textured with the model
					clothes.color  = shirt_colors[rgen.rand()%NUM_SHIRT_COLORS];
					clothes.flags |= RO_FLAG_UNTEXTURED;
				}
				objects.push_back(clothes);
			}
			if (mix_hangers && (rgen.rand()&7) == 0) {hanger.item_flags = rgen.rand();} // switch to a new hanger type occasionally
		} // for i
		objects[hanger_rod_ix].item_flags = uint16_t(objects.size() - hanger_rod_ix); // number of objects hanging on the hanger rod, including hangers and shirts
	}
}

void building_room_geom_t::expand_cabinet(room_object_t const &c) { // called on cabinets, counters, and kitchen sinks
	rand_gen_t rgen(c.create_rgen());
	vect_cube_t &cubes(get_temp_cubes());
	float const wall_thickness(0.04*c.dz()), light_amt(0.25*c.light_amt);
	cube_t interior(c), dishwasher;
	interior.expand_by(-wall_thickness);
	vector3d const c_sz(interior.get_size());
	bool const is_counter(c.type == TYPE_COUNTER || c.type == TYPE_KSINK);
	
	if (c.type == TYPE_KSINK && get_dishwasher_for_ksink(c, dishwasher)) { // avoid placing objects that overlap the dishwasher
		dishwasher.d[c.dim][!c.dir] = c.d[c.dim][!c.dir]; // extend to the back of the cabinet
		dishwasher.expand_by_xy(wall_thickness);
		cubes.push_back(dishwasher);
	}
	unsigned const start_num_cubes(cubes.size()), flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);

	if (is_counter && rgen.rand_bool()) {
		float const tcan_height(c_sz.z*rgen.rand_uniform(0.35, 0.55)), tcan_radius(min(tcan_height/rgen.rand_uniform(1.6, 2.8), 0.4f*min(c_sz.x, c_sz.y)));
		cube_t tcan;
		gen_xy_pos_for_round_obj(tcan, interior, tcan_radius, tcan_height, 1.1*tcan_radius, rgen, 1); // place_at_z1=1
		room_object_t obj(tcan, TYPE_TCAN, c.room_id, c.dim, c.dir, flags, light_amt, (rgen.rand_bool() ? SHAPE_CYLIN : SHAPE_CUBE), tcan_colors[rgen.rand()%NUM_TCAN_COLORS]);
		add_if_not_intersecting(obj, expanded_objs, cubes);
	}
	// maybe add fire extinguisher
	if (is_counter && building_obj_model_loader.is_model_valid(OBJ_MODEL_FIRE_EXT) && rgen.rand_float() < 0.3) { // 30% of the time
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FIRE_EXT)); // D, W, H
		float const fe_height(0.38*c.dz()), fe_radius(fe_height*(0.5*(sz.x + sz.y)/sz.z));

		if (2.2*fe_radius < min(c_sz.x, c_sz.y)) {
			cube_t fire_ext;
			gen_xy_pos_for_round_obj(fire_ext, interior, fe_radius, fe_height, 1.05*fe_radius, rgen, 1); // place_at_z1=1
			room_object_t obj(fire_ext, TYPE_FIRE_EXT, c.room_id, rgen.rand_bool(), rgen.rand_bool(), flags, light_amt, SHAPE_CYLIN);
			add_if_not_intersecting(obj, expanded_objs, cubes);
		}
	}
	// add boxes
	unsigned const num_boxes(rgen.rand()%4); // 0-3
	float const box_sz(0.3*c.get_length());
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
	// add pans
	unsigned const num_pans(rgen.rand()%3); // 0-2
	float const pan_radius(min(sz_scale*rgen.rand_uniform(0.40, 0.45), 0.45f*min(c_sz.x, c_sz.y))), pan_height(rgen.rand_uniform(0.4, 0.5)*pan_radius);

	for (unsigned n = 0; n < num_pans; ++n) {
		cube_t pan;
		gen_xy_pos_for_round_obj(pan, interior, pan_radius, pan_height, 1.1*pan_radius, rgen, 1); // place_at_z1=1
		pan.translate_dim(2, 0.01*pan_height); // shift in +z to prevent z-fighting with the bottom of the pan
		bool const dir(pan.get_center_dim(!c.dim) < c.get_center_dim(!c.dim)); // point toward the cabinet center to avoid the handle clipping through the side/end
		room_object_t obj(pan, TYPE_PAN, c.room_id, c.dim, dir, flags, light_amt, SHAPE_CYLIN, GRAY_BLACK);
		// Note: the pan's handle extends outside its bcube and may clip through other objects, but this isn't very noticeable when viewed from normal head height
		add_if_not_intersecting(obj, expanded_objs, cubes);
	} // for n
	// add bottles
	unsigned const max_bottles(3 + 2*sz_ratio), num_bottles(rgen.rand() % max_bottles); // wider cabinet has more bottles

	for (unsigned n = 0; n < num_bottles; ++n) {
		float const bottle_height(sz_scale*rgen.rand_uniform(0.4, 0.65)), bottle_radius(sz_scale*rgen.rand_uniform(0.07, 0.1));
		if (min(c_sz.x, c_sz.y) < 3.0*bottle_radius) continue; // cabinet not wide/deep enough to add this bottle
		cube_t bottle;
		gen_xy_pos_for_round_obj(bottle, interior, bottle_radius, bottle_height, 1.5*bottle_radius, rgen, 1); // place_at_z1=1
		room_object_t obj(bottle, TYPE_BOTTLE, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN); // vertical
		bool const allow_medicine(rgen.rand_bool()); // medicine is half as common
		obj.set_as_bottle(rgen.rand(), (allow_medicine ? (unsigned)NUM_BOTTLE_TYPES : (unsigned)BOTTLE_TYPE_MEDS)-1, 1); // all bottle types, no_empty=1
		add_if_not_intersecting(obj, expanded_objs, cubes);
	}
	if (cubes.size() > start_num_cubes) {invalidate_small_geom();} // some object was added
}

void building_room_geom_t::expand_med_cab(room_object_t const &c) { // aka house "mirrors"
	rand_gen_t rgen(c.create_rgen());
	// add medicine bottle
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);
	float const wall_thickness(get_med_cab_wall_thickness(c));
	float const height(0.33*c.dz()*rgen.rand_uniform(0.75, 1.0)), radius(0.35*c.get_sz_dim(c.dim)*rgen.rand_uniform(0.75, 1.0));
	cube_t interior(c), bottle;
	interior.expand_by(-wall_thickness);
	if (rgen.rand_bool()) {interior.z1() += 0.5*(interior.dz() + wall_thickness);} // place on middle shelf half the time
	gen_xy_pos_for_round_obj(bottle, interior, radius, height, 0.1*radius, rgen, 1); // place_at_z1=1
	room_object_t obj(bottle, TYPE_BOTTLE, c.room_id, 0, 0, flags, c.light_amt, SHAPE_CYLIN); // vertical
	obj.set_as_bottle(BOTTLE_TYPE_MEDS, BOTTLE_TYPE_MEDS, 1); // medicine, no_empty=1
	expanded_objs.push_back(obj);
	invalidate_small_geom();
}

// Note: see building_room_geom_t::add_breaker_panel()
void building_room_geom_t::expand_breaker_panel(room_object_t const &c, bool has_elevator, bool has_parking_garage) {
	float const box_width(c.get_sz_dim(!c.dim)), box_depth(c.get_sz_dim(c.dim)), box_height(c.dz());
	float const thickness(0.25*box_depth), dir_sign(c.dir ? -1.0 : 1.0);
	cube_t breakers(c);
	breakers.d[c.dim][0] = breakers.d[c.dim][1] = breakers.d[c.dim][!c.dir] - dir_sign*thickness; // move both edges to front of panel
	breakers.d[c.dim][!c.dir] += dir_sign*0.5*thickness; // expand out for correct width
	breakers.expand_in_dim(!c.dim, -0.2*box_width); // shrink the width
	breakers.expand_in_dim(2, -0.1*box_height); // shrink vertically
	float const width(breakers.get_sz_dim(!c.dim)), height(breakers.dz());
	float const breaker_height(0.1*box_depth + 0.025*box_width + 0.05*box_height), breaker_width(2.4*breaker_height);
	unsigned const num_cols(max(1, round_fp(width/breaker_width))), num_rows(max(1, round_fp(height/breaker_height))), num_breakers(num_cols*num_rows);
	float const breaker_dc(width/num_cols), breaker_dr(height/num_rows), breaker_cs(0.1*breaker_dc), breaker_rs(0.06*breaker_dr);
	cube_t breaker(breakers);
	colorRGBA const color(0.03, 0.015, 0.01); // black-brown

	for (unsigned C = 0; C < num_cols; ++C) {
		float const cv(breakers.d[!c.dim][0] + C*breaker_dc);
		breaker.d[!c.dim][0] = cv + breaker_cs;
		breaker.d[!c.dim][1] = cv - breaker_cs + breaker_dc;

		for (unsigned r = 0; r < num_rows; ++r) {
			float const rv(breakers.z1() + r*breaker_dr);
			breaker.z1() = rv + breaker_rs;
			breaker.z2() = rv - breaker_rs + breaker_dr;
			room_object_t obj(breaker, TYPE_BREAKER, c.room_id, c.dim, c.dir, (RO_FLAG_NOCOLL | RO_FLAG_OPEN), c.light_amt, SHAPE_CUBE, color); // starts on
			obj.obj_id     = C*num_rows + r; // store breaker index
			obj.item_flags = num_breakers; // store the number of breakers
			expanded_objs.push_back(obj);
			// add labels
			cube_t label(breaker);
			string label_text;

			if (has_elevator && C == 0 && r == 0) { // first breaker - add elevator label
				label.translate_dim(!c.dim, -0.75*breaker_dc);
				label_text = "Elevator";
			}
			else if (has_parking_garage && C == num_cols-1 && r == num_rows-1) { // last breaker - add parking garage label
				label.translate_dim(!c.dim, 0.75*breaker_dc);
				label_text = "Garage";
			}
			if (!label_text.empty()) {
				label.d[c.dim][!c.dir] = label.d[c.dim][c.dir] + dir_sign*0.05*thickness;
				label.expand_in_dim(!c.dim, -0.2*breaker_dc); // small shrink
				label.expand_in_dim(2,      -0.2*breaker_dr); // small shrink
				expanded_objs.emplace_back(label, TYPE_SIGN, c.room_id, c.dim, !c.dir, RO_FLAG_NOCOLL, c.light_amt, SHAPE_CUBE, BLACK);
				expanded_objs.back().obj_id = register_sign_text(label_text);
			}
		}
	} // for C
	invalidate_small_geom();
}

void building_room_geom_t::expand_dishwasher(room_object_t &c, cube_t const &dishwasher) {
	unsigned const expanded_objs_start(expanded_objs.size());
	c.item_flags = expanded_objs_start; // record the location where we'll add our objects
	//rand_gen_t rgen(c.create_rgen());
	// add TYPE_PLATE, and TYPE_CUP, and TYPE_SILVER
	unsigned const num_plates = 4; // random?
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);
	uint16_t obj_place_index(0);
	// see dishwasher drawing code in building_room_geom_t::add_counter()
	float const dz(c.dz()), width(dishwasher.get_sz_dim(!c.dim)), depth(c.get_depth());
	float const wall_width(0.1*dz), door_thickness(0.035*depth), tray_sz(min(dz, width)), dir_sign(c.dir ? 1.0 : -1.0);
	float const plate_radius(0.3*tray_sz), plate_thickness(0.1*plate_radius), plate_spacing(2.0*plate_thickness), plate_dim_offset(-0.5*num_plates*plate_spacing);
	point tray_center;
	tray_center[ c.dim] = dishwasher.d[c.dim][c.dir] + 0.5*dir_sign*dz;
	tray_center[!c.dim] = dishwasher.get_center_dim(!c.dim);
	tray_center.z       = dishwasher.z1() + door_thickness + wall_width + 0.85*plate_radius; // top of tray
	cube_t plate;
	bool const plate_dim(!c.dim);
	plate.set_from_point(tray_center);
	plate.expand_in_dim( plate_dim, 0.5*plate_thickness); // set thickness
	plate.expand_in_dim(!plate_dim, plate_radius);
	plate.expand_in_dim(2,          plate_radius); // expand in Z
	plate.translate_dim(plate_dim, plate_dim_offset);

	for (unsigned n = 0; n < num_plates; ++n, ++obj_place_index) {
		if (!(c.taken_level & (1<<obj_place_index))) { // plate not yet taken
			expanded_objs.emplace_back(plate, TYPE_PLATE, c.room_id, plate_dim, (c.dir ^ c.dim), flags, c.light_amt, SHAPE_CYLIN, WHITE);
			expanded_objs.back().obj_id = obj_place_index;
		}
		plate.translate_dim(plate_dim, plate_spacing);
	}
	if (expanded_objs.size() > expanded_objs_start) {invalidate_small_geom();} // if something was added
}
void building_room_geom_t::unexpand_dishwasher(room_object_t &c, cube_t const &dishwasher) {
	cube_t door_region(dishwasher);
	door_region.d[c.dim][c.dir] += (c.dir ? 1.0 : -1.0)*c.dz();
	unsigned num_rem(0);
	c.taken_level = 0xFF; // assume everything has been taken unless we find objects to reclaim

	for (unsigned i = c.item_flags; i < expanded_objs.size(); ++i) { // search starting at our expanded_objs insertion point
		room_object_t &obj(expanded_objs[i]);
		if (!door_region.contains_pt(obj.get_cube_center())) break; // break when outside our object range
		if (/*obj.type == TYPE_BLOCKER*/obj.type != TYPE_PLATE && obj.type != TYPE_CUP && obj.type != TYPE_SILVER) continue; // not a dishwasher item
		obj.remove();
		c.taken_level &= ~(1 << (i - c.item_flags)); // mark this object as untaken
		++num_rem;
	}
	// pop removed objects from the end, in case the player repeatedly opens and closes the same dishwasher
	while (!expanded_objs.empty() && expanded_objs.back().type == TYPE_BLOCKER) {expanded_objs.pop_back();}
	c.item_flags = 0xFFFF; // set to some illegal value
	if (num_rem > 0) {invalidate_small_geom();}
}

unsigned building_room_geom_t::get_shelves_for_object(room_object_t const &c, cube_t shelves[4]) {
	unsigned const num_shelves(c.get_num_shelves()); // 2-4 shelves
	float const thickness(0.02*c.dz()), bracket_thickness(0.8*thickness), z_step(c.dz()/(num_shelves + 1)); // include a space at the bottom
	cube_t shelf(c);
	shelf.z2() = shelf.z1() + thickness; // set shelf thickness
	shelf.d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*bracket_thickness; // leave space behind the shelf for brackets

	for (unsigned s = 0; s < num_shelves; ++s) {
		shelf.translate_dim(2, z_step); // move up one step; done first, so bottom shelf is above the floor
		shelves[s] = shelf; // record for later use
	}
	return num_shelves;
}

void building_room_geom_t::get_shelf_objects(room_object_t const &c_in, cube_t const shelves[4], unsigned num_shelves, vect_room_object_t &objects) {
	room_object_t c(c_in);
	c.flags |= RO_FLAG_WAS_EXP;
	bool const is_house(c.is_house());
	vector3d const c_sz(c.get_size());
	float const dz(c_sz.z), width(c_sz[c.dim]), thickness(0.02*dz), bracket_thickness(0.75*thickness), floor_spacing(1.1*dz);
	float const z_step(dz/(num_shelves + 1)), shelf_clearance(z_step - thickness - bracket_thickness), sz_scale(is_house ? 0.7 : 1.0), box_zscale(shelf_clearance*sz_scale);
	rand_gen_t rgen(c.create_rgen());

	// Note: this function doesn't support placement of objects drawn as 3D models such as fire extinguishers
	for (unsigned s = 0; s < num_shelves; ++s) {
		cube_t const &S(shelves[s]);
		vect_cube_t &cubes(get_temp_cubes());
		room_object_t C(c);
		// add crates/boxes
		unsigned const num_boxes(rgen.rand() % (is_house ? 8 : 13)); // 0-12
		cube_t bounds(S);
		bounds.z1() = S.z2(); // place on top of shelf
		add_boxes_to_space(c, objects, bounds, cubes, rgen, num_boxes, 0.42*width*sz_scale*(is_house ? 1.5 : 1.0), 0.4*box_zscale, 0.98*box_zscale, 1, c.flags); // allow_crates=1
		// add computers; what about monitors?
		bool const top_shelf(s+1 == num_shelves);
		float const h_val(0.21*floor_spacing), cheight(0.75*h_val), cwidth(0.44*cheight), cdepth(0.9*cheight); // fixed AR=0.44 to match the texture
		vector3d sz;
		sz[ c.dim] = 0.5*cdepth;
		sz[!c.dim] = 0.5*cwidth;

		if (2.1*sz.x < c_sz.x && 2.1*sz.y < c_sz.y && (top_shelf || cheight < shelf_clearance)) { // if it fits in all dims
			unsigned const num_comps(rgen.rand() % (is_house ? 3 : 6)); // 0-5
			C.dim   = c.dim; C.dir = !c.dir; // reset dim, flip dir
			C.type  = TYPE_COMPUTER;
			C.shape = SHAPE_CUBE;

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
			C.type  = TYPE_KEYBOARD;
			C.shape = SHAPE_CUBE;

			for (unsigned n = 0; n < num_kbds; ++n) {
				gen_xy_pos_for_cube_obj(C, S, sz, kbd_height, rgen);
				C.dir = rgen.rand_bool();
				add_if_not_intersecting(C, objects, cubes);
			}
		}
		// add bottles
		unsigned const num_bottles(((rgen.rand()&3) == 0) ? 0 : (rgen.rand() % 9)); // 0-8, 75% chance
		C.dir   = C.dim = 0;
		C.type  = TYPE_BOTTLE;
		C.shape = SHAPE_CYLIN;

		for (unsigned n = 0; n < num_bottles; ++n) {
			// same as building_t::place_bottle_on_obj()
			float const bottle_height(z_step*rgen.rand_uniform(0.4, 0.7)), bottle_radius(z_step*rgen.rand_uniform(0.07, 0.11));
			if (min(c_sz.x, c_sz.y) < 5.0*bottle_radius) continue; // shelf not wide/deep enough to add this bottle
			gen_xy_pos_for_round_obj(C, S, bottle_radius, bottle_height, 2.0*bottle_radius, rgen);
			C.set_as_bottle(rgen.rand(), BOTTLE_TYPE_MEDS-1, 1); // all bottle types except for medicine, no_empty=1
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
		}
		// add spraypaint cans
		float const spcan_height(0.55*z_step), spcan_radius(0.17*spcan_height); // fixed size

		if (min(c_sz.x, c_sz.y) > 5.0*spcan_radius) { // add if shelf wide/deep enough
			unsigned const num_spcans(((rgen.rand()%5) < 3) ? (rgen.rand() % 4) : 0); // 0-3, 60% chance
			C.dir   = C.dim = 0;
			C.type  = TYPE_SPRAYCAN;
			C.shape = SHAPE_CYLIN;

			for (unsigned n = 0; n < num_spcans; ++n) {
				gen_xy_pos_for_round_obj(C, S, spcan_radius, spcan_height, 1.5*spcan_radius, rgen);
				C.color  = spcan_colors[rgen.rand() % NUM_SPCAN_COLORS]; // random color
				C.obj_id = rgen.rand(); // used to select emissive color
				add_if_not_intersecting(C, objects, cubes);
			}
		}
		// add large balls to houses
		if (is_house) {
			unsigned const num_balls(rgen.rand() % 3); // 0-2
			C.color = WHITE;
			C.type  = TYPE_LG_BALL;
			C.shape = SHAPE_SPHERE;

			for (unsigned n = 0; n < num_balls; ++n) {
				unsigned const btype(rgen.rand() % NUM_BALL_TYPES);
				float const ball_radius(ball_types[btype].radius*floor_spacing/(12*8)); // assumes 8 foot floor spacing, and bt.radius in inches
				if (2.1*ball_radius > min(c_sz.x, c_sz.y)) continue; // shelf is not wide/deep enough for balls of this type
				gen_xy_pos_for_round_obj(C, S, ball_radius, 2.0*ball_radius, ball_radius, rgen);
				C.item_flags = btype;
				add_if_not_intersecting(C, objects, cubes);
			}
			C.item_flags = 0; // reset for next object type
		}
		// add tape rolls
		float const tape_radius(0.22*z_step), tape_height(TAPE_HEIGHT_TO_RADIUS*tape_radius); // fixed size

		if (min(c_sz.x, c_sz.y) > 3.0*tape_radius) { // add if shelf wide/deep enough
			unsigned const num_tapes(((rgen.rand()%4) < 3) ? (rgen.rand() % 4) : 0); // 0-3, 75% chance
			C.dir   = C.dim = 0;
			C.type  = TYPE_TAPE;
			C.shape = SHAPE_CYLIN;

			for (unsigned n = 0; n < num_tapes; ++n) {
				gen_xy_pos_for_round_obj(C, S, tape_radius, tape_height, 1.25*tape_radius, rgen);
				C.color = tape_colors[rgen.rand() % NUM_TAPE_COLORS]; // random color
				add_if_not_intersecting(C, objects, cubes);
			}
		}
		// add flashlight
		float const fl_height(0.1*floor_spacing), fl_radius(0.2*fl_height);

		if (rgen.rand_float() < 0.25 && min(c_sz.x, c_sz.y) > 3.0*fl_radius) { // add 25% of the time if shelf wide/deep enough
			C.dir   = C.dim = 0;
			C.type  = TYPE_FLASHLIGHT;
			C.shape = SHAPE_CYLIN;
			C.color = BLACK;
			gen_xy_pos_for_round_obj(C, S, fl_radius, fl_height, 1.25*fl_radius, rgen);
			add_if_not_intersecting(C, objects, cubes);
		}
	} // for s
}

void building_room_geom_t::expand_shelves(room_object_t const &c) {
	cube_t shelves[4]; // max number of shelves
	unsigned const num_shelves(get_shelves_for_object(c, shelves));
	get_shelf_objects(c, shelves, num_shelves, expanded_objs);
}

void building_room_geom_t::add_wine_rack_bottles(room_object_t const &c, vect_room_object_t &objects) {
	float const height(c.dz()), width(c.get_width()), depth(c.get_depth()), shelf_thick(0.1*depth);
	unsigned const num_rows(max(1, round_fp(2.0*height/depth))), num_cols(max(1, round_fp(2.0*width/depth)));
	float const row_step((height - shelf_thick)/num_rows), col_step((width - shelf_thick)/num_cols);
	float const space_w(col_step - shelf_thick), space_h(row_step - shelf_thick), diameter(min(space_w, space_h) - shelf_thick);
	rand_gen_t rgen(c.create_rgen());
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
void place_money(room_object_t &obj, cube_t const &parent, float length, float width, unsigned room_id, bool dim, bool dir, rand_gen_t &rgen, unsigned max_bills=20) {
	obj = room_object_t(parent, TYPE_MONEY, room_id, dim, dir);
	obj.z2() = obj.z1() + 0.01f*length*((rgen.rand()%max_bills) + 1); // 1-20 bills
	set_rand_pos_for_sz(obj, dim, length, width, rgen);
}
void place_phone(room_object_t &obj, cube_t const &parent, float length, float width, unsigned room_id, bool dim, bool dir, rand_gen_t &rgen) {
	unsigned const NUM_PHONE_COLORS = 7; // for the case
	colorRGBA const phone_colors[NUM_PHONE_COLORS] = {WHITE, GRAY, DK_GRAY, GRAY_BLACK, BLUE, RED, PINK};
	obj = room_object_t(parent, TYPE_PHONE, room_id, dim, dir);
	obj.color = phone_colors[rgen.rand() % NUM_PHONE_COLORS];
	obj.z2()  = obj.z1() + 0.045*length;
	set_rand_pos_for_sz(obj, dim, length, width, rgen);
}
void set_book_id_and_color(room_object_t &obj, rand_gen_t &rgen) {
	obj.obj_id = rgen.rand();
	obj.color  = book_colors[rgen.rand() % NUM_BOOK_COLORS];
}
void place_book(room_object_t &obj, cube_t const &parent, float length, float max_height, unsigned room_id, bool dim, bool dir, rand_gen_t &rgen) {
	float const width(rgen.rand_uniform(0.6, 1.0)*length);
	obj = room_object_t(parent, TYPE_BOOK, room_id, dim, dir);
	set_book_id_and_color(obj, rgen);
	obj.z2() = obj.z1() + min(0.3f*width, rgen.rand_uniform(0.3, 1.0)*max_height);
	set_rand_pos_for_sz(obj, !dim, length, width, rgen);
}

/*static*/ room_object_t building_room_geom_t::get_item_in_drawer(room_object_t const &c, cube_t const &drawer_in, unsigned drawer_ix, unsigned item_ix, float &stack_z1) {
	room_object_t obj; // starts as no item
	if (stack_z1 == drawer_in.z2()) return obj; // already full
	unsigned const per_drawer_ix(123*c.room_id + 17*c.obj_id + 31*drawer_ix); // per-drawer but not per item

	if (c.type == TYPE_FCABINET || c.type == TYPE_DESK) { // supports up to 4 drawers and 4 items
		unsigned const max_items((c.type == TYPE_DESK) ? 2 : 4);
		if (drawer_ix >= 4 || item_ix >= max_items) return obj; // too many drawers or items
		drawer_ix += 4*item_ix; // encode into a single index
	}
	else { // supports up to 16 drawers
		if (drawer_ix >= 16 || item_ix >= 1) return obj; // too many drawers or items; should never exceed 16 drawers with current room objects
	}
	if (c.item_flags & (1U << drawer_ix)) return obj; // item has been taken
	assert(drawer_in.is_strictly_normalized());
	rand_gen_t rgen;
	rgen.set_state((123*drawer_ix + 777*item_ix + 1), (456*c.room_id + 777*c.obj_id + 1));
	rgen.rand_mix();
	if (item_ix > 0 && rgen.rand_bool()) return obj; // no more items
	cube_t drawer(drawer_in); // copy so that we can adjust z1
	unsigned const type_ix(rgen.rand() % 11); // 0-10
	unsigned const types_dresser [11] = {TYPE_BOX, TYPE_PAPER, TYPE_FOLD_SHIRT, TYPE_FOLD_SHIRT, TYPE_BOOK, TYPE_KEY, TYPE_BOTTLE, TYPE_MONEY,  TYPE_PHONE,  TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_desk    [11] = {TYPE_BOX, TYPE_PAPER, TYPE_PEN, TYPE_STAPLER, TYPE_BOOK, TYPE_KEY,           TYPE_BOTTLE, TYPE_MONEY,  TYPE_PHONE,  TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_attic   [11] = {TYPE_BOX, TYPE_PAPER, TYPE_PEN, TYPE_PEN,     TYPE_BOOK, TYPE_KEY,           TYPE_BOTTLE, TYPE_BOX,    TYPE_BOOK,   TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_kcabinet[11] = {TYPE_BOX, TYPE_PAPER, TYPE_PEN, TYPE_PEN,     TYPE_BOOK, TYPE_PLATE,         TYPE_BOTTLE, TYPE_BOTTLE, TYPE_SILVER, TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_fcabinet[11] = {TYPE_BOX, TYPE_PAPER, TYPE_PEN, TYPE_PEN,     TYPE_BOOK, TYPE_STAPLER,       TYPE_PAPER,  TYPE_BOOK,   TYPE_TAPE,   TYPE_STAPLER,  TYPE_TAPE};
	unsigned obj_type(0);
	if (c.in_attic())                 {obj_type = types_attic   [type_ix];} // custom object overrides for attic item drawers
	else if (c.type == TYPE_DESK)     {obj_type = types_desk    [type_ix];}
	else if (c.type == TYPE_COUNTER)  {obj_type = types_kcabinet[type_ix];}
	else if (c.type == TYPE_FCABINET) {obj_type = types_fcabinet[type_ix];}
	else                              {obj_type = types_dresser [type_ix];} // dresser or nightstand
	if (obj_type == TYPE_SILVER && !building_obj_model_loader.is_model_valid(OBJ_MODEL_SILVER)) {obj_type = TYPE_BOOK;} // replace silverware with book
	// if drawer is too small, replace teeshirt with pen
	if (obj_type == TYPE_FOLD_SHIRT && (drawer.get_sz_dim(c.dim) < 0.55*c.dz() || drawer.get_sz_dim(!c.dim) < 0.52*c.dz())) {obj_type = TYPE_PEN;}
	if (obj_type == TYPE_FOLD_SHIRT && !building_obj_model_loader.is_model_valid(OBJ_MODEL_FOLD_SHIRT)) {obj_type = TYPE_TEESHIRT;} // replace folded shirt with teeshirt
	if (obj_type == TYPE_KEY && item_ix > 0) {obj_type = TYPE_BOTTLE;} // key must be first item/no two kes in one drawer
	// object stacking logic
	bool const is_stackable(obj_type == TYPE_BOX || obj_type == TYPE_PAPER || obj_type == TYPE_BOOK || obj_type == TYPE_PLATE || obj_type == TYPE_TAPE || obj_type == TYPE_FOLD_SHIRT);
	bool const is_single_item(obj_type == TYPE_BOTTLE || obj_type == TYPE_SPRAYCAN); // these two don't combine well with other items since they're large horizontal cylinders
	
	if (item_ix == 0) {stack_z1 = drawer.z1();} // base case
	else if (is_stackable) {
		assert(stack_z1 >= drawer.z1());
		drawer.z1() = stack_z1; // shift bottom of drawer to top of stack
		if (drawer.dz() < 0.1*drawer_in.dz()) return obj; // stack too high
	}
	// else should we mark the area covered by the object as used somehow, and not place another object there?
	vector3d const sz(drawer.get_size()); // Note: drawer is the interior area

	switch (obj_type) {
	case TYPE_BOX: // box
	{
		float length(rgen.rand_uniform(0.4, 0.9)*sz[c.dim]), width(rgen.rand_uniform(0.4, 0.9)*sz[!c.dim]);
		min_eq(length, 1.8f*width); min_eq(width, 1.8f*length); // limit aspect ratio to 1.8
		obj = room_object_t(drawer, TYPE_BOX, c.room_id, rgen.rand_bool(), rgen.rand_bool());
		obj.color = gen_box_color(rgen);
		obj.z2()  = obj.z1() + rgen.rand_uniform(0.4, 0.8)*sz.z;
		set_rand_pos_for_sz(obj, c.dim, length, width, rgen);
		break;
	}
	case TYPE_PAPER: // paper
	{
		float const length(2.0*sz.z), width(0.77*length);

		if (length < 0.9*sz[c.dim] && width < 0.9*sz[!c.dim]) { // if it can fit
			obj = room_object_t(drawer, TYPE_PAPER, c.room_id, c.dim, c.dir);
			obj.obj_id = rgen.rand();
			obj.color  = paper_colors[rgen.rand()%NUM_PAPER_COLORS];
			obj.z2()   = obj.z1() + 0.01*sz.z;
			set_rand_pos_for_sz(obj, c.dim, length, width, rgen);
		}
		break;
	}
	case TYPE_PEN: // pen/pencil/marker
	{
		unsigned const types[3] = {TYPE_PEN, TYPE_PENCIL, TYPE_MARKER}, type(types[rgen.rand()%3]);
		bool const dim(!c.dim); // always opposite orient of the drawer
		float const length(min(1.7f*sz.z, 0.9f*sz[dim])), diameter(((type == TYPE_MARKER) ? 0.08 : 0.036)*length);
		obj = room_object_t(drawer, type, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN);
		obj.color = ((obj.type == TYPE_MARKER) ? marker_colors[rgen.rand()&7] : ((obj.type == TYPE_PEN) ? pen_colors[rgen.rand()&3] : pencil_colors[rgen.rand()&1]));
		obj.z2()  = obj.z1() + diameter;
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case TYPE_BOOK: // book
	{
		float const length(rgen.rand_uniform(0.6, 0.9)*min(sz.x, sz.y));
		place_book(obj, drawer, length, 0.35*sz.z, c.room_id, !c.dim, (c.dir ^ c.dim), rgen);
		break;
	}
	case TYPE_PLATE: // plate
	{
		float const diameter(rgen.rand_uniform(0.5, 0.9)*min(sz.x, sz.y));
		obj = room_object_t(drawer, TYPE_PLATE, c.room_id, 0, 0, 0, 1.0, SHAPE_CYLIN);
		obj.obj_id = rgen.rand();
		obj.z2()   = obj.z1() + rgen.rand_uniform(0.5, 1.0)*min(0.2f*diameter, 0.4f*sz.z); // set height
		set_rand_pos_for_sz(obj, 0, diameter, diameter, rgen);
		break;
	}
	case TYPE_KEY: // key; Note: aspect ratio of key depends on aspect ratio of door, but key model is always a constant aspect ratio
	{
		obj = room_object_t(drawer, TYPE_KEY, c.room_id, rgen.rand_bool(), rgen.rand_bool());
		obj.expand_in_dim( obj.dim, -0.40*sz[ obj.dim]); // long  dim
		obj.expand_in_dim(!obj.dim, -0.46*sz[!obj.dim]); // short dim
		for (unsigned d = 0; d < 2; ++d) {obj.translate_dim(d, 0.35*sz[d]*rgen.rand_uniform(-1.0, 1.0));}
		obj.z2() = obj.z1() + 0.05*sz.z;
		break;
	}
	case TYPE_BOTTLE: // bottle
	{
		bool const dim(c.dim ^ bool(per_drawer_ix & 1)); // random orient, but consistent across the items in the drawer
		float const length(rgen.rand_uniform(0.7, 0.9)*min(((c.type == TYPE_COUNTER) ? 2.7f : 1.8f)*sz.z, min(sz.x, sz.y))), diameter(length*rgen.rand_uniform(0.26, 0.34));
		obj = room_object_t(drawer, TYPE_BOTTLE, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN);
		obj.set_as_bottle(rgen.rand()); // can be empty
		obj.z2() = obj.z1() + diameter;
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case TYPE_MONEY: // money
	{
		float const length(0.135*c.dz()), width(2.35*length);
		if (length < 0.9*sz[c.dim] && width < 0.9*sz[!c.dim]) {place_money(obj, drawer, length, width, c.room_id, c.dim, c.dir, rgen);} // if it can fit
		break;
	}
	case TYPE_SILVER: // silverware
	{
		vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SILVER)); // D, W, H
		float const sw_width(0.5*min(sz[!c.dim], (ssz.y/ssz.x)*sz[c.dim])), sw_length(sw_width*ssz.x/ssz.y), sw_height(sw_width*ssz.z/ssz.y);
		obj = room_object_t(drawer, TYPE_SILVER, c.room_id, c.dim, c.dir);
		obj.z2() = obj.z1() + sw_height; // set height
		set_rand_pos_for_sz(obj, c.dim, sw_length, sw_width, rgen);
		break;
	}
	case TYPE_PHONE: // phone
	{
		bool const dim(c.dim ^ rgen.rand_bool()); // random orient
		float const length(0.3*c.dz()), width(0.45*length);
		if (length < 0.9*sz[dim] && width < 0.9*sz[!dim]) {place_phone(obj, drawer, length, width, c.room_id, dim, c.dir, rgen);} // if it can fit
		break;
	}
	case TYPE_SPRAYCAN: // spray paint can
	{
		bool const dim(c.dim ^ bool(per_drawer_ix & 1)); // random orient, but consistent across the items in the drawer
		float const length(rgen.rand_uniform(0.8, 0.9)*min(1.8f*sz.z, min(sz.x, sz.y))), diameter(0.34*length);
		obj = room_object_t(drawer, TYPE_SPRAYCAN, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN, spcan_colors[rgen.rand() % NUM_SPCAN_COLORS]);
		obj.z2()   = obj.z1() + diameter;
		obj.obj_id = rgen.rand(); // used to select emissive color
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case TYPE_TAPE: // tape roll
	{
		float const diameter(0.3*c.dz());

		if (diameter < 0.9*min(sz.x, sz.y)) { // if it can fit
			obj = room_object_t(drawer, TYPE_TAPE, c.room_id, 0, 0, 0, 1.0, SHAPE_CYLIN, tape_colors[rgen.rand() % NUM_TAPE_COLORS]); // dim/dir don't matter, so use 0
			obj.z2() = obj.z1() + 0.5f*TAPE_HEIGHT_TO_RADIUS*diameter; // set height
			set_rand_pos_for_sz(obj, 0, diameter, diameter, rgen);
		}
		break;
	}
	case TYPE_STAPLER: // stapler
	{
		float const length(min(0.5f*min(sz.x, sz.y), 0.2f*(sz.x + sz.y + sz.z))), width(0.2*length), height(0.3*length);
		obj = room_object_t(drawer, TYPE_STAPLER, c.room_id, rgen.rand_bool(), rgen.rand_bool());
		obj.color = stapler_colors[rgen.rand()%NUM_STAPLER_COLORS];
		obj.z2()  = obj.z1() + height;
		set_rand_pos_for_sz(obj, obj.dim, length, width, rgen);
		break;
	}
	case TYPE_TEESHIRT: // T-shirt
	{
		float length(0.9*drawer.get_sz_dim(c.dim)), width(0.9*drawer.get_sz_dim(!c.dim));
		min_eq(length, 1.10f*width );
		min_eq(width,  1.05f*length);
		cube_t shirt(drawer);
		set_rand_pos_for_sz(shirt, c.dim, length, width, rgen);
		shirt.z2() = shirt.z1() + 0.1*drawer.dz(); // set height
		obj = room_object_t(shirt, TYPE_TEESHIRT, c.room_id, c.dim, c.dir, 0, 1.0, SHAPE_CUBE, TSHIRT_COLORS[rgen.rand()%NUM_TSHIRT_COLORS]);
		break;
	}
	case TYPE_FOLD_SHIRT: // folded shirt
	{
		vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FOLD_SHIRT)); // D, W, H
		float const fs_width(0.75*min(sz[!c.dim], (ssz.y/ssz.x)*sz[c.dim])), fs_length(fs_width*ssz.x/ssz.y), fs_height(fs_width*ssz.z/ssz.y);
		obj = room_object_t(drawer, TYPE_FOLD_SHIRT, c.room_id, c.dim, c.dir, 0, 1.0, SHAPE_CUBE, TSHIRT_COLORS[rgen.rand()%NUM_TSHIRT_COLORS]);
		obj.z2() = obj.z1() + fs_height; // set height
		set_rand_pos_for_sz(obj, c.dim, fs_length, fs_width, rgen);
		break;
	}
	default: assert(0);
	} // end switch
	if (obj.type == TYPE_NONE) return obj; // no other updates
	if (is_stackable) {stack_z1 = obj.z2();} // move the stack up
	else if (is_single_item) {stack_z1 = drawer_in.z2();} // explicitly not stackable
	obj.flags    |= (RO_FLAG_WAS_EXP | RO_FLAG_NOCOLL);
	obj.light_amt = c.light_amt;
	return obj;
}
/*static*/ void building_room_geom_t::add_draw_items(room_object_t const &c, cube_t const &drawer, unsigned drawer_ix, vect_room_object_t &objects) {
	float stack_z1(0.0);

	for (unsigned item_ix = 0; item_ix < 16; ++item_ix) { // will likely break before we hit 16 items, though 16 is the max
		room_object_t const obj(get_item_in_drawer(c, drawer, drawer_ix, item_ix, stack_z1));
		if (obj.type == TYPE_NONE) break;
		objects.push_back(obj);
	}
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
	rand_gen_t rgen(box.create_rgen());
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
				set_book_id_and_color(obj, rgen);
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
				objs.back().set_as_bottle(bottle_id, BOTTLE_TYPE_WINE, 1); // 0-3; excludes poison and medicine; no_empty=1
			}
		}
		else if (obj_type == 2) { // ball - only 1
			if (!is_house) continue; // balls in houses only, not office buildings
			unsigned const btype(rgen.rand() % (NUM_BALL_TYPES-2)); // no beach balls or tennis balls
			float const radius(ball_types[btype].radius*floor_spacing/(12*8)); // assumes 8 foot floor spacing, and bt.radius in inches
			if (c.min_len() < 2.0*radius) continue; // can't fit any of this item
			point const center(cube_bot_center(c) + radius*plus_z);
			cube_t ball; ball.set_from_sphere(center, radius);
			objs.emplace_back(ball, TYPE_LG_BALL, room_id, 0, 0, flags, light_amt, SHAPE_SPHERE);
			objs.back().item_flags = btype;
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
				objs.back().color  = spcan_colors[(color_ix++) % NUM_SPCAN_COLORS];
				objs.back().obj_id = rgen.rand(); // used to select emissive color
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
				objs.emplace_back(*i, TYPE_TAPE, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN); // dim/dir don't matter, so use 0
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
		interior->room_geom->invalidate_small_geom();
		break; // if we got here, something was placed in the box
	} // for n
}

room_obj_or_custom_item_t steal_from_car(room_object_t const &car, float floor_spacing, bool do_pickup) {
	rand_gen_t rgen;
	rgen.set_state(car.obj_id+1, car.get_orient()+1);
	rgen.rand_mix();
	bool const is_locked(rgen.rand_bool()); // or was locked
	room_obj_or_custom_item_t ret; // starts off unassigned

	if (is_locked && !car.is_broken()) { // can't open the door
		if (do_pickup) {print_text_onscreen("Door is locked", RED, 1.0, 2.0*TICKS_PER_SECOND, 0);}
		return ret;
	}
	unsigned const item_counts[16] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5}; // 0-4
	unsigned const num_items(item_counts[rgen.rand()%16]);

	if (car.taken_level >= num_items) { // already looted, nothing left
		if (do_pickup) {print_text_onscreen("You found nothing of value", RED, 0.75, 2.0*TICKS_PER_SECOND, 0);}
		return ret;
	}
	for (unsigned n = 0; n < car.taken_level; ++n) {rgen.rand_mix();} // mix up the state so that each item is different

	if (rgen.rand_bool()) { // select a custom item; most of these are of little value and added for humor
		switch (rgen.rand()%10) {
		case 0: ret.item = custom_item_t("bag of chips",              0.5,  0.1); break;
		case 1: ret.item = custom_item_t("sweaty gym clothes",        2.0,  0.4); break;
		case 2: ret.item = custom_item_t("day old pizza",             0.0,  0.6); break;
		case 3: ret.item = custom_item_t("case of CDs",               5.0,  0.2); break;
		case 4: ret.item = custom_item_t("handgun (with no ammo)",  150.0,  1.0); break;
		case 5: ret.item = custom_item_t("locked briefcase",        200.0, 10.0); break;
		case 6: ret.item = custom_item_t("stuffed animal",            5.0,  0.3); break;
		case 7: ret.item = custom_item_t("pair of shoes",            40.0,  2.0); break;
		case 8: ret.item = custom_item_t("cooler full of beer",      60.0, 20.0); break;
		case 9: ret.item = custom_item_t("box of Girl Scout cookies", 8.0,  1.0); break;
		}
	}
	else { // select a standard room item
		float const v(rgen.rand_float()); // 0.0-1.0

		if (v < 0.15) { // 15% chance for money
			float const length(0.025*floor_spacing); // doesn't really matter, cube is unused/irrelevant
			place_money(ret.obj, car, length, 2.35*length, car.room_id, car.dim, car.dir, rgen);
		}
		else if (v < 0.30) { // 15% chance for cell phone
			float const length(0.06*floor_spacing), width(0.45*length);
			place_phone(ret.obj, car, length, width, car.room_id, car.dim, car.dir, rgen);
		}
		else if (v < 0.40) { // 10% chance for laptop, cube is unused/irrelevant
			ret.obj.type = TYPE_LAPTOP;
		}
		else if (v < 0.65) { // 25% chance for bottle, cube is unused/irrelevant
			ret.obj.type = TYPE_BOTTLE;
			ret.obj.set_as_bottle(rgen.rand(), NUM_BOTTLE_TYPES-1, 0, (1 << BOTTLE_TYPE_POISON)); // all bottle types except for poison, no_empty=0
		}
		else if (v < 0.85) { // 20% chance for book
			float const length(rgen.rand_uniform(0.06, 0.09)*floor_spacing);
			place_book(ret.obj, car, length, 0.2*length, car.room_id, car.dim, car.dir, rgen);
		}
		else { // 15% chance for box, cube is unused/irrelevant
			ret.obj.type   = TYPE_BOX;
			ret.obj.obj_id = rgen.rand();
		}
	}
	return ret;
}

bool building_room_geom_t::expand_object(room_object_t &c, building_t const &building) {
	if (c.obj_expanded()) return 0; // already expanded
	switch (c.type) {
	case TYPE_CLOSET:    expand_closet   (c); break;
	case TYPE_SHELVES:   expand_shelves  (c); break;
	case TYPE_WINE_RACK: expand_wine_rack(c); break;
	case TYPE_CABINET: case TYPE_COUNTER: case TYPE_KSINK: expand_cabinet(c); break;
	case TYPE_MIRROR:    expand_med_cab(c); break;
	case TYPE_BRK_PANEL: expand_breaker_panel(c, !building.interior->elevators.empty(), building.has_parking_garage); break;
	default: assert(0); // not a supported expand type
	}
	if (c.type == TYPE_CLOSET) {maybe_spawn_spider_in_drawer(c, c, 0, building.get_window_vspace(), 1);} // spawn spider when first opened
	c.flags |= RO_FLAG_EXPANDED; // flag as expanded
	return 1;
}

