// 3D World - Building Interior Room Geometry Drawing
// by Frank Gennari 7/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

float const TAPE_HEIGHT_TO_RADIUS = 0.6;

extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;

vect_cube_t &get_temp_cubes();
bool get_dishwasher_for_ksink(room_object_t const &c, cube_t &dishwasher);
void set_wall_width(cube_t &wall, float pos, float half_thick, unsigned dim);
float get_med_cab_wall_thickness(room_object_t const &c);
float get_locker_wall_thickness (room_object_t const &c);
float get_radius_for_square_model(unsigned model_id);
bool place_bottle_on_obj(rand_gen_t &rgen, cube_t const &place_on, vect_room_object_t &objs, float vspace,
	unsigned rid, float lamt, unsigned max_type, vect_cube_t const &avoid, bool at_z1, bool allow_transparent=0);
bool place_dcan_on_obj  (rand_gen_t &rgen, cube_t const &place_on, vect_room_object_t &objs, float vspace,
	unsigned rid, float lamt, unsigned max_type, vect_cube_t const &avoid, bool at_z1);

void resize_model_cube_xy(cube_t &cube, float dim_pos, float not_dim_pos, unsigned id, bool dim) {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(id)); // L, W, H
	set_wall_width(cube, dim_pos,     0.5*cube.dz()*sz.x/sz.z,  dim); // length
	set_wall_width(cube, not_dim_pos, 0.5*cube.dz()*sz.y/sz.z, !dim); // width
}

bool is_shirt_model     (unsigned model_id) {return building_obj_model_loader.model_filename_contains(model_id, "shirt", "Shirt");}
bool is_long_shirt_model(unsigned model_id) {return building_obj_model_loader.model_filename_contains(model_id, "long shirt", "Long Shirt");}
bool is_pants_model     (unsigned model_id) {return building_obj_model_loader.model_filename_contains(model_id, "pants", "Pants");}
bool is_bar_hanger_model(unsigned model_id) {return building_obj_model_loader.model_filename_contains(model_id, "bar hanger", "Bar Hanger");}
bool is_shirt_model     (room_object_t const &obj) {return is_shirt_model     (obj.get_model_id());}
bool is_long_shirt_model(room_object_t const &obj) {return is_long_shirt_model(obj.get_model_id());}
bool is_pants_model     (room_object_t const &obj) {return is_pants_model     (obj.get_model_id());}
bool is_bar_hanger_model(room_object_t const &obj) {return is_bar_hanger_model(obj.get_model_id());}

void add_obj_pair(room_object_t const &obj, vect_room_object_t &objects) { // shoes, etc.
	bool const dim(!obj.dim), mirrored(building_obj_model_loader.get_model(obj.get_model_id()).mirrored);
	room_object_t o1(obj), o2(obj);
	o1.d[dim][1] = o2.d[dim][0] = obj.get_center_dim(dim); // abutting
	((obj.dim ^ obj.dir ^ mirrored) ? o2 : o1).flags |= RO_FLAG_ADJ_TOP; // flag as mirrored; use for shoes; flag either o1 or o2 based on whether the shoe is left/right
	o1.flags |= RO_FLAG_ADJ_LO; // first  in pair
	o2.flags |= RO_FLAG_ADJ_HI; // second in pair
	objects.push_back(o1);
	objects.push_back(o2);
}
bool add_if_not_intersecting(room_object_t const &obj, vect_room_object_t &objects, vect_cube_t &cubes, bool add_obj=1, bool is_pair=0) {
	if (has_bcube_int(obj, cubes)) return 0;
	if (!add_obj) {} // nothing
	else if (is_pair) {add_obj_pair(obj, objects);}
	else {objects.push_back(obj);} // single object
	cubes.push_back(obj);
	return 1;
}
void gen_xy_pos_in_cube(point &pos, cube_t const &c, rand_gen_t &rgen) {
	for (unsigned d = 0; d < 2; ++d) {pos[d] = rgen.rand_uniform(c.d[d][0], c.d[d][1]);}
}
void gen_xyz_pos_in_cube(point &pos, cube_t const &c, rand_gen_t &rgen) {
	for (unsigned d = 0; d < 3; ++d) {pos[d] = rgen.rand_uniform(c.d[d][0], c.d[d][1]);}
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
void gen_xy_pos_for_cube_obj(cube_t &C, cube_t const &S, vector3d const &sz, float height, rand_gen_t &rgen, bool place_at_z1) {
	C.set_from_point(gen_xy_pos_in_area(S, sz, rgen));
	float const place_z(place_at_z1 ? S.z1() : S.z2());
	set_cube_zvals(C, place_z, place_z+height);
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

bool add_obj_to_closet(room_object_t const &c, cube_t const &interior, vect_room_object_t &objects, vect_cube_t &cubes, rand_gen_t &rgen,
	vector3d const &sz, unsigned obj_type, unsigned flags, room_obj_shape shape=SHAPE_CUBE, colorRGBA const &color=WHITE, bool against_back=0)
{
	for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts
		point center(gen_xy_pos_in_area(interior, sz, rgen));
		if (against_back) {center[c.dim] = interior.d[c.dim][!c.dir] + (c.dir ? 1.0 : -1.0)*sz[c.dim];}
		cube_t obj;
		obj.set_from_point(center);
		set_cube_zvals(obj, interior.z1(), (interior.z1() + sz.z));
		obj.expand_by_xy(sz);

		if (!has_bcube_int(obj, cubes)) { // check for intersection with boxes
			objects.emplace_back(obj, obj_type, c.room_id, c.dim, c.dir, flags, c.light_amt, shape, color);
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
	rgen.set_state(obj_id+1, ((type == TYPE_POOL_BALL || type == TYPE_MACHINE) ? item_flags : room_id)+1);
	rgen.rand_mix(); // optional, but improves randomness
	return rgen;
}
ball_type_t const &room_object_t::get_ball_type() const {
	if (type == TYPE_LG_BALL  ) {return ball_types[item_flags % NUM_BALL_TYPES];}
	if (type == TYPE_POOL_BALL) {return pool_ball_type;}
	assert(0);
	return pool_ball_type; // never gets here
}

void building_room_geom_t::add_closet_objects(room_object_t const &c, vect_room_object_t &objects) {
	rand_gen_t rgen(c.create_rgen());
	cube_t ccubes[5]; // only used to get interior space
	get_closet_cubes(c, ccubes);
	cube_t const interior(get_closet_interior_space(c, ccubes));
	float const depth(interior.get_sz_dim(c.dim)), box_sz(0.25*depth), window_vspacing(c.dz()*(1.0 + FLOOR_THICK_VAL_HOUSE));
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);
	unsigned const num_boxes((rgen.rand()%3) + (rgen.rand()%4)); // 0-5
	bool const is_hotel(c.has_extra());
	vect_cube_t &cubes(get_temp_cubes());
	add_boxes_to_space(c, objects, interior, cubes, rgen, num_boxes, box_sz, 0.8*box_sz, 1.5*box_sz, 0, flags); // allow_crates=0

	if (!c.is_small_closet()) { // larger closets have more random items
		vector3d sz;

		if (is_hotel || rgen.rand_bool()) { // maybe add a safe; always add for hotels
			float const sheight(0.15*window_vspacing*rgen.rand_uniform(1.0, 1.2)), swidth(1.0*sheight), sdepth(min(1.1f*sheight, 0.67f*depth));
			sz[c.dim] = 0.5*sdepth; sz[!c.dim] = 0.5*swidth; sz.z = sheight;
			add_obj_to_closet(c, interior, objects, cubes, rgen, sz, TYPE_SAFE, flags, SHAPE_CUBE, colorRGBA(0.7, 0.7, 0.7, 1.0), 1); // against_back=1
		}
		if (rgen.rand_bool()) { // maybe add a lamp in the closet
			try_add_lamp(interior, window_vspacing, c.room_id, flags, c.light_amt, cubes, objects, rgen);
		}
		if (!is_hotel && rgen.rand_bool()) { // maybe add an old computer in the closet
			float const height(0.21*window_vspacing*rgen.rand_uniform(1.0, 1.2)), cheight(0.75*height), cwidth(0.44*cheight), cdepth(0.9*cheight);
			sz[c.dim] = 0.5*cdepth; sz[!c.dim] = 0.5*cwidth; sz.z = cheight;
			add_obj_to_closet(c, interior, objects, cubes, rgen, sz, TYPE_COMPUTER, (flags | RO_FLAG_BROKEN));
		}
		if (!is_hotel && rgen.rand_bool()) { // maybe add a keyboard in the closet
			float const kbd_hwidth(0.12*window_vspacing), kbd_depth(0.6*kbd_hwidth), kbd_height(0.06*kbd_hwidth);
			sz[c.dim] = 0.5*kbd_depth; sz[!c.dim] = kbd_hwidth; sz.z = kbd_height;
			add_obj_to_closet(c, interior, objects, cubes, rgen, sz, TYPE_KEYBOARD, flags);
		}
		if (!is_hotel && rgen.rand_bool()) { // maybe add a paint can in the closet
			float const height(0.64*0.2*window_vspacing), radius(0.28*0.2*window_vspacing);
			sz.assign(radius, radius, height);
			add_obj_to_closet(c, interior, objects, cubes, rgen, sz, TYPE_PAINTCAN, flags, SHAPE_CYLIN);
		}
	}
	// add hanger rod
	float const hr_radius(0.007*window_vspacing);
	room_object_t hanger_rod(interior, TYPE_HANGER_ROD, c.room_id, c.dim, c.dir, (RO_FLAG_NOCOLL | RO_FLAG_INTERIOR)); // SHAPE_CUBE, even though it's a horizontal cylinder
	hanger_rod.z1() = c.z1() + 0.8*window_vspacing;
	hanger_rod.z2() = hanger_rod.z1() + 2.0*hr_radius;
	set_wall_width(hanger_rod, (0.45*c.d[c.dim][c.dir] + 0.55*c.d[c.dim][!c.dir]), hr_radius, c.dim); // move slightly toward the back
	objects.push_back(hanger_rod);
	// add hangers and hanging clothes
	unsigned const num_hangers(c.is_small_closet() ? ((rgen.rand()%7) + 2) : ((rgen.rand()%13) + 4)); // 2-8 / 4-16
	add_hangers_and_clothing(window_vspacing, num_hangers, flags, -1, -1, objects, rgen); // hanger_model_id=clothing_model_id=-1
}

void building_room_geom_t::add_hangers_and_clothing(float window_vspacing, unsigned num_hangers, unsigned flags, int hanger_model_id, int clothing_model_id,
	vect_room_object_t &objects, rand_gen_t &rgen)
{
	if (num_hangers == 0) return;
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_HANGER)) return;
	assert(!objects.empty());
	unsigned const hanger_rod_ix(objects.size() - 1);
	room_object_t const hanger_rod(objects.back()); // hanger rod should have been added previously; deep copy to avoid reference invalidation
	bool const dim(hanger_rod.dim);
	float const hr_radius(0.5*hanger_rod.dz()), wire_radius(0.25*hr_radius);
	if (hanger_rod.get_sz_dim(!dim) < 10.0*wire_radius) return;
	room_object_t hanger(hanger_rod);
	hanger.type       = TYPE_HANGER;
	hanger.item_flags = ((hanger_model_id >= 0) ? hanger_model_id : rgen.rand()); // choose a random hanger sub_model_id per-closet if not forced
	set_cube_zvals(hanger, (hanger_rod.z1() - 0.07*window_vspacing), (hanger_rod.z2() + 1.0*wire_radius));
	unsigned const hid(hanger.get_model_id());
	float const edge_spacing(max(4.0f*hr_radius, 0.06f*window_vspacing));
	float const pos_min(hanger_rod.d[!dim][0] + edge_spacing), pos_max(hanger_rod.d[!dim][1] - edge_spacing);
	float const pos_delta(pos_max - pos_min), slot_spacing(pos_delta/63.0);
	uint64_t slots_used(0); // divide the space into 64 slots, initially all empty
	bool const use_model(building_obj_model_loader.is_model_valid(OBJ_MODEL_CLOTHES));
	bool const mix_hangers(hanger_model_id < 0 && rgen.rand_bool());

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
		resize_model_cube_xy(hanger, hanger.get_center_dim(dim), (pos_min + slot_ix*slot_spacing), hid, dim);
		expand_to_nonzero_area(hanger, 0.1*hr_radius, !dim);
		objects.push_back(hanger);
			
		if (use_model && rgen.rand_float() < 0.8) { // maybe add clothing to the hanger
			objects.back().flags |= RO_FLAG_HANGING; // flag the hanger has having clothing hanging on it
			room_object_t clothes(hanger, TYPE_CLOTHES, hanger_rod.room_id, dim, hanger_rod.dir, (flags | RO_FLAG_HANGING), hanger_rod.light_amt, SHAPE_CUBE, WHITE);
			clothes.z2() -= 0.55*hanger.dz(); // top
			clothes.z1() -= 0.3*window_vspacing/(1.0 + FLOOR_THICK_VAL_HOUSE); // bottom
			clothes.item_flags = ((clothing_model_id >= 0) ? clothing_model_id : rgen.rand()); // choose a random clothing sub_model_id if not forced
			if (is_pants_model(clothes) && is_bar_hanger_model(hanger)) {++clothes.item_flags;} // hack to avoid placing pants on bar hangers
			unsigned const cid(clothes.get_model_id());
			float const scale(building_obj_model_loader.get_model_scale(cid));
			if (scale != 0.0) {set_wall_width(clothes, clothes.zc(), 0.5*scale*clothes.dz(), 2);} // rescale zvals around the center
			resize_model_cube_xy(clothes, clothes.get_center_dim(dim), clothes.get_center_dim(!dim), cid, dim);
			bool skip(0);

			for (auto m = objects.begin()+hanger_rod_ix+1; m != objects.end(); ++m) { // skip if this object intersects a previous hanging clothing item
				if (m->type == TYPE_CLOTHES && m->intersects(clothes)) {skip = 1; break;}
			}
			if (skip) continue;

			if (is_shirt_model(clothes)) {
				if (rgen.rand_float() < 0.6) { // 60% of shirts are randomly colored rather than colored + textured with the model
					clothes.color  = shirt_colors[rgen.rand()%NUM_SHIRT_COLORS];
					clothes.flags |= RO_FLAG_UNTEXTURED;
				}
				else if (rgen.rand_float() < 0.7) { // 70% of textured shirts have the store name/logo
					clothes.flags |= RO_FLAG_HAS_EXTRA;
				}
			}
			objects.push_back(clothes);
		}
		if (mix_hangers && (rgen.rand()&7) == 0) {hanger.item_flags = rgen.rand();} // switch to a new hanger type occasionally
	} // for i
	objects[hanger_rod_ix].item_flags = uint16_t(objects.size() - hanger_rod_ix); // number of objects hanging on the hanger rod, including hangers and shirts
}

bool add_cabinet_objects(room_object_t const &c, vect_room_object_t &objects) { // called on cabinets, counters, kitchen sinks, and bathroom vanities
	rand_gen_t rgen(c.create_rgen());
	vect_cube_t &cubes(get_temp_cubes());
	bool const is_vanity(c.type == TYPE_VANITY), in_kitchen(!is_vanity), is_counter(c.type == TYPE_COUNTER || c.type == TYPE_KSINK || is_vanity);
	bool const in_hospital(is_vanity && !c.is_house());
	float const wall_thickness(0.04*(is_counter ? 0.95 : 1.0)*c.dz()); // cabinet height is 95% of total counter; the other 5% is the top surface
	cube_t interior(c), dishwasher;
	if (is_counter) {interior.d[c.dim][c.dir] -= (c.dir ? 1.0 : -1.0)*0.05*c.get_depth();} // subtract front overhang
	interior.expand_by(-wall_thickness);
	vector3d const c_sz(interior.get_size());
	float const light_amt(0.25*c.light_amt), c_min_xy(min(c_sz.x, c_sz.y));
	unsigned const sz_ratio(round_fp(c_sz[!c.dim]/c_sz.z));

	if (c.type == TYPE_KSINK && get_dishwasher_for_ksink(c, dishwasher)) { // avoid placing objects that overlap the dishwasher
		dishwasher.d[c.dim][!c.dir] = c.d[c.dim][!c.dir]; // extend to the back of the cabinet
		dishwasher.expand_by_xy(wall_thickness);
		dishwasher.z1() = c.z1(); // no gap at the bottom
		cubes.push_back(dishwasher);
	}
	unsigned const start_num_cubes(cubes.size()), flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);

	if (is_counter && rgen.rand_float() < (is_vanity ? 0.25 : 0.5)) { // maybe add a trashcan
		float const tcan_height(c_sz.z*rgen.rand_uniform(0.35, 0.55)), tcan_radius(min(tcan_height/rgen.rand_uniform(1.6, 2.8), 0.4f*c_min_xy));
		cube_t tcan;
		gen_xy_pos_for_round_obj(tcan, interior, tcan_radius, tcan_height, 1.1*tcan_radius, rgen, 1); // place_at_z1=1
		tcan.translate_dim(2, 0.01*tcan_height); // move up slightly to prevent Z-fighting
		room_object_t obj(tcan, TYPE_TCAN, c.room_id, c.dim, c.dir, flags, light_amt, (rgen.rand_bool() ? SHAPE_CYLIN : SHAPE_CUBE), tcan_colors[rgen.rand()%NUM_TCAN_COLORS]);
		add_if_not_intersecting(obj, objects, cubes);
	}
	// maybe add fire extinguisher
	if (is_counter && building_obj_model_loader.is_model_valid(OBJ_MODEL_FIRE_EXT) && rgen.rand_float() < 0.3) { // 30% of the time
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FIRE_EXT)); // D, W, H
		float const fe_height(0.38*c.dz()), fe_radius(fe_height*(0.5*(sz.x + sz.y)/sz.z));

		if (2.2*fe_radius < c_min_xy) {
			cube_t fire_ext;
			gen_xy_pos_for_round_obj(fire_ext, interior, fe_radius, fe_height, 1.05*fe_radius, rgen, 1); // place_at_z1=1
			room_object_t obj(fire_ext, TYPE_FIRE_EXT, c.room_id, rgen.rand_bool(), rgen.rand_bool(), flags, light_amt, SHAPE_CYLIN);
			add_if_not_intersecting(obj, objects, cubes);
		}
	}
	// add boxes
	unsigned const num_boxes(rgen.rand()%4); // 0-3
	float const box_sz(0.3*c.get_length()), sz_scale(0.7*c_sz.z);
	room_object_t cb(c);
	cb.light_amt = light_amt;
	add_boxes_to_space(cb, objects, interior, cubes, rgen, num_boxes, box_sz, 0.8*box_sz, 1.5*box_sz, 0, (flags | RO_FLAG_NOCOLL)); // allow_crates=0

	if (!in_hospital && !is_vanity) { // add paint cans (slightly smaller than normal)
		float const pc_height(0.6*sz_scale), pc_radius(0.24*sz_scale);

		if (3*pc_radius < c_min_xy) { // have enough space for for paint cans
			unsigned const num_pcans(rgen.rand()%3); // 0-2

			for (unsigned n = 0; n < num_pcans; ++n) {
				cube_t pcan;
				gen_xy_pos_for_round_obj(pcan, interior, pc_radius, pc_height, 1.2*pc_radius, rgen, 1); // place_at_z1=1
				room_object_t obj(pcan, TYPE_PAINTCAN, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
				add_if_not_intersecting(obj, objects, cubes);
			}
		}
	}
	if (in_kitchen) {
		// add plates
		unsigned const max_plates(0 + 1*sz_ratio), num_plates(rgen.rand() % max_plates); // wider cabinet has more plates
		float const plate_radius(min(sz_scale*rgen.rand_uniform(0.30, 0.35), 0.45f*c_min_xy)), plate_height(0.1*plate_radius);

		for (unsigned n = 0; n < num_plates; ++n) {
			cube_t plate;
			gen_xy_pos_for_round_obj(plate, interior, plate_radius, plate_height, 1.2*plate_radius, rgen, 1); // place_at_z1=1
			room_object_t obj(plate, TYPE_PLATE, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
			if (!add_if_not_intersecting(obj, objects, cubes)) continue; // can't place the bottom plate
			unsigned const stack_height(1 + (rgen.rand()%5)); // 1-6

			for (unsigned s = 1; s < stack_height; ++s) {
				obj.translate_dim(2, plate_height); // shift up in z
				if (obj.z2() + plate_height > interior.z2()) break; // stack is too high, end it here
				objects.push_back(obj);
			}
		} // for n
		// add pans
		unsigned const num_pans(rgen.rand()%3); // 0-2
		float const pan_radius(min(sz_scale*rgen.rand_uniform(0.40, 0.45), 0.45f*c_min_xy)), pan_height(rgen.rand_uniform(0.4, 0.5)*pan_radius);

		for (unsigned n = 0; n < num_pans; ++n) {
			cube_t pan;
			gen_xy_pos_for_round_obj(pan, interior, pan_radius, pan_height, 1.1*pan_radius, rgen, 1); // place_at_z1=1
			pan.translate_dim(2, 0.01*pan_height); // shift in +z to prevent z-fighting with the bottom of the pan
			bool const dir(pan.get_center_dim(!c.dim) < c.get_center_dim(!c.dim)); // point toward the cabinet center to avoid the handle clipping through the side/end
			room_object_t obj(pan, TYPE_PAN, c.room_id, c.dim, dir, flags, light_amt, SHAPE_CYLIN, GRAY_BLACK);
			cube_t const pan_bc(get_pan_bcube_inc_handle(obj));
			if (has_bcube_int(pan_bc, cubes)) continue;
			objects.push_back(obj);
			cubes.push_back(pan_bc);
		} // for n
	}
	if (in_hospital || !is_vanity) { // add bottles; wider cabinet has more; fewer medicine bottles in hospitals
		unsigned const max_bottles(3 + (in_hospital ? 1 : 2)*sz_ratio), num_bottles(rgen.rand() % max_bottles);

		for (unsigned n = 0; n < num_bottles; ++n) {
			float const bottle_height(sz_scale*rgen.rand_uniform(0.45, 0.65)), bottle_radius(sz_scale*rgen.rand_uniform(0.075, 0.1));
			if (c_min_xy < 3.0*bottle_radius) continue; // cabinet not wide/deep enough to add this bottle
			cube_t bottle;
			gen_xy_pos_for_round_obj(bottle, interior, bottle_radius, bottle_height, 1.5*bottle_radius, rgen, 1); // place_at_z1=1
			room_object_t obj(bottle, TYPE_BOTTLE, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN); // vertical

			if (in_hospital) {obj.set_as_bottle(rgen.rand(), BOTTLE_TYPE_MEDS, 0, ~(1 << BOTTLE_TYPE_MEDS));} // medicine only
			else {
				bool const allow_medicine(rgen.rand_bool()); // medicine is half as common
				obj.set_as_bottle(rgen.rand(), (allow_medicine ? (unsigned)NUM_BOTTLE_TYPES : (unsigned)BOTTLE_TYPE_MEDS)-1, 1); // all bottle types, no_empty=1
			}
			add_if_not_intersecting(obj, objects, cubes);
		}
	}
	if (in_hospital) {
		// add tape rolls
		unsigned const num_tape_rolls(rgen.rand()%6); // 0-5
		colorRGBA const h_tape_colors[3] = {WHITE, WHITE, BLUE};

		for (unsigned n = 0; n < num_tape_rolls; ++n) {
			float const tape_height(sz_scale*rgen.rand_uniform(0.08, 0.15)), tape_radius(min(0.4f*c_min_xy, sz_scale*rgen.rand_uniform(0.14, 0.22)));
			cube_t tape;
			gen_xy_pos_for_round_obj(tape, interior, tape_radius, tape_height, 1.2*tape_radius, rgen, 1); // place_at_z1=1
			room_object_t const obj(tape, TYPE_TAPE, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN, h_tape_colors[rgen.rand()%3]); // vertical
			add_if_not_intersecting(obj, objects, cubes);
		}
	}
	else { // non-hospital
		if (is_vanity) { // add toilet paper rolls and paper towels
			unsigned const num_rolls(rgen.rand()%3); // 0-2

			for (unsigned n = 0; n < num_rolls; ++n) {
				bool const paper_towel(rgen.rand_bool());
				float const height((paper_towel ? 0.64 : 0.3)*sz_scale), radius((paper_towel ? 0.25 : 0.5)*height);
				if (3.0*radius > c_min_xy) continue; // cabinet not wide/deep enough to add a roll
				cube_t roll;
				gen_xy_pos_for_round_obj(roll, interior, radius, height, 1.1*radius, rgen, 1); // place_at_z1=1
				room_object_t const obj(roll, TYPE_TPROLL, c.room_id, 0, 0, (flags | (paper_towel ? RO_FLAG_HAS_EXTRA : 0)), light_amt, SHAPE_CYLIN); // vertical
				add_if_not_intersecting(obj, objects, cubes);
			} // for n
		}
		else { // add drink cans
			unsigned const max_cans(2 + 2*sz_ratio), num_cans(rgen.rand() % max_cans); // wider cabinet has more cans
			float const can_height(0.32*sz_scale), can_radius(0.26*can_height);

			if (c_min_xy > 3.0*can_radius) { // cabinet wide/deep enough to add a can
				for (unsigned n = 0; n < num_cans; ++n) {
					cube_t can;
					gen_xy_pos_for_round_obj(can, interior, can_radius, can_height, 1.25*can_radius, rgen, 1); // place_at_z1=1
					room_object_t obj(can, TYPE_DRINK_CAN, c.room_id, 0, 0, flags, light_amt, SHAPE_CYLIN); // vertical
					obj.obj_id = rgen.rand(); // random can type
					add_if_not_intersecting(obj, objects, cubes);
				}
			}
		}
	}
	return (cubes.size() > start_num_cubes); // returns true if some object was added
}
void building_room_geom_t::expand_cabinet(room_object_t const &c) { // called on cabinets, counters, kitchen sinks, and bathroom vanities
	bool any_objs_added(add_cabinet_objects(c, expanded_objs));
	if (any_objs_added) {invalidate_small_geom();} // some object was added
}

void building_room_geom_t::expand_med_cab(room_object_t const &c) {
	rand_gen_t rgen(c.create_rgen());
	// add medicine bottle
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);
	float const wall_thickness(get_med_cab_wall_thickness(c));
	cube_t interior(c), bottle;
	interior.expand_by(-wall_thickness);
	if (rgen.rand_bool()) {interior.z1() += 0.5*(interior.dz() + wall_thickness);} // place on middle shelf half the time
	float const height(0.33*c.dz()*rgen.rand_uniform(0.75, 1.0)), radius(0.4*interior.get_sz_dim(c.dim)*rgen.rand_uniform(0.75, 1.0));
	gen_xy_pos_for_round_obj(bottle, interior, radius, height, 1.1*radius, rgen, 1); // place_at_z1=1
	room_object_t obj(bottle, TYPE_BOTTLE, c.room_id, 0, 0, flags, c.light_amt, SHAPE_CYLIN); // vertical
	obj.set_as_bottle(BOTTLE_TYPE_MEDS, BOTTLE_TYPE_MEDS, 1); // medicine, no_empty=1
	expanded_objs.push_back(obj);
	invalidate_small_geom();
}

void building_room_geom_t::expand_breaker_panel(room_object_t const &c, building_t const &building) {
	float const box_width(c.get_sz_dim(!c.dim)), box_depth(c.get_sz_dim(c.dim)), box_height(c.dz());
	float const thickness(0.25*box_depth), dir_sign(c.dir ? -1.0 : 1.0);
	cube_t breakers(c);
	breakers.d[c.dim][0] = breakers.d[c.dim][1] = breakers.d[c.dim][!c.dir] - dir_sign*thickness; // move both edges to front of panel
	breakers.d[c.dim][!c.dir] += dir_sign*0.5*thickness; // expand out for correct width
	breakers.expand_in_dim(!c.dim, -0.2*box_width); // shrink the width
	breakers.expand_in_dim(2,      -0.1*box_height); // shrink vertically
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
			breaker_zone_t const zone(building.interior->get_circuit_breaker_info(obj.obj_id, num_breakers, building.get_window_vspace()));
			if (zone.invalid() || zone.is_dup) continue; // no zone/label
			float label_xlate(0.0);
			if      (C == 0         ) {label_xlate = -0.75;} // first column/lo side
			else if (C == num_cols-1) {label_xlate =  0.75;} // last column/hi side
			else {continue;} // else middle column; no space to add a label
			cube_t label(breaker);
			label.translate_dim(!c.dim, label_xlate*breaker_dc);
			string label_text(room_names_short[zone.rtype]); // ignores store names
			if (label_text.empty()) continue;

			// special breaker labeling for apartments and hotels;
			// normal residential rooms (to re-label) are before RTYPE_STORAGE and special rooms such as RTYPE_UTILITY (keep the label) are after
			if (building.is_apt_or_hotel() && zone.rtype < RTYPE_STORAGE && zone.pri_room >= 0) {
				room_t const &pri_room(building.get_room(zone.pri_room));

				if (pri_room.is_apt_or_hotel_room()) { // Note: pri_room.unit_id does *not* match the room number on the hallway sign
					label_text = string(building.is_hotel() ? "Rooms" : "Apts") + " " + std::to_string(obj.obj_id);
				}
			}
			label.d[c.dim][!c.dir] = label.d[c.dim][c.dir] + dir_sign*0.05*thickness;
			label.expand_in_dim(!c.dim, -0.2*breaker_dc); // small shrink
			label.expand_in_dim(2,      -0.2*breaker_dr); // small shrink
			expanded_objs.emplace_back(label, TYPE_SIGN, c.room_id, c.dim, !c.dir, RO_FLAG_NOCOLL, c.light_amt, SHAPE_CUBE, BLACK);
			expanded_objs.back().obj_id = register_sign_text(label_text);
		} // for r
	} // for C
	invalidate_small_geom();
}

void building_room_geom_t::expand_dishwasher(room_object_t &c, cube_t const &dishwasher) {
	unsigned const expanded_objs_start(expanded_objs.size());
	c.item_flags = expanded_objs_start; // record the location where we'll add our objects
	//rand_gen_t rgen(c.create_rgen());
	// add TYPE_CUP, and TYPE_SILVER? these are models and are more difficult to expand
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

unsigned room_object_t::get_num_shelves() const {
	assert(type == TYPE_SHELVES);
	if (in_mall() && item_flags == STORE_PETS) return 2; // pet store shelves always have two levels
	if (on_warehouse_floor())                  return 4; // warehouse special case, always 4 shelves if on the floor
	if (in_warehouse())                        return 3; // warehouse special case, always 3 shelves if on a room roof
	return (2 + (room_id % (in_mall() ? 2 : 3))); // 2-4 shelves, 2-3 for mall clothing stores
}
unsigned get_shelves_for_object(room_object_t const &c, cube_t shelves[MAX_SHELVES]) {
	unsigned const num_shelves(c.get_num_shelves()); // 2-4 shelves
	assert(num_shelves <= MAX_SHELVES);
	bool const in_warehouse(c.in_warehouse());
	float const bot_space_ratio(in_warehouse ? 0.25 : ((num_shelves == 2) ? 0.75 : 1.0)); // small space for warehouse shelves; less space for the 2 shelves case
	float const thickness((in_warehouse ? 0.01 : 0.02)*c.dz()), bracket_thickness(0.8*thickness);
	float const z_step(c.dz()/(num_shelves + bot_space_ratio)); // include a space at the bottom
	cube_t shelf(c);
	shelf.z2() = shelf.z1() + thickness; // set shelf thickness
	shelf.d[c.dim][c.dir] += (c.dir ? -1.0 : 1.0)*bracket_thickness; // leave space behind the shelf for brackets

	for (unsigned s = 0; s < num_shelves; ++s) {
		shelf.translate_dim(2, ((s == 0) ? bot_space_ratio : 1.0)*z_step); // move up one step; done first, so bottom shelf is above the floor
		shelves[s] = shelf; // record for later use
	}
	return num_shelves;
}

void set_spraypaint_color(room_object_t &obj, rand_gen_t &rgen, unsigned *color_ix=nullptr) {
	obj.color  = spcan_colors[(color_ix ? (*color_ix)++ : rgen.rand()) % NUM_SPCAN_COLORS];
	obj.obj_id = rgen.rand(); // used to select emissive color
}

void building_room_geom_t::get_shelf_objects(room_object_t const &c_in, cube_t const shelves[MAX_SHELVES], unsigned num_shelves, vect_room_object_t &objects, bool add_models_mode) {
	if (!c_in.is_nonempty()) return; // empty - no objects
	room_object_t c(c_in); // deep copy so that we can set flags and don't invalidate the reference
	if (!add_models_mode) {c.flags |= RO_FLAG_WAS_EXP;} // only expand for non-models mode
	bool const dim(c.dim), dir(c.dir), is_house(c.is_house()), in_mall(c.in_mall());
	bool const in_warehouse(c.in_warehouse()), on_warehouse_floor(c.on_warehouse_floor());
	vector3d const c_sz(c.get_size());
	float const dz(c_sz.z), width(c_sz[dim]), thickness(0.02*dz), bracket_thickness(0.75*thickness), floor_spacing(1.1*dz);
	float const z_step(dz/(num_shelves + 1)), shelf_clearance(z_step - thickness - bracket_thickness), sz_scale(is_house ? 0.7 : (in_warehouse ? 0.8 : 1.0));
	float const box_zscale(shelf_clearance*sz_scale), h_val(0.21*floor_spacing), length_ratio(c_sz[!dim]/width);
	rand_gen_t base_rgen(c.create_rgen());
	bool prev_add_shoe_boxes(0);

	// Note: this function supports placement of objects drawn as 3D models such as fire extinguishers when add_models_mode=1
	for (unsigned s = 0; s < num_shelves; ++s) {
		cube_t const &S(shelves[s]);
		vector3d const s_sz(S.get_size());
		float const s_sz_min_xy(min(s_sz.x, s_sz.y));
		bool const top_shelf(s+1 == num_shelves);
		vect_cube_t &cubes(get_temp_cubes());
		base_rgen.rand_mix();
		rand_gen_t rgen(base_rgen); // create a new rgen so that we can early exit this loop based on add_models_mode and still get deterministic results
		room_object_t C(c);
		vector3d sz;

		if (in_mall) {
			float const shelf_depth(s_sz[dim]), shelf_len(s_sz[!dim]), ld_ratio(shelf_len/shelf_depth);

			if (c.item_flags == STORE_CLOTHING) { // clothing store shelf
				unsigned const num_items(round_fp(rgen.rand_uniform(0.5, 1.0)*ld_ratio)); // 50-100% full
				C.shape = SHAPE_CUBE;
				C.dim   =  dim;
				C.dir   = !dir;

				if (building_obj_model_loader.is_model_valid(OBJ_MODEL_FOLD_SHIRT) && rgen.rand_bool()) { // 50% chance of shirt models
					vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FOLD_SHIRT)); // D, W, H
					float const length(rgen.rand_uniform(0.55, 0.65)*shelf_depth);
					float const width(min(0.5f*shelf_len, length*ssz.x/ssz.y)), height(length*ssz.z/ssz.y);
					C.type = TYPE_FOLD_SHIRT;
					C.z2() = C.z1() + height; // set height
					sz[ dim] = 0.5*length;
					sz[!dim] = 0.5*width;

					for (unsigned n = 0; n < num_items; ++n) {
						// should folded shirts be stacked? this will increase drawing/culling time
						C.color = gen_teeshirt_color(rgen);
						gen_xy_pos_for_cube_obj(C, S, sz, height, rgen);
						add_if_not_intersecting(C, objects, cubes, add_models_mode);
					}
					continue; // done with this shelf
				}
				if (add_models_mode) continue;
				bool const is_teeshirt(rgen.rand_bool());
				float const length((is_teeshirt ? 1.0 : 0.8)*rgen.rand_uniform(0.9, 1.0)*shelf_depth), height(0.02*shelf_depth);
				C.type = (is_teeshirt ? TYPE_TEESHIRT : TYPE_PANTS);
				C.z2() = C.z1() + height; // set height lower
				sz[ dim] = 0.5*length;
				sz[!dim] = 0.5*min(0.98f*length, 0.5f*shelf_len); // width

				for (unsigned n = 0; n < num_items; ++n) {
					C.color = gen_shirt_pants_color(C.type, rgen);
					gen_xy_pos_for_cube_obj(C, S, sz, height, rgen);
					add_if_not_intersecting(C, objects, cubes);
				}
			}
			else if (c.item_flags == STORE_PETS) { // pet store shelf
				if (add_models_mode) { // fishtanks and cages count as models since they have fish models and are added to objs rather than expanded_objs
					unsigned const num_place(round_fp(0.4*rgen.rand_uniform(0.5, 1.0)*ld_ratio));
					vector<unsigned> animal_types{TYPE_FISH, TYPE_SNAKE, TYPE_SPIDER};
					if (building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT      )) {animal_types.push_back(TYPE_RAT );}
					if (building_obj_model_loader.is_model_valid(OBJ_MODEL_BIRD_ANIM)) {animal_types.push_back(TYPE_BIRD);}
					unsigned const animal_type(animal_types[rgen.rand() % animal_types.size()]);
					cube_t tank;

					for (unsigned n = 0; n < num_place; ++n) {
						float const height((top_shelf ? 1.2 : 1.0)*rgen.rand_uniform(0.82, 0.92)*shelf_clearance);
						sz[ dim] = 0.5*rgen.rand_uniform(0.86, 0.98)*shelf_depth; // set depth
						sz[!dim] = 0.5*min(0.5f*shelf_len, rgen.rand_uniform(1.6, 2.8)*shelf_depth); // set length
						gen_xy_pos_for_cube_obj(tank, S, sz, height, rgen);
						if (has_bcube_int(tank, cubes)) continue;
						add_pet_container(tank, C.room_id, C.light_amt, dim, !dir, 1, objects, rgen, animal_type, s); // in_pet_store=1
						tank.expand_in_dim(!dim, 0.05*shelf_depth); // add a bit of extra spacing between tanks
						cubes.push_back(tank);
					} // for n
				}
				// else added during animal update pass
			}
			else if (c.item_flags == STORE_SHOE) { // shoe store shelf; // add shoes displayed sideways
				if (add_models_mode) {
					float length(0.0), spacing(0.0);
					unsigned model_id(0);
					bool add_pairs(0);

					// add shoe boxes rather than shoes, with a preference toward bottom shelves
					if ((s == 0 || prev_add_shoe_boxes) && building_obj_model_loader.is_model_valid(OBJ_MODEL_SHOEBOX) && rgen.rand_float() > (0.5*s + 0.4)) {
						model_id = OBJ_MODEL_SHOEBOX;
						C.type   = TYPE_SHOEBOX;
						C.dir    = dir; // same dir as shelf so that box opens outward
						C.color  = WHITE;
						length   = 0.3*floor_spacing;
						spacing  = 1.1*length;
						prev_add_shoe_boxes = 1;
					}
					else { // add shoes
						assert(building_obj_model_loader.is_model_valid(OBJ_MODEL_SHOE));
						model_id = OBJ_MODEL_SHOE;
						C.type   = TYPE_SHOE;
						C.dir    = rgen.rand_bool(); // random orient, but consistent per shelf
						C.color  = WHITE; // currently always white; most shoe models are a single textured material, and it doesn't look right to color every part
						length   = 0.25*floor_spacing;
						spacing  = 1.25*length;
						add_pairs= (shelf_len < 8.0*floor_spacing); // only for short shelves, not full wall shelves
					}
					vector3d const ssz(building_obj_model_loader.get_model_world_space_size(model_id)); // L, W, H
					float const width(min(0.5*shelf_len, (add_pairs ? 2.0 : 1.0)*length*ssz.y/ssz.x)), height(length*ssz.z/ssz.x); // set max
					C.shape = SHAPE_CUBE;
					C.dim   = !dim;
					C.z1()  = S.z2();
					unsigned const num_slots(shelf_len / spacing);
					float const slot_spacing(shelf_len / num_slots);

					for (unsigned n = 0; n < num_slots; ++n) {
						if (rgen.rand_float() < 0.25) continue; // no shoe in this slot
						C.item_flags = rgen.rand(); // random shoe sub-model
						float const scale(rgen.rand_uniform(0.85, 1.0));
						C.z2() = C.z1() + scale*height; // set height
						set_wall_width(C, c.get_center_dim(dim), 0.5*scale*width, dim);
						set_wall_width(C, (c.d[!dim][0] + (n + 0.5 + 0.05*rgen.signed_rand_float())*slot_spacing), 0.5*scale*length, !dim);
						add_if_not_intersecting(C, objects, cubes, add_models_mode, add_pairs);
					} // for n
				}
			}
			else {assert(0);} // unsupported store type
			continue;
		}
		if (!in_warehouse) { // add fire extinguishers if we have a model and it fits (2 level shelves)
			float const max_height((top_shelf ? 1.0 : 0.8)*shelf_clearance); // more space on the top shelf
			float fe_height(0.0), fe_radius(0.0);

			if (get_fire_ext_height_and_radius(floor_spacing, fe_height, fe_radius) && fe_height < max_height && 2.2*fe_radius < s_sz_min_xy) {
				if (rgen.rand_float() < 0.3) { // 30% of the time
					C.color = WHITE;
					C.type  = TYPE_FIRE_EXT;
					C.shape = SHAPE_CYLIN;
					C.dim   = rgen.rand_bool(); C.dir = rgen.rand_bool(); // random orient
					gen_xy_pos_for_round_obj(C, S, fe_radius, fe_height, 1.1*fe_radius, rgen);
					add_if_not_intersecting(C, objects, cubes, add_models_mode);
				}
			}
		}
		if (add_models_mode) continue; // all items added below are non-models
		unsigned const objs_start(objects.size());

		// add crates/boxes
		unsigned num_boxes(0);

		if (in_warehouse) {
			unsigned const max_boxes(round_fp(1.8*length_ratio));
			num_boxes = rgen.rand_uniform_uint(max_boxes/2, max_boxes);
		}
		else {
			num_boxes = (rgen.rand() % (is_house ? 8 : 13)); // 0-7/12
		}
		cube_t bounds(S);
		bounds.z1() = S.z2(); // place on top of shelf
		// allow_crates=1; copy flags from shelf, except for open flag
		unsigned const box_flags((c.flags & ~RO_FLAG_OPEN) | RO_FLAG_NOCOLL);
		add_boxes_to_space(c, objects, bounds, cubes, rgen, num_boxes, 0.42*width*sz_scale*(is_house ? 1.5 : 1.0), 0.4*box_zscale, 0.98*box_zscale, 1, box_flags);

		if (on_warehouse_floor) { // add pallets; only for shelves on main warehouse floor
			float const pdepth(0.96*s_sz[dim]), pwidth(40.0/48*pdepth), pheight(4.5/48*pdepth); // 48x40x4.5 inches, scaled to shelf depth

			if (pwidth < 3.0*s_sz[!dim]) { // don't place on short shelves
				float const half_width(0.5*pwidth);
				unsigned const max_pallets(round_fp(0.6*length_ratio)), num_pallets(rgen.rand_uniform_uint(max_pallets/2, max_pallets));
				C.dim    = dim; C.dir = 0;
				C.type   = TYPE_PALLET;
				C.shape  = SHAPE_CUBE;
				C.obj_id = objs_start; // store start of shelf objects so that we can more quickly find an object stacked on this pallet (future work)
				set_cube_zvals(C, S.z2(), (S.z2() + pheight));
				set_wall_width(C, S.get_center_dim(dim), 0.5*pdepth, dim); // set depth
				vector<vect_room_object_t::iterator> to_move_up;

				for (unsigned n = 0; n < num_pallets; ++n) {
					for (unsigned N = 0; N < 10; ++N) { // N placement attempts
						set_wall_width(C, rgen.rand_uniform((S.d[!dim][0] + half_width), (S.d[!dim][1] - half_width)), half_width, !dim);
						bool valid(1);
						to_move_up.clear();

						// Note: not using cubes because logic depends on object type
						for (auto i = objects.begin()+objs_start; i != objects.end(); ++i) {
							if (!i->intersects(C)) continue;
							if (i->type != TYPE_BOX && i->type != TYPE_CRATE) {valid = 0; break;} // another pallet, etc.
							cube_t contain_area(C);
							contain_area.expand_by_xy(0.2*i->get_size_xy()); // all a bit of overlap
							if (!contain_area.contains_cube_xy(*i))  {valid = 0; break;} // partial overlap
							if (i->dz() + pheight > shelf_clearance) {valid = 0; break;} // too tall
							to_move_up.push_back(i);
						} // for i
						if (!valid) continue;
						for (auto const &i : to_move_up) {i->translate_dim(2, pheight);} // move up onto pallet
						objects.push_back(C);
						if (!to_move_up.empty()) {objects.back().flags |= RO_FLAG_ADJ_TOP;} // flag as having something on it
						break; // done/success
					} // for N
				} // for n
				C.obj_id = 0; // reset to default
			}
		}
		if (in_warehouse) continue; // done placing objects

		// add computers; what about monitors?
		float const cheight(0.75*h_val), cwidth(0.44*cheight), cdepth(0.9*cheight); // fixed AR=0.44 to match the texture
		sz[ dim] = 0.5*cdepth;
		sz[!dim] = 0.5*cwidth;

		if (2.1*sz.x < s_sz.x && 2.1*sz.y < s_sz.y && (top_shelf || cheight < shelf_clearance)) { // if it fits in all dims
			unsigned const num_comps(rgen.rand() % (is_house ? 3 : 6)); // 0-5
			C.dim    = dim; C.dir = !dir; // reset dim, flip dir
			C.type   = TYPE_COMPUTER;
			C.shape  = SHAPE_CUBE;
			C.flags |= RO_FLAG_BROKEN; // flag as old

			for (unsigned n = 0; n < num_comps; ++n) {
				gen_xy_pos_for_cube_obj(C, S, sz, cheight, rgen);
				add_if_not_intersecting(C, objects, cubes);
			}
			C.flags &= ~RO_FLAG_BROKEN; // unset old flag
		}
		// add keyboards
		float const kbd_hwidth(0.7*0.5*1.1*2.0*h_val), kbd_depth(0.6*kbd_hwidth), kbd_height(0.06*kbd_hwidth);
		sz[ dim] = 0.5*kbd_depth;
		sz[!dim] = kbd_hwidth;

		if (2.1*sz.x < s_sz.x && 2.1*sz.y < s_sz.y && (top_shelf || kbd_height < shelf_clearance)) { // if it fits in all dims
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
			if (s_sz_min_xy < 4.5*bottle_radius) continue; // shelf not wide/deep enough to add this bottle
			gen_xy_pos_for_round_obj(C, S, bottle_radius, bottle_height, 2.0*bottle_radius, rgen);
			C.set_as_bottle(rgen.rand(), BOTTLE_TYPE_MEDS-1, 1); // all bottle types except for medicine, no_empty=1
			add_if_not_intersecting(C, objects, cubes);
		}
		// add paint cans
		float const pc_height(0.64*z_step), pc_radius(0.28*z_step), pc_edge_spacing(1.1*pc_radius);

		if (2.1*pc_edge_spacing < s_sz_min_xy) { // shelf is wide/deep enough for paint cans
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

		if (s_sz_min_xy > 4.5*spcan_radius) { // add if shelf wide/deep enough
			unsigned const num_spcans(((rgen.rand()%5) < 3) ? (rgen.rand() % 4) : 0); // 0-3, 60% chance
			C.dir   = C.dim = 0;
			C.type  = TYPE_SPRAYCAN;
			C.shape = SHAPE_CYLIN;

			for (unsigned n = 0; n < num_spcans; ++n) {
				gen_xy_pos_for_round_obj(C, S, spcan_radius, spcan_height, 1.5*spcan_radius, rgen);
				set_spraypaint_color(C, rgen);
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
				if (2.1*ball_radius > s_sz_min_xy) continue; // shelf is not wide/deep enough for balls of this type
				gen_xy_pos_for_round_obj(C, S, ball_radius, 2.0*ball_radius, ball_radius, rgen);
				C.item_flags = btype;
				add_if_not_intersecting(C, objects, cubes);
			}
			C.item_flags = 0; // reset for next object type
		}
		// add tape rolls
		float const tape_radius(0.22*z_step), tape_height(TAPE_HEIGHT_TO_RADIUS*tape_radius); // fixed size

		if (s_sz_min_xy > 2.5*tape_radius) { // add if shelf wide/deep enough
			unsigned const num_tapes(((rgen.rand()%4) < 3) ? (rgen.rand() % 4) : 0); // 0-3, 75% chance
			C.dir   = C.dim = 0; // vertical
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

		if (rgen.rand_float() < 0.25 && s_sz_min_xy > 2.5*fl_radius) { // add 25% of the time if shelf wide/deep enough
			C.dir   = C.dim = 0;
			C.type  = TYPE_FLASHLIGHT;
			C.shape = SHAPE_CYLIN;
			C.color = BLACK;
			gen_xy_pos_for_round_obj(C, S, fl_radius, fl_height, 1.25*fl_radius, rgen);
			add_if_not_intersecting(C, objects, cubes);
		}
	} // for s
}

void building_room_geom_t::expand_shelves(room_object_t const &c, bool add_models_mode) {
	cube_t shelves[MAX_SHELVES]; // max number of shelves
	unsigned const num_shelves(get_shelves_for_object(c, shelves));
	get_shelf_objects(c, shelves, num_shelves, (add_models_mode ? objs : expanded_objs), add_models_mode); // models expand to objs; everything else expands to expanded_objs
}

void set_book_id_and_color(room_object_t &obj, rand_gen_t &rgen) {
	obj.obj_id = rgen.rand();
	obj.color  = book_colors[rgen.rand() % NUM_BOOK_COLORS];
}
void set_vase_id_and_color(room_object_t &obj, rand_gen_t &rgen) {
	obj.obj_id = rgen.rand();
	obj.color  = gen_vase_color(rgen);
}

void add_rows_of_vcylinders(room_object_t const &c, cube_t const &region, float radius, float height, float spacing_factor,
	unsigned type, unsigned max_cols, unsigned flags, vect_room_object_t &objects, rand_gen_t &rgen, bool dir=0, bool inv_dim=0)
{
	float const length(region.get_sz_dim(!c.dim)), shelf_depth(region.get_sz_dim(c.dim));
	float const space(spacing_factor*radius), stride(2.0*radius + space);
	unsigned const num_rows(length/stride), num_cols(min((1U + (rgen.rand()%max_cols)), max(1U, unsigned(shelf_depth/stride)))); // round down
	float const row_spacing(length/num_rows), col_spacing(shelf_depth/num_cols);
	point pos;
	pos[ c.dim] = region.d[ c.dim][0] + 0.5*col_spacing;
	pos[!c.dim] = region.d[!c.dim][0] + 0.5*row_spacing;
	pos.z = region.z1();
	if (type == TYPE_PAN || type == TYPE_TCAN || type == TYPE_BUCKET || type == TYPE_VASE) {pos.z += 0.002*c.dz();} // shift up slightly to prevent Z-fighting
	cube_t objc(get_cube_height_radius(pos, radius, height));
	bool const is_drink(type == TYPE_BOTTLE || type == TYPE_DRINK_CAN);
	unsigned rand_id(is_drink ? rgen.rand() : 0);
	if (is_drink) {flags |= RO_FLAG_NO_CONS;} // not consumable

	for (unsigned row = 0; row < num_rows; ++row) {
		cube_t row_obj(objc);

		for (unsigned col = 0; col < num_cols; ++col) {
			if (rgen.rand_float() < 0.75) { // 75% chance
				room_object_t obj(row_obj, type, c.room_id, (c.dim ^ inv_dim), dir, flags, c.light_amt, SHAPE_CYLIN);
				if      (type == TYPE_BOTTLE    ) {obj.set_as_bottle(rand_id, 3, 1);} // 0-3; excludes poison and medicine; should we include medicine?; no_empty=1
				else if (type == TYPE_DRINK_CAN ) {obj.obj_id = (uint16_t)(rand_id & 127); obj.color = WHITE;} // strip off empty bit
				else if (type == TYPE_VASE      ) {set_vase_id_and_color(obj, rgen);} // randomize the vase
				else if (type == TYPE_SPRAYCAN  ) {set_spraypaint_color (obj, rgen);}
				else if (type == TYPE_FLASHLIGHT) {obj.color = BLACK;}
				else if (type == TYPE_CANDLE    ) {obj.color = candle_color;}
				else if (type == TYPE_TCAN      ) {obj.color = tcan_colors[rgen.rand()%NUM_TCAN_COLORS];}
				else if (type == TYPE_LAMP      ) {obj.color = lamp_colors[rgen.rand()%NUM_LAMP_COLORS];}
				else if (type == TYPE_BUCKET    ) {obj.color = LT_GRAY;}
				else if (type == TYPE_PAN       ) {obj.color = DK_GRAY;}
				objects.push_back(obj);
			}
			row_obj.translate_dim(c.dim, col_spacing);
		} // for col
		if (is_drink && rgen.rand_float() < 0.1) {rand_id = rgen.rand();} // 10% chance to update bottle/can type
		objc.translate_dim(!c.dim, row_spacing);
	} // for row
}
unsigned add_row_of_cubes(room_object_t const &c, cube_t const &region, float width, float depth, float height, float spacing_factor, unsigned type,
	unsigned flags, vect_room_object_t &objects, rand_gen_t &rgen, bool dir=0, bool inv_dim=0, unsigned max_stack_height=1, bool no_empty=0)
{
	float const length(region.get_sz_dim(!c.dim)), space(spacing_factor*width), stride(width + space);
	unsigned const objects_start(objects.size()), num_rows(length/stride); // round down
	float const row_spacing(length/num_rows), shelf_depth(region.get_sz_dim(c.dim));
	unsigned num_items(0);
	point pos;
	pos[ c.dim] = region.d[ c.dim][0] + 0.5*shelf_depth;
	pos[!c.dim] = region.d[!c.dim][0] + 0.5*row_spacing;
	pos.z = region.z1();
	cube_t objc(pos);
	objc.z2() += height;
	objc.expand_in_dim( c.dim, 0.5*depth);
	objc.expand_in_dim(!c.dim, 0.5*width);

	for (unsigned row = 0; row < num_rows; ++row) {
		if (no_empty || rgen.rand_float() < 0.75) { // 75% chance
			unsigned const stack_height(1 + ((max_stack_height > 1) ? (rgen.rand() % max_stack_height) : 0));
			cube_t objc_stack(objc);
			
			for (unsigned stack = 0; stack < stack_height; ++stack) {
				room_object_t obj(objc_stack, type, c.room_id, (c.dim ^ inv_dim), dir, flags, c.light_amt);

				if (type == TYPE_BOOK) { // reduce size randomly
					obj.z2() -= 0.5*height*rgen.rand_float();
					obj.expand_in_dim( c.dim, -0.1*width*rgen.rand_float()); // length
					obj.expand_in_dim(!c.dim, -0.2*width*rgen.rand_float()); // width
					set_book_id_and_color(obj, rgen);
				}
				else if (type == TYPE_FOOD_BOX  ) {obj.obj_id  = objects_start;} // unique per section
				else if (type == TYPE_TOASTER   ) {obj.color   = toaster_colors[rgen.rand()%NUM_TOASTER_COLORS];} // random color
				else if (type == TYPE_FOLD_SHIRT) {obj.color   = TSHIRT_COLORS [rgen.rand()%NUM_TSHIRT_COLORS ];} // random color
				else if (type == TYPE_MONITOR   ) {obj.obj_id |= 1;} // off by default; set LSB
				objects.push_back(obj);
				++num_items;
				if (stack+1 == stack_height) break; // done with stack
				objc_stack.translate_dim(2, obj.dz()); // shift stack up
				if (objc_stack.z2() > region.z2()) break; // stack is too tall
				
				if (type == TYPE_LAPTOP) { // add a bit of horizontal jitter
					for (unsigned d = 0; d < 2; ++d) {objc_stack.translate_dim(d, 0.05*objc.get_sz_dim(d)*rgen.signed_rand_float());}
				}
			} // for stack
		}
		objc.translate_dim(!c.dim, row_spacing);
	} // for row
	return num_items;
}
void add_row_of_balls(room_object_t const &c, cube_t const &region, float spacing_factor,
	float floor_spacing, unsigned flags, vect_room_object_t &objects, rand_gen_t &rgen)
{
	unsigned const btype(rgen.rand() % NUM_BALL_TYPES);
	float const radius(ball_types[btype].radius*floor_spacing/(12*8)), diameter(2.0*radius), shelf_depth(region.get_sz_dim(c.dim));
	if (diameter > region.dz() || diameter > shelf_depth) return; // not enough space to fit a ball of this type; often happens with beach balls
	float const length(region.get_sz_dim(!c.dim)), space(spacing_factor*diameter), stride(diameter + space);
	unsigned const num_rows(length/stride); // round down
	float const row_spacing(length/num_rows);
	point pos;
	pos[ c.dim] = region.d[ c.dim][0] + 0.5*shelf_depth;
	pos[!c.dim] = region.d[!c.dim][0] + 0.5*row_spacing;
	pos.z = region.z1() + radius;
	cube_t objc;
	objc.set_from_sphere(pos, radius);

	for (unsigned row = 0; row < num_rows; ++row) {
		if (rgen.rand_float() < 0.75) { // 75% chance
			room_object_t obj(objc, TYPE_LG_BALL, c.room_id, 0, 0, (flags | RO_FLAG_RAND_ROT), c.light_amt, SHAPE_SPHERE);
			obj.item_flags = btype; // no dstate
			objects.push_back(obj);
		}
		objc.translate_dim(!c.dim, row_spacing);
	}
}

void building_room_geom_t::get_shelfrack_objects(room_object_t const &c, vect_room_object_t &objects, bool add_models_mode, bool books_only) {
	if (!c.is_nonempty()) return; // empty - no objects
	cube_t back, top, sides[2], shelves[5];
	unsigned const num_shelves(get_shelf_rack_cubes(c, back, top, sides, shelves));
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_WAS_EXP | RO_FLAG_ON_SRACK), unpowered_flags(flags | RO_FLAG_NO_POWER);
	float const floor_spacing(c.dz()/SHELF_RACK_HEIGHT_FS);
	float const top_shelf_z2(top.is_all_zeros() ? c.z2() : top.z1()); // bottom of the top, if present
	bool const add_food_boxes(!global_building_params.food_box_tids.empty());
	rand_gen_t rgen;
	vect_cube_t cubes; // for placed object overlap tests

	for (unsigned dir = 0; dir < 2; ++dir) { // each shelf has two sides/aisles
		unsigned const rack_id((c.item_flags << 1) + dir);
		rgen.set_state(123*c.obj_id, 207*rack_id+1);
		rgen.rand_mix();
		rgen.rand_mix();
		// 0=boxes items, 1=food, 2=household goods, 3=kitchen, 4=electronics
		int const category(((c.drawer_flags > 0) ? (c.drawer_flags-1) : rgen.rand()) % NUM_RETAIL_CAT); // use drawer_flags to select category (for malls), if nonzero

		for (unsigned n = 0; n < num_shelves; ++n) {
			float const ztop(((n+1) == num_shelves) ? top_shelf_z2 : shelves[n+1].z1());
			cube_t shelf(shelves[n]); // this is the valid shelf area where objects can be placed
			shelf.d[c.dim][!dir] = back.d[c.dim][dir]; // excludes the back
			set_cube_zvals(shelf, shelf.z2(), ztop);
			float const length(shelf.get_sz_dim(!c.dim)), depth(shelf.get_sz_dim(c.dim));
			if (length <= depth) continue; // shouldn't happen
			float const height(shelf.dz()), height_val(min(height, c.dz()/5.0f)); // height_val is for use with objects that were tuned for 5 shelf racks
			cubes.clear();
			rgen.rand_mix(); // make sure to change the random values every shelf even if the items are skipped

			if (category == RETAIL_BOXED) { // boxed items
				if (add_models_mode || books_only) continue; // not model or book
				rand_gen_t rgen2(rgen); // local rgen so that we get the same outcome for either value of add_models_model
				unsigned const num_boxes(rgen2.rand() % 51); // 0-50
				// the call below adds boxes randomly; should they be organized in more orderly rows/columns, or have a more consistent size?
				add_boxes_to_space(c, objects, shelf, cubes, rgen2, num_boxes, 0.375*depth, 0.3*height, 0.8*height, 0, flags); // allow_crates=0
				continue;
			}
			if (category == RETAIL_FOOD) { // food
				if (add_models_mode || books_only) continue; // not model or book
				rand_gen_t rgen2(rgen); // local rgen so that we get the same outcome for either value of add_models_model

				if (n == 0) { // bottom shelf, add pizza; box can't be opened by the player
					float const pbox_width(0.8*depth), pbox_height(0.1*pbox_width), pbox_space(0.2*pbox_width), pbox_stride(pbox_width + pbox_space);
					unsigned const num_boxes(length/pbox_stride), max_stack_height(min(4U, unsigned(height/pbox_height))); // round down
					if (num_boxes == 0 || max_stack_height == 0) continue; // none can fit?
					float const spacing(length/num_boxes), stack_xy_rand(1.0*pbox_space/max_stack_height);
					point pos;
					pos[ c.dim] = shelf.get_center_dim(c.dim);
					pos[!c.dim] = shelf.d[!c.dim][0] + 0.5*spacing;
					cube_t pizza;
					pizza.set_from_point(pos);
					set_cube_zvals(pizza, shelf.z1(), shelf.z1()+pbox_height);
					pizza.expand_by_xy(0.5*pbox_width);

					for (unsigned n = 0; n < num_boxes; ++n) {
						unsigned const stack_sz(rgen2.rand() % max_stack_height);
						cube_t pizza_stacked(pizza);

						for (unsigned N = 0; N < stack_sz; ++N) {
							objects.emplace_back(pizza_stacked, TYPE_PIZZA_BOX, c.room_id, c.dim, dir, flags, c.light_amt);
							point xlate;
							xlate.z = pbox_height; // stacked vertically
							for (unsigned d = 0; d < 2; ++d) {xlate[d] = stack_xy_rand*rgen2.signed_rand_float();}
							pizza_stacked.translate(xlate);
							if (pizza_stacked.intersects(back)) break; // leaning over, hit the back of the shelf, no more stacking
						}
						pizza.translate_dim(!c.dim, spacing);
					} // for n
					continue;
				}
				else if (add_food_boxes && n == num_shelves-1) { // add food boxes on the top shelf
					// will fall through to grouped items case below
				}
				else if (rgen2.rand_float() < 0.65) { // add bottles; these aren't consumable by the player because that would be too powerful
					float const bot_height(height_val*rgen2.rand_uniform(0.7, 0.9)), bot_radius(min(0.25f*depth, bot_height*rgen2.rand_uniform(0.12, 0.18)));
					add_rows_of_vcylinders(c, shelf, bot_radius, bot_height, 0.25, TYPE_BOTTLE, 2, flags, objects, rgen2); // 1-2 columns
					continue;
				}
				else { // add drink cans; should these be grouped into six packs?
					float const can_height(0.48*height_val), can_radius(0.13*height_val); // standard height and radius
					add_rows_of_vcylinders(c, shelf, can_radius, can_height, 0.25, TYPE_DRINK_CAN, 2, flags, objects, rgen2); // 1-3 columns
					continue;
				}
			} // end food
			// items grouped into sections
			unsigned const num_sections(min(unsigned(0.75*length/depth), (3U + (rgen.rand()&3)))); // 3-6
			float const section_width(length/num_sections), section_gap(section_width*rgen.rand_uniform(0.01, 0.05));
			float section_offset(0.0);

			for (unsigned s = 0; s < num_sections; ++s) {
				float const lo_edge(shelf.d[!c.dim][0]);
				cube_t section(shelf);
				section.d[!c.dim][0] = lo_edge +  s   *section_width + section_offset + section_gap;
				section_offset = ((s+1 == num_sections) ? 0.0 : min(0.25f*section_width, (section_width - depth))*rgen.signed_rand_float()); // no offset for last section
				section.d[!c.dim][1] = lo_edge + (s+1)*section_width + section_offset - section_gap;
				rgen.rand_mix(); // make sure to change the random values every section even if the items are skipped
				rand_gen_t rgen2(rgen); // local rgen so that we get the same outcome for either value of add_models_model; no use of rgen beyond this point
				if (books_only && category != RETAIL_HOUSE_GOODS) continue; // not books

				if (category == RETAIL_FOOD) { // food boxes (add_models_mode == 0)
					if (add_models_mode) continue; // not model
					float const fheight(height_val*rgen2.rand_uniform(0.7, 0.9)), fdepth(min(depth, fheight)*rgen2.rand_uniform(0.15, 0.25));
					float const fwidth(min(depth, fheight)*rgen2.rand_uniform(0.6, 0.9));
					add_row_of_cubes(c, section, fwidth, fdepth, fheight, 0.2, TYPE_FOOD_BOX, flags, objects, rgen2, dir);
				}
				else if (category == RETAIL_HOUSE_GOODS) { // houshold goods
					unsigned const type_ix(rgen2.rand() % 9);
					if (add_models_mode && type_ix != 8) continue; // not model
					if (books_only      && type_ix != 6) continue; // not books

					if (type_ix == 0) { // paint cans
						float const oheight(height_val*rgen2.rand_uniform(0.6, 0.8)), radius(min(0.4f*depth, 0.44f*oheight));
						add_rows_of_vcylinders(c, section, radius, oheight, 0.1, TYPE_PAINTCAN, 2, flags, objects, rgen2); // 1-2 columns
					}
					else if (type_ix == 1) { // toilet paper rolls or paper towels
						if (rgen2.rand_bool()) { // paper towels
							float const oheight(min(0.9*height, 1.3*height_val)), radius(min(0.4f*depth, 0.25f*oheight));
							add_rows_of_vcylinders(c, section, radius, oheight, 0.2, TYPE_TPROLL, 2, (flags | RO_FLAG_HAS_EXTRA), objects, rgen2); // 1-2 columns
						}
						else { // TP rolls
							float const oheight(min(0.9*height, 0.48*height_val)), radius(min(0.4f*depth, 0.5f*oheight));
							add_rows_of_vcylinders(c, section, radius, oheight, 0.2, TYPE_TPROLL, 3, flags, objects, rgen2); // 1-3 columns
						}
					}
					else if (type_ix == 2) { // spray paint cans
						float const oheight(height_val*rgen2.rand_uniform(0.75, 0.8)), radius(min(0.25f*depth, 0.17f*oheight));
						add_rows_of_vcylinders(c, section, radius, oheight, 0.2, TYPE_SPRAYCAN, 3, flags, objects, rgen2); // 1-3 columns
					}
					else if (type_ix == 3) { // vases
						float const oheight(height*rgen2.rand_uniform(0.6, 0.9)), radius(min(0.4f*depth, oheight*rgen2.rand_uniform(0.3, 0.5)));
						add_rows_of_vcylinders(c, section, radius, oheight, 0.4, TYPE_VASE, 2, flags, objects, rgen2); // 1-2 columns
					}
					else if (type_ix == 4) { // flashlights
						float const oheight(height_val*rgen2.rand_uniform(0.7, 0.85)), radius(min(0.25f*depth, 0.2f*oheight*rgen2.rand_uniform(0.8, 1.25)));
						add_rows_of_vcylinders(c, section, radius, oheight, 0.25, TYPE_FLASHLIGHT, 3, flags, objects, rgen2); // 1-3 columns
					}
					else if (type_ix == 5) { // candles
						float const oheight(height_val*rgen2.rand_uniform(0.55, 0.7)), radius(min(0.2f*depth, 0.16f*oheight*rgen2.rand_uniform(0.8, 1.25)));
						add_rows_of_vcylinders(c, section, radius, oheight, 0.25, TYPE_CANDLE, 3, flags, objects, rgen2); // 1-3 columns
					}
					else if (type_ix == 6) { // books; should these be stacked or upright?
						add_row_of_cubes(c, section, 0.7*depth, 0.9*depth, 0.12*depth, 0.1, TYPE_BOOK, flags, objects, rgen2, (c.dim ^ bool(dir)), 1, 3); // stacked up to 3 high
					}
					else if (type_ix == 7) { // balls; not dynamic objects
						add_row_of_balls(c, section, 0.25, floor_spacing, flags, objects, rgen2);
					}
					else if (type_ix == 8) { // clothing
						if (!add_models_mode) continue; // not adding models
						if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_FOLD_SHIRT)) continue;
						vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FOLD_SHIRT));
						float const sdepth(depth*rgen2.rand_uniform(0.55, 0.7)), sheight(sdepth*sz.z/sz.y), swidth(sheight*sz.x/sz.z);
						add_row_of_cubes(c, section, swidth, sdepth, sheight, 0.2, TYPE_FOLD_SHIRT, flags, objects, rgen2, dir, 0, 3); // stacked up to 3 high
					}
					// else TYPE_TOY?
				} // end household goods
				else if (category == RETAIL_KITCHEN) { // kitchen
					bool can_have_fire_ext(num_shelves <= 4); // if there are more than 4 shelves, there isn't enough space for fire extinguishers
					unsigned const type_ix(rgen2.rand() % (can_have_fire_ext ? 6 : 5));

					if (type_ix == 0) { // cups
						if (!add_models_mode) continue; // not adding models
						if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_CUP)) continue;
						float const oheight(height_val*rgen2.rand_uniform(0.42, 0.44)), radius(oheight*get_radius_for_square_model(OBJ_MODEL_CUP));
						add_rows_of_vcylinders(c, section, radius, oheight, 1.2, TYPE_CUP, 2, flags, objects, rgen2); // 1-2 columns
					}
					else if (type_ix == 1) { // pans
						if (add_models_mode) continue; // not model
						float const radius(depth*rgen2.rand_uniform(0.3, 0.4)), oheight(radius*rgen2.rand_uniform(0.4, 0.6));
						section.expand_in_dim(!c.dim, -1.2*radius); // make room for handles
							
						if (section.get_sz_dim(!c.dim) > depth) { // add if large enough: 1 column, enough spacing for handle, random dir
							add_rows_of_vcylinders(c, section, radius, oheight, 1.25, TYPE_PAN, 1, flags, objects, rgen2, rgen2.rand_bool());
						}
					}
					else if (type_ix == 2) { // trashcans
						if (add_models_mode) continue; // not model
						unsigned const objs_start(objects.size());
						float const oheight(height*rgen2.rand_uniform(0.7, 0.9)), radius(oheight*rgen2.rand_uniform(0.36, 0.6));
						add_rows_of_vcylinders(c, section, radius, oheight, 0.2, TYPE_TCAN, 1, flags, objects, rgen2); // 1 column
						
						for (auto i = objects.begin()+objs_start; i != objects.end(); ++i) {
							if (rgen2.rand_bool()) {i->shape = SHAPE_CUBE;} // 50% cylinders, 50% cubes
						}
					}
					else if (type_ix == 3) { // microwaves
						if (add_models_mode) continue; // not model; should they be, so that they can be opened without taking an object from the rack?
						float const mheight(height*rgen2.rand_uniform(0.9, 0.95)), mwidth(1.7*mheight), mdepth(min(1.2f*mheight, 0.95f*depth)); // AR=1.7 to match texture
						add_row_of_cubes(c, section, mwidth, mdepth, mheight, 0.1, TYPE_MWAVE, unpowered_flags, objects, rgen2, dir);
					}
					else if (type_ix == 4) { // toasters
						if (!add_models_mode) continue; // not adding models
						if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_TOASTER)) continue;
						vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOASTER));
						float const theight(height_val*0.6), twidth(theight*sz.x/sz.z), tdepth(theight*sz.y/sz.z);
						add_row_of_cubes(c, section, twidth, tdepth, theight, 0.2, TYPE_TOASTER, unpowered_flags, objects, rgen2, dir, 1); // inv_dim=1
					}
					else if (type_ix == 5) { // fire extinguishers
						if (!add_models_mode) continue; // not adding models
						if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_FIRE_EXT)) continue;
						float const fheight(height*rgen2.rand_uniform(0.85, 0.95)), radius(fheight*get_radius_for_square_model(OBJ_MODEL_FIRE_EXT));
						add_rows_of_vcylinders(c, section, radius, fheight, 0.25, TYPE_FIRE_EXT, 1, flags, objects, rgen2); // 1 column
					}
				}
				else if (category == RETAIL_ELECTRONICS) { // electronics
					bool can_have_monitors_lamps(num_shelves <= 3); // if there are more than 3 shelves, there isn't enough space for monitors and lamps
					unsigned const type_ix(rgen2.rand() % (can_have_monitors_lamps ? 4 : 2));

					if (type_ix == 0) { // computers
						if (add_models_mode) continue; // not model
						float const cheight(height*rgen2.rand_uniform(0.9, 0.95)), cwidth(0.44*cheight), cdepth(0.9*cheight); // fixed AR=0.44 to match texture
						add_row_of_cubes(c, section, cwidth, cdepth, cheight, 0.35, TYPE_COMPUTER, unpowered_flags, objects, rgen2, dir);
					}
					else if (type_ix == 1) { // laptops
						if (add_models_mode) continue; // not model
						float const lwidth(depth*rgen2.rand_uniform(0.7, 0.75)), ldepth(0.7*lwidth), lheight(0.06*lwidth);
						add_row_of_cubes(c, section, lwidth, ldepth, lheight, 0.15, TYPE_LAPTOP, flags, objects, rgen2, dir, 0, 2); // stacked up to 3 high
					}
					else if (type_ix == 2) { // monitors
						if (!add_models_mode) continue; // not adding models
						if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_TV)) continue;
						vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
						float const mheight(height*rgen2.rand_uniform(0.85, 0.95)), mwidth(mheight*sz.y/sz.z), mdepth(mheight*sz.x/sz.z);
						add_row_of_cubes(c, section, mwidth, mdepth, mheight, 0.2, TYPE_MONITOR, unpowered_flags, objects, rgen2, dir);
					}
					else if (type_ix == 3) { // lamps
						if (!add_models_mode) continue; // not adding models
						if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_LAMP)) continue;
						float const lheight(height*rgen2.rand_uniform(0.85, 0.95)), radius(lheight*get_radius_for_square_model(OBJ_MODEL_LAMP));
						add_rows_of_vcylinders(c, section, radius, lheight, 0.25, TYPE_LAMP, 1, unpowered_flags, objects, rgen2); // 1 column
					}
				}
			} // for s
		} // for n
	} // for dir
}

void building_room_geom_t::expand_shelfrack(room_object_t const &c) {get_shelfrack_objects(c, expanded_objs);}

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
void place_book(room_object_t &obj, cube_t const &parent, float length, float max_height, unsigned room_id, bool dim, bool dir, rand_gen_t &rgen) {
	float const width(rgen.rand_uniform(0.6, 1.0)*length);
	obj = room_object_t(parent, TYPE_BOOK, room_id, dim, dir);
	set_book_id_and_color(obj, rgen);
	obj.z2() = obj.z1() + min(0.3f*width, rgen.rand_uniform(0.3, 1.0)*max_height);
	set_rand_pos_for_sz(obj, !dim, length, width, rgen);
}

void building_room_geom_t::expand_locker(room_object_t const &c) {
	bool const dim(c.dim), dir(c.dir);
	bool const in_hallway(c.state_flags == RTYPE_HALL), in_locker_room(c.state_flags == RTYPE_LOCKER), in_industrial(c.state_flags == RTYPE_OFFICE);
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP);
	float const wall_thickness(get_locker_wall_thickness(c));
	cube_t interior(c);
	interior.expand_by(-wall_thickness);
	interior.z1() += 0.02*c.dz();
	interior.d[dim][!dir] += (dir ? 1.0 : -1.0)*1.5*wall_thickness;
	float const width(interior.get_sz_dim(!dim)), depth(interior.get_sz_dim(dim)), height(interior.dz()), vspace(c.dz()/0.75);
	unsigned const init_objs_sz(expanded_objs.size());
	// only industrial lockers can have beer
	unsigned const max_bottle_type(in_industrial ? BOTTLE_TYPE_BEER : BOTTLE_TYPE_COKE), max_can_type(in_industrial ? DRINK_CAN_TYPE_BEER : DRINK_CAN_TYPE_COKE);
	rand_gen_t rgen(c.create_rgen());
	bool added_phone(0);

	for (unsigned level = 0; level < 2; ++level) { // {bottom, top}
		cube_t place_area(interior);
		if (level == 0) {place_area.z2() -= (1.0 - LOCKER_BOT_SHELF_HEIGHT)*height + 0.5*wall_thickness;} // lower level; matches drawing code
		else            {place_area.z1() +=        LOCKER_BOT_SHELF_HEIGHT *height + 0.5*wall_thickness;} // upper level; matches drawing code

		if (in_hallway) { // add books; stacked up to 4 high; tag with RO_FLAG_USED to indicate school subject books
			if (add_row_of_cubes(c, place_area, width, depth, 0.15*depth, 0.0, TYPE_BOOK, (flags | RO_FLAG_USED), expanded_objs, rgen, !dir, 1, 4) > 0) {
				float const orig_z2(place_area.z2());
				place_area = expanded_objs.back(); // top book
				set_cube_zvals(place_area, place_area.z2(), orig_z2); // space above the book
			}
		}
		unsigned const num_obj_types(6);
		unsigned const obj_types_hall [num_obj_types] = {TYPE_NONE,    TYPE_BOTTLE, TYPE_DRINK_CAN, TYPE_PHONE,   TYPE_TRASH,    TYPE_TEESHIRT};
		unsigned const obj_types_lroom[num_obj_types] = {TYPE_NONE,    TYPE_BOTTLE, TYPE_DRINK_CAN, TYPE_PHONE,   TYPE_TEESHIRT, TYPE_TEESHIRT};
		unsigned const obj_types_ind  [num_obj_types] = {TYPE_HARDHAT, TYPE_BOTTLE, TYPE_DRINK_CAN, TYPE_HARDHAT, TYPE_TEESHIRT, TYPE_TEESHIRT};
		unsigned const *const obj_types(in_hallway ? obj_types_hall : (in_locker_room ? obj_types_lroom : obj_types_ind));
		unsigned num_sel_obj_types(6);

		if (level == 1) { // teeshirts only on the (larger) bottom level
			while (obj_types[num_sel_obj_types-1] == TYPE_TEESHIRT) {--num_sel_obj_types;}
		}
		for (unsigned n = 0; n < 4; ++n) { // 4 tries for valid placement
			bool bad_item(0);

			switch (obj_types[rgen.rand()%num_sel_obj_types]) {
			case TYPE_NONE:      break; // empty
			case TYPE_BOTTLE:    place_bottle_on_obj(rgen, place_area, expanded_objs, vspace, c.room_id, c.light_amt, max_bottle_type, vect_cube_t(), 1, 1); break; // at_z1=1; trans=1
			case TYPE_DRINK_CAN: place_dcan_on_obj  (rgen, place_area, expanded_objs, vspace, c.room_id, c.light_amt, max_can_type,    vect_cube_t(), 1   ); break; // at_z1=1
			case TYPE_PHONE: {
				if (added_phone) {bad_item = 1; break;} // one phone per locker
				float const phone_length(0.085*height), phone_width(0.45*phone_length);
				if (phone_length > 0.9*place_area.get_sz_dim(dim) || phone_width > 0.9*place_area.get_sz_dim(!dim)) break; // doesn't fit
				room_object_t phone;
				place_phone(phone, place_area, phone_length, phone_width, c.room_id, dim, dir, rgen);
				if (phone.z2() < place_area.z2()) {expanded_objs.push_back(phone);} // add if not too tall
				added_phone = 1;
				break;
			}
			case TYPE_TRASH: {
				float const trash_radius(rgen.rand_uniform(0.15, 0.2)*min(place_area.dx(), place_area.dy()));
				if (place_area.dz() < 2.0*trash_radius) {bad_item = 1; break;} // too tall
				cube_t trash;
				gen_xy_pos_for_round_obj(trash, place_area, trash_radius, 2.0*trash_radius, trash_radius, rgen, 1); // place_at_z1=1
				colorRGBA const color(trash_colors[rgen.rand() % NUM_TRASH_COLORS]);
				expanded_objs.emplace_back(trash, TYPE_TRASH, c.room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, c.light_amt, SHAPE_SPHERE, color);
				set_obj_id(expanded_objs);
				break;
			}
			case TYPE_TEESHIRT: { // vertical
				float const shirt_width(sqrt(width*width + depth*depth)*rgen.rand_uniform(0.9, 0.95)), shirt_len(shirt_width*rgen.rand_uniform(1.1, 1.3));
				if (place_area.dz() < shirt_len) {bad_item = 1; break;} // not enough vertical space for a shirt
				cube_t shirt;
				set_cube_zvals(shirt, place_area.z2()-shirt_len, place_area.z2());
				set_wall_width(shirt, c.get_center_dim( dim), 0.01*shirt_width,  dim); // set thickness
				set_wall_width(shirt, c.get_center_dim(!dim), 0.50*shirt_width, !dim); // set width
				expanded_objs.emplace_back(shirt, TYPE_TEESHIRT, c.room_id, c.dim, c.dir, RO_FLAG_HANGING, c.light_amt, SHAPE_CUBE, gen_teeshirt_color(rgen), rgen.rand());
				break;
			}
			case TYPE_HARDHAT: {
				if (level == 0) {bad_item = 1; break;} // only on the top level

				if (rgen.rand_float() < 0.2) { // tophat 20% of the time
					float const th_radius(min(width, depth)*0.5*rgen.rand_uniform(0.65, 0.8)), th_height(rgen.rand_uniform(1.0, 1.8)*th_radius);
					cube_t that(point(c.xc(), c.yc(), place_area.z1()));
					that.z2() += th_height;
					that.expand_by_xy(th_radius);
					expanded_objs.emplace_back(that, TYPE_TOPHAT, c.room_id, c.dim, c.dir, RO_FLAG_NOCOLL, c.light_amt, SHAPE_CYLIN, BKGRAY);
				}
				else { // hard hat
					float const hh_depth(depth*rgen.rand_uniform(0.75, 0.9)), hh_width(0.8*hh_depth), hh_height(0.56*hh_depth);
					cube_t hhat;
					set_cube_zvals(hhat, place_area.z1(), place_area.z1()+hh_height);
					set_wall_width(hhat, c.get_center_dim( dim), 0.5*hh_depth,  dim); // set thickness
					set_wall_width(hhat, c.get_center_dim(!dim), 0.5*hh_width, !dim); // set width
					expanded_objs.emplace_back(hhat, TYPE_HARDHAT, c.room_id, c.dim, c.dir, RO_FLAG_NOCOLL, c.light_amt, SHAPE_CUBE, hardhat_colors[rgen.rand()%NUM_HARDHAT_COLORS]);
				}
				break;
			}
			default: assert(0);
			} // end switch
			if (!bad_item) break; // success
		} // for n
	} // for level
	if (expanded_objs.size() > init_objs_sz) {invalidate_small_geom();}
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
	float const drawer_dz(drawer.dz());
	// TODO: TYPE_GUN, when ready
	unsigned const type_ix(rgen.rand() % 11); // 0-10
	unsigned const types_dresser [11] = {TYPE_FOLD_SHIRT, TYPE_PAPER,  TYPE_BOX,       TYPE_FOLD_SHIRT, TYPE_BOOK, TYPE_KEY,     TYPE_BOTTLE, TYPE_MONEY,  TYPE_PHONE,  TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_desk    [11] = {TYPE_FLASHLIGHT, TYPE_PAPER,  TYPE_DRINK_CAN, TYPE_STAPLER,    TYPE_BOOK, TYPE_KEY,     TYPE_BOTTLE, TYPE_MONEY,  TYPE_PHONE,  TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_attic   [11] = {TYPE_BOX,        TYPE_PAPER,  TYPE_PEN,       TYPE_PEN,        TYPE_BOOK, TYPE_KEY,     TYPE_BOTTLE, TYPE_BOX,    TYPE_BOOK,   TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_kcabinet[11] = {TYPE_FLASHLIGHT, TYPE_BOX,    TYPE_PEN,       TYPE_PEN,        TYPE_BOOK, TYPE_PLATE,   TYPE_BOTTLE, TYPE_BOTTLE, TYPE_SILVER, TYPE_SPRAYCAN, TYPE_TAPE};
	unsigned const types_fcabinet[11] = {TYPE_BOX,        TYPE_PAPER,  TYPE_PEN,       TYPE_PEN,        TYPE_BOOK, TYPE_STAPLER, TYPE_PAPER,  TYPE_BOOK,   TYPE_TAPE,   TYPE_STAPLER,  TYPE_TAPE};
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
	// these don't combine well with others since they're large horiz cylinders
	bool const is_single_item(obj_type == TYPE_BOTTLE || obj_type == TYPE_DRINK_CAN || obj_type == TYPE_SPRAYCAN || obj_type == TYPE_FLASHLIGHT);
	
	if (item_ix == 0) {stack_z1 = drawer.z1();} // base case
	else if (1 || is_stackable) { // any object can be placed on top of a stackable object
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
		float const length(2.0*drawer_dz), width(0.77*length);

		if (length < 0.9*sz[c.dim] && width < 0.9*sz[!c.dim]) { // if it can fit
			obj = room_object_t(drawer, TYPE_PAPER, c.room_id, c.dim, c.dir); // Note: item_flags/btype not set
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
		float const length(min(1.7f*drawer_dz, 0.9f*sz[dim])), diameter(((type == TYPE_MARKER) ? 0.08 : 0.036)*length);

		if (diameter < min(sz.z, sz[dim])) { // should always be true, unless something got broken (such as calling this on closed drawers)
			obj = room_object_t(drawer, type, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN);
			obj.color = ((obj.type == TYPE_MARKER) ? marker_colors[rgen.rand()&7] : ((obj.type == TYPE_PEN) ? pen_colors[rgen.rand()&3] : pencil_colors[rgen.rand()&1]));
			obj.z2()  = obj.z1() + diameter;
			set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		}
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
		obj.z2()   = obj.z1() + 0.05*sz.z;

		if (c.type != TYPE_DESK || c.is_house()) { // only houses have colored keys
			obj.obj_id = rgen.rand() % NUM_LOCK_COLORS; // lock/key color index
			obj.color  = lock_colors[obj.obj_id];
		}
		break;
	}
	case TYPE_BOTTLE: // bottle, on its side
	{
		bool const dim(c.dim ^ bool(per_drawer_ix & 1)); // random orient, but consistent across the items in the drawer
		float const length(rgen.rand_uniform(0.7, 0.9)*min(((c.type == TYPE_COUNTER) ? 2.7f : 1.8f)*drawer_dz, min(sz.x, sz.y)));
		float const diameter(min(0.8f*sz.z, length*rgen.rand_uniform(0.26, 0.34)));
		obj = room_object_t(drawer, TYPE_BOTTLE, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN);
		obj.set_as_bottle(rgen.rand(), NUM_BOTTLE_TYPES-1, 0, 0, 0, 1); // can be empty; allow_transparent=1
		obj.z2() = obj.z1() + diameter;
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case TYPE_DRINK_CAN: // drink can, vertical
	{
		float const height(min(0.25*c.dz(), 0.9*drawer.dz())), diameter(0.53*height);
		obj = room_object_t(drawer, TYPE_DRINK_CAN, c.room_id, 0, 0, RO_FLAG_RAND_ROT, 1.0, SHAPE_CYLIN); // random orient
		obj.obj_id = rgen.rand();
		obj.z2()   = obj.z1() + height;
		set_rand_pos_for_sz(obj, 0, diameter, diameter, rgen);
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
		float const length(rgen.rand_uniform(0.8, 0.9)*min(1.8f*drawer_dz, min(sz.x, sz.y))), diameter(min(0.8f*sz.z, 0.34f*length));
		obj = room_object_t(drawer, TYPE_SPRAYCAN, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN);
		obj.z2() = obj.z1() + diameter;
		set_spraypaint_color(obj, rgen);
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case TYPE_FLASHLIGHT: // flashlight
	{
		bool const dim(c.dim ^ bool(per_drawer_ix & 1)); // random orient, but consistent across the items in the drawer
		float const length(rgen.rand_uniform(0.8, 0.9)*min(1.8f*drawer_dz, min(sz.x, sz.y))), diameter(min(0.8f*sz.z, 0.4f*length));
		obj = room_object_t(drawer, TYPE_FLASHLIGHT, c.room_id, dim, rgen.rand_bool(), 0, 1.0, SHAPE_CYLIN, BLACK);
		obj.z2() = obj.z1() + diameter;
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
		float const length(min(0.5f*min(sz.x, sz.y), 0.2f*(sz.x + sz.y + drawer_dz))), width(0.2*length), height(min(0.8f*sz.z, 0.3f*length));
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
		obj = room_object_t(shirt, TYPE_TEESHIRT, c.room_id, c.dim, c.dir, 0, 1.0, SHAPE_CUBE, gen_teeshirt_color(rgen));
		break;
	}
	case TYPE_FOLD_SHIRT: // folded shirt
	{
		vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FOLD_SHIRT)); // D, W, H
		float const fs_width(0.75*min(sz[!c.dim], (ssz.y/ssz.x)*sz[c.dim])), fs_length(fs_width*ssz.x/ssz.y), fs_height(fs_width*ssz.z/ssz.y);
		obj = room_object_t(drawer, TYPE_FOLD_SHIRT, c.room_id, c.dim, c.dir, 0, 1.0, SHAPE_CUBE, gen_teeshirt_color(rgen));
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

// shrink c to create {height, depth, length} aspect ratio for dims {2, dim, !dim}
void shrink_to_fit_fixed_ar(cube_t &c, float height, float depth, float length, bool dim) { // dim => depth
	float const orig_z1(c.z1());
	vector3d const sz(c.get_size()), ars((dim ? length : depth), (dim ? depth : length), height);
	vector3d const scale_ars(ars/sz), scale_xyz(scale_ars/scale_ars.get_max_val()), shrink_xyz(0.5*sz*(vector3d(1.0, 1.0, 1.0) - scale_xyz));
	c.expand_by(-shrink_xyz);
	c.translate_dim(2, (orig_z1 - c.z1())); // z1 must remain unchanged
	assert(c.is_strictly_normalized());
}

void building_t::add_box_contents(room_object_t const &box_) {
	room_object_t const box(box_); // deep copy since code below may invalidate box_ reference
	rand_gen_t rgen(box.create_rgen());
	cube_t c(box);
	c.expand_by(-0.01*box.get_size()); // shrink to interior area
	vector3d const sz(c.get_size()), sz_xy(sz.x, sz.y, 0.0);
	bool const in_warehouse(box.in_warehouse()), dim(sz.x < sz.y); // long dim
	float const light_amt(box.light_amt), floor_spacing(get_window_vspace()), base_height(0.2*floor_spacing); // avg shelf height
	unsigned const flags(RO_FLAG_WAS_EXP | RO_FLAG_NOCOLL);
	uint8_t const room_id(box.room_id);
	auto &objs(interior->room_geom->expanded_objs);
	vect_cube_t obj_bcubes;
	unsigned house_objs [] = {TYPE_BOOK,    TYPE_BOTTLE, TYPE_LG_BALL,   TYPE_PAINTCAN, TYPE_SPRAYCAN, TYPE_TPROLL, TYPE_TAPE};
	unsigned office_objs[] = {TYPE_BOOK,    TYPE_BOTTLE, TYPE_DRINK_CAN, TYPE_PAINTCAN, TYPE_SPRAYCAN, TYPE_TPROLL, TYPE_TAPE};
	unsigned whouse_objs[] = {TYPE_MACHINE, TYPE_BOTTLE, TYPE_DRINK_CAN, TYPE_COMPUTER, TYPE_SPRAYCAN, TYPE_TPROLL, TYPE_FOOD_BOX};

	for (unsigned n = 0; n < 10; ++n) { // make up to 10 attempts at placing valid item(s) in this box
		//unsigned const obj_type((n == 9 && !in_warehouse) ? 0 : (rgen.rand()%7)); // place book on last iteration unless this is a warehouse
		unsigned obj_type(0);
		if (is_house)          {obj_type = ((n == 9) ? (unsigned)TYPE_BOOK : house_objs [rgen.rand() % (sizeof(house_objs )/sizeof(unsigned))]);} // place book on last iteration
		else if (in_warehouse) {obj_type =                                   whouse_objs[rgen.rand() % (sizeof(whouse_objs)/sizeof(unsigned))] ;} // warehouse shelf
		else                   {obj_type =                                   office_objs[rgen.rand() % (sizeof(office_objs)/sizeof(unsigned))] ;} // office

		if (sz.z < 0.3*(sz.x + sz.y)) { // food and machines can't fit in short boxes
			if      (obj_type == TYPE_FOOD_BOX) {obj_type = TYPE_BOTTLE   ;}
			else if (obj_type == TYPE_MACHINE ) {obj_type = TYPE_DRINK_CAN;}
		}
		if (obj_type == TYPE_MACHINE) {
			unsigned const item_flags(rgen.rand());
			cube_t obj_bc(c);
			obj_bc.expand_by(-0.05*sz_xy);
			obj_bc.z2() -= 0.05*sz.z;
			objs.emplace_back(obj_bc, TYPE_MACHINE, room_id, rgen.rand_bool(), rgen.rand_bool(), flags, light_amt, SHAPE_CUBE, WHITE, item_flags);
			objs.back().obj_id = rgen.rand();
		}
		else if (obj_type == TYPE_BOOK) { // stacked books; can always be placed
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
		else if (obj_type == TYPE_BOTTLE) { // bottles; not consumable, as this would make things too easy for the player
			float const height(base_height*rgen.rand_uniform(0.4, 0.7)), radius(base_height*rgen.rand_uniform(0.07, 0.11));
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			unsigned const bottle_id(rgen.rand()); // same type for all bottles

			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_BOTTLE, room_id, 0, 0, (flags | RO_FLAG_NO_CONS), light_amt, SHAPE_CYLIN);
				objs.back().set_as_bottle(bottle_id, BOTTLE_TYPE_WINE, 1); // 0-3; excludes poison and medicine; no_empty=1
			}
		}
		else if (obj_type == TYPE_DRINK_CAN) { // cans; not consumable
			float const height(0.3*base_height), radius(0.08*base_height); // standard height and radius
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			uint16_t const can_id(rgen.rand() & 127); // strip off empty bit; same type for all cans

			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_DRINK_CAN, room_id, 0, 0, (flags | RO_FLAG_NO_CONS), light_amt, SHAPE_CYLIN);
				objs.back().obj_id = can_id;
			}
		}
		else if (obj_type == TYPE_LG_BALL) { // ball - only 1
			unsigned const btype(rgen.rand() % (NUM_BALL_TYPES-2)); // no beach balls or tennis balls
			float const radius(ball_types[btype].radius*floor_spacing/(12*8)); // assumes 8 foot floor spacing, and bt.radius in inches
			if (c.min_len() < 2.0*radius) continue; // can't fit any of this item
			point const center(cube_bot_center(c) + radius*plus_z);
			cube_t ball; ball.set_from_sphere(center, radius);
			objs.emplace_back(ball, TYPE_LG_BALL, room_id, 0, 0, flags, light_amt, SHAPE_SPHERE, WHITE, btype);
		}
		else if (obj_type == TYPE_BUCKET) { // empty bucket - only 1; currently unused
			cube_t bucket(cube_bot_center(c));
			bucket.expand_by_xy(min(0.35*sz.z, 0.4*min(sz.x, sz.y))); // set radius
			bucket.z1() += 0.01*sz.z; // fix for bottom Z-fighting
			bucket.z2()  = c.z1() + 0.6*sz.z; // set height; leave space for the handle
			objs.emplace_back(bucket, TYPE_BUCKET, room_id, rgen.rand_bool(), 0, flags, light_amt, SHAPE_CYLIN);
		}
		else if (obj_type == TYPE_PAINTCAN) { // paint cans
			float const height(0.64*base_height), radius(0.28*base_height);
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			
			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_PAINTCAN, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
			}
		}
		else if (obj_type == TYPE_SPRAYCAN) { // spraypaint cans
			float const height(0.55*base_height), radius(0.17*height);
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			unsigned color_ix(rgen.rand()); // random starting color
			
			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_SPRAYCAN, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
				set_spraypaint_color(objs.back(), rgen, &color_ix);
			}
		}
		else if (obj_type == TYPE_TPROLL) { // toilet paper rolls; should we add paper towels as well?
			float const height(0.35*0.18*floor_spacing), radius(0.5*height);
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			
			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_TPROLL, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN);
			}
		}
		else if (obj_type == TYPE_TAPE) { // rolls of tape
			float height(0.032*floor_spacing), radius(height/TAPE_HEIGHT_TO_RADIUS);

			for (unsigned m = 0; m < 2; ++m) {
				if (2.0*radius < 0.95*min(sz.x, sz.y)) break; // size is okay
				height *= 0.9; radius *= 0.9; // can't fit any, try making it smaller
			}
			if (!place_objects_in_box(c, obj_bcubes, radius, height)) continue; // can't fit any of this item
			unsigned color_ix(rgen.rand()); // random color

			for (auto i = obj_bcubes.begin(); i != obj_bcubes.end(); ++i) {
				objs.emplace_back(*i, TYPE_TAPE, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN); // dim/dir don't matter, so use 0
				objs.back().color = tape_colors[color_ix % NUM_TAPE_COLORS];
			}
		}
		else if (obj_type == TYPE_FOOD_BOX) { // food box
			if (global_building_params.food_box_tids.empty()) continue; // not available
			bool const fbdim(!dim); // face in short dim
			float const height(min(0.95f*sz.z, base_height*rgen.rand_uniform(0.8, 1.2)));
			float const depth(height*rgen.rand_uniform(0.15, 0.25)), width(height*rgen.rand_uniform(0.6, 0.9));
			unsigned const nx(round_fp(sz.x/(fbdim ? width : depth))), ny(round_fp(sz.y/(fbdim ? depth : width)));
			if (nx == 0 || ny == 0) continue; // box is too small to place this object
			float const xspace(sz.x/nx), yspace(sz.y/ny), xr(0.45*xspace), yr(0.45*yspace);
			uint16_t const obj_id(rgen.rand()); // set box type
			cube_t fbox;
			set_cube_zvals(fbox, c.z1(), c.z1()+height);

			for (unsigned y = 0; y < ny; ++y) {
				set_wall_width(fbox, (c.y1() + (y+0.5)*yspace), yr, 1);

				for (unsigned x = 0; x < nx; ++x) {
					set_wall_width(fbox, (c.x1() + (x+0.5)*xspace), xr, 0);
					objs.emplace_back(fbox, TYPE_FOOD_BOX, room_id, fbdim, 0, flags, light_amt);
					objs.back().obj_id = obj_id;
				}
			} // for y
		}
		else if (obj_type == TYPE_COMPUTER) { // placeholder for electronics - fit to box size
			float const long_sz(sz[dim]), short_sz(sz[!dim]);
			
			if (short_sz > 1.3*sz.z && short_sz < 0.24*floor_spacing) { // short and wide, but not too large: laptop
				cube_t laptop(c);
				laptop.z2() -= 0.4*sz.z; // add some padding at the top and sides
				laptop.expand_by(-0.2*sz_xy);
				shrink_to_fit_fixed_ar(laptop, 0.06, 0.7, 1.0, !dim);
				objs.emplace_back(laptop, TYPE_LAPTOP, room_id, !dim, rgen.rand_bool(), flags, light_amt);
			}
			else if ((long_sz > 1.3*short_sz && sz.z > 1.1*short_sz) || short_sz < 0.26*floor_spacing) { // tall and thin, or small: computer
				cube_t computer(c);
				computer.z2() -= 0.1*sz.z; // add some padding at the top and sides
				computer.expand_by(-0.1*sz_xy);
				shrink_to_fit_fixed_ar(computer, 0.75, 0.9, 0.44, dim);
				objs.emplace_back(computer, TYPE_COMPUTER, room_id, dim, rgen.rand_bool(), flags, light_amt);
			}
			else { // closer to equal dims: microwave
				cube_t mwave(c);
				mwave.z2() -= 0.05*sz.z; // add some padding at the top and sides
				mwave.expand_by(-0.05*sz_xy);
				shrink_to_fit_fixed_ar(mwave, 1.0, 1.2, 1.7, !dim);
				objs.emplace_back(mwave, TYPE_MWAVE, room_id, !dim, rgen.rand_bool(), flags, light_amt);
			}
		}
		else if (obj_type == TYPE_BOX) { // nested box (not currently enabled)
			cube_t box2(c);
			box2.expand_by(-0.04*sz); // shrink by 10% in all dims (original 2% + new 8%)
			objs.emplace_back(box2, TYPE_BOX, room_id, box.dim, box.dir, (flags | (box.flags & ~RO_FLAG_OPEN)), light_amt);
			objs.back().obj_id += rgen.rand();
		}
		else if (obj_type == TYPE_NONE) {continue;} // empty box
		else {assert(0);}
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
		case 0: ret.item = custom_item_t("bag of chips",              0.5,  0.1, 0.1); break; // +10% health
		case 1: ret.item = custom_item_t("sweaty gym clothes",        2.0,  0.4); break;
		case 2: ret.item = custom_item_t("day old pizza",             0.0,  0.6); break;
		case 3: ret.item = custom_item_t("case of CDs",               5.0,  0.2); break;
		case 4: ret.item = custom_item_t("handgun (with no ammo)",  150.0,  1.0); break;
		case 5: ret.item = custom_item_t("locked briefcase",        200.0, 10.0); break;
		case 6: ret.item = custom_item_t("stuffed animal",            5.0,  0.3); break;
		case 7: ret.item = custom_item_t("pair of shoes",            40.0,  2.0); break;
		case 8: ret.item = custom_item_t("cooler full of beer",      60.0, 20.0); break;
		case 9: ret.item = custom_item_t("box of Girl Scout cookies", 8.0,  1.0, 0.2); break; // +20% health
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
	case TYPE_SHELFRACK: expand_shelfrack(c); break;
	case TYPE_WINE_RACK: expand_wine_rack(c); break;
	case TYPE_CABINET: case TYPE_COUNTER: case TYPE_KSINK: case TYPE_VANITY: expand_cabinet(c); break;
	case TYPE_MED_CAB:   expand_med_cab(c); break;
	case TYPE_LOCKER:    expand_locker (c); break;
	case TYPE_BRK_PANEL: expand_breaker_panel(c, building); break;
	//case TYPE_TCAN:      expand_trashcan(c); break;
	default: assert(0); // not a supported expand type
	}
	if (c.type == TYPE_CLOSET) {maybe_spawn_spider_in_drawer(c, c, 0, building.get_window_vspace(), 1);} // spawn spider when first opened
	c.flags |= RO_FLAG_EXPANDED; // flag as expanded
	return 1;
}

