// 3D World - Building Interior Room Geometry Drawing
// by Frank Gennari 7/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"


float get_lamp_width_scale();
vect_cube_t &get_temp_cubes();

cube_t get_closet_interior_space(room_object_t const &c, cube_t const cubes[5]) {
	cube_t interior(c);
	if (!cubes[1].is_all_zeros()) {interior.d[!c.dim][0] = cubes[1].d[!c.dim][1];} // left  side (if wall exists)
	if (!cubes[3].is_all_zeros()) {interior.d[!c.dim][1] = cubes[3].d[!c.dim][0];} // right side (if wall exists)
	interior.d[ c.dim][c.dir] = cubes[2].d[c.dim][!c.dir]; // set to inside of front wall
	interior.expand_by_xy(-0.02*interior.min_len()); // shrink slightly to clip off the area taken by the wall trim where objects shouldn't be placed
	assert(interior.is_strictly_normalized());
	return interior;
}

void building_room_geom_t::add_closet_objects(room_object_t const &c, vector<room_object_t> &objects) {
	cube_t ccubes[5]; // only used to get interior space
	get_closet_cubes(c, ccubes);
	cube_t const interior(get_closet_interior_space(c, ccubes));
	float const depth(interior.get_sz_dim(c.dim)), box_sz(0.25*depth), window_vspacing(c.dz()*(1.0 + FLOOR_THICK_VAL_HOUSE));
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	unsigned const num_boxes((rgen.rand()%3) + (rgen.rand()%4)); // 0-5
	vect_cube_t &cubes(get_temp_cubes());
	room_object_t C(c);
	C.flags = (RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP); // Note: also clears open flag
	C.type  = TYPE_BOX;
	vector3d sz;
	point center;

	for (unsigned n = 0; n < num_boxes; ++n) {
		for (unsigned d = 0; d < 2; ++d) {
			sz    [d] = box_sz*rgen.rand_uniform(0.5, 1.0); // x,y half width
			center[d] = rgen.rand_uniform(interior.d[d][0]+sz[d], interior.d[d][1]-sz[d]); // randomly placed within the bounds of the closet
		}
		C.set_from_point(center);
		set_cube_zvals(C, interior.z1(), (interior.z1() + box_sz*rgen.rand_uniform(0.8, 1.5)));
		C.expand_by_xy(sz);
		if (has_bcube_int(C, cubes)) continue; // intersects - just skip it, don't try another placement
		C.color  = gen_box_color(rgen);
		C.dim    = rgen.rand_bool();
		C.dir    = rgen.rand_bool();
		C.obj_id = n+1; // make it unique so that contents are unique
		objects.push_back(C);
		cubes.push_back(C);
	} // for n
	if (!c.is_small_closet() && rgen.rand_bool()) { // maybe add a lamp in the closet if it's large
		float const height(0.25*window_vspacing), width(height*get_lamp_width_scale()), radius(0.5*width);

		if (width > 0.0 && width < 0.9*min(interior.dx(), interior.dy())) { // check if lamp model is valid and lamp fits in closet
			point center(0.0, 0.0, interior.z1());

			for (unsigned n = 0; n < 4; ++n) { // make up to 4 attempts to place a lamp
				for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(interior.d[d][0]+radius, interior.d[d][1]-radius);}
				cube_t lamp(get_cube_height_radius(center, radius, height));

				if (!has_bcube_int(lamp, cubes)) { // check for intersection with boxes
					objects.emplace_back(lamp, TYPE_LAMP, c.room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_INTERIOR | RO_FLAG_WAS_EXP), 0.0, SHAPE_CYLIN, lamp_colors[rgen.rand()%NUM_LAMP_COLORS]);
					break;
				}
			} // for n
		}
	}
	// add hanger rod
	float const hr_radius(0.015*window_vspacing);
	room_object_t hanger_rod(interior, TYPE_HANGER_ROD, c.room_id, c.dim, c.dir, (RO_FLAG_NOCOLL | RO_FLAG_INTERIOR));
	hanger_rod.z1() = c.z1() + 0.8*window_vspacing;
	hanger_rod.z2() = hanger_rod.z1() + 2.0*hr_radius;
	set_wall_width(hanger_rod, (0.45*c.d[c.dim][c.dir] + 0.55*c.d[c.dim][!c.dir]), hr_radius, c.dim); // move slightly toward the back
	objects.push_back(hanger_rod);
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
void add_if_not_intersecting(room_object_t const &obj, vector<room_object_t> &objects, vect_cube_t &cubes) {
	if (!has_bcube_int(obj, cubes)) {objects.push_back(obj); cubes.push_back(obj);}
}
void gen_xy_pos_for_cube_obj(cube_t &C, cube_t const &S, vector3d const &sz, float height, rand_gen_t &rgen) {
	point center;
	for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform(S.d[d][0]+sz[d], S.d[d][1]-sz[d]);} // randomly placed within the bounds of the shelf
	C.set_from_point(center);
	set_cube_zvals(C, S.z2(), S.z2()+height);
	C.expand_by_xy(sz);
}
void gen_xy_pos_for_round_obj(cube_t &C, cube_t const &S, float radius, float height, float spacing, rand_gen_t &rgen) {
	point center;
	for (unsigned d = 0; d < 2; ++d) {center[d] = rgen.rand_uniform((S.d[d][0] + spacing), (S.d[d][1] - spacing));} // place at least spacing from edge
	C.set_from_sphere(center, radius);
	set_cube_zvals(C, S.z2(), S.z2()+height);
}

void building_room_geom_t::get_shelf_objects(room_object_t const &c_in, cube_t const shelves[4], unsigned num_shelves, vector<room_object_t> &objects) {
	room_object_t c(c_in);
	c.flags |= RO_FLAG_WAS_EXP;
	bool const is_house(c.is_house());
	vector3d const c_sz(c.get_size());
	float const dz(c_sz.z), width(c_sz[c.dim]), thickness(0.02*dz), bracket_thickness(0.75*thickness);
	float const z_step(dz/(num_shelves + 1)), shelf_clearance(z_step - thickness - bracket_thickness), sz_scale(is_house ? 0.5 : 1.0);
	rand_gen_t rgen;
	c.set_rand_gen_state(rgen);
	vect_cube_t &cubes(get_temp_cubes());

	for (unsigned s = 0; s < num_shelves; ++s) {
		cube_t const &S(shelves[s]);
		room_object_t C(c);
		vector3d sz;
		cubes.clear();
		// add crates/boxes
		unsigned const num_boxes(rgen.rand() % (is_house ? 8 : 13)); // 0-12

		for (unsigned n = 0; n < num_boxes; ++n) {
			point center;

			for (unsigned d = 0; d < 2; ++d) {
				sz[d] = 0.5*width*sz_scale*(is_house ? 1.5 : 1.0)*rgen.rand_uniform(0.45, 0.8); // x,y half width
				center[d] = rgen.rand_uniform(S.d[d][0]+sz[d], S.d[d][1]-sz[d]); // randomly placed within the bounds of the shelf
			}
			C.set_from_point(center);
			set_cube_zvals(C, S.z2(), (S.z2() + shelf_clearance*sz_scale*rgen.rand_uniform(0.4, 0.98)));
			C.expand_by_xy(sz);
			if (has_bcube_int(C, cubes)) continue; // intersects - just skip it, don't try another placement
			C.color  = gen_box_color(rgen);
			C.dim    = c.dim ^ bool(rgen.rand()&3) ^ 1; // make the box label face outside 75% of the time
			C.obj_id = rgen.rand(); // used to select crate texture and box contents
			C.type   = ((is_house || rgen.rand_bool()) ? TYPE_BOX : TYPE_CRATE);
			objects.push_back(C);
			cubes.push_back(C);
		} // for n
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
			C.set_as_bottle(rgen.rand() & 127); // no empty bottles - don't set 7th bit
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
	} // for s
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
		obj.set_as_bottle(rgen.rand());
		obj.z2() = (obj.z1() + diameter);
		set_rand_pos_for_sz(obj, dim, length, diameter, rgen);
		break;
	}
	case 7: // money
	{
		float const length(0.5*sz.z), width(2.35*length);

		if (length < 0.9*sz[c.dim] && width < 0.9*sz[!c.dim]) { // if it can fit
			obj = room_object_t(drawer, TYPE_MONEY, c.room_id, c.dim, c.dir);
			obj.z2() = (obj.z1() + 0.01*length*((rgen.rand()%20) + 1)); // 1-20 bills
			set_rand_pos_for_sz(obj, c.dim, length, width, rgen);
		}
		break;
	}
	case 8: // phone
	{
		bool const dim(c.dim ^ rgen.rand_bool()); // random orient
		float const length(1.1*sz.z), width(0.45*length);

		if (length < 0.9*sz[dim] && width < 0.9*sz[!dim]) { // if it can fit
			unsigned const NUM_PHONE_COLORS = 7; // for the case
			colorRGBA const phone_colors[NUM_PHONE_COLORS] = {WHITE, GRAY, DK_GRAY, GRAY_BLACK, BLUE, RED, PINK};
			obj = room_object_t(drawer, TYPE_PHONE, c.room_id, dim, c.dir);
			obj.color = phone_colors[rgen.rand() % NUM_PHONE_COLORS];
			obj.z2()  = (obj.z1() + 0.06*length);
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
	case 10: // empty
		break;
	}
	obj.flags    |= (RO_FLAG_WAS_EXP | RO_FLAG_NOCOLL);
	obj.light_amt = c.light_amt;
	return obj;
}

void building_t::add_box_contents(room_object_t const &box) {
	rand_gen_t rgen;
	box.set_rand_gen_state(rgen);
	cube_t c(box);
	c.expand_by(-0.01*box.get_size()); // shrink to interior area
	vector3d const sz(c.get_size());
	bool const dim(sz.x < sz.y); // long dim
	float const base_height(0.2*get_window_vspace());
	bool was_added(0);

	for (unsigned n = 0; n < 10; ++n) { // make up to 10 attempts at placing valid item(s) in this box
		unsigned const obj_type(/*rgen.rand()%6*/0); // {book, bottles, ball, paint can, spraypaint, toilet paper}

		if (obj_type == 0) { // books
			unsigned const room_id(box.room_id), num_books(1 + (rgen.rand()&3)); // 1-4 books
			float cur_zval(c.z1());

			// Note: the code below may invalidate the reference to box, so we can't use it after this point
			for (unsigned n = 0; n < num_books; ++n) {
				float const length(rgen.rand_uniform(0.7, 0.95)*min(sz[dim], 2.0f*sz[!dim])), width(min(rgen.rand_uniform(0.6, 1.0)*length, 0.95f*sz[!dim]));
				room_object_t obj(c, TYPE_BOOK, room_id, !dim, rgen.rand_bool(), (RO_FLAG_WAS_EXP | RO_FLAG_NOCOLL));
				obj.obj_id = rgen.rand();
				obj.color  = book_colors[rgen.rand() % NUM_BOOK_COLORS];
				set_cube_zvals(obj, cur_zval, (cur_zval + min(0.3f*width, rgen.rand_uniform(0.1, 0.2)*sz.z)));
				if (obj.z2() > c.z2()) break; // book doesn't fit - the stack is too tall
				set_rand_pos_for_sz(obj, dim, length, width, rgen);
				interior->room_geom->expanded_objs.push_back(obj);
				was_added = 1;
				cur_zval  = obj.z2();
			} // for n
			break; // done
		}
		else if (obj_type == 1) { // bottles
			// TODO: 0.4-0.7 base_height
		}
		else if (obj_type == 2) { // ball
			// TODO: radius 0.2*base_height
		}
		else if (obj_type == 3) { // paint can
			// TODO: 0.64*base_height
		}
		else if (obj_type == 4) { // spraypaint
			// TODO: 0.55*base_height
		}
		else if (obj_type == 5) { // toilet paper
			// TODO: ~0.4*base_height
		}
	} // for n
	if (was_added) {interior->room_geom->clear_static_small_vbos();}
}

