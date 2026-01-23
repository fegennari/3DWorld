// 3D World - Building Restaurants
// by Frank Gennari 01/09/2026

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

extern object_model_loader_t building_obj_model_loader;

void create_wall(cube_t &wall, bool dim, float wall_pos, float fc_thick, float wall_half_thick, float wall_edge_spacing);
float get_radius_for_square_model(unsigned model_id);
colorRGBA get_stain_color(rand_gen_t &rgen, bool is_food=0);


void building_t::create_restaurant_floorplan(unsigned part_id, rand_gen_t &rgen) {
	// Note: exterior doors have not yet been placed
	assert(part_id < parts.size());
	cube_t const &part(parts[part_id]);
	vector<room_t> &rooms(interior->rooms);
	assert(rooms.empty()); // must call this first
	// divide into main restaurant area vs. {kitchen, bathrooms, storage}
	bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // room split dim/dir
	float const h_space(get_hspacing_for_part(part, dim)), h_border(get_window_h_border()), part_sz(part.get_sz_dim(dim));
	float const floor_spacing(get_window_vspace()), fc_thick(get_fc_thickness());
	float const wall_thick(get_wall_thickness()), wall_half_thick(0.5*wall_thick), wall_edge_spacing(0.05*wall_thick);
	float split_pos(0.0);
	interior->restaurant_orient = 2*dim + dir;
	
	for (unsigned n = 0; n < 100; ++n) { // 100 tries for a valid split pos
		split_pos = part.d[dim][dir] + rgen.rand_uniform(0.25, 0.35)*part_sz*(dir ? -1.0 : 1.0);
		if (!is_val_inside_window(part, dim, split_pos, h_space, h_border)) break;
	}
	cube_t wall(part);
	create_wall(wall, dim, split_pos, fc_thick, wall_half_thick, wall_edge_spacing);
	cube_t wall_lo(wall), wall_hi(wall);
	wall_lo.z2() = wall_hi.z1() = part.z1() + floor_spacing - fc_thick; // top wall contains the space that would be between the upper ceiling and lower floor
	interior->walls[dim].push_back(wall_hi);
	cube_t main_room(part), side_area(part);
	main_room.d[dim][dir] = side_area.d[dim][!dir] = split_pos;
	add_assigned_room(main_room, part_id, RTYPE_RESTAURANT); // num_lights will be calculated later
	// split side room into {kitchen, men's room, women's room, and maybe storage}
	int const num_side_windows(get_num_windows_on_side(part, !dim)); // in other dim; typically 5-8 total windows and 3-5 kitchen windows
	assert(num_side_windows >= 3);
	bool const add_storage(num_side_windows >= 7), br_side(rgen.rand_bool()), mw_side(rgen.rand_bool());
	float const wspace(part.get_sz_dim(!dim)/num_side_windows), window_step((br_side ? -1.0 : 1.0)*wspace);
	float const br_split(side_area.d[!dim][br_side] + window_step); // split point between men's and women's bathrooms
	float const k_br_split(br_split + window_step); // split point between kitchens and bathrooms
	float const k_s_split(side_area.d[!dim][!br_side] - window_step); // split point between kitchen and storage, or end of kitchen
	cube_t kitchen(side_area), br1(side_area), br2(side_area), storage(side_area);
	br1.d[!dim][!br_side] = br2    .d[!dim][br_side] = br_split;
	br2.d[!dim][!br_side] = kitchen.d[!dim][br_side] = k_br_split;
	if (add_storage) {kitchen.d[!dim][!br_side] = storage.d[!dim][br_side] = k_s_split;}
	add_assigned_room(kitchen, part_id, RTYPE_KITCHEN); // num_lights will be calculated later
	add_assigned_room(br1,     part_id, (mw_side ? RTYPE_MENS : RTYPE_WOMENS));
	add_assigned_room(br2,     part_id, (mw_side ? RTYPE_WOMENS : RTYPE_MENS));
	if (add_storage) {add_assigned_room(storage, part_id, RTYPE_STORAGE);}
	// add doors to wall_lo; all doors are unlocked since restaurants have no keys
	bool storage_conn_to_kitchen(1);
	unsigned const num_wall_doors((add_storage && !storage_conn_to_kitchen) ? 4U : 3U);
	float const doorway_width(get_nominal_doorway_width()), doorway_hwidth(0.5*doorway_width), edge_pad(doorway_hwidth + wall_thick);
	cube_t room_bcs[4] = {br1, br2, kitchen, storage};
	
	for (unsigned n = 0; n < num_wall_doors; ++n) {
		bool const is_br(n < 2);
		float const v1(room_bcs[n].d[!dim][0] + edge_pad), v2(room_bcs[n].d[!dim][1] - edge_pad);
		assert(v1 < v2); // assumes window is wider than door
		float const door_pos(rgen.rand_uniform(v1, v2));
		insert_door_in_wall_and_add_seg(wall_lo, (door_pos - doorway_hwidth), (door_pos + doorway_hwidth), !dim, dir, !br_side, is_br, 1); // opens into room; unlocked
		interior->door_stacks.back().set_mult_floor(); // counts as multi-floor (for drawing top edge)
		if (is_br) {interior->doors.back().rtype = ((bool(n) ^ mw_side) ? RTYPE_MENS : RTYPE_WOMENS);}
	}
	interior->walls[dim].push_back(wall_lo);

	// add walls separating rooms
	for (unsigned n = 0; n < (add_storage ? 3U : 2U); ++n) {
		wall = side_area;
		create_wall(wall, !dim, room_bcs[n].d[!dim][!br_side], fc_thick, wall_half_thick, wall_edge_spacing);

		if (storage_conn_to_kitchen && n == 2) { // storage room door
			cube_t wall_hi(wall);
			wall.z2() = wall_hi.z1() = wall_lo.z2(); // must split into lo and hi parts again to add the door
			interior->walls[!dim].push_back(wall_hi);
			float const v1(storage.d[dim][0] + edge_pad), v2(storage.d[dim][1] - edge_pad);
			assert(v1 < v2); // assumes room is wider than door
			float const door_pos(rgen.rand_uniform(v1, v2));
			insert_door_in_wall_and_add_seg(wall, (door_pos - doorway_hwidth), (door_pos + doorway_hwidth), dim, !br_side, 0, 0, 1); // opens into storage; unlocked
			interior->door_stacks.back().set_mult_floor(); // counts as multi-floor (for drawing top edge)
		}
		interior->walls[!dim].push_back(wall);
	} // for n
	// all main part restaurant rooms are a single floor
	for (room_t &room : rooms) {room.is_single_floor = 1;}
}

void building_t::add_restaurant_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float light_amt, unsigned light_nx, unsigned light_ny) {
	cube_t const room_bounds(get_walkable_room_bounds(room));
	cube_t place_area(room_bounds);
	place_area.expand_by(-0.25*get_wall_thickness()); // common spacing to wall
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size());
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	bool added_desk(0);

	for (tquad_with_ix_t const &door : doors) {
		cube_t bc(door.get_bcube());
		bc.expand_by_xy(wall_thickness); // make nonzero area
		if (!bc.intersects_no_adj(room)) continue; // not for this room
		bool const ddim(bc.dy() < bc.dx()), ddir(bc.get_center_dim(ddim) < room.get_center_dim(ddim)); // dir into room
		float const door_edge(bc.d[ddim][ddir]);

		if (!added_desk) { // place a small desk/table/podium by the front (first) door
			float const dlo(bc.d[!ddim][0] - room.d[!ddim][0]), dhi(room.d[!ddim][1] - bc.d[!ddim][1]);
			bool tside(0);
			if      (dlo < 0.5*dhi) {tside = 1;} // door close to lo wall, put desk on hi side
			else if (dhi < 0.5*dlo) {tside = 0;} // door close to hi wall, put desk on lo side
			else {tside = rgen.rand_bool();} // door not near a wall, choose a random side
			float const table_sz(0.12*floor_spacing);
			cube_t table;
			set_cube_zvals(table, zval, zval+0.4*floor_spacing);
			set_wall_width(table, (door_edge + (ddir ? 1.0 : -1.0)*1.5*table_sz), table_sz, ddim);
			set_wall_width(table, (bc.d[!ddim][tside] + (tside ? 1.0 : -1.0)*1.5*table_sz), table_sz, !ddim);
			objs.emplace_back(table, TYPE_TABLE, room_id, !ddim, !tside, 0, light_amt, SHAPE_TALL, WHITE); // wood

			if (building_obj_model_loader.is_model_valid(OBJ_MODEL_BAR_STOOL)) { // add a stool
				float const chair_height(0.45*floor_spacing), chair_hwidth(chair_height*get_radius_for_square_model(OBJ_MODEL_BAR_STOOL));
				cube_t chair;
				set_cube_zvals(chair, zval, (zval + chair_height));
				set_wall_width(chair, table.get_center_dim(ddim), chair_hwidth, ddim);
				set_wall_width(chair, (table.d[!ddim][tside] + (tside ? 1.0 : -1.0)*1.5*chair_hwidth), chair_hwidth, !ddim);
				objs.emplace_back(chair, TYPE_BAR_STOOL, room_id, !ddim, !tside, 0, light_amt, SHAPE_CUBE);
			}
			added_desk = 1;
		}
		// add rugs (door mats?) by the door(s)
		float const rug_hwidth(rgen.rand_uniform(0.18, 0.22)*floor_spacing), rug_hlen(rgen.rand_uniform(1.5, 1.7)*rug_hwidth);
		cube_t rug;
		set_cube_zvals(rug, zval, (zval + get_rug_thickness()));
		set_wall_width(rug, bc.get_center_dim(!ddim), rug_hlen, !ddim);
		set_wall_width(rug, (door_edge + (ddir ? 1.0 : -1.0)*rgen.rand_uniform(1.1, 1.3)*rug_hwidth), rug_hwidth, ddim); // move away from the door
		objs.emplace_back(rug, TYPE_RUG, room_id, 0, 0, RO_FLAG_NOCOLL, light_amt);
		objs.back().obj_id = uint16_t(11*objs.size() + 17*mat_ix); // determines rug texture
	} // for door
	bool const plastic_tc(0); // custom material?
	fill_room_with_tables_and_chairs(rgen, room, zval, room_id, light_amt, objs_start, plastic_tc);
	unsigned const num_wine_racks(rgen.rand_bool() ? 2 : 1);
	for (unsigned n = 0; n < num_wine_racks; ++n) {add_wine_rack(rgen, room, zval, room_id, light_amt, objs_start);}
	if (rgen.rand_bool()) {add_fishtank_to_room(rgen, room, zval, room_id, light_amt, objs_start, place_area);}
	unsigned const num_plants(6 + (rgen.rand() & 5)); // 6-10
	add_plants_to_room(rgen, room, zval, room_id, light_amt, objs_start, num_plants);
	float const ceil_zval(room.z2() - get_fc_thickness());

	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_CEIL_FAN)) { // add ceiling fans between the ceiling lights
		float const dx(room.dx()), dy(room.dy()), xstep(dx/light_nx), ystep(dy/light_ny);
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_CEIL_FAN)); // D, W, H
		float const diameter(0.75*min(min(xstep, ystep), floor_spacing)), height(diameter*sz.z/sz.y); // assumes width = depth = diameter
		bool skip(1);

		for (unsigned y = 1; y < light_ny; ++y) {
			for (unsigned x = 1; x < light_nx; ++x) {
				skip ^= 1;
				if (skip) continue; // only add every other fan
				cube_t fan(point((room.x1() + x*xstep), (room.y1() + y*ystep), ceil_zval));
				fan.expand_by_xy(0.5*diameter);
				fan.z1() -= height;
				unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_IN_FACTORY); // set factory flag so that fans rotate
				if (rgen.rand_float() < 0.5) {flags |= RO_FLAG_ROTATING;} // make fan rotate when turned on 50% of the time
				objs.emplace_back(fan, TYPE_CEIL_FAN, room_id, 0, 0, flags, light_amt, SHAPE_CYLIN, WHITE);
			} // for x
		} // for y
	}
	// add wooden ceiling beams
	float const beam_hwidth(0.025*floor_spacing), beam_height(3.0*beam_hwidth);
	bool const pri_dim(rgen.rand_bool());

	for (unsigned d = 0; d < 2; ++d) {
		bool const is_pri(bool(d) == pri_dim);
		unsigned const num_beams((d ? light_nx : light_ny) + 2); // includes beams along room edges
		float const beam_spacing(room.get_sz_dim(!d)/(num_beams - 2));
		float const beams_start(room.d[!d][0] - 0.5*beam_spacing);
		float const clamp_lo(room_bounds.d[!d][0] + beam_hwidth), clamp_hi(room_bounds.d[!d][1] - beam_hwidth); // clamp to room walls
		cube_t beam(room_bounds); // spans the entire room, excluding the walls
		set_cube_zvals(beam, (ceil_zval - (is_pri ? 1.2 : 1.0)*beam_height), ceil_zval); // primary dim beams are slightly lower to prevent Z-fighting
		if (!is_pri) {beam.expand_in_dim(0, -2.0*beam_hwidth);} // remove overlaps at ends

		for (unsigned n = 0; n < num_beams; ++n) {
			unsigned skip_faces(EF_Z2 | get_skip_mask_for_xy(d));
			if (n   == 0        ) {skip_faces |= ~get_face_mask(!d, 0);} // skip face against the wall
			if (n+1 == num_beams) {skip_faces |= ~get_face_mask(!d, 1);} // skip face against the wall
			set_wall_width(beam, max(clamp_lo, min(clamp_hi, (beams_start + n*beam_spacing))), beam_hwidth, !d);
			objs.emplace_back(beam, TYPE_IBEAM, room_id, d, 0, (RO_FLAG_NOCOLL | RO_FLAG_IS_HOUSE), light_amt, SHAPE_CUBE, WHITE, skip_faces); // house == wood
		}
	} // for d
	if (rgen.rand_bool()) {add_wall_tv(rgen, room, zval, room_id, light_amt, objs_start);} // maybe TV for sports bar
	// add additional pictures, likely only on the wall separating dining from kitchen and bathrooms
	hang_pictures_whiteboard_chalkboard_in_room(rgen, room, zval, room_id, light_amt, objs_start, 0, 0);
}

void building_t::make_restaurant_light(room_object_t &light) {
	assert(has_room_geom());
	float const radius(0.25*(light.dx() + light.dy())), rod_radius(0.08*radius), rod_len(0.75*light.dz());
	cube_t rod;
	set_cube_zvals(rod, (light.z2() - 1.2*rod_len), light.z2());
	light.translate_dim(2, -rod_len);
	for (unsigned d = 0; d < 2; ++d) {set_wall_width(rod, light.get_center_dim(d), rod_radius, d);}
	interior->room_geom->objs.emplace_back(rod, TYPE_METAL_BAR, light.room_id, 0, 1, RO_FLAG_NOCOLL, light.light_amt, SHAPE_CYLIN, BKGRAY); // vertical
}

template<typename T> void remove_if_intersects(vector<T> &objs, cube_t const &c) {
	auto i(objs.begin()), o(i);

	for (; i != objs.end(); ++i) {
		if (!i->intersects(c)) {*(o++) = *i;} // keep if not this store
	}
	objs.erase(o, objs.end());
}

void building_t::add_mall_restaurant_objs(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, bool dim, bool dir,
	bool no_doorway, float light_amt, cube_t &div_wall, light_ix_assign_t &light_ix_assign)
{
	float const window_vspace(get_window_vspace()), wall_thickness(get_wall_thickness()), wall_hthick(0.5*wall_thickness), fc_thick(get_fc_thickness());
	float const clearance(get_min_front_clearance_inc_people()), trim_thickness(get_trim_thickness()), trim_height(get_trim_height());
	float const r90_center_bias(1.5*get_doorway_width() + wall_thickness);
	bool const can_have_r90(0.45*room.get_sz_dim(!dim) - r90_center_bias > 1.5*window_vspace); // if enough space for kitchen appliances
	rgen.rand_mix();
	// {open with dining, open with no dining, closed with counter, long with dining and no kitchen}; must be closed with counter if no doorway
	int const style(no_doorway ? 2 : (rgen.rand() % (can_have_r90 ? 4 : 3)));
	bool const has_dining(style == 0 || style == 3), is_open(style != 2), has_kitchen(1/*style != 3*/), is_r90(style == 3);
	float const front_wall(room.d[dim][dir]), bot_wall_z1(room.z1() + fc_thick), bot_wall_z2(zval + 0.35*window_vspace);
	colorRGBA const &wall_color(interior->mall_info->mall_wall_color), &trim_color(get_trim_color());
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t front_area(room), windows_area;
	set_wall_width(front_area, front_wall, wall_thickness, dim);

	for (cube_t const &w : interior->mall_info->storefronts) {
		if (w.intersects(front_area)) {windows_area.assign_or_union_with_cube(w);}
	}
	bool const has_windows(!windows_area.is_all_zeros());

	if (!is_open && has_windows) { // windows_area should never be zero area
		// add counter rather than glass window and door; must remove windows and doorways but keep storefronts
		remove_if_intersects(interior->int_windows, front_area);
		remove_if_intersects(interior->mall_info->store_doorways, front_area);
		// add front top and bottom wall sections and counter
		cube_t wall(windows_area);
		set_wall_width(wall, room.d[dim][dir], wall_hthick, dim);
		set_cube_zvals(wall, bot_wall_z1, bot_wall_z2);
		add_restaurant_counter(wall, dim, dir, room_id, light_amt, 1, 1, rgen); // leave_end_gaps=1 (doesn't span the entire room), add_cash_registers=1
		set_cube_zvals(wall, windows_area.z2()-wall_thickness, windows_area.z2()); // narrow strip to fill the bottom edge of the top wall
		objs.emplace_back(wall, TYPE_STAIR_WALL, room_id, dim, 0, RO_FLAG_HANGING, light_amt, SHAPE_CUBE, wall_color); // upper wall; draw bottom
		// add trim around the opening
		float const trim_hwidth(0.5*trim_height);
		cube_t trim(wall);
		trim.expand_in_dim(!dim, trim_hwidth);
		trim.expand_in_dim( dim, 2.0*trim_thickness);
		set_wall_width(trim, wall.z1(), trim_hwidth, 2);
		unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_UNTEXTURED); // not reflective
		objs.emplace_back(trim, TYPE_METAL_BAR, room_id, dim, 0, flags, light_amt, SHAPE_CUBE, trim_color, 0); // draw all sides
		set_cube_zvals(trim, bot_wall_z2+0.1*trim_hwidth, trim.z1());

		for (unsigned d = 0; d < 2; ++d) { // left, right
			set_wall_width(trim, windows_area.d[!dim][d], trim_hwidth, !dim);
			objs.emplace_back(trim, TYPE_METAL_BAR, room_id, dim, 0, flags, light_amt, SHAPE_CUBE, trim_color, EF_Z2); // skip top
		}
	}
	// split between public space in the front and commercial kitchen in the back
	bool const sdim(dim ^ is_r90), pdir(is_r90 ? rgen.rand_bool() : dir); // split dim and public dir
	float const pdsign(pdir ? -1.0 : 1.0), rlen(room.get_sz_dim(sdim));
	float fb_split(0.0); // front/public vs. back/private split pos
	// Note: for R90 restaurants, we ideally want the public area connected to the front door and the private area connected to the back door,
	// but we can't move the doors at this point, and the wall/counter must be axis aligned, so make both doors in the public space;
	// min shift is set to 1.5x doorway width + wall_thickness to allow for counter in the front and the back door to open
	if (is_r90) {fb_split = room.d[sdim][pdir] + pdsign*(0.5 + rgen.rand_uniform(0.0, 0.05))*rlen + pdsign*r90_center_bias;} // not too close to door
	else        {fb_split = front_wall + pdsign*((has_dining ? 0.5 : 0.25) + rgen.rand_uniform(0.0, 0.1))*rlen;}
	// add separator wall on top and bottom, and metal counter
	bool const leave_end_gaps(1);
	float const upper_wall_z1(((is_r90 && has_windows) ? windows_area.z2() : (zval + get_floor_ceil_gap())) + 0.5*trim_height);
	cube_t wall(room);
	set_cube_zvals(wall, bot_wall_z1, bot_wall_z2);
	set_wall_width(wall, fb_split, wall_hthick, sdim);
	cube_t upper_wall(wall);
	if (leave_end_gaps) {wall.expand_in_dim(!sdim, -1.25*clearance);}
	set_cube_zvals(upper_wall, upper_wall_z1, (room.z2() - fc_thick));
	objs.emplace_back(upper_wall, TYPE_STAIR_WALL, room_id, sdim, 0, 0, light_amt, SHAPE_CUBE, wall_color); // draw sides only
	div_wall = upper_wall; // store for use with placing ceiling lights, which happens later
	// add upper wall bottom trim
	cube_t trim(upper_wall);
	set_cube_zvals(trim, upper_wall.z1()-trim_height, upper_wall.z1());
	trim.expand_in_dim(sdim, trim_thickness);
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_UNTEXTURED); // not reflective
	interior->room_geom->objs.emplace_back(trim, TYPE_METAL_BAR, room_id, sdim, 0, flags, light_amt, SHAPE_CUBE, trim_color, get_skip_mask_for_xy(!sdim)); // skip ends
	unsigned const objs_start(objs.size());
	cube_t const counter(add_restaurant_counter(wall, sdim, pdir, room_id, light_amt, leave_end_gaps, is_open, rgen)); // add_cash_registers=is_open

	if (is_open && has_windows) { // add a blocker for the windows and door
		cube_t blocker(windows_area);
		blocker.expand_in_dim(dim, clearance);
		objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, 0, 0, light_amt, SHAPE_CUBE);
	}
	if (leave_end_gaps) { // add end gap blockers
		for (unsigned d = 0; d < 2; ++d) {
			cube_t blocker(room);
			blocker.d[!sdim][!d] = wall.d[!sdim][d];
			set_wall_width(blocker, fb_split, 1.5*clearance, sdim);
			objs.emplace_back(blocker, TYPE_BLOCKER, room_id, 0, 0, 0, 1.0);
		}
	}
	// populate kitchen and public areas
	room_t pub_area(room), kitchen(room);
	pub_area.d[sdim][!pdir] = counter.d[sdim][ pdir];
	kitchen .d[sdim][ pdir] = counter.d[sdim][!pdir];

	if (has_kitchen) {
		add_commercial_kitchen_objs(rgen, kitchen, zval, room_id, 0, light_amt, objs_start, objs_start, light_ix_assign); // floor_ix=0, no lights_start
	}
	else {
		// something more like a coffee shop with a row of smaller appliances along the back wall?
	}
	if (is_open) { // add vending machines
		cube_t const pub_room_bounds(get_walkable_room_bounds(pub_area));
		unsigned const vend_types[2] = {VEND_DRINK, VEND_SNACK}; // one of each type
		for (unsigned n = 0; n < 2; ++n) {add_vending_machine_type(rgen, room, zval, room_id, light_amt, objs_start, pub_room_bounds, vend_types[n]);}
	}
	if (has_dining) {add_cafeteria_objs(rgen, pub_area, zval, room_id, 0, light_amt, objs_start);} // floor_ix=0
	else {add_corner_trashcans(rgen, pub_area, zval, room_id, light_amt, objs_start, !sdim, 1);} // both_ends=1
}

cube_t building_t::add_restaurant_counter(cube_t const &wall, bool dim, bool dir, unsigned room_id, float light_amt, bool leave_end_gaps, bool add_cash_registers, rand_gen_t &rgen) {
	float const window_vspace(get_window_vspace()), wall_thickness(get_wall_thickness()), clearance(get_min_front_clearance_inc_people());
	vect_room_object_t &objs(interior->room_geom->objs);
	add_short_wall_with_trim(wall, dim, room_id, light_amt, interior->mall_info->mall_wall_color); // bottom wall
	cube_t counter(wall);
	set_cube_zvals(counter, wall.z2(), (wall.z2() + wall_thickness));
	counter.expand_in_dim(dim, 2.5*wall_thickness);
	unsigned const skip_faces(leave_end_gaps ? 0 : get_skip_mask_for_xy(!dim)); // skip ends if no gaps
	objs.emplace_back(counter, TYPE_METAL_BAR, room_id, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_CUBE, LT_GRAY, skip_faces);
	// add a collider + blocker around the counter and the windows/entrance
	cube_t blocker(counter);
	blocker.z1() = wall.z1();
	objs.emplace_back(blocker, TYPE_COLLIDER, room_id, dim, RO_FLAG_INVIS,  0, light_amt, SHAPE_CUBE); // for player and people
	blocker.expand_in_dim(dim, clearance);
	objs.emplace_back(blocker, TYPE_BLOCKER,  room_id, dim, RO_FLAG_NOCOLL, 0, light_amt, SHAPE_CUBE); // for objects
	cube_t place_area(counter);
	place_area.expand_in_dim(!dim, -wall_thickness); // skip ends, which may be inside the wall trim
	float const counter_len(place_area.get_sz_dim(!dim));
	vect_cube_t avoid;
	
	if (add_cash_registers && building_obj_model_loader.is_model_valid(OBJ_MODEL_CASHREG)) {
		// add cash registers
		unsigned const num_cr(2 + (rgen.rand()%3)); // 2-4
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_CASHREG)); // D, W, H
		float const height(0.12*window_vspace), hwidth(0.5*height*sz.y/sz.z), hdepth(0.5*height*sz.x/sz.z), edge_space(1.5*hwidth);

		if (hwidth < 0.25*counter_len) { // counter is wide enough for cash registers; should be true
			cube_t cr;
			set_cube_zvals(cr, place_area.z2(), place_area.z2()+height);
			set_wall_width(cr, place_area.get_center_dim(dim), hdepth, dim); // centered on the counter

			for (unsigned n = 0; n < num_cr; ++n) {
				for (unsigned N = 0; N < 10; ++N) {
					set_wall_width(cr, rgen.rand_uniform(place_area.d[!dim][0]+edge_space, place_area.d[!dim][1]-edge_space), hwidth, !dim);
					cube_t cr_exp(cr);
					cr_exp.expand_in_dim(!dim, edge_space); // require a gap between cash registers
					if (has_bcube_int(cr_exp, avoid)) continue; // blocked
					objs.emplace_back(cr, TYPE_CASHREG, room_id, dim, dir, RO_FLAG_IN_MALL, light_amt, SHAPE_CUBE);
					avoid.push_back(cr);
					break; // success
				} // for N
			} // for n
		}
	}
	// place other items between cash registers
	unsigned const num_cont  (rgen.rand() % 6); // 0-5
	unsigned const num_pizza (rgen.rand() % 4); // 0-3
	unsigned const num_drinks(add_cash_registers ? (rgen.rand() % 9) : 0); // 0-8
	unsigned const num_fruit (add_cash_registers ? (rgen.rand() % 7) : 0); // 0-6
	unsigned const num_apple_bowls(rgen.rand() % 3); // 0-2

	for (unsigned n = 0; n < num_cont; ++n) { // containers
		switch (rgen.rand()%3) { // TYPE_SILVER?
		case 0: // cup
			if (place_cup_on_obj  (rgen, place_area, room_id, light_amt, avoid)) {avoid.push_back(objs.back());}
			break;
		case 1: // plate or bowl
			if (place_plate_on_obj(rgen, place_area, room_id, light_amt, avoid, (rgen.rand_float() < 0.35))) {avoid.push_back(objs.back());}
			break;
		case 2: { // tray
			add_cafeteria_tray_to_surface(place_area, dim, room_id, light_amt, avoid, rgen);
			break;
		}
		} // end switch
	}
	for (unsigned n = 0; n < num_pizza; ++n) { // pizza; should these be stacked?
		if (place_pizza_on_obj(rgen, place_area, room_id, light_amt, avoid)) {
			objs.back().flags |= RO_FLAG_IN_MALL;
			objs.back().dim    = dim; // set dim, but leave dir random
			avoid.push_back(objs.back());
		}
	}
	for (unsigned n = 0; n < num_apple_bowls; ++n) { // bowls of apples
		unsigned const bowl_obj_ix(objs.size());
		if (place_bowl_of_apples_on_obj(rgen, place_area, room_id, light_amt, avoid)) {avoid.push_back(objs[bowl_obj_ix]);}
	}
	for (unsigned n = 0; n < num_drinks; ++n) { // drinks; no wine; should they be grouped together?
		if (rgen.rand_bool() ? place_bottle_on_obj(rgen, place_area, room_id, light_amt, avoid, 0, BOTTLE_TYPE_BEER) :
			                   place_dcan_on_obj  (rgen, place_area, room_id, light_amt, avoid, 0, DRINK_CAN_TYPE_BEER)) {avoid.push_back(objs.back());}
	}
	for (unsigned n = 0; n < num_fruit; ++n) { // fruit; apples in a bowl?
		if (rgen.rand_bool() ? place_banana_on_obj(rgen, place_area, room_id, light_amt, avoid) :
			                   place_apple_on_obj (rgen, place_area, room_id, light_amt, avoid)) {avoid.push_back(objs.back());}
	}
	return counter;
}

bool building_t::add_object_to_tray(cube_t const &tray, bool dim, unsigned room_id, float light_amt, bool no_alcohol,
	vect_cube_t const &avoid, rand_gen_t &rgen, unsigned &prev_type)
{
	cube_t place_area(tray);
	place_area.z2() = tray.z1() + 0.1*tray.dz(); // top of tray
	place_area.expand_by_xy(-0.1*tray.get_sz_dim(dim));
	unsigned type_ix(rgen.rand() % 6);
	if (type_ix == prev_type) {type_ix = ((type_ix+1) % 6);} // use a different type
	prev_type = type_ix;

	switch (type_ix) { // simplified version of the logic in add_mall_table_with_chairs()
	case 0: return place_bottle_on_obj(rgen, place_area, room_id, light_amt, avoid, 0, (no_alcohol ? BOTTLE_TYPE_COKE    : BOTTLE_TYPE_WINE   ));
	case 1: return place_dcan_on_obj  (rgen, place_area, room_id, light_amt, avoid, 0, (no_alcohol ? DRINK_CAN_TYPE_COKE : DRINK_CAN_TYPE_BEER));
	case 2: return place_cup_on_obj   (rgen, place_area, room_id, light_amt, avoid);
	case 3: return place_plate_on_obj (rgen, place_area, room_id, light_amt, avoid, (rgen.rand_float() < 0.35)); // plate or (less often) bowl
	case 4: return place_banana_on_obj(rgen, place_area, room_id, light_amt, avoid);
	case 5: return place_apple_on_obj (rgen, place_area, room_id, light_amt, avoid);
	}
	return 0; // never gets here
}
unsigned building_t::add_objects_to_tray(cube_t const &tray, bool dim, unsigned room_id, float light_amt, bool no_alcohol, rand_gen_t &rgen, unsigned max_num) {
	unsigned const num(rgen.rand() % (max_num+1));
	static vect_cube_t avoid;
	unsigned prev_type(255); // start at an invalid type
	avoid.clear();

	for (unsigned n = 0; n < num; ++n) {
		if (add_object_to_tray(tray, dim, room_id, light_amt, no_alcohol, avoid, rgen, prev_type)) {avoid.push_back(interior->room_geom->objs.back());}
	}
	return avoid.size();
}
void building_t::add_cafeteria_tray(cube_t const &tray, bool dim, unsigned room_id, float light_amt, bool no_alcohol, rand_gen_t &rgen) {
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const tray_obj_ix(objs.size());
	objs.emplace_back(tray, TYPE_FOOD_TRAY, room_id, dim, 0, RO_FLAG_NOCOLL, light_amt, SHAPE_ROUNDED_CUBE, LT_GRAY);
	unsigned const num_objs(add_objects_to_tray(tray, dim, room_id, light_amt, no_alcohol, rgen, 2)); // 0-2 objects
	if (num_objs > 0) {objs[tray_obj_ix].flags |= RO_FLAG_ADJ_TOP;} // flag tray as having something on it

	if (rgen.rand_float() < 0.4) { // add a stain on the tray; this won't be removed when the tray is taken, but hopefully that's okay
		float const stain_radius(rgen.rand_uniform(0.25, 0.45)*min(tray.dx(), tray.dy()));
		colorRGBA const color(get_stain_color(rgen, 1)); // is_food=1
		point const pos(tray.xc(), tray.yc(), (tray.z1() + 0.2*tray.dz())); // centered, on top of top surface of tray, under any placed items
		interior->room_geom->decal_manager.add_blood_or_stain(pos, stain_radius, color, 0, 2, 1); // is_blood=0; +z
		objs[tray_obj_ix].flags |= RO_FLAG_BROKEN; // mark as stained
	}
}
void building_t::add_cafeteria_tray_to_surface(cube_t const &surface, bool dim, unsigned room_id, float light_amt, vect_cube_t &avoid, rand_gen_t &rgen) {
	float const width(min(0.22f*get_window_vspace(), 0.8f*surface.get_sz_dim(!dim))), depth(min(0.63*width, 0.95*surface.get_sz_dim(dim))), height(0.03*width);
	cube_t tray;
	set_cube_zvals(tray, surface.z2(), surface.z2()+height);
	vector3d const tray_sz(0.5*(dim ? width : depth), 0.5*(dim ? depth : width), height);

	for (unsigned n = 0; n < 4; ++n) { // 4 attempts
		gen_xy_pos_for_cube_obj(tray, surface, tray_sz, height, rgen);
		if (has_bcube_int(tray, avoid)) continue; // blocked
		add_cafeteria_tray(tray, dim, room_id, light_amt, 0, rgen); // no_alcohol=0
		avoid.push_back(tray);
		break;
	}
}

bool building_interior_t::obj_on_restaurant_counter(room_object_t const &obj) const {
	if (!obj.in_mall()) return 0;
	int const store_id(get_store_id_for_room(obj.room_id));
	if (store_id < 0) return 0;
	return (mall_info->stores[store_id].store_type == STORE_FOOD);
}

