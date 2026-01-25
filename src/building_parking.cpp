// 3D World - Parking Garages/Structures
// by Frank Gennari 5/27/25

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for car_t


void set_car_model_color(car_t &car, unsigned btype);
void subtract_cube_from_floor_ceil(cube_t const &c, vect_cube_t &fs);


float building_t::get_parking_road_width() const {
	return 2.3f*get_parked_car_size().y;
}
float building_t::get_parking_ramp_width() const {
	return (is_parking() ? 3.0f : 2.3f)*get_parked_car_size().y; // wider for parking structure so two cars can pass
}

bool building_t::add_parking_structure_entrance(rand_gen_t rgen) {
	assert(interior && !interior->rooms.empty());
	cube_t const &part(parts.front()); // sets the exterior space
	room_t const &room(interior->rooms.front()); // main above ground room is first; sets the interior space
	float const entrance_width(get_parking_ramp_width()), extend_len(1.0*get_parked_car_size().x), door_width(get_doorway_width());
	// choose a random wall + end to try first; if in the side, use the closest street side
	bool const wdim(street_dir ? ((street_dir-1) >> 1) : rgen.rand_bool()), wdir(street_dir ? ((street_dir-1)&1) : rgen.rand_bool()), wside(rgen.rand_bool());
	cube_t entrance;
	entrance.z1() = room    .z1() + get_fc_thickness();
	entrance.z2() = entrance.z1() + get_floor_ceil_gap();

	for (unsigned i = 0; i < 2; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			for (unsigned k = 0; k < 2; ++k) {
				bool const dim(wdim ^ bool(i)), dir(wdir ^ bool(j)), side(wside ^ bool(k));
				float const ent_wall(part.d[dim][dir]), side_wall(room.d[!dim][side]);
				entrance.d[ dim][ dir ] = ent_wall; // flush with exterior wall
				entrance.d[ dim][!dir ] = ent_wall  + (dir  ? -1.0 : 1.0)*extend_len; // extend into building
				entrance.d[!dim][ side] = side_wall; // flush with interior of exterior wall
				entrance.d[!dim][!side] = side_wall + (side ? -1.0 : 1.0)*entrance_width;
				assert(entrance.is_strictly_normalized());
				if (interior->is_blocked_by_stairs_or_elevator(entrance)) continue;
				cube_t test_cube(entrance);
				test_cube.expand_by_xy(door_width);
				bool int_door(0);

				for (auto const &door : doors) { // check exterior doors
					if (test_cube.intersects(door.get_bcube())) {int_door = 1; break;}
				}
				if (int_door) continue;
				interior->parking_entrance = cube_with_ix_t(entrance, (2*dim + dir));
				// add a driveway
				driveway = entrance;
				set_cube_zvals(driveway, ground_floor_z1, (ground_floor_z1 + get_trim_thickness()));
				driveway.d[dim][!dir]  = ent_wall;
				driveway.d[dim][ dir] += (dir ? 1.0 : -1.0)*1.5*get_window_vspace();
				driveway.expand_in_dim(!dim, 0.5*get_park_struct_wall_thick());
				return 1; // success
			} // for k
		} // for j
	} // for i
	return 0; // failed at all 8 locations
}

bool building_t::add_parking_structure_bathroom(rand_gen_t rgen) {
	assert(interior && !interior->rooms.empty());
	room_t &room(interior->rooms.front()); // main above ground room is first
	bool xside(rgen.rand_bool()), yside(rgen.rand_bool()), door_dim(rgen.rand_bool()); // perferred starting corner/side
	float const floor_spacing(get_window_vspace()), door_width(get_doorway_width()), wall_thick(get_wall_thickness()), bathroom_sz(3.2*door_width);

	for (unsigned xv = 0; xv < 2; ++xv, xside ^= 1) {
		for (unsigned yv = 0; yv < 2; ++yv, yside ^= 1) {
			point const corner(room.d[0][xside], room.d[1][yside], room.z1());
			cube_t bathroom(corner);
			bathroom.z2() += floor_spacing;
			bathroom.d[0][!xside] += (xside ? -1.0 : 1.0)*bathroom_sz;
			bathroom.d[1][!yside] += (yside ? -1.0 : 1.0)*bathroom_sz;

			for (unsigned dd = 0; dd < 2; ++dd, door_dim ^= 1) { // try door on both sides
				bool const door_dir(!(door_dim ? yside : xside));
				cube_t br_with_door(bathroom);
				br_with_door.d[door_dim][door_dir] += (door_dir ? 1.0 : -1.0)*door_width; // add space for door to open
				if (interior->is_blocked_by_stairs_or_elevator(br_with_door, 0.0)) continue; // blocked by stairs, elevator, or ramp
				if (!interior->parking_entrance.is_all_zeros() && br_with_door.intersects(interior->parking_entrance)) continue; // too close to entrance
				room_t orig_room(room); // copy as room reference may be invalidated below
				room.copy_from(bathroom); // first (contained) room must be the bathroom
				room.expand_by_xy(-wall_thick); // subtract off the walls
				room.has_stairs = room.has_elevator = 0; // clear flags for bathroom
				calc_room_ext_sides(room);
				room.assign_all_to(RTYPE_BATH);
				room.set_is_nested();
				room.is_office = 1; // hack to treat this as an office bathroom
				orig_room.set_has_subroom();
				interior->rooms.push_back(orig_room); // re-add full room; invalidates room reference
				// update elevator with the new room; required for elevator light room_id
				if (!interior->elevators.empty() && interior->elevators.front().room_id == 0) {++interior->elevators.front().room_id;}
				// add wall sections and door
				float const fc_thick(get_fc_thickness()), door_hwidth(0.5*door_width);

				for (unsigned wdim = 0; wdim < 2; ++wdim) {
					for (unsigned wdir = 0; wdir < 2; ++wdir) {
						cube_t wall(bathroom);
						clip_wall_to_ceil_floor(wall, fc_thick);
						wall.d[wdim][!wdir] = wall.d[wdim][wdir] + (wdir ? -1.0 : 1.0)*wall_thick; // set wall width

						if (bool(wdim) == door_dim && bool(wdir) == door_dir) { // add door in this wall
							float const door_center(bathroom.get_center_dim(!door_dim));
							// opens into bathroom; keep_high_side=0, is_bathroom=0 (not always closed), make_unlocked=1, make_closed=0
							insert_door_in_wall_and_add_seg(wall, (door_center - door_hwidth), (door_center + door_hwidth), !wdim, !door_dir, 0, 0, 1, 0);
						}
						interior->walls[wdim].push_back(wall);
					} // for wdir
				} // for wdim
				interior->ps_bathroom = br_with_door; // includes space for door
				return 1; // done/success
			} // for dd
		} // for yv
	} // for xv
	return 0;
}

void get_door_cuts(vect_tquad_with_ix_t const &doors, float wall_thick, vect_cube_t &door_cuts) {
	for (tquad_with_ix_t const &d : doors) { // find all doors on the ground floor
		if (d.type == tquad_with_ix_t::TYPE_RDOOR) continue; // roof access door - skip
		cube_t dbc(d.get_bcube());
		bool const dim(dbc.dy() < dbc.dx());
		dbc.expand_in_dim(dim, 2.0*wall_thick);
		door_cuts.push_back(dbc);
	}
}

// index bits: enable dims: 1=x, 2=y, 4=z | disable cube faces: 8=x1, 16=x2, 32=y1, 64=y2, 128=z1, 256=z2
void building_t::get_parking_struct_ext_walls(vect_cube_with_ix_t &walls, bool exterior_surfaces) const {
	if (!exterior_surfaces && interior && !interior->parking_str_walls.empty()) { // use cached value
		walls = interior->parking_str_walls;
		return;
	}
	assert(is_parking());
	assert(real_num_parts == (1 + has_basement()));
	cube_t const &part(parts.front()); // above ground part
	float const floor_spacing(get_window_vspace()), wall_thick(get_park_struct_wall_thick());
	float const lower_wall_height(0.35*floor_spacing), upper_wall_height(0.15*floor_spacing), int_ext_wall_fc_gap(0.0);
	float const gap_z1(part.z1() + lower_wall_height), gap_z2(part.z1() + floor_spacing - upper_wall_height);
	unsigned num_floors(calc_num_floors(part, floor_spacing, get_floor_thickness()));
	assert(num_floors > 0);
	cube_t inner(part);
	inner.expand_by_xy(-wall_thick);
	vect_cube_t wall_parts, door_cuts, temp;
	bool has_entrance(0);
	get_door_cuts(doors, wall_thick, door_cuts);

	if (interior && !interior->parking_entrance.is_all_zeros()) { // handle entrance opening
		door_cuts.push_back(interior->parking_entrance);
		has_entrance = 1;
	}
	// add extra vertical wall segments to either side of the door and entrance cut
	for (cube_t const &c : door_cuts) {
		bool const dim((c == interior->parking_entrance) ? (interior->parking_entrance.ix >> 1) : (c.dy() < c.dx()));
		bool const dir(part.get_center_dim(dim) < c.get_center_dim(dim));
		cube_t wall(part);
		set_cube_zvals(wall, gap_z1, gap_z2);
		wall.d[dim][!dir] = wall.d[dim][dir] + (dir ? -1.0 : 1.0)*wall_thick; // set wall thickness

		for (unsigned side = 0; side < 2; ++side) {
			unsigned const sf[2][2] = {{8, 16}, {32, 64}}; // {dim x dir}
			unsigned face_mask(0);
			float const door_edge(c.d[!dim][side]);
			cube_t w(wall);
			w.d[!dim][!side] = door_edge;
			w.d[!dim][ side] = door_edge + (side ? 1.0 : -1.0)*0.75*floor_spacing;
			w.intersect_with_cube(part); // don't extend outside the part

			if (w.get_sz_dim(!dim) < 1.1*wall_thick) { // entrance corner wall - expand in the other dim
				w.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.50*floor_spacing;
				w.intersect_with_cube(part); // clamp again, just in case
				if (exterior_surfaces) {face_mask = 3 + sf[!dim][!side];} // outside and both ends
				else                   {face_mask = (1 << unsigned(!dim)) + sf[!dim][side];} // inside only
			}
			else {
				if (exterior_surfaces) {face_mask = 3 + sf[dim][!dir] + sf[!dim][!side];} // outside and exposed end
				else                   {face_mask = (1 << unsigned(dim)) + sf[dim][dir];} // inside only
			}
			walls.emplace_back(w, face_mask);
		} // for side
	} // for c
	// one extra floor; each wall slice goes from the ceiling of the level below to the top of the wall of the level above, except for the two ends
	for (unsigned f = 0; f <= num_floors; ++f) {
		float const zval(part.z1() + f*floor_spacing), wall_bot(zval - upper_wall_height), wall_top(zval + lower_wall_height);
		cube_t slice(part);
		set_cube_zvals(slice, max(part.z1(), wall_bot), min(wall_top, part.z2()));
		
		if (exterior_surfaces) { // XY exterior walls around entire perimeter; don't need to handle doors here
			unsigned face_mask(3); // all XY to start
			
			if (has_entrance && f < 2) { // clip out entrance on lower two floors
				cube_with_ix_t const &entrance(interior->parking_entrance);
				bool const dim(entrance.ix >> 1);
				unsigned const entrance_face(8 << entrance.ix), front_face_mask(123 - entrance_face);
				face_mask += entrance_face; // disable entrance face
				
				for (unsigned d = 0; d < 2; ++d) {
					cube_t side(slice);
					side.d[!dim][!d] = entrance.d[!dim][d];
					walls.emplace_back(side, front_face_mask); // draw only entrance face
				}
				if (slice.z2() > entrance.z2()) { // draw top part above entrance
					cube_t upper(slice);
					upper.z1() = entrance.z2();
					for (unsigned d = 0; d < 2; ++d) {upper.d[!dim][d] = entrance.d[!dim][d];}
					walls.emplace_back(upper, front_face_mask); // draw only entrance face
				}
			}
			walls.emplace_back(slice, face_mask);
		}
		cube_t sides[4]; // {-y, +y, center -x, center +x}
		subtract_cube_xy(slice, inner, sides);
		unsigned const face_masks[4] = {127-64, 127-32, 127-16, 127-8}; // enable XYZ but skip all XY but {+y, -y, +x, -x} in XY
		bool split_lu(f < 2); // handle doors on the ground floor lower wall and second floor upper wall

		if (!split_lu && have_walkway_ext_door) { // include walkway doors if any overlap in Z
			for (cube_t const &c : door_cuts) {split_lu |= (c.z1() < slice.z2() && c.z2() > slice.z1());}
		}
		for (unsigned n = 0; n < 4; ++n) { // exterior: top and bottom; interior: inside faces
			unsigned face_mask(exterior_surfaces ? 4 : face_masks[n]);

			if (split_lu) { // cut out slots for doors
				for (unsigned lu = 0; lu < 2; ++lu) { // split into lower and upper sections
					if (f == (lu ? num_floors : 0)) continue; // no lower/upper section for this floor
					cube_t wall(sides[n]);
					if (lu) {wall.z1() = zval + int_ext_wall_fc_gap;} // upper wall
					else    {wall.z2() = zval - int_ext_wall_fc_gap;} // lower wall
					wall_parts.clear();
					subtract_cubes_from_cube(wall, door_cuts, wall_parts, temp, 2); // check zval overlap
					if (exterior_surfaces && f == 0) {face_mask |= 128;} // zkip Z1 for bottom floor
					for (cube_t const &w : wall_parts) {walls.emplace_back(w, face_mask);}
				} // for lu
			}
			else { // add entire wall
				if (exterior_surfaces && f == num_floors) {face_mask |= 256;} // skip Z2 for top floor since it will cause Z-fighting with the roof
				walls.emplace_back(sides[n], face_mask);
			}
		} // for n
	} // for f
	if (!exterior_surfaces && interior) {interior->parking_str_walls = walls;} // cache for future use
}

void building_t::get_parking_garage_wall_openings(vect_cube_with_ix_t &openings) const {
	assert(is_parking());
	cube_t const &part(parts.front()); // above ground part
	float const floor_spacing(get_window_vspace()), wall_thick(get_park_struct_wall_thick());
	float const lower_wall_height(0.35*floor_spacing), upper_wall_height(0.15*floor_spacing);
	unsigned num_floors(calc_num_floors(part, floor_spacing, get_floor_thickness()));
	assert(num_floors > 0);
	cube_t inner(part);
	inner.expand_by_xy(-wall_thick);
	vect_cube_t wall_parts, door_cuts, temp;
	get_door_cuts(doors, wall_thick, door_cuts);

	if (interior && !interior->parking_entrance.is_all_zeros()) { // handle entrance opening
		cube_with_ix_t entrance(interior->parking_entrance);
		bool const dim(entrance.ix >> 1), dir(entrance.ix & 1);
		entrance.d[dim][!dir] = inner.d[dim][dir];
		openings.push_back(entrance);
		entrance.d[dim][!dir] += (dim ? -1.0 : 1.0)*0.5*floor_spacing; // extend back
		door_cuts.push_back(entrance);
	}
	for (cube_t &c : door_cuts) {c.expand_in_dim(!(c.dy() < c.dx()), 0.75*floor_spacing);} // widen

	for (unsigned f = 0; f < num_floors; ++f) {
		float const zval(part.z1() + f*floor_spacing);
		cube_t slice(part);
		set_cube_zvals(slice, (zval + lower_wall_height), (zval + floor_spacing - upper_wall_height));
		cube_t sides[4]; // {-y, +y, center -x, center +x}
		unsigned const ixs[4] = {2, 3, 0, 1};
		subtract_cube_xy(slice, inner, sides);

		for (unsigned n = 0; n < 4; ++n) { // exterior: top and bottom; interior: inside faces
			if (f == 0) { // cut out slots for doors on the ground floor lower wall and second floor upper wall
				cube_t wall(sides[n]);
				wall_parts.clear();
				subtract_cubes_from_cube(wall, door_cuts, wall_parts, temp, 2); // check zval overlap
				for (cube_t const &w : wall_parts) {openings.emplace_back(w, ixs[n]);}
			}
			else {openings.emplace_back(sides[n], ixs[n]);} // add entire wall
		} // for n
	} // for f
}

cube_t building_t::get_parking_structure_roof() const {
	assert(!parts.empty());
	cube_t roof(parts[0]); // roof is on the first part
	roof.expand_by_xy(-get_park_struct_wall_thick());
	return roof;
}
void building_t::get_parking_rooftop_avoid(vect_cube_t &avoid) const {
	float const floor_spacing(get_window_vspace()), doorway_width(get_doorway_width());
	cube_t roof(get_parking_structure_roof());
	set_cube_zvals(roof, roof.z2()-get_floor_thickness(), roof.z2()+floor_spacing); // assume one floor tall

	// currently there are no parking spaces or on the roof, so we only need to avoid the ramp, stairs, and elevator
	if (!interior->pg_ramp.is_all_zeros()) {
		cube_with_ix_t ramp(interior->pg_ramp);
		bool const dim(ramp.ix >> 1), side(roof.get_center_dim(!dim) < ramp.get_center_dim(!dim));
		// extend to cover most of the width of the roof so that cars can access all parking space rows
		for (unsigned d = 0; d < 2; ++d) {ramp.d[dim][d] = roof.d[dim][d] - (d ? 1.0 : -1.0)*doorway_width;}
		ramp.d[!dim][!side] -= (side ? 1.0 : -1.0)*get_parking_road_width(); // add space alongside the ramp
		avoid.push_back(ramp);
	}
	for (stairwell_t const &s : interior->stairwells) {
		if (s.intersects(roof)) {avoid.push_back(get_stairs_bcube_expanded(s, 2.0*doorway_width, 0.1*doorway_width, doorway_width));}
	}
	for (elevator_t const &e : interior->elevators) {
		if (e.intersects(roof)) {avoid.push_back(e.get_bcube_padded(2.0*doorway_width));}
	}
}

void building_t::add_parking_roof_lights() {
	assert(interior);
	roof_lights.clear(); // in case the interior was regenerated
	vect_cube_t avoid;
	get_parking_rooftop_avoid(avoid);
	cube_t const roof(get_parking_structure_roof());
	float const floor_spacing(get_window_vspace()), pole_radius(0.025*floor_spacing), pole_height(1.5*floor_spacing), light_radius(6.0*floor_spacing);
	float const zval(roof.z2()), car_len(get_nom_car_size().x);
	bool const rdim(roof.dx() < roof.dy()); // long dim
	cube_t place_area(roof);
	place_area.expand_by_xy(-3.2*pole_radius);
	unsigned num[2]={};
	float space [2]={};

	for (unsigned d = 0; d < 2; ++d) {
		float const sz(place_area.get_sz_dim(d));

		if (bool(d) == !rdim && car_len > 0.0) { // align lights to rows of parking spaces if there are cars
			float const space_len(1.1*car_len), row_span(2*space_len + get_parking_road_width()), rows_span(place_area.get_sz_dim(!rdim));
			num[d] = rows_span/row_span + 1; // take the floor and add 1, since lights are at the ends of rows
		}
		else {
			num[d] = max(2U, (unsigned)ceil(sz/light_radius));
		}
		space[d] = sz/(num[d] - (bool(d) != rdim));
	} // for d
	unsigned const nrows(num[!rdim]), ncols(num[rdim]); // rows are longer than columns
	float const col_center(place_area.get_center_dim(!rdim));
	cube_t light;
	set_cube_zvals(light, zval, zval+pole_height);

	for (unsigned r = 0; r < nrows; ++r) {
		float const rpos(place_area.d[!rdim][0] + r*space[!rdim]); // along the edges
		bool const ldir(rpos < col_center); // facing toward the center
		unsigned const ix(2*(!rdim) + ldir);
		set_wall_width(light, rpos, pole_radius, !rdim);

		for (unsigned c = 0; c < ncols; ++c) {
			float const cpos(place_area.d[rdim][0] + (c + 0.5)*space[rdim]); // centered
			set_wall_width(light, cpos, pole_radius, rdim);
			if (has_bcube_int(light, avoid)) continue;
			roof_lights.emplace_back(light, ix);
		}
	} // for r
	max_eq(bcube.z2(), light.z2());
}

void building_t::get_rooftop_cars(vector<car_t> &cars) const {
	if (!interior || !is_parking()) return;
	vector3d const nom_car_sz(get_nom_car_size());
	if (nom_car_sz.x == 0.0) return; // no car models loaded
	// Note: parking space logic is similar to add_parking_garage_objs(), except there are no walls, pillars, handicap spots, or half rows;
	// this logic is run before parking spaces are added, so we can't use them for reference
	cube_t const roof(get_parking_structure_roof());
	float const zval(roof.z2()), floor_spacing(get_window_vspace());
	float const road_width(get_parking_road_width()), space_len(1.1*nom_car_sz.x), space_width(1.4*nom_car_sz.y), row_span(2*space_len + road_width);
	bool const rdim(roof.dx() < roof.dy()); // long dim - rows of cars
	cube_t place_area(roof);

	if (!interior->pg_ramp.is_all_zeros()) {
		cube_with_ix_t const &ramp(interior->pg_ramp);
		bool const dim(ramp.ix >> 1), side(roof.get_center_dim(!dim) < ramp.get_center_dim(!dim));
		place_area.d[!dim][side] = ramp.d[!dim][!side] - (side ? 1.0 : -1.0)*road_width; // add space alongside the ramp
	}
	place_area.expand_in_dim(!rdim, -3.2*0.025*floor_spacing); // add space for light poles
	place_area.expand_in_dim( rdim, -0.25*space_width); // add space for roads to the sides; there are no pillars to avoid on the roof
	float const row_len(place_area.get_sz_dim(rdim)), rows_span(place_area.get_sz_dim(!rdim));
	unsigned const num_rows(rows_span/row_span), num_spaces(row_len/space_width); // take the floor
	if (num_rows == 0 || num_spaces == 0) return; // shouldn't get here, but if we do, there are no cars on the roof
	float const row_spacing(rows_span/num_rows), space_spacing(row_len/num_spaces);
	float const row_road_width(row_spacing - 2*space_len), row_offset(0.5*(row_road_width + space_len));
	vect_cube_t avoid;
	get_parking_rooftop_avoid(avoid);
	rand_gen_t rgen;
	rgen.set_state(mat_ix, interior->floors.size());
	point center(0.0, 0.0, zval);
	car_t car;
	// encode building index in cur_city, starting with CITY_BIX_START and incrementing with each new building
	car.cur_city      = ((cars.empty() || cars.back().cur_city < CITY_BIX_START) ? CITY_BIX_START : cars.back().cur_city+1);
	car.cur_road_type = TYPE_BUILDING;
	car.dim           = !rdim;
	vector3d car_sz(nom_car_sz);
	if (!car.dim) {swap(car_sz.x, car_sz.y);}

	for (unsigned row = 0; row < num_rows; ++row) { // each row is a road with parking spaces to either side
		float const row_center(place_area.d[!rdim][0] + (row + 0.5)*row_spacing);

		for (unsigned d = 0; d < 2; ++d) { // for each side of the row
			for (unsigned space = 0; space < num_spaces; ++space) {
				if (rgen.rand_float() > 0.45) continue; // 45% chance of adding a car
				car.dir        = (bool(d) ^ (rgen.rand_float() < 0.2)); // 20% chance of backing in
				car.cur_seg    = rgen.rand(); // sets the model and color
				center[!rdim]  = row_center + (d ? 1.0 : -1.0)*row_offset;
				center[ rdim]  = place_area.d[rdim][0] + (space + 0.5)*space_spacing;
				center[!rdim] += 0.02*nom_car_sz.x*rgen.signed_rand_float(); // small random misalign front/back
				center[ rdim] += 0.03*nom_car_sz.y*rgen.signed_rand_float(); // small random misalign side
				car.set_bcube(point(center.x, center.y, zval), car_sz);
				set_car_model_color(car, btype);
				if (!has_bcube_int(car.bcube, avoid)) {cars.push_back(car);}
			} // for space
		} // for d
	} // for row
}

void building_t::add_parking_struct_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	unsigned num_floors, unsigned &nlights_x, unsigned &nlights_y, float &light_delta_z, light_ix_assign_t &light_ix_assign)
{
	add_parking_garage_objs(rgen, room, zval, room_id, floor_ix, num_floors, nlights_x, nlights_y, light_delta_z, light_ix_assign);
	
	if (floor_ix == 0 && !interior->parking_entrance.is_all_zeros()) { // ground floor entrance
		vect_room_object_t &objs(interior->room_geom->objs);
		cube_with_ix_t const &entrance(interior->parking_entrance);
		bool const dim(entrance.ix >> 1), dir(entrance.ix & 1);
		float const front_pos(entrance.d[dim][dir]), centerline(entrance.get_center_dim(!dim)), entrance_width(entrance.get_sz_dim(!dim));
		// add a double entrance/exit divider line
		cube_t space(entrance);
		set_cube_zvals(space, zval, zval+get_rug_thickness());
		space.d[!dim][0] = centerline + 0.01*entrance_width; // close to half the entrance width

		for (unsigned d = 0; d < 2; ++d) {
			objs.emplace_back(space, TYPE_PARK_SPACE, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_HI), 1.0, SHAPE_CUBE, wall_color);
			if (d == 0) {space.translate_dim(!dim, -0.04*entrance_width);}
		}
		// add parking gate/ticket machine/barrier for each side
		float const floor_spacing(get_window_vspace()), dsign(dir ? -1.0 : 1.0);
		float const gate_height(0.42*floor_spacing), gate_width(0.16*floor_spacing), gate_depth(0.10*floor_spacing);
		bool const first_is_open(rgen.rand_bool());
		cube_t gate;
		set_cube_zvals(gate, zval, zval+gate_height);
		set_wall_width(gate, centerline, 0.5*gate_depth, !dim);
		gate.d[dim][ dir] = front_pos        + dsign*0.2*floor_spacing;
		gate.d[dim][!dir] = gate.d[dim][dir] + dsign*gate_width;

		for (unsigned d = 0; d < 2; ++d) { // shift back and alternate sides
			bool const gdir(bool(d) ^ dim ^ dir), is_open(first_is_open ^ bool(d)); // one is open, the other is closed
			// arm is in the opposite dir because it's for the gate on the opposite side, so that cars have space to pull up before the bar
			bool const arm_side(gdir ^ dim), swap_sides(1);
			unsigned const item_flags((unsigned)arm_side | ((unsigned)swap_sides << 1));
			objs.emplace_back(gate, TYPE_PARK_GATE, room_id, !dim, gdir, (is_open ? RO_FLAG_OPEN : 0), 1.0, SHAPE_CUBE, colorRGBA(0.9, 0.6, 0.0), item_flags);
			if (d == 0) {gate.translate_dim(dim, 2.0*dsign*gate_width);}
		}
	}
}

void building_t::add_parking_garage_ramp(rand_gen_t &rgen) {
	bool const is_parking_str(is_parking());
	assert(interior && !is_house && (has_parking_garage || is_parking_str));
	cube_with_ix_t &ramp(interior->pg_ramp);
	assert(ramp.is_all_zeros()); // must not have been set
	cube_t room(has_basement() ? get_basement() : parts.front()); // basement parking garage or above ground parking structure
	if (is_parking_str) {room.z2() = parts.front().z2();} // extend from basement to top of parking structure
	bool const dim(room.dx() < room.dy()); // long/primary dim
	// see building_t::add_parking_garage_objs(); make sure there's space for a ramp plus both exit dirs within the building width
	float const room_width(room.get_sz_dim(!dim)), ramp_width(min(0.25f*room_width, get_parking_ramp_width())), wall_space((is_parking_str ? 1.2 : 1.0)*ramp_width);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const z1(room.z1() + fc_thick), z2(room.z2() + fc_thick); // bottom level room floor to first floor floor
	bool const ramp_pref_xdir(rgen.rand_bool()), ramp_pref_ydir(rgen.rand_bool());
	//if (is_in_city && is_parking()) {} // don't face toward the walkway/skyway? but the skyway hasn't been placed yet
	bool added_ramp(0), dir(0);

	for (unsigned pass = 0; pass < 2 && !added_ramp; ++pass) {
		for (unsigned xd = 0; xd < 2 && !added_ramp; ++xd) {
			for (unsigned yd = 0; yd < 2; ++yd) {
				bool const xdir(bool(xd) ^ ramp_pref_xdir), ydir(bool(yd) ^ ramp_pref_ydir);
				float const xsz((dim ? 2.0 : 1.0)*ramp_width), ysz((dim ? 1.0 : 2.0)*ramp_width); // longer in !dim
				unsigned const num_ext(unsigned(room.d[0][xdir] == bcube.d[0][xdir]) + unsigned(room.d[1][ydir] == bcube.d[1][ydir]));
				if (num_ext < 2-pass) continue; // must be on the exterior edge of the building in both dims for pass 0, and one dim for pass 1
				dir = (dim ? xdir : ydir);
				point corner(room.d[0][xdir], room.d[1][ydir], z1);
				corner[!dim] += (dir ? -1.0 : 1.0)*wall_space; // shift away from the wall so that cars have space to turn onto the level floor
				point const c1((corner.x - 0.001*(xdir ? 1.0 : -1.0)*xsz), (corner.y - 0.001*(ydir ? 1.0 : -1.0)*ysz), z1); // slight inward shift to prevent z-fighting
				point const c2((corner.x + (xdir ? -1.0 : 1.0)*xsz), (corner.y + (ydir ? -1.0 : 1.0)*ysz), z2);
				cube_t const ramp_cand(c1, c2);
				assert(ramp_cand.is_strictly_normalized());
				cube_t test_cube(ramp_cand);
				test_cube.expand_in_dim(!dim, ramp_width); // extend outward for clearance to enter/exit the ramp (ramp dim is actually !dim)
				if (interior->is_blocked_by_stairs_or_elevator(test_cube)) continue;
				// check for backrooms door, in case it was placed already (but currently it's not)
				if (has_ext_basement() && get_ext_basement_door_blocker().intersects(test_cube)) continue; // blocked by extended basement door
				ramp = cube_with_ix_t(ramp_cand, (((!dim)<<1) + dir)); // encode dim and dir in ramp index field
				added_ramp = 1;
				break; // done
			} // for yd
		} // for xd
	} // for pass
	if (!added_ramp) return; // what if none of the 4 corners work for a ramp?
	// add landings, which are used to draw the vertical edges of the cutout
	unsigned num_floors(calc_num_floors(room, window_vspacing, floor_thickness));
	float z(room.z1() + window_vspacing); // start at upper floor rather than lower floor

	if (!is_parking()) {
		// FIXME: rooms on the ground floor above ramps aren't yet handled (except for parking structures), so clip ramps to avoid disrupting their floors until this is fixed
		ramp.z2() -= 2.0*floor_thickness;
		--num_floors;
		interior->ignore_ramp_placement = 1; // okay to place room objects over ramps because the floor has not been removed
	}
	for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) { // skip first floor - draw pairs of floors and ceilings
		landing_t landing(ramp, 0, f, !dim, dir, 0, SHAPE_RAMP, 0, (f+1 == num_floors), 0, 1); // for_ramp=1
		set_cube_zvals(landing, (z - fc_thick), (z + fc_thick));
		if (is_parking() && f+1 == num_floors) {landing.z2() -= fc_thick;} // parking structure top floor is on the roof and doesn't include a floor above
		interior->landings.push_back(landing);
	}
	// cut out spaces from floors and ceilings
	subtract_cube_from_floor_ceil(ramp, interior->floors  );
	subtract_cube_from_floor_ceil(ramp, interior->ceilings);
	// make rooms over the ramp of type RTYPE_RAMP_EXIT
}
