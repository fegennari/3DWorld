// 3D World - Parking Garages/Structures
// by Frank Gennari 5/27/25

#include "function_registry.h"
#include "buildings.h"
//#include "city_model.h"


bool building_t::add_parking_structure_entrance(rand_gen_t rgen) {
	assert(interior && !interior->rooms.empty());
	cube_t const &part(parts.front()); // sets the exterior space
	room_t const &room(interior->rooms.front()); // main above ground room is first; sets the interior space
	float const entrance_width(2.3*get_parked_car_size().y), extend_len(1.0*get_parked_car_size().x), door_width(get_doorway_width());
	bool const wdim(rgen.rand_bool()), wdir(rgen.rand_bool()), wside(rgen.rand_bool()); // choose a random wall + end to try first
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
				set_cube_zvals(driveway, ground_floor_z1, entrance.z1());
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
				else                   {face_mask = (1 << (!dim)) + sf[!dim][side];} // inside only
			}
			else {
				if (exterior_surfaces) {face_mask = 3 + sf[dim][!dir] + sf[!dim][!side];} // outside and exposed end
				else                   {face_mask = (1 << dim) + sf[dim][dir];} // inside only
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

		for (unsigned n = 0; n < 4; ++n) { // exterior: top and bottom; interior: inside faces
			if (f < 2) { // cut out slots for doors on the ground floor lower wall and second floor upper wall
				for (unsigned lu = 0; lu < 2; ++lu) { // split into lower and upper sections
					if (f == (lu ? num_floors : 0)) continue; // no lower/upper section for this floor
					cube_t wall(sides[n]);
					if (lu) {wall.z1() = zval + int_ext_wall_fc_gap;} // upper wall
					else    {wall.z2() = zval - int_ext_wall_fc_gap;} // lower wall
					wall_parts.clear();
					subtract_cubes_from_cube(wall, door_cuts, wall_parts, temp, 2); // check zval overlap
					unsigned const face_mask(exterior_surfaces ? 4 : face_masks[n]);
					for (cube_t const &w : wall_parts) {walls.emplace_back(w, face_mask);}
				} // for lu
			}
			else {walls.emplace_back(sides[n], (exterior_surfaces ? 4 : face_masks[n]));} // add entire wall
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

void building_t::add_parking_struct_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix,
	unsigned num_floors, unsigned &nlights_x, unsigned &nlights_y, float &light_delta_z, light_ix_assign_t &light_ix_assign)
{
	add_parking_garage_objs(rgen, room, zval, room_id, floor_ix, num_floors, nlights_x, nlights_y, light_delta_z, light_ix_assign);
	//cube_t const &part(parts.front()); // above ground part
	//vect_room_object_t &objs(interior->room_geom->objs);
	
	if (!interior->parking_entrance.is_all_zeros()) {
		// add ticket booth, barrier, etc.
	}
	// TODO: anything else to add?
}
