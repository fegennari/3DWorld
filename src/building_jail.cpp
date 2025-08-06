// 3D World - Jail Cells and Prison Buildings
// by Frank Gennari 7/02/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

float const JAIL_DOOR_WIDTH = 0.8; // relative to regular door width; about the min size that male people/zombies can fit through
float const JAIL_BARS_THICK = 0.4; // relative to wall width

extern object_model_loader_t building_obj_model_loader;

bool has_key_3d_model();


bool building_t::divide_part_into_jail_cells(cube_t const &part, unsigned part_id, unsigned gen_index, rand_gen_t &rgen, bool try_short_dim) {
	vector2d const part_sz(part.get_size_xy());
	float const dx(part_sz.x), dy(part_sz.y);
	bool const long_dim(dx < dy), hall_dim(long_dim ^ try_short_dim), in_basement(part.z1() < ground_floor_z1);
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), wall_hthick(0.5*wall_thickness), fc_thick(get_fc_thickness());
	float const door_width(get_doorway_width()), door_hwidth(0.5*door_width), room_len(part_sz[hall_dim]), room_width(part_sz[!hall_dim]);
	float min_hall_width(2.1*door_width), min_cell_depth(max(floor_spacing, 2.1f*door_width));
	int const num_end_windows(get_num_windows_on_side(part, !hall_dim));
	float const window_width(room_width/num_end_windows);
	if (num_end_windows >= 3) {max_eq(min_cell_depth, window_width);} // if at least three windows, make sure cells contains a full window
	if (gen_index > 2) {min_hall_width += 0.1*floor_spacing*(gen_index - 2);} // increase hall width incrementally as we fail to add stairs
	float const hall_centerline(part.get_center_dim(!hall_dim));
	float min_room_width(2*min_cell_depth + min_hall_width), extra_width(room_width - min_room_width);
	unsigned skip_side(2); // 2=neither
	vector<room_t> &rooms(interior->rooms);
	cube_t hall_room(part);
	if (interior->prison_halls.empty()) {interior->prison_halls = parts;} // start as entire part

	if (extra_width < 0.0) { // too narrow to fit jail cells; should be rare
		if (try_short_dim || room_len > 1.8*room_width) { // consider skipping the cells on one side and creating a wider hallway
			min_room_width -= min_cell_depth;
			extra_width     = room_width - min_room_width;
			// if we have enough space, skip side furthest from edge of bcube and keep exterior cells
			if (extra_width > 0.0) {skip_side = (hall_centerline < bcube.get_center_dim(!hall_dim));}
		}
		if (skip_side == 2) { // if we didn't decide to skip a side above
			if (!try_short_dim) {return divide_part_into_jail_cells(part, part_id, gen_index, rgen, 1);} // try the other dim; try_short_dim=1
			add_room(part, part_id); // add a single room for the entire part
			rooms.back().assign_all_to(RTYPE_JAIL); // maybe shouldn't be a jail in this case?
			return 0;
		}
	}
	float cell_depth(min_cell_depth + min(0.5*min_cell_depth, extra_width/3.0)); // distribute extra width across the hall and cells on either side
	min_eq(cell_depth, 1.2f*max(min_cell_depth, window_width)); // not so wide that we clip a second window
	unsigned added_cells(0); // 1 bit for lo wall, 2 bit for hi wall
	vect_cube_with_ix_t cells;

	if (in_basement) { // basement jail cell
		float const min_cell_len(1.3*min_cell_depth), stairs_ext(door_width + wall_thickness);
		unsigned const num_cells(room_len/min_cell_len); // floor
		float const cell_len(room_len/num_cells);
		bool const dim(!hall_dim);

		for (unsigned d = 0; d < 2; ++d) { // each side
			if (d == skip_side) continue;
			cube_t cell(part);
			cell.d[dim][!d] = hall_room.d[dim][d] = part.d[dim][d] + (d ? -1.0 : 1.0)*cell_depth; // room ends at side of central hallway

			for (unsigned n = 0; n < num_cells; ++n) {
				float const low_edge(part.d[!dim][0] + n*cell_len);
				cell.d[!dim][0] = low_edge;
				cell.d[!dim][1] = ((n+1 == num_cells) ? part.d[!dim][1] : (low_edge + cell_len)); // end exactly at the part
				if (interior->is_cube_close_to_doorway(cell, part, door_width, 1)) continue; // check for extended basement door; inc_open=1
				// check for stairs above so that they have space to extend below
				bool blocked_by_stairs(0);

				for (stairwell_t const &s : interior->stairwells) {
					if (s.z1() > part.z2() + 0.5*floor_spacing) continue; // doesn't extend to basement
					cube_t ext_cube(s);
					if (s.z1() < part.z2() - 1.5*floor_spacing) {ext_cube.expand_in_dim(s.dim, stairs_ext);} // multi-floor basement stairs; expand at both ends
					else {ext_cube.d[s.dim][!s.dir] += (s.dir ? -1.0 : 1.0)*stairs_ext;} // add padding on exit side
					if (ext_cube.intersects_xy(cell)) {blocked_by_stairs = 1; break;}
				}
				if (blocked_by_stairs) continue;
				cells.emplace_back(cell, (2*dim + d));
				added_cells |= (d ? 2 : 1);
			} // for n
		} // for d
	}
	else if (is_rotated()) {
		// exterior walls are rotated with the building - can't add jail cells
	}
	else { // above ground jail cell
		// start by getting the windows associated with this part
		cube_t window_area(part);
		window_area.expand_by_xy(0.5*wall_thickness); // include ext walls exactly on the edges
		window_area.expand_in_z(-fc_thick); // exclude walls on stacked parts above or below
		vect_vnctcc_t const &wall_quad_verts(get_all_drawn_window_verts_as_quads());
		float tx1(0.0), tx2(0.0), tz1(0.0), tz2(0.0);
		cube_t c;

		if (!try_short_dim) { // calculate sum of wall length in each dim and choose the other dim if too unbalanced
			float wall_len[2] = {0,0};

			for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
				if (!get_wall_quad_window_area(wall_quad_verts, i, c, tx1, tx2, tz1, tz2)) continue; // no windows in this wall
				if (!c.intersects(window_area)) continue; // wrong part, skip
				c.intersect_with_cube_xy(window_area); // only windows in this part count
				bool const dim(c.dy() < c.dx());
				wall_len[!dim] += c.get_sz_dim(!dim);
			}
			if (wall_len[long_dim] < 0.67*wall_len[!long_dim]) {return divide_part_into_jail_cells(part, part_id, gen_index, rgen, 1);} // try_short_dim=1
		}
		for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
			if (!get_wall_quad_window_area(wall_quad_verts, i, c, tx1, tx2, tz1, tz2)) continue; // no windows in this wall
			if (!c.intersects(window_area)) continue; // wrong part, skip
			bool const dim(c.dy() < c.dx()), dir(wall_quad_verts[i].get_norm()[dim] > 0.0);
			if (          dim == hall_dim ) continue; // only keep windows on long part edges
			if ((unsigned)dir == skip_side) continue;
			assert(c.get_sz_dim(dim) == 0.0); // must be zero size in one dim (X or Y oriented); could also use the vertex normal
			// here we only care about the window width, not the height, because the rooms will span the same floors as the window
			float const window_width(c.get_sz_dim(!dim)/(tx2 - tx1));
			float cell_edge(part.d[dim][dir] + (dir ? -1.0 : 1.0)*cell_depth); // room ends at side of central hallway
			cube_t cell(part);
			// attempt to move to not intersect window
			float const window_hspacing(part_sz[dim]/get_num_windows_on_side(part, dim)); // for hall end windows
			float const shift_cell_edge(shift_val_to_not_intersect_window(part, cell_edge, window_hspacing, get_window_h_border(), dim));
			float const new_cell_depth(fabs(shift_cell_edge - cell.d[dim][dir]));
			float const new_hall_hwidth(shift_cell_edge - hall_centerline + ((skip_side < 2) ? cell_depth : 0.0)); // empty cell row counts as hallway
			if (new_cell_depth > min_cell_depth && 2.0*new_hall_hwidth > min_hall_width) {cell_edge = shift_cell_edge;} // shift if cell and hall are wide enough
			cell.d[dim][!dir] = cell_edge; // front of cell

			for (float xy = tx1; xy < tx2; xy += 1.0) { // windows along each wall
				float const low_edge(c.d[!dim][0] + (xy - tx1)*window_width);
				cell.d[!dim][0] = low_edge;
				cell.d[!dim][1] = low_edge + window_width;
				if (!window_area.contains_cube_xy(cell) || !part.intersects(cell)) continue; // not contained in this part (outside or adjacent to)
				if (interior->is_cube_close_to_doorway(cell, part, door_width, 1)) continue; // inc_open=1
				// at this point there should be no placed stairs, elevators, or ext doors, so no need to check valid placement
				cell.intersect_with_cube(part); // will assert if slightly outside
				cells.emplace_back(cell, (2*dim + dir));
				hall_room.d[dim][dir] = cell_edge;
				added_cells |= (dir ? 2 : 1);
			} // for xy
		} // for i
		if (!added_cells && !try_short_dim) {return divide_part_into_jail_cells(part, part_id, gen_index, rgen, 1);} // try the other dim; try_short_dim=1
	}
	cube_t extra_room_area;

	if (cells.empty()) {extra_room_area = part;} // all extra room area
	else { // at least one cell was added; make this a cell block hallway
		bool const dim(!hall_dim), two_sides(added_cells == 3);
		float const hall_width(1.25*min_hall_width);
		cube_t cell_block(part);

		// split non-cells side into more rooms if hallway is wide enough
		if (!in_basement && room_len > 4.0*floor_spacing && hall_room.get_sz_dim(dim) > (two_sides ? 2 : 1)*hall_width + 2.5*floor_spacing) {
			extra_room_area = part;

			if (two_sides) { // split into two hallways
				// Note: in this case, hall_room spans both hallways and the extra room, and must be clipped to the parent room before use
				for (unsigned dir = 0; dir < 2; ++dir) {
					cell_block = part;
					cell_block.d[dim][!dir] = extra_room_area.d[dim][dir] = hall_room.d[dim][dir] + (dir ? -1.0 : 1.0)*hall_width;
					add_prison_cells(cells, cell_block, part_id, (1 << (2*dim + (1-dir)))); // skip cells on opposite side
				}
				cells.clear(); // done adding cells
			}
			else { // cells on one side
				assert(added_cells == 1 || added_cells == 2);
				bool const dir(added_cells == 2);
				cell_block.d[dim][!dir] = hall_room.d[dim][!dir] = extra_room_area.d[dim][dir] = hall_room.d[dim][dir] + (dir ? -1.0 : 1.0)*hall_width;
			}
		}
		else if (gen_index < 2) { // skip if we fail to add stairs on the first two iterations
			// look for gaps between cells along interior walls and add room(s) there; these will be sub-rooms of the cell block
			vect_cube_t cands;
			subtract_cube_from_cube(cell_block, hall_room, cands);
			subtract_cubes_from_cubes(cells, cands);
			std::sort(cands.begin(), cands.end(), [](cube_t const &a, cube_t const &b) {return (a.get_area_xy() > b.get_area_xy());}); // largest to smallest area
			
			for (cube_t const &cand : cands) { // Note: should be at least the size of a cell
				if (cand.get_sz_dim(hall_dim) < 4.0*door_width) continue; // too narrow
				bool const dir(hall_centerline < cand.get_center_dim(dim)); // side of hallway
				float const hall_other_side(hall_room.d[dim][!dir]), dsign(dir ? 1.0 : -1.0);
				float wall_pos(cand.d[dim][!dir] + dsign*wall_hthick), window_hspacing(0.0); // align with front of cell walls
				bool at_part_edge[2]={};
				for (unsigned d = 0; d < 2; ++d) {at_part_edge[d] = (cand.d[hall_dim][d] == part.d[hall_dim][d]);}

				if (!in_basement && (at_part_edge[0] || at_part_edge[1])) { // attempt to move to not intersect window
					window_hspacing = part_sz[dim]/get_num_windows_on_side(part, dim); // for hall end windows
					wall_pos        = shift_val_to_not_intersect_window(part, wall_pos, window_hspacing, get_window_h_border(), dim);
				}
				float const new_hall_width(dsign*(wall_pos - hall_other_side)), pref_hall_width(1.25*min_hall_width);
				if (new_hall_width < min_hall_width) continue; // hall too narrow
				
				if (new_hall_width > pref_hall_width && gen_index == 0) { // try expanding into hallway, but only on the first iteration as this may block stairs
					float cand_wall_pos(hall_other_side + dsign*pref_hall_width);
					
					if (window_hspacing > 0.0) {
						cand_wall_pos = shift_val_to_not_intersect_window(part, cand_wall_pos, window_hspacing, get_window_h_border(), dim);
						if (dsign*(cand_wall_pos - hall_other_side) < min_hall_width) {cand_wall_pos = wall_pos;} // too narrow, restore original wall_pos
					}
					wall_pos = cand_wall_pos; // maybe shift the wall outward
				}
				float const old_wall_pos(cand.d[dim][!dir]);
				cube_t sub_room(cand);
				sub_room.d[dim][!dir] = wall_pos;
				assert(sub_room.is_strictly_normalized());
				if (sub_room.get_sz_dim(dim) < min_cell_depth) continue; // room too narrow

				if (dsign*(old_wall_pos - wall_pos) > wall_thickness) { // add horizontal walls
					cube_t wall(sub_room);
					clip_wall_to_ceil_floor(wall, fc_thick);
					wall.d[dim][!dir] = wall_pos - dsign*wall_hthick;
					wall.d[dim][ dir] = old_wall_pos;

					for (unsigned d = 0; d < 2; ++d) {
						if (at_part_edge[d]) continue; // no wall needed
						set_wall_width(wall, sub_room.d[hall_dim][d], wall_hthick, hall_dim);
						interior->walls[hall_dim].push_back(wall);
					}
				}
				add_prison_room(sub_room, part_id, 1, 1); // inc_half_walls=1, is_nested=1
				cube_t wall(sub_room);
				clip_wall_to_ceil_floor(wall, fc_thick);
				set_wall_width(wall, wall_pos, wall_hthick, dim);
				// insert a door; gen_interior_int() won't add a door because this wall is contained in the parent room
				float const door_pos(rgen.rand_uniform(sub_room.d[hall_dim][0]+door_width, sub_room.d[hall_dim][1]-door_width));
				insert_door_in_wall_and_add_seg(wall, (door_pos - door_hwidth), (door_pos + door_hwidth), hall_dim, dir, 0, 0, 0, 0, 1); // opens into hallway; jail_door=1
				interior->walls[dim].push_back(wall); // add remainder
				break; // only add the largest area cand that's valid
			} // for cand
		}
		if (!cells.empty()) {add_prison_cells(cells, cell_block, part_id, 0);} // add cells if not added in split case above
		assert(part_id < interior->prison_halls.size());
		interior->prison_halls[part_id] = hall_room;
	}
	if (!extra_room_area.is_all_zeros()) { // add the extra room(s)
		bool const dim(!hall_dim);
		float const sub_room_width(extra_room_area.get_sz_dim(dim));
		float const target_area(25.0*floor_spacing*floor_spacing), target_len(min(8.0f*floor_spacing, max(4.0f*floor_spacing, target_area/sub_room_width)));
		unsigned num_sub_rooms(room_len / target_len);
		float pref_door_pos(0.0);

		if (num_sub_rooms > 1) { // split further if large
			// if large enough, add a hallway connecting the two ends so that we can have private rooms to the sides
			bool const add_hall(room_len > 8.0f*floor_spacing);
			float const hall_width(max(3.0*(door_width + wall_thickness), 0.04*room_len)), room_len_no_hall(room_len - (add_hall ? hall_width : 0.0));
			float const room_step(room_len_no_hall/num_sub_rooms), rand_amt(0.2*room_step), door_wall_space(1.2*door_hwidth);
			float const room_start(extra_room_area.d[hall_dim][0]);
			unsigned hall_room_ix(num_sub_rooms); // start at an invalid index
			
			if (add_hall) {
				++num_sub_rooms; // hall is an extra room
				hall_room_ix = num_sub_rooms/2; // center room
				if (!(num_sub_rooms & 1) && rgen.rand_bool()) {hall_room_ix -= 1;} // choose a random side if even
			}
			cube_t sub_room(extra_room_area);

			for (unsigned n = 0; n < num_sub_rooms; ++n) {
				bool const is_last(n+1 == num_sub_rooms), is_hall(n == hall_room_ix);
				float &lo_edge(sub_room.d[hall_dim][0]), &hi_edge(sub_room.d[hall_dim][1]);
				if (is_last)      {hi_edge = extra_room_area.d[hall_dim][1];} // end at exactly the high edge
				else if (is_hall) {hi_edge = lo_edge + hall_width;}
				else              {hi_edge = room_start + (n+1)*room_step + rand_amt*rgen.signed_rand_float() - ((n > hall_room_ix) ? (room_step - hall_width) : 0.0);}
				assert(sub_room.is_strictly_normalized());
				add_prison_room(sub_room, part_id, 1, 0, is_hall); // inc_half_walls=1, is_nested=0
				if (is_hall) {pref_door_pos = sub_room.get_center_dim(hall_dim);} // place door in center of hallway
				if (is_last) break;
				lo_edge = hi_edge; // shift to next room; next room starts where the last room ends
				// add wall separating rooms
				cube_t wall(sub_room);
				clip_wall_to_ceil_floor(wall, fc_thick);
				set_wall_width(wall, hi_edge, wall_hthick, hall_dim);
				// insert a door in the wall; this is optional since the room connecting code may add a door anway, but that doesn't always happen for short walls
				float door_pos(sub_room.d[dim][0] + rgen.rand_uniform(0.35, 0.65)*sub_room_width); // center 30% of wall
				max_eq(door_pos, wall.d[dim][0]+door_wall_space); // door must be contained in wall
				min_eq(door_pos, wall.d[dim][1]-door_wall_space);
				insert_door_in_wall_and_add_seg(wall, (door_pos - door_hwidth), (door_pos + door_hwidth), dim, rgen.rand_bool()); // random open_dir; not a jail door
				interior->walls[hall_dim].push_back(wall); // add remainder
			} // for n
		}
		else { // one big room
			add_prison_room(extra_room_area, part_id, 1, 0); // inc_half_walls=1, is_nested=0
		}
		// add side walls; doors will be cut into these walls to connect rooms later
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (!(added_cells & (1 << dir))) continue; // no cells/hall on this side
			cube_t wall(extra_room_area);
			clip_wall_to_ceil_floor(wall, fc_thick);
			set_wall_width(wall, extra_room_area.d[dim][dir], wall_hthick, dim);
			
			if (pref_door_pos != 0.0) { // add hall jail door that opens out of the hallway; doesn't cover the part separator wall
				insert_door_in_wall_and_add_seg(wall, (pref_door_pos - door_hwidth), (pref_door_pos + door_hwidth), !dim, dir, 0, 0, 0, 0, 1);
			}
			interior->walls[dim].push_back(wall);
		} // for dir
	}
	return !cells.empty();
}

void building_t::add_prison_room(cube_t const &room, unsigned part_id, bool inc_half_walls, bool is_nested, bool is_hall) {
	unsigned num_lights(1);
	if (is_hall) {num_lights = max(1, round_fp(0.4*max(room.dx(), room.dy())/min(room.dx(), room.dy())));}
	add_room(room, part_id, num_lights, is_hall);
	if (inc_half_walls) {interior->rooms.back().set_office_floorplan();}
	if (is_nested     ) {interior->rooms.back().set_is_nested();}
}

void building_t::add_prison_cells(vect_cube_with_ix_t const &cells, cube_t const &cell_block, unsigned part_id, unsigned skip_ix_mask) {
	float const wall_thickness(get_wall_thickness()), wall_hthick(0.5*wall_thickness), fc_thick(get_fc_thickness());
	vector<room_t> &rooms(interior->rooms);
	unsigned parent_room_id(rooms.size()); // the room added after all cells below
	assert(cell_block.is_strictly_normalized());

	for (cube_with_ix_t const &cell : cells) {
		if (!(skip_ix_mask & (1 << cell.ix))) {++parent_room_id;}
	}
	for (cube_with_ix_t const &cell : cells) {
		if (skip_ix_mask & (1 << cell.ix)) continue;
		bool const dim(cell.ix >> 1), dir(cell.ix & 1);
		unsigned const room_id(rooms.size());
		add_prison_room(cell, part_id, 0, 1); // inc_half_walls=0, is_nested=1
		rooms.back().assign_all_to(RTYPE_JAIL_CELL);
		rooms.back().mark_open_wall(dim, !dir);
		// add side interior walls
		cube_t wall(cell);
		clip_wall_to_ceil_floor(wall, fc_thick);

		for (unsigned d = 0; d < 2; ++d) {
			float const wall_pos(cell.d[!dim][d]);
			if (fabs(wall_pos - parts[part_id].d[!dim][d]) < wall_thickness) continue; // at edge of part - no wall
			set_wall_width(wall, wall_pos, 0.5*wall_thickness, !dim);
			interior->walls[!dim].push_back(wall);
		}
		// add the door
		bool const hinge_side(room_id & 1);
		float const bars_depth_pos(cell.d[dim][!dir] + (dir ? 1.0 : -1.0)*wall_hthick), bars_hthick(JAIL_BARS_THICK*wall_hthick);
		cube_t bars(cell);
		clip_wall_to_ceil_floor(bars, fc_thick);
		set_wall_width(bars, bars_depth_pos, bars_hthick, dim);
		add_jail_cell_door(bars, room_id, parent_room_id, !dim, dir, hinge_side);
	} // for cell
	add_room(cell_block, part_id); // add a single room for the entire cell block (hallway + cells), last, since cells are sub-rooms
	rooms.back().assign_all_to(RTYPE_JAIL);
	rooms.back().set_has_subroom();
}

void building_t::add_prison_jail_cell_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	// find the open wall dim/dir
	unsigned num_open_walls(0);
	bool dim(0), dir(0);

	for (unsigned wdim = 0; wdim < 2; ++wdim) {
		for (unsigned wdir = 0; wdir < 2; ++wdir) {
			if (!room.has_open_wall(wdim, wdir)) continue;
			dim = bool(wdim);
			dir = bool(wdir);
			++num_open_walls;
		}
	}
	dir ^= 1; // bars dir, not wall dir
	assert(num_open_walls == 1);
	bool const in_basement(zval < ground_floor_z1), sink_on_back_wall(in_basement); // basement has no window, so put the sink on the back wall
	bool const hinge_side(room_id & 1), bed_side(!hinge_side); // door opens to a random side consistent per room
	bool const is_lit(rgen.rand_bool());
	float const wall_hthick(0.5*get_wall_thickness());
	colorRGBA const bar_color(detail_color*0.8); // slightly darker than the roof
	cube_t const &part(parts[room.part_id]);
	cube_t cell(room); // cell wall bounds
	set_cube_zvals(cell, zval, zval+get_floor_ceil_gap());
	cube_t place_area(cell); // used for object placement; exclude window frame area
	cell.d[dim][!dir] += (dir ? 1.0 : -1.0)*wall_hthick; // exclude front/bars only
	place_area.expand_in_dim(dim, -wall_hthick); // exclude front and back
	
	for (unsigned d = 0; d < 2; ++d) { // exclude interior side walls
		bool const at_ext_wall(fabs(room.d[!dim][d] - part.d[!dim][d]) < 3.0*wall_hthick); // Note: tolerance is needed for interior wall room shrink
		if (at_ext_wall && zval < ground_floor_z1) continue; // no windows in basement, no shrink needed
		place_area.d[!dim][d] -= (d ? 1.0 : -1.0)*(at_ext_wall ? 0.8 : 1.0)*wall_hthick; // ext wall needs small shrink to avoid intersecting window bars/trim
		if (!at_ext_wall) {cell.d[!dim][d] -= (d ? 1.0 : -1.0)*wall_hthick;} // cell_area goes right up to the exterior wall
	}
	// Note: open wall/bars dim is inverted because the calls below use the hall dim
	float const bars_depth_pos(room.d[dim][!dir] + (dir ? 1.0 : -1.0)*wall_hthick), bars_hthick(JAIL_BARS_THICK*wall_hthick);
	add_jail_cell_bars(place_area, cube_t(), room_id, tot_light_amt, !dim, dir, hinge_side, bar_color);
	populate_jail_cell(rgen, cell, place_area, zval, room_id, tot_light_amt, !dim, dir, bed_side, sink_on_back_wall, is_lit, bars_hthick, bars_depth_pos);
}

struct room_pref_t {
	room_type rtype=0;
	bool allow_mult   =0  ; // allow multiple rooms of this type
	bool allow_st_el  =0  ; // allow stairs or elevator in this room; hard constraint
	bool is_private   =0  ; // private room; can't be along a path between other rooms
	float min_size    =0.0; // in short dim, relative to floor spacing
	float max_size    =0.0; // in long  dim, relative to floor spacing; 0 is infinite
	float base_weight =0.0; // base priority; higher is more likely
	float gfloor_scale=0.0; // weight factor for ground floor rooms; positive is preferred
	float basement_val=0.0; // weight factor if in basement; higher is preferred, negative is avoided
	float window_val  =0.0; // weight factor if room has windows; higher is preferred, negative is avoided
};
room_pref_t const room_prefs[] = {
	// rtype, allow_mult, allow_st_el, is_private, min_size, max_size, base_weight, floor_scale, basement_val, window_val
	{RTYPE_VISIT,     0, 1, 1, 2.5, 0.0, 0.5, 1.0, -1.0,  0.0}, // flagged as private since it's separated into two parts
	{RTYPE_LAUNDRY,   0, 1, 0, 0.0, 4.0, 0.5, 0.5,  1.0, -0.5},
	{RTYPE_OFFICE,    1, 1, 0, 0.0, 3.0, 0.0, 0.0,  0.0,  1.0},
	{RTYPE_CAFETERIA, 0, 0, 0, 4.0, 0.0, 1.0, 0.0,  0.0,  0.0},
	{RTYPE_GYM,       0, 1, 0, 2.0, 7.0, 1.0, 0.0,  0.5,  0.0},
	{RTYPE_SHOWER,    0, 0, 1, 2.0, 6.0, 1.0, 0.5,  1.0, -1.0},
	{RTYPE_CLASS,     1, 0, 0, 3.0, 8.0, 0.5, 0.5,  0.0,  0.5},
	{RTYPE_BATH,      1, 0, 1, 0.0, 3.0, 0.0, 0.0,  0.0, -0.5}
};
bool building_t::assign_and_fill_prison_room(rand_gen_t rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt,
	unsigned objs_start, unsigned lights_start, unsigned floor_ix, bool is_basement, colorRGBA const &chair_color)
{
	assert(interior);
	unsigned const num_room_pref(sizeof(room_prefs)/sizeof(room_pref_t));
	vector2d const room_sz(room.get_size_xy());
	float const floor_spacing(get_window_vspace()), min_sz(room_sz.get_min_val()/floor_spacing), max_sz(room_sz.get_max_val()/floor_spacing);
	bool const ground_floor(floor_ix == 0 && !is_basement); // Note: prisons don't have stacked parts
	bool const has_st_el(room_has_stairs_or_elevator(room, zval, floor_ix));
	bool const has_ext_wall(!is_basement && count_ext_walls_for_room(room, zval) > 0);
	bool const on_walk_path(/*count_num_int_doors(room) > 1 &&*/ is_room_on_critical_path(room_id, zval));
	vector<pair<float, room_type>> cands;
	//cout << TXT(min_sz) << TXT(max_sz) << TXT(is_basement) << TXT(ground_floor) << TXT(has_st_el) << TXT(has_ext_wall) << TXT(on_walk_path) << endl;

	for (unsigned pass = 0; pass < 2; ++pass) { // first pass is more restrictive
		bool const strict(pass == 0);
		cands.clear();

		for (unsigned i = 0; i < num_room_pref; ++i) {
			room_pref_t const &p(room_prefs[i]);
			if (has_st_el    && !p.allow_st_el) continue;
			if (on_walk_path &&  p.is_private ) continue;

			float weight(p.base_weight);
			
			if (!p.allow_mult && interior->has_room_type(p.rtype)) {
				if (strict) continue;
				weight -= 10.0; // penalize more
			}
			if (min_sz < p.min_size) {
				if (strict) continue;
				weight -= 2.0*(p.min_size - min_sz); // penalize based on undersize
			}
			if (p.max_size > 0.0 && max_sz > p.max_size) {
				if (strict) continue;
				weight -= 1.0*(max_sz - p.max_size); // penalize based on oversize
			}
			if (ground_floor) {weight += p.gfloor_scale;}
			if (is_basement ) {weight += p.basement_val;}
			if (has_ext_wall) {weight += p.window_val  ;} // assumes exterior walls have windows
			//cout << "cand " << room_names[p.rtype] << " weight " << weight << " pass " << pass << endl;
			weight += rgen.rand_float(); // add some randomness
			cands.emplace_back(-weight, p.rtype); // negate weight to get min value
		} // for i
		std::sort(cands.begin(), cands.end());
		unsigned added_bathroom_objs_mask(0); // unused

		for (auto const &cand : cands) {
			room_type const rtype(cand.second);
			//cout << "try " << room_names[rtype] << " weight " << -cand.first << endl;

			switch (rtype) { // split on rtype
			case RTYPE_VISIT:
				if (!add_visit_room_objs(rgen, room, zval, room_id, tot_light_amt, objs_start, lights_start)) continue;
				break;
			case RTYPE_LAUNDRY:
				if (!add_laundry_objs(rgen, room, zval, room_id, tot_light_amt, objs_start, added_bathroom_objs_mask)) continue;
				break;
			case RTYPE_OFFICE: {
				vect_cube_t blockers; // unused
				if (!add_office_objs(rgen, room, blockers, chair_color, zval, room_id, floor_ix, tot_light_amt, objs_start, is_basement)) continue;
				break;
			}
			case RTYPE_CAFETERIA:
				if (!add_cafeteria_objs(rgen, room, zval, room_id, floor_ix, tot_light_amt, objs_start)) continue;
				break;
			case RTYPE_GYM:
				if (!add_gym_objs(rgen, room, zval, room_id, tot_light_amt, objs_start)) continue;
				break;
			case RTYPE_SHOWER:
				if (!add_shower_room_objs(rgen, room, zval, room_id, tot_light_amt, objs_start)) continue;
				break;
			case RTYPE_CLASS: {
				unsigned td_orient(0); // unused
				zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_CARPET);
				if (!add_classroom_objs(rgen, room, zval, room_id, floor_ix, tot_light_amt, objs_start, chair_color, td_orient)) continue;
				break;
			}
			case RTYPE_BATH:
				if (!add_bathroom_objs(rgen, room, zval, room_id, tot_light_amt, lights_start, objs_start, floor_ix, is_basement, 0, added_bathroom_objs_mask)) continue;
				break;
			default: assert(0);
			}
			//cout << "select " << room_names[rtype] << endl;
			room.assign_to(rtype, floor_ix);
			interior->set_room_type(rtype);
			return 1; // success
		} // for cand
	} // for pass
	//cout << "no cand" << endl;
	return 0;
}

bool building_t::add_gym_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	bool const rubber_flooring(rgen.rand_bool());
	unsigned const sub_type(rubber_flooring ? 0 : 1); // select light colored wood
	zval = add_flooring(room, zval, room_id, tot_light_amt, (rubber_flooring ? FLOORING_RUBBER : FLOORING_WOOD), sub_type);
	//vect_room_object_t &objs(interior->room_geom->objs);
	// TODO: TYPE_GYM_WEIGHT, TYPE_BENCH, TYPE_LOCKER, clothing, etc.
	add_door_sign("Gym", room, zval, room_id);
	return 1;
}
bool building_t::add_visit_room_objs(rand_gen_t rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned lights_start) {
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors
	if (doorways.size() < 2) return 0; // must have both public and prisoner doors
	float const floor_spacing(get_window_vspace()), door_width(get_doorway_width());
	float const wall_thick(get_wall_thickness()), wall_hthick(0.5*wall_thick), trim_thick(get_trim_thickness());
	float const min_side_width(1.2*floor_spacing), clearance(0.1*floor_spacing + wall_thick);
	vector2d const room_sz(room.get_size_xy());
	bool dim(room_sz.x < room_sz.y); // long dim
	if (room_sz[!dim] < 2.1*min_side_width + wall_thick) return 0; // too narrow
	vect_room_object_t &objs(interior->room_geom->objs);
	
	// try to find dividing wall
	cube_t room_interior(get_room_wall_bounds(room));
	set_cube_zvals(room_interior, zval, (zval + get_floor_ceil_gap()));
	float wall_pos(0.0);
	bool valid(0);
	
	for (unsigned cand_dim = 0; cand_dim < 2; ++cand_dim) { // consider both dims
		float const lo_pos(room.d[!dim][0] + min_side_width), hi_pos(room.d[!dim][1] - min_side_width);

		for (unsigned n = 0; n < 25; ++n) {
			wall_pos = rgen.rand_uniform(lo_pos, hi_pos);
			cube_t test_cube(room);
			set_wall_width(test_cube, wall_pos, clearance, !dim);
			cube_t test_cube_exp(test_cube);
			test_cube_exp.expand_in_dim(dim, wall_thick); // to intersect doors at room ends
			bool has_door_on_side[2] = {0,0};
			valid = 1;

			for (door_stack_t const &ds : doorways) {
				if (ds.intersects(test_cube_exp)) {valid = 0; break;}
				has_door_on_side[ds.get_center_dim(!dim) < wall_pos] = 1;
			}
			if (!valid) continue;
			cube_t cand_wall(room);
			set_wall_width(cand_wall, wall_pos, wall_hthick, !dim);

			for (auto i = objs.begin()+lights_start; i != objs.end(); ++i) { // check for lights blocking wall
				if (i->intersects(cand_wall)) {valid = 0; break;}
			}
			if (!valid || !has_door_on_side[0] || !has_door_on_side[1])            {valid = 0; continue;} // must have a door on each side
			if (interior->is_blocked_by_stairs_or_elevator(test_cube, door_width)) {valid = 0; continue;}
			if (is_cube_close_to_exterior_doorway(test_cube, door_width, 1))       {valid = 0; continue;}
			break; // valid
		} // for n
		if (valid) break; // success
		if (cand_dim == 1) return 0; // both dims invalid
		if (room_sz[!dim] < 4.0*floor_spacing) return 0; // too short in short dim
		dim ^= 1;
	} // for cand_dim
	zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_CARPET);
	
	// add wall and window separating prisoner from public sides
	float const wall_z2(zval + 0.3*floor_spacing), table_z2(wall_z2 + trim_thick);
	cube_t wall(room_interior);
	wall.z2() = wall_z2;
	set_wall_width(wall, wall_pos, wall_hthick, !dim);

	for (unsigned d = 0; d < 2; ++d) { // add a gap to prevent Z-fighting with exterior walls
		if (classify_room_wall(room, zval, dim, d, 0) == ROOM_WALL_EXT) {wall.d[dim][d] -= (d ? 1.0 : -1.0)*trim_thick;}
	}
	objs.emplace_back(wall, TYPE_STAIR_WALL, room_id, !dim, 0, RO_FLAG_ADJ_TOP, tot_light_amt, SHAPE_CUBE); // draw top
	cube_t window(wall);
	set_cube_zvals(window, wall_z2, room_interior.z2());
	window.expand_in_dim(!dim, -0.3*wall_hthick); // narrower than wall
	window.expand_in_dim( dim,     -wall_hthick); // shorten ends so that trim meets interior room wall
	interior->int_windows.emplace_back(window, room_id);
	room.set_interior_window();
	
	// add tables and dividers; both extend to both sides of the wall/window symmmetrically
	cube_t const walkable_area(get_walkable_room_bounds(room));
	float const seat_width(1.2*door_width), divider_width(0.5*wall_thick), divider_hwidth(0.5*divider_width), room_len(walkable_area.get_sz_dim(dim));
	unsigned const num_seats(room_len/seat_width);
	float const seat_spacing(room_len/num_seats), table_depth(0.15*floor_spacing), divider_depth(0.3*floor_spacing);
	unsigned const NUM_DIV_COLORS = 4;
	colorRGBA const divider_colors[NUM_DIV_COLORS] = {colorRGBA(0.75, 1.0, 0.9, 1.0), colorRGBA(0.7, 0.8, 1.0), WHITE, DK_GRAY}; // blue-green, light blue
	colorRGBA const divider_color(divider_colors[rgen.rand() % NUM_DIV_COLORS]); // random color
	colorRGBA const table_color(LT_GRAY);
	cube_t table(wall), divider(wall);
	set_cube_zvals(table, (table_z2 - 0.4*wall_thick), table_z2);
	divider.z2() += 0.3*floor_spacing; // taller than table
	table  .expand_in_dim(!dim, table_depth  );
	divider.expand_in_dim(!dim, divider_depth);
	table  .d[dim][0] += divider_hwidth; // starts at divider edge
	table  .d[dim][1]  = walkable_area.d[dim][0] + seat_spacing - divider_hwidth; // ends at adjacent divider edge
	divider.d[dim][1]  = walkable_area.d[dim][0] + divider_width;
	colorRGBA const &chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
	vect_cube_t blockers; // unused
	bool has_prev(0);

	for (unsigned n = 0; n < num_seats; ++n) {
		cube_t test_cube(table);
		set_cube_zvals(test_cube, zval, window.z2()); // full room height
		test_cube.expand_in_dim( dim, divider_hwidth); // need space for dividers
		test_cube.expand_in_dim(!dim, (0.5*door_width + divider_depth - table_depth)); // need clearance to the back

		if (is_obj_placement_blocked(test_cube, room, 1)) { // inc_open_doors=1
			divider.translate_dim(dim, seat_spacing);
			has_prev = 0;
		}
		else { // add table + divider(s)
			objs.emplace_back(table, TYPE_METAL_BAR, room_id, !dim, 0, 0, tot_light_amt, SHAPE_CUBE, table_color, get_skip_mask_for_xy(dim)); // skip sides

			if (!has_prev) { // left divider
				cube_t div(divider);
				max_eq(div.d[dim][0], room_interior.d[dim][0]); // don't clip through walls at the ends
				if (div.get_sz_dim(dim) > 0.0) {objs.emplace_back(div, TYPE_STALL, room_id, !dim, 0, 0, tot_light_amt, SHAPE_SHORT, divider_color);}
				has_prev = 1;
			}
			// right divider
			divider.translate_dim(dim, seat_spacing);
			cube_t div(divider);
			min_eq(div.d[dim][1], room_interior.d[dim][1]); // don't clip through walls at the ends
			if (div.get_sz_dim(dim) > 0.0) {objs.emplace_back(divider, TYPE_STALL, room_id, !dim, 0, 0, tot_light_amt, SHAPE_SHORT, divider_color);}
			
			for (unsigned d = 0; d < 2; ++d) { // add chairs and phones on each side
				// chair at table
				point chair_pos(0.0, 0.0, zval);
				chair_pos[!dim] = 0.5*(divider.d[!dim][d] + table.d[!dim][d]);
				chair_pos[ dim] = table.get_center_dim(dim) + 0.12*rgen.signed_rand_float()*seat_spacing; // slightly misaligned
				add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, !dim, !d, tot_light_amt, 0, 0, 0, 0, 1); // office_chair=0, plastic
				// phone on table
				cube_t half_table(table);
				half_table.d[!dim][!d]  = wall.d[!dim][d];
				half_table.d[!dim][ d] += (d ? 1.0 : -1.0)*wall_hthick; // allow phone to hang a bit off the edge of the table
				place_phone_on_obj(rgen, half_table, room_id, tot_light_amt, !dim, d);
			} // for d
		}
		for (unsigned d = 0; d < 2; ++d) { // add chairs against the far walls, even if no table
			point chair_pos(0.0, 0.0, zval);
			chair_pos[!dim] = walkable_area.d[!dim][d] - rgen.rand_uniform(0.22, 0.26)*(d ? 1.0 : -1.0)*floor_spacing; // a bit misaligned
			chair_pos[ dim] = table.get_center_dim(dim) + 0.08*rgen.signed_rand_float()*seat_spacing; // a bit less misaligned
			add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, !dim, !d, tot_light_amt, 0, 0, 0, 1, 1); // office_chair=0, no_push_out=1, plastic
		}
		table.translate_dim(dim, seat_spacing);
	} // for n
	// add chairs for waiting?
	//place_chairs_along_walls(rgen, room, zval, room_id, tot_light_amt, objs_start, chair_color, 1, 3*num_seats/2); // is_plastic=1
	add_door_sign("Visitation", room, zval, room_id);
	return 1;
}
bool building_t::add_shower_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	// TODO: showers, toilets, lockers, benches
	zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_TILE);
	add_door_sign("Shower", room, zval, room_id);
	return 1;
}

void building_t::add_prison_hall_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	cube_t const hall(get_prison_hall_for_room(room));
	vector2d const hall_sz(hall.get_size_xy());
	bool const dim(hall_sz.x < hall_sz.y); // long dim
	float const floor_spacing(get_window_vspace()), door_width(get_doorway_width());
	float const desk_hwidth(0.5*0.9*floor_spacing), desk_hdepth(0.6*desk_hwidth), edge_space(desk_hwidth + 1.5*door_width);
	if (hall_sz[!dim] < 2.0*(desk_hdepth + 1.25*door_width) || hall_sz[dim] < 2.0*edge_space) return; // too small for guard desk
	float const centerline(hall.get_center_dim(!dim));
	bool const dir(centerline - room.get_center_dim(!dim) < 0.1*floor_spacing*rgen.signed_rand_float()); // facing somewhat toward cells, if they're on one side
	cube_t desk;
	set_cube_zvals(desk, zval, zval+0.32*floor_spacing);
	set_wall_width(desk, centerline, desk_hdepth, !dim);

	for (unsigned n = 0; n < 20; ++n) { // 20 tries
		float const val((n == 0) ? room.get_center_dim(dim) : rgen.rand_uniform((room.d[dim][0] + edge_space), (room.d[dim][1] - edge_space)));
		set_wall_width(desk, val, desk_hwidth, dim);
		cube_t desk_exp(desk);
		desk_exp.expand_by_xy(door_width);
		if (intersects_nested_room  (desk_exp, room_id))     continue;
		if (is_cube_close_to_doorway(desk, room, 0.0, 1, 1)) continue; // too close to a doorway
		if (add_reception_desk(rgen, desk, !dim, dir, room_id, tot_light_amt)) break; // add keys on the desk?
	}
}

bool building_t::add_basement_jail_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
	bool is_lit, colorRGBA const &light_color, light_ix_assign_t &light_ix_assign)
{
	float const floor_spacing(get_window_vspace()), dx(room.dx()), dy(room.dy());
	if (min(dx, dy) < 2.4*floor_spacing || max(dx, dy) < 3.0*floor_spacing) return 0; // too small
	bool const dim(dx < dy); // long dim
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // don't need to handle individual doors here
	if (doorways.empty()) return 0; // error?
	float const room_center(room.get_center_dim(dim));
	cube_t end_doors_span;
	bool door_at_end[2] = {0,0};

	for (door_stack_t const &ds : doorways) {
		if (ds.dim != dim) return 0; // not handling doors in long sides of the room yet
		end_doors_span.assign_or_union_with_cube(ds.get_true_bcube());
		door_at_end[room_center < ds.get_center_dim(dim)] = 1; // assumes one door per room end
	}
	assert(!end_doors_span.is_all_zeros()); // no doors for this room?
	if (is_house) {zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_CONCRETE);} // add concrete over the carpet (even if we don't make it a jail)
	float const door_width(get_doorway_width()), wall_hthick(0.5*get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room));
	set_cube_zvals(room_bounds, zval, (zval + get_floor_ceil_gap()));
	float const room_len(room_bounds.get_sz_dim(dim)), min_cell_len(1.25*floor_spacing);
	unsigned const num_cells(room_len/min_cell_len);
	assert(num_cells > 0);
	float const cell_len(room_len/num_cells); // includes walls between cells
	float const bar_lum(rgen.rand_uniform(0.1, 0.5));
	colorRGBA const bar_color(bar_lum, bar_lum, bar_lum);
	vect_room_object_t &objs(interior->room_geom->objs);
	end_doors_span.expand_in_dim(!dim, 0.35*door_width); // add side padding for door frame, etc.
	unsigned const lock_color_ix(rgen.rand_bool() ? 7 : 1); // black or brown, since they look good on walls and cell doors
	bool added_cell(0), added_lock(0);

	for (unsigned dir = 0; dir < 2; ++dir) { // for each side of the room
		float const wall_pos(end_doors_span.d[!dim][dir]);
		cube_t cell_area(room_bounds);
		cell_area.d[!dim][!dir] = wall_pos; // clip off the hallway
		float const cell_depth(cell_area.get_sz_dim(!dim));
		if (cell_depth < 0.9*floor_spacing) continue; // too narrow
		float const bars_depth_pos(wall_pos + (dir ? 1.0 : -1.0)*wall_hthick), bars_hthick(JAIL_BARS_THICK*wall_hthick);
		bool const sink_on_back_wall(cell_depth < 0.67*cell_len); // wide and shallow cell
		cube_t cell(cell_area);

		for (unsigned n = 0; n < num_cells; ++n) {
			float const lo_edge(room_bounds.d[dim][0] + n*cell_len), div_wall_hwidth(sink_on_back_wall ? bars_hthick : wall_hthick);
			cell.d[dim][0] = lo_edge;
			cell.d[dim][1] = lo_edge + cell_len;
			if (n   > 0        ) {cell.d[dim][0] += div_wall_hwidth;} // reserve space for walls/bars
			if (n+1 < num_cells) {cell.d[dim][1] -= div_wall_hwidth;}
			// add bars and door
			unsigned const door_ix(interior->doors.size());
			float const cell_center(cell.get_center_dim(dim));
			bool const hinge_side((room_center < cell_center) ^ bool(dir) ^ 1), bed_side(!hinge_side); // door opens toward hallway center
			cube_t bars(cell);
			set_wall_width(bars, bars_depth_pos, bars_hthick, !dim);
			add_jail_cell_door(bars, room_id, room_id, dim, dir, hinge_side); // parent_room_id=room_id
			add_jail_cell_bars(cell, interior->doors.back(), room_id, tot_light_amt, dim, dir, hinge_side, bar_color);

			if (rgen.rand_bool()) {
				add_padlock_to_door(door_ix, (1 << lock_color_ix), rgen); // force lock_color_ix
				added_lock = 1;
			}
			if (n > 0) { // add divider wall if not the end cell
				cube_t wall(cell);
				set_wall_width(wall, lo_edge, div_wall_hwidth, dim);

				if (sink_on_back_wall) { // add bars between the cells
					wall.d[!dim][!dir] = bars.d[!dim][!dir]; // flush with the bars
					objs.emplace_back(wall, TYPE_JAIL_BARS, room_id, dim, 0, 0, tot_light_amt, SHAPE_CUBE, bar_color, room_id); // dir=0; use room_id as item_flags for material
				}
				else { // sink on side wall; must place a wall to hold the pipes
					unsigned const flags(is_house ? 0 : RO_FLAG_BACKROOM); // flag as backroom for concrete texture in office buildings
					objs.emplace_back(wall, TYPE_PG_WALL, room_id, dim, 0, flags, tot_light_amt, SHAPE_CUBE, WHITE); // dir=0
				}
			}
			// add a small light in each cell
			cube_t light(cube_top_center(cell));
			light.z1() -= 0.01*floor_spacing;
			light.expand_by_xy(0.06*floor_spacing);
			objs.emplace_back(light, TYPE_LIGHT, room_id, dim, 0, (RO_FLAG_NOCOLL | (is_lit ? RO_FLAG_LIT : 0)), 0.0, SHAPE_CYLIN, light_color); // dir=0 (unused)
			objs.back().obj_id = light_ix_assign.get_next_ix();
			populate_jail_cell(rgen, cell, cell, zval, room_id, tot_light_amt, dim, dir, bed_side, sink_on_back_wall, is_lit, bars_hthick, bars_depth_pos);
		} // for n
		added_cell = 1;
	} // for dir
	if (!added_cell) return 0; // not a jail

	if (added_lock && has_key_3d_model()) { // a door lock was added; add a key hanging on the wall opposite the door if there's a single door
		for (unsigned d = 0; d < 2; ++d) {
			if (door_at_end[d]) continue; // blocked by the door; there should be at least one end door
			float const key_sz(0.018*floor_spacing), xlate((d ? -1.0 : 1.0)*0.4*key_sz);
			point key_pos;
			key_pos.z = zval + 0.55*floor_spacing;
			key_pos[ dim] = room_bounds.d[dim][d];
			key_pos[!dim] = end_doors_span.get_center_dim(!dim);
			cube_t key(key_pos);
			key.expand_by(key_sz*vector3d(0.7, 0.7, 2.0)); // make it square in XY since it's small, to avoid all of the orient logic, but make it larger in Z
			key.translate_dim(dim, xlate); // move inside the wall
			objs.emplace_back(key, TYPE_KEY, room_id, dim, d, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), tot_light_amt, SHAPE_CUBE, lock_colors[lock_color_ix]);
			objs.back().obj_id = lock_color_ix;
			// add nail to place the key on
			float const nail_radius(0.14*key_sz);
			key_pos.z += (dim ? 1.96 : 0.06)*key_sz; // offset correctly based on dim, since the swap of dims used in drawing doesn't rotate about the key hole
			cube_t nail(key_pos);
			nail.expand_in_dim(2,    nail_radius);
			nail.expand_in_dim(!dim, nail_radius);
			nail.d[dim][!d] += 2.5*xlate;
			objs.emplace_back(nail, TYPE_METAL_BAR, room_id, dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, DK_GRAY);
		} // for d
	}
	interior->room_geom->jails.push_back(room); // needed for door open logic
	interior->has_jail = 1;
	return 1;
}

unsigned building_t::get_nested_room_parent(unsigned room_id) const {
	room_t const &room(get_room(room_id));

	if (room.is_nested()) { // nested room (prison cell), find the parent; basement jails share the same room
		for (unsigned i = 0; i < interior->rooms.size(); ++i) {
			room_t const &r(get_room(i));
			if (r.has_subroom() && r.contains_cube(room)) return i;
		}
	}
	return room_id; // no parent; return ourself
}

void building_t::add_jail_cell_door(cube_t const &bars, unsigned room_id, unsigned parent_room_id, bool dim, bool dir, bool hinge_side) {
	float const door_width(get_doorway_width());
	float const door_center(bars.get_center_dim(dim) - (hinge_side ? -1.0 : 1.0)*0.1*door_width); // slightly away from bed and room door
	door_t door(bars, !dim, !dir, 0, 0, hinge_side); // open=0, on_stairs=0
	door.for_jail     = 1;
	door.conn_room[0] = parent_room_id; // may be the same room
	door.conn_room[1] = room_id;
	set_wall_width(door, door_center, 0.5*JAIL_DOOR_WIDTH*door_width, dim);
	door.d[!dim][0] = door.d[!dim][1] = bars.get_center_dim(!dim); // shrink to zero width
	add_interior_door(door, 0, 1, 1); // is_bathroom=0, make_unlocked=1, make_closed=1
}
void building_t::add_jail_cell_bars(cube_t const &cell, cube_t door_bc, unsigned room_id, float tot_light_amt, bool dim, bool dir, bool hinge_side, colorRGBA const &color) {
	float const bars_hwidth(JAIL_BARS_THICK*0.5*get_wall_thickness());
	unsigned const parent_room_id(get_nested_room_parent(room_id)); // sets texture/material

	if (door_bc.is_all_zeros()) { // if not specified, find the door associated with this cell and get its bounds
		for (door_stack_t const &ds : interior->door_stacks) {
			if (ds.is_connected_to_room(room_id)) {door_bc = ds; break;}
		}
		assert(!door_bc.is_all_zeros()); // must be found
	}
	for (unsigned d = 0; d < 2; ++d) { // add bars on both sides of the door; dir is facing outside the cell
		cube_t bars(cell);
		set_wall_width(bars, door_bc.get_center_dim(!dim), bars_hwidth, !dim);
		bars.d[dim][!d] = door_bc.d[dim][d];
		assert(bars.is_strictly_normalized());
		interior->room_geom->objs.emplace_back(bars, TYPE_JAIL_BARS, room_id, !dim, dir, 0, tot_light_amt, SHAPE_CUBE, color, parent_room_id);
	}
}

void building_t::populate_jail_cell(rand_gen_t &rgen, cube_t const &cell, cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt,
	bool dim, bool dir, bool bed_side, bool sink_on_back_wall, bool is_lit, float bars_hthick, float bars_depth_pos)
{
	float const floor_spacing(get_window_vspace()), wall_thick(get_wall_thickness()), dsign(dir ? 1.0 : -1.0), bss(bed_side ? 1.0 : -1.0);
	vect_room_object_t &objs(interior->room_geom->objs);
	// add bed; use place_area so that it doesn't clip a window frame
	cube_t bed(place_area);
	bed.expand_by_xy(-get_trim_thickness()); // add a bit of space around the bed
	bed.z2() = zval + 0.32*floor_spacing; // set height
	float const gap_len(bed.get_sz_dim(!dim)), bed_to_bars_gap(max((gap_len - 0.9*floor_spacing), (bars_hthick + 0.5*wall_thick)));
	bed.d[!dim][!dir] = bars_depth_pos + dsign*bed_to_bars_gap; // set length
	bed.d[ dim][!bed_side] = bed.d[dim][bed_side] - bss*0.5*bed.get_sz_dim(!dim); // set width to half length
	assert(bed.is_strictly_normalized());
	objs.emplace_back(bed, TYPE_BED, room_id, !dim, dir, RO_FLAG_IN_JAIL, tot_light_amt);

	if (rgen.rand_bool()) { // make this a bunk bed; set per-room/row?
		room_object_t &bed(objs.back());
		room_object_t bed2(bed);
		bed2.translate_dim(2, bed.dz());
		bed .flags |= RO_FLAG_ADJ_BOT; // flag as bottom bunk
		bed2.flags |= RO_FLAG_ADJ_TOP; // flag as top    bunk
		objs.push_back(bed2); // Note: bed reference is invalidated here
		// add a small ladder
		float const bed_edge(bed2.d[dim][!bed_side]);
		cube_t ladder(bed2);
		ladder.z1()  = zval; // down to the floor
		ladder.z2() -= 0.05*bed2.dz(); // lower to just above mattress level
		ladder.d[ dim][ bed_side] = bed_edge;
		ladder.d[ dim][!bed_side] = bed_edge - bss*0.4*wall_thick; // set depth
		ladder.d[!dim][ dir]  = bed2.d[!dim][!dir] + dsign*0.2*bed2.get_sz_dim(!dim);
		ladder.d[!dim][!dir] += dsign*0.25*wall_thick; // move away from footboard
		unsigned const flags(RO_FLAG_IN_FACTORY | RO_FLAG_NOCOLL);
		objs.emplace_back(ladder, TYPE_INT_LADDER, room_id, dim, bed_side, flags, tot_light_amt, SHAPE_CUBE, GRAY); // metal, like factory ladder
	}
	// add toilet on the far wall and sink to the side or next to the toilet
	cube_t ts_space(place_area); // toilet and sink space
	ts_space.d[dim][bed_side] = bed.d[dim][!bed_side]; // the part not occupied by the bed
	float const space_width(ts_space.get_sz_dim(dim));

	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_TOILET)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
		float const height(0.35*floor_spacing), width(height*sz.y/sz.z), length(height*sz.x/sz.z);
		float const center_pos(ts_space.d[dim][bed_side] - (sink_on_back_wall ? 0.7 : 0.5)*bss*space_width); // further from the bed if there's a sink on the wall
		cube_t toilet(ts_space);
		toilet.z2() = zval + height;
		set_wall_width(toilet, center_pos, 0.5*width, dim);
		toilet.d[!dim][!dir] = ts_space.d[!dim][dir] - dsign*length;
		objs.emplace_back(toilet, TYPE_TOILET, room_id, !dim, !dir, 0, tot_light_amt);
		if (zval < ground_floor_z1) {add_bathroom_plumbing(objs.back());} // only for basement toilets; above ground toilets may be near a window
		float const tp_zval(zval + 0.7*height), tp_length(0.18*height);

		if (sink_on_back_wall) { // on the back wall, not on bars
			add_tp_roll(cell, room_id, tot_light_amt, !dim, dir, tp_length, tp_zval, (toilet.d[dim][!bed_side] - bss*0.4*width));
		}
		else { // on the side wall next to the toilet
			add_tp_roll(cell, room_id, tot_light_amt, dim, !bed_side, tp_length, tp_zval, toilet.get_center_dim(!dim));
		}
	}
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SINK)); // D, W, H
		float const height(0.45*floor_spacing), width(height*sz.y/sz.z), depth(height*sz.x/sz.z);
		cube_t sink(ts_space);
		sink.z2() = zval + height;

		if (sink_on_back_wall) {
			float const center_pos(ts_space.d[dim][bed_side] - 0.3*bss*space_width); // closer to the bed
			set_wall_width(sink, center_pos, 0.5*width, dim);
			sink.d[!dim][!dir] = ts_space.d[!dim][dir] - dsign*depth;
			objs.emplace_back(sink, TYPE_SINK, room_id, !dim, !dir, 0, tot_light_amt);
		}
		else { // sink on side wall
			set_wall_width(sink, ts_space.get_center_dim(!dim), 0.5*width, !dim);
			sink.d[dim][bed_side] = ts_space.d[dim][!bed_side] + bss*depth;
			objs.emplace_back(sink, TYPE_SINK, room_id, dim, bed_side, 0, tot_light_amt);
		}
		bool maybe_against_window(0);

		if (zval >= ground_floor_z1) { // windows are only above ground
			room_t const &room(get_room(room_id));
			cube_t const &part(parts[room.part_id]);
			bool const sdim(sink_on_back_wall ? !dim : dim), sdir(sink_on_back_wall ? dir : !bed_side);
			maybe_against_window = (fabs(ts_space.d[sdim][sdir] - part.d[sdim][sdir]) < wall_thick && classify_room_wall(room, zval, sdim, sdir, 0) == ROOM_WALL_EXT);
		}
		if (!maybe_against_window) {add_bathroom_plumbing(objs.back());}
	}
}

void building_t::add_window_bars(cube_t const &window, bool dim, bool dir, unsigned room_id) {
	if (!has_int_windows()) return; // no interior (or exterior) drawn windows
	cube_t bars(window);
	bars.expand_in_dim(dim, -0.1*window.get_sz_dim(dim)); // small shrink
	interior->room_geom->objs.emplace_back(bars, TYPE_JAIL_BARS, room_id, dim, dir, 0, 1.0, SHAPE_CUBE, GRAY);
}

bool building_t::stairs_or_elevator_blocked_by_nested_room(cube_t const &c, unsigned room_id) const {
	if (!is_prison()) return 0; // currently, only prisons have jail cells that block stairs and elevators
	cube_t c_exp(c);
	c_exp.expand_by_xy(get_doorway_width()); // leave enough room to the sides for players and AI to walk around
	return intersects_nested_room(c_exp, room_id);
}
bool building_t::intersects_nested_room(cube_t const &c, unsigned room_id, cube_t *blocker) const {
	room_t const &parent(get_room(room_id));
	if (!parent.has_subroom()) return 0; // no nested rooms

	for (int r = (int)room_id-1; r >= 0; --r) { // nested rooms are added before the parent room, so iterate backwards
		room_t const &room(get_room(r));
		if (!room.is_nested() || !parent.contains_cube(room)) break; // no more nested rooms for this parent
		if (!c.intersects(room)) continue;
		if (blocker) {*blocker = room;}
		return 1;
	}
	return 0;
}
void building_t::move_cube_to_not_intersect_sub_room(cube_t &c, cube_t const &place_area, unsigned room_id, bool dim) const {
	cube_t blocker;
	if (!intersects_nested_room(c, room_id, &blocker)) return;
	// blocked; move to closer side that has enough space in short dim
	float const r1(place_area.d[dim][0]), r2(place_area.d[dim][1]), b1(blocker.d[dim][0]), b2(blocker.d[dim][1]);
	float const lo_gap(b1 - r1), hi_gap(r2 - b2), lo_pos(0.5*(b1 + r1)), hi_pos(0.5*(b2 + r2));
	float const cur_pos(c.get_center_dim(dim)), light_width(c.get_sz_dim(dim));
	float new_pos(cur_pos);

	if (lo_gap > light_width) { // lo shift is valid
		if (hi_gap > light_width) { // both shifts are valid
			new_pos = ((fabs(cur_pos - lo_pos) > fabs(cur_pos - hi_pos)) ? hi_pos : lo_pos); // choose pos with smaller delta
		} else {new_pos = lo_pos;} // only lo shift is valid
	} else if (hi_gap > light_width) {new_pos = hi_pos;} // only hi shift is valid
	c.translate_dim(dim, (new_pos - cur_pos));
}

bool building_t::is_prison_door_valid(cube_t const &cand, bool dim, bool &open_dir) const { // check if too close to jail cells
	if (!is_prison()) return 1;
	assert(interior);
	float const doorway_width(get_doorway_width());
	cube_t cand_exp(cand);
	cand_exp.expand_in_dim(dim, doorway_width); // clear a path for the door

	for (unsigned d = 0; d < 2; ++d) { // try both open_dir values
		cube_t open_bc(cand);
		open_bc.expand_in_dim(!dim, doorway_width); // increase width on both sides, since we don't know/can't control which direction the door will open
		open_bc.d[dim][ open_dir] += (open_dir ? 1.0 : -1.0)*doorway_width; // extend outward
		open_bc.d[dim][!open_dir] += (open_dir ? 1.0 : -1.0)*get_wall_thickness(); // shrink inward to prevent invalid intersections
		bool valid(1);

		for (room_t const &r : interior->rooms) {
			if ( r.is_ext_basement())                break; // end at extended basement rooms
			if (!r.is_nested() || !r.is_jail_cell()) continue; // not a jail cell
			if ( r.intersects(cand_exp))             return 0; // invalid for either dim
			if ( r.intersects(open_bc )) {valid = 0; break;} // invalid if hits a cell
		}
		if (valid) { // test that door is contained in a part
			valid = 0;

			for (cube_t const &part : parts) {
				if (part.contains_cube(open_bc)) {valid = 1; break;}
			}
		}
		if (valid) return 1;
		open_dir ^= 1;
	} // for d
	return 0; // neither dir is valid
}

cube_t building_t::get_prison_hall_for_room(room_t const &r) const {
	assert(interior);
	assert(r.part_id < interior->prison_halls.size());
	cube_t hall(interior->prison_halls[r.part_id]);
	hall.intersect_with_cube(r); // needed to handle pair of halls in a single part with two different rooms
	return hall;
}
void building_t::get_prison_cell_block_cubes(unsigned room_id, vect_cube_t &out, bool inc_hallway, bool inc_non_cell_subrooms) const {
	assert(is_prison());
	room_t const &parent(get_room(room_id));
	cube_t const hall(get_prison_hall_for_room(parent));
	assert(parent.intersects(hall)); // parent usually contains hall, but may extend slightly outside
	float const wall_thick(get_wall_thickness());
	out.clear();
	if (inc_hallway) {out.push_back(parent);}
	else {subtract_cube_from_cube(cube_t(parent), hall, out);}
	if (!parent.has_subroom()) return; // no cells

	for (int r = (int)room_id-1; r >= 0; --r) { // Note: includes non-cell rooms as well
		room_t const &room(get_room(r));
		if (!room.is_nested() || !parent.contains_cube(room)) break; // no more nested rooms for this parent
		if (inc_non_cell_subrooms && !room.is_jail_cell()) continue;
		cube_t sub(room);
		sub.expand_by_xy(wall_thick); // make sure to remove the walls to the sides of the cells
		subtract_cube_from_cubes(sub, out); // no holes
	}
}

bool building_t::place_stairs_in_prison_room(cube_t &stairs, unsigned room_id, bool stairs_dim, bool &wall_dir, rand_gen_t &rgen) const {
	vect_cube_t cands;
	get_prison_cell_block_cubes(room_id, cands, 0, 1); // inc_hallway=0, inc_non_cell_subrooms=1 (allow stairs in those subrooms)
	if (cands.empty()) return 0; // room is full
	if (cands.size() > 1) {std::shuffle(cands.begin(), cands.end(), rand_gen_wrap_t(rgen));}
	float const doorway_width(get_nominal_doorway_width()), wall_thickness(get_wall_thickness());
	stairs.expand_in_dim(stairs_dim, -0.1*stairs.get_sz_dim(stairs_dim)); // make it smaller and more likely to fit
	vector2d cut_sz(stairs.get_size_xy()), min_sz(cut_sz);
	min_sz[stairs_dim] += 2.0*doorway_width + 2.0*get_window_vspace(); // must have space at the ends to enter and exit and space for a door to an adj part
	room_t const &room(get_room(room_id));

	for (cube_t const &cand : cands) {
		vector2d const cand_sz(cand.get_size_xy());
		if (cand_sz.x <= min_sz.x || cand_sz.y <= min_sz.y) continue; // too small
		bool const edge_dir(room.get_center_dim(!stairs_dim) < cand.get_center_dim(!stairs_dim)), end_dir(rgen.rand_bool());
		float const side_shift_val(cand.d[!stairs_dim][edge_dir] - stairs.d[!stairs_dim][edge_dir] - (edge_dir ? 1.0 : -1.0)*wall_thickness);
		stairs.translate_dim(!stairs_dim, side_shift_val); // shift to the edge of the room
		float const end_shift_val(cand.d[stairs_dim][end_dir] - stairs.d[stairs_dim][end_dir] - (end_dir ? 1.0 : -1.0)*doorway_width);
		stairs.translate_dim(stairs_dim, end_shift_val); // shift to one end to maximize space for a door
		if (is_cube_close_to_doorway(stairs, room)) continue;
		stairs.intersect_with_cube_xy(room); // shouldn't be needed, except for FP error
		wall_dir = edge_dir;
		return 1;
	} // for cand
	return 0; // failed
}

