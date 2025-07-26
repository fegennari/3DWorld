// 3D World - Jail Cells and Prison Buildings
// by Frank Gennari 7/02/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

float const JAIL_DOOR_WIDTH = 0.8; // relative to regular door width; about the min size that male people/zombies can fit through
float const JAIL_BARS_THICK = 0.4; // relative to wall width

extern object_model_loader_t building_obj_model_loader;

bool has_key_3d_model();


bool building_t::divide_part_into_jail_cells(cube_t const &part, unsigned part_id, rand_gen_t &rgen, bool try_short_dim) {
	// TODO: we need to create a room for the stairs, elevator, common room, etc.; maybe the smallest part, or a subset of a part
	vector2d const part_sz(part.get_size_xy());
	float const dx(part_sz.x), dy(part_sz.y);
	bool const long_dim(dx < dy), hall_dim(long_dim ^ try_short_dim), in_basement(part.z1() < ground_floor_z1);
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), wall_hthick(0.5*wall_thickness), fc_thick(get_fc_thickness());
	float const door_width(get_doorway_width()), room_width(part_sz[!hall_dim]), min_hall_width(2.1*door_width);
	float min_cell_depth(max(floor_spacing, 2.1f*door_width));
	int const num_end_windows(get_num_windows_on_side(part, !hall_dim));
	float const window_width(part_sz[!hall_dim]/num_end_windows);
	if (num_end_windows >= 3) {max_eq(min_cell_depth, window_width);} // if at least three windows, make sure cells contains a full window
	float const hall_centerline(part.get_center_dim(!hall_dim));
	float min_room_width(2*min_cell_depth + min_hall_width), extra_width(room_width - min_room_width);
	unsigned skip_side(2); // 2=neither
	vector<room_t> &rooms(interior->rooms);
	cube_t hall_room(part);
	if (interior->prison_halls.empty()) {interior->prison_halls = parts;} // start as entire part

	if (extra_width < 0.0) { // too narrow to fit jail cells; should be rare
		if (try_short_dim || part_sz[hall_dim] > 1.8*part_sz[!hall_dim]) { // consider skipping the cells on one side and creating a wider hallway
			min_room_width -= min_cell_depth;
			extra_width     = room_width - min_room_width;
			// if we have enough space, skip side furthest from edge of bcube and keep exterior cells
			if (extra_width > 0.0) {skip_side = (hall_centerline < bcube.get_center_dim(!hall_dim));}
		}
		if (skip_side == 2) { // if we didn't decide to skip a side above
			if (!try_short_dim) {return divide_part_into_jail_cells(part, part_id, rgen, 1);} // try the other dim; try_short_dim=1
			add_room(part, part_id); // add a single room for the entire part
			rooms.back().assign_all_to(RTYPE_JAIL); // maybe shouldn't be a jail in this case?
			return 0;
		}
	}
	float cell_depth(min_cell_depth + min(0.5*min_cell_depth, extra_width/3.0)); // distribute extra width across the hall and cells on either side
	min_eq(cell_depth, 1.2f*max(min_cell_depth, window_width)); // not so wide that we clip a second window
	vect_cube_with_ix_t cells;

	if (in_basement) { // basement jail cell
		float const room_len(part_sz[hall_dim]), min_cell_len(1.3*min_cell_depth), stairs_ext(door_width + wall_thickness);
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
			} // for n
		} // for d
	}
	else { // above ground jail cell
		// start by getting the windows associated with this part
		cube_t window_area(part);
		window_area.expand_by_xy(0.5*wall_thickness); // include ext walls exactly on the edges
		window_area.expand_in_z(-fc_thick); // exclude walls on stacked parts above or below
		vect_vnctcc_t const &wall_quad_verts(get_all_drawn_window_verts_as_quads());
		bool had_ext_wall(0);
		cube_t c;
		float tx1(0.0), tx2(0.0), tz1(0.0), tz2(0.0);

		if (!try_short_dim) { // calculate sum of wall length in each dim and choose the other dim if too unbalanced
			float wall_len[2] = {0,0};

			for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
				if (!get_wall_quad_window_area(wall_quad_verts, i, c, tx1, tx2, tz1, tz2)) continue; // no windows in this wall
				if (!c.intersects(window_area)) continue; // wrong part, skip
				c.intersect_with_cube_xy(window_area); // only windows in this part count
				bool const dim(c.dy() < c.dx());
				wall_len[!dim] += c.get_sz_dim(!dim);
			}
			if (wall_len[long_dim] < 0.67*wall_len[!long_dim]) {return divide_part_into_jail_cells(part, part_id, rgen, 1);} // try_short_dim=1
		}
		for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
			if (!get_wall_quad_window_area(wall_quad_verts, i, c, tx1, tx2, tz1, tz2)) continue; // no windows in this wall
			if (!c.intersects(window_area)) continue; // wrong part, skip
			bool const dim(c.dy() < c.dx()), dir(wall_quad_verts[i].get_norm()[dim] > 0.0);
			if (          dim == hall_dim ) continue; // only keep windows on long part edges
			if ((unsigned)dir == skip_side) continue;
			assert(c.get_sz_dim(dim) == 0.0); // must be zero size in one dim (X or Y oriented); could also use the vertex normal
			had_ext_wall = 1;
			// here we only care about the window width, not the height, because the rooms will span the same floors as the window
			float const window_width(c.get_sz_dim(!dim)/(tx2 - tx1));
			cube_t cell(part);
			float cell_edge(part.d[dim][dir] + (dir ? -1.0 : 1.0)*cell_depth); // room ends at side of central hallway
			// attempt to move to not intersect window
			float const window_hspacing(part_sz[dim]/get_num_windows_on_side(part, dim)); // for hall end windows
			float const shift_cell_edge(shift_val_to_not_intersect_window(part, cell_edge, window_hspacing, get_window_h_border(), dim));
			float const new_cell_depth(fabs(shift_cell_edge - cell.d[dim][dir]));
			float const new_hall_hwidth(shift_cell_edge - hall_centerline + ((skip_side < 2) ? cell_depth : 0.0)); // empty cell row counts as hallway
			if (new_cell_depth > min_cell_depth && 2.0*new_hall_hwidth > min_hall_width) {cell_edge = shift_cell_edge;} // shift if cell and hall are wide enough
			cell.d[dim][!dir] = hall_room.d[dim][dir] = cell_edge;

			for (float xy = tx1; xy < tx2; xy += 1.0) { // windows along each wall
				float const low_edge(c.d[!dim][0] + (xy - tx1)*window_width);
				cell.d[!dim][0] = low_edge;
				cell.d[!dim][1] = low_edge + window_width;
				if (!window_area.contains_cube_xy(cell) || !part.intersects(cell)) continue; // not contained in this part (outside or adjacent to)
				if (interior->is_cube_close_to_doorway(cell, part, door_width, 1)) continue; // inc_open=1
				// at this point there should be no placed stairs, elevators, or ext doors, so no need to check valid placement
				cell.intersect_with_cube(part); // will assert if slightly outside
				cells.emplace_back(cell, (2*dim + dir));
			} // for xy
		} // for i
		if (!had_ext_wall && !try_short_dim) {return divide_part_into_jail_cells(part, part_id, rgen, 1);} // try the other dim; try_short_dim=1
	}
	unsigned const parent_room_id(rooms.size() + cells.size());

	for (cube_with_ix_t const &cell : cells) {
		bool const dim(cell.ix >> 1), dir(cell.ix & 1);
		unsigned const room_id(rooms.size());
		add_room(cell, part_id);
		rooms.back().assign_all_to(RTYPE_JAIL_CELL);
		rooms.back().set_is_nested();
		rooms.back().mark_open_wall(dim, !dir);
		// add side interior walls
		cube_t wall(cell);
		clip_wall_to_ceil_floor(wall, fc_thick);

		for (unsigned d = 0; d < 2; ++d) {
			float const wall_pos(cell.d[!dim][d]);
			if (fabs(wall_pos - part.d[!dim][d]) < wall_thickness) continue; // at edge of part - no wall
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
	add_room(part, part_id); // add a single room for the entire part, last, since cells are sub-rooms
	rooms.back().assign_all_to(RTYPE_JAIL);
	if (cells.empty()) return 0;
	rooms.back().set_has_subroom();
	assert(part_id < interior->prison_halls.size());
	interior->prison_halls[part_id] = hall_room;
	return 1;
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
	cube_t cell(room);
	set_cube_zvals(cell, zval, zval+get_floor_ceil_gap());
	cell.expand_in_dim(dim, -wall_hthick); // exclude front and back
	cube_t cell_inner(cell);
	cell_inner.expand_in_dim(!dim, -wall_hthick); // sides, even against part edges in case there are windows
	
	for (unsigned d = 0; d < 2; ++d) { // exclude interior side walls
		bool const at_ext_wall(fabs(room.d[!dim][d] - part.d[!dim][d]) < 3.0*wall_hthick); // Note: tolerance is needed for interior wall room shrink
		if (at_ext_wall && zval < ground_floor_z1) continue; // no windows in basement, no shrink needed
		cell.d[!dim][d] -= (d ? 1.0 : -1.0)*(at_ext_wall ? 0.8 : 1.0)*wall_hthick; // ext wall needs small shrink to avoid intersecting window bars/trim
	}
	// Note: open wall/bars dim is inverted because the calls below use the hall dim
	float const bars_depth_pos(room.d[dim][!dir] + (dir ? 1.0 : -1.0)*wall_hthick), bars_hthick(JAIL_BARS_THICK*wall_hthick);
	add_jail_cell_bars(cell, cube_t(), room_id, tot_light_amt, !dim, dir, hinge_side, bar_color);
	populate_jail_cell(rgen, cell, zval, room_id, tot_light_amt, !dim, dir, bed_side, sink_on_back_wall, is_lit, bars_hthick, bars_depth_pos);
}

void building_t::add_prison_main_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	// TODO: desk, keys, stairs; also cafeteria, etc.
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
			populate_jail_cell(rgen, cell, zval, room_id, tot_light_amt, dim, dir, bed_side, sink_on_back_wall, is_lit, bars_hthick, bars_depth_pos);
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
	float const door_width(get_doorway_width()), bars_hwidth(JAIL_BARS_THICK*0.5*get_wall_thickness());
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

void building_t::populate_jail_cell(rand_gen_t &rgen, cube_t const &cell, float zval, unsigned room_id, float tot_light_amt,
	bool dim, bool dir, bool bed_side, bool sink_on_back_wall, bool is_lit, float bars_hthick, float bars_depth_pos)
{
	float const floor_spacing(get_window_vspace()), wall_thick(get_wall_thickness()), dsign(dir ? 1.0 : -1.0), bss(bed_side ? 1.0 : -1.0);
	vect_room_object_t &objs(interior->room_geom->objs);
	// add bed
	cube_t bed(cell);
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
	cube_t ts_space(cell); // toilet and sink space
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
		add_bathroom_plumbing(objs.back());
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
	room_t const &parent(get_room(room_id));
	if (!parent.has_subroom()) return 0; // no nested rooms
	cube_t c_exp(c);
	c_exp.expand_by_xy(get_doorway_width()); // leave enough room to the sides for players and AI to walk around

	// nested rooms are added before the parent room, so iterate backwards
	for (int r = (int)room_id-1; r >= 0; --r) {
		room_t const &room(get_room(r));
		if (!room.is_nested() || !parent.contains_cube(room)) break; // no more nested rooms for this parent
		if (c_exp.intersects(room)) return 1;
	}
	return 0;
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
			if ( r.is_ext_basement()) break; // end at extended basement rooms
			if (!r.is_nested())   continue; // not a jail cell
			if ( r.intersects(cand_exp)) return 0; // invalid for either dim
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

void building_t::get_non_jail_cell_non_hallway_cubes(unsigned room_id, vect_cube_t &out) const {
	assert(is_prison());
	room_t const &parent(get_room(room_id));
	assert(parent.part_id < parts.size());
	assert(parent.part_id < interior->prison_halls.size());
	cube_t const &hall(interior->prison_halls[parent.part_id]);
	assert(parent.intersects(hall)); // parent usually contains hall, but may extend slightly outside
	float const wall_thick(get_wall_thickness());
	subtract_cube_from_cube(cube_t(parent), hall, out, 1); // clear_out=1
	if (!parent.has_subroom()) return; // no cells

	for (int r = (int)room_id-1; r >= 0; --r) {
		room_t const &room(get_room(r));
		if (!room.is_nested() || !parent.contains_cube(room)) break; // no more nested rooms for this parent
		cube_t sub(room);
		sub.expand_by_xy(wall_thick); // make sure to remove the walls to the sides of the cells
		subtract_cube_from_cubes(sub, out); // no holes
	}
}

