// 3D World - Jail Cells and Prison Buildings
// by Frank Gennari 7/02/25

#include "function_registry.h"
#include "buildings.h"
#include "city_model.h"

extern object_model_loader_t building_obj_model_loader;

bool has_key_3d_model();


bool building_t::add_jail_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
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
	float const door_width(get_doorway_width()), wall_hthick(0.5*get_wall_thickness()), trim_thick(get_trim_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room));
	set_cube_zvals(room_bounds, zval, (zval + get_floor_ceil_gap()));
	float const room_len(room_bounds.get_sz_dim(dim)), min_cell_len(1.25*floor_spacing);
	float const jail_door_width(0.8*door_width); // about the min size that male people/zombies can fit through
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
		float const bars_depth_pos(wall_pos + (dir ? 1.0 : -1.0)*wall_hthick), bars_hthick(0.4*wall_hthick);
		bool const sink_on_back_wall(cell_depth < 0.67*cell_len); // wide and shallow cell
		cube_t cell(cell_area);

		for (unsigned n = 0; n < num_cells; ++n) {
			float const lo_edge(room_bounds.d[dim][0] + n*cell_len), div_wall_hwidth(sink_on_back_wall ? bars_hthick : wall_hthick);
			cell.d[dim][0] = lo_edge;
			cell.d[dim][1] = lo_edge + cell_len;
			if (n   > 0        ) {cell.d[dim][0] += div_wall_hwidth;} // reserve space for walls/bars
			if (n+1 < num_cells) {cell.d[dim][1] -= div_wall_hwidth;}
			// add bars and door
			float const cell_center(cell.get_center_dim(dim));
			bool const hinge_side((room_center < cell_center) ^ bool(dir) ^ 1), bed_side(!hinge_side); // door opens toward hallway center
			float const door_center(cell_center - (bed_side ? 1.0 : -1.0)*0.1*door_width); // slightly away from bed and room door
			cube_t bars(cell);
			set_wall_width(bars, bars_depth_pos, bars_hthick, !dim);
			door_t door(bars, !dim, !dir, 0, 0, hinge_side); // open=0, on_stairs=0
			door.for_jail = 1;
			door.conn_room[0] = door.conn_room[1] = room_id; // both sides connect to the same room
			set_wall_width(door, door_center, 0.5*jail_door_width, dim);
			cube_t bar_segs[2] = {bars, bars};
			bar_segs[0].d[dim][1] = door.d[dim][0]; // lo side
			bar_segs[1].d[dim][0] = door.d[dim][1]; // hi side

			for (unsigned d = 0; d < 2; ++d) { // add bars on both sides of the door
				objs.emplace_back(bar_segs[d], TYPE_JAIL_BARS, room_id, !dim, dir, 0, tot_light_amt, SHAPE_CUBE, bar_color); // dir is facing outside the cell
			}
			door.d[!dim][0] = door.d[!dim][1] = bars_depth_pos; // shrink to zero width
			door.set_for_closet(); // flag so that we don't try to add a light switch by this door, etc.
			unsigned const door_ix(interior->doors.size());
			add_interior_door(door, 0, 1, 1); // is_bathroom=0, make_unlocked=1, make_closed=1

			if (rgen.rand_bool()) {
				add_padlock_to_door(door_ix, (1 << lock_color_ix), rgen); // force lock_color_ix
				added_lock = 1;
			}
			// maybe add divider wall
			if (n > 0) {
				cube_t wall(cell);
				set_wall_width(wall, lo_edge, div_wall_hwidth, dim);

				if (sink_on_back_wall) { // add bars between the cells
					wall.d[!dim][!dir] = bars.d[!dim][!dir]; // flush with the bars
					objs.emplace_back(wall, TYPE_JAIL_BARS, room_id, dim, 0, 0, tot_light_amt, SHAPE_CUBE, bar_color); // dir=0
				}
				else { // sink on side wall; must place a wall to hold the pipes
					unsigned const flags(is_house ? 0 : RO_FLAG_BACKROOM); // flag as backroom for concrete texture in office buildings
					objs.emplace_back(wall, TYPE_PG_WALL, room_id, dim, 0, flags, tot_light_amt, SHAPE_CUBE, WHITE); // dir=0
				}
			}
			populate_jail_cell(rgen, cell, zval, room_id, tot_light_amt, dim, dir, bed_side, sink_on_back_wall, is_lit, bars_hthick, bars_depth_pos, light_color, light_ix_assign);
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

void building_t::populate_jail_cell(rand_gen_t &rgen, cube_t const &cell, float zval, unsigned room_id, float tot_light_amt, bool dim, bool dir, bool bed_side,
	bool sink_on_back_wall, bool is_lit, float bars_hthick, float bars_depth_pos, colorRGBA const &light_color, light_ix_assign_t &light_ix_assign)
{
	float const floor_spacing(get_window_vspace()), wall_thick(get_wall_thickness()), dsign(dir ? 1.0 : -1.0), bss(bed_side ? 1.0 : -1.0);
	vect_room_object_t &objs(interior->room_geom->objs);
	// add a small light in each cell
	cube_t light(cube_top_center(cell));
	light.z1() -= 0.01*floor_spacing;
	light.expand_by_xy(0.06*floor_spacing);
	objs.emplace_back(light, TYPE_LIGHT, room_id, dim, 0, (RO_FLAG_NOCOLL | (is_lit ? RO_FLAG_LIT : 0)), 0.0, SHAPE_CYLIN, light_color); // dir=0 (unused)
	objs.back().obj_id = light_ix_assign.get_next_ix();
	// add bed
	cube_t bed(cell);
	bed.expand_by_xy(-get_trim_thickness()); // add a bit of space around the bed
	bed.z2() = zval + 0.32*floor_spacing; // set height
	float const gap_len(bed.get_sz_dim(!dim)), bed_to_bars_gap(max((gap_len - 0.9*floor_spacing), (bars_hthick + 0.5*wall_thick)));
	bed.d[!dim][!dir] = bars_depth_pos + dsign*bed_to_bars_gap; // set length
	bed.d[ dim][!bed_side] = bed.d[dim][bed_side] - bss*0.5*bed.get_sz_dim(!dim); // set width to half length
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
		add_bathroom_plumbing(objs.back());
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


