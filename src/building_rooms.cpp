// 3D World - Building Interior Room Geometry
// by Frank Gennari 4/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#pragma warning(disable : 26812) // prefer enum class over enum


extern int display_mode;

int get_rand_screenshot_texture(unsigned rand_ix);
unsigned get_num_screenshot_tids();


bool building_t::overlaps_other_room_obj(cube_t const &c, unsigned objs_start) const {
	assert(interior && interior->room_geom);
	vector<room_object_t> &objs(interior->room_geom->objs);
	assert(objs_start <= objs.size());

	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (i->intersects(c)) return 1;
	}
	return 0;
}

bool building_t::add_table_and_chairs(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id,
	point const &place_pos, colorRGBA const &chair_color, float rand_place_off, float tot_light_amt, bool is_lit)
{
	float const window_vspacing(get_window_vspace()), room_pad(4.0f*get_wall_thickness());
	uint8_t const obj_flags(is_lit ? RO_FLAG_LIT : 0);
	vector3d const room_sz(room.get_size());
	assert(interior && interior->room_geom);
	vector<room_object_t> &objs(interior->room_geom->objs);
	point table_pos(place_pos);
	vector3d table_sz;
	for (unsigned d = 0; d < 2; ++d) {table_sz [d]  = 0.18*window_vspacing*(1.0 + rgen.rand_float());} // half size relative to window_vspacing
	for (unsigned d = 0; d < 2; ++d) {table_pos[d] += rand_place_off*room_sz[d]*rgen.rand_uniform(-1.0, 1.0);} // near the center of the room
	point llc(table_pos - table_sz), urc(table_pos + table_sz);
	llc.z = table_pos.z; // bottom
	urc.z = table_pos.z + 0.2*window_vspacing;
	cube_t table(llc, urc);
	if (!is_valid_placement_for_room(table, room, blockers, room_pad)) return 0; // check proximity to doors
	objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, obj_flags, tot_light_amt);
	float const chair_sz(0.1*window_vspacing); // half size

	// place some chairs around the table
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (rgen.rand_bool()) continue; // 50% of the time
			point chair_pos(table_pos); // same starting center and z1
			chair_pos[dim] += (dir ? -1.0f : 1.0f)*(rgen.rand_uniform(-0.5, 1.2)*chair_sz + table_sz[dim]);
			cube_t chair(chair_pos, chair_pos);
			chair.z2() += 0.4*window_vspacing; // chair height
			chair.expand_by(vector3d(chair_sz, chair_sz, 0.0));
			if (!is_valid_placement_for_room(chair, room, blockers, room_pad)) continue; // check proximity to doors
			objs.emplace_back(chair, TYPE_CHAIR, room_id, dim, dir, obj_flags, tot_light_amt);
			objs.back().color = chair_color;
		} // for dir
	} // for dim
	return 1;
}

void building_t::add_trashcan_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start, bool check_last_obj) {
	unsigned const NUM_COLORS = 6;
	colorRGBA const colors[NUM_COLORS] = {BLUE, DK_GRAY, LT_GRAY, GRAY, BLUE, WHITE};
	int const rr(rgen.rand()%3), rar(rgen.rand()%3); // three sizes/ARs
	float const radius(0.02f*(3 + rr)*get_window_vspace()), height(0.5f*(3 + rar)*radius);
	vector<room_object_t> &objs(interior->room_geom->objs);
	cube_t room_bounds(get_walkable_room_bounds(room));
	room_bounds.expand_by_xy(-1.1*radius); // leave a slight gap between trashcan and wall
	if (!room_bounds.is_strictly_normalized()) return; // no space for trashcan (likely can't happen)
	point center;
	center.z = zval + 0.001*get_window_vspace(); // slightly above the floor to avoid z-fighting
	unsigned skip_wall(4); // start at an invalid value

	if (check_last_obj) {
		assert(!objs.empty());
		skip_wall = 2*objs.back().dim + (!objs.back().dir); // don't place trashcan on same wall as whiteboard (dir is opposite)
	}
	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a trashcan
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		if ((2*dim + dir) == skip_wall) continue; // don't place a trashcan on this wall
		center[ dim] = room_bounds.d[dim][dir]; // against this wall
		center[!dim] = rgen.rand_uniform(room_bounds.d[!dim][0], room_bounds.d[!dim][1]);
		cube_t c(center, center);
		c.expand_by_xy(radius);
		c.z2() += height;
		if (is_cube_close_to_doorway(c, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(c) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_TCAN, room_id, dim, dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt, room_obj_shape::SHAPE_CYLIN);
		objs.back().color = colors[rgen.rand()%NUM_COLORS];
		return; // done
	} // for n
}

void building_t::add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit, unsigned objs_start) {
	vector<room_object_t> &objs(interior->room_geom->objs);
	cube_t room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return; // room is too small
	float const width(0.35*vspace*rgen.rand_uniform(1.0, 1.2)), depth(0.1*vspace*rgen.rand_uniform(1.0, 1.2)), height(0.7*vspace*rgen.rand_uniform(1.0, 1.2));
	assert(room.part_id < parts.size());
	cube_t const &part(parts[room.part_id]);

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a bookcase
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		if (room.d[dim][dir] == part.d[dim][dir]) continue; // don't place against an exterior wall/window (though check is too strong)
		cube_t c;
		c.z1() = zval;
		c.z2() = zval + height;
		c.d[dim][ dir] = room_bounds.d[dim][dir]; // against this wall
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		c.d[!dim][0] = pos - 0.5*width;
		c.d[!dim][1] = pos + 0.5*width;
		if (is_cube_close_to_doorway(c, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(c) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_BCASE, room_id, dim, !dir, (is_lit ? RO_FLAG_LIT : 0), tot_light_amt); // Note: dir faces into the room, not the wall
		return; // done
	} // for n
}

void building_t::add_rug_to_room(rand_gen_t &rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit) {
	if (!room_object_t::enable_rugs()) return; // disabled
	vector3d const room_sz(room.get_size());
	vector3d center(room.get_cube_center()); // Note: zvals ignored
	bool const min_dim(room_sz.y < room_sz.x);
	float const ar(rgen.rand_uniform(0.65, 0.85)), length(min(0.7f*room_sz[min_dim]/ar, room_sz[!min_dim]*rgen.rand_uniform(0.4, 0.7))), width(length*ar);
	cube_t rug;
	rug.z1() = zval;
	rug.z2() = rug.z1() + 0.001*get_window_vspace(); // almost flat

	for (unsigned d = 0; d < 2; ++d) {
		float const radius(0.5*((bool(d) == min_dim) ? width : length));
		center[d] += 0.05*room_sz[d]*rgen.rand_uniform(-1.0, 1.0); // slight random misalignment
		rug.d[d][0] = center[d] - radius;
		rug.d[d][1] = center[d] + radius;
	}
	uint8_t const obj_flags((is_lit ? RO_FLAG_LIT : 0) | RO_FLAG_NOCOLL);
	vector<room_object_t> &objs(interior->room_geom->objs);
	objs.emplace_back(rug, TYPE_RUG, room_id, 0, 0, obj_flags, tot_light_amt);
	objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix); // determines rug texture
}

bool building_t::hang_pictures_in_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool is_lit) {
	if (!room_object_t::enable_pictures()) return 0; // disabled
	if (!is_house && !room.is_office) return 0; // houses and offices only for now
	if (room.is_sec_bldg) return 0; // no pictures in secondary buildings
	if (room.is_hallway ) return 0; // no pictures in hallways (yet)
	assert(room.part_id < parts.size());
	cube_t const &part(parts[room.part_id]);
	float const floor_height(get_window_vspace()), wall_thickness(get_wall_thickness());
	uint8_t const obj_flags((is_lit ? RO_FLAG_LIT : 0) | RO_FLAG_NOCOLL);
	vector<room_object_t> &objs(interior->room_geom->objs);
	bool was_hung(0);

	if (room.is_office) { // add whiteboards
		if ((rgen.rand() & 3) == 0) return 0; // skip 25% of the time
		bool const pref_dim(rgen.rand_bool()), pref_dir(rgen.rand_bool());

		for (unsigned dim2 = 0; dim2 < 2; ++dim2) {
			for (unsigned dir2 = 0; dir2 < 2; ++dir2) {
				bool const dim(bool(dim2) ^ pref_dim), dir(bool(dir2) ^ pref_dir);
				if (fabs(room.d[dim][dir] - part.d[dim][dir]) < wall_thickness) continue; // on part boundary, likely exterior wall where there may be windows, skip
				float const xy_space(0.2*room.get_sz_dim(!dim));
				cube_t c(room);
				c.z1() = zval + 0.25*floor_height; c.z2() = zval + 0.8*floor_height;
				c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*0.6*wall_thickness; // Note: offset by an additional half wall thickness
				c.d[!dim][0] += xy_space; c.d[!dim][1] -= xy_space;
				if (is_cube_close_to_doorway(c)) continue; // bad placement (TODO: inc_open=1?)
				objs.emplace_back(c, TYPE_WBOARD, room_id, dim, !dir, obj_flags, tot_light_amt); // whiteboard faces dir opposite the wall
				return 1; // done, only need to add one
			} // for dir
		} // for dim
		return 0;
	}
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			float const wall_pos(room.d[dim][dir]);
			if (fabs(room.d[dim][dir] - part.d[dim][dir]) < wall_thickness) continue; // on part boundary, likely exterior wall where there may be windows, skip
			if ((rgen.rand() & 3) == 0) continue; // skip 25% of the time
			float const height(floor_height*rgen.rand_uniform(0.3, 0.6)), width(height*rgen.rand_uniform(1.5, 2.0)); // width > height
			if (width > 0.8*room.get_sz_dim(!dim)) continue; // not enough space
			point center;
			center[ dim] = wall_pos;
			center[!dim] = room.get_center_dim(!dim);
			center.z     = zval + rgen.rand_uniform(0.45, 0.6)*floor_height; // move up

			for (unsigned n = 0; n < 2; ++n) { // make 2 attempts to choose a position along the wall; first iteration is the center
				if (n > 0) {
					float const lo(room.d[!dim][0] + 0.7*width), hi(room.d[!dim][1] - 0.7*width);
					if (hi - lo < width) break; // not enough space to shift, can't place this picture
					center[!dim] = rgen.rand_uniform(lo, hi);
				}
				cube_t c(center, center);
				c.z1() -= 0.5*height; c.z2() += 0.5*height;
				c.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.1*wall_thickness;
				c.d[!dim][0] -= 0.5*width; c.d[!dim][1] += 0.5*width;
				cube_t tc(c);
				tc.d[!dim][0] -= 0.1*width; tc.d[!dim][1] += 0.1*width; // expand slightly to account for frame
				if (is_cube_close_to_doorway(tc) || interior->is_blocked_by_stairs_or_elevator(tc, 4.0*wall_thickness)) continue; // bad placement
				// TODO: check for blocking objects such as bookcases
				objs.emplace_back(c, TYPE_PICTURE, room_id, dim, !dir, obj_flags, tot_light_amt); // picture faces dir opposite the wall
				objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix + 61*dim + 123*dir); // determines picture texture
				was_hung = 1;
				break; // success
			} // for n
		} // for dir
	} // for dim
	return was_hung;
}

void set_light_xy(cube_t &light, point const &center, float light_size, bool light_dim, room_obj_shape light_shape) {
	for (unsigned dim = 0; dim < 2; ++dim) {
		float const sz(((light_shape == SHAPE_CYLIN) ? 1.6 : ((bool(dim) == light_dim) ? 2.2 : 1.0))*light_size);
		light.d[dim][0] = center[dim] - sz;
		light.d[dim][1] = center[dim] + sz;
	}
	light.z1() = light.z2() = center.z; // set so that valid pos can be checked
}

template<typename T> bool has_bcube_int_exp(cube_t const &bcube, vector<T> const &bcubes, float expand) {
	cube_t bcube_exp(bcube);
	bcube_exp.expand_by(expand); // expand in all dirs, including z
	return has_bcube_int(bcube_exp, bcubes);
}

// Note: these three floats can be calculated from mat.get_floor_spacing(), but it's easier to change the constants if we just pass them in
void building_t::gen_room_details(rand_gen_t &rgen, vect_cube_t const &ped_bcubes) {

	assert(interior);
	if (interior->room_geom) return; // already generated?
	//timer_t timer("Gen Room Details");
	interior->room_geom.reset(new building_room_geom_t);
	vector<room_object_t> &objs(interior->room_geom->objs);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	interior->room_geom->obj_scale = window_vspacing; // used to scale room object textures
	unsigned tot_num_rooms(0), num_light_stacks(0);
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {tot_num_rooms += calc_num_floors(*r, window_vspacing, floor_thickness);}
	objs.reserve(tot_num_rooms); // placeholder - there will be more than this many
	room_obj_shape const light_shape(is_house ? SHAPE_CYLIN : SHAPE_CUBE);

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		float const light_amt(window_vspacing*r->get_light_amt()); // multiply perimeter/area by window spacing to make unitless
		float const floor_height(r->is_sec_bldg ? r->dz() : window_vspacing); // secondary buildings are always one floor
		unsigned const num_floors(calc_num_floors(*r, floor_height, floor_thickness));
		unsigned const room_id(r - interior->rooms.begin());
		point room_center(r->get_cube_center());
		// determine light pos and size for this stack of rooms
		bool const room_dim(r->dx() < r->dy()); // longer room dim
		float light_size(floor_thickness); // default size for houses

		if (r->is_office) { // light size varies by office size
			float const room_size(r->dx() + r->dy()); // normalized to office size
			light_size = max(0.015f*room_size, 0.67f*floor_thickness);
		}
		if (r->is_hallway) { // light size varies by hallway size
			float const room_size(min(r->dx(), r->dy())); // normalized to hallway width
			light_size = max(0.06f*room_size, 0.67f*floor_thickness);
		}
		float const light_val(22.0*light_size), room_light_intensity(light_val*light_val/r->get_area_xy()); // average for room, unitless
		cube_t pri_light, sec_light;
		set_light_xy(pri_light, room_center, light_size, room_dim, light_shape);
		if (!r->contains_cube_xy(pri_light)) {pri_light.set_to_zeros();} // disable light if it doesn't fit (small room)
		bool const blocked_by_stairs(!r->is_hallway && interior->is_blocked_by_stairs_or_elevator(pri_light, fc_thick));
		bool use_sec_light(0);
		float z(r->z1());

		if (blocked_by_stairs) { // blocked by stairs - see if we can add a light off to the side in the other orient
			bool const first_dir(rgen.rand_bool());

			for (unsigned d = 0; d < 2; ++d) { // see if we can place it by moving on one direction
				point new_center(room_center);
				new_center[room_dim] += ((bool(d) ^ first_dir) ? -1.0 : 1.0)*0.33*r->get_sz_dim(room_dim);
				set_light_xy(sec_light, new_center, light_size, !room_dim, light_shape); // flip the light dim
				if (!interior->is_blocked_by_stairs_or_elevator(sec_light, fc_thick)) {use_sec_light = 1; break;} // add if not blocked
			}
		}
		// make chair colors consistent for each part by using a few variables for a hash
		colorRGBA chair_colors[12] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN};
		colorRGBA chair_color(chair_colors[(13*r->part_id + 123*tot_num_rooms + 617*mat_ix + 1337*num_floors) % 12]);
		unsigned num_lights_added(0);

		// place objects on each floor for this room
		for (unsigned f = 0; f < num_floors; ++f, z += floor_height) {
			room_center.z = z + fc_thick; // floor height
			bool const top_floor(f+1 == num_floors), check_stairs(!is_house && parts.size() > 1 && top_floor); // top floor of building that may have stairs connecting to upper stack
			bool is_lit(0), light_dim(room_dim), has_stairs(r->has_stairs);

			if (!has_stairs && (f == 0 || top_floor) && interior->stairwells.size() > 1) { // check for stairwells connecting stacked parts
				for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
					if (!r->contains_cube_xy(*s)) continue; // stairs not in this room
					// Note: here we adjust stairs zval by floor_thickness to include stairs in the floor but not in the room above
					if (s->z1() + floor_thickness > r->z2()) continue; // stairs above the room
					if (s->z2() + floor_thickness < r->z1()) continue; // stairs below the room
					has_stairs = 1;
				}
			}
			bool const top_of_stairs(has_stairs && top_floor);
			if (top_of_stairs || check_stairs) {num_light_stacks += num_lights_added;} // these cases may shift the light, so we allocate a new stack
			cube_t light;
			if (!blocked_by_stairs || top_of_stairs) {light = pri_light;}
			else if (use_sec_light) {light = sec_light; light_dim ^= 1;}

			if (!light.is_all_zeros()) { // add a light to the center of the ceiling of this room if there's space (always for top of stairs)
				light.z2() = z + floor_height - fc_thick;
				light.z1() = light.z2() - 0.5*fc_thick;
				is_lit = (r->is_hallway || ((rgen.rand() & (top_of_stairs ? 3 : 1)) != 0)); // 50% of lights are on, 75% for top of stairs, 100% for hallways
				//is_lit |= r->is_sec_bldg; // lit for garages and sheds

				// check ped_bcubes and set is_lit if any are people are in this floor of this room
				for (auto p = ped_bcubes.begin(); p != ped_bcubes.end() && !is_lit; ++p) {
					if (!p->intersects_xy(*r)) continue; // person not in this room
					if (p->z2() < light.z1() && p->z1() + floor_height > light.z2()) {is_lit = 1;} // on this floor
				}
				uint8_t flags(RO_FLAG_NOCOLL); // no collision detection with lights
				if (is_lit)        {flags |= RO_FLAG_LIT;}
				if (top_of_stairs) {flags |= RO_FLAG_TOS;}
				if (has_stairs)    {flags |= RO_FLAG_RSTAIRS;}
				colorRGBA color;
				if (is_house) {color = colorRGBA(1.0, 1.0, 0.85);} // house - yellowish
				else if (r->is_hallway || r->is_office) {color = colorRGBA(0.85, 0.85, 1.0);} // office building - blueish
				else {color = colorRGBA(1.0, 1.0, 1.0);} // white - small office
				unsigned num_lights(r->num_lights);

				if (r->is_hallway && num_lights > 1) { // place a light on each side (of the stairs if they exist), and also between stairs and elevator if there are both
					if (r->has_elevator && r->has_stairs) {num_lights = 3;} // we really should have 3 lights in this case
					float const offset(((num_lights == 3) ? 0.3 : 0.2)*r->get_sz_dim(light_dim)); // closer to the ends in the 3 lights case
					cube_t valid_bounds(*r);
					valid_bounds.expand_by_xy(-0.1*floor_height); // add some padding

					for (unsigned d = 0; d < num_lights; ++d) {
						float const delta((d == 2) ? 0.0 : (d ? -1.0 : 1.0)*offset); // last light is in the center
						cube_t hall_light(light);
						hall_light.translate_dim(delta, light_dim);

						if (check_stairs && has_bcube_int_exp(hall_light, interior->stairwells, fc_thick)) { // keep moving until not blocked by stairs
							cube_t const hall_light_start(hall_light);
							bool is_valid(0);

							for (unsigned shift_dir = 0; shift_dir < 2 && !is_valid; ++shift_dir) {
								hall_light = hall_light_start;

								for (unsigned n = 0; n < 40; ++n) {
									if (!has_bcube_int_exp(hall_light, interior->stairwells, fc_thick)) {is_valid = 1; break;}
									hall_light.translate_dim(0.04*delta*(shift_dir ? -1.0 : 1.0), light_dim);
									if (!valid_bounds.contains_cube_xy(hall_light)) break; // translated outside the hall, give up
								}
							} // for shift_dir
							if (!is_valid) continue; // skip adding this light
						} // end check_stairs
						objs.emplace_back(hall_light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape, color); // dir=0 (unused)
						objs.back().obj_id = num_light_stacks + d;
					} // for d
					num_lights_added = num_lights;
				}
				else { // normal room
					if (check_stairs && has_bcube_int_exp(light, interior->stairwells, fc_thick)) {is_lit = 0;} // disable if blocked by stairs
					else {
						objs.emplace_back(light, TYPE_LIGHT, room_id, light_dim, 0, flags, light_amt, light_shape, color); // dir=0 (unused)
						objs.back().obj_id = num_light_stacks;
					}
					num_lights_added = 1;
				}
				if (is_lit) {r->lit_by_floor |= (1ULL << (f&63));} // flag this floor as being lit (for up to 64 floors)
			} // end light placement
			if (r->no_geom) continue; // no other geometry for this room
			if (has_stairs && !pri_hall.is_all_zeros()) continue; // no other geometry in office building rooms that have stairs
			float tot_light_amt(light_amt); // unitless, somewhere around 1.0
			if (is_lit) {tot_light_amt += room_light_intensity;} // light surface area divided by room surface area with some fudge constant
			unsigned const objs_start(objs.size());

			// place a table and maybe some chairs near the center of the room 95% of the time if it's not a hallway
			bool added_tc(0);

			if (rgen.rand_float() < 0.95) { // 95% of the time
				added_tc = add_table_and_chairs(rgen, *r, ped_bcubes, room_id, room_center, chair_color, 0.1, tot_light_amt, is_lit);
				
				if (!added_tc) {
					// TODO: try to place a desk instead
				}
				if (objs.size() > objs_start) { // an object was placed
					// TODO: maybe add a book on top of objs[objs_start]
				}
			}
			if (is_house) { // place house-specific items
				if (rgen.rand_float() < 0.8) {add_bookcase_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start);}

				if (!has_stairs && (rgen.rand()&3) <= (added_tc ? 0 : 2)) { // maybe add a rug, 25% of the time if there's a table and 75% of the time otherwise
					add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit);
				}
			}
			bool const was_hung(hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit));

			if (rgen.rand_float() < 0.75) { // 75% of the time
				add_trashcan_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, is_lit, objs_start, (was_hung && r->is_office)); // no trashcans on same wall as office whiteboard
			}
			//if (z == bcube.z1()) {} // any special logic that goes on the first floor is here
		} // for f
		num_light_stacks += num_lights_added;
	} // for r
	add_stairs_and_elevators(rgen);
	objs.shrink_to_fit();
	interior->room_geom->light_bcubes.resize(num_light_stacks); // allocate but don't fill un until needed
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness());
	float const stair_dz(window_vspacing/(NUM_STAIRS_PER_FLOOR+1)), stair_height(stair_dz + floor_thickness);
	vector<room_object_t> &objs(interior->room_geom->objs);
	interior->room_geom->stairs_start = objs.size();
	interior->room_geom->has_elevators = (!interior->elevators.empty());

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator) continue; // for elevator, not stairs
		bool const dim(i->dim), dir(i->dir);
		// Note: stairs always start at floor_thickness above the landing z1, ignoring landing z2/height
		float const tot_len(i->get_sz_dim(dim)), floor_z(i->z1() + floor_thickness - window_vspacing), step_len_pos(tot_len/NUM_STAIRS_PER_FLOOR);
		float step_len((dir ? 1.0 : -1.0)*step_len_pos), z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i);

		if (i->shape == SHAPE_STRAIGHT || i->shape == SHAPE_WALLED) { // straight stairs
			for (unsigned n = 0; n < NUM_STAIRS_PER_FLOOR; ++n, z += stair_dz, pos += step_len) {
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // Note: room_id=0, not tracked, unused
			}
		}
		else if (i->shape == SHAPE_U) { // U-shaped stairs
			bool const side(0); // for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
			stair.d[!dim][side] = i->get_center_dim(!dim);
			step_len *= 2.0;

			for (unsigned n = 0; n < NUM_STAIRS_PER_FLOOR; ++n, z += stair_dz, pos += step_len) {
				if (n == NUM_STAIRS_PER_FLOOR/2) { // reverse direction and switch to other side
					step_len *= -1.0;
					stair.d[!dim][ side] = i->d[!dim][side];
					stair.d[!dim][!side] = i->get_center_dim(!dim);
				}
				bool const is_rev(n >= NUM_STAIRS_PER_FLOOR/2);
				stair.d[dim][dir^is_rev^1] = pos; stair.d[dim][dir^is_rev] = pos + step_len;
				stair.z1() = max(floor_z, z); // don't go below the floor
				stair.z2() = z + stair_height;
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, side^is_rev); // Note: room_id=0, not tracked, unused
				objs.back().shape = SHAPE_STAIRS_U;
			} // for n
		}
		else {assert(0);} // unknown stairs shape

		if (i->shape == SHAPE_U || i->shape == SHAPE_WALLED) { // add walls around stairs for this floor
			float const wall_hw(0.15*step_len_pos), half_thick(0.5*floor_thickness);
			stair = *i;
			stair.z2() -= 0.5*floor_thickness; // prevent z-fighting on top floor
			stair.z1()  = max(bcube.z1()+half_thick, floor_z-half_thick); // full height
			set_wall_width(stair, i->d[dim][dir], wall_hw, dim);
			objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // back/end of stairs
			stair.d[dim][!dir] = i->d[dim][!dir];

			for (unsigned d = 0; d < 2; ++d) { // sides of stairs
				set_wall_width(stair, i->d[!dim][d], wall_hw, !dim);
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir);
			}
		}
	} // for i
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		unsigned const elevator_id(i - interior->elevators.begin()); // used for room_object_t::room_id
		cube_t elevator_car(*i);
		elevator_car.z2() = i->z1() + window_vspacing; // currently at the bottom floor
		objs.emplace_back(elevator_car, TYPE_ELEVATOR, elevator_id, i->dim, i->dir, (i->open ? RO_FLAG_OPEN : 0));
	}
}

void building_t::draw_room_geom(shader_t &s, bool shadow_only, bool skip_small_objs) {
	if (interior && interior->room_geom) {interior->room_geom->draw(s, shadow_only, skip_small_objs);}
}
void building_t::gen_and_draw_room_geom(shader_t &s, vect_cube_t &ped_bcubes, unsigned building_ix, int ped_ix, bool shadow_only) {
	if (!interior) return;
	if (is_rotated()) return; // no room geom for rotated buildings

	if (!has_room_geom()) {
		rand_gen_t rgen;
		rgen.set_state(building_ix, parts.size()); // set to something canonical per building
		ped_bcubes.clear();
		if (ped_ix >= 0) {get_ped_bcubes_for_building(ped_ix, building_ix, ped_bcubes);}
		gen_room_details(rgen, ped_bcubes); // generate so that we can draw it
		assert(has_room_geom());
	}
	draw_room_geom(s, shadow_only, 0); // skip_small_objs=0
}

void building_t::clear_room_geom() {
	if (!has_room_geom()) return;
	interior->room_geom->clear(); // free VBO data before deleting the room_geom object
	interior->room_geom.reset();
}

// *** room geometry ***

colorRGBA const WOOD_COLOR(0.9, 0.7, 0.5); // light brown, multiplies wood texture color

// skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2 to match CSG cube flags
void rgeom_mat_t::add_cube_to_verts(cube_t const &c, colorRGBA const &color, unsigned skip_faces, bool swap_tex_st, bool mirror_x) {
	vertex_t v;
	v.set_c4(color);

	// Note: stolen from draw_cube() with tex coord logic, back face culling, etc. removed
	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		unsigned const d[2] = {i, ((i+1)%3)}, n((i+2)%3);

		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			if (skip_faces & (1 << (2*(2-n) + j))) continue; // skip this face
			v.set_ortho_norm(n, j);
			v.v[n] = c.d[n][j];

			for (unsigned s1 = 0; s1 < 2; ++s1) {
				v.v[d[1]] = c.d[d[1]][s1];
				v.t[swap_tex_st] = ((tex.tscale_x == 0.0) ? float(s1) : tex.tscale_x*v.v[d[1]]); // tscale==0.0 => fit texture to cube

				for (unsigned k = 0; k < 2; ++k) { // iterate over vertices
					bool const s2(k^j^s1^1); // need to orient the vertices differently for each side
					v.v[d[0]] = c.d[d[0]][s2];
					v.t[!swap_tex_st] = ((tex.tscale_y == 0.0) ? float(s2) : tex.tscale_y*v.v[d[0]]);
					quad_verts.push_back(v);
					if (mirror_x) {quad_verts.back().t[0] = 1.0 - v.t[0];} // use for pictures
				} // for k
			} // for s1
		} // for j
	} // for i
}

template<typename T> void add_inverted_triangles(T &verts, vector<unsigned> &indices, unsigned verts_start, unsigned ixs_start) {
	unsigned const verts_end(verts.size()), numv(verts_end - verts_start);
	verts.resize(verts_end + numv);

	for (unsigned i = verts_start; i < verts_end; ++i) {
		verts[i+numv] = verts[i];
		verts[i+numv].invert_normal();
	}
	unsigned const ixs_end(indices.size()), numi(ixs_end - ixs_start);
	indices.resize(ixs_end + numi);
	for (unsigned i = 0; i < numi; ++i) {indices[ixs_end + i] = (indices[ixs_end - i - 1] + numv);} // copy in reverse order
}

void rgeom_mat_t::add_vcylin_to_verts(cube_t const &c, colorRGBA const &color, bool draw_bot, bool draw_top, bool two_sided, bool ts_tb, bool inv_tb, float rs_bot, float rs_top) {
	assert(!(ts_tb && inv_tb));
	point const center(c.get_cube_center());
	point const ce[2] = {point(center.x, center.y, c.z1()), point(center.x, center.y, c.z2())};
	unsigned const ndiv(N_CYL_SIDES), num_ends((unsigned)draw_top + (unsigned)draw_bot);
	float const radius(0.5*min(c.dx(), c.dy())), ndiv_inv(1.0/ndiv); // cube X/Y size should be equal/square
	vector3d v12;
	vector_point_norm const &vpn(gen_cylinder_data(ce, radius*rs_bot, radius*rs_top, ndiv, v12));
	color_wrapper const cw(color);
	unsigned itris_start(itri_verts.size()), ixs_start(indices.size()), itix(itris_start), iix(ixs_start);
	itri_verts.resize(itris_start + 2*(ndiv+1));
	indices.resize(ixs_start + 6*ndiv);
	unsigned const ixs_off[6] = {1,2,0, 3,2,1}; // 1 quad = 2 triangles

	for (unsigned i = 0; i <= ndiv; ++i) { // vertex data
		unsigned const s(i%ndiv);
		float const ts(1.0f - i*ndiv_inv);
		norm_comp const normal(0.5*(vpn.n[s] + vpn.n[(i+ndiv-1)%ndiv])); // normalize?
		itri_verts[itix++].assign(vpn.p[(s<<1)+0], normal, ts, 0.0, cw.c);
		itri_verts[itix++].assign(vpn.p[(s<<1)+1], normal, ts, 1.0, cw.c);
	}
	for (unsigned i = 0; i < ndiv; ++i) { // index data
		unsigned const ix0(itris_start + 2*i);
		for (unsigned j = 0; j < 6; ++j) {indices[iix++] = ix0 + ixs_off[j];}
	}
	// room object drawing uses back face culling and single sided lighting; to make lighting two sided, need to add verts with inverted normals/winding dirs
	if (two_sided) {add_inverted_triangles(itri_verts, indices, itris_start, ixs_start);}
	// maybe add top and bottom end cap using triangles, currently using all TCs=0.0
	itris_start = itix = itri_verts.size();
	ixs_start   = iix  = indices.size();
	itri_verts.resize(itris_start + (ndiv + 1)*num_ends);
	indices.resize(ixs_start + 3*ndiv*num_ends);

	for (unsigned bt = 0; bt < 2; ++bt) {
		if (!(bt ? draw_top : draw_bot)) continue; // this disk not drawn
		norm_comp const normal((bool(bt) ^ inv_tb) ? plus_z : -plus_z);
		unsigned const center_ix(itix);
		itri_verts[itix++].assign(ce[bt], normal, 0.0, 0.0, cw.c); // center

		for (unsigned i = 0; i < ndiv; ++i) {
			itri_verts[itix++].assign(vpn.p[(i<<1) + bt], normal, 0.0, 0.0, cw.c);
			indices[iix++] = center_ix; // center
			indices[iix++] = center_ix + i + 1;
			indices[iix++] = center_ix + ((i+1)%ndiv) + 1;
		}
	} // for bt
	if (inv_tb) {std::reverse(indices.begin()+ixs_start, indices.end());} // reverse the order to swap triangle winding order
	if (two_sided) {add_inverted_triangles(itri_verts, indices, itris_start, ixs_start);}
}

class rgeom_alloc_t {
	deque<rgeom_storage_t> free_list; // one per unique texture ID/material
public:
	void alloc(rgeom_storage_t &s) { // attempt to use free_list entry to reuse existing capacity
		if (free_list.empty()) return; // no pre-alloc
		//cout << TXT(free_list.size()) << TXT(free_list.back().get_tot_vert_capacity()) << endl;

		// try to find a free list element with the same tex so that we balance out material memory usage/capacity better
		for (unsigned i = 0; i < free_list.size(); ++i) {
			if (free_list[i].tex.tid != s.tex.tid) continue;
			s.swap(free_list[i]);
			free_list[i].swap(free_list.back());
			free_list.pop_back();
			return; // done
		}
		//s.swap(free_list.back());
		//free_list.pop_back();
	}
	void free(rgeom_storage_t &s) {
		s.clear(); // in case the caller didn't clear it
		free_list.push_back(rgeom_storage_t(s.tex)); // record tex of incoming element
		s.swap(free_list.back()); // transfer existing capacity to free list; clear capacity from s
	}
};

rgeom_alloc_t rgeom_alloc; // static allocator with free list, shared across all buildings; not thread safe

void rgeom_storage_t::clear() {
	quad_verts.clear();
	itri_verts.clear();
	indices.clear();
}
void rgeom_storage_t::swap(rgeom_storage_t &s) {
	quad_verts.swap(s.quad_verts);
	itri_verts.swap(s.itri_verts);
	indices.swap(s.indices);
	std::swap(tex, s.tex);
}

void rgeom_mat_t::clear() {
	vbo.clear();
	delete_and_zero_vbo(ivbo);
	rgeom_storage_t::clear();
	num_qverts = num_itverts = num_ixs = 0;
}

void rgeom_mat_t::create_vbo() {
	num_qverts  = quad_verts.size();
	num_itverts = itri_verts.size();
	num_ixs     = indices.size();
	unsigned qsz(num_qverts*sizeof(vertex_t)), itsz(num_itverts*sizeof(vertex_t));
	vbo.vbo = ::create_vbo();
	check_bind_vbo(vbo.vbo);
	upload_vbo_data(nullptr, get_tot_vert_count()*sizeof(vertex_t));
	upload_vbo_sub_data(quad_verts.data(), 0, qsz);
	upload_vbo_sub_data(itri_verts.data(), qsz, itsz);
	bind_vbo(0);

	if (!indices.empty()) { // we have some indexed quads
		for (auto i = indices.begin(); i != indices.end(); ++i) {*i += num_qverts;} // shift indices to match the new vertex location
		create_vbo_and_upload(ivbo, indices, 1, 1);
	}
	rgeom_alloc.free(*this); // vertex and index data is no longer needed and can be cleared
}

void rgeom_mat_t::draw(shader_t &s, bool shadow_only, bool skip_small_objs) {
	if (shadow_only && !en_shadows)  return; // shadows not enabled for this material (picture, whiteboard, rug, etc.)
	if (shadow_only && tex.emissive) return; // assume this is a light source and shouldn't produce shadows
	assert(vbo.vbo_valid());
	assert(num_qverts > 0 || num_itverts > 0);
	if (!shadow_only) {tex.set_gl(s);} // ignores texture scale for now
	vbo.pre_render();
	vertex_t::set_vbo_arrays();
	if (num_qverts > 0) {draw_quads_as_tris(num_qverts);}

	if (num_itverts > 0 && !skip_small_objs) { // index quads, used for cylinders; treated as small objects
		assert(ivbo > 0);
		bind_vbo(ivbo, 1);
		//glDisable(GL_CULL_FACE); // two sided lighting requires fewer verts (no duplicates), but must be set in the shader
		glDrawRangeElements(GL_TRIANGLES, num_qverts, (num_qverts + num_itverts), num_ixs, GL_UNSIGNED_INT, nullptr);
		//glEnable(GL_CULL_FACE);
		bind_vbo(0, 1);
	}
	if (!shadow_only) {tex.unset_gl(s);}
}

void building_materials_t::clear() {
	for (iterator m = begin(); m != end(); ++m) {m->clear();}
	vector<rgeom_mat_t>::clear();
}
unsigned building_materials_t::count_all_verts() const {
	unsigned num_verts(0);
	for (const_iterator m = begin(); m != end(); ++m) {num_verts += m->get_tot_vert_count();}
	return num_verts;
}
rgeom_mat_t &building_materials_t::get_material(tid_nm_pair_t const &tex, bool inc_shadows) {
	// for now we do a simple linear search because there shouldn't be too many unique materials
	for (iterator m = begin(); m != end(); ++m) {
		if (m->tex != tex) continue;
		if (inc_shadows) {m->enable_shadows();}
		return *m;
	}
	emplace_back(tex); // not found, add a new material
	if (inc_shadows) {back().enable_shadows();}
	rgeom_alloc.alloc(back());
	return back();
}
void building_materials_t::create_vbos() {
	for (iterator m = begin(); m != end(); ++m) {m->create_vbo();}
}
void building_materials_t::draw(shader_t &s, bool shadow_only, bool skip_small_objs) {
	for (iterator m = begin(); m != end(); ++m) {m->draw(s, shadow_only, skip_small_objs);}
}

void get_tc_leg_cubes(cube_t const &c, float width, cube_t cubes[4]) {
	for (unsigned y = 0; y < 2; ++y) {
		for (unsigned x = 0; x < 2; ++x) {
			cube_t leg(c);
			leg.d[0][x] += (1.0f - width)*(x ? -1.0f : 1.0f)*c.dx();
			leg.d[1][y] += (1.0f - width)*(y ? -1.0f : 1.0f)*c.dy();
			cubes[2*y+x] = leg;
		}
	}
}
void building_room_geom_t::add_tc_legs(cube_t const &c, colorRGBA const &color, float width, float tscale) {
	rgeom_mat_t &mat(get_wood_material(tscale));
	cube_t cubes[4];
	get_tc_leg_cubes(c, width, cubes);
	for (unsigned i = 0; i < 4; ++i) {mat.add_cube_to_verts(cubes[i], color, (EF_Z1 | EF_Z2));} // skip top and bottom faces
}

colorRGBA apply_light_color(room_object_t const &o, colorRGBA const &c) {
	if (display_mode & 0x10) return c; // disable this when using indir lighting
	return c * (0.5f + 0.5f*min(sqrt(o.light_amt), 1.5f)); // use c.light_amt as an approximation for ambient lighting due to sun/moon
}
colorRGBA apply_light_color(room_object_t const &o) {return apply_light_color(o, o.color);} // use object color

tid_nm_pair_t const untex_shad_mat(-1, 2.0); // make sure it's different from default tid_nm_pair_t so that it's not grouped with shadowed materials

void building_room_geom_t::add_table(room_object_t const &c, float tscale) { // 6 quads for top + 4 quads per leg = 22 quads = 88 verts
	cube_t top(c), legs_bcube(c);
	top.z1() += 0.85*c.dz(); // 15% of height
	legs_bcube.z2() = top.z1();
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(top, color); // all faces drawn
	add_tc_legs(legs_bcube, color, 0.08, tscale);
}

void building_room_geom_t::add_chair(room_object_t const &c, float tscale) { // 6 quads for seat + 5 quads for back + 4 quads per leg = 27 quads = 108 verts
	float const height(c.dz());
	cube_t seat(c), back(c), legs_bcube(c);
	seat.z1() += 0.32*height;
	seat.z2()  = back.z1() = seat.z1() + 0.07*height;
	legs_bcube.z2() = seat.z1();
	back.d[c.dim][c.dir] += 0.88f*(c.dir ? -1.0f : 1.0f)*c.get_sz_dim(c.dim);
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.2*tscale), 1).add_cube_to_verts(seat, apply_light_color(c)); // all faces drawn
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(back, color, EF_Z1); // skip bottom face
	add_tc_legs(legs_bcube, color, 0.15, tscale);
}

void building_room_geom_t::add_stair(room_object_t const &c, float tscale) {
	get_material(tid_nm_pair_t(MARBLE_TEX, 1.5*tscale), 1).add_cube_to_verts(c, colorRGBA(0.85, 0.85, 0.85)); // all faces drawn
}

unsigned get_face_mask(bool dim, bool dir) {return ~(1 << (2*(2-dim) + dir));} // skip_faces: 1=Z1, 2=Z2, 4=Y1, 8=Y2, 16=X1, 32=X2

tid_nm_pair_t get_tex_auto_nm(int tid, float tscale) {
	return tid_nm_pair_t(tid, get_normal_map_for_bldg_tid(tid), tscale, tscale);
}

void building_room_geom_t::add_elevator(room_object_t const &c, float tscale) {
	// elevator car, all materials are dynamic
	float const thickness(0.051*c.dz());
	cube_t floor(c), ceil(c), back(c);
	floor.z2() = floor.z1() + thickness;
	ceil. z1() = ceil. z2() - thickness;
	floor.expand_by_xy(-0.5f*thickness);
	ceil .expand_by_xy(-0.5f*thickness);
	back.d[c.dim][c.dir] = back.d[c.dim][!c.dir] + (c.dir ? 1.0 : -1.0)*thickness;
	unsigned const front_face_mask(get_face_mask(c.dim, c.dir)), floor_ceil_face_mask(front_face_mask & 60); // +Z faces
	tid_nm_pair_t const paneling(get_tex_auto_nm(PANELING_TEX, 2.0f*tscale));
	get_material(get_tex_auto_nm(TILE_TEX, tscale), 1, 1).add_cube_to_verts(floor, WHITE, floor_ceil_face_mask);
	get_material(get_tex_auto_nm(get_rect_panel_tid(), tscale), 1, 1).add_cube_to_verts(ceil, WHITE, floor_ceil_face_mask);
	get_material(paneling, 1, 1).add_cube_to_verts(back, WHITE, front_face_mask, !c.dim);

	for (unsigned d = 0; d < 2; ++d) { // side walls
		cube_t side(c);
		side.d[!c.dim][!d] = side.d[!c.dim][d] + (d ? -1.0 : 1.0)*thickness;
		get_material(paneling, 1, 1).add_cube_to_verts(side, WHITE, get_face_mask(!c.dim, !d), c.dim);
	}
}

void building_room_geom_t::add_light(room_object_t const &c, float tscale) {
	// Note: need to use a different texture (or -1) for is_on because emissive flag alone does not cause a material change
	bool const is_on(c.is_lit());
	tid_nm_pair_t tp((is_on ? (int)WHITE_TEX : (int)PLASTER_TEX), tscale);
	tp.emissive = is_on;
	rgeom_mat_t &mat(get_material(tp));
	if      (c.shape == SHAPE_CUBE ) {mat.add_cube_to_verts  (c, c.color, EF_Z2);} // untextured, skip top face
	else if (c.shape == SHAPE_CYLIN) {mat.add_vcylin_to_verts(c, c.color, 1, 0);} // bottom only
	else {assert(0);}
}

void building_room_geom_t::add_rug(room_object_t const &c) {
	bool const swap_tex_st(c.dy() < c.dx()); // rug textures are oriented with the long side in X, so swap the coordinates (rotate 90 degrees) if our rug is oriented the other way
	get_material(tid_nm_pair_t(c.get_rug_tid(), 0.0)).add_cube_to_verts(c, WHITE, 61, swap_tex_st); // only draw top/+z face
}

void building_room_geom_t::add_picture(room_object_t const &c) { // also whiteboards
	bool const whiteboard(c.type == TYPE_WBOARD);
	int picture_tid(WHITE_TEX);

	if (!whiteboard) { // picture
		int const user_tid(get_rand_screenshot_texture(c.obj_id));
		picture_tid  = ((user_tid >= 0) ? (unsigned)user_tid : c.get_picture_tid()); // if user texture is valid, use that instead
		num_pic_tids = get_num_screenshot_tids();
		has_pictures = 1;
	}
	unsigned skip_faces(~(1 << (2*(2-c.dim) + c.dir))); // only the face oriented outward
	bool const mirror_x(!whiteboard && !(c.dim ^ c.dir));
	get_material(tid_nm_pair_t(picture_tid, 0.0)).add_cube_to_verts(c, WHITE, skip_faces, !c.dim, mirror_x);
	// add a frame
	cube_t frame(c);
	vector3d exp;
	exp.z = exp[!c.dim] = (whiteboard ? 0.04 : 0.06)*c.dz(); // frame width
	exp[c.dim] = -0.1*c.get_sz_dim(c.dim); // shrink in this dim
	frame.expand_by(exp);
	get_material(tid_nm_pair_t()).add_cube_to_verts(frame, (whiteboard ? GRAY : BLACK), skip_faces, 0);
	
	if (whiteboard) { // add a marker ledge
		cube_t ledge(c);
		ledge.z2() = ledge.z1() + 0.016*c.dz(); // along the bottom edge
		ledge.d[c.dim][c.dir] += (c.dir ? 1.5 : -1.5)*c.get_sz_dim(c.dim); // extrude outward
		get_material(untex_shad_mat, 1).add_cube_to_verts(ledge, GRAY, (1 << (2*(2-c.dim) + !c.dir)), 0); // shadowed
	}
}

void building_room_geom_t::add_book(room_object_t const &c) {
	// TODO - two cubes, use c.color for cover and WHITE for interior pages
	get_material(untex_shad_mat).add_cube_to_verts(c, apply_light_color(c)); // untextured?, all faces
}

void building_room_geom_t::add_bookcase(room_object_t const &c, float tscale) {
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	rgeom_mat_t &mat(get_wood_material(tscale));
	unsigned const skip_faces(!get_face_mask(c.dim, !c.dir)); // skip back face
	cube_t back(c), rest(c);
	back.d[c.dim] [c.dir] += 0.94*(c.dir ? -1.0 : 1.0)*c.get_sz_dim(c.dim);
	rest.d[c.dim][!c.dir]  = back.d[c.dim][c.dir];
	mat.add_cube_to_verts(back, color, skip_faces);
	cube_t top(rest);
	top.z1() += 0.96*c.dz();
	rest.z2() = top.z1();
	mat.add_cube_to_verts(top, color, skip_faces);
	unsigned const num_shelves(3 + (c.obj_id%3)); // 3-5
	cube_t middle(rest);
	
	for (unsigned d = 0; d < 2; ++d) { // left/right sides
		cube_t lr(rest);
		lr.d[!c.dim][d] += 0.94*(d ? -1.0 : 1.0)*c.get_sz_dim(!c.dim);
		mat.add_cube_to_verts(lr, color, skip_faces);
		middle.d[!c.dim][!d] = lr.d[!c.dim][d];
	}
	for (unsigned i = 0; i < num_shelves; ++i) {
		cube_t shelf(middle);
		shelf.z1() += (i+1)*c.dz()/(num_shelves+1);
		shelf.z2()  = shelf.z1() + 0.03*c.dz();
		mat.add_cube_to_verts(shelf, color, skip_faces);
	}
}

void building_room_geom_t::add_desk(room_object_t const &c, float tscale) {
	// desk top and legs, similar to add_table()
	cube_t top(c), legs_bcube(c);
	top.z1() += 0.5*c.dz();
	top.z2()  = top.z1() + 0.1*c.dz();
	legs_bcube.z2() = top.z1();
	colorRGBA const color(apply_light_color(c, WOOD_COLOR));
	get_wood_material(tscale).add_cube_to_verts(top, color); // all faces drawn
	add_tc_legs(legs_bcube, color, 0.10, tscale);
	// TODO - add back/top
	//add_bookcase(c_top_back, tscale);
}

void building_room_geom_t::add_trashcan(room_object_t const &c) {
	rgeom_mat_t &mat(get_material(untex_shad_mat, 1));

	if (c.shape == room_obj_shape::SHAPE_CYLIN) {
		mat.add_vcylin_to_verts(c, apply_light_color(c), 1, 0, 1, 0, 1, 0.7, 1.0); // untextured, bottom only, two_sided cylinder with inverted bottom normal
	}
	else {} // TODO: rounded cube
}

void building_room_geom_t::clear() {
	clear_materials();
	objs.clear();
	light_bcubes.clear();
	has_elevators = 0;
}
void building_room_geom_t::clear_materials() { // can be called to update textures, lighting state, etc.
	materials_s.clear();
	materials_d.clear();
}

rgeom_mat_t &building_room_geom_t::get_wood_material(float tscale) {
	return get_material(get_tex_auto_nm(WOOD2_TEX, tscale), 1); // hard-coded for common material
}

colorRGBA room_object_t::get_color() const {
	switch (type) {
	case TYPE_TABLE:    return WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX));
	case TYPE_CHAIR:    return (color + WOOD_COLOR.modulate_with(texture_color(WOOD2_TEX)))*0.5; // 50% seat color / 50% wood legs color
	case TYPE_STAIR:    return LT_GRAY; // close enough
	case TYPE_ELEVATOR: return LT_BROWN; // ???
	case TYPE_RUG:      return texture_color(get_rug_tid());
	case TYPE_PICTURE:  return texture_color(get_picture_tid());
	case TYPE_WBOARD:   return WHITE;
	case TYPE_BCASE:    return WOOD_COLOR;
	case TYPE_DESK:     return WOOD_COLOR;
	case TYPE_TCAN:     return BLACK;
	default: return color;
	}
	return color; // Note: probably should always set color so that we can return it here
}

void building_room_geom_t::create_static_vbos() {
	//timer_t timer("Gen Room Geom"); // 6.5ms
	float const tscale(2.0/obj_scale);

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible()) continue;
		assert(i->is_strictly_normalized());
		switch (i->type) {
		case TYPE_NONE:  assert(0); // not supported
		case TYPE_TABLE:   add_table   (*i, tscale); break;
		case TYPE_CHAIR:   add_chair   (*i, tscale); break;
		case TYPE_STAIR:   add_stair   (*i, tscale); break;
		case TYPE_LIGHT:   add_light   (*i, tscale); break; // light fixture
		case TYPE_RUG:     add_rug     (*i); break;
		case TYPE_PICTURE: add_picture (*i); break;
		case TYPE_WBOARD:  add_picture (*i); break;
		case TYPE_BOOK:    add_book    (*i); break;
		case TYPE_BCASE:   add_bookcase(*i, tscale); break;
		case TYPE_DESK:    add_desk    (*i, tscale); break;
		case TYPE_TCAN:    add_trashcan(*i); break;
		case TYPE_ELEVATOR: break; // not handled here
		default: assert(0); // undefined type
		}
	} // for i
	// Note: verts are temporary, but cubes are needed for things such as collision detection with the player and ray queries for indir lighting
	//timer_t timer2("Create VBOs"); // 3.8ms
	materials_s.create_vbos();
}

void building_room_geom_t::create_dynamic_vbos() {
	if (!has_elevators) return; // currently only elevators are dynamic, can skip this step if there are no elevators

	for (auto i = objs.begin(); i != objs.end(); ++i) {
		if (!i->is_visible() || i->type != TYPE_ELEVATOR) continue; // only elevators for now
		add_elevator(*i, 2.0/obj_scale);
	}
	materials_d.create_vbos();
}

void building_room_geom_t::draw(shader_t &s, bool shadow_only, bool skip_small_objs) { // non-const because it creates the VBO
	if (empty()) return; // no geom
	unsigned const num_screenshot_tids(get_num_screenshot_tids());

	if (has_pictures && num_pic_tids != num_screenshot_tids) {
		clear_materials(); // user created a new screenshot texture, and this building has pictures - recreate room geom
		num_pic_tids = num_screenshot_tids;
	}
	if (materials_s.empty()) {create_static_vbos ();} // create static  materials if needed
	if (materials_d.empty()) {create_dynamic_vbos();} // create dynamic materials if needed
	enable_blend(); // needed for rugs
	materials_s.draw(s, shadow_only, skip_small_objs);
	materials_d.draw(s, shadow_only, skip_small_objs);
	disable_blend();
	vbo_wrap_t::post_render();
}


