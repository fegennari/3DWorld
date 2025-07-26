// 3D World - Building Interior Room Object Placement
// by Frank Gennari 7/15/2023

#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t

bool const PLACE_LIGHTS_ON_SKYLIGHTS = 1;

extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;
extern bldg_obj_type_t bldg_obj_types[];

bool enable_parked_cars();
car_t car_from_parking_space(room_object_t const &o);
void get_stove_burner_locs(room_object_t const &stove, point locs[4]);
void get_pool_ball_rot_matrix(room_object_t const &c, xform_matrix &rot_matrix);
float get_cockroach_height_from_radius(float radius);
void rotate_obj_cube(cube_t &c, cube_t const &bc, bool in_dim, bool dir);
cube_t get_whiteboard_marker_ledge(room_object_t const &c);


class door_path_checker_t {
	vector<point> door_centers;
public:
	bool check_door_path_blocked(cube_t const &c, room_t const &room, unsigned room_id, float zval, building_t const &building) {
		if (door_centers.empty()) {building.get_all_door_centers_for_room(room, room_id, zval, door_centers);}
		if (door_centers.size() < 2) return 0; // must have at least 2 doors for the path to be blocked

		for (auto p1 = door_centers.begin(); p1 != door_centers.end(); ++p1) {
			for (auto p2 = p1+1; p2 != door_centers.end(); ++p2) {
				if (check_line_clip(*p1, *p2, c.d)) return 1;
			}
		}
		return 0;
	}
	void clear() {door_centers.clear();} // to allow for reuse across rooms
};

bool building_t::is_obj_placement_blocked(cube_t const &c, cube_t const &room, bool inc_open_doors, bool check_open_dir, float dmin) const {
	if (is_cube_close_to_doorway(c, room, 0.0, inc_open_doors, check_open_dir)) return 1; // too close to a doorway
	if (interior && interior->is_blocked_by_stairs_or_elevator(c, dmin))        return 1; // faster to check only one per stairwell, but then we need to store another vector?
	if (is_blocked_by_open_attic_door(c)) return 1; // blocked by attic ladder when open
	if (!check_cube_within_part_sides(c)) return 1; // handle non-cube buildings
	return 0;
}
bool building_t::is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, bool inc_open_doors, vector2d const &room_pad) const {
	cube_t c_exp(c);
	c_exp.expand_by_xy(room_pad);
	if (!room.contains_cube_xy(c_exp))                     return 0; // not contained in interior part of the room
	if (is_obj_placement_blocked(c, room, inc_open_doors)) return 0;
	if (has_bcube_int(c, blockers))                        return 0; // Note: ignores dmin
	if (has_attic() && c.intersects_xy(interior->attic_access) && (c.z2() + get_window_vspace()) > interior->attic_access.z1()) return 0; // blocked by attic access door (when open)
	return 1;
}

float get_radius_for_square_model(unsigned model_id) {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id));
	return 0.5f*0.5f*(sz.x + sz.y)/sz.z; // assume square and take average of xsize and ysize
}
cube_t place_cylin_object(rand_gen_t &rgen, cube_t const &place_on, float radius, float height, float dist_from_edge, bool place_at_z1=0) {
	cube_t c;
	gen_xy_pos_for_round_obj(c, place_on, radius, height, dist_from_edge, rgen, place_at_z1); // place at dist_from_edge from edge
	return c;
}

colorRGBA get_pastic_chair_color(colorRGBA const &color) { // colored plastic chairs are too bright; convert to mostly grayscale
	return blend_color(color*0.8, WHITE*color.get_luminance(), 0.2, 0);
}
vector3d building_t::get_office_chair_size() const {
	float const height(0.425*get_window_vspace()), radius(height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR));
	return vector3d(radius, radius, height);
}
bool building_t::add_chair(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id, point const &place_pos, colorRGBA const &chair_color,
	bool dim, bool dir, float tot_light_amt, bool office_chair, bool enable_rotation, bool bar_stool, bool no_push_out, int wooden_or_plastic, bool reduced_clearance)
{
	assert(!(office_chair && bar_stool)); // can't ask for both
	office_chair &= has_office_chair_model();
	bar_stool    &= building_obj_model_loader.is_model_valid(OBJ_MODEL_BAR_STOOL);
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), dir_sign(dir ? -1.0 : 1.0);
	float chair_height(0.0), chair_hwidth(0.0), min_push_out(0.0), max_push_out(0.0);
	point chair_pos(place_pos); // same starting center and z1
	bool is_plastic(0);

	if (office_chair) {
		vector3d const chair_sz(get_office_chair_size());
		chair_height = chair_sz.z;
		chair_hwidth = chair_sz.x;
		min_push_out = 0.5; // pushed out a bit so that the arms don't intersect the table top, but can push out more
		max_push_out = 1.1;
	}
	else if (bar_stool) {
		chair_height = 0.45*window_vspacing;
		chair_hwidth = chair_height*get_radius_for_square_model(OBJ_MODEL_BAR_STOOL);
		min_push_out = 0.25;
		max_push_out = 0.85;
	}
	else { // normal wooden or plastic chair
		is_plastic   = ((wooden_or_plastic == 2) ? !is_residential() : bool(wooden_or_plastic)); // wooden_or_plastic=2 selects based on building type
		chair_height = (is_plastic ? 0.44 : 0.4)*window_vspacing;
		chair_hwidth = 0.1*window_vspacing; // half width
		min_push_out = -0.5; // variable amount of pushed in/out
		max_push_out = 1.2;
	}
	float const min_wall_dist(fabs(place_pos[dim] - room.d[dim][!dir]) - (reduced_clearance ? 1.0 : 1.33)*get_min_front_clearance_inc_people());

	if (!no_push_out) {
		min_eq(max_push_out, (min_wall_dist - chair_hwidth)/chair_hwidth);
		if (min_push_out >= max_push_out) return 0; // not enough space
		chair_pos[dim] += dir_sign*rgen.rand_uniform(min_push_out, max_push_out)*chair_hwidth;
	}
	vector2d room_pad;
	room_pad[ dim] = 4.0*wall_thickness; // more space in front and back
	room_pad[!dim] = 1.0*wall_thickness; // less space to the sides
	cube_t chair(get_cube_height_radius(chair_pos, chair_hwidth, chair_height));
	if (!is_valid_placement_for_room(chair, room, blockers, 1, room_pad)) return 0; // check proximity to doors; inc_open_doors=1
	vect_room_object_t &objs(interior->room_geom->objs);

	if (office_chair) {
		unsigned const flags(enable_rotation ? RO_FLAG_RAND_ROT : 0);
		float lum(chair_color.get_luminance()); // calculate grayscale luminance
		if (lum > 0.5) {lum = 1.0 - lum;} // not white; clamp to [0.0, 0.5] range
		objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, colorRGBA(lum, lum, lum));
	}
	else if (bar_stool) {
		objs.emplace_back(chair, TYPE_BAR_STOOL, room_id, dim, dir, 0, tot_light_amt, SHAPE_CUBE);
	}
	else {
		unsigned flags(0);

		if (place_pos.z < ground_floor_z1 && !has_mall() && rgen.rand_float() < 0.25) { // fallen basement chair 25% of the time; not for mall stores
			// rotate 90 degrees about back legs bottom, tilting backwards
			cube_t new_chair(chair);
			rotate_obj_cube(new_chair, chair, dim, dir);

			if (is_valid_placement_for_room(new_chair, room, blockers, 1, room_pad)) { // has space to fall
				chair  = new_chair;
				flags |= RO_FLAG_ON_FLOOR;
			}
		}
		colorRGBA color(chair_color);
		if (is_plastic) {color = get_pastic_chair_color(color);}
		objs.emplace_back(chair, TYPE_CHAIR, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, color, (is_plastic ? 1 : 0)); // set item_flags to 1 for plastic chair
	}
	return 1;
}

// Note: must be first placed objects; returns the number of total objects added (table + optional chairs)
unsigned building_t::add_table_and_chairs(rand_gen_t rgen, room_t const &room, vect_cube_t &blockers, unsigned room_id, point const &place_pos,
	colorRGBA const &chair_color, float rand_place_off, float tot_light_amt, unsigned max_chairs, bool use_tall_table, int wooden_or_plastic)
{
	bool const use_bar_stools(use_tall_table), is_store(room.is_store());
	float const table_rscale(use_tall_table ? 0.12 : 0.18), room_dx(room.dx()), room_dy(room.dy());
	float const window_vspacing(get_window_vspace()), room_pad(room.is_store() ? 0.0f : max(4.0f*get_wall_thickness(), get_min_front_clearance_inc_people()));
	// use a long table for a large room 50% of the time
	bool const use_long_table(!use_tall_table && max(room_dx, room_dy) > 3.0*window_vspacing && min(room_dx, room_dy) > 2.0*window_vspacing && rgen.rand_bool());
	bool const long_dim(room_dx < room_dy);
	vect_room_object_t &objs(interior->room_geom->objs);
	point table_pos(place_pos);
	vector3d table_sz; // half size
	for (unsigned d = 0; d < 2; ++d) {table_sz [d]  = table_rscale*window_vspacing*(1.0 + rgen.rand_float());} // half size relative to window_vspacing
	for (unsigned d = 0; d < 2; ++d) {table_pos[d] += rand_place_off*room.get_sz_dim(d)*rgen.rand_uniform(-1.0, 1.0);} // near the center of the room
	if (use_long_table) {table_sz[long_dim] *= 1.5;}
	float const long_edge_len(2.0*max(table_sz.x, table_sz.y));
	bool const is_round(!use_tall_table && !use_long_table && rgen.rand_float() < 0.25); // 25% of the time
	bool const is_plastic((wooden_or_plastic == 2) ? (room.is_office && rgen.rand_bool()) : bool(wooden_or_plastic)); // plastic table for offices 50% of the time unless forced
	float const table_hscale(use_tall_table ? 0.38 : (is_round ? 0.23 : (is_plastic ? 0.24 : 0.22)));
	if (is_round) {table_sz.x = table_sz.y = 0.6f*(table_sz.x + table_sz.y);} // round tables must have square bcubes for now (no oval tables yet); make radius slightly larger
	point llc(table_pos - table_sz), urc(table_pos + table_sz);
	llc.z = table_pos.z; // bottom
	urc.z = table_pos.z + table_hscale*rgen.rand_uniform(1.0, 1.1)*window_vspacing; // top
	cube_t table(llc, urc);
	if (!is_valid_placement_for_room(table, room, blockers, 1, vector2d(room_pad, room_pad))) return 0; // check proximity to doors and collision with blockers
	//if (door_path_checker_t().check_door_path_blocked(table, get_room(room_id), room_id, table_pos.z, *this)) return 0; // optional, but may want to allow for kitchens/dining
	unsigned const item_flags(is_plastic ? 1 : 0);
	unsigned flags(is_house ? RO_FLAG_IS_HOUSE : 0);
	// basements can have broken glass tables 50% of the time; these are short cube tables
	if (!is_round && !use_tall_table && !is_plastic && place_pos.z < ground_floor_z1 && (rgen.rand_float() < 0.5)) {flags |= RO_FLAG_BROKEN;}
	objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, flags, tot_light_amt, (is_round ? SHAPE_CYLIN : SHAPE_CUBE), WHITE, item_flags);
	set_obj_id(objs);
	if (max_chairs == 0) return 1; // table only
	// maybe place some chairs around the table
	unsigned const blockers_sz(blockers.size());
	unsigned num_chairs(0);
	bool prev_not_added(0), pri_dim(rgen.rand_bool()), pri_dir(rgen.rand_bool());

	if (use_tall_table) { // use bar stools
		for (unsigned d = 0; d < 2; ++d) {
			bool const dim(pri_dim ^ bool(d));

			for (unsigned e = 0; e < 2; ++e) {
				bool const dir(pri_dir ^ bool(e)); // does pri_dir actually matter here?
				if (num_chairs == max_chairs) break; // done
				point chair_pos(table_pos); // same starting center and z1
				chair_pos[dim] += (dir ? -1.0f : 1.0f)*table_sz[dim];
				bool const added(add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, 0, 0, use_bar_stools)); // office_chair=0
				if (added) {++num_chairs; blockers.push_back(objs.back());}
			}
			if (num_chairs > 0) break; // done
		} // for d
	}
	else {
		bool const two_chairs_long_side(use_long_table && long_edge_len > 3.0*0.2*window_vspacing);
		if (two_chairs_long_side && max_chairs >= 4) {max_chairs = 6;}

		for (unsigned orient = 0; orient < 4; ++orient) {
			bool const dim(bool(orient >> 1) ^ pri_dim), dir(bool(orient & 1) ^ pri_dir);
			// add 2 chairs on long side of table if there's space
			unsigned const num_this_orient((two_chairs_long_side && dim != long_dim) ? 2 : 1);

			for (unsigned n = 0; n < num_this_orient; ++n) {
				if (num_chairs == max_chairs) break; // done
				if (prev_not_added) {prev_not_added = 0;} // if the previous chair failed to be added, make sure to try the next orient
				else if (is_store) {} // don't skip chairs in furniture stores
				else if (orient == 3 && num_chairs == 0 && n+1 == num_this_orient) {} // make sure to place a chair if we have none and this is our last orient/chair
				else if (rgen.rand_float() > 0.25*(num_this_orient + 1)) continue; // skip 50% of the time for single chair, 25% of the time for double
				point chair_pos(table_pos); // same starting center and z1
				chair_pos[dim] += (dir ? -1.0f : 1.0f)*table_sz[dim];
				if (num_this_orient == 2) {chair_pos[!dim] += ((n==0) ? -1.0 : 1.0)*0.225*long_edge_len;} // offset to sides
				bool const added(add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, 0, 0, 0, 0, wooden_or_plastic)); // office_chair=0
				if (added) {++num_chairs; blockers.push_back(objs.back());}
				else {prev_not_added = 1;}
			} // for n
		} // for orient
	} // use chairs
	blockers.resize(blockers_sz); // remove temp chair blockers
	return num_chairs + 1; // add 1 for the table
}
void building_t::shorten_chairs_in_region(cube_t const &region, unsigned objs_start) {
	for (auto i = interior->room_geom->objs.begin() + objs_start; i != interior->room_geom->objs.end(); ++i) {
		if (i->type != TYPE_CHAIR || !i->intersects(region)) continue;
		i->z2() -= 0.25*i->dz();
		i->shape = SHAPE_SHORT;
	}
}

void building_t::get_doorways_for_room(cube_t const &room, float zval, vect_door_stack_t &doorways, bool all_floors) const { // interior doorways; thread safe
	// find interior doorways connected to this room
	float const floor_thickness(get_floor_thickness());
	cube_t room_exp(room);
	room_exp.expand_by_xy(get_wall_thickness());
	if (!all_floors) {set_cube_zvals(room_exp, (zval + floor_thickness), (zval + get_window_vspace() - floor_thickness));} // clip to z-range of this floor
	doorways.clear();

	for (door_stack_t const &ds : interior->door_stacks) {
		if (ds.not_a_room_separator()) continue; // not a real doorway into the room; should is_bldg_conn doors count? depends on the caller?
		if (ds.intersects_no_adj(room_exp)) {doorways.push_back(ds);} // Note: can't use ds.get_conn_room() because this is called before it's filled in
	}
}
vect_door_stack_t &building_t::get_doorways_for_room(cube_t const &room, float zval, bool all_floors) const { // interior doorways; not thread safe
	static vect_door_stack_t doorways; // reuse across rooms
	get_doorways_for_room(room, zval, doorways, all_floors);
	return doorways;
}
void building_t::get_all_door_centers_for_room(room_t const &room, unsigned room_id, float zval, vector<point> &door_centers) const { // used for door_path_checker_t
	float const floor_spacing(get_window_vspace());
	zval += 0.01*floor_spacing; // shift up so that it intersects objects even if they're placed with their z1 exactly at zval
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors
	for (door_stack_t const &ds : doorways) {door_centers.emplace_back(ds.xc(), ds.yc(), zval);}

	if (zval < (ground_floor_z1 + 0.5*floor_spacing)) { // get exterior doors if on the ground floor
		cube_t room_exp(room);
		room_exp.expand_by_xy(get_wall_thickness());

		for (tquad_with_ix_t const &door : doors) {
			cube_t door_bcube(door.get_bcube());
			if (door_bcube.intersects(room_exp)) {door_centers.emplace_back(door_bcube.xc(), door_bcube.yc(), zval);}
		}
	}
	if (room.has_tunnel_conn()) { // tunnel entrance counts as a door, but not false doors
		for (tunnel_seg_t const &t : interior->tunnels) {
			if (t.room_conn && t.conn_room_ix == room_id) {door_centers.push_back(t.get_room_conn_pt(zval));}
		}
	}
}

bool building_t::is_room_an_exit(cube_t const &room, int room_ix, float zval) const { // for living rooms, etc.
	if (is_room_adjacent_to_ext_door(room, zval, 1)) return 1; // front_door_only=1
	if (!multi_family) return 0; // stairs check is only for multi-family houses
	int const has_stairs(room_or_adj_room_has_stairs(room_ix, zval, 1, 0)); // inc_adj_rooms=1, check_door_open=0
	return (has_stairs == 2); // only if adjacent to stairs
}

void building_t::add_trashcan_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool check_last_obj) {
	int const rr(rgen.rand()%3), rar(rgen.rand()%3); // three sizes/ARs
	float const floor_spacing(get_window_vspace()), radius(0.02f*(3 + rr)*floor_spacing), height(0.55f*(3 + rar)*radius); // radius={0.06, 0.08, 0.10} x AR={1.65, 2.2, 2.75}
	cube_t room_bounds(get_walkable_room_bounds(room));
	room_bounds.expand_by_xy(-1.1*radius); // leave a slight gap between trashcan and wall
	if (!room_bounds.is_strictly_normalized()) return; // no space for trashcan (likely can't happen)
	int const floor_ix(int((zval - room.z1())/floor_spacing));
	bool const cylin(((mat_ix + 13*real_num_parts + 5*hallway_dim + 131*floor_ix) % 7) < 4); // varies per-building, per-floor
	point center;
	center.z = zval + 1.1*get_flooring_thick(); // slightly above the flooring/rug to avoid z-fighting
	unsigned skip_wall(4); // start at an invalid value
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t avoid;

	if (!objs.empty() && objs[objs_start].type == TYPE_TABLE) { // make sure there's enough space for the player and AIs to walk around the table
		avoid = objs[objs_start];
		avoid.expand_by_xy(get_min_front_clearance_inc_people());
	}
	if (check_last_obj) {
		assert(!objs.empty());
		skip_wall = 2*objs.back().dim + (!objs.back().dir); // don't place trashcan on same wall as whiteboard (dir is opposite)
	}
	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a trashcan
		bool dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		if ((2U*dim + dir) == skip_wall) {dir ^= 1;} // don't place a trashcan on this wall, try opposite wall
		center[dim] = room_bounds.d[dim][dir]; // against this wall
		bool is_good(0);

		for (unsigned m = 0; m < 40; ++m) { // try to find a point near a doorway
			center[!dim] = rgen.rand_uniform(room_bounds.d[!dim][0], room_bounds.d[!dim][1]);
			if (doorways.empty()) break; // no doorways, keep this point
				
			for (auto i = doorways.begin(); i != doorways.end(); ++i) {
				float const dmin(radius + i->dx() + i->dy()), dist_sq(p2p_dist_sq(center, i->closest_pt(center)));
				if (dist_sq > 4.0*dmin*dmin) continue; // too far
				if (dist_sq <     dmin*dmin) {is_good = 0; break;} // too close, reject this point
				is_good = 1; // close enough, keep this point
			}
			if (is_good) break; // done; may never get here if no points are good, but the code below will handle that
		} // for m
		cube_t const c(get_cube_height_radius(center, radius, height));
		if (!avoid.is_all_zeros() && c.intersects_xy(avoid)) continue; // bad placement
		if (is_obj_placement_blocked(c, room, !room.is_hallway) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_TCAN, room_id, dim, dir, 0, tot_light_amt, (cylin ? SHAPE_CYLIN : SHAPE_CUBE), tcan_colors[rgen.rand() % NUM_TCAN_COLORS]);

		if (is_house || zval < ground_floor_z1 + 2.0*floor_spacing) { // add trash to the trashcan; not on upper floors of office buildings
			unsigned const trash_start(objs.size());
			unsigned const num_objs(rgen.rand() & 3); // 0-3
			float cur_zval(center.z);

			for (unsigned n = 0; n < num_objs; ++n) {
				float trash_radius(min(0.45*radius, 0.2*height)*rgen.rand_uniform(0.8, 1.2)*(1.0 + 0.1*n));
				min_eq(trash_radius, 0.033f*floor_spacing); // limit to a reasonable size
				point trash_center(center.x, center.y, cur_zval+trash_radius);
				if (trash_center.z + 0.5*trash_radius > c.z2()) break;
				if (n > 0) {trash_center += rgen.signed_rand_vector_spherical_xy(0.4*trash_radius);}
				cube_t trash;
				trash.set_from_sphere(trash_center, trash_radius);
				colorRGBA const color(trash_colors[rgen.rand() % NUM_TRASH_COLORS]);
				objs.emplace_back(trash, TYPE_TRASH, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_SPHERE, color);
				set_obj_id(objs);
				cur_zval += 1.5*trash_radius; // overlap a bit
			} // for n
			reverse(objs.begin()+trash_start, objs.end()); // ordered top down, in the order the player must take the trash
		}
		return; // done
	} // for n
}

cube_t get_book_bcube(rand_gen_t &rgen, point const &pos, float floor_spacing, bool dim, bool dir) {
	float const book_sz(0.07*floor_spacing);
	vector3d book_scale(book_sz*rgen.rand_uniform(0.8, 1.2), book_sz*rgen.rand_uniform(0.8, 1.2), 0.0);
	float const thickness(book_sz*rgen.rand_uniform(0.1, 0.3));
	book_scale[dim] *= 0.8; // slightly smaller in this dim
	cube_t book;
	book.set_from_point(pos);
	book.expand_by(book_scale);
	book.z2() += thickness;
	return book;
}

// Note: no blockers, but does check existing objects
bool building_t::add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement) {
	bool const is_store(room.is_store()); // bookstore or furniture store
	cube_t const room_bounds(get_room_bounds_inside_trim(room));
	float const vspace(get_window_vspace()), wall_thickness(get_wall_thickness());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return 0; // room is too small
	rand_gen_t rgen2;
	rgen2.set_state((room_id + 1), (13*mat_ix + interior->rooms.size() + 1)); // local rgen that's per-building/room; ensures bookcases are all the same size in a library
	float const width(0.4*vspace*rgen2.rand_uniform(1.0, 1.2)), depth(0.12*vspace*rgen2.rand_uniform(1.0, 1.2)), height(0.7*vspace*rgen2.rand_uniform(1.0, 1.2));
	float const clearance(max(0.2f*vspace, get_min_front_clearance_inc_people()));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	set_cube_zvals(c, zval, zval+height);

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a bookcase
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		// don't place against an exterior wall/window, inc. partial ext walls
		if (!is_basement && !is_store && classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) continue;
		c.d[dim][ dir] = room_bounds.d[dim][dir]; // against this wall
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		set_wall_width(c, pos, 0.5*width, !dim);
		cube_t tc(c);
		tc.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance; // increase space to add clearance
		if (is_obj_placement_blocked(tc, room, 1) || overlaps_other_room_obj(tc, objs_start)) continue; // bad placement
		tc.d[dim][ dir] -= (dir ? -1.0 : 1.0)*wall_thickness; // expand to wall to check for interior windows and store doors
		if (has_bcube_int(tc, interior->int_windows) || (is_store && has_bcube_int(tc, interior->mall_info->store_doorways))) continue;
		unsigned flags(room.has_open_wall(dim, dir) ? RO_FLAG_OPEN : 0); // flag as open if along an open wall so that the back is drawn

		if (is_basement && rgen.rand_float() < 0.40) { // fallen bookcase 40% of the time
			// rotate 90 degrees about front edge, tilting forward
			cube_t cand(c);
			rotate_obj_cube(cand, c, dim, dir);

			if (!is_obj_placement_blocked(cand, room, 1) && !overlaps_other_room_obj(cand, objs_start)) { // has space to fall
				c     = cand;
				flags = RO_FLAG_ON_FLOOR; // and clear RO_FLAG_OPEN flag
			}
		}
		objs.emplace_back(c, TYPE_BCASE, room_id, dim, !dir, flags, tot_light_amt); // Note: dir faces into the room, not the wall
		set_obj_id(objs);

		if (is_basement) { // maybe add scattered books on the floor
			unsigned num_books(rgen.rand() % 16); // 0-15, but many will fail to be placed
			if (flags & RO_FLAG_ON_FLOOR) {num_books = max(4U, 2U*num_books);} // more books on floor if bookcase has fallen over
			
			if (num_books > 0) {
				float const place_dist(1.0*height);
				point const center(c.get_cube_center());
				cube_t place_area;
				place_area.set_from_sphere(center, place_dist); // zvals are ignored
				place_area.intersect_with_cube(room_bounds);

				for (unsigned n = 0; n < num_books; ++n) {
					point const pos(gen_xy_pos_in_area(place_area, 0.05*vspace, rgen, zval));
					if (!dist_xy_less_than(pos, center, place_dist)) continue;
					bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
					cube_t const book(get_book_bcube(rgen, pos, vspace, dim, dir));
					float const dx(book.dx()), dy(book.dy()), exp(0.5*sqrt(dx*dx + dy*dy)); // book is randomly rotated, to expand to capture all bcubes
					cube_t const bc(pos.x-exp, pos.x+exp, pos.y-exp, pos.y+exp, book.z1(), book.z2());
					if (!room_bounds.contains_cube_xy(bc) || is_obj_placement_blocked(bc, room, 1) || overlaps_other_room_obj(bc, objs_start, 1)) continue;
					colorRGBA const color(book_colors[rgen.rand() % NUM_BOOK_COLORS]);
					objs.emplace_back(book, TYPE_BOOK, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT | RO_FLAG_ON_FLOOR), tot_light_amt, SHAPE_CUBE, color);
					set_obj_id(objs);
				} // for n
			}
		}
		return 1; // done/success
	} // for n
	return 0; // not placed
}

bool building_t::room_has_stairs_or_elevator(room_t const &room, float zval, unsigned floor) const {
	if (room.has_elevator) return 1; // elevator shafts extend through all rooms in a stack, don't need to check zval
	if (!room.has_stairs_on_floor(floor)) return 0; // no stairs
	assert(interior);
	cube_t c(room);
	set_cube_zvals(c, zval, zval+0.9*get_window_vspace());

	for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
		if (s->intersects(c)) return 1;
	}
	return 0;
}
bool building_t::is_room_office_bathroom(room_t &room, float zval, unsigned floor) const { // Note: may also update room flags
	if (!room.is_office || !is_bathroom(room.get_room_type(floor))) return 0;
	if (!room_has_stairs_or_elevator(room, zval, floor))            return 1;
	room.rtype[wrap_room_floor(floor)] = RTYPE_NOTSET; // not a bathroom; can't call assign_to() because it skips bathrooms
	return 0;
}

bool building_t::add_desk_to_room(rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers, colorRGBA const &chair_color, float zval,
	unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement, unsigned desk_ix, bool no_computer, bool force_computer, bool add_phone)
{
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return 0; // room is too small
	float const width(0.8*vspace*rgen.rand_uniform(1.0, 1.2)), depth(0.38*vspace*rgen.rand_uniform(1.0, 1.2)), height(0.23*vspace*rgen.rand_uniform(1.08, 1.2));
	float const clearance(max(0.5f*depth, get_min_front_clearance_inc_people()));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	set_cube_zvals(c, zval, zval+height);

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a desk
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		float const dsign(dir ? -1.0 : 1.0);
		c.d[dim][ dir] = room_bounds.d[dim][dir] + rgen.rand_uniform(0.1, 1.0)*dsign*get_wall_thickness(); // almost against this wall
		c.d[dim][!dir] = c.d[dim][dir] + dsign*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		set_wall_width(c, pos, 0.5*width, !dim);
		cube_t desk_pad(c);
		desk_pad.d[dim][!dir] += dsign*clearance; // ensure clearance in front of the desk so that a chair can be placed
		if (!is_valid_placement_for_room(desk_pad, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
		if (overlaps_other_room_obj(desk_pad, objs_start))             continue; // check other objects (for bedroom desks or multiple office desks)
		// make short if against an exterior or open wall, in an office, or if there's a complex floorplan (in case there's no back wall)
		bool const is_tall(!room.is_office && !has_complex_floorplan && !room.has_open_wall(dim, dir) && rgen.rand_float() < 0.5 &&
			(is_basement || classify_room_wall(room, zval, dim, dir, 0) != ROOM_WALL_EXT));
		room_object_t desk(c, TYPE_DESK, room_id, dim, !dir, (is_house ? RO_FLAG_IS_HOUSE : 0), tot_light_amt, (is_tall ? SHAPE_TALL : SHAPE_CUBE));
		
		if (is_tall) {
			cube_t back(get_desk_top_back(desk));
			back.d[dim][dir] = room_bounds.d[dim][dir]; // against the wall
			if (overlaps_other_room_obj(back, objs_start)) {desk.shape = SHAPE_CUBE;} // back is blocked, maybe be a wall light
		}
		unsigned const desk_obj_ix(objs.size());
		objs.push_back(desk);
		set_obj_id(objs);
		objs.back().obj_id += 123*desk_ix; // set even more differently per-desk so that they have different drawer contents
		bool const add_computer(!no_computer && (force_computer || rgen.rand_bool()) && add_computer_to_desk(desk, desk_obj_ix, dim, dir, rgen, room_id, tot_light_amt));

		if (!add_computer && !room.is_store()) { // no computer; add paper, pens, and pencils; not for furniture stores
			if (rgen.rand_float() < 0.75) {add_papers_to_surface(c, dim, !dir, 7, rgen, room_id, tot_light_amt);} // add 0-7 sheet(s) of paper 75% of the time
			bool const is_big_office(!is_house && room.is_office && interior->rooms.size() > 40);
			unsigned const max_num_pp(is_big_office ? 2 : 3); // 0-3 for houses, 0-2 for big office buildings
			add_pens_pencils_to_surface(c, dim, dir, max_num_pp, rgen, room_id, tot_light_amt);
		}
		if (add_phone) {
			for (unsigned N = 0; N < 10; ++N) { // N tries to place a phone on the desk
				if (!place_phone_on_obj(rgen, desk, room_id, tot_light_amt, desk.dim, desk.dir)) break;
				room_object_t const &phone(objs.back());
				bool valid(1);

				for (unsigned i = desk_obj_ix+1; i+1 < objs.size(); ++i) {
					if (objs[i].intersects(phone)) {valid = 0; break;}
				}
				if (valid) break; // success
				objs.pop_back(); // remove the phone and maybe try again
			} // for N
		}
		if (desk.desk_has_drawers()) { // place blocker in front of drawers so that they have room to open
			room_object_t &desk_obj(objs[desk_obj_ix]);
			room_object_t drawers(get_desk_drawers_part(desk_obj)); // use the actual object, after adding computer, since drawers side depends on obj_id
			drawers.d[dim][!dir] += dsign*0.55*drawers.get_sz_dim(dim); // apply approximate drawer extend
			objs.emplace_back(drawers, TYPE_BLOCKER, room_id, dim, !dir, RO_FLAG_INVIS);
		}
		if (!(is_house || is_office_bldg()) || rgen.rand_float() > 0.05) { // 5% chance of no chair for office buildings and houses
			point chair_pos;
			chair_pos.z = zval;
			chair_pos[ dim] = c.d[dim][!dir];
			chair_pos[!dim] = pos + rgen.rand_uniform(-0.1, 0.1)*width; // slightly misaligned
			// use office chair models when the desk has a computer monitor; now that occlusion culling works well, it's okay to have a ton of these in office buildings
			bool const office_chair(add_computer);
			add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, office_chair);
		}
		return 1; // done/success
	} // for n
	return 0; // failed
}

bool building_t::add_computer_to_desk(cube_t const &desk, unsigned desk_obj_ix, bool dim, bool dir, rand_gen_t &rgen, unsigned room_id, float tot_light_amt, float sz_scale) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_TV)) return 0;
	// add a computer monitor using the TV model
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(desk_obj_ix < objs.size());
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
	float const desk_height(desk.dz()), desk_depth(desk.get_sz_dim(dim)), dsign(dir ? -1.0 : 1.0);
	float const tv_height(1.1*sz_scale*desk_height), tv_hwidth(0.5*tv_height*sz.y/sz.z), tv_depth(tv_height*sz.x/sz.z), center(desk.get_center_dim(!dim));
	cube_t tv;
	set_cube_zvals(tv, desk.z2(), desk.z2()+tv_height);
	tv.d[dim][ dir] = desk.d[dim][dir] + dsign*0.25*desk_depth; // 25% of the way from the wall
	tv.d[dim][!dir] = tv  .d[dim][dir] + dsign*tv_depth;
	set_wall_width(tv, center, tv_hwidth, !dim);
	objs.emplace_back(tv, TYPE_MONITOR, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_SHORT, BLACK); // monitors are shorter than TVs
	set_obj_id(objs);
	// add a keyboard as well
	float const kbd_hwidth(0.7*tv_hwidth), kbd_depth(0.6*kbd_hwidth), kbd_height(0.06*kbd_hwidth);
	cube_t keyboard;
	set_cube_zvals(keyboard, desk.z2(), desk.z2()+kbd_height);
	keyboard.d[dim][!dir] = desk    .d[dim][!dir] - dsign*0.06*desk_depth; // close to front edge
	keyboard.d[dim][ dir] = keyboard.d[dim][!dir] - dsign*kbd_depth;
	set_wall_width(keyboard, center, kbd_hwidth, !dim);
	objs.emplace_back(keyboard, TYPE_KEYBOARD, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // add as white, will be drawn with gray/black texture
	
	if (rgen.rand_float() < 0.75) { // add a mouse on the right side of the desk 75% of the time, if there's space
		bool const  mp_dir(dim ^ dir ^ 1);
		float const mp_width(1.6*kbd_depth), mp_depth(1.25*kbd_depth), mp_height(1.8*kbd_height); // height includes the mouse
		float const end_space(fabs(keyboard.d[!dim][mp_dir] - desk.d[!dim][mp_dir])), slack(end_space - mp_width);

		if (slack > 0.2*mp_width) { // we have the space
			cube_t mp;
			set_cube_zvals(mp, desk.z2(), desk.z2()+mp_height);
			set_wall_width(mp, keyboard.get_center_dim(dim), 0.5*mp_depth, dim);
			mp.d[!dim][!mp_dir] = keyboard.d[!dim][mp_dir] + (mp_dir ? 1.0 : -1.0)*min(0.1f*mp_width, max(0.5f*mp_width, 0.5f*slack));
			mp.d[!dim][ mp_dir] = mp.d[!dim][!mp_dir] + (mp_dir ? 1.0 : -1.0)*mp_width;
			mp.intersect_with_cube_xy(desk); // don't hang off the front of the desk
			unsigned const NUM_MP_COLORS = 5;
			colorRGBA const mp_colors[NUM_MP_COLORS] = {BKGRAY, GRAY_BLACK, DK_GRAY, DK_BLUE, DK_RED};
			objs.emplace_back(mp, TYPE_COMP_MOUSE, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, mp_colors[rgen.rand() % NUM_MP_COLORS]);
			if (rgen.rand_float() < 0.1) {objs.back().taken_level = 1;} // 10% chance of no mouse
		}
	}
	// add a computer tower under the desk
	float const cheight(0.75*sz_scale*desk_height), cwidth(0.44*cheight), cdepth(0.9*cheight); // fixed AR=0.44 to match the texture
	bool const comp_side(rgen.rand_bool());
	float const pos(desk.d[!dim][comp_side] + (comp_side ? -1.0 : 1.0)*0.8*cwidth);
	cube_t computer;
	set_cube_zvals(computer, desk.z1(), desk.z1()+cheight);
	set_wall_width(computer, pos, 0.5*cwidth, !dim);
	computer.d[dim][ dir] = desk    .d[dim][dir] + dsign*0.5*cdepth;
	computer.d[dim][!dir] = computer.d[dim][dir] + dsign*cdepth;
	objs.emplace_back(computer, TYPE_COMPUTER, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);
	room_object_t &desk_obj(objs[desk_obj_ix]);
	desk_obj.flags |= RO_FLAG_ADJ_TOP; // flag so that we don't place other objects on this desk
	// force even/odd-ness of obj_id based on comp_side so that we know what side to put the drawers on so that they don't intersect the computer
	if (bool(desk_obj.obj_id & 1) == comp_side) {++desk_obj.obj_id;}
	return 1;
}

colorRGBA select_paper_color(rand_gen_t &rgen) {
	return paper_colors[rgen.rand()%NUM_PAPER_COLORS];
}
void building_t::add_papers_to_surface(cube_t const &c, bool dim, bool dir, unsigned max_num, rand_gen_t &rgen, unsigned room_id, float tot_light_amt, cube_t const &avoid) {
	float const floor_spacing(get_window_vspace());
	float const plen(0.115*floor_spacing), pwidth(0.77*plen), espace(0.75*plen), thickness(0.00025*floor_spacing); // 8.5x11
	vect_room_object_t &objs(interior->room_geom->objs);
	if (espace > 0.51*c.get_sz_dim(dim) || pwidth > 0.51*c.get_sz_dim(!dim)) return; // desk/table is too small for papers
	cube_t paper;
	set_cube_zvals(paper, c.z2(), c.z2()+thickness); // very thin
	unsigned const num_papers(rgen.rand() % (max_num+1)); // 0 - max_num

	for (unsigned n = 0; n < num_papers; ++n) { // okay if they overlap
		set_wall_width(paper, rgen.rand_uniform(c.d[ dim][0]+espace, c.d[ dim][1]-espace), 0.5*plen,    dim);
		set_wall_width(paper, rgen.rand_uniform(c.d[!dim][0]+pwidth, c.d[!dim][1]-pwidth), 0.5*pwidth, !dim);
		if (!avoid.is_all_zeros() && paper.intersects_xy(avoid)) continue; // skip this paper
		objs.emplace_back(paper, TYPE_PAPER, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt, SHAPE_CUBE, select_paper_color(rgen), btype);
		set_obj_id(objs);
		paper.z2() += thickness; // to avoid Z-fighting if different colors
	} // for n
}
void building_t::add_pens_pencils_to_surface(cube_t const &c, bool dim, bool dir, unsigned max_num, rand_gen_t &rgen, unsigned room_id, float tot_light_amt, cube_t const &avoid) {
	float const floor_spacing(get_window_vspace());
	float const pp_len(0.077*floor_spacing), pp_dia(0.0028*floor_spacing), edge_space(0.75*pp_len); // ~7.5 inches long
	vect_room_object_t &objs(interior->room_geom->objs);
	if (min(c.dx(), c.dy()) < 2.0*edge_space) return; // desk is too small for pens/pencils
	float const pp_z1(c.z2() + 0.3f*pp_dia); // move above papers, and avoid self shadow from the desk
	cube_t pp_bcube;
	set_cube_zvals(pp_bcube, pp_z1, pp_z1+pp_dia);
	unsigned const num_pp(rgen.rand() % (max_num+1)); // 0 - max_num

	for (unsigned n = 0; n < num_pp; ++n) {
		bool const is_pen(rgen.rand_bool());
		colorRGBA const color(is_pen ? pen_colors[rgen.rand()&3] : pencil_colors[rgen.rand()&1]);
		set_wall_width(pp_bcube, rgen.rand_uniform(c.d[ dim][0]+edge_space, c.d[ dim][1]-edge_space), 0.5*pp_len,  dim);
		set_wall_width(pp_bcube, rgen.rand_uniform(c.d[!dim][0]+edge_space, c.d[!dim][1]-edge_space), 0.5*pp_dia, !dim);
		if (!avoid.is_all_zeros() && pp_bcube.intersects_xy(avoid)) continue; // skip this paper
		// Note: no check for overlap with books and potted plants, but that would be complex to add and this case is rare;
		//       computer monitors/keyboards aren't added in this case, and pencils should float above papers, so we don't need to check those
		if (!pp_bcube.is_strictly_normalized()) continue; // too small, likely due to FP error when far from the origin
		objs.emplace_back(pp_bcube, (is_pen ? TYPE_PEN : TYPE_PENCIL), room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, color);
	} // for n
}

bool building_t::add_filing_cabinet_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const fc_height(rgen.rand_uniform(0.45, 0.6)*get_window_vspace());
	vector3d const fc_sz_scale(rgen.rand_uniform(0.40, 0.45), rgen.rand_uniform(0.25, 0.30), 1.0); // depth, width, height
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness());
	return place_obj_along_wall(TYPE_FCABINET, room, fc_height, fc_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 1); // door clearance
}

// Note: can be for an office building or a house
bool building_t::add_office_objs(rand_gen_t rgen, room_t const &room, vect_cube_t &blockers, colorRGBA const &chair_color,
	float zval, unsigned room_id, unsigned floor, float tot_light_amt, unsigned objs_start, bool is_basement)
{
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const desk_obj_id(objs.size());
	if (!add_desk_to_room(rgen, room, blockers, chair_color, zval, room_id, tot_light_amt, objs_start, is_basement)) return 0;

	if (!is_house && rgen.rand_float() < 0.5 && !room_has_stairs_or_elevator(room, zval, floor)) { // allow two desks in one office
		assert(objs[desk_obj_id].type == TYPE_DESK);
		blockers.push_back(objs[desk_obj_id]); // temporarily add the previous desk as a blocker for the new desk and its chair
		room_object_t const &maybe_chair(objs.back());
		bool const added_chair(maybe_chair.type == TYPE_CHAIR || maybe_chair.type == TYPE_OFF_CHAIR);
		if (added_chair) {blockers.push_back(maybe_chair);}
		add_desk_to_room(rgen, room, blockers, chair_color, zval, room_id, tot_light_amt, objs_start, is_basement, 1); // desk_ix=1
		if (added_chair) {blockers.pop_back();} // remove the chair if it was added
		blockers.pop_back(); // remove the first desk blocker
	}
	if (rgen.rand_float() < (is_house ? 0.25 : 0.75)) { // maybe place a filing cabinet along a wall; more likely for office buildings than houses
		if (add_filing_cabinet_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start)) {blockers.push_back(objs.back());}
	}
	if (!is_basement && is_industrial()) {add_industrial_office_objs(rgen, room, zval, room_id, tot_light_amt, objs_start);} // industrial office
	else if (min(room.dx(), room.dy()) > 3.5*get_window_vspace()) { // large room, add a table and chairs in the center
		add_table_and_chairs(rgen, room, blockers, room_id, point(room.xc(), room.yc(), zval), chair_color, 0.1, tot_light_amt);
	}
	return 1;
}

bool building_t::create_office_cubicles(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) { // assumes no prior placed objects
	if (!room.is_office) return 0; // offices only
	if (!room.interior && (rgen.rand()%3) == 0) return 0; // 66.7% chance for non-interior rooms
	cube_t const room_bounds(get_room_bounds_inside_trim(room)); // fix for Z-fighting of cubicles with exterior walls, and avoids clipping through wall trim
	float const floor_spacing(get_window_vspace());
	// Note: we could choose the primary dim based on door placement like in office building bathrooms, but it seems easier to not place cubes by doors
	bool const long_dim(room.dx() < room.dy());
	float const rlength(room_bounds.get_sz_dim(long_dim)), rwidth(room_bounds.get_sz_dim(!long_dim)), midpoint(room_bounds.get_center_dim(!long_dim));
	if (rwidth < 2.5*floor_spacing || rlength < 3.5*floor_spacing) return 0; // not large enough
	float const front_clearance(get_min_front_clearance_inc_people());
	float const target_width(rgen.rand_uniform(0.8, 0.88)*floor_spacing), max_depth(0.5*(rwidth - 1.2*front_clearance));
	unsigned const num_cubes(floor(rlength/target_width)); // >= 3, but usually >= 4
	float const cube_width(rlength/num_cubes), cube_depth(min(cube_width*rgen.rand_uniform(0.8, 1.2), min(1.2f*target_width, max_depth))); // not quite square
	float const trashcan_z1(zval + 1.1*get_flooring_thick()); // slightly above the flooring/rug to avoid z-fighting
	bool const add_middle_col(rwidth > 4.0*cube_depth + 2.0*get_doorway_width()); // enough to fit 4 rows of cubes and 2 hallways in between
	uint16_t const bldg_id(uint16_t(mat_ix + interior->rooms.size())); // some value that's per-building
	cube_t const &part(get_part_for_room(room));
	vect_room_object_t &objs(interior->room_geom->objs);
	bool const has_office_chair(building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR));
	float lo_pos(room_bounds.d[long_dim][0]), chair_height(0.0), chair_radius(0.0);
	cube_t c;
	set_cube_zvals(c, zval, zval+0.425*floor_spacing);
	bool added_cube(0);

	if (has_office_chair) {
		vector3d const chair_sz(get_office_chair_size());
		chair_height = chair_sz.z;
		chair_radius = chair_sz.x;
	}
	for (unsigned n = 0; n < num_cubes; ++n) {
		float const hi_pos(lo_pos + cube_width);
		c.d[long_dim][0] = lo_pos;
		c.d[long_dim][1] = hi_pos;

		for (unsigned is_middle = 0; is_middle < (add_middle_col ? 2U : 1U); ++is_middle) {
			if (is_middle && (n == 0 || n+1 == num_cubes)) continue; // skip end rows for middle section

			for (unsigned dir = 0; dir < 2; ++dir) {
				float const wall_pos(is_middle ? midpoint : room_bounds.d[!long_dim][dir]), dir_sign(dir ? -1.0 : 1.0);
				c.d[!long_dim][ dir] = wall_pos;
				c.d[!long_dim][!dir] = wall_pos + dir_sign*cube_depth;
				cube_t test_cube(c);
				test_cube.d[!long_dim][!dir] += 1.1*dir_sign*front_clearance; // allow space for people to enter the cubicle
				// technically we should check for pillar coll if (num_sides > 4), but in all the cases I've seen the pillar is between rows of cubicles and looks okay
				if (is_obj_placement_blocked(test_cube, room, 1)) continue; // inc_open_doors=1
				bool const against_window(room.d[!long_dim][dir] == part.d[!long_dim][dir]);
				unsigned const cubicle_obj_id(objs.size());
				objs.emplace_back(c, TYPE_CUBICLE, room_id, !long_dim, dir, 0, tot_light_amt, ((against_window && !is_middle) ? SHAPE_SHORT : SHAPE_CUBE));
				objs.back().obj_id = bldg_id;
				added_cube = 1;
				cube_t sides[2], fronts[2], back, surfaces[3];
				get_cubicle_parts(objs[cubicle_obj_id], sides, fronts, back, surfaces);
				bool has_cup(0), has_trashcan(0);

				for (unsigned n = 0; n < 3; ++n) { // add objects to/under cubicle desk surfaces
					bool const sdim(long_dim ^ (n == 2)), sdir((n == 2) ? !dir : bool(!n));
					cube_t const &surface(surfaces[n]);
					unsigned const pp_start(objs.size());
					if (rgen.rand_float() < 0.6) {add_papers_to_surface      (surface, sdim,  sdir, 4, rgen, room_id, tot_light_amt);} // 0-4
					if (rgen.rand_float() < 0.6) {add_pens_pencils_to_surface(surface, sdim, !sdir, 2, rgen, room_id, tot_light_amt);} // 0-2

					if (!has_cup && rgen.rand_float() < 0.25) { // maybe add a cup; at most one per cubicle
						vect_cube_t const avoid(objs.begin()+pp_start, objs.end()); // add all papers, pens, and pencils
						place_cup_on_obj(rgen, surface, room_id, tot_light_amt, avoid);
						has_cup = 1;
					}
					if (rgen.rand_float() < 0.4) { // maybe add a book
						room_object_t const desk(surface, TYPE_DESK, room_id, sdim, sdir, RO_FLAG_INVIS, tot_light_amt); // fake desk for this surface
						place_book_on_obj(rgen, desk, room_id, tot_light_amt, pp_start, 1, 0, 1); // use_dim_dir=1; skip_if_overlaps=1
					}
					if (!has_trashcan && rgen.rand_float() < 0.25) { // maybe add a trashcan under the cubicle
						float const radius(0.5*rgen.rand_uniform(0.6, 0.8)*min(surface.dx(), surface.dy()));
						cube_t tc;
						tc.set_from_point(gen_xy_pos_in_area(surface, radius, rgen));
						set_cube_zvals(tc, trashcan_z1, zval+rgen.rand_uniform(0.5, 0.75)*(surface.z1() - zval));
						tc.expand_by_xy(radius);
						// no trash, no objects to intersect with
						objs.emplace_back(tc, TYPE_TCAN, room_id, sdim, sdir, RO_FLAG_NOCOLL, tot_light_amt,
							(rgen.rand_bool() ? SHAPE_CYLIN : SHAPE_CUBE), tcan_colors[rgen.rand() % NUM_TCAN_COLORS]);
						has_trashcan = 1;
					}
				} // for n
				if (has_office_chair && (rgen.rand()&3)) { // add office chair 75% of the time
					point center(c.get_cube_center());
					center[!long_dim] += dir_sign*0.2*cube_depth;
					for (unsigned d = 0; d < 2; ++d) {center[d] += 0.15*chair_radius*rgen.signed_rand_float();} // slightly random XY position
					center.z = zval;
					cube_t const chair(get_cube_height_radius(center, chair_radius, chair_height));
					objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, !long_dim, dir, RO_FLAG_RAND_ROT, tot_light_amt, SHAPE_CUBE, GRAY_BLACK);
				}
			} // for d
		} // for col
		lo_pos = hi_pos;
	} // for n
	return added_cube;
}

void building_t::add_office_pillars(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned floor_ix, vect_cube_t const &lights, vect_cube_t &blockers) {
	if (room.has_stairs_on_floor(floor_ix)) return; // no pillar if there are stairs since they're likely to intersect; I guess the stairs are structural supports?
	// add one pillar in the center of the room; room should be empty at this point (except for ceiling lights), so no collision checking is needed
	cube_t pillar;
	pillar.set_from_point(point(room.xc(), room.yc(), zval));
	pillar.z2() += get_floor_ceil_gap();
	pillar.expand_by_xy(1.8*get_wall_thickness());
	if (!check_cube_within_part_sides(pillar)) return; // handle non-cube buildings
	if (has_bcube_int(pillar, lights))         return; // really should all the pillar, then the lights, but this is easier (and usually doesn't fail)
	interior->room_geom->objs.emplace_back(pillar, TYPE_OFF_PILLAR, room_id, 0, 0, 0, 1.0, (use_cylinder_coll() ? SHAPE_CYLIN : SHAPE_CUBE)); // cylin buildings have cylin pillars
	blockers.push_back(pillar);
}

void offset_hanging_tv(room_object_t &obj) {
	if (obj.is_hanging()) {obj.translate_dim(obj.dim, (obj.dir ? -1.0 : 1.0)*0.28*obj.get_depth());} // translate to the wall to account for the missing stand
}

vector2d building_t::get_conf_room_table_length_width(cube_t const &room) const {
	float const doorway_width(get_doorway_width()), end_clearance(1.25*max(doorway_width, get_min_front_clearance_inc_people())), side_clearance(1.25*end_clearance);
	bool const dim(room.dx() < room.dy()); // long dim
	float const table_len(room.get_sz_dim(dim) - 2.0*end_clearance), table_width(1.1*doorway_width);
	if (table_len < 1.5*table_width || room.get_sz_dim(!dim) < table_width + 2.0*side_clearance) return vector2d(); // too small to be a conference room
	return vector2d(table_len, table_width);
}
bool building_t::add_conference_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned floor_ix) {
	if (room.has_stairs_on_floor(floor_ix) || room.has_elevator || room.is_hallway) return 0; // not a valid conference room
	float const floor_spacing(get_window_vspace());
	cube_t room_this_floor(room); // clip to a single floor before checking if blocked by stairs or elevator
	set_cube_zvals(room_this_floor, zval, zval+floor_spacing);
	if (interior->is_blocked_by_stairs_or_elevator(room_this_floor)) return 0;
	cube_t const room_bounds(get_walkable_room_bounds(room));
	bool const dim(room.dx() < room.dy()); // long dim
	vect_room_object_t &objs(interior->room_geom->objs);
	// add conference table
	vector2d const table_lw(get_conf_room_table_length_width(room_bounds));
	if (table_lw == vector2d()) return 0; // too small to be a conference room
	cube_t table;
	set_cube_zvals(table, zval, (zval + 0.30*rgen.rand_uniform(1.0, 1.05)*floor_spacing)); // set height
	set_wall_width(table, room.get_center_dim( dim), 0.5*table_lw.x,  dim); // set length
	set_wall_width(table, room.get_center_dim(!dim), 0.5*table_lw.y, !dim); // set width
	objs.emplace_back(table, TYPE_CONF_TABLE, room_id, dim, 0, 0, tot_light_amt, SHAPE_CUBE, WHITE); // dir=0, flags=0
	cube_t avoid; // avoid placing papers and pens/pencils under the phone
	
	if (rgen.rand_float() < 0.9) { // add conference phone with a random dir 90% of the time
		cube_t table_center(table);
		table_center.expand_in_dim( dim, -0.375*table_lw.x); // center 25%
		table_center.expand_in_dim(!dim, -0.250*table_lw.y); // center 50%
		if (place_phone_on_obj(rgen, table_center, room_id, tot_light_amt, !dim, rgen.rand_bool())) {avoid = objs.back();}
	}
	if (rgen.rand_float() < 0.8) {add_papers_to_surface(table, !dim, rgen.rand_bool(), 5, rgen, room_id, tot_light_amt, avoid);} // add 0-5 sheet(s) of paper 80% of the time

	for (unsigned n = 0; n < 2; ++n) { // add twice, so that dim and dir are different
		add_pens_pencils_to_surface(table, rgen.rand_bool(), rgen.rand_bool(), 1, rgen, room_id, tot_light_amt, avoid); // 0-1, random dim and dir
	}
	// add office chairs along sides of the table
	bool add_rolling_chair(has_office_chair_model() /*&& rgen.rand_float() < 0.5*/); // else pastic chairs
	vector3d const chair_sz(add_rolling_chair ? get_office_chair_size() : vector3d(0.1, 0.1, 0.45)*floor_spacing);
	unsigned const num_chairs(floor(table_lw.x/((add_rolling_chair ? 1.2 : 1.8)*2.0*chair_sz.x))); // more spacing for plastic chairs
	colorRGBA const chair_color(add_rolling_chair ? GRAY_BLACK : WHITE*rgen.rand_uniform(0.1, 1.0));
	float const chair_spacing(table_lw.x/max(1U, num_chairs));
	point chair_pos(0.0, 0.0, zval);
	
	for (unsigned side = 0; side < 2; ++side) {
		for (unsigned n = 0; n < num_chairs; ++n) {
			chair_pos[!dim] = table.d[!dim][side] + (side ? 1.0 : -1.0)*chair_sz.x*rgen.rand_uniform(0.2, 0.8);
			chair_pos[ dim] = table.d[ dim][0] + (n + 0.5)*chair_spacing + 0.1*chair_spacing*rgen.signed_rand_float();
			cube_t chair(get_cube_height_radius(chair_pos, chair_sz.x, chair_sz.z));

			if (add_rolling_chair) {
				objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, !dim, !side, 0, tot_light_amt, SHAPE_CYLIN, chair_color);
			}
			else { // plastic chair
				objs.emplace_back(chair, TYPE_CHAIR, room_id, !dim, !side, 0, tot_light_amt, SHAPE_CUBE, chair_color, 1); // item_flags=1 to specify a plastic chair
			}
		} // for n
	} // for side
	// add a big TV/monitor on the wall
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_TV)) {
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
		float const tv_height(0.5*floor_spacing*rgen.rand_uniform(1.0, 1.2)), tv_hwidth(0.5*tv_height*sz.y/sz.z), tv_depth(tv_height*sz.x/sz.z);
		bool dir(rgen.rand_bool());
		cube_t tv;
		tv.z1() = zval    + 0.25*floor_spacing;
		tv.z2() = tv.z1() + tv_height;
		set_wall_width(tv, room_bounds.get_center_dim(!dim), tv_hwidth, !dim);

		for (unsigned n = 0; n < 2; ++n, dir ^= 1) { // find a valid wall dir
			tv.d[dim][ dir] = room_bounds.d[dim][dir]; // on the wall
			tv.d[dim][!dir] = tv.d[dim][dir] + (dir ? -1.0 : 1.0)*tv_depth;
			if (overlaps_or_adj_int_window(tv)) continue; // check interior windows
			if (is_cube_close_to_doorway(tv, room, 0.0, 1)) continue; // too close to a doorway
			add_tv_to_wall(tv, room_id, tot_light_amt, dim, !dir, 1, 0); // use_monitor_image=1, on_off=0 (off)
			break; // done
		} // for n
	}
	add_door_sign("Conference", room, zval, room_id);
	return 1;
}
void building_t::add_tv_to_wall(cube_t const &tv, unsigned room_id, float light_amt, bool dim, bool dir, bool use_monitor_image, int on_off) {
	vect_room_object_t &objs(interior->room_geom->objs);
	// use_monitor_image: want the name, weight, and value of a TV, but we want the images of a computer monitor, so set as TYPE_TV + SHAPE_SHORT
	objs.emplace_back(tv, TYPE_TV, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), light_amt, (use_monitor_image ? SHAPE_SHORT : SHAPE_CUBE), BLACK);
	offset_hanging_tv(objs.back());
	set_obj_id(objs);
	if      (on_off == 0) {objs.back().obj_id |=  1;} // off; set   LSB
	else if (on_off == 1) {objs.back().obj_id &= ~1;} // on ; clear LSB
}

void add_lounge_blockers(vect_room_object_t const &objs, unsigned objs_start, vect_cube_t &blockers) {
	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (i->no_coll() || i->type == TYPE_BLOCKER) continue; // skip pictures, outlets, blocker padding in front of couches, etc.
		blockers.push_back(*i); // add any previously places tables, chairs, etc.
	}
}
cube_t get_table_blocker(vect_room_object_t const &objs, unsigned objs_start) {
	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (i->type != TYPE_TABLE) continue;
		cube_t table_blocker(*i);
		table_blocker.expand_by_xy(2.05*building_t::get_scaled_player_radius());
		return table_blocker; // there should only be one
	}
	return cube_t();
}
void building_t::add_lounge_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_lobby) {
	cube_t const room_area(get_walkable_room_bounds(room)); // right against the wall
	cube_t place_area(room_area);
	place_area.expand_by(-0.25*get_wall_thickness()); // common spacing to wall
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const room_objs_start(objs.size());
	float const window_vspacing(get_window_vspace());
	vect_cube_t blockers;
	add_lounge_blockers(objs, objs_start, blockers); // add any previously places tables, chairs, etc.
	bool const teacher(is_school());
	bool const add_tall_table(!teacher && !is_lobby && min(place_area.dx(), place_area.dy()) > 1.4*window_vspacing && rgen.rand_float() < 0.75); // 75% if non-teacher lounge
	point const table_pos(room.xc(), room.yc(), zval); // approximate; can be placed 10% away from the room center

	if (teacher || add_tall_table) { // add tall table with two bar stools, or a normal table with chairs for a teacher lounge
		unsigned const max_chairs(add_tall_table ? 2 : 4);
		add_table_and_chairs(rgen, room, blockers, room_id, table_pos, WHITE, 0.1, tot_light_amt, max_chairs, add_tall_table);
	}
	if (teacher) {add_mwave_on_table(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area);} // add microwave on a table
	cube_t const table_blocker(get_table_blocker(objs, objs_start));
	if (!table_blocker.is_all_zeros()) {objs.emplace_back(table_blocker, TYPE_BLOCKER, room_id);} // add as an actual blocker for placement of chairs, fridge, etc.

	if (min(room.dx(), room.dy()) > 1.5*window_vspacing) { // place 1-3 couch(es) along a wall if the room is large enough
		unsigned const couches_start(objs.size());
		unsigned const counts[4] = {1, 2, 2, 3}; // 2 couches is more common
		add_couches_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, counts);
		unsigned const couches_end(objs.size());

		if (couches_start < couches_end) { // at least one couch was placed; add tables in front of couches
			add_lounge_blockers(objs, room_objs_start, blockers);
			if (!table_blocker.is_all_zeros()) {blockers.push_back(table_blocker);}

			for (unsigned i = couches_start; i < couches_end; ++i) {
				room_object_t const &c(objs[i]);
				if (c.type != TYPE_COUCH) continue; // skip blockers
				cube_t table(c);
				table.z2() = zval + 0.18*rgen.rand_uniform(1.0, 1.1)*window_vspacing; // top
				table.translate_dim( c.dim, (c.dir ? 1.0 : -1.0)*1.1*c.get_depth());
				if (!is_valid_placement_for_room(table, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
				table.expand_in_dim(!c.dim, -0.25*c.get_width()); // shrink table length
				table.expand_in_dim( c.dim, -0.20*c.get_depth()); // shrink table width
				objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, 0, tot_light_amt, SHAPE_CUBE);
				set_obj_id(objs);
			} // for i
		}
	}
	if (teacher || rgen.rand_bool()) {add_vending_machine(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area);}
	// add a TV on the wall (not on a table)
	add_wall_tv(rgen, room, zval, room_id, tot_light_amt, objs_start);

	if (!is_lobby) { // place 1-2 bookcases
		unsigned const num_bookcases(1 + (rgen.rand() & 1));
		for (unsigned n = 0; n < num_bookcases; ++n) {add_bookcase_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, 0);} // is_basement=0
	}
	if (is_residential()) {add_fishtank_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area);}

	if (!is_lobby) {
		// place fridge; not at window
		place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, room, 0.75, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.2, 4, 0, WHITE, 1);
		// place 0-3 bar stools; fewer if there's already a table
		unsigned const num_stools((rgen.rand() & 3) + !add_tall_table);

		for (unsigned n = 0; n < num_stools; ++n) {
			place_model_along_wall(OBJ_MODEL_BAR_STOOL, TYPE_BAR_STOOL, room, 0.4, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.0, 4, 0);
		}
	}
	// add 1-4 plants
	unsigned const num_plants(1 + (rgen.rand() & 3));
	add_plants_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, num_plants);
	if (teacher) {add_door_sign("Teacher Lounge", room, zval, room_id);}
}

bool building_t::add_mwave_on_table(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, cube_t const &place_area) {
	float const window_vspacing(get_window_vspace()), table_height(rgen.rand_uniform(0.33, 0.36)*window_vspacing);
	float const mheight(rgen.rand_uniform(1.0, 1.2)*0.14*window_vspacing), mwidth(1.7*mheight), mdepth(1.2*mheight); // fixed AR=1.7 to match the texture
	vector3d const sz_scale(rgen.rand_uniform(1.5, 1.8)*mdepth, rgen.rand_uniform(1.5, 1.8)*mwidth, table_height); // depth, width, height
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const table_obj_ix(objs.size());
	if (!place_obj_along_wall(TYPE_TABLE, room, table_height, sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 1)) return 0;
	room_object_t &table(objs[table_obj_ix]);
	bool const dim(table.dim), dir(table.dir);
	float const pos(rgen.rand_uniform((table.d[!dim][0] + 0.6*mwidth), (table.d[!dim][1] - 0.6*mwidth)));
	table.flags |= RO_FLAG_ADJ_TOP; // flag as having a microwave so that we don't add a book or bottle that could overlap it
	cube_t mwave;
	set_cube_zvals(mwave, table.z2(), table.z2()+mheight);
	set_wall_width(mwave, pos, 0.5*mwidth, !dim);
	mwave.d[dim][!dir] = table.d[dim][!dir] + (dir ? 1.0 : -1.0)*0.05*mdepth;
	mwave.d[dim][ dir] = mwave.d[dim][!dir] + (dir ? 1.0 : -1.0)*mdepth;
	objs.emplace_back(mwave, TYPE_MWAVE, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt);
	return 1;
}
bool building_t::add_vending_machine(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, cube_t const &place_area) {
	if (NUM_VEND_TYPES == 0) return 0; // disabled
	unsigned const vtype_id(rgen.rand() % NUM_VEND_TYPES), obj_id(interior->room_geom->objs.size());
	vending_info_t const &vtype(get_vending_type(vtype_id));
	float const floor_spacing(get_window_vspace()), height(0.75*floor_spacing*(vtype.size.z/72)); // normalized to 72"
	if (!place_obj_along_wall(TYPE_VENDING, room, height, vtype.size, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 0, 4, 0, vtype.color)) return 0;
	room_object_t &vm(interior->room_geom->objs[obj_id]);
	vm.item_flags = vtype_id;

	if (vtype.size.z < 48) { // short objects are placed on a small table; this assumes we have space to raise it up
		float const table_height(0.3*floor_spacing);
		cube_t table(vm);
		table.z2() = zval + table_height;
		vm.translate_dim(2, table_height);
		table.expand_by_xy(0.02*vm.get_width()); // slightly larger; hopefully won't intersect anything
		interior->room_geom->objs.emplace_back(table, TYPE_TABLE, room_id, vm.dim, vm.dir, RO_FLAG_ADJ_TOP, tot_light_amt);
	}
	return 1;
}
bool building_t::add_wall_tv(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const tv_obj_ix(objs.size());
	float const z1(zval + 0.3*get_window_vspace());
	cube_t const room_area(get_walkable_room_bounds(room)); // right against the wall
	if (!place_model_along_wall(OBJ_MODEL_TV, TYPE_TV, room, 0.5, rgen, z1, room_id, tot_light_amt, room_area, objs_start, 4.0, 4, 1, BKGRAY, 1, RO_FLAG_HANGING)) return 0;
	offset_hanging_tv(objs[tv_obj_ix]);
	assert(objs.back().type == TYPE_BLOCKER);
	objs.back().z1() = zval; // extend blocker down to the floor so that we don't place objects such as fishtanks under the TV
	return 1;
}

bool building_t::check_valid_closet_placement(cube_t const &c, room_t const &room, unsigned objs_start, unsigned bed_ix, float min_bed_space) const {
	if (min_bed_space > 0.0) {
		room_object_t const &bed(interior->room_geom->get_room_object_by_index(bed_ix));
		assert(bed.type == TYPE_BED);
		cube_t bed_exp(bed);
		bed_exp.expand_by_xy(min_bed_space);
		if (c.intersects_xy(bed_exp)) return 0; // too close to bed
	}
	if (overlaps_other_room_obj(c, objs_start) || is_cube_close_to_doorway(c, room, 0.0, 1)) return 0;
	if (has_attic() && c.intersects_xy(interior->attic_access) && (c.z2() + get_floor_thickness()) > interior->attic_access.z1()) return 0;
	return 1;
}

float get_lamp_width_scale() {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_LAMP)); // L, W, H
	return ((sz == zero_vector) ? 0.0 : 0.5f*(sz.x + sz.y)/sz.z);
}

bool building_t::add_bedroom_objs(rand_gen_t rgen, room_t &room, vect_cube_t &blockers, colorRGBA const &chair_color, float zval, unsigned room_id,
	unsigned floor, float tot_light_amt, unsigned objs_start, bool room_is_lit, bool is_basement, bool force, light_ix_assign_t &light_ix_assign)
{
	// bedrooms should have at least one window; if windowless/interior, it can't be a bedroom; faster than checking count_ext_walls_for_room(room, zval) > 0
	if (room.interior) return 0;
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const bed_obj_ix(objs.size()); // if placed, it will be this index
	int bed_size_ix(-1); // unset
	room_object_t other_bed; // unset
	if (!add_bed_to_room(rgen, room, blockers, zval, room_id, tot_light_amt, floor, force, bed_size_ix, other_bed)) return 0; // it's only a bedroom if there's bed
	assert(bed_obj_ix < objs.size());
	room_object_t const bed(objs[bed_obj_ix]); // deep copy so that we don't need to worry about invalidating the reference below
	float const doorway_width(get_doorway_width()), front_clearance(max(0.6f*doorway_width, get_min_front_clearance_inc_people()));

	if (is_hotel()) { // maybe add a second bed
		unsigned const sec_bed_obj_ix(objs.size());
		cube_t bed_exp(bed);
		bed_exp.expand_by_xy(front_clearance); // add space around the bed so that two beds aren't placed too close together
		blockers.push_back(bed_exp);
		
		if (add_bed_to_room(rgen, room, blockers, zval, room_id, tot_light_amt, floor, 0, bed_size_ix, bed)) { // force=0; pass in the first bed as other_bed
			assert(sec_bed_obj_ix < objs.size());
			room_object_t &bed2(objs[sec_bed_obj_ix]);
			bed2.obj_id = bed.obj_id; // set second bed to the same obj_id as first bed so that they have the same style
			bed2.color  = bed.color ;
		}
		blockers.pop_back(); // remove the bed
	}
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-get_trim_thickness()); // shrink to leave a small gap
	// closet
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), window_h_border(get_window_h_border());
	float const closet_min_depth(0.65*doorway_width), closet_min_width(1.5*doorway_width), min_dist_to_wall(1.0*doorway_width), min_bed_space(front_clearance);
	unsigned const first_corner(rgen.rand() & 3);
	bool const first_dim(rgen.rand_bool()), is_store(room.is_store());
	cube_t const &part(get_part_for_room(room));
	bool placed_closet(0), placed_lamp(0);
	unsigned closet_obj_id(0);
	bool chk_windows[2][2] = {}; // precompute which walls are exterior and can have windows, {dim}x{dir}

	if (!is_basement && has_int_windows()) { // are bedrooms ever placed in the basement?
		for (unsigned d = 0; d < 4; ++d) {chk_windows[d>>1][d&1] = (classify_room_wall(room, zval, (d>>1), (d&1), 0) == ROOM_WALL_EXT);}
	}
	for (unsigned n = 0; n < 4 && !placed_closet; ++n) { // try 4 room corners
		if (is_store) continue; // no closets in furniture stores
		unsigned const corner_ix((first_corner + n)&3);
		bool const xdir(corner_ix&1), ydir(corner_ix>>1);
		point const corner(room_bounds.d[0][xdir], room_bounds.d[1][ydir], zval);

		for (unsigned d = 0; d < 2 && !placed_closet; ++d) { // try both dims
			bool const dim(bool(d) ^ first_dim), dir(dim ? ydir : xdir), other_dir(dim ? xdir : ydir);
			if (room_bounds.get_sz_dim(!dim) < closet_min_width + min_dist_to_wall) continue; // room is too narrow to add a closet here
			if (chk_windows[dim][dir]) continue; // don't place closets against exterior walls where they would block a window
			float const dir_sign(dir ? -1.0 : 1.0), signed_front_clearance(dir_sign*front_clearance);
			float const window_hspacing(get_hspacing_for_part(part, dim));
			cube_t c(corner, corner);
			c.d[0][!xdir] += (xdir ? -1.0 : 1.0)*(dim ? closet_min_width : closet_min_depth);
			c.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(dim ? closet_min_depth : closet_min_width);
			if (chk_windows[!dim][other_dir] && is_val_inside_window(part, dim, c.d[dim][!dir], window_hspacing, window_h_border)) continue; // check for window intersection
			c.z2() += window_vspacing - floor_thickness;
			c.d[dim][!dir] += signed_front_clearance; // extra padding in front, to avoid placing too close to bed, etc.
			if (!check_valid_closet_placement(c, room, objs_start, bed_obj_ix, min_bed_space)) continue; // bad placement
			// good placement, see if we can make the closet larger
			unsigned const num_steps = 10;
			float const req_dist(chk_windows[!dim][!other_dir] ? (other_dir ? -1.0 : 1.0)*min_dist_to_wall : 0.0); // signed; at least min dist from opposite wall if exterior
			float max_grow((room_bounds.d[!dim][!other_dir] - req_dist) - c.d[!dim][!other_dir]);
			min_eq(max_grow, 1.5f*window_vspacing); // limit to a reasonable length
			float const len_step(max_grow/num_steps), depth_step(dir_sign*0.35*doorway_width/num_steps); // signed

			for (unsigned s1 = 0; s1 < num_steps; ++s1) { // try increasing width
				cube_t c2(c);
				c2.d[!dim][!other_dir] += len_step;
				if (!check_valid_closet_placement(c2, room, objs_start, bed_obj_ix, min_bed_space)) break; // bad placement
				c = c2; // valid placement, update with larger cube
			}
			if (front_clearance < doorway_width && room_object_t(c, TYPE_CLOSET, room_id, dim, !dir).is_small_closet()) {
				// will be a small closet; add extra space in front to make the sure the door can open
				cube_t c_ext(c);
				c_ext.d[dim][!dir] += dir_sign*(doorway_width - front_clearance);
				if (!check_valid_closet_placement(c_ext, room, objs_start, bed_obj_ix, 0.0)) continue; // bad placement; don't need to recheck bed
			}
			for (unsigned s2 = 0; s2 < num_steps; ++s2) { // now try increasing depth
				cube_t c2(c);
				c2.d[dim][!dir] += depth_step;
				if (chk_windows[!dim][other_dir] && is_val_inside_window(part, dim, (c2.d[dim][!dir] - signed_front_clearance),
					window_hspacing, get_window_h_border())) break; // bad placement
				if (!check_valid_closet_placement(c2, room, objs_start, bed_obj_ix, min_bed_space)) break; // bad placement
				c = c2; // valid placement, update with larger cube
			}
			c.d[dim][!dir] -= signed_front_clearance; // subtract off front clearance
			assert(c.is_strictly_normalized());
			unsigned flags(0);
			if (is_house) {flags |= RO_FLAG_IS_HOUSE;}
			if (c.d[!dim][0] == room_bounds.d[!dim][0]) {flags |= RO_FLAG_ADJ_LO;}
			if (c.d[!dim][1] == room_bounds.d[!dim][1]) {flags |= RO_FLAG_ADJ_HI;}
			if (is_hotel()) {flags |= RO_FLAG_HAS_EXTRA;} // flag so that the closet has a safe
			//if ((rgen.rand() % 10) == 0) {flags |= RO_FLAG_OPEN;} // 10% chance of open closet; unclear if this adds any value, but it works
			closet_obj_id = objs.size();
			objs.emplace_back(c, TYPE_CLOSET, room_id, dim, !dir, flags, tot_light_amt, SHAPE_CUBE, wall_color); // closet door is always white; sides should match interior walls
			set_obj_id(objs);
			if (flags & RO_FLAG_OPEN) {interior->room_geom->expand_object(objs.back(), *this);} // expand opened closets immediately
			placed_closet = 1; // done
			// add a light inside the closet
			room_object_t const closet(objs.back()); // deep copy so that we can invalidate the reference
			bool const small_closet(closet.is_small_closet());
			point lpos(cube_top_center(closet));
			lpos[dim] += 0.05*c.get_sz_dim(dim)*(dir ? -1.0 : 1.0); // move slightly toward the front of the closet
			cube_t light(lpos);
			light.z1() -= 0.02*window_vspacing;
			light.expand_by_xy((small_closet ? 0.04 : 0.06)*window_vspacing);
			colorRGBA const color(1.0, 1.0, 0.9); // yellow-ish
			objs.emplace_back(light, TYPE_LIGHT, room_id, dim, 0, (RO_FLAG_NOCOLL | RO_FLAG_IN_CLOSET), 0.0, SHAPE_CYLIN, color); // dir=0 (unused)
			objs.back().obj_id = light_ix_assign.get_next_ix();

			if (small_closet) { // add a blocker in front of the closet to avoid placing furniture that blocks the door from opening
				c.d[dim][!dir] += dir_sign*doorway_width;
				objs.emplace_back(c, TYPE_BLOCKER, room_id, dim, 0, RO_FLAG_INVIS);
				// add closet door
				cube_t cubes[5];
				get_closet_cubes(closet, cubes);
				door_t door(cubes[4], dim, !dir, 0); // open=0
				door.d[dim][0] = door.d[dim][1] = door.get_center_dim(dim); // shrink to zero width
				door.set_for_closet(); // flag so that we don't try to add a light switch by this door, etc.
				add_interior_door(door, 0, 1, 1); // is_bathroom=0, make_unlocked=1, make_closed=1
				interior->doors.back().obj_ix = closet_obj_id;
			}
		} // for d
	} // for n
	// dresser
	float const ds_height(rgen.rand_uniform(0.26, 0.32)*window_vspacing), ds_depth(rgen.rand_uniform(0.20, 0.25)*window_vspacing), ds_width(rgen.rand_uniform(0.6, 0.9)*window_vspacing);
	vector3d const ds_sz_scale(ds_depth/ds_height, ds_width/ds_height, 1.0);
	unsigned const dresser_obj_id(objs.size());
	
	if (place_obj_along_wall(TYPE_DRESSER, room, ds_height, ds_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 1)) { // door clearance
		room_object_t &dresser(objs[dresser_obj_id]);

		// place a mirror on the dresser 25% of the time; skip if against an exterior wall to avoid blocking a window
		if (rgen.rand_float() < 0.25 && classify_room_wall(room, zval, dresser.dim, !dresser.dir, 0) != ROOM_WALL_EXT) {
			room_object_t mirror(dresser);
			mirror.type = TYPE_DRESS_MIR;
			set_cube_zvals(mirror, dresser.z2(), (dresser.z2() + 1.4*dresser.get_height()));
			mirror.d[mirror.dim][mirror.dir] -= (mirror.dir ? 1.0 : -1.0)*0.9*dresser.get_length(); // push it toward the back
			mirror.expand_in_dim(!mirror.dim, -0.02*mirror.get_width()); // shrink slightly
			if (is_residential()) {mirror.flags |= RO_FLAG_IS_HOUSE;} // flag as in a house for reflections logic; should always be true?
			//mirror .flags |= RO_FLAG_NOCOLL; // leave this unset so that light switches aren't blocked, etc.
			dresser.flags |= RO_FLAG_ADJ_TOP; // flag the dresser as having an item on it so that we don't add something else that blocks or intersects the mirror
			objs.push_back(mirror);
			set_obj_id(objs); // for crack texture selection/orient
			room.set_has_mirror();
		}
	}
	// nightstand
	unsigned const pref_orient(2*bed.dim + (!bed.dir)); // prefer the same orient as the bed so that it's placed on the same wall next to the bed
	float const ns_height(rgen.rand_uniform(0.24, 0.26)*window_vspacing), ns_depth(rgen.rand_uniform(0.15, 0.2)*window_vspacing), ns_width(rgen.rand_uniform(1.0, 2.0)*ns_depth);
	vector3d const ns_sz_scale(ns_depth/ns_height, ns_width/ns_height, 1.0);
	place_obj_along_wall(TYPE_NIGHTSTAND, room, ns_height, ns_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 1, pref_orient); // door clearance

	if (placed_closet) { // determine if there's space for the closet doors to fold outward
		room_object_t &closet(objs[closet_obj_id]);

		if (closet.get_sz_dim(!closet.dim) < 1.8*closet.dz()) { // only for medium sized closets
			bool const dim(closet.dim), dir(closet.dir);
			cube_t doors_area(closet);
			doors_area.d[dim][!dir]  = closet.d[dim][dir]; // flush with the front of the closet
			doors_area.d[dim][ dir] += (dir ? 1.0 : -1.0)*0.25*closet.get_sz_dim(!dim); // extend outward by a quarter the closet width
			bool can_fold((room_bounds.d[dim][dir] < doors_area.d[dim][dir]) ^ dir); // should be true, unless closet is very wide and room is very narrow

			for (auto i = objs.begin()+objs_start; i != objs.end() && can_fold; ++i) {
				if (i->type == TYPE_CLOSET || i->type == TYPE_LIGHT) continue; // skip the closet and its light
				can_fold &= !i->intersects(doors_area);
			}
			if (can_fold) { // mark as folding
				closet.flags |= RO_FLAG_HANGING;
				objs.emplace_back(doors_area, TYPE_BLOCKER, room_id, dim, dir, RO_FLAG_INVIS); // prevent adding bookcases/trashcans/balls intersecting open closet doors
			}
		}
	}
	// try to place a lamp on a dresser or nightstand that was added to this room
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_LAMP) && (rgen.rand()&3) != 0) {
		float const height(0.25*window_vspacing), width(height*get_lamp_width_scale());
		point pillow_center(bed.get_cube_center());
		pillow_center[bed.dim] += (bed.dir ? 1.0 : -1.0)*0.5*bed.get_sz_dim(bed.dim); // adjust from bed center to near the pillow(s)
		int obj_id(-1);
		float dmin_sq(0.0);

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) { // choose the dresser or nightstand closest to the bed
			if (i->type != TYPE_DRESSER && i->type != TYPE_NIGHTSTAND) continue; // not a dresser or nightstand
			if (i->flags & RO_FLAG_ADJ_TOP)    continue; // don't add a lamp if there's a mirror on it
			if (min(i->dx(), i->dy()) < width) continue; // too small to place a lamp on
			float const dist_sq(p2p_dist_xy_sq(i->get_cube_center(), pillow_center));
			if (dmin_sq == 0.0 || dist_sq < dmin_sq) {obj_id = (i - objs.begin()); dmin_sq = dist_sq;}
		}
		if (obj_id >= 0) { // found a valid object to place this on
			room_object_t &obj(objs[obj_id]);
			point center(obj.get_cube_center());
			center.z = obj.z2();
			cube_t lamp(get_cube_height_radius(center, 0.5*width, height));
			lamp.translate_dim(obj.dim, (obj.dir ? 1.0 : -1.0)*0.1*width); // move slightly toward the front to avoid clipping through the wall
			float const shift_range(0.4f*(obj.get_sz_dim(!obj.dim) - width)), obj_center(obj.get_center_dim(!obj.dim)), targ_pos(pillow_center[!obj.dim]);
			float shift_val(0.0), dmin(0.0);

			for (unsigned n = 0; n < 4; ++n) { // generate several random positions on the top of the object and choose the one closest to the bed
				float const cand_shift(rgen.rand_uniform(-1.0, 1.0)*shift_range), dist(fabs((obj_center + cand_shift) - targ_pos));
				if (dmin == 0.0 || dist < dmin) {shift_val = cand_shift; dmin = dist;}
			}
			lamp.translate_dim(!obj.dim, shift_val);
			unsigned flags(RO_FLAG_NOCOLL); // no collisions, as an optimization since the player and AI can't get onto the dresser/nightstand anyway
			if (rgen.rand_bool() && !room_is_lit) {flags |= RO_FLAG_LIT;} // 50% chance of being lit if the room is dark (Note: don't let room_is_lit affect rgen)
			obj.flags |= RO_FLAG_ADJ_TOP; // flag this object as having something on it
			objs.emplace_back(lamp, TYPE_LAMP, room_id, obj.dim, obj.dir, flags, tot_light_amt, SHAPE_CYLIN, lamp_colors[rgen.rand()%NUM_LAMP_COLORS]); // Note: invalidates obj ref
			placed_lamp = 1;
		}
	}
	if (!placed_lamp && !is_store && rgen.rand_float() < 0.75) { // attempt to place a lava lamp on a dresser or nightstand 75% of the time
		float const height(0.15*window_vspacing), radius(0.135*height);
		bool const on_dresser(rgen.rand() & 3); // dresser 75% of the time

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (i->type != (on_dresser ? TYPE_DRESSER : TYPE_NIGHTSTAND)) continue; // not a dresser or nightstand
			cube_t place_area(*i);
			if (i->flags & RO_FLAG_ADJ_TOP) {place_area.d[i->dim][!i->dir] = i->get_center_dim(i->dim);} // shrink place area to exclude the mirror
			if (min(place_area.dx(), place_area.dy()) < 3.0*radius) continue; // too small to place a lava lamp on
			cube_t const llamp(place_cylin_object(rgen, place_area, radius, height, 1.0*radius));
			i->flags |= RO_FLAG_ADJ_TOP; // flag this object as having something on it
			objs.emplace_back(llamp, TYPE_LAVALAMP, room_id, i->dim, i->dir, (RO_FLAG_NOCOLL | RO_FLAG_LIT), tot_light_amt, SHAPE_CYLIN, LT_GRAY); // on by default
			placed_lamp = 1; // counts as a lamp
			break; // must break here as we've invalidated the iterator i
		} // for i
	}
	if (min(room_bounds.dx(), room_bounds.dy()) > 2.5*window_vspacing && max(room_bounds.dx(), room_bounds.dy()) > 3.0*window_vspacing) {
		// large room, try to add a desk and chair as well
		add_desk_to_room(rgen, room, blockers, chair_color, zval, room_id, tot_light_amt, objs_start, is_basement);
	}
	if (is_store) return 1; // no flashlight, candle, balls, t-shirt, or jeans in furniture store rooms

	if (is_hotel()) { // hotel specific items
		// add a phone on a dresser, nightstand, or desk
		static vector<unsigned> cands;
		cands.clear();

		for (unsigned i = objs_start; i < objs.size(); ++i) {
			room_object_t const &obj(objs[i]);
			if (obj.flags & RO_FLAG_ADJ_TOP) continue; // skip if an item has already been placed
			unsigned num_slots(0);
			if      (obj.type == TYPE_DRESSER   ) {num_slots = 2;}
			else if (obj.type == TYPE_NIGHTSTAND) {num_slots = 1;}
			else if (obj.type == TYPE_DESK      ) {num_slots = 3;}
			cands.resize(cands.size()+num_slots, i);
		} // for i
		if (!cands.empty()) {
			room_object_t &obj(objs[cands[rgen.rand() % cands.size()]]);
			if (place_phone_on_obj(rgen, obj, room_id, tot_light_amt, obj.dim, obj.dir)) {obj.flags |= RO_FLAG_ADJ_TOP;}
		}
	}
	else { // not a hotel
		// maybe add a flashlight or candle on a dresser, night stand, or desk; or in a drawer?
		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (!((i->type == TYPE_DRESSER || i->type == TYPE_NIGHTSTAND || i->type == TYPE_DESK) && !(i->flags & RO_FLAG_ADJ_TOP))) continue; // not empty dresser/nightstand/desk
			unsigned const rand_val(rgen.rand());
			if (rand_val & 3) continue; // only add 25% of the time
			bool const is_flashlight(rand_val & 4);
			float const height((is_flashlight ? 0.1 : 0.09)*window_vspacing), radius((is_flashlight ? 0.2 : 0.16)*height);
			if (min(i->dx(), i->dy()) < 3.0*radius) continue; // surface is too small
			i->flags |= RO_FLAG_ADJ_TOP;
			cube_t bc;
			gen_xy_pos_for_round_obj(bc, *i, radius, height, 1.2*radius, rgen);
			colorRGBA const color(is_flashlight ? BLACK : candle_color);
			objs.emplace_back(bc, (is_flashlight ? TYPE_FLASHLIGHT : TYPE_CANDLE), room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, color);
			break; // only add one object; adding invalidates the iterator anyway
		} // for i
		if (rgen.rand_float() < 0.3) {add_laundry_basket(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area);} // try to place a laundry basket 25% of the time

		if (rgen.rand_probability(global_building_params.ball_prob)) { // maybe add a ball to the room
			add_ball_to_room(rgen, room, place_area, zval, room_id, tot_light_amt, objs_start);
		}
	}
	cube_t const avoid(placed_closet ? objs[closet_obj_id] : cube_t()); // avoid intersecting the closet, since it meets the ceiling
	
	if (objs_start > 0) { // light is prev placed object; should always be true
		if (!(bed.flags & RO_FLAG_ADJ_BOT)) {replace_light_with_ceiling_fan(rgen, room, avoid, room_id, tot_light_amt, objs_start-1);} // add if no bunk bed
	}
	if (rgen.rand_float() < 0.3) { // maybe add a t-shirt or jeans on the floor
		unsigned const type(rgen.rand_bool() ? TYPE_PANTS : TYPE_TEESHIRT);
		bool already_on_bed(0); // if shirt/pants are already on the bed, don't put them on the floor

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (i->type == type) {already_on_bed = 1; break;}
		}
		if (!already_on_bed) {place_shirt_pants_on_floor(rgen, room, zval, room_id, tot_light_amt, place_area, objs_start, type);}
	}
	return 1; // success
} // end add_bedroom_objs()

void building_t::place_shirt_pants_on_floor(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt,
	cube_t const &place_area, unsigned objs_start, unsigned type)
{
	float const floor_spacing(get_window_vspace());
	float const length(((type == TYPE_TEESHIRT) ? 0.26 : 0.2)*floor_spacing), width(0.98*length), height(0.002*floor_spacing);
	cube_t valid_area(place_area);
	valid_area.expand_by_xy(-0.25*floor_spacing); // not too close to a wall to avoid bookcases, dressers, and nightstands
	bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random orientation
	vector3d size(0.5*length, 0.5*width, height);
	if (dim) {swap(size.x, size.y);}
	if (valid_area.dx() < 2.0*size.x || valid_area.dy() < 2.0*size.y) return; // too small to place; should never be true

	for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to place the object
		cube_t c(gen_xy_pos_in_area(valid_area, size, rgen, zval));
		c.expand_by_xy(size);
		c.z2() += height;
		if (overlaps_other_room_obj(c, objs_start) || is_obj_placement_blocked(c, room, 1)) continue; // bad placement
		interior->room_geom->objs.emplace_back(c, type, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, gen_shirt_pants_color(type, rgen));
		break; // done
	} // for n
}

bool building_t::replace_light_with_ceiling_fan(rand_gen_t &rgen, cube_t const &room, cube_t const &avoid, unsigned room_id, float tot_light_amt, unsigned light_obj_ix) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_CEIL_FAN)) return 0;
	if (rgen.rand_float() > 0.3) return 0;
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(light_obj_ix < objs.size());
	room_object_t &light(objs[light_obj_ix]);
	if (light.type != TYPE_LIGHT) return 0; // error?
	float const floor_spacing(get_window_vspace());
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_CEIL_FAN)); // D, W, H
	float const diameter(min(0.4*min(room.dx(), room.dy()), 0.5*floor_spacing)), height(diameter*sz.z/sz.y); // assumes width = depth = diameter
	cube_t fan(cube_top_center(light)); // center on the light, with z2 on the ceiling
	fan.expand_by_xy(0.5*diameter);
	fan.z1() -= height;
	if (!avoid.is_all_zeros() && avoid.intersects(fan)) return 0; // check for closet intersection

	if (has_attic()) {
		cube_t aa_blocked(interior->attic_access);
		aa_blocked.expand_by(0.25*diameter);
		if (aa_blocked.intersects(fan)) return 0; // too close to attic access
	}
	light.translate_dim(2, -0.895*height); // move near the bottom of the ceiling fan (before invalidating with objs.emplace_back())
	light.flags |= RO_FLAG_INVIS;   // don't draw the light itself; assume the light is part of the bottom of the fan instead
	light.flags |= RO_FLAG_HANGING; // don't draw upward facing light
	unsigned flags(RO_FLAG_NOCOLL);
	if (rgen.rand_float() < 0.65) {flags |= RO_FLAG_ROTATING;} // make fan rotate when turned on 65% of the time
	objs.emplace_back(fan, TYPE_CEIL_FAN, room_id, 0, 0, flags, tot_light_amt, SHAPE_CYLIN, WHITE);
	objs.back().obj_id = light_obj_ix; // store light index in this object
	return 1;
}

colorRGBA get_bed_sheet_color(int tid, rand_gen_t &rgen) {
	unsigned const NUM_COLORS = 8;
	colorRGBA const colors[NUM_COLORS] = {WHITE, WHITE, WHITE, LT_BLUE, LT_BLUE, PINK, PINK, LT_GREEN}; // color of the sheets
	if (tid < 0 || tid == WHITE_TEX || texture_color(tid).get_luminance() > 0.5) {return colors[rgen.rand()%NUM_COLORS];}
	return WHITE; // textured
}

// Note: must be first placed object
bool building_t::add_bed_to_room(rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, float zval,
	unsigned room_id, float tot_light_amt, unsigned floor, bool force, int &bed_size_ix, room_object_t const &other_bed)
{
	cube_t room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace()), wall_thick(get_wall_thickness());
	bool const long_dim(room_bounds.dx() < room_bounds.dy()), is_store(room.is_store());
	bool const dim((is_hotel() && room.get_sz_dim(!long_dim) > 0.8*vspace) ? !long_dim : long_dim); // bed in long room dim, unless in a large hotel room (with 2 beds)
	vector3d expand, bed_sz;
	expand[ dim] = -wall_thick; // small amount of space
	expand[!dim] = -0.3f*vspace; // leave at least some space between the bed and the wall
	room_bounds.expand_by_xy(expand);
	float const room_len(room_bounds.get_sz_dim(dim)), room_width(room_bounds.get_sz_dim(!dim));
	
	if (force || is_store) { // no room is too large
		if (room_len < 1.0*vspace || room_width < 0.55*vspace) return 0; // room is too small to fit a bed
	}
	else if (floor == 0) { // special case for ground floor
		if (room_len < 1.3*vspace || room_width < 0.7*vspace) return 0; // room is too small to fit a bed
		if (room_len > 4.0*vspace || room_width > 2.5*vspace) return 0; // room is too large to be a bedroom
	}
	else { // more relaxed constraints
		if (room_len < 1.1*vspace || room_width < 0.6*vspace) return 0; // room is too small to fit a bed
		if (room_len > 4.5*vspace || room_width > 3.5*vspace) return 0; // room is too large to be a bedroom
	}
	bool first_head_dir(0);
	// place hotel room bed head by an exterior wall so that both beds can be placed there without blocking a door
	if (is_hotel()) {first_head_dir = (room.get_center_dim(dim) < get_part_for_room(room).get_center_dim(dim));}
	else {first_head_dir = rgen.rand_bool();}
	bool const first_wall_dir(rgen.rand_bool()), have_other_bed(other_bed.type == TYPE_BED);
	door_path_checker_t door_path_checker;
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	c.z1() = zval;

	for (unsigned n = 0; n < (force ? 100U : 40U); ++n) { // make 40 attempts to place a bed (100 if forced)
		float const sizes[6][2] = {{38, 75}, {38, 80}, {53, 75}, {60, 80}, {76, 80}, {72, 84}}; // twin, twin XL, full, queen, king, cal king
		unsigned size_ix(0);
		if (bed_size_ix >= 0) {size_ix = bed_size_ix;} // size was selected by the caller (to match pairs of beds in the same room)
		else {size_ix = ((room_width < 0.9*vspace) ? (rgen.rand() % 6) : (2 + (rgen.rand() % 4)));} // select a random size; only add twin beds to narrow rooms
		bed_sz[ dim] = 0.01f*vspace*(sizes[size_ix][1] + 8.0f); // length (mattress + headboard + footboard)
		bed_sz[!dim] = 0.01f*vspace*(sizes[size_ix][0] + 4.0f); // width  (mattress + small gaps)
		if (room_bounds.dx() < 1.5*bed_sz.x || room_bounds.dy() < 1.5*bed_sz.y) continue; // room is too small for a bed of this size
		bool const always_against_wall(room_bounds.dx() < 2.0*bed_sz.x || room_bounds.dy() < 2.0*bed_sz.y); // if room is narrow
		bool against_wall(0), wall_dir(0);
		bed_sz.z = (have_other_bed ? other_bed.dz() : 0.3*vspace*rgen.rand_uniform(1.0, 1.2)); // height
		c.z2()   = zval + bed_sz.z;

		for (unsigned d = 0; d < 2; ++d) { // dim
			float const min_val(room_bounds.d[d][0]), max_val(room_bounds.d[d][1] - bed_sz[d]);

			if (bool(d) == dim && n < 5) { // in the first few iterations, try to place the head of the bed against the wall (maybe not for exterior wall facing window?)
				// if there was another bed, use the same wall/dir (head_dir == !bed_dir); otherwise, alternate head_dir
				bool const head_dir((have_other_bed && other_bed.dim == dim) ? !other_bed.dir : (first_head_dir ^ bool(n&1)));
				c.d[d][0] = (head_dir ? min_val : max_val);
			}
			else if (bool(d) != dim && (always_against_wall || rgen.rand_bool())) { // try to place the bed against the wall sometimes
				wall_dir  = (first_wall_dir ^ bool(n&1) ^ 1);
				c.d[d][0] = (wall_dir ? (max_val + 0.25*vspace) : (min_val - 0.25*vspace));
				against_wall = 1;
			}
			else {
				c.d[d][0] = rgen.rand_uniform(min_val, max_val);
			}
			c.d[d][1] = c.d[d][0] + bed_sz[d];
		} // for d
		if (!is_valid_placement_for_room(c, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
		// prefer not to block the path between doors in the first half of iterations
		if (n < 10 && door_path_checker.check_door_path_blocked(c, room, room_id, zval, *this)) continue;
		bool const dir((room_bounds.d[dim][1] - c.d[dim][1]) < (c.d[dim][0] - room_bounds.d[dim][0])); // head of the bed is closer to the wall; opposite first_head_dir
		unsigned const bed_ix(objs.size());
		objs.emplace_back(c, TYPE_BED, room_id, dim, dir, (is_house ? RO_FLAG_IS_HOUSE : 0), tot_light_amt);
		set_obj_id(objs);
		room_object_t &bed(objs.back());
		// use white color if a texture is assigned that's not close to white
		int const sheet_tid(bed.get_sheet_tid());
		bed.color = get_bed_sheet_color(sheet_tid, rgen);
		unsigned num_beds(1);

		if (is_house && against_wall && size_ix <= 2 && rgen.rand_float() < 0.4) { // maybe add a bunk bed if against a wall and <= twin XL size
			room_object_t bed2(bed);
			bed2.translate_dim(2, c.dz());

			if (is_valid_placement_for_room(bed2, room, blockers, 1)) { // always true?
				bed .flags |= RO_FLAG_ADJ_BOT; // flag as bottom bunk
				bed2.flags |= RO_FLAG_ADJ_TOP; // flag as top    bunk
				objs.push_back(bed2); // same options and style; invalidates bed reference
				++num_beds;
				// add a small ladder; hopefully it shouldn't intersect anything since it's on the inside edge
				float const bed_edge(bed2.d[!dim][!wall_dir]);
				cube_t ladder(bed2);
				ladder.z1()  = zval; // down to the floor
				ladder.z2() -= 0.25*bed2.dz(); // lower to just above mattress level
				ladder.d[!dim][ wall_dir] = bed_edge;
				ladder.d[!dim][!wall_dir] = bed_edge + (wall_dir ? -1.0 : 1.0)*0.75*wall_thick; // set depth
				ladder.d[ dim][ dir]  = bed2.d[dim][!dir] + (dir ? 1.0 : -1.0)*0.22*vspace;
				ladder.d[ dim][!dir] += (dir ? 1.0 : -1.0)*0.25*wall_thick; // move away from footboard
				objs.emplace_back(ladder, TYPE_INT_LADDER, room_id, !dim, wall_dir, (RO_FLAG_IS_HOUSE | RO_FLAG_HANGING), tot_light_amt);
			}
		}
		for (unsigned bix = 0; bix < num_beds; ++bix) {
			room_object_t const &bed(objs[bed_ix + bix]);
			cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
			get_bed_cubes(bed, cubes);
			cube_t const &mattress(cubes[3]), &pillow(cubes[4]);
			float const rand_val(rgen.rand_float());

			if (rand_val < 0.4) { // add a blanket on the bed 40% of the time
				vector3d const mattress_sz(mattress.get_size());
				cube_t blanket(mattress);
				set_cube_zvals(blanket, mattress.z2(), (mattress.z2() + 0.02*mattress_sz.z)); // on top of mattress; set height
				blanket.d[dim][ dir] = pillow.d[dim][!dir] - (dir ? 1.0 : -1.0)*rgen.rand_uniform(0.01, 0.06)*mattress_sz[dim]; // shrink at head
				blanket.d[dim][!dir] += (dir ? 1.0 : -1.0)*rgen.rand_uniform(0.03, 0.08)*mattress_sz[dim]; // shrink at foot
				blanket.expand_in_dim(!dim, -rgen.rand_uniform(0.08, 0.16)*mattress_sz[!dim]); // shrink width
				objs.emplace_back(blanket, TYPE_BLANKET, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt);
				set_obj_id(objs);
			}
			else if (!is_store && rand_val < 0.7) { // add teeshirt or jeans on the bed 30% of the time
				unsigned const type(rgen.rand_bool() ? TYPE_PANTS : TYPE_TEESHIRT);
				float const length(((type == TYPE_TEESHIRT) ? 0.26 : 0.2)*vspace), width(0.98*length), height(0.002*vspace);
				float const clearance(get_min_front_clearance_inc_people());
				cube_t valid_area(mattress);
				valid_area.d[dim][dir] = pillow.d[dim][!dir]; // don't place under the pillow
				vector3d size(0.5*length, 0.5*width, height);
			
				if (valid_area.get_sz_dim(!dim) > 4.0*size[!dim]) { // don't place on the side near a wall if the bed is wide
					if (bed .d[!dim][0] - room.d[!dim][0] < clearance) {valid_area.d[!dim][0] = bed.get_center_dim(!dim);}
					if (room.d[!dim][1] - bed .d[!dim][1] < clearance) {valid_area.d[!dim][1] = bed.get_center_dim(!dim);}
				}
				bool const dim2(rgen.rand_bool()), dir2(rgen.rand_bool()); // choose a random orientation
				if (dim2) {std::swap(size.x, size.y);}

				if (valid_area.dx() > 2.0*size.x && valid_area.dy() > 2.0*size.y) {
					cube_t c(gen_xy_pos_in_area(valid_area, size, rgen, mattress.z2()));
					c.expand_by_xy(size);
					c.z2() += size.z;
					objs.emplace_back(c, type, room_id, dim2, dir2, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, gen_shirt_pants_color(type, rgen));
				}
			}
		} // for bix
		bed_size_ix = size_ix;
		return 1; // done/success
	} // for n
	return 0;
}

bool building_t::add_ball_to_room(rand_gen_t &rgen, room_t const &room, cube_t const &place_area, float zval,
	unsigned room_id, float tot_light_amt, unsigned objs_start, int force_type, cube_t const &avoid_xy, bool in_pool)
{
	unsigned const btype((force_type >= 0) ? force_type : (rgen.rand() % NUM_BALL_TYPES));
	ball_type_t const &bt(ball_types[btype]);
	float const floor_spacing(get_window_vspace()), radius(bt.radius*floor_spacing/(12*8)); // assumes 8 foot floor spacing, and bt.radius in inches
	cube_t ball_area(place_area);
	ball_area.expand_by_xy(-rgen.rand_uniform(radius, max(2.0f*radius, min(10.0f*radius, 0.5f*floor_spacing))));
	vect_room_object_t &objs(interior->room_geom->objs);
	if (!ball_area.is_strictly_normalized()) return 0; // should always be normalized
	bool const is_backrooms(!is_house && has_parking_garage && room.is_ext_basement());
	if (is_backrooms && min(ball_area.dx(), ball_area.dy()) < 4.0*radius) return 0; // too small; shouldn't happen
	float const ceil_zval(zval + get_floor_ceil_gap()); // zval should be valid whether or not in the pool
	
	for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to place the object
		point center(0.0, 0.0, zval);

		if (is_backrooms) { // office backrooms: place anywhere within the room
			center = gen_xy_pos_in_area(ball_area, radius, rgen, center.z);
		}
		else { // house bedroom or swimming pool room: place along a wall
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
			center[ dim] = ball_area.d[dim][dir];
			center[!dim] = rgen.rand_uniform(ball_area.d[!dim][0], ball_area.d[!dim][1]); // random position along the wall
		}
		if (in_pool) {
			if (has_water()) {center.z = (interior->water_zval - 2.0*bt.density*bt.radius*floor_spacing/(12*8));} // floats on the water
			else {get_zval_for_pool_bottom(center, center.z);} // on the bottom of the pool
		}
		center.z += radius;
		set_float_height(center, radius, ceil_zval, bt.density); // maybe floats on water (in backrooms)
		cube_t c(center);
		c.expand_by(radius);
		if (!avoid_xy.is_all_zeros() && c.intersects_xy(avoid_xy)) continue;
		if (overlaps_other_room_obj(c, objs_start) || is_obj_placement_blocked(c, room, 1)) continue; // bad placement
		objs.emplace_back(c, TYPE_LG_BALL, room_id, 0, 0, RO_FLAG_DSTATE, tot_light_amt, SHAPE_SPHERE, WHITE, btype);
		objs.back().obj_id = (uint16_t)interior->room_geom->allocate_dynamic_state(); // allocate a new dynamic state object
		return 1; // done
	} // for n
	return 0;
}

colorRGBA gen_vase_color(rand_gen_t &rgen) { // or urn
	if (rgen.rand_bool()) return WHITE; // will be textured
	return colorRGBA(rgen.rand_float(), rgen.rand_float(), rgen.rand_float(), 1.0); // solid pastel color
}

// Note: modified blockers rather than using it; fireplace must be the first placed object
bool building_t::maybe_add_fireplace_to_room(rand_gen_t &rgen, room_t const &room, vect_cube_t &blockers, float zval, unsigned room_id, float tot_light_amt) {
	// Note: the first part of the code below is run on every first floor room and will duplicate work, so it may be better to factor it out somehow
	cube_t fireplace(get_fireplace()); // make a copy of the exterior fireplace that will be converted to an interior fireplace
	bool dim(0), dir(0);
	if      (fireplace.x1() <= bcube.x1()) {dim = 0; dir = 0;} // Note: may not work on rotated buildings
	else if (fireplace.x2() >= bcube.x2()) {dim = 0; dir = 1;}
	else if (fireplace.y1() <= bcube.y1()) {dim = 1; dir = 0;}
	else if (fireplace.y2() >= bcube.y2()) {dim = 1; dir = 1;}
	else {assert(is_rotated()); return 0;} // can fail on rotated buildings?
	float const depth_signed((dir ? -1.0 : 1.0)*1.0*fireplace.get_sz_dim(dim)), wall_pos(fireplace.d[dim][!dir]), top_gap(0.15*fireplace.dz());
	fireplace.d[dim][ dir] = wall_pos; // flush with the house wall
	fireplace.d[dim][!dir] = wall_pos + depth_signed; // extend out into the room
	fireplace.z2() -= top_gap; // shorten slightly
	if (!room.contains_cube_xy_exp(fireplace, -0.5*get_wall_thickness())) return 0; // fireplace not in this room; allow fireplace to extend slightly into room walls
	// the code below should be run at most once per building
	cube_t fireplace_ext(fireplace);
	fireplace_ext.d[dim][!dir] = fireplace.d[dim][!dir] + 0.4*depth_signed; // extend out into the room for clearance
	if (interior->is_blocked_by_stairs_or_elevator(fireplace_ext)) return 0; // blocked by stairs, don't add (would be more correct to relocate stairs) - should no longer fail
	if (is_cube_close_to_doorway(fireplace_ext, room, 0.0, 1))     return 0; // too close to door
	fireplace.d[dim][dir] = room.d[dim][dir]; // re-align to room to remove any gap between the fireplace and the exterior wall
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(fireplace, TYPE_FPLACE, room_id, dim, dir, 0, tot_light_amt);
	fireplace_ext.d[dim][!dir] += 0.1*depth_signed; // extend out into the room even further for clearance
	cube_t blocker(fireplace_ext);
	blocker.d[dim][dir] = fireplace.d[dim][!dir]; // flush with the front of the fireplace
	objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, dir, RO_FLAG_INVIS);
	blockers.push_back(fireplace_ext); // add as a blocker if it's not already there

	if (rgen.rand_bool()) { // add urn on fireplace
		float const urn_height(rgen.rand_uniform(0.65, 0.95)*top_gap), urn_radius(rgen.rand_uniform(0.35, 0.45)*min(urn_height, fabs(depth_signed)));
		point center(fireplace.get_cube_center());
		center[!dim] += 0.45*fireplace.get_sz_dim(!dim)*rgen.signed_rand_float(); // random placement to the left and right of center
		cube_t urn;
		urn.set_from_sphere(center, urn_radius);
		set_cube_zvals(urn, fireplace.z2(), fireplace.z2()+urn_height); // place on the top of the fireplace
		objs.emplace_back(urn, TYPE_URN, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, gen_vase_color(rgen));
		set_obj_id(objs);
	}
	has_int_fplace = 1;
	return 1;
}

bool building_t::check_if_against_window(cube_t const &c, room_t const &room, bool dim, bool dir) const {
	if (!has_int_windows() || classify_room_wall(room, c.zc(), dim, dir, 0) != ROOM_WALL_EXT) return 0;
	cube_t const &part(get_part_for_room(room));
	float const hspacing(get_hspacing_for_part(part, !dim)), border(get_window_h_border());
	// assume object is no larger than 2x window size and check left, right, and center positions
	return (is_val_inside_window(part, !dim, c.d[!dim][0], hspacing, border) ||
		    is_val_inside_window(part, !dim, c.d[!dim][1], hspacing, border) ||
		    is_val_inside_window(part, !dim, c.get_center_dim(!dim), hspacing, border));
}
// Note: front_clearance is relative to depth, while side_clearance is an absolute distance
bool building_t::place_obj_along_wall(room_object type, room_t const &room, float height, vector3d const &sz_scale, rand_gen_t &rgen, float zval,
	unsigned room_id, float tot_light_amt, cube_t const &place_area, unsigned objs_start, float front_clearance, bool add_door_clearance, unsigned pref_orient,
	bool pref_centered, colorRGBA const &color, bool not_at_window, room_obj_shape shape, float side_clearance, unsigned extra_flags, bool not_ext_wall, bool force_pref)
{
	float const hwidth(0.5*height*sz_scale.y/sz_scale.z), depth(height*sz_scale.x/sz_scale.z);
	float const min_space(max(2.8f*hwidth, 2.1f*(max(hwidth, 0.5f*depth) + get_scaled_player_radius()))); // make sure the player can get around the object
	vector2d const place_area_sz(place_area.get_size_xy());
	if (max(place_area_sz.x, place_area_sz.y) <= min_space) return 0; // can't fit in either dim
	unsigned const force_dim((place_area_sz.x <= min_space) ? 0 : ((place_area_sz.y <= min_space) ? 1 : 2)); // *other* dim; 2=neither
	float const obj_clearance(depth*front_clearance), clearance(max(obj_clearance, get_min_front_clearance_inc_people()));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	set_cube_zvals(c, zval, zval+height);
	bool center_tried[4] = {};

	for (unsigned n = 0; n < 25; ++n) { // make 25 attempts to place the object
		bool const use_pref(pref_orient < 4 && (force_pref || n < 10)); // use pref orient for first 10 tries unless force_pref=1
		bool const dim((force_dim < 2) ? force_dim : (use_pref ? (pref_orient >> 1) : rgen.rand_bool())); // choose a random wall unless forced
		bool const dir(use_pref ? !(pref_orient & 1) : rgen.rand_bool()); // dir is inverted for the model, so we invert pref dir as well
		unsigned const orient(2*dim + dir);
		float center(0.0);
		if (pref_centered && !center_tried[orient]) {center = place_area.get_center_dim(!dim); center_tried[orient] = 1;} // try centered
		else {center = rgen.rand_uniform(place_area.d[!dim][0]+hwidth, place_area.d[!dim][1]-hwidth);} // random position
		c.d[dim][ dir] = place_area.d[dim][dir];
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		set_wall_width(c, center, hwidth, !dim);
		if (!room.contains_cube(c)) continue; // larger than room width?
		if (not_ext_wall  && classify_room_wall(room, c.zc(), dim, dir, 0) == ROOM_WALL_EXT) continue;
		if (not_at_window && check_if_against_window(c, room, dim, dir)) continue;
		if (not_at_window && overlaps_or_adj_int_window(c))              continue; // check interior windows as well, in case this is a conference room
		cube_t c2(c), c3(c); // used for collision tests
		c2.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance; // add front clearance
		if (!room.contains_cube(c2)) continue; // not enough clearance
		cube_t c2b(c2);
		c2b.expand_in_dim(!dim, side_clearance); // side_clearance applies to other objects, stairs, and elevators; used with boxes, medicine cabinets, etc.
		if (overlaps_other_room_obj(c2b, objs_start) || interior->is_blocked_by_stairs_or_elevator(c2b)) continue; // bad placement (Note: not using is_obj_placement_blocked())
		c3.d[dim][!dir] += (dir ? -1.0 : 1.0)*obj_clearance; // smaller clearance value (without player diameter)
		c3.expand_in_dim(!dim, side_clearance);
		if (!place_area.contains_cube_xy(c3)) continue; // not enough clearance in place area

		if (add_door_clearance) {
			if (is_cube_close_to_doorway(c3, room, 0.0, 1)) continue; // bad placement
		}
		else { // we don't need clearance for both door and object; test the object itself against the open door and the object with clearance against the closed door
			if (is_cube_close_to_doorway(c,  room, 0.0, 1)) continue; // bad placement
			if (is_cube_close_to_doorway(c3, room, 0.0, 0)) continue; // bad placement
		}
		if (!check_cube_within_part_sides(c )) continue; // handle non-cube buildings
		if (is_blocked_by_open_attic_door(c3)) continue; // blocked by attic ladder when open
		unsigned flags(extra_flags);
		if (type == TYPE_BOX) {flags |= (RO_FLAG_ADJ_LO << orient);} // set wall edge bit for boxes (what about other dim bit if place in room corner?)
		objs.emplace_back(c, type, room_id, dim, !dir, flags, tot_light_amt, shape, color);
		//if (type == TYPE_BAR_STOOL) {objs.emplace_back(c, TYPE_DBG_SHAPE, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, RED);} // TESTING
		
		if (type == TYPE_TOILET || type == TYPE_SINK || type == TYPE_URINAL || type == TYPE_TUB || type == TYPE_SHOWERTUB) {
			if (type == TYPE_SINK && check_if_against_window(c, room, dim, dir)) {} // no plumbing into window
			else {add_bathroom_plumbing(objs.back());}
		}
		set_obj_id(objs);
		if (front_clearance > 0.0) {objs.emplace_back(c2, TYPE_BLOCKER, room_id, dim, !dir, RO_FLAG_INVIS);} // add blocker cube to ensure no other object overlaps this space
		return 1; // done
	} // for n
	return 0; // failed
}
bool building_t::place_model_along_wall(unsigned model_id, room_object type, room_t const &room, float height, rand_gen_t &rgen, float zval, unsigned room_id,
	float tot_light_amt, cube_t const &place_area, unsigned objs_start, float front_clearance, unsigned pref_orient, bool pref_centered, colorRGBA const &color,
	bool not_at_window, unsigned extra_flags, bool force_pref, bool sideways)
{
	if (place_area.is_all_zeros()) return 0;
	if (!building_obj_model_loader.is_model_valid(model_id)) return 0; // don't have a model of this type
	vect_room_object_t &objs(interior->room_geom->objs);
	vector3d sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
	if (sideways) {swap(sz[0], sz[1]);}
	unsigned const obj_id(objs.size());
	if (!place_obj_along_wall(type, room, height*get_window_vspace(), sz, rgen, zval, room_id, tot_light_amt, place_area,
		objs_start, front_clearance, 0, pref_orient, pref_centered, color, not_at_window, SHAPE_CUBE, 0.0, extra_flags, 0, force_pref)) return 0;
	if (sideways) {objs[obj_id].dim ^= 1;}
	return 1;
}

float building_t::add_flooring(room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned flooring_type) {
	float const new_zval(zval + get_flooring_thick());
	cube_t flooring(get_walkable_room_bounds(room));
	// expand flooring to include half of the walls so that it meets the door and trim is added; not for garages and sheds
	if (!room.is_sec_bldg) {flooring.expand_by_xy(0.5*get_wall_thickness());}

	if (!room.is_ext_basement() && room.part_id < real_num_parts) { // skip if ext basement, including backrooms bathrooms
		if (parts[room.part_id].contains_cube(room)) {flooring.intersect_with_cube_xy(parts[room.part_id]);}
		assert(flooring.is_strictly_normalized());
	}
	set_cube_zvals(flooring, zval, new_zval);
	tot_light_amt = 0.5*tot_light_amt + 0.5; // brighten flooring so that lights shining through doors and flashlights look better
	unsigned flags(RO_FLAG_NOCOLL);
	if (room.open_wall_mask) {flags |= RO_FLAG_OPEN;} // flag flooring as "open" so that color is not adjusted by room light
	// cut elevators out of flooring in case they pass through bathrooms or utility rooms
	static vect_cube_t fparts, temp;
	subtract_cubes_from_cube(flooring, interior->elevators, fparts, temp, 2); // check zval overlap
	for (cube_t const &f : fparts) {interior->room_geom->objs.emplace_back(f, TYPE_FLOORING, room_id, 0, 0, flags, tot_light_amt, SHAPE_CUBE, WHITE, flooring_type);}
	return new_zval;
}

bool building_t::add_bathroom_objs(rand_gen_t rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt,
	unsigned lights_start, unsigned objs_start, unsigned floor, bool is_basement, bool add_shower_tub, unsigned &added_bathroom_objs_mask)
{
	// Note: zval passed by reference
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());

	if (!skylights.empty()) { // check for skylights; should we allow bathroom stalls (with ceilings?) even if there's a skylight?
		cube_t test_cube(room);
		set_cube_zvals(test_cube, zval, zval+floor_spacing);
		if (check_skylight_intersection(test_cube)) return 0;
	}
	cube_t const room_bounds(get_walkable_room_bounds(room));
	cube_t place_area(room_bounds);
	place_area.expand_by(-0.5*wall_thickness);
	vector2d const place_area_sz(place_area.dx(), place_area.dy());
	if (min(place_area_sz.x, place_area_sz.y) < 0.7*floor_spacing) return 0; // room is too small (should be rare)
	bool const have_toilet(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOILET)), have_sink(building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK));
	vect_room_object_t &objs(interior->room_geom->objs);

	if ((have_toilet || have_sink) && is_cube() && !is_industrial()) { // bathroom with at least a toilet or sink; cube shaped parts only; no industrial; add flooring
		int const flooring_type(is_residential() ? (is_basement ? (int)FLOORING_CONCRETE : (int)FLOORING_TILE) : (int)FLOORING_MARBLE);
		if (flooring_type == FLOORING_CONCRETE && get_material().basement_floor_tex.tid == get_concrete_tid()) {} // already concrete
		else { // replace carpet/wood with marble/tile/concrete
			zval = add_flooring(room, zval, room_id, tot_light_amt, flooring_type); // move the effective floor up
		}
	}
	if (have_toilet && room.is_office) { // office bathroom
		if (min(place_area_sz.x, place_area_sz.y) > 1.5*floor_spacing && max(place_area_sz.x, place_area_sz.y) > 2.0*floor_spacing) {
			if (divide_bathroom_into_stalls(rgen, room, zval, room_id, tot_light_amt, floor, lights_start, objs_start)) { // large enough, divide into bathroom stalls
				added_bathroom_objs_mask |= (PLACED_TOILET | PLACED_SINK);
				return 1;
			}
		}
	}
	float const tub_height_factor(0.2); // in units of floor spacing
	bool placed_obj(0), placed_toilet(0), no_tub(0), added_vanity(0);
	
	if (is_house && !is_basement && add_vanity_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start)) { // try to place vanity
		added_vanity = placed_obj = 1;
		added_bathroom_objs_mask |= PLACED_SINK;
	}
	// place toilet first because it's in the corner out of the way and higher priority
	if (have_toilet) { // have a toilet model
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
		float const height(0.35*floor_spacing), width(height*sz.y/sz.z), length(height*sz.x/sz.z); // for toilet
		unsigned const first_corner(rgen.rand() & 3);
		bool const first_dim(rgen.rand_bool());

		for (unsigned n = 0; n < 4 && !placed_toilet; ++n) { // try 4 room corners
			unsigned const corner_ix((first_corner + n)&3);
			bool const xdir(corner_ix&1), ydir(corner_ix>>1);
			point const corner(place_area.d[0][xdir], place_area.d[1][ydir], zval);
			if (!check_pt_within_part_sides(corner)) continue; // invalid corner

			for (unsigned d = 0; d < 2 && !placed_toilet; ++d) { // try both dims
				bool const dim(bool(d) ^ first_dim), dir(dim ? ydir : xdir);
				cube_t c(corner, corner);
				c.d[0][!xdir] += (xdir ? -1.0 : 1.0)*(dim ? width : length);
				c.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(dim ? length : width);
				for (unsigned e = 0; e < 2; ++e) {c.d[!dim][e] += ((dim ? xdir : ydir) ? -1.5 : 1.5)*wall_thickness;} // extra padding on left and right sides
				c.z2() += height;
				cube_t c2(c); // used for placement tests
				c2.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.8*length; // extra padding in front of toilet, to avoid placing other objects there (sink and tub)
				c2.expand_in_dim(!dim, 0.4*width); // more padding on the sides
				if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, room, 0.0, 1)) continue; // bad placement
				objs.emplace_back(c,  TYPE_TOILET,  room_id, dim, !dir, 0, tot_light_amt);
				add_bathroom_plumbing(objs.back());
				objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker cube to ensure no other object overlaps this space
				placed_obj = placed_toilet = 1; // done
				added_bathroom_objs_mask  |= PLACED_TOILET;

				// try to place a roll of toilet paper on the adjacent wall
				bool const tp_dir(dim ? xdir : ydir);
				float const length(0.18*height), wall_pos(c.get_center_dim(dim)), far_edge_pos(wall_pos + (dir ? -1.0 : 1.0)*0.5*length);
				cube_t const &part(get_part_for_room(room));

				// if this wall has windows and bathroom has multiple exterior walls (which means it has non-glass block windows), don't place a TP roll
				if (is_basement || !has_int_windows() || classify_room_wall(room, zval, !dim, tp_dir, 0) != ROOM_WALL_EXT ||
					!is_val_inside_window(part, dim, far_edge_pos, get_hspacing_for_part(part, dim), get_window_h_border()) || count_ext_walls_for_room(room, zval) <= 1)
				{
					add_tp_roll(room_bounds, room_id, tot_light_amt, !dim, tp_dir, length, (c.z1() + 0.7*height), wall_pos);
				}
			} // for d
		} // for n
		if (!placed_toilet) { // if the toilet can't be placed in a corner, allow it to be placed anywhere; needed for small offices
			placed_toilet = place_model_along_wall(OBJ_MODEL_TOILET, TYPE_TOILET, room, 0.35, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8);
			placed_obj   |= placed_toilet;
			added_bathroom_objs_mask |= PLACED_TOILET;

			if (placed_toilet) { // if toilet was placed, try to place a roll of toilet paper on the same wall as the toilet
				room_object_t const &toilet(objs.back()); // okay if this is the blocker
				
				// Note: not calling is_val_inside_window() here because I don't have a test case for that and it may not even be possible to get here when the toilet is next to a window
				if (is_basement || !has_int_windows() || classify_room_wall(room, zval, toilet.dim, !toilet.dir, 0) != ROOM_WALL_EXT) { // check for possible windows
					bool place_dir(rgen.rand_bool()); // pick a random starting side

					for (unsigned d = 0; d < 2; ++d) {
						float const length(0.18*height), wall_pos(toilet.d[!toilet.dim][place_dir] + (place_dir ? 1.0 : -1.0)*0.5*width);
						if (add_tp_roll(room_bounds, room_id, tot_light_amt, toilet.dim, !toilet.dir, length, (toilet.z1() + 0.7*height), wall_pos, 1)) break; // check_valid=1
						place_dir ^= 1; // try the other dir
					} // for d
				}
			}
		}
	}
	unsigned const pre_shower_tub_sink_ix(objs.size());

	for (unsigned n = 0; n < 20; ++n) { // 20 tries to add a tub/shower and sink
		unsigned bathroom_objs_mask(0);

		// try to add a shower; 50% chance if on first floor of a house; not in basements (due to drawing artifacts)
		if (add_shower_tub && (is_apt_or_hotel() || floor > 0 || rgen.rand_bool())) {
			float shower_dx(0.0), shower_dy(0.0), wall_thick(0.0);
			bool hdim(0); // handle dim/long dim
			// use a shower + tub combo for basements due to glass drawing artifacts, and for arpartments and hotels with small bathrooms; requires tub model
			bool const place_shower_tub((is_basement || is_apt_or_hotel() || rgen.rand_float() < 0.25) && building_obj_model_loader.is_model_valid(OBJ_MODEL_TUB));
			float const shower_height(place_shower_tub ? get_floor_ceil_gap() : 0.8*floor_spacing), tub_height(tub_height_factor*floor_spacing);

			if (place_shower_tub) {
				vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TUB)); // D, W, H
				hdim       = rgen.rand_bool();
				shower_dx  = tub_height*(sz.y/sz.z); // width
				shower_dy  = tub_height*(sz.x/sz.z); // depth
				wall_thick = 0.05*shower_dx;
				shower_dx += wall_thick;
				if (hdim) {swap(shower_dx, shower_dy);}
			}
			else {
				for (unsigned d = 0; d < 2; ++d) {(d ? shower_dy : shower_dx) = rgen.rand_uniform(0.4, 0.5)*floor_spacing;}
				hdim = (shower_dx < shower_dy); // larger dim, must match handle/door drawing code
			}
			unsigned const first_corner(rgen.rand() & 3);
			bool placed_shower(0), is_ext_wall[2][2] = {0};
		
			if (!is_basement && has_int_windows()) { // precompute which walls are exterior, {dim}x{dir}; basement walls are not considered exterior because there are no windows
				for (unsigned d = 0; d < 4; ++d) {is_ext_wall[d>>1][d&1] = (classify_room_wall(room, zval, (d>>1), (d&1), 0) == ROOM_WALL_EXT);}
			}
			for (unsigned ar = 0; ar < 2; ++ar) { // try both aspect ratios/door sides
				for (unsigned n = 0; n < 4; ++n) { // try 4 room corners
					unsigned const corner_ix((first_corner + n)&3);
					bool const xdir(corner_ix&1), ydir(corner_ix>>1), dirs[2] = {xdir, ydir};
					point const corner(room_bounds.d[0][xdir], room_bounds.d[1][ydir], zval); // flush against the wall
					cube_t c(corner, corner);
					c.d[0][!xdir] += (xdir ? -1.0 : 1.0)*shower_dx;
					c.d[1][!ydir] += (ydir ? -1.0 : 1.0)*shower_dy;
					c.z2() += shower_height; // set height

					if (!place_shower_tub || has_int_windows()) {
						// exterior walls aren't drawn in the correct order for shower glass alpha blend, so skip any exterior walls;
						// also don't place shower tubs near windows; assume they're long enough to always overlap a window when placed on an exterior wall
						bool is_bad(0);

						for (unsigned d = 0; d < 2; ++d) { // check for window intersection
							if (is_ext_wall[!d][dirs[!d]]) {is_bad = 1; break;}
						}
						if (is_bad) continue;
					}
					cube_t c2(c); // used for placement tests

					if (place_shower_tub) {
						// anything?
					}
					else { // shower: extend out by door width on the side that opens, and a small amount on the other side
						c2.d[0][!xdir] += (xdir ? -1.0 : 1.0)*((!hdim) ? 1.1*shower_dy : 0.2*shower_dx);
						c2.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(  hdim  ? 1.1*shower_dx : 0.2*shower_dy);
					}
					if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, room, 0.0, 1)) continue; // bad placement
					bool const dim(place_shower_tub ? !hdim : xdir), dir(place_shower_tub ? !(hdim ? xdir : ydir) : ydir); // different encodings
					objs.emplace_back(c, (place_shower_tub ? TYPE_SHOWERTUB : TYPE_SHOWER), room_id, dim, dir, 0, tot_light_amt);
					set_obj_id(objs); // selects tile texture/color

					if (place_shower_tub) { // add the tub part as well
						bool const wall_dir(hdim ? ydir : xdir);
						// set flag to indicate which side is the wall for adding the shower head, and make open by default
						unsigned const lo_hi_flag(wall_dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
						objs.back().flags |= lo_hi_flag | RO_FLAG_OPEN;
						objs.back().color  = (is_basement ? WHITE : wall_color); // color of the end wall
						cube_t tub(c);
						tub.z2() = c.z1() + tub_height;
						tub.d[!dim][!wall_dir] -= (wall_dir ? -1.0 : 1.0)*wall_thick; // shrink off the wall
						tub.translate_dim(dim, (dir ? 1.0 : -1.0)*0.1*get_trim_thickness()); // shift away from the wall slightly to prevent Z-fighting
						objs.emplace_back(tub, TYPE_TUB, room_id, dim, dir, lo_hi_flag, tot_light_amt);
						bathroom_objs_mask |= PLACED_TUB;
						no_tub = 1;
					}
					objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker cube to ensure no other object overlaps this space
					placed_obj = placed_shower = 1;
					bathroom_objs_mask  |= PLACED_SHOWER;
					break; // done
				} // for n
				if (placed_shower) break; // done
				swap(shower_dx, shower_dy); // try the other aspect ratio
				hdim ^= 1;
			} // for ar
		}
		// place a tub, but not in office buildings; placed before the sink because it's the largest and the most limited in valid locations
		if (!no_tub && add_shower_tub && (!is_basement || rgen.rand_bool())) { // 50% of the time if in the basement
			cube_t place_area_tub(room_bounds);
			place_area_tub.expand_by(-get_trim_thickness()); // just enough to prevent z-fighting and intersecting the wall trim
		
			if (place_model_along_wall(OBJ_MODEL_TUB, TYPE_TUB, room, tub_height_factor, rgen, zval, room_id, tot_light_amt, place_area_tub, objs_start, 0.4)) {
				placed_obj = 1;
				bathroom_objs_mask |= PLACED_TUB;
			}
		}
		unsigned const sink_obj_ix(objs.size());

		if (added_vanity) {} // added vanity, no need to place a sink
		else if (place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6)) {
			placed_obj = 1;
			bathroom_objs_mask |= PLACED_SINK;
			assert(sink_obj_ix < objs.size());
			room_object_t const &sink(objs[sink_obj_ix]); // sink, not blocker
		
			if (point_in_water_area(sink.get_cube_center())) {} // no medicine cabinet, because the reflection system doesn't support both a mirror and water reflection
			else if (is_parking() && !is_basement) {} // no mirror in parking structure since reflections don't work there
			else if (is_basement || classify_room_wall(room, zval, sink.dim, !sink.dir, 0) != ROOM_WALL_EXT) { // interior wall only
				// add a mirror/medicine cabinet above the sink
				float const mirror_expand(0.1*sink.get_sz_dim(!sink.dim));
				cube_t mirror(sink); // start with the sink left and right position
				mirror.expand_in_dim(!sink.dim, mirror_expand); // make slightly wider
				set_cube_zvals(mirror, sink.z2(), sink.z2()+0.3*floor_spacing);
				mirror.d[sink.dim][!sink.dir] = room_bounds.d[sink.dim][!sink.dir];
				mirror.d[sink.dim][ sink.dir] = mirror.d[sink.dim][!sink.dir] + (sink.dir ? 1.0 : -1.0)*1.0*wall_thickness; // thickness

				if (!overlaps_other_room_obj(mirror, objs_start, 0, &sink_obj_ix)) { // check_all=0; skip sink + blocker
					if (interior->is_cube_close_to_doorway(mirror, room, 0.0, 1, 1)) { // check doorways as well, since mirror is wider thank sink
						mirror.expand_in_dim(!sink.dim, -mirror_expand); // undo expand and try for a narrow mirror
					}
					// this mirror is actually 3D, so we enable collision detection; treat as a house even if it's in an office building
					unsigned flags(RO_FLAG_IS_HOUSE); // Note: not necessarily a house
					if (count_ext_walls_for_room(room, mirror.z1()) == 1) {flags |= RO_FLAG_INTERIOR;} // flag as interior if windows are opaque glass blocks
					objs.emplace_back(mirror, TYPE_MED_CAB, room_id, sink.dim, sink.dir, flags, tot_light_amt);
					set_obj_id(objs); // for crack texture selection/orient
					room.set_has_mirror();
				}
			}
		}
		else if (n+1 < 20) { // not the last try, rip up and retry
			objs.resize(pre_shower_tub_sink_ix);
			continue;
		}
		added_bathroom_objs_mask |= bathroom_objs_mask;
		break;
	} // for n
	if (room.is_office && !room.is_nested()) {add_door_sign("Restroom", room, zval, room_id);} // add office bathroom sign; not for hospital rooms
	return placed_obj;
}

bool building_t::add_vanity_to_room(rand_gen_t &rgen, room_t &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t const room_bounds(get_walkable_room_bounds(room));
	vector2d const room_sz(room_bounds.dx(), room_bounds.dy());
	vect_room_object_t &objs(interior->room_geom->objs);
	// select an interior wall that doesn't have a door
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors; there should be no exterior doors
	vector<pair<float, unsigned>> avail_walls;
	unsigned num_ext_walls(0);

	for (unsigned d = 0; d < 4; ++d) {
		bool const dim(d >> 1), dir(d & 1);
		if (classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) {++num_ext_walls; continue;} // skip exterior wall
		cube_t wall(room);
		set_wall_width(wall, room.d[dim][dir], wall_thickness, dim);
		if (has_bcube_int(wall, doorways)) continue;
		avail_walls.emplace_back(rgen.rand_uniform(1.0, 1.1)*room_sz[dim], d); // prefer long wall, with some randomness
	} // for d
	sort(avail_walls.begin(), avail_walls.end());

	for (auto const &cand : avail_walls) {
		bool const dim(cand.second >> 1), dir(cand.second & 1);
		float const wall_len(room_sz[!dim]), wall_edge(room_bounds.d[dim][dir]); // no trim padding
		float const height(0.34*floor_spacing), depth(0.74*height), length(rgen.rand_uniform(0.8, 1.0)*floor_spacing), dsign(dir ? -1.0 : 1.0);
		unsigned flags(is_house ? RO_FLAG_IS_HOUSE : 0);
		cube_t vanity(room_bounds);
		vanity.expand_in_dim(!dim, -get_trim_thickness()); // avoid Z-fighting with exterior wall
		set_cube_zvals(vanity, zval, zval+height);
		vanity.d[dim][!dir] = wall_edge + dsign*depth;

		if (wall_len <= length) { // short wall: vanity is full width
			flags |= (RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP); // flag to not draw either end
		}
		else if (wall_len < 1.5*length) { // middle size: vanity against 2 walls/corner
			bool const shift_dir(rgen.rand_bool());
			vanity.d[!dim][shift_dir] += (shift_dir ? -1.0 : 1.0)*(wall_len - length);
			flags |= (shift_dir ? RO_FLAG_ADJ_BOT : RO_FLAG_ADJ_TOP); // flag to not draw one end
		}
		else { // long wall: vanity centered
			set_wall_width(vanity, (room_bounds.get_center_dim(!dim) + 0.1*rgen.signed_rand_float()*length), 0.5*length, !dim);
		}
		cube_t blocker(vanity);
		blocker.d[dim][!dir] += 0.9*dsign*depth;
		blocker.expand_in_dim(!dim, 0.05*depth); // include overhang on sides
		if (overlaps_obj_or_placement_blocked(blocker, room, objs_start)) continue; // bad placement; need to check elevators for hospitals
		// Note: has doors but no drawers
		objs.emplace_back(vanity,  TYPE_VANITY,  room_id, dim, !dir, flags, tot_light_amt);
		objs.emplace_back(blocker, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker in front
		// add a mirror/medicine cabinet above the sink; shouldn't need to check for overlaps since vanity is placed first
		cube_t mirror(vanity); // start with the sink left and right position
		mirror.z1() = vanity.z2() + 0.14*floor_spacing;
		mirror.z2() = mirror.z1() + 0.30*floor_spacing;
		set_wall_width(mirror, vanity.get_center_dim(!dim), 0.13*floor_spacing, !dim);
		mirror.d[dim][!dir] = wall_edge + dsign*wall_thickness; // thickness
		unsigned const mirror_flags(RO_FLAG_IS_HOUSE | ((num_ext_walls == 1) ? RO_FLAG_INTERIOR : 0)); // flag as interior if windows are opaque glass blocks
		objs.emplace_back(mirror, TYPE_MED_CAB, room_id, dim, !dir, mirror_flags, tot_light_amt);
		set_obj_id(objs); // for crack texture selection/orient
		room.set_has_mirror();
		return 1; // success/done
	} // for cand
	return 0;
}

void building_t::add_bathroom_plumbing(room_object_t const &obj) { // only water pipes; drains are added when the plumbing fixture is removed
	assert(has_room_geom());
	bool const dim(obj.dim), dir(obj.dir);
	float const wall_thickness(get_wall_thickness());
	float pipe_radius(0.12*wall_thickness);
	unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_HANGING | RO_FLAG_LIT); // RO_FLAG_LIT means the pipe casts shadows
	flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO); // cap the exposed end so that the wall and interior aren't visible
	cube_t pipes[2];
	point pipe_p1;
	pipe_p1[!dim] = obj.get_center_dim(!dim); // starts centered
	pipe_p1[ dim] = obj.d[dim][!dir]; // on the back facing the wall

	if (obj.type == TYPE_TOILET) {
		pipe_p1.z = obj.z1() + 0.52*obj.dz(); // directly runs to the tank, not on the bottom like most toilets
		pipes[0].set_from_point(pipe_p1);
	}
	else if (obj.type == TYPE_SINK) {
		float const spacing(0.2*obj.get_width());
		pipe_p1.z = obj.z1() + 0.78*obj.dz(); // near the top
		pipe_p1[!dim] -= 0.5*spacing;
		pipes[0].set_from_point(pipe_p1);
		pipe_p1[!dim] += 1.0*spacing;
		pipes[1].set_from_point(pipe_p1);
	}
	else if (obj.type == TYPE_URINAL) {
		pipe_p1.z = obj.z1() + 0.75*obj.dz(); // near the top
		pipes[0].set_from_point(pipe_p1);
		pipe_radius *= 1.25; // wider
	}
	else if (obj.type == TYPE_TUB || obj.type == TYPE_SHOWERTUB) {
		// not yet handled because a tub/shower+tub can't be taken
	}
	else {assert(0);} // not a plumbing fixture

	for (unsigned d = 0; d < 2; ++d) { // this loop invalidates obj
		cube_t &pipe(pipes[d]);
		if (pipe.is_all_zeros()) continue;
		pipe.expand_in_dim( dim, 0.5*wall_thickness); // set length
		pipe.expand_in_dim(!dim, pipe_radius);
		pipe.expand_in_dim(2,    pipe_radius);
		interior->room_geom->objs.emplace_back(pipe, TYPE_PIPE, obj.room_id, dim, 0, flags, obj.light_amt, SHAPE_CYLIN, COPPER_C); // horizontal
	} // for d
}

bool building_t::add_tp_roll(cube_t const &room, unsigned room_id, float tot_light_amt, bool dim, bool dir, float length, float zval, float wall_pos, bool check_valid_pos) {
	float const diameter(length);
	cube_t tp;
	set_cube_zvals(tp, zval, (zval + diameter));
	set_wall_width(tp, wall_pos, 0.5*length, !dim); // set length
	tp.d[dim][ dir] = room.d[dim][dir]; // against the wall
	tp.d[dim][!dir] = tp  .d[dim][dir] + (dir ? -1.0 : 1.0)*1.1*diameter; // set the diameter with a bit of extra space for clearance
	// Note: not checked against other bathroom objects because the toilet is placed first
	if (check_valid_pos && (!room.contains_cube(tp) || is_obj_placement_blocked(tp, room, 1))) return 0;
	if (has_small_part && !check_if_placed_on_interior_wall(tp, get_room(room_id), dim, dir))  return 0; // need to check for missing walls
	interior->room_geom->objs.emplace_back(tp, TYPE_TPROLL, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, WHITE);
	set_obj_id(interior->room_geom->objs);
	return 1;
}

bool building_t::divide_bathroom_into_stalls(rand_gen_t &rgen, room_t &room, float zval, unsigned room_id,
	float tot_light_amt, unsigned floor, unsigned lights_start, unsigned lights_end)
{
	// Note: assumes no prior placed objects
	bool const use_sink_model(0);
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	vector3d const tsz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOILET)); // L, W, H
	float const theight(0.35*floor_spacing), twidth(theight*tsz.y/tsz.z), tlength(theight*tsz.x/tsz.z), stall_depth(2.2*tlength);
	float sheight(0), swidth(0), slength(0), uheight(0), uwidth(0), ulength(0);

	if (use_sink_model) {
		vector3d const ssz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SINK)); // L, W, H
		sheight = 0.45*floor_spacing; swidth = sheight*ssz.y/ssz.z; slength = sheight*ssz.x/ssz.z;
	}
	else {
		sheight = 0.36*floor_spacing; swidth = 0.3*floor_spacing; slength = 0.32*floor_spacing;
		//slength = (has_parking_garage ? (tlength + 2.0*wall_thickness) : 0.32*floor_spacing); // align sink drain to toilets for parking garage pipes?
	}
	float stall_width(2.0*twidth), sink_spacing(1.75*swidth);
	bool br_dim(room.dy() < room.dx()), sink_side(0), sink_side_set(0), mens_room(0); // br_dim is the smaller dim
	cube_t place_area(room), avoid;
	place_area.expand_by(-0.5*wall_thickness);
	unsigned br_door_stack_ix(0);

	// check for any stairs or elevators that may partiall overlap a bathroom wall from the adjacent room and avoid them; there should be at most one
	for (stairwell_t const &s : interior->stairwells) {
		if (s.intersects(room) && zval >= s.z1() && zval < s.z2()) {avoid = s; break;}
	}
	for (elevator_t const &e : interior->elevators) {
		if (e.intersects(room) && zval >= e.z1() && zval < e.z2()) {avoid = e; break;}
	}
	// determine men's room vs. women's room
	room_type const init_rtype(room.get_room_type(floor));

	if (room.is_rtype_locked(floor) && (init_rtype == RTYPE_MENS || init_rtype == RTYPE_WOMENS)) { // pre-assigned
		mens_room = (init_rtype == RTYPE_MENS);
	}
	else { // determine if this is the men's or women's room
		point const part_center(get_part_for_room(room).get_cube_center()), room_center(room.get_cube_center());
		mens_room = ((part_center.x < room_center.x) ^ (part_center.y < room_center.y));

		if (!is_cube()) { // alternate between men's and women's restrooms
			if      (interior->room_geom->mens_count < interior->room_geom->womens_count) {mens_room = 1;}
			else if (interior->room_geom->mens_count > interior->room_geom->womens_count) {mens_room = 0;}
		}
		else {
			bool has_second_bathroom(0);

			// if there are two bathrooms (one on each side of the building), assign a gender to each side; if only one, alternate gender per floor
			for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
				if (r->part_id != room.part_id || &(*r) == &room) continue; // different part or same room
				if (is_room_office_bathroom(*r, zval, floor)) {has_second_bathroom = 1; break;}
			}
			if (!has_second_bathroom) {mens_room ^= (floor & 1);}
		}
	}
	bool const add_urinals(mens_room && building_obj_model_loader.is_model_valid(OBJ_MODEL_URINAL));

	if (add_urinals) { // use urinal model
		vector3d const usz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_URINAL)); // L, W, H
		uheight = 0.4*floor_spacing; uwidth = uheight*usz.y/usz.z; ulength = uheight*usz.x/usz.z;
	}
	for (unsigned d = 0; d < 2 && !sink_side_set; ++d) { // try long then short dim
		bool first_sink_side(0);

		if (d == 0 && point_in_mall(room.get_cube_center())) { // place sink on the side closer to the mall concourse
			first_sink_side = (room.get_center_dim(!br_dim) < get_mall_concourse().get_center_dim(!br_dim));
		}
		for (unsigned e = 0; e < 2 && !sink_side_set; ++e) {
			bool const side(bool(e) ^ first_sink_side);
			cube_t c(room);
			set_cube_zvals(c, zval, zval+wall_thickness); // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
			c.d[!br_dim][!side] = c.d[!br_dim][side] + (side ? -1.0 : 1.0)*wall_thickness; // shrink to near zero area in this dim

			for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
				if (i->dim == br_dim) continue; // door in wrong dim
				if (!i->is_connected_to_room(room_id) || !is_cube_close_to_door(c, 0.0, 0, *i, 2)) continue; // check both dirs
				sink_side = side; sink_side_set = 1;
				place_area.d[!br_dim][side] += (sink_side ? -1.0 : 1.0)*(i->get_sz_dim(br_dim) - 0.25*swidth); // add sink clearance for the door to close
				br_door_stack_ix = (i - interior->door_stacks.begin());
				break; // sinks are on the side closest to the door
			}
		} // for e
		if (d == 0 && !sink_side_set) {br_dim ^= 1;} // door not found on long dim - R90 and try short dim
	} // for d
	assert(sink_side_set);
	float const room_len(place_area.get_sz_dim(!br_dim)), room_width(place_area.get_sz_dim(br_dim));
	float const sinks_len(0.4*room_len), stalls_len(room_len - sinks_len), req_depth(2.0f*max(stall_depth, slength));
	if (room_width < req_depth) return 0;
	if (sinks_len < 2.0*sink_spacing) {sink_spacing *= 0.8;} // reduce sink spacing a bit to try and fit at least two
	unsigned const num_stalls(std::floor(stalls_len/stall_width)), num_sinks(std::floor(sinks_len/sink_spacing));
	if (num_stalls < 2 || num_sinks < 1) return 0; // not enough space for 2 stalls and a sink
	stall_width  = stalls_len/num_stalls; // reclaculate to fill the gaps
	sink_spacing = sinks_len/num_sinks;
	bool const two_rows(room_width > 1.5*req_depth), skip_stalls_side(room_id & 1); // put stalls on a side consistent across floors
	float const sink_side_sign(sink_side ? 1.0 : -1.0), stall_step(sink_side_sign*stall_width), sink_step(-sink_side_sign*sink_spacing);
	float const floor_thickness(get_floor_thickness());
	unsigned const NUM_STALL_COLORS = 4;
	colorRGBA const stall_colors[NUM_STALL_COLORS] = {colorRGBA(0.75, 1.0, 0.9, 1.0), colorRGBA(0.7, 0.8, 1.0), WHITE, DK_GRAY}; // blue-green, light blue
	colorRGBA const stall_color(stall_colors[interior->doors.size() % NUM_STALL_COLORS]); // random, but constant for each building
	vect_room_object_t &objs(interior->room_geom->objs);
	room_object_t mirrors[2]; // candidate mirrors for each dir
	++(mens_room ? interior->room_geom->mens_count : interior->room_geom->womens_count);

	for (unsigned dir = 0; dir < 2; ++dir) { // each side of the wall
		if (!two_rows && dir == (unsigned)skip_stalls_side) continue; // no stalls/sinks on this side
		// add stalls
		float const dir_sign(dir ? -1.0 : 1.0), wall_pos(place_area.d[br_dim][dir]), stall_from_wall(wall_pos + dir_sign*(0.5*tlength + wall_thickness));
		float stall_pos(place_area.d[!br_dim][!sink_side] + 0.5*stall_step);

		for (unsigned n = 0; n < num_stalls; ++n, stall_pos += stall_step) {
			point center(stall_from_wall, stall_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t toilet(center), stall(toilet);
			toilet.expand_in_dim( br_dim, 0.5*tlength);
			toilet.expand_in_dim(!br_dim, 0.5*twidth);
			toilet.z2() += theight;
			stall.z2() = stall.z1() + floor_spacing - floor_thickness; // set stall height to room height
			stall.expand_in_dim(!br_dim, 0.5*stall_width);
			stall.d[br_dim][ dir] = wall_pos; // + wall_thickness?
			stall.d[br_dim][!dir] = wall_pos + dir_sign*stall_depth; // extend the front outward from the wall
			if (interior->is_cube_close_to_doorway(stall, room, 0.0, 1)) continue; // skip if close to a door (for rooms with doors at both ends); inc_open=1
			if (!avoid.is_all_zeros() && stall.intersects(avoid))        continue;

			if (!is_cube()) {
				cube_t stall_exp(stall);
				stall_exp.d[br_dim][!dir] += dir_sign*1.0*get_scaled_player_radius(); // make sure the player can fit
				if (!check_cube_within_part_sides(stall_exp)) continue; // outside the building
			}
			bool const is_open(rgen.rand_bool()); // 50% chance of stall door being open
			bool const out_of_order(!is_open && !room.is_ext_basement() && rgen.rand_float() < 0.2); // not for mall bathrooms
			unsigned const flags(out_of_order ? RO_FLAG_BROKEN : 0); // toilet can't be flushed and door can't be opened if out of order
			objs.emplace_back(toilet, TYPE_TOILET, room_id, br_dim, !dir, flags, tot_light_amt);
			add_bathroom_plumbing(objs.back());
			objs.emplace_back(stall,  TYPE_STALL,  room_id, br_dim,  dir, (flags | (is_open ? RO_FLAG_OPEN : 0)), tot_light_amt, SHAPE_CUBE, stall_color);

			if (out_of_order) { // add out-of-order sign
				cube_t stall_clipped(stall);
				stall_clipped.z2() -= 0.33*stall.dz(); // make it shorter; really only need the stall door, but the dim dimension size isn't used anyway
				add_out_or_order_sign(stall_clipped, br_dim, !dir, room_id);
			}
			float const tp_length(0.18*theight), wall_pos(toilet.get_center_dim(br_dim));
			cube_t stall_inner(stall);
			stall_inner.expand_in_dim(!br_dim, -0.0125*stall.dz()); // subtract off stall wall thickness
			add_tp_roll(stall_inner, room_id, tot_light_amt, !br_dim, dir, tp_length, (zval + 0.7*theight), wall_pos);

			// move any lights that intersect stalls, since stalls extend to the ceiling
			for (auto i = objs.begin()+lights_start; i < objs.begin()+lights_end; ++i) {
				if (i->type != TYPE_LIGHT || !i->intersects(stall)) continue;
				float const move_amt(stall.d[br_dim][!dir] - i->d[br_dim][dir] + dir_sign*wall_thickness);
				i->translate_dim(br_dim, move_amt); // shouldn't intersect anything else
			}
		} // for n
		if (add_urinals && dir == (unsigned)skip_stalls_side) continue; // no urinals and sinks are each on one side
		// add sinks
		float const sink_start(place_area.d[!br_dim][sink_side] + 0.5f*sink_step);
		float const sink_from_wall(wall_pos + dir_sign*(0.5f*slength + (use_sink_model ? wall_thickness : 0.0f)));
		float sink_pos(sink_start);
		bool hit_mirror_end(0);
		unsigned last_sink_ix(0);
		cube_t sinks_bcube;

		for (unsigned n = 0; n < num_sinks; ++n, sink_pos += sink_step) {
			point center(sink_from_wall, sink_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t sink(center);
			sink.expand_in_dim(br_dim, 0.5*slength);
			sink.z2() += sheight;
			if (interior->is_cube_close_to_doorway(sink, room, 0.0, 1)) continue; // skip if close to a door, inc_open=1, pre expand
			sink.expand_in_dim(!br_dim, 0.5*(use_sink_model ? swidth : fabs(sink_step))); // tile exactly with the adjacent sink
			if (interior->is_cube_close_to_doorway(sink, room, 0.0, 0)) continue; // skip if close to a door
			if (!check_cube_within_part_sides(sink))                    continue; // outside the building
			if (!avoid.is_all_zeros() && sink.intersects(avoid))        continue;
			if (use_sink_model) {objs.emplace_back(sink, TYPE_SINK,   room_id, br_dim, !dir, 0, tot_light_amt);} // sink 3D model
			else                {objs.emplace_back(sink, TYPE_BRSINK, room_id, br_dim, !dir, 0, tot_light_amt);} // flat basin sink
			if (use_sink_model) {add_bathroom_plumbing(objs.back());}
			// if we started the mirror, but we have a gap with no sink (blocked by a door, etc.), then end the mirror
			hit_mirror_end |= (n > last_sink_ix+1 && !sinks_bcube.is_all_zeros());
			if (!hit_mirror_end) {sinks_bcube.assign_or_union_with_cube(sink);}
			last_sink_ix = n;
		} // for n
		if (add_urinals) { // add urinals opposite the sinks, using same spacing as sinks
			float const u_wall(place_area.d[br_dim][!dir]), u_from_wall(u_wall - dir_sign*(0.5*ulength + 0.01*wall_thickness));
			float u_pos(sink_start);
			cube_t sep_wall;
			set_cube_zvals(sep_wall, zval+0.15*uheight, zval+1.25*uheight);
			sep_wall.d[br_dim][!dir] = u_wall;
			sep_wall.d[br_dim][ dir] = u_wall - dir_sign*0.25*floor_spacing;
			// school bathrooms have a short urinal at one end if there is more than one urinal
			unsigned const short_urinal_ix((is_school() && num_sinks > 1) ? ((room_id & 1) ? 0 : num_sinks-1) : num_sinks);

			for (unsigned n = 0; n < num_sinks; ++n, u_pos += sink_step) {
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				point center(u_from_wall, u_pos, (zval + ((n == short_urinal_ix) ? 0.125 : 0.25)*uheight));
				if (br_dim) {swap(center.x, center.y);} // R90 about z
				cube_t urinal(center);
				urinal.expand_in_dim( br_dim, 0.5*ulength);
				urinal.expand_in_dim(!br_dim, 0.5*uwidth);
				urinal.z2() += uheight;
				if (interior->is_cube_close_to_doorway(urinal, room, 0.0, 1)) continue; // skip if close to a door
				if (!check_cube_within_part_sides(urinal))                    continue; // outside the building
				if (!avoid.is_all_zeros() && urinal.intersects(avoid))        continue;
				
				if (!interior->is_cube_close_to_doorway(sep_wall, room, 0.0, 1)) { // check for doors, when the bathroom door is not centered on the room
					objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, 0, tot_light_amt, SHAPE_SHORT, stall_color);
				}
				objs.emplace_back(urinal, TYPE_URINAL, room_id, br_dim,  dir, 0, tot_light_amt);
				add_bathroom_plumbing(objs.back());
			} // for n
			if (!two_rows) { // skip first wall if adjacent to a stall
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				if (!avoid.is_all_zeros() && sep_wall.intersects(avoid)) {} // skip
				else {objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, 0, tot_light_amt, SHAPE_SHORT, stall_color);}
			}
		}
		if (!sinks_bcube.is_all_zeros()) { // add a long mirror above the sink
			cube_t mirror(sinks_bcube);
			mirror.expand_in_dim(!br_dim, -0.25*wall_thickness); // slightly smaller
			mirror.d[br_dim][ dir] = wall_pos;
			mirror.d[br_dim][!dir] = wall_pos + dir_sign*0.1*wall_thickness;
			mirror.z1() = sinks_bcube.z2() + 0.25*floor_thickness;
			mirror.z2() = zval + 0.9*floor_spacing - floor_thickness;
			if (mirror.is_strictly_normalized()) {mirrors[dir] = room_object_t(mirror, TYPE_MIRROR, room_id, br_dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);}
		}
	} // for dir
	for (unsigned d = 0; d < 2; ++d) { // each candidate mirror
		if (mirrors[d].is_all_zeros()) continue;
		if (ENABLE_MIRROR_REFLECTIONS && d == (unsigned)skip_stalls_side && !mirrors[!d].is_all_zeros()) continue; // select a single side if reflections are enabled
		objs.push_back(mirrors[d]);
		set_obj_id(objs); // for crack texture selection/orient
		room.set_has_mirror();
	}
	room_type const rtype(mens_room ? RTYPE_MENS : RTYPE_WOMENS);
	room.assign_to(rtype, floor);
	
	// make sure doors start closed and unlocked, and flag them as auto_close;
	// if (!is_cube()) we also want to make sure the door opens inward, but we can't change it for only one door in the stack
	for (door_stack_t const &ds : interior->door_stacks) {
		if (!ds.is_connected_to_room(room_id)) continue;
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() > zval || door.z2() < zval) continue; // wrong floor
			door.rtype  = rtype; // tag type so that the correct sign is added
			door.locked = 0; // only needed for non-cube buildings
			door.make_auto_close();
			if (!is_cube() && !door_opens_inward(ds, room)) {door.opens_out_of_br = 1;}
		} // for ds
	} // for ds
	// add a sign outside the bathroom door; should schools use boys/girls? but it doesn't match the text on the signs placed on the doors
	string const sign_text(/*is_school() ? (mens_room ? "Boys" : "Girls") :*/ (mens_room ? "Men" : "Women"));
	add_door_sign(sign_text, room, zval, room_id, 1); // no_check_adj_walls=1

	// make this door/room out of order 10% of the time; only for cube buildings (others need the connectivity), and not for mall bathrooms
	if (is_cube() && !(has_mall() && room.is_ext_basement()) && rgen.rand_float() < 0.1) {
		make_door_out_or_order(room, zval, room_id, br_door_stack_ix);
		room.set_has_out_of_order(); // flag if any floor is out of order
	}
	return 1;
}

void add_door_if_blocker(cube_t const &door, cube_t const &room, bool inc_open, bool dir, bool hinge_side, bool exterior, vect_cube_t &blockers) {
	bool const dim(door.dy() < door.dx()), edir(dim ^ dir ^ hinge_side ^ 1);
	float const width(door.get_sz_dim(!dim));
	cube_t door_exp(door);
	door_exp.expand_in_dim(dim, width); // width
	if (!door_exp.intersects(room)) return; // check against room before expanding along wall to exclude doors in adjacent rooms
	if (exterior) {door_exp.expand_in_dim(dim, width);} // add extra clearance in front
	door_exp.expand_in_dim(!dim, width*0.25); // min expand value
	if (inc_open) {door_exp.d[!dim][edir] += (edir ? 1.0 : -1.0)*0.75*width;} // expand the remainder of the door width in this dir
	blockers.push_back(door_exp);
}
// used for kitchen cabinets
int building_t::gather_room_placement_blockers(cube_t const &room, unsigned objs_start, vect_cube_t &blockers, bool inc_open_doors, bool ignore_chairs) const {
	assert(has_room_geom());
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(objs_start <= objs.size());
	blockers.clear();
	int table_blocker_ix(-1);

	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (ignore_chairs && i->type == TYPE_CHAIR) continue;
		
		if (!i->no_coll() && i->intersects(room)) {
			if (i->type == TYPE_TABLE) {table_blocker_ix = int(blockers.size());} // track which blocker is the table, for use with kitchen counters
			blockers.push_back(*i);
		}
	}
	for (auto const &door : doors) { // exterior doors
		add_door_if_blocker(door.get_bcube(), room, 0, 0, 0, 1, blockers);
	}
	for (door_stack_t const &ds : interior->door_stacks) { // interior doors
		add_door_if_blocker(ds, room, door_opens_inward(ds, room), ds.open_dir, ds.hinge_side, 0, blockers);
	}
	// Note: caller must call is_blocked_by_stairs_or_elevator() to handle stairs and elevators
	return table_blocker_ix;
}

bool building_t::add_kitchen_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool allow_adj_ext_door) {
	// Note: table and chairs have already been placed
	bool const residential(is_residential());
	if (room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be kitchens
	if (!residential && rgen.rand_bool()) return 0; // some office buildings have kitchens, allow it half the time
	// if it has an external door then reject the room half the time; most houses don't have a front door to the kitchen
	if (is_room_adjacent_to_ext_door(room, zval, 1) && (!allow_adj_ext_door || rgen.rand_bool())) return 0; // front_door_only=1
	float const wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall for appliances
	vect_room_object_t &objs(interior->room_geom->objs);
	bool placed_obj(0);
	placed_obj |= place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, room, 0.75, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.2, 4, 0, WHITE, 1); // not at window
	
	if (residential) { // try to place a stove
		unsigned const stove_ix(objs.size()); // can't use objs.back() because there's a blocker
		
		if (place_model_along_wall(OBJ_MODEL_STOVE, TYPE_STOVE, room, 0.46, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0)) {
			assert(stove_ix < objs.size());

			if (building_obj_model_loader.is_model_valid(OBJ_MODEL_HOOD)) { // add hood above the stove
				room_object_t const &stove(objs[stove_ix]);
				vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_HOOD)); // D, W, H
				float const width(stove.get_sz_dim(!stove.dim)), height(width*sz.z/sz.y), depth(width*sz.x/sz.y); // scale to the width of the stove
				float const ceiling_z(zval + get_floor_ceil_gap()), z_top(ceiling_z + get_fc_thickness()); // shift up a bit because it's too low
				cube_t hood(stove);
				set_cube_zvals(hood, z_top-height, z_top);
				hood.d[stove.dim][stove.dir] = stove.d[stove.dim][!stove.dir] + (stove.dir ? 1.0 : -1.0)*depth;

				if (!check_if_against_window(hood, room, stove.dim, !stove.dir)) { // don't place hoods against windows
					objs.emplace_back(hood, TYPE_HOOD, room_id, stove.dim, stove.dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, LT_GRAY);
				
					if (has_attic()) { // add rooftop vent above the hood; should be close to the edge of the house in many cases
						float const attic_floor_zval(get_attic_part().z2()), vent_radius(0.075*(width + depth)); // kitchen can be in a part below the one with the attic
						point const vent_bot_center(hood.xc(), hood.yc(), attic_floor_zval);
						add_attic_roof_vent(vent_bot_center, vent_radius, room_id, 1.0); // light_amt=1.0; room_id is for the kitchen because there's no attic room
					}
				}
			}
			if (!rgen.rand_bool()) { // maybe add a pan on one of the stove burners
				room_object_t const &stove(objs[stove_ix]);
				float const stove_height(stove.dz()), delta_z(0.018*stove_height);
				float const pan_radius(rgen.rand_uniform(0.075, 0.09)*stove_height), pan_height(rgen.rand_uniform(0.035, 0.045)*stove_height);
				point locs[4];
				get_stove_burner_locs(stove, locs);
				unsigned const burner_ix(rgen.rand() & 3); // 0-3
				point &loc(locs[burner_ix]);
				loc.z += delta_z;
				cube_t const pan(get_cube_height_radius(loc, pan_radius, pan_height));
				objs.emplace_back(pan, TYPE_PAN, room_id, stove.dim, stove.dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, DK_GRAY);
			}
			placed_obj = 1;
		}
	}	
	if (residential && placed_obj) { // if we have at least a fridge or stove, try to add countertops
		float const vspace(get_window_vspace()), height(0.345*vspace), depth(0.74*height), floor_thickness(get_floor_thickness());
		float const min_clearance(get_min_front_clearance_inc_people());
		float min_hwidth(0.6*height), front_clearance(max(0.6f*height, min_clearance));
		cube_t cabinet_area(room_bounds);
		cabinet_area.expand_by(-0.05*wall_thickness); // smaller gap than place_area; this is needed to prevent z-fighting with exterior walls
		if (min(cabinet_area.dx(), cabinet_area.dy()) < 4.0*min_hwidth) return placed_obj; // no space for cabinets, room is too small
		unsigned const counters_start(objs.size());
		cube_t c;
		set_cube_zvals(c, zval, zval+height);
		set_cube_zvals(cabinet_area, zval, (zval + vspace - floor_thickness));
		static vect_cube_t blockers;
		int const table_blocker_ix(gather_room_placement_blockers(cabinet_area, objs_start, blockers, 1, 1)); // inc_open_doors=1, ignore_chairs=1
		bool const have_toaster(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOASTER));
		vector3d const toaster_sz(have_toaster ? building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOASTER) : zero_vector); // L, D, H
		bool is_sink(1), placed_mwave(0), placed_toaster(0), had_counter(0);
		unsigned num_paper_towels(0);
		cube_t mwave, toaster;

		for (unsigned n = 0; n < (had_counter ? 50U : 60U); ++n) { // 50 attempts, plus an extra 10 if no counters were placed
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
			if (room.has_open_wall(dim, dir)) continue; // don't place against open walls
			bool const is_ext_wall(classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT); // assumes not in basement
			// only consider exterior walls in the first 20 attempts to prioritize these so that we don't have splits visible through windows; also places kitchen sinks near windows
			if (n < 20 && !is_ext_wall) continue;

			if (n > 30 && !had_counter) { // reduce min width and clearance for later iterations if no counters were placed
				min_hwidth      -= 0.01*height;
				front_clearance -= 0.01*height;
			}
			float const center(rgen.rand_uniform(cabinet_area.d[!dim][0]+min_hwidth, cabinet_area.d[!dim][1]-min_hwidth)); // random position
			float const dir_sign(dir ? -1.0 : 1.0), wall_pos(cabinet_area.d[dim][dir]), front_pos(wall_pos + dir_sign*depth);
			c.d[ dim][ dir] = wall_pos;
			c.d[ dim][!dir] = front_pos + dir_sign*front_clearance;
			c.d[!dim][   0] = center - min_hwidth;
			c.d[!dim][   1] = center + min_hwidth;
			cube_t c_min(c); // min runlength - used for collision tests
			for (unsigned e = 0; e < 2; ++e) {c.d[!dim][e] = cabinet_area.d[!dim][e];} // start at full room width
			bool bad_place(0);

			for (auto i = blockers.begin(); i != blockers.end(); ++i) {
				cube_t b(*i); // expand tables by an extra clearance to allow the player to fit in the diagonal gap between the table and the counter
				if (int(i - blockers.begin()) == table_blocker_ix) {b.expand_in_dim(!dim, min_clearance);}
				if (!b.intersects(c)) continue; // optimization - no cube interaction
				if (b.intersects(c_min)) {bad_place = 1; break;}
				if (b.d[!dim][1] < c_min.d[!dim][0]) {max_eq(c.d[!dim][0], b.d[!dim][1]);} // clip on lo side
				if (b.d[!dim][0] > c_min.d[!dim][1]) {min_eq(c.d[!dim][1], b.d[!dim][0]);} // clip on hi side
			} // for i
			if (bad_place) continue;
			if (interior->is_blocked_by_stairs_or_elevator(c)) continue; // check stairs
			assert(c.contains_cube(c_min));

			for (auto i = objs.begin()+objs_start; i != objs.begin()+counters_start; ++i) { // remove any chairs blocking the counter doors/drawers
				if (i->type == TYPE_CHAIR && i->intersects(c)) {i->remove();}
			}
			c.d[dim][!dir] = front_pos; // remove front clearance
			bool const add_backsplash(!is_ext_wall); // only add to interior walls to avoid windows; assuming not in basement

			for (auto i = objs.begin()+counters_start; i != objs.end(); ++i) { // find adjacencies to previously placed counters and flag to avoid placing doors
				if (i->dim == dim) continue; // not perpendicular
				if (i->d[!i->dim][dir] != wall_pos) continue; // not against the wall on this side
				if (i->d[i->dim][i->dir] != c.d[!dim][0] && i->d[i->dim][i->dir] != c.d[!dim][1]) continue; // not adjacent
				i->flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
				if (add_backsplash) {i->flags |= RO_FLAG_HAS_EXTRA;}
			}
			unsigned const cabinet_id(objs.size());
			objs.emplace_back(c, (is_sink ? TYPE_KSINK : TYPE_COUNTER), room_id, dim, !dir, 0, tot_light_amt);
			set_obj_id(objs);
			had_counter = 1;
			
			if (add_backsplash) {
				objs.back().flags |= (RO_FLAG_ADJ_BOT | RO_FLAG_HAS_EXTRA); // flag back as having a back backsplash
				cube_t bs(c);
				bs.z1()  = c.z2();
				bs.z2() += BACKSPLASH_HEIGHT*c.dz();
				bs.d[dim][!dir] -= (dir ? -1.0 : 1.0)*0.99*depth; // matches building_room_geom_t::add_counter()
				objs.emplace_back(bs, TYPE_BLOCKER, room_id, dim, !dir, RO_FLAG_INVIS); // add blocker to avoid placing light switches here
			}
			// add upper cabinets
			cube_t c2(c);
			set_cube_zvals(c2, (zval + 0.65*vspace), cabinet_area.z2()); // up to the ceiling

			if (is_ext_wall) { // possibly against a window
				max_eq(c2.z1(), (c2.z2() - vspace*get_window_v_border() + 0.5f*floor_thickness)); // increase bottom of cabinet to the top of the window
			}
			if (c2.dz() > 0.1*vspace && !has_bcube_int_no_adj(c2, blockers)) { // add if it's not too short and not blocked
				objs.emplace_back(c2, TYPE_CABINET, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // no collision detection
				set_obj_id(objs);
			}
			blockers.push_back(c); // add to blockers so that later counters don't intersect this one
			bool added_obj(0);

			// place a microwave on a counter 50% of the time
			if (!is_sink && !placed_mwave && c.get_sz_dim(!dim) > 0.5*vspace && rgen.rand_bool()) {
				float const mheight(rgen.rand_uniform(1.0, 1.2)*0.14*vspace), mwidth(1.7*mheight), mdepth(1.2*mheight); // fixed AR=1.7 to match the texture
				float const pos(rgen.rand_uniform((c.d[!dim][0] + 0.6*mwidth), (c.d[!dim][1] - 0.6*mwidth)));
				set_cube_zvals(mwave, c.z2(), c.z2()+mheight);
				set_wall_width(mwave, pos, 0.5*mwidth, !dim);
				mwave.d[dim][ dir] = wall_pos + dir_sign*0.05*mdepth;
				mwave.d[dim][!dir] = mwave.d[dim][dir] + dir_sign*mdepth;
				objs.emplace_back(mwave, TYPE_MWAVE, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);
				objs[cabinet_id].flags |= RO_FLAG_ADJ_TOP; // flag as having a microwave so that we don't add a book or bottle that could overlap it
				placed_mwave = added_obj = 1;
			}
			// place a toaster on a counter 90% of the time
			if (!is_sink && !placed_toaster && have_toaster && rgen.rand_float() < 0.9) {
				float const theight(0.09*vspace), twidth(theight*toaster_sz.x/toaster_sz.z), tdepth(theight*toaster_sz.y/toaster_sz.z);

				if (c.get_sz_dim(!dim) > 1.25*twidth && c.get_sz_dim(dim) > 1.25*tdepth) { // add if it fits
					float const pos_w(rgen.rand_uniform((c.d[!dim][0] + 0.6*twidth), (c.d[!dim][1] - 0.6*twidth)));
					float const pos_d(rgen.rand_uniform((c.d[ dim][0] + 0.6*tdepth), (c.d[ dim][1] - 0.6*tdepth)));
					set_cube_zvals(toaster, c.z2(), c.z2()+theight);
					set_wall_width(toaster, pos_w, 0.5*twidth, !dim);
					set_wall_width(toaster, pos_d, 0.5*tdepth,  dim);

					if (!placed_mwave || !mwave.intersects(toaster)) { // don't overlap the microwave
						objs.emplace_back(toaster, TYPE_TOASTER, room_id, !dim, rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt); // random dir
						objs.back().color = toaster_colors[rgen.rand()%NUM_TOASTER_COLORS];
						objs[cabinet_id].flags |= RO_FLAG_ADJ_TOP; // flag as having an object so that we don't add a book or bottle that could overlap it
						placed_toaster = added_obj = 1;
					}
				}
			}
			// place a paper towel roll on the counter 75% of the time if there's no toaster, microwave, or sink; max 2 per kitchen
			if (!is_sink && !added_obj && num_paper_towels < 2 && rgen.rand_float() < 0.75) {
				float const ptheight(0.115*vspace), ptradius(0.25*ptheight);

				if (min(c.dx(), c.dy()) > 2.5*ptradius) { // add if it fits
					cube_t ptroll;
					gen_xy_pos_for_round_obj(ptroll, c, ptradius, ptheight, 1.2*ptradius, rgen);
					objs.emplace_back(ptroll, TYPE_TPROLL, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_HAS_EXTRA), tot_light_amt);
					objs[cabinet_id].flags |= RO_FLAG_ADJ_TOP; // flag as having an object so that we don't add a book or bottle that could overlap it
					++num_paper_towels;
				}
			}
			if (is_sink) { // kitchen sink; add cups, plates, and cockroaches
				cube_t sink(get_sink_cube(objs[cabinet_id]));
				sink.z2() = sink.z1(); // shrink to zero area at the bottom
				unsigned const objs_start(objs.size()), num_objs(1 + rgen.rand_bool()); // 1-2 objects

				for (unsigned n = 0; n < num_objs; ++n) {
					unsigned const obj_type(rgen.rand()%3);
					static vect_cube_t avoid;
					avoid.clear();
					if (objs.size() > objs_start) {avoid.push_back(objs.back());} // avoid the last object that was placed, if there was one

					if      (obj_type == 0) {place_plate_on_obj(rgen, sink, room_id, tot_light_amt, avoid);} // add a plate
					else if (obj_type == 1) {place_cup_on_obj  (rgen, sink, room_id, tot_light_amt, avoid, 1);} // add a cup; make_empty=1
					else if (obj_type == 2 && building_obj_model_loader.is_model_valid(OBJ_MODEL_ROACH)) { // add a cockroach (upside down?)
						sink.d[dim][!dir] = sink.get_center_dim(dim); // use the half area near the back wall to make sure the roach is visible to the player
						cube_t roach;
						float const radius(sink.get_sz_dim(dim)*rgen.rand_uniform(0.08, 0.12)), height(get_cockroach_height_from_radius(radius));
						gen_xy_pos_for_round_obj(roach, sink, radius, height, 1.1*radius, rgen);
						objs.emplace_back(roach, TYPE_ROACH, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt);
						if (rgen.rand_bool()) {objs.back().flags |= RO_FLAG_BROKEN;} // 50% chance it's dead
					}
				} // for n
			}
			is_sink = 0; // sink is in first placed counter only
		} // for n
		if (had_counter && rgen.rand_float() < 0.75) { // maybe place an island near the room center
			room_object_t const &sink(objs[counters_start]); // first counter is the kitchen sink
			bool const dim(sink.dim);

			if (cabinet_area.get_sz_dim(dim) > 6.0*depth) { // kitchen is wide enough
				float const wall_spacing(2.6*depth);
				cube_t island(sink);
				island.translate_dim(dim, (sink.dir ? 1.0 : -1.0)*3.0*depth);
				max_eq(island.d[!dim][0], cabinet_area.d[!dim][0]+wall_spacing);
				min_eq(island.d[!dim][1], cabinet_area.d[!dim][1]-wall_spacing);
				
				if (island.get_sz_dim(!dim) > depth && cabinet_area.contains_cube(island) &&
					!has_bcube_int(island, blockers) && !interior->is_blocked_by_stairs_or_elevator(island))
				{
					for (auto i = objs.begin()+objs_start; i != objs.begin()+counters_start; ++i) { // remove any chairs intersecting the island
						if (i->type == TYPE_CHAIR && i->intersects(island)) {i->remove();}
					}
					objs.emplace_back(island, TYPE_COUNTER, room_id, dim, !sink.dir, 0, tot_light_amt); // faces opposite the sink
					set_obj_id(objs);
				}
			}
		}
	}
	if (placed_obj && building_obj_model_loader.is_model_valid(OBJ_MODEL_BAN_PEEL) && rgen.rand_bool()) { // maybe place a banana peel on the floor
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_BAN_PEEL));
		float const floor_spacing(get_window_vspace()), length(0.083*floor_spacing), width(length*sz.y/sz.x), height(length*sz.z/sz.x);
		cube_t valid_area(place_area);
		valid_area.expand_by_xy(-0.25*floor_spacing); // not too close to a wall
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random orientation
		vector3d size(0.5*length, 0.5*width, height);
		if (dim) {swap(size.x, size.y);}

		if (valid_area.dx() > 2.0*size.x && valid_area.dy() > 2.0*size.y) { // should always be true
			for (unsigned n = 0; n < 4; ++n) { // make 4 attempts to place the object
				cube_t c(gen_xy_pos_in_area(valid_area, size, rgen, zval));
				c.expand_by_xy(size);
				c.z2() += height;
				if (overlaps_other_room_obj(c, objs_start) || is_obj_placement_blocked(c, room, 1)) continue; // bad placement
				objs.emplace_back(c, TYPE_BAN_PEEL, room_id, dim, dir, RO_FLAG_RAND_ROT, tot_light_amt);
				break; // done
			} // for n
		}
	}
	return placed_obj;
}

bool building_t::add_fishtank_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, cube_t const &place_area) {
	vect_room_object_t &objs(interior->room_geom->objs);
	// first add a tall table
	float const floor_spacing(get_window_vspace());
	float const table_height(rgen.rand_uniform(0.22, 0.24)*floor_spacing);
	vector3d const fc_sz_scale(rgen.rand_uniform(0.7, 0.8), rgen.rand_uniform(1.6, 1.8), 1.0); // depth, width, height
	unsigned const table_obj_ix(objs.size());
	// not_at_a_window=1
	if (!place_obj_along_wall(TYPE_TABLE, room, table_height, fc_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start,
		0.5, 1, 4, 0, WHITE, 0, SHAPE_CUBE, 0.0, 0, 1)) return 0; // not_ext_wall=1
	// then add the fishtank
	float const tank_height(rgen.rand_uniform(0.23, 0.25)*floor_spacing);
	assert(table_obj_ix < objs.size());
	room_object_t &table(objs[table_obj_ix]); // table was placed first, then a blocker
	table.shape  = SHAPE_SHORT;
	table.flags |= RO_FLAG_ADJ_TOP; // flag table as having something on it
	cube_t tank(table);
	tank.expand_by_xy(-0.1*min(table.dx(), table.dy()));
	set_cube_zvals(tank, table.z2(), (table.z2() + tank_height));
	cube_t test_cube(tank);
	test_cube.z1() += 0.1*tank.dz(); // shift up so that it doesn't intersect the table
	// check if fishtank overlaps and object, but not the table; can this happen? maybe if in front of a picture or TV?
	if (overlaps_other_room_obj(test_cube, objs_start)) return 0; // should we remove the table or leave it there?
	interior->room_geom->add_pet_container(tank, room_id, tot_light_amt, table.dim, table.dir, 0, objs, rgen); // in_pet_store=0
	return 1;
}
void building_room_geom_t::add_pet_container(cube_t const &tank, unsigned room_id, float tot_light_amt, bool dim, bool dir, bool in_pet_store,
	vect_room_object_t &objects, rand_gen_t &rgen, unsigned animal_type, unsigned shelf_ix)
{
	bool is_cage(0);
	unsigned flags(RO_FLAG_NOCOLL);
	colorRGBA color(WHITE);

	if (in_pet_store) {
		bool is_lit(0);

		switch (animal_type) {
		case TYPE_FISH  : is_lit = 1; break;
		case TYPE_RAT   : is_lit = 1; is_cage = ((shelf_ix + room_id) & 1); break; // should be a cage, and not lit, but a lit tank looks better visually
		case TYPE_SNAKE : is_lit = (rgen.rand_float() < 0.6); break;
		case TYPE_SPIDER: is_lit = (rgen.rand_float() < 0.4); break;
		case TYPE_BIRD  : is_cage = 1; break; // unlit
		default: assert(0); // unsupported
		}
		if (!is_cage) {flags |= RO_FLAG_ADJ_TOP;} // pet store tanks always have a lid
		if ( is_lit ) {flags |= RO_FLAG_LIT    ;}
		colorRGBA const cage_colors[3] = {WHITE, BLACK, BRASS_C};
		if (is_cage) {color = cage_colors[rgen.rand() % 3];} // select a random cage color
	}
	else if (rgen.rand_float() < 0.80) { // add a lid 80% of the time
		flags |= RO_FLAG_ADJ_TOP;
		if (rgen.rand_float() < 0.85) {flags |= RO_FLAG_LIT;} // light is on 85% of the time
	}
	objects.emplace_back(tank, (is_cage ? TYPE_PET_CAGE : TYPE_FISHTANK), room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, color, animal_type); // item_flags=animal_type
	set_obj_id(objects);
}

colorRGBA get_couch_color(rand_gen_t &rgen) {
	unsigned const NUM_COLORS = 8;
	colorRGBA const colors[NUM_COLORS] = {GRAY_BLACK, WHITE, LT_GRAY, GRAY, DK_GRAY, LT_BROWN, BROWN, DK_BROWN};
	return colors[rgen.rand()%NUM_COLORS];
}
bool building_t::add_livingroom_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!is_residential() || room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be living rooms
	float const wall_thickness(get_wall_thickness());
	cube_t place_area(get_walkable_room_bounds(room)), tv_pref_area;
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall
	vect_room_object_t &objs(interior->room_geom->objs);
	bool placed_couch(0), placed_tv(0);
	// place couches with a variety of colors
	colorRGBA const color(get_couch_color(rgen));
	unsigned tv_pref_orient(4), couch_ix(objs.size()), tv_ix(0);
	
	if (place_model_along_wall(OBJ_MODEL_COUCH, TYPE_COUCH, room, 0.40, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.67, 4, 1, color)) { // pref centered
		room_object_t const &couch(objs[couch_ix]);
		float const tv_overlap(0.5*couch.get_length());
		placed_couch   = 1;
		tv_pref_orient = (2*couch.dim + !couch.dir); // TV should be across from couch
		tv_pref_area   = place_area;
		max_eq(tv_pref_area.d[!couch.dim][0], couch.d[!couch.dim][0]-tv_overlap);
		min_eq(tv_pref_area.d[!couch.dim][1], couch.d[!couch.dim][1]+tv_overlap);
	}
	tv_ix = objs.size();

	// place TV: first attempt to place across from the couch, then prefer centered; maybe should set not_at_window=1, but that seems too restrictive
	if (place_model_along_wall(OBJ_MODEL_TV, TYPE_TV, room, 0.45, rgen, zval, room_id, tot_light_amt, tv_pref_area, objs_start, 4.0, tv_pref_orient, 1, BKGRAY, 0, 0, 1) ||
		place_model_along_wall(OBJ_MODEL_TV, TYPE_TV, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area,   objs_start, 4.0, tv_pref_orient, 1, BKGRAY, 0, 0, 0))
	{
		placed_tv = 1;
		// add a small table to place the TV on so that it's off the floor and not blocked as much by tables and chairs
		room_object_t &tv(objs[tv_ix]);
		float const height(0.4*tv.dz());
		cube_t table(tv); // same XY bounds as the TV
		tv.translate_dim(2, height); // move TV up
		table.z2() = tv.z1();
		objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, (RO_FLAG_IS_HOUSE | RO_FLAG_ADJ_TOP), tot_light_amt, SHAPE_SHORT); // short table; houses only
	}
	if (placed_couch && placed_tv) {
		room_object_t const &couch(objs[couch_ix]), &tv(objs[tv_ix]);

		if (couch.dim == tv.dim && couch.dir != tv.dir) { // placed against opposite walls facing each other
			cube_t region(couch);
			region.union_with_cube(tv);
			shorten_chairs_in_region(region, objs_start); // region represents that space between the couch and the TV
		}
	}
	if (!placed_couch && !placed_tv) return 0; // not a living room

	if (rgen.rand_bool()) { // add rocking chair
		unsigned const chair_ix(objs.size());
		cube_t chair_place_area(place_area);
		chair_place_area.expand_by(-wall_thickness); // move a bit further back from the wall to prevent intersections when rotating

		if (place_model_along_wall(OBJ_MODEL_RCHAIR, TYPE_RCHAIR, room, 0.5, rgen, zval, room_id, tot_light_amt, chair_place_area, objs_start, 1.0)) {
			if (rgen.rand_bool()) { // add a random rotation half the time
				assert(chair_ix < objs.size()); // chair must have been placed
				objs[chair_ix].flags |= RO_FLAG_RAND_ROT; // rotated, but generally facing the center of the room
				objs[chair_ix].shape  = SHAPE_CYLIN; // make it a cylinder since it no longer fits in a tight cube
			}
		}
	}
	if ((is_house || is_apartment()) && (rgen.rand()%3) == 0) { // add fishtank on a tall table 33% of the time
		add_fishtank_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area);
	}
	if (is_house && rgen.rand_float() < 0.75 && is_room_adjacent_to_ext_door(room, zval)) {add_shoes_by_door(rgen, room, zval, room_id, tot_light_amt, objs_start);}
	if (room.is_single_floor && objs_start > 0) {replace_light_with_ceiling_fan(rgen, room, cube_t(), room_id, tot_light_amt, objs_start-1);} // light is prev placed object
	return 1;
}

// Note: this room is decided by the caller and the failure to add objects doesn't make it not a dining room
void building_t::add_diningroom_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	//if (!is_house || room.is_hallway || room.is_sec_bldg || room.is_office) return; // still applies, but unnecessary
	if ((rgen.rand()&3) == 0) return; // no additional objects 25% of the time
	cube_t const room_bounds(get_room_bounds_inside_trim(room));
	float const vspace(get_window_vspace()), clearance(max(0.2f*vspace, get_min_front_clearance_inc_people()));
	vect_room_object_t &objs(interior->room_geom->objs);
	// add a wine rack
	float const width(0.3*vspace*rgen.rand_uniform(1.0, 1.5)), depth(0.16*vspace), height(0.4*vspace*rgen.rand_uniform(1.0, 1.5)); // depth is based on bottle length, which is constant
	cube_t c;
	set_cube_zvals(c, zval, zval+height);

	for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to place a wine rack; similar to placing a bookcase
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		c.d[dim][ dir] = room_bounds.d[dim][dir]; // against this wall
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		set_wall_width(c, pos, 0.5*width, !dim);
		cube_t tc(c);
		tc.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance; // increase space to add clearance
		if (is_obj_placement_blocked(tc, room, 1) || overlaps_other_room_obj(tc, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_WINE_RACK, room_id, dim, !dir, 0, tot_light_amt); // Note: dir faces into the room, not the wall
		set_obj_id(objs);
		break; // done/success
	} // for n
}

void gen_crate_sz(vector3d &sz, rand_gen_t &rgen, float window_vspacing) {
	for (unsigned d = 0; d < 3; ++d) {sz[d] = 0.06*window_vspacing*(1.0 + ((d == 2) ? 1.2 : 2.0)*rgen.rand_float());} // slightly more variation in XY
}

bool building_t::add_storage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement, bool has_stairs) {
	bool const is_garage_or_shed(room.is_garage_or_shed(0)), is_int_garage(room.get_room_type(0) == RTYPE_GARAGE);
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), floor_thickness(get_floor_thickness());
	float const ceil_zval(zval + window_vspacing - floor_thickness), shelf_depth((is_house ? (is_basement ? 0.18 : 0.15) : 0.2)*window_vspacing);
	float shelf_shorten(shelf_depth + 1.0f*wall_thickness);
	// increase shelf shorten for interior garages to account for approx width of exterior door when opened
	if (is_int_garage) {max_eq(shelf_shorten, 0.36f*window_vspacing);}
	cube_t room_bounds(get_walkable_room_bounds(room)), crate_bounds(room_bounds);
	vect_room_object_t &objs(interior->room_geom->objs);
	vector<unsigned> ds_ixs;
	vect_cube_t exclude;
	cube_t test_cube(room);
	set_cube_zvals(test_cube, zval, zval+wall_thickness); // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
	unsigned num_doors(0);

	// first pass to record the doors in this room
	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) { // check both dirs
		if ((i->no_room_conn() || i->get_bldg_conn() || i->is_connected_to_room(room_id)) && is_cube_close_to_door(test_cube, 0.0, 0, *i, 2)) {
			ds_ixs.push_back(i - interior->door_stacks.begin());
			++num_doors;
		}
	}
	if (room.has_tunnel_conn()) { // tunnel entrance counts as a door, but not false doors
		for (tunnel_seg_t const &t : interior->tunnels) {
			if (t.room_conn && t.conn_room_ix == room_id) {
				exclude.push_back(t.bcube_ext);
				exclude.back().expand_in_dim(!t.dim, 0.4*room.get_sz_dim(!t.dim)); // will cover the entire room if there's a door at the opposite end
				++num_doors;
			}
		}
	}
	for (unsigned dsix : ds_ixs) {
		door_stack_t const &ds(interior->door_stacks[dsix]);
		exclude.push_back(ds);
		exclude.back().expand_in_dim(ds.dim, 0.6*room.get_sz_dim(ds.dim));
		
		if (num_doors > 1) { // if there are multiple doors (houses only?), expand the exclude area more in the other dimension to make sure there's a path between doors
			exclude.back().expand_in_dim(!ds.dim, min(1.2f*ds.get_width(), 0.3f*room.get_sz_dim(!ds.dim)));
		}
		exclude.back().union_with_cube(ds.get_open_door_bcube_for_room(room)); // include open door
	}
	// add shelves on walls (avoiding any door(s)), and have crates avoid them
	bool const has_bike(is_int_garage && objs.size() >= 2 && objs[objs.size()-2].type == TYPE_GBIKE); // previously added bike, then blocker
	unsigned const rgen_mod(is_int_garage ? (has_bike ? 5 : 3) : 2); // 50% of the walls, 67% for interior garages, 80% if there's a bike

	for (unsigned dim = 0; dim < 2; ++dim) {
		if (room_bounds.get_sz_dim( dim) < 6.0*shelf_depth  ) continue; // too narrow to add shelves in this dim
		if (room_bounds.get_sz_dim(!dim) < 4.0*shelf_shorten) continue; // too narrow in the other dim

		for (unsigned dir = 0; dir < 2; ++dir) {
			if ((rgen.rand()%rgen_mod) == 0) continue; // only add shelves to some walls
			
			if (is_garage_or_shed) {
				// garage or shed - don't place shelves in front of door, but allow them against windows; basement - don't place against basement door
				cube_t room_wall(room);
				room_wall.d[dim][!dir] = room.d[dim][dir]; // shrink room to zero width along this wall
				if (is_room_adjacent_to_ext_door(room_wall)) continue; // no zval needed
			}
			else if (is_house && !is_basement && has_int_windows() && classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) {
				// don't place shelves against exterior house walls in case there are windows
				cube_t const &part(get_part_for_room(room));
				float const h_spacing(get_hspacing_for_part(part, !dim));
				if (room_bounds.get_sz_dim(!dim) - 2.0*shelf_depth > h_spacing) continue; // shelf width is larger than spacing - likely to intersect a window, don't test center pt
				if (is_val_inside_window(part, !dim, room_bounds.get_center_dim(!dim), h_spacing, get_window_h_border())) continue;
			}
			cube_t shelves(room_bounds);
			set_cube_zvals(shelves, zval, ceil_zval-floor_thickness);
			crate_bounds.d[dim][dir] = shelves.d[dim][!dir] = shelves.d[dim][dir] + (dir ? -1.0 : 1.0)*shelf_depth; // outer edge of shelves, which is also the crate bounds
			shelves.expand_in_dim(!dim, -shelf_shorten); // shorten shelves
			cube_t cands[3] = {shelves, shelves, shelves}; // full, lo half, hi half
			cands[1].d[!dim][1] = cands[2].d[!dim][0] = shelves.get_center_dim(!dim); // split in half
			unsigned const num_cands((cands[1].get_sz_dim(!dim) < 4.0*shelf_shorten) ? 1 : 3); // only check full length if short

			for (unsigned n = 0; n < num_cands; ++n) {
				cube_t const &cand(cands[n]);
				if (has_bcube_int(cand, exclude)) continue; // too close to a doorway
				if (!is_garage_or_shed && interior->is_blocked_by_stairs_or_elevator(cand)) continue;
				if (overlaps_other_room_obj(cand, objs_start)) continue; // can be blocked by bookcase, etc.
				bool const is_empty(rgen.rand_float() < 0.05); // 5% empty
				unsigned const shelf_flags((is_house ? RO_FLAG_IS_HOUSE : 0) | (is_garage_or_shed ? 0 : RO_FLAG_INTERIOR) | (is_empty ? 0 : RO_FLAG_NONEMPTY));
				add_shelves(cand, dim, dir, room_id, tot_light_amt, shelf_flags, 0, rgen); // item_flags=0
				break; // done
			} // for n
		} // for dir
	} // for dim
	if (is_garage_or_shed) return 1; // no chair, crates, or boxes in garages or sheds

	// add a random office chair if there's space
	if (!is_house && min(crate_bounds.dx(), crate_bounds.dy()) > 1.2*window_vspacing && has_office_chair_model()) {
		vector3d const chair_sz(get_office_chair_size());
		float const chair_height(chair_sz.z), chair_radius(chair_sz.x);
		point const pos(gen_xy_pos_in_area(crate_bounds, chair_radius, rgen, zval));
		cube_t chair(get_cube_height_radius(pos, chair_radius, chair_height));
		
		// for now, just make one random attempt; if it fails then there's no chair in this room
		if (!has_bcube_int(chair, exclude) && !is_obj_placement_blocked(chair, room, 1)) {
			objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_RAND_ROT, tot_light_amt, SHAPE_CYLIN, GRAY_BLACK);
		}
	}
	unsigned const num_max(is_house ? (is_basement ? 12 : 5) : 30); // 4-33 for offices, 4-8 for houses, 4-16 for house basements
	add_boxes_and_crates(rgen, room, zval, room_id, tot_light_amt, objs_start, num_max, is_basement, room_bounds, crate_bounds, exclude);
	// add a ladder leaning against the wall if storage room is on the ground floor or basement
	if (zval < ground_floor_z1 + window_vspacing && rgen.rand_bool()) {add_ladder_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);}
	// add office building storage room sign, in a hallway, basement, etc.
	if (!is_house) {add_door_sign((has_stairs ? "Stairs" : "Storage"), room, zval, room_id);}
	return 1; // it's always a storage room, even if it's empty
}

void check_for_blocked_box_flags(vect_room_object_t &objs, unsigned objs_start, unsigned obj_ix) {
	assert(objs_start <= obj_ix && obj_ix < objs.size());
	if (objs_start == obj_ix) return; // no other objects
	room_object_t &box(objs[obj_ix]);
	assert(box.type == TYPE_BOX);

	for (unsigned d = 0; d < 4; ++d) { // each side
		bool const dim(d>>1), dir(d&1);
		cube_t bc(box);
		set_cube_zvals(bc, box.z2(), (box.z2() + 0.5*box.dz())); // flaps extend above the box
		bc.d[dim][!dir]  = box.d[dim][dir]; // flush with the edge
		bc.d[dim][ dir] += (dir ? 1.0 : -1.0)*0.485*box.get_sz_dim(dim); // extend out by about half the box size (see building_room_geom_t::add_box())

		for (unsigned i = objs_start; i < obj_ix; ++i) { // check objects placed before the box
			if (objs[i].intersects(bc)) {box.flags |= (RO_FLAG_ADJ_LO << d); break;}
		}
	} // for d
}
void building_t::add_boxes_and_crates(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start,
	unsigned num_max, bool is_basement, cube_t const &room_bounds, cube_t const &crate_bounds, vect_cube_t const &exclude)
{
	float const window_vspacing(get_window_vspace()), ceil_zval(zval + get_floor_ceil_gap());
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const num_place(4 + (rgen.rand() % num_max)), boxes_start(objs.size());
	door_path_checker_t door_path_checker;
	unsigned num_placed(0);

	for (unsigned n = 0; n < 4*num_place; ++n) { // make up to 4 attempts for every crate/box
		vector3d sz; // half size relative to window_vspacing
		gen_crate_sz(sz, rgen, window_vspacing*(is_house ? (is_basement ? 0.75 : 0.5) : 1.0)); // smaller for houses
		if (crate_bounds.dx() <= 2.0*sz.x || crate_bounds.dy() <= 2.0*sz.y) continue; // too large for this room
		point const pos(gen_xy_pos_in_area(crate_bounds, sz, rgen, zval));
		cube_t crate(get_cube_height_radius(pos, sz, 2.0*sz.z)); // multiply by 2 since this is a size rather than half size/radius
		if (has_bcube_int(crate, exclude)) continue; // don't place crates between the door and the center of the room
		bool bad_placement(0);

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (!i->intersects(crate)) continue;
			// only handle stacking of crates on other crates
			if (i->is_crate_or_box() && i->z1() == zval && (i->z2() + crate.dz() < ceil_zval) && i->contains_pt_xy(pos)) {crate.translate_dim(2, i->dz());}
			else {bad_placement = 1; break;}
		}
		if (bad_placement) continue;
		if (is_obj_placement_blocked(crate, room, 1)) continue;
		if (door_path_checker.check_door_path_blocked(crate, room, room_id, zval, *this)) continue; // don't block the path between doors
		cube_t c2(crate);
		c2.expand_by(vector3d(0.5*c2.dx(), 0.5*c2.dy(), 0.0)); // approx extents of flaps if open
		room_object const type(rgen.rand_bool() ? TYPE_CRATE : TYPE_BOX);
		unsigned flags(0);

		if (type == TYPE_BOX) { // determine which sides are against a wall, for use with box flaps logic
			for (unsigned d = 0; d < 4; ++d) {
				bool const dim(d>>1), dir(d&1);
				if ((c2.d[dim][dir] < room_bounds.d[dim][dir]) ^ dir) {flags |= (RO_FLAG_ADJ_LO << d);}
			}
		}
		objs.emplace_back(crate, type, room_id, rgen.rand_bool(), 0, flags, tot_light_amt, SHAPE_CUBE, gen_box_color(rgen)); // crate or box
		set_obj_id(objs); // used to select texture and box contents
		if (++num_placed == num_place) break; // we're done
	} // for n
	for (unsigned i = boxes_start; i < objs.size(); ++i) { // check for other boxes or crates blocking box flaps
		if (objs[i].type == TYPE_BOX) {check_for_blocked_box_flags(objs, objs_start, i);}
	}
}

void building_t::add_shelves(cube_t const &c, bool dim, bool dir, unsigned room_id, float tot_light_amt, unsigned flags, unsigned item_flags, rand_gen_t &rgen) {
	assert(c.is_strictly_normalized());
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(c, TYPE_SHELVES, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, WHITE, item_flags);
	set_obj_id(objs);
	unsigned const objs_start(objs.size());
	interior->room_geom->expand_shelves(objs.back(), 1); // add_models_mode=1

	if (flags & RO_FLAG_IN_MALL) { // register fishtanks with non-fish animals
		assert(has_mall());

		for (auto c = objs.begin()+objs_start; c != objs.end(); ++c) {
			if (c->is_pet_container() && c->item_flags != TYPE_FISH) {interior->mall_info->pet_tanks.emplace_back(*c, c->item_flags, (c - objs.begin()));}
		}
	}
}

bool building_t::add_ladder_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	float const ladder_height(rgen.rand_uniform(0.77, 0.95)*get_floor_ceil_gap());
	vector3d const ladder_sz(0.25, rgen.rand_uniform(0.2, 0.22), 1.0); // D, W, H
	cube_t const place_area(get_walkable_room_bounds(room));
	return place_obj_along_wall(TYPE_INT_LADDER, room, ladder_height, ladder_sz, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.0, 1, 4, 0);
}

bool building_t::add_interrogation_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	// find the direction to the main (first) door
	bool dim(0), dir(0);
	unsigned const first_door_ix(has_ext_basement() ? interior->ext_basement_door_stack_ix : 0);
	
	for (auto ds = interior->door_stacks.begin()+first_door_ix; ds != interior->door_stacks.end(); ++ds) {
		door_t const &door(get_door(ds->first_door_ix)); // assume a single door stack
		if (!door.is_connected_to_room(room_id)) continue;
		dim = door.dim;
		dir = (room.get_center_dim(dim) < door.get_center_dim(dim));
		break; // found the door - done
	}
	// add chair(s) facing the door; they may have falled over
	unsigned const num_chairs((rgen.rand_float() < 0.25) ? 2 : 1); // add a pair of chairs 25% of the time
	colorRGBA const &chair_color(chair_colors[rgen.rand() % NUM_CHAIR_COLORS]);
	point const room_center(room.xc(), room.yc(), zval);
	point chair_pos[2] = {room_center, room_center};
	bool added_chair(0);

	if (num_chairs == 2) { // add two chairs: pairwise interrogation!
		float const chair_offset(0.2*get_window_vspace());
		chair_pos[0][!dim] -= chair_offset;
		chair_pos[1][!dim] += chair_offset;
	}
	for (unsigned n = 0; n < num_chairs; ++n) {
		if (add_chair(rgen, room, vect_cube_t(), room_id, chair_pos[n], chair_color, dim, dir, tot_light_amt, 0, 0, 0, 1)) { // no_push_out=1
			interior->room_geom->objs.back().flags |= RO_FLAG_PLCOLL; // flag so that player collisions are enabled
			added_chair = 1;
		}
	}
	if (!added_chair) return 0;
	// add bucket(s)
	unsigned const num_buckets((rgen.rand() % 3) + 1); // 1-3
	add_buckets_to_room(rgen, get_walkable_room_bounds(room), zval, room_id, tot_light_amt, objs_start, num_buckets);
	// add a ladder leaning against the wall
	if (rgen.rand_bool()) {add_ladder_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);}
	return 1;
}

void building_t::add_buckets_to_room(rand_gen_t &rgen, cube_t place_area, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num) {
	float const floor_spacing(get_window_vspace()), trim_thickness(get_trim_thickness());
	vect_room_object_t &objs(interior->room_geom->objs);
	place_area.z1() = zval + 1.5*get_flooring_thick(); // slightly above the flooring/rug to avoid z-fighting

	for (unsigned n = 0; n < num; ++n) {
		float const height(rgen.rand_uniform(0.1, 0.15)*floor_spacing), radius(rgen.rand_uniform(0.4, 0.6)*height);

		for (unsigned N = 0; N < 10; ++N) { // 10 place attempts
			cube_t const bucket(place_cylin_object(rgen, place_area, radius, height, (radius + trim_thickness), 1)); // place_at_z1=1
			if (is_obj_placement_blocked(bucket, place_area, 1) || overlaps_other_room_obj(bucket, objs_start)) continue; // bad placement
			objs.emplace_back(bucket, TYPE_BUCKET, room_id, rgen.rand_bool(), 0, 0, tot_light_amt, SHAPE_CYLIN, LT_GRAY); // dir=0 (unused)
			set_obj_id(objs); // used for setting water/liquid properties
			break; // success
		}
	} // for n
}

void building_t::add_garage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	vect_room_object_t &objs(interior->room_geom->objs);

	if (!enable_parked_cars() || (rgen.rand()&3) == 0) { // 75% of garages have cars
		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_BICYCLE) && rgen.rand_bool()) { // if no car, add a bike 50% of the time
			cube_t place_area(get_walkable_room_bounds(room));
			unsigned const bike_id(objs.size());
			
			if (place_model_along_wall(OBJ_MODEL_BICYCLE, TYPE_GBIKE, room, 0.45, rgen, zval, room_id, tot_light_amt,
				place_area, objs.size(), 1.0, 4, 0, WHITE, 0, 0, 0, 1)) // sideways
			{
				objs[bike_id].dir = rgen.rand_bool(); // random dir
			}
		}
		return;
	}
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_USED | RO_FLAG_INVIS); // lines not shown
	bool dim(0), dir(0); // set dir so that cars pull into driveways
	get_garage_dim_dir(room, dim, dir);
	cube_t space(room); // full room, car will be centered here
	set_cube_zvals(space, zval, (zval + get_rug_thickness()));
	room_object_t pspace(space, TYPE_PARK_SPACE, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, WHITE);
	pspace.obj_id = (uint16_t)(objs.size() + rgen.rand()); // will be used for the car model and color
	car_t const car(car_from_parking_space(pspace));
	interior->room_geom->wall_ps_start = objs.size(); // first parking space index
	cube_t collider(car.bcube);
	float const min_spacing(2.1*get_scaled_player_radius()); // space for the player to fit

	for (unsigned d = 0; d < 2; ++d) { // make sure there's enough spacing around the car for the player to walk without getting stuck
		max_eq(collider.d[d][0], (room.d[d][0] + min_spacing));
		min_eq(collider.d[d][1], (room.d[d][1] - min_spacing));
	}
	if (!collider.is_strictly_normalized()) {collider = car.bcube;} // garage is too small for player to fit; shouldn't happen
	objs.push_back(pspace);
	objs.emplace_back(collider, TYPE_COLLIDER, room_id, dim, dir, (RO_FLAG_INVIS | RO_FLAG_FOR_CAR));
	interior->room_geom->has_garage_car = 1;
}

void building_t::add_floor_clutter_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!is_residential()) return; // houses only for now

	if (rgen.rand_float() < 0.10) { // maybe add a toy 10% of the time
		vect_room_object_t &objs(interior->room_geom->objs);

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (i->type == TYPE_TOY) return; // don't place a toy on both a room object and on the floor
		}
		bool const use_model(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOY));
		float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness());
		cube_t place_area(get_walkable_room_bounds(room));
		place_area.expand_by(-1.0*wall_thickness); // add some extra padding
		float const height(0.11*window_vspacing), radius(height*(use_model ? get_radius_for_square_model(OBJ_MODEL_TOY) : 0.5f*0.67f));

		if (radius < 0.1*min(place_area.dx(), place_area.dy())) {
			point const pos(gen_xy_pos_in_area(place_area, radius, rgen, zval));
			cube_t c(get_cube_height_radius(pos, radius, height));

			// for now, just make one random attempt; if it fails then there's no chair in this room
			if (!overlaps_other_room_obj(c, objs_start) && !is_obj_placement_blocked(c, room, 1)) {
				if (use_model) { // symmetric, no dim or dir, but random rotation
					objs.emplace_back(c, TYPE_TOY_MODEL, room_id, 0, 0, (RO_FLAG_RAND_ROT | RO_FLAG_NOCOLL), tot_light_amt);
				}
				else { // random dim/dir; each ring is SHAPE_VERT_TORUS, but the overall toy is SHAPE_CYLIN
					objs.emplace_back(c, TYPE_TOY, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN);
					set_obj_id(objs); // used for color selection
				}
			}
		}
	}
}

cube_t place_cylin_object_maybe_near(rand_gen_t &rgen, cube_t const &place_on, point &prev_pos,
	float radius, float height, float dist_from_edge, float near_prob, float min_dist, float max_dist)
{
	if (prev_pos != all_zeros && rgen.rand_probability(near_prob)) {
		assert(min_dist < max_dist);
		cube_t c;

		for (unsigned n = 0; n < 10; ++n) { // 10 attempts to not intersect the original, otherwise fall through
			c.set_from_sphere((prev_pos + rgen.signed_rand_vector_spherical_xy(max_dist)), radius); // place at least spacing from edge
			if (dist_less_than(c.get_cube_center(), prev_pos, min_dist)) continue; // skip if too close
			if (!place_on.contains_cube_xy(c)) continue; // outside the place area
			set_cube_zvals(c, place_on.z1(), place_on.z1()+height);
			return c;
		}
	}
	cube_t const c(place_cylin_object(rgen, place_on, radius, height, dist_from_edge, 1)); // place_at_z1=1
	prev_pos = c.get_cube_center();
	return c;
}
void building_t::add_basement_clutter_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	bool const add_bottles(rgen.rand_float() < 0.35); // 35% of rooms
	bool const add_trash  (rgen.rand_float() < 0.35); // 35% of rooms
	bool const add_papers (rgen.rand_float() < 0.50); // 50% of rooms
	bool const add_glass  (rgen.rand_float() < 0.65); // 65% of the time
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-get_trim_thickness()); // add some extra padding
	add_floor_clutter_objs(rgen, room, place_area, zval, room_id, tot_light_amt, objs_start, add_bottles, add_trash, add_papers, add_glass);
}
void set_max_drink_can_type(room_object_t &obj, unsigned max_type) {
	if (obj.get_drink_can_type() > max_type) {obj.obj_id = max_type;} // clamp to max_type
}
void building_t::add_floor_clutter_objs(rand_gen_t &rgen, room_t const &room, cube_t place_area, float zval, unsigned room_id,
	float tot_light_amt, unsigned objs_start, bool add_bottles, bool add_trash, bool add_papers, bool add_glass)
{
	vect_room_object_t &objs(interior->room_geom->objs);
	float const floor_spacing(get_window_vspace()), min_place_sz(min(place_area.dx(), place_area.dy()));
	place_area.z1() = zval;
	vect_cube_t avoid;

	if (add_bottles) { // add bottles and cans on the floor
		unsigned const num_bottles((rgen.rand() % 12) + 1); // 1-12
		float const near_prob(0.75), near_dist(0.3*floor_spacing);
		bool const no_alcohol(is_school() && zval >= ground_floor_z1); // alcohol is only allowed in school basements
		unsigned const max_bottle_type(no_alcohol ? BOTTLE_TYPE_COKE : NUM_BOTTLE_TYPES-1);
		point prev_pos;

		for (unsigned n = 0; n < num_bottles; ++n) {
			float const radius(floor_spacing*rgen.rand_uniform(0.012, 0.018));
			if (min_place_sz < 6.0*radius) return; // room is too small to place this bottle; shouldn't get here
			bool const is_can(rgen.rand_float() < 0.33); // 33% cans, 67% bottles
			float const height(is_can ? 3.77*radius : floor_spacing*rgen.rand_uniform(0.075, 0.12)), stain_height(1.5*get_flooring_thick());
			cube_t bottle(place_cylin_object_maybe_near(rgen, place_area, prev_pos, radius, height, max(2.0f*radius, height), near_prob, 2.0*height, near_dist));
			cube_t bc(bottle);
			bool const fallen_over(rgen.rand_float() < 0.75); // make the bottle fallen over 75% of the time
			unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_ON_FLOOR);
			bool dim(0), dir(0);

			if (fallen_over) {
				dim = rgen.rand_bool();
				dir = rgen.rand_bool();
				float const delta(height - 2.0*radius), xy_shift((dir ? 1.0 : -1.0)*delta);
				bottle.z2() -= delta;
				bottle.d[dim][dir] += xy_shift;
				bc.translate_dim(dim, 0.5*xy_shift); // recenter
				bc.expand_by_xy(delta); // account for any rotation
				if (!place_area.contains_cube_xy(bc)) continue;
				flags |= RO_FLAG_RAND_ROT; // rotate into a random orientation
			}
			if (is_obj_placement_blocked(bc, room, 1) || overlaps_other_room_obj(bc, objs_start) || has_bcube_int(bc, avoid)) continue; // bad placement
			avoid.push_back(bc);
			bool const add_stain(fallen_over && rgen.rand_bool());
			colorRGBA stain_color(ALPHA0); // maybe add a stain the same color as the liquid in the bottle

			if (is_can) { // can
				room_object_t obj(bottle, TYPE_DRINK_CAN, room_id, dim, dir, flags, tot_light_amt, SHAPE_CYLIN);
				obj.obj_id = rgen.rand();
				if (no_alcohol) {set_max_drink_can_type(obj, DRINK_CAN_TYPE_COKE);}
				obj.obj_id |= BOTTLE_EMPTY_MASK; // make empty
				if (add_stain ) {stain_color = drink_can_params[obj.get_drink_can_type()].liquid_color;}
				objs.emplace_back(obj);
			}
			else { // bottle
				room_object_t obj(bottle, TYPE_BOTTLE, room_id, dim, dir, flags, tot_light_amt, SHAPE_CYLIN);
				obj.set_as_bottle(rgen.rand(), max_bottle_type, 0, 0, 1, 1); // make_empty=1, allow_transparent=1
				if (add_stain) {stain_color = bottle_params[obj.get_bottle_type()].liquid_color;}
				objs.emplace_back(obj);
			}
			if (stain_color.alpha > 0.0) { // not transparent (water)
				float const stain_radius(10.0*radius*rgen.rand_uniform(0.5, 1.0));
				// stain should be near the open end of the bottle
				point const center(bottle.xc(), bottle.yc(), zval+stain_height);
				vector3d rot_dir(objs.back().get_dir());
				rotate_vector3d(plus_z, objs.back().get_bottle_rot_angle(), rot_dir);
				interior->room_geom->decal_manager.add_blood_or_stain((center - (0.5*height)*rot_dir), stain_radius, stain_color, 0); // is_blood=0
			}
		} // for n
	}
	if (add_trash) { // add trash (paper balls) on the floor
		unsigned const num_trash((rgen.rand() % 6) + 1); // 1-6
		float const near_prob(0.6), near_dist(0.3*floor_spacing);
		point prev_pos;

		for (unsigned n = 0; n < num_trash; ++n) {
			float const radius(rgen.rand_uniform(0.015, 0.03)*min(floor_spacing, min_place_sz));
			cube_t trash(place_cylin_object_maybe_near(rgen, place_area, prev_pos, radius, 2.0*radius, 1.5*radius, near_prob, 2.0*radius, near_dist));
			if (is_obj_placement_blocked(trash, room, 1) || overlaps_other_room_obj(trash, objs_start) || has_bcube_int(trash, avoid)) continue; // bad placement
			avoid.push_back(trash);
			colorRGBA const color(trash_colors[rgen.rand() % NUM_TRASH_COLORS]);
			objs.emplace_back(trash, TYPE_TRASH, room_id, rgen.rand_bool(), rgen.rand_bool(), (RO_FLAG_NOCOLL | RO_FLAG_ON_FLOOR), tot_light_amt, SHAPE_SPHERE, color);
			set_obj_id(objs);
		} // for n
	}
	if (add_papers) { // add sheets of paper on the floor
		// similar to add_papers_to_surface()
		unsigned const num_papers((rgen.rand() % 5) + 1); // 1-5
		unsigned const paper_flags(RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT | RO_FLAG_ON_FLOOR);
		float const plen(0.115*floor_spacing), pwidth(0.77*plen), thickness(0.00025*floor_spacing); // 8.5x11
		point prev_pos;
		cube_t paper;
		set_cube_zvals(paper, zval, (zval + get_rug_thickness() + thickness)); // very thin

		for (unsigned n = 0; n < num_papers; ++n) { // okay if they overlap
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
			set_wall_width(paper, rgen.rand_uniform(place_area.d[ dim][0]+plen,   place_area.d[ dim][1]-plen  ), 0.5*plen,    dim);
			set_wall_width(paper, rgen.rand_uniform(place_area.d[!dim][0]+pwidth, place_area.d[!dim][1]-pwidth), 0.5*pwidth, !dim);
			if (is_obj_placement_blocked(paper, room, 1) || overlaps_other_room_obj(paper, objs_start) || has_bcube_int(paper, avoid)) continue; // bad placement
			objs.emplace_back(paper, TYPE_PAPER, room_id, dim, dir, paper_flags, tot_light_amt, SHAPE_CUBE, select_paper_color(rgen), btype);
			set_obj_id(objs);
			paper.z2() += thickness; // to avoid Z-fighting if different colors
		} // for n
	}
	if (add_glass) { // add broken glass on the floor; often fails to be placed
		float const glass_zval(zval + 1.1*get_flooring_thick()); // slightly above the flooring/rug to avoid z-fighting
		float const radius(rgen.rand_uniform(0.08, 0.12)*min(floor_spacing, min_place_sz));
		cube_t const bc(place_cylin_object(rgen, place_area, radius, 0.1*radius, 1.5*radius, 1)); // place_at_z1=1

		if (!is_obj_placement_blocked(bc, room, 1) && !overlaps_other_room_obj(bc, objs_start) && !has_bcube_int(bc, avoid)) {
			avoid.push_back(bc);
			add_broken_glass_decal(point(bc.xc(), bc.yc(), glass_zval), radius, rgen);
		}
	}
}

void building_t::add_laundry_basket(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, cube_t place_area) {
	float const floor_spacing(get_window_vspace()), radius(rgen.rand_uniform(0.1, 0.12)*floor_spacing), height(rgen.rand_uniform(1.5, 2.25)*radius);
	place_area.expand_by_xy(-radius); // leave a slight gap between laundry basket and wall
	if (!place_area.is_strictly_normalized()) return; // no space for laundry basket (likely can't happen)
	cube_t legal_area(get_part_for_room(room));
	legal_area.expand_by_xy(-(1.0*floor_spacing + radius)); // keep away from part edge/exterior walls to avoid alpha mask drawing problems (unless we use mats_amask)
	point center;
	center.z = zval + 0.002*floor_spacing; // slightly above the floor to avoid z-fighting

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a laundry basket
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		center[ dim] = place_area.d[dim][dir]; // against this wall
		center[!dim] = rgen.rand_uniform(place_area.d[!dim][0], place_area.d[!dim][1]);
		if (!legal_area.contains_pt_xy(center)) continue; // too close to part edge
		cube_t const c(get_cube_height_radius(center, radius, height));
		if (is_obj_placement_blocked(c, room, !room.is_hallway) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
		colorRGBA const colors[4] = {WHITE, LT_BLUE, LT_GREEN, LT_BROWN};
		room_obj_shape const shape(rgen.rand_bool() ? SHAPE_CUBE : SHAPE_CYLIN);
		interior->room_geom->objs.emplace_back(c, TYPE_LBASKET, room_id, dim, dir, 0, tot_light_amt, shape, colors[rgen.rand()%4]);
		set_obj_id(interior->room_geom->objs); // for clothing type and orient
		break; // done
	} // for n
}

bool building_t::add_laundry_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned &added_bathroom_objs_mask) {
	float const front_clearance(get_min_front_clearance_inc_people());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // common spacing to wall for appliances
	vector2d const place_area_sz(place_area.get_size_xy());
	vect_room_object_t &objs(interior->room_geom->objs);
	bool success(0);

	// check for bookcase and add blocker, since it's possible to have both a bookcase and laundry objects placed in the same room
	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (i->type != TYPE_BCASE) continue;
		cube_t blocker(*i);
		blocker.d[i->dim][i->dir] += (i->dir ? 1.0 : -1.0)*2.0*i->get_depth();
		objs.emplace_back(blocker, TYPE_BLOCKER, room_id, i->dim, i->dir, RO_FLAG_INVIS); // add blocker cube
		break; // there should only be one, and we've invalidated the reference
	}
	if (!is_house) { // commercial laundry room
		unsigned const num_floors(round_fp(room.dz()/get_window_vspace())), max_place(num_floors*interior->rooms.size()/8);
		bool const wd_first(rgen.rand_bool());
		bool any_placed(0);

		for (unsigned e = 0; e < 2; ++e) {
			bool const wd(bool(e) ^ wd_first); // {washer, dryer}
			unsigned const obj_ix(objs.size()), model_id(wd ? OBJ_MODEL_DRYER : OBJ_MODEL_WASHER);
			room_object const type(wd ? TYPE_DRYER : TYPE_WASHER);
			float const height(wd ? 0.38 : 0.42);
			if (!place_model_along_wall(model_id, type, room, height, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8)) continue;
			// try to place more objects of this type along the same wall
			room_object_t const obj(objs[obj_ix]); // deep copy
			float const shift_amt(rgen.rand_uniform(1.25, 2.0)*obj.get_width()), clearance(0.8*obj.get_depth());
			unsigned num_place(1);

			for (unsigned dir = 0; dir < 2; ++dir) { // {left, right} of obj
				room_object_t obj2(obj);

				for (unsigned n = 1; n < max_place; ++n) {
					obj2.translate_dim(!obj.dim, (dir ? 1.0 : -1.0)*shift_amt);
					if (!room.contains_cube(obj2)) break; // done
					cube_t blocker(obj2); // used for collision tests
					blocker.d[obj.dim][obj.dir] += (obj.dir ? 1.0 : -1.0)*clearance; // add front clearance
					if (overlaps_obj_or_placement_blocked(blocker, room, objs_start)) continue;
					objs.push_back(obj2);
					objs.emplace_back(blocker, TYPE_BLOCKER, room_id, obj.dim, obj.dir, RO_FLAG_INVIS);
					++num_place;
				} // for n
			} // for dir
			// if no extra objects were placed, try a different wall/pos
			if (num_place == 1) {place_model_along_wall(model_id, type, room, height, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8);}
			any_placed = 1;
		} // for wed
		if (!any_placed) return 0; // failed
		add_door_sign("Laundry", room, zval, room_id);
		success = 1;
	}
	else { // house laundry room
		for (unsigned n = 0; n < 10; ++n) { // 10 attempts to place washer and dryer along the same wall
			unsigned const washer_ix(objs.size());
			bool const placed_washer(place_model_along_wall(OBJ_MODEL_WASHER, TYPE_WASHER, room, 0.42, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8));
			unsigned pref_orient(4); // if washer was placed, prefer to place dryer along the same wall
			if (placed_washer) {pref_orient = objs[washer_ix].get_orient();}
			unsigned const dryer_ix(objs.size());
			bool const placed_dryer(place_model_along_wall(OBJ_MODEL_DRYER, TYPE_DRYER, room, 0.38, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8, pref_orient));
			if (placed_dryer) {objs[dryer_ix].obj_id = (uint16_t)objs.size();} // set dryer obj_id; will determine whether this is gas or electric

			if (placed_washer && placed_dryer && objs[dryer_ix].get_orient() == pref_orient) {success = 1;} // placed both washer and dryer along the same wall
			else if (n+1 == 10) { // last attempt
				if (!(placed_washer || placed_dryer)) return 0; // placed neither washer nor dryer, failed
				if (placed_washer != placed_dryer) {success = 1;} // placed only one of the washer or dryer, allow it
				else if (objs[washer_ix].dim != objs[dryer_ix].dim) {success = 1;} // placed on two adjacent walls, allow it
				// placed on opposite walls; check that there's space for the player to walk between the washer and dryer
				else if (objs[washer_ix].get_sz_dim(objs[washer_ix].dim) + objs[dryer_ix].get_sz_dim(objs[dryer_ix].dim) +
					front_clearance < place_area_sz[objs[washer_ix].dim]) {success = 1;}
			}
			if (success) break;
			objs.resize(objs_start); // remove washer and dryer and try again
		} // for n
	}
	if (!success) return 0; // failed

	// if we've placed at least washer and/or dryer and made this into a laundry room, try to place a sink as well; should this use a different sink model from bathrooms?
	if (place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6)) {
		added_bathroom_objs_mask |= PLACED_SINK;
	}
	add_laundry_basket(rgen, room, zval, room_id, tot_light_amt, objs_start, place_area); // try to place a laundry basket
	return 1; // done
}

void building_t::add_entryway_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	// what about apartment entryway rooms next to entrance interior doors?
	if (is_house && rgen.rand_float() < 0.8) {add_shoes_by_door(rgen, room, zval, room_id, tot_light_amt, objs_start);}
}
void building_t::add_shoes_by_door(rand_gen_t &rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_SHOE)) return;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SHOE)); // L, W, H
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	float const length(0.2*floor_spacing*rgen.rand_uniform(0.75, 1.0)), width(length*sz.y/sz.x), height(length*sz.z/sz.x), hlen(0.5*length), pair_hw(width);
	unsigned const num_shoes(1 + (rgen.rand() % 3)); // 1-3
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t room_exp(room);
	room_exp.expand_by_xy(wall_thickness);

	for (tquad_with_ix_t const &d : doors) { // find the front exterior door
		if (!d.is_exterior_door() || d.type == tquad_with_ix_t::TYPE_RDOOR) continue;
		cube_t dbc(d.get_bcube());
		if (!room_exp.contains_pt(dbc.get_cube_center())) continue;
		bool const dim(dbc.dy() < dbc.dx()); // door dim
		cube_t door_path(dbc);
		door_path.expand_in_dim( dim, 1.0*floor_spacing ); // in front of door
		door_path.expand_in_dim(!dim, 0.5*wall_thickness);
		dbc.expand_in_dim( dim, 0.3*floor_spacing); // near the wall the door is on
		dbc.expand_in_dim(!dim, 0.7*floor_spacing);
		cube_t place_area(get_walkable_room_bounds(room));
		place_area.expand_by(-0.5*wall_thickness);
		place_area.intersect_with_cube(dbc);
		vector2d const place_sz(place_area.get_size_xy());
		cube_t shoes;
		set_cube_zvals(shoes, zval, zval+height);
		bool const dir(rgen.rand_bool());

		for (unsigned n = 0, num_added = 0; n < 20 && num_added < num_shoes; ++n) { // try to place shoes
			if (place_sz[dim] < 1.2*length || place_sz[!dim] < 2.4*pair_hw) continue; // not enough space
			set_wall_width(shoes, rgen.rand_uniform((place_area.d[ dim][0] + hlen   ), (place_area.d[ dim][1] - hlen   )), hlen,     dim);
			set_wall_width(shoes, rgen.rand_uniform((place_area.d[!dim][0] + pair_hw), (place_area.d[!dim][1] - pair_hw)), pair_hw, !dim);
			if (shoes.intersects_xy(door_path) || overlaps_other_room_obj(shoes, objs_start)) continue;
			if (interior->is_cube_close_to_doorway(shoes, room, 0.0, 1, 0) || interior->is_blocked_by_stairs_or_elevator(shoes)) continue; // not checking ext door
			unsigned const item_flags(rgen.rand()); // random shoe sub-model
			add_obj_pair(room_object_t(shoes, TYPE_SHOE, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, WHITE, item_flags), objs);
			++num_added;
		} // for n
		break; // only add to the front door
	} // for d
}

void building_t::add_couches_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned const counts[4]) {
	unsigned const num_couches(counts[rgen.rand() & 3]);
	if (num_couches == 0) return;
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // common spacing to wall
	colorRGBA const color(get_couch_color(rgen)); // make the colors match if there's more than one couch

	for (unsigned n = 0; n < num_couches; ++n) {
		place_model_along_wall(OBJ_MODEL_COUCH, TYPE_COUCH, room, 0.40, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 4, 0, color);
	}
}

cube_t get_pool_table_top_surface(room_object_t const &c) {
	cube_t top(c);
	top.expand_by_xy(-0.12*c.get_width());
	top.z1() = c.z2() - 0.046*c.dz();
	return top;
}
bool check_dist_and_add_center(point const &pos, float diameter_xy, vector<point> &centers) {
	for (point const &p : centers) {
		if (dist_xy_less_than(p, pos, diameter_xy)) return 0;
	}
	centers.push_back(pos);
	return 1;
}
// room with pool table, not swimming pool, that one is below
bool building_t::add_pool_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	float const floor_spacing(get_window_vspace()), sz_in_feet(floor_spacing/8.0), clearance(get_min_front_clearance_inc_people());
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_start(objs.size()), pool_table_obj_ix(objs_start);
	bool const long_dim(room.dx() < room.dy());
	point const pos(room.xc(), room.yc(), zval); // place in the center of the room
	cube_t ptable(pos);
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_POOL_TABLE)) return 0;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_POOL_TABLE)); // L, W, H
	// size is 8ft long x 4ft wide x 2.5 ft tall, with 4.5 feet of side clearance
	float const length(8.0*sz_in_feet), width(length*sz.y/sz.x), height(length*sz.z/sz.x);
	ptable.z2() += height;
	ptable.expand_in_dim( long_dim, 0.5*length);
	ptable.expand_in_dim(!long_dim, 0.5*width);
	cube_t outer_bcube(ptable), pad_bcube(ptable);
	outer_bcube.expand_by_xy(max(clearance, 4.5f*sz_in_feet)); // add spacing for pool cues
	if (!room.contains_cube_xy(outer_bcube)) return 0; // can't fit in this room
	pad_bcube.expand_by_xy(clearance);
	if (is_obj_placement_blocked(pad_bcube, room, 1)) return 0; // inc_open_doors=1
	bool const table_dir(rgen.rand_bool()); // dir doesn't much matter, though the ends of the table are slightly different
	room_object_t const pool_table(ptable, TYPE_POOL_TABLE, room_id, long_dim, table_dir, 0, tot_light_amt, SHAPE_CUBE);
	objs.push_back(pool_table);

	// place pool balls
	float const ball_radius(0.0117*length), ball_diameter(2.0*ball_radius); // pool ball diameter is 2.25", 1.125" radius; length is 8', so radius is ~0.0117*length
	cube_t top(get_pool_table_top_surface(pool_table));
	top.z2() = top.z1() + ball_diameter;
	float const ball_zval(top.zc());
	unsigned const num_balls(16); // including cue ball
	vector<point> centers;

	for (unsigned n = 0; n < num_balls; ++n) {
		for (unsigned N = 0; N < 10; ++N) { // 10 tries to find a valid pos; if we can't find a pos, the ball won't be placed
			point const pos(gen_xy_pos_in_area(top, ball_radius, rgen, ball_zval));
			if (!check_dist_and_add_center(pos, ball_diameter, centers)) continue;
			centers.push_back(pos);
			cube_t ball;
			ball.set_from_sphere(pos, ball_radius);
			objs.emplace_back(ball, TYPE_POOL_BALL, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_DSTATE), tot_light_amt, SHAPE_SPHERE);
			room_object_t &pball(objs.back());
			pball.item_flags  = n; // assign ball number
			pball.obj_id      = (uint16_t)interior->room_geom->allocate_dynamic_state(); // allocate a new dynamic state object
			pball.state_flags = (uint16_t)pool_table_obj_ix; // encode object index in state_flags
			get_pool_ball_rot_matrix(pball, interior->room_geom->get_dstate(pball).rot_matrix);
			break; // success
		} // for n
	} // for n

	// place wall mount and pool cues; these can be on the wall, on the table, or on the floor leaning against the table
	float const cue_len(57*sz_in_feet/12), cue_diameter(1.0*sz_in_feet/12), cue_radius(0.5*cue_diameter); // 57 inches x 1 inch diameter
	bool const cue_dir(rgen.rand_bool());
	cube_t cue;
	// attach wall mount for pool cues
	vector3d const fc_sz_scale(0.15, 3.6, 1.0); // depth, width, height
	cube_t const place_area(get_walkable_room_bounds(room)); // right against the wall
	float const mount_zval(zval + 0.7*cue_len);
	unsigned const wall_mount_obj_ix(objs.size());
	unsigned num_table_cues(0), num_leaning_cues(0);

	if (place_obj_along_wall(TYPE_WALL_MOUNT, room, 0.08*cue_len, fc_sz_scale, rgen, mount_zval, room_id, tot_light_amt, place_area, objs_start, 10.0, 0, 4, 1, WHITE)) {
		assert(objs.size() == wall_mount_obj_ix+2); // must have added wall mount + blocker
		objs.back().z1() = zval; // extend blocker down to the floor to prevent a couch from being placed here
		unsigned const num_cues(4);
		room_object_t const mount(objs[wall_mount_obj_ix]);
		bool const cdim(mount.dim), cdir(mount.dir);
		point pos(0.0, 0.0, mount.zc());
		pos[cdim] = mount.d[cdim][cdir] + (cdir ? 1.0 : -1.0)*cue_radius;
		float const spacing(mount.get_sz_dim(!cdim)/(num_cues + 1)), block_hwidth(0.4*cue_radius);

		for (unsigned n = 0; n < num_cues; ++n) {
			pos[!cdim] = mount.d[!cdim][0] + (n+1)*spacing;
			cue.set_from_point(pos);
			cue.expand_by_xy(cue_radius);
			cue.expand_in_z(0.5*cue_len); // set height/length
			// either add the cue to the holder, on the table, or leaning against the table
			if (rgen.rand_bool()) {++((num_table_cues < 2 && rgen.rand_bool()) ? num_table_cues : num_leaning_cues);}
			else {objs.emplace_back(cue, TYPE_POOL_CUE, room_id, 0, 1, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, WHITE);} // dim=2 (calculated from size)

			for (unsigned d = 0; d < 2; ++d) { // add tiny blocks on each side of the cue, even if not present
				room_object_t block(mount);
				block.expand_in_z(-0.2*mount.dz()); // shrink in Z
				for (unsigned e = 0; e < 2; ++e) {block.d[cdim][e] = cue.d[cdim][e];}
				set_wall_width(block, (pos[!cdim] + (d ? 1.0 : -1.0)*(0.7*cue_radius + block_hwidth)), block_hwidth, !cdim);
				objs.push_back(block);
			}
		} // for n
	}
	else { // failed to place wall mount
		(rgen.rand_bool() ? num_table_cues : num_leaning_cues) = 2; // two cues, either on the table or leaning against it
	}
	if (num_table_cues > 0) { // on the top edges of the pool table
		set_cube_zvals(cue, ptable.z2(), (ptable.z2() + 2.0*cue_radius));
		set_wall_width(cue, ptable.get_center_dim(long_dim), 0.5*cue_len, long_dim);
		bool const first_side(rgen.rand_bool());

		for (unsigned n = 0; n < num_table_cues; ++n) { // for each long side; can be at most 2
			bool const d((n == 1) ^ first_side);
			set_wall_width(cue, (ptable.d[!long_dim][d] + (d ? -1.0 : 1.0)*0.05*width), cue_radius, !long_dim);
			objs.emplace_back(cue, TYPE_POOL_CUE, room_id, long_dim, cue_dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, WHITE);
		}
	}
	if (num_leaning_cues > 0) { // on the floor leaning up against the pool table
		bool const cdim(rgen.rand_bool()); // all on the same or different sides
		centers.clear();

		for (unsigned n = 0; n < num_leaning_cues; ++n) {
			bool const cdir(rgen.rand_bool()); // side to place the pool cues on
			point pos(0.0, 0.0, zval);
			pos[ cdim] = ptable.d[cdim][cdir] + (cdir ? 1.0 : -1.0)*cue_radius;
			pos[!cdim] = rgen.rand_uniform((ptable.d[!cdim][0] + 8.0*cue_radius), (ptable.d[!cdim][1] - 8.0*cue_radius));
			if (!check_dist_and_add_center(pos, cue_diameter, centers)) continue;
			cue.set_from_point(pos);
			cue.expand_by_xy(cue_radius);
			cue.z2() += cue_len;
			objs.emplace_back(cue, TYPE_POOL_CUE, room_id, 0, 1, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, WHITE); // dim=2 (calculated from size)
		} // for n
	}
	
	// maybe place couch(es) along a wall
	unsigned const counts[4] = {0, 1, 1, 2}; // one couch is more common
	add_couches_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, counts);
	// place a mini bar?
	// place two bar stools
	unsigned const bs_id(objs.size());

	if (place_model_along_wall(OBJ_MODEL_BAR_STOOL, TYPE_BAR_STOOL, room, 0.4, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.0, 4, 0)) {
		// if a bar stool was placed, try to place another next to it
		room_object_t const &bs(objs[bs_id]);
		float const translate_dist(rgen.rand_uniform(1.2, 2.0)*bs.get_sz_dim(!bs.dim));
		bool const pref_dir(rgen.rand_bool());
		bool placed(0);

		for (unsigned d = 0; d < 2 && !placed; ++d) {
			bool const dir(bool(d) ^ pref_dir);
			room_object_t bs2(bs);
			bs2.translate_dim(!bs.dim, (dir ? 1.0 : -1.0)*translate_dist);
			if (!place_area.contains_cube_xy(bs2) || overlaps_other_room_obj(bs2, objs_start) || is_obj_placement_blocked(bs2, room, 1)) continue; // inc_open_doors=1
			objs.push_back(bs2);
			placed = 1;
		} // for d
		if (!placed) { // if not placed next to the other bar stool, try placing anywhere along the same wall as the first one
			place_model_along_wall(OBJ_MODEL_BAR_STOOL, TYPE_BAR_STOOL, room, 0.4, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.0, bs.get_orient(), 0);
		}
	}
	return 1;
}

// for indoor pools, currently only in extended basements
void building_t::add_swimming_pool_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	assert(has_pool());
	cube_t const place_area(get_walkable_room_bounds(room));
	indoor_pool_t &pool(interior->pool);
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), tile_thickness(get_flooring_thick());
	float const pool_len(pool.get_sz_dim(pool.dim)), pool_depth(pool.dz());
	float const shallow_depth(0.5*floor_spacing); // about 4 feet
	vect_room_object_t &objs(interior->room_geom->objs);
	
	// add tile around the sides of the pool
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			cube_t side(pool);
			side.z2() += tile_thickness; // fill the gap with the tiles on the floor above
			side.d[dim][!dir]  = pool.d[dim][dir];
			side.d[dim][ dir] += (dir ? 1.0 : -1.0)*tile_thickness;
			expand_to_nonzero_area(side, tile_thickness, dim);
			objs.emplace_back(side, TYPE_POOL_TILE, room_id, dim, dir, RO_FLAG_ADJ_LO); // flag RO_FLAG_ADJ_LO for being inside the pool
		} // for dir
	} // for dim
	// add tile to the bottom of the pool; make it darker for bottomless pools
	cube_t bottom(pool);
	bottom.z2() = pool.z1() + tile_thickness;
	objs.emplace_back(bottom, TYPE_POOL_TILE, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_LO), 1.0, SHAPE_CUBE, (pool.bottomless ? GRAY : WHITE));
	// add ceiling and floor tile
	cube_t ceil(room);
	ceil.z2() = room.z2() - get_fc_thickness();
	ceil.z1() = ceil.z2() - tile_thickness;
	objs.emplace_back(ceil, TYPE_POOL_TILE, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_TOP));
	zval += tile_thickness; // any objects will be placed on the floor tile
	cube_t floor_test(room), drain;
	floor_test.expand_by(-get_wall_thickness());
	cube_t deep_end(pool); // for balls and floats resting on the dry pool bottom
	deep_end.d[pool.dim][pool.dir] = pool.get_center_dim(pool.dim); // start with the half away from the stairs

	for (cube_t const &f : interior->floors) {
		if (!floor_test.intersects(f)) continue;
		cube_t fc(f);
		set_cube_zvals(fc, f.z2(), (f.z2() + tile_thickness));
		objs.emplace_back(fc, TYPE_POOL_TILE, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_BOT));
	}
	// add wall tile
	for (unsigned d = 0; d < 2; ++d) {
		for (cube_t const &wall : interior->walls[d]) {
			if (!wall.intersects_no_adj(room)) continue;
			bool const dir(wall.get_center_dim(d) < room.get_center_dim(d)); // direction facing the center of the room
			cube_t tile(wall);
			tile.d[d][!dir]  = wall.d[d][dir]; // edge of wall facing the room
			tile.d[d][ dir] += (dir ? 1.0 : -1.0)*tile_thickness;
			expand_to_nonzero_area(tile, tile_thickness, d);
			objs.emplace_back(tile, TYPE_POOL_TILE, room_id, d, !dir, RO_FLAG_NOCOLL);
		} // for wall
	} // for d
	// add a sloped ramp at the bottom if deep enough and not bottomless
	if (!pool.bottomless && pool_depth > 1.2*shallow_depth && pool_len > 3.0*floor_spacing) {
		interior->room_geom->pool_ramp_obj_ix = objs.size();
		cube_t slope(pool), upper(pool);
		slope.expand_in_dim(pool.dim, -0.25*pool_len); // shrink to middle 50% (cut 25% off of each end)
		pool.shallow_zval = slope.z2() = upper.z2() = pool.z1() + shallow_depth; // shallow end
		upper   .d[pool.dim][!pool.dir] = slope.d[pool.dim][ pool.dir];
		deep_end.d[pool.dim][ pool.dir] = slope.d[pool.dim][!pool.dir];
		objs.emplace_back(slope, TYPE_POOL_TILE, room_id, pool.dim, pool.dir, (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_BOT), 1.0, SHAPE_ANGLED);
		objs.emplace_back(upper, TYPE_POOL_TILE, room_id, pool.dim, pool.dir, (RO_FLAG_ADJ_LO | RO_FLAG_ADJ_BOT)); // bottom of shallow end
	}
	if (!pool.bottomless) { // add a pool drain at the deep end
		float const drain_radius(min(0.025*pool_len, 0.08*floor_spacing));
		point drain_center;
		drain_center[!pool.dim] = pool.get_center_dim(!pool.dim);
		drain_center[ pool.dim] = pool.d[pool.dim][!pool.dir] + (pool.dir ? 1.0 : -1.0)*5.0*drain_radius; // at deep end
		drain.set_from_sphere(drain_center, drain_radius);
		set_cube_zvals(drain, pool.z1(), pool.z1()+0.08*drain_radius);
		objs.emplace_back(drain, TYPE_DRAIN, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_IN_POOL), tot_light_amt, SHAPE_CYLIN, LT_GRAY);
	}
	unsigned const objs_start(objs.size()); // we can start here, since the pool tile objects placed above don't collide with other placed objects
	cube_t ladder;

	if (pool_len > 5.0*floor_spacing && pool_depth > 0.9*floor_spacing && rgen.rand_float() < 0.75) { // maybe add a diving board if long and deep
		cube_t dboard;
		set_cube_zvals(dboard, zval, (zval + 0.11*floor_spacing));
		set_wall_width(dboard, pool.get_center_dim(!pool.dim), 0.1*floor_spacing, !pool.dim); // width
		set_wall_width(dboard, (pool.d[pool.dim][!pool.dir] + (pool.dir ? 1.0 : -1.0)*0.2*floor_spacing), 0.6*floor_spacing,  pool.dim); // length
		objs.emplace_back(dboard, TYPE_DIV_BOARD, room_id, pool.dim, pool.dir, 0, tot_light_amt, SHAPE_CUBE, colorRGBA(0.75, 1.0, 0.9, 1.0)); // blue-green
	}
	if (pool_len > 3.5*floor_spacing && pool_depth > 0.8*floor_spacing && rgen.rand_float() < 0.75 && building_obj_model_loader.is_model_valid(OBJ_MODEL_POOL_LAD)) {
		// add a side ladder if long and deep enough
		bool const ldir(rgen.rand_bool());
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_POOL_LAD)); // D, W, H
		float const height(0.75*floor_spacing), hdepth(0.5*height*sz.x/sz.z), hwidth(0.5*height*sz.y/sz.z);
		set_cube_zvals(ladder, (pool.z2() - 0.42*height), (pool.z2() + 0.58*height));
		set_wall_width(ladder, (pool.d[ pool.dim][pool.dir] + (pool.dir ? -1.0 : 1.0)*rgen.rand_uniform(0.8, 0.9)*pool_len), hdepth, pool.dim); // pos along pool length
		set_wall_width(ladder, (pool.d[!pool.dim][ldir] + (ldir ? 1.0 : -1.0)*0.16*hwidth), hwidth, !pool.dim); // at pool edge
		objs.emplace_back(ladder, TYPE_POOL_LAD, room_id, !pool.dim, !ldir, RO_FLAG_NOCOLL, tot_light_amt);
	}
	// add benches along the walls
	unsigned const num_benches(2 + (rgen.rand()%3)); // 2-4

	if (num_benches > 0) {
		float const bench_height(0.22*floor_spacing);
		vector3d const bench_sz(rgen.rand_uniform(0.7, 0.9), rgen.rand_uniform(4.0, 6.0), 1.0); // D, W, H
		bool pref_dir(rgen.rand_bool()); // random first dir

		for (unsigned n = 0; n < num_benches; ++n) {
			unsigned const pref_orient(2*(!pool.dim) + pref_dir), bench_obj_id(objs.size());
			if (!place_obj_along_wall(TYPE_BENCH, room, bench_height, bench_sz, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, 1, pref_orient)) break;
			maybe_place_obj_on_bench(objs[bench_obj_id], rgen, 0.3);
			pref_dir ^= 1; // alternate sides
		}
	}
	if (1) { // add a clock on the wall by the door or opposite the door
		bool const digital(rgen.rand_bool());
		float const place_pos(room.get_center_dim(!pool.dim)), clock_z1(max((zval + 0.6*floor_spacing), (room.z2() - 0.4*floor_spacing))); // near the ceiling in tall rooms
		float const clock_height((digital ? 0.12 : 0.2)*floor_spacing), clock_width((digital ? 4.0 : 1.0)*clock_height), clock_depth(0.08*clock_width);
		cube_t clock;
		set_cube_zvals(clock, clock_z1, clock_z1+clock_height);
		set_wall_width(clock, place_pos, 0.5*clock_width, !pool.dim);

		for (unsigned d = 0; d < 2; ++d) {
			bool const dir(pool.dir ^ bool(d));
			float const wall_pos(place_area.d[pool.dim][dir]);
			clock.d[pool.dim][ dir] = wall_pos;
			clock.d[pool.dim][!dir] = wall_pos + (dir ? -1.0 : 1.0)*clock_depth;
			if (is_cube_close_to_doorway(clock, room, wall_thickness, 1)) continue; // inc_open=1
			add_clock(clock, room_id, tot_light_amt, pool.dim, !dir, digital);
			break;
		} // for d
	}
	// add other objects in and around the pool;
	// note that stairs haven't been placed yet, so we shouldn't place objects near the future stairs when there's no water
	bool const pool_has_water(has_water());
	cube_t pool_area(pool_has_water ? pool : deep_end); // limit to deep end when there's no water since balls and floats slide down the ramp and should avoid stairs
	set_cube_zvals(pool_area, (pool_has_water ? interior->water_zval : pool.z1()), pool.z2());
	if (pool_has_water) {pool_area.d[pool.dim][pool.dir] += (pool.dir ? -1.0 : 1.0)*0.6*floor_spacing;} // avoid the railings by the stairs (approximate)
	
	if (pool_len > 3.0*floor_spacing && pool_depth > 0.8*floor_spacing) { // add pool float(s) if pool is large and deep enough
		unsigned const num_floats(rgen.rand() % 3); // 0-2
		float const pf_radius(0.25*floor_spacing), pf_height(0.6*pf_radius);

		for (unsigned n = 0; n < num_floats; ++n) {
			if (min(pool_area.dx(), pool_area.dy()) < 3.0*pf_radius) break; // pool place area is too small
			cube_t pfloat;
			gen_xy_pos_for_round_obj(pfloat, pool_area, pf_radius, pf_height, pf_radius, rgen, 1); // place_at_z1=1
			
			if (!pool_has_water) { // resting on the bottom of the pool if there's no water
				point const center(cube_bot_center(pfloat));
				float float_zval(center.z);
				if (get_zval_for_pool_bottom(center, float_zval)) {pfloat.translate_dim(2, (float_zval - center.z));}
			}
			if (overlaps_other_room_obj(pfloat, objs_start))         continue; // if placement falls, leave it out; should only collide with another float
			if (!ladder.is_all_zeros() && pfloat.intersects(ladder)) continue; // check the pool ladder
			if (!pool_has_water && pfloat.intersects(drain))         continue;
			bool const deflated(!pool_has_water);
			unsigned flags(0);
			if (deflated) {flags |= RO_FLAG_BROKEN; pfloat.z2() -= 0.9*pfloat.dz();}
			colorRGBA const &color(pfloat_colors[rgen.rand() % NUM_PFLOAT_COLORS]);
			objs.emplace_back(pfloat, TYPE_POOL_FLOAT, room_id, 0, 0, flags, tot_light_amt, (deflated ? SHAPE_CYLIN : SHAPE_VERT_TORUS), color);
		} // for n
	}
	if (pool_len > 2.0*floor_spacing) { // add beach ball(s) if pool is large enough
		unsigned const num_balls(rgen.rand() % 3); // 0-2

		for (unsigned n = 0; n < num_balls; ++n) {
			bool const in_pool(rgen.rand_bool());
			cube_t const &avoid_xy (in_pool ? ladder    : pool      );
			cube_t const &ball_area(in_pool ? pool_area : place_area);
			add_ball_to_room(rgen, room, ball_area, zval, room_id, tot_light_amt, objs_start, BALL_TYPE_BEACH, avoid_xy, in_pool);
		}
	}
}

// add objects to rooms connected to walkways, since the walkways themselves aren't actual rooms (and aren't inside the building bcube);
// returns true if there's a real door connecting the walkway with this room
bool building_t::maybe_add_walkway_room_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, light_ix_assign_t &light_ix_assign) {
	assert(has_room_geom());
	float const wall_thickness(get_wall_thickness()), ceil_zval(zval + get_floor_ceil_gap()), floor_spacing(get_window_vspace()), door_thickness(0.25*wall_thickness);
	vect_room_object_t &objs(interior->room_geom->objs);
	bool added(0);

	for (building_walkway_t &w : walkways) {
		bool const dim(w.dim);
		cube_t w_ext(w.bcube);
		w_ext.expand_in_dim(dim, door_thickness); // this will also set the door thickness
		cube_t test_cube(w_ext);
		if (!room.intersects_xy(test_cube))                  continue; // walkway not connected to this room
		if (zval < w.bcube.z1() || ceil_zval > w.bcube.z2()) continue; // wrong floor
		if (w.conn_bldg != nullptr && w.conn_bldg->get_window_vspace() != floor_spacing) continue; // floors not aligned (shouldn't happen?)
		float const door_width(get_office_ext_doorway_width());
		bool const dir(room.get_center_dim(dim) < w.bcube.get_center_dim(dim)), has_door(w.has_ext_door(!dir));
		if (has_door && !is_room_adjacent_to_ext_door(room, zval)) continue; // not the room connected to the walkway
		// update bcube_inc_rooms to include this room's length, so that lights and trim will be drawn when the player is in this room
		// (and possibly other rooms bordering the walkway); applies to both rooms, even if there is no door
		w.bcube_inc_rooms.union_with_cube_xy(room);

		if (w.conn_bldg != nullptr) { // update our neighbor's copy of the walkway as well
			for (building_walkway_t &w2 : w.conn_bldg->walkways) {
				if (w2.bcube == w.bcube) {w2.bcube_inc_rooms.union_with_cube_xy(room); break;}
			}
		}
		if (ADD_WALKWAY_EXT_DOORS && w.is_owner) { // add ceiling light(s) if interior may be accessible (even if no door) and we're the owner
			float const length(w.get_length());

			if (length > 0.4*floor_spacing) { // not too short for a light; should always get here
				unsigned const num_lights(max(1U, unsigned(0.6*length/floor_spacing)));
				float const light_spacing(length/(num_lights + 1));
				unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_EXTERIOR | RO_FLAG_LIT);
				colorRGBA const color(1.0, 1.0, 1.0); // white
				point lpos;
				lpos[ dim] = w.bcube.d[dim][0] + light_spacing;
				lpos[!dim] = w.bcube.get_center_dim(!dim);
				lpos.z     = ceil_zval - 0.02*floor_spacing;

				for (unsigned n = 0; n < num_lights; ++n) {
					cube_t light(lpos);
					light.z1() -= 0.02*floor_spacing;
					light.expand_in_dim( dim, 0.12*floor_spacing);
					light.expand_in_dim(!dim, 0.08*floor_spacing);
					objs.emplace_back(light, TYPE_LIGHT, room_id, dim, 0, flags, 0.0, SHAPE_CUBE, color); // dir=0 (unused); not actually in the room
					objs.back().obj_id = light_ix_assign.get_next_ix();
					lpos[dim] += light_spacing;
				} // for n
			}
		}
		if (has_door) { // has real exterior door
			float const center(w.bcube.get_center_dim(!w.dim));
			if (room.d[!dim][1] < center - 0.5*door_width || room.d[!dim][0] > center + 0.5*door_width) continue; // not overlapping the door (which is centered on the walkway)
			added = 1; // using real exterior doors; don't add false doors; blockers are unnecessary
			continue;
		}
		unsigned const floor_ix((zval - w.bcube.z1())/floor_spacing), door_mask(1 << floor_ix); // within the walkway
		assert(floor_ix < 8);
		if (!room.is_hallway && (w.has_door & door_mask)) continue; // already has a door on this floor; allow multiple doors if this is a hallway
		cube_t entry(w_ext);
		entry.intersect_with_cube(room);
		max_eq(entry.z1(), zval);
		min_eq(entry.z2(), ceil_zval);
		float const entry_width(entry.get_sz_dim(!dim));
		if (entry_width < 1.2*door_width) continue; // too narrow for a doorway; walkway likely ends at a wall separating rooms
		cube_t door(entry);
		door.expand_in_dim(!dim, 0.5*(door_width - entry_width)); // shrink to doorway width
		cube_t blocker(door);
		blocker.d[dim][!dir] += (dir ? -1.0 : 1.0)*1.2*door_width; // add space in the front for the door to open
		
		// skip if blocked by stairs or elevator; maybe the walkway shouldn't be placed here? walkways are added after placing stairs and elevators
		if (interior->is_blocked_by_stairs_or_elevator(blocker)) {
			if (!ADD_WALKWAY_EXT_DOORS) continue; // choose a different door pos
			// if we get here, that means we couldn't add a real door and we can't add a two sided fake door;
			// but a walkway that leads to nowhere is no good, so create a single sided door that exists only in the walkway
			door.d[dim][dir] -= (dir ? -1.0 : 1.0)*(wall_thickness + 0.25*door_thickness); // extend into walkway, enough to prevent Z-fighting
			objs.emplace_back(door, TYPE_FALSE_DOOR, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_WALKWAY | RO_FLAG_INTERIOR), tot_light_amt);
		}
		else {
			door.d[dim][dir] -= (dir ? -1.0 : 1.0)*(wall_thickness + door_thickness); // extend into walkway on the other side of this room so that it's visible from walkway
			objs.emplace_back(door, TYPE_FALSE_DOOR, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_WALKWAY), tot_light_amt);
			objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, dir, RO_FLAG_INVIS); // prevent objects from blocking the door from inside the building
		}
		w.has_door |= door_mask;
	} // for w
	return added;
}

bool get_fire_ext_height_and_radius(float window_vspacing, float &height, float &radius) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_FIRE_EXT)) return 0;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FIRE_EXT)); // D, W, H
	height = 0.24*window_vspacing;
	radius = height*(0.25*(sz.x + sz.y)/sz.z);
	return 1;
}
bool building_t::add_fire_ext_along_wall(cube_t const &room, float zval, unsigned room_id, float tot_light_amt, bool dim, bool dir, rand_gen_t &rgen) {
	float const window_vspacing(get_window_vspace());
	float fe_height(0.0), fe_radius(0.0);
	if (!get_fire_ext_height_and_radius(window_vspacing, fe_height, fe_radius)) return 0;
	float const min_clearance(2.0*fe_radius), wall_pos_lo(room.d[!dim][0] + min_clearance), wall_pos_hi(room.d[!dim][1] - min_clearance);
	if (wall_pos_lo >= wall_pos_hi) return 0; // should never fail?
	float const wall_pos(room.d[dim][dir] + (dir ? -1.0 : 1.0)*0.5*get_wall_thickness());

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts at placing a fire extinguisher
		float const val(rgen.rand_uniform(wall_pos_lo, wall_pos_hi)), cov_lo(val - min_clearance), cov_hi(val + min_clearance);

		if (is_contained_in_wall_range(wall_pos, cov_lo, cov_hi, (zval + 0.5*window_vspacing), dim)) { // shouldn't need to check anything else?
			add_fire_ext(fe_height, fe_radius, zval, wall_pos, val, room_id, tot_light_amt, dim, dir);
			return 1; // done/success
		}
	} // for n
	return 0; // failed to place
}
void building_t::add_fire_ext(float height, float radius, float zval, float wall_edge, float pos_along_wall,
	unsigned room_id, float tot_light_amt, bool dim, bool dir, bool center_mount)
{
	float const window_vspacing(get_window_vspace()), dir_sign(dir ? -1.0 : 1.0), xlate(((dim ^ dir ^ 1) ? 1.0 : -1.0)*0.32*radius);
	point pos(0.0, 0.0, (zval + 0.32*window_vspacing)); // bottom position
	pos[ dim] = wall_edge + dir_sign*1.3*radius; // radius away from the wall + plate thickness
	pos[!dim] = pos_along_wall + (center_mount ? -xlate : 0.0);
	vect_room_object_t &objs(interior->room_geom->objs);
	// add fire extinguisher
	cube_t const fe_bcube(get_cube_height_radius(pos, radius, height));
	objs.emplace_back(fe_bcube, TYPE_FIRE_EXT, room_id, !dim, (dir ^ dim), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN); // mounted sideways
	// add the wall mounting bracket; what about adding a small box with a door that contains the fire extinguisher?
	cube_t wall_mount(fe_bcube);
	wall_mount.expand_in_dim(!dim, -0.3*radius);
	wall_mount.translate_dim(!dim, xlate); // shift to line up with FE body
	wall_mount.d[dim][ dir]  = wall_edge; // extend to touch the wall
	wall_mount.d[dim][!dir] -= dir_sign*0.7*radius; // move inward
	wall_mount.z1() -= 0.02*height; // under the fire extinguisher
	wall_mount.z2() -= 0.50*height;
	objs.emplace_back(wall_mount, TYPE_FEXT_MOUNT, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, GRAY_BLACK);
	// add the sign
	cube_t sign;
	sign.d[dim][ dir] = wall_edge; // extend to touch the wall
	sign.d[dim][!dir] = wall_edge + dir_sign*0.05*radius;
	set_cube_zvals(sign, (zval + 0.65*window_vspacing), (zval + 0.80*window_vspacing));
	set_wall_width(sign, wall_mount.get_center_dim(!dim), 0.5*radius, !dim); // line up with wall bracket
	objs.emplace_back(sign, TYPE_FEXT_SIGN, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);
}

bool building_t::is_contained_in_wall_range(float wall_pos, float cov_lo, float cov_hi, float zval, bool dim) const {
	for (cube_t const &wall : interior->walls[dim]) {
		if (wall.d[ dim][0] > wall_pos || wall.d[ dim][1] < wall_pos) continue; // not on the correct side of this hallway
		if (wall.d[!dim][0] > cov_lo   || wall.d[!dim][1] < cov_hi  ) continue; // range not covered
		if (wall.z1() > zval || wall.z2() < zval) continue; // wrong zval/floor
		return 1;
	}
	return 0;
}
cube_t get_rand_reception_desk_placement_cube(room_object_t const &rdesk, rand_gen_t &rgen, unsigned &sel_ix) {
	unsigned const placed(rdesk.item_flags & 7);
	if (placed == 7) return cube_t(); // all three cubes are occupied
	cube_t cubes[3]; // {front, left, right}
	get_reception_desk_cubes(rdesk, cubes);
	if (placed == 6) return cubes[0]; // only front is available
	if (placed == 5) return cubes[1]; // only left  is available
	if (placed == 3) return cubes[2]; // only right is available
	unsigned rval_to_ix[4] = {0, 0, 1, 2}; // front cube is twice as likely
	if (placed == 1) {rval_to_ix[0] = 1; rval_to_ix[1] = 2;} // front is occupied, redistribute weights to back
	if (placed == 2) {rval_to_ix[2] = 2;} // left  is occupied, redistribute weight to right
	if (placed == 4) {rval_to_ix[3] = 1;} // right is occupied, redistribute weight to left
	sel_ix = rval_to_ix[rgen.rand()&3];
	return cubes[sel_ix];
}
void building_t::add_pri_hall_objs(rand_gen_t rgen, rand_gen_t room_rgen, room_t const &room, float zval,
	unsigned room_id, float tot_light_amt, unsigned floor_ix, unsigned objs_start)
{
	bool const long_dim(room.dx() < room.dy()), first_floor(floor_ix == 0);
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), hall_width(room.get_sz_dim(!long_dim));
	float const hall_side_clearance(max(get_doorway_width(), get_min_front_clearance_inc_people())); // front clearance is generally slightly larger
	vect_room_object_t &objs(interior->room_geom->objs);

	// place ground floor/lobby objects + hospital per-floor desks; not for schools
	if ((first_floor && room.z1() == ground_floor_z1 && !is_school()) || is_hospital()) {
		// add lobby reception desks
		float const nom_desk_width(0.9*window_vspacing);
		
		if (hall_width > (nom_desk_width + 1.6*hall_side_clearance)) { // hallway is wide enough for a reception desk
			float const desk_width(min(nom_desk_width, 0.5f*(hall_width - hall_side_clearance))); // shrink to make sure there's clearance on each side
			float const centerline(room.get_center_dim(!long_dim)), desk_depth(0.6*nom_desk_width), hall_len(room.get_sz_dim(long_dim));
			cube_t desk;
			set_cube_zvals(desk, zval, zval+0.32*window_vspacing);
			set_wall_width(desk, centerline, 0.5*desk_width, !long_dim);

			for (unsigned dir = 0; dir < 2; ++dir) { // add a reception desk at each entrance
				float const hall_start(room.d[long_dim][dir]), dir_sign(dir ? -1.0 : 1.0);

				if (first_floor) { // use a range of reasonable desk placements along the hall
					float const val1(hall_start + max(0.1f*hall_len, window_vspacing)*dir_sign), val2(hall_start + 0.3*hall_len*dir_sign);

					for (unsigned n = 0; n < 10; ++n) { // try to find a valid placement near the door; make 10 random attempts
						float const val(rgen.rand_uniform(min(val1, val2), max(val1, val2)));
						set_wall_width(desk, val, 0.5*desk_depth, long_dim);
						if (add_reception_desk(rgen, desk, long_dim, dir, room_id, tot_light_amt)) break;
					}
				}
				else { // place near the far end of the hallway, since there's no door to block
					if (dir == ((interior->rooms.size() & 1) ? 1 : 0)) continue; // place only at one end of the hallway, consistent per floor
					set_wall_width(desk, (hall_start + 0.8*window_vspacing*dir_sign), 0.5*desk_depth, long_dim);
					add_reception_desk(rgen, desk, long_dim, !dir, room_id, tot_light_amt); // place facing the stairs/elevator
				}
			} // for dir
		} // end reception desk
	}
	// maybe add a clock and US flag on the back of the stairs, if this is the lobby or the floor directly above the retail area
	if (first_floor && is_ground_floor_excluding_retail(room.z1()) && room.has_stairs) {
		for (stairwell_t const &s : interior->stairwells) {
			if (!s.is_u_shape() && s.shape != SHAPE_WALLED_SIDES) continue;
			if (s.extends_below && !s.is_u_shape())               continue; // skip stairs extending down to the basement/parking garage
			if (s.z1() > zval || !room.contains_cube(s))          continue; // stairs not on ground floor, or not contained
			cube_t place_cube(s);
			place_cube.d[s.dim][s.dir] += (s.dir ? 1.0 : -1.0)*0.3*wall_thickness; // account for stairs wall (approximate)
			bool const digital(rgen.rand_bool());
			add_clock_to_cube(place_cube, zval, room_id, tot_light_amt, s.dim, s.dir, digital);

			if (rgen.rand_float() < 0.75) { // place US flag
				// Note: flags are two sided, so lighting doesn't look correct on the unlit side
				bool const side(s.dir ^ s.dim); // side of clock; always place on the right because the left side is lit with correct normals
				float const flag_pos(s.d[!s.dim][side] + (side ? -1.0 : 1.0)*max(0.1*window_vspacing, 0.1*s.get_width())), wall_pos(place_cube.d[s.dim][s.dir]);
				add_wall_us_flag(wall_pos, flag_pos, zval, s.dim, s.dir, room_id, tot_light_amt);
			}
		} // for sb
	}
	// add a fire extinguisher on the wall
	bool const pri_dir(room_rgen.rand_bool()); // random, but the same across all floors
	add_fire_ext_along_wall(room, zval, room_id, tot_light_amt, !long_dim, pri_dir, room_rgen);

	// add a water fountain to the center of the hall (likely between stairs and elevator)
	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_WFOUNTAIN)) {
		float used_space_center(0.0); // calculate occupied hallway width near water fountain, assuming the first stairwell and elevator are in the primary hallway
		if (!interior->stairwells.empty()) {max_eq(used_space_center, interior->stairwells.front().get_width());}
		if (!interior->elevators .empty()) {max_eq(used_space_center, interior->elevators .front().get_width());}
		bool const dim(!long_dim), dir(!pri_dir); // opposite of fire extinguisher so as not to overlap
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_WFOUNTAIN)); // D, W, H
		float const height(0.25*window_vspacing), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z);
		float const side_space(0.5*(hall_width - used_space_center));

		if (side_space - depth > hall_side_clearance) { // add if there's enough clearance for people to pass by
			float const z1(zval+0.18*window_vspacing), wall_pos(room.d[dim][dir] + (dir ? -1.0 : 1.0)*0.5*wall_thickness);
			float const hall_center(room.get_center_dim(!dim)), pref_offset((room_rgen.rand_bool() ? 1.0 : -1.0)*0.8*window_vspacing);
			float const test_offsets[3] = {0.0, pref_offset, -pref_offset};

			for (unsigned n = 0; n < 3; ++n) { // try center plus one offset in each direction (in case there's a central door or hallway)
				cube_t wf;
				set_wall_width(wf, (hall_center + test_offsets[n]), hwidth, !dim);

				if (is_contained_in_wall_range(wall_pos, wf.d[!dim][0], wf.d[!dim][1], z1, dim)) { // shouldn't need to check anything else?
					set_cube_zvals(wf, z1, z1+height);
					wf.d[dim][ dir] = wall_pos;
					wf.d[dim][!dir] = wall_pos - (dir ? 1.0 : -1.0)*depth;
					objs.emplace_back(wf, TYPE_WFOUNTAIN, room_id, dim, dir, 0, tot_light_amt, SHAPE_CUBE);
					break;
				}
			} // for n
		}
	}
	// add trashcans
	bool const both_ends(floor_ix == 0); // floors above the ground floor have only one trashcan
	add_corner_trashcans(rgen, room, zval, room_id, tot_light_amt, long_dim, both_ends);
	// add cameras to each end of the hallway
	add_cameras_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start);
}

void building_t::add_corner_trashcans(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, bool dim, bool both_ends) {
	// add one or more trashcan to a corner; use a large cylinder, the same as a mall trashcan
	float const window_vspacing(get_window_vspace());
	float const tcan_height(0.26*window_vspacing), tcan_radius(0.1*window_vspacing), wall_spacing(1.2*tcan_radius);
	cube_t tcan;
	set_cube_zvals(tcan, zval, zval+tcan_height); // set height
	bool add_to_ends[2] = {1, 1}; // defaults to both ends
	if (!both_ends) {add_to_ends[rgen.rand_bool()] = 0;} // add only one trashcan

	for (unsigned end_dir = 0; end_dir < 2; ++end_dir) {
		if (!add_to_ends[end_dir]) continue;
		bool const side_dir(rgen.rand_bool());
		float const side_pos(room.d[!dim][side_dir] + (side_dir ? -1.0 : 1.0)*wall_spacing);
		float const end_pos (room.d[ dim][end_dir ] + (end_dir  ? -1.0 : 1.0)*wall_spacing);
		set_wall_width(tcan, side_pos, tcan_radius, !dim);
		set_wall_width(tcan, end_pos,  tcan_radius,  dim);
		if (is_obj_placement_blocked(tcan, room, 1)) continue; // inc_open=1
		room_object_t const tcan_obj(tcan, TYPE_TCAN, room_id, dim, end_dir, RO_FLAG_IN_HALLWAY, tot_light_amt, SHAPE_CYLIN, LT_GRAY);
		interior->room_geom->objs.push_back(tcan_obj);
		add_large_trashcan_contents(rgen, tcan_obj, room_id, tot_light_amt);
	} // for d
}

void building_t::add_wall_us_flag(float wall_pos, float flag_pos, float zval, bool dim, bool dir, unsigned room_id, float tot_light_amt) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_FLAG)) return;
	// Note: flags are two sided, so lighting doesn't look correct on the unlit side
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FLAG)); // W, D, H
	float const window_vspacing(get_window_vspace());
	float const flag_height(0.3*window_vspacing), flag_hwidth(0.5*flag_height*sz.x/sz.z), flag_depth(flag_height*sz.y/sz.z);
	cube_t flag;
	flag.z1() = zval + 0.53*window_vspacing;
	flag.z2() = flag.z1() + flag_height;
	flag.d[dim][!dir] = wall_pos;
	flag.d[dim][ dir] = wall_pos + (dir ? 1.0 : -1.0)*flag_depth;
	set_wall_width(flag, flag_pos, flag_hwidth, !dim);
	interior->room_geom->objs.emplace_back(flag, TYPE_US_FLAG, room_id, !dim, (dim ^ dir), RO_FLAG_NOCOLL, tot_light_amt);
}

bool building_t::add_reception_desk(rand_gen_t &rgen, cube_t const &desk, bool dim, bool dir, unsigned room_id, float tot_light_amt) {
	if (interior->is_blocked_by_stairs_or_elevator(desk)) return 0; // bad location, try a new one
	vect_room_object_t &objs(interior->room_geom->objs);

	if (has_office_chair_model()) { // add an office chair behind the desk
		vector3d const chair_sz(get_office_chair_size());
		point pos;
		pos.z = desk.z1();
		pos[!dim] = desk.get_center_dim(!dim);
		pos[ dim] = desk.get_center_dim( dim) + (dir ? -1.0 : 1.0)*(-0.05*desk.get_sz_dim(dim) + chair_sz.x); // push the chair into the cutout of the desk
		cube_t const chair(get_cube_height_radius(pos, chair_sz.x, chair_sz.z));
		if (interior->is_blocked_by_stairs_or_elevator(chair)) return 0; // bad location, try a new one
		objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, dim, dir, 0, tot_light_amt, SHAPE_CYLIN, GRAY_BLACK);
	}
	unsigned const rdesk_ix(objs.size());
	objs.emplace_back(desk, TYPE_RDESK, room_id, dim, dir, 0, tot_light_amt, SHAPE_CUBE);

	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_PHONE) && rgen.rand_float() < 0.75) { // add phone 75% of the time
		unsigned sel_ix(0);
		cube_t const place_area(get_rand_reception_desk_placement_cube(objs.back(), rgen, sel_ix));

		if (place_area.is_strictly_normalized()) {
			bool pdim(0), pdir(0); // facing the chair
			if      (sel_ix == 0) {pdim =  dim; pdir = !dir;} // front
			else if (sel_ix == 1) {pdim = !dim; pdir = 1;} // left
			else if (sel_ix == 2) {pdim = !dim; pdir = 0;} // right

			if (place_phone_on_obj(rgen, place_area, room_id, tot_light_amt, pdim, pdir, 0.2)) { // allow a bit of overhang
				objs[rdesk_ix].item_flags |= (1 << sel_ix); // mark this surface as occupied
			}
		}
	}
	return 1;
}

void building_t::add_cameras_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	// add cameras at each end of the room in the long dim (for hallways, retail, etc.)
	bool const long_dim(room.dx() < room.dy());
	float const window_vspacing(get_window_vspace()), doorway_width(get_doorway_width());
	float const ceil_zval(room.is_single_floor ? (room.z2() - get_fc_thickness()) : (zval + get_floor_ceil_gap())); // on upper floor of retail
	float const length(0.09*window_vspacing), width(0.4*length), height(0.5*length);
	cube_t const place_area(get_walkable_room_bounds(room));
	bool const camera_side(rgen.rand_bool()), is_ground_floor(zval >= ground_floor_z1 && zval < ground_floor_z1 + window_vspacing);
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t camera;
	set_cube_zvals(camera, (ceil_zval - height), ceil_zval);

	for (unsigned dir = 0; dir < 2; ++dir) {
		float const wall_pos(place_area.d[long_dim][!dir]), signed_len((dir ? 1.0 : -1.0)*length);
		float pos(room.get_center_dim(!long_dim));
		bool offset(is_ground_floor);
		
		if (!walkways.empty() && !offset) { // check for walkway doors
			float const check_radius(width + get_wall_thickness());
			point test_pt;
			test_pt[ long_dim] = wall_pos;
			test_pt[!long_dim] = pos;
			test_pt.z = 0.5*(zval + ceil_zval);

			for (auto const &door : doors) {
				if (door.get_bcube().contains_pt_exp(test_pt, check_radius)) {offset = 1; break;}
			}
		}
		if (offset) {pos += 0.65*doorway_width*((bool(dir) ^ camera_side) ? 1.0 : -1.0);} // place off to the side to avoid blocking doorway and exit sign
		set_wall_width(camera, pos, 0.5*width, !long_dim);
		camera.d[long_dim][!dir] = wall_pos + 0.25*signed_len; // near the wall
		camera.d[long_dim][ dir] = wall_pos + 1.25*signed_len; // extend away from the wall
		bool blocked(0);

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) { // check if blocked by a walkway false door, etc.
			if (i->type == TYPE_BLOCKER && i->intersects(camera)) {blocked = 1; break;}
		}
		if (!blocked) {objs.emplace_back(camera, TYPE_CAMERA, room_id, long_dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, WHITE);}
		// Note: if we want to have cameras on the bottom retail floor, they need brackets extending from the wall
	} // for dir
}

bool building_t::add_server_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) { // for office buildings
	float const window_vspacing(get_window_vspace());
	float const server_height(0.7*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const server_width (0.3*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hwidth(0.5*server_width);
	float const server_depth (0.4*window_vspacing*rgen.rand_uniform(0.9, 1.1)), server_hdepth(0.5*server_depth);
	float const comp_height  (0.2*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const min_spacing  (0.1*window_vspacing*rgen.rand_uniform(0.9, 1.1));
	float const comp_hwidth(0.5*0.44*comp_height), comp_hdepth(0.5*0.9*comp_height); // fixed AR=0.44 to match the texture
	float const server_period(server_width + min_spacing);
	bool const long_dim(room.dx() < room.dy());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // server spacing from walls
	zval = add_flooring(room, zval, room_id, tot_light_amt, FLOORING_CONCRETE); // add concrete and move the effective floor up
	cube_t server, computer;
	set_cube_zvals(server,   zval, (zval + server_height));
	set_cube_zvals(computer, zval, (zval + comp_height  ));
	point center;
	unsigned num_servers(0), num_comps(0);
	vect_room_object_t &objs(interior->room_geom->objs);

	// try to line servers up against each wall wherever they fit
	for (unsigned D = 0; D < 2; ++D) {
		bool const dim(bool(D) ^ long_dim); // place along walls in long dim first
		float const room_len(place_area.get_sz_dim(dim));
		unsigned const num(room_len/server_period); // take the floor
		if (num == 0) continue; // not enough space for a server in this dim
		float const server_spacing(room_len/num);
		center[dim] = place_area.d[dim][0] + 0.5*server_spacing; // first position at half spacing

		for (unsigned n = 0; n < num; ++n, center[dim] += server_spacing) {
			set_wall_width(server, center[dim], server_hwidth, dim); // position along the wall

			for (unsigned dir = 0; dir < 2; ++dir) {
				float const dir_sign(dir ? -1.0 : 1.0);
				center[!dim] = place_area.d[!dim][dir] + dir_sign*server_hdepth;
				set_wall_width(server, center[!dim], server_hdepth, !dim); // position from the wall
				cube_t server_exp(server);
				server_exp.expand_in_dim(dim, server_hwidth); // check for more side/width spacing for doors
				
				// Note: overlaps_other_room_obj includes previously placed servers, so we don't have to check for intersections at the corners of rooms
				if (is_obj_placement_blocked(server_exp, room, 1) || overlaps_other_room_obj(server, objs_start)) { // no space for server; try computer instead
					set_wall_width(computer,  center[ dim], comp_hwidth, dim); // position along the wall
					set_wall_width(computer, (place_area.d[!dim][dir] + 1.2*dir_sign*comp_hdepth), comp_hdepth, !dim); // position from the wall
					if (is_obj_placement_blocked(computer, room, 1) || overlaps_other_room_obj(computer, objs_start)) continue;
					objs.emplace_back(computer, TYPE_COMPUTER, room_id, !dim, !dir, 0, tot_light_amt);
					++num_comps;
					continue;
				}
				objs.emplace_back(server, TYPE_SERVER, room_id, !dim, !dir, 0, tot_light_amt);
				cube_t blocker(server);
				blocker.d[!dim][ dir]  = server.d[!dim][!dir]; // front of server
				blocker.d[!dim][!dir] += dir_sign*server_width; // add space in the front for the door to open (don't block with another server)
				objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, 0, RO_FLAG_INVIS);
				++num_servers;
			} // for dir
		} // for n
	} // for dim
	if (num_servers == 0 && num_comps == 0) return 0; // both servers and computers count
	
	if (num_servers > 0) { // add a keyboard to the master server
		unsigned const master_server(rgen.rand() % num_servers);

		for (unsigned i = objs_start, server_ix = 0; i < objs.size(); ++i) {
			room_object_t const &server(objs[i]);
			if (server.type != TYPE_SERVER  ) continue;
			if (server_ix++ != master_server) continue;
			float const kbd_hwidth(0.8*server_hwidth), kbd_depth(0.6*kbd_hwidth), kbd_height(0.04*kbd_hwidth); // slightly flatter than regular keyboards
			bool const dim(server.dim), dir(server.dir);
			float const kbd_z1(server.z1() + 0.57*server.dz()), server_front(server.d[dim][dir]);
			cube_t keyboard;
			set_cube_zvals(keyboard, kbd_z1, kbd_z1+kbd_height);
			keyboard.d[dim][!dir] = server_front; // at front of server
			keyboard.d[dim][ dir] = server_front + (dir ? 1.0 : -1.0)*kbd_depth; // sticks out of the front
			set_wall_width(keyboard, server.get_center_dim(!dim), kbd_hwidth, !dim);
			if (is_obj_placement_blocked(keyboard, room, 1)) break; // Note: not checking overlaps_other_room_obj() because it will overlap server blockers
			objs.emplace_back(keyboard, TYPE_KEYBOARD, room_id, dim, dir, RO_FLAG_HANGING, tot_light_amt); // add as white, will be drawn with gray/black texture
			break;
		} // for i
	}
	// maybe add laptops on top of some servers, to reward the player for finding this room
	for (unsigned i = objs_start; i < objs.size(); ++i) {
		room_object_t const &server(objs[i]);
		if (server.type != TYPE_SERVER) continue;
		if (rgen.rand_float() > 0.2)    continue; // place laptops 20% of the time
		bool const dim(server.dim), dir(server.dir);
		float const server_front(server.d[dim][dir]); // copy before reference is invalidated
		if (!place_laptop_on_obj(rgen, server, room_id, tot_light_amt)) continue; // no avoid, use_dim_dir=0
		// make the laptop hang over the edge of the front of the server so that the player can see and take it
		room_object_t &laptop(objs.back());
		float const xlate(server_front - laptop.d[dim][dir] + (dir ? 1.0 : -1.0)*rgen.rand_uniform(0.05, 0.35)*laptop.get_sz_dim(dim));
		laptop.translate_dim(dim, xlate);
		laptop.flags |= RO_FLAG_HANGING;
	} // for i
	add_door_sign("Server Room", room, zval, room_id);
	return 1;
}

bool building_t::add_security_room_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt, unsigned objs_start) { // for office buildings
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_TV)) return 0; // no TV/monitor model, can't create a security room
	float const floor_spacing(get_window_vspace()), ceil_zval(zval + get_floor_ceil_gap());
	float const start_zval(zval + 0.3*floor_spacing), place_z_range(ceil_zval - start_zval); // should be above the top of the desk
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
	float const tv_height(0.23*floor_spacing*rgen.rand_uniform(1.0, 1.2)); // common height for all monitors
	float const tv_hwidth(0.5*tv_height*sz.y/sz.z), tv_depth(tv_height*sz.x/sz.z);
	float const vert_space(1.25*tv_height), horiz_space(2.5*tv_hwidth);
	unsigned const num_rows(round_fp(place_z_range/vert_space));
	float const row_spacing(place_z_range/num_rows);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t breaker_panel, tv;
	vect_cube_t blockers;
	interior->security_room_ix = room_id;
	
	if (1) { // add breaker panel
		vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // get interior doors
		assert(!doorways.empty());
		door_stack_t const &door(doorways.front()); // choose the first door (there is likely only one)
		bool const side(!door.get_check_dirs()), dim(door.dim), dir(door.get_center_dim(dim) > room.get_center_dim(dim)); // the wall the door is on
		float const door_edge(door.d[!dim][side]), wall_edge(room_bounds.d[!dim][side]);
		float const wall_len(fabs(door_edge - wall_edge)), wall_center(0.5*(door_edge + wall_edge)), wall_pos(room_bounds.d[dim][dir]);
		float const width(min(0.5f*wall_len, rgen.rand_uniform(0.25, 0.35)*floor_spacing)), depth(0.04*floor_spacing);
		set_cube_zvals(breaker_panel, (ceil_zval - 0.75*floor_spacing), (ceil_zval - rgen.rand_uniform(0.25, 0.3)*floor_spacing));
		set_wall_width(breaker_panel, wall_center, 0.5*width, !dim);
		breaker_panel.d[dim][ dir] = wall_pos;
		breaker_panel.d[dim][!dir] = wall_pos + (dir ? -1.0 : 1.0)*depth;
		add_breaker_panel(rgen, breaker_panel, ceil_zval, door.dim, dir, room_id, tot_light_amt);
		breaker_panel.d[dim][!dir] += (dir ? -1.0 : 1.0)*width; // add padding for desk placement
		blockers.push_back(breaker_panel);
	}
	add_desk_to_room(rgen, room, blockers, DK_GRAY, zval, room_id, tot_light_amt, objs_start, 0, 0, 1); // is_basement=0, desk_ix=0, no_computer=1

	// add computer monitors along all walls that don't have doors
	for (unsigned dim = 0; dim < 2; ++dim) {
		float const wall_len(room_bounds.get_sz_dim(!dim));
		unsigned const num_cols(floor(wall_len/horiz_space));
		float const col_spacing(wall_len/num_cols), col_start(room_bounds.d[!dim][0] + 0.5*col_spacing);

		for (unsigned dir = 0; dir < 2; ++dir) {
			tv.d[dim][ dir] = room_bounds.d[dim][dir]; // on the wall
			tv.d[dim][!dir] = tv.d[dim][dir] + (dir ? -1.0 : 1.0)*tv_depth;

			for (unsigned col = 0; col < num_cols; ++col) { // XY
				set_wall_width(tv, (col_start + ((dim ^ dir) ? (num_cols-col-1) : col)*col_spacing), tv_hwidth, !dim); // ordered left to right
				if (!breaker_panel.is_all_zeros() && breaker_panel.intersects_xy(tv)) continue; // ignore zvals so that we don't put a monitor above the breaker panel

				for (unsigned row = 0; row < num_rows; ++row) { // Z
					float const z1(start_zval + row*row_spacing);
					set_cube_zvals(tv, z1, z1+tv_height);
					if (is_obj_placement_blocked(tv, room, 1)) continue;
					//if (overlaps_other_room_obj(tv, objs_start)) continue; // not needed since there are no objects placed first?
					objs.emplace_back(tv, TYPE_MONITOR, room_id, dim, !dir, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), tot_light_amt, SHAPE_SHORT, BLACK); // monitors are shorter than TVs
					offset_hanging_tv(objs.back());
					set_obj_id(objs);
					objs.back().obj_id &= ~1; // on by default; strip off LSB
				} // for row
			} // for col
		} // for dir
	} // for dim
	add_door_sign("Security Room", room, zval, room_id);
	return 1;
}

bool building_t::place_book_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id,
	float tot_light_amt, unsigned objs_start, bool use_dim_dir, unsigned extra_flags, bool skip_if_overlaps)
{
	point center(place_on.get_cube_center());
	for (unsigned d = 0; d < 2; ++d) {center[d] += 0.1*place_on.get_sz_dim(d)*rgen.rand_uniform(-1.0, 1.0);} // add a slight random shift
	cube_t const room_bounds(get_walkable_room_bounds(get_room(room_id)));
	assert(room_bounds.contains_pt(center));
	// book is randomly oriented for tables and rotated 90 degrees from desk orient
	bool const dim(use_dim_dir ? !place_on.dim : rgen.rand_bool()), dir(use_dim_dir ? (place_on.dir^place_on.dim) : rgen.rand_bool());
	cube_t book(get_book_bcube(rgen, point(center.x, center.y, place_on.z2()), get_window_vspace(), dim, dir));
	book.intersect_with_cube(room_bounds); // clip to room bounds
	vect_room_object_t &objs(interior->room_geom->objs);

	// check if there's anything in the way; only handling pens and pencils here; paper is ignored, and larger objects should already be handled
	for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
		if (i->type != TYPE_PEN && i->type != TYPE_PENCIL) {
			if (skip_if_overlaps && *i != place_on) return 0; // something else like a cup
			continue;
		}
		if (!i->intersects(book)) continue;
		set_cube_zvals(book, i->z2(), i->z2()+book.dz()); // place book on top of object; maybe the book should be tilted?
	}
	colorRGBA const color(book_colors[rgen.rand() % NUM_BOOK_COLORS]);
	unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT | extra_flags);
	if (place_on.is_glass_table()) {flags |= RO_FLAG_HAS_EXTRA;} // flag so that shadows are enabled
	objs.emplace_back(book, TYPE_BOOK, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, color); // Note: invalidates place_on reference
	set_obj_id(objs);
	return 1;
}

bool place_bottle_on_obj(rand_gen_t &rgen, cube_t const &place_on, vect_room_object_t &objs, float vspace,
	unsigned rid, float lamt, unsigned max_type, vect_cube_t const &avoid, bool at_z1, bool allow_transparent=0)
{
	float const height(vspace*rgen.rand_uniform(0.075, 0.12)), radius(vspace*rgen.rand_uniform(0.012, 0.018));
	if (min(place_on.dx(), place_on.dy()) < 5.0*radius) return 0; // surface is too small to place this bottle
	if (at_z1 && height > place_on.dz())                return 0; // too tall to place in object
	cube_t const bottle(place_cylin_object(rgen, place_on, radius, height, 2.0*radius, at_z1));
	if (has_bcube_int(bottle, avoid)) return 0; // only make one attempt
	objs.emplace_back(bottle, TYPE_BOTTLE, rid, 0, 0, RO_FLAG_NOCOLL, lamt, SHAPE_CYLIN);
	objs.back().set_as_bottle(rgen.rand(), max_type, 0, 0, 0, allow_transparent); // can be empty
	return 1;
}
bool place_dcan_on_obj(rand_gen_t &rgen, cube_t const &place_on, vect_room_object_t &objs, float vspace,
	unsigned rid, float lamt, unsigned max_type, vect_cube_t const &avoid, bool at_z1)
{
	float const height(0.06*vspace), radius(0.016*vspace);
	if (min(place_on.dx(), place_on.dy()) < 5.0*radius) return 0; // surface is too small to place this can
	if (at_z1 && height > place_on.dz())                return 0; // too tall to place in object
	cube_t const can(place_cylin_object(rgen, place_on, radius, height, 2.0*radius, at_z1));
	if (has_bcube_int(can, avoid)) return 0; // only make one attempt
	objs.emplace_back(can, TYPE_DRINK_CAN, rid, 0, 0, RO_FLAG_NOCOLL, lamt, SHAPE_CYLIN, LT_GRAY);
	objs.back().obj_id = (uint16_t)(rgen.rand() & 127); // strip off empty bit
	set_max_drink_can_type(objs.back(), max_type);
	return 1;
}
bool building_t::place_bottle_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt,
	vect_cube_t const &avoid, bool place_at_z1, unsigned max_type)
{
	return ::place_bottle_on_obj(rgen, place_on, interior->room_geom->objs, get_window_vspace(), room_id, tot_light_amt, BOTTLE_TYPE_WINE, avoid, place_at_z1);
}
bool building_t::place_dcan_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt,
	vect_cube_t const &avoid, bool place_at_z1, unsigned max_type) // drink can
{
	return ::place_dcan_on_obj(rgen, place_on, interior->room_geom->objs, get_window_vspace(), room_id, tot_light_amt, DRINK_CAN_TYPE_BEER, avoid, place_at_z1);
}

colorRGBA choose_pot_color(rand_gen_t &rgen) {
	unsigned const num_colors = 8;
	colorRGBA const pot_colors[num_colors] = {LT_GRAY, GRAY, DK_GRAY, BKGRAY, WHITE, LT_BROWN, RED, colorRGBA(1.0, 0.35, 0.18)};
	return pot_colors[rgen.rand() % num_colors];
}
bool building_t::place_plant_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, float sz_scale, vect_cube_t const &avoid) {
	float const window_vspacing(get_window_vspace()), height(rgen.rand_uniform(0.25, 0.4)*sz_scale*window_vspacing), max_radius(min(place_on.dx(), place_on.dy())/3.0f);
	vect_room_object_t &objs(interior->room_geom->objs);

	if (building_obj_model_loader.is_model_valid(OBJ_MODEL_PLANT)) { // prefer to place potted plant models, if any exist
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_PLANT)); // D, W, H
		float const radius_to_height(0.25*(sz.x + sz.y)/sz.z), radius(min(radius_to_height*height*0.8f, max_radius)); // cylindrical, 80% the size of regular plants
		cube_t const plant(place_cylin_object(rgen, place_on, radius, radius/radius_to_height, 1.2*radius)); // recompute height from radius
		
		if (!has_bcube_int(plant, avoid)) { // only make one attempt
			unsigned const item_flags(rgen.rand()); // choose a random potted plant model if there are more than one
			objs.emplace_back(plant, TYPE_PLANT_MODEL, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_BOT), tot_light_amt, SHAPE_CYLIN, WHITE, item_flags);
			return 1;
		} // else try to place a non-model plant
	}
	float const radius(min(rgen.rand_uniform(0.06, 0.08)*sz_scale*window_vspacing, max_radius));
	cube_t const plant(place_cylin_object(rgen, place_on, radius, height, 1.2*radius));
	if (has_bcube_int(plant, avoid)) return 0; // only make one attempt
	objs.emplace_back(plant, TYPE_PLANT, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_BOT), tot_light_amt, SHAPE_CYLIN, choose_pot_color(rgen)); // flag for full player pickup
	set_obj_id(objs);
	return 1;
}

bool building_t::place_laptop_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid, bool use_dim_dir) {
	point center(place_on.get_cube_center());
	for (unsigned d = 0; d < 2; ++d) {center[d] += 0.1*place_on.get_sz_dim(d)*rgen.rand_uniform(-1.0, 1.0);} // add a slight random shift
	bool const dim(use_dim_dir ? place_on.dim : rgen.rand_bool()), dir(use_dim_dir ? (place_on.dir^place_on.dim^1) : rgen.rand_bool()); // Note: dir is inverted
	float const width(0.136*get_window_vspace());
	vector3d sz;
	sz[!dim] = width;
	sz[ dim] = 0.7*width;  // depth
	sz.z     = 0.06*width; // height
	point const llc(center.x, center.y, place_on.z2());
	cube_t laptop(llc, (llc + sz));
	if (has_bcube_int(laptop, avoid)) return 0; // only make one attempt
	if (!get_walkable_room_bounds(get_room(room_id)).contains_cube(laptop)) return 0; // may clip through a wall
	interior->room_geom->objs.emplace_back(laptop, TYPE_LAPTOP, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt); // Note: invalidates place_on reference
	return 1;
}

bool building_t::place_pizza_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid) {
	float const width(0.15*get_window_vspace());
	if (min(place_on.dx(), place_on.dy()) < 1.2*width) return 0; // place_on is too small
	cube_t pizza;
	gen_xy_pos_for_cube_obj(pizza, place_on, vector3d(0.5*width, 0.5*width, 0.0), 0.1*width, rgen);
	bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
	if (has_bcube_int(pizza, avoid)) return 0; // only make one attempt
	interior->room_geom->objs.emplace_back(pizza, TYPE_PIZZA_BOX, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt); // Note: invalidates place_on reference
	set_obj_id(interior->room_geom->objs);
	return 1;
}

float get_plate_radius(rand_gen_t &rgen, cube_t const &place_on, float window_vspacing) {
	return min(rgen.rand_uniform(0.05, 0.07)*window_vspacing, 0.25f*min(place_on.dx(), place_on.dy()));
}

bool building_t::place_plate_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid) {
	float const radius(get_plate_radius(rgen, place_on, get_window_vspace()));
	cube_t const plate(place_cylin_object(rgen, place_on, radius, 0.1*radius, 1.1*radius));
	if (has_bcube_int(plate, avoid)) return 0; // only make one attempt
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(plate, TYPE_PLATE, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN);
	set_obj_id(objs);
	return 1;
}

bool building_t::place_cup_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid, bool make_empty) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_CUP)) return 0;
	float const height(0.06*get_window_vspace()), radius(height*get_radius_for_square_model(OBJ_MODEL_CUP)); // almost square
	if (min(place_on.dx(), place_on.dy()) < 2.5*radius) return 0; // surface is too small to place this cup
	cube_t const cup(place_cylin_object(rgen, place_on, radius, height, 1.2*radius));
	if (has_bcube_int(cup, avoid)) return 0; // only make one attempt
	unsigned flags(RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT);
	if (!make_empty && rgen.rand_bool()) {flags |= RO_FLAG_NONEMPTY;} // add coffee to the cup
	// random dim/dir, plus more randomness on top
	interior->room_geom->objs.emplace_back(cup, TYPE_CUP, room_id, rgen.rand_bool(), rgen.rand_bool(), flags, tot_light_amt, SHAPE_CYLIN);
	return 1;
}

bool building_t::place_toy_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid) {
	float const height(0.11*get_window_vspace()), radius(0.5f*height*0.67f);
	if (min(place_on.dx(), place_on.dy()) < 2.5*radius) return 0; // surface is too small to place this toy
	cube_t const toy(place_cylin_object(rgen, place_on, radius, height, 1.1*radius));
	if (has_bcube_int(toy, avoid)) return 0; // only make one attempt
	interior->room_geom->objs.emplace_back(toy, TYPE_TOY, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN);
	set_obj_id(interior->room_geom->objs); // used for color selection
	return 1;
}

bool building_t::place_banana_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_BANANA)) return 0;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_BANANA));
	float const length(0.083*get_window_vspace()), width(length*sz.y/sz.x), height(length*sz.z/sz.x);
	if (min(place_on.dx(), place_on.dy()) < 1.1*length) return 0; // surface is too small to place this banana
	point bsz(0.5*length, 0.5*width, height);
	bool const dim(rgen.rand_bool());
	if (dim) {swap(bsz.x, bsz.y);}
	cube_t bbc;
	gen_xy_pos_for_cube_obj(bbc, place_on, bsz, height, rgen);
	if (has_bcube_int(bbc, avoid)) return 0; // only make one attempt
	interior->room_geom->objs.emplace_back(bbc, TYPE_BANANA, room_id, dim, rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt);
	return 1;
}
bool building_t::place_apple_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, vect_cube_t const &avoid) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_APPLE)) return 0;
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_APPLE));
	float const radius(rgen.rand_uniform(0.02, 0.025)*get_window_vspace()), height(4.0*radius*sz.z/(sz.x + sz.y)); // assumes xsize == ysize
	if (min(place_on.dx(), place_on.dy()) < 2.5*radius) return 0; // surface is too small to place this apple
	cube_t bbc;
	gen_xy_pos_for_round_obj(bbc, place_on, radius, height, 1.2*radius, rgen);
	if (has_bcube_int(bbc, avoid)) return 0; // only make one attempt
	unsigned const item_flags(rgen.rand()); // random apple sub-model
	interior->room_geom->objs.emplace_back(bbc, TYPE_APPLE, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_SPHERE, WHITE, item_flags);
	return 1;
}

bool building_t::place_phone_on_obj(rand_gen_t &rgen, cube_t const &place_on, unsigned room_id, float tot_light_amt, bool dim, bool dir, float overhang_amt) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_PHONE)) return 0;
	float const height(0.03*get_window_vspace()), radius((1.0 - overhang_amt)*height*get_radius_for_square_model(OBJ_MODEL_PHONE)); // almost square
	if (min(place_on.dx(), place_on.dy()) < 2.1*radius) return 0; // surface is too small to place this phone
	cube_t const phone(place_cylin_object(rgen, place_on, radius, height, radius));
	interior->room_geom->objs.emplace_back(phone, TYPE_CONF_PHONE, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt);
	return 1;
}

bool building_t::add_rug_to_room(rand_gen_t rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!room_object_t::enable_rugs()) return 0; // disabled
	vector2d const room_sz(room.get_size_xy());
	bool const min_dim(room_sz.y < room_sz.x);
	float const ar(rgen.rand_uniform(0.65, 0.85)), length(min(0.7f*room_sz[min_dim]/ar, room_sz[!min_dim]*rgen.rand_uniform(0.4, 0.7))), width(length*ar);
	cube_t rug;
	set_cube_zvals(rug, zval, (zval + get_rug_thickness())); // almost flat
	vect_room_object_t &objs(interior->room_geom->objs);
	float sz_scale(1.0);

	for (unsigned n = 0; n < 10; ++n) { // make 10 attempts at choosing a valid alignment
		vector3d center(room.get_cube_center()); // Note: zvals ignored
		bool valid_placement(1);

		for (unsigned d = 0; d < 2; ++d) {
			float const radius(0.5*((bool(d) == min_dim) ? width : length)), scaled_radius(radius*sz_scale);
			center[d] += (0.05f*room_sz[d] + (radius - scaled_radius))*rgen.rand_uniform(-1.0, 1.0); // slight random misalignment, increases with decreasing sz_scale
			rug.d[d][0] = center[d] - radius;
			rug.d[d][1] = center[d] + radius;
		}
		for (auto i = objs.begin() + objs_start; i != objs.end() && valid_placement; ++i) { // check for objects overlapping the rug
			if (i->type == TYPE_FLOORING) continue; // allow placing rugs over flooring
			if (!i->intersects(rug))      continue;

			if (bldg_obj_types[i->type].attached) { // rugs can't overlap these object types; first, see if we can shrink the rug on one side and get it to fit
				float max_area(0.0);
				cube_t best_cand;

				for (unsigned dim = 0; dim < 2; ++dim) {
					for (unsigned dir = 0; dir < 2; ++dir) {
						cube_t cand(rug);
						cand.d[dim][dir] = i->d[dim][!dir] + (dir ? -1.0 : 1.0)*0.025*rug.get_sz_dim(dim); // leave a small gap
						float const area(cand.dx()*cand.dy());
						if (area > max_area) {best_cand = cand; max_area = area;}
					}
				}
				if (max_area > 0.8*rug.dx()*rug.dy()) {rug = best_cand;} // good enough
				else {valid_placement = 0;} // shrink is not enough, try again
			}
			else if (i->type == TYPE_TABLE || i->type == TYPE_DESK || i->type == TYPE_FCABINET || i->type == TYPE_CONF_TABLE) { // rugs can't partially overlap these object types
				valid_placement = rug.contains_cube_xy(*i); // don't expand as that could cause the rug to intersect a previous object
				// maybe beds should be included as well, but then rugs are unlikely to be placed in bedrooms
			}
		} // for i
		if (valid_placement && interior->is_blocked_by_stairs_or_elevator(rug)) {valid_placement = 0;} // check stairs (required for ext basement rooms); no need to check doors

		if (valid_placement) {
			cube_t place_area(room);
			place_area.expand_by_xy(-0.1*get_wall_thickness()); // add small border to avoid alpha blending artifacts if the rug intersects the wall
			rug.intersect_with_cube_xy(place_area); // make sure the rug stays within the room bounds

			if (rug.is_strictly_normalized()) {
				objs.emplace_back(rug, TYPE_RUG, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt);
				room_object_t &rug_obj(objs.back());
				rug_obj.obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix); // determines rug texture
				
				// don't use the same texture as a blanket because that looks odd
				for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
					if (i->type == TYPE_BLANKET && i->get_rug_tid() == rug_obj.get_rug_tid()) {++rug_obj.obj_id;} // select a different texture
				}
				return 1;
			}
		}
		sz_scale *= 0.9; // decrease rug size and try again
	} // for n
	return 0;
}

// return value: 0=invalid, 1=valid and good, 2=valid but could be better
int building_t::check_valid_picture_placement(room_t const &room, cube_t const &c, float width, float zval, bool dim, bool dir, unsigned objs_start) const {
	assert(interior != nullptr);
	float const wall_thickness(get_wall_thickness()), clearance(4.0*wall_thickness), side_clearance(1.0*wall_thickness);
	cube_t tc(c), keepout(c);
	tc.expand_in_dim(!dim, 0.1*width); // expand slightly to account for frame
	//keepout.z1() = zval; // extend to the floor
	keepout.z1() -= 0.10*c.dz(); // more padding on the bottom
	keepout.z2() += 0.05*c.dz(); // small padding on the top for the frame
	keepout.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance;
	keepout.d[dim][ dir] -= (dir ? -1.0 : 1.0)*0.5*wall_thickness; // move into wall to capture thin objects such as wall vents
	keepout.expand_in_dim(!dim, side_clearance); // make sure there's space for the frame
	if (overlaps_other_room_obj(keepout, objs_start, 1))   return 0; // check_all=1, to include outlets, vents, etc.
	bool const inc_open(!is_house && !room.is_office);
	if (is_cube_close_to_doorway(tc, room, 0.0, inc_open)) return 0; // bad placement
	// Note: it's not legal to guard the below check with (room.has_stairs || room.has_elevator) because room.has_stairs may not be set for stack connector stairs that split a wall
	float const se_clearance(2.0*get_min_front_clearance_inc_people() + 4.0*wall_thickness); // enough space to walk in front, and then some
	if (interior->is_blocked_by_stairs_or_elevator(tc, se_clearance)) return 0; // check stairs and elevators
	if (!check_cube_within_part_sides(tc)) return 0; // handle non-cube buildings
	if (overlaps_or_adj_int_window   (tc)) return 0; // check interior windows
	if (!inc_open && !room.is_hallway && is_cube_close_to_doorway(tc, room, 0.0, 1)) return 2; // success, but could be better (doors never open into hallway)

	if (has_complex_floorplan && c.z1() > ground_floor_z1) { // check for office building whiteboards placed on room sides that aren't true walls; skip basements
		cube_t test_cube(c);
		test_cube.expand_by_xy(2.0*wall_thickness); // max sure it extends through the wall
		unsigned num_parts_int(0);

		for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
			if (p->intersects(test_cube)) {++num_parts_int;}
		}
		assert(num_parts_int > 0);

		if (num_parts_int > 1) { // on the border between two parts, check if there's a wall between them
			cube_t wall_mount(c);
			wall_mount.d[dim][0] = wall_mount.d[dim][1] = c.d[dim][dir] + (dir ? 1.0 : -1.0)*0.5*wall_thickness; // should be in the center of the wall
			bool found_wall(0);

			for (auto const &w: interior->walls[dim]) {
				if (w.contains_cube(wall_mount)) {found_wall = 1; break;}
			}
			if (!found_wall) return 0;
		}
	}
	return 1; // success
}

bool building_t::hang_pictures_whiteboard_chalkboard_in_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id,
	float tot_light_amt, unsigned objs_start, unsigned floor_ix, bool is_basement, unsigned pref_orient)
{
	if (!room_object_t::enable_pictures()) return 0; // disabled
	
	if (!is_house && !room.is_office) {
		if (room.is_hallway) return 0; // no pictures or whiteboards in office building hallways (what about rooms with stairs?)
		// room in a commercial building or hotel/apartment - add whiteboard/picture when there is a full wall to use
	}
	if (room.is_sec_bldg) return 0; // no pictures in secondary buildings
	if (room.get_room_type(0) == RTYPE_STORAGE) return 0; // no pictures or whiteboards in storage rooms (always first floor)
	cube_t const &part(get_part_for_room(room));
	float const floor_height(get_window_vspace()), wall_thickness(get_wall_thickness());
	bool const no_ext_walls(!is_basement && (has_int_windows() || !is_cube())); // don't place on ext walls with windows or non-square orients
	vect_room_object_t &objs(interior->room_geom->objs);
	bool was_hung(0);

	if (!is_residential() || room.is_office) { // add whiteboards
		room_type const rtype(room.get_room_type(floor_ix));
		bool const is_conference(rtype == RTYPE_CONF); // conference rooms always have a whiteboard
		if (!is_conference && pref_orient == 4 && rgen.rand_float() < 0.1) return 0; // skip 10% of the time; don't skip if pref_orient was set (classroom)
		bool const pref_dim((pref_orient < 4) ? (pref_orient >> 1) : rgen.rand_bool()), pref_dir((pref_orient < 4) ? (pref_orient & 1) : rgen.rand_bool());
		bool const use_blackboards(is_school());
		float const floor_thick(get_floor_thickness()); // blackboard for schools, whiteboard for offices
		colorRGBA const color(use_blackboards ? GRAY_BLACK : WHITE);

		for (unsigned dim2 = 0; dim2 < 2; ++dim2) {
			for (unsigned dir2 = 0; dir2 < 2; ++dir2) {
				bool const dim(bool(dim2) ^ pref_dim), dir(bool(dir2) ^ pref_dir);
				if (no_ext_walls && fabs(room.d[dim][dir] - part.d[dim][dir]) < 1.1*wall_thickness) continue; // on part boundary, likely ext wall where there may be windows, skip
				cube_t c(room);
				set_cube_zvals(c, zval+0.25*floor_height, zval+0.9*floor_height-floor_thick);
				c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*0.6*wall_thickness;
				// translate by half wall thickness if not interior hallway or office wall
				if (!(room.inc_half_walls() && classify_room_wall(room, zval, dim, dir, 0) != ROOM_WALL_EXT)) {c.translate_dim(dim, (dir ? 1.0 : -1.0)*0.5*wall_thickness);}
				float const room_len(room.get_sz_dim(!dim));
				if (room_len < 1.0*floor_height) continue; // wall too short for a whiteboard
				c.expand_in_dim(!dim, -0.2*room_len); // xy_space
				float const wb_len(c.get_sz_dim(!dim)), wb_max_len(3.0*floor_height);
				if (wb_len > wb_max_len) {c.expand_in_dim(!dim, -0.5*(wb_len - wb_max_len));} // shrink to max length if needed
				
				if (!check_valid_picture_placement(room, c, 0.6*room_len, zval, dim, dir, objs_start)) { // fails wide/tall placement
					cube_t c_prev(c);
					c.expand_in_dim(!dim, -0.167*c.get_sz_dim(!dim)); // shrink width a bit and try again
					
					if (!check_valid_picture_placement(room, c, 0.4*room_len, zval, dim, dir, objs_start)) { // fails narrow/tall placement
						c = c_prev;
						c.z2() -= 0.15*c.dz(); // shrink height and try again with wide placement
						if (!check_valid_picture_placement(room, c, 0.6*room_len, zval, dim, dir, objs_start)) continue; // give up/fail
					}
				}
				assert(c.is_strictly_normalized());
				room_object_t wboard(c, TYPE_WBOARD, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, color); // whiteboard faces dir opposite the wall
				if (use_blackboards && rgen.rand_float() < 0.25) {wboard.item_flags = 1;} // flag as having math writing
				objs.push_back(wboard);
				cube_t blocker(wboard);
				blocker.d[dim][!dir] += (dir ? -1.0 : 1.0)*get_min_front_clearance_inc_people();
				objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, !dir, RO_FLAG_INVIS);

				if (rgen.rand_float() < 0.8) { // add marker/chalk and maybe eraser
					float const marker_hlen(0.5*0.06*floor_height), marker_radius(0.005*floor_height);
					cube_t const ledge(get_whiteboard_marker_ledge(wboard));

					if (marker_hlen < 0.1*wb_len && marker_radius < 0.4*ledge.get_sz_dim(dim)) { // have space for a marker
						float const end_space(1.1*marker_hlen), side_space(1.1*marker_radius);
						cube_t marker;
						set_cube_zvals(marker, ledge.z2(), ledge.z2()+2.0*marker_radius);
						set_wall_width(marker, rgen.rand_uniform(ledge.d[!dim][0]+end_space,  ledge.d[!dim][1]-end_space ), marker_hlen,  !dim);
						set_wall_width(marker, rgen.rand_uniform(ledge.d[ dim][0]+side_space, ledge.d[ dim][1]-side_space), marker_radius, dim);
						colorRGBA const marker_color(use_blackboards ? WHITE : marker_colors[rgen.rand() & 7]); // white chalk for school blackboards, colored markers for whiteboards
						objs.emplace_back(marker, TYPE_MARKER, room_id, !dim, rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, marker_color);
						float const eraser_hlen(0.5*0.055*floor_height), eraser_hwidth(0.5*0.025*floor_height), eraser_height(0.012*floor_height);

						if (eraser_hlen < 0.1*wb_len && eraser_hwidth < 0.48*ledge.get_sz_dim(dim) && rgen.rand_float() < 0.7) { // add eraser
							cube_t eraser;
							set_cube_zvals(eraser, ledge.z2(), ledge.z2()+eraser_height);
							set_wall_width(eraser, rgen.rand_uniform(ledge.d[!dim][0]+eraser_hlen,   ledge.d[!dim][1]-eraser_hlen  ), eraser_hlen,  !dim);
							set_wall_width(eraser, rgen.rand_uniform(ledge.d[ dim][0]+eraser_hwidth, ledge.d[ dim][1]-eraser_hwidth), eraser_hwidth, dim);

							if (!eraser.intersects(marker)) { // skip if eraser intersects marker
								objs.emplace_back(eraser, TYPE_ERASER, room_id, !dim, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, BKGRAY);
							}
						}
					}
				}
				if (use_blackboards && rtype == RTYPE_CLASS) {add_objects_next_to_classroom_chalkboard(rgen, wboard, room, zval, objs_start);} // add classroom flag, clock, etc.
				if (!is_conference || rgen.rand_bool()) return 1; // only need to add one, except for conference rooms
			} // for dir
		} // for dim
		return 0;
	}
	// add pictures
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			float const wall_pos(room.d[dim][dir]);
			if (no_ext_walls && fabs(room.d[dim][dir] - part.d[dim][dir]) < 1.1*wall_thickness) continue; // on part boundary, likely ext wall where there may be windows, skip
			if (!room.is_hallway && rgen.rand_float() < 0.2) continue; // skip 20% of the time unless it's a hallway
			float const height(floor_height*rgen.rand_uniform(0.3, 0.6)*(is_basement ? 0.8 : 1.0)); // smaller pictures in basement to avoid the pipes
			float const width(height*rgen.rand_uniform(1.5, 2.0)); // width > height
			if (width > 0.8*room.get_sz_dim(!dim)) continue; // not enough space
			float const base_shift((dir ? -1.0 : 1.0)*0.5*wall_thickness); // half a wall's thickness in dir
			point center;
			center[ dim] = wall_pos;
			center[!dim] = room.get_center_dim(!dim);
			center.z     = zval + rgen.rand_uniform(0.45, 0.55)*floor_height; // move up
			float const lo(room.d[!dim][0] + 0.7*width), hi(room.d[!dim][1] - 0.7*width);
			cube_t best_pos;

			for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to choose a position along the wall; first iteration is the center
				if (n > 0) { // try centered first, then non-centered
					if (hi - lo < width) break; // not enough space to shift, can't place this picture
					center[!dim] = rgen.rand_uniform(lo, hi);
				}
				cube_t c(center);
				c.expand_in_z(0.5*height);
				c.d[dim][!dir] += 0.2*base_shift; // move out to prevent z-fighting
				// add an additional half wall thickness for interior hallway and office walls
				if (room.inc_half_walls() && classify_room_wall(room, zval, dim, dir, 0) != ROOM_WALL_EXT) {c.translate_dim(dim, base_shift);}
				c.expand_in_dim(!dim, 0.5*width);
				int const ret(check_valid_picture_placement(room, c, width, zval, dim, dir, objs_start));
				if (ret == 0) continue; // invalid, retry
				best_pos = c;
				if (ret == 1) break; // valid and good - keep this pos
			} // for n
			if (best_pos.is_all_zeros()) continue; // failed placement
			assert(best_pos.is_strictly_normalized());
			unsigned flags(RO_FLAG_NOCOLL);
			if (rgen.rand_bool()) {flags |= RO_FLAG_HAS_EXTRA;} // flag as being abstract (shader) art
			objs.emplace_back(best_pos, TYPE_PICTURE, room_id, dim, !dir, flags, tot_light_amt); // picture faces dir opposite the wall
			objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 17*floor_ix + 31*mat_ix + 61*dim + 123*dir); // determines picture texture
			was_hung = 1;
		} // for dir
	} // for dim
	return was_hung;
}

void building_t::add_plants_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num) {
	if (num == 0) return;
	float const window_vspacing(get_window_vspace());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-get_trim_thickness()); // shrink to leave a small gap
	zval += 1.1*get_flooring_thick(); // move up slightly to avoid z-fighting of bottom when the dirt is taken
	
	for (unsigned n = 0; n < num; ++n) {
		float const height(rgen.rand_uniform(0.6, 0.9)*window_vspacing), width(rgen.rand_uniform(0.15, 0.35)*window_vspacing);
		vector3d const sz_scale(width/height, width/height, 1.0);
		place_obj_along_wall(TYPE_PLANT, room, height, sz_scale, rgen, zval, room_id, tot_light_amt,
			place_area, objs_start, 0.0, 0, 4, 0, choose_pot_color(rgen), 0, SHAPE_CYLIN); // no clearance, pref_orient, or color
	}
}

void building_t::add_boxes_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned max_num) {
	if (max_num == 0) return; // why did we call this?
	float const window_vspacing(get_window_vspace());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // shrink to leave a small gap for open flaps
	unsigned const num(rgen.rand() % (max_num+1));
	bool const allow_crates(!is_residential() && room.is_ext_basement()); // backrooms
	vect_room_object_t &objs(interior->room_geom->objs);

	for (unsigned n = 0; n < num; ++n) {
		vector3d sz;
		gen_crate_sz(sz, rgen, window_vspacing);
		sz *= 1.5; // make larger than storage room boxes
		room_object const type((allow_crates && rgen.rand_bool()) ? TYPE_CRATE : TYPE_BOX);
		unsigned const obj_ix(objs.size());
		float const side_clearance(0.05*max(sz.x, sz.y)); // small extra clearance on sides for box flaps
		if (!place_obj_along_wall(type, room, sz.z, sz, rgen, zval, room_id, tot_light_amt, place_area, objs_start,
			0.0, 0, 4, 0, gen_box_color(rgen), 0, SHAPE_CUBE, side_clearance)) continue;
		if (type == TYPE_BOX) {check_for_blocked_box_flags(objs, objs_start, obj_ix);} // check for nearby objects that would block the box flaps from opening
	} // for n
}

colorRGBA get_stain_color(rand_gen_t &rgen) {
	colorRGBA color(BLACK); // color.B = 0.0
	color.R = rgen.rand_uniform(0.0, 0.5);
	color.G = rgen.rand_uniform(0.0, 0.5);
	return color;
}
void building_t::add_stains_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	assert(has_room_geom());
	bool const backrooms(room.is_backrooms()), parking_garage(room.is_parking());
	unsigned const max_stains((backrooms ? 20 : (parking_garage ? 5 : 2))); // 2 for regular rooms, 20 for backrooms, and 5 for parking garages
	unsigned const num_floor_stains(rgen.rand() % (max_stains+1));
	unsigned const num_wall_stains((backrooms || parking_garage) ? 0 : (rgen.rand() % (2*max_stains+1))); // no backrooms or parking garage wall stains - no interior walls
	if (num_floor_stains == 0 && num_wall_stains == 0) return;
	float const window_vspacing(get_window_vspace()), height(1.5*get_flooring_thick()), rmax(0.2*window_vspacing);
	cube_t const place_area(get_walkable_room_bounds(room));
	float const place_area_xymin(min(place_area.dx(), place_area.dy()));
	if (rmax > 0.25*place_area_xymin) return; // room is too small
	// stains on the floor
	add_floor_stains(rgen, place_area, zval, room_id, tot_light_amt, objs_start, num_floor_stains, rmax);

	// black mold stains on the ceiling if wood or tile
	if (!room.is_ext_basement() && !has_parking_garage && (!is_house || has_basement_pipes)) {
		float const bot_zval(zval + get_floor_ceil_gap() - height), ceil_rmax(min(0.5f*window_vspacing, 0.4f*place_area_xymin)); // larger radius
		unsigned const num_ceil_stains(rgen.rand() % (max_stains+1));

		for (unsigned n = 0; n < num_ceil_stains; ++n) {
			float const radius(ceil_rmax*rgen.rand_uniform(0.5, 1.0));
			point const pos(gen_xy_pos_in_area(place_area, radius, rgen, bot_zval));
			cube_t const c(get_cube_height_radius(pos, radius, height));
			if (is_obj_placement_blocked(c, room, 0)) continue; // for now, just make one random attempt; okay if overlaps a ceiling light or other object
			interior->room_geom->decal_manager.add_blood_or_stain(pos, radius, BLACK, 0, 2, 0); // is_blood=0; -z
		}
	}
	// stains on the walls
	if (num_wall_stains == 0) return;
	float const fc_gap(get_floor_ceil_gap()), wall_thickness(get_wall_thickness());
	if (rmax > 0.3*fc_gap) return; // room is not tall enough (should not happen)
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval));

	for (unsigned n = 0; n < num_wall_stains; ++n) {
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool());
		if (classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) continue; // stains don't blend against exterior walls
		float const radius(rmax*rgen.rand_uniform(0.5, 1.0)), min_wall_spacing(1.2*radius);
		point pos;
		pos[ dim] = place_area.d[dim][dir] + (dir ? -1.0 : 1.0)*height;
		pos[!dim] = rgen.rand_uniform((place_area.d[!dim][0] + min_wall_spacing), (place_area.d[!dim][1] - min_wall_spacing));
		pos.z     = rgen.rand_uniform(zval+min_wall_spacing, zval+0.75*fc_gap-min_wall_spacing); // not too high
		cube_t c; c.set_from_sphere(pos, radius);
		c.d[dim][ dir] = place_area.d[dim][dir];
		c.d[dim][!dir] = pos[dim];
		c.expand_in_dim(dim, 0.5*wall_thickness); // required to intersect doors
		if (overlaps_other_room_obj(c, objs_start, 1))     continue; // check for things like closets; check_all=1 to include blinds
		if (interior->is_blocked_by_stairs_or_elevator(c)) continue; // check stairs and elevators
		bool bad_place(0);

		for (auto const &d : doorways) {
			if (d.get_true_bcube().intersects(c)) {bad_place = 1; break;}
		}
		if (bad_place) continue;
		interior->room_geom->decal_manager.add_blood_or_stain(pos, radius, get_stain_color(rgen), 0, dim, !dir); // is_blood=0
	} // for n
}
void building_t::add_floor_stains(rand_gen_t &rgen, cube_t const &place_area, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num, float rmax) {
	float const height(1.5*get_flooring_thick());

	for (unsigned n = 0; n < num; ++n) {
		float const radius(rmax*rgen.rand_uniform(0.5, 1.0));
		point const pos(gen_xy_pos_in_area(place_area, radius, rgen, zval));
		cube_t const c(get_cube_height_radius(pos, radius, height));
		if (overlaps_other_room_obj(c, objs_start) || is_obj_placement_blocked(c, place_area, 0)) continue; // for now, just make one random attempt
		interior->room_geom->decal_manager.add_blood_or_stain(point(pos.x, pos.y, zval+height), radius, get_stain_color(rgen), 0, 2, 1); // is_blood=0; +z
	}
}

room_object_t get_conduit(bool dim, bool dir, float radius, float wall_pos_dim, float wall_pos_not_dim, float z1, float z2, unsigned room_id) {
	cube_t conduit;
	set_wall_width(conduit, wall_pos_not_dim, radius, !dim);
	conduit.d[dim][ dir] = wall_pos_dim; // flush with wall
	conduit.d[dim][!dir] = conduit.d[dim][dir] + (dir ? -1.0 : 1.0)*2.0*radius;
	set_cube_zvals(conduit, z1, z2);
	return room_object_t(conduit, TYPE_PIPE, room_id, 0, 1, RO_FLAG_NOCOLL, 1.0, SHAPE_CYLIN, LT_GRAY); // vertical
}

cube_t building_t::get_light_switch_bounds(float floor_zval, float wall_edge, float wall_pos, bool dim, bool dir) const {
	float const wall_thickness(get_wall_thickness()), switch_thickness(0.2*wall_thickness);
	cube_t c;
	c.z1() = floor_zval + 0.38*get_window_vspace(); // same for every switch
	c.z2() = c.z1() + 1.8*wall_thickness; // set height
	c.d[dim][ dir] = wall_edge; // flush with wall
	c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*switch_thickness; // expand out a bit to set thickness
	expand_to_nonzero_area(c, switch_thickness, dim);
	set_wall_width(c, wall_pos, 0.5*wall_thickness, !dim);
	return c;
}
void building_t::add_light_switches_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start,
	bool is_ground_floor, bool is_basement, bool is_jail)
{
	float const wall_thickness(get_wall_thickness()), switch_hwidth(0.5*wall_thickness);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	if (min(room_bounds.dx(), room_bounds.dy()) < 8.0*switch_hwidth) return; // room is too small; shouldn't happen
	float const ceil_zval(zval + get_floor_ceil_gap());
	vect_door_stack_t &doorways(get_doorways_for_room(room, zval)); // place light switch next to a door
	std::shuffle(doorways.begin(), doorways.end(), rand_gen_wrap_t(rgen));
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_end(objs.size());
	bool const first_side(rgen.rand_bool());
	vect_door_stack_t ext_doors; // not really door stacks, but we can fill in the data to treat them as such
	unsigned garage_door_mask(0);

	if (is_ground_floor || (!walkways.empty() && !DRAW_CITY_INT_WINDOWS)) { // handle exterior doors on the ground floor, and walkway doors if there are no windows
		cube_t room_exp(room);
		room_exp.expand_by(wall_thickness, wall_thickness, -wall_thickness); // expand in XY and shrink in Z
		set_cube_zvals(room_exp, zval, ceil_zval); // limit to this floor of the room

		for (auto d = doors.begin(); d != doors.end(); ++d) {
			if (!d->is_exterior_door() || d->type == tquad_with_ix_t::TYPE_RDOOR) continue;
			cube_t bc(d->get_bcube());
			if (!room_exp.contains_pt(bc.get_cube_center())) continue;
			if (d->type == tquad_with_ix_t::TYPE_GDOOR) {garage_door_mask |= (1 << ext_doors.size());}
			bool const dim(bc.dy() < bc.dx());
			bc.expand_in_dim(dim, 0.4*wall_thickness); // expand slightly to make it nonzero area
			ext_doors.emplace_back(door_t(bc, dim, 0), 0); // dir=0, first_door_ix=0 because it's unused
		}
	}
	for (unsigned ei = 0; ei < 2; ++ei) { // exterior, interior
		vect_door_stack_t const &cands(ei ? doorways : ext_doors);
		unsigned const max_ls(is_residential() ? 2 : 1); // place up to 2 light switches in this room if it's a residential, otherwise place only 1
		unsigned num_ls(0);

		for (auto i = cands.begin(); i != cands.end() && num_ls < max_ls; ++i) {
			if (!is_house && room.is_ext_basement() && room_bounds.contains_cube_xy(*i)) continue; // skip interior backrooms doors
			if (i->for_jail == 1) continue; // skip jail cell doors
			// check for windows if (real_num_parts > 1)? is it actually possible for doors to be within far_spacing of a window? yes, for office building walkways
			bool const dim(i->dim), dir(i->get_center_dim(dim) > room.get_center_dim(dim));
			float const dir_sign(dir ? -1.0 : 1.0), door_width(i->get_width());
			assert(door_width > 0.0);
			bool const is_gdoor(ei == 0 && (garage_door_mask & (1 << (i - cands.begin()))));
			float const nom_width(is_gdoor ? 0.25*get_doorway_width() : (is_jail ? 0.75 : 1.0)*door_width); // use small door width for garage doors and closer dist for jails
			float const near_spacing(0.25*nom_width), far_spacing(1.25*nom_width); // off to side of door when open
			float const min_wall_spacing(switch_hwidth + (is_gdoor ? 0.5 : 2.0)*wall_thickness); // reduced spacing for garage doors
			cube_t const &wall_bounds(ei ? room_bounds : room); // exterior door should use the original room, not room_bounds
			bool done(0);

			for (unsigned Side = 0; Side < 2 && !done; ++Side) { // try both sides of the doorway
				bool const side(bool(Side) ^ first_side);

				for (unsigned nf = 0; nf < 2; ++nf) { // {near, far}
					float const spacing(nf ? far_spacing : near_spacing), wall_pos(i->d[!dim][side] + (side ? 1.0 : -1.0)*spacing);
					if (wall_pos < room_bounds.d[!dim][0] + min_wall_spacing || wall_pos > room_bounds.d[!dim][1] - min_wall_spacing) continue; // too close to the adjacent wall
					cube_t c(get_light_switch_bounds(zval, wall_bounds.d[dim][dir], wall_pos, dim, dir)), c_test(c); // should have enough thickness for pool tile
					c_test.d[dim][!dir] += dir_sign*wall_thickness; // expand out more so that it's guaranteed to intersect appliances placed near the wall
					if (overlaps_other_room_obj(c_test, objs_start))                     continue;
					if (!is_gdoor && is_obj_placement_blocked(c_test, room, (ei==1), 1)) continue; // inc_open_doors=1/check_open_dir=1 for inside, to avoid placing behind open door
					if (!check_if_placed_on_interior_wall(c, room, dim, dir))            continue; // ensure the switch is on a wall
					if (overlaps_or_adj_int_window(c))                                   continue; // check interior windows
					// if is_basement, and this is an exterior wall, use a non-recessed light switch? but the basement ext wall will never have a doorway; next to basement stairs?
					unsigned flags(RO_FLAG_NOCOLL);

					if (is_house && is_basement && classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) { // house exterior basement wall; non-recessed
						room_object_t const conduit(get_conduit(dim, dir, 0.25*switch_hwidth, c.d[dim][dir], wall_pos, c.z2(), ceil_zval, room_id));

						if (!overlaps_other_room_obj(conduit, objs_start)) {
							objs.push_back(conduit);
							c.d[dim][!dir] += dir_sign*1.0*switch_hwidth; // shift front outward more
							flags |= RO_FLAG_HANGING;
						}
					}
					objs.emplace_back(c, TYPE_SWITCH, room_id, dim, dir, flags, 1.0); // dim/dir matches wall; fully lit
					done = 1; // done, only need to add one for this door
					++num_ls;
					break;
				} // for nf
			} // for side
		} // for i
	} // for ei
	if (!is_residential() || is_basement) return; // no closets - done

	// add closet light switches
	for (unsigned i = objs_start; i < objs_end; ++i) { // can't iterate over objs while modifying it
		room_object_t const &obj(objs[i]);
		if (obj.type != TYPE_CLOSET) continue;
		cube_t cubes[5]; // front left, left side, front right, right side, door
		get_closet_cubes(obj, cubes); // for_collision=0
		bool const dim(obj.dim), dir(!obj.dir);
		bool side_of_door(0);
		if (obj.is_small_closet()) {side_of_door = 1;} // same side as door handle
		else { // large closet, put the light switch on the side closer to the center of the room
			float const room_center(room.get_center_dim(!dim));
			side_of_door = (fabs(cubes[2].get_center_dim(!dim) - room_center) < fabs(cubes[0].get_center_dim(!dim) - room_center));
		}
		cube_t const &target_wall(cubes[2*side_of_door]); // front left or front right
		cube_t c(get_light_switch_bounds(zval, target_wall.d[dim][!dir], target_wall.get_center_dim(!dim), dim, dir));
		// since nothing is placed against the exterior wall of the closet near the door (to avoid blocking it), we don't need to check for collisions with room objects
		objs.emplace_back(c, TYPE_SWITCH, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_IN_CLOSET), 1.0); // dim/dir matches wall; fully lit; flag for closet
		//break; // there can be only one closet per room; done (unless I add multiple closets later?)
	} // for i
}

void building_t::add_outlets_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id,
	unsigned objs_start, bool is_ground_floor, bool is_basement, bool is_kitchen)
{
	float const wall_thickness(get_wall_thickness()), plate_height(1.8*wall_thickness), plate_hwidth(0.5*wall_thickness), min_wall_spacing(4.0*plate_hwidth);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	if (min(room_bounds.dx(), room_bounds.dy()) < 3.0*min_wall_spacing) return; // room is too small; shouldn't happen
	float const outlet_z1(zval + get_trim_height() + 0.4*plate_height); // wall trim height + some extra padding; same for every outlet
	float plate_thickness(0.03*wall_thickness);
	if ((int)room_id == interior->pool.room_ix) {plate_thickness += 0.5*get_flooring_thick();} // add over pool tile, somewhat recessed
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval));
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const num_outlets_per_wall(is_kitchen ? 2 : 1); // more outlets in kitchen

	// try to add an outlet to each wall, down near the floor so that they don't intersect objects such as pictures
	for (unsigned wall = 0; wall < 4; ++wall) {
		bool const dim(wall >> 1), dir(wall & 1);
		if (!is_commercial() && room.get_sz_dim(!dim) < room.get_sz_dim(dim)) continue; // only add outlets to the long walls of office building rooms
		bool const is_exterior_wall(classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT); // includes basement
		if (is_exterior_wall && !is_cube()) continue; // don't place on ext wall if it's not X/Y aligned
		bool const might_have_windows(!is_basement && has_int_windows() && is_exterior_wall);
		if (might_have_windows && is_industrial()) continue; // industrial buildings have irregular sub-room window placement, so skip these exterior walls
		cube_t const &wall_bounds(is_exterior_wall ? room : room_bounds); // exterior wall should use the original room, not room_bounds
		float const wall_face(wall_bounds.d[dim][dir]), dir_sign(dir ? -1.0 : 1.0);

		for (unsigned n = 0; n < num_outlets_per_wall; ++n) {
			float const wall_pos(rgen.rand_uniform((room_bounds.d[!dim][0] + min_wall_spacing), (room_bounds.d[!dim][1] - min_wall_spacing)));
			cube_t c;
			set_cube_zvals(c, outlet_z1, (outlet_z1 + plate_height));
			set_wall_width(c, wall_pos, plate_hwidth, !dim);
			c.d[dim][ dir] = wall_face; // flush with wall
			c.d[dim][!dir] = wall_face + dir_sign*plate_thickness; // expand out a bit

			if (might_have_windows) { // check for window intersection
				cube_t const &part(get_part_for_room(room));
				float const window_hspacing(get_hspacing_for_part(part, !dim)), window_h_border(get_window_h_border());
				// expand by the width of the window trim, plus some padded wall plate width, then check to the left and right;
				// 2*xy_expand should be smaller than a window so we can't have a window fit in between the left and right sides
				float const xy_expand(get_trim_thickness() + 1.2f*plate_hwidth);
				if (is_val_inside_window(part, !dim, (wall_pos - xy_expand), window_hspacing, window_h_border) ||
					is_val_inside_window(part, !dim, (wall_pos + xy_expand), window_hspacing, window_h_border)) continue;
			}
			cube_t c_exp(c);
			c_exp.expand_by_xy(wall_thickness);
			bool hit_cabinet(0);
		
			if (overlaps_other_room_obj(c_exp, objs_start, 1)) { // check for things like closets; check_all=1 to include blinds
				if (!is_kitchen) continue;

				for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) { // check for kitchen cabinets
					if (i->type != TYPE_COUNTER && i->type != TYPE_KSINK) continue;
					if (!i->intersects(c_exp)) continue;
					if (i->type == TYPE_KSINK && get_sink_cube(*i).intersects_xy(c_exp)) break; // don't place behind the kitchen sink
					float new_zval(i->z2());
					if (i->has_extra() && (i->flags & RO_FLAG_ADJ_BOT)) {new_zval += BACKSPLASH_HEIGHT*i->dz();} // place above the backsplash if present
					c.z1() = c_exp.z1() = new_zval + 0.1*plate_height; // slightly above counter
					c.z2() = c_exp.z2() = c.z1() + plate_height;
					hit_cabinet = 1;
					break;
				} // for i
				if (!hit_cabinet) continue; // hit something else
				if (overlaps_other_room_obj(c_exp, objs_start, 1)) continue;
			}
			if (interior->is_blocked_by_stairs_or_elevator(c_exp)) continue; // check stairs and elevators
			if (!check_cube_within_part_sides(c_exp))              continue; // handle non-cube buildings
			if (overlaps_or_adj_int_window(c_exp))                 continue; // check interior windows
			bool bad_place(0);

			if (is_ground_floor || !walkways.empty()) { // handle exterior doors
				for (auto d = doors.begin(); d != doors.end(); ++d) {
					if (!d->is_exterior_door() || d->type == tquad_with_ix_t::TYPE_RDOOR) continue;
					cube_t bc(d->get_bcube());
					bc.expand_in_dim(dim, wall_thickness); // make sure it's nonzero area
					if (bc.intersects(c_exp)) {bad_place = 1; break;}
				}
				if (bad_place) continue;
			}
			for (auto const &d : doorways) {
				if (d.get_true_bcube().intersects(c_exp)) {bad_place = 1; break;}
			}
			if (bad_place) continue;
			if (!hit_cabinet && !check_if_placed_on_interior_wall(c, room, dim, dir)) continue; // ensure the outlet is on a wall; skip for cabinets/backsplashes
			unsigned flags(RO_FLAG_NOCOLL);

			if (is_house && is_basement && is_exterior_wall) { // house exterior basement wall; non-recessed
				room_object_t const conduit(get_conduit(dim, dir, 0.25*plate_hwidth, c.d[dim][dir], wall_pos, c.z2(), (zval + get_floor_ceil_gap()), room_id));

				if (!overlaps_other_room_obj(conduit, objs_start)) {
					objs.push_back(conduit);
					c.d[dim][!dir] += dir_sign*1.2*plate_hwidth; // shift front outward more
					flags |= RO_FLAG_HANGING;
				}
			}
			expand_to_nonzero_area(c, plate_thickness, dim);
			objs.emplace_back(c, TYPE_OUTLET, room_id, dim, dir, flags, 1.0); // dim/dir matches wall; fully lit
		} // for n
	} // for wall
}

cube_t door_base_t::get_open_door_bcube_for_room(cube_t const &room) const {
	bool const dir(get_check_dirs());
	cube_t bcube(get_true_bcube());
	if (door_opens_inward(*this, room)) {bcube.d[!dim][dir] += (dir ? 1.0 : -1.0)*get_width();} // include door fully open position
	return bcube;
}
unsigned door_base_t::get_conn_room(unsigned room_id) const {
	assert(!no_room_conn());
	if (room_id == conn_room[0]) return conn_room[1];
	if (room_id == conn_room[1]) return conn_room[0];
	cout << TXT(room_id) << TXT(conn_room[0]) << TXT(conn_room[1]) << TXT(str()) << endl; // TESTING
	assert(0); // invalid room
	return room_id;
}

bool building_t::add_wall_vent_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool check_for_ducts) {
	float const wall_thickness(get_wall_thickness()), ceiling_zval(zval + get_floor_ceil_gap());
	float height(2.5*wall_thickness), hwidth(2.0*wall_thickness);
	if (room.is_store()) {height *= 2.0; hwidth *= 2.0;} // larger wall vents for mall stores
	float const thickness(0.1*wall_thickness), min_wall_spacing(1.5*hwidth);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	if (min(room_bounds.dx(), room_bounds.dy()) < 3.0*min_wall_spacing) return 0; // room is too small; shouldn't happen
	bool const pref_dim(room.dx() < room.dy()); // shorter dim, to make it less likely to conflict with whiteboards
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	c.z2() = ceiling_zval - 0.1*height;
	c.z1() = c.z2() - height;

	for (unsigned n = 0; n < 100; ++n) { // 100 tries
		bool const dim((n < 10) ? pref_dim : rgen.rand_bool()), dir(rgen.rand_bool());
		if (classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) continue; // skip exterior walls
		float const wall_pos(rgen.rand_uniform((room_bounds.d[!dim][0] + min_wall_spacing), (room_bounds.d[!dim][1] - min_wall_spacing)));
		float const wall_face(room_bounds.d[dim][dir]);
		c.d[dim][ dir] = wall_face; // flush with wall
		c.d[dim][!dir] = wall_face + (dir ? -1.0 : 1.0)*thickness; // expand out a bit
		set_wall_width(c, wall_pos, hwidth, !dim);
		cube_t c_exp(c);
		c_exp.expand_by_xy(0.5*wall_thickness);
		c_exp.d[dim][!dir] += (dir ? -1.0 : 1.0)*hwidth; // add some clearance in front
		if (overlaps_other_room_obj(c_exp, objs_start, 1))     continue; // check for objects; check_all=1 to inc whiteboards; excludes picture frames
		if (interior->is_blocked_by_stairs_or_elevator(c_exp)) continue; // check stairs and elevators
		if (overlaps_or_adj_int_window(c_exp))                 continue; // check interior windows
		cube_t door_test_cube(c_exp);
		door_test_cube.expand_in_dim(!dim, 0.25*hwidth); // not too close to doors
		bool bad_place(0);

		for (auto const &d : doorways) {
			if (d.get_open_door_bcube_for_room(room).intersects(door_test_cube)) {bad_place = 1; break;}
		}
		if (bad_place) continue;
		if (!check_if_placed_on_interior_wall(c, room, dim, dir)) continue; // ensure the vent is on a wall
		if (!check_cube_within_part_sides(c)) continue; // handle non-cube buildings

		if (check_for_ducts) { // if this is a utility room, check to see if we can connect the vent to a furnace with a duct
			assert(objs_start <= objs.size());

			for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
				//if (i->type == TYPE_DUCT) {} // maybe this can be implemented once we add ducts to rooms (rather than only attics)
				if (i->type != TYPE_FURNACE) continue;
				if (i->dim != dim || i->dir == dir) continue; // wrong wall
				bool const side(i->get_center_dim(!dim) < c.get_center_dim(!dim));
				float const duct_wall_shift((dir ? -1.0 : 1.0)*0.6*height), duct_end_shift((side ? -1.0 : 1.0)*wall_thickness);
				cube_t duct(c);
				duct.d[ dim][!dir ] = wall_face + duct_wall_shift; // expand outward
				duct.d[!dim][!side] = i->d[!dim][side]; // flush with the side of the furnace
				cube_t test_cube(duct);
				test_cube.d[!dim][!side] -= duct_end_shift; // shrink slightly so as not to overlap the furnace
				duct     .d[!dim][!side] += duct_end_shift; // extend duct into furnace duct (since there's a gap on the edge of the furnace)
				if (overlaps_other_room_obj(test_cube, objs_start, 1)) continue; // check for objects
				if (is_obj_placement_blocked(duct, room, 1))           continue; // too close to a doorway/stairs/elevator; inc_open=1
				objs.emplace_back(duct, TYPE_DUCT, room_id, !dim, 0, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, DUCT_COLOR);
				c.translate_dim(dim, duct_wall_shift); // place vent on the duct
				break; // only connect to one furnace
			} // for i
		}
		objs.emplace_back(c, TYPE_VENT, room_id, dim, dir, RO_FLAG_NOCOLL, 1.0); // dim/dir matches wall; fully lit
		return 1; // done
	} // for n
	return 0; // failed
}

// what about floor vent?
bool building_t::add_ceil_vent_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start) {
	float const wall_thickness(get_wall_thickness()), thickness(0.1*wall_thickness), hlen(2.0*wall_thickness), hwid(1.25*wall_thickness);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	if (min(room_bounds.dx(), room_bounds.dy()) < 4.0*hlen) return 0; // room is too small; shouldn't happen
	float const ceiling_zval(room.is_single_floor ? (room.z2() - get_fc_thickness()) : (zval + get_floor_ceil_gap()));
	cube_t c, attic_access;
	set_cube_zvals(c, ceiling_zval-thickness, ceiling_zval);

	if (has_attic()) {
		attic_access = interior->attic_access;
		attic_access.z1() -= get_floor_thickness(); // expand downward a bit do that it would intersect a vent below
	}
	for (unsigned n = 0; n < 10; ++n) { // 10 tries
		bool const dim(rgen.rand_bool());
		point sz;
		sz[ dim] = hlen;
		sz[!dim] = hwid;
		point const center(gen_xy_pos_in_area(room_bounds, sz, rgen));
		set_wall_width(c, center[ dim], hlen,  dim);
		set_wall_width(c, center[!dim], hwid, !dim);
		cube_t c_exp(c);
		c_exp.expand_by_xy(0.5*wall_thickness); // add a bit of padding
		if (has_bcube_int(c_exp, interior->stairwells) || has_bcube_int(c_exp, interior->elevators)) continue; // check stairs and elevators
		if (overlaps_other_room_obj(c_exp, objs_start, 1)) continue; // check for things like closets; check_all=1 to inc whiteboards; excludes picture frames
		if (is_cube_close_to_doorway(c, room, 0.0, 1, 1))  continue;
		if (vent_in_attic_test(c, dim) == 2)               continue; // not enough clearance in attic for duct
		if (has_attic() && c.intersects(attic_access))     continue; // check attic access door
		interior->room_geom->objs.emplace_back(c, TYPE_VENT, room_id, dim, 0, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), 1.0); // dir=0; fully lit
		return 1; // done
	} // for n
	return 0; // failed
}

bool building_t::check_if_placed_on_interior_wall(cube_t const &c, room_t const &room, bool dim, bool dir) const {
	bool const mall_or_store(room.is_mall_or_store());
	if (!has_small_part && !mall_or_store && !room.has_open_wall(dim, dir) && !room_might_have_clipped_wall(room)) return 1; // check not needed, any non-door loc is a wall
	float const wall_thickness(get_wall_thickness()), wall_face(c.d[dim][dir]);
	cube_t test_cube(c);
	test_cube.d[dim][0] = test_cube.d[dim][1] = wall_face - (dir ? -1.0 : 1.0)*0.5*wall_thickness; // move inward
	test_cube.expand_in_dim(!dim, 0.5*wall_thickness);

	if (mall_or_store) { // special logic to handle mall store windows and doors
		if (has_bcube_int(test_cube, interior->int_windows              )) return 0; // on window, not wall
		if (has_bcube_int(test_cube, interior->mall_info->store_doorways)) return 0; // on store door, not wall
		return 1;
	}
	// check for exterior wall
	bool intersects_part(0);

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		if (p->intersects(test_cube)) {intersects_part = 1; break;}
	}
	if (!intersects_part) return 1; // not contained in a part, must be an exterior wall or extended basement
	// check for interior wall
	for (auto const &w : interior->walls[dim]) {
		if (w.contains_cube(test_cube)) return 1;
	}
	return 0;
}

bool building_t::place_eating_items_on_table(rand_gen_t &rgen, unsigned table_obj_id) {
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(table_obj_id < objs.size());
	room_object_t const table(objs[table_obj_id]); // deep copy to avoid invalidating the reference
	float const floor_spacing(get_window_vspace()), plate_radius(get_plate_radius(rgen, table, floor_spacing)), plate_height(0.1*plate_radius), spacing(1.33*plate_radius);
	unsigned const objs_size(objs.size());
	bool added_obj(0);

	for (unsigned i = (table_obj_id + 1); i < objs_size; ++i) { // iterate over chairs (by index, since we're adding to objs)
		if (objs[i].type != TYPE_CHAIR) break; // done with chairs for this table
		point const chair_center(objs[i].get_cube_center()), table_center(table.get_cube_center());
		point pos;

		if (table.shape == SHAPE_CYLIN) { // circular
			float const dist(table.get_radius() - spacing);
			pos = table_center + dist*(chair_center - table_center).get_norm();
		}
		else { // rectangular
			cube_t place_bounds(table);
			place_bounds.expand_by_xy(-spacing);
			pos = place_bounds.closest_pt(chair_center);
		}
		unsigned const plate_flags(RO_FLAG_NOCOLL | (table.is_glass_table() ? RO_FLAG_ADJ_BOT : 0)); // flag for bottom to be drawn if on a glass table
		cube_t plate;
		plate.set_from_sphere(pos, plate_radius);
		set_cube_zvals(plate, table.z2(), table.z2()+plate_height); // place on the table
		objs.emplace_back(plate, TYPE_PLATE, table.room_id, 0, 0, plate_flags, table.light_amt, SHAPE_CYLIN);
		set_obj_id(objs);

		if (building_obj_model_loader.is_model_valid(OBJ_MODEL_SILVER)) {
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_SILVER)); // D, W, H
			float const sw_height(0.0075*floor_spacing), sw_hwidth(0.5*sw_height*sz.x/sz.z), sw_hlen(0.5*sw_height*sz.y/sz.z); // Note: x/y swapped
			vector3d const offset(pos - table_center);
			bool const dim(objs[i].dim), dir(!objs[i].dir); // use chair dim/dir
			cube_t sw_bc;
			set_cube_zvals(sw_bc, table.z2()+0.1*sw_height, table.z2()+sw_height);
			set_wall_width(sw_bc, pos[!dim] + ((dim ^ dir) ? 1.0 : -1.0)*1.2*(plate_radius + sw_hlen), sw_hlen, !dim);
			set_wall_width(sw_bc, pos[ dim], sw_hwidth, dim);
			objs.emplace_back(sw_bc, TYPE_SILVER, table.room_id, dim, dir, RO_FLAG_NOCOLL, table.light_amt, SHAPE_CUBE, GRAY);
		}
		added_obj = 1;
	} // for i
	if (added_obj && !is_school()) { // place a vase in the center of the table, but not for schools
		float const vase_radius(rgen.rand_uniform(0.35, 0.6)*plate_radius), vase_height(rgen.rand_uniform(2.0, 6.0)*vase_radius);
		cube_t vase;
		vase.set_from_sphere(table.get_cube_center(), vase_radius);
		set_cube_zvals(vase, table.z2(), table.z2()+vase_height); // place on the table
		objs.emplace_back(vase, TYPE_VASE, table.room_id, 0, 0, RO_FLAG_NOCOLL, table.light_amt, SHAPE_CYLIN, gen_vase_color(rgen));
		set_obj_id(objs);
	}
	return added_obj;
}

void building_t::place_objects_onto_surfaces(rand_gen_t rgen, room_t const &room, unsigned room_id, float tot_light_amt,
	unsigned objs_start, unsigned floor, bool is_basement, bool not_private)
{
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(objs.size() > objs_start);
	bool const is_library   (room.get_room_type(floor) == RTYPE_LIBRARY);
	bool const is_kitchen   (room.get_room_type(floor) == RTYPE_KITCHEN);
	bool const is_conference(room.get_room_type(floor) == RTYPE_CONF   );
	bool const sparse_place(floor > 0 && !is_conference && interior->rooms.size() > 40); // fewer objects on upper floors of large office buildings as an optimization
	float const place_book_prob(( is_house ? 1.0 : 0.5)*(room.is_office ? 0.80 : 1.00)*(sparse_place ? 0.75 : 1.0));
	float const place_bottle_prob(is_house ? 1.0 :      (room.is_office ? 0.80 : 0.50)*(sparse_place ? 0.50 : 1.0));
	float const place_cup_prob   (is_house ? 1.0 :      (room.is_office ? 0.50 : 0.25)*(sparse_place ? 0.50 : 1.0));
	float const place_plant_prob (is_house ? 1.0 :      (room.is_office ? 0.25 : 0.15)*(sparse_place ? 0.75 : 1.0));
	float const place_laptop_prob(is_house ? 0.4 :      (room.is_office ? 0.60 : 0.50)*(sparse_place ? 0.80 : 1.0));
	float const place_pizza_prob (is_house ? 1.0 :      (room.is_office ? 0.30 : 0.15)*(sparse_place ? 0.75 : 1.0));
	float const place_banana_prob(is_house ? 1.0 :      (room.is_office ? 0.40 : 0.20)*(sparse_place ? 0.50 : 1.0));
	unsigned const objs_end(objs.size());
	bool placed_book_on_counter(0);

	// see if we can place objects on any room object top surfaces
	for (unsigned i = objs_start; i < objs_end; ++i) { // can't iterate over objs because we modify it
		room_object_t const &obj(objs[i]);
		if (obj.is_on_floor()) continue; // can't place on fallen over object
		// add place settings to kitchen and dining room tables 50% of the time
		bool const is_table(obj.type == TYPE_TABLE || obj.type == TYPE_CONF_TABLE); // for eating, and directionless (meaning objects can have any orient)
		bool const is_eating_table(is_table && (room.get_room_type(floor) == RTYPE_KITCHEN || room.get_room_type(floor) == RTYPE_DINING) && rgen.rand_bool());
		if (is_eating_table && place_eating_items_on_table(rgen, i)) continue; // no other items to place
		float book_prob(0.0), bottle_prob(0.0), cup_prob(0.0), plant_prob(0.0), laptop_prob(0.0), pizza_prob(0.0), toy_prob(0.0), banana_prob(0.0);
		static vect_cube_t avoid; // reuse across buildings
		avoid.clear();

		if (obj.type == TYPE_TABLE && i == objs_start) { // only first table (not TV table)
			book_prob   = 0.4*place_book_prob;
			bottle_prob = 0.6*place_bottle_prob;
			cup_prob    = 0.5*place_cup_prob;
			plant_prob  = 0.6*place_plant_prob;
			laptop_prob = 0.3*place_laptop_prob;
			pizza_prob  = 0.8*place_pizza_prob;
			banana_prob = 0.7*place_banana_prob;
			if ((is_house || (is_apartment() && !not_private)) && !is_kitchen) {toy_prob = 0.5;} // toys are in houses and private apartments rooms only; not on kitchen tables
		}
		else if (obj.type == TYPE_CONF_TABLE) {
			book_prob   = 0.3*place_book_prob;
			bottle_prob = 0.4*place_bottle_prob;
			cup_prob    = 0.6*place_cup_prob;
			laptop_prob = 0.8*place_laptop_prob;
			pizza_prob  = 1.0*place_pizza_prob;
			if (i+1 < objs_end && objs[i+1].type == TYPE_CONF_PHONE) {avoid.push_back(objs[i+1]);} // avoid overlapping the conference phone
		}
		else if (obj.type == TYPE_DESK && !(obj.flags & RO_FLAG_ADJ_TOP)) { // desk with no computer monitor
			book_prob   = 0.8*place_book_prob;
			bottle_prob = 0.4*place_bottle_prob;
			cup_prob    = 0.3*place_cup_prob;
			plant_prob  = 0.3*place_plant_prob;
			laptop_prob = 0.7*place_laptop_prob;
			pizza_prob  = 0.4*place_pizza_prob;
			banana_prob = 0.2*place_banana_prob;
		}
		else if (obj.type == TYPE_RDESK) { // reception desk
			book_prob   = 0.4*place_book_prob;
			bottle_prob = 0.6*place_bottle_prob;
			cup_prob    = 0.7*place_cup_prob;
			plant_prob  = 1.0*place_plant_prob;
		}
		else if (obj.type == TYPE_COUNTER && !(obj.flags & RO_FLAG_ADJ_TOP)) { // counter without a microwave
			book_prob   = (placed_book_on_counter ? 0.0 : 0.5); // only place one book per counter
			bottle_prob = 0.25*place_bottle_prob;
			cup_prob    = 0.30*place_cup_prob;
			plant_prob  = 0.10*place_plant_prob;
			laptop_prob = 0.05*place_laptop_prob;
			pizza_prob  = 0.50*place_pizza_prob;
			banana_prob = 0.90*place_banana_prob;
		}
		else if ((obj.type == TYPE_DRESSER || obj.type == TYPE_NIGHTSTAND) && !(obj.flags & RO_FLAG_ADJ_TOP)) { // dresser or nightstand with nothing on it yet; no pizza
			book_prob   = 0.25*place_book_prob;
			bottle_prob = 0.15*place_bottle_prob;
			cup_prob    = 0.15*place_cup_prob;
			plant_prob  = 0.1*place_plant_prob;
			laptop_prob = 0.1*place_laptop_prob;
			toy_prob    = 0.15;
		}
		else {
			continue;
		}
		if (is_library) {book_prob *= 2.5;} // higher probability of books placed in a library
		if (is_kitchen) {cup_prob  *= 2.0; pizza_prob *= 2.0;} // higher probability of cups and pizza if placed in a kitchen
		room_object_t surface(obj); // deep copy to allow modification and avoid using an invalidated reference
		unsigned sel_ix(0); // for reception desks
		
		if (obj.shape == SHAPE_CYLIN) { // find max contained XY rectangle (simpler than testing distance to center vs. radius)
			for (unsigned d = 0; d < 2; ++d) {surface.expand_in_dim(d, -0.5*(1.0 - SQRTOFTWOINV)*surface.get_sz_dim(d));}
		}
		else if (obj.type == TYPE_RDESK) { // select a random cube
			surface.copy_from(get_rand_reception_desk_placement_cube(obj, rgen, sel_ix));
			if (surface.is_all_zeros()) continue; // no surface to place on
		}
		// Note: after this point, the obj reference is invalid
		if (is_eating_table) { // table in a room for eating, add a plate
			if (place_plate_on_obj(rgen, surface, room_id, tot_light_amt, avoid)) {avoid.push_back(objs.back());}
		}
		if (avoid.empty() && rgen.rand_probability(book_prob)) { // place book if it's the first item (no plate)
			placed_book_on_counter |= (surface.type == TYPE_COUNTER);
			place_book_on_obj(rgen, surface, room_id, tot_light_amt, objs_start, !is_table);
			avoid.push_back(objs.back());
		}
		if (surface.type == TYPE_DESK) { // if this is a desk, try to avoid placing an object that overlaps a pen or pencil; papers are okay to overlap since they're flat
			for (unsigned j = i+1; j < objs_end; ++j) {
				room_object_t const &obj2(objs[j]);
				if (obj2.type == TYPE_PEN || obj2.type == TYPE_PENCIL) {avoid.push_back(obj2);}
			}
		}
		unsigned const num_obj_types = 7;
		unsigned const obj_type_start(rgen.rand() % num_obj_types); // select a random starting point to remove bias toward objects checked first
		bool placed(0);

		for (unsigned n = 0; n < num_obj_types && !placed; ++n) { // place a single object; ***Note: obj is invalidated by this loop and can't be used***
			switch ((n + obj_type_start) % num_obj_types) {
			case 0: placed = (rgen.rand_probability(bottle_prob) && place_bottle_on_obj(rgen, surface, room_id, tot_light_amt, avoid)); break;
			case 1: placed = (rgen.rand_probability(cup_prob   ) && place_cup_on_obj   (rgen, surface, room_id, tot_light_amt, avoid)); break;
			case 2: placed = (rgen.rand_probability(laptop_prob) && place_laptop_on_obj(rgen, surface, room_id, tot_light_amt, avoid, !is_table)); break;
			case 3: placed = (rgen.rand_probability(pizza_prob ) && place_pizza_on_obj (rgen, surface, room_id, tot_light_amt, avoid)); break;
			case 4: placed = (rgen.rand_probability(banana_prob) && place_banana_on_obj(rgen, surface, room_id, tot_light_amt, avoid)); break;
			case 5: placed = (!is_basement && rgen.rand_probability(plant_prob) && place_plant_on_obj(rgen, surface, room_id, tot_light_amt, 0.7, avoid)); break; // sz_scale=0.7
			case 6: placed = (rgen.rand_probability(toy_prob)    && place_toy_on_obj   (rgen, surface, room_id, tot_light_amt, avoid)); break;
			}
		} // for n
		if (placed && objs[i].type == TYPE_RDESK) {objs[i].item_flags |= (1 << sel_ix);} // mark this surface of the reception desk as occupied
	} // for i
}

template<typename T> bool any_cube_contains(cube_t const &cube, T const &cubes) {
	for (auto const &c : cubes) {if (c.contains_cube(cube)) return 1;}
	return 0;
}
bool building_t::is_light_placement_valid(cube_t const &light, room_t const &room, float pad) const {
	if (room.is_industrial()) return 1; // industrial lights hang from the ceiling; assume they are placed correctly
	cube_t light_ext(light);
	light_ext.expand_by_xy(pad);
	if (!room.contains_cube(light_ext))             return 0; // room too small?
	if (has_bcube_int(light, interior->elevators )) return 0;
	if (has_bcube_int(light, interior->escalators)) return 0; // conservative; is this needed?
	if (!check_cube_within_part_sides(light))       return 0; // handle non-cube buildings
	if (!interior->elevator_equip_room.is_all_zeros() && interior->elevator_equip_room.intersects(light)) return 0;

	// the fc_occluders test below will handle stairs cutouts, but we still need to handle the ceiling above stairs;
	// check stairs with walled sides because these may clip through ceiling lights; lights completely contained in the stairs look okay and are allowed
	for (stairwell_t const &s : interior->stairwells) {
		if (s.intersects(light) && !s.contains_cube_xy(light) && (s.has_walled_sides() || s.is_u_shape())) return 0;
	}
	unsigned const pg_wall_start(interior->room_geom->wall_ps_start);

	// check for intersection with low pipes such as sprinkler pipes that have been previously placed; only works for top level of parking garage; skip for backrooms
	if (light.z1() < ground_floor_z1 && has_parking_garage && get_basement().contains_cube(light) && pg_wall_start > 0) {
		vect_room_object_t const &objs(interior->room_geom->objs);
		assert(pg_wall_start < objs.size());

		for (auto i = objs.begin()+pg_wall_start; i != objs.end(); ++i) {
			if (i->type == TYPE_PIPE && i->intersects(light)) return 0; // check for pipe intersections (in particular horizontal sprinkler pipes)
		}
	}
	if (has_tall_retail() && room.is_retail()) { // handle lights blocked by retail stairs
		for (stairwell_t const &s : interior->stairwells) {
			if (s.z2() < room.z2() || !s.intersects(room)) continue; // wrong stairs
			float const landing_width(s.get_retail_landing_width(get_window_vspace()));
			cube_t stairs_ext(s);
			stairs_ext.d[s.dim][!s.dir] += (s.dir ? -1.0 : 1.0)*(landing_width + get_wall_thickness()); // extend outward to include the landing
			if (stairs_ext.intersects(light)) return 0;
		} // for s
	}
	// handle parking structure bathroom; skip if this room is the bathroom
	if (!interior->ps_bathroom.is_all_zeros() && light_ext.intersects(interior->ps_bathroom) && !interior->ps_bathroom.contains_cube(room)) return 0;
	light_ext.z1() = light_ext.z1() = light.z2() + get_fc_thickness(); // shift in between the ceiling and floor so that we can do a cube contains check
	if (any_cube_contains(light_ext, interior->fc_occluders)) return 1; // Note: don't need to check skylights because fc_occluders excludes skylights
	if (PLACE_LIGHTS_ON_SKYLIGHTS && any_cube_contains(light_ext, skylights)) return 1; // place on a skylight
	return 0;
}

void building_t::try_place_light_on_ceiling(cube_t const &light, room_t const &room, unsigned room_id, bool room_dim, float pad,
	bool allow_rot, bool allow_mult, unsigned nx, unsigned ny, unsigned check_coll_start, vect_cube_t &lights, rand_gen_t &rgen) const
{
	assert(has_room_geom());
	float const window_vspacing(get_window_vspace());
	int light_placed(0); // 0=no, 1=at center, 2=at alternate location

	if (is_hospital() && room.has_subroom()) { // avoid placing light inside bathroom
		// hopefully a quadratic iteration over rooms is okay here, since this is only for the main building; there should be less than 100 rooms in total
		for (auto r = interior->rooms.begin(); r != interior->rooms.begin()+room_id; ++r) {
			if (!r->is_nested() || !room.contains_cube(*r)) continue;
			cube_t sub_room_inc_walls(*r);
			sub_room_inc_walls.expand_by_xy(0.5*get_wall_thickness());
			if (light.intersects_xy(sub_room_inc_walls)) return; // skip this light
			break; // okay for now, since there's only one nested room
		}
	}
	if (is_light_placement_valid(light, room, pad) && !overlaps_other_room_obj(light, check_coll_start)) {
		lights.push_back(light); // valid placement, done
		light_placed = 1;
	}
	else {
		point room_center(room.get_cube_center());
		bool const first_dir(rgen.rand_bool()); // randomize shift direction
		cube_t light_cand(light); // Note: same logic for cube and cylinder light shape - cylinder uses bounding cube
		unsigned const num_shifts = 10;

		if (allow_rot) { // flip aspect ratio
			float const sz_diff(0.5*(light.dx() - light.dy()));
			light_cand.expand_in_x(-sz_diff);
			light_cand.expand_in_y( sz_diff);
		}
		for (unsigned D = 0; D < 2 && !light_placed; ++D) { // try both dims
			bool const dim(room_dim ^ bool(D));
			unsigned const num(room_dim ? ny : nx);
			float const shift_step((0.5*(room.get_sz_dim(dim) - light_cand.get_sz_dim(dim)))/(num*num_shifts)); // shift within the bounds of placement grid based on num

			for (unsigned d = 0; d < 2; ++d) { // dir: see if we can place it by moving on one direction
				for (unsigned n = 1; n <= num_shifts; ++n) { // try different shift values
					cube_t cand(light_cand);
					cand.translate_dim(dim, ((bool(d) ^ first_dir) ? -1.0 : 1.0)*n*shift_step);
					if (!is_light_placement_valid(cand, room, pad))        continue;
					if (overlaps_other_room_obj(cand, check_coll_start))   continue; // intersects wall, pillar, etc.
					if (lights.size() <= 8 && has_bcube_int(cand, lights)) continue; // check for intersecting other lights if there aren't too many
					lights.push_back(cand);
					light_placed = 2;
					break;
				} // for n
				if (!allow_mult && light_placed) break;
			} // for d
		} // for D
	}
	if (light_placed) {
		cube_t &cur_light(lights.back());
		
		// check doors if placed off-center or centered but close to the room bounds; only needed for non-centered lights, lights in small rooms, and backrooms
		if (light_placed == 2 || room.is_backrooms() || !room.contains_cube_xy_exp(cur_light, get_doorway_width())) {
			cube_t test_cube(cur_light);
			test_cube.z1() -= 0.4*window_vspacing; // lower Z1 so that it's guaranteed to overlap a door

			// maybe should exclude basement doors, since they don't show as open? but then it would be wrong if I later draw basement doors;
			// note that this test is conservative for cylindrical house lights
			if (is_cube_close_to_doorway(test_cube, room, 0.0, 1, 1)) { // inc_open=1, check_open_dir=1
				float const orig_z1(cur_light.z1());
				cur_light.z1() += 0.98*cur_light.dz(); // if light intersects door, move it up into the ceiling rather than letting it hang down into the room
				if (cur_light.z1() == cur_light.z2()) {cur_light.z1() = orig_z1;} // fix to avoid zero area cube assert due to FP error
			}
		}
	}
	// else place light on a wall instead? do we ever get here?
}

void building_t::try_place_light_on_wall(cube_t const &light, room_t const &room, bool room_dim, float zval, vect_cube_t &lights, rand_gen_t &rgen) const {
	float const wall_thickness(get_wall_thickness()), window_vspacing(get_window_vspace());
	float const length(light.dz()), radius(0.25*min(light.dx(), light.dy())), min_wall_spacing(2.0*radius);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	if (min(room_bounds.dx(), room_bounds.dy()) < 3.0*min_wall_spacing) return; // room is too small; shouldn't happen
	bool const pref_dim(!room_dim); // shorter dim
	float const light_ws_dmax = 20.0; // max distance in world space from the origin; used to prevent shadow acne when values are large
	point const lpos(light.get_cube_center());
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval)); // must use floor zval
	cube_t c;
	c.z2() = light.z2() - 0.1*window_vspacing;
	c.z1() = c.z2() - 2.0*radius;

	for (unsigned n = 0; n < 100; ++n) { // 100 tries
		bool const dim((n < 10) ? pref_dim : rgen.rand_bool()), dir(rgen.rand_bool());
		if (dir ? (lpos[dim] < -light_ws_dmax) : (lpos[dim] > light_ws_dmax)) continue; // too far from the origin
		float const wall_edge_spacing(max(min_wall_spacing, 0.25f*room_bounds.get_sz_dim(!dim))); // not near a corner of the room
		float const wall_pos(rgen.rand_uniform((room_bounds.d[!dim][0] + wall_edge_spacing), (room_bounds.d[!dim][1] - wall_edge_spacing)));
		float const wall_face(room_bounds.d[dim][dir]);
		c.d[dim][ dir] = wall_face; // flush with wall
		c.d[dim][!dir] = wall_face + (dir ? -1.0 : 1.0)*length; // expand outward
		set_wall_width(c, wall_pos, radius, !dim);
		cube_t c_exp(c);
		c_exp.expand_by_xy(0.5*wall_thickness);
		c_exp.d[dim][!dir] += (dir ? -1.0 : 1.0)*2.0*(length + radius); // add some clearance in front
		if (interior->is_blocked_by_stairs_or_elevator(c_exp)) continue; // check stairs and elevators; no other room objects have been placed yet
		cube_t door_test_cube(c_exp);
		door_test_cube.expand_in_dim(!dim, 1.0*radius); // not too close to doors
		bool bad_place(0);

		for (auto const &d : doorways) {
			if (d.get_open_door_bcube_for_room(room).intersects(door_test_cube)) {bad_place = 1; break;}
		}
		if (bad_place) continue;
		if (!check_if_placed_on_interior_wall(c, room, dim, dir)) continue; // ensure the light is on a wall; is this really needed?
		if (!check_cube_within_part_sides(c)) continue; // handle non-cube buildings
		lights.push_back(c);
		break; // success/done
	} // for n
}

bool building_t::add_padlock_to_door(unsigned door_ix, unsigned lock_color_mask, rand_gen_t &rgen) {
	if (lock_color_mask == 0 || !has_room_geom() || !building_obj_model_loader.is_model_valid(OBJ_MODEL_PADLOCK)) return 0;
	door_t &door(get_door(door_ix));
	if (door.open) return 0; // door is not closed
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned color_ix(0);

	if (door.is_padlocked()) {color_ix = door.get_padlock_color_ix();} // already has a padlock (from a previous room geom gen), use the same color
	else { // select a lock from the colors available from keys found in this building
		assert(door.obj_ix < 0); // not yet assigned
		static vector<unsigned> avail_colors;
		avail_colors.clear();

		for (unsigned n = 0; n < NUM_LOCK_COLORS; ++n) {
			if (lock_color_mask & (1 << n)) {avail_colors.push_back(n);}
		}
		color_ix = avail_colors[(avail_colors.size() == 1) ? 0 : (rgen.rand() % avail_colors.size())];
	}
	door.obj_ix = objs.size();
	door.set_padlock_color_ix(color_ix); // padlocked
	cube_t const door_bc(door.get_true_bcube());
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_PADLOCK)); // D, W, H
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	float const door_width(door_bc.get_sz_dim(!door.dim)), door_height(door_bc.dz());
	float const height((door.for_jail ? 0.09 : 0.078)*door_height), hwidth(0.5*height*sz.y/sz.z), depth(height*sz.x/sz.z);
	float const edge_dist(door.for_jail ? 0.088 : 0.062);
	bool const side(door.get_check_dirs());
	cube_t lock;
	lock.z1() = door.z1() + (door.for_jail ? 0.365 : 0.41)*door_height;
	lock.z2() = lock.z1() + height;
	set_wall_width(lock, (door_bc.d[!door.dim][0] + (side ? edge_dist : (1.0 - edge_dist))*door_width), hwidth, !door.dim);
	colorRGBA const color(lock_colors[color_ix]);

	// since we don't know which side of the door the player will be on, add the padlock to both sides
	for (unsigned d = 0; d < 2; ++d) {
		if (door.for_jail && bool(d) == door.open_dir) continue; // no padlock on the inside/cell side of the door
		float const pos(door_bc.d[door.dim][!d]);
		lock.d[door.dim][ d] = pos;
		lock.d[door.dim][!d] = pos + (d ? -1.0 : 1.0)*depth;
		point center(lock.get_cube_center());
		center[door.dim] += (d ? -1.0 : 1.0)*wall_thickness; // make sure the point is fully inside the room and not inside the wall area/doorway
		int const room_id(get_room_containing_pt(center));
		assert(room_id >= 0); // must be found
		// if the room has a single interior door, no exterior door, and no stairs, then the player can't be inside it unless it's unlocked, and we can skip the interior lock
		room_t const &room(get_room(room_id));
		unsigned const floor_ix(room.get_floor_containing_zval(lock.z1(), floor_spacing));
		if (!room.has_stairs_on_floor(floor_ix) && !(floor_ix == 0 && is_room_adjacent_to_ext_door(room)) && count_num_int_doors(room) == 1) continue;
		objs.emplace_back(lock, TYPE_PADLOCK, ((room_id < 0) ? 0 : room_id), door.dim, d, (RO_FLAG_NOCOLL | RO_FLAG_IS_ACTIVE), 1.0, SHAPE_CUBE, color); // starts attached
	} // for d
	interior->room_geom->invalidate_model_geom();
	return 1;
}
bool building_t::remove_padlock_from_door(unsigned door_ix, point const &remove_pos) {
	if (!has_room_geom() || !building_obj_model_loader.is_model_valid(OBJ_MODEL_PADLOCK)) return 0;
	door_t &door(get_door(door_ix));
	if (!door.is_padlocked()) return 0; // no padlock
	//assert(!door.open); // too strong?
	door.locked = 0; // unlocked
	bool const dim(door.dim);
	float const center(door.get_center_dim(dim));
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(door.obj_ix >= 0 && unsigned(door.obj_ix+1) < objs.size()); // must be space for two locks
	bool saw_padlock(0);

	for (unsigned d = 0; d < 2; ++d) {
		room_object_t &obj(objs[door.obj_ix + d]);
		if (obj.type != TYPE_PADLOCK) continue; // may have been skipped above if the door has a single entry
		saw_padlock = 1;
		bool const keep(((remove_pos[dim] - center)*(obj.get_center_dim(dim) - center)) > 0.0); // keep if it's on the side that was opened

		if (keep) {
			obj.flags &= ~RO_FLAG_IS_ACTIVE; // no longer attached
			obj.translate_dim(2, (door.z1() - obj.z1())); // move to the floor at the base of the door
		}
		else {obj.remove();}
	} // for d
	assert(saw_padlock); // at least one padlock must have been removed
	interior->room_geom->invalidate_model_geom();
	return 1;
}

