// 3D World - Building Interior Room Geometry Placement
// by Frank Gennari 4/30/2020

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city.h" // for object_model_loader_t
#include "profiler.h"

enum {PLACED_TOILET=1, PLACED_SINK=2, PLACED_TUB=4, PLACED_SHOWER=8}; // for bathroom objects

extern bool camera_in_building;
extern int display_mode;
extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;
extern bldg_obj_type_t bldg_obj_types[];

void setup_bldg_obj_types();
bool enable_parked_cars();
car_t car_from_parking_space(room_object_t const &o);


class light_ix_assign_t {
	vector<pair<point2d<float>, unsigned>> cur;
	unsigned next_ix;
public:
	light_ix_assign_t() : next_ix(0) {}
	void next_room() {cur.clear();}
	unsigned get_next_ix() {return next_ix++;}

	unsigned get_ix_for_light(cube_t const &c) {
		point2d<float> const pos(c.x1(), c.y1());

		for (auto i = cur.begin(); i != cur.end(); ++i) {
			if (i->first == pos) return i->second; // existing light is part of the same stack and is valid to return
		}
		cur.emplace_back(pos, get_next_ix()); // allocate a new light
		return cur.back().second;
	}
};

bool building_t::is_valid_placement_for_room(cube_t const &c, cube_t const &room, vect_cube_t const &blockers, bool inc_open_doors, float room_pad) const {
	cube_t place_area(room);
	if (room_pad != 0.0f) {place_area.expand_by_xy(-room_pad);} // shrink by dmin
	if (!place_area.contains_cube_xy(c)) return 0; // not contained in interior part of the room
	if (is_cube_close_to_doorway(c, room, 0.0, inc_open_doors)) return 0; // too close to a doorway
	if (interior && interior->is_blocked_by_stairs_or_elevator(c)) return 0; // faster to check only one per stairwell, but then we need to store another vector?
	if (has_bcube_int(c, blockers)) return 0; // Note: ignores dmin
	return 1;
}

float get_radius_for_square_model(unsigned model_id) {
	vector3d const chair_sz(building_obj_model_loader.get_model_world_space_size(model_id));
	return 0.5f*(chair_sz.x + chair_sz.y)/chair_sz.z; // assume square and take average of xsize and ysize
}

bool building_t::add_chair(rand_gen_t &rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id, point const &place_pos,
	colorRGBA const &chair_color, bool dim, bool dir, float tot_light_amt, bool office_chair_model)
{
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR)) {office_chair_model = 0;}
	float const window_vspacing(get_window_vspace()), room_pad(4.0f*get_wall_thickness()), chair_height(0.4*window_vspacing);
	float chair_hwidth(0.0), push_out(0.0);
	point chair_pos(place_pos); // same starting center and z1

	if (office_chair_model) {
		chair_hwidth = 0.5f*chair_height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR);
		push_out     = 0.5; // pushed out a bit so that the arms don't intersect the table top
	}
	else {
		chair_hwidth = 0.1*window_vspacing; // half width
		push_out     = rgen.rand_uniform(-0.5, 1.2); // varible amount of pushed in/out
	}
	chair_pos[dim] += (dir ? -1.0f : 1.0f)*push_out*chair_hwidth;
	cube_t chair(get_cube_height_radius(chair_pos, chair_hwidth, chair_height));
	
	if (!is_valid_placement_for_room(chair, room, blockers, 0, room_pad)) { // check proximity to doors
		if (office_chair_model) return 0; // can't push office chair in any more
		float const max_push_in((dir ? -1.0f : 1.0f)*(-0.5 - push_out)*chair_hwidth);
		chair.translate_dim(dim, max_push_in*rgen.rand_uniform(0.5, 1.0)); // push the chair mostly in and try again
		if (!is_valid_placement_for_room(chair, room, blockers, 0, room_pad)) return 0;
	}
	vect_room_object_t &objs(interior->room_geom->objs);

	if (office_chair_model) {
		float const lum(chair_color.get_weighted_luminance()); // calculate grayscale luminance
		objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, dim, dir, 0, tot_light_amt, SHAPE_CUBE, colorRGBA(lum, lum, lum));
	}
	else {
		objs.emplace_back(chair, TYPE_CHAIR, room_id, dim, dir, 0, tot_light_amt, SHAPE_CUBE, chair_color);
	}
	return 1;
}

// Note: must be first placed objects; returns the number of total objects added (table + optional chairs)
unsigned building_t::add_table_and_chairs(rand_gen_t rgen, cube_t const &room, vect_cube_t const &blockers, unsigned room_id,
	point const &place_pos, colorRGBA const &chair_color, float rand_place_off, float tot_light_amt)
{
	float const window_vspacing(get_window_vspace()), room_pad(max(4.0f*get_wall_thickness(), get_min_front_clearance()));
	vector3d const room_sz(room.get_size());
	vect_room_object_t &objs(interior->room_geom->objs);
	point table_pos(place_pos);
	vector3d table_sz;
	for (unsigned d = 0; d < 2; ++d) {table_sz [d]  = 0.18*window_vspacing*(1.0 + rgen.rand_float());} // half size relative to window_vspacing
	for (unsigned d = 0; d < 2; ++d) {table_pos[d] += rand_place_off*room_sz[d]*rgen.rand_uniform(-1.0, 1.0);} // near the center of the room
	bool const is_round((rgen.rand()&3) == 0); // 25% of the time
	if (is_round) {table_sz.x = table_sz.y = 0.6f*(table_sz.x + table_sz.y);} // round tables must have square bcubes for now (no oval tables yet); make radius slightly larger
	point llc(table_pos - table_sz), urc(table_pos + table_sz);
	llc.z = table_pos.z; // bottom
	urc.z = table_pos.z + rgen.rand_uniform(0.20, 0.22)*window_vspacing; // top
	cube_t table(llc, urc);
	if (!is_valid_placement_for_room(table, room, blockers, 0, room_pad)) return 0; // check proximity to doors and collision with blockers
	objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, (is_house ? RO_FLAG_IS_HOUSE : 0), tot_light_amt, (is_round ? SHAPE_CYLIN : SHAPE_CUBE));
	set_obj_id(objs);
	unsigned num_added(1); // start with the table

	// place some chairs around the table
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			if (rgen.rand_bool()) continue; // 50% of the time
			point chair_pos(table_pos); // same starting center and z1
			chair_pos[dim] += (dir ? -1.0f : 1.0f)*table_sz[dim];
			num_added += add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, 0); // office_chair_model=0
		}
	}
	return num_added;
}
void building_t::shorten_chairs_in_region(cube_t const &region, unsigned objs_start) {
	for (auto i = interior->room_geom->objs.begin() + objs_start; i != interior->room_geom->objs.end(); ++i) {
		if (i->type != TYPE_CHAIR || !i->intersects(region)) continue;
		i->z2() -= 0.25*i->dz();
		i->shape = SHAPE_SHORT;
	}
}

vect_door_stack_t &building_t::get_doorways_for_room(room_t const &room, float zval) const { // interior doorways
	// find interior doorways connected to this room
	float const floor_thickness(get_floor_thickness());
	cube_t room_exp(room);
	room_exp.expand_by_xy(get_wall_thickness());
	set_cube_zvals(room_exp, (zval + floor_thickness), (zval + get_window_vspace() - floor_thickness)); // clip to z-range of this floor
	static vect_door_stack_t doorways; // reuse across rooms
	doorways.clear();

	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
		if (i->on_stairs) continue; // skip basement door
		if (i->intersects(room_exp)) {doorways.push_back(*i);}
	}
	return doorways;
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
	center.z = zval + 0.0012*floor_spacing; // slightly above the floor/rug to avoid z-fighting
	unsigned skip_wall(4); // start at an invalid value
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t avoid;

	if (!objs.empty() && objs[objs_start].type == TYPE_TABLE) { // make sure there's enough space for the player to walk around the table
		avoid = objs[objs_start];
		avoid.expand_by_xy(get_min_front_clearance());
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
		if (is_cube_close_to_doorway(c, room, 0.0, !room.is_hallway) || interior->is_blocked_by_stairs_or_elevator(c) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_TCAN, room_id, dim, dir, 0, tot_light_amt, (cylin ? SHAPE_CYLIN : SHAPE_CUBE), tcan_colors[rgen.rand()%NUM_TCAN_COLORS]);
		return; // done
	} // for n
}

// Note: no blockers, but does check existing objects
bool building_t::add_bookcase_to_room(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement) {
	cube_t room_bounds(get_walkable_room_bounds(room));
	room_bounds.expand_by_xy(-get_trim_thickness());
	float const vspace(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return 0; // room is too small
	rand_gen_t rgen2;
	rgen2.set_state((room_id + 1), (13*mat_ix + interior->rooms.size() + 1)); // local rgen that's per-building/room; ensures bookcases are all the same size in a library
	float const width(0.4*vspace*rgen2.rand_uniform(1.0, 1.2)), depth(0.12*vspace*rgen2.rand_uniform(1.0, 1.2)), height(0.7*vspace*rgen2.rand_uniform(1.0, 1.2));
	float const clearance(max(0.2f*vspace, get_min_front_clearance()));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	set_cube_zvals(c, zval, zval+height);

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a bookcase
		bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
		if (!is_basement && classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) continue; // don't place against an exterior wall/window, inc. partial ext walls
		c.d[dim][ dir] = room_bounds.d[dim][dir]; // against this wall
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		float const pos(rgen.rand_uniform(room_bounds.d[!dim][0]+0.5*width, room_bounds.d[!dim][1]-0.5*width));
		set_wall_width(c, pos, 0.5*width, !dim);
		cube_t tc(c);
		tc.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance; // increase space to add clearance
		if (is_cube_close_to_doorway(tc, room, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(tc) || overlaps_other_room_obj(tc, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_BCASE, room_id, dim, !dir, 0, tot_light_amt); // Note: dir faces into the room, not the wall
		set_obj_id(objs);
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

bool building_t::is_room_office_bathroom(room_t const &room, float zval, unsigned floor) const {
	return room.is_office && room.get_room_type(floor) == RTYPE_BATH && !room_has_stairs_or_elevator(room, zval, floor);
}

// Note: must be first placed object
bool building_t::add_desk_to_room(rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers, colorRGBA const &chair_color,
	float zval, unsigned room_id, unsigned floor, float tot_light_amt, bool is_basement)
{
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace());
	if (min(room_bounds.dx(), room_bounds.dy()) < 1.0*vspace) return 0; // room is too small
	float const width(0.8*vspace*rgen.rand_uniform(1.0, 1.2)), depth(0.38*vspace*rgen.rand_uniform(1.0, 1.2)), height(0.21*vspace*rgen.rand_uniform(1.0, 1.2));
	float const clearance(max(0.5f*depth, get_min_front_clearance()));
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned num_placed(0);
	cube_t c, placed_desk;
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
		if (num_placed > 0 && desk_pad.intersects(placed_desk)) continue; // intersects previously placed desk
		if (!is_valid_placement_for_room(desk_pad, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
		// make short if against an exterior wall or in an office
		bool const is_tall(!room.is_office && rgen.rand_float() < 0.5 && (is_basement || classify_room_wall(room, zval, dim, dir, 0) != ROOM_WALL_EXT));
		unsigned const desk_obj_ix(objs.size());
		objs.emplace_back(c, TYPE_DESK, room_id, dim, !dir, 0, tot_light_amt, (is_tall ? SHAPE_TALL : SHAPE_CUBE));
		set_obj_id(objs);
		cube_t bc(c);
		bool const add_computer(building_obj_model_loader.is_model_valid(OBJ_MODEL_TV) && rgen.rand_bool());

		if (add_computer) {
			// add a computer monitor using the TV model
			vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TV)); // D, W, H
			float const tv_height(1.1*height), tv_hwidth(0.5*tv_height*sz.y/sz.z), tv_depth(tv_height*sz.x/sz.z), center(c.get_center_dim(!dim));
			cube_t tv;
			set_cube_zvals(tv, c.z2(), c.z2()+tv_height);
			tv.d[dim][ dir] = c. d[dim][dir] + dsign*0.25*depth; // 25% of the way from the wall
			tv.d[dim][!dir] = tv.d[dim][dir] + dsign*tv_depth;
			set_wall_width(tv, center, tv_hwidth, !dim);
			objs.emplace_back(tv, TYPE_MONITOR, room_id, dim, !dir, 0, tot_light_amt, SHAPE_SHORT, BLACK); // monitors are shorter than TVs
			set_obj_id(objs);
			// add a keyboard as well
			float const kbd_hwidth(0.7*tv_hwidth), kbd_depth(0.6*kbd_hwidth), kbd_height(0.06*kbd_hwidth);
			cube_t keyboard;
			set_cube_zvals(keyboard, c.z2(), c.z2()+kbd_height);
			keyboard.d[dim][!dir] = c.d[dim][!dir] - dsign*0.06*depth; // close to front edge
			keyboard.d[dim][ dir] = keyboard.d[dim][!dir] - dsign*kbd_depth;
			set_wall_width(keyboard, center, kbd_hwidth, !dim);
			objs.emplace_back(keyboard, TYPE_KEYBOARD, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // add as white, will be drawn with gray/black texture
			// add a computer tower under the desk
			float const cheight(0.75*height), cwidth(0.44*cheight), cdepth(0.9*cheight); // fixed AR=0.44 to match the texture
			bool const comp_side(rgen.rand_bool());
			float const pos(c.d[!dim][comp_side] + (comp_side ? -1.0 : 1.0)*0.8*cwidth);
			cube_t computer;
			set_cube_zvals(computer, c.z1(), c.z1()+cheight);
			set_wall_width(computer, pos, 0.5*cwidth, !dim);
			computer.d[dim][ dir] = c.d[dim][dir] + dsign*0.5*cdepth;
			computer.d[dim][!dir] = computer.d[dim][dir] + dsign*cdepth;
			objs.emplace_back(computer, TYPE_COMPUTER, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);
			// force even/odd-ness of obj_id based on comp_side so that we know what side to put the drawers on so that they don't intersect the computer
			if (bool(objs[desk_obj_ix].obj_id & 1) == comp_side) {++objs[desk_obj_ix].obj_id;}
		}
		else { // no computer
			if ((rgen.rand()%3) != 0) { // add sheet(s) of paper 75% of the time
				float const pheight(0.115*vspace), pwidth(0.77*pheight), thickness(0.00025*vspace); // 8.5x11

				if (pheight < 0.5*c.get_sz_dim(dim) && pwidth < 0.5*c.get_sz_dim(!dim)) { // desk is large enough for papers
					cube_t paper;
					set_cube_zvals(paper, c.z2(), c.z2()+thickness); // very thin
					unsigned const num_papers(rgen.rand() % 8); // 0-7

					for (unsigned n = 0; n < num_papers; ++n) { // okay if they overlap
						set_wall_width(paper, rgen.rand_uniform(c.d[ dim][0]+pheight, c.d[ dim][1]-pheight), 0.5*pheight,  dim);
						set_wall_width(paper, rgen.rand_uniform(c.d[!dim][0]+pwidth,  c.d[!dim][1]-pwidth),  0.5*pwidth,  !dim);
						objs.emplace_back(paper, TYPE_PAPER, room_id, dim, !dir, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt, SHAPE_CUBE, paper_colors[rgen.rand()%NUM_PAPER_COLORS]);
						set_obj_id(objs);
						paper.z2() += thickness; // to avoid Z-fighting if different colors
					} // for n
				}
			}
			float const pp_len(0.077*vspace), pp_dia(0.0028*vspace), edge_space(0.75*pp_len); // ~7.5 inches long

			if (edge_space < 0.25*min(c.dx(), c.dy())) { // desk is large enough for pens/pencils
				float const pp_z1(c.z2() + 0.3f*pp_dia); // move above papers, and avoid self shadow from the desk
				cube_t pp_bcube;
				set_cube_zvals(pp_bcube, pp_z1, pp_z1+pp_dia);
				bool const is_big_office(!is_house && room.is_office && interior->rooms.size() > 40);
				unsigned const num_pp(rgen.rand()&(is_big_office ? 2 : 3)); // 0-3 for houses, 0-2 for big office buildings

				for (unsigned n = 0; n < num_pp; ++n) {
					bool const is_pen(rgen.rand_bool());
					colorRGBA const color(is_pen ? pen_colors[rgen.rand()&3] : pencil_colors[rgen.rand()&1]);
					set_wall_width(pp_bcube, rgen.rand_uniform(c.d[ dim][0]+edge_space, c.d[ dim][1]-edge_space), 0.5*pp_len,  dim);
					set_wall_width(pp_bcube, rgen.rand_uniform(c.d[!dim][0]+edge_space, c.d[!dim][1]-edge_space), 0.5*pp_dia, !dim);
					// Note: no check for overlap with books and potted plants, but that would be complex to add and this case is rare;
					//       computer monitors/keyboards aren't added in this case, and pencils should float above papers, so we don't need to check those
					objs.emplace_back(pp_bcube, (is_pen ? TYPE_PEN : TYPE_PENCIL), room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, color);
				} // for n
			}
		}
		if (rgen.rand_float() > 0.05) { // 5% chance of no chair
			point chair_pos;
			chair_pos.z = zval;
			chair_pos[dim]  = c.d[dim][!dir];
			chair_pos[!dim] = pos + rgen.rand_uniform(-0.1, 0.1)*width; // slightly misaligned
			// there are too many desks in office buildings, and they have office chairs in cubicles anyway, so only use chair models for desks in houses with computer monitors
			bool const office_chair_model(add_computer && is_house);
			
			if (add_chair(rgen, room, blockers, room_id, chair_pos, chair_color, dim, dir, tot_light_amt, office_chair_model)) {
				cube_t const &chair(objs.back());
				if (num_placed > 0 && chair.intersects(placed_desk)) {objs.pop_back();} // intersects previously placed desk, remove it
				else {bc.union_with_cube(chair);} // include the chair
			}
		}
		++num_placed;
		if (room.is_office && num_placed == 1 && rgen.rand_float() < 0.5 && !room_has_stairs_or_elevator(room, zval, floor)) {placed_desk = bc; continue;} // allow two desks in one office
		break; // done/success
	} // for n
	return (num_placed > 0);
}

bool building_t::create_office_cubicles(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) { // assumes no prior placed objects
	if (!room.is_office) return 0; // offices only
	if (!room.interior && (rgen.rand()%3) == 0) return 0; // 66.7% chance for non-interior rooms
	cube_t const room_bounds(get_walkable_room_bounds(room));
	float const floor_spacing(get_window_vspace());
	// Note: we could choose the primary dim based on door placement like in office building bathrooms, but it seems easier to not place cubes by doors
	bool const long_dim(room.dx() < room.dy());
	float const rlength(room_bounds.get_sz_dim(long_dim)), rwidth(room_bounds.get_sz_dim(!long_dim)), midpoint(room_bounds.get_center_dim(!long_dim));
	if (rwidth < 2.5*floor_spacing || rlength < 3.5*floor_spacing) return 0; // not large enough
	unsigned const num_cubes(round_fp(rlength/(rgen.rand_uniform(0.75, 0.9)*floor_spacing))); // >= 4
	float const cube_width(rlength/num_cubes), cube_depth(cube_width*rgen.rand_uniform(0.8, 1.2)); // not quite square
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
		chair_height = 0.425*floor_spacing;
		chair_radius = 0.5f*chair_height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR);
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
				test_cube.d[!long_dim][!dir] += dir_sign*0.5*cube_depth; // allow space for people to enter the cubicle
				if (interior->is_cube_close_to_doorway(test_cube, room, 0.0, 1)) continue; // too close to a doorway; inc_open=1
				if (interior->is_blocked_by_stairs_or_elevator(test_cube)) continue;
				bool const against_window(room.d[!long_dim][dir] == part.d[!long_dim][dir]);
				objs.emplace_back(c, TYPE_CUBICLE, room_id, !long_dim, dir, 0, tot_light_amt, ((against_window && !is_middle) ? SHAPE_SHORT : SHAPE_CUBE));
				objs.back().obj_id = bldg_id;
				added_cube = 1;
				// add colliders to allow the player to enter the cubicle but not cross the side walls
				cube_t c2(c), c3(c), c4(c);
				c2.d[long_dim][0] = hi_pos - 0.06*cube_width;
				c3.d[long_dim][1] = lo_pos + 0.06*cube_width;
				c4.d[!long_dim][!dir] = wall_pos + dir_sign*0.12*cube_depth;
				objs.emplace_back(c2, TYPE_COLLIDER, room_id, !long_dim, dir, RO_FLAG_INVIS, tot_light_amt); // side1
				objs.emplace_back(c3, TYPE_COLLIDER, room_id, !long_dim, dir, RO_FLAG_INVIS, tot_light_amt); // side2
				objs.emplace_back(c4, TYPE_COLLIDER, room_id, !long_dim, dir, RO_FLAG_INVIS, tot_light_amt); // back (against wall)

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

// Note: applies to both houses and office buildings, but only houses will have bedrooms
bool building_t::can_be_bedroom_or_bathroom(room_t const &room, unsigned floor, bool skip_conn_check) const { // check room type and existence of exterior door
	if (room.has_stairs_on_floor(floor) || room.has_elevator || room.is_hallway || room.is_office || room.is_sec_bldg) return 0; // no bed/bath in these cases
	
	if (floor == 0) { // run special logic for bedrooms and bathrooms (private rooms) on the first floor of a house
		if (is_room_adjacent_to_ext_door(room)) return 0; // door to house does not open into a bedroom/bathroom
		if (skip_conn_check) return 1;

		if (room.dz() > 1.5*get_window_vspace()) { // more than one floor
			if (interior->stairwells.empty()) return 1; // failed to place stairs in this house, maybe because it was too small; I guess we just return 0 here
			// determine if this room is on the shortest path from an exterior door to the stairs; if so, it can't be a bedroom or bathroom;
			// okay, that's not easy/fast to do, so determine if there is any path from the exterior door to the stairs that doesn't go through this room;
			// this won't work when there are two paths from the door to the stairs and this room is only on one of the paths, so we could put a BR/BR on both paths
			int cur_room(-1);
			vector<unsigned> door_rooms, stairs_rooms;

			for (unsigned i = 0; i < interior->rooms.size(); ++i) {
				room_t const &r(interior->rooms[i]);
				if (r == room) {cur_room = i; continue;} // this room; we know it can't have stairs or an exterior door
				if (r.is_sec_bldg || r.z2() <= ground_floor_z1) continue; // skip basement rooms, garages, and sheds
				if (r.has_stairs_on_floor(floor))  {stairs_rooms.push_back(i);}
				if (is_room_adjacent_to_ext_door(r)) {door_rooms.push_back(i);}
			}
			if (cur_room < 0 || stairs_rooms.empty()) {cout << TXT(bcube.str());}
			assert(cur_room >= 0); // must be found
			assert(!stairs_rooms.empty());

			if (!is_rotated() && is_simple_cube() && !has_complex_floorplan) { // too strong for rotated or non-cube buildings, where door placement can sometimes fail
				assert(!doors.empty());
				assert(!door_rooms.empty());
			}
			for (auto d = door_rooms.begin(); d != door_rooms.end(); ++d) {
				for (auto s = stairs_rooms.begin(); s != stairs_rooms.end(); ++s) {
					if (!are_rooms_connected_without_using_room(*d, *s, cur_room)) return 0;
				}
			}
		}
	}
	return 1;
}
bool building_t::can_be_bathroom(room_t const &room) const { // Note: assumes caller has checked can_be_bedroom_or_bathroom()
	float const vspace(get_window_vspace());
	return (min(room.dx(), room.dy()) < 2.0*vspace && max(room.dx(), room.dy()) < 3.0*vspace && count_num_int_doors(room) == 1);
}

unsigned building_t::count_num_int_doors(room_t const &room) const {
	cube_t room_exp(room);
	float const wall_thickness(get_wall_thickness());
	room_exp.expand_by(wall_thickness, wall_thickness, -wall_thickness); // expand in XY and shrink in Z
	unsigned num(0);
	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {num += i->intersects(room_exp);}
	return num;
}

bool building_t::check_valid_closet_placement(cube_t const &c, room_t const &room, unsigned objs_start, unsigned bed_ix, float min_bed_space) const {
	if (min_bed_space > 0.0) {
		room_object_t const &bed(interior->room_geom->get_room_object_by_index(bed_ix));
		assert(bed.type == TYPE_BED);
		cube_t bed_exp(bed);
		bed_exp.expand_by_xy(min_bed_space);
		if (c.intersects_xy(bed_exp)) return 0; // too close to bed
	}
	return (!overlaps_other_room_obj(c, objs_start) && !is_cube_close_to_doorway(c, room, 0.0, 1));
}

float get_lamp_width_scale() {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_LAMP)); // L, W, H
	return ((sz == zero_vector) ? 0.0 : 0.5f*(sz.x + sz.y)/sz.z);
}

bool building_t::add_bedroom_objs(rand_gen_t rgen, room_t const &room, vect_cube_t const &blockers, float zval, unsigned room_id,
	unsigned floor, float tot_light_amt, unsigned objs_start, bool room_is_lit, bool is_basement, light_ix_assign_t &light_ix_assign)
{
	if (room.interior) return 0; // bedrooms should have at least one window; if windowless/interior, it can't be a bedroom
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const bed_obj_ix(objs.size()); // if placed, it will be this index
	if (!add_bed_to_room(rgen, room, blockers, zval, room_id, tot_light_amt, floor)) return 0; // it's only a bedroom if there's bed
	assert(bed_obj_ix < objs.size());
	room_object_t const bed(objs[bed_obj_ix]); // deep copy so that we don't need to worry about invalidating the reference below
	float const window_vspacing(get_window_vspace());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-get_trim_thickness()); // shrink to leave a small gap
	// closet
	float const doorway_width(get_doorway_width()), floor_thickness(get_floor_thickness()), front_clearance(max(0.6f*doorway_width, get_min_front_clearance()));
	float const closet_min_depth(0.65*doorway_width), closet_min_width(1.5*doorway_width), min_dist_to_wall(1.0*doorway_width), min_bed_space(front_clearance);
	unsigned const first_corner(rgen.rand() & 3);
	bool const first_dim(rgen.rand_bool());
	cube_t const part(get_part_for_room(room));
	bool placed_closet(0);
	unsigned closet_obj_id(0);
	bool chk_windows[2][2] = {0}; // precompute which walls are exterior and can have windows, {dim}x{dir}

	if (!is_basement && has_windows()) { // are bedrooms ever plaed in the basement?
		for (unsigned d = 0; d < 4; ++d) {chk_windows[d>>1][d&1] = (classify_room_wall(room, zval, (d>>1), (d&1), 0) == ROOM_WALL_EXT);}
	}
	for (unsigned n = 0; n < 4 && !placed_closet; ++n) { // try 4 room corners
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
			if (chk_windows[!dim][other_dir] && is_val_inside_window(part, dim, c.d[dim][!dir], window_hspacing, get_window_h_border())) continue; // check for window intersection
			c.z2() += window_vspacing - floor_thickness;
			c.d[dim][!dir] += signed_front_clearance; // extra padding in front, to avoid placing too close to bed
			if (!check_valid_closet_placement(c, room, objs_start, bed_obj_ix, min_bed_space)) continue; // bad placement
			// good placement, see if we can make the closet larger
			unsigned const num_steps = 10;
			float const req_dist(chk_windows[!dim][!other_dir] ? (other_dir ? -1.0 : 1.0)*min_dist_to_wall : 0.0); // signed; at least min dist from the opposite wall if it's exterior
			float const max_grow((room_bounds.d[!dim][!other_dir] - req_dist) - c.d[!dim][!other_dir]);
			float const len_step(max_grow/num_steps), depth_step(dir_sign*0.35*doorway_width/num_steps); // signed

			for (unsigned s1 = 0; s1 < num_steps; ++s1) { // try increasing width
				cube_t c2(c);
				c2.d[!dim][!other_dir] += len_step;
				if (!check_valid_closet_placement(c2, room, objs_start, bed_obj_ix, min_bed_space)) break; // bad placement
				c = c2; // valid placement, update with larger cube
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
			if (c.d[!dim][0] == room_bounds.d[!dim][0]) {flags |= RO_FLAG_ADJ_LO;}
			if (c.d[!dim][1] == room_bounds.d[!dim][1]) {flags |= RO_FLAG_ADJ_HI;}
			//if ((rgen.rand() % 10) == 0) {flags |= RO_FLAG_OPEN;} // 10% chance of open closet; unclear if this adds any value, but it works
			closet_obj_id = objs.size();
			objs.emplace_back(c, TYPE_CLOSET, room_id, dim, !dir, flags, tot_light_amt, SHAPE_CUBE, wall_color); // closet door is always white; sides should match interior walls
			set_obj_id(objs);
			if (flags & RO_FLAG_OPEN) {interior->room_geom->expand_object(objs.back());} // expand opened closets
			placed_closet = 1; // done
			// add a light inside the closet
			room_object_t const &closet(objs.back());
			point lpos(closet.xc(), closet.yc(), closet.z2());
			lpos[dim] += 0.05*c.get_sz_dim(dim)*(dir ? -1.0 : 1.0); // move slightly toward the front of the closet
			cube_t light(lpos);
			light.z1() -= 0.02*window_vspacing;
			light.expand_by_xy((closet.is_small_closet() ? 0.04 : 0.06)*window_vspacing);
			colorRGBA const color(1.0, 1.0, 0.9); // yellow-ish
			objs.emplace_back(light, TYPE_LIGHT, room_id, dim, 0, (RO_FLAG_NOCOLL | RO_FLAG_IN_CLOSET), 0.0, SHAPE_CYLIN, color); // dir=0 (unused)
			objs.back().obj_id = light_ix_assign.get_next_ix();

			if (closet.is_small_closet()) { // add a blocker in front of the closet to avoid placing furniture that blocks the door from opening
				c.d[dim][!dir] += dir_sign*doorway_width;
				objs.emplace_back(c, TYPE_BLOCKER, room_id, dim, 0, RO_FLAG_INVIS);
			}
		} // for d
	} // for n
	// dresser
	float const ds_height(rgen.rand_uniform(0.26, 0.32)*window_vspacing), ds_depth(rgen.rand_uniform(0.20, 0.25)*window_vspacing), ds_width(rgen.rand_uniform(0.6, 0.9)*window_vspacing);
	vector3d const ds_sz_scale(ds_depth/ds_height, ds_width/ds_height, 1.0);
	place_obj_along_wall(TYPE_DRESSER, room, ds_height, ds_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0);
	// nightstand
	unsigned const pref_orient(2*bed.dim + (!bed.dir)); // prefer the same orient as the bed so that it's placed on the same wall next to the bed
	float const ns_height(rgen.rand_uniform(0.24, 0.26)*window_vspacing), ns_depth(rgen.rand_uniform(0.15, 0.2)*window_vspacing), ns_width(rgen.rand_uniform(1.0, 2.0)*ns_depth);
	vector3d const ns_sz_scale(ns_depth/ns_height, ns_width/ns_height, 1.0);
	place_obj_along_wall(TYPE_NIGHTSTAND, room, ns_height, ns_sz_scale, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0, pref_orient);

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

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) { // choose the dresser or nightstand closest to be bed
			if (i->type != TYPE_DRESSER && i->type != TYPE_NIGHTSTAND) continue; // not a dresser or nightstand
			float const dist_sq(p2p_dist_xy_sq(i->get_cube_center(), pillow_center));
			if (dmin_sq == 0.0 || dist_sq < dmin_sq) {obj_id = (i - objs.begin()); dmin_sq = dist_sq;}
		}
		if (obj_id >= 0) { // found a valid object to place this on
			room_object_t const &obj(objs[obj_id]);
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
			objs.emplace_back(lamp, TYPE_LAMP, room_id, obj.dim, obj.dir, flags, tot_light_amt, SHAPE_CYLIN, lamp_colors[rgen.rand()%NUM_LAMP_COLORS]); // Note: invalidates obj ref
		}
	}
	if (rgen.rand_float() < global_building_params.ball_prob) { // maybe add a ball to the room
		float const radius(0.048*window_vspacing); // 4.7 inches
		cube_t ball_area(place_area);
		ball_area.expand_by_xy(-radius*rgen.rand_uniform(1.0, 10.0));

		if (ball_area.is_strictly_normalized()) {
			for (unsigned n = 0; n < 10; ++n) { // make 10 attempts to place the object
				bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
				point center(0.0, 0.0, (zval + radius));
				center[ dim] = ball_area.d[dim][dir];
				center[!dim] = rgen.rand_uniform(ball_area.d[!dim][0], ball_area.d[!dim][1]); // random position along the wall
				cube_t c(center);
				c.expand_by(radius);
				if (overlaps_other_room_obj(c, objs_start) || interior->is_blocked_by_stairs_or_elevator(c) || is_cube_close_to_doorway(c, room, 0.0, 1)) continue; // bad placement
				objs.emplace_back(c, TYPE_LG_BALL, room_id, 0, 0, RO_FLAG_DSTATE, tot_light_amt, SHAPE_SPHERE, WHITE);
				objs.back().obj_id     = (uint16_t)interior->room_geom->allocate_dynamic_state(); // allocate a new dynamic state object
				objs.back().item_flags = rgen.rand_bool(); // selects ball type
				break; // done
			} // for n
		}
	}
	return 1; // success
}

// Note: must be first placed object
bool building_t::add_bed_to_room(rand_gen_t &rgen, room_t const &room, vect_cube_t const &blockers, float zval, unsigned room_id, float tot_light_amt, unsigned floor) {
	unsigned const NUM_COLORS = 8;
	colorRGBA const colors[NUM_COLORS] = {WHITE, WHITE, WHITE, LT_BLUE, LT_BLUE, PINK, PINK, LT_GREEN}; // color of the sheets
	cube_t room_bounds(get_walkable_room_bounds(room));
	float const vspace(get_window_vspace()), wall_thick(get_wall_thickness());
	bool const dim(room_bounds.dx() < room_bounds.dy()); // longer dim
	vector3d expand, bed_sz;
	expand[ dim] = -wall_thick; // small amount of space
	expand[!dim] = -0.3f*vspace; // leave at least some space between the bed and the wall
	room_bounds.expand_by_xy(expand);
	float const room_len(room_bounds.get_sz_dim(dim)), room_width(room_bounds.get_sz_dim(!dim));
	
	if (floor == 0) { // special case for ground floor
		if (room_len < 1.3*vspace || room_width < 0.7*vspace) return 0; // room is too small to fit a bed
		if (room_len > 4.0*vspace || room_width > 2.5*vspace) return 0; // room is too large to be a bedroom
	}
	else { // more relaxed constraints
		if (room_len < 1.1*vspace || room_width < 0.6*vspace) return 0; // room is too small to fit a bed
		if (room_len > 4.5*vspace || room_width > 3.5*vspace) return 0; // room is too large to be a bedroom
	}
	bool const first_head_dir(rgen.rand_bool()), first_wall_dir(rgen.rand_bool());
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	c.z1() = zval;

	for (unsigned n = 0; n < 20; ++n) { // make 20 attempts to place a bed
		float const sizes[6][2] = {{38, 75}, {38, 80}, {53, 75}, {60, 80}, {76, 80}, {72, 84}}; // twin, twin XL, full, queen, king, cal king
		unsigned const size_ix((room_width < 0.9*vspace) ? (rgen.rand() % 6) : (2 + (rgen.rand() % 4))); // only add twin beds to narrow rooms
		bed_sz[ dim] = 0.01f*vspace*(sizes[size_ix][1] + 8.0f); // length (mattress + headboard + footboard)
		bed_sz[!dim] = 0.01f*vspace*(sizes[size_ix][0] + 4.0f); // width  (mattress + small gaps)
		if (room_bounds.dx() < 1.5*bed_sz.x || room_bounds.dy() < 1.5*bed_sz.y) continue; // room is too small for a bed of this size
		bed_sz.z = 0.3*vspace*rgen.rand_uniform(1.0, 1.2); // height
		c.z2()   = zval + bed_sz.z;

		for (unsigned d = 0; d < 2; ++d) {
			float const min_val(room_bounds.d[d][0]), max_val(room_bounds.d[d][1] - bed_sz[d]);

			if (bool(d) == dim && n < 5) { // in the first few iterations, try to place the head of the bed against the wall (maybe not for exterior wall facing window?)
				c.d[d][0] = ((first_head_dir ^ bool(n&1)) ? min_val : max_val);
			}
			else if (bool(d) != dim && rgen.rand_bool()) { // try to place the bed against the wall sometimes
				c.d[d][0] = ((first_wall_dir ^ bool(n&1)) ? (min_val - 0.25*vspace) : (max_val + 0.25*vspace));
			}
			else {
				c.d[d][0] = rgen.rand_uniform(min_val, max_val);
			}
			c.d[d][1] = c.d[d][0] + bed_sz[d];
		} // for d
		if (!is_valid_placement_for_room(c, room, blockers, 1)) continue; // check proximity to doors and collision with blockers
		bool const dir((room_bounds.d[dim][1] - c.d[dim][1]) < (c.d[dim][0] - room_bounds.d[dim][0])); // head of the bed is closer to the wall
		objs.emplace_back(c, TYPE_BED, room_id, dim, dir, 0, tot_light_amt);
		room_object_t &bed(objs.back());
		bed.obj_id = (uint16_t)objs.size();
		// use white color if a texture is assigned that's not close to white
		int const sheet_tid(bed.get_sheet_tid());
		if (sheet_tid < 0 || sheet_tid == WHITE_TEX || texture_color(sheet_tid).get_luminance() > 0.5) {bed.color = colors[rgen.rand()%NUM_COLORS];}
		return 1; // done/success
	} // for n
	return 0;
}

// Note: modified blockers rather than using it; fireplace must be the first placed object
bool building_t::maybe_add_fireplace_to_room(room_t const &room, vect_cube_t &blockers, float zval, unsigned room_id, float tot_light_amt) {
	if (has_int_fplace) return 0; // already added an interior fireplace
	// Note: the first part of the code below is run on every first floor room and will duplicate work, so it may be better to factor it out somehow
	cube_t fireplace(get_fireplace()); // make a copy of the exterior fireplace that will be converted to an interior fireplace
	bool dim(0), dir(0);
	if      (fireplace.x1() <= bcube.x1()) {dim = 0; dir = 0;} // Note: may not work on rotated buildings
	else if (fireplace.x2() >= bcube.x2()) {dim = 0; dir = 1;}
	else if (fireplace.y1() <= bcube.y1()) {dim = 1; dir = 0;}
	else if (fireplace.y2() >= bcube.y2()) {dim = 1; dir = 1;}
	else {assert(is_rotated()); return 0;} // can fail on rotated buildings?
	float const depth_signed((dir ? -1.0 : 1.0)*1.0*fireplace.get_sz_dim(dim)), wall_pos(fireplace.d[dim][!dir]);
	fireplace.d[dim][ dir] = wall_pos; // flush with the house wall
	fireplace.d[dim][!dir] = wall_pos + depth_signed; // extend out into the room
	fireplace.z2() -= 0.15*fireplace.dz(); // shorten slightly
	cube_t room_exp(room);
	room_exp.expand_by_xy(0.5*get_wall_thickness()); // allow fireplace to extend slightly into room walls
	if (!room_exp.contains_cube_xy(fireplace)) return 0; // fireplace not in this room
	// the code below should be run at most once per building
	cube_t fireplace_ext(fireplace);
	fireplace_ext.d[dim][!dir] = fireplace.d[dim][!dir] + 0.5*depth_signed; // extend out into the room even further for clearance
	if (interior->is_blocked_by_stairs_or_elevator(fireplace_ext)) return 0; // blocked by stairs, don't add (would be more correct to relocate stairs) - should no longer fail
	fireplace.d[dim][dir] = room.d[dim][dir]; // re-align to room to remove any gap between the fireplace and the exterior wall
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(fireplace, TYPE_FPLACE, room_id, dim, dir, 0, tot_light_amt);
	cube_t blocker(fireplace_ext);
	blocker.d[dim][ dir] = fireplace.d[dim][!dir]; // flush with the front of the fireplace
	objs.emplace_back(blocker, TYPE_BLOCKER, room_id, dim, dir, RO_FLAG_INVIS);
	blockers.push_back(fireplace_ext); // add as a blocker if it's not already there
	has_int_fplace = 1;
	return 1;
}

bool building_t::place_obj_along_wall(room_object type, room_t const &room, float height, vector3d const &sz_scale, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
	cube_t const &place_area, unsigned objs_start, float front_clearance, unsigned pref_orient, bool pref_centered, colorRGBA const &color, bool not_at_window, room_obj_shape shape)
{
	float const hwidth(0.5*height*sz_scale.y/sz_scale.z), depth(height*sz_scale.x/sz_scale.z), min_space(2.8*hwidth);
	vector3d const place_area_sz(place_area.get_size());
	if (max(place_area_sz.x, place_area_sz.y) <= min_space) return 0; // can't fit in either dim
	unsigned const force_dim((place_area_sz.x <= min_space) ? 0 : ((place_area_sz.y <= min_space) ? 1 : 2)); // *other* dim; 2=neither
	float const obj_clearance(depth*front_clearance), clearance(max(obj_clearance, get_min_front_clearance()));
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c;
	set_cube_zvals(c, zval, zval+height);
	bool center_tried[4] = {};

	for (unsigned n = 0; n < 25; ++n) { // make 25 attempts to place the object
		bool const use_pref(pref_orient < 4 && n < 10); // use pref orient for first 10 tries
		bool const dim((force_dim < 2) ? force_dim : (use_pref ? (pref_orient >> 1) : rgen.rand_bool())); // choose a random wall unless forced
		bool const dir(use_pref ? !(pref_orient & 1) : rgen.rand_bool()); // dir is inverted for the model, so we invert pref dir as well
		unsigned const orient(2*dim + dir);
		float center(0.0);
		if (pref_centered && !center_tried[orient]) {center = place_area.get_center_dim(!dim); center_tried[orient] = 1;} // try centered
		else {center = rgen.rand_uniform(place_area.d[!dim][0]+hwidth, place_area.d[!dim][1]-hwidth);} // random position
		c.d[ dim][ dir] = place_area.d[dim][dir];
		c.d[ dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*depth;
		c.d[!dim][   0] = center - hwidth;
		c.d[!dim][   1] = center + hwidth;

		if (not_at_window && classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) {
			cube_t const part(get_part_for_room(room));
			float const hspacing(get_hspacing_for_part(part, !dim)), border(get_window_h_border());
			// assume object is no larger than 2x window size and check left, right, and center positions
			if (is_val_inside_window(part, !dim, c.d[!dim][0], hspacing, border) ||
				is_val_inside_window(part, !dim, c.d[!dim][1], hspacing, border) ||
				is_val_inside_window(part, !dim, c.get_center_dim(!dim), hspacing, border)) continue;
		}
		cube_t c2(c); // used for collision tests
		c2.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance;
		if (overlaps_other_room_obj(c2, objs_start) || interior->is_blocked_by_stairs_or_elevator(c2)) continue; // bad placement
		// we don't need clearance for both the door and the object; test the object itself against the open door and the object with clearance against the closed door
		if (is_cube_close_to_doorway(c, room, 0.0, 1)) continue; // bad placement
		cube_t c3(c); // used for collision tests
		c3.d[dim][!dir] += (dir ? -1.0 : 1.0)*obj_clearance; // smaller clearance value (without player diameter)
		if (is_cube_close_to_doorway(c3, room, 0.0, 0)) continue; // bad placement
		unsigned const flags((type == TYPE_BOX) ? (RO_FLAG_ADJ_LO << orient) : 0); // set wall edge bit for boxes (what about other dim bit if place in room corner?)
		objs.emplace_back(c, type, room_id, dim, !dir, flags, tot_light_amt, shape, color);
		set_obj_id(objs);
		if (front_clearance > 0.0) {objs.emplace_back(c2, TYPE_BLOCKER, room_id, dim, !dir, RO_FLAG_INVIS);} // add blocker cube to ensure no other object overlaps this space
		return 1; // done
	} // for n
	return 0; // failed
}
bool building_t::place_model_along_wall(unsigned model_id, room_object type, room_t const &room, float height, rand_gen_t &rgen, float zval, unsigned room_id, float tot_light_amt,
	cube_t const &place_area, unsigned objs_start, float front_clearance, unsigned pref_orient, bool pref_centered, colorRGBA const &color, bool not_at_window)
{
	if (!building_obj_model_loader.is_model_valid(model_id)) return 0; // don't have a model of this type
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(model_id)); // D, W, H
	return place_obj_along_wall(type, room, height*get_window_vspace(), sz, rgen, zval, room_id, tot_light_amt,
		place_area, objs_start, front_clearance, pref_orient, pref_centered, color, not_at_window);
}

float building_t::add_flooring(room_t const &room, float &zval, unsigned room_id, float tot_light_amt) {
	float const new_zval(zval + 0.0012*get_window_vspace());
	cube_t floor(get_walkable_room_bounds(room));
	set_cube_zvals(floor, zval, new_zval);
	interior->room_geom->objs.emplace_back(floor, TYPE_FLOORING, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt);
	return new_zval;
}

bool building_t::add_bathroom_objs(rand_gen_t rgen, room_t const &room, float &zval, unsigned room_id, float tot_light_amt,
	unsigned objs_start, unsigned floor, bool is_basement, unsigned &added_bathroom_objs_mask)
{
	// Note: zval passed by reference
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.5*wall_thickness);
	if (min(place_area.dx(), place_area.dy()) < 0.7*floor_spacing) return 0; // room is too small (should be rare)
	bool const have_toilet(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOILET)), have_sink(building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK));
	vect_room_object_t &objs(interior->room_geom->objs);

	if (!is_house && (have_toilet || have_sink)) { // office with at least a toilet or sink - replace carpet with tile
		zval       = add_flooring(room, zval, room_id, tot_light_amt); // move the effective floor up
		objs_start = objs.size(); // exclude this from collision checks
	}
	if (have_toilet && room.is_office && min(place_area.dx(), place_area.dy()) > 1.5*floor_spacing && max(place_area.dx(), place_area.dy()) > 2.0*floor_spacing) {
		if (divide_bathroom_into_stalls(rgen, room, zval, room_id, tot_light_amt, floor)) { // large enough, try to divide into bathroom stalls
			added_bathroom_objs_mask |= (PLACED_TOILET | PLACED_SINK);
			return 1;
		}
	}
	bool placed_obj(0), placed_toilet(0);
	
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
				objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker cube to ensure no other object overlaps this space
				placed_obj = placed_toilet = 1; // done
				added_bathroom_objs_mask  |= PLACED_TOILET;

				// try to place a roll of toilet paper on the adjacent wall
				bool const tp_dir(dim ? xdir : ydir);
				float const length(0.18*height), wall_pos(c.get_center_dim(dim)), far_edge_pos(wall_pos + (dir ? -1.0 : 1.0)*0.5*length);
				cube_t const part(get_part_for_room(room));

				// if this wall has windows and bathroom has multiple exterior walls (which means it has non-glass block windows), don't place a TP roll
				if (is_basement || !has_windows() || classify_room_wall(room, zval, !dim, tp_dir, 0) != ROOM_WALL_EXT ||
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
				if (is_basement || !has_windows() || classify_room_wall(room, zval, toilet.dim, !toilet.dir, 0) != ROOM_WALL_EXT) { // check for possible windows
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
	if (is_house && !is_basement && (floor > 0 || rgen.rand_bool())) { // try to add a shower; 50% chance if on first floor; not in basements (due to drawing artifacts)
		float const shower_height(0.8*floor_spacing);
		float shower_dx(rgen.rand_uniform(0.4, 0.5)*floor_spacing), shower_dy(rgen.rand_uniform(0.4, 0.5)*floor_spacing);
		bool hdim(shower_dx < shower_dy); // larger dim, ust match handle/door drawing code
		unsigned const first_corner(rgen.rand() & 3);
		//cube_t const part(get_part_for_room(room));
		bool placed_shower(0), is_ext_wall[2][2] = {0};
		
		if (!is_basement && has_windows()) { // precompute which walls are exterior, {dim}x{dir}; basement walls are not considered exterior because there are no windows
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
				bool is_bad(0);

				for (unsigned d = 0; d < 2; ++d) { // check for window intersection
					// Update: exterior walls aren't drawn in the correct order for glass alpha blend, so skip any exterior walls
					if (is_ext_wall[!d][dirs[!d]] /*&& is_val_inside_window(part, d, c.d[d][!dirs[d]], get_hspacing_for_part(part, d), get_window_h_border())*/) {is_bad = 1; break;}
				}
				if (is_bad) continue;
				cube_t c2(c); // used for placement tests; extend out by door width on the side that opens, and a small amount on the other side
				c2.d[0][!xdir] += (xdir ? -1.0 : 1.0)*((!hdim) ? 1.1*shower_dy : 0.2*shower_dx);
				c2.d[1][!ydir] += (ydir ? -1.0 : 1.0)*(  hdim  ? 1.1*shower_dx : 0.2*shower_dy);
				if (overlaps_other_room_obj(c2, objs_start) || is_cube_close_to_doorway(c2, room, 0.0, 1)) continue; // bad placement
				objs.emplace_back(c,  TYPE_SHOWER,  room_id, xdir, ydir, 0, tot_light_amt);
				objs.emplace_back(c2, TYPE_BLOCKER, room_id, 0, 0, RO_FLAG_INVIS); // add blocker cube to ensure no other object overlaps this space
				placed_obj = placed_shower = 1;
				added_bathroom_objs_mask  |= PLACED_SHOWER;
				break; // done
			} // for n
			if (placed_shower) break; // done
			swap(shower_dx, shower_dy); // try the other aspect ratio
			hdim ^= 1;
		} // for ar
	}
	if (is_house && (!is_basement || rgen.rand_bool())) { // 50% of the time if in the basement
		// place a tub, but not in office buildings; placed before the sink because it's the largest and the most limited in valid locations
		cube_t place_area_tub(room_bounds);
		place_area_tub.expand_by(-get_trim_thickness()); // just enough to prevent z-fighting and intersecting the wall trim
		
		if (place_model_along_wall(OBJ_MODEL_TUB, TYPE_TUB, room, 0.2, rgen, zval, room_id, tot_light_amt, place_area_tub, objs_start, 0.4)) {
			placed_obj = 1;
			added_bathroom_objs_mask |= PLACED_TUB;
		}
	}
	if (place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6)) {
		placed_obj = 1;
		added_bathroom_objs_mask |= PLACED_SINK;
		room_object_t const &sink((objs.back().type == TYPE_SINK) ? objs.back() : objs[objs.size()-2]); // find sink, skip blocker
		
		if (is_basement || classify_room_wall(room, zval, sink.dim, !sink.dir, 0) != ROOM_WALL_EXT) { // interior wall only
			// add a mirror/medicine cabinet above the sink; could later make into medicine cabinet
			cube_t mirror(sink); // start with the sink left and right position
			mirror.expand_in_dim(!sink.dim, 0.1*mirror.get_sz_dim(!sink.dim)); // make slightly wider
			set_cube_zvals(mirror, sink.z2(), sink.z2()+0.3*floor_spacing);
			mirror.d[sink.dim][!sink.dir] = room_bounds.d[sink.dim][!sink.dir];
			mirror.d[sink.dim][ sink.dir] = mirror.d[sink.dim][!sink.dir] + (sink.dir ? 1.0 : -1.0)*1.0*wall_thickness; // thickness
			// this mirror is actually 3D, so we enable collision detection; treat as a house even if it's in an office building
			objs.emplace_back(mirror, TYPE_MIRROR, room_id, sink.dim, sink.dir, RO_FLAG_IS_HOUSE, tot_light_amt);
		}
	}
	return placed_obj;
}

bool building_t::add_tp_roll(cube_t const &room, unsigned room_id, float tot_light_amt, bool dim, bool dir, float length, float zval, float wall_pos, bool check_valid_pos) {
	float const diameter(length);
	cube_t tp;
	set_cube_zvals(tp, zval, (zval + diameter));
	set_wall_width(tp, wall_pos, 0.5*length, !dim); // set length
	tp.d[dim][ dir] = room.d[dim][dir]; // against the wall
	tp.d[dim][!dir] = tp  .d[dim][dir] + (dir ? -1.0 : 1.0)*diameter; // set the diameter
	// Note: not checked against other bathroom objects because the toilet is placed first
	if (check_valid_pos && (!room.contains_cube(tp) || is_cube_close_to_doorway(tp, room, 0.0, 1))) return 0;
	interior->room_geom->objs.emplace_back(tp, TYPE_TPROLL, room_id, dim, dir, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN, WHITE);
	set_obj_id(interior->room_geom->objs);
	return 1;
}

bool building_t::divide_bathroom_into_stalls(rand_gen_t &rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned floor) {
	// Note: assumes no prior placed objects
	bool const use_sink_model(0 && building_obj_model_loader.is_model_valid(OBJ_MODEL_SINK)); // not using sink models
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
	bool br_dim(room.dy() < room.dx()), sink_side(0), sink_side_set(0); // br_dim is the smaller dim
	cube_t place_area(room), br_door;
	place_area.expand_by(-0.5*wall_thickness);

	// determine men's room vs. women's room
	point const part_center(get_part_for_room(room).get_cube_center()), room_center(room.get_cube_center());
	bool mens_room((part_center.x < room_center.x) ^ (part_center.y < room_center.y)), has_second_bathroom(0);

	// if there are two bathrooms (one on each side of the building), assign a gender to each side; if only one, alternate gender per floor
	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (r->part_id != room.part_id || &(*r) == &room) continue; // different part or same room
		if (is_room_office_bathroom(*r, zval, floor)) {has_second_bathroom = 1; break;}
	}
	if (!has_second_bathroom) {mens_room ^= (floor & 1);}
	bool const add_urinals(mens_room && building_obj_model_loader.is_model_valid(OBJ_MODEL_URINAL));

	if (add_urinals) { // use urinal model
		vector3d const usz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_URINAL)); // L, W, H
		uheight = 0.4*floor_spacing; uwidth = uheight*usz.y/usz.z; ulength = uheight*usz.x/usz.z;
	}
	for (unsigned d = 0; d < 2 && !sink_side_set; ++d) {
		for (unsigned side = 0; side < 2 && !sink_side_set; ++side) {
			cube_t c(room);
			set_cube_zvals(c, zval, zval+wall_thickness); // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
			c.d[!br_dim][!side] = c.d[!br_dim][side] + (side ? -1.0 : 1.0)*wall_thickness; // shrink to near zero area in this dim

			for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
				if ((i->dy() < i->dx()) == br_dim) continue; // door in wrong dim
				if (!is_cube_close_to_door(c, 0.0, 0, *i, 2)) continue; // check both dirs
				sink_side = side; sink_side_set = 1;
				place_area.d[!br_dim][side] += (sink_side ? -1.0 : 1.0)*(i->get_sz_dim(br_dim) - 0.25*swidth); // add sink clearance for the door to close
				br_door = *i;
				break; // sinks are on the side closest to the door
			}
		} // for side
		if (d == 0 && !sink_side_set) {br_dim ^= 1;} // door not found on long dim - R90 and try short dim
	} // for d
	assert(sink_side_set);
	float const room_len(place_area.get_sz_dim(!br_dim)), room_width(place_area.get_sz_dim(br_dim));
	float const sinks_len(0.4*room_len), stalls_len(room_len - sinks_len), req_depth(2.0f*max(stall_depth, slength));
	if (room_width < req_depth) return 0;
	unsigned const num_stalls(std::floor(stalls_len/stall_width)), num_sinks(std::floor(sinks_len/sink_spacing));
	if (num_stalls < 2 || num_sinks < 2) return 0; // not enough space for 2 stalls and 2 sinks
	stall_width  = stalls_len/num_stalls; // reclaculate to fill the gaps
	sink_spacing = sinks_len/num_sinks;
	bool const two_rows(room_width > 1.5*req_depth), skip_stalls_side(room_id & 1); // put stalls on a side consistent across floors
	float const sink_side_sign(sink_side ? 1.0 : -1.0), stall_step(sink_side_sign*stall_width), sink_step(-sink_side_sign*sink_spacing);
	float const floor_thickness(get_floor_thickness());
	unsigned const NUM_STALL_COLORS = 4;
	colorRGBA const stall_colors[NUM_STALL_COLORS] = {colorRGBA(0.75, 1.0, 0.9, 1.0), colorRGBA(0.7, 0.8, 1.0), WHITE, DK_GRAY}; // blue-green, light blue
	colorRGBA const stall_color(stall_colors[interior->doors.size() % NUM_STALL_COLORS]); // random, but constant for each building
	vect_room_object_t &objs(interior->room_geom->objs);

	for (unsigned dir = 0; dir < 2; ++dir) { // each side of the wall
		if (!two_rows && dir == (unsigned)skip_stalls_side) continue; // no stalls/sinks on this side
		// add stalls
		float const dir_sign(dir ? -1.0 : 1.0), wall_pos(place_area.d[br_dim][dir]), stall_from_wall(wall_pos + dir_sign*(0.5*tlength + wall_thickness));
		float stall_pos(place_area.d[!br_dim][!sink_side] + 0.5*stall_step);

		for (unsigned n = 0; n < num_stalls; ++n) {
			point center(stall_from_wall, stall_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t toilet(center, center), stall(toilet);
			toilet.expand_in_dim( br_dim, 0.5*tlength);
			toilet.expand_in_dim(!br_dim, 0.5*twidth);
			toilet.z2() += theight;
			stall.z2() = stall.z1() + floor_spacing - floor_thickness; // set stall height to room height
			stall.expand_in_dim(!br_dim, 0.5*stall_width);
			stall.d[br_dim][ dir] = wall_pos; // + wall_thickness?
			stall.d[br_dim][!dir] = wall_pos + dir_sign*stall_depth;
			
			if (!interior->is_cube_close_to_doorway(stall, room, 0.0, 1)) { // skip if close to a door (for rooms with doors at both ends); inc_open=1
				bool const is_open(rgen.rand_bool()); // 50% chance of stall door being open
				objs.emplace_back(toilet, TYPE_TOILET, room_id, br_dim, !dir, 0, tot_light_amt);
				objs.emplace_back(stall,  TYPE_STALL,  room_id, br_dim,  dir, (is_open ? RO_FLAG_OPEN : 0), tot_light_amt, SHAPE_CUBE, stall_color);
				float const tp_length(0.18*theight), wall_pos(toilet.get_center_dim(br_dim));
				cube_t stall_inner(stall);
				stall_inner.expand_in_dim(!br_dim, -0.0125*stall.dz()); // subtract off stall wall thickness
				add_tp_roll(stall_inner, room_id, tot_light_amt, !br_dim, dir, tp_length, (zval + 0.7*theight), wall_pos);
			}
			stall_pos += stall_step;
		} // for n
		if (add_urinals && dir == (unsigned)skip_stalls_side) continue; // no urinals and sinks are each on one side
		// add sinks
		float const sink_start(place_area.d[!br_dim][sink_side] + 0.5f*sink_step);
		float const sink_from_wall(wall_pos + dir_sign*(0.5f*slength + (use_sink_model ? wall_thickness : 0.0f)));
		float sink_pos(sink_start);
		cube_t sinks_bcube;

		for (unsigned n = 0; n < num_sinks; ++n) {
			point center(sink_from_wall, sink_pos, zval);
			if (br_dim) {swap(center.x, center.y);} // R90 about z
			cube_t sink(center, center);
			sink.expand_in_dim(br_dim, 0.5*slength);
			sink.z2() += sheight;

			if (use_sink_model) { // sink 3D model
				sink.expand_in_dim(!br_dim, 0.5*swidth);
				objs.emplace_back(sink, TYPE_SINK, room_id, br_dim, !dir, 0, tot_light_amt);
			}
			else { // flat basin sink
				sink.expand_in_dim(!br_dim, 0.5*fabs(sink_step)); // tile exactly with the adjacent sink
				objs.emplace_back(sink, TYPE_BRSINK, room_id, br_dim, !dir, 0, tot_light_amt);
			}
			sinks_bcube.assign_or_union_with_cube(sink);
			sink_pos += sink_step;
		} // for n
		if (add_urinals) { // add urinals opposite the sinks, using same spacing as sinks
			float const u_wall(place_area.d[br_dim][!dir]), u_from_wall(u_wall - dir_sign*(0.5*ulength + 0.01*wall_thickness));
			float u_pos(sink_start);
			cube_t sep_wall;
			set_cube_zvals(sep_wall, zval+0.15*uheight, zval+1.25*uheight);
			sep_wall.d[br_dim][!dir] = u_wall;
			sep_wall.d[br_dim][ dir] = u_wall - dir_sign*0.25*floor_spacing;

			for (unsigned n = 0; n < num_sinks; ++n) {
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, 0, tot_light_amt, SHAPE_SHORT, stall_color);
				point center(u_from_wall, u_pos, (zval + 0.2*uheight));
				if (br_dim) {swap(center.x, center.y);} // R90 about z
				cube_t urinal(center, center);
				urinal.expand_in_dim( br_dim, 0.5*ulength);
				urinal.expand_in_dim(!br_dim, 0.5*uwidth);
				urinal.z2() += uheight;
				objs.emplace_back(urinal, TYPE_URINAL, room_id, br_dim, dir, 0, tot_light_amt);
				u_pos += sink_step;
			} // for n
			if (!two_rows) { // skip first wall if adjacent to a stall
				set_wall_width(sep_wall, (u_pos - 0.5*sink_step), 0.2*wall_thickness, !br_dim);
				objs.emplace_back(sep_wall, TYPE_STALL, room_id, br_dim, !dir, 0, tot_light_amt, SHAPE_SHORT, stall_color);
			}
		}
		if (!sinks_bcube.is_all_zeros()) { // add a long mirror above the sink
			if (!ENABLE_MIRROR_REFLECTIONS || dir != (unsigned)skip_stalls_side) { // don't add mirrors to both sides if reflections are enabled
				cube_t mirror(sinks_bcube);
				mirror.expand_in_dim(!br_dim, -0.25*wall_thickness); // slightly smaller
				mirror.d[br_dim][ dir] = wall_pos;
				mirror.d[br_dim][!dir] = wall_pos + dir_sign*0.1*wall_thickness;
				mirror.z1() = sinks_bcube.z2() + 0.25*floor_thickness;
				mirror.z2() = zval + 0.9*floor_spacing - floor_thickness;
				if (mirror.is_strictly_normalized()) {objs.emplace_back(mirror, TYPE_MIRROR, room_id, br_dim, !dir, RO_FLAG_NOCOLL, tot_light_amt);}
			}
		}
	} // for dir
	// add a sign outside the bathroom door
	bool const shift_dir(room_center[br_dim] < part_center[br_dim]); // put the sign toward the outside of the building because there's more space and more light
	float const door_width(br_door.get_sz_dim(br_dim));
	cube_t sign(br_door);
	set_cube_zvals(sign, zval+0.50*floor_spacing, zval+0.55*floor_spacing);
	sign.translate_dim( br_dim, (shift_dir ? -1.0 : 1.0)*0.8*door_width);
	sign.expand_in_dim( br_dim, -(mens_room ? 0.35 : 0.25)*door_width); // shrink a bit
	sign.translate_dim(!br_dim, sink_side_sign*0.5*wall_thickness); // move to outside wall
	sign.d[!br_dim][sink_side] += sink_side_sign*0.1*wall_thickness; // make nonzero area
	objs.emplace_back(sign, TYPE_SIGN, room_id, !br_dim, sink_side, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CUBE, DK_BLUE); // technically should use hallway room_id
	string const sign_text(mens_room ? "Men" : "Women");
	objs.back().obj_id = register_sign_text(sign_text);
	return 1;
}

void add_door_if_blocker(cube_t const &door, cube_t const &room, bool inc_open, bool dir, bool hinge_side, vect_cube_t &blockers) {
	bool const dim(door.dy() < door.dx()), edir(dim ^ dir ^ hinge_side ^ 1);
	float const width(door.get_sz_dim(!dim));
	cube_t door_exp(door);
	door_exp.expand_in_dim(dim, width);
	if (!door_exp.intersects(room)) return; // check against room before expanding along wall to exclude doors in adjacent rooms
	door_exp.expand_in_dim(!dim, width*0.25); // min expand value
	if (inc_open) {door_exp.d[!dim][edir] += (edir ? 1.0 : -1.0)*0.75*width;} // expand the remainder of the door width in this dir
	blockers.push_back(door_exp);
}
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
	for (auto i = doors.begin(); i != doors.end(); ++i) {add_door_if_blocker(i->get_bcube(), room, 0, 0, 0, blockers);} // exterior doors, inc_open=0

	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) { // interior doors
		add_door_if_blocker(*i, room, door_opens_inward(*i, room), i->open_dir, i->hinge_side, blockers);
	}
	float const doorway_width(get_doorway_width());

	for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
		cube_t tc(*s);
		// expand only in stairs entrance dim for the first floor (could do the opposite for top floor)
		bool const first_floor(room.z1() <= s->z1() + get_floor_thickness()); // for these stairs, not for the building
		if (first_floor) {tc.d[s->dim][!s->dir] += (s->dir ? -1.0 : 1.0);}
		else {tc.expand_in_dim(s->dim, doorway_width);} // add extra space at both ends of stairs
		if (tc.intersects(bcube)) {blockers.push_back(tc);}
	}
	for (auto e = interior->elevators.begin(); e != interior->elevators.end(); ++e) {
		cube_t tc(*e);
		tc.d[e->dim][e->dir] += doorway_width*(e->dir ? 1.0 : -1.0); // add extra space in front of the elevator
		if (tc.intersects(bcube)) {blockers.push_back(tc);}
	}
	return table_blocker_ix;
}

bool building_t::add_kitchen_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool allow_adj_ext_door) {
	// Note: table and chairs have already been placed
	if (room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be kitchens
	if (!is_house && rgen.rand_bool()) return 0; // some office buildings have kitchens, allow it half the time
	// if it has an external door then reject the room half the time; most houses don't have a front door to the kitchen
	if (is_room_adjacent_to_ext_door(room, 1) && (!allow_adj_ext_door || rgen.rand_bool())) return 0; // front_door_only=1
	float const wall_thickness(get_wall_thickness());
	cube_t room_bounds(get_walkable_room_bounds(room)), place_area(room_bounds);
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall for appliances
	bool placed_obj(0);
	placed_obj |= place_model_along_wall(OBJ_MODEL_FRIDGE, TYPE_FRIDGE, room, 0.75, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.2, 4, 0, WHITE, 1); // not at window
	if (is_house) {placed_obj |= place_model_along_wall(OBJ_MODEL_STOVE, TYPE_STOVE, room, 0.46, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 1.0);}
		
	if (is_house && placed_obj) { // if we have at least a fridge or stove, try to add countertops
		float const vspace(get_window_vspace()), height(0.345*vspace), depth(0.74*height), min_hwidth(0.6*height), floor_thickness(get_floor_thickness());
		float const min_clearance(get_min_front_clearance()), front_clearance(max(0.6f*height, min_clearance));
		cube_t cabinet_area(room_bounds);
		cabinet_area.expand_by(-0.05*wall_thickness); // smaller gap than place_area; this is needed to prevent z-fighting with exterior walls
		if (min(cabinet_area.dx(), cabinet_area.dy()) < 4.0*min_hwidth) return placed_obj; // no space for cabinets, room is too small
		vect_room_object_t &objs(interior->room_geom->objs);
		unsigned const counters_start(objs.size());
		cube_t c;
		set_cube_zvals(c, zval, zval+height);
		set_cube_zvals(cabinet_area, zval, (zval + vspace - floor_thickness));
		static vect_cube_t blockers;
		int const table_blocker_ix(gather_room_placement_blockers(cabinet_area, objs_start, blockers, 1, 1)); // inc_open_doors=1, ignore_chairs=1
		bool const have_toaster(building_obj_model_loader.is_model_valid(OBJ_MODEL_TOASTER));
		vector3d const toaster_sz(have_toaster ? building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_TOASTER) : zero_vector); // L, D, H
		bool is_sink(1), placed_mwave(0), placed_toaster(0);
		cube_t mwave, toaster;

		for (unsigned n = 0; n < 50; ++n) { // 50 attempts
			bool const dim(rgen.rand_bool()), dir(rgen.rand_bool()); // choose a random wall
			bool const is_ext_wall(classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT); // assumes not in basement
			// only consider exterior walls in the first 20 attempts to prioritize these so that we don't have splits visible through windows; also places kitchen sinks near windows
			if (n < 20 && !is_ext_wall) continue;
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
			assert(c.contains_cube(c_min));
			c.d[dim][!dir] = front_pos; // remove front clearance
			bool const add_backsplash(!is_ext_wall); // only add to interior walls to avoid windows; assuming not in basement

			for (auto i = objs.begin()+counters_start; i != objs.end(); ++i) { // find adjacencies to previously placed counters and flag to avoid placing doors
				if (i->dim == dim) continue; // not perpendicular
				if (i->d[!i->dim][dir] != wall_pos) continue; // not against the wall on this side
				if (i->d[i->dim][i->dir] != c.d[!dim][0] && i->d[i->dim][i->dir] != c.d[!dim][1]) continue; // not adjacent
				i->flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
				if (add_backsplash) {i->item_flags |= (1 << (dir+1));} // flag side as having a backsplash
			}
			unsigned const cabinet_id(objs.size());
			objs.emplace_back(c, (is_sink ? TYPE_KSINK : TYPE_COUNTER), room_id, dim, !dir, 0, tot_light_amt);
			set_obj_id(objs);
			
			if (add_backsplash) {
				objs.back().item_flags |= 1; // flag back as having a backsplash
				cube_t bs(c);
				bs.z1()  = c.z2();
				bs.z2() += 0.33*c.dz();
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
				placed_mwave = 1;
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
						unsigned const NUM_TOASTER_COLORS = 7;
						colorRGBA const toaster_colors[NUM_TOASTER_COLORS] = {WHITE, LT_GRAY, GRAY, DK_GRAY, GRAY_BLACK, colorRGBA(0.0, 0.0, 0.5), colorRGBA(0.5, 0.0, 0.0)};
						objs.emplace_back(toaster, TYPE_TOASTER, room_id, !dim, rgen.rand_bool(), RO_FLAG_NOCOLL, tot_light_amt); // random dir
						objs.back().color = toaster_colors[rgen.rand()%NUM_TOASTER_COLORS];
						objs[cabinet_id].flags |= RO_FLAG_ADJ_TOP; // flag as having a toaster so that we don't add a book or bottle that could overlap it
						placed_toaster = 1;
					}
				}
			}
			is_sink = 0; // sink is in first placed counter only
		} // for n
	}
	return placed_obj;
}

bool building_t::add_livingroom_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!is_house || room.is_hallway || room.is_sec_bldg || room.is_office) return 0; // these can't be living rooms
	float const wall_thickness(get_wall_thickness());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*wall_thickness); // common spacing to wall for appliances
	vect_room_object_t &objs(interior->room_geom->objs);
	bool placed_couch(0), placed_tv(0);
	// place couches with a variety of colors
	unsigned const NUM_COLORS = 8;
	colorRGBA const colors[NUM_COLORS] = {GRAY_BLACK, WHITE, LT_GRAY, GRAY, DK_GRAY, LT_BROWN, BROWN, DK_BROWN};
	colorRGBA const &couch_color(colors[rgen.rand()%NUM_COLORS]);
	unsigned tv_pref_orient(4), couch_ix(objs.size()), tv_ix(0);
	
	if (place_model_along_wall(OBJ_MODEL_COUCH, TYPE_COUCH, room, 0.40, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.67, 4, 1, couch_color)) { // pref centered
		placed_couch   = 1;
		tv_pref_orient = (2*objs[couch_ix].dim + !objs[couch_ix].dir); // TV should be across from couch
	}
	tv_ix = objs.size();

	// place TV: pref centered; maybe should set not_at_window=1, but that seems too restrictive
	if (place_model_along_wall(OBJ_MODEL_TV, TYPE_TV, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 4.0, tv_pref_orient, 1, BKGRAY, 0)) {
		placed_tv = 1;
		// add a small table to place the TV on so that it's off the floor and not blocked as much by tables and chairs
		room_object_t &tv(objs[tv_ix]);
		float const height(0.4*tv.dz());
		cube_t table(tv); // same XY bounds as the TV
		tv.translate_dim(2, height); // move TV up
		table.z2() = tv.z1();
		objs.emplace_back(table, TYPE_TABLE, room_id, 0, 0, RO_FLAG_IS_HOUSE, tot_light_amt, SHAPE_SHORT); // short table; houses only
	}
	if (placed_couch && placed_tv) {
		room_object_t const &couch(objs[couch_ix]), &tv(objs[tv_ix]);

		if (couch.dim == tv.dim && couch.dir != tv.dir) { // placed against opposite walls facing each other
			cube_t region(couch);
			region.union_with_cube(tv);
			shorten_chairs_in_region(region, objs_start); // region represents that space between the couch and the TV
		}
	}
	return (placed_couch || placed_tv);
}

// Note: this room is decided by the caller and the failure to add objects doesn't make it not a dining room
void building_t::add_diningroom_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	//if (!is_house || room.is_hallway || room.is_sec_bldg || room.is_office) return; // still applies, but unnecessary
	if ((rgen.rand()&3) == 0) return; // no additional objects 25% of the time
	cube_t room_bounds(get_walkable_room_bounds(room));
	room_bounds.expand_by_xy(-get_trim_thickness());
	float const vspace(get_window_vspace()), clearance(max(0.2f*vspace, get_min_front_clearance()));
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
		if (is_cube_close_to_doorway(tc, room, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(tc) || overlaps_other_room_obj(tc, objs_start)) continue; // bad placement
		objs.emplace_back(c, TYPE_WINE_RACK, room_id, dim, !dir, 0, tot_light_amt); // Note: dir faces into the room, not the wall
		set_obj_id(objs);
		break; // done/success
	} // for n
}

bool building_t::add_library_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement) {
	if (room.is_hallway || room.is_sec_bldg) return 0; // these can't be libraries
	unsigned num_added(0);

	for (unsigned n = 0; n < 8; ++n) { // place up to 8 bookcases
		bool const added(add_bookcase_to_room(rgen, room, zval, room_id, tot_light_amt, objs_start, is_basement));
		if (added) {++num_added;} else {break;}
	}
	return (num_added > 0);
}

void gen_crate_sz(vector3d &sz, rand_gen_t &rgen, float window_vspacing) {
	for (unsigned d = 0; d < 3; ++d) {sz[d] = 0.06*window_vspacing*(1.0 + ((d == 2) ? 1.2 : 2.0)*rgen.rand_float());} // slightly more variation in XY
}

bool building_t::add_storage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement) {
	bool const is_garage_or_shed(room.is_garage_or_shed(0)), is_int_garage(room.get_room_type(0) == RTYPE_GARAGE);
	float const window_vspacing(get_window_vspace()), wall_thickness(get_wall_thickness()), floor_thickness(get_floor_thickness());
	float const ceil_zval(zval + window_vspacing - floor_thickness), shelf_depth((is_house ? (is_basement ? 0.18 : 0.15) : 0.2)*window_vspacing);
	float shelf_shorten(shelf_depth + 1.0f*wall_thickness);
	// increase shelf shorten for interior garages to account for approx width of exterior door when opened
	if (is_int_garage) {max_eq(shelf_shorten, 0.36f*window_vspacing);}
	cube_t room_bounds(get_walkable_room_bounds(room)), crate_bounds(room_bounds);
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const num_crates(4 + (rgen.rand() % (is_house ? (is_basement ? 12 : 5) : 30))); // 4-33 for offices, 4-8 for houses, 4-16 for house basements
	vect_cube_t exclude;
	cube_t test_cube(room);
	set_cube_zvals(test_cube, zval, zval+wall_thickness); // reduce to a small z strip for this floor to avoid picking up doors on floors above or below
	unsigned num_placed(0), num_doors(0);

	// first pass to count the number of doors in this room
	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
		num_doors += is_cube_close_to_door(test_cube, 0.0, 0, *i, 2); // check both dirs
	}
	for (auto i = interior->door_stacks.begin(); i != interior->door_stacks.end(); ++i) {
		if (!is_cube_close_to_door(test_cube, 0.0, 0, *i, 2)) continue; // check both dirs
		exclude.push_back(*i);
		exclude.back().expand_in_dim( i->dim, 0.6*room.get_sz_dim(i->dim));
		// if there are multiple doors (houses only?), expand the exclude area more in the other dimension to make sure there's a path between doors
		exclude.back().expand_in_dim(!i->dim, max(0.1*i->get_width(), ((num_doors > 1) ? 0.3*room.get_sz_dim(!i->dim) : 0.0)));
	}
	// add shelves on walls (avoiding any door(s)), and have crates avoid them
	for (unsigned dim = 0; dim < 2; ++dim) {
		if (room_bounds.get_sz_dim( dim) < 6.0*shelf_depth  ) continue; // too narrow to add shelves in this dim
		if (room_bounds.get_sz_dim(!dim) < 4.0*shelf_shorten) continue; // too narrow in the other dim

		for (unsigned dir = 0; dir < 2; ++dir) {
			if (is_int_garage ? ((rgen.rand()%3) == 0) : rgen.rand_bool()) continue; // only add shelves to 50% of the walls, 67% for interior garages
			
			if (is_garage_or_shed) { // garage or shed - don't place shelves in front of door, but allow them against windows
				cube_t wall(room);
				wall.d[dim][!dir] = wall.d[dim][dir]; // shrink room to zero width along this wall
				if (is_room_adjacent_to_ext_door(wall)) continue;
			}
			else if (is_house && !is_basement && has_windows() && classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT) {
				// don't place shelves against exterior house walls in case there are windows
				cube_t const part(get_part_for_room(room));
				float const h_spacing(get_hspacing_for_part(part, !dim));
				if (room_bounds.get_sz_dim(!dim) - 2.0*shelf_depth > h_spacing) continue; // shelf width is larger than spacing - likely to intersect a window, don't test center pt
				if (is_val_inside_window(part, !dim, room_bounds.get_center_dim(!dim), h_spacing, get_window_h_border())) continue;
			}
			cube_t shelves(room_bounds);
			set_cube_zvals(shelves, zval, ceil_zval-floor_thickness);
			crate_bounds.d[dim][dir] = shelves.d[dim][!dir] = shelves.d[dim][dir] + (dir ? -1.0 : 1.0)*shelf_depth; // outer edge of shelves, which is also the crate bounds
			shelves.expand_in_dim(!dim, -shelf_shorten); // shorten shelves
			if (has_bcube_int(shelves, exclude)) continue; // too close to a doorway
			if (!is_garage_or_shed && interior->is_blocked_by_stairs_or_elevator(shelves)) continue;
			unsigned const shelf_flags((is_house ? RO_FLAG_IS_HOUSE : 0) | (is_garage_or_shed ? 0 : RO_FLAG_INTERIOR));
			objs.emplace_back(shelves, TYPE_SHELVES, room_id, dim, dir, shelf_flags, tot_light_amt);
			set_obj_id(objs);
		} // for dir
	} // for dim
	if (is_garage_or_shed) return 1; // no chair, crates, or boxes in garages or sheds

	// add a random office chair if there's space
	if (!is_house && min(crate_bounds.dx(), crate_bounds.dy()) > 1.2*window_vspacing && building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR)) {
		float const chair_height(0.425*window_vspacing), chair_radius(0.5f*chair_height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR));
		point pos(gen_xy_pos_in_area(crate_bounds, chair_radius, rgen));
		pos.z = zval;
		cube_t chair(get_cube_height_radius(pos, chair_radius, chair_height));
		
		// for now, just make one random attempt; if it fails then there's no chair in this room
		if (!has_bcube_int(chair, exclude) && !is_cube_close_to_doorway(chair, room, 0.0, 1) && !interior->is_blocked_by_stairs_or_elevator(chair)) {
			objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, rgen.rand_bool(), rgen.rand_bool(), RO_FLAG_RAND_ROT, tot_light_amt, SHAPE_CYLIN, GRAY_BLACK);
		}
	}
	for (unsigned n = 0; n < 4*num_crates; ++n) { // make up to 4 attempts for every crate/box
		vector3d sz; // half size relative to window_vspacing
		gen_crate_sz(sz, rgen, window_vspacing*(is_house ? (is_basement ? 0.75 : 0.5) : 1.0)); // smaller for houses
		if (crate_bounds.dx() <= 2.0*sz.x || crate_bounds.dy() <= 2.0*sz.y) continue; // too large for this room
		point pos(gen_xy_pos_in_area(crate_bounds, sz, rgen));
		pos.z = zval;
		cube_t crate(get_cube_height_radius(pos, sz, 2.0*sz.z)); // multiply by 2 since this is a size rather than half size/radius
		if (has_bcube_int(crate, exclude)) continue; // don't place crates between the door and the center of the room
		bool bad_placement(0);

		for (auto i = objs.begin()+objs_start; i != objs.end(); ++i) {
			if (!i->intersects(crate)) continue;
			// only handle stacking of crates on other crates
			if ((i->type == TYPE_CRATE || i->type == TYPE_BOX) && i->z1() == zval && (i->z2() + crate.dz() < ceil_zval) && i->contains_pt_xy(pos)) {crate.translate_dim(2, i->dz());}
			else {bad_placement = 1; break;}
		}
		if (bad_placement) continue;
		if (is_cube_close_to_doorway(crate, room, 0.0, 1) || interior->is_blocked_by_stairs_or_elevator(crate)) continue;
		cube_t c2(crate);
		c2.expand_by(vector3d(0.5*c2.dx(), 0.5*c2.dy(), 0.0)); // approx extents of flaps if open
		unsigned flags(0);
		
		for (unsigned d = 0; d < 4; ++d) { // determine which sides are against a wall
			bool const dim(d>>1), dir(d&1);
			if ((c2.d[dim][dir] < room_bounds.d[dim][dir]) ^ dir) {flags |= (RO_FLAG_ADJ_LO << d);}
		}
		objs.emplace_back(crate, (rgen.rand_bool() ? TYPE_CRATE : TYPE_BOX), room_id, rgen.rand_bool(), 0, flags, tot_light_amt, SHAPE_CUBE, gen_box_color(rgen)); // crate or box
		set_obj_id(objs); // used to select texture and box contents
		if (++num_placed == num_crates) break; // we're done
	} // for n
	return 1; // it's always a storage room, even if it's empty
}

void building_t::add_garage_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	if (!enable_parked_cars() || (rgen.rand()&3) == 0) return; // 75% of garages have cars
	unsigned const flags(RO_FLAG_NOCOLL | RO_FLAG_USED | RO_FLAG_INVIS); // lines not shown
	bool const dim(room.dx() < room.dy()); // long dim
	bool dir(0); // set dir so that cars pull into driveways
	if (street_dir > 0 && bool((street_dir-1)>>1) == dim) {dir = !((street_dir-1)&1);} // use street_dir if it's set and dims agree
	else {dir = (room.get_center_dim(dim) < bcube.get_center_dim(dim));} // assumes the garage is at an exterior wall and doesn't occupy the entire house width
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t space(room); // full room, car will be centered here
	set_cube_zvals(space, zval, (zval + 0.001*get_window_vspace()));
	room_object_t pspace(space, TYPE_PARK_SPACE, room_id, dim, dir, flags, tot_light_amt, SHAPE_CUBE, WHITE);
	pspace.obj_id = (uint16_t)(objs.size() + rgen.rand()); // will be used for the car model and color
	car_t car(car_from_parking_space(pspace));
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

bool building_t::add_office_utility_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	// TODO
	return 0;
}

bool building_t::add_laundry_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned &added_bathroom_objs_mask) {
	float const front_clearance(get_min_front_clearance());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // common spacing to wall for appliances
	vector3d const place_area_sz(place_area.get_size());
	vect_room_object_t &objs(interior->room_geom->objs);

	for (unsigned n = 0; n < 10; ++n) { // 10 attempts to place washer and dryer along the same wall
		unsigned const washer_ix(objs.size());
		bool const placed_washer(place_model_along_wall(OBJ_MODEL_WASHER, TYPE_WASHER, room, 0.42, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8));
		unsigned pref_orient(4); // if washer was placed, prefer to place dryer along the same wall
		if (placed_washer) {pref_orient = objs[washer_ix].get_orient();}
		unsigned const dryer_ix(objs.size());
		bool const placed_dryer(place_model_along_wall(OBJ_MODEL_DRYER, TYPE_DRYER, room, 0.38, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.8, pref_orient));
		bool success(0);
		if (placed_washer && placed_dryer && objs[dryer_ix].get_orient() == pref_orient) {success = 1;} // placed both washer and dryer along the same wall
		else if (n+1 == 10) { // last attempt
			if (!(placed_washer || placed_dryer)) return 0; // placed neither washer nor dryer, failed
			if (placed_washer != placed_dryer) {success = 1;} // placed only one of the washer or dryer, allow it
			else if (objs[washer_ix].dim != objs[dryer_ix].dim) {success = 1;} // placed on two adjacent walls, allow it
			// placed on opposite walls; check that there's space for the player to walk between the washer and dryer
			else if (objs[washer_ix].get_sz_dim(objs[washer_ix].dim) + objs[dryer_ix].get_sz_dim(objs[dryer_ix].dim) + front_clearance < place_area_sz[objs[washer_ix].dim]) {success = 1;}
		}
		if (success) {
			// if we've placed a washer and/or dryer and made this into a laundry room, try to place a sink as well; should this use a different sink model from bathrooms?
			if (place_model_along_wall(OBJ_MODEL_SINK, TYPE_SINK, room, 0.45, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.6)) {
				added_bathroom_objs_mask |= PLACED_SINK;
			}
			// try to place a laundry basket
			float const floor_spacing(get_window_vspace()), radius(rgen.rand_uniform(0.1, 0.12)*floor_spacing), height(rgen.rand_uniform(1.5, 2.2)*radius);
			place_area.expand_by_xy(-radius); // leave a slight gap between laundry basket and wall
			if (!place_area.is_strictly_normalized()) return 1; // no space for laundry basket (likely can't happen)
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
				if (is_cube_close_to_doorway(c, room, 0.0, !room.is_hallway) || interior->is_blocked_by_stairs_or_elevator(c) || overlaps_other_room_obj(c, objs_start)) continue; // bad placement
				colorRGBA const colors[4] = {WHITE, LT_BLUE, LT_GREEN, LT_BROWN};
				objs.emplace_back(c, TYPE_LBASKET, room_id, dim, dir, 0, tot_light_amt, SHAPE_CYLIN, colors[rgen.rand()%4]);
				break; // done
			} // for n
			return 1; // done
		}
		objs.resize(objs_start); // remove washer and dryer and try again
	} // for n
	return 0; // failed
}

void building_t::add_pri_hall_objs(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt) {
	float const window_vspacing(get_window_vspace()), desk_width(0.9*window_vspacing);
	bool const long_dim(room.dx() < room.dy());
	if (room.get_sz_dim(!long_dim) < (desk_width + 1.6*get_doorway_width())) return; // hallway is too narrow
	float const centerline(room.get_center_dim(!long_dim)), desk_depth(0.6*desk_width);
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t desk;
	set_cube_zvals(desk, zval, zval+0.32*window_vspacing);
	set_wall_width(desk, centerline, 0.5*desk_width, !long_dim);

	for (unsigned dir = 0; dir < 2; ++dir) { // add a reception desk at each entrance
		float const hall_len(room.get_sz_dim(long_dim)), hall_start(room.d[long_dim][dir]), dir_sign(dir ? -1.0 : 1.0);
		float const val1(hall_start + max(0.1f*hall_len, window_vspacing)*dir_sign), val2(hall_start + 0.3*hall_len*dir_sign); // range of reasonable desk placements along the hall

		for (unsigned n = 0; n < 10; ++n) { // try to find the closest valid placement to the door, make 10 random attempts
			float const val(rgen.rand_uniform(min(val1, val2), max(val1, val2)));
			set_wall_width(desk, val, 0.5*desk_depth, long_dim);
			if (interior->is_blocked_by_stairs_or_elevator(desk)) continue; // bad location, try a new one

			if (building_obj_model_loader.is_model_valid(OBJ_MODEL_OFFICE_CHAIR)) {
				float const chair_height(0.425*window_vspacing), chair_radius(0.5f*chair_height*get_radius_for_square_model(OBJ_MODEL_OFFICE_CHAIR));
				point pos;
				pos.z = zval;
				pos[!long_dim] = centerline;
				pos[ long_dim] = val + dir_sign*(-0.05*desk_depth + chair_radius); // push the chair into the cutout of the desk
				cube_t const chair(get_cube_height_radius(pos, chair_radius, chair_height));
				if (interior->is_blocked_by_stairs_or_elevator(chair)) continue; // bad location, try a new one
				objs.emplace_back(chair, TYPE_OFF_CHAIR, room_id, long_dim, dir, 0, tot_light_amt, SHAPE_CYLIN, GRAY_BLACK);
			}
			objs.emplace_back(desk, TYPE_RDESK, room_id, long_dim, dir, 0, tot_light_amt, SHAPE_CUBE);
			break; // done
		} // for n
	} // for dir
}

colorRGBA choose_pot_color(rand_gen_t &rgen) {
	unsigned const num_colors = 8;
	colorRGBA const pot_colors[num_colors] = {LT_GRAY, GRAY, DK_GRAY, BKGRAY, WHITE, LT_BROWN, RED, colorRGBA(1.0, 0.35, 0.18)};
	return pot_colors[rgen.rand() % num_colors];
}

void building_t::place_book_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, bool use_dim_dir) {
	point center(place_on.get_cube_center());
	for (unsigned d = 0; d < 2; ++d) {center[d] += 0.1*place_on.get_sz_dim(d)*rgen.rand_uniform(-1.0, 1.0);} // add a slight random shift
	float const book_sz(0.07*get_window_vspace());
	// book is randomly oriented for tables and rotated 90 degrees from desk orient
	bool const dim(use_dim_dir ? !place_on.dim : rgen.rand_bool()), dir(use_dim_dir ? (place_on.dir^place_on.dim) : rgen.rand_bool());
	cube_t book;
	vector3d book_scale(book_sz*rgen.rand_uniform(0.8, 1.2), book_sz*rgen.rand_uniform(0.8, 1.2), 0.0);
	book_scale[dim] *= 0.8; // slightly smaller in this dim
	book.set_from_point(point(center.x, center.y, place_on.z2()));
	book.expand_by(book_scale);
	book.z2() += book_sz*rgen.rand_uniform(0.1, 0.3);
	colorRGBA const color(book_colors[rgen.rand() % NUM_BOOK_COLORS]);
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(book, TYPE_BOOK, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt, SHAPE_CUBE, color); // Note: invalidates place_on reference
	set_obj_id(objs);
}

cube_t place_cylin_object(rand_gen_t rgen, cube_t const &place_on, float radius, float height, float dist_from_edge) {
	cube_t c;
	gen_xy_pos_for_round_obj(c, place_on, radius, height, dist_from_edge, rgen); // place at dist_from_edge from edge
	return c;
}

bool building_t::place_bottle_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid) {
	float const window_vspacing(get_window_vspace());
	float const height(window_vspacing*rgen.rand_uniform(0.075, 0.12)), radius(window_vspacing*rgen.rand_uniform(0.012, 0.018));
	if (min(place_on.dx(), place_on.dy()) < 6.0*radius) return 0; // surface is too small to place this bottle
	cube_t const bottle(place_cylin_object(rgen, place_on, radius, height, 2.0*radius));
	if (!avoid.is_all_zeros() && bottle.intersects(avoid)) return 0; // only make one attempt
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(bottle, TYPE_BOTTLE, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN);
	objs.back().set_as_bottle(rgen.rand(), 3); // 0-3; excludes poison and medicine
	return 1;
}

bool building_t::place_plant_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid) {
	float const window_vspacing(get_window_vspace()), height(rgen.rand_uniform(0.25, 0.4)*window_vspacing);
	float const radius(min(rgen.rand_uniform(0.06, 0.08)*window_vspacing, min(place_on.dx(), place_on.dy())/3.0f));
	cube_t const plant(place_cylin_object(rgen, place_on, radius, height, 1.2*radius));
	if (!avoid.is_all_zeros() && plant.intersects(avoid)) return 0; // only make one attempt
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(plant, TYPE_PLANT, room_id, 0, 0, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_BOT), tot_light_amt, SHAPE_CYLIN, choose_pot_color(rgen));
	set_obj_id(objs);
	return 1;
}

bool building_t::place_laptop_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid, bool use_dim_dir) {
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
	if (!avoid.is_all_zeros() && laptop.intersects(avoid)) return 0; // only make one attempt
	interior->room_geom->objs.emplace_back(laptop, TYPE_LAPTOP, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt); // Note: invalidates place_on reference
	return 1;
}

float get_plate_radius(rand_gen_t &rgen, room_object_t const &place_on, float window_vspacing) {
	return min(rgen.rand_uniform(0.05, 0.07)*window_vspacing, 0.25f*min(place_on.dx(), place_on.dy()));
}

bool building_t::place_plate_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid) {
	float const radius(get_plate_radius(rgen, place_on, get_window_vspace()));
	cube_t const plate(place_cylin_object(rgen, place_on, radius, 0.1*radius, 1.1*radius));
	if (!avoid.is_all_zeros() && plate.intersects(avoid)) return 0; // only make one attempt
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(plate, TYPE_PLATE, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt, SHAPE_CYLIN);
	set_obj_id(objs);
	return 1;
}

bool building_t::place_cup_on_obj(rand_gen_t &rgen, room_object_t const &place_on, unsigned room_id, float tot_light_amt, cube_t const &avoid) {
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_CUP)) return 0;
	float const window_vspacing(get_window_vspace()), height(0.06*window_vspacing), radius(0.5f*height*get_radius_for_square_model(OBJ_MODEL_CUP)); // almost square
	cube_t const cup(place_cylin_object(rgen, place_on, radius, height, 1.2*radius));
	if (!avoid.is_all_zeros() && cup.intersects(avoid)) return 0; // only make one attempt
	// random dim/dir, plus more randomness on top
	interior->room_geom->objs.emplace_back(cup, TYPE_CUP, room_id, rgen.rand_bool(), rgen.rand_bool(), (RO_FLAG_NOCOLL | RO_FLAG_RAND_ROT), tot_light_amt, SHAPE_CYLIN);
	return 1;
}

bool building_t::add_rug_to_room(rand_gen_t rgen, cube_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start) {
	if (!room_object_t::enable_rugs()) return 0; // disabled
	vector3d const room_sz(room.get_size());
	bool const min_dim(room_sz.y < room_sz.x);
	float const ar(rgen.rand_uniform(0.65, 0.85)), length(min(0.7f*room_sz[min_dim]/ar, room_sz[!min_dim]*rgen.rand_uniform(0.4, 0.7))), width(length*ar);
	cube_t rug;
	set_cube_zvals(rug, zval, zval+0.001*get_window_vspace()); // almost flat
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
			if (!i->intersects(rug)) continue;

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
				break;
			}
			else if (i->type == TYPE_TABLE || i->type == TYPE_DESK) { // rugs can't partially overlap these object types
				valid_placement = rug.contains_cube_xy(*i); // don't expand as that could cause the rug to intersect a previous object
				break;
				// maybe beds should be included as well, but then rugs are unlikely to be placed in bedrooms
			}
		} // for i
		if (valid_placement) {
			rug.intersect_with_cube(room); // make sure the rug stays within the room bounds
			objs.emplace_back(rug, TYPE_RUG, room_id, 0, 0, RO_FLAG_NOCOLL, tot_light_amt);
			objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix); // determines rug texture
			return 1;
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
	keepout.z1() = zval; // extend to the floor
	keepout.d[dim][!dir] += (dir ? -1.0 : 1.0)*clearance;
	keepout.expand_in_dim(!dim, side_clearance); // make sure there's space for the frame
	if (overlaps_other_room_obj(keepout, objs_start)) return 0;
	bool const inc_open(!is_house && !room.is_office);
	if (is_cube_close_to_doorway(tc, room, 0.0, inc_open)) return 0; // bad placement
	// Note: it's not legal to guard the below check with (room.has_stairs || room.has_elevator) because room.has_stairs may not be set for stack connector stairs that split a wall
	if (interior->is_blocked_by_stairs_or_elevator(tc, 4.0*wall_thickness)) return 0; // check stairs and elevators
	if (!inc_open && !room.is_hallway && is_cube_close_to_doorway(tc, room, 0.0, 1)) return 2; // success, but could be better (doors never open into hallway)

	if (has_complex_floorplan) { // check for office building whiteboards placed on room sides that aren't true walls
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

bool building_t::hang_pictures_in_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, bool is_basement) {
	if (!room_object_t::enable_pictures()) return 0; // disabled
	
	if (!is_house && !room.is_office) {
		if (room.is_hallway) return 0; // no pictures or whiteboards in office building hallways (what about rooms with stairs?)
		// room in a commercial building - add whiteboard when there is a full wall to use
	}
	if (room.is_sec_bldg) return 0; // no pictures in secondary buildings
	if (room.get_room_type(0) == RTYPE_STORAGE) return 0; // no pictures or whiteboards in storage rooms (always first floor)
	cube_t const &part(get_part_for_room(room));
	float const floor_height(get_window_vspace()), wall_thickness(get_wall_thickness());
	bool const check_for_windows(!is_basement && has_windows());
	vect_room_object_t &objs(interior->room_geom->objs);
	bool was_hung(0);

	if (!is_house || room.is_office) { // add whiteboards
		if (rgen.rand_float() < 0.2) return 0; // skip 20% of the time
		bool const pref_dim(rgen.rand_bool()), pref_dir(rgen.rand_bool());
		float const floor_thick(get_floor_thickness());

		for (unsigned dim2 = 0; dim2 < 2; ++dim2) {
			for (unsigned dir2 = 0; dir2 < 2; ++dir2) {
				bool const dim(bool(dim2) ^ pref_dim), dir(bool(dir2) ^ pref_dir);
				if (check_for_windows && fabs(room.d[dim][dir] - part.d[dim][dir]) < 1.1*wall_thickness) continue; // on part boundary, likely exterior wall where there may be windows, skip
				cube_t c(room);
				set_cube_zvals(c, zval+0.25*floor_height, zval+0.9*floor_height-floor_thick);
				c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*0.6*wall_thickness; // Note: offset by an additional half wall thickness
				c.expand_in_dim(!dim, -0.2*room.get_sz_dim(!dim)); // xy_space
				if (!check_valid_picture_placement(room, c, 0.6*room.get_sz_dim(!dim), zval, dim, dir, objs_start)) continue;
				objs.emplace_back(c, TYPE_WBOARD, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // whiteboard faces dir opposite the wall
				return 1; // done, only need to add one
			} // for dir
		} // for dim
		return 0;
	}
	for (unsigned dim = 0; dim < 2; ++dim) {
		for (unsigned dir = 0; dir < 2; ++dir) {
			float const wall_pos(room.d[dim][dir]);
			if (check_for_windows && fabs(room.d[dim][dir] - part.d[dim][dir]) < 1.1*wall_thickness) continue; // on part boundary, likely exterior wall where there may be windows, skip
			if (!room.is_hallway && rgen.rand_float() < 0.2) continue; // skip 20% of the time unless it's a hallway
			float const height(floor_height*rgen.rand_uniform(0.3, 0.6)), width(height*rgen.rand_uniform(1.5, 2.0)); // width > height
			if (width > 0.8*room.get_sz_dim(!dim)) continue; // not enough space
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
				cube_t c(center, center);
				c.expand_in_dim(2, 0.5*height);
				c.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.1*wall_thickness; // move out to prevent z-fighting
				if (room.is_hallway) {c.translate_dim(dim, (dir ? -1.0 : 1.0)*0.5*wall_thickness);} // add an additional half wall thickness for hallways
				c.expand_in_dim(!dim, 0.5*width);
				int const ret(check_valid_picture_placement(room, c, width, zval, dim, dir, objs_start));
				if (ret == 0) continue; // invalid, retry
				best_pos = c;
				if (ret == 1) break; // valid and good - keep this pos
			} // for n
			if (best_pos.is_all_zeros()) continue; // failed placement
			objs.emplace_back(best_pos, TYPE_PICTURE, room_id, dim, !dir, RO_FLAG_NOCOLL, tot_light_amt); // picture faces dir opposite the wall
			objs.back().obj_id = uint16_t(objs.size() + 13*room_id + 31*mat_ix + 61*dim + 123*dir); // determines picture texture
			was_hung = 1;
		} // for dir
	} // for dim
	return was_hung;
}

void building_t::add_plants_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned num) {
	float const window_vspacing(get_window_vspace());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-get_trim_thickness()); // shrink to leave a small gap
	zval += 0.01*get_floor_thickness(); // move up slightly to avoid z-fithing of bottom when the dirt is taken
	
	for (unsigned n = 0; n < num; ++n) {
		float const height(rgen.rand_uniform(0.6, 0.9)*window_vspacing), width(rgen.rand_uniform(0.15, 0.35)*window_vspacing);
		vector3d const sz_scale(width/height, width/height, 1.0);
		place_obj_along_wall(TYPE_PLANT, room, height, sz_scale, rgen, zval, room_id, tot_light_amt,
			place_area, objs_start, 0.0, 4, 0, choose_pot_color(rgen), 0, SHAPE_CYLIN); // no clearance, pref_orient, or color
	}
}

void building_t::add_boxes_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned max_num) {
	if (max_num == 0) return; // why did we call this?
	float const window_vspacing(get_window_vspace());
	cube_t place_area(get_walkable_room_bounds(room));
	place_area.expand_by(-0.25*get_wall_thickness()); // shrink to leave a small gap
	unsigned const num(rgen.rand() % (max_num+1));

	for (unsigned n = 0; n < num; ++n) {
		vector3d sz;
		gen_crate_sz(sz, rgen, window_vspacing);
		sz *= 1.5; // make larger than storage room boxes
		place_obj_along_wall(TYPE_BOX, room, sz.z, sz, rgen, zval, room_id, tot_light_amt, place_area, objs_start, 0.0, 4, 0, gen_box_color(rgen));
	} // for n
}

void building_t::add_light_switches_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool is_ground_floor) {
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), switch_thickness(0.2*wall_thickness);
	float const switch_height(1.8*wall_thickness), switch_hwidth(0.5*wall_thickness), min_wall_spacing(switch_hwidth + 2.0*wall_thickness);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	if (min(room_bounds.dx(), room_bounds.dy()) < 8.0*switch_hwidth) return; // room is too small; shouldn't happen
	vect_door_stack_t &doorways(get_doorways_for_room(room, zval)); // place light switch next to a door
	if (doorways.size() > 1 && rgen.rand_bool()) {std::reverse(doorways.begin(), doorways.end());} // random permute if more than 2 doorways?
	vect_room_object_t &objs(interior->room_geom->objs);
	unsigned const objs_end(objs.size());
	bool const first_side(rgen.rand_bool());
	vect_door_stack_t ext_doors; // not really door stacks, but we can fill in the data to treat them as such
	cube_t c;
	c.z1() = zval + 0.38*floor_spacing; // same for every switch
	c.z2() = c.z1() + switch_height;

	if (is_ground_floor) { // handle exterior doors
		cube_t room_exp(room);
		room_exp.expand_by(wall_thickness, wall_thickness, -wall_thickness); // expand in XY and shrink in Z

		for (auto d = doors.begin(); d != doors.end(); ++d) {
			if (!d->is_exterior_door() || d->type == tquad_with_ix_t::TYPE_RDOOR) continue;
			cube_t bc(d->get_bcube());
			if (!room_exp.contains_pt(bc.get_cube_center())) continue;
			bool const dim(bc.dy() < bc.dx());
			bc.expand_in_dim(dim, 0.4*wall_thickness); // expand slightly to make it nonzero area
			ext_doors.emplace_back(door_t(bc, dim, 0), 0); // dir=0, first_door_ix=0 because it's unused
		}
	}
	for (unsigned ei = 0; ei < 2; ++ei) { // exterior, interior
		vect_door_stack_t const &cands(ei ? doorways : ext_doors);
		unsigned const max_ls(is_house ? 2 : 1); // place up to 2 light switches in this room if it's a house, otherwise place only 1
		unsigned num_ls(0);

		for (auto i = cands.begin(); i != cands.end() && num_ls < max_ls; ++i) {
			// check for windows if (real_num_parts > 1)? is it actually possible for doors to be within far_spacing of a window?
			bool const dim(i->dim), dir(i->get_center_dim(dim) > room.get_center_dim(dim));
			float const door_width(i->get_width()), near_spacing(0.25*door_width), far_spacing(1.25*door_width); // off to the side of the door when open
			assert(door_width > 0.0);
			cube_t const &wall_bounds(ei ? room_bounds : room); // exterior door should use the original room, not room_bounds
			c.d[dim][ dir] = wall_bounds.d[dim][dir]; // flush with wall
			c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*switch_thickness; // expand out a bit
			bool done(0);

			for (unsigned Side = 0; Side < 2 && !done; ++Side) { // try both sides of the doorway
				bool const side(bool(Side) ^ first_side);

				for (unsigned nf = 0; nf < 2; ++nf) { // {near, far}
					float const spacing(nf ? far_spacing : near_spacing), wall_pos(i->d[!dim][side] + (side ? 1.0 : -1.0)*spacing);
					if (wall_pos < room_bounds.d[!dim][0] + min_wall_spacing || wall_pos > room_bounds.d[!dim][1] - min_wall_spacing) continue; // too close to the adjacent wall
					set_wall_width(c, wall_pos, switch_hwidth, !dim);
					cube_t c_test(c);
					c_test.d[dim][!dir] += (dir ? -1.0 : 1.0)*wall_thickness; // expand out more so that it's guaranteed to intersect appliances placed near the wall
					if (overlaps_other_room_obj(c_test, objs_start))        continue;
					if (is_cube_close_to_doorway(c, room, 0.0, (ei==1), 1)) continue; // inc_open=1/check_open_dir=1 for inside, to avoid placing light switch behind an open door
					if (interior->is_blocked_by_stairs_or_elevator(c))      continue; // check stairs and elevators
					if (!check_of_placed_on_interior_wall(c, room, dim, dir)) continue; // ensure the switch is on a wall
					objs.emplace_back(c, TYPE_SWITCH, room_id, dim, dir, RO_FLAG_NOCOLL, 1.0); // dim/dir matches wall; fully lit
					done = 1; // done, only need to add one for this door
					++num_ls;
					break;
				} // for nf
			} // for side
		} // for i
	} // for ei
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
		c.d[dim][ dir] = target_wall.d[dim][!dir]; // flush with wall
		c.d[dim][!dir] = c.d[dim][dir] + (dir ? -1.0 : 1.0)*switch_thickness; // expand out a bit
		set_wall_width(c, target_wall.get_center_dim(!dim), switch_hwidth, !dim);
		// since nothing is placed against the exterior wall of the closet near the door (to avoid blocking it), we don't need to check for collisions with room objects
		objs.emplace_back(c, TYPE_SWITCH, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_IN_CLOSET), 1.0); // dim/dir matches wall; fully lit; flag for closet
		//break; // there can be only one closet per room; done (unless I add multiple closets later?)
	} // for i
}

void building_t::add_outlets_to_room(rand_gen_t rgen, room_t const &room, float zval, unsigned room_id, unsigned objs_start, bool is_ground_floor, bool is_basement) {
	float const wall_thickness(get_wall_thickness());
	float const plate_thickness(0.03*wall_thickness), plate_height(1.8*wall_thickness), plate_hwidth(0.5*wall_thickness), min_wall_spacing(4.0*plate_hwidth);
	cube_t const room_bounds(get_walkable_room_bounds(room));
	if (min(room_bounds.dx(), room_bounds.dy()) < 3.0*min_wall_spacing) return; // room is too small; shouldn't happen
	vect_door_stack_t const &doorways(get_doorways_for_room(room, zval));
	cube_t c;
	c.z1() = zval + get_trim_height() + 0.4*plate_height; // wall trim height + some extra padding; same for every outlet
	c.z2() = c.z1() + plate_height;

	// try to add an outlet to each wall, down near the floor so that they don't intersect objects such as pictures
	for (unsigned wall = 0; wall < 4; ++wall) {
		bool const dim(wall >> 1), dir(wall & 1);
		if (!is_house && room.get_sz_dim(!dim) < room.get_sz_dim(dim)) continue; // only add outlets to the long walls of office building rooms
		bool const is_exterior_wall(classify_room_wall(room, zval, dim, dir, 0) == ROOM_WALL_EXT); // includes basement
		cube_t const &wall_bounds(is_exterior_wall ? room : room_bounds); // exterior wall should use the original room, not room_bounds
		float const wall_pos(rgen.rand_uniform((room_bounds.d[!dim][0] + min_wall_spacing), (room_bounds.d[!dim][1] - min_wall_spacing)));
		float const wall_face(wall_bounds.d[dim][dir]);
		c.d[dim][ dir] = wall_face; // flush with wall
		c.d[dim][!dir] = wall_face + (dir ? -1.0 : 1.0)*plate_thickness; // expand out a bit
		set_wall_width(c, wall_pos, plate_hwidth, !dim);

		if (!is_basement && has_windows() && is_exterior_wall) { // check for window intersection
			cube_t const part(get_part_for_room(room));
			float const window_hspacing(get_hspacing_for_part(part, !dim)), window_h_border(get_window_h_border());
			// expand by the width of the window trim, plus some padded wall plate width, then check to the left and right;
			// 2*xy_expand should be smaller than a window so we can't have a window fit in between the left and right sides
			float const xy_expand(get_trim_thickness() + 1.2f*plate_hwidth);
			if (is_val_inside_window(part, !dim, (wall_pos - xy_expand), window_hspacing, window_h_border) ||
				is_val_inside_window(part, !dim, (wall_pos + xy_expand), window_hspacing, window_h_border)) continue;
		}
		cube_t c_exp(c);
		c_exp.expand_by_xy(0.5*wall_thickness);
		// Note: outlets can still be partially blocked by picture frames; I guess this is okay?
		if (overlaps_other_room_obj(c_exp, objs_start, 1))     continue; // check for things like closets; check_all=1 to include blinds
		if (interior->is_blocked_by_stairs_or_elevator(c_exp)) continue; // check stairs and elevators
		bool bad_place(0);

		if (is_ground_floor) { // handle exterior doors
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
		if (!check_of_placed_on_interior_wall(c, room, dim, dir)) continue; // ensure the outlet is on a wall
		// Note: it may be more efficient to have outlets be stored as static quads rather than objects, since there's currently no player or AI interaction with them;
		// in fact, maybe wall trim should be the same way since there's so much of it?
		interior->room_geom->objs.emplace_back(c, TYPE_OUTLET, room_id, dim, dir, RO_FLAG_NOCOLL, 1.0); // dim/dir matches wall; fully lit
	} // for wall
}

bool building_t::check_of_placed_on_interior_wall(cube_t const &c, room_t const &room, bool dim, bool dir) const {
	if (!has_small_part && (is_house || !room.is_hallway)) return 1; // check not needed in this case, any non-door location is a wall
	float const wall_thickness(get_wall_thickness()), wall_face(c.d[dim][dir]);
	cube_t test_cube(c);
	test_cube.d[dim][0] = test_cube.d[dim][1] = wall_face - (dir ? -1.0 : 1.0)*0.5*wall_thickness; // move inward
	test_cube.expand_in_dim(!dim, 0.5*wall_thickness);
	// check for exterior wall
	bool intersects_part(0);

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		if (p->intersects(test_cube)) {intersects_part = 1; break;}
	}
	if (!intersects_part) return 1; // not contained in a part, must be an exterior wall
	// check for interior wall
	for (auto const &w : interior->walls[dim]) {
		if (w.contains_cube(test_cube)) return 1;
	}
	return 0;
}

bool building_t::place_eating_items_on_table(rand_gen_t &rgen, unsigned table_obj_id) {
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(table_obj_id < objs.size());
	room_object_t const &table(objs[table_obj_id]);
	float const radius(get_plate_radius(rgen, table, get_window_vspace())), height(0.1*radius), spacing(1.33*radius);
	bool added_obj(0);

	for (auto i = (objs.begin() + table_obj_id + 1); i != objs.end(); ++i) {
		if (i->type != TYPE_CHAIR) break; // done with chairs for this table
		point const chair_center(i->get_cube_center());
		point pos;

		if (table.shape == SHAPE_CYLIN) { // circular
			float const dist(table.get_radius() - spacing);
			point const table_center(table.get_cube_center());
			pos = table_center + dist*(chair_center - table_center).get_norm();
		}
		else { // rectangular
			cube_t place_bounds(table);
			place_bounds.expand_by_xy(-spacing);
			pos = place_bounds.closest_pt(chair_center);
		}
		cube_t plate;
		plate.set_from_sphere(pos, radius);
		set_cube_zvals(plate, table.z2(), table.z2()+height); // place on the table
		objs.emplace_back(plate, TYPE_PLATE, table.room_id, 0, 0, RO_FLAG_NOCOLL, table.light_amt, SHAPE_CYLIN);
		set_obj_id(objs);
		added_obj = 1;
	} // for i
	return added_obj;
}

void building_t::place_objects_onto_surfaces(rand_gen_t rgen, room_t const &room, unsigned room_id, float tot_light_amt, unsigned objs_start, unsigned floor, bool is_basement) {
	if (room.is_hallway) return; // no objects placed in hallways, but there shouldn't be any surfaces either (except for reception desk?)
	vect_room_object_t &objs(interior->room_geom->objs);
	assert(objs.size() > objs_start);
	bool const is_library(room.get_room_type(floor) == RTYPE_LIBRARY);
	bool const is_kitchen(room.get_room_type(floor) == RTYPE_KITCHEN);
	bool const sparse_place(floor > 0 && interior->rooms.size() > 40); // fewer objects on upper floors of large office buildings as an optimization
	float const place_book_prob(( is_house ? 1.0 : 0.5)*(room.is_office ? 0.80 : 1.00)*(sparse_place ? 0.75 : 1.0));
	float const place_bottle_prob(is_house ? 1.0 :      (room.is_office ? 0.80 : 0.50)*(sparse_place ? 0.50 : 1.0));
	float const place_cup_prob   (is_house ? 1.0 :      (room.is_office ? 0.50 : 0.25)*(sparse_place ? 0.50 : 1.0));
	float const place_plant_prob (is_house ? 1.0 :      (room.is_office ? 0.25 : 0.15)*(sparse_place ? 0.75 : 1.0));
	float const place_laptop_prob(is_house ? 0.4 :      (room.is_office ? 0.60 : 0.50)*(sparse_place ? 0.80 : 1.0));
	unsigned const objs_end(objs.size());
	bool placed_book_on_counter(0);

	// see if we can place objects on any room object top surfaces
	for (unsigned i = objs_start; i < objs_end; ++i) { // can't iterate over objs because we modify it
		room_object_t const &obj(objs[i]);
		// add place settings to kitchen and dining room tables 50% of the time
		bool const is_eating_table(obj.type == TYPE_TABLE && (room.get_room_type(floor) == RTYPE_KITCHEN || room.get_room_type(floor) == RTYPE_DINING) && rgen.rand_bool());
		if (is_eating_table && place_eating_items_on_table(rgen, i)) continue; // no other items to place
		float book_prob(0.0), bottle_prob(0.0), cup_prob(0.0), plant_prob(0.0), laptop_prob(0.0);
		cube_t avoid;

		if (obj.type == TYPE_TABLE && i == objs_start) { // only first table (not TV table)
			book_prob   = 0.4*place_book_prob;
			bottle_prob = 0.6*place_bottle_prob;
			cup_prob    = 0.5*place_cup_prob;
			plant_prob  = 0.6*place_plant_prob;
			laptop_prob = 0.3*place_laptop_prob;
		}
		else if (obj.type == TYPE_DESK && (i+1 == objs_end || objs[i+1].type != TYPE_MONITOR)) { // desk with no computer monitor
			book_prob   = 0.8*place_book_prob;
			bottle_prob = 0.4*place_bottle_prob;
			cup_prob    = 0.3*place_cup_prob;
			plant_prob  = 0.3*place_plant_prob;
			laptop_prob = 0.7*place_laptop_prob;
		}
		else if (obj.type == TYPE_COUNTER && !(obj.flags & RO_FLAG_ADJ_TOP)) { // counter without a microwave
			book_prob   = (placed_book_on_counter ? 0.0 : 0.5); // only place one book per counter
			bottle_prob = 0.25*place_bottle_prob;
			cup_prob    = 0.30*place_cup_prob;
			plant_prob  = 0.10*place_plant_prob;
			laptop_prob = 0.05*place_laptop_prob;
		}
		else {
			continue;
		}
		if (is_library) {book_prob *= 2.5;} // higher probability of books placed in a library
		if (is_kitchen) {cup_prob  *= 2.0;} // higher probability of cups  placed in a kitchen
		room_object_t surface(obj); // deep copy to allow modification and avoid using an invalidated reference
		
		if (obj.shape == SHAPE_CYLIN) { // find max contained XY rectangle (simpler than testing distance to center vs. radius)
			for (unsigned d = 0; d < 2; ++d) {surface.expand_in_dim(d, -0.5*(1.0 - SQRTOFTWOINV)*surface.get_sz_dim(d));}
		}
		if (is_eating_table) { // table in a room for eating, add a plate
			if (place_plate_on_obj(rgen, surface, room_id, tot_light_amt, avoid)) {avoid = objs.back();}
		}
		if (avoid.is_all_zeros() && book_prob > 0.0 && rgen.rand_float() < book_prob) { // place book if it's the first item (no plate)
			placed_book_on_counter |= (obj.type == TYPE_COUNTER);
			place_book_on_obj(rgen, surface, room_id, tot_light_amt, (obj.type != TYPE_TABLE));
			avoid = objs.back();
		}
		if (avoid.is_all_zeros() && obj.type == TYPE_DESK) {
			// if we have no other avoid object, and this is a desk, try to avoid placing an object that overlaps a pen or pencil
			for (unsigned j = i+1; j < objs_end; ++j) {
				room_object_t const &obj2(objs[j]);
				if (obj2.type == TYPE_PEN || obj2.type == TYPE_PENCIL) {avoid = obj2; break;} // we can only use the first one
			}
		}
		if      (bottle_prob > 0.0 && rgen.rand_float() < bottle_prob && place_bottle_on_obj(rgen, surface, room_id, tot_light_amt, avoid)) {}
		else if (cup_prob    > 0.0 && rgen.rand_float() < cup_prob    && place_cup_on_obj   (rgen, surface, room_id, tot_light_amt, avoid)) {}
		else if (laptop_prob > 0.0 && rgen.rand_float() < laptop_prob && place_laptop_on_obj(rgen, surface, room_id, tot_light_amt, avoid, (obj.type != TYPE_TABLE))) {}
		// don't add both a plant and a bottle; don't add plants in the basement
		else if (!is_basement && plant_prob > 0.0 && rgen.rand_float() < plant_prob && place_plant_on_obj(rgen, surface, room_id, tot_light_amt, avoid)) {}
	} // for i
}

void set_light_xy(cube_t &light, point const &center, float light_size, bool light_dim, room_obj_shape light_shape) {
	for (unsigned dim = 0; dim < 2; ++dim) {
		float const sz(((light_shape == SHAPE_CYLIN || light_shape == SHAPE_SPHERE) ? 1.6 : ((bool(dim) == light_dim) ? 2.2 : 1.0))*light_size);
		light.d[dim][0] = center[dim] - sz;
		light.d[dim][1] = center[dim] + sz;
	}
}

bool has_bcube_int_stairs_exp(cube_t const &bcube, vect_stairwell_t const &stairs, float expand, bool stacked_only) {
	cube_t bcube_exp(bcube);
	bcube_exp.expand_by(expand); // expand in all dirs, including z
	
	for (auto s = stairs.begin(); s != stairs.end(); ++s) {
		if (stacked_only && !s->stack_conn) continue; // not stacked stairs
		if (s->intersects(bcube)) return 1;
	}
	return 0;
}

unsigned calc_num_floors_room(room_t const &r, float window_vspacing, float floor_thickness) {
	return (r.is_sec_bldg ? 1 : calc_num_floors(r, window_vspacing, floor_thickness));
}

bool any_cube_contains(cube_t const &cube, vect_cube_t const &cubes) {
	for (cube_t const &c : cubes) {if (c.contains_cube(cube)) return 1;}
	return 0;
}
bool building_t::is_light_placement_valid(cube_t const &light, room_t const &room, float pad) const {
	cube_t light_ext(light);
	light_ext.expand_by_xy(pad);
	if (!room.contains_cube(light_ext)) return 0; // room too small?
	light_ext.z1() = light_ext.z1() = light.z2() + get_fc_thickness(); // shift in between the ceiling and floor so that we can do a cube contains check
	return any_cube_contains(light_ext, interior->fc_occluders);
}
void building_t::try_place_light_on_ceiling(cube_t const &light, room_t const &room, bool room_dim, float pad, bool allow_rot, bool allow_mult, vect_cube_t &lights, rand_gen_t &rgen) {
	assert(has_room_geom());
	if (is_light_placement_valid(light, room, pad)) {lights.push_back(light); return;} // contained = done
	point room_center(room.get_cube_center());
	bool const first_dir(rgen.rand_bool());
	float const window_vspacing(get_window_vspace()), light_width(light.get_sz_dim(!room_dim));
	cube_t light_cand(light);

	if (allow_rot) { // flip aspect ratio
		float const sz_diff(0.5*(light.dx() - light.dy()));
		light_cand.expand_in_dim(0, -sz_diff);
		light_cand.expand_in_dim(1,  sz_diff);
	}
	for (unsigned d = 0; d < 2; ++d) { // see if we can place it by moving on one direction
		for (unsigned n = 0; n < 10; ++n) { // try 10 different shift values
			cube_t cand(light_cand);
			cand.translate_dim(room_dim, ((bool(d) ^ first_dir) ? -1.0 : 1.0)*n*light_width);
			if (!is_light_placement_valid(cand, room, pad)) continue;
			cube_t test_cube(cand);
			test_cube.expand_in_dim(2, 0.4*window_vspacing); // expand to cover nearly an entire floor so that it's guaranteed to overlap a door
			
			// maybe should exclude basement doors, since they don't show as open? but then it would be wrong if I later draw basement doors
			if (is_cube_close_to_doorway(test_cube, room, 0.0, 1, 1)) { // inc_open=1, check_open_dir=1
				cand.z1() += 0.99*cand.dz(); // if light intersects door, move it up into the ceiling rather than letting it hang down into the room
			}
			lights.push_back(cand);
			break;
		} // for n
		if (!allow_mult) break;
	} // for d
}

colorRGBA get_light_color_temp(float t) {
	// 0.0: 1.0 1.0 0.5
	// 0.5: 1.0 1.0 1.0
	// 1.0: 0.5 0.5 1.0
	if (t > 0.5) {return colorRGBA(1.5-t, 1.5-t, 1.0  );} // high temp blue spectrum
	else         {return colorRGBA(1.0,   1.0,   t+0.5);} // low temp yellow spectrum
}
colorRGBA get_light_color_temp_range(float tmin, float tmax, rand_gen_t &rgen) {
	return get_light_color_temp(((tmin == tmax) ? tmin : rgen.rand_uniform(tmin, tmax)));
}

// Note: these three floats can be calculated from get_window_vspace(), but it's easier to change the constants if we just pass them in
void building_t::gen_room_details(rand_gen_t &rgen, unsigned building_ix) {

	assert(interior);
	if (interior->room_geom) return; // already generated?
	setup_bldg_obj_types(); // initialize object types if not already done
	//highres_timer_t timer("Gen Room Details");
	// Note: people move from room to room, so using their current positions for room object generation is both nondeterministic and unnecessary
	vect_cube_t blockers, valid_lights; // blockers are used for fireplaces
	interior->room_geom.reset(new building_room_geom_t(bcube.get_llc()));
	vect_room_object_t &objs(interior->room_geom->objs);
	vector<room_t> &rooms(interior->rooms);
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness);
	float const light_thick(0.025*window_vspacing), def_light_size(0.1*window_vspacing);
	interior->room_geom->obj_scale = window_vspacing; // used to scale room object textures
	unsigned tot_num_rooms(0), num_bathrooms(0);
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {tot_num_rooms += calc_num_floors_room(*r, window_vspacing, floor_thickness);}
	objs.reserve(tot_num_rooms); // placeholder - there will be more than this many
	float const extra_bathroom_prob((is_house ? 2.0 : 1.0)*0.02*min((int(tot_num_rooms) - 4), 20));
	room_obj_shape const light_shape(is_house ? SHAPE_CYLIN : SHAPE_CUBE);
	unsigned cand_bathroom(rooms.size()); // start at an invalid value
	unsigned added_kitchen_mask(0); // per-floor
	unsigned added_bathroom_objs_mask(0);
	bool added_bedroom(0), added_living(0), added_library(0), added_dining(0), added_laundry(0), added_basement_utility(0);
	light_ix_assign_t light_ix_assign;
	interior->create_fc_occluders(); // not really part of room geom, but needed for generating and drawing room geom, so we create them here

	if (rooms.size() > 1) { // choose best room assignments for required rooms; if a single room, skip this step
		float min_score(0.0);

		// Note: assigning cand_bathroom when has_pri_hall() is not strictly necessary, but may help add a bathroom to an upper stacked part
		for (auto r = rooms.begin(); r != rooms.end(); ++r) {
			if (r->is_sec_bldg) continue; // garage/shed excluded - not a normal room
			if (has_basement() && r->part_id == (int)basement_part_ix) continue; // skip the basement
			unsigned const num_floors(calc_num_floors_room(*r, window_vspacing, floor_thickness));

			// find best bathroom with no hard size constraints;
			// use the top floor for the test since it's less restrictive than the ground floor; will be checked per-floor later
			if (can_be_bedroom_or_bathroom(*r, (num_floors-1), 0)) { // skip_conn_check=0
				if (has_chimney == 2 && num_floors == 1) { // can't be a bathroom if there's a fireplace
					cube_t test_cube(*r);
					test_cube.expand_by_xy(floor_thickness);
					if (test_cube.intersects(get_fireplace())) continue;
				}
				float score(r->dx() + r->dy()); // starts as half the perimeter
				score *= (1.0 + 10.0*(max(count_num_int_doors(*r), 1U) - 1U)); // multiply by a large value if there are mult doors so we only choose this if there are no alternatives
				if (min_score == 0.0 || score < min_score) {cand_bathroom = (r - rooms.begin()); min_score = score;}
			}
		} // for r
	}
	for (auto r = rooms.begin(); r != rooms.end(); ++r) {
		bool const is_basement(has_basement() && r->part_id == (int)basement_part_ix);
		float const light_amt(is_basement ? 0.0f : window_vspacing*r->get_light_amt()); // exterior light: multiply perimeter/area by window spacing to make unitless; none for basement rooms
		float const floor_height(r->is_sec_bldg ? r->dz() : window_vspacing); // secondary buildings are always one floor
		unsigned const num_floors(calc_num_floors_room(*r, floor_height, floor_thickness)), room_id(r - rooms.begin());
		point room_center(r->get_cube_center());

		// determine light pos and size for this stack of rooms
		bool const room_dim(r->dx() < r->dy()); // longer room dim
		bool const must_be_bathroom(room_id == cand_bathroom && num_bathrooms == 0); // cand bathroom, and bathroom not already placed
		bool const is_parking_garage(r->get_room_type(0) == RTYPE_PARKING); // all floors should be parking garage
		float light_size(def_light_size); // default size for houses
		unsigned const room_objs_start(objs.size());
		unsigned nx(1), ny(1); // number of lights in X and Y for this room

		if (r->is_office) { // more lights for large offices; parking garages are handled later
			nx = max(1U, unsigned(0.5*r->dx()/window_vspacing));
			ny = max(1U, unsigned(0.5*r->dy()/window_vspacing));
		}
		if (r->is_sec_bldg) {
			if    (has_garage) {r->assign_all_to(RTYPE_GARAGE);}
			else if (has_shed) {r->assign_all_to(RTYPE_SHED);}
		}
		if (r->is_office) { // light size varies by office size
			float const room_size(r->dx() + r->dy()); // normalized to office size
			light_size = max(0.015f*room_size, 0.67f*def_light_size);
		}
		if (r->is_hallway) { // light size varies by hallway size
			float const room_size(min(r->dx(), r->dy())); // normalized to hallway width
			light_size = max(0.06f*room_size, 0.67f*def_light_size);
		}
		float const light_val(22.0*light_size);
		r->light_intensity = light_val*light_val/r->get_area_xy(); // average for room, unitless; light surface area divided by room surface area with some fudge constant
		cube_t light;
		set_light_xy(light, room_center, light_size, room_dim, light_shape);
		bool added_bathroom(0);
		float z(r->z1());
		r->interior = (is_basement || get_part_for_room(*r).contains_cube_xy_no_adj(*r)); // AKA windowless
		// make chair colors consistent for each part by using a few variables for a hash
		colorRGBA chair_colors[12] = {WHITE, WHITE, GRAY, DK_GRAY, LT_GRAY, BLUE, DK_BLUE, LT_BLUE, YELLOW, RED, DK_GREEN, LT_BROWN};
		colorRGBA chair_color(chair_colors[(13*r->part_id + 123*tot_num_rooms + 617*mat_ix + 1367*num_floors) % 12]);
		light_ix_assign.next_room();
		// select light color for this room
		colorRGBA color;
		if      (is_house)          {color = get_light_color_temp(0.4);} // house - yellowish
		else if (is_parking_garage) {color = get_light_color_temp_range(0.2, 0.5, rgen);} // parking garage - yellow-white
		else if (r->is_office)      {color = get_light_color_temp(0.6);} // office - blueish
		else if (r->is_hallway)     {color = get_light_color_temp(0.6);} // hallway - blueish
		else                        {color = get_light_color_temp(0.5);} // small office - white

		// place objects on each floor for this room
		for (unsigned f = 0; f < num_floors; ++f, z += floor_height) {
			room_center.z = z + fc_thick; // floor height
			// top floor may have stairs connecting to upper stack
			bool const top_floor(f+1 == num_floors);
			bool const has_stairs_this_floor(r->has_stairs_on_floor(f));
			bool is_lit(0), has_light(1), light_dim(room_dim), has_stairs(has_stairs_this_floor), top_of_stairs(has_stairs && top_floor);
			float light_delta_z(0.0);

			if (is_parking_garage) { // parking garage; added first because this sets the number of lights
				r->interior = 1;
				add_parking_garage_objs(rgen, *r, room_center.z, room_id, f, num_floors, nx, ny, light_delta_z);
				for (auto i = objs.begin() + room_objs_start; i != objs.end(); ++i) {i->flags |= RO_FLAG_INTERIOR;}
			}
			if ((!has_stairs && (f == 0 || top_floor) && interior->stairwells.size() > 1) || top_of_stairs) { // should this be outside the loop?
				// check for stairwells connecting stacked parts (is this still needed?); check for roof access stairs and set top_of_stairs=0
				for (auto s = interior->stairwells.begin(); s != interior->stairwells.end(); ++s) {
					if (!r->contains_cube_xy(*s)) continue; // stairs not in this room
					// Note: here we adjust stairs zval by floor_thickness to include stairs in the floor but not in the room above
					if (s->z1() + floor_thickness > r->z2()) continue; // stairs above the room
					if (s->z2() + floor_thickness < r->z1()) continue; // stairs below the room
					if (s->roof_access) {top_of_stairs = 0;}
					has_stairs = 1;
				} // for s
			}
			int light_obj_ix(-1);
			unsigned num_lights(r->num_lights);
			float const light_z2(z + floor_height - fc_thick + light_delta_z);
			// 50% of lights are on, 75% for top of stairs, 100% for hallways, 100% for parking garages
			is_lit = (r->is_hallway || is_parking_garage || ((rgen.rand() & (top_of_stairs ? 3 : 1)) != 0));

			if (!is_lit) { // check people and set is_lit if anyone is in this floor of this room
				for (person_t const &p : interior->people) {
					cube_t const bc(p.get_bcube());
					if (!bc.intersects_xy(*r)) continue; // person not in this room
					if (bc.z2() < light_z2 && bc.z1() + floor_height > light_z2) {is_lit = 1; break;} // on this floor
				}
			}
			unsigned flags(RO_FLAG_NOCOLL); // no collision detection with lights
			if (is_lit)     {flags |= RO_FLAG_LIT | RO_FLAG_EMISSIVE;}
			if (has_stairs) {flags |= RO_FLAG_RSTAIRS;}
			// add a light to the ceiling of this room if there's space (always for top of stairs);
			set_cube_zvals(light, (light_z2 - light_thick), light_z2);
			valid_lights.clear();

			if (r->is_hallway && num_lights > 1) { // hallway: place a light on each side (of the stairs if they exist), and also between stairs and elevator if there are both
				if (r->has_elevator && r->has_stairs == 255) {num_lights = 3;} // main hallway with elevator + stairs on all floors: we really should have 3 lights in this case
				float const offset(((num_lights == 3) ? 0.3 : 0.2)*r->get_sz_dim(light_dim)); // closer to the ends in the 3 lights case

				for (unsigned d = 0; d < num_lights; ++d) {
					float const delta((d == 2) ? 0.0 : (d ? -1.0 : 1.0)*offset); // last light is in the center
					cube_t hall_light(light);
					hall_light.translate_dim(light_dim, delta);
					try_place_light_on_ceiling(hall_light, *r, room_dim, fc_thick, 0, 0, valid_lights, rgen); // allow_rot=0, allow_mult=0
				}
			}
			else if (nx > 1 || ny > 1) { // office or parking garage with multiple lights
				float const dx(r->dx()), dy(r->dy()), xstep(dx/nx), ystep(dy/ny);
				vector3d const shrink(0.5*light.dx()*sqrt((nx - 1)/nx), 0.5*light.dy()*sqrt((ny - 1)/ny), 0.0);

				for (unsigned y = 0; y < ny; ++y) {
					for (unsigned x = 0; x < nx; ++x) {
						cube_t cur_light(light);
						cur_light.expand_by_xy(-shrink);
						cur_light.translate(point((-0.5f*dx + (x + 0.5)*xstep), (-0.5f*dy + (y + 0.5)*ystep), 0.0));
						try_place_light_on_ceiling(cur_light, *r, room_dim, fc_thick, 0, 0, valid_lights, rgen); // allow_rot=0, allow_mult=0
					}
				} // for y
			}
			else { // normal room with a single light
				try_place_light_on_ceiling(light, *r, room_dim, fc_thick, 1, 1, valid_lights, rgen); // allow_rot=1, allow_mult=1
				if (!valid_lights.empty()) {light_obj_ix = objs.size();} // this will be the index of the light to be added later
			}
			rand_gen_t rgen_lights(rgen); // copy state so that we don't modify rgen

			for (cube_t const &l : valid_lights) {
				objs.emplace_back(l, TYPE_LIGHT, room_id, (light.dx() < light.dy()), 0, flags, light_amt, light_shape, color); // reclaculate dim; dir=0 (unused)
				objs.back().obj_id = light_ix_assign.get_ix_for_light(l);
				if (is_parking_garage && (rgen_lights.rand()%50 == 13)) {objs.back().flags |= RO_FLAG_BROKEN;} // 2% chance of a flickering light
			}
			if (is_lit) {r->lit_by_floor |= (1ULL << (f&63));} // flag this floor as being lit (for up to 64 floors)
			if (is_parking_garage) continue; // generated above, done; no outlets or light switches
			float tot_light_amt(light_amt); // unitless, somewhere around 1.0
			if (is_lit) {tot_light_amt += r->light_intensity;}
			bool const is_ground_floor(f == 0 && !is_basement), is_garage_or_shed(r->is_garage_or_shed(f));
			rgen.rand_mix();

			if (r->no_geom || is_garage_or_shed) {
				if (is_garage_or_shed) {
					if (r->get_room_type(0) == RTYPE_GARAGE) {
						room_center.z = add_flooring(*r, room_center.z, room_id, tot_light_amt);
						add_garage_objs(rgen, *r, room_center.z, room_id, tot_light_amt);
					}
					// is there enough clearance between shelves and a car parked in the garage? there seems to be in all the cases I've seen
					add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs.size(), is_basement);
				}
				add_outlets_to_room(rgen, *r, room_center.z, room_id, objs.size(), is_ground_floor, is_basement);
				if (has_light) {add_light_switches_to_room(rgen, *r, room_center.z, room_id, objs.size(), is_ground_floor);} // shed, garage, or hallway

				if (is_house && r->is_hallway) { // allow pictures, rugs, and light switches in the hallways of houses
					hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs.size(), is_basement);
					if (rgen.rand_bool()) {add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs.size());} // 50% of the time; not all rugs will be placed
				}
				if (!is_house && r->is_hallway && f == 0 && *r == pri_hall) { // first floor primary hallway, make it the lobby
					add_pri_hall_objs(rgen, *r, room_center.z, room_id, tot_light_amt);
					r->assign_to(RTYPE_LOBBY, f);
				}
				continue; // no other geometry for this room
			}
			//if (has_stairs && !pri_hall.is_all_zeros()) continue; // no other geometry in office building base part rooms that have stairs
			unsigned const objs_start(objs.size()), floor_mask(1<<f);
			bool added_tc(0), added_desk(0), added_obj(0), can_place_onto(0), is_bathroom(0), is_bedroom(0), is_kitchen(0), is_living(0), is_dining(0), is_storage(0), no_whiteboard(0);
			unsigned num_chairs(0);

			// place room objects
			bool const allow_br(!is_house || must_be_bathroom || f > 0 || num_floors == 1 || (rgen.rand_float() < 0.33f*(added_living + (added_kitchen_mask&1) + 1))); // bed/bath
			bool is_office_bathroom(is_room_office_bathroom(*r, room_center.z, f)), has_fireplace(0);
			blockers.clear(); // clear for this new room
			
			if (has_chimney == 2 && !is_basement && f == 0) { // handle fireplaces on the first floor
				has_fireplace = maybe_add_fireplace_to_room(*r, blockers, room_center.z, room_id, tot_light_amt);
			}
			if (is_office_bathroom) { // bathroom is already assigned
				added_obj = is_bathroom = added_bathroom = no_whiteboard =
					add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement, added_bathroom_objs_mask); // add bathroom
			}
			else if (!is_house && f == 0 && r->get_room_type(f) == RTYPE_UTILITY) { // office building utility room; currently first floor only
				added_obj = add_office_utility_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
			}
			// bedroom or bathroom case; need to check first floor even if must_be_bathroom
			if (!added_obj && allow_br && can_be_bedroom_or_bathroom(*r, f)) {
				// place a bedroom 75% of the time unless this must be a bathroom; if we got to the second floor and haven't placed a bedroom, always place it; houses only
				if (is_house && !must_be_bathroom && !is_basement && ((f > 0 && !added_bedroom) || rgen.rand_float() < 0.75)) {
					added_obj = added_bedroom = is_bedroom =
						add_bedroom_objs(rgen, *r, blockers, room_center.z, room_id, f, tot_light_amt, objs_start, is_lit, is_basement, light_ix_assign);
					if (is_bedroom) {r->assign_to(RTYPE_BED, f);}
					// Note: can't really mark room type as bedroom because it varies per floor; for example, there may be a bedroom over a living room connected to an exterior door
				}
				if (!added_obj && !has_fireplace && (must_be_bathroom || (can_be_bathroom(*r) && (num_bathrooms == 0 || rgen.rand_float() < extra_bathroom_prob)))) {
					// bathrooms can be in both houses and office buildings
					added_obj = is_bathroom = added_bathroom =
						add_bathroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, f, is_basement, added_bathroom_objs_mask); // add bathroom
					if (is_bathroom) {r->assign_to(RTYPE_BATH, f);}
				}
			}
			if (!added_obj && r->is_office) { // add cubicles if this is a large office
				added_obj = no_whiteboard = create_office_cubicles(rgen, *r, room_center.z, room_id, tot_light_amt);
			}
			if (!added_obj && rgen.rand_float() < (is_basement ? 0.4 : (r->is_office ? 0.6 : (is_house ? 0.95 : 0.5)))) {
				// place a table and maybe some chairs near the center of the room if it's not a hallway;
				// 60% of the time for offices, 95% of the time for houses, and 50% for other buildings
				unsigned const num_tcs(add_table_and_chairs(rgen, *r, blockers, room_id, room_center, chair_color, 0.1, tot_light_amt));
				if (num_tcs > 0) {added_tc = added_obj = can_place_onto = 1; num_chairs = num_tcs - 1;}
				// on ground floor, try to make this a kitchen; not all houses will have a kitchen with this logic - maybe we need fewer bedrooms?
				if (!(added_kitchen_mask & floor_mask) && (!is_house || f == 0) && !is_basement) { // office buildings can also have kitchens, even on non-ground floors
					if (added_tc || (is_house && (r+1) == rooms.end())) { // make it a kitchen if it's the last room in a house, even if there's no table
						if (add_kitchen_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, added_living)) {
							r->assign_to(RTYPE_KITCHEN, f);
							added_kitchen_mask |= floor_mask;
							is_kitchen = added_obj = 1;
						}
					}
				}
			}
			if (!added_obj && (is_basement || (r->is_office && r->interior && f == 0 /*&& r->z1() == ground_floor_z1*/)) && rgen.rand_bool()) {
				// if we haven't added any objects yet, and this room is an interior office on the first floor or basement, make it a storage room 50% of the time
				added_obj = no_whiteboard = is_storage = add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
				if (added_obj) {r->assign_to(RTYPE_STORAGE, f);}
			}
			if (!added_obj && (!is_basement || rgen.rand_bool())) { // try to place a desk if there's no table, bed, etc.
				added_obj = can_place_onto = added_desk = add_desk_to_room(rgen, *r, blockers, chair_color, room_center.z, room_id, f, tot_light_amt, is_basement);
				if (added_obj && !has_stairs_this_floor) {r->assign_to((is_house ? (room_type)RTYPE_STUDY : (room_type)RTYPE_OFFICE), f);} // or other room type - may overwrite below
			}
			if (is_house && (added_tc || added_desk) && !is_kitchen && f == 0) { // don't add second living room unless we added a kitchen and have enough rooms
				if ((!added_living && !r->has_center_stairs && rooms.size() >= 8 && (added_kitchen_mask || rgen.rand_bool())) || is_room_adjacent_to_ext_door(*r, 1)) { // front_door_only=1
					// add a living room on the ground floor if it has a table or desk but isn't a kitchen
					added_living = is_living = add_livingroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
					if (is_living) {r->assign_to(RTYPE_LIVING, f);}
				}
			}
			if (is_house && added_tc && num_chairs > 0 && !is_living && !is_kitchen) { // room with table and chair that's not a kitchen
				if (f == 0 && !is_basement) { // dining room, must be on the first floor
					if (light_obj_ix >= 0) { // handle dining room light (assume there is only one): extend downward and make it a sphere
						assert((unsigned)light_obj_ix < objs.size());
						room_object_t &light(objs[light_obj_ix]);
						light.shape = SHAPE_SPHERE;
						light.z2() += 0.5f*light.dz();
						light.z1() -= 0.22f*(light.dx() + light.dy());
					}
					if (!added_dining) {add_diningroom_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);} // only one room is the primary dining room
					r->assign_to(RTYPE_DINING, f);
					is_dining = added_dining = 1;
				}
				else if (!added_library && add_library_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement)) { // add library, at most one
					r->assign_to(RTYPE_LIBRARY, f);
					added_library = 1;
				}
			}
			if (!is_house && r->is_office && !no_whiteboard && (rgen.rand() % (pri_hall.is_all_zeros() ? 30U : max(50U, (unsigned)interior->rooms.size()))) == 0) {
				// office, no cubicles or bathroom - try to make it a library (in rare cases)
				if (add_library_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement)) {r->assign_to(RTYPE_LIBRARY, f);}
			}
			if (can_place_onto) { // an object was placed (table, desk, counter, etc.), maybe add a book or bottle on top of it
				place_objects_onto_surfaces(rgen, *r, room_id, tot_light_amt, objs_start, f, is_basement);
			}
			if (is_house) { // place house-specific items
				if (!is_bathroom && !is_kitchen && rgen.rand_float() < (is_basement ? 0.25 : 0.8)) { // place bookcase 80% of the time, but not in bathrooms or kitchens
					rand_gen_t rgen2(rgen); // copy so that rgen isn't updated in the call below
					add_bookcase_to_room(rgen2, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
				}
				if (!has_stairs && (rgen.rand()&3) <= (added_tc ? 0 : 2) && !is_kitchen) { // maybe add a rug, 25% of the time if there's a table and 75% of the time otherwise
					add_rug_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				}
			}
			bool const room_type_was_not_set(r->get_room_type(f) == RTYPE_NOTSET);

			if (room_type_was_not_set) { // attempt to assign it with an optional room type
				if (!is_basement && f == 0 && is_room_adjacent_to_ext_door(*r)) { // entryway/lobby if on first floor, has exterior door, and unassigned
					r->assign_to((is_house ? (room_type)RTYPE_ENTRY : (room_type)RTYPE_LOBBY), f); // office building lobby can have a whiteboard - is that okay?
				}
				else if (!is_house) {r->assign_to(RTYPE_OFFICE, f);} // any unset room in an office building is an office
				// else house
				else if (has_stairs && !is_basement) {} // will be marked as RTYPE_STAIRS below
				else if ((!added_obj || is_basement) && f == 0 && !added_laundry && !has_fireplace &&
					add_laundry_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, added_bathroom_objs_mask))
				{
					r->assign_to(RTYPE_LAUNDRY, f);
					added_laundry = 1;
				}
				else if (!added_obj && !has_fireplace) { // make it a storage room until we add some other room type that it can be
					add_storage_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement);
					r->assign_to(RTYPE_STORAGE, f);
					is_storage = 1; // mark it as a storage room whether or not we've added anything to it
				}
				else if (is_basement) {r->assign_to(RTYPE_CARD, f);} // basement card room
				else { // unassigned room of house on upper floor with added object/table
					// this case is relatively rare, and we've already added a table, so it's too late to make this a bedroom/bathroom if can_be_bedroom_or_bathroom(*r, f)
					r->assign_to((rgen.rand_bool() ? (room_type)RTYPE_PLAY : (room_type)RTYPE_ART), f); // play room or art room
				}
			}
			if (is_house && is_basement && !added_basement_utility && !has_stairs && (is_storage || room_type_was_not_set) && rgen.rand_bool()) {
				// basement laundry, storage, or card room; should this be placed before adding boxes to the floor of storage rooms?
				added_basement_utility = add_basement_utility_objs(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start);
				if (added_basement_utility) {r->assign_to(RTYPE_UTILITY, f);}
			}
			if (!is_bathroom && !is_bedroom && !is_kitchen && !is_storage && !is_basement) { // add potted plants to some room types
				// 0-2 for living/dining rooms, 50% chance for houses, 25% (first floor) / 10% (other floors) chance for offices
				unsigned const num(is_house ? (rgen.rand() % ((is_living || is_dining) ? 3 : 2)) : ((rgen.rand()%((f == 0) ? 4 : 10)) == 0));
				if (num > 0) {add_plants_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, num);}
			}
			add_outlets_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor, is_basement);
			if (has_light) {add_light_switches_to_room(rgen, *r, room_center.z, room_id, objs_start, is_ground_floor);} // add a light switch if this room has a light
			// pictures and whiteboards must not be placed behind anything, excluding trashcans; so we add them here
			bool const can_hang((is_house || !(is_bathroom || is_kitchen || no_whiteboard)) && !is_storage); // no whiteboards in office bathrooms or kitchens
			bool const was_hung(can_hang && hang_pictures_in_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, is_basement));

			if (is_bathroom || is_kitchen || rgen.rand_float() < 0.8) { // 80% of the time, always in bathrooms and kitchens
				add_trashcan_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, (was_hung && !is_house)); // no trashcans on same wall as office whiteboard
			}
			if (is_house && !(is_bathroom || is_kitchen || is_storage) && rgen.rand_float() < ((f > 0) ? 0.15 : 0.25)) {
				unsigned const max_num(is_bedroom ? 1 : 2);
				add_boxes_to_room(rgen, *r, room_center.z, room_id, tot_light_amt, objs_start, max_num); // place boxes in this room
			}
			if (has_stairs_this_floor && r->get_room_type(f) == RTYPE_NOTSET) {r->assign_to(RTYPE_STAIRS, f);} // default to stairs if not set above
		} // for f (floor)
		if (added_bathroom) {++num_bathrooms;}

		if (r->interior) { // tag objects as interior if room is interior
			for (auto i = objs.begin() + room_objs_start; i != objs.end(); ++i) {i->flags |= RO_FLAG_INTERIOR;}
		}
	} // for r (room)
	if (is_rotated()) {} // skip for rotated buildings, since toilets, etc. may not be placed
	else if (num_bathrooms == 0) { // can happen, but very rare
		cout << "no bathroom in building " << bcube.xc() << " " << bcube.yc() << endl;
		if (cand_bathroom < rooms.size()) {cout << "cand bathroom was at " << rooms[cand_bathroom].str() << endl;}
	}
	else {
		if (!(added_bathroom_objs_mask & PLACED_TOILET)) {cout << "no toilet in building " << bcube.xc() << " " << bcube.yc() << endl;}
		if (!(added_bathroom_objs_mask & PLACED_SINK  )) {cout << "no sink in building "   << bcube.xc() << " " << bcube.yc() << endl;}
		//if (is_house && !(added_bathroom_objs_mask & (PLACED_TUB | PLACED_SHOWER))) {cout << "no bathtub or shower in building " << bcube.xc() << " " << bcube.yc() << endl;} // common
	}
	if (is_house && has_basement()) {add_basement_electrical_house(rgen);}
	maybe_add_fire_escape(rgen);
	add_extra_obj_slots(); // needed to handle balls taken from one building and brought to another
	add_stairs_and_elevators(rgen); // the room objects - stairs and elevators have already been placed within a room
	add_exterior_door_signs (rgen);
	objs.shrink_to_fit(); // Note: currently up to around 15K objs max for large office buildings
	interior->room_geom->light_bcubes.resize(light_ix_assign.get_next_ix()); // allocate but don't fill un until needed
	// randomly vary wood color for this building
	colorRGBA &wood_color(interior->room_geom->wood_color);
	float const luminance(rgen.rand_uniform(0.4, 1.6));
	for (unsigned i = 0; i < 3; ++i) {wood_color[i] = luminance*WOOD_COLOR[i]*rgen.rand_uniform(0.9, 1.1);}
	wood_color.set_valid_color();
	max_eq(wood_color.R, max(wood_color.G, wood_color.B)); // make sure wood isn't blue or green tinted
}

void building_t::maybe_add_fire_escape(rand_gen_t &rgen) {
	if (!is_house) return; // houses only for now
	// our hard-coded fire escape model is designed for a 5 story building; but the max number of floors for a 'house' is 5-6 anyway, which makes them relatively rare
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_FESCAPE)) return;
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fe_height(4.25*window_vspacing);

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) {
		unsigned const num_floors(calc_num_floors(*p, window_vspacing, floor_thickness));
		if (num_floors != 5 && num_floors != 6) continue; // not 5-6 stories
		unsigned const pref_dim_dir(rgen.rand() & 3);
		// it's uncommon to get here, so we only check if the model size here
		vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_FESCAPE)); // D, W, H
		float const fe_hwidth(0.5*fe_height*sz.y/sz.z), fe_depth(fe_height*sz.x/sz.z);

		for (unsigned d = 0; d < 4; ++d) {
			unsigned const dd((d + pref_dim_dir) & 3);
			bool const dim(dd >> 1), dir(dd & 1);
			if (p->d[dim][dir] != bcube.d[dim][dir]) continue; // not on the building bcube - could intersect another part, porch, etc.
			if (p->get_sz_dim(!dim) < 3.0*fe_hwidth) continue; // wall is too narrow
			cube_t fe_bc;
			set_cube_zvals(fe_bc, p->z1(), (p->z1() + fe_height));
			set_wall_width(fe_bc, rgen.rand_uniform((p->d[!dim][0] + 1.2*fe_hwidth), (p->d[!dim][1] - 1.2*fe_hwidth)), fe_hwidth, !dim);
			fe_bc.d[dim][0] = fe_bc.d[dim][1] = p->d[dim][dir];
			fe_bc.d[dim][dir] += (dir ? 1.0 : -1.0)*fe_depth;
			if (has_bcube_int_no_adj(fe_bc, parts)) continue; // check for intersection with other parts, in particular the chimney and fireplace
			if (has_driveway() && fe_bc.intersects_xy(driveway)) continue; // skip if intersects driveway or garage
			interior->room_geom->objs.emplace_back(fe_bc, TYPE_FESCAPE, 0, dim, dir, 0, 1.0, SHAPE_CUBE, BLACK); // room_id=0
			break; // success/done
		} // for d
	} // for p
}

void building_t::add_extra_obj_slots() {
	assert(has_room_geom());
	vect_room_object_t &objs(interior->room_geom->objs);
	if (objs.empty()) return; // if there are no objects (empty building), don't allocate any extra slots
	unsigned num_slots(0);
	for (auto i = objs.begin(); i != objs.end(); ++i) {num_slots += (i->type == TYPE_BLOCKER);}
	if (num_slots >= 10) return;
	// make sure there are at least 10 blockers that will create free slots when adding dynamic objects
	float const v(0.01*get_wall_thickness()); // some tiny number
	point const llc(bcube.get_llc());
	cube_t const c(llc, llc+vector3d(v, v, v));
	for (unsigned n = num_slots; n < 20; ++n) {objs.emplace_back(c, TYPE_BLOCKER, 0, 0, 0, (RO_FLAG_INVIS | RO_FLAG_NOCOLL));}
}

void building_t::add_wall_and_door_trim_if_needed() {
	if (!is_cube()) return; // not yet supported for non-cube buildings
	if (!interior->room_geom->trim_objs.empty()) return; // trim already generated
	add_wall_and_door_trim();
	interior->room_geom->trim_objs.shrink_to_fit();
}

void building_t::add_wall_and_door_trim() { // and window trim
	//highres_timer_t timer("Add Wall And Door Trim");
	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), fc_thick(0.5*floor_thickness), wall_thickness(get_wall_thickness());
	float const trim_height(get_trim_height()), trim_thickness(get_trim_thickness()), expand_val(2.0*trim_thickness);
	float const door_trim_exp(2.0*trim_thickness + 0.5*wall_thickness), door_trim_width(0.5*wall_thickness), floor_to_ceil_height(window_vspacing - floor_thickness);
	float const trim_toler(0.1*trim_thickness); // required to handle wall intersections that were calculated with FP math and may misalign due to FP rounding error
	float const ext_wall_toler(0.01*trim_thickness); // required to prevent z-fighting when AA is disabled
	unsigned const flags(RO_FLAG_NOCOLL);
	// ceiling trim disabled for large office buildings with outside corners because there's a lot of trim to add, and outside corners don't join correctly;
	// ceiling trim also disabled for non-houses (all office buildings), because it doesn't really work with acoustic paneling
	bool const has_outside_corners(!is_house && !pri_hall.is_all_zeros()), has_ceil_trim(!has_outside_corners && is_house);
	colorRGBA const &trim_color(is_house ? WHITE : DK_GRAY);
	vect_room_object_t &objs(interior->room_geom->trim_objs);
	vect_cube_t trim_cubes;

	for (auto d = interior->door_stacks.begin(); d != interior->door_stacks.end(); ++d) { // vertical strips on each side + strip on top of interior doors
		if (d->on_stairs) continue; // no frame for stairs door, skip
		cube_t trim(*d);
		trim.expand_in_dim(d->dim, door_trim_exp);

		for (unsigned side = 0; side < 2; ++side) { // left/right of door
			trim.d[!d->dim][0] = d->d[!d->dim][side] - (side ? trim_thickness : door_trim_width);
			trim.d[!d->dim][1] = d->d[!d->dim][side] + (side ? door_trim_width : trim_thickness);
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, d->dim, side, (flags | RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP), 1.0, SHAPE_TALL, trim_color); // abuse tall flag
		}
		// add trim at top of door
		unsigned const num_floors(calc_num_floors(*d, window_vspacing, floor_thickness));
		float z(d->z1() + floor_to_ceil_height);
		trim.d[!d->dim][0] = d->d[!d->dim][0] + trim_thickness;
		trim.d[!d->dim][1] = d->d[!d->dim][1] - trim_thickness;

		for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
			set_cube_zvals(trim, z-trim_thickness, z); // z2=ceil height
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, d->dim, 0, (flags | RO_FLAG_ADJ_TOP), 1.0, SHAPE_SHORT, trim_color);
		}
	} // for d
	for (auto d = doors.begin(); d != doors.end(); ++d) { // exterior doors
		if (d->type == tquad_with_ix_t::TYPE_RDOOR) { // roof access door
			continue; // this requires a completely different approach to trim and has not yet been implemented
		}
		bool const garage_door(d->type == tquad_with_ix_t::TYPE_GDOOR);
		cube_t door(d->get_bcube());
		bool const dim(door.dy() < door.dx());
		//door.expand_in_dim(!dim, -0.1*door_trim_width); // shrink slightly so that the edge of the wall is contained in the trim
		cube_t trim(door);
		trim.expand_in_dim(dim, door_trim_exp);
		bool dir(0);
		unsigned ext_flags(flags);
		colorRGBA const &ext_trim_color(garage_door ? WHITE : door_color); // garage doors are always white
		float const trim_width(((garage_door && has_int_garage) ? 1.5 : 1.0)*door_trim_width); // interior garage doors have thicker trim

		for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
			if (!i->intersects_no_adj(trim)) continue;
			trim.intersect_with_cube_xy(*i); // clip to containing part
			dir = (i->get_center_dim(dim) < trim.get_center_dim(dim));
			trim.d[dim][dir] -= (dir ? -1.0 : 1.0)*0.025*window_vspacing; // move to to same offset for door
			ext_flags |= (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO);
			break;
		}
		for (unsigned side = 0; side < 2; ++side) { // left/right of door
			trim.d[!dim][0] = door.d[!dim][side] - (side ? trim_width : 0.0);
			trim.d[!dim][1] = door.d[!dim][side] + (side ? 0.0 : trim_width);
			objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, side, (ext_flags | RO_FLAG_ADJ_BOT | RO_FLAG_ADJ_TOP), 1.0, SHAPE_TALL, ext_trim_color); // abuse tall flag
		}
		// add trim at bottom of door for threshold
		trim.d[!dim][0] = door.d[!dim][0];
		trim.d[!dim][1] = door.d[!dim][1];
		set_cube_zvals(trim, door.z1()+fc_thick, door.z1()+fc_thick+2.0*trim_thickness); // floor height
		objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, (ext_flags | RO_FLAG_ADJ_BOT), 1.0, SHAPE_SHORT, ext_trim_color);

		if (d->type == tquad_with_ix_t::TYPE_HDOOR || d->is_building_door() || garage_door) { // add trim at top of exterior door, houses and office buildings
			set_cube_zvals(trim, door.z2()-0.03*door.dz(), door.z2()); // ends at top of door texture; see logic in clip_door_to_interior()
		}
		if (d->is_building_door()) { // different logic for building doors
			ext_flags = flags; // unlike hdoors, need to draw the back face to hide the gap betweeen ceiling and floor above
			trim.d[dim][dir] += (dir ? -1.0 : 1.0)*0.005*window_vspacing; // minor shift back toward building to prevent z-fighting
		}
		objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_SHORT, ext_trim_color); // top of door
	} // for d
	for (unsigned dim = 0; dim < 2; ++dim) { // add horizontal strips along each wall at each floor/ceiling
		for (auto w = interior->walls[dim].begin(); w != interior->walls[dim].end(); ++w) {
			cube_t trim(*w);
			trim.expand_in_dim(dim, trim_thickness);

			if (has_outside_corners) { // handle outside corners of office building hallway intersections
				for (auto W = interior->walls[!dim].begin(); W != interior->walls[!dim].end(); ++W) { // check walls in other dim for an outside corner
					for (unsigned d = 0; d < 2; ++d) {
						if (W->z1() > w->z2() || W->z2() < w->z1()) continue; // no z overlap, wrong stack
						if (W->d[!dim][0] > w->d[!dim][d]+trim_toler || W->d[!dim][1] < w->d[!dim][d]-trim_toler) continue; // not adjacent/overlapping
						if (W->d[ dim][0] < w->d[ dim][0]-trim_toler && W->d[ dim][1] > w->d[ dim][1]+trim_toler) continue; // skip T junctions
						trim.d[!dim][d] = W->d[!dim][d] + (d ? 1.0 : -1.0)*trim_thickness; // expand to cover gap at outside corners of hallway walls
					}
				} // for W
			}
			unsigned const num_floors(calc_num_floors(*w, window_vspacing, floor_thickness));
			// snap to the nearest floor to handle short walls due to cut out stairs
			float const ground_wall_z1(bcube.z1() + fc_thick);
			float z(ground_wall_z1 + window_vspacing*round_fp((w->z1() - ground_wall_z1)/window_vspacing));

			for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
				set_cube_zvals(trim, z, z+trim_height); // starts at floor height
				objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, 0, flags, 1.0, SHAPE_CUBE, trim_color); // floor trim
				if (!has_ceil_trim) continue;
				trim.z2() = z + floor_to_ceil_height; // ceil height
				trim.z1() = trim.z2() - trim_height;

				for (unsigned dir = 0; dir < 2; ++dir) { // for each side of wall
					cube_t ceil_trim(trim);
					ceil_trim.d[dim][!dir] = w->d[dim][dir];
					objs.emplace_back(ceil_trim, TYPE_WALL_TRIM, 0, dim, dir, flags, 1.0, SHAPE_ANGLED, trim_color); // ceiling trim
				}
			} // for f
		} // for w
	} // for d
	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) { // add trim for exterior walls
		if (is_basement(i)) continue; // skip basement walls because they're bare concrete
		bool const is_sec_bldg(i == get_real_parts_end());
		unsigned const num_floors(is_sec_bldg ? 1 : calc_num_floors(*i, window_vspacing, floor_thickness));

		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				cube_t trim(*i);
				trim.d[dim][!dir]  = i->d[dim][dir] + (dir ? -1.0 : 1.0)*trim_thickness;
				trim.d[dim][ dir] += (dir ? -1.0 : 1.0)*ext_wall_toler; // slight bias away from the exterior wall
				unsigned const ext_flags(flags | (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));
				float z(i->z1() + fc_thick);

				for (unsigned f = 0; f < num_floors; ++f, z += window_vspacing) {
					set_cube_zvals(trim, z, z+trim_height); // starts at floor height
					trim_cubes.clear();
					trim_cubes.push_back(trim); // start with entire length
					float const room_height(is_sec_bldg ? (i->dz() - floor_thickness) : floor_to_ceil_height);
					float const ceil_trim_z2(z + room_height), ceil_trim_z1(ceil_trim_z2 - trim_height); // ceil height

					for (auto j = parts.begin(); j != get_real_parts_end(); ++j) { // clip against other parts
						if (j == i) continue; // skip self
						cube_t clip_cube(*j);
						clip_cube.expand_in_dim(dim, expand_val); // expand to clip trim on the other side of the split wall
						subtract_cube_from_cubes(clip_cube, trim_cubes); // subtract this part from current trim cubes by clipping in XY
					}
					if (has_ceil_trim && is_house) { // houses have shorter doors and ceiling trim extends above the door, so draw full range
						for (auto c = trim_cubes.begin(); c != trim_cubes.end(); ++c) {
							cube_t trim(*c); // copy so that we can modify it
							set_cube_zvals(trim, ceil_trim_z1, ceil_trim_z2);
							objs.emplace_back(trim, TYPE_WALL_TRIM, 0, dim, !dir, flags, 1.0, SHAPE_ANGLED, trim_color); // ceiling trim
						}
					}
					if (f == 0) { // first floor, cut out areas for exterior doors
						for (auto d = doors.begin(); d != doors.end(); ++d) {
							cube_t door(d->get_bcube());
							bool const door_dim(door.dy() < door.dx());
							if (door_dim != bool(dim)) continue;
							door.expand_in_dim(door_dim, (expand_val + wall_thickness)); // expand to nonzero area; use a larger expand to account for distance door is offset away from ext wall
							subtract_cube_from_cubes(door, trim_cubes); // subtract this door from current trim cubes by clipping in XY
						}
					}
					for (auto c = trim_cubes.begin(); c != trim_cubes.end(); ++c) {
						objs.emplace_back(*c, TYPE_WALL_TRIM, 0, dim, 0, ext_flags, 1.0, SHAPE_CUBE, trim_color); // floor trim
						if (!has_ceil_trim || is_house) continue;
						set_cube_zvals(*c, ceil_trim_z1, ceil_trim_z2); // okay to edit in-place here
						objs.emplace_back(*c, TYPE_WALL_TRIM, 0, dim, !dir, flags, 1.0, SHAPE_ANGLED, trim_color); // ceiling trim
					}
				} // for f
			} // for dir
		} // for dim
	} // for i
	// add window trim
	if (is_rotated())   return; // not yet working for rotated buildings
	if (!has_windows()) return; // no windows
	float const border_mult(0.94); // account for the frame part of the window texture, which is included in the interior cutout of the window
	float const window_h_border(border_mult*get_window_h_border()), window_v_border(border_mult*get_window_v_border()); // (0, 1) range
	// Note: depth must be small to avoid object intersections; this applies to the windowsill as well
	float const window_trim_width(0.75*wall_thickness), window_trim_depth(1.0*trim_thickness), windowsill_depth(1.0*trim_thickness);
	float const window_offset(0.01*window_vspacing); // must match building_draw_t::add_section()
	static vect_vnctcc_t wall_quad_verts;
	wall_quad_verts.clear();
	get_all_drawn_window_verts_as_quads(wall_quad_verts);
	assert((wall_quad_verts.size() & 3) == 0); // must be a multiple of 4

	for (unsigned i = 0; i < wall_quad_verts.size(); i += 4) { // iterate over each quad
		auto const &v0(wall_quad_verts[i]);
		cube_t c(v0.v);
		float tx1(v0.t[0]), tx2(tx1), tz1(v0.t[1]), tz2(tz1); // tex coord ranges (xy, z); should generally be whole integers

		for (unsigned j = 1; j < 4; ++j) {
			auto const &vj(wall_quad_verts[i + j]);
			c.union_with_pt(vj.v);
			min_eq(tx1, vj.t[0]);
			max_eq(tx2, vj.t[0]);
			min_eq(tz1, vj.t[1]);
			max_eq(tz2, vj.t[1]);
		}
		if (tx1 == tx2 || tz1 == tz2) continue; // wall is too small to contain a window
		assert(tx2 - tx1 < 1000.0f && tz2 - tz1 < 1000.0f); // sanity check - less than 1000 windows in each dim
		assert(c.dz() > 0.0);
		bool const dim(c.dy() < c.dx()), dir(v0.get_norm()[dim] > 0.0);
		assert(c.get_sz_dim(dim) == 0.0); // must be zero size in one dim (X or Y oriented); could also use the vertex normal
		float const d_tx_inv(1.0f/(tx2 - tx1)), d_tz_inv(1.0f/(tz2 - tz1));
		float const window_width(c.get_sz_dim(!dim)*d_tx_inv), window_height(c.dz()*d_tz_inv); // window_height should be equal to window_vspacing
		float const border_xy(window_width*window_h_border), border_z(window_height*window_v_border), dscale(dir ? -1.0 : 1.0);
		cube_t window(c); // copy dim <dim>
		window.translate_dim(dim, dscale*window_offset);
		window.d[dim][!dir] += dscale*window_trim_depth; // add thickness on interior of building
		window.d[dim][ dir] += dscale*ext_wall_toler; // slight bias away from the exterior wall
		unsigned ext_flags(flags | (dir ? RO_FLAG_ADJ_HI : RO_FLAG_ADJ_LO));

		for (float z = tz1; z < tz2; z += 1.0) { // each floor
			float const bot_edge(c.z1() + (z - tz1)*window_height);
			set_cube_zvals(window, bot_edge+border_z, bot_edge+window_height-border_z);

			for (float xy = tx1; xy < tx2; xy += 1.0) { // windows along each wall
				float const low_edge(c.d[!dim][0] + (xy - tx1)*window_width);
				window.d[!dim][0] = low_edge + border_xy;
				window.d[!dim][1] = low_edge + window_width - border_xy;
				float const window_ar(window.get_sz_dim(!dim)/window.dz());
				float const side_trim_width(window_trim_width*((window_ar > 1.5) ? (window_ar - 0.5) : 1.0)); // widen for very wide windows to cover any holes at stretched edges
				cube_t top(window), bot(window), side(window);
				top.z1()  = window.z2();
				top.z2() += window_trim_width;
				bot.z2()  = window.z1();
				bot.z1() -= window_trim_width;
				bot.d[dim][!dir] += dscale*(windowsill_depth - window_trim_depth); // shift out further for windowsill
				top.expand_in_dim(!dim, side_trim_width);
				bot.expand_in_dim(!dim, side_trim_width);
				objs.emplace_back(top, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL, trim_color);
				objs.emplace_back(bot, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL, trim_color);

				for (unsigned s = 0; s < 2; ++s) { // left/right sides
					side.d[!dim][ s] = window.d[!dim][s] - (s ? -1.0 : 1.0)*side_trim_width;
					side.d[!dim][!s] = window.d[!dim][s];
					objs.emplace_back(side, TYPE_WALL_TRIM, 0, dim, dir, ext_flags, 1.0, SHAPE_TALL, trim_color);
				}
				add_window_coverings(window, dim, dir);
			} // for xy
		} // for z
	} // for i
}

void building_t::add_window_coverings(cube_t const &window, bool dim, bool dir) {
	// add blinds to some windows based on the containing room type for this floor
	bool is_split(0);
	int const room_id(get_room_id_for_window(window, dim, dir, is_split));
	if (room_id < 0) return; // room not found - should this be an error?
	if (is_split)    return; // window split across multiple rooms - how do we handle this? for now skip it
	room_t const &room(get_room(room_id));
	unsigned const floor(room.get_floor_containing_zval(window.zc(), get_window_vspace()));
	room_type const rtype(room.get_room_type(floor));

	switch (rtype) {
	case RTYPE_BED:  add_window_blinds  (window, dim, dir, room_id, floor); break; // bedroom
	case RTYPE_BATH: add_bathroom_window(window, dim, dir, room_id, floor); break; // bathroom
	} // end switch
}

void building_t::add_window_blinds(cube_t const &window, bool dim, bool dir, unsigned room_id, unsigned floor) {
	float const floor_spacing(get_window_vspace()), wall_thickness(get_wall_thickness()), extend(0.9*wall_thickness); // extend is 15% larger than window trim width
	vect_room_object_t &objs(interior->room_geom->objs);
	bool vertical((mat_ix + interior->rooms.size() + parts.size()) & 1); // something that's per-building
	
	if (vertical) { // check for horizontal wall clearance
		room_t const &room(get_room(room_id));
		if ((window.d[!dim][0] - 2.0*wall_thickness) < room.d[!dim][0] || (window.d[!dim][1] + 2.0*wall_thickness) > room.d[!dim][1]) {vertical = 0;} // not enough space for vertical blinds
	}
	rand_gen_t rgen;
	rgen.set_state((123*room_id + 211*interior->rooms.size()), (777*floor + 1));
	// open_amt is a mix of 50% room-based and 50% window-based to get somewhat consistent levels per room
	bool const full_open(rgen.rand_float() < 0.75);
	float const open_amt(0.9*(full_open ? 1.0 : 0.5*(rgen.rand_float() + fract(1123.7*objs.size())))); // 0-90% 25% the time, 90% for the rest
	cube_t c(window);
	c.d[dim][ dir] += (dir ? -1.0 : 1.0)*0.01*wall_thickness; // slight gap for wall trim
	c.d[dim][!dir] += (dir ? -1.0 : 1.0)*0.15*wall_thickness*(vertical ? 0.05 : (open_amt + 0.025)); // vertical blinds have no furniture clearance and can't bunch up
	c.expand_in_dim(2, extend); // extend in Z to cover window trim

	if (vertical) {
		c.expand_in_dim(!dim, 1.5*wall_thickness); // larger expand value (beyond the wall trim)
		float const center(c.get_center_dim(!dim)), half_width(0.5*c.get_sz_dim(!dim));
		float const shift_val(1.44*max(0.0f, (open_amt - 0.4f))*half_width); // more likely to be fully closed

		for (unsigned d = 0; d < 2; ++d) { // left, right
			cube_t c2(c);
			c2.d[!dim][!d] = center - (d ? -1.0 : 1.0)*shift_val;
			objs.emplace_back(c2, TYPE_BLINDS, room_id, dim, dir, (RO_FLAG_NOCOLL | (d ? RO_FLAG_ADJ_LO : RO_FLAG_ADJ_HI)), 1.0); // always fully lit
		}
	}
	else {
		c.expand_in_dim(!dim, extend); // expand width to cover trim +15% WT
		c.z2() += extend + 0.05*floor_spacing; // expand height to allow space for it to bunch up at the top
		c.z1() += open_amt*window.dz(); // raise amount is random per-room
		objs.emplace_back(c, TYPE_BLINDS, room_id, dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), 1.0); // always fully lit
	}
}

void building_t::add_bathroom_window(cube_t const &window, bool dim, bool dir, unsigned room_id, unsigned floor) { // frosted window blocks
	room_t const &room(get_room(room_id));
	if (count_ext_walls_for_room(room, window.z1()) != 1) return; // it looks odd to have window block walls at the corner of a building, so only enable this for single exterior walls
	vect_room_object_t &objs(interior->room_geom->objs);
	cube_t c(window);
	c.translate_dim(dim, (dir ? 1.0 : -1.0)*get_trim_thickness());
	unsigned const flags(RO_FLAG_NOCOLL | (room.is_lit_on_floor(floor) ? RO_FLAG_LIT : 0));
	objs.emplace_back(c, TYPE_WINDOW, room_id, dim, dir, flags, 1.0, SHAPE_CUBE, WHITE); // always lit
}

int building_t::get_room_id_for_window(cube_t const &window, bool dim, bool dir, bool &is_split) const {
	assert(interior);
	float const wall_thickness(get_wall_thickness());
	point const center(window.get_cube_center());

	for (auto r = interior->rooms.begin(); r != interior->rooms.end(); ++r) {
		if (center[!dim] < r->d[!dim][0] || center[!dim] >= r->d[!dim][1]) continue; // test center point for windows that straddle two rooms
		if (center.z < r->z1() || center.z > r->z2()) continue;
		if (fabs(center[dim] - r->d[dim][dir]) > wall_thickness) continue; // wrong wall
		is_split = (window.d[!dim][0] < r->d[!dim][0] || window.d[!dim][1] > r->d[!dim][1]);
		return (r - interior->rooms.begin()); // found
	} // for r
	return -1; // not found
}

void add_elevator_button(point const &pos, float button_radius, bool dim, bool dir, unsigned elevator_id, unsigned floor_id, bool inside, vect_room_object_t &objs) {
	cube_t c; c.set_from_point(pos);
	c.expand_in_dim(!dim, button_radius);
	c.expand_in_dim(2, button_radius); // Z
	c.d[dim][dir] += (dir ? 1.0 : -1.0)*0.25*button_radius;
	objs.emplace_back(c, TYPE_BUTTON, elevator_id, dim, dir, (RO_FLAG_NOCOLL | (inside ? RO_FLAG_IN_ELEV : 0)), 1.0, SHAPE_CYLIN, colorRGBA(1.0, 0.9, 0.5)); // room_id=elevator_id
	objs.back().obj_id = floor_id; // encode floor index as obj_id
}
void add_floor_number(unsigned floor_ix, unsigned floor_offset, bool has_parking_garage, ostringstream &oss) { // Note: floor_ix=1 is ground floor
	oss.str("");
	int const adj_floor_ix(int(floor_ix) - int(floor_offset));
	if (adj_floor_ix <= 0) {oss << (has_parking_garage ? "P" : "B") << (1 - adj_floor_ix);} // basement floors
	else {oss << adj_floor_ix;} // above ground floors
}
void set_floor_text_for_sign(room_object_t &sign, unsigned floor_ix, unsigned floor_offset, bool has_parking_garage, ostringstream &oss) { // Note: floor_ix=1 is ground floor
	add_floor_number(floor_ix, floor_offset, has_parking_garage, oss);
	sign.obj_id = register_sign_text(oss.str());
	int const adj_floor_ix(int(floor_ix) - int(floor_offset));
	float width_adj(0.0);
	if      (adj_floor_ix <= 0 ) {width_adj =  ((adj_floor_ix <  0) ? 0.1 : 0.0);} // basement floors; widen if lower than B1
	else if (adj_floor_ix <  10) {width_adj = -((adj_floor_ix == 1) ? 0.2 : 0.1);} // 1-10: make narrow
	else if (adj_floor_ix >= 20) {width_adj =  0.1;} // 20+: widen
	if (width_adj != 0.0) {sign.expand_in_dim(!sign.dim, width_adj*sign.get_sz_dim(!sign.dim));}
}
unsigned building_t::calc_floor_offset(float zval) const { // for basements
	return ((zval < ground_floor_z1) ? round_fp((ground_floor_z1 - zval)/get_window_vspace()) : 0);
}

void building_t::add_stairs_and_elevators(rand_gen_t &rgen) {

	float const window_vspacing(get_window_vspace()), floor_thickness(get_floor_thickness()), half_thick(0.5*floor_thickness);
	float const wall_thickness(get_wall_thickness()), elevator_car_z1_add(0.05*floor_thickness), fc_thick_scale(get_elevator_fc_thick_scale());
	vect_room_object_t &objs(interior->room_geom->objs);
	ostringstream oss; // reused across elevators/floors

	// add floor signs for U-shaped stairs
	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator || i->for_ramp || i->shape != SHAPE_U) continue; // not U-shaped stairs
		// stacked conn stairs start at floor 0 but are really the top floor of the part below; i->floor is not a global index and can't be used
		unsigned const floor_offset(calc_floor_offset(bcube.z1())); // use building z1 - should return number of underground levels
		unsigned const real_floor(round_fp((i->z1() - bcube.z1())/get_window_vspace()));
		point center;
		center[ i->dim] = i->d[i->dim][!i->dir]; // front of stairs
		center[!i->dim] = i->get_center_dim(!i->dim);
		if (has_parking_garage && i->z1() < ground_floor_z1) {center[!i->dim] += 0.25*i->get_sz_dim(!i->dim);} // shift to the side for parking garages to avoid center beams
		center.z = i->z1();
		cube_t sign(center);
		sign.d[i->dim][!i->dir] += (i->dir ? -1.0 : 1.0)*0.25*wall_thickness; // set sign thickness
		sign.expand_in_dim(!i->dim, 1.0*wall_thickness); // set sign width
		sign.z1() -= 2.5*wall_thickness; // set sign height
		objs.emplace_back(sign, TYPE_SIGN, 0, i->dim, !i->dir, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), 1.0, SHAPE_CUBE, DK_BLUE); // no room_id
		set_floor_text_for_sign(objs.back(), real_floor, floor_offset, has_parking_garage, oss);

		// if this is the top landing, we need to add a floor sign on the ceiling above it for the top floor
		if (i->is_at_top && !i->roof_access) {
			sign.translate_dim(2, window_vspacing); // move up one floor
			objs.emplace_back(sign, TYPE_SIGN, 0, i->dim, !i->dir, (RO_FLAG_NOCOLL | RO_FLAG_HANGING), 1.0, SHAPE_CUBE, DK_BLUE); // no room_id
			set_floor_text_for_sign(objs.back(), (real_floor + 1), floor_offset, has_parking_garage, oss);
		}
	} // for i

	// add elevator lights, and signs on each floor; must be done before setting buttons_start
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		// add light
		i->light_obj_id = objs.size();
		float const light_zval(max(ground_floor_z1, i->z1()) + elevator_car_z1_add + (1.0 - fc_thick_scale)*window_vspacing); // starts on the ground floor
		cube_t light(point(i->xc(), i->yc(), light_zval));
		light.z1() -= 0.02*window_vspacing;
		light.expand_by_xy(0.06*window_vspacing);
		objs.emplace_back(light, TYPE_LIGHT, i->room_id, i->dim, i->dir, (RO_FLAG_NOCOLL | RO_FLAG_IN_ELEV | RO_FLAG_LIT), 0.0, SHAPE_CYLIN, WHITE);
		objs.back().obj_id = uint16_t(i - interior->elevators.begin()); // encode elevator index as obj_id
		// add floor signs
		unsigned const num_floors(calc_num_floors(*i, window_vspacing, floor_thickness));
		unsigned const floor_offset(calc_floor_offset(i->z1()));
		float const ewidth(i->get_width());
		cube_t sign;
		sign.d[i->dim][0] = sign.d[i->dim][1] = i->d[i->dim][i->dir];
		sign.d[i->dim][i->dir] += (i->dir ? 1.0 : -1.0)*0.1*wall_thickness; // front of sign
		set_wall_width(sign, (i->d[!i->dim][1] - 0.1*ewidth), 0.04*ewidth, !i->dim); // to the high side, opposite the call button

		for (unsigned f = 0; f < num_floors; ++f) { // Note: floor number starts at 1 even if the elevator doesn't extend to the ground floor
			sign.z1() = i->z1()   + (f + 0.5)*window_vspacing;
			sign.z2() = sign.z1() + 0.1*ewidth;
			objs.emplace_back(sign, TYPE_SIGN, i->room_id, i->dim, i->dir, RO_FLAG_NOCOLL, 1.0, SHAPE_CUBE, DK_BLUE);
			set_floor_text_for_sign(objs.back(), f+1, floor_offset, has_parking_garage, oss);
		}
	} // for e
	interior->room_geom->buttons_start = objs.size();

	// add elevator buttons for each floor; must be done before setting stairs_start
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		float const button_radius(0.3*wall_thickness), ewidth(i->get_width());
		unsigned const num_floors(calc_num_floors(*i, window_vspacing, floor_thickness)), elevator_id(i - interior->elevators.begin());
		assert(num_floors > 1); // otherwise, why have an elevator?
		i->button_id_start = objs.size();

		// call buttons on each floor outside the elevator
		for (unsigned f = 0; f < num_floors; ++f) {
			point pos;
			pos[ i->dim] = i->d[i->dim][i->dir]; // front of the elevator
			pos[!i->dim] = i->d[!i->dim][0] + 0.1*ewidth; // to the low side
			pos.z = i->z1() + (f + 0.45)*window_vspacing;
			add_elevator_button(pos, button_radius, i->dim, i->dir, elevator_id, f, 0, objs);
		}
		// call buttons for each floor inside the elevator car; first find the panel location for the starting elevator car position
		cube_t elevator_car(*i);
		max_eq(elevator_car.z1(), ground_floor_z1); // always starts on the ground floor, not the bottom of the basement
		elevator_car.z1() += elevator_car_z1_add;
		elevator_car.z2()  = elevator_car.z1() + window_vspacing; // currently at the bottom floor
		cube_t const panel(get_elevator_car_panel(room_object_t(elevator_car, TYPE_ELEVATOR, elevator_id, i->dim, i->dir, 0), fc_thick_scale));
		float const dz(panel.dz()), button_spacing(dz/(num_floors + 1)); // add extra spacing on bottom and top of panel
		float const inner_button_radius(min(button_radius, min(0.35f*button_spacing, 0.25f*panel.get_sz_dim(!i->dim)))); // may need to be smaller
		point pos;
		pos[ i->dim] = panel.d[i->dim][!i->dir]; // front face of inside panel
		pos[!i->dim] = panel.get_center_dim(!i->dim) + 0.8*inner_button_radius; // a bit right of center to make room for floor number text
		
		for (unsigned f = 0; f < num_floors; ++f) {
			pos.z = panel.z1() + (f + 1)*button_spacing;
			add_elevator_button(pos, inner_button_radius, i->dim, !i->dir, elevator_id, f, 1, objs); // inside, pointing in opposite dir
		}
		i->button_id_end = objs.size();
	} // for e
	interior->room_geom->stairs_start  = objs.size();
	interior->room_geom->has_elevators = (!interior->elevators.empty());
	colorRGBA const railing_colors[3] = {GOLD, LT_GRAY, BLACK};
	colorRGBA const railing_color(railing_colors[rgen.rand()%3]); // set per-building

	for (auto i = interior->landings.begin(); i != interior->landings.end(); ++i) {
		if (i->for_elevator || i->for_ramp) continue; // for elevator or ramp, not stairs
		unsigned const num_stairs(i->get_num_stairs());
		float const stair_dz(window_vspacing/(num_stairs+1)), stair_height(stair_dz + floor_thickness);
		bool const dim(i->dim), dir(i->dir), has_side_walls(i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES || i->shape == SHAPE_U);
		bool const side(dir); // for U-shaped stairs; for now this needs to be consistent for the entire stairwell, can't use rgen.rand_bool()
		// Note: stairs always start at floor_thickness above the landing z1, ignoring landing z2/height
		float const tot_len(i->get_sz_dim(dim)), floor_z(i->z1() + floor_thickness - window_vspacing), step_len_pos(tot_len/num_stairs);
		float const wall_hw(min(max(0.15*step_len_pos, 0.15*stair_dz), 0.25*stair_dz));
		float step_len((dir ? 1.0 : -1.0)*step_len_pos), z(floor_z - floor_thickness), pos(i->d[dim][!dir]);
		cube_t stair(*i);

		if (i->shape != SHAPE_U) { // straight stairs
			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				stair.d[dim][!dir] = pos; stair.d[dim][dir] = pos + step_len;
				set_cube_zvals(stair, max(bcube.z1(), (z + 0.5f*half_thick)), z+stair_height); // don't go below the floor
				assert(stair.z1() < stair.z2());
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, dir); // Note: room_id=0, not tracked, unused
			}
		}
		else { // U-shaped stairs
			float const mid(i->get_center_dim(!dim));
			stair.d[!dim][side] = mid;
			step_len *= 2.0;

			for (unsigned n = 0; n < num_stairs; ++n, z += stair_dz, pos += step_len) {
				if (n == num_stairs/2) { // reverse direction and switch to other side
					step_len *= -1.0;
					stair.d[!dim][ side] = i->d[!dim][side];
					stair.d[!dim][!side] = mid;
				}
				assert(!(num_stairs & 1)); // require num_stairs to be an even number
				bool const is_rev(n >= num_stairs/2), stairs_dir(dir^is_rev);
				stair.d[dim][!stairs_dir] = pos; stair.d[dim][stairs_dir] = pos + step_len;
				set_cube_zvals(stair, max(floor_z, z), z+stair_height);
				objs.emplace_back(stair, TYPE_STAIR, 0, dim, side^is_rev); // Note: room_id=0, not tracked, unused
				objs.back().shape = SHAPE_STAIRS_U;
			} // for n
		}
		// add walls and railings
		bool const extend_walls_up(i->is_at_top && !i->roof_access); // space above is open, add a wall so that people can't fall down the stairs
		float const railing_z2(i->z2() + (i->roof_access ? 0.025*i->dz() : 0.0)); // capture z2 before we change it; move roof access railing up a bit to offset the shrink resize
		float const wall_bottom(floor_z - half_thick), railing_side_dz(0.5*stair_dz); // for U-shaped stairs
		cube_t wall(*i);
		if (extend_walls_up) {wall.z2() += window_vspacing - floor_thickness;}
		else {wall.z2() -= 0.5*floor_thickness;} // prevent z-fighting on top floor
		wall.z1() = max(bcube.z1()+half_thick, wall_bottom); // full height
		set_wall_width(wall, i->d[dim][dir], wall_hw, dim);

		if ((i->shape == SHAPE_WALLED && !(i->against_wall[0] || i->against_wall[1]) && (!i->stack_conn || !i->is_at_top)) || i->shape == SHAPE_U) {
			objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir); // add wall at back/end of stairs
		}
		else if ((i->shape == SHAPE_WALLED || i->shape == SHAPE_WALLED_SIDES) && extend_walls_up) { // add upper section only
			cube_t wall_upper(wall);
			set_wall_width(wall_upper, (i->d[dim][!dir] + (dir ? 1.0 : -1.0)*wall_hw), wall_hw, dim); // move to the other side
			wall_upper.z1() = railing_z2;

			for (unsigned d = 0; d < 2; ++d) {
				// if there's no wall, extend to cover the gap where the wall would be; slightly smaller to avoid z-fighting on exterior walls (happens to be the trim thickness);
				// we can't just skip this face for exterior walls because it may be visible through a window
				if (i->against_wall[d]) {wall_upper.d[!dim][d] += (d ? 1.0 : -1.0)*0.9*wall_thickness;}
			}
			objs.emplace_back(wall_upper, TYPE_STAIR_WALL, 0, dim, dir); // add wall at back/end of stairs
		}
		wall.d[dim][!dir] = i->d[dim][!dir];

		for (unsigned d = 0; d < 2; ++d) { // sides of stairs
			set_wall_width(wall, i->d[!dim][d], wall_hw, !dim);
			wall.expand_in_dim(dim, 0.01*wall_hw); // just enough to avoid z-fighting with stairs
			bool const add_wall(has_side_walls && !i->against_wall[d]); // don't add a wall if the stairs are already against a wall
			if (add_wall) {objs.emplace_back(wall, TYPE_STAIR_WALL, 0, dim, dir);} // add walls around stairs for this floor

			if (i->has_railing) { // add railings
				bool railing_dir(dir);
				cube_t railing(wall);
				unsigned flags(add_wall ? RO_FLAG_NOCOLL : 0);
				if (!has_side_walls) {flags |= RO_FLAG_OPEN;} // use this flag to indicate no walls, need balusters
				railing.z2() = railing_z2;

				if (add_wall || i->roof_access) {
					railing.translate_dim(!dim, (d ? -1.0 : 1.0)*2.0*wall_hw); // shift railing inside of walls
					railing.expand_in_dim( dim, -(i->roof_access ? 2.0 : 1.0)*wall_hw); // shrink slightly to avoid clipping through an end wall
				}
				if (i->shape == SHAPE_U) { // adjust railing height/angle to match stairs
					float const z_split(railing.zc());
					if (bool(d) == side) {railing.z1() = z_split + railing_side_dz; flags |= RO_FLAG_ADJ_HI; railing_dir ^= 1;}
					else                 {railing.z2() = z_split - railing_side_dz; flags |= RO_FLAG_ADJ_LO;}
				}
				objs.emplace_back(railing, TYPE_RAILING, 0, dim, railing_dir, flags, 1.0, SHAPE_CUBE, railing_color);
			}
		} // for d
		if (i->has_railing && i->shape == SHAPE_U) { // add a railing for the back wall of U-shaped stairs
			float const railing_zc(wall_bottom + 0.819*window_vspacing); // determined experimentally
			cube_t railing(*i);
			set_wall_width(railing, (i->d[dim][dir] + (dir ? -1.0 : 1.0)*2.0*wall_hw), wall_hw, dim);
			set_wall_width(railing, railing_zc, 1.4*railing_side_dz, 2); // set zvals
			objs.emplace_back(railing, TYPE_RAILING, 0, !dim, dir, (RO_FLAG_NOCOLL | RO_FLAG_ADJ_HI | RO_FLAG_ADJ_LO | RO_FLAG_ADJ_BOT), 1.0, SHAPE_CUBE, railing_color); // no ends
		}
		else if (i->has_railing && (i->stack_conn || (extend_walls_up && i->shape == SHAPE_STRAIGHT))) {
			// add railings around the top if: straight + top floor with no roof access, connector stairs, or basement stairs
			room_object_t railing(*i, TYPE_RAILING, 0, !dim, dir, (RO_FLAG_TOS | RO_FLAG_ADJ_BOT), 1.0, SHAPE_CUBE, railing_color); // flag to skip drawing ends
			railing.z1()  = railing.z2(); // starts at the floor
			railing.z2() += window_vspacing - floor_thickness;
			set_wall_width(railing, (i->d[dim][!dir] + (dir ? -1.0 : 1.0)*wall_hw), wall_hw, dim); // no overlap with stairs cutout

			for (unsigned d = 0; d < 2; ++d) {
				if (has_side_walls && !i->against_wall[d]) {railing.d[!dim][d] += (d ? -1.0 : 1.0)*2.0*wall_hw;} // shift railing inside of walls
			}
			// if the stairs extend to a wall, the back railing can be omitted; this is rare but happens in one house near the starting area
			cube_t railing_center(railing);
			railing_center.d[dim][dir] = railing_center.d[dim][!dir]; // shrink to zero area in this dim
			bool in_wall(0);

			for (auto const &w : interior->walls[dim]) {
				if (w.contains_cube(railing_center)) {in_wall = 1;}
			}
			if (!in_wall) {objs.emplace_back(railing);} // back railing
			railing.d[dim][dir] = i->d[dim][dir]; // extend to the front of the stairs
			railing.dim  ^= 1;
			railing.flags = RO_FLAG_TOS | RO_FLAG_ADJ_TOP; // flag so that no vertical pole is added

			for (unsigned d = 0; d < 2; ++d) { // sides of stairs
				railing.dir = bool(d);
				set_wall_width(railing, i->d[!dim][d], wall_hw, !dim);
				if (has_side_walls && !i->against_wall[d]) {railing.translate_dim(!dim, (d ? -1.0 : 1.0)*2.0*wall_hw);} // shift railing inside of walls
				objs.emplace_back(railing);
			}
		}
	} // for i (landings)
	for (auto i = interior->elevators.begin(); i != interior->elevators.end(); ++i) {
		unsigned const elevator_id(i - interior->elevators.begin()); // used for room_object_t::room_id
		cube_t elevator_car(*i);
		max_eq(elevator_car.z1(), ground_floor_z1); // always starts on the ground floor, not the bottom of the basement
		elevator_car.z1() += elevator_car_z1_add; // to prevent z-fighting when looking at the building from the bottom
		elevator_car.z2()  = elevator_car.z1() + window_vspacing; // one floor of height
		i->car_obj_id = objs.size();
		objs.emplace_back(elevator_car, TYPE_ELEVATOR, elevator_id, i->dim, i->dir, RO_FLAG_DYNAMIC);
		objs.back().drawer_flags = (uint16_t)calc_num_floors(*i, window_vspacing, floor_thickness); // store the number of floors in drawer_flags; used for drawing
		objs.back().item_flags   = (uint16_t)calc_floor_offset(i->z1()); // use correct starting floor index
	} // for i
}

int building_t::get_ext_door_dir(cube_t const &door_bcube, bool dim) const { // erturn value of 2 means 'not found'
	float const width(door_bcube.get_sz_dim(!dim));

	for (auto p = parts.begin(); p != get_real_parts_end(); ++p) { // find part containing this door so that we can get the correct dir
		if (is_basement(p)) continue; // skip the basement
		if (p->z1() != ground_floor_z1) continue; // not ground floor
		if (p->d[!dim][1] < door_bcube.d[!dim][1] || p->d[!dim][0] > door_bcube.d[!dim][0]) {continue;} // not contained in this dim
		if      (fabs(p->d[dim][0] - door_bcube.d[dim][0]) < 0.1*width) return 0;
		else if (fabs(p->d[dim][1] - door_bcube.d[dim][1]) < 0.1*width) return 1;
	} // for p
	cout << "Warning: Failed to find building exterior door: " << TXT(bcube.str()) << TXT(door_bcube.str()) << TXT(is_house) << endl; // debug printout
	//assert(0); // never gets here (too strong?)
	return 2; // not found
}

void building_t::add_sign_by_door(tquad_with_ix_t const &door, bool outside, std::string const &text, colorRGBA const &color, bool emissive) {
	cube_t const door_bcube(door.get_bcube());
	bool const dim(door_bcube.dy() < door_bcube.dx());
	int const dir_ret(get_ext_door_dir(door_bcube, dim));
	if (dir_ret > 1) return; // not found, skip sign
	bool dir(dir_ret != 0);
	float const width(door_bcube.get_sz_dim(!dim)), height(door_bcube.dz());
	cube_t c(door_bcube);

	if (outside) { // outside, place above the door
		c.z2() = door_bcube.z2() + 0.1*height;
	}
	else { // inside, place hanging near the top of the door
		c.z2() = door_bcube.z1() + get_window_vspace() - get_fc_thickness(); // right against the ceiling
	}
	c.z1() = c.z2() - 0.05*height;
	float const sign_width(0.8*text.size()*c.dz()), shrink(0.5f*(width - sign_width));
	c.expand_in_dim(!dim, -shrink);
	if (!outside) {dir ^= 1; c.translate_dim(dim, (dir ? 1.0 : -1.0)*0.1*height);} // move inside the building
	c.d[dim][dir] += (dir ? 1.0 : -1.0)*0.01*height;

	if (outside) {
		for (auto p2 = get_real_parts_end_inc_sec(); p2 != parts.end(); ++p2) {
			if (p2->intersects(c)) return; // sign intersects porch roof, skip this building
		}
	}
	unsigned flags(RO_FLAG_LIT | RO_FLAG_NOCOLL | (emissive ? RO_FLAG_EMISSIVE : 0) | (outside ? 0 : RO_FLAG_HANGING));
	vect_room_object_t &objs(interior->room_geom->objs);
	objs.emplace_back(c, TYPE_SIGN, 0, dim, dir, flags, 1.0, SHAPE_CUBE, color); // always lit; room_id is not valid
	objs.back().obj_id = register_sign_text(text);
}

void building_t::add_doorbell(tquad_with_ix_t const &door) {
	cube_t const door_bcube(door.get_bcube());
	bool const dim(door_bcube.dy() < door_bcube.dx());
	int const dir_ret(get_ext_door_dir(door_bcube, dim));
	if (dir_ret > 1) return; // not found, skip doorbell
	bool dir(dir_ret != 0);
	bool const side(dir ^ dim); // currently always to the right, which matches the door handle side
	float const door_width(door_bcube.get_sz_dim(!dim)), half_width(0.016*door_width), half_height(1.8*half_width);
	float const zval(door_bcube.z1() + 0.55*door_bcube.dz());
	float const pos(door_bcube.d[!dim][side] + (side ? 1.0 : -1.0)*5.0*half_width);
	cube_t c;
	c.d[dim][0  ]  = c.d[dim][1] = door_bcube.d[dim][dir] - 0.02*(dir ? 1.0 : -1.0)*get_window_vspace(); // slightly in front of exterior wall
	c.d[dim][dir] += (dir ? 1.0 : -1.0)*0.1*half_width;
	set_cube_zvals(c, (zval - half_height), (zval + half_height));
	set_wall_width(c, pos, half_width, !dim);
	interior->room_geom->objs.emplace_back(c, TYPE_BUTTON, 0, dim, dir, (RO_FLAG_LIT | RO_FLAG_NOCOLL), 1.0, SHAPE_CYLIN); // always lit; room_id is not valid
}

void building_t::add_exterior_door_signs(rand_gen_t &rgen) {
	if (is_house) { // maybe add welcome sign and add doorbell
		assert(!doors.empty());
		add_doorbell(doors.front());
		if (rgen.rand() & 3) return; // only 25% of houses have a welcome sign
		add_sign_by_door(doors.front(), 1, "Welcome", DK_BROWN, 0); // front door only, outside
	}
	else { // add exit signs
		if (pri_hall.is_all_zeros() && rgen.rand_bool()) return; // place exit signs on buildings with primary hallways and 50% of other buildings
		colorRGBA const exit_color(rgen.rand_bool() ? RED : GREEN);
		
		for (auto d = doors.begin(); d != doors.end(); ++d) {
			if (has_courtyard && (d+1) == doors.end()) break; // courtyard door is not an exit
			if (d->is_building_door()) {add_sign_by_door(*d, 0, "Exit", exit_color, 1);} // inside, emissive
		}
	}
}

room_t::room_t(cube_t const &c, unsigned p, unsigned nl, bool is_hallway_, bool is_office_, bool is_sec_bldg_) :
	cube_t(c), has_stairs(0), has_elevator(0), has_center_stairs(0), no_geom(is_hallway_), is_hallway(is_hallway_), is_office(is_office_), // no geom in hallways
	is_sec_bldg(is_sec_bldg_), interior(0), ext_sides(0), part_id(p), num_lights(nl), lit_by_floor(0), light_intensity(0.0)
{
	if      (is_sec_bldg) {assign_all_to(RTYPE_GARAGE);} // or RTYPE_SHED - will be set later
	else if (is_hallway)  {assign_all_to(RTYPE_HALL);}
	else if (is_office)   {assign_all_to(RTYPE_OFFICE);}
	else if (has_stairs)  {assign_all_to(RTYPE_STAIRS);} // not really correct since has_stairs is now a per-floor bit flag, but this will likely be overwritten later anyway
	else                  {assign_all_to(RTYPE_NOTSET);}
}
void room_t::assign_all_to(room_type rt) {
	for (unsigned n = 0; n < NUM_RTYPE_SLOTS; ++n) {rtype[n] = rt;}
}
void room_t::assign_to(room_type rt, unsigned floor) {
	min_eq(floor, NUM_RTYPE_SLOTS-1U); // room types are only tracked up to the 4th floor, and every floor above that has the same type as the 4th floor; good enough for houses at least
	if (rtype[floor] != RTYPE_BATH) {rtype[floor] = rt;} // assign unless already set to a bathroom, since we need that for has_bathroom()
}
bool room_t::has_bathroom() const {
	for (unsigned n = 0; n < NUM_RTYPE_SLOTS; ++n) {
		if (rtype[n] == RTYPE_BATH) return 1;
	}
	return 0;
}

