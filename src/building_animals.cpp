// 3D World - Building Animals (rats, etc.)
// by Frank Gennari 1/16/22

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h" // for object_model_loader_t

float const RAT_FOV_DEG     = 60.0; // field of view in degrees
float const RAT_VIEW_FLOORS = 4.0; // view distance in floors
float const RAT_FOV_DP(cos(0.5*RAT_FOV_DEG*TO_RADIANS));

extern int animate2;
extern float fticks;
extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;
extern bldg_obj_type_t bldg_obj_types[];

bool maybe_inside_room_object(room_object_t const &obj, point const &pos, float radius);


float rat_t::get_width () const {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT)); // L, W, H
	return radius*sz.y/sz.x; // scale radius by ratio of width to length
}
float rat_t::get_height() const {
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT)); // L=3878, W=861, H=801
	return 2.0*radius*sz.z/max(sz.x, sz.y); // use max of x/y size; the x/y size represents the bcube across rotations
}
cube_t rat_t::get_bcube() const {
	cube_t bcube(pos, pos);
	bcube.expand_by_xy(radius);
	bcube.z2() += get_height();
	return bcube;
}
cube_t rat_t::get_bcube_with_dir() const {
	bool const pri_dim(fabs(dir.x) < fabs(dir.y));
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT));
	cube_t bcube(pos, pos);
	bcube.expand_in_dim( pri_dim, radius); // larger dim
	bcube.expand_in_dim(!pri_dim, radius*min(sz.x, sz.y)/max(sz.x, sz.y)); // smaller dim
	bcube.z2() += get_height();
	return bcube;
}

void building_t::update_animals(unsigned building_ix) {
	if (global_building_params.num_rats_max == 0 || !animate2) return;
	if (is_rotated() || !has_room_geom() || interior->rooms.empty()) return;
	vect_rat_t &rats(interior->room_geom->rats);
	if (rats.placed && rats.empty()) return; // no rats placed in this building
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) return; // no rat model
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, mat_ix+1); // unique per building

	if (!rats.placed) { // new building - place rats
		float const base_radius(0.1*get_window_vspace());
		unsigned const rmin(global_building_params.num_rats_min), rmax(global_building_params.num_rats_max);
		unsigned const num(rmin + ((rmin == rmax) ? 0 : (rgen.rand() % (rmax - rmin + 1))));
		rats.reserve(num);

		for (unsigned n = 0; n < num; ++n) {
			float const radius(base_radius*rgen.rand_uniform(0.8, 1.2));
			point const pos(gen_rat_pos(radius, rgen));
			if (pos == all_zeros) continue; // bad pos? skip this rat
			rats.emplace_back(pos, radius, rgen.signed_rand_vector_xy().get_norm());
		}
		rats.placed = 1; // even if there were no rats placed
	}
	for (rat_t &rat : rats) {update_rat(rat, rgen);}
}

point building_t::gen_rat_pos(float radius, rand_gen_t &rgen) const {
	for (unsigned n = 0; n < 100; ++n) { // make up to 100 tries
		unsigned const room_ix(rgen.rand() % interior->rooms.size());
		room_t const &room(interior->rooms[room_ix]);
		if (room.z1() > ground_floor_z1) continue; // not on the ground floor or basement
		cube_t place_area(room); // will represent the usable floor area
		place_area.expand_by_xy(-(radius + get_wall_thickness()));
		point pos(gen_xy_pos_in_area(place_area, radius, rgen));
		pos.z = place_area.z1() + get_fc_thickness(); // on top of the floor
		if (is_valid_ai_placement(pos, radius)) {return pos;} // check room objects; start in the open, not under something
	} // for n
	return all_zeros; // failed
}

bool line_int_cube_exp(point const &p1, point const &p2, cube_t const &cube, vector3d const &expand) {
	cube_t tc(cube);
	tc.expand_by(expand);
	return tc.line_intersects(p1, p2);
}
bool line_int_cubes_exp(point const &p1, point const &p2, cube_t const *const cubes, unsigned num_cubes, vector3d const &expand) {
	for (unsigned n = 0; n < num_cubes; ++n) {
		if (line_int_cube_exp(p1, p2, cubes[n], expand)) return 1;
	}
	return 0;
}
template<typename T> bool line_int_cubes_exp(point const &p1, point const &p2, vector<T> const &cubes, vector3d const &expand) {
	for (auto const &c : cubes) {
		if (line_int_cube_exp(p1, p2, c, expand)) return 1;
	}
	return 0;
}
bool can_hide_under(room_object_t const &c) {
	if (c.shape == SHAPE_CYLIN) return 0; // skip round tables for now
	if (c.shape != SHAPE_CUBE ) return 0; // cubes only for now
	return (c.type == TYPE_BED || c.type == TYPE_DESK || c.type == TYPE_DRESSER || c.type == TYPE_NIGHTSTAND || c.type == TYPE_TABLE || c.type == TYPE_CHAIR);
}

// p1 and p2 are line end points; radius applies in X and Y, hheight is half height and applies in +/- z
bool building_t::check_line_of_sight(point const &p1, point const &p2, float radius, float hheight) const {
	assert(p1.z == p2.z); // assumes same zval; can relax this later if needed
	assert(interior != nullptr);
	vector3d const expand(radius, radius, hheight);
	
	// check interior walls, doors, stairwells, and elevators
	for (unsigned d = 0; d < 2; ++d) {
		if (line_int_cubes_exp(p1, p2, interior->walls[d], expand)) return 0;
	}
	for (auto const &door : interior->doors) {
		if (door.open) {
			cube_t door_bounds(door);
			door_bounds.expand_by_xy(door.get_width());
			if (!line_int_cube_exp(p1, p2, door_bounds, expand)) continue; // optimization
			tquad_with_ix_t const door_tq(set_interior_door_from_cube(door));
			door_bounds = door_tq.get_bcube(); // somewhat more accurate
			door_bounds.expand_by_xy(door.get_thickness()); // conservative
			if (line_int_cube_exp(p1, p2, door_bounds, expand)) return 0;
		}
		else if (line_int_cube_exp(p1, p2, door, expand)) return 0;
	} // for door
	if (line_int_cubes_exp(p1, p2, interior->stairwells, expand)) return 0;
	if (line_int_cubes_exp(p1, p2, interior->elevators,  expand)) return 0;
	// check exterior walls
	vector3d cnorm; // unused
	float t(0.0); // unused
	if (ray_cast_exterior_walls(p1, p2, cnorm, t)) return 0;
	if (!has_room_geom()) return 1; // done (but really shouldn't get here)
	// check room objects; ignore expanded objects for now
	auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
	float const obj_z1(p1.z - hheight), obj_z2(p1.z + hheight);

	for (auto c = interior->room_geom->objs.begin(); c != objs_end; ++c) {
		if (c->no_coll() || !bldg_obj_types[c->type].ai_coll) continue; // skip non-colliding objects
		if (c->z1() > obj_z2 || c->z2() < obj_z1)   continue; // wrong floor
		if (!line_int_cube_exp(p1, p2, *c, expand)) continue;

		if (c->shape == SHAPE_CYLIN) { // vertical cylinder
			cylinder_3dw cylin(c->get_cylinder());
			cylin.p1.z -= hheight; cylin.p2.z += hheight; // extend top and bottom
			cylin.r1   += radius ; cylin.r2   += radius;
			//return !line_int_cylinder(p1, p2, cylin.p1, cylin.p2, cylin.r1, cylin.r2, 0, t); // check_ends=0
			return !line_intersect_cylinder_with_t(p1, p2, cylin, 0, t); // check_ends=0
		}
		if (c->type == TYPE_CLOSET) {
			cube_t cubes[5];
			get_closet_cubes(*c, cubes, 1); // get cubes for walls and door; for_collision=1
			// skip check of open doors for large closets since this case is more complex
			return !line_int_cubes_exp(p1, p2, cubes, ((c->is_open() && !c->is_small_closet()) ? 4U : 5U), expand);
		}
		else if (c->type == TYPE_BED) {
			cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
			get_bed_cubes(*c, cubes);
			if (line_int_cube_exp(p1, p2, cubes[0], expand)) return 0; // check bed frame (in case p1.z is high enough)
			get_tc_leg_cubes(cubes[5], 0.04, cubes); // head_width=0.04
			return !line_int_cubes_exp(p1, p2, cubes, 4, expand); // check legs
		}
		else if (c->type == TYPE_DESK || c->type == TYPE_DRESSER || c->type == TYPE_NIGHTSTAND || c->type == TYPE_TABLE) {
			cube_t cubes[5];
			get_table_cubes(*c, cubes); // body and legs
			return !line_int_cubes_exp(p1, p2, cubes, 5, expand);
		}
		else if (c->type == TYPE_CHAIR) {
			cube_t cubes[3], leg_cubes[4]; // seat, back, legs_bcube
			get_chair_cubes(*c, cubes);
			if (line_int_cube_exp(p1, p2, cubes[0], expand)) return 0; // check seat
			get_tc_leg_cubes(cubes[2], 0.15, leg_cubes); // width=0.15
			return !line_int_cubes_exp(p1, p2, leg_cubes, 4, expand); // check legs
		}
		else if (c->type == TYPE_STALL && maybe_inside_room_object(*c, p2, radius)) {} // TODO - inside test only applied to end point
		return 0; // intersection: no line of sight
	} // for c
	return 1;
}

void building_t::update_rat(rat_t &rat, rand_gen_t &rgen) const {
	float const floor_spacing(get_window_vspace()), view_dist(RAT_VIEW_FLOORS*floor_spacing);
	float const length(rat.get_length()), width(rat.get_width()), height(rat.get_height()), hheight(0.5*height);
	float const coll_radius(1.2*width); // maybe should use length so that the rat doesn't collide when turning?
	float const line_project_dist(0.5*1.1*(length - width)); // extra space in front of the target destination
	float const move_dist(fticks*rat.speed);

	if (rat.speed > 0.0) {
		// TODO: collision detection with dynamic objects: doors, balls, the player? people? other animals?
		rat.pos += move_dist*rat.dir; // apply movement
	}
	if (rat.fear > 0.0) { // find hiding spot (in opposite direction from fear_pos?)
		auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators

		for (auto c = interior->room_geom->objs.begin(); c != objs_end; ++c) {
			if (!can_hide_under(*c)) continue;
			if (c->z1() > (rat.pos.z + height) || c->z2() < rat.pos.z) continue; // wrong floor
			// TODO: consider this object as a hiding spot
		} // for c
		rat.fear = max(0.0f, (rat.fear - 0.1f*(fticks/TICKS_PER_SECOND))); // reduce fear over 10s
	}
	else if (dist_less_than(rat.pos, rat.dest, max(0.05f*rat.radius, 1.5f*move_dist))) { // no dest/at dest - choose a new dest
		vector3d const center_dz(0.0, 0.0, hheight);
		rat.speed = 0.0; // stop until we've found a valid destination
		float target_fov_dp(RAT_FOV_DP), target_max_dist(view_dist); // start at nominal/max values

		for (unsigned n = 0; n < 200; ++n) { // make 100 tries
			if (n > 50) { // we've been at this for a while, maybe we need to relax our constraints, maybe follow the walls?
				target_fov_dp   -= 0.02; // allow for turns outside our field of view
				target_max_dist *= 0.96;  // decrease the max distance considered
			}
			vector3d const vdir(rgen.signed_rand_vector_xy().get_norm()); // random XY direction
			if (dot_product(rat.dir, vdir) < target_fov_dp) continue; // not in field of view, use a new direction
			float const dist(rgen.rand_uniform(0.1, 1.0)*target_max_dist); // random distance out to max view dist
			point const cand(rat.pos + dist*vdir);
			if (!bcube.contains_pt_xy_exp(cand, -coll_radius)) continue; // check for end point inside building bcube
			point const p1(rat.pos + center_dz), p2(cand + line_project_dist*vdir + center_dz); // extend in vdir so that the head doesn't collide
			if (!check_line_of_sight(p1, p2, coll_radius, hheight)) continue;
			rat.dest  = cand;
			rat.speed = global_building_params.rat_speed*rgen.rand_uniform(0.5, 1.0); // random speed
			break; // success
		} // for n
	}
	if (rat.dest != rat.pos) {
		float const delta_dir(min(1.0f, 1.5f*(1.0f - pow(0.7f, fticks))));
		vector3d const new_dir((rat.dest - rat.pos).get_norm());
		rat.dir = (delta_dir*new_dir + (1.0 - delta_dir)*rat.dir).get_norm();
	}
	assert(rat.dir.z == 0.0); // must be in XY plane
}

void building_t::scare_animals(point const &scare_pos, float sight_amt, float sound_amt) {
	if (!animate2) return; // don't accumulate fear when AI is not updating
	vector<rat_t> &rats(interior->room_geom->rats);
	if (rats.empty()) return; // nothing to do
	float const amount(min((sight_amt + sound_amt), 1.0f)); // for now we don't differentiate between sight and sound for rats
	assert(amount > 0.0);
	float const max_scare_dist(3.0*get_window_vspace()), scare_dist(max_scare_dist*amount);
	int const scare_room(get_room_containing_pt(scare_pos));
	if (scare_room < 0) return; // error?

	for (rat_t &rat : rats) {
		if (rat.fear > 0.99) continue; // already max fearful (optimization)
		float const dist(p2p_dist(rat.pos, scare_pos));
		if (dist < scare_dist) continue; // optimization
		int const rat_room(get_room_containing_pt(rat.pos));
		assert(rat_room >= 0);
		if (rat_room != scare_room) continue; // only scared if in the same room
		float fear(amount);

		if (sight_amt > 0.0 && !check_line_of_sight(rat.pos, scare_pos, 0.0, 0.0)) { // check line of sight with zero radius/height line
			if (sound_amt == 0.0) continue; // not visible and no sound = no fear
			fear = min(sound_amt, 1.0f); // fear from sound only
		}
		fear = fear*max_scare_dist - dist;
		if (fear <= 0.0) continue;
		rat.fear     = min(1.0f, (rat.fear + fear));
		rat.fear_pos = scare_pos;
	} // for rat
}
