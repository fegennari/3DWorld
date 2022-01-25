// 3D World - Building Animals (rats, etc.)
// by Frank Gennari 1/16/22

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h" // for object_model_loader_t

float const RAT_FOV_DEG     = 60.0; // field of view in degrees
float const RAT_VIEW_FLOORS = 4.0; // view distance in floors
float const RAT_FOV_DP(cos(0.5*RAT_FOV_DEG*TO_RADIANS));

extern int animate2, camera_surf_collide, frame_counter;
extern float fticks;
extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;

float get_closest_building_sound(point const &at_pos, point &sound_pos, float floor_spacing);
sphere_t get_cur_frame_loudest_sound();


float rat_t::get_hwidth() const {
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

void building_t::update_animals(point const &camera_bs, unsigned building_ix) { // 0.01ms for 2 rats
	if (global_building_params.num_rats_max == 0 || !animate2) return;
	if (is_rotated() || !has_room_geom() || interior->rooms.empty()) return;
	vect_rat_t &rats(interior->room_geom->rats);
	if (rats.placed && rats.empty()) return; // no rats placed in this building
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) return; // no rat model

	if (!rats.placed) { // new building - place rats
		rand_gen_t rgen;
		rgen.set_state(building_ix+1, mat_ix+1); // unique per building
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
	// update rats
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame

	for (rat_t &rat : rats) {
		rgen.rand_mix(); // make sure it's different per rat
		update_rat(rat, camera_bs, rgen);
	}
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

bool can_hide_under(room_object_t const &c, float &zbot) {
	if (c.shape != SHAPE_CUBE ) return 0; // cubes only for now

	if (c.type == TYPE_CLOSET && c.is_open() && c.is_small_closet()) { // open small closet
		zbot = c.z2(); // closet ceiling
		return 1;
	}
	else if (c.type == TYPE_BED) {
		cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
		get_bed_cubes(c, cubes);
		zbot = cubes[0].z1(); // frame bottom
		return 1;
	}
	else if (c.type == TYPE_DESK || c.type == TYPE_DRESSER || c.type == TYPE_NIGHTSTAND || c.type == TYPE_TABLE) {
		cube_t cubes[5];
		get_table_cubes(c, cubes); // body and legs
		zbot = cubes[0].z1(); // body bottom
		return 1;
	}
	else if (c.type == TYPE_CHAIR) {
		cube_t cubes[3]; // seat, back, legs_bcube
		get_chair_cubes(c, cubes);
		zbot = cubes[0].z1(); // seat bottom
		return 1;
	}
	return 0;
}

void building_t::update_rat(rat_t &rat, point const &camera_bs, rand_gen_t &rgen) const {
	float const floor_spacing(get_window_vspace()), trim_thickness(get_trim_thickness()), view_dist(RAT_VIEW_FLOORS*floor_spacing);
	float const hlength(rat.get_hlength()), hwidth(rat.get_hwidth()), height(rat.get_height()), hheight(0.5*height);
	float const coll_radius(1.2f*hwidth); // slightly larger than half-width; maybe should use length so that the rat doesn't collide when turning?
	float const line_project_dist(max(1.1f*(hlength - coll_radius), 0.0f)); // extra space in front of the target destination
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	// set dist_thresh based on the distance we can move this frame; if set too low, we may spin in circles trying to turn to stop on the right spot
	float const move_dist(timestep*rat.speed), dist_thresh(3.0f*timestep*max(rat.speed, global_building_params.rat_speed));
	assert(hwidth <= hlength); // otherwise the model is probably in the wrong orientation
	bool collided(0);

	// move the rat
	if (rat.speed == 0.0) {rat.anim_time = 0.0;} // reset animation to rest pos
	else {
		rat.pos       += move_dist*rat.dir; // apply movement
		rat.anim_time += timestep*rat.speed;
		collided       = check_and_handle_dynamic_obj_coll(rat.pos, rat.radius, height, camera_bs);
	}
	vector3d const center_dz(0.0, 0.0, hheight);
	point const p1(rat.pos + center_dz);
	vector3d dir_to_fear;
	bool has_fear_dest(0);
	// apply scare logic
	bool const was_scared(rat.fear > 0.0);
	if (camera_surf_collide) {scare_rat(rat, camera_bs, 0.5, 1);} // the sight of the player walking in the building scares the rats
	sphere_t const cur_sound(get_cur_frame_loudest_sound());
	if (cur_sound.radius > 0.0) {scare_rat(rat, cur_sound.pos, 4.0*cur_sound.radius, 0);}
	// what about fear from sudden light changes? maybe it's enough that the light switch makes a click sound
	bool const is_scared(rat.fear > 0.0), newly_scared(is_scared && !was_scared);

	// determine destination
	if (is_scared) {
		// find hiding spot (pref in opposite direction from fear_pos);
		// we must check this each frame in case the player took or moved the object we were hiding under
		auto objs_end(interior->room_geom->get_std_objs_end()); // skip buttons/stairs/elevators
		float const rat_z1(rat.pos.z), rat_z2(rat.pos.z + height);
		point best_dest;
		float best_score(0.0), zbot(0.0);
		dir_to_fear   = (rat.fear_pos - rat.pos);
		dir_to_fear.z = 0.0; // XY plane only
		dir_to_fear.normalize();

		for (auto c = interior->room_geom->objs.begin(); c != objs_end; ++c) {
			if (c->z1() > rat_z2 || c->z2() < rat_z1) continue; // wrong floor
			if (!can_hide_under(*c, zbot)) continue;
			float const top_gap(zbot - (rat_z2 - 0.1*height)); // space between top of rat and bottom of object; allow the rat to squish a bit
			if (top_gap < 0.0) continue; // rat can't fit under this object
			point center(c->xc(), c->yc(), p1.z);
			float misalign(0.0);
			bool is_occupied(0);

			for (rat_t &other_rat : interior->room_geom->rats) {
				if (&other_rat == &rat) continue; // skip ourself
				float const move_dist(0.8f*(rat.radius + other_rat.radius) - p2p_dist_xy(center, other_rat.pos)); // smaller dist (head can overlap tail)
				
				if (move_dist > 0.0) { // another rat is in this spot
					center     += (p1 - center).get_norm()*move_dist; // move our target in front of this other rat
					misalign   += move_dist;
					is_occupied = 1;
				}
			} // for other_rat
			float const side_coverage(0.5f*min(c->dx(), c->dy()) - hlength - misalign); // amount of overhang of the object around the rat's extents
			float const dist(p2p_dist(p1, center)*(is_occupied ? 1.5 : 1.0)); // less desirable if occupied
			if (dist < dist_thresh) {has_fear_dest = 1; rat.speed = 0.0; break;} // already at this location, and it's valid, so stay there
			// Note: I tried using the dot product between this vector and dir_to_fear, but that causes instability when the rat is between two objects
			float const dist_to_fear(p2p_dist(rat.fear_pos, center));
			float const score((side_coverage - 0.5f*top_gap + 0.2f*dist_to_fear)/max(dist, dist_thresh)); // can be positive or negative
			if (best_score != 0.0 && score <= best_score) continue;
			if (check_line_coll_expand(p1, center, coll_radius, hheight)) continue;
			best_dest  = point(center.x, center.y, rat.pos.z); // keep zval on the floor
			best_score = score;
			if (center.x == rat.dest.x && center.y == rat.dest.y) break; // keep the same dest (optimization)
		} // for c
		if (!has_fear_dest && best_score != 0.0) { // found a valid hiding place; score can be positive or negative
			rat.dest = best_dest;
			if (dist_less_than(rat.pos, rat.dest, dist_thresh)) {rat.speed = 0.0;} // close enough - stop
			else {rat.speed = global_building_params.rat_speed*rgen.rand_uniform(1.0, 1.5);} // high speed if not at dest
			has_fear_dest = 1; // set to avoid triggering the code below if close to dest
			assert(rat.pos.z == rat.dest.z);
		}
		rat.fear = max(0.0f, (rat.fear - 0.2f*(timestep/TICKS_PER_SECOND))); // reduce fear over 5s
	}
	if (!has_fear_dest && (rat.speed == 0.0 || newly_scared || collided || dist_less_than(rat.pos, rat.dest, dist_thresh))) {
		// stopped, no dest, at dest, collided, or newly scared - choose a new dest
		cube_t valid_area(bcube);
		valid_area.expand_by_xy(-(hlength + trim_thickness));
		float target_fov_dp(RAT_FOV_DP), target_max_dist(view_dist); // start at nominal/max values
		float const dist_upper_bound((rat.fear > 0.8) ? 0.2 : 1.0); // shorten the distance if very scared, so that we can avoid more easily
		rat.speed = 0.0; // stop until we've found a valid destination

		for (unsigned n = 0; n < 200; ++n) { // make 100 tries
			if (n > 50) { // we've been at this for a while, maybe we need to relax our constraints, maybe follow the walls?
				target_fov_dp   -= 0.02; // allow for turns outside our field of view
				target_max_dist *= 0.96;  // decrease the max distance considered
			}
			vector3d const vdir(rgen.signed_rand_vector_xy().get_norm()); // random XY direction
			if (is_scared && n <= 100 && dot_product(dir_to_fear, vdir) > 0.0) continue; // don't move toward danger; may make the rat back into a corner
			if (dot_product(rat.dir, vdir) < target_fov_dp) continue; // not in field of view, use a new direction
			float const dist(rgen.rand_uniform(0.1, dist_upper_bound)*target_max_dist); // random distance out to max view dist
			if (dist <= dist_thresh) continue; // distance is too short
			point const cand(rat.pos + dist*vdir);
			if (!valid_area.contains_pt_xy(cand)) continue; // check for end point inside building bcube
			point const p2(cand + line_project_dist*vdir + center_dz); // extend in vdir so that the head doesn't collide
			point const p1_ext(p1 + coll_radius*vdir); // move the line slightly toward the dest to prevent collisions at the initial location
			if (check_line_coll_expand(p1_ext, p2, coll_radius, hheight)) continue;
			rat.dest  = cand;
			rat.speed = global_building_params.rat_speed*rgen.rand_uniform(0.5, 1.0)*(is_scared ? 1.5 : 1.0); // random speed
			break; // success
		} // for n
		assert(rat.pos.z == rat.dest.z);
	}
	// update direction
	bool const is_close(dist_less_than(rat.pos, rat.dest, dist_thresh));
	vector3d new_dir;
	if      (!is_close) {new_dir = (rat.dest - rat.pos).get_norm();} // point toward our destination
	else if (is_scared) {new_dir = dir_to_fear;} // point toward what we fear
	// else dir is unchanged

	if (new_dir != zero_vector) { // update dir if new_dir was set above
		float const delta_dir((is_scared ? 2.0 : 1.0)*min(1.0f, 1.5f*(1.0f - pow(0.7f, timestep)))); // higher turning rate when scared
		rat.dir = (delta_dir*new_dir + (1.0 - delta_dir)*rat.dir).get_norm();
	}
	if (rat.dir == zero_vector) {rat.dir = rgen.signed_rand_vector_xy().get_norm();} // dir must always be valid
	assert(rat.dir.z == 0.0); // must be in XY plane
}

void building_t::scare_rat(rat_t &rat, point const &scare_pos, float amount, bool by_sight) const {
	assert(amount > 0.0);
	if (rat.fear > 0.99 && dist_less_than(rat.fear_pos, scare_pos, rat.radius)) return; // already max fearful of this location (optimization)
	int const scare_room(get_room_containing_pt(scare_pos)), rat_room(get_room_containing_pt(rat.pos));
	assert(rat_room >= 0);
	if (rat_room != scare_room) {amount *= 0.67;} // less fearful if in a different room
	float const max_scare_dist(RAT_VIEW_FLOORS*get_window_vspace()), scare_dist(max_scare_dist*min(amount, 1.0f));
	float const fear((scare_dist - p2p_dist(rat.pos, scare_pos))/max_scare_dist);
	if (fear <= 0.0) return;
	if (by_sight && !check_line_of_sight_large_objs(rat.pos, scare_pos)) return; // check line of sight
	rat.fear     = min(1.0f, (rat.fear + fear));
	rat.fear_pos = scare_pos;
}

