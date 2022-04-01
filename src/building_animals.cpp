// 3D World - Building Animals (rats, etc.)
// by Frank Gennari 1/16/22

#include "3DWorld.h"
#include "function_registry.h"
#include "buildings.h"
#include "city_model.h" // for object_model_loader_t
#include "openal_wrap.h"

float const RAT_FOV_DEG      = 60.0; // field of view in degrees
float const RAT_VIEW_FLOORS  = 4.0; // view distance in floors
float const RAT_FEAR_SPEED   = 1.3; // multiplier
float const RAT_ATTACK_SPEED = 1.2; // multiplier
float const RAT_FOV_DP(cos(0.5*RAT_FOV_DEG*TO_RADIANS));
float const SPIDER_VIEW_FLOORS = 4.0; // view distance in floors

extern int animate2, camera_surf_collide, frame_counter, display_mode;
extern float fticks;
extern double tfticks;
extern building_params_t global_building_params;
extern object_model_loader_t building_obj_model_loader;

float get_closest_building_sound(point const &at_pos, point &sound_pos, float floor_spacing);
sphere_t get_cur_frame_loudest_sound();
bool in_building_gameplay_mode();
bool player_take_damage(float damage_scale, bool poisoned=0, bool *has_key=nullptr);
cube_t get_true_room_obj_bcube(room_object_t const &c);


void building_animal_t::sleep_for(float time_secs_min, float time_secs_max) {
	wake_time = (float)tfticks + rand_uniform(time_secs_min, time_secs_max)*TICKS_PER_SECOND;
	dist_since_sleep = 0.0; // reset the counter
}
void building_animal_t::move(float timestep) {
	// update animation time using position change; note that we can't just do the update in the rat movement code below because pos may be reset in case of collision
	anim_time += p2p_dist(pos, last_pos)/radius; // scale with size so that small rat/spider legs move faster
	last_pos   = pos;
	if (anim_time > 100000.0) {anim_time = 0.0;} // reset animation after awhile to avoid precision problems; should this be done when the player isn't looking?

	if (is_sleeping()) {
		if ((float)tfticks > wake_time) {wake_time = speed = 0.0;} // time to wake up
	}
	else if (speed == 0.0) {
		anim_time = 0.0; // reset animation to rest pos
	}
	else { // apply movement and check for collisions with dynamic objects
		vector3d const dest_dir((dest - pos).get_norm());

		if (dot_product(dest_dir, dir) > 0.75) { // only move if we're facing our dest, to avoid walking through an object
			float const move_dist(timestep*speed);
			pos               = pos + move_dist*dir;
			dist_since_sleep += move_dist;
		}
	}
}

void building_t::update_animals(point const &camera_bs, unsigned building_ix, int ped_ix) {
	if (!animate2 || is_rotated() || !has_room_geom() || interior->rooms.empty()) return;
	update_rats   (camera_bs, building_ix, ped_ix);
	update_spiders(camera_bs, building_ix, ped_ix);
}

template<typename T> void building_t::add_animals_on_floor(T &animals, unsigned building_ix, unsigned num_min, unsigned num_max, float sz_min, float sz_max) const {
	if (animals.placed) return; // already placed
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, mat_ix+1); // unique per building
	float const floor_spacing(get_window_vspace());
	unsigned const num(num_min + ((num_min == num_max) ? 0 : (rgen.rand() % (num_max - num_min + 1))));
	animals.reserve(num);

	for (unsigned n = 0; n < num; ++n) {
		// there's no error check for min_sz <= max_sz, so just use min_sz in that case
		float const sz_scale((sz_min >= sz_max) ? sz_min : rgen.rand_uniform(sz_min, sz_max)), radius(0.5f*floor_spacing*sz_scale);
		point const pos(gen_animal_floor_pos(radius, rgen));
		if (pos == all_zeros) continue; // bad pos? skip this animal
		animals.add(typename T::value_type(pos, radius, rgen.signed_rand_vector_xy().get_norm()));
	}
	animals.placed = 1; // even if there were no animals placed
}

point building_t::gen_animal_floor_pos(float radius, rand_gen_t &rgen) const {
	for (unsigned n = 0; n < 100; ++n) { // make up to 100 tries
		unsigned const room_ix(rgen.rand() % interior->rooms.size());
		room_t const &room(interior->rooms[room_ix]);
		if (room.z1() > ground_floor_z1) continue; // not on the ground floor or basement
		cube_t place_area(room); // will represent the usable floor area
		place_area.expand_by_xy(-(radius + get_wall_thickness()));
		if (min(place_area.dx(), place_area.dy()) < 4.0*radius) continue; // room too small (can happen for has_complex_floorplan office buildings)
		point pos(gen_xy_pos_in_area(place_area, radius, rgen));
		pos.z = place_area.z1() + get_fc_thickness(); // on top of the floor
		if (is_valid_ai_placement(pos, radius)) {return pos;} // check room objects; start in the open, not under something
	} // for n
	return all_zeros; // failed
}

bool building_t::is_pos_inside_building(point const &pos, float xy_pad, float hheight) const {
	if (!bcube.contains_pt_xy_exp(pos, -xy_pad)) return 0; // check for end point inside building bcube
	cube_t req_area(pos, pos);
	req_area.expand_by_xy(xy_pad);
	req_area.z2() += hheight;
	return is_cube_contained_in_parts(req_area);
}

void update_dir_incremental(vector3d &dir, vector3d const &new_dir, float turn_rate, float timestep, rand_gen_t &rgen) {
	if (new_dir != zero_vector) { // update dir if new_dir was set above
		float const delta_dir(turn_rate*min(1.0f, 1.5f*(1.0f - pow(0.7f, timestep))));
		dir = (delta_dir*new_dir + (1.0 - delta_dir)*dir).get_norm();
	}
	if (dir == zero_vector) {dir = rgen.signed_rand_vector_xy().get_norm();} // dir must always be valid
}
void try_resolve_coll(vector3d const &dir, vector3d const &upv, vector3d const &coll_dir, vector3d &new_dir, bool try_tangent) {
	if (dot_product(coll_dir, new_dir) > 0.0) { // must move away from the collision direction
		if (try_tangent) { // try to preserve direction by moving in a tangent to the collider
			new_dir = cross_product(coll_dir, upv).get_norm(); // should be orthogonal, but we can be safe and normalize
			if (dot_product(new_dir, dir) < 0.0) {new_dir.negate();} // there are two solutions; choose the one closer to our current dir
		}
		else {new_dir.negate();} // otherwise reverse direction
	}
}

bool play_attack_sound(point const &pos, float gain, float pitch, rand_gen_t &rgen) {
	static double last_sound_time(0.0);

	if (tfticks - last_sound_time > 0.4*TICKS_PER_SECOND) {
		gen_sound_thread_safe(SOUND_SQUISH, pos, gain, pitch);
		last_sound_time = tfticks + 0.2f*TICKS_PER_SECOND*rgen.rand_float(); // add some randomness
		return 1;
	}
	return 0;
}


// *** Rats ***

rat_t::rat_t(point const &pos_, float radius_, vector3d const &dir_) : building_animal_t(pos_, radius_, dir_),
fear(0.0), is_hiding(0), near_player(0), attacking(0)
{
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT)); // L=3878, W=861, H=801
	hwidth = radius*sz.y/sz.x; // scale radius by ratio of width to length
	height = 2.0*radius*sz.z/max(sz.x, sz.y); // use max of x/y size; the x/y size represents the bcube across rotations
}
cube_t rat_t::get_bcube() const {
	cube_t bcube(pos, pos);
	bcube.expand_by_xy(radius);
	bcube.z2() += height;
	return bcube;
}
cube_t rat_t::get_bcube_with_dir() const {
	bool const pri_dim(fabs(dir.x) < fabs(dir.y));
	vector3d const sz(building_obj_model_loader.get_model_world_space_size(OBJ_MODEL_RAT));
	cube_t bcube(pos, pos);
	bcube.expand_in_dim( pri_dim, radius); // larger dim
	bcube.expand_in_dim(!pri_dim, radius*min(sz.x, sz.y)/max(sz.x, sz.y)); // smaller dim
	bcube.z2() += height;
	return bcube;
}

bool building_t::add_rat(point const &pos, float hlength, vector3d const &dir, point const &placed_from) {
	point rat_pos(pos);
	if (!get_zval_of_floor(pos, hlength, rat_pos.z)) return 0; // place on the floor, skip if there's no floor here
	rat_t rat(rat_pos, hlength, vector3d(dir.x, dir.y, 0.0).get_norm()); // dir in XY plane
	
	if (check_line_coll_expand(pos, rat_pos, hlength, rat.height)) { // something is in the way
		point const test_pos(rat_pos + vector3d(0.0, 0.0, rat.height));
		auto objs_end(interior->room_geom->get_placed_objs_end()); // skip buttons/stairs/elevators

		for (auto i = interior->room_geom->objs.begin(); i != objs_end; ++i) { // check for a toilet that we can drop the rat into
			if (i->type == TYPE_TOILET && i->contains_pt(test_pos) && i->room_id == get_room_containing_pt(placed_from)) {
				point const sound_origin(i->get_cube_center());
				gen_sound_thread_safe(SOUND_FLUSH, local_to_camera_space(sound_origin));
				register_building_sound(sound_origin, 0.5);
				register_achievement("Sleep with the Fishes");
				return 1;
			}
		}
		return 0; // can't place the rat here
	}
	rat.fear_pos = placed_from;
	rat.fear     = 1.0; // starts off with max fear
	interior->room_geom->rats.add(rat);
	interior->room_geom->modified_by_player = 1;
	return 1;
}

void building_t::update_rats(point const &camera_bs, unsigned building_ix, int ped_ix) {
	if (global_building_params.num_rats_max == 0) return;
	vect_rat_t &rats(interior->room_geom->rats);
	if (rats.placed && rats.empty()) return; // no rats placed in this building
	if (!building_obj_model_loader.is_model_valid(OBJ_MODEL_RAT)) return; // no rat model
	//timer_t timer("Update Rats"); // multi-part: 1.1ms, open door 1.7ms; office building 1.7ms
	add_animals_on_floor(rats, building_ix, global_building_params.num_rats_min, global_building_params.num_rats_max,
		global_building_params.rat_size_min, global_building_params.rat_size_max);
	// update rats
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	unsigned num_near_player(0);
	point rat_alert_pos;

	for (rat_t &rat : rats) { // must be done before sorting
		rat.move(timestep);
		num_near_player += rat.near_player;
		if (num_near_player == global_building_params.min_attack_rats) {rat_alert_pos = rat.pos;}
	}
	static bool prev_can_attack_player(0);
	bool const can_attack_player(num_near_player >= global_building_params.min_attack_rats);
	
	if (can_attack_player && !prev_can_attack_player) { // play sound on first attack
		gen_sound_thread_safe(SOUND_RAT_SQUEAK, local_to_camera_space(rat_alert_pos), 1.0, 1.2); // high pitch
	}
	prev_can_attack_player = can_attack_player;
	rats.do_sort();
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame

	if (frame_counter & 1) { // reverse iteration, to avoid directional bias
		for (auto r = rats.rbegin(); r != rats.rend(); ++r) {update_rat(*r, camera_bs, ped_ix, timestep, rats.max_xmove, can_attack_player, rgen);}
	}
	else { // forward iteration; ~0.004ms per rat
		for (auto r = rats. begin(); r != rats. end(); ++r) {update_rat(*r, camera_bs, ped_ix, timestep, rats.max_xmove, can_attack_player, rgen);}
	}
}

bool can_hide_under(room_object_t const &c, cube_t &hide_area) {
	cube_t dishwasher; // used below

	if (c.type == TYPE_CLOSET && c.is_open() && c.is_small_closet()) { // open small closet
		hide_area = c;
		hide_area.expand_by(-get_closet_wall_thickness(c)); // we want the inside of the closet, excluding the walls
		hide_area.z1() += 0.5*hide_area.dz(); // use the halfway point; somewhat arbitrary, but will affect the score
		return 1;
	}
	else if (c.type == TYPE_BED) {
		cube_t cubes[6]; // frame, head, foot, mattress, pillow, legs_bcube
		get_bed_cubes(c, cubes);
		hide_area = cubes[0]; // frame
		return 1;
	}
	else if (c.type == TYPE_DESK || c.type == TYPE_TABLE) {
		cube_t cubes[5];
		get_table_cubes(c, cubes); // body and legs
		hide_area = cubes[0]; // body
		return 1;
	}
	else if (c.type == TYPE_DRESSER || c.type == TYPE_NIGHTSTAND) {
		hide_area = get_dresser_middle(c);
		return 1;
	}
	else if (c.type == TYPE_CHAIR) {
		cube_t cubes[3]; // seat, back, legs_bcube
		get_chair_cubes(c, cubes);
		hide_area = cubes[0]; // seat
		return 1;
	}
	else if (c.type == TYPE_BCASE) {
		cube_t top, middle, back, lr[2];
		get_bookcase_cubes(c, top, middle, back, lr);
		hide_area = middle;
		return 1;
	}
	else if (c.type == TYPE_KSINK && get_dishwasher_for_ksink(c, dishwasher)) {
		hide_area = dishwasher;
		hide_area.d[c.dim][!c.dir] = c.d[c.dim][!c.dir]; // use the back of the cabinet, not the back of the dishwasher door
		return 1;
	}
	else if (c.type == TYPE_COUCH) {
		hide_area = c;
		hide_area.z1() += 0.06*c.dz(); // there's space under the couch
		return 1;
	}
	else if (c.type == TYPE_BRSINK) { // office building bathroom sink
		// not a very good hiding spot, but there aren't many in office buildings;
		// this is just a placeholder anyway, since sinks don't extend down to the floor and won't pass the rat zval test
		hide_area = c;
		return 1;
	}
	else if (c.type == TYPE_SHELVES) {
		if (c.get_num_shelves() < 4) return 0; // only 4 shelves provide a small enough space at the bottom to hide under
		hide_area = c;
		hide_area.z1() += 0.05*c.dz();
		return 1;
	}
	else if (c.type == TYPE_COLLIDER && (c.flags & RO_FLAG_FOR_CAR)) { // parked car
		hide_area = c;
		hide_area.z1() += 0.12*c.dz(); // there's space under the car
		return 1;
	}
	//else if (c.type == TYPE_STAIR) {} // what about hiding under the stairs? only for SHAPE_WALLED_SIDES?
	return 0;
}

class dir_gen_t {
	vector<vector3d> dirs;
	unsigned dir_ix;

	void gen_dirs() {
		rand_gen_t rgen;
		dirs.resize(1000);
		for (auto &dir : dirs) {dir = rgen.signed_rand_vector_spherical_xy().get_norm();}
	}
public:
	dir_gen_t() : dir_ix(0) {}

	vector3d const &gen_dir() {
		if (dirs.empty()) {gen_dirs();}
		vector3d const &dir(dirs[dir_ix++]);
		if (dir_ix == dirs.size()) {dir_ix = 0;} // wrap around
		return dir;
	}
};
dir_gen_t dir_gen;

void building_t::update_rat(rat_t &rat, point const &camera_bs, int ped_ix, float timestep, float &max_xmove, bool can_attack_player, rand_gen_t &rgen) const {

	float const floor_spacing(get_window_vspace()), trim_thickness(get_trim_thickness()), view_dist(RAT_VIEW_FLOORS*floor_spacing);
	float const hlength(rat.get_hlength()), hwidth(rat.hwidth), height(rat.height), hheight(0.5*height);
	float const squish_hheight(0.75*hheight); // rats can squish to get under low objects and walk onto small steps
	float const coll_radius(1.2f*hwidth); // slightly larger than half-width; maybe should use length so that the rat doesn't collide when turning?
	float const line_project_dist(max(1.1f*(hlength - coll_radius), 0.0f)); // extra space in front of the target destination
	// set dist_thresh based on the distance we can move this frame; if set too low, we may spin in circles trying to turn to stop on the right spot
	float const dist_thresh(2.0f*timestep*max(rat.speed, global_building_params.rat_speed));
	float const xy_pad(hlength + trim_thickness);
	vector3d const center_dz(0.0, 0.0, hheight); // or squish_hheight?
	assert(hwidth <= hlength); // otherwise the model is probably in the wrong orientation
	bool update_path(0);
	vector3d coll_dir;
	point const prev_pos(rat.pos); // capture the pre-collision point
	rgen.rand_mix(); // make sure it's different per rat

	if (rat.is_sleeping() && rat.fear == 0.0) {} // peacefully sleeping, no collision needed
	else if (check_and_handle_dynamic_obj_coll(rat.pos, rat.radius, rat.pos.z, (rat.pos.z + height), camera_bs, 0)) { // check for collisions; for_spider=0
		coll_dir = (prev_pos - rat.pos).get_norm(); // points toward the collider in the XY plane

		// check if new pos is valid, and has a path to dest
		if (!is_pos_inside_building(rat.pos, xy_pad, hheight)) {
			rat.pos = prev_pos; // restore previous pos before collision
			rat.sleep_for(0.1, 0.2); // wait 0.1-0.2s so that we don't immediately collide and get pushed out again
		}
		else if (check_line_coll_expand((rat.pos + center_dz), (rat.dest + center_dz), coll_radius, squish_hheight)) {
			rat.pos = prev_pos; // restore previous pos before collision
		}
		else {
			max_eq(max_xmove, fabs(rat.pos.x - prev_pos.x));
			// update the path about every 30 frames of colliding; this prevents the rat from being stuck while also avoiding jittering due to frequent dest updates;
			// using randomness rather than updating every 30 frames makes this process less regular and mechanical
			update_path = ((rgen.rand()%30) == 0);
		}
	}
	point const p1(rat.pos + center_dz);
	vector3d dir_to_fear;
	bool has_fear_dest(0);
	// apply scare logic
	bool const was_scared(rat.fear > 0.0);
	scare_rat(rat, camera_bs, ped_ix);
	rat.attacking = (rat.near_player && can_attack_player);
	if (rat.attacking) {rat.fear = 0.0;} // no fear when attacking
	bool const is_scared(rat.fear > 0.0), newly_scared(is_scared && !was_scared);

	// determine destination
	if (rat.attacking) {
		float const player_radius(CAMERA_RADIUS*global_building_params.player_coll_radius_scale), min_dist(player_radius + hlength);
		point target(camera_bs.x, camera_bs.y, rat.pos.z);
		vector3d const vdir((target - rat.pos).get_norm());
		target -= vdir*(1.01*min_dist); // get within attacking range, but not at the center of the player
		
		if (is_pos_inside_building(target, xy_pad, hheight)) { // check if outside the valid area
			point const p1_ext(p1 + coll_radius*vdir); // move the line slightly toward the dest to prevent collisions at the initial location
			
			if (!check_line_coll_expand(p1_ext, target, coll_radius, squish_hheight)) {
				rat.dest      = target;
				rat.speed     = RAT_ATTACK_SPEED*global_building_params.rat_speed;
				rat.wake_time = 0.0; // wake up
				update_path   = 0;

				if (dist_xy_less_than(rat.pos, target, 0.05*min_dist)) { // do damage when nearly colliding with the player
					play_attack_sound(local_to_camera_space(rat.pos), 1.0, 1.0, rgen);
					if (player_take_damage(0.004, 0)) {register_achievement("Rat Food");} // damage over time; achievement if the player dies
				}
			}
		}
	}
	if (is_scared) {
		// find hiding spot (pref in opposite direction from fear_pos);
		// we must check this each frame in case the player took or moved the object we were hiding under
		float const rat_z1(rat.pos.z), rat_z2(rat.pos.z + height), rat_squish_z2(p1.z + squish_hheight);
		point best_dest;
		float best_score(0.0);
		dir_to_fear   = (rat.fear_pos - rat.pos);
		dir_to_fear.z = 0.0; // XY plane only
		dir_to_fear.normalize();
		rat.wake_time = 0.0; // wake up
		vect_room_object_t::const_iterator b, e;
		get_begin_end_room_objs_on_ground_floor(rat_z2, 0, b, e); // for_spider=0

		for (auto c = b; c != e; ++c) {
			if (c->z1() > rat_z2 || c->z2() < rat_z1) continue; // wrong floor, or object not on the floor
			cube_t hide_area; // will be a subset of the object
			if (c->shape != SHAPE_CUBE || !can_hide_under(*c, hide_area)) continue; // only cubes for now
			float const top_gap(hide_area.z1() - rat_squish_z2); // space between top of rat and bottom of object
			if (top_gap < 0.0) continue; // rat can't fit under this object
			if (!dist_xy_less_than(hide_area.get_cube_center(), p1, view_dist)) continue; // too far away to see
			// select our destination under this hiding spot;
			// this must be unique per rat so that rats don't compete for the exact same spot, and must be the same across calls for stability;
			// also use the obj_id to mix things up between objects, and throw in the type in case obj_id is left at 0;
			// we can't use the obj vector position because it may change if the player takes or drops objects
			rand_gen_t my_rgen;
			my_rgen.set_state((rat.id + 1), (c->obj_id + (c->type << 16) + 1));
			cube_t safe_area(hide_area);
			point cand_dest(0.0, 0.0, p1.z); // x/y will be set below

			for (unsigned d = 0; d < 2; ++d) {
				// shrink by half length so that any point inside this cube will be covered; add an extra factor of 1.5 to avoid table/chair/desk legs, etc.;
				// make sure safe_area is strictly normalized
				safe_area.expand_in_dim(d, -min(1.5f*hlength, 0.49f*safe_area.get_sz_dim(d)));
				cand_dest[d] = my_rgen.rand_uniform(safe_area.d[d][0], safe_area.d[d][1]);
			}
			float const dist(p2p_dist(p1, cand_dest));
			
			if (dist < dist_thresh) { // already at this location
				if (check_line_coll_expand(p1, cand_dest, coll_radius, squish_hheight)) {update_path = 1; continue;} // location is invalid, need to update the path below
				has_fear_dest = 1; // it's valid, so stay there
				rat.speed     = 0.0;
				break;
			}
			float side_cov(0.5f*min(hide_area.dx(), hide_area.dy()) - hlength); // amount of overhang of the object around the rat's extents
			float const dist_to_fear(p2p_dist(rat.fear_pos, cand_dest));
			float score(side_cov - 0.5f*top_gap + 0.2f*dist_to_fear - 0.1f*max(dist, dist_thresh)); // can be positive or negative
			if (best_score != 0.0 && score <= best_score) continue; // check score before iterating over other rats; it can only decrease below
			if (check_line_coll_expand(p1, cand_dest, coll_radius, squish_hheight)) continue; // use center before checking other rats so that the entire path is valid
			float tot_mdist(0.0);
			bool skip(0);
			float const radius_scale = 0.8; // smaller dist (head can overlap tail)
			vect_rat_t const &rats(interior->room_geom->rats);
			float const rsum_max(radius_scale*(rat.radius + rats.max_radius) + rats.max_xmove), coll_x1(cand_dest.x - rsum_max), coll_x2(cand_dest.x + rsum_max);
			auto start(rats.get_first_with_xv_gt(coll_x1)); // use a binary search to speed up iteration

			for (auto r = start; r != rats.end(); ++r) {
				if (r->pos.x > coll_x2) break; // no rat after this can overlap - done
				if (&(*r) == &rat) continue; // skip ourself
				float const r_sum(radius_scale*(rat.radius + r->radius)); // smaller dist (head can overlap tail)
				if (!dist_xy_less_than(cand_dest, r->pos, r_sum)) continue; // no rat in this spot
				float const move_dist(1.01*r_sum - p2p_dist_xy(cand_dest, r->pos)); // slightly larger than r_sum to prevent collisions
				cand_dest += (p1 - cand_dest).get_norm()*move_dist; // move our target in front of this other rat
				side_cov  -= move_dist; // moving to this misaligned position loses side coverage
				score      = 4.0*side_cov - 0.5f*top_gap + 0.25f*dist_to_fear - 0.1*max(dist, dist_thresh); // update score
				score     -= 0.2*dist; // less desirable when occupied
				score     -= 2.0*move_dist; // use prev accumulated move dist; even less desirable if there are many rats in the way
				tot_mdist += move_dist;
				if (tot_mdist > 4.0*rat.radius || !hide_area.contains_pt_xy(cand_dest)) {skip = 1; break;} // moved too far, too many other rats at this dest, skip it
				if (best_score != 0.0 && score <= best_score) {skip = 1; break;} // score dropped too low
				r = start-1; // go back and test the other rats for collisions with this new position; maybe creates a bit more determinism/less chaos
			} // for r
			if (skip) continue;
			if (tot_mdist > 0.0 && !is_pos_inside_building(cand_dest, xy_pad, hheight)) continue; // check if outside the valid area if center was moved
			best_dest  = point(cand_dest.x, cand_dest.y, rat.pos.z); // keep zval on the floor
			best_score = score;
			if (cand_dest.x == rat.dest.x && cand_dest.y == rat.dest.y) break; // keep the same dest (optimization)
		} // for c
		if (!has_fear_dest && best_score != 0.0) { // found a valid hiding place; score can be positive or negative
			rat.dest = best_dest;
			if (dist_less_than(rat.pos, rat.dest, dist_thresh)) {rat.speed = 0.0;} // close enough - stop
			else {rat.speed = RAT_FEAR_SPEED*global_building_params.rat_speed;} // high speed if not at dest
			has_fear_dest = 1; // set to avoid triggering the code below if close to dest
			assert(rat.pos.z == rat.dest.z);
		}
		rat.fear = max(0.0f, (rat.fear - 0.2f*(timestep/TICKS_PER_SECOND))); // reduce fear over 5s
	}
	bool const is_at_dest(dist_less_than(rat.pos, rat.dest, dist_thresh));

	if (!is_scared && !rat.is_sleeping() && is_at_dest && rat.dist_since_sleep > 1.5*floor_spacing && rgen.rand_bool()) { // 50% chance of taking a rest
		rat.sleep_for(0.0, 4.0); // 0-4s
		rat.speed = 0.0; // will reset anim_time in the next frame
	}
	else if (!has_fear_dest && !rat.is_sleeping() &&
		(rat.speed == 0.0 || newly_scared || update_path || is_at_dest || check_line_coll_expand(rat.pos, rat.dest, coll_radius, hheight)))
	{
		// stopped, no dest, at dest, collided, or newly scared - choose a new dest
		float target_fov_dp(RAT_FOV_DP), target_max_dist(view_dist); // start at nominal/max values
		float const dist_upper_bound(0.12 + 0.88*(1.0 - rat.fear)); // shorten the distance based on the amount of fear to evade more easily
		float const min_step(min(dist_thresh, 0.05f*rat.radius));
		rat.speed = 0.0; // stop until we've found a valid destination

		for (unsigned n = 0; n < 200; ++n) { // make 200 tries
			if (n > 50) { // we've been at this for a while, maybe we need to relax our constraints, maybe follow the walls?
				target_fov_dp   -= 0.02; // allow for turns outside our field of view
				target_max_dist *= 0.96;  // decrease the max distance considered
			}
			vector3d vdir(dir_gen.gen_dir()); // random XY direction

			if (coll_dir != zero_vector) { // resolve the collision; target_fov_dp is ignored in this case
				try_resolve_coll(rat.dir, plus_z, coll_dir, vdir, (n <= 10)); // if earlier in the iteration, try moving in a tangent
			}
			else { // not colliding; check if the new direction is close enough to our current direction
				float dp(dot_product(rat.dir, vdir));
				if (n < 180 && dp < 0.0) {vdir.negate(); dp = -dp;} // only allow switching directions in the last 20 iterations; can still fail the test below
				if (dp < target_fov_dp) continue; // not in field of view, use a new direction
			}
			if (is_scared && n <= 100 && dot_product(dir_to_fear, vdir) > 0.0) continue; // don't move toward danger; may make the rat move into a corner
			float dist(rgen.rand_uniform(0.1, dist_upper_bound)*target_max_dist); // random distance out to max view dist
			max_eq(dist, min_step); // make sure distance isn't too short
			point const cand(rat.pos + dist*vdir);
			if (!is_pos_inside_building(cand, xy_pad, hheight)) continue; // check if outside the valid area
			point const p2(cand + line_project_dist*vdir + center_dz); // extend in vdir so that the head doesn't collide
			point const p1_ext(p1 + coll_radius*vdir); // move the line slightly toward the dest to prevent collisions at the initial location
			if (check_line_coll_expand(p1_ext, p2, coll_radius, squish_hheight)) continue;
			rat.dest  = cand;
			rat.speed = global_building_params.rat_speed*rgen.rand_uniform(0.5, 1.0)*(is_scared ? 1.5 : 1.0); // random speed
			break; // success
		} // for n
		assert(rat.pos.z == rat.dest.z);
	}
	// update direction
	vector3d new_dir;

	if (!dist_less_than(rat.pos, rat.dest, dist_thresh)) {
		new_dir = (rat.dest - rat.pos).get_norm(); // point toward our destination
	}
	else if (has_fear_dest) { // stop, rest, and point toward what we fear
		max_eq(max_xmove, fabs(rat.pos.x - rat.dest.x));
		new_dir   = dir_to_fear;
		rat.speed = rat.dist_since_sleep = 0.0;
		rat.pos   = rat.dest; // just move it to the dest to prevent instability
	}
	// else dir is unchanged
	rat.is_hiding = (has_fear_dest && dist_less_than(rat.pos, rat.dest, 2.0*dist_thresh)); // close to fear_dest
	float const turn_rate(is_scared ? 1.1 : 1.0); // higher turning rate when scared
	update_dir_incremental(rat.dir, new_dir, turn_rate, timestep, rgen);
	assert(rat.dir.z == 0.0); // must be in XY plane
}

void building_t::scare_rat(rat_t &rat, point const &camera_bs, int ped_ix) const {
	// Note: later calls to scare_rat_at_pos() have priority and will set rat.fear_pos, but all calls will accumulate fear
	float const sight_scare_amt = 0.5;
	vect_cube_t ped_bcubes;
	if (ped_ix >= 0) {get_ped_bcubes_for_building(ped_ix, ped_bcubes, 1);} // moving_only=1

	for (cube_t const &c : ped_bcubes) { // only the cube center is needed
		scare_rat_at_pos(rat, point(c.xc(), c.yc(), c.z1()), sight_scare_amt, 1); // other people in the building scare the rats
	}
	rat.near_player = 0;

	if (camera_surf_collide) {
		if (global_building_params.min_attack_rats > 0 && in_building_gameplay_mode()) { // rat attacks are enabled in gameplay mode
			// determine if the player is close and visible for attack strength; can't use a return value of scare_rat_at_pos() due to early termination
			if (fabs(rat.pos.z - camera_bs.z) < get_window_vspace()) { // same floor
				if (dist_less_than(rat.pos, camera_bs, RAT_VIEW_FLOORS*get_window_vspace())) { // close enough; doesn't have to be in the same room
					if (check_line_of_sight_large_objs(rat.pos, camera_bs)) {rat.near_player = 1;}
				}
			}
		}
		scare_rat_at_pos(rat, camera_bs, sight_scare_amt, 1); // the sight of the player walking in the building scares the rats
	}
	sphere_t const cur_sound(get_cur_frame_loudest_sound());
	if (cur_sound.radius > 0.0) {scare_rat_at_pos(rat, cur_sound.pos, 4.0*cur_sound.radius, 0);}
}

void building_t::scare_rat_at_pos(rat_t &rat, point const &scare_pos, float amount, bool by_sight) const {
	assert(amount > 0.0);
	if (fabs(rat.pos.z - scare_pos.z) > get_window_vspace()) return; // on a different floor, ignore
	if (rat.fear > 0.99 && dist_less_than(rat.fear_pos, scare_pos, rat.radius)) return; // already max fearful of this location (optimization)
	point const pos(rat.get_center()); // use center zval, not floor zval
	int const scare_room(get_room_containing_pt(scare_pos)), rat_room(get_room_containing_pt(pos));
	//assert(rat_room >= 0); // this generally doesn't fail, unless rats collide with each other and push each other outside of a room
	if (rat_room != scare_room) {amount *= 0.67;} // less fearful if in a different room
	float const max_scare_dist(RAT_VIEW_FLOORS*get_window_vspace()), scare_dist(max_scare_dist*min(amount, 1.0f));
	float const fear((scare_dist - p2p_dist(pos, scare_pos))/max_scare_dist);
	if (fear <= 0.0) return;
	if (by_sight && !check_line_of_sight_large_objs(pos, scare_pos)) return; // check line of sight
	rat.fear     = min(1.0f, (rat.fear + fear));
	rat.fear_pos = scare_pos;
}


// *** Spiders ***

// Note: radius is really used for height and body size; legs can extend out to ~2x radius, so we decrease radius to half the user-specified value to account for this
spider_t::spider_t(point const &pos_, float radius_, vector3d const &dir_) : building_animal_t(pos_, 0.5*radius_, dir_), upv(plus_z), update_time(0.0) {
	pos.z += radius; // shift upward so that the center is off the ground
}
cube_t spider_t::get_bcube() const {
	cube_t bcube(pos);
	bcube.expand_by(vector3d(2.0*radius, 2.0*radius, radius)); // conservative?
	return bcube;
}
void spider_t::choose_new_dir(rand_gen_t &rgen) {
	if (upv == zero_vector) {upv = plus_z;} // can this ever happen?
	dir = cross_product(rgen.signed_rand_vector(), upv).get_norm(); // must be orthogonal to upv
}

void building_t::update_spiders(point const &camera_bs, unsigned building_ix, int ped_ix) {
	if (global_building_params.num_spiders_max == 0) return;
	vect_spider_t &spiders(interior->room_geom->spiders);
	if (spiders.placed && spiders.empty()) return; // no spiders placed in this building
	//timer_t timer("Update Spiders");
	add_animals_on_floor(spiders, building_ix, global_building_params.num_spiders_min, global_building_params.num_spiders_max,
		global_building_params.spider_size_min, global_building_params.spider_size_max);
	// update spiders
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	for (spider_t &spider : spiders) {spider.move(timestep);}
	spiders.do_sort();
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame
	for (spider_t &spider : spiders) {update_spider(spider, camera_bs, timestep, spiders.max_xmove, rgen);}
}

class surface_orienter_t {
	point pos;
	float r_inner, r_outer;
	vect_cube_t colliders; // share across calls?
public:
	surface_orienter_t(point const &pos_, float rxy, float rz) : pos(pos_), r_inner(min(rxy, rz)), r_outer(max(rxy, rz)) {}

	void register_cube(cube_t const &c, unsigned enable_faces) {
		if (c.contains_pt_exp(pos, r_outer)) {colliders.push_back(c);} // record for later processing
	}
	void register_cubes(vect_cube_t const &cubes, unsigned enable_faces) {
		for (cube_t const &c : cubes) {register_cube(c, enable_faces);}
	}
	bool align_to_surfaces(point &pos, vector3d &forward, vector3d &up, point const &last_pos, float hheight, float speed, float timestep, point const &camera_bs, rand_gen_t &rgen) {
		if (colliders.empty()) return 0;
		// Note: assumes last_pos is valid and non-intersecting; may not hold for initial placement or when objects are moved
		float const dist(p2p_dist(pos, last_pos));
		unsigned const NUM_DIRS = 64; // generate this many random directions
		vector3d best_dir(forward), best_up(up);
		point best_pos(pos);
		float best_score(0.0); // distance squared

		// generate a number of possible new vectors in the forward direction, perform sphere-cube intersection with each one, and choose the new location closest to the target
		for (unsigned n = 0; n <= NUM_DIRS; ++n) { // include forward on first iteration
			vector3d dir;

			if (n == 0) {dir = forward;}
			else { // maybe this should be more consistent, such as all multiples of 45 degrees from "forward"
				dir = rgen.signed_rand_vector().get_norm();
				if (dot_product(dir, forward) < 0.0) {dir.negate();} // face forward
			}
			point cand(last_pos + dist*dir);
			//if (dmin_sq > 0.0 && p2p_dist_sq(pos, cand) >= dmin_sq) continue; // not closer, skip
			bool had_coll(0);
			for (auto const &c : colliders) {had_coll |= sphere_cube_intersect(cand, r_outer, c);}
			if (!had_coll) continue; // floating in space, not a valid dir
			had_coll = 0;
			vector3d coll_norm(up);
			for (auto const &c : colliders) {had_coll |= sphere_cube_int_update_pos(cand, r_inner, c, last_pos, 1, 0, &coll_norm);}

			if (n == 0 && !had_coll) { // forward vector, no coll
				// rerun with outer radius to ensure coll_norm is correct
				for (auto const &c : colliders) {point tmp(cand); sphere_cube_int_update_pos(tmp, r_outer, c, last_pos, 1, 0, &coll_norm);}
				//return; // forward direction is still valid, done (early termination optimization)
			}
			float const score(p2p_dist_sq(cand, last_pos)*dot_product(dir, forward)); // squared movement distance, with higher weight for moving in the correct direction
			if (score > best_score) {best_dir = dir; best_up = coll_norm; best_pos = cand; best_score = score;}
		} // for n
		if (best_score == 0.0) return 0; // no valid directions, must be floating in space
		//float const delta_dir(min(1.0f, 1.5f*(1.0f - pow(0.7f, timestep))));
		//forward = ((delta_dir*best_dir + (1.0 - delta_dir)*forward)).get_norm(); // slowly reorient
		forward = best_dir;
		up      = best_up;
		pos     = best_pos;
		orthogonalize_dir(forward, up, forward, 1);
		// TODO: stuck against each other
		return 1;
	}
}; // surface_orienter_t

class obj_avoid_t {
	point &pos;
	point const &p_last;
	float radius;
public:
	bool had_coll;
	obj_avoid_t(point &pos_, point const &p_last_, float radius_) : pos(pos_), p_last(p_last_), radius(radius_), had_coll(0) {}
	void register_avoid_cube(cube_t const &c) {had_coll |= sphere_cube_int_update_pos(pos, radius, c, p_last);}
};

void building_t::update_spider_pos_orient(spider_t &spider, point const &camera_bs, float timestep, bool &on_surface, bool &had_coll, rand_gen_t &rgen) const {
	assert(interior);
	unsigned const EF_XY12(EF_X12 | EF_Y12);
	float const trim_thickness(get_trim_thickness());
	float const coll_radius(2.0f*spider.radius);
	surface_orienter_t surface_orienter(spider.pos, spider.radius, coll_radius);
	// Note: we can almost use fc_occluders, except this doesn't contain the very bottom floor because it's not an occluder
	//surface_orienter.register_cubes(interior->fc_occluders, EF_Z12); // Z surface only
	surface_orienter.register_cubes(interior->floors,   EF_Z2); // Z2 surface only
	surface_orienter.register_cubes(interior->ceilings, EF_Z1); // Z1 surface only
	for (unsigned d = 0; d < 2; ++d) {surface_orienter.register_cubes(interior->walls[d], EF_XY12);} // XY walls
	// check doors
	cube_t tc(spider.pos);
	tc.expand_by_xy(coll_radius);
	tc.expand_in_dim(2, 1.5*spider.radius); // smaller expand in Z
	obj_avoid_t obj_avoid(spider.pos, spider.last_pos, coll_radius);

	for (auto const &ds : interior->door_stacks) {
		if (ds.z1() > tc.z2() || ds.z2() < tc.z1()) continue; // wrong floor for this stack/part
		// calculate door bounds for bcube test, assuming it's open
		cube_t door_bounds(ds);
		door_bounds.expand_by_xy(ds.get_width());
		if (!tc.intersects(door_bounds)) continue; // optimization
		assert(ds.first_door_ix < interior->doors.size());

		for (unsigned dix = ds.first_door_ix; dix < interior->doors.size(); ++dix) {
			door_t const &door(interior->doors[dix]);
			if (!ds.is_same_stack(door)) break; // moved to a different stack, done
			if (door.z1() > tc.z2() || door.z2() < tc.z1()) continue; // wrong floor

			if (door.open) { // how to handle open doors? they're not cubes; avoid them entirely? use their bcubes?
				tquad_with_ix_t const door_tq(set_interior_door_from_cube(door));
				cube_t tight_door_bounds(door_tq.get_bcube()); // somewhat more accurate
				tight_door_bounds.expand_by_xy(door.get_thickness()); // conservative
				if (tc.intersects(tight_door_bounds)) {obj_avoid.register_avoid_cube(tight_door_bounds);}
			}
			else if (tc.intersects(door)) { // closed door
				surface_orienter.register_cube(door, (door.dim ? EF_Y12 : EF_X12));
			}
		}
	} // for door_stacks
	// check interior objects
	static vect_cube_t cubes, avoid;
	vect_room_object_t::const_iterator b, e;
	get_begin_end_room_objs_on_ground_floor(tc.z2(), 1, b, e); // for_spider=1

	for (auto i = b; i != e; ++i) {
		if (i->z1() > tc.z2() || i->z2() < tc.z1()) continue;
		if (!tc.intersects((i->type == TYPE_CLOSET) ? get_true_room_obj_bcube(*i) : *i)) continue; // no intersection with this object
		if (!i->is_spider_collidable()) continue;
		if (i->get_max_extent() < spider.radius) continue; // too small, skip
		get_room_obj_cubes(*i, spider.pos, cubes, avoid, avoid); // climb on large objects and avoid small and non-cube objects
		unsigned faces(EF_XY12 | EF_Z2); // all sides except for Z1
		if (i->type == TYPE_LIGHT || i->type == TYPE_MIRROR) {faces |= EF_Z1;} // hanging objects
		surface_orienter.register_cubes(cubes, faces);
		for (cube_t const &c : avoid) {obj_avoid.register_avoid_cube(c);}
		cubes.clear();
		avoid.clear();
	} // for i
	// check elevators
	for (elevator_t const &e : interior->elevators) {surface_orienter.register_cube(e, EF_XY12);} // XY surfaces; should we avoid open elevators?
	// check exterior walls; exterior doors are ignored for now (meaning the spider can walk on them)
	for (auto i = parts.begin(); i != get_real_parts_end_inc_sec(); ++i) {
		for (unsigned dim = 0; dim < 2; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				cube_t cube(*i);
				cube.d[dim][!dir] = i->d[dim][dir] + (dir ? -1.0 : 1.0)*trim_thickness; // shrink to trim thickness
				cubes.clear();
				cubes.push_back(cube); // start with entire length

				for (auto j = parts.begin(); j != get_real_parts_end(); ++j) { // clip against other parts
					if (j == i) continue; // skip self
					subtract_cube_from_cubes(*j, cubes); // subtract this part from current cubes by clipping in XY
				}
				unsigned const face(dim ? (dir ? EF_Y1 : EF_Y2) : (dir ? EF_X1 : EF_X2));
				for (cube_t const &c : cubes) {surface_orienter.register_cube(c, face);}
			} // for dir
		} // for dim
	} // for i
	on_surface = surface_orienter.align_to_surfaces(spider.pos, spider.dir, spider.upv, spider.last_pos, spider.radius, spider.speed, timestep, camera_bs, rgen);
	had_coll   = obj_avoid.had_coll;
}

void building_t::update_spider(spider_t &spider, point const &camera_bs, float timestep, float &max_xmove, rand_gen_t &rgen) const {
	
	if (spider.is_sleeping()) return; // peacefully sleeping
	rgen.rand_mix(); // make sure it's different per spider
	// Note: we use a small xy_pad to allow the spider to climb on exterior walls
	float const radius(spider.radius), height(2.0*radius), coll_radius(2.0f*radius), xy_pad(radius/*coll_radius + get_trim_thickness()*/);
	if (spider.speed == 0.0) {spider.speed = global_building_params.spider_speed*rgen.rand_uniform(0.5, 1.0);} // random speed

	if (!is_pos_inside_building(spider.pos, xy_pad, radius)) {
		spider.pos = spider.last_pos; // restore previous pos before collision
		spider.dir = zero_vector; // will be set to a valid value on the next frame
		return;
	}
	bool on_surface(0), had_coll(0);
	update_spider_pos_orient(spider, camera_bs, timestep, on_surface, had_coll, rgen);

	if (had_coll || spider.dir.mag() < 0.5) {spider.choose_new_dir(rgen);} // regenerate dir if collided, zero, or otherwise bad
	else if (on_surface && (float)tfticks > spider.update_time) { // direction change or sleep
		if (spider.dist_since_sleep > 2.0*get_window_vspace() && rgen.rand_bool()) { // 50% chance of taking a rest
			spider.sleep_for(0.0, 4.0); // 0-4s
			spider.speed = 0.0; // will reset anim_time in the next frame
		}
		else {
			spider.update_time = (float)tfticks + rand_uniform(4.0, 10.0)*TICKS_PER_SECOND; // 4-10s
			vector3d const prev_dir(spider.dir);
			spider.choose_new_dir(rgen);
			spider.dir = (spider.dir + prev_dir).get_norm(); // 50% mix of prev and new dir to avoid sharp turns
		}
	}
	vector3d coll_dir;
	point const prev_pos(spider.pos); // capture the pre-collision point

	if (check_and_handle_dynamic_obj_coll(spider.pos, coll_radius, (spider.pos.z - height), (spider.pos.z + height), camera_bs, 1)) { // check for collisions; for_spider=1
		coll_dir = (prev_pos - spider.pos).get_norm(); // points toward the collider in the XY plane

		// check if new pos is valid, and has a path to dest
		if (!is_pos_inside_building(spider.pos, xy_pad, radius)) {
			spider.pos = prev_pos; // restore previous pos before collision
			spider.sleep_for(0.1, 0.2); // wait 0.1-0.2s so that we don't immediately collide and get pushed out again
		}
		else {
			max_eq(max_xmove, fabs(spider.pos.x - prev_pos.x));
			// update the path about every 30 frames of colliding; this prevents the spider from being stuck while also avoiding jittering due to frequent dest updates;
			// using randomness rather than updating every 30 frames makes this process less regular and mechanical
			if ((rgen.rand()%30) == 0) {
				vector3d const old_dir(spider.dir);
				spider.choose_new_dir(rgen);
				try_resolve_coll(old_dir, spider.upv, coll_dir, spider.dir, rgen.rand_bool()); // try_tangent=random
			}
		}
	}
	spider.dest = spider.pos + radius*spider.dir; // put our destination in front of us
	//update_dir_incremental(spider.dir, new_dir, 1.0, timestep, rgen);

	// handle biting the player
	if (!in_building_gameplay_mode()) return;
	if (dot_product_ptv(spider.dir, camera_bs, spider.pos) < 0.0) return; // facing the wrong direction
	if (get_floor_for_zval(camera_bs.z) != get_floor_for_zval(spider.pos.z)) return; // wrong floor
	float const player_radius(CAMERA_RADIUS*global_building_params.player_coll_radius_scale), min_dist(player_radius + coll_radius);

	if (dist_xy_less_than(spider.pos, camera_bs, (player_radius + coll_radius))) { // do damage when nearly colliding with the player
		bool const played_sound(play_attack_sound(local_to_camera_space(spider.pos), 0.5, 1.5, rgen)); // quieter and higher pitch than rats
		if (played_sound) {player_take_damage(0.1, 1);} // poisoned every so often
	}
}

