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
void apply_building_gravity(float &vz, float dt_ticks);
void apply_fc_cube_max_merge_xy(vect_cube_t &cubes);


void building_animal_t::sleep_for(float time_secs_min, float time_secs_max) {
	wake_time = (float)tfticks + rand_uniform(time_secs_min, time_secs_max)*TICKS_PER_SECOND;
	dist_since_sleep = 0.0; // reset the counter
}
float building_animal_t::move(float timestep, bool can_move_forward) { // returns movement distance
	// update animation time using position change; note that we can't just do the update in the rat movement code below because pos may be reset in case of collision
	anim_time += p2p_dist(pos, last_pos)/radius; // scale with size so that small rat/spider legs move faster
	last_pos   = pos;
	if (anim_time > 100000.0) {anim_time = 0.0;} // reset animation after awhile to avoid precision problems; should this be done when the player isn't looking?
	float move_dist(0.0);

	if (is_sleeping()) {
		if ((float)tfticks > wake_time) {wake_time = speed = 0.0;} // time to wake up
	}
	else if (speed == 0.0) {
		anim_time = 0.0; // reset animation to rest pos
	}
	else if (can_move_forward) { // apply movement
		move_dist = timestep*speed;
		pos = pos + move_dist*dir;
		dist_since_sleep += move_dist;
	}
	return move_dist;
}

void building_t::update_animals(point const &camera_bs, unsigned building_ix) {
	if (!animate2 || is_rotated() || !has_room_geom() || interior->rooms.empty()) return;
	update_rats   (camera_bs, building_ix);
	update_spiders(camera_bs, building_ix);
	update_snakes (camera_bs, building_ix);
}

template<typename T> void building_t::add_animals_on_floor(T &animals, unsigned building_ix, unsigned num_min, unsigned num_max, float sz_min, float sz_max) const {
	if (animals.placed) return; // already placed
	animals.placed = 1; // even if there were no animals placed
	if (num_max == 0) return; // none to place
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, mat_ix+1); // unique per building
	float const floor_spacing(get_window_vspace());
	unsigned const num(num_min + ((num_min == num_max) ? 0 : (rgen.rand() % (num_max - num_min + 1))));
	animals.reserve(num);

	for (unsigned n = 0; n < num; ++n) {
		// there's no error check for min_sz <= max_sz, so just use min_sz in that case
		float const sz_scale((sz_min >= sz_max) ? sz_min : rgen.rand_uniform(sz_min, sz_max)), radius(0.5f*floor_spacing*sz_scale);
		point const pos(gen_animal_floor_pos(radius, T::value_type::allow_in_attic(), rgen));
		if (pos == all_zeros) continue; // bad pos? skip this animal
		animals.add(typename T::value_type(pos, radius, rgen.signed_rand_vector_xy().get_norm(), n));
	}
}

point building_t::gen_animal_floor_pos(float radius, bool place_in_attic, rand_gen_t &rgen) const {
	for (unsigned n = 0; n < 100; ++n) { // make up to 100 tries
		cube_t place_area;

		if (place_in_attic && has_attic() && (rgen.rand()%10) == 0) { // 10% of the time in the attic
			place_area = get_attic_part();
			place_area.z1() = place_area.z2();
		}
		else {
			unsigned const room_ix(rgen.rand() % interior->rooms.size());
			room_t const &room(interior->rooms[room_ix]);
			if (room.z1() > ground_floor_z1) continue; // not on the ground floor or basement
			place_area = room; // will represent the usable floor area
			place_area.z1() += get_fc_thickness(); // on top of the floor
		}
		place_area.expand_by_xy(-(radius + get_wall_thickness()));
		if (min(place_area.dx(), place_area.dy()) < 4.0*radius) continue; // room too small (can happen for has_complex_floorplan office buildings)
		point const pos(gen_xy_pos_in_area(place_area, radius, rgen, place_area.z1()));
		if (is_valid_ai_placement(pos, radius)) {return pos;} // check room objects; start in the open, not under something
	} // for n
	return all_zeros; // failed
}

bool building_t::is_pos_inside_building(point const &pos, float xy_pad, float hheight, bool inc_attic) const {
	float bcube_pad(xy_pad);
	if (inc_attic && has_attic() && pos.z >= interior->attic_access.z2()) {bcube_pad += get_attic_beam_depth();} // add extra spacing for attic beams (approximate)
	if (!bcube.contains_pt_xy_exp(pos, -bcube_pad)) return 0; // check for end point inside building bcube
	cube_t req_area(pos, pos);
	req_area.expand_by_xy(xy_pad);
	req_area.z2() += hheight;
	return (is_cube_contained_in_parts(req_area) || (inc_attic && cube_in_attic(req_area)));
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

rat_t::rat_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_) : building_animal_t(pos_, radius_, dir_, id_), dest(pos) {
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
bool rat_t::is_facing_dest() const {
	return (dot_product((dest - pos).get_norm(), dir) > 0.75); // only move if we're facing our dest, to avoid walking through an object
}

bool building_t::add_rat(point const &pos, float hlength, vector3d const &dir, point const &placed_from) {
	if (!rat_t::allow_in_attic() && point_in_attic(pos)) return 0;
	point rat_pos(pos);
	if (!get_zval_of_floor(pos, hlength, rat_pos.z)) return 0; // place on the floor, skip if there's no floor here
	rat_t rat(rat_pos, hlength, vector3d(dir.x, dir.y, 0.0).get_norm(), interior->room_geom->rats.size()); // dir in XY plane
	
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

void building_t::update_rats(point const &camera_bs, unsigned building_ix) {
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
		rat.move(timestep, rat.is_facing_dest());
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
		for (auto r = rats.rbegin(); r != rats.rend(); ++r) {update_rat(*r, camera_bs, timestep, rats.max_xmove, can_attack_player, rgen);}
	}
	else { // forward iteration; ~0.004ms per rat
		for (auto r = rats. begin(); r != rats. end(); ++r) {update_rat(*r, camera_bs, timestep, rats.max_xmove, can_attack_player, rgen);}
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
	else if (c.type == TYPE_RCHAIR) {
		return 0; // not a hiding spot - yet
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

template<bool xy_only> class dir_gen_t {
	vector<vector3d> dirs;
	unsigned dir_ix;

	void gen_dirs() {
		rand_gen_t rgen;
		dirs.resize(1009); // make it a prime number
		for (auto &dir : dirs) {dir = (xy_only ? rgen.signed_rand_vector_spherical_xy() : rgen.signed_rand_vector_spherical()).get_norm();}
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
dir_gen_t<1> dir_gen_xy;
dir_gen_t<0> dir_gen_xyz;

void building_t::update_rat(rat_t &rat, point const &camera_bs, float timestep, float &max_xmove, bool can_attack_player, rand_gen_t &rgen) const {

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
	else if (check_and_handle_dynamic_obj_coll(rat.pos, rat.pos, rat.radius, rat.pos.z, (rat.pos.z + height), camera_bs, 0)) { // check for collisions; for_spider=0
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
	scare_rat(rat, camera_bs);
	rat.attacking = (rat.near_player && can_attack_player);
	if (rat.attacking) {rat.fear = 0.0;} // no fear when attacking
	bool const is_scared(rat.fear > 0.0), newly_scared(is_scared && !was_scared);

	// determine destination
	if (rat.attacking) {
		float const player_radius(get_scaled_player_radius()), min_dist(player_radius + hlength);
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
			vector3d vdir(dir_gen_xy.gen_dir()); // random XY direction

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

void building_t::scare_rat(rat_t &rat, point const &camera_bs) const {
	// Note: later calls to scare_rat_at_pos() have priority and will set rat.fear_pos, but all calls will accumulate fear
	assert(interior);
	float const sight_scare_amt = 0.5;
	rat.near_player = 0;

	for (person_t const &p : interior->people) { // other people in the building scare the rats
		if (p.is_waiting_or_stopped()) continue; // only scare if moving
		scare_rat_at_pos(rat, point(p.pos.x, p.pos.y, p.get_z1()), sight_scare_amt, 1);
	}
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
spider_t::spider_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_) : building_animal_t(pos_, 0.5*radius_, dir_, id_), upv(plus_z) {
	pos.z += radius; // shift upward so that the center is off the ground
}
cube_t spider_t::get_bcube() const {
	cube_t bcube(pos);
	bcube.expand_by(vector3d(2.0*radius, 2.0*radius, radius)); // conservative?
	return bcube;
}
vector3d spider_t::get_size() const {
	float const xy_radius(get_xy_radius());
	vector3d sz(xy_radius, xy_radius, xy_radius);
	sz[get_max_dim(upv)] = radius; // use height for primary up dir, which should at least work for walking on cube surfaces
	return sz;
}
void spider_t::choose_new_dir(rand_gen_t &rgen) {
	if (upv == zero_vector) {upv = plus_z;} // can this ever happen?
	dir = cross_product(rgen.signed_rand_vector(), upv).get_norm(); // must be orthogonal to upv
}
void spider_t::jump(float vel) {
	assert(vel > 0.0);
	jump_vel_z = vel; // overwrite if already set
	speed      = 0.5*vel; // set forward speed lower for a steeper jump; dir must be set (and should probably be in the XY plane)
}
void spider_t::end_jump() {
	if (is_jumping()) {jump_vel_z = speed = 0.0;}
}

bool building_room_geom_t::maybe_spawn_spider_in_drawer(room_object_t const &c, cube_t const &drawer, unsigned drawer_id, float floor_spacing, bool is_door) {
	if (global_building_params.spider_drawer_prob == 0.0) return 0; // no spiders
	if (!spider_t::allow_in_attic() && c.in_attic()) return 0; // no spiders in the attic
	rand_gen_t rgen;
	rgen.set_state((unsigned(c.obj_id) << 8)+drawer_id+1, c.room_id+1);
	if (rgen.rand_float() >= global_building_params.spider_drawer_prob) return 0; // no spider
	float radius(0.5f*floor_spacing*rgen.rand_uniform(global_building_params.spider_size_min, global_building_params.spider_size_max));
	min_eq(radius, rgen.rand_uniform(0.6, 0.9)*min(drawer.dz(), 0.5f*min(drawer.dx(), drawer.dy()))); // make sure it fits in the drawer
	vector3d dir;
	dir[ c.dim] = (c.dir ? 1.0 : -1.0); // face the outside of the drawer
	dir[!c.dim] = 0.25*rgen.signed_rand_float(); // not straight out
	dir.normalize();
	point const pos(drawer.xc(), drawer.yc(), drawer.z1());
	spiders.emplace_back(pos, radius, dir, spiders.size());
	if (!is_door) {spiders.back().jump(0.002*rgen.rand_uniform(1.0, 1.4));} // jump out of drawers
	return 1;
}

void building_t::update_spiders(point const &camera_bs, unsigned building_ix) {
	vect_spider_t &spiders(interior->room_geom->spiders);
	if (spiders.placed && spiders.empty()) return; // no spiders placed in this building
	//timer_t timer("Update Spiders");
	add_animals_on_floor(spiders, building_ix, global_building_params.num_spiders_min, global_building_params.num_spiders_max,
		global_building_params.spider_size_min, global_building_params.spider_size_max);
	// update spiders
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	
	for (spider_t &spider : spiders) {
		if (!spider.squished) {spider.move(timestep);}
	}
	spiders.do_sort();
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame
	for (spider_t &spider : spiders) {update_spider(spider, camera_bs, timestep, spiders.max_xmove, rgen);}
}

class surface_orienter_t {
	vector3d size;
	cube_t check_cube;
	vect_cube_t colliders;
public:
	void init(point const &p1, point const &p2, vector3d const &size_) {
		size = size_;
		check_cube.set_from_point(p1); // include both cur and last points
		check_cube.union_with_pt (p2);
		check_cube.expand_by(size.get_max_val()); // we could expand by size, but expanding by the max dim is more conservative and possibly better
		colliders.clear();
	}
	void register_cube(cube_t const &c, bool check_z_merge=0) {
		if (!c.intersects(check_cube)) return;
		
		if (check_z_merge) { // attempt to merge ceilings into floors
			for (cube_t &c2 : colliders) {
				if (c2.z1() == c.z2()) {c2.z1() = c.z1();return;}
			}
		}
		colliders.push_back(c); // record for later processing in align_to_surfaces()
	}
	void register_cubes(vect_cube_t const &cubes, bool check_z_merge=0) {
		for (cube_t const &c : cubes) {register_cube(c, check_z_merge);}
	}
	void clip_and_max_expand_cubes() {
		if (colliders.size() <= 1) return; // not needed
		for (cube_t &c : colliders) {c.intersect_with_cube_xy(check_cube);}
		apply_fc_cube_max_merge_xy(colliders);
		sort_and_unique(colliders);
	}
	bool align_to_surfaces(spider_t &s, float delta_dir, point const &camera_bs, rand_gen_t &rgen) {
		if (colliders.empty()) return 0; // floating in midair
		// Note: assumes last_pos is valid and non-intersecting; may not hold for initial placement or when objects are moved
		float const dist(p2p_dist(s.pos, s.last_pos)), r_inner(size.get_min_val()), r_outer(size.get_max_val());
		vector3d dir, best_dir(s.dir), best_up(s.upv);
		point best_pos(s.pos);
		float best_score(0.0); // distance squared

		// generate a number of possible new vectors in the forward direction, perform sphere-cube intersection with each one, and choose the new location closest to the target
		for (unsigned n = 0; n <= 50; ++n) { // include forward on first iteration
			if (n == 0) {dir = s.dir;}
			else { // maybe this should be more consistent, such as all multiples of 45 degrees from "forward"
				dir = dir_gen_xyz.gen_dir(); // random XYZ direction
				if (dot_product(dir, s.dir) < 0.0) {dir.negate();} // face forward
			}
			point cand(s.last_pos + dist*dir);
			bool had_coll(0);
			for (auto const &c : colliders) {had_coll |= ellipse_cube_intersect(cand, 1.05*size, c);} // slightly larger radius to pick up nearby objects
			if (!had_coll) continue; // floating in space, not a valid dir
			had_coll = 0;
			vector3d coll_norm(s.upv);
			for (auto const &c : colliders) {had_coll |= sphere_cube_int_update_pos(cand, r_inner, c, s.last_pos, 1, 0, &coll_norm);}

			if (n == 0 && !had_coll) { // forward vector, no coll
				// rerun with outer radius to ensure coll_norm is correct
				for (auto const &c : colliders) {point tmp(cand); sphere_cube_int_update_pos(tmp, r_outer, c, s.last_pos, 1, 0, &coll_norm);}
			}
			// squared movement distance, with higher weight for moving in the correct direction, and a preference for moving orthogonal in upv
			float const dir_dp(dot_product(dir, s.dir)), up_dp(fabs(dot_product(dir, s.upv)));
			float const score(p2p_dist_sq(cand, s.last_pos)*(dir_dp + 0.25*up_dp));
			if (score > best_score) {best_dir = dir; best_up = coll_norm; best_pos = cand; best_score = score;}
			if (n == 0 && !had_coll) break; // forward direction is still valid, done (early termination optimization)
		} // for n
		if (best_score == 0.0) { // no valid directions, must be floating in space
			if (colliders.size() > 1) return 0; // multiple colliders is too complex to handle here
			point const contact_pt(colliders.back().closest_pt(s.pos));
			best_up  = (s.pos - contact_pt).get_norm(); // point away from the collider
			best_pos = contact_pt + r_inner*best_up;
			orthogonalize_dir(best_dir, best_up, best_dir, 1);
			delta_dir *= 0.25; // update more smoothly for stability
		}
		vector3d const delta(best_pos - s.pos);
		float const delta_mag(delta.mag());
		s.pos = ((delta_mag > s.radius) ? (s.pos + (s.radius/delta_mag)*delta) : best_pos); // limit movement distance to radius
		s.dir = ((delta_dir*best_dir + (1.0 - delta_dir)*s.dir)).get_norm(); // slowly reorient
		s.upv = ((delta_dir*best_up  + (1.0 - delta_dir)*s.upv)).get_norm(); // slowly reorient
		orthogonalize_dir(s.dir, s.upv, s.dir, 1);
		return 1;
	}
}; // surface_orienter_t

surface_orienter_t surface_orienter; // reused across spiders

class obj_avoid_t {
	point &pos;
	point const &p_last;
	float radius;
public:
	bool had_coll;
	obj_avoid_t(point &pos_, point const &p_last_, float radius_) : pos(pos_), p_last(p_last_), radius(radius_), had_coll(0) {}
	void register_avoid_cube(cube_t const &c) {had_coll |= sphere_cube_int_update_pos(pos, radius, c, p_last);}
};

bool building_t::update_spider_pos_orient(spider_t &spider, point const &camera_bs, float timestep, rand_gen_t &rgen) const {
	assert(interior);
	float const trim_thickness(get_trim_thickness()), coll_radius(spider.get_xy_radius());
	vector3d const size(spider.get_size());
	bool const in_attic(point_in_attic(spider.pos));
	obj_avoid_t obj_avoid(spider.pos, spider.last_pos, coll_radius); // use xy_radius for all dims
	surface_orienter.init(spider.pos, spider.last_pos, size);
	// Note: we can almost use fc_occluders, except this doesn't contain the very bottom floor because it's not an occluder, and maybe the overlaps would cause problems
	surface_orienter.register_cubes(interior->floors);
	surface_orienter.register_cubes(interior->ceilings, 1); // check_z_merge=1
	surface_orienter.clip_and_max_expand_cubes(); // required to remove false splits between ceilings and floors
	
	if (!in_attic) {
		for (unsigned d = 0; d < 2; ++d) {surface_orienter.register_cubes(interior->walls[d]);} // XY walls
	}
	cube_t tc(spider.pos);
	tc.expand_by_xy(size); // use xy_radius for all dims; okay to be convervative
	tc.expand_in_dim(2, 1.5*spider.radius); // smaller expand in Z

	if (!in_attic) { // check doors
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
					cube_t door_bcube(get_door_bounding_cube(door));
					bool const dir(door.get_check_dirs()); // side of the door frame the door opens to
					door_bcube.d[!door.dim][!dir] += (dir ? 1.0 : -1.0)*spider.radius; // shift edge away from door frame to allow spider to walk on inside of the frame
				
					if (tc.intersects(door_bcube)) {
						obj_avoid.register_avoid_cube(door_bcube);
						// while the extruded polygon collision check below is more accurate, it can result in spiders getting stuck behind open doors
						//if (sphere_ext_poly_intersect(door_tq.pts, 4, normal, spider.pos, coll_radius, door.get_thickness(), 0.0)) {obj_avoid.had_coll = 1;}
					}
				}
				else if (tc.intersects(door)) {surface_orienter.register_cube(door);} // closed door
			}
		} // for door_stacks
	}
	// check interior objects
	static vect_cube_t cubes, avoid;
	vect_room_object_t::const_iterator b, e;
	get_begin_end_room_objs_on_ground_floor(tc.z2(), 1, b, e); // for_spider=1

	for (auto i = b; i != e; ++i) {
		if (i->z1() > tc.z2()) continue; // object is too high
		if (i->z2() < tc.z1() && !(i->type == TYPE_DESK && i->shape == SHAPE_TALL)) continue; // sigh, tall desks are special
		if (!tc.intersects_xy((i->type == TYPE_CLOSET) ? get_true_room_obj_bcube(*i) : (cube_t)*i)) continue; // no intersection with this object
		if (!i->is_spider_collidable()) continue;
		if (i->get_max_extent() < spider.radius) continue; // too small, skip
		get_room_obj_cubes(*i, spider.pos, cubes, cubes, avoid); // climb on large and small objects and avoid non-cube objects
		surface_orienter.register_cubes(cubes);
		for (cube_t const &c : avoid) {obj_avoid.register_avoid_cube(c);}
		cubes.clear();
		avoid.clear();
	} // for i
	// check elevators
	for (elevator_t const &e : interior->elevators) {surface_orienter.register_cube(e);} // should we avoid open elevators?
	
	if (!in_attic) { // check exterior walls; exterior doors are ignored for now (meaning the spider can walk on them)
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
					for (cube_t const &c : cubes) {surface_orienter.register_cube(c);}
				} // for dir
			} // for dim
		} // for i
	}
	else { // in_attic case
		// the attic roof is not a cube we can walk on and the beams aren't real objects;
		// also, it's not really possible to move from the attic floor to the roof without getting stuck in/on the corner
	}
	if (obj_avoid.had_coll) {spider.end_jump();}
	float const delta_dir(min(1.0f, 1.5f*(1.0f - pow(0.7f, timestep))));

	if (!surface_orienter.align_to_surfaces(spider, delta_dir, camera_bs, rgen)) { // not on a surface
		if (!spider.is_jumping()) { // if jumping, we continue the jump; otherwise, drop to a surface below
			if (!spider.on_web) {
				spider.web_start_zval = max(spider.pos.z, spider.last_pos.z) + spider.get_xy_radius();
				spider.on_web = 1;
			}
			spider.pos    = spider.last_pos;
			spider.dir    = (delta_dir*-plus_z + (1.0 - delta_dir)*spider.dir).get_norm(); // slowly reorient to -z
			spider.pos.z -= 0.5*timestep*spider.speed; // drop at half speed
			orthogonalize_dir(spider.upv, spider.dir, spider.upv, 1);
		}
	}
	else { // on a surface
		spider.end_jump();
		spider.on_web = 0;
	}
	return obj_avoid.had_coll;
}

void building_t::maybe_bite_and_poison_player(point const &pos, point const &camera_bs, vector3d const &dir, float coll_radius, float damage, bool poison, rand_gen_t &rgen) const {
	if (!in_building_gameplay_mode()) return;
	if (dot_product_ptv(dir, camera_bs, pos) < 0.0) return; // facing the wrong direction
	if (get_floor_for_zval(camera_bs.z) != get_floor_for_zval(pos.z)) return; // wrong floor

	if (dist_xy_less_than(pos, camera_bs, (get_scaled_player_radius() + coll_radius))) { // do damage when nearly colliding with the player
		bool const played_sound(play_attack_sound(local_to_camera_space(pos), 0.5, 1.5, rgen)); // quieter and higher pitch than rats
		if (played_sound) {player_take_damage(damage, poison);} // bite and maybe poisoned every so often
	}
}

void building_t::update_spider(spider_t &spider, point const &camera_bs, float timestep, float &max_xmove, rand_gen_t &rgen) const {
	
	if (spider.squished || spider.is_sleeping()) return; // horribly squished or peacefully sleeping
	float const radius(spider.radius), height(2.0*radius), coll_radius(2.0f*radius);
	float const spider_z1(spider.pos.z - height), spider_z2(spider.pos.z + height);

	if (spider.is_jumping()) { // jumping
		spider.jump_dist += p2p_dist_xy(spider.pos, spider.last_pos);
		spider.pos.z     += timestep*spider.jump_vel_z;
		apply_building_gravity(spider.jump_vel_z, timestep);
		// if we're still in the initial phase of our jump, skip collision and movement logic below to avoid initial coll with our starting object
		if (spider.jump_dist < 2.0*spider.radius && is_pos_inside_building(spider.pos, radius, radius)) return;
	}
	rgen.rand_mix(); // make sure it's different per spider
	// Note: we use a small xy_pad to allow the spider to climb on exterior walls
	if (spider.speed == 0.0) {spider.speed = global_building_params.spider_speed*rgen.rand_uniform(0.5, 1.0);} // random speed

	if (!is_pos_inside_building(spider.pos, radius, radius)) {
		spider.end_jump();
		spider.pos = spider.last_pos; // restore previous pos before collision

		if (!is_pos_inside_building(spider.pos, radius, radius)) { // still not valid
			if (spider.last_valid_pos == all_zeros) { // bad spawn pos - retry
				gen_animal_floor_pos(radius, spider_t::allow_in_attic(), rgen);
				return;
			}
			spider.pos = spider.last_valid_pos; // restore to prev frame pos
		}
		assert(is_pos_inside_building(spider.pos, radius, radius));
		spider.choose_new_dir(rgen);
		return; // or could continue below?
	}
	spider.last_valid_pos = spider.pos;
	bool const had_coll(update_spider_pos_orient(spider, camera_bs, timestep, rgen));
	if (had_coll) {spider.end_jump();}

	// regenerate dir if collided and not on a web, or if dir is somehow bad
	if ((had_coll && !spider.on_web) || spider.dir.mag() < 0.5) {spider.choose_new_dir(rgen);}
	else if ((float)tfticks > spider.update_time) { // direction change or sleep
		if (spider.on_web) {
			spider.update_time = (float)tfticks + 1.0*TICKS_PER_SECOND; // wait another 1s before updating
		}
		else if (spider.dist_since_sleep > 2.0*get_window_vspace() && rgen.rand_bool()) { // 50% chance of taking a rest
			spider.sleep_for(0.1, 5.0); // 0.1-5s
			spider.speed = 0.0; // will reset anim_time in the next frame
		}
		else {
			spider.update_time = (float)tfticks + rand_uniform(5.0, 15.0)*TICKS_PER_SECOND; // 5-15s
			vector3d const prev_dir(spider.dir);
			spider.choose_new_dir(rgen);
			spider.dir = (spider.dir + prev_dir).get_norm(); // 50% mix of prev and new dir to avoid sharp turns
		}
	}
	if (spider.on_web) { // rotate around the vertical axis
		bool const rotate_dir(spider.id & 1);
		float const angle(0.04*timestep*(rotate_dir ? 1.0 : -1.0));
		rotate_vector3d(plus_z, angle, spider.upv);
		spider.dir.normalize();
	}
	vector3d coll_dir;
	point const prev_pos(spider.pos); // capture the pre-collision point

	if (!spider.is_jumping() && check_and_handle_dynamic_obj_coll(spider.pos, spider.pos, coll_radius, spider_z1, spider_z2, camera_bs, 1)) { // check for collisions; for_spider=1
		spider.end_jump();
		coll_dir = (prev_pos - spider.pos).get_norm(); // points toward the collider in the XY plane

		// check if new pos is valid, and has a path to dest
		if (!is_pos_inside_building(spider.pos, radius, radius)) {
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
	maybe_bite_and_poison_player(spider.pos, camera_bs, spider.dir, coll_radius, 0.1, 1, rgen); // handle biting the player: 0.1 damage with poison
}

bool building_t::maybe_squish_spider(room_object_t const &obj) {
	assert(has_room_geom());
	bool any_squished(0);

	for (spider_t &spider : interior->room_geom->spiders) {
		if (spider.squished) continue; // already squished
		if (!!obj.contains_pt_xy(spider.pos) || !obj.intersects(spider.get_bcube())) continue;
		if (obj.get_size().get_max_val() < spider.get_xy_radius()) continue; // object is too small to squish this spider
		add_blood_decal(spider.pos, 1.5*spider.get_xy_radius());
		spider.pos.z -= 0.4*spider.get_height(); // move it near the ground since it will be drawn flattened
		any_squished  = spider.squished = 1;
		register_achievement("Splat the Spider");
	} // for spider
	return any_squished;
}


// *** Snakes ***

snake_t::snake_t(point const &pos_, float radius_, vector3d const &dir_, unsigned id_) : building_animal_t(pos_, radius_, dir_, id_) {
	rand_gen_t rgen;
	rgen.set_state(id+1, 3*id+7);
	rgen.rand_mix();
	length  = 2.0*radius; // input radius is half length
	radius *= 0.04;
	color   = WHITE*(1.0 - rgen.rand_float()*rgen.rand_float()); // random color shade, weighted toward lighter
	has_rattle = rgen.rand_bool();
	unsigned const NUM_SEGS = 20; // head + 18 segments + tail
	float const seg_length(length/NUM_SEGS);
	vector3d const seg_step(-seg_length*dir); // head -> tail
	segments.resize(NUM_SEGS, pos);
	for (unsigned n = 0; n < NUM_SEGS; ++n) {segments[n] += seg_step*(n - NUM_SEGS/2.0 + 0.5);} // set segment centers
}
float snake_t::get_xy_radius() const {
	float xy_radius_sq(0.0);
	for (point const &p : segments) {max_eq(xy_radius_sq, p2p_dist_xy_sq(pos, p));}
	return (sqrt(xy_radius_sq) + radius); // don't forget to add the body radius
}
cube_t snake_t::get_bcube() const {
	cube_t bcube;
	bcube.set_from_points(segments);
	bcube.expand_by_xy(radius);
	bcube.z2() += radius;
	return bcube;
}
// Note: seg_ix is a float so that we can specify fractional segments; actual segment index is offset by 0.5 for the segment center
float snake_t::get_seg_radius(float seg_ix) const {
	assert(!segments.empty());
	float const t(seg_ix/segments.size()); // varies from 0.0 at tip of nose to 1.0 at tip of tail
	
	if (t < 0.4) { // head
		float const v(t/0.4); // 0-1
		return 0.5*(1.0 + v*v)*radius; // 0.5 => 1.0
	}
	if (t > 0.6) { // tail
		float const v((t - 0.6)/0.4); // 0-1
		return sqrt(1.0 - v)*radius; // 1.0 => 0.0
	}
	return radius; // middle of body
}
void snake_t::move_segments(float dist) {
	if (dist == 0.0) return;
	assert(dist > 0.0);
	assert(segments.size() >= 2);
	float const seg_len(get_seg_length());
	segments[0] += dist*dir; // move the head forward incrementally

	// work backwards toward the tail and apply the new segment directions
	for (unsigned n = 1; n < segments.size(); ++n) {
		vector3d const seg_delta(segments[n-1] - segments[n]);
		segments[n] = segments[n-1] - seg_delta*(seg_len/seg_delta.mag());
	}
	// update pos as new bcube center; pos.z remains unchanged
	cube_t const bcube(get_bcube());
	pos.x = bcube.xc();
	pos.y = bcube.yc();
}

void building_t::update_snakes(point const &camera_bs, unsigned building_ix) {
	vect_snake_t &snakes(interior->room_geom->snakes);
	if (snakes.placed && snakes.empty()) return; // no snakes placed in this building
	//timer_t timer("Update Snakes");
	add_animals_on_floor(snakes, building_ix, global_building_params.num_snakes_min, global_building_params.num_snakes_max,
		global_building_params.snake_size_min, global_building_params.snake_size_max);
	// update snakes
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	snakes.do_sort(); // is this necessary?
	rand_gen_t rgen;
	rgen.set_state(building_ix+1, frame_counter+1); // unique per building and per frame
	for (snake_t &snake : snakes) {update_snake(snake, camera_bs, timestep, snakes.max_xmove, rgen);}
}

void building_t::update_snake(snake_t &snake, point const &camera_bs, float timestep, float &max_xmove, rand_gen_t &rgen) const {
	if (snake.speed == 0.0) {snake.speed = global_building_params.snake_speed*rgen.rand_uniform(0.5, 1.0);} // random speed
	float const radius(snake.radius), height(snake.get_height()), hheight(0.5*height);
	vector3d const center_dz(0.0, 0.0, hheight);
	point const &old_head_pos(snake.get_head_pos()); // only the head moves
	point head_pos(old_head_pos + snake.dir*(timestep*snake.speed)); // move the head one timestep
	point const pre_head_pos(head_pos);
	vector3d coll_dir;
	bool change_dir(0);
	
	if (snake.is_sleeping()) {} // peacefully sleeping, no collision needed
	else if (check_and_handle_dynamic_obj_coll(head_pos, snake.pos, radius, head_pos.z, (head_pos.z + height), camera_bs, 0)) { // check for collisions; for_spider=0
		coll_dir   = (pre_head_pos - head_pos).get_norm(); // points toward the collider in the XY plane
		change_dir = 1;
	}
	// check if pos is valid
	// TODO: should look ahead and avoid obstacles
	// TODO: return coll normal and use that to select a mew direction away from the object
	else if (!is_pos_inside_building(head_pos, radius, hheight)) {change_dir = 1;}
	else if (check_line_coll_expand((old_head_pos + center_dz), (head_pos + center_dz), radius, hheight)) {change_dir = 1;}
	else {max_eq(max_xmove, fabs(head_pos.x - old_head_pos.x));}
	
	if (change_dir) { // collision, change direction rather than moving
		//vector3d const new_dir_hemisphere((coll_dir == zero_vector) ? snake.dir : -coll_dir); // use coll_dir if available
		vector3d const new_dir_hemisphere(snake.dir);
		snake.dir = rgen.signed_rand_vector_xy().get_norm();
		// keep turn angle less than a bit more than 180 degrees
		// add a stuck counter that allows a sharp turn if snake has been stuck for too many frames? or just allow any angle every 64 frames?
		if (dot_product(snake.dir, new_dir_hemisphere) < -0.1 && (rgen.rand() & 63)) {snake.dir.negate();}
	}
	else { // move forward
		snake.last_valid_dir = snake.dir;
		maybe_bite_and_poison_player((head_pos + center_dz), camera_bs, snake.dir, 2.0*snake.radius, 0.5, snake.has_rattle, rgen); // 0.5 damage, poison if has a rattle
		float const move_dist(snake.move(timestep));
		snake.move_segments(move_dist);
		
		if (move_dist > 0.0) { // move head in a winding motion if moving
			vector3d const side_dir(cross_product(snake.dir, plus_z));
			float const speed_factor(snake.speed/global_building_params.snake_speed); // [0.5, 1.0]
			float const rot_amt(sin(0.1*snake.anim_time)); // rotation amount should be independent of speed
			rotate_vector3d(plus_z, 0.02*fticks*PI*rot_amt*speed_factor, snake.dir);
			snake.dir.normalize(); // is this needed?
		}
	}
}

