// 3D World - City Birds Implementation
// by Frank Gennari
// 09/21/23

#include "city_objects.h"

extern int animate2;
extern float fticks;
extern double tfticks;
extern object_model_loader_t building_obj_model_loader;

enum {BIRD_STATE_FLYING=0, BIRD_STATE_GLIDING, BIRD_STATE_LANDING, BIRD_STATE_STANDING, BIRD_STATE_TAKEOFF, NUM_BIRD_STATES};


city_bird_t::city_bird_t(point const &pos_, float height, vector3d const &init_dir, rand_gen_t &rgen) : city_bird_base_t(pos_, height, init_dir, OBJ_MODEL_BIRD_ANIM) {
	state     = BIRD_STATE_STANDING;
	anim_time = 1.0*TICKS_PER_SECOND*rgen.rand_float(); // 1s random variation so that birds aren't all in sync
}

void city_bird_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	// animations: 0=flying, 1=gliding, 2=landing, 3=standing, 4=takeoff
	animation_state_t anim_state(1, ANIM_ID_SKELETAL, 0.025*anim_time/TICKS_PER_SECOND, get_model_anim_id()); // enabled=1
	building_obj_model_loader.draw_model(dstate.s, pos, bcube, dir, WHITE, dstate.xlate, OBJ_MODEL_BIRD_ANIM, shadow_only, 0, &anim_state);
}
/*static*/ void city_bird_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	animation_state_t anim_state(1); // enabled=1
	anim_state.clear_animation_id(dstate.s); // clear animations
	model3d::bind_default_flat_normal_map();
}

void city_bird_t::set_takeoff_time(rand_gen_t &rgen) {
	takeoff_time = tfticks + rgen.rand_uniform(10.0, 30.0)*TICKS_PER_SECOND; // wait 10-30s
}

void city_bird_t::next_frame(float timestep, bool &tile_changed, rand_gen_t &rgen) {
	// update state
	uint8_t const prev_state(state);
	float const new_anim_time(anim_time + timestep);

	switch (state) {
	case BIRD_STATE_FLYING:
		if (velocity.z < 0.0 && is_anim_cycle_complete(new_anim_time)) {state = BIRD_STATE_GLIDING;} // maybe switch to gliding
		// fall through
	case BIRD_STATE_GLIDING:
		if (velocity.z > 0.0) {state = BIRD_STATE_FLYING;} // should we call is_anim_cycle_complete()?
		// check if close enough to dest, then switch to landing
		//state = BIRD_STATE_LANDING;
		break;
	case BIRD_STATE_LANDING:
		if (is_anim_cycle_complete(new_anim_time)) { // wait for animation to complete
			state    = BIRD_STATE_STANDING;
			velocity = zero_vector; // stopped
			set_takeoff_time(rgen);
		}
		break;
	case BIRD_STATE_STANDING:
		if (takeoff_time > tfticks) {
			state = BIRD_STATE_TAKEOFF;
			// TODO: choose dest
		}
		break;
	case BIRD_STATE_TAKEOFF:
		if (is_anim_cycle_complete(new_anim_time)) {state = BIRD_STATE_FLYING;} // wait for animation to complete
		break;
	default: assert(0);
	}
	if (state != prev_state) {anim_time = 0.0;} // reset animation time on state change
	else {anim_time = new_anim_time;}

	if (state != BIRD_STATE_STANDING) {
		// TODO: update velocity
	}
	if (velocity != zero_vector) { // update direction and apply movement
		dir    = vector3d(velocity.x, velocity.y, 0.0).get_norm(); // always point in XY direction of velocity
		point const prev_pos(pos);
		pos   += velocity*timestep;
		bcube += (pos - bcube.get_cube_center()); // translate bcube to match - keep it centered on pos
		tile_changed |= (get_tile_id_containing_point_no_xyoff(pos) != get_tile_id_containing_point_no_xyoff(prev_pos));
	}
}

bool city_bird_t::is_anim_cycle_complete() const {
	return 0; // TODO
}

// this is here because only birds are updated each frame
void city_obj_placer_t::next_frame() {
	if (!animate2 || birds.empty()) return;
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	bool tile_changed(0);
	for (city_bird_t &bird : birds) {bird.next_frame(timestep, tile_changed, bird_rgen);}
	
	if (tile_changed) { // update bird_groups; is there a more efficient way than rebuilding bird_groups each frame?
		bird_groups.clear();
		for (unsigned i = 0; i < birds.size(); ++i) {bird_groups.insert_obj_ix(birds[i].bcube, i);}
		bird_groups.create_groups(birds, all_objs_bcube);
	}
}

void vect_bird_place_t::add_placement(cube_t const &obj, bool dim, bool dir, bool orient_dir, float spacing, rand_gen_t &rgen) {
	point pos(0.0, 0.0, obj.z2());
	pos[ dim] = obj.d[dim][dir];
	pos[!dim] = rgen.rand_uniform(obj.d[!dim][0]+spacing, obj.d[!dim][1]-spacing); // random position along the edge
	emplace_back(pos, dim, orient_dir);
}
void vect_bird_place_t::add_placement_rand_dim_dir(cube_t const &obj, float spacing, rand_gen_t &rgen) {
	bool const dir(rgen.rand_bool());
	add_placement(obj, rgen.rand_bool(), dir, dir, spacing, rgen);
}
void vect_bird_place_t::add_placement_top_center(cube_t const &obj, rand_gen_t &rgen) {
	emplace_back(cube_top_center(obj), rgen.rand_bool(), rgen.rand_bool());
}

