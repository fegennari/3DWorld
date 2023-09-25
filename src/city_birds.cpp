// 3D World - City Birds Implementation
// by Frank Gennari
// 09/21/23

#include "city_objects.h"

extern int animate2;
extern float fticks;
extern double tfticks;
extern object_model_loader_t building_obj_model_loader;

enum {BIRD_STATE_FLYING=0, BIRD_STATE_GLIDING, BIRD_STATE_LANDING, BIRD_STATE_STANDING, BIRD_STATE_TAKEOFF, NUM_BIRD_STATES};


city_bird_t::city_bird_t(point const &pos_, float height, vector3d const &init_dir) : city_bird_base_t(pos_, height, init_dir, OBJ_MODEL_BIRD_ANIM) {
	state = BIRD_STATE_STANDING;
}

void city_bird_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	// animations: 0=flying, 1=gliding, 2=landing, 3=standing, 4=takeoff
	bool const anim_enabled(1 || animate2);
	unsigned const anim_id(0), model_anim_id(state);
	animation_state_t anim_state(anim_enabled, anim_id, anim_time, model_anim_id);
	building_obj_model_loader.draw_model(dstate.s, pos, bcube, dir, WHITE, dstate.xlate, OBJ_MODEL_BIRD_ANIM, shadow_only, 0, &anim_state);
}

void city_bird_t::next_frame(float timestep, bool &tile_changed, rand_gen_t &rgen) {
	// update state
	uint8_t const prev_state(state);

	switch (state) {
	case BIRD_STATE_FLYING:
		if (velocity.z < 0.0) {
			//state = BIRD_STATE_GLIDING; // maybe switch
		}
		// fall through
	case BIRD_STATE_GLIDING:
		if (velocity.z > 0.0) {
			state = BIRD_STATE_FLYING;
		}
		// check if close enough to dest
		break;
	case BIRD_STATE_LANDING:
		if (is_anim_cycle_complete()) { // wait for animation to complete
			state = BIRD_STATE_STANDING;
			takeoff_time = tfticks + rgen.rand_uniform(10.0, 30.0)*TICKS_PER_SECOND; // wait 10-30s
		}
		break;
	case BIRD_STATE_STANDING:
		if (takeoff_time > tfticks) {
			state = BIRD_STATE_TAKEOFF;
			// TODO: choose dest
		}
		break;
	case BIRD_STATE_TAKEOFF:
		if (is_anim_cycle_complete()) {state = BIRD_STATE_FLYING;} // wait for animation to complete
		break;
	default: assert(0);
	}
	if (state != prev_state) {anim_time = 0.0;} // reset animation time on state change
	else {anim_time += timestep;}
	// TODO: update velocity, dir
	if (velocity != zero_vector) {dir = vector3d(velocity.x, velocity.y, 0.0).get_norm();} // always point in XY direction of velocity

	if (velocity != zero_vector) { // apply movement
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

