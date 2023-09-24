// 3D World - City Birds Implementation
// by Frank Gennari
// 09/21/23

#include "city_objects.h"

extern int animate2;
extern float fticks;
extern object_model_loader_t building_obj_model_loader;

enum {BIRD_STATE_FLYING=0, BIRD_STATE_GLIDING, BIRD_STATE_LANDING, BIRD_STATE_STANDING, BIRD_STATE_TAKEOFF};


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

void city_bird_t::next_frame(float timestep) {
	// TODO: update velocity, dir, pos, etc.
	// update state
	// choose dest
	anim_time += timestep;
}

// this is here because only birds are updated each frame
void city_obj_placer_t::next_frame() {
	if (!animate2 || birds.empty()) return;
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	for (city_bird_t &bird : birds) {bird.next_frame(timestep);}
	// update bird_groups; is there a more efficient way than rebuilding bird_groups each frame?
	bird_groups.clear();
	for (unsigned i = 0; i < birds.size(); ++i) {bird_groups.insert_obj_ix(birds[i].bcube, i);}
	bird_groups.create_groups(birds, all_objs_bcube);
}

