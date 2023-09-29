// 3D World - City Birds Implementation
// by Frank Gennari
// 09/21/23

#include "city_objects.h"

float const anim_time_scale(1.0/TICKS_PER_SECOND);

extern int animate2;
extern float fticks;
extern double tfticks;
extern object_model_loader_t building_obj_model_loader;

enum {BIRD_STATE_FLYING=0, BIRD_STATE_GLIDING, BIRD_STATE_LANDING, BIRD_STATE_STANDING, BIRD_STATE_TAKEOFF, NUM_BIRD_STATES};


city_bird_t::city_bird_t(point const &pos_, float height, vector3d const &init_dir, unsigned loc_ix_, rand_gen_t &rgen) :
	city_bird_base_t(pos_, height, init_dir, OBJ_MODEL_BIRD_ANIM), state(BIRD_STATE_STANDING), loc_ix(loc_ix_)
{
	anim_time = 1.0*TICKS_PER_SECOND*rgen.rand_float(); // 1s random variation so that birds aren't all in sync
}

void city_bird_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!dstate.check_cube_visible(bcube, dist_scale)) return;
	// animations: 0=flying, 1=gliding, 2=landing, 3=standing, 4=takeoff
	float const model_anim_time(anim_time_scale*anim_time/SKELETAL_ANIM_TIME_CONST); // divide by constant to cancel out multiply in draw_model()
	animation_state_t anim_state(1, ANIM_ID_SKELETAL, model_anim_time, get_model_anim_id()); // enabled=1
	building_obj_model_loader.draw_model(dstate.s, pos, bcube, dir, WHITE, dstate.xlate, OBJ_MODEL_BIRD_ANIM, shadow_only, 0, &anim_state);

	if (0 && dest_valid()) { // debug drawing
		post_draw(dstate, shadow_only); // clear animations
		vector<vert_color> pts;
		pts.emplace_back(pos,  YELLOW);
		pts.emplace_back(dest, YELLOW);
		select_texture(WHITE_TEX);
		draw_verts(pts, GL_LINES);
	}
}
/*static*/ void city_bird_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	animation_state_t anim_state(1); // enabled=1
	anim_state.clear_animation_id(dstate.s); // clear animations
	model3d::bind_default_flat_normal_map();
}

void city_bird_t::set_takeoff_time(rand_gen_t &rgen) {
	takeoff_time = tfticks + rgen.rand_uniform(10.0, 30.0)*TICKS_PER_SECOND; // wait 10-30s
}
bool city_bird_t::is_anim_cycle_complete(float new_anim_time) const {
	if (anim_time == 0.0) return 0; // anim_time was just reset
	return building_obj_model_loader.check_anim_wrapped(OBJ_MODEL_BIRD_ANIM, get_model_anim_id(), anim_time_scale*anim_time, anim_time_scale*new_anim_time);
}
bool city_bird_t::in_landing_dist() const {
	assert(dest_valid());
	assert(state == BIRD_STATE_FLYING || state == BIRD_STATE_GLIDING);
	if (velocity == zero_vector) return 0; // not moving
	float const land_delay_secs(building_obj_model_loader.get_anim_duration(OBJ_MODEL_BIRD_ANIM, get_model_anim_id()));
	float const land_delay_ticks(land_delay_secs/anim_time_scale);
	vector3d const delta(dest - pos);
	float const dest_time(delta.mag_sq()/dot_product(velocity, delta)); // distance/velocity_to_dest = delta.mag()/dot_product(velocity, delta/deta.mag())
	return (dest_time < land_delay_ticks);
}

void city_bird_t::next_frame(float timestep, bool &tile_changed, city_obj_placer_t &placer, rand_gen_t &rgen) { // timestep is in ticks
	// update state
	uint8_t const prev_state(state);
	float const new_anim_time(anim_time + timestep);

	switch (state) {
	case BIRD_STATE_FLYING:
		if (velocity.z < 0.0 && is_anim_cycle_complete(new_anim_time)) {state = BIRD_STATE_GLIDING;} // maybe switch to gliding
		// fall through
	case BIRD_STATE_GLIDING:
		if (velocity.z > 0.0)  {state = BIRD_STATE_FLYING ;} // should we call is_anim_cycle_complete()?
		if (in_landing_dist()) {state = BIRD_STATE_LANDING;} // check if close enough to dest, then switch to landing
		break;
	case BIRD_STATE_LANDING:
		if (is_anim_cycle_complete(new_anim_time)) { // wait for animation to complete
			state    = BIRD_STATE_STANDING;
			velocity = zero_vector; // stopped
			set_takeoff_time(rgen);
		}
		break;
	case BIRD_STATE_STANDING:
		if (takeoff_time == 0.0) { // set initial takeoff time
			set_takeoff_time(rgen);
		}
		else if (tfticks > takeoff_time && is_anim_cycle_complete(new_anim_time)) {
			if (placer.choose_bird_dest(radius, loc_ix, dest, dest_dir)) {state = BIRD_STATE_TAKEOFF;}
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
		//assert(dest_valid());
		//vector3d const dest_dir((dest - pos).get_norm());
		// TODO: update velocity
	}
	if (velocity != zero_vector) { // update direction and apply movement
		point const prev_pos(pos);
		dir    = vector3d(velocity.x, velocity.y, 0.0).get_norm(); // always point in XY direction of velocity
		pos   += velocity*timestep;
		bcube += (pos - bcube.get_cube_center()); // translate bcube to match - keep it centered on pos
		tile_changed |= (get_tile_id_containing_point_no_xyoff(pos) != get_tile_id_containing_point_no_xyoff(prev_pos));
	}
}

// this is here because only birds are updated each frame
void city_obj_placer_t::next_frame() {
	if (!animate2 || birds.empty()) return;
	float const enable_birds_dist(0.5f*(X_SCENE_SIZE + Y_SCENE_SIZE)); // half the pedestrian AI distance
	point const camera_bs(get_camera_pos() - get_tiled_terrain_model_xlate());
	if (!all_objs_bcube.closest_dist_less_than(camera_bs, enable_birds_dist)) return; // too far from the player
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	bool tile_changed(0);
	for (city_bird_t &bird : birds) {bird.next_frame(timestep, tile_changed, *this, bird_rgen);}
	
	if (tile_changed) { // update bird_groups; is there a more efficient way than rebuilding bird_groups each frame?
		bird_groups.clear();
		for (unsigned i = 0; i < birds.size(); ++i) {bird_groups.insert_obj_ix(birds[i].bcube, i);}
		bird_groups.create_groups(birds, all_objs_bcube);
	}
}

bool city_obj_placer_t::choose_bird_dest(float radius, unsigned &loc_ix, point &dest_pos, vector3d &dest_dir) {
	dest_pos = all_zeros; // reset
	if (bird_locs.size() == birds.size()) return 0; // all locs in use
	assert(loc_ix < bird_locs.size());
	bird_place_t &loc(bird_locs[loc_ix]);
	assert(loc.in_use);
	float const xlate_dist(2.0*radius); // move away from the object to avoid intersecting it
	
	for (unsigned n = 0; n < 10; ++n) { // make up to 10 attempts to avoid long runtime
		unsigned const new_loc_ix(bird_rgen.rand() % bird_locs.size());
		if (new_loc_ix == loc_ix) continue; // same
		bird_place_t &new_loc(bird_locs[new_loc_ix]);
		if (new_loc.in_use) continue;
		vector3d const dir((new_loc.pos - loc.pos).get_norm());
		point const start_pos(loc.pos + xlate_dist*dir), end_pos(new_loc.pos - xlate_dist*dir);
		float t(0.0); // unused
		if (line_intersect(start_pos, end_pos, t)) continue;
		bool valid(1);
		// TODO: check for buildings, etc.
		if (!valid) continue;
		dest_pos = new_loc.pos;
		dest_dir = new_loc.orient;
		loc    .in_use = 0;
		new_loc.in_use = 1;
		return 1; // success
	} // for n
	return 0; // failed, try again next frame or next animation cycle
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

