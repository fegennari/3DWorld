// 3D World - City Birds Implementation
// by Frank Gennari
// 09/21/23

#include "city_objects.h"

float const BIRD_ACCEL    = 0.0005;
float const BIRD_MAX_VEL  = 0.002;
float const BIRD_ZV_SCALE = 0.5; // Z vs. XY velocity/acceleration
float const BIRD_MAX_ALT  = 2.0; // above destination; in multiples of road width
float const anim_time_scale(1.0/TICKS_PER_SECOND);

extern int animate2;
extern float fticks;
extern double tfticks;
extern city_params_t city_params;
extern object_model_loader_t building_obj_model_loader;

enum {BIRD_STATE_FLYING=0, BIRD_STATE_GLIDING, BIRD_STATE_LANDING, BIRD_STATE_STANDING, BIRD_STATE_TAKEOFF, NUM_BIRD_STATES};


city_bird_t::city_bird_t(point const &pos_, float height, vector3d const &init_dir, unsigned loc_ix_, rand_gen_t &rgen) :
	city_bird_base_t(pos_, height, init_dir, OBJ_MODEL_BIRD_ANIM), state(BIRD_STATE_STANDING), loc_ix(loc_ix_)
{
	anim_time = 1.0*TICKS_PER_SECOND*rgen.rand_float(); // 1s random variation so that birds aren't all in sync
}

void city_bird_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (dstate.check_cube_visible(bcube, dist_scale)) {
		// animations: 0=flying, 1=gliding, 2=landing, 3=standing, 4=takeoff
		float const model_anim_time(anim_time_scale*anim_time/SKELETAL_ANIM_TIME_CONST); // divide by constant to cancel out multiply in draw_model()
		animation_state_t anim_state(1, ANIM_ID_SKELETAL, model_anim_time, get_model_anim_id()); // enabled=1
		building_obj_model_loader.draw_model(dstate.s, pos, bcube, dir, WHITE, dstate.xlate, OBJ_MODEL_BIRD_ANIM, shadow_only, 0, &anim_state);
	}
	if (0 && dest_valid()) { // debug drawing, even if bcube not visible
		post_draw(dstate, shadow_only); // clear animations
		vector<vert_color> pts;
		pts.emplace_back(pos,  BLUE);
		pts.emplace_back(dest, BLUE);
		select_texture(WHITE_TEX);
		draw_verts(pts, GL_LINES);
	}
}
/*static*/ void city_bird_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	animation_state_t anim_state(1); // enabled=1
	anim_state.clear_animation_id(dstate.s); // clear animations
	model3d::bind_default_flat_normal_map();
}

bool city_bird_t::is_close_to_player() const {
	return dist_less_than(pos, (get_camera_pos() - get_tiled_terrain_model_xlate()), 1.0*city_params.road_width);
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
	if (velocity.x == 0.0 && velocity.y == 0.0) return 0; // not moving in XY
	float const land_delay_secs(building_obj_model_loader.get_anim_duration(OBJ_MODEL_BIRD_ANIM, get_model_anim_id()));
	float const land_delay_ticks(land_delay_secs/anim_time_scale);
	vector3d const delta_xy(dest.x-pos.x, dest.y-pos.y, 0.0), v_xy(velocity.x, velocity.y, 0.0); // XY only
	float const dest_time(delta_xy.mag_sq()/dot_product(velocity, delta_xy)); // distance/velocity_to_dest = delta.mag()/dot_product(velocity, delta/deta.mag())
	if (dest_time < 0.0) return 0; // flying away from dest

	if (is_close_to_player()) {
		float const dest_time_secs(dest_time/TICKS_PER_SECOND);
		float const dp(dot_product(v_xy, delta_xy)/(v_xy.mag()*delta_xy.mag()));
		cout << TXT(land_delay_secs) << TXT(dest_time_secs) << endl;
	}
	return (dest_time < land_delay_ticks);
}
float city_bird_t::get_path_progress() const {
	assert(dest_valid());
	if (init_dest_dist == 0.0) return 0.0; // error?
	return CLIP_TO_01(1.0f - p2p_dist_xy(pos, dest)/init_dest_dist); // 0.0 at start, 1.0 at end
}

// timestep is in ticks
void city_bird_t::next_frame(float timestep, float delta_dir, bool &tile_changed, bool &bird_moved, city_obj_placer_t &placer, rand_gen_t &rgen) {
	// update state
	uint8_t const prev_state(state);
	float const new_anim_time(anim_time + timestep);

	switch (state) {
		// Note: if flying level (velocity.z == 0.0), maintain the current flying vs. gliding state
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
			if (placer.choose_bird_dest(radius, loc_ix, dest, dest_dir)) {
				init_dest_dist = p2p_dist_xy(pos, dest);
				state = BIRD_STATE_TAKEOFF;
			}
		}
		break;
	case BIRD_STATE_TAKEOFF:
		if (is_anim_cycle_complete(new_anim_time)) {state = BIRD_STATE_FLYING;} // wait for animation to complete
		break;
	default: assert(0);
	}
	if (state != prev_state) {anim_time = 0.0;} // reset animation time on state change
	else {anim_time = new_anim_time;}

	if (state != BIRD_STATE_STANDING) { // update direction and velocity
		vector3d dest_dir_xy(vector3d(dest.x-pos.x, dest.y-pos.y, 0.0).get_norm()); // XY only

		// update direction
		if (dot_product(dir, dest_dir_xy) < 0.999) { // not oriented in dir
			// if directions are nearly opposite, pick a side to turn using the cross product to get an orthogonal vector
			if (dot_product(dest_dir_xy, dir) < -0.9) {dest_dir_xy = cross_product(dir, plus_z)*((loc_ix & 1) ? -1.0 : 1.0);} // random turn direction (CW/CCW)
			dir = delta_dir*dest_dir_xy + (1.0 - delta_dir)*dir; // merge new_dir into dir gradually for smooth turning
			dir.normalize();
		}
		// handle vertical velocity component
		float const path_progress(get_path_progress());

		if (path_progress < 0.5) { // first half of path
			if (pos.z > BIRD_MAX_ALT*city_params.road_width + dest.z) {velocity.z = 0.0;} // too high; level off (should we decelerate smoothly?)
			else {velocity.z = min( BIRD_ZV_SCALE*BIRD_MAX_VEL, (velocity.z + BIRD_ZV_SCALE*BIRD_ACCEL));} // rise in the air
		}
		else { // second half of path
			if (pos.z <= dest.z) {velocity.z = 0.0;} // don't fall below the destination zval
			else {velocity.z = max(-BIRD_ZV_SCALE*BIRD_MAX_VEL, (velocity.z - BIRD_ZV_SCALE*BIRD_ACCEL));} // glide down
		}
		// handle horizontal velocity component
		float const accel_scale((state == BIRD_STATE_TAKEOFF) ? 0.5 : 1.0); // slower during takeoff
		velocity += BIRD_ACCEL*accel_scale*dir; // accelerate in XY; acceleration always follows direction
		float const xy_mag(velocity.xy_mag());
		
		if (xy_mag > BIRD_MAX_VEL) { // apply velocity cap
			float const xy_scale(BIRD_MAX_VEL/xy_mag);
			velocity.x *= xy_scale;
			velocity.y *= xy_scale;
		}
	}
	if (velocity != zero_vector) { // apply movement
		point const prev_pos(pos);
		pos   += velocity*timestep;
		bcube += (pos - bcube.get_cube_center()); // translate bcube to match - keep it centered on pos
		tile_changed |= (get_tile_id_containing_point_no_xyoff(pos) != get_tile_id_containing_point_no_xyoff(prev_pos));
		bird_moved    = 1;
	}
}

// this is here because only birds are updated each frame
void city_obj_placer_t::next_frame() {
	if (!animate2 || birds.empty()) return;
	float const enable_birds_dist(0.5f*(X_SCENE_SIZE + Y_SCENE_SIZE)); // half the pedestrian AI distance
	point const camera_bs(get_camera_pos() - get_tiled_terrain_model_xlate());
	if (!all_objs_bcube.closest_dist_less_than(camera_bs, enable_birds_dist)) return; // too far from the player
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	float const delta_dir(0.1*(1.0 - pow(0.7f, fticks))); // controls bird turning rate
	bool tile_changed(0), bird_moved(0);
	for (city_bird_t &bird : birds) {bird.next_frame(timestep, delta_dir, tile_changed, bird_moved, *this, bird_rgen);}
	
	if (tile_changed) { // update bird_groups; is there a more efficient way than rebuilding bird_groups each frame?
		//cout << TXT(birds.size()) << TXT(bird_groups.size()) << endl; // TESTING
		bird_groups.clear();
		for (unsigned i = 0; i < birds.size(); ++i) {bird_groups.insert_obj_ix(birds[i].bcube, i);}
		bird_groups.create_groups(birds, all_objs_bcube);
	}
	else if (bird_moved) { // incrementally update group bcubes and all_objs_bcube
		unsigned start_ix(0);

		for (auto &g : bird_groups) {
			cube_t const prev(g);
			g.set_to_zeros(); // clear bcube for this pass
			assert(start_ix <= g.ix && g.ix <= birds.size());
			for (auto i = birds.begin()+start_ix; i != birds.begin()+g.ix; ++i) {g.assign_or_union_with_cube(i->bcube);}
			all_objs_bcube.union_with_cube(g);
			start_ix = g.ix;
		} // for g
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
		if (check_path_segment_coll(start_pos, end_pos, radius)) continue;
		loc_ix   = new_loc_ix;
		dest_pos = new_loc.pos;
		dest_dir = new_loc.orient;
		loc    .in_use = 0;
		new_loc.in_use = 1;
		return 1; // success
	} // for n
	return 0; // failed, try again next frame or next animation cycle
}
bool city_obj_placer_t::check_path_segment_coll(point const &p1, point const &p2, float radius) const {
	float t(0.0); // unused
	if (line_intersect(p1, p2, t) || check_city_building_line_coll_bs_any(p1, p2)) return 1;

	if (radius > 0.0) { // cylinder case: check 4 points a distance radius from the center
		vector3d const dir((p2 - p1).get_norm()), v1(cross_product(dir, plus_z).get_norm()), v2(cross_product(v1, dir).get_norm()); // orthogonalize_dir?
		vector3d const offs[4] = {v1, -v1, v2, -v2};

		for (unsigned n = 0; n < 4; ++n) {
			vector3d const off(radius*offs[n]);
			point const p1o(p1 + off), p2o(p2 + off);
			if (line_intersect(p1o, p2o, t) || check_city_building_line_coll_bs_any(p1o, p2o)) return 1;
		}
	}
	return 0;
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

