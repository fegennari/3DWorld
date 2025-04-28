// 3D World - City Birds Implementation
// by Frank Gennari
// 09/21/23

#include "city_objects.h"
#include "openal_wrap.h"
//#include "profiler.h"

float const BIRD_ACCEL       = 0.00025;
float const BIRD_MAX_VEL     = 0.002;
float const BIRD_ZV_RISE     = 0.4; // Z vs. XY velocity/acceleration for ascent
float const BIRD_ZV_FALL     = 0.8; // Z vs. XY velocity/acceleration for descent
float const BIRD_MAX_ALT_RES = 0.8; // above destination; in multiples of road width; for residential cities
float const BIRD_MAX_ALT_COM = 1.8; // above destination; in multiples of road width; for commercial  cities
float const anim_time_scale(1.25/TICKS_PER_SECOND);

extern bool camera_in_building, player_in_walkway;
extern int frame_counter;
extern float fticks;
extern double tfticks;
extern vector3d wind;
extern city_params_t city_params;
extern object_model_loader_t building_obj_model_loader;


void bird_poop_manager_t::next_frame(float fticks_stable) {
	if (poops.empty()) return;
	assert(!city_bounds.is_all_zeros()); // must call init() first
	float const z_ground(city_bounds.z2());
	float const gravity  = 0.0002; // unsigned magnitude
	float const term_vel = 0.1;

	for (poop_t &p : poops) { // what about wind?
		assert(p.radius > 0.0);
		p.vel.z -= gravity*fticks_stable; // apply gravitational acceleration
		max_eq(p.vel.z, -term_vel);
		p.pos += fticks_stable*p.vel;

		if (p.pos.z < z_ground || !city_bounds.contains_pt_xy(p.pos)) { // outside city bounds
			if (p.pos.z < z_ground) { // leave a splat
				point const spos(p.pos.x, p.pos.y, (z_ground + 0.25*p.radius)); // slightly above the ground
				splats.emplace_back(spos, 4.0*p.radius, rgen.signed_rand_vector_spherical_xy_norm());
				if (splats.size() > 40) {remove_oldest_splat();} // limit the max number of splats to 40
				gen_sound_thread_safe(SOUND_SM_SPLAT, (spos + get_camera_coord_space_xlate()));
			}
			p.radius = 0.0; // will be removed below
		}
	} // for p
	poops.erase(std::remove_if(poops.begin(), poops.end(), [](poop_t const &p) {return (p.radius == 0.0);}), poops.end());
}
void bird_poop_manager_t::remove_oldest_splat() {
	if (splats.empty()) return; // error?
	float tmin(0.0);
	unsigned oldest(0);

	for (auto s = splats.begin(); s != splats.end(); ++s) {
		if (tmin == 0.0 || s->time < tmin) {tmin = s->time; oldest = (s - splats.begin());}
	}
	swap(splats[oldest], splats.back());
	splats.pop_back();
}
void bird_poop_manager_t::draw(shader_t &s, vector3d const &xlate) {
	if (poops.empty() && splats.empty()) return;
	bind_default_flat_normal_map();

	if (!poops.empty()) {
		s.set_cur_color(WHITE);
		select_texture(WHITE_TEX);
		begin_sphere_draw(0); // textured=0
		unsigned const ndiv = 16;

		for (poop_t const &p : poops) {
			if (!camera_pdu.sphere_visible_test((p.pos + xlate), p.radius)) continue; // no distance check since poop should be dropped above the player
			draw_sphere_vbo(p.pos, p.radius, ndiv, 0); // textured=0
		}
		end_sphere_draw();
	}
	cube_t splat_range;

	for (splat_t const &s : splats) {
		point const pos_cs(s.pos + xlate);
		if (!dist_less_than(camera_pdu.pos, pos_cs, 400.0*s.radius) || !camera_pdu.sphere_visible_test(pos_cs, s.radius)) continue;
		s.add_quad(splat_qbd);
		splat_range.assign_or_union_with_pt(s.pos);
	}
	if (!splat_qbd.empty()) {
		select_texture(get_texture_by_name("splatter.png"));
		s.add_uniform_float("min_alpha", 0.9); // background color is black, so doesn't blend properly
		try_bind_tile_smap_at_point((splat_range.get_cube_center() + xlate), s); // view dist is low, so assume all splats are in the same tile
		splat_qbd.draw_and_clear();
		s.add_uniform_float("min_alpha", DEF_CITY_MIN_ALPHA); // reset
	}
}


city_bird_t::city_bird_t(point const &pos_, float height_, vector3d const &init_dir, colorRGBA const &color_, unsigned loc_ix_, rand_gen_t &rgen) :
	city_bird_base_t(pos_, height_, init_dir, OBJ_MODEL_BIRD_ANIM), state(BIRD_STATE_STANDING), loc_ix(loc_ix_), height(height_), prev_frame_pos(pos), color(color_)
{
	anim_time = 1.0*TICKS_PER_SECOND*rgen.rand_float(); // 1s random variation so that birds aren't all in sync
}

void city_bird_t::draw(draw_state_t &dstate, city_draw_qbds_t &qbds, float dist_scale, bool shadow_only) const {
	if (!shadow_only && dist_scale > 0.0 && !dist_less_than(dstate.camera_bs, pos, dist_scale*dstate.draw_tile_dist)) return;
	cube_t draw_bcube(bcube);
	draw_bcube.translate_dim(2, -0.28*radius); // required for GLB bird model

	if (dstate.is_visible_and_unoccluded(draw_bcube, dist_scale)) {
		// animations: 0=flying, 1=gliding, 2=landing, 3=standing, 4=takeoff
		float const model_anim_time(anim_time_scale*anim_time/SKELETAL_ANIM_TIME_CONST); // divide by constant to cancel out multiply in draw_model()
		animation_state_t anim_state(1, ANIM_ID_SKELETAL, model_anim_time, get_model_anim_id()); // enabled=1
		building_obj_model_loader.draw_model(dstate.s, pos, draw_bcube, dir, color, dstate.xlate, OBJ_MODEL_BIRD_ANIM, shadow_only, 0, &anim_state);
	}
	if (0 && !shadow_only && dest_valid() && is_close_to_player()) { // debug drawing, even if bcube not visible
		post_draw(dstate, shadow_only); // clear animations
		vector<vert_color> pts;
		pts.emplace_back(pos,  BLUE);
		pts.emplace_back(dest, BLUE);
		select_texture(WHITE_TEX);
		draw_verts(pts, GL_LINES);
		dstate.s.set_cur_color(BLUE);
		draw_sphere_vbo(dest, radius, N_SPHERE_DIV, 0); // draw destination
	}
}
/*static*/ void city_bird_t::post_draw(draw_state_t &dstate, bool shadow_only) {
	animation_state_t anim_state(1); // enabled=1
	anim_state.clear_animation_id(dstate.s); // clear animations
	bind_default_flat_normal_map();
}

bool city_bird_t::is_close_to_player() const {
	return dist_less_than(pos, get_camera_building_space(), 1.0*city_params.road_width);
}
void city_bird_t::set_takeoff_time(rand_gen_t &rgen) {
	takeoff_time = tfticks + rgen.rand_uniform(10.0, 30.0)*TICKS_PER_SECOND; // wait 10-30s
}
void city_bird_t::adjust_new_dest_zval() {
	dest.z += (0.5 + BIRD_ZVAL_ADJ)*height; // place feet at dest, not bird center
}
bool city_bird_t::is_anim_cycle_complete(float new_anim_time) const {
	if (anim_time == 0.0) return 0; // anim_time was just reset
	return building_obj_model_loader.check_anim_wrapped(OBJ_MODEL_BIRD_ANIM, get_model_anim_id(), anim_time_scale*anim_time, anim_time_scale*new_anim_time);
}
bool city_bird_t::in_landing_dist() const {
	assert(dest_valid());
	assert(state == BIRD_STATE_FLYING || state == BIRD_STATE_GLIDING);
	// play landing animation when dest is reached
	float const frame_dist(p2p_dist(pos, prev_frame_pos)), dist_thresh(frame_dist + 0.05*radius); // include previous frame distance to avoid overshoot
	return dist_less_than(pos, dest, dist_thresh);
}
bool city_bird_t::check_for_mid_flight_coll(float dir_dp, city_obj_placer_t &placer, rand_gen_t &rgen) {
	if (state != BIRD_STATE_FLYING && state != BIRD_STATE_GLIDING)     return 0; // only needed when flying or gliding
	if (((loc_ix + frame_counter) & 15) != 0)                          return 0; // check for pending collisions every 16 frames
	// here we assume there are no overhangs, so if we can see the point in front of and below us it must be reachable;
	// the descent approach is a straight line and should also follow this post=>dest path, except for the landing step at the end, which is too late to change dest;
	// we may even be able to get away with checking velocity.z >= 0.0, but that likely isn't necessary because descent shouldn't happen until dir is aligned;
	// however, this is only legal when descending/gliding due to other dir changes that can be made while flying (in particular when ascending)
	if (state == BIRD_STATE_GLIDING && dir_dp > 0.99)                  return 0;
	if (!placer.check_path_segment_coll(pos, dest, radius))            return 0;
	if (!placer.choose_bird_dest(pos, radius, loc_ix, dest, dest_dir)) return 0; // if this fails, continue to original dest
	adjust_new_dest_zval();
	max_eq(start_end_zmax, dest.z);
	// no state update since we're already flying
	return 1;
}

// timestep is in ticks
void city_bird_t::next_frame(float timestep, float delta_dir, point const &camera_bs, bool &tile_changed, bool &bird_moved, city_obj_placer_t &placer, rand_gen_t &rgen) {
	// update state
	uint8_t const prev_state(state);
	float const new_anim_time(anim_time + timestep);

	switch (state) {
		// Note: if flying level (velocity.z == 0.0), maintain the current flying vs. gliding state
	case BIRD_STATE_FLYING:
		if (velocity.z < 0.0 && is_anim_cycle_complete(new_anim_time)) {state = BIRD_STATE_GLIDING;} // maybe switch to gliding
		// fall through
	case BIRD_STATE_GLIDING:
		if (velocity.z > 0.0) {state = BIRD_STATE_FLYING;} // should we call is_anim_cycle_complete()?
		
		if (in_landing_dist()) { // check if close enough to dest, then switch to landing
			state = BIRD_STATE_LANDING;
			pos   = dest; // snap to dest; is this a good idea? move more slowly?
			dest  = all_zeros; // clear for next cycle
		}
		break;
	case BIRD_STATE_LANDING:
		if (is_anim_cycle_complete(new_anim_time)) { // wait for animation to complete
			state       = BIRD_STATE_STANDING;
			velocity    = zero_vector; // make sure to stop
			hit_min_alt = 0; // reset for next cycle
			set_takeoff_time(rgen);
		}
		break;
	case BIRD_STATE_STANDING:
		if (takeoff_time == 0.0) { // set initial takeoff time
			set_takeoff_time(rgen);
		}
		else if (tfticks > takeoff_time && is_anim_cycle_complete(new_anim_time)) {
			if (placer.choose_bird_dest(pos, radius, loc_ix, dest, dest_dir)) {
				adjust_new_dest_zval();
				start_end_zmax = max(dest.z, pos.z);
				state          = BIRD_STATE_TAKEOFF;
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

	if (state == BIRD_STATE_STANDING || state == BIRD_STATE_LANDING || pos == dest) { // stationary
		if (dest_dir != zero_vector) { // maybe update direction to match intended direction of our new destination
			float const dp(dot_product(dir, dest_dir));

			if (dp < 0.999) { // not oriented in dir
				if (dp < 0.0) {dest_dir.negate();} // dest_dir is typically something like a wall, so facing the opposite direction should be okay
				dir = delta_dir*dest_dir + (1.0 - delta_dir)*dir; // merge new_dir into dir gradually for smooth turning
				dir.normalize();
			}
		}
		velocity = zero_vector; // stopped
	}
	else { // update direction and velocity
		vector3d dest_dir_xy(dest.x-pos.x, dest.y-pos.y, 0.0); // XY only
		float const dist_xy(dest_dir_xy.mag());
		if (dist_xy > TOLERANCE) {dest_dir_xy /= dist_xy;} else {dest_dir_xy = dir;} // normalize if nonzero, otherwise use cur dir
		float const accel(((state == BIRD_STATE_TAKEOFF) ? 0.1 : 1.0)*BIRD_ACCEL); // slower during takeoff
		float const dir_dp(dot_product(dir, dest_dir_xy));

		// update direction
		if (dir_dp < 0.999) { // not oriented in dir
			delta_dir *= max(1.0f, 4.0f*radius/dist_xy); // faster correction when close to destination to prevent circling
			//if (dir_dp < 0.9 && state == BIRD_STATE_FLYING && ((loc_ix + frame_counter + 1) & 7) == 0) {} // check for pending collision?
			// if directions are nearly opposite, pick a side to turn using the cross product to get an orthogonal vector
			if (dir_dp < -0.9) {dest_dir_xy = cross_product(dir, plus_z)*((loc_ix & 1) ? -1.0 : 1.0);} // random turn direction (CW/CCW)
			dir = delta_dir*dest_dir_xy + (1.0 - delta_dir)*dir; // merge new_dir into dir gradually for smooth turning
			dir.normalize();
		}
		// handle vertical velocity component; this assumes the start and end points are similar zvals
		// check dir is aligned and approach angle is steep enough to glide down;
		// the extra radius factor allows for smooth landings when close; closer to radius is more smooth, while omitting this results in flappy landings
		bool const can_descend(hit_min_alt && dir_dp > 0.5 && BIRD_ZV_FALL*dist_xy < (pos.z - dest.z - 0.1*radius));

		if (!can_descend) { // ascent stage
			float const max_alt_mult(placer.has_residential() ? BIRD_MAX_ALT_RES : BIRD_MAX_ALT_COM);
			float const max_alt(max_alt_mult*(city_params.road_width + 4.0*radius)); // larger birds can fly a bit higher
			hit_min_alt |= (pos.z > start_end_zmax + 4.0*radius); // lift off at least 4x radius to clear the starting object
			if (pos.z > max_alt + start_end_zmax) {velocity.z = 0.0;} // too high; level off (should we decelerate smoothly?)
			else {velocity.z = min(BIRD_ZV_RISE*BIRD_MAX_VEL, (velocity.z + BIRD_ZV_RISE*accel));} // rise in the air
		}
		else { // descent stage
			float const z_clearance(min(1.0f*radius, 0.5f*dist_xy)); // approach angle not too shallow
			
			if (pos.z <= dest.z + z_clearance) { // don't fall below the destination zval
				velocity.z = 0.0;
				max_eq(pos.z, dest.z); // can't go below dest
			}
			else {velocity.z = max(-BIRD_ZV_FALL*BIRD_MAX_VEL, (velocity.z - BIRD_ZV_FALL*accel));} // glide down
		}
		// handle horizontal velocity component
		float const frame_dist_xy(p2p_dist_xy(pos, prev_frame_pos)), dist_thresh(frame_dist_xy + 0.03*radius);
		
		if (dist_xy_less_than(pos, dest, dist_thresh)) { // close enough in XY, but may need to descend in Z
			velocity.x = velocity.y = 0.0;
		}
		else {
			velocity += accel*dir; // accelerate in XY; acceleration always follows direction
			float const max_vel(BIRD_MAX_VEL*max(0.2f, min(1.0f, 0.05f*dist_xy/radius))); // prevent overshoot and circling by slowing when close
			float const xy_mag(velocity.xy_mag());
		
			if (xy_mag > max_vel) { // apply velocity cap
				float const xy_scale(max_vel/xy_mag);
				velocity.x *= xy_scale;
				velocity.y *= xy_scale;
			}
			check_for_mid_flight_coll(dir_dp, placer, rgen);
		}
		if (dist_xy < radius) { // special case to pull bird in when close
			vector3d const delta(dest - pos);
			float const dist(delta.mag()), vmag(max(min(dist/timestep, min(velocity.mag(), 0.2f*BIRD_MAX_VEL)), 0.01f*BIRD_MAX_VEL)); // limit vmag to not overshoot
			velocity = delta*(vmag/dist);
		}
	}
	if (velocity != zero_vector) { // apply movement
		prev_frame_pos = pos;
		pos           += velocity*timestep;
		bcube         += (pos - bcube.get_cube_center()); // translate bcube to match - keep it centered on pos
		tile_changed  |= (get_tile_id_containing_point_no_xyoff(pos) != get_tile_id_containing_point_no_xyoff(prev_frame_pos));
		bird_moved     = 1;
		
		// poop on the player if above and the player is in the open
		if (!camera_in_building && !player_in_walkway && pos.z > camera_bs.z && camera_bs.z > placer.city_zval &&
			dist_xy_less_than(pos, camera_bs, 2.0*CAMERA_RADIUS) && tfticks > next_poop_time && !placer.player_under_roof(camera_bs))
		{
			placer.add_bird_poop(pos, 0.18*radius, (velocity + 0.001*wind)); // use bird's initial velocity and add a small amount of wind; should wind apply acceleration?
			next_poop_time = tfticks + rgen.rand_uniform(2.0, 5.0)*TICKS_PER_SECOND; // wait 2-5s before pooping again
		}
	}
}

// this template function is here rather than in city_obj_placer.cpp because it's currently only used for birds
template<typename T> void city_obj_groups_t::update_obj_pos(vector<T> const &objs, cube_t &all_objs_bcube) {
	unsigned start_ix(0);

	for (auto &g : *this) {
		cube_t const prev(g);
		g.set_to_zeros(); // clear bcube for this pass
		assert(start_ix <= g.ix && g.ix <= objs.size());
		for (auto i = objs.begin()+start_ix; i != objs.begin()+g.ix; ++i) {g.assign_or_union_with_cube(i->bcube);}
		bcube.union_with_cube(g);
		start_ix = g.ix;
	} // for g
	all_objs_bcube.union_with_cube(bcube);
}

void city_obj_placer_t::next_frame_birds(point const &camera_bs, float fticks_stable) {
	if (birds.empty()) return;
	float const enable_birds_dist(0.5f*(X_SCENE_SIZE + Y_SCENE_SIZE)); // half the pedestrian AI distance
	if (!all_objs_bcube.closest_dist_less_than(camera_bs, enable_birds_dist)) return; // too far from the player
	//highres_timer_t timer("Update Birds");
	float const timestep(min(fticks, 4.0f)); // clamp fticks to 100ms
	float const delta_dir(0.1*(1.0 - pow(0.7f, fticks))); // controls bird turning rate
	bool tile_changed(0), bird_moved(0);
	for (city_bird_t &bird : birds) {bird.next_frame(timestep, delta_dir, camera_bs, tile_changed, bird_moved, *this, bird_rgen);}
	
	if (tile_changed) { // update bird_groups; is there a more efficient way than rebuilding bird_groups each frame?
		bird_groups.rebuild(birds, all_objs_bcube);
	}
	else if (bird_moved) { // incrementally update group bcubes and all_objs_bcube
		bird_groups.update_obj_pos(birds, all_objs_bcube);
	}
	bird_poop_manager.next_frame(fticks_stable);
}

bool city_obj_placer_t::choose_bird_dest(point const &pos, float radius, unsigned &loc_ix, point &dest_pos, vector3d &dest_dir) {
	if (bird_locs.size() == birds.size()) return 0; // all locs in use
	assert(loc_ix < bird_locs.size());
	bird_place_t &old_loc(bird_locs[loc_ix]);
	assert(old_loc.in_use);
	float const xlate_dist(2.0*radius); // move away from the object to avoid intersecting it
	bool const find_closest = 0; // makes for easier debugging
	vector<pair<float, unsigned>> pref_locs;

	if (find_closest) {
		for (unsigned i = 0; i < bird_locs.size(); ++i) {
			if (i == loc_ix || bird_locs[i].in_use) continue;
			pref_locs.emplace_back(p2p_dist_xy(pos, bird_locs[i].pos), i);
		}
		sort(pref_locs.begin(), pref_locs.end()); // sort closes to furthest
	}
	// else find a random visible/reachable dest
	for (unsigned n = 0; n < 10; ++n) { // make up to 10 attempts to avoid long runtime
		unsigned const new_loc_ix((n < pref_locs.size()) ? pref_locs[n].second : (bird_rgen.rand() % bird_locs.size()));
		if (new_loc_ix == loc_ix) continue; // same
		bird_place_t &new_loc(bird_locs[new_loc_ix]);
		if (new_loc.in_use) continue;
		assert(new_loc.pos != pos);
		if (new_loc.pos.z > pos.z + BIRD_ZV_RISE*p2p_dist_xy(new_loc.pos, pos)) continue; // too steep a rise
		if (new_loc.pos.z < pos.z - BIRD_ZV_FALL*p2p_dist_xy(new_loc.pos, pos)) continue; // too steep a drop
		vector3d const dir((new_loc.pos - pos).get_norm());
		point const start_pos(pos + xlate_dist*dir), end_pos(new_loc.pos - xlate_dist*dir);
		if (check_path_segment_coll(start_pos, end_pos, radius)) continue;
		loc_ix   = new_loc_ix;
		dest_pos = new_loc.pos;
		if (new_loc.use_orient) {dest_dir = new_loc.orient;} // prefer this orient, for example if standing on a wall
		old_loc.in_use = 0;
		new_loc.in_use = 1;
		return 1; // success
	} // for n
	return 0; // failed, try again next frame or next animation cycle
}

// returns: 0=no coll, 1=city object coll, 2=building coll
int city_obj_placer_t::check_path_segment_coll(point const &p1, point const &p2, float radius) const {
	float t(0.0); // unused
	if (line_intersect(p1, p2, t)) return 1;
	if (check_city_building_line_coll_bs_any(p1, p2)) return 2;

	if (radius > 0.0) { // projected cylinder case: check 4 points a distance radius from the center
		vector3d const dir((p2 - p1).get_norm()), v_side(cross_product(dir, plus_z).get_norm()); // orthogonalize_dir?
		vector3d const z_off(0.0, 0.0, 2.0*radius); // move upward to clear any low-lying obstacles such as the source and dest objects, since we'll be flying upward anyway
		vector3d const offs[3] = {v_side, -v_side, -plus_z}; // check both sides and the point below at the feet

		for (unsigned n = 0; n < 3; ++n) {
			vector3d const off(radius*offs[n] + z_off);
			point const p1o(p1 + off), p2o(p2 + off);
			if (line_intersect(p1o, p2o, t)) return 1; // doesn't include objects such as power lines

			for (walkway_t const &w : walkways) { // special handling of walkways so that we don't try to land on power lines just under them
				cube_t bc_ext(w.bcube);
				bc_ext.z1() -= w.floor_spacing; // extend down by one floor
				if (bc_ext.line_intersects(p1o, p2o)) return 1;
			}
			if (skyway.valid) { // check skyway intersection
				cube_t bc_ext(skyway.bcube);
				bc_ext.z1() -= 0.5*skyway.bcube.dz(); // extend down by half height
				if (bc_ext.line_intersects(p1o, p2o)) return 1;
			}
			if (check_city_building_line_coll_bs_any(p1o, p2o)) return 2; // doesn't include objects such as building rooftop signs?
		} // for n
	}
	return 0;
}

bool city_obj_placer_t::check_bird_walkway_clearance(cube_t const &bc) const { // and skyway, and elevators
	for (walkway_t const &w : walkways) {
		cube_t bc_ext(w.bcube);
		bc_ext.z1() -= w.floor_spacing; // extend lower edge by one floor for clearance
		if (bc_ext.intersects_xy(bc)) return 0;
	}
	for (ww_elevator_t const &e : elevators) {
		if (e.bcube.intersects_xy(bc)) return 0;
	}
	if (skyway.valid) { // check skyway intersection
		cube_t bc_ext(skyway.bcube);
		bc_ext.z1() -= 0.5*skyway.bcube.dz(); // extend down by half height
		if (bc_ext.intersects_xy(bc)) return 0;
	}
	return 1;
}

template<typename T> bool check_obj_under(vector<T> const &objs, point const &pos) {
	for (T const &i : objs) {
		if (pos.z < i.bcube.z2() && i.bcube.contains_pt_xy(pos)) return 1;
	}
	return 0;
}
bool city_obj_placer_t::player_under_roof(point const &camera_bs) const { // called for bird poop logic
	return (check_obj_under(p_solars, camera_bs) || check_obj_under(walkways, camera_bs)); // what about pool deck roofs?
}

void city_obj_placer_t::add_bird_poop(point const &pos, float radius, vector3d const &init_vel) {
	if (skyway.valid && pos.z > skyway.bcube.z1() && skyway.bcube.contains_pt_xy(pos)) return; // don't poop over skyway

	for (walkway_t const &w : walkways) {
		if (pos.z > w.bcube.z1() && w.bcube.contains_pt_xy(pos)) return; // don't poop over walkways
	}
	bird_poop_manager.add(pos, radius, init_vel);
}

void vect_bird_place_t::add_placement(cube_t const &obj, bool dim, bool dir, bool orient_dir, float spacing, rand_gen_t &rgen) {
	point pos(0.0, 0.0, obj.z2());
	pos[ dim] = obj.d[dim][dir];
	pos[!dim] = rgen.rand_uniform(obj.d[!dim][0]+spacing, obj.d[!dim][1]-spacing); // random position along the edge
	emplace_back(pos, dim, orient_dir);
}
void vect_bird_place_t::add_placement_centerline(cube_t const &obj, bool dim, bool dir, rand_gen_t &rgen) {
	float const spacing(0.1*obj.get_sz_dim(!dim)); // 10% from either end
	point pos(0.0, 0.0, obj.z2());
	pos[ dim] = obj.get_center_dim(dim);
	pos[!dim] = rgen.rand_uniform(obj.d[!dim][0]+spacing, obj.d[!dim][1]-spacing); // random position along the centerline
	emplace_back(pos, dim, dir);
}
void vect_bird_place_t::add_placement_rand_dim_dir(cube_t const &obj, float spacing, rand_gen_t &rgen) {
	bool const dir(rgen.rand_bool());
	add_placement(obj, rgen.rand_bool(), dir, dir, spacing, rgen);
}
void vect_bird_place_t::add_placement_top_center(cube_t const &obj, rand_gen_t &rgen) {
	emplace_back(cube_top_center(obj), rgen.rand_bool(), rgen.rand_bool(), 0); // use_orient=0
}

