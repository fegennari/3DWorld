// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#include "animals.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "model3d.h"
#include "tiled_mesh.h"
#include "format_text.h"

float const FISH_RADIUS = 0.05;
float const BIRD_RADIUS = 0.1;
float const BFLY_RADIUS = 0.015;
float const FISH_SPEED  = 0.002;
float const BIRD_SPEED  = 0.05;
float const BFLY_SPEED  = 0.0045;

extern bool water_is_lava;
extern int window_width, animate2, display_mode;
extern float fticks, water_plane_z, temperature, atmosphere, ocean_wave_height, model_mat_lod_thresh, zmax_est, tree_scale;
extern colorRGBA cur_fog_color;


void do_xy_rotate_normal(float rot_sin, float rot_cos, point &pos);
bool choose_pt_in_city_park(point const &pos, point &park_pos, rand_gen_t &rgen);

bool birds_active     () {return (light_factor >= 0.4);} // birds are only active whe the sun is out
bool debug_animal_draw() {return (0/*display_mode & 0x10*/);}

enum {BF_BODY=0, BF_LWING, BF_RWING, BF_NUM_PARTS};

class animal_model_loader_t : public model3ds { // currently for fish only

	struct model_info_t {
		bool invalid;
		unsigned id;
		model_info_t() : invalid(0), id(0) {}
		bool is_loaded() const {return (id > 0);}
		bool try_load () const {return (!is_loaded() && !invalid);}
		void set_id(unsigned id_) {id = id_; invalid = (id == 0);}
	};
	model_info_t fish_info, bfly_info[BF_NUM_PARTS];

	model3d &get_model(unsigned id) {
		assert(id > 0);
		assert(id-1 < size()); // Note: id is vector index offset by 1
		return operator[](id-1);
	}
	unsigned load_model(string const &fn, colorRGBA const &def_color=WHITE, int def_tid=-1) {
		unsigned const id(size());
		bool const write_file    = 0;
		int const recalc_normals = 0; // okay for loading model3d

		if (!load_model_file(fn, *this, geom_xform_t(), "", def_tid, def_color, 0, 0.0, 1.0, recalc_normals, 0, write_file, 1)) { // no animations in this model
			cerr << format_red("Error: Failed to read model file '" + fn + "'") << endl;
			return 0;
		}
		return id + 1;
	}
	void draw_model(unsigned id, shader_t &s, vector3d const &pos, float radius, vector3d const &dir, rotation_t const &local_rotate,
		bool cancel_local_translate, colorRGBA const &color=WHITE, bool use_custom_tid=0, bool is_shadow_pass=0, float lod_mult=1.0)
	{
		if (color != WHITE) {s.add_uniform_color("color_modulate", color);}
		model3d &model(get_model(id));
		model.bind_all_used_tids();
		cube_t const &bcube(model.get_bcube());
		float const bcube_radius(0.5*bcube.max_len()), sz_scale(radius / bcube_radius);
		float const lod_dist(p2p_dist(get_camera_pos(), pos));
		bool const camera_pdu_valid(camera_pdu.valid);
		camera_pdu.valid = 0; // disable VFC, since we're doing custom transforms here
		lod_mult *= sz_scale/model_mat_lod_thresh; // model_mat_lod_thresh doesn't apply here, so divide to cancel it out
		int orig_tid(-1);
		unsigned const custom_tex_mat_id(0); // currently hard-coded to 0
		if (use_custom_tid) {orig_tid = model.set_texture_for_material(custom_tex_mat_id, -2);} // use a special tid of -2 to skip model texture bind
		fgPushMatrix();
		translate_to(pos);
		rotate_to_plus_x(dir);
		local_rotate.apply_gl();
		uniform_scale(sz_scale);
		if (cancel_local_translate) {translate_to(-bcube.get_cube_center());}
		model.render_materials(s, is_shadow_pass, 0, 0, 1, 3, 3, model.get_unbound_material(), rotation_t(), nullptr, nullptr, 0, lod_mult, lod_dist, 0, 1); // scaled
		if (use_custom_tid) {model.set_texture_for_material(custom_tex_mat_id, orig_tid);} // restore original value
		if (color != WHITE) {s.add_uniform_color("color_modulate", WHITE);} // reset
		fgPopMatrix();
		camera_pdu.valid = camera_pdu_valid;
	}
public:
	bool load_fish_model() {
		if (fish_info.try_load()) {fish_info.set_id(load_model("model_data/fish/fishOBJ.model3d"));}
		return fish_info.is_loaded();
	}
	bool load_butterfly_model() {
		string const part_fns[BF_NUM_PARTS] = {"Butterfly Body", "Butterfly Left Wing", "Butterfly Right Wing"};
		cube_t all_bcube;

		for (unsigned n = 0; n < BF_NUM_PARTS; ++n) {
			model_info_t &minfo(bfly_info[n]);
			if (minfo.try_load()) {minfo.set_id(load_model("../models/butterfly/" + part_fns[n] + ".model3d"));}
			if (!minfo.is_loaded()) return 0; // failed to load this part
			all_bcube.assign_or_union_with_cube(back().get_bcube());
		}
		for (unsigned n = 0; n < BF_NUM_PARTS; ++n) {get_model(bfly_info[n].id).union_bcube_with(all_bcube);} // the entire model has one unified bcube
		return 1; // success
	}
	void draw_fish_model(shader_t &s, vector3d const &pos, float radius, vector3d const &dir, float anim_time, colorRGBA const &color=WHITE) {
		rotation_t const local_rotate(plus_y, -45.0); // fish model is angled 45 degree upward, so need to rotate it back down
		// invert dir - fish model faces in -x, and we want it to be in +x
		s.add_uniform_float("animation_time", anim_time);
		draw_model(fish_info.id, s, pos, radius, -dir, local_rotate, 1, color); // not shadow pass, cancel local translate
	}
	void draw_butterfly_model(shader_t &s, vector3d const &pos, float radius, vector3d const &dir, float rot_time, bool draw_body, colorRGBA const &color_in=WHITE) {
		float const flap_angle(5.0*TWO_PI*rot_time), rot_angle(45.0f*(sin(flap_angle) + 0.5)); // 5 flaps per second; more positive angle
		point pos_mod(pos);
		pos_mod.z += 0.25*radius*sin(flap_angle - 0.25f*PI); // add some vertical bounce as the wings move; 22.5 degree lag
		colorRGBA const color(color_in, 1.0); // make alpha 1.0
		bool const use_custom_tid(color_in != WHITE);

		if (use_custom_tid) { // not a Monarch butterfly
			select_texture(get_texture_by_name((color == ALPHA0) ? "../models/butterfly/TextureAtlasWhite.png" : "../models/butterfly/TextureAtlasRecolor.png"));
		}
		for (unsigned n = BF_LWING; n <= BF_RWING; ++n) { // draw wings
			float const wing_sign((n == BF_LWING) ? 1.0 : -1.0);
			rotation_t const wing_rotate(plus_x, (90.0f + wing_sign*rot_angle));
			// this isn't perfect because in reality there are two wings on each side, and each of the four wings rotates about a slightly different point;
			// the result has a bit of clipping of wings through the body and each other, but maybe it's close enough
			draw_model(bfly_info[n].id, s, pos_mod, radius, dir, wing_rotate, 0, color, use_custom_tid);
		}
		if (draw_body) {
			rotation_t const body_rotate(plus_x, 90.0); // model has up=y, and we want up=z
			draw_model(bfly_info[BF_BODY].id, s, pos_mod, radius, dir, body_rotate, 0, WHITE, 0, 0, 1.5); // not shadow pass, custom lod_mult for legs
		}
	}
};

animal_model_loader_t animal_model_loader;

void free_animal_context() {animal_model_loader.free_context();}

// static functions for animal types
bool fish_t     ::type_enabled() {return animal_model_loader.load_fish_model();}
bool butterfly_t::type_enabled() {return animal_model_loader.load_butterfly_model();}

bool fish_t     ::can_place_in_tile(tile_t const *const tile) {return tile->has_water();}
bool butterfly_t::can_place_in_tile(tile_t const *const tile) {return (!tile->all_water() && tile->has_grass());} // only spawn over grass

bool have_fish_model() {return fish_t::type_enabled();}


void animal_t::gen_dir_vel(rand_gen_t &rgen, float speed) {
	dir      = rgen.signed_rand_vector_xy().get_norm(); // in xy plane
	velocity = (speed*rgen.rand_uniform(0.6, 1.0))*dir;
}

float animal_t::get_mesh_zval_at_pos(tile_t const *const tile) const {
	// use tile to interpolate z-value if present; faster, but should return a similar value to interpolate_mesh_zval()
	if (tile) {return tile->get_zval_at(pos.x, pos.y, 1);} // in global space
	return interpolate_mesh_zval((pos.x - DX_VAL*xoff2), (pos.y - DY_VAL*yoff2), 0.0, 1, 1); // in local camera space
}

bool fish_t::gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile) {

	assert(range.is_strictly_normalized());
	enabled = 0;
	if (water_is_lava || temperature <= W_FREEZE_POINT || temperature >= WATER_MAX_TEMP) return 0; // too hot/cold for fish
	pos = rgen.gen_rand_cube_point_xy(range);
	float const mesh_height(get_mesh_zval_at_pos(tile)), depth(water_plane_z - mesh_height);
	if (depth < 0.1) return 0; // no water
	radius  = FISH_RADIUS*rgen.rand_uniform(0.4, 1.0);
	float const fzmin(mesh_height + 1.6*get_half_height()), fzmax(water_plane_z - ocean_wave_height - get_tess_wave_height() - 2.0*get_half_height());
	if (fzmin > fzmax) return 0; // water is too shallow for this size of fish
	pos.z   = rgen.rand_uniform(fzmin, fzmax); // random depth within the valid range
	gen_dir_vel(rgen, FISH_SPEED);
	color   = WHITE;
	enabled = 1;
	return 1;
}

bool bird_t::gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile) {

	assert(range.is_strictly_normalized());
	if (atmosphere < 0.5) {enabled = 0; return 0;} // no atmosphere, no clouds, no birds
	pos     = rgen.gen_rand_cube_point(range);
	radius  = BIRD_RADIUS*rgen.rand_uniform(0.6, 1.0);
	gen_dir_vel(rgen, BIRD_SPEED);
	color   = BLACK;
	time    = rgen.rand_uniform(0.0, 100.0); // start at random time offsets
	enabled = 1;
	flocking= 0;
	return 1;
}

float get_butterfly_max_alt() {return 0.10f*(X_SCENE_SIZE + Y_SCENE_SIZE);}
float get_butterfly_min_alt() {return 0.10f*get_butterfly_max_alt();}

bool butterfly_t::gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile) { // Note: tile is unused

	assert(range.is_strictly_normalized());
	enabled = 0;
	if (atmosphere < 0.5) return 0; // no atmosphere, no clouds, no butterflies
	pos     = rgen.gen_rand_cube_point_xy(range);
	float const mesh_height(get_mesh_zval_at_pos(tile));
	if (mesh_height < water_plane_z) return 0; // no butterflies over water
	if (get_rel_height(mesh_height, -zmax_est, zmax_est) > 0.6) return 0; // maybe not necessary, since we know the tile must contain some grass
	pos.z   = mesh_height + rgen.rand_uniform(get_butterfly_min_alt(), get_butterfly_max_alt()); // random amount above the mesh
	radius  = BFLY_RADIUS*rgen.rand_uniform(0.8, 1.0)/tree_scale; // scale with tree size
	gen_dir_vel(rgen, BFLY_SPEED/tree_scale);
	unsigned const NUM_COLORS = 10;
	colorRGBA const colors[NUM_COLORS] = {WHITE, WHITE, WHITE, WHITE, ALPHA0, ALPHA0, ORANGE, YELLOW, LT_BLUE, colorRGBA(0.8, 0.6, 1.0)}; // Note: WHITE==default texture
	color   = colors[rgen.rand()%NUM_COLORS];
	//color   = WHITE; // textured, not colored
	time    = rgen.rand_uniform(0.0, 100.0); // start at random time offsets
	gender  = rgen.rand_bool(); // 50% male, 50% female
	enabled = 1;
	return 1;
}


bool fish_t::update(rand_gen_t &rgen, tile_t const *const tile) {

	if (!enabled || !animate2) return 0;
	point const camera(get_camera_pos()), pos_(get_camera_space_pos());
	if (!dist_less_than(pos_, camera, 200.0*radius)) return 1; // to far away to simulate (optimization)
	if (pos.z - 1.1*get_half_height() < get_mesh_zval_at_pos(tile)) {enabled = 0; anim_time = 0.0; return 0;}
	bool const chased(dist_less_than(pos_, camera, 15.0*radius));
	float const speed(velocity.mag());

	if (chased) { // scared by the player, swim away
		dir   = (pos_ - camera);
		dir.z = 0.0; // only swims in the xy plane
		dir.normalize();
		velocity = (10.0*FISH_SPEED) * dir; // swim away at a constant velocity
	}
	else if (!dist_less_than(pos_, camera, 20.0*radius)) { // far enough away from the player
		if (speed > FISH_SPEED) { // moving fast
			velocity *= pow(0.96f, fticks); // slow down
		}
		else if ((rand() & 127) == 0) { // randomly update direction
			dir += rgen.signed_rand_vector_xy(0.25); // 25% change max
			dir.normalize();
			velocity = dir * speed; // always flies in the direction it's pointed in
		}
	}
	point const orig_pos(pos);
	pos += velocity*fticks; // try to move
	
	if (pos.z - 1.5*get_half_height() < get_mesh_zval_at_pos(tile)) { // water is too shallow
		pos = orig_pos;
		if (chased) {velocity = zero_vector; anim_time = 0.0;} // stop
		else {gen_dir_vel(rgen, FISH_SPEED);} // pick a new direction randomly
		// or swim up?
	}
	anim_time += speed*fticks;
	if (anim_time > 10000) {anim_time -= 10000;} // reset after many ticks to avoid FP error with large numbers
	return 1;
}

bool bird_t::update(rand_gen_t &rgen, tile_t const *const tile) { // Note: tile is unused

	if (!enabled || !animate2 || !birds_active()) return 0;

	if (!flocking && (rand() & 1) == 0) { // randomly update direction
		float const speed(velocity.mag());
		dir += rgen.signed_rand_vector_xy(0.05); // 5% change max
		dir.normalize();
		velocity = dir * speed; // always flies in the direction it's pointed in
	}
	flocking = 0; // reset for next frame
	pos  += velocity*fticks; // always moving
	time += fticks;
	if (time > 600*TICKS_PER_SECOND) {time = 0.0;} // reset every 10 min.
	return 1;
}

void bird_t::apply_force_xy_const_vel(vector3d const &force) {

	float const vmag(velocity.mag());
	velocity.x += force.x; velocity.y += force.y;
	velocity *= vmag/velocity.mag(); // re-normalize
	dir = velocity/vmag; // normalized
	flocking = 1;
}

void adj_tiles_t::ensure_valid(tile_t const *const tile) {
	if (valid) return;
	assert(tile != nullptr);
	tile_xy_pair const tp(tile->get_tile_xy_pair());
	unsigned ix(0);

	for (int dy = -1; dy <= 1; ++dy) {
		for (int dx = -1; dx <= 1; ++dx) {adj[ix++] = get_tile_from_xy(tile_xy_pair(tp.x + dx, tp.y + dy));}
	}
	valid = 1;
}

void vect_bird_t::flock(tile_t const *const tile, adj_tiles_t &adj_tiles) { // boids, called per-tile

	// see https://www.blog.drewcutchins.com/blog/2018-8-16-flocking
	if (!animate2 || this->empty()) return;
	float const neighbor_dist(0.5*get_tile_width()), nd_sq(neighbor_dist*neighbor_dist);
	float const sep_dist_sq(0.2*nd_sq), cohesion_dist_sq(0.3*nd_sq), align_dist_sq(0.25*nd_sq);
	float const dmax(sqrt(max(sep_dist_sq, max(cohesion_dist_sq, align_dist_sq))));
	float const mass(100.0), sep_strength(0.05), cohesion_strength(0.05), align_strength(0.5);
	adj_tiles.ensure_valid(tile);

	for (auto i = this->begin(); i != this->end(); ++i) {
		if (!i->is_enabled()) continue;
		point const cs_pos(i->get_camera_space_pos());
		vector3d avg_pos(zero_vector), avg_vel(zero_vector), tot_force(zero_vector);
		unsigned pcount(0), vcount(0);

		for (unsigned adj_ix = 0; adj_ix < 9; ++adj_ix) {
			tile_t *const adj_tile(adj_tiles.adj[adj_ix]);
			if (!adj_tile) continue;
			vect_bird_t &birds(adj_tile->get_birds());
			if (!birds.empty() && adj_tile != tile && !adj_tile->get_mesh_bcube().closest_dist_xy_less_than(cs_pos, dmax)) continue; // tile is too far away

			for (auto j = birds.begin(); j != birds.end(); ++j) {
				if (!j->is_enabled()) continue;
				if (i == j) continue; // skip self
				float const dxy_sq(p2p_dist_xy_sq(i->pos, j->pos)); // Note: ignores zval
				if (dxy_sq < sep_dist_sq     ) {tot_force += (i->pos - j->pos)*(sep_strength/dxy_sq);} // separation force decreases with distance
				if (dxy_sq < cohesion_dist_sq) {avg_pos   += j->pos;      ++pcount;}
				if (dxy_sq < align_dist_sq   ) {avg_vel   += j->velocity; ++vcount;}
			} // for j
		} // for adj_ix
		if (pcount > 0) {tot_force += (avg_pos/pcount - i->pos)*cohesion_strength;} // cohesion
		if (vcount > 0) {tot_force += avg_vel*(align_strength/vcount);} // alignment
		if (tot_force != zero_vector) {i->apply_force_xy_const_vel(tot_force/mass);}
	} // for i
}

void update_accel(float &accel, rand_gen_t &rgen) {accel = CLIP_TO_pm1(accel + 0.25f*fticks*rgen.signed_rand_float());}

bool butterfly_t::update(rand_gen_t &rgen, tile_t const *const tile) {

	if (!enabled || !animate2 || !birds_active()) return 0;

	if (rest_time > 0.0) { // resting
		rest_time   -= fticks;
		if (rest_time > 0.0) return 1; // still resting
		rest_time    = 0.0; // done resting
		alt_change   = 0.5; // make sure we lift off into the air
		explore_time = TICKS_PER_SECOND*rgen.rand_uniform(15.0, 30.0); // free roam for 15-30s; this allows us to cross into a different tile
	}
	mate_time = max(0.0f, (mate_time - fticks));
	point const prev_pos(pos), prev_cs_pos(get_camera_space_pos());
	float const update_factor(0.01f*fticks);
	update_accel(fwd_accel, rgen);
	update_accel(rot_accel, rgen);
	update_accel(alt_accel, rgen);
	speed_factor = min(1.5f, max( 0.5f, (speed_factor + update_factor*fwd_accel)));
	rot_rate     = min(1.0f, max(-1.0f, (rot_rate     + update_factor*rot_accel)));
	alt_change   = min(1.0f, max(-1.0f, (alt_change   + update_factor*alt_accel)));
	float const delta_t(speed_factor*fticks), vmag(velocity.mag());
	float const rot_angle(0.0005*TWO_PI*delta_t*rot_rate);
	vector3d const prev_dir(dir);
	do_xy_rotate_normal(sin(rot_angle), cos(rot_angle), dir);
	float min_alt_factor(1.0);

	if (dest_valid) { // if we have a destination, bias our direction and altitude toward it based on how far away it is
		vector3d dest_dir(cur_dest - pos);
		float const rel_xy_dist(dest_dir.xy_mag()/radius);
		dest_alignment = 1.0 - min(0.005f*rel_xy_dist, 1.0f); // [0.0, 1.0]
		min_eq(min_alt_factor, 0.1f*rel_xy_dist); // allow low altitude flight, assuming destination is near the ground
		dest_dir.normalize();
		float const prev_dp(dot_product(prev_dir, dest_dir)), cur_dp(dot_product(dir, dest_dir));
		float const min_dp(0.5 + 0.5*dest_alignment); // [0.5, 1.0]
		
		if (cur_dp < min_dp) { // face us toward our destination
			if (prev_dp >= min_dp) {dir = prev_dir;} // revert to previous dir
			else if (dest_dir != zero_vector) { // slowly merge to dest_dir, if valid
				float const blend(pow(0.1, delta_t)*min_dp);
				dir = ((1.0 - blend)*dir + blend*dest_dir).get_norm();
			}
		}
		if (dest_alignment > 0.5) { // getting within range, start to adjust altitude
			float const v(dest_alignment - 0.5);
			alt_change = CLIP_TO_pm1(alt_change + 2.0f*v*v*SIGN(dest_dir.z)*min(fabs(dest_dir.z)/radius, 1.0f));
		}
	}
	else {
		dest_alignment = 0.0;
	}
	float const delta_z(0.4*alt_change*delta_t*radius), dir_mag(dir.mag());
	if (dir_mag < TOLERANCE) {dir = plus_x; velocity = vmag*dir;} // avoid divide-by-zero; unclear if this can happen
	else {velocity = (vmag/dir_mag)*dir;}
	pos   += velocity*delta_t;
	time  += delta_t; // controls wing speed
	if (dest_alignment > 0.75 && fabs(pos.z - cur_dest.z) < fabs(delta_z)) {pos.z = cur_dest.z;} // use exact alignment if possible
	else {pos.z += delta_z;}
	if (time > 600*TICKS_PER_SECOND && !is_visible(get_camera_space_pos(), 0.4)) {time = 0.0;} // reset every 10 min. if not visible
	float const coll_radius(2.0*radius); // use a larger radius for a buffer
	float const mesh_height(get_mesh_zval_at_pos(tile));
	float const zmin_val(max(mesh_height, water_plane_z) + min_alt_factor*coll_radius); // keep within the correct altitude range
	float const max_dz(5.0*vmag*fticks); // 5x nominal forward velocity
	max_eq(pos.z, zmin_val); // required min altitude
	max_eq(pos.z, min((pos.z + max_dz), (zmin_val + min_alt_factor*get_butterfly_min_alt())));
	min_eq(pos.z, max((pos.z - max_dz), (zmin_val + get_butterfly_max_alt())));
	point cs_pos(get_camera_space_pos());
	vector3d cnorm;
	// skip collision checks 7/8 times when not visible to the player (optimization); prev_cs_pos won't be accurate, but collisions should be close enough;
	// we could cache the prev_pos from the most recent collision check, but there would be a jump when the butterfly moves to a new tile or the player's origin is shifted
	unsigned const update_ix(gender + 2*(dir.x < 0) + 4*(dir.y < 0)); // 0-8 based on somewhat random but uniformly distributed state
	bool const skip_coll_check(((update_ix + frame_counter)&7) && !is_visible(cs_pos, 0.4));
	// skip tile scenery intersection check if close to the destination or within the dest plant bsphere
	bool const skip_tile_int_check(dest_valid && !is_mating && (dist_less_than(pos, cur_dest, 2.0*radius) ||
		(dest_bsphere.radius > 0.0 && dist_less_than(pos, dest_bsphere.pos, (coll_radius + dest_bsphere.radius)))));
	bool check_tree_coll(dist_less_than(cs_pos, get_camera_pos(), (X_SCENE_SIZE + Y_SCENE_SIZE)));
	
	// skip building interiors (shouldn't be there) and cars (too slow)
	if (!skip_coll_check && proc_city_sphere_coll(cs_pos, prev_cs_pos, coll_radius, prev_cs_pos.z, 0, &cnorm, 0)) {
		pos = cs_pos - get_camera_coord_space_xlate(); // back to world space
		calc_reflection_angle(dir, dir, cnorm); // reflect
		dir.normalize();
		velocity     = dir*vmag; // change direction but preserve velocity
		dest_valid   = is_mating = 0; // choose a new destination, in case this one is blocked by a building
		explore_time = TICKS_PER_SECOND*rgen.rand_uniform(2.0, 5.0); // explore a bit more to get out from between the buildings
	}
	else if ((mesh_height < water_plane_z - 0.5*get_butterfly_max_alt()) || // over deep water
		// check collision with trees (maybe) and scenery
		(!skip_tile_int_check && !skip_coll_check && tile && tile->check_sphere_collision(cs_pos, coll_radius, check_tree_coll, check_tree_coll, 1)))
	{
		dir.negate(); // just negate the direction because we don't have the collision normal
		velocity.negate();
	}
	else if (pos.x != prev_pos.x || pos.y != prev_pos.y) {
		dir   = (pos - prev_pos); // align direction to velocity vector
		dir.z = 0.0; // always level
		dir.normalize();
	}
	update_dest(rgen, tile);
	return 1;
}

void butterfly_t::update_dest(rand_gen_t &rgen, tile_t const *const tile) {
	if (explore_time > 0.0) { // exploring
		explore_time -= fticks;
		if (explore_time > 0.0) return; // still exploring
		explore_time = 0.0; // done exploring, find a destination
	}
	if (dest_valid) {
		if (dist_less_than(pos, cur_dest, radius)) { // at destination, including mating
			if (is_mating) {mate_time = TICKS_PER_SECOND*rgen.rand_uniform(15.0, 30.0);} // done mating, wait before mating again
			else           {rest_time = TICKS_PER_SECOND*rgen.rand_uniform( 5.0, 10.0);} // rest for 5-10s if not mating
			prev_dest  = cur_dest;
			cur_dest   = pos;
			dest_valid = is_mating = 0;
		}
		return; // no dest update
	}
	dest_bsphere.radius = 0.0;
	if (rgen.rand_float() < 0.9) return; // only look for a destination 10% of the time as an optimization
	if (have_cities() && rgen.rand_bool() && choose_pt_in_city_park(pos, cur_dest, rgen)) {} // choose a city park
	else if (tile && tile->choose_butterfly_dest(cur_dest, dest_bsphere, rgen)) { // choose a new destination within this tile
		cur_dest.z -= 0.5*radius; // shift slightly downward since the body is below the butterfly center
	}
	else {return;} // no dest
	// Note: we're picky about the dest here; many dests will be rejected, and we'll get into this code again in a later frame; this gives butterflies more time to explore
	if (cur_dest == prev_dest) return; // same as previous dest, ignore
	if (dot_product_ptv(dir, cur_dest, pos) < 0.0) return; // behind us, ignore
	float t(0.0); // unused
	int xpos(0), ypos(0); // unused
	if (tile && tile->line_intersect_mesh(get_camera_space_pos(), get_camera_space_dest(), t, xpos, ypos, 0)) return; // blocked by mesh; exclude pine trees
	dest_valid = 1; // success
}

bool butterfly_t::can_mate_with(butterfly_t const &b) const {return (gender != b.gender && color == b.color);} // different gender and same color

void vect_butterfly_t::run_mating(tile_t const *const tile, adj_tiles_t &adj_tiles) {
	if (!animate2 || this->empty()) return;
	adj_tiles.ensure_valid(tile);
	float const mate_dmax(0.75f*(X_SCENE_SIZE + Y_SCENE_SIZE)); // 0.75 tile

	for (auto i = this->begin(); i != this->end(); ++i) {
		if (!i->is_enabled() || i->mate_time > 0.0 || i->rest_time > 0 || i->explore_time > 0) continue; // skip if waiting to mate, resting, or exploring
		if (!debug_animal_draw() && rgen.rand_float() < 0.9) continue; // only run 10% of the time; skip this optimization for smooth debug draw
		if (i->is_mating) {i->dest_valid = i->is_mating = 0;} // reset for this iteration
		point const cs_pos(i->get_camera_space_pos());
		float dmin_sq(mate_dmax*mate_dmax);

		for (unsigned adj_ix = 0; adj_ix < 9; ++adj_ix) {
			tile_t *const adj_tile(adj_tiles.adj[adj_ix]);
			if (!adj_tile) continue;
			vect_butterfly_t &bflies(adj_tile->get_bflies());
			if (!bflies.empty() && adj_tile != tile && !adj_tile->get_mesh_bcube().closest_dist_xy_less_than(cs_pos, mate_dmax)) continue; // tile is too far away

			for (auto j = bflies.begin(); j != bflies.end(); ++j) {
				if (!j->is_enabled() || i == j || !i->can_mate_with(*j)) continue; // skip self or same gender
				if (j->rest_time > 0) continue; // skip another butterfly resting on a plant
				float const dxy_sq(p2p_dist_xy_sq(i->pos, j->pos)); // Note: ignores zval
				if (dxy_sq > dmin_sq) continue; // not closer
				if ((j->cur_dest != i->pos) && dot_product_ptv(i->dir, j->pos, i->pos) < 0.0) continue; // behind us, ignore unless this one wants to mate us (works half the time)
				point const cs_pos(i->get_camera_space_pos()), cs_dest(j->get_camera_space_pos());
				float t(0.0); // unused
				int xpos(0), ypos(0); // unused
				// check both tiles (approximate)
				if (                        tile->line_intersect_mesh(cs_pos, cs_dest, t, xpos, ypos, 0)) continue; // blocked by mesh; exclude pine trees
				if (adj_tile != tile && adj_tile->line_intersect_mesh(cs_pos, cs_dest, t, xpos, ypos, 0)) continue; // blocked by mesh; exclude pine trees
				i->cur_dest  = j->pos;
				i->is_mating = i->dest_valid = 1;
				dmin_sq      = dxy_sq;
			} // for j
		} // for adj_ix
	} // for i
}


bool animal_t::distance_check(point const &pos_, float vis_dist_scale) const {
	return dist_less_than(pos_, get_camera_pos(), 1000.0*vis_dist_scale*radius);
}
bool animal_t::is_visible(point const &pos_, float vis_dist_scale) const {
	if (!enabled) return 0;
	if (!distance_check(pos_, vis_dist_scale)) return 0;
	return sphere_in_camera_view(pos_, radius, 0);
}

int animal_t::get_ndiv(point const &pos_) const {
	return min(N_SPHERE_DIV, max(3, int(4.0*sqrt(radius*window_width/distance_to_camera(pos_)))));
}

void draw_animated_fish_model(shader_t &s, vector3d const &pos, float radius, vector3d const &dir, float anim_time, colorRGBA const &color) {
	if (!s.is_setup()) { // shader not setup; not needed if using lighting shader for fish in buildings
		s.set_vert_shader("fish_animate");
		s.set_frag_shader("simple_texture");
		s.begin_shader();
		s.add_uniform_float("min_alpha", 0.0);
		s.add_uniform_int("tex0", 0);
	}
	animal_model_loader.draw_fish_model(s, pos, radius, dir, anim_time, color);
}

void fish_t::draw(shader_t &s, tile_t const *const tile, bool &first_draw) const { // Note: tile is unused, but could be used for shadows

	point const pos_(get_camera_space_pos());
	if (!is_visible(pos_, 0.15)) return;
	colorRGBA draw_color(color);
	water_color_atten_at_pos(draw_color, pos_ - vector3d(0.0, 0.0, radius)); // move down slightly
	point const camera(get_camera_pos());
	float const t(CLIP_TO_01((water_plane_z - pos_.z)/max(1.0E-6f, fabs(camera.z - pos_.z))));
	point const p_int(pos_ + (camera - pos_)*t); // ray-water intersection point
	float const dist(p2p_dist(p_int, camera) + 2.0*p2p_dist(p_int, pos_)); // water is 2x as optically dense as air
	draw_color.alpha = min(1.0f, 1.5f*exp(-0.028f*dist/radius));
	if (draw_color.alpha < 0.01) return;
	if (draw_color.alpha < 0.1) {glDepthMask(GL_FALSE);} // disable depth writing to avoid alpha blend order problems
	//draw_color = lerp(cur_fog_color, draw_color, alpha);
	draw_animated_fish_model(s, pos_, radius, dir, anim_time, draw_color);
	if (draw_color.alpha < 0.1) {glDepthMask(GL_TRUE);}
}

void bird_t::draw(shader_t &s, tile_t const *const tile, bool &first_draw) const { // Note: tile is unused

	point const pos_(get_camera_space_pos());
	if (!is_visible(pos_, 1.0)) return;
	// Note: birds use distance-based transparency rather than fog, because they may be against the blue sky above rather than the distant gray fog on the horizon;
	// also, fog is modeled to be lower to the ground, and doesn't affect things like the sky, clouds, and sun
	float const dist(p2p_dist(pos_, get_camera_pos()));
	float const alpha(CLIP_TO_01(2000.0f*radius/dist - 2.0f)); // 1.0 under half clip distance, after that linear falloff to zero
	if (alpha < 0.01) return;
	
	if (!s.is_setup()) {
		s.begin_color_only_shader();
		bind_draw_sphere_vbo(0, 0); // no textures or normals
	}
	s.set_cur_color(colorRGBA(color, alpha));
	int const ndiv(get_ndiv(pos_));
	fgPushMatrix();
	translate_to(pos_);
	rotate_to_plus_x(dir);
	scale_by(vector3d(1.0, 0.2, 0.14)*radius); // x = length, y = width, z = height
	draw_sphere_vbo_pre_bound(ndiv, 0); // body
	fgTranslate(0.1, 0.0, 0.0);
	float const angle(50.0*(sin(0.2*time) + 0.3));

	for (unsigned d = 0; d < 2; ++d) {
		float const v(d ? 1.0 : -1.0);
		fgPushMatrix();
		fgRotate(v*angle, 1.0, 0.0, 0.0); // rotate about x-axis to flap wings
		fgScale(0.3, 5.0, 0.5);
		fgTranslate(0.0, 0.8*v, 0.0);
		draw_sphere_vbo_pre_bound(ndiv, 0); // wings
		fgPopMatrix();
	}
	fgPopMatrix();
}

void butterfly_t::draw(shader_t &s, tile_t const *const tile, bool &first_draw) const {

	if (!enabled) return;
	point const pos_(get_camera_space_pos());
	bool const visible(is_visible(pos_, 0.4));
	bool const debug_lines(debug_animal_draw()), debug_spheres(0/*display_mode & 0x20*/);
	if (!debug_lines && !debug_spheres && !visible) return;
	
	if (!s.is_setup()) {
		s.set_prefix("#define TWO_SIDED_LIGHTING", 1); // FS
		tile_draw_t::tree_branch_shader_setup(s, shadow_map_enabled(), 0, 0, 0);
	}
	if (first_draw && tile != nullptr) {
		tile->bind_and_setup_shadow_map(s);
		first_draw = 0;
	}
	if (visible) { // draw butterfly
		bool const draw_body(distance_check(pos_, 0.15)); // draw body if close enough to the player
		animal_model_loader.draw_butterfly_model(s, pos_, radius, dir, time/TICKS_PER_SECOND, draw_body, color);
	}
	if (debug_lines) { // debug draw line to destination
		vector<vert_norm_color> line_pts;
		select_no_texture();

		if (dest_valid) {
			colorRGBA const color(blend_color((is_mating ? YELLOW : RED), (is_mating ? ORANGE : BLUE), dest_alignment, 0));
			point const cs_dest(get_camera_space_dest());
			line_pts.emplace_back(pos_,    plus_z, color);
			line_pts.emplace_back(cs_dest, plus_z, color);

			if (!is_mating) { // draw destination point
				s.set_cur_color(GREEN);
				draw_sphere_vbo(cs_dest, 0.5*radius, 16, 0);
			}
		}
		if (rest_time == 0) { // not resting; draw line pointing in current direction
			line_pts.emplace_back( pos_, plus_z, GREEN);
			line_pts.emplace_back((pos_ + 4.0*radius*dir), plus_z, GREEN);
		}
		draw_verts(line_pts, GL_LINES);
	}
	if (debug_spheres && distance_check(pos_, 1.0)) { // debug draw path as spheres
		if (path.empty() || !dist_less_than(pos, path.back(), 2.0*radius)) {path.push_back(pos);}
		select_no_texture();
		s.set_cur_color(WHITE);
		begin_sphere_draw(0); // untextured
		for (auto const &p : path) {draw_sphere_vbo((p + get_camera_coord_space_xlate()), 0.25*radius, 8, 0);}
		end_sphere_draw();
	}
	else {path.clear();}
}


template<typename A> void animal_group_t<A>::gen(unsigned num, cube_t const &range, tile_t const *const tile) { // Note: okay if nonempty

	assert(tile);
	generated = 1;
	if (!range.is_strictly_normalized()) return; // if the range is empty or denormalized in z, mark as generated but skip
	if (!A::type_enabled() || num == 0)  return; // model not valid/loaded, or no objects enabled
	if (!A::can_place_in_tile(tile))     return; // tile not valid for this object
	rgen.set_state(long(1000*range.d[0][0]), long(1000*range.d[1][0])); // seed with x1, y1
	this->reserve(this->size() + num); // optional

	for (unsigned n = 0; n < num; ++n) {
		A animal{};
		if (animal.gen(rgen, range, tile)) {this->push_back(animal);} // only add if generation was successful
	}
	bcube = range; // initial value (approx)
}

template<typename A> void animal_group_t<A>::update(tile_t const *const tile) {

	if (!animate2) return;
	bool bcube_set(0);
	bcube.set_from_point(all_zeros);

	for (auto i = this->begin(); i != this->end(); ++i) {
		i->update(rgen, tile);
		if (!i->is_enabled()) continue;
		if (!bcube_set) {bcube.set_from_sphere(*i); bcube_set = 1;}
		else {bcube.union_with_sphere(*i);}
	}
}

template<typename A> void animal_group_t<A>::remove(unsigned ix) {
	assert(ix < this->size());
	std::swap(this->operator[](ix), this->back());
	this->pop_back();
}

template<typename A> void animal_group_t<A>::draw_animals(shader_t &s, tile_t const *const tile) const {
	if (this->empty() || bcube.is_zero_area()) return;
	if (!debug_animal_draw() && !camera_pdu.cube_visible(bcube + get_camera_coord_space_xlate())) return;
	bool first_draw(1);
	for (auto i = this->begin(); i != this->end(); ++i) {i->draw(s, tile, first_draw);}
	bind_vbo(0); // needed for birds, okay for others
}

// explicit instantiations
template class animal_group_t<fish_t>;
template class animal_group_t<bird_t>;
template class animal_group_t<butterfly_t>;

