// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#include "animals.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "model3d.h"
#include "tiled_mesh.h"

float const FISH_RADIUS = 0.05;
float const BIRD_RADIUS = 0.1;
float const FISH_SPEED  = 0.002;
float const BIRD_SPEED  = 0.05;

extern bool water_is_lava;
extern int window_width, animate2, display_mode;
extern float fticks, water_plane_z, temperature, atmosphere, ocean_wave_height;
extern colorRGBA cur_fog_color;


bool birds_active() {return (light_factor >= 0.4);} // birds are only active whe the sun is out


class animal_model_loader_t : public model3ds { // currently for fish only

	unsigned fish_id;

public:
	animal_model_loader_t() : fish_id(0) {}

	unsigned load_model(string const &fn, int recalc_normals=1, colorRGBA const &def_color=WHITE, int def_tid=-1) {
		unsigned const id(size());

		if (!load_model_file(fn, *this, geom_xform_t(), def_tid, def_color, 0, 0.0, recalc_normals, 0, 0, 1)) {
			cerr << "Error: Failed to read model file '" << fn << "'" << endl;
			exit(1);
		}
		return id + 1;
	}
	void draw_model(unsigned id, shader_t &s, vector3d const &pos, float radius, vector3d const &dir, rotation_t const &local_rotate=rotation_t(), colorRGBA const &color=WHITE, bool is_shadow_pass=0) {
		assert(id > 0 && id-1 < size()); // Note: id is vector index offset by 1
		bool const camera_pdu_valid(camera_pdu.valid);
		camera_pdu.valid = 0; // disable VFC, since we're doing custom transforms here
		s.add_uniform_color("color_modulate", color);
		model3d &model(operator[](id-1));
		model.bind_all_used_tids();
		cube_t const &bcube(model.get_bcube());
		fgPushMatrix();
		translate_to(pos);
		rotate_to_plus_x(dir);
		local_rotate.apply_gl();
		uniform_scale(radius / (0.5*bcube.max_len()));
		translate_to(-bcube.get_cube_center()); // cancel out model local translate
		model.render_materials(s, is_shadow_pass, 0, 0, 1, 3, 3, model.get_unbound_material(), rotation_t(), nullptr);
		fgPopMatrix();
		camera_pdu.valid = camera_pdu_valid;
	}
	void draw_fish_model(shader_t &s, vector3d const &pos, float radius, vector3d const &dir, colorRGBA const &color=WHITE) {
		if (fish_id == 0) { // fish model not yet loaded
			int const recalc_normals = 0; // okay for loading model3d
			string const fish_model_fn( "model_data/fish/fishOBJ.model3d" ); // provide in config file?
			fish_id = load_model(fish_model_fn, recalc_normals);
		}
		assert(fish_id > 0);
		vector3d const fish_dir(-dir); // fish model faces in -x, and we want it to be in +x
		rotation_t const local_rotate(plus_y, -45.0); // fish model is angled 45 degree upward, so need to rotate it back down
		draw_model(fish_id, s, pos, radius, fish_dir, local_rotate, color, 0); // not shadow pass
	}
};

animal_model_loader_t animal_model_loader;

void free_animal_context() {animal_model_loader.free_context();}


void animal_t::gen_dir_vel(rand_gen_t &rgen, float speed) {
	dir      = rgen.signed_rand_vector_xy().get_norm(); // in xy plane
	velocity = (speed*rgen.rand_uniform(0.6, 1.0))*dir;
}

float fish_t::get_mesh_zval_at_pos(tile_t const *const tile) const {
	// use tile to interpolate z-value if present; faster, but should return a similar value to interpolate_mesh_zval()
	if (tile) {return tile->get_zval_at(pos.x, pos.y, 1);} // in global space
	return interpolate_mesh_zval((pos.x - DX_VAL*xoff2), (pos.y - DY_VAL*yoff2), 0.0, 1, 1); // in local camera space
}

bool fish_t::gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile) {

	assert(range.is_strictly_normalized());
	enabled = 0;
	if (water_is_lava || temperature <= W_FREEZE_POINT || temperature >= WATER_MAX_TEMP) return 0; // too hot/cold for fish
	pos = rgen.gen_rand_cube_point(range);
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

bool bird_t::gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile) { // Note: tile is unused

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

bool butterfly_t::gen(rand_gen_t &rgen, cube_t const &range, tile_t const *const tile) { // Note: tile is unused

	assert(range.is_strictly_normalized());
	if (atmosphere < 0.5) {enabled = 0; return 0;} // no atmosphere, no clouds, no butterflies
	//pos     = ;
	//radius  = ;
	//gen_dir_vel(rgen, ?);
	color   = WHITE;
	time    = rgen.rand_uniform(0.0, 100.0); // start at random time offsets
	enabled = 1;
	return 1;
}


bool fish_t::update(rand_gen_t &rgen, tile_t const *const tile) {

	if (!enabled || !animate2) return 0;
	point const camera(get_camera_pos()), pos_(get_draw_pos());
	if (!dist_less_than(pos_, camera, 200.0*radius)) return 1; // to far away to simulate (optimization)
	if (pos.z - 1.1*get_half_height() < get_mesh_zval_at_pos(tile)) {enabled = 0; return 0;}
	bool const chased(dist_less_than(pos_, camera, 15.0*radius));

	if (chased) { // scared by the player, swim away
		dir   = (pos_ - camera);
		dir.z = 0.0; // only swims in the xy plane
		dir.normalize();
		velocity = (10.0*FISH_SPEED) * dir; // swim away at a constant velocity
	}
	else if (!dist_less_than(pos_, camera, 20.0*radius)) { // far enough away from the player
		float const speed(velocity.mag());
		
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
		if (chased) {velocity = zero_vector;} // stop
		else {gen_dir_vel(rgen, FISH_SPEED);} // pick a new direction randomly
		// or swim up?
	}
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
	return 1;
}

void bird_t::apply_force_xy_const_vel(vector3d const &force) {

	float const vmag(velocity.mag());
	velocity.x += force.x; velocity.y += force.y;
	velocity *= vmag/velocity.mag(); // re-normalize
	dir = velocity/vmag; // normalized
	flocking = 1;
}

void vect_bird_t::flock(tile_t const *const tile) { // boids, called per-tile

	// see https://www.blog.drewcutchins.com/blog/2018-8-16-flocking
	if (!animate2 || this->empty()) return;
	float const neighbor_dist(0.5*get_tile_width()), nd_sq(neighbor_dist*neighbor_dist);
	float const sep_dist_sq(0.2*nd_sq), cohesion_dist_sq(0.3*nd_sq), align_dist_sq(0.25*nd_sq);
	float const mass(100.0), sep_strength(0.05), cohesion_strength(0.05), align_strength(0.5);
	tile_xy_pair const tp(tile->get_tile_xy_pair());
	tile_t *adj_tiles[9] = {0};
	unsigned ix(0);

	for (int dy = -1; dy <= 1; ++dy) {
		for (int dx = -1; dx <= 1; ++dx) {adj_tiles[ix++] = get_tile_from_xy(tile_xy_pair(tp.x + dx, tp.y + dy));}
	}
	for (auto i = this->begin(); i != this->end(); ++i) {
		if (!i->is_enabled()) continue;
		vector3d avg_pos(zero_vector), avg_vel(zero_vector), tot_force(zero_vector);
		unsigned pcount(0), vcount(0);

		for (unsigned adj_ix = 0; adj_ix < 9; ++adj_ix) {
			tile_t *const adj_tile(adj_tiles[adj_ix]);
			if (!adj_tile) continue;
			vect_bird_t &birds(adj_tile->get_birds());

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

bool butterfly_t::update(rand_gen_t &rgen, tile_t const *const tile) { // Note: tile is unused

	if (!enabled || !animate2 || !birds_active()) return 0;
	pos  += velocity*fticks;
	time += fticks;
	return 1;
}


bool animal_t::is_visible(point const &pos_, float vis_dist_scale) const {

	if (!enabled) return 0;
	if (!dist_less_than(pos_, get_camera_pos(), 1000.0*vis_dist_scale*radius)) return 0;
	return sphere_in_camera_view(pos_, radius, 0);
}

int animal_t::get_ndiv(point const &pos_) const {
	return min(N_SPHERE_DIV, max(3, int(4.0*sqrt(radius*window_width/distance_to_camera(pos_)))));
}

point animal_t::get_draw_pos() const {return (pos + get_camera_coord_space_xlate());}

void fish_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
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
	animal_model_loader.draw_fish_model(s, pos_, radius, dir, draw_color);
	if (draw_color.alpha < 0.1) {glDepthMask(GL_TRUE);}
}

void bird_t::draw(shader_t &s) const {

	if (!birds_active()) return;
	point const pos_(get_draw_pos());
	if (!is_visible(pos_, 1.0)) return;
	// Note: birds use distance-based transparency rather than fog, because they may be against the blue sky above rather than the distant gray fog on the horizon;
	// also, fog is modeled to be lower to the ground, and doesn't affect things like the sky, clouds, and sun
	float const dist(p2p_dist(pos_, get_camera_pos()));
	float const alpha(CLIP_TO_01(2000.0f*radius/dist - 2.0f)); // 1.0 under half clip distance, after that linear falloff to zero
	if (alpha < 0.01) return;
	s.set_cur_color(colorRGBA(color, alpha));
	int const ndiv(get_ndiv(pos_));
	bind_draw_sphere_vbo(0, 0); // no textures or normals
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
	bind_vbo(0);
}

void butterfly_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
	if (!is_visible(pos_, 0.15)) return;
	animal_model_loader.draw_fish_model(s, pos_, radius, dir, color);
}


/*static*/ void vect_fish_t::begin_draw(shader_t &s) {
	s.begin_simple_textured_shader(); // no lighting
	enable_blend(); // for distance fog
}
/*static*/ void vect_bird_t::begin_draw(shader_t &s) {
	s.begin_color_only_shader();
	enable_blend(); // for distance fog
}
/*static*/ void vect_butterfly_t::begin_draw(shader_t &s) {
	s.begin_simple_textured_shader(); // no lighting
	enable_blend(); // for distance fog
}

/*static*/ void vect_fish_t::end_draw(shader_t &s) {
	disable_blend(); // for distance fog
	s.add_uniform_color("color_modulate", WHITE); // reset
	s.end_shader();
}
/*static*/ void vect_bird_t::end_draw(shader_t &s) {
	disable_blend(); // for distance fog
	s.end_shader();
}
/*static*/ void vect_butterfly_t::end_draw(shader_t &s) {
	disable_blend(); // for distance fog
	s.end_shader();
}


template<typename A> void animal_group_t<A>::gen(unsigned num, cube_t const &range, tile_t const *const tile) { // Note: okay if nonempty

	generated = 1;
	if (!range.is_strictly_normalized()) return; // if the range is empty or denormalized in z, mark as generated but skip
	rgen.set_state(long(1000*range.d[0][0]), long(1000*range.d[1][0])); // seed with x1, y1
	this->reserve(this->size() + num); // optional

	for (unsigned n = 0; n < num; ++n) {
		A animal;
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

template<typename A> void animal_group_t<A>::remove_disabled() {
	auto i(this->begin()), o(i);
	for (; i != this->end(); ++i) {if (i->is_enabled()) {*(o++) = *i;}}
	this->erase(o, this->end());
}

template<typename A> void animal_group_t<A>::draw_animals(shader_t &s) const {
	if (this->empty() || bcube.is_zero_area() || !camera_pdu.cube_visible(bcube + get_camera_coord_space_xlate())) return;
	for (auto i = this->begin(); i != this->end(); ++i) {i->draw(s);}
}

// explicit instantiations
template class animal_group_t<fish_t>;
template class animal_group_t<bird_t>;
template class animal_group_t<butterfly_t>;

