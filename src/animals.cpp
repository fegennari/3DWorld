// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#include "animals.h"
#include "shaders.h"

float const FISH_RADIUS = 0.05;
float const BIRD_RADIUS = 0.1;
float const FISH_SPEED  = 0.002;
float const BIRD_SPEED  = 0.05;

extern bool water_is_lava;
extern int window_width, animate2;
extern float fticks, water_plane_z, temperature, atmosphere, ocean_wave_height;


void animal_t::gen_dir_vel(rand_gen_t &rgen, float speed) {
	dir      = rgen.signed_rand_vector_xy().get_norm(); // in xy plane
	velocity = (speed*rgen.rand_uniform(0.6, 1.0))*dir;
}

// FIXME: slow in cases where a noise function is evaluated - can we use the tile zvals instead?
float fish_t::get_mesh_zval_at_pos() const {
	return interpolate_mesh_zval((pos.x - DX_VAL*xoff2), (pos.y - DY_VAL*yoff2), 0.0, 1, 1); // in local camera space
}

bool fish_t::gen(rand_gen_t &rgen, cube_t const &range) {

	assert(range.is_strictly_normalized());
	enabled = 0;
	if (water_is_lava || temperature <= W_FREEZE_POINT || temperature >= WATER_MAX_TEMP) return 0; // too hot/cold for fish
	pos = rgen.gen_rand_cube_point(range);
	float const mesh_height(get_mesh_zval_at_pos()), depth(water_plane_z - mesh_height);
	if (depth < 0.1) return 0; // no water
	radius  = FISH_RADIUS*rgen.rand_uniform(0.4, 1.0);
	scale   = vector3d(1.0, 0.24*rgen.rand_uniform(0.6, 1.0), 0.4*rgen.rand_uniform(0.8, 1.0)); // x = length, y = width, z = height
	float const fzmin(mesh_height + 1.6*radius*scale.z), fzmax(water_plane_z - ocean_wave_height - get_tess_wave_height() - 2.0*radius*scale.z);
	if (fzmin > fzmax) return 0; // water is too shallow for this size of fish
	pos.z   = rgen.rand_uniform(fzmin, fzmax); // random depth within the valid range
	gen_dir_vel(rgen, FISH_SPEED);
	color   = colorRGBA(rgen.rand_uniform(0.02, 0.03), rgen.rand_uniform(0.05, 0.07), rgen.rand_uniform(0.08, 0.10), 1.0);
	enabled = 1;
	//tile_off.set_from_xyoff2(); // Note: TT specific
	return 1;
}

bool bird_t::gen(rand_gen_t &rgen, cube_t const &range) {

	assert(range.is_strictly_normalized());
	enabled = 0;
	if (atmosphere < 0.5) return 0; // no atmosphere, no clouds, no birds
	pos     = rgen.gen_rand_cube_point(range);
	radius  = BIRD_RADIUS*rgen.rand_uniform(0.6, 1.0);
	gen_dir_vel(rgen, BIRD_SPEED);
	color   = BLACK;
	time    = rgen.rand_uniform(0.0, 100.0); // start at random time offsets
	enabled = 1;
	//tile_off.set_from_xyoff2(); // Note: TT specific
	return 1;
}


bool fish_t::update(rand_gen_t &rgen) {

	if (!enabled || !animate2) return 0;
	point const camera(get_camera_pos()), pos_(get_draw_pos());
	if (!dist_less_than(pos_, camera, 200.0*radius)) return 1; // to far away to simulate (optimization)
	if (pos.z - 1.1*radius*scale.z < get_mesh_zval_at_pos()) {enabled = 0; return 0;}
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
	
	if (pos.z - 1.5*radius*scale.z < get_mesh_zval_at_pos()) { // water is too shallow
		pos = orig_pos;
		if (chased) {velocity = zero_vector;} // stop
		else {gen_dir_vel(rgen, FISH_SPEED);} // pick a new direction randomly
		// FIXME: swim up?
	}
	return 1;
}

bool bird_t::update(rand_gen_t &rgen) {

	if (!enabled || !animate2) return 0;

	if ((rand() & 3) == 0) { // randomly update direction
		float const speed(velocity.mag());
		dir += rgen.signed_rand_vector_xy(0.1); // 10% change max
		dir.normalize();
		velocity = dir * speed; // always flies in the direction it's pointed in
	}
	pos  += velocity*fticks; // always moving
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

vector3d get_pos_offset() {return vector3d((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0);}
point animal_t::get_draw_pos() const {return (pos + get_pos_offset());}

float animal_t::get_alpha_at_cur_dist(point const &pos_, float vis_dist_scale) const {
	return CLIP_TO_01(2000.0f*vis_dist_scale*radius/p2p_dist(pos_, get_camera_pos()) - 2.0f); // 1.0 under half clip distance, after that linear falloff to zero
}

void rotate_to_plus_x(vector3d const &dir) {
	rotate_about(TO_DEG*get_norm_angle(dir, plus_x), vector3d(0.0, 0.0, dir.y));
}

void fish_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
	if (!is_visible(pos_, 0.06)) return;
	colorRGBA draw_color(color);
	water_color_atten_at_pos(draw_color, pos_);
	float const alpha(get_alpha_at_cur_dist(pos_, 0.06)); // FIXME: closer to underwater fog
	if (alpha < 0.1) return;
	s.set_cur_color(colorRGBA(draw_color, alpha));
	fgPushMatrix();
	translate_to(pos_);
	rotate_to_plus_x(dir);
	scale_by(radius*scale);
	draw_sphere_vbo_back_to_front(all_zeros, 1.0, get_ndiv(pos_), 0);
	fgPopMatrix();
}

void bird_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
	if (!is_visible(pos_, 1.0)) return;
	float const alpha(get_alpha_at_cur_dist(pos_, 1.0));
	if (alpha < 0.05) return;
	s.set_cur_color(colorRGBA(color, alpha)); // FIXME: fog
	int const ndiv(get_ndiv(pos_));
	fgPushMatrix();
	translate_to(pos_);
	rotate_to_plus_x(dir);
	uniform_scale(radius);
	fgScale(1.0, 0.2, 0.14); // x = length, y = width, z = height
	draw_sphere_vbo(all_zeros, 1.0, ndiv, 0); // body
	// wings - FIXME: flap up and down in a v-shape using time
	fgTranslate(0.1, 0.0, 0.0);
	fgScale(0.3, 5.0, 0.5);
	fgTranslate(0.0,  0.8, 0.0);
	draw_sphere_vbo(all_zeros, 1.0, ndiv, 0); // wing
	fgTranslate(0.0, -1.6, 0.0);
	draw_sphere_vbo(all_zeros, 1.0, ndiv, 0); // wing
	fgPopMatrix();
}


/*static*/ void animal_group_base_t::begin_draw(shader_t &s) {
	s.begin_color_only_shader(); // FIXME: placeholder
	enable_blend(); // for distance fog
}
/*static*/ void animal_group_base_t::end_draw(shader_t &s) {
	disable_blend(); // for distance fog
	s.end_shader();
}


template<typename A> void animal_group_t<A>::gen(unsigned num, cube_t const &range) { // Note: okay if nonempty

	generated = 1;
	if (!range.is_strictly_normalized()) return; // if the range is empty or denormalized in z, mark as generated but skip
	rgen.set_state(long(1000*range.d[0][0]), long(1000*range.d[1][0])); // seed with x1, y1
	reserve(size() + num); // optional

	for (unsigned n = 0; n < num; ++n) {
		A animal;
		if (animal.gen(rgen, range)) {push_back(animal);} // only add if generation was successful
	}
}

template<typename A> void animal_group_t<A>::update() {

	if (!animate2) return;
	bool bcube_set(0);
	bcube.set_from_point(all_zeros);

	for (iterator i = begin(); i != end(); ++i) {
		i->update(rgen);
		if (!i->is_enabled()) continue;
		if (!bcube_set) {bcube.set_from_sphere(*i); bcube_set = 1;}
		else {bcube.union_with_sphere(*i);}
	}
}

template<typename A> void animal_group_t<A>::remove(unsigned ix) {
	assert(ix < size());
	std::swap(operator[](ix), back());
	pop_back();
}

template<typename A> void animal_group_t<A>::remove_disabled() {
	iterator i(begin()), o(i);
	for (; i != end(); ++i) {if (i->is_enabled()) {*(o++) = *i;}}
	erase(o, end());
}

template<typename A> void animal_group_t<A>::draw_animals(shader_t &s) const {
	if (bcube.is_zero_area() || !camera_pdu.cube_visible(bcube + get_pos_offset())) return;
	for (const_iterator i = begin(); i != end(); ++i) {i->draw(s);}
}


void vect_fish_t::draw() const {

	shader_t s;
	begin_draw(s);
	draw_animals(s);
	end_draw(s);
}

void vect_bird_t::draw() const {

	shader_t s;
	begin_draw(s);
	draw_animals(s);
	end_draw(s);
}


// explicit instantiations
template class animal_group_t<fish_t>;
template class animal_group_t<bird_t>;

