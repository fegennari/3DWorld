// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#include "animals.h"
#include "shaders.h"

extern bool water_is_lava;
extern int window_width, animate2;
extern float fticks, water_plane_z, temperature, atmosphere;


bool fish_t::gen(rand_gen_t &rgen, cube_t const &range) {

	// FIXME: shallow water fish?
	assert(range.is_strictly_normalized());
	enabled = 0;
	if (water_is_lava || temperature <= W_FREEZE_POINT || temperature >= WATER_MAX_TEMP) return 0; // too hot/cold for fish
	pos = rgen.gen_rand_cube_point(range);
	float const mesh_height(interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1)), depth(water_plane_z - mesh_height);
	if (depth < 0.1) return 0; // no water
	radius  = 0.05*rgen.rand_uniform(0.4, 1.0);
	scale   = vector3d(0.2*rgen.rand_uniform(0.6, 1.0), 0.4*rgen.rand_uniform(0.8, 1.0), 1.0); // x = width, y = height, z = length
	float const fzmin(mesh_height + 1.6*radius*scale.y), fzmax(water_plane_z - 2.0*radius*scale.y);
	if (fzmin > fzmax) return 0; // water is too shallow for this size of fish
	pos.z   = min(fzmax, max(fzmin, pos.z)); // clamp to valid depth range
	dir     = vector3d(rgen.signed_rand_float(), rgen.signed_rand_float(), 0.0).get_norm(); // in xy plane
	velocity= zero_vector;
	color   = colorRGBA(rgen.rand_uniform(0.27, 0.33), rgen.rand_uniform(0.20, 0.24), rgen.rand_uniform(0.14, 0.16), 1.0);
	enabled = 1;
	//tile_off.set_from_xyoff2(); // Note: TT specific
	return 1;
}

bool bird_t::gen(rand_gen_t &rgen, cube_t const &range) {

	assert(range.is_strictly_normalized());
	enabled = 0;
	if (atmosphere < 0.5) return 0; // no atmosphere, no clouds, no birds
	pos     = rgen.gen_rand_cube_point(range);
	radius  = 0.1*rgen.rand_uniform(0.6, 1.0);
	dir     = vector3d(rgen.signed_rand_float(), rgen.signed_rand_float(), 0.0).get_norm(); // in xy plane
	velocity= 0.05*dir;
	color   = BLACK;
	enabled = 1;
	//tile_off.set_from_xyoff2(); // Note: TT specific
	return 1;
}


bool fish_t::update(rand_gen_t &rgen) {

	if (!enabled || !animate2) return 0;
	if (pos.z - 1.1*radius*scale.y < interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1)) {enabled = 0; return 0;}
	point const camera(get_camera_pos()), pos_(get_draw_pos());

	if (dist_less_than(pos_, camera, 15.0*radius)) { // scared by the player, swim away
		dir   = (pos_ - camera);
		dir.z = 0.0; // only swims in the xy plane
		dir.normalize();
		velocity = 0.02 * dir; // swim away at a constant velocity
	}
	if (velocity == zero_vector) return 1; // stationary, done

	if (!dist_less_than(pos_, camera, 20.0*radius)) { // far enough away from the player
		velocity *= pow(0.95f, fticks); // slow down
		if (velocity.mag() < 0.001) {velocity = zero_vector;}
	}
	point const orig_pos(pos);
	pos += velocity*fticks;
	if (pos.z - 1.5*radius*scale.y < interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1)) {pos = orig_pos; velocity = zero_vector;} // water is too shallow, stop
	return 1;
}

bool bird_t::update(rand_gen_t &rgen) {

	if (!enabled || !animate2) return 0;

	if (rand() & 127) { // randomly update direction
		float const speed(velocity.mag());
		dir  += rgen.signed_rand_vector(0.1); // 10% change max
		dir.z = 0.0; // only flies in the xy plane
		dir.normalize();
		velocity = dir * speed; // always flies in the direction it's pointed in
	}
	pos += velocity*fticks; // always moving
	return 1;
}


bool animal_t::is_visible(point const &pos_) const {

	if (!enabled) return 0;
	if (!dist_less_than(pos_, get_camera_pos(), 1000.0*radius)) return 0;
	return sphere_in_camera_view(pos_, radius, 0);
}

int animal_t::get_ndiv(point const &pos_) const {
	return min(N_SPHERE_DIV, max(3, int(4.0*sqrt(radius*window_width/distance_to_camera(pos_)))));
}

point animal_t::get_draw_pos() const {
	return (pos + vector3d((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0));
}

void fish_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
	if (!is_visible(pos_)) return;
	colorRGBA draw_color(color);
	water_color_atten_at_pos(draw_color, pos);
	s.set_cur_color(draw_color);
	fgPushMatrix();
	translate_to(pos_);
	rotate_into_plus_z(dir); // FIXME: rotate about z axis (in xy plane)
	scale_by(radius*scale);
	draw_sphere_vbo(all_zeros, 1.0, get_ndiv(pos_), 0);
	fgPopMatrix();
}

void bird_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
	if (!is_visible(pos_)) return;
	float const alpha(CLIP_TO_01(2000.0f*radius/p2p_dist(pos_, get_camera_pos()) - 2.0f)); // 1.0 under half clip distance, after that linear falloff to zero
	s.set_cur_color(colorRGBA(color, alpha)); // FIXME: fog
	fgPushMatrix();
	translate_to(pos_);
	uniform_scale(0.5*radius); // FIXME: adjust radius better in gen()
	draw_sphere_vbo(all_zeros, 1.0, get_ndiv(pos_), 0); // FIXME: placeholder
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
	for (iterator i = begin(); i != end(); ++i) {i->update(rgen);}
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
	// FIXME: VFC using bounding cube?
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

