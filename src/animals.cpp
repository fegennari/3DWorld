// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#include "animals.h"
#include "shaders.h"
#include "gl_ext_arb.h"

float const FISH_RADIUS = 0.05;
float const BIRD_RADIUS = 0.1;
float const FISH_SPEED  = 0.002;
float const BIRD_SPEED  = 0.05;

extern bool water_is_lava;
extern int window_width, animate2;
extern float fticks, water_plane_z, temperature, atmosphere, ocean_wave_height;
extern colorRGBA cur_fog_color;


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

	if ((rand() & 1) == 0) { // randomly update direction
		float const speed(velocity.mag());
		dir += rgen.signed_rand_vector_xy(0.05); // 5% change max
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

void rotate_to_plus_x(vector3d const &dir) {
	rotate_about(TO_DEG*get_norm_angle(dir, plus_x), vector3d(0.0, 0.0, dir.y));
}

void fish_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
	if (!is_visible(pos_, 0.15)) return;
	colorRGBA draw_color(color);
	water_color_atten_at_pos(draw_color, pos_ );
	point const camera(get_camera_pos());
	float const t(CLIP_TO_01((water_plane_z - pos_.z)/max(1.0E-6f, fabs(camera.z - pos_.z))));
	point const p_int(pos_ + (camera - pos_)*t); // ray-water intersection point
	float const dist(p2p_dist(p_int, camera) + 2.0*p2p_dist(p_int, pos_)); // water is 2x as optically dense as air
	float const alpha(min(1.0f, 1.5f*exp(-0.028f*dist/radius)));
	if (alpha < 0.01) return;
	if (alpha < 0.1) {glDepthMask(GL_FALSE);} // disable depth writing to avoid alpha blend order problems
	//draw_color = lerp(cur_fog_color, draw_color, alpha);
	s.set_cur_color(colorRGBA(draw_color, alpha));
	fgPushMatrix();
	translate_to(pos_);
	rotate_to_plus_x(dir);
	scale_by(radius*scale);
	draw_sphere_vbo_back_to_front(all_zeros, 1.0, get_ndiv(pos_), 0);
	if (alpha < 0.1) {glDepthMask(GL_TRUE);}
	fgPopMatrix();
}

void bird_t::draw(shader_t &s) const {

	point const pos_(get_draw_pos());
	if (!is_visible(pos_, 1.0)) return;
	float const dist(p2p_dist(pos_, get_camera_pos()));
	float const alpha(CLIP_TO_01(2000.0f*radius/dist - 2.0f)); // 1.0 under half clip distance, after that linear falloff to zero
	if (alpha < 0.01) return;
	s.set_cur_color(colorRGBA(color, alpha)); // FIXME: fog
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
	bcube = range; // initial value (approx)
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
	if (empty() || bcube.is_zero_area() || !camera_pdu.cube_visible(bcube + get_pos_offset())) return;
	for (const_iterator i = begin(); i != end(); ++i) {i->draw(s);}
}


void vect_fish_t::draw() const {

	if (empty()) return;
	shader_t s;
	begin_draw(s);
	draw_animals(s);
	end_draw(s);
}

void vect_bird_t::draw() const {

	if (empty()) return;
	shader_t s;
	begin_draw(s);
	draw_animals(s);
	end_draw(s);
}


// explicit instantiations
template class animal_group_t<fish_t>;
template class animal_group_t<bird_t>;

