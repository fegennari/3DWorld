// 3D World - Classes for procedural animals
// by Frank Gennari
// 3-17-16

#include "animals.h"
#include "shaders.h"
#include "function_registry.h"
#include "inlines.h"

extern bool water_is_lava;
extern int window_width;
extern float fticks, water_plane_z, temperature, atmosphere;


bool fish_t::gen(rand_gen_t &rgen, cube_t const &range) {

	enabled = 0;
	if (water_is_lava || temperature <= W_FREEZE_POINT || temperature >= WATER_MAX_TEMP) return 0; // too hot/cold for fish
	pos = rgen.gen_rand_cube_point(range);
	float const mesh_height(interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1)), depth(water_plane_z - mesh_height);
	if (depth < 0.1) return 0; // no water
	radius  = 0.1*rgen.rand_uniform(0.4, 1.0);
	scale   = vector3d(rgen.rand_uniform(0.6, 1.0), rgen.rand_uniform(0.8, 1.0), 1.0); // x = width, y = height, z = length
	float const fzmin(mesh_height + 1.6*radius*scale.y), fzmax(water_plane_z - 2.0*radius*scale.y);
	if (fzmin > fzmax) return 0; // water is too shallow for this size of fish
	pos.z   = min(fzmax, max(fzmin, pos.z)); // clamp to valid depth range
	color   = colorRGBA(rgen.rand_uniform(0.27, 0.33), rgen.rand_uniform(0.20, 0.24), rgen.rand_uniform(0.14, 0.16), 1.0);
	enabled = 1;
	return 1;
}

bool bird_t::gen(rand_gen_t &rgen, cube_t const &range) {

	enabled = 0;
	if (atmosphere < 0.5) return 0; // no atmosphere, no clouds, no birds
	pos     = rgen.gen_rand_cube_point(range);
	// FIXME: stuff
	radius  = 0.1*rgen.rand_uniform(0.6, 1.0);
	color   = BLACK;
	enabled = 1;
	return 1;
}


bool fish_t::update(rand_gen_t &rgen) {

	if (!enabled) return 0;
	if (pos.z - 1.1*radius*scale.y < interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1)) {enabled = 0; return 0;}
	point const camera(get_camera_pos());

	if (dist_less_than(pos, camera, 100.0*radius)) { // scared by the player, swim away
		dir   = (pos - camera);
		dir.z = 0.0; // only swims in the xy plane
		dir.normalize();
		velocity = 0.2 * dir; // swim away at a constant velocity
	}
	if (velocity == zero_vector) return 1; // stationary, done

	if (!dist_less_than(pos, camera, 120.0*radius)) { // far enough away from the player
		velocity *= pow(0.95f, fticks); // slow down
		if (velocity.mag() < 0.001) {velocity = zero_vector;}
	}
	point const orig_pos(pos);
	pos += velocity*fticks;
	if (pos.z - 1.5*radius*scale.y < interpolate_mesh_zval(pos.x, pos.y, 0.0, 1, 1)) {pos = orig_pos; velocity = zero_vector;} // water is too shallow, stop
	return 1;
}

bool bird_t::update(rand_gen_t &rgen) {

	if (!enabled) return 0;

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


bool animal_t::is_visible() const {
	return (enabled && sphere_in_camera_view(pos, radius, 0));
}

int animal_t::get_ndiv() const {
	return min(N_SPHERE_DIV, max(3, int(4.0*sqrt(radius*window_width/distance_to_camera(pos)))));
}

void fish_t::draw(shader_t &s) const {

	if (!is_visible()) return;
	s.set_cur_color(color);
	fgPushMatrix();
	translate_to(pos);
	rotate_into_plus_z(dir);
	scale_by(radius*scale);
	draw_sphere_vbo(all_zeros, 1.0, get_ndiv(), 0);
	fgPopMatrix();
}

void bird_t::draw(shader_t &s) const {

	if (!is_visible()) return;
	s.set_cur_color(color);
	fgPushMatrix();
	translate_to(pos);
	uniform_scale(radius);
	draw_sphere_vbo(all_zeros, 1.0, get_ndiv(), 0); // FIXME: placeholder
	fgPopMatrix();
}


void vect_fish_t::draw() const {

	shader_t s;
	s.begin_color_only_shader(); // FIXME: placeholder
	draw_animals(*this, s);
	s.end_shader();
}

void vect_bird_t::draw() const {

	shader_t s;
	s.begin_color_only_shader(); // FIXME: placeholder
	draw_animals(*this, s);
	s.end_shader();
}


