// 3D World - Dynamic Particle class definition
// by Frank Gennari
// 7/17/06

#include "3DWorld.h"
#include "dynamic_particle.h"
#include "mesh.h"
#include "physics_objects.h"


bool     const DEBUG_TIME     = 0;
bool     const ADD_DP_COBJS   = 0;
unsigned const NUM_COLL_STEPS = 4;
float    const TERMINAL_VEL   = 100.0;
float    const MAX_D_HEIGHT   = 0.1;


dynamic_particle_system d_part_sys;


extern int window_width, iticks, begin_motion, animate2, display_mode, frame_counter;
extern float zbottom, ztop, czmin, czmax, fticks, base_gravity, TIMESTEP, XY_SCENE_SIZE;
extern texture_t textures[];
extern obj_type object_types[];



// ************ dynamic_particle ************


dynamic_particle::dynamic_particle() : sphere_t(all_zeros, rand_uniform(0.03, 0.07)), moves(1), shadows(1), lighted(1),
	collides(1), chdir(0), gravity(0), tid(-1), cid(-1), intensity(rand_uniform(0.6, 1.3)*0.2*XY_SCENE_SIZE),
	bwidth(1.0), velocity(signed_rand_vector(rand_uniform(0.6, 3.0)))
{
	//radius = 24.0*HALF_DXY*intensity*intensity;
	colorRGBA const colors[] = {WHITE, RED, GREEN, BLUE, YELLOW};
	color = colors[rand() % (sizeof(colors)/sizeof(colorRGBA))];
	gen_pos();
}


dynamic_particle::dynamic_particle(point const &p, float r, vector3d const &v, colorRGBA const &c, float i,
		float bw, int t, bool m, bool s, bool l, bool coll, bool cd, bool g)
		: sphere_t(p, r), moves(m), shadows(s), lighted(l), collides(coll), chdir(cd), gravity(g),
		tid(t), cid(-1), intensity(i), bwidth(bw), velocity(v), color(c)
{
	if (pos == all_zeros) gen_pos();
}


void dynamic_particle::gen_pos() {
	
	do {
		rand_xy_point(rand_uniform(zbottom, (MAX_D_HEIGHT + max(ztop, czmax))), pos, 0);
	} while (point_inside_voxel_terrain(pos));
}


void dynamic_particle::draw() const { // lights, color, texture, shadowed, FIXME SHADERS: uses fixed function pipeline

	bool const texture(tid >= 0 && glIsTexture(tid));
	if (texture) select_texture(tid);
	
	if (lighted) { // set emission to color, ambient and diffuse to black (could disable lighting)
		set_color(BLACK);
		set_color_e(color);
	}
	else { // this case is unused
		bool is_shadowed(!is_visible_to_light_cobj(pos, get_light(), radius, -1, 0));
		colorRGBA color_a(color);
		if (color_a != BLACK) {get_indir_light(color_a, (pos + vector3d(0.0, 0.0, 0.01)));}
		set_color_a(color_a);
		set_color_d(is_shadowed ? colorRGBA(0.0, 0.0, 0.0, color.alpha) : color);
	}
	int const ndiv(min(N_SPHERE_DIV, max(3, int(3.0f*sqrt(radius*window_width/distance_to_camera(pos))))));
	bool const bfc(!texture || !textures[tid].has_alpha());
	draw_subdiv_sphere(pos, radius, ndiv, texture, bfc); // point if far away?
	if (lighted) set_color_e(BLACK);
	if (texture) glDisable(GL_TEXTURE_2D);
}


// multiple steps?
void dynamic_particle::apply_physics(float stepsize, int index) { // begin_motion, move, random dir change, collision (mesh and cobjs), forces applied to?

	if (!begin_motion || !animate2) return;

	while (1) {
		if (!is_over_mesh(pos) || pos.z > (MAX_D_HEIGHT + max(ztop, czmax)) || pos.z < zbottom) {
			gen_pos(); // keep within simulation area
			continue;
		}
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
		
		if (point_outside_mesh(xpos, ypos)) { // what about water/ice? stuck in cobj?
			gen_pos();
			continue;
		}
		if (!is_mesh_disabled(xpos, ypos)) {
			float const zval(interpolate_mesh_zval(pos.x, pos.y, radius, 0, 0));

			if ((pos.z - radius) < zval) { // bounce off the surface of the mesh
				pos.z = zval + radius;
				vector3d bounce_v;
				calc_reflection_angle(velocity, bounce_v, surface_normals[ypos][xpos]);
				velocity = bounce_v;
			}
		}
		break;
	}
	if (moves) {
		float const timestep(TIMESTEP*fticks*stepsize);

		if (gravity) {
			float const vz(-min(TERMINAL_VEL, -(velocity.z - base_gravity*GRAVITY*timestep)));
			if (vz < velocity.z) velocity.z = vz;
		}
		if (chdir && (rand() % (100*int(NUM_COLL_STEPS))) < iticks) {
			float const vmag(velocity.mag());
			velocity = signed_rand_vector_norm()*vmag; // same magnitude
		}
		pos += velocity*timestep;
	}
	if (collides) { // hack - check this for correctness
		dwobject obj(DYNAM_PART, pos, velocity, 1, 10000.0); // make a DYNAM_PART object for collision detection
		object_types[DYNAM_PART].radius = radius;
		//obj.multistep_coll(last_pos, index, NUM_COLL_STEPS);
		obj.check_vert_collision(index, 0, 0); // ignoring return value
		pos = obj.pos;
		float const vmag(obj.velocity.mag());
		if (vmag > TOLERANCE) velocity = obj.velocity*(velocity.mag()/vmag); // same magnitude
	}
}


void dynamic_particle::add_light() const { // dynamic lights
	if (lighted) add_dynamic_light(intensity, pos, color, velocity, bwidth); // beam in direction of velocity
}

void dynamic_particle::add_cobj_shadows() const { // cobjs, dynamic objects
	if (shadows) add_shadow_obj(pos, radius, -1);
}

void dynamic_particle::add_mesh_shadows() const {
	if (shadows) dynamic_sphere_shadow(pos, radius, CHECK_ALL_SHADOW, 0); // sphere_shadow? quality?
}

void dynamic_particle::add_cobj() {
	if (ADD_DP_COBJS) cid = add_coll_sphere(pos, radius, cobj_params(0.7, color, 0, 1));
}

void dynamic_particle::remove_cobj() {
	if (ADD_DP_COBJS) remove_coll_object(cid);
	cid = -1;
}


// ************ dynamic_particle_system ************


void dynamic_particle_system::create_particles(unsigned num, bool only_if_empty) {

	if (only_if_empty && size() > 0) return;
	RESET_TIME;
	clear();
	particles.reserve(num);

	for (unsigned i = 0; i < num; ++i) {
		add_particle(dynamic_particle());
	}
	if (DEBUG_TIME && size() > 0) PRINT_TIME("DPS Create");
}


void dynamic_particle_system::draw() const {
	
	RESET_TIME;

	for (unsigned i = 0; i < size(); ++i) {
		particles[i].draw();
	}
	if (DEBUG_TIME && size() > 0) PRINT_TIME("DPS Draw");
}


void dynamic_particle_system::apply_physics(float stepsize) {
	
	RESET_TIME;
	for (unsigned i = 0; i < size(); ++i) {
		particles[i].remove_cobj();

		for (unsigned s = 0; s < NUM_COLL_STEPS; ++s) {
			particles[i].apply_physics(stepsize/NUM_COLL_STEPS, i);
		}
		particles[i].add_cobj();
	}
	if (DEBUG_TIME && size() > 0) PRINT_TIME("DPS Physics");
}


void dynamic_particle_system::add_light() const {
	
	RESET_TIME;
	for (unsigned i = 0; i < size(); ++i) {
		particles[i].add_light();
	}
	if (DEBUG_TIME && size() > 0) PRINT_TIME("DPS Add DLight");
}


void dynamic_particle_system::add_cobj_shadows() const {
	
	RESET_TIME;
	for (unsigned i = 0; i < size(); ++i) {
		particles[i].add_cobj_shadows();
	}
	if (DEBUG_TIME && size() > 0) PRINT_TIME("DPS Cobj Shadows");
}


void dynamic_particle_system::add_mesh_shadows() const {
	
	RESET_TIME;
	for (unsigned i = 0; i < size(); ++i) {
		particles[i].add_mesh_shadows();
	}
	if (DEBUG_TIME && size() > 0) PRINT_TIME("DPS Mesh Shadows");
}


void dynamic_particle_system::build_lookup_matrix() {

	bins.clear();
	bins.resize(XY_MULT_SIZE);

	for (unsigned i = 0; i < size(); ++i) {
		point const &pos(particles[i].get_pos());
		int const xpos(get_xpos(pos.x)), ypos(get_ypos(pos.y));
		if (!point_outside_mesh(xpos, ypos)) bins[xpos + MESH_X_SIZE*ypos].push_back(i);
	}
}

