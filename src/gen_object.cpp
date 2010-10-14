// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/10/02

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"


float    const SMOKE_ZVEL      = 3.0;
unsigned const NUM_STARS       = 2500;
unsigned const MAX_BUBBLES     = 2000;
unsigned const MAX_PART_CLOUDS = 200;
unsigned const MAX_FIRES       = 50;


// Global Variables
int nstars(0);
long rseed1(1), rseed2(1);
vector<star> stars(NUM_STARS);
obj_vector_t<bubble> bubbles(MAX_BUBBLES);
obj_vector_t<particle_cloud> part_clouds(MAX_PART_CLOUDS);
obj_vector_t<particle_cloud> cloud_volumes;
obj_vector_t<fire> fires(MAX_FIRES);
float gauss_rand_arr[N_RAND_DIST+2];


extern int star_init, begin_motion, animate2;
extern float zmax_est, zmax, ztop;
extern int coll_id[];
extern point leaf_points[];
extern obj_group obj_groups[];
extern obj_type object_types[];


int valid_star_pos(point &pos);



void gen_stars(int rand_seed, float alpha, int half_sphere) {

	float const cmin[3] = {0.7, 0.9, 0.6};
	int const num_stars(stars.size());

	if (!star_init) {
		if (rand_seed > 1) srand(rand_seed);
		nstars    = num_stars - rand()%max(1, num_stars/4);
		star_init = 1;

		for (int i = 0; i < nstars; ++i) {
			gen_star(stars[i], half_sphere);
		}
	}
	else {
		nstars += rand()%max(1, num_stars/200);
		int const old_num(nstars);
		nstars  = min(nstars, num_stars);

		for (int i = old_num; i < nstars; ++i) {
			gen_star(stars[i], half_sphere);
		}
		if ((rand()&3) == 0) {
			for (int i = 0; i < nstars; ++i) {
				int const rnum(rand());

				if (rnum%21 == 0) {
					stars[i].intensity *= 1.0 + 0.1*signed_rand_float();
					stars[i].intensity = min(stars[i].intensity, 1.0f);
				}
				if (rnum%51 == 0) {
					for (unsigned j = 0; j < 3; ++j) {
						stars[i].color[j] *= 1.0 + 0.05*signed_rand_float();
						stars[i].color[j]  = min(max(stars[i].color[j], cmin[j]), 1.0f);
					}
				}
				if (rnum%101 == 0) {
					point const orig_pos(stars[i].pos);
					vadd_rand(stars[i].pos, 0.002);
					for (unsigned d = 0; d < 3; ++d) stars[i].pos[d] *= SCENE_SIZE[d];
					if (!valid_star_pos(stars[i].pos)) stars[i].pos = orig_pos;
				}
				if (rnum%1001 == 0) gen_star(stars[i], half_sphere);
			}
		}
	}
	draw_stars(alpha);
}


void gen_star(star &star1, int half_sphere) {

	float const radius(0.7*(FAR_CLIP+X_SCENE_SIZE)), theta(rand_uniform(0.0, TWO_PI));
	float phi(gen_rand_phi<rand_uniform>());
	if (!half_sphere && (rand()&1) == 0) phi = PI - phi;
	star1.pos       = rtp_to_xyz(radius, theta, phi);
	star1.intensity = rand_uniform(0.1, 1.0);
	star1.color.assign(rand_uniform(0.7, 1.0), rand_uniform(0.9, 1.0), rand_uniform(0.6, 1.0));
}


inline int valid_star_pos(point &pos) {

	return !(pos.z < 0.0 || (fabs(pos.x) < 20.0*X_SCENE_SIZE && fabs(pos.y) < 20.0*Y_SCENE_SIZE && pos.z < 40.0*Z_SCENE_SIZE));
}


void rand_xy_point(float zval, point &pt, unsigned flags) {

	for (unsigned i = 0; i < 2; ++i) {
		pt[i] = ((flags & FALL_EVERYWHERE) ? 1.5 : 1.0)*SCENE_SIZE[i]*signed_rand_float();
	}
	pt.z = zval;
}


void gen_object_pos(point &position, unsigned flags) {

	rand_xy_point((CLOUD_CEILING + ztop)*(1.0 + rand_uniform(-0.1, 0.1)), position, flags);
}


void basic_physics_obj::init_gen_rand(point const &p2, float rxy, float rz) {

	status = 1;
	time   = 0;
	pos.assign(rand_uniform(-rxy, rxy), rand_uniform(-rxy, rxy), rand_uniform(-rz, rz));
	pos   += p2;
}


void dwobject::print_and_terminate() const { // only called when there is an error

	cout << "pos = "; pos.print();
	cout << ", vel = "; velocity.print();
	cout << ", type = " << int(type) << ", status = " << int(status) << endl;
	assert(0);
}


void bubble::gen(point const &p, float r, colorRGBA const &c) {

	init_gen_rand(p, 0.005, 0.01);
	radius   = ((r == 0.0) ? rand_uniform(0.001, 0.004) : r);
	velocity = rand_uniform(0.15, 0.2);
	pos.z   -= 0.01;
	color    = c;
}


void particle_cloud::gen(point const &p, colorRGBA const &bc, vector3d const &iv, float r,
						 float den, float dark, float dam, int src, bool as, bool use_parts)
{
	init_gen_rand(p, 0.005, 0.025);
	acc_smoke  = as;
	source     = src;
	base_color = bc;
	init_vel   = iv;
	radius     = r;
	darkness   = dark;
	density    = den;
	damage     = dam;

	if (use_parts) {
		parts.resize((rand()&3) + 2); // 2-5 parts

		for (unsigned i = 0; i < parts.size(); ++i) {
			parts[i].pos    = signed_rand_vector(1.0);
			parts[i].pos.z  = fabs(parts[i].pos.z); // always +z
			parts[i].radius = rand_uniform(0.5, 1.0);
			parts[i].status = 1; // always 1
		}
	}
}


void fire::gen(point const &p, float size, int src) {

	init_gen_rand(p, 0.07, 0.01);
	radius   = size*rand_uniform(0.02, 0.05);
	heat     = rand_uniform(0.5, 0.8);
	velocity = zero_vector;
	source   = src;
	cval     = rand_uniform(0.3, 0.7);
	inten    = rand_uniform(0.3, 0.7);
	float const zval(interpolate_mesh_zval(pos.x, pos.y, radius, 0, 0));
	if (fabs(pos.z - zval) > 2.0*radius) status = 2; // above the ground
}


void gen_bubble(point const &pos, float r, colorRGBA const &c) {

	bubbles[bubbles.choose_element()].gen(pos, r, c);
}


void gen_line_of_bubbles(point const &p1, point const &p2, float r, colorRGBA const &c) {

	//RESET_TIME;
	point cur(p1);
	vector3d const dir(p2 - p1);
	float const step_sz(0.02), dist(dir.mag());
	vector3d const step(dir*(step_sz/dist));
	unsigned const nsteps(min(bubbles.size(), unsigned(dist/step_sz)));
	vector<unsigned> ixs;
	bubbles.choose_elements(ixs, nsteps);
	assert(ixs.size() == nsteps);

	for (unsigned i = 0, cur_ix = 0; i < nsteps; ++i) {
		cur += step;
		if (!is_over_mesh(cur)) break;
		int xpos(get_xpos(cur.x)), ypos(get_ypos(cur.y));
		if (point_outside_mesh(xpos, ypos) || cur.z < mesh_height[ypos][xpos]) break;
		if (is_underwater(cur)) bubbles[ixs[cur_ix++]].gen(cur, r, c);
	}
	//PRINT_TIME("Bubble Gen");
}


bool gen_arb_smoke(point const &pos, colorRGBA const &bc, vector3d const &iv,
				   float r, float den, float dark, float dam, int src, bool as)
{
	if (!animate2 || is_underwater(pos)) return 0;
	unsigned const chosen();
	part_clouds[part_clouds.choose_element()].gen(pos, bc, iv, r, den, dark, dam, src, as);
	return 1;
}


void gen_smoke(point const &pos) {

	gen_arb_smoke(pos, WHITE, vector3d(0.0, 0.0, SMOKE_ZVEL),
		rand_uniform(0.01, 0.025), rand_uniform(0.7, 0.9), rand_uniform(0.75, 0.95), 0.0, -2, 1);
}


bool gen_fire(point const &pos, float size, int source) {

	assert(size > 0.0);
	if (!is_over_mesh(pos) || is_underwater(pos)) return 0; // off the mesh or under water/ice
	fires[fires.choose_element()].gen(pos, size, source);
	return 1;
}


void gen_particles(point const &pos, unsigned num, float lt_scale, bool fade) { // lt_scale: 0.0 = full lt, 1.0 = no lt

	obj_group &objg(obj_groups[coll_id[PARTICLE]]);

	for (unsigned o = 0; o < num; ++o) {
		int const i(objg.choose_object());
		objg.create_object_at(i, pos);
		objg.get_obj(i).velocity = gen_rand_vector((fade ? 1.5 : 1.0)*rand_uniform(3.0, 5.0), (fade ? 1.5 : 1.8), 0.75*PI);
		if (lt_scale != 1.0) objg.get_obj(i).time = lt_scale*object_types[PARTICLE].lifetime;
		if (fade) objg.get_obj(i).flags |= TYPE_FLAG;
	}
}


void gen_leaf_at(point const *const points, vector3d const &normal, int type, colorRGB const &color) {

	if (!begin_motion) return;
	int const cid(coll_id[LEAF]);
	if (cid < 0) return;
	obj_group &objg(obj_groups[cid]);
	if (objg.max_objs == 0) return;
	int const max_t_i(objg.choose_object());
	point const pos(get_center(points, 4));
	vector3d const delta(points[1] - points[0]), ldelta(leaf_points[1] - leaf_points[0]);
	float const leaf_size(delta.mag()/ldelta.mag());
	assert(leaf_size > 0.0);
	objg.create_object_at(max_t_i, pos);
	dwobject &obj(objg.get_obj(max_t_i));
	obj.init_dir.z = leaf_size; // sets leaf size
	obj.init_dir.x = 2.0*PI*rand_float(); // angle rotated about z-axis
	//obj.init_dir.x = get_angle(delta.get_norm(), ldelta.get_norm());
	obj.source     = (short)type;

	for (unsigned i = 0; i < 3; ++i) {
		obj.vdeform[i] = color[i]; // stuff color into vdeform
	}
	obj.set_orient_for_coll(&normal);
}


void gen_cloud_volumes() { // 3D cloud puffs

	unsigned const NCLOUDS = 10;
	unsigned const NPARTS  = 1000;
	cloud_volumes.clear();
	srand(123);

	for (unsigned c = 0; c < NCLOUDS; ++c) {
		point const center(4.0*X_SCENE_SIZE*signed_rand_float(),
			               4.0*Y_SCENE_SIZE*signed_rand_float(),
						   (ztop + CLOUD_CEILING + Z_SCENE_SIZE*rand_uniform(0.25, 0.75)));
		point const size(X_SCENE_SIZE*rand_uniform(1.0, 2.0),
			             Y_SCENE_SIZE*rand_uniform(1.0, 2.0),
						 Z_SCENE_SIZE*rand_uniform(0.4, 0.8));
		unsigned const nparts(rand()%(NPARTS/2) + NPARTS/2);
		unsigned const ix(cloud_volumes.size());
		cloud_volumes.resize(ix + nparts);

		for (unsigned p = 0; p < nparts; ++p) {
			point pos(signed_rand_vector_spherical(1.0, 0));

			for (unsigned i = 0; i < 3; ++i) {
				pos[i] *= size[i];
			}
			if (pos.z < 0.0) pos.z *= 0.5; // compressed on the bottom
			pos += center;
			float const radius(0.045*(X_SCENE_SIZE + Y_SCENE_SIZE)*rand_uniform(0.5, 1.0));
			float const density(rand_uniform(0.05, 0.12));
			cloud_volumes[ix + p].gen(pos, WHITE, zero_vector, radius, density, 0.0, 0.0, -((int)c+2), 0, 0);
		}
	}
}


void gen_gauss_rand_arr() {

	float const mconst(2.0E-4*RG_NORM), aconst(((float)N_RAND_GAUSS)*RG_NORM);

	for (int i = 0; i < N_RAND_DIST+2; ++i) {
		float val(0.0);

		for (int j = 0; j < N_RAND_GAUSS; ++j) {
			val += rand()%10000;
		}
		gauss_rand_arr[i] = mconst*val - aconst;
	}
}







