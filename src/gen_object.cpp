// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/10/02

#include "3DWorld.h"
#include "mesh.h"
#include "physics_objects.h"
#include "tree_leaf.h"


float    const SMOKE_ZVEL      = 3.0;
unsigned const NUM_STARS       = 5000;
unsigned const MAX_BUBBLES     = 2500;
unsigned const MAX_PART_CLOUDS = 250;
unsigned const MAX_FIRES       = 75;
unsigned const MAX_DECALS      = 3000;


// Global Variables
vector<star> stars;
obj_vector_t<bubble> bubbles(MAX_BUBBLES);
obj_vector_t<particle_cloud> part_clouds(MAX_PART_CLOUDS);
obj_vector_t<fire> fires(MAX_FIRES);
obj_vector_t<decal_obj> decals(MAX_DECALS);
water_particle_manager water_part_man;
float gauss_rand_arr[N_RAND_DIST+2];
rand_gen_t global_rand_gen;


extern int begin_motion, animate2, show_fog;
extern float zmax, ztop, water_plane_z, FAR_CLIP;
extern int coll_id[];
extern obj_group obj_groups[];
extern obj_type object_types[];



void gen_stars(float alpha, int half_sphere) {

	if (world_mode == WMODE_GROUND && show_fog) return;

	if (stars.empty()) {
		stars.resize(NUM_STARS);
		for (unsigned i = 0; i < stars.size(); ++i) {gen_star(stars[i], half_sphere);}
	}
	if ((rand()&15) == 0) {
		float const cmin[3] = {0.7, 0.9, 0.6};

		for (unsigned i = 0; i < stars.size(); ++i) {
			int const rnum(rand());

			if ((rnum & 15) == 0) { // change intensity
				stars[i].intensity *= 1.0 + 0.1*signed_rand_float();
				stars[i].intensity  = min(stars[i].intensity, 1.0f);
			}
			if ((rnum & 31) == 0) { // change color
				for (unsigned j = 0; j < 3; ++j) {
					stars[i].color[j] *= 1.0 + 0.1*signed_rand_float();
					stars[i].color[j]  = min(max(stars[i].color[j], cmin[j]), 1.0f);
				}
			}
			if (rnum%2000 == 0) {gen_star(stars[i], half_sphere);} // create a new star and destroy the old
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


void rand_xy_point(float zval, point &pt, unsigned flags) {

	for (unsigned i = 0; i < 2; ++i) {
		pt[i] = ((flags & FALL_EVERYWHERE) ? 1.5 : 1.0)*SCENE_SIZE[i]*signed_rand_float();
	}
	pt.z = zval;
}


void gen_object_pos(point &position, unsigned flags) {

	rand_xy_point((CLOUD_CEILING + ztop)*(1.0 + rand_uniform(-0.1, 0.1)), position, flags);
}


void basic_physics_obj::init(point const &p) {

	status = 1;
	time   = 0;
	pos    = p;
}


void basic_physics_obj::init_gen_rand(point const &p, float rxy, float rz) {

	init(p);
	pos += point(rand_uniform(-rxy, rxy), rand_uniform(-rxy, rxy), rand_uniform(-rz, rz));
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
						 float den, float dark, float dam, int src, int dt, bool as, bool use_parts, bool nl)
{
	init_gen_rand(p, 0.005, 0.025);
	acc_smoke  = as;
	source     = src;
	damage_type= dt;
	base_color = bc;
	init_vel   = iv;
	radius     = r;
	init_radius= r;
	darkness   = dark;
	density    = den;
	damage     = dam;
	no_lighting= nl;
	red_only   = 0;

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


void fire::gen(point const &p, float size, float intensity, int src, bool is_static_, float light_bw) {

	if (is_static_) {
		init(p);
		radius = 0.04*size;
		heat   = intensity;
		cval   = 0.5;
		inten  = 0.5;
	}
	else {
		init_gen_rand(p, 0.07, 0.01);
		radius = size*rand_uniform(0.02, 0.05);
		heat   = intensity*rand_uniform(0.5, 0.8);
		cval   = rand_uniform(0.3, 0.7);
		inten  = rand_uniform(0.3, 0.7);
	}
	is_static    = is_static_;
	light_bwidth = light_bw;
	velocity     = zero_vector;
	source       = src;
	float const zval(interpolate_mesh_zval(pos.x, pos.y, radius, 0, 0));
	if (fabs(pos.z - zval) > 2.0*radius) status = 2; // above the ground
}


void decal_obj::gen(point const &p, float r, float ang, vector3d const &o, int lt, int tid_, int cid_, colorRGBA const &color_, bool is_glass_, tex_range_t const &tr) {

	assert(r > 0.0 && color.alpha > 0.0 && lt > 0);
	cid       = cid_; // must be set first
	tid       = tid_;
	lifetime  = lt;
	ipos      = p;
	ipos     -= get_platform_delta(); // make relative to the at-rest platform pos
	init(ipos);
	radius    = r;
	rot_angle = ang;
	color     = color_;
	alpha     = color.alpha;
	orient    = o.get_norm(); // normal of attached surface at collision/anchor point
	pos      += min(0.1*radius, 1.5*DECAL_OFFSET)*orient; // move away from the object it's attached to
	is_glass  = is_glass_;
	tex_range = tr;
}


void gen_bubble(point const &pos, float r, colorRGBA const &c) {

	if (animate2 && begin_motion) {bubbles[bubbles.choose_element()].gen(pos, r, c);}
}


void gen_line_of_bubbles(point const &p1, point const &p2, float r, colorRGBA const &c) {

	//RESET_TIME;
	if (!animate2) return;
	point cur(p1);
	vector3d const dir(p2 - p1);
	float const step_sz(0.02), dist(dir.mag());
	vector3d const step(dir*(step_sz/dist));
	unsigned const nsteps(min((unsigned)bubbles.size(), unsigned(dist/step_sz)));
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
				   float r, float den, float dark, float dam, int src, int dt, bool as)
{
	if (!animate2 || is_underwater(pos) || is_under_mesh(pos)) return 0;
	// Note: we scale by 0.62 since we're using BLUR_CENT_TEX rather than BLUR_TEX to draw smoke (to reduce fill rate)
	part_clouds[part_clouds.choose_element()].gen(pos, bc, iv, 0.62*r, den, dark, dam, src, dt, as);
	return 1;
}


void gen_smoke(point const &pos, float zvel_scale) {

	gen_arb_smoke(pos, WHITE, vector3d(0.0, 0.0, SMOKE_ZVEL*zvel_scale),
		rand_uniform(0.01, 0.025), rand_uniform(0.7, 0.9), rand_uniform(0.75, 0.95), 0.0, NO_SOURCE, SMOKE, 1);
}


bool gen_fire(point const &pos, float size, int source, bool allow_close, bool is_static, float light_bwidth, float intensity) {

	assert(size > 0.0);
	if (!is_over_mesh(pos) || is_underwater(pos)) return 0; // off the mesh or under water/ice

	if (!allow_close) { // check if there are any existing fires near this location and if so, skip this one
		for (unsigned i = 0; i < fires.size(); ++i) {
			if (fires[i].status != 0 && dist_less_than(pos, fires[i].pos, 2.0*fires[i].radius)) {
				fires[i].radius = max(fires[i].radius, size*rand_uniform(0.02, 0.05)); // make it larger
				return 0;
			}
		}
	}
	fires[fires.choose_element()].gen(pos, size, intensity, source, is_static, light_bwidth);
	return 1;
}


void gen_decal(point const &pos, float radius, vector3d const &orient, int tid, int cid, colorRGBA const &color,
	bool is_glass, bool rand_angle, int lifetime, tex_range_t const &tr)
{
	static point last_pos(all_zeros);
	if (dist_less_than(pos, last_pos, 0.5*radius)) return; // skip duplicate/close locations
	last_pos = pos;
	float const rot_angle(rand_angle ? rand_uniform(0.0, TWO_PI) : 0.0);
	decal_obj decal;
	decal.gen(pos, radius, rot_angle, orient, lifetime, tid, cid, color, is_glass, tr);
	if (decal.is_on_cobj(cid)) {decals[decals.choose_element()] = decal;}
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


int gen_fragment(point const &pos, vector3d const &velocity, float size_mult, float time_mult,
	colorRGBA const &color, int tid, float tscale, int source, bool tri_fragment, float hotness)
{
	obj_group &objg(obj_groups[coll_id[FRAGMENT]]);
	int const ix(objg.choose_object());
	objg.create_object_at(ix, pos);
	dwobject &obj(objg.get_obj(ix));
	UNROLL_3X(obj.init_dir[i_] = color[i_];)
	obj.coll_id     = -(tid + 2); // < 0;
	assert(obj.coll_id < 0);
	obj.velocity    = (velocity + gen_rand_vector(rand_uniform(0.3, 0.7), 1.0, PI))*rand_uniform(10.0, 15.0);
	obj.angle       = 360.0*rand_float();
	obj.orientation = signed_rand_vector_norm();
	obj.vdeform.x   = size_mult*(0.5 + rand_float()); // size
	obj.vdeform.y   = color.alpha;
	obj.vdeform.z   = fabs(tscale);
	obj.time        = int(time_mult*object_types[FRAGMENT].lifetime);
	obj.source      = source;
	obj.direction   = (unsigned char)(255.0*CLIP_TO_01(hotness));
	if (tri_fragment) {obj.flags |= TYPE_FLAG;}
	return ix;
}


void gen_leaf_at(point const *const points, vector3d const &normal, int type, colorRGB const &color) {

	if (!begin_motion) return;
	int const cid(coll_id[LEAF]);
	if (cid < 0) return;
	obj_group &objg(obj_groups[cid]);
	if (objg.max_objs == 0) return;
	int const max_t_i(objg.choose_object());
	point const pos(get_center(points, 4));
	float const leaf_size();
	objg.create_object_at(max_t_i, pos);
	dwobject &obj(objg.get_obj(max_t_i));
	obj.init_dir.z = 0.5*p2p_dist(points[0], points[1]); // sets leaf size
	assert(obj.init_dir.z > 0.0);
	obj.init_dir.x = 2.0*PI*rand_float(); // angle rotated about z-axis
	//obj.init_dir.x = get_angle(delta.get_norm(), ldelta.get_norm());
	obj.source     = (short)type;
	for (unsigned i = 0; i < 3; ++i) {obj.vdeform[i] = color[i];} // stuff color into vdeform
	obj.set_orient_for_coll(&normal);
}


void water_particle_manager::gen_particles(point const &pos, vector3d const &vadd, float vmag, float gen_radius, colorRGBA const &color, unsigned num) {

	if (!is_pos_valid(pos)) return; // origin invalid

	for (unsigned i = 0; i < num; ++i) {
		point ppos;
		do {ppos = pos + signed_rand_vector_spherical(gen_radius);} while (!is_pos_valid(ppos)); // find a valid particle starting pos
		vector3d pvel(vadd + signed_rand_vector_spherical(vmag));
		if (pvel.z < 0.0) {pvel.z *= -1.0;} // make sure it's going up
		parts.push_back(part_t(ppos, pvel, color));
	}
}


void add_water_particles(point const &pos, vector3d const &vadd, float vmag, float gen_radius, float mud_mix, float blood_mix, unsigned num) {
	water_part_man.gen_particles(pos, vadd, vmag, gen_radius, water_part_man.calc_color(mud_mix, blood_mix), num);
}


void gen_gauss_rand_arr() {

	float const RG_NORM(sqrt(3.0/N_RAND_GAUSS)), mconst(2.0E-4*RG_NORM), aconst(((float)N_RAND_GAUSS)*RG_NORM);

	for (int i = 0; i < N_RAND_DIST+2; ++i) {
		float val(0.0);
		for (int j = 0; j < N_RAND_GAUSS; ++j) {val += rand()%10000;}
		gauss_rand_arr[i] = mconst*val - aconst;
	}
}


double rgen_core_t::randd() { // FIXME: for some reason, inlining this function changes the values returned
	double rand_num;
	randome_int(rand_num);
	return rand_num/2147483563.;
}

void rgen_pregen_t::pregen_floats(unsigned num) {
	pregen_rand_reals.resize(num);
	for (unsigned i = 0; i < num; ++i) {pregen_rand_reals[i] = rgen_core_t::randd();}
	cur_pos = 0;
}

double rgen_pregen_t::randd() {
	if (pregen_rand_reals.empty()) {return rgen_core_t::randd();}
	float const val(pregen_rand_reals[cur_pos++]);
	if (cur_pos == pregen_rand_reals.size()) {cur_pos = 0;}
	return val;
}

template<typename base> vector3d rand_gen_template_t<base>::signed_rand_vector(float scale) {
	assert(scale > 0.0);
	return vector3d(scale*signed_rand_float(), scale*signed_rand_float(), scale*signed_rand_float());
}

template<typename base> vector3d rand_gen_template_t<base>::signed_rand_vector_norm(float scale) {
	assert(scale > 0.0);

	while (1) {
		vector3d const v(signed_rand_vector(scale));
		if (v.mag_sq() > scale*TOLERANCE) return v.get_norm();
	}
	return zero_vector; // never gets here
}

template<typename base> vector3d rand_gen_template_t<base>::signed_rand_vector_spherical(float scale) {
	assert(scale > 0.0);

	while (1) {
		vector3d const v(signed_rand_vector(scale));
		if (v.mag_sq() < scale*scale) return v;
	}
	return zero_vector; // never gets here
}

template<typename base> point rand_gen_template_t<base>::gen_rand_cube_point(cube_t const &c) {
	point pt;
	UNROLL_3X(pt[i_] = rand_uniform(c.d[i_][0], c.d[i_][1]););
	return pt;
}


// explicit template instantiations
template class rand_gen_template_t<rgen_core_t>;
template class rand_gen_template_t<rgen_pregen_t>;

