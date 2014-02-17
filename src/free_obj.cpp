// 3D World - free_obj and derived classes for universe mode
// by Frank Gennari
// 10/18/05

#include "ship.h"
#include "ship_util.h"
#include "explosion.h"
#include "draw_utils.h"
#include "shaders.h"


bool const PART_COLL_DESTROY = 0;
bool const LG_OBJ_SHADOWS    = 1; // somewhat slow
bool const STENCIL_SHADOWS   = 1; // slightly slow
bool const SELF_SHADOW       = 1; // slow

unsigned const SEEK_CTIME    = unsigned(0.50*TICKS_PER_SECOND);
unsigned const PROJ_ARM_T    = unsigned(0.25*TICKS_PER_SECOND);
unsigned const PART_DELAY    = unsigned(0.40*TICKS_PER_SECOND);
unsigned const ao_s_time     = unsigned(2.00*TICKS_PER_SECOND);
unsigned const ao_e_time     = unsigned(4.00*TICKS_PER_SECOND);

float const GRAVITY_FACTOR   = 1.0E-8; // G?
float const SOLAR_WIND_PRES  = 2.5E-6;
float const COLLISION_ENERGY = 5.0E7;
float const PROJ_ROT_ATTEN   = 0.92;
float const ROT_FORCE_CONST  = 0.5;
float const MAX_ROT_RATE     = 0.03;
float const BLACK_HOLE_GRAV  = 2000.0;
float const ASTEROID_DENSITY = 0.1; // lower density/gravity than planets/moons/stars due to incorrect relative size
float const DEF_AMBIENT_SCALE= 3.2; // hack to increase light on ships (isn't ambient reset below in some cases?)
float const ao_d_time(float(ao_e_time - ao_s_time));

vector3d const init_dir(1.0, 0.0, 0.0);
vector3d const init_up(0.0, 1.0, 0.0); // must be orthogonal to init_dir
vector3d const dir0(1.0, 0.0, 0.0);
vector3d const upv0(0.0, 1.0, 0.0);


extern bool univ_stencil_shadows;
extern int iticks, begin_motion, display_mode;
extern float fticks;
extern point player_death_pos;
extern pos_dir_up player_pdu;
extern vector<us_weapon> us_weapons;
extern exp_type_params et_params[];
extern pt_line_drawer particle_pld, glow_pld;
extern pt_line_drawer_no_lighting_t emissive_pld;


// ************ FREE_OBJ ************


void free_obj::reset() {

	status      = 0; // no longer destroyed
	time        = 0;
	reset_timer = 0;
	shadow_val  = 0;
	sobj_coll_tid = -1;
	temperature = 0.0;
	extra_mass  = 0.0;
	rot_rate    = 0.0;
	sobj_dist   = 0.0;
	ambient_scale = DEF_AMBIENT_SCALE;
	draw_rscale = 1.0;
	target_obj  = NULL;
	parent      = NULL;
	pos         = reset_pos;
	velocity    = zero_vector;
	dvel        = zero_vector;
	upv         = init_up;
	dir         = init_dir.get_norm(); // must be orthogonal to up_vector
	rot_axis    = plus_z;
	gvect       = zero_vector;
	near_b_hole = 0;
	ra1 = ra2   = 0.0;
	invalidate_rotv();
	reset_lights(); // should already be reset
}


free_obj const *free_obj::get_root_parent() const {
	
	free_obj const *p1(this);
	
	while (p1->parent) {
		if (VERIFY_REFS || 1) p1->verify_status(); // this one is important
		p1 = p1->parent;
	}
	return p1;
}


void free_obj::check_ref_objs() {
	
	if (status != 0) return;

	while (target_obj != NULL && target_obj->not_a_target()) { // try to find a valid target from our current target's parent/ancestors
		if (target_obj->parent != NULL && target_obj->invalid() && target_valid(target_obj->parent) &&
			target_obj->parent->get_align() != alignment) // almost correct
		{
			target_obj = target_obj->parent; // now go after the ship's parent
		}
		else {
			target_obj = NULL;
		}
	}
	while (parent != NULL && parent->invalid()) { // find a valid parent/ancestor to be our new parent
		if (is_ship()) {
			parent = parent->parent; // take orders from parent's parent now, can't add to new parent's fighter list
		}
		else {
			parent = NULL;
		}
	}
	if (VERIFY_REFS) {
		assert(!target_obj || (!target_obj->invalid() && !target_obj->to_be_removed()));
		assert(!parent     || (!parent->invalid()     && !parent->to_be_removed()));
		if (target_obj != NULL) target_obj->verify_status();
		if (parent     != NULL) parent->verify_status();
	}
}


void free_obj::fix_upv() {

	vector3d const cp(cross_product(upv, dir));
	upv = dir;
	rotate_vector3d_norm(cp, PI_TWO, upv); // force upv to be orthogonal to dir
	invalidate_rotv();
}


void free_obj::accelerate(float speed, float accel) {

	if (speed == 0.0 || accel == 0.0) return;
	float const vmag0(velocity.mag());
	velocity += dir*(accel*speed); // small change in velocity
	velocity.set_max_mag(max(vmag0, float(fabs(speed)))); // normalize to speed = max speed
	assert(!is_nan(velocity));
}
	

void free_obj::decelerate(float speed, float decel) {

	if (speed == 0.0 || decel == 0.0) return;
	float const vmag(velocity.mag()), vmax(max(vmag, float(fabs(speed))));

	if (vmag < decel*vmax) {
		velocity = zero_vector; // stopped
	}
	else if (vmag > TOLERANCE) {
		velocity *= ((vmag - decel*vmax)/vmag); // slow down
		assert(!is_nan(velocity));
	}
}


void free_obj::arb_rotate_about(float rangle, vector3d const &raxis) {

	rotate_vector3d_x2(raxis, rangle, upv, dir);
	upv.normalize();
	dir.normalize();
	invalidate_rotv();
}


void free_obj::do_rotate(float pitch, float yaw) { // upv rotates as you move the mouse in xy circles

	vector3d const vy(cross_product(dir, upv));
	rotate_vector3d(upv, yaw, dir); // x-rotation
	arb_rotate_about(pitch, vy); // y-rotation, rotate up vector
}


void free_obj::tilt(float val) {
	
	rotate_vector3d_norm(dir, val, upv);
	invalidate_rotv();
}


void free_obj::add_gravity_swp(vector3d const &gravity, vector3d const &swp, float gscale, bool near_bh) {

	near_b_hole = near_bh;
	gvect       = gravity;
	velocity   += swp*SOLAR_WIND_PRES*gscale; // apply solar wind pressure (usually negligible)
	if (gravity == zero_vector) return;
	velocity   += gravity*GRAVITY_FACTOR*gscale; // velocity scale due to gravitational pull
	assert(!is_nan(velocity));
	
	if (near_b_hole && gravity.mag() > 0.95*BLACK_HOLE_GRAV) { // black hole does collision damage
		damage(1000.0, DAMAGE_COLL, pos, this, WCLASS_COLLISION);
	}
}


void free_obj::apply_torque_force(point const &fpos, vector3d const &fdir, float fmag, float mass) {

	assert(mass > 0.0);
	if (fpos == pos) return;
	float const dist(p2p_dist(fpos, pos)), fdir_mag(fdir.mag()), fm(mass*fdir_mag);
	if (fm < TOLERANCE) return;
	vector3d new_rot_axis(cross_product(fdir, vector3d(fpos, pos)));
	new_rot_axis *= ROT_FORCE_CONST*fmag/fm; // apply force (independent of radius)
	if (dist > c_radius) new_rot_axis *= c_radius/dist; // difficult to do since we don't know the real point of intersection
	rot_axis     *= rot_rate;       // apply magnitude to old rotation axis
	rot_axis     += new_rot_axis;   // sum the forces (torques)
	rot_rate      = rot_axis.mag(); // new rotation speed

	if (rot_rate < TOLERANCE) {
		rot_axis = plus_z;
		rot_rate = 0.0;
	}
	else {
		rot_axis /= rot_rate;
		rot_rate  = min(MAX_ROT_RATE, rot_rate);
	}
}


// could be in math3d.cpp, but modifies pos and indirectly velocity and requires a ton of parameters
float free_obj::coll_physics(point const &copos, vector3d const &vcoll, float obj_mass, float obj_radius,
							 free_obj const *const source, float elastic, float *dscale)
{
	// elastic collision: m1 = mass0, m2 = obj_mass, v1 = velocity, v2 = vcoll
	// conservation of momentum: m1*v1   + m2*v2   = m1*v1'   + m2*v2'
	// conservation of energy  : m1*v1^2 + m2*v2^2 = m1*v1'^2 + m2*v2'^2 (1/2 factor has been removed)
	// solve                   : v1' = v1*(m1 - m2)/(m1 + m2) + v2*2*m2/(m1 + m2)
	// 3D:
	//  c = n.(v1i - v2i)
	//  v1f = v1i - ((m2c)/(m1 + m2))*(1 + e)*n
	float const mass0(get_mass());
	assert(radius > 0.0 && obj_mass > 0.0 && mass0 > 0.0);
	vector3d const vtot(get_tot_vel_at(copos)), vd(vtot, vcoll);
	vector3d norm(pos, copos); // normal from collided object to current object
	bool const stopped(vd.mag() < TOLERANCE), move_pos((stopped || obj_mass >= 0.5*mass0) && obj_radius > 0.0 && c_radius > 0.0);
	point new_pos(copos);

	if (norm.mag_sq() < TOLERANCE*TOLERANCE) { // points on top of each other, just pick a direction
		norm = signed_rand_vector_norm();
	}
	else {
		norm.normalize();
	}
	if (stopped && !move_pos) { // non-movable objects starting on top of each other - what to do?
		dvel += norm*1.0E-6; // give it a small velocity
		return 0.0;
	}
	if (move_pos) { // start with the intersection point on the edge of the source collision object's bounding sphere
		float const rval(obj_radius + c_radius + 0.01*min(obj_radius, c_radius)); // slightly larger than rsum
		//new_pos += norm*rval; // move the smaller object by its radius from the collision point so as not to intersect
		get_sphere_mov_sphere_int_pt(copos, pos, vd, rval, new_pos); // velocity with or without rotations?
			
		if (obj_mass >= 2.0*mass0) { // does this help?
			float const dist(p2p_dist(copos, new_pos));
	
			if (dist < rval) {
				norm    = (new_pos - copos)/max(dist, TOLERANCE);
				new_pos = copos + norm*rval;
			}
		}
	}
	if (source != NULL) {
		// Technically either 1.0/NUM_TIMESTEPS or fticks/NUM_TIMESTEPS is correct, but extending a little further shoudn't hurt
		// and this way we can make it more robust in cases where fticks is highly variable and/or far from 1.0.
		// Also, should really take into account linear velocity at the collision point
		// due to angular velocity of the object (from turning, rolling, and spinning from torque)
		float const timestep(1.0 + fticks);
		point const last_pos(pos - velocity*timestep); // close but not entirely correct (delta velocity?)
		intersect_params ip(last_pos, new_pos, norm);

		if (source->has_detailed_coll(this) && source->obj_int_obj(this, ip)) { // will update new_pos and norm if necessary
			new_pos = ip.p_int;
			norm    = ip.norm;
			if (dscale) *dscale = ip.dscale;
		}
		if (move_pos && !is_particle()) {
			vector3d const dpos(new_pos - pos);
			float const move_dist(dpos.mag()), frame_dist(max(0.05f*c_radius, p2p_dist(last_pos, pos)));

			if (move_dist > max(TOLERANCE, frame_dist)) { // moved too far - clamp max dist to frame_dist
				new_pos = pos + dpos*(frame_dist/move_dist); // doesn't work if one object starts inside of another
			}
		}
	}
	if (move_pos) { // technically should move both objects, but more difficult to ensure no collision afterwards
		assert(!is_nan(new_pos));
		pos = new_pos; // technically should move both objects, but more difficult to ensure no collision afterwards
	}
	if (stopped) return 0.0;
	
	// compute collision force
	float const c(-dot_product(norm, vd)); // norm is not really correct in all cases
	//if (c <= 0.0) return 0.0; // ???
	vector3d const deltav(norm*(c*(obj_mass/(mass0 + obj_mass))*(1.0 + elastic)));
	assert(!is_nan(deltav));
	dvel += deltav; // accumulate velocity change and add into velocity after all collisions are processed
	return COLLISION_ENERGY*get_coll_energy(velocity, (velocity + deltav), mass0); // too bad we don't know the other half of the energy
}


void free_obj::check_distant() {
	
	if (is_distant(pos, get_draw_radius())) {flags |= OBJ_FLAGS_DIST;} else {flags &= ~OBJ_FLAGS_DIST;}
}


void free_obj::apply_physics() {

	if (reset_timer > 0) {
		unsigned const ticks(iticks); //max(1, iticks)

		if (reset_timer <= ticks) { // dead
			bool const player(is_player_ship());
			if (player) player_death_pos = get_player_pos2();
			reset();
			if (player) reset_player_universe(); // hack
			reset_timer = 0;
		}
		else {
			reset_timer -= ticks;
		}
	}
	time += iticks;
	check_distant();
}


void free_obj::advance_time(float timestep) {

	if (dvel != zero_vector) {
		velocity += dvel; // collision velocity
		dvel      = zero_vector;
		assert(!is_nan(velocity));
	}
	UNROLL_3X(if (fabs(velocity[i_]) < TOLERANCE) {velocity[i_] = 0.0;}) // fix for denormalized numbers?
	//if (is_ship() && !is_player_ship()) {velocity = zero_vector;} // testing
	pos += velocity*timestep;
	assert(!is_nan(pos));

	if (rot_rate != 0.0) { // apply rotation
		if (fabs(rot_rate) < 1.0E-6) {
			rot_rate = 0.0;
		}
		else {
			assert(rot_rate > 0.0);
			arb_rotate_about(fticks*rot_rate, rot_axis);
		}
	}
	sobj_coll_tid = -1;
}


void free_obj::set_max_speed(float max_speed) { // Note: if max speed is 0 then won't move, no effect of collision or gravity

	if (max_speed < TOLERANCE) {
		velocity = zero_vector;
		return;
	}
	if (velocity.mag_sq() > max_speed*max_speed) {
		velocity.normalize();
		velocity *= max_speed;
		assert(!is_nan(velocity));
	}
}


void free_obj::set_temp(float temp, point const &tcenter, free_obj const *source) {
	
	temperature = temp;

	if (is_burning()) {
		float const otf(get_over_temp_factor());
		
		if (otf > 0.0) {
			if (source == NULL) source = this;
			float const dscale(40.0*max(0.0001f, radius)); // larger objects take more temperature damage since they have more surface area
			damage(dscale*otf, DAMAGE_HEAT, pos, source, WCLASS_HEAT); // burn over time
		}
	}
}


void free_obj::add_light(unsigned index) {

	if (num_exp_lights < NUM_EXP_LIGHTS) {
		exp_lights[num_exp_lights++] = index;
		return;
	}
	assert(num_exp_lights <= NUM_EXP_LIGHTS);
	unsigned old_i(0);

	for (unsigned i = 1; i < num_exp_lights; ++i) { // find lowest priority light
		if (higher_priority(exp_lights[old_i], exp_lights[i])) old_i = i;
	}
	if (higher_priority(index, exp_lights[old_i])) {
		exp_lights[old_i] = index; // index is higher priority than oldest light
	}
}


vector3d free_obj::get_orient() const {

	switch (auto_orient()) {
	case 0: return dir; // should already be normalized
	case 1: return velocity.get_norm();
	case 2:
		{
			if (time <= ao_s_time) return dir;
			float const vmag(velocity.mag()); // should already be normalized
			if (vmag < TOLERANCE ) return dir; // stopped
			if (time >= ao_e_time) return velocity/vmag;
			float const val((float(time - ao_s_time))/ao_d_time);
			return velocity*(val/vmag) + dir*(1.0 - val); // should be normalized
		}
	default: assert(0);
	}
	return zero_vector; // shouldn't get here
}


void free_obj::calc_rotation_vectors() const {

	if (rv1 != zero_vector) return; // already cached

	if (!is_ship() && !calc_rvs()) {
		ra1 = ra2 = 0.0;
		rv1 = rv2 = dir;
		return;
	}
	force_calc_rotation_vectors();
}


void free_obj::force_calc_rotation_vectors() const {

	// rotate in direction <dir>
	vector3d const orient(get_orient()); // normalized
	cross_product(upv0, dir0, rv2); // constant, could be cached?
	cross_product(rv2, orient, rv1);
	if (rv1.mag() < TOLERANCE) {rv1 = plus_z;} // otherwise this will assertion fail
	ra1 = get_angle(rv2.get_norm(), orient);

	// rotate up to <upv>
	vector3d upvr;
	rotate_vector3d(upv, rv1, ra1, upvr);
	upvr.normalize();
	ra2 = get_angle(upvr, dir0);
	if (dot_product(upvr, upv0) > 0.0) ra2 = -ra2;
}


template<typename T> void free_obj::rotate_point(pointT<T> &pt) const { // rotate from global coordinate system into object coordinate system

	calc_rotation_vectors();
	rotate_vector3d(pointT<T>(rv1), ra1, pt); // inverse rotate in direction <dir>
	rotate_vector3d(pointT<T>(rv2), ra2, pt); // inverse rotate up to <upv>
}


template<typename T> void free_obj::xform_point(pointT<T> &pt) const { // transform from global coordinate system into object coordinate system

	pt -= pos; // translate to (0,0,0)
	rotate_point(pt);
	pt /= radius; // scale by inverse of radius
}


template<typename T> void free_obj::rotate_point_inv(pointT<T> &pt) const { // rotate from object coordinate system into global coordinate system

	calc_rotation_vectors();
	rotate_vector3d(pointT<T>(rv2), -ra2, pt); // rotate up to <upv>
	rotate_vector3d(pointT<T>(rv1), -ra1, pt); // rotate in direction <dir>
}


template<typename T> void free_obj::xform_point_inv(pointT<T> &pt) const { // transform from object coordinate system into global coordinate system

	rotate_point_inv(pt);
	pt *= T(radius); // scale by radius
	pt += pos; // translate to pos
}


void free_obj::xform_point_x2(point &p1, point &p2) const { // transform from global coordinate system into object coordinate system

	p1 -= pos; // translate to (0,0,0)
	p2 -= pos; // translate to (0,0,0)
	calc_rotation_vectors();
	rotate_vector3d_x2(rv1, ra1, p1, p2); // inverse rotate in direction <dir>
	rotate_vector3d_x2(rv2, ra2, p1, p2); // inverse rotate up to <upv>
	p1 /= radius; // scale by inverse of radius
	p2 /= radius; // scale by inverse of radius
}


void free_obj::xform_point_inv_multi(upos_point_type *pts, unsigned npts) const {

	calc_rotation_vectors();
	rotate_vector3d_multi(upos_point_type(rv2), -ra2, pts, npts); // rotate up to <upv>
	rotate_vector3d_multi(upos_point_type(rv1), -ra1, pts, npts); // rotate in direction <dir>

	for (unsigned i = 0; i < npts; ++i) {
		pts[i] *= radius; // scale by radius
		pts[i] += pos; // translate to pos
	}
}

template void free_obj::rotate_point    (point   &pt) const;
template void free_obj::rotate_point    (point_d &pt) const;
template void free_obj::xform_point     (point   &pt) const;
template void free_obj::xform_point     (point_d &pt) const;
template void free_obj::rotate_point_inv(point   &pt) const;
template void free_obj::rotate_point_inv(point_d &pt) const;
template void free_obj::xform_point_inv (point   &pt) const;
template void free_obj::xform_point_inv (point_d &pt) const;


vector3d free_obj::calc_angular_vel(point const &cpos, vector3d const &axis, float rate) const {

	if (rate == 0.0) return zero_vector;
	float const dist_sq(p2p_dist_sq(cpos, pos));
	if (dist_sq > c_radius*c_radius) rate *= c_radius/sqrt(dist_sq); // normalize to max dist
	return cross_product((cpos - pos), axis)*TWO_PI*rate;
}


vector3d free_obj::get_tot_vel_at(point const &cpos) const {
	
	return (velocity + calc_angular_vel(cpos, rot_axis, rot_rate));
}


float free_obj::damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {
	
	if (val <= 0.0) return 0.0; // no damage
	def_explode((10.0 + 0.1*min(100.0f, val)), ETYPE_ANIM_FIRE, dir, WCLASS_EXPLODE, alignment, 0, parent); // no shields or armor
	return 0.0;
}


// get_all currently only works with free_objs, not uobjects
void free_obj::large_obj_visibility_test(point const &lpos, float exp_radius, vector<uobject const*> &sobjs, bool get_all) const {

	float const dist(p2p_dist(lpos, pos));
	if (dist < 0.1*radius) return;
	free_obj *fobj;
	line_int_data li_data(pos, (lpos - pos), min(50.0f*radius, dist), this, NULL, 0, 0, exp_radius);
	li_data.visible_only = 1;
	li_data.use_lpos     = 1;
	li_data.lpos         = lpos;
	if (get_all) li_data.sobjs = &sobjs;
	uobject const *const obj(line_intersect_objects(li_data, fobj, OBJ_TYPE_LARGE));
	if (!get_all && obj != NULL) sobjs.push_back(obj);
}


void free_obj::rotate() const {

	rotate_about(TO_DEG*ra1, rv1); // rotate in direction <dir>
	rotate_about(TO_DEG*ra2, rv2); // rotate up to <upv>
}


void free_obj::inverse_rotate() const {

	rotate_about(-TO_DEG*ra2, rv2);
	rotate_about(-TO_DEG*ra1, rv1);
}


void free_obj::transform_and_draw_obj(uobj_draw_data &udd, bool specular, bool first_pass, bool final_pass) const {

	glPushMatrix();
	global_translate(pos);
	uniform_scale(radius);

	if (near_b_hole && gvect.mag() > 0.05*BLACK_HOLE_GRAV) { // stretch the object
		vector3d gscale;
		UNROLL_3X(gscale[i_] = fabs(gvect[i_]) + 0.1*BLACK_HOLE_GRAV;)
		gscale.normalize();
		float volume(1.0);
		UNROLL_3X(volume *= gscale[i_];)
		gscale *= pow(1.0/volume, 1.0/3.0);
		scale_by(gscale);
	}
	if (!is_particle()) glPushMatrix();
	if (!udd.draw_as_pt() && calc_rvs()) rotate();
	udd.specular_en = specular;
	udd.first_pass  = first_pass;
	udd.final_pass  = final_pass;
	draw_obj(udd);
	if (!is_particle()) glPopMatrix();
	glPopMatrix();
}


void free_obj::draw(shader_t shader[2]) const { // view culling has already been performed

	//RESET_TIME;
	if (!is_ok()) return; // dead
	bool const lg_obj_type(is_ship() || is_stationary());
	float const dist(p2p_dist(pos, get_player_pos2())), type_scale(lg_obj_type ? 0.75 : (is_particle() ? 1.2 : 1.0));
	float const dscale(type_scale*NDIV_SCALE_U*(get_draw_radius()/(dist + 0.1*radius + TOLERANCE)));
	if (dscale < 0.5 || (is_particle() && dscale < 1.0)) return; // too far/small - clip it
	int ndiv(max(3, min((int)FREE_OBJ_MAX_NDIV, int(sqrt(10.0*dscale)))));
	if (ndiv > 8 && (ndiv&1)) ++ndiv; // an even size is divisible by 2, so hemispheres can be created exactly
	bool const stencil_shadows(STENCIL_SHADOWS && univ_stencil_shadows && lg_obj_type), no_lighting(no_light());
	uobject const *sobj(NULL);
	static vector<uobject const *> sobjs;
	sobjs.resize(0);
	point sun_pos; // should be the same as lpos?
	int shadowed(0);
	bool partial_shadow(0);
	
	if (LG_OBJ_SHADOWS && !no_lighting && lg_obj_type && ndiv > 5 && get_universe_sun_pos(pos, sun_pos)) {
		large_obj_visibility_test(sun_pos, c_radius, sobjs, 1); // radius or c_radius? seems to get decoy flares as well as ships
		
		for (unsigned i = 0; i < sobjs.size(); ++i) {
			assert(sobjs[i]);
			shadowed = max(shadowed, (1 + (sobjs[i]->get_radius() > 1.5*radius)));
		}
		if (SELF_SHADOW && ndiv >= 16 && self_shadow()) {
			cobj_vector_t const &cobjs(get_cobjs());

			if (cobjs.size() > 1 || (cobjs.size() == 1 && cobjs[0]->is_concave())) { // maybe don't have to check this?
				sobjs.push_back(this); // self shadow
				shadowed = max(shadowed, 1);
			}
		}
	}
	bool const known_shadowed(shadowed == 2 || (stencil_shadows && shadowed));
	int const shadow_thresh(stencil_shadows ? 0 : 1);
	int const light_val(no_lighting ? 0 : set_uobj_color(pos, c_radius, known_shadowed, shadow_thresh, sun_pos, sobj, ambient_scale, ambient_scale));
	shadow_val = max(shadowed, light_val); // only updated if drawn - close enough?
	if (is_player_ship()) return; // don't draw player ship
	if (light_val > 0 && sobj != NULL) sobjs.push_back(sobj);
	assert(num_exp_lights <= NUM_EXP_LIGHTS);
	unsigned const nlights(min((unsigned)dscale, num_exp_lights)); // hack to avoid adding too many lights to tiny objects

	if (stencil_shadows && light_val != 3) {
		if (shadowed) {
			assert(light_val >= 0); // sun_pos is set
			partial_shadow = 1; // shadow model?
		}
		else if (light_val == 1 || light_val == 2) {
			partial_shadow = 1;
		}
	}
	set_fill_mode();
	calc_rotation_vectors();
	unsigned const npasses(partial_shadow ? get_num_draw_passes() : 1); // will be slow if > 1
	bool const specular(!known_shadowed && (light_val == 0 || (!stencil_shadows && light_val == 1))); // less than half shadowed
	uobj_draw_data udd(this, &shader[0], ndiv, time, powered(), specular, 0, pos, velocity, dir, upv,
		dist, radius, c_radius/radius, (nlights > 0), 1, !partial_shadow, 1, (npasses == 1));

	if (ndiv > 3) {
		for (unsigned i = 0; i < nlights; ++i) {
			setup_br_light(exp_lights[i], pos, (EXPLOSION_LIGHT + i));
		}
	}
	for (unsigned pass = 0; pass < npasses; ++pass) {
		if (pass > 0) {
			set_uobj_color(pos, c_radius, known_shadowed, shadow_thresh, sun_pos, sobj, ambient_scale, ambient_scale);
			udd.phase1 = 0;
			udd.phase2 = 1;
		}
		transform_and_draw_obj(udd, specular, 1, !partial_shadow);

		if (partial_shadow) { // partially shadowed - draw the sun's light with a stencil pass
			// http://www.gamasutra.com/features/20021011/lengyel_05.htm
			if (shader[1].is_setup()) {shader[1].enable(); udd.shader = &shader[1];}
			assert(!sobjs.empty());
			glPushMatrix();
			global_translate(pos);
			glClear(GL_STENCIL_BUFFER_BIT);
			glEnable(GL_STENCIL_TEST);
			glStencilFunc(GL_ALWAYS, 0, ~0);
			glEnable(GL_DEPTH_TEST);
			glDepthFunc(GL_LESS);
			glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
			glDepthMask(GL_FALSE);
			glEnable(GL_CULL_FACE);

			for (unsigned d = 0; d < sobjs.size(); ++d) { // *** planet/moon exact shadow? ***
				draw_shadow_volumes_from(sobjs[d], sun_pos, dscale, ndiv, 0);
			}
			glPopMatrix();
			glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
			glDepthFunc(GL_EQUAL);
			glStencilFunc(GL_EQUAL, 0, ~0);
			glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
			set_additive_blend_mode(); //glBlendFunc(GL_ONE, GL_ONE);
			glCullFace(GL_BACK);
			glDisable(GL_CULL_FACE);

			set_uobj_color(pos, c_radius, 0, 2, sun_pos, sobj, 0.0, 0.0); // enable diffuse/specular only
			transform_and_draw_obj(udd, 1, 0, 1);

			glDepthMask(GL_TRUE);
			glDepthFunc(GL_LEQUAL);
			glDisable(GL_STENCIL_TEST);
			set_std_blend_mode();
			if (shader[0].is_setup()) {shader[0].enable(); udd.shader = &shader[0];}

			if (display_mode & 0x10) { // testing
				set_emissive_color(GREEN); // will be reset
				glPushMatrix();
				global_translate(pos);

				for (unsigned d = 0; d < sobjs.size(); ++d) {
					if (sobjs[d] != &player_ship()) draw_shadow_volumes_from(sobjs[d], sun_pos, dscale, ndiv, 1);
				}
				glPopMatrix();
				clear_emissive_color();
			}
			glDepthFunc(GL_LESS);
		} // partial_shadow
	} // pass
	//if (GET_DELTA_TIME > 10) cout << get_name() << ": " << GET_DELTA_TIME << endl;
}


// ************ STATIONARY_OBJ ************


stationary_obj::stationary_obj(unsigned type_, point const &pos_, float radius_, unsigned lt) : type(type_), lifetime(lt) {

	assert(type < NUM_SO_TYPES);
	flags    = OBJ_FLAGS_STAT; // so that targeting works
	pos      = pos_;
	radius   = radius_;
	c_radius = radius;
	dir      = signed_rand_vector_norm();
	
	switch (type) {
		case SO_BLACK_HOLE:
			flags |= (OBJ_FLAGS_NOLT | OBJ_FLAGS_NCOL);
			break;
		case SO_ASTEROID: // do nothing
			break;
		default:
			assert(0);
	}
}


float stationary_obj::damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {
	
	if (type == DAMAGE_DESTROY) def_explode(2.5, ETYPE_ANIM_FIRE, dir, WCLASS_EXPLODE, alignment, 0, NULL); // can't be damaged by anything else
	return 0.0;
}


void stationary_obj::apply_physics() {

	time += iticks;
	if (lifetime > 0 && time > lifetime) status = 1; // lifetime expired, destroy it
}


// similar to uobj_solid::add_gravity_vector()
int stationary_obj::get_gravity(vector3d &vgravity, point const &mpos) const {

	switch (type) {
		case SO_BLACK_HOLE:
			add_gravity_vector_base(vgravity, mpos, 0.1*BLACK_HOLE_GRAV*radius, BLACK_HOLE_GRAV);
			return 2;
		case SO_ASTEROID: // not spherical
			add_gravity_vector_base(vgravity, mpos, ASTEROID_DENSITY*radius, MAX_SOBJ_GRAVITY); // gravity = density*radius
			return 1;
		default: assert(0);
	}
	return 0;
}


void stationary_obj::draw_obj(uobj_draw_data &ddata) const {

	switch (type) {
		case SO_BLACK_HOLE:
			ddata.draw_black_hole();
			break;
		default:
			assert(0);
	}
}


// ************ UPARTICLE ************


void uparticle::set_params(unsigned ptype_, point const &pos_, vector3d const &vel, vector3d const &d, float radius_,
						   colorRGBA const &c1, colorRGBA const &c2, unsigned lt, float damage_, unsigned align, bool coll_, int tid) {
	
	flags = OBJ_FLAGS_PART;
	
	if (!coll_) {
		add_flag(OBJ_FLAGS_NCOL);
		add_flag(OBJ_FLAGS_NEXD);
	}
	ptype     = ptype_;
	color1    = c1;
	color2    = c2;
	lifetime  = lt;
	pos       = pos_;
	velocity  = vel;
	dir       = d;
	radius    = radius_;
	c_radius  = radius;
	damage_v  = damage_;
	alignment = align;
	texture_id= tid;
	angle     = 0.0;
	rrate     = 2.0*rand_float();
	axis      = signed_rand_vector_norm();
	draw_rscale = ((ptype == PTYPE_GLOW) ? 1.5 : 1.0);
	invalidate_rotv();
}


bool uparticle::dec_ref() {
	
	if (alloc_block == NULL) return 0;
	alloc_block->free_obj();
	return 1;
}


vector3d uparticle::get_tot_vel_at(point const &cpos) const {

	return free_obj::get_tot_vel_at(cpos) + calc_angular_vel(cpos, axis, 0.1*rrate);
}


void uparticle::apply_physics() {

	if (time > lifetime) {status = 1; return;}
	free_obj::apply_physics();
	angle += rrate*fticks; // increase rotation
}


bool uparticle::collision(point const &copos, vector3d const &vcoll, float obj_mass, float obj_radius, free_obj *source, float elastic) {

	if (damage_v > 1.0E-3 && source != NULL) {
		source->damage(damage_v, DAMAGE_HEAT, pos, this, WCLASS_HEAT); // damage the target (hot particle)
		damage_v *= 0.95; // particle loses energy
	}
	if (PART_COLL_DESTROY) {
		damage(1.0, DAMAGE_COLL, copos, source, SWCLASS_UNDEF);
	}
	else {
		coll_physics(copos, vcoll, obj_mass, obj_radius, source, elastic);
	}
	return 1;
}


float uparticle::damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {
	
	if (val <= 0.0) return 0.0; // no damage
	if (time > PART_DELAY && (ptype != PTYPE_GLOW || val > 1000.0)) status = 1; // destroy if damaged, no explosion, it just goes away
	return 0.0;
}


void uparticle::draw_obj(uobj_draw_data &ddata) const {

	colorRGBA color(color1);
	if (color1 != color2) {blend_color(color, color1, color2, CLIP_TO_01(1.0f - ((float)time)/((float)lifetime)), 1);}

	switch (ptype) {
	case PTYPE_GLOW:
		if (ddata.draw_as_pt()) {
			emissive_pld.add_pt(make_pt_global(pos), color); // Note: may not be in correct back to front ordering for alpha blending
		}
		else if (60.0*radius < ddata.dist) {
			glow_pld.add_pt(make_pt_global(pos), vector3d(2.0*radius, 0.0, 0.0), color); // FIXME: radius encoded as normal.x
		}
		else {
			ddata.setup_colors_draw_flare(pos, all_zeros, 2.0, 2.0, color); // Note: draw order isn't always correct
		}
		break;

	case PTYPE_SPHERE: // low resolution particles (ship pieces)
		if (ddata.draw_as_pt(0.8)) {
			if (texture_id >= 0) {color = color.modulate_with(texture_color(texture_id));}
			particle_pld.add_pt(make_pt_global(pos), (get_player_pos2() - pos), color);
		}
		else {
			color.do_glColor();
			if (texture_id >= 0) select_texture(texture_id);
			draw_sphere_vbo(all_zeros, 1.0, min((no_coll() ? ddata.ndiv : max(3, 3*ddata.ndiv/4)), N_SPHERE_DIV/2), (texture_id > 0)); // fewer ndiv/more irregular?
			if (texture_id >= 0) end_texture();
		}
		break;

	case PTYPE_TRIANGLE:
		if (ddata.draw_as_pt(0.6)) {
			vector3d normal(plus_z);
			rotate_vector3d(axis, angle, normal);
			if (dot_product(normal, vector3d(get_player_pos2() - pos)) < 0.0) {normal = -normal;} // side facing the player
			particle_pld.add_pt(make_pt_global(pos), normal, color);
		}
		// medium distance: vector of untextured triangles to draw at end?
		else {
			color.do_glColor();
			ddata.draw_one_triangle(axis, angle); // rotate around some random axis
		}
		break;
	default:
		assert(0);
	}
}


// ************ UPARTICLE_CLOUD ************


shader_t upc_shader; // FIXME: some way to pass this through uobj_draw_data so it's not a global?

void end_part_cloud_draw() {upc_shader.end_shader();}


void add_uparticle_cloud(point const &pos, float rmin, float rmax, colorRGBA const &ci1, colorRGBA const &co1,
	colorRGBA const &ci2, colorRGBA const &co2, unsigned lt, float damage, float expand_exp, float noise_scale)
{
	uparticle_cloud *upc(new uparticle_cloud(pos, rmin, rmax, ci1, co1, ci2, co2, lt, damage, expand_exp, noise_scale));
	bool const coll(add_uobj(upc, 0));
	assert(!coll);
}


uparticle_cloud::uparticle_cloud(point const &pos_, float rmin_, float rmax_, colorRGBA const &ci1, colorRGBA const &co1,
	colorRGBA const &ci2, colorRGBA const &co2, unsigned lt, float damage_, float expand_exp_, float noise_scale_)
	: lifetime(lt), rmin(rmin_), rmax(rmax_), damage_v(damage_), expand_exp(expand_exp_), noise_scale(noise_scale_)
{
	assert(rmin > 0.0 && rmin < rmax);
	free_obj::reset();
	colors[0][0] = ci1;
	colors[1][0] = co1;
	colors[0][1] = ci2;
	colors[1][1] = co2;
	flags       = OBJ_FLAGS_NCOL | OBJ_FLAGS_NOPC | OBJ_FLAGS_NEXD | OBJ_FLAGS_NOLT | OBJ_FLAGS_PARC; // not sure about OBJ_FLAGS_NOLT
	pos         = reset_pos = pos_;
	radius      = c_radius = rmin_;
	alignment   = ALIGN_NEUTRAL;
	draw_rscale = 1.0; // ?
	hashval     = 1000*pos.x;
	gen_pts(1.0); // generate the points using a radius of 1.0 and scale them to the current radius during rendering
}


void uparticle_cloud::apply_physics() {

	if (time > lifetime) {status = 1; return;}
	free_obj::apply_physics();
	radius = c_radius = rmin + (rmax - rmin)*pow(get_lt_scale(), expand_exp);
}


void uparticle_cloud::draw_obj(uobj_draw_data &ddata) const { // Note: assumes GL_BLEND is already enabled

	float const lt_scale(get_lt_scale());
	colorRGBA cur_colors[2]; // {inner, outer}
	for (unsigned d = 0; d < 2; ++d) {blend_color(cur_colors[d], colors[d][1], colors[d][0], lt_scale, 1);}
	shader_t &s(upc_shader);
	shader_setup(s, 1); // grayscale noise
	s.enable();
	s.add_uniform_float("noise_scale", noise_scale);
	s.add_uniform_color("color1i",     cur_colors[0]);
	s.add_uniform_color("color1o",     cur_colors[1]);
	s.add_uniform_float("radius",      1.0); // vertex will be scaled by radius
	s.add_uniform_float("offset",      hashval); // used as a hash
	s.add_uniform_vector3d("view_dir", (get_camera_pos() - pos).get_norm()); // local object space
	cur_colors[0].do_glColor();
	draw_quads();
	bool const use_shaders((display_mode & 0x08) != 0);
	if (use_shaders && ddata.shader) {ddata.shader->enable();}
}


// ************ US_PROJECTILE ************


us_projectile::us_projectile(unsigned type) : tup_time(0), alloc_block(NULL) {

	flags = (OBJ_FLAGS_TARG | OBJ_FLAGS_PROJ);
	set_type(type);
}


void us_projectile::set_type(unsigned type) {

	assert(type < us_weapons.size());
	wclass   = type;
	radius   = specs().radius;
	c_radius = specs().c_radius;
	armor    = specs().armor;
	if (specs().no_coll)    flags |= OBJ_FLAGS_NCOL;
	if (specs().no_exp_dam) flags |= OBJ_FLAGS_NEXD;
	if (specs().is_decoy)   flags |= OBJ_FLAGS_DECY;
	if (specs().no_light)   flags |= OBJ_FLAGS_NOLT;
	if (!us_weapons[type].hit_proj) flags |= OBJ_FLAGS_NOPC; else flags &= ~OBJ_FLAGS_NOPC;
	if ( us_weapons[type].c2_flag)  flags |= OBJ_FLAGS_NOC2; else flags &= ~OBJ_FLAGS_NOC2;
}


bool us_projectile::dec_ref() {
	
	if (alloc_block == NULL) return 0;
	alloc_block->free_obj();
	return 1;
}


us_weapon const &us_projectile::specs() const {

	assert(wclass < NUM_UWEAP);
	return us_weapons[wclass];
}


inline unsigned us_projectile::get_eflags() const {
	
	return (specs().no_ffire ? EXP_FLAGS_NO_FFIRE : 0);
}


void us_projectile::ai_action() {

	if (!is_ok() || !specs().seeking || time < PROJ_ARM_T || !begin_motion) return;
	if (rot_rate != 0.0) rot_rate *= pow(PROJ_ROT_ATTEN, fticks); // stabilize
	float const max_dist(specs().seek_dist), target_dist((target_obj == NULL) ? 0.0 : p2p_dist(pos, target_obj->get_pos()));
	free_obj const *const ptarg((parent == NULL || !target_valid(parent->get_target())) ? NULL : parent->get_target());
	assert(max_dist > 0.0);
	bool missile_lock(0);

	if (ptarg != NULL && (target_obj == NULL || target_obj == ptarg || ptarg->is_decoy()) && dist_less_than(pos, ptarg->get_pos(), max_dist)) {
		target_obj   = ptarg; // keep using parent's target or decoy since it's in seeking range
		missile_lock = 1;
	}
	else {
		if (target_obj != NULL && target_dist > 2.0*max_dist) target_obj = NULL; // too far - loose target
		
		if (target_obj == NULL || (!target_obj->is_decoy() && time > (tup_time + SEEK_CTIME)) || target_dist > max_dist) { // execute seeking code
			float seek_dist(max_dist);
			if (target_obj != NULL) seek_dist = min(seek_dist, 0.7f*target_dist); // hysteresis to keep current target
			tup_time   = time;
			target_obj = get_closest_ship(pos, 0.0, seek_dist, 1, 0, 0, 1);
			if (target_obj != NULL) missile_lock = 1;
			bool const decoy(target_obj != NULL && target_obj->is_decoy()); // Note: Decoy will only work if fighting enemy teams

			if (ptarg != NULL && ptarg != target_obj && !decoy) {
				if (target_obj == NULL || (p2p_dist_sq(pos, ptarg->get_pos()) < p2p_dist_sq(pos, target_obj->get_pos()))) {
					target_obj = ptarg; // seek out your parent's current target (even if it's not an enemy)
				}
			}
		}
	}
	if (target_obj != NULL) { // target acquired
		if (target_obj->is_player_ship() && missile_lock) {
			send_warning_message(string("Missile Lock Warning: ") + get_name());
		}
		vector3d const seek_dir(target_obj->get_pos(), pos);
		
		if (dot_product(seek_dir, dir) > 0.0) { // only seek if in the same direction
			float const smag_sq(seek_dir.mag_sq());

			if (smag_sq > TOLERANCE && smag_sq < max_dist*max_dist) {
				float const smag(sqrt(smag_sq)), vmag(velocity.mag()), ss(smag/max_dist);
				velocity += seek_dir*(0.025*(1.0 - ss) + 0.1*(1.0 - ss*ss) + 0.4*(1.0 - ss*ss*ss));
				velocity.set_max_mag(max(specs().speed, vmag)); // normalize to original velocity
			}
			dir = velocity.get_norm(); // force dir to align with velocity when seeking
			invalidate_rotv();
		}
	}
}


void us_projectile::apply_physics() {

	if (!is_ok()) return;
	unsigned const exp_type(specs().exp_type);

	if (time > specs().lifetime) {
		if (specs().bradius > 0.0 && !specs().no_exp_dam) {
			explode(0.4*specs().damage, 0.4*specs().bradius, exp_type,
				dir, 0, wclass, alignment, get_eflags(), parent); // smaller explosion than collision
		}
		status = 1;
	}
	else {
		if (specs().const_dam) {
			float const tval(((float)time+1)/((float)specs().lifetime+1));
			radius   = (1.0 - tval)*specs().radius + tval*specs().bradius;
			c_radius = (specs().c_radius/specs().radius)*radius;
			float const damage((0.25 + 0.75*fticks)*(1.0 - tval)*specs().damage); // hard to model this correctly since some damage is absorbed
			float const exp_size(1.5*radius);
			register_explosion(pos, exp_size, damage, (get_eflags() | EXP_FLAGS_NO_PART), wclass, this, parent);

			if (exp_type != ETYPE_NONE) {
				add_blastr(pos, dir, exp_size, 0.0, int(0.2*TICKS_PER_SECOND), alignment,
					et_params[exp_type].c1, et_params[exp_type].c2, exp_type, parent);
			}
		}
		free_obj::apply_physics();
	}
}


bool us_projectile::collision(point const &copos, vector3d const &vcoll, float obj_mass, float obj_radius, free_obj *source, float elastic) {
	
	assert(!specs().is_fighter);
	if (source != NULL && source == parent && time < W_SELF_ARM_T) return 0; // not yet armed
	float dscale(1.0);
	float const energy(coll_physics(copos, vcoll, obj_mass, obj_radius, source, elastic, &dscale));

	if (source != NULL && specs().bradius == 0.0 && dscale > 0.0 && specs().damage > 0.0 &&
		(parent == NULL || !specs().no_ffire || !parent->is_related(source)))
	{
		source->damage(dscale*specs().damage, DAMAGE_PROJ, pos, this, wclass); // damage the target
	}
	damage(energy, DAMAGE_COLL, copos, source, SWCLASS_UNDEF); // explode on contact, no bounce (in most cases)
	return 1;
}


float us_projectile::damage(float val, int type, point const &hit_pos, free_obj const *source, int wc) {

	if (val <= 0.0 || specs().no_coll) return 0.0; // no damage or can't be destroyed

	if (armor > val && ((type != DAMAGE_COLL || specs().damage == 0.0) ||
		(source != NULL && source != target_obj && source->is_proj() && !source->is_ship() && !source->is_decoy())))
	{
		armor -= val;
		return val; // hasn't collided with its target, it survives
	}
	unsigned const etype(specs().exp_type);

	if (specs().bradius > 0.0) {
		float const damage_done(specs().no_exp_dam ? 0.0 : specs().damage);
		explode(damage_done, specs().bradius, etype, dir, specs().btime, wclass, alignment, get_eflags(), parent);

		if (sobj_coll_tid >= 0 && damage_done >= 20.0) { // exploded on a planet, moon, or asteroid
			point const hit_pos2(pos + (pos - hit_pos).get_norm()*radius);
			// FIXME: get_fragment_tid() called on the wrong object - should be called on the uobject that was collided with
			gen_moving_fragments(hit_pos2, min(25U, max(1U, unsigned(min(10.0f, 0.4f*damage_done)))), sobj_coll_tid, 3.2, 6.0);
		}
	}
	else if (etype != ETYPE_NONE) {
		explode(0.0, 4.0*radius, etype, dir, 0, wclass, alignment, get_eflags(), parent);
	}
	status = 1;
	return armor;
}


void us_projectile::draw_obj(uobj_draw_data &ddata) const { // front is in -z

	if (specs().auto_orient != 0) {ddata.dir = get_orient();}
	ddata.lifetime = (unsigned)specs().lifetime; // why is lifetime a float?

	switch (wclass) {
	case UWEAP_ROCKET:  ddata.draw_usw_rocket();   break;
	case UWEAP_NUKEDEV: ddata.draw_usw_nukedev();  break;
	case UWEAP_TORPEDO: ddata.draw_usw_torpedo();  break;
	case UWEAP_ENERGY:  ddata.draw_usw_energy();   break;
	case UWEAP_ATOMIC:  ddata.draw_usw_atomic();   break;
	case UWEAP_EMP:     ddata.draw_usw_emp();      break;
	case UWEAP_DFLARE:  ddata.draw_usw_dflare();   break;
	case UWEAP_CHAFF:   ddata.draw_usw_chaff();    break;
	case UWEAP_RFIRE:   ddata.draw_usw_rfire();    break;
	case UWEAP_SHIELDD: ddata.draw_usw_shieldd();  break;
	case UWEAP_THUNDER: ddata.draw_usw_thunder();  break;
	case UWEAP_STAR:    ddata.draw_usw_star();     break;
	case UWEAP_SEIGEC:  ddata.draw_usw_seige();    break;
	default: assert(0); // otherwise can't draw here - not a real projectile (beams, etc.)
	}
}



