// 3D World - Member function definitions for platform classes
// by Frank Gennari
// 9/25/09

#include "collision_detect.h"


platform_cont platforms;

extern int animate2, enable_shadow_maps;
extern float fticks;
extern vector<coll_obj> coll_objects;



// shadow_mode: 0 = no dynamic shadows, 1 = only on object, 2 = on object and targets
platform::platform(float fs, float rs, float sd, float rd, float dst, float ad,
				   point const &o, vector3d const &dir_, int sm, bool c)
				   : cont(c), shadow_mode(sm), fspeed(fs), rspeed(rs), sdelay(sd), rdelay(rd), ext_dist(dst),
				   act_dist(ad), origin(o), dir(dir_.get_norm()), delta(all_zeros), s_d_chg(0)
{
	assert(shadow_mode <= 2);
	assert(dir_ != all_zeros);
	assert(fspeed > 0.0 && rspeed > 0.0 && sdelay >= 0.0 && ext_dist > 0.0 && act_dist >= 0.0);
	reset();
}


void platform::reset() {

	state   = ST_NOACT;
	ns_time = 0.0;
	pos     = origin;
	s_d_chg = !cont;
}


void platform::activate() {

	assert(state == ST_NOACT);
	assert(pos == origin);
	assert(ns_time == 0.0);
	state   = ST_WAIT; // activated
	ns_time = sdelay;
	s_d_chg = !cont;
}


bool platform::check_activate(point const &p, float radius) {

	if (cont || state != ST_NOACT || cobjs.empty()) return 1; // continuous, already activated, or no cobjs
	
	if (dist_less_than(p, pos, (act_dist + radius))) {
		activate();
		return 1;
	}
	return 0;
}


void platform::next_frame() {
	
	cobjs.clear();
	delta = all_zeros;
}


void platform::advance_timestep() {

	if (fticks == 0.0 || cobjs.empty()) return; // no progress or no cobjs
	
	if (state == ST_NOACT) { // not activated
		assert(pos == origin);
		assert(ns_time == 0.0);
		if (!cont) return;
		activate();
	}
	ns_time -= fticks;
	point const last_pos(pos);

	while (ns_time < 0.0) { // guaranteed to terminate as long as the assertions in the constructor are met
		switch (state) {
			case ST_WAIT: // activated, waiting to start moving forward
				assert(pos == origin);
				state = ST_FWD; // wait has ended, start moving forward (fallthrough)
			case ST_FWD: // moving forward
				{
					float dist_traveled(-fspeed*ns_time), cur_dist(p2p_dist(pos, origin)); // dist is pos
					assert(dist_traveled > 0.0);
					
					if (dist_traveled + cur_dist > ext_dist) { // traveled past the end
						dist_traveled = ext_dist - cur_dist;
						ns_time      += dist_traveled/fspeed; // ns_time will generally still be neg at this step
						ns_time      += max(0.0f, rdelay); // add in reverse delay (ns_time can be pos or neg)
						state         = ST_CHDIR;
					}
					else { // keep moving forward
						ns_time = 0.0;
					}
					pos += dir*dist_traveled;
				}
				break;

			case ST_CHDIR: // waiting to change direction
				if (rdelay < 0.0) return; // no reverse phase - stay in this state forever
				state = ST_REV; // time is up, reverse (fallthrough)
			case ST_REV: // moving in reverse
				{
					float dist_traveled(rspeed*ns_time), cur_dist(p2p_dist(pos, origin)); // dist is neg
					assert(dist_traveled < 0.0);
					
					if (dist_traveled + cur_dist < 0.0) { // traveled past the start
						// we might have some time left over this frame: (ns_time - cur_dist/fspeed)
						// if in continuous mode, we can have multiple iterations per frame if fticks is large
						// but we don't want to do this for efficiency/complexity reasons
						reset(); // reset time to start a new iteration
					}
					else { // keep moving in reverse
						ns_time = 0.0;
						pos    += dir*dist_traveled;
					}
				}
				break;

			case ST_NOACT: // not activated
			default:
				assert(0);
		}
	}
	if (pos != last_pos) {
		delta = (pos - last_pos); // can accumulate error, fix?

		for (vector<unsigned>::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
			assert(*i < coll_objects.size());
			coll_obj &cobj(coll_objects[*i]);
			// need to update collision structure when there is an x/y delta by removing/adding to coll_cells (except for cubes)
			bool const update_colls(cobj.type != COLL_CUBE && (delta.x != 0.0 || delta.y != 0.0));

			if (shadow_mode > 0) {
				if (shadow_mode > 1 && s_d_chg) { // static/dynamic state change
					cobj.clear_dependent_cobjs_lightmaps(coll_objects, *i);
				}
				cobj.clear_lightmap(0, 0, 1);
			}
			if (update_colls) remove_coll_object(*i, 0);
			cobj.shift_by(delta); // move object
			if (update_colls) cobj.re_add_coll_cobj(*i, 0);
			// squish player or stop when hit player?
			// what about coll_cell::cvz?
		}
		s_d_chg = 0;
	} // pos != last_pos
}


vector3d platform::get_velocity() const {

	if (state == ST_FWD) {
		return dir*fspeed;
	}
	else if (state == ST_REV) {
		return dir*-rspeed;
	}
	return zero_vector;
}


void platform::add_cobj(unsigned cobj) {
	
	assert(cobj < coll_objects.size());
	cobjs.push_back(cobj);
}


void platform::shift_by(vector3d const &val) {
	
	origin += val;
	pos    += val;
}


bool platform_cont::add_from_file(FILE *fp) {

	assert(fp);
	float fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist;
	point origin;
	vector3d dir;
	int cont, shadow_mode;
	
	if (fscanf(fp, "%f%f%f%f%f%f%f%f%f%f%f%f%i%i", &fspeed, &rspeed, &sdelay, &rdelay, &ext_dist,
		&act_dist, &origin.x, &origin.y, &origin.z, &dir.x, &dir.y, &dir.z, &shadow_mode, &cont) != 14)
	{
		return 0;
	}
	if (cont != 0 && cont != 1) return 0; // not a bool
	push_back(platform(fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist, origin, dir, shadow_mode, (cont != 0)));
	return 1;
}


void platform_cont::check_activate(point const &p, float radius) {

	for (unsigned i = 0; i < size(); ++i) {
		operator[](i).check_activate(p, radius);
	}
}


void platform_cont::shift_by(vector3d const &val) {

	for (unsigned i = 0; i < size(); ++i) {
		operator[](i).shift_by(val);
	}
}


void platform_cont::advance_timestep() {

	for (unsigned i = 0; i < size(); ++i) { // cache this?
		operator[](i).next_frame();
	}
	for (unsigned i = 0; i < coll_objects.size(); ++i) {
		int const pid(coll_objects[i].platform_id);
		if (pid < 0) continue;
		assert(pid < (int)size());
		operator[](pid).add_cobj(i);
	}
	for (unsigned i = 0; i < size(); ++i) {
		operator[](i).advance_timestep();
	}
}


void coll_obj::add_to_platform() const {

	if (platform_id < 0) return; // no platform
	assert(size_t(platform_id) < platforms.size());
	platforms[platform_id].add_cobj(id);
}


bool coll_obj::dynamic_shadows_only() const {

	if (cp.color.alpha < 0.5) return 0; // semi-transparent
	if (platform_id    < 0)   return 0;
	assert(platform_id < (int)platforms.size());
	return (enable_shadow_maps == 2 || platforms[platform_id].has_dynamic_shadows());
}



