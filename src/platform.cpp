// 3D World - Member function definitions for platform classes
// by Frank Gennari
// 9/25/09

#include "collision_detect.h"
#include "openal_wrap.h"
#include "model3d.h" //for geom_xform_t


platform_cont platforms;

extern bool user_action_key;
extern int animate2;
extern float fticks, CAMERA_RADIUS;
extern coll_obj_group coll_objects;


bool trigger_t::register_player_pos(point const &p, float act_radius, int activator, bool clicks) {

	// Note: since only the camera/player can issue an action, we assume requires_action implies player_only
	if ((player_only || requires_action) && activator != CAMERA_ID) return 0; // not activated by player

	if (requires_action) {
		if (!user_action_key) return 0; // check action key
		if (!use_act_region && !camera_pdu.point_visible_test(act_pos)) return 0; // player not looking at the activation pos
	}
	if (use_act_region) {return act_region.contains_pt(p);} // check active region containment
	else if (act_dist == 0.0) {return 0;} // act_dist of 0 disables this trigger
	if (!dist_less_than(p, act_pos, (act_dist + act_radius))) {return 0;} // too far to activate
	if (requires_action && clicks) {gen_sound(SOUND_CLICK, act_pos, 1.0);}
	return 1;
}


platform::platform(float fs, float rs, float sd, float rd, float dst, float ad,
				   point const &o, vector3d const &dir_, bool c)
				   : cont(c), fspeed(fs), rspeed(rs), sdelay(sd), rdelay(rd), ext_dist(dst),
				   act_dist(ad), origin(o), dir(dir_.get_norm()), delta(all_zeros), trigger(origin, act_dist, 0)
{
	assert(dir_ != all_zeros);
	assert(fspeed > 0.0 && rspeed > 0.0 && sdelay >= 0.0 && ext_dist > 0.0 && act_dist >= 0.0);
	reset();
}


void platform::reset() {

	state   = ST_NOACT;
	ns_time = 0.0;
	pos     = origin;
}


void platform::activate() {

	assert(state == ST_NOACT);
	assert(pos == origin);
	assert(ns_time == 0.0);
	state   = ST_WAIT; // activated
	ns_time = sdelay;
}


bool platform::check_activate(point const &p, float radius, int activator) {

	if (cont || state != ST_NOACT || cobjs.empty()) return 1; // continuous, already activated, or no cobjs
	if (!trigger.register_player_pos(p, radius, activator, 1)) return 0; // not yet triggered
	activate();
	return 1;
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
			if (update_colls) remove_coll_object(*i, 0);
			cobj.shift_by(delta); // move object
			if (update_colls) cobj.re_add_coll_cobj(*i, 0);
			// squish player or stop when hit player?
		}
	} // pos != last_pos
}


vector3d platform::get_velocity() const {

	if (state == ST_FWD) {return dir*fspeed;}
	else if (state == ST_REV) {return dir*-rspeed;}
	return zero_vector;
}


void platform::add_cobj(unsigned cobj) {
	
	assert(cobj < coll_objects.size());
	cobjs.push_back(cobj);
}


void platform::shift_by(vector3d const &val) {
	
	origin += val;
	pos    += val;
	trigger.shift_by(val);
}


bool platform_cont::add_from_file(FILE *fp, geom_xform_t const &xf, trigger_t const &trigger) {

	assert(fp);
	float fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist; // in seconds/units-per-second
	point origin;
	vector3d dir;
	int cont;
	
	if (fscanf(fp, "%f%f%f%f%f%f%f%f%f%f%f%f%i", &fspeed, &rspeed, &sdelay, &rdelay, &ext_dist,
		&act_dist, &origin.x, &origin.y, &origin.z, &dir.x, &dir.y, &dir.z, &cont) != 13)
	{
		return 0;
	}
	if (cont != 0 && cont != 1) return 0; // not a bool
	sdelay *= TICKS_PER_SECOND;
	rdelay *= TICKS_PER_SECOND;
	fspeed /= TICKS_PER_SECOND;
	rspeed /= TICKS_PER_SECOND;
	xf.xform_pos(origin);
	xf.xform_pos_rm(dir);
	push_back(platform(fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist, origin, dir, (cont != 0)));
	if (trigger.is_active()) {back().set_trigger(trigger);}
	return 1;
}


void platform_cont::check_activate(point const &p, float radius, int activator) {
	for (auto i = begin(); i != end(); ++i) {i->check_activate(p, radius, activator);}
}

void platform_cont::shift_by(vector3d const &val) {
	for (auto i = begin(); i != end(); ++i) {i->shift_by(val);}
}


void platform_cont::advance_timestep() {

	for (auto i = begin(); i != end(); ++i) {i->next_frame();} // cache this?

	for (cobj_id_set_t::const_iterator i = coll_objects.platform_ids.begin(); i != coll_objects.platform_ids.end(); ++i) {
		assert(*i < coll_objects.size());
		int const pid(coll_objects[*i].platform_id);
		assert(pid >= 0 && pid < (int)size());
		operator[](pid).add_cobj(*i);
	}
	for (auto i = begin(); i != end(); ++i) {i->advance_timestep();}
}


void coll_obj::add_to_platform() const {

	if (platform_id < 0) return; // no platform
	assert(size_t(platform_id) < platforms.size());
	platforms[platform_id].add_cobj(id);
}



