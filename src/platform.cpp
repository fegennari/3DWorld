// 3D World - Member function definitions for platform classes
// by Frank Gennari
// 9/25/09

#include "collision_detect.h"
#include "openal_wrap.h"
#include "model3d.h"  // for geom_xform_t
#include "lightmap.h" // for light_source_trig


platform_cont platforms;

extern bool user_action_key;
extern int animate2;
extern float fticks, CAMERA_RADIUS;
extern coll_obj_group coll_objects;
extern vector<light_source_trig> light_sources_d;


unsigned trigger_t::register_player_pos(point const &p, float act_radius, int activator, bool clicks) {

	// Note: since only the camera/player can issue an action, we assume requires_action implies player_only
	if ((player_only || requires_action) && activator != CAMERA_ID) return 0; // not activated by player

	if (requires_action) {
		if (!user_action_key) return 0; // check action key
		if (!use_act_region && !camera_pdu.point_visible_test(act_pos)) return 0; // player not looking at the activation pos
	}
	if (use_act_region) {return (act_region.contains_pt(p) ? (requires_action ? 2 : 1) : 0);} // check active region containment
	else if (act_dist == 0.0) {return 0;} // act_dist of 0 disables this trigger
	if (!dist_less_than(p, act_pos, (act_dist + act_radius))) {return 0;} // too far to activate
	if (requires_action && clicks) {gen_sound(SOUND_CLICK, act_pos, 1.0);}
	return (requires_action ? 2 : 1); // triggered by proximity = 1, triggered by player action = 2
}


unsigned multi_trigger_t::register_player_pos(point const &p, float act_radius, int activator, bool clicks) {
	unsigned ret(0);
	for (iterator i = begin(); i != end(); ++i) {ret |= i->register_player_pos(p, act_radius, activator, clicks);}
	return ret;
}

void multi_trigger_t::shift_by(vector3d const &val) {
	for (iterator i = begin(); i != end(); ++i) {i->shift_by(val);}
}

float multi_trigger_t::get_auto_on_time() const { // the first trigger to activate activates the group
	float min_aot(0.0);
	for (const_iterator i = begin(); i != end(); ++i) {
		if (i->auto_on_time > 0.0) {min_aot = ((min_aot == 0.0) ? i->auto_on_time : min(min_aot, i->auto_on_time));}
	}
	return min_aot;
}

float multi_trigger_t::get_auto_off_time() const { // the last trigger to deactivate deactivates the group
	float max_aot(0.0);
	for (const_iterator i = begin(); i != end(); ++i) {max_aot = max(max_aot, i->auto_off_time);}
	return max_aot;
}


platform::platform(float fs, float rs, float sd, float rd, float dst, float ad, point const &o, vector3d const &dir_, bool c, bool ir)
				   : cont(c), is_rot(ir), fspeed(fs), rspeed(rs), sdelay(sd), rdelay(rd), ext_dist(dst), act_dist(ad), origin(o), dir(dir_.get_norm()), delta(all_zeros)
{
	assert(dir_ != all_zeros);
	assert(fspeed > 0.0 && rspeed > 0.0 && sdelay >= 0.0 && ext_dist > 0.0 && act_dist >= 0.0);
	if (act_dist > 0.0) {triggers.push_back(trigger_t(origin, act_dist));}
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

	if (cont || state != ST_NOACT || empty()) return 1; // continuous, already activated, or no cobjs/lights
	if (!triggers.register_player_pos(p, radius, activator, 1)) return 0; // not yet triggered
	activate();
	return 1;
}


void platform::next_frame() {
	
	cobjs.clear(); // lights is constant
	delta = all_zeros;
}


void platform::move_platform(float dist_traveled) {
	if (is_rot) {assert(0);} // FIXME: WRITE - difficult because rotations may change cobj type (such as cube -> extruded polygon)
	else {pos += dir*dist_traveled;}
}


void platform::advance_timestep() {

	if (fticks == 0.0 || empty()) return; // no progress or no cobjs/lights
	
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
					else {ns_time = 0.0;} // keep moving forward
					move_platform(dist_traveled);
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
						move_platform(dist_traveled);
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
			cobj.move_cobj(delta, update_colls); // move object
			// squish player or stop when hit player?
		}
		for (vector<unsigned>::const_iterator i = lights.begin(); i != lights.end(); ++i) {
			assert(*i < light_sources_d.size());
			light_sources_d[*i].shift_by(delta);
		}
	} // pos != last_pos
}


vector3d platform::get_velocity() const { // Note: linear velocity, or angular velocity if is_rot==1

	if      (state == ST_FWD) {return dir* fspeed;}
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
	triggers.shift_by(val);
}


bool platform_cont::add_from_file(FILE *fp, geom_xform_t const &xf, multi_trigger_t const &triggers) {

	assert(fp);
	float fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist; // in seconds/units-per-second
	point origin;
	vector3d dir; // or rot_axis
	int cont(0), is_rotation(0);
	
	if (fscanf(fp, "%f%f%f%f%f%f%f%f%f%f%f%f%i", &fspeed, &rspeed, &sdelay, &rdelay, &ext_dist,
		&act_dist, &origin.x, &origin.y, &origin.z, &dir.x, &dir.y, &dir.z, &cont) != 13)
	{
		return 0;
	}
	if (cont != 0 && cont != 1) return 0; // not a bool
	fscanf(fp, "%i", &is_rotation); // try to read is_rotation - if fails, leave at 0
	if (is_rotation != 0 && is_rotation != 1) return 0; // not a bool
	sdelay *= TICKS_PER_SECOND;
	rdelay *= TICKS_PER_SECOND;
	fspeed /= TICKS_PER_SECOND;
	rspeed /= TICKS_PER_SECOND;
	xf.xform_pos(origin);
	xf.xform_pos_rm(dir);
	push_back(platform(fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist, origin, dir, (cont != 0)));
	if (!triggers.empty()) {back().add_triggers(triggers);} // if a custom trigger is used, reset any built-in trigger
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

bool platform_cont::any_active() const {

	for (auto i = begin(); i != end(); ++i) {if (i->is_moving()) return 1;}
	return 0;
}

bool platform_cont::any_moving_platforms_in_view(pos_dir_up const &pdu) const { // Note: untested

	for (auto i = begin(); i != end(); ++i) {
		if (!i->is_moving()) continue;

		for (auto j = i->cobjs.begin(); j != i->cobjs.end(); ++j) {
			coll_obj const &cobj(coll_objects[*j]);
			if (!cobj.disabled() && pdu.cube_visible(cobj)) return 1; // test cube in view frustum
		}
	}
	return 0;
}


void coll_obj::add_to_platform() const {

	if (platform_id < 0) return; // no platform
	assert(size_t(platform_id) < platforms.size());
	platforms[platform_id].add_cobj(id);
}


