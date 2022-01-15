// 3D World - Member function definitions for platform classes
// by Frank Gennari
// 9/25/09

#include "collision_detect.h"
#include "openal_wrap.h"
#include "model3d.h"  // for geom_xform_t
#include "lightmap.h" // for light_source_trig


platform_cont platforms;

extern bool user_action_key;
extern int display_mode;
extern float fticks, CAMERA_RADIUS, temperature;
extern coll_obj_group coll_objects;
extern set<unsigned> moving_cobjs;
extern vector<light_source_trig> light_sources_d;
extern vector<sphere_t> cur_frame_explosions;
extern obj_vector_t<fire> fires;


// ***** triggers *****

// Note: activator can be player (CAMERA_ID), smileys (0-N), lights (0-N), or explosions (NO_SOURCE)
unsigned trigger_t::register_activator_pos(point const &p, float act_radius, int activator, bool clicks) const {

	// Note: since only the camera/player can issue an action, we assume requires_action implies player_only
	if (player_only && activator != CAMERA_ID) return 0; // not activated by player
	if (req_keycard_id >= 0 && !has_keycard_id(activator, req_keycard_id)) return 0; // activator doesn't have the required keycard
	bool const is_explosion_action(activator == NO_SOURCE);

	if (requires_action && !is_explosion_action) { // requires explicit player action
		if (activator != CAMERA_ID) return 0; // not activated by player
		if (!user_action_key) return 0; // check action key
		if (!use_act_region && !camera_pdu.point_visible_test(act_pos)) return 0; // player not looking at the activation pos
	}
	if (use_act_region) {return (act_region.contains_pt(p) ? (requires_action ? 2 : 1) : 0);} // check active region containment
	else if (act_dist == 0.0) {return 0;} // act_dist of 0 disables this trigger
	if (!dist_less_than(p, act_pos, (act_dist + act_radius))) {return 0;} // too far to activate
	if (requires_action && clicks) {gen_sound(SOUND_CLICK, act_pos, 1.0);}
	return (requires_action ? 2 : 1); // triggered by proximity = 1, triggered by player action = 2
}

void trigger_t::write_to_cobj_file(std::ostream &out) const {
	out << "trigger " << act_pos.raw_str() << " " << act_dist << " " << auto_on_time << " " << auto_off_time
		<< " " << (player_only != 0) << " " << (requires_action != 0) << " " << (player_only ? req_keycard_id : obj_type_id);
	if (use_act_region) {out << " " << act_region.raw_str();}
	out << endl;
}


unsigned multi_trigger_t::register_activator_pos(point const &p, float act_radius, int activator, bool clicks) {
	unsigned ret(0);
	for (iterator i = begin(); i != end(); ++i) {ret |= i->register_activator_pos(p, act_radius, activator, clicks);}
	return ret;
}

unsigned multi_trigger_t::check_for_activate_this_frame() {
	unsigned ret(0);
	for (iterator i = begin(); i != end(); ++i) {ret |= i->check_for_activate_this_frame();}
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

void multi_trigger_t::write_to_cobj_file(std::ostream &out) const {
	for (const_iterator i = begin(); i != end(); ++i) {i->write_to_cobj_file(out);}
}
void multi_trigger_t::write_end_triggers_cobj_file(std::ostream &out) const {
	if (!empty()) {out << "trigger" << endl;} // empty/end trigger
}


// ***** sensors *****

bool check_for_light(point const &pos, float thresh) {
	if (is_visible_to_any_dir_light(pos, 0.0, -1, 0)) return 1; // check sun and moon; skip_dynamic=0
	if (thresh > 1.0) return 0; // only sun/moon
	if (display_mode & 0x0200) return 0; // skip checking of dynamic lights when dynamic spheres are enabled, since they uselessly trigger light sensors
	return is_any_dlight_visible(pos); // Note: thresh is unused
}

bool check_for_heat(point const &pos, float radius, float thresh) {

	if (temperature > thresh) return 1;

	for (auto i = cur_frame_explosions.begin(); i != cur_frame_explosions.end(); ++i) { // check explosions (one frame)
		if (dist_less_than(pos, i->pos, (radius + i->radius))) return 1; // "infinitely" hot
	}
	for (auto i = fires.begin(); i != fires.end(); ++i) { // check fires
		if (i->enabled() && dist_less_than(pos, i->pos, (radius + i->radius))) return 1; // "infinitely" hot
	}
	if (get_ground_fire_intensity(pos, radius) > 0.0) return 1; // "infinitely" hot
	return 0;
}

bool check_for_metal(point const &pos, float radius) {

	cube_t bcube;
	bcube.set_from_sphere(pos, radius);
	vector<unsigned> cobjs;

	for (unsigned dynamic = 0; dynamic < 2; ++dynamic) { // check both dynamic (material spheres) and static (movable) cobjs
		cobjs.clear();
		get_intersecting_cobjs_tree(bcube, cobjs, -1, 0.0, (dynamic != 0), 0);

		for (vector<unsigned>::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
			coll_obj const &cobj(coll_objects.get_cobj(*i));
			if (cobj.cp.metalness == 0.0) continue;
			if (!cobj.may_be_dynamic())   continue; // assume fully static cobjs don't count, otherwise the result will always be the same
			if (cobj.sphere_intersects(pos, radius)) return 1; // most expensive test last
		}
	}
	return 0;
}

bool check_for_pressure(point const &pos, float radius) {

	if (check_player_proximity(pos, radius, 1)) return 1; // player/smiley interaction (using player bottom sphere)
	
	for (auto i = moving_cobjs.begin(); i != moving_cobjs.end(); ++i) {
		coll_obj const &cobj(coll_objects.get_cobj(*i));
		if (cobj.sphere_intersects(pos, radius)) return 1; // movable cobj interaction
	}
	// what about dynamic objects (balls, material spheres, etc.)?
	return 0;
}

bool sensor_t::check_active_int() const {

	switch (type) {
	case SENSOR_DISABLED:   assert(0); // shouldn't get here
	case SENSOR_ALWAYS_OFF: return 0;
	case SENSOR_ALWAYS_ON:  return 1;
	case SENSOR_LIGHT:      return check_for_light(pos, thresh);
	case SENSOR_SOUND:      return check_for_active_sound(pos, radius, thresh);
	case SENSOR_HEAT:       return check_for_heat(pos, radius, thresh);
	case SENSOR_METAL:      return check_for_metal(pos, radius);
	case SENSOR_WATER:      return is_underwater(pos);
	case SENSOR_PRESSURE:   return check_for_pressure(pos, radius);
	case SENSOR_SMOKE:      return (get_smoke_at_pos(pos) > thresh);
	default: assert(0);
	}
	return 0;
}

//SENSOR_DISABLED=0, SENSOR_ALWAYS_OFF, SENSOR_ALWAYS_ON, SENSOR_LIGHT, SENSOR_SOUND, SENSOR_HEAT, SENSOR_METAL, SENSOR_WATER, SENSOR_PRESSURE, SENSOR_SMOKE, NUM_SENSOR_TYPES
string const sensor_type_names[NUM_SENSOR_TYPES] = {"disabled", "off", "on", "light", "sound", "heat", "metal", "water", "pressure", "smoke"};

bool sensor_t::read_from_file(FILE *fp, geom_xform_t const &xf) {
	// sensor type [pos.x pos.y pos.z [invert [radius [thresh]]]]
	int inv(0);
	char str[MAX_CHARS] = {0};
	if (fscanf(fp, "%255s", str) < 1) return 0;
	int const ix(atoi(str));
	if (ix > 0) {type = ix;} // a number was specified
	else {
		string const type_str(str);
		if (type_str == "0") {type = ix;} // a number was specified
		else {
			type = -1;
			for (unsigned i = 0; i < NUM_SENSOR_TYPES; ++i) {
				if (type_str == sensor_type_names[i]) {type = i; break;}
			}
			if (type < 0) {cerr << "Error: Unknown sensor type: '" << type_str << "'" << endl; return 0;}
		}
	}
	if (type < SENSOR_DISABLED || type >= NUM_SENSOR_TYPES) {cerr << "Error: Invalid sensor type: " << type << endl; return 0;}
	if (type >= SENSOR_LIGHT && fscanf(fp, "%f%f%f%i%f%f", &pos.x, &pos.y, &pos.z, &inv, &radius, &thresh) < 3) return 0;
	invert = (inv != 0);
	xf.xform_pos(pos);
	radius *= xf.scale;
	return 1;
}

void sensor_t::write_to_cobj_file(std::ostream &out) const {
	if (type == SENSOR_DISABLED) return; // nothing to write
	out << "sensor " << type;
	if (type >= SENSOR_LIGHT) {out << " " << pos.raw_str() << " " << invert << " " << radius << " " << thresh;}
	out << endl;
}
void sensor_t::write_end_sensor_to_cobj_file(std::ostream &out) const {
	if (type == SENSOR_DISABLED) return; // nothing to write
	out << "sensor disabled" << endl;
}

bool multi_sensor_t::enabled() const {
	for (auto i = begin(); i != end(); ++i) {
		if (i->enabled()) return 1;
	}
	return 0;
}
bool multi_sensor_t::check_active() const {
	for (auto i = begin(); i != end(); ++i) {
		if (i->check_active()) return 1;
	}
	return 0;
}
void multi_sensor_t::write_to_cobj_file(std::ostream &out) const { // one per line
	for (auto i = begin(); i != end(); ++i) {i->write_to_cobj_file(out);}
	//out << "sensor disabled" << endl; // should there be an option for this?
}


// ***** platforms *****

platform::platform(float fs, float rs, float sd, float rd, float dst, float ad, point const &o, vector3d const &dir_, bool c, bool ir, bool ul, bool destroys_, int sid, sensor_t const &cur_sensor) :
	cont(c), is_rot(ir), update_light(ul), destroys(destroys_), fspeed(fs), rspeed(rs), sdelay(sd), rdelay(rd), ext_dist(dst), act_dist(ad),
	origin(o), dir(dir_.get_norm()), sound_id(sid), delta(all_zeros), sensor(cur_sensor)
{
	assert(dir_ != all_zeros);
	assert(fspeed > 0.0 && sdelay >= 0.0 && act_dist >= 0.0);
	assert(is_rot || ext_dist > 0.0);
	if (act_dist > 0.0) {triggers.push_back(trigger_t(origin, act_dist));}
	reset();
}

void platform::reset() {

	state      = ST_NOACT;
	is_stopped = 0;
	is_paused  = 0;
	ns_time    = active_time = 0.0;
	pos        = origin;
	cur_angle  = 0.0;
}

void platform::activate() {

	assert(state == ST_NOACT);
	assert(pos == origin && cur_angle == 0.0);
	assert(ns_time == 0.0);
	state   = ST_WAIT; // activated
	is_stopped = 0;
	ns_time = sdelay;
	active_time = 0.0;
}

void platform::pause() {

	assert(!is_paused);
	state      = ST_NOACT;
	is_stopped = 0;
	is_paused  = 1;
}

void platform::unpause() {

	assert(is_paused);
	state      = ST_FWD; // Note: currently can only pause while moving forward
	is_stopped = 0;
	is_paused  = 0;
}

bool platform::check_activate(point const &p, float radius, int activator) {

	if (cont || empty()) return 1; // continuous or no cobjs/lights

	if (state != ST_NOACT && !is_stopped && !is_paused) { // already activated
		if (ext_dist == 0.0) { // unbounded rotation
			assert(is_rot);
			if (triggers.register_activator_pos(p, radius, activator, 1) == 2) {pause(); check_play_sound(); return 0;} // deactivate/pause
		}
		return 1;
	}
	if (!triggers.register_activator_pos(p, radius, activator, 1)) return 0; // not yet triggered
	register_activate();
	return 1;
}

void platform::register_activate() {
	if (is_stopped) {is_stopped = 0; check_play_sound();}
	else if (is_paused) {unpause(); check_play_sound();}
	else {activate();}
}

void platform::move_platform(float dist_traveled) { // linear distance or rotation angle
	
	if (is_rot) {cur_angle += dist_traveled;} // rotate
	else        {pos       += dir*dist_traveled; return;} // translate

	for (auto i = cobjs.begin(); i != cobjs.end(); ++i) { // handle rotation
		coll_objects.get_cobj(*i).rotate_about(origin, dir, dist_traveled, 1);
	}
}

void platform::check_play_sound() const {
	if (sound_id >= 0) {gen_sound(sound_id, pos, 1.0);}
}

void platform::advance_timestep() {

	if (fticks == 0.0 || empty() || is_paused) return; // no progress, no cobjs/lights, or paused
	if (state == ST_NOACT && triggers.check_for_activate_this_frame()) {register_activate();} // not active - check if activated by an object
	if (is_stopped && !is_sensor_active()) return; // stopped and waiting for trigger re-activate
	is_stopped = 0;
	
	if (state == ST_NOACT) { // not activated
		assert(pos == origin && cur_angle == 0.0);
		assert(ns_time == 0.0);
		if (!cont && !is_sensor_active()) return;
		activate();
	}
	float const auto_off_time(triggers.get_auto_off_time());

	if (auto_off_time > 0.0) { // automatically turns off after a fixed period of time
		if (active_time > TICKS_PER_SECOND*auto_off_time) {active_time = 0.0; is_stopped = 1; return;} // reached auto off time, stop
		active_time += fticks;
	}
	ns_time -= fticks;
	point const last_pos(pos);

	while (ns_time < 0.0) { // guaranteed to terminate as long as the assertions in the constructor are met
		switch (state) {
			case ST_WAIT: // activated, waiting to start moving forward
				assert(pos == origin && cur_angle == 0.0);
				state = ST_FWD; // wait has ended, start moving forward (fallthrough)
				check_play_sound();
			case ST_FWD: // moving forward
				{
					float dist_traveled(-fspeed*ns_time), cur_dist(get_dist_traveled()); // dist is pos
					assert(dist_traveled > 0.0);
					
					if (ext_dist > 0.0 && dist_traveled + cur_dist > ext_dist) { // traveled past the end
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
				check_play_sound();
			case ST_REV: // moving in reverse
				{
					if (rspeed == 0.0) { // no reverse state
						if (cont && is_rot) {cur_angle = 0.0; state = ST_FWD;} // back to forward state - intinite loop
						else {} // wait in this state forever
						ns_time = 0.0;
						break;
					}
					float dist_traveled(rspeed*ns_time), cur_dist(get_dist_traveled()); // dist is neg
					assert(dist_traveled < 0.0);
					
					if (dist_traveled + cur_dist < 0.0f) { // traveled past the start
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
			coll_obj &cobj(coll_objects.get_cobj(*i));
			// need to update collision structure when there is an x/y delta by removing/adding to coll_cells (except for cubes)
			bool const update_colls(cobj.type != COLL_CUBE && (delta.x != 0.0 || delta.y != 0.0));
			cobj.move_cobj(delta, update_colls); // move object
			// may crush the player (see vert_coll_detector::check_cobj_intersect()); could also maybe make the platform stop when hitting the player
			if (destroys) {destroy_coll_objs(cobj.get_center_pt(), 1000.0, NO_SOURCE, IMPACT, cobj.get_bsphere_radius());}
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

void platform_cont::read_sound_filename(string const &name) {cur_sound_id = read_sound_file(name);}

bool platform_cont::add_from_file(FILE *fp, geom_xform_t const &xf, multi_trigger_t const &triggers, sensor_t const &cur_sensor) {

	assert(fp);
	float fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist; // in seconds/units-per-second
	point origin;
	vector3d dir; // or rot_axis
	int cont(0), is_rotation(0), update_light(0), destroys(0);
	// Q enabled [fspeed rspeed sdelay rdelay ext_dist act_dist origin<x,y,z> dir<x,y,z> cont [is_rotation=0 [update_light=0 [destroys=0]]]]
	if (fscanf(fp, "%f%f%f%f%f%f%f%f%f%f%f%f%i%i%i%i", &fspeed, &rspeed, &sdelay, &rdelay, &ext_dist,
		&act_dist, &origin.x, &origin.y, &origin.z, &dir.x, &dir.y, &dir.z, &cont, &is_rotation, &update_light, &destroys) < 13)  {return 0;}
	if ((cont != 0 && cont != 1) || (is_rotation != 0 && is_rotation != 1) || (update_light != 0 && update_light != 1)) return 0; // not a bool
	sdelay *= TICKS_PER_SECOND;
	rdelay *= TICKS_PER_SECOND;
	fspeed /= TICKS_PER_SECOND;
	rspeed /= TICKS_PER_SECOND;
	xf.xform_pos(origin);
	xf.xform_pos_rm(dir);
	push_back(platform(fspeed, rspeed, sdelay, rdelay, ext_dist, act_dist, origin, dir, (cont != 0), (is_rotation != 0), (update_light != 0), (destroys != 0), cur_sound_id, cur_sensor));
	cur_sound_id = -1; // reset to null after use
	if (!triggers.empty()) {back().add_triggers(triggers);} // if a custom trigger is used, reset any built-in trigger
	return 1;
}

void platform::write_to_cobj_file(std::ostream &out) const {

	triggers.write_to_cobj_file(out);
	sensor.write_to_cobj_file(out);

	if (sound_id >= 0) {
		string const sound_name(get_sound_name(sound_id));
		if (!sound_name.empty()) {out << "sound_file " << sound_name << endl;} // skip if sounds not enabled
	}
	// 'Q': // platform: enabled [fspeed rspeed sdelay rdelay ext_dist|rot_angle act_dist origin<x,y,z> dir|rot_axis<x,y,z> cont [is_rotation=0 [update_light=0]]]
	out << "Q 1 " << fspeed*TICKS_PER_SECOND << " " << rspeed*TICKS_PER_SECOND << " " << sdelay/TICKS_PER_SECOND << " " << rdelay/TICKS_PER_SECOND << " " << ext_dist << " " << act_dist
		<< " " << origin.raw_str() << " " << dir.raw_str() << " " << cont << " " << is_rotation() << " " << update_light << " " << destroys << endl; // always enabled
	
	for (auto i = lights.begin(); i != lights.end(); ++i) {
		assert(*i < light_sources_d.size());
		//light_sources_d[*i].write_to_cobj_file(out, 1);
		// FIXME: see 'L' command - need to output lights bound to this platform here
	}
	sensor.write_end_sensor_to_cobj_file(out);
	triggers.write_end_triggers_cobj_file(out);
}


void platform_cont::check_activate(point const &p, float radius, int activator) {
	for (auto i = begin(); i != end(); ++i) {i->check_activate(p, radius, activator);}
}

void platform_cont::shift_by(vector3d const &val) {
	for (auto i = begin(); i != end(); ++i) {i->shift_by(val);}
}

void platform_cont::add_current_cobjs() {
	for (auto i = begin(); i != end(); ++i) {i->clear_cobjs();}

	for (cobj_id_set_t::const_iterator i = coll_objects.platform_ids.begin(); i != coll_objects.platform_ids.end(); ++i) {
		get_cobj_platform(coll_objects.get_cobj(*i)).add_cobj(*i);
	}
}

void platform_cont::advance_timestep() {
	for (auto i = begin(); i != end(); ++i) {i->next_frame();}
	add_current_cobjs();
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

