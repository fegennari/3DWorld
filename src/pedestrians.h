// 3D World - Buildings Interface
// by Frank Gennari
// 4-5-18
#pragma once

#include "3DWorld.h"

float const PED_WIDTH_SCALE  = 0.5; // ratio of collision radius to model radius (x/y)
float const PED_HEIGHT_SCALE = 2.5; // ratio of collision radius to model height (z)

extern double tfticks;

struct waiting_obj_t {
	float waiting_start;
	waiting_obj_t() : waiting_start(0.0) {}
	void reset_waiting() {waiting_start = tfticks;}
	float get_wait_time_secs() const {return (float(tfticks) - waiting_start)/TICKS_PER_SECOND;} // Note: only meaningful for cars stopped at lights or peds stopped at roads
};

struct person_base_t : public waiting_obj_t {

	point target_pos;
	vector3d dir, vel;
	point pos;
	float radius, speed, anim_time;
	unsigned short model_id, ssn;
	bool in_building, is_stopped, destroyed, is_female;

	person_base_t(float radius_) : target_pos(all_zeros), dir(zero_vector), vel(zero_vector), pos(all_zeros), radius(radius_),
		speed(0.0), anim_time(0.0), model_id(0), ssn(0), in_building(0), is_stopped(0), destroyed(0), is_female(0) {}
	std::string get_name() const;
	unsigned get_unique_id() const {return ssn;} // technically only unique if there are <= 65536 people
	float get_height () const {return PED_HEIGHT_SCALE*radius;}
	float get_width  () const {return PED_WIDTH_SCALE *radius;}
	float get_z1     () const {return (pos.z - radius);}
	float get_z2     () const {return (get_z1() + get_height());}
	cube_t get_bcube () const;
	point get_eye_pos() const {return (pos + vector3d(0.0, 0.0, (0.9*get_height() - radius)));}
	bool target_valid() const {return (target_pos != all_zeros);}
	bool is_waiting_or_stopped() const {return (speed == 0.0 || waiting_start > 0);}
	void set_velocity(vector3d const &v) {vel = v*(speed/v.mag());} // normalize to original velocity
	void stop();
	void go();
	void wait_for(float seconds);
	void destroy() {destroyed = 1;} // that's it, no other effects
	bool is_close_to_player() const;
};

struct person_t : public person_base_t { // building person
	float retreat_time;
	int cur_bldg;
	unsigned short cur_rseed;
	bool following_player, is_on_stairs, has_key;

	person_t(float radius_) : person_base_t(radius_), retreat_time(0.0), cur_bldg(-1), cur_rseed(1), following_player(0), is_on_stairs(0), has_key(0) {in_building = 1;}
	bool on_stairs() const {return is_on_stairs;}
};

