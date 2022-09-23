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
	bool in_building, is_stopped, is_female;

	person_base_t(float radius_) : target_pos(all_zeros), dir(zero_vector), vel(zero_vector), pos(all_zeros), radius(radius_),
		speed(0.0), anim_time(0.0), model_id(0), ssn(0), in_building(0), is_stopped(0), is_female(0) {}
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
	bool is_close_to_player() const;
};

enum {AI_STOP=0, AI_WAITING, AI_NEXT_PT, AI_BEGIN_PATH, AI_AT_DEST, AI_MOVING, AI_TO_REMOVE};
enum {GOAL_TYPE_NONE=0, GOAL_TYPE_ROOM, GOAL_TYPE_PLAYER, GOAL_TYPE_SOUND};

struct person_t : public person_base_t { // building person
	float retreat_time;
	int cur_bldg, cur_room, dest_room; // Note: -1 is unassigned
	unsigned short cur_rseed;
	uint8_t goal_type;
	bool following_player, is_on_stairs, has_key, is_first_path, on_new_path_seg;
	vector<point> path; // stored backwards, next point on path is path.back()

	person_t(float radius_) : person_base_t(radius_), retreat_time(0.0), cur_bldg(-1), cur_room(-1), dest_room(-1), cur_rseed(1), goal_type(GOAL_TYPE_NONE),
		following_player(0), is_on_stairs(0), has_key(0), is_first_path(1), on_new_path_seg(0) {in_building = 1;}
	bool on_stairs() const {return is_on_stairs;}
	void next_path_pt(bool starting_path);
};

