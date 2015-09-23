// 3D World - trigger classes for lights and platforms
// by Frank Gennari
// 5/2/15

#ifndef _TRIGGER_H_
#define _TRIGGER_H_

#include "3DWorld.h"

struct trigger_t {

	point act_pos;
	float act_dist, auto_off_time, auto_on_time;
	bool player_only, use_act_region, requires_action;
	cube_t act_region;

	trigger_t(point const &ap=all_zeros, float ad=0.0, float off_t=0.0, float on_t=0.0, bool po=0) :
	act_pos(ap), act_dist(ad), auto_off_time(off_t), auto_on_time(on_t), player_only(po), use_act_region(0), requires_action(0) {}
	void set_act_region(cube_t const ar) {act_region = ar; use_act_region = 1; act_dist = 0.0;}
	unsigned register_player_pos(point const &p, float act_radius, int activator, bool clicks=0);
	bool is_active() const {return (act_dist > 0.0 || use_act_region || auto_on_time > 0.0);}
	void shift_by(vector3d const &val) {act_pos += val; if (use_act_region) {act_region.translate(val);}}
	void write_to_cobj_file(std::ostream &out) const;
};


struct multi_trigger_t : public vector<trigger_t> {
	void add_triggers(multi_trigger_t const &t) {copy(t.begin(), t.end(), back_inserter(*this));}
	unsigned register_player_pos(point const &p, float act_radius, int activator, bool clicks=0);
	bool is_active() const {return !empty();}
	void shift_by(vector3d const &val);
	float get_auto_on_time() const;
	float get_auto_off_time() const;
	void write_to_cobj_file(std::ostream &out) const;
	void write_end_triggers_cobj_file(std::ostream &out) const;
};


#endif // _TRIGGER_H_

