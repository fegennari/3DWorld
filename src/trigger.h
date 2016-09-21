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
	unsigned register_activator_pos(point const &p, float act_radius, int activator, bool clicks=0) const;
	bool is_active() const {return (act_dist > 0.0 || use_act_region || auto_on_time > 0.0);}
	void shift_by(vector3d const &val) {act_pos += val; if (use_act_region) {act_region.translate(val);}}
	void write_to_cobj_file(std::ostream &out) const;
};


struct multi_trigger_t : public vector<trigger_t> {
	void add_triggers(multi_trigger_t const &t) {copy(t.begin(), t.end(), back_inserter(*this));}
	unsigned register_activator_pos(point const &p, float act_radius, int activator, bool clicks=0);
	bool is_active() const {return !empty();}
	void shift_by(vector3d const &val);
	float get_auto_on_time() const;
	float get_auto_off_time() const;
	void write_to_cobj_file(std::ostream &out) const;
	void write_end_triggers_cobj_file(std::ostream &out) const;
};


enum {SENSOR_DISABLED=0, SENSOR_ALWAYS_OFF, SENSOR_ALWAYS_ON, SENSOR_LIGHT, SENSOR_SOUND, SENSOR_HEAT, SENSOR_METAL, SENSOR_WATER, SENSOR_PRESSURE, SENSOR_SMOKE, NUM_SENSOR_TYPES};

struct geom_xform_t;

class sensor_t {
	point pos;
	float radius; // only used by some sensor types (sound, metal, pressure)
	float thresh; // only used by some sensor types (light?, sound?, heat, smoke)
	int type;
	bool invert;

	bool check_active_int() const;

public:
	sensor_t() : pos(all_zeros), radius(0.0), thresh(0.0), type(SENSOR_DISABLED), invert(0) {}
	sensor_t(point const &pos_, int type_, bool invert_=0, float radius_=0.0, float thresh_=0.0) : pos(pos_), radius(radius_), thresh(thresh_), type(type_), invert(invert_) {
		assert(type >= SENSOR_DISABLED && type < NUM_SENSOR_TYPES);
	}
	bool enabled() const {return (type != SENSOR_DISABLED);}
	bool check_active() const {return (check_active_int() ^ invert);}
	bool read_from_file(FILE *fp, geom_xform_t const &xf);
	void write_to_cobj_file(std::ostream &out) const;
};

struct multi_sensor_t : public vector<sensor_t> {
	// WRITE
};


#endif // _TRIGGER_H_

