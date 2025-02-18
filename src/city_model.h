// 3D World - Loading and Drawing of 3D Models for Cities
// by Frank Gennari
// 6/5/2020

#pragma once

#include "3DWorld.h"
#include "model3d.h"


float const SKELETAL_ANIM_TIME_CONST = 32.0; // maybe should rework these constants; this should be for people only, not all models

enum {ANIM_ID_NONE=0, ANIM_ID_WALK, ANIM_ID_BHOP, ANIM_ID_FLIP, ANIM_ID_TWIRL, ANIM_ID_MARCH, ANIM_ID_ALIEN,
	  ANIM_ID_RAT, ANIM_ID_SPIDER, ANIM_ID_SKELETAL, ANIM_ID_HCOPTER, ANIM_ID_SWINGS, ANIM_ID_FISH_TAIL, ANIM_ID_FLAP_ARMS};
// Note: MODEL_ANIM_STAIRS is for future use when walking up and down stairs; MODEL_ANIM_ATTACK is for future use with zombies
enum {MODEL_ANIM_WALK=0, MODEL_ANIM_IDLE, MODEL_ANIM_STAIRS_UP, MODEL_ANIM_STAIRS_DOWN, MODEL_ANIM_ATTACK, MODEL_ANIM_CROUCH, NUM_MODEL_ANIMS};
string const animation_names[NUM_MODEL_ANIMS] = {"walking", "idle", "stairs_up", "stairs_down", "attack", "crouch"};

struct animation_state_t {
	bool enabled=0, fixed_anim_speed=0, cache_animations=0;
	unsigned anim_id=0, model_anim_id=0, model_anim_id2=0;
	float anim_time=0.0, anim_time2=0.0, blend_factor=0.0;
	bone_transform_data_t *cached=nullptr;

	animation_state_t(bool enabled_=0, unsigned anim_id_=ANIM_ID_NONE, float anim_time_=0.0, unsigned model_anim_id_=MODEL_ANIM_WALK) :
		enabled(enabled_), anim_id(anim_id_), model_anim_id(model_anim_id_), model_anim_id2(model_anim_id_), anim_time(anim_time_) {}
	int get_anim_id_for_setup_bone_transforms () const;
	int get_anim_id2_for_setup_bone_transforms() const;
	void set_animation_id_and_time(shader_t &s, bool has_bone_animations=0, float anim_speed=1.0) const;
	void clear_animation_id(shader_t &s) const;
};


struct city_model_t {

	string fn;
	bool valid=0, swap_xz=0, swap_yz=1, two_sided=0, mirrored=0, is_zombie=0, tried_to_load=0;
	int body_mat_id=-1, fixed_color_id=-1, recalc_normals=1, centered=0; // recalc_normals: 0=no, 1=yes, 2=face_weight_avg
	int blade_mat_id=-1; // for helicopters
	int model3d_id=-1; // index into model3ds vector; -1 is not set
	uint64_t rev_winding_mask=0;
	float xy_rot=0.0, lod_mult=1.0, scale=1.0, model_anim_scale=1.0, anim_speed=1.0; // xy_rot in degrees
	point rotate_about;
	colorRGBA custom_color;
	string default_anim_name;
	vector<unsigned> shadow_mat_ids;

	struct model_anim_t {
		string fn, anim_name;
		model_anim_t() {}
		model_anim_t(string const &fn_, string const &name) : fn(fn_), anim_name(name) {}
	};
	vector<model_anim_t> anim_fns; // for ped models, etc.

	city_model_t() {}
	city_model_t(string const &fn_, int bmid, int fcid, float rot, float dz_, float lm, vector<unsigned> const &smids) :
		fn(fn_), body_mat_id(bmid), fixed_color_id(fcid), xy_rot(rot), lod_mult(lm), shadow_mat_ids(smids) {}
	bool is_loaded() const {return (model3d_id >= 0);}
	bool read(FILE *fp, bool is_helicopter=0, bool is_person=0);
	bool check_filename();
	bool has_animation(string const &anim_name) const;
};


class city_model_loader_t : public model3ds {
protected:
	model3d &get_model3d(unsigned id);
public:
	virtual ~city_model_loader_t() {}
	virtual bool has_low_poly_model() {return 0;}
	virtual bool can_skip_model(unsigned id) const {return 0;}
	virtual unsigned num_models() const = 0;
	virtual unsigned get_num_sub_models(unsigned id) const {return 1;}
	virtual city_model_t const &get_model(unsigned id) const = 0;
	virtual city_model_t &get_model(unsigned id) = 0;
	float get_model_scale(unsigned id) const {return get_model(id).scale;}
	vector3d get_model_world_space_size(unsigned id);
	colorRGBA get_avg_color(unsigned id, bool area_weighted=1);
	bool model_filename_contains(unsigned id, string const &str, string const &str2="") const;
	bool is_model_valid(unsigned id);
	void load_model_id(unsigned id);
	bool check_anim_wrapped(unsigned model_id, unsigned model_anim_id, float old_time, float new_time);
	float get_anim_duration(unsigned model_id, unsigned model_anim_id);
	void draw_model(shader_t &s, vector3d const &pos, cube_t const &obj_bcube, vector3d const &dir, colorRGBA const &color, vector3d const &xlate,
		unsigned model_id, bool is_shadow_pass=0, bool low_detail=0, animation_state_t *anim_state=nullptr, unsigned skip_mat_mask=0,
		bool untextured=0, bool force_high_detail=0, bool upside_down=0, bool emissive=0, bool do_local_rotate=0, int mirror_dim=3);
	static void rotate_model_from_plus_x_to_dir(vector3d const &dir);
};

class car_model_loader_t : public city_model_loader_t {
public:
	virtual bool has_low_poly_model() {return 1;}
	virtual unsigned num_models() const;
	virtual city_model_t const &get_model(unsigned id) const;
	virtual city_model_t       &get_model(unsigned id);
};

class helicopter_model_loader_t : public city_model_loader_t {
public:
	virtual unsigned num_models() const;
	virtual city_model_t const &get_model(unsigned id) const;
	virtual city_model_t       &get_model(unsigned id);
};

class ped_model_loader_t : public city_model_loader_t {
	vector<unsigned> people_models, zombie_models;
public:
	virtual unsigned num_models() const;
	virtual city_model_t const &get_model(unsigned id) const;
	virtual city_model_t       &get_model(unsigned id);
	int select_random_model(int rand_val, bool choose_zombie, unsigned pref_gender);
	bool has_mix_of_model_types() const {return (!people_models.empty() && !zombie_models.empty());} // Note: must call select_random_model() first
};

class object_model_loader_t : public city_model_loader_t {
	city_model_t null_model;
	int get_valid_sub_model_id(unsigned id, vector<city_model_t> const &models) const;
public:
	virtual unsigned num_models() const;
	virtual unsigned get_num_sub_models(unsigned id) const;
	virtual bool can_skip_model(unsigned id) const;
	virtual city_model_t const &get_model(unsigned id) const;
	virtual city_model_t       &get_model(unsigned id);
};

