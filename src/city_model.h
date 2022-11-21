// 3D World - Loading and Drawing of 3D Models for Cities
// by Frank Gennari
// 6/5/2020

#pragma once

#include "3DWorld.h"
#include "model3d.h"


struct animation_state_t {
	bool enabled;
	unsigned anim_id, model_anim_id;
	float anim_time;

	animation_state_t(bool enabled_=0, unsigned anim_id_=0, float anim_time_=0.0, unsigned model_anim_id_=0) :
		enabled(enabled_), anim_id(anim_id_), model_anim_id(model_anim_id_), anim_time(anim_time_) {}
	void set_animation_time(float anim_time_) {anim_time = anim_time_;}
	void set_animation_id_and_time(shader_t &s, bool has_bone_animations=0) const;
	void clear_animation_id(shader_t &s) const;
};


struct city_model_t {

	string fn;
	bool valid, swap_xz, swap_yz, two_sided;
	int body_mat_id, fixed_color_id, recalc_normals, centered; // recalc_normals: 0=no, 1=yes, 2=face_weight_avg
	int blade_mat_id; // for helicopters
	int model3d_id; // index into model3ds vector; -1 is not set
	float xy_rot, lod_mult, scale; // xy_rot in degrees
	colorRGBA custom_color;
	vector<unsigned> shadow_mat_ids;

	city_model_t() : valid(0), swap_xz(0), swap_yz(1), two_sided(0), body_mat_id(-1), fixed_color_id(-1), recalc_normals(1),
		centered(0), blade_mat_id(-1), model3d_id(-1), xy_rot(0.0), lod_mult(1.0), scale(1.0) {}
	city_model_t(string const &fn_, int bmid, int fcid, float rot, float dz_, float lm, vector<unsigned> const &smids) :
		fn(fn_), valid(0), swap_xz(0), swap_yz(1), two_sided(0), body_mat_id(bmid), fixed_color_id(fcid), recalc_normals(1),
		centered(0), blade_mat_id(-1), model3d_id(-1), xy_rot(rot), lod_mult(lm), scale(1.0), shadow_mat_ids(smids) {}
	bool is_loaded() const {return (model3d_id >= 0);}
	bool read(FILE *fp, bool is_helicopter=0);
	bool check_filename();
};


class city_model_loader_t : public model3ds {
protected:
	void ensure_models_loaded() {if (empty()) {load_models();}}
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
	void load_models();
	void load_model_id(unsigned id);
	void draw_model(shader_t &s, vector3d const &pos, cube_t const &obj_bcube, vector3d const &dir, colorRGBA const &color,
		vector3d const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail=0, animation_state_t *anim_state=nullptr,
		unsigned skip_mat_mask=0, bool untextured=0, bool force_high_detail=0, bool upside_down=0, bool emissive=0);
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
public:
	virtual unsigned num_models() const;
	virtual city_model_t const &get_model(unsigned id) const;
	virtual city_model_t       &get_model(unsigned id);
};

class object_model_loader_t : public city_model_loader_t {
	city_model_t null_model;
public:
	virtual unsigned num_models() const;
	virtual unsigned get_num_sub_models(unsigned id) const;
	virtual bool can_skip_model(unsigned id) const;
	virtual city_model_t const &get_model(unsigned id) const;
	virtual city_model_t       &get_model(unsigned id);
};

