// 3D World - Loading and Drawing of 3D Models for Cities
// by Frank Gennari
// 6/5/2020

#pragma once

#include "3DWorld.h"
#include "model3d.h"


struct city_model_t {

	string fn;
	bool valid, loaded, swap_xz, swap_yz;
	int body_mat_id, fixed_color_id, recalc_normals; // recalc_normals: 0=no, 1=yes, 2=face_weight_avg
	int blade_mat_id; // for helicopters
	float xy_rot, lod_mult, scale; // xy_rot in degrees
	vector<unsigned> shadow_mat_ids;

	city_model_t() : valid(0), loaded(0), swap_xz(0), swap_yz(1), body_mat_id(-1), fixed_color_id(-1), recalc_normals(1),
		blade_mat_id(-1), xy_rot(0.0), lod_mult(1.0), scale(1.0) {}
	city_model_t(string const &fn_, int bmid, int fcid, float rot, float dz_, float lm, vector<unsigned> const &smids) :
		fn(fn_), valid(0), loaded(0), swap_xz(0), swap_yz(1), body_mat_id(bmid), fixed_color_id(fcid), recalc_normals(1),
		blade_mat_id(-1), xy_rot(rot), lod_mult(lm), scale(1.0), shadow_mat_ids(smids) {}
	bool read(FILE *fp, bool is_helicopter=0);
	bool check_filename();
};


class city_model_loader_t : public model3ds {
protected:
	void ensure_models_loaded() {if (empty()) {load_models();}}
public:
	virtual ~city_model_loader_t() {}
	virtual bool has_low_poly_model() {return 0;}
	virtual unsigned num_models() const = 0;
	virtual city_model_t const &get_model(unsigned id) const = 0;
	virtual city_model_t &get_model(unsigned id) = 0;
	vector3d get_model_world_space_size(unsigned id);
	colorRGBA get_avg_color(unsigned id);
	bool is_model_valid(unsigned id);
	void load_models();
	bool load_model_id(unsigned id);
	void draw_model(shader_t &s, vector3d const &pos, cube_t const &obj_bcube, vector3d const &dir, colorRGBA const &color,
		vector3d const &xlate, unsigned model_id, bool is_shadow_pass, bool low_detail=0, bool enable_animations=0, unsigned skip_mat_mask=0);
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
	virtual city_model_t const &get_model(unsigned id) const;
	virtual city_model_t       &get_model(unsigned id);
};

