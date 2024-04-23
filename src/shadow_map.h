// 3D World - Shared shadow map related classes
// by Frank Gennari
// 5/23/14
#pragma once

#include "transform_obj.h" // for xform_matrix

unsigned const DEF_LOCAL_SMAP_SZ      = 1024;
unsigned const LOCAL_SMAP_START_TU_ID = 16;
unsigned const GLOBAL_SMAP_START_TU_ID= 6; // for ground mode and tiled terrain mode
unsigned const MAX_DLIGHT_SMAPS       = 120; // must agree with the value used in dynamic_lighting.part


class smap_texture_array_t {

	unsigned num_layers=0, num_layers_used=0;
public:
	unsigned tid=0, gen_id=1, gpu_mem=0; // gen_id starts at 1

	bool is_allocated() const {return (tid > 0);}
	void free_gl_state();
	void ensure_tid(unsigned xsize, unsigned ysize);
	void reserve_num_layers(unsigned num);
	unsigned new_layer();
	void clear() {num_layers = num_layers_used = gpu_mem = 0; free_gl_state();}
};

class smap_data_state_t {

protected:
	smap_texture_array_t *tex_arr=nullptr;
	unsigned fbo_id=0, local_tid=0, gen_id=0, layer_id=0;
public:
	bool is_csm=0;

	bool is_arrayed() const {return (tex_arr != nullptr);}
	void bind_tex_array(smap_texture_array_t *tex_arr_);
	unsigned get_tid   () const {return (is_arrayed() ? tex_arr->tid : local_tid);}
	bool is_allocated  () const {return (get_tid() > 0 && (!is_arrayed() || gen_id == tex_arr->gen_id));}
	unsigned *get_layer() {return (tex_arr ? &layer_id : nullptr);}
	void free_gl_state ();
	void disown() {assert(!is_arrayed()); local_tid = gen_id = fbo_id = 0;}
};

struct smap_data_t : public smap_data_state_t { // used for all types of lights: ground mode, tiled terrain, and model3d

	unsigned tu_id, smap_sz;
	pos_dir_up pdu;
	point last_lpos;
	xform_matrix texture_matrix;
	vector<xform_matrix> cascade_matrices; // light space matrices, for CSMs

	smap_data_t(unsigned tu_id_, unsigned smap_sz_, smap_data_state_t const &init_state=smap_data_state_t())
	  : smap_data_state_t(init_state), tu_id(tu_id_), smap_sz(smap_sz_), last_lpos(all_zeros), texture_matrix(glm::mat4(1.0)) {}
	virtual ~smap_data_t() {} // free_gl_state()?
	bool set_smap_shader_for_light(shader_t &s, int light, xform_matrix const *const mvm=nullptr) const;
	bool bind_smap_texture(bool light_valid=1) const;
	void set_csm_matrices(shader_t &s) const;
	void create_shadow_map_for_light(point const &lpos, cube_t const *const bounds=nullptr, bool use_world_space=0, bool no_update=0, bool force_update=0);
	unsigned get_gpu_mem() const;
	virtual void render_scene_shadow_pass(point const &lpos) = 0;
	virtual bool needs_update(point const &lpos);
	virtual bool is_local() const {return 0;} // for debugging only
};

struct cached_dynamic_smap_data_t : public smap_data_t { // used for all types of lights (ground mode directional sun/moon, point, spotlight)

	bool last_has_dynamic;
	cached_dynamic_smap_data_t(unsigned tu_id_, unsigned smap_sz_) : smap_data_t(tu_id_, smap_sz_), last_has_dynamic(0) {}
};

struct local_smap_data_t : public cached_dynamic_smap_data_t { // for point/spot lights that may be dynamic; CSMs not supported

	bool used, outdoor_shadows;
	unsigned user_smap_id;

	local_smap_data_t(unsigned tu_id_, unsigned smap_sz_=DEF_LOCAL_SMAP_SZ, bool outdoor_shadows_=0)
		: cached_dynamic_smap_data_t(tu_id_, smap_sz_), used(0), outdoor_shadows(outdoor_shadows_), user_smap_id(0) {}
	bool set_smap_shader_for_light(shader_t &s, bool &arr_tex_set) const;
	virtual void render_scene_shadow_pass(point const &lpos);
	virtual bool needs_update(point const &lpos);
	virtual bool is_local() const {return 1;} // for debugging only
};

struct local_cube_map_smap_data_t : public local_smap_data_t { // unused; to be implemented/used later
	local_cube_map_smap_data_t(unsigned tu_id_, unsigned smap_sz_=DEF_LOCAL_SMAP_SZ) : local_smap_data_t(tu_id_, smap_sz_) {}
};


struct rotation_t {
	vector3d axis;
	float angle; // in degrees

	rotation_t() : axis(zero_vector), angle(0.0) {}
	rotation_t(vector3d const &axis_, float angle_) : axis(axis_), angle(angle_) {}
	bool operator==(rotation_t const &r) const {return (r.axis == axis && r.angle == angle);}
	bool operator< (rotation_t const &r) const {return ((angle == r.angle) ? (axis < r.axis) : (angle < r.angle));}
	void apply_gl() const {if (angle != 0.0) {rotate_about(angle, axis);}}
	void rotate_point(point &pos, float sign) const {
		if (angle != 0.0) {rotate_vector3d(axis, sign*double(TO_RADIANS)*angle, pos);}
	}
};


template<class SD> struct vect_smap_t : public vector<SD> { // one per light source (sun, moon); used by tiled terrain
	void set_for_all_lights(shader_t &s, xform_matrix const *const mvm) const {
		for (unsigned i = 0; i < vector<SD>::size(); ++i) {vector<SD>::operator[](i).set_smap_shader_for_light(s, i, mvm);}
	}
	void clear() {
		for (auto i = vector<SD>::begin(); i != vector<SD>::end(); ++i) {i->free_gl_state();}
		vector<SD>::clear();
	}
	void create_if_needed(cube_t const &bcube, rotation_t const *const inv_rot=nullptr) {
		for (unsigned i = 0; i < vector<SD>::size(); ++i) {
			point lpos;
			if (!light_valid_and_enabled(i, lpos)) continue;
			if (inv_rot != nullptr) {inv_rot->rotate_point(lpos, 1.0);} // inverse rotation from shadow casting object
			vector<SD>::operator[](i).create_shadow_map_for_light(lpos, &bcube);
		}
	}
};

unsigned get_empty_smap_tid(bool is_csm=0);
void bind_default_sun_moon_smap_textures();

