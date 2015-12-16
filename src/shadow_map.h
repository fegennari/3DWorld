// 3D World - Shared shadow map related classes
// by Frank Gennari
// 5/23/14

#ifndef _SHADOW_MAP_H_
#define _SHADOW_MAP_H_

#include "transform_obj.h" // for xform_matrix

unsigned const DEF_LOCAL_SMAP_SZ      = 1024;
unsigned const LOCAL_SMAP_START_TU_ID = 16;
unsigned const MAX_DLIGHT_SMAPS       = 16;


class smap_texture_array_t {

	unsigned num_layers, num_layers_used;

public:
	unsigned tid, gen_id;

	smap_texture_array_t() : num_layers(0), num_layers_used(0), tid(0), gen_id(1) {} // gen_id starts at 1
	bool is_allocated() const {return (tid > 0);}
	void free_gl_state();
	void ensure_tid(unsigned xsize, unsigned ysize);
	void reserve_num_layers(unsigned num);
	unsigned new_layer();
	void clear() {num_layers = num_layers_used = 0; free_gl_state();}
};

class smap_data_state_t {

protected:
	smap_texture_array_t *tex_arr;
	unsigned fbo_id, local_tid, gen_id, layer_id;

public:
	smap_data_state_t() : tex_arr(nullptr), fbo_id(0), local_tid(0), gen_id(0), layer_id(0) {}
	bool is_arrayed() const {return (tex_arr != nullptr);}
	void bind_tex_array(smap_texture_array_t *tex_arr_);
	unsigned get_tid() const {return (is_arrayed() ? tex_arr->tid : local_tid);}
	bool is_allocated() const {return (get_tid() > 0 && (!is_arrayed() || gen_id == tex_arr->gen_id));}
	unsigned *get_layer() {return (tex_arr ? &layer_id : nullptr);}
	void free_gl_state();
	void disown() {assert(!is_arrayed()); local_tid = gen_id = fbo_id = 0;}
};

struct smap_data_t : public smap_data_state_t {

	unsigned tu_id, smap_sz;
	pos_dir_up pdu;
	point last_lpos;
	xform_matrix texture_matrix;

	smap_data_t(unsigned tu_id_, unsigned smap_sz_, smap_data_state_t const &init_state=smap_data_state_t())
		: smap_data_state_t(init_state), tu_id(tu_id_), smap_sz(smap_sz_), last_lpos(all_zeros) {}
	virtual ~smap_data_t() {} // free_gl_state()?
	bool set_smap_shader_for_light(shader_t &s, int light, xform_matrix const *const mvm=nullptr) const;
	void bind_smap_texture(bool light_valid=1) const;
	void create_shadow_map_for_light(point const &lpos, cube_t const *const bounds=nullptr, bool use_world_space=0);
	virtual void render_scene_shadow_pass(point const &lpos) = 0;
	virtual bool needs_update(point const &lpos);
	virtual bool is_local() const {return 0;} // for debugging only
};

struct cached_dynamic_smap_data_t : public smap_data_t {

	bool last_has_dynamic;
	cached_dynamic_smap_data_t(unsigned tu_id_, unsigned smap_sz_) : smap_data_t(tu_id_, smap_sz_), last_has_dynamic(0) {}
};

struct local_smap_data_t : public cached_dynamic_smap_data_t {

	bool used;

	local_smap_data_t(unsigned tu_id_, unsigned smap_sz_=DEF_LOCAL_SMAP_SZ) : cached_dynamic_smap_data_t(tu_id_, smap_sz_), used(0) {}
	bool set_smap_shader_for_light(shader_t &s, bool &arr_tex_set) const;
	virtual void render_scene_shadow_pass(point const &lpos);
	virtual bool needs_update(point const &lpos);
	virtual bool is_local() const {return 1;} // for debugging only
};

struct local_cube_map_smap_data_t : public local_smap_data_t { // to be implemented/used later
	local_cube_map_smap_data_t(unsigned tu_id_, unsigned smap_sz_=DEF_LOCAL_SMAP_SZ) : local_smap_data_t(tu_id_, smap_sz_) {}
};


template<class SD> struct vect_smap_t : public vector<SD> { // one per light source (sun, moon)
	void set_for_all_lights(shader_t &s, xform_matrix const *const mvm) const {
		for (unsigned i = 0; i < size(); ++i) {operator[](i).set_smap_shader_for_light(s, i, mvm);}
	}
	void clear() {
		for (iterator i = begin(); i != end(); ++i) {i->free_gl_state();}
		vector<SD>::clear();
	}
	void create_if_needed(cube_t const &bcube) {
		for (unsigned i = 0; i < size(); ++i) {
			point lpos;
			if (!light_valid_and_enabled(i, lpos)) continue;
			operator[](i).create_shadow_map_for_light(lpos, &bcube);
		}
	}
};


#endif // _SHADOW_MAP_H_


