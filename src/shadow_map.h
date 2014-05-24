// 3D World - Shared shadow map related classes
// by Frank Gennari
// 5/23/14

#ifndef _SHADOW_MAP_H_
#define _SHADOW_MAP_H_

#include "transform_obj.h" // for xform_matrix


struct smap_data_t {

	unsigned tid, tu_id, fbo_id, base_tu_id;
	pos_dir_up pdu;
	point last_lpos;
	xform_matrix texture_matrix;

	smap_data_t(unsigned base_tu_id_) : tid(0), tu_id(0), fbo_id(0), base_tu_id(base_tu_id_), last_lpos(all_zeros) {}
	virtual ~smap_data_t() {} // free_gl_state()?

	void free_gl_state() {
		free_texture(tid);
		free_fbo(fbo_id);
	}
	bool set_smap_shader_for_light(shader_t &s, int light, float z_bias) const;
	void create_shadow_map_for_light(int light, point const &lpos, cube_t const &bounds);
	virtual void render_scene_shadow_pass(point const &lpos) = 0;
	virtual bool needs_update(int light, point const &lpos);
};


#endif // _SHADOW_MAP_H_


