// 3D World
// by Frank Gennari
// object transformation/deformation class definitions
// 6/11/06
#ifndef _TRANSFORM_OBJ_H_
#define _TRANSFORM_OBJ_H_

#include "3DWorld.h"
#include "mesh2d.h"
#include <glm/mat4x4.hpp>


class xform_matrix {

protected:
	float m[16];

public:
	xform_matrix() {load_identity();}
	xform_matrix(float const *const m_) {assign(m_);}
	void assign(float const *const m_) {for(unsigned i = 0; i < 16; ++i) {m[i] = m_[i];}}
	void assign(float v0, float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8, float v9, float v10, float v11, float v12, float v13, float v14, float v15) {
		m[0] = v0; m[1] = v1; m[2] = v2; m[3] = v3; m[4] = v4; m[5] = v5; m[6] = v6; m[7] = v7; m[8] = v8; m[9] = v9; m[10] = v10; m[11] = v11; m[12] = v12; m[13] = v13; m[14] = v14; m[15] = v15;
	}
	void apply() const;
	void assign_mv_from_gl();
	void assign_pj_from_gl();
	void normalize();
	void load_identity();
	void rotate(float angle, vector3d const &rot);
	float *get_ptr() {return m;}
};


struct xform_matrix_glm : public glm::mat4 {

	xform_matrix_glm() {}
	xform_matrix_glm(glm::mat4 const &m) : glm::mat4(m) {}
	void apply() const;
	void assign_mv_from_gl();
	void assign_pj_from_gl();
	float *get_ptr();
	float const *get_ptr() const;
};


struct transform_data {

	vector<xform_matrix> matrices;
	vector<mesh2d> perturb_maps;

	void resize(unsigned sz) {
		matrices.resize(sz);
		perturb_maps.resize(sz);
	}
	xform_matrix &get_matrix(unsigned i) {
		assert(i < matrices.size());
		return matrices[i];
	}
	mesh2d &get_mesh(unsigned i) {
		assert(i < perturb_maps.size());
		return perturb_maps[i];
	}
	void set_perturb_size(unsigned i, unsigned sz);
	void add_perturb_at(unsigned s, unsigned t, unsigned i, float val, float min_mag, float max_mag);
	void reset_perturb_if_set(unsigned i);
	~transform_data();
};


class instance_render_t { // is this a base class of shader_t?

	vector<xform_matrix> inst_xforms;
	int loc;

public:
	instance_render_t(int loc_=-1) : loc(loc_) {}
	void set_loc(int loc_) {loc = loc_;}
	void add_cur_inst();
	void add_inst(xform_matrix const &xf) {inst_xforms.push_back(xf);}
	void draw_and_clear(int prim_type, unsigned count, unsigned cur_vbo=0, int index_type=GL_NONE, void *indices=NULL);
	unsigned size () const {return inst_xforms.size();}
	bool     empty() const {return inst_xforms.empty();}
};


// function prototypes
void apply_obj_mesh_roll(xform_matrix &matrix, point const &pos, point const &lpos, float radius, float a_add=0.0, float a_mult=1.0);


#endif

