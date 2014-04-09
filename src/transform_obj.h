// 3D World
// by Frank Gennari
// object transformation/deformation class definitions
// 6/11/06
#ifndef _TRANSFORM_OBJ_H_
#define _TRANSFORM_OBJ_H_

#include "3DWorld.h"
#include "mesh2d.h"

#define GLM_FORCE_RADIANS
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>


inline glm::vec3 vec3_from_vector3d(vector3d const &v) {return glm::vec3(v.x, v.y, v.z);}


struct xform_matrix : public glm::mat4 {

	xform_matrix() {}
	xform_matrix(glm::mat4 const &m) : glm::mat4(m) {}
	void assign_mv_from_gl();
	void assign_pj_from_gl();
	float *get_ptr();
	float const *get_ptr() const;
	void normalize();
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


class matrix_stack_t {

	vector<xform_matrix> m;

public:
	matrix_stack_t() {m.push_back(xform_matrix());} // will be identity
	void push() {assert(!m.empty()); m.push_back(m.back());} // duplicate top element
	void pop()  {assert(!m.empty()); m.pop_back(); assert(!m.empty());} // can't start or end as an empty matrix stack
	xform_matrix const &top() const {assert(!m.empty()); return m.back();}
	void assign(xform_matrix const &v) {assert(!m.empty()); m.back() = v;}
	void identity() {assert(!m.empty()); m.back() = xform_matrix();}
};


class instance_render_t {

	vector<xform_matrix> inst_xforms;
	int loc;

public:
	instance_render_t(int loc_=-1) : loc(loc_) {}
	void set_loc(int loc_) {loc = loc_;}
	void add_cur_inst();
	void add_inst(xform_matrix const &xf) {inst_xforms.push_back(xf);}
	void draw_and_clear(int prim_type, unsigned count, unsigned cur_vbo=0, int index_type=GL_NONE, void *indices=NULL, unsigned first=0);
	unsigned size () const {return inst_xforms.size();}
	bool     empty() const {return inst_xforms.empty();}
};


// function prototypes
void apply_obj_mesh_roll(xform_matrix &matrix, shader_t &shader, point const &pos, point const &lpos, float radius, float a_add=0.0, float a_mult=1.0);


#endif

