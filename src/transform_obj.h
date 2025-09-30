// 3D World
// by Frank Gennari
// object transformation/deformation class definitions
// 6/11/06
#pragma once

#include "3DWorld.h"
#include "mesh2d.h"

#define GLM_FORCE_RADIANS
#include <glm/mat4x4.hpp>
#include <glm/vec3.hpp>


inline glm::vec3 vec3_from_vector3d (vector3d const &v) {return glm::vec3(v.x, v.y, v.z);}
inline vector3d  vector3d_from_vec3(glm::vec3 const &v) {return vector3d (v.x, v.y, v.z);}

void print_matrix(float const *const m, std::string const &prefix=std::string(), std::ostream &out=std::cout);
glm::mat4 get_rotation_matrix(vector3d const &vrot, float angle);


struct xform_matrix : public glm::mat4 { // Note: maybe better to use glm::gtx::simd_mat4?

	xform_matrix() : glm::mat4(1.0) {} // identity
	xform_matrix(glm::mat4 const &m) : glm::mat4(m) {}
	xform_matrix inverse() const; // defined in shaders.cpp
	float *get_ptr();
	float const *get_ptr() const;
	void get_as_doubles(double md[16]) const;
	void normalize();
	void check_valid(const char *msg_str) const;
	void apply_to_vector3d(vector3d &v) const;
	void apply_rotate(float angle_degrees, float x, float y, float z);
	void print(std::string const &prefix=std::string(), std::ostream &out=std::cout) const {print_matrix(get_ptr(), prefix, out);}
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
};


class matrix_stack_t {

	vector<xform_matrix> m;

public:
	matrix_stack_t() {m.push_back(xform_matrix());} // will be identity
	void push() {assert(!m.empty()); m.push_back(m.back());} // duplicate top element
	void pop()  {assert(!m.empty()); m.pop_back(); assert(!m.empty());} // can't start or end as an empty matrix stack
	void push_identity() {m.push_back(xform_matrix());}
	xform_matrix const &top() const {assert(!m.empty()); return m.back();}
	void assign(xform_matrix const &v) {assert(!m.empty()); m.back() = v;}
	void identity() {assert(!m.empty()); m.back() = glm::mat4(1.0);}
};


// function prototypes
xform_matrix const &fgGetMVM();
xform_matrix const &fgGetPJM();
void apply_roll_to_matrix(xform_matrix &matrix, point const &pos, point const &lpos, vector3d const &ground_normal, float radius, float a_add=0.0, float a_mult=1.0);
void apply_obj_mesh_roll(xform_matrix &matrix, point const &pos, point const &lpos, float radius, float a_add=0.0, float a_mult=1.0);


class instance_render_t {

	vector<xform_matrix> inst_xforms;
	int loc;

public:
	instance_render_t(int loc_=-1) : loc(loc_) {}
	void set_loc(int loc_) {loc = loc_;}
	void add_cur_inst() {add_inst(fgGetMVM());}
	void add_inst(xform_matrix const &xf) {inst_xforms.push_back(xf);}
	void draw_and_clear(int prim_type, unsigned count, unsigned cur_vbo=0, int index_type=GL_NONE, void const *const indices=NULL, unsigned first=0, unsigned cur_vao=0);
	unsigned size () const {return inst_xforms.size();}
	bool     empty() const {return inst_xforms.empty();}
};

struct bone_transform_data_t {
	int anim_id=-2; // -2 is unset/invalid
	float anim_time=0.0;
	vector<xform_matrix> transforms;
	void clear() {anim_id = -2; anim_time = 0.0; clear_container(transforms);} // and free the memory
};

