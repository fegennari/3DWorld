// 3D World
// by Frank Gennari
// object transformation/deformation class definitions
// 6/11/06
#ifndef _TRANSFORM_OBJ_H_
#define _TRANSFORM_OBJ_H_

#include "3DWorld.h"


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
	void normalize();
	void apply() const {glMultMatrixf(m);}
	void assign_mv_from_gl() {glGetFloatv(GL_MODELVIEW_MATRIX,  m);}
	void assign_pj_from_gl() {glGetFloatv(GL_PROJECTION_MATRIX, m);}
	void load_identity();
	void rotate(float angle, vector3d const &rot);
	void translate(vector3d const &t);
	void scale(vector3d const &s);
	void scale(float s) {scale(vector3d(s,s,s));}
	float *get_ptr() {return m;}
};


struct mesh2d {

	float *pmap; // perturbation map
	bool *rmap;  // render map
	float *emap; // expand map
	point *ptsh; // point shift map
	unsigned size;

	unsigned get_index(unsigned s, unsigned t) const {assert(s < size && t <= size); return (s*(size+1) + t);}

public:
	float expand;

	mesh2d() : pmap(NULL), rmap(NULL), emap(NULL), ptsh(NULL), size(0), expand(0.0) {}
	//~mesh2d() {clear();}
	void clear();
	unsigned get_num()     const {assert(size > 0); return size*(size+1);} // square, with an extra row
	unsigned choose_rand() const {return (rand() % get_num());}
	void set_size(unsigned sz);
	template<typename T> void alloc_ptr(T *&p, T const val);
	void alloc_pmap();
	void alloc_rmap();
	void alloc_emap();
	void alloc_ptsh();
	void reset_pmap();
	void add_random(float mag, float min_mag, float max_mag, unsigned skipval=0);
	void mult_by(float val);
	void unset_rand_rmap(unsigned num_remove);
	void set_rand_expand(float mag, unsigned num_exp);
	void set_rand_translate(point const &tp, unsigned num_trans);
	void set_val(unsigned s, unsigned t, float val)      {assert(pmap); pmap[get_index(s, t)] = val;}
	float get_val(unsigned s, unsigned t) const          {assert(pmap); return pmap[get_index(s, t)];}
	void  set_rm(unsigned s, unsigned t, bool val)       {assert(rmap); rmap[get_index(s, t)] = val;}
	bool  get_rm(unsigned s, unsigned t) const           {assert(rmap); return rmap[get_index(s, t)];}
	void  set_em(unsigned s, unsigned t, float val)      {assert(emap); emap[get_index(s, t)] = val;}
	float get_em(unsigned s, unsigned t) const           {assert(emap); return emap[get_index(s, t)];}
	void  set_pt(unsigned s, unsigned t, point const &p) {assert(ptsh); ptsh[get_index(s, t)] = p;}
	point get_pt(unsigned s, unsigned t) const           {assert(ptsh); return ptsh[get_index(s, t)];}
	unsigned get_size() const {return size;}
	void draw_perturbed_sphere(point const &pos, float radius, int ndiv, bool tex_coord) const;
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
	void add_rand_perturb(unsigned i, float mag, float min_mag, float max_mag);
	void add_perturb_at(unsigned s, unsigned t, unsigned i, float val, float min_mag, float max_mag);
	void reset_perturb_if_set(unsigned i);
	~transform_data();
};


// function prototypes
void apply_obj_mesh_roll(xform_matrix &matrix, point const &pos, point const &lpos, float radius, float a_add=0.0, float a_mult=1.0);
void deform_obj(dwobject &obj, vector3d const &norm, vector3d const &v0);
void update_deformation(dwobject &obj);


#endif

