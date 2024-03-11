// 3D World
// by Frank Gennari
// 1/8/06
#pragma once

#include "3DWorld.h"
#include "gl_ext_arb.h"

class sd_sphere_d; // forward reference
class sd_sphere_vbo_d; // forward reference
class upsurface;
class coll_obj_group;
class instance_render_t;


class sphere_point_norm { // size = 12

protected:
	unsigned ndiv;
	bool is_full;
	float radius;
	point center;
	point **points;
	vector3d **norms;

public:
	friend class sd_sphere_d;
	sphere_point_norm() : ndiv(0), is_full(0), radius(0.0), points(NULL), norms(NULL) {}
	void alloc(unsigned ndiv_);
	void set_pointer_stride(unsigned ndiv_);
	void free_data();
	point    **get_points() const {return points;}
	vector3d **get_norms () const {return norms ;}
	unsigned get_ndiv    () const {return ndiv  ;}
};


class sd_sphere_d : public sphere_point_norm { // size = 40

protected:
	point pos;
	float radius, def_pert;
	float const *perturb_map;
	upsurface const *surf;

public:
	typedef vert_norm_tc vertex_type_t;
	//typedef vert_norm_comp_tc vertex_type_t; // lower resolution but less memory used, seems to make no real difference
	typedef unsigned index_type_t;

	sd_sphere_d() : pos(all_zeros), radius(0.0), def_pert(0.0), perturb_map(NULL), surf(NULL) {}
	sd_sphere_d(point const &p, float r, int n, float const *pm=NULL, float dp=0.0, upsurface const *const s=NULL) {set_data(p, r, n, pm, dp, s);}
	void gen_points_norms_static(float s_beg=0.0, float s_end=1.0, float t_beg=0.0, float t_end=1.0);
	void gen_points_norms(sphere_point_norm &cur_spn, float s_beg=0.0, float s_end=1.0, float t_beg=0.0, float t_end=1.0);
	float get_rmax() const;
	void draw_subdiv_sphere(point const &vfrom, int texture, bool disable_bfc, unsigned char const *const render_map=NULL,
		float const *const exp_map=NULL, point const *const pt_shift=NULL, float expand=0.0,
		float s_beg=0.0, float s_end=1.0, float t_beg=0.0, float t_end=1.0, bool back_to_front=0) const;
	void get_quad_points(vector<vert_norm_tc> &quad_pts, vector<unsigned> *indices=nullptr, bool use_tri_strip=0, float s_beg=0.0, float s_end=1.0, float t_beg=0.0, float t_end=1.0) const;
	void get_triangles(vector<vert_wrap_t> &verts) const;
	void get_triangle_strip_pow2(vector<vertex_type_t> &verts, unsigned skip=1, float s_beg=0.0, float s_end=1.0, float t_beg=0.0, float t_end=1.0) const;
	void get_triangle_vertex_list(vector<vertex_type_t> &verts) const;
	void get_triangle_index_list_pow2(vector<index_type_t> &indices, unsigned skip=1) const;
	void get_faceted_triangles(vector<vertex_type_t> &verts) const;
	void set_data(point const &p, float r, int n, float const *pm, float dp=0.0, upsurface const *const s=NULL);
	bool equal(point const &p, float r, int n) const {return (p == pos && r == radius && n == (int)ndiv && points != NULL);}
};


class sd_sphere_vbo_d : public sd_sphere_d, public indexed_vao_manager_t {

	vector<unsigned> ix_offsets;
	bool faceted;

	void ensure_vbos();
	unsigned draw_setup(unsigned draw_ndiv);
	unsigned get_index_type_enum()          const {return ((sizeof(index_type_t) == 4) ? GL_UNSIGNED_INT : GL_UNSIGNED_SHORT);}
	unsigned get_count       (unsigned lod) const {assert(lod+1 < ix_offsets.size()); return (ix_offsets[lod+1] - ix_offsets[lod]);}
	void const* get_index_ptr(unsigned lod) const {assert(lod+1 < ix_offsets.size()); return (void const *)(ix_offsets[lod]*sizeof(index_type_t));}

public:
	sd_sphere_vbo_d() : faceted(0) {}
	sd_sphere_vbo_d(point const &p, float r, int n, float const *pm=NULL, float dp=0.0, upsurface const *const s=NULL)
		: sd_sphere_d(p, r, n, pm, dp, s), faceted(0) {}
	void make_faceted() {faceted = 1;}
	void clear_vbos();
	void draw_ndiv_pow2_vbo(unsigned draw_ndiv);
	void draw_instances(unsigned draw_ndiv, instance_render_t &inst_render);
};

