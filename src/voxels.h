// 3D World - Voxel Header
// by Frank Gennari
// 2/25/12
#ifndef _VOXELS_H_
#define _VOXELS_H_

#include "3DWorld.h"
#include "model3d.h"


struct voxel_render_params_t {

	unsigned tids[2];
	colorRGBA colors[2];
	colorRGBA base_color;
	voxel_render_params_t() {tids[0] = tids[1] = 0; colors[0] = colors[1] = base_color = WHITE;}

	voxel_render_params_t(unsigned t1, unsigned t2, colorRGBA const &c1, colorRGBA const &c2, colorRGBA const &bc) : base_color(bc)
	{
		tids  [0] = t1; tids  [1] = t2;
		colors[0] = c1; colors[1] = c2;
	}
};


struct voxel_params_t {

	float isolevel, elasticity;
	bool make_closed_surface, invert, remove_unconnected, keep_at_scene_edge, remove_under_mesh;
	voxel_render_params_t rp;

	voxel_params_t(float il=0.0, float e=0.5, bool mcs=0, bool inv=0, bool ru=0, bool kase=0, bool rum=0)
		: isolevel(il), elasticity(e), make_closed_surface(mcs), invert(inv), remove_unconnected(ru),
		keep_at_scene_edge(kase), remove_under_mesh(rum) {}
};


class noise_texture_manager_t {

	unsigned noise_tid, tsize;

public:
	noise_texture_manager_t() : noise_tid(0), tsize(0) {}
	void setup(unsigned size, float mag=1.0, float freq=1.0, vector3d const &offset=zero_vector);
	void bind_texture(unsigned tu_id) const;
	void clear();
};


// stored internally in yxz order
template<typename V> class voxel_grid : public vector<V> {
public:
	unsigned nx, ny, nz;
	vector3d vsz; // size of a voxel in x,y,z
	point center, lo_pos;

	voxel_grid() : nx(0), ny(0), nz(0), vsz(zero_vector) {}
	void init(unsigned nx_, unsigned ny_, unsigned nz_, vector3d const &vsz_, point const &center_);

	bool get_ix(point const &p, unsigned &ix) const { // returns whether or not the point was inside the voxel volume
		int i[3]; // x,y,z
		UNROLL_3X(i[i_] = int((p[i_] - lo_pos[i_])/vsz[i_]);); // convert to voxel space
		if (i[0] < 0 || i[1] < 0 || i[2] < 0 || i[0] >= (int)nx || i[1] >= (int)ny || i[2] >= (int)nz) return 0;
		ix = i[2] + (i[0] + i[1]*nx)*nz;
		return 1;
	}
	unsigned get_ix(unsigned x, unsigned y, unsigned z) const {
		assert(x < nx && y < ny && z < nz);
		return (z + (x + y*nx)*nz);
	}
	V const &get(unsigned x, unsigned y, unsigned z) const {
		unsigned const ix(get_ix(x, y, z));
		assert(ix < size());
		return operator[](ix);
	}
	V &get_ref(unsigned x, unsigned y, unsigned z) {
		unsigned const ix(get_ix(x, y, z));
		assert(ix < size());
		return operator[](ix);
	}
	void set(unsigned x, unsigned y, unsigned z, V const &val) {
		unsigned const ix(get_ix(x, y, z));
		assert(ix < size());
		operator[](ix) = val;
	}
};

typedef voxel_grid<float> float_voxel_grid;


class voxel_manager : public float_voxel_grid {

	voxel_grid<unsigned char> outside;

	point interpolate_pt(float isolevel, point const &pt1, point const &pt2, float const val1, float const val2) const;

protected:
	void get_triangles_for_voxel(vector<triangle> &triangles, voxel_params_t const &vp, unsigned x, unsigned y, unsigned z) const;

public:
	void clear() {outside.clear(); float_voxel_grid::clear();}
	void create_procedural(float mag, float freq, vector3d const &offset, bool normalize_to_1, int rseed1, int rseed2);
	void atten_at_edges(float val);
	void atten_at_top_only(float val);
	void determine_voxels_outside(voxel_params_t const &vp);
	void remove_unconnected_outside(bool keep_at_scene_edge);
	void create_triangles(vector<triangle> &triangles, voxel_params_t const &vp) const;
	void get_triangles(vector<triangle> &triangles, voxel_params_t const &vp);
	bool point_inside_volume(point const &pos) const;
};


struct coll_tquad;


class voxel_model : public voxel_manager {

	voxel_render_params_t params;
	typedef vert_norm vertex_type_t;
	typedef vntc_vect_block_t<vertex_type_t> tri_data_t;
	tri_data_t tri_data;
	noise_texture_manager_t noise_tex_gen;

	struct pt_ix_t {
		point pt;
		unsigned ix;
		pt_ix_t(point const &pt_=all_zeros, unsigned const ix_=0) : pt(pt_), ix(ix_) {}
	};

	vector<pt_ix_t> pt_to_ix;

	struct comp_by_dist {
		point const p;
		comp_by_dist(point const &p_) : p(p_) {}
		bool operator()(pt_ix_t const &a, pt_ix_t const &b) const {return (p2p_dist_sq(a.pt, p) < p2p_dist_sq(b.pt, p));}
	};

public:
	void clear();
	void build(voxel_params_t const &vp, bool add_cobjs);
	void render(bool is_shadow_pass);
	void free_context();
};

#endif

