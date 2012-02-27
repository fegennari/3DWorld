// 3D World - Voxel Header
// by Frank Gennari
// 2/25/12
#ifndef _VOXELS_H_
#define _VOXELS_H_

#include "3DWorld.h"


struct float_voxel_t {
	float v; // right now just a floating-point number
	float_voxel_t(float v_=0.0) : v(v_) {}
};


template<typename V> class voxel_grid : public vector<V> {
public:
	unsigned nx, ny, nz;
	vector3d vsz; // size of a voxel in x,y,z
	point center, lo_pos;

	voxel_grid() : nx(0), ny(0), nz(0), vsz(zero_vector) {}
	void init(unsigned nx_, unsigned ny_, unsigned nz_, vector3d const &vsz_, point const &center_);

	unsigned get_ix(unsigned x, unsigned y, unsigned z) const {
		assert(x < nx && y < ny && z < nz);
		return (x + (y + z*ny)*nx);
	}
	V const &get(unsigned x, unsigned y, unsigned z) const {
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

typedef voxel_grid<float_voxel_t> float_voxel_grid;


class voxel_manager : public float_voxel_grid {

	point interpolate_pt(float isolevel, point const &pt1, point const &pt2, float const val1, float const val2) const;

public:
	void create_procedural(float mag, float freq);
	void add_triangles_from_voxel(vector<triangle> &triangles, unsigned x, unsigned y, unsigned z, float isolevel=0.0) const;
	void voxel_manager::get_triangles(vector<triangle> &triangles, float isolevel=0.0) const;
};

#endif

