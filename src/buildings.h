// 3D World - Buildings Interface
// by Frank Gennari
// 4-5-18

#ifndef _BUILDING_H_
#define _BUILDING_H_

#include "3DWorld.h"

struct building_occlusion_state_t {
	point pos;
	vector3d xlate;
	vector<unsigned> building_ids;
	vector<point> temp_points;

	void init(point const &pos_, vector3d const &xlate_) {
		pos   = pos_;
		xlate = xlate_;
		building_ids.clear();
	}
};

struct cube_with_zval_t : public cube_t {
	float zval;
	cube_with_zval_t() : zval(0.0) {}
	cube_with_zval_t(cube_t const &c, float zval_=0.0) : cube_t(c), zval(zval_) {}
};

typedef vector<cube_with_zval_t> vect_cube_with_zval_t;

void get_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state);
bool check_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t &state);
bool has_bcube_int_xy(cube_t const &bcube, vect_cube_t const &bcubes, float pad_dist=0.0);

#endif // _BUILDING_H_
