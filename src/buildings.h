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

void get_building_occluders(pos_dir_up const &pdu, building_occlusion_state_t &state);
bool check_pts_occluded(point const *const pts, unsigned npts, building_occlusion_state_t &state);

#endif // _BUILDING_H_
