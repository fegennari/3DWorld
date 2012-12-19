// 3D World - Cloud and Nebula header
// by Frank Gennari
// 12/16/12

#ifndef _CLOUDS_H_
#define _CLOUDS_H_

#include "universe_base.h"


class unebula : public uobject_base {

	colorRGBA color[2];
	typedef vert_norm_comp vert_type_t;
	vector<vert_type_t> points;

public:
	static void begin_render(shader_t &s);
	static void end_render(shader_t &s);
	void gen(float range, ellipsoid_t const &bounds);
	void draw(point_d pos_, point const &camera, float max_dist, shader_t &s) const;
	void free() {points.clear();}
};


#endif // _CLOUDS_H_


