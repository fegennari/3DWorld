// 3D World - Universe Asteroid Header
// by Frank Gennari
// 12/15/12

#ifndef _ASTEROID_H_
#define _ASTEROID_H_

#include "universe.h"

class uasteroid : public uobject_base, public rotated_obj {
	
	unsigned inst_id;
	vector3d scale;

public:
	uasteroid() : inst_id(0) {}
	void gen(float max_dist, float max_radius);
	void draw(point_d const &pos_, point const &camera, shader_t &s) const;
	vector3d const &get_scale() const {return scale;}
};


class uasteroid_field : public uobject_base, public vector<uasteroid> {

public:
	static void begin_render(shader_t &shader);
	static void end_render(shader_t &shader);
	void gen_asteroids();
	void draw(point_d const &pos_, point const &camera, shader_t &s);
	void free() {clear();}
};


#endif // _ASTEROID_H_


