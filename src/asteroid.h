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
	void draw(point_d const &pos_) const;
};


class uasteroid_field : public uobject_base, public vector<uasteroid> {
	
	float max_aradius;

public:
	uasteroid_field() : max_aradius(0.0) {}
	void gen_asteroids(unsigned num);
	void draw(point_d const &pos_, point const &camera) const;
	void free() {clear();}
};


#endif // _ASTEROID_H_


