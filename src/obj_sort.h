// 3D World - Object performance wrappers/caching and sorting code for universe model
// by Frank Gennari
// 11/01/05

#ifndef _OBJ_SORT_H_
#define _OBJ_SORT_H_

#include "ship.h"

unsigned const LEFT_EDGE_BIT = 0x80000000U;


struct cached_obj : public sphere_t {

	free_obj *obj;
	unsigned flags;

	cached_obj() : obj(NULL), flags(0) {}
	cached_obj(free_obj *ptr) {set_obj(ptr);}

	void set_obj(free_obj *ptr) {
		assert(ptr != NULL);
		if (VERIFY_REFS) ptr->verify_status();
		obj    = ptr;
		radius = ptr->get_c_radius();
		pos    = ptr->get_pos();
		flags  = ptr->get_flags();
		if (ptr->status     == 1) flags |= OBJ_FLAGS_BAD_;
		if (ptr->get_time() == 0) flags |= OBJ_FLAGS_NEW_;
	}
	inline void refresh() {set_obj(obj);}
};


struct interval {

	float val;
	unsigned ix;

	interval() {}
	interval(float v, unsigned i, bool l) : val(v), ix(i) {if (l) ix |= LEFT_EDGE_BIT;}
	bool operator<(interval const &iv) const {return (val < iv.val);}
};


struct comp_co_fast_x {
	bool operator()(cached_obj const &o1, cached_obj const &o2) {
		return (o1.pos.x < o2.pos.x);
	}
};



#endif

