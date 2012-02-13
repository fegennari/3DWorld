// 3D World
// by Frank Gennari
// CSG class definitions
// 1/14/06
#ifndef _CSG_H_
#define _CSG_H_

#include "3DWorld.h"

// too large and there can be "holes" in the geometry (usually due to object destruction)
// too small and there will be performance problems due to FP error
float const TOLER = 1.0E-6;


class rect {

public:
	float d[2][2];

	rect() {}

	rect(float const r[2][2]) {
		d[0][0] = r[0][0]; d[0][1] = r[0][1]; d[1][0] = r[1][0]; d[1][1] = r[1][1];
		assert(nonzero());
	}
	rect(float const r[3][2], unsigned d0, unsigned d1);
	float area()   const {return (d[0][1] - d[0][0])*( d[1][1] - d[1][0]);}
	bool nonzero() const {return (d[0][1] > d[0][0] && d[1][1] > d[1][0]);}
	
	float clipped_area(float const c[2]) const {
		return ((d[0][0] >= c[1] || d[0][1] <= c[0]) ? 0.0 : (min(c[1], d[0][1]) - max(c[0], d[0][0]))*(d[1][1] - d[1][0]));
	}
	bool is_near_zero_area() const {
		return (fabs(d[0][0] - d[0][1]) < TOLER || fabs(d[1][0] - d[1][1]) < TOLER);
	}
	bool equal(float const r[2][2]) const {
		return (d[0][0] == r[0][0] && d[0][1] == r[0][1] && d[1][0] == r[1][0] && d[1][1] == r[1][1]);
	}
	bool contains_pt(float const pt[2]) const {
		return (d[0][0] <= pt[0]   && d[0][1] >= pt[0]   && d[1][0] <= pt[1]   && d[1][1] >= pt[1]);
	}
	bool contains(float const r[2][2]) const {
		return (d[0][0] <= r[0][0] && d[0][1] >= r[0][1] && d[1][0] <= r[1][0] && d[1][1] >= r[1][1]);
	}
	bool overlaps(float const r[2][2]) const {
		return (d[0][0] < r[0][1]  && d[0][1] > r[0][0]  && d[1][0] < r[1][1]  && d[1][1] > r[1][0]);
	}
	void clip_to(float const c[2][2]);
	void subtract_from(rect const &rr, deque<rect> &new_rects) const;
	bool merge_with(rect const &r);
	void print() const;
}; // class rect


class csg_cube : public cube_t {

	unsigned char eflags;

	csg_cube(unsigned char eflags0) : eflags(eflags0) {} // eflags constructor (internal)
	bool subtract_from_internal(const csg_cube &cube, vector<csg_cube> &output) const;

public:
	csg_cube() : eflags(0) {}
	csg_cube(float x1, float x2, float y1, float y2, float z1, float z2)
		: cube_t(x1, x2, y1, y2, z1, z2), eflags(0) {} // float constructor
	csg_cube(const coll_obj &cobj, bool use_bounding_cube=0);
	csg_cube(cube_t const &cube, unsigned char eflags0=0) : cube_t(cube), eflags(eflags0) {}
	void write_to_cobj(coll_obj &cobj) const;
	bool cube_intersection(const csg_cube &cube, csg_cube &res) const;
	bool subtract_from_cube(coll_obj_group &new_cobjs, coll_obj const &cobj) const;
	bool subtract_from_cylinder(coll_obj_group &new_cobjs, coll_obj &cobj) const;
	bool subtract_from_polygon(coll_obj_group &new_cobjs, coll_obj const &cobj) const;
	bool subtract_from_thick_polygon(coll_obj_group &new_cobjs, coll_obj const &cobj) const;
	bool cube_merge(csg_cube &cube); // const cube?
	void unset_adjacent_edge_flags(coll_obj &cobj) const;
	void unset_intersecting_edge_flags(coll_obj &cobj) const;
	float get_d(unsigned dim, bool dir) const {assert(dim >= 0 && dim < 3); return d[dim][dir];}
}; // class csg_cube


class r_profile {

	bool filled;
	rect bb;
	float tot_area, avg_alpha;
	vector<rect> rects;
	deque<rect> pend;

	void add_rect_int(rect const &r);

public:
	r_profile() : filled(0), tot_area(0.0), avg_alpha(1.0) {}

	r_profile(float const bb_[2][2]) : filled(0), bb(bb_), tot_area(bb.area()), avg_alpha(1.0) {	
		assert(tot_area > 0.0);
	}
	bool is_filled() const {return filled;}
	void reset_bbox(float const bb_[2][2]);
	bool add_rect(float const d[3][2], unsigned d0, unsigned d1, float alpha);
	float den_inv() {return clipped_den_inv(bb.d[0]);}
	float clipped_den_inv(float const c[2]) const; // clip by first dimension
	void clear();
	void clear_within(float const c[2]);
};


#endif

