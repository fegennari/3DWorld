// 3D World - Mesh Intersection classes header
// by Frank Gennari
// 6/30/06
#ifndef _MESH_INTERSECT_H_
#define _MESH_INTERSECT_H_

#include "3DWorld.h"


struct bsp_tree_node {
	cube_t c;
};

struct mesh_query_ret {
	int xpos, ypos;
	float zval;
	mesh_query_ret() : xpos(0), ypos(0), zval(0.0) {}
};


class mesh_bsp_tree; // forward reference


class mesh_intersector {

	point v1, v2;
	vector3d v2_v1;
	int fast;
	mesh_query_ret ret;

	bool line_int_surface_cached();
	bool line_intersect_surface();
	bool line_intersect_surface_fast();
	bool check_iter_clip(int fast2);
	bool line_intersect_plane(int x1, int x2, int y1, int y2);
	bool intersect_plane(point const &mesh_pt, vector3d const &norm, bool sign);

public:
	mesh_intersector(point const &v1_, point const &v2_, int fast_) :
		v1(v1_), v2(v2_), v2_v1(v2 - v1), fast(fast_) {}
	bool get_intersection(int &xpos_, int &ypos_, float &zval_, bool cached);
	bool get_intersection();
	bool get_any_non_intersection(point const *const pts, unsigned npts);
	bool intersect_mesh_quad(int x, int y);
	mesh_query_ret const &get_ret() const {return ret;}
};


class mesh_bsp_tree {

	bool dir0; // first subdivision direction
	unsigned nlevels; // inclusive
	vector<bsp_tree_node> bsp_data;
	vector<bsp_tree_node *> tree; // {level, y/x}

	mesh_bsp_tree(mesh_bsp_tree const &); // forbidden
	void operator=(mesh_bsp_tree const &); // forbidden

	bool search_recur(point v1, point v2, unsigned x, unsigned y, unsigned level, mesh_query_ret &ret) const; // Note: thread safe

public:
	mesh_bsp_tree();
	bool search(point const &v1, point const &v2, mesh_query_ret &ret) const {return search_recur(v1, v2, 0, 0, 0, ret);} // Note: thread safe
	bsp_tree_node const &get_root() const {assert(!tree.empty()); return tree[0][0];}
};


#endif _MESH_INTERSECT_H_


