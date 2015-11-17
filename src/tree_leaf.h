// 3D World - Tree leaf texture/color definitions
// by Frank Gennari
// 7/5/06

#ifndef _TREE_LEAF_H_
#define _TREE_LEAF_H_

#include "3DWorld.h"


enum {TREE_MAPLE = 0, TREE_LIVE_OAK, TREE_A, TREE_B, PAPAYA, NUM_TREE_TYPES};

point const leaf_points[4] = {point(-1,0,0), point(-1,0,2), point(1,0,2), point(1,0,0)};

// should have this for small trees as well
struct tree_type {

	int bark_tex, leaf_tex;
	float branch_size, branch_radius, leaf_size, leaf_x_ar, height_scale, branch_break_off, branch_tscale, branch_color_var, bush_prob;
	colorRGBA barkc, leafc;

	tree_type(int bt, int lt, float bsz, float br, float lsz, float ar, float hs, float bbo, float bts, float bcv, float bp, colorRGBA const &bc, colorRGBA const &lc)
		: bark_tex(bt), leaf_tex(lt), branch_size(bsz), branch_radius(br), leaf_size(lsz), leaf_x_ar(ar), height_scale(hs),
		branch_break_off(bbo), branch_tscale(bts), branch_color_var(bcv), bush_prob(bp), barkc(bc), leafc(lc) {}
};

extern tree_type tree_types[];


#endif // _TREE_LEAF_H_

