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
	float branch_size, branch_radius, leaf_size, leaf_x_ar, height_scale, branch_break_off, branch_tscale, branch_color_var;
	colorRGBA barkc, leafc;

	tree_type(int bt, int lt, float bsz, float br, float lsz, float ar, float hs, float bbo, float bts, float bcv, colorRGBA const &bc, colorRGBA const &lc)
		: bark_tex(bt), leaf_tex(lt), branch_size(bsz), branch_radius(br), leaf_size(lsz), leaf_x_ar(ar), height_scale(hs),
		branch_break_off(bbo), branch_tscale(bts), branch_color_var(bcv), barkc(bc), leafc(lc) {}
};

// bark_tex, leaf_tex, branch_size, branch_radius, leaf_size, leaf_x_ar, height_scale, branch_break_off, branch_tscale, branch_color_var, barkc, leafc
tree_type const tree_types[NUM_TREE_TYPES] = {
	tree_type(BARK3_TEX, LEAF_TEX,     1.0, 0.7, 1.0, 1.00, 1.0, 1.0, 1.0, 0.1,  colorRGBA(0.7, 0.7,  0.5,  1.0), colorRGBA(0.2, 1.0, 0.2, 1.0)),
	tree_type(BARK4_TEX, LIVE_OAK_TEX, 1.0, 1.0, 1.0, 0.63, 1.0, 1.0, 1.5, 0.1,  colorRGBA(1.0, 0.9,  0.8,  1.0), WHITE),
	tree_type(BARK1_TEX, LEAF2_TEX,    1.0, 1.0, 1.0, 0.82, 2.0, 0.5, 1.0, 0.1,  colorRGBA(0.8, 0.5,  0.3,  1.0), WHITE),
	tree_type(BARK5_TEX, LEAF3_TEX,    1.0, 0.7, 1.5, 0.81, 1.0, 1.0, 0.3, 0.01, colorRGBA(0.8, 0.75, 0.65, 1.0), WHITE), // birch bark
	tree_type(BARK6_TEX, PAPAYA_TEX,   1.0, 1.0, 1.0, 1.00, 2.0, 2.0, 0.5, 0.1,  colorRGBA(0.7, 0.6,  0.5,  1.0), WHITE)
};


#endif // _TREE_LEAF_H_

