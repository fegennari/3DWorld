// 3D World - Tree leaf texture/color definitions
// by Frank Gennari
// 7/5/06

#ifndef _TREE_LEAF_H_
#define _TREE_LEAF_H_

#include "3DWorld.h"


enum {TREE_MAPLE = 0, TREE_LIVE_OAK, TREE_A, TREE_B, PAPAYA, NUM_TREE_TYPES};

// should have this for small trees as well
struct tree_type {

	int bark_tex, leaf_tex;
	float size;
	colorRGBA barkc, leafc;

	tree_type(int bt, int lt, float sz, colorRGBA const &bc, colorRGBA const &lc)
		: bark_tex(bt), leaf_tex(lt), size(sz), barkc(bc), leafc(lc) {}
};

// bark_tex, leaf_tex, barkc, leafc
tree_type const tree_types[NUM_TREE_TYPES] = {
	tree_type(WOOD_TEX, LEAF_TEX,     1.0, colorRGBA(0.7, 0.7, 0.7, 1.0), colorRGBA(0.2, 1.0, 0.2, 1.0)),
	tree_type(WOOD_TEX, LIVE_OAK_TEX, 1.0, colorRGBA(0.7, 0.6, 0.5, 1.0), WHITE),
	tree_type(WOOD_TEX, LEAF2_TEX,    1.0, colorRGBA(0.7, 0.4, 0.0, 1.0), WHITE),
	tree_type(WOOD_TEX, LEAF3_TEX,    1.0, colorRGBA(0.8, 0.6, 0.4, 1.0), WHITE),
	tree_type(WOOD_TEX, PAPAYA_TEX,   1.0, colorRGBA(0.6, 0.5, 0.3, 1.0), WHITE)
};


#endif // _TREE_LEAF_H_

