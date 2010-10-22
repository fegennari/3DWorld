// 3D World - OpenGL CS184 Computer Graphics Project - collision detection BSP Tree
// by Frank Gennari
// 10/16/10

#include "3DWorld.h"
#include "physics_objects.h"


bool const BUILD_COBJ_TREE = 1;


int last_update_frame(1); // first drawn frame is 1

extern int display_mode, frame_counter;
extern vector<coll_obj> coll_objects;


// 3: BSP Tree, 8: OctTree
template<unsigned NUM> class cobj_tree_t {
protected:
	struct tree_node : public cube_t { // size = 36
		unsigned start, end; // index into cixs for leaves
		unsigned next_node_id;

		tree_node(unsigned s=0, unsigned e=0) : start(s), end(e), next_node_id(0) {
			UNROLL_3X(d[i_][0] = d[i_][1] = 0.0;)
		}
	};

	vector<coll_obj> const &cobjs;
	vector<unsigned> cixs;
	vector<tree_node> nodes;

	coll_obj const &get_cobj(unsigned ix) const {
		assert(ix < cixs.size() && cixs[ix] < cobjs.size());
		return cobjs[cixs[ix]];
	}
	void mark_as_bin(unsigned ix, unsigned bix) {cixs[ix] |= (bix << 29);}
	void unmark_as_bin(unsigned ix)             {cixs[ix] &= 0x07FFFFFF;}
	unsigned get_bin_ix(unsigned ix)            {return (cixs[ix] >> 29);}

	// Note: start by adding all cobjs to [start, end], then calc bbox, then split into branches and leaves
	void calc_node_bbox(tree_node &n) const {
		assert(n.start < n.end);
		n.copy_from(get_cobj(n.start));

		for (unsigned i = n.start+1; i < n.end; ++i) { // bbox union
			coll_obj const &cobj(get_cobj(i));
			UNROLL_3X(n.d[i_][0] = min(n.d[i_][0], cobj.d[i_][0]);)
			UNROLL_3X(n.d[i_][1] = max(n.d[i_][1], cobj.d[i_][1]);)
		}
	}

	void build_tree(unsigned nix, unsigned skip_dims) {assert(0);}


public:
	cobj_tree_t(vector<coll_obj> const &cobjs_) : cobjs(cobjs_) {}

	void clear() {
		nodes.clear();
		cixs.clear();
	}

	void add_cobjs() {
		RESET_TIME;
		clear();
		cixs.reserve(cobjs.size());

		for (vector<coll_obj>::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
			if (i->status != COLL_STATIC) continue;
			cixs.push_back(i - cobjs.begin());
		}
		if (cixs.empty()) return; // nothing to be done
		assert(cixs.size() < (1 << 29));
		nodes.reserve(6*cixs.size()/5); // conservative
		nodes.push_back(tree_node(0, cixs.size()));
		assert(nodes.size() == 1);
		build_tree(0, 0);
		nodes[0].next_node_id = nodes.size();
		PRINT_TIME("Cobj BSP Tree Create");
		cout << "cobjs: " << cobjs.size() << ", leaves: " << cixs.size() << ", nodes: " << nodes.size() << endl; // testing
	}

	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj, bool exact) const {
		cindex = -1;
		bool ret(0);
		float t(0.0), tmin(0.0), tmax(1.0);
		vector3d dinv(p2 - p1);
		dinv.invert(0, 1);
		
		for (unsigned nix = 0; nix < nodes.size();) {
			tree_node const &n(nodes[nix]);
			assert(n.start <= n.end);

			if (!get_line_clip2(p1, dinv, n.d)) {
				assert(n.next_node_id > nix);
				nix = n.next_node_id; // failed the bbox test
				assert(nix > 0);
				continue;
			}
			for (unsigned i = n.start; i < n.end; ++i) { // check leaves
				// Note: we test cobj against the original (unclipped) p1 and p2 so that t is correct
				// Note: we probably don't need to return cnorm and cpos in inexact mode, but it shouldn't be too expensive to do so
				if ((int)cixs[i] == ignore_cobj) continue;
				coll_obj const &cobj(get_cobj(i));
				bool const coll(cobj.status == COLL_STATIC && cobj.line_int_exact(p1, p2, t, cnorm, tmin, tmax));
				
				if (coll) {
					cindex = cixs[i];
					cpos   = p1 + (p2 - p1)*t;
					if (!exact) return 1; // return first hit
					dinv   = vector3d(cpos - p1);
					dinv.invert(0, 1);
					tmax   = t;
					ret    = 1;
				}
			}
			++nix;
		}
		return ret;
	}
}; // cobj_tree_t


// BSP Tree (left, mid, right) kids
template <> void cobj_tree_t<3>::build_tree(unsigned nix, unsigned skip_dims) {
	
	assert(nix < nodes.size());
	tree_node &n(nodes[nix]);
	assert(n.start < n.end);
	calc_node_bbox(n);
	unsigned const num(n.end - n.start);
	if (num <= 2 || skip_dims == 7) return; // base case
	
	// determine split dimension and value
	unsigned dim(0);
	float max_sz(0);

	for (unsigned d = 0; d < 3; ++d) {
		if (skip_dims & (1 << d)) continue;
		float const dim_sz(n.d[d][1] - n.d[d][0]);
		assert(dim_sz >= 0.0);
		
		if (dim_sz > max_sz) {
			max_sz = dim_sz;
			dim    = d;
		}
	}
	assert(!(skip_dims & (1 << dim)));
	float const sval(0.5*(n.d[dim][1] + n.d[dim][0])); // center point (mean seems to work better than median)

	// split in this dimension: use upper 2 bits of cixs for storing bin index
	for (unsigned i = n.start; i < n.end; ++i) {
		unsigned bix(2);
		coll_obj const &cobj(get_cobj(i));
		assert(cobj.d[dim][0] <= cobj.d[dim][1]);
		if (cobj.d[dim][1] <= sval) bix = 0; // ends   before the split, put in bin 0
		if (cobj.d[dim][0] >= sval) bix = 1; // starts after  the split, put in bin 1
		mark_as_bin(i, bix);
	}
	sort((cixs.begin() + n.start), (cixs.begin() + n.end)); // sort by bix then by ix
	unsigned bin_count[3] = {0};

	for (unsigned i = n.start; i < n.end; ++i) {
		unsigned const bix(get_bin_ix(i));
		assert(bix < 3);
		++bin_count[bix];
		unmark_as_bin(i);
	}

	// check that dataset has been subdivided (not all in one bin)
	assert(bin_count[0] < num && bin_count[1] < num);

	if (bin_count[2] == num) {
		build_tree(nix, (skip_dims | (1 << dim))); // single bin, rebin with a different dim
		return;
	}

	// create child nodes and call recursively
	unsigned cur(n.start);

	for (unsigned bix = 0; bix < 3; ++bix) { // Note: this loop will invalidate the reference to 'n'
		unsigned const count(bin_count[bix]);
		if (count == 0) continue; // empty bin
		unsigned const kid(nodes.size());
		nodes.push_back(tree_node(cur, cur+count));
		build_tree(kid, skip_dims);
		nodes[kid].next_node_id = nodes.size();
		cur += count;
	}
	tree_node &n2(nodes[nix]); // create a new reference
	assert(cur == n2.end);
	n2.start = n2.end = 0; // branch node has no leaves
}


// OctTree
template <> void cobj_tree_t<8>::build_tree(unsigned nix, unsigned skip_dims) {

	assert(nix < nodes.size());
	tree_node &n(nodes[nix]);
	assert(n.start < n.end);
	calc_node_bbox(n);
	unsigned const num(n.end - n.start);
	if (num <= 2 || skip_dims) return; // base case
	
	// determine split values
	point const sval(n.get_center()); // center point

	// split in this dimension: use upper 3 bits of cixs for storing bin index
	for (unsigned i = n.start; i < n.end; ++i) {
		coll_obj const &cobj(get_cobj(i));
		point const center(cobj.get_center());
		unsigned bix(0);

		for (unsigned d = 0; d < 3; ++d) {
			if (center[d] > sval[d]) bix |= (1 << d);
		}
		mark_as_bin(i, bix);
	}
	sort((cixs.begin() + n.start), (cixs.begin() + n.end)); // sort by bix then by ix
	unsigned bin_count[8] = {0};

	for (unsigned i = n.start; i < n.end; ++i) {
		unsigned const bix(get_bin_ix(i));
		assert(bix < 8);
		++bin_count[bix];
		unmark_as_bin(i);
	}

	// create child nodes and call recursively
	unsigned cur(n.start);

	for (unsigned bix = 0; bix < 8; ++bix) { // Note: this loop will invalidate the reference to 'n'
		unsigned const count(bin_count[bix]);
		if (count == 0) continue; // empty bin
		unsigned const kid(nodes.size());
		nodes.push_back(tree_node(cur, cur+count));
		build_tree(kid, (count == num)); // if all in one bin, make that bin a leaf
		nodes[kid].next_node_id = nodes.size();
		cur += count;
	}
	tree_node &n2(nodes[nix]); // create a new reference
	assert(cur == n2.end);
	n2.start = n2.end = 0; // branch node has no leaves
}


cobj_tree_t<8> cobj_tree(coll_objects); // 3: BSP Tree, 8: Octtree


void build_cobj_tree() {

	if (BUILD_COBJ_TREE) cobj_tree.add_cobjs();
	last_update_frame = max(last_update_frame, frame_counter);
}

void update_cobj_tree() {

	if (last_update_frame != frame_counter) build_cobj_tree();
}

// can use with ray trace lighting, snow collision?, maybe water reflections
bool check_coll_line_exact_tree(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj) {

	return cobj_tree.check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 1);
}

// can use with snow shadows, grass shadows, tree leaf shadows
bool check_coll_line_tree(point const &p1, point const &p2, int &cindex, int ignore_cobj) {

	vector3d cnorm; // unused
	point cpos; // unused
	return cobj_tree.check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 0);
}



