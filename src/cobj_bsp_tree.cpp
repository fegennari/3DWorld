// 3D World - OpenGL CS184 Computer Graphics Project - collision detection BSP Tree
// by Frank Gennari
// 10/16/10

#include "3DWorld.h"
#include "physics_objects.h"


bool const BUILD_COBJ_TREE = 0;


extern int display_mode;
extern vector<coll_obj> coll_objects;


// 3: BSP Tree, 8: OctTree
template<unsigned NUM> class cobj_tree_t {
protected:
	template<unsigned NUM> struct tree_node : public cube_t { // size = 20 - 40
		unsigned start, end; // index into cixs for leaves
		unsigned b[NUM]; // branches

		tree_node(unsigned s=0, unsigned e=0) : start(s), end(e) {
			UNROLL_3X(d[i_][0] = d[i_][1] = 0.0;)
			for (unsigned i = 0; i < NUM; ++i) b[i] = 0; // kids of 0 means empty kids
		}
	};

	struct coll_state {
		point const &p1, p2;
		point &cpos;
		vector3d &cnorm;
		int &cindex;
		int ignore_cobj;
		bool exact;
		float tmin, tmax;

		coll_state(point const &p1_, point const &p2_, point &cp, vector3d &cn, int &ci, int ic, bool ex)
			: p1(p1_), p2(p2_), cpos(cp), cnorm(cn), cindex(ci), ignore_cobj(ic), exact(ex), tmin(0.0), tmax(1.0) {}
		void update_cpos() {cpos = p1 + (p2 - p1)*tmax;}
	};

	vector<coll_obj> const &cobjs;
	vector<unsigned> cixs;
	vector<tree_node<NUM> > nodes;
	unsigned root_node;

	coll_obj const &get_cobj(unsigned ix) const {
		assert(ix < cixs.size() && cixs[ix] < cobjs.size());
		return cobjs[cixs[ix]];
	}
	void mark_as_bin(unsigned ix, unsigned bix) {cixs[ix] |= (bix << 29);}
	void unmark_as_bin(unsigned ix)             {cixs[ix] &= 0x07FFFFFF;}
	unsigned get_bin_ix(unsigned ix)            {return (cixs[ix] >> 29);}

	// Note: start by adding all cobjs to [start, end], then calc bbox, then split into branches and leaves
	void calc_node_bbox(tree_node<NUM> &n) const {
		assert(n.start < n.end);
		n.copy_from(get_cobj(n.start));

		for (unsigned i = n.start+1; i < n.end; ++i) { // bbox union
			coll_obj const &cobj(get_cobj(i));
			UNROLL_3X(n.d[i_][0] = min(n.d[i_][0], cobj.d[i_][0]);)
			UNROLL_3X(n.d[i_][1] = max(n.d[i_][1], cobj.d[i_][1]);)
		}
	}

	void build_tree(unsigned nix, unsigned skip_dims) {assert(0);}

	bool test_cobj(coll_state &state, unsigned ix) const {
		float t(0.0);
		// Note: we test cobj against the original (unclipped) p1 and p2 so that t is correct
		bool const ret(get_cobj(ix).line_int_exact(state.p1, state.p2, t, state.cnorm, state.tmin, state.tmax));
		
		if (ret) {
			state.tmax   = t;
			state.cindex = cixs[ix];
			state.update_cpos();
		}
		return ret;
	}

	bool search_tree(coll_state &state, unsigned nix) const {
		assert(nix < nodes.size());
		tree_node<NUM> const &n(nodes[nix]);
		assert(n.start <= n.end);
		if (!check_line_clip(state.p1, state.cpos, n.d)) return 0;
		bool ret(0);

		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			ret |= test_cobj(state, i);
			if (!state.exact && ret) return 1;
		}

		// interate in correct order based on (p2 - p1)?
		for (unsigned i = 0; i < NUM; ++i) { // check branches
			if (n.b[i] != 0) ret |= search_tree(state, n.b[i]);
			if (!state.exact && ret) return 1;
		}
		return ret;
	}


public:
	cobj_tree_t(vector<coll_obj> const &cobjs_) : cobjs(cobjs_), root_node(0) {}

	void clear() {
		nodes.clear();
		cixs.clear();
		root_node = 0;
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
		root_node = nodes.size();
		nodes.push_back(tree_node<NUM>(0, cixs.size()));
		build_tree(root_node, 0);
		PRINT_TIME("Cobj BSP Tree Create");
		cout << "cobjs: " << cobjs.size() << ", leaves: " << cixs.size() << ", nodes: " << nodes.size() << endl; // testing
	}

	bool check_coll_line(point p1, point p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj, bool exact) const {
		cindex = -1;
		if (nodes.empty()) return 0;
		cpos = p2;
		coll_state state(p1, p2, cpos, cnorm, cindex, ignore_cobj, exact);
		return search_tree(state, root_node);
	}
}; // cobj_tree_t


// BSP Tree (left, mid, right) kids
template <> void cobj_tree_t<3>::build_tree(unsigned nix, unsigned skip_dims) {
	
	assert(nix < nodes.size());
	tree_node<3> &n(nodes[nix]);
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
		nodes[nix].b[bix] = kid;
		nodes.push_back(tree_node<3>(cur, cur+count));
		build_tree(kid, skip_dims);
		cur += count;
	}
	tree_node<3> &n2(nodes[nix]); // create a new reference
	assert(cur == n2.end);
	n2.start = n2.end = 0; // branch node has no leaves
}


// OctTree
template <> void cobj_tree_t<8>::build_tree(unsigned nix, unsigned skip_dims) {

	assert(nix < nodes.size());
	tree_node<8> &n(nodes[nix]);
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
		nodes[nix].b[bix] = kid;
		nodes.push_back(tree_node<8>(cur, cur+count));
		build_tree(kid, (count == num)); // if all in one bin, make that bin a leaf
		cur += count;
	}
	tree_node<8> &n2(nodes[nix]); // create a new reference
	assert(cur == n2.end);
	n2.start = n2.end = 0; // branch node has no leaves
}


cobj_tree_t<8> cobj_tree(coll_objects); // 3: BSP Tree, 8: Octtree


void build_cobj_tree() {

	if (BUILD_COBJ_TREE) cobj_tree.add_cobjs();
}

// can use with ray trace lighting, snow collision, maybe water reflections
bool check_coll_line_exact_tree(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj) {

	return cobj_tree.check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 1);
}

// can use with snow shadows, grass shadows, tree leaf shadows
bool check_coll_line_tree(point const &p1, point const &p2, int &cindex, int ignore_cobj) {

	vector3d cnorm; // unused
	point cpos; // unused
	return cobj_tree.check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 0);
}



