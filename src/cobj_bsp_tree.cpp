// 3D World - OpenGL CS184 Computer Graphics Project - collision detection BSP Tree
// by Frank Gennari
// 10/16/10

#include "3DWorld.h"
#include "physics_objects.h"


bool const BUILD_COBJ_TREE = 1;


bool cobj_tree_valid(0);

extern int display_mode, frame_counter, cobj_counter;
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
	bool is_static, is_dynamic;

	coll_obj const &get_cobj(unsigned ix) const {
		//assert(ix < cixs.size() && cixs[ix] < cobjs.size());
		return cobjs[cixs[ix]];
	}
	inline void mark_as_bin(unsigned ix, unsigned bix) {cixs[ix] |= (bix << 29);}
	inline void unmark_as_bin(unsigned ix)             {cixs[ix] &= 0x07FFFFFF;}
	inline unsigned get_bin_ix(unsigned ix)            {return (cixs[ix] >> 29);}

	// Note: start by adding all cobjs to [start, end], then calc bbox, then split into branches and leaves
	void calc_node_bbox(tree_node &n) const {
		assert(n.start < n.end);
		n.copy_from(get_cobj(n.start));

		for (unsigned i = n.start+1; i < n.end; ++i) { // bbox union
			n.union_with_cube(get_cobj(i));
		}
	}

	bool create_cixs() {
		if (is_static) cixs.reserve(cobjs.size());

		for (vector<coll_obj>::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
			if (obj_ok(*i)) cixs.push_back(i - cobjs.begin());
		}
		assert(cixs.size() < (1 << 29));
		return !cixs.empty();
	}

	void build_tree(unsigned nix, unsigned skip_dims) {assert(0);}

	bool obj_ok(coll_obj const &cobj) const {
		return ((is_static && cobj.status == COLL_STATIC) || (is_dynamic && cobj.status == COLL_DYNAMIC));
	}


public:
	cobj_tree_t(vector<coll_obj> const &cobjs_, bool s, bool d) : cobjs(cobjs_), is_static(s), is_dynamic(d) {}

	void clear() {
		nodes.resize(0);
		cixs.resize(0);
	}

	void add_cobjs(bool verbose) {
		RESET_TIME;
		clear();
		if (!create_cixs()) return; // nothing to be done
		nodes.reserve(6*cixs.size()/5); // conservative
		nodes.push_back(tree_node(0, cixs.size()));
		assert(nodes.size() == 1);
		build_tree(0, 0);
		nodes[0].next_node_id = nodes.size();

		if (verbose) {
			PRINT_TIME(" Cobj Tree Create");
			cout << "cobjs: " << cobjs.size() << ", leaves: " << cixs.size() << ", nodes: " << nodes.size() << endl;
		}
	}

	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj, bool exact, int test_alpha) const {
		cindex = -1;
		if (nodes.empty()) return 0;
		bool ret(0);
		float t(0.0), tmin(0.0), tmax(1.0);
		vector3d dinv(p2 - p1);
		dinv.invert(0, 1);
		assert(test_alpha != 2); // don't support this mode
		
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
				if (!obj_ok(cobj))               continue;
				if (test_alpha == 1 && cobj.is_semi_trans())                   continue; // semi-transparent, can see through
				if (test_alpha == 3 && cobj.cp.color.alpha < MIN_SHADOW_ALPHA) continue; // less than min alpha
				if (!cobj.line_int_exact(p1, p2, t, cnorm, tmin, tmax))        continue;
				cindex = cixs[i];
				cpos   = p1 + (p2 - p1)*t;
				if (!exact) return 1; // return first hit
				dinv   = vector3d(cpos - p1);
				dinv.invert(0, 1);
				tmax   = t;
				ret    = 1;
			}
			++nix;
		}
		return ret;
	}

	void get_intersecting_cobjs(cube_t const &cube, vector<unsigned> &cobjs, int ignore_cobj, float toler, bool check_ccounter, int id_for_cobj_int) const {
		for (unsigned nix = 0; nix < nodes.size();) {
			tree_node const &n(nodes[nix]);
			assert(n.start <= n.end);

			if (!cube.intersects(n, toler)) {
				assert(n.next_node_id > nix);
				nix = n.next_node_id; // failed the bbox test
				assert(nix > 0);
				continue;
			}
			for (unsigned i = n.start; i < n.end; ++i) { // check leaves
				if ((int)cixs[i] == ignore_cobj) continue;
				coll_obj const &cobj(get_cobj(i));
				if (check_ccounter && cobj.counter == cobj_counter) continue;
				if (!obj_ok(cobj) || !cube.intersects(cobj, toler)) continue;
				if (id_for_cobj_int >= 0 && coll_objects[id_for_cobj_int].intersects_cobj(cobj, toler) != 1) continue;
				cobjs.push_back(cixs[i]);
			}
			++nix;
		}
	}
}; // cobj_tree_t


// BSP Tree (left, right, mid) kids
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
	point const sval(n.get_cube_center()); // center point

	// split in this dimension: use upper 3 bits of cixs for storing bin index
	for (unsigned i = n.start; i < n.end; ++i) {
		point const center(get_cobj(i).get_cube_center());
		unsigned bix(0);
		UNROLL_3X(if (center[i_] > sval[i_]) bix |= (1 << i_);)
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


// 3: BSP Tree, 8: Octtree
typedef cobj_tree_t<8> cobj_tree_type;
cobj_tree_type cobj_tree_static (coll_objects, 1, 0);
cobj_tree_type cobj_tree_dynamic(coll_objects, 0, 1);
int last_update_frame[2] = {1, 1}; // first drawn frame is 1


cobj_tree_type &get_tree(bool dynamic) {
	return (dynamic ? cobj_tree_dynamic : cobj_tree_static);
}

void build_cobj_tree(bool dynamic, bool verbose) {
	if (BUILD_COBJ_TREE) {
		get_tree(dynamic).add_cobjs(verbose);
		cobj_tree_valid = 1;
	}
}

void update_cobj_tree(bool dynamic, bool verbose) {

	if (last_update_frame[dynamic] < frame_counter && (dynamic || !cobj_tree_valid)) {
		last_update_frame[dynamic] = frame_counter;
		build_cobj_tree(dynamic, verbose);
	}
}

// can use with ray trace lighting, snow collision?, maybe water reflections
bool check_coll_line_exact_tree(point const &p1, point const &p2, point &cpos,
								vector3d &cnorm, int &cindex, int ignore_cobj, bool dynamic, int test_alpha)
{
	return get_tree(dynamic).check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 1, test_alpha);
}

// can use with snow shadows, grass shadows, tree leaf shadows
bool check_coll_line_tree(point const &p1, point const &p2, int &cindex, int ignore_cobj, bool dynamic, int test_alpha) {

	vector3d cnorm; // unused
	point cpos; // unused
	return get_tree(dynamic).check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 0, test_alpha);
}


void get_intersecting_cobjs_tree(cube_t const &cube, vector<unsigned> &cobjs, int ignore_cobj, float toler,
	bool dynamic, bool check_ccounter, int id_for_cobj_int)
{
	get_tree(dynamic).get_intersecting_cobjs(cube, cobjs, ignore_cobj, toler, check_ccounter, id_for_cobj_int);
}



