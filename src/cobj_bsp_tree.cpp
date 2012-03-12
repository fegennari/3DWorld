// 3D World - OpenGL CS184 Computer Graphics Project - collision detection BSP/KD/Oct Tree
// by Frank Gennari
// 10/16/10

#include "3DWorld.h"
#include "cobj_bsp_tree.h"


unsigned const MAX_LEAF_SIZE = 2;
float const POLY_TOLER       = 1.0E-6;
float const OVERLAP_AMT      = 0.02;


vector<unsigned> dynamic_ranges;

extern int display_mode, frame_counter, cobj_counter, begin_motion;
extern coll_obj_group coll_objects;
extern vector<unsigned> falling_cobjs;
extern platform_cont platforms;


coll_tquad::coll_tquad(coll_obj const &c) : tquad_t(c.npoints), normal(c.norm), cid(c.id) {

	assert(is_cobj_valid(c));
	for (unsigned i = 0; i < npts; ++i) {pts[i] = c.points[i];}
	if (npts == 3) pts[3] = pts[2]; // duplicate the last point so that it's valid
}


coll_tquad::coll_tquad(polygon_t const &p) : tquad_t((unsigned)p.size()) {

	assert(npts == 3 || npts == 4);
	color.set_c4(p.color);
	for (unsigned i = 0; i < npts; ++i) {pts[i]  = p[i].v;}
	if (npts == 3) pts[3] = pts[2]; // duplicate the last point so that it's valid
	get_normal(pts[0], pts[1], pts[2], normal, 1);
}


coll_tquad::coll_tquad(triangle const &t, colorRGBA const &c) {

	npts = 3;
	UNROLL_3X(pts[i_] = t.pts[i_];);
	get_normal(pts[0], pts[1], pts[2], normal, 1);
	color.set_c4(c);
}


bool tquad_t::is_valid() const {return (npts >= 3 && is_triangle_valid(pts[0], pts[1], pts[2]));}


#define UPDATE_CUBE(i) {if (pts[i][i_] < c.d[i_][0]) c.d[i_][0] = pts[i][i_]; if (pts[i][i_] > c.d[i_][1]) c.d[i_][1] = pts[i][i_];}


void tquad_t::update_bcube(cube_t &c) const {

	UNROLL_3X(UPDATE_CUBE(0));
	UNROLL_3X(UPDATE_CUBE(1));
	UNROLL_3X(UPDATE_CUBE(2));
	if (npts == 4) UNROLL_3X(UPDATE_CUBE(3));
}


cube_t tquad_t::get_bcube() const {

	cube_t c(pts[0], pts[1]);
	UNROLL_3X(UPDATE_CUBE(2));
	if (npts == 4) UNROLL_3X(UPDATE_CUBE(3));
	return c;
}


unsigned cobj_tree_base::tree_node::get_split_dim(float &max_sz, float &sval, unsigned skip_dims) const {

	unsigned dim(0);
	max_sz = 0;

	for (unsigned i = 0; i < 3; ++i) {
		if (skip_dims & (1 << i)) continue;
		float const dim_sz(d[i][1] - d[i][0]);
		assert(dim_sz >= 0.0);
		
		if (dim_sz > max_sz) {
			max_sz = dim_sz;
			dim    = i;
		}
	}
	if (max_sz > 0.0) {
		sval = get_cube_center()[dim]; // center point (mean seems to work better than median)
		assert(!(skip_dims & (1 << dim)));
	}
	return dim;
}


void cobj_tree_tquads_t::calc_node_bbox(tree_node &n) const {

	assert(n.start < n.end);
	cube_t &c(n);
	c = cube_t(X_SCENE_SIZE, -X_SCENE_SIZE, Y_SCENE_SIZE, -Y_SCENE_SIZE, czmax, czmin);

	for (unsigned i = n.start; i < n.end; ++i) { // bbox union
		tquads[i].update_bcube(c);
	}
	c.expand_by(POLY_TOLER);
}


bool cobj_tree_base::check_for_leaf(unsigned num, unsigned skip_dims) {

	if (num <= MAX_LEAF_SIZE || skip_dims == 7) { // base case
		register_leaf(num);
		return 1;
	}
	return 0;
}


bool cobj_tree_base::node_ix_mgr::check_node(unsigned &nix) const {

	tree_node const &n(nodes[nix]);

	if (!get_line_clip2(p1, dinv, n.d)) {
		assert(n.next_node_id > nix);
		nix = n.next_node_id; // failed the bbox test
		return 0;
	}
	++nix;
	return 1;
}


void cobj_tree_tquads_t::build_tree(unsigned nix, unsigned skip_dims, unsigned depth) {
	assert(nix < nodes.size());
	tree_node &n(nodes[nix]);
	calc_node_bbox(n);
	unsigned const num(n.end - n.start);
	max_depth = max(max_depth, depth);
	if (check_for_leaf(num, skip_dims)) return; // base case

	// determine split dimension and value
	float max_sz(0), sval(0);
	unsigned const dim(n.get_split_dim(max_sz, sval, skip_dims));

	if (max_sz == 0) { // can't split
		register_leaf(num);
		return;
	}
	float const sval_lo(sval+OVERLAP_AMT*max_sz), sval_hi(sval-OVERLAP_AMT*max_sz);
	unsigned pos(n.start), bin_count[3] = {0};

	// split in this dimension
	for (unsigned i = n.start; i < n.end; ++i) {
		unsigned bix(2);
		coll_tquad const &t(tquads[i]);
		float vlo(min(min(t.pts[0][dim], t.pts[1][dim]), t.pts[2][dim])), vhi(max(max(t.pts[0][dim], t.pts[1][dim]), t.pts[2][dim]));
		if (t.npts == 4) {vlo = min(vlo, t.pts[3][dim]); vhi = max(vhi, t.pts[3][dim]);}
		if (vhi <= sval_lo) bix =  (depth&1); // ends   before the split, put in bin 0
		if (vlo >= sval_hi) bix = !(depth&1); // starts after  the split, put in bin 1
		++bin_count[bix];
		temp_bins[bix].push_back(tquads[i]);
	}
	for (unsigned d = 0; d < 3; ++d) {
		for (unsigned i = 0; i < temp_bins[d].size(); ++i) {
			tquads[pos++] = temp_bins[d][i];
		}
		temp_bins[d].resize(0);
	}
	assert(pos == n.end);

	// check that dataset has been subdivided (not all in one bin)
	if (bin_count[0] == num || bin_count[1] == num || bin_count[2] == num) {
		build_tree(nix, (skip_dims | (1 << dim)), depth); // single bin, rebin with a different dim
		return;
	}

	// create child nodes and call recursively
	unsigned cur(n.start);

	for (unsigned bix = 0; bix < 3; ++bix) { // Note: this loop will invalidate the reference to 'n'
		unsigned const count(bin_count[bix]);
		if (count == 0) continue; // empty bin
		unsigned const kid((unsigned)nodes.size());
		nodes.push_back(tree_node(cur, cur+count));
		build_tree(kid, skip_dims, depth+1);
		nodes[kid].next_node_id = (unsigned)nodes.size();
		cur += count;
	}
	assert(cur == nodes[nix].end);
	nodes[nix].start = nodes[nix].end = 0; // branch node has no leaves
}


void cobj_tree_tquads_t::build_tree_top(bool verbose) {

	nodes.reserve(get_conservative_num_nodes(tquads.size()));
	nodes.push_back(tree_node(0, (unsigned)tquads.size()));
	assert(nodes.size() == 1);
	max_depth = max_leaf_count = num_leaf_nodes = 0;
	build_tree(0, 0, 0);
	nodes[0].next_node_id = (unsigned)nodes.size();
	for (unsigned i = 0; i < 3; ++i) {vector<coll_tquad>().swap(temp_bins[i]);}

	if (verbose) {
		cout << "tquads: " << tquads.size() << ", nodes: " << nodes.size()  << ", depth: "
				<< max_depth << ", max_leaf: " << max_leaf_count << ", leaf_nodes: " << num_leaf_nodes << endl;
	}
}


void cobj_tree_tquads_t::add_cobjs(coll_obj_group const &cobjs, bool verbose) {

	RESET_TIME;
	clear();
	tquads.reserve(cobjs.size()); // is this a good idea?
		
	for (coll_obj_group::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
		if (i->status != COLL_STATIC) continue;
		assert(i->type == COLL_POLYGON && i->thickness <= MIN_POLY_THICK);
		tquads.push_back(coll_tquad(*i));
	}
	build_tree_top(verbose);
	PRINT_TIME(" Cobj Tree Triangles Create (from Cobjs)");
	//exit(0); // TESTING
}


void cobj_tree_tquads_t::add_polygons(vector<polygon_t> const &polygons, bool verbose) {

	RESET_TIME;
	clear();
	tquads.reserve(polygons.size());
		
	for (vector<polygon_t>::const_iterator i = polygons.begin(); i != polygons.end(); ++i) {
		tquads.push_back(coll_tquad(*i));
	}
	build_tree_top(verbose);
	PRINT_TIME(" Cobj Tree Triangles Create (from Polygons)");
	//exit(0); // TESTING
}


bool cobj_tree_tquads_t::check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA *color, int *cindex, int ignore_cobj, bool exact) const {

	if (nodes.empty()) return 0;
	bool ret(0);
	float t(0.0), tmin(0.0), tmax(1.0);
	node_ix_mgr nixm(nodes, p1, p2);
	unsigned const num_nodes((unsigned)nodes.size());

	for (unsigned nix = 0; nix < num_nodes;) {
		tree_node const &n(nodes[nix]);
		if (!nixm.check_node(nix)) continue;

		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			// Note: we test cobj against the original (unclipped) p1 and p2 so that t is correct
			// Note: we probably don't need to return cnorm and cpos in inexact mode, but it shouldn't be too expensive to do so
			if (cindex && (int)tquads[i].cid == ignore_cobj)             continue;
			if (!tquads[i].line_int_exact(p1, p2, t, cnorm, tmin, tmax)) continue;
			if (cindex) *cindex = tquads[i].cid;
			if (color ) *color  = tquads[i].color.get_c4();
			cpos = p1 + (p2 - p1)*t;
			if (!exact) return 1; // return first hit
			nixm.dinv = vector3d(cpos - p1);
			nixm.dinv.invert(0, 1);
			tmax = t;
			ret  = 1;
		}
	}
	return ret;
}


template<unsigned NUM> bool cobj_tree_t<NUM>::create_cixs() {

	if (moving_only) {
		for (platform_cont::const_iterator i = platforms.begin(); i != platforms.end(); ++i) {
			for (vector<unsigned>::const_iterator j = i->cobjs.begin(); j != i->cobjs.end(); ++j) {
				add_cobj(*j);
			}
		}
		for (unsigned i = 0; i < falling_cobjs.size(); ++i) {
			add_cobj(falling_cobjs[i]);
		}
	}
	else if (is_dynamic && !is_static && !dynamic_ranges.empty()) { // use dynamic_ranges here
		assert(dynamic_ranges.size() & 1);

		for (unsigned r = 0; r < dynamic_ranges.size(); r += 2) {
			unsigned const end((r+1 < dynamic_ranges.size()) ? dynamic_ranges[r+1] : cobjs.size());
			for (unsigned i = dynamic_ranges[r]; i < end; ++i) add_cobj(i);
		}
	}
	else {
		bool const normal_static_mode(is_static && !occluders_only && !cubes_only);
		bool in_static_range(1);

		if (normal_static_mode) {
			cixs.reserve(cobjs.size());
			dynamic_ranges.resize(0); // recompute
		}
		for (unsigned i = 0; i < cobjs.size(); ++i) {
			if (normal_static_mode && cobjs[i].truly_static() != in_static_range) {
				dynamic_ranges.push_back(i);
				in_static_range ^= 1;
			}
			add_cobj(i);
		}
		if (normal_static_mode && in_static_range) dynamic_ranges.push_back(cobjs.size());
	}
	assert(cixs.size() < (1 << 29));
	return !cixs.empty();
}


template<unsigned NUM> void cobj_tree_t<NUM>::calc_node_bbox(tree_node &n) const {

	assert(n.start < n.end);
	n.copy_from(get_cobj(n.start));

	for (unsigned i = n.start+1; i < n.end; ++i) { // bbox union
		n.union_with_cube(get_cobj(i));
	}
}


template<unsigned NUM> void cobj_tree_t<NUM>::add_cobjs(bool verbose) {

	RESET_TIME;
	clear();
	if (!create_cixs()) return; // nothing to be done
	nodes.reserve(get_conservative_num_nodes(cixs.size()));
	nodes.push_back(tree_node(0, (unsigned)cixs.size()));
	assert(nodes.size() == 1);
	max_depth = max_leaf_count = num_leaf_nodes = 0;
	build_tree(0, 0, 0);
	nodes[0].next_node_id = (unsigned)nodes.size();
	for (unsigned i = 0; i < NUM; ++i) {vector<unsigned>().swap(temp_bins[i]);}

	if (verbose) {
		PRINT_TIME(" Cobj Tree Create");
		cout << "cobjs: " << cobjs.size() << ", leaves: " << cixs.size() << ", nodes: " << nodes.size()
				<< ", depth: " << max_depth << ", max_leaves: " << max_leaf_count << ", leaf_nodes: " << num_leaf_nodes << endl;
	}
}


// test_alpha: 0 = allow any alpha value, 1 = require alpha = 1.0, 2 = get intersected cobj with max alpha, 3 = require alpha >= MIN_SHADOW_ALPHA
template<unsigned NUM> bool cobj_tree_t<NUM>::check_coll_line(point const &p1, point const &p2, point &cpos,
	vector3d &cnorm, int &cindex, int ignore_cobj, bool exact, int test_alpha, bool skip_non_drawn) const
{
	cindex = -1;
	if (nodes.empty()) return 0;
	bool ret(0);
	float t(0.0), tmin(0.0), tmax(1.0), max_alpha(0.0);
	node_ix_mgr nixm(nodes, p1, p2);
	unsigned const num_nodes((unsigned)nodes.size());

	for (unsigned nix = 0; nix < num_nodes;) {
		tree_node const &n(nodes[nix]);
		if (!nixm.check_node(nix)) continue;

		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			// Note: we test cobj against the original (unclipped) p1 and p2 so that t is correct
			// Note: we probably don't need to return cnorm and cpos in inexact mode, but it shouldn't be too expensive to do so
			if ((int)cixs[i] == ignore_cobj)  continue;
			coll_obj const &c(get_cobj(i));
			if (!obj_ok(c))                   continue;
			if (skip_non_drawn && !c.cp.might_be_drawn())               continue;
			if (test_alpha == 1 && c.is_semi_trans())                   continue; // semi-transparent, can see through
			if (test_alpha == 2 && c.cp.color.alpha <= max_alpha)       continue; // lower alpha than an earlier object
			if (test_alpha == 3 && c.cp.color.alpha < MIN_SHADOW_ALPHA) continue; // less than min alpha
			if (!c.line_int_exact(p1, p2, t, cnorm, tmin, tmax))        continue;
			cindex = cixs[i];
			cpos   = p1 + (p2 - p1)*t;
			//if (c.type == COLL_POLYGON && dot_product((p2 - p1), c.norm) < 0.0) {} // back-facing polygon test
			if (!exact && test_alpha != 2) return 1; // return first hit
			max_alpha = c.cp.color.alpha; // we need all intersections to find the max alpha
			nixm.dinv = vector3d(cpos - p1);
			nixm.dinv.invert(0, 1);
			tmax = t;
			ret  = 1;
		}
	}
	return ret;
}


template<unsigned NUM> void cobj_tree_t<NUM>::get_intersecting_cobjs(cube_t const &cube, vector<unsigned> &cobjs,
	int ignore_cobj, float toler, bool check_ccounter, int id_for_cobj_int) const
{
	unsigned const num_nodes((unsigned)nodes.size());

	for (unsigned nix = 0; nix < num_nodes;) {
		tree_node const &n(nodes[nix]);
		assert(n.start <= n.end);

		if (!cube.intersects(n, toler)) {
			assert(n.next_node_id > nix);
			nix = n.next_node_id; // failed the bbox test
			continue;
		}
		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			if ((int)cixs[i] == ignore_cobj) continue;
			coll_obj const &c(get_cobj(i));
			if (check_ccounter && c.counter == cobj_counter) continue;
			// get_intersecting_cobjs_tree() calls this on both static and moving cobj_trees, so we want to check to make sure we don't double include it
			if (!moving_only && c.maybe_is_moving())         continue;
			if (!obj_ok(c) || !cube.intersects(c, toler))    continue;
			if (id_for_cobj_int >= 0 && coll_objects[id_for_cobj_int].intersects_cobj(c, toler) != 1) continue;
			cobjs.push_back(cixs[i]);
		}
		++nix;
	}
}


template<unsigned NUM> bool cobj_tree_t<NUM>::is_cobj_contained(point const &p1, point const &p2, point const &viewer,
	point const *const pts, unsigned npts, int ignore_cobj, int &cobj) const
{
	if (nodes.empty()) return 0;
	node_ix_mgr nixm(nodes, p1, p2);
	unsigned const num_nodes((unsigned)nodes.size());

	for (unsigned nix = 0; nix < num_nodes;) {
		tree_node const &n(nodes[nix]);
		if (!nixm.check_node(nix)) continue;

		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			if ((int)cixs[i] == ignore_cobj) continue;
			coll_obj const &c(get_cobj(i));
				
			if (obj_ok(c) && c.is_occluder() && c.intersects_all_pts(viewer, pts, npts)) {
				cobj = cixs[i];
				return 1;
			}
		}
	}
	return 0;
}


template<unsigned NUM> void cobj_tree_t<NUM>::get_coll_line_cobjs(point const &pos1, point const &pos2,
	int ignore_cobj, vector<int> *cobjs, cobj_query_callback *cqc, bool occlude) const
{
	assert(cobjs || cqc);
	if (nodes.empty()) return;
	node_ix_mgr nixm(nodes, pos1, pos2);
	unsigned const num_nodes((unsigned)nodes.size());

	for (unsigned nix = 0; nix < num_nodes;) {
		tree_node const &n(nodes[nix]);
		if (!nixm.check_node(nix)) continue;
			
		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			if ((int)cixs[i] == ignore_cobj) continue;
			coll_obj const &c(get_cobj(i));
			if (!obj_ok(c)) continue;
			if (occlude && !(c.is_big_occluder() && check_line_clip_expand(pos1, pos2, c.d, GET_OCC_EXPAND))) continue;
			if (cobjs) cobjs->push_back(cixs[i]);
			if (cqc  ) cqc->register_cobj(c);
		}
	}
}


// Note: actually, this only returns sphere intersection candidates
template<unsigned NUM> void cobj_tree_t<NUM>::get_coll_sphere_cobjs(point const &center, float radius, int ignore_cobj, vert_coll_detector &vcd) const {

	if (nodes.empty()) return;
	unsigned const num_nodes((unsigned)nodes.size());
	cube_t bcube(center, center);
	bcube.expand_by(radius);

	for (unsigned nix = 0; nix < num_nodes;) {
		tree_node const &n(nodes[nix]);

		if (!n.intersects(bcube)/* && !sphere_cube_intersect(center, radius, n)*/) {
			assert(n.next_node_id > nix);
			nix = n.next_node_id; // failed the bbox test
			continue;
		}
		++nix;
		
		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			if ((int)cixs[i] != ignore_cobj && get_cobj(i).intersects(bcube)) vcd.check_cobj(cixs[i]);
		}
	}
}


// BSP Tree/KD-Tree (left, right, mid) kids
template <> void cobj_tree_t<3>::build_tree(unsigned nix, unsigned skip_dims, unsigned depth) {
	
	assert(nix < nodes.size());
	tree_node &n(nodes[nix]);
	calc_node_bbox(n);
	unsigned const num(n.end - n.start);
	max_depth = max(max_depth, depth);
	if (check_for_leaf(num, skip_dims)) return; // base case
	
	// determine split dimension and value
	float max_sz(0), sval(0);
	unsigned const dim(n.get_split_dim(max_sz, sval, skip_dims));

	if (max_sz == 0) { // can't split
		register_leaf(num);
		return;
	}
	float const sval_lo(sval+OVERLAP_AMT*max_sz), sval_hi(sval-OVERLAP_AMT*max_sz);
	unsigned pos(n.start), bin_count[3] = {0};

	// split in this dimension: use upper 2 bits of cixs for storing bin index
	for (unsigned i = n.start; i < n.end; ++i) {
		unsigned bix(2);
		float const *vals(get_cobj(i).d[dim]);
		assert(vals[0] <= vals[1]);
		if (vals[1] <= sval_lo) bix =  (depth&1); // ends   before the split, put in bin 0
		if (vals[0] >= sval_hi) bix = !(depth&1); // starts after  the split, put in bin 1
		++bin_count[bix];
		temp_bins[bix].push_back(cixs[i]);
	}
	for (unsigned d = 0; d < 3; ++d) {
		for (unsigned i = 0; i < temp_bins[d].size(); ++i) {
			cixs[pos++] = temp_bins[d][i];
		}
		temp_bins[d].resize(0);
	}
	assert(pos == n.end);

	// check that dataset has been subdivided (not all in one bin)
	if (bin_count[0] == num || bin_count[1] == num || bin_count[2] == num) {
		build_tree(nix, (skip_dims | (1 << dim)), depth); // single bin, rebin with a different dim
		return;
	}

	// create child nodes and call recursively
	unsigned cur(n.start);

	for (unsigned bix = 0; bix < 3; ++bix) { // Note: this loop will invalidate the reference to 'n'
		unsigned const count(bin_count[bix]);
		if (count == 0) continue; // empty bin
		unsigned const kid((unsigned)nodes.size());
		nodes.push_back(tree_node(cur, cur+count));
		build_tree(kid, skip_dims, depth+1);
		nodes[kid].next_node_id = (unsigned)nodes.size();
		cur += count;
	}
	assert(cur == nodes[nix].end);
	nodes[nix].start = nodes[nix].end = 0; // branch node has no leaves
}


// OctTree
template <> void cobj_tree_t<8>::build_tree(unsigned nix, unsigned skip_dims, unsigned depth) {

	assert(nix < nodes.size());
	tree_node &n(nodes[nix]);
	calc_node_bbox(n);
	unsigned const num(n.end - n.start);
	max_depth = max(max_depth, depth);
	if (check_for_leaf(num, skip_dims)) return; // base case
	
	// determine split values
	point const sval(n.get_cube_center()); // center point
	unsigned pos(n.start), bin_count[8] = {0};

	// split in this dimension: use upper 3 bits of cixs for storing bin index
	for (unsigned i = n.start; i < n.end; ++i) {
		point const center(get_cobj(i).get_cube_center());
		unsigned bix(0);
		UNROLL_3X(if (center[i_] > sval[i_]) bix |= (1 << i_);)
		++bin_count[bix];
		temp_bins[bix].push_back(cixs[i]);
	}
	for (unsigned d = 0; d < 8; ++d) {
		for (unsigned i = 0; i < temp_bins[d].size(); ++i) {
			cixs[pos++] = temp_bins[d][i];
		}
		temp_bins[d].resize(0);
	}
	assert(pos == n.end);

	// create child nodes and call recursively
	unsigned cur(n.start);

	for (unsigned bix = 0; bix < 8; ++bix) { // Note: this loop will invalidate the reference to 'n'
		unsigned const count(bin_count[bix]);
		if (count == 0) continue; // empty bin
		unsigned const kid((unsigned)nodes.size());
		nodes.push_back(tree_node(cur, cur+count));
		build_tree(kid, ((count == num) ? 7 : 0), depth+1); // if all in one bin, make that bin a leaf
		nodes[kid].next_node_id = (unsigned)nodes.size();
		cur += count;
	}
	assert(cur == nodes[nix].end);
	nodes[nix].start = nodes[nix].end = 0; // branch node has no leaves
}


cobj_tree_type cobj_tree_static (coll_objects, 1, 0, 0, 0);
cobj_tree_type cobj_tree_dynamic(coll_objects, 0, 1, 0, 0);
cobj_tree_type cobj_tree_occlude(coll_objects, 1, 0, 1, 0);
cobj_tree_type cobj_tree_moving (coll_objects, 1, 0, 0, 1);
cobj_tree_tquads_t cobj_tree_triangles;


cobj_tree_type &get_tree(bool dynamic) {
	return (dynamic ? cobj_tree_dynamic : cobj_tree_static);
}

void build_cobj_tree(bool dynamic, bool verbose) {

	get_tree(dynamic).add_cobjs(verbose);
	if (!dynamic) cobj_tree_occlude.add_cobjs(verbose);
	//if (!dynamic) cobj_tree_triangles.add_cobjs(coll_objects, verbose);
}

void build_moving_cobj_tree() {
	cobj_tree_moving.add_cobjs(0);
}

// can use with ray trace lighting, snow collision?, maybe water reflections
bool check_coll_line_exact_tree(point const &p1, point const &p2, point &cpos,
								vector3d &cnorm, int &cindex, int ignore_cobj, bool dynamic, int test_alpha, bool skip_non_drawn)
{
	//return cobj_tree_triangles.check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 1);
	return get_tree(dynamic).check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 1, test_alpha, skip_non_drawn);
}

// can use with snow shadows, grass shadows, tree leaf shadows
bool check_coll_line_tree(point const &p1, point const &p2, int &cindex, int ignore_cobj, bool dynamic, int test_alpha, bool skip_non_drawn) {

	vector3d cnorm; // unused
	point cpos; // unused
	return get_tree(dynamic).check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 0, test_alpha, skip_non_drawn);
}


void get_intersecting_cobjs_tree(cube_t const &cube, vector<unsigned> &cobjs, int ignore_cobj, float toler,
	bool dynamic, bool check_ccounter, int id_for_cobj_int)
{
	get_tree(dynamic).get_intersecting_cobjs(cube, cobjs, ignore_cobj, toler, check_ccounter, id_for_cobj_int);
	cobj_tree_moving.get_intersecting_cobjs (cube, cobjs, ignore_cobj, toler, check_ccounter, id_for_cobj_int);
}


bool cobj_contained_tree(point const &p1, point const &p2, point const &viewer, point const *const pts,
	unsigned npts, int ignore_cobj, int &cobj)
{
	return cobj_tree_occlude.is_cobj_contained(p1, p2, viewer, pts, npts, ignore_cobj, cobj);
}


void get_coll_line_cobjs_tree(point const &pos1, point const &pos2, int ignore_cobj,
	vector<int> *cobjs, cobj_query_callback *cqc, bool dynamic, bool occlude)
{
	(occlude ? cobj_tree_occlude : get_tree(dynamic)).get_coll_line_cobjs(pos1, pos2, ignore_cobj, cobjs, cqc, occlude);
}


void get_coll_sphere_cobjs_tree(point const &center, float radius, int cobj, vert_coll_detector &vcd, bool dynamic) {
	get_tree(dynamic).get_coll_sphere_cobjs(center, radius, cobj, vcd);
}


bool have_occluders() {
	return !cobj_tree_occlude.is_empty();
}


