// 3D World - OpenGL CS184 Computer Graphics Project - collision detection BSP/KD/Oct Tree
// by Frank Gennari
// 10/16/10

#include "3DWorld.h"
#include "cobj_bsp_tree.h"


unsigned const MAX_LEAF_SIZE = 2;
float const POLY_TOLER       = 1.0E-6;
float const OVERLAP_AMT      = 0.02;


extern bool mt_cobj_tree_build;
extern int display_mode, frame_counter, cobj_counter, begin_motion;
extern coll_obj_group coll_objects;
extern vector<unsigned> falling_cobjs;
extern platform_cont platforms;


// *** coll_tquad / tquad_t ***


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
	update_normal();
}


coll_tquad::coll_tquad(triangle const &t, colorRGBA const &c) {

	npts = 3;
	UNROLL_3X(pts[i_] = t.pts[i_];);
	update_normal();
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


// *** cobj_tree_base ***


bool cobj_tree_base::get_root_bcube(cube_t &bc) const {
	
	if (nodes.empty()) {return 0;}
	bc = nodes[0];
	return 1;
}


bool cobj_tree_base::check_for_leaf(unsigned num, unsigned skip_dims) {

	if (num <= MAX_LEAF_SIZE || skip_dims == 7) { // base case
		register_leaf(num);
		return 1;
	}
	return 0;
}


// performance critical
template<bool xneg, bool yneg, bool zneg> bool get_line_clip(point const &p1, vector3d const &dinv, float const d[3][2]) {

	float tmin(0.0), tmax(1.0);
	float const t1((d[0][xneg] - p1.x)*dinv.x), t2((d[0][!xneg] - p1.x)*dinv.x);
	if (t2 < tmax) {tmax = t2;} if (t1 > tmin) {tmin = t1;}

	if (tmin < tmax) {
		float const t1((d[1][yneg] - p1.y)*dinv.y), t2((d[1][!yneg] - p1.y)*dinv.y);
		if (t2 < tmax) {tmax = t2;} if (t1 > tmin) {tmin = t1;}

		if (tmin < tmax) {
			float const t1((d[2][zneg] - p1.z)*dinv.z), t2((d[2][!zneg] - p1.z)*dinv.z);
			if (t2 < tmax) {tmax = t2;} if (t1 > tmin) {tmin = t1;}
			if (tmin < tmax) {return 1;}
		}
	}
	return 0;
}

cobj_tree_base::node_ix_mgr::node_ix_mgr(vector<tree_node> const &nodes_, point const &p1_, point const &p2_)
			: nodes(nodes_), p1(p1_), p2(p2_), dinv(p2 - p1)
{
	dinv.invert();
	if (dinv.x < 0.0) {
		if (dinv.y < 0.0) {
			if (dinv.z < 0.0) {get_line_clip_func = get_line_clip<1,1,1>;}
			else              {get_line_clip_func = get_line_clip<1,1,0>;}} else {
			if (dinv.z < 0.0) {get_line_clip_func = get_line_clip<1,0,1>;}
			else              {get_line_clip_func = get_line_clip<1,0,0>;}}} else {
		if (dinv.y < 0.0) {
			if (dinv.z < 0.0) {get_line_clip_func = get_line_clip<0,1,1>;}
			else              {get_line_clip_func = get_line_clip<0,1,0>;}} else {
			if (dinv.z < 0.0) {get_line_clip_func = get_line_clip<0,0,1>;}
			else              {get_line_clip_func = get_line_clip<0,0,0>;}}}
}

bool cobj_tree_base::node_ix_mgr::check_node(unsigned &nix) const {

	tree_node const &n(nodes[nix]);

	if (!get_line_clip_func(p1, dinv, n.d)) {
		assert(n.next_node_id > nix);
		nix = n.next_node_id; // failed the bbox test
		return 0;
	}
	++nix;
	return 1;
}


// *** cobj_tree_simple_type_t ***


inline float get_vlo(coll_tquad const &t, unsigned dim) {
	float vlo(min(min(t.pts[0][dim], t.pts[1][dim]), t.pts[2][dim]));
	if (t.npts == 4) {vlo = min(vlo, t.pts[3][dim]);}
	return vlo;
}
inline float get_vhi(coll_tquad const &t, unsigned dim) {
	float vhi(max(max(t.pts[0][dim], t.pts[1][dim]), t.pts[2][dim]));
	if (t.npts == 4) {vhi = max(vhi, t.pts[3][dim]);}
	return vhi;
}

inline float get_vlo(sphere_t const &s, unsigned dim) {
	return (s.pos[dim] - s.radius);
}
inline float get_vhi(sphere_t const &s, unsigned dim) {
	return (s.pos[dim] + s.radius);
}


template<typename T> void cobj_tree_simple_type_t<T>::build_tree(unsigned nix, unsigned skip_dims, unsigned depth) {

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
	unsigned pos(n.start), bin_count[3];

	// split in this dimension
	for (unsigned i = n.start; i < n.end; ++i) {
		unsigned bix(2);
		T const &obj(objects[i]);
		if (get_vhi(obj, dim) <= sval_lo) {bix =  (depth&1);} // ends   before the split, put in bin 0
		if (get_vlo(obj, dim) >= sval_hi) {bix = !(depth&1);} // starts after  the split, put in bin 1
		temp_bins[bix].push_back(obj);
	}
	for (unsigned d = 0; d < 3; ++d) {
		bin_count[d] = temp_bins[d].size();
		for (unsigned i = 0; i < bin_count[d]; ++i) {objects[pos++] = temp_bins[d][i];}
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


template<typename T> void cobj_tree_simple_type_t<T>::build_tree_top(bool verbose) {

	nodes.reserve(get_conservative_num_nodes(objects.size()));
	nodes.push_back(tree_node(0, (unsigned)objects.size()));
	assert(nodes.size() == 1);
	max_depth = max_leaf_count = num_leaf_nodes = 0;
	if (!objects.empty()) {build_tree(0, 0, 0);}
	nodes[0].next_node_id = (unsigned)nodes.size();
	for (unsigned i = 0; i < 3; ++i) {vector<T>().swap(temp_bins[i]);}

	if (verbose) {
		cout << "objects: " << objects.size() << ", nodes: " << nodes.size()  << ", depth: "
			 << max_depth << ", max_leaf: " << max_leaf_count << ", leaf_nodes: " << num_leaf_nodes << endl;
	}
}


// *** cobj_tree_tquads_t ***


void cobj_tree_tquads_t::calc_node_bbox(tree_node &n) const {

	assert(n.start < n.end);
	cube_t &c(n);
	c = cube_t(X_SCENE_SIZE, -X_SCENE_SIZE, Y_SCENE_SIZE, -Y_SCENE_SIZE, czmax, czmin);
	for (unsigned i = n.start; i < n.end; ++i) {objects[i].update_bcube(c);} // bbox union
	c.expand_by(POLY_TOLER);
}


void cobj_tree_tquads_t::add_cobjs(coll_obj_group const &cobjs, bool verbose) {

	RESET_TIME;
	clear();
	objects.reserve(cobjs.size()); // is this a good idea?
		
	for (coll_obj_group::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
		if (i->status != COLL_STATIC) continue;
		assert(i->type == COLL_POLYGON && i->thickness <= MIN_POLY_THICK);
		objects.push_back(coll_tquad(*i));
	}
	build_tree_top(verbose);
	PRINT_TIME(" Cobj Tree Triangles Create (from Cobjs)");
}


void cobj_tree_tquads_t::add_polygons(vector<polygon_t> const &polygons, bool verbose) {

	RESET_TIME;
	clear();
	objects.reserve(polygons.size());
	for (vector<polygon_t>::const_iterator i = polygons.begin(); i != polygons.end(); ++i) {objects.push_back(coll_tquad(*i));}
	build_tree_top(verbose);
	PRINT_TIME(" Cobj Tree Triangles Create (from Polygons)");
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
			if (cindex && (int)objects[i].cid == ignore_cobj)             continue;
			if (!objects[i].line_int_exact(p1, p2, t, cnorm, tmin, tmax)) continue;
			if (cindex) *cindex = objects[i].cid;
			if (color ) *color  = objects[i].color.get_c4();
			cpos = p1 + (p2 - p1)*t;
			if (!exact) return 1; // return first hit
			nixm.dinv = vector3d(cpos - p1);
			nixm.dinv.invert();
			tmax = t;
			ret  = 1;
		}
	}
	return ret;
}


// *** cobj_tree_sphere_t ***


void cobj_tree_sphere_t::calc_node_bbox(tree_node &n) const {

	assert(n.start < n.end);
	cube_t &c(n);

	for (unsigned i = n.start; i < n.end; ++i) { // bbox union
		if (i == n.start) {c.set_from_sphere(objects[i]);} else {c.union_with_sphere(objects[i]);}
	}
}


void cobj_tree_sphere_t::add_spheres(vector<sphere_with_id_t> &spheres_, bool verbose) {

	clear();
	objects.swap(spheres_); // copy, destroy input
	build_tree_top(verbose);
}


void cobj_tree_sphere_t::get_ids_int_sphere(point const &center, float radius, vector<unsigned> &ids) const {

	if (objects.empty()) return;
	unsigned const num_nodes((unsigned)nodes.size());

	for (unsigned nix = 0; nix < num_nodes;) {
		tree_node const &n(nodes[nix]);
		assert(n.start <= n.end);

		if (!sphere_cube_intersect(center, radius, n)) {
			assert(n.next_node_id > nix);
			nix = n.next_node_id; // failed the bounding sphere test
			continue;
		}
		for (unsigned i = n.start; i < n.end; ++i) { // check leaves
			if (dist_less_than(center, objects[i].pos, (radius + objects[i].radius))) {ids.push_back(objects[i].id);}
		}
		++nix;
	}
}


// *** cobj_bvh_tree ***


bool cobj_bvh_tree::create_cixs() {

	if (is_dynamic && !is_static) { // use dynamic_ids
		for (cobj_id_set_t::const_iterator i = cobjs->dynamic_ids.begin(); i != cobjs->dynamic_ids.end(); ++i) {
			assert(*i < cobjs->size());
			assert((*cobjs)[*i].status == COLL_DYNAMIC);
			add_cobj(*i);
		}
	}
	else {
		if (is_static && !occluders_only && !cubes_only) {cixs.reserve(cobjs->size());} // normal static mode
		for (unsigned i = 0; i < cobjs->size(); ++i) {add_cobj(i);}
	}
	assert(cixs.size() < (1 << 29));
	return !cixs.empty();
}


void cobj_bvh_tree::calc_node_bbox(tree_node &n) const {

	assert(n.start < n.end);
	n.copy_from(get_cobj(n.start));
	for (unsigned i = n.start+1; i < n.end; ++i) {n.union_with_cube(get_cobj(i));} // bbox union
}


void cobj_bvh_tree::clear() {

	cobj_tree_base::clear();
	cixs.resize(0);
}


void cobj_bvh_tree::add_cobjs(bool verbose) {

	RESET_TIME;
	clear();
	if (!create_cixs()) return; // nothing to be done
	bool const do_mt_build(mt_cobj_tree_build && cixs.size() > 10000);
	build_tree_from_cixs(do_mt_build);

	if (verbose) {
		PRINT_TIME(" Cobj Tree Create");
		cout << "cobjs: " << cobjs->size() << ", leaves: " << cixs.size() << ", nodes: " << nodes.size()
				<< ", depth: " << max_depth << ", max_leaves: " << max_leaf_count << ", leaf_nodes: " << num_leaf_nodes << endl;
	}
}


// to be called from within add_cobjs() or after a call to add_cobj_ids()
void cobj_bvh_tree::build_tree_from_cixs(bool do_mt_build) {

	max_depth = max_leaf_count = num_leaf_nodes = 0;
	nodes.resize(get_conservative_num_nodes(cixs.size()) + 64*do_mt_build); // add 8 extra nodes for each of 8 top level splits
	unsigned const root(0);
	nodes[root] = tree_node(0, (unsigned)cixs.size());

	if (do_mt_build) { // 2x faster build time, 10% slower traversal
		build_tree_top_level_omp();
	}
	else {
		per_thread_data ptd(1, nodes.size(), 1);
		build_tree(root, 0, 0, ptd);
		nodes.resize(ptd.get_next_node_ix());
	}
	nodes[root].next_node_id = (unsigned)nodes.size();
}


// test_alpha: 0 = allow any alpha value, 1 = require alpha = 1.0, 2 = get intersected cobj with max alpha, 3 = require alpha >= MIN_SHADOW_ALPHA
bool cobj_bvh_tree::check_coll_line(point const &p1, point const &p2, point &cpos,
	vector3d &cnorm, int &cindex, int ignore_cobj, bool exact, int test_alpha, bool skip_non_drawn) const
{
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
			nixm.dinv.invert();
			tmax = t;
			ret  = 1;
		}
	}
	return ret;
}


void cobj_bvh_tree::get_intersecting_cobjs(cube_t const &cube, vector<unsigned> &cobjs,
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
			if (!obj_ok(c) || !cube.intersects(c, toler))    continue;
			if (id_for_cobj_int >= 0 && coll_objects[id_for_cobj_int].intersects_cobj(c, toler) != 1) continue;
			cobjs.push_back(cixs[i]);
		}
		++nix;
	}
}


bool cobj_bvh_tree::is_cobj_contained(point const &p1, point const &p2, point const &viewer,
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


void cobj_bvh_tree::get_coll_line_cobjs(point const &pos1, point const &pos2,
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
void cobj_bvh_tree::get_coll_sphere_cobjs(point const &center, float radius, int ignore_cobj, vert_coll_detector &vcd) const {

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


void cobj_bvh_tree::build_tree_top_level_omp() { // single octtree level

	vector<unsigned> top_temp_bins[8];
	unsigned const nix(0);
	tree_node &n(nodes[nix]);
	unsigned const num(n.end - n.start);

	// calculate bbox and determine mean values
	point sval(all_zeros);
	n.copy_from(get_cobj(n.start));

	for (unsigned i = n.start; i < n.end; ++i) {
		coll_obj const &cobj(get_cobj(i));
		n.union_with_cube(cobj);
		sval += cobj.get_cube_center();
	}
	sval /= num;
	unsigned pos(n.start);

	// split in this dimension
	for (unsigned i = n.start; i < n.end; ++i) {
		point const center(get_cobj(i).get_cube_center());
		unsigned bix(0);
		UNROLL_3X(if (center[i_] > sval[i_]) bix |= (1 << i_);)
		top_temp_bins[bix].push_back(cixs[i]);
	}
	for (unsigned d = 0; d < 8; ++d) {
		memcpy(&cixs[pos], &top_temp_bins[d].front(), top_temp_bins[d].size()*sizeof(unsigned));
		pos += top_temp_bins[d].size();
	}
	assert(pos == n.end);

	// create child nodes and call recursively
	unsigned cur(n.start), cur_nix(1);
	unsigned curs[8], cur_nixs[8];
	
	for (int bix = 0; bix < 8; ++bix) {
		unsigned const count(top_temp_bins[bix].size());
		if (count == 0) continue; // empty bin
		curs[bix]     = cur;
		cur_nixs[bix] = cur_nix;
		cur     += count;
		cur_nix += get_conservative_num_nodes(count);
		assert(cur_nix <= nodes.size());
	}

	#pragma omp parallel for schedule(static,1)
	for (int bix = 0; bix < 8; ++bix) {
		unsigned const count(top_temp_bins[bix].size());
		if (count == 0) continue; // empty bin
		unsigned const kid(cur_nixs[bix]), alloc_sz(get_conservative_num_nodes(count)), end_nix(cur_nixs[bix] + alloc_sz);
		nodes[kid] = tree_node(curs[bix], curs[bix]+count);
		per_thread_data ptd(cur_nixs[bix]+1, end_nix, 0);
		build_tree(kid, ((count == num) ? 7 : 0), 1, ptd); // if all in one bin, make that bin a leaf
		unsigned const next_kid(ptd.get_next_node_ix());
		assert(next_kid <= end_nix);
		if (next_kid < end_nix) {nodes[next_kid].next_node_id = end_nix;} // close the gap of unused nodes
		nodes[kid].next_node_id = end_nix;
	}
	nodes.resize(cur_nix);
	assert(cur == n.end);
	n.start = n.end = 0; // branch node has no leaves
}


// BVH (left, right, mid) kids
void cobj_bvh_tree::build_tree(unsigned nix, unsigned skip_dims, unsigned depth, per_thread_data &ptd) {
	
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
	unsigned pos(n.start), bin_count[3];

	// split in this dimension: use upper 2 bits of cixs for storing bin index
	for (unsigned i = n.start; i < n.end; ++i) {
		unsigned bix(2);
		float const *vals(get_cobj(i).d[dim]);
		assert(vals[0] <= vals[1]);
		if (vals[1] <= sval_lo) {bix =  (depth&1);} // ends   before the split, put in bin 0
		if (vals[0] >= sval_hi) {bix = !(depth&1);} // starts after  the split, put in bin 1
		if (bix == 0) {cixs[pos++] = cixs[i];} else {ptd.temp_bins[bix].push_back(cixs[i]);}
	}
	bin_count[0] = (pos - n.start);

	for (unsigned d = 1; d < 3; ++d) {
		bin_count[d] = ptd.temp_bins[d].size();
		for (unsigned i = 0; i < bin_count[d]; ++i) {cixs[pos++] = ptd.temp_bins[d][i];}
		ptd.temp_bins[d].resize(0);
	}
	assert(pos == n.end);

	// check that dataset has been subdivided (not all in one bin)
	if (bin_count[0] == num || bin_count[1] == num || bin_count[2] == num) {
		build_tree(nix, (skip_dims | (1 << dim)), depth, ptd); // single bin, rebin with a different dim
		return;
	}
	// create child nodes and call recursively
	unsigned cur(n.start);

	for (unsigned bix = 0; bix < 3; ++bix) {
		unsigned const count(bin_count[bix]);
		if (count == 0) continue; // empty bin
		unsigned const kid(ptd.get_next_node_ix());
		ptd.increment_node_ix();

		if (ptd.at_node_end()) {
			assert(ptd.can_be_resized);
			unsigned const old_nodes_size(nodes.size());
			nodes.resize(5*old_nodes_size/4); // increase by 25% (will invalidate n reference)
			cout << "Warning: Resizing cobj_bvh_tree nodes from " << old_nodes_size << " to " << nodes.size() << endl;
			ptd.advance_end_range(nodes.size());
		}
		nodes[kid] = tree_node(cur, cur+count);
		build_tree(kid, skip_dims, depth+1, ptd);
		nodes[kid].next_node_id = ptd.get_next_node_ix();
		cur += count;
	}
	assert(cur == nodes[nix].end);
	nodes[nix].start = nodes[nix].end = 0; // branch node has no leaves
}


// is_static is_dynamic occluders_only cubes_only inc_voxel_cobjs
cobj_bvh_tree cobj_tree_static (&coll_objects, 1, 0, 0, 0, 0); // does not include voxels
cobj_bvh_tree cobj_tree_dynamic(&coll_objects, 0, 1, 0, 0, 0);
cobj_bvh_tree cobj_tree_occlude(&coll_objects, 1, 0, 1, 0, 0);
cobj_bvh_tree cobj_tree_static_moving(&coll_objects, 1, 0, 0, 0, 0);
cobj_tree_tquads_t cobj_tree_triangles;


cobj_bvh_tree &get_tree(bool dynamic) {
	return (dynamic ? cobj_tree_dynamic : cobj_tree_static);
}

void build_cobj_tree(bool dynamic, bool verbose) {

	get_tree(dynamic).add_cobjs(verbose);
	
	if (!dynamic) { // static
		cobj_tree_occlude.add_cobjs(verbose);
		//cobj_tree_triangles.add_cobjs(coll_objects, verbose);
	}
	else { // dynamic
		cobj_tree_static_moving.clear();
		vector<unsigned> moving_cids(falling_cobjs);

		for (platform_cont::const_iterator i = platforms.begin(); i != platforms.end(); ++i) {
			copy(i->cobjs.begin(), i->cobjs.end(), back_inserter(moving_cids));
		}
		if (!moving_cids.empty()) {
			cobj_tree_static_moving.add_cobj_ids(moving_cids);
			cobj_tree_static_moving.build_tree_from_cixs(0);
		}
	}
}

// can use with ray trace lighting, snow collision?, maybe water reflections
bool check_coll_line_exact_tree(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex,
	int ignore_cobj, bool dynamic, int test_alpha, bool skip_non_drawn, bool include_voxels)
{
	cindex = -1;
	//return cobj_tree_triangles.check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 1);
	bool ret(get_tree(dynamic).check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 1, test_alpha, skip_non_drawn));
	if (!dynamic) {ret |= cobj_tree_static_moving.check_coll_line(p1, (ret ? cpos : p2), cpos, cnorm, cindex, ignore_cobj, 1, test_alpha, skip_non_drawn);}
	if (!dynamic && include_voxels) {ret |= check_voxel_coll_line(p1, (ret ? cpos : p2), cpos, cnorm, cindex, ignore_cobj, 1);}
	return ret;
}

// can use with snow shadows, grass shadows, tree leaf shadows
bool check_coll_line_tree(point const &p1, point const &p2, int &cindex, int ignore_cobj,
	bool dynamic, int test_alpha, bool skip_non_drawn, bool include_voxels)
{
	vector3d cnorm; // unused
	point cpos; // unused
	cindex = -1;
	if (get_tree(dynamic).check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 0, test_alpha, skip_non_drawn)) return 1;
	if (!dynamic && cobj_tree_static_moving.check_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 0, test_alpha, skip_non_drawn)) return 1;
	if (!dynamic && include_voxels && check_voxel_coll_line(p1, p2, cpos, cnorm, cindex, ignore_cobj, 0)) return 1;
	return 0;
}

// used in destroy_cobj for cobj destroy/modification and connected/anchoring tests
void get_intersecting_cobjs_tree(cube_t const &cube, vector<unsigned> &cobjs, int ignore_cobj, float toler,
	bool dynamic, bool check_ccounter, int id_for_cobj_int)
{
	get_tree(dynamic).get_intersecting_cobjs(cube, cobjs, ignore_cobj, toler, check_ccounter, id_for_cobj_int);
	if (!dynamic) {cobj_tree_static_moving.get_intersecting_cobjs(cube, cobjs, ignore_cobj, toler, check_ccounter, id_for_cobj_int);}
}

// used in cobj_contained_ref() for grass occlusion
bool cobj_contained_tree(point const &p1, point const &p2, point const &viewer, point const *const pts,
	unsigned npts, int ignore_cobj, int &cobj)
{
	return cobj_tree_occlude.is_cobj_contained(p1, p2, viewer, pts, npts, ignore_cobj, cobj);
}

// used in get_occluders() for occlusion culling
void get_coll_line_cobjs_tree(point const &pos1, point const &pos2, int ignore_cobj,
	vector<int> *cobjs, cobj_query_callback *cqc, bool dynamic, bool occlude)
{
	(occlude ? cobj_tree_occlude : get_tree(dynamic)) .get_coll_line_cobjs(pos1, pos2, ignore_cobj, cobjs, cqc, occlude);
	if (!dynamic && !occlude) {cobj_tree_static_moving.get_coll_line_cobjs(pos1, pos2, ignore_cobj, cobjs, cqc, occlude);}
}

// used in vert_coll_detector for object collision detection
void get_coll_sphere_cobjs_tree(point const &center, float radius, int cobj, vert_coll_detector &vcd, bool dynamic) {
	get_tree(dynamic).get_coll_sphere_cobjs(center, radius, cobj, vcd);
	if (!dynamic) {cobj_tree_static_moving.get_coll_sphere_cobjs(center, radius, cobj, vcd);}
	if (!dynamic) {get_voxel_coll_sphere_cobjs(center, radius, cobj, vcd);}
}


bool have_occluders() {
	return !cobj_tree_occlude.is_empty();
}


