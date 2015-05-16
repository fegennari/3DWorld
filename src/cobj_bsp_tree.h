// 3D World Collision Object BSP/KD/Oct Tree
// by Frank Gennari
// 11/15/11
#ifndef _COBJ_BSP_TREE_H_
#define _COBJ_BSP_TREE_H_

#include "physics_objects.h"


class cobj_tree_base {

protected:
	struct tree_node : public cube_t { // size = 36
		unsigned start, end; // index into cixs for leaves
		unsigned next_node_id;

		tree_node(unsigned s=0, unsigned e=0) : start(s), end(e), next_node_id(0) {
			UNROLL_3X(d[i_][0] = d[i_][1] = 0.0;)
		}
		tree_node(unsigned s, unsigned e, cube_t const &cube) : cube_t(cube), start(s), end(e), next_node_id(0) {}
	};

	vector<tree_node> nodes;
	unsigned max_depth, max_leaf_count, num_leaf_nodes;

	inline void register_leaf(unsigned num) {
		++num_leaf_nodes;
		max_leaf_count = max(max_leaf_count, num);
	}
	bool check_for_leaf(unsigned num, unsigned skip_dims);
	unsigned get_conservative_num_nodes(unsigned num) const {return (3*num/2 + 8);}

	struct node_ix_mgr {
		point const p1, p2;
		vector3d dinv;
		vector<tree_node> const &nodes;

		node_ix_mgr(vector<tree_node> const &nodes_, point const &p1_, point const &p2_);
		bool check_node(unsigned &nix) const;
		bool (* get_line_clip_func) (point const &p1, vector3d const &dinv, float const d[3][2]); // function pointer
	};

public:
	cobj_tree_base() : max_depth(0), max_leaf_count(0), num_leaf_nodes(0) {}
	bool is_empty() const {return nodes.empty();}
	void clear() {nodes.resize(0);}
	bool get_root_bcube(cube_t &bc) const;
};


template<typename T> class cobj_tree_simple_type_t : public cobj_tree_base {

protected:
	vector<T> objects, temp_bins[3];

	virtual void calc_node_bbox(tree_node &n) const = 0;
	void build_tree(unsigned nix, unsigned skip_dims, unsigned depth);

public:
	virtual ~cobj_tree_simple_type_t() {}
	
	void clear() {
		cobj_tree_base::clear();
		objects.clear(); // reserve(0)?
	}
	void build_tree_top(bool verbose);
};


class cobj_tree_tquads_t : public cobj_tree_simple_type_t<coll_tquad> {

	virtual void calc_node_bbox(tree_node &n) const;

public:
	vector<coll_tquad> &get_tquads_ref() {return objects;}
	void add_cobjs(coll_obj_group const &cobjs, bool verbose);
	void add_polygons(vector<polygon_t> const &polygons, bool verbose);
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA *color, int *cindex, int ignore_cobj, bool exact) const;

	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj, bool exact) const {
		cindex = -1;
		return check_coll_line(p1, p2, cpos, cnorm, NULL, &cindex, ignore_cobj, exact);
	}
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) const {
		return check_coll_line(p1, p2, cpos, cnorm, &color, NULL, 0, exact);
	}
};


struct sphere_with_id_t : public sphere_t {
	unsigned id;
	sphere_with_id_t() : id(0) {}
	sphere_with_id_t(point const &p, float r, unsigned id_) : sphere_t(p, r), id(id_) {}
};


class cobj_tree_sphere_t : public cobj_tree_simple_type_t<sphere_with_id_t> {

	virtual void calc_node_bbox(tree_node &n) const;

public:
	vector<unsigned> ids; // for using in get_ids_int_sphere()
	void add_spheres(vector<sphere_with_id_t> &spheres_, bool verbose);
	void get_ids_int_sphere(point const &center, float radius, vector<unsigned> &ids) const;
};


class cobj_bvh_tree : public cobj_tree_base {

	coll_obj_group const *cobjs;
	vector<unsigned> cixs;
	bool is_static, is_dynamic, occluders_only, cubes_only, inc_voxel_cobjs;

	struct per_thread_data {
		vector<unsigned> temp_bins[3];
		unsigned start_nix, end_nix, cur_nix;
		bool can_be_resized;

		per_thread_data(unsigned start, unsigned end, bool cbr) : start_nix(start), end_nix(end), cur_nix(start), can_be_resized(cbr) {}
		bool at_node_end() const {return (cur_nix == end_nix);}
		void advance_end_range(unsigned new_end_nix) {assert(new_end_nix >= end_nix); end_nix = new_end_nix;}
		unsigned get_next_node_ix() const {assert(cur_nix < end_nix); return cur_nix;}
		void increment_node_ix() {assert(cur_nix >= start_nix); cur_nix++;}
	};

	void add_cobj(unsigned ix) {if (obj_ok((*cobjs)[ix])) cixs.push_back(ix);}
	coll_obj const &get_cobj(unsigned ix) const {return (*cobjs)[cixs[ix]];}
	bool create_cixs();
	void calc_node_bbox(tree_node &n) const;
	void build_tree_top_level_omp();
	void build_tree(unsigned nix, unsigned skip_dims, unsigned depth, per_thread_data &ptd);

	bool obj_ok(coll_obj const &c) const {
		return (((is_static && c.status == COLL_STATIC) || (is_dynamic && c.status == COLL_DYNAMIC) || (!is_static && !is_dynamic)) &&
			(!occluders_only || c.is_occluder()) && !c.cp.no_coll && (!cubes_only || c.type == COLL_CUBE) &&
			(inc_voxel_cobjs || c.cp.cobj_type != COBJ_TYPE_VOX_TERRAIN));
	}

public:
	cobj_bvh_tree(coll_obj_group const *cobjs_, bool s, bool d, bool o, bool c, bool v)
		: cobjs(cobjs_), is_static(s), is_dynamic(d), occluders_only(o), cubes_only(c), inc_voxel_cobjs(v) {assert(cobjs);}

	unsigned get_num_objs() const {return cixs.size();}
	void clear();
	void add_cobj_ids(vector<unsigned> const &cids) {assert(cixs.empty() && !cids.empty()); cixs = cids;}
	void add_cobjs(bool verbose);
	void build_tree_from_cixs(bool do_mt_build);
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, int &cindex, int ignore_cobj,
		bool exact, int test_alpha, bool skip_non_drawn, bool skip_init_colls) const;
	bool check_point_contained(point const &p, int &cindex) const;
	void get_intersecting_cobjs(cube_t const &cube, vector<unsigned> &cobjs, int ignore_cobj, float toler, bool check_ccounter, int id_for_cobj_int) const;
	bool is_cobj_contained(point const &p1, point const &p2, point const &viewer, point const *const pts, unsigned npts, int ignore_cobj, int &cobj) const;
	void get_coll_line_cobjs(point const &pos1, point const &pos2, int ignore_cobj, vector<int> *cobjs, cobj_query_callback *cqc, bool occlude) const;
	void get_coll_sphere_cobjs(point const &center, float radius, int ignore_cobj, vert_coll_detector &vcd) const;
};


#endif // _COBJ_BSP_TREE_H_


