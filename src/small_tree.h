// 3D World - Small Tree (Low-Detail Tree) Code
// by Frank Gennari
// 6/16/12

#ifndef _SMALL_TREE_H_
#define _SMALL_TREE_H_

#include "tree_3dw.h"
#include "draw_utils.h"


enum {TREE_NONE = -1, T_PINE, T_DECID, T_TDECID, T_BUSH, T_PALM, T_SH_PINE, NUM_ST_TYPES};


class small_tree { // size = 81 (82)

	char type; // 0 = pine, 1 = decidious, 2 = tall, 3 = bush, 4 = palm, 5 = short pine
	vector<int> coll_id;
	int vbo_mgr_ix[2]; // {high, low} detail
	float height, width, r_angle, rx, ry, rv[3];
	point pos;
	colorRGBA color;

public:
	small_tree() : type(-1) {clear_vbo_mgr_ix();}
	small_tree(point const &p, float h, float w, int t, bool calc_z);
	void setup_rotation();
	vector3d get_rot_dir() const;
	void add_cobjs(cobj_params &cp, cobj_params &cp_trunk);
	void remove_cobjs();
	void clear_vbo_mgr_ix(int which=3);
	void calc_points(vbo_quad_block_manager_t &vbo_manager, bool low_detail);
	void add_trunk_as_line(vector<point> &points) const;
	colorRGBA get_bark_color() const;
	void draw(int mode, bool shadow_only, bool do_cull, vbo_quad_block_manager_t const &vbo_manager,
		bool low_detail=0, vector3d const xlate=zero_vector, vector<point> *points=NULL) const;
	void translate_by(vector3d const &vd) {pos += vd;}
	bool operator<(small_tree const &t) const {return (type < t.type);} // sort by type
	point get_pos() const {return pos;}
	int get_type () const {return type;}
	float get_pine_tree_radius() const;
	float get_zmax() const;

	struct comp_by_type_dist {
		point pos;
		comp_by_type_dist(point const &pos_) : pos(pos_) {}

		bool operator()(small_tree const &a, small_tree const &b) const {
			return ((a.type != b.type) ? (a.type < b.type) : (p2p_dist_sq(a.pos, pos) < p2p_dist_sq(b.pos, pos)));
		}
	};
};


struct small_tree_group : public vector<small_tree> {

	vbo_quad_block_manager_t vbo_manager[2]; // {high, low} detail
	vector<point> trunk_pts;
	bool generated;
	unsigned num_pine_trees;
	
	small_tree_group() : generated(0), num_pine_trees(0) {}
	void sort_by_type() {sort(begin(), end());}

	void sort_by_dist_to_camera() {
		sort(begin(), end(), small_tree::comp_by_type_dist(get_camera_pos()));
	}
	void add_tree(small_tree &st);
	void finalize(bool low_detail);
	void add_trunk_pts(point const &xlate, vector<point> &pts) const;
	void clear_vbos();
	void clear_vbo_manager(int which=3);
	void clear_vbo_manager_and_ids(int which=3);
	void clear_all();
	void add_cobjs();
	void remove_cobjs();
	void translate_by(vector3d const &vd);
	void draw_branches(bool shadow_only, vector3d const xlate=zero_vector, vector<point> *points=NULL) const;
	void draw_leaves(bool shadow_only, bool low_detail=0, bool draw_all_pine=0, vector3d const xlate=zero_vector) const;
	void gen_trees(int x1, int y1, int x2, int y2);
};


float calc_tree_size();


#endif // _SMALL_TREE_H_

