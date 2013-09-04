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
	int vbo_mgr_ix; // high detail
	float height, width, r_angle, rx, ry, rv[3];
	point pos;
	colorRGBA color;

public:
	small_tree() : type(-1) {clear_vbo_mgr_ix();}
	small_tree(point const &p, float h, float w, int t, bool calc_z, rand_gen_t &rgen);
	void setup_rotation(rand_gen_t &rgen);
	vector3d get_rot_dir() const;
	cylinder_3dw get_trunk_cylin() const;
	void add_cobjs(cobj_params &cp, cobj_params &cp_trunk);
	void remove_cobjs();
	bool check_sphere_coll(point &center, float radius) const;
	void clear_vbo_mgr_ix() {vbo_mgr_ix = -1;}
	void calc_points(vbo_vnc_block_manager_t &vbo_manager, bool low_detail, bool update_mode=0);
	void update_points_vbo(vbo_vnc_block_manager_t &vbo_manager, bool low_detail);
	void add_trunk_as_line(vector<point> &points) const;
	colorRGBA get_bark_color() const;
	void draw_pine(vbo_vnc_block_manager_t const &vbo_manager) const;
	void draw(int mode, bool shadow_only, vbo_vnc_block_manager_t const &vbo_manager, vector3d const &xlate=zero_vector, vector<point> *points=NULL) const;
	void translate_by(vector3d const &vd) {pos += vd;}
	bool operator<(small_tree const &t) const {return (type < t.type);} // sort by type
	point get_pos()     const {return pos;}
	float get_height()  const {return height;}
	int get_type ()     const {return type;}
	bool is_pine_tree() const {return (type == T_PINE || type == T_SH_PINE);}
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

	vbo_vnc_block_manager_t vbo_manager[2]; // {high, low} detail
	vector<point> trunk_pts;
	rand_gen_t rgen;
	bool generated;
	unsigned num_pine_trees;
	float max_pt_radius;
	
	small_tree_group() : generated(0), num_pine_trees(0), max_pt_radius(0.0) {}
	void sort_by_type() {stable_sort(begin(), end());}

	void sort_by_dist_to_camera() {
		sort(begin(), end(), small_tree::comp_by_type_dist(get_camera_pos()));
	}
	void add_tree(small_tree &st);
	void calc_trunk_pts();
	void finalize(bool low_detail);
	void finalize_upload_and_clear_pts(bool low_detail);
	void add_trunk_pts(point const &xlate, vector<point> &pts) const;
	void clear_vbos();
	void clear_vbo_manager(int which=3);
	void clear_vbo_manager_and_ids(int which=3);
	void clear_vbo_and_ids_if_needed(bool low_detail);
	void clear_all();
	void add_cobjs_range(iterator b, iterator e);
	void add_cobjs() {add_cobjs_range(begin(), end());}
	void remove_cobjs();
	bool check_sphere_coll(point &center, float radius) const;
	void translate_by(vector3d const &vd);
	void draw_branches(bool shadow_only, vector3d const &xlate=zero_vector, vector<point> *points=NULL) const;
	void draw_pine_leaves(bool shadow_only, bool low_detail=0, bool draw_all_pine=0, vector3d const &xlate=zero_vector) const;
	void draw_non_pine_leaves(bool shadow_only, vector3d const &xlate=zero_vector) const;
	void gen_trees(int x1, int y1, int x2, int y2, float vegetation_);
	unsigned get_gpu_mem() const {return (vbo_manager[0].get_gpu_mem() + vbo_manager[1].get_gpu_mem());}
	bool is_uploaded(bool low_detail) const {return vbo_manager[low_detail].is_uploaded();}
	void update_zmax(float &tzmax) const;
	bool update_zvals(int x1, int y1, int x2, int y2);
	float get_rmax() const {return max_pt_radius;}
};


float calc_tree_size();


#endif // _SMALL_TREE_H_

