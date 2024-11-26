// 3D World - Small Tree (Low-Detail Tree) Code
// by Frank Gennari
// 6/16/12
#pragma once

#include "tree_3dw.h"
#include "draw_utils.h"

enum {TREE_NONE = -1, T_PINE, T_DECID, T_TDECID, T_BUSH, T_PALM, T_SH_PINE, NUM_ST_TYPES};


class small_tree { // size = 116

	short type=-1; // -1 = unset, 0 = pine, 1 = decidious, 2 = tall, 3 = bush, 4 = palm, 5 = short pine
	short inst_id=-1; // for instancing
	int vbo_mgr_ix=-1; // high detail
	float height=0.0, width=0.0, r_angle=0.0, rx=0.0, ry=0.0;
	float branch_xy_scale=1.0; // should only be se to a value != 1.0 for pine trees
	point pos;
	colorRGB leaf_color, bark_color;
	cylinder_3dw trunk_cylin;
	int coll_id[2];

	struct palm_verts_t {
		vector<vert_norm_comp_color> v;
		vector<int> coll_id;
		mutable vbo_wrap_t vbo; // created dynamically during drawing
		unsigned get_gpu_mem() const {return v.size()*sizeof(vert_norm_comp_color);}
	};
	std::shared_ptr<palm_verts_t> palm_verts; // for palm trees only

public:
	small_tree() {init();}
	small_tree(point const &p, unsigned instance_id);
	small_tree(point const &p, float h, float w, int t, bool calc_z, rand_gen_t &rgen, bool allow_rotation=0, float bxys=1.0);
	void init() {coll_id[0] = coll_id[1] = -1; clear_vbo_mgr_ix();}
	void setup_rotation(rand_gen_t &rgen);
	vector3d get_rot_dir() const;
	cylinder_3dw get_trunk_cylin() const;
	unsigned get_palm_mem() const {return (palm_verts ? palm_verts->get_gpu_mem() : 0);}
	void add_cobjs(cobj_params &cp, cobj_params &cp_trunk);
	void remove_cobjs();
	void add_bounds_to_bcube(cube_t &bcube) const;
	void clear_vbo() {if (palm_verts != nullptr) {palm_verts->vbo.clear_vbo();}}
	bool check_sphere_coll(point &center, float radius) const;
	bool line_intersect(point const &p1, point const &p2, float *t=NULL) const;
	void clear_vbo_mgr_ix() {vbo_mgr_ix = -1;}
	void calc_points(vbo_vnc_block_manager_t &vbo_manager, bool low_detail, bool update_mode=0);
	void calc_palm_tree_points();
	void update_points_vbo(vbo_vnc_block_manager_t &vbo_manager, bool low_detail);
	void add_trunk_as_line(vector<point> &points) const;
	vector<vert_norm_comp_color> const &get_palm_verts() const;
	colorRGBA get_leaf_color() const {return leaf_color;}
	void draw_pine(vbo_vnc_block_manager_t const &vbo_manager, unsigned num_instances=1) const;
	bool are_leaves_visible(vector3d const &xlate) const;
	void draw_pine_leaves(vbo_vnc_block_manager_t const &vbo_manager, vector3d const &xlate) const;
	void get_palm_trunk_verts(vector<vert_norm_comp_tc_color> &verts, unsigned nsides) const;
	bool draw_trunk(bool shadow_only, bool all_visible, bool skip_lines=0, vector3d const &xlate=zero_vector, vector<vert_norm_tc> *cylin_verts=nullptr) const;
	void draw_palm_leaves(unsigned num_instances=1) const;
	void draw_leaves(bool shadow_only, int xlate_loc, int scale_loc, vector3d const &xlate=zero_vector) const;
	void translate_by(vector3d const &vd) {pos += vd;}
	bool operator<(small_tree const &t) const {return (type < t.type);} // sort by type
	point get_pos()     const {return pos;}
	point get_center()  const {return pos;}
	float get_height()  const {return height;}
	float get_width()   const {return width;}
	int get_type ()     const {return type;}
	bool is_pine_tree() const {return (type == T_PINE || type == T_SH_PINE);}
	unsigned get_inst_id() const {assert(inst_id >= 0); return inst_id;}
	float get_pine_tree_radius() const;
	float get_radius() const {return (is_pine_tree() ? branch_xy_scale*get_pine_tree_radius() : width);} // approximate
	float get_zmax() const;
	float get_xy_radius() const {return max(1.5*branch_xy_scale*width, 0.5*height);}
	float get_trunk_bsphere_radius() const {return (trunk_cylin.r1 + 0.5*((r_angle == 0.0) ? fabs(trunk_cylin.p1.z - trunk_cylin.p2.z) : trunk_cylin.get_length()));}
	void write_to_cobj_file(std::ostream &out) const;

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
	vbo_wrap_t trunk_pts_vbo;
	vector<point> inst_pts;
	rand_gen_t rgen;
	bool generated=0, instanced=0;
	unsigned num_pine_trees=0, num_palm_trees=0, num_trunk_pts=0, palm_vbo_mem=0;
	float max_tree_radius=0.0;
	point last_cpos;
	cube_t all_bcube;

	struct tree_inst_t {
		unsigned id;
		point pt;

		tree_inst_t(unsigned id_, point const &pt_) : id(id_), pt(pt_) {}
		bool operator<(tree_inst_t const &i) const {return (id < i.id);}
	};
	vector<tree_inst_t> tree_insts[2]; // pine trees, palm trees

	void enable_instanced() {instanced |= ((num_pine_trees + num_palm_trees) == size());} // only if all are pine/palm trees
	void sort_by_type() {stable_sort(begin(), end());}
	void sort_by_dist_to_camera();
	void add_tree(small_tree const &st);
	void finalize(bool low_detail);
	void finalize_upload_and_clear_pts(bool low_detail);
	void draw_trunk_pts();
	void draw(bool shadow_only, int reflection_pass);
	void clear_vbos();
	void clear_vbo_manager(int which=3);
	void clear_vbo_manager_and_ids(int which=3);
	void clear_vbo_and_ids_if_needed(bool low_detail);
	void clear_all();
	void add_cobjs_range(iterator b, iterator e);
	void add_cobjs() {add_cobjs_range(begin(), end());}
	void remove_cobjs();
	void calc_bcube();
	bool check_sphere_coll(point &center, float radius) const;
	bool line_intersect(point const &p1, point const &p2, float *t=NULL) const;
	void translate_by(vector3d const &vd);
	void get_back_to_front_ordering(vector<pair<float, unsigned> > &to_draw, vector3d const &xlate) const;
	bool draw_trunks(bool shadow_only, bool all_visible=0, bool skip_lines=0, vector3d const &xlate=zero_vector) const;
	void draw_tree_insts(shader_t &s, bool draw_all, vector3d const &xlate, int xlate_loc, vector<tree_inst_t> &insts, bool is_pine);
	void draw_pine_insts(shader_t &s, bool draw_all, vector3d const &xlate, int xlate_loc) {draw_tree_insts(s, draw_all, xlate, xlate_loc, tree_insts[0], 1);}
	void draw_palm_insts(shader_t &s, bool draw_all, vector3d const &xlate, int xlate_loc) {draw_tree_insts(s, draw_all, xlate, xlate_loc, tree_insts[1], 0);}
	void draw_pine_leaves(shader_t &s, bool shadow_only, bool low_detail=0, bool draw_all=0, bool sort_front_to_back=0, vector3d const &xlate=zero_vector, int xlate_loc=-1);
	void draw_non_pine_leaves(bool shadow_only, bool draw_palm, bool draw_non_palm, int xlate_loc=-1, int scale_loc=-1, vector3d const &xlate=zero_vector) const;
	int get_ntrees_for_mesh_xy(int i, int j, float ntrees_mult_density);
	void maybe_add_tree(int i, int j, float zpos_in, float tsize, int skip_val, bool check_hmap_normal);
	void gen_trees(int x1, int y1, int x2, int y2, float const density[4]);
	void gen_trees_tt_within_radius(int x1, int y1, int x2, int y2, point const &pos, float radius, bool is_square, float mesh_dz=-1.0, tile_t const *const cur_tile=nullptr);
	unsigned get_gpu_mem() const {return (palm_vbo_mem + vbo_manager[0].get_gpu_mem() + vbo_manager[1].get_gpu_mem());}
	bool is_uploaded(bool low_detail) const {return vbo_manager[low_detail].is_uploaded();}
	void update_zmax(float &tzmax) const;
	bool update_zvals(int x1, int y1, int x2, int y2);
	float get_rmax() const {return max_tree_radius;}
	small_tree const &get_tree(unsigned ix) const {assert(ix < size()); return operator[](ix);}
};


float calc_tree_size();
bool can_have_pine_palm_trees_in_zrange(float z_min, float z_max, bool skip_range_check_if_manually_placed=0);
bool can_have_decid_trees_in_zrange(float z_min, float z_max, bool skip_range_check_if_manually_placed=0);

