// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 5/12/02

#ifndef _TREE_H_
#define _TREE_H_

#include "3DWorld.h"
#include "gl_ext_arb.h" // for indexed_vbo_manager_t


float const TREE_DIST_SCALE = 100.0;
float const TREE_DEN_THRESH = 0.55;


struct blastr; // forward reference

// small tree classes
enum {TREE_CLASS_NONE=0, TREE_CLASS_PINE, TREE_CLASS_DECID, TREE_CLASS_PALM, TREE_CLASS_DETAILED, NUM_TREE_CLASSES};


class tree_lod_render_t {

	struct leaf_inst_data_t {
		point pos;
		float radius;
		leaf_inst_data_t(point const &pos_=all_zeros, float radius_=0.0) : pos(pos_), radius(radius_) {}
	};

	typedef map<texture_pair_t, vector<leaf_inst_data_t> > leaf_map_t;
	leaf_map_t leaf_map;
	bool enabled;

public:
	tree_lod_render_t(bool enabled_) : enabled(enabled_) {}
	bool is_enabled() const {return enabled;}
	bool empty() {return leaf_map.empty();}
	void clear() {leaf_map.clear();}
	void add_leaves(texture_pair_t const &tp, point const &pos, float radius) ;
	void render_quads_facing_camera() const;
};


struct tree_leaf { // size = 28 + 48 = 76

	int shadow_bits;
	float color, lred, lgreen;
	vector3d norm;
	point pts[4];

	tree_leaf() : shadow_bits(0) {}
	void create_init_color(bool deterministic);
	colorRGB calc_leaf_color(colorRGBA const &leaf_color, colorRGBA const &base_color) const;
	float get_norm_scale(unsigned pt_ix) const;
	point get_center() const {return 0.25*(pts[0] + pts[1] + pts[2] + pts[3]);} // average of all 4 leaf points
};


inline bool comp_leaf(const tree_leaf &A, const tree_leaf &B) {
	return (A.pts[0].mag_sq() < B.pts[0].mag_sq());
}


struct draw_cylin : public cylinder_3dw { // size = 35 (36)

	unsigned char level;
	unsigned short branch_id;

	draw_cylin() : level(0), branch_id(0) {}
	unsigned get_num_div() const {return (N_CYL_SIDES >> 1) - ((level - 1) << 2);}
	bool can_merge(draw_cylin const &c) const {return (level == c.level && branch_id == c.branch_id);}
};


struct tree_cylin : public draw_cylin { // size = 55 (56)

	float length, deg_rotate;
	vector3d rotate;

	void assign_params(unsigned char lev, unsigned short bid, float r1_, float r2_, float len, float drot) {
		level = lev; branch_id = bid; r1 = r1_; r2 = r2_; length = len; deg_rotate = drot;
	}
};


struct tree_branch { // size = 12

	tree_cylin *cylin;
	float total_length;
	short num_cylins, num_branches;

	void clear_num() {num_cylins = num_branches = 0;}
};


class tree_builder_t {

	static vector<tree_cylin >   cylin_cache;
	static vector<tree_branch>   branch_cache;
	static vector<tree_branch *> branch_ptr_cache;

	tree_branch base, roots, *branches_34[2], **branches;
	int base_num_cylins, root_num_cylins, ncib, num_1_branches, num_big_branches_min, num_big_branches_max;
	int num_2_branches_min, num_2_branches_max, num_34_branches[2], num_3_branches_min, num_3_branches_max;
	int tree_slimness, tree_wideness, base_break_off;
	float base_radius, base_length_min, base_length_max, base_curveness, num_leaves_per_occ;
	float branch_curveness, branch_upwardness, branch_distribution, branch_1_distribution, num_cylin_factor, base_cylin_factor;
	float branch_1_var, branch_1_rad_var, branch_1_start, branch_2_var, branch_2_rad_var, branch_2_start, branch_4_max_radius, rotate_factor;
	float angle_rotate, branch_min_angle, branch_max_angle, branch_4_length;
	float max_2_angle_rotate, max_3_angle_rotate;  //max angle to rotate 3rd order branches around from the 2nd order branch

	float gen_bc_size(float branch_var);
	float gen_bc_size2(float branch_var);
	void gen_next_cylin(tree_cylin &cylin, tree_cylin &lcylin, float var, float rad_var, int level, int branch_id, bool rad_var_test);
	void gen_first_cylin(tree_cylin &cylin, tree_cylin &src_cylin, float bstart, float rad_var, float rotate_start, int level, int branch_id);
	void create_1_order_branch(int base_cylin_num, float rotate_start, int branch_num);
	void create_2nd_order_branch(int i, int j, int cylin_num, bool branch_deflected, int rotation);
	void create_3rd_order_branch(int i, int j, int cylin_num, int branch_num, bool branch_deflected, int rotation);
	void gen_b4(tree_branch &branch, int &branch_num, int num_4_branches, int i, int k);
	void create_4th_order_branches(int nbranches);
	void generate_4th_order_branch(tree_branch &src_branch, int j, float rotate_start, float temp_deg, int branch_num);
	int generate_next_cylin(int cylin_num, int ncib, bool branch_just_created, bool &branch_deflected);
	void add_leaves_to_cylin(tree_cylin const &cylin, float tsize, float deadness, vector<tree_leaf> &leaves) const;
	void process_cylins(tree_cylin *const cylins, unsigned num, int tree_type, float deadness,
		vector<draw_cylin> &all_cylins, vector<tree_leaf> &leaves);

public:
	float create_tree_branches(int tree_type, int size, float tree_depth, colorRGBA &base_color);
	void create_all_cylins_and_leaves(int tree_type, float deadness, vector<draw_cylin> &all_cylins, vector<tree_leaf> &leaves);
};


class tree_data_t {

	typedef vert_norm_comp_color leaf_vert_type_t;
	typedef vert_norm_comp_tc branch_vert_type_t;
	typedef unsigned short branch_index_t;

	indexed_vbo_manager_t branch_vbo_manager;
	unsigned leaf_vbo, num_branch_quads, num_unique_pts;
	int tree_type;
	colorRGBA base_color, leaf_color;
	vector<leaf_vert_type_t> leaf_data;
	vector<draw_cylin> all_cylins;
	vector<tree_leaf> leaves;
	texture_pair_t render_leaf_texture;
	int last_update_frame;
	//unsigned ref_count;
	bool leaves_changed, reset_leaves;

	void create_indexed_quads_for_branches(vector<branch_vert_type_t> &data, vector<branch_index_t> &idata) const;
	void clear_vbo_ixs();

public:
	float base_radius, sphere_radius, sphere_center_zoff;

	tree_data_t(bool priv=1) : leaf_vbo(0), num_branch_quads(0), num_unique_pts(0), tree_type(-1), last_update_frame(0),
		leaves_changed(0), reset_leaves(0), base_radius(0.0), sphere_radius(0.0), sphere_center_zoff(0.0) {}
	vector<draw_cylin> const &get_all_cylins() const {return all_cylins;}
	vector<tree_leaf>  const &get_leaves    () const {return leaves;}
	vector<tree_leaf>        &get_leaves    ()       {return leaves;}
	void make_private_copy(tree_data_t &dest) const;
	void gen_tree_data(int tree_type_, int size, float tree_depth);
	void gen_leaf_color();
	void update_all_leaf_colors();
	void update_leaf_color(unsigned i, bool no_mark_changed=0);
	colorRGB get_leaf_color(unsigned i) const;
	bool leaf_data_allocated() const {return !leaf_data.empty();}
	bool is_created() const {return !all_cylins.empty();} // as good a check as any
	bool check_if_needs_updated();
	void remove_leaf_ix(unsigned i, bool update_data);
	void bend_leaf(unsigned i, float angle);
	void draw_tree_shadow_only(bool draw_branches, bool draw_leaves) const;
	void draw_branches(float size_scale);
	void draw_leaves(float size_scale);
	texture_pair_t get_render_leaf_texture() const {return render_leaf_texture;}
	bool leaf_draw_setup(bool leaf_dynamic_en);
	void check_leaf_render_texture();
	void update_normal_for_leaf(unsigned i);
	void reset_leaf_pos_norm();
	void alloc_leaf_data() {leaf_data.resize(4*leaves.size());}
	void clear_data();
	void clear_context();
	unsigned get_gpu_mem() const;
	int get_tree_type() const {return tree_type;}
	point get_center() const {return point(0.0, 0.0, sphere_center_zoff);}

	static void pre_draw (bool branches_or_leaves, bool shadow_only);
	static void post_draw(bool branches_or_leaves, bool shadow_only);
};


class tree
{
	tree_data_t priv_tree_data; // by pointer?
	tree_data_t *tree_data; // by index?
	void make_private_tdata_copy();
	tree_data_t const &tdata() const {return (tree_data ? *tree_data : priv_tree_data);}
	tree_data_t       &tdata()       {return (tree_data ? *tree_data : priv_tree_data);}
	bool td_is_private() const {return (tree_data == NULL);}

	int type, created; // should type be a member of tree_data_t?
	bool no_delete, not_visible, leaf_orients_valid;
	point tree_center;
	float damage, damage_scale;
	colorRGBA tree_color, bcolor;
	vector<int> branch_cobjs, leaf_cobjs;

	coll_obj &get_leaf_cobj(unsigned i) const;
	void update_leaf_orients();
	bool has_leaf_data()   const {return tdata().leaf_data_allocated();}
	bool physics_enabled() const;
	void get_abs_leaf_pts(point pts[4], unsigned ix) const;
	void create_leaf_obj(unsigned ix) const;

	bool is_over_mesh() const;
	bool is_visible_to_camera(vector3d const &xlate) const;
	void burn_leaves();
	void blast_damage(blastr const *const blast_radius);
	void lightning_damage(point const &ltpos);
	void drop_leaves();
	void remove_leaf(unsigned i, bool update_data);
	bool damage_leaf(unsigned i, float damage_done);
	void draw_tree_branches(shader_t const &s, float size_scale, vector3d const &xlate, int shader_loc);
	void draw_tree_leaves(shader_t const &s, float size_scale, vector3d const &xlate);
	void update_leaf_cobj_color(unsigned i);
	void copy_color(unsigned i, bool no_mark_changed=0);

public:
	tree() : tree_data(NULL), created(0), no_delete(0), not_visible(0), leaf_orients_valid(0) {}
	void bind_to_td(tree_data_t *td);
	void gen_tree(point const &pos, int size, int ttype, int calc_z, bool add_cobjs, bool user_placed);
	void calc_leaf_shadows();
	void gen_tree_shadows(unsigned light_sources);
	void add_tree_collision_objects();
	void remove_collision_objects();
	bool check_sphere_coll(point &center, float radius) const;
	void draw_tree(shader_t const &s, tree_lod_render_t &lod_renderer, bool draw_branches, bool draw_leaves,
		bool shadow_only, vector3d const &xlate, int shader_loc);
	void shift_tree(vector3d const &vd) {tree_center += vd;}
	void clear_context();
	int delete_tree();
	int get_type()            const {return type;}
	float get_radius()        const {return tdata().sphere_radius;}
	point sphere_center()     const {return (tree_center + tdata().get_center());}
	point const &get_center() const {return tree_center;}
	unsigned get_gpu_mem()    const {return (td_is_private() ? tdata().get_gpu_mem() : 0);}
	bool get_no_delete()      const {return no_delete;}
	void set_no_delete(bool no_delete_) {no_delete = no_delete_;}
	bool operator<(tree const &t) const {return ((type != t.type) ? (type < t.type) : (tree_data < t.tree_data));}
	void check_leaf_render_texture() {tdata().check_leaf_render_texture();}
};


class tree_data_manager_t : public vector<tree_data_t> {

	float last_tree_scale;
	int last_rgi;

public:
	tree_data_manager_t() : last_tree_scale(1.0), last_rgi(0) {}
	void ensure_init();
	void clear_context();
	unsigned get_gpu_mem() const;
};


class tree_cont_t : public vector<tree> {

	tree_data_manager_t &shared_tree_data;
	bool generated;

public:
	tree_cont_t(tree_data_manager_t &tds) : shared_tree_data(tds), generated(0) {}
	bool was_generated() const {return generated;}
	void remove_cobjs();
	bool check_sphere_coll(point &center, float radius) const;
	void draw_branches_and_leaves(shader_t const &s, tree_lod_render_t &lod_renderer, bool draw_branches, bool draw_leaves, bool shadow_only, vector3d const &xlate);
	void check_leaf_shadow_change();
	static void pre_leaf_draw(shader_t &shader);
	static void post_leaf_draw(shader_t &shader);
	void draw(bool shadow_only);
	unsigned delete_all();
	unsigned scroll_trees(int ext_x1, int ext_x2, int ext_y1, int ext_y2);
	void post_scroll_remove();
	void gen_deterministic(int x1, int y1, int x2, int y2, float vegetation_);
	void shift_by(vector3d const &vd);
	void add_cobjs();
	void clear_context();
	void clear() {delete_all(); vector<tree>::clear();}
	unsigned get_gpu_mem() const;
	float get_rmax() const;
	void update_zmax(float &tzmax) const;
	void check_leaf_render_textures();
};


// function prototypes - tree
float get_tree_z_bottom(float z, point const &pos);
void remove_tree_cobjs();
void draw_trees(bool shadow_only=0);
void delete_trees();
void regen_trees(bool recalc_shadows, bool keep_old);
void shift_trees(vector3d const &vd);
void add_tree_cobjs();
void clear_tree_context();

// function prototypes - small trees
int add_small_tree(point const &pos, float height, float width, int tree_type, bool calc_z);
void add_small_tree_coll_objs();
void remove_small_tree_cobjs();
void gen_small_trees();
void draw_small_trees(bool shadow_only);
void shift_small_trees(vector3d const &vd);


#endif // _TREE_H_

