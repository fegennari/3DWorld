// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 5/12/02

#ifndef _TREE_H_
#define _TREE_H_

#include "3DWorld.h"
#include "gl_ext_arb.h" // for indexed_vbo_manager_t


float const TREE_DIST_SCALE = 100.0;
float const TREE_DEPTH      = 0.1;


struct blastr; // forward reference
class tree_data_t;
class cobj_bvh_tree;

// small tree classes
enum {TREE_CLASS_NONE=0, TREE_CLASS_PINE, TREE_CLASS_DECID, TREE_CLASS_PALM, TREE_CLASS_DETAILED, NUM_TREE_CLASSES};


class tree_lod_render_t {

	struct entry_t {
		tree_data_t const *td;
		point pos;
		color_wrapper cw;

		entry_t() : td(NULL) {}
		entry_t(tree_data_t const *td_, point const &pos_, colorRGBA const &color) : td(td_), pos(pos_) {assert(td); cw.set_c4(color);}
		bool operator<(entry_t const &e) const {return (td < e.td);} // compare tree data pointer values
	};

	vector<entry_t> leaf_vect, branch_vect;
	bool enabled;

public:
	int leaf_opacity_loc, branch_opacity_loc;

	tree_lod_render_t(bool enabled_) : enabled(enabled_), leaf_opacity_loc(-1), branch_opacity_loc(-1) {}
	void set_enabled(bool enabled_) {enabled = enabled_;} // to be called before use, not during rendering
	bool is_enabled()   const {return enabled;}
	bool has_leaves()   const {return !leaf_vect.empty();}
	bool has_branches() const {return !branch_vect.empty();}
	bool empty()        const {return (!has_leaves() && !has_branches());}
	void clear()       {leaf_vect.clear(); branch_vect.clear();}
	void resize_zero() {leaf_vect.resize(0); branch_vect.resize(0);}

	void add_leaves(tree_data_t const *td, point const &pos, float opacity) {
		leaf_vect.push_back(entry_t(td, pos, colorRGBA(1, 1, 1, opacity)));
	}
	void add_branches(tree_data_t const *td, point const &pos, float opacity, colorRGBA const &bcolor) {
		branch_vect.push_back(entry_t(td, pos, colorRGBA(bcolor, opacity)));
	}
	void finalize();
	void render_billboards(bool render_branches) const;
};


struct tree_leaf { // size = 28 + 48 = 76

	int shadow_bits;
	float color, lred, lgreen;
	vector3d norm;
	point pts[4];

	tree_leaf() : shadow_bits(0), color(0.0), lred(0.0), lgreen(0.0) {}
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
	unsigned get_num_div() const {return (N_CYL_SIDES >> 1) - ((level - 1) << 2);} // 0:20 1:16 2:12 3:8 4:4
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


struct tree_xform_t {

	float tree_temp_matrix[3], i_matrix[9], re_matrix[3];
	float langle, sin_term, cos_term;

	tree_xform_t() : langle(0.0), sin_term(0.0), cos_term(1.0) {}

	void add_rotation(point &dest, point const &src, float mult) {
		UNROLL_3X(dest[i_] = src[i_] + mult*re_matrix[i_];)
	}
	void setup_rotate(vector3d &rotate, float rotate_start, float temp_deg) {
		float const angle(rotate_start/TO_DEG + temp_deg);
		rotate.assign(cosf(angle), sinf(angle), 0.0);
	}
	void rotate_around_axis(tree_cylin const &c) {
		rotate_all(c.rotate, c.deg_rotate/TO_DEG, 0.0, 0.0, c.length);
	}
	void rotate_pts_around_axis(point const &p, point &rotation_v, float deg_rotate) {
		rotate_all(rotation_v, deg_rotate/TO_DEG, p.x, p.y, p.z);
	}
	void rotate_cylin(tree_cylin &c);
	void rotate_it();
	void rotate_itv2();
	void turn_into_i(float *m);
	void rotate_all(point const &rotate, float angle, float x, float y, float z);
	void gen_cylin_rotate(vector3d &rotate, vector3d &lrotate, float rotate_start);
};


class tree_builder_t : public tree_xform_t {

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
	float angle_rotate, branch_min_angle, branch_max_angle, branch_4_length, leaf_acc;
	float max_2_angle_rotate, max_3_angle_rotate;  //max angle to rotate 3rd order branches around from the 2nd order branch
	cube_t const *clip_cube;

	float gen_bc_size(float branch_var);
	float gen_bc_size2(float branch_var);
	void gen_next_cylin(tree_cylin &cylin, tree_cylin &lcylin, float var, float rad_var, int level, int branch_id, bool rad_var_test);
	void gen_first_cylin(tree_cylin &cylin, tree_cylin &src_cylin, float bstart, float rad_var, float rotate_start, int level, int branch_id);
	void create_1_order_branch(int base_cylin_num, float rotate_start, int branch_num);
	void create_2nd_order_branch(int i, int j, int cylin_num, bool branch_deflected, int rotation);
	void create_3rd_order_branch(int i, int j, int cylin_num, int branch_num, bool branch_deflected, int rotation);
	void gen_b4(tree_branch &branch, int &branch_num, int num_4_branches, int i, int k);
	void create_4th_order_branches(int nbranches, float branch_scale);
	void generate_4th_order_branch(tree_branch &src_branch, int j, float rotate_start, float temp_deg, int branch_num);
	int generate_next_cylin(int cylin_num, int ncib, bool branch_just_created, bool &branch_deflected);
	void add_leaves_to_cylin(unsigned cylin_ix, int tree_type, float rel_leaf_size, float deadness, vector<tree_leaf> &leaves);

public:
	tree_builder_t(cube_t const *clip_cube_) : branches(NULL), clip_cube(clip_cube_) {branches_34[0] = branches_34[1] = NULL;}
	float create_tree_branches(int tree_type, int size, float tree_depth, colorRGBA &base_color,
		float height_scale, float br_scale, float nl_scale, float bbo_scale, bool has_4th_branches, bool create_bush);
	void create_all_cylins_and_leaves(vector<draw_cylin> &all_cylins, vector<tree_leaf> &leaves,
		int tree_type, float deadness, float br_scale, float nl_scale, bool has_4th_branches);
};


//#define USE_TREE_BB_TEX_ATLAS

#ifdef USE_TREE_BB_TEX_ATLAS
typedef texture_atlas_t tree_bb_tex_t;
#else
typedef texture_pair_t tree_bb_tex_t;
#endif

bool const TREE_BILLBOARD_MULTISAMPLE = 0;


class tree_data_t {

	typedef vert_norm_comp_color leaf_vert_type_t;
	typedef vert_norm_comp_tc branch_vert_type_t;

	indexed_vbo_manager_t branch_manager;
	unsigned leaf_vbo, num_branch_quads, num_unique_pts, branch_index_bytes;
	int tree_type;
	colorRGBA base_color, leaf_color;
	vector<leaf_vert_type_t> leaf_data;
	vector<draw_cylin> all_cylins;
	vector<tree_leaf> leaves;
	tree_bb_tex_t render_leaf_texture, render_branch_texture;
	int last_update_frame;
	unsigned leaf_change_start, leaf_change_end;
	bool reset_leaves, has_4th_branches;

	void clear_vbo_ixs();
	template<typename branch_index_t> void create_branch_vbo();

public:
	float base_radius, sphere_radius, sphere_center_zoff, br_scale, b_tex_scale;
	float lr_z_cent, lr_x, lr_y, lr_z, br_x, br_y, br_z; // bounding cylinder data for leaves and branches

	tree_data_t(bool priv=1) : leaf_vbo(0), num_branch_quads(0), num_unique_pts(0), branch_index_bytes(0), tree_type(-1), last_update_frame(0),
		render_leaf_texture(TREE_BILLBOARD_MULTISAMPLE), render_branch_texture(TREE_BILLBOARD_MULTISAMPLE),
		leaf_change_start(0), leaf_change_end(0), reset_leaves(0), has_4th_branches(0), base_radius(0.0), sphere_radius(0.0), sphere_center_zoff(0.0),
		br_scale(1.0), b_tex_scale(1.0), lr_z_cent(0.0), lr_x(0.0), lr_y(0.0), lr_z(0.0), br_x(0.0), br_y(0.0), br_z(0.0) {}
	vector<draw_cylin> const &get_all_cylins() const {return all_cylins;}
	vector<tree_leaf>  const &get_leaves    () const {return leaves;}
	vector<tree_leaf>        &get_leaves    ()       {return leaves;}
	void make_private_copy(tree_data_t &dest) const;
	void gen_tree_data(int tree_type_, int size, float tree_depth, float height_scale, float br_scale_mult, float nl_scale,
		float bbo_scale, bool has_4th_branches_, cube_t const *clip_cube, bool create_bush);
	void mark_leaf_changed(unsigned ix);
	void gen_leaf_color();
	void update_all_leaf_colors();
	void update_leaf_color(unsigned i, bool no_mark_changed=0);
	colorRGB get_leaf_color(unsigned i) const;
	bool leaf_data_allocated() const {return !leaf_data.empty();}
	bool is_created() const {return !all_cylins.empty();} // as good a check as any
	bool leaf_vbo_valid() const {return (leaf_vbo > 0);}
	bool get_has_4th_branches() const {return has_4th_branches;}
	float get_size_scale_mult() const;
	bool check_if_needs_updated();
	void remove_leaf_ix(unsigned i, bool update_data);
	bool spraypaint_leaves(point const &pos, float radius, colorRGBA const &color, bool check_only);
	void bend_leaf(unsigned i, float angle);
	void draw_leaf_quads_from_vbo(unsigned max_leaves) const;
	void draw_leaves_shadow_only(float size_scale);
	void ensure_branch_vbo();
	void draw_branches(shader_t &s, float size_scale, bool force_low_detail);
	void ensure_leaf_vbo();
	void draw_leaves(float size_scale);
	tree_bb_tex_t const &get_render_leaf_texture  () const {return render_leaf_texture  ;}
	tree_bb_tex_t const &get_render_branch_texture() const {return render_branch_texture;}
	bool leaf_draw_setup(bool no_leaf_reset);
	void check_render_textures();
	void update_normal_for_leaf(unsigned i);
	void reset_leaf_pos_norm();
	void alloc_leaf_data() {leaf_data.resize(4*leaves.size());}
	void clear_data();
	void clear_context();
	unsigned get_gpu_mem() const;
	int get_tree_type() const {return tree_type;}
	point get_center() const {return point(0.0, 0.0, sphere_center_zoff);}

	static void pre_branch_draw(shader_t &s, bool shadow_only);
	static void pre_leaf_draw(shader_t &s);
	static void post_branch_draw(bool shadow_only);
	static void post_leaf_draw();
};


class tree {

	tree_data_t priv_tree_data; // by pointer?
	tree_data_t *tree_data; // by index?
	void make_private_tdata_copy();
	tree_data_t const &tdata() const {return (tree_data ? *tree_data : priv_tree_data);}
	tree_data_t       &tdata()       {return (tree_data ? *tree_data : priv_tree_data);}
	bool td_is_private() const {return (tree_data == NULL);}

	int type, created; // should type be a member of tree_data_t?
	bool no_delete, not_visible, leaf_orients_valid, enable_leaf_wind, use_clip_cube;
	point tree_center;
	float damage, damage_scale, last_size_scale, tree_nl_scale;
	colorRGBA tree_color;
	vector<int> branch_cobjs, leaf_cobjs;
	cube_t clip_cube;

	coll_obj &get_leaf_cobj(unsigned i) const;
	void update_leaf_orients();
	bool has_leaf_data()   const {return tdata().leaf_data_allocated();}
	bool physics_enabled() const;
	void get_abs_leaf_pts(point pts[4], unsigned ix) const;
	void create_leaf_obj(unsigned ix) const;

	bool is_over_mesh() const;
	bool is_visible_to_camera(vector3d const &xlate) const;
	void burn_leaves();
	void lightning_damage(point const &ltpos);
	void drop_leaves();
	void remove_leaf(unsigned i, bool update_data);
	bool damage_leaf(unsigned i, float damage_done);
	void update_leaf_cobj_color(unsigned i);
	void copy_color(unsigned i, bool no_mark_changed=0);

public:
	tree(bool en_lw=1) : type(-1), tree_data(NULL), created(0), no_delete(0), not_visible(0), leaf_orients_valid(0),
		enable_leaf_wind(en_lw), use_clip_cube(0), damage(0.0), damage_scale(0.0), last_size_scale(0.0), tree_nl_scale(1.0) {}
	void enable_clip_cube(cube_t const &cc) {clip_cube = cc; use_clip_cube = 1;}
	void bind_to_td(tree_data_t *td);
	void gen_tree(point const &pos, int size, int ttype, int calc_z, bool add_cobjs, bool user_placed,
		float height_scale=1.0, float br_scale_mult=1.0, float nl_scale=1.0, bool has_4th_branches=0);
	void calc_leaf_shadows();
	void add_tree_collision_objects();
	void remove_collision_objects();
	bool check_sphere_coll(point &center, float radius) const;
	float calc_size_scale(point const &draw_pos) const;
	void draw_branches_top(shader_t &s, tree_lod_render_t &lod_renderer, bool shadow_only, bool reflection_pass, vector3d const &xlate, int wsoff_loc);
	void draw_leaves_top(shader_t &s, tree_lod_render_t &lod_renderer, bool shadow_only, vector3d const &xlate, int wsoff_loc, int tex0_loc);
	void shift_tree(vector3d const &vd) {tree_center += vd;}
	void clear_context();
	int delete_tree();
	int get_type()            const {return type;}
	float get_radius()        const {return tdata().sphere_radius;}
	float get_height_offset() const {return tdata().sphere_center_zoff;}
	point sphere_center()     const {return (tree_center + tdata().get_center());}
	point const &get_center() const {return tree_center;}
	unsigned get_gpu_mem()    const {return (td_is_private() ? tdata().get_gpu_mem() : 0);}
	bool get_no_delete()      const {return no_delete;}
	void set_no_delete(bool no_delete_) {no_delete = no_delete_;}
	bool operator<(tree const &t) const {return ((type != t.type) ? (type < t.type) : (tree_data < t.tree_data));}
	void check_render_textures() {tdata().check_render_textures();}
	bool spraypaint_leaves(point const &pos, float radius, colorRGBA const &color);
	void blast_damage(blastr const *const blast_radius);
	void write_to_cobj_file(std::ostream &out) const;
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
	void draw_branches_and_leaves(shader_t &s, tree_lod_render_t &lod_renderer, bool draw_branches, bool draw_leaves,
		bool shadow_only, bool reflection_pass, vector3d const &xlate);
	void check_leaf_shadow_change();
	static void pre_leaf_draw(shader_t &shader, bool enable_opacity, bool shadow_only=0, bool use_fs_smap=0);
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
	bool update_zvals(int x1, int y1, int x2, int y2);
	void spraypaint_leaves(point const &pos, float radius, colorRGBA const &color);
	void check_render_textures();
	void apply_exp_damage(point const &epos, float damage, float bradius, int type);
};


// function prototypes - tree
float get_tree_z_bottom(float z, point const &pos);
void remove_tree_cobjs();
void draw_trees(bool shadow_only=0);
void delete_trees();
void regen_trees(bool keep_old);
void shift_trees(vector3d const &vd);
void add_tree_cobjs();
void clear_tree_context();
void spraypaint_tree_leaves(point const &pos, float radius, colorRGBA const &color);

// function prototypes - small trees
int add_small_tree(point const &pos, float height, float width, int tree_type, bool calc_z);
void add_small_tree_coll_objs();
void remove_small_tree_cobjs();
void gen_small_trees();
void draw_small_trees(bool shadow_only);
void shift_small_trees(vector3d const &vd);


#endif // _TREE_H_

