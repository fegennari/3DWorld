// 3D World - 3D Model Rendering Classes
// by Frank Gennari
// 8/17/11

#ifndef _MODEL3D_H_
#define _MODEL3D_H_

#include "3DWorld.h"
#include "collision_detect.h" // for polygon_t
#include "cobj_bsp_tree.h" // for cobj_tree_tquads_t
#include "shadow_map.h" // for smap_data_t

using namespace std;

typedef map<string, unsigned> string_map_t;

unsigned const MAX_VMAP_SIZE     = (1 << 18); // 256K
float const POLY_COPLANAR_THRESH = 0.98;


struct geom_xform_t { // should be packed, can read/write as POD

	vector3d tv;
	float scale;
	bool mirror[3], swap_dim[3][3];

	geom_xform_t(vector3d const &tv_=zero_vector, float scale_=1.0) : tv(tv_), scale(scale_) {
		restore_mirror_and_swap();
	}
	void restore_mirror_and_swap() {
		for (unsigned i = 0; i < 3; ++i) {
			UNROLL_3X(swap_dim[i][i_] = 0;)
			mirror[i] = 0;
		}
	}
	void xform_pos_rm(point &pos) const {
		UNROLL_3X(if (mirror[i_]) {pos[i_] = -pos[i_];})
		
		for (unsigned i = 0; i < 3; ++i) {
			UNROLL_3X(if (swap_dim[i][i_]) {swap(pos[i], pos[i_]);})
		}
	}
	void inv_xform_pos_rm(point &pos) const { // Note: unused/untested
		for (unsigned i = 0; i < 3; ++i) {
			UNROLL_3X(if (swap_dim[2-i][2-i_]) {swap(pos[2-i], pos[2-i_]);})
		}
		UNROLL_3X(if (mirror[i_]) {pos[i_] = -pos[i_];}) // order-independent, inverse has no effect
	}
	void xform_pos_rms(point &pos) const {
		xform_pos_rm(pos);
		pos *= scale;
	}
	void inv_xform_pos_rms(point &pos) const {
		assert(scale != 0.0);
		pos /= scale;
		inv_xform_pos_rm(pos);
	}
	void xform_pos(point &pos) const {
		xform_pos_rms(pos);
		pos += tv;
	}
	void inv_xform_pos(point &pos) const {
		pos -= tv;
		inv_xform_pos_rms(pos);
	}
	cube_t get_xformed_cube_ts(cube_t const &cube) const {return cube*scale + tv;} // Note: RM ignored
	void xform_vect(vector<point> &v) const {
		for (vector<point>::iterator i = v.begin(); i != v.end(); ++i) {xform_pos(*i);}
	}
	bool operator==(geom_xform_t const &x) const;
	bool operator!=(geom_xform_t const &x) const {return !operator==(x);}
};


struct model3d_xform_t : public geom_xform_t { // should be packed, can read/write as POD

	vector3d axis;
	float angle;
	colorRGBA color;
	int tid;

	model3d_xform_t(vector3d const &tv_=zero_vector, float scale_=1.0) : geom_xform_t(tv_, scale_), axis(zero_vector), angle(0.0), color(ALPHA0), tid(-1) {}
	cube_t get_xformed_cube(cube_t const &cube) const;
	void apply_inv_xform_to_pdu(pos_dir_up &pdu) const;
	void apply_gl() const;

	bool eq_xforms(model3d_xform_t const &x) const {
		return (axis == x.axis && angle == x.angle && geom_xform_t::operator==(x));
	}
	bool operator==(model3d_xform_t const &x) const {
		return (eq_xforms(x) && x.color == color && x.tid == tid);
	}
	bool operator!=(model3d_xform_t const &x) const {return !operator==(x);}
	bool is_identity() const {return eq_xforms(model3d_xform_t());}
};


struct vntc_ix_t {
	unsigned vix, nix, tix;
	vntc_ix_t(unsigned vix_=0, unsigned nix_=0, unsigned tix_=0) : vix(vix_), nix(nix_), tix(tix_) {}
};


struct poly_header_t {
	unsigned npts, obj_id;
	int mat_id;
	vector3d n;

	poly_header_t(int mat_id_=-1, unsigned obj_id_=0) : npts(0), mat_id(mat_id_), obj_id(obj_id_), n(zero_vector) {}
};


struct poly_data_block {
	vector<poly_header_t> polys;
	vector<vntc_ix_t> pts;
};


struct model3d_stats_t {
	unsigned verts, quads, tris, blocks, mats, transforms;
	model3d_stats_t() : verts(0), quads(0), tris(0), blocks(0), mats(0), transforms(0) {}
	void print() const;
};


template<typename T> class vertex_map_t : public map<T, unsigned> {

	int last_mat_id;
	unsigned last_obj_id;
	bool average_normals;

public:
	vertex_map_t(bool average_normals_=0) : last_mat_id(-1), last_obj_id(0), average_normals(average_normals_) {}
	bool get_average_normals() const {return average_normals;}
	
	void check_for_clear(int mat_id) {
		if (mat_id != last_mat_id || size() >= MAX_VMAP_SIZE) {
			last_mat_id = mat_id;
			clear();
		}
	}
};

typedef vertex_map_t<vert_norm_tc> vntc_map_t;
typedef vertex_map_t<vert_norm_tc_tan> vntct_map_t;


template<typename T> void clear_cont(T &cont) {T().swap(cont);}


template<typename T> class vntc_vect_t : public vector<T> {

protected:
	bool has_tangents, finalized;
	unsigned vbo, ivbo;
	sphere_t bsphere;
	cube_t bcube;

public:
	unsigned obj_id;

	vntc_vect_t(unsigned obj_id_=0) : has_tangents(0), finalized(0), vbo(0), ivbo(0), obj_id(obj_id_) {}
	void render(shader_t &shader, bool is_shadow_pass, unsigned npts);
	void clear();
	void free_vbos();
	void make_private_copy() {vbo = ivbo = 0;} // Note: to be called *only* after a deep copy
	void add_poly(vntc_vect_t const &poly);
	void calc_bounding_volumes();
	cube_t get_bcube () const {return get_polygon_bbox(*this);}
	point get_center () const {return bsphere.pos;}
	float get_bradius() const {return bsphere.radius;}
	void optimize(unsigned npts) {remove_excess_cap();}
	void remove_excess_cap() {if (20*size() < 19*capacity()) vector<value_type>(*this).swap(*this);}
	void write(ostream &out) const;
	void read(istream &in);
};


template<typename T> class indexed_vntc_vect_t : public vntc_vect_t<T> {

	vector<unsigned> indices;
	bool need_normalize;

	struct geom_block_t {
		unsigned start_ix, num;
		cube_t bcube;
		geom_block_t() : start_ix(0), num(0) {}
		geom_block_t(unsigned s, unsigned n, cube_t const &bc) : start_ix(s), num(n), bcube(bc) {}
	};
	vector<geom_block_t> blocks;

public:
	indexed_vntc_vect_t(unsigned obj_id_=0) : vntc_vect_t(obj_id_), need_normalize(0) {}
	void calc_tangents(unsigned npts) {assert(0);}
	void render(shader_t &shader, bool is_shadow_pass, unsigned npts, bool no_vfc=0);
	void reserve_for_num_verts(unsigned num_verts);
	void add_poly(polygon_t const &poly, vertex_map_t<T> &vmap);
	void add_triangle(triangle const &t, vertex_map_t<T> &vmap);
	unsigned add_vertex(T const &v, vertex_map_t<T> &vmap);
	void add_index(unsigned ix) {assert(ix < size()); indices.push_back(ix);}
	void subdiv_recur(vector<unsigned> const &ixs, unsigned npts, unsigned skip_dims);
	void optimize(unsigned npts);
	void finalize(unsigned npts);
	void simplify(vector<unsigned> &out, float target) const;
	void clear();
	unsigned num_verts() const {return unsigned(indices.empty() ? size() : indices.size());}
	T       &get_vert(unsigned i)       {return (*this)[indices.empty() ? i : indices[i]];}
	T const &get_vert(unsigned i) const {return (*this)[indices.empty() ? i : indices[i]];}
	unsigned get_ix  (unsigned i) const {assert(i < indices.size()); return indices[i];}
	void get_polygons(vector<coll_tquad> &polygons, colorRGBA const &color, unsigned npts, bool quads_only) const;
	void write(ostream &out) const;
	void read(istream &in);
	bool indexing_enabled() const {return !indices.empty();}
	void mark_need_normalize() {need_normalize = 1;}
};


template<typename T> struct vntc_vect_block_t : public deque<indexed_vntc_vect_t<T> > {

	void optimize(unsigned npts);
	void clear() {free_vbos(); deque<indexed_vntc_vect_t<T> >::clear();}
	void free_vbos();
	cube_t get_bcube() const;
	float calc_draw_order_score() const;
	unsigned num_verts() const;
	unsigned num_unique_verts() const;
	void get_stats(model3d_stats_t &stats) const {stats.blocks += (unsigned)size(); stats.verts += num_unique_verts();}
	void get_polygons(vector<coll_tquad> &polygons, colorRGBA const &color, unsigned npts, bool quads_only) const;
	bool write(ostream &out) const;
	bool read(istream &in);
};


template<typename T> struct geometry_t {

	vntc_vect_block_t<T> triangles, quads;

	void calc_tangents_blocks(vntc_vect_block_t<T> &blocks, unsigned npts) {assert(0);}
	void calc_tangents();
	void render_blocks(shader_t &shader, bool is_shadow_pass, vntc_vect_block_t<T> &blocks, unsigned npts);
	void render(shader_t &shader, bool is_shadow_pass);
	bool empty() const {return (triangles.empty() && quads.empty());}
	void add_poly_to_polys(polygon_t const &poly, vntc_vect_block_t<T> &v, vertex_map_t<T> &vmap, unsigned obj_id=0) const;
	void add_poly(polygon_t const &poly, vertex_map_t<T> vmap[2], unsigned obj_id=0);
	void get_polygons(vector<coll_tquad> &polygons, colorRGBA const &color, bool quads_only) const;
	cube_t get_bcube() const;
	void optimize()  {triangles.optimize(3); quads.optimize(4);}
	void free_vbos() {triangles.free_vbos(); quads.free_vbos();}
	void clear();
	void get_stats(model3d_stats_t &stats) const;
	bool write(ostream &out) const {return (triangles.write(out) && quads.write(out));}
	bool read(istream &in)         {return (triangles.read (in ) && quads.read (in ));}
};


class texture_manager {

protected:
	deque<texture_t> textures;
	string_map_t tex_map; // maps texture filenames to texture indexes

public:
	unsigned create_texture(string const &fn, bool is_alpha_mask, bool verbose);
	void clear();
	void free_tids();
	void free_textures();
	void ensure_texture_loaded(texture_t &t, int tid, bool is_bump);
	void bind_alpha_channel_to_texture(int tid, int alpha_tid);
	void ensure_tid_loaded(int tid, bool is_bump) {if (tid >= 0) {ensure_texture_loaded(get_texture(tid), tid, is_bump);}}
	void ensure_tid_bound(int tid) {if (tid >= 0) {get_texture(tid).check_init();}} // if allocated
	void bind_texture(int tid) const {get_texture(tid).bind_gl();}
	colorRGBA get_tex_avg_color(int tid) const {return get_texture(tid).get_avg_color();}
	bool has_binary_alpha(int tid) const {return get_texture(tid).has_binary_alpha;}
	bool might_have_alpha_comp(int tid) const {return (tid >= 0 && get_texture(tid).ncolors == 4);}

	texture_t const &get_texture(int tid) const {
		assert((unsigned)tid < textures.size());
		return textures[tid];
	}
	texture_t &get_texture(int tid) {
		assert((unsigned)tid < textures.size());
		return textures[tid];
	}
};


struct material_params_t {

	colorRGB ka, kd, ks, ke, tf;
	float ns, ni, alpha, tr;
	unsigned illum;
	bool skip, is_used, unused1, unused2; // unused bools to pad the struct

	material_params_t() : ka(WHITE), kd(WHITE), ks(BLACK), ke(BLACK), tf(WHITE),
		ns(1.0), ni(1.0), alpha(1.0), tr(0.0), illum(2), skip(0), is_used(0), unused1(0), unused2(0) {}
}; // must be padded


struct material_t : public material_params_t {

	bool might_have_alpha_comp;
	int a_tid, d_tid, s_tid, alpha_tid, bump_tid, refl_tid;
	float draw_order_score;
	string name, filename;

	geometry_t<vert_norm_tc> geom;
	geometry_t<vert_norm_tc_tan> geom_tan;

	material_t(string const &name_=string(), string const &fn=string())
		: might_have_alpha_comp(0), a_tid(-1), d_tid(-1), s_tid(-1), alpha_tid(-1), bump_tid(-1), refl_tid(-1),
		draw_order_score(0.0), name(name_), filename(fn) {}
	bool add_poly(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], unsigned obj_id=0);
	void mark_as_used() {is_used = 1;}
	bool mat_is_used () const {return is_used;}
	bool use_bump_map() const;
	bool use_spec_map() const;
	void optimize() {geom.optimize(); geom_tan.optimize();}
	int get_render_texture() const {return ((d_tid >= 0) ? d_tid : a_tid);}
	bool get_needs_alpha_test() const {return (alpha_tid >= 0 || might_have_alpha_comp);}
	bool is_partial_transparent() const {return (alpha < 1.0 || get_needs_alpha_test());}
	void init_textures(texture_manager &tmgr);
	void render(shader_t &shader, texture_manager const &tmgr, int default_tid, bool is_shadow_pass, bool enable_alpha_mask);
	colorRGBA get_ad_color() const;
	colorRGBA get_avg_color(texture_manager const &tmgr, int default_tid=-1) const;
	bool write(ostream &out) const;
	bool read(istream &in);
};


struct voxel_params_t; // forward declaration
class voxel_manager; // forward declaration

class model3d {

	// geometry
	geometry_t<vert_norm_tc> unbound_geom;
	int unbound_tid;
	colorRGBA unbound_color;
	vector<polygon_t> split_polygons_buffer;
	cube_t bcube;
	bool from_model3d_file, has_cobjs, needs_alpha_test, needs_bump_maps;

	// materials
	deque<material_t> materials;
	string_map_t mat_map; // maps material names to materials indexes
	set<string> undef_materials; // to reduce warning messages
	cobj_tree_tquads_t coll_tree;

	// transforms
	vector<model3d_xform_t> transforms;

	// shadows
	struct model_smap_data_t : public smap_data_t {
		model3d *model;

		model_smap_data_t(unsigned tu_id_, model3d *model_) : smap_data_t(tu_id_), model(model_) {assert(model);}
		virtual void render_scene_shadow_pass(point const &lpos);
		//virtual bool needs_update(point const &lpos);
	};
	vect_smap_t<model_smap_data_t> smap_data;

	void update_bbox(polygon_t const &poly);

public:
	// textures
	texture_manager &tmgr;

	model3d(texture_manager &tmgr_, int def_tid=-1, colorRGBA const &def_c=WHITE, bool ignore_a=0)
		: tmgr(tmgr_), unbound_tid((def_tid >= 0) ? def_tid : WHITE_TEX), unbound_color(def_c), bcube(all_zeros_cube),
		from_model3d_file(0), has_cobjs(0), needs_alpha_test(0), needs_bump_maps(0) {}
	~model3d() {clear();}
	size_t num_materials(void) const {return materials.size();}

	material_t &get_material(int mat_id) {
		assert(mat_id >= 0 && (unsigned)mat_id < materials.size());
		return materials[mat_id];
	}

	// creation and query
	void set_has_cobjs() {has_cobjs = 1;}
	void add_transform(model3d_xform_t const &xf) {transforms.push_back(xf);}
	unsigned add_triangles(vector<triangle> const &triangles, colorRGBA const &color, int mat_id=-1, unsigned obj_id=0);
	unsigned add_polygon(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], int mat_id=-1, unsigned obj_id=0);
	void add_triangle(polygon_t const &tri, vntc_map_t &vmap, int mat_id=-1, unsigned obj_id=0);
	void get_polygons(vector<coll_tquad> &polygons, bool quads_only=0) const;
	void get_cubes(vector<cube_t> &cubes, float spacing) const;
	int get_material_ix(string const &material_name, string const &fn);
	int find_material(string const &material_name);
	void mark_mat_as_used(int mat_id);
	void optimize();
	void clear();
	void free_context();
	void clear_smaps() {smap_data.clear();} // frees GL state
	void load_all_used_tids();
	void bind_all_used_tids();
	void render_materials_def(shader_t &shader, bool is_shadow_pass, bool enable_alpha_mask, unsigned bmap_pass_mask, xform_matrix const *const mvm=nullptr) {
		render_materials(shader, is_shadow_pass, enable_alpha_mask, bmap_pass_mask, unbound_color, unbound_tid, mvm);
	}
	void render_materials(shader_t &shader, bool is_shadow_pass, bool enable_alpha_mask, unsigned bmap_pass_mask,
		colorRGBA const &cur_ub_color, int cur_ub_tid, xform_matrix const *const mvm=nullptr);
	void render(shader_t &shader, bool is_shadow_pass, bool enable_alpha_mask, unsigned bmap_pass_mask, vector3d const &xlate);
	void setup_shadow_maps();
	bool has_any_transforms() const {return !transforms.empty();}
	cube_t const &get_bcube() const {return bcube;}
	void build_cobj_tree(bool verbose);
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) const;
	bool get_needs_alpha_test() const {return needs_alpha_test;}
	bool get_needs_bump_maps () const {return needs_bump_maps;}
	void get_stats(model3d_stats_t &stats) const;
	void show_stats() const;
	void get_all_mat_lib_fns(set<std::string> &mat_lib_fns) const;
	bool write_to_disk (string const &fn) const;
	bool read_from_disk(string const &fn);
	static void proc_counted_normals(vector<counted_normal> &cn, float nmag_thresh);
};


struct model3ds : public deque<model3d> {

	texture_manager tmgr;

	void clear();
	void free_context();
	void render(bool is_shadow_pass, vector3d const &xlate); // non-const
	bool has_any_transforms() const;
	cube_t get_bcube() const;
	void build_cobj_trees(bool verbose);
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) const;
};


template<typename T> bool split_polygon(polygon_t const &poly, vector<T> &ppts, float coplanar_thresh);

void coll_tquads_from_triangles(vector<triangle> const &triangles, vector<coll_tquad> &ppts, colorRGBA const &color);
void free_model_context();
void render_models(bool shadow_pass, vector3d const &xlate=zero_vector);
void add_transform_for_cur_model(model3d_xform_t const &xf);

bool read_model_file(string const &filename, vector<coll_tquad> *ppts, vector<cube_t> *cubes, cube_t &model_bcube,
	geom_xform_t const &xf, int def_tid, colorRGBA const &def_c, float voxel_xy_spacing, bool load_model_file,
	bool recalc_normals, bool write_file, bool verbose);


#endif // _MODEL3D_H_
