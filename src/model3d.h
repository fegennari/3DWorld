// 3D World - 3D Model Rendering Classes
// by Frank Gennari
// 8/17/11

#ifndef _MODEL3D_H_
#define _MODEL3D_H_

#include "3DWorld.h"
#include "collision_detect.h" // for polygon_t

using namespace std;

typedef map<string, unsigned> string_map_t;

unsigned const MAX_VMAP_SIZE     = (1 << 18); // 256K
float const POLY_COPLANAR_THRESH = 0.98;
colorRGB const def_color(0.0, 0.0, 0.0);


struct geom_xform_t {
	vector3d tv;
	float scale;
	bool mirror[3], swap_dim[3][3];

	geom_xform_t() : tv(zero_vector), scale(1.0) {
		for (unsigned i = 0; i < 3; ++i) {
			UNROLL_3X(swap_dim[i][i_] = 0;)
			mirror[i] = 0;
		}
	}
	void xform_pos_rm(point &pos) const {
		UNROLL_3X(if (mirror[i_]) pos[i_] = -pos[i_];)
		
		for (unsigned i = 0; i < 3; ++i) {
			UNROLL_3X(if (swap_dim[i][i_]) swap(pos[i], pos[i_]);)
		}
	}
	void xform_pos_rms(point &pos) const {
		xform_pos_rm(pos);
		pos *= scale;
	}
	void xform_pos(point &pos) const {
		xform_pos_rms(pos);
		pos += tv;
	}
	void xform_vect(vector<point> &v) const {
		for (vector<point>::iterator i = v.begin(); i != v.end(); ++i) {
			xform_pos(*i);
		}
	}
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
	unsigned verts, quads, tris, blocks, mats;
	model3d_stats_t() : verts(0), quads(0), tris(0), blocks(0), mats(0) {}
	void print() const {cout << "verts: " << verts << ", quads: " << quads << ", tris: " << tris << ", blocks: " << blocks << ", mats: " << mats << endl;}
};


struct vert_norm_tc_tan : public vert_norm_tc { // size = 48
	vector4d tangent;

	vert_norm_tc_tan() {}
	vert_norm_tc_tan(vert_norm_tc const &vntc, vector4d const &tangent_=vector4d(0,0,0,0))
		: vert_norm_tc(vntc), tangent(tangent_) {}
	vert_norm_tc_tan(point const &v_, vector3d const &n_, float ts, float tt, vector4d const &tangent_=vector4d(0,0,0,0))
		: vert_norm_tc(v_, n_, ts, tt), tangent(tangent_) {}

	bool operator<(vert_norm_tc_tan const &p) const {
		if (v < p.v) return 1;
		if (p.v < v) return 0;
		if (n < p.n) return 1;
		if (p.n < n) return 0;
		if (t[0] < p.t[0]) return 1;
		if (p.t[0] < t[0]) return 0;
		if (t[1] < p.t[1]) return 1;
		if (p.t[1] < t[1]) return 0;
		return (tangent < p.tangent);
	}
};


template<typename T> class vertex_map_t : public map<T, unsigned> {

	int last_mat_id;
	unsigned last_obj_id;

public:
	vertex_map_t() : last_mat_id(-1), last_obj_id(0) {}
	
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
	bool has_tangents;
	unsigned vbo, ivbo;
	sphere_t bsphere;
	cube_t bcube;

public:
	unsigned obj_id;

	vntc_vect_t(unsigned obj_id_=0) : has_tangents(0), vbo(0), ivbo(0), obj_id(obj_id_) {}
	void render(shader_t &shader, bool is_shadow_pass, int prim_type);
	void free_vbos();
	void add_poly(vntc_vect_t const &poly);
	void calc_bounding_volumes();
	cube_t get_bbox() const {return get_polygon_bbox(*this);}
	void remove_excess_cap() {if (20*size() < 19*capacity()) vector<value_type>(*this).swap(*this);}
	void write(ostream &out) const;
	void read(istream &in);
};


template<typename T> class indexed_vntc_vect_t : public vntc_vect_t<T> {

	vector<unsigned> indices;

public:
	indexed_vntc_vect_t(unsigned obj_id_=0) : vntc_vect_t(obj_id_) {}
	void calc_tangents(unsigned npts) {assert(0);}
	void render(shader_t &shader, bool is_shadow_pass, int prim_type);
	void add_poly(polygon_t const &poly, vertex_map_t<T> &vmap);
	void add_vertex(T const &v, vertex_map_t<T> &vmap);
	void clear() {vntc_vect_t<T>::clear(); indices.clear();}
	unsigned num_verts() const {return (indices.empty() ? size() : indices.size());}
	T       &get_vert(unsigned i)       {return (*this)[indices.empty() ? i : indices[i]];}
	T const &get_vert(unsigned i) const {return (*this)[indices.empty() ? i : indices[i]];}
	void get_polygons(vector<polygon_t> &polygons, colorRGBA const &color, unsigned npts) const;
	void write(ostream &out) const;
	void read(istream &in);
};


template<typename T> struct vntc_vect_block_t : public deque<indexed_vntc_vect_t<T> > {

	void remove_excess_cap();
	void free_vbos();
	cube_t get_bbox() const;
	unsigned num_verts() const;
	unsigned num_unique_verts() const;
	void get_stats(model3d_stats_t &stats) const {stats.blocks += size(); stats.verts += num_unique_verts();}
	void get_polygons(vector<polygon_t> &polygons, colorRGBA const &color, unsigned npts) const;
	bool write(ostream &out) const;
	bool read(istream &in);
};


template<typename T> struct geometry_t {

	vntc_vect_block_t<T> triangles, quads;

	void calc_tangents() {assert(0);}
	void render(shader_t &shader, bool is_shadow_pass);
	bool empty() const {return (triangles.empty() && quads.empty());}
	void add_poly_to_polys(polygon_t const &poly, vntc_vect_block_t<T> &v, vertex_map_t<T> &vmap, unsigned obj_id=0) const;
	void add_poly(polygon_t const &poly, vertex_map_t<T> vmap[2], unsigned obj_id=0);
	void get_polygons(vector<polygon_t> &polygons, colorRGBA const &color) const;
	cube_t get_bbox() const;
	void remove_excess_cap() {triangles.remove_excess_cap(); quads.remove_excess_cap();}
	void free_vbos()         {triangles.free_vbos(); quads.free_vbos();}
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
	void ensure_tid_loaded(int tid, bool is_bump);
	void ensure_tid_bound(int tid);
	void bind_texture(int tid) const;
	colorRGBA get_tex_avg_color(int tid) const;
};


struct material_params_t {

	colorRGB ka, kd, ks, ke, tf;
	float ns, ni, alpha, tr;
	unsigned illum;
	bool skip, is_used;

	material_params_t() : ka(def_color), kd(def_color), ks(def_color), ke(def_color), tf(def_color),
		ns(1.0), ni(1.0), alpha(1.0), tr(0.0), illum(2), skip(0), is_used(0) {}
};


struct material_t : public material_params_t {

	int a_tid, d_tid, s_tid, alpha_tid, bump_tid, refl_tid;
	string name, filename;

	geometry_t<vert_norm_tc> geom;
	geometry_t<vert_norm_tc_tan> geom_tan;

	material_t(string const &name_=string(), string const &fn=string())
		: a_tid(-1), d_tid(-1), s_tid(-1), alpha_tid(-1), bump_tid(-1), refl_tid(-1), name(name_), filename(fn) {}
	bool add_poly(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], unsigned obj_id=0);
	void mark_as_used() {is_used = 1;}
	bool mat_is_used () const {return is_used;}
	bool use_bump_map() const;
	bool use_spec_map() const;
	int get_render_texture() const {return ((d_tid >= 0) ? d_tid : a_tid);}
	bool is_partial_transparent() const {return (alpha < 1.0 || alpha_tid >= 0);}
	void render(shader_t &shader, texture_manager const &tmgr, int default_tid, bool ignore_ambient, bool is_shadow_pass);
	colorRGBA get_ad_color() const;
	colorRGBA get_avg_color(texture_manager const &tmgr, int default_tid=-1) const;
	bool write(ostream &out) const;
	bool read(istream &in);
};


class cobj_tree_tquads_t;


class model3d {

	// geometry
	geometry_t<vert_norm_tc> unbound_geom;
	int unbound_tid;
	colorRGBA unbound_color;
	vector<polygon_t> split_polygons_buffer;
	cube_t bbox;
	bool from_model3d_file, ignore_ambient, has_cobjs;

	// materials
	deque<material_t> materials;
	string_map_t mat_map; // maps material names to materials indexes
	set<string> undef_materials; // to reduce warning messages
	cobj_tree_tquads_t *coll_tree;

public:
	// textures
	texture_manager &tmgr;

	model3d(texture_manager &tmgr_, int def_tid=-1, colorRGBA const &def_c=WHITE, bool ignore_a=0) : tmgr(tmgr_),
		unbound_tid((def_tid >= 0) ? def_tid : WHITE_TEX), unbound_color(def_c), bbox(all_zeros_cube),
		from_model3d_file(0), ignore_ambient(ignore_a), has_cobjs(0), coll_tree(NULL) {}
	~model3d() {clear();}
	unsigned num_materials(void) const {return materials.size();}

	material_t &get_material(int mat_id) {
		assert(mat_id >= 0 && (unsigned)mat_id < materials.size());
		return materials[mat_id];
	}

	// creation and query
	void set_has_cobjs() {has_cobjs = 1;}
	unsigned add_polygon(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], int mat_id, unsigned obj_id=0);
	void get_polygons(vector<polygon_t> &polygons) const;
	int get_material_ix(string const &material_name, string const &fn);
	int find_material(string const &material_name);
	void mark_mat_as_used(int mat_id);
	void remove_excess_cap();
	void clear();
	void free_context();
	void load_all_used_tids();
	void bind_all_used_tids();
	void render(shader_t &shader, bool is_shadow_pass, bool bmap_pass); // const?
	cube_t const &get_bbox() const {return bbox;}
	void build_cobj_tree(bool verbose);
	void free_cobj_tree();
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) const;
	void get_stats(model3d_stats_t &stats) const;
	void show_stats() const;
	void get_all_mat_lib_fns(set<std::string> &mat_lib_fns) const;
	bool write_to_disk (string const &fn) const;
	bool read_from_disk(string const &fn);
};


struct model3ds : public deque<model3d> {

	texture_manager tmgr;

	void clear();
	void free_context();
	void render(bool is_shadow_pass); // const?
	cube_t get_bbox() const;
	void build_cobj_trees(bool verbose);
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) const;
};


bool split_polygon(polygon_t const &poly, vector<polygon_t> &ppts, float coplanar_thresh);

void free_model_context();
void render_models(bool shadow_pass);

bool read_object_file(string const &filename, vector<polygon_t> *ppts, geom_xform_t const &xf, int def_tid,
	colorRGBA const &def_c, bool load_model_file, bool recalc_normals, bool write_file, bool ignore_ambient, bool verbose);


#endif // _MODEL3D_H_
