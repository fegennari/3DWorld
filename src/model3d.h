// 3D World - 3D Model Rendering Classes
// by Frank Gennari
// 8/17/11
#pragma once

#include "3DWorld.h"
#include "collision_detect.h" // for polygon_t
#include "cobj_bsp_tree.h" // for cobj_tree_tquads_t
#include "shadow_map.h" // for smap_data_t and rotation_t
#include "gl_ext_arb.h"

#include <unordered_map>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/quaternion.hpp>

using namespace std;

typedef map<string, unsigned> string_map_t;

unsigned const MAX_VMAP_SIZE     = (1 << 18); // 256K
unsigned const BUILTIN_TID_START = (1 << 16); // 65K
unsigned const MAX_NUM_BONES_PER_VERTEX = 4;
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
	void xform_pos(point &pos) const { // rotate, mirror, scale, translate
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
	xform_matrix create_xform_matrix() const;
	bool operator==(geom_xform_t const &x) const;
	bool operator!=(geom_xform_t const &x) const {return !operator==(x);}
};


struct model3d_xform_t : public geom_xform_t, public rotation_t { // should be packed, can read/write as POD

	base_mat_t material;
	int group_cobjs_level=0;
	float voxel_spacing=0.0;
	cube_t bcube_xf;

	model3d_xform_t(vector3d const &tv_=zero_vector, float scale_=1.0) : geom_xform_t(tv_, scale_) {}
	model3d_xform_t(geom_xform_t const &xf) : geom_xform_t(xf) {}
	cube_t get_xformed_cube(cube_t const &cube) const;
	cube_t const &get_xformed_bcube(cube_t const &bcube);
	void clear_bcube() {bcube_xf.set_to_zeros();}
	void apply_inv_xform_to_pdu(pos_dir_up &pdu) const;
	void apply_to_tquad(coll_tquad &tquad) const;
	void apply_gl() const;

	bool eq_xforms(model3d_xform_t const &x) const {return (rotation_t::operator==(x) && geom_xform_t::operator==(x));}
	bool operator==(model3d_xform_t const &x) const {return (eq_xforms(x) && material == x.material);}
	bool operator!=(model3d_xform_t const &x) const {return !operator==(x);}
	bool is_identity() const {return eq_xforms(model3d_xform_t());}

	void xform_pos(point &pos) const { // rotate, mirror, scale, arb_rotate, translate
		xform_pos_rms(pos);
		rotate_point(pos, -1.0); // negative rotate?
		pos += tv;
	}
	void inv_xform_pos(point &pos) const {
		pos -= tv;
		rotate_point(pos, 1.0); // negative rotate?
		inv_xform_pos_rms(pos);
	}
	void apply_material_override(base_mat_t &mat) const {
		if (material.tid >= 0) {mat.tid = material.tid;}
		if (material.shine > 0.0) {mat.shine = material.shine;}
		if (material.color != ALPHA0) {mat.color = material.color;}
		if (material.spec_color != BLACK) {mat.spec_color = material.spec_color;}
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

	poly_header_t(int mat_id_=-1, unsigned obj_id_=0) : npts(0), obj_id(obj_id_), mat_id(mat_id_) {}
};

struct poly_data_block {
	vector<poly_header_t> polys;
	vector<vntc_ix_t> pts;
};

struct model3d_stats_t {
	unsigned verts=0, quads=0, tris=0, blocks=0, mats=0, transforms=0;
	void print() const;
};

// for computing vertex normals from face normals
struct counted_normal : public vector3d { // size = 16
	unsigned count=0;

	counted_normal() {}
	counted_normal(vector3d const &n) : vector3d(n), count(1) {}
	void add_normal(vector3d const &n) {*this += n; ++count;}
	bool is_valid() const {return (count > 0);}
};

// unused
struct weighted_normal : public vector3d { // size = 16
	float weight=0.0;

	weighted_normal() {}
	weighted_normal(vector3d const &n, float w=1.0) : vector3d(w*n), weight(w) {}
	void add_normal(vector3d const &n, float w=1.0) {*this += w*n; weight += w;}
	void add_normal(vector3d const &n, point const *const poly_pts, unsigned npts) {add_normal(n, polygon_area(poly_pts, npts));}
	bool is_valid() const {return (weight > 0.0);}
};

//template<typename T> class vertex_map_t : public unordered_map<T, unsigned, hash_by_bytes<T>> {
template<typename T> class vertex_map_t : public map<T, unsigned> {
	int last_mat_id=-1;
	unsigned last_obj_id=0;
	bool average_normals=0;

public:
	vertex_map_t(bool average_normals_=0) : average_normals(average_normals_) {}
	bool get_average_normals() const {return average_normals;}
	
	void check_for_clear(int mat_id) {
		if (mat_id != last_mat_id || this->size() >= MAX_VMAP_SIZE) {
			last_mat_id = mat_id;
			this->clear();
		}
	}
};

typedef vertex_map_t<vert_norm_tc> vntc_map_t;
typedef vertex_map_t<vert_norm_tc_tan> vntct_map_t;


struct get_polygon_args_t {
	vector<coll_tquad> &polygons;
	bool quads_only=0;
	unsigned lod_level=0;

	get_polygon_args_t(vector<coll_tquad> &polygons_, bool quads_only_=0, unsigned lod_level_=0)
		: polygons(polygons_), quads_only(quads_only_), lod_level(lod_level_) {}
};


template<typename T> cube_t get_polygon_bbox(vector<T> const &p) {
	if (p.empty()) return all_zeros_cube;
	cube_t bbox(p.front().v, p.front().v);
	for (unsigned i = 1; i < p.size(); ++i) {bbox.union_with_pt(p[i].v);}
	return bbox;
}


struct vertex_bone_data_t { // Note: must be packed
	unsigned ids    [MAX_NUM_BONES_PER_VERTEX] = {};
	float    weights[MAX_NUM_BONES_PER_VERTEX] = {};
	void add(unsigned id, float weight, bool &had_vertex_error);
	void normalize();
};
struct mesh_bone_data_t {
	vector<vertex_bone_data_t> vertex_to_bones;
};

struct model_anim_t {
	unordered_map<string, unsigned> bone_name_to_index_map;
	vector<xform_matrix> bone_transforms, bone_offset_matrices;
	xform_matrix global_inverse_transform, root_transform;
	string model_name; // for debug printouts
	mutable bool had_anim_id_error=0;

	struct anim_node_t {
		int bone_index; // cached to avoid bone_name_to_index_map lookup; -1 is no bone
		mutable bool no_anim_data=0;
		string name;
		xform_matrix transform;
		vector<unsigned> children; // indexes into anim_nodes
		anim_node_t(string const &name_, xform_matrix const &transform_, unsigned bone_index_) : bone_index(bone_index_), name(name_), transform(transform_) {}
	};
	vector<anim_node_t> anim_nodes;

	template<typename T> struct anim_val_t {
		float time=0.0;
		T v;
		anim_val_t() {} // for merged resize()
		anim_val_t(float time_, T const &v_) : time(time_), v(v_) {}
	};
	struct merged_anim_data_t {
		vector3d pos, scale;
		glm::quat q=glm::quat(1.0, 0.0, 0.0, 0.0); // initialize to identity
		void set(vector3d const &pos_, vector3d const &scale_, glm::quat const &q_) {pos = pos_; scale = scale_; q = q_;}
	};
	struct anim_data_t {
		bool uses_scale=0;
		vector<anim_val_t<vector3d>> pos, scale;
		vector<anim_val_t<glm::quat>> rot;
		vector<anim_val_t<merged_anim_data_t>> merged;
		void init(unsigned np, unsigned nr, unsigned ns);
	};
	struct animation_t {
		float ticks_per_sec=25.0, duration=1.0; // duration is in ticks
		string name;
		unordered_map<string, anim_data_t> anim_data; // per bone
		animation_t(string const &name_="") : name(name_) {}
	};
	vector<animation_t> animations;

	unsigned get_bone_id(string const &bone_name);
	void transform_node_hierarchy_recur(float anim_time, animation_t const &animation, unsigned node_ix, xform_matrix const &parent_transform);
	void get_bone_transforms(unsigned anim_id, float cur_time);
	bool check_anim_wrapped(unsigned anim_id, float old_time, float new_time) const;
	float get_anim_duration(unsigned anim_id) const;
private:
	xform_matrix apply_anim_transform(float anim_time, animation_t const &animation, anim_node_t const &node) const;
public:
	void blend_animations_simple(unsigned anim_id1, unsigned anim_id2, float blend_factor, float cur_time1, float cur_time2);
	void blend_animations(unsigned anim_id1, unsigned anim_id2, float blend_factor, float delta_time, float &cur_time1, float &cur_time2);
	void get_blended_bone_transforms(float anim_time1, float anim_time2, animation_t const &animation1, animation_t const &animation2,
		unsigned node_ix, xform_matrix const &parent_transform, float blend_factor);
	void merge_anim_transforms();
	void merge_from(model_anim_t const &anim);
	int get_animation_id_by_name(string const &anim_name) const;
};


template<typename T> class vntc_vect_t : public vector<T>, public indexed_vao_manager_with_shadow_t {

protected:
	bool has_tangents=0, finalized=0;
	sphere_t bsphere;
	cube_t bcube;
public:
	using vector<T>::empty;
	using vector<T>::size;
	unsigned obj_id;

	vntc_vect_t(unsigned obj_id_=0) : obj_id(obj_id_) {bcube.set_to_zeros();}
	void clear();
	void make_private_copy() {vbo = ivbo = 0;} // Note: to be called *only* after a deep copy
	void calc_bounding_volumes();
	void ensure_bounding_volumes() {if (bsphere.radius == 0.0) {calc_bounding_volumes();}}
	cube_t get_bcube () const {return get_polygon_bbox(*this);}
	point get_center () const {return bsphere.pos;}
	float get_bradius() const {return bsphere.radius;}
	size_t get_gpu_mem() const {return (vbo_valid() ? size()*sizeof(T) : 0);}
	void optimize(unsigned npts) {remove_excess_cap();}
	void remove_excess_cap() {if (20*vector<T>::size() < 19*vector<T>::capacity()) {vector<T>::shrink_to_fit();}}
	void write(ostream &out) const;
	void read(istream &in);
};


template<typename T> class indexed_vntc_vect_t : public vntc_vect_t<T> {
public:
	typedef unsigned index_type_t;
	vector<unsigned> indices; // needs to be public for merging operation
	mesh_bone_data_t bone_data;
	bool has_bones() const {return !bone_data.vertex_to_bones.empty();}
private:
	bool need_normalize, optimized, prev_ucc;
	float avg_area_per_tri, amin, amax;

	struct geom_block_t {
		unsigned start_ix=0, num=0;
		cube_t bcube;
		geom_block_t() {}
		geom_block_t(unsigned s, unsigned n, cube_t const &bc) : start_ix(s), num(n), bcube(bc) {}
	};
	vector<geom_block_t> blocks;

	struct lod_block_t {
		unsigned start_ix=0, num=0;
		float tri_area=0.0;
		lod_block_t() {}
		lod_block_t(unsigned s, unsigned n, float a) : start_ix(s), num(n), tri_area(a) {}
		unsigned get_end_ix() const {return (start_ix + num);}
	};
	vector<lod_block_t> lod_blocks;
	unsigned get_block_ix(float area) const;

public:
	using vntc_vect_t<T>::size;
	using vntc_vect_t<T>::empty;
	using vntc_vect_t<T>::begin;
	using vntc_vect_t<T>::end;
	using vntc_vect_t<T>::at;
	using vntc_vect_t<T>::operator[];
	using vntc_vect_t<T>::finalized;
	using vntc_vect_t<T>::bcube;
	using vntc_vect_t<T>::bsphere;
	
	indexed_vntc_vect_t(unsigned obj_id_=0) : vntc_vect_t<T>(obj_id_), need_normalize(0), optimized(0), prev_ucc(0), avg_area_per_tri(0.0), amin(0.0), amax(0.0) {}
	void calc_tangents(unsigned npts) {assert(0);}
	void setup_bones(shader_t &shader, bool is_shadow_pass);
	void unset_bone_attrs();
	void render(shader_t &shader, bool is_shadow_pass, point const *const xlate, unsigned npts, bool no_vfc=0);
	void reserve_for_num_verts(unsigned num_verts);
	void add_poly(polygon_t const &poly, vertex_map_t<T> &vmap);
	void add_triangle(triangle const &t, vertex_map_t<T> &vmap);
	unsigned add_vertex(T const &v, vertex_map_t<T> &vmap);
	void add_index(unsigned ix) {assert(ix < size()); indices.push_back(ix);}
	void subdiv_recur(vector<unsigned> const &ixs, unsigned npts, unsigned skip_dims, cube_t *bcube_in=nullptr);
	void optimize(unsigned npts);
	void gen_lod_blocks(unsigned npts);
	void finalize(unsigned npts);
	void finalize_lod_blocks(unsigned npts);
	void simplify(vector<unsigned> &out, float target) const;
	void simplify_meshoptimizer(vector<unsigned> &out, float target) const;
	void simplify_indices(float reduce_target);
	void reverse_winding_order(unsigned npts);
	void clear();
	void clear_blocks() {blocks.clear(); lod_blocks.clear();}
	unsigned num_verts() const {return unsigned(indices.empty() ? size() : indices.size());}
	T       &get_vert(unsigned i)       {return (*this)[indices.empty() ? i : indices[i]];}
	T const &get_vert(unsigned i) const {return (*this)[indices.empty() ? i : indices[i]];}
	unsigned get_ix  (unsigned i) const {assert(i < indices.size()); return indices[i];}
	float get_prim_area(unsigned i, unsigned npts) const;
	float calc_area(unsigned npts);
	void get_polygons(get_polygon_args_t &args, unsigned npts) const;
	size_t get_gpu_mem() const {return (vntc_vect_t<T>::get_gpu_mem() + (this->ivbo_valid() ? indices.size()*sizeof(unsigned) : 0));}
	void invert_tcy();
	void write(ostream &out) const;
	void read(istream &in, unsigned npts);
	void write_to_obj_file(ostream &out, unsigned &cur_vert_ix, unsigned npts) const;
	bool indexing_enabled() const {return !indices.empty();}
	void mark_need_normalize() {need_normalize = 1;}
};


template<typename T> struct vntc_vect_block_t : public deque<indexed_vntc_vect_t<T> > {

	using deque<indexed_vntc_vect_t<T> >::begin;
	using deque<indexed_vntc_vect_t<T> >::end;
	
	void finalize(unsigned npts);
	void clear() {free_vbos(); deque<indexed_vntc_vect_t<T> >::clear();}
	void free_vbos();
	cube_t get_bcube() const;
	size_t get_gpu_mem() const;
	float calc_draw_order_score() const;
	unsigned num_verts() const;
	unsigned num_unique_verts() const;
	void get_stats(model3d_stats_t &stats) const {stats.blocks += (unsigned)this->size(); stats.verts += num_unique_verts();}
	float calc_area(unsigned npts);
	void get_polygons(get_polygon_args_t &args, unsigned npts) const;
	void invert_tcy();
	void simplify_indices(float reduce_target);
	void reverse_winding_order(unsigned npts);
	void merge_into_single_vector();
	bool write(ostream &out) const;
	bool read(istream &in, unsigned npts);
	bool write_to_obj_file(ostream &out, unsigned &cur_vert_ix, unsigned npts) const;
};


template<typename T> struct geometry_t {

	vntc_vect_block_t<T> triangles, quads;

	void calc_tangents_blocks(vntc_vect_block_t<T> &blocks, unsigned npts) {assert(0);}
	void calc_tangents();
	void render_blocks(shader_t &shader, bool is_shadow_pass, point const *const xlate, vntc_vect_block_t<T> &blocks, unsigned npts);
	void render(shader_t &shader, bool is_shadow_pass, point const *const xlate);
	bool empty() const {return (triangles.empty() && quads.empty());}
	size_t get_gpu_mem() const {return (triangles.get_gpu_mem() + quads.get_gpu_mem());}
	unsigned add_triangles(vector<vert_norm_tc> const &verts, vector<unsigned> const &indices, bool add_new_block);
	void add_poly_to_polys(polygon_t const &poly, vntc_vect_block_t<T> &v, vertex_map_t<T> &vmap, unsigned obj_id=0) const;
	void add_poly(polygon_t const &poly, vertex_map_t<T> vmap[2], unsigned obj_id=0);
	void get_polygons(get_polygon_args_t &args) const;
	cube_t get_bcube() const;
	void invert_tcy() {triangles.invert_tcy(); quads.invert_tcy();}
	void finalize  () {triangles.finalize(3); quads.finalize(4);}
	void free_vbos () {triangles.free_vbos(); quads.free_vbos();}
	void clear();
	void get_stats(model3d_stats_t &stats) const;
	void calc_area(float &area, unsigned &ntris);
	void simplify_indices(float reduce_target);
	void reverse_winding_order();
	bool write(ostream &out) const {return (triangles.write(out)  && quads.write(out)) ;}
	bool read(istream &in)         {return (triangles.read(in, 3) && quads.read(in, 4));}
	bool write_to_obj_file(ostream &out, unsigned &cur_vert_ix) const {return (triangles.write_to_obj_file(out, cur_vert_ix, 3) && quads.write_to_obj_file(out, cur_vert_ix, 4));}
};


class texture_manager {
	struct tex_work_item_t {
		unsigned tid;
		bool is_nm;
		tex_work_item_t(unsigned tid_, bool is_nm_) : tid(tid_), is_nm(is_nm_) {}
		// we really shouldn't have the same tid with different is_nm; should we just ingore is_nm when comparing
		bool operator< (tex_work_item_t const &w) const {return ((tid == w.tid) ? (is_nm < w.is_nm) : (tid < w.tid));}
		bool operator==(tex_work_item_t const &w) const {return (tid == w.tid && is_nm == w.is_nm);}
	};
	static size_t tot_textures, tot_gpu_mem, tot_cpu_mem; // summed across all texture managers
protected:
	deque<texture_t> textures;
	string_map_t tex_map; // maps texture filenames to texture indexes
	vector<tex_work_item_t> to_load;
public:
	unsigned create_texture(string const &fn, bool is_alpha_mask, bool verbose,
		bool invert_alpha=0, bool wrap=1, bool mirror=0, bool force_grayscale=0, bool is_nm=0, bool invert_y=0, bool no_cache=0, bool load_now=0);
	bool empty() const {return textures.empty();}
	void clear();
	void free_tids();
	void free_textures();
	void free_client_mem();
	void remove_last_texture();
	bool ensure_texture_loaded(int tid, bool is_bump);
	void post_load_texture_from_memory(int tid);
	void bind_alpha_channel_to_texture(int tid, int alpha_tid);
	bool ensure_tid_loaded(int tid, bool is_bump) {return ((tid >= 0) ? ensure_texture_loaded(tid, is_bump) : 0);}
	void ensure_tid_bound(int tid);
	void add_work_item(int tid, bool is_nm);
	void load_work_items_mt();
	void bind_texture(int tid) const {get_texture(tid).bind_gl();}
	void bind_texture_tu_or_white_tex(int tid, unsigned tu_id) const;
	colorRGBA get_tex_avg_color(int tid) const {return get_texture(tid).get_avg_color();}
	bool has_binary_alpha(int tid) const {return get_texture(tid).has_binary_alpha;}
	bool might_have_alpha_comp(int tid) const {return (tid >= 0 && get_texture(tid).ncolors == 4);}
	texture_t const &get_texture(int tid) const;
	texture_t &get_texture(int tid);
	size_t get_cpu_mem() const;
	size_t get_gpu_mem() const;
	static size_t get_tot_textures() {return tot_textures;}
	static size_t get_tot_gpu_mem () {return tot_gpu_mem ;}
	static size_t get_tot_cpu_mem () {return tot_cpu_mem ;}
};


struct material_params_t { // Warning: changing this struct will invalidate the model3d file format

	colorRGB ka=WHITE, kd=WHITE, ks=BLACK, ke=BLACK, tf=WHITE;
	float ns=1.0, ni=1.0, alpha=1.0, tr=0.0;
	unsigned illum=2;
	bool skip=0, is_used=0, no_blend=0, unused_field=0; // unused bool to pad the struct
}; // must be padded


struct material_t : public material_params_t {

	bool might_have_alpha_comp=0, tcs_checked=0, no_lod_cull=0;
	int a_tid=-1, d_tid=-1, s_tid=-1, ns_tid=-1, alpha_tid=-1, bump_tid=-1, refl_tid=-1; // ambient, diffuse, specular, shininess, transparency, normal, reflection (unused)
	float draw_order_score=0.0, avg_area_per_tri=0.0, tot_tri_area=0.0;
	float metalness=-1.0; // < 0 disables; should go into material_params_t, but that would invalidate the model3d file format
	string name, filename;

	geometry_t<vert_norm_tc> geom;
	geometry_t<vert_norm_tc_tan> geom_tan;

	material_t(string const &name_="", string const &fn="") : name(name_), filename(fn) {}
	bool empty() const {return (geom.empty() && geom_tan.empty());}
	mesh_bone_data_t &get_bone_data_for_last_added_tri_mesh();
	unsigned add_triangles(vector<vert_norm_tc> const &verts, vector<unsigned> const &indices, bool add_new_block); // Note: no quads or tangents
	bool add_poly(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], unsigned obj_id=0);
	void mark_as_used() {is_used = 1;}
	bool mat_is_used () const {return is_used;}
	bool use_bump_map() const;
	bool use_spec_map() const;
	unsigned get_gpu_mem() const {return (geom.get_gpu_mem() + geom_tan.get_gpu_mem());}
	void finalize() {geom.finalize(); geom_tan.finalize();}
	int get_render_texture() const {return ((d_tid >= 0 || a_tid < 0) ? d_tid : a_tid);} // return diffuse texture unless ambient texture is specified but diffuse texture is not
	bool get_needs_alpha_test  () const {return (alpha_tid >= 0 || might_have_alpha_comp);}
	bool is_partial_transparent() const {return ((alpha < 1.0 || get_needs_alpha_test()) && !no_blend);}
	bool has_alpha_mask        () const {return (alpha_tid >= 0 && get_render_texture() >= 0);}
	void compute_area_per_tri();
	void simplify_indices(float reduce_target);
	void reverse_winding_order();
	void ensure_textures_loaded(texture_manager &tmgr);
	void init_textures(texture_manager &tmgr);
	void queue_textures_to_load(texture_manager &tmgr);
	void check_for_tc_invert_y(texture_manager &tmgr);
	void render(shader_t &shader, texture_manager const &tmgr, int default_tid, bool is_shadow_pass, bool is_z_prepass,
		int enable_alpha_mask, bool is_bmap_pass, point const *const xlate, bool no_set_min_alpha=0);
	colorRGBA get_ad_color() const;
	colorRGBA get_avg_color(texture_manager const &tmgr, int default_tid=-1) const;
	bool write(ostream &out) const;
	bool read(istream &in);
	bool write_to_obj_file(ostream &out, unsigned &cur_vert_ix) const;
	void write_mtllib_entry(ostream &out, texture_manager const &tmgr) const;
};


struct voxel_params_t; // forward declaration
class voxel_manager; // forward declaration

class model3d {
	// read/write options
	string filename;
	int recalc_normals=0, group_cobjs_level=0;

	// geometry
	geometry_t<vert_norm_tc> unbound_geom;
	base_mat_t unbound_mat;
	vector<polygon_t> split_polygons_buffer;
	cube_t bcube, bcube_all_xf, occlusion_cube;
	unsigned model_refl_tid=0, model_refl_tsize=0, model_refl_last_tsize=0, model_indir_tid=0;
	int reflective=0; // reflective: 0=none, 1=planar, 2=cube map
	int indoors=2; // 0=no/outdoors, 1=yes/indoors, 2=unknown
	bool from_model3d_file=0, has_cobjs=0, needs_alpha_test=0, needs_bump_maps=0, has_spec_maps=0, has_gloss_maps=0, xform_zvals_set=0, needs_trans_pass=0, has_alpha_mask=0;
	float metalness=0.0; // should be per-material, but not part of the material file and specified per-object instead
	float lod_scale=1.0;

	// materials
	deque<material_t> materials;
	string_map_t mat_map; // maps material names to materials indexes
	set<string> undef_materials; // to reduce warning messages
	cobj_tree_tquads_t coll_tree;
	colorRGBA cached_avg_color=ALPHA0; // used by get_and_cache_avg_color()
	bool textures_loaded=0;

	// transforms
	vector<model3d_xform_t> transforms;

	// shadows
	struct model_smap_data_t : public smap_data_t {
		model3d *model;

		model_smap_data_t(unsigned tu_id_, unsigned smap_sz_, model3d *model_) : smap_data_t(tu_id_, smap_sz_), model(model_) {assert(model);}
		virtual void render_scene_shadow_pass(point const &lpos);
		//virtual bool needs_update(point const &lpos);
	};
	typedef vect_smap_t<model_smap_data_t> per_model_smap_data;
	map<rotation_t, per_model_smap_data> smap_data;

	// lighting
	string sky_lighting_fn;
	unsigned sky_lighting_sz[3];
	float sky_lighting_weight=0.0;
	//lmap_manager_t local_lmap_manager;

	// temporaries to be reused
	vector<pair<float, unsigned> > to_draw, to_draw_xf;

	void update_bbox(polygon_t const &poly);
	void create_indir_texture();
public:
	texture_manager &tmgr; // stores all textures
	model_anim_t model_anim_data;

	model3d(string const &filename_, texture_manager &tmgr_, int def_tid=-1, colorRGBA const &def_c=WHITE, int reflective_=0,
		float metalness_=0.0, float lod_scale_=1.0, int recalc_normals_=0, int group_cobjs_level_=0)
		: filename(filename_), recalc_normals(recalc_normals_), group_cobjs_level(group_cobjs_level_), unbound_mat(((def_tid >= 0) ? def_tid : WHITE_TEX), def_c),
		reflective(reflective_), metalness(metalness_), lod_scale(lod_scale_), tmgr(tmgr_)
	{UNROLL_3X(sky_lighting_sz[i_] = 0;)}
	~model3d() {clear();}
	size_t num_materials() const {return materials.size();}
	material_t &get_material(int mat_id, bool alloc_if_needed=0);
	base_mat_t const &get_unbound_material() const {return unbound_mat;}

	// creation and query
	bool empty              () const {return (materials.empty() && unbound_geom.empty());}
	bool are_textures_loaded() const {return textures_loaded;}
	string const &get_filename() const {return filename;}
	void set_has_cobjs() {has_cobjs = 1;}
	void add_transform(model3d_xform_t const &xf) {transforms.push_back(xf);}
	unsigned add_polygon(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], int mat_id=-1, unsigned obj_id=0);
	void add_triangle(polygon_t const &tri, vntc_map_t &vmap, int mat_id=-1, unsigned obj_id=0);
	void get_polygons(vector<coll_tquad> &polygons, bool quads_only=0, bool apply_transforms=0, unsigned lod_level=0) const;
	void get_transformed_bcubes(vector<cube_t> &bcubes) const;
	void get_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf) const;
	colorRGBA get_avg_color() const;
	colorRGBA get_area_weighted_avg_color();
	colorRGBA get_and_cache_avg_color(bool area_weighted=0);
	size_t get_gpu_mem() const;
	int get_material_ix(string const &material_name, string const &fn, bool okay_if_exists=0);
	int find_material(string const &material_name);
	void mark_mat_as_used(int mat_id);
	void set_xform_zval_from_tt_height(bool flatten_mesh);
	void finalize();
	void clear();
	void free_context();
	void clear_smaps(); // frees GL state
	void load_all_used_tids();
	void bind_all_used_tids();
	void calc_tangent_vectors();
	void simplify_indices(float reduce_target);
	void reverse_winding_order(uint64_t mats_mask=~uint64_t(0));
	void set_sky_lighting_file(string const &fn, float weight, unsigned sz[3]);
	void set_occlusion_cube(cube_t const &cube) {occlusion_cube = cube;}
	geom_xform_t fit_to_scene();
	void set_target_translate_scale(point const &target_pos, float target_radius, geom_xform_t &xf) const;
	void render_materials_def(shader_t &shader, bool is_shadow_pass, int reflection_pass, bool is_z_prepass, int enable_alpha_mask,
		unsigned bmap_pass_mask, int trans_op_mask, point const *const xlate, xform_matrix const *const mvm=nullptr)
	{
		render_materials(shader, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, trans_op_mask, unbound_mat, rotation_t(), xlate, mvm);
	}
	void render_materials(shader_t &shader, bool is_shadow_pass, int reflection_pass, bool is_z_prepass, int enable_alpha_mask, unsigned bmap_pass_mask,
		int trans_op_mask, base_mat_t const &unbound_mat, rotation_t const &rot, point const *const xlate=nullptr, xform_matrix const *const mvm=nullptr,
		bool force_lod=0, float model_lod_mult=1.0, float fixed_lod_dist=0.0, bool skip_cull_face=0, bool is_scaled=0, bool no_set_min_alpha=0, unsigned skip_mat_mask=0);
	void render_material(shader_t &shader, unsigned mat_id, bool is_shadow_pass, bool is_z_prepass=0, int enable_alpha_mask=0, bool is_bmap_pass=0,
		point const *const xlate=nullptr, bool no_set_min_alpha=0);
	void render_with_xform(shader_t &shader, model3d_xform_t &xf, xform_matrix const &mvm, bool is_shadow_pass,
		int reflection_pass, bool is_z_prepass, int enable_alpha_mask, unsigned bmap_pass_mask, int reflect_mode, int trans_op_mask);
	void render(shader_t &shader, bool is_shadow_pass, int reflection_pass, bool is_z_prepass, int enable_alpha_mask,
		unsigned bmap_pass_mask, int reflect_mode, int trans_op_mask, vector3d const &xlate);
	material_t *get_material_by_name(string const &name);
	colorRGBA set_color_for_material(unsigned mat_id, colorRGBA const &color);
	int set_texture_for_material(unsigned mat_id, int tid);
	void set_material_emissive_color(unsigned mat_id, colorRGBA const &color);
	void ensure_reflection_cube_map();
	cube_t get_single_transformed_bcube(vector3d const &xlate=zero_vector) const;
	void setup_shadow_maps();
	void set_local_model_scene_bounds(shader_t &s);
	bool has_any_transforms() const {return !transforms.empty();}
	cube_t const &get_bcube() const {return bcube;}
	cube_t calc_bcube_including_transforms();
	void union_bcube_with(cube_t const &c) {bcube.assign_or_union_with_cube(c);}
	void build_cobj_tree(bool verbose);
	bool check_coll_line_cur_xf(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact);
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact, bool build_bvh_if_needed=0);
	bool get_needs_alpha_test() const {return needs_alpha_test;}
	bool get_needs_trans_pass() const {return needs_trans_pass;}
	bool get_has_alpha_mask  () const {return has_alpha_mask;}
	bool get_needs_bump_maps () const {return needs_bump_maps;}
	bool uses_spec_map()        const {return has_spec_maps;}
	bool uses_gloss_map()       const {return has_gloss_maps;}
	bool is_planar_reflective() const {return (reflective == 1);}
	bool is_cube_map_reflective() const {return (reflective == 2);}
	bool is_reflective()        const {return (reflective != 0);}
	void compute_area_per_tri();
	void get_stats(model3d_stats_t &stats) const;
	void show_stats() const;
	void get_all_mat_lib_fns(set<std::string> &mat_lib_fns) const;
	bool write_to_disk (string const &fn) const;
	bool read_from_disk(string const &fn);
	bool write_as_obj_file(string const &fn);
	static void proc_model_normals(vector<counted_normal> &cn, int recalc_normals, float nmag_thresh=0.7);
	static void proc_model_normals(vector<weighted_normal> &wn, int recalc_normals, float nmag_thresh=0.7);
	void write_to_cobj_file(std::ostream &out) const;

	// animations
	unsigned num_animations() const {return model_anim_data.animations.size();}
	bool has_animations    () const {return (num_animations() > 0);}
	void setup_bone_transforms(shader_t &shader, float anim_time, int anim_id=-1);
	void setup_bone_transforms_cached(bone_transform_data_t &cached, shader_t &shader, float anim_time, int anim_id=-1);
	void setup_bone_transforms_blended(shader_t &shader, float anim_time1, float anim_time2, float blend_factor, int anim_id1=-1, int anim_id2=-1);
	bool check_anim_wrapped(unsigned anim_id, float old_time, float new_time) const;
	float get_anim_duration(unsigned anim_id) const;
	void merge_animation_from(model3d const &anim_model) {model_anim_data.merge_from(anim_model.model_anim_data);}
protected:
	unsigned get_anim_id(shader_t &shader, string const &prop_name, int anim_id=-1) const;
	void add_bone_transforms_to_shader(shader_t &shader, vector<xform_matrix> const &bone_transforms) const;
	void add_bone_transforms_to_shader(shader_t &shader) const {add_bone_transforms_to_shader(shader, model_anim_data.bone_transforms);}
};


struct model3ds : public deque<model3d> {
	texture_manager tmgr;

	void clear();
	void free_context();
	void render(bool is_shadow_pass, int reflection_pass, int trans_op_mask, vector3d const &xlate); // non-const
	void ensure_reflection_cube_maps();
	void set_xform_zval_from_tt_height(bool flatten_mesh);
	bool has_any_transforms() const;
	bool has_any_animations() const;
	cube_t calc_and_return_bcube(bool only_reflective);
	void get_all_model_bcubes(vector<cube_t> &bcubes) const;
	size_t get_gpu_mem() const;
	void build_cobj_trees(bool verbose);
	bool check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact, bool build_bvh_if_needed=0);
	void write_to_cobj_file(std::ostream &out) const;
};


class model_from_file_t {
	string rel_path;
protected:
	model3d &model;

public:
	model_from_file_t(string const &fn, model3d &model_) : model(model_) {rel_path = get_path(fn);}
	string open_include_file(string const &fn, string const &type, ifstream &in_inc) const;
	static string get_path(string const &fn);
	int get_texture(string const &fn, bool is_alpha_mask, bool verbose, bool invert_alpha=0, bool wrap=1, bool mirror=0, bool force_grayscale=0);
	void check_and_bind(int &tid, string const &tfn, bool is_alpha_mask, bool verbose, bool invert_alpha=0, bool wrap=1, bool mirror=0);
};

template<typename T> bool split_polygon(polygon_t const &poly, vector<T> &ppts, float coplanar_thresh, bool allow_quads=1);

bool use_model3d_bump_maps();
void coll_tquads_from_triangles(vector<triangle> const &triangles, vector<coll_tquad> &ppts, colorRGBA const &color);
void free_model_context();
void render_models(int shadow_pass, int reflection_pass, int trans_op_mask=3, vector3d const &xlate=zero_vector);
void ensure_model_reflection_cube_maps();
void auto_calc_model_zvals();
void get_cur_model_polygons(vector<coll_tquad> &ppts, model3d_xform_t const &xf=model3d_xform_t(), unsigned lod_level=0);
size_t get_loaded_models_gpu_mem();
void get_cur_model_edges_as_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf);
void get_cur_model_as_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf);
bool add_transform_for_cur_model(model3d_xform_t const &xf);
void set_sky_lighting_file_for_cur_model(string const &fn, float weight, unsigned sz[3]);
void set_occlusion_cube_for_cur_model(cube_t const &cube);
geom_xform_t fit_cur_model_to_scene();
bool have_cur_model();
cube_t calc_and_return_all_models_bcube(bool only_reflective=0);
void get_all_model_bcubes(vector<cube_t> &bcubes);
void write_models_to_cobj_file(std::ostream &out);
void adjust_zval_for_model_coll(point &pos, float radius, float mesh_zval, float step_height=0.0);
void check_legal_movement_using_model_coll(point const &prev, point &cur, float radius=0.0);

bool load_model_file(string const &filename, model3ds &models, geom_xform_t const &xf, string const &anim_name, int def_tid, colorRGBA const &def_c,
	int reflective, float metalness, float lod_scale, int recalc_normals, int group_cobjs_level, bool write_file, bool verbose, uint64_t rev_winding_mask=0);
bool read_model_file(string const &filename, vector<coll_tquad> *ppts, geom_xform_t const &xf, int def_tid, colorRGBA const &def_c,
	int reflective, float metalness, float lod_scale, bool load_model_file, int recalc_normals, int group_cobjs_level, bool write_file, bool verbose);

