// 3D World - 3D Model Rendering Code
// by Frank Gennari
// 8/17/11

#include "model3d.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "voxels.h"
#include "vertex_opt.h"
#include "voxels.h" // for get_cur_model_edges_as_cubes
#include "csg.h" // for clip_polygon_to_cube
#include "lightmap.h" // for lmap_manager_t
#include <fstream>
#include <queue>
#include "meshoptimizer.h"
#include "profiler.h"
#include "format_text.h"
#include "binary_file_io.h"

#include <glm/gtc/matrix_transform.hpp>

bool const ENABLE_BUMP_MAPS  = 1;
bool const ENABLE_SPEC_MAPS  = 1;
bool const ENABLE_INTER_REFLECTIONS = 1;
bool const SHOW_MODEL_BCUBE_CENTER  = 0;
bool const ENABLE_ANIMATION_SHADOWS = 1;
bool const USE_ANIM_MODEL_TANGENTS  = 1;
unsigned const MAGIC_NUMBER      = 42987143; // arbitrary file signature; no   animations
unsigned const MAGIC_NUMBER_ANIM = 42987144; // arbitrary file signature; with animations
unsigned const BLOCK_SIZE        = 32768; // in vertex indices
unsigned const BONE_IDS_LOC      = 4;
unsigned const BONE_WEIGHTS_LOC  = 5;

bool model_calc_tan_vect(1); // slower and more memory but sometimes better quality/smoother transitions

extern bool group_back_face_cull, enable_model3d_tex_comp, disable_shader_effects, texture_alpha_in_red_comp, use_model3d_tex_mipmaps, enable_model3d_bump_maps;
extern bool two_sided_lighting, have_indir_smoke_tex, use_core_context, model3d_wn_normal, invert_model_nmap_bscale, use_z_prepass, all_model3d_ref_update;
extern bool use_interior_cube_map_refl, enable_model3d_custom_mipmaps, enable_tt_model_indir, no_subdiv_model, auto_calc_tt_model_zvals, use_model_lod_blocks;
extern bool flatten_tt_mesh_under_models, no_store_model_textures_in_memory, disable_model_textures, allow_model3d_quads, merge_model_objects, invert_model3d_faces;
extern unsigned shadow_map_sz, reflection_tid;
extern int display_mode, animate2, default_anim_id;
extern float model3d_alpha_thresh, model3d_texture_anisotropy, model_triplanar_tc_scale, model_mat_lod_thresh, cobj_z_bias, model_hemi_lighting_scale, light_int_scale[];
extern double tfticks;
extern pos_dir_up orig_camera_pdu;
extern bool vert_opt_flags[3];
extern vector<texture_t> textures;


model3ds all_models;


bool enable_bump_map() {return (ENABLE_BUMP_MAPS && !disable_shader_effects && enable_model3d_bump_maps);} // enabled by default
bool enable_spec_map() {return (ENABLE_SPEC_MAPS && !disable_shader_effects);}
bool no_sparse_smap_update();
bool enable_reflection_dynamic_updates();
string texture_str(int tid);
bool use_model3d_bump_maps() {return enable_bump_map();} // global function export
void draw_transparent_city_bldg_geom(int reflection_pass, vector3d const &xlate);


// ************ texture_manager ************

size_t texture_manager::tot_textures = 0;
size_t texture_manager::tot_gpu_mem  = 0;
size_t texture_manager::tot_cpu_mem  = 0;

unsigned texture_manager::create_texture(string const &fn, bool is_alpha_mask, bool verbose, bool invert_alpha,
	bool wrap, bool mirror, bool force_grayscale, bool is_nm, bool invert_y, bool no_cache, bool load_now)
{
	assert(!(wrap && mirror)); // can't both be set

	if (!no_cache) { // temp images aren't in the texture cache
		string_map_t::const_iterator it(tex_map.find(fn));

		if (it != tex_map.end()) { // found (already loaded)
			assert(it->second < textures.size());
			return it->second; // check invert_alpha?
		}
	}
	unsigned const tid((unsigned)textures.size());
	if (!no_cache) {tex_map[fn] = tid;}
	if (verbose) cout << "creating texture " << fn << endl;
	bool const compress(!is_alpha_mask && enable_model3d_tex_comp);
	// Note: custom mipmap code is also present in ensure_texture_loaded(), but that may not be called early enough if load_now=0
	int const use_mipmaps((use_model3d_tex_mipmaps && !is_alpha_mask) ? (enable_model3d_custom_mipmaps ? 4 : 1) : 0);
	unsigned ncolors((is_alpha_mask || force_grayscale) ? 1 : 3);
	// type=read_from_file format=auto width height wrap_mir ncolors use_mipmaps name [do_compress]
	// always RGB wrapped+mipmap (normal map flag set later)
	textures.emplace_back(0, IMG_FMT_AUTO, 0, 0, (mirror ? 2 : (wrap ? 1 : 0)), ncolors, use_mipmaps, fn, invert_y, compress, model3d_texture_anisotropy, 1.0, is_nm);
	textures.back().invert_alpha = invert_alpha;
	if (load_now) {ensure_texture_loaded(tid, is_nm);} // must load temp images now
	++tot_textures;
	return tid; // can't fail
}

void texture_manager::clear() {
	tot_textures -= textures.size();
	free_textures();
	textures.clear();
	tex_map.clear();
}
void texture_manager::free_tids() {
	for (auto &t : textures) {tot_gpu_mem -= t.get_gpu_mem(); t.gl_delete();}
}
void texture_manager::free_textures() {
	free_tids();
	free_client_mem();
}
void texture_manager::free_client_mem() { // Note: should not be called if model textures can overlap with predefined textures
	for (auto &t : textures) {tot_cpu_mem -= t.get_cpu_mem(); t.free_client_mem();}
}

void texture_manager::remove_last_texture() {
	assert(!textures.empty());
	texture_t &t(textures.back());
	assert(!t.is_bound());
	tot_cpu_mem -= t.get_cpu_mem();
	--tot_textures;
	t.free_client_mem(); // should not be allocated?
	tex_map.erase(t.name);
	textures.pop_back();
}

bool texture_manager::ensure_texture_loaded(int tid, bool is_bump) {

	texture_t &t(get_texture(tid));
	if (t.is_loaded()) return 0; // already loaded from disk
	if (t.is_bound() ) return 0; // already bound to a texture/sent to the GPU, no need to reload
	//if (is_bump) {t.do_compress = 0;} // don't compress normal maps
	// Note: it's incorrect to call t.has_alpha() here because that uses color, which hasn't been computed yet (t.init() is called later);
	// but that's okay, do_gl_init() will disable custom mipmaps for textures with color.A == 1.0
	if (use_model3d_tex_mipmaps && enable_model3d_custom_mipmaps /*&& t.has_alpha()*/) {t.use_mipmaps = 4;}
	t.load(-1);
	assert(t.is_loaded());
		
	if (t.alpha_tid >= 0 && t.alpha_tid != tid) { // if alpha is the same texture then the alpha channel should already be set
		ensure_tid_loaded(t.alpha_tid, 0);
		t.copy_alpha_from_texture(get_texture(t.alpha_tid), texture_alpha_in_red_comp);
	}
	if (is_bump && t.ncolors == 1) {t.make_normal_map();} // make RGB normal map from grayscale bump map
	t.init(); // must be after alpha copy
	tot_cpu_mem += t.get_cpu_mem(); // not thread safe?
	return 1;
}

void texture_manager::post_load_texture_from_memory(int tid) { // called on internal texture formats; what about external textures, do we need expand_grayscale_to_rgb()?
	texture_t &t(get_texture(tid));
	assert(t.is_loaded());
	t.expand_grayscale_to_rgb();
	t.fix_word_alignment(); // untested
	t.init(); // calls calc_color()
	tot_cpu_mem += t.get_cpu_mem(); // not thread safe?
}

void texture_manager::bind_alpha_channel_to_texture(int tid, int alpha_tid) {

	if (tid < 0 || alpha_tid < 0) return; // no texture
	assert((unsigned)alpha_tid < textures.size());
	texture_t &t(get_texture(tid));
	assert(t.ncolors == 3 || t.ncolors == 4);
	if (t.alpha_tid == alpha_tid) return; // already bound
	assert(tid < (int)BUILTIN_TID_START); // can't modify builtin textures
	assert(t.alpha_tid < 0); // can't rebind to a different value
	assert(!t.is_allocated()); // must not yet be loaded
	t.alpha_tid = alpha_tid;
	t.ncolors   = 4; // add alpha channel
	if (t.use_mipmaps) {t.use_mipmaps = 3;} // generate custom alpha mipmaps
}

void texture_manager::ensure_tid_bound(int tid) { // if allocated
	if (tid < 0) return; // not setup
	texture_t &t(get_texture(tid));
	if (t.is_bound()) return; // already updated
	t.check_init();
	tot_gpu_mem += t.get_gpu_mem();
}

void texture_manager::add_work_item(int tid, bool is_nm) {
	if (tid < 0) return;
	if (get_texture(tid).is_loaded() || get_texture(tid).is_bound()) return; // already loaded or bound, nothing to do
	to_load.emplace_back(tid, is_nm);
}
void texture_manager::load_work_items_mt() {
	if (to_load.empty()) return; // nothing to do
	sort_and_unique(to_load);
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < (int)to_load.size(); ++i) {ensure_texture_loaded(to_load[i].tid, to_load[i].is_nm);}
	to_load.clear();
}

texture_t &get_builtin_texture(int tid) {
	assert((unsigned)tid < textures.size());
	return textures[tid];
}
texture_t const &texture_manager::get_texture(int tid) const {
	if (tid >= (int)BUILTIN_TID_START) {return get_builtin_texture(tid - BUILTIN_TID_START);} // global textures lookup
	assert((unsigned)tid < textures.size());
	return textures[tid]; // local textures lookup
}
texture_t &texture_manager::get_texture(int tid) {
	if (tid >= (int)BUILTIN_TID_START) {return get_builtin_texture(tid - BUILTIN_TID_START);} // global textures lookup
	assert((unsigned)tid < textures.size());
	return textures[tid]; // local textures lookup
}

size_t texture_manager::get_cpu_mem() const {
	size_t mem(0);
	for (auto t = textures.begin(); t != textures.end(); ++t) {mem += t->get_cpu_mem();}
	return mem;
}
size_t texture_manager::get_gpu_mem() const {
	size_t mem(0);
	for (auto t = textures.begin(); t != textures.end(); ++t) {mem += t->get_gpu_mem();}
	return mem;
}


// ************ vntc_vect_t/indexed_vntc_vect_t ************

// explicit template instantiations of vert_norm case, used for voxel_model, where tc=0.0
template class indexed_vntc_vect_t<vert_norm>;


template<typename T> void vntc_vect_t<T>::clear() {
	
	clear_vbos();
	vector<T>::clear();
	finalized = has_tangents = 0;
	bsphere.radius = 0.0;
}


template<typename T> void vntc_vect_t<T>::calc_bounding_volumes() {

	assert(!empty());
	bsphere.pos = zero_vector;
	for (unsigned i = 0; i < size(); ++i) {bsphere.pos += (*this)[i].v;}
	bsphere.pos /= size();
	bsphere.radius = 0.0;
	bcube = cube_t(bsphere.pos, bsphere.pos);
	
	for (unsigned i = 0; i < size(); ++i) {
		bsphere.radius = max(bsphere.radius, p2p_dist_sq(bsphere.pos, (*this)[i].v));
		bcube.union_with_pt((*this)[i].v);
	}
	bsphere.radius = sqrt(bsphere.radius);
}


template<typename T> void indexed_vntc_vect_t<T>::subdiv_recur(vector<unsigned> const &ixs, unsigned npts, unsigned skip_dims, cube_t *bcube_in) {

	unsigned const num(ixs.size());
	assert(num > 0 && (num % npts) == 0);
	cube_t bc;
	
	if (bcube_in != nullptr) {bc = *bcube_in;} // use cached bcube
	else { // compute bcube
		bc.set_from_point(at(ixs.front()).v);

		for (vector<unsigned>::const_iterator i = ixs.begin()+1; i != ixs.end(); ++i) { // slow due to cache misses (indices are not in order)
			bc.union_with_pt(at(*i).v); // update bounding cube
		}
	}
	if (num > BLOCK_SIZE) { // subdiv case
		float max_sz(0), sval(0);
		unsigned const dim(bc.get_split_dim(max_sz, sval, skip_dims));

		if (max_sz > 0) { // can split
			vector<unsigned> bins[2];
			for (unsigned i = 0; i < 2; ++i) {bins[i].reserve(3*num/5);} // reserve to 60%

			for (unsigned i = 0; i < num; i += npts) {
				vector<unsigned> &dest(bins[(at(ixs[i]).v[dim] > sval)]); // use the first point to determine the bin
				for (unsigned j = i; j < i+npts; ++j) {dest.push_back(ixs[j]);}
			}
			if (bins[0].empty() || bins[1].empty()) {skip_dims |= (1 << dim);}

			for (unsigned i = 0; i < 2; ++i) {
				if (!bins[i].empty()) {subdiv_recur(bins[i], npts, skip_dims);}
			}
			return;
		}
	}
	blocks.emplace_back(indices.size(), num, bc);
	copy(ixs.begin(), ixs.end(), back_inserter(indices)); // make leaf
}


template<typename T> void indexed_vntc_vect_t<T>::optimize(unsigned npts) {

	if (optimized) return;
	optimized = 1;
	vntc_vect_t<T>::optimize(npts);

	if (vert_opt_flags[0]) { // only if not subdivided?
		vert_optimizer optimizer(indices, size(), npts);
		optimizer.run(vert_opt_flags[1], vert_opt_flags[2]);
	}
}


unsigned get_area_pow2(float area, float amin) {return unsigned(log2(max(area/amin, 1.0f)));} // truncate

template<typename T> unsigned indexed_vntc_vect_t<T>::get_block_ix(float area) const {
	unsigned const ix(get_area_pow2(area, amin));
	assert(ix < lod_blocks.size());
	return (lod_blocks.size() - ix - 1);
}

template<typename T> void indexed_vntc_vect_t<T>::gen_lod_blocks(unsigned npts) {

	//timer_t timer("Gen LOD Blocks");
	unsigned const num(indices.size()), num_prims(num/npts);
	assert(num > 0 && (num % npts) == 0);

	// compute min/max area, and use this to determine the number of LOD blocks
	for (unsigned i = 0; i < num_prims; ++i) {
		float const area(get_prim_area(i*npts, npts));
		if (i == 0) {amin = amax = area;} else {amin = min(amin, area); amax = max(amax, area);}
	}
	amin = max(amin, amax/1024.0f); // limit to a reasonable number of blocks
	assert(amin > 0.0);
	unsigned const num_blocks(get_area_pow2(amax, amin) + 1);
	//cout << TXT(amin) << TXT(amax) << TXT(num_blocks) << endl;
	if (num_blocks == 1) return; // all triangles have similar area, don't subdivide
	lod_blocks.resize(num_blocks);

	// count and record the number of triangles in each block and start index for each block
	for (unsigned i = 0; i < num_prims; ++i) {
		lod_blocks[get_block_ix(get_prim_area(i*npts, npts))].num += npts;
	}
	for (unsigned i = 0; i < num_blocks-1; ++i) {lod_blocks[i+1].start_ix = lod_blocks[i].get_end_ix();}
	//cout << "start: "; for (unsigned i = 0; i < num_blocks; ++i) {cout << lod_blocks[i].start_ix << " ";} cout << endl;
	//cout << "num  : "; for (unsigned i = 0; i < num_blocks; ++i) {cout << lod_blocks[i].num << " ";} cout << endl;
	assert(lod_blocks.back().get_end_ix() == indices.size());
	for (unsigned i = 0; i < num_blocks; ++i) {lod_blocks[i].num = 0;} // reset for next pass

	// reorder indices based on LOD blocks
	vector<unsigned> ixs(indices.size());

	for (unsigned i = 0; i < num_prims; ++i) {
		unsigned const ix(get_block_ix(get_prim_area(i*npts, npts)));
		unsigned const cur_pos(lod_blocks[ix].get_end_ix());
		assert(cur_pos + npts <= ((ix+1 == num_blocks) ? indices.size() : lod_blocks[ix+1].start_ix));
		for (unsigned n = 0; n < npts; ++n) {ixs[cur_pos + n] = indices[i*npts + n];} // copy indices
		lod_blocks[ix].num += npts;
	}
	indices.swap(ixs);
}

template<typename T> void indexed_vntc_vect_t<T>::finalize(unsigned npts) { // Note: called when reading obj files, not model3d files

	optimize(npts);

	if (need_normalize) {
		for (auto i = begin(); i != end(); ++i) {i->n.normalize();}
		need_normalize = 0;
	}
	if (!empty()) {this->ensure_bounding_volumes();}
	if (indices.empty() || finalized) return; // nothing to do
#if 0 // TESTING, maybe this doesn't really go here
	if (npts == 3) { // Note: triangles only
		unsigned const init_num_indices(indices.size());
		for (unsigned n = 0; n < 3; ++n) {
			unsigned const prev_num(indices.size());
			simplify_indices(0.5);
			if (indices.size() == prev_num) break; // no change, done
		}
		out << "*** SIMPLIFY: input=" << init_num_indices << " output=" << indices.size() << endl;
	}
#endif
	finalized = 1;
	finalize_lod_blocks(npts);
}

template<typename T> void indexed_vntc_vect_t<T>::finalize_lod_blocks(unsigned npts) {

	assert((num_verts() % npts) == 0); // triangles or quads
	assert(blocks.empty() && lod_blocks.empty());

	if (use_model_lod_blocks && indices.size() > 1024) {
		gen_lod_blocks(npts);
	}
	else if (!no_subdiv_model && num_verts() > 2*BLOCK_SIZE) { // subdivide large buffers
		//timer_t timer("Subdivide Model");
		vector<unsigned> ixs;
		ixs.swap(indices);
		subdiv_recur(ixs, npts, 0, &bcube);
	}
}


template<unsigned N> struct vert_to_tri_t {
	unsigned t[N] = {0}, n=0; // if this vertex is used in more than N triangles we give up and never remove it

	void add(unsigned ix) {if (n < N) {t[n] = ix;} ++n;} // only add if it fits, but always increment n
	void remove(unsigned tix) {assert(tix < min(n, N)); t[tix] = t[n-1]; --n;} // move last element to position tix
	unsigned get_first_index_ix(unsigned tix) const {assert(tix < min(n, N)); return 3*t[tix];} // multiply by 3 to convert from triangle to index
	bool ix_overflow() const {return (n > N);}
};

struct merge_entry_t {
	unsigned vix;
	float val;

	merge_entry_t(unsigned vix_=0, float val_=0.0) : vix(vix_), val(val_) {}
	bool operator<(merge_entry_t const &e) const {return (val < e.val);}
};

struct vertex_remap_t {

	vector<unsigned> remap;

	vertex_remap_t(unsigned num_verts) {
		remap.resize(num_verts);
		for (unsigned i = 0; i < num_verts; ++i) {remap[i] = i;}
	}
	unsigned get_remapped_val(unsigned ix) { // union find with path compression
		if (is_remapped(ix)) {remap[ix] = get_remapped_val(remap[ix]);} // follow the chain
		return remap[ix];
	}
	bool is_remapped(unsigned ix) const {assert(ix < remap.size()); return (remap[ix] != ix);}
	void remap_vertex(unsigned from, unsigned to) {assert(from != to); assert(!is_remapped(from)); remap[from] = to;}
};

struct mesh_edge_t {
	unsigned a, b; // a < b
	mesh_edge_t(unsigned a_=0, unsigned b_=0) : a(min(a_, b_)), b(max(a_, b_)) {assert(a != b);}
	bool operator==(mesh_edge_t const &e) const {return (a == e.a && b == e.b);}
};


// target = ratio of output to input vertices in (0.0, 1.0)
// Note: works on triangles only (not quads), intended for 2-manifold meshes
template<typename T> void indexed_vntc_vect_t<T>::simplify(vector<unsigned> &out, float target) const {

	RESET_TIME;
	assert(target < 1.0 && target > 0.0);
	out.clear();
	unsigned const num_verts(size()), num_ixs(indices.size()), target_num_verts(unsigned(target*num_verts));
	if (target_num_verts <= 3) {out = indices; return;} // can't simplify

	// build vertex to face/triangle mapping
	assert((num_ixs % 3) == 0); // must be triangles; num_tris = num_ixs/3
	vector<vert_to_tri_t<8>> vert_to_tri(num_verts);

	for (unsigned i = 0; i < num_ixs; ++i) {
		assert(indices[i] < num_verts);
		vert_to_tri[indices[i]].add(i);
	}
	// determine which verts/edges can be removed
	std::priority_queue<merge_entry_t> merge_queue; // of vertices

	for (unsigned i = 0; i < num_verts; ++i) {
		auto const &vt(vert_to_tri[i]);
		if (vt.ix_overflow()) continue; // don't remove this vertex
		bool on_mesh_edge(0);
		counted_normal normal_sum;

		for (unsigned t = 0; t < vt.n && !on_mesh_edge; ++t) { // iterate over triangles at this vertex
			assert(vt.t[t] < num_ixs);
			unsigned const tix(vt.t[t]/3), six(3*tix);
			assert(six+3 <= num_ixs);
			
			for (unsigned j = 0; j < 3; ++j) {
				normal_sum.add_normal(operator[](indices[six+j]).n);
				mesh_edge_t const e1(indices[six+j], indices[six+((j+1)%3)]);
				if (e1.a != i && e1.b != i) continue; // edge doesn't contain the vertex of interest
				bool found(0);

				for (unsigned t2 = 0; t2 < vt.n && !found; ++t2) {
					if (t2 == t) continue; // same triangle
					unsigned const tix2(vt.t[t2]/3), six2(3*tix2);
					if (e1 == mesh_edge_t(indices[six2+0], indices[six2+1]) ||
						e1 == mesh_edge_t(indices[six2+1], indices[six2+2]) ||
						e1 == mesh_edge_t(indices[six2+2], indices[six2+0])) {found = 1; break;}
				} // for t2
				if (!found) {on_mesh_edge = 1; break;}
			} // for j
		} // for t
		if (on_mesh_edge) continue; // can't remove this vertex
		float const val(normal_sum.mag()/normal_sum.count);
		merge_queue.push(merge_entry_t(i, val));
	}
	// mapping from orig vertex to collapsed (new) vertex
	unsigned const cand_verts(merge_queue.size());
	vertex_remap_t remap(num_verts);
	unsigned num_valid_verts(num_verts);

	while (!merge_queue.empty() && num_valid_verts > target_num_verts) {
		unsigned const src_ix(merge_queue.top().vix); // remove this one (half edge collapse)
		merge_queue.pop();
		assert(src_ix < num_verts);
		unsigned dest_ix(src_ix);
		float min_edge_len_sq(0.0);
		
		for (unsigned t = 0; t < vert_to_tri[src_ix].n; ++t) { // iterate over triangles at this vertex
			unsigned const tix(vert_to_tri[src_ix].t[t]/3), six(3*tix);

			for (unsigned j = 0; j < 3; ++j) {
				if (remap.is_remapped(indices[six+j])) continue; // already remapped
				float const edge_len_sq(p2p_dist_sq(operator[](src_ix).v, operator[](indices[six+j]).v));

				if (min_edge_len_sq == 0.0 || edge_len_sq < min_edge_len_sq) { // first or shortest edge
					dest_ix = indices[six+j]; // new vertex
					min_edge_len_sq = edge_len_sq;
				}
			}
		} // for i
		if (min_edge_len_sq == 0.0) continue; // can't remove this one
		remap.remap_vertex(src_ix, dest_ix);
		assert(num_valid_verts > 0);
		--num_valid_verts;
	} // while

	// generate output
	out.reserve(unsigned(target*num_ixs));

	for (unsigned i = 0; i < num_ixs; i += 3) { // iterate by triangle
		unsigned new_ixs[3] = {};
		for (unsigned n = 0; n < 3; ++n) {new_ixs[n] = remap.get_remapped_val(indices[i+n]);}
		if (new_ixs[0] == new_ixs[1] || new_ixs[1] == new_ixs[2] || new_ixs[2] == new_ixs[0]) continue; // duplicate vertices, degenerate triangle, skip
		UNROLL_3X(out.push_back(new_ixs[i_]);)
	}
	cout << TXT(num_verts) << TXT(num_ixs) << TXT(target_num_verts) << TXT(cand_verts) << TXT(num_valid_verts) << TXT(out.size()) << endl;
	PRINT_TIME("Simplify");
}

template<typename T> void indexed_vntc_vect_t<T>::simplify_meshoptimizer(vector<unsigned> &out, float target) const { // triangles only

	timer_t timer("Meshoptimizer Simplify");
	assert(target < 1.0 && target > 0.0);
	float const target_error = 0.01;
	unsigned const num_verts(size()), num_ixs(indices.size()), target_num_ixs(max(3U, unsigned(target*num_ixs)));
	if (num_ixs < 3) return; // no triangles to simplify
	out.resize(num_ixs); // allocate space
	size_t const num_ixs_out(meshopt_simplify(out.data(), indices.data(), num_ixs, &this->front().v.x, num_verts, sizeof(T), target_num_ixs, target_error));
	cout << TXT(num_ixs) << TXT(target_num_ixs) << TXT(num_ixs_out) << endl;
	assert(num_ixs_out > 0 && num_ixs_out <= num_ixs);
	out.resize(num_ixs_out); // truncate to correct size
}

template<typename T> void indexed_vntc_vect_t<T>::simplify_indices(float reduce_target) {
	vector<unsigned> simplified_indices;
	simplify_meshoptimizer(simplified_indices, reduce_target);
	indices.swap(simplified_indices);
}

template<typename T> void indexed_vntc_vect_t<T>::reverse_winding_order(unsigned npts) {
	if (indices.empty()) return; // unsupported (error?)
	unsigned const nverts(num_verts());
	assert((nverts%npts) == 0);
	for (unsigned i = 0; i < nverts; i += npts) {reverse(indices.begin()+i, indices.begin()+i+npts);}
	for (auto &v : *this) {v.invert_normal();}
}

template<typename T> void indexed_vntc_vect_t<T>::clear() {
	vntc_vect_t<T>::clear();
	indices.clear();
	clear_blocks();
	need_normalize = 0;
}

void ensure_valid_tangent(vector4d &tangent) {
	if ((vector3d)tangent == zero_vector) {tangent.assign(0.0, 0.0, 1.0, tangent.w);}
}

template<> void indexed_vntc_vect_t<vert_norm_tc_tan>::calc_tangents(unsigned npts) {

	if (has_tangents) return; // already computed
	has_tangents = 1;
	assert(npts >= 3); // at least triangles
	unsigned const nverts(num_verts());
	assert((nverts%npts) == 0);

	for (unsigned i = 0; i < nverts; i += npts) {
		vert_norm_tc const &A(get_vert(i)), &B(get_vert(i+1)), &C(get_vert(i+2));
		vector3d const v1(A.v - B.v), v2(C.v - B.v);
		float const t1(A.t[1] - B.t[1]), t2(C.t[1] - B.t[1]), s1(A.t[0] - B.t[0]), s2(C.t[0] - B.t[0]);
		float const val(s1*t2 - s2*t1), w((val < 0.0) ? -1.0 : 1.0);
		vector4d tangent((v1*t2 - v2*t1).get_norm(), w);
		ensure_valid_tangent(tangent);
		for (unsigned j = i; j < i+npts; ++j) {get_vert(j).tangent += tangent;}
	}
	for (iterator i = begin(); i != end(); ++i) { // need to renormalize tangents
		i->tangent.normalize();
		i->tangent.w = ((i->tangent.w < 0.0) ? -1.0 : 1.0);
		ensure_valid_tangent(i->tangent);
	}
}


template<typename T> void vntc_vect_t<T>::write(ostream &out) const {
	write_vector(out, *this);
}

template<typename T> void vntc_vect_t<T>::read(istream &in) {
	// Note: it would be nice to write/read without the tangent vectors and recalculate them later,
	// but knowing which materials require tangents requires loading the material file first, but that requires the model,
	// so we would have to read the model3d material headers, then read the material file, then read the polygon data into the correct geometry type,
	// which would also require writing out the model3d file in two passes and smaller blocks of data at a time
	read_vector(in, *this);
	has_tangents = (sizeof(T) == sizeof(vert_norm_tc_tan)); // HACK to get the type
	calc_bounding_volumes();
}

void vertex_bone_data_t::add(unsigned id, float weight, bool &had_vertex_error) {
	assert(weight > 0.0);
	float min_weight(weight);
	unsigned min_weight_ix(0);

	for (unsigned i = 0; i < MAX_NUM_BONES_PER_VERTEX; ++i) {
		if (weights[i] < min_weight) {
			min_weight    = weights[i];
			min_weight_ix = i;
			if (min_weight == 0.0) break; // zero weight is the min - done
		}
	}
	if (min_weight > 0.0 && !had_vertex_error) {
		cerr << "Warning: too many weights/bones for a single vertex; using the " << MAX_NUM_BONES_PER_VERTEX << " largest weights" << endl;
		had_vertex_error = 1;
	}
	if (min_weight == weight) return; // no slot for this weight
	ids    [min_weight_ix] = id;
	weights[min_weight_ix] = weight;
	//printf("bone %d weight %f index %i\n", id, weight, min_weight_ix);
}
void vertex_bone_data_t::normalize() { // make sure all weights sum to 1.0
	float w_sum(0.0);
	for (unsigned i = 0; i < MAX_NUM_BONES_PER_VERTEX; ++i) {w_sum += weights[i];}
	assert(w_sum > 0.0); // or just skip in this case?
	for (unsigned i = 0; i < MAX_NUM_BONES_PER_VERTEX; ++i) {weights[i] /= w_sum;}
}

class model_error_logger_t {
	set<string> errors_printed;
public:
	void log_error(string const &error_str) {
		if (errors_printed.insert(error_str).second) {cerr << format_red(error_str) << endl;} // print if unique
	}
};
model_error_logger_t model_error_logger;

unsigned model3d::get_anim_id(shader_t &shader, string const &prop_name, int anim_id) const {
	assert(has_animations());
	if (anim_id >= 0) return anim_id; // already specified
	string const anim_name(shader.get_property(prop_name));

	if (anim_name.empty()) { // no named animation, use the first one; could also use "default" for the name as that matches the default name
		model_error_logger.log_error("Warning: No animation name or ID specified for model '" + filename + "'; Using the first animation");
		return 0; // use first animation
	}
	anim_id = model_anim_data.get_animation_id_by_name(anim_name);

	if (anim_id < 0) {
		model_error_logger.log_error("Error: Animation '" + anim_name + "' not found in model '" + filename + "'; Using the first animation");
		return 0;
	}
	return anim_id;
}
void model3d::setup_bone_transforms_cached(bone_transform_data_t &cached, shader_t &shader, float anim_time, int anim_id) {
	if (!has_animations()) return;

	if (cached.anim_id == anim_id && cached.anim_time == anim_time) {
		assert(cached.transforms.size() == model_anim_data.bone_transforms.size());
		add_bone_transforms_to_shader(shader, cached.transforms);
	}
	else {
		setup_bone_transforms(shader, anim_time, anim_id);
		unsigned const orig_sz(model_anim_data.bone_transforms.size());
		cached.transforms.swap(model_anim_data.bone_transforms);
		cached.anim_id   = anim_id;
		cached.anim_time = anim_time;
		model_anim_data.bone_transforms.resize(orig_sz); // must be the correct size
	}
}
void model3d::setup_bone_transforms(shader_t &shader, float anim_time, int anim_id) {
	if (!has_animations()) return;
	//highres_timer_t timer("Setup Bone Transforms"); // 0.021ms
	model_anim_data.get_bone_transforms(get_anim_id(shader, "animation_name", anim_id), anim_time);
	add_bone_transforms_to_shader(shader);
}
// do we need a cached version of this one as well?
void model3d::setup_bone_transforms_blended(shader_t &shader, float anim_time1, float anim_time2, float blend_factor, int anim_id1, int anim_id2) {
	if (!has_animations()) return;
	//highres_timer_t timer("Setup Bone Transforms Blend");
	unsigned const id1(get_anim_id(shader, "animation_name" , anim_id1));
	unsigned const id2(get_anim_id(shader, "animation_name2", anim_id2));
	
	if (id1 == id2) { // same animations (maybe there's only one loaded?)
		model_anim_data.get_bone_transforms(id1, anim_time1); // unclear which time to use, so arbitrarily choose time1
	}
	else { // two different animations
		//model_anim_data.blend_animations(id1, id2, blend_factor, delta_time, anim_time1, anim_time2); // requires passing in delta_time (which we don't have), and modifying anim_time1/2
		model_anim_data.blend_animations_simple(id1, id2, blend_factor, anim_time1, anim_time2);
	}
	add_bone_transforms_to_shader(shader);
}
bool model3d::check_anim_wrapped(unsigned anim_id, float old_time, float new_time) const {
	assert(has_animations());
	return model_anim_data.check_anim_wrapped(anim_id, old_time, new_time);
}
float model3d::get_anim_duration(unsigned anim_id) const { // in seconds
	assert(has_animations());
	return model_anim_data.get_anim_duration(anim_id);
}
void model3d::add_bone_transforms_to_shader(shader_t &shader, vector<xform_matrix> const &bone_transforms) const {
	unsigned const MAX_MODEL_BONES = 200; // must agree with shader code
	unsigned num_bones(bone_transforms.size());
	assert(num_bones > 0);

	if (num_bones > MAX_MODEL_BONES) {
		cerr << "Error: Too many bones for model: " << num_bones << "; Max is " << MAX_MODEL_BONES << endl;
		assert(0); // or return? or ignore some bones?
		num_bones = MAX_MODEL_BONES;
	}
	if (!shader.add_uniform_matrix_4x4("bones", bone_transforms[0].get_ptr(), 0, num_bones)) { // transpose=0
		//assert(0); // too strong, as this can trigger when debugging and when using a shader that was setup before the model was loaded
	}
}

template<typename T> void indexed_vntc_vect_t<T>::setup_bones(shader_t &shader, bool is_shadow_pass) {

	if (this->vaos[is_shadow_pass].vao) return; // already set
	this->vaos[is_shadow_pass].ensure_vao_bound();
	vector<T> const &data(*this);
	size_t const vert_mem(data.size()*sizeof(T));

	if (this->vbo) {indexed_vbo_manager_t::pre_render(!indices.empty());}
	else {
		this->vbo = create_vbo();
		check_bind_vbo(this->vbo);
		unsigned const bone_mem(bone_data.vertex_to_bones.size()*sizeof(vertex_bone_data_t)), tot_mem(vert_mem + bone_mem);
		upload_vbo_data(nullptr, tot_mem); // allocate space
		upload_vbo_sub_data(data.data(), 0, vert_mem); // vertex data
		upload_vbo_sub_data(bone_data.vertex_to_bones.data(), vert_mem, bone_mem); // bone data
		this->gpu_mem += tot_mem;

		if (!this->ivbo && !indices.empty()) {
			create_vbo_and_upload(this->ivbo, indices, 1); // is_index=1
			this->gpu_mem += indices.size()*sizeof(index_type_t);
		}
	}
	T::set_vbo_arrays();
	unsigned const stride(sizeof(vertex_bone_data_t));
	glEnableVertexAttribArray(BONE_IDS_LOC);
	glEnableVertexAttribArray(BONE_WEIGHTS_LOC);
	glVertexAttribIPointer(BONE_IDS_LOC,     MAX_NUM_BONES_PER_VERTEX, GL_UNSIGNED_INT, stride, (void *)vert_mem); // bone_ids
	glVertexAttribPointer (BONE_WEIGHTS_LOC, MAX_NUM_BONES_PER_VERTEX, GL_FLOAT, 0, stride, (void *)(vert_mem + MAX_NUM_BONES_PER_VERTEX*sizeof(unsigned))); // bone_weights
}
template<typename T> void indexed_vntc_vect_t<T>::unset_bone_attrs() {
	glDisableVertexAttribArray(BONE_IDS_LOC);
	glDisableVertexAttribArray(BONE_WEIGHTS_LOC);
}

// Note: non-const due to VBO caching
template<typename T> void indexed_vntc_vect_t<T>::render(shader_t &shader, bool is_shadow_pass, point const *const xlate, unsigned npts, bool no_vfc) {

	if (empty()) return;
	assert(npts == 3 || npts == 4);
	//if (is_shadow_pass && this->vbo == 0 && world_mode == WMODE_GROUND) return; // don't create the vbo on the shadow pass (voxel terrain problems - works now?)
	no_vfc |= has_bones(); // disable VFC if animated because animations may move verts outside of the original bcube (should we use a conservative bcube?)

	if (no_vfc) {
		// do nothing
	}
	else if (is_shadow_pass) { // Note: makes shadow map caching more difficult
		if (no_sparse_smap_update() && !orig_camera_pdu.projected_cube_visible(bcube, camera_pdu.pos)) return; // light_pos == camera_pdu.pos for the shadow pass
	}
	else if (this->vbo) { // don't cull if vbo hasn't yet been allocated because this will cause it to be skipped in the shadow pass
		if (!camera_pdu.sphere_and_cube_visible_test(bsphere.pos, bsphere.radius, bcube)) return; // view frustum culling
		
		// Note: null xlate implies there are transforms other than translate, so skip occlusion culling
		if (world_mode == WMODE_GROUND && indices.size() >= 100 && xlate != nullptr && (display_mode & 0x08) != 0) {
			if (cube_cobj_occluded((camera_pdu.pos + *xlate), (bcube + *xlate))) return; // occlusion culling
		}
	}
	assert(!indices.empty()); // now always using indexed drawing
	int prim_type(GL_TRIANGLES);
	unsigned ixn(1), ixd(1), end_ix(indices.size());

	if (!is_shadow_pass && !lod_blocks.empty()) { // block LOD
		float const dmin(2.0*bsphere.radius), dist(p2p_dist(camera_pdu.pos, bsphere.pos));

		if (dist > dmin) { // no LOD if within the bounding sphere
			float const area_thresh((dist - dmin)*(dist - dmin)/(1.0E5f*model_mat_lod_thresh));
			if      (area_thresh > amax) {return;} // draw none
			else if (area_thresh > amin) {end_ix = lod_blocks[get_block_ix(area_thresh)].get_end_ix();}
			assert(end_ix <= indices.size());
		}
	}
	if (npts == 4 && prev_ucc != use_core_context) { // need to rebuild VBOs on core context mode change
		this->clear_vbos();
		prev_ucc = use_core_context;
	}
	assert(!has_bones() || npts == 3); // bones only supported for triangle data

	if (use_core_context && npts == 4) {
		if (!this->ivbo || !this->is_vao_setup(is_shadow_pass)) { // have to setup IVBO once (okay to redo for shadow pass), and VAO for both passes
			vector<unsigned> tixs;
			convert_quad_ixs_to_tri_ixs(indices, tixs);
			this->create_and_upload(*this, tixs, is_shadow_pass, 0, 1); // dynamic_level=0, setup_pointers=1
		}
		ixn = 6; ixd = 4; // convert quads to 2 triangles
	}
	else {
		if (npts == 4) {prim_type = GL_QUADS;}
		if (has_bones()) {setup_bones(shader, is_shadow_pass);}
		else {this->create_and_upload(*this, indices, is_shadow_pass, 0, 1);} // dynamic_level=0, setup_pointers=1
	}
	this->pre_render(is_shadow_pass);
	check_mvm_update();
	
	if (is_shadow_pass || blocks.empty() || no_vfc || camera_pdu.sphere_completely_visible_test(bsphere.pos, bsphere.radius)) { // draw the entire range
		glDrawRangeElements(prim_type, 0, (unsigned)size(), (unsigned)(ixn*end_ix/ixd), GL_UNSIGNED_INT, 0);
		++num_frame_draw_calls;
	}
	else { // draw each block independently
		// could use glDrawElementsIndirect(), but the draw calls don't seem to add any significant overhead for the current set of models
		for (auto i = blocks.begin(); i != blocks.end(); ++i) {
			if (camera_pdu.cube_visible(i->bcube)) {
				glDrawRangeElements(prim_type, 0, (unsigned)size(), (ixn*i->num/ixd), GL_UNSIGNED_INT, (void *)((ixn*i->start_ix/ixd)*sizeof(unsigned)));
				++num_frame_draw_calls;
			}
		}
	}
	this->post_render();
	T::unset_attrs();
	if (has_bones()) {unset_bone_attrs();}
}


template<typename T> void indexed_vntc_vect_t<T>::reserve_for_num_verts(unsigned num_verts) {
	if (empty()) {indices.reserve(num_verts);}
}

template<typename T> void indexed_vntc_vect_t<T>::add_poly(polygon_t const &poly, vertex_map_t<T> &vmap) {
	for (unsigned i = 0; i < poly.size(); ++i) {add_vertex(poly[i], vmap);}
}

template<typename T> void indexed_vntc_vect_t<T>::add_triangle(triangle const &t, vertex_map_t<T> &vmap) {

	vector3d const normal(t.get_normal());
	//vector3d const normal(t.get_normal(!vmap.get_average_normals())); // weight by triangle area
	UNROLL_3X(add_vertex(T(t.pts[i_], normal), vmap);)
}


template<typename T> unsigned indexed_vntc_vect_t<T>::add_vertex(T const &v, vertex_map_t<T> &vmap) {

	T v2(v);
	if (vmap.get_average_normals()) {v2.n = zero_vector;}
	auto it(vmap.find(v2));
	unsigned ix;

	if (it == vmap.end()) { // not found
		ix = (unsigned)size();
		this->push_back(v);
		vmap[v2] = ix;
	}
	else { // found
		ix = it->second;
		assert(ix < size());

		if (vmap.get_average_normals()) {
			operator[](ix).n += v.n; // sum the normals
			need_normalize = 1;
		}
	}
	indices.push_back(ix);
	return ix;
}


template<typename T> float indexed_vntc_vect_t<T>::get_prim_area(unsigned i, unsigned npts) const {

	assert(i+npts <= num_verts());
	float area(triangle_area(get_vert(i).v, get_vert(i+1).v, get_vert(i+2).v)); // first triangle
	if (npts == 4) {area += triangle_area(get_vert(i).v, get_vert(i+2).v, get_vert(i+3).v);} // second triangle (for quads)
	return area;
}

template<typename T> float indexed_vntc_vect_t<T>::calc_area(unsigned npts) {
	
	assert(npts > 0);
	unsigned const nv(num_verts());
	if (nv == 0) return 0.0;
	float area(0.0);
	for (unsigned i = 0; i < nv; i += npts) {area += get_prim_area(i, npts);}
	avg_area_per_tri = area/float(nv/npts);
	return area;
}


struct shared_vertex_t {
	unsigned ai=0, bi=0;
	bool shared=0;
	shared_vertex_t() {}
	shared_vertex_t(unsigned ai_, unsigned bi_) : ai(ai_), bi(bi_), shared(1) {}
};

template<typename T> void indexed_vntc_vect_t<T>::get_polygons(get_polygon_args_t &args, unsigned npts) const {

	if (args.lod_level > 1 && !indices.empty()) {
		indexed_vntc_vect_t<T> simplified_this;
		simplified_this.insert(simplified_this.begin(), begin(), end()); // copy only vertex data; indices will be filled in below, and other fields are unused
		//simplify(simplified_this.indices, 1.0/args.lod_level);
		simplify_meshoptimizer(simplified_this.indices, 1.0/args.lod_level);
		get_polygon_args_t args2(args);
		args2.lod_level = 0;
		simplified_this.get_polygons(args2, npts);
		return;
	}
	unsigned const nv(num_verts());
	if (nv == 0) return;
	assert(npts ==3 || npts == 4);
	assert((nv % npts) == 0);
	polygon_t poly, quad_poly;
	poly.resize(npts);
	quad_poly.resize(4);

	for (unsigned i = 0; i < nv; i += npts) {
		if (npts == 3 && (i+npts) < nv) { // attempt to merge two adjacent triangles into quads
			shared_vertex_t shared1, shared2;

			for (unsigned a = 0; a < 3 && !shared2.shared; ++a) {
				for (unsigned b = 0; b < 3; ++b) {
					if (get_vert(i+a) == get_vert(i+b+3)) {
						(shared1.shared ? shared2 : shared1) = shared_vertex_t(a, b);
						break;
					}
				}
			}
			if (shared2.shared) { // merge two triangles into a single quad
				unsigned nsa(3), nsb(3); // non-sshared

				for (unsigned j = 0; j < 3; ++j) {
					if (shared1.ai != j && shared2.ai != j) nsa = j;
					if (shared1.bi != j && shared2.bi != j) nsb = j;
				}
				assert(nsa < 3 && nsb < 3);
				quad_poly[0] = get_vert(i+shared1.ai);
				quad_poly[1] = get_vert(i+nsa);
				quad_poly[2] = get_vert(i+shared2.ai);
				quad_poly[3] = get_vert(i+nsb+3);

				if (quad_poly.is_coplanar(POLY_COPLANAR_THRESH) && quad_poly.is_convex()) {
					args.polygons.push_back(quad_poly);
					i += npts;
					continue;
				}
			}
		}
		for (unsigned p = 0; p < npts; ++p) {poly[p] = get_vert(i+p);}

		if (args.quads_only) {
			if (npts != 4) {cerr << "Error: model3d mode only works on quads, not triangles" << endl;}
			assert(npts == 4);
			args.polygons.push_back(poly);
		}
		else {
			split_polygon(poly, args.polygons, POLY_COPLANAR_THRESH);
		}
	}
}

template<typename T> void invert_vert_tcy(T &vert) {vert.t[1] = 1.0 - vert.t[1];}
void invert_vert_tcy(vert_norm &v) {} // do nothing (no tcs)

template<typename T> void indexed_vntc_vect_t<T>::invert_tcy() {
	for (auto i = begin(); i != end(); ++i) {invert_vert_tcy(*i);}
}

template<typename T> void indexed_vntc_vect_t<T>::write(ostream &out, bool inc_animations) const {
	vntc_vect_t<T>::write(out);
	write_vector(out, indices);
	if (inc_animations) {write_vector(out, bone_data.vertex_to_bones);}
}
template<typename T> void indexed_vntc_vect_t<T>::read(istream &in, unsigned npts, bool inc_animations) {
	vntc_vect_t<T>::read(in);
	read_vector(in, indices);
	finalize_lod_blocks(npts);
	if (inc_animations) {read_vector(in, bone_data.vertex_to_bones);}
}

// Note: will also match vert_norm_tc_tan, but we don't write the tangent
void write_vertex_to_obj_file(vert_norm_tc const &v, ostream &out) {
	out << "v "  << v.v.raw_str() << endl; // vertex
	out << "vn " << v.n.raw_str() << endl; // normal
	out << "vt " << v.t[0] << " " << v.t[1] << endl; // tex coords (always? do we skip these when they're invalid?)
}

template<typename T> void indexed_vntc_vect_t<T>::write_to_obj_file(ostream &out, unsigned &cur_vert_ix, unsigned npts) const {
	unsigned const nv(num_verts());
	assert((nv % npts) == 0);
	unsigned const start_vert_ix(cur_vert_ix);
	for (auto v = begin(); v != end(); ++v) {write_vertex_to_obj_file(*v, out);}

	for (unsigned i = 0; i < indices.size(); i += npts) {
		out << "f";
		
		for (unsigned n = 0; n < npts; ++n) {
			unsigned const ix(start_vert_ix + indices[i + n] + 1); // always the same index for v/vn/vt; starts at 1
			out << " " << ix << "/" << ix << "/" << ix;
		}
		out << endl;
	} // for i
	cur_vert_ix += size();
}


// ************ polygon_t ************

void polygon_t::from_triangle(triangle const &t) {
	resize(3);
	float const tc[2] = {0.0, 0.0}; // all zero?
	vector3d const normal(t.get_normal());
	UNROLL_3X(operator[](i_) = vert_norm_tc(t.pts[i_], normal, tc);)
}

bool polygon_t::is_convex() const {
	unsigned const npts((unsigned)size());
	assert(npts >= 3);
	if (npts == 3) return 1;
	unsigned counts[2] = {0};
	vector3d const norm(get_planar_normal());

	for (unsigned i = 0; i < npts; ++i) {
		unsigned const ip((i+npts-1)%npts), in((i+1)%npts);
		++counts[dot_product(norm, cross_product((*this)[i].v-(*this)[ip].v, (*this)[in].v-(*this)[i].v)) < 0.0];
	}
	return !(counts[0] && counts[1]);
}

bool polygon_t::is_coplanar(float thresh) const {
	assert(size() >= 3);
	if (size() == 3 || thresh == 0.0) return 1;
	vector3d n2;
	get_normal((*this)[0].v, (*this)[2].v, (*this)[3].v, n2, 1);
	return (dot_product(get_planar_normal(), n2) > thresh);
}

vector3d polygon_t::get_planar_normal() const {
	assert(size() >= 3);
	vector3d norm;
	get_normal((*this)[0].v, (*this)[1].v, (*this)[2].v, norm, 1);
	return norm;
}

void polygon_t::from_points(vector<point> const &pts) {
	resize(pts.size());
	for (unsigned i = 0; i < size(); ++i) {(*this)[i].v = pts[i];}
}


// ************ vntc_vect_block_t ************

// explicit template instantiations of all used vntc_vect_block_t types
template struct vntc_vect_block_t<vert_norm>;
template struct vntc_vect_block_t<vert_norm_tc>;
template struct vntc_vect_block_t<vert_norm_tc_tan>;


template<typename T> void vntc_vect_block_t<T>::finalize(unsigned npts) {
	for (auto i = begin(); i != end(); ++i) {i->finalize(npts);}
}

template<typename T> void vntc_vect_block_t<T>::free_vbos() {
	for (auto i = begin(); i != end(); ++i) {i->clear_vbos();}
}

template<typename T> cube_t vntc_vect_block_t<T>::get_bcube() const {
	if (this->empty()) return all_zeros_cube;
	cube_t bcube(this->front().get_bcube());
	for (auto i = begin()+1; i != end(); ++i) {bcube.union_with_cube(i->get_bcube());}
	return bcube;
}

template<typename T> size_t vntc_vect_block_t<T>::get_gpu_mem() const {
	size_t mem(0);
	for (auto i = begin(); i != end(); ++i) {mem += i->get_gpu_mem();}
	return mem;
}

template<typename T> unsigned vntc_vect_block_t<T>::num_verts() const {
	unsigned s(0);
	for (auto i = begin(); i != end(); ++i) {s += i->num_verts();}
	return s;
}

template<typename T> unsigned vntc_vect_block_t<T>::num_unique_verts() const {
	unsigned s(0);
	for (auto i = begin(); i != end(); ++i) {s += (unsigned)i->size();}
	return s;
}

template<typename T> float vntc_vect_block_t<T>::calc_draw_order_score() const {

	float area(0.0);
	unsigned count(0);

	for (auto i = begin(); i != end(); ++i) {
		count += i->size();
		area  += i->get_bradius()*i->get_bradius();
	}
	return ((count == 0) ? 0.0 : -area/count);
}

template<typename T> float vntc_vect_block_t<T>::calc_area(unsigned npts) {
	float area(0.0);
	for (auto i = begin(); i != end(); ++i) {area += i->calc_area(npts);}
	return area;
}

template<typename T> void vntc_vect_block_t<T>::get_polygons(get_polygon_args_t &args, unsigned npts) const {
	for (auto i = begin(); i != end(); ++i) {i->get_polygons(args, npts);}
}

template<typename T> void vntc_vect_block_t<T>::invert_tcy() {
	for (auto i = begin(); i != end(); ++i) {i->invert_tcy();}
}

template<typename T> void vntc_vect_block_t<T>::simplify_indices(float reduce_target) {
	for (auto i = begin(); i != end(); ++i) {i->simplify_indices(reduce_target);}
}

template<typename T> void vntc_vect_block_t<T>::reverse_winding_order(unsigned npts) {
	for (auto i = begin(); i != end(); ++i) {i->reverse_winding_order(npts);}
}

template<typename T> void vntc_vect_block_t<T>::merge_into_single_vector() {
	if (this->size() <= 1) return; // nothing to merge
	unsigned tot_verts(0), tot_ixs(0);

	for (auto i = begin(); i != end(); ++i) {
		for (auto j = i->indices.begin(); j != i->indices.end(); ++j) {*j += tot_verts;} // offset indices by current vertex offset
		tot_verts += i->size();
		tot_ixs   += i->indices.size();
	}
	auto &dest(this->front());
	dest.reserve(tot_verts);
	dest.indices.reserve(tot_ixs);

	for (auto i = begin()+1; i != end(); ++i) { // skip first block
		vector_add_to(*i, dest); // merge verts
		vector_add_to(i->indices, dest.indices); // merge indices
	}
	dest.calc_bounding_volumes(); // can be optimized
	dest.clear_blocks(); // no longer valid
	//dest.finalize_lod_blocks(npts); // is this needed? it doesn't really make sense to merge blocks and then re-split, so I guess not
	this->resize(1); // remove all but the first block
}

template<typename T> bool vntc_vect_block_t<T>::write(ostream &out, bool inc_animations) const {
	write_uint(out, (unsigned)this->size());
	for (auto i = begin(); i != end(); ++i) {i->write(out, inc_animations);}
	return out.good();
}
template<typename T> bool vntc_vect_block_t<T>::read(istream &in, unsigned npts, bool inc_animations) {
	this->clear();
	this->resize(read_uint(in));
	for (auto i = begin(); i != end(); ++i) {i->read(in, npts, inc_animations);}
	if (merge_model_objects) {merge_into_single_vector();} // model was split per object, and we don't want that; merge into a single vector
	return in.good();
}

template<typename T> bool vntc_vect_block_t<T>::write_to_obj_file(ostream &out, unsigned &cur_vert_ix, unsigned npts) const {
	for (auto i = begin(); i != end(); ++i) {i->write_to_obj_file(out, cur_vert_ix, npts);}
	return out.good();
}


// ************ geometry_t ************


template<> void geometry_t<vert_norm_tc_tan>::calc_tangents_blocks(vntc_vect_block_t<vert_norm_tc_tan> &blocks, unsigned npts) {
	for (auto i = blocks.begin(); i != blocks.end(); ++i) {i->calc_tangents(npts);}
}

template<typename T> void geometry_t<T>::calc_tangents() {
	calc_tangents_blocks(triangles, 3);
	calc_tangents_blocks(quads,     4);
}

template<typename T> void geometry_t<T>::render_blocks(shader_t &shader, bool is_shadow_pass, point const *const xlate, vntc_vect_block_t<T> &blocks, unsigned npts) {
	for (auto i = blocks.begin(); i != blocks.end(); ++i) {i->render(shader, is_shadow_pass, xlate, npts);}
}

template<typename T> void geometry_t<T>::render(shader_t &shader, bool is_shadow_pass, point const *const xlate) {
	render_blocks(shader, is_shadow_pass, xlate, triangles, 3);
	render_blocks(shader, is_shadow_pass, xlate, quads,     4);
}

template<typename T> void geometry_t<T>::add_poly_to_polys(polygon_t const &poly, vntc_vect_block_t<T> &v, vertex_map_t<T> &vmap, unsigned obj_id) const {

	if (v.empty() || v.back().size() > MAX_VMAP_SIZE || (!merge_model_objects && obj_id > v.back().obj_id)) {
		vmap.clear();
		if (!merge_model_objects || v.empty()) {v.push_back(indexed_vntc_vect_t<T>(obj_id));}
	}
	v.back().add_poly(poly, vmap);
}

template<typename T> void geometry_t<T>::add_poly(polygon_t const &poly, vertex_map_t<T> vmap[2], unsigned obj_id) {
	
	if (poly.size() == 3) { // triangle
		add_poly_to_polys(poly, triangles, vmap[0], obj_id);
	}
	else if (poly.size() == 4) {
		add_poly_to_polys(poly, quads, vmap[1], obj_id);
	}
	else {
		assert(0); // shouldn't get here
	}
}

template<typename T> void geometry_t<T>::get_polygons(get_polygon_args_t &args) const {
	triangles.get_polygons(args, 3); // should be empty in quads_only mode (will be checked)
	quads    .get_polygons(args, 4);
}

template<typename T> cube_t geometry_t<T>::get_bcube() const {

	cube_t bcube; // will return this if empty

	if (!triangles.empty()) {
		bcube = triangles.get_bcube();
		if (!quads.empty()) {bcube.union_with_cube(quads.get_bcube());}
	}
	else if (!quads.empty()) {
		bcube = quads.get_bcube();
	}
	return bcube;
}

template<typename T> void geometry_t<T>::clear() {

	free_vbos();
	triangles.clear();
	quads    .clear();
}

template<typename T> void geometry_t<T>::get_stats(model3d_stats_t &stats) const {
	
	stats.tris  += triangles.num_verts()/3;
	stats.quads += quads    .num_verts()/4;
	triangles.get_stats(stats);
	quads    .get_stats(stats);
}

template<typename T> void geometry_t<T>::calc_area(float &area, unsigned &ntris) {
	area  += triangles.calc_area(3) + quads.calc_area(4);
	ntris += (triangles.num_verts()/3) + 2*(quads.num_verts()/4); // quads count as 2 triangles
}

template<typename T> void geometry_t<T>::simplify_indices(float reduce_target) {
	triangles.simplify_indices(reduce_target); // mesh simplification only applies to triangles, not quads
}

template<typename T> void geometry_t<T>::reverse_winding_order() {
	triangles.reverse_winding_order(3);
	quads    .reverse_winding_order(4);
}


// ************ material_t ************


template<typename T> void update_score(vntc_vect_block_t<T> const &v, float &score, unsigned &num_nonempty) {
	if (v.empty()) return;
	score += v.calc_draw_order_score();
	++num_nonempty;
}


void material_t::compute_area_per_tri() {

	if (avg_area_per_tri > 0) return; // already computed
	unsigned tris(0);
	tot_tri_area = 0;
	geom    .calc_area(tot_tri_area, tris);
	geom_tan.calc_area(tot_tri_area, tris);
	if (tris > 0) {avg_area_per_tri = alpha*tot_tri_area/tris;}
	//cout << "name: " << name << " " << TXT(tris) << TXT(tot_tri_area) << TXT(alpha) << "value: " << (1.0E6*avg_area_per_tri) << endl;
}

void material_t::simplify_indices(float reduce_target) {
	geom    .simplify_indices(reduce_target);
	geom_tan.simplify_indices(reduce_target);
}

void material_t::reverse_winding_order() {
	geom    .reverse_winding_order();
	geom_tan.reverse_winding_order();
}

void material_t::ensure_textures_loaded(texture_manager &tmgr) {
	tmgr.ensure_tid_loaded(get_render_texture(), 0); // only one tid for now
	// if bump_tid is set, but bump maps are disabled, then clear bump_tid because either a) we won't use it, or b) it won't be loaded later when we try to use it
	if (use_bump_map()) {tmgr.ensure_tid_loaded(bump_tid, 1);} else {bump_tid = -1;}
	if (use_spec_map()) {tmgr.ensure_tid_loaded( s_tid,   0);} else {s_tid    = -1;}
	if (use_spec_map()) {tmgr.ensure_tid_loaded(ns_tid,   0);} else {ns_tid   = -1;}
}

void maybe_upload_and_free(texture_manager &tmgr, unsigned tid) {
	if (tid < BUILTIN_TID_START) {tmgr.ensure_tid_bound(tid);} // upload to GPU and free if not a built-in texture
}

void material_t::init_textures(texture_manager &tmgr) {

	if (!mat_is_used()) return;
	int const tid(get_render_texture());
	tmgr.bind_alpha_channel_to_texture(tid, alpha_tid);
	ensure_textures_loaded(tmgr);
	might_have_alpha_comp |= tmgr.might_have_alpha_comp(tid);
	
	if (no_store_model_textures_in_memory) { // now that textures have been loaded, send the data to the GPU so that they can be freed early
		maybe_upload_and_free(tmgr, get_render_texture());
		if (use_bump_map()) {maybe_upload_and_free(tmgr, bump_tid);}
		if (use_spec_map()) {maybe_upload_and_free(tmgr, s_tid);}
		if (use_spec_map()) {maybe_upload_and_free(tmgr, ns_tid);}
	}
}

void material_t::queue_textures_to_load(texture_manager &tmgr) {

	if (!mat_is_used()) return;
	int const tid(get_render_texture());
	tmgr.bind_alpha_channel_to_texture(tid, alpha_tid);
	tmgr.add_work_item(tid, 0);
	if (use_bump_map()) {tmgr.add_work_item(bump_tid, 1);} else {bump_tid = -1;}
	if (use_spec_map()) {tmgr.add_work_item( s_tid,   0);} else {s_tid    = -1;}
	if (use_spec_map()) {tmgr.add_work_item(ns_tid,   0);} else {ns_tid   = -1;}
}

void material_t::check_for_tc_invert_y(texture_manager &tmgr) {

	if (tcs_checked) return; // already done
	int const tid(get_render_texture());
	if (tid < 0) return; // no texture
	texture_t &texture(tmgr.get_texture(tid));

	if (texture.is_inverted_y_type() && !texture.invert_y) { // compressed DDS texture, need to invert tex coord in Y
		geom.invert_tcy();
		geom_tan.invert_tcy();
		texture.invert_y ^= 1; // already inverted, don't try to invert again (FIXME: doesn't work if used in multiple materials)
	}
	tcs_checked = 1;
}


void texture_manager::bind_texture_tu_or_white_tex(int tid, unsigned tu_id) const {
	bind_texture_tu_def_white_tex(((tid >= 0) ? get_texture(tid).get_tid() : 0), tu_id);
}

// enable_alpha_mask: 0=non-alpha mask only, 1=alpha mask only, 2=both
void material_t::render(shader_t &shader, texture_manager const &tmgr, int default_tid, bool is_shadow_pass,
	bool is_z_prepass, int enable_alpha_mask, bool is_bmap_pass, point const *const xlate, bool no_set_min_alpha)
{
	if (empty() || skip || alpha == 0.0) return; // empty or transparent
	if (is_shadow_pass && alpha < MIN_SHADOW_ALPHA) return;

	if (draw_order_score == 0) {
		unsigned num_nonempty(0);
		update_score(geom.triangles,     draw_order_score, num_nonempty);
		update_score(geom.quads,         draw_order_score, num_nonempty);
		update_score(geom_tan.triangles, draw_order_score, num_nonempty);
		update_score(geom_tan.quads,     draw_order_score, num_nonempty);
		assert(num_nonempty > 0);
		draw_order_score /= num_nonempty; // take the average
	}
	int const tex_id(get_render_texture());

	if (is_z_prepass) { // no textures
		if (alpha < 1.0 || has_alpha_mask()) return; // partially transparent or has alpha mask
		geom.render(shader, 0, xlate);
		geom_tan.render(shader, 0, xlate);
	}
	else if (is_shadow_pass) {
		bool const has_amask(has_alpha_mask());
		if (enable_alpha_mask != 2 && has_amask != bool(enable_alpha_mask)) return; // incorrect pass
		if (has_amask) {tmgr.bind_texture(tex_id);} // enable alpha mask texture
		geom.render(shader, 1, xlate);
		geom_tan.render(shader, 1, xlate);
		if (has_amask) {select_texture(WHITE_TEX);} // back to a default white texture
	}
	else {
		int bump_map_mag_ix(-1), metalness_ix(-1), refract_ix_ix(-1);
		
		if (disable_model_textures) {
			select_texture(WHITE_TEX);
		}
		else if (tex_id >= 0) {
			tmgr.bind_texture(tex_id);
		}
		else if (tex_id != -2) { // Note: the special tid of -2 is used to indicate the caller has already bound the correct texture, so we should not bind a texture here
			select_texture((default_tid >= 0) ? default_tid : WHITE_TEX); // no texture specified - use white texture
		}
		if (use_bump_map()) {
			bind_texture_tu(tmgr.get_texture(bump_tid).get_tid(), 5);
		}
		else if (is_bmap_pass) {
			if (enable_bump_map()) {bind_default_flat_normal_map();} // use default normal map in this case instead of leaving it unbound, or bound to the previous material
			// disable bump map, including TBN matrix transform of eye_pos and light_dir; needed when there are no TCs or tangent vectors
			bump_map_mag_ix = shader.get_uniform_loc("bump_map_mag");
			if (bump_map_mag_ix >= 0) {shader.set_uniform_float(bump_map_mag_ix, 0.0);}
		}
		if (enable_spec_map()) { // all white/specular if no specular map texture
			tmgr.bind_texture_tu_or_white_tex(s_tid,  8); // specular map
			tmgr.bind_texture_tu_or_white_tex(ns_tid, 9); // gloss map (Note: unclear how to interpret map_ns in object files)
		}
		if (!disable_shader_effects && metalness >= 0.0) {metalness_ix = shader.get_uniform_loc("metalness");}
		if (metalness_ix >= 0) {shader.set_uniform_float(metalness_ix, metalness);} // set metalness if specified/valid; may or may not be used
		// set index of refraction - may not actually be used (ENABLE_CUBE_MAP_REFLECT/ENABLE_REFLECTIONS)
		if (!disable_shader_effects /*&& alpha < 1.0*/ && ni != 1.0 && ni > 0.0) {refract_ix_ix = shader.get_uniform_loc("refract_ix");}
		if (refract_ix_ix >= 0) {shader.set_uniform_float(refract_ix_ix, ni);}
		bool const need_blend(is_partial_transparent()); // conservative, but should be okay
		if (need_blend) {enable_blend(); /*glDepthMask(GL_FALSE);*/ /*glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE);*/}

		if (!no_set_min_alpha && !shader.get_user_flag(SHADER_FLAG_NO_ALPHA_TEST)) { // only set min_alpha if alpha test is enabled
			bool const has_binary_alpha(!disable_model_textures && tex_id >= 0 && tmgr.has_binary_alpha(tex_id));
			float const min_alpha(min(0.99*alpha, (get_needs_alpha_test() ? (has_binary_alpha ? 0.9 : model3d_alpha_thresh) : 0.0)));
			shader.add_uniform_float("min_alpha", min_alpha);
		}
		if (ns > 0.0)    {shader.set_specular_color(ks, ns);} // ns<=0 is undefined?
		if (ke != BLACK) {shader.set_color_e(colorRGBA(ke, alpha));}
		// Note: ka is ignored here because it represents a "fake" lighting model;
		// 3DWorld uses a more realistic lighting model where ambient comes from indirect lighting that's computed independently from the material;
		// however, it might make sense to use ka instead of ke when ke is not specified?
		shader.set_cur_color(get_ad_color());
		geom    .render(shader, 0, xlate);
		geom_tan.render(shader, 0, xlate);
		shader.clear_color_e(); // Note: caller can cheat and override the emissive of the first material because it won't be set above but will be cleared here
		if (ns > 0.0)    {shader.clear_specular();}
		if (need_blend)  {disable_blend(); /*glDepthMask(GL_TRUE);*/ /*glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);*/}
		if (refract_ix_ix   >= 0) {shader.set_uniform_float(refract_ix_ix,   1.0);}
		if (metalness_ix    >= 0) {shader.set_uniform_float(metalness_ix,    0.0);} // if metalness was specified, reset to the default of 0.0 for the next material
		if (bump_map_mag_ix >= 0) {shader.set_uniform_float(bump_map_mag_ix, 1.0);} // reset back to default
	}
}


bool material_t::use_bump_map() const {return (enable_bump_map() && bump_tid >= 0);}
bool material_t::use_spec_map() const {return (enable_spec_map() && (s_tid >= 0 || ns_tid >= 0));}

colorRGBA material_t::get_ad_color() const {
	colorRGBA c(kd, alpha);
	c.set_valid_color();
	return c;
}

colorRGBA material_t::get_avg_color(texture_manager const &tmgr, int default_tid) const {
	colorRGBA avg_color(get_ad_color());
	int tex_id(get_render_texture());
	if (tex_id >= 0) {return avg_color.modulate_with(tmgr.get_tex_avg_color(tex_id));}
	else if (default_tid >= 0) {return avg_color.modulate_with(texture_color(default_tid));}
	return avg_color;
}

mesh_bone_data_t &material_t::get_bone_data_for_last_added_tri_mesh() {
	if (USE_ANIM_MODEL_TANGENTS && model_calc_tan_vect && use_bump_map()) {
		assert(!geom_tan.triangles.empty());
		return geom_tan.triangles.back().bone_data;
	}
	assert(!geom.triangles.empty());
	return geom.triangles.back().bone_data;
}

template<typename T> unsigned geometry_t<T>::add_triangles(vector<vert_norm_tc> const &verts, vector<unsigned> const &indices, bool add_new_block) {
	if (add_new_block || triangles.empty()) { // new block; can directly copy indices
		triangles.push_back(indexed_vntc_vect_t<T>());
		vector_add_to(verts, triangles.back());
		triangles.back().indices = indices;
		return 0;
	}
	auto &dest(triangles.back());
	assert(indices.empty() == dest.indices.empty()); // can't mix indexed with non-indexed triangles
	unsigned const ixs_off(dest.size());
	vector_add_to(verts, dest);
	for (unsigned ix : indices) {dest.indices.push_back(ix + ixs_off);} // adjust indices based on existing vertices
	return ixs_off;
}

// Note: no quads; indices are optional; returns the index offset of the first vertex
unsigned material_t::add_triangles(vector<vert_norm_tc> const &verts, vector<unsigned> const &indices, bool add_new_block) {
	if (verts.empty()) {assert(indices.empty()); return 0;} // no triangles? error?
	mark_as_used();
	if (USE_ANIM_MODEL_TANGENTS && model_calc_tan_vect && use_bump_map()) {return geom_tan.add_triangles(verts, indices, add_new_block);} // with tangents for normal maps
	return geom.add_triangles(verts, indices, add_new_block); // without tangents
}

bool material_t::add_poly(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], unsigned obj_id) {
	if (skip) return 0;
	if (model_calc_tan_vect && use_bump_map()) {geom_tan.add_poly(poly, vmap_tan, obj_id);}
	else {geom.add_poly(poly, vmap, obj_id);}
	mark_as_used();
	return 1;
}


bool material_t::write(ostream &out, bool inc_animations) const {
	out.write((char const *)this, sizeof(material_params_t));
	write_vector(out, name);
	write_vector(out, filename);
	return (geom.write(out, inc_animations) && geom_tan.write(out, inc_animations));
}
bool material_t::read(istream &in, bool inc_animations) {
	in.read((char *)this, sizeof(material_params_t));
	read_vector(in, name);
	read_vector(in, filename);
	return (geom.read(in, inc_animations) && geom_tan.read(in, inc_animations));
}

bool material_t::write_to_obj_file(ostream &out, unsigned &cur_vert_ix) const {
	out << "usemtl " << name << endl;
	return (geom.write_to_obj_file(out, cur_vert_ix) && geom_tan.write_to_obj_file(out, cur_vert_ix));
}

void write_mtllib_texture(string const &opt_name, int tid, texture_manager const &tmgr, ostream &out) {
	if (tid >= 0) {out << "\t" << opt_name << " " << tmgr.get_texture(tid).name << endl;}
}
void material_t::write_mtllib_entry(ostream &out, texture_manager const &tmgr) const {
	out << "newmtl "  << name  << endl;
	out << "\tNs "    << ns    << endl;
	out << "\tNi "    << ni    << endl;
	out << "\td "     << alpha << endl;
	out << "\tTr "    << tr    << endl;
	out << "\tTf "    << tf.raw_str() << endl;
	out << "\tillum " << illum << endl;
	out << "\tKa "    << ka.raw_str() << endl;
	out << "\tKd "    << kd.raw_str() << endl;
	out << "\tKs "    << ks.raw_str() << endl;
	out << "\tKe "    << ke.raw_str() << endl;
	write_mtllib_texture("map_Ka",   a_tid,     tmgr, out);
	write_mtllib_texture("map_Kd",   d_tid,     tmgr, out);
	write_mtllib_texture("map_Ks",   s_tid,     tmgr, out);
	write_mtllib_texture("map_ns",   ns_tid,    tmgr, out);
	write_mtllib_texture("map_d",    alpha_tid, tmgr, out);
	write_mtllib_texture("map_bump", bump_tid,  tmgr, out);
	write_mtllib_texture("map_refl", refl_tid,  tmgr, out);
	// what about metalness?
}


// ************ model3d ************


material_t &model3d::get_material(int mat_id, bool alloc_if_needed) {
	if (alloc_if_needed && mat_id >= (int)materials.size()) {materials.resize(mat_id+1);} // allocate additional material(s) if needed
	assert(mat_id >= 0);
	assert(mat_id < (int)materials.size());
	return materials[mat_id];
}

void coll_tquads_from_triangles(vector<triangle> const &triangles, vector<coll_tquad> &ppts, colorRGBA const &color) {
	ppts.reserve(ppts.capacity() + triangles.size());
	for (unsigned i = 0; i < triangles.size(); ++i) ppts.push_back(coll_tquad(triangles[i], color));
}

unsigned model3d::add_polygon(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], int mat_id, unsigned obj_id) {
	
	for (unsigned d = 0; d < 2; ++d) {
		vmap[d].check_for_clear(mat_id);
		vmap_tan[d].check_for_clear(mat_id);
	}
	split_polygons_buffer.resize(0);
	split_polygon(poly, split_polygons_buffer, 0.0, allow_model3d_quads);

	for (polygon_t &p : split_polygons_buffer) {
		if (mat_id < 0) {unbound_geom.add_poly(p, vmap, obj_id);}
		else {
			assert((unsigned)mat_id < materials.size());
			materials[mat_id].add_poly(p, vmap, vmap_tan, obj_id);
		}
	}
	if (mat_id < 0 || !materials[mat_id].skip) {update_bbox(poly);} // don't include skipped materials in the bbox
	return (unsigned)split_polygons_buffer.size();
}

void model3d::add_triangle(polygon_t const &tri, vntc_map_t &vmap, int mat_id, unsigned obj_id) {

	assert(tri.size() == 3);
	vmap.check_for_clear(mat_id);

	if (mat_id < 0) {
		unbound_geom.add_poly_to_polys(tri, unbound_geom.triangles, vmap, obj_id);
	}
	else {
		assert((unsigned)mat_id < materials.size());
		materials[mat_id].geom.add_poly_to_polys(tri, materials[mat_id].geom.triangles, vmap, obj_id);
		materials[mat_id].mark_as_used();
	}
	if (mat_id < 0 || !materials[mat_id].skip) {update_bbox(tri);} // don't include skipped materials in the bbox
}


void model3d::update_bbox(polygon_t const &poly) {
	cube_t const bb(get_polygon_bbox(poly));
	if (bcube == all_zeros_cube) {bcube = bb;} else {bcube.union_with_cube(bb);}
}

void model3d::get_transformed_bcubes(vector<cube_t> &bcubes) const {
	if (transforms.empty()) {bcubes.push_back(bcube); return;} // no transforms
	for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {bcubes.push_back(xf->get_xformed_cube(bcube));} // Note: const, non-cached call
}


void model3d::get_polygons(vector<coll_tquad> &polygons, bool quads_only, bool apply_transforms, unsigned lod_level) const {

	unsigned const start_pix(polygons.size());

	if (start_pix == 0) { // Note: we count quads as 1.5 polygons because some of them may be split into triangles
		model3d_stats_t stats;
		get_stats(stats);
		unsigned const num_copies((!apply_transforms || transforms.empty()) ? 1 : transforms.size());
		polygons.reserve(num_copies*(quads_only ? stats.quads : (stats.tris + 1.5*stats.quads)));
	}
	get_polygon_args_t args(polygons, quads_only, lod_level);
	unbound_geom.get_polygons(args);

	for (material_t const &m : materials) {
		m.geom.get_polygons    (args);
		m.geom_tan.get_polygons(args);
	}
	if (apply_transforms && !transforms.empty()) { // handle transforms
		// first clone the polygons for each transform; first transform is done already
		unsigned const num_polys(polygons.size() - start_pix);

		for (unsigned i = 1; i < transforms.size(); ++i) { // N-1 block copies
			for (unsigned p = 0; p < num_polys; ++p) {polygons.push_back(polygons[start_pix + p]);}
		}
		assert(polygons.size() == start_pix + transforms.size()*num_polys); // can be removed later
		unsigned pix(start_pix);

		for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {
			for (unsigned p = 0; p < num_polys; ++p) {xf->apply_to_tquad(polygons[pix++]);}
		}
		assert(pix == polygons.size());
	}
	//::remove_excess_cap(polygons); // slightly slower, but slightly less memory usage
}


void calc_bounds(cube_t const &c, int bounds[2][2], float spacing) {

	for (unsigned d = 0; d < 2; ++d) {
		for (unsigned e = 0; e < 2; ++e) {
			bounds[d][e] = round_fp(c.d[d][e]/spacing);
		}
	}
}


struct float_plus_dir {
	float f=0.0;
	bool d=0;
	float_plus_dir() {}
	float_plus_dir(float f_, bool d_) : f(f_), d(d_) {}
	bool operator<(float_plus_dir const &fd) const {return ((f == fd.f) ? (d < fd.d) : (f < fd.f));}
};

void xform_polygons(vector<coll_tquad> &ppts, model3d_xform_t const &xf, unsigned start_ix=0) {
	if (xf.is_identity()) return;
#pragma omp parallel for schedule(static,1)
	for (int i = start_ix; i < (int)ppts.size(); ++i) {xf.apply_to_tquad(ppts[i]);}
}

template<typename T> unsigned add_polygons_to_voxel_grid(vector<coll_tquad> &polygons, T const &cont,
	vector<vector<float_plus_dir> > &zvals, int bounds[2][2], int num_xy[2], float spacing, unsigned &nhq, model3d_xform_t const &xf)
{
	model3d_stats_t stats;
	cont.get_stats(stats);
	polygons.clear();
	polygons.reserve(stats.quads);
	get_polygon_args_t args(polygons, 1); // quads_only=1
	cont.get_polygons(args);
	xform_polygons(polygons, xf, 0);
	
	for (coll_tquad const &i : polygons) { // Note: colors are unused
		assert(i.npts == 4);
		if (fabs(i.normal.z) < 0.99)    continue; // only keep top/bottom cube sides
		cube_t const bcube(i.get_bcube());
		if ((bcube.dz()) > 0.5*spacing) continue; // can this happen?
		int cbounds[2][2];
		calc_bounds(bcube, cbounds, spacing);
		bool const is_top(i.normal.z > 0.0);
		float_plus_dir const fd(bcube.d[2][0], is_top);
		++nhq;
		
		for (int y = cbounds[1][0]; y < cbounds[1][1]; ++y) {
			for (int x = cbounds[0][0]; x < cbounds[0][1]; ++x) {
				int const xv(x - bounds[0][0]), yv(y - bounds[1][0]);
				assert(xv >= 0 && yv >= 0 && xv < num_xy[0] && yv < num_xy[1]);
				zvals[xv + num_xy[0]*yv].push_back(fd);
			}
		}
	}
	return (unsigned)polygons.size();
}


// Note: ignores model transforms, which is why xf is passed in
void model3d::get_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf) const {

	RESET_TIME;
	float const spacing(xf.voxel_spacing);
	assert(spacing > 0.0);

	// calculate scene voxel bounds
	int bounds[2][2] = {}, num_xy[2] = {}; // {x,y}x{lo,hi}
	calc_bounds(bcube, bounds, spacing);
	cout << "bounds: ";

	for (unsigned d = 0; d < 2; ++d) {
		cout << bounds[d][0] << " " << bounds[d][1] << " ";
		num_xy[d] = (bounds[d][1] - bounds[d][0]);
	}
	unsigned const num_tot(num_xy[0]*num_xy[1]);
	vector<vector<float_plus_dir> > zvals(num_tot);
	unsigned num_horiz_quads(0), num_polys(0), num_pre_merged_cubes(0);
	cout << ", size: " << num_xy[0] << "x" << num_xy[1] << " = " << num_tot << endl;

	// create z ranges for each xy voxel column from polygons
	{
		// we technically only want the horizontal quads, but it's difficult to filter them out earlier
		vector<coll_tquad> polygons;
		num_polys += add_polygons_to_voxel_grid(polygons, unbound_geom, zvals, bounds, num_xy, spacing, num_horiz_quads, xf);

		for (material_t const &m : materials) {
			num_polys += add_polygons_to_voxel_grid(polygons, m.geom,     zvals, bounds, num_xy, spacing, num_horiz_quads, xf);
			num_polys += add_polygons_to_voxel_grid(polygons, m.geom_tan, zvals, bounds, num_xy, spacing, num_horiz_quads, xf);
		}
	}

	// convert voxel columns to cubes
	unsigned const MERGE_LENGTH = 16;
	cubes.reserve(cubes.size() + num_horiz_quads/2);

	for (int y = bounds[1][0]; y < bounds[1][1]; ++y) {
		for (int x = bounds[0][0]; x < bounds[0][1]; ++x) {
			int const xv(x - bounds[0][0]), yv(y - bounds[1][0]);
			assert(xv >= 0 && yv >= 0 && xv < num_xy[0] && yv < num_xy[1]);
			vector<float_plus_dir> &zv(zvals[xv + num_xy[0]*yv]);
			sort(zv.begin(), zv.end());
			cube_t cube(x*spacing, (x+1)*spacing, y*spacing, (y+1)*spacing, 0.0, 0.0);
			//bool in_cube(0);

			for (unsigned j = 0; j < zv.size(); ++j) {
				float const val(zv[j].f);

				if (j+1 < zv.size() && fabs(zv[j+1].f - val) < TOLERANCE) {
					assert(zv[j].d != zv[j+1].d);
					++j;
					continue; // a canceling pair forming a zero length segment - skip both
				}
				if (!zv[j].d) { // bottom
					//assert(!in_cube);
					cube.d[2][0] = val;
				}
				else { // top
					//assert(in_cube);
					assert(j > 0); // bottom must be set
					cube.d[2][1] = val;
					assert(cube.d[2][0] < cube.d[2][1]); // no zero height cubes (should be ok, but can remove later if too strict)
					bool merged(0);
					unsigned num(0);

					for (vector<cube_t>::reverse_iterator i = cubes.rbegin(); i != cubes.rend() && !merged && num < MERGE_LENGTH; ++i, ++num) {
						merged = i->cube_merge(cube);
					}
					if (!merged) {cubes.push_back(cube);}
					cube.d[2][0] = val; // next segment starts here in case we get two top edges in a row
					++num_pre_merged_cubes;
				}
				//in_cube ^= 1;
			}
			//assert(!in_cube);
		}
	}
	PRINT_TIME("Get Cubes");
	cout << "polygons: " << num_polys << ", hquads: " << num_horiz_quads << ", pre_merged_cubes: " << num_pre_merged_cubes << ", cubes: " << cubes.size() << endl;
}


colorRGBA model3d::get_avg_color() const {
	colorRGBA avg_color(BLACK);
	for (auto m = materials.begin(); m != materials.end(); ++m) {avg_color += m->get_avg_color(tmgr, unbound_mat.tid);}
	if (!unbound_geom.empty()) {avg_color += unbound_mat.color.modulate_with(texture_color(unbound_mat.tid));}
	if (avg_color.alpha > 0.0) {avg_color = avg_color/avg_color.alpha; avg_color.alpha = 1.0;} // normalize to alpha of 1.0
	return avg_color;
}
colorRGBA color_mult_inc_alpha(colorRGBA const &c, float v) {return colorRGBA(c.R*v, c.G*v, c.B*v, c.A*v);}

colorRGBA model3d::get_area_weighted_avg_color() {
	colorRGBA avg_color(BLACK);

	for (auto m = materials.begin(); m != materials.end(); ++m) {
		m->compute_area_per_tri(); // Note: will cache, which is why this function isn't const
		avg_color += color_mult_inc_alpha(m->get_avg_color(tmgr, unbound_mat.tid), m->tot_tri_area);
	}
	if (!unbound_geom.empty()) {
		float area(0.0);
		unsigned ntris(0); // unused
		unbound_geom.calc_area(area, ntris);
		avg_color += color_mult_inc_alpha(unbound_mat.color.modulate_with(texture_color(unbound_mat.tid)), area);
	}
	if (avg_color.alpha > 0.0) {avg_color = avg_color/avg_color.alpha; avg_color.alpha = 1.0;} // normalize to alpha of 1.0
	return avg_color;
}
colorRGBA model3d::get_and_cache_avg_color(bool area_weighted) {
	if (cached_avg_color.A == 0.0) {cached_avg_color = (area_weighted ? get_area_weighted_avg_color() : get_avg_color());}
	return cached_avg_color;
}

size_t model3d::get_gpu_mem() const {
	size_t mem(unbound_geom.get_gpu_mem());
	for (auto m = materials.begin(); m != materials.end(); ++m) {mem += m->get_gpu_mem();}
	return mem;
}


int model3d::get_material_ix(string const &material_name, string const &fn, bool okay_if_exists) {

	unsigned mat_id(0);
	string_map_t::const_iterator it(mat_map.find(material_name));

	if (it == mat_map.end()) {
		mat_id = (unsigned)materials.size();
		mat_map[material_name] = mat_id;
		materials.push_back(material_t(material_name, fn));
	}
	else {
		if (!from_model3d_file && !okay_if_exists) {cerr << "Warning: Redefinition of material " << material_name << " in file " << fn << endl;}
		mat_id = it->second;
	}
	assert(mat_id < materials.size());
	return mat_id;
}


int model3d::find_material(string const &material_name) {

	string_map_t::const_iterator it(mat_map.find(material_name));

	if (it == mat_map.end()) {
		if (undef_materials.find(material_name) == undef_materials.end()) {
			if (material_name != "null" && material_name != "<null>") {
				cerr << "Error: Material " << material_name << " not found in any included material libraries" << endl;
			}
			undef_materials.insert(material_name);
		}
		return -1; // return -1 on failure
	}
	assert(it->second < materials.size());
	return it->second;
}


void model3d::mark_mat_as_used(int mat_id) {

	if (mat_id < 0) return;
	assert((unsigned)mat_id < materials.size());
	materials[mat_id].mark_as_used();
}


void model3d::finalize() {

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < (int)materials.size(); ++i) {materials[i].finalize();}
	unbound_geom.finalize();
}


void model3d::clear() {

	free_context();
	unbound_geom.clear();
	materials.clear();
	undef_materials.clear();
	mat_map.clear();
	coll_tree.clear();
	textures_loaded = 0;
}

void model3d::free_context() {

	for (material_t &m : materials) {
		m.geom.free_vbos();
		m.geom_tan.free_vbos();
	}
	unbound_geom.free_vbos();
	clear_smaps();
	free_texture(model_refl_tid);
	free_texture(model_indir_tid);
	if (no_store_model_textures_in_memory) {textures_loaded = 0;} // must reload textures
}

void model3d::clear_smaps() { // frees GL state
	
	for (auto i = smap_data.begin(); i != smap_data.end(); ++i) {
		for (auto j = i->second.begin(); j != i->second.end(); ++j) {j->free_gl_state();}
	}
	smap_data.clear();
}

void model3d::load_all_used_tids() {

	if (textures_loaded) return; // is this safe to skip?
	timer_t timer("Model3d Texture Load", (!tmgr.empty() && !materials.empty()));
#if 0
	// MT loading flow; will fail with an error if any textures require calling texture_t::resize() due to nested OpenGL calls;
	// while this mostly works, it's only ~1s faster because of poor load balancing, and there's a lot of time taken later to compress, build mipmaps, and send to the GPU
	// build a worklist of files to load from disk
	for (auto &m : materials) {m.queue_textures_to_load(tmgr);}
	// load the textures from disk in parallel; requires sorting and uniquing textures, which may be shared across multiple materials
	tmgr.load_work_items_mt();
	// run serial post load steps and make any required OpenGL calls
#endif
	for (auto &m : materials) {m.init_textures(tmgr);}
	textures_loaded = 1;
	if (no_store_model_textures_in_memory) {tmgr.free_client_mem();}
}


void model3d::bind_all_used_tids() {

	load_all_used_tids();
		
	for (material_t &m : materials) {
		if (!m.mat_is_used()) continue;
		m.check_for_tc_invert_y(tmgr);
		tmgr.ensure_tid_bound(m.get_render_texture()); // only one tid for now
		
		if (m.use_bump_map()) {
			if (model_calc_tan_vect && !m.geom.empty()) {
				cerr << "Error loading model3d material " << m.name << ": Geometry is missing tangent vectors, so bump map cannot be enabled." << endl;
				m.bump_tid = -1; // disable bump map
			}
			else {
				tmgr.ensure_tid_bound(m.bump_tid);
			}
			needs_bump_maps = 1;
		}
		if (m.use_spec_map()) {
			tmgr.ensure_tid_bound(m.s_tid);
			tmgr.ensure_tid_bound(m.ns_tid);
			has_spec_maps  |= (m.s_tid  >= 0);
			has_gloss_maps |= (m.ns_tid >= 0);
		}
		needs_alpha_test |= m.get_needs_alpha_test();
		needs_trans_pass |= m.is_partial_transparent();
		has_alpha_mask   |= m.has_alpha_mask();
	} // for m
	calc_tangent_vectors();
}

void model3d::calc_tangent_vectors() {

	for (material_t &m : materials) {
		if (!m.mat_is_used() || !m.use_bump_map()) continue;
		m.geom_tan.calc_tangents();
	}
}

void model3d::simplify_indices(float reduce_target) {
	for (material_t &m : materials) {m.simplify_indices(reduce_target);}
	unbound_geom.simplify_indices(reduce_target);
}

void model3d::reverse_winding_order(uint64_t mats_mask) { // Note: only handles up to 64 materials
	if (mats_mask == 0) return; // nothing to do

	for (unsigned m = 0; m < materials.size(); ++m) {
		if (mats_mask & (uint64_t(1) << m)) {materials[m].reverse_winding_order();}
	}
	if (mats_mask & (uint64_t(1) << materials.size())) {unbound_geom.reverse_winding_order();} // use the last material slot
}


void set_def_spec_map() {
	if (enable_spec_map()) {select_texture(WHITE_TEX, 8);} // all white/specular (no specular map texture)
}
bool cull_front_faces(int reflection_pass) {
	// the reflection pass uses a mirror, which changes the winding direction, so we cull the front faces instead
	return ((reflection_pass == 1) ^ invert_model3d_faces);
}

void model3d::render_materials(shader_t &shader, bool is_shadow_pass, int reflection_pass, bool is_z_prepass, int enable_alpha_mask, unsigned bmap_pass_mask,
	int trans_op_mask, base_mat_t const &unbound_mat, rotation_t const &rot, point const *const xlate, xform_matrix const *const mvm, bool force_lod,
	float model_lod_mult, float fixed_lod_dist, bool skip_cull_face, bool is_scaled, bool no_set_min_alpha, unsigned skip_mat_mask)
{
	bool const is_normal_pass(!is_shadow_pass && !is_z_prepass), is_bmap_pass((bmap_pass_mask & 2) != 0);
	if (is_normal_pass) {smap_data[rot].set_for_all_lights(shader, mvm);} // choose correct shadow map based on rotation

	if (group_back_face_cull && reflection_pass != 2 && !skip_cull_face) { // okay enable culling if is_shadow_pass on some scenes
		if (cull_front_faces(reflection_pass)) {glCullFace(GL_FRONT);}
		glEnable(GL_CULL_FACE);
	}
	// render geom that was not bound to a material; skip if skip_mat_mask was set since this material is non-indexable
	if ((bmap_pass_mask & 1) && unbound_mat.color.alpha > 0.0 && (trans_op_mask & 1) && !unbound_geom.empty() && skip_mat_mask == 0) { // not in bump map only pass; assume opaque
		if (is_normal_pass) { // cur_ub_tid texture shouldn't have an alpha mask, so we don't need to use it in the shadow pass
			assert(unbound_mat.tid >= 0);
			select_texture(unbound_mat.tid);
			shader.set_material(unbound_mat);
			if (!shader.get_user_flag(SHADER_FLAG_NO_ALPHA_TEST)) {shader.add_uniform_float("min_alpha", 0.0);}
			set_def_spec_map();
		}
		if (is_normal_pass || enable_alpha_mask != 1) {unbound_geom.render(shader, is_shadow_pass, xlate);} // skip shadow + alpha mask only pass
		if (is_normal_pass) {shader.clear_specular();}
	}
	bool check_lod(force_lod);
	point center(all_zeros);
	float const lod_thresh(1.0E6*model_mat_lod_thresh*model_lod_mult*lod_scale);

	if ((world_mode == WMODE_INF_TERRAIN || use_model_lod_blocks) && !is_shadow_pass) { // setup LOD/distance culling
		point pts[2] = {bcube.get_llc(), bcube.get_urc()};
		rot.rotate_point(pts[0], -1.0); rot.rotate_point(pts[1], -1.0);
		cube_t const bcube_rot(pts[0], pts[1]);
		// Note: this code assumes the model is in camera space and only rot is used to transform it;
		// if the model is translated, it's up to the caller to translate camera_pdu.pos into model space for LOD to work;
		// if the model is scaled, the containment test below is likely invalid, so we assume the camera is not inside the model and set check_lod=1
		check_lod |= (is_scaled || !bcube_rot.contains_pt(camera_pdu.pos));
		if (check_lod) {center = bcube_rot.get_cube_center();}
	}
	// render all materials (opaque then transparent)
	for (unsigned pass = 0; pass < (is_z_prepass ? 1U : 2U); ++pass) { // opaque, transparent
		if (!(trans_op_mask & (1<<pass))) continue; // wrong opaque vs. transparent pass

		for (unsigned i = 0; i < materials.size(); ++i) {
			if (skip_mat_mask & (1<<i)) continue; // skip this material
			material_t const &mat(materials[i]);

			if (mat.is_partial_transparent() == (pass != 0) && (bmap_pass_mask & (1 << unsigned(mat.use_bump_map())))) {
				if (check_lod && mat.avg_area_per_tri > 0.0 && !mat.no_lod_cull) {
					if ((fixed_lod_dist ? fixed_lod_dist : p2p_dist(camera_pdu.pos, center)) > lod_thresh*mat.avg_area_per_tri) continue; // LOD/dist culling
				}
				to_draw.push_back(make_pair(mat.draw_order_score, i));
			}
		} // for i
		sort(to_draw.begin(), to_draw.end());

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			materials[to_draw[i].second].render(shader, tmgr, unbound_mat.tid, is_shadow_pass, is_z_prepass, enable_alpha_mask, is_bmap_pass, xlate, no_set_min_alpha);
		}
		to_draw.clear();
	} // for pass
	if (group_back_face_cull && reflection_pass != 2 && !skip_cull_face) { // okay enable culling if is_shadow_pass on some scenes
		if (cull_front_faces(reflection_pass)) {glCullFace(GL_BACK);} // restore the default
		glDisable(GL_CULL_FACE);
	}
}

void model3d::render_material(shader_t &shader, unsigned mat_id, bool is_shadow_pass, bool is_z_prepass,
	int enable_alpha_mask, bool is_bmap_pass, point const *const xlate, bool no_set_min_alpha)
{
	if (mat_id == materials.size()) {unbound_geom.render(shader, is_shadow_pass, xlate);} // unbound geom is material ID materials.size() (one past the end)
	else {get_material(mat_id).render(shader, tmgr, unbound_mat.tid, is_shadow_pass, is_z_prepass, enable_alpha_mask, is_bmap_pass, xlate, no_set_min_alpha);}
}

material_t *model3d::get_material_by_name(string const &name) {
	int const mat_ix(find_material(name));
	return ((mat_ix < 0) ? nullptr : &materials[mat_ix]);
}

colorRGBA model3d::set_color_for_material(unsigned mat_id, colorRGBA const &color) { // returns the original color
	colorRGBA orig_color;

	if (mat_id == materials.size()) { // unbound geom is material ID materials.size() (one past the end)
		orig_color = unbound_mat.color;
		unbound_mat.color = color;
	}
	else {
		material_t &mat(get_material(mat_id));
		orig_color = mat.kd;
		mat.ka = mat.kd = color;
	}
	return orig_color;
}

int model3d::set_texture_for_material(unsigned mat_id, int tid) { // returns the original texture ID
	int orig_tid(tid);
	if (mat_id == materials.size()) {swap(orig_tid, unbound_mat.tid);} // unbound geom is material ID materials.size() (one past the end)
	else {swap(orig_tid, get_material(mat_id).d_tid);}
	return orig_tid;
}

void model3d::set_material_emissive_color(unsigned mat_id, colorRGBA const &color) {
	get_material(mat_id).ke = color; // Note: mat_id cannot map to unbound_mat
}


bool geom_xform_t::operator==(geom_xform_t const &x) const {
	if (tv != x.tv || scale != x.scale) return 0;
	UNROLL_3X(if (mirror[i_] != x.mirror[i_]) return 0;)

	for (unsigned i = 0; i < 3; ++i) {
		UNROLL_3X(if (swap_dim[i][i_] != x.swap_dim[i][i_]) return 0;)
	}
	return 1;
}

xform_matrix geom_xform_t::create_xform_matrix() const { // mirror - swap - scale - translate
	glm::mat4 const M_mirror(glm::scale(glm::mat4(1.0), glm::vec3((mirror[0] ? -1.0 : 1.0), (mirror[1] ? -1.0 : 1.0), (mirror[2] ? -1.0 : 1.0))));
	glm::mat4 M_rotate(1.0); // starts as identity
	
	// Warning: the behavior of swap_dim is different from the other functions such as xform_pos_rm();
	// this is intentional, because this function is applied as a pre-transform for animated models rather than for instancing static models,
	// and the rotation behavior makes more sense in that context since the caller will be responsible for object placement and orientation;
	// also, the order of i vs. j matters here as it determines the rotation direction
	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 3; ++j) {
			if (swap_dim[i][j] && i != j) { // skip identity swap
				vector3d A, B;
				A[i] = 1.0;
				B[j] = 1.0;
				M_rotate = glm::rotate(M_rotate, PI_TWO, vec3_from_vector3d(cross_product(A, B)));
			}
		}
	}
	glm::mat4 const M_scale(glm::scale(glm::mat4(1.0), glm::vec3(scale, scale, scale)));
	glm::mat4 const M_translate(glm::translate(glm::mat4(1.0), glm::vec3(tv.x, tv.y, tv.z)));
	return M_mirror * M_rotate * M_scale * M_translate;
}

void model3d_xform_t::apply_inv_xform_to_pdu(pos_dir_up &pdu) const { // Note: RM ignored
	// Note: since pdu's don't have an xform matrix, and don't track applied xforms, we must do the translate first
	pdu.translate(-tv);
	//pdu.rotate(axis, -angle); // FIXME: incorrect - we want to rotate about the model's origin, not the frustum/camera origin
	assert(scale != 0.0);
	//assert(scale > 0.0); // negative scales represent an XYZ mirror, which isn't something we can support for a PDU, so instead we take the abs()
	pdu.scale(1.0/fabs(scale));
	if (angle != 0.0) {pdu.valid = 0;} // since we can't transform the pdu correctly, we give up and disable using it for VFC
}

void model3d_xform_t::apply_to_tquad(coll_tquad &tquad) const {
	for (unsigned i = 0; i < tquad.npts; ++i) {xform_pos_rms(tquad.pts[i]);}
	if (angle != 0.0) {rotate_vector3d_multi(axis, -TO_RADIANS*(double)angle, tquad.pts, tquad.npts);} // negative rotate?
	for (unsigned i = 0; i < tquad.npts; ++i) {tquad.pts[i] += tv;}
	tquad.update_normal(); // simplest to recalculate it
}

cube_t model3d_xform_t::get_xformed_cube(cube_t const &cube) const { // Note: RM ignored
	if (angle == 0.0) {return cube*scale + tv;} // optimization
	return rotate_cube(cube*scale, axis, -TO_RADIANS*angle) + tv; // negative rotate?
}
cube_t const &model3d_xform_t::get_xformed_bcube(cube_t const &bcube) {
	if (bcube_xf.is_all_zeros()) {bcube_xf = get_xformed_cube(bcube);}
	return bcube_xf;
}

void model3d_xform_t::apply_gl() const {
	assert(scale != 0.0);
	translate_to(tv);
	rotation_t::apply_gl();
	
	// Note: unused/untested, but has the same rotation semantics as in create_xform_matrix(); maybe we should just call create_xform_matrix() here?
	for (unsigned i = 0; i < 3; ++i) {
		for (unsigned j = 0; j < 3; ++j) {
			if (swap_dim[i][j] && i != j) { // skip identity swap; the order of i vs. j doesn't matter
				vector3d A, B;
				A[i] = 1.0;
				B[j] = 1.0;
				rotate_about_radians(PI_TWO, cross_product(A, B));
			}
		}
	}
	vector3d scale_xyz(scale, scale, scale);
	UNROLL_3X(if (mirror[i_]) {scale_xyz[i_] = -scale_xyz[i_];}) // Note: untested
	scale_by(scale_xyz);
}

struct camera_pdu_transform_wrapper {

	pos_dir_up prev_pdu, prev_orig_pdu;
	bool active;

	camera_pdu_transform_wrapper(model3d_xform_t const &xf) : active(!xf.is_identity()) {
		if (!active) return;
		prev_pdu      = camera_pdu;
		prev_orig_pdu = orig_camera_pdu;
		xf.apply_inv_xform_to_pdu(camera_pdu);
		xf.apply_inv_xform_to_pdu(orig_camera_pdu);
		fgPushMatrix();
		xf.apply_gl();
	}
	~camera_pdu_transform_wrapper() {
		if (!active) return;
		fgPopMatrix();
		camera_pdu      = prev_pdu;
		orig_camera_pdu = prev_orig_pdu;
	}
};


bool is_cube_visible_to_camera(cube_t const &cube, bool is_shadow_pass, bool animation_enabled) {
	if (animation_enabled) return 1; // TODO: or use a conservative bcube that includes all animations?
	if (!camera_pdu.cube_visible(cube)) return 0;
	if (!(display_mode & 0x08) && !is_shadow_pass) return 1; // check occlusion culling, but allow occlusion culling during the shadow pass
	return !cube_cobj_occluded(camera_pdu.pos, cube);
}


geom_xform_t model3d::fit_to_scene() {
	model3d_xform_t xf;
	if (!transforms.empty()) {cerr << "Error: Can't fit model3d to scene when transforms have been added" << endl; return xf;}
	cube_t scene(get_scene_bounds());
	max_eq(scene.z2(), (scene.z1() + Z_SCENE_SIZE)); // make sure delta Z is at least Z_SCENE_SIZE
	vector3d const model_sz(bcube.get_size()), scene_sz(scene.get_size());
	xf.scale = min(scene_sz.x/model_sz.x, min(scene_sz.y/model_sz.y, scene_sz.z/model_sz.z)); // make sure it fits in all dims
	xf.tv    = scene.get_cube_center() - xf.scale*bcube.get_cube_center();
	transforms.push_back(xf);
	return xf;
}

void model3d::set_target_translate_scale(point const &target_pos, float target_radius, geom_xform_t &xf) const {
	xf.scale = target_radius / (0.5*bcube.max_len());
	xf.tv    = target_pos - xf.scale*bcube.get_cube_center(); // scale is applied before translate
}

void model3d::render_with_xform(shader_t &shader, model3d_xform_t &xf, xform_matrix const &mvm, bool is_shadow_pass,
	int reflection_pass, bool is_z_prepass, int enable_alpha_mask, unsigned bmap_pass_mask, int reflect_mode, int trans_op_mask)
{
	if (!is_cube_visible_to_camera(xf.get_xformed_bcube(bcube), is_shadow_pass, num_animations())) return; // Note: xlate has already been applied to camera_pdu
	// Note: it's simpler and more efficient to inverse transfrom the camera frustum rather than transforming the geom/bcubes
	// Note: currently, only translate is supported (and somewhat scale)
	camera_pdu_transform_wrapper cptw2(xf);
	base_mat_t ub_mat(unbound_mat);
	xf.apply_material_override(ub_mat);
	//point xlate2(xlate); // complex transforms, occlusion culling disabled
	render_materials(shader, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, trans_op_mask, ub_mat, xf, nullptr, &mvm);
	// cptw2 dtor called here
}

// non-const due to vbo caching, normal computation, bcube caching, etc.
void model3d::render(shader_t &shader, bool is_shadow_pass, int reflection_pass, bool is_z_prepass, int enable_alpha_mask,
	unsigned bmap_pass_mask, int reflect_mode, int trans_op_mask, vector3d const &xlate)
{
	assert(trans_op_mask > 0 && trans_op_mask <= 3); // 1 bit = draw opaque, 2 bit = draw transparent
	if (!needs_trans_pass && !(trans_op_mask & 1)) return; // transparent only pass, but no transparent materials
	if (is_shadow_pass && has_animations() && !ENABLE_ANIMATION_SHADOWS) return; // if animated, and no shadow animations, then don't draw the shadow map at all
	if (transforms.empty() && !is_cube_visible_to_camera(bcube+xlate, is_shadow_pass, num_animations())) return;
	
	if (enable_tt_model_indir && world_mode == WMODE_INF_TERRAIN && !is_shadow_pass) {
		if (model_indir_tid == 0) {create_indir_texture();}
		if (model_indir_tid != 0) {
			bind_texture_tu(model_indir_tid, 1); // indir texture uses TU_ID=1
			shader.add_uniform_color("const_indir_color", colorRGB(0,0,0)); // set black indir color - assumes all models will get here, so not reset
		}
		set_local_model_scene_bounds(shader);
	}
	if (reflect_mode) {
		shader.add_uniform_float("metalness", metalness); // may or may not be used
		shader.add_uniform_float("cube_map_reflect_mipmap_level", 0.0); // may or may not be used; should actually be per-material, based on specular exponent
	}
	if (reflective == 2 && !is_shadow_pass && !is_z_prepass) { // cube map reflections
		cube_t const bcube_xf(get_single_transformed_bcube(xlate));

		if (reflection_pass) { // creating the reflection texture
			if (reflection_pass == 2 && !use_interior_cube_map_refl && bcube_xf.get_cube_center() == camera_pdu.pos) return; // skip self reflections
		}
		else if (reflect_mode == 2 && model_refl_tid) { // using the reflection texture
			setup_shader_cube_map_params(shader, bcube_xf, model_refl_tid, model_refl_tsize); // Note: xlate should be all zeros
#if 0 // TESTING
			select_texture(WHITE_TEX);
			set_def_spec_map();
			shader.set_cur_color(WHITE); // or BLACK
			shader.set_specular_color(WHITE, 60.0);
			shader.set_color_e(BLACK);
			//draw_cube(bcube_xf.get_cube_center(), 0.1, 0.1, 0.1, 1);
			draw_subdiv_sphere(bcube_xf.get_cube_center(), 0.1, N_SPHERE_DIV, 1, 1);
			return; // TESTING
#endif
		}
	}
	if (has_animations()) {
		static float cur_tfticks(0.0);
		if (animate2) {cur_tfticks = tfticks;}
		setup_bone_transforms(shader, cur_tfticks/TICKS_PER_SECOND, default_anim_id);
	}
	xform_matrix const mvm(fgGetMVM());
	model3d_xform_t const xlate_xf(xlate);
	camera_pdu_transform_wrapper cptw(xlate_xf);

	// we need the vbo to be created here even in the shadow pass,
	// and the textures are needed for determining whether or not we need to build the tanget_vectors for bump mapping
	bind_all_used_tids();
	if (is_shadow_pass) {select_texture(WHITE_TEX);} // make sure a valid texture is enabled for the shadow pass

	if (transforms.empty()) { // no transforms case
		render_materials_def(shader, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, trans_op_mask, &xlate, &mvm);
	}
	else if (world_mode == WMODE_INF_TERRAIN) {
		//timer_t timer("Draw Models");
		float const view_dist(get_inf_terrain_fog_dist() + bcube.get_max_dim_sz()), view_dist_sq(view_dist*view_dist); // or get_draw_tile_dist()?
		to_draw_xf.clear();
		
		for (unsigned i = 0; i < transforms.size(); ++i) {
			float const dist_sq(distance_to_camera_sq(transforms[i].tv + xlate)); // only use translate; assumes models are approx centered and rotated about their centers
			if (dist_sq < view_dist_sq) {to_draw_xf.emplace_back(dist_sq, i);} // add if not too far away
		}
		if (trans_op_mask < 3) {sort(to_draw_xf.begin(), to_draw_xf.end());} // drawing opaque and trans in separate passes, sort by dist
		vector<cube_t> occluders;
		point const camera(get_camera_pos() - xlate);

		if (!is_shadow_pass && (display_mode & 0x08) != 0 && !occlusion_cube.is_all_zeros()) { // enable occlusion culling
			for (auto i = to_draw_xf.begin(); i != to_draw_xf.end(); ++i) {
				if (transforms[i->second].get_xformed_bcube(bcube).contains_pt(camera)) continue; // skip cases where the camera is inside the cube
				cube_t const occ_cube(transforms[i->second].get_xformed_cube(occlusion_cube));
				if (!camera_pdu.cube_visible(occ_cube)) continue;
				occluders.push_back(occ_cube);
				if (occluders.size() >= 5) break; // at most 5 occluders, starting from closest model, which is largest in screen space
			}
		}
		if (trans_op_mask == 2) {reverse(to_draw_xf.begin(), to_draw_xf.end());} // opaque: front-to-back, transparent: back-to-front

		for (auto i = to_draw_xf.begin(); i != to_draw_xf.end(); ++i) {
			bool skip(0);

			if (!occluders.empty()) {
				cube_t const bc(transforms[i->second].get_xformed_bcube(bcube));
				point pts[8];
				unsigned const ncorners(get_cube_corners(bc.d, pts, camera, 0)); // 8 corners allocated, but only 6 used

				for (cube_t const &j : occluders) {
					bool int_all(1);

					for (unsigned c = 0; c < ncorners; ++c) {
						if (!check_line_clip(camera, pts[c], j.d)) {int_all = 0; break;}
					}
					if (int_all) {skip = 1; break;}
				}
			}
			if (!skip) {render_with_xform(shader, transforms[i->second], mvm, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, reflect_mode, trans_op_mask);}
		}
	}
	else { // ground mode, no sorting or distance culling
		for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {
			render_with_xform(shader, *xf, mvm, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, reflect_mode, trans_op_mask);
		}
	}
	// cptw dtor called here
}

void model3d::set_xform_zval_from_tt_height(bool flatten_mesh) { // set zval to place bottom center of bcube at the mesh surface

	if (xform_zvals_set) return; // already set
	xform_zvals_set = 1;
	timer_t timer("Calc Zvals");
	vector3d const xlate(-xoff2*DX_VAL, -yoff2*DY_VAL, 0.0); // cancel out xoff2/yoff2 translate
	if (transforms.empty()) {transforms.push_back(model3d_xform_t());} // no transforms case - insert identity transform

#pragma omp parallel for schedule(static,1)
	for (int i = 0; i < (int)transforms.size(); ++i) {
		model3d_xform_t &xf(transforms[i]);
		cube_t const &bcube_xf(xf.get_xformed_bcube(bcube));
		point const center(bcube_xf.get_cube_center() + xlate);
		xf.tv.z += get_exact_zval(center.x, center.y) - bcube_xf.d[2][0];
		if (flatten_mesh) {flatten_hmap_region(bcube_xf);} // flatten the mesh under the bcube to a height of mesh_zval
		xf.clear_bcube(); // invalidate and force recompute
	}
}

void model3d::set_sky_lighting_file(string const &fn, float weight, unsigned sz[3]) {
	sky_lighting_fn = fn;
	sky_lighting_weight = weight;
	assert(weight > 0.0);
	UNROLL_3X(assert(sz[i_] > 0); sky_lighting_sz[i_] = sz[i_];)
}

void model3d::create_indir_texture() {

	unsigned const xsize(sky_lighting_sz[0]), ysize(sky_lighting_sz[1]), zsize(sky_lighting_sz[2]), tot_sz(xsize*ysize*zsize);
	
	if (tot_sz == 0) {
		static bool showed_warning(0);
		if (!showed_warning) {cout << "Failed to allocate indir lighting texture: " << TXT(xsize) << TXT(ysize) << TXT(zsize) << endl;}
		showed_warning = 1;
		return; // nothing to do
	}
	timer_t timer("Create Indir Texture");
	lmap_manager_t local_lmap_manager; // store in the model3d and cache for reuse on context change (at the cost of more CPU memory usage)? only matters when ray tracing (below)?
	lmcell init_lmcell;
	unsigned char **need_lmcell = nullptr; // not used - dense mode
	local_lmap_manager.alloc(tot_sz, xsize, ysize, zsize, need_lmcell, init_lmcell);
	float const init_weight(light_int_scale[LIGHTING_SKY]); // record orig value

	if (!sky_lighting_fn.empty() && local_lmap_manager.read_data_from_file(sky_lighting_fn.c_str(), LIGHTING_SKY)) {
		assert(sky_lighting_weight > 0.0);
		light_int_scale[LIGHTING_SKY] = sky_lighting_weight;
	}
	else {
		// TODO: run raytracing to fill in local_lmap_manager, use bcube for bounds (scale to get_scene_bounds_bcube())
	}
	vector<unsigned char> tex_data;
	indir_light_tex_from_lmap(model_indir_tid, local_lmap_manager, tex_data, xsize, ysize, zsize);
	light_int_scale[LIGHTING_SKY] = init_weight; // restore orig value
}

void model3d::set_local_model_scene_bounds(shader_t &s) { // tight bounds
	s.setup_scene_bounds_from_bcube(bcube);
}

void model3d::ensure_reflection_cube_map() {

	if (reflective != 2) return; // no cube map reflections
	bool const dynamic_update(enable_reflection_dynamic_updates());
	if (model_refl_tid && !dynamic_update) return; // reflection texture is valid and scene is static
	cube_t const bcube_xf(get_single_transformed_bcube());
	if (model_refl_tid && !camera_pdu.cube_visible(bcube_xf)) return; // reflection texture is valid and model is not in view
	if (model_refl_tid && (display_mode & 0x08) != 0 && cube_cobj_occluded(get_camera_pos(), bcube_xf)) return; // occlusion culling

	if (indoors == 2) { // not yet known - calculate it (very approximate but should be okay for simple/easy cases)
		int cindex(-1); // unused
		point const test_pt(bcube_xf.get_cube_center());
		indoors = ::check_coll_line(point(test_pt.x, test_pt.y, bcube_xf.d[2][1]+SMALL_NUMBER), point(test_pt.x, test_pt.y, czmax), cindex, -1, 1, 0, 1, 0);
	}
	create_cube_map_reflection(model_refl_tid, model_refl_tsize, model_refl_last_tsize, -1, bcube_xf, (model_refl_tid != 0 && !all_model3d_ref_update), (indoors == 1));
}

cube_t model3d::get_single_transformed_bcube(vector3d const &xlate) const {

	cube_t bcube_xf(bcube + xlate);

	if (!transforms.empty()) {
		if (transforms.size() > 1) { // instancing not supported with a single cube map reflection texture
			cerr << "Error: Instanced models cannot be used with reflection cube maps" << endl;
			exit(1);
		}
		bcube_xf = transforms[0].get_xformed_cube(bcube_xf);
	}
	return bcube_xf;
}


void model3d::model_smap_data_t::render_scene_shadow_pass(point const &lpos) {

	model->bind_all_used_tids();

	for (unsigned sam_pass = 0; sam_pass < (model->get_has_alpha_mask() ? 2U : 1U); ++sam_pass) { // {normal, alpha mask}
		shader_t s;
		s.begin_shadow_map_shader(sam_pass == 1);
		model->render_materials_def(s, 1, 0, 0, (sam_pass == 1), 3, 3, &zero_vector); // no transforms; both opaque and partially transparent
		s.end_shader();
	}
}


void model3d::setup_shadow_maps() {

	if (!shadow_map_enabled()) return; // disabled

	if (smap_data.empty()) { // allocate new shadow maps
		if (transforms.empty()) {smap_data[rotation_t()];} // no transforms case, insert def rotation in map
		for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {smap_data[*xf];}
	}
	for (auto m = smap_data.begin(); m != smap_data.end(); ++m) {
		if (m->second.empty()) {
			for (unsigned i = 0; i < NUM_LIGHT_SRC; ++i) {m->second.push_back(model_smap_data_t(6+i, shadow_map_sz, this));} // uses tu_id 6 and 7
		}
		m->second.create_if_needed(get_bcube(), &m->first); // inverse rotate light source
	}
}


cube_t model3d::calc_bcube_including_transforms() { // non-const because bcube_xf is cached

	if (transforms.empty()) return bcube; // no transforms case
	if (bcube_all_xf != all_zeros_cube) return bcube_all_xf; // already calculated
	
	for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {
		cube_t const &bc(xf->get_xformed_bcube(bcube));
		if (bcube_all_xf == all_zeros_cube) {bcube_all_xf = bc;} else {bcube_all_xf.union_with_cube(bc);}
	}
	return bcube_all_xf;
}


void model3d::build_cobj_tree(bool verbose) {

	if (!coll_tree.is_empty() || has_cobjs) return; // already built or not needed because cobjs will be used instead
	RESET_TIME;
	get_polygons(coll_tree.get_tquads_ref());
	PRINT_TIME(" Get Model3d Polygons");
	coll_tree.build_tree_top(verbose);
	PRINT_TIME(" Cobj Tree Create (from model3d)");
}

bool model3d::check_coll_line_cur_xf(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) {
	if (!check_line_clip(p1, p2, bcube.d)) return 0;
	build_cobj_tree(1);
	return coll_tree.check_coll_line(p1, p2, cpos, cnorm, color, exact);
}

// Note: const as long as build_bvh_if_needed=0
bool model3d::check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact, bool build_bvh_if_needed) {

	if (!build_bvh_if_needed && coll_tree.is_empty()) return 0;
	//timer_t timer("Check Coll Line");
	if (transforms.empty()) {return check_coll_line_cur_xf(p1, p2, cpos, cnorm, color, exact);}
	if (bcube_all_xf != all_zeros_cube && !check_line_clip(p1, p2, bcube_all_xf.d)) return 0; // use bcube of transformed bcubes
	bool coll(0);
	point cur(p2);

	for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {
		if (!check_line_clip(p1, p2, xf->get_xformed_bcube(bcube).d)) continue;
		point p1x(p1), p2x(cur);
		xf->inv_xform_pos(p1x);
		xf->inv_xform_pos(p2x);

		if (check_coll_line_cur_xf(p1x, p2x, cpos, cnorm, color, exact)) { // Note: only modifies cnorm and color if a collision is found
			xf->xform_pos(cpos);
			xf->xform_pos_rm(cnorm);
			coll = 1;
			cur  = cpos; // closer intersection point - shorten the segment
		}
	}
	return coll;
}


void model3d::get_all_mat_lib_fns(set<string> &mat_lib_fns) const {
	for (material_t const &m : materials) {mat_lib_fns.insert(m.filename);}
}

void model3d::compute_area_per_tri() {
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < (int)materials.size(); ++i) {materials[i].compute_area_per_tri();}

	float max_area_per_tri(0.0);
	for (material_t const &m : materials) {max_eq(max_area_per_tri, m.avg_area_per_tri);}
	for (material_t &m : materials) {m.no_lod_cull = (m.avg_area_per_tri == max_area_per_tri);} // material with max area per tri is never culled
}

void model3d::get_stats(model3d_stats_t &stats) const {

	stats.transforms += transforms.size();
	unbound_geom.get_stats(stats);
	
	for (material_t const &m : materials) {
		m.geom.get_stats(stats);
		m.geom_tan.get_stats(stats);
		++stats.mats;
	}
}

void model3d::show_stats() const {

	if (SHOW_MODEL_BCUBE_CENTER) { // more verbose
		cout << "Model stats:" << endl;
		cout << "bcube: " << bcube.str() << endl;
		cout << "center: " << bcube.get_cube_center().str() << ", size: " << bcube.get_size().str() << endl;
	}
	else { // less verbose
		cout << "Model stats: "; // no newline
	}
	model3d_stats_t stats;
	get_stats(stats);
	stats.print();
}


// ************ read/write of animations and models ************

bool model_anim_t::write(ostream &out) const {
	// bone_name_to_index_map
	write_uint(out, bone_name_to_index_map.size());

	for (auto const &kv : bone_name_to_index_map) {
		write_string(out, kv.first );
		write_uint  (out, kv.second);
	}
	if (!out.good()) return 0;
	// transforms
	write_vector(out, bone_transforms); // better be packed
	write_vector(out, bone_offset_matrices);
	write_val   (out, global_inverse_transform);
	write_val   (out, root_transform);
	write_string(out, model_name);
	if (!out.good()) return 0;
	// anim_nodes
	write_uint(out, anim_nodes.size());

	for (anim_node_t const &n : anim_nodes) { // Note: no_anim_data is not written
		write_val   (out, n.bone_index);
		write_string(out, n.name);
		write_val   (out, n.transform);
		write_vector(out, n.children);
	}
	if (!out.good()) return 0;
	// animations
	write_uint(out, animations.size());

	for (animation_t const &a : animations) {
		write_val   (out, a.ticks_per_sec);
		write_val   (out, a.duration);
		write_string(out, a.name);
		write_uint  (out, a.anim_data.size());

		for (auto const &kv : a.anim_data) {
			write_string(out, kv.first );
			anim_data_t const &d(kv.second);
			write_val   (out, d.uses_scale);
			write_vector(out, d.pos);
			write_vector(out, d.scale); // if (uses_scale)?
			write_vector(out, d.rot);
			// Note: not writing d.merged
		} // for kv
	} // for a
	return out.good();
}

bool model_anim_t::read(istream &in) {
	// bone_name_to_index_map
	unsigned const bone_name_to_index_map_sz(read_uint(in));
	string str;

	for (unsigned i = 0; i < bone_name_to_index_map_sz; ++i) {
		read_string(in, str);
		bool const did_ins(bone_name_to_index_map.insert(make_pair(str, read_uint(in))).second);
		assert(did_ins);
	}
	if (!in.good()) return 0;
	// transforms
	read_vector(in, bone_transforms);
	read_vector(in, bone_offset_matrices);
	read_val   (in, global_inverse_transform);
	read_val   (in, root_transform);
	read_string(in, model_name);
	if (!in.good()) return 0;
	// anim_nodes
	anim_nodes.resize(read_uint(in));

	for (anim_node_t &n : anim_nodes) {
		read_val   (in, n.bone_index);
		read_string(in, n.name);
		read_val   (in, n.transform);
		read_vector(in, n.children);
	}
	if (!in.good()) return 0;
	// animations
	animations.resize(read_uint(in));

	for (animation_t &a : animations) {
		read_val   (in, a.ticks_per_sec);
		read_val   (in, a.duration);
		read_string(in, a.name);
		unsigned const anim_data_sz(read_uint(in));

		for (unsigned i = 0; i < anim_data_sz; ++i) {
			read_string(in, str);
			anim_data_t &d(a.anim_data[str]);
			assert(d.pos.empty() && d.scale.empty() && d.rot.empty()); // must be a new entry
			read_val   (in, d.uses_scale); // or calculate from scale?
			read_vector(in, d.pos);
			read_vector(in, d.scale); // if (uses_scale)?
			read_vector(in, d.rot);
			// Note: not reading d.merged
		} // for i
	} // for a
	merge_anim_transforms(); // need to recalculate merged because it's not written or read
	return in.good();
}


bool model3d::write_to_disk(string const &fn) const { // as model3d file; Note: transforms not written

	ofstream out(fn, ios::out | ios::binary);
	
	if (!out.good()) {
		cerr << format_red("Error opening model3d file for write: " + fn) << endl;
		return 0;
	}
	cout << "Writing model3d file " << fn << endl;
	bool const inc_animations(has_animations());
	write_uint(out, (inc_animations ? MAGIC_NUMBER_ANIM : MAGIC_NUMBER));
	out.write((char const *)&bcube, sizeof(cube_t));
	if (!unbound_geom.write(out, inc_animations)) return 0;
	write_uint(out, (unsigned)materials.size());
	
	for (material_t const &m : materials) {
		if (!m.write(out, inc_animations)) {
			cerr << "Error writing material " << m.name << endl;
			return 0;
		}
	}
	if (inc_animations && !model_anim_data.write(out)) {cerr << "Error writing animation data" << endl; return 0;}
	return out.good();
}

bool model3d::read_from_disk(string const &fn) { // as model3d file; Note: transforms not read

	ifstream in(fn, ios::in | ios::binary);
	
	if (!in.good()) {
		cerr << format_red("Error opening model3d file for read: " + fn) << endl;
		return 0;
	}
	clear(); // may not be needed
	unsigned const magic_number_comp(read_uint(in));
	bool const inc_animations(magic_number_comp == MAGIC_NUMBER_ANIM);

	if (!inc_animations && magic_number_comp != MAGIC_NUMBER) {
		cerr << "Error reading model3d file " << fn << ": Invalid file format (magic number check failed)." << endl;
		return 0;
	}
	cout << "Reading model3d file " << fn << endl;
	from_model3d_file = 1;
	in.read((char *)&bcube, sizeof(cube_t));
	if (!unbound_geom.read(in, inc_animations)) return 0;
	materials.resize(read_uint(in));
	
	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->read(in, inc_animations)) {
			cerr << "Error reading material" << endl;
			return 0;
		}
		mat_map[m->name] = (m - materials.begin());
	}
	if (inc_animations && !model_anim_data.read(in)) {cerr << "Error reading animation data" << endl; return 0;}
	//simplify_indices(0.1); // TESTING
	return in.good();
}

bool model3d::write_as_obj_file(string const &fn) {

	ofstream out(fn, ios::out);

	if (!out.good()) {
		cerr << "Error opening obj file for write: " << fn << endl;
		return 0;
	}
	string const mtllib_fn((endswith(fn, ".obj") ? fn.substr(0, fn.size()-4) : fn) + ".mtl"); // replace ".obj" with ".mtl" or append ".mtl"
	cout << "Writing obj file " << fn << endl;
	out << "# Created by 3DWorld, by Frank Gennari 2022" << endl;
	out << "mtllib " << mtllib_fn << endl;
	unsigned cur_vert_ix(0);
	if (!unbound_geom.write_to_obj_file(out, cur_vert_ix)) return 0; // no usemtl

	for (material_t const &m : materials) {
		if (!m.write_to_obj_file(out, cur_vert_ix)) {
			cerr << "Error writing material " << m.name << endl;
			return 0;
		}
	} // for m
	if (!out.good()) return 0;
	// Note: model_anim_data not written
	out.close();

	// write mtllib file
	out.open(mtllib_fn, ios::out);

	if (!out.good()) {
		cerr << "Error opening mtllib file for write: " << mtllib_fn << endl;
		return 0;
	}
	cout << "Writing mtllib file " << mtllib_fn << endl;
	out << "# Created by 3DWorld, by Frank Gennari 2022" << endl;
	for (material_t const &m : materials) {m.write_mtllib_entry(out, tmgr);}
	return out.good();
}


void model3d::proc_model_normals(vector<counted_normal> &cn, int recalc_normals, float nmag_thresh) {
	for (counted_normal &i : cn) {
		if (!i.is_valid()) continue; // invalid, remains invalid
		float const mag(i.mag());
		if (mag < TOLERANCE) {i.count = 0; continue;} // make invalid
		i /= mag; // normalize
		i.count = (recalc_normals > 1 || mag > nmag_thresh*i.count); // stores the 'valid' state of the normal
	}
}

void model3d::proc_model_normals(vector<weighted_normal> &wn, int recalc_normals, float nmag_thresh) { // unused
	for (weighted_normal &i : wn) {
		if (!i.is_valid()) continue; // invalid, remains invalid
		float const mag(i.mag());
		if (mag < TOLERANCE) {i.weight = 0.0; continue;} // make invalid
		i /= mag; // normalize
	}
}


void model3d::write_to_cobj_file(ostream &out) const {

	// 'O' <filename> <group_cobjs_level> <recalc_normals/use_vertex_normals> <write_file> [<voxel_xy_spacing>]
	out << "metalness " << metalness << endl;
	out << "cube_map_ref " << (reflective == 2) << endl;
	out << "l " << 0.5 << " " << unbound_mat.color.raw_str() << " " << texture_str(unbound_mat.tid) << endl; // Note: elastic is hard-coded as 0.5
	out << "r 1.0 " << unbound_mat.shine << " " << unbound_mat.spec_color.raw_str() << endl;
	out << "O " << filename << " " << group_cobjs_level << " " << recalc_normals << " " << 0 << endl; // write_file=0

	for (auto i = transforms.begin(); i != transforms.end(); ++i) {
		out << "l " << 0.5 << " " << i->material.color.raw_str() << " " << texture_str(i->material.tid) << endl; // Note: elastic is hard-coded as 0.5
		out << "r 1.0 " << i->material.shine << " " << i->material.spec_color.raw_str() << endl;
		out << "Z " << i->group_cobjs_level << " " << i->tv.x << " " << i->tv.y << " " << i->tv.z << " " << i->scale << " "
			<< i->axis.x << " " << i->axis.y << " " << i->axis.z << " " << i->angle << " " << i->voxel_spacing << endl;
	}
}


// ************ model3ds ************

void model3ds::clear() {
	for (iterator m = begin(); m != end(); ++m) {m->clear();}
	deque<model3d>::clear();
	tmgr.clear();
}

void model3ds::free_context() {
	for (iterator m = begin(); m != end(); ++m) {m->free_context();}
	tmgr.free_tids();
}


void model3ds::render(bool is_shadow_pass, int reflection_pass, int trans_op_mask, vector3d const &xlate) { // Note: xlate is only used in tiled terrain mode
	
	if (empty()) return;
	bool const tt_mode(world_mode == WMODE_INF_TERRAIN);
	bool const shader_effects(!disable_shader_effects && !is_shadow_pass);
	bool const use_custom_smaps(shader_effects && shadow_map_enabled() && tt_mode);
	bool const enable_any_reflections(shader_effects && !is_shadow_pass && (reflection_pass == 0 || ENABLE_INTER_REFLECTIONS));
	// Note: planar reflections are disabled during the cube map reflection creation pass because they don't work (wrong point is reflected)
	bool const enable_planar_reflections(reflection_pass != 2 && enable_any_reflections && reflection_tid > 0 && use_reflection_plane());
	bool const enable_cube_map_reflections(enable_any_reflections && enable_all_reflections());
	bool const allow_animations(!disable_shader_effects && (!is_shadow_pass || ENABLE_ANIMATION_SHADOWS));
	// Note: in ground mode, lighting is global, so transforms are included in vpos with use_mvm=1; in TT mode, lighting is relative to each model instance
	bool const use_mvm(!tt_mode && has_any_transforms()), v(!tt_mode), use_smap(1 || v);
	bool needs_alpha_test(0), needs_bump_maps(0), any_planar_reflective(0), any_cube_map_reflective(0), any_non_reflective(0);
	bool use_spec_map(0), use_gloss_map(0), needs_trans_pass(0), has_alpha_mask(0);
	unsigned num_models_with_animations(0);

	for (iterator m = begin(); m != end(); ++m) {
		needs_alpha_test |= m->get_needs_alpha_test();
		needs_trans_pass |= m->get_needs_trans_pass();
		has_alpha_mask   |= m->get_has_alpha_mask();
		use_spec_map     |= (enable_spec_map() && m->uses_spec_map());
		use_gloss_map    |= (enable_spec_map() && m->uses_gloss_map());
		if (allow_animations) {num_models_with_animations += m->has_animations();}
		if      (enable_planar_reflections   && m->is_planar_reflective  ()) {any_planar_reflective   = 1;}
		else if (enable_cube_map_reflections && m->is_cube_map_reflective()) {any_cube_map_reflective = 1;}
		else                                                                 {any_non_reflective      = 1;}
		if (shader_effects  ) {needs_bump_maps |= m->get_needs_bump_maps();} // optimization, makes little difference
		if (use_custom_smaps) {m->setup_shadow_maps();} else if (!is_shadow_pass) {m->clear_smaps();}
	}
	if (any_planar_reflective + any_cube_map_reflective > 1) {
		cerr << "Error: Cannot mix planar reflections and cube map reflections for model3ds" << endl;
		exit(1); // unsupported by shaders, but there are currently no cases that get here/need this
	}
	if (!needs_trans_pass && !(trans_op_mask & 1)) return; // transparent only pass, but no transparent materials
	bool const any_animated(num_models_with_animations > 0), all_animated(num_models_with_animations == size());
	shader_t s;
	set_fill_mode();

	if (use_z_prepass && !is_shadow_pass && reflection_pass == 0 && (trans_op_mask & 1) && !all_animated) { // check use_mvm?
		// faster for scenes with high depth complexity and slow fragment shaders; slower when vertex/transform limited
		s.set_prefix("#define POS_FROM_EPOS_MULT", 0); // VS - needed to make transformed vertices agree with the normal rendering flow
		s.begin_color_only_shader(BLACK); // don't even need colors, only need depth
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color rendering, we only want to write to the Z-Buffer
		
		for (iterator m = begin(); m != end(); ++m) {
			if (m->has_animations()) continue; // not supported in z-prepass
			m->render(s, 0, 0, 1, 0, 3, 0, trans_op_mask, xlate);
		}
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		s.end_shader();
		set_std_depth_func_with_eq();
	}
	int const reflect_mode(any_planar_reflective ? 1 : (any_cube_map_reflective ? 2 : 0));
	assert(!reflect_mode || xlate == all_zeros); // xlate not supported for reflections (and not used anyway)
	float const min_alpha(needs_alpha_test ? 0.5 : 0.0); // will be reset per-material, but this variable is used to enable alpha testing

	// the bump map pass is first and the regular pass is second; this way, transparent objects such as glass that don't have bump maps are drawn last
	for (int bmap_pass = (needs_bump_maps ? 1 : 0); bmap_pass >= 0; --bmap_pass) {
		for (unsigned sam_pass = 0; sam_pass < ((is_shadow_pass && has_alpha_mask) ? 2U : 1U); ++sam_pass) {
			for (unsigned ref_pass = (any_non_reflective ? 0U : 1U); ref_pass < (reflect_mode ? 2U : 1U); ++ref_pass) {
				for (unsigned anim_pass = 0; anim_pass < 2; ++anim_pass) {
					if ( anim_pass && !any_animated) continue;
					if (!anim_pass &&  all_animated) continue;
					int const cur_reflect_mode(ref_pass ? reflect_mode : 0);
					bool reset_bscale(0);
					if (anim_pass) {s.add_property("animation_shader", "model_animation.part+");}

					if (is_shadow_pass) {
						if (ENABLE_ANIMATION_SHADOWS && anim_pass) { // Note: should count as a dynamic object for shadow update, but currently does not
							// need to use smoke shader for animations, but disable all shading options
							setup_smoke_shaders(s, min_alpha, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, use_mvm);
						}
						else {s.begin_shadow_map_shader(sam_pass == 1);}
					}
					else if (shader_effects) {
						int const use_bmap((bmap_pass == 0) ? 0 : (model_calc_tan_vect ? 2 : 1));
						int const is_outside((is_shadow_pass || reflection_pass == 1) ? 0 : 2); // enable wet effect coverage mask
						if (model3d_wn_normal) {s.set_prefix("#define USE_WINDING_RULE_FOR_NORMAL", 1);} // FS
						setup_smoke_shaders(s, min_alpha, 0, 0, (enable_tt_model_indir || v), 1, v, v, 0, (use_smap ? 2 : 1), use_bmap, use_spec_map, use_mvm,
							two_sided_lighting, 0.0, model_triplanar_tc_scale, 0, cur_reflect_mode, is_outside, 1, 0, use_gloss_map);
						if (use_custom_smaps) {s.add_uniform_float("z_bias", cobj_z_bias);} // unnecessary?
						if (use_bmap && invert_model_nmap_bscale) {s.add_uniform_float("bump_b_scale", 1.0); reset_bscale = 1;}
						if (ref_pass && any_planar_reflective) {bind_texture_tu(reflection_tid, 14);}
						if (model3d_wn_normal) {s.add_uniform_float("winding_normal_sign", (cull_front_faces(reflection_pass) ? -1.0 : 1.0));}
						s.add_uniform_float("hemi_lighting_scale", model_hemi_lighting_scale);
					}
					else {
						s.begin_simple_textured_shader(0.0, 1); // with lighting
						s.clear_specular();
					}
					for (iterator m = begin(); m != end(); ++m) { // non-const
						if (any_non_reflective && (reflect_mode != 0) && (ref_pass != 0) != m->is_reflective()) continue; // wrong reflection pass for this object
						if (allow_animations && m->has_animations() != (anim_pass != 0)) continue; // wrong animation pass
						m->render(s, is_shadow_pass, reflection_pass, 0, (sam_pass == 1), (shader_effects ? (1 << bmap_pass) : 3), cur_reflect_mode, trans_op_mask, xlate);
					}
					if (reset_bscale) {s.add_uniform_float("bump_b_scale", -1.0);} // may be unnecessary
					if (anim_pass   ) {s.remove_property("animation_shader");}
					s.clear_specular(); // may be unnecessary
					s.end_shader();
				} // anim_pass
			} // ref_pass
		} // sam_pass
	} // bmap_pass
	if (use_z_prepass) {set_std_depth_func();} // reset to default
}

void model3ds::ensure_reflection_cube_maps() {
	if (!enable_all_reflections()) return;
	for (iterator m = begin(); m != end(); ++m) {m->ensure_reflection_cube_map();}
}

void model3ds::set_xform_zval_from_tt_height(bool flatten_mesh) {
	if (!auto_calc_tt_model_zvals) return;
	for (iterator m = begin(); m != end(); ++m) {m->set_xform_zval_from_tt_height(flatten_mesh);}
}

bool model3ds::has_any_transforms() const {
	for (const_iterator m = begin(); m != end(); ++m) {if (m->has_any_transforms()) return 1;}
	return 0;
}
bool model3ds::has_any_animations() const {
	for (const_iterator m = begin(); m != end(); ++m) {if (m->has_animations()) return 1;}
	return 0;
}


cube_t model3ds::calc_and_return_bcube(bool only_reflective) { // Note: calculates bcubes, so non-const

	cube_t bcube(all_zeros_cube); // will return this if empty()
	bool bcube_set(0);

	for (iterator m = begin(); m != end(); ++m) {
		cube_t const bb(m->calc_bcube_including_transforms());
		if (only_reflective && !m->is_planar_reflective()) continue;
		if (!bcube_set) {bcube = bb; bcube_set = 1;} else {bcube.union_with_cube(bb);}
	}
	return bcube;
}

void model3ds::get_all_model_bcubes(vector<cube_t> &bcubes) const {
	for (const_iterator m = begin(); m != end(); ++m) {m->get_transformed_bcubes(bcubes);}
}

size_t model3ds::get_gpu_mem() const {
	size_t mem(tmgr.get_gpu_mem());
	for (const_iterator m = begin(); m != end(); ++m) {mem += m->get_gpu_mem();}
	return mem;
}

void model3ds::build_cobj_trees(bool verbose) {
	for (iterator m = begin(); m != end(); ++m) {m->build_cobj_tree(verbose);}
}


bool model3ds::check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact, bool build_bvh_if_needed) {

	bool ret(0);
	point end_pos(p2);

	for (iterator m = begin(); m != end(); ++m) { // Note: const as long as build_bvh_if_needed=0
		if (m->check_coll_line(p1, end_pos, cpos, cnorm, color, exact, build_bvh_if_needed)) {
			end_pos = cpos; // advance so that we get the closest intersection point to p1
			ret = 1;
		}
	}
	return ret;
}


void model3ds::write_to_cobj_file(ostream &out) const {
	for (const_iterator m = begin(); m != end(); ++m) {m->write_to_cobj_file(out);}
}


void model3d_stats_t::print() const {
	
	cout << "verts: " << verts << ", quads: " << quads << ", tris: " << tris << ", blocks: " << blocks << ", mats: " << mats;
	if (transforms) {cout << ", transforms: " << transforms;}
	cout << endl;
}


// ************ Free Functions ************


void free_model_context() {all_models.free_context();}

void render_models(int shadow_pass, int reflection_pass, int trans_op_mask, vector3d const &xlate) { // shadow_only: 0=non-shadow pass, 1=sun/moon shadow, 2=dynamic shadow
	all_models.render((shadow_pass != 0), reflection_pass, trans_op_mask, xlate);
	if (trans_op_mask & 1) {draw_buildings(shadow_pass, 0, xlate);} // opaque pass (first); Note: not passing reflection_pass (which is for water plane, not mirrors)
	
	if (world_mode == WMODE_INF_TERRAIN) {
		if ((trans_op_mask & 2) && !shadow_pass) {
			draw_building_lights(xlate); // transparent pass (second); not drawn in the shadow pass
			draw_transparent_city_bldg_geom(reflection_pass, xlate);
		}
		draw_cities(shadow_pass, reflection_pass, trans_op_mask, xlate);
	}
}
void ensure_model_reflection_cube_maps() {all_models.ensure_reflection_cube_maps();}
void auto_calc_model_zvals() {all_models.set_xform_zval_from_tt_height(flatten_tt_mesh_under_models);}

model3d &get_cur_model(string const &operation) {

	if (all_models.empty()) {
		cerr << "No current model to " << operation << endl;
		exit(1);
	}
	return all_models.back();
}

void get_cur_model_polygons(vector<coll_tquad> &ppts, model3d_xform_t const &xf, unsigned lod_level) {
	RESET_TIME;
	unsigned const start_ix(ppts.size());
	model3d &cur_model(get_cur_model("extract polygons from"));
	cur_model.get_polygons(ppts, 0, 0, lod_level);
	cur_model.set_has_cobjs();
	xform_polygons(ppts, xf, start_ix);
	PRINT_TIME("Create and Xform Model3d Polygons");
}

size_t get_loaded_models_gpu_mem() {return all_models.get_gpu_mem();}

cube_t get_polygons_bcube(vector<coll_tquad> const &ppts) {

	cube_t bcube(all_zeros);

	for (auto i = ppts.begin(); i != ppts.end(); ++i) { // rasterize ppts to cubes in {x,y,z}
		if (i == ppts.begin()) {bcube = i->get_bcube();} else {bcube.union_with_cube(i->get_bcube());}
	}
	return bcube;
}

void get_cur_model_edges_as_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf) {

	float const grid_spacing(xf.voxel_spacing);
	assert(grid_spacing > 0.0);
	vector<coll_tquad> ppts;
	get_cur_model_polygons(ppts, xf);
	//cube_t const bcube(xf.get_xformed_cube(get_cur_model("get bcube").get_bcube()));
	cube_t const bcube(get_polygons_bcube(ppts));
	vector3d const csz(bcube.get_size());
	unsigned ndiv[3] = {};
	uint64_t tot_grid(1);
	for (unsigned i = 0; i < 3; ++i) {ndiv[i] = max(2U, min(1024U, unsigned(csz[i]/grid_spacing))); tot_grid *= ndiv[i];} // clamp to [2,1024] range
	cout << "polygons: " << ppts.size() << ", grid: " << ndiv[0] << "x" << ndiv[1] << "x" << ndiv[2] << endl;
	if (tot_grid > (1ULL<<20)) {cerr << "Error: Too many model3d grid voxels: " << tot_grid << endl; assert(0);}
	RESET_TIME;
	voxel_grid<cube_t> grid;
	grid.init(ndiv[0], ndiv[1], ndiv[2], bcube, all_zeros_cube);
	vector<point> pts_out;

	for (auto i = ppts.begin(); i != ppts.end(); ++i) { // rasterize ppts to cubes in {x,y,z}
		cube_t const c(i->get_bcube());
		int llc[3], urc[3]; // {low, high} indices
		grid.get_bcube_ix_bounds(c, llc, urc);

		for (int y = llc[1]; y <= urc[1]; ++y) {
			for (int x = llc[0]; x <= urc[0]; ++x) {
				for (int z = llc[2]; z <= urc[2]; ++z) {
					point const llc(grid.get_pt_at(x, y, z));
					cube_t const gc_max(llc, llc+grid.vsz);
					cube_t &gc(grid.get_ref(x, y, z));
					if (gc == gc_max) continue; // already at max

					if (gc_max.contains_cube(c)) { // optimization for contained case
						if (gc == all_zeros_cube) {gc = c;} else {gc.union_with_cube(c);}
						break;
					}
					clip_polygon_to_cube(gc_max, i->pts, i->npts, c, pts_out);
					if (pts_out.empty()) continue; // no intersection

					for (auto p = pts_out.begin(); p != pts_out.end(); ++p) {
						point pt(*p);
						gc_max.clamp_pt(pt); // not required, but needed for FP precision to avoid the assertion below
						if (gc == all_zeros_cube) {gc.set_from_point(pt);} else {gc.union_with_pt(pt);}
					}
					assert(gc_max.contains_cube(gc));
				} // for z
			} // for x
		} // for y
	} // for i
	for (auto i = grid.begin(); i != grid.end(); ++i) {
		if (i->is_near_zero_area()) continue; // skip zero area (volume?) cubes
		if (cubes.empty() || !cubes.back().cube_merge(*i)) {cubes.push_back(*i);}
	}
	PRINT_TIME("Model3d Polygons to Cubes");
	cout << "grid size: " << grid.size() << ", cubes out: " << cubes.size() << endl;
}

//void get_cur_model_edges_as_spheres(vector<sphere_t> &spheres, model3d_xform_t const &xf, float grid_spacing) {}

void get_cur_model_as_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf) { // Note: only xf.scale is used
	RESET_TIME;
	model3d &cur_model(get_cur_model("extract cubes from"));
	cur_model.get_cubes(cubes, xf);
	//cur_model.set_has_cobjs(); // billboard cobjs are not added, and the colors/textures are missing
	PRINT_TIME("Create Model3d Cubes");
}

bool add_transform_for_cur_model(model3d_xform_t const &xf) {
	if (all_models.empty()) {cerr << "Error: No model to add transform to" << endl; return 0;} // nonfatal
	get_cur_model("transform").add_transform(xf);
	return 1;
}
void set_sky_lighting_file_for_cur_model(string const &fn, float weight, unsigned sz[3]) {
	get_cur_model("sky_lighting_file").set_sky_lighting_file(fn, weight, sz);
}
void set_occlusion_cube_for_cur_model(cube_t const &cube) {
	get_cur_model("model_occlusion_cube").set_occlusion_cube(cube);
}
geom_xform_t fit_cur_model_to_scene() {
	return get_cur_model("fit_to_scene").fit_to_scene();
}
bool have_cur_model() {return (!all_models.empty());}

cube_t calc_and_return_all_models_bcube(bool only_reflective) {return all_models.calc_and_return_bcube(only_reflective);}
void get_all_model_bcubes(vector<cube_t> &bcubes) {all_models.get_all_model_bcubes(bcubes);}

void write_models_to_cobj_file(ostream &out) {all_models.write_to_cobj_file(out);}

void adjust_zval_for_model_coll(point &pos, float radius, float mesh_zval, float step_height) { // used for player

	if (pos.z > mesh_zval) { // above the mesh
		assert(step_height >= 0.0);
		vector3d const xlate(get_tiled_terrain_model_xlate());
		assert(xlate.z == 0.0);
		point p1(pos - xlate), p2(p1);
		p1.z += step_height;
		p2.z  = mesh_zval;
		point cpos(pos);
		vector3d cnorm; // unused
		colorRGBA color; // unused
		
		// ray cast from step height above current point down to mesh
		if (all_models.check_coll_line(p1, p2, cpos, cnorm, color, 1, 1)) { // exact=1 build_bvh_if_needed=1
			pos.z = cpos.z; // only update zval

			if (radius > 0.0) { // test 4 points around the center to avoid falling into a small hole
				float const test_dist(radius/SQRT2);

				for (unsigned dim = 0; dim < 2; ++dim) {
					for (unsigned dir = 0; dir < 2; ++dir) {
						point p1b(p1), p2b(p2);
						p1b[dim] += (dir ? -1.0 : 1.0)*test_dist;
						p2b[dim] += (dir ? -1.0 : 1.0)*test_dist;
						if (all_models.check_coll_line(p1b, p2b, cpos, cnorm, color, 1, 1)) {max_eq(pos.z, cpos.z);} // take highest point
					}
				}
			}
			return;
		}
	}
	pos.z = mesh_zval; // no model collision, place on mesh
}

void check_legal_movement_using_model_coll(point const &prev, point &cur, float radius) {
	
	if (prev == cur) return; // no change
	vector3d const dir((cur - prev).get_norm());
	vector3d const xlate(get_tiled_terrain_model_xlate());
	point const p2(cur + dir*radius - xlate); // extend outward by radius in movement direction
	point cpos(cur);
	vector3d cnorm; // unused
	colorRGBA color; // unused

	if (all_models.check_coll_line((prev - xlate), p2, cpos, cnorm, color, 1, 1)) { // exact=1 build_bvh_if_needed=1
		//cur = cpos - dir*radius + xlate; // move radius away from collision point in reverse movement dir
		UNROLL_3X(cur[i_] += fabs(cnorm[i_])*((cpos[i_] + xlate[i_] - dir[i_]*radius) - cur[i_]);) // allow for sliding along walls
	}
}


