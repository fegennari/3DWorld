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
#include <fstream>
#include <queue>

bool const ENABLE_BUMP_MAPS  = 1;
bool const ENABLE_SPEC_MAPS  = 1;
bool const ENABLE_INTER_REFLECTIONS = 1;
unsigned const MAGIC_NUMBER  = 42987143; // arbitrary file signature
unsigned const BLOCK_SIZE    = 32768; // in vertex indices

bool model_calc_tan_vect(1); // slower and more memory but sometimes better quality/smoother transitions

extern bool group_back_face_cull, enable_model3d_tex_comp, disable_shader_effects, texture_alpha_in_red_comp, use_model2d_tex_mipmaps;
extern bool two_sided_lighting, have_indir_smoke_tex, use_core_context, model3d_wn_normal, invert_model_nmap_bscale, use_z_prepass;
extern unsigned shadow_map_sz, reflection_tid;
extern int display_mode, begin_motion;
extern float model3d_alpha_thresh, model3d_texture_anisotropy, model_triplanar_tc_scale, cobj_z_bias;
extern pos_dir_up orig_camera_pdu;
extern bool vert_opt_flags[3];
extern vector<texture_t> textures;


model3ds all_models;


bool enable_bump_map() {return (ENABLE_BUMP_MAPS && !disable_shader_effects && (display_mode & 0x20) == 0);} // enabled by default
bool enable_spec_map() {return (ENABLE_SPEC_MAPS && !disable_shader_effects);}
bool no_sparse_smap_update();


// ************ texture_manager ************

unsigned texture_manager::create_texture(string const &fn, bool is_alpha_mask, bool verbose, bool invert_alpha, bool wrap, bool mirror) {

	assert(!(wrap && mirror)); // can't both be set
	string_map_t::const_iterator it(tex_map.find(fn));

	if (it != tex_map.end()) { // found (already loaded)
		assert(it->second < textures.size());
		return it->second; // check invert_alpha?
	}
	unsigned const tid((unsigned)textures.size());
	tex_map[fn] = tid;
	if (verbose) cout << "creating texture " << fn << endl;
	bool const compress(!is_alpha_mask && enable_model3d_tex_comp);
	// type=read_from_file format=auto width height wrap ncolors use_mipmaps name [do_compress]
	textures.push_back(texture_t(0, 7, 0, 0, wrap, (is_alpha_mask ? 1 : 3), (use_model2d_tex_mipmaps && !is_alpha_mask),
		fn, 0, compress, model3d_texture_anisotropy)); // always RGB wrapped+mipmap (normal map flag set later)
	textures.back().invert_alpha = invert_alpha;
	textures.back().mirror = mirror;
	return tid; // can't fail
}

void texture_manager::clear() {

	free_textures();
	textures.clear();
	tex_map.clear();
}

void texture_manager::free_tids() {
	for (deque<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {t->gl_delete();}
}
void texture_manager::free_textures() {
	for (deque<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {t->free_data();}
}

void texture_manager::ensure_texture_loaded(texture_t &t, int tid, bool is_bump) {

	if (t.is_loaded()) return;
	//if (is_bump) {t.do_compress = 0;} // don't compress normal maps
	t.load(-1);
		
	if (t.alpha_tid >= 0 && t.alpha_tid != tid) { // if alpha is the same texture then the alpha channel should already be set
		ensure_tid_loaded(t.alpha_tid, 0);
		t.copy_alpha_from_texture(get_texture(t.alpha_tid), texture_alpha_in_red_comp);
	}
	if (is_bump) {t.make_normal_map();}
	t.init(); // must be after alpha copy
	assert(t.is_loaded());
}

void texture_manager::bind_alpha_channel_to_texture(int tid, int alpha_tid) {

	if (tid < 0 || alpha_tid < 0) return; // no texture
	assert((unsigned)alpha_tid < textures.size());
	texture_t &t(get_texture(tid));
	assert(t.ncolors == 3 || t.ncolors == 4);
	if (t.alpha_tid == alpha_tid) return; // already bound
	assert(tid < BUILTIN_TID_START); // can't modify builtin textures
	assert(t.alpha_tid < 0); // can't rebind to a different value
	assert(!t.is_allocated()); // must not yet be loaded
	t.alpha_tid = alpha_tid;
	t.ncolors   = 4; // add alpha channel
	if (t.use_mipmaps) {t.use_mipmaps = 3;} // generate custom alpha mipmaps
}

texture_t &get_builtin_texture(int tid) {
	assert((unsigned)tid < textures.size());
	return textures[tid];
}
texture_t const &texture_manager::get_texture(int tid) const {
	if (tid >= BUILTIN_TID_START) {return get_builtin_texture(tid - BUILTIN_TID_START);} // global textures lookup
	assert((unsigned)tid < textures.size());
	return textures[tid]; // local textures lookup
}
texture_t &texture_manager::get_texture(int tid) {
	if (tid >= BUILTIN_TID_START) {return get_builtin_texture(tid - BUILTIN_TID_START);} // global textures lookup
	assert((unsigned)tid < textures.size());
	return textures[tid]; // local textures lookup
}


// ************ read/write code ************

void write_uint(ostream &out, unsigned val) {
	out.write((const char *)&val, sizeof(unsigned));
}

unsigned read_uint(istream &in) {
	unsigned val;
	in.read((char *)&val, sizeof(unsigned));
	return val;
}

template<typename V> void write_vector(ostream &out, V const &v) {
	write_uint(out, (unsigned)v.size());
	out.write((const char *)&v.front(), v.size()*sizeof(V::value_type));
}

template<typename V> void read_vector(istream &in, V &v) {
	v.clear();
	v.resize(read_uint(in));
	in.read((char *)&v.front(), v.size()*sizeof(V::value_type));
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


template<typename T> void indexed_vntc_vect_t<T>::subdiv_recur(vector<unsigned> const &ixs, unsigned npts, unsigned skip_dims) {

	unsigned const num(ixs.size());
	assert(num > 0 && (num % npts) == 0);
	point const &v0(at(ixs.front()).v);
	cube_t bc(v0, v0);
	
	for (vector<unsigned>::const_iterator i = ixs.begin()+1; i != ixs.end(); ++i) {
		bc.union_with_pt(at(*i).v); // update bounding cube
	}
	if (num > BLOCK_SIZE) { // subdiv case
		float max_sz(0), sval(0);
		unsigned const dim(bc.get_split_dim(max_sz, sval, skip_dims));

		if (max_sz > 0) { // can split
			vector<unsigned> bins[2];
			for (unsigned i = 0; i < 2; ++i) {bins[i].reserve(num/2);}

			for (unsigned i = 0; i < num; i += npts) {
				bool const bix(at(ixs[i]).v[dim] > sval); // use the first point to determine the bin
				for (unsigned j = i; j < i+npts; ++j) {bins[bix].push_back(ixs[j]);}
			}
			if (bins[0].empty() || bins[1].empty()) {skip_dims |= (1 << dim);}

			for (unsigned i = 0; i < 2; ++i) {
				if (!bins[i].empty()) {subdiv_recur(bins[i], npts, skip_dims);}
			}
			return;
		}
	}
	blocks.push_back(geom_block_t(indices.size(), num, bc));
	copy(ixs.begin(), ixs.end(), back_inserter(indices)); // make leaf
}


template<typename T> void indexed_vntc_vect_t<T>::optimize(unsigned npts) {

	vntc_vect_t<T>::optimize(npts);

	if (vert_opt_flags[0]) { // only if not subdivided?
		vert_optimizer optimizer(indices, size(), npts);
		optimizer.run(vert_opt_flags[1], vert_opt_flags[2]);
	}
}


template<typename T> void indexed_vntc_vect_t<T>::finalize(unsigned npts) {

	if (need_normalize) {
		for (iterator i = begin(); i != end(); ++i) {i->n.normalize();}
		need_normalize = 0;
	}
	if (indices.empty() || finalized) return; // nothing to do

	bool const do_simplify = 0; // TESTING, maybe this doesn't really go here
	if (do_simplify && npts == 3) {
		for (unsigned n = 0; n < 1; ++n) {
			vector<unsigned> v;
			simplify(v, 0.5);
			indices.swap(v);
		}
	}
	finalized = 1;
	//optimize(npts); // now optimized in the object file loading phase
	assert((num_verts() % npts) == 0); // triangles or quads
	assert(blocks.empty());

	if (num_verts() > 2*BLOCK_SIZE) { // subdivide large buffers
		vector<unsigned> ixs;
		ixs.swap(indices);
		subdiv_recur(ixs, npts, 0);
	}
}


template<unsigned N> struct vert_to_tri_t {

	unsigned t[N], n; // if this vertex is used in more than N triangles we give up and never remove it

	vert_to_tri_t() : n(0) {}
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
		unsigned new_ixs[3];
		for (unsigned n = 0; n < 3; ++n) {new_ixs[n] = remap.get_remapped_val(indices[i+n]);}
		if (new_ixs[0] == new_ixs[1] || new_ixs[1] == new_ixs[2] || new_ixs[2] == new_ixs[0]) continue; // duplicate vertices, degenerate triangle, skip
		UNROLL_3X(out.push_back(new_ixs[i_]);)
	}
	cout << TXT(num_verts) << TXT(num_ixs) << TXT(target_num_verts) << TXT(cand_verts) << TXT(num_valid_verts) << TXT(out.size()) << endl;
	PRINT_TIME("Simplify");
}


template<typename T> void indexed_vntc_vect_t<T>::clear() {
	
	vntc_vect_t<T>::clear();
	indices.clear();
	blocks.clear();
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


template<typename T> void indexed_vntc_vect_t<T>::render(shader_t &shader, bool is_shadow_pass, point const *const xlate, unsigned npts, bool no_vfc) {

	if (empty()) return;
	assert(npts == 3 || npts == 4);
	finalize(npts);
	if (bsphere.radius == 0.0) {calc_bounding_volumes();}
	//if (is_shadow_pass && vbo == 0 && world_mode == WMODE_GROUND) return; // don't create the vbo on the shadow pass (voxel terrain problems - works now?)

	if (no_vfc) {
		// do nothing
	}
	else if (is_shadow_pass) { // Note: makes shadow map caching more difficult
		if (no_sparse_smap_update() && !orig_camera_pdu.projected_cube_visible(bcube, camera_pdu.pos)) return; // light_pos == camera_pdu.pos for the shadow pass
	}
	else if (vbo) { // don't cull if vbo hasn't yet been allocated because this will cause it to be skipped in the shadow pass
		if (!camera_pdu.sphere_and_cube_visible_test(bsphere.pos, bsphere.radius, bcube)) return; // view frustum culling
		
		if (indices.size() >= 100 && xlate != nullptr && (display_mode & 0x08) != 0) { // Note: null xlate implies there are transforms other than translate, so skip occlusion culling
			if (cube_cobj_occluded((camera_pdu.pos + *xlate), (bcube + *xlate))) return; // occlusion culling
		}
	}
	assert(!indices.empty()); // now always using indexed drawing
	int prim_type(GL_TRIANGLES);
	unsigned ixn(1), ixd(1);

	if (use_core_context && npts == 4) {
		if (!ivbo) {
			vector<unsigned> tixs;
			convert_quad_ixs_to_tri_ixs(indices, tixs); // Note: can use geometry shader, see http://github.prideout.net/quad-meshes
			create_and_upload(*this, tixs);
		}
		ixn = 6; ixd = 4; // convert quads to 2 triangles
	}
	else {
		if (npts == 4) {prim_type = GL_QUADS;}
		create_and_upload(*this, indices);
	}
	pre_render();
	// FIXME: we need this call here because we don't know if the VAO was created with the same enables/locations: consider normal vs. shadow pass
	T::set_vbo_arrays(); // calls check_mvm_update()

	if (is_shadow_pass || blocks.empty() || no_vfc || camera_pdu.sphere_completely_visible_test(bsphere.pos, bsphere.radius)) { // draw the entire range
		glDrawRangeElements(prim_type, 0, (unsigned)size(), (unsigned)(ixn*indices.size()/ixd), GL_UNSIGNED_INT, 0);
	}
	else { // draw each block independently
		// could use glDrawElementsIndirect(), but the draw calls don't seem to add any significant overhead for the current set of models
		for (vector<geom_block_t>::const_iterator i = blocks.begin(); i != blocks.end(); ++i) {
			if (camera_pdu.cube_visible(i->bcube)) {
				glDrawRangeElements(prim_type, 0, (unsigned)size(), ixn*i->num/ixd, GL_UNSIGNED_INT, (void *)((ixn*i->start_ix/ixd)*sizeof(unsigned)));
			}
		}
	}
	post_render();
	T::unset_attrs();
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
	vertex_map_t<T>::const_iterator it(vmap.find(v2));
	unsigned ix;

	if (it == vmap.end()) { // not found
		ix = (unsigned)size();
		push_back(v);
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


struct shared_vertex_t {
	unsigned ai, bi;
	bool shared;
	shared_vertex_t() : ai(0), bi(0), shared(0) {}
	shared_vertex_t(unsigned ai_, unsigned bi_) : ai(ai_), bi(bi_), shared(1) {}
};


template<typename T> void indexed_vntc_vect_t<T>::get_polygons(get_polygon_args_t &args, unsigned npts) const {

	if (args.lod_level > 1 && !indices.empty()) {
		indexed_vntc_vect_t<T> simplified_this(*this); // FIXME: inefficient to copy everything
		simplify(simplified_this.indices, 1.0/args.lod_level);
		get_polygon_args_t args2(args);
		args2.lod_level = 0;
		simplified_this.get_polygons(args2, npts);
		return;
	}
	unsigned const nv(num_verts());
	if (nv == 0) return;
	assert((nv % npts) == 0);
	polygon_t poly(args.color), quad_poly(args.color);
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
			assert(npts == 4);
			args.polygons.push_back(poly);
		}
		else {
			split_polygon(poly, args.polygons, POLY_COPLANAR_THRESH);
		}
	}
}


template<typename T> void indexed_vntc_vect_t<T>::write(ostream &out) const {

	vntc_vect_t<T>::write(out);
	write_vector(out, indices);
}


template<typename T> void indexed_vntc_vect_t<T>::read(istream &in) {

	vntc_vect_t<T>::read(in);
	read_vector(in, indices);
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


template<typename T> void vntc_vect_block_t<T>::optimize(unsigned npts) {
	for (iterator i = begin(); i != end(); ++i) {i->optimize(npts);}
}

template<typename T> void vntc_vect_block_t<T>::free_vbos() {
	for (iterator i = begin(); i != end(); ++i) {i->clear_vbos();}
}


template<typename T> cube_t vntc_vect_block_t<T>::get_bcube() const {

	if (empty()) return all_zeros_cube;
	cube_t bcube(front().get_bcube());
	for (const_iterator i = begin()+1; i != end(); ++i) {bcube.union_with_cube(i->get_bcube());}
	return bcube;
}


template<typename T> unsigned vntc_vect_block_t<T>::num_verts() const {

	unsigned s(0);
	for (const_iterator i = begin(); i != end(); ++i) {s += i->num_verts();}
	return s;
}


template<typename T> unsigned vntc_vect_block_t<T>::num_unique_verts() const {

	unsigned s(0);
	for (const_iterator i = begin(); i != end(); ++i) {s += (unsigned)i->size();}
	return s;
}


template<typename T> float vntc_vect_block_t<T>::calc_draw_order_score() const {

	float area(0.0);
	unsigned count(0);

	for (const_iterator i = begin(); i != end(); ++i) {
		count += i->size();
		area  += i->get_bradius()*i->get_bradius();
	}
	return ((count == 0) ? 0.0 : -area/count);
}


template<typename T> void vntc_vect_block_t<T>::get_polygons(get_polygon_args_t &args, unsigned npts) const {
	for (const_iterator i = begin(); i != end(); ++i) {i->get_polygons(args, npts);}
}


template<typename T> bool vntc_vect_block_t<T>::write(ostream &out) const {

	write_uint(out, (unsigned)size());
	for (const_iterator i = begin(); i != end(); ++i) {i->write(out);}
	return 1;
}


template<typename T> bool vntc_vect_block_t<T>::read(istream &in) {

	clear();
	resize(read_uint(in));
	for (iterator i = begin(); i != end(); ++i) {i->read(in);}
	return 1;
}


// ************ geometry_t ************


template<> void geometry_t<vert_norm_tc_tan>::calc_tangents_blocks(vntc_vect_block_t<vert_norm_tc_tan> &blocks, unsigned npts) {

	for (vntc_vect_block_t<vert_norm_tc_tan>::iterator i = blocks.begin(); i != blocks.end(); ++i) {
		i->calc_tangents(npts);
	}
}

template<typename T> void geometry_t<T>::calc_tangents() {

	calc_tangents_blocks(triangles, 3);
	calc_tangents_blocks(quads,     4);
}


template<typename T> void geometry_t<T>::render_blocks(shader_t &shader, bool is_shadow_pass, point const *const xlate, vntc_vect_block_t<T> &blocks, unsigned npts) {

	for (vntc_vect_block_t<T>::iterator i = blocks.begin(); i != blocks.end(); ++i) {
		i->render(shader, is_shadow_pass, xlate, npts);
	}
}


template<typename T> void geometry_t<T>::render(shader_t &shader, bool is_shadow_pass, point const *const xlate) {

	render_blocks(shader, is_shadow_pass, xlate, triangles, 3);
	render_blocks(shader, is_shadow_pass, xlate, quads,     4);
}


template<typename T> void geometry_t<T>::add_poly_to_polys(polygon_t const &poly, vntc_vect_block_t<T> &v, vertex_map_t<T> &vmap, unsigned obj_id) const {

	if (v.empty() || v.back().size() > MAX_VMAP_SIZE || obj_id > v.back().obj_id) {
		vmap.clear();
		v.push_back(indexed_vntc_vect_t<T>(obj_id));
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
	quads.get_polygons    (args, 4);
}


template<typename T> cube_t geometry_t<T>::get_bcube() const {

	cube_t bcube(all_zeros_cube); // will return this if empty

	if (!triangles.empty()) {
		bcube = triangles.get_bcube();
		if (!quads.empty()) bcube.union_with_cube(quads.get_bcube());
	}
	else if (!quads.empty()) {
		bcube = quads.get_bcube();
	}
	return bcube;
}


template<typename T> void geometry_t<T>::clear() {

	free_vbos();
	triangles.clear();
	quads.clear();
}


template<typename T> void geometry_t<T>::get_stats(model3d_stats_t &stats) const {
	
	stats.tris  += triangles.num_verts()/3;
	stats.quads += quads.num_verts()/4;
	triangles.get_stats(stats);
	quads.get_stats(stats);
}


// ************ material_t ************


template<typename T> void update_score(vntc_vect_block_t<T> const &v, float &score, unsigned &num_nonempty) {
	if (v.empty()) return;
	score += v.calc_draw_order_score();
	++num_nonempty;
}


void material_t::init_textures(texture_manager &tmgr) {

	if (!mat_is_used()) return;
	int const tid(get_render_texture());
	tmgr.bind_alpha_channel_to_texture(tid, alpha_tid);
	tmgr.ensure_tid_loaded(tid, 0); // only one tid for now
	if (use_bump_map()) tmgr.ensure_tid_loaded(bump_tid, 1);
	if (use_spec_map()) tmgr.ensure_tid_loaded(s_tid, 0);
	might_have_alpha_comp |= tmgr.might_have_alpha_comp(tid);
}


void material_t::render(shader_t &shader, texture_manager const &tmgr, int default_tid, bool is_shadow_pass, bool is_z_prepass, bool enable_alpha_mask, point const *const xlate) {

	if ((geom.empty() && geom_tan.empty()) || skip || alpha == 0.0) return; // empty or transparent
	if (is_shadow_pass && alpha < MIN_SHADOW_ALPHA) return;

	if (draw_order_score == 0) {
		unsigned num_nonempty;
		update_score(geom.triangles,     draw_order_score, num_nonempty);
		update_score(geom.quads,         draw_order_score, num_nonempty);
		update_score(geom_tan.triangles, draw_order_score, num_nonempty);
		update_score(geom_tan.quads,     draw_order_score, num_nonempty);
		assert(num_nonempty > 0);
		draw_order_score /= num_nonempty; // take the average
	}
	int const tex_id(get_render_texture());

	if (is_z_prepass) { // no textures
		if (alpha < 1.0 || (tex_id >= 0 && alpha_tid >= 0)) return; // partially transparent or has alpha mask
		geom.render(shader, 0, xlate);
		geom_tan.render(shader, 0, xlate);
	}
	else if (is_shadow_pass) {
		bool const has_alpha_mask(tex_id >= 0 && alpha_tid >= 0);
		if (has_alpha_mask != enable_alpha_mask) return; // incorrect pass
		if (has_alpha_mask) {tmgr.bind_texture(tex_id);} // enable alpha mask texture
		geom.render(shader, 1, xlate);
		geom_tan.render(shader, 1, xlate);
		if (has_alpha_mask) {select_texture(WHITE_TEX);} // back to a default white texture
	}
	else {
		bool has_binary_alpha(1);
		
		if (tex_id >= 0) {
			tmgr.bind_texture(tex_id);
			has_binary_alpha = tmgr.has_binary_alpha(tex_id);
		}
		else {
			select_texture((default_tid >= 0) ? default_tid : WHITE_TEX); // no texture specified - use white texture
		}
		if (use_bump_map()) {
			set_active_texture(5);
			tmgr.bind_texture(bump_tid);
			set_active_texture(0);
		}
		if (enable_spec_map()) { // all white/specular if no specular map texture
			set_active_texture(8);
			if (s_tid >= 0) {tmgr.bind_texture(s_tid);} else {select_texture(WHITE_TEX);}
			set_active_texture(0);
		}
		//if (!disable_shader_effects && alpha < 1.0 && ni != 1.0) {shader.add_uniform_float("refract_ix", ni);} // FIXME: set index of refraction
		bool const need_blend(is_partial_transparent()); // conservative, but should be okay
		if (need_blend) {enable_blend();}
		float const min_alpha(min(0.99*alpha, ((alpha_tid >= 0) ? (has_binary_alpha ? 0.9 : model3d_alpha_thresh) : 0.0)));
		shader.add_uniform_float("min_alpha", min_alpha);
		if (ns > 0.0) {shader.set_specular_color(ks, ns);} // ns<=0 is undefined?
		shader.set_color_e(colorRGBA(ke, alpha));
		// Note: ka is ignored here because it represents a "fake" lighting model;
		// 3DWorld uses a more realistic lighting model where ambient comes from indirect lighting that's computed independently from the material;
		// however, it might make sense to use ka instead of ke when ke is not specified?
		shader.set_cur_color(get_ad_color());
		geom.render(shader, 0, xlate);
		geom_tan.render(shader, 0, xlate);
		shader.clear_color_e();
		if (ns > 0.0) {shader.clear_specular();}
		if (need_blend) {disable_blend();}
		//if (!disable_shader_effects && alpha < 1.0 && ni != 1.0) {shader.add_uniform_float("refract_ix", 1.0);}
	}
}


bool material_t::use_bump_map() const {

	return (enable_bump_map() && bump_tid >= 0);
}


bool material_t::use_spec_map() const {

	return (enable_spec_map() && s_tid >= 0);
}


colorRGBA material_t::get_ad_color() const {

	colorRGBA c(kd, alpha);
	c.set_valid_color();
	return c;
}


colorRGBA material_t::get_avg_color(texture_manager const &tmgr, int default_tid) const {

	colorRGBA avg_color(get_ad_color());
	int tex_id(get_render_texture());
	
	if (tex_id >= 0) {
		return avg_color.modulate_with(tmgr.get_tex_avg_color(tex_id));
	}
	else if (default_tid >= 0) {
		return avg_color.modulate_with(texture_color(default_tid));
	}
	return avg_color;
}


bool material_t::add_poly(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], unsigned obj_id) {
	
	if (skip) return 0;

	if (model_calc_tan_vect && use_bump_map()) {
		geom_tan.add_poly(poly, vmap_tan, obj_id);
	}
	else {
		geom.add_poly(poly, vmap, obj_id);
	}
	mark_as_used();
	return 1;
}


bool material_t::write(ostream &out) const {

	out.write((char const *)this, sizeof(material_params_t));
	write_vector(out, name);
	write_vector(out, filename);
	return (geom.write(out) && geom_tan.write(out));
}


bool material_t::read(istream &in) {

	in.read((char *)this, sizeof(material_params_t));
	read_vector(in, name);
	read_vector(in, filename);
	return (geom.read(in) && geom_tan.read(in));
}


// ************ model3d ************


void coll_tquads_from_triangles(vector<triangle> const &triangles, vector<coll_tquad> &ppts, colorRGBA const &color) {

	ppts.reserve(ppts.capacity() + triangles.size());
	for (unsigned i = 0; i < triangles.size(); ++i) ppts.push_back(coll_tquad(triangles[i], color));
}


unsigned model3d::add_triangles(vector<triangle> const &triangles, colorRGBA const &color, int mat_id, unsigned obj_id) {

	// average_normals=1 should turn most of these face normals into vertex normals
	vntc_map_t  vmap    [2] = {vntc_map_t (1), vntc_map_t (1)};
	vntct_map_t vmap_tan[2] = {vntct_map_t(1), vntct_map_t(1)};
	unsigned tot_added(0);
	polygon_t poly(color);

	for (vector<triangle>::const_iterator i = triangles.begin(); i != triangles.end(); ++i) {
		poly.from_triangle(*i);
		tot_added += add_polygon(poly, vmap, vmap_tan, mat_id, obj_id);
	}
	return tot_added;
}


unsigned model3d::add_polygon(polygon_t const &poly, vntc_map_t vmap[2], vntct_map_t vmap_tan[2], int mat_id, unsigned obj_id) {
	
	for (unsigned d = 0; d < 2; ++d) {
		vmap[d].check_for_clear(mat_id);
		vmap_tan[d].check_for_clear(mat_id);
	}
	split_polygons_buffer.resize(0);
	split_polygon(poly, split_polygons_buffer, 0.0);

	for (vector<polygon_t>::iterator i = split_polygons_buffer.begin(); i != split_polygons_buffer.end(); ++i) {
		if (mat_id < 0) {
			unbound_geom.add_poly(*i, vmap, obj_id);
		}
		else {
			assert((unsigned)mat_id < materials.size());
			materials[mat_id].add_poly(*i, vmap, vmap_tan, obj_id);
		}
	}
	update_bbox(poly);
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
	update_bbox(tri);
}


void model3d::update_bbox(polygon_t const &poly) {
	cube_t const bb(get_polygon_bbox(poly));
	if (bcube == all_zeros_cube) {bcube = bb;} else {bcube.union_with_cube(bb);}
}


void model3d::get_polygons(vector<coll_tquad> &polygons, bool quads_only, bool apply_transforms, unsigned lod_level) const {

	unsigned const start_pix(polygons.size());

	if (start_pix == 0) { // Note: we count quads as 1.5 polygons because some of them may be split into triangles
		model3d_stats_t stats;
		get_stats(stats);
		unsigned const num_copies((!apply_transforms || transforms.empty()) ? 1 : transforms.size());
		polygons.reserve(num_copies*(quads_only ? stats.quads : (stats.tris + 1.5*stats.quads)));
	}
	colorRGBA def_color(WHITE);
	if (unbound_mat.tid >= 0) {def_color.modulate_with(texture_color(unbound_mat.tid));}
	get_polygon_args_t args(polygons, def_color, quads_only, lod_level);
	unbound_geom.get_polygons(args);

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		args.color = m->get_avg_color(tmgr, unbound_mat.tid);
		m->geom.get_polygons    (args);
		m->geom_tan.get_polygons(args);
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
	float f;
	bool d;
	float_plus_dir() {}
	float_plus_dir(float f_, bool d_) : f(f_), d(d_) {}
	bool operator<(float_plus_dir const &fd) const {return ((f == fd.f) ? (d < fd.d) : (f < fd.f));}
};


template<typename T> unsigned add_polygons_to_voxel_grid(vector<coll_tquad> &polygons, T const &cont,
	vector<vector<float_plus_dir> > &zvals, int bounds[2][2], int num_xy[2], float spacing, unsigned &nhq, model3d_xform_t const &xf)
{
	model3d_stats_t stats;
	cont.get_stats(stats);
	polygons.clear();
	polygons.reserve(stats.quads);
	get_polygon_args_t args(polygons, WHITE, 1); // quads_only=1
	cont.get_polygons(args);
	xform_polygons(polygons, xf, 0);
	
	for (vector<coll_tquad>::const_iterator i = polygons.begin(); i != polygons.end(); ++i) {
		assert(i->npts == 4);
		if (fabs(i->normal.z) < 0.99) continue; // only keep top/bottom cube sides
		cube_t const bcube(i->get_bcube());
		if ((bcube.d[2][1] - bcube.d[2][0]) > 0.5*spacing) continue; // can this happen?
		int cbounds[2][2];
		calc_bounds(bcube, cbounds, spacing);
		bool const is_top(i->normal.z > 0.0);
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
void model3d::get_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf, float spacing) const {

	RESET_TIME;
	assert(spacing > 0.0);

	// calculate scene voxel bounds
	int bounds[2][2], num_xy[2]; // {x,y}x{lo,hi}
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

		for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
			num_polys += add_polygons_to_voxel_grid(polygons, m->geom,     zvals, bounds, num_xy, spacing, num_horiz_quads, xf);
			num_polys += add_polygons_to_voxel_grid(polygons, m->geom_tan, zvals, bounds, num_xy, spacing, num_horiz_quads, xf);
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
			cerr << "Error: Material " << material_name << " not found in any included material libraries" << endl;
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


void model3d::optimize() {
	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {m->optimize();}
	unbound_geom.optimize();
}


void model3d::clear() {

	free_context();
	unbound_geom.clear();
	materials.clear();
	undef_materials.clear();
	mat_map.clear();
	coll_tree.clear();
	smap_data.clear(); // unnecessary
}


void model3d::free_context() {

	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		m->geom.free_vbos();
		m->geom_tan.free_vbos();
	}
	unbound_geom.free_vbos();
	clear_smaps();
	free_texture(model_refl_tid);
}


void model3d::load_all_used_tids() {
	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {m->init_textures(tmgr);}
}


void model3d::bind_all_used_tids() {

	load_all_used_tids();
		
	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->mat_is_used()) continue;
		tmgr.ensure_tid_bound(m->get_render_texture()); // only one tid for now
		
		if (m->use_bump_map()) {
			if (model_calc_tan_vect && !m->geom.empty()) {
				cerr << "Error loading model3d material " << m->name << ": Geometry is missing tangent vectors, so bump map cannot be enabled." << endl;
				m->bump_tid = -1; // disable bump map
			}
			else {
				tmgr.ensure_tid_bound(m->bump_tid);
				m->geom_tan.calc_tangents();
			}
			needs_bump_maps = 1;
		}
		if (m->use_spec_map()) {
			tmgr.ensure_tid_bound(m->s_tid);
			has_spec_maps = 1;
		}
		needs_alpha_test |= m->get_needs_alpha_test();
	}
}


void set_def_spec_map() {
	if (enable_spec_map()) {select_multitex(WHITE_TEX, 8);} // all white/specular (no specular map texture)
}

void model3d::render_materials(shader_t &shader, bool is_shadow_pass, int reflection_pass, bool is_z_prepass, bool enable_alpha_mask,
	unsigned bmap_pass_mask, base_mat_t const &unbound_mat, point const *const xlate, xform_matrix const *const mvm)
{
	bool const is_normal_pass(!is_shadow_pass && !is_z_prepass);
	if (is_normal_pass) {smap_data.set_for_all_lights(shader, mvm);}

	if (group_back_face_cull && reflection_pass != 2) { // okay enable culling if is_shadow_pass on some scenes
		if (reflection_pass == 1) {glCullFace(GL_FRONT);} // the reflection pass uses a mirror, which changes the winding direction, so we cull the front faces instead
		glEnable(GL_CULL_FACE);
	}

	// render geom that was not bound to a material
	if ((bmap_pass_mask & 1) && unbound_mat.color.alpha > 0.0) { // enabled, not in bump map only pass
		if (is_normal_pass) { // cur_ub_tid texture shouldn't have an alpha mask, so we don't need to use it in the shadow pass
			assert(unbound_mat.tid >= 0);
			select_texture(unbound_mat.tid);
			shader.set_material(unbound_mat);
			shader.add_uniform_float("min_alpha", 0.0);
			set_def_spec_map();
		}
		if (is_normal_pass || !enable_alpha_mask) {unbound_geom.render(shader, is_shadow_pass, xlate);} // skip shadow + alpha mask pass
		if (is_normal_pass) {shader.clear_specular();}
	}
	
	// render all materials (opaque then transparent)
	vector<pair<float, unsigned> > to_draw;

	for (unsigned pass = 0; pass < (is_z_prepass ? 1U : 2U); ++pass) { // opaque, transparent
		for (unsigned i = 0; i < materials.size(); ++i) {
			if (materials[i].is_partial_transparent() == (pass != 0) && (bmap_pass_mask & (1 << unsigned(materials[i].use_bump_map())))) {
				to_draw.push_back(make_pair(materials[i].draw_order_score, i));
			}
		}
		sort(to_draw.begin(), to_draw.end());

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			materials[to_draw[i].second].render(shader, tmgr, unbound_mat.tid, is_shadow_pass, is_z_prepass, enable_alpha_mask, xlate);
		}
		to_draw.clear();
	}
	if (group_back_face_cull && reflection_pass != 2) { // okay enable culling if is_shadow_pass on some scenes
		if (reflection_pass == 1) {glCullFace(GL_BACK);} // restore the default
		glDisable(GL_CULL_FACE);
	}
}


bool geom_xform_t::operator==(geom_xform_t const &x) const {
	if (tv != x.tv || scale != x.scale) return 0;
	UNROLL_3X(if (mirror[i_] != x.mirror[i_]) return 0;)

	for (unsigned i = 0; i < 3; ++i) {
		UNROLL_3X(if (swap_dim[i][i_] != x.swap_dim[i][i_]) return 0;)
	}
	return 1;
}

void model3d_xform_t::apply_inv_xform_to_pdu(pos_dir_up &pdu) const { // Note: RM ignored
	// Note: since pdu's don't have an xform matrix, and don't track applied xforms, we must do the translate first
	assert(scale != 0.0);
	pdu.translate(-tv);
	//pdu.rotate(axis, -angle); // FIXME: incorrect - we want to rotate about the model's origin, not the frustum/camera origin
	pdu.scale(1.0/fabs(scale)); // FIXME: what to do about negative scales?
	if (angle != 0.0) {pdu.valid = 0;} // since we can't transform the pdu correctly, we give up and disable using it for VFC
}

void model3d_xform_t::apply_to_tquad(coll_tquad &tquad) const {
	for (unsigned i = 0; i < tquad.npts; ++i) {xform_pos_rms(tquad.pts[i]);}
	if (angle != 0.0) {rotate_vector3d_multi(axis, -TO_RADIANS*angle, tquad.pts, tquad.npts);} // negative rotate?
	for (unsigned i = 0; i < tquad.npts; ++i) {tquad.pts[i] += tv;}
	tquad.update_normal(); // simplest to recalculate it
}

cube_t model3d_xform_t::get_xformed_cube(cube_t const &cube) const { // Note: RM ignored
	if (angle == 0.0) {return cube*scale + tv;} // optimization
	return rotate_cube(cube*scale, axis, -TO_RADIANS*angle) + tv; // negative rotate?
}

void model3d_xform_t::apply_gl() const {
	assert(scale != 0.0);
	translate_to(tv);
	rotation_t::apply_gl();
	for (unsigned i = 0; i < 3; ++i) {UNROLL_3X(assert(!swap_dim[i][i_]);)} // swap not supported
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


bool is_cube_visible_to_camera(cube_t const &cube, bool is_shadow_pass) {
	if (!camera_pdu.cube_visible(cube)) return 0;
	if (!(display_mode & 0x08) && !is_shadow_pass) return 1; // check occlusion culling, but allow occlusion culling during the shadow pass
	return !cube_cobj_occluded(camera_pdu.pos, cube);
}


void model3d::set_target_translate_scale(point const &target_pos, float target_radius, geom_xform_t &xf) const {
	xf.scale = target_radius / (0.5*bcube.max_len());
	xf.tv    = target_pos - xf.scale*bcube.get_cube_center(); // scale is applied before translate
}

void model3d::render_with_xform(shader_t &shader, model3d_xform_t const &xf, xform_matrix const &mvm, bool is_shadow_pass,
	int reflection_pass, bool is_z_prepass, bool enable_alpha_mask, unsigned bmap_pass_mask, int reflect_mode)
{
	if (!is_cube_visible_to_camera(xf.get_xformed_cube(bcube), is_shadow_pass)) return; // Note: xlate has already been applied to camera_pdu
	// Note: it's simpler and more efficient to inverse transfrom the camera frustum rather than transforming the geom/bcubes
	// Note: currently, only translate is supported (and somewhat scale)
	camera_pdu_transform_wrapper cptw2(xf);
	base_mat_t ub_mat(unbound_mat);
	xf.apply_material_override(ub_mat);
	//point xlate2(xlate); // complex transforms, occlusion culling disabled
	render_materials(shader, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, ub_mat, nullptr, &mvm);
	// cptw2 dtor called here
}

// non-const due to vbo caching, normal computation, etc.
void model3d::render(shader_t &shader, bool is_shadow_pass, int reflection_pass, bool is_z_prepass, bool enable_alpha_mask, unsigned bmap_pass_mask, int reflect_mode, vector3d const &xlate) {

	if (transforms.empty() && !is_cube_visible_to_camera(bcube+xlate, is_shadow_pass)) return;
	
	if (reflect_mode) {
		shader.add_uniform_float("metalness", metalness); // may or may not be used
		shader.add_uniform_float("cube_map_reflect_mipmap_level", 0.0); // may or may not be used; should actually be per-material, based on specular exponent
	}
	if (reflective == 2 && !is_shadow_pass && !is_z_prepass) { // cube map reflections
		cube_t const bcube_xf(get_single_transformed_bcube(xlate));

		if (reflection_pass) { // creating the reflection texture
			if (reflection_pass == 2 && bcube_xf.get_cube_center() == camera_pdu.pos) return; // skip self reflections
		}
		else if (reflect_mode == 2 && model_refl_tid) { // using the reflection texture
			setup_shader_cube_map_params(shader, bcube_xf, model_refl_tid); // Note: xlate should be all zeros
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
	xform_matrix const mvm(fgGetMVM());
	model3d_xform_t const xlate_xf(xlate);
	camera_pdu_transform_wrapper cptw(xlate_xf);

	// we need the vbo to be created here even in the shadow pass,
	// and the textures are needed for determining whether or not we need to build the tanget_vectors for bump mapping
	bind_all_used_tids();

	if (transforms.empty()) { // no transforms case
		render_materials_def(shader, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, &xlate, &mvm);
	}
	for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {
		render_with_xform(shader, *xf, mvm, is_shadow_pass, reflection_pass, is_z_prepass, enable_alpha_mask, bmap_pass_mask, reflect_mode);
	}
	// cptw dtor called here
}

void model3d::ensure_reflection_cube_map() {

	if (reflective != 2) return; // no cube map reflections
	bool const dynamic_update(begin_motion != 0); // FIXME: do something better
	if (model_refl_tid && !dynamic_update) return; // reflection texture is valid and scene is static
	cube_t const bcube_xf(get_single_transformed_bcube());
	if (model_refl_tid && !camera_pdu.cube_visible(bcube_xf)) return; // reflection texture is valid and model is not in view
	if (model_refl_tid && (display_mode & 0x08) != 0 && cube_cobj_occluded(get_camera_pos(), bcube_xf)) return; // occlusion culling

	if (indoors == 2) { // not yet known - calculate it (very approximate but should be okay for simple/easy cases)
		int cindex(-1); // unused
		point const test_pt(bcube_xf.get_cube_center());
		indoors = ::check_coll_line(point(test_pt.x, test_pt.y, bcube_xf.d[2][1]+SMALL_NUMBER), point(test_pt.x, test_pt.y, czmax), cindex, -1, 1, 0, 1, 0);
	}
	create_cube_map_reflection(model_refl_tid, -1, bcube_xf, (model_refl_tid != 0), (indoors == 1));
}

cube_t model3d::get_single_transformed_bcube(vector3d const &xlate) const {

	cube_t bcube_xf(bcube + xlate);

	if (!transforms.empty()) {
		assert(transforms.size() == 1); // FIXME: instancing not supported with a single cube map refelction texture
		bcube_xf = transforms[0].get_xformed_cube(bcube_xf);
	}
	return bcube_xf;
}


void setup_smap_shader(shader_t &s, bool sam_pass) {

	if (sam_pass == 1) {
		s.begin_simple_textured_shader(MIN_SHADOW_ALPHA);
		select_texture(WHITE_TEX);
	}
	else {
		s.begin_color_only_shader(); // really don't even need colors
	}
}


void model3d::model_smap_data_t::render_scene_shadow_pass(point const &lpos) {

	model->bind_all_used_tids();

	for (unsigned sam_pass = 0; sam_pass < 2U; ++sam_pass) {
		shader_t s;
		setup_smap_shader(s, (sam_pass != 0));
		model->render_materials_def(s, 1, 0, 0, (sam_pass == 1), 3, &zero_vector); // no transforms
		s.end_shader();
	}
}


void model3d::setup_shadow_maps() {

	if (!shadow_map_enabled()) return; // disabled

	if (smap_data.empty()) { // allocate new shadow maps
		for (unsigned i = 0; i < NUM_LIGHT_SRC; ++i) {smap_data.push_back(model_smap_data_t(6+i, shadow_map_sz, this));} // uses tu_id 6 and 7
	}
	smap_data.create_if_needed(get_bcube());
}


cube_t model3d::calc_bcube_including_transforms() const {

	if (transforms.empty()) {return bcube;} // no transforms case
	cube_t bcube_xf(all_zeros_cube); // will return this if transforms.empty()
	
	for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {
		cube_t const bc(xf->get_xformed_cube(bcube));
		if (bcube_xf == all_zeros_cube) {bcube_xf = bc;} else {bcube_xf.union_with_cube(bc);}
	}
	return bcube_xf;
}


void model3d::build_cobj_tree(bool verbose) {

	if (!coll_tree.is_empty() || has_cobjs) return; // already built or not needed because cobjs will be used instead
	RESET_TIME;
	get_polygons(coll_tree.get_tquads_ref());
	coll_tree.build_tree_top(verbose);
	PRINT_TIME(" Cobj Tree Create (from model3d)");
}


bool model3d::check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) const {

	if (transforms.empty()) {return coll_tree.check_coll_line(p1, p2, cpos, cnorm, color, exact);}
	bool coll(0);
	point cur(p2);

	// Note: this case unused/untested
	for (auto xf = transforms.begin(); xf != transforms.end(); ++xf) {
		point p1x(p1), p2x(cur);
		xf->inv_xform_pos(p1x);
		xf->inv_xform_pos(p2x);

		if (coll_tree.check_coll_line(p1x, p2x, cpos, cnorm, color, exact)) { // Note: only modifies cnorm and color if a collision is found
			xf->xform_pos(cpos);
			xf->xform_pos_rm(cnorm);
			coll = 1;
			cur  = cpos; // closer intersection point - shorten the segment
		}
	}
	return coll;
}


void model3d::get_all_mat_lib_fns(set<string> &mat_lib_fns) const {

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		mat_lib_fns.insert(m->filename);
	}
}


void model3d::get_stats(model3d_stats_t &stats) const {

	stats.transforms += transforms.size();
	unbound_geom.get_stats(stats);
	
	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		m->geom.get_stats(stats);
		m->geom_tan.get_stats(stats);
		++stats.mats;
	}
}


void model3d::show_stats() const {

	model3d_stats_t stats;
	get_stats(stats);
	stats.print();
}


bool model3d::write_to_disk(string const &fn) const { // Note: transforms not written

	ofstream out(fn, ios::out | ios::binary);
	
	if (!out.good()) {
		cerr << "Error opening model3d file for write: " << fn << endl;
		return 0;
	}
	cout << "Writing model3d file " << fn << endl;
	write_uint(out, MAGIC_NUMBER);
	out.write((char const *)&bcube, sizeof(cube_t));
	if (!unbound_geom.write(out)) return 0;
	write_uint(out, (unsigned)materials.size());
	
	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->write(out)) {
			cerr << "Error writing material" << endl;
			return 0;
		}
	}
	return out.good();
}


bool model3d::read_from_disk(string const &fn) { // Note: transforms not read

	ifstream in(fn, ios::in | ios::binary);
	
	if (!in.good()) {
		cerr << "Error opening model3d file for read: " << fn << endl;
		return 0;
	}
	clear(); // ???
	unsigned const magic_number_comp(read_uint(in));

	if (magic_number_comp != MAGIC_NUMBER) {
		cerr << "Error reading model3d file " << fn << ": Invalid file format (magic number check failed)." << endl;
		return 0;
	}
	cout << "Reading model3d file " << fn << endl;
	from_model3d_file = 1;
	in.read((char *)&bcube, sizeof(cube_t));
	if (!unbound_geom.read(in)) return 0;
	materials.resize(read_uint(in));
	
	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->read(in)) {
			cerr << "Error reading material" << endl;
			return 0;
		}
		mat_map[m->name] = (m - materials.begin());
	}
	return in.good();
}


void model3d::proc_counted_normals(vector<counted_normal> &cn, float nmag_thresh) {

	for (vector<counted_normal>::iterator i = cn.begin(); i != cn.end(); ++i) { // if recalc_normals
		if (!i->is_valid()) continue; // invalid, remains invalid
		// Note: could assign higher weights to triangles with smaller area
		*i /= (float)i->count;
		float const mag(i->mag());
		if (mag < 1E-6) {i->count = 0; continue;} // invalid
		assert(mag < 1.001);
		*i /= mag; // normalize
		i->count = (mag > nmag_thresh); // stores the 'valid' state of the normal
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


void model3ds::render(bool is_shadow_pass, int reflection_pass, vector3d const &xlate) { // Note: xlate is only used in tiled terrain mode
	
	if (empty()) return;
	bool const shader_effects(!disable_shader_effects && !is_shadow_pass);
	bool const use_custom_smaps(shader_effects && shadow_map_enabled() && world_mode == WMODE_INF_TERRAIN);
	bool const enable_any_reflections(shader_effects && !is_shadow_pass && (reflection_pass == 0 || ENABLE_INTER_REFLECTIONS));
	// Note: planar reflections are disabled during the cube map reflection creation pass because they don't work (wrong point is reflected)
	bool const enable_planar_reflections(reflection_pass != 2 && enable_any_reflections && reflection_tid > 0 && use_reflection_plane());
	bool const enable_cube_map_reflections(enable_any_reflections && enable_all_reflections());
	bool const use_mvm(has_any_transforms()), v(world_mode == WMODE_GROUND), use_smap(1 || v);
	bool needs_alpha_test(0), needs_bump_maps(0), any_planar_reflective(0), any_cube_map_reflective(0), any_non_reflective(0), use_spec_map(0);
	shader_t s;
	set_fill_mode();

	for (iterator m = begin(); m != end(); ++m) {
		needs_alpha_test |= m->get_needs_alpha_test();
		use_spec_map     |= (enable_spec_map() && m->uses_spec_map());
		if      (enable_planar_reflections   && m->is_planar_reflective  ()) {any_planar_reflective   = 1;}
		else if (enable_cube_map_reflections && m->is_cube_map_reflective()) {any_cube_map_reflective = 1;}
		else                                                                 {any_non_reflective      = 1;}
		if (shader_effects  ) {needs_bump_maps |= m->get_needs_bump_maps();} // optimization, makes little difference
		if (use_custom_smaps) {m->setup_shadow_maps();} else if (!is_shadow_pass) {m->clear_smaps();}
	}
	if (any_planar_reflective + any_cube_map_reflective > 1) {
		cerr << "Error: Cannot mix planar reflections and cube map reflections for model3ds" << endl;
		exit(1); // FIXME: better/earlier error? make this work?
	}
	if (use_z_prepass && !is_shadow_pass && reflection_pass == 0) { // check use_mvm?
		// faster for scenes with high depth complexity and slow fragment shaders; slower when vertex/transform limited
		s.set_prefix("#define POS_FROM_EPOS_MULT", 0); // VS - needed to make transformed vertices agree with the normal rendering flow
		s.begin_color_only_shader(BLACK); // don't even need colors, only need depth
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); // Disable color rendering, we only want to write to the Z-Buffer
		for (iterator m = begin(); m != end(); ++m) {m->render(s, 0, 0, 1, 0, 3, 0, xlate);}
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		s.end_shader();
		glDepthFunc(GL_LEQUAL);
	}
	int const reflect_mode(any_planar_reflective ? 1 : (any_cube_map_reflective ? 2 : 0));
	assert(!reflect_mode || xlate == all_zeros); // xlate not supported for reflections (and not used anyway)

	// the bump map pass is first and the regular pass is second; this way, transparent objects such as glass that don't have bump maps are drawn last
	for (int bmap_pass = (needs_bump_maps ? 1 : 0); bmap_pass >= 0; --bmap_pass) {
		for (unsigned sam_pass = 0; sam_pass < (is_shadow_pass ? 2U : 1U); ++sam_pass) {
			for (unsigned ref_pass = (any_non_reflective ? 0U : 1U); ref_pass < (reflect_mode ? 2U : 1U); ++ref_pass) {
				int const cur_reflect_mode(ref_pass ? reflect_mode : 0);
				bool reset_bscale(0);

				if (is_shadow_pass) {
					setup_smap_shader(s, (sam_pass != 0));
				}
				else if (shader_effects) {
					int const use_bmap((bmap_pass == 0) ? 0 : (model_calc_tan_vect ? 2 : 1));
					float const min_alpha(needs_alpha_test ? 0.5 : 0.0); // will be reset per-material, but this variable is used to enable alpha testing
					int const is_outside((is_shadow_pass || reflection_pass == 1) ? 0 : 2); // enable wet effect coverage mask
					if (model3d_wn_normal) {s.set_prefix("#define USE_WINDING_RULE_FOR_NORMAL", 1);} // FS
					setup_smoke_shaders(s, min_alpha, 0, 0, v, 1, v, v, 0, use_smap, use_bmap, use_spec_map, use_mvm, two_sided_lighting,
						0.0, model_triplanar_tc_scale, 0, cur_reflect_mode, is_outside);
					if (use_custom_smaps) {s.add_uniform_float("z_bias", cobj_z_bias);} // unnecessary?
					if (use_bmap && invert_model_nmap_bscale) {s.add_uniform_float("bump_b_scale", 1.0); reset_bscale = 1;}
					if (ref_pass && any_planar_reflective) {bind_texture_tu(reflection_tid, 14);}
					if (model3d_wn_normal) {s.add_uniform_float("winding_normal_sign", ((reflection_pass == 1) ? -1.0 : 1.0));}
				}
				else {
					s.begin_simple_textured_shader(0.0, 1); // with lighting
					s.clear_specular();
				}
				for (iterator m = begin(); m != end(); ++m) { // non-const
					if (any_non_reflective && (reflect_mode != 0) && (ref_pass != 0) != m->is_reflective()) continue; // wrong reflection pass for this object
					m->render(s, is_shadow_pass, reflection_pass, 0, (sam_pass == 1), (shader_effects ? (1 << bmap_pass) : 3), cur_reflect_mode, xlate);
				}
				if (reset_bscale) {s.add_uniform_float("bump_b_scale", -1.0);} // may be unnecessary
				s.clear_specular(); // may be unnecessary
				s.end_shader();
			} // ref_pass
		} // sam_pass
	} // bmap_pass
	if (use_z_prepass) {glDepthFunc(GL_LESS);} // reset to default
}

void model3ds::ensure_reflection_cube_maps() {
	if (!enable_all_reflections()) return;
	for (iterator m = begin(); m != end(); ++m) {m->ensure_reflection_cube_map();}
}


bool model3ds::has_any_transforms() const {
	for (const_iterator m = begin(); m != end(); ++m) {if (m->has_any_transforms()) return 1;}
	return 0;
}


cube_t model3ds::get_bcube(bool only_reflective) const {

	cube_t bcube(all_zeros_cube); // will return this if empty()
	bool bcube_set(0);

	for (const_iterator m = begin(); m != end(); ++m) {
		cube_t const bb(m->calc_bcube_including_transforms());
		if (only_reflective && !m->is_planar_reflective()) continue;
		if (!bcube_set) {bcube = bb; bcube_set = 1;} else {bcube.union_with_cube(bb);}
	}
	return bcube;
}


void model3ds::build_cobj_trees(bool verbose) {
	for (iterator m = begin(); m != end(); ++m) {m->build_cobj_tree(verbose);}
}


bool model3ds::check_coll_line(point const &p1, point const &p2, point &cpos, vector3d &cnorm, colorRGBA &color, bool exact) const {

	bool ret(0);
	point end_pos(p2);

	for (const_iterator m = begin(); m != end(); ++m) {
		if (m->check_coll_line(p1, end_pos, cpos, cnorm, color, exact)) {
			end_pos = cpos; // advance so that we get the closest intersection point to p1
			ret = 1;
		}
	}
	return ret;
}


void model3d_stats_t::print() const {
	
	cout << "verts: " << verts << ", quads: " << quads << ", tris: " << tris << ", blocks: " << blocks << ", mats: " << mats;
	if (transforms) {cout << ", transforms: " << transforms;}
	cout << endl;
}


// ************ Free Functions ************


void free_model_context() {
	all_models.free_context();
}
void render_models(bool shadow_pass, int reflection_pass, vector3d const &xlate) {
	all_models.render(shadow_pass, reflection_pass, xlate);
}
void ensure_model_reflection_cube_maps() {
	all_models.ensure_reflection_cube_maps();
}

model3d &get_cur_model(string const &operation) {

	if (all_models.empty()) {
		cerr << "No current model to " << operation << endl;
		exit(1);
	}
	return all_models.back();
}

void xform_polygons(vector<coll_tquad> &ppts, model3d_xform_t const &xf, unsigned start_ix=0) {
	if (xf.is_identity()) return;
	#pragma omp parallel for schedule(static,1)
	for (int i = start_ix; i < (int)ppts.size(); ++i) {xf.apply_to_tquad(ppts[i]);}
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

cube_t get_polygons_bcube(vector<coll_tquad> const &ppts) {

	cube_t bcube(all_zeros, all_zeros);

	for (auto i = ppts.begin(); i != ppts.end(); ++i) { // rasterize ppts to cubes in {x,y,z}
		if (i == ppts.begin()) {bcube = i->get_bcube();} else {bcube.union_with_cube(i->get_bcube());}
	}
	return bcube;
}

void get_cur_model_edges_as_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf, float grid_spacing) {

	assert(grid_spacing > 0.0);
	vector<coll_tquad> ppts;
	get_cur_model_polygons(ppts, xf);
	//cube_t const bcube(xf.get_xformed_cube(get_cur_model("get bcube").get_bcube()));
	cube_t const bcube(get_polygons_bcube(ppts));
	vector3d const csz(bcube.get_size());
	unsigned ndiv[3];
	for (unsigned i = 0; i < 3; ++i) {ndiv[i] = max(2U, min(1024U, unsigned(csz[i]/grid_spacing)));} // clamp to [2,1024] range
	cout << "polygons: " << ppts.size() << ", grid: " << ndiv[0] << "x" << ndiv[1] << "x" << ndiv[2] << endl;
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

void get_cur_model_as_cubes(vector<cube_t> &cubes, model3d_xform_t const &xf, float voxel_xy_spacing) { // Note: only xf.scale is used
	RESET_TIME;
	model3d &cur_model(get_cur_model("extract cubes from"));
	cur_model.get_cubes(cubes, xf, voxel_xy_spacing);
	//cur_model.set_has_cobjs(); // billboard cobjs are not added, and the colors/textures are missing
	PRINT_TIME("Create Model3d Cubes");
}

void add_transform_for_cur_model(model3d_xform_t const &xf) {
	get_cur_model("transform").add_transform(xf);
}

cube_t get_all_models_bcube(bool only_reflective) {return all_models.get_bcube(only_reflective);}


