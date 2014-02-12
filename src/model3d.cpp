// 3D World - 3D Model Rendering Code
// by Frank Gennari
// 8/17/11

#include "model3d.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include "voxels.h"
#include "vertex_opt.h"
#include <fstream>

bool const ENABLE_BUMP_MAPS  = 1;
bool const ENABLE_SPEC_MAPS  = 1;
bool const USE_INDEXED_VERTS = 1;
bool const CALC_TANGENT_VECT = 1; // slower and more memory but sometimes better quality/smoother transitions
unsigned const MAGIC_NUMBER  = 42987143; // arbitrary file signature
unsigned const BLOCK_SIZE    = 32768; // in vertex indices

extern bool group_back_face_cull, enable_model3d_tex_comp, disable_shaders, texture_alpha_in_red_comp, use_model2d_tex_mipmaps;
extern bool two_sided_lighting, have_indir_smoke_tex;
extern int display_mode;
extern float model3d_alpha_thresh;
extern pos_dir_up orig_camera_pdu;
extern bool vert_opt_flags[3];


model3ds all_models;


bool enable_bump_map() {return (ENABLE_BUMP_MAPS && !disable_shaders && (display_mode & 0x20) == 0);} // enabled by default
bool enable_spec_map() {return (ENABLE_SPEC_MAPS && !disable_shaders);}
bool no_sparse_smap_update();


// ************ texture_manager ************

unsigned texture_manager::create_texture(string const &fn, bool is_alpha_mask, bool verbose) {

	string_map_t::const_iterator it(tex_map.find(fn));

	if (it != tex_map.end()) { // found (already loaded)
		assert(it->second < textures.size());
		return it->second;
	}
	unsigned const tid((unsigned)textures.size());
	tex_map[fn] = tid;
	if (verbose) cout << "creating texture " << fn << endl;
	bool const compress(!is_alpha_mask && enable_model3d_tex_comp);
	// type=read_from_file format=auto width height wrap ncolors use_mipmaps name [do_compress]
	textures.push_back(texture_t(0, 7, 0, 0, 1, (is_alpha_mask ? 1 : 3), (use_model2d_tex_mipmaps && !is_alpha_mask), fn, 0, compress)); // always RGB wrapped+mipmap
	return tid; // can't fail
}


void texture_manager::clear() {

	free_textures();
	textures.clear();
	tex_map.clear();
}


void texture_manager::free_tids() {

	for (deque<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {
		t->gl_delete();
	}
}

void texture_manager::free_textures() {

	for (deque<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {
		t->free_data();
	}
}


void texture_manager::ensure_texture_loaded(texture_t &t, int tid, bool is_bump) {

	if (!t.is_allocated()) {
		t.load(-1);
		
		if (t.alpha_tid >= 0 && t.alpha_tid != tid) { // if alpha is the same texture then the alpha channel should already be set
			ensure_tid_loaded(t.alpha_tid, 0);
			assert((unsigned)t.alpha_tid < textures.size());
			t.copy_alpha_from_texture(textures[t.alpha_tid], texture_alpha_in_red_comp);
		}
		if (is_bump) t.make_normal_map();
		t.init(); // must be after alpha copy
	}
	assert(t.is_allocated());
}


void texture_manager::bind_alpha_channel_to_texture(int tid, int alpha_tid) {

	if (tid < 0 || alpha_tid < 0) return; // no texture
	assert((unsigned)tid < textures.size() && (unsigned)alpha_tid < textures.size());
	texture_t &t(textures[tid]);
	assert(t.ncolors == 3 || t.ncolors == 4);
	if (t.alpha_tid == alpha_tid) return; // already bound
	assert(t.alpha_tid < 0); // can't rebind to a different value
	assert(!t.is_allocated()); // must not yet be loaded
	t.alpha_tid = alpha_tid;
	t.ncolors   = 4; // add alpha channel
	if (t.use_mipmaps) t.use_mipmaps = 3; // generate custom alpha mipmaps
}


void texture_manager::ensure_tid_loaded(int tid, bool is_bump) {

	if (tid < 0) return; // not allocated
	assert((unsigned)tid < textures.size());
	ensure_texture_loaded(textures[tid], tid, is_bump);
}


void texture_manager::ensure_tid_bound(int tid) {

	if (tid < 0) return; // not allocated
	assert((unsigned)tid < textures.size());
	textures[tid].check_init();
}


void texture_manager::bind_texture(int tid) const {

	assert((unsigned)tid < textures.size());
	textures[tid].bind_gl();
}


colorRGBA texture_manager::get_tex_avg_color(int tid) const {

	assert((unsigned)tid < textures.size());
	return textures[tid].get_avg_color();
}


bool texture_manager::has_binary_alpha(int tid) const {

	assert((unsigned)tid < textures.size());
	return textures[tid].has_binary_alpha;
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
	
	free_vbos();
	vector<T>::clear();
	finalized = has_tangents = 0;
	bsphere.radius = 0.0;
}


template<typename T> void vntc_vect_t<T>::free_vbos() {

	delete_and_zero_vbo(vbo);
	delete_and_zero_vbo(ivbo);
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
	
	for (vector<unsigned>::const_iterator i = ixs.begin(); i != ixs.end(); ++i) {
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
				if (!bins[i].empty()) subdiv_recur(bins[i], npts, skip_dims);
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


template<typename T> void indexed_vntc_vect_t<T>::finalize(int prim_type) {

	if (need_normalize) {
		for (iterator i = begin(); i != end(); ++i) i->n.normalize();
		need_normalize = 0;
	}
	if (indices.empty() || finalized) return; // nothing to do
	finalized = 1;
	//optimize(prim_type); // now optimized in the object file loading phase
	unsigned const npts((prim_type == GL_TRIANGLES) ? 3 : 4), nverts(num_verts()); // triangles or quads
	assert((nverts % npts) == 0);
	assert(blocks.empty());

	if (nverts > 2*BLOCK_SIZE) { // subdivide large buffers
		//RESET_TIME;
		vector<unsigned> ixs;
		ixs.reserve(indices.size());
		ixs.swap(indices);
		subdiv_recur(ixs, npts, 0);
		//PRINT_TIME("Subdiv");
	}

#if 0
	vector<pair<unsigned, unsigned> > area_ix_pairs;
	area_ix_pairs.resize(nverts/npts);
	point pts[4];

	for (unsigned i = 0, ix = 0; i < nverts; i += npts, ++ix) {
		for (unsigned j = 0; j < npts; ++j) {pts[j] = get_vert(i+j).v;}
		//area_ix_pairs[ix] = make_pair(-polygon_area(pts, npts), i); // sort by triangle/quad size, largest to smallest
		area_ix_pairs[ix] = make_pair(indices[i], i); // sort by first vertex index
	}
	sort(area_ix_pairs.begin(), area_ix_pairs.end());
	vector<unsigned> new_indices;
	new_indices.resize(nverts);
	
	for (unsigned i = 0, ix = 0; i < nverts; i += npts, ++ix) {
		unsigned const vix(area_ix_pairs[ix].second);
		assert(vix+npts <= indices.size());
		for (unsigned j = 0; j < npts; ++j) {new_indices[i+j] = indices[vix+j];}
	}
	indices.swap(new_indices);
#endif
}


template<typename T> void indexed_vntc_vect_t<T>::clear() {
	
	vntc_vect_t<T>::clear();
	indices.clear();
	blocks.clear();
	need_normalize = 0;
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
		vector4d const tangent((v1*t2 - v2*t1).get_norm(), w);
		for (unsigned j = i; j < i+npts; ++j) {get_vert(j).tangent += tangent;}
	}
	if (!indices.empty()) { // using index array, need to renormalilze tangents
		for (iterator i = begin(); i != end(); ++i) {
			i->tangent.normalize();
			i->tangent.w = ((i->tangent.w < 0.0) ? -1.0 : 1.0);
		}
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


template<typename T> void indexed_vntc_vect_t<T>::render(shader_t &shader, bool is_shadow_pass, int prim_type, bool no_vfc) {

	if (empty()) return;
	finalize(prim_type);
	if (bsphere.radius == 0.0) calc_bounding_volumes();
	if (is_shadow_pass && vbo == 0) return; // don't create the vbo on the shadow pass (voxel terrain problems)

	if (no_vfc) {
		// do nothing
	}
	else if (is_shadow_pass) { // Note: makes shadow map caching more difficult
		if (no_sparse_smap_update() && !orig_camera_pdu.projected_cube_visible(bcube, camera_pdu.pos)) return; // light_pos == camera_pdu.pos for the shadow pass
	}
	else if (vbo) { // don't cull if vbo hasn't yet been allocated because this will cause it to be skipped in the shadow pass
		if (!camera_pdu.sphere_visible_test(bsphere.pos, bsphere.radius) || !camera_pdu.cube_visible(bcube)) return; // view frustum culling
	}
	unsigned const stride(sizeof(T));
	bool const have_normals(stride >= sizeof(vert_norm) && !is_shadow_pass), have_tex_coords(stride >= sizeof(vert_norm_tc) && !is_shadow_pass);
	set_array_client_state(1, have_tex_coords, have_normals, 0);
	int loc(-1);
	create_bind_vbo_and_upload(vbo, *this, 0);

	if (CALC_TANGENT_VECT && enable_bump_map() && !is_shadow_pass && has_tangents && shader.is_setup()) { // Note: if we get here, T must be a vert_norm_tc_tan
		assert(stride == sizeof(vert_norm_tc_tan));
		loc = shader.get_attrib_loc("tangent");
		assert(loc > 0);
		glEnableVertexAttribArray(loc);
		glVertexAttribPointer(loc, 4, GL_FLOAT, GL_FALSE, stride, (void *)sizeof(vert_norm_tc)); // stuff in at the end
	}
	glVertexPointer(3, GL_FLOAT, stride, 0);
	if (have_normals)    {glNormalPointer(     GL_FLOAT, stride, (void *)sizeof(point));}
	if (have_tex_coords) {glTexCoordPointer(2, GL_FLOAT, stride, (void *)sizeof(vert_norm));}

	if (indices.empty()) { // draw regular arrays
		glDrawArrays(prim_type, 0, (unsigned)size());
	}
	else { // draw indexed arrays
		create_bind_vbo_and_upload(ivbo, indices, 1);

		if (is_shadow_pass || blocks.empty() || no_vfc || camera_pdu.sphere_completely_visible_test(bsphere.pos, bsphere.radius)) { // draw the entire range
			glDrawRangeElements(prim_type, 0, (unsigned)size(), (unsigned)indices.size(), GL_UNSIGNED_INT, 0);
		}
		else { // draw each block independently
			for (vector<geom_block_t>::const_iterator i = blocks.begin(); i != blocks.end(); ++i) {
				if (camera_pdu.cube_visible(i->bcube)) {
					glDrawRangeElements(prim_type, 0, (unsigned)size(), i->num, GL_UNSIGNED_INT, (void *)(i->start_ix*sizeof(unsigned)));
				}
			}
		}
		bind_vbo(0, 1);
	}
	bind_vbo(0);
	if (loc >= 0) {glDisableVertexAttribArray(loc);}
}


template<typename T> void indexed_vntc_vect_t<T>::reserve_for_num_verts(unsigned num_verts) {

	if (!empty()) {
		// do nothing
	}
	else if (USE_INDEXED_VERTS) {
		indices.reserve(num_verts);
	}
	else {
		reserve(num_verts);
	}
}


template<typename T> void indexed_vntc_vect_t<T>::add_poly(polygon_t const &poly, vertex_map_t<T> &vmap) {

	for (unsigned i = 0; i < poly.size(); ++i) {add_vertex(poly[i], vmap);}
}


template<typename T> void indexed_vntc_vect_t<T>::add_triangle(triangle const &t, vertex_map_t<T> &vmap) {

	vector3d const normal(t.get_normal());
	//vector3d const normal(t.get_normal(!(USE_INDEXED_VERTS && vmap.get_average_normals()))); // weight by triangle area
	UNROLL_3X(add_vertex(T(t.pts[i_], normal), vmap);)
}


template<typename T> void indexed_vntc_vect_t<T>::add_vertex(T const &v, vertex_map_t<T> &vmap) {

	if (USE_INDEXED_VERTS) {
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
	}
	else {
		assert(vmap.empty());
		push_back(v);
	}
}


struct shared_vertex_t {
	unsigned ai, bi;
	bool shared;
	shared_vertex_t() : shared(0) {}
	shared_vertex_t(unsigned ai_, unsigned bi_) : ai(ai_), bi(bi_), shared(1) {}
};


template<typename T> void indexed_vntc_vect_t<T>::get_polygons(vector<coll_tquad> &polygons,
	colorRGBA const &color, unsigned npts, bool quads_only) const
{
	unsigned const nv(num_verts());
	assert((nv % npts) == 0);
	polygon_t poly(color), quad_poly(color);
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
					polygons.push_back(quad_poly);
					i += npts;
					continue;
				}
			}
		}
		for (unsigned p = 0; p < npts; ++p) {poly[p] = get_vert(i+p);}

		if (quads_only) {
			assert(npts == 4);
			polygons.push_back(poly);
		}
		else {
			split_polygon(poly, polygons, POLY_COPLANAR_THRESH);
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

	for (iterator i = begin(); i != end(); ++i) {i->free_vbos();}
}


template<typename T> cube_t vntc_vect_block_t<T>::get_bbox() const {

	if (empty()) return all_zeros_cube;
	cube_t bbox(front().get_bbox());
	for (const_iterator i = begin()+1; i != end(); ++i) {bbox.union_with_cube(i->get_bbox());}
	return bbox;
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


template<typename T> void vntc_vect_block_t<T>::get_polygons(vector<coll_tquad> &polygons,
	colorRGBA const &color, unsigned npts, bool quads_only) const
{
	for (const_iterator i = begin(); i != end(); ++i) {
		i->get_polygons(polygons, color, npts, quads_only);
	}
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


template<typename T> void geometry_t<T>::render_blocks(shader_t &shader, bool is_shadow_pass, vntc_vect_block_t<T> &blocks, int prim_type) {

	for (vntc_vect_block_t<T>::iterator i = blocks.begin(); i != blocks.end(); ++i) {
		i->render(shader, is_shadow_pass, prim_type);
	}
}


template<typename T> void geometry_t<T>::render(shader_t &shader, bool is_shadow_pass) {

	render_blocks(shader, is_shadow_pass, triangles, GL_TRIANGLES);
	render_blocks(shader, is_shadow_pass, quads,     GL_QUADS);
}


template<typename T> void geometry_t<T>::add_poly_to_polys(polygon_t const &poly, vntc_vect_block_t<T> &v, vertex_map_t<T> &vmap, unsigned obj_id) const {

	unsigned const max_entries(1 << 18); // 256K

	if (v.empty() || v.back().size() > max_entries || obj_id > v.back().obj_id) {
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


template<typename T> void geometry_t<T>::get_polygons(vector<coll_tquad> &polygons, colorRGBA const &color, bool quads_only) const {

	if (!quads_only) triangles.get_polygons(polygons, color, 3, 0);
	quads.get_polygons(polygons, color, 4, quads_only);
}


template<typename T> cube_t geometry_t<T>::get_bbox() const {

	cube_t bbox(all_zeros_cube); // will return this if empty

	if (!triangles.empty()) {
		bbox = triangles.get_bbox();
		if (!quads.empty()) bbox.union_with_cube(quads.get_bbox());
	}
	else if (!quads.empty()) {
		bbox = quads.get_bbox();
	}
	return bbox;
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


void material_t::render(shader_t &shader, texture_manager const &tmgr, int default_tid, bool is_shadow_pass) {

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
	if (is_shadow_pass) {
		geom.render(shader, 1);
		geom_tan.render(shader, 1);
	}
	else {
		int const tex_id(get_render_texture());
		bool has_binary_alpha(1);
		
		if (tex_id >= 0) {
			tmgr.bind_texture(tex_id);
			has_binary_alpha = tmgr.has_binary_alpha(tex_id);
		}
		else {
			select_texture(((default_tid >= 0) ? default_tid : WHITE_TEX), 0); // no texture specified - use white texture
		}
		if (use_bump_map()) {
			set_active_texture(5);
			tmgr.bind_texture(bump_tid);
			set_active_texture(0);
		}
		if (enable_spec_map()) { // all white/specular if no specular map texture
			set_active_texture(8);
			if (s_tid >= 0) {tmgr.bind_texture(s_tid);} else {select_texture(WHITE_TEX, 0);}
			set_active_texture(0);
		}
		if (alpha < 1.0 && ni != 1.0) {
			//shader.add_uniform_float("refract_index", ni); // FIXME: set index of refraction (and reset it at the end)
		}
		if (alpha_tid >= 0) enable_blend();
		float const spec_val((ks.R + ks.G + ks.B)/3.0);
		float const min_alpha((alpha_tid >= 0) ? (has_binary_alpha ? 0.9 : model3d_alpha_thresh) : 0.0);
		if (ns > 0.0) {set_specular(spec_val, ns);} // ns<=0 is undefined?
		set_color_e(colorRGBA(ke, alpha));

		if (shader.is_setup()) {
			set_color_d(get_ad_color());
			shader.add_uniform_float("min_alpha", min_alpha);
		}
		else {
			glAlphaFunc(GL_GREATER, min_alpha);
			set_color_a(colorRGBA((ignore_ambient ? kd : ka), alpha));
			set_color_d(colorRGBA(kd, alpha));
		}
		geom.render(shader, 0);
		geom_tan.render(shader, 0);
		set_color_e(BLACK);
		if (ns > 0.0) {set_specular(0.0, 1.0);}
		if (alpha_tid >= 0) disable_blend();
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
	if (!ignore_ambient) c += colorRGBA(ka, 0.0);
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

	if (CALC_TANGENT_VECT && use_bump_map()) {
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
	cube_t const bb(get_polygon_bbox(poly));
	if (bbox == all_zeros_cube) {bbox = bb;} else {bbox.union_with_cube(bb);}
	return (unsigned)split_polygons_buffer.size();
}


void model3d::get_polygons(vector<coll_tquad> &polygons, bool quads_only) const {

	if (polygons.empty()) { // Note: we count quads as 1.5 polygons because some of them may be split into triangles
		model3d_stats_t stats;
		get_stats(stats);
		polygons.reserve((quads_only ? stats.quads : (stats.tris + 1.5*stats.quads)));
	}
	colorRGBA unbound_color(WHITE);
	if (unbound_tid >= 0) unbound_color.modulate_with(texture_color(unbound_tid));
	unbound_geom.get_polygons(polygons, unbound_color, quads_only);

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		colorRGBA const color(m->get_avg_color(tmgr, unbound_tid));
		m->geom.get_polygons    (polygons, color, quads_only);
		m->geom_tan.get_polygons(polygons, color, quads_only);
	}
	::remove_excess_cap(polygons); // probably a good idea
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
	vector<vector<float_plus_dir> > &zvals, int bounds[2][2], int num_xy[2], float spacing, unsigned &nhq)
{
	model3d_stats_t stats;
	cont.get_stats(stats);
	polygons.resize(0);
	polygons.reserve(stats.quads);
	cont.get_polygons(polygons, WHITE, 1);
	
	for (vector<coll_tquad>::const_iterator i = polygons.begin(); i != polygons.end(); ++i) {
		assert(i->npts == 4);
		if (fabs(i->normal.z) < 0.99) continue; // only keep top/bottom cube sides
		cube_t const bcube(i->get_bcube());
		if ((bcube.d[2][1] - bcube.d[2][0]) > 0.5*spacing) continue; // can this happen?
		int cbounds[2][2];
		calc_bounds(bcube, cbounds, spacing);
		bool const is_top(i->normal.z > 0.0);
		++nhq;
		
		for (int y = cbounds[1][0]; y < cbounds[1][1]; ++y) {
			for (int x = cbounds[0][0]; x < cbounds[0][1]; ++x) {
				int const xv(x - bounds[0][0]), yv(y - bounds[1][0]);
				assert(xv >= 0 && yv >= 0 && xv < num_xy[0] && yv < num_xy[1]);
				zvals[xv + num_xy[0]*yv].push_back(float_plus_dir(bcube.d[2][0], is_top));
			}
		}
	}
	return (unsigned)polygons.size();
}


void model3d::get_cubes(vector<cube_t> &cubes, float spacing) const {

	assert(spacing > 0.0);

	// calculate scene voxel bounds
	int bounds[2][2], num_xy[2]; // {x,y}x{lo,hi}
	calc_bounds(bbox, bounds, spacing);
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
		num_polys += add_polygons_to_voxel_grid(polygons, unbound_geom, zvals, bounds, num_xy, spacing, num_horiz_quads);

		for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
			num_polys += add_polygons_to_voxel_grid(polygons, m->geom,     zvals, bounds, num_xy, spacing, num_horiz_quads);
			num_polys += add_polygons_to_voxel_grid(polygons, m->geom_tan, zvals, bounds, num_xy, spacing, num_horiz_quads);
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
					if (!merged) cubes.push_back(cube);
					cube.d[2][0] = val; // next segment starts here in case we get two top edges in a row
					++num_pre_merged_cubes;
				}
				//in_cube ^= 1;
			}
			//assert(!in_cube);
		}
	}
	cout << "polygons: " << num_polys << ", hquads: " << num_horiz_quads << ", pre_merged_cubes: " << num_pre_merged_cubes << ", cubes: " << cubes.size() << endl;
}


int model3d::get_material_ix(string const &material_name, string const &fn) {

	unsigned mat_id(0);
	string_map_t::const_iterator it(mat_map.find(material_name));

	if (it == mat_map.end()) {
		mat_id = (unsigned)materials.size();
		mat_map[material_name] = mat_id;
		materials.push_back(material_t(material_name, fn, ignore_ambient));
	}
	else {
		if (!from_model3d_file) cerr << "Warning: Redefinition of material " << material_name << " in file " << fn << endl;
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

	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		m->optimize();
	}
	unbound_geom.optimize();
}


void model3d::clear() {

	free_context();
	unbound_geom.clear();
	materials.clear();
	undef_materials.clear();
	mat_map.clear();
	coll_tree.clear();
}


void model3d::free_context() {

	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		m->geom.free_vbos();
		m->geom_tan.free_vbos();
	}
	unbound_geom.free_vbos();
}


void model3d::load_all_used_tids() {

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->mat_is_used()) continue;
		int const tid(m->get_render_texture());
		tmgr.bind_alpha_channel_to_texture(tid, m->alpha_tid);
		tmgr.ensure_tid_loaded(tid, 0); // only one tid for now
		if (m->use_bump_map()) tmgr.ensure_tid_loaded(m->bump_tid, 1);
		if (m->use_spec_map()) tmgr.ensure_tid_loaded(m->s_tid, 0);
	}
}


void model3d::bind_all_used_tids() {

	load_all_used_tids();
		
	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->mat_is_used()) continue;
		tmgr.ensure_tid_bound(m->get_render_texture()); // only one tid for now
		
		if (m->use_bump_map()) {
			if (CALC_TANGENT_VECT && !m->geom.empty()) {
				cerr << "Error loading model3d material " << m->name << ": Geometry is missing tangent vectors, so bump map cannot be enabled." << endl;
				m->bump_tid = -1; // disable bump map
			}
			else {
				tmgr.ensure_tid_bound(m->bump_tid);
				m->geom_tan.calc_tangents();
			}
			needs_bump_maps = 1;
		}
		if (m->use_spec_map()) tmgr.ensure_tid_bound(m->s_tid);
		needs_alpha_test |= m->get_needs_alpha_test();
	}
}


void model3d::render(shader_t &shader, bool is_shadow_pass, unsigned bmap_pass_mask) { // const?

	// we need the vbo to be created here even in the shadow pass,
	// and the textures are needed for determining whether or not we need to build the tanget_vectors for bump mapping
	bind_all_used_tids();
	if (group_back_face_cull || is_shadow_pass) glEnable(GL_CULL_FACE);

	// render geom that was not bound to a material
	if ((bmap_pass_mask & 1) && unbound_color.alpha > 0.0) { // enabled, not in bump map only pass
		if (!is_shadow_pass) {
			assert(unbound_tid >= 0);
			select_texture(unbound_tid, 0);
			set_color_d(unbound_color);
			if (shader.is_setup()) shader.add_uniform_float("min_alpha", 0.0);
		}
		unbound_geom.render(shader, is_shadow_pass);
	}
	
	// render all materials (opaque then transparent)
	for (unsigned pass = 0; pass < 2; ++pass) { // opaque, transparent
		vector<pair<float, unsigned> > to_draw;

		for (unsigned i = 0; i < materials.size(); ++i) {
			if (materials[i].is_partial_transparent() == (pass != 0) && (bmap_pass_mask & (1 << unsigned(materials[i].use_bump_map())))) {
				to_draw.push_back(make_pair(materials[i].draw_order_score, i));
			}
		}
		sort(to_draw.begin(), to_draw.end());

		for (unsigned i = 0; i < to_draw.size(); ++i) {
			materials[to_draw[i].second].render(shader, tmgr, unbound_tid, is_shadow_pass);
		}
	}
	if (group_back_face_cull || is_shadow_pass) glDisable(GL_CULL_FACE);
}


void model3d::build_cobj_tree(bool verbose) {

	if (!coll_tree.is_empty() || has_cobjs) return; // already built or not needed because cobjs will be used instead
	RESET_TIME;
	get_polygons(coll_tree.get_tquads_ref());
	coll_tree.build_tree_top(verbose);
	PRINT_TIME(" Cobj Tree Create (from model3d)");
}


void model3d::get_all_mat_lib_fns(set<string> &mat_lib_fns) const {

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		mat_lib_fns.insert(m->filename);
	}
}


void model3d::get_stats(model3d_stats_t &stats) const {

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


bool model3d::write_to_disk(string const &fn) const {

	ofstream out(fn, ios::out | ios::binary);
	
	if (!out.good()) {
		cerr << "Error opening model3d file for write: " << fn << endl;
		return 0;
	}
	cout << "Writing model3d file " << fn << endl;
	write_uint(out, MAGIC_NUMBER);
	out.write((char const *)&bbox, sizeof(cube_t));
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


bool model3d::read_from_disk(string const &fn) {

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
	in.read((char *)&bbox, sizeof(cube_t));
	if (!unbound_geom.read(in)) return 0;
	materials.resize(read_uint(in));
	
	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->read(in)) {
			cerr << "Error reading material" << endl;
			return 0;
		}
		mat_map[m->name]  = (m - materials.begin());
		m->ignore_ambient = ignore_ambient;
	}
	return in.good();
}


// ************ model3ds ************

void model3ds::clear() {

	for (iterator m = begin(); m != end(); ++m) {
		m->clear();
	}
	deque<model3d>::clear();
	tmgr.clear();
}


void model3ds::free_context() {

	for (iterator m = begin(); m != end(); ++m) {
		m->free_context();
	}
	tmgr.free_tids();
}


void model3ds::render(bool is_shadow_pass) {
	
	if (empty()) return;
	bool const shader_effects(!disable_shaders && !is_shadow_pass);
	set_lighted_sides(2);
	set_fill_mode();
	
	if (is_shadow_pass) {
		// FIXME: textuing is disabled, so alpha mask textures won't work
		// FIXME SHADERS: uses fixed function pipeline
	}
	else if (shader_effects && have_indir_smoke_tex) {
		set_color_a(BLACK); // ambient will be set by indirect lighting in the shader
	}
	else {
		set_color_a(WHITE);
	}
	BLACK.do_glColor();
	set_specular(0.0, 1.0);
	bool needs_alpha_test(0), needs_bump_maps(0);

	for (const_iterator m = begin(); m != end(); ++m) {
		needs_alpha_test |= m->get_needs_alpha_test();
		if (shader_effects) {needs_bump_maps |= m->get_needs_bump_maps();} // optimization, makes little difference
	}
	float const min_alpha(needs_alpha_test ? 0.5 : 0.0); // will be reset per-material, but this variable is used to enable alpha testing

	for (unsigned bmap_pass = 0; bmap_pass < (needs_bump_maps ? 2U : 1U); ++bmap_pass) {
		shader_t s;

		if (is_shadow_pass) {

		}
		else if (shader_effects) {
			int const use_bmap((bmap_pass == 0) ? 0 : (CALC_TANGENT_VECT ? 2 : 1));
			setup_smoke_shaders(s, min_alpha, 0, 0, 1, 1, 1, 1, 0, 1, use_bmap, enable_spec_map(), 0, two_sided_lighting);
		}
		else {
			s.setup_enabled_lights(2, 1); // sun and moon VS lighting
			s.set_vert_shader("ads_lighting.part*+two_lights_texture");
			s.set_frag_shader("simple_texture");
			s.begin_shader();
			s.add_uniform_int("tex0", 0);
		}
		for (iterator m = begin(); m != end(); ++m) { // non-const
			m->render(s, is_shadow_pass, (shader_effects ? (1 << bmap_pass) : 3));
		}
		if (!is_shadow_pass) {s.end_shader();}
	}
	set_lighted_sides(1);
	set_specular(0.0, 1.0);
}


cube_t model3ds::get_bbox() const {

	cube_t bbox(all_zeros_cube); // will return this if empty()

	for (const_iterator m = begin(); m != end(); ++m) {
		cube_t const &bb(m->get_bbox());
		if (m == begin()) {bbox = bb;} else {bbox.union_with_cube(bb);}
	}
	return bbox;
}


void model3ds::build_cobj_trees(bool verbose) {

	for (iterator m = begin(); m != end(); ++m) {
		m->build_cobj_tree(verbose);
	}
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


// ************ Free Functions ************


void free_model_context() {
	all_models.free_context();
}

void render_models(bool shadow_pass) {
	all_models.render(shadow_pass);
}


