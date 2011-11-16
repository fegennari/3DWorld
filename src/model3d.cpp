// 3D World - 3D Model Rendering Code
// by Frank Gennari
// 8/17/11

#include "model3d.h"
#include "shaders.h"
#include "gl_ext_arb.h"
#include <fstream>

bool const USE_SHADERS       = 1;
bool const ENABLE_BUMP_MAPS  = 1;
bool const ENABLE_SPEC_MAPS  = 1;
bool const USE_INDEXED_VERTS = 1;
unsigned const MAGIC_NUMBER  = 42987143; // arbitrary file signature

extern bool group_back_face_cull, enable_model3d_tex_comp, disable_shaders;
extern int display_mode;


model3ds all_models;


bool get_use_shaders() {return (USE_SHADERS && !disable_shaders);}
bool enable_bump_map() {return (ENABLE_BUMP_MAPS && get_use_shaders() && (display_mode & 0x20) == 0);} // enabled by default
bool enable_spec_map() {return (ENABLE_SPEC_MAPS && get_use_shaders());}


// ************ texture_manager ************

unsigned texture_manager::create_texture(string const &fn, bool is_alpha_mask, bool verbose) {

	string_map_t::const_iterator it(tex_map.find(fn));

	if (it != tex_map.end()) { // found (already loaded)
		assert(it->second < textures.size());
		return it->second;
	}
	unsigned const tid(textures.size());
	tex_map[fn] = tid;
	if (verbose) cout << "creating texture " << fn << endl;
	bool const compress(!is_alpha_mask && enable_model3d_tex_comp);
	// type=read_from_file format=auto width height wrap ncolors use_mipmaps name [do_compress]
	textures.push_back(texture_t(0, 7, 0, 0, 1, (is_alpha_mask ? 1 : 3), !is_alpha_mask, fn, compress)); // always RGB targa wrapped+mipmap
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
		t->free();
	}
}


void texture_manager::ensure_texture_loaded(texture_t &t, int tid, bool is_bump) {

	if (!t.is_allocated()) {
		t.load(-1);
		
		if (t.alpha_tid >= 0 && t.alpha_tid != tid) { // if alpha is the same texture then the alpha channel should already be set
			ensure_tid_loaded(t.alpha_tid, 0);
			assert((unsigned)t.alpha_tid < textures.size());
			t.copy_alpha_from_texture(textures[t.alpha_tid]);
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
	write_uint(out, v.size());
	out.write((const char *)&v.front(), v.size()*sizeof(V::value_type));
}

template<typename V> void read_vector(istream &in, V &v) {
	v.clear();
	v.resize(read_uint(in));
	in.read((char *)&v.front(), v.size()*sizeof(V::value_type));
}


// ************ vntc_vect_t/indexed_vntc_vect_t ************

template<typename T> void vntc_vect_t<T>::free_vbos() {

	delete_vbo(vbo);
	delete_vbo(ivbo);
	vbo = ivbo = 0;
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
		for (unsigned j = i; j < i+npts; ++j) get_vert(j).tangent += tangent;
	}
	if (!indices.empty()) { // using index array, need to renormalilze tangents
		for (iterator i = begin(); i != end(); ++i) {
			i->tangent.normalize();
			i->tangent.w = ((i->tangent.w < 0.0) ? -1.0 : 1.0); // FIXME: what if 0.0?
		}
	}
}


template<typename T> void vntc_vect_t<T>::write(ostream &out) const {

	write_vector(out, *this);
}


template<typename T> void vntc_vect_t<T>::read(istream &in) {

	read_vector(in, *this);
	has_tangents = (sizeof(T) == sizeof(vert_norm_tc_tan)); // HACK to get the type
	calc_bounding_volumes();
}


template<typename T> void indexed_vntc_vect_t<T>::render(shader_t &shader, bool is_shadow_pass, int prim_type) {

	if (empty()) return;
	if (bsphere.radius == 0.0) calc_bounding_volumes();
	if (!is_shadow_pass && !camera_pdu.sphere_visible_test(bsphere.pos, bsphere.radius)) return; // view frustum culling
	if (!is_shadow_pass && !camera_pdu.cube_visible(bcube)) return; // test the bounding cube as well
	set_array_client_state(1, !is_shadow_pass, !is_shadow_pass, 0);
	unsigned const stride(sizeof(T));
	int loc(-1);

	if (vbo == 0) {
		vbo = create_vbo();
		assert(vbo > 0);
		bind_vbo(vbo);
		upload_vbo_data(&front(), size()*stride);
	}
	else {
		bind_vbo(vbo);
	}
	if (enable_bump_map() && !is_shadow_pass && has_tangents) { // Note: if we get here, T must be a vert_norm_tc_tan
		assert(stride == sizeof(vert_norm_tc_tan));
		int const loc(shader.get_attrib_loc("tangent"));
		assert(loc > 0);
		glEnableVertexAttribArray(loc);
		glVertexAttribPointer(loc, 4, GL_FLOAT, GL_FALSE, stride, (void *)sizeof(vert_norm_tc)); // stuff in at the end
	}
	glVertexPointer(  3, GL_FLOAT, stride, 0);
	glNormalPointer(     GL_FLOAT, stride, (void *)sizeof(point));
	glTexCoordPointer(2, GL_FLOAT, stride, (void *)sizeof(vert_norm));

	if (indices.empty()) { // draw regular arrays
		glDrawArrays(prim_type, 0, size());
	}
	else { // draw indexed arrays
		if (ivbo == 0) {
			ivbo = create_vbo();
			assert(ivbo > 0);
			bind_vbo(ivbo, 1);
			upload_vbo_data(&indices.front(), indices.size()*sizeof(unsigned), 1);
		}
		else {
			bind_vbo(ivbo, 1);
		}
		glDrawRangeElements(prim_type, 0, size(), indices.size(), GL_UNSIGNED_INT, 0);
		bind_vbo(0, 1);
	}
	bind_vbo(0);
	if (loc >= 0) glDisableVertexAttribArray(loc);
}


template<typename T> void indexed_vntc_vect_t<T>::add_poly(polygon_t const &poly, vertex_map_t<T> &vmap) {

	for (unsigned i = 0; i < poly.size(); ++i) {add_vertex(poly[i], vmap);}
}


template<typename T> void indexed_vntc_vect_t<T>::add_vertex(T const &v, vertex_map_t<T> &vmap) {

	if (USE_INDEXED_VERTS) {
		vertex_map_t<T>::const_iterator it(vmap.find(v));
		unsigned ix;

		if (it == vmap.end()) { // not found
			ix = size();
			push_back(v);
			vmap[v] = ix;
		}
		else { // found
			ix = it->second;
		}
		assert(ix < size());
		indices.push_back(ix);
	}
	else {
		assert(vmap.empty());
		push_back(v);
	}
}


template<typename T> void indexed_vntc_vect_t<T>::get_polygons(vector<polygon_t> &polygons, colorRGBA const &color, unsigned npts) const {

	unsigned const nv(num_verts());
	assert((nv % npts) == 0);

	for (unsigned i = 0; i < nv; i += npts) {
		polygon_t poly(color);
		poly.resize(npts);
		for (unsigned p = 0; p < npts; ++p) {poly[p] = get_vert(i+p);}
		split_polygon(poly, polygons, POLY_COPLANAR_THRESH);
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

bool polygon_t::is_convex() const {

	unsigned const npts(size());
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
	unsigned counts[2] = {0};
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

template<typename T> void vntc_vect_block_t<T>::remove_excess_cap() {

	for (iterator i = begin(); i != end(); ++i) {i->remove_excess_cap();}
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
	for (const_iterator i = begin(); i != end(); ++i) {s += i->size();}
	return s;
}


template<typename T> void vntc_vect_block_t<T>::get_polygons(vector<polygon_t> &polygons, colorRGBA const &color, unsigned npts) const {

	for (const_iterator i = begin(); i != end(); ++i) {
		i->get_polygons(polygons, color, npts);
	}
}


template<typename T> bool vntc_vect_block_t<T>::write(ostream &out) const {

	write_uint(out, size());
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

template<> void geometry_t<vert_norm_tc_tan>::calc_tangents() {

	for (vntc_vect_block_t<vert_norm_tc_tan>::iterator i = triangles.begin(); i != triangles.end(); ++i) {
		i->calc_tangents(3);
	}
	for (vntc_vect_block_t<vert_norm_tc_tan>::iterator i = quads.begin(); i != quads.end(); ++i) {
		i->calc_tangents(4);
	}
}


template<typename T> void geometry_t<T>::render(shader_t &shader, bool is_shadow_pass) {

	for (vntc_vect_block_t<T>::iterator i = triangles.begin(); i != triangles.end(); ++i) {
		i->render(shader, is_shadow_pass, GL_TRIANGLES);
	}
	for (vntc_vect_block_t<T>::iterator i = quads.begin(); i != quads.end(); ++i) {
		i->render(shader, is_shadow_pass, GL_QUADS);
	}
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


template<typename T> void geometry_t<T>::get_polygons(vector<polygon_t> &polygons, colorRGBA const &color) const {

	triangles.get_polygons(polygons, color, 3);
	quads.get_polygons    (polygons, color, 4);
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


void material_t::render(shader_t &shader, texture_manager const &tmgr, int default_tid, bool is_shadow_pass) {

	if ((geom.empty() && geom_tan.empty()) || skip || alpha == 0.0) return; // empty or transparent
	if (is_shadow_pass && alpha < MIN_SHADOW_ALPHA) return;

	if (is_shadow_pass) {
		geom.render(shader, 1);
		geom_tan.render(shader, 1);
	}
	else {
		int const tex_id(get_render_texture());
		
		if (tex_id >= 0) {
			tmgr.bind_texture(tex_id);
		}
		else {
			select_texture(((default_tid >= 0) ? default_tid : WHITE_TEX), 0); // no texture specified - use white texture
		}
		if (use_bump_map()) {
			set_multitex(5);
			tmgr.bind_texture(bump_tid);
		}
		if (enable_spec_map()) { // all white/specular if no specular map texture
			set_multitex(8);
			if (s_tid >= 0) tmgr.bind_texture(s_tid); else select_texture(WHITE_TEX);
		}
		if (alpha < 1.0 && ni != 1.0) {
			// set index of refraction (and reset it at the end)
		}
		if (alpha_tid >= 0) enable_blend();
		float const spec_val((ks.R + ks.G + ks.B)/3.0);
		set_specular(spec_val, ns);
		set_color_e(colorRGBA(ke, alpha));

		if (shader.is_setup()) {
			set_color_d(get_ad_color());
			shader.add_uniform_float("min_alpha", ((alpha_tid >= 0) ? 0.9 : 0.0)); // FIXME: check has_binary_alpha?
		}
		else {
			set_color_a(colorRGBA((ignore_ambient ? kd : ka), alpha));
			set_color_d(colorRGBA(kd, alpha));
		}
		geom.render(shader, 0);
		geom_tan.render(shader, 0);
		if (use_bump_map())    disable_multitex(5, 1);
		if (enable_spec_map()) disable_multitex(8, 1);
		set_color_e(BLACK);
		set_specular(0.0, 1.0);
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

	if (use_bump_map()) {
		geom_tan.add_poly(poly, vmap_tan, obj_id);
	}
	else {
		geom.add_poly(poly, vmap, obj_id);
	}
	mark_as_used();
	return 1;
}


bool material_t::write(ostream &out) const {

	out.write((char const *)this, sizeof(material_params_t)); // FIXME: can write garbage bytes
	write_vector(out, name);
	write_vector(out, filename);
	return (geom.write(out) && geom_tan.write(out));
}


bool material_t::read(istream &in) {

	in.read((char *)this, sizeof(material_params_t)); // FIXME: can read garbage bytes
	read_vector(in, name);
	read_vector(in, filename);
	return (geom.read(in) && geom_tan.read(in));
}


// ************ model3d ************


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
	return split_polygons_buffer.size();
}


void model3d::get_polygons(vector<polygon_t> &polygons) const {

	if (polygons.empty()) {
		model3d_stats_t stats;
		get_stats(stats);
		polygons.reserve(stats.tris + 1.5*stats.quads); // Note: we count quads as 1.5 polygons because some of them may be split into triangles
	}
	colorRGBA unbound_color(def_color);
	if (unbound_tid >= 0) unbound_color.modulate_with(texture_color(unbound_tid));
	unbound_geom.get_polygons(polygons, unbound_color);

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		colorRGBA const color(m->get_avg_color(tmgr, unbound_tid));
		m->geom.get_polygons(polygons, color);
		m->geom_tan.get_polygons(polygons, color);
	}
	::remove_excess_cap(polygons); // probably a good idea
}


int model3d::get_material_ix(string const &material_name, string const &fn) {

	unsigned mat_id(0);
	string_map_t::const_iterator it(mat_map.find(material_name));

	if (it == mat_map.end()) {
		mat_id = materials.size();
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


void model3d::remove_excess_cap() {

	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		m->geom.remove_excess_cap();
		m->geom_tan.remove_excess_cap();
	}
	unbound_geom.remove_excess_cap();
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
			if (!m->geom.empty()) {
				cerr << "Error loading model3d material " << m->name << ": Geometry is missing tangent vectors, so bump map cannot be enabled." << endl;
				m->bump_tid = -1; // disable bump map
			}
			else {
				tmgr.ensure_tid_bound(m->bump_tid);
				m->geom_tan.calc_tangents();
			}
		}
		if (m->use_spec_map()) tmgr.ensure_tid_bound(m->s_tid);
	}
}


void model3d::render(shader_t &shader, bool is_shadow_pass, bool bmap_pass) { // const?

	// we need the vbo to be created here even in the shadow pass,
	// and the textures are needed for determining whether or not we need to build the tanget_vectors for bump mapping
	bind_all_used_tids();
	bool const do_cull(group_back_face_cull && !is_shadow_pass);
	if (do_cull) glEnable(GL_CULL_FACE);

	// render geom that was not bound to a material
	if (!bmap_pass && unbound_color.alpha > 0.0) { // enabled, not in bump map pass
		if (!is_shadow_pass) {
			assert(unbound_tid >= 0);
			select_texture(unbound_tid, 0);
			set_color_d(unbound_color);
			if (shader.is_setup()) shader.add_uniform_float("min_alpha", 0.0);
		}
		unbound_geom.render(shader, is_shadow_pass);
	}
	
	// render all materials (opaque then transparen)
	for (unsigned pass = 0; pass < 2; ++pass) { // opaque, transparent
		for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
			if (m->is_partial_transparent() == (pass != 0) && m->use_bump_map() == bmap_pass) {
				m->render(shader, tmgr, unbound_tid, is_shadow_pass);
			}
		}
	}
	if (do_cull) glDisable(GL_CULL_FACE);
}


void model3d::build_cobj_tree(bool verbose) {

	if (!coll_tree.is_empty() || has_cobjs) return; // already built or not needed because cobjs will be used instead
	vector<polygon_t> polygons;
	get_polygons(polygons);
	coll_tree.add_polygons(polygons, verbose);
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
	write_uint(out, materials.size());
	
	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->write(out)) {
			cerr << "Error writing material" << endl;
			return 0;
		}
	}
	return out.good();
}


bool model3d::read_from_disk(string const &fn) {

	cout << "fn: " << fn << endl;
	ifstream in(fn, ios::in | ios::binary);
	
	if (!in.good()) {
		cerr << "Error opening model3d file for read: " << fn << endl;
		return 0;
	}
	clear(); // ???
	unsigned const magic_number_comp(read_uint(in));
	cout << "number: " << magic_number_comp << endl;

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

	bool const use_shaders(get_use_shaders() && !is_shadow_pass);
	set_lighted_sides(2);
	set_fill_mode();
	
	if (is_shadow_pass) {
		glDisable(GL_LIGHTING);
	}
	else if (use_shaders) {
		glDisable(GL_LIGHTING); // custom lighting calculations from this point on
		set_color_a(BLACK); // ambient will be set by indirect lighting in the shader
	}
	else {
		set_color_a(WHITE); // ambient will be set by indirect lighting in the shader
		glEnable(GL_TEXTURE_2D);
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.5);
	}
	BLACK.do_glColor();
	set_specular(0.0, 1.0);
	float const min_alpha(0.5); // will be reset per-material

	for (unsigned bmap_pass = 0; bmap_pass < 2; ++bmap_pass) {
		shader_t s;
		colorRGBA orig_fog_color;
		
		if (use_shaders) {
			orig_fog_color = setup_smoke_shaders(s, min_alpha, 0, 0, 1, 1, 1, 1, 0, shadow_map_enabled(), (bmap_pass != 0), enable_spec_map());
		}
		bool const render_if_bmap(0); // set this later when bump maps are supported

		for (iterator m = begin(); m != end(); ++m) {
			m->render(s, is_shadow_pass, (bmap_pass != 0));
		}
		if (use_shaders) end_smoke_shaders(s, orig_fog_color);
	}
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
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


void free_model_context() {
	all_models.free_context();
}

void render_models(bool shadow_pass) {
	all_models.render(shadow_pass);
}


