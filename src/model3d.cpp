// 3D World - 3D Model Rendering Code
// by Frank Gennari
// 8/17/11

#include "model3d.h"
#include "shaders.h"
#include "gl_ext_arb.h"

extern bool group_back_face_cull, enable_model3d_tex_comp;


model3ds all_models;


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
	// type format width height wrap ncolors use_mipmaps name [bump_name [id [color]]]
	textures.push_back(texture_t(0, 4, 0, 0, 1, (is_alpha_mask ? 1 : 3), !is_alpha_mask, fn)); // always RGB targa wrapped+mipmap
	textures.back().do_compress = (!is_alpha_mask && enable_model3d_tex_comp);
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


void texture_manager::ensure_texture_loaded(texture_t &t) {

	if (!t.data) {
		t.load(-1);
		
		if (t.alpha_tid >= 0) {
			ensure_tid_loaded(t.alpha_tid);
			assert((unsigned)t.alpha_tid < textures.size());
			texture_t &at(textures[t.alpha_tid]);
			assert(at.data && t.data); // check that data is allocated in both textures
			assert(at.width == t.width && at.height == t.height); // check for matching sizes
			assert(t.ncolors == 4); // check for alpha channel
			assert(t.tid == 0); // check that texture isn't already bound
			assert(at.ncolors == 1);
			unsigned const npixels(t.num_pixels());

			for (unsigned i = 0; i < npixels; ++i) {
				t.data[4*i+3] = at.data[i]; // copy alpha values
			}
		}
		t.init(); // must be after alpha copy
	}
	assert(t.data);
}


void texture_manager::bind_alpha_channel_to_texture(int tid, int alpha_tid) {

	if (tid < 0 || alpha_tid < 0) return; // no texture
	assert((unsigned)tid < textures.size() && (unsigned)alpha_tid < textures.size());
	texture_t &t(textures[tid]);
	assert(t.ncolors == 3 || t.ncolors == 4);
	if (t.alpha_tid == alpha_tid) return; // already bound
	assert(t.alpha_tid < 0); // can't rebind to a different value
	assert(!t.data); // must not yet be loaded
	t.alpha_tid = alpha_tid;
	t.ncolors   = 4; // add alpha channel
	if (t.use_mipmaps) t.use_mipmaps = 3; // generate custom alpha mipmaps
}


void texture_manager::ensure_tid_loaded(int tid) {

	if (tid < 0) return; // not allocated
	assert((unsigned)tid < textures.size());
	ensure_texture_loaded(textures[tid]);
}


void texture_manager::ensure_tid_bound(int tid) {

	if (tid < 0) return; // not allocated
	assert((unsigned)tid < textures.size());
	textures[tid].check_init();
}


void texture_manager::bind_texture(int tid) const {

	//select_texture(WHITE_TEX, 0); return; // TESTING
	assert((unsigned)tid < textures.size());
	assert(textures[tid].tid > 0);
	glBindTexture(GL_TEXTURE_2D, textures[tid].tid);
}


colorRGBA texture_manager::get_tex_avg_color(int tid) const {

	assert((unsigned)tid < textures.size());
	return textures[tid].color;
}


// ************ vntc_vect_t ************

void vntc_vect_t::render(bool is_shadow_pass) const {

	assert(size() >= 3);

	for (const_iterator v = begin(); v != end(); ++v) {
		if (!is_shadow_pass) {
			glTexCoord2fv(v->t);
			v->n.do_glNormal();
		}
		v->v.do_glVertex();
	}
}


void vntc_vect_t::render_array(bool is_shadow_pass, int prim_type) {

	if (empty()) return;
	set_array_client_state(1, !is_shadow_pass, !is_shadow_pass, 0);

	if (vbo == 0) {
		vbo = create_vbo();
		assert(vbo > 0);
		bind_vbo(vbo);
		upload_vbo_data(&front(), size()*sizeof(vert_norm_tc));
	}
	else {
		bind_vbo(vbo);
	}
	glVertexPointer(  3, GL_FLOAT, sizeof(vert_norm_tc), 0);
	glNormalPointer(     GL_FLOAT, sizeof(vert_norm_tc), (void *)sizeof(point));
	glTexCoordPointer(2, GL_FLOAT, sizeof(vert_norm_tc), (void *)sizeof(vert_norm));
	glDrawArrays(prim_type, 0, size());
	bind_vbo(0);
}


void vntc_vect_t::free_vbo() {

	delete_vbo(vbo);
	vbo = 0;
}


bool vntc_vect_t::is_convex() const {

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


vector3d vntc_vect_t::get_planar_normal() const {

	assert(size() >= 3);
	vector3d norm;
	get_normal((*this)[0].v, (*this)[1].v, (*this)[2].v, norm, 1);
	return norm;
}


void vntc_vect_t::from_points(vector<point> const &pts) {

	resize(pts.size());
	for (unsigned i = 0; i < size(); ++i) (*this)[i].v = pts[i];
}


// ************ geometry_t ************

void geometry_t::render(bool is_shadow_pass) {

	triangles.render_array(is_shadow_pass, GL_TRIANGLES);
	quads.render_array    (is_shadow_pass, GL_QUADS);
}


void geometry_t::add_poly(vntc_vect_t const &poly) {
	
	if (poly.size() == 3) { // triangle
		for (unsigned i = 0; i < 3; ++i) {triangles.push_back(poly[i]);}
		return;
	}
	if (poly.size() == 4) {
		for (unsigned i = 0; i < 4; ++i) {quads.push_back(poly[i]);}
		return;
	}
	assert(0); // shouldn't get here
}


void geometry_t::free_vbos() {

	triangles.free_vbo();
	quads.free_vbo();
}


void geometry_t::clear() {

	free_vbos();
	triangles.clear();
	quads.clear();
}


// ************ material_t ************

void material_t::render(texture_manager const &tmgr, int default_tid, bool is_shadow_pass) {

	if (geom.empty() || skip || alpha == 0.0)       return; // empty or transparent
	if (is_shadow_pass && alpha < MIN_SHADOW_ALPHA) return;

	if (!is_shadow_pass) {
		int const tex_id(get_render_texture());
		
		if (tex_id >= 0) {
			tmgr.bind_texture(tex_id);
		}
		else {
			select_texture(((default_tid >= 0) ? default_tid : WHITE_TEX), 0); // no texture specified - use white texture
		}
		if (alpha < 1.0 && ni != 1.0) {
			// set index of refraction (and reset it at the end)
		}
		if (alpha_tid >= 0) enable_blend();
		float const spec_val((ks.R + ks.G + ks.B)/3.0);
		set_specular(spec_val, ns);
		//set_color_a(colorRGBA(ka, alpha));
		//set_color_d(colorRGBA(kd, alpha));
		set_color_d(get_ad_color());
		set_color_e(colorRGBA(ke, alpha));
	}
	geom.render(is_shadow_pass);

	if (!is_shadow_pass) {
		set_color_e(BLACK);
		set_specular(0.0, 1.0);
		if (alpha_tid >= 0) disable_blend();
	}
}


colorRGBA material_t::get_ad_color() const {

	colorRGBA c(colorRGBA(kd, alpha) + colorRGBA(ka, 0.0));
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


bool material_t::add_poly(vntc_vect_t const &poly) {
	
	if (skip) return 0;
	geom.add_poly(poly);
	mark_as_used();
	return 1;
}


// ************ model3d ************

void model3d::add_polygon(vntc_vect_t const &poly, int mat_id, vector<polygon_t> *ppts) {

	//assert(mat_id >= 0); // must be set/valid - too strict?
	split_polygons_buffer.resize(0);
	split_polygon(poly, split_polygons_buffer);

	for (vector<polygon_t>::const_iterator i = split_polygons_buffer.begin(); i != split_polygons_buffer.end(); ++i) {
		if (mat_id < 0) {
			unbound_geom.add_poly(*i);
			if (ppts) ppts->push_back(*i);
		}
		else {
			assert((unsigned)mat_id < materials.size());
		
			if (materials[mat_id].add_poly(*i)) {
				if (ppts) {
					ppts->push_back(*i);
					ppts->back().color = materials[mat_id].get_avg_color(tmgr, unbound_tid);
				}
			}
		}
	}
}


int model3d::get_material_ix(string const &material_name, string const &fn) {

	unsigned mat_id(0);
	string_map_t::const_iterator it(mat_map.find(material_name));

	if (it == mat_map.end()) {
		mat_id = materials.size();
		mat_map[material_name] = mat_id;
		materials.push_back(material_t());
	}
	else {
		cerr << "Warning: Redefinition of material " << material_name << " in file " << fn << endl;
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


void model3d::clear() {

	free_context();
	unbound_geom.clear();
	materials.clear();
	undef_materials.clear();
	mat_map.clear();
}


void model3d::free_context() {

	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		m->geom.free_vbos();
	}
	unbound_geom.free_vbos();
}


void model3d::load_all_used_tids() {

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->mat_is_used()) continue;
		int const tid(m->get_render_texture());
		tmgr.bind_alpha_channel_to_texture(tid, m->alpha_tid);
		tmgr.ensure_tid_loaded(tid); // only one tid for now
	}
}


void model3d::bind_all_used_tids() {

	load_all_used_tids();
		
	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (!m->mat_is_used()) continue;
		tmgr.ensure_tid_bound(m->get_render_texture()); // only one tid for now
	}
}


void model3d::render(bool is_shadow_pass) { // const?

	if (!is_shadow_pass) bind_all_used_tids();
	bool const do_cull(group_back_face_cull && !is_shadow_pass);
	if (do_cull) glEnable(GL_CULL_FACE);

	// render geom that was not bound to a material
	if (unbound_color.alpha > 0.0) { // enabled
		assert(unbound_tid >= 0);
		select_texture(unbound_tid, 0);
		set_color_d(unbound_color);
		unbound_geom.render(is_shadow_pass);
	}
	
	// render all materials (opaque then transparen)
	for (unsigned pass = 0; pass < 2; ++pass) { // opaque, transparent
		for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
			if ((unsigned)m->is_partial_transparent() == pass) m->render(tmgr, unbound_tid, is_shadow_pass);
		}
	}
	if (do_cull) glDisable(GL_CULL_FACE);
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

	set_lighted_sides(2);
	set_fill_mode();
	glDisable(GL_LIGHTING); // custom lighting calculations from this point on
	BLACK.do_glColor();
	set_color_a(BLACK); // ambient will be set by indirect lighting in the shader
	set_specular(0.0, 1.0);
	shader_t s;
	colorRGBA orig_fog_color;
	float const min_alpha(0.5); // since we're using alpha masks we must set min_alpha > 0.0
	if (!is_shadow_pass) orig_fog_color = setup_smoke_shaders(s, min_alpha, 0, 0, 1, 1, 1, 1, 0, shadow_map_enabled());

	for (iterator m = begin(); m != end(); ++m) {
		m->render(is_shadow_pass);
	}
	if (!is_shadow_pass) end_smoke_shaders(s, orig_fog_color);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	set_lighted_sides(1);
	set_specular(0.0, 1.0);
}


void free_model_context() {
	all_models.free_context();
}

void render_models(bool shadow_pass) {
	all_models.render(shadow_pass);
}


