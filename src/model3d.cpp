// 3D World - 3D Model Rendering Code
// by Frank Gennari
// 8/17/11

#include "model3d.h"
#include "shaders.h"

extern bool group_back_face_cull;


model3ds all_models;


void polygon_t::render(bool textured) const {

	assert(size() >= 3);

	for (const_iterator v = begin(); v != end(); ++v) {
		if (textured) glTexCoord2fv(v->t);
		v->n.do_glNormal();
		v->v.do_glVertex();
	}
}


void poly_group::render(bool textured) const {

	// FIXME: very inefficient - use glDrawArrays() or vbo with polygons split into triangles
	for (vector<polygon_t>::const_iterator p = polygons.begin(); p != polygons.end(); ++p) {
		glBegin(GL_POLYGON);
		p->render(textured);
		glEnd();
	}
}


void material_t::render(deque<texture_t> const &textures, bool is_shadow_pass) const {

	if (geom.empty()) return; // nothing to do
	bool textured(0);

	if (!is_shadow_pass) {
		int const tex_id(get_render_texture());
		textured = (tex_id >= 0);

		if (textured) {
			assert((unsigned)tex_id < textures.size());
			assert(textures[tex_id].tid > 0);
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, textures[tex_id].tid);
		}
		if (alpha < 1.0 && ni != 1.0) {
			// set index of refraction (and reset it at the end)
		}
		float const spec_val((ks.R + ks.G + ks.B)/3.0);
		set_specular(spec_val, ns);
		set_color_a(colorRGBA(ka, alpha));
		set_color_d(colorRGBA(kd, alpha));
		set_color_e(colorRGBA(kd, alpha));
	}
	// FIXME: sort back to front if is_partial_transparent()?
	geom.render(textured);

	if (!is_shadow_pass) {
		set_color_e(BLACK);
		set_specular(0.0, 1.0);
		if (textured) glDisable(GL_TEXTURE_2D);
	}
}


// creation and query
void model3d::add_polygon(polygon_t const &poly, int mat_id) {

	//assert(mat_id >= 0); // must be set/valid - FIXME: too strict?

	if (mat_id < 0) {
		unbound_geom.add_poly(poly);
	}
	else {
		assert((unsigned)mat_id < materials.size());
		materials[mat_id].geom.add_poly(poly);
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


unsigned model3d::create_texture(string const &fn, bool verbose) {

	string_map_t::const_iterator it(tex_map.find(fn));

	if (it != tex_map.end()) { // found (already loaded)
		assert(it->second < textures.size());
		return it->second;
	}
	unsigned const tid(textures.size());
	tex_map[fn] = tid;
	if (verbose) cout << "loading texture " << fn << endl;
	// type format width height wrap ncolors use_mipmaps name [bump_name [id [color]]]
	textures.push_back(texture_t(0, 4, 0, 0, 1, 3, 1, fn)); // always RGB targa wrapped+mipmap
	return tid; // can't fail
}


// clear/free
void model3d::clear() {

	unbound_geom.clear();
	materials.clear();
	undef_materials.clear();
	mat_map.clear();
	free_textures();
	textures.clear();
	tex_map.clear();
}


void model3d::free_tids() {

	for (deque<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {
		t->gl_delete();
	}
}

void model3d::free_textures() {

	for (deque<texture_t>::iterator t = textures.begin(); t != textures.end(); ++t) {
		t->free();
	}
}


// texture loading
void model3d::ensure_texture_loaded(texture_t &t) const {
	if (!t.data) t.load(-1);
	assert(t.data);
}

void model3d::ensure_tid_loaded(int tid) {

	if (tid < 0) return; // not allocated
	assert((unsigned)tid < textures.size());
	ensure_texture_loaded(textures[tid]);
}


void model3d::ensure_tid_bound(int tid) {

	if (tid < 0) return; // not allocated
	assert((unsigned)tid < textures.size());
	textures[tid].check_init();
}


void model3d::load_all_used_tids() {

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (m->geom.empty()) continue;
		ensure_tid_loaded(m->get_render_texture()); // only one tid for now
		ensure_tid_loaded(m->alpha_tid);
		// FIXME: use alpha_tid
	}
}


void model3d::bind_all_used_tids() {

	load_all_used_tids();
		
	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (m->geom.empty()) continue;
		ensure_tid_bound(m->get_render_texture()); // only one tid for now
	}
}


// rendering
void model3d::render(bool is_shadow_pass) { // const?

	bind_all_used_tids();
	if (group_back_face_cull) glEnable(GL_CULL_FACE);
	unbound_geom.render(0);
	
	for (unsigned pass = 0; pass < 2; ++pass) { // opaque, transparent
		for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
			if ((unsigned)m->is_partial_transparent() == pass) m->render(textures, is_shadow_pass);
		}
	}
	if (group_back_face_cull) glDisable(GL_CULL_FACE);
}


void model3ds::clear() {

	for (iterator m = begin(); m != end(); ++m) {
		m->clear();
	}
}


void model3ds::free_tids() {

	for (iterator m = begin(); m != end(); ++m) {
		m->free_tids();
	}
}


void model3ds::render(bool is_shadow_pass) {
	shader_t s;
	colorRGBA const orig_fog_color(setup_smoke_shaders(s, 0.0, 0, 0, 1, 1, 1, 1, 0, shadow_map_enabled()));

	for (iterator m = begin(); m != end(); ++m) {
		m->render(is_shadow_pass);
	}
	end_smoke_shaders(s, orig_fog_color);
}


void free_model_textures() {
	all_models.free_tids();
}

void render_models(bool shadow_pass) {
	all_models.render(shadow_pass);
}


