// 3D World - 3D Model Rendering Code
// by Frank Gennari
// 8/17/11

#include "model3d.h"
#include "shaders.h"
#include "gl_ext_arb.h"

extern bool group_back_face_cull, enable_model3d_tex_comp;


model3ds all_models;


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


void vntc_vect_t::render_array(bool is_shadow_pass) {

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
	glDrawArrays(GL_TRIANGLES, 0, size());
	bind_vbo(0);
}


void vntc_vect_t::free_vbos() {

	delete_vbo(vbo);
	vbo = 0;
}


void geom_data_t::add_poly(vntc_vect_t const &poly) {
	
	polygons.push_back(poly);

	if (poly.size() == 3) { // triangle
		for (unsigned i = 0; i < 3; ++i) {triangles.push_back(poly[i]);}
		return;
	}
	if (poly.size() == 4) {
		static vector<point> poly_pts;
		poly_pts.resize(4);
		for (unsigned d = 0; d < 4; ++d) {poly_pts[d] = poly[d].v;}

		if (is_poly_convex(poly_pts)) { // quad
			unsigned const ixs[6] = {0,1,2,0,2,3};
			for (unsigned i = 0; i < 6; ++i) {triangles.push_back(poly[ixs[i]]);}
			return;
		}
	}
	// split to triangles - WRITE
}


void geom_data_t::render_polygons(bool is_shadow_pass) const {

	for (vector<vntc_vect_t>::const_iterator p = polygons.begin(); p != polygons.end(); ++p) {
		glBegin(GL_POLYGON);
		p->render(is_shadow_pass);
		glEnd();
	}
}


void material_t::render(texture_manager const &tm, int default_tid, bool is_shadow_pass) {

	if (geom.empty() || alpha == 0.0) return; // empty or transparent

	if (!is_shadow_pass) {
		int const tex_id(get_render_texture());
		
		if (tex_id >= 0) {
			tm.bind_texture(tex_id);
		}
		else {
			select_texture(((default_tid >= 0) ? default_tid : WHITE_TEX), 0); // no texture specified - use white texture
		}
		if (alpha < 1.0 && ni != 1.0) {
			// set index of refraction (and reset it at the end)
		}
		float const spec_val((ks.R + ks.G + ks.B)/3.0);
		set_specular(spec_val, ns);
		//set_color_a(colorRGBA(ka, alpha));
		//set_color_d(colorRGBA(kd, alpha));
		set_color_d(colorRGBA(kd, alpha) + colorRGBA(ka, 0.0));
		set_color_e(colorRGBA(ke, alpha));
	}
	geom.render_array(is_shadow_pass);

	if (!is_shadow_pass) {
		set_color_e(BLACK);
		set_specular(0.0, 1.0);
	}
}


// creation and query
void model3d::add_polygon(vntc_vect_t const &poly, int mat_id, vector<vector<point> > *ppts) {

	//assert(mat_id >= 0); // must be set/valid - FIXME: too strict?

	if (mat_id < 0) {
		unbound_geom.add_poly(poly);
	}
	else {
		assert((unsigned)mat_id < materials.size());
		materials[mat_id].geom.add_poly(poly);
	}
	if (ppts) { // FIXME: split polyg if needed
		ppts->push_back(vector<point>());

		for (vntc_vect_t::const_iterator i = poly.begin(); i != poly.end(); ++i) {
			ppts->back().push_back(i->v);
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


unsigned texture_manager::create_texture(string const &fn, bool verbose) {

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
	textures.back().do_compress = enable_model3d_tex_comp;
	return tid; // can't fail
}


// clear/free
void model3d::clear() {

	free_context();
	unbound_geom.clear();
	materials.clear();
	undef_materials.clear();
	mat_map.clear();
}


void model3d::free_context() {

	for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
		m->geom.free_context();
	}
	unbound_geom.free_context();
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


// texture loading
void texture_manager::ensure_texture_loaded(texture_t &t) const {
	if (!t.data) t.load(-1);
	assert(t.data);
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

	assert((unsigned)tid < textures.size());
	assert(textures[tid].tid > 0);
	glBindTexture(GL_TEXTURE_2D, textures[tid].tid);
}


void model3d::load_all_used_tids() {

	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (m->geom.empty()) continue;
		tm.ensure_tid_loaded(m->get_render_texture()); // only one tid for now
		tm.ensure_tid_loaded(m->alpha_tid);
		// FIXME: use alpha_tid
	}
}


void model3d::bind_all_used_tids() {

	load_all_used_tids();
		
	for (deque<material_t>::const_iterator m = materials.begin(); m != materials.end(); ++m) {
		if (m->geom.empty()) continue;
		tm.ensure_tid_bound(m->get_render_texture()); // only one tid for now
	}
}


// rendering
void model3d::render(bool is_shadow_pass) { // const?

	if (!is_shadow_pass) bind_all_used_tids();
	bool const do_cull(group_back_face_cull && !is_shadow_pass);
	if (do_cull) glEnable(GL_CULL_FACE);

	// render geom that was not bound to a material
	if (unbound_color.alpha > 0.0) { // enabled
		assert(unbound_tid >= 0);
		select_texture(unbound_tid, 0);
		set_color_d(unbound_color);
		unbound_geom.render_array(is_shadow_pass);
	}
	
	// render all materials (opaque then transparen)
	for (unsigned pass = 0; pass < 2; ++pass) { // opaque, transparent
		for (deque<material_t>::iterator m = materials.begin(); m != materials.end(); ++m) {
			if ((unsigned)m->is_partial_transparent() == pass) m->render(tm, unbound_tid, is_shadow_pass);
		}
	}
	if (do_cull) glDisable(GL_CULL_FACE);
}


void model3ds::clear() {

	for (iterator m = begin(); m != end(); ++m) {
		m->clear();
	}
	deque<model3d>::clear();
	tm.clear();
}


void model3ds::free_context() {

	for (iterator m = begin(); m != end(); ++m) {
		m->free_context();
	}
	tm.free_tids();
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
	if (!is_shadow_pass) orig_fog_color = setup_smoke_shaders(s, 0.0, 0, 0, 1, 1, 1, 1, 0, shadow_map_enabled());

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


