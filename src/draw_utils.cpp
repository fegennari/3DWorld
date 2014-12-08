// 3D World - Drawing Utility Classes
// by Frank Gennari
// 6/15/12
#include "draw_utils.h"
#include "function_registry.h"
#include "gl_ext_arb.h"
#include "shaders.h"


extern int window_height, display_mode;
extern bool use_core_context;
extern shader_t *cur_shader;


void colorRGBA::set_for_cur_shader() const {

	assert(cur_shader != NULL);
	cur_shader->set_cur_color(*this);
}


void set_array_client_state(bool va, bool tca, bool na, bool ca, bool actually_set_state) {

	assert(cur_shader != NULL);
	if (actually_set_state) {cur_shader->enable_vnct_atribs(va, tca, na, ca);}
	else {check_mvm_update();}
}


void const *ptr_add(void const *p, unsigned off) {return (unsigned char const *)(unsigned(p) + off);}

void set_vn_ptrs(unsigned stride, bool comp, void const *vbo_ptr_offset) {
	assert(cur_shader);
	cur_shader->set_vertex_ptr(stride, vbo_ptr_offset);
	cur_shader->set_normal_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(point)), comp);
}

void vert_wrap_t::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 0, 0, set_state);
	cur_shader->set_vertex_ptr(sizeof(point), vbo_ptr_offset);
}

void vert_tc_t::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 0, 0, set_state);
	unsigned const stride(sizeof(vert_tc_t));
	cur_shader->set_vertex_ptr(stride, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_wrap_t)), 0);
}

void vert_norm::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 1, 0, set_state);
	set_vn_ptrs(sizeof(vert_norm), 0, vbo_ptr_offset);
}

void vert_norm_comp::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 1, 0, set_state);
	set_vn_ptrs(sizeof(vert_norm_comp), 1, vbo_ptr_offset);
}

void vert_norm_comp_tc::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 1, 0, set_state);
	unsigned const stride(sizeof(vert_norm_comp_tc));
	set_vn_ptrs(stride, 1, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_comp)), 0);
}

void vert_norm_comp_tc_comp::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 1, 0, set_state);
	unsigned const stride(sizeof(vert_norm_comp_tc_comp));
	set_vn_ptrs(stride, 1, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_comp)), 1);
}

void vert_norm_tc::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 1, 0, set_state);
	unsigned const stride(sizeof(vert_norm_tc));
	set_vn_ptrs(stride, 0, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm)), 0);
}

void vert_norm_tc_tan::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 1, 0, set_state);
	unsigned const stride(sizeof(vert_norm_tc_tan));
	set_vn_ptrs(stride, 0, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm)), 0);
	// Note: would be cleaner but less efficient to use cur_shader->get_attrib_loc("tangent")
	int const loc(cur_shader->attrib_loc_by_ix(TANGENT_ATTR, 1)); // okay if fails
	if (loc >= 0) {
		glEnableVertexAttribArray(loc);
		glVertexAttribPointer(loc, 4, GL_FLOAT, GL_FALSE, stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_tc)));
	}
}

void vert_norm_tc_tan::unset_attrs() {
	assert(cur_shader);
	int const loc(cur_shader->attrib_loc_by_ix(TANGENT_ATTR, 1)); // okay if fails
	if (loc >= 0) {glDisableVertexAttribArray(loc);}
}

void vert_color::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 0, 1, set_state);
	unsigned const stride(sizeof(vert_color));
	cur_shader->set_vertex_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(color_wrapper)));
	cur_shader->set_color4_ptr(stride, vbo_ptr_offset, 1);
}

void vert_norm_color::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 1, 1, set_state);
	unsigned const stride(sizeof(vert_norm_color));
	set_vn_ptrs(stride, 0, vbo_ptr_offset);
	cur_shader->set_color4_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm)), 1);
}

void vert_norm_comp_color::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 1, 1, set_state);
	unsigned const stride(sizeof(vert_norm_comp_color));
	set_vn_ptrs(stride, 1, vbo_ptr_offset);
	cur_shader->set_color4_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_comp)), 1);
}

void vert_norm_tc_color::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 1, 1, set_state);
	unsigned const stride(sizeof(vert_norm_tc_color));
	set_vn_ptrs(stride, 0, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm)), 0);
	cur_shader->set_color4_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_tc)), 1);
}

void vert_tc_color::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 0, 1, set_state);
	unsigned const stride(sizeof(vert_tc_color));
	cur_shader->set_vertex_ptr(stride, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(point)), 0);
	cur_shader->set_color4_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_tc_t)), 1);
}

void vert_norm_comp_tc_color::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 1, 1, set_state);
	unsigned const stride(sizeof(vert_norm_comp_tc_color));
	set_vn_ptrs(stride, 1, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_comp)), 0);
	cur_shader->set_color4_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_comp_tc)), 1);
}

void vert_norm_comp_tc_comp_color::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 1, 1, 1, set_state);
	unsigned const stride(sizeof(vert_norm_comp_tc_comp_color));
	set_vn_ptrs(stride, 1, vbo_ptr_offset);
	cur_shader->set_tcoord_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_comp)), 1);
	cur_shader->set_color4_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm_comp_tc_comp)), 1);
}

void vert_norm_texp::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 1, 0, set_state);
	unsigned const stride(sizeof(vert_norm_texp));
	set_vn_ptrs(stride, 0, vbo_ptr_offset);
	
	for (unsigned d = 0; d < 2; ++d) {
		unsigned const loc(cur_shader->attrib_loc_by_ix(d ? TEX0_T_ATTR : TEX0_S_ATTR));
		if (loc >= 0) {glEnableVertexAttribArray(loc);}
		if (loc >= 0) {glVertexAttribPointer(loc, 4, GL_FLOAT, GL_FALSE, stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm)+4*d*sizeof(float)));}
	}
}
void vert_norm_texp::unset_attrs() {
	for (unsigned d = 0; d < 2; ++d) {
		unsigned const loc(cur_shader->attrib_loc_by_ix(d ? TEX0_T_ATTR : TEX0_S_ATTR));
		if (loc >= 0) {glDisableVertexAttribArray(loc);}
	}
}


// required when using a core context, but slow
vbo_ring_buffer_t vbo_ring_buffer(1 << 23); // 8MB

void clear_vbo_ring_buffer() {vbo_ring_buffer.free_vbo();}

void const *get_dynamic_vbo_ptr(void const *const verts, unsigned size_bytes) {
	assert(verts != NULL && size_bytes > 0);
	return vbo_ring_buffer.add_verts_bind_vbo(verts, size_bytes);
}

bool bind_temp_vbo_from_verts(void const *const verts, unsigned count, unsigned vert_size, void const *&vbo_ptr_offset) {
	if (!use_core_context) return 0; // not needed
	vbo_ptr_offset = get_dynamic_vbo_ptr(verts, count*vert_size);
	return 1;
}

void unbind_temp_vbo() {
	if (!use_core_context) return; // not needed
	bind_vbo(0);
}


void pt_line_drawer::add_textured_pt(point const &v, colorRGBA c, int tid) {

	if (tid >= 0) c = c.modulate_with(texture_color(tid));
	vector3d const view_dir(get_camera_pos(), v);
	add_pt(v, view_dir, c);
}

void pt_line_drawer::add_textured_line(point const &v1, point const &v2, colorRGBA c, int tid) {

	if (tid >= 0) c = c.modulate_with(texture_color(tid));
	vector3d view_dir(get_camera_pos(), (v1 + v2)*0.5);
	orthogonalize_dir(view_dir, (v2 - v1), view_dir, 0);
	add_line(v1, view_dir, c, v2, view_dir, c);
}

void pt_line_drawer::draw() const {
	
	assert(!(lines.size() & 1));
	if (!points.empty()) {draw_verts(points, GL_POINTS);}
	if (!lines.empty ()) {draw_verts(lines,  GL_LINES );}
	//cout << "mem: " << get_mem() << endl;
}


void pt_line_drawer_no_lighting_t::draw() const {

	assert(!(lines.size() & 1));
	draw_verts(points, GL_POINTS);
	draw_verts(lines,  GL_LINES);
}

void pt_line_drawer_no_lighting_t::draw_vbo() {

	if (!points.empty()) {
		create_bind_vbo_and_upload(vbo[0], points);
		draw_verts<vert_color>(NULL, points.size(), GL_POINTS);
	}
	if (!lines.empty()) {
		assert(!(lines.size() & 1));
		create_bind_vbo_and_upload(vbo[1], lines);
		draw_verts<vert_color>(NULL, lines.size(), GL_LINES);
	}
	bind_vbo(0);
}

void pt_line_drawer_no_lighting_t::free_vbo() {
	for (unsigned d = 0; d < 2; ++d) {delete_and_zero_vbo(vbo[d]);}
}


template<class vert_type_t> void sized_vert_t<vert_type_t>::set_size_attr(unsigned stride, void const *vbo_ptr_offset) {
	int const loc(cur_shader->get_attrib_loc("point_size"));
	assert(loc > 0);
	glEnableVertexAttribArray(loc);
	glVertexAttribPointer(loc, 4, GL_FLOAT, GL_FALSE, stride, ptr_add(vbo_ptr_offset, sizeof(vert_type_t)));
}
template<class vert_type_t> void sized_vert_t<vert_type_t>::unset_attrs() {
	int const loc(cur_shader->get_attrib_loc("point_size"));
	assert(loc > 0);
	glDisableVertexAttribArray(loc);
}

template<> void sized_vert_t<vert_norm_color>::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 1, 1, set_state);
	unsigned const stride(sizeof(sized_vert_t<vert_norm_color>));
	set_vn_ptrs(stride, 0, vbo_ptr_offset);
	cur_shader->set_color4_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(vert_norm)), 1);
	set_size_attr(stride, vbo_ptr_offset);
}
template<> void sized_vert_t<vert_color>::set_vbo_arrays(bool set_state, void const *vbo_ptr_offset) {
	set_array_client_state(1, 0, 0, 1, set_state);
	unsigned const stride(sizeof(sized_vert_t<vert_color>));
	cur_shader->set_vertex_ptr(stride, ptr_add(vbo_ptr_offset, sizeof(color_wrapper))); // color is first in the struct
	cur_shader->set_color4_ptr(stride, vbo_ptr_offset, 1);
	set_size_attr(stride, vbo_ptr_offset);
}


template<class vert_type_t> void point_sprite_drawer_t<vert_type_t>::draw(int tid, float const_point_size, bool enable_lighting) const {

	if (empty()) return;
	shader_t s;
#if 1 // point sprite variant
	if (const_point_size) {s.set_prefix("#define CONSTANT_PT_SIZE", 0);} // VS

	if (enable_lighting) {
		s.setup_enabled_lights(2, 1); // sun and moon VS lighting
		s.set_prefix("#define ENABLE_LIGHTING", 0); // VS
		s.set_vert_shader("ads_lighting.part*+point_sprite"); // no fog
	}
	else {
		s.set_vert_shader("point_sprite");
	}
	s.set_frag_shader("point_sprite_texture");
	s.begin_shader();
	s.add_uniform_float("point_scale", 2.0*(const_point_size ? const_point_size : window_height)); // diameter = 2*radius
#else // geometry shader variant - doesn't support const_point_size or lighting
	s.set_prefix("#define SIZE_FROM_ATTRIB", 2); // GS
	s.set_vert_shader("particle_draw");
	s.set_frag_shader("simple_texture");
	s.set_geom_shader("pt_billboard_tri"); // point => 1 triangle
	s.begin_shader();
#endif
	s.add_uniform_int("tex0", 0);
	s.add_uniform_float("min_alpha", 0.0);
	set_point_sprite_mode(1);
	select_texture(tid);

	if (!const_point_size) { // use variable attribute point size
		assert(points.size() == sizes.size());
		vector<sized_vert_t<vert_type_t> > verts(points.size());
		for (unsigned i = 0; i < points.size(); ++i) {verts[i] = sized_vert_t<vert_type_t>(points[i], sizes[i]);}
		draw_verts(verts, GL_POINTS);
	}
	else {
		draw_verts(points, GL_POINTS);
	}
	set_point_sprite_mode(0);
	s.end_shader();
}

template class point_sprite_drawer_t<vert_color     >;
template class point_sprite_drawer_t<vert_norm_color>;


void quad_batch_draw::add_quad_pts(point const pts[4], colorRGBA const &c, vector3d const &n, tex_range_t const &tr) {

	float const t[4][2] = {{tr.x1,tr.y1}, {tr.x2,tr.y1}, {tr.x2,tr.y2}, {tr.x1,tr.y2}};
	unsigned const v[6] = {0,2,1, 0,3,2}; // Note: reversed from quad_to_tris_ixs
	color_wrapper cw;
	cw.set_c4(c);

	for (unsigned i = 0; i < 6; ++i) {
		verts.push_back(vert_norm_tc_color(pts[v[i]], n, t[v[i]][0], t[v[i]][1], cw.c, 1));
	}
}

void quad_batch_draw::add_quad_dirs(point const &pos, vector3d const &dx, vector3d const &dy,
	colorRGBA const &c, vector3d const &n, tex_range_t const &tr)
{
	point const pts[4] = {(pos - dx - dy), (pos + dx - dy), (pos + dx + dy), (pos - dx + dy)};
	add_quad_pts(pts, c, n, tr);
}

// unused
void quad_batch_draw::add_quad_dirs_single_tri(point const &pos, vector3d const &dx, vector3d const &dy, colorRGBA const &c, vector3d const &n) {

	// add a single triangle with an inscribed square, using tex coords outside the [0,1] range
	color_wrapper cw;
	cw.set_c4(c);
	verts.push_back(vert_norm_tc_color(pos-2*dx-dy, n, -0.5, 0.0, cw.c, 1));
	verts.push_back(vert_norm_tc_color(pos+2*dx-dy, n,  1.5, 0.0, cw.c, 1));
	verts.push_back(vert_norm_tc_color(pos+3*dy,    n,  0.0, 2.0, cw.c, 1));
}

void quad_batch_draw::add_xlated_billboard(point const &pos, point const &xlate, point const &viewer, vector3d const &up_dir,
	colorRGBA const &c, float xsize, float ysize, tex_range_t const &tr, bool minimize_fill)
{
	vector3d const vdir(viewer - pos); // z
	vector3d const v1((cross_product(vdir, up_dir).get_norm())*xsize); // x (what if colinear?)
	vector3d const v2(cross_product(v1, vdir).get_norm()*ysize); // y
	vector3d const normal(vdir.get_norm());

	if (minimize_fill) { // draw as octagon
		assert(tr.x1 == 0 && tr.y1 == 0 && tr.x2 == 1 && tr.y2 == 1);
		float const p[8][2] = {{0.7,0.0}, {1.0,0.3}, {1.0,0.7}, {0.7,1.0}, {0.3,1.0}, {0.0,0.7}, {0.0,0.3}, {0.3,0.0}};
		unsigned const v[18] = {0,1,7 ,1,6,7, 1,2,6, 2,5,6, 2,3,5, 3,4,5};
		color_wrapper cw;
		cw.set_c4(c);

		for (unsigned i = 0; i < 18; ++i) {
			float const tcx(p[v[i]][0]), tcy(p[v[i]][1]);
			point const p(xlate + v1*(2.0*tcx - 1.0) + v2*(2.0*tcy - 1.0));
			verts.push_back(vert_norm_tc_color(p, normal, tcx, tcy, cw.c, 1));
		}
	}
	else { // draw as quad (2 triangles)
		add_quad_dirs(xlate, v1, v2, c, normal, tr);
	}
}

void quad_batch_draw::add_animated_billboard(point const &pos, point const &viewer, vector3d const &up_dir, colorRGBA const &c, float xsize, float ysize, float timescale) {

	// fixed 4x4 animation
	int const frame_id(max(0, min(15, int(16*timescale)))), tx(frame_id&3), ty(frame_id>>2);
	point const gpos(make_pt_global(pos));
	add_billboard(gpos, (viewer + gpos - pos), up_dir, c, xsize, ysize, tex_range_t::from_atlas(tx, ty, 4, 4)); // upside down
}


void quad_batch_draw::draw_as_flares_and_clear(int flare_tex) { // Note: used in flows where texturing is always enabled

	glDepthMask(GL_FALSE);
	select_texture(flare_tex);
	draw_and_clear();
	glDepthMask(GL_TRUE);
}


template<typename T> void indexed_mesh_draw<T>::clear() {

	verts.clear();
	free_context();
	nx = ny = ivbo_size = 0;
}

template<typename T> void indexed_mesh_draw<T>::init(unsigned nx_, unsigned ny_) {

	if (nx == nx_ && ny == ny_ && ivbo) return; // already setup
	nx = nx_; ny = ny_;
	assert(nx > 0 && ny > 0);
	verts.resize((nx+1)*(ny+1));
	vector<unsigned> indices(6*nx*ny);
	ivbo_size = indices.size();
		
	for (unsigned y = 0; y < ny; ++y) {
		for (unsigned x = 0; x < nx; ++x) {
			unsigned const iix(6*(y*nx + x)), vix(get_vert_ix(x, y)); // 6 verts / 2 tris / 1 quad
			indices[iix+0] = indices[iix+3] = vix; // 0
			indices[iix+1] = vix + 1; // 1
			indices[iix+2] = indices[iix+4] = vix + (nx+1) + 1; // 2
			indices[iix+5] = vix + (nx+1); // 3
		}
	}
	free_context();
	create_vbo_and_upload(ivbo, indices, 1, 1);
}

template<typename T> void indexed_mesh_draw<T>::render() const {

	if (verts.empty()) return;
	assert(ivbo && ivbo_size > 0);
	bind_vbo(ivbo, 1);
	set_ptr_state(&verts.front(), verts.size());
	glDrawRangeElements(GL_TRIANGLES, 0, verts.size(), ivbo_size, GL_UNSIGNED_INT, NULL);
	bind_vbo(0, 1);
	unset_ptr_state(&verts.front());
}

template<typename T> void indexed_mesh_draw<T>::render_z_plane(float x1, float y1, float x2, float y2, float zval, unsigned nx_, unsigned ny_) {

	assert(x1 < x2 && y1 < y2);
	init(nx_, ny_);
	float const xinc((x2 - x1)/nx), yinc((y2 - y1)/ny);
	float yval(y1);

	for (unsigned y = 0; y <= ny; ++y) {
		float xval(x1);

		for (unsigned x = 0; x <= nx; ++x) {
			set_vert(x, y, point(xval, yval, zval));
			xval += xinc;
		}
		yval += yinc;
	}
	render();
}

template class indexed_mesh_draw<vert_wrap_t>;


// Note: these classes are actually declared in gl_ext_arb.h
// simplified and optimized version of draw_cube_mapped_sphere() (no center, tex coords, or normals, and radius=1.0)
cube_map_sphere_drawer_t::cube_map_sphere_drawer_t(unsigned ndiv) {

	assert(!vbo);
	assert(ndiv > 0);
	nverts   = 6*(ndiv+1)*(ndiv+1);
	nindices = 6*ndiv*(2*(ndiv+1)+1);
	vector<vert_wrap_t> verts; verts.reserve(nverts);
	vector<unsigned> indices; indices.reserve(nindices);
	float const vstep(2.0/ndiv);

	for (unsigned i = 0; i < 3; ++i) { // iterate over dimensions
		for (unsigned j = 0; j < 2; ++j) { // iterate over opposing sides, min then max
			unsigned const d1(i), d2((i+1)%3), dn((i+2)%3), start_vix(verts.size());
			point pt;
			pt[dn] = (j ? 1.0 : -1.0);
				
			for (unsigned s = 0; s <= ndiv; ++s) { // create vertex data
				pt[d1] = -1.0 + s*vstep;

				for (unsigned T = 0; T <= ndiv; ++T) {
					pt[d2] = -1.0 + T*vstep;
					verts.push_back(pt.get_norm()); // Note: could also normalize in the shader
				}
			} // for s
			for (unsigned s = 0; s < ndiv; ++s) { // create indices
				for (unsigned T = 0; T <= ndiv; ++T) {
					unsigned const ix_off(start_vix + s*(ndiv+1) + (j ? T : ndiv-T)); // reverse between sides
					indices.push_back(ix_off); // current s row
					indices.push_back(ix_off + ndiv+1); // next s row
				}
				indices.push_back(PRIMITIVE_RESTART_IX); // restart the strip
			} // for s
		} // for j
	} // for i
	assert(verts.size()   == nverts);
	assert(indices.size() == nindices);
	create_and_upload(verts, indices);
}

void cube_map_sphere_drawer_t::draw() const {

	pre_render();
	vert_wrap_t::set_vbo_arrays();
	glEnable(GL_PRIMITIVE_RESTART);
	glPrimitiveRestartIndex(PRIMITIVE_RESTART_IX);
	glDrawRangeElements(GL_TRIANGLE_STRIP, 0, nverts, nindices, GL_UNSIGNED_INT, nullptr);
	glDisable(GL_PRIMITIVE_RESTART);
	post_render();
}

void cube_map_sphere_manager_t::draw_sphere(unsigned ndiv) {
	ndiv_sphere_map_t::iterator it(cached.find(ndiv));
	if (it == cached.end()) {it = cached.insert(make_pair(ndiv, cube_map_sphere_drawer_t(ndiv))).first;} // create a new one
	it->second.draw();
}
void cube_map_sphere_manager_t::clear() {
	for (ndiv_sphere_map_t::iterator i = cached.begin(); i != cached.end(); ++i) {i->second.clear_vbos();}
	cached.clear();
}


template< typename vert_type_t >
unsigned vbo_block_manager_t<vert_type_t>::get_offset_for_last_points_added() {

	if (offsets.empty()) {offsets.push_back(0);} // start at 0
	unsigned const next_ix(offsets.size() - 1);
	offsets.push_back(pts.size()); // range will be [start_ix, start_ix+p.size()]
	return next_ix;
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::add_points_int(vector<vert_type_t> &dest, typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color) {

	assert(p != NULL && npts > 0);
	color_wrapper cw;
	cw.set_c4(color);
	for (unsigned i = 0; i < npts; ++i) {dest.push_back(vert_type_t(p[i], cw));}
}

template<>
void vbo_block_manager_t<vert_norm_tc>::add_points_int(vector<vert_norm_tc> &dest, vert_norm_tc const *const p, unsigned npts, colorRGBA const &color) { // color is ignored
	assert(p != NULL && npts > 0);
	copy(p, p+npts, back_inserter(dest));
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::render_range(unsigned six, unsigned eix, unsigned num_instances) const {

	assert(six < eix && eix < offsets.size());
	unsigned const count(offsets[eix] - offsets[six]);
	// Note: currently always used to render quads, but can be made more general in the future
	check_mvm_update();

	if (use_core_context) {
		draw_quads_as_tris(count, offsets[six], num_instances);
	}
	else {
		glDrawArraysInstanced(GL_QUADS, offsets[six], count, num_instances);
	}
}

template< typename vert_type_t >
bool vbo_block_manager_t<vert_type_t>::upload() {

	if (vbo || !has_data()) return 0; // already uploaded or empty
	assert(!pts.empty());
	if (offsets.empty()) {offsets.push_back(0);} // start at 0
	if (offsets.back() != pts.size()) {offsets.push_back(pts.size());} // add terminator
	create_vbo_and_upload(vbo, pts, 0, 1);
	return 1;
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::update_range(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color, unsigned six, unsigned eix) {

	if (!vbo) return; // vbo not uploaded (assertion?)
	assert(six < eix && eix < offsets.size());
	unsigned const start(offsets[six]), update_size(offsets[eix] - start);
	assert(npts == update_size);
	vector<vert_type_t> update_verts;
	add_points_int(update_verts, p, npts, color);
	pre_render();
	upload_vbo_sub_data(&update_verts.front(), start*sizeof(vert_type_t), update_size*sizeof(vert_type_t));
	post_render();
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::begin_render() const {
	if (!has_data()) return;
	pre_render();
	vert_type_t::set_vbo_arrays();
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::clear() {
	clear_points();
	temp_points.clear();
	offsets.clear();
	clear_vbo();
}

// explicit template instantiations
template class vbo_block_manager_t<vert_color>;
template class vbo_block_manager_t<vert_norm_comp_color>;
template class vbo_block_manager_t<vert_norm_tc_color>;
template class vbo_block_manager_t<vert_norm_tc>;


class quad_ix_buffer_t {

	unsigned ivbo_16, ivbo_32;
	unsigned size_16, size_32;

	template< typename T > static void ensure_quad_ixs(unsigned &ivbo, unsigned size) {
		if (ivbo != 0) return; // vbo already valid
		assert((size % 6) == 0); // must be in groups of 2 tris = 1 quad
		unsigned const num_quads(size/6);
		vector<T> ixs(size);

		// Note: quad is split along a different axis from GL quads, so interpolation is different
		for (unsigned q = 0; q < num_quads; ++q) {
			for (unsigned i = 0; i < 6; ++i) {ixs[6*q+i] = 4*q + quad_to_tris_ixs[i];}
		}
		create_vbo_and_upload(ivbo, ixs, 1);
	}

public:
	quad_ix_buffer_t() : ivbo_16(0), ivbo_32(0), size_16(0), size_32(0) {}
	
	void free_context() {
		delete_and_zero_vbo(ivbo_16);
		delete_and_zero_vbo(ivbo_32);
		size_16 = size_32 = 0;
	}

	bool bind_quads_as_tris_ivbo(unsigned num_quad_verts) {
		assert((num_quad_verts & 3) == 0); // must be a multiple of 4
		unsigned const num_tri_verts(6*(num_quad_verts/4));
		unsigned const max_quad_verts = 65532; // largest multiple of 4 and 6 smaller than 2^16
		bool const use_32_bit(num_quad_verts > max_quad_verts);
		unsigned &ivbo    (use_32_bit ? ivbo_32 : ivbo_16);
		unsigned &cur_size(use_32_bit ? size_32 : size_16);
		
		if (num_tri_verts > cur_size) { // increase the size
			delete_vbo(ivbo); ivbo = 0;
			cur_size = max(96U, max(num_tri_verts, 2U*cur_size)); // at least double
		}
		assert(num_tri_verts <= cur_size);

		if (use_32_bit) { // use 32-bit verts
			ensure_quad_ixs<unsigned>(ivbo, cur_size);
		}
		else { // use 16-bit verts
			ensure_quad_ixs<unsigned short>(ivbo, cur_size);
		}
		assert(ivbo != 0);
		bind_vbo(ivbo, 1);
		return use_32_bit;
	}

	void draw_quads_as_tris(unsigned num_quad_verts, unsigned start_quad_vert, unsigned num_instances) { // # vertices
		if (num_quad_verts == 0) return; // nothing to do
		assert((num_quad_verts & 3) == 0 && (start_quad_vert & 3) == 0); // must be a multiple of 4
		unsigned const end_quad_vert(start_quad_vert + num_quad_verts);
		unsigned const num_tri_verts(6*(num_quad_verts/4)), start_tri_vert(6*(start_quad_vert/4)); // # indices
		bool const use_32_bit(bind_quads_as_tris_ivbo(end_quad_vert));
		int const index_type(use_32_bit ? GL_UNSIGNED_INT : GL_UNSIGNED_SHORT);
		unsigned const bytes_offset((use_32_bit ? sizeof(unsigned) : sizeof(unsigned short))*start_tri_vert);

		if (num_instances > 1) { // instanced drawing not supported/efficient on some cards, so only enabled it when needed
			glDrawElementsInstanced(GL_TRIANGLES, num_tri_verts, index_type, (void *)bytes_offset, num_instances);
		}
		else {
			glDrawElements(GL_TRIANGLES, num_tri_verts, index_type, (void *)bytes_offset);
		}
		bind_vbo(0, 1);
	}
};

quad_ix_buffer_t quad_ix_buffer; // singleton

void clear_quad_ix_buffer_context() {quad_ix_buffer.free_context();}
void draw_quads_as_tris(unsigned num_quad_verts, unsigned start_quad_vert, unsigned num_instances) {
	quad_ix_buffer.draw_quads_as_tris(num_quad_verts, start_quad_vert, num_instances);
}
bool bind_quads_as_tris_ivbo(unsigned num_quad_verts) {return quad_ix_buffer.bind_quads_as_tris_ivbo(num_quad_verts);}


void convert_quad_ixs_to_tri_ixs(vector<unsigned> const &qixs, vector<unsigned> &tixs) { // what about 16-bit indices?

	tixs.resize(6*qixs.size()/4);

	for (unsigned i = 0, j = 0; i < qixs.size(); i += 4) { // step a quad at a time
		UNROLL_4X(tixs[j++] = qixs[i+i_];) // copy quad verts
		tixs[j++] = qixs[i+0];
		tixs[j++] = qixs[i+2];
	}
}


void lt_atten_manager_t::enable() {
	const char *lt_atten_uniform_strs[5] = {"light_atten", "refract_ix", "cube_bb", "sphere_center", "sphere_radius"};

	for (unsigned i = 0; i < 5; ++i) {
		ulocs[i] = shader.get_uniform_loc(lt_atten_uniform_strs[i]);
		//assert(ulocs[i] >= 0); // not all need to be available
	}
}

void lt_atten_manager_t::next_object(float light_atten, float refract_ix) {
	if (ulocs[0] >= 0 && light_atten != last_light_atten) {
		shader.set_uniform_float(ulocs[0], light_atten);
		last_light_atten = light_atten;
	}
	if (ulocs[1] >= 0 && refract_ix != last_refract_ix) {
		shader.set_uniform_float(ulocs[1], refract_ix);
		last_refract_ix = refract_ix;
	}
}

void lt_atten_manager_t::next_cube(float light_atten, float refract_ix, cube_t const &cube) {
	if (light_atten > 0.0) {shader.set_uniform_float_array(ulocs[2], (float const *)cube.d, 6);} // per-cube data
	next_object(light_atten, refract_ix);
}

void lt_atten_manager_t::next_sphere(float light_atten, float refract_ix, point const &pos, float radius) {
	if (light_atten > 0.0) {
		shader.set_uniform_vector3d(ulocs[3], pos);
		shader.set_uniform_float(ulocs[4], radius);
	}
	next_object(light_atten, refract_ix);
}


