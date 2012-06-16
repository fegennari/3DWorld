// 3D World - Drawing Utility Classes
// by Frank Gennari
// 6/15/12
#include "draw_utils.h"
#include "gl_ext_arb.h"



template class pt_line_drawer_t<color_wrapper      >;
template class pt_line_drawer_t<color_wrapper_float>;


template<typename cwt> void pt_line_drawer_t<cwt>::add_textured_pt(point const &v, colorRGBA c, int tid) {

	if (tid >= 0) c = c.modulate_with(texture_color(tid));
	vector3d const view_dir(get_camera_pos(), v);
	add_pt(v, view_dir, c);
}

template<typename cwt> void pt_line_drawer_t<cwt>::add_textured_line(point const &v1, point const &v2, colorRGBA c, int tid) {

	if (tid >= 0) c = c.modulate_with(texture_color(tid));
	vector3d view_dir(get_camera_pos(), (v1 + v2)*0.5);
	orthogonalize_dir(view_dir, (v2 - v1), view_dir, 0);
	add_line(v1, view_dir, c, v2, view_dir, c);
}

template<typename cwt> void pt_line_drawer_t<cwt>::vnc_cont::draw(int type) const {
	
	if (empty()) return; // nothing to do
	glVertexPointer(3, GL_FLOAT,     sizeof(vnc), &(front().v));
	glNormalPointer(   GL_FLOAT,     sizeof(vnc), &(front().n));
	glColorPointer( 4, cwt::gl_type, sizeof(vnc), &(front().c));
	glDrawArrays(type, 0, (unsigned)size());
}

template<typename cwt> void pt_line_drawer_t<cwt>::draw() const {
		
	if (points.empty() && lines.empty()) return;
	GLboolean const col_mat_en(glIsEnabled(GL_COLOR_MATERIAL));
	assert(!(lines.size() & 1));
	assert((triangles.size() % 3) == 0);
	if (!col_mat_en) glEnable(GL_COLOR_MATERIAL);
	set_array_client_state(1, 0, 1, 1);
	points.draw(GL_POINTS);
	lines.draw(GL_LINES);
	triangles.draw(GL_TRIANGLES);
	if (!col_mat_en) glDisable(GL_COLOR_MATERIAL);
	//cout << "mem: " << get_mem() << endl;
}


void quad_batch_draw::add_quad_vect(vector<vert_norm> const &points, colorRGBA const &color) {
	
	assert(!(points.size() & 3)); // must be a multiple of 4
	float const tcx[4] = {0,1,1,0}, tcy[4] = {0,0,1,1}; // 00 10 11 01
	color_wrapper cw;
	cw.set_c3(color);

	for (unsigned i = 0; i < points.size(); ++i) {
		verts.push_back(vert_norm_tc_color(points[i].v, points[i].n, tcx[i&3], tcy[i&3], cw.c));
	}
	unsigned const batch_size(4096);
	if (size() > batch_size) draw_and_clear();
}

void quad_batch_draw::draw() const {
	
	if (verts.empty()) return;
	assert(!(verts.size() & 3)); // must be a multiple of 4
	verts.front().set_state();
	glDrawArrays(GL_QUADS, 0, (unsigned)size());
}


unsigned vbo_quad_block_manager_t::add_points(vector<vert_norm> const &p) {

	assert(!p.empty());
	assert((p.size()&3) == 0); // must be quads
	unsigned const num_quads(p.size()/4), start_ix(pts.size());
		
	for (vector<vert_norm>::const_iterator i = p.begin(); i != p.end(); ++i) {
		pts.push_back(vert_norm_tc(*i));
	}
	gen_quad_tex_coords(pts[start_ix].t, num_quads, sizeof(vert_norm_tc)/sizeof(float));
	assert(!offsets.empty());
	unsigned const next_ix(offsets.size() - 1);
	offsets.push_back(pts.size()); // range will be [start_ix, start_ix+p.size()]
	return next_ix;
}

void vbo_quad_block_manager_t::render_range(unsigned six, unsigned eix) const {

	assert(six < eix && eix < offsets.size());
	assert(offsets[eix] <= pts.size());
	glDrawArrays(GL_QUADS, offsets[six], offsets[eix]-offsets[six]);
}

void vbo_quad_block_manager_t::upload() {

	if (vbo || empty()) return; // already uploaded or empty
	vbo = create_vbo();
	bind_vbo(vbo);
	upload_vbo_data(&pts.front(), pts.size()*sizeof(vert_norm_tc));
	bind_vbo(0);
}

void vbo_quad_block_manager_t::begin_render() const {

	if (empty()) return;
	assert(vbo);
	bind_vbo(vbo);
	vert_norm_tc::set_vbo_arrays();
}

void vbo_quad_block_manager_t::end_render() const {

	bind_vbo(0);
}

void vbo_quad_block_manager_t::clear_vbo() {
	
	delete_vbo(vbo);
	vbo = 0;
}

void vbo_quad_block_manager_t::clear() {

	pts.clear();
	offsets.clear();
	offsets.push_back(0); // start at 0
	clear_vbo();
}



