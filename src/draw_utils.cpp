// 3D World - Drawing Utility Classes
// by Frank Gennari
// 6/15/12
#include "draw_utils.h"
#include "function_registry.h"
#include "gl_ext_arb.h"


void set_array_client_state(bool va, bool tca, bool na, bool ca) {

	bool const enables[4] = {va, tca, na, ca};
	int  const arrays [4] = {GL_VERTEX_ARRAY, GL_TEXTURE_COORD_ARRAY, GL_NORMAL_ARRAY, GL_COLOR_ARRAY};

	for (unsigned i = 0; i < 4; ++i) {
		if (enables[i]) {
			glEnableClientState(arrays[i]);
		}
		else {
			glDisableClientState(arrays[i]);
		}
	}
}


void set_vn_ptrs(unsigned stride, bool comp) {
	glVertexPointer(3, GL_FLOAT, stride, (void *)(0));
	glNormalPointer((comp ? GL_BYTE : GL_FLOAT), stride, (void *)(sizeof(point)));
}

void vert_wrap_t::set_vbo_arrays(bool set_state) {
	if (set_state) {set_array_client_state(1, 0, 0, 0);}
	glVertexPointer(3, GL_FLOAT, sizeof(point), (void *)(0));
}

void vert_tc_t::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 0, 0);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_tc_t));
	glVertexPointer(3, GL_FLOAT, sizeof(point), (void *)(0));
	glTexCoordPointer(2, GL_FLOAT, stride, (void *)(sizeof(vert_wrap_t)));
}

void vert_norm::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 0, 1, 0);}
	set_vn_ptrs((force_stride ? force_stride : sizeof(vert_norm)), 0);
}

void vert_norm_comp::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 0, 1, 0);}
	set_vn_ptrs(force_stride ? force_stride : sizeof(vert_norm_comp), 1);
}

void vert_norm_comp_tc::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 0);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_comp_tc));
	set_vn_ptrs(stride, 1);
	glTexCoordPointer(2, GL_FLOAT, stride, (void *)(sizeof(vert_norm_comp)));
}

void vert_norm_comp_tc_comp::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 0);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_comp_tc_comp));
	set_vn_ptrs(stride, 1);
	glTexCoordPointer(2, GL_SHORT, stride, (void *)(sizeof(vert_norm_comp)));
}

void vert_norm_tc::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 0);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_tc));
	set_vn_ptrs(stride, 0);
	glTexCoordPointer(2, GL_FLOAT, stride, (void *)(sizeof(vert_norm)));
}

void vert_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 0, 0, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_color));
	glVertexPointer(3, GL_FLOAT, stride, (void *)(0));
	glColorPointer(4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(point)));
}

void vert_norm_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 0, 1, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_color));
	set_vn_ptrs(stride, 0);
	glColorPointer(4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(vert_norm)));
}

void vert_norm_comp_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 0, 1, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_comp_color));
	set_vn_ptrs(stride, 1);
	glColorPointer(4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(vert_norm_comp)));
}

void vert_norm_tc_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_tc_color));
	set_vn_ptrs(stride, 0);
	glTexCoordPointer(2, GL_FLOAT,         stride, (void *)(sizeof(vert_norm)));
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(vert_norm_tc)));
}

void vert_tc_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 0, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_tc_color));
	glVertexPointer  (3, GL_FLOAT, stride, (void *)(0));
	glTexCoordPointer(2, GL_FLOAT,         stride, (void *)(sizeof(point)));
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(vert_tc_t)));
}

void vert_norm_comp_tc_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_comp_tc_color));
	set_vn_ptrs(stride, 1);
	glTexCoordPointer(2, GL_FLOAT,         stride, (void *)(sizeof(vert_norm_comp)));
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(vert_norm_comp_tc)));
}

void vert_norm_comp_tc_comp_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_comp_tc_comp_color));
	set_vn_ptrs(stride, 1);
	glTexCoordPointer(2, GL_SHORT,         stride, (void *)(sizeof(vert_norm_comp)));
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(vert_norm_comp_tc_comp)));
}


// set_state() functions for all vertex classes - typically called on element 0
void vert_wrap_t::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 0, 0);
	glVertexPointer(3, GL_FLOAT, stride, &v);
}

void vert_tc_t::set_state() const {
	set_array_client_state(1, 1, 0, 0);
	unsigned const stride(sizeof(*this));
	glVertexPointer  (3, GL_FLOAT, stride, &v);
	glTexCoordPointer(2, GL_FLOAT, stride, &t);
}

void vert_norm::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 1, 0);
	glVertexPointer(3, GL_FLOAT, stride, &v);
	glNormalPointer(   GL_FLOAT, stride, &n);
}

void vert_norm_comp::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 1, 0);
	glVertexPointer(3, GL_FLOAT, stride, &v);
	glNormalPointer(GL_BYTE,     stride, &n);
}

void vert_norm_comp_tc::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 1, 1, 0);
	glVertexPointer(3, GL_FLOAT,   stride, &v);
	glNormalPointer(GL_BYTE,       stride, &n);
	glTexCoordPointer(2, GL_FLOAT, stride, &t);
}

void vert_norm_comp_tc_comp::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 1, 1, 0);
	glVertexPointer(3, GL_FLOAT,   stride, &v);
	glNormalPointer(GL_BYTE,       stride, &n);
	glTexCoordPointer(2, GL_SHORT, stride, &t);
}

void vert_norm_tc::set_state() const {
	set_array_client_state(1, 1, 1, 0);
	unsigned const stride(sizeof(*this));
	glVertexPointer  (3, GL_FLOAT, stride, &v);
	glNormalPointer  (   GL_FLOAT, stride, &n);
	glTexCoordPointer(2, GL_FLOAT, stride, &t);
}

void vert_norm_tc_color::set_state() const {
	set_array_client_state(1, 1, 1, 1);
	unsigned const stride(sizeof(*this));
	glVertexPointer  (3, GL_FLOAT,         stride, &v);
	glNormalPointer  (   GL_FLOAT,         stride, &n);
	glTexCoordPointer(2, GL_FLOAT,         stride, &t);
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_tc_color::set_state() const {
	set_array_client_state(1, 1, 0, 1);
	unsigned const stride(sizeof(*this));
	glVertexPointer  (3, GL_FLOAT,         stride, &v);
	glTexCoordPointer(2, GL_FLOAT,         stride, &t);
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_color::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 0, 1);
	glVertexPointer(3, GL_FLOAT,         stride, &v);
	glColorPointer (4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_norm_color::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 1, 1);
	glVertexPointer(3, GL_FLOAT,         stride, &v);
	glNormalPointer(   GL_FLOAT,         stride, &n);
	glColorPointer (4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_norm_comp_color::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 1, 1);
	glVertexPointer(3, GL_FLOAT,        stride, &v);
	glNormalPointer(GL_BYTE,            stride, &n);
	glColorPointer(4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_norm_comp_tc_color::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 1, 1, 1);
	glVertexPointer(3, GL_FLOAT,        stride, &v);
	glNormalPointer(GL_BYTE,            stride, &n);
	glTexCoordPointer(2, GL_FLOAT,      stride, &t);
	glColorPointer(4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_norm_comp_tc_comp_color::set_state() const {
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 1, 1, 1);
	glVertexPointer(3, GL_FLOAT,        stride, &v);
	glNormalPointer(GL_BYTE,            stride, &n);
	glTexCoordPointer(2, GL_SHORT,      stride, &t);
	glColorPointer(4, GL_UNSIGNED_BYTE, stride, &c);
}


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
		
	if (empty()) return;
	assert(!(lines.size() & 1));
	assert((triangles.size() % 3) == 0);
	set_array_client_state(1, 0, 1, 1);
	points.draw(GL_POINTS);
	lines.draw(GL_LINES);
	triangles.draw(GL_TRIANGLES);
	//cout << "mem: " << get_mem() << endl;
}


void pt_line_drawer_no_lighting_t::draw() const {

	assert(!(lines.size() & 1));
	draw_verts(points, GL_POINTS);
	draw_verts(lines,  GL_LINES);
}


void quad_batch_draw::add_quad_pts(point const pts[4], colorRGBA const &c, vector3d const &n, tex_range_t const &tr) {

	float const t[4][2] = {{tr.x1,tr.y1}, {tr.x2,tr.y1}, {tr.x2,tr.y2}, {tr.x1,tr.y2}};
	unsigned const v[6] = {0,2,1, 0,3,2};
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
	draw_and_clear(GL_TRIANGLES);
	glDepthMask(GL_TRUE);
}


template<typename T> void indexed_mesh_draw<T>::clear() {

	verts.clear();
	indices.clear();
	nx = ny = 0;
}

template<typename T> void indexed_mesh_draw<T>::init(unsigned nx_, unsigned ny_) {

	if (nx == nx_ && ny == ny_) return; // already setup
	nx = nx_; ny = ny_;
	assert(nx > 0 && ny > 0);
	verts.resize((nx+1)*(ny+1));
	indices.resize(4*nx*ny);
		
	for (unsigned y = 0; y < ny; ++y) {
		for (unsigned x = 0; x < nx; ++x) {
			unsigned const iix(4*(y*nx + x)), vix(get_vert_ix(x, y));
			indices[iix+0] = vix;
			indices[iix+1] = vix + 1;
			indices[iix+2] = vix + (nx+1) + 1;
			indices[iix+3] = vix + (nx+1);
		}
	}
}

template<typename T> void indexed_mesh_draw<T>::render() const {

	if (verts.empty()) return;
	verts.front().set_state();
	glDrawRangeElements(GL_QUADS, 0, verts.size(), indices.size(), GL_UNSIGNED_INT, &indices.front());
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
void vbo_block_manager_t<vert_type_t>::render_range(int gl_type, unsigned six, unsigned eix, unsigned num_instances) const {

	assert(six < eix && eix < offsets.size());

	if (num_instances > 0) { // instanced rendering
		glDrawArraysInstanced(gl_type, offsets[six], offsets[eix]-offsets[six], num_instances);
	}
	else { // normal rendering
		glDrawArrays(gl_type, offsets[six], offsets[eix]-offsets[six]);
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
	bind_cur_vbo();
	upload_vbo_sub_data(&update_verts.front(), start*sizeof(vert_type_t), update_size*sizeof(vert_type_t));
	bind_vbo(0);
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::begin_render() const {

	if (!has_data()) return;
	set_color(BLACK);
	bind_cur_vbo();
	vert_type_t::set_vbo_arrays();
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::end_render() const {
	bind_vbo(0);
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::bind_cur_vbo() const {
	assert(vbo);
	bind_vbo(vbo);
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::clear() {

	clear_points();
	temp_points.clear();
	offsets.clear();
	clear_vbo();
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::clear_vbo() {
	delete_and_zero_vbo(vbo);
}

// explicit template instantiations
template class vbo_block_manager_t<vert_color>;
template class vbo_block_manager_t<vert_norm_comp_color>;
template class vbo_block_manager_t<vert_norm_tc_color>;
template class vbo_block_manager_t<vert_norm_tc>;



