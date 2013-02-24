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

void vert_norm::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 0, 1, 0);}
	set_vn_ptrs((force_stride ? force_stride : sizeof(vert_norm)), 0);
}

void vert_norm_comp_tc::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 0);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_comp_tc));
	set_vn_ptrs(stride, 1);
	glTexCoordPointer(2, GL_FLOAT, stride, (void *)(sizeof(vert_norm_comp)));
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

void vert_norm_comp_tc_color::set_vbo_arrays(unsigned force_stride, bool set_state) {
	if (set_state) {set_array_client_state(1, 1, 1, 1);}
	unsigned const stride(force_stride ? force_stride : sizeof(vert_norm_comp_tc_color));
	set_vn_ptrs(stride, 1);
	glTexCoordPointer(2, GL_FLOAT,         stride, (void *)(sizeof(vert_norm_comp)));
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, (void *)(sizeof(vert_norm_comp_tc)));
}


void vert_wrap_t::set_state() const { // typically called on element 0
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 0, 0);
	glVertexPointer(3, GL_FLOAT, stride, &v);
}

void vert_norm_comp::set_state() const { // typically called on element 0
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 1, 0);
	glVertexPointer(3, GL_FLOAT, stride, &v);
	glNormalPointer(GL_BYTE,     stride, &n);
}

void vert_norm_tc::set_state() const {
	set_array_client_state(1, 1, 1, 0);
	unsigned const stride(sizeof(*this));
	glVertexPointer  (3, GL_FLOAT, stride, &v);
	glNormalPointer  (   GL_FLOAT, stride, &n);
	glTexCoordPointer(2, GL_FLOAT, stride, &t);
}

void vert_norm_tc_color::set_state() const { // typically called on element 0
	set_array_client_state(1, 1, 1, 1);
	unsigned const stride(sizeof(*this));
	glVertexPointer  (3, GL_FLOAT,         stride, &v);
	glNormalPointer  (   GL_FLOAT,         stride, &n);
	glTexCoordPointer(2, GL_FLOAT,         stride, &t);
	glColorPointer   (4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_color::set_state() const { // typically called on element 0
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 0, 1);
	glVertexPointer(3, GL_FLOAT,         stride, &v);
	glColorPointer (4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_norm_color::set_state() const { // typically called on element 0
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 1, 1);
	glVertexPointer(3, GL_FLOAT,         stride, &v);
	glNormalPointer(   GL_FLOAT,         stride, &n);
	glColorPointer (4, GL_UNSIGNED_BYTE, stride, &c);
}

void vert_norm_comp_color::set_state() const { // typically called on element 0
	unsigned const stride(sizeof(*this));
	set_array_client_state(1, 0, 1, 1);
	glVertexPointer(3, GL_FLOAT,        stride, &v);
	glNormalPointer(GL_BYTE,            stride, &n);
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
void vbo_block_manager_t<vert_type_t>::add_points(vector<typename vert_type_t::non_color_class> const &p, colorRGBA const &color) {

	assert(!p.empty());
	color_wrapper cw;
	cw.set_c4(color);
		
	for (vector<vert_type_t::non_color_class>::const_iterator i = p.begin(); i != p.end(); ++i) {
		pts.push_back(vert_type_t(*i, cw));
	}
}

template< typename vert_type_t >
unsigned vbo_block_manager_t<vert_type_t>::add_points_with_offset(vector<typename vert_type_t::non_color_class> const &p, colorRGBA const &color) {

	add_points(p, color);
	if (offsets.empty()) {offsets.push_back(0);} // start at 0
	unsigned const next_ix(offsets.size() - 1);
	offsets.push_back(pts.size()); // range will be [start_ix, start_ix+p.size()]
	return next_ix;
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::render_range(int gl_type, unsigned six, unsigned eix) const {

	assert(six < eix && eix < offsets.size());
	glDrawArrays(gl_type, offsets[six], offsets[eix]-offsets[six]);
}


template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::draw_no_vbos(int gl_type) const { // unused

	if (pts.empty()) return;
	pts.front().set_state();
	glDrawArrays(gl_type, 0, (unsigned)pts.size());
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
void vbo_block_manager_t<vert_type_t>::begin_render(bool color_mat) const {

	if (!has_data()) return;
	if (color_mat) {glEnable(GL_COLOR_MATERIAL);}
	set_color(BLACK);
	assert(vbo);
	bind_vbo(vbo);
	vert_type_t::set_vbo_arrays();
}

template< typename vert_type_t >
void vbo_block_manager_t<vert_type_t>::end_render() const {

	glDisable(GL_COLOR_MATERIAL);
	bind_vbo(0);
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
	
	delete_vbo(vbo);
	vbo = 0;
}

// explicit template instantiations
template class vbo_block_manager_t<vert_color>;
template class vbo_block_manager_t<vert_norm_comp_color>;
template class vbo_block_manager_t<vert_norm_tc_color  >;



