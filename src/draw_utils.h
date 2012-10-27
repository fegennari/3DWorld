// 3D World - Drawing Utility Classes
// by Frank Gennari
// 6/15/12

#ifndef _DRAW_UTILS_H_
#define _DRAW_UTILS_H_

#include "3DWorld.h"


template<typename cwt> class pt_line_drawer_t { // and triangles too!

	struct vnc : public vert_norm, public cwt { // size = 28
		vnc() {}
		vnc(point const &v_, vector3d const &n_, colorRGBA const &c_) : vert_norm(v_, n_) {set_c4(c_);}
	};

	struct vnc_cont : public vector<vnc> {
		void draw(int type) const;
	};

	vnc_cont points, lines, triangles;

public:
	void clear() {
		points.resize(0);
		lines.resize(0);
		triangles.resize(0);
	}
	void free_mem() {
		points.swap(vnc_cont());
		lines.swap(vnc_cont());
		triangles.swap(vnc_cont());
	}
	void add_pt(point const &v, vector3d const &n, colorRGBA const &c) {
		points.push_back(vnc(v, n, c));
	}
	void add_line(point const &v1, vector3d const &n1, colorRGBA const &c1, point const &v2, vector3d const &n2, colorRGBA const &c2) {
		lines.push_back(vnc(v1, n1, c1));
		lines.push_back(vnc(v2, n2, c2));
	}
	void add_triangle(point const &v1, point const &v2, point const &v3, vector3d const &n, colorRGBA const &c) {
		points.push_back(vnc(v1, n, c));
		points.push_back(vnc(v2, n, c));
		points.push_back(vnc(v3, n, c));
	}
	void add_textured_pt(point const &v, colorRGBA c, int tid);
	void add_textured_line(point const &v1, point const &v2, colorRGBA c, int tid);
	void draw() const;
	void draw_and_clear() {draw(); clear();}
	size_t get_mem() const {return (points.capacity() + lines.capacity() + triangles.capacity())*sizeof(vnc);}
	bool empty() const {return (points.empty() && lines.empty() && triangles.empty());}
};


typedef pt_line_drawer_t<color_wrapper      > pt_line_drawer;
typedef pt_line_drawer_t<color_wrapper_float> pt_line_drawer_hdr;


template<typename T> class indexed_mesh_draw { // quads

	unsigned nx, ny; // in quads
	vector<T> verts;
	vector<unsigned> indices;

public:
	indexed_mesh_draw() : nx(0), ny(0) {}
	void clear();
	void init(unsigned nx_, unsigned ny_);

	void set_vert(unsigned x, unsigned y, T const &v) {
		assert(x <= nx && y <= ny);
		verts[y*(nx+1) + x] = v;
	}
	void render() const;
	void render_z_plane(float x1, float y1, float x2, float y2, float zval, unsigned nx_, unsigned ny_);
};


template< typename vert_type_t > class vbo_block_manager_t {

	vector<vert_type_t> pts;
	vector<unsigned> offsets;
	unsigned vbo;

	bool has_data() const {return (!pts.empty() || offsets.size() > 1);}

public:
	vector<vert_norm> temp_points;

	vbo_block_manager_t() : vbo(0) {clear();}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~vbo_block_manager_t() {clear_vbo();}
	bool is_uploaded() const {return (vbo != 0);}
	void reserve_pts(unsigned num) {assert(pts.empty()); pts.reserve(num);}
	void add_points(vector<typename vert_type_t::non_color_class> const &p, colorRGBA const &color);
	unsigned add_points_with_offset(vector<typename vert_type_t::non_color_class> const &p, colorRGBA const &color);
	void render_range(int gl_type, unsigned six, unsigned eix) const;
	void render_all(int gl_type) const {if (has_data()) {render_range(gl_type, 0, offsets.size()-1);}}
	void draw_no_vbos(int gl_type) const;
	bool upload();
	void begin_render(bool color_mat) const;
	void end_render() const;
	void clear_points() {pts.swap(vector<vert_type_t>());}
	void clear_vbo();
	void clear();
	void upload_and_clear_points() {upload(); clear_points();}
	unsigned get_gpu_mem() const {return ((vbo && has_data()) ? offsets.back()*sizeof(vert_type_t) : 0);}
};


typedef vbo_block_manager_t<vert_color> vbo_vc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_comp_color> vbo_vnc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_tc_color  > vbo_vntc_block_manager_t;


#endif // _DRAW_UTILS_H_

