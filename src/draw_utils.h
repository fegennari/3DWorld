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
	void add_quad(point const v[4], vector3d const &n, colorRGBA const &c) {
		add_triangle(v[0], v[2], v[1], n, c);
		add_triangle(v[0], v[3], v[2], n, c);
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


class pt_line_drawer_no_lighting_t {

	vector<vert_color> points, lines;

public:
	void clear() {points.resize(0); lines.resize(0);}
	void add_pt(point const &v, colorRGBA const &c) {
		points.push_back(vert_color(v, c));
	}
	void add_line(point const &v1, colorRGBA const &c1, point const &v2, colorRGBA const &c2) {
		lines.push_back(vert_color(v1, c1));
		lines.push_back(vert_color(v2, c2));
	}
	void draw() const;
	void draw_and_clear() {draw(); clear();}
	bool empty() const {return (points.empty() && lines.empty());}
};


struct quad_batch_draw { // Note: might want an indexed version of this

	vector<vert_norm_tc_color> verts;

	void add_quad_pts(point const pts[4], colorRGBA const &c, vector3d const &n=plus_z, float tx1=0.0, float ty1=0.0, float tx2=1.0, float ty2=1.0);
	void add_quad_dirs(point const &pos, vector3d const &dx, vector3d const &dy, colorRGBA const &c, vector3d const &n=plus_z,
		float tx1=0.0, float ty1=0.0, float tx2=1.0, float ty2=1.0);
	void add_xlated_billboard(point const &pos, point const &xlate, point const &viewer, vector3d const &up_dir, colorRGBA const &c,
		float xsize, float ysize, float tx1=0.0, float ty1=0.0, float tx2=1.0, float ty2=1.0, bool minimize_fill=0);
	void add_billboard(point const &pos, point const &viewer, vector3d const &up_dir, colorRGBA const &c,
		float xsize, float ysize, float tx1=0.0, float ty1=0.0, float tx2=1.0, float ty2=1.0, bool minimize_fill=0) {
			add_xlated_billboard(pos, pos, viewer, up_dir, c, xsize, ysize, tx1, ty1, tx2, ty2, minimize_fill);
	}
	void add_animated_billboard(point const &pos, point const &viewer, vector3d const &up_dir, colorRGBA const &c, float xsize, float ysize, float timescale);
	void draw(int prim_type=GL_TRIANGLES) const {draw_verts(verts, prim_type);} // GL_QUADS or GL_TRIANGLES
	void draw_and_clear(int prim_type=GL_TRIANGLES) {draw(prim_type); verts.clear();}
	void draw_as_flares_and_clear(int flare_tex=BLUR_TEX);
};


class line_tquad_draw_t {

	vector<vert_tc_color> verts;

public:
	void add_line_tquad  (point const &p1, point const &p2, float w1, float w2, colorRGBA const &color1, colorRGBA const &color2,
		point const* const prev=NULL, point const *const next=NULL);
	void add_line_as_tris(point const &p1, point const &p2, float w1, float w2, colorRGBA const &color1, colorRGBA const &color2,
		point const* const prev=NULL, point const *const next=NULL, bool make_global=0);
	void draw(int prim_type) const;
};


template<typename T> class indexed_mesh_draw { // quads

protected:
	unsigned nx, ny; // in quads
	vector<T> verts;
	vector<unsigned> indices;

public:
	indexed_mesh_draw() : nx(0), ny(0) {}
	void clear();
	void init(unsigned nx_, unsigned ny_);

	unsigned get_vert_ix(unsigned x, unsigned y) const {
		assert(!verts.empty() && x <= nx && y <= ny);
		return (y*(nx+1) + x);
	}
	void set_vert(unsigned x, unsigned y, T const &v) {verts[get_vert_ix(x, y)] = v;}
	void render() const;
	void render_z_plane(float x1, float y1, float x2, float y2, float zval, unsigned nx_, unsigned ny_);
};


template< typename vert_type_t > class vbo_block_manager_t {

	vector<vert_type_t> pts;
	vector<unsigned> offsets;
	unsigned vbo;

	bool has_data() const {return (!pts.empty() || offsets.size() > 1);}
	void add_points_int(vector<vert_type_t> &dest, typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color);

public:
	vector<vert_norm> temp_points;

	vbo_block_manager_t() : vbo(0) {clear();}
	// can't free in the destructor because the gl context may be destroyed before this point
	//~vbo_block_manager_t() {clear_vbo();}
	bool is_uploaded() const {return (vbo != 0);}
	void reserve_pts(unsigned num) {assert(pts.empty()); pts.reserve(num);}
	void add_points(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color) {add_points_int(pts, p, npts, color);}
	void add_points(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color) {add_points_int(pts, &v.front(), v.size(), color);}
	unsigned add_points_with_offset(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color);
	unsigned add_points_with_offset(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color) {return add_points_with_offset(&v.front(), v.size(), color);}
	void render_range(int gl_type, unsigned six, unsigned eix, unsigned num_instances=0) const;
	void render_all(int gl_type, unsigned num_instances=0) const {if (has_data()) {render_range(gl_type, 0, offsets.size()-1, num_instances);}}
	void draw_no_vbos(int gl_type) const {draw_verts(pts, gl_type);} // unused
	bool upload();
	void update_range(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color, unsigned six, unsigned eix);
	void update_range(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color, unsigned six, unsigned eix) {update_range(&v.front(), v.size(), color, six, eix);}
	void begin_render(bool color_mat) const;
	void end_render() const;
	void bind_cur_vbo() const;
	void clear_points() {pts.swap(vector<vert_type_t>());}
	void clear_vbo();
	void clear();
	void upload_and_clear_points() {upload(); clear_points();}
	unsigned get_gpu_mem() const {return ((vbo && has_data()) ? offsets.back()*sizeof(vert_type_t) : 0);}
};


typedef vbo_block_manager_t<vert_color> vbo_vc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_comp_color> vbo_vnc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_tc_color  > vbo_vntc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_tc        > vbo_vnt_block_manager_t;


#endif // _DRAW_UTILS_H_

