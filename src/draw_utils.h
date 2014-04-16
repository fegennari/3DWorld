// 3D World - Drawing Utility Classes
// by Frank Gennari
// 6/15/12

#ifndef _DRAW_UTILS_H_
#define _DRAW_UTILS_H_

#include "3DWorld.h"


class pt_line_drawer {

	vector<vert_norm_color> points, lines;

public:
	void clear() {
		points.resize(0);
		lines.resize(0);
	}
	void free_mem() {
		points.swap(vector<vert_norm_color>());
		lines.swap(vector<vert_norm_color>());
	}
	void add_pt(point const &v, vector3d const &n, colorRGBA const &c) {
		points.push_back(vert_norm_color(v, n, c));
	}
	void add_line(point const &v1, vector3d const &n1, colorRGBA const &c1, point const &v2, vector3d const &n2, colorRGBA const &c2) {
		lines.push_back(vert_norm_color(v1, n1, c1));
		lines.push_back(vert_norm_color(v2, n2, c2));
	}
	void add_textured_pt(point const &v, colorRGBA c, int tid);
	void add_textured_line(point const &v1, point const &v2, colorRGBA c, int tid);
	void draw() const;
	void draw_and_clear() {draw(); clear();}
	size_t get_mem() const {return (points.capacity() + lines.capacity())*sizeof(vert_norm_color);}
	bool empty() const {return (points.empty() && lines.empty());}
};


class pt_line_drawer_no_lighting_t { // uses vbo for static points/lines; used for universe mode stars

	vector<vert_color> points, lines;
	unsigned vbo[2];

public:
	pt_line_drawer_no_lighting_t() {vbo[0] = vbo[1] = 0;}
	void clear() {points.resize(0); lines.resize(0); free_vbo();}
	void add_pt(point const &v, colorRGBA const &c) {
		points.push_back(vert_color(v, c));
	}
	void add_line(point const &v1, colorRGBA const &c1, point const &v2, colorRGBA const &c2) {
		lines.push_back(vert_color(v1, c1));
		lines.push_back(vert_color(v2, c2));
	}
	void draw() const;
	void draw_vbo();
	void draw_and_clear() {draw(); clear();}
	bool empty() const {return (points.empty() && lines.empty());}
	void free_vbo();
};


// Note: we may want versions with and without normals
template<class vert_type_t> class point_sprite_drawer_t {

	vector<vert_type_t> points;
	vector<float> sizes;

public:
	void clear() {points.resize(0); sizes.resize(0);}
	void reserve_pts(unsigned sz) {points.reserve(sz);}
	void add_pt(vert_type_t const &v) {points.push_back(v);} // constant size version
	void add_pt(vert_type_t const &v, float s) {points.push_back(v); sizes.push_back(s);} // variable size version
	void draw(int tid, float const_point_size=0.0, bool enable_lighting=0) const;
	void draw_and_clear(int tid, float const_point_size=0.0, bool enable_lighting=0) {draw(tid, const_point_size, enable_lighting); clear();}
	bool empty() const {return points.empty();}
};

typedef point_sprite_drawer_t<vert_color     > point_sprite_drawer;
typedef point_sprite_drawer_t<vert_norm_color> point_sprite_drawer_norm;


struct quad_batch_draw { // Note: might want an indexed version of this

	vector<vert_norm_tc_color> verts;

	void add_quad_pts(point const pts[4], colorRGBA const &c, vector3d const &n=plus_z, tex_range_t const &tr=tex_range_t());
	void add_quad_dirs(point const &pos, vector3d const &dx, vector3d const &dy, colorRGBA const &c, vector3d const &n=plus_z, tex_range_t const &tr=tex_range_t());
	void add_xlated_billboard(point const &pos, point const &xlate, point const &viewer, vector3d const &up_dir, colorRGBA const &c,
		float xsize, float ysize, tex_range_t const &tr=tex_range_t(), bool minimize_fill=0);
	void add_billboard(point const &pos, point const &viewer, vector3d const &up_dir, colorRGBA const &c,
		float xsize, float ysize, tex_range_t const &tr=tex_range_t(), bool minimize_fill=0) {
			add_xlated_billboard(pos, pos, viewer, up_dir, c, xsize, ysize, tr, minimize_fill);
	}
	bool empty() const {return verts.empty();}
	void clear() {verts.clear();}
	void add_animated_billboard(point const &pos, point const &viewer, vector3d const &up_dir, colorRGBA const &c, float xsize, float ysize, float timescale);
	void draw() const {draw_verts(verts, GL_TRIANGLES);}
	void draw_and_clear() {draw(); clear();}
	void draw_and_clear_quads() {draw_quad_verts_as_tris(verts); clear();}
	void draw_as_flares_and_clear(int flare_tex=BLUR_TEX);
};


class line_tquad_draw_t {

	vector<vert_tc_color> verts;

public:
	void add_line_as_tris(point const &p1, point const &p2, float w1, float w2, colorRGBA const &color1, colorRGBA const &color2,
		point const* const prev=NULL, point const *const next=NULL, bool make_global=0);
	void draw() const;
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
	unsigned get_offset_for_last_points_added();
	unsigned add_points_with_offset(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color) {
		add_points(p, npts, color);
		return get_offset_for_last_points_added();
	}
	unsigned add_points_with_offset(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color) {return add_points_with_offset(&v.front(), v.size(), color);}
	void render_range(unsigned six, unsigned eix, unsigned num_instances=1) const;
	void render_all(unsigned num_instances=1) const {if (has_data()) {render_range(0, offsets.size()-1, num_instances);}}
	bool upload();
	void update_range(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color, unsigned six, unsigned eix);
	void update_range(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color, unsigned six, unsigned eix) {update_range(&v.front(), v.size(), color, six, eix);}
	void begin_render() const;
	void end_render() const;
	void bind_cur_vbo() const;
	void clear_points() {pts.swap(vector<vert_type_t>());}
	vector<vert_type_t> &get_pts_vector_for_adding() {return pts;}
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

