// 3D World - Drawing Utility Classes
// by Frank Gennari
// 6/15/12
#pragma once

#include "3DWorld.h"
#include "gl_ext_arb.h"


template<class vert_type_t> struct sized_vert_t : public vert_type_t { // size = 20-32
	float size;

	sized_vert_t() : size(0.0f) {}
	sized_vert_t(vert_type_t const &v, float size_) : vert_type_t(v), size(size_) {}
	static void set_vbo_arrays(bool set_state=1, void const *vbo_ptr_offset=NULL);
	static void set_size_attr(unsigned stride, void const *vbo_ptr_offset);
	static void unset_attrs();
};


class pt_line_drawer {

	vector<vert_norm_color> points, lines;

public:
	void clear() {
		points.resize(0);
		lines.resize(0);
	}
	void free_mem() {
		clear_cont(points);
		clear_cont(lines);
	}
	void add_pt(point const &v, vector3d const &n, colorRGBA const &c) {
		points.emplace_back(v, n, c);
	}
	void add_line(point const &v1, vector3d const &n1, colorRGBA const &c1, point const &v2, vector3d const &n2, colorRGBA const &c2) {
		lines.emplace_back(v1, n1, c1);
		lines.emplace_back(v2, n2, c2);
	}
	void add_textured_pt(point const &v, colorRGBA c, int tid);
	void add_textured_line(point const &v1, point const &v2, colorRGBA c, int tid);
	void draw() const;
	void draw_and_clear() {draw(); clear();}
	size_t get_mem() const {return (get_cont_mem_usage(points) + get_cont_mem_usage(lines));}
	bool empty() const {return (points.empty() && lines.empty());}
};


class pt_line_drawer_no_lighting_t { // uses vbo for static points/lines; used for universe mode stars

	vector<vert_color> points, lines;
	unsigned vbo[2]={};

public:
	void clear() {points.clear(); lines.clear(); free_vbo();}
	void add_pt  (point const &v,  colorRGBA const &c) {points.emplace_back(v, c);}
	void add_line(point const &v1, colorRGBA const &c1, point const &v2, colorRGBA const &c2) {
		lines.emplace_back(v1, c1);
		lines.emplace_back(v2, c2);
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

public:
	void clear() {points.clear();}
	void reserve_pts(unsigned sz) {points.reserve(sz);}
	void add_pt(vert_type_t const &v) {points.push_back(v);}
	void sort_back_to_front();
	void draw(int tid, float const_point_size=0.0, bool enable_lighting=0, bool use_geom_shader=0, float min_alpha=0.0) const;
	void draw_and_clear(int tid, float const_point_size=0.0, bool enable_lighting=0, bool use_geom_shader=0, float min_alpha=0.0) {
		draw(tid, const_point_size, enable_lighting, use_geom_shader, min_alpha); clear();
	}
	bool empty() const {return points.empty();}
};

typedef point_sprite_drawer_t<vert_color                   > point_sprite_drawer;
typedef point_sprite_drawer_t<sized_vert_t<vert_color>     > point_sprite_drawer_sized;
typedef point_sprite_drawer_t<vert_norm_color              > point_sprite_drawer_norm;
typedef point_sprite_drawer_t<sized_vert_t<vert_norm_color>> point_sprite_drawer_norm_sized;


struct quad_batch_draw { // Note: might want an indexed version of this

	vector<vert_norm_tc_color> verts;

	void add_quad_pts(point const pts[4], color_wrapper const &cw, vector3d const &n=plus_z, tex_range_t const &tr=tex_range_t());
	void add_quad_pts_vert_norms(vert_norm const pts[4], color_wrapper const &cw, tex_range_t const &tr=tex_range_t());
	void add_quad_pts_vert_norms(point const pts[4], vector3d const n[4], color_wrapper const &cw, tex_range_t const &tr=tex_range_t());
	void add_quad_dirs(point const &pos, vector3d const &dx, vector3d const &dy, colorRGBA const &c, vector3d const &n=plus_z, tex_range_t const &tr=tex_range_t());
	void add_quad_dirs_single_tri(point const &pos, vector3d const &dx, vector3d const &dy, colorRGBA const &c, vector3d const &n);
	void add_xlated_billboard(point const &pos, point const &xlate, point const &viewer, vector3d const &up_dir, colorRGBA const &c,
		float xsize, float ysize, tex_range_t const &tr=tex_range_t(), bool minimize_fill=0, vector3d const *const normal_=nullptr);
	void add_billboard(point const &pos, point const &viewer, vector3d const &up_dir, colorRGBA const &c,
		float xsize, float ysize, tex_range_t const &tr=tex_range_t(), bool minimize_fill=0, vector3d const *const normal_=nullptr) {
			add_xlated_billboard(pos, pos, viewer, up_dir, c, xsize, ysize, tr, minimize_fill, normal_);
	}
	bool empty() const {return verts.empty();}
	void clear() {verts.clear();}
	void add_animated_billboard(point const &pos, point const &viewer, vector3d const &up_dir, colorRGBA const &c, float xsize, float ysize, float timescale);
	void draw() const {draw_verts(verts, GL_TRIANGLES);}
	void draw_and_clear() {draw(); clear();}
	void draw_and_clear_quads() {draw_quad_verts_as_tris(verts); clear();}
	void draw_as_flares_and_clear(int flare_tex=BLUR_TEX);
	void add_quads(quad_batch_draw const &qbd) {verts.insert(verts.end(), qbd.verts.begin(), qbd.verts.end());}
};


class line_tquad_draw_t {

	vector<vert_tc_color> verts;

public:
	bool empty() const {return verts.empty();}
	size_t size() {return verts.size();}
	void clear() {verts.clear();}
	void reserve_verts(unsigned size) {verts.reserve(size);}
	void add_line_as_tris(point const &p1, point const &p2, float w1, float w2, colorRGBA const &color1, colorRGBA const &color2,
		point const* const prev=NULL, point const *const next=NULL, bool make_global=0, bool use_2d_tex_coords=0);
	void draw(float noise_scale=0.0) const;
	void draw_and_clear(float noise_scale=0.0) {draw(noise_scale); clear();}
	void draw_tri_verts() const {draw_verts(verts, GL_TRIANGLES);}
};


template<typename T> class indexed_mesh_draw { // quads

protected:
	unsigned ivbo, ivbo_size, nx, ny; // in quads
	vector<T> verts;

public:
	indexed_mesh_draw() : ivbo(0), ivbo_size(0), nx(0), ny(0) {}
	void clear();
	void init(unsigned nx_, unsigned ny_);

	unsigned get_vert_ix(unsigned x, unsigned y) const {
		assert(!verts.empty() && x <= nx && y <= ny);
		return (y*(nx+1) + x);
	}
	void set_vert(unsigned x, unsigned y, T const &v) {verts[get_vert_ix(x, y)] = v;}
	void render() const;
	void render_z_plane(float x1, float y1, float x2, float y2, float zval, unsigned nx_, unsigned ny_);
	void free_context() {delete_and_zero_vbo(ivbo);}
};


template< typename vert_type_t > class vbo_block_manager_t : public vbo_wrap_t {

	vector<vert_type_t> pts;
	vector<unsigned> offsets;
	int prim_type;

	bool has_data() const {return (!pts.empty() || offsets.size() > 1);}
	void add_points_int(vector<vert_type_t> &dest, typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color);
public:
	// can't free in the destructor because the gl context may be destroyed before this point
	vbo_block_manager_t() : prim_type(GL_QUADS) {} // default is quads
	//~vbo_block_manager_t() {clear_vbo();}
	void set_prim_type(int type) {prim_type = type;}
	bool is_uploaded() const {return (vbo != 0);}
	void reserve_pts(unsigned num) {assert(pts.empty()); pts.reserve(num);}
	void reserve_offsets(unsigned num) {assert(offsets.empty()); offsets.reserve(num);}
	void add_point(vert_type_t const &pt) {pts.push_back(pt);}
	void add_points(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color) {add_points_int(pts, p, npts, color);}
	void add_points(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color) {add_points_int(pts, &v.front(), v.size(), color);}
	unsigned get_offset_for_last_points_added();
	unsigned add_points_with_offset(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color) {
		add_points(p, npts, color);
		return get_offset_for_last_points_added();
	}
	unsigned add_points_with_offset(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color) {
		assert(!v.empty());
		return add_points_with_offset(&v.front(), v.size(), color);
	}
	unsigned alloc_points_with_offset(unsigned npts) {
		pts.resize(pts.size() + npts);
		return get_offset_for_last_points_added();
	}
	unsigned get_fill_start(unsigned npts, unsigned six) const;
	void fill_pts_from(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color, unsigned six);
	void render_range (unsigned six, unsigned eix, unsigned num_instances=1) const;
	void render_single(unsigned ix, unsigned num_instances=1) const {render_range(ix, ix+1, num_instances);}
	void render_all(unsigned num_instances=1) const {if (has_data()) {render_range(0, offsets.size()-1, num_instances);}}
	bool upload();
	void update_range(typename vert_type_t::non_color_class const *const p, unsigned npts, colorRGBA const &color, unsigned six, unsigned eix);
	void update_range(vector<typename vert_type_t::non_color_class> const &v, colorRGBA const &color, unsigned six, unsigned eix) {
		assert(!v.empty());
		update_range(&v.front(), v.size(), color, six, eix);
	}
	void begin_render() const;
	void end_render() const {post_render();}
	void clear_points() {clear_cont(pts);}
	void swap_points(vector<vert_type_t> &other_pts) {pts.swap(other_pts);}
	vector<vert_type_t> &get_pts_vector_for_adding() {return pts;}
	void clear(bool free_pts_mem=1);
	size_t get_gpu_mem() const {return ((vbo_valid() && has_data()) ? offsets.back()*sizeof(vert_type_t) : 0);}
};


typedef vbo_block_manager_t<vert_color> vbo_vc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_comp_color> vbo_vnc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_tc_color  > vbo_vntc_block_manager_t;
typedef vbo_block_manager_t<vert_norm_tc        > vbo_vnt_block_manager_t;
typedef vbo_block_manager_t<vert_norm_comp_tc   > vbo_vnct_block_manager_t;


class lt_atten_manager_t {

	shader_t &shader;
	int ulocs[5];
	float last_light_atten, last_refract_ix;

public:
	lt_atten_manager_t(shader_t &shader_) : shader(shader_), last_light_atten(-2.0), last_refract_ix(0.0) // set to invalid values to start
	{ulocs[0] = ulocs[1] = ulocs[2] = ulocs[3] = ulocs[4] = -1;}
	void enable();
	void next_object(float light_atten, float refract_ix);
	void next_cube(float light_atten, float refract_ix, cube_t const &cube);
	void next_sphere(float light_atten, float refract_ix, point const &pos, float radius);
	bool is_enabled() const {return (ulocs[0] >= 0);} // check that first uniform loc is assigned
};


class reflect_plane_selector {

	vector<cube_t> bcubes;
	int sel_cube;

public:
	reflect_plane_selector() : sel_cube(-1) {}
	bool empty() const {return bcubes.empty();}
	bool enabled() const {return (!empty() && sel_cube >= 0);}
	void add(cube_t const &c) {bcubes.push_back(c);}
	cube_t const &get_selected() const {assert(enabled()); assert(sel_cube < (int)bcubes.size()); return bcubes[sel_cube];}
	float get_refl_plane() const {cube_t const &c(get_selected()); return 0.5f*(c.d[2][0] + c.d[2][1]);}
	void select_best_reflection_plane();
};


struct fire_drawer_t {

	quad_batch_draw qbd;

	void add_fire(point const &pos, float radius, int frame_ix, float alpha=1.0);
	void draw(shader_t &s);
};

