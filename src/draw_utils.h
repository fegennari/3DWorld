// 3D World - Drawing Utility Classes
// by Frank Gennari
// 6/15/12

#ifndef _DRAW_UTILS_H_
#define _DRAW_UTILS_H_

#include "3DWorld.h"


template<typename cwt> class pt_line_drawer_t { // and triangles too!

	struct vnc : public cwt { // size = 28
		point v;
		vector3d n;

		vnc() {}
		vnc(point const &v_, vector3d const &n_, colorRGBA const &c_) : v(v_), n(n_) {set_c4(c_);}
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


class quad_batch_draw { // unused, but could possibly use for pine trees and plants
	vector<vert_norm_tc_color> verts;

public:
	void add_quad_vect(vector<vert_norm> const &points, colorRGBA const &color);
	void draw() const;
	void draw_and_clear() {draw(); verts.resize(0);}
	size_t size() const {return verts.size();}
	void reserve(size_t sz) {verts.reserve(sz);}
};


class vbo_quad_block_manager_t {

	typedef vert_norm_tc_color vert_type_t;
	vector<vert_type_t> pts;
	vector<unsigned> offsets;
	unsigned vbo;

public:
	vbo_quad_block_manager_t() : vbo(0) {clear();}
	~vbo_quad_block_manager_t() {clear_vbo();}
	bool empty() const {return pts.empty();}
	unsigned add_points(vector<vert_norm> const &p, colorRGBA const &color);
	void render_range(unsigned six, unsigned eix) const;
	void render_all() const {if (!empty()) {render_range(0, offsets.size()-1);}}
	void upload();
	void begin_render() const;
	void end_render() const;
	void clear_vbo();
	void clear();
};


#endif // _DRAW_UTILS_H_

