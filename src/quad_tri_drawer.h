// 3D World - Quad/Triangle Batched Drawing Code
// by Frank Gennari
// 12/4/10

#ifndef _QUAD_TRI_DRAW_H_
#define _QUAD_TRI_DRAW_H_

#include "3DWorld.h"


template<typename T> class quad_tri_drawer {

	vector<T> data;
	int in_type, out_type;
	bool in_strip;
	unsigned dpos;

public:
	// Note: in_type starts as GL_POINTS, which is something known to be invalid here
	quad_tri_drawer(int out_type_) : in_type(GL_POINTS), out_type(out_type_), in_strip(0), dpos(0) {
		assert(out_type == GL_TRIANGLES || out_type == GL_QUADS);
	}
	void clear() {data.resize(0); end();}
	bool  empty() const {return data.empty();}
	size_t size() const {return data.size();}
	size_t get_num_prim() const;
	void begin(int type);
	void end();
	void add_raw(T const &v);
	void add_vert(T const &v);
	void draw(unsigned vbo) const;
	void draw_and_clear(unsigned vbo) {draw(vbo); clear();}

	void add_internal(T const &v) {
		data.push_back(v);
		++dpos;
	}
};

#endif // _QUAD_TRI_DRAW_H_
