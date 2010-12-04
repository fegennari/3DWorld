// 3D World - Quad/Triangle Batched Drawing Code
// by Frank Gennari
// 12/4/10

#include "quad_tri_drawer.h"


template<typename T> size_t quad_tri_drawer<T>::get_num_prim() const {

	if (out_type == GL_TRIANGLES) {
		assert((size()%3) == 0);
		return (size()/3);
	}
	else if (out_type == GL_QUADS) {
		assert((size()%4) == 0);
		return (size()/4);
	}
	return 0;
}


template<typename T> void quad_tri_drawer<T>::begin(int type) {

	assert(!in_strip && !dpos);

	switch (type) {
			case GL_TRIANGLE_STRIP:
				in_strip = 1;
			case GL_TRIANGLES:
				in_type  = GL_TRIANGLES;
				break;
			case GL_QUAD_STRIP:
				in_strip = 1;
			case GL_QUADS:
				in_type  = GL_QUADS;
				break;
			default:
				assert(0);
	}
}


template<typename T> void quad_tri_drawer<T>::end() {

	in_type  = GL_POINTS;
	in_strip = 0;
	dpos     = 0;
}


template<typename T> void quad_tri_drawer<T>::add_raw(T const &v) {

	assert(in_type != GL_POINTS); // begin() never called
	add_internal(v);

	if (in_type == GL_TRIANGLES && out_type == GL_QUADS) {
		// convert 1 triangle into 1 quad by duplicating the last point
		assert(!in_strip); // not yet supported
		if ((dpos&3) == 3) add_internal(v); // last vertex of a triangle, duplicate the last point
	}
	else if (in_type == GL_QUADS && out_type == GL_TRIANGLES) {
		// convert 1 quad into 2 triangles by duplicating 2 intermediate points and reordering
		assert(0); // not yet supported
	}
	else {
		assert(in_type == out_type);
	}
}


template<typename T> void quad_tri_drawer<T>::add_vert(T const &v) {

	if (in_strip) {
		unsigned const sz(size());

		if (in_type == GL_TRIANGLES && dpos >= 3) {
			add_internal(data[sz-2]); // duplicate the last 2 points for each triangle after the first one
			add_internal(data[sz-1]);
		}
		else if (in_type == GL_QUADS && dpos >= 4 && (dpos&3) == 0) {
			add_internal(data[sz-1]); // duplicate the last 2 points for each quad after the first one
			add_internal(data[sz-2]); // note that they have been swapped so we need to swap them back
		}
	}
	add_raw(v);

	if (in_strip && in_type == GL_QUADS) {
		if ((dpos&3) == 0) swap(data[size()-2], data[size()-1]); // swap the last two vertices
	}
}


template<typename T> void quad_tri_drawer<T>::draw(unsigned vbo) const { // vbo can be 0

	assert(!in_strip && !dpos);
	if (empty() || get_num_prim() == 0) return; // will do error checking
	data[0].set_state(vbo);
	glDrawArrays(out_type, 0, size());
}


template class quad_tri_drawer<vert_color>;
