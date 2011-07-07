// 3D World - GLUT Polygon Tessellation code
// by Frank Gennari
// 5/29/06

#include "3DWorld.h"
#include "collision_detect.h"


#ifdef _WIN32
#define fgCALLBACK CALLBACK
#else
#define fgCALLBACK GLAPIENTRY
#endif


triangle cur_triangle;
deque<triangle> triangles;
deque<point> added_pts;
bool mode_valid(0), self_int(0), has_tess_error(0);
int mode(0);
unsigned vertex(0);



void fgCALLBACK tess_error(GLenum errno) {

	cout << "Polygon tesselation error " << errno << " has occurred: " << gluErrorString(errno) << endl;
	has_tess_error = 1;
}


void fgCALLBACK do_combine(GLdouble c[3], void *[4], GLfloat [4], void **out) { // self-intersecting polygon

	self_int = 1; // do nothing
	added_pts.push_back(point(c[0], c[1], c[2]));
	*out = (void *)(&added_pts.back());
}


void fgCALLBACK do_coord(void *coord) {

	assert(coord);
	assert(mode_valid);
	assert(vertex <= 3);
	cur_triangle.pts[vertex] = *((point *)coord);
	bool const emit_tri(vertex == 2);
	if (emit_tri) triangles.push_back(cur_triangle);

	switch (mode) {
	case GL_TRIANGLES:
		if (emit_tri) vertex = 0;
		break;
	case GL_TRIANGLE_STRIP:
		if (emit_tri) {
			for (unsigned i = 0; i < 2; ++i) cur_triangle.pts[i] = cur_triangle.pts[i+1];
		}
		break;
	case GL_TRIANGLE_FAN:
		if (emit_tri) cur_triangle.pts[1] = cur_triangle.pts[2];
		break;
	default:
		assert(0);
	}
	if (!emit_tri) ++vertex;
}


void fgCALLBACK do_begin(int type) {

	assert(!mode_valid);
	mode       = type;
	mode_valid = 1;
	vertex     = 0;
}


void fgCALLBACK do_end() {

	assert(mode_valid);
	assert(!triangles.empty());
	mode_valid = 0;

	switch (mode) {
	case GL_TRIANGLES:
		assert(vertex == 0);
		break;
	case GL_TRIANGLE_STRIP:
	case GL_TRIANGLE_FAN:
		assert(vertex == 2);
		break;
	default:
		assert(0);
	}
}


GLUtesselator *init_tess() {

	GLUtesselator *tobj(gluNewTess());
	assert(tobj != NULL);
	gluTessCallback(tobj, GLU_TESS_VERTEX,  (void (fgCALLBACK *)(void))do_coord);
	gluTessCallback(tobj, GLU_TESS_BEGIN,   (void (fgCALLBACK *)(void))do_begin);
	gluTessCallback(tobj, GLU_TESS_END,     (void (fgCALLBACK *)(void))do_end);
	gluTessCallback(tobj, GLU_TESS_ERROR,   (void (fgCALLBACK *)(void))tess_error);
	gluTessCallback(tobj, GLU_TESS_COMBINE, (void (fgCALLBACK *)(void))do_combine);
	return tobj;
}


void tessellate_polygon(vector<point> const &poly_pts) {

	assert(!has_tess_error);
	unsigned const size(poly_pts.size());
	assert(size >= 3);
	static GLUtesselator *tobj = NULL;
	if (tobj == NULL) tobj = init_tess();
	mode_valid = 0;
	self_int   = 0;
	vertex     = 0;
	added_pts.clear();
	vector3d const norm(get_poly_norm(&poly_pts.front()));
	gluTessNormal(tobj, norm.x, norm.y, norm.z);
	gluTessBeginPolygon(tobj, (void *)(&poly_pts.front()));
	gluTessBeginContour(tobj);
	double coord[3];

	for (unsigned i = 0; i < size; ++i) {
		for (unsigned j = 0; j < 3; ++j) coord[j] = poly_pts[i][j];
		gluTessVertex(tobj, coord, (void *)&poly_pts[i]);
	}
	gluTessEndContour(tobj);
	gluTessEndPolygon(tobj);
	assert(self_int == !added_pts.empty());

	if (has_tess_error) { // must delete and recreate tess object
		has_tess_error = 0;
		gluDeleteTess(tobj);
		tobj = NULL;
	}
	else {
		assert(!mode_valid);
	}
	if (self_int) cout << "* Warning: Self-intersecting polygon." << endl;
}


// unused, but could be useful
void split_polygon_to_tris(vector<triangle> &triangles_out, vector<point> const &poly_pts) {

	unsigned const npts(poly_pts.size());
	assert(npts >= 3);
	
	if (npts == 3) {
		triangles_out.push_back(triangle(poly_pts[0], poly_pts[1], poly_pts[2]));
	}
	else {
		tessellate_polygon(poly_pts);
		copy(triangles.begin(), triangles.end(), back_inserter(triangles_out));
		triangles.clear();
	}
}


bool is_poly_convex(vector<point> const &points) { // untested

	unsigned const npts(points.size());
	assert(npts >= 3);
	if (npts == 3) return 1;
	unsigned counts[2] = {0};
	vector3d const norm(get_poly_norm(&points.front()));

	for (unsigned i = 0; i < npts; ++i) {
		unsigned const ip((i+npts-1)%npts), in((i+1)%npts);
		++counts[dot_product(norm, cross_product(points[i]-points[ip], points[in]-points[i])) < 0.0];
	}
	return !(counts[0] && counts[1]);
}


bool split_polygon_to_cobjs(coll_obj const &cobj, vector<coll_obj> &split_polygons, vector<point> const &poly_pts, bool split_quads) {

	unsigned const npts(poly_pts.size());
	assert(npts >= 3);
	
	if (npts <= N_COLL_POLY_PTS && !(split_quads && npts > 3) && is_poly_convex(poly_pts)) { // triangle or convex quad
		if (!is_poly_valid(&poly_pts.front())) return 0; // invalid zero area polygon - skip
		split_polygons.push_back(cobj);
		split_polygons.back().npoints = npts;

		for (unsigned i = 0; i < npts; ++i) {
			split_polygons.back().points[i] = poly_pts[i];
		}
		return 0;
	}
	tessellate_polygon(poly_pts);
	
	for (unsigned i = 0; i < triangles.size(); ++i) {
		if (!is_poly_valid(triangles[i].pts)) continue; // invalid zero area triangle - skip
		split_polygons.push_back(cobj);
		split_polygons.back().npoints = 3; // triangles
		UNROLL_3X(split_polygons.back().points[i_] = triangles[i].pts[i_];)
	}
	triangles.clear();
	assert(!split_polygons.empty());
	return 1;
}





