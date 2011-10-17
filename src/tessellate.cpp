// 3D World - GLUT Polygon Tessellation code
// by Frank Gennari
// 5/29/06

#include "3DWorld.h"
#include "collision_detect.h"
#include "model3d.h"
#include "subdiv.h"


#ifdef _WIN32
#define fgCALLBACK CALLBACK
#else
#define fgCALLBACK GLAPIENTRY
#endif


triangle_vntc cur_triangle;
deque<triangle_vntc> triangles;
deque<vert_norm_tc> added_pts;
bool mode_valid(0), self_int(0), has_tess_error(0);
int mode(0);
unsigned vertex(0);



void fgCALLBACK tess_error(GLenum errno) {

	cout << "Polygon tesselation error " << errno << " has occurred: " << gluErrorString(errno) << endl;
	has_tess_error = 1;
}


void fgCALLBACK do_combine(GLdouble c[3], vert_norm_tc *vertex_data[4], GLfloat weight[4], vert_norm_tc **out) { // self-intersecting polygon

	self_int = 1;
	vector3d norm(zero_vector);
	float tc[2] = {0.0, 0.0};

	for (unsigned i = 0; i < 4; ++i) {
		if (!vertex_data[i]) { // the spec says this can't happen, but it does seem to occur
			assert(weight[i] == 0.0);
			continue;
		}
		norm  += vertex_data[i]->n    * weight[i];
		tc[0] += vertex_data[i]->t[0] * weight[i];
		tc[1] += vertex_data[i]->t[1] * weight[i];
	}
	added_pts.push_back(vert_norm_tc(point(c[0], c[1], c[2]), norm, tc));
	*out = &added_pts.back();
}


void fgCALLBACK do_coord(vert_norm_tc *coord) {

	assert(coord);
	assert(mode_valid);
	assert(vertex <= 3);
	cur_triangle.pts[vertex] = *coord;
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


void tessellate_polygon(polygon_t const &poly) {

	assert(!has_tess_error);
	unsigned const size(poly.size());
	assert(size >= 3);
	static GLUtesselator *tobj = NULL;
	if (tobj == NULL) tobj = init_tess();
	mode_valid = 0;
	self_int   = 0;
	vertex     = 0;
	added_pts.clear();
	vector3d const norm(poly.get_planar_normal());
	gluTessNormal(tobj, norm.x, norm.y, norm.z);
	gluTessBeginPolygon(tobj, (void *)(&poly.front()));
	gluTessBeginContour(tobj);
	double coord[3];

	for (unsigned i = 0; i < size; ++i) {
		UNROLL_3X(coord[i_] = poly[i].v[i_];)
		gluTessVertex(tobj, coord, (void *)&poly[i]);
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
	if (self_int) {
		static bool had_self_int_warning(0);

		if (!had_self_int_warning) {
			cout << "* Warning: Self-intersecting polygon." << endl;
			had_self_int_warning = 1;
		}
	}
}


bool split_polygon(polygon_t const &poly, vector<polygon_t> &ppts) {

	unsigned const npts(poly.size());
	assert(npts >= 3);
	
	if (npts <= 4 && (npts == 3 || poly.is_convex())) { // triangle or convex quad
		if (!poly.is_valid()) return 0; // invalid zero area polygon - skip
		ppts.push_back(poly);
		return 1;
	}
	tessellate_polygon(poly);

	// calculate polygon normal (assuming planar polygon)
	vector3d n(poly.get_planar_normal()), cp_sum(zero_vector);

	for (unsigned i = 0; i < npts; ++i) {
		cp_sum += cross_product(poly[i].v, poly[(i+1)%npts].v);
	}
	if (dot_product(n, cp_sum) < 0.0) n *= -1.0;
	static polygon_t new_poly;
	new_poly.resize(3);
	new_poly.color = poly.color;

	// triangles can be empty if they're all small fragments that get dropped
	for (unsigned i = 0; i < triangles.size(); ++i) {
		UNROLL_3X(new_poly[i_] = triangles[i].pts[i_];)
		if (!new_poly.is_valid()) continue; // invalid zero area triangle - skip
		if (dot_product(new_poly.get_planar_normal(), n) < 0.0) swap(new_poly[0], new_poly[2]); // invert draw order
		ppts.push_back(new_poly);
	}
	// triangles and split_polygons can be empty here if they're all small fragments that get dropped
	triangles.clear();
	return 1;
}


void split_polygon_to_cobjs(coll_obj const &cobj, vector<coll_obj> &split_polygons, vector<point> const &poly_pts, bool split_quads) {

	if (poly_pts.size() == 3 && is_poly_valid(&poly_pts.front())) { // optimization
		split_polygons.push_back(cobj);
		split_polygons.back().set_from_pts(&poly_pts.front(), 3);
		return;
	}
	static polygon_t poly;
	static vector<polygon_t> ppts;
	ppts.resize(0);
	poly.from_points(poly_pts);
	split_polygon(poly, ppts);

	for (vector<polygon_t>::const_iterator i = ppts.begin(); i != ppts.end(); ++i) {
		split_polygons.push_back(cobj);
		copy_polygon_to_cobj(*i, split_polygons.back()); // Note: poly color is ignored
	}
}


void cobj_triangle_visitor::proc_poly(point const *p, unsigned npts) {

	if (npts == 3) {
		proc_tri(p);
	}
	else if (npts == 4) {
		proc_quad(p);
	}
	else {
		assert(0);
	}
}


void cobj_triangle_visitor::proc_cobj(coll_obj const &c) {

	if (c.status != COLL_STATIC && (skip_dynamic || c.status != COLL_DYNAMIC)) return;

	switch (c.type) {
	case COLL_CUBE:
		for (unsigned dim = 0; dim < 3; ++dim) {
			for (unsigned dir = 0; dir < 2; ++dir) {
				if (c.cp.surfs & EFLAGS[dim][dir]) continue;
				unsigned const d0((dim+1)%3), d1((dim+2)%3);
				point pts[4], p;
				p[dim] = c.d[dim][dir];
				p[d0 ] = c.d[d0][0];
				p[d1 ] = c.d[d1][0]; pts[0] = p;
				p[d0 ] = c.d[d0][1]; pts[1] = p;
				p[d1 ] = c.d[d1][1]; pts[2] = p;
				p[d0 ] = c.d[d0][0]; pts[3] = p;
				proc_quad(pts);
			}
		}
		break;

	case COLL_SPHERE:
		{
			unsigned const ndiv(N_SPHERE_DIV);
			sd_sphere_d sd(c.points[0], c.radius, ndiv);
			sd.gen_points_norms();
			point **points = sd.get_points();
	
			for (unsigned s = 0; s < ndiv; ++s) {
				unsigned const sn((s+1)%ndiv);

				for (unsigned t = 0; t < ndiv; ++t) {
					point const pts[4] = {points[s][t], points[sn][t], points[sn][t+1], points[s][t+1]};
					proc_quad(pts);
				}
			}
		}
		break;

	case COLL_CYLINDER:
	case COLL_CYLINDER_ROT:
		{
			unsigned const ndiv(N_CYL_SIDES);
			assert(c.radius > 0.0 || c.radius2 > 0.0);
			bool const draw_ends(!(c.cp.surfs & 1));
			point const ce[2] = {c.points[0], c.points[1]};
			vector3d v12; // (ce[1] - ce[0]).get_norm()
			vector_point_norm const &vpn(gen_cylinder_data(ce, c.radius, c.radius2, ndiv, v12));

			for (unsigned S = 0; S < ndiv; ++S) { // ndiv can change
				unsigned const s[2] = {S, (S+ndiv-1)%ndiv};
				point pts[4];

				for (unsigned i = 0; i < 4; ++i) {
					unsigned const ss(s[i>>1]);
					pts[i] = vpn.p[(ss<<1) + (i==1||i==2)];
				}
				proc_quad(pts);
			}
			if (draw_ends) {
				float const r[2] = {c.radius, c.radius2};

				for (unsigned i = 0; i < 2; ++i) {
					if (r[i] == 0.0) continue;
					point pts[3];
					pts[0] = ce[i];

					for (unsigned s = 0; s < ndiv; ++s) { // ndiv can change
						pts[1] = vpn.p[(s<<1)+i];
						pts[1] = vpn.p[(((s+1)%ndiv)<<1)+i];
						proc_tri(pts);
					}
				}
			}
		}
		break;

	case COLL_POLYGON:
		if (c.thickness > MIN_POLY_THICK) { // thick polygon
			vector3d const norm(get_poly_norm(c.points));
			vector<point> pts[2];
			gen_poly_planes(c.points, c.npoints, norm, fabs(c.thickness), pts);
			for (unsigned d = 0; d < 2; ++d) proc_poly(&pts[d].front(), c.npoints); // top/bottom
	
			for (int i = 0; i < c.npoints; ++i) { // sides
				int const ii((i+1)%c.npoints);
				point const side_pts[4] = {pts[0][i], pts[0][ii], pts[1][ii], pts[1][i]};
				proc_quad(side_pts);
			}
		}
		else {
			proc_poly(c.points, c.npoints);
		}
		break;

	default:
		assert(0);
	}
}


void cobj_triangle_visitor::proc_cobjs(vector<coll_obj> const &cobjs) {

	RESET_TIME;
	tris_visited = 0;

	for (vector<coll_obj>::const_iterator i = cobjs.begin(); i != cobjs.end(); ++i) {
		proc_cobj(*i);
	}
	PRINT_TIME("Visit Cobj Triangles");
	cout << "triangles visited: " << tris_visited << endl;
}



