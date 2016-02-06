// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/24/02

#include "function_registry.h"
#include "sinf.h"


extern float orig_timestep, base_gravity;
extern int display_mode; // for debugging


bool line_intersect_torus(double ax, double ay, double az, double bx, double by,
						  double bz, double R, double r, double vlength, float &t);


// ************ BASIC VECTOR MATH, ETC. ************


float fix_angle(float angle) { // not sure if this is really necessary since sin/cos functions should be able to handle large angles

	if      (angle > TWO_PI) {angle -= TWO_PI;}
	else if (angle < 0.0)    {angle += TWO_PI;}
	if (angle == -0.0)       {angle = 0.0;} // stupid -0
	return angle;
}


void calc_reflection_angle(vector3d const &v_inc, vector3d &v_ref, vector3d const &norm) {

	float const cos_t1(-dot_product(norm, v_inc));
	v_ref = v_inc + norm*(2.0*cos_t1);
}


bool calc_refraction_angle(vector3d const &v_inc, vector3d &v_ref, vector3d const &norm, float n1, float n2) {

	assert(n2 != 0.0);
	float const cos_t1(-dot_product(norm, v_inc)), n_ratio(n1/n2);
	float const arg(1.0 - n_ratio*n_ratio*(1.0 - cos_t1*cos_t1));
	if (arg < 0.0) return 0; // total internal reflection - no refraction angle
	float const cos_t2(sqrt(arg));
	v_ref = v_inc*n_ratio + norm*(n_ratio*cos_t1 - fabs(cos_t2));
	return 1;
}


float get_fresnel_reflection(vector3d const &v_inc, vector3d const &norm, float n1, float n2) { // vectors must be normalized

	// http://en.wikipedia.org/wiki/Fresnel_equations
	float const cos_theta_i(dot_product(v_inc, norm)), sin_theta_i(cross_product(v_inc, norm).mag());
	float const val((n1/n2)*sin_theta_i), cos_theta_t(sqrt(1.0 - val*val));
	float const rs_sqrt((n1*cos_theta_i - n2*cos_theta_t)/(n1*cos_theta_i + n2*cos_theta_t));
	float const rp_sqrt((n1*cos_theta_t - n2*cos_theta_i)/(n1*cos_theta_t + n2*cos_theta_i));
	float const r(0.5*(rs_sqrt*rs_sqrt + rp_sqrt*rp_sqrt)); // average of rs and rp
	return r;
}


float get_reflected_weight(float fresnel_ref, float alpha) {
	return (alpha + (1.0 - alpha)*CLIP_TO_01(fresnel_ref));
}

float get_coll_energy(vector3d const &v1, vector3d const &v2, float mass) {

	if (v1 == v2) return 0.0;
	float const t(orig_timestep/DEF_TIMESTEP), vsq(fabs(v1.mag_sq() - v2.mag_sq()));
	return ((vsq < TOLERANCE) ? 0.0 : 0.5*mass*vsq*t*t); // 0.5*M*V^2 (should t be squared?)
}

point triangle_centroid(point const &p1, point const &p2, point const &p3) {return (p1 + p2 + p3)/3.0;}

float triangle_area(point const &p1, point const &p2, point const &p3) {

	float const a(p2p_dist(p1, p2)), b(p2p_dist(p2, p3)), c(p2p_dist(p3, p1));
	return 0.25*sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c));
}

float polygon_area(point const *const points, unsigned npoints) {

	assert(npoints == 3 || npoints == 4);
	float area(triangle_area(points[0], points[1], points[2]));
	if (npoints == 4) {area += triangle_area(points[2], points[3], points[0]);} // other triangle
	return area;
}

float get_closest_pt_on_line_t(point const &pos, point const &l1, point const &l2) { // use for line lights as well?
	vector3d const L(l2 - l1);
	float const L_mag(L.mag_sq());
	return ((L_mag < TOLERANCE) ? 0.0 : CLIP_TO_01(dot_product((pos - l1), L)/L_mag)); // clipped: line segment, not infinite line
}
point get_closest_pt_on_line(point const &pos, point const &l1, point const &l2) {
	return l1 + get_closest_pt_on_line_t(pos, l1, l2)*(l2 - l1);
}


// ************ SHAPE INTERSECTION ************


float const UV_TOLER(1.0E-6);

inline bool test_0_1(float v) {

	return (v >= (0.0 + UV_TOLER) && v <= (1.0 - UV_TOLER));
}


// from Graphics Gems V - actually could use quadrangle version but it seems to have problems
bool planar_contour_intersect(point const *points, unsigned npoints, point const &pos, vector3d const &norm) {

	if (npoints < 3) return 0;
	assert(npoints <= 4); // use get_max_dim()?
	unsigned const dx((fabs(norm.x) > max(fabs(norm.y), fabs(norm.z))) ? 1 : 0); // x is largest
	unsigned const dy((fabs(norm.z) > max(fabs(norm.x), fabs(norm.y))) ? 1 : 2); // z is largest
	point2d<double> const A(points[0][dx], points[0][dy]);
	point2d<double> const B(points[1][dx], points[1][dy]);
	point2d<double> const C(points[2][dx], points[2][dy]);
	point2d<double> const M(pos[dx], pos[dy]), AM(M, A), AB(B, A), AC(C, A);
	double const d(AM.cp_mag(AC)), a(AB.cp_mag(AC));

	if (fabs(a) > TOLERANCE) {
		double const u(d/a);
		if (test_0_1(u)) {
			double const v(AB.cp_mag(AM)/a);
			if (test_0_1(v) && test_0_1(u+v)) return 1;
		}
	}
	if (npoints == 4) { // quad: test other triangle
		point2d<double> const D(points[3][dx], points[3][dy]);
		point2d<double> const AD(D, A);
		double const a(AD.cp_mag(AC));

		if (fabs(a) > TOLERANCE) {
			double const u(d/a);
			if (test_0_1(u)) {
				double const v(AD.cp_mag(AM)/a);
				return (test_0_1(v) && test_0_1(u+v));
			}
		}
	}
	return 0;
}

bool point_in_polygon_2d(float sval, float tval, const point *points, int npts, int ds, int dt) {

	assert(ds >= 0 && ds <= 3 && dt >= 0 && dt <= 3 && ds != dt);
	assert(points != NULL && npts > 2);
	int nint(0);
	float s2(points[npts-1][ds]), t2(points[npts-1][dt]);

	for (int i = 0; i < npts; ++i) { // check if line (s1,t1), (s2,t2) intersects (sval,tval), (sout,tval)
		float const s1(points[i][ds]), t1(points[i][dt]);

		if ((t1 < tval) ^ (t2 < tval)) { // edge crosses t value
			//if ((s1 + (s2 - s1)*(tval - t1)/(t2 - t1)) < sval) { // s-intersection
			if (((s2 - s1)*(tval - t1) < (sval - s1)*(t2 - t1)) ^ (t2 > t1)) { // s-intersection
				if (t2 > t1) ++nint; else --nint;
			}
		}
		s2 = s1;
		t2 = t1;
	}
	return (nint != 0);
}

bool point_in_convex_planar_polygon(vector<point> const &pts, point const &normal, point const &pt) {

	if (pts.size() < 3) return 0; // not a polygon
	// Note: use the 2D projection method; could also use the normal and cross products
	int const proj_dim(get_min_dim(normal)), d1((proj_dim+1)%3), d2((proj_dim+2)%3);
	return point_in_polygon_2d(pt[d1], pt[d2], &pts.front(), pts.size(), d1, d2);
}


// z1 and z2 must be initialized (z1 > z2 to start)
bool get_poly_zminmax(point const *const pts, unsigned npts, vector3d const &norm, float dval, cube_t const &cube, float &z1, float &z2) {

	assert(pts && npts >= 3);
	unsigned num_inside[2] = {0, 0};

	for (unsigned i = 0; i < npts; ++i) { // test the points
		if (cube.contains_pt_xy(pts[i])) {
			z1 = min(z1, pts[i].z);
			z2 = max(z2, pts[i].z);
			++num_inside[0];
		}
	}
	if (num_inside[0] == npts) return 1;

	for (unsigned i = 0; i < npts; ++i) { // test the edges
		point edge[2] = {pts[i], pts[(i+1)%npts]};
		
		if (do_line_clip(edge[0], edge[1], cube.d)) {
			for (unsigned j = 0; j < 2; ++j) {
				z1 = min(z1, edge[j].z);
				z2 = max(z2, edge[j].z);
			}
			++num_inside[1];
		}
	}
	if (num_inside[1] == npts) return 1;
	point const center(cube.get_cube_center());

	if (point_in_polygon_2d(center.x, center.y, pts, npts, 0, 1)) {
		float const zv((fabs(norm.z) > 0.001) ? -(center.x*norm.x + center.y*norm.y + dval)/norm.z : pts[0].z);
		z1 = min(z1, zv);
		z2 = max(z2, zv);
		return 1;
	}
	return (num_inside[0] || num_inside[1]);
}


bool get_poly_zvals(vector<tquad_t> const &pts, float xv, float yv, float &z1, float &z2) {

	bool coll(0);

	for (unsigned p = 0; p < pts.size(); ++p) {
		if (point_in_polygon_2d(xv, yv, pts[p].pts, pts[p].npts, 0, 1)) coll = 1;
		vector3d const norm2(pts[p].get_norm());

		if (fabs(norm2.z) > TOLERANCE) {
			float const zv(-(xv*norm2.x + yv*norm2.y - dot_product(norm2, pts[p][0]))/norm2.z);
			z1 = min(z1, zv);
			z2 = max(z2, zv);
		}
	}
	return coll;
}


void gen_poly_planes(point const *const points, unsigned npoints, vector3d const &norm, float thick, point pts[2][4]) {

	for (unsigned i = 0; i < 2; ++i) { // back face cull?
		float const tv(0.5*(i ? -thick : thick));
		for (unsigned j = 0; j < npoints; ++j) {pts[i][j] = points[j] + norm*tv;}
	}
}


void thick_poly_to_sides(point const *const points, unsigned npoints, vector3d const &norm, float thick, vector<tquad_t> &sides) {

	assert(npoints >= 3);
	sides.resize(npoints + 2);

	for (unsigned i = 0; i < 2; ++i) { // same as gen_poly_planes()
		float const tv(0.5*(i ? -thick : thick));
		for (unsigned j = 0; j < npoints; ++j) {sides[i][j] = points[j] + norm*tv;}
		sides[i].npts = npoints;
	}
	for (unsigned i = 0; i < npoints; ++i) { // create the <npoints> sides
		unsigned const inext((i+1)%npoints);
		sides[i+2].npts = 4;
		sides[i+2][0] = sides[0][i];
		sides[i+2][1] = sides[1][i];
		sides[i+2][2] = sides[1][inext];
		sides[i+2][3] = sides[0][inext];
	}
	std::reverse(sides[1].pts, sides[1].pts+sides[1].npts); // reverse point order of bottom side
}


bool line_int_plane(point const &p1, point const &p2, point const &pp0, vector3d const &norm, point &p_int, float &t, bool ignore_t) {

	vector3d const v1(p2, p1);
	float const denom(dot_product(norm, v1)); // doesn't require norm to be normalized
	if (fabs(denom) < TOLERANCE) return 0;
	t = dot_product_ptv(norm, pp0, p1)/denom;
	if (!ignore_t && (t < 0.0 || t > 1.0)) return 0;
	p_int = p1 + v1*t;
	return 1;
}

bool line_int_plane_test_only(point const &p1, point const &p2, point const &pp0, vector3d const &norm) { // unused

	vector3d const v1(p2, p1);
	float const denom(dot_product(norm, v1)); // doesn't require norm to be normalized
	float const dp(dot_product_ptv(norm, pp0, p1));
	return (SIGN(dp) == SIGN(denom) && fabs(dp) <= fabs(denom));
}


// Note: this test considers the extruded polygon to be a collection of N planar polygons, not a volume,
// which means that lines completely inside the polygon volume (p1 and p2 contained) are not intersecting;
// if we need to treat this case as intersecting, we can use something like sphere_ext_poly_intersect()
bool thick_poly_intersect(vector3d const &v1, point const &p1, vector3d const &norm,
						  point const pts[2][4], bool test_side, unsigned npoints)
{ // test extruded (3D) polygon
	assert(npoints == 3 || npoints == 4);
	if (line_poly_intersect(v1, p1, pts[test_side], npoints, norm)) return 1;

	for (unsigned j = 0; j < npoints; ++j) { // now test the <npoints> sides
		unsigned const jnext((j+1)%npoints);
		point const side_pts[4] = {pts[0][j], pts[0][jnext], pts[1][jnext], pts[1][j]};
		if (line_poly_intersect(v1, p1, side_pts, 4, 1)) return 1;
	}
	return 0;
}


bool sphere_intersect_poly_sides(vector<tquad_t> const &pts, point const &center, float radius, float &dist, vector3d &norm, bool strict) {

	bool found(0);
	dist = FAR_DISTANCE;

	for (unsigned i = 0; i < pts.size(); ++i) { // test the <npoints> sides
		vector3d const side_norm(pts[i].get_norm());
		float tdist(radius - dot_product_ptv(side_norm, center, pts[i][0]));
		if (strict && tdist < 0.0) return 0; // outside of the shape
		
		if (!found || fabs(tdist) < fabs(dist)) {
			dist  = tdist;
			norm  = side_norm;
			found = 1;
		}
	}
	return found;
}


bool pt_line_seg_dist_less_than(point const &P, point const &L1, point const &L2, float dist) {

	if (dot_product(P-L1, P-L2) > 0.0) return 0; // pt not between s1 and s2
	return pt_line_dist_less_than(P, L1, L2, dist);
}

// Note: pt should lie in the plane of the polygon and be contained within the polygon
float min_dist_from_pt_to_polygon_edge(point const &pt, point const *const pts, unsigned npts) {

	assert(pts != nullptr && npts >= 3);
	float dmin(0.0);

	for (unsigned i = 0; i < npts; ++i) {
		float const dist(pt_line_dist(pt, pts[i], pts[(i+1 == npts) ? 0 : i+1]));
		dmin = ((i == 0) ? dist : min(dist, dmin));
	}
	return dmin;
}


bool sphere_poly_intersect(const point *points, unsigned npoints, point const &pos, vector3d const &norm, float rdist, float radius) {

	// test the points (point to point distance)
	for (unsigned i = 0; i < npoints; ++i) {
		if (dist_less_than(points[i], pos, radius)) return 1;
	}

	// test the edges (point to line distance)
	for (unsigned i = 0; i < npoints; ++i) {
		if (pt_line_seg_dist_less_than(pos, points[i], points[(i+1 == npoints) ? 0 : i+1], radius)) return 1;
	}

	// test for sphere center projected onto the polygon's plane
	return planar_contour_intersect(points, npoints, (pos - norm*rdist), norm);
}


// so many parameters for such a short function, but I want to make sure this isn't repeated and can be changed in only one place
bool sphere_ext_poly_int_base(point const &pt, vector3d const &norm, point const &pos, float radius,
							  float thickness, float &thick, float &rdist)
{
	thick = 0.5*thickness + radius;
	rdist = dot_product_ptv(norm, pos, pt); // not quite right, maybe should test each polygon face?
	return (fabs(rdist) <= thick);
}


bool sphere_ext_poly_intersect(point const * const points, unsigned npoints, vector3d const &norm,
							   point const &pos, float radius, float thickness, float t_adj)
{
	float thick, rdist;
	if (!sphere_ext_poly_int_base(points[0], norm, pos, radius, thickness, thick, rdist)) return 0;
	if (thickness <= MIN_POLY_THICK) return (sphere_poly_intersect(points, npoints, pos, norm, rdist, max(0.0f, (thick - t_adj))));
	vector<tquad_t> pts;
	thick_poly_to_sides(points, npoints, norm, thickness, pts); // slow

	for (unsigned i = 0; i < pts.size(); ++i) { // adapted from sphere_intersect_poly_sides()
		if ((radius - dot_product_ptv(pts[i].get_norm(), pos, pts[i][0])) < 0.0) return 0;
	}
	return 1;
}


// similar to sphere_test_comp(), but needs some intermediate values and assumes the line extends infinitely
// assumes infinite line - v12 should be normalized
bool line_intersect_sphere(point const &p1, vector3d const &v12, point const &sc, float radius, float &rad, float &dist, float &t) {

	vector3d const vsp(sc, p1);
	dist = vsp.mag();
	float const dotp(dot_product(v12, vsp));
	if (dotp < 0.0 && dist > radius) return 0; // pointing away from sphere and not inside it
	t    = dotp/dist; // t must be <= 1.0, assuming v12 is normalized
	rad  = sin(safe_acosf(t))*dist; // rad in (-dist, dist)
	return (radius > rad);
}


// p2 = starting point of the line, p1 = sphere center, v1 = line direction (inverted?), r2sq = square of the sphere radius
bool sphere_test_comp(point const &p2, point const &p1, vector3d const &v1, float r2sq, float &t) { // line segment/sphere intersection

	vector3d const v2(p2, p1);
	if (v2.mag_sq() <= r2sq) {t = 0.0; return 1;} // starting point is inside of the sphere
	float const dotp(dot_product(v1, v2));
	if (dotp < 0.0)        return 0; // line is pointing away from sphere
	float const denom(v1.mag_sq());
	if (denom < TOLERANCE) return 0; // point/sphere intersect
	if (dotp  > denom)     return (p2p_dist_sq(v1, v2) <= r2sq); // test to see if the end point is inside of the sphere
	t = dotp/denom;
	return ((v2 - v1*t).mag_sq() <= r2sq);
}


// v1 is the line dir, p1 is the line start point, {center, radius} define the sphere, and lsint is the intersect point
// similar to sphere_test_comp(), but assumes infinite line, so t test is invalid, and also needs the intersection point
bool line_sphere_int(vector3d const &v1, point const &p1, point const &center, float radius, point &lsint, bool test_neg_t) {

	lsint = center;
	vector3d v2(p1, center); // dir from target sphere to object
	float const t(-dot_product(v1, v2)); // v1 should be normalized
	if (test_neg_t && t < 0.0) return 0;
	v2 += v1*t;
	float const dsq(v2.mag_sq()); // length of perpendicular to sphere center
	if (dsq >= radius*radius)  return 0;
	lsint += v2 - v1*sqrt(radius*radius - dsq); // distance along line to outside surface of sphere
	return 1;
}


bool sphere_vert_cylin_intersect(point &center, float radius, cylinder_3dw const &c) {

	float const rsum(radius + max(c.r1, c.r2));
	assert(rsum > radius); // cylin not degenerate
	if (center.z < min(c.p1.z, c.p2.z) || center.z > max(c.p1.z, c.p2.z)) return 0; // not enough z overlap (could bias by up to radius)
	if (!dist_xy_less_than(c.p1, center, rsum)) return 0;
	vector3d const normal(vector3d(center.x-c.p1.x, center.y-c.p1.y, 0.0).get_norm()); // Note: assumes (near) vertical cylinder, not correct for all trees
	center = point(c.p1.x, c.p1.y, center.z) + normal*(1.001*rsum); // move center out along cylin normal so it doesn't intersect
	return 1;
}


void get_sphere_border_pts(point *qp, point const &pos, point const &viewed_from, float radius, unsigned num_pts) {

	assert(qp && num_pts <= 5);
	vector3d vortho;

	if (num_pts > 2) {
		cross_product((viewed_from - pos), plus_z, vortho); // xy direction
		vortho.normalize();
	}
	for (unsigned i = 0; i < num_pts; ++i) { // center, z+, xy-, xy+, z-
		qp[i] = pos;
		
		if      (i == 1 || i == 4) {qp[i].z += ((i == 1) ? radius : -radius);}
		else if (i == 2 || i == 3) {qp[i] += vortho*((i == 2) ? radius : -radius);}
	}
}


void get_sphere_points(point const &pos, float radius, point *pts, unsigned npts, vector3d const &dir) {

	vector3d const v1(cross_product(plus_z, dir).get_norm());
	vector3d const v2(cross_product(v1,     dir).get_norm());

	for (unsigned i = 0; i < npts; ++i) {
		float const theta(TWO_PI*i/npts);
		pts[i] = pos + v1*(sinf(theta)*radius) + v2*(cosf(theta)*radius);
	}
}


// dir must be normalized
void dir_to_sphere_s_t(vector3d const &dir, vector3d const &sdir, double &s, double &t) {

	double const stheta(atan2(sdir.y, sdir.x)), sphi(safe_acosf(-sdir.z));
	double const angle(-stheta + PI_TWO), sinv(sin(angle)), cosv(cos(angle));
	double const xval (-dir.x*cosv + dir.y*sinv), yval0(-dir.x*sinv - dir.y*cosv);
	double const angle2(-sphi + PI_TWO), sinv2(sin(angle2)), cosv2(cos(angle2));
	double const yval(yval0*cosv2 + dir.z*sinv2), zval(yval0*sinv2 - dir.z*cosv2);
	double const theta(atan2(yval, xval)), phi(safe_acosf(zval));
	s = (theta + PI_TWO)/TWO_PI + 1.0;
	t = phi/PI + 1.0;
}


// assumes infinite line - untransformed sdir = (1.0, 0.0, 0.0) ?
bool line_sphere_intersect_s_t(point const &p1, point const &p2, point const &sc, float radius,
							   vector3d const &sdir, double &s, double &t)
{
	point p_int;
	vector3d const v1(vector3d(p2, p1).get_norm());
	if (!line_sphere_int(v1, p1, sc, radius, p_int, 0)) return 0;
	vector3d const dir(vector3d(p_int, sc).get_norm());
	dir_to_sphere_s_t(dir, sdir, s, t);
	return 1;
}


// shortest distance between (p1b - p1a) and (p2b - p2a)
float line_line_dist(point const &p1a, point const &p1b, point const &p2a, point const &p2b) {

	vector3d const a(p1b, p1a), b(p2b, p2a), cp(cross_product(a, b));
	float const cp_mag(cp.mag());

	if (fabs(cp_mag) < TOLERANCE) { // lines are parallel
		vector3d const w(p2a - p1a), v_para(a*(dot_product(a, w)/a.mag_sq())), v_perp(w - v_para);
		return v_perp.mag();
	}
	return fabs(dot_product_ptv(cp, p2a, p1a))/cp_mag;
}


// similar to get_closest_pt_on_line_t(), but without the [0,1] clamp
float get_cylinder_params(point const &cp1, point const &cp2, point const &pos, vector3d &v1, vector3d &v2) {

	v1 = cp1 - cp2; // vector along length of cylinder
	v2 = cp1 - pos; // cylinder->object vector
	float const c_len(v1.mag_sq());
	assert(c_len > 0.0/*TOLERANCE*/);
	return dot_product(v1, v2)/c_len; // position of pos along cylinder center line often in [0,1]
}


// radius == 0.0 at cp1
int line_intersect_trunc_cone(point const &p1, point const &p2, point const &cp1, point const &cp2,
							  float r1, float r2, bool check_ends, float &t, bool swap_ends)
{
	// P = p1, V = cp1
	// M = A*A'*g*g*I, c2 = D'*M*D, c1 = D'*M*d, c0 = d'*M*d
	// c2*t^2 + 2*c1*t + c0 = 0
	// t = (-c1 +/- sqrt(c1*c1 - c2*c0))/c0
	// X = P + t*D
	assert(r2 > 0.0 && r2 > r1);
	point V(cp1);
	vector3d dir(cp2, cp1);
	if (r1 > 0.0) V -= dir*(r1/(r2 - r1)); // extend to the cone's apex
	vector3d A(cp2, V), D(p2 - p1), d(p1, V);
	float const g(cosf(atan2f(r2, A.mag())));
	A.normalize();
	double M[3][3];

	for (unsigned i = 0; i < 3; ++i) {
		UNROLL_3X(M[i][i_] = A[i]*A[i_];)
		M[i][i] -= g*g;
	}
	vector3d Md, MD;
	matrix_mult(d, Md, M);
	matrix_mult(D, MD, M);
	float c0(0.0), c1(0.0), c2(0.0);

	for (unsigned i = 0; i < 3; ++i) {
		c0 += D[i]*MD[i];
		c1 += D[i]*Md[i];
		c2 += d[i]*Md[i];
	}
	float num(c1*c1 - c2*c0);
	int t_set(0);

	if (num >= 0.0) {
		float const len(dir.mag());
		num = sqrt(num);

		for (unsigned i = 0; i < 2; ++i) {
			float const ti((-c1 + (1-2*i)*num)/c0);

			if (ti >= 0.0 && ti <= 1.0 && (!t_set || ti < t)) { // between p1 and p2
				float const dp(dot_product_ptv(A, (p1 + D*ti), cp1));

				if (dp >= 0.0 && dp <= len) { // test off beginning and end of cone (between cp1 and cp2)
					t     = ti;
					t_set = 1;
				}
			}
		}
	}
	if (check_ends) { // test end circle(s)
		float const r[2]  = {r1,  r2};
		point const cp[2] = {cp1, cp2};

		for (unsigned i = 0; i < 2; ++i) {
			float ti;

			if (r[i] > 0.0 && circle_test_comp(p1, cp[i], D, A, r[i]*r[i], ti)) {
				if (ti >= 0.0 && ti <= 1.0 && (!t_set || ti < t)) {
					t     = ti;
					t_set = (i ^ unsigned(swap_ends)) + 2;
				}
			}
		}
	}
	if (!t_set && check_ends && point_in_cylinder(cp1, cp2, p1, r1, r2)) { // what about p2 inside cylinder?
		t     = 0.0; // intersects at the starting pos
		t_set = 1; // sides or end?
	}
	return t_set;
}


// This one supports different radii at the ends and tests against a solid cylinder
bool line_intersect_cylinder(point const &p1, point const &p2, cylinder_3dw const &c, bool check_ends) {

	if (line_line_dist(p1, p2, c.p1, c.p2) > max(c.r1, c.r2)) return 0;
	int npts(0);
	point pts[4];
	vector3d const v1(p1, p2); // backwards
	vector3d v2, norm;
	cylinder_quad_projection(pts, c, v1, v2, npts);
	assert(npts > 2 && npts <= 4);
	get_normal(c.p1, c.p2, pts[c.r1 == 0.0], norm, 1);
	float const denom(dot_product(norm, v1));

	if (fabs(denom) > TOLERANCE) {
		float const t(dot_product_ptv(norm, c.p1, p2)/denom); // check if point is inside cylinder - close to correct

		if ((t >= 0.0 && t <= 1.0) || point_in_cylinder(c.p1, c.p2, p1, c.r1, c.r2) || point_in_cylinder(c.p1, c.p2, p2, c.r1, c.r2)) {
			vector3d const norm2(get_poly_norm(pts));
			if (planar_contour_intersect(pts, npts, (p2 + v1*t), norm2)) return 1;
		}
	}
	if (check_ends) {
		float t; // unused
		point const pt[2] = {p1, p2};
		point const cp[2] = {c.p1, c.p2};
		float const cr[2] = {c.r1, c.r2};

		for (unsigned d = 0; d < 2; ++d) {
			if (cr[d] > 0.0 && circle_test_comp(p2, cp[d], v1, v2, cr[d]*cr[d], t)) return 1; // end pt
			if (point_in_cylinder(c.p1, c.p2, pt[d], c.r1, c.r2)) return 1;
		}
	}
	return 0;
}


// This one calculates the intersection point and tests against a hollow cylinder
int line_int_thick_cylinder(point const &p1, point const &p2, point const &cp1, point const &cp2,
							float ri1, float ri2, float ro1, float ro2, bool check_ends, float &t)
{
	if (line_line_dist(p1, p2, cp1, cp2) > max(ro1, ro2)) return 0;

	if (ri1 == 0.0 && ri2 == 0.0 && (ro1 != ro2)) {
		if (ro1 < ro2) {
			return line_intersect_trunc_cone(p1, p2, cp1, cp2, ro1, ro2, check_ends, t, 0);
		}
		else {
			return line_intersect_trunc_cone(p1, p2, cp2, cp1, ro2, ro1, check_ends, t, 1); // reverse the ends
		}
	}
	assert(ri1 == ri2 && ro1 == ro2);
	assert(ro1 > 0.0 && ro2 > 0.0 && ri1 >= 0.0 && ri2 >= 0.0);
	assert(ro1 >= ri1 && ro2 >= ri2);
	//if (point_in_cylinder(cp1, cp2, p1, ro1, ro2) && !point_in_cylinder(cp1, cp2, p1, ri1, ri2)) return 3; // t = ???
	point v1(p1, cp1), v2(p2, cp1), c2(cp2, cp1); // translate so that cp1 is at (0,0,0), swap if r1 > r2?
	float const len(c2.mag());
	vector3d const cv(c2/len); // normalize
	rotate_vector3d_by_vr(cv, plus_z, v1); // rotate both points and cylinder so that the cylinder is oriented in +z
	rotate_vector3d_by_vr(cv, plus_z, v2);
	float const dz(v2.z - v1.z);
	float ta((0.0 - v1.z)/dz), tb((len - v1.z)/dz);
	bool const swapped(tb < ta);
	if (swapped) {swap(ta, tb);}
	if (ta > 1.0 || tb < 0.0) return 0; // doesn't cross between the cylinder ends
	float const dx(v2.x - v1.x), dy(v2.y - v1.y), dr2(dx*dx + dy*dy);
	
	if (ta >= 0.0) { // test intersection with the closest cylinder end
		if (check_ends) {
			float const xval(v1.x + ta*dx), yval(v1.y + ta*dy), dist_sq(xval*xval + yval*yval);
			if (dist_sq <= (swapped ? ro2*ro2 : ro1*ro1) && dist_sq >= (swapped ? ri2*ri2 : ri1*ri1)) {
				t = ta;
				return (swapped ? 3 : 2); // intersects the end circle
			}
		}
	}
	else {
		ta = 0.0;
	}
	if (dr2 < TOLERANCE) return 0; // line is parallel to cylinder
	t  = 2.0; // set > 1.0
	tb = min(1.0f, tb);
	float const D(v1.x*v2.y - v2.x*v1.y);
	unsigned const niter(1 + (ri1 != ro1 || ri2 != ro2));

	for (unsigned r = 0; r < niter; ++r) {
		float const disc((r ? ri1*ri1 : ro1*ro1)*dr2 - D*D); // can we generalize this if r1 != r2?
		if (disc < 0.0) continue; // no circle intersection (t should be 2.0)
		float const val(float(fabs(dy)*sqrt(disc)));

		for (unsigned i = 0; i < 2; ++i) { // 2 roots of quadratic equation
			float const yy((-D*dx + (i ? val : -val))/dr2), tt((yy - v1.y)/dy);
			if (tt >= ta && tt <= tb && tt < t) t = tt;
		}
	}
	return (t <= 1.0);
}


// spheres and vertical cylinders are rotationally invariant about z, so we can rotate a non-axis aligned cylinder about z
// to make it axis aligned, and compute a tighter bounding cube to use with a separating axis test
bool cylin_proj_circle_z_SAT_test(point const &cc, float cr, point const &cp1, point const &cp2, float r1, float r2) {
	point pts[2] = {cp1, cp2};
	vector3d const dir(cp2 - cp1);

	for (unsigned d = 0; d < 2; ++d) {
		pts[d] -= cc; // translate to circle center
		rotate_vector3d_by_vr(dir, plus_x, pts[d]); // rotate around vert cylinder
		pts[d] += cc; // translate back
	}
	cube_t bcube;
	cylinder_3dw(pts[0], pts[1], r1, r2).calc_bcube(bcube);
	return circle_rect_intersect(cc, cr, bcube, 2); // in z
}


// inaccurate if fabs(r2 - r1) > p2p_dist(cp1, cp2) == radius difference > cylinder length
bool sphere_int_cylinder_pretest(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2,
								 bool check_ends, vector3d &v1, vector3d &v2, float &t, float &rad)
{
	if (!cylin_proj_circle_z_SAT_test(sc, sr, cp1, cp2, r1, r2)) return 0;
	t   = get_cylinder_params(cp1, cp2, sc, v1, v2); // v1 = cylinder vector, v2 = cylinder_p1-sphere vector
	float const t_clamped(CLIP_TO_01(t));
	rad = (r1 + t_clamped*(r2 - r1)); // radius of cylinder at closest point to sphere

	if (cp1.z == cp2.z) {
		float const closest_z(cp1.z + t_clamped*(cp2.z - cp1.z)), sphere_dist(fabs(closest_z - sc.z));
		if (sphere_dist < sr) {rad += sqrt(sr*sr - sphere_dist*sphere_dist);}
	}
	else {
		rad += sr; // add sphere radius (FIXME: inaccurate)
	}
	if (check_ends || (t >= 0.0 && t <= 1.0)) {
		v2 -= v1*t; // vector from point to collision point along cylinder axis
		if (v2.mag_sq() <= rad*rad) return 1;
	}
	if (r1 == r2) return 0;
	vector3d cpb;
	orthogonalize_dir(v2, v1, cpb, 0);
	point l1(cp1), l2(cp2);
	float const mag(cpb.mag());
	
	if (mag > TOLERANCE) {
		if (r1 != 0.0) {l1 -= cpb*(r1/mag);}
		if (r2 != 0.0) {l2 -= cpb*(r2/mag);}
	}
	return line_sphere_intersect(l1, l2, sc, sr); // r1 != r2 - do a line/sphere intersection
}


bool sphere_intersect_cylinder_ipt(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2,
								   bool check_ends, point &p_int, vector3d &norm, bool calc_int)
{
	float t, rad;
	vector3d v1, v2; // v1: cp1-cp2, v2: cp1-sc
	if (!sphere_int_cylinder_pretest(sc, sr, cp1, cp2, r1, r2, check_ends, v1, v2, t, rad)) return 0;
	int const tok(t >= 0.0 && t <= 1.0);
	if (!calc_int && tok) return 1;
	unsigned npos(0);
	float dmin(0.0);
	point cpos[3];
	vector3d norms[3];
	float const len(v1.mag()), rdist(v2.mag());
	float const toler(0.0001);
	assert(len > 0.0);

	if (tok && rdist < rad) { // collision with side
		float const val(rad - rdist + toler);

		if (rdist < min(TOLERANCE, toler*rad)) { // center along cylinder centerline (rarely occurs)
			norm  = all_zeros;
			norm[get_min_dim(v1)] = 1.0; // move out of the way in an orthogonal direction
			p_int = sc + norm*val;
			return 1;
		}
		cpos[npos]  = sc;
		norms[npos] = v2;
		//if (check_ends || !point_in_cylinder(cp1, cp2, sc, r1, r2)) norms[npos].negate();
		norms[npos].negate();
		cpos[npos] += norms[npos]*(val/rdist); // r1 != r2 and !check_ends doesn't work well
		++npos;
	}
	if (check_ends) {
		for (unsigned d = 0; d < 2; ++d) {
			float const t_clamped(CLIP_TO_01(t)), tv(d ? (1.0f - t) : t), tv_clamped(d ? (1.0f - t_clamped) : t_clamped);

			if (((d ? r2 : r1) > 0.0) && (fabs(tv_clamped)*len < min(sr, rdist))) { // collision with p1/p2
				if (!calc_int) return 1;
				cpos[npos]  = sc;
				norms[npos] = v1;
				if (d) {norms[npos].negate();}
				
				if (len > TOLERANCE) {
					float const adj(tv + (sr + toler)/len);
					if (adj < 0.0) continue; // not a real intersection (due to fp error, etc.)
					cpos[npos] += norms[npos]*adj;
				}
				++npos;
			}
		}
	}
	if (npos == 0) return 0;

	for (unsigned p = 0; p < npos; ++p) {
		float const pd(p2p_dist(sc, cpos[p]));

		if (p == 0 || pd < dmin) {
			dmin  = pd;
			p_int = cpos[p];
			norm  = norms[p];
		}
	}
	norm.normalize();
	assert(!is_nan(p_int));
	assert(!is_nan(norm));
	return 1;
}


// line/sphere-scaled sphere/cylinder ???

// torus is oriented in the z direction (for now)
bool line_torus_intersect(point const &p1, point const &p2, point const &tc, float ri, float ro, float &t) {

	assert(ri >= 0.0 && ro > 0.0 && ri <= ro);
	float const rsum(ro + ri), dist[3] = {rsum, rsum, ri};

	for (unsigned i = 0; i < 3; ++i) {
		if (p1[i] > (tc[i] + dist[i]) && p2[i] > (tc[i] + dist[i])) return 0;
		if (p1[i] < (tc[i] - dist[i]) && p2[i] < (tc[i] - dist[i])) return 0;
	}
	//if (!sphere_test_comp(p1, tc, v1, rsum*rsum)) return 0; // clip (optional), could test bounding cube or cylinder as well
	vector3d v1(p2, p1); // line direction
	double const vmag(v1.mag());
	vector3d l1(p1 - tc); // translate to torus center
	
	if (p1 == tc) { // l1 is zero vector, this case doesn't work correctly
		l1 = p2 - tc;
		v1.negate();
	}
	return line_intersect_torus(l1.x, l1.y, l1.z, v1.x/vmag, v1.y/vmag, v1.z/vmag, ro, ri, vmag, t);
}


// torus is oriented in the z direction (for now)
bool sphere_torus_intersect(point const &sc, float sr, point const &tc, float ri, float ro, point &p_int, vector3d &norm, bool calc_int) {

	assert(sr >= 0.0 && ri >= 0.0 && ro > 0.0 && ri <= ro); // sphere radius can be 0
	float const r2s_sq((sr + ro + ri)*(sr + ro + ri));
	vector3d const t2s(sc, tc);
	if (t2s.mag_sq() > r2s_sq) return 0; // too far away
	float const dxy_sq(t2s.x*t2s.x + t2s.y*t2s.y);
	if (dxy_sq < TOLERANCE)    return 0; // TOLERANCE in here to avoid later divide-by-zero
	if (dxy_sq > r2s_sq)       return 0; // too far away in the x-y plane
	float const r1s_sq((-sr + ro - ri)*(-sr + ro - ri));
	if (dxy_sq < r1s_sq)       return 0; // too close to the center in the x-y plane
	float const dxy(sqrt(dxy_sq)), drxy(fabs(dxy - ro));
	float const rcs_sq((sr + ri)*(sr + ri)), dist_sq(drxy*drxy + t2s.z*t2s.z);
	if (dist_sq > rcs_sq)      return 0; // no intersection
	if (!calc_int)             return 1; // intersection, no p_int/norm required
	p_int    = tc;
	p_int.x += ro*t2s.x/dxy;
	p_int.y += ro*t2s.y/dxy; // intersection point along rcent in the z=tc.z plane
	norm     = sc - p_int;
	norm.normalize();
	p_int += norm*(sr + ri);
	return 1;
}


#define DMIN_CHECK(i) {if      (pos[i] < cube.d[i][0]) {float const dist(pos[i] - cube.d[i][0]); dmin += dist*dist;} \
					   else if (pos[i] > cube.d[i][1]) {float const dist(pos[i] - cube.d[i][1]); dmin += dist*dist;} \
					   if (dmin > r2) return 0;}

bool circle_rect_intersect(point const &pos, float radius, cube_t const &cube, int dim) {

	float dmin(0.0);
	float const r2(radius*radius);
	UNROLL_3X(if (dim != i_) {DMIN_CHECK(i_);})
	return 1;
}

bool sphere_cube_intersect(point const &pos, float radius, cube_t const &cube) {

	float dmin(0.0);
	float const r2(radius*radius);
	UNROLL_3X(DMIN_CHECK(i_));
	return 1;
}


void cylinder_3dw::calc_bcube(cube_t &bcube) const {

	vector3d const norm(get_norm_dir_vect());
			
	for (unsigned i = 0; i < 3; ++i) {
		float const ni(sqrt(1.0 - norm[i]*norm[i]));
		bcube.d[i][0] = min((p1[i] - ni*r1), (p2[i] - ni*r2));
		bcube.d[i][1] = max((p1[i] + ni*r1), (p2[i] + ni*r2));
	}
}


// Note: assumes a planar convex polygon, and assumes the cylinder is closed
bool approx_poly_cylin_int(point const *const pts, unsigned npts, cylinder_3dw const &cylin) {

	cube_t cbcube;
	cylin.calc_bcube(cbcube);
	if (!cube_t(pts, npts).intersects(cbcube)) return 0; // optimization
	float t;
	vector3d const poly_norm(get_poly_norm(pts)); // could factor out and pass in
	if (line_poly_intersect(cylin.p1, cylin.p2, pts, npts, poly_norm, t)) return 1;

	for (unsigned i = 0; i < npts; ++i) {
		if (line_intersect_cylinder(pts[i], pts[(i+1)%npts], cylin, 1)) return 1; // check_ends=1
	}
	return 0; // inexact, but close
}


// p_int is the location of the sphere at intersection, norm is the normal of the intersected surface
bool sphere_cube_intersect(point const &pos, float radius, cube_t const &cube, point const &p_last, point &p_int,
						   vector3d &norm, unsigned &cdir, bool check_int, bool skip_z)
{
	if (check_int && !sphere_cube_intersect(pos, radius, cube)) return 0;
	float min_dist(0.0);
	bool found(0);

	// first iteration:  find closest side where object has crossed a face of the cube, and if not found
	// second iteration: find closest side with no constraints
	// Note: still not exactly correct, but close, and OK in most cases
	for (unsigned iter = (pos == p_last); iter < 2 && !found; ++iter) {
		for (unsigned i = 0; i < unsigned(2 + !skip_z); ++i) {
			for (unsigned j = 0; j < 2; ++j) {
				//if (iter == 0 && pos[i] != p_last[i] && ((pos[i] - p_last[i]) < 0) ^ j) continue; // ignore back-facing sides (is this correct?)
				float const delta(j ? 1.0 : -1.0), side_pos(cube.d[i][j] + delta*radius);
				if (iter == 0 && !((p_last[i] < side_pos) ^ j) && ((pos[i] >= side_pos) ^ j)) continue;
				float const dist(fabs(pos[i] - side_pos));

				if (!found || dist < min_dist) {
					min_dist = dist;
					p_int    = pos;
					p_int[i] = side_pos;
					norm     = zero_vector;
					norm[i]  = delta;
					cdir     = (i << 1) + j; // x0 x1 y0 y1 z0 z1
					found    = 1;
				}
			}
		}
	}
	return found;
}


#define TEST_CLIP_T(reg, va, vb, vd, vc) \
	if (region3 & (reg)) { \
		assert((vd) != 0.0); /* this assertion may be removed as an optimization */ \
		float const t(((va) - (vb))/(vd)); \
		if ((vc) > 0.0) {if (t > tmin) tmin = t;} else {if (t < tmax) tmax = t;} \
		if (tmin >= tmax) return 0; \
	}

// performance critical
bool get_line_clip(point const &v1, point const &v2, float const d[3][2], float &tmin, float &tmax) {

	int const region1(get_region(v1, d)), region2(get_region(v2, d));
	if (region1 & region2) return 0; // line outside
	int const region3(region1 | region2);
	tmax = 1.0;
	tmin = 0.0;
	if (region3 == 0) return 1; // both points inside => entire line inside
	vector3d const dv(v2, v1);
	TEST_CLIP_T(0x01, d[0][0], v1.x, dv.x,  dv.x); // -x plane
	TEST_CLIP_T(0x02, d[0][1], v1.x, dv.x, -dv.x); // +x plane
	TEST_CLIP_T(0x04, d[1][0], v1.y, dv.y,  dv.y); // -y plane
	TEST_CLIP_T(0x08, d[1][1], v1.y, dv.y, -dv.y); // +y plane
	TEST_CLIP_T(0x10, d[2][0], v1.z, dv.z,  dv.z); // -z plane
	TEST_CLIP_T(0x20, d[2][1], v1.z, dv.z, -dv.z); // +z plane
	return 1;
}


// performance critical: return 1 if line intersects the cube
bool do_line_clip(point &v1, point &v2, float const d[3][2]) {

	int const region1(get_region(v1, d)), region2(get_region(v2, d));
	if (region1 & region2) return 0; // line outside
	int const region3(region1 | region2);
	if (region3 == 0) return 1; // both points inside => entire line inside
	float tmin(0.0), tmax(1.0);
	vector3d const dv(v2, v1);
	TEST_CLIP_T(0x01, d[0][0], v1.x, dv.x,  dv.x); // -x plane
	TEST_CLIP_T(0x02, d[0][1], v1.x, dv.x, -dv.x); // +x plane
	TEST_CLIP_T(0x04, d[1][0], v1.y, dv.y,  dv.y); // -y plane
	TEST_CLIP_T(0x08, d[1][1], v1.y, dv.y, -dv.y); // +y plane
	TEST_CLIP_T(0x10, d[2][0], v1.z, dv.z,  dv.z); // -z plane
	TEST_CLIP_T(0x20, d[2][1], v1.z, dv.z, -dv.z); // +z plane
	if (tmax > TOLERANCE)         {v2  = v1 + dv*tmax;}
	if (tmin < (1.0 - TOLERANCE)) {v1 += dv*tmin;}
	return 1;
}


bool check_line_clip_expand(point const &v1, point const &v2, float const d[3][2], float expand) {

	float de[3][2];
	UNROLL_3X(de[i_][0]=d[i_][0]-expand; de[i_][1]=d[i_][1]+expand;)
	return check_line_clip(v1, v2, de);
}


// ************ PROJECTIONS, CENTER, BOUNDING VOLUME, AND TRANSFORMS ************


void cylinder_quad_projection(point *pts, cylinder_3dw const &c, vector3d const &v1, vector3d &v2, int &npts) {

	assert(pts != NULL);
	npts = 0;
	v2   = c.p2 - c.p1;
	vector3d const crossp(cross_product(v1, v2).get_norm());
	gen_cylin_pts(pts, npts, c.p1, c.r1, crossp);
	gen_cylin_pts(pts, npts, c.p2, c.r2, crossp);
	if (npts == 4) {swap(pts[2], pts[3]);}
}


template<typename T> pointT<T> get_center_arb(pointT<T> const *const pts, int npts) {

	assert(pts != NULL && npts > 0);
	point pos(pts[0]);

	for (int i = 1; i < npts; ++i) {
		pos += pts[i]; // average to get the center
	}
	return pos/(float)npts;
}

template point   get_center_arb(point   const *const pts, int npts);
template point_d get_center_arb(point_d const *const pts, int npts);


// only_vis_outline => 6 corners, otherwise 8
unsigned get_cube_corners(float const d[3][2], point corners[8], point const &viewed_from, bool all_corners) {

	for (unsigned i = 0; i < 2; ++i) {
		for (unsigned j = 0; j < 2; ++j) {
			for (unsigned k = 0; k < 2; ++k) {
				unsigned const ix((((i<<1)+j)<<1)+k);
				corners[ix][0] = d[0][i];
				corners[ix][1] = d[1][j];
				corners[ix][2] = d[2][k];
			}
		}
	}
	if (all_corners) return 8;
	float dmin(0.0), dmax(0.0);
	unsigned cmin(0), cmax(0);

	for (unsigned i = 0; i < 8; ++i) {
		float const dist_sq(p2p_dist_sq(viewed_from, corners[i]));

		if (i == 0 || dist_sq < dmin) {
			dmin = dist_sq;
			cmin = i;
		}
		if (i == 0 || dist_sq >= dmax) { // >= required so that different points are chosen if all dists are equal
			dmax = dist_sq;
			cmax = i;
		}
	}
	assert(cmax != cmin);
	swap(corners[cmin], corners[7]);
	if (cmax == 7) cmax = cmin; // account for the above swap
	swap(corners[cmax], corners[6]);
	return 6;
}


void get_closest_cube_norm(float const d[3][2], point const &p, vector3d &norm) {

	unsigned dim(2), dir(0);
	float dmin(fabs(p[2] - d[2][0]));
	UNROLL_2X({float const dist(fabs(p[i_] - d[i_][0])); if (dist < dmin) {dmin = dist; dim = i_;}})
	UNROLL_3X({float const dist(fabs(p[i_] - d[i_][1])); if (dist < dmin) {dmin = dist; dim = i_; dir = 1;}})
	norm      = all_zeros;
	norm[dim] = (dir ? 1.0 : -1.0);
}


void cylinder_bounding_sphere(point const *const p, float r1, float r2, point &center, float &radius) {

	center = get_center_n2(p);
	radius = sqrt(p2p_dist_sq(center, p[0]) + 2.0*max(r1, r2)*max(r1, r2));
}


void polygon_bounding_sphere(const point *pts, int npts, float thick, point &center, float &radius) {

	center = get_center(pts, npts);
	radius = 0.0;
	for (int i = 0; i < npts; ++i) {radius = max(radius, p2p_dist_sq(center, pts[i]));}
	radius = sqrt(radius) + thick*SQRT3; // not sure this thickness thing is quite right...
}


void add_rotated_quad_pts(vert_norm *points, unsigned &ix, float theta, float z, point const &pos, float xscale1, float xscale2, float yscale, float zscale) {

	point pts[4];
	vector3d const v1(SINF(theta), COSF(theta), 0.0), v2(cross_product(v1, plus_z)); // should be normalized
	pts[0] = 2.0*yscale*v1 + xscale2*v2; pts[0].z += z - zscale;
	pts[1] = 2.0*yscale*v1 - xscale2*v2; pts[1].z += z - zscale;
	pts[2] = -xscale1*v2; pts[2].z += z + zscale;
	pts[3] =  xscale1*v2; pts[3].z += z + zscale;
	vector3d const norm(cross_product((pts[1] - pts[0]), (pts[2] - pts[1])).get_norm());
	for (unsigned i = 0; i < 4; ++i) {points[ix++] = vert_norm((pts[i] + pos), norm);}
}


// n is already normalized (unused)
void vproj_plane(vector3d const &vin, vector3d const &n, vector3d &vout) { // project vector v onto plane with normal n

	double const m[3][3] = {
		{ n.y*n.y+n.z*n.z, -n.x*n.y,         -n.x*n.z},
		{-n.y*n.x,          n.x*n.x+n.z*n.z, -n.y*n.z},
		{-n.z*n.x,         -n.z*n.y,          n.x*n.x+n.y*n.y}
	};
	matrix_mult(vin, vout, m);
}


// OK if vrot is zero_vector, sin = sqrt(1 - cos^2) ???
// Note: angle is in radians
#define CREATE_ROT_MATRIX(vrot, angle) \
	double const mag(vrot.mag()); assert(mag != 0.0); \
	double const X(vrot.x/mag), Y(vrot.y/mag), Z(vrot.z/mag), c(cos(angle)), s(sin(angle)); \
	double const t(1.0 - c), tX(t*X), tY(t*Y); \
	double const m[3][3] = { \
		{tX*X + c,    tX*Y + s*Z,  tX*Z  - s*Y}, \
		{tX*Y - s*Z,  tY*Y + c,    tY*Z  + s*X}, \
		{tX*Z + s*Y,  tY*Z - s*X,  t*Z*Z + c  }, \
	};


// Note: making vin by const reference is bad (fails ship weapon intersection assertion)
template<typename T> void rotate_vector3d(pointT<T> vin, pointT<T> const &vrot, double angle, pointT<T> &vout) { // rotate vin by angle about vrot to get vout

	if (angle == 0.0) return;
	CREATE_ROT_MATRIX(vrot, angle);
	matrix_mult(vin, vout, m);
}


template<typename T> void rotate_vector3d_multi(pointT<T> const &vrot, double angle, pointT<T> *vout, unsigned nv) {

	assert(vout != NULL);
	if (angle == 0.0) return;
	CREATE_ROT_MATRIX(vrot, angle);
	
	for (unsigned i = 0; i < nv; ++i) {
		pointT<T> const vin(vout[i]); // have to cache this
		matrix_mult(vin, vout[i], m);
	}
}


// vrot must be normalized
void rotate_vector3d_x2(point const &vrot, double angle, point &vout1, point &vout2) { // rotate by angle about vrot

	CREATE_ROT_MATRIX(vrot, angle);
	vector3d const vin1(vout1), vin2(vout2); // have to cache these
	matrix_mult(vin1, vout1, m);
	matrix_mult(vin2, vout2, m);
}

template void rotate_vector3d(vector3d vin, vector3d const &vrot, double angle, vector3d &vout);
template void rotate_vector3d(vector3d_d vin, vector3d_d const &vrot, double angle, vector3d_d &vout);
template void rotate_vector3d_multi(vector3d const &vrot, double angle, vector3d *vout, unsigned nv);
template void rotate_vector3d_multi(vector3d_d const &vrot, double angle, vector3d_d *vout, unsigned nv);


// apply the same rotation to vout that is required to rotate v1 to v2
void rotate_vector3d_by_vr(vector3d v1, vector3d v2, vector3d &vout) { // v1 rotated by vout = v2

	v1.normalize();
	v2.normalize();
	vector3d const v(cross_product(v2, v1));
	double const c(dot_product(v1, v2));
	if (abs(c + 1.0) < TOLERANCE) return; // v1 and v2 are parallel
	double const t(1.0/(1.0+c)), tX(t*v.x), tY(t*v.y);
	double const m[3][3] = {
		{tX*v.x + c,    tX*v.y + v.z,  tX*v.z    - v.y},
		{tX*v.y - v.z,  tY*v.y + c,    tY*v.z    + v.x},
		{tX*v.z + v.y,  tY*v.z - v.x,  t*v.z*v.z + c},
	};
	vector3d const vin(vout);
	matrix_mult(vin, vout, m);
}

cube_t rotate_cube(cube_t const &cube, vector3d const &axis, float angle_in_radians) {

	point pts[2] = {cube.get_cube_center(), 0.5*cube.get_size()}; // {center, extent}
	rotate_vector3d_multi(axis, angle_in_radians, pts, 2);
	return cube_t(pts[0]-pts[1], pts[0]+pts[1]);
}


// unused
float angle_of_projected_vectors(vector3d const &v1, vector3d const &v2, vector3d n) {

	vector3d vp1, vp2;
	n.normalize();
	vproj_plane(v1, n, vp1);
	vproj_plane(v2, n, vp2);
	return get_norm_angle(vp1, vp2); // angle from 'v1' to 'v2' about 'n'
}


vector3d rtp_to_xyz(float radius, double theta, double phi) {

	double const msin_phi(radius*sin(phi));
	return vector3d(cos(theta)*msin_phi, sin(theta)*msin_phi, radius*cos(phi));
}


// ************ MISC ************


template <rand_func rfunc> vector3d gen_rand_vector_template(float mag, float zscale, float phi_term) {

	float phi;

	if (phi_term == PI || phi_term == TWO_PI) {
		phi = gen_rand_phi<rfunc>();
		if (phi_term == PI) phi = abs(phi);
	}
	else {
		phi = rfunc(0.0, phi_term);
	}
	vector3d v(rtp_to_xyz(mag, rfunc(0.0, TWO_PI), phi));
	v.z *= zscale;
	return v;
}


vector3d gen_rand_vector_uniform(float mag) {

	return rtp_to_xyz(mag*rand_float(), rand_uniform(0.0, TWO_PI), gen_rand_phi<rand_uniform>());
}


vector3d gen_rand_vector(float mag, float zscale, float phi_term) {

	return gen_rand_vector_template<rand_uniform>(mag, zscale, phi_term);
}


vector3d gen_rand_vector2(float mag, float zscale, float phi_term) {

	return gen_rand_vector_template<rand_uniform2>(mag, zscale, phi_term);
}


vector3d lead_target(point const &ps, point const &pt, vector3d const &vs, vector3d const &vt, float vweap) {

	point pt0(pt, ps); // make ps the reference orgin
	vector3d const vt0(vt, vs); // make vs the reference velocity
	if (vt0.mag_sq() < TOLERANCE*TOLERANCE) return pt0.get_norm(); // target is stopped, aim right at it

	// quadratic equation - solve for time t
	double const a(vt0.mag_sq() - vweap*vweap);

	if (fabs(a) < TOLERANCE || fabs(a) < 0.01*vweap*vweap) {
		// velocities are too close - prediction will be inaccurate and unstable
		// Note: when used for shooter movement, can lead to instability in AI
		//       due to switching target prediction on and off as speed difference changes
		return pt0.get_norm();
	}
	double const b(2.0*dot_product(pt0, vt0)), c(pt0.mag_sq()), val(b*b - 4.0*a*c);
	if (val < 0.0) return zero_vector; // can't hit the target
	double const t((-b - sqrt(val))/(2.0*a)); // choose minus
	if (t <= 0.0) return zero_vector; // can't hit the target

	// now use to to calculate the collision position pc0, then the direction vector
	return (pt0 + vt0*t).get_norm(); // pc0 - (0,0,0)
}


vector3d get_firing_dir(vector3d const &src, vector3d const &dest, float fvel, float gravity_scale) {

	if (fvel <= 0.0) return zero_vector; // edge cases that we don't want to deal with
	if (dest.x == src.x && dest.y == src.y) return vector3d(0.0, 0.0, (dest.z - src.z)/fabs(dest.z - src.z)); // straight up or down
	if (gravity_scale == 0.0 || base_gravity == 0.0) return (dest - src).get_norm(); // no gravity, aim directly at target
	float const gravity(gravity_scale*base_gravity*GRAVITY), a_2(-0.5*gravity); // in -z
	float const dist(p2p_dist_xy(src, dest)), height(dest.z - src.z); // ok if height is negative
	float const dsq(dist*dist), hsq(height*height), vsq(fvel*fvel), val(2*a_2*height*dsq);
	float const a(-dsq - hsq), b(vsq*dsq + 2*vsq*hsq - val); // Note: a is always negative
	float const c(-hsq*vsq*vsq + val*vsq - a_2*a_2*dsq*dsq);
	float const in_sqrt(b*b - 4*a*c); // quadratic formula
	assert(a < 0.0);
	if (in_sqrt < 0.0) return zero_vector; // no solution
	float const sqrt_part(sqrt(in_sqrt)), vz_sq1((-b - sqrt_part)/(2*a)), vz_sq2((-b + sqrt_part)/(2*a));
	float vz_sq;

	if (vz_sq1 >= 0.0 && vz_sq2 >= 0.0) { // two solutions
		vz_sq = min(vz_sq1, vz_sq2); // choose the smaller value to minimize time
	}
	else if (vz_sq1 >= 0.0) { // one solution
		vz_sq = vz_sq1;
	}
	else if (vz_sq2 >= 0.0) { // one solution
		vz_sq = vz_sq2;
	}
	else {
		return zero_vector; // no solution
	}
	assert(vz_sq >= 0.0 && vz_sq <= vsq+TOLERANCE);
	float const vz(sqrt(vz_sq)), v_xy(sqrt(vsq - vz_sq));
	//cout << "fvel: " << fvel << ", dist: " << dist << ", height: " << height << ", g: " << gravity << ", a: " << a << ", b: " << b << ", c: " << c << ", vz: " << vz << ", t: " << dist/v_xy << endl;
	vector3d dir_xy(vector3d(dest.x-src.x, dest.y-src.y, 0.0).get_norm());
	return (dir_xy*v_xy + vector3d(0.0, 0.0, vz)).get_norm();
}


