// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 4/9/02
#ifndef _INLINES_H_
#define _INLINES_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gl_includes.h"


extern int MESH_X_SIZE, MESH_Y_SIZE, world_mode, do_zoom, xoff2, yoff2;
extern float X_SCENE_SIZE, Y_SCENE_SIZE, Z_SCENE_SIZE, DX_VAL, DY_VAL;
extern float light_factor, relh_adj_tex, glaciate_exp_inv, cview_radius, czmin, czmax, zbottom, ztop;
extern point cview_dir, camera_origin, camera_pos;
extern upos_point_type cur_origin;
extern vector3d up_vector;
extern colorRGB uw_atten_max, uw_atten_scale;
extern pos_dir_up camera_pdu, player_pdu;
extern char **mesh_draw;
extern float SCENE_SIZE[];
extern rand_gen_t global_rand_gen;


// ***************** MATH FUNCTIONS ********************


inline float SIGN(float v) {return (v < 0.0) ? -1.0 : 1.0;}

inline float safe_acosf(float val) {return acos(max(-1.0f, min(1.0f, val)));}

// fast 1/sqrt(x), accurate to ~0.17% error
inline float InvSqrt(float x) {
	
	union {
		float f;
		int i;
	} tmp;

	tmp.f = x;
	tmp.i = 0x5f3759df - (tmp.i >> 1);
	float const y(tmp.f);
	return y*(1.5f - 0.5f*x*y*y); // repeat this iteration for more accuracy
}

#define INTERP_1D(v, s, t, n, i) t*(s*v[2]i + (1.0-s)*v[n-1]i) + (1.0-t)*(s*v[1]i + (1.0-s)*v[0]i)

// s => n1 - n0, t => n3 - n0
template<typename T> inline T interpolate_3d(T const *v, unsigned npts, float s, float t) {
	return T(INTERP_1D(v, s, t, npts, [0]), INTERP_1D(v, s, t, npts, [1]), INTERP_1D(v, s, t, npts, [2]));
}

// these interpolation functions take and return x in (0.0, 1.0)
inline float cubic_interpolate  (float x) {return ((3 - 2*x)*x*x);}
inline float quintic_interpolate(float x) {return ((5*x*x - 15*x + 10)*x*x*x);}

inline int round_fp(double val) {return ((val > 0.0) ? int(val + 0.5) : int(val - 0.5));}


// ***************** RANDOM NUMBER GENERATION ********************


typedef float (*rand_func)(float, float);

inline int rand2()                {return global_rand_gen.rand();}
inline int rand2_seed_mix()       {return global_rand_gen.rand_seed_mix();}
inline void rand2_mix()           {return global_rand_gen.rand_mix();}
inline double rand2d()            {return global_rand_gen.randd();}
inline float rand_float2()        {return global_rand_gen.rand_float();} // uniform 0 to 1
inline float signed_rand_float2() {return global_rand_gen.signed_rand_float();}
inline float rand_uniform2(float val1, float val2) {return global_rand_gen.rand_uniform(val1, val2);}
inline unsigned rand_uniform_uint2(unsigned min_val, unsigned max_val) {return global_rand_gen.rand_uniform_uint(min_val, max_val);}
inline float rgauss2() {return global_rand_gen.rgauss();} // mean = 0.0, std_dev = 1.0
inline float rand_gaussian2(float mean, float std_dev) {return global_rand_gen.rand_gaussian(mean, std_dev);}
inline void set_rand2_state(long rs1, long rs2)    {global_rand_gen.set_state(rs1, rs2);}
inline vector3d signed_rand_vector2(float scale=1.0)           {return global_rand_gen.signed_rand_vector(scale);}
inline vector3d signed_rand_vector2_norm(float scale=1.0)      {return global_rand_gen.signed_rand_vector_norm(scale);}
inline vector3d signed_rand_vector2_spherical(float scale=1.0) {return global_rand_gen.signed_rand_vector_spherical(scale);}

inline float rand_float()        {return 0.0001*(rand()%10000);} // uniform 0 to 1 (only 16-bit random numbers)
inline float signed_rand_float() {return 2.0*rand_float() - 1.0;}
inline float rand_uniform(float val1, float val2) {return 0.5*((val1 + val2) + fabs(val2 - val1)*signed_rand_float());}
inline float rgauss() {return gauss_rand_arr[rand ()%N_RAND_DIST];} // mean = 0.0, std_dev = 1.0
inline float rand_gaussian(float mean, float std_dev) {return mean + std_dev*rgauss();}


inline point gen_rand_scene_pos() {
	return point(X_SCENE_SIZE*signed_rand_float(),
		         Y_SCENE_SIZE*signed_rand_float(), 
				 rand_uniform(min(zbottom, czmin), max(ztop, czmax)));
}


template <rand_func rfunc> float gen_rand_phi() {
	return safe_acosf(2.0*rfunc(0.0, 1.0) - 1.0);
}


inline vector3d signed_rand_vector(float scale=1.0) {
	assert(scale > 0.0);
	return vector3d(scale*signed_rand_float(), scale*signed_rand_float(), scale*signed_rand_float());
}


inline vector3d signed_rand_vector_norm(float scale=1.0) {
	assert(scale > 0.0);

	while (1) {
		vector3d const v(signed_rand_vector());
		if (v.mag_sq() > TOLERANCE) return scale*v.get_norm();
	}
	return zero_vector; // never gets here
}


inline vector3d signed_rand_vector_spherical(float scale=1.0) {
	assert(scale > 0.0);

	while (1) {
		vector3d const v(signed_rand_vector(scale));
		if (v.mag_sq() < scale*scale) return v;
	}
	return zero_vector; // never gets here
}


inline void vadd_rand(vector3d &v, float rval, bool spherical=0) {v += (spherical ? signed_rand_vector_spherical(rval) : signed_rand_vector(rval));}


// ***************** VECTOR MATH ********************


template<typename T> inline void cross_product(pointT<T> const &v1, pointT<T> const &v2, pointT<T> &v3) {
	v3.x = v1.y*v2.z - v1.z*v2.y;
	v3.y = v1.z*v2.x - v1.x*v2.z;
	v3.z = v1.x*v2.y - v1.y*v2.x;
}


template<typename T> inline pointT<T> cross_product(pointT<T> const &v1, pointT<T> const &v2) {
	pointT<T> cp;
	cross_product(v1, v2, cp);
	return cp;
}


inline float vector_determinant(vector3d const &v1, vector3d const &v2, vector3d const &v3) { // 3x3
	return v1.x*(v2.y*v3.z - v3.y*v2.z) - v2.x*(v1.y*v3.z - v3.y*v1.z) + v3.x*(v1.y*v2.z - v2.y*v1.z);
}

inline float get_water_coll_angle(vector3d const &v) { // normal n = (0, 0, 1.0)
	return safe_acosf(-v.z/v.mag());
}

inline float pt_line_dist(point const &P, point const &L1, point const &L2) {
	vector3d const L(L2 - L1);
	float const L_mag(L.mag());
	return ((L_mag < TOLERANCE) ? p2p_dist(L1, P) : cross_product(L, (L1 - P)).mag()/L_mag);
}
inline bool pt_line_dist_less_than(point const &P, point const &L1, point const &L2, float dist) {
	vector3d const L(L2 - L1), cp(cross_product(L, (L1 - P)));
	return (cp.mag_sq() < dist*dist*L.mag_sq());
}
inline bool pt_line_dir_dist_less_than(point const &P, point const &L1, vector3d const &Ldir, float dist) { // Ldir is normalized
	return (cross_product(Ldir, (L1 - P)).mag_sq() < dist*dist);
}

template<typename T, typename S> inline float p2p_dist_sq(const pointT<T> &pt1, const pointT<S> &pt2) {
	return (pt1.x-pt2.x)*(pt1.x-pt2.x) + (pt1.y-pt2.y)*(pt1.y-pt2.y) + (pt1.z-pt2.z)*(pt1.z-pt2.z);
}
template<typename T, typename S> inline float p2p_dist(const pointT<T> &pt1, const pointT<S> &pt2) {
	return sqrt(p2p_dist_sq(pt1, pt2));
}

template<typename T, typename S> inline float p2p_dist_xy_sq(const pointT<T> &pt1, const pointT<S> &pt2) {
	return (pt1.x-pt2.x)*(pt1.x-pt2.x) + (pt1.y-pt2.y)*(pt1.y-pt2.y);
}
template<typename T, typename S> inline float p2p_dist_xy(const pointT<T> &pt1, const pointT<S> &pt2) {
	return sqrt(p2p_dist_xy_sq(pt1, pt2));
}

template<typename T, typename S, typename V> inline bool dist_less_than(pointT<T> const &pt1, pointT<S> const &pt2, V dval) {
	return (p2p_dist_sq(pt1, pt2) < dval*dval);
}
template<typename T, typename S, typename V> inline bool dist_xy_less_than(pointT<T> const &pt1, pointT<S> const &pt2, V dval) {
	return (p2p_dist_xy_sq(pt1, pt2) < dval*dval);
}

template<typename T> inline T dot_product(pointT<T> const &A, pointT<T> const &B, pointT<T> const &C) {
	return ((B.x - A.x)*(C.x - A.x) + (B.y - A.y)*(C.y - A.y) + (B.z - A.z)*(C.z - A.z));
}
template<typename T> inline T dot_product(pointT<T> const &A, pointT<T> const &B) {
	return (A.x*B.x + A.y*B.y + A.z*B.z);
}
template<typename S, typename T> inline double dot_product(pointT<S> const &A, pointT<T> const &B) {
	return (A.x*B.x + A.y*B.y + A.z*B.z);
}

template<typename T> inline T dot_product_ptv(pointT<T> const &V, pointT<T> const &A, pointT<T> const &B) {
	return (V.x*(A.x - B.x) + V.y*(A.y - B.y) + V.z*(A.z - B.z));
}
template<typename T> inline T dot_product_ptv_norm(pointT<T> const &V, pointT<T> const &A, pointT<T> const &B) {
	return (V.x*(A.x - B.x) + V.y*(A.y - B.y) + V.z*(A.z - B.z))/(V.mag()*p2p_dist(A, B));
}

template<typename T> inline void matrix_mult(pointT<T> const &vin, pointT<T> &vout, double const m[3][3]) {

	// basic V[3] = M[3][3]xV[3]
	vout[0] = T(vin[0]*m[0][0] + vin[1]*m[0][1] + vin[2]*m[0][2]);
	vout[1] = T(vin[0]*m[1][0] + vin[1]*m[1][1] + vin[2]*m[1][2]);
	vout[2] = T(vin[0]*m[2][0] + vin[1]*m[2][1] + vin[2]*m[2][2]);
}


inline vector3d get_poly_dir_norm(vector3d const &norm, point const &p1, vector3d const &v1, float t) {
	return ((dot_product_ptv(norm, p1, (p1 + v1*t)) < 0.0) ? -norm : norm);
}

inline point get_center_n2(point const *const pts) {
	return ((pts[0] + pts[1])*0.5);
}


inline vector3d get_norm_rand(vector3d const &v) {

	float const vmag(v.mag());
	if (vmag < TOLERANCE) {return signed_rand_vector_norm();}
	return v/vmag;
}


template<typename T> inline void get_normal(pointT<T> const &v1, pointT<T> const &v2, pointT<T> const &v3, pointT<T> &norm, bool normalize) {
	cross_product((v2 - v1), (v3 - v2), norm);
	if (normalize) {norm.normalize();}
}

inline void orthogonalize_dir(vector3d const &vin, vector3d const &dir, vector3d &vortho, bool normalize) {
	cross_product(dir, cross_product(vin, dir), vortho);
	if (normalize) {vortho.normalize();}
}


inline vector3d get_poly_norm(point const *const points, bool normalize) { // requires at least 3 points

	assert(points != NULL);
	vector3d norm;
	get_normal(points[0], points[1], points[2], norm, normalize);
	return norm;
}


inline bool is_axis_aligned(vector3d const &n) { // n must be normalized
	return (fabs(n.x) > 0.99 || fabs(n.y) > 0.99 || fabs(n.z) > 0.99);
}

inline bool is_triangle_valid(point const &p1, point const &p2, point const &p3) {
	return (!dist_less_than(p1, p2, TOLERANCE) && !dist_less_than(p2, p3, TOLERANCE) && !dist_less_than(p3, p1, TOLERANCE));
}

// only applies to the first 3 points (first triangle) since this corresponds to the points used in get_poly_norm()
inline bool is_poly_valid(point const *const p) {
	return is_triangle_valid(p[0], p[1], p[2]);
}


inline bool line_intersect_sphere(point const &p1, vector3d const &v12, point const &sc, float radius) {

	float rad, t, dist; // unused
	return line_intersect_sphere(p1, v12, sc, radius, rad, dist, t);
}

inline bool sphere_int_cylinder_sides(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2) {

	float t, rad; // unused
	vector3d v1, v2; // unused
	return sphere_int_cylinder_pretest(sc, sr, cp1, cp2, r1, r2, 0, v1, v2, t, rad);
}

inline bool sphere_intersect_cylinder(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2) {

	point p_int; // unused
	vector3d norm; // unused
	return sphere_intersect_cylinder_ipt(sc, sr, cp1, cp2, r1, r2, 1, p_int, norm, 0);
}

// p2 = line starting point, p1 = circle center, v1 = line direction, norm = plane normal, and r2sq = square of the circle radius
inline bool circle_test_comp(point const &p2, point const &p1, vector3d const &v1, vector3d norm, float r2sq, float &t) {

	norm.normalize();
	point pos;
	return (line_int_plane(p2, (v1 + p2), p1, norm, pos, t, 0) && p2p_dist_sq(p1, pos) < r2sq);
}

inline bool sphere_test_comp(point const &p2, point const &p1, vector3d const &v1, float r2sq) {
	float t;
	return sphere_test_comp(p2, p1, v1, r2sq, t);
}

inline bool line_sphere_intersect(point const &p1, point const &p2, point const &c, float r) {

	vector3d const v1(p1, p2);
	return sphere_test_comp(p1, c, v1, r*r);
}

inline bool line_sphere_int_cont(point const &p1, point const &p2, point const &c, float r) {

	if (dist_less_than(p1, c, r) || dist_less_than(p2, c, r)) return 1; // p1 or p2 inside sphere (c, r)
	return line_sphere_intersect(p1, p2, c, r);
}


inline int line_int_cylinder(point const &p1, point const &p2, point const &cp1, point const &cp2, float r1, float r2, bool check_ends, float &t) {
	return line_int_thick_cylinder(p1, p2, cp1, cp2, 0.0, 0.0, r1, r2, check_ends, t);
}

inline bool point_in_cylinder(point const &cp1, point const &cp2, point const &pos, float r1, float r2) {
	return sphere_intersect_cylinder(pos, 0.0, cp1, cp2, r1, r2);
}


inline bool line_poly_intersect(point const &p1, point const &p2, point const *points, unsigned npts, vector3d const &norm, float &t) {
	point pos;
	return (line_int_plane(p1, p2, points[0], norm, pos, t, 0) && planar_contour_intersect(points, npts, pos, norm));
}

inline bool line_poly_intersect(vector3d const &v1, point const &p1, point const *points, unsigned npts, vector3d const &norm) {
	float t;
	return line_poly_intersect(p1, (p1 + v1), points, npts, norm, t);
}

inline bool line_poly_intersect(vector3d const &v1, point const &p1, point const *points, unsigned npts, bool bfc=0) {

	vector3d norm;
	get_normal(points[0], points[1], points[2], norm, 0); // doesn't require norm to be normalized
	if (bfc && dot_product(v1, norm) < 0.0) return 0;
	return line_poly_intersect(v1, p1, points, npts, norm);
}


inline float get_angle(vector3d const &v1, vector3d const &v2) {
	return safe_acosf(dot_product(v1, v2));
}

inline float get_norm_angle(vector3d const &v1, vector3d const &v2) {
	return get_angle(v1.get_norm(), v2.get_norm());
}


template<typename T> inline pointT<T> get_center(const pointT<T> *pts, int npts) {

	assert(pts != NULL);
	if (npts == 3) return (pts[0] + pts[1] + pts[2]) * 0.3333333; // 1/3
	if (npts == 4) return (pts[0] + pts[1] + pts[2] + pts[3]) * 0.25; // 1/4
	return get_center_arb(pts, npts);
}


inline bool line_is_axis_aligned(point const &p1, point const &p2) {

	unsigned eq_dims(0);
	UNROLL_3X(if(fabs(p1[i_] - p2[i_]) < TOLERANCE) ++eq_dims;)
	return (eq_dims >= 2);
}


// Note: v is velocity of sphere at p2
template <typename T> void get_sphere_mov_sphere_int_pt(pointT<T> const &p1, pointT<T> const &p2, vector3d const &v, float rsum, pointT<T> &cpos) {

	pointT<T> const norm(p2, p1);
	double const vmag_sq(v.mag_sq());
	if (vmag_sq < TOLERANCE) {cpos = p2; return;}
	double const t1(-dot_product(norm, v)); // negate?
	pointT<T> const q(p2 + v*(t1/vmag_sq));
	double const d_sq(p2p_dist_sq(p1, q)), ival(rsum*rsum - d_sq), t2((ival > 0.0) ? sqrt(ival) : 0.0);
	cpos = q - v*(t2/sqrt(vmag_sq));
}


// *********************** CAMERA STUFF ************************


inline point get_camera_pos() {
	return camera_pos;
}

inline float distance_to_camera(point const &pos) {
	return p2p_dist(get_camera_pos(), pos);
}

inline float distance_to_camera_sq(point const &pos) {
	return p2p_dist_sq(get_camera_pos(), pos);
}

inline point get_camera_all() {
	return ((world_mode == WMODE_UNIVERSE) ? get_player_pos() : get_camera_pos());
}

inline vector3d get_vdir_all() {
	return ((world_mode == WMODE_UNIVERSE) ? get_player_dir() : cview_dir);
}

inline vector3d get_upv_all() {
	return ((world_mode == WMODE_UNIVERSE) ? get_player_up() : up_vector);
}

inline float get_zoom_scale() {
	return (do_zoom ? SQRT_ZOOMF : 1.0);
}

inline bool sphere_in_camera_view(point const &pos, float radius, int max_level) {
	return sphere_in_view(camera_pdu, pos, radius, max_level);
}

inline bool univ_sphere_vis(point const &pos, float radius) {
	return player_pdu.sphere_visible_test(pos, radius);
}

template<typename T> inline pointT<T> make_pt_global(pointT<T> const &pos) {
	return ((cur_origin == all_zeros) ? pos : (pos - cur_origin));
}

inline vector3d get_norm_camera_orient(vector3d const &normal, point const &center) {
	bool const inv_norm(dot_product_ptv(normal, get_camera_pos(), center) < 0.0);
	return normal*(inv_norm ? -1.0 : 1.0);
}


struct cmp_back_to_front {
	bool operator()(point const &a, point const &b) const {
		return (distance_to_camera_sq(a) > distance_to_camera_sq(b));
	}
	template<typename vert_t> bool operator()(vert_t const &a, vert_t const &b) const {
		return (distance_to_camera_sq(a.v) > distance_to_camera_sq(b.v));
	}
};


// *********************** FLOATING POINT ************************


inline int is_nan(float f) {
	return ((!(f >= 0.0) && !(f < 0.0)) || f > 1.0E30 || f < -1.0E30);
}

inline bool is_nan(vector3d const &v) {
	return (is_nan(v.x) || is_nan(v.y) || is_nan(v.z));
}

inline void check_fp_val(float val) {
	assert(val == 0.0 || fabs(val) > TOLERANCE);
}

inline void fix_fp_mag(float &v) {
	if (fabs(v) < TOLERANCE) v = 0.0;
}


// *********************** SCENE CLIP ************************


inline bool point_outside_mesh(int xpos, int ypos) {
	return (xpos < 0 || ypos < 0 || xpos >= MESH_X_SIZE || ypos >= MESH_Y_SIZE);
}

inline bool point_interior_to_mesh(int xpos, int ypos) {
	return (xpos > 0 && ypos > 0 && xpos < MESH_X_SIZE-1 && ypos < MESH_Y_SIZE-1);
}

inline bool is_over_mesh(point const &pos) {
	return (pos.x > -X_SCENE_SIZE && pos.x < X_SCENE_SIZE && pos.y > -Y_SCENE_SIZE && pos.y < Y_SCENE_SIZE);
}


inline bool check_line_clip(point const &v1, point const &v2, float const d[3][2]) {

	float tmin, tmax;
	return get_line_clip(v1, v2, d, tmin, tmax);
}

inline bool do_line_clip_scene(point &v1, point &v2, float z1, float z2) {

	float const d[3][2] = {{-X_SCENE_SIZE, X_SCENE_SIZE}, {-Y_SCENE_SIZE, Y_SCENE_SIZE}, {z1, z2}};
	return do_line_clip(v1, v2, d);
}


inline int get_region(point const &v, float const d[3][2]) {

	int region(0);
	if (v.x < d[0][0]) region |= 0x01; else if (v.x >= d[0][1]) region |= 0x02;
	if (v.y < d[1][0]) region |= 0x04; else if (v.y >= d[1][1]) region |= 0x08;
	if (v.z < d[2][0]) region |= 0x10; else if (v.z >= d[2][1]) region |= 0x20;
	return region;
}


// ****************** MISC GL ************************


inline void blend_color(colorRGB &C, const colorRGB &A, const colorRGB &B, float mix) {
	UNROLL_3X(C[i_] = mix*A[i_] + (1.0 - mix)*B[i_];);
}

inline void blend_color(colorRGBA &C, const colorRGBA &A, const colorRGBA &B, float mix, int calc_alpha) {

	UNROLL_3X(C[i_] = mix*A[i_] + (1.0 - mix)*B[i_];);
	if (calc_alpha) C[3] = mix*A[3] + (1.0 - mix)*B[3];
}

inline colorRGBA blend_color(const colorRGBA &A, const colorRGBA &B, float mix, int calc_alpha) {

	colorRGBA C(0.0, 0.0, 0.0, A.alpha);
	blend_color(C, A, B, mix, calc_alpha);
	return C;
}


inline void translate_to(point const &p) {
	fgTranslate(p.x, p.y, p.z);
}

template<typename T> inline void global_translate(pointT<T> const &pos) {
	translate_to(make_pt_global(pos));
}

inline void rotate_about(float angle, vector3d const &v) { // Note: angle is in degrees
	fgRotate(angle, v.x, v.y, v.z);
}

inline void rotate_by_vector(vector3d const &dir, float vadd=0.0) {
	fgRotate((-TO_DEG*safe_acosf(-dir.z) + vadd), -dir.y, dir.x, 0.0);
}

inline void scale_by(vector3d const &scale) {
	fgScale(scale.x, scale.y, scale.z);
}

inline void uniform_scale(float scale) {
	fgScale(scale, scale, scale);
}


template<typename T> inline void rotate_vector3d(pointT<T> const &vrot, double angle, pointT<T> &vout) { // rotate vout by angle (radians) about vrot
	rotate_vector3d(vout, vrot, angle, vout);
}

template<typename T> inline void rotate_vector3d_norm(pointT<T> const &vrot, double angle, pointT<T> &vout) { // rotate vout by angle (radians) about vrot

	rotate_vector3d(vout, vrot, angle, vout);
	vout.normalize();
}

inline void rotate_from_v2v(vector3d const &v1, vector3d const &v2) {
	rotate_about(TO_DEG*get_norm_angle(v1, v2), cross_product(v2, v1));
}

inline void rotate_into_plus_z(vector3d const &v) {
	rotate_about(TO_DEG*safe_acosf(v.z/v.mag()), vector3d(-v.y, v.x, 0.0)); // rotate_from_v2v(v, plus_z)
}

inline void rotate_to_dir(vector3d const &dir, float vadd, float vmult=1.0) {
	fgRotate(TO_DEG*vmult*atan2(dir.y, dir.x) + vadd, 0.0, 0.0, 1.0);
}


inline void fix_nsides(int &nsides) { // use discrete steps

	assert(nsides >= 2);
	if (nsides > 4)  nsides &= ~1;
	if (nsides > 8)  nsides &= ~3;
	if (nsides > 16) nsides &= ~7;
}


inline int get_min_dim(vector3d const &v) {
	return ((fabs(v[0]) < fabs(v[1])) ? ((fabs(v[0]) < fabs(v[2])) ? 0:2) : ((fabs(v[1]) < fabs(v[2])) ? 1:2));
}
inline int get_max_dim(vector3d const &v) {
	return ((fabs(v[0]) > fabs(v[1])) ? ((fabs(v[0]) > fabs(v[2])) ? 0:2) : ((fabs(v[1]) > fabs(v[2])) ? 1:2));
}
inline int get_min_dim(cube_t const &c) {return get_min_dim(c.get_size());}


inline int get_light()          {return ((light_factor >= 0.5) ? LIGHT_SUN : LIGHT_MOON);}
inline int get_specular_light() {return ((light_factor >= 0.4) ? LIGHT_SUN : LIGHT_MOON);} // sun takes priority

inline point get_light_pos(int light=-1) {

	point lpos;
	get_light_pos(lpos, ((light >= 0) ? light : get_light()));
	return lpos;
}


inline void water_color_atten(float *c, int x, int y, point const &pos) {
	water_color_atten_pt(c, x, y, pos, get_camera_pos(), get_light_pos());
}


inline void atten_by_water_depth(float *c, float dist) {
	UNROLL_3X(c[i_] *= (1.0 - min(uw_atten_max[i_], uw_atten_scale[i_]*dist));)
	//UNROLL_3X(c[i_] *= max(1.0f-m[i_], exp(-s[i_]*dist));)
}


inline float get_rel_height(float zval, float zmin0, float zmax0) {

	float const zv(relh_adj_tex + (zval - zmin0)/(zmax0 - zmin0));
	return ((zv > 0.0) ? pow(zv, glaciate_exp_inv) : 0.0);
}


inline bool is_mesh_disabled(int xpos, int ypos) {

	int const x(xpos + xoff2), y(ypos + yoff2);
	return (mesh_draw != NULL && !point_outside_mesh(x, y) && !mesh_draw[y][x]);
}


inline int add_coll_cylinder(cylinder_3dw const &c, cobj_params const &cparams, int platform_id=-1, int dhcm=0) {
	return add_coll_cylinder(c.p1, c.p2, c.r1, c.r2, cparams, platform_id, dhcm);
}


inline void gen_cylin_pts(point *pts, int &npts, point const &p, float radius, vector3d const &v) {

	if (radius == 0.0) {
		pts[npts++] = p;
	}
	else {
		pts[npts++] = p + v*radius;
		pts[npts++] = p - v*radius;
	}
}


// ****************** matrix allocation/deletion/clearing ************************


template<typename T> void matrix_base_alloc_2d(T **&data, unsigned nx, unsigned ny) {
	data    = new T*[ny];
	data[0] = new T[nx*ny];
}

template<typename T> void matrix_ptr_fill_2d(T **&data, unsigned nx, unsigned ny) {
	for (unsigned i = 0; i < ny; ++i) {data[i] = data[0] + i*nx;}
}

template<typename T> void matrix_gen_2d(T **&data, unsigned nx, unsigned ny) {

	//matrix_delete_2d(data); // make sure to delete before re-allocating (if only everything was initialized to NULL)
	matrix_base_alloc_2d(data, nx, ny);
	matrix_ptr_fill_2d(data, nx, ny);
}

template<typename T> void matrix_gen_2d(T **&data) {
	matrix_gen_2d(data, MESH_X_SIZE, MESH_Y_SIZE);
}

template<typename T> void matrix_gen_3d(T ***&data, unsigned nz) {
	data = new T**[nz];
	for (unsigned i = 0; i < nz; ++i) {matrix_gen_2d(data[i]);}
}

template<typename T> void matrix_delete_2d(T **&data) {

	if (data) {
		delete [] data[0];
		delete [] data;
		data = NULL;
	}
}

template<typename T> void matrix_delete_3d(T ***&data, unsigned nz) {

	if (data) {
		for (unsigned i = 0; i < nz; ++i) {matrix_delete_2d(data[i]);}
		delete [] data;
	}
}

template<typename T> void matrix_clear_1d(T *data) {
	memset(data, 0, XY_MULT_SIZE*sizeof(T));
}

template<typename T> void matrix_clear_2d(T **data) {

	assert(data);
	matrix_clear_1d(data[0]);
}

template<typename T> void remove_excess_cap(vector<T> &v) {
	if (v.size() < v.capacity()) {vector<T>(v).swap(v);}
}


// string converters

template<typename T> inline std::string make_string(T const val) {
	std::ostringstream oss;
	oss << val;
	return oss.str();
}


template<typename T> void set_bit_flag_to(T &flags, unsigned mask, bool val) {
	if (val) {flags |= mask;} else {flags &= ~mask;}
}

#endif

