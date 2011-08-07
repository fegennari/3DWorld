// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/10/02

#ifndef _3DWORLD_H_
#define _3DWORLD_H_

#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>

 // STL include
#include <vector>
#include <deque>
#include <algorithm>
#include <set>
#include <map>
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
using std::vector;
using std::deque;
using std::set;
using std::map;
using std::swap;
using std::pair;
using std::make_pair;
using std::cout;
using std::endl;
using std::min;
using std::max;

#ifndef PI
#define PI 3.141592654
#endif

int      const MIX_BLOOD_WATER  = 1;
int      const CAMERA_ID        = -1;
int      const NO_SOURCE        = -2;
unsigned const INIT_CCELL_SIZE  = 16;
float    const LARGE_OBJ_RAD    = 0.01;
float    const TOLERANCE        = 1.0E-12;
float    const CAMERA_RADIUS    = 0.06;
float    const ABSOLUTE_ZERO    = -273; // in degrees C
float    const MAX_SPLASH_DEPTH = 0.1;
float    const WATER_INDEX_REFRACT = 1.333;
float    const WATER_COL_ATTEN  = 0.6;

unsigned const TICKS_PER_SECOND = 40;

int      const MESH_BORDER_OBJ  = 1;

float    const MESH_HEIGHT0     = 0.1;
float    const OCEAN_DEPTH0     = 0.05;
float    const HEIGHT_SCALE     = 100.0;
float    const ISLAND_BOT_OFF   = 0.15;
float    const SMALL_NUMBER     = 0.001;
float    const MIN_POLY_THICK   = 0.001;
float    const LARGE_DIST       = 200.0;
unsigned const MAX_CHARS        = 256;

float    const LIGHT_ROT_AMT    = 0.05;
float    const WIND_ADJUST      = 0.2;

float    const DEF_TIMESTEP     = 0.007;
int      const TIMESCALE        = 1; // (integer) larger makes surface movement faster
int      const MAX_I_TIMESCALE  = 8; // (integer) larger makes surface movement slower (max inverse of TIMESCALE)
float    const GRAVITY          = 300.0;
float    const PIXEL_WIND_VEL   = 10.0;

float    const CLOUD_CEILING0   = 1.5;
int      const LITNING_TIME     = 50;
float    const MEL_EVAP_RATIO   = 0.75;
float    const LITNING_LINEAR_I = 1.0;

float    const STICK_THRESHOLD  = 1.0;
float    const COLL_DAMAGE      = 0.75;
float    const WATER_DAMAGE     = 5.0;
float    const WATER_FRAME_DAM  = 1.0;
float    const TEMP_INCREMENT   = 10.0;
float    const RAIN_MIN_TEMP    = 2.0; // degrees C
float    const SNOW_MAX_TEMP    = -2.0;
float    const WATER_MAX_TEMP   = 100.0; // boiling point
float    const DEF_TEMPERATURE  = 20.0;

float    const SNOW_ACC         = 10.0;
float    const W_FREEZE_POINT   = -0.1;
float    const WATER_DENSITY    = 1.0;
float    const SPLASH_BASE_SZ   = 0.01;

float    const DEF_AMBIENT      = 0.5;
float    const DEF_DIFFUSE      = 0.9;
float    const SKY_REPEATS      = 2.0;
float    const SKY_HEIGHT       = 12.0;
float    const RAIN_TAIL_MIN_V  = -1.0;

float    const WATER_ALPHA      = 0.75;
float    const ICE_ALPHA        = 0.88;
float    const SPLASH_ALPHA     = 1.0;
float    const SNOW_ALPHA       = 1.0;
float    const TRANS_AMOUNT     = 0.4;
float    const SHADOW_COLOR     = 0.5;

int      const N_RAND_GAUSS     = 10;
int      const N_RAND_DIST      = 10000;
int      const N_STAR_POINTS    = 5;
int      const N_CYL_SIDES      = 32;
int      const N_SPHERE_DIV     = 32;
int      const N_COLL_POLY_PTS  = 4;
int      const MAX_SPLASH_DROP  = 100;
int      const SMILEY_NCHUNKS   = 25;
int      const NUM_CHUNK_BLOCKS = 4;

float    const GROUND_SPEED     = 0.0175;
float    const SIDESTEP_SPEED   = 0.8;
float    const BACKWARD_SPEED   = 0.85;
float    const C_STEP_HEIGHT    = 0.6;
float    const NEAR_CLIP        = 0.01;
float    const FAR_CLIP         = 100.0;
float    const PERSP_ANGLE      = 60.0;
float    const ZOOM_FACTOR      = 5.0;

float    const RES_STEP         = 2.0;
float    const STEP_SIZE        = 1.0;
float    const PROJ_MAX_DSTEPS  = 2000.0;

float    const MIN_SHADOW_ALPHA = 0.5;
float    const GET_OCC_EXPAND   = 0.02;

unsigned const CELL_Z_DIVS      = 8;
unsigned const SMALL_NDIV       = 8;
int const FAST_VISIBILITY_CALC  = 1;

float const RG_NORM         = sqrt(3.0/N_RAND_GAUSS);
float const TWO_PI          = 2.0*PI;
float const PI_TWO          = PI/2.0;
float const PI_INV          = 1.0/PI;
float const SQRT2           = sqrt(2.0);
float const SQRT3           = sqrt(3.0);
float const TWO_SQRT2       = 2.0*SQRT2;
float const SQRTOFTWOINV    = 1.0/SQRT2;
float const TO_DEG          = 180.0/PI;
float const TO_RADIANS      = PI/180.0;
float const SQRT_ZOOMF      = sqrt(ZOOM_FACTOR);
float const SQRT_ZOOMF_INV  = 1.0/SQRT_ZOOMF;


#define CLIP_TO_01(x)  max( 0.0f, min(1.0f, (x)))
#define CLIP_TO_pm1(x) max(-1.0f, min(1.0f, (x)))

#define BITSHIFT_CEIL(num, bs) (((num-1) >> bs) + 1)

#define UNROLL_3X(expr) {{unsigned const i_(0); expr} {unsigned const i_(1); expr} {unsigned const i_(2); expr}}


enum {CAM_FILT_DAMAGE=0, CAM_FILT_FOG, CAM_FILT_BURN, CAMERA_FILT_BKG, CAM_FILT_UWATER};


template<typename T> struct point2d { // size = 8

	T x, y;

	point2d() {}
	point2d(T x_, T y_) : x(x_), y(y_) {}
	point2d(point2d const &a, point2d const &b) : x(a.x - b.x), y(a.y - b.y) {}
	T mag_sq() const {return (x*x + y*y);}
	T mag()    const {return sqrt(mag_sq());}
	T cp_mag(point2d const &p) const {return (x*p.y - y*p.x);}

	void normalize() {
		T const d(mag());
		if (d >= TOLERANCE) {x /= d; y /= d;}
	}
};


template<typename T> struct pointT { // size = 12 (float), 24(double)

	T x, y, z;

	pointT() {}
	pointT(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}
	pointT(pointT const &p1, pointT const &p2) : x(p1.x-p2.x), y(p1.y-p2.y), z(p1.z-p2.z) {} // take the difference (vector)
	template<typename S> pointT(S const &p) : x(p.x), y(p.y), z(p.z) {}

	void print() const {cout << x << ", " << y << ", " << z;}
	template<typename S> void operator=(pointT<S> const &p) {x = p.x; y = p.y; z = p.z;}
	bool operator==(const pointT &p) const {return (p.x == x && p.y == y && p.z == z);}
	bool operator!=(const pointT &p) const {return !operator==(p);}

	const T &operator[](unsigned i) const {
		switch(i) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: assert(0);
		}
		return x; // never gets here
	}
	T &operator[](unsigned i) {
		switch(i) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: assert(0);
		}
		return x; // never gets here
	}
	void assign(T x_, T y_, T z_)    {x = x_; y = y_; z = z_;}
	void operator+=(pointT const &p) {x += p.x; y += p.y; z += p.z;}
	void operator-=(pointT const &p) {x -= p.x; y -= p.y; z -= p.z;}
	void operator*=(double m)        {x *= m; y *= m; z *= m;}

	void operator/=(double d) {
		assert(d != 0.0);
		T const m(1.0/d);
		x *= m; y *= m; z *= m;
	}
	void normalize() {
		T const d(mag());
		if (d >= TOLERANCE) {
			T const m(1.0/d);
			x *= m; y *= m; z *= m;
		}
	}
	bool normalize_test() {
		T const d(mag());
		if (d < TOLERANCE) return 0;
		T const m(1.0/d);
		x *= m; y *= m; z *= m;
		return 1;
	}
	void negate() {x = -x; y = -y; z = -z;}

	void invert(bool err_div_0, bool fix_div_0) {
		if (err_div_0) {
			assert(x != 0.0 && y != 0.0 && z != 0.0);
		}
		else if (fix_div_0) {
			if (x == 0.0) x = TOLERANCE;
			if (y == 0.0) y = TOLERANCE;
			if (z == 0.0) z = TOLERANCE;
		}
		x = 1.0/x; y = 1.0/y; z = 1.0/z;
	}
	pointT get_norm() const {
		T const vmag(mag());
		return ((vmag < TOLERANCE) ? *this : pointT(x/vmag, y/vmag, z/vmag));
	}
	void set_max_mag(T vmax) {
		T const vmag(mag());
		if (vmag > TOLERANCE && vmag > vmax) operator*=(vmax/vmag);
	}
	pointT operator+(pointT const &p)  const {return pointT((x+p.x), (y+p.y), (z+p.z));}
	pointT operator-(pointT const &p)  const {return pointT((x-p.x), (y-p.y), (z-p.z));}
	pointT operator*(T      const val) const {return pointT(x*val, y*val, z*val);}
	
	pointT operator/(T      const val) const {
		assert(val != 0.0);
		T const val_inv(1.0/val);
		return pointT(x*val_inv, y*val_inv, z*val_inv);
	}
	T mag_sq()    const {return (x*x + y*y + z*z);}
	T mag()       const {return sqrt(mag_sq());}
	T xy_mag_sq() const {return (x*x + y*y);}
	T xy_mag()    const {return sqrt(xy_mag_sq());}
	bool is_nonzero() const {return (x != 0.0 || y != 0.0 || z != 0.0);}
	void do_glVertex() const;
	void do_glNormal() const;

	bool operator<(pointT const &p) const { // needed for maps and stuff
		if (z > p.z) return 1; // greater than operation?
		if (z < p.z) return 0;
		if (y < p.y) return 1;
		if (y > p.y) return 0;
		return (x < p.x);
	}
};

template<> inline void pointT<float >::do_glVertex() const {glVertex3fv((float *)this);}
template<> inline void pointT<float >::do_glNormal() const {glNormal3fv((float *)this);}

template<> inline void pointT<double>::do_glVertex() const {glVertex3dv((double *)this);}
template<> inline void pointT<double>::do_glNormal() const {glNormal3dv((double *)this);}


typedef pointT<float>  point;
typedef pointT<float>  vector3d;
typedef pointT<double> point_d;
typedef pointT<double> vector3d_d;
typedef point_d        upos_point_type;


// fixed settings
point    const all_zeros(0, 0, 0);
vector3d const plus_z(0, 0, 1);
vector3d const zero_vector(0, 0, 0);
vector3d const all_ones(1, 1, 1);


template<typename T> struct plane_t : public pointT<T> { // unused
	T d;

	plane_t() {}
	plane_t(T x_, T y_, T z_, T d_) : pointT(x_, y_, z_), d(d_) {}
	plane_t(pointT<T> const &p, T d_) : pointT(p), d(d_) {}
	void print() const {cout << x << ", " << y << ", " << z << ", " << d;}
	bool operator==(const plane_t &p) const {return (p.x == x && p.y == y && p.z == z && p.d == d);}
	bool operator!=(const plane_t &p) const {return !operator==(p);}
	void assign(T x_, T y_, T z_) {x = x_; y = y_; z = z_;}
};


struct triangle {

	point pts[3];

	triangle() {}
	triangle(point const &p1, point const &p2, point const &p3) {pts[0] = p1; pts[1] = p2; pts[2] = p3;}
};


struct sphere_t {

	point pos;
	float radius;
	
	sphere_t(point const &p=all_zeros, float r=0.0) : pos(p), radius(r) {}
	bool operator==(sphere_t const &s) const {return (pos == s.pos && radius == s.radius);}
	bool operator!=(sphere_t const &s) const {return (pos != s.pos || radius != s.radius);}
	point const &get_pos() const {return pos;}
	float get_radius()     const {return radius;}
};


struct cube_t { // size = 24

	float d[3][2]; // {x,y,z},{min,max}

	cube_t() {}
	
	cube_t(float x1, float x2, float y1, float y2, float z1, float z2) {
		d[0][0] = x1; d[0][1] = x2;
		d[1][0] = y1; d[1][1] = y2;
		d[2][0] = z1; d[2][1] = z2;
	}
	void copy_from(cube_t const &c) {
		UNROLL_3X(d[i_][0] = c.d[i_][0];)
		UNROLL_3X(d[i_][1] = c.d[i_][1];)
	}
	void translate(point const &p) {
		UNROLL_3X(d[i_][0] += p[i_]; d[i_][1] += p[i_];)
	}
	void set_from_points(point const *const pts, unsigned npts);
	void print() const;
	bool is_near_zero_area() const;

	void union_with_pt(point const &pt) {
		UNROLL_3X(d[i_][0] = min(d[i_][0], pt[i_]); d[i_][1] = max(d[i_][1], pt[i_]);)
	}
	void union_with_cube(cube_t const &c) {
		UNROLL_3X(d[i_][0] = min(d[i_][0], c.d[i_][0]); d[i_][1] = max(d[i_][1], c.d[i_][1]);)
	}
	void normalize() {
		UNROLL_3X(if (d[i_][1] < d[i_][0]) swap(d[i_][0], d[i_][1]);)
	}
	bool is_zero_area() const {
		UNROLL_3X(if (d[i_][0] == d[i_][1]) return 1;)
		return 0;
	}
	bool intersects(const cube_t &cube, float toler) const {
		UNROLL_3X(if (cube.d[i_][1] < (d[i_][0] + toler) || cube.d[i_][0] > (d[i_][1] - toler)) return 0;)
		return 1;
	}
	bool contains_cube(const cube_t &cube) const {
		UNROLL_3X(if (cube.d[i_][0] < d[i_][0] || cube.d[i_][1] > d[i_][1]) return 0;)
		return 1;
	}
	bool contains_pt(point const &pt) const {
		UNROLL_3X(if (pt[i_] < d[i_][0] || pt[i_] > d[i_][1]) return 0;)
		return 1;
	}
	bool contains_pt_xy(point const &pt) const {
		return (pt.x > d[0][0] && pt.x < d[0][1] && pt.y > d[1][0] && pt.y < d[1][1]);
	}
	bool quick_intersect_test(const cube_t &cube) const {
		UNROLL_3X(if (cube.d[i_][0] >= d[i_][1] || cube.d[i_][1] <= d[i_][0]) return 0;)
		return 1;
	}
	float get_volume() const {
		return fabs(d[0][1] - d[0][0])*fabs(d[1][1] - d[1][0])*fabs(d[2][1] - d[2][0]);
	}
	float max_len() const {
		return max((d[0][1] - d[0][0]), max((d[1][1] - d[1][0]), (d[2][1] - d[2][0])));
	}
	float min_len() const {
		return min((d[0][1] - d[0][0]), min((d[1][1] - d[1][0]), (d[2][1] - d[2][0])));
	}
	float second_largest_len() const {
		return min(max((d[0][1] - d[0][0]), (d[1][1] - d[1][0])),
			   min(max((d[1][1] - d[1][0]), (d[2][1] - d[2][0])),
			       max((d[2][1] - d[2][0]), (d[0][1] - d[0][0]))));
	}
	point get_cube_center() const {
		return point(0.5*(d[0][0]+d[0][1]), 0.5*(d[1][0]+d[1][1]), 0.5*(d[2][0]+d[2][1]));
	}
	float get_bsphere_radius() const {
		return 0.5*sqrt((d[0][1]-d[0][0])*(d[0][1]-d[0][0]) + (d[1][1]-d[1][0])*(d[1][1]-d[1][0]) + (d[2][1]-d[2][0])*(d[2][1]-d[2][0]));
	}
	void expand_by(float val) {UNROLL_3X(d[i_][0] -= val; d[i_][1] += val;)}
	bool cube_intersection(const cube_t &cube, cube_t &res) const;
	float get_overlap_volume(const cube_t &cube) const;
	vector3d closest_side_dir(point const &pos) const;
	point gen_rand_pt_in_cube() const;
};


struct line_3dw {

	point p1, p2;

	line_3dw() {}
	line_3dw(point const &p1_, point const &p2_) : p1(p1_), p2(p2_) {assert(p1 != p2);}
	vector3d get_norm_dir_vect() const {return (p2 - p1).get_norm();}
};


struct vector_point_norm {

	vector<point>    p;
	vector<vector3d> n;
};


struct pos_dir_up { // defines a view frustum

	point pos;
	vector3d dir, upv, upv_, cp;
	float tterm, sterm, tterm_sq2_inv, near_, far_;
	double A;
	bool valid;

	pos_dir_up(void) : valid(0) {}
	pos_dir_up(point const &p, vector3d const &d, vector3d const &u, float t, float s, float n, float f, float a=0.0);
	bool sphere_visible_test(point const &pos_, float radius) const;
	bool cube_visible(cube_t const &cube) const;
	void draw_frustum() const;
};


struct cylinder_3dw : public line_3dw { // size = 32

	float r1, r2;

	cylinder_3dw() {}
	cylinder_3dw(point const &p1_, point const &p2_, float r1_, float r2_) : line_3dw(p1_, p2_), r1(r1_), r2(r2_) {}
};


struct colorRGB { // size = 12

	float R, G, B;
	colorRGB() {}
	colorRGB(float r, float g, float b) : R(r), G(g), B(b) {}
	void set_to_val (float val) {R = G = B = val;}

	const float &operator[](unsigned i) const {
		switch(i) {
			case 0: return R;
			case 1: return G;
			case 2: return B;
			default: assert(0);
		}
		return R; // never gets here
	}
	float &operator[](unsigned i) {
		switch(i) {
			case 0: return R;
			case 1: return G;
			case 2: return B;
			default: assert(0);
		}
		return R; // never gets here
	}
	void operator*=(float val) {
		R *= val; G *= val; B *= val;
	}
	void set_valid_color() {
		R = CLIP_TO_01(R);
		G = CLIP_TO_01(G);
		B = CLIP_TO_01(B);
	}
	void from_normal(vector3d const &n) { // for normal maps
		R = 0.5*(n.x + 1.0);
		G = 0.5*(n.y + 1.0);
		B = 0.5*(n.z + 1.0);
	}
	void to_normal(vector3d &n) const { // for normal maps
		n.x = 2.0*R - 1.0;
		n.y = 2.0*G - 1.0;
		n.z = 2.0*B - 1.0;
	}
	void do_glColor() const {glColor3fv((float *)this);}
};


struct colorRGBA { // size = 16

	float red, green, blue, alpha;

	colorRGBA() {}
	colorRGBA(float r, float g, float b, float a=1.0) : red(r), green(g), blue(b), alpha(a) {}
	colorRGBA(colorRGB const &color, float a=1.0) : red(color.R), green(color.G), blue(color.B), alpha(a) {}

	bool operator==(const colorRGBA &c) const {
		return (c.red == red && c.green == green && c.blue == blue && c.alpha == alpha);
	}
	bool operator!=(const colorRGBA &c) const {
		return !operator==(c);
	}
	const float &operator[](unsigned i) const {
		switch(i) {
			case 0: return red;
			case 1: return green;
			case 2: return blue;
			case 3: return alpha;
			default: assert(0);
		}
		return red; // never gets here
	}
	float &operator[](unsigned i) {
		switch(i) {
			case 0: return red;
			case 1: return green;
			case 2: return blue;
			case 3: return alpha;
			default: assert(0);
		}
		return red; // never gets here
	}
	void assign(float R, float G, float B, float A=1.0) {
		red = R; green = G; blue = B; alpha = A;
	}
	bool operator<(const colorRGBA &c) const { // greater than operation?
		if (alpha > c.alpha) return 1; // note: alpha less than so that low alpha colors are last
		if (alpha < c.alpha) return 0;
		if (red   < c.red)   return 1;
		if (red   > c.red)   return 0;
		if (green < c.green) return 1;
		if (green > c.green) return 0;
		return (blue < c.blue);
	}
	colorRGBA operator*(float val) const {
		return colorRGBA(red*val, green*val, blue*val, alpha);
	}
	void operator*=(float val) {
		red *= val; green *= val; blue *= val;
	}
	colorRGBA operator+(colorRGBA const &c) const {
		return colorRGBA(red+c.red, green+c.green, blue+c.blue, alpha+c.alpha);
	}
	void operator+=(colorRGBA const &c) {
		red += c.red; green += c.green; blue += c.blue; alpha += c.alpha;
	}
	colorRGBA modulate_with(colorRGBA const &c) const {
		return colorRGBA(red*c.red, green*c.green, blue*c.blue, alpha*c.alpha);
	}
	void set_valid_color() {
		red   = CLIP_TO_01(red);
		green = CLIP_TO_01(green);
		blue  = CLIP_TO_01(blue);
		alpha = CLIP_TO_01(alpha);
	}
	void normalize_to_alpha_1() {
		if (alpha == 1.0) return;
		red  *= alpha; green *= alpha; blue *= alpha;
		alpha = 1.0;
	}
	bool within_thresh_of_rgb(float thresh, colorRGBA const &c) const { // no alpha check
		return ((fabs(red-c.red) + fabs(green-c.green) + fabs(blue-c.blue)) < thresh);
	}
	bool within_thresh_of_rgba(float thresh, colorRGBA const &c) const { // no alpha check
		return ((fabs(red-c.red) + fabs(green-c.green) + fabs(blue-c.blue) + fabs(alpha-c.alpha)) < thresh);
	}
	float get_luminance() const {return (red + green + blue)/3.0;}
	colorRGB rgb()  const {return colorRGB(red, green, blue);}
	bool is_valid() const {return (red >= 0 && green >= 0 && blue >= 0 && alpha >= 0 && red <= 1 && green <= 1 && blue <= 1 && alpha <= 1);}
	void print() const {cout << "R: " << red << ", G: " << green << ", B: " << blue << ", A: " << alpha;}
	void do_glColor() const {glColor4fv((float *)this);}
	void do_glColor4ubv() const;
};


struct vert_norm { // size = 24
	point v;
	vector3d n;

	vert_norm() {}
	vert_norm(point const &v_, vector3d const &n_) : v(v_), n(n_) {}
	void assign(point const &v_, vector3d const &n_) {v = v_; n = n_;}
};


struct vert_norm_tc : public vert_norm { // size = 32
	float t[2];

	vert_norm_tc() {}
	vert_norm_tc(point const &v_, vector3d const &n_, float ts, float tt) : vert_norm(v_, n_) {t[0] = ts; t[1] = tt;}
};


struct color_wrapper { // size = 4
	static int gl_type;
	unsigned char c[4]; // Note: c[3] (alpha component) is not used in all cases

	template<typename T> void set_c3(T const &c_) {UNROLL_3X(c[i_] = (unsigned char)(255.0*CLIP_TO_01(c_[i_]));)}
	void set_c4(colorRGBA const &c_) {set_c3(c_); c[3] = (unsigned char)(255.0*CLIP_TO_01(c_[3]));}
	colorRGB  get_c3() const {return colorRGB(c[0]/255.0, c[1]/255.0, c[2]/255.0);}
	colorRGBA get_c4() const {return colorRGBA(get_c3(), c[3]/255.0);}
};


struct color_wrapper_float { // size = 16
	static int gl_type;
	colorRGBA c; // Note: c[3] (alpha component) is not used in all cases

	template<typename T> void set_c3(T const &c_) {c = c_;}
	void set_c4(colorRGBA const &c_) {c = c_;}
	colorRGB  get_c3() const {return colorRGB(c.red, c.green, c.blue);}
	colorRGBA get_c4() const {return c;}
};


inline void colorRGBA::do_glColor4ubv() const {
	color_wrapper cw;
	cw.set_c4(*this);
	glColor4ubv(cw.c);
}


struct vert_color : public color_wrapper { // size = 16
	point v;

	vert_color() {}
	vert_color(point const &v_, colorRGBA const &c_)     : v(v_) {set_c4(c_);}
	vert_color(point const &v_, unsigned char const *c_) : v(v_) {c[0]=c_[0]; c[1]=c_[1]; c[2]=c_[2]; c[3]=c_[3];}
	void set_state(unsigned vbo) const;
};


struct vert_norm_tc_color : public vert_norm_tc, public color_wrapper { // size = 36
	vert_norm_tc_color() {}
	vert_norm_tc_color(point const &v_, vector3d const &n_, float ts, float tt, colorRGB const &c_)
		: vert_norm_tc(v_, n_, ts, tt) {set_c3(c_);}
	vert_norm_tc_color(point const &v_, vector3d const &n_, float ts, float tt, colorRGBA const &c_)
		: vert_norm_tc(v_, n_, ts, tt) {set_c4(c_);}
	vert_norm_tc_color(point const &v_, vector3d const &n_, float ts, float tt, unsigned char const *const c_, bool has_alpha=0)
		: vert_norm_tc(v_, n_, ts, tt) {c[0] = c_[0]; c[1] = c_[1]; c[2] = c_[2]; if (has_alpha) c[3] = c_[3];}
	void assign(point const &v_, vector3d const &n_, float ts, float tt, unsigned char const *const c_, bool has_alpha=0) {
		v = v_; n = n_; t[0] = ts; t[1] = tt; c[0] = c_[0]; c[1] = c_[1]; c[2] = c_[2]; if (has_alpha) c[3] = c_[3];
	}
	static void set_vbo_arrays(unsigned stride_mult=1);
	void set_state(unsigned stride_mult=1) const;
};


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
	unsigned get_mem() const {return (points.capacity() + lines.capacity() + triangles.capacity())*sizeof(vnc);}
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


struct valley { // size = 70

	struct spill_func { // size = 16

		short index, i, j, si, sj, spill;
		float z_over;

		spill_func() : index(-1), i(0), j(0), si(0), sj(0), spill(0), z_over(0.0) {}
		spill_func(short ix, short i_, short j_, short si_, short sj_, short s, float z)
			: index(ix), i(i_), j(j_), si(si_), sj(sj_), spill(s), z_over(z) {}
	};

	short x, y, spill_index;
	bool has_spilled;
	float w_volume, spill_vol, lwv, zval, min_zval, dz, area, fvol, depth, blood_mix, mud_mix, spill_integral;
	spill_func sf;

	valley(short x_=0, short y_=0) : x(x_), y(y_), spill_index(-1), has_spilled(0),
		w_volume(0.0), spill_vol(0.0), lwv(0.0), zval(0.0), min_zval(0.0), dz(0.0),
		area(0.0), fvol(0.0), depth(0.0), blood_mix(0.0), mud_mix(0.0), spill_integral(0.0) {}
	void copy_state_from(valley const &v);
	void create(int wsi);
	float get_volume() const;
};


struct valley_w { // size = 8

	short wsi, x, y, inside8;
};


struct surf_adv { // size = 4

	short x, y;
	void assign(int x_, int y_) {x = x_; y = y_;}
};


class line3d { // size = 28

public:
	float width;
	vector<point> points;
	colorRGBA color;

	line3d() {}
	void draw() const;
	void destroy();
};


struct lightning { // size = 40

	int time, enabled;
	point start, end;
	vector<line3d> path;
	typedef unsigned long long cell_ix_t;
	set<cell_ix_t> cells_seen;

	void gen();
	void gen_recur(point &start, float strength, int xpos, int ypos, int zpos, float zval, int l_frame_counter);
	void draw() const;
	
	cell_ix_t get_cell_ix(unsigned const x, unsigned const y, unsigned const z) const {
		return (x + (cell_ix_t(y) << 16) + (cell_ix_t(z) << 32)); // x, y, z < 2^16
	}
};


colorRGBA const DEF_TEX_COLOR(0.0, 0.0, 0.0, 0.0); // black


struct texture { // size = 78 (80)

	char type, format, use_mipmaps;
	bool wrap;
	int width, height, ncolors;
	unsigned char *data, *orig_data, *colored_data, *mm_data;
	std::string name;
	GLuint tid;
	colorRGBA color;
	vector<unsigned> mm_offsets;

	texture(char t, char f, int w, int h, bool wra, int nc, int um, std::string const &n,
		GLuint tex=0, colorRGBA const &c=DEF_TEX_COLOR) : type(t), format(f), use_mipmaps(um), wrap(wra),
		width(w), height(h), ncolors(nc), data(NULL), orig_data(NULL), colored_data(NULL), mm_data(0),
		name(n), tid(tex), color(c) {}
	void init();
	void do_gl_init();
	GLenum calc_internal_format() const;
	GLenum calc_format() const;
	void calc_color();
	void build_mipmaps();
	void create_custom_mipmaps();
	unsigned char const *get_mipmap_data(unsigned level) const;
	void set_to_color(colorRGBA const &c);
	void alloc();
	void free_mm_data();
	void free();
	void gl_delete();
};


struct camera_filter {

	int tid;
	unsigned time;
	colorRGBA color;

	camera_filter() : time(0) {}
	camera_filter(colorRGBA const &c, unsigned t, int tid_) : tid(tid_), time(t), color(c) {}
	void draw();
};


struct portal {

	point pts[4]; // quads only, for now
	void draw() const;
	point get_center_pt() const {return (pts[0] + pts[1] + pts[2] + pts[3])*0.25;}
	bool is_visible() const;
};


struct user_waypt_t {

	int type;
	point pos;
	user_waypt_t(int type_=0, point const &pos_=all_zeros) : type(type_), pos(pos_) {}
};


// colors
colorRGBA const RED      (1.0,  0.0,  0.0,  1.0);
colorRGBA const GREEN    (0.0,  1.0,  0.0,  1.0);
colorRGBA const BLUE     (0.0,  0.0,  1.0,  1.0);
colorRGBA const BLACK    (0.0,  0.0,  0.0,  1.0);
colorRGBA const WHITE    (1.0,  1.0,  1.0,  1.0);
colorRGBA const CYAN     (0.0,  1.0,  1.0,  1.0);
colorRGBA const MAGENTA  (1.0,  0.0,  1.0,  1.0);
colorRGBA const YELLOW   (1.0,  1.0,  0.0,  1.0);
colorRGBA const PURPLE   (0.5,  0.0,  0.6,  1.0);
colorRGBA const TRANS    (1.0,  1.0,  1.0,  TRANS_AMOUNT);
colorRGBA const ALPHA0   (1.0,  1.0,  1.0,  0.0);
colorRGBA const RGBA0    (0.0,  0.0,  0.0,  0.0);
colorRGBA const ALPHA0_5 (1.0,  1.0,  1.0,  0.5);
colorRGBA const BROWN    (0.6,  0.25, 0.1,  1.0);
colorRGBA const DK_BROWN (0.3,  0.15, 0.08, 1.0);
colorRGBA const LT_BROWN (0.6,  0.4,  0.2,  1.0);
colorRGBA const LT_RED   (1.0,  0.58, 0.58, 1.0);
colorRGBA const DK_RED   (0.7,  0.0,  0.0,  1.0);
colorRGBA const LT_GREEN (0.58, 1.0,  0.58, 1.0);
colorRGBA const DK_GREEN (0.0,  0.7,  0.0,  1.0);
colorRGBA const OLIVE    (0.3,  0.4,  0.2,  1.0);
colorRGBA const LT_BLUE  (0.58, 0.58, 1.0,  1.0);
colorRGBA const DK_BLUE  (0.0,  0.0,  0.7,  1.0);
colorRGBA const ORANGE   (1.0,  0.5,  0.0,  1.0);
colorRGBA const PINK     (1.0,  0.5,  0.5,  1.0);
colorRGBA const GRAY     (0.5,  0.5,  0.5,  1.0);
colorRGBA const LT_GRAY  (0.75, 0.75, 0.75, 1.0);
colorRGBA const DK_GRAY  (0.25, 0.25, 0.25, 1.0);
colorRGBA const GRAY_BLACK(0.1, 0.1,  0.1,  1.0);
colorRGBA const BKGRAY   (0.05, 0.05, 0.05, 1.0);
colorRGBA const SILVER   (0.15, 0.15, 0.15, 1.0);
colorRGBA const WATER    (0.35, 0.35, 1.0,  1.0);
colorRGBA const SUN_C    (1.0,  0.8,  0.5,  1.0);
colorRGBA const SUN_LT_C (1.0,  1.0,  0.85, 1.0);
colorRGBA const MOON_C   (0.9,  0.8,  0.7,  1.0);
colorRGBA const LEAF_C   (0.0,  0.6,  0.0,  1.0);
colorRGBA const LITN_C   (0.8,  0.9,  1.0,  1.0);
colorRGBA const WATER_C  (0.4,  0.4,  1.0,  WATER_ALPHA);
colorRGBA const ICE_C    (0.65, 0.65, 1.0,  ICE_ALPHA);
colorRGBA const SPLASH_C (0.7,  0.7,  1.0,  SPLASH_ALPHA);
colorRGBA const BLOOD_C  (0.4,  0.0,  0.0,  1.0);
colorRGBA const CLOUD_C  (0.9,  0.9,  0.9,  1.0);
colorRGBA const TREE_C   (0.7,  0.4,  0.2,  1.0);
colorRGBA const PTREE_C  (0.4,  0.3,  0.25, 1.0);
colorRGBA const GOLD     (0.9,  0.65, 0.05, 1.0); // more like yellow
colorRGBA const MUD_C    (0.3,  0.18, 0.09, 1.0);
colorRGBA const MUD_S_C  (0.7,  0.35, 0.25, 1.0);
colorRGBA const BRASS_C  (1.0,  1.0,  0.8,  1.0);
colorRGBA const BRONZE_C (1.0,  1.0,  0.6,  1.0);
colorRGBA const MED_GREEN(0.1,  0.7,  0.1,  1.0);
colorRGBA const SNOW_COLOR(1.0, 1.0,  1.4,  1.0); // slightly bluish (saturated)
colorRGBA const BACKGROUND_DAY(0.2, 0.3, 0.8, 1.0);
colorRGBA const BACKGROUND_NIGHT(BLACK);


// macros
#define RESET_TIME       int const timer1(glutGet(GLUT_ELAPSED_TIME));
#define GET_DELTA_TIME   (glutGet(GLUT_ELAPSED_TIME) - timer1)
#define PRINT_TIME(str) {cout << str << " time = " << GET_DELTA_TIME << endl;}


inline void atten_by_water_depth(float *c, float dist) {
	float const m[3] = {0.98, 0.97, 0.95};
	float const s[3] = {1.5,  0.9,  0.5 };
	UNROLL_3X(c[i_] *= (1.0 - min(m[i_], s[i_]*dist));)
	//UNROLL_3X(c[i_] *= max(1.0f-m[i_], exp(-s[i_]*dist));)
}

#define set_color(Color)   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, (const float *)&(Color))
#define set_color_a(Color) glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,             (const float *)&(Color))
#define set_color_d(Color) glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,             (const float *)&(Color))
#define set_color_e(Color) glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,            (const float *)&(Color))
#define set_obj_color(obj) set_color(object_types[obj.type].color)


// world modes
enum {WMODE_GROUND=0, WMODE_UNIVERSE, WMODE_INF_TERRAIN, NUM_WMODE};

// movement directions
enum {MOVE_FRONT, MOVE_BACK, MOVE_LEFT, MOVE_RIGHT, MOVE_STOP};

// gameplay powerups
enum {PU_DAMAGE = 0, PU_REGEN, PU_SHIELD, PU_SPEED, PU_FLIGHT, PU_INVISIBILITY, NUM_POWERUPS};

// object definitions
enum {SF_EYE = 0, SF_NOSE, SF_TONGUE, NUM_SMILEY_PARTS};

enum {RAIN = 0, SNOW, HAIL, LEAF, BALL, S_BALL, SMILEY, BLOOD, CHARRED, CHUNK,
	SFPART, ROCKET, LANDMINE, SEEK_D, STAR5, PLASMA, GRENADE, CGRENADE, SHRAPNEL, SHELLC,
	PROJC, DROPLET, WDROPLET, SAND, DIRT, ROCK, FRAGMENT, PARTICLE, HEALTH, SHIELD,
	POWERUP, WEAPON, AMMO, WA_PACK, CAMERA, PRECIP, BLAST_RADIUS, PROJECTILE, BEAM, IMPACT,
	PLASMA_LT_D, LASER, DROWNED, BURNED, FIRE, FELL, FROZEN, SUFFOCATED, CRUSHED, GASSED,
	WAYPOINT, SMOKE, DYNAM_PART, SKULL, GRASS, NUM_TOT_OBJS};

// landscapes
enum {PEAK = 0, RIDGE, MINIMA, VALLEY, PLATEAU, HILLSIDE, BOUNDARY};

// weapons
enum {W_UNARMED = 0, W_BBBAT, W_BALL, W_SBALL, W_ROCKET, W_LANDMINE, W_SEEK_D, W_STAR5, W_M16, W_SHOTGUN, W_GRENADE,
	W_LASER, W_PLASMA, W_BLADE, W_GASSER, /* non-selectable*/ W_CGRENADE};

// textures
enum {GROUND_TEX = 0, DARK_ROCK_TEX, WATER_TEX, WATER2_TEX, CLOUD_TEX, SUN_TEX, MOON_TEX, EARTH_TEX, MARBLE_TEX, SNOW_TEX,
	LEAF_TEX, WOOD_TEX, SAND_TEX, DIRT_TEX, CAMOFLAGE_TEX, GRASS_TEX, BRICK_TEX, MANHOLE_TEX, PALM_TEX, SMOKE_TEX,
	PLASMA_TEX, GEN_TEX, LANDSCAPE_TEX, TREE_END_TEX, TREE_HEMI_TEX, SHINGLE_TEX, PANELING_TEX, CBLOCK_TEX, MJ_LEAF_TEX, LIVE_OAK_TEX,
	LEAF2_TEX, LEAF3_TEX, PLANT1_TEX, PLANT2_TEX, PLANT3_TEX, PLANT4_TEX, FENCE_TEX, SKULL_TEX, RADIATION_TEX, YUCK_TEX,
	SAW_TEX, SAW_B_TEX, BLUR_TEX, SBLUR_TEX, PINE_TEX, NOISE_TEX, WOOD2_TEX, HBBRICK_TEX, PARTB_TEX, PLASTER_TEX,
	TILE_TEX, LOGO_TEX, DISINT_TEX, BLUR_TEX_INV, HSTRIPE_TEX, VSTRIPE_TEX, BCUBE_TEX, EXPLOSION_TEX, SHIP_HULL_TEX, BCUBE2_TEX,
	BCUBE_T_TEX, ROCK_SPHERE_TEX, PAPAYA_TEX, COFFEE_TEX, SMILEY_SKULL_TEX, ICE_TEX, ROCK_TEX, BLACK_TEX, WHITE_TEX, FIRE_TEX,
	CLOUD_RAW_TEX, SNOWFLAKE_TEX, BLUR_CENT_TEX, GRADIENT_TEX, GRASS_BLADE_TEX, WIND_TEX
};

// collision object destroyability
enum {NON_DEST=0, DESTROYABLE, SHATTERABLE, SHATTER_TO_PORTAL, EXPLODEABLE};

int const dodgeball_tids[] = {SKULL_TEX, RADIATION_TEX, YUCK_TEX};
unsigned const NUM_DB_TIDS(sizeof(dodgeball_tids)/sizeof(int));

// shadow mask bits
#define MESH_SHADOW      0x02
#define DYNAMIC_SHADOW   0x04 // such as objects
#define OBJECT_SHADOW    0x40

#define SHADOWED_ALL     0xCF

// light sources
enum {LIGHT_SUN = 0, LIGHT_MOON, NUM_LIGHT_SRC};

#define SUN_SHADOW       0x01
#define MOON_SHADOW      0x02
#define TREE_ONLY        0x08
#define CHECK_ALL_SHADOW 0x03


#define EF_Z1  0x01
#define EF_Z2  0x02
#define EF_Y1  0x04
#define EF_Y2  0x08
#define EF_X1  0x10
#define EF_X2  0x20
#define EF_ALL 0x3F

unsigned char const EFLAGS[3][2] = {{EF_X1, EF_X2}, {EF_Y1, EF_Y2}, {EF_Z1, EF_Z2}};

#define PRIM_DISABLED -1
#define PRIM_UNSET    -2

// forward references
class  sd_sphere_d;
struct dwobject;
class  obj_group;
struct star;
class  cobj_params;
class  coll_obj;
class  shape3d;
struct lightning;
struct color_tid_vol;
class shader_t;


// function prototypes - main (3DWorld.cpp, etc.)
bool check_gl_error(unsigned loc_id);
void enable_blend();
void disable_blend();
void set_array_client_state(bool va, bool tca, bool na, bool ca);
void set_lighted_sides(int num);
void enable_point_specular();
void disable_point_specular();
void reset_fog();
void set_light_atten(int light, float attenuation=1.0);
void set_perspective(float fovy, float nc_scale);
float get_moon_light_factor();
void setup_basic_fog();
void check_zoom();
void reset_camera_pos();
void update_cpos();
void advance_camera(int dir);
bool open_file(FILE *&fp, char const *const fn, std::string const &file_type, char const *const mode="r");
void fire_weapon();
bool has_extension(std::string const &ext);

// function prototypes - visibility
void calc_mesh_shadows(unsigned l, point const &lpos, float const *const mh, unsigned char *smask, int xsize, int ysize,
					   float const *sh_in_x=NULL, float const *sh_in_y=NULL, float *sh_out_x=NULL, float *sh_out_y=NULL);
void calc_visibility(unsigned light_sources);
bool is_visible_to_light_cobj(point const &pos, int light, float radius, int cobj, int skip_dynamic, int *cobj_ix=NULL);
bool coll_pt_vis_test(point pos, point pos2, float dist, int &index, int cobj, int skip_dynamic, int test_alpha);
void calc_view_test_terms(float &tterm, float &sterm, bool is_zoomed);
void calc_viewing_cone();
void set_camera_pdu();
bool sphere_cobj_occluded(point const &viewer, point const &sc, float radius);
bool sphere_in_view(pos_dir_up const &pdu, point const &pos, float radius, int max_level, bool no_frustum_test=0);
int  get_light_pos(point &lpos, int light);
void update_sun_shadows();
void create_shadows();
void update_sun_and_moon();
int light_valid(unsigned light_sources, int l, point &lpos);

// function prototypes - shadows
void add_cobj_shadows(unsigned light_sources);
int  camera_shadow(point const &camera);
int  get_shape_shadow_bb(point const *points, int npoints, int l, int quality, point const &lpos,
	int &xmin, int &ymin, int &xmax, int &ymax, int &ret_val, unsigned char stype);
void get_sphere_points(point const &pos, float radius, point *pts, unsigned npts, vector3d const &dir);
int  enable_shadow_envelope(point const &pos, float radius, unsigned light_sources, int is_dynamic);
void disable_shadow_envelope(unsigned light_sources);
int  sphere_shadow2(point const &pos, float radius, unsigned light_sources, int is_dynamic, int quality);
int  sphere_shadow(point const &pos, float radius, unsigned light_sources, int is_dynamic, int quality);
int  cylinder_shadow(point p1, point p2, float radius1, float radius2, unsigned light_sources, int shadow_ends, int is_dynamic, int quality);
int  polygon_shadow(point const *points, vector3d const &norm, int npoints, float thick, unsigned light_sources,
					int is_dynamic, int quality, int is_cube, int tid=-1);
int  cube_shadow(cube_t const &cube, unsigned light_sources, int is_dynamic, int quality);
void reset_shadows(unsigned char type);

// function prototypes - mesh_intersect
bool sphere_visible_to_pt(point const &pt, point const &center, float radius);
bool is_visible_from_light(point const &pos, point const &lpos, int fast);
bool line_intersect_surface_cached(point const &v1, point const &v2, int &xpos, int &ypos, float &zval, int fast=0);
bool line_intersect_mesh(point const &v1, point const &v2, int &xpos, int &ypos, float &zval, int fast=0, bool cached=0);
bool line_intersect_mesh(point const &v1, point const &v2, int fast=0);
void gen_mesh_bsp_tree();
bool line_int_mesh_bsp(point const &v1, point const &v2);

// function prototypes - build
void create_object_groups();
void shift_all_objs(vector3d const &vd);
void process_platforms();
void process_groups();
void gen_scene(int generate_mesh, int gen_trees, int keep_sin_table, int update_zvals, int rgt_only, bool cobjs_re_add=0);
void init_models();
void free_models();

// function prototypes - display_world
point get_sun_pos();
point get_moon_pos();
colorRGBA get_bkg_color(point const &p1, vector3d const &v12);

// function prototypes - draw_world
void set_fill_mode();
int get_universe_ambient_light();
void set_gl_light_pos(int light, point const &pos, float w);
void set_colors_and_enable_light(int light, float const ambient[4], float const diffuse[4]);
int get_light();
void draw_solid_object_groups();
void draw_transparent_object_groups();
void set_shadowed_color(colorRGBA const &color, point const &pos, bool is_shadowed, bool precip=0, bool no_dynamic=0);
bool pt_is_shadowed(point const &pos, int light, int status, float radius, int cid, int fast);
void draw_select_groups(int solid);
void draw_group(obj_group &objg);
colorRGBA get_powerup_color(int powerup);
void draw_shadow_volume(point const &pos, point const &lpos, float radius, int &inverts);
int  draw_shadowed_objects(int light);
void check_drawing_flags(unsigned flags, int init_draw);
void set_specular(float specularity, float shininess);
colorRGBA get_glowing_obj_color(point const &pos, int time, int lifetime, float &stime, bool shrapnel_cscale, bool fade);
colorRGBA const &get_landmine_light_color(int time);
float get_landmine_sensor_height(float radius, int time);
colorRGBA get_plasma_color(float size);
void get_enabled_lights();
void set_dlights_booleans(shader_t &s, bool enable, int shader_type);
colorRGBA setup_smoke_shaders(shader_t &s, float min_alpha, int use_texgen, bool keep_alpha, bool indir_lighting,
	bool direct_lighting, bool dlights, bool smoke_en, bool has_lt_atten=0, bool use_smap=0);
void end_smoke_shaders(shader_t &s, colorRGBA const &orig_fog_color);
void setup_object_render_data();
void draw_coll_surfaces(bool draw_solid, bool draw_trans);
void draw_stars(float alpha);
void draw_sun();
void draw_moon();
void draw_earth();
void draw_stationary_earth(float radius);
void apply_red_sky(colorRGBA &color);
colorRGBA get_cloud_color();
float get_cloud_density(point const &pt, vector3d const &dir);
void draw_sky(int order);
void draw_stationary_sky(float radius, float density);
void compute_brightness();
void draw_water_plane(float zval, unsigned reflection_tid, int const *const hole_bounds);
void draw_bubbles();
void draw_smoke();
void draw_fires();
void draw_scorches();
void add_camera_filter(colorRGBA const &color, unsigned time, int tid, unsigned ix);
void draw_camera_filters(vector<camera_filter> &cfs);
void draw_projectile_effects();
void draw_env_other();
void mouse_draw_on_ground(int x, int y);
void draw_splash(float x, float y, float z, float size, colorRGBA color=WATER_C);
void draw_text(float x, float y, float z, char const *text, float tsize=1.0, bool bitmap_font=0);
void draw_framerate(float val);
void draw_compass_and_alt();
void exec_universe_text(std::string const &text);

// function prototypes - draw shapes
bool is_above_mesh(point const &pos);
bool check_face_containment(cube_t const &cube, int dim, int dir, int cobj);
float get_mesh_zmax(point const *const pts, unsigned npts);
void add_shadow_obj(point const &pos, float radius, int coll_id);
void add_coll_shadow_objs();
void get_occluders();

// function prototypes - draw primitives
void get_ortho_vectors(vector3d const &v12, vector3d *vab, int force_dim=-1);
vector_point_norm const &gen_cylinder_data(point const ce[2], float radius1, float radius2, unsigned ndiv, vector3d &v12,
										   float const *const perturb_map=NULL, float s_beg=0.0, float s_end=1.0, int force_dim=-1);
void draw_cylinder(float length, float radius1, float radius2, int ndiv, int nstacks, bool draw_ends, bool first_end_only=0, bool last_end_only=0);
void draw_cylinder_nstacks(float len, float r1, float r2, int ndiv, int nstacks, bool texture);
void draw_cylinder(point const &p1, float length, float radius1, float radius2, int ndiv, int nstacks, bool draw_ends);
void draw_circle_normal(float r_inner, float r_outer, int ndiv, int invert_normals);
void draw_circle_normal_at_z(float z, float r_inner, float r_outer, int ndiv, int invert_normals);
void draw_rotated_cylinder(point const &p1, point const &p2, float radius1, float radius2, int ndiv, int nstacks,
						   bool draw_ends, vector3d const &scale=zero_vector);
void draw_fast_cylinder(point const &p1, point const &p2, float radius1, float radius2, int ndiv, bool texture, int draw_sides_ends=0,
						float const *const perturb_map=NULL, bool const *const render_map=NULL, float const *const exp_map=NULL,
						point const *const pt_shift=NULL, float expand=0.0, float s_beg=0.0, float s_end=1.0);
void draw_cylindrical_section(point const &pos, float length, float r_inner, float r_outer, int ndiv, bool texture=0);
void draw_trunc_cone(point pos, vector3d v1, float length, float radius, float radius2, bool is_camera=0);
void draw_sphere_at(point const &pos, float radius, int ndiv);
void draw_sphere_at_tc(point const &pos, float radius, int ndiv, bool texture, bool cull);
void draw_subdiv_sphere(point const &pos, float radius, int ndiv, point const &vfrom, float const *perturb_map,
						int texture, bool disable_bfc, bool const *const render_map=NULL, float const *const exp_map=NULL,
						point const *const pt_shift=NULL, float expand=0.0, float s_beg=0.0, float s_end=1.0, float t_beg=0.0, float t_end=1.0);
void draw_subdiv_sphere(point const &pos, float radius, int ndiv, int texture, bool disable_bfc);
void draw_subdiv_sphere_section(point const &pos, float radius, int ndiv, int texture,
								float s_beg, float s_end, float t_beg, float t_end);
void rotate_sphere_tex_to_dir(vector3d const &dir);
void draw_cube_map_sphere(point const &pos, float radius, int ndiv, bool texture, bool disable_bfc=0);
void draw_torus(float ri, float ro, unsigned ndivi, unsigned ndivo, bool do_tex);
void rotate_towards_camera(point const &pos);
void draw_textured_square(float size, float z, int tid);
void draw_textured_square_alpha_test(float size, float z, int tid);
void draw_flare(point const &pos, point const &xlate, float xsize, float ysize);
void enable_flares(colorRGBA const &color, bool zoomed=0);
void disable_flares();
void draw_textured_quad(float xsize, float ysize, float z, int tid);
void draw_tquad(float xsize, float ysize, float z, bool texture);
void draw_one_tquad(float x1, float y1, float x2, float y2, float z, bool texture, float tx1=0.0, float ty1=0.0, float tx2=1.0, float ty2=1.0);
void draw_one_mult_tex_quad(float x1, float y1, float x2, float y2, float z, float tx1=0.0, float ty1=0.0, float tx2=1.0, float ty2=1.0);
void draw_billboard(point const &pos, point const &viewer, vector3d const &up_dir,
					float xsize, float ysize, float tx1=0.0, float ty1=0.0, float tx2=1.0, float ty2=1.0);
void draw_line_tquad(point const &p1, point const &p2, float w1, float w2, colorRGBA const &color1, colorRGBA const &color2, bool globalize);
void begin_line_tquad_draw();
void end_line_tquad_draw();
void draw_animated_billboard(point const &pos, float size, float timescale);
int  draw_simple_cube(cube_t const &c, bool texture, int in_cur_prim=PRIM_DISABLED, bool no_normals=0, int eflags=0,
	float texture_scale=1.0, vector3d const *const view_dir=NULL);
void draw_cube(point const &pos, float sx, float sy, float sz, bool texture, unsigned ndiv, bool scale_ndiv=0,
			   float texture_scale=1.0, bool proportional_texture=0, vector3d const *const view_dir=NULL);
int draw_cylin_quad_proj(cylinder_3dw const &cylin, vector3d const &view_dir, int in_cur_prim=PRIM_DISABLED, bool no_normals=0);
int  draw_simple_polygon(point const *const points, int npoints, vector3d const &norm, int in_cur_prim=PRIM_DISABLED, bool no_normals=0);
int  draw_simple_extruded_polygon(float thick, point const *const points, int npoints, int in_cur_prim=PRIM_DISABLED, bool no_normals=0);
void gen_quad_tex_coords(float *tdata, unsigned num, unsigned stride);
void gen_quad_tri_tex_coords(float *tdata, unsigned num, unsigned stride);
void draw_quads_from_pts(vector<vert_norm> const &points, unsigned draw_num=0);
void free_dlists();
void setup_dlists();
void draw_cylin_fast(float r1, float r2, float l, int ndiv, bool texture, bool restore_matrix, bool r_approx=0);
void draw_sphere_dlist_raw(int ndiv, bool textured, bool half=0);
void draw_sphere_dlist(point const &pos, float radius, int ndiv, bool textured, bool half=0, bool bfc=0);
void draw_sphere_dlist_back_to_front(point const &pos, float radius, int ndiv, bool textured, bool half=0);
void draw_rotated_cylinder_dlist(point const &p1, point const &p2, float r, int ndiv, vector3d const &scale=zero_vector);

// function prototypes - draw mesh
float integrate_water_dist(point const &targ_pos, point const &src_pos, float const water_z);
void water_color_atten_pt(float *c, int x, int y, point const &pos, point const &p1, point const &p2);
float get_cloud_shadow_atten(int x, int y);
colorRGBA setup_mesh_lighting();
void run_post_mesh_draw();
void set_landscape_texgen(float tex_scale, int xoffset, int yoffset, int xsize, int ysize, bool use_detail_tex=1);
void display_mesh();
float display_mesh3(int const *const hole_bounds, float wpz);
void draw_water_sides(int check_zvals);
float get_inf_terrain_fog_dist();

// function prototypes - tiled mesh
int get_tile_radius();
float draw_tiled_terrain(bool add_hole, float wpz);
void clear_tiled_terrain();
void reset_tiled_terrain_state();

// function prototypes - map_view
void draw_overhead_map();

// function prototypes - gen_obj
void gen_stars(float alpha, int half_sphere);
void gen_star(star &star1, int half_sphere);
void rand_xy_point(float zval, point &pt, unsigned flags);
void gen_object_pos(point &position, unsigned flags);
void gen_bubble(point const &pos, float r=0.0, colorRGBA const &c=WATER_C);
void gen_line_of_bubbles(point const &p1, point const &p2, float r=0.0, colorRGBA const &c=WATER_C);
bool gen_arb_smoke(point const &pos, colorRGBA const &bc, vector3d const &iv,
				   float r, float den, float dark, float dam, int src, int dt, bool as);
void gen_smoke(point const &pos);
bool gen_fire(point const &pos, float size, int source);
void gen_scorch_mark(point const &pos, float radius, vector3d const &orient, int cid=-1, float init_alpha=1.0, float rgb_val=0.0);
void gen_particles(point const &pos, unsigned num, float lt_scale=1.0, bool fade=0);
int gen_fragment(point const &pos, vector3d const &velocity, float size_mult, float time_mult,
	colorRGBA const &color, int tid, float tscale, int source, bool tri_fragment);
void gen_leaf_at(point const *const points, vector3d const &normal, int type, colorRGB const &color);
void gen_cloud_volumes();
float rgauss();
void gen_gauss_rand_arr();
int  rand2();
double rand2d();

// function prototypes - mesh_gen
bool bmp_to_chars(char *fname, char **&data);
bool verify_bmp_header(FILE *&fp, bool grayscale);
bool open_image_file(char *filename, FILE *&fp, bool is_bmp, bool is_grayscale);
void gen_mesh_sine_table(float **matrix, float **jterms, int x_offset, int y_offset, int xsize, int ysize, int make_island);
void gen_mesh(int surface_type, int make_island, int keep_sin_table, int update_zvals);
float get_glaciated_zval(float zval);
float calc_glaciated_rel_value(float value);
void init_terrain_mesh();
void build_xy_mesh_arrays(float *xv, float *yv, int nx, int ny);
float fast_eval_from_index(int x, int y, bool use_cache=1, bool glaciate=1);
float eval_one_surface_point(float xval, float yval);
void reset_offsets();
bool is_mesh_disabled(int xpos, int ypos);
float get_water_z_height();
void update_mesh(float dms, bool do_regen_trees);
bool is_under_mesh(point const &p);
bool read_mesh(const char *filename, float zmm=0.0);
bool write_mesh(const char *filename);
bool load_state(const char *filename);
bool save_state(const char *filename);

// function prototypes - physics
float get_max_t(int obj_type);
void init_objects();
void set_coll_rmax(float rmax);
void change_timestep(float mult_factor);
vector3d get_local_wind(point const &pt);
void reanimate_objects();
void accumulate_object(point const &pos, int type);
void shift_other_objs(vector3d const &vd);
void advance_physics_objects();
void reset_other_objects_status();
void auto_advance_time();

// function prototypes - ai
void advance_smiley(dwobject &obj, int smiley_id);
void shift_player_state(vector3d const &vd, int smiley_id);
void smiley_action(int smiley_id);

// function prototypes - matrix
void set_scene_constants();
void alloc_matrices();
void delete_matrices();
void compute_matrices();
void update_matrix_element(int i, int j);
void update_mesh_height(int xpos, int ypos, int rad, float scale, float offset, int mode);
vector3d get_matrix_surf_norm(float **matrix, char **enabled, int xsize, int ysize, int i, int j);
void calc_matrix_normal_at(float **matrix, vector3d **vn, vector3d **sn, char **enabled, int xsize, int ysize, int i, int j);
void calc_matrix_normals(float **matrix, vector3d **vn, vector3d **sn, char **enabled, int xsize, int ysize);
void get_matrix_point(int xpos, int ypos, point &pt);
int  is_in_ice(int xpos, int ypos);
float interpolate_mesh_zval(float xval, float yval, float rad, int use_real_equation, int ignore_ice);
float int_mesh_zval_pt_off(point const &pos, int use_real_equation, int ignore_ice);
void calc_motion_direction();
float lowest_mesh_point(point const &pt, float radius);
float highest_mesh_point(point const &pt, float radius);

// function prototypes - collision detection
int  add_coll_cube(cube_t &cube, cobj_params const &cparams, int platform_id=-1, int dhcm=0);
int  add_coll_cylinder(float x1, float y1, float z1, float x2, float y2, float z2,
					   float radius, float radius2, cobj_params const &cparams, int platform_id=-1, int dhcm=0);
int  add_coll_sphere(point const &pt, float radius, cobj_params const &cparams, int platform_id=-1, int dhcm=0);
int  add_coll_polygon(const point *points, int npoints, cobj_params const &cparams,
	float thickness, point const &xlate=all_zeros, int platform_id=-1, int dhcm=0);
int  remove_coll_object(int index, bool reset_draw=1);
int  remove_reset_coll_obj(int &index);
void purge_coll_freed(bool force);
void remove_all_coll_obj();
void cobj_optimize();
int  collision_detect_large_sphere(point &pos, float radius, unsigned flags);
int  check_legal_move(int x_new, int y_new, float zval, float radius, int &cindex);

// function prototypes - coll_cell_search
void build_cobj_tree( bool dynamic=0, bool verbose=1);
void update_cobj_tree(bool dynamic=0, bool verbose=1);
void build_moving_cobj_tree();
bool check_coll_line_exact_tree(point const &p1, point const &p2, point &cpos, vector3d &cnorm,
	int &cindex, int ignore_cobj, bool dynamic=0, int test_alpha=0);
bool check_coll_line_tree(point const &p1, point const &p2, int &cindex, int ignore_cobj, bool dynamic=0, int test_alpha=0);
bool cobj_contained_tree(point const &p1, point const &p2, point const &viewer, point const *const pts, unsigned npts,
	int ignore_cobj, int &cobj);
void get_coll_line_cobjs_tree(point const &pos1, point const &pos2, int ignore_cobj, vector<int> &cobjs);
void get_intersecting_cobjs_tree(cube_t const &cube, vector<unsigned> &cobjs, int ignore_cobj, float toler,
	bool dynamic, bool check_ccounter, int id_for_cobj_int);
bool is_contained(point const &pos, point const *const pts, unsigned npts, float const d[3][2]);
bool check_coll_line(point pos1, point pos2, int &cindex, int c_obj, int skip_dynamic, int test_alpha, bool no_tree=0);
bool check_coll_line_exact(point pos1, point pos2, point &cpos, vector3d &coll_norm, int &cindex, float splash_val=0.0,
						   int ignore_cobj=-1, bool fast=0, bool test_alpha=0, bool skip_dynamic=0, bool no_tree=0);
bool cobj_contained(point pos1, point center, const point *pts, unsigned npts, int cobj);
bool get_coll_line_cobjs(point pos1, point pos2, int cobj, vector<int> &cobjs);
bool coll_pt_vis_test_large(point pos1, point pos2, vector<int> &cobjs, int cobj, float radius, int skip_dynamic);
bool is_occluded(vector<int> const &occluders, point const *const pts0, int npts, point const &camera);
void add_camera_cobj(point const &pos);
void force_onto_surface_mesh(point &pos);
int  set_true_obj_height(point &pos, point const &lpos, float step_height, float &zvel, int type, int id,
	bool flight, bool on_snow, bool skip_dynamic=0, bool test_only=0);

// function prototypes - math3d
float fix_angle(float angle);
vector3d get_poly_norm(point const *points);
void get_face_normal(shape3d &shape, int face_id);
void calc_reflection_angle(vector3d const &v_inc, vector3d &v_ref, vector3d const &norm);
bool calc_refraction_angle(vector3d const &v_inc, vector3d &v_ref, vector3d const &norm, float n1, float n2);
float get_fresnel_reflection(vector3d const &v_inc, vector3d const &norm, float n1, float n2);
float get_reflected_weight(float fresnel_ref, float alpha);
float get_coll_energy(vector3d const &v1, vector3d const &v2, float mass);
float triangle_area(point const *const points);
bool planar_contour_intersect(const point *points, unsigned npoints, point const &pos, vector3d const &norm);
bool point_in_polygon_2d(float xval, float yval, const point *points, int npts, int dx, int dy);
bool get_poly_zminmax(point const *const pts, unsigned npts, vector3d const &norm, float dval,
					  cube_t const &cube, float &z1, float &z2);
void grow_poly_about_center(point *pts, unsigned npts, float scale);
bool get_poly_zvals(vector<vector<point> > const &pts, float xv, float yv, float &z1, float &z2);
void gen_poly_planes(point const *const points, unsigned npoints, vector3d const &norm, float thick, vector<point> pts[2]);
vector<vector<point> > const &thick_poly_to_sides(point const *const points, unsigned npoints, vector3d const &norm, float thick);
bool line_int_plane(point const &p1, point const &p2, point const &pp0, vector3d const &norm,
					point &p_int, float &t, bool ignore_t);
bool thick_poly_intersect(vector3d const &v1, point const &p1, vector3d const &norm,
						  vector<point> const pts[2], bool test_side, unsigned npoints);
vector3d get_poly_dir_norm(vector3d const &norm, point const &p1, vector3d const &v1, float t);
bool sphere_intersect_poly_sides(vector<vector<point> > const &pts, point const &center,
								 float radius, float &dist, vector3d &norm, bool strict);
bool pt_line_seg_dist_less_than(point const &P, point const &L1, point const &L2, float dist);
bool sphere_poly_intersect(const point *points, unsigned npoints, point const &pos, vector3d const &norm, float rdist, float radius);
bool sphere_ext_poly_int_base(point const &pt, vector3d const &norm, point const &pos, float radius,
							  float thickness, float &thick, float &rdist);
bool sphere_ext_poly_intersect(point const *const points, unsigned npoints, vector3d const &norm,
							   point const &pos, float radius, float thickness, float t_adj);
void get_sphere_mov_sphere_int_pt(point const &p1, point const &p2, vector3d const &v, float rsum, point &cpos);
bool sphere_test_comp(point const &p2, point const &p1, vector3d const &v1, float r2sq, float &t);
bool line_sphere_intersect(point const &p1, point const &p2, point const &c, float r);
bool circle_test_comp(point const &p2, point const &p1, vector3d const &v1, vector3d norm, float r2sq, float &t);
void dir_to_sphere_s_t(vector3d const &dir, vector3d const &sdir, double &s, double &t);
bool line_sphere_intersect_s_t(point const &p1, point const &p2, point const &sc, float radius,
							   vector3d const &sdir, double &s, double &t);
bool line_sphere_int(vector3d const &v1, point const &p1, point const &center, float radius, point &lsint, bool test_neg_t);
bool line_intersect_sphere(point const &p1, vector3d const &v12, point const &sc, float radius, float &rad, float &dist, float &t);
void get_sphere_border_pts(point *qp, point const &pos, point const &viewed_from, float radius, unsigned num_pts);
bool line_torus_intersect(point const &p1, point const &p2, point const &tc, float ri, float ro, float &t);
bool sphere_torus_intersect(point const &sc, float sr, point const &tc, float ri, float ro, point &p_int, vector3d &norm, bool calc_int);
bool circle_rect_intersect(point const &pos, float radius, cube_t const &cube);
bool sphere_cube_intersect(point const &pos, float radius, cube_t const &cube);
bool sphere_cube_intersect(point const &pos, float radius, cube_t const &cube, point const &p_last,
						   point &p_int, vector3d &norm, unsigned &cdir, bool check_int, bool skip_z=0);
bool do_line_clip(point &v1, point &v2, float const d[3][2]);
bool get_line_clip(point const &v1, point const &v2, float const d[3][2], float &tmin, float &tmax);
bool get_line_clip2(point const &v1, vector3d const &dinv, float const d[3][2]);
bool check_line_clip_expand(point const &v1, point const &v2, float const d[3][2], float expand);
float line_line_dist(point const &p1a, point const &p1b, point const &p2a, point const &p2b);
float get_cylinder_params(point const &cp1, point const &cp2, point const &pos, vector3d &v1, vector3d &v2);
int  line_intersect_trunc_cone(point const &p1, point const &p2, point const &cp1, point const &cp2,
							   float r1, float r2, bool check_ends, float &t);
bool line_intersect_cylinder(point const &p1, point const &p2, cylinder_3dw const &c, bool check_ends);
int  line_int_thick_cylinder(point const &p1, point const &p2, point const &cp1, point const &cp2,
							 float ri1, float ri2, float ro1, float ro2, bool check_ends, float &t);
bool sphere_int_cylinder_pretest(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2,
								 bool check_ends, vector3d &v1, vector3d &v2, float &t, float &rad);
bool sphere_intersect_cylinder_ipt(point const &sc, float sr, point const &cp1, point const &cp2, float r1, float r2,
							   bool check_ends, point &p_int, vector3d &norm, bool calc_int);
void cylinder_quad_projection(point *pts, cylinder_3dw const &c, vector3d const &v1, vector3d &v2, int &npts);
template<typename T> pointT<T> get_center_arb(pointT<T> const *const pts, int npts);
unsigned get_cube_corners(float const d[3][2], point corners[8], point const &viewed_from=all_zeros, bool all_corners=1);
void get_closest_cube_norm(float const d[3][2], point const &p, vector3d &norm);
void cylinder_bounding_sphere(point const *const pts, float r1, float r2, point &center, float &radius);
void polygon_bounding_sphere(const point *pts, int npts, float thick, point &center, float &radius);
void add_rotated_quad_pts(vector<vert_norm> &points, float theta, float rd, float z, point const &pos, vector3d const &scale);
void vproj_plane(vector3d const &vin, vector3d const &n, vector3d &vout);
template<typename T> void rotate_vector3d(pointT<T> vin, pointT<T> const &vrot, double angle, pointT<T> &vout);
template<typename T> void rotate_vector3d_multi(pointT<T> const &vrot, double angle, pointT<T> *vout, unsigned nv);
void rotate_vector3d_x2(point const &vrot, double angle, point &vout1, point &vout2);
float angle_of_projected_vectors(vector3d const &v1, vector3d const &v2, vector3d n);
vector3d rtp_to_xyz(float radius, double theta, double phi);
vector3d gen_rand_vector_uniform(float mag);
vector3d gen_rand_vector(float mag, float zscale=1.0, float phi_term=PI);
vector3d gen_rand_vector2(float mag, float zscale=1.0, float phi_term=PI);
vector3d lead_target(point const &ps, point const &pt, vector3d const &vs, vector3d const &vt, float vweap);
vector3d get_firing_dir(vector3d const &src, vector3d const &dest, float fvel, float gravity_scale);

// function prototypes - water
bool get_water_enabled(int x, int y);
bool has_water(int x, int y);
bool mesh_is_underwater(int x, int y);
void select_water_ice_texture(colorRGBA &color, float *use_this_temp=NULL);
void draw_water();
void add_splash(int xpos, int ypos, float energy, float radius);
bool add_water_section(float x1, float y1, float x2, float y2, float zval, float wvol);
void float_downstream(point &pos, float radius);
void calc_watershed();
bool is_underwater(point const &pos, int check_bottom=0, float *depth=NULL);
void select_liquid_color(colorRGBA &color, int xpos, int ypos);
void select_liquid_color(colorRGBA &color, point const &pos);
void add_water_spring(point const &pos, vector3d const &vel, float rate, float diff, int calc_z, int gen_vel);
void shift_water_springs(vector3d const &vd);
void update_water_zvals(int x1, int y1, int x2, int y2);

// function prototypes - lightning
void compute_volume_matrix();

// function prototypes - textures
void load_textures();
int get_texture_by_name(std::string const &name);
bool select_texture(int id, bool enable=1);
float get_tex_ar(int id);
void bind_1d_texture(unsigned tid);
void bind_2d_texture(unsigned tid);
void setup_texture(unsigned &tid, int type, bool mipmap, bool wrap_s, bool wrap_t, bool mirror_s=0, bool mirror_t=0, bool nearest=0);
void free_textures();
void reset_textures();
void free_texture(unsigned &tid);
void setup_landscape_tex_colors(colorRGBA const &c1, colorRGBA const &c2);
colorRGBA texture_color(int tid);
unsigned get_texture_size(int tid, bool dim);
void get_lum_alpha(colorRGBA const &color, int tid, float &luminance, float &alpha);
colorRGBA get_landscape_texture_color(int xpos, int ypos);
void update_lttex_ix(int &ix);
void get_tids(float relh, int NTEXm1, float const *const h_tex, int &k1, int &k2, float &t);
void create_landscape_texture();
float add_crater_to_landscape_texture(float xval, float yval, float radius);
void add_cutout_to_landscape_texture(int x1, int y1, int x2, int y2);
void add_color_to_landscape_texture(colorRGBA const &color, float xval, float yval, float radius, int check_unique);
void add_snow_to_landscape_texture(point const &pos, float acc);
void update_landscape_texture();
void gen_tex_height_tables();
void set_texgen_vec4(float const v[4], bool s_or_t, bool enable_and_set_mode, shader_t *shader=NULL);
void setup_texgen_full(float sx, float sy, float sz, float sw, float tx, float ty, float tz, float tw, int mode=GL_EYE_LINEAR, shader_t *shader=NULL);
void setup_texgen(float xscale, float yscale, float tx, float ty, float z_off=0.0, int mode=GL_EYE_LINEAR);
void disable_texgen();
void disable_textures_texgen();
void setup_polygon_texgen(vector3d const &norm, float const scale[2], float const xlate[2], vector3d const &offset, bool swap_txy=0, shader_t *shader=NULL);
void get_tex_coord(vector3d const &dir, vector3d const &sdir, unsigned txsize, unsigned tysize, int &tx, int &ty, bool invert);
float get_texture_component(unsigned tid, float xval, float yval, int comp);
bool is_billboard_texture_transparent(point const *const points, point const &pos, int tid);
void set_texture_specular(bool val);

// function prototypes - sun flares
void DoFlares(point const &from, point const &at, point const &light, float near_clip, float size);
void load_flare_textures();
void free_flare_textures();

// function prototypes - ocean
void draw_ocean();
double gamma (double x);
double MitsuyasuDistribution (double f, double theta);
double EPM(double f);
double EnergyDistribution(double f, double theta);
double Amplitude(double f, double theta, double k);
void AnimateWater();
void InitWater();
void calc_ocean_normals();

// function prototypes - gameplay/ai
void camera_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void smiley_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void landmine_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void health_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void shield_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void powerup_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void weapon_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void ammo_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void pack_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void rock_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void sball_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void dodgeball_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
void skull_collision(int index, int obj_index, vector3d const &velocity, point const &position, float energy, int type);
int get_smiley_hit(vector3d &hdir, int index);
void blast_radius(point const &pos, int type, int obj_index, int shooter, int chain_level);
void create_explosion(point const &pos, int shooter, int chain_level, float damage, float size, int type, bool cview);
void do_area_effect_damage(point &pos, float effect_radius, float damage, int index, int source, int type);
void switch_player_weapon(int val);
void draw_lasers();
void show_blood_on_camera();
void update_weapon_cobjs();
int select_dodgeball_texture(int shooter);
void draw_weapon_simple(point const &pos, vector3d const &dir, float radius, int cid, int wid, float scale);
void draw_weapon_in_hand(int shooter);
bool weap_has_transparent(int shooter);
void draw_scheduled_weapons();
void add_weapon_lights(int shooter);
void show_crosshair(int do_zoom);
void show_user_stats();
void show_other_messages();
void print_text_onscreen(std::string const &text, colorRGBA const &color, float size, int time, int priority=0, bool bitmap=0);
void print_weapon(int weapon_id);
bool check_underwater(int who, float &depth);
void player_fall(int id);
void update_camera_velocity(vector3d const &v);
void init_game_state();
void gamemode_rand_appear();
bool has_invisibility(int id);
void init_smileys();
void init_game_mode();
void update_game_frame();
void change_game_mode();
void free_dodgeballs(bool camera, bool smileys);
int gen_smiley_or_player_pos(point &pos, int index);
colorRGBA get_smiley_team_color(int smiley_id);
void select_smiley_texture(int smiley_id);
void free_smiley_textures();
int get_ammo_or_obj(int wid);
int wid_need_weapon(int wid);

// function prototypes - explosion
void update_blasts();
void draw_blasts();

// function prototypes - scenery
void gen_scenery();
void draw_scenery(bool draw_opaque, bool draw_transparent, bool shadow_only=0);
void update_scenery_zvals(int x1, int y1, int x2, int y2);
void free_scenery();
void add_scenery_cobjs();
void shift_scenery(vector3d const &vd);
void add_plant(point const &pos, float height, float radius, int type, int calc_z);

// function prototypes - grass
void setup_wind_for_shader(shader_t &s);
bool no_grass();
void gen_grass(bool full_regen);
void update_grass_vbos();
void draw_grass();
void modify_grass_at(point const &pos, float radius, bool crush=0, bool burn=0, bool cut=0, bool update_mh=0, bool check_uw=0);
bool place_obj_on_grass(point &pos, float radius);
float get_grass_density(point const &pos);

// function prototypes - draw mech
void draw_hmv();
void build_hmv_shape();
void delete_hmv_shape();
void add_shape_coll_objs();
void shift_hmv(vector3d const &vd);

// function prototypes - tree + sm_tree (see also tree_3dw.h)
void mult_leaf_points_by(float val);
int get_tree_type_from_height(float zpos);
void set_leaf_shader(shader_t &s, float min_alpha, bool use_wind=0);

// function prototypes - csg
void get_cube_points(const float d[3][2], point pts[8]);
void check_cubes(vector<coll_obj> &cobjs);
void merge_cubes(vector<coll_obj> &cobjs);
void process_negative_shapes(vector<coll_obj> &cobjs);
void remove_overlapping_cubes(vector<coll_obj> &cobjs);
void subdiv_cubes(vector<coll_obj> &cobjs);
void sort_cobjs_for_rendering(vector<coll_obj> &cobjs);

// function prototypes - ship
upos_point_type const &get_player_pos();
vector3d const &get_player_dir();
vector3d const &get_player_up();
void set_player_pos(point const &pos_);
void set_player_dir(vector3d const &dir_);
void set_player_up(vector3d const &upv_);
void init_universe_display();
void set_univ_pdu();
void setup_current_system();
void apply_univ_physics();
void draw_universe(bool static_only=0, bool skip_closest=0, bool no_distant=0);
void draw_universe_stats();

// function prototypes - lightmap
void update_flow_for_voxels(cube_t const &cube);
void shift_light_sources(vector3d const &vd);
void shift_lightmap(vector3d const &vd);
void regen_lightmap();
void clear_lightmap();
void build_lightmap(bool verbose);
void add_smoke(point const &pos, float val);
void distribute_smoke();
float get_smoke_at_pos(point const &pos);
void add_line_light(point const &p1, point const &p2, colorRGBA const &color, float size, float intensity=1.0);
void add_dynamic_light(float sz, point const &p, colorRGBA const &c=WHITE, vector3d const &d=plus_z, float bw=1.0);
colorRGBA gen_fire_color(float &cval, float &inten);
void clear_dynamic_lights();
void add_dynamic_lights();
void get_xyz(point const &p, int v[3]);
void get_xyz_v2(point &p, int const v[3]);
bool upload_smoke_3d_texture();
void upload_dlights_textures();
void setup_dlight_textures(shader_t &s);
bool is_shadowed_lightmap(point const &p);
bool is_in_darkness(point const &pos, float radius, int cobj);
bool get_dynamic_light(int x, int y, int z, point const &p, float lightscale, float *ls, vector3d const *const norm, float const *const spec);
void get_sd_light(int x, int y, int z, point const &p, bool no_dynamic, float lightscale, float *ls, vector3d const *const norm, float const *const spec);
float get_indir_light(colorRGBA &a, colorRGBA cscale, point const &p, bool no_dynamic, bool shadowed, vector3d const *const norm, float const *const spec);
unsigned enable_dynamic_lights(point const center=all_zeros, float radius=0.0);
void disable_dynamic_lights(unsigned num_dlights);

// function prototypes - tessellate
void split_polygon_to_tris(vector<triangle> &triangles_out, vector<point> const &poly_pts);
bool split_polygon_to_cobjs(coll_obj const &cobj, vector<coll_obj> &split_polygons, vector<point> const &poly_pt, bool split_quads);

// function prototypes - shaders
char const *append_array_ix(std::string &s, unsigned i);
bool setup_shaders();
void clear_shaders();

// function prototypes - snow
bool snow_enabled();
void gen_snow_coverage();
void reset_snow_vbos();
void draw_snow();
bool get_snow_height(point const &p, float radius, float &zval, vector3d &norm, bool crush_snow);

// function prototypes - waypoints
void create_waypoints(vector<user_waypt_t> const &user_waypoints);
void shift_waypoints(vector3d const &vd);
void draw_waypoints();

// function prototypes - destroy_cobj
void destroy_coll_objs(point const &pos, float damage, int shooter, bool big);
void check_falling_cobjs();

// function prototypes - shadow_map
bool shadow_map_enabled();
unsigned get_shadow_map_tu_id(int light);
unsigned get_shadow_map_tid(int light);
int get_smap_ndiv(float radius);
void set_smap_shader_for_light(shader_t &s, int light, float z_bias);
void set_smap_shader_for_all_lights(shader_t &s, float z_bias=0.0005);
void draw_scene_bounds_and_light_frustum(point const &lpos);
void create_shadow_map();
void free_shadow_map_textures();

// function prototypes - screenshot (these are C functions)
#ifdef ENABLE_JPEG
extern "C" int screenshot(unsigned window_width, unsigned window_height, char *file_path);
extern "C" int write_jpeg(unsigned window_width, unsigned window_height, char *file_path);
#else // why compile as C when libjpeg isn't needed?
int screenshot(unsigned window_width, unsigned window_height, char *file_path);
int write_jpeg(unsigned window_width, unsigned window_height, char *file_path);
#endif


#include "inlines.h"


#endif // _3DWORLD_H_

