// 3D World - OpenGL CS184 Computer Graphics Project
// by Frank Gennari
// 3/10/02

#pragma once

// timer_t is used in types.h on linux, and also in 3DWorld, so we have to typedef it as something else for these includes
#define timer_t stdlib_timer_t

#include "globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gl_includes.h"
#include "rand_gen.h"

// STL include (others come from rand_gen.h)
#include <deque>
#include <algorithm>
#include <set>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>

#undef timer_t

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
using std::string;

#ifndef PI
#define PI 3.141592654f
#endif

int      const CAMERA_ID        = -1;
int      const NO_SOURCE        = -2;
unsigned const INIT_CCELL_SIZE  = 16;
float    const LARGE_OBJ_RAD    = 0.01;
float    const TOLERANCE        = 1.0E-12;
float    const ABSOLUTE_ZERO    = -273; // in degrees C
float    const MAX_SPLASH_DEPTH = 0.1;
float    const WATER_INDEX_REFRACT = 1.333;
float    const ICE_INDEX_REFRACT   = 1.309;
float    const WATER_COL_ATTEN  = 0.6;

float    const MESH_BOT_QUAD_DZ = 0.050; // measured from zbottom
float    const MIN_WATER_DZ     = 0.045;
float    const MESH_LOWEST_DZ   = 0.040;

unsigned const TICKS_PER_SECOND = 40;
float    const SMALL_NUMBER     = 0.001;
float    const MIN_POLY_THICK   = 0.001;
unsigned const MAX_CHARS        = 256;

float    const LIGHT_ROT_AMT    = 0.05;
float    const WIND_ADJUST      = 0.2;

float    const DEF_TIMESTEP     = 0.007;
int      const TIMESCALE        = 1; // (integer) larger makes surface movement faster
int      const MAX_I_TIMESCALE  = 8; // (integer) larger makes surface movement slower (max inverse of TIMESCALE)
float    const GRAVITY          = 300.0;

float    const CLOUD_CEILING0   = 1.5;
int      const LITNING_TIME     = 50; // in ticks
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
float    const RAIN_TAIL_MIN_V  = -1.0;

float    const WATER_ALPHA      = 0.75;
float    const ICE_ALPHA        = 0.88;
float    const SPLASH_ALPHA     = 1.0;
float    const SNOW_ALPHA       = 1.0;

int      const N_RAND_GAUSS     = 10;
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
float    const DEF_NEAR_CLIP    = 0.01;
float    const DEF_FAR_CLIP     = 100.0;
float    const FAR_DISTANCE     = 100.0;
float    const PERSP_ANGLE      = 60.0;
float    const ZOOM_FACTOR      = 5.0;
float    const UNIV_NCLIP_SCALE = 0.02;

float    const MIN_SHADOW_ALPHA = 0.5;
float    const GET_OCC_EXPAND   = 0.02;
float    const DEF_Z_BIAS       = 0.0005;

unsigned const CELL_Z_DIVS      = 8;
unsigned const SMALL_NDIV       = 8;
int const FAST_VISIBILITY_CALC  = 1;

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

unsigned const MAX_SHADER_LIGHTS = 8;

unsigned const quad_to_tris_ixs [6] = {0,1,2, 0,2,3};
unsigned const cube_dim_table[2][3] = {{1, 2, 0}, {2, 0, 1}};


#define CLIP_TO_01(x)  max( 0.0f, min(1.0f, (x)))
#define CLIP_TO_pm1(x) max(-1.0f, min(1.0f, (x)))

#define BITSHIFT_CEIL(num, bs) (((num-1) >> bs) + 1)

#define UNROLL_2X(expr) {{unsigned const i_(0); expr} {unsigned const i_(1); expr}}
#define UNROLL_3X(expr) {UNROLL_2X(expr) {unsigned const i_(2); expr}}
#define UNROLL_4X(expr) {UNROLL_3X(expr) {unsigned const i_(3); expr}}

template<typename T> inline void min_eq(T &A, T const B) {A = min(A, B);}
template<typename T> inline void max_eq(T &A, T const B) {A = max(A, B);}

enum {CAM_FILT_DAMAGE=0, CAM_FILT_FOG, CAM_FILT_BURN, CAMERA_FILT_BKG, CAM_FILT_UWATER, CAM_FILT_TELEPORT, CAM_FILT_FROZEN, CAM_FILT_END};
enum {FG_PROJECTION=0, FG_MODELVIEW};
enum {GAME_MODE_NONE=0, GAME_MODE_FPS, GAME_MODE_DODGEBALL};
unsigned const GAME_MODE_BUILDINGS = 2; // same as GAME_MODE_DODGEBALL


template<typename T> struct point2d { // size = 8

	T x, y;

	point2d() : x(0.0), y(0.0) {}
	point2d(T x_, T y_) : x(x_), y(y_) {}
	point2d(point2d const &a, point2d const &b) : x(a.x - b.x), y(a.y - b.y) {}
	void assign(T x_, T y_) {x = x_; y = y_;}
	bool operator==(point2d const &p) const {return (x == p.x && y == p.y);}
	bool operator!=(point2d const &p) const {return (x != p.x || y != p.y);}
	T mag_sq() const {return (x*x + y*y);}
	T mag()    const {return sqrt(mag_sq());}
	T cp_mag(point2d const &p) const {return (x*p.y - y*p.x);}
	T get_max_val() const {return max(x, y);}
	T get_min_val() const {return min(x, y);}
	void operator+=(point2d const &p) {x += p.x; y += p.y;}
	void operator-=(point2d const &p) {x -= p.x; y -= p.y;}
	void operator*=(point2d const &p) {x *= p.x; y *= p.y;} // component multiply
	void operator+=(T const &v)      {x += v; y += v;}
	void operator-=(T const &v)      {x -= v; y -= v;}
	void operator*=(T m)             {x *= m; y *= m;}
	point2d operator+(point2d const &p) const {return point2d((x+p.x), (y+p.y));}
	point2d operator-(point2d const &p) const {return point2d((x-p.x), (y-p.y));}
	point2d operator+(T const &v)       const {return point2d((x+v), (y+v));}
	point2d operator-(T const &v)       const {return point2d((x-v), (y-v));}
	point2d operator*(T      const val) const {return point2d(x*val, y*val);}
	point2d operator*(point2d const &p) const {return point2d(x*p.x, y*p.y);} // component multiply
	point2d operator-()                 const {return point2d(-x, -y);}

	const T &operator[](unsigned i) const {
		switch(i) {
			case 0: return x;
			case 1: return y;
			default: assert(0);
		}
		return x; // never gets here
	}
	T &operator[](unsigned i) {
		switch(i) {
			case 0: return x;
			case 1: return y;
			default: assert(0);
		}
		return x; // never gets here
	}
	void negate() {x = -x; y = -y;}
	void normalize() {
		T const d(mag());
		if (d >= TOLERANCE) {x /= d; y /= d;}
	}
	point2d get_norm() const {
		T const vmag(mag());
		return ((vmag < TOLERANCE) ? *this : point2d(x/vmag, y/vmag));
	}
};

typedef point2d<float> vector2d;


template<typename T> struct pointT { // size = 12 (float), 24(double)

	typedef T value_type;
	T x, y, z;

	pointT() : x(0.0), y(0.0), z(0.0) {}
	//pointT(T v) : x(v), y(v), z(v) {} // unsafe?
	pointT(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}
	pointT(pointT const &p1, pointT const &p2) : x(p1.x-p2.x), y(p1.y-p2.y), z(p1.z-p2.z) {} // take the difference (vector)
	template<typename S> pointT(pointT<S> const &p) : x(p.x), y(p.y), z(p.z) {}

	string str() const {std::ostringstream oss; oss << x << ", " << y << ", " << z; return oss.str();}
	string raw_str() const {std::ostringstream oss; oss << x << " " << y << " " << z; return oss.str();}
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
	void operator*=(pointT const &p) {x *= p.x; y *= p.y; z *= p.z;} // component multiply
	void operator+=(T const &v)      {x += v; y += v; z += v;}
	void operator-=(T const &v)      {x -= v; y -= v; z -= v;}
	void operator*=(T m)             {x *= m; y *= m; z *= m;}

	void operator/=(T d) {
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

	void invert() {
		if (x == 0.0) {x = TOLERANCE;}
		if (y == 0.0) {y = TOLERANCE;}
		if (z == 0.0) {z = TOLERANCE;}
		x = 1.0/x; y = 1.0/y; z = 1.0/z;
	}
	pointT inverse() const {return pointT(1.0/x, 1.0/y, 1.0/z);} // no divide by zero check

	pointT get_norm() const {
		T const vmag(mag());
		return ((vmag < TOLERANCE) ? *this : pointT(x/vmag, y/vmag, z/vmag));
	}
	void set_max_mag(T vmax) {
		T const vmag(mag());
		if (vmag > TOLERANCE && vmag > vmax) operator*=(vmax/vmag);
	}
	T sum() const {return (x + y + z);}
	pointT operator+(pointT const &p)  const {return pointT((x+p.x), (y+p.y), (z+p.z));}
	pointT operator-(pointT const &p)  const {return pointT((x-p.x), (y-p.y), (z-p.z));}
	pointT operator+(T const &v)       const {return pointT((x+v), (y+v), (z+v));}
	pointT operator-(T const &v)       const {return pointT((x-v), (y-v), (z-v));}
	pointT operator*(T      const val) const {return pointT(x*val, y*val, z*val);}
	pointT operator*(pointT const &p)  const {return pointT(x*p.x, y*p.y, z*p.z);} // component multiply
	pointT operator/(pointT const &p)  const {return pointT(x/p.x, y/p.y, z/p.z);} // component division
	pointT operator-()                 const {return pointT(-x, -y, -z);}
	
	pointT operator/(T      const val) const {
		assert(val != 0.0);
		T const val_inv(1.0/val);
		return pointT(x*val_inv, y*val_inv, z*val_inv);
	}
	float  dot  (pointT const &v) const {return (x*v.x + y*v.y + z*v.z);}
	pointT cross(pointT const &v) const {return pointT((y*v.z - z*v.y), (z*v.x - x*v.z), (x*v.y - y*v.x));}
	pointT min  (pointT const &v) const {return pointT(std::min(x, v.x), std::min(y, v.y), std::min(z, v.z));}
	pointT max  (pointT const &v) const {return pointT(std::max(x, v.x), std::max(y, v.y), std::max(z, v.z));}
	T mag_sq()    const {return (x*x + y*y + z*z);}
	T mag()       const {return sqrt(mag_sq());}
	T xy_mag_sq() const {return (x*x + y*y);}
	T xy_mag()    const {return sqrt(xy_mag_sq());}
	T get_min_val()   const {return std::min(x, std::min(y, z));}
	T get_max_val()   const {return std::max(x, std::max(y, z));}
	bool is_nonzero() const {return (x != 0.0 || y != 0.0 || z != 0.0);}

	bool operator<(pointT const &p) const { // needed for maps and stuff
		if (z > p.z) return 1;
		if (z < p.z) return 0; // greater than operation?
		if (y < p.y) return 1;
		if (y > p.y) return 0;
		return (x < p.x);
	}
};

// premultiply a pointT by a scalar
template<typename S, typename T> point2d<T> inline operator*(S const v, point2d<T> const &p) {return point2d<T>(v*p.x, v*p.y);}
template<typename S, typename T> pointT <T> inline operator*(S const v, pointT <T> const &p) {return pointT <T>(v*p.x, v*p.y, v*p.z);}


typedef pointT<float>  point;
typedef pointT<float>  vector3d;
typedef pointT<double> point_d;
typedef pointT<double> vector3d_d;
typedef point_d        upos_point_type;


// constants
point    const all_zeros(0, 0, 0);
vector3d const plus_x(1, 0, 0);
vector3d const plus_y(0, 1, 0);
vector3d const plus_z(0, 0, 1);
vector3d const zero_vector(0, 0, 0);
vector3d const all_ones(1, 1, 1);


template<typename T> uint32_t jenkins_one_at_a_time_hash(const T* key, size_t length) { // T is an unsigned integer type
	size_t i = 0;
	uint32_t hash = 0;
	while (i != length) {hash += key[i++]; hash += hash << 10; hash ^= hash >> 6;}
	hash += hash << 3;
	hash ^= hash >> 11;
	hash += hash << 15;
	return hash;
}

template<typename T> struct hash_by_bytes { // should work with all packed vertex types
	uint32_t operator()(T const &v) const {return jenkins_one_at_a_time_hash((const uint8_t*)&v, sizeof(T));} // slower but better quality hash
	//uint32_t operator()(T const &v) const {return jenkins_one_at_a_time_hash((const uint32_t*)&v, sizeof(T)>>2);} // faster but lower quality hash
};
inline unsigned hash_point(point const &p) {return hash_by_bytes<point>()(p);}

inline void hash_mix_point(point const &p, unsigned &hv) {
	hv += hash_point(p);
	hv += hv << 10;
	hv ^= hv >> 6;
}
template<typename T> unsigned hash_vect_as_int(vector<T> const &v) {
	assert((sizeof(T) % sizeof(int)) == 0); // must be a multiple of 4 bytes
	return jenkins_one_at_a_time_hash((int const*)v.data(), sizeof(T)*v.size()/sizeof(int));
}


struct vector4d : public vector3d { // size = 16
	float w;

	vector4d() : w(0.0) {}
	vector4d(float x_, float y_, float z_, float w_) : vector3d(x_, y_, z_), w(w_) {}
	vector4d(vector3d const &v, float w_) : vector3d(v), w(w_) {}
	void assign(float x_, float y_, float z_, float w_)    {x = x_; y = y_; z = z_; w = w_;}
	string str() const {std::ostringstream oss; oss << x << ", " << y << ", " << z << ", " << w; return oss.str();}
	vector4d operator+ (vector4d const &p) const {return vector4d((x+p.x), (y+p.y), (z+p.z), (w+p.w));}
	vector4d operator- (vector4d const &p) const {return vector4d((x-p.x), (y-p.y), (z-p.z), (w-p.w));}
	void     operator+=(vector4d const &p) {x += p.x; y += p.y; z += p.z; w += p.w;}
	void     operator-=(vector4d const &p) {x -= p.x; y -= p.y; z -= p.z; w -= p.w;}
	vector4d operator-()                   const {return vector4d(-x, -y, -z, -w);}
	bool     operator==(vector4d const &v) const {return (v.x == x && v.y == y && v.z == z && v.w == w);}
	bool     operator!=(vector4d const &v) const {return !operator==(v);}
	bool     operator< (vector4d const &p) const {return ((w == p.w) ? vector3d::operator<(p) : (w < p.w));}
};


struct sphere_t {

	point pos;
	float radius;
	
	sphere_t(point const &p=all_zeros, float r=0.0) : pos(p), radius(r) {}
	bool operator==(sphere_t const &s) const {return (pos == s.pos && radius == s.radius);}
	bool operator!=(sphere_t const &s) const {return (pos != s.pos || radius != s.radius);}
	point const &get_pos() const {return pos;}
	float get_radius()     const {return radius;}
	float get_volume()     const {return (4.0/3.0)*PI*radius*radius*radius;}
	float get_surf_area()  const {return 4.0*PI*radius*radius;}
	bool contains_point(point const &p) const;
};


struct cube_t { // size = 24; Note: AABB, not actually a cube

	float d[3][2]; // {x,y,z},{min,max}

	cube_t() {set_to_zeros();}
	
	cube_t(float x1_, float x2_, float y1_, float y2_, float z1_, float z2_) {
		x1() = x1_; x2() = x2_; y1() = y1_; y2() = y2_; z1() = z1_; z2() = z2_;
	}
	cube_t(point const &p1, point const &p2) {UNROLL_3X(d[i_][0] = min(p1[i_], p2[i_]); d[i_][1] = max(p1[i_], p2[i_]);)}
	cube_t(point const &pt) {set_from_point(pt);}
	cube_t(point const *const pts, unsigned npts) {set_from_points(pts, npts);}
	void set_to_zeros() {set_from_point(all_zeros);}
	void copy_from(cube_t const &c) {
		UNROLL_3X(d[i_][0] = c.d[i_][0]; d[i_][1] = c.d[i_][1];)
	}
	void set_from_point(point const &pt) {UNROLL_3X(d[i_][0] = d[i_][1] = pt[i_];)}
	void set_from_sphere(point const &pt, float radius) {
		UNROLL_3X(d[i_][0] = pt[i_]-radius; d[i_][1] = pt[i_]+radius;)
	}
	void set_from_sphere(sphere_t const &s) {set_from_sphere(s.pos, s.radius);}
	bool operator==(cube_t const &c) const {
		UNROLL_3X(if (d[i_][0] != c.d[i_][0]) return 0;)
		UNROLL_3X(if (d[i_][1] != c.d[i_][1]) return 0;)
		return 1;
	}
	bool operator<(cube_t const &c) const { // for use in sorting/uniquing
		UNROLL_3X(if (d[i_][0] != c.d[i_][0]) return (d[i_][0] < c.d[i_][0]);)
		UNROLL_3X(if (d[i_][1] != c.d[i_][1]) return (d[i_][1] < c.d[i_][1]);)
		return 0;
	}
	float x1() const {return d[0][0];}
	float x2() const {return d[0][1];}
	float y1() const {return d[1][0];}
	float y2() const {return d[1][1];}
	float z1() const {return d[2][0];}
	float z2() const {return d[2][1];}
	float &x1() {return d[0][0];}
	float &x2() {return d[0][1];}
	float &y1() {return d[1][0];}
	float &y2() {return d[1][1];}
	float &z1() {return d[2][0];}
	float &z2() {return d[2][1];}
	float xc() const {return 0.5f*(x1() + x2());}
	float yc() const {return 0.5f*(y1() + y2());}
	float zc() const {return 0.5f*(z1() + z2());}
	float dx() const {return (x2() - x1());}
	float dy() const {return (y2() - y1());}
	float dz() const {return (z2() - z1());}
	bool operator!=(cube_t const &c) const {return !operator==(c);}
	cube_t operator+ (vector3d const &p) const {cube_t c(*this); c += p; return c;}
	cube_t operator- (vector3d const &p) const {cube_t c(*this); c -= p; return c;}
	cube_t operator* (vector3d const &p) const {cube_t c(*this); c *= p; return c;}
	cube_t operator* (float scale      ) const {cube_t c(*this); c *= scale; return c;}
	void   operator+=(vector3d const &p) {translate( p);}
	void   operator-=(vector3d const &p) {translate(-p);}
	void   operator*=(vector3d const &p) {UNROLL_3X(d[i_][0] *= p[i_]; d[i_][1] *= p[i_];)}
	void   operator*=(float scale      ) {UNROLL_3X(d[i_][0] *= scale; d[i_][1] *= scale;)}

	void translate(point const &p) {UNROLL_3X(d[i_][0] += p[i_]; d[i_][1] += p[i_];)}
	void translate_dim(unsigned dim, float v) {assert(dim < 3); d[dim][0] += v; d[dim][1] += v;}
	void swap_dims(unsigned d1, unsigned d2) {assert(d1 < 3 && d2 < 3); swap(d[d1][0], d[d2][0]); swap(d[d1][1], d[d2][1]);}
	void set_from_points(point const *const pts, unsigned npts);
	void set_from_points(vector<point> const &pts) {set_from_points(pts.data(), pts.size());}
	string str() const;
	string raw_str() const;
	bool is_near_zero_area() const;
	bool is_all_zeros() const {return (x1() == 0 && x2() == 0 && y1() == 0 && y2() == 0 && z1() == 0 && z2() == 0);}

	void union_with_pt(point const &pt) {
		UNROLL_3X(min_eq(d[i_][0], pt[i_]); max_eq(d[i_][1], pt[i_]);)
	}
	void assign_or_union_with_pt(point const &pt) {
		if (is_all_zeros()) {set_from_point(pt);} else {union_with_pt(pt);} // Note: won't work if pt == (0,0,0)
	}
	void assign_or_union_with_sphere(point const &pt, float radius) {
		if (is_zero_area()) {set_from_sphere(pt, radius);} else {union_with_sphere(pt, radius);} // Note: won't work if pt == (0,0,0)
	}
	void union_with_sphere(point const &pt, float radius) {
		UNROLL_3X(min_eq(d[i_][0], pt[i_]-radius); max_eq(d[i_][1], pt[i_]+radius);)
	}
	void union_with_sphere(sphere_t const &s) {union_with_sphere(s.pos, s.radius);}
	void union_with_cube(cube_t const &c) {
		UNROLL_3X(min_eq(d[i_][0], c.d[i_][0]); max_eq(d[i_][1], c.d[i_][1]);)
	}
	void union_with_cube_xy(cube_t const &c) {
		UNROLL_2X(min_eq(d[i_][0], c.d[i_][0]); max_eq(d[i_][1], c.d[i_][1]);)
	}
	void assign_or_union_with_cube(cube_t const &c) {
		if (c.is_zero_area()) return;
		if (is_zero_area()) {copy_from(c);} else {union_with_cube(c);}
	}
	void intersect_with_cube(cube_t const &c) { // Note: cube and *this must overlap
		UNROLL_3X(max_eq(d[i_][0], c.d[i_][0]); min_eq(d[i_][1], c.d[i_][1]);)
	}
	void intersect_with_cube_xy(cube_t const &c) { // Note: cube and *this must overlap
		UNROLL_2X(max_eq(d[i_][0], c.d[i_][0]); min_eq(d[i_][1], c.d[i_][1]);)
	}
	void normalize() {UNROLL_3X(if (d[i_][1] < d[i_][0]) swap(d[i_][0], d[i_][1]);)}

	bool is_zero_area() const {
		UNROLL_3X(if (d[i_][0] == d[i_][1]) return 1;)
		return 0;
	}
	bool is_normalized() const {
		UNROLL_3X(if (d[i_][0] > d[i_][1]) return 0;)
		return 1;
	}
	bool is_strictly_normalized() const {
		UNROLL_3X(if (d[i_][0] >= d[i_][1]) return 0;)
		return 1;
	}
	bool intersects(const cube_t &cube) const { // includes adjacency
		UNROLL_3X(if (cube.d[i_][1] < d[i_][0] || cube.d[i_][0] > d[i_][1]) return 0;)
		return 1;
	}
	bool intersects_no_adj(const cube_t &cube) const { // excludes adjacency
		UNROLL_3X(if (cube.d[i_][1] <= d[i_][0] || cube.d[i_][0] >= d[i_][1]) return 0;)
		return 1;
	}
	bool intersects_xy(const cube_t &cube) const {
		UNROLL_2X(if (cube.d[i_][1] < d[i_][0] || cube.d[i_][0] > d[i_][1]) return 0;)
		return 1;
	}
	bool intersects_xy_no_adj(const cube_t &cube) const { // excludes adjacency
		UNROLL_2X(if (cube.d[i_][1] <= d[i_][0] || cube.d[i_][0] >= d[i_][1]) return 0;)
		return 1;
	}
	bool intersects(const cube_t &cube, float toler) const { // Note: toler > 0 makes adjacent cubes *not* intersect
		UNROLL_3X(if (cube.d[i_][1] < (d[i_][0] + toler) || cube.d[i_][0] > (d[i_][1] - toler)) return 0;)
		return 1;
	}
	bool contains_cube(const cube_t &cube) const {
		UNROLL_3X(if (cube.d[i_][0] < d[i_][0] || cube.d[i_][1] > d[i_][1]) return 0;)
		return 1;
	}
	bool contains_cube_xy(const cube_t &cube) const {
		UNROLL_2X(if (cube.d[i_][0] < d[i_][0] || cube.d[i_][1] > d[i_][1]) return 0;)
		return 1;
	}
	bool contains_cube_xy_exp(const cube_t &cube, float exp) const {
		UNROLL_2X(if (cube.d[i_][0]-exp < d[i_][0] || cube.d[i_][1]+exp > d[i_][1]) return 0;)
			return 1;
	}
	bool contains_cube_xy_no_adj(const cube_t &cube) const {
		UNROLL_2X(if (cube.d[i_][0] <= d[i_][0] || cube.d[i_][1] >= d[i_][1]) return 0;)
		return 1;
	}
	bool contains_cube_xy_overlaps_z(const cube_t &cube) const { // no adj in Z
		return (contains_cube_xy(cube) && cube.z1() < z2() && cube.z2() > z1());
	}
	bool contains_pt(point const &pt) const { // includes points on the edge
		UNROLL_3X(if (pt[i_] < d[i_][0] || pt[i_] > d[i_][1]) return 0;)
		return 1;
	}
	bool contains_pt_xy             (point const &pt) const {return (pt.x >  x1() && pt.x <  x2() && pt.y >  y1() && pt.y <  y2());} // excludes points on the edge
	bool contains_pt_xy_inc_low_edge(point const &pt) const {return (pt.x >= x1() && pt.x <  x2() && pt.y >= y1() && pt.y <  y2());} // includes points on the lower edges
	bool contains_pt_xy_inclusive   (point const &pt) const {return (pt.x >= x1() && pt.x <= x2() && pt.y >= y1() && pt.y <= y2());} // includes points on the edge
	bool contains_pt_xy_exp         (point const &pt, float exp) const {return (pt.x > x1()-exp && pt.x < x2()+exp && pt.y > y1()-exp && pt.y < y2()+exp);}
	bool contains_pt_exp            (point const &pt, float exp) const {return (pt.x > x1()-exp && pt.x < x2()+exp && pt.y > y1()-exp && pt.y < y2()+exp && pt.z > z1()-exp && pt.z < z2()+exp);}
	bool contains_pt_exp_xy_only    (point const &pt, float exp) const {return (pt.x > x1()-exp && pt.x < x2()+exp && pt.y > y1()-exp && pt.y < y2()+exp && pt.z > z1() && pt.z < z2());}

	bool quick_intersect_test(const cube_t &cube) const {
		UNROLL_3X(if (cube.d[i_][0] >= d[i_][1] || cube.d[i_][1] <= d[i_][0]) return 0;)
		return 1;
	}
	bool line_intersects(point const &p1, point const &p2) const;
	void clamp_pt   (point &pt) const {UNROLL_3X(pt[i_] = min(d[i_][1], max(d[i_][0], pt[i_]));)}
	void clamp_pt_xy(point &pt) const {UNROLL_2X(pt[i_] = min(d[i_][1], max(d[i_][0], pt[i_]));)}
	float get_volume() const {return fabs(x2() - x1())*fabs(y2() - y1())*fabs(z2() - z1());}
	float get_area  () const {return 2.0f*(fabs(x2() - x1())*fabs(y2() - y1()) + fabs(y2() - y1())*fabs(z2() - z1()) + fabs(z2() - z1())*fabs(x2() - x1()));}
	float get_area_xy()const {return dx()*dy();}
	float max_len   () const {return max((x2() - x1()), max((y2() - y1()), (z2() - z1())));}
	float min_len   () const {return min((x2() - x1()), min((y2() - y1()), (z2() - z1())));}

	float second_largest_len() const {
		return min(max((x2() - x1()), (y2() - y1())), min(max((y2() - y1()), (z2() - z1())),  max((z2() - z1()), (x2() - x1()))));
	}
	point get_cube_center() const {
		return point(0.5f*(x1()+x2()), 0.5f*(y1()+y2()), 0.5f*(z1()+z2()));
	}
	float get_bsphere_radius() const {
		return 0.5f*sqrt((x2()-x1())*(x2()-x1()) + (y2()-y1())*(y2()-y1()) + (z2()-z1())*(z2()-z1()));
	}
	float get_xy_bsphere_radius() const {
		return 0.5f*sqrt((x2()-x1())*(x2()-x1()) + (y2()-y1())*(y2()-y1()));
	}
	sphere_t get_bsphere() const {return sphere_t(get_cube_center(), get_bsphere_radius());}
	sphere_t get_bcylin () const {return sphere_t(get_cube_center(), get_xy_bsphere_radius());}
	point get_llc() const {return point(x1(), y1(), z1());}
	point get_urc() const {return point(x2(), y2(), z2());}
	vector3d get_size   () const {return vector3d((x2()-x1()), (y2()-y1()), (z2()-z1()));}
	vector2d get_size_xy() const {return vector2d((x2()-x1()), (y2()-y1()));}
	float get_center_dim(unsigned dim) const {assert(dim < 3); return 0.5f*(d[dim][0] + d[dim][1]);}
	float get_sz_dim    (unsigned dim) const {assert(dim < 3); return (d[dim][1] - d[dim][0]);}
	void expand_by(float val) {UNROLL_3X(d[i_][0] -= val; d[i_][1] += val;)}
	void expand_by(float x, float y, float z) {expand_by(vector3d(x, y, z));}
	void expand_by(vector3d const &val) {UNROLL_3X(d[i_][0] -= val[i_]; d[i_][1] += val[i_];)}
	void expand_by_xy(float val) {UNROLL_2X(d[i_][0] -= val; d[i_][1] += val;)}
	void expand_by_xy(vector3d const &val) {UNROLL_2X(d[i_][0] -= val[i_]; d[i_][1] += val[i_];)}
	void expand_by_xy(vector2d const &val) {UNROLL_2X(d[i_][0] -= val[i_]; d[i_][1] += val[i_];)}
	void expand_in_dim(unsigned dim, float val) {assert(dim < 3); d[dim][0] -= val; d[dim][1] += val;}
	void expand_in_x(float val) {x1() -= val; x2() += val;}
	void expand_in_y(float val) {y1() -= val; y2() += val;}
	void expand_in_z(float val) {z1() -= val; z2() += val;}
	unsigned get_split_dim(float &max_sz, float &sval, unsigned skip_dims) const;
	bool cube_intersection(const cube_t &cube, cube_t &res) const;
	float get_overlap_volume(const cube_t &cube) const;
	vector3d closest_side_dir(point const &pos, unsigned skip_dims=0) const;
	bool closest_dist_less_than(point const &pos, float dist) const;
	bool closest_dist_xy_less_than(point const &pos, float dist) const;
	
	point closest_pt(point const &pos) const { // closest point inside this cube
		point pt(pos);
		clamp_pt(pt);
		return pt;
	}
	float get_max_extent() const { // from (0,0,0)
		float mextent(0.0);
		UNROLL_3X(mextent = max(mextent, max(-d[i_][0], d[i_][1]));)
		return mextent;
	}
	float get_max_dim_sz() const {return max(dz(), max(dx(), dy()));}
	float get_min_dim_sz() const {return min(dz(), min(dx(), dy()));}

	float furthest_dist_to_pt(point const &pos) const {
		vector3d dmax;
		UNROLL_3X(dmax[i_] = max((pos[i_] - d[i_][0]), (d[i_][1] - pos[i_]));)
		return dmax.mag();
	}
	int closest_face(point const &pos) const;
	bool cube_merge(cube_t const &cube);
	void get_points(point pts[8]) const;
};

typedef vector<cube_t> vect_cube_t;
cube_t const all_zeros_cube(0,0,0,0,0,0);

struct cube_with_ix_t : public cube_t {
	unsigned ix;
	cube_with_ix_t(unsigned ix_=0) : ix(ix_) {}
	cube_with_ix_t(cube_t const &c, unsigned ix_=0) : cube_t(c), ix(ix_) {}
};
typedef vector<cube_with_ix_t> vect_cube_with_ix_t;


vector3d get_poly_norm(point const *const points, bool normalize=1);

struct tquad_t { // size = 52

	point pts[4];
	unsigned npts;
	
	tquad_t(unsigned npts_=0) : npts(npts_) {}
	bool is_valid() const;
	void update_bcube(cube_t &c) const;
	cube_t get_bcube() const;
	vector3d get_norm(bool normalize=1) const {return get_poly_norm(pts, normalize);}
	point const &operator[](unsigned i) const {return pts[i];}
	point       &operator[](unsigned i)       {return pts[i];}
};


struct line_3dw {

	point p1, p2;

	line_3dw() : p1(all_zeros), p2(all_zeros) {}
	line_3dw(point const &p1_, point const &p2_) : p1(p1_), p2(p2_) {assert(p1 != p2);}
	vector3d get_norm_dir_vect() const {return (p2 - p1).get_norm();}
	float get_length() const {return (p1 - p2).mag();}
	void translate(point const &p) {p1 += p; p2 += p;}
};


struct vector_point_norm {
	vector<point>    p;
	vector<vector3d> n;
};


struct pos_dir_up { // defines a view frustum

	point pos;
	vector3d dir, upv, upv_, cp;
	float angle, tterm, sterm, x_sterm, behind_sphere_mult, near_, far_;
	double A; // aspect ratio x/y
	bool valid;

	pos_dir_up(void) : angle(0.0f), tterm(0.0f), sterm(0.0f), x_sterm(0.0f), behind_sphere_mult(0.0f), near_(0.0f), far_(0.0f), A(0.0), valid(0) {}
	pos_dir_up(point const &p, vector3d const &d, vector3d const &u, float angle_, float n, float f, float a=0.0, bool no_zoom=0);
	void orthogonalize_up_dir();
	bool point_visible_test(point const &pos_) const;
	bool line_visible_test(point const &p1, point const &p2) const;
	bool sphere_visible_test(point const &pos_, float radius) const;
	bool sphere_visible_test_no_inside_test(point const &pos_, float radius) const;
	bool sphere_completely_visible_test(point const &pos_, float radius) const {return sphere_visible_test(pos_, -radius);}
	template<unsigned N> bool pt_set_visible(point const *const pts) const;
	bool cube_visible(cube_t const &c) const;
	bool cube_completely_visible(cube_t const &c) const;
	bool cube_visible_likely(cube_t const &c) const {return (!valid || point_visible_test(c.get_cube_center()) || cube_visible(c));}
	bool cube_visible_for_light_cone(cube_t const &c) const;
	bool projected_cube_visible(cube_t const &cube, point const &proj_pt) const;
	bool sphere_and_cube_visible_test(point const &pos_, float radius, cube_t const &cube) const;
	void get_frustum_corners(point pts[8]) const;
	point get_frustum_center() const;
	void draw_frustum() const;
	void translate(vector3d const &tv) {pos += tv;}
	void scale(float s) {pos *= s; near_ *= s; far_ *= s;}
	void rotate(vector3d const &axis, float angle);
	void apply_z_mirror(float zval) {apply_dim_mirror(2, zval);}
	void apply_dim_mirror(unsigned dim, float val);
};


struct cylinder_3dw : public line_3dw { // size = 32

	float r1, r2;

	cylinder_3dw() : r1(0.0), r2(0.0) {}
	cylinder_3dw(point const &p1_, point const &p2_, float r1_, float r2_) : line_3dw(p1_, p2_), r1(r1_), r2(r2_) {}
	void calc_bcube(cube_t &bcube) const;
	float get_volume() const {return PI*(r1*r1 + r1*r2 + r2*r2)*get_length()/3.0f;}
	float get_surface_area() const;
	point get_center() const {return 0.5f*(p1 + p2);}
	float get_avg_radius() const {return 0.5f*(r1 + r2);}
	float get_bounding_radius() const;
};


struct colorRGB { // size = 12

	float R, G, B;
	colorRGB() : R(0.0f), G(0.0f), B(0.0f) {}
	colorRGB(float r, float g, float b) : R(r), G(g), B(b) {}
	void assign(float r, float g, float b) {R = r; G = g; B = b;}
	void set_to_val(float val) {R = G = B = val;}
	bool operator==(const colorRGB &c) const {return (c.R == R && c.G == G && c.B == B);}
	bool operator!=(const colorRGB &c) const {return !operator==(c);}

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
	bool operator<(const colorRGB &c) const { // greater than operation?
		if (R < c.R) return 1;
		if (R > c.R) return 0;
		if (G < c.G) return 1;
		if (G > c.G) return 0;
		return (B < c.B);
	}
	colorRGB operator+ (colorRGB const &c) const {return colorRGB(R+c.R, G+c.G, B+c.B);}
	void     operator+=(colorRGB const &c)       {R += c.R; G += c.G; B += c.B;}
	colorRGB operator*(float val) const {return colorRGB(R*val, G*val, B*val);}
	void operator*=(float val) {R *= val; G *= val; B *= val;}
	colorRGB modulate_with(colorRGB const &c) const {return colorRGB(R*c.R, G*c.G, B*c.B);}

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
	void normalize_to_max_comp() {
		float const max_comp(get_max_component());
		if (max_comp > TOLERANCE) {R /= max_comp; G /= max_comp; B /= max_comp;}
	}
	string str    () const {std::ostringstream oss; oss << "R: " << R << ", G: " << G << ", B: " << B; return oss.str();}
	string raw_str() const {std::ostringstream oss; oss << R << " " << G << " " << B; return oss.str();}
	float get_luminance() const {return (R + G + B)/3.0f;}
	float get_weighted_luminance() const {return (0.2126*R + 0.7152*G + 0.0722*B);} // see https://www.w3.org/WAI/GL/wiki/Relative_luminance
	float get_max_component() const {return max(R, max(G, B));}
	void set_for_cur_shader() const;
};


struct colorRGBA : public colorRGB { // size = 16

	union {float A; float alpha;}; // A and alpha are both valid components

	colorRGBA() : alpha(1.0) {}
	colorRGBA(float r, float g, float b, float a=1.0) : colorRGB(r, g, b), A(a) {}
	colorRGBA(colorRGB const &color, float a=1.0) : colorRGB(color), A(a) {}
	void assign(float r, float g, float b, float a=1.0) {R = r; G = g; B = b; A = a;}
	bool operator==(const colorRGBA &c) const {return (c.R == R && c.G == G && c.B == B && c.A == A);}
	bool operator!=(const colorRGBA &c) const {return !operator==(c);}

	const float &operator[](unsigned i) const {
		switch(i) {
			case 0: return R;
			case 1: return G;
			case 2: return B;
			case 3: return A;
			default: assert(0);
		}
		return R; // never gets here
	}
	float &operator[](unsigned i) {
		switch(i) {
			case 0: return R;
			case 1: return G;
			case 2: return B;
			case 3: return A;
			default: assert(0);
		}
		return R; // never gets here
	}
	bool operator<(const colorRGBA &c) const { // greater than operation?
		if (A > c.A) return 1; // note: alpha less than so that low alpha colors are last
		if (A < c.A) return 0;
		return colorRGB::operator<(c);
	}
	colorRGBA operator* (float val) const             {return colorRGBA(R*val, G*val, B*val, A);}
	colorRGBA operator/ (float val) const             {return colorRGBA(R/val, G/val, B/val, A);}
	colorRGBA operator+ (colorRGBA const &c) const    {return colorRGBA(R+c.R, G+c.G, B+c.B, A+c.A);}
	void      operator+=(colorRGBA const &c)          {R += c.R; G += c.G; B += c.B; A += c.A;}
	colorRGBA modulate_with(colorRGBA const &c) const {return colorRGBA(R*c.R, G*c.G, B*c.B, A*c.A);}

	void set_valid_color() {
		colorRGB::set_valid_color();
		A = CLIP_TO_01(A);
	}
	void normalize_to_alpha_1() {
		if (A == 1.0) return;
		R *= A; G *= A; B *= A;
		A = 1.0;
	}
	bool within_thresh_of_rgb(float thresh, colorRGBA const &c) const { // no alpha check
		return ((fabs(R-c.R) + fabs(G-c.G) + fabs(B-c.B)) < thresh);
	}
	bool within_thresh_of_rgba(float thresh, colorRGBA const &c) const { // no alpha check
		return ((fabs(R-c.R) + fabs(G-c.G) + fabs(B-c.B) + fabs(A-c.A)) < thresh);
	}
	bool is_valid() const {return (R >= 0 && G >= 0 && B >= 0 && A >= 0 && R <= 1 && G <= 1 && B <= 1 && A <= 1);}
	string str    () const {std::ostringstream oss; oss << "R: " << R << ", G: " << G << ", B: " << B << ", A: " << A; return oss.str();}
	string raw_str() const {std::ostringstream oss; oss << R << " " << G << " " << B << " " << A; return oss.str();}
	void set_for_cur_shader() const;
};


struct tex_range_t {

	float x1, y1, x2, y2;
	bool clip_quad, swap_xy;

	tex_range_t() : x1(0.0), y1(0.0), x2(1.0), y2(1.0), clip_quad(0), swap_xy(0) {}
	tex_range_t(float x1_, float y1_, float x2_, float y2_, bool clip_quad_=0, bool swap_xy_=0) : x1(x1_), y1(y1_), x2(x2_), y2(y2_), clip_quad(clip_quad_), swap_xy(swap_xy_) {}
	void mirror_x() {swap(x1, x2);}
	void mirror_y() {swap(y1, y2);}

	static tex_range_t from_atlas(unsigned xv, unsigned yv, unsigned nx, unsigned ny) {
		assert(nx > 0 && ny > 0 && xv < nx && yv < ny);
		return tex_range_t(xv/float(nx), yv/float(ny), (xv+1)/float(nx), (yv+1)/float(ny));
	}
};

#include "vertex_types.h" // must be included here


bool bind_temp_vbo_from_verts(void const *const verts, unsigned count, unsigned vert_size, void const *&vbo_ptr_offset);
void unbind_temp_vbo();

template< typename T> void set_ptr_state(T const *const verts, unsigned count, unsigned start_ix=0, bool set_array_client_state=1) {
	void const *ptr_offset = NULL;
	if (verts && !bind_temp_vbo_from_verts(verts+start_ix, count, sizeof(T), ptr_offset)) {ptr_offset = verts + start_ix;}
	T::set_vbo_arrays(set_array_client_state, ptr_offset);
}
template <typename T> void unset_ptr_state(T const *const verts) {
	T::unset_attrs();
	if (verts) {unbind_temp_vbo();}
}

extern unsigned num_frame_draw_calls;
template <typename T> void draw_verts(T const *const verts, unsigned count, int gl_type, unsigned start_ix=0, bool set_array_client_state=1) {
	assert(count > 0);
	set_ptr_state(verts, count, start_ix, set_array_client_state);
	glDrawArrays(gl_type, start_ix, count);
	++num_frame_draw_calls;
	unset_ptr_state(verts);
}
template <typename T> void draw_verts(vector<T> const &verts, int gl_type, unsigned start_ix=0, bool set_array_client_state=1) {
	if (!verts.empty()) {draw_verts(&verts.front(), verts.size(), gl_type, start_ix, set_array_client_state);}
}

template <typename T> void draw_and_clear_verts(vector<T> &verts, int gl_type) {
	draw_verts(verts, gl_type);
	verts.resize(0); // clear()?
}

void draw_quads_as_tris(unsigned num_quad_verts, unsigned start_quad_vert=0, unsigned num_instances=1);
bool bind_quads_as_tris_ivbo(unsigned num_quad_verts);
void convert_quad_ixs_to_tri_ixs(vector<unsigned> const &qixs, vector<unsigned> &tixs);

template <typename T> void draw_quad_verts_as_tris(T const *const verts, unsigned count, unsigned start_ix=0, unsigned num_instances=1, bool set_array_client_state=1) {
	assert(count > 0);
	set_ptr_state(verts, count, start_ix, set_array_client_state);
	draw_quads_as_tris(count, start_ix, num_instances);
	unset_ptr_state(verts);
}
template <typename T> void draw_quad_verts_as_tris(vector<T> const &verts, unsigned start_ix=0, unsigned num_instances=1, bool set_array_client_state=1) {
	if (!verts.empty()) {draw_quad_verts_as_tris(&verts.front(), verts.size(), start_ix, num_instances, set_array_client_state);}
}

template <typename T> void draw_quad_verts_as_tris_and_clear(vector<T> &verts) {
	draw_quad_verts_as_tris(verts); verts.clear();
}

extern bool use_core_context;
template<typename T> void draw_vect_quads(T const &verts) {
	if (use_core_context) {draw_quad_verts_as_tris(verts);} else {draw_verts(verts, GL_QUADS);}
}


template <typename T> void translate_verts(vector<T> &verts, vector3d const &xlate) {
	for (auto i = verts.begin(); i != verts.end(); ++i) {i->v += xlate;}
}

template <typename T> void scale_verts(vector<T> &verts, vector3d const &scale) {
	for (auto i = verts.begin(); i != verts.end(); ++i) {i->v *= scale;}
}

template <typename T> void tri_strip_push(vector<T> &v) {
	assert(v.size() >= 3);
	v.push_back(v[v.size()-2]);
	v.push_back(v[v.size()-2]);
}


class draw_call_counter {
	string name;
	unsigned start_num_draw_calls;
public:
	draw_call_counter(string const &name_) : name(name_), start_num_draw_calls(num_frame_draw_calls) {}
	~draw_call_counter() {std::cout << name << ": " << (num_frame_draw_calls - start_num_draw_calls) << std::endl;}
};


template<typename T> struct triangle_t {

	T pts[3];

	triangle_t() {}
	triangle_t(T const &p1, T const &p2, T const &p3) {pts[0] = p1; pts[1] = p2; pts[2] = p3;}
	triangle_t(T const *const p) {pts[0] = p[0]; pts[1] = p[1]; pts[2] = p[2];}
	vector3d get_normal(bool normalize=1) const {return get_poly_norm(pts, normalize);}
	void operator+=(T const &p) {pts[0] += p; pts[1] += p; pts[2] += p;}

	cube_t get_bbox(vector<T> const &p) const { // Note: only works on some types of triangles
		cube_t bbox(pts[0], pts[0]);
		bbox.union_with_pt(pts[1]);
		bbox.union_with_pt(pts[2]);
		return bbox;
	}
};

typedef triangle_t<point>        triangle;
typedef triangle_t<point_d>      triangle_d; // unused
typedef triangle_t<vert_norm_tc> triangle_vntc;


struct ray3d { // size = 40

	point pts[2];
	colorRGBA color;
	
	ray3d(point const &pt0, point const &pt1, colorRGBA const &c) : color(c) {pts[0] = pt0; pts[1] = pt1;}
	ray3d() {}
};


class line_tquad_draw_t;

struct beam3d : public ray3d { // size = 48

	bool distant;
	short shooter;
	float intensity;

	beam3d(bool dist, int shoot, point const &pt0, point const &pt1, colorRGBA const &c, float int_=1.0)
		: ray3d(pt0, pt1, c), distant(dist), shooter(shoot), intensity(int_) {}
	beam3d() : distant(0), shooter(0), intensity(0.0f) {}
	void draw(line_tquad_draw_t &drawer) const;
};


class line3d { // size = 28

public:
	float width;
	vector<point> points;
	colorRGBA color;

	line3d() : width(0.0), color(0,0,0,0) {}
	void draw_lines(bool fade_ends, bool no_end_draw=0) const;
	bool empty() const {return points.empty();}
};


colorRGBA const DEF_TEX_COLOR(0.0, 0.0, 0.0, 0.0); // black with alpha of 0.0

// format: 0: RAW, 1: BMP, 2: RAW (upside down), 3: RAW (alpha channel), 4: targa (*tga), 5: jpeg, 6: png, 7: auto, 8: tiff, 10: DDS, 11:ppm, 12: tex2d
enum {IMG_FMT_RAW_RGB=0, IMG_FMT_BMP, IMG_FMT_RAW_INVY, IMG_FMT_RAW_RGBA, IMG_FMT_TGA, IMG_FMT_JPG, IMG_FMT_PNG, IMG_FMT_AUTO,
	IMG_FMT_TIFF, IMG_FMT_GEN, IMG_FMT_DDS, IMG_FMT_PPM, IMG_FMT_TEX2D, IMG_FMT_OTHER};


class texture_t { // size >= 116

public:
	char type=0, format=0, use_mipmaps=0, defer_load_type=DEFER_TYPE_NONE;
	bool wrap=0, mirror=0, invert_y=0, do_compress=0, has_binary_alpha=0, is_16_bit_gray=0, no_avg_color_alpha_fill=0, invert_alpha=0, normal_map=0;
	int width=0, height=0, ncolors=0, bump_tid=-1, alpha_tid=-1;
	float anisotropy=1.0, mipmap_alpha_weight=1.0;
	string name;

protected:
	unsigned char *data=nullptr, *orig_data=nullptr, *colored_data=nullptr;
	unsigned tid=0;
	colorRGBA color=DEF_TEX_COLOR;
	enum {DEFER_TYPE_NONE=0, DEFER_TYPE_DDS, DEFER_TYPE_TEX2D, NUM_DEFER_TYPE};

	void maybe_swap_rb(unsigned char *ptr) const;

public:
	texture_t() {}
	texture_t(char t, char f, int w, int h, int wrap_mir, int nc, char um, string const &n, bool inv=0, bool do_comp=1, float a=1.0, float maw=1.0, bool nm=0)
		: type(t), format(f), use_mipmaps(um), wrap(wrap_mir != 0), mirror(wrap_mir == 2), invert_y(inv), do_compress(do_comp),
		normal_map(nm), width(w), height(h), ncolors(nc), anisotropy(a), mipmap_alpha_weight(maw), name(n) {}
	bool is_inverted_y_type() const {return (defer_load_type == DEFER_TYPE_DDS);}
	void set_existing_tid(unsigned tid_, colorRGBA const &color_) {tid = tid_; color = color_;}
	void set_16_bit_grayscale();
	void init() {calc_color();}
	void do_gl_init(bool free_after_upload=0);
	void compress_and_send_texture_with_mipmaps();
	void write_texture2d_binary(string const &fn="") const;
	void read_texture2d_binary();
	void upload_cube_map_face(unsigned ix);
	bool is_texture_compressed() const;
	GLenum calc_internal_format() const;
	GLenum calc_format() const;
	GLenum get_data_format() const {return (is_16_bit_gray ? GL_UNSIGNED_SHORT : GL_UNSIGNED_BYTE);}
	void calc_color();
	void copy_alpha_from_texture(texture_t const &at, bool alpha_in_red_comp);
	void merge_in_alpha_channel(texture_t const &at);
	void set_to_color(colorRGBA const &c);
	void maybe_assign_normal_map_tid(int nm_tid) {if (nm_tid >= 0 && bump_tid < 0) {bump_tid = nm_tid;}}
	void alloc();
	void bind_gl(unsigned tu_id=0) const;
	GLuint64 get_bindless_handle(bool make_resident=1) const;
	void free_client_mem();
	void free_data() {gl_delete(); free_client_mem();}
	void gl_delete();
	void load(int index, bool allow_diff_width_height=0, bool allow_two_byte_grayscale=0, bool ignore_word_alignment=0);
	void set_image_size(int w, int h, bool allow_diff_width_height);
	void load_raw_bmp(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale);
	void load_targa(int index, bool allow_diff_width_height);
	void load_jpeg(int index, bool allow_diff_width_height);
	void load_png(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale);
	void load_tiff(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale);
	bool load_stb_image(int index, bool allow_diff_width_height, bool allow_two_byte_grayscale=0, unsigned char const *const load_from_data=nullptr, unsigned load_from_size=0);
	void load_dds();
	void deferred_load_dds();
	void load_ppm(int index, bool allow_diff_width_height);
	void set_default_white_texture();
	void auto_insert_alpha_channel(int index);
	void fill_to_grayscale_color(unsigned char color_val);
	void fill_transparent_with_avg_color();
	void do_invert_y();
	void fix_word_alignment();
	void add_alpha_channel();
	void expand_grayscale_to_rgb();
	void resize(int new_w, int new_h);
	bool try_compact_to_lum();
	void make_normal_map();
	void gen_rand_texture(unsigned char val, unsigned char a_add=0, unsigned a_rand=256);
	void load_from_gl();
	void deferred_load_and_bind();
	void update_texture_data(int x1, int y1, int x2, int y2);
	int write_to_jpg(string const &fn) const;
	int write_to_bmp(string const &fn) const;
	int write_to_png(string const &fn) const;
	unsigned get_texel_ix(float u, float v) const;
	// assumes width and height are a power of 2
	unsigned get_texel_ix_fast_pow2(float u, float v) const {return (width*(int(height*v) & (height-1)) + (int(width*u) & (width-1)));}
	colorRGBA get_texel(unsigned ix) const;
	colorRGBA get_texel(float u, float v) const {return get_texel(get_texel_ix(u, v));}
	float get_component(float u, float v, int comp) const;
	float get_component_grayscale_pow2(float u, float v) const {return data[get_texel_ix_fast_pow2(u, v)]/255.0;}
	void check_init(bool free_after_upload=0) {if (tid == 0) do_gl_init(free_after_upload);}
	unsigned num_pixels() const {return unsigned(width*height);}
	unsigned num_bytes()  const {return ncolors*num_pixels();}
	unsigned bytes_per_channel() const {return (is_16_bit_gray ? 2U : 1U);}
	unsigned get_cpu_mem() const {return (is_allocated() ? num_bytes() : 0);} // Note: ignores other data; excludes deferred load/DDS textures
	unsigned get_gpu_mem() const;
	unsigned get_tid() const {return tid;} // for passing into bind_texture_tu() calls
	void set_color_alpha_to_one() {color.alpha = 1.0;} // to make has_alpha() return 0
	bool has_alpha()    const {return (color.alpha < 1.0 || alpha_tid >= 0);}
	bool is_bound()     const {return (tid > 0);}
	bool is_allocated() const {return (data != nullptr);}
	bool defer_load()   const {return (defer_load_type != DEFER_TYPE_NONE);}
	bool is_loaded()    const {return (is_allocated() || defer_load());}
	colorRGBA get_avg_color() const {return color;}
	unsigned char *get_data() {assert(data); return data;}
	unsigned char const *get_data() const {assert(data); return data;}
	void write_pixel_16_bits(unsigned ix, float val);
};


struct camera_filter {

	int tid;
	unsigned time, init_time; // in ticks
	bool fades;
	colorRGBA color;

	camera_filter() : tid(-1), time(0), init_time(0), fades(0) {}
	camera_filter(colorRGBA const &c, unsigned t, int tid_, bool fades_) : tid(tid_), time(t), init_time(t), fades(fades_), color(c) {}
	void draw(bool apply_texture=1);
};


struct portal {

	point pts[4]; // quads only, for now
	vector3d normal; // for back face determination

	portal() : normal(zero_vector) {}
	static void pre_draw(vector<vert_wrap_t> &verts);
	static void post_draw(vector<vert_wrap_t> &verts);
	void draw(vector<vert_wrap_t> &verts) const;
	point get_center_pt() const {return (pts[0] + pts[1] + pts[2] + pts[3])*0.25;}
	bool is_visible(int reflection_pass) const;
};


struct fire_elem_t {

	float hp, fuel, burn_amt;

	fire_elem_t() : hp(0.0), fuel(0.0), burn_amt(0.0) {}
	bool burn(float val);
	void next_frame(float burn_rate, float consume_rate, float die_rate=1.0);
	static float get_burn_rate();
};


class shader_t;
class vpc_shader_t;

class volume_part_cloud {

public:
	typedef vert_norm_comp vert_type_t;

protected:
	static vector<vert_type_t> unscaled_points[2];
	vector<vert_type_t> points;

public:
	static colorRGBA gen_color(rand_gen_t &rgen);
	static void calc_unscaled_points(bool simplified);
	vector<vert_type_t> const &get_points() const {return points;}
	void gen_pts(vector3d const &size, point const &pos=all_zeros, bool simplified=0);
	void gen_pts(float radius, point const &pos=all_zeros, bool simplified=0) {gen_pts(vector3d(radius, radius, radius), pos, simplified);}
	static void shader_setup(vpc_shader_t &s, unsigned noise_ncomp, bool ridged=1, float alpha_bias=-0.4, float dist_bias=0.0,
		unsigned num_octaves=5, bool enable_lighting=0, bool use_cloud_mode=0, bool irregular_shape=0);
	void draw_quads(bool depth_map_already_disabled=0) const;
	bool enabled() const {return !points.empty();}
	void clear() {points.clear();}
};


struct water_params_t {

	float alpha, mud, algae, bright, reflect, green, wave_amp;
	water_params_t() {set_def_water();}
	void set_def_water();
	void set_def_lava();
};


struct text_string_t {

	string str;
	point pos;
	float size;
	colorRGBA color;

	text_string_t() : size(0.0), color(0,0,0,0) {}
	text_string_t(string const &s, point const &p, float sz, colorRGBA const &c) : str(s), pos(p), size(sz), color(c) {}
};

class popup_text_t : public text_string_t {
	float dist, time, tfticks_last_drawn;
	unsigned mode; // 0=one time, 1=on enter, 2=continuous
	bool any_active, prev_active;

public:
	popup_text_t() : dist(0.0), time(1.0), tfticks_last_drawn(0.0), mode(0), any_active(0), prev_active(0) {}
	bool read(FILE *fp, unsigned &line_num);
	void write(std::ostream &out) const;
	void check_player_prox();
	void draw() const;
};

struct text_drawer_t {
	vector<text_string_t> strs;
	void draw() const;
};

template<typename T> void clear_cont(T &cont) {T().swap(cont);}
template<typename T> unsigned get_cont_mem_usage(vector<T> const &v) {return v.capacity()*sizeof(T);}


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
colorRGBA const TRANS    (1.0,  1.0,  1.0,  0.4);
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
colorRGBA const PTREE_C  (0.8,  0.6,  0.5,  1.0);
colorRGBA const MUD_C    (0.3,  0.18, 0.09, 1.0);
colorRGBA const MUD_S_C  (0.7,  0.35, 0.25, 1.0);
colorRGBA const GOLD     (0.7,  0.45, 0.05, 1.0);
colorRGBA const BRASS_C  (0.7,  0.65, 0.25, 1.0);
colorRGBA const BRONZE_C (0.52, 0.23, 0.17, 1.0);
colorRGBA const COPPER_C (0.7,  0.3,  0.05,  1.0);
colorRGBA const MED_GREEN(0.2,  0.7,  0.2,  1.0);
colorRGBA const LT_YELLOW(1.0,  1.0,  0.1,  1.0);
colorRGBA const SNOW_COLOR(1.0, 1.0,  1.4,  1.0); // slightly bluish (saturated)
colorRGBA const LAVA_COLOR(1.0, 0.15, 0.05, 1.0); // red-orange
colorRGBA const FREEZE_COLOR(0.3, 0.3, 1.0, 1.0);
colorRGBA const BACKGROUND_DAY(0.2, 0.3, 0.8, 1.0);
colorRGBA const BACKGROUND_NIGHT(BLACK);


void register_timing_value(const char *str, int delta_time, bool no_loading_screen=0);
void toggle_timing_profiler();
void timing_profiler_stats();

// macros
#define GET_TIME_MS()    glutGet(GLUT_ELAPSED_TIME)
#define RESET_TIME       int const timer1(GET_TIME_MS());
#define GET_DELTA_TIME   (GET_TIME_MS() - timer1)
#define PRINT_TIME(str) {register_timing_value(str, GET_DELTA_TIME);}
#define PRINT_TIME_STR(str) {PRINT_TIME((str).c_str());}
#define PRINT_TIME_ONSCREEN(str) {std::ostringstream oss; oss << str << " time = " << GET_DELTA_TIME; print_debug_text(oss);}
#define TXT(x) #x"=" << x << " "
#define TXTn(x) #x"=" << x << "\n"
// these two are intended to be used with char/unsigned char where we want to print them as integers rather than characters
#define TXTi(x) #x"=" << (int)x << " "
#define TXTin(x) #x"=" << (int)x << "\n"

#if 0 // still only ms accuracy in windows
#include <time.h>
#define GET_TIME_CLOCK()  clock()
#define RESET_TIME2       clock_t const timer1(GET_TIME_CLOCK());
#define GET_DELTA_TIME2   double(GET_TIME_CLOCK() - timer1)/double(CLOCKS_PER_SEC)
#define PRINT_TIME2(str) {cout << str << " time = " << GET_DELTA_TIME << endl;}
#endif

class timer_t {
	string name;
	int timer1;
	bool enabled, no_loading_screen;
public:
	timer_t(char const *const  name_, bool enabled_=1, bool nls=0) : name(name_), timer1(GET_TIME_MS()), enabled(enabled_), no_loading_screen(nls) {}
	timer_t(string      const &name_, bool enabled_=1, bool nls=0) : name(name_), timer1(GET_TIME_MS()), enabled(enabled_), no_loading_screen(nls) {}
	~timer_t() {end();}
	void end() {if (enabled && !name.empty()) {register_timing_value(name.c_str(), GET_DELTA_TIME, no_loading_screen); name.clear();}}
};

struct status_bar_t {
	colorRGBA color;
	float val;
	unsigned icon_id;
	status_bar_t(colorRGBA const &c, float v, unsigned id=0) : color(c), val(v), icon_id(id) {}
};

// status bar icons
enum {ICON_HEALTH=0, ICON_SHIELD, ICON_POWER, ICON_DRUNK, ICON_TOILET, ICON_WATER, ICON_OXYGEN, ICON_CARRY, NUM_ICONS};


// world modes
enum {WMODE_GROUND=0, WMODE_UNIVERSE, WMODE_INF_TERRAIN, NUM_WMODE};

// movement directions
enum {MOVE_FRONT, MOVE_BACK, MOVE_LEFT, MOVE_RIGHT, MOVE_STOP};

// gameplay powerups
enum {PU_NONE=-1, PU_DAMAGE, PU_REGEN, PU_SHIELD, PU_SPEED, PU_FLIGHT, PU_INVISIBILITY, NUM_POWERUPS};

// object definitions
enum {RAIN = 0, SNOW, HAIL, LEAF, BALL, S_BALL, SMILEY, BLOOD, CHARRED, CHUNK,
	SFPART, ROCKET, LANDMINE, SEEK_D, STAR5, PLASMA, GRENADE, CGRENADE, SHRAPNEL, SHELLC,
	PROJC, DROPLET, WDROPLET, SAND, DIRT, ROCK, FRAGMENT, PARTICLE, HEALTH, SHIELD,
	POWERUP, WEAPON, AMMO, WA_PACK, CAMERA, PRECIP, BLAST_RADIUS, PROJECTILE, BEAM, IMPACT,
	PLASMA_LT_D, LASER, DROWNED, BURNED, FIRE, FELL, FROZEN, SUFFOCATED, CRUSHED, GASSED,
	WAYPOINT, SMOKE, DYNAM_PART, SKULL, GRASS, TELEFRAG, SAWBLADE, MAT_SPHERE, COLLISION, RAPT_PROJ,
	FREEZE_BOMB, XLOCATOR, XLOCATOR_DEATH, JUMP_PAD, TELEPORTER, KEYCARD, NUM_TOT_OBJS};

inline bool is_rocket_type(int type) {return (type == ROCKET || type == SEEK_D || type == RAPT_PROJ || type == FREEZE_BOMB);}

// textures
enum {GROUND_TEX = 0, DARK_ROCK_TEX, WATER_TEX, STUCCO_TEX, CLOUD_TEX, BRICK2_TEX, MOON_TEX, EARTH_TEX, MARBLE_TEX, SNOW_TEX,
	LEAF_TEX, WOOD_TEX, SAND_TEX, DIRT_TEX, CAMOFLAGE_TEX, HEDGE_TEX, BRICK_TEX, MANHOLE_TEX, PALM_FROND_TEX, SMOKE_TEX,
	PLASMA_TEX, GEN_TEX, LANDSCAPE_TEX, TREE_END_TEX, TREE_HEMI_TEX, SHINGLE_TEX, PANELING_TEX, CBLOCK_TEX, MJ_LEAF_TEX, LIVE_OAK_TEX,
	LEAF2_TEX, LEAF3_TEX, PLANT1_TEX, PLANT2_TEX, PLANT3_TEX, PLANT4_TEX, FENCE_TEX, SKULL_TEX, RADIATION_TEX, YUCK_TEX,
	SAW_TEX, SAW_B_TEX, BLUR_TEX, SBLUR_TEX, PINE_TEX, NOISE_TEX, WOOD2_TEX, HBBRICK_TEX, PARTB_TEX, PLASTER_TEX,
	TILE_TEX, LOGO_TEX, DISINT_TEX, BLUR_TEX_INV, HSTRIPE_TEX, VSTRIPE_TEX, BCUBE_TEX, EXPLOSION_TEX, SHIP_HULL_TEX, BCUBE2_TEX,
	BCUBE_T_TEX, ROCK_SPHERE_TEX, PAPAYA_TEX, COFFEE_TEX, SMILEY_SKULL_TEX, ICE_TEX, ROCK_TEX, BLACK_TEX, WHITE_TEX, FIRE_TEX,
	CLOUD_RAW_TEX, SNOWFLAKE_TEX, BLUR_CENT_TEX, GRADIENT_TEX, GRASS_BLADE_TEX, WIND_TEX, MOSSY_ROCK_TEX, BARK1_TEX, BARK2_TEX, BARK2_NORMAL_TEX,
	BARK3_TEX, BARK4_TEX, WATER_NORMAL_TEX, OCEAN_WATER_NORMAL_TEX, WATER_CAUSTIC_TEX, PS_NOISE_TEX, NOISE_GEN_TEX, NOISE_GEN_MIPMAP_TEX, SPARSE_NOISE_TEX, PLAYER_BBB_TEX,
	PINE_TREE_TEX, FLARE1_TEX, FLARE2_TEX, FLARE3_TEX, FLARE4_TEX, FLARE5_TEX, FOAM_TEX, BULLET_D_TEX, BULLET_A_TEX, BULLET_N_TEX,
	ROCK_NORMAL_TEX, RAINDROP_TEX, SPACESHIP1_TEX, SPACESHIP2_TEX, BLOOD_SPLAT_TEX, LICHEN_TEX, PALM_BARK_TEX, DAISY_TEX, LAVA_TEX, SMOKE_PUFF_TEX,
	BARK5_TEX, BARK6_TEX, RIPPLE_MAP_TEX, STARBURST_TEX, ROCK1_NORMAL_TEX, ROCK2_NORMAL_TEX, ROCK3_NORMAL_TEX, DIRT_NORMAL_TEX, FLAT_NMAP_TEX, RED_TEX,
	HAZARD_TEX, BLDG_WINDOW_TEX, BLDG_WIND_TRANS_TEX, KEYCARD_TEX, NUM_PREDEF_TEXTURES
};

// lighting files/types
enum {LIGHTING_SKY=0, LIGHTING_GLOBAL, LIGHTING_LOCAL, LIGHTING_COBJ_ACCUM, LIGHTING_DYNAMIC /*must be last*/, NUM_LIGHTING_TYPES};

// heightmap generation modes
enum {MGEN_SINE=0, MGEN_SIMPLEX, MGEN_PERLIN, MGEN_SIMPLEX_GPU, MGEN_DWARP_GPU, MGEN_END};


// shadow mask bits
#define MESH_SHADOW      0x02
#define SHADOWED_ALL     0xCF

// light sources
enum {LIGHT_SUN = 0, LIGHT_MOON, NUM_LIGHT_SRC};

#define SUN_SHADOW  0x01
#define MOON_SHADOW 0x02

#define EF_Z1  0x01
#define EF_Z2  0x02
#define EF_Y1  0x04
#define EF_Y2  0x08
#define EF_X1  0x10
#define EF_X2  0x20
#define EF_Z12 0x03
#define EF_Y12 0x0c
#define EF_X12 0x30
#define EF_ALL 0x3F

unsigned char const EFLAGS[3][2] = {{EF_X1, EF_X2}, {EF_Y1, EF_Y2}, {EF_Z1, EF_Z2}};

// forward references
struct base_mat_t;
class  sd_sphere_d;
struct dwobject;
class  obj_group;
struct cobj_params;
class  coll_obj;
class  coll_obj_group;
class  shape3d;
struct color_tid_vol;
class  vert_coll_detector;
struct cobj_query_callback;
struct user_waypt_t;
class  voxel_model;

