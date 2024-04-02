// 3D World - Random Number Generators
// by Frank Gennari
// 6/23/20
#pragma once

#include <vector>
#include <memory>
#include <assert.h>

int const N_RAND_DIST = 10000;

template<typename T> struct pointT;
typedef pointT<float> point;
typedef pointT<float> vector3d;
struct cube_t;

extern float gauss_rand_arr[];

class rgen_core_t {
protected:
	// this is a good random number generator written by Stephen E. Derenzo
	template<typename T> inline void randome_int(T &ranptr) {
		if ((rseed1 = 40014*(rseed1%53668) - 12211*(rseed1/53668)) < 0) rseed1 += 2147483563;
		if ((rseed2 = 40692*(rseed2%52774) - 3791 *(rseed2/52774)) < 0) rseed2 += 2147483399;
		if ((ranptr = (T)rseed1 - (T)rseed2) < 1) ranptr += 2147483562;
	}

public:
	long rseed1, rseed2;

	rgen_core_t() {set_state(1,1);}
	void set_state(long rs1, long rs2) {rseed1 = rs1; rseed2 = rs2;}
	double randd();
};

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
typedef struct { uint64_t state; uint64_t inc; } pcg32_random_t;

inline uint32_t pcg32_random_r(pcg32_random_t* rng) {
	uint64_t oldstate = rng->state;
	// Advance internal state
	rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
	// Calculate output function (XSH RR), uses old state for max ILP
	uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
	uint32_t rot = oldstate >> 59u;
	return (xorshifted >> rot) | (xorshifted << ((-(int32_t)rot) & 31));
}

class rgen_pregen_t : public rgen_core_t {

	std::shared_ptr<std::vector<double>> pregen_rand_reals;
	unsigned cur_pos;

public:
	rgen_pregen_t() : cur_pos(0) {}
	void pregen_floats(unsigned num);
	double randd();
};

template<typename base> class rand_gen_template_t : public base {

public:
	using rgen_core_t::rseed1;
	using rgen_core_t::rseed2;

	int rand() {
		int rand_num;
		base::randome_int(rand_num);
		return rand_num;
	}
	int rand_fast() { // faster but lower quality, only uses rseed1
		rseed1 *= 16807;
		return rseed1;
	}
	float rand_float_fast() { // faster but lower quality, only uses rseed1 (uniform 0 to 1)
		union {float fres; unsigned int ires;};
		rseed1 *= 16807;
		ires = ((((unsigned long)rseed1)>>9) | 0x3f800000);
		return fres - 1.0f;
	}
	int rand_seed_mix() {
	  int val1(rand()); std::swap(rseed1, rseed2); return (val1 + rand()); // more random
		//return (rseed1 ^ (rseed2 >> 8)); // faster (should call rand2_mix() after)
	}
	void rand_mix() {rand(); std::swap(rseed1, rseed2);}
	float rand_float() {return 0.000001*(rand()%1000000);} // uniform 0 to 1
	float signed_rand_float() {return 2.0*float(base::randd()) - 1.0;}
	bool rand_bool() {return ((rand()&1) != 0);}
	float rand_uniform(float val1, float val2) {assert(val1 <= val2); return val1 + (val2 - val1)*float(base::randd());}
	unsigned rand_uniform_uint(unsigned min_val, unsigned max_val) {assert(min_val <= max_val); return (min_val + (rand() % (max_val - min_val + 1)));}
	int rand_int(int start, int end) {return (rand()%(end - start + 1) + start);} // used for trees; start and end should be positive
	float rgauss() {return gauss_rand_arr[rand()%N_RAND_DIST];} // mean = 0.0, std_dev = 1.0
	float rand_gaussian(float mean, float std_dev) {return mean + std_dev*rgauss();}
	bool rand_probability(float prob) {return (prob >= 1.0 || (prob > 0.0 && rand_float() < prob));}
	vector3d rand_vector(float scale=1.0);
	vector3d signed_rand_vector(float scale=1.0);
	vector3d signed_rand_vector_xy(float scale=1.0);
	vector3d signed_rand_vector_norm(float scale=1.0);
	vector3d signed_rand_vector_spherical(float scale=1.0);
	vector3d signed_rand_vector_spherical_xy(float scale=1.0);
	vector3d signed_rand_vector_spherical_noloop(float scale=1.0);
	vector3d signed_rand_vector_spherical_xy_norm();
	point gen_rand_cube_point(cube_t const &c);
	point gen_rand_cube_point_xy(cube_t const &c, float z=0.0);
};

typedef rand_gen_template_t<rgen_core_t> rand_gen_t;
typedef rand_gen_template_t<rgen_pregen_t> rand_gen_pregen_t;

// based on xxHash
// see: https://blogs.unity3d.com/2015/01/07/a-primer-on-repeatable-random-numbers/
inline unsigned RotateLeft(unsigned value, unsigned count) {return (value << count) | (value >> (32 - count));}

inline unsigned xxHash_uint(unsigned buf, unsigned seed=0) {
	unsigned h32 = seed + 374761393U;
	h32 += 4U;
	h32 += buf * 3266489917U;
	h32 = RotateLeft(h32, 17) * 668265263U;
	h32 ^= h32 >> 15;
	h32 *= 2246822519U;
	h32 ^= h32 >> 13;
	h32 *= 3266489917U;
	h32 ^= h32 >> 16;
	return h32;
}
