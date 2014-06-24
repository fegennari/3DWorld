// 3D World - upsurface class, used for planets, moons, asteroids, etc.
// by Frank Gennari
// 4/5/07

#ifndef _UPSURFACE_H_
#define _UPSURFACE_H_

#include "3DWorld.h"
#include "subdiv.h"


unsigned const MAX_TEXTURE_SIZE  = 256; // must be a power of 2
unsigned const SINES_PER_FREQ    = 12;
unsigned const MAX_FREQ_BINS     = 5;
unsigned const N_RAND_MAG_TESTS  = 100;

unsigned const TOT_NUM_SINES   = SINES_PER_FREQ*MAX_FREQ_BINS;
unsigned const NUM_SINE_PARAMS = 2*3+1; // 2*NUM_DIMENSIONS+1
unsigned const SINE_DATA_SIZE  = NUM_SINE_PARAMS*TOT_NUM_SINES;
unsigned const SUBDIV_SECTS    = 8;


class color_gen_class {

public:
	virtual void get_surface_color(unsigned char *data, float val, float phi) const = 0;
};


class ref_counted_obj {

	unsigned ref_count;

public:
	ref_counted_obj() : ref_count(0) {} // not 1
	void inc_ref() {++ref_count;}
	void dec_ref() {assert(ref_count > 0); --ref_count;}
	bool unrefed() const {return (ref_count == 0);}
};


class noise_gen_3d {
public:
	unsigned num_sines;
	float rdata[SINE_DATA_SIZE];
	rand_gen_t rgen;

	noise_gen_3d() : num_sines(0) {}
	void set_rand_seeds(int rs1, int rs2) {rgen.set_state(rs1, rs2);}
	void gen_sines(float mag, float freq);
	void gen_xyz_vals(point const &start, vector3d const &step, unsigned const xyz_num[3], vector<float> xyz_vals[3]);
	float get_val(unsigned x, unsigned y, unsigned z, vector<float> const xyz_vals[3]) const;
	float get_val(point const &pt) const;
};


class upsurface : public ref_counted_obj, public noise_gen_3d { // size = 104 + 4*336 = 1784 (+cache)

private:
	struct cache_entry {
		point p;
		float val;

		cache_entry(point const &pt=all_zeros) : p(pt), val(0.0) {}
		size_t hash() const {return (10831*(*(int *)(&p[0])) + 15601*(*(int *)(&p[1])) + 21401*(*(int *)(&p[2])));}
	};

	mutable vector<cache_entry> val_cache;

public:
	int type;
	unsigned ssize;
	float max_mag, rmax, min_cutoff;
	vector<float> heightmap;
	sphere_point_norm spn;
	sd_sphere_vbo_d sd;

	upsurface(int type_=0) : type(type_), ssize(0), max_mag(0.0), rmax(0.0), min_cutoff(0.0) {}
	~upsurface();
	void gen(float mag, float freq, unsigned ntests=N_RAND_MAG_TESTS, float mm_scale=1.0);
	void setup(unsigned size, float mcut, bool alloc_hmap);
	float get_one_minus_cutoff() const {return 1.0/max(0.01, (1.0 - min_cutoff));} // avoid div-by-zero
	float get_height_at(point const &pt, bool use_cache=0) const;
	void setup_draw_sphere(point const &pos, float radius, float dp, int ndiv, float const *const pmap);
	void calc_rmax() {rmax = sd.get_rmax();}
	void free_context() {sd.clear_vbos();}
	void clear_cache() {vector<cache_entry>().swap(val_cache);}
	bool has_heightmap() const {return (!heightmap.empty());}
};

typedef std::shared_ptr<upsurface> p_upsurface;


#endif // _UPSURFACE_H_


