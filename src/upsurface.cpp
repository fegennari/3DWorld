// 3D World - upsurface class, used for planets, moons, asteroids, etc.
// by Frank Gennari
// 4/5/07

#include "function_registry.h"
#include "upsurface.h"
#include "universe.h"
#include "sinf.h"
#include "textures_3dw.h"


unsigned const ND_TEST     = 32;
unsigned const ND_CUBE     = 16;
float const M_ATTEN_FACTOR = 0.5;
float const F_ATTEN_FACTOR = 0.4;


extern int display_mode;
extern texture_t textures[];



void noise_gen_3d::gen_sines(float mag, float freq) {

	assert(SINES_PER_FREQ >= 2);
	assert(mag > 0.0 && freq > 0.0);

	for (unsigned i = 0; i < MAX_FREQ_BINS; ++i) { // low frequencies first
		unsigned const offset2(SINES_PER_FREQ*i);

		for (unsigned j = 0; j < SINES_PER_FREQ; ++j) {
			unsigned const offset(NUM_SINE_PARAMS*(offset2 + j));
			rdata[offset+0] = rgen.rand_uniform(0.2, 1.0)*mag;  // magnitude
			rdata[offset+1] = rgen.rand_uniform(0.1, 1.0)*freq; // x frequency
			rdata[offset+2] = rgen.randd()*TWO_PI; // x phase
			rdata[offset+3] = rgen.rand_uniform(0.1, 1.0)*freq; // y frequency
			rdata[offset+4] = rgen.randd()*TWO_PI; // y phase
			rdata[offset+5] = rgen.rand_uniform(0.1, 1.0)*freq; // z frequency
			rdata[offset+6] = rgen.randd()*TWO_PI; // z phase
		}
		mag  *= M_ATTEN_FACTOR;
		freq /= F_ATTEN_FACTOR;
	}
	num_sines = TOT_NUM_SINES; // initial value; may decrease later
}


void noise_gen_3d::gen_xyz_vals(point const &start, vector3d const &step, unsigned const xyz_num[3], vector<float> xyz_vals[3]) {

	for (unsigned d = 0; d < 3; ++d) {
		xyz_vals[d].resize(num_sines*xyz_num[d]);
		float val(start[d]);

		for (unsigned i = 0; i < xyz_num[d]; ++i) {
			for (unsigned k = 0; k < num_sines; ++k) {
				unsigned const index2(NUM_SINE_PARAMS*k + 2*d);
				xyz_vals[d][i*num_sines + k] = SINF(rdata[index2+1]*val + rdata[index2+2]);
			}
			val += step[d];
		}
	}
}


float noise_gen_3d::get_val(unsigned x, unsigned y, unsigned z, vector<float> const xyz_vals[3]) const {

	float val(0.0);
	unsigned const xyz[3] = {x, y, z};
	UNROLL_3X(assert(num_sines*xyz[i_]+num_sines <= xyz_vals[i_].size()););

	for (unsigned k = 0; k < num_sines; ++k) { // performance critical
		val += rdata[NUM_SINE_PARAMS*k]*(xyz_vals[0][x*num_sines + k])*(xyz_vals[1][y*num_sines + k])*(xyz_vals[2][z*num_sines + k]);
	}
	return val;
}


float noise_gen_3d::get_val(point const &pt) const {

	float val(0.0);

	for (unsigned k = 0; k < num_sines; ++k) { // performance critical
		unsigned const index2(NUM_SINE_PARAMS*k);
		float const x(SINF(rdata[index2+1]*pt.x + rdata[index2+2])); // faster sinf() calls
		float const y(SINF(rdata[index2+3]*pt.y + rdata[index2+4]));
		float const z(SINF(rdata[index2+5]*pt.z + rdata[index2+6]));
		val += rdata[index2]*x*y*z;
	}
	return val;
}


void upsurface::gen(float mag, float freq, unsigned ntests, float mm_scale) {

	max_mag = 0.0;
	pair<float, unsigned> hf_comps[SINES_PER_FREQ];
	gen_sines(mag, freq);
	
	for (unsigned i = 0; i < MAX_FREQ_BINS; ++i) { // low frequencies first
		unsigned const offset2(SINES_PER_FREQ*i);

		for (unsigned j = 0; j < SINES_PER_FREQ; ++j) {
			unsigned const offset(NUM_SINE_PARAMS*(offset2 + j));
			float const fmin(min(min(rdata[offset+1], rdata[offset+3]), rdata[offset+5]));
			hf_comps[j] = make_pair(fmin*rdata[offset+0], j);
		}
		sort(hf_comps, (hf_comps + SINES_PER_FREQ));
		float const largest(hf_comps[SINES_PER_FREQ-1].first), next_largest(hf_comps[SINES_PER_FREQ-2].first);
		
		if (largest > 1.5*next_largest) { // scale down the dominant high frequency
			rdata[NUM_SINE_PARAMS*(offset2 + hf_comps[SINES_PER_FREQ-1].second)+0] *= 1.5*next_largest/largest;
		}
	}
	for (unsigned i = 0; i < ntests; ++i) { // choose random test points to approximate max value (procedural way to do this?)
		float val(0.0);

		for (unsigned j = 0; j < TOT_NUM_SINES; ++j) {
			unsigned const offset(NUM_SINE_PARAMS*j);
			float localval(rdata[offset]);
			UNROLL_3X(localval *= SINF(rdata[offset+(i_<<1)+1]*rgen.randd() + rdata[offset+(i_<<1)+2]);)
			val += fabs(localval);
		}
		max_mag = max(max_mag, val);
	}
	max_mag /= mm_scale;
	val_cache.clear();
}


void upsurface::setup(unsigned size, float mcut, bool alloc_hmap) {

	ssize      = size;
	min_cutoff = mcut;
	if (alloc_hmap) heightmap.resize(ssize*ssize);
	unsigned max_freq(MAX_FREQ_BINS - 4);

	for (unsigned i = 8; i <= MAX_TEXTURE_SIZE; i <<= 1) {
		if (ssize <= i) break;
		++max_freq;
	}
	max_freq  = max(1u, min(MAX_FREQ_BINS, max_freq));
	num_sines = max_freq*SINES_PER_FREQ;
}


float upsurface::get_height_at(point const &pt, bool use_cache) const {

	unsigned const CACHE_SIZE(1<<18); // 4MB (only allocated if use_cache is enabled) (or use prime number 104729)
	size_t cache_index(0);
	cache_entry ce(pt);

	if (use_cache) {
		if (val_cache.empty()) {val_cache.resize(CACHE_SIZE);}
		assert(val_cache.size() == CACHE_SIZE);
		cache_index = ce.hash()%CACHE_SIZE;
		if (val_cache[cache_index].p == pt) {return val_cache[cache_index].val;}
	}
	float val(get_val(pt)); // performance critical
	val = 0.5*(max(-1.0f, min(1.0f, (1.5f/max_mag)*val)) + 1.0); // duplicate code
	
	if (use_cache) {
		ce.val = val;
		val_cache[cache_index] = ce;
	}
	return val;
}


void upsurface::setup_draw_sphere(point const &pos, float radius, float dp, int ndiv, float const *const pmap) {

	sd.set_data(pos, radius, ndiv, pmap, dp, (pmap ? NULL : this));
	sd.gen_points_norms(spn);
}


upsurface::~upsurface() {

	spn.free_data();
	free_context();
}


void urev_body::gen_surface() {

	set_rseeds();
	surface.reset(new upsurface(type)); // may delete a previous surface
	float mag(SURFACE_HEIGHT*radius), freq(((type == UTYPE_MOON) ? 1.5 : 1.0)*INITIAL_FREQ*TWO_PI);
	surface->set_rand_seeds(urseed1, urseed2);
	surface->gen(mag, freq);
}


// Note: many planet/sphere renderers use a texture with width = 2*height, which yields square regions at the equator
// here we use a square texture for simplicity, so that this code can be shared with (and be similar to)
// the rest of the 3DWorld sphere generation and drawing code; it also produces more uniform regions near the poles
void urev_body::gen_texture_data_and_heightmap(unsigned char *data, unsigned size) {

	//RESET_TIME;
	get_colors(a, b);
	unsigned size_p2(0);
	for (unsigned sz = size; sz > 1; sz >>= 1, ++size_p2);
	assert((1U<<size_p2) == size); // size must be a power of 2
	assert(surface != nullptr);
	unsigned const table_size(MAX_TEXTURE_SIZE << 1); // larger is more accurate
	static float xtable[TOT_NUM_SINES*table_size], ytable[TOT_NUM_SINES*table_size];
	surface->setup(size, max(water, lava), 1); // use_heightmap=1
	unsigned const num_sines(surface->num_sines);
	float const *const rdata(surface->rdata);
	float const sizef(size/TWO_PI), mt2(0.5*(table_size-1)), scale(1.5/surface->max_mag);
	float const delta(TWO_PI/size), sin_ds(sin(delta)), cos_ds(cos(delta));
	unsigned const pole_thresh(size>>3);
	wr_scale = 1.0/max(0.01, (1.0 - water));

	for (unsigned i = 0; i < table_size; ++i) { // build sin table
		unsigned const offset(i*num_sines);
		float const sarg(i/mt2 - 1.0);

		for (unsigned k = 0; k < num_sines; ++k) { // create x and y tables
			unsigned const index2(NUM_SINE_PARAMS*k);
			xtable[offset+k] = SINF(rdata[index2+1]*sarg + rdata[index2+2]);
			ytable[offset+k] = SINF(rdata[index2+3]*sarg + rdata[index2+4]);
		}
	}

	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < (int)size; ++i) { // phi values
		unsigned const hmoff(i*size), ti(size-i-1), texoff(ti*size);
		float const phi((float(i)/(size-1))*PI);
		float const sin_phi((i == size-1) ? 0.0 : sinf(phi)), zval((i == size-1) ? -1.0 : cosf(phi));
		float sin_s(0.0), cos_s(1.0);
		float ztable[TOT_NUM_SINES];

		for (unsigned k = 0; k < num_sines; ++k) { // create z table
			unsigned const index2(NUM_SINE_PARAMS*k);
			ztable[k] = rdata[index2]*SINF(rdata[index2+5]*zval + rdata[index2+6]);
		}
		for (unsigned j = 0; j < size; ++j) { // theta values, Note: x and y are swapped because theta is out of phase by 90 degrees to match tex coords
			float const s(sin_s), c(cos_s), xval(sin_phi*s), yval(sin_phi*c);
			unsigned const tj(size-j-1), index(3*(texoff + tj));
			unsigned const ox1((unsigned((xval+1.0)*mt2))*num_sines), oy1((unsigned((yval+1.0)*mt2))*num_sines);
			float val(0.0);

			if (i <= (int)pole_thresh || i >= int(size-pole_thresh-1)) { // slower version near the poles
				for (unsigned k = 0; k < num_sines; ++k) {
					unsigned const index2(NUM_SINE_PARAMS*k);
					val += ztable[k]*SINF(rdata[index2+1]*xval + rdata[index2+2])*SINF(rdata[index2+3]*yval + rdata[index2+4]);
				}
			}
			else {
				// Note: chooses the closest precomputed grid point for efficiency -
				// no interpolation, so has artifacts closer to the poles
				for (unsigned k = 0; k < num_sines; ++k) {val += ztable[k]*xtable[ox1+k]*ytable[oy1+k];}
			}
			val = 0.5*(max(-1.0f, min(1.0f, scale*val)) + 1.0);
			surface->heightmap[hmoff + j] = val;
			get_surface_color((data + index), val, phi);
			sin_s = s*cos_ds + c*sin_ds;
			cos_s = c*cos_ds - s*sin_ds;
		} // for j
	} // for i
	//if (size >= MAX_TEXTURE_SIZE) PRINT_TIME("Gen");
}


bool urev_body::surface_test(float rad, point const &p, float &coll_r, bool simple) const {

	// not quite right - should take into consideration peaks in surrounding geometry that also intersect the sphere
	if (has_heightmap()) {
		if (!dist_less_than(p, pos, (radius*(1.0 + get_hmap_scale()*0.5) + rad))) return 0; // test rmax
		coll_r = get_radius_at(p, !simple);
	}
	return 1;
}


float urev_body::get_radius_at(point const &p, bool exact) const {

	if (surface == NULL) return radius;
	double const val(get_dheight_at(p, exact)), cutoff(surface->min_cutoff);
	return radius*(1.0 + get_hmap_scale()*((max(cutoff, val) - cutoff)*surface->get_one_minus_cutoff() - 0.5));
}


float urev_body::get_dheight_at(point const &p, bool exact) const {

	if (!has_heightmap()) return 0.0;
	point coll_from(p - pos);
	rotate_vector(coll_from); // in radians
	coll_from.normalize();
	if (exact) {return surface->get_height_at(coll_from, 0);} // slower but more accurate

	unsigned const tsize(surface->ssize);
	unsigned const tx(unsigned(tsize*((atan2(coll_from.y, coll_from.x)/TWO_PI) + 0.5) + 0.25*tsize)%tsize);
	unsigned const ty(unsigned(tsize*(1.0-safe_acosf(coll_from.z)/PI)));
	assert(tx < tsize && ty < tsize);
	return surface->heightmap[(tsize-tx-1) + tsize*(tsize-ty-1)];
}


bool urev_body::pt_over_land(point const &p) const {

	if (water < 0.1 || !has_heightmap()) return 1;
	return (get_dheight_at(p, 1) > surface->min_cutoff);
}

