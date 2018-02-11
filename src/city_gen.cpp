// 3D World - City Generation
// by Frank Gennari
// 2/10/18

#include "3DWorld.h"
#include "mesh.h"
#include "heightmap.h"
#include "file_utils.h"

using std::string;

bool const CHECK_HEIGHT_BORDER_ONLY = 1; // choose building site to minimize edge discontinuity rather than amount of land that needs to be modified


extern int rand_gen_index, display_mode;
extern float water_plane_z;


class city_plot_gen_t {

protected:
	struct rect_t {
		unsigned x1, y1, x2, y2;
		rect_t() : x1(0), y1(0), x2(0), y2(0) {}
		rect_t(unsigned x1_, unsigned y1_, unsigned x2_, unsigned y2_) : x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
		bool is_valid() const {return (x1 < x2 && y1 < y2);}
		unsigned get_area() const {return (x2 - x1)*(y2 - y1);}
		bool operator== (rect_t const &r) const {return (x1 == r.x1 && y1 == r.y1 && x2 == r.x2 && y2 == r.y2);}
		bool has_overlap(rect_t const &r) const {return (x1 < r.x2 && y1 < r.y2 && r.x1 < x2 && r.y1 < y2);}
	};

	float *heightmap;
	unsigned xsize, ysize;
	int last_rgi;
	rand_gen_t rgen;
	vector<rect_t> used;

	bool is_valid_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		return (x1 < x2 && y1 < y2 && x2 <= xsize && y2 <= ysize);
	}
	bool overlaps_used(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		rect_t const cur(x1, y1, x2, y2);
		for (vector<rect_t>::const_iterator i = used.begin(); i != used.end(); ++i) {if (i->has_overlap(cur)) return 1;} // simple linear iteration
		return 0;
	}
	void mark_used(unsigned x1, unsigned y1, unsigned x2, unsigned y2) {used.emplace_back(x1, y1, x2, y2);}

	float any_underwater(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				if (heightmap[y*xsize + x] < water_plane_z) return 1;
			}
		}
		return 0;
	}
	float get_avg_height(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		assert(is_valid_region(x1, y1, x2, y2));
		float sum(0.0), denom(0.0);

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				sum   += heightmap[y*xsize + x];
				denom += 1.0;
			}
		}
		return sum/denom;
	}
	float get_rms_height_diff(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		float const avg(get_avg_height(x1, y1, x2, y2));
		float diff(0.0);

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				float const delta(heightmap[y*xsize + x] - avg);
				diff += delta*delta; // square the difference
			}
		}
		return diff;
	}
public:
	city_plot_gen_t() : heightmap(nullptr), xsize(0), ysize(0), last_rgi(0) {}

	void init(float *heightmap_, unsigned xsize_, unsigned ysize_) {
		heightmap = heightmap_; xsize = xsize_; ysize = ysize_;
		assert(heightmap != nullptr);
		assert(xsize > 0 && ysize > 0); // any size is okay
		if (rand_gen_index != last_rgi) {rgen.set_state(rand_gen_index, 12345); last_rgi = rand_gen_index;} // only when rand_gen_index changes
	}
	bool find_best_city_location(unsigned width, unsigned height, unsigned border, unsigned num_samples, unsigned &x_llc, unsigned &y_llc) {
		cout << TXT(xsize) << TXT(ysize) << TXT(width) << TXT(height) << TXT(border) << endl;
		assert(num_samples > 0);
		assert((width + 2*border) < xsize && (height + 2*border) < ysize); // otherwise the city can't fit in the map
		unsigned const num_iters(100*num_samples); // upper bound
		unsigned xend(xsize - width - 2*border + 1), yend(ysize - width - 2*border + 1); // max rect LLC, inclusive
		unsigned num_cands(0);
		float best_diff(0.0);

		for (unsigned n = 0; n < num_iters; ++n) { // find min RMS height change across N samples
			unsigned const x1(border + (rgen.rand()%xend)), y1(border + (rgen.rand()%yend));
			unsigned const x2(x1 + width), y2(y1 + height);
			if (overlaps_used (x1, y1, x2, y2)) continue; // skip
			if (any_underwater(x1, y1, x2, y2)) continue; // skip
			float const diff(get_rms_height_diff(x1, y1, x2, y2));
			if (num_cands == 0 || diff < best_diff) {x_llc = x1; y_llc = y1; best_diff = diff;}
			if (++num_cands == num_samples) break; // done
		} // for n
		if (num_cands == 0) return 0;
		cout << "cands: " << num_cands << ", diff: " << best_diff << ", loc: " << x_llc << "," << y_llc << endl;
		mark_used(x_llc, y_llc, x_llc+width, y_llc+height);
		return 1; // success
	}
	void flatten_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float const *const height=nullptr) {
		assert(is_valid_region(x1, y1, x2, y2));
		float const delta_h = 0.0; // for debugging in map view
		float const h(height ? *height : (get_avg_height(x1, y1, x2, y2) + delta_h));

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {heightmap[y*xsize + x] = h;}
		}
	}
};


struct city_params_t {

	unsigned num_cities, num_samples, city_size, city_border;

	city_params_t() : num_cities(0), num_samples(100), city_size(0), city_border(0) {}
	bool enabled() const {return (num_cities > 0 && city_size > 0);}
	static bool read_error(string const &str) {cout << "Error reading city config option " << str << "." << endl; return 0;}

	bool read_option(FILE *fp) {
		char strc[MAX_CHARS] = {0};
		if (!read_str(fp, strc)) return 0;
		string const str(strc);

		if (str == "num_cities") {
			if (!read_uint(fp, num_cities)) {return read_error(str);}
		}
		else if (str == "num_samples") {
			if (!read_uint(fp, num_samples) || num_samples == 0) {return read_error(str);}
		}
		else if (str == "city_size") {
			if (!read_uint(fp, city_size)) {return read_error(str);}
		}
		else if (str == "city_border") {
			if (!read_uint(fp, city_border)) {return read_error(str);}
		}
		else {
			cout << "Unrecognized city keyword in input file: " << str << endl;
			return 0;
		}
		return 1;
	}
};

class city_gen_t : public city_plot_gen_t {
	bool gen_city(unsigned city_size, unsigned border, unsigned num_samples) {
		timer_t t("Choose City Location");
		unsigned x1(0), y1(0);
		if (!find_best_city_location(city_size, city_size, border, num_samples, x1, y1)) return 0;
		unsigned const x2(x1 + city_size), y2(y1 + city_size);
		flatten_region(x1, y1, x2, y2);
		cube_t pos_range(all_zeros);
		int const dx(-int(xsize)/2), dy(-int(ysize)/2); // convert from center to LLC
		pos_range.d[0][0] = get_xval(x1 + dx);
		pos_range.d[0][1] = get_xval(x2 + dx);
		pos_range.d[1][0] = get_yval(y1 + dy);
		pos_range.d[1][1] = get_yval(y2 + dy);
		set_buildings_pos_range(pos_range);
		return 1;
	}
public:
	void gen_cities(city_params_t const &params) {
		for (unsigned n = 0; n < params.num_cities; ++n) {gen_city(params.city_size, params.city_border, params.num_samples);}
	}
};


city_params_t city_params;
city_gen_t city_gen;


bool parse_city_option(FILE *fp) {return city_params.read_option(fp);}
bool have_cities() {return city_params.enabled();}

void gen_cities(float *heightmap, unsigned xsize, unsigned ysize) {
	if (!have_cities()) return; // nothing to do
	city_gen.init(heightmap, xsize, ysize); // only need to call once for any given heightmap
	city_gen.gen_cities(city_params);
}

