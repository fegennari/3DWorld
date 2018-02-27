// 3D World - City Generation
// by Frank Gennari
// 2/10/18

#include "3DWorld.h"
#include "mesh.h"
#include "heightmap.h"
#include "file_utils.h"
#include "draw_utils.h"
#include "shaders.h"

using std::string;

bool const CHECK_HEIGHT_BORDER_ONLY = 1; // choose building site to minimize edge discontinuity rather than amount of land that needs to be modified
float const ROAD_HEIGHT             = 0.002;
float const OUTSIDE_TERRAIN_HEIGHT  = 0.0;
colorRGBA const road_color          = WHITE; // all road parts are the same color, to make the textures match

enum {TID_SIDEWLAK=0, TID_STRAIGHT, TID_BEND_90, TID_3WAY,   TID_4WAY,   NUM_RD_TIDS };
enum {TYPE_PLOT   =0, TYPE_RSEG,    TYPE_ISEC2,  TYPE_ISEC3, TYPE_ISEC4, NUM_RD_TYPES};


extern int rand_gen_index, display_mode;
extern float water_plane_z, shadow_map_pcf_offset, cobj_z_bias, fticks;


struct city_params_t {

	unsigned num_cities, num_samples, num_conn_tries, city_size_min, city_size_max, city_border, road_border, slope_width;
	float road_width, road_spacing, conn_road_seg_len, max_road_slope;
	// cars
	unsigned num_cars;
	float car_speed;

	city_params_t() : num_cities(0), num_samples(100), num_conn_tries(50), city_size_min(0), city_size_max(0), city_border(0), road_border(0),
		slope_width(0), road_width(0.0), road_spacing(0.0), conn_road_seg_len(1000.0), max_road_slope(1.0), num_cars(0), car_speed(0.0) {}
	bool enabled() const {return (num_cities > 0 && city_size_min > 0);}
	bool roads_enabled() const {return (road_width > 0.0 && road_spacing > 0.0);}
	float get_road_ar() const {return nearbyint(road_spacing/road_width);} // round to nearest texture multiple
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
		else if (str == "num_conn_tries") {
			if (!read_uint(fp, num_conn_tries) || num_conn_tries == 0) {return read_error(str);}
		}
		else if (str == "city_size_min") {
			if (!read_uint(fp, city_size_min)) {return read_error(str);}
			if (city_size_max == 0) {city_size_max = city_size_min;}
			if (city_size_max < city_size_min) {return read_error(str);}
		}
		else if (str == "city_size_max") {
			if (!read_uint(fp, city_size_max)) {return read_error(str);}
			if (city_size_min == 0) {city_size_min = city_size_max;}
			if (city_size_max < city_size_min) {return read_error(str);}
		}
		else if (str == "city_border") {
			if (!read_uint(fp, city_border)) {return read_error(str);}
		}
		else if (str == "road_border") {
			if (!read_uint(fp, road_border)) {return read_error(str);}
		}
		else if (str == "slope_width") {
			if (!read_uint(fp, slope_width)) {return read_error(str);}
		}
		else if (str == "road_width") {
			if (!read_float(fp, road_width) || road_width < 0.0) {return read_error(str);}
		}
		else if (str == "road_spacing") {
			if (!read_float(fp, road_spacing) || road_spacing < 0.0) {return read_error(str);}
		}
		else if (str == "conn_road_seg_len") {
			if (!read_float(fp, conn_road_seg_len) || conn_road_seg_len <= 0.0) {return read_error(str);}
		}
		else if (str == "max_road_slope") {
			if (!read_float(fp, max_road_slope) || max_road_slope <= 0.0) {return read_error(str);}
		}
		else if (str == "num_cars") {
			if (!read_uint(fp, num_cars)) {return read_error(str);}
		}
		else if (str == "car_speed") {
			if (!read_float(fp, car_speed) || car_speed < 0.0) {return read_error(str);}
		}
		else {
			cout << "Unrecognized city keyword in input file: " << str << endl;
			return 0;
		}
		return 1;
	}
}; // city_params_t

city_params_t city_params;


vector3d const get_query_xlate() {
	return vector3d((world_mode == WMODE_INF_TERRAIN) ? vector3d((xoff - xoff2)*DX_VAL, (yoff - yoff2)*DY_VAL, 0.0) : zero_vector);
}

bool check_bcube_sphere_coll(cube_t const &bcube, point const &sc, float radius, bool xy_only) {
	return (xy_only ? sphere_cube_intersect_xy(sc, radius, bcube) : sphere_cube_intersect(sc, radius, bcube));
}
template<typename T> bool check_bcubes_sphere_coll(vector<T> const &bcubes, point const &sc, float radius, bool xy_only) {
	for (auto i = bcubes.begin(); i != bcubes.end(); ++i) {
		if (check_bcube_sphere_coll(*i, sc, radius, xy_only)) return 1;
	}
	return 0;
}

template<typename T> void get_all_bcubes(vector<T> const &v, vector<cube_t> &bcubes) {
	for (auto i = v.begin(); i != v.end(); ++i) {bcubes.push_back(*i);}
}

template<typename T> static void add_flat_road_quad(T const &r, quad_batch_draw &qbd, float ar) { // z1 == z2
	float const z(r.z1()); // z1
	point const pts[4] = {point(r.x1(), r.y1(), z), point(r.x2(), r.y1(), z), point(r.x2(), r.y2(), z), point(r.x1(), r.y2(), z)};
	qbd.add_quad_pts(pts, road_color, plus_z, r.get_tex_range(ar));
}

float smooth_interp(float a, float b, float mix) {
	mix = mix * mix * (3.0 - 2.0 * mix); // cubic Hermite interoplation (smoothstep)
	return mix*a + (1.0 - mix)*b;
}

struct rect_t {
	unsigned x1, y1, x2, y2;
	rect_t() : x1(0), y1(0), x2(0), y2(0) {}
	rect_t(unsigned x1_, unsigned y1_, unsigned x2_, unsigned y2_) : x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
	bool is_valid() const {return (x1 < x2 && y1 < y2);}
	unsigned get_area() const {return (x2 - x1)*(y2 - y1);}
	bool operator== (rect_t const &r) const {return (x1 == r.x1 && y1 == r.y1 && x2 == r.x2 && y2 == r.y2);}
	bool has_overlap(rect_t const &r) const {return (x1 < r.x2 && y1 < r.y2 && r.x1 < x2 && r.y1 < y2);}
};
struct flatten_op_t : public rect_t {
	float z1, z2;
	bool dim;
	unsigned border;
	flatten_op_t() : z1(0.0), z2(0.0), dim(0), border(0) {}
	flatten_op_t(unsigned x1_, unsigned y1_, unsigned x2_, unsigned y2_, float z1_, float z2_, bool dim_, unsigned border_) :
		rect_t(x1_, y1_, x2_, y2_), z1(z1_), z2(z2_), dim(dim_), border(border_) {}
};

struct road_t : public cube_t {
	bool dim; // dim the road runs in
	bool slope; // 0: z1 applies to first (lower) point; 1: z1 applies to second (upper) point

	road_t() : dim(0), slope(0) {}
	road_t(cube_t const &c, bool dim_, bool slope_=0) : cube_t(c), dim(dim_), slope(slope_) {}
	road_t(point const &s, point const &e, float width, bool dim_, bool slope_=0) : dim(dim_), slope(slope_) {
		assert(s != e);
		assert(width > 0.0);
		vector3d const dw(0.5*width*cross_product((e - s), plus_z).get_norm());
		point const pts[4] = {(s - dw), (s + dw), (e + dw), (e - dw)};
		set_from_points(pts, 4);
	}
	float get_length   () const {return (d[dim][1] - d[dim][0]);}
	float get_slope_val() const {return get_dz()/get_length();}
};
struct road_seg_t : public road_t {
	road_seg_t(road_t const &r) : road_t(r) {}
	road_seg_t(cube_t const &c, bool dim_, bool slope_=0) : road_t(c, dim_, slope_) {}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, -ar, (dim ? -1.0 : 1.0), 0, dim);}

	void add_road_quad(quad_batch_draw &qbd, float ar) const { // specialized here for sloped roads
		if (z1() == z2()) {add_flat_road_quad(*this, qbd, ar); return;}
		bool const s(slope ^ dim);
		point pts[4] = {point(x1(), y1(), d[2][!s]), point(x2(), y1(), d[2][!s]), point(x2(), y2(), d[2][ s]), point(x1(), y2(), d[2][ s])};
		if (!dim) {swap(pts[0].z, pts[2].z);}
		vector3d const normal(cross_product((pts[2] - pts[1]), (pts[0] - pts[1])).get_norm());
		qbd.add_quad_pts(pts, road_color, normal, get_tex_range(ar));
	}
};
struct road_isec_t : public cube_t {
	uint8_t conn; // connected roads in {-x, +x, -y, +y}
	road_isec_t() : conn(15) {}
	road_isec_t(cube_t const &c, uint8_t conn_=15) : cube_t(c), conn(conn_) {}

	tex_range_t get_tex_range(float ar) const {
		switch (conn) {
		case 5 : return tex_range_t(0.0, 0.0, -1.0,  1.0, 0, 0); // 2-way: MX
		case 6 : return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 2-way: R0
		case 9 : return tex_range_t(0.0, 0.0, -1.0, -1.0, 0, 0); // 2-way: MXMY
		case 10: return tex_range_t(0.0, 0.0,  1.0, -1.0, 0, 0); // 2-way: MY
		case 7 : return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 3-way: R0
		case 11: return tex_range_t(0.0, 0.0, -1.0, -1.0, 0, 0); // 3-way: MY
		case 13: return tex_range_t(0.0, 0.0,  1.0, -1.0, 0, 1); // 3-way: R90MY
		case 14: return tex_range_t(0.0, 0.0, -1.0,  1.0, 0, 1); // 3-way: R90MX
		case 15: return tex_range_t(0.0, 0.0,  1.0,  1.0, 0, 0); // 4-way: R0
		default: assert(0);
		}
		return tex_range_t(0.0, 0.0, 1.0, 1.0); // never gets here
	}
};
struct road_plot_t : public cube_t {
	road_plot_t() {}
	road_plot_t(cube_t const &c) : cube_t(c) {}
	tex_range_t get_tex_range(float ar) const {return tex_range_t(0.0, 0.0, ar, ar);}
};

struct car_t {
	cube_t bcube;
	bool dim, dir;
	unsigned char cur_road_type, color_id;
	unsigned short cur_city, cur_road, cur_seg;
	float dz;
	car_t() : bcube(all_zeros), dim(0), dir(0), cur_road_type(0), color_id(0), cur_city(0), cur_road(0), cur_seg(0), dz(0.0) {}
	bool is_valid() const {return !bcube.is_all_zeros();}
	point get_center() const {return bcube.get_cube_center();}
};

class heightmap_query_t {
protected:
	float *heightmap;
	unsigned xsize, ysize;

public:
	flatten_op_t last_flatten_op;

	heightmap_query_t() : heightmap(nullptr), xsize(0), ysize(0) {}
	heightmap_query_t(float *hmap, unsigned xsize_, unsigned ysize_) : heightmap(hmap), xsize(xsize_), ysize(ysize_) {}
	float get_x_value(int x) const {return get_xval(x - int(xsize)/2);} // convert from center to LLC
	float get_y_value(int y) const {return get_yval(y - int(ysize)/2);}
	int get_x_pos(float x) const {return (get_xpos(x) + int(xsize)/2);}
	int get_y_pos(float y) const {return (get_ypos(y) + int(ysize)/2);}
	float  get_height(unsigned x, unsigned y) const {return heightmap[y*xsize + x];} // Note: not bounds checked
	float &get_height(unsigned x, unsigned y)       {return heightmap[y*xsize + x];} // Note: not bounds checked
	bool is_normalized_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <  x2 && y1 <  y2 && x2 <= xsize && y2 <= ysize);}
	bool is_valid_region     (unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {return (x1 <= x2 && y1 <= y2 && x2 <= xsize && y2 <= ysize);}
	bool is_inside_terrain(int x, int y) const {return (x >= 0 && y >= 0 && x < (int)xsize && y < (int)ysize);}
	
	float get_height_at(float xval, float yval) const {
		int const x(get_x_pos(xval)), y(get_y_pos(yval));
		return (is_inside_terrain(x, y) ? get_height(x, y) : OUTSIDE_TERRAIN_HEIGHT);
	}
	float any_underwater(unsigned x1, unsigned y1, unsigned x2, unsigned y2, bool check_border=0) const {
		min_eq(x2, xsize); min_eq(y2, ysize); // clamp upper bound
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (check_border && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				if (get_height(x, y) < water_plane_z) return 1;
			}
		}
		return 0;
	}
	void flatten_region_to(cube_t const c, unsigned slope_width, bool decrease_only=0) {
		flatten_region_to(get_x_pos(c.x1()), get_y_pos(c.y1()), get_x_pos(c.x2()), get_y_pos(c.y2()), slope_width, (c.z1() - ROAD_HEIGHT), decrease_only);
	}
	void flatten_region_to(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float elevation, bool decrease_only=0) {
		assert(is_valid_region(x1, y1, x2, y2));

		for (unsigned y = max((int)y1-(int)slope_width, 0); y < min(y2+slope_width, ysize); ++y) {
			for (unsigned x = max((int)x1-(int)slope_width, 0); x < min(x2+slope_width, xsize); ++x) {
				float &h(get_height(x, y));
				if (decrease_only && h < elevation) continue; // don't increase

				if (slope_width > 0) {
					float const dx(max(0, max(((int)x1 - (int)x), ((int)x - (int)x2 + 1))));
					float const dy(max(0, max(((int)y1 - (int)y), ((int)y - (int)y2 + 1))));
					h = smooth_interp(h, elevation, min(1.0f, sqrt(dx*dx + dy*dy)/slope_width));
				} else {h = elevation;}
			} // for x
		} // for y
	}
	float flatten_sloped_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float z1, float z2, bool dim, unsigned border, bool stats_only=0, bool decrease_only=0) {
		if (!stats_only) {last_flatten_op = flatten_op_t(x1, y1, x2, y2, z1, z2, dim, border);} // cache for later replay
		assert(is_valid_region(x1, y1, x2, y2));
		if (x1 == x2 || y1 == y2) return 0.0; // zero area
		float const run_len(dim ? (y2 - y1) : (x2 - x1)), denom(1.0f/max(run_len, 1.0f)), dz(z2 - z1), border_inv(1.0/border);
		int const pad(border + 1U); // pad an extra 1 texel to handle roads misaligned with the texture
		unsigned px1(x1), py1(y1), px2(x2), py2(y2);
		float tot_dz(0.0);
		
		if (dim) {
			px1 = max((int)x1-pad, 0);
			px2 = min(x2+pad, xsize);
			py1 = max((int)y1-1, 0); // pad by 1 in road dim as well to blend with edge of city
			py2 = min(y2+1, ysize);
		}
		else {
			py1 = max((int)y1-pad, 0);
			py2 = min(y2+pad, ysize);
			px1 = max((int)x1-1, 0);
			px2 = min(x2+1, xsize);
		}
		for (unsigned y = py1; y < py2; ++y) {
			for (unsigned x = px1; x < px2; ++x) {
				float const t(((dim ? int(y - y1) : int(x - x1)) + ((dz < 0.0) ? 1 : -1))*denom); // bias toward the lower zval
				float const road_z(z1 + dz*t - ROAD_HEIGHT);
				float &h(get_height(x, y));
				if (decrease_only && h < road_z) continue; // don't increase
				float new_h;

				if (border > 0) {
					unsigned const dist(dim ? max(0, max(((int)x1 - (int)x - 1), ((int)x - (int)x2))) : max(0, max(((int)y1 - (int)y - 1), ((int)y - (int)y2))));
					new_h = smooth_interp(h, road_z, dist*border_inv);
				} else {new_h = road_z;}
				tot_dz += fabs(h - new_h);
				if (!stats_only) {h = new_h;} // apply the height change
			} // for x
		} // for y
		return tot_dz;
	}
	float flatten_for_road(road_t const &road, unsigned border, bool stats_only=0, bool decrease_only=0) {
		float const z_adj(ROAD_HEIGHT + 0.5*road.get_slope_val()*(road.dim ? DY_VAL : DX_VAL)); // account for a half texel of error for sloped roads
		unsigned const rx1(get_x_pos(road.x1())), ry1(get_y_pos(road.y1())), rx2(get_x_pos(road.x2())), ry2(get_y_pos(road.y2()));
		return flatten_sloped_region(rx1, ry1, rx2, ry2, road.d[2][road.slope]-z_adj, road.d[2][!road.slope]-z_adj, road.dim, border, stats_only, decrease_only);
	}
};

class city_plot_gen_t : public heightmap_query_t {

protected:
	int last_rgi;
	rand_gen_t rgen;
	vector<rect_t> used;
	vector<cube_t> plots; // same size as used
	cube_t bcube;

	bool overlaps_used(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		rect_t const cur(x1, y1, x2, y2);
		for (vector<rect_t>::const_iterator i = used.begin(); i != used.end(); ++i) {if (i->has_overlap(cur)) return 1;} // simple linear iteration
		return 0;
	}
	cube_t add_plot(unsigned x1, unsigned y1, unsigned x2, unsigned y2, float elevation) {
		cube_t c;
		c.x1() = get_x_value(x1);
		c.x2() = get_x_value(x2);
		c.y1() = get_y_value(y1);
		c.y2() = get_y_value(y2);
		c.z1() = c.z2() = elevation;
		if (plots.empty()) {bcube = c;} else {bcube.union_with_cube(c);}
		plots.push_back(c);
		used.emplace_back(x1, y1, x2, y2);
		return c;
	}
	float get_avg_height(unsigned x1, unsigned y1, unsigned x2, unsigned y2) const {
		assert(is_normalized_region(x1, y1, x2, y2));
		float sum(0.0), denom(0.0);

		for (unsigned y = y1; y < y2; ++y) {
			for (unsigned x = x1; x < x2; ++x) {
				if (CHECK_HEIGHT_BORDER_ONLY && y != y1 && y != y2-1 && x == x1+1) {x = x2-1;} // jump to right edge
				sum   += get_height(x, y);
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
				float const delta(get_height(x, y) - avg);
				diff += delta*delta; // square the difference
			}
		}
		return diff;
	}
public:
	city_plot_gen_t() : last_rgi(0), bcube(all_zeros) {}

	void init(float *heightmap_, unsigned xsize_, unsigned ysize_) {
		heightmap = heightmap_; xsize = xsize_; ysize = ysize_;
		assert(heightmap != nullptr);
		assert(xsize > 0 && ysize > 0); // any size is okay
		if (rand_gen_index != last_rgi) {rgen.set_state(rand_gen_index, 12345); last_rgi = rand_gen_index;} // only when rand_gen_index changes
	}
	bool find_best_city_location(unsigned wmin, unsigned hmin, unsigned wmax, unsigned hmax, unsigned border, unsigned slope_width, unsigned num_samples,
		unsigned &cx1, unsigned &cy1, unsigned &cx2, unsigned &cy2)
	{
		assert(num_samples > 0);
		assert((wmax + 2*border) < xsize && (hmax + 2*border) < ysize); // otherwise the city can't fit in the map
		unsigned const num_iters(100*num_samples); // upper bound
		unsigned xend(xsize - wmax - 2*border + 1), yend(ysize - hmax - 2*border + 1); // max rect LLC, inclusive
		unsigned num_cands(0);
		float best_diff(0.0);

		for (unsigned n = 0; n < num_iters; ++n) { // find min RMS height change across N samples
			unsigned const x1(border + (rgen.rand()%xend)), y1(border + (rgen.rand()%yend));
			unsigned const x2(x1 + ((wmin == wmax) ? wmin : rgen.rand_int(wmin, wmax)));
			unsigned const y2(y1 + ((hmin == hmax) ? hmin : rgen.rand_int(hmin, hmax)));
			if (overlaps_used (x1-slope_width, y1-slope_width, x2+slope_width, y2+slope_width)) continue; // skip if plot expanded by slope_width overlaps an existing city
			if (any_underwater(x1, y1, x2, y2, CHECK_HEIGHT_BORDER_ONLY))   continue; // skip
			float const diff(get_rms_height_diff(x1, y1, x2, y2));
			if (num_cands == 0 || diff < best_diff) {cx1 = x1; cy1 = y1; cx2 = x2; cy2 = y2; best_diff = diff;}
			if (++num_cands == num_samples) break; // done
		} // for n
		if (num_cands == 0) return 0;
		//cout << "City cands: " << num_cands << ", diff: " << best_diff << ", loc: " << (cx1+cx2)/2 << "," << (cy1+cy2)/2 << endl;
		return 1; // success
	}
	float flatten_region(unsigned x1, unsigned y1, unsigned x2, unsigned y2, unsigned slope_width, float const *const height=nullptr) {
		float const elevation(height ? *height : get_avg_height(x1, y1, x2, y2));
		flatten_region_to(x1, y1, x2, y2, slope_width, elevation);
		return elevation;
	}
	bool check_plot_sphere_coll(point const &pos, float radius, bool xy_only=1) const {
		if (plots.empty()) return 0;
		point const query_pos(pos - get_query_xlate());
		if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
		return check_bcubes_sphere_coll(plots, query_pos, radius, xy_only);
	}
}; // city_plot_gen_t


class road_mat_mrg_t {

	bool inited;
	unsigned tids[NUM_RD_TIDS];

public:
	road_mat_mrg_t() : inited(0) {}

	void ensure_road_textures() {
		if (inited) return;
		timer_t timer("Load Road Textures");
		string const img_names[NUM_RD_TIDS] = {"sidewalk.jpg", "straight_road.jpg", "bend_90.jpg", "int_3_way.jpg", "int_4_way.jpg"};
		float const aniso[NUM_RD_TIDS] = {4.0, 16.0, 8.0, 8.0, 8.0};
		for (unsigned i = 0; i < NUM_RD_TIDS; ++i) {tids[i] = get_texture_by_name(("roads/" + img_names[i]), 0, 0, 1, aniso[i]);}
		inited = 1;
	}
	void set_texture(unsigned type) {
		assert(type < NUM_RD_TYPES);
		ensure_road_textures();
		select_texture(tids[type]);
	}
};

road_mat_mrg_t road_mat_mrg;


struct draw_state_t {
	shader_t s;
	vector3d xlate;
	bool use_smap, use_bmap;
protected:
	bool emit_now;

public:
	draw_state_t() : xlate(zero_vector), use_smap(0), use_bmap(0), emit_now(0) {}
	virtual void draw_unshadowed() {}
	void begin_tile(point const &pos) {emit_now = (use_smap && try_bind_tile_smap_at_point((pos + xlate), s));}

	void pre_draw(vector3d const &xlate_) {
		xlate = xlate_;
		use_smap = shadow_map_enabled();
		if (!use_smap) return;
		setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 1, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
		s.add_uniform_float("z_bias", cobj_z_bias);
		s.add_uniform_float("pcf_offset", 10.0*shadow_map_pcf_offset);
	}
	void post_draw() {
		emit_now = 0;
		if (use_smap) {s.end_shader();}
		setup_smoke_shaders(s, 0.0, 0, 0, 0, 1, 0, 0, 0, 0, use_bmap, 0, 0, 0, 0.0, 0.0, 0, 0, 1); // is_outside=1
		draw_unshadowed();
		s.end_shader();
	}
	bool check_cube_visible(cube_t const &bc) const {
		cube_t const bcx(bc + xlate);
		return (camera_pdu.cube_visible(bcx) && dist_less_than(camera_pdu.pos, bcx.closest_pt(camera_pdu.pos), get_draw_tile_dist()));
	}
}; // draw_state_t


class city_road_gen_t {

	struct range_pair_t {
		unsigned s, e; // Note: e is one past the end
		range_pair_t() : s(0), e(0) {}
		void update(unsigned v) {
			if (s == 0 && e == 0) {s = v;} // first insert
			else {assert(s < e && v >= e);} // v must strictly increase
			e = v+1; // one past the end
		}
	};

	class road_draw_state_t : public draw_state_t {
		quad_batch_draw qbd_batched[NUM_RD_TYPES];
		float ar;

	public:
		road_draw_state_t() : ar(1.0) {}

		void pre_draw(vector3d const &xlate_) {
			draw_state_t::pre_draw(xlate_);
			ar = city_params.get_road_ar();
		}
		virtual void draw_unshadowed() {
			for (unsigned i = 0; i < NUM_RD_TYPES; ++i) { // only unshadowed blocks
				road_mat_mrg.set_texture(i);
				qbd_batched[i].draw_and_clear();
			}
		}
		template<typename T> void add_road_quad(T const &r, quad_batch_draw &qbd) const {add_flat_road_quad(r, qbd, ar);} // generic flat road case
		template<> void add_road_quad(road_seg_t  const &r, quad_batch_draw &qbd) const {r.add_road_quad(qbd, ar);}

		template<typename T> void draw_road_region(vector<T> const &v, range_pair_t const &rp, quad_batch_draw &cache, unsigned type_ix) {
			assert(rp.s <= rp.e && rp.e <= v.size());
			assert(type_ix < NUM_RD_TYPES);
			
			if (cache.empty()) { // generate and cache quads
				for (unsigned i = rp.s; i < rp.e; ++i) {add_road_quad(v[i], cache);}
			}
			if (emit_now) { // draw shadow blocks directly
				road_mat_mrg.set_texture(type_ix);
				cache.draw();
			} else {qbd_batched[type_ix].add_quads(cache);} // add non-shadow blocks for drawing later
		}
	}; // road_draw_state_t

	class road_network_t {
		vector<road_t> roads; // full overlapping roads, for collisions, etc.
		vector<road_seg_t> segs; // non-overlapping road segments, for drawing with textures
		vector<road_isec_t> isecs[3]; // for drawing with textures: {4-way, 3-way, 2-way}
		vector<road_plot_t> plots; // plots of land that can hold buildings
		cube_t bcube;
		vector<road_t> segments; // reused temporary
		set<unsigned> connected_to; // vector?
		unsigned cluster_id;
		//string city_name; // future work

		static uint64_t get_tile_id_for_cube(cube_t const &c) {return get_tile_id_containing_point_no_xyoff(c.get_cube_center());}

		struct cmp_by_tile { // not the most efficient solution, but no memory overhead
			bool operator()(cube_t const &a, cube_t const &b) const {return (get_tile_id_for_cube(a) < get_tile_id_for_cube(b));}
		};
		struct tile_block_t { // collection of road parts for a given tile
			range_pair_t ranges[NUM_RD_TYPES]; // {plot, seg, isec2, isec3, isec4}
			quad_batch_draw quads[NUM_RD_TYPES];
			cube_t bcube;
			tile_block_t(cube_t const &bcube_) : bcube(bcube_) {}
		};
		vector<tile_block_t> tile_blocks;

		template<typename T> void add_tile_blocks(vector<T> &v, map<uint64_t, unsigned> &tile_to_block_map, unsigned type_ix) {
			assert(type_ix < NUM_RD_TYPES);
			sort(v.begin(), v.end(), cmp_by_tile());

			for (unsigned i = 0; i < v.size(); ++i) {
				uint64_t const tile_id(get_tile_id_for_cube(v[i]));
				auto it(tile_to_block_map.find(tile_id));
				unsigned block_id(0);
			
				if (it == tile_to_block_map.end()) { // not found, add new block
					tile_to_block_map[tile_id] = block_id = tile_blocks.size();
					tile_blocks.push_back(tile_block_t(v[i]));
				}
				else {block_id = it->second;}
				assert(block_id < tile_blocks.size());
				tile_blocks[block_id].ranges[type_ix].update(i);
				tile_blocks[block_id].bcube.union_with_cube(v[i]);
			} // for i
		}
	public:
		road_network_t() : bcube(all_zeros), cluster_id(0) {}
		road_network_t(cube_t const &bcube_) : bcube(bcube_), cluster_id(0) {bcube.d[2][1] += ROAD_HEIGHT;} // make it nonzero size
		cube_t const &get_bcube() const {return bcube;}
		void set_bcube(cube_t const &bcube_) {bcube = bcube_;}
		unsigned num_roads() const {return roads.size();}
		bool empty() const {return roads.empty();}
		void set_cluster(unsigned id) {cluster_id = id;}
		void register_connected_city(unsigned id) {connected_to.insert(id);}
		set<unsigned> const &get_connected() const {return connected_to;}

		void clear() {
			roads.clear();
			segs.clear();
			plots.clear();
			for (unsigned i = 0; i < 3; ++i) {isecs[i].clear();}
			tile_blocks.clear();
		}
		bool gen_road_grid(float road_width, float road_spacing) {
			cube_t const &region(bcube); // use our bcube as the region to process
			vector3d const size(region.get_size());
			assert(size.x > 0.0 && size.y > 0.0);
			float const half_width(0.5*road_width), zval(region.z1() + ROAD_HEIGHT);
			float const rx1(region.x1() + half_width), rx2(region.x2() - half_width), ry1(region.y1() + half_width), ry2(region.y2() - half_width); // shrink to include centerlines
			float road_pitch_x(road_width + road_spacing), road_pitch_y(road_pitch_x);
			int const num_x_roads((rx2 - rx1)/road_pitch_x), num_y_roads((ry2 - ry1)/road_pitch_y);
			road_pitch_x = 0.9999*(rx2 - rx1)/num_x_roads; // auto-calculate, round down slightly to avoid FP error
			road_pitch_y = 0.9999*(ry2 - ry1)/num_y_roads;

			// create a grid, for now; crossing roads will overlap
			for (float x = rx1; x < rx2; x += road_pitch_x) {
				roads.emplace_back(point(x, region.y1(), zval), point(x, region.y2(), zval), road_width, false);
			}
			unsigned const num_x(roads.size());

			for (float y = ry1; y < ry2; y += road_pitch_y) {
				roads.emplace_back(point(region.x1(), y, zval), point(region.x2(), y, zval), road_width, true);
			}
			unsigned const num_r(roads.size()), num_y(num_r - num_x);
			if (num_x <= 1 || num_y <= 1) {clear(); return 0;} // not enough space for roads
			bcube.x1() = roads[0      ].x1(); // actual bcube x1 from first x road
			bcube.x2() = roads[num_x-1].x2(); // actual bcube x2 from last  x road
			bcube.y1() = roads[num_x  ].y1(); // actual bcube y1 from first y road
			bcube.y2() = roads[num_r-1].y2(); // actual bcube y2 from last  y road

			// create road segments and intersections
			segs .reserve(num_x*(num_y-1) + (num_x-1)*num_y + 4); // X + Y segments, allocate one extra per side for connectors
			plots.reserve((num_x-1)*(num_y-1));

			if (num_x > 2 && num_y > 2) {
				isecs[0].reserve(4); // 2-way, always exactly 4 at each corner
				isecs[1].reserve(2*((num_x-2) + (num_y-2)) + 4); // 3-way, allocate one extra per side for connectors
				isecs[2].reserve((num_x-2)*(num_y-2)); // 4-way
			}
			for (unsigned x = 0; x < num_x; ++x) {
				for (unsigned y = num_x; y < num_r; ++y) {
					bool const FX(x == 0), FY(y == num_x), LX(x+1 == num_x), LY(y+1 == num_r);
					cube_t const &rx(roads[x]), &ry(roads[y]);
					unsigned const num_conn((!FX) + (!LX) + (!FY) + (!LY));
					if (num_conn < 2) continue; // error?
					uint8_t const conn(((!FX) << 0) | ((!LX) << 1) | ((!FY) << 2) | ((!LY) << 3)); // 1-15
					isecs[num_conn - 2].emplace_back(cube_t(rx.x1(), rx.x2(), ry.y1(), ry.y2(), zval, zval), conn); // intersections
					
					if (!LX) { // skip last y segment
						cube_t const &rxn(roads[x+1]);
						segs.emplace_back(cube_t(rx.x2(), rxn.x1(), ry.y1(), ry.y2(), zval, zval), false); // y-segments
					}
					if (!LY) { // skip last x segment
						cube_t const &ryn(roads[y+1]);
						segs.emplace_back(cube_t(rx.x1(), rx.x2(), ry.y2(), ryn.y1(), zval, zval), true); // x-segments

						if (!LX) { // skip last y segment
							cube_t const &rxn(roads[x+1]);
							plots.push_back(cube_t(rx.x2(), rxn.x1(), ry.y2(), ryn.y1(), zval, zval)); // plots between roads
						}
					}
				} // for y
			} // for x
			return 1;
		}
		void calc_bcube_from_roads() {
			if (roads.empty()) return; // no roads
			bcube = roads.front(); // first road
			for (auto r = roads.begin()+1; r != roads.end(); ++r) {bcube.union_with_cube(*r);} // skip first road
		}
	private:
		int find_conn_int_seg(cube_t const &c, bool dim, bool dir) const {
			for (unsigned i = 0; i < segs.size(); ++i) {
				road_seg_t const &s(segs[i]);
				if (s.dim == dim) continue; // not perp dim
				if (s.d[dim][dir] != bcube.d[dim][dir]) continue; // not on edge of road grid
				if (s.d[!dim][1] < c.d[!dim][0] || s.d[!dim][0] > c.d[!dim][1]) continue; // no overlap/projection in other dim
				if (c.d[!dim][0] > s.d[!dim][0] && c.d[!dim][1] < s.d[!dim][1]) return i; // c contained in segment in other dim, this is the one we want
				return -1; // partial overlap in other dim, can't split, fail
			} // for i
			return -1; // not found
		}
		int find_conn_intersection(cube_t const &c, bool dim, bool dir, unsigned num_conn) const {
			assert(num_conn == 2 || num_conn == 3); // 2 or 3 way intersections; 4 way can't have more connections
			vector<road_isec_t> const &v(isecs[num_conn-2]);
			for (unsigned i = 0; i < v.size(); ++i) {if (c.intersects_xy(v[i])) return i;}
			return -1; // not found
		}
		void make_4way_int(unsigned int3_ix) { // turn a 3-way intersection into a 4-way intersection for a connector road
			assert(int3_ix < isecs[1].size());
			road_isec_t &isec(isecs[1][int3_ix]); // bbox doesn't change, only conn changes
			isec.conn = 15; // all connected
			isecs[2].push_back(isec); // add as 4-way intersection
			isecs[1][int3_ix] = isecs[1].back(); // remove original 3-way intersection
			isecs[1].pop_back();
		}
	public:
		bool check_valid_conn_intersection(cube_t const &c, bool dim, bool dir) const {return (find_conn_int_seg(c, dim, dir) >= 0);}
		void insert_conn_intersection(cube_t const &c, bool dim, bool dir) {
			int const seg_id(find_conn_int_seg(c, dim, dir));
			assert(seg_id >= 0 && (unsigned)seg_id < segs.size());
			segs.push_back(segs[seg_id]); // clone the segment first
			segs[seg_id].d[!dim][1] = c.d[!dim][0]; // low part
			segs.back() .d[!dim][0] = c.d[!dim][1]; // high part
			cube_t ibc(segs[seg_id]); // intersection bcube
			ibc.d[!dim][0] = c.d[!dim][0]; // copy width from c
			ibc.d[!dim][1] = c.d[!dim][1];
			uint8_t const conns[4] = {7, 11, 13, 14};
			isecs[1].emplace_back(ibc, conns[2*(!dim) + dir]);
		}
		float create_connector_road(cube_t const &bcube1, cube_t const &bcube2, vector<cube_t> &blockers, road_network_t *rn1, road_network_t *rn2, heightmap_query_t &hq,
			float road_width, float conn_pos, bool dim, bool check_only=0)
		{
			bool const dir(bcube1.d[dim][0] < bcube2.d[dim][0]);
			point p1, p2;
			p1.z = bcube1.d[2][1];
			p2.z = bcube2.d[2][1];
			p1[!dim] = p2[!dim] = conn_pos;
			p1[ dim] = bcube1.d[dim][ dir];
			p2[ dim] = bcube2.d[dim][!dir];
			bool const slope((p1.z < p2.z) ^ dir);
			road_t const road(p1, p2, road_width, dim, slope);
			float const road_len(road.get_length()), delta_z(road.get_dz()), max_slope(city_params.max_road_slope);
			assert(road_len > 0.0 && delta_z >= 0.0);
			if (delta_z/road_len > max_slope) return -1.0; // slope is too high (split segments will have even higher slopes)
			unsigned const x1(hq.get_x_pos(road.x1())), y1(hq.get_y_pos(road.y1())), x2(hq.get_x_pos(road.x2())), y2(hq.get_y_pos(road.y2()));

			if (check_only) { // only need to do these checks in this case
				// FIXME: use find_conn_intersection(4)/make_4way_int() to create a new 4-way intersection on one city if one of these fails?
				if (rn1 && !rn1->check_valid_conn_intersection(road, dim,  dir)) return -1.0; // invalid, don't make any changes
				if (rn2 && !rn2->check_valid_conn_intersection(road, dim, !dir)) return -1.0;

				for (auto b = blockers.begin(); b != blockers.end(); ++b) {
					if ((rn1 && b->contains_cube(bcube1)) || (rn2 && b->contains_cube(bcube2))) continue; // skip current cities
					// create an intersection if blocker is a road, and happens to be the same elevation?
					if (b->intersects_xy(road)) return -1.0; // bad intersection, fail
				}
				if (hq.any_underwater(x1, y1, x2+1, y2+1)) return -1.0; // underwater (Note: bounds check is done here)
			}
			if (!check_only) { // create intersections and add blocker
				if (rn1) {rn1->insert_conn_intersection(road, dim,  dir);}
				if (rn2) {rn2->insert_conn_intersection(road, dim, !dir);}
				float const blocker_padding(max(city_params.road_spacing, 2.0f*city_params.road_border*max(DX_VAL, DY_VAL))); // use road_spacing?
				blockers.push_back(road);
				blockers.back().expand_by(blocker_padding); // add extra padding
			}
			if (road_len <= city_params.conn_road_seg_len) { // simple single road segment case
				if (!check_only) {roads.push_back(road);}
				return hq.flatten_sloped_region(x1, y1, x2, y2, road.d[2][slope]-ROAD_HEIGHT, road.d[2][!slope]-ROAD_HEIGHT, dim, city_params.road_border, check_only);
			}
			unsigned const num_segs(ceil(road_len/city_params.conn_road_seg_len));
			assert(num_segs > 0 && num_segs < 1000); // sanity check
			float const seg_len(road_len/num_segs);
			assert(seg_len <= city_params.conn_road_seg_len);
			road_t rs(road); // keep d[!dim][0], d[!dim][1] and dim
			rs.z1() = road.d[2][slope];
			segments.clear();
			float tot_dz(0.0);

			for (unsigned n = 0; n < num_segs; ++n) {
				rs.d[dim][1] = rs.d[dim][0] + seg_len;
				point pos;
				pos[ dim] = rs.d[dim][1];
				pos[!dim] = conn_pos;
				rs.z2()   = hq.get_height_at(pos.x, pos.y) + ROAD_HEIGHT; // terrain height at end of segment
				rs.slope  = (rs.z2() < rs.z1());
				
				if (fabs(rs.get_slope_val()) > max_slope) { // slope is too high, clamp z2 to max allowed value
					if (n+1 == num_segs) return -1.0;
					rs.z2() = rs.z1() + seg_len*max_slope*SIGN(rs.get_dz());
				}
				segments.push_back(rs);
				rs.d[dim][0] = rs.d[dim][1]; rs.z1() = rs.z2(); // shift segment end point
			} // for n
			for (auto s = segments.begin(); s != segments.end(); ++s) {
				if (s->z2() < s->z1()) {swap(s->z2(), s->z1());} // swap zvals if needed
				assert(s->is_normalized());
				tot_dz += hq.flatten_for_road(*s, city_params.road_border, check_only);
				if (!check_only) {roads.push_back(*s);}
			}
			if (!check_only) { // post-flatten pass to fix up dirt at road joints - doesn't help much
				for (auto s = segments.begin(); s != segments.end(); ++s) {hq.flatten_for_road(*s, city_params.road_border, 0, 1);} // decrease_only=1
			}
			return tot_dz; // success
		}
		void create_connector_bend(cube_t const &int_bcube, bool dx, bool dy, heightmap_query_t &hq) {
			hq.flatten_region_to(int_bcube, city_params.road_border);
			uint8_t const conns[4] = {6, 5, 10, 9};
			isecs[0].emplace_back(int_bcube, conns[2*dy + dx]);
			//blockers.push_back(int_bcube); // ???
		}
		void split_connector_roads(float road_spacing) {
			// Note: here we use segs, maybe 2-way isecs for bends, but not plots
			for (auto r = roads.begin(); r != roads.end(); ++r) {
				bool const d(r->dim), slope(r->slope);
				float const len(r->get_length());
				if (len <= road_spacing) {segs.push_back(*r); continue;} // single segment road
				assert(len > 0.0);
				unsigned const num_segs(ceil(len/road_spacing));
				float const seg_len(len/num_segs), z1(r->d[2][slope]), z2(r->d[2][!slope]); // use fixed-length segments
				assert(seg_len <= road_spacing);
				cube_t c(*r); // start by copying the road's bcube
				
				for (unsigned n = 0; n < num_segs; ++n) {
					c.d[d][1] = c.d[d][0] + seg_len;
					for (unsigned e = 0; e < 2; ++e) {c.d[2][e] = z1 + (z2 - z1)*((c.d[d][e] - r->d[d][0])/len);} // interpolate road height across segments
					if (c.z2() < c.z1()) {swap(c.z2(), c.z1());} // swap zvals if needed
					assert(c.is_normalized());
					segs.emplace_back(c, d, r->slope);
					c.d[d][0] = c.d[d][1]; // shift segment end point
				} // for n
			} // for r
		}
		void gen_tile_blocks() {
			tile_blocks.clear(); // should already be empty?
			map<uint64_t, unsigned> tile_to_block_map;
			add_tile_blocks(segs,  tile_to_block_map, TYPE_RSEG);
			add_tile_blocks(plots, tile_to_block_map, TYPE_PLOT);
			for (unsigned i = 0; i < 3; ++i) {add_tile_blocks(isecs[i], tile_to_block_map, (TYPE_ISEC2 + i));}
			//cout << "tile_to_block_map: " << tile_to_block_map.size() << ", tile_blocks: " << tile_blocks.size() << endl;
		}
		void get_road_bcubes(vector<cube_t> &bcubes) const {get_all_bcubes(roads, bcubes);}
		void get_plot_bcubes(vector<cube_t> &bcubes) const {get_all_bcubes(plots, bcubes);}

		bool check_road_sphere_coll(point const &pos, float radius, bool include_intersections, bool xy_only=1) const {
			if (roads.empty()) return 0;
			point const query_pos(pos - get_query_xlate());
			if (!check_bcube_sphere_coll(bcube, query_pos, radius, xy_only)) return 0;
			if (check_bcubes_sphere_coll(roads, query_pos, radius, xy_only)) return 1;

			if (include_intersections) { // used for global road network
				for (unsigned i = 0; i < 3; ++i) {
					if (check_bcubes_sphere_coll(isecs[i], query_pos, radius, xy_only)) return 1;
				}
			}
			return 0;
		}
		void draw(road_draw_state_t &dstate) {
			if (empty()) return;
			if (!dstate.check_cube_visible(bcube)) return; // VFC/too far

			for (auto b = tile_blocks.begin(); b != tile_blocks.end(); ++b) {
				if (!dstate.check_cube_visible(b->bcube)) continue; // VFC/too far
				dstate.begin_tile(b->bcube.get_cube_center());
				dstate.draw_road_region(segs,  b->ranges[TYPE_RSEG], b->quads[TYPE_RSEG], TYPE_RSEG); // road segments
				dstate.draw_road_region(plots, b->ranges[TYPE_PLOT], b->quads[TYPE_PLOT], TYPE_PLOT); // plots
				for (unsigned i = 0; i < 3; ++i) {dstate.draw_road_region(isecs[i], b->ranges[TYPE_ISEC2 + i], b->quads[TYPE_ISEC2 + i], (TYPE_ISEC2 + i));} // intersections
			} // for b
		}

		// cars
		bool add_car(car_t &car, rand_gen_t &rgen) const {
			if (segs.empty()) return 0; // no segments to place car on
			unsigned const seg_ix(rgen.rand()%segs.size());
			road_seg_t const &seg(segs[seg_ix]); // chose a random segment
			car.dim   = seg.dim;
			car.dir   = (rgen.rand()&1);
			car.dz    = 0.0; // flat
			car.cur_road = 0; // FIXME
			car.cur_seg  = seg_ix;
			car.cur_road_type = TYPE_RSEG;
			vector3d car_sz(vector3d(0.5*rgen.rand_uniform(0.9, 1.1), 0.18*rgen.rand_uniform(0.9, 1.1), 0.1*rgen.rand_uniform(0.9, 1.1))*city_params.road_width); // {length, width, height}
			point pos(seg.get_cube_center());
			pos[!seg.dim] += (car.dir ? -1.0 : 1.0)*(0.15*city_params.road_width); // place in right lane
			pos.z += 0.5*car_sz.z; // place above road surface
			if (seg.dim) {swap(car_sz.x, car_sz.y);}
			car.bcube.set_from_point(pos);
			car.bcube.expand_by(0.5*car_sz);
			return 1;
		}
		void update_car(car_t &car, float speed, rand_gen_t &rgen) const {
			// WRITE: move
			// WRITE: collision detection
		}
	}; // road_network_t

	vector<road_network_t> road_networks; // one per city
	road_network_t global_rn; // connects cities together; no plots
	road_draw_state_t dstate;
	rand_gen_t rgen;

	static float rgen_uniform(float val1, float val2, rand_gen_t &rgen) {return (val1 + (val2 - val1)*rgen.rand_float());}

	void assign_city_clusters() {
		vector<unsigned char> used(road_networks.size(), 0);
		vector<unsigned> pend;
		unsigned cluster_id(0);

		for (unsigned i = 0; i < used.size(); ++i) {
			if (used[i]) continue; // already used
			pend.push_back(i);

			while (!pend.empty()) {
				unsigned const id(pend.back());
				pend.pop_back();
				road_networks[id].set_cluster(cluster_id);
				used[id] = 1;
				set<unsigned> const &conn(road_networks[id].get_connected());

				for (auto c = conn.begin(); c != conn.end(); ++c) {
					assert(*c < used.size());
					if (!used[*c]) {pend.push_back(*c);}
				}
			} // while
			++cluster_id;
		} // for i
		cout << "City clusters: " << cluster_id << endl;
	}

public:
	void gen_roads(cube_t const &region, float road_width, float road_spacing) {
		//timer_t timer("Gen Roads"); // ~0.5ms
		road_networks.push_back(road_network_t(region));
		if (!road_networks.back().gen_road_grid(road_width, road_spacing)) {road_networks.pop_back();}
		//else {cout << "Roads: " << road_networks.back().num_roads() << endl;}
	}
	bool connect_two_cities(unsigned city1, unsigned city2, vector<cube_t> &blockers, heightmap_query_t &hq, float road_width) {
		assert(city1 < road_networks.size() && city2 < road_networks.size());
		assert(city1 != city2); // check for self reference
		//cout << "Connect city " << city1 << " and " << city2 << endl;
		road_network_t &rn1(road_networks[city1]), &rn2(road_networks[city2]);
		cube_t const &bcube1(rn1.get_bcube()), &bcube2(rn2.get_bcube());
		assert(!bcube1.intersects_xy(bcube2));
		float const min_edge_dist(4.0*road_width), min_jog(2.0*road_width), half_width(0.5*road_width);
		// Note: cost function should include road length, number of jogs, total elevation change, and max slope

		for (unsigned d = 0; d < 2; ++d) { // try for single segment
			float const shared_min(max(bcube1.d[d][0], bcube2.d[d][0])), shared_max(min(bcube1.d[d][1], bcube2.d[d][1]));
			
			if (shared_max - shared_min > min_edge_dist) { // can connect with single road segment in dim !d, if the terrain in between is passable
				float const val1(shared_min + 0.5*min_edge_dist), val2(shared_max - 0.5*min_edge_dist);
				float best_conn_pos(0.0), best_cost(-1.0);

				for (unsigned n = 0; n < city_params.num_conn_tries; ++n) { // make up to num_tries attempts at connecting the cities with a straight line
					float const conn_pos(rgen_uniform(val1, val2, rgen)); // chose a random new connection point and try it
					float const cost(global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, hq, road_width, conn_pos, !d, 1)); // check_only=1
					if (cost >= 0.0 && (best_cost < 0.0 || cost < best_cost)) {best_conn_pos = conn_pos; best_cost = cost;}
				}
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Single segment dim: << "d " << cost: " << best_cost << endl;
					global_rn.create_connector_road(bcube1, bcube2, blockers, &rn1, &rn2, hq, road_width, best_conn_pos, !d, 0); // check_only=0; actually make the change
					return 1;
				}
			}
		} // for d
		point const center1(bcube1.get_cube_center()), center2(bcube2.get_cube_center());
		bool const dx(center1.x < center2.x), dy(center1.y < center2.y), slope(dx ^ dy);
		cube_t const bc[2] = {bcube1, bcube2};
		
		if ((bc[dx].x1() - bc[!dx].x1() > min_jog) && (bc[dy].y1() - bc[!dy].y1() > min_jog)) {
			// connect with two road segments using a jog: Note: assumes cities are all the same size
			bool const inv_dim(rgen.rand()&1);
			//cout << "Try connect using jog in dim " << inv_dim << endl;
			cube_t bc1_conn(bcube1), bc2_conn(bcube2);
			bc1_conn.d[0][ dx] = bcube2.d[0][!dx];
			bc2_conn.d[0][!dx] = bcube1.d[0][ dx];
			bc1_conn.d[1][ dy] = bcube2.d[1][!dy];
			bc2_conn.d[1][!dy] = bcube1.d[1][ dy];

			for (unsigned d = 0; d < 2; ++d) { // x-then-y vs. y-then-x
				bool const fdim((d != 0) ^ inv_dim); // first segment dim: 0=x first, 1=y first
				bool const range_dir1(fdim ? dy : dx), range_dir2(fdim ? dx : dy);
				cube_t region1(bc1_conn), region2(bc2_conn); // regions of valid jog connectors
				region1.d[ fdim][!range_dir1] = region1.d[ fdim][ range_dir1]; // fixed edge
				region2.d[!fdim][ range_dir2] = region2.d[!fdim][!range_dir2]; // fixed edge
				region1.d[!fdim][0] += min_edge_dist; region1.d[!fdim][1] -= min_edge_dist; // variable edges
				region2.d[ fdim][0] += min_edge_dist; region1.d[ fdim][1] -= min_edge_dist; // variable edges
				float const xmin(region1.d[!fdim][0]), xmax(region1.d[!fdim][1]), ymin(region2.d[fdim][0]), ymax(region2.d[fdim][1]);
				float best_xval(0.0), best_yval(0.0), best_cost(-1.0);
				cube_t best_int_cube;

				for (unsigned n = 0; n < city_params.num_conn_tries; ++n) { // make up to num_tries attempts at connecting the cities with a single jog
					float xval(rgen_uniform(xmin, xmax, rgen)), yval(rgen_uniform(ymin, ymax, rgen));
					if (!fdim) {swap(xval, yval);}
					float const height(hq.get_height_at(xval, yval) + ROAD_HEIGHT);
					cube_t const int_cube(xval-half_width, xval+half_width, yval-half_width, yval+half_width, height, height); // the candidate intersection point
					bool has_int(0);

					for (auto b = blockers.begin(); b != blockers.end(); ++b) {
						if (b->intersects_xy(int_cube)) {has_int = 1; break;} // bad intersection, fail
					}
					//cout << TXT(dx) << TXT(dy) << TXT(slope) << TXT(fdim) << TXT(range_dir1) << TXT(range_dir2) << TXT(xval) << TXT(yval) << TXT(height) << TXT(has_int) << endl;
					if (has_int) continue; // bad intersection, fail
					float const cost1(global_rn.create_connector_road(bcube1, int_cube, blockers, &rn1, nullptr, hq, road_width, (fdim ? xval : yval),  fdim, 1)); // check_only=1
					if (cost1 < 0.0) continue; // bad segment
					if (best_cost > 0.0 && cost1 > best_cost) continue; // bound - early terminate
					float const cost2(global_rn.create_connector_road(int_cube, bcube2, blockers, nullptr, &rn2, hq, road_width, (fdim ? yval : xval), !fdim, 1)); // check_only=1
					if (cost2 < 0.0) continue; // bad segment
					float const cost(cost1 + cost2); // Note: cost function will prefer shorter routes
					if (best_cost < 0.0 || cost < best_cost) {best_xval = xval; best_yval = yval; best_int_cube = int_cube; best_cost = cost;}
				} // for n
				if (best_cost >= 0.0) { // found a candidate - use connector with lowest cost
					//cout << "Double segment cost: " << best_cost << " " << TXT(best_xval) << TXT(best_yval) << TXT(fdim) << ", int_cube: " << best_int_cube.str() << endl;
					global_rn.create_connector_bend(best_int_cube, (dx ^ fdim), (dy ^ fdim), hq); // do this first to improve flattening
					global_rn.create_connector_road(bcube1, best_int_cube, blockers, &rn1, nullptr, hq, road_width, (fdim ? best_xval : best_yval),  fdim, 0); // check_only=0
					flatten_op_t const fop(hq.last_flatten_op); // cache for reuse later during decrease_only pass
					global_rn.create_connector_road(best_int_cube, bcube2, blockers, nullptr, &rn2, hq, road_width, (fdim ? best_yval : best_xval), !fdim, 0); // check_only=0
					hq.flatten_sloped_region(fop.x1, fop.y1, fop.x2, fop.y2, fop.z1, fop.z2, fop.dim, fop.border, 0, 1); // decrease_only=1; remove any dirt that the prev road added
					hq.flatten_region_to(best_int_cube, city_params.road_border, 1); // one more pass to fix mesh that was raised above the intersection by a sloped road segment
					return 1;
				}
			} // for d
		}
		return 0;
	}
	void connect_all_cities(float *heightmap, unsigned xsize, unsigned ysize, float road_width, float road_spacing) {
		if (road_width == 0.0 || road_spacing == 0.0) return; // no roads
		unsigned const num_cities(road_networks.size());
		if (num_cities < 2) return; // not cities to connect
		timer_t timer("Connect Cities");
		heightmap_query_t hq(heightmap, xsize, ysize);
		vector<unsigned> is_conn(num_cities, 0); // start with all cities unconnected (0=unconnected, 1=connected, 2=done/connect failed
		unsigned cur_city(0), num_conn(0), num_done(0);
		vector<pair<float, unsigned>> cands;
		vector<cube_t> blockers; // existing cities and connector roads that we want to avoid intersecting

		// gather city blockers
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {
			blockers.push_back(i->get_bcube());
			blockers.back().expand_by(road_spacing); // separate roads by at least this value
		}
		// full cross-product connectivity
		for (unsigned i = 0; i < num_cities; ++i) {
			for (unsigned j = i+1; j < num_cities; ++j) {
				bool const success(connect_two_cities(i, j, blockers, hq, road_width));
				//cout << "Trying to connect city " << i << " to city " << j << ": " << success << endl;
				if (!success) continue;
				road_networks[i].register_connected_city(j);
				road_networks[j].register_connected_city(i);
			} // for j
		} // for i
		assign_city_clusters();
		global_rn.calc_bcube_from_roads();
		global_rn.split_connector_roads(road_spacing);
	}
	void gen_tile_blocks() {
		for (auto i = road_networks.begin(); i != road_networks.end(); ++i) {i->gen_tile_blocks();}
		global_rn.gen_tile_blocks();
	}
	void get_all_road_bcubes(vector<cube_t> &bcubes) const {
		global_rn.get_road_bcubes(bcubes); // not sure if this should be included
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_road_bcubes(bcubes);}
	}
	void get_all_plot_bcubes(vector<cube_t> &bcubes) const { // Note: no global_rn
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->get_plot_bcubes(bcubes);}
	}
	bool check_road_sphere_coll(point const &pos, float radius, bool xy_only) const {return global_rn.check_road_sphere_coll(pos, radius, xy_only);}

	void draw(vector3d const &xlate) { // non-const because qbd is modified
		if (road_networks.empty() && global_rn.empty()) return;
		//timer_t timer("Draw Roads");
		fgPushMatrix();
		translate_to(xlate);
		glDepthFunc(GL_LEQUAL); // helps prevent Z-fighting
		dstate.pre_draw(xlate);
		for (auto r = road_networks.begin(); r != road_networks.end(); ++r) {r->draw(dstate);}
		global_rn.draw(dstate);
		dstate.post_draw();
		glDepthFunc(GL_LESS);
		fgPopMatrix();
	}

	// cars
	bool add_car(car_t &car, rand_gen_t &rgen) const {
		if (road_networks.empty()) return 0; // no cities to add cars to
		unsigned const city(rgen.rand()%road_networks.size());
		car.cur_city = city;
		return road_networks[city].add_car(car, rgen);
	}
	void update_car(car_t &car, float speed, rand_gen_t &rgen) const {

	}
}; // city_road_gen_t


unsigned const NUM_CAR_COLORS = 8;
colorRGBA const car_colors[NUM_CAR_COLORS] = {WHITE, BLACK, DK_GRAY, RED, DK_BLUE, DK_GREEN, YELLOW, BROWN};

class car_manager_t {

	class car_draw_state_t : public draw_state_t {
		quad_batch_draw qbds[2]; // unshadowed, shadowed
	public:
		car_draw_state_t() {}

		void pre_draw(vector3d const &xlate_) {
			draw_state_t::pre_draw(xlate_);
			select_texture(WHITE_TEX);
		}
		virtual void draw_unshadowed() {qbds[0].draw_and_clear();}

		void draw_car(car_t const &car) { // Note: all quads, CCW
			if (!check_cube_visible(car.bcube)) return;
			begin_tile(car.get_center()); // enable shadows
			quad_batch_draw &qbd(qbds[emit_now]);
			assert(car.color_id < NUM_CAR_COLORS);
			colorRGBA const &color(car_colors[car.color_id]);
			cube_t const &c(car.bcube);
			bool const d(car.dim), D(car.dir);
			float fz1, fz2, bz1, bz2; // front/back top/bottom
			if (car.dz > 0.0) {fz1 = c.z1() + car.dz; fz2 = c.z2(); bz1 = c.z1(); bz2 = c.z2() - car.dz;} // going uphill
			else              {fz1 = c.z1(); fz2 = c.z2() + car.dz; bz1 = c.z1() - car.dz; bz2 = c.z2();} // going downhill or level
			point p[8];
			p[0][!d] = p[4][!d] = c.d[!d][1]; p[0][d] = p[4][d] = c.d[d][ D]; p[0].z = fz1; p[4].z = fz2; // front right
			p[1][!d] = p[5][!d] = c.d[!d][0]; p[1][d] = p[5][d] = c.d[d][ D]; p[1].z = fz1; p[5].z = fz2; // front left
			p[2][!d] = p[6][!d] = c.d[!d][0]; p[2][d] = p[6][d] = c.d[d][!D]; p[2].z = bz1; p[6].z = bz2; // back left
			p[3][!d] = p[7][!d] = c.d[!d][1]; p[3][d] = p[7][d] = c.d[d][!D]; p[3].z = bz1; p[7].z = bz2; // back right
			float const sign((d^D) ? -1.0 : 1.0);
			vector3d const top_n  (cross_product((p[2] - p[1]), (p[0] - p[1])).get_norm()*sign);
			vector3d const front_n(cross_product((p[5] - p[1]), (p[0] - p[1])).get_norm()*sign);
			vector3d right_n(all_zeros); right_n[!d] = -1.0;
			//qbd.add_quad_pts(p+0, color, -top_n); // bottom - not actually drawn
			qbd.add_quad_pts(p+4, color,  top_n); // top
			{point const pts[4] = {p[0], p[1], p[5], p[4]}; qbd.add_quad_pts(pts, color,  front_n);} // front
			{point const pts[4] = {p[2], p[3], p[7], p[6]}; qbd.add_quad_pts(pts, color, -front_n);} // back
			{point const pts[4] = {p[1], p[2], p[6], p[5]}; qbd.add_quad_pts(pts, color,  right_n);} // right
			{point const pts[4] = {p[3], p[0], p[4], p[7]}; qbd.add_quad_pts(pts, color, -right_n);} // left
			if (emit_now) {qbds[1].draw_and_clear();} // shadowed
		}
	}; // car_draw_state_t

	city_road_gen_t const &road_gen;
	vector<car_t> cars;
	car_draw_state_t dstate;
	rand_gen_t rgen;

public:
	car_manager_t(city_road_gen_t const &road_gen_) : road_gen(road_gen_) {}
	void clear() {cars.clear();}

	void init_cars(unsigned num) {
		timer_t timer("Init Cars");
		cars.reserve(num);
		
		for (unsigned n = 0; n < num; ++n) {
			car_t car;
			
			if (road_gen.add_car(car, rgen)) {
				car.color_id = rgen.rand() % NUM_CAR_COLORS;
				assert(car.is_valid());
				cars.push_back(car);
			}
		} // for i
		cout << "Cars: " << cars.size() << endl;
	}
	void next_frame(float car_speed) {
		float const speed(car_speed*fticks);
		for (auto i = cars.begin(); i != cars.end(); ++i) {road_gen.update_car(*i, speed, rgen);}
	}
	void draw(vector3d const &xlate) {
		if (cars.empty()) return;
		//timer_t timer("Draw Cars");
		dstate.xlate = xlate;
		fgPushMatrix();
		translate_to(xlate);
		dstate.pre_draw(xlate);
		for (auto i = cars.begin(); i != cars.end(); ++i) {dstate.draw_car(*i);}
		dstate.post_draw();
		fgPopMatrix();
	}
};


class city_gen_t : public city_plot_gen_t {

	city_road_gen_t road_gen;
	car_manager_t car_manager;

public:
	city_gen_t() : car_manager(road_gen) {}

	bool gen_city(city_params_t const &params, cube_t &cities_bcube) {
		unsigned x1(0), y1(0), x2(0), y2(0);
		if (!find_best_city_location(params.city_size_min, params.city_size_min, params.city_size_max, params.city_size_max,
			params.city_border, params.slope_width, params.num_samples, x1, y1, x2, y2)) return 0;
		float const elevation(flatten_region(x1, y1, x2, y2, params.slope_width));
		cube_t const pos_range(add_plot(x1, y1, x2, y2, elevation));
		if (cities_bcube.is_all_zeros()) {cities_bcube = pos_range;} else {cities_bcube.union_with_cube(pos_range);}
		if (params.roads_enabled()) {road_gen.gen_roads(pos_range, params.road_width, params.road_spacing);}
		return 1;
	}
	void gen_cities(city_params_t const &params) {
		if (params.num_cities == 0) return;
		cube_t cities_bcube(all_zeros);
		{ // open a scope
			timer_t t("Choose City Location");
			for (unsigned n = 0; n < params.num_cities; ++n) {gen_city(params, cities_bcube);}
		}
		bool const is_const_zval(cities_bcube.z1() == cities_bcube.z2());
		if (!cities_bcube.is_all_zeros()) {set_buildings_pos_range(cities_bcube, is_const_zval);}
		road_gen.connect_all_cities(heightmap, xsize, ysize, params.road_width, params.road_spacing);
		road_gen.gen_tile_blocks();
		car_manager.init_cars(city_params.num_cars);
	}
	void get_all_road_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_road_bcubes(bcubes);}
	void get_all_plot_bcubes(vector<cube_t> &bcubes) const {road_gen.get_all_plot_bcubes(bcubes);}

	bool check_city_sphere_coll(point const &pos, float radius, bool xy_only=1) const {
		if (check_plot_sphere_coll(pos, radius, xy_only)) return 1;
		return road_gen.check_road_sphere_coll(pos, radius, xy_only);
	}
	void next_frame() {car_manager.next_frame(city_params.car_speed);}

	void draw(bool shadow_only, int reflection_pass, vector3d const &xlate) { // for now, there are only roads
		if (!shadow_only && reflection_pass == 0) {road_gen.draw(xlate);} // roads don't cast shadows and aren't reflected in water
		car_manager.draw(xlate);
		// buildings are drawn through draw_buildings()
	}
}; // city_gen_t

city_gen_t city_gen;


bool parse_city_option(FILE *fp) {return city_params.read_option(fp);}
bool have_cities() {return city_params.enabled();}
float get_road_max_len() {return city_params.road_spacing;}

void gen_cities(float *heightmap, unsigned xsize, unsigned ysize) {
	if (!have_cities()) return; // nothing to do
	city_gen.init(heightmap, xsize, ysize); // only need to call once for any given heightmap
	city_gen.gen_cities(city_params);
}
void get_city_road_bcubes(vector<cube_t> &bcubes) {city_gen.get_all_road_bcubes(bcubes);}
void get_city_plot_bcubes(vector<cube_t> &bcubes) {city_gen.get_all_plot_bcubes(bcubes);}
void next_city_frame() {city_gen.next_frame();}
void draw_cities(bool shadow_only, int reflection_pass, vector3d const &xlate) {city_gen.draw(shadow_only, reflection_pass, xlate);}

bool check_city_sphere_coll(point const &pos, float radius) {
	if (!have_cities()) return 0;
	point center(pos);
	if (world_mode == WMODE_INF_TERRAIN) {center += vector3d(xoff*DX_VAL, yoff*DY_VAL, 0.0);} // apply xlate for all static objects
	return city_gen.check_city_sphere_coll(center, radius);
}
bool check_valid_scenery_pos(point const &pos, float radius) {
	if (check_buildings_sphere_coll(pos, radius, 1, 1)) return 0; // apply_tt_xlate=1, xy_only=1
	if (check_city_sphere_coll(pos, radius)) return 0;
	return 1;
}

